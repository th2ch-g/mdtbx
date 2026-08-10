Agent workflows
===============

``mdtbx`` exposes the complete command surface as JSON and adds immutable,
resource-aware execution plans for Codex, Claude Code, and other CLI agents.
Existing command behavior is unchanged unless ``--json`` or ``--dry-run`` is
used.

Machine-readable commands
-------------------------

Use ``agent_schema`` to discover command arguments and the formal version-2
request, profile, plan, state, and result JSON Schemas:

.. code-block:: console

   $ pixi run mdtbx agent_schema run_fep

Every non-Agent command accepts ``--json`` and ``--dry-run``. ``--json`` emits
one version-1 command-result envelope. ``--dry-run`` validates arguments and
emits an immutable Agent version-2 plan without executing the command:

.. code-block:: console

   $ pixi run mdtbx show_npy results.npy --dry-run

Approved execution
------------------

Create a version-2 request JSON containing a non-empty ``steps`` list.
``cwd`` and ``cluster_profile`` are optional. Workflow-level ``execution`` is
``auto``, ``local``, or ``batch``. Each step has an ``id``, ``command``,
``arguments``, and optional ``depends_on``, ``resources``, ``evidence``, and
``confidence`` fields. Per-step execution modes are rejected so one workflow
cannot silently mix login-node and compute-node execution.

.. code-block:: json

   {
     "schema_version": 2,
     "name": "analysis",
     "cwd": "/work/study",
     "cluster_profile": "/home/user/.config/mdtbx/clusters/site.json",
     "execution": "auto",
     "steps": [
       {
         "id": "tica",
         "command": "tica",
         "arguments": {
           "input": ["features-a.npy", "features-b.npy"],
           "lagtime": 10,
           "n_components": 3,
           "output": "tica.npz"
         }
       }
     ]
   }

.. code-block:: console

   $ pixi run mdtbx agent_plan --request request.json --output plan.json
   $ pixi run mdtbx agent_run --plan plan.json --approve PLAN_ID
   $ pixi run mdtbx agent_status --run RUN_DIRECTORY
   $ pixi run mdtbx agent_collect --run RUN_DIRECTORY
   $ pixi run mdtbx agent_cancel --run RUN_DIRECTORY --approve PLAN_ID

The SHA-256 ``plan_id`` covers commands, arguments, dependencies, resource
allocations, paths, filesystem snapshots, and the cluster-profile fingerprint.
Any change invalidates approval. Step identifiers, resource values, profile
capacities, and dependency graphs are strictly validated.

Existing files, directories, symlinks, and broken symlinks at output paths are
listed in the plan and require ``--approve-overwrite``. If an output or
destructive target appears or changes after planning, execution stops and a
new plan must be created. Unsafe and destructive approval flags also apply
only to the same immutable ``plan_id``.

Run state is atomically replaced after every local step and each scheduler
submission. Partial submission failures retain already-issued job IDs. mdtbx
does not automatically cancel, resubmit, or alter a scientific protocol.
``agent_collect`` persists normalized scheduler states and records result validity
and artifact size/hash metadata. Run IDs are resolved from ``.mdtbx/runs`` in the
current directory or any parent directory, so status and collection commands work
from nested workflow directories. Slurm scheduler stdout and stderr are stored next
to each generated job script instead of the submission directory.

Cluster profiles
----------------

Cluster capabilities are external version-2 JSON data. Generic incomplete
templates for Slurm, Altair Grid Engine, and PJM are stored in
``agent-profiles/``. Do not execute them before replacing all
``configure-me`` values.

On a login node, ``agent_probe`` performs read-only scheduler inspection and
returns a draft plus ``draft_id``. Saving that draft requires approval of its
exact hash:

.. code-block:: console

   $ pixi run mdtbx agent_probe --scheduler slurm
   $ pixi run mdtbx agent_profile_save \
       --draft cluster.json --output ~/.config/mdtbx/clusters/site.json \
       --approve DRAFT_ID

Profiles are not replaced unless ``agent_profile_save --replace`` is explicitly
approved. The recommended location is
``~/.config/mdtbx/clusters/NAME.json``.

Set ``MDTBX_CLUSTER_PROFILE`` to select the approved profile. Commands in
configured batch classes must run through an approved compute-node plan.
Low-confidence pilot-capable commands generate a pilot plan; production
requires a new plan and separate approval.

Shared Agent guidance
---------------------

The canonical ``md-research`` Skill is in ``agent-skills/md-research`` and is
linked into both ``.agents/skills`` and ``.claude/skills``. Codex rules and a
Claude ``PreToolUse`` hook prompt for profile persistence and execution while
blocking raw scheduler submission.
