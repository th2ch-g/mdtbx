Agent workflows
===============

``mdtbx`` exposes the complete command surface as JSON and adds immutable,
resource-aware execution plans for Codex, Claude Code, and other CLI agents.
Existing command behavior is unchanged unless ``--json`` or ``--dry-run`` is
used.

Machine-readable commands
-------------------------

Use ``agent_schema`` to discover command arguments:

.. code-block:: console

   $ pixi run mdtbx agent_schema run_fep

Every non-Agent command accepts ``--json`` and ``--dry-run``. ``--json`` emits
one result envelope. ``--dry-run`` validates the arguments and emits an
immutable plan without executing the command:

.. code-block:: console

   $ pixi run mdtbx show_npy results.npy --dry-run

Approved execution
------------------

Create a version-1 request JSON containing ``cwd``, ``cluster_profile``, and a
non-empty ``steps`` list. Each step has an ``id``, ``command``, ``arguments``,
and optional ``depends_on``, ``resources``, ``execution``, ``evidence``, and
``confidence`` fields.

.. code-block:: console

   $ pixi run mdtbx agent_plan --request request.json --output plan.json
   $ pixi run mdtbx agent_run --plan plan.json --approve PLAN_ID
   $ pixi run mdtbx agent_status --run RUN_DIRECTORY
   $ pixi run mdtbx agent_collect --run RUN_DIRECTORY

The SHA-256 ``plan_id`` covers all commands, arguments, dependencies, resource
allocations, paths, and the external cluster-profile fingerprint. Any change
invalidates approval.

Existing output files are listed in the plan and require
``--approve-overwrite``. Unsafe, destructive, and overwrite approval flags
apply only to the same immutable ``plan_id``.

Cluster profiles
----------------

Cluster capabilities are external JSON data. Generic incomplete templates for
Slurm, Altair Grid Engine, and PJM are stored in ``agent-profiles/``. Do not
execute them before replacing all ``configure-me`` values.

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

Set ``MDTBX_CLUSTER_PROFILE`` to select the approved profile. Heavy commands
must be submitted through an approved compute-node plan. Low-confidence
pilot-capable commands generate a pilot plan; production requires a new plan
and a separate approval.

Shared Agent guidance
---------------------

The canonical ``md-research`` Skill is in ``agent-skills/md-research`` and is
linked into both ``.agents/skills`` and ``.claude/skills``. Codex rules and a
Claude ``PreToolUse`` hook prompt for profile persistence and execution while
blocking raw scheduler submission.
