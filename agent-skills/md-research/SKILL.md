---
name: md-research
description: Plan, approve, execute, and collect reproducible molecular-dynamics research workflows with mdtbx. Use for MD system preparation, pilot simulations, production calculations, trajectory analysis, free-energy work, cluster resource discovery, scheduler submission, or any request that should run through Slurm, Altair Grid Engine, or PJM.
---

# MD Research

Use the typed mdtbx Agent CLI to turn a scientific request into an immutable,
resource-aware plan. Keep scientific protocol decisions and compute allocation
visible to the user before any execution or submission.

## Workflow

1. Read `AGENTS.md`, `CLAUDE.md`, and the relevant project files.
2. Run `mdtbx agent_schema [command]` before constructing command arguments.
3. Inspect input sizes, topology, trajectory metadata, prior run results, and
   relevant primary literature or official software documentation. Ask before
   downloading or installing executables, models, force fields, or datasets.
   Do not run heavy calculations during inspection.
4. Resolve a cluster profile:
   - Use `MDTBX_CLUSTER_PROFILE` or an explicitly supplied profile when present.
   - Otherwise run `mdtbx agent_probe` on a login node.
   - Treat `draft_profile.resources[].incomplete=true` as unresolved.
   - Present the proposed profile and `draft_id`; save it only after the user
     approves that exact hash and output path with `mdtbx agent_profile_save`.
   - Prefer `~/.config/mdtbx/clusters/NAME.json`. Never replace an existing
     profile unless the user approves `--replace`.
5. Create a schema-version-2 request and run `mdtbx agent_plan`.
6. Present the plan kind, commands, dependencies, inputs, artifacts, resource
   allocation, evidence, confidence, risks, and exact `plan_id`.
7. Ask the user to approve that exact plan. Do not infer approval from an
   earlier request or approval of a different plan.
8. After approval, run only
   `mdtbx agent_run --plan PLAN --approve PLAN_ID`. Add unsafe or destructive
   approval flags only when separately approved. If the plan lists existing
   artifacts, require separate approval for `--approve-overwrite`.
9. Poll with `mdtbx agent_status`; obey the profile polling interval. Collect
   with `mdtbx agent_collect`.
10. Summarize provenance, results, limitations, and the next proposed plan.

## Scientific and Execution Boundaries

- Run heavy calculations only through an approved batch plan on a compute node.
- Keep local execution to lightweight inspection, preparation, and analysis
  allowed by the selected profile.
- Select the minimum resource that satisfies the evidence-based estimate.
- When confidence is below profile policy and the command supports a pilot,
  propose a separately approved pilot. Collect it before proposing production.
- Never convert pilot approval into production approval.
- Treat an unknown cluster or unbenchmarked production MD calculation as low
  confidence and require a separately approved pilot.
- Never alter the protocol, resources, dependencies, output paths, or command
  arguments after `plan_id` creation. Re-plan and request approval.
- Never submit with raw `sbatch`, `qsub`, or `pjsub`.
- Never cancel, resubmit, expand sampling, or recover automatically.
- Never write an incomplete or unapproved cluster profile.

Read [references/resource-planning.md](references/resource-planning.md) when
estimating resources or preparing a profile. Read
[references/approval-contract.md](references/approval-contract.md) before
execution, scheduler submission, or profile persistence.
