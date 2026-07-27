# Approval Contract

Approval is bound to a SHA-256 identifier and the command shown in the Agent
UI.

## Cluster profile

Show the complete proposed profile, output path, and `draft_id`. Persist it only
with:

```text
mdtbx agent_profile_save --draft PROFILE.json --output PATH --approve DRAFT_ID
```

Any content edit changes the hash and requires new approval. Refuse to replace
an existing profile unless the user separately approves the same command with
`--replace`.

## Execution plan

Show at least:

- `plan_id` and `plan_kind`;
- every step, dependency, and exact command arguments;
- input, output, existing-artifact, and destructive-target paths;
- scheduler, resource, nodes, CPUs, GPUs, memory, and walltime;
- confidence and evidence;
- unsafe, destructive, or overwrite risk flags.

Execute only with:

```text
mdtbx agent_run --plan PLAN.json --approve PLAN_ID
```

Use `--approve-unsafe`, `--approve-destructive`, or `--approve-overwrite` only
after the user explicitly approves that risk for the same plan.

An approval is invalid after any plan or profile change. Scheduler state changes
do not authorize cancellation, resubmission, resource changes, or protocol
changes.
