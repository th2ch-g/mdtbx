# Resource Planning

## Evidence order

Use evidence in this order:

1. Explicit user constraints.
2. Successful prior runs with the same command, system scale, and software.
3. An approved pilot result.
4. Input file sizes and parsed system metadata.
5. Command resource class and cluster profile capabilities.
6. Conservative defaults with safety factors.

Record evidence in each request step with a source such as `approved_pilot`,
`prior_run`, or `benchmark`. Set `confidence` between zero and one only when the
evidence supports it. Never claim benchmark evidence that was not collected.

## Minimum sufficient allocation

Choose the smallest profile resource satisfying:

- required GPU count and type;
- memory estimate multiplied by the profile safety factor;
- CPU or MPI requirements;
- walltime estimate multiplied by the profile safety factor;
- node, queue, group, account, and scheduler limits.

Treat an unknown cluster or unbenchmarked production MD calculation as low
confidence. Use a separately approved pilot to estimate performance and memory.
A pilot must be scientifically useful for estimation and must not silently
become a production trajectory.

## Cluster profile discovery

Run only read-only scheduler inspection commands through `agent_probe`.
Profiles are external JSON data and must not be hard-coded in mdtbx.

For Slurm, verify partitions, node capabilities, generic resources, limits,
accounts, and required environment setup.

For Altair Grid Engine, verify queues, resource names, parallel environments,
groups, limits, and required environment setup.

For PJM, verify resource groups, node limits, elapsed-time limits, groups,
dependency syntax, and required environment setup.

If probing cannot establish a value, mark it unresolved and ask the user. Do
not manufacture cluster-specific values.
