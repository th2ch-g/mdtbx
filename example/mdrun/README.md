# mdrun

## T-REMD MDP templates

- `example/mdp/solution/tremd.mdp`
- `example/mdp/membrane/tremd.mdp`

## Distance REUS with PLUMED

`example/remd/setup_reus_dist.sh` prepares one replica directory per line in
`example/remd/reus_distances.txt`. Run it from a directory containing
`gmx.top`, `index.ndx`, and the listed structures:

```bash
bash /path/to/mdtbx/example/remd/setup_reus_dist.sh
```

The `index.ndx` file must define `TARGET1` and `TARGET2`. Their center-of-mass
distance is restrained by `example/plumed/reus_dist.dat`; the target distance
is read in nm from `reus_distances.txt`. The generated `plumed.dat` writes the
distance and restraint bias to `COLVAR`.

Override `TARGETS_FILE`, `TEMPLATE_MDP`, `TEMPLATE_PLUMED`, `SUBMIT_OUTPUT`, or
`REPLEX` through environment variables when using custom inputs.

## Reference

- [A Guide to CUDA Graphs in GROMACS 2023](https://developer.nvidia.com/blog/a-guide-to-cuda-graphs-in-gromacs-2023/)
- [Massively Improved Multi-node NVIDIA GPU Scalability with GROMAC](https://developer.nvidia.com/blog/massively-improved-multi-node-nvidia-gpu-scalability-with-gromacs/)
