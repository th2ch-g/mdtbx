# Kinetic analysis

`mdtbx` preserves independent trajectory boundaries in each NPZ archive. Input
feature arrays must be finite, two-dimensional `.npy` files with shape
`(frames, features)`. Discrete `.npy` trajectories must contain non-negative
integer labels.

```bash
pixi run mdtbx tica -i feature_replica1.npy feature_replica2.npy \
  --lagtime 10 --n-components 3 -o tica.npz
pixi run mdtbx cluster -i tica.npz --n-clusters 100 \
  --seed 42 --max-iter 500 --n-jobs 1 -o clusters.npz
pixi run mdtbx msm -i clusters.npz --lagtime 10 \
  --count-mode effective -o msm.npz
```

The archives use named NumPy arrays and load with `allow_pickle=False`.
`tICA.npz` contains projections, trajectory lengths, components, singular
values, and timescales. `clusters.npz` contains discrete trajectories, lengths,
centers, and convergence metadata. `msm.npz` contains transition/count
matrices, the stationary distribution, eigenvalues, timescales, and retained
state symbols.

For batch execution, place these commands in an Agent schema-v2 workflow and
select an approved cluster profile. See the project Agent workflow
documentation.

## Reference

- [deeptime](https://deeptime-ml.github.io/)
- [PaCS-Toolkit](https://github.com/Kitaolab/PaCS-Toolkit)
