#!/usr/bin/env bash
set -euo pipefail

pixi run mdtbx tica -i feature_replica1.npy feature_replica2.npy \
    --lagtime 10 --n-components 3 -o tica.npz
pixi run mdtbx cluster -i tica.npz --n-clusters 100 \
    --seed 42 --max-iter 500 --n-jobs 1 -o clusters.npz
pixi run mdtbx msm -i clusters.npz --lagtime 10 \
    --count-mode effective -o msm.npz
