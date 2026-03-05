#!/usr/bin/env python3

import numpy as np
from pathlib import Path

n_frame = 100

# (directory, last_dim) pairs: comdist is scalar (1), comvec is 3D vector (3)
targets = [
    ("cvs/comdist/", 1),
    ("cvs/comvec/", 3),
]

for cv_dir, ndim in targets:
    for npy in Path(cv_dir).glob("*.npy"):
        a = np.load(npy)
        print(f"{npy}: {a.shape}", end=" -> ")

        total_frame = a.shape[0] - 1  # drop first frame
        a = a[1:].reshape(total_frame // n_frame, n_frame, ndim)

        print(a.shape)
        np.save(f"{npy.stem}_reshaped.npy", a)
