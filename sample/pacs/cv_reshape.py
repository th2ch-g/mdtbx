#!/usr/bin/env python3

import numpy as np
from pathlib import Path

n_frame = 100

for npy in Path("cvs/comdist/").glob("*.npy"):
    a = np.load(npy)

    print(a.shape)

    total_frame = a.shape[0]
    total_frame = total_frame - 1

    a = a[1:].reshape(total_frame // n_frame, n_frame, 1)

    print(a.shape)

    prefix = npy.stem
    np.save(f"{prefix}_reshaped.npy", a)


for npy in Path("cvs/comvec/").glob("*.npy"):
    a = np.load(npy)

    print(a.shape)

    total_frame = a.shape[0]
    total_frame = total_frame - 1

    a = a[1:].reshape(total_frame // n_frame, n_frame, 3)

    print(a.shape)

    prefix = npy.stem
    np.save(f"{prefix}_reshaped.npy", a)
