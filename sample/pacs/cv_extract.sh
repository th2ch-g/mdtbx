#!/bin/bash
set -e

mkdir -p cvs/comdist
mkdir -p cvs/comvec

for i in {1..15};
do
    mdtbx comdist \
        -p trial$i/trial001/rmmol_top.gro \
        -t trial$i/trial001/prd_all.xtc \
        -s1 "protein" \
        -s2 "resname NKP" \
        -o cvs/comdist/trial$i.npy

    mdtbx comvec \
        -p trial$i/trial001/rmmol_top.gro \
        -t trial$i/trial001/prd_all.xtc \
        -s1 "protein" \
        -s2 "resname NKP" \
        -o cvs/comvec/trial$i.npy
done

echo done
