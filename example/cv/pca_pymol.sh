#!/bin/bash
# pca_pymol.sh
# Backbone PCA を計算し、PyMOL で主成分ベクトルを矢印表示する
#
# 使用例:
#   bash pca_pymol.sh

set -e

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${REPO_ROOT}"

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
REFERENCE="ref.gro"

mkdir -p cvs

mdtbx pca \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "protein and backbone" \
    -scr "protein and backbone" \
    -n 3 \
    -o cvs/pca.npy \
    -oz cvs/pca_backbone.npz \
    -oa cvs/pca_average.pdb

pymol <<EOF
load cvs/pca_average.pdb, avg
hide everything, avg
show cartoon, avg
color gray80, avg
show_pca_mode_from_npz cvs/pca_backbone.npz, 1, 5.0, 5, CA, 0.10, cyan red
zoom avg
EOF
