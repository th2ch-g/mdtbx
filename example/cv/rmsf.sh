#!/bin/bash
# rmsf.sh
# 残基ごと・原子ごとの揺らぎ (RMSF) を計算する
#
# B-factor に対応する指標; 柔軟な領域や剛直な領域の同定に使用する
# --resolution residue: 残基ごとに集約 (デフォルト)
# --resolution atom   : 原子ごとに出力
#
# 出力: .npy (shape: [n_residues] or [n_atoms], 単位: nm)
#
# 使用例:
#   bash rmsf.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# 残基単位の RMSF (B-factor 相当)
# -----------------------------------------------------------------------
mdtbx rmsf \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    --selection "protein" \
    --resolution residue \
    -o cvs/rmsf_residue.npy

echo "rmsf residue done -> cvs/rmsf_residue.npy"

# -----------------------------------------------------------------------
# 原子単位の RMSF (backbone のみ)
# -----------------------------------------------------------------------
mdtbx rmsf \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    --selection "protein and backbone" \
    --resolution atom \
    -o cvs/rmsf_atom.npy

echo "rmsf atom done -> cvs/rmsf_atom.npy"

# -----------------------------------------------------------------------
# Gromacs gmx rmsf を使う場合 (--gmx)
# -----------------------------------------------------------------------
# mdtbx rmsf \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     --selection "Backbone" \
#     --resolution residue \
#     --gmx \
#     -o cvs/rmsf_gmx.npy

echo "All done."
