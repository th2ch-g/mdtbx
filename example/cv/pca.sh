#!/bin/bash
# pca.sh
# 主成分分析 (PCA) で集団運動モードを抽出する
#
# 自由エネルギー地形の可視化や、PACs MD の多次元 CV として使用する
# フィット選択 (-sft/-sfr) と PCA 計算選択 (-sct/-scr) を分けられる
#
# 出力:
#   - .npy: scores (shape: [n_frames, n_components])
#   - .npz: PCA metadata for PyMOL visualization
#   - .pdb: average structure after fitting
#
# 使用例:
#   bash pca.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
REFERENCE="ref.gro"

mkdir -p cvs

# -----------------------------------------------------------------------
# 基本ケース: MDtraj + scikit-learn による PCA
# backbone で重ね合わせ後、backbone の主成分を抽出
# -----------------------------------------------------------------------
mdtbx pca \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "protein and backbone" \
    -scr "protein and backbone" \
    -n 10 \
    -o cvs/pca.npy \
    -oz cvs/pca_backbone.npz \
    -oa cvs/pca_average.pdb

echo "pca done -> cvs/pca.npy  (shape: [n_frames, 10])"
echo "metadata done -> cvs/pca_backbone.npz"
echo "average structure done -> cvs/pca_average.pdb"

# -----------------------------------------------------------------------
# Gromacs gmx covar/anaeig を使う場合 (--gmx)
# 大規模系や gmx との一貫性が必要な場合に使用する
# -----------------------------------------------------------------------
# mdtbx pca \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -sft "Backbone" \
#     -sct "Backbone" \
#     --gmx \
#     -idx gmx.ndx \
#     -n 10 \
#     -o cvs/pca_gmx.npy \
#     -oz cvs/pca_gmx.npz \
#     -oa cvs/pca_gmx_average.pdb

echo "All done."
