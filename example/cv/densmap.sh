#!/bin/bash
# densmap.sh
# 特定の原子群の 2D 密度マップ (ヒストグラム) を計算する
#
# 膜系でのリガンド分布や、特定平面への投影密度の可視化に使用する
# --axis: 投影する平面 (xy / xz / yz)
# --bins: 各軸のビン数
#
# 出力: .npy (object array [counts, edges0, edges1])
#   counts  : shape [bins, bins]  - 各セルの頻度
#   edges0/1: shape [bins+1]      - ビンの境界
#
# 使用例:
#   bash densmap.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# xy 平面への投影 (膜面内の分布)
# -----------------------------------------------------------------------
mdtbx densmap \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    --axis xy \
    --bins 100 \
    -o cvs/densmap_xy.npy

echo "densmap xy done -> cvs/densmap_xy.npy"

# -----------------------------------------------------------------------
# xz 平面への投影 (膜の厚さ方向を含む断面)
# -----------------------------------------------------------------------
mdtbx densmap \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    --axis xz \
    --bins 100 \
    -o cvs/densmap_xz.npy

echo "densmap xz done -> cvs/densmap_xz.npy"

# -----------------------------------------------------------------------
# Gromacs gmx densmap を使う場合 (--gmx)
# -----------------------------------------------------------------------
# mdtbx densmap \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/densmap_gmx.npy

echo "All done."
