#!/bin/bash
# mindist.sh
# 2つの原子グループ間の全ペアから最小距離を計算する
#
# comdist (重心距離) と異なり、2グループの接触の有無を反映しやすい
# 結合部位での接触距離や、タンパク質-タンパク質界面の解析に使用する
#
# 出力: .npy (shape: [n_frames], 単位: nm)
#
# 使用例:
#   bash mindist.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# 基本ケース: 活性部位残基とリガンド間の最小距離
# -----------------------------------------------------------------------
mdtbx mindist \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "resid 50 to 70" \
    -s2 "resname LIG" \
    -o cvs/mindist.npy

echo "mindist done -> cvs/mindist.npy"

# -----------------------------------------------------------------------
# タンパク質-タンパク質界面の最小距離
# -----------------------------------------------------------------------
# mdtbx mindist \
#     -p ${TOPOLOGY} \
#     -t ${TRAJECTORY} \
#     -s1 "chainid 0" \
#     -s2 "chainid 1" \
#     -o cvs/mindist_pp.npy

echo "All done."
