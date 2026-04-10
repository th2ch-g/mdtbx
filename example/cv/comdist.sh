#!/bin/bash
# comdist.sh
# 2つの原子グループ間の重心距離 (COM distance) を計算する
#
# 出力: .npy (shape: [n_frames], 単位: nm)
#
# 使用例:
#   bash comdist.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# 基本ケース: MDtraj による計算
# protein と リガンド間の重心距離
# -----------------------------------------------------------------------
mdtbx comdist \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "protein" \
    -s2 "resname LIG" \
    -o cvs/comdist.npy

echo "comdist done -> cvs/comdist.npy"

# -----------------------------------------------------------------------
# Gromacs インターフェース (--gmx): 大規模系で高速
# -s1/-s2 に ndx グループ名を指定する
# -----------------------------------------------------------------------
# mdtbx comdist \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s1 "Protein" \
#     -s2 "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/comdist_gmx.npy

echo "All done."
