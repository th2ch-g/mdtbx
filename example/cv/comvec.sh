#!/bin/bash
# comvec.sh
# 2つの原子グループ間の重心ベクトル (COM vector) を計算する
#
# comdist がスカラー(距離)を返すのに対し、comvec は3次元ベクトルを返す
# チャネル透過や膜挿入のような方向性のある運動の解析に適している
#
# 出力: .npy (shape: [n_frames, 3], 単位: nm)
#
# 使用例:
#   bash comvec.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# 基本ケース: MDtraj による計算
# -----------------------------------------------------------------------
mdtbx comvec \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "protein" \
    -s2 "resname LIG" \
    -o cvs/comvec.npy

echo "comvec done -> cvs/comvec.npy"

# -----------------------------------------------------------------------
# Gromacs インターフェース (--gmx)
# -----------------------------------------------------------------------
# mdtbx comvec \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s1 "Protein" \
#     -s2 "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/comvec_gmx.npy

echo "All done."
