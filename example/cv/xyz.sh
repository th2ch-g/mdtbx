#!/bin/bash
# xyz.sh
# 特定の原子群の XYZ 座標時系列を抽出する
#
# スカラー CV では失われる空間情報が必要な場合や、
# カスタム CV のための前処理として使用する
#
# 出力: .npy (shape: [n_frames, n_atoms, 3], 単位: nm)
#
# 使用例:
#   bash xyz.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# リガンドの全原子座標
# -----------------------------------------------------------------------
mdtbx xyz \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    -o cvs/lig_xyz.npy

echo "xyz done -> cvs/lig_xyz.npy"

# -----------------------------------------------------------------------
# 活性部位 Cα 座標 (特定残基)
# -----------------------------------------------------------------------
mdtbx xyz \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resid 50 to 70 and name CA" \
    -o cvs/active_site_ca_xyz.npy

echo "xyz done -> cvs/active_site_ca_xyz.npy"

echo "All done."
