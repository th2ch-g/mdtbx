#!/bin/bash
# rmsd.sh
# reference 構造に対する RMSD を計算する
#
# フィット選択 (-sft/-sfr) と RMSD 計算選択 (-sct/-scr) を分けられる
# 例: backbone で重ね合わせたうえで活性部位ループの RMSD を算出
#
# 出力: .npy (shape: [n_frames], 単位: nm)
#
# 使用例:
#   bash rmsd.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
REFERENCE="ref.gro"

mkdir -p cvs

# -----------------------------------------------------------------------
# backbone 全体の RMSD
# -----------------------------------------------------------------------
mdtbx rmsd \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "protein and backbone" \
    -scr "protein and backbone" \
    -o cvs/rmsd_backbone.npy

echo "rmsd backbone done -> cvs/rmsd_backbone.npy"

# -----------------------------------------------------------------------
# 活性部位ループの RMSD (backbone で重ね合わせ後にループ領域を計算)
# -sft/-sfr: 重ね合わせに使う選択 (グローバルフィット)
# -sct/-scr: RMSD を計算する選択 (局所的な変化を検出)
# -----------------------------------------------------------------------
mdtbx rmsd \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "resid 100 to 120 and backbone" \
    -scr "resid 100 to 120 and backbone" \
    -o cvs/rmsd_loop.npy

echo "rmsd loop done -> cvs/rmsd_loop.npy"

echo "All done."
