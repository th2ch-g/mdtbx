#!/bin/bash
# run.sh
# 3D-RISM を使って結晶水・保存水サイトを推定し、PDB に配置する
#
# 前提: leap (build_solution) で生成した parm7/rst7 がある
#
# 使用例:
#   bash run.sh

set -e

PRMTOP="leap.parm7"
COORD="leap.rst7"
OUTPUT="placed_water.pdb"

# -----------------------------------------------------------------------
# 基本ケース: 1D-RISM から xvv を自動生成してそのまま 3D-RISM を実行
# -----------------------------------------------------------------------
mdtbx place_solvent \
    -p ${PRMTOP} \
    -x ${COORD} \
    -o ${OUTPUT} \
    --solvent-model SPC \
    --temperature 300.0 \
    --closure kh \
    --grdspc 0.5 \
    --buffer 14.0 \
    --solvcut 14.0 \
    --threshold 1.5 \
    --exclusion-radius 2.6

echo "Basic run done -> ${OUTPUT}"

# -----------------------------------------------------------------------
# xvv 再利用ケース: 同一溶媒モデル・温度の計算を複数構造に適用する場合
# 1D-RISM は一度だけ実行して xvv を保存しておくと時間を節約できる
# -----------------------------------------------------------------------
# (初回に --keepfiles で xvv を取得しておく場合の例)
# mdtbx place_solvent \
#     -p ${PRMTOP} -x ${COORD} \
#     --keepfiles \
#     --solvent-model SPC --temperature 300.0
# # 生成された xvv を保存
# cp /tmp/rism3d_*/SPC_300.00.xvv ./SPC_300.xvv

XVV="SPC_300.xvv"
if [ -f "${XVV}" ]; then
    mdtbx place_solvent \
        -p ${PRMTOP} \
        -x ${COORD} \
        -o placed_water_reuse.pdb \
        --xvv ${XVV} \
        --temperature 300.0 \
        --closure kh \
        --threshold 1.5 \
        --exclusion-radius 2.6
    echo "xvv reuse run done -> placed_water_reuse.pdb"
fi

# -----------------------------------------------------------------------
# 上位 N サイトのみ配置するケース
# タンパク質活性部位周辺の主要な保存水だけを選びたい場合に有効
# -----------------------------------------------------------------------
mdtbx place_solvent \
    -p ${PRMTOP} \
    -x ${COORD} \
    -o placed_water_top20.pdb \
    --solvent-model SPC \
    --temperature 300.0 \
    --closure kh \
    --threshold 1.5 \
    --exclusion-radius 2.6 \
    --max-sites 20

echo "Top-20 sites run done -> placed_water_top20.pdb"

# -----------------------------------------------------------------------
# 高精度ケース: グリッドを細かくし、閾値を下げて弱い水和サイトも検出
# 計算コストは grdspc の 3 乗に比例するので注意
# -----------------------------------------------------------------------
mdtbx place_solvent \
    -p ${PRMTOP} \
    -x ${COORD} \
    -o placed_water_hires.pdb \
    --solvent-model SPC \
    --temperature 300.0 \
    --closure kh \
    --grdspc 0.3 \
    --buffer 14.0 \
    --solvcut 14.0 \
    --tolerance 1e-6 \
    --threshold 1.0 \
    --exclusion-radius 2.6

echo "High-resolution run done -> placed_water_hires.pdb"

echo "All done."
