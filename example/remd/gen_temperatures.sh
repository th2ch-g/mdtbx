#!/bin/bash
# gen_temperatures.sh
# REMD シミュレーション用の温度リストを生成する
#
# アルゴリズム: de Pablo らの方法 (virtualchemistry.org/remd-temperature-generator)
# 系の自由度・拘束条件・溶媒モデルから交換確率が目標値になる温度を逐次計算する
#
# 前提: gmx.gro / gmx.top から原子数・水分子数を確認しておく
#   水分子数:  grep -c "SOL" gmx.gro (3原子/分子なので /3 する)
#   タンパク原子数: grep -c "Protein" gmx.ndx 等で確認
#
# 使用例:
#   bash gen_temperatures.sh

set -e

# -----------------------------------------------------------------------
# 水溶液系タンパク質 (ff14SB + TIP3P/SPC 相当) の標準設定
# --pc 1 : タンパク質中の水素結合のみ拘束 (typical LINCS/SHAKE 設定)
# --wc 3 : 水は剛体モデル (TIP3P/SPC/SPC-E)
# --hff 0 : 全水素原子を含む力場 (ff14SB)
# -----------------------------------------------------------------------
NW=10000    # 水分子数
NP=3000     # タンパク質原子数
TLOW=300.0  # 最低温度 [K]
THIGH=400.0 # 最高温度 [K]
PDES=0.25   # 目標交換確率 (0.2-0.3 が一般的)

echo "=== Standard protein-water REMD ==="
mdtbx gen_temperatures \
    --pdes ${PDES} \
    --tlow ${TLOW} \
    --thigh ${THIGH} \
    --nw ${NW} \
    --np ${NP} \
    --pc 1 \
    --wc 3 \
    --hff 0

# -----------------------------------------------------------------------
# 小分子・ペプチド系: 原子数が少なく交換確率が高くなりやすい
# レプリカ数が少なくなるので温度範囲を絞るか pdes を下げて調整する
# -----------------------------------------------------------------------
echo ""
echo "=== Small peptide REMD ==="
mdtbx gen_temperatures \
    --pdes 0.20 \
    --tlow 280.0 \
    --thigh 380.0 \
    --nw 3000 \
    --np 300 \
    --pc 1 \
    --wc 3 \
    --hff 0

# -----------------------------------------------------------------------
# 全原子フレキシブル (拘束なし) 設定
# --pc 0 : タンパク質拘束なし
# --wc 0 : 水も完全フレキシブル (SPC/Ef 等)
# 自由度が増えるためレプリカ数が増える傾向にある
# -----------------------------------------------------------------------
echo ""
echo "=== Fully flexible (no constraints) ==="
mdtbx gen_temperatures \
    --pdes 0.25 \
    --tlow 300.0 \
    --thigh 400.0 \
    --nw 10000 \
    --np 3000 \
    --pc 0 \
    --wc 0 \
    --hff 0

echo "All done."
