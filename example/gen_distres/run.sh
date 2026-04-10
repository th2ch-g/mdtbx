#!/bin/bash
# run.sh
# Gromacs 用距離拘束ファイル (distres.itp) を生成し、topology.top に組み込む
#
# 前提: amb2gro または centering_gro で gmx.gro / gmx.top が生成済みであること
#
# 出力:
#   distres.itp  - [ intermolecular_interactions ] を含む拘束定義
#   gmx.top      - distres.itp の #include が末尾に追記される
#
# NOTE:
#   - 距離の単位は nm (0.26 nm = 2.6 Å)
#   - up2 > up1 > lo の順に設定する
#   - 各ペアに対して MDtraj atom selection 言語で1原子を一意に指定する
#     例: "resname MGX and name O2 and resid 200"
#   - [ intermolecular_interactions ] は異なる分子グループ間にのみ適用される
#     (同一分子内拘束は [ bonds ] に記述する)
#
# 使用例:
#   bash run.sh

set -e

GRO="gmx.gro"
TOP="gmx.top"

# -----------------------------------------------------------------------
# 基本ケース: 金属イオン(MGX)と活性部位残基の距離拘束
# タンパク質-金属間の配位結合を維持するための拘束
# -----------------------------------------------------------------------

# 拘束するペアを配列で管理すると可読性が上がる
# 各ペア: (lo, up1, up2, sel1, sel2) で記述
# lo=0.0 は下限なし (距離は非負なのでデフォルト 0.0 でよい)

# sel1 と sel2 にセットされる MDtraj atom selection 文字列
# カンマ区切りで複数ペアを指定できる
MGX_O1="resname MGX and name O1"
MGX_O2="resname MGX and name O2"
GLY_N="resname GLY and name N and resid 10"
GLY_N2="resname GLY and name N and resid 20"
HID_ND1="resname HID and name ND1"
HID_N="resname HID and name N"
GLH_OE1="resname GLH and name OE1"
GLH_OE2="resname GLH and name OE2"

mdtbx gen_distres \
    -g ${GRO} \
    -p ${TOP} \
    -lo 0.0 \
    -up1 0.26 0.26 0.28 0.30 0.30 \
    -up2 0.28 0.28 0.30 0.32 0.32 \
    -o distres \
    -s1 "${MGX_O2},${MGX_O2}, ${MGX_O1}, ${MGX_O2},${GLY_N2}" \
    -s2 "${GLY_N}, ${HID_ND1},${GLH_OE2},${HID_N}, ${GLH_OE1}"

echo "distres.itp generated and appended to ${TOP}"

# -----------------------------------------------------------------------
# mdp でのアクティベーション方法
# mdp ファイルに以下を追記して拘束を有効化する:
#   define = -DDISTRES -DDISTRES_FC=1000
# DISTRES_FC は力定数 (kJ/mol/nm^2)
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# ユニフォームバウンドケース: 全ペアに同じ上下限を適用
# -lo, -up1, -up2 に単一値を渡すと全ペアに適用される
# -----------------------------------------------------------------------
# mdtbx gen_distres \
#     -g ${GRO} \
#     -p ${TOP} \
#     -lo 0.0 \
#     -up1 0.28 \
#     -up2 0.30 \
#     -o distres_uniform \
#     -s1 "${MGX_O2},${MGX_O2}" \
#     -s2 "${GLY_N}, ${HID_ND1}"

echo "All done."
