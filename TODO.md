# TODO

## 実装タスク一覧

- [x] **[A] `find_bond`: SS-bond検出時にCYS→CYMリネーム機能追加**
- [x] **[B] `mdtbx addh` を実装する**
- [x] **[C] `mdtbx build_vacuum` を実装する**
- [x] **[D] `example/build/all_atom/protein_water.sh` を実装する（ss-bond対応）**
- [x] **[E] `example/build/all_atom/protein.sh` を実装する（vacuum対応）**
- [x] **[F] `example/build/all_atom/protein_water_membrane.sh` を実装する**
- [x] **[G] `example/build/all_atom/` でglycanに対応する（membrane系2ファイル）**
- [x] **[H] `mdtbx place_solvent` を実装する（3D-RISM使用）**

---

## 詳細仕様

### [A] `find_bond`: CYS→CYMリネーム機能追加

**対象ファイル**: `src/build/find_bond.py`

**変更内容**:
- `--output-pdb` (`-op`) オプションを新規追加
  - SS-bond検出時に対象のCYSをCYMにrename（`resn CYS` → `resn CYM`）したPDBを出力
  - PyMOLの `cmd.alter()` で `resn` を変更し `cmd.save()` で出力
- 既存の `-o` オプション（bond text出力）はそのまま維持

**出力**:
- `-o bonds.txt`: `bond SYS.1.SG SYS.2.SG` 形式のtleapコマンド（既存）
- `-op cym.pdb`: CYS→CYMにrenameしたPDB（新規）

**使用例**:
```bash
mdtbx find_bond -s input.pdb -o bonds.txt -op cym.pdb
```

---

### [B] `mdtbx addh` を実装する

**対象ファイル**: `src/build/addh.py`（現在ほぼ空）

**サブコマンド登録**: `cli.py` に import と `addh.add_subcmd(subparsers)` を追加

**引数**:
- `-s/--structure` (required): 入力構造ファイル（PDB等）
- `-o/--output_prefix` (default: `out_h`): 出力ファイルプレフィックス
- `--method` (default: `reduce`, choices: `reduce`, `pymol`): 水素付加ツール
  - `reduce`: `reduce -build {input} > {output_prefix}.pdb` を subprocess で実行
  - `pymol`: `cmd.load()` → `cmd.h_add()` → `cmd.save(f"{output_prefix}.pdb")`

**出力**: `{output_prefix}.pdb`

---

### [C] `mdtbx build_vacuum` を実装する

**対象ファイル**: `src/build/build_vacuum.py`（新規作成）

**サブコマンド登録**: `cli.py` に import と `build_vacuum.add_subcmd(subparsers)` を追加

**概要**: 水・イオン追加なし、boxsize指定なしのtleap実行（vacuum系の構築）

**引数**:
- `-i/--input` (required): 入力PDBファイル
- `-o/--outdir` (default: `./`): 出力ディレクトリ
- `--ligparam` (optional): リガンドパラメータ `FRCMOD:LIB` 形式（build_solutionと同様）
- `--addprecmd` (optional): loadpdb前に追加するtleapコマンド（例: `source leaprc.GLYCAM_06j-1`）
- `--addpostcmd` (optional): loadpdb後に追加するtleapコマンド（例: ss-bond設定）
- `--keepfiles` (flag): 中間ファイル（tleap.in, leap.log）を保持

**tleapテンプレート**: Pythonコード内にテンプレート文字列として埋め込む（`template_tleap.in`方式は使わない）。
最小構成:
```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
ADDPRECMD
LIGAND_PARAMS
LOADPDB
ADDPOSTCMD
saveamberparm SYSTEM_NAME OUT_DIR/leap.parm7 OUT_DIR/leap.rst7
savepdb SYSTEM_NAME OUT_DIR/leap.pdb
quit
```

**出力**: `{outdir}/leap.parm7`, `{outdir}/leap.rst7`, `{outdir}/leap.pdb`

---

### [D] `example/build/all_atom/protein_water.sh` を実装する

**参考**: `protein_water_ligand.sh`（リガンドなし版）

**処理フロー**:
1. `mdtbx addh -s input.pdb -o input_h`（水素付加）
2. `mdtbx addace -s input_h.pdb -o ace`
3. `mdtbx addnme -s ace.pdb -o ace_nme`
4. `mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb`（SS-bond検出＋CYMリネーム）
5. SS-bondが検出された場合は `cym.pdb` を使用し、`--addpostcmd "$(cat bonds.txt)"` を渡す
6. `mdtbx build_solution -i {input}.pdb -o ./ --ion_conc 0.15 --cation Na+ --anion Cl- --boxsize 100 100 100 [--addpostcmd "$(cat bonds.txt)"]`
7. `mdtbx amb2gro -p leap.parm7 -x leap.rst7 --type parmed`
8. `mdtbx add_ndx -g gmx.gro`
9. `mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein`
10. `mdtbx gen_posres -p gmx.top -s "protein and backbone" -o posres`
11. `mdtbx rmfile` + ファイル整理

**SS-bond有無の分岐**: `find_bond` の出力ファイル `bonds.txt` が空かどうかで分岐する

---

### [E] `example/build/all_atom/protein.sh` を実装する

**概要**: build_vacuumを使ったprotein-only（vacuum）系の構築

**処理フロー**:
1. `mdtbx addh -s input.pdb -o input_h`
2. `mdtbx addace -s input_h.pdb -o ace`
3. `mdtbx addnme -s ace.pdb -o ace_nme`
4. `mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb`
5. SS-bond分岐（[D]と同様）
6. `mdtbx build_vacuum -i {input}.pdb -o ./ [--addpostcmd "$(cat bonds.txt)"]`
7. クリーンアップ

**出力**: AMBERトポロジー (`leap.parm7`, `leap.rst7`, `leap.pdb`) のみ（Gromacsへの変換は不要）

---

### [F] `example/build/all_atom/protein_water_membrane.sh` を実装する

**参考**: `protein_water_ligand_membrane.sh`（リガンドなし版）

**処理フロー**:
1. `mdtbx addh -s input.pdb -o input_h`
2. `mdtbx addace -s input_h.pdb -o ace`
3. `mdtbx addnme -s ace.pdb -o ace_nme`
4. `mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb`（SS-bond検出＋CYMリネーム）
5. SS-bond分岐（cym.pdbを使用）
6. `mdtbx cmd packmol-memgen` で膜系構築（`protein_water_ligand_membrane.sh`から`--ligand_param`を除いたもの）
7. `mdtbx amb2gro`, `mdtbx add_ndx`, `mdtbx centering_gro`, `mdtbx gen_posres`, `mdtbx rmfile` + ファイル整理

---

### [G] `example/build/all_atom/` glycan対応

**対象ファイル**:
- `protein_water_membrane_glycan.sh`
- `protein_water_ligand_membrane_glycan.sh`

**概要**: グリカンがタンパク質と同じPDBに含まれている前提で、膜系を構築する。
GLYCAM06-jはtleap内でロードする。packmol-memgenの`--keepligs`でグリカンを保持する。

**処理フロー** (`protein_water_membrane_glycan.sh`):
1. `mdtbx addh -s input.pdb -o input_h`
2. `mdtbx addace -s input_h.pdb -o ace`
3. `mdtbx addnme -s ace.pdb -o ace_nme`
4. `mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb`
5. `mdtbx cmd packmol-memgen` で膜系構築
   - `--keepligs` でグリカン座標を保持（タンパク質と同じPDBに含まれている前提）
   - `--parametrize --gaff2` は**使わない**（GLYCAM06-jを使うため）
6. 後処理は[F]と同様

**tleapでGLYCAMをロードするための追加コマンド**:
- build_solutionの `--addprecmd "source leaprc.GLYCAM_06j-1"` を渡す

**`protein_water_ligand_membrane_glycan.sh`**:
- [G] `protein_water_membrane_glycan.sh` にリガンドパラメータ (`--ligand_param`) を追加したもの

---

### [H] `mdtbx place_solvent` を実装する

**Ref**: https://ambermd.org/tutorials/advanced/tutorial34/index.html

**対象ファイル**: `src/build/place_solvent.py`（現在ほぼ空）

**サブコマンド登録**: `cli.py` に import と `place_solvent.add_subcmd(subparsers)` を追加

**概要**: 3D-RISMによる溶媒配置。`rism3d.snglpnt` を実行し、溶媒分子の確率密度から配置候補を出力する。

**引数**:
- `-p/--prmtop` (required): AMBERトポロジーファイル (.parm7)
- `-x/--coord` (required): AMBERコーディネートファイル (.rst7)
- `-o/--output` (default: `solvent_placed.pdb`): 出力PDBファイル
- `--solvent` (default: `water`, choices: `water`): 溶媒タイプ
- `--temperature` (default: 300.0): 温度 [K]
- `--keepfiles` (flag): 中間ファイルを保持

**処理フロー**:
1. `rism3d.snglpnt` コマンドを subprocess で実行
2. 出力された3D密度分布からピーク位置を抽出してPDBとして出力

**出力**: `{output}` (PDB形式の溶媒配置構造)
