# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code)
when working with code in this repository.

## Overview

`mdtbx` はMDシミュレーション用のツールボックス。系の構築・シミュレーション実行・軌跡解析・自由エネルギー計算をサポートするCLIツール。

依存ツール: AMBER, PyMOL, OpenBabel, Gromacs, Gaussian16
力場: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## 設計哲学

各サブコマンドは**1つの機能のみ**を担当し、なるべく小さく保つ。複雑なタスクは複数サブコマンドを連携させて実現する。

```bash
# Good: combine small subcommands to build a pipeline
pixi run mdtbx addh --input protein.pdb --output protein_h.pdb
pixi run mdtbx gen_am1bcc --input ligand.mol2 --output ligand.mol2
pixi run mdtbx build_solution --receptor protein_h.pdb --ligand ligand.mol2
pixi run mdtbx fit --trajectory md.xtc --output fitted.xtc
pixi run mdtbx rmsd --trajectory fitted.xtc --output rmsd.npy
```

新しいサブコマンドを追加する際は「既存コマンドの組み合わせで実現できないか」を先に検討すること。

## 開発コマンド

```bash
# 環境構築
pixi install

# CLIの実行
pixi run mdtbx <subcommand>
pixi run gmx ...           # gromacs gmx バイナリを直接実行(pixi環境のgromacsを使用)

# テスト
pixi run test              # 全テスト実行
pixi run test-fast         # 最初の失敗で停止 (-x)

# コードフォーマット・Lint
pixi run r                 # ruff format + ruff check を一括実行
pixi run ruff-format       # フォーマットのみ
pixi run ruff-lint         # Lintのみ

# 更新
pixi run update            # git pull && pixi install

# PyMOL設定
pixi run pymolrc           # ~/.pymolrcを生成

# JupyterLab (リモート)
pixi run jupyter_remote
```

## アーキテクチャ

```text
src/
  __main__.py    # エントリポイント: main() -> cli()
  cli.py         # argparseサブコマンドの登録・ディスパッチ
  config.py      # グローバル定数(水密度、Gaussian設定、MAXWARN等)
  logger.py      # ロガー生成ユーティリティ
  utils/         # 汎用ユーティリティ(mod_mdp, convert, rmfile, cmd, shell_hook, show_mdtraj, show_npy, partial_tempering)
               # ※ atom_selection_parser.py, parse_top.py はサブコマンドでなくライブラリユーティリティ
  build/         # 系構築サブコマンド(addace, addh, addnme, add_ndx, mv_crds_mol2, calc_ion_conc, centering_gro,
               #   find_bond, gen_am1bcc, gen_resp, gen_modres_am1bcc, gen_modres_resp, gen_posres,
               #   gen_distres, modeling_cf, amb2gro, build_solution, build_vacuum, place_solvent, gen_temperatures)
               # ※ mutate.py は未登録(cli.pyへの追加が必要)
  trajectory/    # 軌跡処理サブコマンド(fit, trjcat, pacs_trjcat, print_perf)
               # ※ opt_perf.py は未登録
  analysis/      # 解析サブコマンド(extract_str, extract_ave_str)
               # ※ contactmap.py, distmat.py は未登録
  cv/            # Collective Variable計算(comdist, comvec, densmap, mindist, rmsd, rmsf, xyz, pca)

tests/
  conftest.py    # 共有fixture・PyMOLモック設定
  fixtures/      # テストデータ(sample.mdp, sample.top, sample.pdb)
  test_utils/    # src/utils/ のテスト
  test_build/    # src/build/ のテスト
  test_trajectory/ # src/trajectory/ のテスト
  test_analysis/ # src/analysis/ のテスト
  test_cv/       # src/cv/ のテスト
  test_cli.py    # 全サブコマンドのCLI登録確認

pymol-plugins/
  pymol_plugins/ # PyMOLプラグイン(builder, visualizer, selector等)

example/         # 用途別のサンプルノートブック・スクリプト
install_scripts/ # Gromacs/PLUMED等の手動インストールスクリプト
```

### サブコマンドの追加パターン

各モジュール(`src/build/*.py`, `src/trajectory/*.py`, `src/analysis/*.py`,
`src/cv/*.py`, `src/utils/*.py`)は以下の2関数を実装する:

```python
def add_subcmd(subparsers):
    # argparse サブコマンドの引数定義

def run(args):
    # 実装本体
```

`cli.py` に以下の2箇所を追加して登録する:

```python
# 1. importブロック (カテゴリに応じて選択)
from .build import <name>       # 系構築
from .trajectory import <name>  # 軌跡処理
from .analysis import <name>    # 解析
from .cv import <name>          # CV計算
from .utils import <name>       # 汎用

# 2. add_subcmdの呼び出し
<name>.add_subcmd(subparsers)
```

モジュール内で `src/utils/` のパーサーを使う場合は `..utils.` で参照する:

```python
from ..utils.atom_selection_parser import AtomSelector
from ..utils.parse_top import GromacsTopologyParser
```

### 設定 (`src/config.py`)

- `MAXWARN`: grompp の最大警告数
- `GAUSSIAN_CMD`, `STRUCTURE_OPTIMIZATION`, `SINGLE_POINT_CALCULATION`: Gaussian設定
- 各水モデル(TIP3P/TIP4P/TIP5P/OPC)の密度・体積定数
- 起動時に `.pixi/envs/default/bin` をPATHに追加する

## 環境管理

- パッケージ管理: `pixi`(conda + pip の混在)
- Python: 3.10固定
- `pixi.lock` で再現性を保証
- Docker対応: `Dockerfile` でコンテナビルドも可能
