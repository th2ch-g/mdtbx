# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`mdtbx` はMDシミュレーション用のツールボックス。系の構築・シミュレーション実行・軌跡解析・自由エネルギー計算をサポートするCLIツール。

依存ツール: AMBER, PyMOL, OpenBabel, Gromacs, Gaussian16
力場: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## 開発コマンド

```bash
# 環境構築
pixi install

# CLIの実行
pixi run mdtbx <subcommand>
pixi run gmx ...           # mdtbx cmd gmx ... の短縮形

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

```
src/
  __main__.py    # エントリポイント: main() -> cli()
  cli.py         # argparseサブコマンドの登録・ディスパッチ
  config.py      # グローバル定数(水密度、Gaussian設定、MAXWARN等)
  logger.py      # ロガー生成ユーティリティ
  utils/         # 各サブコマンドの実装(build/analysis/general)
  cv/            # Collective Variable計算(comdist, rmsd, pca等)

pymol-plugins/
  pymol_plugins/ # PyMOLプラグイン(builder, visualizer, selector等)

example/         # 用途別のサンプルノートブック・スクリプト
install_scripts/ # Gromacs/PLUMED等の手動インストールスクリプト
```

### サブコマンドの追加パターン

各モジュール(`src/utils/*.py`, `src/cv/*.py`)は以下の2関数を実装する:

```python
def add_subcmd(subparsers):
    # argparse サブコマンドの引数定義

def run(args):
    # 実装本体
```

`cli.py` に `add_subcmd` 呼び出しと `elif sys.argv[1] == "..."` の分岐を追加して登録する。

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
