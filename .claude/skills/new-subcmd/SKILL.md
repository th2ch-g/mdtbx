---
name: new-subcmd
description: mdtbxに新しいargparseサブコマンドを追加する
---

# new-subcmd スキル

mdtbxに新しいargparseサブコマンドを追加するスキル。

## ワークフロー

### ステップ1: ユーザーへの確認

以下をユーザーに確認する:
1. **サブコマンド名** (例: `calc_rmsd`, `gen_topology`) — CLIで `mdtbx <名前>` として呼ぶ
2. **配置場所**: 機能に応じて選択
   - `src/build/` — 系構築 (addace, amb2gro, gen_posres等と同種)
   - `src/trajectory/` — 軌跡処理 (fit, trjcat等と同種)
   - `src/analysis/` — 解析 (extract_str等と同種)
   - `src/cv/` — Collective Variable計算
   - `src/utils/` — 汎用ユーティリティ (mod_mdp, convert等と同種)
3. **引数**: 必須引数・オプション引数の名前と型
4. **処理内容**: 何をするコマンドか

### ステップ2: サブコマンドファイルの作成

`${SKILL_ROOT}/template.py` を参考に新ファイルを作成する。

必須パターン:
- `add_subcmd(subparsers)` の末尾に必ず `parser.set_defaults(func=run)` を置く
- `run(args)` が実装本体
- ロガーは `from ..logger import generate_logger` / `LOGGER = generate_logger(__name__)`
- `argparse.ArgumentDefaultsHelpFormatter` を使う
- `src/utils/` のパーサーを使う場合は `from ..utils.atom_selection_parser import AtomSelector` のように参照する

### ステップ3: cli.py への登録

`src/cli.py` に2箇所追加:

1. **importブロック** (対応するカテゴリのブロックに追加):
   ```python
   from .build import <name>       # build の場合
   from .trajectory import <name>  # trajectory の場合
   from .analysis import <name>    # analysis の場合
   from .cv import <name>          # cv の場合
   from .utils import <name>       # utils の場合
   ```

2. **add_subcmdの呼び出し** (同カテゴリのブロックに追加):
   ```python
   <name>.add_subcmd(subparsers)
   ```

### ステップ4: テストファイルの作成

`tests/test_<category>/test_<name>.py` にユニットテストを追加する:
- 外部ツール不要な純粋計算関数は直接テスト
- ファイルI/Oは `tmp_path` fixture を使用
- 軌跡が必要な場合は `conftest.py` の `trajectory_files` fixture を使用
- PyMOL/subprocess 依存は `unittest.mock.patch` でモック

### ステップ5: Lintチェック・テスト実行

```bash
pixi run r     # ruff format + lint
pixi run test  # 全テスト
```

エラーがあれば修正する。
