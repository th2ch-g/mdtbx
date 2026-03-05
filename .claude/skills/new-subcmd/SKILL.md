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
2. **配置場所**: `src/*/<名前>.py`
3. **引数**: 必須引数・オプション引数の名前と型
4. **処理内容**: 何をするコマンドか

### ステップ2: サブコマンドファイルの作成

`${SKILL_ROOT}/template.py` を参考に、以下のパスに新ファイルを作成する:
- utils の場合: `$TOOLS/mdtbx/src/utils/<name>.py`
- cv の場合: `$TOOLS/mdtbx/src/cv/<name>.py`

必須パターン:
- `add_subcmd(subparsers)` の末尾に必ず `parser.set_defaults(func=run)` を置く
- `run(args)` が実装本体
- ロガーは `from ..logger import generate_logger` / `LOGGER = generate_logger(__name__)`
- `argparse.ArgumentDefaultsHelpFormatter` を使う

### ステップ3: cli.py への登録

`$TOOLS/mdtbx/src/cli.py` に2箇所追加:

1. **importブロック** (utils/cvの対応するブロックに追加):
   ```python
   from .utils import <name>   # utils の場合
   from .cv import <name>      # cv の場合
   ```

2. **add_subcmdの呼び出し** (同カテゴリのブロックに追加):
   ```python
   <name>.add_subcmd(subparsers)
   ```

### ステップ4: Lintチェック

```bash
cd $TOOLS/mdtbx && pixi run r
```

エラーがあれば修正する。
