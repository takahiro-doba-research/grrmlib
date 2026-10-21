# grrmlib

GRRM / Gaussian の計算で扱う分子構造や反応経路のデータを読み書き・加工する Python パッケージ。
AFIR 探索の結果から分子を分離・連結し、Gaussian SP → GRRM MIN の計算を流す workflow も含む。

## 作業のルール

- **`src/grrmlib` 以下を変更したら、`git add` の前に必ずユーザーに確認する**（変更点の要約か差分を見せて承認を得てから add / commit）
- このリポジトリでは git の user.name / user.email が未設定。コミットは過去と同じ作者で行う:
  `git -c user.name=takahiro-doba-research -c user.email=takahiro.doba.research@gmail.com commit ...`
- README.md の API の説明（`to_gv`, `EQs`, `GroupNetwork` など）は古く、今のコードと合っていない。コードを正とする

## 実行環境

```bash
docker compose build                      # requirements.txt を変えたら再ビルド
docker compose run --rm app python xxx.py # スクリプトやテストはコンテナ内で実行
```

- サービス名は `app`。リポジトリを `/work` にマウントし、`PYTHONPATH=/work/src` で読む（ソースの変更に再ビルドは不要）
- torch は CPU 版（このマシンに GPU はない）。numpy と torch は fairchem-core の制約（`numpy<2.5`, `torch~=2.13.0`）に合わせて固定している
- テストスイートはまだない。確認用スクリプトはリポジトリ直下に一時的に置いて実行し、終わったら消す

## 構成（src/grrmlib）

| パッケージ | 内容 |
|---|---|
| `core/` | `Molecule`, `Molecules`（キー → Molecule の dict）, `GroupedMolecules`, `ReactionPathNetwork` と幾何計算（`separate`, `overlay` など） |
| `io/` | ファイルの読み書き。形式ごとに `gaussian.py`, `xyz.py`, `connectable.py`, `grrm/`。`to_dataframe` で polars に変換 |
| `operations/` | 分子の連結（connectable）、電荷・多重度の組み合わせ |
| `calculators/` | ASE calculator（fairchem UMA）による SP・構造最適化 |
| `jobs/` | `GaussianJob`, `ASEJobRunner`（並列実行、結果を 1 件ずつ保存して再開可能） |
| `schedulers/` | 京大スパコン用の Slurm スケジューラなど |
| `workflows/` | `SeparationWorkflow`, `ConnectionWorkflow` |

公開する名前は `grrmlib/__init__.py` に集めている。新しい関数を足したらそこにも追加する。

## コードの約束事

- **読み書きは関数にする**: `read_xxx(path)` / `write_xxx(mol, path)`、複数ファイルは `read_xxxs(folder, basename)` / `write_xxxs(mols, folder, basename)`。クラスは状態が必要な内部の解析（`io/grrm/min.py` の `_MINOutputParser`）だけに使う
- **キーとフォルダ**: 複数ファイルの読み書きは `io/_files.py` の共通処理を使う
  - `prefixes` なし → `folder/1/2/basename`、あり → `folder/SG=1/SEQ=2/basename`（workflows はこちら）
  - 読むときは `"None"` → None、整数 → int、それ以外は str。キーが 1 段ならタプルでなく値そのもの
  - 計算結果の読み込みでは `ignore_errors=True` で、読めないファイルを空の `Molecule` にして読み進める
- **Molecule の属性**: 名前は cclib に合わせる（`scfenergy`, `atomcoords`, `vibfreqs` など）。基本属性以外は `Molecule._allowed_extra` に登録しないと `TypeError` になる
- **notes**: 原子ごとの追加情報（連結の指定など）。座標の後ろの列として読み書きし、値は int / float。書式は `operations/connectable.py` の `read_connectivity` を参照。Gaussian / GRRM の入力に書き出すときは workflows では `with_notes=False`
- **単位**: grrmlib のエネルギーは Hartree、座標は Å。ASE（eV）との変換は `calculators.HARTREE_IN_EV`
- コメント・docstring は日本語。既存のコードのスタイル（4 スペース、引数ごとの改行、`# Functions for reading.` の区切り）に合わせる
- DataFrame は polars を使う（pandas は使わない）

## 既知の問題

- `io/grrm/lup_ts.py` と `_ircirc.py` は存在しない `_optopt` を import しており動かない（公開していない）
- `operations/connectable.py` の `connect_all_and_write_as_gaussian_input` は未定義の `GaussianInputWriter` を使っており動かない
