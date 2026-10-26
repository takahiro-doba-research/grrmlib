# grrmlib

GRRM / Gaussian の計算で扱う分子構造や反応経路のデータを読み書き・加工する Python パッケージ。
AFIR 探索の結果から分子を分離・連結して xyz に書き、そこから Gaussian SP → GRRM MIN の入力を作る pipeline も含む。

- **パブリックリポジトリ**（GitHub の `takahiro-doba-research/grrmlib`）
- 開発はプライベートの `alphachem`（`~/calculation/alphachem`。grrmlib を履歴ごと引き継いで作った）で行い、
  公開してよい部分をユーザーに確認してから grrmlib に移す。UMA / ASE による計算（alphachem の `calculators/`、
  `run_ase`、`read_ase_result(s)`）は grrmlib には入れない（v0.3.0 で削除した。v0.3.1 の workflows も計算のループ `run_sp` / `run_opt` を除いて移した）

## 作業のルール

- **`src/grrmlib` 以下を変更したら、`git add` の前に必ずユーザーに確認する**（変更点の要約か差分を見せて承認を得てから add / commit）
- このリポジトリでは git の user.name / user.email が未設定。コミットは過去と同じ作者で行う:
  `git -c user.name=takahiro-doba-research -c user.email=takahiro.doba.research@gmail.com commit ...`
- README.md と食い違うときはコードを正とする。API を変えたら README.md も直す

## 実行環境

```bash
docker compose build                      # requirements.txt を変えたら再ビルド
docker compose run --rm app python xxx.py # スクリプトやテストはコンテナ内で実行
```

- サービス名は `app`。リポジトリを `/work` にマウントし、`PYTHONPATH=/work/src` で読む（ソースの変更に再ビルドは不要）

## テスト

```bash
docker compose run --rm app pytest                # 全部
docker compose run --rm app pytest tests/io       # 一部だけ
```

- `tests/` は `src/grrmlib` と同じ構成（`io/`, `core/`, `operations/`, `pipelines/`, `workflows/`）。コードを変えたら対応するテストも直す
- `tests/data/` には公開してよい実際の計算の入力・出力を置く。大きい log（数十 MB の EQ_list など）は先頭の数件だけに縮めて置き、元のファイルは `tests/data_private/` に置く
- **未公開の入力・出力ファイルは必ず `tests/data_private/` に置く**。リポジトリはパブリックなので、`.gitignore` で除外して GitHub に載せない
  - `tests/data/` には公開してよいファイルだけを置く。未公開のファイルを `tests/data/` に移したり写したりしない
  - `git add -f` で `tests/data_private/` のファイルを追加しない。コミットの前に `git status` で入っていないことを確かめる
  - テストからは fixture で `data_private("gaussian/sp.log")` のように使う。ファイルがない環境（GitHub から clone した人など）ではそのテストを skip する
- 共通の fixture（`data`, `data_private`, `water`, `eqs`, `pts`, `folder_keys`）は `tests/conftest.py`
- `--import-mode=importlib` なので、別のフォルダに同じ名前のテストファイルを置ける。テストどうしで import はできないので、共通の処理は conftest.py に書く
- 既知のバグは `pytest.mark.xfail(strict=True)` で残す。直したら xfail を外す

## 構成（src/grrmlib）

| パッケージ | 内容 |
|---|---|
| `core/` | `Molecule`, `Molecules`（キー → Molecule の dict）, `GroupedMolecules`, `ReactionPathNetwork` と幾何計算（`separate`, `overlay` など） |
| `io/` | ファイルの読み書き。形式ごとに `gaussian.py`, `xyz.py`, `connectable.py`, `grrm/`。`to_dataframe` で polars に変換 |
| `operations/` | core → core の変換。分子の連結（connectable）、電荷・多重度の組み合わせ、EQ の分離（`separate_eqs`）、最安定の選択（`select_lowest`） |
| `pipelines/` | ファイル → ファイルの変換（read → operations → write）。`eq_list_to_xyzs`, `xyzs_to_gaussian_inputs`, `connectables_to_gaussian_sp`, `gaussian_sp_to_grrm_min`, `collect_*` |
| `workflows/` | 計算の手順（SC-AFIR の結果からエネルギー記述子を作る、論文 SI の Figure S8）の 1 段ずつを関数にしたもの。`make_separation_sp`, `select_sp`, `make_opt`, `select_opt`, `make_connection_sp`, `sum_energies`。計算そのものは含まない |

公開する名前は `grrmlib/__init__.py` に集めている。新しい関数を足したらそこにも追加する。

## 設計の方針

- 関数は 5 種類に分ける: ファイル → core（`io` の `read_xxx`）、core → ファイル（`io` の `write_xxx`）、core → core（`core`, `operations`）、ファイル → ファイル（`pipelines`。ファイル → core → ファイルをつないだもの）、計算の手順の 1 段（`workflows`。pipelines と operations をつなぎ、多数のファイルをまとめて扱う）
- pipelines のファイル名は、変換する 2 つの形式の名前を単数形にし、アルファベット順に並べて付ける
  （`eq_list_to_xyzs` → `eq_list_xyz.py`、`xyzs_to_gaussian_inputs` → `gaussian_input_xyz.py`、`collect_gaussian_outputs` → `csv_gaussian_output.py`）
  - `io/` のファイル名（`xyz.py` など）と同じく形式の名前にし、変換の向きは関数名で表す。逆向きの変換も同じファイルに置く
  - ファイル名を関数名と同じにしない（`grrmlib.pipelines.xxx` がファイルではなく関数を指してしまう）
- workflow のような状態を持つクラスは作らない。pipelines と workflows の関数は入出力のパスを全て引数で受け取り、カレントディレクトリに決まった名前で書かない
- workflows（v0.1 の `SeparationWorkflow` などのクラスとは別物）
  - 呼ぶ向きは core → io / operations → pipelines → workflows。下の層から workflows を呼ばない
  - 1 つの関数が、プロジェクトのスクリプト 1 つ分（例：`select_sp`）。プロジェクトに決まったもの（パス、結果ファイルの名前、置換基の名前）は引数で受け取る。
    進み具合や件数は print で表示する
  - 計算の結果は、入力と同じフォルダの xyz（タイトル行の JSON に charge, mult, scfenergy, success, optsteps）を `read_xyzs` で読む。計算そのもの（SP・最適化のループ）は入れない
  - 共通の処理（電荷の選択 `select_charges`、誘導体への展開 `expand_derivatives`、notes の写し `copy_notes`）は operations に置く
- 計算の実行（Gaussian, GRRM, 並列化, ジョブの投入・監視, ジョブ一覧の作成, エラーの後始末）は bash / Slurm に任せ、grrmlib には入れない
- bash / Slurm の実例はユーザーが実際に動かしたものを置く（Claude が雛形を作らない）

## コードの約束事

- **読み書きは関数にする**: `read_xxx(path)` / `write_xxx(mol, path)`、複数ファイルは `read_xxxs(folder, basename)` / `write_xxxs(mols, folder, basename)`。クラスは状態が必要な内部の解析（`io/grrm/min.py` の `_MINOutputParser`）だけに使う
- **キーとフォルダ**: 複数ファイルの読み書きは `io/_files.py` の共通処理を使う
  - `prefixes` なし → `folder/1/2/basename`、あり → `folder/SG=1/SEQ=2/basename`（pipelines はこちら）
  - 読むときは `"None"` → None、整数 → int、それ以外は str。キーが 1 段ならタプルでなく値そのもの
  - 計算結果の読み込みでは `ignore_errors=True` で、読めないファイルを空の `Molecule` にして読み進める
- **xyz のタイトル行**: `write_xyz(s)` の `record=True` で charge, mult, labels を JSON で書く。`read_xyz` はタイトル行が JSON の dict なら属性として読み、そうでなければ `title` にする（scfenergy など charge, mult, labels 以外の値も属性として読む）
- **Molecule の属性**: 名前は cclib に合わせる（`scfenergy`, `atomcoords`, `vibfreqs` など）。基本属性以外は `Molecule._allowed_extra` に登録しないと `TypeError` になる
- **notes**: 原子ごとの追加情報（連結の指定など）。座標の後ろの列として読み書きし、値は int / float。書式は `operations/connectable.py` の `read_connectivity` を参照。Gaussian / GRRM の入力に書き出すときは pipelines では `with_notes=False`
- **単位**: grrmlib のエネルギーは Hartree、座標は Å
- **GRRM の出力**: 多重度は Spin(**2) = S(S+1) から `io/grrm/_spin.py` の `spin2_to_mult` で求める。ヘッセ行列（`hessian`）は下三角だけに値が入り、上三角は NaN（仕様）
- コメント・docstring は日本語。既存のコードのスタイル（4 スペース、引数ごとの改行、`# Functions for reading.` の区切り）に合わせる
- DataFrame は polars を使う（pandas は使わない）

## 既知の問題

- `io/grrm/lup_ts.py` の `read_lup_ts` は動くが、まだ公開していない（`from grrmlib.io.grrm.lup_ts import read_lup_ts`）
- `operations/connectable.py` の `connect_all_and_write_as_gaussian_input` は未定義の `GaussianInputWriter` を使っており動かない
