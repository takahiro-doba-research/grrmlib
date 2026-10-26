import json
import shutil
from pathlib import Path

import polars as pl
import pytest

import grrmlib as gl


def write_results(folder: Path, basename_in: str, basename_out: str, optsteps: int | None = None) -> None:
    # 計算の代わりに、タイトル行の JSON に scfenergy と success を書いた結果ファイルを置く
    # （エネルギーは原子数と電荷で決め、電荷が 0 のものを最も低くする）
    for path in folder.rglob(basename_in):
        lines = path.read_text().splitlines()
        record = json.loads(lines[1])
        record.update(
            scfenergy=-float(len(lines) - 2) + 0.1 * abs(record["charge"]),
            success=True,
            optsteps=optsteps,
        )
        lines[1] = json.dumps(record)
        path.with_name(basename_out).write_text("\n".join(lines) + "\n")


def test_separation(tmp_path, data):
    # 分離 → SP の結果 → 選択 → 最適化の結果 → 電荷の選択
    folder = tmp_path / "separation"
    table = gl.make_separation_sp(
        data / "grrm/td_c1_193_SCAFIR6_EQ_list.log",
        folder / "sp",
        folder / "separation.csv",
        charges=[-1, 0, 1],
        mults=1,
    )
    write_results(folder / "sp", "init.xyz", "result.xyz")
    gl.select_sp(
        folder / "sp", "result.xyz", folder / "sp.csv", folder / "sp_selected.csv",
        key_names=gl.SEPARATION_KEYS, group_by=["SG", "charge", "mult"],
    )
    gl.make_opt(folder / "sp_selected.csv", folder / "sp", folder / "opt", key_names=gl.SEPARATION_KEYS)
    write_results(folder / "opt", "init.xyz", "result.xyz", optsteps=1)
    selected = gl.select_opt(
        folder / "opt", "result.xyz", folder / "separation.csv", folder / "opt.csv", folder / "opt_selected.csv",
        key_names=gl.SEPARATION_KEYS,
    )

    # EQ 0〜2 はどれも 3 つの中性の分子に分かれる
    assert selected.columns == ["G", "fragment", *gl.SEPARATION_KEYS, "scfenergy"]
    assert selected.height == 3 * table["G"].n_unique()
    assert set(selected["charge"]) == {0}
    assert pl.read_csv(folder / "opt.csv").columns == [*gl.SEPARATION_KEYS, "scfenergy", "success", "optsteps"]


@pytest.fixture
def connection(tmp_path, data):
    # SG 0 の最適化構造として core.com を置き、置換基 sub（3 つの角度）をつなぐ
    core = gl.read_connectable(data / "connectable/core.com")
    core.charge, core.mult = 0, 1
    gl.write_xyz(core, tmp_path / "opt/SG=0/SEQ=0/charge=0/mult=1/result.xyz", record=True)
    (tmp_path / "sub").mkdir()
    shutil.copy(data / "connectable/sub.com", tmp_path / "sub/0.com")
    pl.DataFrame({"SG": [0], "SEQ": [0], "charge": [0], "mult": [1]}).write_csv(tmp_path / "selected.csv")
    pl.DataFrame({
        "G": [0], "EQ": [0], "SEQ": [0], "SG": [0], "fragment": [0], "labels": ["1_2_3_4_5"],
    }).write_csv(tmp_path / "separation.csv")
    return tmp_path


def test_connection(connection, data):
    folder = connection / "connection"
    key_names = gl.connection_key_names(["sub"])
    n = gl.make_connection_sp(
        connection / "selected.csv", connection / "opt", "result.xyz", data / "connectable/core.com",
        [connection / "sub"], folder / "sp",
    )
    assert n == 3
    assert key_names == [*gl.SEPARATION_KEYS, "sub", "sub_angle"]

    write_results(folder / "sp", "init.xyz", "result.xyz")
    selected = gl.select_sp(
        folder / "sp", "result.xyz", folder / "sp.csv", folder / "sp_selected.csv",
        key_names=key_names, group_by=[*gl.SEPARATION_KEYS, "sub"],
    )
    # 角度の違う 3 つの配座から 1 つを選ぶ
    assert selected.height == 1

    gl.make_opt(folder / "sp_selected.csv", folder / "sp", folder / "opt", key_names=key_names)
    write_results(folder / "opt", "init.xyz", "result.xyz", optsteps=1)
    gl.select_opt(
        folder / "opt", "result.xyz", connection / "separation.csv", folder / "opt.csv", folder / "opt_selected.csv",
        key_names=key_names, substituents={"sub": [0]},
    )
    table = gl.sum_energies(
        folder / "opt_selected.csv", connection / "energy_long.csv", connection / "energy.csv", substituents=["sub"],
    )
    assert table.columns == ["sub", "G=0"]
    assert table["G=0"].to_list() == [-5.0]
