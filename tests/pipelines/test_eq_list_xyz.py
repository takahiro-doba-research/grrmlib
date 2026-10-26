import polars as pl
import pytest

import grrmlib as gl


@pytest.fixture
def paths(tmp_path, data):
    return {
        "eq_list": data / "grrm/td_c1_193_SCAFIR6_EQ_list.log",
        "folder": tmp_path / "separation",
        "table": tmp_path / "separation.csv",
    }


def _eq_list_to_xyzs(paths, **kwargs):
    return gl.eq_list_to_xyzs(
        paths["eq_list"],
        paths["folder"],
        paths["table"],
        **kwargs
    )


def test_eq_list_to_xyzs(paths, folder_keys):
    _eq_list_to_xyzs(paths, charges=[0, 1], mults=[1, 2])

    # SG ごとの SEQ は tests/operations/test_separation.py を参照。
    # それぞれ (charge, mult) = (0, 1), (1, 2) の 2 通り
    seqs = [(0, 0), (0, 6), (1, 1), (1, 4), (1, 7), (2, 2), (2, 5), (2, 8), (3, 3)]
    assert folder_keys(paths["folder"], "seq.xyz") == sorted(
        (sg, seq, charge, mult) for sg, seq in seqs for charge, mult in [(0, 1), (1, 2)]
    )

    # タイトル行の JSON から charge, mult と元の EQ のラベルを読める
    mol = gl.read_xyz(paths["folder"] / "SG=1/SEQ=4/charge=1/mult=2/seq.xyz")
    assert (mol.charge, mol.mult) == (1, 2)
    assert len(mol) == 12
    table = pl.read_csv(paths["table"])
    labels = table.filter(pl.col("SEQ") == 4)["labels"].item()
    assert "_".join(map(str, mol.labels)) == labels


def test_eq_list_to_xyzs_table(paths):
    table = _eq_list_to_xyzs(paths, charges=0, mults=1)

    assert pl.read_csv(paths["table"]).equals(table)
    assert table.columns == ["G", "EQ", "SEQ", "SG", "fragment", "labels"]
    assert table.height == 9


def test_eq_list_to_xyzs_overwrite(paths):
    _eq_list_to_xyzs(paths, charges=0, mults=1)
    with pytest.raises(FileExistsError):
        _eq_list_to_xyzs(paths, charges=0, mults=1)
    _eq_list_to_xyzs(paths, charges=0, mults=1, overwrite=True)
