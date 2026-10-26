import json

import pytest

import grrmlib as gl


def test_write_and_read_ase_result(tmp_path, water):
    water.labels = water.labels + 10
    water.notes = [[1, 2, 3, 120.0], [1], [2]]
    water.scfenergy = -1.0
    gl.write_ase_result(water, tmp_path / "sp.xyz")
    mol = gl.read_ase_result(tmp_path / "sp.xyz")

    assert list(mol.labels) == [11, 12, 13]
    assert mol.notes == [[1, 2, 3, 120.0], [1], [2]]
    assert mol.scfenergy == -1.0
    assert (mol.charge, mol.mult) == (0, 1)
    assert mol.atomcoords == pytest.approx(water.atomcoords)
    assert not (tmp_path / "sp.xyz.tmp").exists()


def test_write_ase_result_overwrite(tmp_path, water):
    gl.write_ase_result(water, tmp_path / "sp.xyz")
    with pytest.raises(FileExistsError):
        gl.write_ase_result(water, tmp_path / "sp.xyz")
    gl.write_ase_result(water, tmp_path / "sp.xyz", overwrite=True)


def test_read_ase_result_old_format(tmp_path, water):
    # 古い結果ファイルは key を持ち、labels がない
    gl.write_ase_result(water, tmp_path / "sp.xyz")
    lines = (tmp_path / "sp.xyz").read_text().splitlines()
    record = json.loads(lines[1])
    del record["labels"]
    lines[1] = json.dumps({"key": [0, 1]} | record)
    (tmp_path / "sp.xyz").write_text("\n".join(lines))

    mol = gl.read_ase_result(tmp_path / "sp.xyz")
    assert list(mol.labels) == [1, 2, 3]


def test_read_ase_results(tmp_path, water):
    gl.write_ase_result(water, tmp_path / "0/1/sp.xyz")
    gl.write_ase_result(water, tmp_path / "1/0/sp.xyz")
    (tmp_path / "1/1").mkdir(parents=True)
    (tmp_path / "1/1/sp.xyz").write_text("broken\n")

    mols = gl.read_ase_results(tmp_path, "sp.xyz", ignore_errors=True)
    assert list(mols) == [(0, 1), (1, 0), (1, 1)]
    # 読めないファイルは空の Molecule
    assert mols[(1, 1)].atomcoords is None
