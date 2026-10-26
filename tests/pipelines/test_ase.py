import pytest
from ase.calculators.emt import EMT

import grrmlib as gl


@pytest.fixture
def path_com(tmp_path, water):
    path = tmp_path / "0/gaussian_sp.com"
    gl.write_gaussian_input(water, path)
    return path


def test_run_ase_opt_then_sp(tmp_path, path_com):
    opt = gl.run_ase(path_com, tmp_path / "0/opt.xyz", "opt", EMT(), fmax=0.05)
    assert opt.success and opt.status == "Converged"

    # opt の結果ファイルから続けて sp を流せる
    sp = gl.run_ase(tmp_path / "0/opt.xyz", tmp_path / "0/sp.xyz", "sp", EMT())
    assert sp.status == "SP"
    assert sp.scfenergy == pytest.approx(opt.scfenergy)
    assert (sp.charge, sp.mult) == (0, 1)
    assert list(gl.read_ase_results(tmp_path, "sp.xyz")) == [0]


def test_run_ase_resumes(tmp_path, path_com):
    path_out = tmp_path / "0/sp.xyz"
    sp = gl.run_ase(path_com, path_out, "sp", EMT())
    mtime = path_out.stat().st_mtime_ns

    # 結果があれば計算せずに読む
    sp_again = gl.run_ase(path_com, path_out, "sp", EMT())
    assert path_out.stat().st_mtime_ns == mtime
    assert sp_again.scfenergy == pytest.approx(sp.scfenergy)


def test_run_ase_records_errors(tmp_path, path_com):
    mol = gl.run_ase(path_com, tmp_path / "0/xxx.xyz", "xxx", EMT())

    assert mol.success is False
    assert mol.status == "Error: ValueError: Unknown task: xxx"
    assert gl.read_ase_result(tmp_path / "0/xxx.xyz").success is False


def test_run_ase_from_xyz(tmp_path, water):
    # write_xyz の record=True で書いた xyz を入力にできる
    water.charge, water.mult = 1, 2
    gl.write_xyz(water, tmp_path / "0/mol.xyz", record=True)
    sp = gl.run_ase(tmp_path / "0/mol.xyz", tmp_path / "0/sp.xyz", "sp", EMT())

    assert sp.status == "SP"
    assert (sp.charge, sp.mult) == (1, 2)
