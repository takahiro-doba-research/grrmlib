from __future__ import annotations

import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ..core import Molecule
from ..io import read_ase_result, read_gaussian_input, write_ase_result

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator


def run_ase(
    path_in: str | Path,
    path_out: str | Path,
    task: str,
    calc: Calculator,
    **opt_kwargs: Any
) -> Molecule:
    """
    path_in の分子を ASE calculator で計算し、結果を path_out に書いて返す。

    - task は "sp" または "opt"。opt_kwargs は optimize に渡す（fmax, steps, optimizer）
    - path_in は Gaussian の入力（.com / .gjf）か、タイトル行が JSON の xyz
      （write_xyz の record=True で書いたもの、または run_ase の結果ファイル。opt → sp など）。
      charge と mult が必要
    - path_out が既にあれば計算せずにそれを読んで返す（止まっても同じコマンドで再開できる）
    - 計算中のエラーは success=False と status に記録する

    並列化は bash / Slurm に任せる。calc の読み込み（UMA など）は時間がかかるので、
    1 つのプロセスで複数のファイルを流すときは同じ calc を使い回す。
    """
    # ase / fairchem は重いので必要になったときだけ import する
    from ..calculators import optimize, single_point

    path_in = Path(path_in)
    path_out = Path(path_out)

    if path_out.exists():
        return read_ase_result(path_out)

    if path_in.suffix in (".com", ".gjf"):
        mol = read_gaussian_input(path_in)
    else:
        mol = read_ase_result(path_in)

    t = time.perf_counter()

    try:
        if task == "sp":
            mol_new = single_point(mol, calc)
        elif task == "opt":
            mol_new = optimize(mol, calc, **opt_kwargs)
        else:
            raise ValueError(f"Unknown task: {task}")
    except Exception as e:
        mol_new = mol.copy()
        mol_new.scfenergy = None
        mol_new.success = False
        mol_new.status = f"Error: {type(e).__name__}: {e}"
        mol_new.elapsed = time.perf_counter() - t

    write_ase_result(mol_new, path_out)
    return mol_new
