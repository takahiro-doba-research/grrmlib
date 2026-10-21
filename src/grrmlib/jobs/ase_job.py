from __future__ import annotations

import json
import os
import time
from collections.abc import Callable, Hashable
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from ..core import Molecule, Molecules
from ..calculators import single_point, optimize
from ..io import key_to_folder

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator


# 結果ファイルの xyz のタイトル行に JSON で保存する属性
_RECORD_ATTRS = [
    "charge",
    "mult",
    "scfenergy",
    "success",
    "status",
    "optsteps",
    "elapsed",
]

# ワーカープロセスごとに 1 つだけ保持する calculator
_calc: Calculator | None = None


def write_result(mol: Molecule, key: Hashable, path: Path) -> None:
    record = {"key": list(key) if isinstance(key, tuple) else key}
    record |= {a: getattr(mol, a, None) for a in _RECORD_ATTRS}
    # labels と notes も残さないと読み直したときに元の分子と対応がとれない
    record["labels"] = [int(l) for l in mol.labels]
    # notes には回転角などの float も入るので、numpy の数値を Python の数値に直すだけにする
    record["notes"] = (
        [[np.asarray(n).item() for n in note] for note in mol.notes]
        if mol.notes is not None
        else None
    )

    lines = [str(len(mol)), json.dumps(record)]
    for s, (x, y, z), _ in mol.iter_atoms():
        lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}")
    lines.append("")

    # 途中で止まっても壊れたファイルが残らないように rename で書き込む
    path.parent.mkdir(parents=True, exist_ok=True)
    path_tmp = path.with_name(path.name + ".tmp")
    path_tmp.write_text("\n".join(lines))
    os.replace(path_tmp, path)


def read_result(path: str | Path) -> tuple[Hashable, Molecule]:
    lines = Path(path).read_text().splitlines()
    n_atoms = int(lines[0])
    record = json.loads(lines[1])
    key = record.pop("key")
    key = tuple(key) if isinstance(key, list) else key
    # labels を保存していない古い結果ファイルは 1 から振り直す
    labels = record.pop("labels", None)
    labels = np.arange(1, n_atoms + 1) if labels is None else np.array(labels)

    lines_coord = lines[2:2 + n_atoms]
    mol = Molecule(
        labels=labels,
        symbols=[l.split()[0] for l in lines_coord],
        atomcoords=np.array([list(map(float, l.split()[1:4])) for l in lines_coord]),
        **record,
    )
    return key, mol


def read_results(folder: str | Path, basename: str) -> Molecules:
    mols = Molecules()
    for path in sorted(Path(folder).rglob(basename)):
        key, mol = read_result(path)
        mols[key] = mol
    return mols


def _init_worker(calc_factory: Callable[[], Calculator], nthreads: int) -> None:
    global _calc
    import torch
    torch.set_num_threads(nthreads)
    _calc = calc_factory()


def _run_one(
    key: Hashable,
    mol: Molecule,
    path: Path,
    task: str,
    opt_kwargs: dict[str, Any]
) -> tuple[Hashable, bool]:
    t = time.perf_counter()

    try:
        if task == "sp":
            mol_new = single_point(mol, _calc)
        elif task == "opt":
            mol_new = optimize(mol, _calc, **opt_kwargs)
        else:
            raise ValueError(f"Unknown task: {task}")
    except Exception as e:
        mol_new = mol.copy()
        mol_new.scfenergy = None
        mol_new.success = False
        mol_new.status = f"Error: {type(e).__name__}: {e}"
        mol_new.elapsed = time.perf_counter() - t

    # labels と notes は計算前の分子のものを残す
    mol_new.labels = mol.labels
    write_result(mol_new, key, path)
    return key, mol_new.success


class ASEJobRunner:
    """
    ASE calculator で SP / 構造最適化をまとめて実行する。
    結果は 1 件ずつ folder/key/basename に保存し、既に結果があるものは飛ばす
    （途中で止めても同じコマンドで再開できる）。
    """

    def __init__(
        self,
        calc_factory: Callable[[], Calculator],
        *,
        nprocs: int = 1,
        nthreads: int = 1,
    ) -> None:
        self.calc_factory = calc_factory
        self.nprocs = nprocs
        self.nthreads = nthreads

    def run(
        self,
        mols: Molecules,
        folder: str | Path,
        task: str,
        *,
        basename: str | None = None,
        **opt_kwargs: Any
    ) -> Molecules:
        basename = basename or f"{task}.xyz"
        jobs = [
            (key, mol, key_to_folder(folder, key) / basename)
            for key, mol in mols.items()
        ]
        jobs_todo = [job for job in jobs if not job[2].exists()]
        # 大きい分子から先に流して並列の待ち時間を減らす
        jobs_todo.sort(key=lambda job: len(job[1]), reverse=True)
        print(f"{task}: {len(jobs) - len(jobs_todo)} done, {len(jobs_todo)} to run", flush=True)

        if jobs_todo:
            self._run_jobs(jobs_todo, task, opt_kwargs)

        results = Molecules()
        for key, _, path in jobs:
            results[key] = read_result(path)[1]
        return results

    def _run_jobs(self, jobs, task, opt_kwargs) -> None:
        t = time.perf_counter()
        n_failed = 0

        with ProcessPoolExecutor(
            max_workers=self.nprocs,
            mp_context=get_context("spawn"),
            initializer=_init_worker,
            initargs=(self.calc_factory, self.nthreads),
        ) as executor:
            futures = [
                executor.submit(_run_one, key, mol, path, task, opt_kwargs)
                for key, mol, path in jobs
            ]

            for i, future in enumerate(as_completed(futures), 1):
                key, success = future.result()
                n_failed += not success

                if i % 10 == 0 or i == len(jobs):
                    print(
                        f"{task}: {i}/{len(jobs)} finished "
                        f"({n_failed} failed, {time.perf_counter() - t:.0f} s)",
                        flush=True
                    )
