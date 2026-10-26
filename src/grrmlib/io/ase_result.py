import json
import os
from pathlib import Path

import numpy as np

from ..core import Molecule, Molecules
from ._files import build_atomcoords, read_files


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


# ----------------------
# Functions for reading.
# ----------------------
def read_ase_result(path: str | Path) -> Molecule:
    lines = Path(path).read_text().splitlines()
    n_atoms = int(lines[0])
    record = json.loads(lines[1])
    # 古い結果ファイルには key も入っているが、key はフォルダから決める
    record.pop("key", None)
    # labels を保存していない古い結果ファイルは 1 から振り直す
    labels = record.pop("labels", None)
    labels = np.arange(1, n_atoms + 1) if labels is None else np.array(labels)

    lines_coord = lines[2:2 + n_atoms]
    return Molecule(
        labels=labels,
        symbols=[l.split()[0] for l in lines_coord],
        atomcoords=np.array([list(map(float, l.split()[1:4])) for l in lines_coord]),
        **record,
    )


def read_ase_results(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(
        folder,
        basename,
        read_ase_result,
        ignore_errors=ignore_errors
    )


# ----------------------
# Functions for writing.
# ----------------------
def write_ase_result(
    mol: Molecule,
    path: str | Path,
    *,
    overwrite: bool = False
) -> None:
    """
    xyz の形式で、タイトル行に charge, mult, scfenergy などを JSON で書く。
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists")

    record = {a: getattr(mol, a, None) for a in _RECORD_ATTRS}
    # labels と notes も残さないと読み直したときに元の分子と対応がとれない
    record["labels"] = [int(l) for l in mol.labels]
    # notes には回転角などの float も入るので、numpy の数値を Python の数値に直すだけにする
    record["notes"] = (
        [[np.asarray(n).item() for n in note] for note in mol.notes]
        if mol.notes is not None
        else None
    )

    lines = [str(len(mol)), json.dumps(record)]
    lines += build_atomcoords(mol, with_notes=False)
    lines.append("")

    # 途中で止まっても壊れたファイルが残らないように rename で書き込む
    path.parent.mkdir(parents=True, exist_ok=True)
    path_tmp = path.with_name(path.name + ".tmp")
    path_tmp.write_text("\n".join(lines))
    os.replace(path_tmp, path)
