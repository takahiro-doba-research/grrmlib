import json
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from ..core import Molecule, Molecules
from ._files import build_atomcoords, parse_atomcoords, read_files, write_files


# record=True でタイトル行に JSON で書く属性
_RECORD_ATTRS = [
    "charge",
    "mult",
]


# ----------------------
# Functions for reading.
# ----------------------
def _parse_record(title: str) -> dict | None:
    try:
        record = json.loads(title)
    except json.JSONDecodeError:
        return None
    return record if isinstance(record, dict) else None


def read_xyz(path: str | Path) -> Molecule:
    """
    タイトル行が JSON の dict なら、その値を属性として読む（write_xyz の record=True）。
    そうでなければタイトル行を title として読む。
    """
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    n_atoms = int(lines[0])
    labels, symbols, atomcoords, _ = parse_atomcoords(lines[2:2 + n_atoms])

    record = _parse_record(lines[1])
    if record is None:
        record = {"title": lines[1]}
    # labels がなければ 1 から振る
    record_labels = record.pop("labels", None)
    if record_labels is not None:
        labels = np.array(record_labels)

    return Molecule(
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        **record,
    )


def read_xyzs(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(folder, basename, read_xyz, ignore_errors=ignore_errors)


# ----------------------
# Functions for writing.
# ----------------------
def write_xyz(
    mol: Molecule,
    path: str | Path,
    *,
    title: str = "",
    record: bool = False,
    overwrite: bool = False
) -> None:
    """
    record=True のときは、title の代わりにタイトル行に charge, mult, labels を JSON で書く。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if record:
        values = {a: getattr(mol, a, None) for a in _RECORD_ATTRS}
        # labels を残さないと元の分子（EQ など）の原子と対応がとれない
        values["labels"] = [int(l) for l in mol.labels]
        title = json.dumps(values)

    lines = []
    lines.append(str(len(mol)))
    lines.append(title)
    lines += build_atomcoords(mol, with_notes=False)
    lines.append("")

    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_xyzs(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    prefixes: str | Iterable[str] | None = None,
    title: str = "",
    record: bool = False,
    overwrite: bool = False
) -> None:
    write_files(
        mols,
        folder,
        basename,
        write_xyz,
        prefixes=prefixes,
        title=title,
        record=record,
        overwrite=overwrite
    )
