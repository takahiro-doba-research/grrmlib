"""
connectable ファイル: Gaussian 入力と同じ形式で、座標の後ろに結合の情報 (notes) を書いたもの。
notes の書式は operations.connectable.read_connectivity を参照。
"""
import re
from collections.abc import Iterable
from pathlib import Path

from ..core import Molecule, Molecules
from ._files import build_atomcoords, parse_atomcoords, read_files, write_files


# ----------------------
# Functions for reading.
# ----------------------
def _parse_charge_mult(line: str) -> tuple[int, int]:
    charge, mult = re.search(r"\s*(-?\d+)\s+(\d+)\s*", line).groups()
    return int(charge), int(mult)


def read_connectable(path: str | Path) -> Molecule:
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    indices_blank = [i for i, line in enumerate(lines) if line.strip() == ""]

    title = lines[indices_blank[0] + 1:indices_blank[1]]

    lines_charge_mult = lines[indices_blank[1] + 1]
    charge, mult = _parse_charge_mult(lines_charge_mult)

    lines_atomcoords = lines[indices_blank[1] + 2:indices_blank[2]]
    labels, symbols, atomcoords, notes = parse_atomcoords(lines_atomcoords)

    return Molecule(
        title=title,
        charge=charge,
        mult=mult,
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes,
    )


def read_connectables(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(
        folder,
        basename,
        read_connectable,
        ignore_errors=ignore_errors
    )


# ----------------------
# Functions for writing.
# ----------------------
def write_connectable(
    mol: Molecule,
    path: str | Path,
    *,
    title: list[str] | None = None,
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    charge = mol.charge if mol.charge is not None else 0
    mult = mol.mult if mol.mult is not None else 1

    lines = []
    lines.append("#")
    lines.append("")
    lines += title if title is not None else ["None"]
    lines.append("")
    lines.append(f"{charge} {mult}")
    lines += build_atomcoords(mol, with_notes=True)
    lines.append("")
    lines.append("")

    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_connectables(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    prefixes: str | Iterable[str] | None = None,
    title: list[str] | None = None,
    overwrite: bool = False
) -> None:
    write_files(
        mols,
        folder,
        basename,
        write_connectable,
        prefixes=prefixes,
        title=title,
        overwrite=overwrite
    )
