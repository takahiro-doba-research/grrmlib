import re
from collections.abc import Iterable
from pathlib import Path

import cclib
import numpy as np
from periodictable import elements

from ..core import Molecule, Molecules
from ._files import build_atomcoords, parse_atomcoords, read_files, write_files


# ----------------------
# Functions for reading.
# ----------------------
def _parse_charge_mult(line: str) -> tuple[int, int]:
    charge, mult = re.search(r"\s*(-?\d+)\s+(\d+)\s*", line).groups()
    return int(charge), int(mult)


def read_gaussian_input(path: str | Path) -> Molecule:
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    # parse lines
    indices_blank = [i for i, line in enumerate(lines) if line.strip() == ""]

    lines_link0_route = lines[:indices_blank[0]]
    link0 = [l for l in lines_link0_route if l.startswith("%")]

    index_route = [i for i, l in enumerate(lines_link0_route) if l.startswith("#")][0]
    route = lines_link0_route[index_route:]

    title = lines[indices_blank[0] + 1:indices_blank[1]]

    lines_charge_mult = lines[indices_blank[1] + 1]
    charge, mult = _parse_charge_mult(lines_charge_mult)

    lines_atomcoords = lines[indices_blank[1] + 2:indices_blank[2]]
    labels, symbols, atomcoords, notes = parse_atomcoords(lines_atomcoords)

    return Molecule(
        link0=link0,
        route=route,
        title=title,
        charge=charge,
        mult=mult,
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes
    )


def read_gaussian_inputs(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(
        folder,
        basename,
        read_gaussian_input,
        ignore_errors=ignore_errors
    )


def read_gaussian_output(path: str | Path) -> Molecule:
    path = Path(path)
    parser = cclib.io.ccopen(path)
    data = parser.parse()
    success = data.metadata["success"]

    if success:
        symbols = [str(elements[n]) for n in data.atomnos]
        labels = np.arange(1, len(symbols) + 1)
        return Molecule(
            charge=data.charge,
            mult=data.mult,
            labels=labels,
            symbols=symbols,
            atomcoords=data.atomcoords[-1],
            scfenergy=data.scfenergies[-1],
            success=success,
        )

    else:
        return Molecule(success=success)


def read_gaussian_outputs(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(
        folder,
        basename,
        read_gaussian_output,
        ignore_errors=ignore_errors
    )


# ----------------------
# Functions for writing.
# ----------------------
def _build_charge_mult(mol: Molecule) -> list[str]:
    charge = mol.charge if mol.charge is not None else 0
    mult = mol.mult if mol.mult is not None else 1
    return [f"{charge} {mult}"]


def write_gaussian_input(
    mol: Molecule,
    path: str | Path,
    *,
    link0: list[str] | None = None,
    route: list[str] | None = None,
    title: list[str] | None = None,
    with_notes: bool = True,
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines += link0 if link0 is not None else []
    lines += route if route is not None else ["#"]
    lines.append("")
    lines += title if title is not None else ["None"]
    lines.append("")
    lines += _build_charge_mult(mol)
    lines += build_atomcoords(mol, with_notes=with_notes)
    lines.append("")
    lines.append("")

    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_gaussian_inputs(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    prefixes: str | Iterable[str] | None = None,
    link0: list[str] | None = None,
    route: list[str] | None = None,
    title: list[str] | None = None,
    with_notes: bool = True,
    overwrite: bool = False
) -> None:
    write_files(
        mols,
        folder,
        basename,
        write_gaussian_input,
        prefixes=prefixes,
        link0=link0,
        route=route,
        title=title,
        with_notes=with_notes,
        overwrite=overwrite
    )


def write_gaussian_scan(
    mols: Molecules,
    path: str | Path,
    *,
    overwrite: bool = False
) -> None:
    """
    GaussView で連続した構造として見るための疑似的な scan の出力ファイルを書く。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    num = len(mols)
    lines = [" #p\n", " \n"]

    for i, mol in enumerate(mols.values(), start=1):
        lines += [
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            "                          Input orientation:                          \n",
            " ---------------------------------------------------------------------\n",
            " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
            " Number     Number       Type             X           Y           Z   \n",
            " ---------------------------------------------------------------------\n",
            " ---------------------------------------------------------------------\n",
            " \n",
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            f" Step number   1 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
            " \n",
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            "                          Input orientation:                          \n",
            " ---------------------------------------------------------------------\n",
            " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
            " Number     Number       Type             X           Y           Z   \n",
            " ---------------------------------------------------------------------\n",
            *[
                f"{j:7d} {elements.symbol(s).number:10d}           0     {x:11.6f} {y:11.6f} {z:11.6f}\n"
                for j, (s, (x, y, z), _) in enumerate(mol.iter_atoms(), start=1)
            ],
            " ---------------------------------------------------------------------\n",
            f" SCF Done:  E = {mol.scfenergy:15.12f}     A.U.\n",
            " \n",
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            f" Step number   2 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
            " \n",
        ]

    lines += [
        " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
        " Normal termination of Gaussian 16\n"
    ]

    with path.open("w" if overwrite else "x") as f:
        f.writelines(lines)
