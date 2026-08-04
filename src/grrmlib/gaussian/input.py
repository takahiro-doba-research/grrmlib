import re
from pathlib import Path

import numpy as np

from ..core import Molecule, Molecules


# ----------------------
# Functions for reading.
# ----------------------
def _parse_charge_mult(line: str) -> tuple[int, int]:
    charge, mult = re.search(r"\s*(-?\d+)\s+(\d+)\s*", line).groups()
    return int(charge), int(mult)


def _parse_atomcoords(lines: list[str]) -> tuple[
    np.ndarray,
    list[str],
    np.ndarray,
    list[list[str]] | None
]:
    labels = np.arange(1, len(lines) + 1)
    symbols = [l.split()[0] for l in lines]
    atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
    notes = [l.split()[4:] for l in lines]
    
    if all(not note for note in notes):
        notes = None
    
    return labels, symbols, atomcoords, notes


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
    labels, symbols, atomcoords, notes = _parse_atomcoords(lines_atomcoords)
    
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
    basename: str | Path
) -> Molecules:
    folder = Path(folder)
    basename = str(basename)
    mols = Molecules()
    
    for path in sorted(folder.rglob(basename)):
        key = tuple(map(int, path.parent.relative_to(folder).parts))
        key = key[0] if len(key) == 1 else key
        mols[key] = read_gaussian_input(path)
    
    return mols


# ----------------------
# Functions for writing.
# ----------------------
def _build_charge_mult(mol: Molecule) -> list[str]:
    charge = mol.charge if mol.charge is not None else 0
    mult = mol.mult if mol.mult is not None else 1
    return [f"{charge} {mult}"]


def _build_atomcoords(mol: Molecule) -> list[str]:
    lines = []
    for s, (x, y, z), n in mol.iter_atoms():
        if n is None:
            lines.append(
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
            )
        else:
            lines.append(
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
                f" {' '.join(map(str, n))}"
            )
    return lines


def write_gaussian_input(
    mol: Molecule,
    path: str | Path,
    *,
    link0: list[str] | None = None,
    route: list[str] | None = None,
    title: list[str] | None = None,
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
    lines += _build_atomcoords(mol)
    lines.append("")
    lines.append("")
    
    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_gaussian_inputs(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    link0: list[str] | None = None,
    route: list[str] | None = None,
    title: list[str] | None = None,
    overwrite: bool = False
) -> None:
    folder = Path(folder)
    basename = Path(basename)
    
    for key, mol in mols.items():
        key = (key, ) if isinstance(key, int) else key
        path = folder.joinpath(*map(str, key)) / basename
        
        write_gaussian_input(
            mol,
            path,
            link0=link0,
            route=route,
            title=title,
            overwrite=overwrite
        )