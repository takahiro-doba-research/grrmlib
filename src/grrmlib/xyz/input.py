from pathlib import Path

import numpy as np

from ..core import Molecule, Molecules


# ----------------------
# Functions for reading.
# ----------------------
def _parse_atomcoords(lines: list[str]) -> tuple[np.ndarray, list[str], np.ndarray]:
    labels = np.arange(1, len(lines) + 1)
    symbols = [l.split()[0] for l in lines]
    atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
    return labels, symbols, atomcoords


def read_xyz(path: str | Path) -> Molecule:
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()
    
    n_atoms = int(lines[0])
    title = lines[1]
    labels, symbols, atomcoords = _parse_atomcoords(lines[2:2 + n_atoms])
    
    return Molecule(
        title=title,
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
    )


def read_xyzs(
    folder: str | Path,
    basename: str | Path
) -> Molecules:
    folder = Path(folder)
    basename = str(basename)
    mols = Molecules()
    
    for path in sorted(folder.rglob(basename)):
        key = tuple(map(int, path.parent.relative_to(folder).parts))
        key = key[0] if len(key) == 1 else key
        mols[key] = read_xyz(path)
    
    return mols


# ----------------------
# Functions for writing.
# ----------------------
def _build_atomcoords(mol: Molecule) -> list[str]:
    lines = []
    for s, (x, y, z), _ in mol.iter_atoms():
        lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}")
    return lines


def write_xyz(
    mol: Molecule,
    path: str | Path,
    *,
    title: str = "",
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    lines = []
    lines.append(str(len(mol)))
    lines.append(title)
    lines += _build_atomcoords(mol)
    lines.append("")
    
    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_xyzs(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    title: str = "",
    overwrite: bool = False
) -> None:
    folder = Path(folder)
    basename = Path(basename)
    
    for key, mol in mols.items():
        key = (key, ) if isinstance(key, int) else key
        path = folder.joinpath(*map(str, key)) / basename
        
        write_xyz(
            mol,
            path,
            title=title,
            overwrite=overwrite
        )