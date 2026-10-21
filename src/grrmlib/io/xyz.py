from collections.abc import Iterable
from pathlib import Path

from ..core import Molecule, Molecules
from ._files import build_atomcoords, parse_atomcoords, read_files, write_files


# ----------------------
# Functions for reading.
# ----------------------
def read_xyz(path: str | Path) -> Molecule:
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    n_atoms = int(lines[0])
    title = lines[1]
    labels, symbols, atomcoords, _ = parse_atomcoords(lines[2:2 + n_atoms])

    return Molecule(
        title=title,
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
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
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

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
    overwrite: bool = False
) -> None:
    write_files(
        mols,
        folder,
        basename,
        write_xyz,
        prefixes=prefixes,
        title=title,
        overwrite=overwrite
    )
