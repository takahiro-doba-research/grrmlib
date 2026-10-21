import re
from collections.abc import Iterable
from pathlib import Path

from ...core import Molecule, Molecules
from .._files import build_atomcoords, parse_atomcoords, read_files, write_files


# ----------------------
# Functions for reading.
# ----------------------
def _parse_link0(lines: list[str]) -> str | None:
    infile = None

    for line in lines:
        if line.lower().startswith("%infile"):
            infile = re.search(r"=\s*(.+)", line).group(1)

    return infile


def _parse_charge_mult(line: str) -> tuple[int, int]:
    charge, mult = re.search(r"\s*(-?\d+)\s+(\d+)\s*", line).groups()
    return int(charge), int(mult)


def read_grrm_input(path: str | Path) -> Molecule:
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    lines_link0 = [l for l in lines if l.startswith("%")]
    infile = _parse_link0(lines_link0)

    index_route = [i for i, l in enumerate(lines) if l.startswith("#")][0]
    route = lines[index_route]

    line_title = lines[index_route + 1]
    title = line_title if line_title != "" else None

    lines_charge_mult = lines[index_route + 2]
    charge, mult = _parse_charge_mult(lines_charge_mult)

    lines_atomcoords_options = lines[index_route + 3:]
    index_options = [
        i for i, l in enumerate(lines_atomcoords_options)
        if l.startswith("Options")
    ]
    if index_options:
        lines_atomcoords = lines_atomcoords_options[:index_options[0]]
        lines_options = lines_atomcoords_options[index_options[0] + 1:]
        options = [l for l in lines_options if l.strip()]
    else:
        lines_atomcoords = [l for l in lines_atomcoords_options if l.strip()]
        options = None
    labels, symbols, atomcoords, notes = parse_atomcoords(lines_atomcoords)

    return Molecule(
        infile=infile,
        route=route,
        title=title,
        charge=charge,
        mult=mult,
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes,
        options=options
    )


def read_grrm_inputs(
    folder: str | Path,
    basename: str | Path,
    *,
    ignore_errors: bool = False
) -> Molecules:
    return read_files(
        folder,
        basename,
        read_grrm_input,
        ignore_errors=ignore_errors
    )


# ----------------------
# Functions for writing.
# ----------------------
def write_grrm_input(
    mol: Molecule,
    path: str | Path,
    *,
    infile: str | None = None,
    route: str | None = None,
    title: str | None = None,
    options: list[str] | None = None,
    with_notes: bool = True,
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    charge = mol.charge if mol.charge is not None else 0
    mult = mol.mult if mol.mult is not None else 1

    lines = []
    lines += [f"%infile={infile}"] if infile is not None else []
    lines.append(route if route is not None else "#")
    lines.append(title if title is not None else "")
    lines.append(f"{charge} {mult}")
    lines += build_atomcoords(mol, with_notes=with_notes)
    lines += ["Options"] + options if options is not None else []
    lines.append("")

    with path.open("w" if overwrite else "x") as f:
        f.write("\n".join(lines))


def write_grrm_inputs(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    *,
    prefixes: str | Iterable[str] | None = None,
    infile: str | None = None,
    route: str | None = None,
    title: str | None = None,
    options: list[str] | None = None,
    with_notes: bool = True,
    overwrite: bool = False
) -> None:
    write_files(
        mols,
        folder,
        basename,
        write_grrm_input,
        prefixes=prefixes,
        infile=infile,
        route=route,
        title=title,
        options=options,
        with_notes=with_notes,
        overwrite=overwrite
    )
