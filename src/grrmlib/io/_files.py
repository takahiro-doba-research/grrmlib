from collections.abc import Callable, Hashable, Iterable
from pathlib import Path
from typing import Any

import numpy as np

from ..core import Molecule, Molecules


def key_to_folder(
    folder: str | Path,
    key: Hashable,
    prefixes: str | Iterable[str] | None = None
) -> Path:
    """
    key の各要素を 1 階層ずつのフォルダにする。
    prefixes を渡すと "prefix=value" の形にする（例: SG=1/SEQ=2）。
    """
    key = key if isinstance(key, tuple) else (key, )

    if prefixes is None:
        parts = [str(k) for k in key]
    else:
        prefixes = (prefixes, ) if isinstance(prefixes, str) else tuple(prefixes)
        if len(prefixes) != len(key):
            raise ValueError("prefixes and key must have the same length")
        parts = [f"{p}={k}" for p, k in zip(prefixes, key)]

    return Path(folder).joinpath(*parts)


def _parse_part(part: str) -> Hashable:
    value = part.split("=", 1)[-1]

    if value == "None":
        return None

    try:
        return int(value)
    except ValueError:
        return value


def folder_to_key(folder: str | Path, path: str | Path) -> Hashable:
    """
    key_to_folder の逆。"prefix=value" の形なら value の部分を使う。
    """
    parts = Path(path).parent.relative_to(folder).parts
    key = tuple(_parse_part(p) for p in parts)
    return key[0] if len(key) == 1 else key


def read_files(
    folder: str | Path,
    basename: str | Path,
    read: Callable[[Path], Molecule],
    *,
    ignore_errors: bool = False
) -> Molecules:
    folder = Path(folder)
    mols = Molecules()

    for path in sorted(folder.rglob(str(basename))):
        key = folder_to_key(folder, path)

        try:
            mols[key] = read(path)
        except Exception as e:
            if not ignore_errors:
                raise
            # 計算が途中で落ちたファイルなどは空の Molecule にして読み進める
            print(f"{path}: {type(e).__name__}: {e}")
            mols[key] = Molecule()

    return mols


def write_files(
    mols: Molecules,
    folder: str | Path,
    basename: str | Path,
    write: Callable[..., None],
    *,
    prefixes: str | Iterable[str] | None = None,
    **kwargs: Any
) -> None:
    for key, mol in mols.items():
        path = key_to_folder(folder, key, prefixes) / basename
        write(mol, path, **kwargs)


def parse_note(tokens: list[str]) -> list[int | float]:
    return [int(t) if t.lstrip("-").isdigit() else float(t) for t in tokens]


def build_atomcoords(mol: Molecule, *, with_notes: bool = True) -> list[str]:
    lines = []
    for s, (x, y, z), n in mol.iter_atoms():
        line = f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
        if with_notes and n is not None:
            line += f" {' '.join(map(str, n))}"
        lines.append(line)
    return lines


def parse_atomcoords(lines: list[str]) -> tuple[
    np.ndarray,
    list[str],
    np.ndarray,
    list[list[int | float]] | None
]:
    labels = np.arange(1, len(lines) + 1)
    symbols = [l.split()[0] for l in lines]
    atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
    notes = [parse_note(l.split()[4:]) for l in lines]

    if all(not note for note in notes):
        notes = None

    return labels, symbols, atomcoords, notes
