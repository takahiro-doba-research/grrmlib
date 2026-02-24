import re
from pathlib import Path

import numpy as np

from ..core import Molecule


class ConnectableReader:
    
    def read(self, path: str | Path) -> Molecule:
        path = Path(path)
        text = path.read_text()
        lines = text.splitlines(keepends=True)
        return self.parse(lines)
    
    def _parse_charge_mult(self, lines: str) -> tuple[int, int]:
        charge, mult = re.search(r"\s*(-?\d+)\s+(\d+)\s*", lines).groups()
        return int(charge), int(mult)
    
    def _parse_atomcoords(
        self,
        lines: list[str]
    ) -> tuple[
        np.ndarray,
        list[str],
        np.ndarray,
        list[list[int]] | None
    ]:
        labels = np.arange(1, len(lines) + 1)
        symbols = [l.split()[0] for l in lines]
        atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
        notes = [list(map(int, l.split()[4:])) for l in lines]
        return labels, symbols, atomcoords, notes
    
    def parse(self, lines: list[str]) -> Molecule:
        indices_blank = [i for i, line in enumerate(lines) if line.strip() == ""]
        
        title = lines[indices_blank[0] + 1:indices_blank[1]]
        
        lines_charge_mult = lines[indices_blank[1] + 1]
        charge, mult = self._parse_charge_mult(lines_charge_mult)
        
        lines_atomcoords = lines[indices_blank[1] + 2:indices_blank[2]]
        labels, symbols, atomcoords, notes = self._parse_atomcoords(lines_atomcoords)
        
        return Molecule(
            title=title,
            charge=charge,
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            notes=notes,
        )


def read_connectable(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_blank = [i for i, line in enumerate(lines) if line == "\n"]
    header = lines[:indices_blank[1] + 2]
    footer = lines[indices_blank[2]:]
    lines_coord = lines[indices_blank[1] + 2:indices_blank[2]]
    labels = np.arange(1, len(lines_coord) + 1)
    symbols = [line.split()[0] for line in lines_coord]
    atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
    notes = [list(map(int, line.split()[4:])) for line in lines_coord]
    
    return Molecule(
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes,
        header=header,
        footer=footer,
    )