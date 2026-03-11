import re
from pathlib import Path

import numpy as np

from ..core import Molecule


class GRRMInputReader:
    
    def read(self, path: str | Path) -> Molecule:
        path = Path(path)
        text = path.read_text()
        lines = text.splitlines(keepends=False)
        return self.parse(lines)
    
    def _parse_link0(self, lines: list[str]) -> str | None:
        infile = None
        
        for line in lines:
            if line.lower().startswith("%infile"):
                infile = re.search(r"=\s*(.+)", line).group(1)
        
        return infile
    
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
        if all(not note for note in notes):
            notes = None
        
        return labels, symbols, atomcoords, notes
    
    def parse(self, lines: list[str]) -> Molecule:
        lines_link0 = [l for l in lines if l.startswith("%")]
        infile = self._parse_link0(lines_link0)
        
        index_route = [i for i, l in enumerate(lines) if l.startswith("#")][0]
        route = lines[index_route]
        
        line_title = lines[index_route + 1]
        title = line_title if line_title != "" else None
        
        lines_charge_mult = lines[index_route + 2]
        charge, mult = self._parse_charge_mult(lines_charge_mult)
        
        lines_atomcoords_options = lines[index_route + 3:]
        index_options = [
            i for i, l in enumerate(lines_atomcoords_options)
            if l.startswith("Options")
        ]
        if index_options:
            lines_atomcoords = lines_atomcoords_options[:index_options[0]]
            labels, symbols, atomcoords, notes = self._parse_atomcoords(lines_atomcoords)
            lines_options = lines_atomcoords_options[index_options[0] + 1:]
            options = [l for l in lines_options if l.strip()]
        else:
            lines_atomcoords = [l for l in lines_atomcoords_options if l.strip()]
            labels, symbols, atomcoords, notes = self._parse_atomcoords(lines_atomcoords)
            options = None
        
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