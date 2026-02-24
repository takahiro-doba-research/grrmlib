from pathlib import Path

import numpy as np

from ..core import Molecule


class XYZReader:
    
    def read(self, path: str | Path) -> Molecule:
        path = Path(path)
        text = path.read_text()
        lines = text.splitlines(keepends=True)
        return self.parse(lines)
    
    def _parse_n_atoms(self, line: str) -> int:
        return int(line.strip())
    
    def _parse_title(self, line: str) -> str:
        return line.rstrip("\r\n")
    
    def _parse_atomcoords(
        self,
        lines: list[str]
    ) -> tuple[np.ndarray, list[str], np.ndarray]:
        labels = np.arange(1, len(lines) + 1)
        symbols = [l.split()[0] for l in lines]
        atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
        return labels, symbols, atomcoords
    
    def parse(self, lines: list[str]) -> Molecule:
        n_atoms = self._parse_n_atoms(lines[0])
        title = self._parse_title(lines[1])
        labels, symbols, atomcoords = self._parse_atomcoords(lines[2:2 + n_atoms])
        
        return Molecule(
            title=title,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
        )