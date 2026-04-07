import re
from pathlib import Path

import cclib
import numpy as np
from periodictable import elements

from ..core import Molecule, Molecules


class GaussianInputReader:
    
    def read_mols(
        self,
        folder: str | Path,
        basename: str = "gaussian.com"
    ) -> Molecules:
        paths = sorted(Path(folder).rglob(basename))
        mols = Molecules()
        
        for path in paths:
            key = tuple(p.split("=")[1] for p in path.parts[1:-1])
            key = tuple(None if k == "None" else int(k) for k in key)
            try:
                mols[key] = self.read(path)
            except Exception as e:
                mols[key] = Molecule()
                print(key, e)
        
        return mols
    
    def read(self, path: str | Path) -> Molecule:
        path = Path(path)
        text = path.read_text()
        lines = text.splitlines()
        return self.parse(lines)
    
    def _parse_link0(
        self,
        lines: list[str]
    ) -> tuple[int | None, str | None, str | None]:
        nprocshared = None
        mem = None
        chk = None
        
        for line in lines:
            if line.lower().startswith("%nprocshared"):
                nprocshared = int(re.search(r"=\s*(\d+)", line).group(1))
            elif line.lower().startswith("%mem"):
                mem = re.search(r"=\s*(.+)", line).group(1)
            elif line.lower().startswith("%chk"):
                chk = re.search(r"=\s*(.+)", line).group(1)
        
        return nprocshared, mem, chk
    
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
        indices_blank = [i for i, line in enumerate(lines) if line.strip() == ""]
        
        lines_link0_route = lines[:indices_blank[0]]
        lines_link0 = [l for l in lines_link0_route if l.startswith("%")]
        nprocshared, mem, chk = self._parse_link0(lines_link0)
        
        index_route = [i for i, l in enumerate(lines_link0_route) if l.startswith("#")][0]
        route = lines_link0_route[index_route:]
        
        title = lines[indices_blank[0] + 1:indices_blank[1]]
        
        lines_charge_mult = lines[indices_blank[1] + 1]
        charge, mult = self._parse_charge_mult(lines_charge_mult)
        
        lines_atomcoords = lines[indices_blank[1] + 2:indices_blank[2]]
        labels, symbols, atomcoords, notes = self._parse_atomcoords(lines_atomcoords)
        
        return Molecule(
            nprocshared=nprocshared,
            mem=mem,
            chk=chk,
            route=route,
            title=title,
            charge=charge,
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            notes=notes,
        )


class GaussianOutputReader:
    
    def read(self, path: str | Path) -> Molecule:
        path = Path(path)
        parser = cclib.io.ccopen(path)
        data = parser.parse()
        
        symbols = [str(elements[n]) for n in data.atomnos]
        labels = np.arange(1, len(symbols) + 1)
        
        return Molecule(
            charge=data.charge,
            mult=data.mult,
            labels=labels,
            symbols=symbols,
            atomcoords=data.atomcoords[-1],
            scfenergy=data.scfenergies[-1],
            success=data.metadata["success"],
        )
    
    def read_mols(
        self,
        folder: str | Path,
        basename: str = "gaussian.log"
    ) -> Molecules:
        paths = sorted(Path(folder).rglob(basename))
        mols = Molecules()
        
        for path in paths:
            key = tuple(p.split("=")[1] for p in path.parts[1:-1])
            key = tuple(None if k == "None" else int(k) for k in key)
            try:
                mol = self.read(path)
                mols[key] = mol
            except Exception as e:
                mols[key] = Molecule()
                print(key, e)
        
        return mols


def read_gaussian_input(path):
    
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
    if all(not note for note in notes):
        notes = None
    
    return Molecule(
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes,
        header=header,
        footer=footer,
    )