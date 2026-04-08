from pathlib import Path
from typing import Iterable

from ..core import Molecule, Molecules


class GRRMInputWriter:
    
    def __init__(
        self,
        *,
        infile: str | None = None,
        route: str | None = None,
        title: str | None = None,
        options: list[str] | None = None,
        with_notes: bool = False,
    ) -> None:
        self.infile = infile
        self.route = route
        self.title = title
        self.options = options
        self.with_notes = with_notes
    
    def _build_link0(self) -> list[str]:
        if self.infile is not None:
            return [f"%infile={self.infile}"]
        else:
            return []
    
    def _build_route(self) -> list[str]:
        if self.route is not None:
            return [self.route]
        else:
            return ["#"]
    
    def _build_title(self) -> list[str]:
        if self.title is not None:
            return [self.title]
        else:
            return [""]
    
    def _build_charge_mult(self, mol: Molecule) -> list[str]:
        charge = mol.charge if mol.charge is not None else 0
        mult = mol.mult if mol.mult is not None else 1
        return [f"{charge} {mult}"]
    
    def _build_atomcoords(self, mol: Molecule) -> list[str]:
        lines = []
        
        if self.with_notes:
            for s, (x, y, z), n in mol.iter_atoms(with_notes=True):
                lines.append(
                    f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
                    f" {' '.join(map(str, n))}"
                )
        else:
            for s, (x, y, z) in mol.iter_atoms():
                lines.append(
                    f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
                )
        
        return lines
    
    def _build_options(self) -> list[str]:
        if self.options is not None:
            return ["Options"] + self.options
        else:
            return []
    
    def build(self, mol: Molecule) -> str:
        lines = []
        lines += self._build_link0()
        lines += self._build_route()
        lines += self._build_title()
        lines += self._build_charge_mult(mol)
        lines += self._build_atomcoords(mol)
        lines += self._build_options()
        lines.append("")
        return "\n".join(lines)
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "grrm.com",
        *,
        overwrite: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        text = self.build(mol)
        
        with path.open("w" if overwrite else "x") as f:
            f.write(text)
        
        return path
    
    def write_mols(
        self,
        mols: Molecules,
        folder: str | Path,
        prefix: str | Iterable[str] = "name",
        basename: str | Path = "grrm.com",
        *,
        overwrite: bool = False
    ) -> None:
        folder = Path(folder)
        basename = Path(basename)
        
        prefix_tuple = (prefix, ) if isinstance(prefix, str) else tuple(prefix)
        
        for name, mol in mols.items():
            name_tuple = name if isinstance(name, tuple) else (name, )
            
            if len(prefix_tuple) != len(name_tuple):
                raise ValueError("prefix and name must have the same length")
            
            folder_new = folder
            for prefix, name in zip(prefix_tuple, name_tuple):
                folder_new = folder_new / f"{prefix}={name}"
            
            self.write(mol, folder_new / basename, overwrite=overwrite)