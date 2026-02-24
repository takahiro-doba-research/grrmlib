from pathlib import Path

from ..core import Molecule


class XYZWriter:

    def __init__(self, title: str | None = None) -> None:
        self.title = title  # does not contain \n
    
    def _build_n_atoms(self, mol: Molecule) -> list[str]:
        return [f"{len(mol)}\n"]
    
    def _build_title(self) -> list[str]:
        if self.title is not None:
            return [f"{self.title}\n"]
        else:
            return ["\n"]
    
    def _build_atomcoords(self, mol: Molecule) -> list[str]:
        lines = []
        
        for s, (x, y, z) in mol.iter_atoms():
            lines.append(
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}\n"
            )
        
        return lines
    
    def build(self, mol: Molecule) -> str:
        lines = []
        lines += self._build_n_atoms(mol)
        lines += self._build_title()
        lines += self._build_atomcoords(mol)
        return "".join(lines)
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "molecule.xyz",
        *,
        overwrite: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        text = self.build(mol)
        
        with path.open("w" if overwrite else "x") as f:
            f.write(text)
        
        return path