from pathlib import Path

from ..core import Molecule


class XYZWriter:
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "molecule.xyz",
        *,
        exist_ok: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        lines = [
            f"{len(mol)}\n",
            f"{getattr(mol, 'comments', '') or ''}\n"
        ]
        
        for s, (x, y, z) in mol.iter_atoms():
            lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}\n")
        
        with path.open("w" if exist_ok else "x") as f:
            f.writelines(lines)
        
        return path