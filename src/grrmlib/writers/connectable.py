from pathlib import Path

from ..core import Molecule, GroupedMolecules


class ConnectableWriter:
    
    def __init__(
        self,
        *,
        title: list[str] | None = None,
    ) -> None:
        self.title = title
    
    def _build_route(self) -> list[str]:
        return ["#"]
    
    def _build_title(self) -> list[str]:
        if self.title is not None:
            return self.title
        else:
            return ["None"]
    
    def _build_charge_mult(self, mol: Molecule) -> list[str]:
        charge = mol.charge if mol.charge is not None else 0
        mult = mol.mult if mol.mult is not None else 1
        return [f"{charge} {mult}"]
    
    def _build_atomcoords(self, mol: Molecule) -> list[str]:
        lines = []
        
        for s, (x, y, z), n in mol.iter_atoms(with_notes=True):
            lines.append(
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
                f" {' '.join(map(str, n))}"
            )
        
        return lines
    
    def build(self, mol: Molecule) -> str:
        lines = []
        lines += self._build_route()
        lines.append("")
        lines += self._build_title()
        lines.append("")
        lines += self._build_charge_mult(mol)
        lines += self._build_atomcoords(mol)
        lines.append("")
        return "\n".join(lines)
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "connectable.com",
        *,
        overwrite: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        text = self.build(mol)
        
        with path.open("w" if overwrite else "x") as f:
            f.write(text)
        
        return path
    
    def write_grouped(
        self,
        gmols: GroupedMolecules,
        folder: str | Path,
        prefix_group: str | Path = "group",
        prefix_mol: str | Path = "molecule",
        basename: str | Path = "connectable.com",
        *,
        overwrite: bool = False
    ) -> None:
        folder = Path(folder)
        
        for group, mols in gmols.items():
            for name, mol in mols.items():
                folder_new = (
                    folder
                    / f"{prefix_group}={group}"
                    / f"{prefix_mol}={name}"
                )
                self.write(mol, folder_new / basename, overwrite=overwrite)