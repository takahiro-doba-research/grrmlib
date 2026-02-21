from pathlib import Path

from ..core import (
    Molecule,
    Molecules,
    GroupedMolecules,
    atomic_number
)


class GaussianInputWriter:
    
    def __init__(
        self,
        *,
        nprocshared: int | None = None,
        mem: str | None = None,
        chk: str | None = None,
        route: list[str] | None = None,
        title: list[str] | None = None,
        with_notes: bool = False,
    ) -> None:
        self.nprocshared = nprocshared
        self.mem = mem
        self.chk = chk
        self.route = route
        self.title = title
        self.with_notes = with_notes
    
    def _build_link0(self) -> list[str]:
        lines = []
        
        if self.nprocshared is not None:
            lines.append(f"%NProcShared={self.nprocshared}\n")
        
        if self.mem is not None:
            lines.append(f"%Mem={self.mem}\n")
        
        if self.chk is not None:
            lines.append(f"%Chk={self.chk}\n")
        
        return lines
    
    def _build_route(self) -> list[str]:
        if self.route is not None:
            return self.route
        else:
            return ["#\n"]
    
    def _build_title(self) -> list[str]:
        if self.title is not None:
            return self.title
        else:
            return ["None\n"]
    
    def _build_charge_mult(self, mol: Molecule) -> list[str]:
        charge = mol.charge if mol.charge is not None else 0
        mult = mol.mult if mol.mult is not None else 1
        return [f"{charge} {mult}\n"]
    
    def _build_atomcoords(self, mol: Molecule) -> list[str]:
        lines = []
        
        if self.with_notes:
            for s, (x, y, z), n in mol.iter_atoms(with_notes=True):
                lines.append(
                    f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}"
                    f" {' '.join(map(str, n))}\n"
                )
        else:
            for s, (x, y, z) in mol.iter_atoms():
                lines.append(
                    f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}\n"
                )
        
        return lines
    
    def build(self, mol: Molecule) -> str:
        lines = []
        lines += self._build_link0()
        lines += self._build_route()
        lines.append("\n")
        lines += self._build_title()
        lines.append("\n")
        lines += self._build_charge_mult(mol)
        lines += self._build_atomcoords(mol)
        lines.append("\n")
        return "".join(lines)
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "gaussian_input.com",
        *,
        exist_ok: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        text = self.build(mol)
        
        with path.open("w" if exist_ok else "x") as f:
            f.write(text)
        
        return path
    
    def write_grouped(
        self,
        gmols: GroupedMolecules,
        folder: str | Path,
        prefix_group: str | Path = "group",
        prefix_mol: str | Path = "molecule",
        basename: str | Path = "gaussian_input.com",
        *,
        exist_ok: bool = False
    ) -> None:
        folder = Path(folder)
        
        for group, mols in gmols.items():
            for name, mol in mols.items():
                folder_new = (
                    folder
                    / f"{prefix_group}{group}"
                    / f"{prefix_mol}{name}"
                )
                self.write(mol, folder_new / basename, exist_ok=exist_ok)


class GaussianOutputWriter:
    
    def write_scan(
        self,
        mols: Molecules,
        path: str | Path = "gaussian_scan.log",
        *,
        exist_ok: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        num = len(mols)
        lines = [" #p\n", " \n"]
        
        for i, mol in enumerate(mols.values(), start=1):
            lines += [
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                " ---------------------------------------------------------------------\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   1 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                *[
                    f"{j:7d} {atomic_number(s):10d}           0     {x:11.6f} {y:11.6f} {z:11.6f}\n"
                    for j, (s, (x, y, z)) in enumerate(mol.iter_atoms(), start=1)
                ],
                " ---------------------------------------------------------------------\n",
                f" SCF Done:  E({mol.functional or 'B3LYP'}) = {mol.scfenergy:15.12f}     A.U.\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   2 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
                " \n",
            ]
        
        lines += [
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            " Normal termination of Gaussian 16\n"
        ]
        
        with path.open("w" if exist_ok else "x") as f:
            f.writelines(lines)
        
        return path