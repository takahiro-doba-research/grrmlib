import re

import numpy as np

from .data import atomic_number
from .molecule import Molecule


class EQList:

    def __init__(self, name=None, molecules=None):
        self.name = name
        self.molecules = molecules
    
    def __len__(self):
        return len(self.molecules)
    
    def __getitem__(self, item):
        return self.molecules[item]
    
    def to_gv(self, path):
        num = len(self.molecules)
        lines = [" #p\n", " \n"]
        
        for i, molecule in enumerate(self.molecules):
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
                f" Step number   1 out of a maximum of   2 on scan point {i+1:5d} out of {num:5d}\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                *[
                    f"{i+1:7d} {atomic_number(sym):10d}           0     {coord[0]:11.6f} {coord[1]:11.6f} {coord[2]:11.6f}\n"
                    for i, (sym, coord) in enumerate(zip(molecule.symbols, molecule.atomcoords))
                ],
                " ---------------------------------------------------------------------\n",
                f" SCF Done:  E({molecule.functional}) = {molecule.scfenergy:15.12f}     A.U.\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   2 out of a maximum of   2 on scan point {i+1:5d} out of {num:5d}\n",
                " \n",
            ]
        
        lines += [
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            " Normal termination of Gaussian 16\n"
        ]
        
        with open(path, "w") as f:
            f.writelines(lines)


def read_eq_list(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_num = [i for i, line in enumerate(lines) if line.startswith("#")]
    indices_num.append(len(lines))
    molecules = []
    
    for i, j in zip(indices_num, indices_num[1:]):
        lines_eq = lines[i:j]
        name = "EQ" + re.search(r"EQ (\d+),", lines_eq[0]).group(1)
        index_energy = [i for i, line in enumerate(lines_eq) if line.startswith("Energy")][0]
        lines_coord = lines_eq[1:index_energy]
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy = float(re.search(r"\(\s*(-?\d+\.?\d+)\s*:", lines_eq[index_energy]).group(1))
        afirenergy = float(re.search(r"=\s*(-?\d+\.?\d+)\s*\(", lines_eq[index_energy]).group(1))
        mult = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 1]).group(1))
        zpve = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 2]).group(1))
        lines_nmeigen = lines_eq[index_energy + 4:]
        nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
        molecule = Molecule(
            name=name,
            mult=mult,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            afirenergy=afirenergy,
            zpve=zpve,
            nmeigen=nmeigen,
        )
        molecules.append(molecule)
    
    return EQList(name=path, molecules=molecules)