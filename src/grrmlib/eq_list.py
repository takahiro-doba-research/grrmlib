import re

import numpy as np

from .molecule import EQ
from .molecules import Molecules


def read_eq_list(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_num = [i for i, line in enumerate(lines) if line.startswith("#")]
    indices_num.append(len(lines))
    mols = Molecules()
    
    for i, j in zip(indices_num, indices_num[1:]):
        lines_eq = lines[i:j]
        name = int(re.search(r"EQ (\d+),", lines_eq[0]).group(1))
        index_energy = [i for i, line in enumerate(lines_eq) if line.startswith("Energy")][0]
        lines_coord = lines_eq[1:index_energy]
        labels = np.arange(1, len(lines_coord) + 1)
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy = float(re.search(r"\(\s*(-?\d+\.?\d+)\s*:", lines_eq[index_energy]).group(1))
        afirenergy = float(re.search(r"=\s*(-?\d+\.?\d+)\s*\(", lines_eq[index_energy]).group(1))
        mult = int(float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 1]).group(1)))
        zpve = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 2]).group(1))
        lines_nmeigen = lines_eq[index_energy + 4:]
        nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
        eq = EQ(
            name=name,
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            afirenergy=afirenergy,
            zpve=zpve,
            nmeigen=nmeigen,
        )
        mols[name] = eq
    
    return mols