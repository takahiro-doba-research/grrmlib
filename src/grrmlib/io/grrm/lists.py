import re

import numpy as np

from ...core import Molecule, Molecules
from ._spin import spin2_to_mult


def _parse_energy(line: str) -> tuple[float, float | None]:
    """
    (scfenergy, afirenergy) を返す。Energy の行は 2 通り:
    - AFIR を使った計算: "Energy = <AFIR エネルギー> (<SCF エネルギー> : <AFIR の力の項>)"
    - 使わない計算:     "Energy = <SCF エネルギー>"  （afirenergy は None）
    """
    m = re.search(r"=\s*(-?\d+\.?\d+)\s*\(\s*(-?\d+\.?\d+)\s*:", line)
    
    if m:
        afirenergy, scfenergy = map(float, m.groups())
        return scfenergy, afirenergy
    
    scfenergy = float(re.search(r"=\s*(-?\d+\.?\d+)", line).group(1))
    return scfenergy, None


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
        scfenergy, afirenergy = _parse_energy(lines_eq[index_energy])
        mult = spin2_to_mult(float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 1]).group(1)))
        zpve = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_eq[index_energy + 2]).group(1))
        lines_nmeigen = lines_eq[index_energy + 4:]
        nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
        mol = Molecule(
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
        mols[name] = mol
    
    return mols

def read_pt_list(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_num = [i for i, line in enumerate(lines) if line.startswith("#")]
    indices_num.append(len(lines))
    mols = Molecules()
    
    for i, j in zip(indices_num, indices_num[1:]):
        lines_pt = lines[i:j]
        name = int(re.search(r"TS (\d+),", lines_pt[0]).group(1))
        index_energy = [i for i, line in enumerate(lines_pt) if line.startswith("Energy")][0]
        lines_coord = lines_pt[1:index_energy]
        labels = np.arange(1, len(lines_coord) + 1)
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy, afirenergy = _parse_energy(lines_pt[index_energy])
        mult = spin2_to_mult(float(re.search(r"=\s*(-?\d+\.?\d+)", lines_pt[index_energy + 1]).group(1)))
        zpve = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_pt[index_energy + 2]).group(1))
        index_connection = [i for i, line in enumerate(lines_pt) if line.startswith("CONNECTION")][0]
        lines_nmeigen = lines_pt[index_energy + 4:index_connection]
        nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
        connection = re.search(r":\s*(\d+|\?\?)\s*-\s*(\d+|\?\?)", lines_pt[index_connection]).groups()
        connection = tuple(n if n == "??" else int(n) for n in connection)
        mol = Molecule(
            name=name,
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            afirenergy=afirenergy,
            zpve=zpve,
            nmeigen=nmeigen,
            connection=connection,
        )
        mols[name] = mol
    
    return mols