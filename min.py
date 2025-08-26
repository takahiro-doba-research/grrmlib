import re

import numpy as np

from .data import atomic_number
from .molecule import Molecule


class MIN:

    def __init__(self, name=None, itrs=None, optimized=None):
        self.name = name
        self.itrs = itrs
        self.optimized = optimized

    def to_gv(self, path):
        num = len(self.itrs)
        lines = [" #p\n", " \n"]
        
        for i, molecule in enumerate(self.itrs):
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


def read_min(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_optopt = [i for i, line in enumerate(lines) if line.startswith("OPTOPT")]
    lines_optopt = lines[indices_optopt[0]:indices_optopt[1] + 1]
    itrs, optimized = _read_optopt(lines_optopt)
    min_ = MIN(
        itrs=itrs,
        optimized=optimized
    )
    return min_


def _read_optopt(lines):
    """
    Read iterations.
    """
    indices_optopt = [i for i, line in enumerate(lines) if line.startswith("OPTOPT")]
    index_lineline = [i for i, line in enumerate(lines) if line.startswith("======")][0]
    lines_opt = lines[indices_optopt[0]:index_lineline]
    indices_itr = [i for i, line in enumerate(lines_opt) if line.startswith("#")]
    indices_itr.append(len(lines_opt))
    molecules_itr = []
    
    for i, j in zip(indices_itr, indices_itr[1:]):
        lines_itr = lines_opt[i:j]
        name = "ITR" + re.search(r"ITR. (\d+)", lines_itr[0]).group(1)
        index_item = [i for i, line in enumerate(lines_itr) if "Item" in line][0]
        lines_coord = lines_itr[1:index_item]
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy = float(re.search(r"\(\s*(-?\d+\.?\d+)\s*:", lines_itr[index_item + 1]).group(1))
        afirenergy = float(re.search(r"ENERGY\s*(-?\d+\.?\d+)\s*\(", lines_itr[index_item + 1]).group(1))
        mult = float(re.search(r"Spin\(\*\*2\)\s*(-?\d+\.?\d+)", lines_itr[index_item + 2]).group(1))
        index_nmeigen = [i for i, line in enumerate(lines_itr) if line.startswith("NORMAL MODE EIGENVALUE")][0]
        lines_nmeigen = lines_itr[index_nmeigen + 1:]
        nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
        molecule = Molecule(
            name=name,
            mult=mult,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            afirenergy=afirenergy,
            nmeigen=nmeigen,
        )
        molecules_itr.append(molecule)
    
    """
    Read the optimized structure.
    """
    lines_optimized = lines[index_lineline:indices_optopt[1]]
    index_energy = [i for i, line in enumerate(lines_optimized) if line.startswith("ENERGY")][0]
    lines_coord = lines_optimized[2:index_energy]
    symbols = [line.split()[0] for line in lines_coord]
    atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
    afirenergy = float(re.search(r"=\s*(-?\d+\.?\d+)\s*", lines_optimized[index_energy]).group(1))
    mult = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_optimized[index_energy + 1]).group(1))
    zpve = float(re.search(r"=\s*(-?\d+\.?\d+)", lines_optimized[index_energy + 2]).group(1))
    index_grad = [i for i, line in enumerate(lines_optimized) if line.startswith("GRADIENT VECTOR")][0]
    index_hess = [i for i, line in enumerate(lines_optimized) if line.startswith("HESSIAN MATRIX")][0]
    lines_grad = lines_optimized[index_grad + 1:index_hess]
    grads = np.array([float(line.split()[0]) for line in lines_grad])
    index_nmeigen = [i for i, line in enumerate(lines_optimized) if line.startswith("NORMAL MODE EIGENVALUE")][0]
    
    # read hessian
    lines_hess = lines_optimized[index_hess + 1: index_nmeigen]
    indices_hess_chunk = [i for i, line in enumerate(lines_hess) if len(line) <= 15]
    indices_hess_chunk.append(len(lines_hess))
    hesses = []
    
    for i, j in zip(indices_hess_chunk, indices_hess_chunk[1:]):
        lines_hess_chunk = lines_hess[i:j]
        lines_hess_chunk = [list(map(float, line.split())) for line in lines_hess_chunk]
        
        for k in range(len(lines_hess_chunk[-1])):
            hess = [line[k] for line in lines_hess_chunk[k:]]
            hesses.append(hess)
    
    hess_1d = np.concatenate(hesses)
    dim = len(hesses[0])
    hess_2d = np.full((dim, dim), np.nan)
    hess_2d[np.tril_indices(dim)] = hess_1d[::-1]
    hess_2d = hess_2d[::-1, ::-1].T
    
    lines_nmeigen = lines_optimized[index_nmeigen + 1:-2]
    nmeigen = np.concatenate([list(map(float, line.split())) for line in lines_nmeigen])
    status = lines_optimized[-1].strip()
    molecule_optimized = Molecule(
        name="Optimized structure",
        mult=mult,
        symbols=symbols,
        atomcoords=atomcoords,
        afirenergy=afirenergy,
        zpve=zpve,
        grads=grads,
        hessian=hess_2d,
        nmeigen=nmeigen,
        status=status,
    )
    return molecules_itr, molecule_optimized