from __future__ import annotations

import time
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator

    from ..core import Molecule


# ASE は eV、grrmlib は Hartree を使う
HARTREE_IN_EV = 27.211386245988

OPTIMIZERS = {
    "LBFGS": LBFGS,
    "BFGS": BFGS,
    "FIRE": FIRE,
}


def load_uma(
    model_name: str = "uma-s-1p1",
    task_name: str = "omol",
    device: str = "cpu"
) -> Calculator:
    # fairchem は重いので必要になったときだけ import する
    from fairchem.core import FAIRChemCalculator, pretrained_mlip

    predictor = pretrained_mlip.get_predict_unit(model_name, device=device)
    return FAIRChemCalculator(predictor, task_name=task_name)


def to_atoms(mol: Molecule) -> Atoms:
    mol.validate()

    if mol.charge is None or mol.mult is None:
        raise ValueError("charge and mult must be set")

    atoms = Atoms(symbols=list(mol.symbols), positions=mol.atomcoords)
    # OMol25 系のモデルは info の charge と spin (= 多重度) を参照する
    atoms.info["charge"] = mol.charge
    atoms.info["spin"] = mol.mult
    return atoms


def single_point(mol: Molecule, calc: Calculator) -> Molecule:
    t = time.perf_counter()
    atoms = to_atoms(mol)
    atoms.calc = calc

    mol_new = mol.copy()
    mol_new.scfenergy = atoms.get_potential_energy() / HARTREE_IN_EV
    mol_new.success = True
    mol_new.status = "SP"
    mol_new.elapsed = time.perf_counter() - t
    return mol_new


def optimize(
    mol: Molecule,
    calc: Calculator,
    *,
    fmax: float = 0.02,
    steps: int = 2000,
    optimizer: str = "LBFGS"
) -> Molecule:
    """
    fmax は eV/Å 単位（0.02 eV/Å ≈ 3.9e-4 Hartree/Bohr。Gaussian の既定の
    最大力の閾値 4.5e-4 Hartree/Bohr と同程度）。
    """
    t = time.perf_counter()
    atoms = to_atoms(mol)
    atoms.calc = calc

    opt = OPTIMIZERS[optimizer](atoms, logfile=None)
    converged = bool(opt.run(fmax=fmax, steps=steps))

    mol_new = mol.copy()
    mol_new.atomcoords = np.array(atoms.get_positions())
    mol_new.scfenergy = atoms.get_potential_energy() / HARTREE_IN_EV
    mol_new.success = converged
    mol_new.status = "Converged" if converged else "Not converged"
    mol_new.optsteps = opt.nsteps
    mol_new.elapsed = time.perf_counter() - t
    return mol_new
