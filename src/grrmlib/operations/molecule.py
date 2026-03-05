from ..core import Molecule, Molecules


def expand_charge_mult(
    mol: Molecule,
    charges: list[int] | int,
    mults: list[int] | int
) -> Molecules:
    if isinstance(charges, int):
        charges = [charges]
    
    if isinstance(mults, int):
        mults = [mults]
    
    mols = Molecules()
    
    for charge in charges:
        for mult in mults:
            mol_new = mol.copy()
            mol_new.charge = charge
            mol_new.mult = mult
            if mol_new.charge_mult_is_valid():
                mols[(charge, mult)] = mol_new
    
    return mols