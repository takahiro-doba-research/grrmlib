from ..core import Molecule, Molecules


def product_charge_mult(
    mol: Molecule,
    charges: list[int] | int,
    mults: list[int] | int
) -> Molecules:
    
    charges = [charges] if isinstance(charges, int) else charges
    mults = [mults] if isinstance(mults, int) else mults
    
    mols = Molecules()
    
    for charge in charges:
        for mult in mults:
            mol_new = mol.copy()
            mol_new.charge = charge
            mol_new.mult = mult
            if mol_new.charge_mult_is_valid():
                mols[(charge, mult)] = mol_new
    
    return mols