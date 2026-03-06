import networkx as nx

from ..core import Molecule, Molecules


def separate(mol: Molecule) -> Molecules:
    mol.validate()
    A = mol.get_adj_matrix()
    G = nx.from_numpy_array(A)
    #components = sorted(nx.connected_components(G), key=len, reverse=True)
    components = sorted(nx.connected_components(G), key=min)
    
    mols = Molecules()
    for i, component in enumerate(components):
        indices = sorted(component)
        fragment = mol._select_by_indices(indices)
        mols[i] = fragment
    
    return mols


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