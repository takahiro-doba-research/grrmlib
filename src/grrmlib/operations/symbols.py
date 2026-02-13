import copy


def with_symbols_from(mol0, mol1):
    mol0_new = mol0.copy()
    mol0_new.symbols = copy.copy(mol1.symbols)
    return mol0_new