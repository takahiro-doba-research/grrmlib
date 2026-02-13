import copy


def with_header_from(mol0, mol1):
    mol0_new = mol0.copy()
    mol0_new.header = copy.copy(mol1.header)
    return mol0_new