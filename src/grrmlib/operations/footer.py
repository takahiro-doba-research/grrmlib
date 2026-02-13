import copy


def with_footer_from(mol0, mol1):
    mol0_new = mol0.copy()
    mol0_new.footer = copy.copy(mol1.footer)
    return mol0_new