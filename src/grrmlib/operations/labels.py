import copy


def with_labels_from(mol0, mol1):
    mol0_new = mol0.copy()
    mol0_new.labels = copy.copy(mol1.labels)
    return mol0_new