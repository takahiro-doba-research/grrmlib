import numpy as np


def is_identical(mol0, mol1):
    """
    True if molecules are strictly identical in:
    - atom labels
    - atom order and symbols
    - bonding topology (adjacency matrix)
    """
    return (
        np.array_equal(mol0.labels, mol1.labels)
        and (mol0.symbols == mol1.symbols)
        and np.array_equal(mol0.get_adj_matrix(), mol1.get_adj_matrix())
    )