from pathlib import Path

from ..operations import FragmentNotFoundError, connect
from ..writers import GaussianInputWriter


def connection_factory(mols_tuple, gmols_tuples):
    """
    Example
    -------
    mols_tuple = ("SEQ", Molecules({0: Molecule(), 1: Molecule(), ...}))
    gmols_tuples = [
        ("arene", GroupedMolecules({0: Molecules(), 1: Molecules(), ...}),
        ("alkene", GroupedMolecules({0: Molecules(), 1: Molecules(), ...})),
        ("pg", GroupedMolecules({0: Molecules(), 1: Molecules(), ...})),
        ("backbone", GroupedMolecules({0: Molecules(), 1: Molecules(), ...})),
        ("pyridone", GroupedMolecules({0: Molecules(), 1: Molecules(), ...})),
    ]
    """
    
    def connect_recursive(mol, gmols_tuples, names):
        if not gmols_tuples:
            yield mol, names
            return
        
        key, gmols = gmols_tuples[0]
        connected = False
        
        for name, mols in gmols.items():
            for angle, mol_next in mols.items():
                try:
                    mol_new = connect(mol, mol_next)
                    connected = True
                    yield from connect_recursive(
                        mol_new,
                        gmols_tuples[1:],
                        names + [(key, name), ("angle", angle)]
                    )
                except FragmentNotFoundError:
                    break
            if not connected:
                break
        
        if not connected:
            yield from connect_recursive(
                mol,
                gmols_tuples[1:],
                names
            )
    
    names = []
    key, mols = mols_tuple
    
    for name, mol in mols.items():
        yield from connect_recursive(
            mol,
            gmols_tuples,
            names + [(key, name)]
        )


def build_connection_folder(folder, names):
    """
    Example
    -------
    names = [
        ('SEQ', 7286),
        ('arene', 0), ('angle', 0.0),
        ('alkene', 0), ('angle', 240.0),
        ('pg', 0), ('angle', 0.0),
        ('backbone', 5), ('angle', 240.0),
        ('pyridone', 9), ('angle', 0.0)
    ]
    """
    folder_new = Path(folder)
    
    for key, value in names:
        value = round(value) if key == "angle" else value
        folder_new /= f"{key}{value}"
    
    return folder_new


def connect_all_and_write_as_gaussian_input(
    mols_tuple,
    gmols_tuples,
    mol_method,
    folder="connection",
    basename="gaussian_input.com"
):
    writer = GaussianInputWriter(
        header=mol_method.header,
        footer=mol_method.footer,
    )
    
    for mol, names in connection_factory(mols_tuple, gmols_tuples):
        folder_new = build_connection_folder(folder, names)
        folder_new.mkdir(parents=True, exist_ok=True)
        
        try:
            writer.write(mol, folder_new / basename)
        except FileExistsError as e:
            print(e)