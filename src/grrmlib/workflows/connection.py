from pathlib import Path

from ..core import Molecules
from ..operations import FragmentNotFoundError, connect
from ..writers import GaussianInputWriter


class ConnectionWorkflow:
    
    def __init__(
        self,
        *,
        folder_seq: str,
        folders_sub: list[str],
        path_gaussian_sp: str | Path,
        path_grrm_min: str | Path,
    ) -> None:
        self.folder_seq = folder_seq
        self.folders_sub = folders_sub
        self.path_gaussian_sp = Path(path_gaussian_sp)
        self.path_grrm_min = Path(path_grrm_min)
    
    def write_gaussian_sp(
        self,
        charges: int | list[int],
        mults: int | list[int],
    ) -> None:
        paths_seq = sorted(Path(self.folder_seq).rglob("*.com"))
        paths_sub_list = [sorted(Path(f).rglob("*.com")) for f in self.folders_sub]
        
        reader = ConnectableReader()
        seqs = reader.read_mols(paths_seq)
        subs_list = [reader.read_mols(p) for p in paths_sub_list]
        
        mol_method = GaussianInputReader().read(self.path_gaussian_sp)
        writer = GaussianInputWriter(
            nprocshared=mol_method.nprocshared,
            mem=mol_method.mem,
            route=mol_method.route,
            title=mol_method.title,
        )
        writer.write_mols(
            seqs,
            "separation",
            ["SG", "SEQ", "charge", "mult"],
            "gaussian_sp.com"
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