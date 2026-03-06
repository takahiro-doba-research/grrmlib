from pathlib import Path

from ..readers import (
    GaussianInputReader,
    ConnectableReader,
    read_eq_list,
)
from ..core import is_identical
from ..operations import (
    separate,
    product_charge_mult,
)
from ..writers import (
    GaussianInputWriter,
    ConnectableWriter,
)


class SeparationWorkflow:
    
    def __init__(
        self,
        *,
        folder: str | Path,
        path_eq_list: str | Path,
        path_gaussian_sp: str | Path,
        path_grrm_min: str | Path,
    ) -> None:
        self.folder = Path(folder)
        self.path_eq_list = Path(path_eq_list)
        self.path_gaussian_sp = Path(path_gaussian_sp)
        self.path_grrm_min = Path(path_grrm_min)
        
    def write_gaussian_sp(
        self,
        nprocshared: int,
        mem: str,
        charges: int | list[int],
        mults: int | list[int],
    ) -> None:
        eqs = read_eq_list(self.path_eq_list)
        seqs = eqs.expand(lambda eq: separate(eq)).flatten().reset_keys()
        seqs = seqs.cluster(is_identical).flatten()
        seqs = seqs.expand(lambda seq: product_charge_mult(seq, charges, mults)).flatten()
        
        mol_method = GaussianInputReader().read(self.path_gaussian_sp)
        writer = GaussianInputWriter(
            nprocshared=nprocshared,
            mem=mem,
            route=mol_method.route,
            title=mol_method.title,
        )
        writer.write_mols(
            seqs,
            self.folder,
            ["SG", "SEQ", "charge", "mult"],
            "gaussian_sp.com"
        )
    
    def list_gaussian_sp(self) -> None:
        paths = sorted(self.folder.rglob("gaussian_sp.com"))
        Path("gaussian_sp.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} jobs listed.")
    
    def run_gaussian_sp():
        pass

    def write_grrm_min():
        pass

    def run_grrm_min():
        pass