from pathlib import Path

import polars as pl

from ..readers import (
    GaussianInputReader,
    GaussianOutputReader,
    GRRMInputReader,
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
    GRRMInputWriter,
    ConnectableWriter,
)


class SeparationWorkflow:
    
    def __init__(
        self,
        *,
        path_eq_list: str | Path,
        path_gaussian_sp: str | Path,
        path_grrm_min: str | Path,
    ) -> None:
        self.path_eq_list = Path(path_eq_list)
        self.path_gaussian_sp = Path(path_gaussian_sp)
        self.path_grrm_min = Path(path_grrm_min)
    
    def write_gaussian_sp(
        self,
        charges: int | list[int],
        mults: int | list[int],
    ) -> None:
        eqs = read_eq_list(self.path_eq_list)
        seqs = eqs.expand(lambda eq: separate(eq)).flatten().reset_keys()
        seqs = seqs.cluster(is_identical).flatten()
        seqs = seqs.expand(lambda seq: product_charge_mult(seq, charges, mults)).flatten()
        
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
    
    def list_gaussian_sp(self) -> None:
        paths = sorted(Path("separation").rglob("gaussian_sp.com"))
        Path("gaussian_sp.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} gaussian_sp jobs listed.")
    
    def analyze_gaussian_sp(self) -> None:
        reader = GaussianOutputReader()
        seqs = reader.read_mols("separation", "gaussian_sp.log")
        df = pl.DataFrame(
            [(*name, seq.scfenergy, seq.success) for name, seq in seqs.items()],
            schema=["SG", "SEQ", "charge", "mult", "scfenergy", "success"],
            orient="row"
        )
        df.write_csv("gaussian_sp.csv")
        
        df_error = df.filter(pl.col("success") == False)
        print(f"{df_error.height} gaussian_sp jobs failed.")
    
    def write_grrm_min(self) -> None:
        df = (
            pl.read_csv("gaussian_sp.csv")
            .group_by(["SG", "charge", "mult"], maintain_order=True)
            .agg(pl.all().sort_by("scfenergy").first())
            .select(["SG", "SEQ", "charge", "mult", "scfenergy", "success"])
        )
        df.write_csv("gaussian_sp_summary.csv")
        
        paths = sorted(Path("separation").rglob("gaussian_sp.com"))
        folders = [
            path.parent
            for path in paths
            for parts in path.parts
            if parts.split("=")[0] == "SEQ"
            and int(parts.split("=")[1]) in df["SEQ"]
        ]
        reader = GaussianInputReader()
        mol_method = GRRMInputReader().read(self.path_grrm_min)
        writer = GRRMInputWriter(
            route=mol_method.route,
            options=mol_method.options
        )
        
        for folder in folders:
            mol = reader.read(folder / "gaussian_sp.com")
            writer.write(mol, folder / "grrm_min.com", overwrite=True)
    
    def list_grrm_min(self) -> None:
        paths = sorted(Path("separation").rglob("grrm_min.com"))
        Path("grrm_min.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} grrm_min jobs listed.")