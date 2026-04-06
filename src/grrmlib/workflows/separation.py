from pathlib import Path

import polars as pl

from ..readers import (
    GaussianInputReader,
    GaussianOutputReader,
    GRRMInputReader,
    GRRMMINOutputReader,
    ConnectableReader,
    read_eq_list,
)
from ..core import Molecules, is_identical
from ..operations import (
    separate,
    product_charge_mult,
)
from ..writers import (
    GaussianInputWriter,
    GRRMInputWriter,
    ConnectableWriter,
    PolarsExporter,
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
        *,
        overwrite: bool = False
    ) -> None:
        eqs = read_eq_list(self.path_eq_list)
        seqs = eqs.expand(lambda eq: separate(eq)).flatten().reset_keys()
        seqs = seqs.cluster(is_identical).flatten()
        seqs = seqs.expand(lambda seq: product_charge_mult(seq, charges, mults)).flatten()
        
        reader = GaussianInputReader()
        mol_method = reader.read(self.path_gaussian_sp)
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
            "gaussian_sp.com",
            overwrite=overwrite
        )
    
    def list_gaussian_sp(self) -> None:
        paths = sorted(Path("separation").rglob("gaussian_sp.com"))
        Path("gaussian_sp.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} gaussian_sp jobs listed.")
    
    def analyze_gaussian_sp(self) -> None:
        reader = GaussianOutputReader()
        seqs = reader.read_mols("separation", "gaussian_sp.log")
        exporter = PolarsExporter()
        df = exporter.export(
            seqs,
            ["SG", "SEQ", "charge", "mult"],
            ["scfenergy", "success"]
        )
        df.write_csv("gaussian_sp.csv")
        
        df_error = df.filter(pl.col("success") == False)
        print(f"{df_error.height} gaussian_sp jobs failed.")
    
    def write_grrm_min(self, *, overwrite: bool = False) -> None:
        df = (
            pl.read_csv("gaussian_sp.csv")
            .group_by(["SG", "charge", "mult"], maintain_order=True)
            .agg(pl.all().sort_by("scfenergy").first())
            .select(["SG", "SEQ", "charge", "mult", "scfenergy", "success"])
        )
        df.write_csv("gaussian_sp_selected.csv")
        
        reader = GaussianInputReader()
        seqs = reader.read_mols("separation", "gaussian_sp.com")
        seqs_selected = Molecules({
            row: seqs[row]
            for row in df.select(["SG", "SEQ", "charge", "mult"]).iter_rows()
        })
        
        reader = GRRMInputReader()
        mol_method = reader.read(self.path_grrm_min)
        writer = GRRMInputWriter(
            route=mol_method.route,
            options=mol_method.options
        )
        writer.write_mols(
            seqs_selected,
            "separation",
            ["SG", "SEQ", "charge", "mult"],
            "grrm_min.com",
            overwrite=overwrite
        )
    
    def list_grrm_min(self) -> None:
        paths = sorted(Path("separation").rglob("grrm_min.com"))
        Path("grrm_min.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} grrm_min jobs listed.")
    
    def analyze_grrm_min(self) -> None:
        reader = GRRMMINOutputReader()
        seqs = reader.read_mols("separation", "grrm_min.log")
        exporter = PolarsExporter()
        df = exporter.export(
            seqs,
            ["SG", "SEQ", "charge", "mult"],
            ["scfenergy", "status"]
        )
        df.write_csv("grrm_min.csv")
        
        df_error = df.filter(pl.col("status") != "Minimum point was found")
        print(f"{df_error.height} grrm_min jobs failed.")