from pathlib import Path

import numpy as np
import polars as pl

from ..core import Molecules
from ..operations import step_to_angles, connect_all_iter
from ..readers import (
    GaussianInputReader,
    GaussianOutputReader,
    ConnectableReader
)
from ..writers import (
    GaussianInputWriter,
    PolarsExporter,
)


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
    
    def write_gaussian_sp(self, old=False) -> None:
        reader = ConnectableReader()
        paths_seq = sorted(Path(self.folder_seq).glob("*.com"))
        seqs = Molecules({int(p.stem): reader.read(p) for p in paths_seq})
        
        if old:
            df_labels = pl.read_csv("separated_group_info.csv")
            df_labels = df_labels.with_columns(
                pl.col("labels")
                .str.split("_")
                .list.eval(pl.element().cast(pl.Int64))
            )
            df_sgroup = pl.read_csv("separated_group.csv")
            
            for name, seq in seqs.items():
                sgroup = df_sgroup.filter(pl.col("separated_EQ") == name)["separated_group"][0]
                labels = df_labels.filter(pl.col("separated_group") == sgroup)["labels"].to_list()[0]
                seq.labels = np.array(labels)
        
        paths_sub_list = [sorted(Path(f).glob("*.com")) for f in self.folders_sub]
        subs_list = [
            Molecules({int(p.stem): reader.read(p) for p in paths_sub})
            for paths_sub in paths_sub_list
        ]
        subs_list = [
            subs.expand(step_to_angles).flatten()
            for subs in subs_list
        ]
        
        mols = Molecules({
            tuple(int(n) if n is not None else n for n in names): seq
            for names, seq in connect_all_iter(seqs, subs_list)
        })
        
        reader = GaussianInputReader()
        mol_method = reader.read(self.path_gaussian_sp)
        writer = GaussianInputWriter(
            nprocshared=mol_method.nprocshared,
            mem=mol_method.mem,
            route=mol_method.route,
            title=mol_method.title,
        )
        folder = [self.folder_seq] + [y for x in self.folders_sub for y in [x, "angle"]]
        writer.write_mols(
            mols,
            "connection",
            folder,
            "gaussian_sp.com"
        )
    
    def list_gaussian_sp(self) -> None:
        paths = sorted(Path("connection").rglob("gaussian_sp.com"))
        Path("gaussian_sp.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} gaussian_sp jobs listed.")
    
    def analyze_gaussian_sp(self) -> None:
        reader = GaussianOutputReader()
        seqs = reader.read_mols("connection", "gaussian_sp.log")
        exporter = PolarsExporter()
        folder = [self.folder_seq] + [y for x in self.folders_sub for y in [x, f"{x}_angle"]]
        df = exporter.export(
            seqs,
            folder,
            ["scfenergy", "success"]
        )
        df.write_csv("gaussian_sp.csv")
        
        df_error = df.filter(pl.col("success") == False)
        print(f"{df_error.height} gaussian_sp jobs failed.")
    
    def write_grrm_min(self) -> None:
        pass
    
    def list_grrm_min(self) -> None:
        paths = sorted(Path("connection").rglob("grrm_min.com"))
        Path("grrm_min.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} grrm_min jobs listed.")
    
    def analyze_grrm_min(self) -> None:
        pass