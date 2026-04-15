from pathlib import Path

import numpy as np
import polars as pl

from ..core import Molecules
from ..operations import step_to_angles, connect_all_iter
from ..readers import (
    GaussianInputReader,
    GaussianOutputReader,
    GRRMInputReader,
    GRRMMINOutputReader,
    ConnectableReader
)
from ..writers import (
    GaussianInputWriter,
    GRRMInputWriter,
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
        
        self.folders = [folder_seq] + [y for x in folders_sub for y in [x, f"{x}_angle"]]
        
        if not Path(folder_seq).is_dir():
            raise FileNotFoundError(f"{folder_seq} not found.")
        
        for folder in folders_sub:
            if not Path(folder).is_dir():
                raise FileNotFoundError(f"{folder} not found.")
    
    def write_gaussian_sp(self, old=False) -> None:
        reader = ConnectableReader()
        paths_seq = sorted(Path(self.folder_seq).glob("*.com"))
        seqs = Molecules({int(p.stem): reader.read(p) for p in paths_seq})
        
        # For the old system in which the notes of the separated_EQs were not relabeled.
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
        writer.write_mols(
            mols,
            "connection",
            self.folders,
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
        df = exporter.export(
            seqs,
            self.folders,
            ["scfenergy", "success"]
        )
        df.write_csv("gaussian_sp.csv")
        
        df_error = df.filter(pl.col("success") == False)
        print(f"{df_error.height} gaussian_sp jobs failed.")
        
        df_selected = (
            pl.read_csv("gaussian_sp.csv")
            .group_by([self.folder_seq] + self.folders_sub, maintain_order=True)
            .agg(pl.all().sort_by("scfenergy", nulls_last=True, maintain_order=True).first())
        )
        df_selected.write_csv("gaussian_sp_selected.csv")
        
        df_selected_error = df_selected.filter(pl.col("success") == False)
        print(f"{df_selected_error.height} gaussian_sp selected jobs failed.")
    
    def write_grrm_min(self, *, overwrite: bool = False) -> None:
        df = pl.read_csv("gaussian_sp_selected.csv")
        
        reader = GaussianInputReader()
        mols = reader.read_mols("connection", "gaussian_sp.com")
        mols_selected = Molecules({
            row: mols[row]
            for row in df.select(self.folders).iter_rows()
        })
        
        reader = GRRMInputReader()
        mol_method = reader.read(self.path_grrm_min)
        writer = GRRMInputWriter(
            route=mol_method.route,
            options=mol_method.options
        )
        writer.write_mols(
            mols_selected,
            "connection",
            self.folders,
            "grrm_min.com",
            overwrite=overwrite
        )
    
    def list_grrm_min(self) -> None:
        paths = sorted(Path("connection").rglob("grrm_min.com"))
        Path("grrm_min.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} grrm_min jobs listed.")
    
    def list_grrm_min_error(self) -> None:
        paths = sorted(Path("connection").rglob("grrm_min_message_ERROR.rrm"))
        paths_new = [str(p.with_name("grrm_min.com")) for p in paths]
        Path("grrm_min_error.txt").write_text("\n".join(paths_new))
        print(f"{len(paths_new)} grrm_min error jobs listed.")
    
    def write_grrm_min_continue(self) -> None:
        paths = sorted(Path("connection").rglob("grrm_min_message_ERROR.rrm"))
        count = 0
        
        for path in paths:
            path_new = path.with_name("grrm_min_message_CONTINUE.rrm")
            path_new.touch()
            count += 1
        
        print(f"{count} grrm_min_message_CONTINUE.rrm files created")
    
    def list_grrm_min_unprocessed(self) -> None:
        paths = sorted(Path("connection").rglob("grrm_min.com"))
        paths_new = []
        
        for path in paths:
            files = [f.name for f in path.parent.iterdir()]
            
            if (
                "grrm_min_message_ERROR.rrm" not in files
                and "grrm_min.log" not in files
            ):
                paths_new.append(str(path))
        
        Path("grrm_min_unprocessed.txt").write_text("\n".join(paths_new))
        print(f"{len(paths_new)} grrm_min unprocessed jobs listed.")
    
    def analyze_grrm_min(self, version: str = "GRRM23") -> None:
        reader = GRRMMINOutputReader(version=version)
        mols = reader.read_mols("connection", "grrm_min.log")
        exporter = PolarsExporter()
        df = exporter.export(
            mols,
            self.folders,
            ["scfenergy", "status"]
        )
        df.write_csv("grrm_min.csv")
        
        df_error = df.filter(pl.col("status") != "Minimum point was found")
        print(f"{df_error.height} grrm_min jobs failed.")