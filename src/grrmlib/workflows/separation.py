from pathlib import Path

import polars as pl

from ..core import GroupedMolecules, Molecules, is_identical, separate
from ..io import (
    read_eq_list,
    read_gaussian_input,
    read_gaussian_inputs,
    read_gaussian_outputs,
    read_grrm_input,
    read_grrm_min_outputs,
    to_dataframe,
    write_gaussian_inputs,
    write_grrm_inputs,
)
from ..operations import product_charge_mult


# separation/SG=1/SEQ=2/charge=0/mult=1/ のようにフォルダを分ける
KEYS = ["SG", "SEQ", "charge", "mult"]


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
        # キーは (EQ, fragment)。その順に SEQ の番号を振る
        frags = eqs.expand(lambda eq: separate(eq)).flatten()
        sgs = frags.reset_keys().cluster(is_identical)
        df = self._build_table(eqs, frags, sgs)
        
        seqs = sgs.flatten()
        seqs = seqs.expand(lambda seq: product_charge_mult(seq, charges, mults)).flatten()
        
        mol_method = read_gaussian_input(self.path_gaussian_sp)
        write_gaussian_inputs(
            seqs,
            "separation",
            "gaussian_sp.com",
            prefixes=KEYS,
            # chk は各ジョブで要らないので引き継がない
            link0=[l for l in mol_method.link0 if not l.lower().startswith("%chk")],
            route=mol_method.route,
            title=mol_method.title,
            with_notes=False,
            overwrite=overwrite
        )
        df.write_csv("separation.csv")
    
    @staticmethod
    def _build_table(
        eqs: Molecules,
        frags: Molecules,
        sgs: GroupedMolecules
    ) -> pl.DataFrame:
        """
        G, EQ, SEQ, SG の対応表を作る。1 行が 1 つの SEQ。
        G は EQ を、SG は SEQ を is_identical でまとめたグループ。
        labels は SEQ の原子の元の EQ でのラベルを "_" でつないだもの。
        """
        eq_to_g = {
            eq: g
            for g, mols in eqs.cluster(is_identical).items()
            for eq in mols
        }
        seq_to_sg = {
            seq: sg
            for sg, mols in sgs.items()
            for seq in mols
        }
        rows = [
            {
                "G": eq_to_g[eq],
                "EQ": eq,
                "SEQ": seq,
                "SG": seq_to_sg[seq],
                "fragment": fragment,
                "labels": "_".join(str(l) for l in frag.labels),
            }
            for seq, ((eq, fragment), frag) in enumerate(frags.items())
        ]
        return pl.DataFrame(rows)
    
    def list_gaussian_sp(self) -> None:
        paths = sorted(Path("separation").rglob("gaussian_sp.com"))
        Path("gaussian_sp.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} gaussian_sp jobs listed.")
    
    def analyze_gaussian_sp(self) -> None:
        seqs = read_gaussian_outputs("separation", "gaussian_sp.log", ignore_errors=True)
        df = to_dataframe(seqs, KEYS, ["scfenergy", "success"])
        df.write_csv("gaussian_sp.csv")
        
        # log を読めなかったジョブは null なので、それも失敗に数える
        df_error = df.filter(pl.col("success").ne_missing(True))
        print(f"{df_error.height} gaussian_sp jobs failed.")
    
    def write_grrm_min(self, *, overwrite: bool = False) -> None:
        df = (
            pl.read_csv("gaussian_sp.csv")
            .group_by(["SG", "charge", "mult"], maintain_order=True)
            .agg(pl.all().sort_by("scfenergy").first())
            .select(["SG", "SEQ", "charge", "mult", "scfenergy", "success"])
        )
        df.write_csv("gaussian_sp_selected.csv")
        
        seqs = read_gaussian_inputs("separation", "gaussian_sp.com", ignore_errors=True)
        seqs_selected = Molecules({
            row: seqs[row]
            for row in df.select(KEYS).iter_rows()
        })
        
        mol_method = read_grrm_input(self.path_grrm_min)
        write_grrm_inputs(
            seqs_selected,
            "separation",
            "grrm_min.com",
            prefixes=KEYS,
            route=mol_method.route,
            options=mol_method.options,
            with_notes=False,
            overwrite=overwrite
        )
    
    def list_grrm_min(self) -> None:
        paths = sorted(Path("separation").rglob("grrm_min.com"))
        Path("grrm_min.txt").write_text("\n".join(str(p) for p in paths))
        print(f"{len(paths)} grrm_min jobs listed.")
    
    def analyze_grrm_min(self) -> None:
        seqs = read_grrm_min_outputs("separation", "grrm_min.log", ignore_errors=True)
        df = to_dataframe(seqs, KEYS, ["scfenergy", "status"])
        df.write_csv("grrm_min.csv")
        
        # log を読めなかったジョブは null なので、それも失敗に数える
        df_error = df.filter(pl.col("status").ne_missing("Minimum point was found"))
        print(f"{df_error.height} grrm_min jobs failed.")