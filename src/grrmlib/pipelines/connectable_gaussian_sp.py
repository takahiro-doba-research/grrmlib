from pathlib import Path

from ..core import Molecules
from ..io import read_connectable, write_gaussian_inputs
from ..operations import connect_all, step_to_angles
from ._templates import read_gaussian_template


def connection_keys(
    folder_seq: str | Path,
    folders_sub: list[str | Path],
    *,
    with_angles: bool = True
) -> list[str]:
    """
    connectables_to_gaussian_sp が書くフォルダの key の名前。
    例: ["seq", "arene", "arene_angle", "alkene", "alkene_angle"]
    with_angles=False では角度の列を除く（置換基の組み合わせごとに選ぶときの group_by）。
    """
    names = [Path(f).name for f in folders_sub]
    names = [y for x in names for y in ([x, f"{x}_angle"] if with_angles else [x])]
    return [Path(folder_seq).name] + names


def _read_connectables(folder: str | Path) -> Molecules:
    # folder/1.com, folder/2.com, ... を key が 1, 2, ... の Molecules にする
    paths = sorted(Path(folder).glob("*.com"))
    return Molecules({int(p.stem): read_connectable(p) for p in paths})


def connectables_to_gaussian_sp(
    folder_seq: str | Path,
    folders_sub: list[str | Path],
    path_template: str | Path,
    folder: str | Path,
    *,
    basename: str = "gaussian_sp.com",
    overwrite: bool = False
) -> None:
    """
    folder_seq の骨格に folders_sub の置換基を順につなぎ、全ての組み合わせと
    回転角について Gaussian の入力を書く。

    - folder_seq, folders_sub には connectable の入力を 1.com, 2.com, ... の名前で置く
    - 入力は folder/seq=1/arene=0/arene_angle=120/.../basename のように書く
      （名前は connection_keys）
    - link0（%chk を除く）, route, title は path_template の Gaussian の入力から写す
    """
    for f in [folder_seq, *folders_sub]:
        if not Path(f).is_dir():
            raise FileNotFoundError(f"{f} not found.")

    seqs = _read_connectables(folder_seq)
    subs_list = [
        _read_connectables(f).expand(step_to_angles).flatten()
        for f in folders_sub
    ]

    write_gaussian_inputs(
        connect_all(seqs, subs_list),
        folder,
        basename,
        prefixes=connection_keys(folder_seq, folders_sub),
        with_notes=False,
        overwrite=overwrite,
        **read_gaussian_template(path_template)
    )
