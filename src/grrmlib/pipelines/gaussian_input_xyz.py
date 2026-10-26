from pathlib import Path

from ..io import read_xyz, write_gaussian_input
from ._templates import read_gaussian_template


def xyzs_to_gaussian_inputs(
    folder: str | Path,
    basename_xyz: str,
    path_template: str | Path,
    *,
    basename_com: str = "gaussian_sp.com",
    overwrite: bool = False
) -> int:
    """
    folder の下の全ての basename_xyz から、同じフォルダの basename_com に Gaussian の入力を書く。

    - xyz のタイトル行に charge と mult が必要（write_xyz の record=True）
    - link0（%chk を除く）, route, title は path_template の Gaussian の入力から写す
    - 書いた入力の数を返す
    """
    template = read_gaussian_template(path_template)

    paths = sorted(Path(folder).rglob(basename_xyz))
    for path in paths:
        write_gaussian_input(
            read_xyz(path),
            path.with_name(basename_com),
            with_notes=False,
            overwrite=overwrite,
            **template
        )
    return len(paths)
