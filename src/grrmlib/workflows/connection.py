from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ..core import Molecules
from ..io import key_to_folder, read_connectable, read_xyz, write_xyzs
from ..operations import connect_all_iter, copy_notes, step_to_angles
from ..pipelines import SEPARATION_KEYS, connection_keys


def connection_key_names(substituents: Iterable[str]) -> list[str]:
    """
    置換基を導入した構造のフォルダの key の名前。
    例：["SG", "SEQ", "charge", "mult", "arene", "arene_angle", "alkene", "alkene_angle"]
    """
    return SEPARATION_KEYS + connection_keys("", list(substituents))[1:]


def _read_substituents(folder: Path) -> Molecules:
    # folder/0.com, folder/1.com, ... を key が (置換基の番号, 角度) の Molecules にする
    paths = sorted(Path(folder).glob("*.com"), key=lambda p: int(p.stem))
    mols = Molecules({int(p.stem): read_connectable(p) for p in paths})
    return mols.expand(step_to_angles).flatten()


def make_connection_sp(
    path_selected: str | Path,
    folder_opt: str | Path,
    basename_opt: str,
    path_core: str | Path,
    folders_sub: Iterable[str | Path],
    folder: str | Path,
    *,
    basename: str = "init.xyz"
) -> int:
    """
    path_selected の (SG, SEQ, charge, mult) の最適化構造（folder_opt の下の basename_opt）に置換基をつなぎ、
    回転させた配座を folder の下に書く。

    - path_core：EQ と同じ原子番号（labels）を持つ構造に、置換基を付け替える位置を notes で指定したもの（note.com）。
      SG の原子には labels で notes を写す（copy_notes）
    - folders_sub：置換基のテンプレート（<番号>.com）のフォルダ。フォルダ名が置換基の名前になる
    - 結果：folder/SG=*/SEQ=*/charge=*/mult=*/<置換基>=*/<置換基>_angle=*/.../basename。
      SG に含まれない置換基と角度は None。タイトル行に charge, mult, labels（つないだ後の 1, 2, ...）を書く
    - 初期構造は入力から決まるので、書き直しても中身は変わらない（overwrite=True）。書いた数を返す
    """
    folders_sub = [Path(f) for f in folders_sub]
    key_names = connection_key_names(f.name for f in folders_sub)
    core = read_connectable(path_core)
    subs_list = [_read_substituents(f) for f in folders_sub]
    selected = pl.read_csv(path_selected)

    n = 0
    for key in selected.select(SEPARATION_KEYS).iter_rows():
        mol = read_xyz(key_to_folder(folder_opt, key, SEPARATION_KEYS) / basename_opt)
        mol = copy_notes(mol, core)

        mols = Molecules()
        for names, mol_new in connect_all_iter(Molecules({key: mol}), subs_list):
            # 角度は int にする（connect_all と同じ）
            subs = tuple(int(x) if x is not None else None for x in names[1:])
            mols[(*names[0], *subs)] = mol_new
        write_xyzs(mols, folder, basename, prefixes=key_names, record=True, overwrite=True)
        n += len(mols)

    print(f"connection inputs: {n} written from {selected.height} (SG, charge, mult)")
    return n
