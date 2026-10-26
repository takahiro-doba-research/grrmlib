from pathlib import Path
from typing import Any

from ..io import read_gaussian_input, read_grrm_input


def read_gaussian_template(path: str | Path) -> dict[str, Any]:
    """
    Gaussian の入力ファイルから write_gaussian_input(s) に渡す link0, route, title を読む。
    """
    mol = read_gaussian_input(path)
    return {
        # chk は各ジョブで要らないので引き継がない
        "link0": [l for l in mol.link0 if not l.lower().startswith("%chk")],
        "route": mol.route,
        "title": mol.title,
    }


def read_grrm_template(path: str | Path) -> dict[str, Any]:
    """
    GRRM の入力ファイルから write_grrm_input(s) に渡す route, options を読む。
    """
    mol = read_grrm_input(path)
    return {
        "route": mol.route,
        "options": mol.options,
    }
