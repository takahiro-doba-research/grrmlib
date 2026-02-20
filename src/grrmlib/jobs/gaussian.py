import subprocess
from pathlib import Path


class GaussianJob:
    
    def __init__(
        self,
        path_input: str | Path,
        *,
        nprocshared: int | None = None,
        mem: str | None = None,
        chk: str | None = None,
    ) -> None:
        self.path_input = Path(path_input)
        self.path_output = self.path_input.with_suffix(".log")
        self.cwd = self.path_input.parent
        self.nprocshared = nprocshared
        self.mem = mem
        self.chk = chk
    
    def run(self) -> None:
        with (
            open(self.path_input, "r") as f_in,
            open(self.path_output, "w") as f_out
        ):
            subprocess.run(
                ["g16"],
                stdin=f_in,
                stdout=f_out,
                cwd=self.cwd,
                check=True,
            )