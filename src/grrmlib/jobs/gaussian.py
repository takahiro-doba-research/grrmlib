import subprocess
from pathlib import Path


class GaussianJob:
    
    def __init__(self, path_input: str | Path) -> None:
        self.path_input = Path(path_input)
        self.cwd = self.path_input.parent
    
    def run(self) -> None:
        subprocess.run(
            ["g16", self.path_input.name],
            cwd=self.cwd,
        )