import subprocess
from pathlib import Path


class GaussianJob:
    
    def __init__(
        self,
        path_input: str | Path,
        command: str = "g16"
    ) -> None:
        self.path_input = Path(path_input)
        self.command = command
        self.cwd = self.path_input.parent
    
    def run(self) -> None:
        subprocess.run(
            [self.command, self.path_input.name],
            cwd=self.cwd,
        )
    
    def schedule(self, scheduler) -> None:
        scheduler.submit(self)
    
    def status(self) -> str:
        pass