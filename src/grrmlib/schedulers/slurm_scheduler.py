import re
import subprocess
import time
from pathlib import Path


class KyotoUnivSlurmScheduler:
    
    def __init__(
        self,
        job_name: str,
        output: str,
        *,
        partition: str,
        time: str = "14-00:00:00",
        ntasks: int = 1,
        thread_spec: int = 8,
        cpus_per_task: int = 8,
        mem: str = "8G",
        array: str | None = None,
        commands: list[str] | None = None,
    ) -> None:
        self.job_name = job_name
        self.output = output
        self.partition = partition
        self.time = time
        self.ntasks = ntasks
        self.thread_spec = thread_spec
        self.cpus_per_task = cpus_per_task
        self.mem = mem
        self.array = array
        self.commands = commands
    
    def _build_rsc(self) -> str:
        # Build --rsc for Kyoto University Slurm
        return (
            "#SBATCH --rsc "
            f"p={self.ntasks}"
            f":t={self.thread_spec}"
            f":c={self.cpus_per_task}"
            f":m={self.mem}"
        )
        
    def build(self) -> str:
        lines = [
            "#!/bin/bash",
            f"#SBATCH -J {self.job_name}",
            f"#SBATCH -o {self.output}",
            f"#SBATCH -p {self.partition}",
            f"#SBATCH -t {self.time}",
        ]
        
        lines.append(self._build_rsc())
        
        if self.array is not None:
            lines.append(f"#SBATCH --array={self.array}")
        
        if self.commands is not None:
            lines.append("")
            lines.extend(self.commands)
        
        return "\n".join(lines)
    
    def write(
        self,
        path: str | Path,
        *,
        overwrite: bool = False
    ) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        text = self.build()
        
        with path.open("w" if overwrite else "x") as f:
            f.write(text)
        
        return path
    
    def sbatch(self, path: str | Path) -> str:
        result = subprocess.run(
            ["sbatch", str(path)],
            capture_output=True,
            text=True,
            check=True,
        )
        match = re.search(r"\d+", result.stdout)
        
        if not match:
            raise RuntimeError("Failed to parse job ID")
            
        return match.group(0)
    
    def status(self, job_id: str) -> str:
        result = subprocess.run(
            ["sacct", "-j", str(job_id), "--format=JobID,State,Elapsed"],
            capture_output=True,
            text=True,
        )
        return result.stdout
    
    def cancel(self, job_id: str) -> None:
        subprocess.run(["scancel", str(job_id)], check=True)

    def wait(self, job_id: str, interval: int = 30) -> None:
        while True:
            status = self.status(job_id)
            
            if "COMPLETED" in status or "FAILED" in status or "CANCELLED" in status:
                return
            
            time.sleep(interval)