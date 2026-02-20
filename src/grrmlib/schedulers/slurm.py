from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class SlurmScript:
    job_name: str
    partition: str
    time: str                  # "02:00:00"
    cpus_per_task: int
    nodes: int = 1
    ntasks: int = 1
    account: str | None = None
    qos: str | None = None
    constraint: str | None = None
    gres: str | None = None
    output: str | None = None
    error: str | None = None
    modules: list[str] = field(default_factory=list)
    commands: list[str] = field(default_factory=list)

    def add_command(self, cmd: str) -> None:
        self.commands.append(cmd)

    def add_commands(self, cmds: Iterable[str]) -> None:
        self.commands.extend(cmds)

    def render(self) -> str:
        lines = [
            "#!/bin/bash",
            f"#SBATCH -J {self.job_name}",
            f"#SBATCH -p {self.partition}",
            f"#SBATCH -t {self.time}",
            f"#SBATCH -N {self.nodes}",
            f"#SBATCH -n {self.ntasks}",
            f"#SBATCH -c {self.cpus_per_task}",
        ]

        if self.account:
            lines.append(f"#SBATCH -A {self.account}")
        if self.qos:
            lines.append(f"#SBATCH --qos={self.qos}")
        if self.constraint:
            lines.append(f"#SBATCH --constraint={self.constraint}")
        if self.gres:
            lines.append(f"#SBATCH --gres={self.gres}")
        if self.output:
            lines.append(f"#SBATCH -o {self.output}")
        if self.error:
            lines.append(f"#SBATCH -e {self.error}")

        lines.append("")
        lines.append("set -e")
        lines.append("")

        for m in self.modules:
            lines.append(f"module load {m}")

        if self.modules:
            lines.append("")

        lines.extend(self.commands)

        return "\n".join(lines)

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(self.render())
        return path

    def sbatch(self, path: str | Path) -> None:
        import subprocess
        subprocess.run(["sbatch", str(path)], check=True)


@dataclass
class SlurmArrayScript(SlurmScript):
    array: str | None = None        # "0-999"
    array_parallelism: int | None = None  # 20

    def render(self) -> str:
        lines = super().render().split("\n")

        insert_index = 1  # #!/bin/bash の次

        if self.array:
            if self.array_parallelism:
                lines.insert(
                    insert_index,
                    f"#SBATCH --array={self.array}%{self.array_parallelism}"
                )
            else:
                lines.insert(
                    insert_index,
                    f"#SBATCH --array={self.array}"
                )

        return "\n".join(lines)