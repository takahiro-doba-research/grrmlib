import re
from pathlib import Path

import numpy as np

from ..core import Molecule, Molecules, Trajectory


class GRRMMINOutputReader:
    
    def __init__(self, version: str = "GRRM23") -> None:
        self.version = version
    
    def read_mols(
        self,
        folder: str | Path,
        basename: str = "grrm.log"
    ) -> Molecules:
        paths = sorted(Path(folder).rglob(basename))
        mols = Molecules()
        
        for path in paths:
            key = tuple(p.split("=")[1] for p in path.parts[1:-1])
            key = tuple(None if k == "None" else int(k) for k in key)
            try:
                _, opt, _ = self.read(path)
                mols[key] = opt
            except Exception as e:
                mols[key] = Molecule()
                print(key, e)
        
        return mols
    
    def read(self, path: str | Path) -> tuple[Trajectory, Molecule, Molecule | None]:
        path = Path(path)
        text = path.read_text()
        lines = text.splitlines()
        return self.parse(lines)
    
    def parse(self, lines: list[str]) -> tuple[Trajectory, Molecule, Molecule | None]:
        indices_optopt = [i for i, l in enumerate(lines) if l.startswith("OPTOPT")]
        lines_optopt = lines[indices_optopt[0] + 1: indices_optopt[1]]
        itrs, opt = self._parse_optopt(lines_optopt)
        
        indices_freqfreq = [i for i, l in enumerate(lines) if l.startswith("FREQFREQ")]
        if indices_freqfreq:
            lines_freqfreq = lines[indices_freqfreq[0] + 1: indices_freqfreq[1]]
            if self.version == "GRRM23":
                freq = self._parse_freqfreq(lines_freqfreq)
            elif self.version == "GRRMdev":
                freq = self._parse_freqfreq_dev(lines_freqfreq)
            else:
                freq = self._parse_freqfreq(lines_freqfreq)
        else:
            freq = None
        
        return itrs, opt, freq
    
    def _parse_optopt(self, lines: list[str]) -> tuple[Trajectory, Molecule]:
        indices_eqeq = [i for i, l in enumerate(lines) if l.startswith("======")]
        
        if indices_eqeq:
            index_eqeq = indices_eqeq[0]
            lines_itrs = lines[0:index_eqeq]
            itrs = self._parse_itrs(lines_itrs)
            
            lines_opt = lines[index_eqeq + 1:]
            opt = self._parse_opt(lines_opt)
            opt.status = lines[-1].strip()
        
        else:
            lines_itrs = lines[0:-1]
            itrs = self._parse_itrs(lines_itrs)
            
            opt = Molecule()
            opt.status = lines[-1].strip()
        
        return itrs, opt
    
    def _parse_itrs(self, lines: list[str]) -> Trajectory:
        indices_itr = [i for i, l in enumerate(lines) if l.startswith("#")]
        indices_itr.append(len(lines))
        traj = Trajectory()
        
        for i, j in zip(indices_itr, indices_itr[1:]):
            lines_itr = lines[i:j]
            itr = self._parse_itr(lines_itr)
            traj.append(itr)
        
        return traj
    
    def _parse_itr(self, lines: list[str]) -> Molecule:
        name = int(re.search(r"ITR. (\d+)", lines[0]).group(1))
        index_item = [i for i, l in enumerate(lines) if "Item" in l][0]
        lines_atomcoords = lines[1:index_item]
        labels, symbols, atomcoords = self._parse_atomcoords(lines_atomcoords)
        scfenergy = float(re.search(r"\(\s*(-?\d+\.?\d+)\s*:", lines[index_item + 1]).group(1))
        afirenergy = float(re.search(r"ENERGY\s*(-?\d+\.?\d+)\s*\(", lines[index_item + 1]).group(1))
        mult = float(re.search(r"Spin\(\*\*2\)\s*(-?\d+\.?\d+)", lines[index_item + 2]).group(1)) * 2 + 1
        index_nmeigen = [i for i, l in enumerate(lines) if l.startswith("NORMAL MODE EIGENVALUE")][0]
        lines_nmeigen = lines[index_nmeigen + 1:]
        nmeigen = np.concatenate([list(map(float, l.split())) for l in lines_nmeigen])
        return Molecule(
            name=name,
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            afirenergy=afirenergy,
            nmeigen=nmeigen,
        )
    
    def _parse_opt(self, lines: list[str]) -> Molecule:
        index_energy = [i for i, l in enumerate(lines) if l.startswith("ENERGY")][0]
        lines_atomcoords = lines[1:index_energy]
        labels, symbols, atomcoords = self._parse_atomcoords(lines_atomcoords)
        scfenergy = float(re.search(r"=\s*(-?\d+\.?\d+)\s*", lines[index_energy]).group(1))
        mult = float(re.search(r"=\s*(-?\d+\.?\d+)", lines[index_energy + 1]).group(1)) * 2 + 1
        zpve = float(m.group(1)) if (m := re.search(r"=\s*(-?\d+\.?\d+)", lines[index_energy + 2])) else None
        index_grad = [i for i, l in enumerate(lines) if l.startswith("GRADIENT VECTOR")][0]
        index_hess = [i for i, l in enumerate(lines) if l.startswith("HESSIAN MATRIX")][0]
        lines_grad = lines[index_grad + 1:index_hess]
        grads = np.array([float(l.split()[0]) for l in lines_grad])
        
        index_nmeigen = [i for i, l in enumerate(lines) if l.startswith("NORMAL MODE EIGENVALUE")][0]
        lines_hess = lines[index_hess + 1: index_nmeigen]
        hessian = self._parse_hessian(lines_hess)
        
        lines_nmeigen = lines[index_nmeigen + 1:-2]
        nmeigen = np.concatenate([list(map(float, l.split())) for l in lines_nmeigen])
        
        return Molecule(
            mult=mult,
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
            zpve=zpve,
            grads=grads,
            hessian=hessian,
            nmeigen=nmeigen,
        )
    
    def _parse_hessian(self, lines: list[str]) -> np.ndarray:
        indices_chunk = [i for i, l in enumerate(lines) if len(l) <= 15]
        indices_chunk.append(len(lines))
        hesses = []
        
        for i, j in zip(indices_chunk, indices_chunk[1:]):
            lines_chunk = lines[i:j]
            lines_chunk = [list(map(float, l.split())) for l in lines_chunk]
            
            for k in range(len(lines_chunk[-1])):
                hess = [l[k] for l in lines_chunk[k:]]
                hesses.append(hess)
        
        hess_1d = np.concatenate(hesses)
        dim = len(hesses[0])
        hess_2d = np.full((dim, dim), np.nan)
        hess_2d[np.tril_indices(dim)] = hess_1d[::-1]
        hess_2d = hess_2d[::-1, ::-1].T
        return hess_2d
    
    def _parse_freqfreq(self, lines: list[str]) -> Molecule:
        indices_freq = [i for i, l in enumerate(lines) if l.startswith("Freq.")]
        lines_atomcoords = lines[1:indices_freq[0] - 2]
        labels, symbols, atomcoords = self._parse_atomcoords(lines_atomcoords)
        vibfreqs = np.array([float(freq) for i in indices_freq for freq in lines[i].split()[2:5]])
        
        index_thermo_all = [i for i, l in enumerate(lines) if l.startswith("Thermochemistry (use all")][0]
        lines_thermo_all = lines[index_thermo_all + 1: index_thermo_all + 15]
        (
            scfenergy,
            zpve,
            scfzpvenergy,
            transenergy,
            rotenergy,
            vibenergy,
            enthalpy_correction,
            enthalpy,
            scfentropy,
            transentropy,
            rotentropy,
            vibentropy,
            freeenergy_correction,
            freeenergy,
        ) = [float(re.search(r"=\s*(-?\d+\.?\d+)", l).group(1)) for l in lines_thermo_all]
        
        index_thermo_repl = [i for i, l in enumerate(lines) if l.startswith("Thermochemistry (after")][0]
        lines_thermo_repl = lines[index_thermo_repl + 1: index_thermo_repl + 15]
        (
            scfenergy_repl,
            zpve_repl,
            scfzpvenergy_repl,
            transenergy_repl,
            rotenergy_repl,
            vibenergy_repl,
            enthalpy_correction_repl,
            enthalpy_repl,
            scfentropy_repl,
            transentropy_repl,
            rotentropy_repl,
            vibentropy_repl,
            freeenergy_correction_repl,
            freeenergy_repl,
        ) = [float(re.search(r"=\s*(-?\d+\.?\d+)", l).group(1)) for l in lines_thermo_repl]
        
        return Molecule(
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            vibfreqs=vibfreqs,
            scfenergy=scfenergy,
            zpve=zpve,
            scfzpvenergy=scfzpvenergy,
            transenergy=transenergy,
            rotenergy=rotenergy,
            vibenergy=vibenergy,
            enthalpy_correction=enthalpy_correction,
            enthalpy=enthalpy,
            scfentropy=scfentropy,
            transentropy=transentropy,
            rotentropy=rotentropy,
            vibentropy=vibentropy,
            freeenergy_correction=freeenergy_correction,
            freeenergy=freeenergy,
            scfenergy_repl=scfenergy_repl,
            zpve_repl=zpve_repl,
            scfzpvenergy_repl=scfzpvenergy_repl,
            transenergy_repl=transenergy_repl,
            rotenergy_repl=rotenergy_repl,
            vibenergy_repl=vibenergy_repl,
            enthalpy_correction_repl=enthalpy_correction_repl,
            enthalpy_repl=enthalpy_repl,
            scfentropy_repl=scfentropy_repl,
            transentropy_repl=transentropy_repl,
            rotentropy_repl=rotentropy_repl,
            vibentropy_repl=vibentropy_repl,
            freeenergy_correction_repl=freeenergy_correction_repl,
            freeenergy_repl=freeenergy_repl,
        )
    
    def _parse_freqfreq_dev(self, lines: list[str]) -> Molecule:
        indices_freq = [i for i, l in enumerate(lines) if l.startswith("Freq.")]
        lines_atomcoords = lines[1:indices_freq[0] - 2]
        labels, symbols, atomcoords = self._parse_atomcoords(lines_atomcoords)
        vibfreqs = np.array([float(freq) for i in indices_freq for freq in lines[i].split()[2:5]])
        
        index_thermo_repl = [i for i, l in enumerate(lines) if l.startswith("Thermochemistry at")][0]
        lines_thermo_repl = lines[index_thermo_repl + 1: index_thermo_repl + 15]
        (
            scfenergy_repl,
            zpve_repl,
            scfzpvenergy_repl,
            transenergy_repl,
            rotenergy_repl,
            vibenergy_repl,
            enthalpy_correction_repl,
            enthalpy_repl,
            scfentropy_repl,
            transentropy_repl,
            rotentropy_repl,
            vibentropy_repl,
            freeenergy_correction_repl,
            freeenergy_repl,
        ) = [float(re.search(r"=\s*(-?\d+\.?\d+)", l).group(1)) for l in lines_thermo_repl]
        
        return Molecule(
            labels=labels,
            symbols=symbols,
            atomcoords=atomcoords,
            vibfreqs=vibfreqs,
            scfenergy_repl=scfenergy_repl,
            zpve_repl=zpve_repl,
            scfzpvenergy_repl=scfzpvenergy_repl,
            transenergy_repl=transenergy_repl,
            rotenergy_repl=rotenergy_repl,
            vibenergy_repl=vibenergy_repl,
            enthalpy_correction_repl=enthalpy_correction_repl,
            enthalpy_repl=enthalpy_repl,
            scfentropy_repl=scfentropy_repl,
            transentropy_repl=transentropy_repl,
            rotentropy_repl=rotentropy_repl,
            vibentropy_repl=vibentropy_repl,
            freeenergy_correction_repl=freeenergy_correction_repl,
            freeenergy_repl=freeenergy_repl,
        )
    
    def _parse_atomcoords(self, lines: list[str]) -> tuple[np.ndarray, list[str], np.ndarray]:
        labels = np.arange(1, len(lines) + 1)
        symbols = [l.split()[0] for l in lines]
        atomcoords = np.array([list(map(float, l.split()[1:4])) for l in lines])
        return labels, symbols, atomcoords