# 02_relax_mof.py — relax atoms and cell of empty HKUST-1
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from common import get_calc

atoms = read("hkust1.cif")
atoms.calc = get_calc()
print(f"Start: a = {atoms.cell.cellpar()[0]:.3f} Å, E = {atoms.get_potential_energy():.4f} eV")

# FrechetCellFilter lets the optimizer change the cell as well as positions
opt = BFGS(FrechetCellFilter(atoms), trajectory="mof_relax.traj", logfile="mof_relax.log")
opt.run(fmax=0.02, steps=1000)            # stop when every force < 0.02 eV/Å

print(f"Converged: {opt.converged()}  after {opt.nsteps} steps")
print(f"Final: a = {atoms.cell.cellpar()[0]:.3f} Å, E = {atoms.get_potential_energy():.4f} eV")
write("mof_relaxed.extxyz", atoms)