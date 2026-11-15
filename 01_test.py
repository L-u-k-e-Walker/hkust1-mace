# 01_test.py — inspect the structure and time one MACE calculation
import time
import numpy as np
from ase.io import read
from common import get_calc

atoms = read("hkust1.cif")
print("Atoms:", len(atoms), "  formula:", atoms.get_chemical_formula())
print("Cell a, b, c (Å), angles (°):", atoms.cell.cellpar().round(3))
print("Periodic:", atoms.pbc)

atoms.calc = get_calc()
t = time.time()
e = atoms.get_potential_energy()          # eV, whole cell
f = atoms.get_forces()                     # eV/Å, one vector per atom
print(f"Energy: {e:.4f} eV   ({time.time() - t:.1f} s, includes model loading)")
print(f"Largest force: {np.abs(f).max():.3f} eV/Å")

atoms.rattle(0.001)                        # tiny nudge so it recalculates
t = time.time()
atoms.get_potential_energy()
print(f"Second calculation: {time.time() - t:.2f} s")