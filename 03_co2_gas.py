# 03_co2_gas.py — relax an isolated CO2 molecule
from ase.build import molecule
from ase.io import write
from ase.optimize import BFGS
from common import get_calc

co2 = molecule("CO2")                      # ASE's built-in geometry
co2.center(vacuum=10.0)
co2.pbc = False                            # a true isolated molecule
co2.calc = get_calc()

BFGS(co2, logfile="co2_relax.log").run(fmax=0.01)

print("Symbols:", co2.get_chemical_symbols())
print("Distances (Å):\n", co2.get_all_distances().round(3))
print("O-C-O angle (°):", round(co2.get_angle(1, 0, 2), 1))
print(f"E(CO2) = {co2.get_potential_energy():.4f} eV")
write("co2_relaxed.extxyz", co2)