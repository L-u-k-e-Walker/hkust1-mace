# 05_binding.py — binding energy with and without D3 dispersion
from ase.io import read
from common import get_calc

EV_TO_KJMOL = 96.485
structs = {
    "mof":  read("mof_relaxed.extxyz"),
    "co2":  read("co2_relaxed.extxyz"),
    "cplx": read("complex_tilt30_relaxed.extxyz"),
}

for disp in (True, False):
    calc = get_calc(dispersion=disp)
    E = {}
    for name, atoms in structs.items():
        atoms.calc = calc
        E[name] = atoms.get_potential_energy()
    eb = E["cplx"] - E["mof"] - E["co2"]
    label = "MACE-MP + D3" if disp else "MACE-MP only"
    print(f"{label:14s}  E_bind = {eb:+.3f} eV = {eb * EV_TO_KJMOL:+.1f} kJ/mol")