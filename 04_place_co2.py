# 04_place_co2.py — put CO2 at an open Cu site (tilted start) and relax
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from common import get_calc

TILT_DEG = 30          # angle between the CO2 axis and the Cu-Cu line at the start
CU_O_START = 2.40      # starting Cu...O distance (Å)
TAG = f"tilt{TILT_DEG}"

mof = read("mof_relaxed.extxyz")
cu = [a.index for a in mof if a.symbol == "Cu"]
i = cu[0]                                    # the Cu we will bind to

# Find the paddlewheel partner: nearest other Cu, using periodic images
vecs = mof.get_distances(i, cu, mic=True, vector=True)
dists = np.linalg.norm(vecs, axis=1)
dists[cu.index(i)] = np.inf
k = int(np.argmin(dists))
print(f"Cu {i} partner: Cu {cu[k]} at {dists[k]:.3f} Å")

# u: unit vector from the partner Cu through Cu i, pointing into the pore (open site)
u = -vecs[k] / dists[k]

# w: any unit vector perpendicular to u, used to tilt the CO2
w = np.cross(u, [1.0, 0.0, 0.0])
if np.linalg.norm(w) < 0.1:
    w = np.cross(u, [0.0, 1.0, 0.0])
w /= np.linalg.norm(w)

t = np.radians(TILT_DEG)
axis = np.cos(t) * u + np.sin(t) * w         # CO2 molecular axis direction

p = mof.positions[i]
O1 = p + CU_O_START * u                      # binding O, straight out from Cu
C = O1 + 1.17 * axis
O2 = C + 1.17 * axis
co2 = Atoms("OCO", positions=[O1, C, O2])

complex_ = mof + co2                         # CO2 atoms are the last three
n = len(complex_)
iO1, iC, iO2 = n - 3, n - 2, n - 1

# Clash check with periodic images: each CO2 atom vs all framework atoms
for name, idx in (("O1", iO1), ("C", iC), ("O2", iO2)):
    d = complex_.get_distances(idx, list(range(len(mof))), mic=True)
    print(f"{name}: nearest framework atom {d.min():.2f} Å")

print(f"Start: Cu-O-C = {complex_.get_angle(i, iO1, iC, mic=True):.1f}°")

complex_.calc = get_calc()
opt = BFGS(complex_, trajectory=f"complex_{TAG}_relax.traj",
           logfile=f"complex_{TAG}_relax.log")
opt.run(fmax=0.02, steps=1000)               # cell stays fixed

print(f"Converged: {opt.converged()}  after {opt.nsteps} steps")
print(f"Cu-O distance: {complex_.get_distance(i, iO1, mic=True):.3f} Å")
print(f"Cu-O-C angle: {complex_.get_angle(i, iO1, iC, mic=True):.1f}°")
print(f"O-C-O angle:  {complex_.get_angle(iO1, iC, iO2, mic=True):.1f}°")
print(f"E(MOF+CO2) = {complex_.get_potential_energy():.4f} eV   (straight start was -4486.2068 eV)")
write(f"complex_{TAG}_relaxed.extxyz", complex_)