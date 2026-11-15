# CO2 binding at an open Cu site in HKUST-1 with MACE-MP-0

This project calculates the binding energy of one CO2 molecule at an open
(coordinatively unsaturated) Cu site in HKUST-1 (Cu-BTC). It uses the
**MACE-MP-0 medium** foundation model plus **D3 dispersion**, run through ASE
on the GPU (CUDA, float64 for all optimisations).

Binding energy is defined as

    E_bind = E(MOF + CO2) - E(MOF) - E(CO2, gas)

The structures are relaxed with BFGS to fmax = 0.02 eV/Å (0.01 eV/Å for gas-phase CO2).

## Environment

- Conda env `ml` (mace-torch 0.3.17, ASE 3.29.0)
- Run with `conda run -n ml python <script>.py`

## Scripts (in run order)

| Script | What it does |
|---|---|
| `common.py` | Shared helper that returns the MACE-MP-0 medium calculator (CUDA, D3 optional). |
| `01_test.py` | Reads `hkust1.cif`, prints the cell, and times one energy/force call. |
| `02_relax_mof.py` | Relaxes the empty framework's atoms and cell together (FrechetCellFilter), writing `mof_relaxed.extxyz`. |
| `03_co2_gas.py` | Relaxes an isolated CO2 molecule, writing `co2_relaxed.extxyz`. |
| `04_place_co2.py` | Puts CO2 end-on at a Cu open site with a tilted start, relaxes it with the cell fixed, and writes `complex_tilt<N>_relaxed.extxyz`. |
| `05_binding.py` | Computes the binding energy with and without D3. |

## Files

- `hkust1.cif`: starting structure (Fm-3m, a = 26.343 Å, 624 atoms)
- `*_relaxed.extxyz`: final relaxed structures (MOF, CO2, tilted complex)
- `*.log`: BFGS optimiser logs
- `co2_relax.gif`: animation of the CO2 relaxation at the Cu site

The trajectories (`*.traj`) and the extxyz files converted from them
(`*_relax.extxyz`) are not tracked. Rerun the scripts to regenerate them.

## Results

**Relaxed framework** (`mof_relaxed.extxyz`): cubic, a = 26.464 Å (the CIF
starts at 26.343 Å, so the cell expands by 0.5%). E = −4463.3201 eV.

**Gas-phase CO2:** E = −22.8050 eV.

**Straight end-on start** (Cu–O–C = 180°; files not tracked): stayed at
180°, E = −4486.2068 eV. This turned out to be a **saddle point**: the
symmetric start has no force pushing it to bend, so the optimiser never left it.

**Tilted start** (30° off the Cu–Cu axis, `complex_tilt30_relaxed.extxyz`):
relaxed to a bent end-on geometry that is 1.4 kJ/mol lower in energy:

| Quantity | Value |
|---|---|
| Cu–O | 2.496 Å |
| Cu–O–C | 144.4° |
| O–C–O | 179.6° |
| E(MOF + CO2) | −4486.2213 eV |

**Binding energy** (tilted minimum):

| Model | E_bind |
|---|---|
| MACE-MP-0 + D3 | **−9.3 kJ/mol** |
| MACE-MP-0 without D3 | +0.7 kJ/mol |

### Interpretation

All of the binding comes from the D3 dispersion. Without D3 the interaction
is slightly repulsive. The result is about 2–3× weaker than typical
experimental heats of adsorption for CO2 in HKUST-1. MACE-MP-0 does not
capture the electrostatic/orbital interaction at the open Cu site, so it
treats the Cu site like any other pore surface.

## Next steps

- Compare against the GCMC isosteric heat (see `../cubtc_isotherm`) and literature values
- Try a model trained on ODAC (Open DAC) MOF–adsorbate data
- Run DFT on a Cu paddlewheel cluster to benchmark the Cu···OCO interaction
- Look at other binding sites (small tetrahedral cages, cage windows)
