# common.py — shared MACE-MP calculator settings
from mace.calculators import mace_mp

def get_calc(dtype="float64", dispersion=True):
    """MACE-MP foundation model on the GPU.
    float64 for geometry optimization (precise forces);
    float32 for MD (about twice as fast)."""
    return mace_mp(model="medium", device="cuda",
                   default_dtype=dtype, dispersion=dispersion)