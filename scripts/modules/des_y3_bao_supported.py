"""DES-Y3 BAO likelihood with explicit table-support handling."""

from pathlib import Path

import numpy as np
from cosmosis.datablock import option_section


DATA = (
    Path(__file__).resolve().parents[3]
    / "likelihood/des-y3-bao/chi2profile_dvdesy3_cosmoplanck18_covcosmolike.csv"
)

# Fiducial distances used to define the DES-Y3 BAO alpha measurement.
REDSHIFT = 0.835
D_A_FIDUCIAL_MPC = 1616.9
R_S_FIDUCIAL_MPC = 147.6


def setup(options):
    strict = options.get_bool(option_section, "strict_support", default=True)
    chi2_column = options.get_int(option_section, "chi2_column", default=1)
    rd_scale = options.get_double(option_section, "rd_scale", default=1.0)
    if not np.isfinite(rd_scale) or rd_scale <= 0:
        raise ValueError("rd_scale must be a finite positive calibration factor")
    data = np.loadtxt(DATA, delimiter=",", skiprows=1)
    if chi2_column < 1 or chi2_column >= data.shape[1]:
        raise ValueError("Invalid DES-Y3 BAO chi-square column")
    alpha, chi2 = data[:, 0], data[:, chi2_column]
    if np.any(np.diff(alpha) <= 0):
        raise ValueError("DES-Y3 BAO alpha grid must increase strictly")
    return alpha, chi2, strict, rd_scale


def execute(block, config):
    grid, chi2_grid, strict, rd_scale = config
    z = block["distances", "z"]
    d_m = block["distances", "D_M"]
    r_s = block["distances", "rs_zdrag"]
    if z[1] < z[0]:
        z, d_m = z[::-1], d_m[::-1]
    if not z[0] <= REDSHIFT <= z[-1]:
        raise ValueError(f"Theory distance grid must cover DES-Y3 BAO redshift {REDSHIFT}")
    alpha = float(
        np.interp(REDSHIFT, z, d_m)
        / (D_A_FIDUCIAL_MPC * (1 + REDSHIFT))
        * R_S_FIDUCIAL_MPC
        / (r_s * rd_scale)
    )
    supported = bool(np.isfinite(alpha) and grid[0] <= alpha <= grid[-1])
    block["des_y3_bao", "alpha"] = alpha
    block["des_y3_bao", "in_support"] = float(supported)
    if not supported and strict:
        block["likelihoods", "des_y3_bao_like"] = -np.inf
        return 0
    chi2 = float(np.interp(alpha, grid, chi2_grid))
    block["des_y3_bao", "chi2"] = chi2
    block["likelihoods", "des_y3_bao_like"] = -0.5 * chi2
    return 0
