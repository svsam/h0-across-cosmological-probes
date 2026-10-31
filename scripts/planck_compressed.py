"""Conservative, diagonal Planck 2018 CMB summary likelihood.

This is a compressed approximation, not the Planck power-spectrum likelihood.
Its three inputs are 100 theta_star, omega_b h^2, and omega_c h^2.
"""

from cosmosis.datablock import names


# The acoustic-scale mean is Planck TT,TE,EE+lowE (2018), eq. (9) of
# Planck Collaboration VI. The physical-density centers are from the same
# paper's abstract. Errors are deliberately broadened to account partly for
# model extension and omitted off-diagonal covariance; they are choices of
# this approximate likelihood, not a published joint covariance.
CENTERS = (1.04109, 0.0224, 0.120)
SIGMAS = (0.0006, 0.0003, 0.003)


def setup(options):
    return True


def execute(block, config):
    cosmo = names.cosmological_parameters
    theta = block[names.distances, "thetastar"]
    h = block[cosmo, "h0"]
    omega_b = block[cosmo, "omega_b"]
    omega_m = block[cosmo, "omega_m"]
    values = (theta, omega_b * h * h, (omega_m - omega_b) * h * h)
    chi2 = sum(((value - center) / sigma) ** 2
               for value, center, sigma in zip(values, CENTERS, SIGMAS))
    block[names.likelihoods, "planck_compressed_like"] = -0.5 * chi2
    block[names.data_vector, "planck_compressed_chi2"] = chi2
    return 0


def cleanup(config):
    pass
