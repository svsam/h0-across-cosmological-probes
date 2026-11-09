"""Curved wCDM Astropy distances with explicit standard radiation closure.

This project-local module keeps the shared CosmoSIS background module untouched.
Omega_DE is derived from Omega_m + Omega_r + Omega_DE + Omega_k = 1, with
photons and massless neutrinos calculated from the fixed standard settings.
"""

import astropy.constants
import astropy.cosmology
import astropy.units as u
import numpy as np
from cosmosis.datablock import names, option_section


def setup(options):
    return {
        "zmax": options.get_double(option_section, "zmax", default=3.0),
        "nz": options.get_int(option_section, "nz", default=500),
        "tcmb": options.get_double(option_section, "tcmb", default=2.7255),
        "neff": options.get_double(option_section, "neff", default=3.046),
    }


def build_cosmology(hubble, omega_m, omega_k, w, tcmb, neff):
    """Build the curved model with radiation included in density closure."""
    hubble = hubble * u.km / u.s / u.Mpc
    tcmb = tcmb * u.K
    m_nu = np.zeros(3) * u.eV

    # Radiation density does not depend on Ode0. Build a probe cosmology to
    # calculate it, then derive dark energy from exact present-day closure.
    probe = astropy.cosmology.wCDM(
        H0=hubble,
        Om0=omega_m,
        Ode0=0.0,
        w0=w,
        Tcmb0=tcmb,
        Neff=neff,
        m_nu=m_nu,
    )
    omega_r = float(probe.Ogamma0 + probe.Onu0)
    omega_de = 1.0 - omega_m - omega_r - omega_k
    if omega_m <= 0 or omega_de <= 0 or not np.isfinite(omega_de):
        raise ValueError("Unphysical density closure: Omega_m and Omega_DE must be positive")

    model = astropy.cosmology.wCDM(
        H0=hubble,
        Om0=omega_m,
        Ode0=omega_de,
        w0=w,
        Tcmb0=tcmb,
        Neff=neff,
        m_nu=m_nu,
    )
    if not np.isclose(float(model.Ok0), omega_k, atol=2e-8, rtol=0):
        raise ValueError(f"Astropy curvature closure mismatch: requested {omega_k}, got {model.Ok0}")
    return model, omega_r, omega_de


def make_cosmology(block, config):
    cp = names.cosmological_parameters
    hubble = float(block[cp, "h0"]) * 100.0 if block.has_value(cp, "h0") else float(block[cp, "hubble"])
    return build_cosmology(
        hubble, float(block[cp, "omega_m"]), float(block[cp, "omega_k"]),
        float(block[cp, "w"]), config["tcmb"], config["neff"],
    )


def execute(block, config):
    model, omega_r, omega_de = make_cosmology(block, config)
    z = np.linspace(0.0, config["zmax"], config["nz"])
    a = 1.0 / (1.0 + z)
    d_c = model.comoving_distance(z).to_value(u.Mpc)
    omega_k = float(model.Ok0)
    if np.isclose(omega_k, 0.0, atol=1e-14):
        d_m = d_c
    else:
        root = np.sqrt(abs(omega_k))
        dh = model.hubble_distance.to_value(u.Mpc)
        if omega_k > 0:
            d_m = dh / root * np.sinh(root * d_c / dh)
        else:
            d_m = dh / root * np.sin(root * d_c / dh)

    d_l = d_m * (1.0 + z)
    d_a = d_m / (1.0 + z)
    mu = np.full_like(z, -np.inf)
    mu[1:] = 5.0 * np.log10(d_l[1:]) + 25.0
    h_z = (model.H(z) / astropy.constants.c).to_value(1 / u.Mpc)
    d_v = ((1.0 + z) ** 2 * z * d_a**2 / h_z) ** (1.0 / 3.0)

    dist = names.distances
    for key, value in {
        "z": z, "a": a, "mu": mu, "D_L": d_l, "D_A": d_a,
        "D_V": d_v, "D_M": d_m, "D_C": d_c, "H": h_z,
    }.items():
        block[dist, key] = value
    for key in ("D_L", "D_A", "D_M", "D_C", "D_V"):
        block.put_metadata(dist, key, "unit", "Mpc")
    block.put_metadata(dist, "H", "unit", "1.0/Mpc")

    cp = names.cosmological_parameters
    block[cp, "omega_r"] = omega_r
    block[cp, "omega_lambda"] = omega_de
    block[cp, "omega_k"] = omega_k
    block[cp, "closure_sum"] = float(model.Om0 + model.Ogamma0 + model.Onu0 + model.Ode0 + model.Ok0)
    block[cp, "q0"] = float(0.5 * model.Om0 + omega_r + 0.5 * (1.0 + 3.0 * model.w0) * model.Ode0)
    return 0
