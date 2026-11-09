"""Render GetDist corner plots from the saved continuous CosmoSIS chains.

The four likelihoods were run independently. The comparison therefore overlays
only their shared parameters; individual figures include every varied parameter.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
from getdist import MCSamples, plots


ROOT = Path(__file__).resolve().parent.parent
CHAINS = ROOT / ".run-data" / "continuous"
FIGURES = ROOT / "output" / "continuous" / "figures"
BURN_FRACTION = 0.30
THIN = 10
COMMON = ("h0", "omega_m", "omega_k", "w")
# Planck 2018 base-Lambda-CDM benchmarks; M is the configured nuisance
# fiducial, not a universal constant. Values use the plotted units.
REFERENCE_VALUES = {"h0": 67.4, "omega_m": 0.315, "omega_k": 0.0,
                    "w": -1.0, "omega_b": 0.0493, "m": -19.35}
PARAMETERS = {
    "cosmological_parameters--h0": ("h0", r"H_0\,[\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{Mpc}^{-1}]", 100),
    "cosmological_parameters--omega_m": ("omega_m", r"\Omega_m", 1),
    "cosmological_parameters--omega_k": ("omega_k", r"\Omega_k", 1),
    "cosmological_parameters--w": ("w", "w", 1),
    "supernova_params--m": ("m", "M", 1),
    "cosmological_parameters--omega_b": ("omega_b", r"\Omega_b", 1),
}
RUNS = (
    ("des_y3_bao", "DES-Y3 BAO", "#0876ac"),
    ("pantheonplussh0es", "Pantheon+SH0ES", "#d17b13"),
    ("desi_dr2_bao", "DESI DR2 BAO", "#198a68"),
    ("planck_compressed", "Planck 2018 compressed", "#6a4ca2"),
)


def load_chain(stem, label):
    """Load one CosmoSIS chain as labeled GetDist samples for corner plots.

    Drop the prior/post columns, discard the initial burn fraction, thin the
    remaining rows, and convert the stored h0 fraction to H0 in km/s/Mpc.
    """
    path = CHAINS / f"{stem}.txt"
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().removeprefix("#").strip().split("\t")
    if len(header) < 3 or header[-2:] != ["prior", "post"]:
        raise ValueError(f"Unexpected CosmoSIS header in {path}")
    unknown = set(header[:-2]) - PARAMETERS.keys()
    if unknown:
        raise ValueError(f"Unrecognised varied parameters in {path}: {sorted(unknown)}")
    values = np.loadtxt(path, comments="#", delimiter="\t", ndmin=2)
    if values.shape[1] != len(header) or len(values) < 1000:
        raise ValueError(f"Incomplete or malformed chain in {path}: {values.shape}")
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite chain entries in {path}")
    values = values[int(len(values) * BURN_FRACTION)::THIN, :len(header) - 2]
    names, labels, scales = zip(*(PARAMETERS[name] for name in header[:-2]))
    samples = values * np.asarray(scales)
    return MCSamples(samples=samples, names=list(names), labels=list(labels),
                     label=label, settings={"ignore_rows": 0})


def render(samples, params, colors, labels, destination, width):
    """Export a filled triangle plot for the chosen samples and parameters."""
    plotter = plots.get_subplot_plotter(width_inch=width)
    plotter.settings.num_plot_contours = 2
    plotter.settings.axes_fontsize = 10
    plotter.settings.lab_fontsize = 12
    plotter.settings.legend_fontsize = 10
    plotter.triangle_plot(samples, params=list(params), filled=True,
                         contour_colors=list(colors), legend_labels=list(labels),
                         legend_loc="upper right")
    for row, y_name in enumerate(params):
        for column, x_name in enumerate(params[:row + 1]):
            ax = plotter.subplots[row, column]
            # Include references outside the posterior's automatic axis limits.
            x_limits = ax.get_xlim()
            x_ref = REFERENCE_VALUES[x_name]
            ax.set_xlim(min(x_limits[0], x_ref), max(x_limits[1], x_ref))
            ax.axvline(x_ref, color="#222222", ls="--", lw=1.1, zorder=10)
            if row != column:
                y_limits = ax.get_ylim()
                y_ref = REFERENCE_VALUES[y_name]
                ax.set_ylim(min(y_limits[0], y_ref), max(y_limits[1], y_ref))
                ax.axhline(y_ref, color="#222222", ls="--", lw=1.1, zorder=10)
    references = ", ".join(
        f"${PARAMETERS[key][1]} = {REFERENCE_VALUES[name]:g}$"
        for name in params for key in PARAMETERS if PARAMETERS[key][0] == name)
    plotter.fig.text(0.5, 0.005, "Dashed reference lines: " + references,
                     ha="center", va="bottom", fontsize=8)
    plotter.export(str(destination))


def main():
    """Write the shared-parameter overlay and each chain's full corner plot."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    chains = [load_chain(stem, label) for stem, label, _ in RUNS]
    render(chains, COMMON, [color for _, _, color in RUNS],
           [label for _, label, _ in RUNS], FIGURES / "continuous_corner.png", 12)
    for chain, (stem, label, color) in zip(chains, RUNS):
        render([chain], chain.getParamNames().list(), [color], [label],
               FIGURES / f"continuous_corner_{stem}.png", 10)


if __name__ == "__main__":
    main()
