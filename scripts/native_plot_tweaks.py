"""Presentation labels for plots made by CosmoSIS's postprocessor.

The solid posterior curves are drawn by CosmoSIS. This supported tweak hook
converts its h axis to H0, adds model labels, and draws the tabulated Gaussian
summaries as explicitly identified approximations.
"""

import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from cosmosis.postprocessing.plots import Tweaks


SUMMARY = Path(__file__).resolve().parent.parent / ".run-data/curved/curved_summary.csv"
GROUPS = {
    "omega_k_minus_0p01": ("Ωₖ = −0.01", -0.01, ("Pantheon+SH0ES", "DES-Y3 BAO")),
    "omega_k_zero": ("Ωₖ = 0.00", 0.0, ("Pantheon+SH0ES", "DES-Y3 BAO")),
    "omega_k_plus_0p01": ("Ωₖ = +0.01", 0.01, ("Pantheon+SH0ES", "DES-Y3 BAO")),
    "pantheon_curvature": ("Pantheon+SH0ES", None, ("Ωₖ = −0.01", "Ωₖ = 0.00", "Ωₖ = +0.01")),
    "bao_curvature": ("DES-Y3 BAO", None, ("Ωₖ = −0.01", "Ωₖ = 0.00", "Ωₖ = +0.01")),
}
COLORS = ("#bb3f4b", "#277da1", "#8661aa")


class CurvedPosteriorPresentation(Tweaks):
    filename = ["cosmological_parameters--h0", "cosmological_parameters--w"]

    def run(self):
        group = os.environ["COSMOSIS_FIGURE_GROUP"]
        title, omega_k, labels = GROUPS[group]
        parameter = self.info[0]
        is_h0 = parameter.endswith("--h0")
        fig = plt.gcf()
        ax = fig.axes[0]
        native_lines = list(ax.lines)
        if len(native_lines) != len(labels):
            raise RuntimeError(f"Expected {len(labels)} native curves in {group}, found {len(native_lines)}")

        for line, label, color in zip(native_lines, labels, COLORS):
            if is_h0:
                line.set_xdata(np.asarray(line.get_xdata()) * 100.0)
            line.set_color(color)
            line.set_linewidth(2.4)
            line.set_label(f"{label} · sampled posterior")

        if is_h0 and omega_k is not None:
            with SUMMARY.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            x = np.linspace(60.0, 95.0, 700)
            for label, color in zip(labels, COLORS):
                row = next(row for row in rows
                           if row["dataset"] == label and float(row["omega_k"]) == omega_k)
                mean = float(row["gaussian_peak_H0"])
                sigma = float(row["gaussian_sigma_H0"])
                gaussian = np.exp(-0.5 * ((x - mean) / sigma) ** 2)
                ax.plot(x, gaussian, "--", color=color, linewidth=1.8,
                        label=f"{label} · Gaussian {mean:.2f} ± {sigma:.2f}")
                ax.axvline(mean, color=color, linestyle=":", linewidth=1.2)
                ax.text(mean, 1.035, f"{mean:.2f}", ha="center", va="bottom",
                        fontsize=9, color=color, fontweight="bold")

        fig.set_size_inches(11.4, 6.3)
        ax.set_position([0.09, 0.19, 0.58, 0.63])
        ax.set_ylim(0, 1.16)
        ax.set_ylabel("Relative posterior density (each peak = 1)", fontsize=11)
        if is_h0:
            ax.set_xlim(60, 95)
            ax.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", fontsize=11)
            subject = "Hubble constant"
        else:
            ax.set_xlim(-1.5, -1.3)
            ax.set_xlabel(r"$w$ (dimensionless)", fontsize=11)
            subject = "Dark-energy equation of state"
        ax.grid(alpha=0.17)
        ax.tick_params(labelsize=10)
        old_legend = ax.get_legend()
        if old_legend is not None:
            old_legend.remove()
        ax.legend(loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False,
                  fontsize=9, borderaxespad=0)

        fig.text(0.09, 0.955, f"{subject} posterior · {title}",
                 ha="left", va="top", fontsize=16, fontweight="bold")
        fig.text(0.09, 0.905,
                 "Ωₘ = 0.30 fixed   ·   w ∈ [−1.5, −1.3] (phantom)   ·   standard radiation   ·   ΩDE from closure",
                 ha="left", va="top", fontsize=9.5)
        if is_h0:
            footer = ("Solid: CosmoSIS sampled posterior   ·   Dashed: Gaussian summary   ·   "
                      "Dotted: Gaussian mean") if omega_k is not None else "Solid: CosmoSIS sampled posterior"
            footer += "   ·   BAO upper H₀ prior = 95"
        else:
            footer = "Solid: CosmoSIS sampled posterior   ·   w prior = [−1.5, −1.3]"
        fig.text(0.09, 0.065, footer, ha="left", va="bottom", fontsize=8.5, color="#444444")
