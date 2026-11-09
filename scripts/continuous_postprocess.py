"""CosmoSIS postprocessor extension for a four-parameter corner comparison."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from cosmosis.plotting.kde import KDE
from cosmosis.postprocessing.elements import PostProcessorElement
from cosmosis.postprocessing.outputs import PostprocessPlot


ROOT = Path(__file__).resolve().parent.parent
SUMMARY = ROOT / ".run-data/continuous/summary.json"
PARAMETERS = (
    ("cosmological_parameters--h0", r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", (45, 110), 100),
    ("cosmological_parameters--omega_m", r"$\Omega_m$", (0.20, 0.40), 1),
    ("cosmological_parameters--omega_k", r"$\Omega_k$", (-0.01, 0.01), 1),
    ("cosmological_parameters--w", r"$w$", (-1.0, -0.333334), 1),
)
DATASETS = ("DES-Y3 BAO", "Pantheon+SH0ES", "DESI DR2 BAO", "Planck 2018 compressed")
COLORS = ("#0876ac", "#e49c21", "#23a782", "#6a4ca2")


def contour_levels(density, xs, ys, gx, gy):
    """Return 95% and 68% KDE levels using sample counts on the KDE grid.

    Ordering cells by density and accumulating their sample counts yields
    contour thresholds for the empirical credible regions.
    """
    dx, dy = gx[1] - gx[0], gy[1] - gy[0]
    counts, _, _ = np.histogram2d(xs, ys,
                                  bins=[np.r_[gx - dx / 2, gx[-1] + dx / 2],
                                        np.r_[gy - dy / 2, gy[-1] + dy / 2]])
    order = np.argsort(density.ravel())[::-1]
    cumulative = np.cumsum(counts.ravel()[order]) / max(counts.sum(), 1)
    sorted_density = density.ravel()[order]
    outer = sorted_density[min(np.searchsorted(cumulative, 0.95), len(order) - 1)]
    inner = sorted_density[min(np.searchsorted(cumulative, 0.68), len(order) - 1)]
    return float(outer), float(inner)


class ContinuousCorner(PostProcessorElement):
    """Accumulate each independent chain into the comparison plot products."""

    def run(self):
        """Draw one dataset's histograms, H0 summary, and corner contours."""
        index = self.source.index
        if index >= len(DATASETS):
            raise ValueError("Expected one chain for each dataset")
        product = self.get_output("continuous_corner")
        if product is None:
            filename = self.filename(self.options.get("file_type", "png"), "continuous_corner")
            product = PostprocessPlot("continuous_corner", filename, self.new_figure())
            self.set_output("continuous_corner", product)
        fig = product.value
        axes = fig._study_axes
        color = COLORS[index]
        vectors = []
        for name, _, _, scale in PARAMETERS:
            values = np.asarray(self.reduced_col(name)) * scale
            vectors.append(values)
        valid = np.logical_and.reduce([np.isfinite(vector) for vector in vectors])
        vectors = [vector[valid] for vector in vectors]
        if len(vectors[0]) < 1000:
            raise RuntimeError("Too few post-burn samples for the corner plot")
        summary = next(row for row in json.loads(SUMMARY.read_text())
                       if row["dataset"] == DATASETS[index])

        histogram_product = self.get_output("continuous_h0_histograms")
        if histogram_product is None:
            histogram_name = self.filename(self.options.get("file_type", "png"),
                                           "continuous_h0_histograms")
            histogram_product = PostprocessPlot("continuous_h0_histograms", histogram_name,
                                                self.new_histogram_figure())
            self.set_output("continuous_h0_histograms", histogram_product)
        h0_hist_ax = histogram_product.value.axes[index]
        h0_edges = np.linspace(vectors[0].min(), vectors[0].max(), 41)
        h0_hist_ax.set_xlim(h0_edges[0], h0_edges[-1])
        h0_hist_ax.hist(vectors[0], bins=h0_edges, density=True, color=color,
                        edgecolor="white", linewidth=0.4, alpha=0.86)
        h0_hist_ax.axvline(summary["gaussian_peak_h0"], color="#222", lw=1.4,
                           ls="--", label=f"Mean {summary['gaussian_peak_h0']:.2f}")
        h0_hist_ax.legend(frameon=False, fontsize=9, loc="upper right")
        h0_hist_ax.set_title(DATASETS[index], fontsize=11, fontweight="bold", loc="left")

        shared_product = self.get_output("continuous_h0_histograms_50_80")
        if shared_product is None:
            shared_name = self.filename(self.options.get("file_type", "png"),
                                        "continuous_h0_histograms_50_80")
            shared_product = PostprocessPlot("continuous_h0_histograms_50_80",
                                             shared_name,
                                             self.new_histogram_figure(shared=True))
            self.set_output("continuous_h0_histograms_50_80", shared_product)
        shared_ax = shared_product.value.axes[index]
        shared_edges = np.linspace(50, 80, 61)
        bin_width = shared_edges[1] - shared_edges[0]
        # Normalize against the full chain so the visible area shows its fraction in range.
        weights = np.full(len(vectors[0]), 1.0 / (len(vectors[0]) * bin_width))
        shared_ax.hist(vectors[0], bins=shared_edges, weights=weights,
                       color=color, edgecolor="white", linewidth=0.3, alpha=0.86)
        shared_ax.axvline(summary["gaussian_peak_h0"], color="#222", lw=1.4,
                          ls="--", label=f"Mean {summary['gaussian_peak_h0']:.2f}")
        fraction = np.mean((vectors[0] >= 50) & (vectors[0] <= 80))
        shared_ax.text(0.02, 0.94, f"{fraction:.0%} of chain in view",
                       transform=shared_ax.transAxes, ha="left", va="top",
                       fontsize=8.5, color="#444")
        shared_ax.legend(frameon=False, fontsize=9, loc="upper right")
        shared_ax.set_title(DATASETS[index], fontsize=11, fontweight="bold", loc="left")

        parameter_product = self.get_output("continuous_parameter_histograms")
        if parameter_product is None:
            parameter_name = self.filename(self.options.get("file_type", "png"),
                                           "continuous_parameter_histograms")
            parameter_product = PostprocessPlot("continuous_parameter_histograms",
                                                parameter_name,
                                                self.new_parameter_histogram_figure())
            self.set_output("continuous_parameter_histograms", parameter_product)
        for column, vector in enumerate(vectors):
            ax = parameter_product.value.axes[index * len(PARAMETERS) + column]
            edges = h0_edges if column == 0 else np.linspace(*PARAMETERS[column][2], 31)
            if column == 0:
                ax.set_xlim(h0_edges[0], h0_edges[-1])
            ax.hist(vector, bins=edges, density=True, color=color,
                    edgecolor="white", linewidth=0.25, alpha=0.86)

        h0_product = self.get_output("continuous_h0")
        if h0_product is None:
            h0_name = self.filename(self.options.get("file_type", "png"), "continuous_h0")
            h0_product = PostprocessPlot("continuous_h0", h0_name, self.new_h0_figure())
            self.set_output("continuous_h0", h0_product)
        h0_ax = h0_product.value.axes[0]
        h0_x = np.linspace(*PARAMETERS[0][2], 400)
        h0_mean, h0_sigma = summary["gaussian_peak_h0"], summary["gaussian_sigma_h0"]
        h0_ax.axvspan(h0_mean - h0_sigma, h0_mean + h0_sigma,
                      color=color, alpha=0.075)
        h0_ax.plot(h0_x, np.exp(-0.5 * ((h0_x - h0_mean) / h0_sigma) ** 2),
                   color=color, lw=2.3,
                   label=f"{DATASETS[index]}: {h0_mean:.2f} ± {h0_sigma:.2f}")
        h0_ax.axvline(h0_mean, color=color, lw=1.3, ls=":")
        h0_ax.text(h0_mean, 1.065 - index * 0.12, f"{h0_mean:.2f}", color=color,
                   ha="center", va="top", fontsize=10, fontweight="bold",
                   bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1})

        for column, vector in enumerate(vectors):
            ax = axes[column, column]
            x, density = KDE([vector], factor=1.5).grid_evaluate(220, PARAMETERS[column][2])
            ax.plot(x, density / density.max(), color=color, lw=1.9)
            if column == 0:
                mean, sigma = summary["gaussian_peak_h0"], summary["gaussian_sigma_h0"]
                gaussian = np.exp(-0.5 * ((x - mean) / sigma) ** 2)
                ax.plot(x, gaussian, ls="--", color=color, lw=1.5)
                ax.axvline(mean, ls=":", color=color, lw=1.1)
                ax.text(mean, 1.08 - index * 0.15, f"{mean:.1f} ± {sigma:.1f}",
                        ha="center", va="top", color=color, fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1})

        for row in range(1, len(PARAMETERS)):
            for column in range(row):
                ax = axes[row, column]
                xs, ys = vectors[column], vectors[row]
                (gx, gy), density = KDE([xs, ys], factor=1.5).grid_evaluate(
                    84, [PARAMETERS[column][2], PARAMETERS[row][2]])
                outer, inner = contour_levels(density, xs, ys, gx, gy)
                upper = float(density.max()) * 1.001
                # Nearly flat or coincident levels cannot define nested contours.
                if 0 < outer < inner < upper:
                    ax.contourf(gx, gy, density.T, levels=[outer, upper], colors=[color], alpha=0.14)
                    ax.contourf(gx, gy, density.T, levels=[inner, upper], colors=[color], alpha=0.26)
                    ax.contour(gx, gy, density.T, levels=[outer, inner], colors=[color],
                               linewidths=[0.65, 1.15], alpha=0.85)
        if index == len(DATASETS) - 1:
            h0_ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
        return [product.filename, h0_product.filename,
                histogram_product.filename, shared_product.filename,
                parameter_product.filename]

    @staticmethod
    def new_histogram_figure(shared=False):
        """Prepare four H0 histograms, optionally with a common 50-80 axis."""
        fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.8), facecolor="white")
        fig.subplots_adjust(left=0.09, right=0.97, bottom=0.13, top=0.80,
                            hspace=0.34, wspace=0.23)
        for ax in axes.flat:
            if shared:
                ax.set_xlim(50, 80)
            ax.set_xlabel(PARAMETERS[0][1])
            ax.set_ylabel("Posterior density")
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(axis="y", alpha=0.14)
        title = ("CosmoSIS chain histograms · common $H_0$ axis"
                 if shared else "CosmoSIS chain histograms · marginalized $H_0$")
        fig.text(0.09, 0.96, title,
                 fontsize=16, fontweight="bold", va="top")
        fig.text(0.09, 0.91,
                 r"Sampled: $\Omega_m\in[0.20,0.40]$, $\Omega_k\in[-0.01,0.01]$, "
                 r"$w\in[-1,-1/3)$; $\Omega_{DE}$ by closure",
                 fontsize=10, va="top")
        note = ("Bars: 0.5 km/s/Mpc bins on a shared 50–80 axis · "
                "dashed: full-chain mean · out-of-range samples omitted"
                if shared else
                "Bars: post-burn chain frequencies (40 bins over each observed range) · "
                "dashed: chain mean · note different H0 axes")
        fig.text(0.09, 0.855, note,
                 fontsize=9, va="top", color="#444")
        footer = ("Bars retain each chain's full-sample normalization; displayed area is the "
                  "fraction inside 50–80. All chains remain provisional."
                  if shared else
                  "Density normalization accounts for different chain lengths. "
                  "All four chains remain provisional: convergence criterion not met.")
        fig.text(0.09, 0.04, footer,
                 fontsize=9, color="#8b3030")
        return fig

    @staticmethod
    def new_parameter_histogram_figure():
        """Prepare dataset rows and shared-parameter columns for chain histograms."""
        fig, axes = plt.subplots(4, 4, figsize=(13.2, 11.4), facecolor="white")
        fig.subplots_adjust(left=0.16, right=0.97, bottom=0.10, top=0.84,
                            hspace=0.38, wspace=0.29)
        for row, label in enumerate(DATASETS):
            for column, (_, axis_label, limits, _) in enumerate(PARAMETERS):
                ax = axes[row, column]
                ax.set_xlim(*limits)
                ax.tick_params(direction="in", top=True, right=True, labelsize=8)
                ax.grid(axis="y", alpha=0.12)
                if row == 3:
                    ax.set_xlabel(axis_label, fontsize=10)
                if column == 0:
                    ax.set_ylabel("Posterior density", fontsize=9)
            fig.text(0.015, 0.745 - row * 0.198, label, rotation=90,
                     va="center", ha="left", fontsize=10, fontweight="bold",
                     color=COLORS[row])
        fig.text(0.16, 0.97, "CosmoSIS chain histograms · shared cosmological parameters",
                 fontsize=16, fontweight="bold", va="top")
        fig.text(0.16, 0.915,
                 r"Columns: $H_0$, $\Omega_m$, $\Omega_k$, $w$ · "
                 r"$\Omega_{DE}$ derived by closure · 30 bins per density parameter",
                 fontsize=10, va="top")
        fig.text(0.16, 0.055,
                 "Post-burn, thinned chains · density normalized within each panel · "
                 "all results provisional pending convergence",
                 fontsize=9, color="#8b3030")
        return fig

    @staticmethod
    def new_h0_figure():
        """Prepare the axis for overlaid Gaussian H0 summaries from each chain."""
        fig, ax = plt.subplots(figsize=(14.0, 6.4), facecolor="white")
        fig.subplots_adjust(left=0.08, right=0.60, bottom=0.24, top=0.78)
        ax.set_xlim(*PARAMETERS[0][2])
        ax.set_ylim(0, 1.14)
        ax.set_xlabel(PARAMETERS[0][1], fontsize=11)
        ax.set_ylabel("Relative Gaussian density (each peak = 1)", fontsize=11)
        ax.grid(alpha=0.15)
        ax.tick_params(direction="in", top=True, right=True)
        fig.text(0.09, 0.96, "Marginalized $H_0$ · continuous curved $w$CDM",
                 ha="left", va="top", fontsize=16, fontweight="bold")
        fig.text(0.09, 0.895,
                 r"Sampled: $\Omega_m\in[0.20,0.40]$, $\Omega_k\in[-0.01,0.01]$, "
                 r"$w\in[-1,-1/3)$; $\Omega_{DE}$ by closure",
                 ha="left", va="top", fontsize=10)
        fig.text(0.09, 0.065,
                 "Solid: Gaussian summary  ·  Dotted: Gaussian peak (mean)  ·  "
                 "Pale band: Gaussian mean ± 1σ",
                 ha="left", va="bottom", fontsize=9, color="#444")
        fig.text(0.09, 0.025,
                 "Provisional chain diagnostics: none of the four chains passed the agreed autocorrelation rule.",
                 ha="left", va="bottom", fontsize=9, color="#8b3030")
        return fig

    @staticmethod
    def new_figure():
        """Prepare the four-parameter corner grid and cache its axes for reuse."""
        fig, axes = plt.subplots(4, 4, figsize=(13.2, 12.0), facecolor="white")
        fig.subplots_adjust(left=0.10, right=0.94, bottom=0.10, top=0.81,
                            hspace=0.09, wspace=0.10)
        for row in range(4):
            for column in range(4):
                ax = axes[row, column]
                if row < column:
                    ax.set_visible(False)
                    continue
                ax.set_xlim(*PARAMETERS[column][2])
                if row == column:
                    ax.set_ylim(0, 1.14)
                    ax.set_yticks([0, 0.5, 1.0] if column == 0 else [])
                else:
                    ax.set_ylim(*PARAMETERS[row][2])
                ax.tick_params(direction="in", top=True, right=True, labelsize=8)
                ax.grid(alpha=0.10)
                if row < 3:
                    ax.tick_params(labelbottom=False)
                else:
                    ax.set_xlabel(PARAMETERS[column][1], fontsize=11)
                    ax.set_xticks(((50, 70, 90), (0.22, 0.30, 0.38),
                                   (-0.008, 0.0, 0.008), (-0.9, -0.7, -0.5))[column])
                if column > 0:
                    ax.tick_params(labelleft=False)
                elif row > 0:
                    ax.set_ylabel(PARAMETERS[row][1], fontsize=11)
        axes[0, 0].set_ylabel("Relative density", fontsize=9)
        fig._study_axes = axes
        fig.text(0.10, 0.965, "Continuous curved $w$CDM · separate dataset posteriors",
                 fontsize=17, fontweight="bold", ha="left", va="top")
        fig.text(0.10, 0.921,
                 r"Sampled: $H_0\in[45,110]$, $\Omega_m\in[0.20,0.40]$, "
                 r"$\Omega_k\in[-0.01,0.01]$, $w\in[-1,-1/3)$  ·  $\Omega_{DE}$ by closure",
                 fontsize=10, ha="left", va="top")
        fig.text(0.10, 0.884,
                 "BAO sound horizon calibrated by CAMB  ·  Pantheon magnitude and Planck baryon density sampled",
                 fontsize=9.5, ha="left", va="top")
        handles = [Patch(facecolor=color, edgecolor=color, alpha=0.65, label=label)
                   for label, color in zip(DATASETS, COLORS)]
        handles += [Line2D([0], [0], color="#555", lw=1.8, label="Sampled posterior"),
                    Line2D([0], [0], color="#555", lw=1.8, ls="--", label="$H_0$ Gaussian summary")]
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.10, 0.856),
                   ncol=4, frameon=False, fontsize=9)
        fig.text(0.10, 0.045,
                 "Diagonal: marginalized posteriors  ·  Lower triangle: 68%/95% joint contours  ·  "
                 "Dotted $H_0$ lines: Gaussian peaks (mean ± 1σ)",
                 fontsize=9, color="#444", ha="left")
        fig.text(0.10, 0.021,
                 "Provisional chain diagnostics: convergence criterion not met by any dataset.",
                 fontsize=9, color="#8b3030", ha="left")
        return fig
