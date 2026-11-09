"""CosmoSIS postprocessing extension for the curved H0–w corner comparison.

Omega_k is fixed separately in each chain, so it is represented by color rather
than by a continuous corner-plot axis. Only H0 and w were sampled continuously.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from cosmosis.plotting.kde import KDE
from cosmosis.postprocessing.elements import PostProcessorElement
from cosmosis.postprocessing.outputs import PostprocessPlot


CURVATURE = (-0.01, 0.0, 0.01)
COLORS = ("#0876ac", "#e49c21", "#24a989")
DATASETS = ("Pantheon+SH0ES", "DES-Y3 BAO")
H0_RANGE = (60.0, 95.0)
W_RANGE = (-1.5, -1.3)


def credible_levels(density, h0, w, h0_grid, w_grid):
    """Return outer/inner KDE thresholds enclosing 95%/68% of the samples.

    Histogram counts on the KDE grid supply the enclosed probability; the
    density only determines which cells are visited first.
    """
    h0_step = h0_grid[1] - h0_grid[0]
    w_step = w_grid[1] - w_grid[0]
    h0_edges = np.r_[h0_grid - h0_step / 2, h0_grid[-1] + h0_step / 2]
    w_edges = np.r_[w_grid - w_step / 2, w_grid[-1] + w_step / 2]
    counts, _, _ = np.histogram2d(h0, w, bins=[h0_edges, w_edges])
    order = np.argsort(density.ravel())[::-1]
    cumulative = np.cumsum(counts.ravel()[order]) / counts.sum()
    flat = density.ravel()[order]
    inner = flat[np.searchsorted(cumulative, 0.68)]
    outer = flat[np.searchsorted(cumulative, 0.95)]
    return float(outer), float(inner)


class CurvedCorner(PostProcessorElement):
    """One CosmoSIS plot product, filled by six independent chain inputs."""

    def run(self):
        """Add one chain's marginals and joint contours to its dataset panel."""
        index = self.source.index
        if index >= 6:
            raise ValueError("Expected three fixed-curvature chains per dataset")
        dataset_index, curvature_index = divmod(index, 3)
        color = COLORS[curvature_index]
        h0 = np.asarray(self.reduced_col("cosmological_parameters--h0")) * 100.0
        w = np.asarray(self.reduced_col("cosmological_parameters--w"))
        valid = np.isfinite(h0) & np.isfinite(w)
        h0, w = h0[valid], w[valid]
        if len(h0) < 1000:
            raise RuntimeError("Too few post-burn samples for the corner plot")

        product = self.get_output("h0_w_corner")
        if product is None:
            fig = self._new_figure()
            filename = self.filename(self.options.get("file_type", "png"), "h0_w_corner")
            product = PostprocessPlot("h0_w_corner", filename, fig)
            self.set_output("h0_w_corner", product)
        fig = product.value
        top, joint, w_diagonal = fig._study_axes[dataset_index]

        x_h0, h0_density = KDE([h0], factor=1.6).grid_evaluate(240, H0_RANGE)
        x_w, w_density = KDE([w], factor=1.6).grid_evaluate(240, W_RANGE)
        top.plot(x_h0, h0_density / h0_density.max(), color=color, lw=1.8)
        w_diagonal.plot(x_w, w_density / w_density.max(), color=color, lw=1.8)

        (h0_grid, w_grid), density = KDE([h0, w], factor=1.6).grid_evaluate(
            80, [H0_RANGE, W_RANGE])
        outer, inner = credible_levels(density, h0, w, h0_grid, w_grid)
        z = density.T
        joint.contourf(h0_grid, w_grid, z, levels=[outer, z.max() * 1.001],
                       colors=[color], alpha=0.16)
        joint.contourf(h0_grid, w_grid, z, levels=[inner, z.max() * 1.001],
                       colors=[color], alpha=0.31)
        joint.contour(h0_grid, w_grid, z, levels=[outer, inner],
                      colors=[color], linewidths=[0.7, 1.25], alpha=0.85)
        return [product.filename]

    @staticmethod
    def _new_figure():
        """Create paired H0-w corner panels and cache their axes for later chains."""
        fig = plt.figure(figsize=(14.2, 6.7), facecolor="white")
        layout = fig.add_gridspec(2, 5, left=0.065, right=0.965, bottom=0.13,
                                 top=0.78, hspace=0.04, wspace=0.07,
                                 width_ratios=[1, 1, 0.34, 1, 1])
        groups = []
        for dataset_index, start in enumerate((0, 3)):
            top = fig.add_subplot(layout[0, start])
            joint = fig.add_subplot(layout[1, start], sharex=top)
            w_diagonal = fig.add_subplot(layout[1, start + 1])
            groups.append((top, joint, w_diagonal))
            fig.text((0.29 if dataset_index == 0 else 0.77), 0.815,
                     DATASETS[dataset_index], ha="center", va="bottom",
                     fontsize=13, fontweight="bold")
            top.set_ylabel("Relative posterior density")
            joint.set_ylabel(r"$w$ (dimensionless)")
            joint.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]")
            w_diagonal.set_xlabel(r"$w$ (dimensionless)")
            top.tick_params(axis="x", labelbottom=False)
            w_diagonal.tick_params(axis="y", labelleft=False)
            joint.set_xticks([60, 70, 80, 90])
            joint.set_yticks([-1.5, -1.4, -1.3])
            w_diagonal.set_xticks([-1.5, -1.4, -1.3])
            for axis in (top, w_diagonal):
                axis.set_ylim(0, 1.07)
                axis.set_yticks([0, 0.5, 1.0])
            for axis in (top, joint):
                axis.set_xlim(*H0_RANGE)
            # The Lambda-CDM reference is outside these historical phantom-w
            # priors; show it explicitly without changing the saved samples.
            joint.set_ylim(W_RANGE[0], -0.98)
            w_diagonal.set_xlim(W_RANGE[0], -0.98)
            top.axvline(67.4, color="#222", ls="--", lw=1.1, zorder=10)
            joint.axvline(67.4, color="#222", ls="--", lw=1.1, zorder=10)
            joint.axhline(-1.0, color="#222", ls="--", lw=1.1, zorder=10)
            w_diagonal.axvline(-1.0, color="#222", ls="--", lw=1.1, zorder=10)
            for axis in (top, joint, w_diagonal):
                axis.tick_params(direction="in", top=True, right=True, labelsize=9)
                axis.grid(alpha=0.10)
        fig._study_axes = groups
        fig.text(0.065, 0.965, "Curved phantom wCDM · marginalized $H_0$–$w$ posteriors",
                 ha="left", va="top", fontsize=16, fontweight="bold")
        fig.text(0.065, 0.91,
                 "Ωₘ = 0.30 fixed  ·  Ωₖ fixed per color  ·  w ∈ [−1.5, −1.3] sampled  ·  standard radiation  ·  ΩDE by closure",
                 ha="left", va="top", fontsize=10)
        handles = [Patch(facecolor=color, edgecolor=color, alpha=0.65,
                         label=f"Ωₖ = {value:+.2f}")
                   for value, color in zip(CURVATURE, COLORS)]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.89),
                   ncol=3, frameon=False, fontsize=10)
        fig.text(0.065, 0.055,
                 "Diagonal: marginalized posterior  ·  Inner/outer shaded contours: 68%/95% credible regions  ·  BAO $H_0$ prior ≤ 95",
                 ha="left", va="bottom", fontsize=9, color="#444444")
        fig.text(0.065, 0.02,
                 "Dashed references: H₀ = 67.4 km/s/Mpc; w = −1 (outside the sampled phantom-w prior).",
                 ha="left", va="bottom", fontsize=9, color="#222222")
        return fig
