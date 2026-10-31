#!/usr/bin/env python3
"""Run and summarize the three curved phantom-wCDM comparisons."""

from __future__ import annotations

import configparser
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import tempfile
from pathlib import Path

import numpy as np

import hubble_analysis as base


HERE = Path(__file__).resolve().parent.parent
CONFIGS = HERE / "configs"
SCRIPTS = HERE / "scripts"
OUT = HERE / ".run-data" / "curved"
K_VALUES = (-0.01, 0.0, 0.01)
W_PRIOR = (-1.5, -1.3)
BASELINES = {
    "Pantheon+SH0ES": ("pantheon_wcdm.ini", "pantheon_wcdm_values.ini",
                       ("cosmological_parameters--h0", "cosmological_parameters--w",
                        "supernova_params--m")),
    "DES-Y3 BAO": ("des_y3_bao_wcdm.ini", "des_y3_bao_wcdm_values.ini",
                   ("cosmological_parameters--h0", "cosmological_parameters--w")),
}


def read_config(path: Path):
    config = configparser.ConfigParser(interpolation=None)
    config.optionxform = str
    config.read(path)
    return config


def make_run_files(label: str, omega_k: float, temp: Path):
    ini_name, values_name, _ = BASELINES[label]
    cfg = read_config(CONFIGS / ini_name)
    values = read_config(CONFIGS / values_name)
    cp = values["cosmological_parameters"]
    cp["w"] = f"{W_PRIOR[0]} -1.4 {W_PRIOR[1]}"
    cp["omega_k"] = f"{omega_k:.2f}"
    cp["omega_m"] = "0.3"

    values_path = temp / f"{base.slug(label)}_values.ini"
    values_path.write_text(base._config_text(values))
    cfg["pipeline"]["values"] = str(values_path)
    chain_path = OUT / f"{base.slug(label)}_ok_{omega_k:+.2f}.txt"
    cfg["output"]["filename"] = str(chain_path)

    if label == "Pantheon+SH0ES":
        cfg["pipeline"]["modules"] = "consistency curved_astropy_background pantheon_plus"
        cfg["curved_astropy_background"] = {
            "file": str(SCRIPTS / "curved_astropy_background.py"),
            "zmax": "3.0", "nz": "500", "tcmb": "2.7255", "neff": "3.046",
        }
    cfg_path = temp / f"{base.slug(label)}_ok_{omega_k:+.2f}.ini"
    cfg_path.write_text(base._config_text(cfg))
    return cfg_path, chain_path


def validation_grid():
    """Check closure, positive H(z)^2, and acceleration across prior corners."""
    from curved_astropy_background import build_cosmology

    rows = []
    z = np.linspace(0.0, 3.0, 301)
    for omega_k in K_VALUES:
        for w in W_PRIOR:
            model, omega_r, omega_de = build_cosmology(
                70.0, 0.3, omega_k, w, 2.7255, 3.046)
            closure = float(model.Om0 + model.Ogamma0 + model.Onu0 + model.Ode0 + model.Ok0)
            e2 = np.asarray(model.efunc(z)) ** 2
            q0 = float(0.5 * model.Om0 + omega_r + 0.5 * (1.0 + 3.0 * w) * model.Ode0)
            rows.append({"omega_k": omega_k, "w": w, "omega_m": 0.3,
                         "omega_r": omega_r, "omega_de": omega_de,
                         "closure_sum": closure, "min_E2_z_0_to_3": float(e2.min()),
                         "q0": q0, "accelerating_today": q0 < 0})
            if not np.isclose(closure, 1.0, atol=2e-8) or e2.min() <= 0 or q0 >= 0:
                raise RuntimeError(f"Physical-prior validation failed: {rows[-1]}")
    return rows


def posterior_samples(header, chain):
    burn = int(0.3 * chain.shape[0])
    h0 = chain[burn:, :, header.index("cosmological_parameters--h0")].reshape(-1) * 100.0
    w = chain[burn:, :, header.index("cosmological_parameters--w")].reshape(-1)
    valid = np.isfinite(h0) & np.isfinite(w)
    if np.count_nonzero(valid) < 100:
        raise RuntimeError("Too few finite post-burn samples")
    return h0[valid], w[valid]


def summarize(samples):
    summary = base.density_summary(samples, plot_points=500)
    return {
        "gaussian_peak_H0": summary["mean_and_gaussian_peak"],
        "gaussian_sigma_H0": summary["gaussian_sigma"],
        "median_H0": summary["median"],
        "q16_H0": summary["q16"],
        "q84_H0": summary["q84"],
        "post_burn_samples": summary["n_post_burn"],
        "kde_gaussian_l1": summary["kde_gaussian_l1_on_99_6pct_interval"],
        "skewness_H0": summary["skewness"],
    }


def write_report(results, convergence):
    lines = [
        "# Curved phantom wCDM: DES-Y3 BAO and Pantheon+SH0ES", "",
        "This study fits DES-Y3 BAO and Pantheon+SH0ES independently at fixed "
        "curvature values Ωk = −0.01, 0, +0.01. Each run samples H0 and constant "
        "w ∈ [−1.5, −1.3], fixes Ωm = 0.30, and derives ΩDE from density closure. "
        "This full w range is phantom and is presented as phenomenological, not as "
        "a canonical scalar-field model.", "",
        "Pantheon distances use photons at Tcmb = 2.7255 K and Neff = 3.046 "
        "massless neutrinos, with Ωm + Ωr + ΩDE + Ωk = 1. DES-Y3 uses CAMB with "
        "the study's fixed early-universe settings to calculate rd; its H0 posterior "
        "therefore remains conditional on that sound-horizon calibration.", "",
        "The tabulated Gaussian peak is the mean of marginalized H0 samples; its "
        "1σ is their standard deviation. It need not coincide with the mode of "
        "the sampled posterior. In the CosmoSIS figures, solid curves show smoothed "
        "sampled posteriors and dashed curves show Gaussian summaries; the horizontal "
        "H0 axis uses km/s/Mpc. Because w is sampled "
        "continuously over [−1.5, −1.3], each table reports its posterior mean "
        "rather than a fixed w value.", "",
        "## Results", "",
    ]
    for omega_k in K_VALUES:
        lines += [f"### Ωk = {omega_k:+.2f}", "",
                  "| Dataset | Posterior mean w | Gaussian peak H0 (km/s/Mpc) | Gaussian 1σ (km/s/Mpc) |",
                  "|---|---:|---:|---:|"]
        for label in BASELINES:
            key = f"{label}|{omega_k:+.2f}"
            s = results[key]["summary"]
            lines.append(f"| {label} | {s['w_posterior_mean']:.3f} | "
                         f"{s['gaussian_peak_H0']:.2f} | {s['gaussian_sigma_H0']:.2f} |")
        lines.append("")
    lines += ["## Convergence", ""]
    for record in convergence:
        final = record["history"][-1]
        lines.append(f"- {record['dataset']}: converged at {final['steps_per_walker']:,} steps/walker; "
                     f"steps/τ = {', '.join(f'{x:.1f}' for x in final['steps_per_tau'])}; "
                     f"ESS = {', '.join(f'{x:.0f}' for x in final['effective_sample_size'])}.")
    lines += ["", "## Physics checks", "",
              "A prior-corner grid checked the present-day density closure, positive "
              "E(z)² from z=0 to 3, and acceleration today (q0 < 0). The validation "
              "values are in `physics_validation.json`.", "",
              "## Figures and data", "",
              "- `native_postprocess/corner/h0_w_corner.png`: H0–w corner comparison for both datasets; color distinguishes the three fixed Ωk runs. The contours use samples thinned by ten after burn-in.",
              "- `native_postprocess/`: CosmoSIS-native H0 and w posterior PNG curves, with separate comparisons for each Ωk.",
              "- `curved_summary.json` and `curved_summary.csv`: machine-readable posterior summaries.",
              "- The six CosmoSIS text chains and convergence JSON files are retained alongside the figures.", "",
              "The DES-Y3 BAO sampled H0 posterior rises toward its upper prior bound, "
              "H0 = 95 km/s/Mpc. "
              "Its Gaussian mean and standard deviation summarize a markedly non-Gaussian, "
              "prior-limited distribution; they should not be interpreted as the actual "
              "posterior mode and width.", ""]
    (OUT / "curved_wcdm_report.md").write_text("\n".join(lines))


def main():
    if not base.COSMOSIS.exists():
        raise FileNotFoundError(f"Existing CosmoSIS executable not found: {base.COSMOSIS}")
    OUT.mkdir(parents=True, exist_ok=True)
    validation = validation_grid()
    (OUT / "physics_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    jobs = [(label, omega_k, sampled)
            for label, (_, _, sampled) in BASELINES.items()
            for omega_k in K_VALUES]

    def run_one(job):
        label, omega_k, sampled = job
        with tempfile.TemporaryDirectory(prefix="curved_wcdm_") as temp_name:
            cfg_path, chain_path = make_run_files(label, omega_k, Path(temp_name))
            run_name = f"{label} Ωk={omega_k:+.2f}"
            record = base.run_fit(run_name, "", chain_path.name, sampled,
                                  output_subdir=".run-data/curved", ini_path=cfg_path,
                                  output_path=chain_path)
        header, chain = base.read_chain(chain_path)
        for param, (lower, upper) in {
            "cosmological_parameters--h0": (0.50, 0.95),
            "cosmological_parameters--w": W_PRIOR,
        }.items():
            values = chain[:, :, header.index(param)]
            if np.nanmin(values) < lower - 1e-10 or np.nanmax(values) > upper + 1e-10:
                raise RuntimeError(f"Samples escaped configured prior {param} in {chain_path}")
        k_column = "cosmological_parameters--omega_k"
        if k_column in header and not np.allclose(chain[:, :, header.index(k_column)], omega_k):
            raise RuntimeError(f"Sampled curvature differs from configured value in {chain_path}")
        h0_samples, w_samples = posterior_samples(header, chain)
        summary = summarize(h0_samples)
        summary["w_posterior_mean"] = float(np.mean(w_samples))
        return (f"{label}|{omega_k:+.2f}",
                {"dataset": label, "omega_k": omega_k,
                 "summary": summary,
                 "chain": str(chain_path)}, record)

    results, convergence = {}, []
    # The likelihood evaluations are independent. Three concurrent CosmoSIS
    # processes use the available CPU budget without oversubscribing the host.
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(run_one, job) for job in jobs]
        for future in as_completed(futures):
            key, result, record = future.result()
            results[key] = result
            convergence.append(record)
    convergence.sort(key=lambda record: record["dataset"])

    write_report(results, convergence)
    (OUT / "convergence_summary.json").write_text(json.dumps(convergence, indent=2) + "\n")
    serializable = {}
    csv_rows = []
    for omega_k in K_VALUES:
        for label in BASELINES:
            key = f"{label}|{omega_k:+.2f}"
            summary = results[key]["summary"]
            row = {"dataset": label, "omega_k": omega_k,
                   **{name: value for name, value in summary.items()
                      if not name.startswith("plot_")}}
            serializable[key] = row
            csv_rows.append(row)
    (OUT / "curved_summary.json").write_text(json.dumps(serializable, indent=2) + "\n")
    with (OUT / "curved_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Curved analysis complete. Report: {OUT / 'curved_wcdm_report.md'}")


if __name__ == "__main__":
    main()
