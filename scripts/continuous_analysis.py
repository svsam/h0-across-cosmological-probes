#!/usr/bin/env python3
"""Two continuous curved-wCDM fits, with at most three convergence checks each."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
import subprocess
from pathlib import Path
from itertools import product

import numpy as np

import hubble_analysis as base
from curved_astropy_background import build_cosmology


HERE = Path(__file__).resolve().parent.parent
OUT = HERE / ".run-data" / "continuous"
TARGETS = (200, 800, 1500)
PARAMETERS = ("cosmological_parameters--h0", "cosmological_parameters--omega_m",
              "cosmological_parameters--omega_k", "cosmological_parameters--w")
RUNS = {
    "DES-Y3 BAO": ("continuous_des_y3_bao.ini", "des_y3_bao.txt", PARAMETERS),
    "Pantheon+SH0ES": ("continuous_pantheon.ini", "pantheonplussh0es.txt",
                       (*PARAMETERS, "supernova_params--m")),
}


def validate_physics():
    rows = []
    for h0, omega_m, omega_k, w in product(
        (0.45, 1.10), (0.20, 0.40), (-0.01, 0.01), (-1.0, -0.333334)
    ):
        model, omega_r, omega_de = build_cosmology(
            h0 * 100, omega_m, omega_k, w, 2.7255, 3.046
        )
        redshift = np.linspace(0, 3, 301)
        e2 = np.asarray(model.efunc(redshift)) ** 2
        closure = float(model.Om0 + model.Ogamma0 + model.Onu0 + model.Ode0 + model.Ok0)
        q0 = 0.5 * omega_m + omega_r + 0.5 * (1 + 3 * w) * omega_de
        row = {"h0": h0, "omega_m": omega_m, "omega_k": omega_k, "w": w,
               "omega_r": omega_r, "omega_de": omega_de, "closure_sum": closure,
               "min_E2_z0_to_3": float(e2.min()), "q0": q0}
        if not np.isclose(closure, 1, atol=2e-8) or e2.min() <= 0:
            raise RuntimeError(f"Unphysical prior corner: {row}")
        rows.append(row)
    (OUT / "physics_validation.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def fit(label, setting):
    ini_name, chain_name, sampled = setting
    ini = HERE / "configs" / ini_name
    chain_path = OUT / chain_name
    audit_path = OUT / f"{chain_path.stem}_convergence.json"
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {"dataset": label, "checks": []}
    checks = audit["checks"]
    if chain_path.exists():
        base.validate_chain_settings(chain_path, ini)
        existing, _ = base._fit_chain(chain_path, label, sampled)
        current_steps = existing.shape[0]
    else:
        current_steps = 0
    if len(checks) > len(TARGETS):
        raise RuntimeError(f"More than three previous autocorrelation checks for {label}")
    for target in TARGETS:
        if any(record["converged"] for record in checks):
            break
        if len(checks) >= 3 or any(record["target_steps"] == target for record in checks):
            continue
        if current_steps > target:
            raise RuntimeError(f"Existing {label} chain has unrecorded steps; inspect it before continuing")
        if current_steps < target:
            command = [str(base.COSMOSIS), str(ini), "-p", f"emcee.samples={target}"]
            print(f"{label}: extending the same chain to {target:,} steps/walker", flush=True)
            subprocess.run(command, cwd=base.ROOT, env={**os.environ, "PWD": str(base.ROOT)}, check=True)
        data, indices = base._fit_chain(chain_path, label, sampled)
        current_steps = data.shape[0]
        record = base.convergence_record(data, indices)
        record["target_steps"] = target
        record["steps_per_walker"] = current_steps
        checks.append(record)
        audit_path.write_text(json.dumps(audit, indent=2) + "\n")
        print(f"{label}: converged={record['converged']}; "
              f"minimum steps/tau={min(record['steps_per_tau'] or [0]):.1f}", flush=True)
        if record["converged"]:
            break
    if not checks:
        raise RuntimeError(f"No autocorrelation check was recorded for {label}")
    # Each resumed stage contains every earlier sample. The longest checked
    # chain uses the most information, irrespective of its inferred H0 value.
    best = checks[-1]
    audit["best_target_steps"] = best["target_steps"]
    audit["converged"] = best["converged"]
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def summarize(label, chain_name, audit):
    header, chain = base.read_chain(OUT / chain_name)
    samples = chain[int(0.3 * len(chain)):].reshape(-1, len(header))
    get = lambda name: samples[:, header.index("cosmological_parameters--" + name)]
    h0, om, ok, w = (get(name) for name in ("h0", "omega_m", "omega_k", "w"))
    h0_kms = h0 * 100.0
    summary = base.density_summary(h0_kms)
    # Radiation has fixed temperature and massless Neff; its density fraction
    # varies as H0^-2 when H0 is sampled.
    _, radiation_at_70, _ = build_cosmology(70.0, 0.3, 0.0, -1.0, 2.7255, 3.046)
    radiation = radiation_at_70 * (0.7 / h0) ** 2
    dark_energy = 1.0 - om - ok - radiation
    q0 = 0.5 * om + radiation + 0.5 * (1.0 + 3.0 * w) * dark_energy
    if np.any(dark_energy <= 0) or np.max(np.abs(om + ok + radiation + dark_energy - 1)) > 1e-10:
        raise RuntimeError(f"Density closure failed for {label}")
    limits = {"h0": (0.45, 1.10), "omega_m": (0.20, 0.40),
              "omega_k": (-0.01, 0.01), "w": (-1.0, -0.333334)}
    for name, values in zip(limits, (h0, om, ok, w)):
        low, high = limits[name]
        if values.min() < low - 1e-8 or values.max() > high + 1e-8:
            raise RuntimeError(f"{label} {name} chain escaped its prior")
    return {
        "dataset": label,
        "converged": audit["converged"],
        "checks": len(audit["checks"]),
        "steps_per_walker": len(chain),
        "gaussian_peak_h0": summary["mean_and_gaussian_peak"],
        "gaussian_sigma_h0": summary["gaussian_sigma"],
        "sampled_kde_mode_h0": float(summary["plot_x"][np.argmax(summary["plot_kde"])]),
        "h0_credible_16": summary["q16"], "h0_credible_84": summary["q84"],
        "h0_skewness": summary["skewness"],
        "h0_kde_gaussian_l1": summary["kde_gaussian_l1_on_99_6pct_interval"],
        "mean_w": float(np.mean(w)), "mean_omega_m": float(np.mean(om)),
        "mean_omega_k": float(np.mean(ok)), "mean_omega_de": float(np.mean(dark_energy)),
        "mean_omega_r": float(np.mean(radiation)),
        "fraction_accelerating_today": float(np.mean(q0 < 0)),
        "min_effective_sample_size": float(min(audit["checks"][-1]["effective_sample_size"] or [0])),
    }


def write_report(results, audits):
    lines = [
        "# Continuous curved wCDM: DES-Y3 BAO and Pantheon+SH0ES", "",
        "Each dataset has one independently sampled CosmoSIS chain. The continuous "
        "uniform priors are H0 = 45–110 km/s/Mpc, Ωm = 0.20–0.40, Ωk = −0.01–+0.01, "
        "and constant w = −1 to −0.333334. Pantheon+SH0ES also samples the "
        "supernova magnitude M = −21 to −18. Ωb = 0.048 is held fixed and is "
        "part of Ωm, so Ωc = Ωm − Ωb for massless neutrinos.", "",
        "At the present epoch, Ωr + Ωm + Ωk + ΩDE = 1. Thus ΩDE is derived "
        "as 1 − Ωr − Ωm − Ωk for each proposal, rather than sampled separately. "
        "The radiation temperature and massless-neutrino effective number are "
        "fixed at Tcmb = 2.7255 K and Neff = 3.046; the fractional Ωr changes "
        "slightly with sampled H0 because the critical density changes. "
        "CosmoSIS's consistency module writes a variable named Ωlambda from "
        "1 − Ωm − Ωk; that label omits the tiny radiation term and should be "
        "read as a bookkeeping approximation when w differs from −1. The "
        "Pantheon distance module uses exact radiation closure, while CAMB "
        "handles radiation internally for DES-Y3 BAO.", "",
        "The normalized spatial-curvature sign k has the discrete choices "
        "+1, 0, −1, while Ωk = −k c²/(a0² H0²) is continuous in magnitude. "
        "Negative Ωk means positively curved/closed spatial geometry; positive "
        "Ωk means negatively curved/open geometry. The zero value is flat. "
        "Changing Ωk and Ωm adjusts ΩDE through closure. It does not force a "
        "fixed ΩDE value.", "",
        "The condition w < −1/3 gives negative enough pressure for the dark "
        "energy component, but it does not alone guarantee that the total "
        "universe accelerates today. We calculate q0 = Ωm/2 + Ωr + "
        "(1+3w)ΩDE/2 for every saved sample and report its accelerating "
        "fraction. Some prior volume near w = −1/3 decelerates. Expansion "
        "accelerates today only for roughly w < −0.42 at Ωm = 0.20 and "
        "w < −0.56 at Ωm = 0.40 (for negligible radiation and curvature). "
        "Expansion remains real across the positive-density prior box for the redshifts "
        "used here.", "",
        "## Posterior summaries", "",
        "The Gaussian is an approximation built from the marginalized H0 chain: "
        "peak = sample mean and 1σ = sample standard deviation. The figure "
        "overlays this Gaussian on the smoothed sampled posterior. The Gaussian "
        "peak need not be the sampled posterior's mode.", "",
        "| Dataset | Gaussian H0 peak ± 1σ (km/s/Mpc) | 68% equal-tail H0 interval | "
        "Mean w | Mean Ωm | Mean Ωk | Accelerating fraction | Converged? |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in results:
        lines.append(
            f"| {row['dataset']} | {row['gaussian_peak_h0']:.2f} ± "
            f"{row['gaussian_sigma_h0']:.2f} | "
            f"[{row['h0_credible_16']:.2f}, {row['h0_credible_84']:.2f}] | "
            f"{row['mean_w']:.3f} | {row['mean_omega_m']:.3f} | "
            f"{row['mean_omega_k']:+.4f} | "
            f"{row['fraction_accelerating_today']:.1%} | "
            f"{'yes' if row['converged'] else 'no; provisional'} |"
        )
    lines.extend(["", "The sampled KDE modes are " +
                  ", ".join(f"{row['dataset']}: {row['sampled_kde_mode_h0']:.2f} km/s/Mpc"
                            for row in results) + ". "
                  "The BAO posterior is broad and its KDE mode differs substantially "
                  "from the Gaussian mean/peak; its Gaussian should be read only "
                  "as a descriptive approximation.", ""])
    lines.extend(["## Chain checks", ""])
    for label, audit in audits.items():
        lines.append(f"- {label}: {len(audit['checks'])} autocorrelation checks, "
                     f"best at {audit['best_target_steps']:,} steps/walker; "
                     f"criterion is at least 50 steps per integrated "
                     f"autocorrelation time and effective sample size of 1,000 "
                     f"per parameter. Final records are in "
                     f"`{RUNS[label][1].removesuffix('.txt')}_convergence.json`.")
    lines.extend(["", "## Interpretation and files", "",
                  "DES-Y3 BAO uses a sound horizon calculated by CAMB for "
                  "each sampled cosmology, with fixed Ωb, Tcmb, Neff, and "
                  "massless-neutrino settings. Physical densities Ωb h² and "
                  "Ωm h² change as H0 and Ωm vary. Its absolute H0 inference is "
                  "conditional on those early-Universe assumptions. DES-Y3 BAO alone has "
                  "substantial degeneracy among H0, Ωm, Ωk, and w, so assess "
                  "the sampled shape and prior sensitivity before treating "
                  "its Gaussian summary as a precise measurement. The two "
                  "datasets are not multiplied into a joint estimate. This "
                  "two-chain study does not measure sensitivity to alternate "
                  "prior bounds or early-Universe calibration choices.", "",
                  "- `figures/continuous_corner.png`: CosmoSIS "
                  "postprocessor corner figure with sampled marginals and "
                  "joint contours, plus dashed H0 Gaussian approximations.",
                  "- `figures/continuous_h0.png`: enlarged H0 sampled "
                  "posterior and Gaussian comparison with labelled peak lines.",
                  "- `summary.json` and `.csv`: numerical summaries.",
                  "- `physics_validation.json`: all 16 prior-corner checks for "
                  "density closure and positive E(z)² over 0 ≤ z ≤ 3.",
                  "- Two `.txt` files: the independent CosmoSIS chains.", ""])
    (OUT / "legacy_two_dataset_report.md").write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    validate_physics()
    if not base.COSMOSIS.exists():
        raise FileNotFoundError(f"Existing CosmoSIS executable unavailable: {base.COSMOSIS}")
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {label: pool.submit(fit, label, setting) for label, setting in RUNS.items()}
        audits = {label: job.result() for label, job in jobs.items()}
    results = [summarize(label, setting[1], audits[label]) for label, setting in RUNS.items()]
    (OUT / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    with (OUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    write_report(results, audits)
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
