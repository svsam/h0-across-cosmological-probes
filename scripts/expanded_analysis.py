#!/usr/bin/env python3
"""Audit independent DESI BAO and compressed Planck chains."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import continuous_analysis as study
import hubble_analysis as base


OUT = study.OUT
NEW = {
    "DESI DR2 BAO": ("continuous_desi_dr2_bao.ini", "desi_dr2_bao.txt", study.PARAMETERS),
    "Planck 2018 compressed": (
        "continuous_planck_compressed.ini", "planck_compressed.txt",
        (*study.PARAMETERS, "cosmological_parameters--omega_b"),
    ),
}
ALL = {**study.RUNS, **NEW}
TARGETS = {
    "DESI DR2 BAO": (100, 300, 600),
    "Planck 2018 compressed": (100, 300, 3000),
}


def audit_fit(label, setting, *, extend):
    targets = TARGETS[label]
    ini_name, chain_name, sampled = setting
    chain_path = OUT / chain_name
    ini = study.HERE / "configs" / ini_name
    audit_path = OUT / f"{chain_path.stem}_convergence.json"
    if not chain_path.exists() or chain_path.stat().st_size == 0:
        raise FileNotFoundError(f"No chain to audit: {chain_path}")
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {
        "dataset": label, "checks": []}
    checks = audit["checks"]
    if len(checks) > 3:
        raise RuntimeError(f"More than three autocorrelation checks for {label}")
    base.validate_chain_settings(chain_path, ini)
    chain, indices = base._fit_chain(chain_path, label, sampled)
    steps = len(chain)
    if steps not in targets and not (checks and checks[-1]["target_steps"] < steps < targets[-1]):
        raise RuntimeError(f"Unexpected chain length {steps} for {label}; inspect it")
    if steps in targets and not any(item["target_steps"] == steps for item in checks):
        record = base.convergence_record(chain, indices)
        record.update(target_steps=steps, steps_per_walker=steps)
        checks.append(record)
        audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    for target in targets:
        if not extend or checks[-1]["converged"] or len(checks) >= 3:
            break
        if target <= steps:
            continue
        command = [sys.executable, str(study.HERE / "scripts" / "run_cosmosis_fork.py"),
                   str(ini), "--smp", "4", "-p", f"emcee.samples={target}"]
        log_dir = OUT / "logs"
        log_dir.mkdir(exist_ok=True)
        with (log_dir / f"{chain_path.stem}_run.log").open("a") as log:
            subprocess.run(command, cwd=base.ROOT,
                           env={**os.environ, "PWD": str(base.ROOT)},
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        chain, indices = base._fit_chain(chain_path, label, sampled)
        steps = len(chain)
        if steps != target:
            raise RuntimeError(f"Expected {target} steps/walker, found {steps}")
        record = base.convergence_record(chain, indices)
        record.update(target_steps=target, steps_per_walker=target)
        checks.append(record)
        audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    best = checks[-1]
    audit.update(best_target_steps=best["target_steps"], converged=best["converged"])
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def write_report(rows, audits):
    lines = [
        "# Independent BAO, supernova, and Planck CMB comparison", "",
        "The four independent CosmoSIS likelihood runs share uniform priors "
        "H0 = 45–110 km s⁻¹ Mpc⁻¹, Ωm = 0.20–0.40, Ωk = −0.01–0.01, "
        "and constant w = −1 to −0.333334. None of their likelihoods are multiplied "
        "together. ΩDE follows Friedmann closure: 1 − Ωr − Ωm − Ωk. "
        "Negative Ωk denotes closed spatial geometry and positive Ωk open geometry. "
        "Although w < −1/3 gives negative-pressure dark energy, some prior "
        "combinations need not accelerate today; the machine-readable summary "
        "records the accelerating fraction of each chain.", "",
        "DESI DR2 BAO uses seven BAO redshift bins (0.295 to 2.330). It contains "
        "no BBN observational likelihood. Its absolute H0 scale still depends on "
        "CAMB's calculated drag-epoch sound horizon r_d, with Ωb = 0.048, "
        "Tcmb = 2.7255 K, Neff = 3.046 and massless neutrinos. DES-Y3 BAO uses "
        "the same calibration. Fixed calibration settings are a modeling assumption, "
        "not additional BAO data. These BAO-only runs can leave H0 broad.", "",
        "Pantheon+SH0ES marginalizes over supernova magnitude M = −21 to −18. "
        "The Planck comparison uses an approximate compressed CMB likelihood "
        "on 100θ* = 1.04109 ± 0.0006, Ωb h² = 0.0224 ± 0.0003, and "
        "Ωc h² = 0.120 ± 0.003. These centers come from Planck 2018 "
        "TT,TE,EE+lowE; the broadened diagonal errors are study choices and "
        "do not reproduce the full covariance. It samples Ωb = 0.03–0.07 "
        "in addition to the four shared parameters. This is not a full Planck "
        "power-spectrum posterior; the full Planck Lite configuration is "
        "available separately but its chain was not used in this comparison. "
        "Neutrino mass and Neff remain fixed in all runs.", "",
        "## Marginalized H0 summaries", "",
        "**No H0 measurement is claimed:** none of these chains met the agreed "
        "autocorrelation criterion in at most three checks. The numbers below "
        "are provisional diagnostics of the current samples. The Gaussian peak "
        "is the mean of the marginalized H0 chain; its 1σ width is the sample "
        "standard deviation. The Gaussian may miss skewness or multiple modes.", "",
        "| Independent dataset | Gaussian H0 peak ± 1σ (km s⁻¹ Mpc⁻¹) | Mean w | Convergence |",
        "|---|---:|---:|---|",
    ]
    for row in rows:
        lines.append(f"| {row['dataset']} | {row['gaussian_peak_h0']:.2f} ± "
                     f"{row['gaussian_sigma_h0']:.2f} | {row['mean_w']:.3f} | "
                     f"{'passed' if row['converged'] else 'not passed; provisional'} |")
    lines.extend(["", "## Autocorrelation checks", ""])
    for label, audit in audits.items():
        record = audit["checks"][-1]
        ratio = min(record.get("steps_per_tau") or [0])
        ess = min(record.get("effective_sample_size") or [0])
        lines.append(f"- {label}: {len(audit['checks'])}/3 checks; "
                     f"{audit['best_target_steps']} steps/walker; "
                     f"minimum steps/τ = {ratio:.1f}, minimum effective samples "
                     f"= {ess:.0f}. Required: ≥50 and ≥1000.")
    lines.extend(["", "The enlarged H0 figure shows Gaussian curves only, "
                  "with vertical peak markers and mean ± 1σ bands; its Gaussian "
                  "parameters are calculated from the marginalized CosmoSIS chains. "
                  "The corner figure retains sampled marginals and joint contours "
                  "to show their shape. All parameter labels and ranges appear "
                  "on the figures. Numerical data are in `summary.json` and "
                  "`summary.csv`; chain and audit files are alongside them. "
                  "No new prior-bound or sound-horizon-calibration sensitivity "
                  "run was made for DESI DR2, so those sensitivities remain "
                  "unquantified here.", "",
                  "## Likelihood sources", "",
                  "- DESI DR2 BAO: https://arxiv.org/abs/2503.14738",
                  "- Planck 2018 acoustic scale and physical densities: "
                  "https://arxiv.org/abs/1807.06209",
                  "- CosmoSIS Planck Lite module: "
                  "https://cosmosis.readthedocs.io/en/latest/reference/standard_library/planck_py.html", ""])
    (OUT / "report.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extend", action="store_true", help="resume up to three checks")
    args = parser.parse_args()
    audits = {label: json.loads((OUT / f"{setting[1][:-4]}_convergence.json").read_text())
              for label, setting in study.RUNS.items()}
    for label, setting in NEW.items():
        audits[label] = audit_fit(label, setting, extend=args.extend)
    rows = [study.summarize(label, setting[1], audits[label])
            for label, setting in ALL.items()]
    (OUT / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    with (OUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_report(rows, audits)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
