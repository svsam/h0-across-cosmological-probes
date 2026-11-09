#!/usr/bin/env python3
"""Run the two H0 fits, check convergence, and write comparison artifacts."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde, norm, skew


HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
CONFIGS = HERE / "configs"
_ENV_COSMOSIS = Path(sys.executable).with_name("cosmosis")
COSMOSIS = _ENV_COSMOSIS if _ENV_COSMOSIS.is_file() else Path(shutil.which("cosmosis") or "cosmosis")
IGNORED_RUN_SETTINGS = {
    ("runtime", "root"), ("pipeline", "values"), ("pipeline", "priors"),
    ("output", "filename"), ("test", "save_dir"),
    ("emcee", "samples"), ("emcee", "nsteps"),
}
MAX_STEPS = 200_000
GROWTH_STEPS = (2_000, 2_500, 3_000, 3_500, 4_000, 5_000, 5_500, 6_000, 6_500, 7_000,
                7_500, 10_000, 10_500, 11_000, 12_000, 20_000, 40_000, 80_000,
                120_000, 160_000, 200_000)
PARAMS = {
    "Pantheon+SH0ES": ("pantheon_wcdm.ini", "pantheon_wcdm.txt",
                       ("cosmological_parameters--h0", "cosmological_parameters--w",
                        "supernova_params--m")),
    "DES-Y3 BAO": ("des_y3_bao_wcdm.ini", "des_y3_bao_wcdm.txt",
                   ("cosmological_parameters--h0", "cosmological_parameters--w")),
}


def _chain_preamble(path: Path):
    """Read CosmoSIS metadata and its embedded input files before the samples."""
    columns = None
    metadata = {}
    blocks = {name: [] for name in ("PARAMS", "VALUES", "PRIORS")}
    active_block = None
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            content = line.lstrip("# ").strip()
            if content.startswith("START_OF_") and content.endswith("_INI"):
                active_block = content.removeprefix("START_OF_").removesuffix("_INI")
            elif content.startswith("END_OF_") and content.endswith("_INI"):
                active_block = None
            elif active_block in blocks:
                blocks[active_block].append(line.removeprefix("##").strip())
            elif line.startswith("##"):
                continue
            elif "=" in content:
                key, value = content.split("=", 1)
                metadata[key] = value
            else:
                candidate = content.split()
                if "post" in candidate and any("--" in name for name in candidate):
                    columns = candidate
    if columns is None:
        raise RuntimeError(f"No CosmoSIS chain header in {path}")
    return columns, metadata, blocks


def read_chain(path: Path):
    header, metadata, _ = _chain_preamble(path)
    try:
        walkers = int(metadata["walkers"])
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"Missing or invalid walker count in {path}") from error
    if walkers < 1:
        raise RuntimeError(f"Invalid walker count in {path}: {walkers}")
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != len(header) or values.shape[0] % walkers:
        raise RuntimeError(f"Malformed chain or row count not divisible by {walkers}: {path}")
    return header, values.reshape((-1, walkers, len(header)))


def _read_config(source: Path | list[str], *, preserve_case: bool = False):
    config = configparser.ConfigParser(interpolation=None)
    config.optionxform = str if preserve_case else str.lower
    if isinstance(source, Path):
        if not config.read(source):
            raise FileNotFoundError(source)
    else:
        config.read_string("\n".join(source) or "[empty]\n")
        config.remove_section("empty")
    return config


def _resolved_config_path(value: str):
    """CosmoSIS expands ${PWD} relative to the repository run directory."""
    path = Path(value.replace("${PWD}", str(ROOT))).expanduser()
    return path if path.is_absolute() else ROOT / path


def _settings(config, *, omit=()):
    ignored = set(omit)
    settings = {}
    for section in config.sections():
        for option, raw_value in config.items(section):
            key = (section.lower(), option.lower())
            if key in ignored:
                continue
            value = " ".join(raw_value.split())
            if option.lower() == "file":
                path = _resolved_config_path(value)
                value = str(path.resolve())
            settings[key] = value
    return settings


def _expected_chain_settings(ini_path: Path, overrides: dict[str, str] | None):
    params = _read_config(ini_path)
    for key, value in (overrides or {}).items():
        section, option = key.split(".", 1)
        params[section][option] = value
    values = _read_config(_resolved_config_path(params["pipeline"]["values"]))
    priors_path = params["pipeline"].get("priors")
    priors = _read_config(_resolved_config_path(priors_path)) if priors_path else _read_config([])
    return params, values, priors


def validate_chain_settings(path: Path, ini_path: Path, overrides: dict[str, str] | None = None):
    """Reject old chains before CosmoSIS can resume them under a new model."""
    _, metadata, blocks = _chain_preamble(path)
    expected = _expected_chain_settings(ini_path, overrides)
    embedded = tuple(_read_config(blocks[name]) for name in ("PARAMS", "VALUES", "PRIORS"))
    for name, current, previous in zip(("run", "values", "priors"), expected, embedded):
        ignored = IGNORED_RUN_SETTINGS if name == "run" else ()
        current_settings = _settings(current, omit=ignored)
        previous_settings = _settings(previous, omit=ignored)
        if current_settings != previous_settings:
            differences = sorted(key for key in current_settings.keys() | previous_settings.keys()
                                 if current_settings.get(key) != previous_settings.get(key))
            changed = ", ".join(f"{section}.{option}" for section, option in differences)
            raise RuntimeError(f"Existing chain {path} has different {name} settings ({changed}); move it aside before rerunning")
    if metadata.get("sampler", "").lower() != expected[0]["runtime"]["sampler"].lower():
        raise RuntimeError(f"Existing chain {path} uses a different sampler; move it aside before rerunning")


def _file_hash(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_fingerprint(ini_path: Path, overrides: dict[str, str] | None):
    """Track code and likelihood data that a chain header cannot identify."""
    params, values, priors = _expected_chain_settings(ini_path, overrides)
    settings = [_settings(params, omit=IGNORED_RUN_SETTINGS), _settings(values), _settings(priors)]
    modules = params["pipeline"]["modules"].split()
    source_files = {_resolved_config_path(params[module]["file"]) for module in modules}
    if "pantheon_plus" in modules:
        pantheon = ROOT / "likelihood" / "pantheon_plus"
        options = params["pantheon_plus"]
        source_files.update((
            pantheon / "pantheon_covariance_io.py",
            _resolved_config_path(options.get("data_file", str(pantheon / "Pantheon+SH0ES.dat"))),
            _resolved_config_path(options.get("covmat_file", str(pantheon / "Pantheon+SH0ES_STAT+SYS.cov_compressed.gz"))),
        ))
    if "des_y3_bao" in modules:
        source_files.add(ROOT / "likelihood" / "des-y3-bao" /
                         "chi2profile_dvdesy3_cosmoplanck18_covcosmolike.csv")
    payload = {
        "settings": [{f"{section}.{option}": value for (section, option), value in group.items()}
                     for group in settings],
        "files": {str(path.resolve()): _file_hash(path) for path in sorted(source_files)},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def convergence_record(chain: np.ndarray, columns: list[int]):
    import emcee
    from emcee.autocorr import AutocorrError, integrated_time

    try:
        chain = chain[:, :, columns]
        tau = np.asarray(integrated_time(chain, tol=0, quiet=True), dtype=float)
    except (AutocorrError, ValueError) as error:
        return {"converged": False, "reason": str(error), "tau": None,
                "steps_per_tau": None, "effective_sample_size": None}
    tau = np.atleast_1d(tau)
    steps_per_tau = chain.shape[0] / tau
    post_burn = chain[int(0.3 * len(chain)):]
    ess = (post_burn.shape[0] * chain.shape[1]) / tau
    enough = bool(np.all(np.isfinite(tau)) and np.all(tau > 0)
                  and np.all(steps_per_tau >= 50) and np.all(ess >= 1000))
    return {"converged": enough,
            "emcee_version": emcee.__version__,
            "tau": tau.tolist(),
            "steps_per_tau": steps_per_tau.tolist(),
            "effective_sample_size": ess.tolist(),
            "criterion": "at least 50 steps per integrated autocorrelation time and ESS >= 1000 per parameter"}


def _fit_chain(path: Path, label: str, sampled_params: tuple[str, ...]):
    columns, chain = read_chain(path)
    missing = [name for name in (*sampled_params, "post") if name not in columns]
    if missing:
        raise RuntimeError(f"{label} chain lacks columns: {missing}")
    indices = [columns.index(name) for name in sampled_params]
    if not np.isfinite(chain[:, :, [*indices, columns.index("post")]]).all():
        raise RuntimeError(f"Non-finite samples or posterior values in {label} chain")
    return chain, indices


def run_fit(label: str, ini_name: str, output_name: str, sampled_params: tuple[str, ...],
            overrides: dict[str, str] | None = None, output_subdir: str = ".run-data",
            ini_path: Path | None = None, output_path: Path | None = None):
    out_dir = HERE / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    output = output_path or (out_dir / output_name)
    base_ini = ini_path or (CONFIGS / ini_name)
    inputs_file = output.with_suffix(output.suffix + ".inputs.json")
    source_fingerprint = _source_fingerprint(base_ini, overrides)
    new_chain = not output.exists()
    history = []
    converged = False
    existing_steps = 0
    if not new_chain:
        validate_chain_settings(output, base_ini, overrides)
        if inputs_file.exists():
            saved = json.loads(inputs_file.read_text())
            if saved.get("sha256") != source_fingerprint:
                raise RuntimeError(f"Source code or likelihood data changed since {output} was generated; move the chain aside before rerunning")
        previous_chain, indices = _fit_chain(output, label, sampled_params)
        existing_steps = previous_chain.shape[0]
        walkers = previous_chain.shape[1]
        if existing_steps >= 2000:
            current = convergence_record(previous_chain, indices)
            current["steps_per_walker"] = int(existing_steps)
            history.append(current)
            converged = current["converged"]

    for target in (() if converged else GROWTH_STEPS):
        if target > MAX_STEPS:
            break
        if target <= existing_steps:
            continue
        command = [str(COSMOSIS), str(base_ini), "-p", f"emcee.samples={target}"]
        if overrides:
            for key, value in overrides.items():
                command.extend(["-p", f"{key}={value}"])
        subprocess.run(command, cwd=ROOT, check=True)
        validate_chain_settings(output, base_ini, overrides)
        if _source_fingerprint(base_ini, overrides) != source_fingerprint:
            raise RuntimeError(f"Source code or likelihood data changed while generating {output}")
        chain, indices = _fit_chain(output, label, sampled_params)
        walkers = chain.shape[1]
        if new_chain:
            inputs_file.write_text(json.dumps({"sha256": source_fingerprint}, indent=2) + "\n")
            new_chain = False
        record = convergence_record(chain, indices)
        record["steps_per_walker"] = int(chain.shape[0])
        history.append(record)
        print(f"{label}: {chain.shape[0]} steps/walker; convergence={record['converged']}", flush=True)
        if record["converged"]:
            converged = True
            break

    result = {"dataset": label, "chain": str(output), "converged": converged,
              "walkers": walkers, "history": history}
    (out_dir / f"{output.stem}_convergence.json").write_text(json.dumps(result, indent=2) + "\n")
    if not converged:
        raise RuntimeError(f"{label} did not meet convergence checks by {MAX_STEPS} steps per walker")
    return result


def posterior_samples(path: Path):
    """Return finite post-burn H0 samples in km/s/Mpc from a saved chain."""
    header, chain = read_chain(path)
    burn = int(chain.shape[0] * 0.3)
    h0 = chain[burn:, :, header.index("cosmological_parameters--h0")].reshape(-1) * 100.0
    return h0[np.isfinite(h0)]


def density_summary(samples: np.ndarray, *, plot_points: int = 512):
    """Summarize H0 samples and tabulate KDE/Gaussian curves for plotting.

    The plot grid spans the central 99.6% of samples. The Gaussian uses the
    sample mean and standard deviation as an approximation to the posterior.
    """
    if len(samples) < 100:
        raise RuntimeError("Too few post-burn-in samples for a posterior summary")
    mu, sd = float(np.mean(samples)), float(np.std(samples, ddof=1))
    q16, median, q84 = map(float, np.quantile(samples, [0.16, 0.5, 0.84]))
    x = np.linspace(float(np.quantile(samples, .002)), float(np.quantile(samples, .998)), plot_points)
    kde = gaussian_kde(samples)
    gaussian = norm.pdf(x, loc=mu, scale=sd)
    empirical = kde(x)
    return {"mean_and_gaussian_peak": mu, "gaussian_sigma": sd, "median": median,
            "q16": q16, "q84": q84, "n_post_burn": int(len(samples)),
            "skewness": float(skew(samples)),
            "kde_gaussian_l1_on_99_6pct_interval": float(np.trapezoid(np.abs(empirical-gaussian), x)),
            "plot_x": x, "plot_kde": empirical, "plot_gaussian": gaussian}


def write_report(summaries: dict[str, dict], convergence: list[dict], out_dir: Path):
    lines = ["# DES-Y3 BAO and Pantheon+SH0ES: Hubble-constant comparison", "",
             "Both likelihoods are fit independently in flat constant-\\(w\\) CDM with \\Omega_m=0.3 fixed. The sampled parameters are \\(H_0\\) and \\(w\\); Pantheon+SH0ES also samples its absolute-magnitude nuisance parameter. DES-Y3 uses the standard early-Universe CAMB sound-horizon calculation with fixed input densities, so its absolute \\(H_0\\) posterior is conditional on that ruler calibration.", "",
             "Inputs are read from the shared standard-library DES-Y3 BAO chi-square profile and Pantheon+SH0ES data/covariance; no likelihood data are copied into this project. The CAMB calibration fixes \\(\\omega_b=0.048\\), \\(A_s=2.19\\times10^{-9}\\), \\(n_s=0.97\\), \\(N_{eff}=3.046\\), and the standard photon temperature.", "",
             "The Gaussian is an approximation to each marginalized posterior, assessed with skewness and KDE-to-Gaussian distance.", "",
             "## Main physical-prior results", "",
             "| Dataset | Gaussian peak H₀ | 68% posterior interval | Gaussian σ | Skewness | KDE/Gaussian L1 | Samples |", "|---|---:|---:|---:|---:|---:|---:|"]
    for name, summary in summaries.items():
        lines.append(f"| {name} | {summary['mean_and_gaussian_peak']:.2f} | {summary['q16']:.2f}–{summary['q84']:.2f} (median {summary['median']:.2f}) | {summary['gaussian_sigma']:.2f} | {summary['skewness']:.3f} | {summary['kde_gaussian_l1_on_99_6pct_interval']:.3f} | {summary['n_post_burn']:,} |")
    lines += ["", "## Assumptions and interpretation", "",
              "- Main prior: \\( -1 \\leq w < -1/3 \\). A separate sensitivity run uses \\( -1.5 \\leq w < -1/3 \\); its phantom part is phenomenological and is not treated as a canonical scalar-field model. The sampled \\(H_0\\) prior is 50–95 km s⁻¹ Mpc⁻¹.",
              "- DES-Y3 BAO constrains a distance-to-sound-horizon ratio. The quoted BAO \\(H_0\\) result is conditional on the fixed standard early-Universe inputs used to calculate \\(r_d\\). It is not an independent, ruler-free absolute measurement.",
              "- Pantheon+SH0ES includes the SH0ES Cepheid calibration; the fitted constraint is not an independent replication of SH0ES.",
              "- These are separate posteriors, not a combined estimate. The fixed \\(\\Omega_m\\) choice and the \\(w\\) prior affect their widths and locations.", "",
              "## Workflow tree", "",
              "```mermaid",
              "flowchart TD",
              "    A[Shared standard-library data] --> B{Independent likelihood fits}",
              "    B --> C[DES-Y3 BAO chi-square profile + CAMB distances and sound horizon]",
              "    B --> D[Pantheon+SH0ES data + compressed statistical/systematic covariance]",
              "    C --> E[Flat wCDM: sample H0 and w; fix Omega_m = 0.3]",
              "    D --> E",
              "    E --> F[Check integrated autocorrelation and ESS]",
              "    F --> G[Marginalize H0; summarize mean, sigma, and 68% interval]",
              "    G --> H[Report separate posterior summaries]",
              "    G --> I[Repeat broad-w and +/-1% rd sensitivity fits]",
              "    H --> J[Report separate constraints and calibration dependence]",
              "    I --> J",
              "```", "",
              "## Sampling checks", ""]
    for item in convergence:
        final = item["history"][-1]
        lines.append(f"- {item['dataset']}: converged with {final['steps_per_walker']:,} steps/walker; steps/\\(\\tau\\) = {', '.join(f'{x:.1f}' for x in final['steps_per_tau'])}; ESS = {', '.join(f'{x:.0f}' for x in final['effective_sample_size'])}.")
    lines += ["", "## Gaussian approximation and sensitivity", "",
              "The DES-Y3 BAO posterior is much broader than Pantheon+SH0ES and its KDE/Gaussian L1 distance is comparatively large, so its Gaussian summary is not a complete description. Use the separate 68% interval when reading its uncertainty. Flat-model broad-\\(w\\) and ±1% \\(r_d\\) sensitivity results are in `../sensitivity/README.md`; the phantom-inclusive prior is phenomenological, and the ruler shifts quantify the BAO calibration dependence.", ""]
    (out_dir / "README.md").write_text("\n".join(lines))


def run_sensitivity():
    """Run a broad-w prior and +/-1% DES ruler calibration stress check."""
    runs = []
    sensitivity_summaries = {}
    for label, (ini, _, sampled) in PARAMS.items():
        cfg = _read_config(CONFIGS / ini, preserve_case=True)
        values_file = CONFIGS / ("pantheon_wcdm_values.ini" if label == "Pantheon+SH0ES" else "des_y3_bao_wcdm_values.ini")
        vals = _read_config(values_file, preserve_case=True)
        vals["cosmological_parameters"]["w"] = "-1.5 -0.7 -0.333333333333"
        with tempfile.TemporaryDirectory(prefix="h0_sensitivity_") as temp_name:
            temp = Path(temp_name)
            vals_path = temp / "values.ini"
            vals_path.write_text(_config_text(vals))
            out = HERE / ".run-data" / "sensitivity"
            out.mkdir(parents=True, exist_ok=True)
            chain_path = out / f"{slug(label)}_phantom_w.txt"
            cfg["pipeline"]["values"] = str(vals_path)
            cfg["output"]["filename"] = str(chain_path)
            cfg_path = temp / "model.ini"
            cfg_path.write_text(_config_text(cfg))
            result = run_fit(label + " broad-w sensitivity", ini, chain_path.name, sampled,
                             output_subdir=".run-data/sensitivity", ini_path=cfg_path,
                             output_path=chain_path)
            runs.append(result)
            sensitivity_summaries[label + " broad-w"] = density_summary(posterior_samples(chain_path))

    # DES ruler stress checks use the same samples/configuration with rd_scale varied.
    # The BAO likelihood is smooth over this small perturbation, so these are
    # rerun as separate posterior chains and never mixed into the main overlay.
    for scale in (0.99, 1.01):
        with tempfile.TemporaryDirectory(prefix="h0_rd_sensitivity_") as temp_name:
            temp = Path(temp_name)
            cfg = _read_config(CONFIGS / "des_y3_bao_wcdm.ini", preserve_case=True)
            out = HERE / ".run-data" / "sensitivity"
            out.mkdir(parents=True, exist_ok=True)
            chain_path = out / f"des_y3_bao_rdscale_{scale:.2f}.txt"
            cfg["output"]["filename"] = str(chain_path)
            cfg["des_y3_bao"]["rd_scale"] = str(scale)
            cfg_path = temp / "des_y3_bao.ini"
            cfg_path.write_text(_config_text(cfg))
            result = run_fit(f"DES-Y3 BAO r_d scale {scale:.2f}", "des_y3_bao_wcdm.ini",
                             chain_path.name, PARAMS["DES-Y3 BAO"][2],
                             output_subdir=".run-data/sensitivity", ini_path=cfg_path,
                             output_path=chain_path)
            runs.append(result)
            sensitivity_summaries[f"DES-Y3 BAO r_d × {scale:.2f}"] = density_summary(posterior_samples(chain_path))
    serializable = {name: {k: v for k, v in summary.items() if not k.startswith("plot_")}
                    for name, summary in sensitivity_summaries.items()}
    out_path = HERE / ".run-data" / "sensitivity" / "sensitivity_summary.json"
    out_path.write_text(json.dumps({"runs": runs, "summaries": serializable}, indent=2) + "\n")
    base = json.loads((HERE / ".run-data" / "flat_wcdm" / "posterior_summary.json").read_text())
    report = ["", "## Sensitivity runs", "",
              "| Variant | Gaussian peak H₀ | 68% interval | Shift from corresponding baseline |", "|---|---:|---:|---:|"]
    for name, summary in sensitivity_summaries.items():
        base_name = "Pantheon+SH0ES" if name.startswith("Pantheon") else "DES-Y3 BAO"
        delta = summary["mean_and_gaussian_peak"] - base[base_name]["mean_and_gaussian_peak"]
        report.append(f"| {name} | {summary['mean_and_gaussian_peak']:.2f} | {summary['q16']:.2f}–{summary['q84']:.2f} | {delta:+.2f} |")
    report.append("")
    with (HERE / ".run-data" / "sensitivity" / "README.md").open("w") as stream:
        stream.write("\n".join(report))


def _config_text(config):
    from io import StringIO
    stream = StringIO()
    config.write(stream)
    return stream.getvalue()


def slug(label):
    return label.lower().replace("+", "plus").replace(" ", "_").replace("-", "_")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sensitivity", action="store_true", help="also run broad-w and r_d sensitivity chains")
    args = parser.parse_args()
    if not COSMOSIS.exists():
        raise FileNotFoundError(f"Existing CosmoSIS executable not found: {COSMOSIS}")
    results = []
    summaries = {}
    for label, (ini, chain, sampled) in PARAMS.items():
        results.append(run_fit(label, ini, chain, sampled,
                               output_subdir=".run-data/flat_wcdm"))
        summaries[label] = density_summary(posterior_samples(HERE / ".run-data" / "flat_wcdm" / chain))
    out_dir = HERE / ".run-data" / "flat_wcdm"
    write_report(summaries, results, out_dir)
    serializable = {name: {k: v for k, v in summary.items() if not k.startswith("plot_")}
                    for name, summary in summaries.items()}
    (out_dir / "posterior_summary.json").write_text(json.dumps(serializable, indent=2) + "\n")
    if args.sensitivity:
        run_sensitivity()
    print(f"Wrote {out_dir / 'README.md'}")


if __name__ == "__main__":
    main()
