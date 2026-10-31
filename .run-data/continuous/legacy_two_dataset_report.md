# Continuous curved wCDM: DES-Y3 BAO and Pantheon+SH0ES

Each dataset has one independently sampled CosmoSIS chain. The continuous uniform priors are H0 = 45–110 km/s/Mpc, Ωm = 0.20–0.40, Ωk = −0.01–+0.01, and constant w = −1 to −0.333334. Pantheon+SH0ES also samples the supernova magnitude M = −21 to −18. Ωb = 0.048 is held fixed and is part of Ωm, so Ωc = Ωm − Ωb for massless neutrinos.

At the present epoch, Ωr + Ωm + Ωk + ΩDE = 1. Thus ΩDE is derived as 1 − Ωr − Ωm − Ωk for each proposal, rather than sampled separately. The radiation temperature and massless-neutrino effective number are fixed at Tcmb = 2.7255 K and Neff = 3.046; the fractional Ωr changes slightly with sampled H0 because the critical density changes. CosmoSIS's consistency module writes a variable named Ωlambda from 1 − Ωm − Ωk; that label omits the tiny radiation term and should be read as a bookkeeping approximation when w differs from −1. The Pantheon distance module uses exact radiation closure, while CAMB handles radiation internally for DES-Y3 BAO.

The normalized spatial-curvature sign k has the discrete choices +1, 0, −1, while Ωk = −k c²/(a0² H0²) is continuous in magnitude. Negative Ωk means positively curved/closed spatial geometry; positive Ωk means negatively curved/open geometry. The zero value is flat. Changing Ωk and Ωm adjusts ΩDE through closure. It does not force a fixed ΩDE value.

The condition w < −1/3 gives negative enough pressure for the dark energy component, but it does not alone guarantee that the total universe accelerates today. We calculate q0 = Ωm/2 + Ωr + (1+3w)ΩDE/2 for every saved sample and report its accelerating fraction. Some prior volume near w = −1/3 decelerates. Expansion accelerates today only for roughly w < −0.42 at Ωm = 0.20 and w < −0.56 at Ωm = 0.40 (for negligible radiation and curvature). Expansion remains real across the positive-density prior box for the redshifts used here.

## Posterior summaries

The Gaussian is an approximation built from the marginalized H0 chain: peak = sample mean and 1σ = sample standard deviation. The figure overlays this Gaussian on the smoothed sampled posterior. The Gaussian peak need not be the sampled posterior's mode.

| Dataset | Gaussian H0 peak ± 1σ (km/s/Mpc) | 68% equal-tail H0 interval | Mean w | Mean Ωm | Mean Ωk | Accelerating fraction | Converged? |
|---|---:|---:|---:|---:|---:|---:|---|
| DES-Y3 BAO | 77.23 ± 16.96 | [58.52, 97.44] | -0.727 | 0.311 | +0.0001 | 87.3% | no; provisional |
| Pantheon+SH0ES | 73.36 ± 1.00 | [72.36, 74.35] | -0.877 | 0.280 | +0.0004 | 100.0% | no; provisional |

The sampled KDE modes are DES-Y3 BAO: 65.67 km/s/Mpc, Pantheon+SH0ES: 73.13 km/s/Mpc. The BAO posterior is broad and its KDE mode differs substantially from the Gaussian mean/peak; its Gaussian should be read only as a descriptive approximation.

## Chain checks

- DES-Y3 BAO: 3 autocorrelation checks, best at 1,500 steps/walker; criterion is at least 50 steps per integrated autocorrelation time and effective sample size of 1,000 per parameter. Final records are in `des_y3_bao_convergence.json`.
- Pantheon+SH0ES: 3 autocorrelation checks, best at 1,500 steps/walker; criterion is at least 50 steps per integrated autocorrelation time and effective sample size of 1,000 per parameter. Final records are in `pantheonplussh0es_convergence.json`.

## Interpretation and files

DES-Y3 BAO uses a sound horizon calculated by CAMB for each sampled cosmology, with fixed Ωb, Tcmb, Neff, and massless-neutrino settings. Physical densities Ωb h² and Ωm h² change as H0 and Ωm vary. Its absolute H0 inference is conditional on those early-Universe assumptions. DES-Y3 BAO alone has substantial degeneracy among H0, Ωm, Ωk, and w, so assess the sampled shape and prior sensitivity before treating its Gaussian summary as a precise measurement. The two datasets are not multiplied into a joint estimate. This two-chain study does not measure sensitivity to alternate prior bounds or early-Universe calibration choices.

- `figures/continuous_corner.png`: CosmoSIS postprocessor corner figure with sampled marginals and joint contours, plus dashed H0 Gaussian approximations.
- `figures/continuous_h0.png`: enlarged H0 sampled posterior and Gaussian comparison with labelled peak lines.
- `summary.json` and `.csv`: numerical summaries.
- `physics_validation.json`: all 16 prior-corner checks for density closure and positive E(z)² over 0 ≤ z ≤ 3.
- Two `.txt` files: the independent CosmoSIS chains.
