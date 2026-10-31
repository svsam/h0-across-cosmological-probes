# Curved phantom wCDM: DES-Y3 BAO and Pantheon+SH0ES

This study fits DES-Y3 BAO and Pantheon+SH0ES independently at fixed curvature values Ωk = −0.01, 0, +0.01. Each run samples H0 and constant w ∈ [−1.5, −1.3], fixes Ωm = 0.30, and derives ΩDE from density closure. This full w range is phantom and is presented as phenomenological, not as a canonical scalar-field model.

Pantheon distances use photons at Tcmb = 2.7255 K and Neff = 3.046 massless neutrinos, with Ωm + Ωr + ΩDE + Ωk = 1. DES-Y3 uses CAMB with the study's fixed early-universe settings to calculate rd; its H0 posterior therefore remains conditional on that sound-horizon calibration.

The tabulated Gaussian peak is the mean of marginalized H0 samples; its 1σ is their standard deviation. It need not coincide with the mode of the sampled posterior. In the CosmoSIS figures, solid curves show smoothed sampled posteriors and dashed curves show Gaussian summaries; the horizontal H0 axis uses km/s/Mpc. Because w is sampled continuously over [−1.5, −1.3], each table reports its posterior mean rather than a fixed w value.

## Results

### Ωk = -0.01

| Dataset | Posterior mean w | Gaussian peak H0 (km/s/Mpc) | Gaussian 1σ (km/s/Mpc) |
|---|---:|---:|---:|
| Pantheon+SH0ES | -1.306 | 75.10 | 1.03 |
| DES-Y3 BAO | -1.354 | 88.64 | 5.02 |

### Ωk = +0.00

| Dataset | Posterior mean w | Gaussian peak H0 (km/s/Mpc) | Gaussian 1σ (km/s/Mpc) |
|---|---:|---:|---:|
| Pantheon+SH0ES | -1.306 | 75.05 | 1.04 |
| DES-Y3 BAO | -1.357 | 88.65 | 5.03 |

### Ωk = +0.01

| Dataset | Posterior mean w | Gaussian peak H0 (km/s/Mpc) | Gaussian 1σ (km/s/Mpc) |
|---|---:|---:|---:|
| Pantheon+SH0ES | -1.306 | 74.98 | 1.04 |
| DES-Y3 BAO | -1.358 | 88.37 | 5.13 |

## Convergence

- Pantheon+SH0ES Ωk=+0.00: converged at 10,000 steps/walker; steps/τ = 232.2, 156.1, 230.9; split-Rhat = 1.0058, 1.0070, 1.0057; ESS = 5202, 3496, 5172.
- Pantheon+SH0ES Ωk=+0.01: converged at 10,000 steps/walker; steps/τ = 218.1, 168.2, 218.6; split-Rhat = 1.0076, 1.0067, 1.0074; ESS = 4886, 3769, 4897.
- Pantheon+SH0ES Ωk=-0.01: converged at 12,000 steps/walker; steps/τ = 265.8, 159.2, 262.8; split-Rhat = 1.0035, 1.0080, 1.0036; ESS = 5955, 3566, 5886.
- DES-Y3 BAO Ωk=+0.00: converged at 7,500 steps/walker; steps/τ = 179.5, 169.9; split-Rhat = 1.0095, 1.0071; ESS = 4020, 3805.
- DES-Y3 BAO Ωk=-0.01: converged at 10,000 steps/walker; steps/τ = 233.1, 223.8; split-Rhat = 1.0032, 1.0070; ESS = 5222, 5012.
- DES-Y3 BAO Ωk=+0.01: converged at 10,000 steps/walker; steps/τ = 218.0, 216.1; split-Rhat = 1.0049, 1.0084; ESS = 4884, 4841.

## Physics checks

A prior-corner grid checked the present-day density closure, positive E(z)² from z=0 to 3, and acceleration today (q0 < 0). The validation values are in `physics_validation.json`.

## Figures and data

- [`native_postprocess/corner/h0_w_corner.png`](native_postprocess/corner/h0_w_corner.png): H0–w corner comparison for both datasets; color distinguishes the three fixed Ωk runs. The contours use samples thinned by ten after burn-in.
- [`native_postprocess/`](native_postprocess/README.md): CosmoSIS-native H0 and w posterior PNG curves, with separate comparisons for each Ωk.
- `curved_summary.json` and `curved_summary.csv`: machine-readable posterior summaries.
- The six CosmoSIS text chains and convergence JSON files are retained alongside the figures.

The DES-Y3 BAO sampled H0 posterior rises toward its upper prior bound, H0 = 95 km/s/Mpc. Its Gaussian mean and standard deviation summarize a markedly non-Gaussian, prior-limited distribution; they should not be interpreted as the actual posterior mode and width.
