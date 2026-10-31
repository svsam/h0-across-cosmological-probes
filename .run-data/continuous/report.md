# Independent BAO, supernova, and Planck CMB comparison

The four independent CosmoSIS likelihood runs share uniform priors H0 = 45–110 km s⁻¹ Mpc⁻¹, Ωm = 0.20–0.40, Ωk = −0.01–0.01, and constant w = −1 to −0.333334. None of their likelihoods are multiplied together. ΩDE follows Friedmann closure: 1 − Ωr − Ωm − Ωk. Negative Ωk denotes closed spatial geometry and positive Ωk open geometry. Although w < −1/3 gives negative-pressure dark energy, some prior combinations need not accelerate today; the machine-readable summary records the accelerating fraction of each chain.

DESI DR2 BAO uses seven BAO redshift bins (0.295 to 2.330). It contains no BBN observational likelihood. Its absolute H0 scale still depends on CAMB's calculated drag-epoch sound horizon r_d, with Ωb = 0.048, Tcmb = 2.7255 K, Neff = 3.046 and massless neutrinos. DES-Y3 BAO uses the same calibration. Fixed calibration settings are a modeling assumption, not additional BAO data. These BAO-only runs can leave H0 broad.

Pantheon+SH0ES marginalizes over supernova magnitude M = −21 to −18. The Planck comparison uses an approximate compressed CMB likelihood on 100θ* = 1.04109 ± 0.0006, Ωb h² = 0.0224 ± 0.0003, and Ωc h² = 0.120 ± 0.003. These centers come from Planck 2018 TT,TE,EE+lowE; the broadened diagonal errors are study choices and do not reproduce the full covariance. It samples Ωb = 0.03–0.07 in addition to the four shared parameters. This is not a full Planck power-spectrum posterior; the full Planck Lite configuration is available separately but its chain was not used in this comparison. Neutrino mass and Neff remain fixed in all runs.

## Marginalized H0 summaries

**No H0 measurement is claimed:** none of these chains met the agreed autocorrelation criterion in at most three checks. The numbers below are provisional diagnostics of the current samples. The Gaussian peak is the mean of the marginalized H0 chain; its 1σ width is the sample standard deviation. The Gaussian may miss skewness or multiple modes.

| Independent dataset | Gaussian H0 peak ± 1σ (km s⁻¹ Mpc⁻¹) | Mean w | Convergence |
|---|---:|---:|---|
| DES-Y3 BAO | 77.23 ± 16.96 | -0.727 | not passed; provisional |
| Pantheon+SH0ES | 73.36 ± 1.00 | -0.877 | not passed; provisional |
| DESI DR2 BAO | 64.66 ± 3.01 | -0.911 | not passed; provisional |
| Planck 2018 compressed | 63.70 ± 2.84 | -0.851 | not passed; provisional |

## Autocorrelation checks

- DES-Y3 BAO: 3/3 checks; 1500 steps/walker; minimum steps/τ = 21.8, minimum effective samples = 732. Required: ≥50 and ≥1000.
- Pantheon+SH0ES: 3/3 checks; 1500 steps/walker; minimum steps/τ = 24.3, minimum effective samples = 817. Required: ≥50 and ≥1000.
- DESI DR2 BAO: 3/3 checks; 600 steps/walker; minimum steps/τ = 12.4, minimum effective samples = 208. Required: ≥50 and ≥1000.
- Planck 2018 compressed: 3/3 checks; 3000 steps/walker; minimum steps/τ = 18.3, minimum effective samples = 307. Required: ≥50 and ≥1000.

The enlarged H0 figure shows Gaussian curves only, with vertical peak markers and mean ± 1σ bands; its Gaussian parameters are calculated from the marginalized CosmoSIS chains. The corner figure retains sampled marginals and joint contours to show their shape. All parameter labels and ranges appear on the figures. Numerical data are in `summary.json` and `summary.csv`; chain and audit files are alongside them. No new prior-bound or sound-horizon-calibration sensitivity run was made for DESI DR2, so those sensitivities remain unquantified here.

Separate H0 histograms show binned frequencies from each saved CosmoSIS chain; a second histogram grid shows H0, Ωm, Ωk, and w for every pipeline. These use 30% burn and thinning by 10, with density normalization so chains of different lengths can be compared. Each H0 histogram uses 40 bins across that chain's observed range, so its horizontal scale differs by dataset; use the Gaussian comparison for a common H0 axis. Bars show the actual sampled frequencies, while the Gaussian-only comparison remains a distinct figure. The original text chains, convergence audits, summaries, corner plot, and figures are available in `output/continuous/`; standard CosmoSIS marginal statistics and extrema are in `output/continuous/postprocess/`. The native CosmoSIS text tables report its `h0` parameter in units of 100 km s⁻¹ Mpc⁻¹; the figures and study summary convert this to H0 in km s⁻¹ Mpc⁻¹. Histograms and all numerical summaries remain provisional because the chains did not pass convergence.

The additional `continuous_h0_histograms_50_80` figure uses the same 50–80 km s⁻¹ Mpc⁻¹ horizontal axis and 0.5 km s⁻¹ Mpc⁻¹ bins for all four datasets. Each panel states what fraction of its chain lies in this window. Histogram heights are normalized against the entire post-burn chain, so truncating DES-Y3's broad tails does not artificially inflate its visible bars.

## Likelihood sources

- DESI DR2 BAO: https://arxiv.org/abs/2503.14738
- Planck 2018 acoustic scale and physical densities: https://arxiv.org/abs/1807.06209
- CosmoSIS Planck Lite module: https://cosmosis.readthedocs.io/en/latest/reference/standard_library/planck_py.html
