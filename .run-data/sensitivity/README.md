# Flat-model sensitivity results

These runs extend the flat \(w\)CDM comparison summarized in [`../flat_wcdm/README.md`](../flat_wcdm/README.md). They were fit independently; the BAO results remain conditional on the sound-horizon calibration. The broad prior includes phantom values and is phenomenological.

| Variant | Gaussian peak \(H_0\) | 68% interval | Gaussian σ | Converged |
|---|---:|---:|---:|:---:|
| Pantheon+SH0ES, \(w\in[-1.5,-1/3]\) | 73.48 | 72.45–74.50 | 1.02 | yes |
| DES-Y3 BAO, \(w\in[-1.5,-1/3]\) | 74.61 | 59.66–88.93 | 12.52 | yes |
| DES-Y3 BAO, \(r_d\times0.99\) | 72.97 | 59.46–86.81 | 11.92 | yes |
| DES-Y3 BAO, \(r_d\times1.01\) | 71.28 | 57.54–85.47 | 12.14 | yes |

Machine-readable values and convergence records are in [`sensitivity_summary.json`](sensitivity_summary.json); the chain files are retained in this folder.
