# BidiagonalSVD sweep statistics (2026-06-24)
Post-processed from `docs/sweeps_2026-06-18/*.log` — no reruns.
New per-row column **`ratio = t_eval / t_dbdsqr`** (solver wall time vs DBDSQR σ-only reference).
Per-row CSVs (with the new ratio column) live in `docs/sweeps_2026-06-18/stats/<suite>_<method>.csv`.

Stats key per column: **min / mean / max** over rows where the value is finite. Two row sets per suite × method:
- **ALL** rows that matched the row schema (includes FAIL inf/nan, which are excluded from each column's min/mean/max via a finite mask).
- **PASS** rows only.

**Caveat — synth coverage**: the orchestrator only printed one in every 100 PASS rows for synth (full periodic print would have made the log ~ 18 k lines × 3 methods). All 5 038 / 691 / 3 FAILs are present in the synth advisor / selfcontained / dbdsqr CSVs (parsed from the trailing `Failing (…)` block), so synth PASS stats below are a ~1 % sample but synth FAIL stats are exact. The other three suites are 100 % complete.

## pract

### advisor  — rows=127  PASS=64

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.251e+03 | 9.941e+03 | 127 | 0 |
| n | PASS | 3.000e+00 | 2.205e+02 | 1.919e+03 | 64 | 0 |
| res | ALL | 0.000e+00 | 1.982e-01 | 1.170e+01 | 113 | 14 |
| res | PASS | 0.000e+00 | 1.188e-01 | 1.000e+00 | 64 | 0 |
| ortU | ALL | 0.000e+00 | 7.434e+04 | 8.307e+06 | 113 | 14 |
| ortU | PASS | 0.000e+00 | 1.684e+00 | 5.000e+00 | 64 | 0 |
| ortV | ALL | 0.000e+00 | 7.434e+04 | 8.307e+06 | 113 | 14 |
| ortV | PASS | 0.000e+00 | 1.630e+00 | 5.000e+00 | 64 | 0 |
| drift | ALL | 0.000e+00 | 2.655e-02 | 1.000e+00 | 113 | 14 |
| drift | PASS | 0.000e+00 | 4.688e-02 | 1.000e+00 | 64 | 0 |
| t_eval | ALL | 4.192e-06 | 3.106e-01 | 1.560e+00 | 113 | 14 |
| t_eval | PASS | 4.192e-06 | 4.822e-02 | 1.398e+00 | 64 | 0 |
| t_dbdsqr | ALL | 8.950e-07 | 1.657e-02 | 8.787e-02 | 113 | 14 |
| t_dbdsqr | PASS | 8.950e-07 | 2.737e-03 | 6.425e-02 | 64 | 0 |
| ratio | ALL | 1.117e+00 | 1.051e+02 | 7.331e+03 | 113 | 14 |
| ratio | PASS | 1.117e+00 | 2.856e+01 | 1.290e+03 | 64 | 0 |

### selfcontained  — rows=127  PASS=114

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.251e+03 | 9.941e+03 | 127 | 0 |
| n | PASS | 3.000e+00 | 9.798e+02 | 9.941e+03 | 114 | 0 |
| res | ALL | 0.000e+00 | 8.339e+14 | 1.051e+17 | 126 | 1 |
| res | PASS | 0.000e+00 | 2.799e-01 | 2.963e+00 | 114 | 0 |
| ortU | ALL | 0.000e+00 | 1.191e+13 | 1.501e+15 | 126 | 1 |
| ortU | PASS | 0.000e+00 | 1.277e+00 | 4.726e+00 | 114 | 0 |
| ortV | ALL | 0.000e+00 | 1.191e+13 | 1.501e+15 | 126 | 1 |
| ortV | PASS | 0.000e+00 | 1.271e+00 | 4.724e+00 | 114 | 0 |
| drift | ALL | 0.000e+00 | 8.653e+09 | 1.060e+12 | 123 | 4 |
| drift | PASS | 0.000e+00 | 9.589e+09 | 1.060e+12 | 111 | 3 |
| t_eval | ALL | 2.921e-03 | 4.103e+00 | 7.310e+01 | 126 | 1 |
| t_eval | PASS | 2.921e-03 | 2.711e+00 | 7.310e+01 | 114 | 0 |
| t_dbdsqr | ALL | 3.440e-07 | 5.273e-02 | 1.137e+00 | 123 | 4 |
| t_dbdsqr | PASS | 3.440e-07 | 3.599e-02 | 1.137e+00 | 111 | 3 |
| ratio | ALL | 1.962e+01 | 3.757e+02 | 8.491e+03 | 123 | 4 |
| ratio | PASS | 1.962e+01 | 3.967e+02 | 8.491e+03 | 111 | 3 |

### dbdsqr  — rows=127  PASS=114

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.251e+03 | 9.941e+03 | 127 | 0 |
| n | PASS | 3.000e+00 | 7.663e+02 | 2.873e+03 | 114 | 0 |
| res | ALL | 0.000e+00 | 3.263e-01 | 4.100e+00 | 114 | 13 |
| res | PASS | 0.000e+00 | 3.263e-01 | 4.100e+00 | 114 | 0 |
| ortU | ALL | 0.000e+00 | 1.781e-01 | 8.000e-01 | 114 | 13 |
| ortU | PASS | 0.000e+00 | 1.781e-01 | 8.000e-01 | 114 | 0 |
| ortV | ALL | 0.000e+00 | 2.044e-01 | 1.400e+00 | 114 | 13 |
| ortV | PASS | 0.000e+00 | 2.044e-01 | 1.400e+00 | 114 | 0 |
| drift | ALL | 0.000e+00 | 8.736e-15 | 1.260e-13 | 114 | 13 |
| drift | PASS | 0.000e+00 | 8.736e-15 | 1.260e-13 | 114 | 0 |
| t_eval | ALL | 1.941e-06 | 1.253e+01 | 9.605e+01 | 114 | 13 |
| t_eval | PASS | 1.941e-06 | 1.253e+01 | 9.605e+01 | 114 | 0 |
| t_dbdsqr | ALL | 8.300e-07 | 1.652e-02 | 9.261e-02 | 114 | 13 |
| t_dbdsqr | PASS | 8.300e-07 | 1.652e-02 | 9.261e-02 | 114 | 0 |
| ratio | ALL | 1.071e+00 | 2.541e+02 | 1.145e+03 | 114 | 13 |
| ratio | PASS | 1.071e+00 | 2.541e+02 | 1.145e+03 | 114 | 0 |

## synth

### advisor  — rows=5176  PASS=138

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.174e+02 | 5.000e+02 | 5176 | 0 |
| n | PASS | 4.000e+00 | 9.408e+01 | 2.610e+02 | 138 | 0 |
| res | ALL | 0.000e+00 | 3.092e+03 | 2.182e+06 | 5176 | 0 |
| res | PASS | 0.000e+00 | 1.043e-01 | 6.000e-01 | 138 | 0 |
| ortU | ALL | 1.000e-01 | 2.952e+10 | 4.289e+13 | 5176 | 0 |
| ortU | PASS | 1.000e-01 | 1.586e+00 | 5.000e+00 | 138 | 0 |
| ortV | ALL | 1.000e-01 | 2.952e+10 | 4.289e+13 | 5176 | 0 |
| ortV | PASS | 1.000e-01 | 1.581e+00 | 5.000e+00 | 138 | 0 |
| drift | ALL | 0.000e+00 | 9.706e-18 | 1.150e-15 | 5176 | 0 |
| drift | PASS | 0.000e+00 | 9.703e-17 | 1.080e-15 | 138 | 0 |
| t_eval | ALL | 4.315e-06 | 7.711e-03 | 4.310e-02 | 5176 | 0 |
| t_eval | PASS | 4.315e-06 | 2.780e-03 | 3.168e-02 | 138 | 0 |
| t_dbdsqr | ALL | 1.560e-06 | 4.438e-04 | 5.200e-03 | 5176 | 0 |
| t_dbdsqr | PASS | 1.780e-06 | 1.582e-04 | 1.314e-03 | 138 | 0 |
| ratio | ALL | 1.568e+00 | 2.451e+01 | 1.094e+03 | 5176 | 0 |
| ratio | PASS | 1.720e+00 | 2.285e+01 | 2.270e+02 | 138 | 0 |

### selfcontained  — rows=870  PASS=179

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 2.000e+00 | 1.085e+02 | 3.000e+02 | 870 | 0 |
| n | PASS | 4.000e+00 | 9.992e+01 | 2.670e+02 | 179 | 0 |
| res | ALL | 1.400e-02 | 2.878e-01 | 5.315e+00 | 870 | 0 |
| res | PASS | 2.500e-02 | 2.889e-01 | 2.333e+00 | 179 | 0 |
| ortU | ALL | 4.000e-03 | 7.568e+00 | 1.642e+02 | 870 | 0 |
| ortU | PASS | 4.000e-03 | 7.564e-01 | 4.245e+00 | 179 | 0 |
| ortV | ALL | 4.000e-03 | 7.570e+00 | 1.642e+02 | 870 | 0 |
| ortV | PASS | 4.000e-03 | 7.558e-01 | 4.193e+00 | 179 | 0 |
| drift | ALL | 4.390e-16 | 3.332e-15 | 5.110e-14 | 870 | 0 |
| drift | PASS | 6.390e-16 | 4.616e-15 | 5.110e-14 | 179 | 0 |
| t_eval | ALL | 2.962e-03 | 3.065e-02 | 1.470e-01 | 870 | 0 |
| t_eval | PASS | 2.962e-03 | 1.536e-02 | 1.036e-01 | 179 | 0 |
| t_dbdsqr | ALL | 4.910e-07 | 4.673e-04 | 2.700e-03 | 870 | 0 |
| t_dbdsqr | PASS | 1.574e-06 | 2.149e-04 | 2.286e-03 | 179 | 0 |
| ratio | ALL | 2.005e+01 | 2.055e+02 | 6.314e+03 | 870 | 0 |
| ratio | PASS | 2.005e+01 | 3.141e+02 | 2.034e+03 | 179 | 0 |

### dbdsqr  — rows=187  PASS=184

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 9.779e+01 | 2.670e+02 | 187 | 0 |
| n | PASS | 4.000e+00 | 9.931e+01 | 2.670e+02 | 184 | 0 |
| res | ALL | 0.000e+00 | 5.225e-01 | 1.000e+01 | 187 | 0 |
| res | PASS | 0.000e+00 | 3.859e-01 | 3.200e+00 | 184 | 0 |
| ortU | ALL | 0.000e+00 | 2.037e-01 | 9.000e-01 | 187 | 0 |
| ortU | PASS | 0.000e+00 | 1.957e-01 | 9.000e-01 | 184 | 0 |
| ortV | ALL | 0.000e+00 | 2.064e-01 | 1.000e+00 | 187 | 0 |
| ortV | PASS | 0.000e+00 | 2.027e-01 | 1.000e+00 | 184 | 0 |
| drift | ALL | 0.000e+00 | 2.805e-15 | 1.390e-14 | 187 | 0 |
| drift | PASS | 0.000e+00 | 2.845e-15 | 1.390e-14 | 184 | 0 |
| t_eval | ALL | 2.125e-06 | 4.605e-03 | 5.896e-02 | 187 | 0 |
| t_eval | PASS | 2.125e-06 | 4.680e-03 | 5.896e-02 | 184 | 0 |
| t_dbdsqr | ALL | 1.755e-06 | 2.365e-04 | 2.539e-03 | 187 | 0 |
| t_dbdsqr | PASS | 1.755e-06 | 2.403e-04 | 2.539e-03 | 184 | 0 |
| ratio | ALL | 7.005e-01 | 1.238e+01 | 4.762e+01 | 187 | 0 |
| ratio | PASS | 7.005e-01 | 1.254e+01 | 4.762e+01 | 184 | 0 |

## 379

### advisor  — rows=375  PASS=317

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.675e+02 | 4.290e+02 | 375 | 0 |
| n | PASS | 3.000e+00 | 1.569e+02 | 4.000e+02 | 317 | 0 |
| res | ALL | 0.000e+00 | 3.554e+11 | 3.688e+13 | 373 | 2 |
| res | PASS | 0.000e+00 | 7.319e-02 | 6.000e-01 | 317 | 0 |
| ortU | ALL | 0.000e+00 | 3.828e+11 | 3.999e+13 | 373 | 2 |
| ortU | PASS | 0.000e+00 | 9.861e-01 | 5.000e+00 | 317 | 0 |
| ortV | ALL | 0.000e+00 | 5.674e+11 | 4.504e+13 | 373 | 2 |
| ortV | PASS | 0.000e+00 | 9.801e-01 | 5.000e+00 | 317 | 0 |
| drift | ALL | 0.000e+00 | 3.430e-01 | 4.150e+01 | 373 | 2 |
| drift | PASS | 0.000e+00 | 3.942e-01 | 4.150e+01 | 317 | 0 |
| t_eval | ALL | 2.571e-06 | 3.762e-03 | 5.573e-02 | 373 | 2 |
| t_eval | PASS | 2.571e-06 | 3.221e-03 | 5.226e-02 | 317 | 0 |
| t_dbdsqr | ALL | 4.650e-07 | 5.898e-04 | 8.441e-03 | 373 | 2 |
| t_dbdsqr | PASS | 4.650e-07 | 4.927e-04 | 8.441e-03 | 317 | 0 |
| ratio | ALL | 5.063e-01 | 2.631e+01 | 7.705e+02 | 373 | 2 |
| ratio | PASS | 5.063e-01 | 3.003e+01 | 7.705e+02 | 317 | 0 |

### selfcontained  — rows=375  PASS=368

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.675e+02 | 4.290e+02 | 375 | 0 |
| n | PASS | 3.000e+00 | 1.666e+02 | 4.290e+02 | 368 | 0 |
| res | ALL | 0.000e+00 | 3.076e-01 | 2.006e+00 | 375 | 0 |
| res | PASS | 0.000e+00 | 3.100e-01 | 2.006e+00 | 368 | 0 |
| ortU | ALL | 0.000e+00 | 1.501e+11 | 2.252e+13 | 375 | 0 |
| ortU | PASS | 0.000e+00 | 4.716e-01 | 3.961e+00 | 368 | 0 |
| ortV | ALL | 0.000e+00 | 3.002e+10 | 1.126e+13 | 375 | 0 |
| ortV | PASS | 0.000e+00 | 4.737e-01 | 3.969e+00 | 368 | 0 |
| drift | ALL | 0.000e+00 | 1.198e+53 | 4.290e+55 | 358 | 17 |
| drift | PASS | 0.000e+00 | 3.024e+09 | 1.060e+12 | 352 | 16 |
| t_eval | ALL | 2.864e-03 | 5.342e-02 | 3.172e-01 | 375 | 0 |
| t_eval | PASS | 2.864e-03 | 5.231e-02 | 3.172e-01 | 368 | 0 |
| t_dbdsqr | ALL | 4.160e-07 | 5.326e-04 | 7.576e-03 | 358 | 17 |
| t_dbdsqr | PASS | 4.160e-07 | 5.277e-04 | 7.576e-03 | 352 | 16 |
| ratio | ALL | 1.111e+01 | 9.372e+02 | 1.035e+04 | 358 | 17 |
| ratio | PASS | 1.111e+01 | 9.495e+02 | 1.035e+04 | 352 | 16 |

### dbdsqr  — rows=375  PASS=375

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 3.000e+00 | 1.675e+02 | 4.290e+02 | 375 | 0 |
| n | PASS | 3.000e+00 | 1.675e+02 | 4.290e+02 | 375 | 0 |
| res | ALL | 0.000e+00 | 3.309e-01 | 4.500e+00 | 375 | 0 |
| res | PASS | 0.000e+00 | 3.309e-01 | 4.500e+00 | 375 | 0 |
| ortU | ALL | 0.000e+00 | 2.349e-01 | 1.000e+00 | 375 | 0 |
| ortU | PASS | 0.000e+00 | 2.349e-01 | 1.000e+00 | 375 | 0 |
| ortV | ALL | 0.000e+00 | 2.435e-01 | 1.400e+00 | 375 | 0 |
| ortV | PASS | 0.000e+00 | 2.435e-01 | 1.400e+00 | 375 | 0 |
| drift | ALL | 0.000e+00 | 2.672e-05 | 5.010e-03 | 375 | 0 |
| drift | PASS | 0.000e+00 | 2.672e-05 | 5.010e-03 | 375 | 0 |
| t_eval | ALL | 1.694e-06 | 4.532e-02 | 4.705e-01 | 375 | 0 |
| t_eval | PASS | 1.694e-06 | 4.532e-02 | 4.705e-01 | 375 | 0 |
| t_dbdsqr | ALL | 6.630e-07 | 5.937e-04 | 8.451e-03 | 375 | 0 |
| t_dbdsqr | PASS | 6.630e-07 | 5.937e-04 | 8.451e-03 | 375 | 0 |
| ratio | ALL | 8.658e-01 | 3.971e+01 | 2.471e+02 | 375 | 0 |
| ratio | PASS | 8.658e-01 | 3.971e+01 | 2.471e+02 | 375 | 0 |

## dense_to_bidiag

### advisor  — rows=224  PASS=179

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 8.000e+00 | 1.721e+02 | 1.005e+03 | 224 | 0 |
| n | PASS | 8.000e+00 | 1.617e+02 | 4.000e+02 | 179 | 0 |
| res | ALL | 0.000e+00 | 2.011e+12 | 4.504e+14 | 224 | 0 |
| res | PASS | 0.000e+00 | 8.659e-02 | 5.000e-01 | 179 | 0 |
| ortU | ALL | 0.000e+00 | 5.144e+00 | 5.524e+02 | 224 | 0 |
| ortU | PASS | 0.000e+00 | 8.564e-01 | 5.000e+00 | 179 | 0 |
| ortV | ALL | 0.000e+00 | 2.192e+12 | 4.504e+14 | 224 | 0 |
| ortV | PASS | 0.000e+00 | 8.615e-01 | 5.000e+00 | 179 | 0 |
| drift | ALL | 0.000e+00 | 3.406e-01 | 1.000e+00 | 224 | 0 |
| drift | PASS | 0.000e+00 | 4.039e-01 | 1.000e+00 | 179 | 0 |
| t_eval | ALL | 4.043e-06 | 5.858e-03 | 4.559e-01 | 224 | 0 |
| t_eval | PASS | 4.043e-06 | 2.840e-03 | 5.483e-02 | 179 | 0 |
| t_dbdsqr | ALL | 2.308e-06 | 8.478e-04 | 2.447e-02 | 224 | 0 |
| t_dbdsqr | PASS | 2.308e-06 | 6.143e-04 | 5.950e-03 | 179 | 0 |
| ratio | ALL | 2.139e-01 | 1.381e+01 | 3.422e+02 | 224 | 0 |
| ratio | PASS | 2.139e-01 | 1.530e+01 | 3.422e+02 | 179 | 0 |

### selfcontained  — rows=224  PASS=220

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 8.000e+00 | 1.721e+02 | 1.005e+03 | 224 | 0 |
| n | PASS | 8.000e+00 | 1.716e+02 | 1.005e+03 | 220 | 0 |
| res | ALL | 0.000e+00 | 3.323e-01 | 1.254e+00 | 224 | 0 |
| res | PASS | 0.000e+00 | 3.375e-01 | 1.254e+00 | 220 | 0 |
| ortU | ALL | 0.000e+00 | 7.623e-01 | 3.116e+01 | 224 | 0 |
| ortU | PASS | 0.000e+00 | 5.090e-01 | 3.564e+00 | 220 | 0 |
| ortV | ALL | 0.000e+00 | 7.657e-01 | 3.115e+01 | 224 | 0 |
| ortV | PASS | 0.000e+00 | 5.124e-01 | 3.586e+00 | 220 | 0 |
| drift | ALL | 0.000e+00 | 1.393e+29 | 2.870e+31 | 206 | 18 |
| drift | PASS | 0.000e+00 | 1.421e+29 | 2.870e+31 | 202 | 18 |
| t_eval | ALL | 2.995e-03 | 4.938e-02 | 1.040e+00 | 224 | 0 |
| t_eval | PASS | 2.995e-03 | 4.850e-02 | 1.040e+00 | 220 | 0 |
| t_dbdsqr | ALL | 1.390e-06 | 7.760e-04 | 2.191e-02 | 206 | 18 |
| t_dbdsqr | PASS | 1.390e-06 | 7.652e-04 | 2.191e-02 | 202 | 18 |
| ratio | ALL | 6.042e+00 | 3.868e+02 | 9.781e+03 | 206 | 18 |
| ratio | PASS | 6.042e+00 | 3.931e+02 | 9.781e+03 | 202 | 18 |

### dbdsqr  — rows=224  PASS=224

| column | set | min | mean | max | n_finite | n_nonfinite |
|---|---|---:|---:|---:|---:|---:|
| n | ALL | 8.000e+00 | 1.721e+02 | 1.005e+03 | 224 | 0 |
| n | PASS | 8.000e+00 | 1.721e+02 | 1.005e+03 | 224 | 0 |
| res | ALL | 0.000e+00 | 2.625e-01 | 4.000e+00 | 224 | 0 |
| res | PASS | 0.000e+00 | 2.625e-01 | 4.000e+00 | 224 | 0 |
| ortU | ALL | 0.000e+00 | 2.670e-01 | 1.500e+00 | 224 | 0 |
| ortU | PASS | 0.000e+00 | 2.670e-01 | 1.500e+00 | 224 | 0 |
| ortV | ALL | 0.000e+00 | 2.424e-01 | 1.300e+00 | 224 | 0 |
| ortV | PASS | 0.000e+00 | 2.424e-01 | 1.300e+00 | 224 | 0 |
| drift | ALL | 0.000e+00 | 4.018e-02 | 9.000e+00 | 224 | 0 |
| drift | PASS | 0.000e+00 | 4.018e-02 | 9.000e+00 | 224 | 0 |
| t_eval | ALL | 1.503e-06 | 5.507e-02 | 2.549e+00 | 224 | 0 |
| t_eval | PASS | 1.503e-06 | 5.507e-02 | 2.549e+00 | 224 | 0 |
| t_dbdsqr | ALL | 2.055e-06 | 8.488e-04 | 2.441e-02 | 224 | 0 |
| t_dbdsqr | PASS | 2.055e-06 | 8.488e-04 | 2.441e-02 | 224 | 0 |
| ratio | ALL | 5.512e-01 | 2.981e+01 | 2.464e+02 | 224 | 0 |
| ratio | PASS | 5.512e-01 | 2.981e+01 | 2.464e+02 | 224 | 0 |
