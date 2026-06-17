# DBDSVDMR3 failures across all 4 suites — ranked by severity

Total failures parsed from logs: **176** (pract 54, 379 40, dense_to_bidiag 32, synth top-50).

> Note: synth ran 18438 matrices with 4478 failures total, but the original orchestrator capped the "Failing" block at top-50 by max(ortU, ortV). The ranker has been updated to dump every failure; the full synth list will appear once the sweep is re-run. The 50 captured here are already the worst cases.

Severity = `max(res, ortU, ortV)` in units of `n·eps`. Thresholds: res ≤ 7, ortU/ortV ≤ 5.

Buckets (decreasing severity):
- **catastrophic** (TIMEOUT / parse-fail / INF): 18
- **huge** (≥ 1e6 n·ε): 28
- **big** (1e3 .. 1e6 n·ε): 43
- **mild** (< 1e3 n·ε): 87

Each row carries a one-liner to reproduce it via the orchestrator.

```
cd python_fortran
bash BidiagonalSVD_TGK/build.sh   # once
```

## Catastrophic (TIMEOUT / parse-fail / INF) (18)

| # | suite | matrix | n | res | ortU | ortV | note | repro |
|---|---|---|---|---|---|---|---|---|
| 1 | pract | `T_Alemdar_1` | 6245 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_Alemdar_1'` |
| 2 | pract | `T_TSC_OPF_300` | 9774 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_TSC_OPF_300'` |
| 3 | pract | `T_bcsstkm10_4` | 4344 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm10_4'` |
| 4 | pract | `T_bcsstkm11_3` | 4419 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm11_3'` |
| 5 | pract | `T_bcsstkm11_4` | 5892 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm11_4'` |
| 6 | pract | `T_bcsstkm12_3` | 4419 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm12_3'` |
| 7 | pract | `T_bcsstkm13_3` | 6009 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm13_3'` |
| 8 | pract | `T_c-40` | 9941 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_c-40'` |
| 9 | pract | `T_nasa1824_3` | 5472 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa1824_3'` |
| 10 | pract | `T_nasa4704_1` | 4704 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa4704_1'` |
| 11 | pract | `T_sts4098_1` | 4098 | inf | inf | inf | TIMEOUT | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_sts4098_1'` |
| 12 | synth | `ev2_ec1_n35_s0_gm` | 105 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev2_ec1_n35_s0_gm'` |
| 13 | synth | `ev2_ec1_n75_s2_gm` | 225 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev2_ec1_n75_s2_gm'` |
| 14 | synth | `ev2_ec4_n87_s1_gm` | 261 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev2_ec4_n87_s1_gm'` |
| 15 | synth | `ev2_ec4_n92_s3_gm` | 276 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev2_ec4_n92_s3_gm'` |
| 16 | synth | `ev5_ec4_n64_s2_gm` | 192 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev5_ec4_n64_s2_gm'` |
| 17 | synth | `ev5_ec4_n88_s3_gm` | 264 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev5_ec4_n88_s3_gm'` |
| 18 | synth | `ev5_ec4_n95_s2_gm` | 285 | inf | inf | inf | parse-fail | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'ev5_ec4_n95_s2_gm'` |

## Huge (>= 1e6 n.eps) (28)

| # | suite | matrix | n | res | ortU | ortV | note | repro |
|---|---|---|---|---|---|---|---|---|
| 1 | dense_to_bidiag | `dense:dense_three_clusters@10` | 10 | 0.288 | 0.5 | 4.11e+14 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_three_clusters@10'` |
| 2 | 379 | `adv:demmel_S2pe@100` | 100 | 0.649 | 2.6 | 4.5e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:demmel_S2pe@100'` |
| 3 | dense_to_bidiag | `dense:dense_repeated_sv@100` | 100 | 0.829 | 2.49 | 2.99e+13 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_repeated_sv@100'` |
| 4 | 379 | `adv:gl_wilkp@200` | 200 | 1.2 | 1.43 | 2.25e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_wilkp@200'` |
| 5 | 379 | `adv:wilkinson_exact@200` | 200 | 0.957 | 1.36 | 2.25e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:wilkinson_exact@200'` |
| 6 | 379 | `adv:demmel_S2pe@200` | 200 | 0.399 | 2.36 | 2.25e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:demmel_S2pe@200'` |
| 7 | dense_to_bidiag | `dense:dense_three_clusters@200` | 200 | 0.101 | 1.75 | 1.89e+13 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_three_clusters@200'` |
| 8 | 379 | `adv:gl_bwilkp@100` | 100 | 2.91 | 1.46e+13 | 1.46e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_bwilkp@100'` |
| 9 | 379 | `adv:gl_wilkm@100` | 100 | 1.04 | 1.46e+13 | 6.26 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_wilkm@100'` |
| 10 | 379 | `adv:demmel_S2pe@400` | 400 | 0.243 | 1.96 | 1.12e+13 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:demmel_S2pe@400'` |
| 11 | dense_to_bidiag | `dense:dense_repeated_sv@200` | 200 | 0.234 | 0.856 | 9.39e+12 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_repeated_sv@200'` |
| 12 | 379 | `adv:gl_bwilkp@200` | 200 | 1.24 | 6.15e+12 | 6.15e+12 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_bwilkp@200'` |
| 13 | 379 | `adv:gl_wilkm@200` | 200 | 0.43 | 6.15e+12 | 1.7e+12 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_wilkm@200'` |
| 14 | dense_to_bidiag | `dense:dense_repeated_sv@400` | 400 | 0.406 | 3.27 | 4.94e+12 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_repeated_sv@400'` |
| 15 | 379 | `adv:gl_bwilkp@400` | 400 | 0.525 | 2.59e+12 | 2.59e+12 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_bwilkp@400'` |
| 16 | 379 | `adv:gl_wilkm@400` | 400 | 0.395 | 2.59e+12 | 6.08e+11 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_wilkm@400'` |
| 17 | pract | `T_W21_g_1e+12` | 2100 | 1.98e+06 | 8.3e+06 | 8.3e+06 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e+12'` |
| 18 | synth | `tri2_n23_gm` | 69 | 3.34e+05 | 5.49e+06 | 5.49e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n23_gm'` |
| 19 | synth | `tri2_n65_gm` | 195 | 1.49e+05 | 4.59e+06 | 4.59e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n65_gm'` |
| 20 | synth | `tri4_n3_gm` | 9 | 1.5e+06 | 2.64e+06 | 2.02e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n3_gm'` |
| 21 | synth | `tri2_n37_gm` | 111 | 1.22e+05 | 2.44e+06 | 2.44e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n37_gm'` |
| 22 | synth | `tri3_n11_gm` | 33 | 4.48e+05 | 1.7e+06 | 1.7e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri3_n11_gm'` |
| 23 | synth | `tri2_n27_gm` | 81 | 9.4e+04 | 1.68e+06 | 1.68e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n27_gm'` |
| 24 | synth | `tri2_n17_gm` | 51 | 8.37e+04 | 1.39e+06 | 1.39e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n17_gm'` |
| 25 | synth | `tri2_n13_gm` | 39 | 2.02e+05 | 1.34e+06 | 1.34e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n13_gm'` |
| 26 | synth | `tri2_n9_gm` | 27 | 2.56e+05 | 1.3e+06 | 1.3e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n9_gm'` |
| 27 | synth | `tri2_n25_gm` | 75 | 6.37e+04 | 9.99e+05 | 1.11e+06 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n25_gm'` |
| 28 | synth | `tri3_n7_gm` | 21 | 3.7e+05 | 1.02e+06 | 9.81e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri3_n7_gm'` |

## Big (1e3 .. 1e6 n.eps) (43)

| # | suite | matrix | n | res | ortU | ortV | note | repro |
|---|---|---|---|---|---|---|---|---|
| 1 | synth | `tri4_n57_gm` | 171 | 5.48e+04 | 9.57e+05 | 9.57e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n57_gm'` |
| 2 | synth | `tri4_n9_gm` | 27 | 2.95e+05 | 8.75e+05 | 8.33e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n9_gm'` |
| 3 | synth | `tri2_n73_gm` | 219 | 2.62e+04 | 8.62e+05 | 8.56e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n73_gm'` |
| 4 | synth | `tri2_n3_gm` | 9 | 3.7e+05 | 8.5e+05 | 6.51e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n3_gm'` |
| 5 | synth | `tri3_n9_gm` | 27 | 2.97e+05 | 8.35e+05 | 8.35e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri3_n9_gm'` |
| 6 | synth | `tri2_n55_gm` | 165 | 4.74e+04 | 8.24e+05 | 8.04e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n55_gm'` |
| 7 | synth | `tri2_n35_gm` | 105 | 4.65e+04 | 7.72e+05 | 7.74e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n35_gm'` |
| 8 | synth | `tri2_n7_gm` | 21 | 2.59e+05 | 7.7e+05 | 7.7e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n7_gm'` |
| 9 | synth | `tri4_n55_gm` | 165 | 7.73e+04 | 7.19e+05 | 7.19e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n55_gm'` |
| 10 | synth | `tri2_n61_gm` | 183 | 2.3e+04 | 7.14e+05 | 7.14e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n61_gm'` |
| 11 | synth | `tri4_n19_gm` | 57 | 1.74e+05 | 7.04e+05 | 7.04e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n19_gm'` |
| 12 | synth | `tri4_n43_gm` | 129 | 9.77e+04 | 6.88e+05 | 6.88e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n43_gm'` |
| 13 | synth | `tri2_n29_gm` | 87 | 2.7e+04 | 6.15e+05 | 6.84e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n29_gm'` |
| 14 | synth | `tri2_n93_gm` | 279 | 2.21e+04 | 6.51e+05 | 6.51e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n93_gm'` |
| 15 | synth | `tri4_n21_gm` | 63 | 9.58e+04 | 6.46e+05 | 6.46e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n21_gm'` |
| 16 | synth | `tri2_n33_gm` | 99 | 4.76e+04 | 6.11e+05 | 6.11e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n33_gm'` |
| 17 | synth | `tri2_n21_gm` | 63 | 7.64e+04 | 5.93e+05 | 5.93e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n21_gm'` |
| 18 | synth | `tri3_n5_gm` | 15 | 2.42e+05 | 5.65e+05 | 5.65e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri3_n5_gm'` |
| 19 | synth | `tri2_n57_gm` | 171 | 2.69e+04 | 5.36e+05 | 5.3e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n57_gm'` |
| 20 | synth | `tri2_n47_gm` | 141 | 2.74e+04 | 5.05e+05 | 5.05e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n47_gm'` |
| 21 | synth | `tri4_n35_gm` | 105 | 6.11e+04 | 4.76e+05 | 4.63e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n35_gm'` |
| 22 | synth | `tri4_n41_gm` | 123 | 9.38e+04 | 4.67e+05 | 4.57e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n41_gm'` |
| 23 | synth | `tri2_n51_gm` | 153 | 1.87e+04 | 4.64e+05 | 4.64e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n51_gm'` |
| 24 | synth | `tri4_n31_gm` | 93 | 8.02e+04 | 4.34e+05 | 4.34e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n31_gm'` |
| 25 | synth | `tri2_n19_gm` | 57 | 6.64e+04 | 4.32e+05 | 4.32e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n19_gm'` |
| 26 | synth | `tri2_n95_gm` | 285 | 1.46e+04 | 4.14e+05 | 4.1e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n95_gm'` |
| 27 | synth | `tri2_n53_gm` | 159 | 1.94e+04 | 4.08e+05 | 4.04e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n53_gm'` |
| 28 | synth | `tri2_n31_gm` | 93 | 5.68e+04 | 4.02e+05 | 4.02e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n31_gm'` |
| 29 | synth | `tri2_n45_gm` | 135 | 3.09e+04 | 3.94e+05 | 3.9e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n45_gm'` |
| 30 | synth | `tri3_n13_gm` | 39 | 9.88e+04 | 3.89e+05 | 3.89e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri3_n13_gm'` |
| 31 | synth | `tri4_n15_gm` | 45 | 7.92e+04 | 3.85e+05 | 3.85e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri4_n15_gm'` |
| 32 | synth | `tri2_n5_gm` | 15 | 9.28e+04 | 3.76e+05 | 3.63e+05 |  | `python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper --only 'tri2_n5_gm'` |
| 33 | 379 | `adv:three_clusters@400` | 400 | 0.9 | 1.98e+05 | 1.98e+05 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:three_clusters@400'` |
| 34 | 379 | `adv:spike@100` | 100 | 5 | 1.46e+05 | 1.47e+05 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:spike@100'` |
| 35 | 379 | `adv:spike@400` | 400 | 1.3 | 7.53e+04 | 7.52e+04 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:spike@400'` |
| 36 | pract | `T_W21_g_1e+02` | 2100 | 885 | 6.76e+04 | 6.76e+04 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e+02'` |
| 37 | 379 | `adv:spike@200` | 200 | 4.9 | 6.03e+04 | 6.05e+04 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:spike@200'` |
| 38 | 379 | `adv:saw_tooth@400` | 400 | 0 | 2.32e+04 | 2.32e+04 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:saw_tooth@400'` |
| 39 | pract | `T_W21_g_1e+06` | 2100 | 2.09e+03 | 8.77e+03 | 8.77e+03 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e+06'` |
| 40 | pract | `T_W21_g_1e+14` | 2100 | 3.55e+03 | 4.02e+03 | 4.02e+03 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e+14'` |
| 41 | 379 | `adv:demmel_S2ps@400` | 400 | 0.2 | 2.73e+03 | 2.73e+03 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:demmel_S2ps@400'` |
| 42 | 379 | `stcoll:B_Kimura_429` | 429 | 252 | 1.02e+03 | 1.02e+03 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'stcoll:B_Kimura_429'` |
| 43 | pract | `STColl_B_Kimura_429` | 429 | 252 | 1.02e+03 | 1.02e+03 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'STColl_B_Kimura_429'` |

## Mild (< 1e3 n.eps) (87)

| # | suite | matrix | n | res | ortU | ortV | note | repro |
|---|---|---|---|---|---|---|---|---|
| 1 | 379 | `adv:saw_tooth@200` | 200 | 0 | 522 | 522 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:saw_tooth@200'` |
| 2 | pract | `T_W21_g_1e+04` | 2100 | 10.6 | 519 | 519 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e+04'` |
| 3 | pract | `T_0007a` | 7 | 84.5 | 288 | 288 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_0007a'` |
| 4 | pract | `T_SkewW21gve+3` | 2100 | 251 | 161 | 161 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_SkewW21gve+3'` |
| 5 | dense_to_bidiag | `paper:glued_wilk_5x21_sqrteps` | 105 | 49.1 | 240 | 240 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:glued_wilk_5x21_sqrteps'` |
| 6 | pract | `Lipshitz_4` | 1088 | 12.8 | 206 | 206 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'Lipshitz_4'` |
| 7 | pract | `T_bcsstkm10_2` | 2172 | 21.1 | 95.1 | 95.1 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm10_2'` |
| 8 | pract | `T_nasa1824_1` | 1824 | 8.7 | 82.4 | 82.4 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa1824_1'` |
| 9 | pract | `T_W21_g_1e-14` | 2100 | 2.4 | 53.3 | 53.3 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_W21_g_1e-14'` |
| 10 | pract | `T_0016_smalleig` | 16 | 52.4 | 51.7 | 51.5 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_0016_smalleig'` |
| 11 | pract | `T_nasa2910` | 2910 | 4.9 | 46.3 | 46.3 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa2910'` |
| 12 | dense_to_bidiag | `paper:glued_wilk_3x21_sqrteps` | 63 | 5.4 | 34 | 34 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:glued_wilk_3x21_sqrteps'` |
| 13 | pract | `T_bcsstkm10_3` | 3258 | 5.8 | 30 | 30 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm10_3'` |
| 14 | pract | `T_0010_stexrfailure_TGK` | 20 | 9.6 | 28.2 | 28.2 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_0010_stexrfailure_TGK'` |
| 15 | pract | `T_matlab_ud_2000` | 2000 | 2.3 | 28.2 | 28.2 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_ud_2000'` |
| 16 | pract | `T_matlab_ud_2250` | 2250 | 1.1 | 21.9 | 21.9 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_ud_2250'` |
| 17 | dense_to_bidiag | `paper:all_equal_nontrivial_glued_small@100` | 99 | 1.3 | 21.2 | 21.2 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:all_equal_nontrivial_glued_small@100'` |
| 18 | pract | `T_matlab_ud_1750` | 1750 | 0.7 | 18.2 | 18.2 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_ud_1750'` |
| 19 | dense_to_bidiag | `dense:dense_random@400` | 400 | 0.8 | 17 | 17 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_random@400'` |
| 20 | pract | `T_matlab_ud_1250` | 1250 | 1.6 | 16.3 | 16.3 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_ud_1250'` |
| 21 | 379 | `adv:gl_abcon3@400` | 400 | 0 | 15.3 | 15.3 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_abcon3@400'` |
| 22 | pract | `T_bcsstkm12_1` | 1473 | 2.1 | 14.1 | 14.1 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm12_1'` |
| 23 | 379 | `adv:wilkinson_like@10` | 10 | 3.9 | 13.9 | 14 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:wilkinson_like@10'` |
| 24 | pract | `T_matlab_nd_1250` | 1250 | 0.8 | 13.9 | 13.9 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_nd_1250'` |
| 25 | dense_to_bidiag | `dense:dense_moler_like@200` | 200 | 2.8 | 13.6 | 13.6 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_moler_like@200'` |
| 26 | 379 | `adv:two_clusters@10` | 10 | 0.2 | 13.5 | 13.5 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:two_clusters@10'` |
| 27 | pract | `T_bcsstkm01_3` | 144 | 2.7 | 12.8 | 12.8 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_bcsstkm01_3'` |
| 28 | dense_to_bidiag | `dense:dense_wilkinson_sv@200` | 199 | 0.1 | 12.7 | 12.7 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_wilkinson_sv@200'` |
| 29 | dense_to_bidiag | `paper:glued_wilk_3x21_1e14` | 63 | 1.7 | 12.2 | 12.2 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:glued_wilk_3x21_1e14'` |
| 30 | dense_to_bidiag | `paper:random_zero_e_10pct@10` | 10 | 3.6 | 12 | 12 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:random_zero_e_10pct@10'` |
| 31 | dense_to_bidiag | `paper:glued_wilk201_5x_sqrteps` | 1005 | 3.5 | 11.1 | 11.1 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'paper:glued_wilk201_5x_sqrteps'` |
| 32 | pract | `T_matlab_nd_0750` | 750 | 0.6 | 11.1 | 11.1 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_nd_0750'` |
| 33 | pract | `T_nasa1824` | 1824 | 0.8 | 11.1 | 11.1 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa1824'` |
| 34 | 379 | `adv:step_function@200` | 200 | 0.8 | 10.8 | 10.8 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:step_function@200'` |
| 35 | dense_to_bidiag | `dense:dense_moler_like@400` | 400 | 1.3 | 10.4 | 10.4 |  | `python3 eval_dbdsvdmr3.py --suite dense_to_bidiag --only 'dense:dense_moler_like@400'` |
| 36 | pract | `T_nasa2146` | 2146 | 0.4 | 10.2 | 10.2 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_nasa2146'` |
| 37 | 379 | `adv:gl_wilkw@100` | 100 | 3.8 | 9.8 | 9.8 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:gl_wilkw@100'` |
| 38 | pract | `T_matlab_nd_1750` | 1750 | 0.7 | 9.8 | 9.8 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_nd_1750'` |
| 39 | 379 | `adv:random_dense_clusters@400` | 400 | 4.4 | 9.6 | 9.6 |  | `python3 eval_dbdsvdmr3.py --suite 379 --only 'adv:random_dense_clusters@400'` |
| 40 | pract | `T_matlab_ud_1000` | 1000 | 0.8 | 9.4 | 9.4 |  | `python3 eval_dbdsvdmr3.py --suite pract --only 'T_matlab_ud_1000'` |

_…47 more rows omitted; full list derivable from `docs/*_dbdsvdmr3_2026-06-12.log`._

## Reproducing the entire sweep

```
cd python_fortran
bash BidiagonalSVD_TGK/build.sh
python3 eval_dbdsvdmr3.py --suite pract           2>&1 | tee docs/pract_dbdsvdmr3.log
python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper 2>&1 | tee docs/synth_dbdsvdmr3.log
python3 eval_dbdsvdmr3.py --suite 379             2>&1 | tee docs/lapack_dbdsvr_379_dbdsvdmr3.log
python3 eval_dbdsvdmr3.py --suite dense_to_bidiag 2>&1 | tee docs/dense_to_bidiag_dbdsvdmr3.log
```
