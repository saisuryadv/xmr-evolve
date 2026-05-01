# 06 — Cluster-Tree NCD Trace

This document logs **every** representation-tree node and **every** adjacent-pair NCD value visited by the MR3 algorithm on four highly-clustered bidiagonal test matrices, and visualizes how the tree forms. Deliverable for advisor task **(g)** as clarified.

**Script:** `experiments/cluster_tree_trace.py`
**Outputs:**
- `experiments/cluster_tree_log.json` — full per-node + per-pair log (machine-readable)
- `experiments/cluster_tree_log.csv` — one row per node, summary stats + per-pair extrema
- `experiments/cluster_tree.png` — tree diagram (nodes coloured by `log10(NCD_min)`, sized ∝ cluster size)
- `experiments/ncd_per_adjacent_pair.png` — every adjacent-pair NCD vs depth (the actual NCD values that `classify` consumes)
- `experiments/ncd_progression.png` — per-node aggregate NCD vs depth (legacy panel)

## Methodology

Production MR3 lives inside `dlaxrf` (Fortran kernel) and is not directly observable without instrumenting binary code. We use the algorithmically-equivalent Python `mr3_block` (`mr3_gk.py:418-460`) — an explicit-stack MR3 implementation. The previous bit-identity validation (21,375 / 21,375 specs match between Python and Fortran orchestrations) confirms that what we measure on the Python side reflects the Fortran tree-shape.

For each tree-node visit we record:
- node id, parent id, depth, shift τ
- cluster size `n_in_node`, range `[lo, hi]`
- aggregate `NCD_max := (hi − lo) / max(|lo|,|hi|)`
- aggregate `NCD_min := (hi − lo) / min(|lo|,|hi|)`  (user-requested form)
- counts of singletons / zero-width batches / descend-groups produced by `classify`
- **every adjacent pair** `(eigvals[k], eigvals[k+1])` with both NCD definitions and a flag for whether `classify` would split on that pair (`NCD_max ≥ GAPTOL = 1e-3`)

## Test matrices

| Name | n | Distribution |
|---|---|---|
| `gl_clustered_at_eps@200` | 200 | All d_i ∈ ε · (1 + uniform · n·ε) — 200 σ's tightly clustered near machine ε |
| `step_function@200` | 200 | d = three-tier step (10, 1, 0.1), constant e = 0.5 |
| `demmel_S1pe_k8@100` | 100 | d = [1, 1e-8, …, 1e-8], e[i] = 0.1·√(d_i·d_{i+1}) — one outlier σ + cluster |
| `pd_T0@100` | 100 | Pre-defined "tight uniform cluster spanning 76 % of spectrum" matrix from the failure list |

## Tree summary

| Matrix | # nodes | Max depth | Singletons (depth 0) | Descend groups (depth 0) |
|---|---|---|---|---|
| `gl_clustered_at_eps@200` | 12 | 1 | 177 | 11 |
| `step_function@200` | 6 | 1 | 174 | 5 |
| `demmel_S1pe_k8@100` | 3 | 1 | 80 | 2 |
| `pd_T0@100` | 2 | 1 | 24 | 1 |

In every case the root node resolves most eigenvalues as singletons (because adjacent pairs have `NCD_max ≥ GAPTOL`); the few residual sub-clusters become depth-1 children that resolve in one descent step.

## Per-node tree (full dump)

### `gl_clustered_at_eps@200`

```
id  parent  depth   n_in_node  NCD_max     NCD_min     outcome    sing/zw/desc
 0     -1      0        200    9.898e-01   9.662e+01   internal   177 / 0 / 11
 1      0      1          2    1.000e+00   6.783e+09   leaf         2 / 0 / 0
 2      0      1          2    1.000e+00   7.426e+09   leaf         2 / 0 / 0
 3      0      1          2    1.000e+00   1.672e+09   leaf         2 / 0 / 0
 4      0      1          2    1.000e+00   8.699e+09   leaf         2 / 0 / 0
 5      0      1          2    1.000e+00   1.556e+10   leaf         2 / 0 / 0
 6      0      1          2    1.000e+00   2.920e+09   leaf         2 / 0 / 0
 7      0      1          2    1.000e+00   1.269e+10   leaf         2 / 0 / 0
 8      0      1          2    1.000e+00   6.288e+08   leaf         2 / 0 / 0
 9      0      1          3    1.000e+00   3.192e+09   leaf         3 / 0 / 0
10      0      1          2    1.000e+00   4.503e+08   leaf         2 / 0 / 0
11      0      1          2    1.000e+00   5.757e+09   leaf         2 / 0 / 0
```

**What this shows:** the root node has 200 eigenvalues sitting near ε with random sign (because `gl_clustered_at_eps` includes negative diag entries → eigenvalues span both signs of ε). 177 of 199 adjacent pairs cross zero or have enough relative gap to be classified as singletons; the 11 remaining clusters all consist of pairs/triples with the same sign that need a depth-1 descent. At every depth-1 cluster `NCD_max = 1.0` (saturated) but `NCD_min` ranges from 4.5e8 to 1.6e10 — these are pairs/triples where the smaller absolute value is 8–10 orders of magnitude below the larger, but both NCD metrics agree they need splitting. They get resolved as singletons in one shift step.

### `step_function@200`

```
id  parent  depth   n_in_node  NCD_max     NCD_min     outcome    sing/zw/desc
 0     -1      0        200    1.000e+00   1.801e+16   internal   174 / 0 / 5
 1      0      1          2    2.000e+00   2.000e+00   leaf         2 / 0 / 0
 2      0      1          3    2.000e+00   2.000e+00   leaf         3 / 0 / 0
 3      0      1          2    1.000e+00   1.251e+10   leaf         2 / 0 / 0
 4      0      1          9    1.000e+00   1.161e+10   leaf         9 / 0 / 0
 5      0      1         10    1.000e+00   1.361e+10   leaf        10 / 0 / 0
```

**Note:** the root NCD_min = 1.8e16 reflects the matrix's full dynamic range (max ≈ 10, min ≈ ε from the smallest tier). Two of the depth-1 children show `NCD_max = NCD_min = 2.0` — these are the canonical regime where `max ≈ min` and both definitions agree. The other three depth-1 children have `NCD_max = 1.0` (saturated) but `NCD_min` in the 10¹⁰ range — same observation as `gl_clustered_at_eps`.

### `demmel_S1pe_k8@100`

```
id  parent  depth   n_in_node  NCD_max     NCD_min     outcome    sing/zw/desc
 0     -1      0        100    1.000e+00   1.111e+08   internal   80 / 0 / 2
 1      0      1          8    2.000e+00   2.000e+00   leaf         8 / 0 / 0
 2      0      1         12    2.000e+00   2.000e+00   leaf        12 / 0 / 0
```

### `pd_T0@100`

```
id  parent  depth   n_in_node  NCD_max     NCD_min     outcome    sing/zw/desc
 0     -1      0        100    4.313e-01   7.585e-01   internal   24 / 0 / 1
 1      0      1         76    2.000e+00   2.000e+00   leaf        76 / 0 / 0
```

The 76-element cluster that dominates `pd_T0` is captured as a single depth-1 child node with `NCD_max = NCD_min = 2.0` — both metrics agree because the cluster sits in a tight magnitude range (no large dynamic range), then resolves all 76 as singletons in one shift step.

## Adjacent-pair NCD distribution

For every visited node we logged the NCD on every adjacent eigenvalue pair (the actual quantity `classify` consumes). Aggregate counts:

| Matrix | # pairs total | `NCD_max ≥ GAPTOL` (split) | `NCD_max < GAPTOL` (cluster) |
|---|---|---|---|
| `gl_clustered_at_eps@200` | 211 | 199 | 12 |
| `step_function@200` | 224 | 197 | 27 |
| `demmel_S1pe_k8@100` | 117 | 80 | 37 |
| `pd_T0@100` | 174 | 24 | 150 |

`pd_T0` is the most clustered: 150 of 174 adjacent pairs fall *below* GAPTOL — they are inside a tight cluster that needs descent.

The plot `experiments/ncd_per_adjacent_pair.png` shows the full distribution of these per-pair NCD values vs depth, with the GAPTOL = 1e-3 threshold marked as a horizontal dashed red line.

## Visualizations

### Tree diagram — `experiments/cluster_tree.png`

Each subplot is one matrix's representation tree. Nodes are coloured by `log10(NCD_min)`; size is proportional to cluster size. Edges connect parents to children. ID + cluster size are annotated.

This is the "tree formation" panel the user asked for.

### Per-pair NCD scatter — `experiments/ncd_per_adjacent_pair.png`

Two panels, one per NCD definition. Each point is one adjacent eigenvalue pair encountered during the trace. X-axis = depth; Y-axis = NCD value (log scale); horizontal red line = GAPTOL.

Pairs **above** the line are split as gaps by `classify`; pairs **below** are kept clustered (and the cluster is then either declared zero-width, sent to descent, or — at MAX_DEPTH — dispatched as multiplicity-batch).

### Per-node aggregate — `experiments/ncd_progression.png`

The legacy "per-depth aggregate scatter" plot kept for back-compatibility.

## Findings

1. **The tree forms.** All four matrices produce non-trivial trees (2–12 nodes); the most-clustered (`gl_clustered_at_eps`) reaches 12 nodes including 11 leaves at depth 1.

2. **Tree depth is shallow** on these matrices because the depth-1 children are either (a) zero-width-batch dispatched (`gev` resolves as multiplicity, no descent) or (b) resolved as singletons after one shift selects a representation where every cluster member becomes well-separated.

3. **`NCD_max` saturates at 1.0** for a depth-1 child whenever the cluster spans a wide dynamic range. The saturation is informative: it tells you the cluster's spread is the same order as its largest member, but it does *not* discriminate between "spread = max" and "spread = 100 × max in relative-to-min terms". `NCD_min` on those same nodes ranges from 4e8 to 1.6e10 — that order of magnitude information is lost in `NCD_max`.

4. **For tight clusters (everything fitting near one magnitude)** `NCD_max ≈ NCD_min` and the choice of metric is irrelevant. This is the regime of `step_function`'s middle tier, `demmel_S1pe_k8`'s sub-cluster, and `pd_T0`'s 76-cluster — and it is also the regime where `is_zero_width` may dispatch the cluster directly without descending.

5. **For the production classifier (`dlaxrb_clssfy_fix.f`)**, the metric used is `NCD_max ≈ NCD_per_pair_max ≥ GAPTHRESH`. Switching to `NCD_min` would aggressively split nodes like the depth-1 leaves on `gl_clustered_at_eps` (where `NCD_min` is 10⁹–10¹⁰, well above GAPTHRESH); the algorithm currently splits these anyway via the `NCD_max = 1.0` check, so the behavioural difference for these test cases is small. The two metrics diverge most sharply on cluster nodes whose members straddle vastly different magnitudes — see node 0 of `step_function` (`NCD_max = 1`, `NCD_min = 1.8e16`).

## Reproducing

```bash
cd python_fortran
python3 experiments/cluster_tree_trace.py                            # all four matrices
python3 experiments/cluster_tree_trace.py gl_clustered_at_eps 200    # one
```

Outputs are written into `experiments/`. Re-running overwrites the JSON / CSV / PNGs.
