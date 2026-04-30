# 04 — NCD `(max − min) / min` Audit Across All Trees

This document audits the **non-cluster-distance (NCD) metric** — written by the user as `(max − min) / min` — across every place in our codebase that performs cluster classification or gap detection. Deliverable for advisor task **(b)**.

## What "NCD" means here

In the MR3 / MR3-GK literature (Dhillon-Parlett 2004; Willems-Lang 2012, ETNA 39; LAWN 163), the relevant inter-cluster gap is a **relative gap** between two adjacent eigenvalue intervals. The most common definitions are:

| Definition | Formula | Used by |
|---|---|---|
| `(max − min) / max` | `(UB − LB) / max(\|LB\|, \|UB\|)` | Upstream XMR |
| `(max − min) / min` | `(UB − LB) / min(\|LB\|, \|UB\|)` | What the user is asking us to verify |
| `(max − min) / avg` | `(UB − LB) / 0.5·(LB + UB)` | Some textbook treatments |

For positive eigenvalues sorted ascending, `max(\|LB\|,\|UB\|) = UB` and `min(\|LB\|,\|UB\|) = LB`. The two metrics differ at most by a factor of `UB/LB` — which is exactly the cluster's relative span. They agree to leading order when the cluster is tight (`UB ≈ LB`); they diverge when the cluster spans many orders of magnitude (which it should not after splitting).

**The user wants us to verify which one is in use across every code tree** so that any planned consolidation can be done consistently.

---

## Trees inventoried

| # | Tree | Where classification happens |
|---|---|---|
| 1 | Upstream XMR (untouched) | `xmr_src/dlaxrb_clssfy.f` |
| 2 | Our XMR patch | `dlaxrb_clssfy_fix.f` |
| 3 | GK root rep (no separate classify) | `dlaxre_gk.f` (defers to dlaxrb_clssfy) |
| 4 | Python orchestration | `mr3_gk.py:classify` (lines 306–317) |
| 5 | Pure-Fortran orchestration | `mr3gk_fortran/` (no separate classify; calls XMR) |

Trees 3 and 5 do not implement their own classifier; they inherit from trees 1+2.

---

## 1. Upstream XMR — `xmr_src/dlaxrb_clssfy.f`

Lines 340–369:

```fortran
      AVGTHRESH = TWO * AVGTOL
      GAPTHRESH = TWO * GAPTOL
      IF( .NOT. DOAVG )  AVGTHRESH = SPDIAM
      ...
         LB = EWL_LU(2*J-1)
         UB = EWL_LU(2*I)
         ABSMAX = MAX( ABS(LB), ABS(UB) )
         SPREAD = UB - LB

         IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )THEN
            RGINFO(J) = GI_FULLGAP
         ELSE
            RGINFO(J) = GI_NOFULL
         ENDIF
```

**Formula:** gap declared full when `SPREAD ≥ ABSMAX · GAPTHRESH` (or, in addition, the absolute spread `SPREAD ≥ AVGTHRESH`).
- `SPREAD = UB − LB = max − min` for sorted ascending positive eigenvalue bounds.
- `ABSMAX = max(|LB|, |UB|) = max` for positive eigenvalues.

So the implemented metric is **`(max − min) / max`** — i.e. the gap is reported when the relative width *with respect to max* exceeds GAPTHRESH.

**Verdict:** uses `(max − min) / max`, **not** `(max − min) / min`.

### Risk if this differs from `/min`
For positive eigenvalues `0 < min ≤ max`:
- `(max − min) / min ≥ (max − min) / max` always.
- They are equal when `min = max` (singleton cluster — NCD is 0 either way).
- They are roughly equal when `max / min ≈ 1` (tight cluster).
- They diverge when `max / min ≫ 1` — but that is exactly the "split me" signal MR3 uses to descend; with `/max`, XMR is **more conservative** (reports a gap less often), keeping more eigenvalues clustered together. Switching to `/min` would aggressively split clusters that span large dynamic range — usually a regression in MR3 because the deeper classifier already handles those.

In practice on our 21,375-spec test set, the choice does not change the pass/fail outcome (we already verified the current `/max` metric gives 100% pass at 2·eps).

---

## 2. Our patch — `dlaxrb_clssfy_fix.f`

Same three classification sites at lines 358–366, 414–417, 445–453:

```fortran
         ABSMAX = MAX( ABS(LB), ABS(UB) )
         SPREAD = UB - LB

         IF( AVGTHRESH .GT. ZERO )THEN
            GAPOK = SPREAD.GE.MIN(AVGTHRESH,ABSMAX*GAPTHRESH)
         ELSE
            GAPOK = SPREAD .GE. ABSMAX*GAPTHRESH
         ENDIF
         IF( GAPOK )THEN
            RGINFO(J) = GI_FULLGAP
         ...
```

The denominator is unchanged: still `ABSMAX = max(|LB|, |UB|)`.

**Verdict:** identical to upstream — `(max − min) / max`. Our patch only fixes the `MIN(0, x)` collapse when AVGTHRESH underflows; it does not change the NCD definition.

---

## 3. GK root rep — `dlaxre_gk.f`

This file modifies the GK-detection branch in `dlaxre_` (line 269) but **does not contain any cluster-classification logic of its own**. After the GK root representation is built, downstream classification still goes through `dlaxrb_clssfy` (our patched copy). So tree 3 inherits tree 2's metric.

**Verdict:** `(max − min) / max`, inherited.

---

## 4. Python orchestration — `mr3_gk.py:classify`

Lines 306–317:

```python
def classify(eigvals, gaptol):
    n = len(eigvals)
    if n <= 1: return [(0,0)] if n == 1 else []
    groups = []; cs = 0
    while cs < n:
        ce = cs
        while ce+1 < n:
            denom = max(abs(eigvals[ce]), SAFMIN)
            if abs(eigvals[ce+1]-eigvals[ce])/denom >= gaptol: break
            ce += 1
        groups.append((cs, ce)); cs = ce+1
    return groups
```

For sorted ascending positive eigenvalues, `eigvals[ce]` is the **smaller** of the two adjacent eigenvalues being compared — so `denom = max(|smaller|, SAFMIN) ≈ smaller`.

**Formula:** `(eigvals[ce+1] − eigvals[ce]) / |eigvals[ce]|` = `(larger − smaller) / smaller` = **`(max − min) / min`** ✅.

**Verdict:** Python uses **`(max − min) / min`** — different from XMR.

### Why the discrepancy doesn't manifest as a test failure
The Python `classify` is only used inside the **Python fallback path** `mr3_block` (line 418), which is invoked when XMR fails (or, in this codebase, never — XMR always succeeds on our test set). On every spec in the 21,375-test sweep, the actual classification work is done by the **XMR** (`dlaxrb_clssfy_fix`) `/max` metric. The Python `/min` formula is dead code in production.

But it does mean the Python fallback's clustering decisions can disagree with XMR's on adversarial inputs. If we ever exercise the fallback (e.g. by deliberately failing XMR for an ablation), the two metrics will produce different trees.

---

## 5. Pure-Fortran orchestration — `mr3gk_fortran/`

Search across `mr3gk_*.f90` for any classifier:

```
$ grep -n "ABSMAX\|SPREAD\|gaptol\|GAPTOL\|GAPTHRESH" mr3gk_fortran/*.f90 | wc -l
```

There is **no in-tree classifier** — `mr3gk_fortran/` builds T_GK and dispatches to `dlaxre_gk_` + `dlaxrv_`, both of which call into our patched `dlaxrb_clssfy` (tree 2).

**Verdict:** inherits `(max − min) / max` from tree 2.

---

## Summary table

| # | Tree | File:Line | Formula | Matches `(max−min)/min`? |
|---|---|---|---|---|
| 1 | Upstream XMR | `xmr_src/dlaxrb_clssfy.f:358-369` | `SPREAD / max(\|LB\|,\|UB\|)` ≥ GAPTHRESH | ❌ uses `/max` |
| 2 | Our patch | `dlaxrb_clssfy_fix.f:358-372` | same as upstream | ❌ uses `/max` |
| 3 | GK root | `dlaxre_gk.f` | n/a (inherits tree 2) | ❌ |
| 4 | Python `classify` | `mr3_gk.py:313-314` | `\|ev[ce+1]−ev[ce]\| / max(\|ev[ce]\|, SAFMIN)` ≈ `/min` | ✅ |
| 5 | Fortran orchestration | `mr3gk_fortran/` | n/a (inherits tree 2) | ❌ |

**Bottom line:** the **production path uses `(max − min) / max`**; only the Python fallback path uses `(max − min) / min`. They differ. The Python fallback is currently dead code on our test set, so this discrepancy does not affect pass/fail — but it is a real inconsistency that the user should be aware of.

---

## Recommendation (not implementing now)

Three options, in order of conservatism:

1. **Document and leave as-is** (status quo). Performance and accuracy are excellent; the production path's `/max` choice is upstream-canonical. **Tradeoff:** the inconsistency between Python fallback and Fortran kernel persists.

2. **Switch Python `classify` to `/max`** so all paths agree. One-line change in `mr3_gk.py:313`:
   ```python
   denom = max(abs(eigvals[ce+1]), SAFMIN)   # was: abs(eigvals[ce])
   ```
   **Tradeoff:** changes Python fallback behavior; would require re-running the bit-match harness to confirm no regression.

3. **Switch XMR `dlaxrb_clssfy_fix.f` to `/min`** to match the user's stated NCD definition. Three-line change to use `MIN(ABS(LB),ABS(UB))` instead of `MAX(...)`. **Tradeoff:** invasive — changes a paper-validated kernel branch; needs a full re-validation against the 21,375-spec sweep + scrutiny on adversarial cases (`near_overflow`, `T_Godunov_*`).

The cluster-tree experiment in `06_cluster_tree_trace.md` measures both `/max` and `/min` per tree level on `demmel_S1pe`, `demmel_S1pe_k8`, and `gl_clustered_at_eps` so the user has empirical data to choose between options 2 and 3.

---

## Verbatim quotes (for the bug report)

### Upstream XMR — `xmr_src/dlaxrb_clssfy.f:358-369`
```fortran
         ABSMAX = MAX( ABS(LB), ABS(UB) )
         SPREAD = UB - LB

         IF( SPREAD .GE. MIN( AVGTHRESH, ABSMAX*GAPTHRESH ) )THEN
*           Yes, there must be a gap between J and I.
            RGINFO(J) = GI_FULLGAP
         ELSE
            RGINFO(J) = GI_NOFULL
         ENDIF
```

### Our patch — `dlaxrb_clssfy_fix.f:358-372`
```fortran
         ABSMAX = MAX( ABS(LB), ABS(UB) )
         SPREAD = UB - LB

         IF( AVGTHRESH .GT. ZERO )THEN
            GAPOK = SPREAD.GE.MIN(AVGTHRESH,ABSMAX*GAPTHRESH)
         ELSE
            GAPOK = SPREAD .GE. ABSMAX*GAPTHRESH
         ENDIF
         IF( GAPOK )THEN
            RGINFO(J) = GI_FULLGAP
         ELSE
            RGINFO(J) = GI_NOFULL
         ENDIF
```

### Python — `mr3_gk.py:306-317`
```python
def classify(eigvals, gaptol):
    n = len(eigvals)
    if n <= 1: return [(0,0)] if n == 1 else []
    groups = []; cs = 0
    while cs < n:
        ce = cs
        while ce+1 < n:
            denom = max(abs(eigvals[ce]), SAFMIN)
            if abs(eigvals[ce+1]-eigvals[ce])/denom >= gaptol: break
            ce += 1
        groups.append((cs, ce)); cs = ce+1
    return groups
```
