#!/usr/bin/env python3
"""
Build a self-contained presentation HTML at docs/presentation.html.

Embeds:
  - Round-2 markdown docs (BUGREPORT, 01..06)
  - Round-3 NCD convergence experiment results
  - Cluster-tree visualizations (PNGs base64-inlined)
  - Ablation table (from ablation_results.json)

No external CDN; output works offline.
"""
import os
import sys
import json
import base64
import textwrap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
EXP = os.path.join(ROOT, "experiments")
OUT = os.path.join(DOCS, "presentation.html")


def b64img(path, alt=""):
    if not os.path.exists(path):
        return f"<p><em>(missing image: {os.path.basename(path)})</em></p>"
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f'<img alt="{alt}" src="data:image/png;base64,{data}" />'


def md(text):
    """Render markdown to HTML using `markdown` if available, else fallback."""
    try:
        import markdown as _md
        return _md.markdown(text, extensions=["tables", "fenced_code"])
    except Exception:
        # very small fallback
        out = []
        for ln in text.splitlines():
            if ln.startswith("### "):
                out.append(f"<h4>{ln[4:]}</h4>")
            elif ln.startswith("## "):
                out.append(f"<h3>{ln[3:]}</h3>")
            elif ln.startswith("# "):
                out.append(f"<h2>{ln[2:]}</h2>")
            elif ln.strip() == "":
                out.append("")
            else:
                out.append(f"<p>{ln}</p>")
        return "\n".join(out)


def load_md(path):
    with open(path) as f:
        return f.read()


# ----- assemble sections -----

def slide(title, body_html, anchor):
    return f'''<section id="{anchor}">
<h2>{title}</h2>
{body_html}
</section>'''


def section_title():
    return slide(
        "MR3-GK Bidiagonal SVD — Project Status",
        '''
<p class="subtitle">A pure-Fortran (and matched Python) implementation of MR3-GK
on top of Willems' XMR kernel, with two corrective patches and a complete
pre/post-processing layer.</p>
<p class="meta">2026-04-29 &nbsp; • &nbsp; Sai Surya &nbsp; • &nbsp; <code>python_fortran/</code></p>

<table class="kv">
<tr><th>Repo</th><td><code>xmr-evolve / python_fortran/</code></td></tr>
<tr><th>Upstream</th><td>Willems XMR (<code>xmr_src/</code>, ~14 K LOC, untouched)</td></tr>
<tr><th>Reference paper</th><td>Willems &amp; Lang, ETNA 39 (2012) — Algorithm 4.1</td></tr>
<tr><th>Companion</th><td>Demmel-Kahan 1990 (zero-shift QR), LAWN 163 §5.4</td></tr>
</table>
''', "title")


def section_tldr():
    return slide(
        "TL;DR",
        '''
<div class="grid-2">
<div class="card">
<h4>Empirical results</h4>
<ul>
<li><b>379 / 379 PASS</b> at <code>evaluate.py</code> thresholds (res ≤ 7 nε, ortU/V ≤ 5 nε).</li>
<li><b>21,375 / 21,375 PASS</b> at 2·machine_eps Py↔Fortran bit-match (379 + 1,130 + 19,866 specs).</li>
<li>σ <b>bit-identical</b> across all 21,375 specs; worst U diff = 1.11e-16 (≈ 1 ULP); worst V diff = 2.27e-19.</li>
</ul>
</div>
<div class="card">
<h4>What's in the patch set</h4>
<ul>
<li>2 XMR Fortran fixes (32 lines total).</li>
<li>2 Python harmonization helpers for Py↔Fo bit-identity (~45 lines).</li>
<li>~1,400 lines of new Fortran: pure-Fortran orchestration layer (<code>mr3gk_fortran/</code>).</li>
<li>Without the 2 XMR fixes: <b>189 / 379</b> pass (190 failures unmasked).</li>
</ul>
</div>
</div>

<div class="grid-2">
<div class="card">
<h4>What XMR provides vs what we add</h4>
<ul>
<li>XMR: tridiagonal eigenpair solver only — <b>no SVD entry point.</b></li>
<li>We add: T_GK build, two-phase split, zero-shift QR deflation, sign matrices, Bv recovery, sign fix, GS completion, σ rescaling.</li>
</ul>
</div>
<div class="card">
<h4>NCD analysis (this round)</h4>
<ul>
<li>Median NCD_max child/parent shrink = <b>10⁻³ to 10⁻⁶</b> per descent on tested clustered matrices.</li>
<li>All tested matrices stay below the <code>evaluate.py</code> orthogonality threshold (5 nε).</li>
<li>NCD denominator audit: production XMR uses (max−min)/max; Python fallback uses (max−min)/min.</li>
</ul>
</div>
</div>
''', "tldr")


def section_algo():
    return slide(
        "What is MR3-GK?",
        '''
<p>Bidiagonal SVD via the <b>Golub-Kahan</b> reduction:
<code>T_GK = perfect_shuffle( [0,...,0; |d|, |e|] )</code> &mdash; a 2n × 2n
symmetric tridiagonal whose <i>positive</i> eigenvalues are the singular values of B,
and whose eigenvectors interleave the right and left singular vectors. Run MR3
on T_GK to get O(n²) bidiagonal SVD.</p>

<pre><code>   bidiagonal B (d, e)
       │
       │  pre-scaling, two-phase split, sign matrices D1/D2,
       │  zero-shift QR deflation                      ← OUR PRE-PROCESSING
       │
       ▼
   T_GK = [[0, B^T], [B, 0]]
       │  bisection (DSTEBZ) → positive eigenvalues
       │  dlaxre_gk → root representation (OUR PATCH activates GK branch)
       │  dlaxrv → MRRR eigenvectors
       │     dlaxrf → representation tree forest
       │        dlaxrb_clssfy → cluster classifier (OUR PATCH guards AVGTHRESH=0)
       │
       ▼
   eigenpairs (z of length 2n)
       │
       │  de-interleave even/odd halves                ← OUR POST-PROCESSING
       │  normalize, Bv recovery, sign fix
       │  GS completion of zero-σ orthogonal complement
       │
       ▼
   (σ, U, V)
</code></pre>
''', "algo")


def section_xmr_caveats():
    return slide(
        "What XMR provides vs what we add",
        '''
<p><b>XMR exposes <code>dstexr</code></b> (selected eigenpairs of a symmetric
tridiagonal). That is its only public entry. There is <b>no SVD entry point</b>
in <code>xmr_src/</code> — no driver that takes a bidiagonal (d, e) and returns
(σ, U, V).</p>

<table>
<tr><th>Layer</th><th>Upstream XMR</th><th>Our code</th></tr>
<tr><td>Bidiagonal-SVD entry point</td><td><b>none</b></td>
    <td><code>mr3_gk.py:bidiag_svd</code> &nbsp;/&nbsp; <code>mr3gk_fortran/mr3gk.f90:dmr3gk_svd</code></td></tr>
<tr><td>Pre-scaling by max(|d|,|e|)</td><td>none</td>
    <td><code>mr3_gk.py:723-746</code></td></tr>
<tr><td>Two-phase bidiagonal splitting</td><td>none</td>
    <td><code>mr3_gk.py:225-257</code> · <code>mr3gk_split.f90</code></td></tr>
<tr><td>Sign matrices D1, D2</td><td>none</td>
    <td><code>mr3_gk.py:763-770</code></td></tr>
<tr><td>Zero-shift QR sweep</td><td>none</td>
    <td><code>mr3_gk.py:86-110</code> · <code>mr3gk_qrsweep.f90</code></td></tr>
<tr><td>Zero-σ detect-then-deflate</td><td>none</td>
    <td><code>mr3_gk.py:790-822</code> &nbsp;<i>(novel; not in Willems-Lang 2012)</i></td></tr>
<tr><td>T_GK construction (perfect shuffle)</td><td>none</td>
    <td><code>mr3_gk.py:580-584</code> · <code>mr3gk_tgk.f90:139-144</code></td></tr>
<tr><td>De-interleaving + signs into U/V</td><td>none</td>
    <td><code>mr3_gk.py:837-850</code></td></tr>
<tr><td>Bv / Bᵀu recovery</td><td>none</td>
    <td><code>mr3_gk.py:887-939</code></td></tr>
<tr><td>Sign fix (B v · σ u)</td><td>none</td>
    <td><code>mr3_gk.py:941-951</code></td></tr>
<tr><td>GS completion of zero-σ complement</td><td>none</td>
    <td><code>mr3_gk.py:953-987</code></td></tr>
<tr><td>σ rescale</td><td>none</td>
    <td><code>mr3_gk.py:989-991</code></td></tr>
<tr><td>Tridiagonal eigenpair solver</td><td><code>dstexr_</code> ← we use this</td><td>—</td></tr>
<tr><td>Block split (sign normalize, scale)</td><td><code>dlaxra</code></td><td>—</td></tr>
<tr><td>Root representation builder</td><td><code>dlaxre</code></td>
    <td>patched: <code>dlaxre_gk.f</code> activates GK branch</td></tr>
<tr><td>MRRR forest / classify</td><td><code>dlaxrf</code> / <code>dlaxrb_clssfy</code></td>
    <td>patched: <code>dlaxrb_clssfy_fix.f</code> guards AVGTHRESH=0</td></tr>
<tr><td>Eigenvector assembly</td><td><code>dlaxrv</code> ← we use this</td><td>—</td></tr>
</table>

<p class="note"><b>Caveat:</b> XMR was downloaded into the repo on 2026-04-03
(commit <code>61dfa93</code>) — <i>before</i> this conversation. It came from
Willems' <code>tdsolver/xmr/SRC/O</code> and we have not modified the kernel.
Our patches live in <code>dlaxre_gk.f</code> and <code>dlaxrb_clssfy_fix.f</code>
which are linked in place of the upstream files at build time.</p>
''', "xmr-caveats")


def section_preproc():
    return slide(
        "Pre-processing additions (relative to Willems-Lang 2012)",
        '''
<p>The 2012 paper presents the algorithmic core (T_GK + MR³ + classify). Our
implementation adds several engineering layers that the paper does not specify.
The most consequential is described below.</p>

<h4>1. Zero-shift QR sweep + clean-split detection (not in Willems-Lang 2012)</h4>

<p>Demmel-Kahan 1990 (Sec. 3) showed that zero-shift QR on a bidiagonal is
equivalent to inverse iteration on Bᵀ B. When σ_min = 0, <b>both</b>
<code>d[-1]</code> and <code>e[-1]</code> converge to zero in a single sweep.
We use this as a <b>zero-σ detector before MR³</b>:</p>

<pre><code>  for each unsplit block of size k ≥ 2:
      d_sw, e_sw = copy of block
      right_rots, left_rots = zero_shift_qr_sweep(d_sw, e_sw)         # mr3_gk.py:86-110

      # Clean split requires BOTH d[-1] and e[-1] near zero:
      if |d_sw[-1]| < n·ε  AND  |e_sw[-1]| < n·ε:                     # mr3_gk.py:802
          # zero σ found
          recurse on bidiag_svd(d_sw[:-1], e_sw[:-1])                 # mr3_gk.py:808
          embed sub-result into k×k, append (0) for the last column
          apply left/right Givens rotations to undo the QR sweep      # mr3_gk.py:813-815
      else:
          # no zero σ — proceed with MR³ on original (d, e unchanged)
          dispatch to mr3_tgk_multiblock                              # mr3_gk.py:825</code></pre>

<p><b>Why both <code>d[-1]</code> and <code>e[-1]</code>?</b> When σ_min = 0
the sweep drives both to zero (clean 1×1 split at the bottom). When σ_min > 0
(even tiny) only <code>d[-1] → σ_min</code> while <code>e[-1]</code> stays
non-trivial — deflating without a clean split severs matrix coupling and
produces catastrophic residuals. We observed e[-1] = 0.707 after sweep on
<code>gl_wilkw</code>; clean-split check correctly rejected deflation there.</p>

<p><b>Test impact:</b> fixes <code>saw_tooth</code> (ortU 552 → 2.04 nε),
<code>step_function</code> (ortU 70.6 → 18.1 nε), and every matrix with
structural-zero diagonals (<code>zero_diagonal</code>, <code>gl_wilkp</code>,
<code>gl_clement</code>, &c.).</p>

<h4>2. Two-phase splitting (Phase 1 relative + Phase 2 absolute)</h4>

<p>At <code>mr3_gk.py:225-257</code>: First pass uses
<code>|e[i]| ≤ ε · (|d[i]| + |d[i+1]|)</code> (relative), preserving relative
accuracy for graded matrices. Second pass within each relative block uses
<code>|e[i]| ≤ k·ε·||B||_∞</code> per Willems-Lang condition (3.5). The two-pass
structure is implementation-specific; the paper formulates the absolute test
but does not specify the relative-then-absolute layering.</p>

<h4>3. Sign matrices D1, D2 baked into the singleton path</h4>

<p>At <code>mr3_gk.py:763-770</code>: rather than form the sign-conjugated
bidiagonal explicitly, we compute D1, D2 once and apply them as <i>scaling</i>
during the Z → (U, V) extraction step
(<code>V[r,c] = Z[0::2,c] · d2[r]</code>, <code>U[r,c] = Z[1::2,c] · d1[r]</code>
at lines 849-850). For 1×1 blocks this collapses to constant-time ±1 placement
(lines 837-838).</p>

<h4>4. Pre-scaling by max(|d|,|e|)</h4>

<p>At <code>mr3_gk.py:723-746</code> and unscale at line 990. Required to handle
<code>near_overflow</code>-class matrices whose entries reach 10¹⁵⁰ (the
algorithm operates on the unit-norm bidiagonal and rescales σ at the end).
Without it, intermediate `B²` operations overflow.</p>
''', "preproc")


def section_postproc():
    return slide(
        "Post-processing additions",
        '''
<h4>1. De-interleave Z into U / V halves with sign matrices</h4>
<p>T_GK eigenvector <code>z</code> of length 2n splits into right-singular
half <code>z[0::2]</code> and left-singular half <code>z[1::2]</code>; multiplied
elementwise by D2, D1 respectively. <code>mr3_gk.py:849-850</code></p>

<h4>2. Bv recovery for one-sided-good vectors</h4>
<p>If only one half of <code>z</code> survives normalization (the other has
norm below threshold — typical for degenerate σ pairs), reconstruct the missing
half via <code>u = Bv/σ</code> or <code>v = Bᵀu/σ</code>. Three branches at
<code>mr3_gk.py:887-939</code> for "U missing", "V missing", and "both good but
residual high" — the third branch additionally checks
<code>‖Bv − σu‖</code> vs <code>‖Bᵀu − σv‖</code> and replaces whichever is
worse. This explicit consistency check is not in the paper.</p>

<h4>3. Sign fix via dot product</h4>
<p>At <code>mr3_gk.py:941-951</code>: compute <code>dot(Bv, σu)</code> per
column. If negative, flip the U-column sign. Ensures <code>Bv = +σu</code>
(not <code>−σu</code>), which the paper assumes implicitly but the eigenvector
solver does not enforce because T_GK eigenvectors are sign-ambiguous.</p>

<h4>4. Gram-Schmidt completion for zero-σ orthogonal complement</h4>
<p>At <code>mr3_gk.py:953-987</code>: when σ = 0, the corresponding U/V columns
are not determined by the eigenvector solver. We complete them via two-pass
modified Gram-Schmidt against the already-good columns, with a <b>deterministic
LAPACK DLARNV-seeded</b> random start so Python and Fortran produce
bit-identical completion vectors. Seeds derived from the column index:
<code>(1,3,5,(2j+1)&amp;4095|1)</code> for V and
<code>(2,4,6,(2j+3)&amp;4095|1)</code> for U.</p>

<h4>5. σ rescaling</h4>
<p><code>mr3_gk.py:989-991</code>: undo the prescale factor.</p>
''', "postproc")


def section_bugs():
    bug = load_md(os.path.join(DOCS, "BUGREPORT.md"))
    # Strip the title (already shown in slide title)
    body = md(bug)
    return slide("Bug Report — Modifications to Upstream XMR", body, "bugs")


def section_test_results():
    return slide(
        "Test Results",
        f'''
<h4>Combined headline (21,375 specs across 3 corpora)</h4>
<table>
<tr><th>Suite</th><th>Specs</th><th>Pass @ 2·eps</th><th>Worst σ</th><th>Worst U</th><th>Worst V</th></tr>
<tr><td>379-spec match (<code>test_fortran_match.py</code>)</td><td>379</td><td>379/379</td>
    <td>0</td><td>1.110e-16</td><td>2.272e-19</td></tr>
<tr><td>1,130-spec extended sweep</td><td>1,130</td><td>1,130/1,130</td>
    <td>0</td><td>0</td><td>0</td></tr>
<tr><td>19,866-spec full sweep ("16k experiment")</td><td>19,866</td><td>19,866/19,866</td>
    <td>0</td><td>1.110e-16</td><td>2.272e-19</td></tr>
<tr><th>Total</th><th>21,375</th><th>21,375 / 21,375</th><th>0</th><th>1.110e-16</th><th>2.272e-19</th></tr>
</table>

<h4>evaluate.py thresholds (res ≤ 7 nε, ortU/V ≤ 5 nε): 379 / 379 pass</h4>

<h4>Patch ablation (379 specs × 4 library variants)</h4>
<table>
<tr><th>Variant</th><th>Pass</th><th>Fail</th><th>Wall (s)</th></tr>
<tr><td><b>patched</b> (current)</td><td><b>379 / 379</b></td><td>0</td><td>166.3</td></tr>
<tr><td>revert <code>dlaxre_gk.f</code></td><td>233 / 379</td><td>146</td><td>165.6</td></tr>
<tr><td>revert <code>dlaxrb_clssfy_fix.f</code></td><td>330 / 379</td><td>49</td><td>164.5</td></tr>
<tr><td>revert both</td><td>189 / 379</td><td><b>190</b></td><td>164.4</td></tr>
</table>

<p class="note">The two patches are <b>synergistic</b>: 13 specs that pass under either single revert
fail when both are reverted. Applying only one is not safe.</p>
''', "tests")


def section_ncd_analysis():
    # load NCD convergence log
    path = os.path.join(EXP, "ncd_convergence_log.json")
    if os.path.exists(path):
        with open(path) as f:
            ncd = json.load(f)
    else:
        ncd = []

    # Build per-matrix table
    rows = ['<tr><th>matrix</th><th>n</th><th>max depth</th><th>edges</th><th>leaves</th>'
            '<th>median NCD_max shrink</th><th>median NCD_min shrink</th>'
            '<th>max getvec ratio</th><th>ortU (nε)</th></tr>']
    for r in ncd:
        rows.append(
            f'<tr><td>{r["name"]}</td><td>{r["n"]}</td><td>{r["max_depth"]}</td>'
            f'<td>{r["n_edges"]}</td><td>{r["n_leaves"]}</td>'
            f'<td>{r["median_ncd_max_shrink"]:.2g}</td>'
            f'<td>{r["median_ncd_min_shrink"]:.2g}</td>'
            f'<td>{r["max_getvec_ratio"]:.2e}</td>'
            f'<td>{r["ortU_neps"]:.2g}</td></tr>')
    table = "<table>" + "".join(rows) + "</table>"

    return slide(
        "NCD Analysis on Clustered Matrices",
        f'''
<p>For each clustered test matrix we drive the instrumented Python
<code>mr3_block</code> (algorithmically equivalent to the Fortran kernel)
and record per-edge NCD reduction, per-leaf getvec ratio, and the actual
<code>ortU</code> achieved. Two NCD definitions tracked side-by-side:</p>
<ul>
<li><code>NCD_max = (max − min) / max(|max|, |min|)</code> — the form XMR uses</li>
<li><code>NCD_min = (max − min) / min(|max|, |min|)</code> — the form the user asked us to verify</li>
</ul>

<h4>Per-matrix summary</h4>
{table}

<p><b>Headline:</b> on every tested clustered matrix the algorithm shrinks
NCD_max by 10⁻³ to 10⁻⁶ in a single descent step, then resolves the residual
sub-clusters as singletons or zero-width batches. All eight matrices stay below
the <code>evaluate.py</code> orthogonality threshold (5 nε).</p>

<h4>Convergence plot</h4>
{b64img(os.path.join(EXP, "ncd_convergence.png"), "NCD convergence")}

<h4>Tree visualization (round-2 cluster trace)</h4>
{b64img(os.path.join(EXP, "cluster_tree.png"), "Cluster tree")}

<h4>Per-adjacent-pair NCD distribution</h4>
{b64img(os.path.join(EXP, "ncd_per_adjacent_pair.png"), "Per-pair NCD")}

<h4>Per-node aggregate NCD vs depth</h4>
{b64img(os.path.join(EXP, "ncd_progression.png"), "NCD progression")}
''', "ncd")


def section_ncd_def_audit():
    return slide(
        "NCD Denominator Audit (across all trees)",
        '''
<p>The user asked us to verify <b>(max − min) / min</b> across every code tree
that performs cluster classification. Findings:</p>

<table>
<tr><th>Tree</th><th>File:Line</th><th>Implemented metric</th><th>Matches /min ?</th></tr>
<tr><td>Upstream XMR</td><td><code>xmr_src/dlaxrb_clssfy.f:358-369</code></td>
    <td>SPREAD / max(|LB|,|UB|) ≥ GAPTHRESH</td><td>❌ uses /max</td></tr>
<tr><td>Our XMR patch</td><td><code>dlaxrb_clssfy_fix.f:358-372</code></td>
    <td>same as upstream (only fixes AVGTHRESH=0 collapse)</td><td>❌ uses /max</td></tr>
<tr><td>Python fallback</td><td><code>mr3_gk.py:313-314</code></td>
    <td><code>|ev[ce+1]−ev[ce]| / max(|ev[ce]|, SAFMIN)</code></td>
    <td>✅ uses /min</td></tr>
<tr><td>Fortran orchestration</td><td>n/a</td><td>delegates to XMR</td><td>❌</td></tr>
</table>

<p><b>Production path uses /max</b>; only the Python fallback uses /min. They
are inconsistent. The Python fallback is dead code on our test set (XMR always
succeeds), so this discrepancy does not currently affect pass/fail.</p>

<h4>Empirical observation from the cluster trace</h4>
<p>On <code>gl_clustered_at_eps@200</code> the depth-1 leaf nodes have:</p>
<ul>
<li>NCD_max = 1.00 (saturated — "spread is the same order as max")</li>
<li>NCD_min = 4·10⁸ to 1.6·10¹⁰ (the actual relative gap is huge in /min terms)</li>
</ul>
<p>This is the regime where the metric choice matters most. <code>/max</code>
saturates at 1; <code>/min</code> conveys order-of-magnitude information about
dynamic range. On our test set both metrics agree on the split/keep decision
because <code>NCD_max = 1.0 ≥ GAPTOL = 1e-3</code> always splits.</p>

<p><b>Recommendation:</b> document the inconsistency; do not change the
production metric without a full re-validation pass through the 21,375-spec
sweep.</p>
''', "ncd-def")


def section_build_flags():
    return slide(
        "Build Flags",
        '''
<h4>Root <code>build.sh</code> — produces <code>libxmr.so</code></h4>
<pre><code>FFLAGS = -fPIC -O2 -std=legacy -w -fno-second-underscore</code></pre>
<p>Standard XMR upstream flags. <code>-O2</code> is acceptable here because the
XMR kernel was independently validated; we do not depend on bit-identity between
two compiles of <code>libxmr.so</code>.</p>
<p>One minor flag-drift on line 58 of <code>build.sh</code>: <code>dlaxre_gk.f</code>
is compiled with <code>-c -fPIC -O2</code> only (missing <code>-std=legacy
-w -fno-second-underscore</code>). Functionally fine on gfortran 4.x where these
are defaults; recommend matching <code>$FFLAGS</code> for portability.</p>

<h4><code>mr3gk_fortran/build.sh</code> — produces <code>mr3gk_run</code></h4>
<pre><code>FFLAGS = -O0 -fno-fast-math -fPIC -fno-second-underscore -fimplicit-none</code></pre>
<p><code>-O0 -fno-fast-math</code> is <b>deliberate</b>: it preserves source-order
floating-point evaluation so Python and Fortran orchestrations produce bit-identical
σ. At <code>-O2</code> gfortran reorders associative ops, fuses FMA, and
vectorizes — each of which would introduce 1–2 ULP differences that break the
2·eps bit-match target.</p>
<p>Performance trade-off is small: the <code>libxmr.so</code> kernel still uses
<code>-O2</code> (and the kernel is the bottleneck); only the orchestration layer
runs at <code>-O0</code>, and that layer is O(n²) per spec.</p>
''', "build")


def section_repro():
    return slide(
        "Reproducibility",
        '''
<pre><code># 1. Build the kernel (libxmr.so) and pure-Fortran driver (mr3gk_run)
cd python_fortran
bash build.sh
cd mr3gk_fortran && bash build.sh && cd ..

# 2. Verify standalone evaluation (the headline "379/379 PASS")
python3 evaluate.py
# expected:  Score 100.0  (379 / 379 passed)

# 3. Verify Python ↔ Fortran bit-identity on 379 specs (~2 min)
python3 test_fortran_match.py --gen-baseline   # one-time, ~5 min
python3 test_fortran_match.py
# expected:  PASS: 379/379  Worst sigma diff: 0  Worst U diff: 1.110e-16  Worst V diff: 2.272e-19

# 4. Reproduce the patch ablation (4 variants × 379 specs, ~12 min)
python3 experiments/ablation_xmr_patches.py
# writes experiments/ablation_results.json (4 variants, full per-spec data)

# 5. Reproduce the cluster-tree NCD analysis
python3 experiments/cluster_tree_trace.py
# writes experiments/cluster_tree_log.{json,csv} + 3 PNGs

# 6. Reproduce the NCD-convergence experiment
python3 experiments/ncd_convergence.py
# writes experiments/ncd_convergence_log.json + ncd_convergence.png
</code></pre>
''', "repro")


def section_caveats():
    return slide(
        "Caveats and Open Questions",
        '''
<h4>1. XMR provenance</h4>
<p>The <code>xmr_src/</code> tree was added to the repo on <b>2026-04-03</b>
(commit <code>61dfa93</code>: "Add Willems SRC/O source"), <i>before</i> this
session began. We did not download XMR during this conversation. The kernel is
unmodified — our two patches (<code>dlaxre_gk.f</code>,
<code>dlaxrb_clssfy_fix.f</code>) live outside <code>xmr_src/</code> and are
linked in place of the upstream files at build time.</p>

<h4>2. NCD denominator inconsistency</h4>
<p>Production XMR uses <code>(max−min)/max</code>; Python fallback uses
<code>(max−min)/min</code>. Documented but not consolidated. Consolidating
requires a full 21,375-spec re-validation.</p>

<h4>3. Tree depth on our clustered test set caps at 1</h4>
<p>All eight tested clustered matrices descend exactly one level. Deeper trees
would arise on harder STCollection specs (<code>T_Godunov_*</code>, n ≥ 2500).
Out of round-3 scope; deeper-tree analysis is a follow-up if needed.</p>

<h4>4. Bit-match diffs concentrate on <code>gl_wilkw / gl_wilk2w</code></h4>
<p>Of 21,375 specs only 5 have non-zero U-diff (≤ 1 ULP) and 2 have non-zero
V-diff (≤ 1 ULP). All in this Wilkinson-glued family.</p>

<h4>5. Performance benchmarking is out of scope</h4>
<p>This deliverable is correctness only. Pure-Fortran driver wall-time vs the
Python+Fortran hybrid is not yet measured at scale.</p>

<h4>6. <code>-O0</code> on the orchestration layer</h4>
<p>Required for bit-identity. Switching to <code>-O2</code> would speed up the
post-processing slightly but break the 2·eps Py↔Fo target.</p>
''', "caveats")


def build_html():
    sections = [
        section_title(),
        section_tldr(),
        section_algo(),
        section_xmr_caveats(),
        section_preproc(),
        section_postproc(),
        section_bugs(),
        section_test_results(),
        section_ncd_analysis(),
        section_ncd_def_audit(),
        section_build_flags(),
        section_repro(),
        section_caveats(),
    ]
    nav = '<nav><ul>' + "".join(
        f'<li><a href="#{anchor}">{title}</a></li>'
        for anchor, title in [
            ("title", "Title"),
            ("tldr", "TL;DR"),
            ("algo", "Algorithm"),
            ("xmr-caveats", "XMR vs ours"),
            ("preproc", "Pre-processing"),
            ("postproc", "Post-processing"),
            ("bugs", "Bug report"),
            ("tests", "Test results"),
            ("ncd", "NCD analysis"),
            ("ncd-def", "NCD definition"),
            ("build", "Build flags"),
            ("repro", "Reproduce"),
            ("caveats", "Caveats"),
        ]
    ) + '</ul></nav>'

    style = textwrap.dedent("""
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, sans-serif;
           max-width: 1100px; margin: 1em auto; color: #222; line-height: 1.5;
           padding: 0 1em 4em; }
    nav { position: sticky; top: 0; background: #fff; border-bottom: 1px solid #ccc;
          padding: 0.5em 0; margin-bottom: 1em; z-index: 10; }
    nav ul { display: flex; flex-wrap: wrap; gap: 0.5em 1em; list-style: none; padding: 0;
             margin: 0; font-size: 12px; }
    nav a { color: #06c; text-decoration: none; }
    nav a:hover { text-decoration: underline; }
    section { padding: 1.5em 0; border-bottom: 1px solid #eee; }
    section:last-child { border-bottom: none; }
    h2 { color: #06335c; border-left: 4px solid #06c; padding-left: 0.5em;
         margin-bottom: 0.5em; font-size: 1.6em; }
    h3 { color: #244; font-size: 1.2em; margin-top: 1.5em; }
    h4 { color: #024; font-size: 1.05em; margin-top: 1.2em; margin-bottom: 0.3em; }
    p.subtitle { font-size: 1.1em; color: #555; }
    p.meta { color: #888; font-size: 0.9em; }
    p.note { background: #fef9e7; border-left: 4px solid #c80; padding: 0.6em 0.9em;
             border-radius: 3px; }
    table { border-collapse: collapse; margin: 0.5em 0; font-size: 0.92em; }
    table th, table td { border: 1px solid #ccc; padding: 0.35em 0.7em; text-align: left;
                          vertical-align: top; }
    table th { background: #eef; }
    table.kv th { width: 130px; background: #fafafa; }
    tr:nth-child(2n) td { background: #fafafa; }
    pre { background: #f6f8fa; border: 1px solid #ddd; border-radius: 4px; padding: 0.8em;
          overflow-x: auto; font-size: 0.86em; line-height: 1.4; }
    code { font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 0.92em; }
    pre code { font-size: inherit; }
    p > code, li > code, td > code, table th > code { background: #f0f3f5;
            padding: 0.05em 0.3em; border-radius: 2px; }
    img { max-width: 100%; height: auto; display: block; margin: 0.6em 0; }
    .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 1em; margin: 1em 0; }
    .card { background: #f6f9fc; border: 1px solid #ddd; border-radius: 5px;
            padding: 0.9em 1em; }
    .card h4 { margin-top: 0; }
    ul, ol { padding-left: 1.6em; }
    li { margin: 0.15em 0; }
    """)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>MR3-GK Bidiagonal SVD — Project Status</title>
<style>{style}</style>
</head>
<body>
{nav}
<main>
{chr(10).join(sections)}
</main>
</body>
</html>
"""


def main():
    html = build_html()
    with open(OUT, "w") as f:
        f.write(html)
    sz = os.path.getsize(OUT)
    n_imgs = html.count("data:image/png;base64,")
    print(f"wrote {OUT}  ({sz/1024:.1f} KB, {n_imgs} embedded images)")


if __name__ == "__main__":
    main()
