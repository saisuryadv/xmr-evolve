"""Post-process the 12 paper-norm sweep logs (no reruns) and emit:

  1. Per-row CSV per (suite, method) with an extra column  ratio = t_eval / t_dbdsqr.
  2. Markdown stats table per suite: min / mean / max for every numeric column,
     broken out per method (advisor / selfcontained / dbdsqr) and per status set
     (ALL rows then PASS-only).
"""
import os, re, glob, math, csv
import statistics as S

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP_DIR = os.path.join(HERE, 'docs', 'sweeps_2026-06-18')
OUT_DIR = os.path.join(HERE, 'docs', 'sweeps_2026-06-18', 'stats')
os.makedirs(OUT_DIR, exist_ok=True)

NUM = r"(?:[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?|inf|nan|-inf|-nan)"

# Matrix-name allows internal spaces (formatted with {:<48}). Pin to the
# numeric tail to distinguish from prelude / banner lines.
ROW = re.compile(
    r"^\s*(?P<name>\S(?:.*?\S)?)\s+"
    r"(?P<n>\d+)\s+"
    r"(?P<res>" + NUM + r")\s+"
    r"(?P<ortU>" + NUM + r")\s+"
    r"(?P<ortV>" + NUM + r")\s+"
    r"(?P<drift>" + NUM + r")\s+"
    r"(?P<t_eval>" + NUM + r")\s+"
    r"(?P<t_dbdsqr>" + NUM + r")\s+"
    r"(?P<reps_e>\d+)/(?P<reps_r>\d+)\s+"
    r"(?P<status>PASS|FAIL|TIMEOUT)"
    r"(?:\s+(?P<note>.*))?$"
)

# Failing-block row (printed for every FAIL even when periodic-print is off):
#   <name>  n=<n>  res=<res>  ortU=<ou>  ortV=<ov>  <note>  drift=<d>  t_eval=<e>  t_dbdsqr=<r>
FAIL_ROW = re.compile(
    r"^\s+(?P<name>\S(?:.*?\S)?)\s+"
    r"n=\s*(?P<n>\d+)\s+"
    r"res=\s*(?P<res>" + NUM + r")\s+"
    r"ortU=\s*(?P<ortU>" + NUM + r")\s+"
    r"ortV=\s*(?P<ortV>" + NUM + r")\s+"
    r"(?P<note>.*?)\s+"
    r"drift=\s*(?P<drift>" + NUM + r")\s+"
    r"t_eval=\s*(?P<t_eval>" + NUM + r")\s+"
    r"t_dbdsqr=\s*(?P<t_dbdsqr>" + NUM + r")\s*$"
)

SUITES = ['pract', 'synth', '379', 'dense_to_bidiag']
METHODS = ['advisor', 'selfcontained', 'dbdsqr']
COLS = ['n', 'res', 'ortU', 'ortV', 'drift', 't_eval', 't_dbdsqr', 'ratio']


def to_float(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')


def _finalize(row):
    te, tr = row['t_eval'], row['t_dbdsqr']
    row['ratio'] = (te / tr) if (math.isfinite(te) and
                                  math.isfinite(tr) and tr > 0) \
                              else float('nan')
    return row


def parse(path):
    rows = []
    seen = set()
    in_failing = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('Failing ('):
                in_failing = True
                continue
            # First try the in-flight columnar row format.
            m = ROW.match(line)
            if m:
                d = m.groupdict()
                row = _finalize({
                    'name':   d['name'].strip(),
                    'n':      int(d['n']),
                    'res':    to_float(d['res']),
                    'ortU':   to_float(d['ortU']),
                    'ortV':   to_float(d['ortV']),
                    'drift':  to_float(d['drift']),
                    't_eval': to_float(d['t_eval']),
                    't_dbdsqr': to_float(d['t_dbdsqr']),
                    'reps_e': int(d['reps_e']),
                    'reps_r': int(d['reps_r']),
                    'status': d['status'],
                    'note':   (d['note'] or '').strip(),
                })
                if row['name'] not in seen:
                    seen.add(row['name'])
                    rows.append(row)
                continue
            if not in_failing:
                continue
            m = FAIL_ROW.match(line)
            if not m:
                continue
            d = m.groupdict()
            row = _finalize({
                'name':   d['name'].strip(),
                'n':      int(d['n']),
                'res':    to_float(d['res']),
                'ortU':   to_float(d['ortU']),
                'ortV':   to_float(d['ortV']),
                'drift':  to_float(d['drift']),
                't_eval': to_float(d['t_eval']),
                't_dbdsqr': to_float(d['t_dbdsqr']),
                'reps_e': -1,
                'reps_r': -1,
                'status': 'FAIL',
                'note':   d['note'].strip(),
            })
            if row['name'] in seen:
                continue
            seen.add(row['name'])
            rows.append(row)
    return rows


def stats(values):
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return None
    return (min(finite), S.fmean(finite), max(finite), len(finite),
            len(values) - len(finite))


def fmt(s):
    if s is None:
        return ('       -', '       -', '       -')
    mn, mean, mx, _, _ = s
    return (f"{mn:>10.3e}", f"{mean:>10.3e}", f"{mx:>10.3e}")


def write_csv(rows, path):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'n', 'res', 'ortU', 'ortV', 'sv_drift',
                    't_eval', 't_dbdsqr', 'ratio', 'reps_eval',
                    'reps_ref', 'status', 'note'])
        for r in rows:
            w.writerow([r['name'], r['n'], r['res'], r['ortU'], r['ortV'],
                        r['drift'], r['t_eval'], r['t_dbdsqr'], r['ratio'],
                        r['reps_e'], r['reps_r'], r['status'], r['note']])


def main():
    md = []
    md.append("# BidiagonalSVD sweep statistics (2026-06-24)\n")
    md.append("Post-processed from `docs/sweeps_2026-06-18/*.log` — no reruns.\n")
    md.append("New per-row column **`ratio = t_eval / t_dbdsqr`** "
              "(solver wall time vs DBDSQR σ-only reference).\n")
    md.append("Per-row CSVs (with the new ratio column) live in "
              "`docs/sweeps_2026-06-18/stats/<suite>_<method>.csv`.\n")
    md.append("\nStats key per column: **min / mean / max** over rows where the value "
              "is finite. Two row sets per suite × method:\n")
    md.append("- **ALL** rows that matched the row schema (includes FAIL inf/nan, which "
              "are excluded from each column's min/mean/max via a finite mask).\n")
    md.append("- **PASS** rows only.\n")
    md.append("\n**Caveat — synth coverage**: the orchestrator only printed one "
              "in every 100 PASS rows for synth (full periodic print would have "
              "made the log ~ 18 k lines × 3 methods). All 5 038 / 691 / 3 FAILs "
              "are present in the synth advisor / selfcontained / dbdsqr CSVs "
              "(parsed from the trailing `Failing (…)` block), so synth PASS "
              "stats below are a ~1 % sample but synth FAIL stats are exact. "
              "The other three suites are 100 % complete.\n")

    for suite in SUITES:
        md.append(f"\n## {suite}\n")
        for method in METHODS:
            log = os.path.join(SWEEP_DIR, f'{suite}_{method}_paper_norms.log')
            if not os.path.exists(log):
                md.append(f"\n### {method}\n_no log_\n")
                continue
            rows = parse(log)
            write_csv(rows, os.path.join(OUT_DIR, f'{suite}_{method}.csv'))
            passes = [r for r in rows if r['status'] == 'PASS']
            md.append(f"\n### {method}  — rows={len(rows)}  PASS={len(passes)}\n")
            md.append("\n| column | set | min | mean | max | n_finite | n_nonfinite |")
            md.append("\n|---|---|---:|---:|---:|---:|---:|")
            for col in COLS:
                for label, src in (('ALL', rows), ('PASS', passes)):
                    vals = [r[col] for r in src]
                    s = stats(vals)
                    if s is None:
                        md.append(f"\n| {col} | {label} | - | - | - | 0 | "
                                  f"{len(vals)} |")
                    else:
                        mn, mean, mx, nf, nn = s
                        md.append(f"\n| {col} | {label} | "
                                  f"{mn:.3e} | {mean:.3e} | {mx:.3e} | "
                                  f"{nf} | {nn} |")
            md.append("\n")

    out_md = os.path.join(OUT_DIR, 'stats_summary.md')
    with open(out_md, 'w') as fh:
        fh.write(''.join(md))
    print(f"wrote {out_md}")
    for f in sorted(os.listdir(OUT_DIR)):
        if f.endswith('.csv'):
            print(f"  csv:  {os.path.join(OUT_DIR, f)}")


if __name__ == '__main__':
    main()
