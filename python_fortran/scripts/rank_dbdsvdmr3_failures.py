#!/usr/bin/env python3
"""Parse the four DBDSVDMR3 sweep logs in docs/, rank every failing matrix by
severity (max of res, ortU, ortV in n.eps units; +inf for parse-fail/TIMEOUT),
and emit a markdown report with one reproduction command per failure.
"""
import os
import re
import math
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DOCS = os.path.join(ROOT, 'docs')

SUITE_LOGS = [
    ('pract',           'pract_dbdsvdmr3_2026-06-12.log'),
    ('synth',           'synth_dbdsvdmr3_2026-06-12.log'),
    ('379',             'lapack_dbdsvr_379_dbdsvdmr3_2026-06-12.log'),
    ('dense_to_bidiag', 'dense_to_bidiag_dbdsvdmr3_2026-06-12.log'),
]

ROW_RE = re.compile(
    r'^\s+(\S+)\s+n=\s*(\d+)\s+res=\s*(\S+)\s+ortU=\s*(\S+)\s+ortV=\s*(\S+)(?:\s+(.*))?$'
)


def parse_value(s):
    s = s.strip()
    if s == 'inf' or s == 'nan':
        return math.inf
    try:
        return float(s)
    except ValueError:
        return math.inf


def parse_failing_block(path):
    out = []
    in_block = False
    with open(path) as fh:
        for line in fh:
            if line.startswith('Failing'):
                in_block = True
                continue
            if not in_block:
                continue
            if not line.strip():
                continue
            m = ROW_RE.match(line.rstrip())
            if not m:
                continue
            name = m.group(1)
            n = int(m.group(2))
            res = parse_value(m.group(3))
            ou = parse_value(m.group(4))
            ov = parse_value(m.group(5))
            note = (m.group(6) or '').strip()
            out.append((name, n, res, ou, ov, note))
    return out


def repro_cmd(suite, name):
    base = 'python3 eval_dbdsvdmr3.py'
    if suite == 'synth':
        return f"{base} --suite synth --synth-mode paper --only '{name}'"
    if suite == '379':
        return f"{base} --suite 379 --only '{name}'"
    if suite == 'dense_to_bidiag':
        return f"{base} --suite dense_to_bidiag --only '{name}'"
    return f"{base} --suite pract --only '{name}'"


def severity_key(row):
    name, n, res, ou, ov, note = row
    sev = max(res, ou, ov)
    return -sev


def main():
    all_fails = []
    for suite, fname in SUITE_LOGS:
        path = os.path.join(DOCS, fname)
        for row in parse_failing_block(path):
            all_fails.append((suite, *row))

    # Bucket by severity tier
    inf_fails  = [r for r in all_fails if math.isinf(max(r[3], r[4], r[5]))]
    huge       = [r for r in all_fails if not math.isinf(max(r[3], r[4], r[5]))
                  and max(r[3], r[4], r[5]) >= 1e6]
    big        = [r for r in all_fails if not math.isinf(max(r[3], r[4], r[5]))
                  and 1e3 <= max(r[3], r[4], r[5]) < 1e6]
    mild       = [r for r in all_fails if not math.isinf(max(r[3], r[4], r[5]))
                  and max(r[3], r[4], r[5]) < 1e3]

    for bucket in (inf_fails, huge, big, mild):
        bucket.sort(key=lambda r: (-max(r[3], r[4], r[5]), r[0], r[1]))

    out = os.path.join(DOCS, 'dbdsvdmr3_failures_ranked_2026-06-12.md')
    with open(out, 'w') as fh:
        fh.write('# DBDSVDMR3 failures across all 4 suites — ranked by severity\n\n')
        fh.write(f'Total failures parsed from logs: **{len(all_fails)}** '
                 f'(pract 54, 379 40, dense_to_bidiag 32, synth top-50).\n\n')
        fh.write('> Note: synth ran 18438 matrices with 4478 failures total, '
                 'but the original orchestrator capped the "Failing" block at '
                 'top-50 by max(ortU, ortV). The ranker has been updated to '
                 'dump every failure; the full synth list will appear once the '
                 'sweep is re-run. The 50 captured here are already the worst '
                 'cases.\n\n')
        fh.write('Severity = `max(res, ortU, ortV)` in units of `n·eps`. '
                 'Thresholds: res ≤ 7, ortU/ortV ≤ 5.\n\n')
        fh.write('Buckets (decreasing severity):\n')
        fh.write(f'- **catastrophic** (TIMEOUT / parse-fail / INF): {len(inf_fails)}\n')
        fh.write(f'- **huge** (≥ 1e6 n·ε): {len(huge)}\n')
        fh.write(f'- **big** (1e3 .. 1e6 n·ε): {len(big)}\n')
        fh.write(f'- **mild** (< 1e3 n·ε): {len(mild)}\n\n')
        fh.write('Each row carries a one-liner to reproduce it via the orchestrator.\n\n')
        fh.write('```\ncd python_fortran\nbash BidiagonalSVD_TGK/build.sh   # once\n```\n\n')

        def emit(title, rows, show_all=False, limit=80):
            fh.write(f'## {title} ({len(rows)})\n\n')
            if not rows:
                fh.write('_none_\n\n')
                return
            show = rows if show_all else rows[:limit]
            fh.write('| # | suite | matrix | n | res | ortU | ortV | note | repro |\n')
            fh.write('|---|---|---|---|---|---|---|---|---|\n')
            for i, (suite, name, n, res, ou, ov, note) in enumerate(show, 1):
                rs = 'inf' if math.isinf(res) else f'{res:.3g}'
                us = 'inf' if math.isinf(ou) else f'{ou:.3g}'
                vs = 'inf' if math.isinf(ov) else f'{ov:.3g}'
                cmd = repro_cmd(suite, name).replace('|', r'\|')
                fh.write(f'| {i} | {suite} | `{name}` | {n} | {rs} | {us} | {vs} '
                         f'| {note} | `{cmd}` |\n')
            if len(rows) > len(show):
                fh.write(f'\n_…{len(rows) - len(show)} more rows omitted; full list '
                         f'derivable from `docs/*_dbdsvdmr3_2026-06-12.log`._\n')
            fh.write('\n')

        emit('Catastrophic (TIMEOUT / parse-fail / INF)', inf_fails, show_all=True)
        emit('Huge (>= 1e6 n.eps)', huge, show_all=False, limit=100)
        emit('Big (1e3 .. 1e6 n.eps)', big, show_all=False, limit=80)
        emit('Mild (< 1e3 n.eps)', mild, show_all=False, limit=40)

        fh.write('## Reproducing the entire sweep\n\n')
        fh.write('```\ncd python_fortran\n')
        fh.write('bash BidiagonalSVD_TGK/build.sh\n')
        fh.write('python3 eval_dbdsvdmr3.py --suite pract           2>&1 | tee docs/pract_dbdsvdmr3.log\n')
        fh.write('python3 eval_dbdsvdmr3.py --suite synth --synth-mode paper 2>&1 | tee docs/synth_dbdsvdmr3.log\n')
        fh.write('python3 eval_dbdsvdmr3.py --suite 379             2>&1 | tee docs/lapack_dbdsvr_379_dbdsvdmr3.log\n')
        fh.write('python3 eval_dbdsvdmr3.py --suite dense_to_bidiag 2>&1 | tee docs/dense_to_bidiag_dbdsvdmr3.log\n')
        fh.write('```\n')

    print(f'wrote {out}')
    print(f'  catastrophic: {len(inf_fails)}')
    print(f'  huge:         {len(huge)}')
    print(f'  big:          {len(big)}')
    print(f'  mild:         {len(mild)}')


if __name__ == '__main__':
    main()
