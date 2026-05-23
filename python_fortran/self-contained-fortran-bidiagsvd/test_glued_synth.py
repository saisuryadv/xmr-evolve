#!/usr/bin/env python3
"""
Systematic glued variant test suite for MR3-GK bidiagonal SVD.

Mirrors the Willems-Lang 2013 Synth methodology: ALL 90 base patterns
tested in glued-small and glued-medium variants, plus larger glued
Wilkinson matrices matching the paper's dimension range.

Reference: Willems & Lang, "A Framework for the MR3 Algorithm," SIAM J.
Sci. Comput. 35:2 (2013), Section 8.1. Their 116,874 Synth testset uses
each base type in 3 variants: as-is, glued-small (~n*eps*||T||),
glued-medium (~n*sqrt(eps)*||T||).
"""
import numpy as np
import sys
import os
import time
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_eval import make, adv_names
from test_dense_to_bidiag import (glue_bidiagonals, test_one, wilkinson_tridiag,
                                   tridiag_to_bidiag)

EPS = np.finfo(np.float64).eps
SQRT_EPS = np.sqrt(EPS)
RES_THRESH = 7.0
ORTHO_THRESH = 5.0

# Patterns to skip (don't make sense to glue)
SKIP_PATTERNS = {'single_element', 'diagonal_only', 'chkbd_4', 'chkbd_16'}


def make_glued(base_name, n_block, glue_type):
    """Generate a glued bidiagonal from a base pattern.
    glue_type: 'small' (3 copies, gamma=n*eps*||B||) or 'medium' (2 copies, gamma=n*sqrt(eps)*||B||)
    """
    d_base, e_base = make(base_name, n_block)
    bnorm = max(np.max(np.abs(d_base)), 1e-300)
    if len(e_base) > 0:
        bnorm = max(bnorm, np.max(np.abs(e_base)))

    if glue_type == 'small':
        copies = 3
        gamma = n_block * EPS * bnorm
    else:  # medium
        copies = 2
        gamma = n_block * SQRT_EPS * bnorm

    d_list = [d_base.copy() for _ in range(copies)]
    e_list = [e_base.copy() for _ in range(copies)]
    return glue_bidiagonals(d_list, e_list, gamma)


def run_glued_synth(block_sizes, timeout_s=60):
    """Run all 90 base patterns in glued-small and glued-medium variants."""
    results = {}
    total_pass = 0
    total = 0

    glueable = [p for p in adv_names if p not in SKIP_PATTERNS]

    for n_block in block_sizes:
        for glue_type in ['small', 'medium']:
            copies = 3 if glue_type == 'small' else 2
            tag = f"glued_{glue_type}"
            print(f"\n--- {tag}, n_block={n_block} (total n~{n_block*copies}) ---")
            sz_pass = 0
            sz_total = 0

            for base_name in glueable:
                try:
                    d, e = make_glued(base_name, n_block, glue_type)
                except Exception as ex:
                    key = f"{base_name}_{tag}@{n_block}"
                    results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                    total += 1
                    sz_total += 1
                    continue

                n = len(d)
                key = f"{base_name}_{tag}@{n}"
                total += 1
                sz_total += 1

                try:
                    ok, res, ou, ov, dt = test_one(d, e)
                    results[key] = (ok, res, ou, ov, dt)
                    if ok:
                        total_pass += 1
                        sz_pass += 1
                    status = 'PASS' if ok else 'FAIL'
                    if not ok or ou > 3.0:  # Print failures and near-failures
                        print(f"  {base_name:35s} n={n:4d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  {status}")
                except Exception as ex:
                    results[key] = (False, float('inf'), float('inf'), float('inf'), 0)
                    print(f"  {base_name:35s} n={n:4d}  ERROR: {ex}")

            print(f"  subtotal: {sz_pass}/{sz_total} passed")

    return results, total_pass, total


def run_large_wilkinson(timeout_s=120):
    """Run larger glued Wilkinson matrices matching Willems-Lang dim range."""
    configs = [
        ('glued_wilk_3x101_sqrteps', 50, 3, SQRT_EPS),
        ('glued_wilk_3x201_sqrteps', 100, 3, SQRT_EPS),
        ('glued_wilk_5x101_sqrteps', 50, 5, SQRT_EPS),
        ('glued_wilk_3x101_1e14',    50, 3, 1e-14),
        ('glued_wilk_3x201_1e14',    100, 3, 1e-14),
        ('glued_wilk_5x201_sqrteps', 100, 5, SQRT_EPS),
    ]

    results = {}
    total_pass = 0
    total = 0

    print("\n--- Large Glued Wilkinson ---")
    for name, m, copies, gamma in configs:
        d_w, e_w = wilkinson_tridiag(m)
        d_b, e_b = tridiag_to_bidiag(d_w, e_w)
        d_list = [d_b.copy() for _ in range(copies)]
        e_list = [e_b.copy() for _ in range(copies)]
        d, e = glue_bidiagonals(d_list, e_list, gamma)
        n = len(d)
        total += 1

        try:
            ok, res, ou, ov, dt = test_one(d, e)
            results[name] = (ok, res, ou, ov, dt)
            if ok:
                total_pass += 1
            status = 'PASS' if ok else 'FAIL'
            print(f"  {name:40s} n={n:5d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  t={dt:.2f}s  {status}")
        except Exception as ex:
            results[name] = (False, float('inf'), float('inf'), float('inf'), 0)
            print(f"  {name:40s} n={n:5d}  ERROR: {ex}")

    return results, total_pass, total


def main():
    parser = argparse.ArgumentParser(description='Systematic Glued Synth Test Suite')
    parser.add_argument('--quick', action='store_true', help='Quick mode: n_block=50 only')
    parser.add_argument('--no-large', action='store_true', help='Skip large Wilkinson tests')
    parser.add_argument('--timeout', type=int, default=60, help='Per-test timeout')
    args = parser.parse_args()

    if args.quick:
        block_sizes = [50]
    else:
        block_sizes = [50, 100]

    all_results = {}
    grand_pass = 0
    grand_total = 0

    # Part A: Systematic glued variants
    print("=" * 80)
    print("PART A: Systematic Glued Variants of All 90 Base Patterns")
    print(f"  Willems-Lang 2013 Synth methodology: glued-small (3 copies, gamma=n*eps*||B||)")
    print(f"                                       glued-medium (2 copies, gamma=n*sqrt(eps)*||B||)")
    print(f"  Block sizes: {block_sizes}")
    n_glueable = len([p for p in adv_names if p not in SKIP_PATTERNS])
    print(f"  {n_glueable} base patterns x 2 glue types x {len(block_sizes)} sizes = {n_glueable * 2 * len(block_sizes)} tests")
    print("=" * 80)

    r, p, t = run_glued_synth(block_sizes, args.timeout)
    all_results.update(r)
    grand_pass += p
    grand_total += t

    # Part B: Large glued Wilkinson
    if not args.no_large:
        print("\n" + "=" * 80)
        print("PART B: Large Glued Wilkinson (matching Willems-Lang dim range)")
        print("=" * 80)

        r, p, t = run_large_wilkinson(args.timeout)
        all_results.update(r)
        grand_pass += p
        grand_total += t

    # Summary
    print("\n" + "=" * 80)
    print(f"GRAND TOTAL: {grand_pass}/{grand_total} passed ({100 * grand_pass / max(grand_total, 1):.1f}%)")

    failed = [(k, v) for k, v in sorted(all_results.items()) if not v[0]]
    if failed:
        print(f"\nFailing ({len(failed)}):")
        for key, (ok, res, ou, ov, dt) in failed:
            print(f"  {key:55s}  res={res:10.3f}  ortU={ou:10.3f}  ortV={ov:10.3f}")

        # Show top 20 worst ortU (including passing)
        print(f"\nTop 20 worst ortU (including passing):")
        by_ortu = sorted(all_results.items(), key=lambda x: -x[1][2])
        for key, (ok, res, ou, ov, dt) in by_ortu[:20]:
            status = 'PASS' if ok else 'FAIL'
            print(f"  {key:55s}  ortU={ou:10.3f}  {status}")
    else:
        print("\nAll tests passed!")
    print("=" * 80)


if __name__ == '__main__':
    main()
