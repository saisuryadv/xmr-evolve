#!/usr/bin/env python3
"""Compare raw S,U,V dumps from C and Fortran hgbsvd libraries."""
import sys
import numpy as np

def parse_dump(path):
    """Parse dump file into dict of {(name, n): {'info': int, 'sigma': array, 'U': array, 'V': array}}"""
    results = {}
    with open(path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MATRIX"):
            parts = line.split()
            name = parts[1]
            n = int(parts[2].split("=")[1])
            info = int(parts[3].split("=")[1])

            # Parse SIGMA
            i += 1
            sigma_parts = lines[i].strip().split()[1:]  # skip "SIGMA"
            sigma = np.array([float(x) for x in sigma_parts])

            # Parse U
            i += 1
            u_parts = lines[i].strip().split()[1:]  # skip "U"
            U = np.array([float(x) for x in u_parts])

            # Parse V
            i += 1
            v_parts = lines[i].strip().split()[1:]  # skip "V"
            V = np.array([float(x) for x in v_parts])

            results[(name, n)] = {
                'info': info,
                'sigma': sigma,
                'U': U.reshape(n, n) if len(U) == n*n else U,
                'V': V.reshape(n, n) if len(V) == n*n else V,
            }
        i += 1
    return results

def compare(c_path, f_path):
    c_data = parse_dump(c_path)
    f_data = parse_dump(f_path)

    all_keys = sorted(set(c_data.keys()) | set(f_data.keys()))

    print(f"C library: {len(c_data)} matrices")
    print(f"Fortran library: {len(f_data)} matrices")
    print()

    # Check for matrices present in one but not the other
    only_c = set(c_data.keys()) - set(f_data.keys())
    only_f = set(f_data.keys()) - set(c_data.keys())
    if only_c:
        print(f"ONLY in C ({len(only_c)}):")
        for k in sorted(only_c):
            print(f"  {k[0]} n={k[1]} info={c_data[k]['info']}")
        print()
    if only_f:
        print(f"ONLY in Fortran ({len(only_f)}):")
        for k in sorted(only_f):
            print(f"  {k[0]} n={k[1]} info={f_data[k]['info']}")
        print()

    # Compare common matrices
    common = sorted(set(c_data.keys()) & set(f_data.keys()))
    print(f"Common matrices: {len(common)}")
    print()

    # Info mismatches
    info_mismatches = []
    for k in common:
        if c_data[k]['info'] != f_data[k]['info']:
            info_mismatches.append((k, c_data[k]['info'], f_data[k]['info']))

    if info_mismatches:
        print(f"INFO MISMATCHES ({len(info_mismatches)}):")
        for k, ci, fi in info_mismatches:
            print(f"  {k[0]} n={k[1]}: C={ci}, F={fi}")
        print()
    else:
        print("No info mismatches.")
        print()

    # Detailed comparison for matrices where both succeeded
    print(f"{'Matrix':<30} {'n':>3} | {'MaxΔσ':>12} {'RelΔσ':>12} | {'MaxΔU':>12} {'MaxΔV':>12} | {'Note'}")
    print("-" * 120)

    worst_sigma_abs = 0
    worst_sigma_rel = 0
    worst_U = 0
    worst_V = 0
    worst_sigma_name = ""
    worst_U_name = ""
    worst_V_name = ""

    for k in common:
        cd = c_data[k]
        fd = f_data[k]

        if cd['info'] != 0 or fd['info'] != 0:
            note = f"info: C={cd['info']} F={fd['info']}"
            print(f"{k[0]:<30} {k[1]:>3} | {'N/A':>12} {'N/A':>12} | {'N/A':>12} {'N/A':>12} | {note}")
            continue

        # Sigma comparison
        s_c = cd['sigma']
        s_f = fd['sigma']
        n = min(len(s_c), len(s_f))

        sigma_abs_diff = np.max(np.abs(s_c[:n] - s_f[:n]))
        sigma_max = np.max(np.abs(np.concatenate([s_c[:n], s_f[:n]])))
        sigma_rel_diff = sigma_abs_diff / sigma_max if sigma_max > 0 else 0

        # U comparison — need to handle sign ambiguity per column
        # Each singular vector is determined up to ±1 sign
        U_c = cd['U']
        U_f = fd['U']
        u_diff = 0
        if U_c.shape == U_f.shape:
            nn = U_c.shape[0]
            for col in range(nn):
                # Try both signs
                d1 = np.max(np.abs(U_c[:, col] - U_f[:, col]))
                d2 = np.max(np.abs(U_c[:, col] + U_f[:, col]))
                u_diff = max(u_diff, min(d1, d2))

        # V comparison — same sign ambiguity
        V_c = cd['V']
        V_f = fd['V']
        v_diff = 0
        if V_c.shape == V_f.shape:
            nn = V_c.shape[0]
            for col in range(nn):
                d1 = np.max(np.abs(V_c[:, col] - V_f[:, col]))
                d2 = np.max(np.abs(V_c[:, col] + V_f[:, col]))
                v_diff = max(v_diff, min(d1, d2))

        note = ""
        if sigma_abs_diff > 1e-10 or u_diff > 1e-10 or v_diff > 1e-10:
            note = "LARGE DIFF"
        elif sigma_abs_diff > 0 or u_diff > 0 or v_diff > 0:
            note = "tiny diff"
        else:
            note = "exact"

        print(f"{k[0]:<30} {k[1]:>3} | {sigma_abs_diff:>12.4e} {sigma_rel_diff:>12.4e} | {u_diff:>12.4e} {v_diff:>12.4e} | {note}")

        if sigma_abs_diff > worst_sigma_abs:
            worst_sigma_abs = sigma_abs_diff
            worst_sigma_name = f"{k[0]} n={k[1]}"
        if sigma_rel_diff > worst_sigma_rel:
            worst_sigma_rel = sigma_rel_diff
        if u_diff > worst_U:
            worst_U = u_diff
            worst_U_name = f"{k[0]} n={k[1]}"
        if v_diff > worst_V:
            worst_V = v_diff
            worst_V_name = f"{k[0]} n={k[1]}"

    print()
    print("=== SUMMARY ===")
    print(f"Worst sigma abs diff: {worst_sigma_abs:.4e}  ({worst_sigma_name})")
    print(f"Worst sigma rel diff: {worst_sigma_rel:.4e}")
    print(f"Worst U diff:         {worst_U:.4e}  ({worst_U_name})")
    print(f"Worst V diff:         {worst_V:.4e}  ({worst_V_name})")

if __name__ == "__main__":
    c_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/hgb_dump_c.txt"
    f_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/hgb_dump_f.txt"
    compare(c_path, f_path)
