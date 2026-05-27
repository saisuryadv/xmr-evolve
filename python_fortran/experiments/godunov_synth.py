#!/usr/bin/env python3
"""Generate interactive MR3 tree visualization for the Godunov problem.

Creates docs/godunov_tree.html — a self-contained interactive visualization
showing the MR3 representation tree for a synthetic Godunov bidiagonal,
comparing with/without additive perturbation.

Usage:
    python3 experiments/godunov_synth.py [--n N] [--real 1e-5]

    --n N       Use synthetic Godunov of size N (default: 20)
    --real X    Use real T_Godunov_X from STCollection instead of synthetic
"""

import sys, os, io, subprocess, json, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
EPS = np.finfo(np.float64).eps


def make_synthetic_godunov(n):
    """Synthetic bidiagonal with Godunov period-4 structure."""
    d = np.array([1.0 if i % 2 == 0 else 0.04468810 for i in range(n)])
    e = np.array([0.99900100 if i % 2 == 0 else 2.48388500e-7 for i in range(n - 1)])
    return d, e


def load_real_godunov(param):
    """Load real Godunov from STCollection."""
    path = f'/tmp/STCollection/DATA/T_Godunov_{param}.dat'
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found. Clone STCollection first.')
    with open(path) as f:
        n = int(f.readline().strip())
        d_tri = np.zeros(n)
        e_tri = np.zeros(max(n - 1, 0))
        for i in range(n):
            parts = f.readline().split()
            d_tri[i] = float(parts[1])
            if i < n - 1 and len(parts) > 2:
                e_tri[i] = float(parts[2])
    # Suppress side-effect prints from test_dense_to_bidiag import
    _out = sys.stdout
    sys.stdout = io.StringIO()
    from test_dense_to_bidiag import tridiag_to_bidiag
    sys.stdout = _out
    return tridiag_to_bidiag(d_tri, e_tri)


def run_variant(d, e, libpath):
    """Run bidiag_svd with a specific library, capturing trace from stderr."""
    import tempfile
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Write arrays to temp file to avoid code-string size issues
    tmpf = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
    json.dump({'d': d.tolist(), 'e': e.tolist()}, tmpf)
    tmpf.close()

    code = f'''
import sys, os, io, json, numpy as np
os.environ["XMR_LIB_PATH"] = "{libpath}"
sys.path.insert(0, ".")
sys.stdout = io.StringIO()
try:
    from test_dense_to_bidiag import tridiag_to_bidiag
except: pass
sys.stdout = sys.__stdout__
from mr3_gk import bidiag_svd
EPS = np.finfo(np.float64).eps
with open("{tmpf.name}") as f:
    data = json.load(f)
d = np.array(data["d"])
e = np.array(data["e"])
sigma, U, V, info = bidiag_svd(d, e)
n = len(d)
B = np.diag(d) + np.diag(e, 1)
bnorm = max(np.max(np.abs(d) + np.append(np.abs(e), 0)), 1e-300)
ou = np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS)
ov = np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS)
res = np.max(np.abs(B - U @ np.diag(sigma) @ V.T)) / (bnorm * n * EPS)
print(json.dumps({{"ortU": ou, "ortV": ov, "res": res, "sigma": sorted(sigma.tolist())}}))
'''
    proc = subprocess.run(
        ['python3', '-c', code],
        capture_output=True, text=True, timeout=300, cwd=base
    )
    os.unlink(tmpf.name)
    metrics = json.loads(proc.stdout.strip().split('\n')[-1])
    trace = parse_trace(proc.stderr)
    return metrics, trace


def parse_trace(stderr_text):
    """Parse dbg_ trace lines into structured data."""
    nodes = []
    singletons = []
    clusters = []
    cluster_fails = []
    singleton_rqis = []

    for raw_line in stderr_text.split('\n'):
        line = raw_line.strip()
        if '[dlaxr' not in line:
            continue
        # Parse key=value pairs from the line
        kv = {}
        for token in line.split():
            if '=' in token and not token.startswith('['):
                k, v = token.split('=', 1)
                kv[k] = v
        try:
            if 'IL' in kv and 'SINGLETON' not in line and 'CLUSTER' not in line:
                nodes.append({
                    'depth': int(kv['depth']),
                    'il': int(kv['IL']),
                    'iu': int(kv['IU']),
                    'taubar': float(kv['taubar'])
                })
            elif 'SINGLETON depth' in line and 'RQI' not in line:
                singletons.append({
                    'depth': int(kv['depth']),
                    'i': int(kv['i']),
                    'lambda': float(kv.get('lambda', '0')),
                    'taubar': float(kv['taubar']),
                    'w': float(kv['w'])
                })
            elif 'CLUSTER SHIFT' in line:
                clusters.append({
                    'depth': int(kv['depth']),
                    'i': int(kv['i']),
                    'j': int(kv['j']),
                    'tau': float(kv['tau']),
                    'new_taubar': float(kv['new_taubar'])
                })
            elif 'CLUSTER FAIL' in line:
                cluster_fails.append({
                    'depth': int(kv['depth']),
                    'i': int(kv['i']),
                    'j': int(kv['j']),
                    'status': int(kv.get('xrfsta', '0'))
                })
            elif 'SINGLETON_RQI' in line:
                singleton_rqis.append({
                    'idx': int(kv['idx']),
                    'nrqi': int(kv['nrqi']),
                    'nbis': int(kv['nbis']),
                    'resid': float(kv['resid'])
                })
        except (KeyError, ValueError):
            continue

    # Filter to first XMR call only.
    # Each dlaxrv call starts with a depth=0 node. Take everything from
    # the first depth=0 node up to (but not including) the second.
    first_d0 = next((i for i, nd in enumerate(nodes) if nd['depth'] == 0), 0)
    second_d0 = next((i for i, nd in enumerate(nodes) if nd['depth'] == 0 and i > first_d0), len(nodes))
    nodes = nodes[first_d0:second_d0]

    # For singletons/clusters/rqis, use order: take only those that appeared
    # between the first and second depth-0 node markers in the raw trace.
    # Since we parsed in order, count how many total events precede the second call.
    # Simpler: just take the first N_pos singletons where N_pos = IU - IL + 1
    if nodes:
        n_evals = nodes[0]['iu'] - nodes[0]['il'] + 1
        singletons = singletons[:n_evals]
        singleton_rqis = singleton_rqis[:n_evals]
        # Clusters: take only those within the depth range of first call
        max_depth = max(nd['depth'] for nd in nodes) if nodes else 0
        n_clusters = sum(1 for c in clusters if c['depth'] <= max_depth)
        clusters = clusters[:n_clusters]
        cluster_fails = cluster_fails[:len(cluster_fails)]  # keep all (usually few)

    return {
        'nodes': nodes,
        'singletons': singletons,
        'clusters': clusters,
        'cluster_fails': cluster_fails,
        'singleton_rqis': singleton_rqis,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=20)
    parser.add_argument('--real', type=str, default=None)
    args = parser.parse_args()

    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if args.real:
        print(f'Loading real Godunov {args.real}...')
        d, e = load_real_godunov(args.real)
    else:
        print(f'Creating synthetic Godunov n={args.n}...')
        d, e = make_synthetic_godunov(args.n)

    n = len(d)
    print(f'  n={n}, d range=[{np.min(np.abs(d)):.3e}, {np.max(np.abs(d)):.3e}]')
    print(f'  e range=[{np.min(np.abs(e)):.3e}, {np.max(np.abs(e)):.3e}]')

    # True SVD for reference
    B = np.diag(d) + np.diag(e, 1)
    _, s_true, _ = np.linalg.svd(B)
    s_true = np.sort(s_true)

    # Check libraries exist
    lib_fix = os.path.join(base, 'libxmr_traced.so')
    lib_noadd = os.path.join(base, 'libxmr_traced_noadd.so')
    for lib in [lib_fix, lib_noadd]:
        if not os.path.exists(lib):
            print(f'ERROR: {lib} not found. Run: bash build_traced_tree.sh')
            sys.exit(1)

    print('\nRunning WITH additive fix...')
    metrics_fix, trace_fix = run_variant(d, e, lib_fix)
    print(f'  ortU={metrics_fix["ortU"]:.3f}  res={metrics_fix["res"]:.3f}')
    print(f'  nodes={len(trace_fix["nodes"])} singletons={len(trace_fix["singletons"])} clusters={len(trace_fix["clusters"])}')

    print('\nRunning WITHOUT additive fix...')
    metrics_noadd, trace_noadd = run_variant(d, e, lib_noadd)
    print(f'  ortU={metrics_noadd["ortU"]:.3f}  res={metrics_noadd["res"]:.3f}')
    print(f'  nodes={len(trace_noadd["nodes"])} singletons={len(trace_noadd["singletons"])} clusters={len(trace_noadd["clusters"])}')

    # Build data for HTML
    data = {
        'matrix': {
            'n': n,
            'type': f'real_godunov_{args.real}' if args.real else f'synthetic_godunov_n{n}',
            'd': d.tolist(),
            'e': e.tolist(),
        },
        'true_sigma': s_true.tolist(),
        'variants': {
            'additive_fix': {'metrics': metrics_fix, 'trace': trace_fix},
            'no_fix': {'metrics': metrics_noadd, 'trace': trace_noadd},
        }
    }

    outpath = os.path.join(base, 'docs', 'godunov_tree_data.json')
    with open(outpath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'\nSaved data to {outpath} ({os.path.getsize(outpath)} bytes)')


if __name__ == '__main__':
    main()
