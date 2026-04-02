"""
Diagnostic mode for MR3-GK bidiagonal SVD.
Prints everything: input, splits, T_GK, shifts, eigenvalues, eigenvectors, singular vectors.

Usage: python3 diagnose.py <test_name> <n>
Example: python3 diagnose.py two_clusters 10
"""
import sys
import numpy as np
from full_eval import make
from mr3_gk import split_bidiag, bidiag_svd, EPS, SAFMIN

np.set_printoptions(precision=10, linewidth=200, suppress=False)


def diagnose(name, n):
    d_orig, e_orig = make(name, n)
    d = d_orig.copy()
    e = e_orig.copy()

    print("=" * 80)
    print(f"DIAGNOSTIC: {name} n={n}")
    print("=" * 80)

    # --- 1. Input ---
    print("\n--- INPUT BIDIAGONAL ---")
    print(f"d ({len(d)}): {d}")
    print(f"e ({len(e)}): {e}")
    print(f"max|d|={np.max(np.abs(d)):.6e}  max|e|={np.max(np.abs(e)):.6e}")

    # --- 2. Reference SVD ---
    B = np.diag(d) + np.diag(e, 1)
    s_ref = np.linalg.svd(B, compute_uv=False)
    s_ref_sorted = np.sort(s_ref)[::-1]
    print(f"\nReference singular values (descending): {s_ref_sorted}")

    # --- 3. Pre-scaling ---
    bnorm_scale = max(np.max(np.abs(d)), np.max(np.abs(e)) if n > 1 else 0.0)
    if bnorm_scale > 0.0:
        scale_factor = 1.0 / bnorm_scale
        d_scaled = d * scale_factor
        e_scaled = e * scale_factor
    else:
        scale_factor = 1.0
        d_scaled = d.copy()
        e_scaled = e.copy()
    print(f"\n--- PRE-SCALING ---")
    print(f"bnorm_scale={bnorm_scale:.6e}  scale_factor={scale_factor:.6e}")
    if scale_factor != 1.0:
        print(f"d_scaled: {d_scaled}")
        print(f"e_scaled: {e_scaled}")

    # --- 4. Bidiagonal splitting ---
    blocks = split_bidiag(d_scaled, e_scaled, n)
    print(f"\n--- BIDIAGONAL SPLITS ---")
    print(f"Number of blocks: {len(blocks)}")
    for i, (bbeg, bend) in enumerate(blocks):
        bk = bend - bbeg + 1
        kind = "singleton" if bk == 1 else "multi"
        print(f"  block[{i}]: [{bbeg},{bend}] size={bk} ({kind})")
        if bk == 1:
            print(f"    d[{bbeg}]={d_scaled[bbeg]:.10e}  sigma=|d|={abs(d_scaled[bbeg]):.10e}")
        else:
            print(f"    d[{bbeg}:{bend+1}] = {d_scaled[bbeg:bend+1]}")
            print(f"    e[{bbeg}:{bend}]   = {e_scaled[bbeg:bend]}")

    # --- 5. D1, D2 sign matrices ---
    d1 = np.ones(n)
    d2 = np.ones(n)
    d1[0] = 1.0
    d2[0] = np.sign(d_scaled[0]) if d_scaled[0] != 0 else 1.0
    for i in range(n - 1):
        s_e = np.sign(e_scaled[i]) if e_scaled[i] != 0 else 1.0
        d2[i + 1] = s_e / d1[i] if abs(d1[i]) > 0 else s_e
        s_d = np.sign(d_scaled[i + 1]) if d_scaled[i + 1] != 0 else 1.0
        d1[i + 1] = s_d / d2[i + 1] if abs(d2[i + 1]) > 0 else s_d
    print(f"\n--- SIGN MATRICES ---")
    print(f"D1: {d1}")
    print(f"D2: {d2}")

    # --- 6. For each multi-block, build T_GK and analyze ---
    from xmr_ctypes import xmr_eigenvectors, build_repr_from_tridiag
    from scipy.linalg.lapack import _flapack

    for bi, (bbeg, bend) in enumerate(blocks):
        bk = bend - bbeg + 1
        if bk <= 1:
            continue

        print(f"\n{'=' * 60}")
        print(f"MULTI-BLOCK [{bbeg},{bend}] bk={bk}")
        print(f"{'=' * 60}")

        m2k = 2 * bk
        bd_abs = np.abs(d_scaled[bbeg:bend + 1])
        be = np.abs(e_scaled[bbeg:bend])

        # Build T_GK
        bc = np.zeros(m2k)
        bte = np.zeros(m2k - 1)
        bte[0::2] = bd_abs
        if bk > 1:
            bte[1::2] = be[:bk - 1]

        print(f"\n--- T_GK CONSTRUCTION ---")
        print(f"T_GK size: {m2k}x{m2k}")
        print(f"diagonal (bc): all zeros")
        print(f"off-diagonal (bte): {bte}")

        # Build explicit T_GK matrix for reference
        T_GK = np.zeros((m2k, m2k))
        for i in range(m2k - 1):
            T_GK[i, i + 1] = bte[i]
            T_GK[i + 1, i] = bte[i]

        # Reference eigenvalues
        ref_evals = np.sort(np.linalg.eigvalsh(T_GK))
        print(f"\nReference T_GK eigenvalues (all {m2k}):")
        print(f"  {ref_evals}")

        n_nonneg = (m2k + 1) // 2
        ref_nonneg = ref_evals[m2k - n_nonneg:]
        print(f"Non-negative eigenvalues ({n_nonneg}): {ref_nonneg}")

        # Reference eigenvectors
        ref_evals_full, ref_evecs_full = np.linalg.eigh(T_GK)
        print(f"\nReference eigenvectors (non-negative half):")
        for j in range(n_nonneg):
            idx = m2k - n_nonneg + j
            ev = ref_evals_full[idx]
            z = ref_evecs_full[:, idx]
            even_norm = np.linalg.norm(z[0::2])
            odd_norm = np.linalg.norm(z[1::2])
            print(f"  eval={ev:+.10e}  |z_even|={even_norm:.6f}  |z_odd|={odd_norm:.6f}")

        # --- Bisection eigenvalues ---
        wil = m2k - n_nonneg + 1  # 1-based
        wiu = m2k
        m_bis, w_bis, _, _, info_bis = _flapack.dstebz(
            bc[:m2k].astype(np.float64), bte[:m2k - 1].astype(np.float64),
            2, 0.0, 0.0, int(wil), int(wiu), 0.0, b'E')
        bisect_evals = np.sort(w_bis[:m_bis])
        print(f"\nBisection eigenvalues (wil={wil} wiu={wiu}):")
        print(f"  {bisect_evals}")
        print(f"  vs reference: max|diff|={np.max(np.abs(bisect_evals - ref_nonneg)):.3e}")

        # All eigenvalues for spdiam
        m_all, w_all, _, _, _ = _flapack.dstebz(
            bc[:m2k].astype(np.float64), bte[:m2k - 1].astype(np.float64),
            2, 0.0, 0.0, 1, m2k, 0.0, b'E')
        all_evals = np.sort(w_all[:m_all])
        spdiam = all_evals[-1] - all_evals[0] if len(all_evals) > 0 else 0.0

        print(f"\nSpectral diameter: {spdiam:.10e}")

        # --- XMR eigenvectors ---
        print(f"\n--- XMR EIGENVECTOR COMPUTATION ---")
        root = build_repr_from_tridiag(bc, bte, m2k, k=m2k)
        w_xmr, Z_xmr, info_xmr = xmr_eigenvectors(
            root, bte, all_evals,
            wil=wil, wiu=wiu,
            spdiam=spdiam, gaptol=1e-3)

        print(f"XMR info: {info_xmr}")

        if info_xmr != 0:
            print("XMR FAILED, falling back to Python MR3")
            from mr3_gk import mr3_block, GAPTOL
            pos_evals = all_evals[wil - 1:]
            w_xmr, Z_xmr = mr3_block(bc, bte, m2k, pos_evals, spdiam)
            print(f"Python MR3 computed {len(w_xmr)} eigenvalues")

        print(f"\nComputed eigenvalues ({len(w_xmr)}):")
        print(f"  {w_xmr}")
        print(f"  vs reference: max|diff|={np.max(np.abs(np.sort(w_xmr) - ref_nonneg)):.3e}")

        # Eigenvector analysis
        print(f"\nEigenvector analysis:")
        for j in range(min(len(w_xmr), Z_xmr.shape[1])):
            z = Z_xmr[:m2k, j]
            even_norm = np.linalg.norm(z[0::2])
            odd_norm = np.linalg.norm(z[1::2])
            total_norm = np.linalg.norm(z)

            # Residual: T_GK * z - lambda * z
            Tz = T_GK @ z
            res = np.linalg.norm(Tz - w_xmr[j] * z)

            print(f"  [{j}] eval={w_xmr[j]:+.10e}  |z|={total_norm:.6f}  "
                  f"|z_even|(V)={even_norm:.6f}  |z_odd|(U)={odd_norm:.6f}  "
                  f"res={res:.3e}")

        # Orthogonality of eigenvectors
        ZZ = Z_xmr[:m2k, :len(w_xmr)].T @ Z_xmr[:m2k, :len(w_xmr)]
        orth_err = np.max(np.abs(ZZ - np.eye(len(w_xmr))))
        print(f"\nEigenvector orthogonality: max|Z^T Z - I| = {orth_err:.6e}")

        # Extract V (even rows) and U (odd rows) BEFORE D1/D2
        V_raw = Z_xmr[0::2, :len(w_xmr)]
        U_raw = Z_xmr[1::2, :len(w_xmr)]

        print(f"\nExtracted V (even rows) orthogonality: max|V^T V - I| = "
              f"{np.max(np.abs(V_raw.T @ V_raw - np.eye(len(w_xmr)) * 0.5)) / 0.5:.6e}")
        print(f"Extracted U (odd rows) orthogonality: max|U^T U - I| = "
              f"{np.max(np.abs(U_raw.T @ U_raw - np.eye(len(w_xmr)) * 0.5)) / 0.5:.6e}")

        # Per-pair orthogonality for V and U
        nev = len(w_xmr)
        if nev > 1:
            VV = V_raw.T @ V_raw
            UU = U_raw.T @ U_raw
            worst_v_i, worst_v_j = -1, -1
            worst_v_val = 0
            worst_u_i, worst_u_j = -1, -1
            worst_u_val = 0
            for ii in range(nev):
                for jj in range(ii + 1, nev):
                    v_dot = abs(VV[ii, jj])
                    u_dot = abs(UU[ii, jj])
                    if v_dot > worst_v_val:
                        worst_v_val = v_dot
                        worst_v_i, worst_v_j = ii, jj
                    if u_dot > worst_u_val:
                        worst_u_val = u_dot
                        worst_u_i, worst_u_j = ii, jj

            print(f"\nWorst V pair: ({worst_v_i},{worst_v_j}) |v_i^T v_j|={worst_v_val:.6e} "
                  f"evals=({w_xmr[worst_v_i]:.6e}, {w_xmr[worst_v_j]:.6e})")
            print(f"Worst U pair: ({worst_u_i},{worst_u_j}) |u_i^T u_j|={worst_u_val:.6e} "
                  f"evals=({w_xmr[worst_u_i]:.6e}, {w_xmr[worst_u_j]:.6e})")

    # --- 7. Run full bidiag_svd ---
    print(f"\n{'=' * 60}")
    print("FULL BIDIAG_SVD OUTPUT")
    print(f"{'=' * 60}")
    sigma, U, V, info = bidiag_svd(d_orig, e_orig)

    print(f"\nComputed singular values: {sigma}")
    print(f"Reference singular values: {s_ref_sorted}")
    print(f"max|sigma - sigma_ref|: {np.max(np.abs(np.sort(sigma)[::-1] - s_ref_sorted)):.3e}")

    # Quality metrics
    B_orig = np.diag(d_orig) + np.diag(e_orig, 1)
    bnorm = np.max(np.abs(B_orig))
    res_mat = B_orig - U @ np.diag(sigma) @ V.T
    res = np.max(np.abs(res_mat)) / (bnorm * n * EPS)
    ortU = np.max(np.abs(U.T @ U - np.eye(n))) / (n * EPS)
    ortV = np.max(np.abs(V.T @ V - np.eye(n))) / (n * EPS)

    print(f"\nres={res:.3f}  ortU={ortU:.3f}  ortV={ortV:.3f}")
    passed = res <= 5 and ortU <= 5 and ortV <= 5
    print(f"PASS: {passed}")

    # Per-column quality
    print(f"\nPer-column analysis:")
    for j in range(n):
        Bvj = B_orig @ V[:, j]
        res_j = np.linalg.norm(Bvj - sigma[j] * U[:, j])
        u_norm = np.linalg.norm(U[:, j])
        v_norm = np.linalg.norm(V[:, j])
        print(f"  [{j}] sigma={sigma[j]:.6e}  |Bv-su|={res_j:.3e}  "
              f"|u|={u_norm:.6f}  |v|={v_norm:.6f}")

    # Worst orthogonality pairs
    if n > 1:
        UU = U.T @ U - np.eye(n)
        VV = V.T @ V - np.eye(n)
        for label, mat in [("U", UU), ("V", VV)]:
            worst_i, worst_j = np.unravel_index(np.argmax(np.abs(mat)), mat.shape)
            if worst_i != worst_j:
                print(f"\n  Worst {label} pair: ({worst_i},{worst_j}) "
                      f"|{label}_i^T {label}_j|={abs(mat[worst_i,worst_j]):.6e} "
                      f"= {abs(mat[worst_i,worst_j])/(n*EPS):.1f}*n*eps  "
                      f"sigma=({sigma[worst_i]:.6e}, {sigma[worst_j]:.6e})")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 diagnose.py <test_name> <n>")
        sys.exit(1)
    name = sys.argv[1]
    n = int(sys.argv[2])
    diagnose(name, n)
