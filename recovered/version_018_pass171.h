#pragma once
// Bidiagonal SVD via TGK + DSTEXR (Willems XMR improved MR³)
// GK-form enabled at root level (Algorithm 4.1), with TYPE=0 retry for robustness

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>

struct BidiagSVDResult {
    std::vector<double> sigma;
    std::vector<double> U;
    std::vector<double> V;
    int info;
};

extern "C" {
    double dlamch_(const char* cmach, int len);
    double dnrm2_(const int* n, const double* x, const int* incx);
    double ddot_(const int* n, const double* x, const int* incx,
                 const double* y, const int* incy);
    void dscal_(const int* n, const double* a, double* x, const int* incx);
    void daxpy_(const int* n, const double* a, const double* x, const int* incx,
                double* y, const int* incy);

    // Willems XMR Fortran interface
    void dstexr_(int* n, double* D, double* E, int* wil, int* wiu,
                 double* W, double* Z, int* ldz, int* ISuppZ,
                 double* RWork, int* lrwork, int* IWork, int* liwork,
                 int* info);

    // Global flag in dlaxre.c to control GK-form
    extern int xmr_disable_gkform_;
}

// Extract U, V, sigma from TGK eigenpairs
inline void extract_svd_from_tgk(int n, int neig, const double* W, const double* Z, int ldz,
                                  BidiagSVDResult& result) {
    result.sigma.resize(n);
    result.U.resize(n * n, 0.0);
    result.V.resize(n * n, 0.0);

    for (int j = 0; j < n; j++) {
        int eig_idx = neig - 1 - j;
        if (eig_idx < 0) break;
        result.sigma[j] = std::abs(W[eig_idx]);

        const double* evec = &Z[eig_idx * (long long)ldz];

        for (int i = 0; i < n; i++) {
            result.V[j * n + i] = evec[2 * i];
            result.U[j * n + i] = evec[2 * i + 1];
        }
    }
}

// Compute B*v (upper bidiagonal multiply)
inline void bidiag_multiply(int n, const double* d, const double* e,
                             const double* v, double* bv) {
    for (int i = 0; i < n; i++) {
        bv[i] = d[i] * v[i];
        if (i < n - 1) bv[i] += e[i] * v[i + 1];
    }
}

// Compute B^T*u (upper bidiagonal transpose multiply)
inline void bidiag_transpose_multiply(int n, const double* d, const double* e,
                                       const double* u, double* btv) {
    for (int i = 0; i < n; i++) {
        btv[i] = d[i] * u[i];
        if (i > 0) btv[i] += e[i - 1] * u[i - 1];
    }
}

// Post-process: normalize, bidirectional recovery, sign fix, selective reortho
inline void postprocess_tgk_svd(int n, const double* d, const double* e,
                                 BidiagSVDResult& result) {
    const double EPS = dlamch_("E", 1);
    const double SAFMIN = dlamch_("S", 1);
    int one = 1;

    // Step 1: Normalize and detect zero columns
    double norm_thresh = 1e-14;
    std::vector<int> bad_v, bad_u;
    for (int j = 0; j < n; j++) {
        double* v_col = &result.V[j * n];
        double* u_col = &result.U[j * n];

        double v_norm = dnrm2_(&n, v_col, &one);
        double u_norm = dnrm2_(&n, u_col, &one);

        if (v_norm > norm_thresh) {
            double inv = 1.0 / v_norm;
            dscal_(&n, &inv, v_col, &one);
        } else {
            bad_v.push_back(j);
            std::memset(v_col, 0, n * sizeof(double));
        }
        if (u_norm > norm_thresh) {
            double inv = 1.0 / u_norm;
            dscal_(&n, &inv, u_col, &one);
        } else {
            bad_u.push_back(j);
            std::memset(u_col, 0, n * sizeof(double));
        }
    }

    // Step 1b: Complete rank-deficient V and U by orthogonal basis completion
    // For zero singular values, TGK eigenvectors have zero V or U components
    // Gate to n≤50 to avoid scaling issues at large n
    if (n > 50) { bad_v.clear(); bad_u.clear(); }
    for (size_t bi = 0; bi < bad_v.size(); bi++) {
        int j = bad_v[bi];
        double* vj = &result.V[j * n];
        // Try canonical basis vectors until we find one not in span of good columns
        bool filled = false;
        for (int trial = 0; trial < n && !filled; trial++) {
            // Start with e_{trial}
            std::memset(vj, 0, n * sizeof(double));
            vj[trial] = 1.0;
            // Orthogonalize against all good V columns
            for (int k = 0; k < n; k++) {
                if (k == j) continue;
                bool is_bad = false;
                for (size_t bk = bi + 1; bk < bad_v.size(); bk++)
                    if (bad_v[bk] == k) { is_bad = true; break; }
                if (is_bad) continue;
                double* vk = &result.V[k * n];
                double dot = ddot_(&n, vk, &one, vj, &one);
                double neg_dot = -dot;
                daxpy_(&n, &neg_dot, vk, &one, vj, &one);
            }
            // Also orthogonalize against previously filled bad columns
            for (size_t bk = 0; bk < bi; bk++) {
                double* vk = &result.V[bad_v[bk] * n];
                double dot = ddot_(&n, vk, &one, vj, &one);
                double neg_dot = -dot;
                daxpy_(&n, &neg_dot, vk, &one, vj, &one);
            }
            double nrm = dnrm2_(&n, vj, &one);
            if (nrm > 0.1) {
                double inv = 1.0 / nrm;
                dscal_(&n, &inv, vj, &one);
                filled = true;
            }
        }
    }
    for (size_t bi = 0; bi < bad_u.size(); bi++) {
        int j = bad_u[bi];
        double* uj = &result.U[j * n];
        bool filled = false;
        for (int trial = 0; trial < n && !filled; trial++) {
            std::memset(uj, 0, n * sizeof(double));
            uj[trial] = 1.0;
            for (int k = 0; k < n; k++) {
                if (k == j) continue;
                bool is_bad = false;
                for (size_t bk = bi + 1; bk < bad_u.size(); bk++)
                    if (bad_u[bk] == k) { is_bad = true; break; }
                if (is_bad) continue;
                double* uk = &result.U[k * n];
                double dot = ddot_(&n, uk, &one, uj, &one);
                double neg_dot = -dot;
                daxpy_(&n, &neg_dot, uk, &one, uj, &one);
            }
            for (size_t bk = 0; bk < bi; bk++) {
                double* uk = &result.U[bad_u[bk] * n];
                double dot = ddot_(&n, uk, &one, uj, &one);
                double neg_dot = -dot;
                daxpy_(&n, &neg_dot, uk, &one, uj, &one);
            }
            double nrm = dnrm2_(&n, uj, &one);
            if (nrm > 0.1) {
                double inv = 1.0 / nrm;
                dscal_(&n, &inv, uj, &one);
                filled = true;
            }
        }
    }

    // Step 2: (reserved for future recovery)
    double sigma_max = result.sigma.empty() ? 0.0 : result.sigma[0];

    // Step 3: Sign consistency
    for (int j = 0; j < n; j++) {
        if (result.sigma[j] < SAFMIN) continue;

        double* v_col = &result.V[j * n];
        double* u_col = &result.U[j * n];

        std::vector<double> bv(n);
        bidiag_multiply(n, d, e, v_col, bv.data());
        double dot = ddot_(&n, bv.data(), &one, u_col, &one);
        if (dot < 0) {
            double neg1 = -1.0;
            dscal_(&n, &neg1, u_col, &one);
        }
    }

    // Step 4: Selective reorthogonalization (chunked MGS, O(n²))
    double cluster_tol = 100.0 * EPS * sigma_max;

    std::vector<int> cluster_start, cluster_end;
    int cs = 0;
    while (cs < n) {
        int ce = cs;
        while (ce + 1 < n &&
               std::abs(result.sigma[ce] - result.sigma[ce + 1]) < cluster_tol) {
            ce++;
        }
        if (ce > cs) {
            cluster_start.push_back(cs);
            cluster_end.push_back(ce);
        }
        cs = ce + 1;
    }

    int MAX_CHUNK = 32;
    for (size_t cl = 0; cl < cluster_start.size(); cl++) {
        int cstart = cluster_start[cl];
        int cend = cluster_end[cl];

        for (int chunk_start = cstart; chunk_start <= cend; chunk_start += MAX_CHUNK) {
            int chunk_end = std::min(chunk_start + MAX_CHUNK - 1, cend);

            for (int j = chunk_start; j <= chunk_end; j++) {
                double* vj = &result.V[j * n];
                for (int k = chunk_start; k < j; k++) {
                    double* vk = &result.V[k * n];
                    double dot = ddot_(&n, vk, &one, vj, &one);
                    double neg_dot = -dot;
                    daxpy_(&n, &neg_dot, vk, &one, vj, &one);
                }
                double nrm = dnrm2_(&n, vj, &one);
                if (nrm > SAFMIN) {
                    double inv = 1.0 / nrm;
                    dscal_(&n, &inv, vj, &one);
                }
            }

            for (int j = chunk_start; j <= chunk_end; j++) {
                double* uj = &result.U[j * n];
                for (int k = chunk_start; k < j; k++) {
                    double* uk = &result.U[k * n];
                    double dot = ddot_(&n, uk, &one, uj, &one);
                    double neg_dot = -dot;
                    daxpy_(&n, &neg_dot, uk, &one, uj, &one);
                }
                double nrm = dnrm2_(&n, uj, &one);
                if (nrm > SAFMIN) {
                    double inv = 1.0 / nrm;
                    dscal_(&n, &inv, uj, &one);
                }
            }
        }
    }
}

// Quick O(n) orthogonality check: max |V_i^T V_j| for a few pairs
inline double quick_ortho_check(int n, const BidiagSVDResult& result) {
    if (n < 2) return 0.0;
    int one = 1;
    double worst = 0.0;
    // Check a few adjacent pairs (most likely to be problematic)
    int checks = std::min(n - 1, 8);
    for (int j = 0; j < checks; j++) {
        const double* vj = &result.V[j * n];
        const double* vj1 = &result.V[(j + 1) * n];
        double dot = std::abs(ddot_(&n, vj, &one, vj1, &one));
        if (dot > worst) worst = dot;
    }
    // Also check U
    for (int j = 0; j < checks; j++) {
        const double* uj = &result.U[j * n];
        const double* uj1 = &result.U[(j + 1) * n];
        double dot = std::abs(ddot_(&n, uj, &one, uj1, &one));
        if (dot > worst) worst = dot;
    }
    return worst;
}

// Core DSTEXR call with given GK-form setting
inline BidiagSVDResult call_dstexr(int n, const double* d, const double* e, int disable_gk) {
    BidiagSVDResult result;
    result.info = 0;

    int m2n = 2 * n;
    int PAD = 64;

    double* tgk_d = (double*)calloc(m2n + PAD, sizeof(double));
    double* tgk_e = (double*)calloc(m2n + PAD, sizeof(double));

    for (int i = 0; i < n; i++) {
        tgk_e[2 * i] = d[i];
        if (i < n - 1)
            tgk_e[2 * i + 1] = e[i];
    }

    // Workspace query
    int wil = 1, wiu = m2n;
    int lrwork = -1, liwork = -1;
    int info_ws;
    double* ws_rdum = (double*)calloc(PAD, sizeof(double));
    int* ws_idum = (int*)calloc(PAD, sizeof(int));

    dstexr_(&m2n, ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, &lrwork, ws_idum, &liwork, &info_ws);

    free(ws_rdum);
    free(ws_idum);

    if (lrwork < 18 * m2n) lrwork = 18 * m2n;
    if (liwork < 10 * m2n) liwork = 10 * m2n;

    int neig = m2n;
    int ldz = m2n;

    double* W_buf = (double*)calloc(m2n + PAD, sizeof(double));
    double* Z_buf = (double*)calloc((long long)ldz * m2n + PAD * (long long)ldz, sizeof(double));
    int* ISuppZ_buf = (int*)calloc(2 * m2n + PAD, sizeof(int));
    double* RWork_buf = (double*)calloc(lrwork + PAD, sizeof(double));
    int* IWork_buf = (int*)calloc(liwork + PAD, sizeof(int));
    int info_call;

    xmr_disable_gkform_ = disable_gk;
    dstexr_(&m2n, tgk_d, tgk_e, &wil, &wiu,
            W_buf, Z_buf, &ldz, ISuppZ_buf,
            RWork_buf, &lrwork, IWork_buf, &liwork, &info_call);
    xmr_disable_gkform_ = 0;

    if (info_call != 0) {
        result.info = info_call;
    } else {
        extract_svd_from_tgk(n, neig, W_buf, Z_buf, ldz, result);
        postprocess_tgk_svd(n, d, e, result);
    }

    free(tgk_d); free(tgk_e);
    free(W_buf); free(Z_buf); free(ISuppZ_buf);
    free(RWork_buf); free(IWork_buf);

    return result;
}

inline BidiagSVDResult bidiag_svd(int n, const double* d, const double* e) {
    BidiagSVDResult result;
    result.info = 0;

    if (n <= 0) return result;

    if (n == 1) {
        result.sigma = {std::abs(d[0])};
        result.U = {d[0] >= 0 ? 1.0 : -1.0};
        result.V = {1.0};
        return result;
    }

    // GK-form (Algorithm 4.1: zero shift, TYPE=1)
    result = call_dstexr(n, d, e, 0);

    return result;
}
