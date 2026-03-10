#pragma once
// Hybrid bidiagonal SVD: HGBSVD (coupling-based O(n²)) with TGK+DSTEMR fallback
// Self-contained — no external header dependencies beyond standard library + LAPACK/BLAS

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>

// === Section 1: Result struct ===

struct BidiagSVDResult {
    std::vector<double> sigma;  // singular values (descending)
    std::vector<double> U;      // left singular vectors (n x n, column-major)
    std::vector<double> V;      // right singular vectors (n x n, column-major)
    int info;
};

// === Section 2: extern "C" declarations ===

extern "C" {
    // HGBSVD coupling-based SVD (Großer-Lang)
    // f2c-generated: requires ftnlen args for character*1 params
    void dbdsgr_(const char* jobz, const char* range,
                 const int* n, double* d, double* e,
                 const double* vl, const double* vu,
                 const int* il, const int* iu,
                 const double* abstol, int* m, double* w,
                 double* z, const int* ldz, int* isuppz,
                 double* work, const int* lwork,
                 int* iwork, const int* liwork, int* info,
                 int jobz_len, int range_len);

    // LAPACK MR³ eigensolver
    void dstemr_(const char* jobz, const char* range,
                 const int* n, double* d, double* e,
                 const double* vl, const double* vu,
                 const int* il, const int* iu,
                 int* m, double* w, double* z, const int* ldz,
                 const int* nzc, int* isuppz,
                 int* tryrac, double* work, const int* lwork,
                 int* iwork, const int* liwork, int* info);

    // Machine constants and BLAS
    double dlamch_(const char* cmach, int len);
    double dnrm2_(const int* n, const double* x, const int* incx);
    double ddot_(const int* n, const double* x, const int* incx,
                 const double* y, const int* incy);
    void dscal_(const int* n, const double* a, double* x, const int* incx);
    void daxpy_(const int* n, const double* a, const double* x, const int* incx,
                double* y, const int* incy);
}

// === Section 3: TGK construction ===

// Build TGK tridiagonal: diagonal = 0, off-diag = d[0], e[0], d[1], e[1], ..., d[n-1]
inline void build_tgk(int n, const double* d, const double* e,
                       std::vector<double>& tgk_d, std::vector<double>& tgk_e) {
    int m2n = 2 * n;
    tgk_d.assign(m2n, 0.0);
    tgk_e.assign(m2n - 1, 0.0);

    for (int i = 0; i < n; i++) {
        tgk_e[2 * i] = d[i];
        if (i < n - 1)
            tgk_e[2 * i + 1] = e[i];
    }
}

// === Section 4: U/V extraction from TGK eigenvectors ===

inline void extract_svd_from_tgk(int n, const double* W, const double* Z, int ldz,
                                  BidiagSVDResult& result) {
    int m2n = 2 * n;
    result.sigma.resize(n);
    result.U.resize(n * n, 0.0);
    result.V.resize(n * n, 0.0);

    for (int j = 0; j < n; j++) {
        int eig_idx = m2n - 1 - j;  // largest positive eigenvalue first
        result.sigma[j] = std::abs(W[eig_idx]);

        const double* evec = &Z[eig_idx * (long long)ldz];

        // V from even rows (0,2,4,...), U from odd rows (1,3,5,...)
        for (int i = 0; i < n; i++) {
            result.V[j * n + i] = evec[2 * i];
            result.U[j * n + i] = evec[2 * i + 1];
        }
    }
}

// === Section 5: Post-processing (normalize, sign, onesided, reortho) ===

inline void postprocess_tgk_svd(int n, const double* d, const double* e,
                                 BidiagSVDResult& result) {
    const double EPS = dlamch_("E", 1);
    const double SAFMIN = dlamch_("S", 1);
    int one = 1;

    // Normalize U and V columns
    for (int j = 0; j < n; j++) {
        double* v_col = &result.V[j * n];
        double* u_col = &result.U[j * n];

        double v_norm = dnrm2_(&n, v_col, &one);
        double u_norm = dnrm2_(&n, u_col, &one);

        if (v_norm > SAFMIN) {
            double inv = 1.0 / v_norm;
            dscal_(&n, &inv, v_col, &one);
        }
        if (u_norm > SAFMIN) {
            double inv = 1.0 / u_norm;
            dscal_(&n, &inv, u_col, &one);
        }
    }

    // Sign consistency: ensure B*v ≈ σ*u (not -σ*u)
    for (int j = 0; j < n; j++) {
        if (result.sigma[j] < SAFMIN) continue;

        double* v_col = &result.V[j * n];
        double* u_col = &result.U[j * n];

        std::vector<double> bv(n, 0.0);
        for (int i = 0; i < n; i++) {
            bv[i] = d[i] * v_col[i];
            if (i < n - 1) bv[i] += e[i] * v_col[i + 1];
        }

        double dot = ddot_(&n, bv.data(), &one, u_col, &one);
        if (dot < 0) {
            double neg1 = -1.0;
            dscal_(&n, &neg1, u_col, &one);
        }
    }

    // One-sided recovery for small singular values
    double sigma_max = result.sigma.empty() ? 0.0 : result.sigma[0];
    double onesided_thresh = std::sqrt(n * EPS) * sigma_max;

    for (int j = 0; j < n; j++) {
        if (result.sigma[j] < onesided_thresh && result.sigma[j] > SAFMIN) {
            double* v_col = &result.V[j * n];
            double* u_col = &result.U[j * n];

            for (int i = 0; i < n; i++) {
                u_col[i] = d[i] * v_col[i];
                if (i < n - 1) u_col[i] += e[i] * v_col[i + 1];
                u_col[i] /= result.sigma[j];
            }

            double u_norm = dnrm2_(&n, u_col, &one);
            if (u_norm > SAFMIN) {
                double inv = 1.0 / u_norm;
                dscal_(&n, &inv, u_col, &one);
            }
        }
    }

    // Selective reorthogonalization (chunked MGS, O(n²))
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

// === Section 6: HGBSVD path ===

inline bool try_hgbsvd(int n, const double* d, const double* e, BidiagSVDResult& result) {
    int PAD = 64;

    // Copy d and e (dbdsgr_ may modify them)
    double* d_buf = (double*)calloc(n + PAD, sizeof(double));
    double* e_buf = (double*)calloc(n + PAD, sizeof(double));
    for (int i = 0; i < n; i++) d_buf[i] = d[i];
    for (int i = 0; i < n - 1; i++) e_buf[i] = e[i];

    // JOBZ='A' requires LDZ >= 2*N
    int ldz = 2 * n;
    char jobz = 'A';
    char range = 'A';
    double vl = 0.0, vu = 0.0;
    int il = 1, iu = n;
    double abstol = 0.0;
    int m_found;

    // Workspace query
    int lwork = -1, liwork = -1;
    double work_query;
    int iwork_query;
    int info_ws;
    double* ws_w = (double*)calloc(n + PAD, sizeof(double));
    int* ws_isuppz = (int*)calloc(2 * ldz + PAD, sizeof(int));

    dbdsgr_(&jobz, &range, &n, d_buf, e_buf,
            &vl, &vu, &il, &iu, &abstol,
            &m_found, ws_w, nullptr, &ldz, ws_isuppz,
            &work_query, &lwork, &iwork_query, &liwork, &info_ws, 1, 1);

    lwork = (int)work_query + 1;
    liwork = iwork_query + 1;
    if (lwork < 20 * n) lwork = 20 * n;
    if (liwork < 10 * n) liwork = 10 * n;

    // Re-copy d and e (workspace query may have modified them)
    for (int i = 0; i < n; i++) d_buf[i] = d[i];
    for (int i = 0; i < n - 1; i++) e_buf[i] = e[i];

    // Allocate
    double* W_buf = (double*)calloc(n + PAD, sizeof(double));
    double* Z_buf = (double*)calloc((long long)ldz * n + PAD * (long long)ldz, sizeof(double));
    int* isuppz_buf = (int*)calloc(2 * ldz + PAD, sizeof(int));
    double* work_buf = (double*)calloc(lwork + PAD, sizeof(double));
    int* iwork_buf = (int*)calloc(liwork + PAD, sizeof(int));
    int info_call;

    dbdsgr_(&jobz, &range, &n, d_buf, e_buf,
            &vl, &vu, &il, &iu, &abstol,
            &m_found, W_buf, Z_buf, &ldz, isuppz_buf,
            work_buf, &lwork, iwork_buf, &liwork, &info_call, 1, 1);

    free(ws_w);
    free(ws_isuppz);

    if (info_call != 0) {
        free(d_buf); free(e_buf);
        free(W_buf); free(Z_buf); free(isuppz_buf);
        free(work_buf); free(iwork_buf);
        return false;
    }

    // Extract results
    // W_buf: m_found singular values
    // Z_buf layout (column-major, ldz rows):
    //   rows 0..n-1: V (right singular vectors)
    //   rows ldz/2..ldz/2+n-1: U (left singular vectors)
    int half = ldz / 2;

    result.sigma.resize(n, 0.0);
    result.U.assign(n * n, 0.0);
    result.V.assign(n * n, 0.0);

    for (int i = 0; i < m_found && i < n; i++) {
        result.sigma[i] = W_buf[i];

        for (int r = 0; r < n; r++) {
            result.V[i * n + r] = Z_buf[i * ldz + r];
        }
        for (int r = 0; r < n; r++) {
            result.U[i * n + r] = Z_buf[i * ldz + half + r];
        }
    }

    free(d_buf); free(e_buf);
    free(W_buf); free(Z_buf); free(isuppz_buf);
    free(work_buf); free(iwork_buf);
    return true;
}

// === Section 7: TGK+DSTEMR path ===

inline bool try_tgk_dstemr(int n, const double* d, const double* e, BidiagSVDResult& result) {
    int m2n = 2 * n;

    // Build TGK
    std::vector<double> tgk_d, tgk_e;
    build_tgk(n, d, e, tgk_d, tgk_e);

    // Workspace query
    int m_found;
    double vl = 0.0, vu = 0.0;
    int il = 1, iu = m2n;
    int tryrac = 0;
    int info_ws;
    double work_query;
    int iwork_query;
    int lwork_q = -1, liwork_q = -1;
    int ldz = m2n, nzc = m2n;

    dstemr_("V", "A", &m2n, tgk_d.data(), tgk_e.data(),
            &vl, &vu, &il, &iu,
            &m_found, nullptr, nullptr, &ldz, &nzc, nullptr,
            &tryrac, &work_query, &lwork_q, &iwork_query, &liwork_q, &info_ws);

    int lwork = (int)work_query + 1;
    int liwork = iwork_query + 1;
    if (lwork < 18 * m2n) lwork = 18 * m2n;
    if (liwork < 10 * m2n) liwork = 10 * m2n;

    // Re-build TGK (dstemr may have modified arrays)
    build_tgk(n, d, e, tgk_d, tgk_e);

    std::vector<double> W(m2n);
    std::vector<double> Z((long long)m2n * m2n);
    std::vector<int> isuppz(2 * m2n);
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    tryrac = 0;
    int info_call;

    dstemr_("V", "A", &m2n, tgk_d.data(), tgk_e.data(),
            &vl, &vu, &il, &iu,
            &m_found, W.data(), Z.data(), &ldz, &nzc, isuppz.data(),
            &tryrac, work.data(), &lwork, iwork.data(), &liwork, &info_call);

    if (info_call != 0) {
        return false;
    }

    extract_svd_from_tgk(n, W.data(), Z.data(), ldz, result);
    postprocess_tgk_svd(n, d, e, result);
    return true;
}

// === Section 8: Main entry point ===

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

    // Try HGBSVD first (coupling-based O(n²), best accuracy)
    if (try_hgbsvd(n, d, e, result)) {
        return result;
    }

    // Fallback: TGK + DSTEMR with post-processing
    if (try_tgk_dstemr(n, d, e, result)) {
        return result;
    }

    // Both failed
    result.info = -999;
    return result;
}
