#pragma once
// Common TGK construction + U/V extraction + post-processing
// Used by both DSTEMR and DSTEXR variants

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
}

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

// Extract U, V, sigma from TGK eigenpairs
// eigenvalues W[0..2n-1] ascending, eigenvectors Z column-major m2n x m2n
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

// Post-process: normalize, sign consistency, one-sided recovery, selective reortho
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
