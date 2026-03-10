#pragma once
// Bidiagonal SVD via DBDSGR (Großer-Lang hgbsvd, coupling-based O(n²) approach)

#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

extern "C" {
    void dbdsgr_(const char* jobz, const char* range,
                 const int* n, double* d, double* e,
                 const double* vl, const double* vu,
                 const int* il, const int* iu,
                 const double* abstol, int* m, double* w,
                 double* z, const int* ldz, int* isuppz,
                 double* work, const int* lwork,
                 int* iwork, const int* liwork, int* info,
                 int jobz_len, int range_len);  // f2c ftnlen args

    double dlamch_(const char* cmach, int len);
}

struct BidiagSVDResult {
    std::vector<double> sigma;  // singular values (descending)
    std::vector<double> U;      // left singular vectors (n x n, column-major)
    std::vector<double> V;      // right singular vectors (n x n, column-major)
    int info;
};

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

    int PAD = 64;

    // Copy d and e (dbdsgr_ may modify them)
    double* d_buf = (double*)calloc(n + PAD, sizeof(double));
    double* e_buf = (double*)calloc(n + PAD, sizeof(double));
    for (int i = 0; i < n; i++) d_buf[i] = d[i];
    for (int i = 0; i < n - 1; i++) e_buf[i] = e[i];

    // JOBZ='A' requires LDZ >= 2*N
    // V in rows 0..N-1, U in rows LDZ/2..LDZ/2+N-1 (0-indexed)
    int ldz = 2 * n;
    char jobz = 'A';
    char range = 'A';
    double vl = 0.0, vu = 0.0;
    int il = 1, iu = n;
    double abstol = 0.0;  // use default
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
        result.info = info_call;
        free(d_buf); free(e_buf);
        free(W_buf); free(Z_buf); free(isuppz_buf);
        free(work_buf); free(iwork_buf);
        return result;
    }

    // Extract results
    // W_buf contains m_found singular values in ascending order
    // Z_buf layout (column-major, ldz rows):
    //   rows 0..n-1: V (right singular vectors), columns 0..m_found-1
    //   rows ldz/2..ldz/2+n-1: U (left singular vectors), columns 0..m_found-1
    int half = ldz / 2;

    result.sigma.resize(n, 0.0);
    result.U.assign(n * n, 0.0);
    result.V.assign(n * n, 0.0);

    // Reverse order: dbdsgr returns descending, we keep descending
    // Storage convention: column-major with leading dimension n
    // V[col*n + row] = V_{row,col},  U[col*n + row] = U_{row,col}
    for (int i = 0; i < m_found && i < n; i++) {
        result.sigma[i] = W_buf[i];

        // V: rows 0..n-1 of Z column i (column-major: V[col*n+row])
        for (int r = 0; r < n; r++) {
            result.V[i * n + r] = Z_buf[i * ldz + r];
        }
        // U: rows half..half+n-1 of Z column i
        for (int r = 0; r < n; r++) {
            result.U[i * n + r] = Z_buf[i * ldz + half + r];
        }
    }

    free(d_buf); free(e_buf);
    free(W_buf); free(Z_buf); free(isuppz_buf);
    free(work_buf); free(iwork_buf);

    return result;
}
