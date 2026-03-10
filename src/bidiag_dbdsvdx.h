#pragma once
// Bidiagonal SVD via DBDSVDX (bisection+inverse iteration on TGK)
// Pure C LAPACK (f2c-generated)

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
    // f2c uses int (ftnlen) for hidden CHARACTER string lengths
    void dbdsvdx_(const char* uplo, const char* jobz, const char* range,
                  const int* n, double* d, double* e,
                  const double* vl, const double* vu,
                  const int* il, const int* iu,
                  int* ns, double* s, double* z, const int* ldz,
                  double* work, int* iwork, int* info,
                  int uplo_len, int jobz_len, int range_len);

    double dlamch_(const char* cmach, int len);
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

    int ldz = 2 * n;
    int ns;

    // All buffers via calloc with generous padding to accommodate
    // f2c's 1-based pointer adjustments (--d, --e, z -= offset, etc.)
    int PAD = 64;
    double* d_buf = (double*)calloc(n + PAD, sizeof(double));
    double* e_buf = (double*)calloc(n + PAD, sizeof(double));
    double* s_buf = (double*)calloc(n + PAD, sizeof(double));
    double* z_buf = (double*)calloc((long long)ldz * 2 * n + PAD * ldz, sizeof(double));
    double* work  = (double*)calloc(14 * n + PAD, sizeof(double));
    int*    iwork = (int*)calloc(12 * n + PAD, sizeof(int));

    std::copy(d, d + n, d_buf);
    if (n > 1) std::copy(e, e + n - 1, e_buf);

    double vl = 0.0, vu = 0.0;
    int il = 1, iu = n;
    int info;

    dbdsvdx_((char*)"U", (char*)"V", (char*)"A",
             &n, d_buf, e_buf,
             &vl, &vu, &il, &iu,
             &ns, s_buf, z_buf, &ldz,
             work, iwork, &info,
             1, 1, 1);

    if (info == 0) {
        // DBDSVDX returns singular values in ASCENDING order
        // Z rows 0..n-1 = U columns, rows n..2n-1 = V columns
        result.sigma.resize(ns);
        result.U.resize(n * n, 0.0);
        result.V.resize(n * n, 0.0);

        for (int j = 0; j < ns; j++) {
            // Reverse to descending order
            int src = ns - 1 - j;
            result.sigma[j] = s_buf[src];

            for (int i = 0; i < n; i++) {
                result.U[j * n + i] = z_buf[src * ldz + i];       // rows 0..n-1
                result.V[j * n + i] = z_buf[src * ldz + n + i];   // rows n..2n-1
            }
        }
    } else {
        result.info = info;
    }

    free(d_buf);
    free(e_buf);
    free(s_buf);
    free(z_buf);
    free(work);
    free(iwork);

    return result;
}
