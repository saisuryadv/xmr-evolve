#pragma once
// Bidiagonal SVD via DBDSQR (QR iteration baseline)

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>

struct BidiagSVDResult {
    std::vector<double> sigma;
    std::vector<double> U;
    std::vector<double> V;
    int info;
};

extern "C" {
    void dbdsqr_(const char* uplo, const int* n, const int* ncvt,
                 const int* nru, const int* ncc,
                 double* d, double* e,
                 double* vt, const int* ldvt,
                 double* u, const int* ldu,
                 double* c, const int* ldc,
                 double* work, int* info,
                 int uplo_len);

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

    result.sigma.resize(n);
    std::copy(d, d + n, result.sigma.data());
    std::vector<double> e_copy(n - 1);
    std::copy(e, e + n - 1, e_copy.data());

    result.U.resize(n * n, 0.0);
    for (int i = 0; i < n; i++) result.U[i * n + i] = 1.0;

    std::vector<double> VT(n * n, 0.0);
    for (int i = 0; i < n; i++) VT[i * n + i] = 1.0;

    std::vector<double> work(4 * n);
    int info;
    int ncvt = n, nru = n, ncc = 0;
    int ldvt = n, ldu = n, ldc = 1;

    dbdsqr_("U", &n, &ncvt, &nru, &ncc,
            result.sigma.data(), e_copy.data(),
            VT.data(), &ldvt,
            result.U.data(), &ldu,
            nullptr, &ldc,
            work.data(), &info,
            1);

    // Transpose VT -> V
    result.V.resize(n * n);
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            result.V[k * n + j] = VT[j * n + k];

    result.info = info;
    return result;
}
