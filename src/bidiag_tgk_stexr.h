#pragma once
// Bidiagonal SVD via TGK + DSTEXR (Willems XMR improved MR³)

#include "bidiag_tgk_common.h"

extern "C" {
    // Willems XMR Fortran interface
    void dstexr_(int* n, double* D, double* E, int* wil, int* wiu,
                 double* W, double* Z, int* ldz, int* ISuppZ,
                 double* RWork, int* lrwork, int* IWork, int* liwork,
                 int* info);
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

    int m2n = 2 * n;
    int PAD = 64;

    // Build TGK matrix with padding for f2c pointer adjustments
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
    int idum = -1;
    double rdum = 0.0;
    // Use padded buffers for workspace query too
    double* ws_rdum = (double*)calloc(PAD, sizeof(double));
    int* ws_idum = (int*)calloc(PAD, sizeof(int));
    ws_idum[0] = -1;

    dstexr_(&m2n, ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, &lrwork, ws_idum, &liwork, &info_ws);

    free(ws_rdum);
    free(ws_idum);

    // Safety minimums
    if (lrwork < 18 * m2n) lrwork = 18 * m2n;
    if (liwork < 10 * m2n) liwork = 10 * m2n;

    int wlen = m2n;
    int ldz = m2n;

    // All buffers via calloc with padding to accommodate
    // f2c's 1-based pointer adjustments (--d, --e, z -= offset, etc.)
    double* W_buf = (double*)calloc(wlen + PAD, sizeof(double));
    double* Z_buf = (double*)calloc((long long)ldz * wlen + PAD * (long long)ldz, sizeof(double));
    int* ISuppZ_buf = (int*)calloc(2 * wlen + PAD, sizeof(int));
    double* RWork_buf = (double*)calloc(lrwork + PAD, sizeof(double));
    int* IWork_buf = (int*)calloc(liwork + PAD, sizeof(int));
    int info_call;

    dstexr_(&m2n, tgk_d, tgk_e, &wil, &wiu,
            W_buf, Z_buf, &ldz, ISuppZ_buf,
            RWork_buf, &lrwork, IWork_buf, &liwork, &info_call);

    if (info_call != 0) {
        result.info = info_call;
        free(tgk_d); free(tgk_e);
        free(W_buf); free(Z_buf); free(ISuppZ_buf);
        free(RWork_buf); free(IWork_buf);
        return result;
    }

    extract_svd_from_tgk(n, W_buf, Z_buf, ldz, result);
    postprocess_tgk_svd(n, d, e, result);

    free(tgk_d); free(tgk_e);
    free(W_buf); free(Z_buf); free(ISuppZ_buf);
    free(RWork_buf); free(IWork_buf);

    return result;
}
