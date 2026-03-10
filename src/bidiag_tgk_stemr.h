#pragma once
// Bidiagonal SVD via TGK + DSTEMR (LAPACK MR³)

#include "bidiag_tgk_common.h"

extern "C" {
    void dstemr_(const char* jobz, const char* range,
                 const int* n, double* d, double* e,
                 const double* vl, const double* vu,
                 const int* il, const int* iu,
                 int* m, double* w, double* z, const int* ldz,
                 const int* nzc, int* isuppz,
                 int* tryrac, double* work, const int* lwork,
                 int* iwork, const int* liwork, int* info,
                 int jobz_len, int range_len);
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
            &tryrac, &work_query, &lwork_q, &iwork_query, &liwork_q, &info_ws,
            1, 1);

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
            &tryrac, work.data(), &lwork, iwork.data(), &liwork, &info_call,
            1, 1);

    if (info_call != 0) {
        result.info = info_call;
        return result;
    }

    extract_svd_from_tgk(n, W.data(), Z.data(), ldz, result);
    postprocess_tgk_svd(n, d, e, result);

    return result;
}
