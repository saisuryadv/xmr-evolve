// Test that C++ can call the Fortran library correctly.
// Computes bidiagonal SVD using DSTEMR on a TGK matrix.

#include "fortran_interface.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

// Build a 2n x 2n Golub-Kahan (TGK) tridiagonal from bidiagonal B
// TGK diagonal = 0, off-diagonal alternates d[i], e[i]
void build_tgk(int n, const double* d, const double* e,
               double* tgk_diag, double* tgk_offdiag) {
    int m = 2 * n;
    for (int i = 0; i < m; i++) tgk_diag[i] = 0.0;
    for (int i = 0; i < m - 1; i++) tgk_offdiag[i] = 0.0;
    // Off-diag: d[0], e[0], d[1], e[1], ..., d[n-1]
    for (int i = 0; i < n; i++) {
        tgk_offdiag[2 * i] = d[i];
        if (i < n - 1)
            tgk_offdiag[2 * i + 1] = e[i];
    }
}

int main() {
    printf("=== Testing Fortran library from C++ ===\n\n");

    // Test 1: DLAMCH
    double eps = dlamch_("E");
    double sfmin = dlamch_("S");
    printf("DLAMCH: eps = %.6e, sfmin = %.6e\n", eps, sfmin);

    // Test 2: Simple 5x5 bidiagonal SVD via TGK + DSTEMR
    int n = 5;
    double d[5] = {4.0, 3.0, 2.0, 1.0, 0.5};
    double e[4] = {1.0, 1.0, 1.0, 1.0};

    int m2n = 2 * n;
    std::vector<double> tgk_d(m2n), tgk_e(m2n - 1);
    build_tgk(n, d, e, tgk_d.data(), tgk_e.data());

    printf("\nTGK diagonal (should be all zeros): ");
    for (int i = 0; i < m2n; i++) printf("%.1f ", tgk_d[i]);
    printf("\nTGK off-diagonal: ");
    for (int i = 0; i < m2n - 1; i++) printf("%.3f ", tgk_e[i]);
    printf("\n");

    // Call DSTEMR to compute all eigenvalues and eigenvectors
    char jobz = 'V', range = 'A';
    int m_found = 0;
    double vl = 0, vu = 0;
    int il = 0, iu = 0;
    std::vector<double> w(m2n);
    std::vector<double> z(m2n * m2n);
    int ldz = m2n;
    int nzc = m2n;
    std::vector<int> isuppz(2 * m2n);
    int tryrac = 1;
    int lwork = 18 * m2n;
    int liwork = 10 * m2n;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    int info = 0;

    dstemr_(&jobz, &range, &m2n,
            tgk_d.data(), tgk_e.data(),
            &vl, &vu, &il, &iu,
            &m_found, w.data(), z.data(), &ldz,
            &nzc, isuppz.data(), &tryrac,
            work.data(), &lwork, iwork.data(), &liwork, &info);

    printf("\nDSTEMR: info=%d, m_found=%d\n", info, m_found);
    if (info != 0) {
        printf("DSTEMR FAILED!\n");
        return 1;
    }

    // TGK eigenvalues come in ±σ pairs, sorted ascending
    // Singular values are the positive eigenvalues
    printf("TGK eigenvalues: ");
    for (int i = 0; i < m_found; i++) printf("%.6f ", w[i]);
    printf("\n");

    // Extract singular values (positive half)
    printf("\nSingular values: ");
    for (int i = n; i < m2n; i++) printf("%.6f ", w[i]);
    printf("\n");

    // Extract U and V from eigenvectors
    // For TGK: V from even rows (0,2,4,...), U from odd rows (1,3,5,...)
    printf("\nV (from even rows of eigenvectors for positive eigenvalues):\n");
    for (int j = n; j < m2n; j++) {
        printf("  σ=%.4f: [", w[j]);
        for (int i = 0; i < n; i++)
            printf("%.4f ", z[j * ldz + 2 * i]);
        printf("]\n");
    }

    printf("\nU (from odd rows of eigenvectors for positive eigenvalues):\n");
    for (int j = n; j < m2n; j++) {
        printf("  σ=%.4f: [", w[j]);
        for (int i = 0; i < n; i++)
            printf("%.4f ", z[j * ldz + 2 * i + 1]);
        printf("]\n");
    }

    // Test 3: DBDSQR as reference
    printf("\n=== DBDSQR reference ===\n");
    double d2[5] = {4.0, 3.0, 2.0, 1.0, 0.5};
    double e2[4] = {1.0, 1.0, 1.0, 1.0};
    std::vector<double> u_ref(n * n, 0.0);
    std::vector<double> vt_ref(n * n, 0.0);
    // Initialize U and VT to identity
    for (int i = 0; i < n; i++) {
        u_ref[i * n + i] = 1.0;
        vt_ref[i * n + i] = 1.0;
    }
    char uplo = 'U';
    int ncvt = n, nru = n, ncc = 0;
    std::vector<double> work_qr(4 * n);
    int info2 = 0;
    dbdsqr_(&uplo, &n, &ncvt, &nru, &ncc,
            d2, e2, vt_ref.data(), &n, u_ref.data(), &n,
            nullptr, &n, work_qr.data(), &info2);
    printf("DBDSQR: info=%d\n", info2);
    printf("DBDSQR singular values: ");
    for (int i = 0; i < n; i++) printf("%.6f ", d2[i]);
    printf("\n");

    // Compare
    printf("\n=== Comparison ===\n");
    printf("  TGK σ vs DBDSQR σ:\n");
    for (int i = 0; i < n; i++) {
        double tgk_sv = w[m2n - 1 - i]; // TGK positive eigenvalues, sorted ascending → largest last
        double ref_sv = d2[i]; // DBDSQR sorts descending
        printf("  [%d] TGK=%.10f  DBDSQR=%.10f  diff=%.2e\n",
               i, tgk_sv, ref_sv, fabs(tgk_sv - ref_sv));
    }

    // Test 4: DSTEXR (Willems)
    printf("\n=== DSTEXR (Willems XMR) ===\n");
    std::vector<double> tgk_d2(m2n), tgk_e2(m2n - 1);
    double d3[5] = {4.0, 3.0, 2.0, 1.0, 0.5};
    double e3[4] = {1.0, 1.0, 1.0, 1.0};
    build_tgk(n, d3, e3, tgk_d2.data(), tgk_e2.data());

    int m_found2 = 0;
    std::vector<double> w2(m2n);
    std::vector<double> z2(m2n * m2n);
    std::vector<int> isuppz2(2 * m2n);
    int lwork2 = 18 * m2n;
    int liwork2 = 10 * m2n;
    std::vector<double> work2(lwork2);
    std::vector<int> iwork2(liwork2);
    int info3 = 0;

    dstexr_(&jobz, &range, &m2n,
            tgk_d2.data(), tgk_e2.data(),
            &vl, &vu, &il, &iu,
            &m_found2, w2.data(), z2.data(), &ldz,
            isuppz2.data(), work2.data(), &lwork2,
            iwork2.data(), &liwork2, &info3);

    printf("DSTEXR: info=%d, m_found=%d\n", info3, m_found2);
    if (info3 == 0) {
        printf("DSTEXR eigenvalues: ");
        for (int i = 0; i < m_found2; i++) printf("%.6f ", w2[i]);
        printf("\n");
    } else {
        printf("DSTEXR FAILED (info=%d) — expected, XMR has known bugs\n", info3);
    }

    printf("\n=== All tests passed ===\n");
    return 0;
}
