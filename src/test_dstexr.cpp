#include <cstdio>
#include <cmath>
#include <vector>
#include "clapack/f2c.h"

extern "C" {
    int dstexr_(integer *n, doublereal *d, doublereal *e,
                integer *wil, integer *wiu,
                doublereal *w, doublereal *z, integer *ldz,
                integer *isuppz,
                doublereal *rwork, integer *lrwork,
                integer *iwork, integer *liwork,
                integer *info);

    doublereal dlamch_(char *cmach, ftnlen cmach_len);
}

int main() {
    // Simple 3x3 tridiagonal: d=[2,2,2], e=[1,1]
    // Eigenvalues: 2-sqrt(2), 2, 2+sqrt(2)
    integer n = 3;
    doublereal d[3] = {2.0, 2.0, 2.0};
    doublereal e[2] = {1.0, 1.0};

    // Workspace query
    integer wil = 1, wiu = 3;
    integer lrwork = -1, liwork = -1;
    integer info;
    integer idum = -1;
    doublereal rdum = 0.0;

    dstexr_(&n, &rdum, &rdum, &idum, &idum,
            &rdum, &rdum, &idum, &idum,
            &rdum, &lrwork, &idum, &liwork, &info);

    printf("Workspace query: lrwork=%d, liwork=%d, info=%d\n", (int)lrwork, (int)liwork, (int)info);

    if (lrwork < 20*n) lrwork = 20*n;
    if (liwork < 12*n) liwork = 12*n;

    integer ldz = n;
    std::vector<doublereal> w(n);
    std::vector<doublereal> z(n * n);
    std::vector<integer> isuppz(2 * n);
    std::vector<doublereal> rwork(lrwork);
    std::vector<integer> iwork(liwork);

    dstexr_(&n, d, e, &wil, &wiu,
            w.data(), z.data(), &ldz, isuppz.data(),
            rwork.data(), &lrwork, iwork.data(), &liwork, &info);

    printf("dstexr returned info=%d\n", (int)info);

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("  w[%d] = %.15f\n", i, w[i]);
    }

    double expected[3] = {2.0 - sqrt(2.0), 2.0, 2.0 + sqrt(2.0)};
    printf("Expected:\n");
    for (int i = 0; i < n; i++) {
        printf("  expected[%d] = %.15f\n", i, expected[i]);
    }

    printf("\nEigenvectors (columns):\n");
    for (int j = 0; j < n; j++) {
        printf("  v[%d] = [", j);
        for (int i = 0; i < n; i++) {
            printf(" %.6f", z[j*n + i]);
        }
        printf(" ]\n");
    }

    // Check orthogonality
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double dot = 0;
            for (int k = 0; k < n; k++)
                dot += z[i*n + k] * z[j*n + k];
            if (i == j) {
                printf("  ||v[%d]|| = %.15f\n", i, sqrt(dot));
            } else {
                printf("  v[%d].v[%d] = %.2e\n", i, j, dot);
            }
        }
    }

    return 0;
}
