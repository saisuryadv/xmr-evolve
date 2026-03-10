#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "clapack/f2c.h"

extern "C" {
    int dbdsvdx_(char *uplo, char *jobz, char *range,
                 integer *n, doublereal *d, doublereal *e,
                 doublereal *vl, doublereal *vu,
                 integer *il, integer *iu,
                 integer *ns, doublereal *s, doublereal *z, integer *ldz,
                 doublereal *work, integer *iwork, integer *info,
                 ftnlen uplo_len, ftnlen jobz_len, ftnlen range_len);
}

int main() {
    // saw_tooth n=10
    integer n = 10;
    std::vector<doublereal> d(n), e(n-1);

    // Simple graded matrix
    for (int i = 0; i < n; i++) d[i] = (double)(i + 1);
    for (int i = 0; i < n-1; i++) e[i] = 0.5;

    integer ldz = 2 * n;
    integer ns;
    std::vector<doublereal> s(n);
    std::vector<doublereal> z((long long)ldz * 2 * n);
    std::vector<doublereal> work(14 * n);
    std::vector<integer> iwork(12 * n);
    doublereal vl = 0, vu = 0;
    integer il = 1, iu = n;
    integer info;

    printf("Calling dbdsvdx_...\n");
    dbdsvdx_((char*)"U", (char*)"V", (char*)"A",
             &n, d.data(), e.data(),
             &vl, &vu, &il, &iu,
             &ns, s.data(), z.data(), &ldz,
             work.data(), iwork.data(), &info,
             (ftnlen)1, (ftnlen)1, (ftnlen)1);

    printf("info=%d, ns=%d\n", (int)info, (int)ns);
    for (int i = 0; i < ns && i < 5; i++)
        printf("  s[%d] = %.10f\n", i, s[i]);

    // Try 20 calls to provoke the corruption
    for (int trial = 0; trial < 20; trial++) {
        for (int i = 0; i < n; i++) d[i] = (double)(i + 1 + trial * 0.1);
        for (int i = 0; i < n-1; i++) e[i] = 0.5 + trial * 0.01;

        dbdsvdx_((char*)"U", (char*)"V", (char*)"A",
                 &n, d.data(), e.data(),
                 &vl, &vu, &il, &iu,
                 &ns, s.data(), z.data(), &ldz,
                 work.data(), iwork.data(), &info,
                 (ftnlen)1, (ftnlen)1, (ftnlen)1);
        printf("Trial %d: info=%d, ns=%d\n", trial, (int)info, (int)ns);
    }

    printf("All done.\n");
    return 0;
}
