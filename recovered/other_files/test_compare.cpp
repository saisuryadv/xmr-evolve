// Compare GK-form vs TYPE=0 on a few specific failing matrices
#include "bidiag_svd.h"
#include <cstdio>
#include <cmath>
#include <cstring>

double compute_ortho_v(int n, const BidiagSVDResult& r) {
    if (n < 2) return 0;
    int one = 1;
    double worst = 0;
    for (int j = 0; j < n; j++) {
        for (int k = j + 1; k < n; k++) {
            double dot = std::abs(ddot_(&n, &r.V[j*n], &one, &r.V[k*n], &one));
            if (dot > worst) worst = dot;
        }
    }
    return worst;
}

void test_matrix(const char* name, int n, const double* d, const double* e) {
    printf("%-30s n=%d\n", name, n);

    // GK-form (disable_gk=0)
    BidiagSVDResult r1 = call_dstexr(n, d, e, 0);
    double ortho1 = compute_ortho_v(n, r1);
    printf("  GK-form:  info=%d orthoV=%.4e sigma[0]=%.6f\n",
           r1.info, ortho1, r1.sigma.empty() ? 0.0 : r1.sigma[0]);

    // TYPE=0 (disable_gk=1)
    BidiagSVDResult r2 = call_dstexr(n, d, e, 1);
    double ortho2 = compute_ortho_v(n, r2);
    printf("  TYPE=0:   info=%d orthoV=%.4e sigma[0]=%.6f\n",
           r2.info, ortho2, r2.sigma.empty() ? 0.0 : r2.sigma[0]);

    printf("  Winner: %s\n\n", ortho1 < ortho2 ? "GK-form" : "TYPE=0");
}

int main() {
    // B_bug414: n=4
    {
        double d[] = {1.0, 2.0, 3.0, 4.0};
        double e[] = {0.5, 0.5, 0.5};
        test_matrix("simple_4x4", 4, d, e);
    }

    // stemr_killer-like: graded spectrum
    {
        int n = 10;
        double d[10], e[9];
        double eps = dlamch_("E", 1);
        double c = 1.0 / std::sqrt(eps);
        for (int i = 0; i < n; i++) d[i] = std::pow(c, -(double)(i)/(n-1));
        for (int i = 0; i < n-1; i++) e[i] = d[i] * 0.1;
        test_matrix("graded_10", n, d, e);
    }

    // two clusters
    {
        int n = 10;
        double d[10], e[9];
        for (int i = 0; i < 5; i++) d[i] = 1.0 + i * 1e-14;
        for (int i = 5; i < 10; i++) d[i] = 2.0 + (i-5) * 1e-14;
        for (int i = 0; i < 9; i++) e[i] = 0.001;
        test_matrix("two_clusters_10", n, d, e);
    }

    // Identity-like
    {
        int n = 5;
        double d[] = {1, 1, 1, 1, 1};
        double e[] = {1e-15, 1e-15, 1e-15, 1e-15};
        test_matrix("near_identity_5", n, d, e);
    }

    return 0;
}
