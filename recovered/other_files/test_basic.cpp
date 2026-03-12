#include "bidiag_svd.h"
#include <cstdio>
#include <cmath>

int main() {
    printf("Starting basic test...\n");
    fflush(stdout);

    // Simple 3x3 bidiagonal: d=[3,2,1], e=[1,1]
    int n = 3;
    double d[] = {3.0, 2.0, 1.0};
    double e[] = {1.0, 1.0};

    printf("Calling bidiag_svd with n=%d\n", n);
    fflush(stdout);

    BidiagSVDResult result = bidiag_svd(n, d, e);

    printf("info = %d\n", result.info);
    printf("sigma:");
    for (int i = 0; i < (int)result.sigma.size(); i++)
        printf(" %.6f", result.sigma[i]);
    printf("\n");

    return 0;
}
