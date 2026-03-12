#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" {
    double dlamch_(const char* cmach, int len);
    void dstexr_(int* n, double* D, double* E, int* wil, int* wiu,
                 double* W, double* Z, int* ldz, int* ISuppZ,
                 double* RWork, int* lrwork, int* IWork, int* liwork,
                 int* info);
}

int main() {
    int n = 3;
    double d[] = {3.0, 2.0, 1.0};
    double e[] = {1.0, 1.0};

    int m2n = 2 * n;
    int PAD = 64;

    double* tgk_d = (double*)calloc(m2n + PAD, sizeof(double));
    double* tgk_e = (double*)calloc(m2n + PAD, sizeof(double));

    for (int i = 0; i < n; i++) {
        tgk_e[2 * i] = d[i];
        if (i < n - 1)
            tgk_e[2 * i + 1] = e[i];
    }
    printf("TGK diagonal (should be zero): ");
    for (int i = 0; i < m2n; i++) printf("%.2f ", tgk_d[i]);
    printf("\nTGK off-diagonal: ");
    for (int i = 0; i < m2n - 1; i++) printf("%.2f ", tgk_e[i]);
    printf("\n\n");

    // Workspace query
    int lrwork = -1, liwork = -1;
    int info_ws;
    double* ws_rdum = (double*)calloc(PAD, sizeof(double));
    int* ws_idum = (int*)calloc(PAD, sizeof(int));
    dstexr_(&m2n, ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, ws_rdum, ws_idum, ws_idum,
            ws_rdum, &lrwork, ws_idum, &liwork, &info_ws);
    printf("Workspace query: lrwork=%d liwork=%d info=%d\n", lrwork, liwork, info_ws);
    if (lrwork < 18 * m2n) lrwork = 18 * m2n;
    if (liwork < 10 * m2n) liwork = 10 * m2n;

    int ldz = m2n;

    // Test 1: ALL eigenvalues
    {
        int wil = 1, wiu = m2n;
        double* W = (double*)calloc(m2n + PAD, sizeof(double));
        double* Z = (double*)calloc((long long)ldz * m2n + PAD * (long long)ldz, sizeof(double));
        int* ISuppZ = (int*)calloc(2 * m2n + PAD, sizeof(int));
        double* RWork = (double*)calloc(lrwork + PAD, sizeof(double));
        int* IWork = (int*)calloc(liwork + PAD, sizeof(int));
        int info_call;

        dstexr_(&m2n, tgk_d, tgk_e, &wil, &wiu,
                W, Z, &ldz, ISuppZ,
                RWork, &lrwork, IWork, &liwork, &info_call);

        printf("ALL eigenvalues (wil=%d, wiu=%d): info=%d\n", wil, wiu, info_call);
        printf("  W: ");
        for (int i = 0; i < m2n; i++) printf("%.6f ", W[i]);
        printf("\n");

        free(W); free(Z); free(ISuppZ); free(RWork); free(IWork);
    }

    // Test 2: Upper half only
    {
        int wil = n + 1, wiu = m2n;
        double* W = (double*)calloc(m2n + PAD, sizeof(double));
        double* Z = (double*)calloc((long long)ldz * m2n + PAD * (long long)ldz, sizeof(double));
        int* ISuppZ = (int*)calloc(2 * m2n + PAD, sizeof(int));
        double* RWork = (double*)calloc(lrwork + PAD, sizeof(double));
        int* IWork = (int*)calloc(liwork + PAD, sizeof(int));
        int info_call;

        dstexr_(&m2n, tgk_d, tgk_e, &wil, &wiu,
                W, Z, &ldz, ISuppZ,
                RWork, &lrwork, IWork, &liwork, &info_call);

        printf("UPPER HALF (wil=%d, wiu=%d): info=%d\n", wil, wiu, info_call);
        printf("  W: ");
        for (int i = 0; i < n; i++) printf("%.6f ", W[i]);
        printf("\n");

        free(W); free(Z); free(ISuppZ); free(RWork); free(IWork);
    }

    free(tgk_d); free(tgk_e);
    return 0;
}
