#include "bidiag_svd.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

// Read a .dat file (format: n on first line, then rows: index d_i e_i)
bool read_dat(const char* path, std::vector<double>& d, std::vector<double>& e, int& n) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    std::getline(f, line);
    std::istringstream iss0(line);
    iss0 >> n;
    d.resize(n);
    e.resize(n > 1 ? n - 1 : 0);
    for (int i = 0; i < n; i++) {
        std::getline(f, line);
        std::istringstream iss(line);
        int idx;
        double di, ei;
        iss >> idx >> di >> ei;
        d[i] = di;
        if (i < n - 1) e[i] = ei;
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <matrix.dat>\n", argv[0]);
        return 1;
    }

    std::vector<double> d, e;
    int n;
    if (!read_dat(argv[1], d, e, n)) {
        printf("Failed to read %s\n", argv[1]);
        return 1;
    }

    printf("Matrix: %s, n=%d\n", argv[1], n);
    printf("d: "); for (int i = 0; i < std::min(n, 10); i++) printf("%.4e ", d[i]); printf("\n");
    printf("e: "); for (int i = 0; i < std::min(n-1, 10); i++) printf("%.4e ", e[i]); printf("\n\n");

    BidiagSVDResult result = bidiag_svd(n, d.data(), e.data());
    printf("info=%d\n", result.info);
    printf("sigma: "); for (int i = 0; i < std::min(n, 10); i++) printf("%.6e ", result.sigma[i]); printf("\n\n");

    // Compute V^T * V (column orthogonality)
    printf("=== V^T * V (column dot products) ===\n");
    double max_vtv_offdiag = 0.0;
    for (int j = 0; j < n; j++) {
        for (int k = j; k < n; k++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++)
                dot += result.V[j * n + i] * result.V[k * n + i];
            double target = (j == k) ? 1.0 : 0.0;
            double err = std::abs(dot - target);
            if (j != k) max_vtv_offdiag = std::max(max_vtv_offdiag, err);
            if (err > 1e-10 || j == k)
                if (n <= 20 || err > 1e-10)
                    printf("  V^T*V[%d,%d] = %+.6e (err=%.2e)\n", j, k, dot, err);
        }
    }
    printf("  max off-diag |V^T*V - I|: %.2e\n\n", max_vtv_offdiag);

    // Compute V * V^T (row orthogonality) - this is what evaluate checks
    printf("=== V * V^T (row-based, evaluate metric) ===\n");
    double max_vvt_err = 0.0;
    int worst_i = -1, worst_j = -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double qqt = 0.0;
            for (int k = 0; k < n; k++)
                qqt += result.V[k * n + i] * result.V[k * n + j];
            double target = (i == j) ? 1.0 : 0.0;
            double err = std::abs(qqt - target);
            if (err > max_vvt_err) {
                max_vvt_err = err;
                worst_i = i;
                worst_j = j;
            }
            if (err > 1e-10)
                printf("  V*V^T[%d,%d] = %+.6e (err=%.2e)\n", i, j, qqt, err);
        }
    }
    double eps = 2.2204460492503131e-16;
    printf("  max |V*V^T - I|: %.2e at (%d,%d)\n", max_vvt_err, worst_i, worst_j);
    printf("  ortV metric: %.2e\n\n", max_vvt_err / (n * eps));

    // Same for U
    printf("=== U * U^T (evaluate metric) ===\n");
    double max_uut_err = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double qqt = 0.0;
            for (int k = 0; k < n; k++)
                qqt += result.U[k * n + i] * result.U[k * n + j];
            double target = (i == j) ? 1.0 : 0.0;
            double err = std::abs(qqt - target);
            max_uut_err = std::max(max_uut_err, err);
            if (err > 1e-10)
                printf("  U*U^T[%d,%d] = %+.6e (err=%.2e)\n", i, j, qqt, err);
        }
    }
    printf("  max |U*U^T - I|: %.2e\n", max_uut_err);
    printf("  ortU metric: %.2e\n", max_uut_err / (n * eps));

    return 0;
}
