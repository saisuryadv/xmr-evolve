// Evaluation framework for bidiagonal SVD
// Computes residual, orthogonality metrics matching LAPACK's DBDCHKRSLT
// Reads test matrices from STCollection .dat files
// Exhaustive test suite: STCollection + papers (Demmel 2008, Grosser-Lang 2001,
//   Marques 2020, Willems-Lang 2012/2013, Dhillon-Parlett-Vomel 2005)
// Output: metrics for OpenEvolve

#include "bidiag_svd.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <map>
#include <tuple>
#include <set>
#include <cfloat>
#include <csignal>

// Global cumulative timeout: 30 seconds total for the entire evaluation.
// If exceeded, print score=0 and exit immediately.
static const double CUMULATIVE_TIMEOUT_SEC = 30.0;
// Per-test timeout: 2 seconds per individual test.
// Any single test exceeding this is auto-FAIL with huge metrics.
static const double PER_TEST_TIMEOUT_SEC = 2.0;
static std::chrono::high_resolution_clock::time_point g_eval_start;

namespace fs = std::filesystem;

static const double EPS = 2.2204460492503131e-16;
static const double SQRT_EPS = 1.4901161193847656e-08;

// ============================================================
// Metric computation (matches LAPACK DBDCHKRSLT)
// ============================================================

struct SVDMetrics {
    double residual;    // ||B - U*Σ*V'||_max / (||B|| * n * eps)
    double ortho_u;     // ||I - U*U'||_max / (n * eps)
    double ortho_v;     // ||I - V*V'||_max / (n * eps)
    double time_sec;
    int info;
    std::string name;
};

double compute_bnorm(int n, const double* d, const double* e) {
    double bnorm = std::abs(d[n - 1]);
    for (int i = 0; i < n - 1; i++) {
        double val = std::abs(d[i]) + std::abs(e[i]);
        bnorm = std::max(bnorm, val);
    }
    return bnorm;
}

double compute_residual(int n, const double* d, const double* e,
                        const double* sigma, const double* U, const double* V) {
    double eps = dlamch_("E", 1);
    double bnorm = compute_bnorm(n, d, e);
    if (bnorm == 0.0) bnorm = 1.0;

    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double b_ij = 0.0;
            if (i == j) b_ij = d[i];
            else if (j == i + 1) b_ij = e[i];

            double usv_ij = 0.0;
            for (int k = 0; k < n; k++) {
                usv_ij += U[k * n + i] * sigma[k] * V[k * n + j];
            }

            double err = std::abs(b_ij - usv_ij);
            max_err = std::max(max_err, err);
        }
    }

    return max_err / (bnorm * n * eps);
}

double compute_ortho(int n, const double* Q) {
    double eps = dlamch_("E", 1);
    double max_err = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double qqt = 0.0;
            for (int k = 0; k < n; k++) {
                qqt += Q[k * n + i] * Q[k * n + j];
            }
            double target = (i == j) ? 1.0 : 0.0;
            double err = std::abs(target - qqt);
            max_err = std::max(max_err, err);
        }
    }

    return max_err / (n * eps);
}

// Check cumulative timeout — if exceeded, print score=0 and hard-exit
static void check_cumulative_timeout() {
    auto now = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(now - g_eval_start).count();
    if (elapsed > CUMULATIVE_TIMEOUT_SEC) {
        printf("\n*** TIMEOUT: evaluation exceeded %.0f seconds (%.1fs elapsed) ***\n",
               CUMULATIVE_TIMEOUT_SEC, elapsed);
        printf("*** Score = 0. Your algorithm has a pathological performance bug. ***\n");
        printf("\n=== OPENEVOLVE METRICS ===\n");
        printf("pass_rate=0.0000\n");
        printf("avg_residual=1000000.0000\n");
        printf("avg_ortho_u=1000000.0000\n");
        printf("avg_ortho_v=1000000.0000\n");
        printf("max_residual=1000000.0000\n");
        printf("max_ortho_u=1000000.0000\n");
        printf("max_ortho_v=1000000.0000\n");
        printf("pass_avg_residual=0.0000\n");
        printf("pass_avg_ortho_u=0.0000\n");
        printf("pass_avg_ortho_v=0.0000\n");
        printf("pass_max_residual=0.0000\n");
        printf("pass_max_ortho_u=0.0000\n");
        printf("pass_max_ortho_v=0.0000\n");
        printf("pass_worst_scaling=100.0000\n");
        printf("composite_score=0.0000\n");
        std::exit(1);
    }
}

SVDMetrics evaluate_matrix(int n, const double* d, const double* e, const std::string& name) {
    SVDMetrics metrics;
    metrics.name = name;

    // Check cumulative timeout before starting this test
    check_cumulative_timeout();

    auto t0 = std::chrono::high_resolution_clock::now();
    BidiagSVDResult result = bidiag_svd(n, d, e);
    auto t1 = std::chrono::high_resolution_clock::now();

    metrics.time_sec = std::chrono::duration<double>(t1 - t0).count();
    metrics.info = result.info;

    // Per-test timeout: if this single test took too long, auto-FAIL
    if (metrics.time_sec > PER_TEST_TIMEOUT_SEC) {
        printf("  *** PER-TEST TIMEOUT: %s n=%d took %.2fs (limit %.1fs) ***\n",
               name.c_str(), n, metrics.time_sec, PER_TEST_TIMEOUT_SEC);
        metrics.residual = 1e10;
        metrics.ortho_u = 1e10;
        metrics.ortho_v = 1e10;
        return metrics;
    }

    if (result.info != 0 || result.sigma.empty()) {
        metrics.residual = 1e10;
        metrics.ortho_u = 1e10;
        metrics.ortho_v = 1e10;
        return metrics;
    }

    metrics.residual = compute_residual(n, d, e, result.sigma.data(),
                                         result.U.data(), result.V.data());
    metrics.ortho_u = compute_ortho(n, result.U.data());
    metrics.ortho_v = compute_ortho(n, result.V.data());

    return metrics;
}

// ============================================================
// Test matrix loading from STCollection .dat files
// ============================================================

struct TestMatrix {
    std::string name;
    int n;
    std::vector<double> d, e;
};

TestMatrix load_stcollection(const std::string& path) {
    TestMatrix tm;
    tm.name = fs::path(path).stem().string();

    std::ifstream f(path);
    if (!f.is_open()) {
        fprintf(stderr, "Cannot open %s\n", path.c_str());
        tm.n = 0;
        return tm;
    }

    std::string line;
    if (std::getline(f, line)) {
        std::istringstream iss(line);
        iss >> tm.n;
    }

    if (tm.n <= 0) {
        tm.n = 0;
        return tm;
    }

    tm.d.resize(tm.n);
    tm.e.resize(tm.n > 1 ? tm.n - 1 : 0);

    for (int i = 0; i < tm.n && std::getline(f, line); i++) {
        std::istringstream iss(line);
        int idx;
        double di, ei;
        iss >> idx >> di >> ei;
        tm.d[i] = di;
        if (i < tm.n - 1) {
            tm.e[i] = ei;
        }
    }

    return tm;
}

// ============================================================
// Deterministic RNG (to avoid srand/rand state issues)
// ============================================================

struct DetRNG {
    uint32_t state;
    DetRNG(uint32_t seed = 42) : state(seed) {}
    uint32_t next() {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }
    double uniform() { return (double)(next() & 0x7FFFFFFF) / 0x7FFFFFFF; }
};

// ============================================================
// Exhaustive synthetic test matrices
// Sources: Demmel 2008, Grosser-Lang 2001 (dmatgen.f IDs 110-244),
//          Marques 2020, Willems-Lang 2012/2013, Dhillon-Parlett-Vomel 2005,
//          test_scaling_parallel.py (22 patterns)
// ============================================================

TestMatrix make_adversarial(const std::string& name, int n) {
    TestMatrix tm;
    tm.name = name;
    tm.n = n;
    tm.d.resize(n);
    tm.e.resize(n > 1 ? n - 1 : 0);

    DetRNG rng(42);

    // ================================================================
    // GROUP 1: Original adversarial patterns (from test_scaling_parallel.py)
    // ================================================================

    if (name == "exponential_graded") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 16.0 / n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = tm.d[i] * 0.5;
    } else if (name == "glued_repeated") {
        int block = 10;
        std::fill(tm.d.begin(), tm.d.end(), 0.0);
        std::fill(tm.e.begin(), tm.e.end(), 0.0);
        for (int b = 0; b < n / block; b++) {
            int s = b * block;
            for (int i = s; i < std::min(s + block, n); i++)
                tm.d[i] = 1.0 + b * 1e-14;
            for (int i = s; i < std::min(s + block - 1, n - 1); i++)
                tm.e[i] = 0.5;
            if (s > 0 && s - 1 < n - 1)
                tm.e[s - 1] = 1e-15;
        }
    } else if (name == "saw_tooth") {
        for (int i = 0; i < n; i++)
            tm.d[i] = (i % 2 == 0) ? 1.0 : 1e-8;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.5;
    } else if (name == "stemr_killer") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 20.0 / n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]);
    } else if (name == "huge_condition") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 15.0 / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.9;
    } else if (name == "spike") {
        for (int i = 0; i < n; i++)
            tm.d[i] = 0.01;
        tm.d[n / 2] = 100.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.01;
    } else if (name == "wilkinson_like") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::abs(i - n / 2) + 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 1.0;
    } else if (name == "two_clusters") {
        for (int i = 0; i < n; i++)
            tm.d[i] = (i < n / 2) ? 1.0 : 1e-8;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 1e-10;
        if (n / 2 - 1 >= 0 && n / 2 - 1 < n - 1)
            tm.e[n / 2 - 1] = 1e-16;
    } else if (name == "random_uniform") {
        for (int i = 0; i < n; i++)
            tm.d[i] = 0.1 + rng.uniform() * 9.9;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.01 + rng.uniform() * 0.99;
    } else if (name == "diagonal_only") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 8.0 / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.0;
    } else if (name == "constant") {
        for (int i = 0; i < n; i++)
            tm.d[i] = 3.14;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.0;
    } else if (name == "all_equal_nontrivial") {
        for (int i = 0; i < n; i++) tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "one_big_cluster") {
        for (int i = 0; i < n; i++) tm.d[i] = 1.0 + rng.uniform() * 1e-12;
        for (int i = 0; i < n - 1; i++) tm.e[i] = rng.uniform() * 1e-12;
    } else if (name == "arithmetic_progression") {
        for (int i = 0; i < n; i++) tm.d[i] = 1.0 - (double)i / (n - 1) * (1.0 - 1e-8);
        for (int i = 0; i < n - 1; i++) tm.e[i] = (tm.d[i] + tm.d[i + 1]) / 4.0;
    } else if (name == "many_near_zero") {
        for (int i = 0; i < n; i++) tm.d[i] = 1e-15;
        for (int i = 0; i < std::min(5, n); i++)
            tm.d[i] = std::pow(10.0, -(double)i);
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1e-16;
        for (int i = 0; i < std::min(4, n - 1); i++) tm.e[i] = tm.d[i] * 0.1;
    } else if (name == "random_dense_clusters") {
        for (int i = 0; i < n; i++) tm.d[i] = 0.1 + rng.uniform() * 9.9;
        std::sort(tm.d.begin(), tm.d.end(), std::greater<double>());
        for (int i = 0; i < n - 10; i += 10) {
            double base = tm.d[i];
            for (int j = i; j < std::min(i + 10, n); j++)
                tm.d[j] = base + rng.uniform() * 1e-13;
        }
        for (int i = 0; i < n - 1; i++) tm.e[i] = 0.01 + rng.uniform() * 0.09;
    } else if (name == "constant_d_graded_e") {
        for (int i = 0; i < n; i++) tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::pow(10.0, -(double)i * 16.0 / (n - 1));
    } else if (name == "random_clustered_5") {
        double centers[5];
        for (int c = 0; c < 5; c++) centers[c] = 0.5 + rng.uniform() * 4.5;
        for (int i = 0; i < n; i++) tm.d[i] = centers[i % 5] + rng.uniform() * 2e-12 - 1e-12;
        std::sort(tm.d.begin(), tm.d.end(), std::greater<double>());
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1e-13 + rng.uniform() * (1e-12 - 1e-13);
    } else if (name == "alternating_sign") {
        for (int i = 0; i < n; i++) tm.d[i] = (i % 2 == 0) ? 1.0 : -1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 0.3;
    } else if (name == "step_function") {
        for (int i = 0; i < n; i++) {
            if (i < n / 3) tm.d[i] = 10.0;
            else if (i < 2 * n / 3) tm.d[i] = 1.0;
            else tm.d[i] = 0.1;
        }
        for (int i = 0; i < n - 1; i++) tm.e[i] = 0.5;
    } else if (name == "three_clusters") {
        for (int i = 0; i < n; i++) {
            if (i < n / 3) tm.d[i] = 1.0;
            else if (i < 2 * n / 3) tm.d[i] = 1e-4;
            else tm.d[i] = 1e-8;
        }
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1e-10;
    } else if (name == "random_sparse_e") {
        for (int i = 0; i < n; i++) tm.d[i] = 0.5 + rng.uniform() * 1.5;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 0.0;
        int nz = std::max(1, n / 10);
        for (int k = 0; k < nz; k++) {
            int idx = rng.next() % (n - 1);
            tm.e[idx] = 0.1 + rng.uniform() * 0.9;
        }

    // ================================================================
    // GROUP 2: Marques 2020 matrices
    // ================================================================

    } else if (name == "chkbd") {
        // CHKBD: d[i] = 10^(-(2i+1)), e[i] = 10^(-(2(i+1)))
        // Triggers DBDSVDX reorthogonalization bug. Cond ~ 10^(2n-2)
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)(2 * i + 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::pow(10.0, -(double)(2 * (i + 1)));
    } else if (name == "marques_graded") {
        // Bidiagonal with singular values spanning 1 to sqrt(eps)
        // Inspired by Marques 2020 Appendix B failure (lambda_j = c^((j-1)/(n-1)))
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(SQRT_EPS, (double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.3;

    // ================================================================
    // GROUP 3: Grosser-Lang entry-based matrices (dmatgen.f IDs 200-244)
    // ================================================================

    } else if (name == "gl_abcon0") {
        // ID 200: Heat equation tridiag, alpha=2, beta=1. d=alpha, e=beta
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon1") {
        // ID 201: ABCON with d[0] = alpha-beta
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 1.0; // alpha - beta = 2 - 1
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon2") {
        // ID 202: ABCON with d[0]=alpha-beta, d[n-1]=alpha+beta
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 1.0;    // alpha - beta
        tm.d[n-1] = 3.0;  // alpha + beta
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon3") {
        // ID 203: ABCON with d[0]=d[n-1]=alpha+beta
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 3.0;    // alpha + beta
        tm.d[n-1] = 3.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_random") {
        // ID 210: Random bidiagonal (DLARAN with fixed seed)
        for (int i = 0; i < n; i++) tm.d[i] = rng.uniform();
        for (int i = 0; i < n - 1; i++) tm.e[i] = rng.uniform();
    } else if (name == "gl_gradp") {
        // ID 220: Graded+ (geometric from right, d[n-1]=1, d[i]=beta*d[i+1])
        double beta = 0.5;
        tm.d[n-1] = 1.0;
        for (int i = n - 2; i >= 0; i--) {
            tm.d[i] = beta * tm.d[i + 1];
            tm.e[i] = tm.d[i];
        }
    } else if (name == "gl_gradm") {
        // ID 221: Graded- (geometric from left, d[0]=1, d[i+1]=beta*d[i])
        double beta = 0.5;
        tm.d[0] = 1.0;
        for (int i = 0; i < n - 1; i++) {
            tm.d[i + 1] = beta * tm.d[i];
            tm.e[i] = tm.d[i + 1];
        }
    } else if (name == "gl_wilkp") {
        // ID 222: Wilkinson W_{2L+1}^+ (exact from dmatgen.f)
        int L = (n - 1) / 2;
        double temp = (double)L;
        for (int i = 0; i < L; i++) {
            tm.d[i] = temp;
            tm.d[n - 1 - i] = temp;
            temp -= 1.0;
        }
        tm.d[L] = 0.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_wilkm") {
        // ID 223: Wilkinson W^- (antisymmetric)
        int L = (n - 1) / 2;
        double temp = (double)L;
        for (int i = 0; i < L; i++) {
            tm.d[i] = temp;
            tm.d[n - 1 - i] = -temp;
            temp -= 1.0;
        }
        tm.d[L] = 0.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_wilkw") {
        // ID 224: Fernando-Parlett Wilkinson variant
        double beta = 1.0 - (double)n / 2.0;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::abs((double)(i + 1) + beta);
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_wilk2w") {
        // ID 225: Two copies of WILKW
        double beta = 1.0 - (double)n / 4.0;
        int half = n / 2;
        for (int i = 0; i < half; i++)
            tm.d[i] = std::abs((double)(i + 1) + beta);
        for (int i = half; i < n; i++)
            tm.d[i] = tm.d[i - half];
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_clement") {
        // ID 230: Clement matrix. d=0, e[i]=sqrt(i*(n-i)) (exact from dmatgen.f)
        for (int i = 0; i < n; i++) tm.d[i] = 0.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt((double)(i + 1) * (double)(n - 1 - i));
    } else if (name == "gl_gro0") {
        // ID 240: GRO0 — d[0]=1, rest d=alpha, all e=alpha
        double alpha = 0.01;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro1") {
        // ID 241: GRO1 — d[0]=d[1]=1, rest d=alpha, all e=alpha
        double alpha = 0.01;
        tm.d[0] = 1.0;
        tm.d[1] = 1.0;
        for (int i = 2; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro2") {
        // ID 242: GRO2 — d[0..3]=1, rest d=alpha, all e=alpha
        double alpha = 0.01;
        for (int i = 0; i < std::min(4, n); i++) tm.d[i] = 1.0;
        for (int i = 4; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro3") {
        // ID 243: GRO3 — d[0..3]=1 e[0..3]=1, rest d=alpha e=alpha
        double alpha = 0.01;
        for (int i = 0; i < std::min(4, n); i++) tm.d[i] = 1.0;
        for (int i = 4; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < std::min(4, n - 1); i++) tm.e[i] = 1.0;
        for (int i = 4; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_bwilkp") {
        // ID 244: bWILK+ (Wilkinson+ shifted by +1)
        int L = (n - 1) / 2;
        double temp = (double)L;
        for (int i = 0; i < L; i++) {
            tm.d[i] = temp;
            tm.d[n - 1 - i] = temp;
            temp -= 1.0;
        }
        tm.d[L] = 0.0;
        // Shift by +1
        for (int i = 0; i < n; i++) tm.d[i] += 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;

    // ================================================================
    // GROUP 4: Grosser-Lang spectrum-based (IDs 110-121)
    // Prescribed singular values as bidiagonal d, with geometric-mean coupling e
    // ================================================================

    } else if (name == "gl_ones") {
        // ID 110: All singular values = 1
        for (int i = 0; i < n; i++) tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = EPS;
    } else if (name == "gl_uniform_eps") {
        // ID 111: Uniform eps-apart: sv[j] = j*eps for j=1..n-1, sv[n]=1
        for (int i = 0; i < n - 1; i++) tm.d[i] = (double)(i + 1) * EPS;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_uniform_sqrteps") {
        // ID 112: sv[0]=eps, sv[j]=1+j*eps^(1/4) for j=1..n-2, sv[n-1]=2
        double eps14 = std::pow(EPS, 0.25);
        tm.d[0] = EPS;
        for (int i = 1; i < n - 1; i++) tm.d[i] = 1.0 + (double)i * eps14;
        tm.d[n - 1] = 2.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_uniform_eps_to_1") {
        // ID 113: Uniform distribution eps to 1
        double step = (1.0 - EPS) / (n - 1);
        for (int i = 0; i < n; i++) tm.d[i] = EPS + (double)i * step;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_geometric_eps_to_1") {
        // ID 115: Geometric distribution eps to 1
        double base = std::pow(EPS, 1.0 / (n - 1));
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(base, (double)(n - 1 - i));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_random_spectrum") {
        // ID 117: Random singular values in (0,1)
        for (int i = 0; i < n; i++) tm.d[i] = rng.uniform();
        std::sort(tm.d.begin(), tm.d.end(), std::greater<double>());
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_clustered_at_1") {
        // ID 118: sv[0]=eps^2, sv[j]~1 with tiny perturbation
        tm.d[0] = EPS * EPS;
        for (int i = 1; i < n; i++)
            tm.d[i] = 1.0 * (10.0 * n + rng.uniform()) / (10.0 * n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_pm1") {
        // ID 119: sv[0]=eps, rest are ±1 randomly
        tm.d[0] = EPS;
        for (int i = 1; i < n; i++)
            tm.d[i] = (rng.uniform() > 0.5) ? 1.0 : -1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_eps") {
        // ID 120: sv[j]~eps^2 with tiny perturbation, sv[n-1]=1
        for (int i = 0; i < n - 1; i++)
            tm.d[i] = EPS * EPS * (10.0 * n + rng.uniform()) / (10.0 * n);
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_pmeps") {
        // ID 121: sv[j]=±eps randomly, sv[n-1]=1
        for (int i = 0; i < n - 1; i++)
            tm.d[i] = (rng.uniform() > 0.5) ? EPS : -EPS;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;

    // ================================================================
    // GROUP 5: Demmel 2008 — Strongly clustered (S1, S2)
    // Bidiagonal with prescribed singular value distributions
    // kappa = 1/eps for 'e' variants, 1/sqrt(eps) for 's' variants
    // ================================================================

    } else if (name == "demmel_S1pe") {
        // S1 MODE 1: d[0]=1, d[1..n-1]=1/kappa. kappa=1/eps, positive
        double kappa = 1.0 / EPS;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S1ps") {
        // S1 MODE 1: kappa=1/sqrt(eps)
        double kappa = 1.0 / SQRT_EPS;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S2pe") {
        // S2 MODE 2: d[0..n-2]=1, d[n-1]=1/kappa. kappa=1/eps
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = 1.0;
        tm.d[n - 1] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S2ps") {
        // S2 MODE 2: kappa=1/sqrt(eps)
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = 1.0;
        tm.d[n - 1] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;

    // ================================================================
    // GROUP 6: Demmel 2008 — Weakly clustered (W1, W2, W3)
    // ================================================================

    } else if (name == "demmel_W1") {
        // W1: d[i]=(i+1)/kappa for i=0..n-2, d[n-1]=1
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = (double)(i + 1) / kappa;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_W2") {
        // W2: d[0]=1/kappa, d[n-1]=2, d[1..n-2]=1+j/sqrt(kappa)
        double kappa = 1.0 / EPS;
        double sqk = std::sqrt(kappa);
        tm.d[0] = 1.0 / kappa;
        for (int i = 1; i < n - 1; i++) tm.d[i] = 1.0 + (double)i / sqk;
        tm.d[n - 1] = 2.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_W3") {
        // W3: d[i]=1+100*(i+1)/kappa (tight cluster near 1)
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++) tm.d[i] = 1.0 + 100.0 * (double)(i + 1) / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;

    // ================================================================
    // GROUP 7: Demmel 2008 — Geometric (G1, G2) and Uniform (U1, U2)
    // ================================================================

    } else if (name == "demmel_G1") {
        // G1 MODE 3: d[i] = kappa^(-(i)/(n-1)), geometric spacing
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_G1s") {
        // G1 with kappa=1/sqrt(eps)
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_U1") {
        // U1 MODE 4: d[i] = 1 - i/(n-1) * (1-1/kappa), arithmetic
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = 1.0 - (double)i / (n - 1) * (1.0 - 1.0 / kappa);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_U1s") {
        // U1 with kappa=1/sqrt(eps)
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = 1.0 - (double)i / (n - 1) * (1.0 - 1.0 / kappa);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;

    // ================================================================
    // GROUP 8: Exact Wilkinson and Glued Wilkinson (Dhillon-Parlett-Vomel 2005)
    // ================================================================

    } else if (name == "wilkinson_exact") {
        // Exact W_{2m+1}^+ as bidiagonal: d[i]=|m-(i)|, e[i]=1
        // Requires n odd; if even, use n-1
        int m = (n - 1) / 2;
        for (int i = 0; i < n; i++)
            tm.d[i] = (double)std::abs(m - i);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 1.0;
    } else if (name == "glued_wilkinson") {
        // p copies of W_{21}^+ glued with gamma=sqrt(eps)
        // Each block is W_{21}^+ (m=10, size 21), total n = p*21
        double gamma = SQRT_EPS;
        int block_size = 21;
        int p = std::max(1, n / block_size);
        tm.n = p * block_size;
        tm.d.resize(tm.n);
        tm.e.resize(tm.n - 1);
        for (int b = 0; b < p; b++) {
            int s = b * block_size;
            int m = (block_size - 1) / 2; // m=10
            for (int i = 0; i < block_size; i++)
                tm.d[s + i] = (double)std::abs(m - i);
            for (int i = 0; i < block_size - 1; i++)
                tm.e[s + i] = 1.0;
            // Glue between blocks
            if (b < p - 1)
                tm.e[s + block_size - 1] = gamma;
        }
    } else if (name == "glued_wilkinson_tight") {
        // 100 copies of W_{21}^+ glued with gamma=1e-11 (Marques case 966)
        double gamma = 1e-11;
        int block_size = 21;
        int p = std::max(1, n / block_size);
        tm.n = p * block_size;
        tm.d.resize(tm.n);
        tm.e.resize(tm.n - 1);
        for (int b = 0; b < p; b++) {
            int s = b * block_size;
            int m = (block_size - 1) / 2;
            for (int i = 0; i < block_size; i++)
                tm.d[s + i] = (double)std::abs(m - i);
            for (int i = 0; i < block_size - 1; i++)
                tm.e[s + i] = 1.0;
            if (b < p - 1)
                tm.e[s + block_size - 1] = gamma;
        }

    // ================================================================
    // GROUP 9: Willems-Lang 2012 specific matrices
    // ================================================================

    } else if (name == "wl_example48") {
        // Example 4.8: B=[1,1;0,alpha] with alpha~eps (element growth test)
        // For larger n, extend pattern: d graded from 1 to eps, e=1 (non-trivial coupling)
        if (n <= 2) {
            tm.n = 2;
            tm.d.resize(2);
            tm.e.resize(1);
            tm.d[0] = 1.0;
            tm.d[1] = EPS;
            tm.e[0] = 1.0;
        } else {
            for (int i = 0; i < n; i++)
                tm.d[i] = std::pow(EPS, (double)i / (n - 1));
            for (int i = 0; i < n - 1; i++)
                tm.e[i] = 1.0;
        }

    // ================================================================
    // GROUP 10: Parlett-Dhillon 2000 diagnostic matrices (as bidiagonal)
    // ================================================================

    } else if (name == "pd_T0") {
        // T0 (3x3, benign): diag=[1, 1+eps, 1], off=[0.5, 0.5]
        tm.n = 3;
        tm.d.resize(3);
        tm.e.resize(2);
        tm.d[0] = 1.0; tm.d[1] = 1.0 + EPS; tm.d[2] = 1.0;
        tm.e[0] = 0.5; tm.e[1] = 0.5;
    } else if (name == "pd_T1") {
        // T1 (3x3, hard): diag=[1, 1+eps, 1], off=[1, 1]
        tm.n = 3;
        tm.d.resize(3);
        tm.e.resize(2);
        tm.d[0] = 1.0; tm.d[1] = 1.0 + EPS; tm.d[2] = 1.0;
        tm.e[0] = 1.0; tm.e[1] = 1.0;

    // ================================================================
    // GROUP 11: Additional stress tests from papers
    // ================================================================

    } else if (name == "zero_diagonal") {
        // All-zero diagonal, tests handling of d=0
        for (int i = 0; i < n; i++) tm.d[i] = 0.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "single_element") {
        // n=1 edge case
        tm.n = 1;
        tm.d.resize(1);
        tm.e.resize(0);
        tm.d[0] = 3.14;
    } else if (name == "near_overflow") {
        // Entries near overflow/underflow boundary
        double big = 1e150;
        for (int i = 0; i < n; i++) tm.d[i] = big * std::pow(0.9, (double)i);
        for (int i = 0; i < n - 1; i++) tm.e[i] = big * std::pow(0.9, (double)i) * 0.5;
    } else if (name == "near_underflow") {
        // Entries near underflow
        double tiny = 1e-150;
        for (int i = 0; i < n; i++) tm.d[i] = tiny * std::pow(1.1, (double)i);
        for (int i = 0; i < n - 1; i++) tm.e[i] = tiny * std::pow(1.1, (double)i) * 0.5;
    } else if (name == "mixed_signs") {
        // Random mix of positive and negative d entries with large e
        for (int i = 0; i < n; i++)
            tm.d[i] = (rng.uniform() > 0.5 ? 1.0 : -1.0) * (0.1 + rng.uniform() * 9.9);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.5 + rng.uniform() * 4.5;
    } else if (name == "checkerboard") {
        // Alternating large/small with strong coupling
        for (int i = 0; i < n; i++)
            tm.d[i] = (i % 2 == 0) ? 1.0 : EPS;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1]));
    // ================================================================
    // GROUP 12: Condition number variants
    // Same structural patterns but with different condition numbers (κ)
    // Matches Synth methodology: vary κ across {1e4, 1e8, 1e12, 1e16}
    // ================================================================

    } else if (name == "exponential_graded_k4") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 4.0 / n);
        for (int i = 0; i < n - 1; i++) tm.e[i] = tm.d[i] * 0.5;
    } else if (name == "exponential_graded_k8") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 8.0 / n);
        for (int i = 0; i < n - 1; i++) tm.e[i] = tm.d[i] * 0.5;
    } else if (name == "exponential_graded_k12") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 12.0 / n);
        for (int i = 0; i < n - 1; i++) tm.e[i] = tm.d[i] * 0.5;
    } else if (name == "stemr_killer_k5") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 5.0 / n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]);
    } else if (name == "stemr_killer_k10") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 10.0 / n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]);
    } else if (name == "huge_condition_k5") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 5.0 / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.9;
    } else if (name == "huge_condition_k10") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)i * 10.0 / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.9;
    } else if (name == "demmel_G1_k4") {
        double kappa = 1e4;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_G1_k8") {
        double kappa = 1e8;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_G1_k12") {
        double kappa = 1e12;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_S1pe_k4") {
        double kappa = 1e4;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S1pe_k8") {
        double kappa = 1e8;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "chkbd_4") {
        // CHKBD at n=4 (hardcoded): κ ~ 10^6
        tm.n = 4;
        tm.d.resize(4);
        tm.e.resize(3);
        for (int i = 0; i < 4; i++)
            tm.d[i] = std::pow(10.0, -(double)(2 * i + 1));
        for (int i = 0; i < 3; i++)
            tm.e[i] = std::pow(10.0, -(double)(2 * (i + 1)));
    } else if (name == "chkbd_16") {
        // CHKBD at n=16 (hardcoded): κ ~ 10^30
        tm.n = 16;
        tm.d.resize(16);
        tm.e.resize(15);
        for (int i = 0; i < 16; i++)
            tm.d[i] = std::pow(10.0, -(double)(2 * i + 1));
        for (int i = 0; i < 15; i++)
            tm.e[i] = std::pow(10.0, -(double)(2 * (i + 1)));
    } else if (name == "marques_graded_k4") {
        double c = 1e-4;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(c, (double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.3;
    } else if (name == "marques_graded_k8") {
        double c = 1e-8;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(c, (double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.3;
    } else {
        // Default: identity-like
        for (int i = 0; i < n; i++)
            tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.0;
    }

    return tm;
}

// ============================================================
// Main evaluation driver
// ============================================================

int main(int argc, char** argv) {
    // Start the cumulative timeout clock
    g_eval_start = std::chrono::high_resolution_clock::now();

    std::string stcoll_dir = "";
    int adv_size = 100;

    if (argc > 1) stcoll_dir = argv[1];
    if (argc > 2) adv_size = std::atoi(argv[2]);

    std::vector<SVDMetrics> all_metrics;
    int pass_count = 0, total_count = 0;

    const double RES_THRESH = 7.0;
    const double ORTHO_THRESH = 5.0;

    // --- STCollection matrices ---
    if (!stcoll_dir.empty() && fs::exists(stcoll_dir)) {
        printf("=== STCollection Matrices ===\n");
        for (auto& entry : fs::directory_iterator(stcoll_dir)) {
            if (entry.path().extension() == ".dat" &&
                entry.path().filename().string().substr(0, 2) == "B_") {
                TestMatrix tm = load_stcollection(entry.path().string());
                if (tm.n == 0) continue;

                SVDMetrics m = evaluate_matrix(tm.n, tm.d.data(), tm.e.data(), tm.name);
                all_metrics.push_back(m);
                total_count++;

                bool pass = (m.info == 0 && m.residual <= RES_THRESH &&
                             m.ortho_u <= ORTHO_THRESH && m.ortho_v <= ORTHO_THRESH);
                if (pass) pass_count++;

                printf("  %-35s n=%4d  res=%8.3f  ortU=%8.3f  ortV=%8.3f  t=%.4fs  %s\n",
                       tm.name.c_str(), tm.n, m.residual, m.ortho_u, m.ortho_v,
                       m.time_sec, pass ? "PASS" : "FAIL");
            }
        }
    }

    // --- Adversarial matrices (exhaustive, multi-size) ---
    // Multi-size testing catches size-dependent bugs (Wilkinson FP underflow
    // for m>=50, element growth compounding, cluster scaling with n).
    // Matches Willems-Lang 2012 Synth methodology: same spectrum type at
    // multiple sizes. Build sorted, deduplicated size list.
    std::vector<int> test_sizes;
    {
        std::set<int> sz_set = {10, 100, adv_size, adv_size * 2};
        test_sizes.assign(sz_set.begin(), sz_set.end()); // sorted, unique
    }

    std::vector<std::string> adv_names = {
        // Group 1: Original patterns (22 from test_scaling_parallel.py)
        "exponential_graded", "glued_repeated", "saw_tooth",
        "stemr_killer", "huge_condition", "spike",
        "wilkinson_like", "two_clusters", "random_uniform",
        "diagonal_only", "constant",
        "all_equal_nontrivial", "one_big_cluster", "arithmetic_progression",
        "many_near_zero", "random_dense_clusters", "constant_d_graded_e",
        "random_clustered_5", "alternating_sign", "step_function",
        "three_clusters", "random_sparse_e",
        // Group 2: Marques 2020
        "chkbd", "marques_graded",
        // Group 3: Grosser-Lang entry-based (IDs 200-244)
        "gl_abcon0", "gl_abcon1", "gl_abcon2", "gl_abcon3",
        "gl_random", "gl_gradp", "gl_gradm",
        "gl_wilkp", "gl_wilkm", "gl_wilkw", "gl_wilk2w",
        "gl_clement",
        "gl_gro0", "gl_gro1", "gl_gro2", "gl_gro3",
        "gl_bwilkp",
        // Group 4: Grosser-Lang spectrum-based (IDs 110-121)
        "gl_ones", "gl_uniform_eps", "gl_uniform_sqrteps",
        "gl_uniform_eps_to_1", "gl_geometric_eps_to_1",
        "gl_random_spectrum",
        "gl_clustered_at_1", "gl_clustered_at_pm1",
        "gl_clustered_at_eps", "gl_clustered_at_pmeps",
        // Group 5: Demmel 2008 strongly clustered
        "demmel_S1pe", "demmel_S1ps", "demmel_S2pe", "demmel_S2ps",
        // Group 6: Demmel 2008 weakly clustered
        "demmel_W1", "demmel_W2", "demmel_W3",
        // Group 7: Demmel 2008 geometric/uniform
        "demmel_G1", "demmel_G1s", "demmel_U1", "demmel_U1s",
        // Group 8: Exact Wilkinson / Glued Wilkinson
        "wilkinson_exact", "glued_wilkinson", "glued_wilkinson_tight",
        // Group 9: Willems-Lang
        "wl_example48",
        // Group 10: Parlett-Dhillon diagnostics
        "pd_T0", "pd_T1",
        // Group 11: Stress tests
        "zero_diagonal", "single_element",
        "near_overflow", "near_underflow",
        "mixed_signs", "checkerboard",
        // Group 12: Condition number variants (κ = 1e4, 1e8, 1e12, 1e16)
        // Each key pattern tested at multiple κ to match Synth methodology
        "exponential_graded_k4", "exponential_graded_k8", "exponential_graded_k12",
        "stemr_killer_k5", "stemr_killer_k10",
        "huge_condition_k5", "huge_condition_k10",
        "demmel_G1_k4", "demmel_G1_k8", "demmel_G1_k12",
        "demmel_S1pe_k4", "demmel_S1pe_k8",
        "chkbd_4", "chkbd_16",
        "marques_graded_k4", "marques_graded_k8",
    };

    // Track per-pattern timings and pass/fail at each size for scaling computation
    // Map: pattern_name -> {size -> time}
    std::map<std::string, std::map<int, double>> pattern_times;
    // Map: pattern_name -> {size -> passed}
    std::map<std::string, std::map<int, bool>> pattern_passed;

    // Track patterns that catastrophically fail at small n — skip them at larger n
    // "Catastrophic" = any metric > 1000 (way beyond threshold, no point testing larger)
    std::set<std::string> catastrophic_patterns;

    for (int sz : test_sizes) {
        printf("\n=== Adversarial Matrices (n=%d) ===\n", sz);
        for (auto& name : adv_names) {
            // Early abort: skip larger sizes for patterns that catastrophically failed
            if (catastrophic_patterns.count(name)) {
                // Still count as a test (failed) for scoring
                SVDMetrics skip_m;
                skip_m.name = name + "_n" + std::to_string(sz);
                skip_m.residual = 1e10;
                skip_m.ortho_u = 1e10;
                skip_m.ortho_v = 1e10;
                skip_m.time_sec = 0.0;
                skip_m.info = -1;
                all_metrics.push_back(skip_m);
                total_count++;
                pattern_times[name][sz] = 0.0;
                pattern_passed[name][sz] = false;
                printf("  %-35s n=%4d  SKIPPED (catastrophic fail at smaller n)\n",
                       name.c_str(), sz);
                continue;
            }

            TestMatrix tm = make_adversarial(name, sz);
            SVDMetrics m = evaluate_matrix(tm.n, tm.d.data(), tm.e.data(), tm.name);
            all_metrics.push_back(m);
            total_count++;

            bool pass = (m.info == 0 && m.residual <= RES_THRESH &&
                         m.ortho_u <= ORTHO_THRESH && m.ortho_v <= ORTHO_THRESH);
            if (pass) pass_count++;

            pattern_times[name][sz] = m.time_sec;
            pattern_passed[name][sz] = pass;

            // Mark catastrophic failures for early abort at larger sizes
            if (m.residual > 1000.0 || m.ortho_u > 1000.0 || m.ortho_v > 1000.0 ||
                m.time_sec > PER_TEST_TIMEOUT_SEC) {
                catastrophic_patterns.insert(name);
            }

            printf("  %-35s n=%4d  res=%8.3f  ortU=%8.3f  ortV=%8.3f  t=%.4fs  %s\n",
                   tm.name.c_str(), tm.n, m.residual, m.ortho_u, m.ortho_v,
                   m.time_sec, pass ? "PASS" : "FAIL");
        }
    }

    // --- Scaling analysis (passing matrices only) ---
    // Only include patterns that PASSED at ALL test sizes.
    // Compute consecutive doubling ratios from sizes >= 100.
    printf("\n=== Scaling Test (passing matrices, %zu patterns) ===\n", adv_names.size());
    std::vector<int> scale_sizes;
    for (int s : test_sizes) {
        if (s >= 100) scale_sizes.push_back(s);
    }

    std::vector<std::tuple<double, std::string, int, int>> all_ratios;
    int patterns_used = 0;
    for (auto& name : adv_names) {
        auto& times = pattern_times[name];
        auto& passed = pattern_passed[name];
        // Only include patterns that passed at every test size
        bool passed_all = true;
        for (int s : test_sizes) {
            if (!passed[s]) { passed_all = false; break; }
        }
        if (!passed_all) continue;
        patterns_used++;
        // Only use the largest size transition (most stable timing)
        int s_prev = scale_sizes[scale_sizes.size() - 2];
        int s_curr = scale_sizes[scale_sizes.size() - 1];
        double t_prev = times[s_prev];
        double t_curr = times[s_curr];
        if (t_prev <= 0.0) continue;
        double ratio = t_curr / t_prev;
        all_ratios.push_back({ratio, name, s_prev, s_curr});
    }
    std::sort(all_ratios.begin(), all_ratios.end(),
              [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });

    double worst_ratio = 0.0;
    std::string worst_pattern;
    int worst_from = 0, worst_to = 0;
    if (!all_ratios.empty()) {
        auto& [r, nm, sf, st] = all_ratios[0];
        worst_ratio = r;
        worst_pattern = nm;
        worst_from = sf;
        worst_to = st;
    }

    printf("  %d/%zu patterns pass all sizes\n", patterns_used, adv_names.size());
    printf("  Top-10 worst pass scaling ratios:\n");
    for (size_t i = 0; i < std::min((size_t)10, all_ratios.size()); i++) {
        auto& [r, nm, sf, st] = all_ratios[i];
        printf("    %.2fx  %-35s  n=%d->%d\n", r, nm.c_str(), sf, st);
    }
    printf("  pass_worst_scaling=%.2f (%s n=%d->%d, expect <=5 for O(n^2))\n",
           worst_ratio, worst_pattern.c_str(), worst_from, worst_to);

    // --- Summary ---
    printf("\n=== SUMMARY ===\n");
    printf("Pass: %d/%d\n", pass_count, total_count);

    // Averages over ALL examples (capped at 1e6 per value)
    double avg_res = 0, avg_ortU = 0, avg_ortV = 0;
    double max_res = 0, max_ortU = 0, max_ortV = 0;
    // Averages over PASSED examples only
    double pavg_res = 0, pavg_ortU = 0, pavg_ortV = 0;
    double pmax_res = 0, pmax_ortU = 0, pmax_ortV = 0;
    int pcount = 0;

    for (auto& m : all_metrics) {
        double r = std::min(m.residual, 1e6);
        double ou = std::min(m.ortho_u, 1e6);
        double ov = std::min(m.ortho_v, 1e6);
        avg_res += r;
        avg_ortU += ou;
        avg_ortV += ov;
        max_res = std::max(max_res, r);
        max_ortU = std::max(max_ortU, ou);
        max_ortV = std::max(max_ortV, ov);

        bool p = (m.info == 0 && m.residual <= RES_THRESH &&
                  m.ortho_u <= ORTHO_THRESH && m.ortho_v <= ORTHO_THRESH);
        if (p) {
            pavg_res += m.residual;
            pavg_ortU += m.ortho_u;
            pavg_ortV += m.ortho_v;
            pmax_res = std::max(pmax_res, m.residual);
            pmax_ortU = std::max(pmax_ortU, m.ortho_u);
            pmax_ortV = std::max(pmax_ortV, m.ortho_v);
            pcount++;
        }
    }
    if (!all_metrics.empty()) {
        avg_res /= all_metrics.size();
        avg_ortU /= all_metrics.size();
        avg_ortV /= all_metrics.size();
    }
    if (pcount > 0) {
        pavg_res /= pcount;
        pavg_ortU /= pcount;
        pavg_ortV /= pcount;
    }

    printf("  All  avg: res=%.2f ortU=%.2f ortV=%.2f\n", avg_res, avg_ortU, avg_ortV);
    printf("  All  max: res=%.2f ortU=%.2f ortV=%.2f\n", max_res, max_ortU, max_ortV);
    printf("  Pass avg: res=%.4f ortU=%.4f ortV=%.4f\n", pavg_res, pavg_ortU, pavg_ortV);
    printf("  Pass max: res=%.4f ortU=%.4f ortV=%.4f\n", pmax_res, pmax_ortU, pmax_ortV);

    double pass_rate = total_count > 0 ? (double)pass_count / total_count : 0;
    double score = pass_rate * 50.0
                   + std::max(0.0, 5.0 - pavg_res) * 2.0
                   + std::max(0.0, 5.0 - pavg_ortU) * 2.0
                   + std::max(0.0, 5.0 - pavg_ortV) * 2.0;

    // O(n²) is a HARD CONSTRAINT, not a soft bonus.
    // If worst scaling ratio > 5x (n doubles → time should ≤ 4x for O(n²), allow 5x for margin),
    // the algorithm is NOT O(n²) and the score is capped at 5 (compilation only).
    // This prevents gaming by using O(n³) fallbacks (e.g., DBDSQR) to boost pass rate.
    score += 5.0; // compilation bonus (we got here)
    if (worst_ratio > 5.0) {
        printf("\n*** O(n²) VIOLATION: worst scaling ratio %.2fx > 5.0x ***\n", worst_ratio);
        printf("*** Score capped at 5. An O(n³) fallback (e.g., DBDSQR) is NOT acceptable. ***\n");
        score = 5.0;  // only compilation credit
    } else {
        score += 10.0;  // full scaling bonus
    }

    printf("\n=== OPENEVOLVE METRICS ===\n");
    printf("pass_rate=%.4f\n", pass_rate);
    printf("avg_residual=%.4f\n", avg_res);
    printf("avg_ortho_u=%.4f\n", avg_ortU);
    printf("avg_ortho_v=%.4f\n", avg_ortV);
    printf("max_residual=%.4f\n", max_res);
    printf("max_ortho_u=%.4f\n", max_ortU);
    printf("max_ortho_v=%.4f\n", max_ortV);
    printf("pass_avg_residual=%.4f\n", pavg_res);
    printf("pass_avg_ortho_u=%.4f\n", pavg_ortU);
    printf("pass_avg_ortho_v=%.4f\n", pavg_ortV);
    printf("pass_max_residual=%.4f\n", pmax_res);
    printf("pass_max_ortho_u=%.4f\n", pmax_ortU);
    printf("pass_max_ortho_v=%.4f\n", pmax_ortV);
    printf("pass_worst_scaling=%.4f\n", worst_ratio);
    printf("composite_score=%.4f\n", score);

    return 0;
}
