// Dump raw S, U, V from dbdsgr_ (hgbsvd) for C-vs-Fortran comparison.
// Uses ALL 90 adversarial patterns at sizes {10, 50, 100} matching evaluator.
// For n>=100 the output files get huge (100x100 = 10K entries per matrix),
// so we use n<=100 but cover every pattern.
//
// Build (C lib):
//   g++ -std=c++17 -O2 -Isrc/clapack -o /tmp/dump_hgb_c tools/dump_hgbsvd.cpp -Llib -lxmr_c -lm
// Build (Fortran lib):
//   g++ -std=c++17 -O2 -o /tmp/dump_hgb_f tools/dump_hgbsvd.cpp -Llib -lxmr -lxmrlapack -framework Accelerate -L/opt/homebrew/lib/gcc/15 -lgfortran

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>

extern "C" {
    void dbdsgr_(const char* jobz, const char* range,
                 const int* n, double* d, double* e,
                 const double* vl, const double* vu,
                 const int* il, const int* iu,
                 const double* abstol, int* m, double* w,
                 double* z, const int* ldz, int* isuppz,
                 double* work, const int* lwork,
                 int* iwork, const int* liwork, int* info);

    double dlamch_(const char* cmach, int len);
}

static const double EPS = 2.2204460492503131e-16;
static const double SQRT_EPS = 1.4901161193847656e-08;

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

struct TestMatrix {
    std::string name;
    int n;
    std::vector<double> d, e;
};

// ============================================================
// EXACT COPY of make_adversarial from evaluate.cpp
// ============================================================
TestMatrix make_adversarial(const std::string& name, int n) {
    TestMatrix tm;
    tm.name = name;
    tm.n = n;
    tm.d.resize(n);
    tm.e.resize(n > 1 ? n - 1 : 0);

    DetRNG rng(42);

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
    } else if (name == "chkbd") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(10.0, -(double)(2 * i + 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::pow(10.0, -(double)(2 * (i + 1)));
    } else if (name == "marques_graded") {
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(SQRT_EPS, (double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.3;
    } else if (name == "gl_abcon0") {
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon1") {
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon2") {
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 1.0;
        tm.d[n-1] = 3.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_abcon3") {
        for (int i = 0; i < n; i++) tm.d[i] = 2.0;
        tm.d[0] = 3.0;
        tm.d[n-1] = 3.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_random") {
        for (int i = 0; i < n; i++) tm.d[i] = rng.uniform();
        for (int i = 0; i < n - 1; i++) tm.e[i] = rng.uniform();
    } else if (name == "gl_gradp") {
        double beta = 0.5;
        tm.d[n-1] = 1.0;
        for (int i = n - 2; i >= 0; i--) {
            tm.d[i] = beta * tm.d[i + 1];
            tm.e[i] = tm.d[i];
        }
    } else if (name == "gl_gradm") {
        double beta = 0.5;
        tm.d[0] = 1.0;
        for (int i = 0; i < n - 1; i++) {
            tm.d[i + 1] = beta * tm.d[i];
            tm.e[i] = tm.d[i + 1];
        }
    } else if (name == "gl_wilkp") {
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
        double beta = 1.0 - (double)n / 2.0;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::abs((double)(i + 1) + beta);
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_wilk2w") {
        double beta = 1.0 - (double)n / 4.0;
        int half = n / 2;
        for (int i = 0; i < half; i++)
            tm.d[i] = std::abs((double)(i + 1) + beta);
        for (int i = half; i < n; i++)
            tm.d[i] = tm.d[i - half];
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_clement") {
        for (int i = 0; i < n; i++) tm.d[i] = 0.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt((double)(i + 1) * (double)(n - 1 - i));
    } else if (name == "gl_gro0") {
        double alpha = 0.01;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro1") {
        double alpha = 0.01;
        tm.d[0] = 1.0;
        tm.d[1] = 1.0;
        for (int i = 2; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro2") {
        double alpha = 0.01;
        for (int i = 0; i < std::min(4, n); i++) tm.d[i] = 1.0;
        for (int i = 4; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_gro3") {
        double alpha = 0.01;
        for (int i = 0; i < std::min(4, n); i++) tm.d[i] = 1.0;
        for (int i = 4; i < n; i++) tm.d[i] = alpha;
        for (int i = 0; i < std::min(4, n - 1); i++) tm.e[i] = 1.0;
        for (int i = 4; i < n - 1; i++) tm.e[i] = alpha;
    } else if (name == "gl_bwilkp") {
        int L = (n - 1) / 2;
        double temp = (double)L;
        for (int i = 0; i < L; i++) {
            tm.d[i] = temp;
            tm.d[n - 1 - i] = temp;
            temp -= 1.0;
        }
        tm.d[L] = 0.0;
        for (int i = 0; i < n; i++) tm.d[i] += 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "gl_ones") {
        for (int i = 0; i < n; i++) tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = EPS;
    } else if (name == "gl_uniform_eps") {
        for (int i = 0; i < n - 1; i++) tm.d[i] = (double)(i + 1) * EPS;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_uniform_sqrteps") {
        double eps14 = std::pow(EPS, 0.25);
        tm.d[0] = EPS;
        for (int i = 1; i < n - 1; i++) tm.d[i] = 1.0 + (double)i * eps14;
        tm.d[n - 1] = 2.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_uniform_eps_to_1") {
        double step = (1.0 - EPS) / (n - 1);
        for (int i = 0; i < n; i++) tm.d[i] = EPS + (double)i * step;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_geometric_eps_to_1") {
        double base = std::pow(EPS, 1.0 / (n - 1));
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(base, (double)(n - 1 - i));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_random_spectrum") {
        for (int i = 0; i < n; i++) tm.d[i] = rng.uniform();
        std::sort(tm.d.begin(), tm.d.end(), std::greater<double>());
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "gl_clustered_at_1") {
        tm.d[0] = EPS * EPS;
        for (int i = 1; i < n; i++)
            tm.d[i] = 1.0 * (10.0 * n + rng.uniform()) / (10.0 * n);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_pm1") {
        tm.d[0] = EPS;
        for (int i = 1; i < n; i++)
            tm.d[i] = (rng.uniform() > 0.5) ? 1.0 : -1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_eps") {
        for (int i = 0; i < n - 1; i++)
            tm.d[i] = EPS * EPS * (10.0 * n + rng.uniform()) / (10.0 * n);
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "gl_clustered_at_pmeps") {
        for (int i = 0; i < n - 1; i++)
            tm.d[i] = (rng.uniform() > 0.5) ? EPS : -EPS;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S1pe") {
        double kappa = 1.0 / EPS;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S1ps") {
        double kappa = 1.0 / SQRT_EPS;
        tm.d[0] = 1.0;
        for (int i = 1; i < n; i++) tm.d[i] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S2pe") {
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = 1.0;
        tm.d[n - 1] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_S2ps") {
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = 1.0;
        tm.d[n - 1] = 1.0 / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_W1") {
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n - 1; i++) tm.d[i] = (double)(i + 1) / kappa;
        tm.d[n - 1] = 1.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_W2") {
        double kappa = 1.0 / EPS;
        double sqk = std::sqrt(kappa);
        tm.d[0] = 1.0 / kappa;
        for (int i = 1; i < n - 1; i++) tm.d[i] = 1.0 + (double)i / sqk;
        tm.d[n - 1] = 2.0;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_W3") {
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++) tm.d[i] = 1.0 + 100.0 * (double)(i + 1) / kappa;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_G1") {
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_G1s") {
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = std::pow(kappa, -(double)i / (n - 1));
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(tm.d[i] * tm.d[i + 1]) * 0.1;
    } else if (name == "demmel_U1") {
        double kappa = 1.0 / EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = 1.0 - (double)i / (n - 1) * (1.0 - 1.0 / kappa);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "demmel_U1s") {
        double kappa = 1.0 / SQRT_EPS;
        for (int i = 0; i < n; i++)
            tm.d[i] = 1.0 - (double)i / (n - 1) * (1.0 - 1.0 / kappa);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1])) * 0.1;
    } else if (name == "wilkinson_exact") {
        int m = (n - 1) / 2;
        for (int i = 0; i < n; i++)
            tm.d[i] = (double)std::abs(m - i);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 1.0;
    } else if (name == "glued_wilkinson") {
        double gamma = SQRT_EPS;
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
    } else if (name == "glued_wilkinson_tight") {
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
    } else if (name == "wl_example48") {
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
    } else if (name == "pd_T0") {
        tm.n = 3;
        tm.d.resize(3);
        tm.e.resize(2);
        tm.d[0] = 1.0; tm.d[1] = 1.0 + EPS; tm.d[2] = 1.0;
        tm.e[0] = 0.5; tm.e[1] = 0.5;
    } else if (name == "pd_T1") {
        tm.n = 3;
        tm.d.resize(3);
        tm.e.resize(2);
        tm.d[0] = 1.0; tm.d[1] = 1.0 + EPS; tm.d[2] = 1.0;
        tm.e[0] = 1.0; tm.e[1] = 1.0;
    } else if (name == "zero_diagonal") {
        for (int i = 0; i < n; i++) tm.d[i] = 0.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 1.0;
    } else if (name == "single_element") {
        tm.n = 1;
        tm.d.resize(1);
        tm.e.resize(0);
        tm.d[0] = 3.14;
    } else if (name == "near_overflow") {
        double big = 1e150;
        for (int i = 0; i < n; i++) tm.d[i] = big * std::pow(0.9, (double)i);
        for (int i = 0; i < n - 1; i++) tm.e[i] = big * std::pow(0.9, (double)i) * 0.5;
    } else if (name == "near_underflow") {
        double tiny = 1e-150;
        for (int i = 0; i < n; i++) tm.d[i] = tiny * std::pow(1.1, (double)i);
        for (int i = 0; i < n - 1; i++) tm.e[i] = tiny * std::pow(1.1, (double)i) * 0.5;
    } else if (name == "mixed_signs") {
        for (int i = 0; i < n; i++)
            tm.d[i] = (rng.uniform() > 0.5 ? 1.0 : -1.0) * (0.1 + rng.uniform() * 9.9);
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = 0.5 + rng.uniform() * 4.5;
    } else if (name == "checkerboard") {
        for (int i = 0; i < n; i++)
            tm.d[i] = (i % 2 == 0) ? 1.0 : EPS;
        for (int i = 0; i < n - 1; i++)
            tm.e[i] = std::sqrt(std::abs(tm.d[i] * tm.d[i + 1]));
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
        tm.n = 4;
        tm.d.resize(4);
        tm.e.resize(3);
        for (int i = 0; i < 4; i++)
            tm.d[i] = std::pow(10.0, -(double)(2 * i + 1));
        for (int i = 0; i < 3; i++)
            tm.e[i] = std::pow(10.0, -(double)(2 * (i + 1)));
    } else if (name == "chkbd_16") {
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
        for (int i = 0; i < n; i++) tm.d[i] = 1.0;
        for (int i = 0; i < n - 1; i++) tm.e[i] = 0.0;
    }

    return tm;
}

void run_hgbsvd(int n, const double* d, const double* e,
                double* sigma, double* U, double* V, int* info_out) {
    int PAD = 64;
    double* d_buf = (double*)calloc(n + PAD, sizeof(double));
    double* e_buf = (double*)calloc(n + PAD, sizeof(double));
    for (int i = 0; i < n; i++) d_buf[i] = d[i];
    for (int i = 0; i < n - 1; i++) e_buf[i] = e[i];

    int ldz = 2 * n;
    char jobz = 'A', range = 'A';
    double vl = 0, vu = 0, abstol = 0;
    int il = 1, iu = n, m_found;

    int lwork = -1, liwork = -1;
    double work_query;
    int iwork_query, info_ws;
    double* ws_w = (double*)calloc(n + PAD, sizeof(double));
    int* ws_isuppz = (int*)calloc(2 * ldz + PAD, sizeof(int));

    dbdsgr_(&jobz, &range, &n, d_buf, e_buf,
            &vl, &vu, &il, &iu, &abstol,
            &m_found, ws_w, nullptr, &ldz, ws_isuppz,
            &work_query, &lwork, &iwork_query, &liwork, &info_ws);

    lwork = (int)work_query + 1;
    liwork = iwork_query + 1;
    if (lwork < 20 * n) lwork = 20 * n;
    if (liwork < 10 * n) liwork = 10 * n;

    for (int i = 0; i < n; i++) d_buf[i] = d[i];
    for (int i = 0; i < n - 1; i++) e_buf[i] = e[i];

    double* W_buf = (double*)calloc(n + PAD, sizeof(double));
    double* Z_buf = (double*)calloc((long long)ldz * n + PAD * (long long)ldz, sizeof(double));
    int* isuppz_buf = (int*)calloc(2 * ldz + PAD, sizeof(int));
    double* work_buf = (double*)calloc(lwork + PAD, sizeof(double));
    int* iwork_buf = (int*)calloc(liwork + PAD, sizeof(int));
    int info_call;

    dbdsgr_(&jobz, &range, &n, d_buf, e_buf,
            &vl, &vu, &il, &iu, &abstol,
            &m_found, W_buf, Z_buf, &ldz, isuppz_buf,
            work_buf, &lwork, iwork_buf, &liwork, &info_call);

    free(ws_w); free(ws_isuppz);
    *info_out = info_call;

    if (info_call != 0) {
        for (int i = 0; i < n; i++) sigma[i] = 0;
        for (int i = 0; i < n * n; i++) { U[i] = 0; V[i] = 0; }
    } else {
        int half = ldz / 2;
        for (int i = 0; i < m_found && i < n; i++) {
            sigma[i] = W_buf[i];
            for (int r = 0; r < n; r++) {
                V[i * n + r] = Z_buf[i * ldz + r];
                U[i * n + r] = Z_buf[i * ldz + half + r];
            }
        }
    }

    free(d_buf); free(e_buf);
    free(W_buf); free(Z_buf); free(isuppz_buf);
    free(work_buf); free(iwork_buf);
}

int main() {
    // ALL 90 patterns from evaluate.cpp
    std::vector<std::string> adv_names = {
        "exponential_graded", "glued_repeated", "saw_tooth",
        "stemr_killer", "huge_condition", "spike",
        "wilkinson_like", "two_clusters", "random_uniform",
        "diagonal_only", "constant",
        "all_equal_nontrivial", "one_big_cluster", "arithmetic_progression",
        "many_near_zero", "random_dense_clusters", "constant_d_graded_e",
        "random_clustered_5", "alternating_sign", "step_function",
        "three_clusters", "random_sparse_e",
        "chkbd", "marques_graded",
        "gl_abcon0", "gl_abcon1", "gl_abcon2", "gl_abcon3",
        "gl_random", "gl_gradp", "gl_gradm",
        "gl_wilkp", "gl_wilkm", "gl_wilkw", "gl_wilk2w",
        "gl_clement",
        "gl_gro0", "gl_gro1", "gl_gro2", "gl_gro3",
        "gl_bwilkp",
        "gl_ones", "gl_uniform_eps", "gl_uniform_sqrteps",
        "gl_uniform_eps_to_1", "gl_geometric_eps_to_1",
        "gl_random_spectrum",
        "gl_clustered_at_1", "gl_clustered_at_pm1",
        "gl_clustered_at_eps", "gl_clustered_at_pmeps",
        "demmel_S1pe", "demmel_S1ps", "demmel_S2pe", "demmel_S2ps",
        "demmel_W1", "demmel_W2", "demmel_W3",
        "demmel_G1", "demmel_G1s", "demmel_U1", "demmel_U1s",
        "wilkinson_exact", "glued_wilkinson", "glued_wilkinson_tight",
        "wl_example48",
        "pd_T0", "pd_T1",
        "zero_diagonal", "single_element",
        "near_overflow", "near_underflow",
        "mixed_signs", "checkerboard",
        "exponential_graded_k4", "exponential_graded_k8", "exponential_graded_k12",
        "stemr_killer_k5", "stemr_killer_k10",
        "huge_condition_k5", "huge_condition_k10",
        "demmel_G1_k4", "demmel_G1_k8", "demmel_G1_k12",
        "demmel_S1pe_k4", "demmel_S1pe_k8",
        "chkbd_4", "chkbd_16",
        "marques_graded_k4", "marques_graded_k8",
    };

    // Same sizes as evaluator stage 1: {10, 50, 100}
    // n=100 means 10K U entries per matrix — manageable
    std::vector<int> sizes = {10, 50, 100};

    int total = 0;
    for (auto& pat : adv_names) {
        for (int sz : sizes) {
            TestMatrix tm = make_adversarial(pat, sz);
            int n = tm.n;

            std::vector<double> sigma(n, 0.0);
            std::vector<double> U(n * n, 0.0);
            std::vector<double> V(n * n, 0.0);
            int info;

            run_hgbsvd(n, tm.d.data(), tm.e.data(),
                       sigma.data(), U.data(), V.data(), &info);

            printf("MATRIX %s n=%d info=%d\n", pat.c_str(), n, info);

            printf("SIGMA");
            for (int i = 0; i < n; i++) printf(" %.17e", sigma[i]);
            printf("\n");

            printf("U");
            for (int i = 0; i < n * n; i++) printf(" %.17e", U[i]);
            printf("\n");

            printf("V");
            for (int i = 0; i < n * n; i++) printf(" %.17e", V[i]);
            printf("\n");

            total++;
        }
    }
    fprintf(stderr, "Dumped %d matrices\n", total);

    return 0;
}
