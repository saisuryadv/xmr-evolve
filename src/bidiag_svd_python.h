#pragma once
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>

struct BidiagSVDResult {
    std::vector<double> sigma;
    std::vector<double> U;
    std::vector<double> V;
    int info;
};

extern "C" {
    double dlamch_(const char*, int);
}

static FILE* g_py_in = nullptr;
static FILE* g_py_out = nullptr;
static pid_t g_py_pid = 0;

inline void init_python_bridge() {
    if (g_py_in) return;
    int pipe_to_py[2], pipe_from_py[2];
    pipe(pipe_to_py);
    pipe(pipe_from_py);
    g_py_pid = fork();
    if (g_py_pid == 0) {
        close(pipe_to_py[1]);
        close(pipe_from_py[0]);
        dup2(pipe_to_py[0], STDIN_FILENO);
        dup2(pipe_from_py[1], STDOUT_FILENO);
        close(pipe_to_py[0]);
        close(pipe_from_py[1]);
        execlp("python3", "python3", "/home/claude/bsvd_server.py", (char*)nullptr);
        _exit(1);
    }
    close(pipe_to_py[0]);
    close(pipe_from_py[1]);
    g_py_in = fdopen(pipe_to_py[1], "wb");
    g_py_out = fdopen(pipe_from_py[0], "rb");
}

inline BidiagSVDResult bidiag_svd(int n, const double* d, const double* e) {
    BidiagSVDResult result;
    result.sigma.resize(n);
    result.U.resize(n * n, 0.0);
    result.V.resize(n * n, 0.0);
    result.info = 0;
    if (n <= 0) return result;
    
    init_python_bridge();
    
    fwrite(&n, sizeof(int), 1, g_py_in);
    fwrite(d, sizeof(double), n, g_py_in);
    int ne = n > 1 ? n - 1 : 0;
    if (ne > 0) fwrite(e, sizeof(double), ne, g_py_in);
    fflush(g_py_in);
    
    fread(result.sigma.data(), sizeof(double), n, g_py_out);
    
    std::vector<double> tmp(n * n);
    fread(tmp.data(), sizeof(double), n * n, g_py_out);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            result.U[k * n + i] = tmp[i * n + k];
    
    fread(tmp.data(), sizeof(double), n * n, g_py_out);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            result.V[k * n + i] = tmp[i * n + k];
    
    return result;
}
