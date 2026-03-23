/*
 * C wrapper for XMR block factorization routines.
 * Provides clean interfaces callable from Python via ctypes.
 */

#include <string.h>
#include <stdlib.h>

/* Fortran COMMON block /XMRSTATS/ - must match the Fortran declaration exactly */
struct {
    double xtime1, xtime2, xtime3;
    double xddddd;
    int xnblcks;
    int xnnodes, xmaxdepth, xnumfn, xnumft, xnumgv, xnumgv0;
    int xnumfs_2, xnumfs_k;
    int xnbis_init;
    int xnbis_cob, xnbis_iib, xnbis_class, xnbis_sng, xnbis_clb;
    int xnbis_wasted;
    int xnrqi, xnrqibis;
    int xmaxnrqi, xmaxnrqibis;
    int xnenvgv, xnenvtf;
    int xiiiii;
    int xstealthmode;  /* Fortran LOGICAL = int */
} xmrstats_;

/* Fortran subroutine declarations (name mangling: lowercase + underscore) */
extern void dlaxrr_(int *n, int *k, int *type, double *e, double *pivbase,
                    double *repr, int *repi);

extern void dlaxrs_(int *n, int *ia, int *ie, double *repr, int *repi,
                    double *tau, double *shfprt,
                    double *dplus, int *omgadp, double *rplus, int *omgarp,
                    double *gammap, int *twistok,
                    double *rwork);

extern void dlaxrt_(int *n, int *j1, int *j2, double *repr, int *repi,
                    double *tau,
                    double *dplus, double *rplus, double *gammap, int *twistok,
                    double *rwork);

extern int dlaxrn_(int *n, double *repr, int *repi, double *tau);

/* Initialize the COMMON block */
void xmr_init(void) {
    memset(&xmrstats_, 0, sizeof(xmrstats_));
    xmrstats_.xstealthmode = 0;  /* FALSE */
}

/*
 * Build representation from G, OMEGA, E arrays.
 *
 * n: matrix size
 * k: twist index (1-based)
 * g[1:n]: pivots (1-based, so g[0] unused)
 * omega[0:n+1]: block flags
 * e[1:n-1]: off-diagonals (1-based)
 * pivbase: pivot threshold
 *
 * Output:
 *   repr[4*n+3]: representation data
 *   repi[6+n+n/2]: integer representation data
 */
void xmr_build_repr(int n, int k, double *g, int *omega, double *e,
                     double pivbase, double *repr, int *repi) {
    int type = 0;
    int i;
    /* IYOMGA=4 in Fortran 1-based = index 3 in C 0-based */
    int iyomga_c = 3;

    /* Copy G into REPR(1:N) → C repr[0:N-1] */
    for (i = 0; i < n; i++) {
        repr[i] = g[i + 1];  /* g is 1-based */
    }

    /* Copy OMEGA into REPI(IYOMGA:IYOMGA+N+1) → C repi[3:3+N+1] */
    for (i = 0; i <= n + 1; i++) {
        repi[iyomga_c + i] = omega[i];
    }

    dlaxrr_(&n, &k, &type, &e[1], &pivbase, repr, repi);
}

/*
 * Shift transformation with blocks.
 *
 * All arrays are allocated by the caller (Python).
 * repr/repi: from xmr_build_repr
 * tau: shift value
 * shfprt[1:n]: shift partition (usually all 1.0)
 *
 * Output:
 *   dplus[1:n-1], rplus[2:n]: shifted pivots
 *   omgadp[1:n], omgarp[1:n]: shifted block flags
 *   gammap[ia:ie]: gamma values
 *   twistok[ia:ie]: valid twist flags
 */
void xmr_shift(int n, int ia, int ie, double *repr, int *repi,
               double tau, double *shfprt,
               double *dplus, int *omgadp, double *rplus, int *omgarp,
               double *gammap, int *twistok, double *rwork) {
    dlaxrs_(&n, &ia, &ie, repr, repi,
            &tau, shfprt,
            dplus, omgadp, rplus, omgarp,
            gammap, twistok, rwork);
}

/*
 * Twisted factorization for eigenvector extraction.
 *
 * repr/repi: representation
 * tau: eigenvalue approximation
 *
 * Output:
 *   dplus[1:n-1], rplus[2:n]: twisted pivots
 *   gammap[j1:j2]: gamma values
 *   twistok[j1:j2]: valid twist flags
 */
void xmr_twist(int n, int j1, int j2, double *repr, int *repi,
               double tau,
               double *dplus, double *rplus, double *gammap, int *twistok,
               double *rwork) {
    dlaxrt_(&n, &j1, &j2, repr, repi,
            &tau,
            dplus, rplus, gammap, twistok,
            rwork);
}

/*
 * Inertia count using the representation (block-aware).
 * Returns xi = 2*negcount + issing.
 *   xi even: tau is in interval (ew xi/2, ew xi/2+1)
 *   xi odd:  tau is on eigenvalue (xi+1)/2
 */
int xmr_negcount(int n, double *repr, int *repi, double tau) {
    return dlaxrn_(&n, repr, repi, &tau);
}

/*
 * Find child representation via dlaxrf.
 *
 * father_repr/father_repi: parent representation
 * e: off-diagonals (1-based, e[1..N-1])
 * icbeg, icend: cluster eigenvalue index range (1-based)
 * ewl_lu: eigenvalue bounds [2*icbeg-1 .. 2*icend] (interleaved lower/upper)
 * lgap, ugap: gaps to left and right of cluster
 * gaptol, spdiam: parameters
 * taubar: accumulated shift so far
 *
 * Output:
 *   son_repr/son_repi: child representation
 *   tau_out: chosen shift
 *   returns INFO (0=success)
 */
extern void dlaxrf_(int *n, double *e, double *gaptol, double *spdiam,
                    int *depth, double *taubar, int *icbeg, int *icend,
                    double *frepr, int *frepi, double *flgap, double *fugap,
                    int *fewl_ae, double *fewl_lu, int *frginfo,
                    int *dogrt, double *frepr_grt, double *shfmod, double *env, int *mode,
                    double *srepr, int *srepi, double *slgap, double *sugap,
                    int *sewl_ae, double *sewl_lu, double *tau,
                    double *rwork, int *lrwork, int *iwork, int *liwork, int *info);

extern void wsreq_xrf_(int *n, int *icbeg, int *icend, int *reqr, int *reqi);

int xmr_find_child_repr(
    int n, double *e, double *father_repr, int *father_repi,
    int icbeg, int icend, double *evals, double lgap, double ugap,
    double gaptol, double spdiam, double taubar,
    double *son_repr, int *son_repi, double *tau_out)
{
    int iclen = icend - icbeg + 1;
    int ewl_size = 2 * icend;  /* indices 2*icbeg-1 .. 2*icend */

    /* Allocate eigenvalue list arrays (1-indexed in Fortran) */
    int *fewl_ae = (int *)calloc(ewl_size + 1, sizeof(int));
    double *fewl_lu = (double *)calloc(ewl_size + 1, sizeof(double));
    int *sewl_ae = (int *)calloc(ewl_size + 1, sizeof(int));
    double *sewl_lu = (double *)calloc(ewl_size + 1, sizeof(double));

    /* FRGINFO: gap info for each eigenvalue boundary.
       FRGINFO(i) for i=icbeg-1..icend.
       Value 2 = full gap (at cluster boundaries)
       Value 0 = no gap (inside cluster)
       Array index 0 = FRGINFO(icbeg-1) = left boundary
       Array index iclen = FRGINFO(icend) = right boundary */
    int frginfo_size = iclen + 2;
    int *frginfo = (int *)calloc(frginfo_size, sizeof(int));
    frginfo[0] = 2;           /* left boundary: full gap */
    for (int i = 1; i < iclen; i++) frginfo[i] = 0;  /* internal: no gap */
    frginfo[iclen] = 2;       /* right boundary: full gap */

    /* Set up eigenvalue list: each eigenvalue i has bounds [EWL_LU(2i-1), EWL_LU(2i)] */
    /* EWL_AE tracks which interval an eigenvalue belongs to */
    /* For a cluster, all eigenvalues initially share one interval: [icbeg, icend] */
    /* But individual bounds should bracket each eigenvalue tightly */
    for (int j = 0; j < iclen; j++) {
        int i = icbeg + j;  /* 1-based eigenvalue index */
        double ev = evals[j];
        /* Tight bounds: ±n*eps*|ev| or ±eps for tiny eigenvalues */
        double err = n * 2.22e-16 * (ev > 0 ? ev : -ev);
        if (err < 2.22e-16) err = 2.22e-16;
        fewl_lu[2*i - 2] = ev - err;  /* lower bound */
        fewl_lu[2*i - 1] = ev + err;  /* upper bound */
        /* All eigenvalues in cluster share the same interval [icbeg, icend] */
        fewl_ae[2*i - 2] = icbeg;
        fewl_ae[2*i - 1] = icend;
    }

    /* SHFMOD: all ones (no shift perturbation) */
    double *shfmod = (double *)calloc(3 * n, sizeof(double));
    for (int i = 0; i < 3 * n; i++) shfmod[i] = 1.0;

    /* ENV: envelope, set to spdiam (permissive) */
    double *env = (double *)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) env[i] = spdiam;

    /* Workspace query */
    int reqr = -1, reqi = -1;
    wsreq_xrf_(&n, &icbeg, &icend, &reqr, &reqi);

    double *rwork = (double *)calloc(reqr + 1, sizeof(double));
    int *iwork = (int *)calloc(reqi + 1, sizeof(int));

    int depth = 0;
    int dogrt = 0;  /* no gap retries */
    int mode = 1;   /* try both outside and inside */
    int info = 0;
    double slgap = 0.0, sugap = 0.0;

    /* frepr_grt not used when dogrt=0, but must point somewhere valid */
    double *frepr_grt = father_repr;

    dlaxrf_(&n, &e[1], &gaptol, &spdiam, &depth, &taubar, &icbeg, &icend,
            father_repr, father_repi,
            &lgap, &ugap,
            &fewl_ae[2*icbeg - 2], &fewl_lu[2*icbeg - 2],
            frginfo,
            &dogrt, frepr_grt, shfmod, env, &mode,
            son_repr, son_repi, &slgap, &sugap,
            &sewl_ae[2*icbeg - 2], &sewl_lu[2*icbeg - 2], tau_out,
            rwork, &reqr, iwork, &reqi, &info);

    free(fewl_ae); free(fewl_lu);
    free(sewl_ae); free(sewl_lu);
    free(frginfo); free(shfmod); free(env);
    free(rwork); free(iwork);

    return info;
}

/*
 * Full XMR eigenvector computation via dlaxrv.
 *
 * Computes eigenpairs wil:wiu of the tridiagonal given by root representation.
 *
 * n: matrix dimension
 * e: off-diagonals, 1-based e[1..n-1] (e[0] unused)
 * rootr/rooti: root representation
 * evals_in: initial eigenvalue approximations, 0-based [0..n-1] for ALL n eigenvalues
 * wil, wiu: wanted eigenvalue range (1-based)
 * spdiam, gaptol: parameters
 *
 * Output:
 *   w_out: computed eigenvalues [0..wiu-wil], 0-based
 *   z_out: eigenvectors column-major [n * (wiu-wil+1)], 0-based
 *   returns INFO
 */
extern void dlaxrv_(int *n, double *e, double *rootr, int *rooti,
                    int *ewl_ae, double *ewl_lu,
                    int *wil, int *wiu, double *spdiam, double *gaptol,
                    double *w, double *z, int *ldz, int *isuppz,
                    double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
extern void wsreq_xrv_(int *n, int *reqr, int *reqi);

extern void dlaxre_(int *n, double *d, double *e, double *gl, double *gu,
                    double *abserr, double *gaptol,
                    double *repr, int *repi, double *tau,
                    int *ewl_ae, double *ewl_lu, char *emode,
                    double *rwork, int *info, int emode_len);

int xmr_eigenvectors(
    int n, double *e_1based, double *rootr, int *rooti,
    double *evals_in, int wil, int wiu,
    double spdiam, double gaptol,
    double *w_out, double *z_out)
{
    int nwant = wiu - wil + 1;

    /* Step 1: Build root representation and EWL using dlaxre.
       dlaxre takes the tridiagonal entries D(1:N) and E(1:N) where
       E(N) should be 0. It builds a proper root representation and
       initializes the eigenvalue list. */

    double *d_arr = (double *)calloc(n, sizeof(double));
    double *e_arr = (double *)calloc(n, sizeof(double));
    /* D = diagonal (zeros for T_GK), E = off-diag + trailing zero */
    for (int i = 0; i < n; i++) d_arr[i] = 0.0;  /* T_GK has zero diagonal */
    for (int i = 0; i < n - 1; i++) e_arr[i] = e_1based[i + 1];
    e_arr[n - 1] = 0.0;

    /* Gerschgorin bounds */
    double gl = evals_in[0] - 0.01 * (evals_in[0] > 0 ? evals_in[0] : -evals_in[0]) - 1e-20;
    double gu = evals_in[n-1] + 0.01 * (evals_in[n-1] > 0 ? evals_in[n-1] : -evals_in[n-1]) + 1e-20;
    double abserr = 0.0;

    double *repr = (double *)calloc(4 * n + 3, sizeof(double));
    int *repi = (int *)calloc(6 + n + n / 2, sizeof(int));
    double tau_re = 0.0;

    int *ewl_ae = (int *)calloc(2 * n, sizeof(int));
    double *ewl_lu = (double *)calloc(2 * n, sizeof(double));
    double *rwork_re = (double *)calloc(6 * n + 10, sizeof(double));
    char emode = 'o';  /* Gerschgorin-based initialization */
    int info_re = 0;

    dlaxre_(&n, d_arr, e_arr, &gl, &gu, &abserr, &gaptol,
            repr, repi, &tau_re,
            ewl_ae, ewl_lu, &emode,
            rwork_re, &info_re, 1);

    free(d_arr); free(rwork_re);

    if (info_re != 0) {
        free(e_arr); free(repr); free(repi);
        free(ewl_ae); free(ewl_lu);
        return -10 - info_re;
    }

    /* Step 2: Call dlaxrv with the proper root representation and EWL */
    int reqr = -1, reqi = -1;
    wsreq_xrv_(&n, &reqr, &reqi);

    double *rwork = (double *)calloc(reqr + 1, sizeof(double));
    int *iwork = (int *)calloc(reqi + 1, sizeof(int));

    double *w = (double *)calloc(nwant, sizeof(double));
    double *z = (double *)calloc((size_t)n * nwant, sizeof(double));
    int *isuppz = (int *)calloc(2 * nwant, sizeof(int));
    int ldz = n;
    int info = 0;

    dlaxrv_(&n, e_arr, repr, repi,
            ewl_ae, ewl_lu,
            &wil, &wiu, &spdiam, &gaptol,
            w, z, &ldz, isuppz,
            rwork, &reqr, iwork, &reqi, &info);

    /* Copy results and add back the dlaxre shift */
    for (int j = 0; j < nwant; j++) {
        w_out[j] = w[j] + tau_re;
        for (int i = 0; i < n; i++) {
            z_out[i + n * j] = z[i + n * j];
        }
    }

    free(e_arr); free(repr); free(repi);
    free(ewl_ae); free(ewl_lu);
    free(rwork); free(iwork); free(isuppz);
    free(w); free(z);

    return info;
}
