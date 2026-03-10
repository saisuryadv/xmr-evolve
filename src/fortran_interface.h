#pragma once
// C++ interface to Fortran routines compiled in libxmrlapack.a
// All Fortran routines use lowercase names with trailing underscore (gfortran default)
// Arrays are column-major (Fortran order)

#include <cstdint>

extern "C" {

// ============================================================
// LAPACK MR³ eigensolver chain
// ============================================================

// DSTEMR: MR³ tridiagonal eigensolver (the core O(n²) solver)
void dstemr_(const char* jobz, const char* range,
             const int* n, double* d, double* e,
             const double* vl, const double* vu,
             const int* il, const int* iu,
             int* m, double* w, double* z, const int* ldz,
             const int* nzc, int* isuppz,
             int* tryrac, double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);

// DLARRE: Find eigenvalues of symmetric tridiagonal (used by DSTEMR)
void dlarre_(const char* range, const int* n,
             double* vl, double* vu,
             const int* il, const int* iu,
             double* d, double* e, double* e2,
             const double* rtol1, const double* rtol2,
             const double* spltol, int* nsplit,
             int* isplit, int* m, double* w,
             double* werr, double* wgap, int* iblock,
             int* indexw, double* gers,
             double* pivmin, double* work, int* iwork, int* info);

// DLARRV: Compute eigenvectors from eigenvalues (used by DSTEMR)
void dlarrv_(const int* n, const double* vl, const double* vu,
             double* d, double* l, const double* pivmin,
             const int* isplit, const int* m,
             const int* dol, const int* dou,
             const double* minrgp, const double* rtol1, const double* rtol2,
             double* w, double* werr, double* wgap,
             const int* iblock, const int* indexw, const double* gers,
             double* z, const int* ldz, int* isuppz,
             double* work, int* iwork, int* info);

// DLAR1V: Scaled inverse iteration vector
void dlar1v_(const int* n, const int* b1, const int* bn,
             const double* lambda, const double* d, const double* l,
             const double* ld, const double* lld,
             const double* pivmin, const double* gaptol,
             double* z, const int* wantnc, int* negcnt,
             double* ztz, double* mingma, int* r,
             int* isuppz, double* nrminv, double* resid,
             double* rqcorr, double* work);

// DLARRF: Find child representation
void dlarrf_(const int* n, const double* d, const double* l, const double* ld,
             const int* clstrt, const int* clend,
             const double* w, double* wgap, const double* werr,
             const double* spdiam, const double* clgapl, const double* clgapr,
             const double* pivmin,
             double* sigma, double* dplus, double* lplus, double* work,
             int* info);

// ============================================================
// LAPACK Bidiagonal SVD drivers
// ============================================================

// DBDSQR: QR iteration bidiagonal SVD
void dbdsqr_(const char* uplo, const int* n, const int* ncvt,
             const int* nru, const int* ncc,
             double* d, double* e,
             double* vt, const int* ldvt,
             double* u, const int* ldu,
             double* c, const int* ldc,
             double* work, int* info);

// DBDSVDX: Bisection-based bidiagonal SVD via TGK
void dbdsvdx_(const char* uplo, const char* jobz, const char* range,
              const int* n, const double* d, const double* e,
              const double* vl, const double* vu,
              const int* il, const int* iu,
              int* ns, double* s, double* z, const int* ldz,
              double* work, int* iwork, int* info);

// DBDSDC: Divide-and-conquer bidiagonal SVD
void dbdsdc_(const char* uplo, const char* compq, const int* n,
             double* d, double* e, double* u, const int* ldu,
             double* vt, const int* ldvt, double* q, int* iq,
             double* work, int* iwork, int* info);

// ============================================================
// XMR (Willems' MR³ implementation)
// ============================================================

// DSTEXR: Willems' XMR extended MR³ tridiagonal eigensolver
// NOTE: XMR interface differs from DSTEMR — no jobz/range/vl/vu/m params,
// uses wil/wiu for index range directly
void dstexr_(int* n, double* d, double* e,
             int* wil, int* wiu,
             double* w, double* z, int* ldz, int* isuppz,
             double* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

// ============================================================
// hgbsvd (Großer-Lang bidiagonal SVD)
// ============================================================

// DBDSGR: Großer-Lang bidiagonal SVD driver (DSTEMR-like interface)
// JOBZ: 'N'=values only, 'L'=left, 'R'=right, 'A'=both
// RANGE: 'A'=all, 'V'=value range, 'I'=index range (only 'A' supported)
// For JOBZ='A': V in Z rows 1:N, U in Z rows LDZ/2+1:LDZ/2+N
void dbdsgr_(const char* jobz, const char* range,
             const int* n, double* d, double* e,
             const double* vl, const double* vu,
             const int* il, const int* iu,
             const double* abstol, int* m, double* w,
             double* z, const int* ldz, int* isuppz,
             double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);

// DLARRI: Eigenvalue refinement for bidiagonal SVD
void dlarri_(const int* n, const double* d, double* e, double* e2,
             const double* pivmin, int* nsplit, int* isplit,
             int* m, double* w, double* werr,
             const double* rtol1, const double* rtol2,
             int* iblock, int* indexw, double* work, int* iwork,
             int* info);

// ============================================================
// LAPACK utility routines
// ============================================================

double dlamch_(const char* cmach);
double dlapy2_(const double* x, const double* y);
double dnrm2_(const int* n, const double* x, const int* incx);
double ddot_(const int* n, const double* x, const int* incx,
             const double* y, const int* incy);
void dcopy_(const int* n, const double* x, const int* incx,
            double* y, const int* incy);
void dscal_(const int* n, const double* a, double* x, const int* incx);
void daxpy_(const int* n, const double* a, const double* x, const int* incx,
            double* y, const int* incy);
void dswap_(const int* n, double* x, const int* incx, double* y, const int* incy);
void dlacpy_(const char* uplo, const int* m, const int* n,
             const double* a, const int* lda, double* b, const int* ldb);
void dlaset_(const char* uplo, const int* m, const int* n,
             const double* alpha, const double* beta, double* a, const int* lda);
void dlasrt_(const char* id, const int* n, double* d, int* info);
int ilaenv_(const int* ispec, const char* name, const char* opts,
            const int* n1, const int* n2, const int* n3, const int* n4);
void xerbla_(const char* srname, const int* info);

} // extern "C"
