/* Minimal f2c runtime stubs — replaces libf2c.a
   Only the functions actually referenced by our CLAPACK/XMR/hgbsvd code. */

#include <math.h>
#include <string.h>
#include "f2c.h"

/* d_sign: return |a| * sign(b) */
double d_sign(doublereal *a, doublereal *b)
{
    double x = (*a >= 0 ? *a : -*a);
    return *b >= 0 ? x : -x;
}

/* pow_dd: a^b (double, double) */
double pow_dd(doublereal *ap, doublereal *bp)
{
    return pow(*ap, *bp);
}

/* pow_di: a^i (double, integer) */
double pow_di(doublereal *ap, integer *bp)
{
    double pow_, x;
    integer n;
    unsigned long u;

    pow_ = 1;
    x = *ap;
    n = *bp;

    if (n < 0) {
        n = -n;
        x = 1.0 / x;
    }
    u = (unsigned long)n;
    for ( ; u; u >>= 1) {
        if (u & 1)
            pow_ *= x;
        x *= x;
    }
    return pow_;
}

/* d_lg10: log10(x) */
double d_lg10(doublereal *x)
{
    return log10(*x);
}

/* i_nint: nearest integer */
integer i_nint(real *x)
{
    return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

/* i_dnnt: nearest integer (double) */
integer i_dnnt(doublereal *x)
{
    return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

/* d_int: truncation to double */
double d_int(doublereal *x)
{
    return (*x > 0) ? floor(*x) : -floor(-*x);
}

/* s_cmp: compare strings a[0..la-1] and b[0..lb-1], blank-padded */
integer s_cmp(const char *a, const char *b, ftnlen la, ftnlen lb)
{
    const char *aend, *bend;
    aend = a + la;
    bend = b + lb;

    /* compare common part */
    if (la <= lb) {
        while (a < aend) {
            if (*a != *b) return (*a < *b ? -1 : 1);
            ++a; ++b;
        }
        /* a exhausted, check if b has trailing non-blanks */
        while (b < bend) {
            if (*b != ' ') return (-1);
            ++b;
        }
    } else {
        while (b < bend) {
            if (*a != *b) return (*a < *b ? -1 : 1);
            ++a; ++b;
        }
        while (a < aend) {
            if (*a != ' ') return (1);
            ++a;
        }
    }
    return 0;
}

/* s_copy: copy string b to a, blank-pad if la > lb */
void s_copy(char *a, const char *b, ftnlen la, ftnlen lb)
{
    char *aend;
    const char *bend;

    aend = a + la;
    if (la <= lb) {
        while (a < aend)
            *a++ = *b++;
    } else {
        bend = b + lb;
        while (b < bend)
            *a++ = *b++;
        while (a < aend)
            *a++ = ' ';
    }
}

/* i_len: length of string (Fortran LEN intrinsic) */
integer i_len(const char *s, ftnlen n)
{
    return n;
}

/* i_sign: sign transfer for integers (Fortran ISIGN) */
integer i_sign(integer *a, integer *b)
{
    integer x = *a >= 0 ? *a : -*a;
    return *b >= 0 ? x : -x;
}

/* pow_ii: integer^integer */
integer pow_ii(integer *ap, integer *bp)
{
    integer pow_, x, n;
    unsigned u;

    x = *ap;
    n = *bp;

    if (n <= 0) {
        if (n == 0 || x == 1) return 1;
        if (x != -1) return x == 0 ? 1 : 0;
        n = -n;
    }
    u = (unsigned)n;
    pow_ = 1;
    for ( ; ; ) {
        if (u & 1) pow_ *= x;
        u >>= 1;
        if (!u) break;
        x *= x;
    }
    return pow_;
}

/* I/O stubs — used by hgbsvd debug prints and xerbla.
   We don't need Fortran I/O; these are no-ops. */
integer s_wsfe(cilist *a) { return 0; }
integer e_wsfe(void) { return 0; }
integer do_fio(integer *number, char *ptr, ftnlen len) { return 0; }
integer s_wsle(cilist *a) { return 0; }
integer e_wsle(void) { return 0; }
integer do_lio(integer *type, integer *number, char *ptr, ftnlen len) { return 0; }
