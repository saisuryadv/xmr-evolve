/* dlaxrf_seltw_part_f77.f -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int dlaxrf_seltw_part__(integer *n, integer *ia, integer *ie,
	 integer *dir, doublereal *g, doublereal *gn, integer *omega, integer 
	*twistok, doublereal *s, doublereal *arcsq, doublereal *abeta, 
	doublereal *agmax, doublereal *agsmax, doublereal *abrcfac)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    doublereal f;
    integer i__, j;
    doublereal z__, beta, gmax, rcsq, gsmax;
    integer iprev, inext;
    doublereal invns, trpns, brcfac;

/*     IMPLICIT NONE */



/*  Purpose */
/*  ======= */

/*     Compute auxiliary quantities to evaluate element growth and */
/*     relative condition for a bidiagonal factorization, depending */
/*     on IA and IE either */
/*      (L)  a top-to-bottom one if 1 = IA <= IE < N, DIR=+1 */
/*      (U)  a bottom-to-top one if N = IA >= IE > 1, DIR=-1. */
/*     For the following assume case (L), that is, G = D and GN = LD for */
/*     an LDL (partial) factor. */

/*     The quantity of main interest is an estimate for the relative */
/*     condition number */
/*       RC := || X L' s || * || invX invL s || */
/*     for a nonsingular diagonal scaling matrix X which minimizes it. */
/*     The latter is equivalent to saying */
/*       RC := ||r||^2 with r(i) := sqrt[ (L'*s)(i) * (invL*s)(i) ] */

/*     The result arrays are changed as follows by this routine (again */
/*     only for the LDL case): */
/*        ARCSQ(i) += ||r(1:i-1)||^2,     i = 2:IE+1. */
/*        ABETA(i) += (invL*s)(i) - s(i), i = 2:IE+1. */
/*        AGMAX(i)  is maxed to Max{ |D(i)| : 1 <= j < i } */
/*        AGSMAX(i) is maxed to Max{ |D(i)*s(i)| : 1 <= j < i }. */
/*        ABRCFAC(i) gives an additional factor to multiply the final */
/*          RC estimate for twist i with. This is to account for the */
/*          fact that the standard RC-formula needs to be modified */
/*          if blocks are present. The entries are maxed as well. */

/*       CC   ARCSQ_NOS is analogous, except that we bound RC by */
/*       CC      || X *L' || * || invX * invL ||. */
/*       CC    cancelled for now: need separate beta return as well */

/*     For all cases, the entry 1 is unchanged, as the corresponding */
/*     quantity would be zero. Besides that, only entries where a twist */
/*     is allowed are touched (TWISTOK(i)!=0), this should include all */
/*     indices where a no (non-final) block ends. */

/*     Notes. */
/*     a) */
/*         (L'*s)(i) = s(i) + */
/*            k(i)*s(i+2), if i starts a block and i != n-1 */
/*            l(i)*s(i+1), otherwise and i < n, */
/*            0, i=n or i=n-1 and i starts a block. */
/*     b) */
/*         (invL*s)(i) = s(i) + */
/*            0,   if i=1 or i ends a block, */
/*           -l(i-1) * (invL*s)(i-1),   if i-1 does not end a block, */
/*           -k(i-2) * (invL*s)(i-2) - l(i-1)*s(i-1),   otherwise. */

/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Local Variables .. */



/*  ===== Executable Statements ========================================== */


    /* Parameter adjustments */
    --abrcfac;
    --agsmax;
    --agmax;
    --abeta;
    --arcsq;
    --s;
    --twistok;
    --omega;
    --gn;
    --g;

    /* Function Body */
    i__ = *ia;
    j = *ie + *dir;

    beta = 0.;
    rcsq = 0.;
    gmax = 0.;
    gsmax = 0.;
L90001:
/*        Here we have: */
/*          I != J, result quantities for I have been set and are */
/*          given by the corresponding locals. */

    inext = i__ + *dir;

    invns = beta + s[i__];
/* Computing MAX */
    d__2 = gmax, d__3 = (d__1 = g[i__], abs(d__1));
    gmax = max(d__2,d__3);
/* Computing MAX */
    d__2 = gsmax, d__3 = (d__1 = g[i__] * s[i__], abs(d__1));
    gsmax = max(d__2,d__3);

    if (omega[inext] == 0) {
	f = gn[i__] / g[i__];
	trpns = s[i__] + f * s[inext];
	rcsq += (d__1 = invns * trpns, abs(d__1));
	beta = -f * invns;
    } else if (inext == 1 || inext == *n) {
	rcsq += (d__1 = invns * s[i__], abs(d__1));
	beta = 0.;
    } else {
/* Computing 2nd power */
	d__1 = gn[i__];
	z__ = gn[inext] / (g[i__] * g[inext] - d__1 * d__1);
	f = -gn[i__] * z__;
/*           now F=k(i) */
	trpns = s[i__] + f * s[inext + *dir];
	rcsq += (d__1 = invns * trpns, abs(d__1));
/* Computing MAX */
	d__2 = gmax, d__3 = (d__1 = g[i__], abs(d__1));
	gmax = max(d__2,d__3);
/* Computing MAX */
	d__2 = gsmax, d__3 = (d__1 = g[i__] * s[i__], abs(d__1));
	gsmax = max(d__2,d__3);
	beta = -f * invns;
	iprev = i__;
	i__ = inext;
	inext += *dir;
	f = g[iprev] * z__;
/*           now F = N(iprev) */
	trpns = s[i__] + f * s[inext];
	rcsq += (d__1 = s[i__] * trpns, abs(d__1));
	beta -= f * s[i__];
    }
    if (twistok[inext] != 0) {
	arcsq[inext] += rcsq;
	abeta[inext] += beta;
/* Computing MAX */
	d__1 = agmax[inext];
	agmax[inext] = max(d__1,gmax);
/* Computing MAX */
	d__1 = agsmax[inext];
	agsmax[inext] = max(d__1,gsmax);
    }
    if (inext == j) {
	goto L90002;
    }
    i__ = inext;
    inext += *dir;
    goto L90001;
L90002:

    brcfac = 1.;
    i__ = *ia;
L90003:
/* Computing MAX */
    d__1 = abrcfac[i__];
    abrcfac[i__] = max(d__1,brcfac);
    if (i__ == *ie + *dir) {
	goto L90004;
    }
    if (omega[i__] != 0) {
/* Computing MAX */
	d__3 = brcfac, d__4 = (d__1 = g[i__ - *dir] / gn[i__ - *dir], abs(
		d__1)), d__3 = max(d__3,d__4), d__4 = (d__2 = g[i__] / gn[i__ 
		- *dir], abs(d__2));
	brcfac = max(d__3,d__4);
    }
    i__ += *dir;
    goto L90003;
L90004:

    return 0;
} /* dlaxrf_seltw_part__ */

