/* Real implementations for XMR computational kernels.
   Hand-translated from Willems' Fortran 90 originals.
   These replace the no-op stubs that were causing incorrect eigenvalues. */

#include "f2c.h"
#include <math.h>

/* ======================================================================
   DLAXRT_STAT: Stationary twisted factorization

   Computes stationary factorization from 1 to J or from N down to J.

   I0=1 and I1=N-1:
     assumes that the twist of the source is >= J
     produces GPLUS(1),...,GPLUS(J-1) and S(1),...,S(J)
     --> only S(1)=0 if J=1

   I0=2 and I1=N:
     assumes that the twist of the source is <= J
     produces GPLUS(J+1),...,GPLUS(N) and S(J),...,S(N)
     --> only S(N)=0 if J=N
   ====================================================================== */
int dlaxrt_stat__(integer *n, integer *j, integer *i0, integer *i1,
    doublereal *g, integer *omega, doublereal *gnsq, doublereal *ngn,
    doublereal *bdet, doublereal *pivmin, doublereal *tau,
    doublereal *gplus, doublereal *s)
{
    /* PRBRK: Determines branch to choose for breaking a block. Default: 8 */
    const integer PRBRK = 8;

    integer i, iprev, inext, dir;
    doublereal sbound, smt, oldsmt, x;

    /* Parameter adjustments: Fortran 1-based to C 0-based */
    --g; --gnsq; --ngn; --bdet; --s;
    /* omega is 0:N+1 indexed, no adjustment needed since caller passes &repi[iyomga] */
    /* gplus is I0:I1 indexed — adjust by I0 not 1 */
    gplus -= *i0;

    if (*i0 == 1 && *i1 == *n - 1) {
        dir = 1;
        i = 1;
    } else if (*i0 == 2 && *i1 == *n) {
        dir = -1;
        i = *n;
    } else {
        dir = 0;
        i = -1;
    }
    inext = i + dir;

    sbound = (PRBRK < *n ? PRBRK : *n) * fabs(*tau);
    s[i] = 0.0;

    for (;;) {
        /* Inner loop: process until we hit J or a block boundary (omega != 0) */
        for (;;) {
            if (i == *j || omega[inext] != 0) break;
            smt = s[i] - *tau;
            gplus[i] = g[i] + smt;
            if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
            s[inext] = ngn[i] * (smt / gplus[i]);
            i = inext;
            inext = inext + dir;
        }
        if (i == *j) break;

        /* Block boundary: break the block */
        smt = s[i] - *tau;
        gplus[i] = g[i] + smt;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
        s[inext] = -gnsq[i] / gplus[i];
        if (inext == *j) break;

        iprev = i;
        i = inext;
        inext = inext + dir;

        oldsmt = smt;
        smt = s[i] - *tau;
        gplus[i] = g[i] + smt;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);

        /* Compute next s after a block was broken
           Note: Recall that NGN(IPREV) is undefined */
        if (g[iprev] == 0.0) {
            s[inext] = -gnsq[i] / gplus[i];
        } else {
            if (fabs(s[i]) <= sbound ||
                (g[iprev] >= 0.0 ? 1.0 : -1.0) != (gplus[iprev] >= 0.0 ? 1.0 : -1.0))
            {
                x = g[iprev] * smt + gnsq[iprev];
            } else {
                x = -s[i] * oldsmt - g[iprev] * (*tau);
            }
            s[inext] = (gnsq[i] * x) / (bdet[iprev] * gplus[i]);
        }
        i = inext;
        inext = inext + dir;
    }

    return 0;
}


/* ======================================================================
   DLAXRT_PROG: Progressive twisted factorization

   Computes a progressive factorization from J1 to J2.

   J1 < J2 and I0 = 1 and I1 = N-1:
     assumes J1 is twist in source
     produces GPLUS(J1),...,GPLUS(J2-1) and P(J1+1),...,P(J2)

   J1 > J2 and I0 = 2 and I1 = N:
     assumes J2 is twist in source
     produces GPLUS(J1),...,GPLUS(J2+1) and P(J1-1),...,P(J2)
   ====================================================================== */
int dlaxrt_prog__(integer *n, integer *j1, integer *j2,
    integer *i0, integer *i1, doublereal *g, integer *omega,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, doublereal *off,
    doublereal *gplus, doublereal *p)
{
    const integer PRBRK = 8;

    integer i, iprev, inext, dir, j;
    doublereal pbound, aux;

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet; --p;
    /* gplus is I0:I1 indexed — adjust by I0 not 1 */
    gplus -= *i0;

    i = *j1;
    j = *j2;

    if (*j1 < *j2 && *i0 == 1 && *i1 == *n - 1) {
        dir = 1;
    } else if (*j2 < *j1 && *i0 == 2 && *i1 == *n) {
        dir = -1;
    } else {
        i = -1;
        dir = 0;
    }

    pbound = (PRBRK < *n ? PRBRK : *n) * fabs(*tau);
    inext = i + dir;

    aux = *off + (g[i] - *tau);
    if (omega[i] == 1) {
        gplus[i] = aux;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
        p[inext] = g[inext] - gnsq[inext] / gplus[i];
        i = inext;
        inext = inext + dir;
        aux = p[i] - *tau;
    }

    for (;;) {
        /* Inner loop */
        for (;;) {
            if (i == j || omega[inext] != 0) break;
            gplus[i] = ngn[inext] + aux;
            if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
            p[inext] = (g[inext] / gplus[i]) * aux;
            i = inext;
            inext = inext + dir;
            aux = p[i] - *tau;
        }
        if (i == j) break;

        /* Block boundary */
        gplus[i] = ngn[inext] + aux;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
        p[inext] = g[inext] - gnsq[inext] / gplus[i];
        if (inext == j) break;

        iprev = i;
        i = inext;
        inext = inext + dir;

        gplus[i] = p[i] - *tau;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);

        if (fabs(p[i]) <= pbound || g[inext] == 0.0) {
            p[inext] = g[inext] - gnsq[inext] / gplus[i];
        } else {
            p[inext] = (bdet[inext] * (aux / gplus[iprev])
                        - g[inext] * (*tau))
                       / gplus[i];
        }
        i = inext;
        inext = inext + dir;
        aux = p[i] - *tau;
    }

    return 0;
}


/* ======================================================================
   DLAXRN_STAT: Negcount computation kernel (stationary)

   Used by dlaxrn_ for the tau != 0 case to count negative pivots
   in a twisted factorization.
   ====================================================================== */
int dlaxrn_stat__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbegk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *ncount,
    doublereal *s_out)
{
    const integer PRBRK = 8;

    doublereal sbound, negpiv;
    doublereal smt, gp, oldgp, oldsmt, x, s;
    integer i, iprev;

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet;
    /* lbbegk: Fortran LBBEGK(*) is 1-based, but the C pointer starts at element 0.
       We decrement so that lbbegk[IBB] (1-based) maps correctly. */
    --lbbegk;

    sbound = (PRBRK < *n ? PRBRK : *n) * fabs(*tau);
    negpiv = -(*pivmin);

    if (*dir == 1) {
        i = 1;
    } else {
        i = *n;
    }

    s = 0.0;
    *ncount = 0;

    for (;;) {
        /* Inner loop: process until we hit a block boundary */
        for (;;) {
            if (i == lbbegk[*ibb]) break;
            smt = s - *tau;
            gp = g[i] + smt;
            if (fabs(gp) < *pivmin) gp = negpiv;
            if (gp < 0.0) ++(*ncount);
            s = ngn[i] * (smt / gp);
            i = i + *dir;
        }
        if (i == *k) break;

        /* Block boundary */
        *ibb = *ibb + *dir;

        smt = s - *tau;
        gp = g[i] + smt;
        if (fabs(gp) < *pivmin) gp = negpiv;
        if (gp < 0.0) ++(*ncount);
        s = -gnsq[i] / gp;
        oldgp = gp;
        oldsmt = smt;
        iprev = i;
        i = i + *dir;
        if (i == *k) break;

        smt = s - *tau;
        gp = g[i] + smt;
        if (fabs(gp) < *pivmin) gp = negpiv;
        if (gp < 0.0) ++(*ncount);

        /* Compute next s after a block was broken
           Notes:
           (1) Recall that NGN(I-DIR) is undefined
           (2) SIGN(A,B) = Abs(A) if B>=0, -Abs(A) ow */
        if (g[iprev] == 0.0) {
            s = -gnsq[i] / gp;
        } else {
            doublereal sgngip = g[iprev] >= 0.0 ? 1.0 : -1.0;
            if (fabs(s) <= sbound || sgngip != (oldgp >= 0.0 ? 1.0 : -1.0)) {
                x = g[iprev] * smt + gnsq[iprev];
            } else {
                x = -s * oldsmt - g[iprev] * (*tau);
            }
            s = (gnsq[i] * x) / (bdet[iprev] * gp);
        }
        i = i + *dir;
    }

    *s_out = s;
    return 0;
}


/* ======================================================================
   DLAXRS_STAT: Shift selection with stationary factorization

   Computes stationary factorization with shift proportionality (SHFPRT).
   Much more complex than DLAXRT_STAT due to block creation/breaking logic.
   ====================================================================== */
int dlaxrs_stat__(integer *n, integer *j, integer *i0, integer *i1,
    doublereal *g, integer *omega, doublereal *gnsq, doublereal *ngn,
    doublereal *bdet, doublereal *pivmin, doublereal *tau,
    doublereal *shfprt, doublereal *gplus, integer *oplus,
    doublereal *s)
{
    /* Parameters from Fortran */
    const doublereal KBLOCK = 1.0 / 8.0;
    const doublereal KBLMOD = KBLOCK * 0.999;
    const doublereal PK1 = KBLMOD;
    const doublereal PK2 = KBLMOD / 3.01;
    const doublereal PRBRK = 5.0;
    const doublereal PRCRE = 5.0;
    const doublereal PROSQ = 0.25;
    const integer PMAXOSEQLEN = 1;
    const logical USEBLOCKS = TRUE_;

    doublereal abstau, oldsmt, smt, aux, cplus, bdetp, x;
    integer i, dir, iprev, inext, bcount;
    logical bres;

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet;
    --shfprt;
    /* gplus is I0:I1 indexed — adjust by I0 not 1 */
    gplus -= *i0;
    --oplus;
    --s;

    abstau = fabs(*tau);

    if (*i0 == 1 && *i1 == *n - 1) {
        dir = 1;
        i = 1;
    } else if (*i0 == 2 && *i1 == *n) {
        dir = -1;
        i = *n;
    } else {
        dir = 0;
        i = -1;
    }
    inext = i + dir;

    /* Main Loop */
    aux = 0.0;
    s[i] = 0.0;
    oplus[i] = 0;

    for (;;) {
        if (i == *j) break;

        smt = aux - (*tau) * shfprt[i];
        gplus[i] = g[i] + smt;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = 0.0;
        oplus[inext] = 0;

        if (omega[inext] != 0) {
            /* Check if the block should be kept */
            cplus = g[inext] - (*tau) * shfprt[inext];
            bres = (fabs(gplus[i] * cplus) < PK1 * gnsq[i]);
            if (bres && (inext != *j || *j == *n || *j == 1)) {
                /* -- Keep a block -- */
                aux = 0.0;
                oplus[inext] = 1;
                if (inext != *j) {
                    s[inext] = 0.0;
                    gplus[inext] = cplus;

                    bdetp = gplus[i] * cplus - gnsq[i];
                    if (g[i] == 0.0) {
                        aux = -(gnsq[inext] / bdetp) * gplus[i];
                    } else {
                        /* recall that NGN(i) is not defined YET */
                        x = (gnsq[i] / g[i]) * smt
                            - gplus[i] * ((*tau) * shfprt[inext]);
                        aux = ngn[inext] * (x / bdetp);
                    }
                    i = inext;
                    inext = inext + dir;
                    oplus[inext] = 0;
                }
            } else {
                /* -- Initiate breaking the block -- */
                aux = -gnsq[i] / gplus[i];
            }
        } else {
            /* check if a new block should be created */
            if (!USEBLOCKS || (inext == *j && *j != 1 && *j != *n)) {
                bres = FALSE_;
            } else if (inext == *j || omega[inext + dir] == 0) {
                /* laxer criterion, includes creating block at end */
                cplus = g[inext] + (ngn[i] - (*tau) * shfprt[inext]);
                bres = (fabs(gplus[i] * cplus) < PK1 * gnsq[i]);
            } else {
                doublereal max_val;
                doublereal a1 = fabs((*tau) * shfprt[inext]);
                doublereal a2 = fabs(g[inext]);
                max_val = a1 > a2 ? a1 : a2;
                bres = (gplus[i] == 0.0 || PMAXOSEQLEN > 0)
                    && (fabs(gplus[i]) * max_val < PK2 * gnsq[i])
                    && (fabs(gplus[i]) < PK2 * fabs(g[i]));
            }

            if (bres) {
                /* -- Initiate creation of a block -- */
                oplus[inext] = 1;
                aux = ngn[i];
            } else {
                /* -- Standard dstqds -- */
                if (gplus[i] == 0.0) gplus[i] = -(*pivmin);
                aux = ngn[i] * (smt / gplus[i]);
            }
        }
        oldsmt = smt;
        s[inext] = aux;

        iprev = i;
        i = inext;
        inext = inext + dir;

        /* =======================================================
           ==  Block structure change and overlap-control loop  ==
           ======================================================= */
        bcount = 0;
        if (oplus[i] != 0 && omega[inext] != 0) bcount = 1;

        for (;;) {
            if (i == *j || ((omega[i] == 0) == (oplus[i] == 0))) break;

            smt = aux - (*tau) * shfprt[i];
            gplus[i] = g[i] + smt;
            if (fabs(gplus[i]) < *pivmin) gplus[i] = 0.0;
            oplus[inext] = 0;

            /* reset overlap-counter if perturbations can be attributed to the shift */
            if (fabs(aux) < PROSQ * abstau) {
                bcount = 0;
            }

            if (oplus[i] != 0) {
                cplus = gplus[i];
                bdetp = gplus[iprev] * cplus - gnsq[iprev];
                if (omega[inext] != 0) {
                    /* = - ngnp(i) */
                    aux = -(gnsq[i] / bdetp) * gplus[iprev];
                } else {
                    /* -- End by create or clean create -- */
                    if (gplus[iprev] == 0.0) {
                        aux = ngn[i];
                    } else {
                        doublereal sgn_gprev = gplus[iprev] >= 0.0 ? 1.0 : -1.0;
                        doublereal sgn_g = g[iprev] >= 0.0 ? 1.0 : -1.0;
                        if (omega[iprev] != 0 ||
                            fabs(aux) <= PRCRE * abstau ||
                            sgn_gprev != sgn_g)
                        {
                            x = gplus[iprev] * smt - gnsq[iprev];
                        } else {
                            x = aux * oldsmt - gplus[iprev] * (*tau) * shfprt[i];
                        }
                        aux = ngn[i] * (x / bdetp);
                    }
                }
            } else {
                /* check if a new block (with overlap) may be created */
                if (!USEBLOCKS || (inext == *j && *j != 1 && *j != *n)) {
                    bres = FALSE_;
                } else {
                    bres = (inext == *j || bcount < PMAXOSEQLEN ||
                            fabs(gplus[i]) * gnsq[inext]
                            < (1.0 - KBLMOD) * (PROSQ * abstau) * gnsq[i]);
                    {
                        doublereal a1 = fabs((*tau) * shfprt[inext]);
                        doublereal a2 = fabs(g[inext]);
                        doublereal max_val = a1 > a2 ? a1 : a2;
                        bres = bres &&
                            (fabs(gplus[i]) * max_val < PK2 * gnsq[i]) &&
                            (fabs(gplus[i] * g[iprev]) < PK2 * fabs(bdet[iprev]));
                    }
                }
                if (bres) {
                    /* -- Create block, continuing the sequence -- */
                    oplus[inext] = 1;
                    aux = ngn[i];
                    bcount = bcount + 1;
                } else {
                    /* -- End by break or clean break -- */
                    if (gplus[i] == 0.0) gplus[i] = -(*pivmin);

                    if (g[iprev] == 0.0) {
                        aux = -gnsq[i] / gplus[i];
                    } else {
                        doublereal sgn_gprev = gplus[iprev] >= 0.0 ? 1.0 : -1.0;
                        doublereal sgn_g = g[iprev] >= 0.0 ? 1.0 : -1.0;
                        if (oplus[iprev] != 0 ||
                            fabs(aux) <= PRBRK * fabs((*tau) * shfprt[i]) ||
                            sgn_gprev != sgn_g)
                        {
                            x = g[iprev] * smt + gnsq[iprev];
                        } else {
                            x = -aux * oldsmt - g[iprev] * (*tau) * shfprt[i];
                        }
                        aux = (gnsq[i] * x) / (bdet[iprev] * gplus[i]);
                    }
                }
            }
            s[inext] = aux;
            oldsmt = smt;
            iprev = i;
            i = inext;
            inext = inext + dir;
        }
    }

    return 0;
}


/* ======================================================================
   DLAXRS_PROG: Progressive shift selection factorization

   Computes a progressive factorization from J1 to J2 with shift
   proportionality (SHFPRT). Non-blocked version.
   ====================================================================== */
int dlaxrs_prog__(integer *n, integer *j1, integer *j2,
    integer *i0, integer *i1, doublereal *g, integer *omega,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, doublereal *shfprt,
    doublereal *off, doublereal *gplus, integer *oplus,
    doublereal *p)
{
    doublereal aux;
    integer i, inext, j, dir;

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet;
    --shfprt;
    /* gplus is I0:I1 indexed — adjust by I0 not 1 */
    gplus -= *i0;
    --oplus;
    --p;

    i = *j1;
    j = *j2;

    if (*j1 < *j2 && *i0 == 1 && *i1 == *n - 1) {
        dir = 1;
    } else if (*j2 < *j1 && *i0 == 2 && *i1 == *n) {
        dir = -1;
    } else {
        i = -1;
        dir = 0;
    }
    inext = i + dir;

    aux = *off + (g[i] - (*tau) * shfprt[i]);
    for (;;) {
        gplus[i] = ngn[inext] + aux;
        if (fabs(gplus[i]) < *pivmin) gplus[i] = -(*pivmin);
        p[inext] = (g[inext] / gplus[i]) * aux;
        oplus[inext] = 0;
        if (inext == j) break;
        i = inext;
        inext = inext + dir;
        aux = p[i] - (*tau) * shfprt[i];
    }

    return 0;
}


/* ======================================================================
   DLAXRM_STAT2: Multi-negcount (2 shifts simultaneously)
   ====================================================================== */
int dlaxrm_stat2__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    const integer PRBRK = 8;

    doublereal sgngip, negpiv;
    integer i, iprev;

    doublereal aux1, gpl1, smt1, oldsmt1, sbound1;
    doublereal aux2, gpl2, smt2, oldsmt2, sbound2;
    integer negc1, negc2;
    logical brnchi1, brnchi2;

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet;
    --tau; --anegc; --aaux;
    --lbbk; /* Fortran 1-based LBBEGK(*) */

    negpiv = -(*pivmin);
    if (*dir == 1) {
        i = 1;
    } else {
        i = *n;
    }
    sbound1 = (PRBRK < *n ? PRBRK : *n) * fabs(tau[1]);
    sbound2 = (PRBRK < *n ? PRBRK : *n) * fabs(tau[2]);
    negc1 = 0; negc2 = 0;
    aux1 = 0.0; aux2 = 0.0;

    for (;;) {
        for (;;) {
            if (i == lbbk[*ibb]) break;
            aux1 = aux1 - tau[1];
            aux2 = aux2 - tau[2];
            gpl1 = g[i] + aux1;
            gpl2 = g[i] + aux2;
            if (fabs(gpl1) < *pivmin) gpl1 = negpiv;
            if (fabs(gpl2) < *pivmin) gpl2 = negpiv;
            if (gpl1 < 0.0) negc1++;
            if (gpl2 < 0.0) negc2++;
            aux1 = ngn[i] * (aux1 / gpl1);
            aux2 = ngn[i] * (aux2 / gpl2);
            i = i + *dir;
        }
        if (i == *k) break;

        *ibb = *ibb + *dir;

        smt1 = aux1 - tau[1];
        smt2 = aux2 - tau[2];
        gpl1 = g[i] + smt1;
        gpl2 = g[i] + smt2;
        if (fabs(gpl1) < *pivmin) gpl1 = negpiv;
        if (fabs(gpl2) < *pivmin) gpl2 = negpiv;
        if (gpl1 < 0.0) negc1++;
        if (gpl2 < 0.0) negc2++;
        aux1 = -gnsq[i] / gpl1;
        aux2 = -gnsq[i] / gpl2;
        iprev = i;
        i = i + *dir;
        if (i == *k) break;

        sgngip = g[iprev] >= 0.0 ? 1.0 : -1.0;
        oldsmt1 = smt1;
        brnchi1 = (fabs(aux1) <= sbound1 || sgngip != (gpl1 >= 0.0 ? 1.0 : -1.0));
        oldsmt2 = smt2;
        brnchi2 = (fabs(aux2) <= sbound2 || sgngip != (gpl2 >= 0.0 ? 1.0 : -1.0));
        smt1 = aux1 - tau[1];
        smt2 = aux2 - tau[2];
        gpl1 = g[i] + smt1;
        gpl2 = g[i] + smt2;
        if (fabs(gpl1) < *pivmin) gpl1 = negpiv;
        if (fabs(gpl2) < *pivmin) gpl2 = negpiv;
        if (gpl1 < 0.0) negc1++;
        if (gpl2 < 0.0) negc2++;
        if (g[iprev] == 0.0) {
            aux1 = -gnsq[i] / gpl1;
            aux2 = -gnsq[i] / gpl2;
        } else {
            if (brnchi1) {
                aux1 = g[iprev] * smt1 + gnsq[iprev];
            } else {
                aux1 = -aux1 * oldsmt1 - g[iprev] * tau[1];
            }
            if (brnchi2) {
                aux2 = g[iprev] * smt2 + gnsq[iprev];
            } else {
                aux2 = -aux2 * oldsmt2 - g[iprev] * tau[2];
            }
            aux1 = (gnsq[i] * aux1) / (bdet[iprev] * gpl1);
            aux2 = (gnsq[i] * aux2) / (bdet[iprev] * gpl2);
        }
        i = i + *dir;
    }

    anegc[1] = negc1;
    anegc[2] = negc2;
    aaux[1] = aux1;
    aaux[2] = aux2;
    return 0;
}


/* ======================================================================
   Generic dlaxrm_stat implementation for M simultaneous shifts
   Used to implement dlaxrm_stat4 through dlaxrm_stat64
   ====================================================================== */
static int dlaxrm_stat_generic_(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux, integer m)
{
    const integer PRBRK = 8;

    doublereal sgngip, negpiv;
    integer i, iprev, q;

    doublereal aux[64], gpl[64], smt[64], oldsmt[64], sbound[64];
    integer negc[64];
    logical brnchi[64];

    /* Parameter adjustments */
    --g; --gnsq; --ngn; --bdet;
    --tau; --anegc; --aaux;
    --lbbk; /* Fortran 1-based LBBEGK(*) */

    negpiv = -(*pivmin);
    if (*dir == 1) {
        i = 1;
    } else {
        i = *n;
    }
    for (q = 0; q < m; q++) {
        sbound[q] = (PRBRK < *n ? PRBRK : *n) * fabs(tau[q+1]);
        negc[q] = 0;
        aux[q] = 0.0;
    }

    for (;;) {
        for (;;) {
            if (i == lbbk[*ibb]) break;
            for (q = 0; q < m; q++) {
                aux[q] = aux[q] - tau[q+1];
                gpl[q] = g[i] + aux[q];
                if (fabs(gpl[q]) < *pivmin) gpl[q] = negpiv;
                if (gpl[q] < 0.0) negc[q]++;
                aux[q] = ngn[i] * (aux[q] / gpl[q]);
            }
            i = i + *dir;
        }
        if (i == *k) break;

        *ibb = *ibb + *dir;

        for (q = 0; q < m; q++) {
            smt[q] = aux[q] - tau[q+1];
            gpl[q] = g[i] + smt[q];
            if (fabs(gpl[q]) < *pivmin) gpl[q] = negpiv;
            if (gpl[q] < 0.0) negc[q]++;
            aux[q] = -gnsq[i] / gpl[q];
        }
        iprev = i;
        i = i + *dir;
        if (i == *k) break;

        sgngip = g[iprev] >= 0.0 ? 1.0 : -1.0;
        for (q = 0; q < m; q++) {
            oldsmt[q] = smt[q];
            brnchi[q] = (fabs(aux[q]) <= sbound[q] ||
                         sgngip != (gpl[q] >= 0.0 ? 1.0 : -1.0));
            smt[q] = aux[q] - tau[q+1];
            gpl[q] = g[i] + smt[q];
            if (fabs(gpl[q]) < *pivmin) gpl[q] = negpiv;
            if (gpl[q] < 0.0) negc[q]++;
        }
        if (g[iprev] == 0.0) {
            for (q = 0; q < m; q++) {
                aux[q] = -gnsq[i] / gpl[q];
            }
        } else {
            for (q = 0; q < m; q++) {
                if (brnchi[q]) {
                    aux[q] = g[iprev] * smt[q] + gnsq[iprev];
                } else {
                    aux[q] = -aux[q] * oldsmt[q] - g[iprev] * tau[q+1];
                }
                aux[q] = (gnsq[i] * aux[q]) / (bdet[iprev] * gpl[q]);
            }
        }
        i = i + *dir;
    }

    for (q = 0; q < m; q++) {
        anegc[q+1] = negc[q];
        aaux[q+1] = aux[q];
    }
    return 0;
}

int dlaxrm_stat4__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    return dlaxrm_stat_generic_(n, k, dir, g, ibb, lbbk, gnsq, ngn, bdet,
                                pivmin, tau, anegc, aaux, 4);
}

int dlaxrm_stat8__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    return dlaxrm_stat_generic_(n, k, dir, g, ibb, lbbk, gnsq, ngn, bdet,
                                pivmin, tau, anegc, aaux, 8);
}

int dlaxrm_stat16__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    return dlaxrm_stat_generic_(n, k, dir, g, ibb, lbbk, gnsq, ngn, bdet,
                                pivmin, tau, anegc, aaux, 16);
}

int dlaxrm_stat32__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    return dlaxrm_stat_generic_(n, k, dir, g, ibb, lbbk, gnsq, ngn, bdet,
                                pivmin, tau, anegc, aaux, 32);
}

int dlaxrm_stat64__(integer *n, integer *k, integer *dir,
    doublereal *g, integer *ibb, integer *lbbk,
    doublereal *gnsq, doublereal *ngn, doublereal *bdet,
    doublereal *pivmin, doublereal *tau, integer *anegc,
    doublereal *aaux)
{
    return dlaxrm_stat_generic_(n, k, dir, g, ibb, lbbk, gnsq, ngn, bdet,
                                pivmin, tau, anegc, aaux, 64);
}
