/* dlaxrf_selshf_f77.f -- translated by f2c (version 20240504).
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

/* Table of constant values */

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int dlaxrf_selshf__(integer *n, doublereal *repr, integer *
	repi, integer *icbeg, integer *icend, doublereal *lgap, doublereal *
	ugap, integer *ewl_ae__, doublereal *ewl_lu__, integer *rginfo, 
	doublereal *taubar, doublereal *gaptol, integer *maxcpo, integer *
	maxcpi, integer *maxnc, integer *maxnb, integer *ncand, integer *
	acloc, doublereal *actau, integer *nbatch, integer *abend, integer *
	iwork)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer i_sign(integer *, integer *), pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    logical goinside;
    integer i__, j, m, i0, j0;
    logical gooutside;
    doublereal lb, ub;
    integer kl, ku;
    doublereal off;
    integer dir;
    logical iok, jok;
    doublereal eps, off0;
    integer flag__, bloc;
    doublereal gapl;
    integer iloc;
    doublereal prec, gapu;
    integer ktry;
    doublereal raclb;
    integer bbtch;
    doublereal delta;
    integer iclen;
    doublereal avgap, ingap;
    integer ilbin, nlbin;
    doublereal bound;
    integer nscnd, index, ixf77a, nlocs;
    doublereal width;
    integer bslot, npart, nprim, itmp77, islot, nlout, ntrny;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibatch;
    doublereal mingap;
    integer ningap;
    doublereal offmax;
    integer prefer, nlocph;
    doublereal srnglb, ewsign;
    integer ilbout;
    doublereal srngub;
    integer nlbout, trymax, maxnloc;

/*     IMPLICIT NONE */




/*  Purpose */
/*  ======= */

/*     Produces shift candidates to try around a cluster. */

/*     For candidate 1 <= i <= NCAND, ACTAU(i) is the shift (absolute), */
/*     and the location is encoded in ACLOC(i): */
/*        ABS( ACLOC(i) )  is the eigenvalue we are shifting to, */
/*        SIGN( ACLOC(i) ) is the direction of the offset/gap, that is, */
/*                         +1 if shifting to the right, -1 if left. */
/*     As example, ACLOC(1)=-5 means the first shift goes to the left */
/*     of eigenvalue 5. This implies that all shifts to this location */
/*     are <= the lower bound of ew 5. */
/*       A special kind of candidates are those which are placed blindly */
/*     and not inside a gap, thos do not have a location but instead */
/*     ACLOC(i) is set to ICEND+1. */

/*     The shifts are grouped into up to NBATCH <= 2*MAXCPL batches */
/*                   1, ... , ABEND(1), */
/*          ABEND(1)+1, ... , ABEND(2) */
/*                      ... */
/*          ABEND(NBATCH-1), ... , ABEND(NBATCH), */
/*     where ABEND(NBATCH)=NCAND */

/*     Out: */
/*       NCAND=0,NBATCH=0 is possible. */
/*       Batches contain at least one candidate, and each candidate is */
/*       contained in exactly one batch (see above). */

/*  Arguments */
/*  ========= */

/*  MAXCPO  (input) INTEGER */
/*          Maximal number of candidates for both outside locations. */
/*          Set to 0 to deactivate outside shifts. */

/*  MAXCPI  (input) INTEGER */
/*          Maximal number of candidates for inside locations. */
/*          Set to 0 to deactivate inside shifts. */

/*  MAXNC   (input) INTEGER */
/*          Maximal number of candidates. */
/*          This is not a parameter, but must be set to */
/*            MAXNC  :=  2*MAXCPO + 2*(ICEND-ICBEG)*MAXCPI + 1. */
/*          It is checked at runtime that the correct value was set. */
/*          The purpose of this argument is just to make declaring the */
/*          array in the interface and at the caller easier. */

/*  MAXNB   (input) INTEGER */
/*          MAximal number of batches. */
/*          This is not a parameter, but must be set to */
/*            MAXNB  :=  3 * MAX(MAXCPO,MAXCPI) */
/*          It is checked at runtime that the correct value was set. */
/*          The purpose of this argument is just to make declaring the */
/*          array in the interface and at the caller easier. */


/*  ====================================================================== */

/*     .. Declarations .. */


/*     .. Constants .. */


/*     .. Parameters .. */


/*     .. Parameters .. */

/*         Cluster boundaries will be refined to relative accuracy */
/*         RACLBFAC*N*EPS */
/*     Maximal number of inside shifts for primary batches */
/*     For an inside shift to become primary, we require that the */
/*     smaller group of eigenvalues left/right has at least */
/*         MAX( 2, MININPRIMPART * (#ews in cluster) ) */
/*     elements. */
/*     Activate the last resort option to just throw a shift somewhere in the */
/*     middle of a tight cluster. */

/*     .. Locals .. */


/*  === Executable Statements ============================================ */

    /* Parameter adjustments */
    --repi;
    --repr;
    --iwork;
    --ewl_ae__;
    --ewl_lu__;
    --rginfo;
    --actau;
    --acloc;
    --abend;

    /* Function Body */
    iclen = *icend - *icbeg + 1;
    maxnloc = iclen << 1;
    eps = dlamch_("Epsilon", (ftnlen)7);
    prec = dlamch_("Precision", (ftnlen)9);
    lb = ewl_lu__[(*icbeg << 1) - 1];
    ub = ewl_lu__[*icend * 2];
    ewsign = d_sign(&c_b4, &lb);
    raclb = (*n << 2) * eps;
    goinside = *maxcpi > 0;
    gooutside = *maxcpo > 0;
    i0 = -1;
    j0 = -1;
/* Computing MIN */
    d__1 = *lgap * .25, d__2 = *gaptol * abs(lb);
    srnglb = lb - min(d__1,d__2);
/* Computing MIN */
    d__1 = *ugap * .25, d__2 = *gaptol * abs(ub);
    srngub = ub + min(d__1,d__2);
    if (*taubar > 0. && lb > 0.) {
/* Computing MIN */
	d__1 = srngub, d__2 = *taubar + lb * 2;
	srngub = min(d__1,d__2);
    } else if (*taubar < 0. && ub < 0.) {
/* Computing MAX */
	d__1 = srnglb, d__2 = *taubar + ub * 2;
	srnglb = max(d__1,d__2);
    }

/*     =================================================================== */
/*                         Select Shift Locations */
/*     =================================================================== */

/*     More precisely, we order them with the most promising first. */
/*     The locations are stored in IWORK. */
/*     Set PREFER TO +1 to emphasize shifts in the right half, */
/*     and to -1 for left. */
/*     TODO: prefer side "closer to min/max ew" (Parlett). Can */
/*           detect using taubar and spdiam */
/*     At the moment: prefer near side */
    if (lb > 0.) {
	prefer = -1;
    } else {
	prefer = 1;
    }
    ningap = 0;
    if (goinside) {
	m = (*icbeg + *icend) / 2;
	i0 = ewl_ae__[(m + (*icbeg + *icend) % 2 << 1) - 1];
	j0 = ewl_ae__[m * 2];
/*        For gaps <= I0 we normally shift only to left ends of subgroups */
/*        except as a last resort, similar for J0. */
/*        Loop once over all inner locations and determine where we want */
/*        them: */
/*          1 - eligible even for primary */
/*          2 - eligible for secondary */
/*          3 - only last resort */
/*        We store these flags temporarily in ACLOC(1:2*ILEN), where */
/*          ACLOC(2*(I-ICBEG) + 1)  is for shift left of I */
/*          ACLOC(2*(I-ICBEG) + 2)  is for shift right of I. */
/*        These flags are later set to zero as soon as the location was */
/*        chosen. The outer ones are always used, we set their flags to */
/*        0 as stoppers. */
	acloc[1] = 0;
	i__ = *icbeg;
L90001:
L90003:
	j = ewl_ae__[i__ * 2];
	i__ = j + 1;
	if (rginfo[j] > 0) {
	    goto L90004;
	}
	goto L90003;
L90004:
	if (j == *icend) {
	    goto L90002;
	}
/*           Inner gap ]J,I[ */
	++ningap;
/* Computing MIN */
	i__1 = i__ - *icbeg, i__2 = *icend - j;
	npart = min(i__1,i__2);
	if (npart < 2) {
	    acloc[(i__ - *icbeg << 1) + 1] = 3;
	    acloc[(j - *icbeg << 1) + 2] = 3;
	} else {
	    if ((doublereal) npart >= iclen * .4) {
		acloc[(i__ - *icbeg << 1) + 1] = 1;
		acloc[(j - *icbeg << 1) + 2] = 1;
	    } else {
		acloc[(i__ - *icbeg << 1) + 1] = 2;
		acloc[(j - *icbeg << 1) + 2] = 2;
	    }
	    if (i__ > i0) {
		acloc[(i__ - *icbeg << 1) + 1] = 3;
	    }
	    if (j < j0) {
		acloc[(j - *icbeg << 1) + 2] = 3;
	    }
	}
	goto L90001;
L90002:
	acloc[iclen * 2] = 0;
    }
    nlocs = 0;
/*     First pass: Inside shifts for primary batches. */
    if (goinside && ningap > 0) {
	i__ = i0;
	j = j0;
L90005:
	if (i__ == *icbeg && j == *icend || nlocs == 2) {
	    goto L90006;
	}
L90007:
	if (rginfo[j] > 0) {
	    goto L90008;
	}
	j = ewl_ae__[(j + 1) * 2];
	goto L90007;
L90008:
L90009:
	if (rginfo[i__ - 1] > 0) {
	    goto L90010;
	}
	i__ = ewl_ae__[(i__ - 1 << 1) - 1];
	goto L90009;
L90010:
/*           Now there are gaps left of I and right of J */
	iok = acloc[(i__ - *icbeg << 1) + 1] == 1;
	jok = acloc[(j - *icbeg << 1) + 2] == 1;
	if (iok && (! jok || prefer != 1)) {
	    ++nlocs;
	    iwork[nlocs] = -i__;
	    acloc[(i__ - *icbeg << 1) + 1] = 0;
	    --i__;
	} else if (jok) {
	    ++nlocs;
	    iwork[nlocs] = j;
	    acloc[(j - *icbeg << 1) + 2] = 0;
	    ++j;
	}
	if (! iok && i__ > *icbeg) {
	    --i__;
	}
	if (! jok && j < *icend) {
	    ++j;
	}
	goto L90005;
L90006:
	;
    }
/*     Outside locations */
    nlout = 0;
    if (gooutside) {
	if (prefer == -1) {
	    iwork[nlocs + 1] = -(*icbeg);
	    iwork[nlocs + 2] = *icend;
	} else {
	    iwork[nlocs + 1] = *icend;
	    iwork[nlocs + 2] = -(*icbeg);
	}
	nlocs += 2;
	nlout += 2;
    }
    nprim = nlocs;
/*     Secondary & Ternary: Remaining inside shifts, first those on the */
/*     right gap sides, then all others */
    nscnd = 0;
    ntrny = 0;
    if (goinside && ningap > 0) {
	for (flag__ = 2; flag__ <= 3; ++flag__) {
	    i__ = i0;
	    j = j0;
L90013:
	    if (i__ == *icbeg && j == *icend) {
		goto L90014;
	    }
L90015:
	    if (rginfo[j] > 0) {
		goto L90016;
	    }
	    j = ewl_ae__[(j + 1) * 2];
	    goto L90015;
L90016:
L90017:
	    if (rginfo[i__ - 1] > 0) {
		goto L90018;
	    }
	    i__ = ewl_ae__[(i__ - 1 << 1) - 1];
	    goto L90017;
L90018:
/*              Now there are gaps left of I and right of J */
	    iok = acloc[(i__ - *icbeg << 1) + 1] == flag__;
	    jok = acloc[(j - *icbeg << 1) + 2] == flag__;
	    if (iok && (! jok || prefer != 1)) {
		++nlocs;
		iwork[nlocs] = -i__;
		acloc[(i__ - *icbeg << 1) + 1] = 0;
		--i__;
	    } else if (jok) {
		++nlocs;
		iwork[nlocs] = j;
		acloc[(j - *icbeg << 1) + 2] = 0;
		++j;
	    }
	    if (! iok && i__ > *icbeg) {
		--i__;
	    }
	    if (! jok && j < *icend) {
		++j;
	    }
	    goto L90013;
L90014:
	    if (flag__ == 2) {
		nscnd = nlocs - nprim;
	    } else {
		ntrny = nlocs - nscnd - nprim;
	    }
/* L90011: */
	}
/* L90012: */
    }

/*     =================================================================== */
/*                         Place the Shift Candidates */
/*     =================================================================== */

    i__1 = *maxnb;
    for (itmp77 = 1; itmp77 <= i__1; ++itmp77) {
	abend[itmp77] = 0;
/* L90019: */
    }
/* L90020: */
    i__1 = *maxnc;
    for (itmp77 = 1; itmp77 <= i__1; ++itmp77) {
	acloc[itmp77] = 0;
/* L90021: */
    }
/* L90022: */
/*     We fill the candidate data fields ACLOC and ACTAU by scanning the */
/*     locations once from left to right. This leaves some entries unset */
/*     (ACLOC(i)=0). The fields will be compressed later. */
/*     At the same time we count in the entries of ABEND the number of */
/*     elements per batch. */
/*     The first slot of candidates within the current phase (primary, */
/*     secondary, ternary). */
    bslot = 1;
/*     Number of the first batch within the current phase */
    bbtch = 1;
/*     Number of locations in the current phase */
    nlocph = nprim;
/*     First location in current phase */
    bloc = 1;
/*     The current phase has NLBOUT outside and NLBIN inside locations. */
    nlbout = nlout;
    nlbin = nprim - nlbout;
    ilbout = 0;
    ilbin = 0;
/*      write(*,*) '================================================' */
/*      write(*,*) 'primary: nlocph = ',NLOCPH */
    *ncand = 0;
    i__1 = nlocs;
    for (iloc = 1; iloc <= i__1; ++iloc) {
	if (iloc == nprim + 1) {
/*           We now move to secondary locations. */
	    bslot = bslot + nlout * *maxcpo + (nprim - nlout) * *maxcpi;
	    bbtch += max(*maxcpo,*maxcpi);
	    nlocph = nscnd;
	    bloc = iloc;
/*            write(*,*) */
/*     $         'secondary: bslot=',bslot,'bbtch=',bbtch, */
/*     $         'nlocph=',nlocph,'bloc=',bloc */
	    ilbout = 0;
	    ilbin = 0;
	    nlbout = 0;
	    nlbin = nscnd;
	}
/*        Note that we might have no secondary or ternary locations */
/*        at all. */
	if (iloc == nprim + nscnd + 1) {
/*           We now move to ternary locations. */
	    bslot += nscnd * *maxcpi;
	    bbtch += *maxcpi;
	    nlocph = ntrny;
	    bloc = iloc;
/*            write(*,*) */
/*     $         'ternary: bslot=',bslot,'bbtch=',bbtch, */
/*     $         'nlocph=',nlocph,'bloc=',bloc */
	    ilbout = 0;
	    ilbin = 0;
	    nlbout = 0;
	    nlbin = ntrny;
	}
	index = (i__2 = iwork[iloc], abs(i__2));
	dir = i_sign(&c__1, &iwork[iloc]);
/*        We might be interested in just a chunk of the cluster */
	if (dir == 1) {
	    kl = *icbeg;
	    ku = index;
	    bound = ewl_lu__[index * 2];
	} else {
	    kl = index;
	    ku = *icend;
	    bound = ewl_lu__[(index << 1) - 1];
	}
/*        Determine minimal absolute gap separating the chunk */
	if (kl == *icbeg) {
	    gapl = *lgap;
	} else {
	    gapl = ewl_lu__[(kl << 1) - 1] - ewl_lu__[(kl << 1) - 2];
	}
	if (ku == *icend) {
	    gapu = *ugap;
	} else {
	    gapu = ewl_lu__[(kl << 1) + 1] - ewl_lu__[kl * 2];
	}
	mingap = min(gapl,gapu);
/*        See if we have an internal gap inside the chunk, on the */
/*        other side of the outermost interval */
	ingap = 0.;
	if (dir == 1) {
	    i__ = ewl_ae__[(index << 1) - 1];
	    if (i__ > kl) {
		ingap = ewl_lu__[(i__ << 1) - 1] - ewl_lu__[(i__ << 1) - 2];
	    }
	} else {
	    j = ewl_ae__[index * 2];
	    if (j < ku) {
		ingap = ewl_lu__[(j << 1) + 1] - ewl_lu__[j * 2];
	    }
	}
	width = ewl_lu__[ku * 2] - ewl_lu__[(kl << 1) - 1];
	avgap = 0.;
	if (kl != ku) {
	    avgap = width / (ku - kl);
	}
/*        Set maximal allowed offset */
	if (dir == 1) {
	    offmax = srngub - bound;
	} else {
	    offmax = bound - srnglb;
	}
/* Computing MIN */
	d__1 = offmax, d__2 = mingap * .25, d__1 = min(d__1,d__2), d__2 = abs(
		bound) * .25;
	offmax = min(d__1,d__2);
/*        now offmax is negative if bound is outside SRNG */
/*        Set maximal number of tries at this location, depending on */
/*        if we are inside or not. */
	if (iwork[iloc] == -(*icbeg) || iwork[iloc] == *icend) {
	    ++ilbout;
	    trymax = *maxcpo;
	} else {
	    ++ilbin;
	    trymax = *maxcpi;
	}
/*        The current location is the ILBOUT'th outer location and the */
/*        ILBIN'th inner location within the current phase. */
/*        Place the shifts at this location */
	ktry = 0;
	off0 = abs(bound) * (eps * 8);
/* Computing MAX */
	d__1 = off0, d__2 = max(avgap,ingap) / pow_ii(&c__2, &trymax);
	delta = max(d__1,d__2);
	if (off0 > offmax) {
	    off0 = 0.;
	}
L90025:
	off = off0 + (pow_ii(&c__2, &ktry) - 1) * delta;
	if (ktry >= trymax || off > offmax) {
	    goto L90026;
	}
	ibatch = bbtch + ktry;
/*            ISLOT  = BSLOT + KTRY*NLOCPH + ILOC-BLOC */
	islot = bslot - 1 + min(ktry,*maxcpo) * nlbout + min(ktry,*maxcpi) * 
		nlbin;
	if (ktry + 1 <= *maxcpo) {
	    islot += ilbout;
	}
	if (ktry + 1 <= *maxcpi) {
	    islot += ilbin;
	}
	actau[islot] = bound + dir * off;
	acloc[islot] = iwork[iloc];
/*            write(*,*) 'batch',IBATCH,'slot',ISLOT,'loc',IWORK(ILOC) */
/*     $        'bound',bound,'dir',dir,'off',off,'tau',actau(islot) */
	++ktry;
	++(*ncand);
	++abend[ibatch];
	goto L90025;
L90026:
/* L90023: */
	;
    }
/* L90024: */

/*     .. Decode Batch Boundaries .. */

    *nbatch = 0;
    j = 0;
    ibatch = 1;
    i__1 = *maxnb;
    for (ibatch = 1; ibatch <= i__1; ++ibatch) {
	m = abend[ibatch];
/*        Now M is the number of candidates in batch IBATCH */
	if (m > 0) {
	    ++(*nbatch);
	    j += m;
	    abend[*nbatch] = j;
	}
/* L90027: */
    }
/* L90028: */
    i__1 = *maxnb;
    for (ixf77a = *nbatch + 1; ixf77a <= i__1; ++ixf77a) {
	abend[ixf77a] = 0;
/* L90029: */
    }
/* L90030: */

/*     Special Backup: If there are no inside gaps (tight cluster) */
/*     we just throw a last resort shift somewhere inside in the */
/*     ternary batch. */
/*     These shifts could be regarded as outer ones, since they do */
/*     only depend on the outer bounds of the eigenvalues. However, */
/*     we think they should only be used as last resort, so they are */
/*     only placed if no inside shifts would be possible and in the last batch */
/* Computing MAX */
    d__1 = abs(lb), d__2 = abs(ub);
    if (ningap == 0 && *ncand < *maxnc && ub - lb <= sqrt(prec) * max(d__1,
	    d__2)) {
	++(*ncand);
/*        mark as special */
	acloc[*ncand] = *icend + 1;
/*        just take the midpoint */
	actau[*ncand] = (lb + ub) * .5;
	++abend[*nbatch];
    }

/*     .. Compress the candidate fields .. */

    i__ = 1;
L90031:
    if (i__ > *ncand) {
	goto L90032;
    }
    if (acloc[i__] == 0) {
	j = i__ + 1;
L90033:
	if (acloc[j] != 0) {
	    goto L90034;
	}
	++j;
	goto L90033;
L90034:
	acloc[i__] = acloc[j];
	actau[i__] = actau[j];
	acloc[j] = 0;
    }
    ++i__;
    goto L90031;
L90032:
    return 0;
} /* dlaxrf_selshf__ */

