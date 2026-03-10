      SUBROUTINE WSREQ_XRF(N, ICBEG, ICEND, REQR, REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)   ::  N, ICBEG, ICEND
      INTEGER, INTENT(OUT)  ::  REQR, REQI
      EXTERNAL DLAXRF
      DOUBLE PRECISION  R, A(2)
      INTEGER           I, T(2)
      LOGICAL           L
      REQR = -1
      REQI = -1
      CALL DLAXRF( N, A, R, R, I, R, ICBEG, ICEND,
     $             A, T, R, R, T, A, T,
     $             L, A, A, A, 1,
     $             A, T, R, R, T, A, R,
     $             A, REQR, T, REQI, I )
      END SUBROUTINE WSREQ_XRF
*
************************************************************************
*
      SUBROUTINE DLAXRF(
     $            N, E, GAPTOL, SPDIAM, DEPTH, TAUBAR, ICBEG, ICEND,
     $            FREPR, FREPI, FLGAP, FUGAP, FEWL_AE, FEWL_LU, FRGINFO,
     $            DOGRT, FREPR_GRT, SHFMOD, ENV, MODE,
     $            SREPR, SREPI, SLGAP, SUGAP, SEWL_AE, SEWL_LU, TAU,
     $            RWORK, LRWORK, IWORK, LIWORK, INFO
     $           )
      IMPLICIT NONE
* note: we want to be able to call this routine both within a df or a
* bf tree traversion strategy. For bf, the EWI-List of the son would
* be copied over the father afterwards and cached quants discarded,
* for df however we can directly use them.
*
      INTEGER,          INTENT(IN)  ::  N, DEPTH, ICBEG, ICEND, MODE
      INTEGER,          INTENT(IN)  ::  FREPI(6+N+N/2)
      INTEGER,          INTENT(IN)  ::  FRGINFO(ICBEG-1:ICEND)
      DOUBLE PRECISION, INTENT(IN)  ::  GAPTOL, SPDIAM, TAUBAR
      DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), FREPR(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  FREPR_GRT(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFMOD(3*N), ENV(N)
      LOGICAL,          INTENT(IN)  ::  DOGRT
*
      INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      INTEGER,          INTENT(INOUT)  ::  FEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  FLGAP, FUGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  FEWL_LU(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK )
*
      INTEGER,          INTENT(OUT)  ::  INFO
      INTEGER,          INTENT(OUT)  ::  SREPI(6+N+N/2)
      INTEGER,          INTENT(OUT)  ::  SEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP, TAU
      DOUBLE PRECISION, INTENT(OUT)  ::  SREPR(4*N+3)
      DOUBLE PRECISION, INTENT(OUT)  ::  SEWL_LU(2*ICBEG-1:2*ICEND)
*#
*  Purpose
*  =======
*
*     Find a shift such that the computed child representation is an RRR
*     (if possible).
*
*     Modes
*      1 - try both outside and inside locations
*      2 - place only shifts on the outside, or more precisely, that
*          depend on the outer bounds
*      3 - place only true inside shifts, that is, exactly those
*          not done in mode 2.
*
*
*     Return status
*      INFO=0  we found a shift TAU and the sons ewl-list and repdata are
*              set
*
*     For INFO != 0, we did not find an acceptable rep
*      INFO=1  no candidate found in select shift
*      INFO=2  none of the candidates we tried was acceptable
*
*  Arguments
*  =========
*
*  SHFMOD   (input) DOUBLE PRECISION array, dimension(3*N)
*           Perturbations (1*alpha) used to modify the shift, with the
*           intent of making emergence of gaps in children more likely.
*           Three sets should be provided:
*             SHFMOD(1:N)     <= 1, for shifting left of a bound
*             SHFMOD(N+1:2N)   = 1, for no shift perturbation
*             SHFMOD(2N+1:3N) >= 1, for shifting right of a bound
*           These can all be one if perturbing the shift is deactivated.
*
*  DOGRT    (input) LOGICAL
*           Specifies if for children without a gap or at least wide
*           enough local spectrum, a gap retry is to be initiated, using
*           the same shift but a perturbed father representation.
*
*  FREPR_GRT  (input) DOUBLE PRECISION array, dimension( (4*N+3) )
*           Perturbed real data of the father rep, to be used for gap
*           retrues.
*           Only accessed if DOGRT = .TRUE.
*
*  LRWORK   (input/output) INTEGER
*  LIWORK   (input/output) INTEGER
*           The lengths of the workspace arrays RWORK and IWORK,
*           respectively. They support a workspace query, if one or both
*           of them is -1, the routine just computes the required
*           workspace and sets LRWORK and LIWORK to it.
*           A ws query needs N, ICBEG, ICEND and MODE to be set,
*           no other argument is referenced.
*           Recommended: Don't use this directly, call WSREQ_XRF.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRS(
     $             N, IA, IE, REPR, REPI,
     $             TAU, SHFPRT,
     $             DPLUS, OMGADP, RPLUS, OMGARP, GAMMAP, TWISTOK,
     $             RWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, IA, IE
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  TAU
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), SHFPRT(N)
*
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK(2*N)
*
      INTEGER,          INTENT(OUT)  ::  OMGADP(N), OMGARP(N)
      INTEGER,          INTENT(OUT)  ::  TWISTOK(IA:IE)
      DOUBLE PRECISION, INTENT(OUT)  ::  DPLUS(1:N-1), RPLUS(2:N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GAMMAP(IA:IE)
      END SUBROUTINE DLAXRS
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRF_SELSHF(
     $             N, REPR, REPI,
     $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU, RGINFO,
     $             TAUBAR, GAPTOL, MAXCPO, MAXCPI, MAXNC, MAXNB,
     $             NCAND, ACLOC, ACTAU, NBATCH, ABEND,
     $             IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND
      INTEGER,          INTENT(IN)  ::  MAXCPO, MAXCPI, MAXNC, MAXNB
      INTEGER,          INTENT(IN)  ::  RGINFO(ICBEG-1:ICEND)
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  TAUBAR, GAPTOL
*
      INTEGER,          INTENT(INOUT) :: EWL_AE(2*ICBEG-1:2*ICEND)
      INTEGER,          INTENT(INOUT) :: IWORK( 2*(ICEND-ICBEG+1) )
      DOUBLE PRECISION, INTENT(INOUT) :: LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT) :: EWL_LU(2*ICBEG-1:2*ICEND)
*
      INTEGER,          INTENT(OUT) :: NCAND, NBATCH
      INTEGER,          INTENT(OUT) :: ACLOC( MAXNC )
      INTEGER,          INTENT(OUT) :: ABEND( MAXNB )
      DOUBLE PRECISION, INTENT(OUT) :: ACTAU( MAXNC )
      END SUBROUTINE DLAXRF_SELSHF
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRF_SELTW(
     $             N, KMINUS, D, OMEGAD, R, OMEGAR, E, GAMMA, TWISTOK,
     $             ENV, SPDIAM, MINGAP,
     $             DIR, LOC, TAU, WINOK,
     $             K, EVAL, XIZERO,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      LOGICAL,          INTENT(IN)     ::  WINOK
      INTEGER,          INTENT(IN)     ::  N, KMINUS, DIR, LOC
      INTEGER,          INTENT(IN)     ::  OMEGAD(N), OMEGAR(N)
      INTEGER,          INTENT(IN)     ::  TWISTOK( N )
      DOUBLE PRECISION, INTENT(IN)     ::  SPDIAM, MINGAP, TAU
      DOUBLE PRECISION, INTENT(IN)     ::  D(1:N-1), R(2:N), E(N-1)
      DOUBLE PRECISION, INTENT(IN)     ::  ENV(N)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( N+2 )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 5*N )
      DOUBLE PRECISION, INTENT(INOUT)  ::  GAMMA(N)
*
      INTEGER,          INTENT(OUT)    ::  K, XIZERO
      DOUBLE PRECISION, INTENT(OUT)    ::  EVAL
      END SUBROUTINE DLAXRf_SELTW
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRF_COB(
     $             N, FLGAP, FLB, FUB, FUGAP, REPR, REPI, TAU, GAPTOL,
     $             ICBEG, ICEND, XIZERO, SLGAP, SUGAP, EWL_AE, EWL_LU,
     $             STATUS
     $         )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND, XIZERO
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  FLGAP, FLB, FUB, FUGAP
      DOUBLE PRECISION, INTENT(IN)  ::  TAU, GAPTOL
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
*
      INTEGER,          INTENT(OUT)  ::  STATUS
      INTEGER,          INTENT(OUT)  ::  EWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP
      DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*ICBEG-1:2*ICEND)
      END SUBROUTINE DLAXRF_COB
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRF_IIB(
     $             N, REPR, REPI, TAU, ICBEG, ICEND,
     $             FEWL_AE, FEWL_LU,
     $             SLGAP, SUGAP, SEWL_AE, SEWL_LU
     $         )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      INTEGER,          INTENT(IN)  ::  FEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(IN)  ::  FEWL_LU(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(IN)  ::  TAU
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  SEWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  SLGAP, SUGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  SEWL_LU(2*ICBEG-1:2*ICEND)
      END SUBROUTINE DLAXRF_IIB
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRR( N, K, TYPE, E, PIVBASE, REPR, REPI )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, TYPE
      DOUBLE PRECISION, INTENT(IN)  ::  PIVBASE
      DOUBLE PRECISION, INTENT(IN)  ::  E(1:N-1)
*
      INTEGER,          INTENT(INOUT)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(INOUT)  ::  REPR(4*N+3)
      END SUBROUTINE DLAXRR
      END INTERFACE
      INTERFACE
      FUNCTION DLAXRN(N, REPR, REPI, TAU)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN
      INTEGER,          INTENT(IN)  ::  N
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU
      END FUNCTION DLAXRN
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRL_REFINE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, I, J, LAMBDA, XI
     $           )
      IMPLICIT NONE
*
      INTEGER         , INTENT(IN)  ::  IL, IU, I, J, XI
      DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      END SUBROUTINE DLAXRL_REFINE
      END INTERFACE
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0

*
*     .. Parameters ..
*


*
*     .. Parameters ..
*
*     The evaluation determines a real value EVAL for each candidate.
*     Based on this we categorize them broadly as follows
*        EVAL <= 1        passes a priori, will take it directly
*        EVAL <= EVALOK   is ok
*        EVAL >  EVALOK   only take as last resort
      DOUBLE PRECISION, PARAMETER  ::  EVALOK = 5.0D0

*     To transform a father bound between bisection and shifted zero
*     inertia semantics, we allow up to MAXNFBTRF tries, where the
*     relative backoff at the K'th try is increased by
*       FBTRFGRAN**K * N*EPS,
*     thus
      LOGICAL, PARAMETER  ::  WANTINEROK = .TRUE.
      INTEGER, PARAMETER  ::  FBTRFGRAN = 4
      INTEGER, PARAMETER  ::  MAXNFBTRF = 4

*     For eigenvalue underflow with a consistent negcount, fudge the
*     shift by |tau|*FUDGEFAC*Eps away in the right direction and try
*     again.
*     Set either one to 0 to deactivate fudging.
      INTEGER, PARAMETER  ::  FUDGEFAC = 2
      INTEGER, PARAMETER  ::  MAXFUDGE = 1

*     If true, the eigenvalues of a candidate that has passed all other
*     tests are refined until we can see that the cluster is wide enough.
*     If that is not possible, then either a gap retry is initiated, if
*     that is activated, or the candidate is discarded.
      LOGICAL, PARAMETER  ::  FORCEGAP = .TRUE.

*     Separate max number of candidates for inside locations, effective
*     will be minimum of this and MAXCPL
      INTEGER, PARAMETER  ::  MAXCPL_OUTSIDE = 3
      INTEGER, PARAMETER  ::  MAXCPL_INSIDE  = 2

*
*     .. Constants ..
*
      INTEGER, PARAMETER  ::  MAXNBATCH =
     $                            3*MAX( MAXCPL_OUTSIDE, MAXCPL_INSIDE )
*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, PREC, PIVBASE
      DOUBLE PRECISION  CLB, CUB, RCLWID, EVAL, BBESTEVAL
      DOUBLE PRECISION  WNTRW
      DOUBLE PRECISION  XA, XB, XC, XD, MID

      INTEGER           JXG, JYOMGA, WSREQR, WSREQI, KMINUS
      INTEGER           IXDP, IXRP, IXGMAP, IXBLUP, IXWORK
      INTEGER           IYTWOK, IYOMDP, IYOMRP, IYWORK
      INTEGER           IYCLOC, IYCTWI, IYCXIZ, IYCFDG, IYCGRT
      INTEGER           IXCTAU, IXCEVL
      INTEGER           NCAND, NBATCH, ICLEN
      INTEGER           SHPSEL, SHPOFF
      INTEGER           ICAND, IBATCH, BBEG, BEND, INDFBB
      INTEGER           NIBWEVAL, NIBACTIVE
      INTEGER           I, J, KL, KU, DIR, LOC, TWIST, XIZERO, XI
      INTEGER           NFBTRF, NFUDGE, COBSTA, IIBSTA
      INTEGER           ABEND( MAXNBATCH )
      INTEGER           MAXCPO, MAXCPI, MAXNCAND, MAXNB

      LOGICAL           LASTHOPE, HADEVAL, DOKILL, INEROK, SUCCESS
      LOGICAL           LIBNONEVAL, CONT_A, CONT_B, CONT_C, ISGAPRETRY

      INTEGER BISTRYCNT
      COMMON /XMRSTATS/
     $       XTIME1, XTIME2, XTIME3,
     $       XDDDDD,
     $       XNBLCKS,
     $       XNNODES, XMAXDEPTH, XNUMFN, XNUMFT, XNUMGV, XNUMGV0,
     $       XNUMFS_2, XNUMFS_K,
     $       XNBIS_INIT,
     $       XNBIS_COB, XNBIS_IIB, XNBIS_CLASS, XNBIS_SNG, XNBIS_CLB,
     $       XNBIS_WASTED,
     $       XNRQI, XNRQIBIS, XMAXNRQI, XMAXNRQIBIS,
     $       XNENVGV, XNENVTF,
     $       XIIIII,
     $       XSTEALTHMODE
*
*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*           !!!! SYNCHRONIZE ANY CHANGES HERE WITH xmr.h !!!!!
*           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*     Accumulative part: All these are accumulated over multiple
*     calls of stexr/laxrv.

*        Times for stages 1-3 (before/during/after block loop)
         DOUBLE PRECISION  XTIME1, XTIME2, XTIME3

*        Marker to verify data integrity
         DOUBLE PRECISION  XDDDDD

         INTEGER XNBLCKS
*
         INTEGER XNNODES, XMAXDEPTH

*        Number of calls to dlaxrn, dlaxrt, dlaxrg, and
*        dlaxrs (all twists/single twist)
         INTEGER XNUMFN, XNUMFT, XNUMGV, XNUMGV0, XNUMFS_2, XNUMFS_K

*        Number of bisection steps in various parts.
         INTEGER XNBIS_INIT
         INTEGER XNBIS_COB, XNBIS_IIB, XNBIS_CLASS, XNBIS_SNG, XNBIS_CLB

         INTEGER XNBIS_WASTED

*        Count/max rqi steps
         INTEGER XNRQI, XNRQIBIS
         INTEGER XMAXNRQI, XMAXNRQIBIS

*        For computing envelopes, number of full twisted factos
*        and gv-calls
         INTEGER XNENVGV, XNENVTF

*        Marker to verify data integrity
         INTEGER XIIIII
      INTEGER IXF77A, IXF77B, IXF77C

*    Temporaries and Configuration
*        If true, calls to internal routines dlaxr[n,t,s] are not
*        counted. Higher level routines like DLAXRB are not affected.
*        We need this to hide bisections steps done within aux
*        routines like xmr_estrc.
         LOGICAL XSTEALTHMODE

      SAVE /XMRSTATS/

*
*  ===== Executable Statements ==========================================
*
      ICLEN = ICEND - ICBEG + 1
      KMINUS = FREPI(2)

      MAXCPO = MAXCPL_OUTSIDE
      MAXCPI = MAXCPL_INSIDE
      IF( MODE .EQ. 2 )THEN
         MAXCPI = 0
      ELSEIF( MODE .EQ. 3 )THEN
         MAXCPO = 0
      ENDIF
C One idea may be to limit the backoffs for deeper levels, to prefer close
C ones to slightly better (but still bad) a priori criteria.
c Better solution:  Introduce relative offset from base bound into eval, maybe
C   like eval  <-  eval^2 * reloff
C      IF( DEPTH.GT.0 )THEN
C         MAXCPO = MIN( MAXCPO, 1 )
C         MAXCPI = MIN( MAXCPI, 1 )
C      ENDIF
CC~
      MAXNCAND = 2*MAXCPO + 2*(ICEND-ICBEG)*MAXCPI + 1
      MAXNB = 3 * MAX( MAXCPO, MAXCPI )
*
*     ----------------------
*      Workspace Allocation
*     ----------------------
*
      IXWORK = 1
      IYWORK = 1

*     Buffer arrays for DLAXRS
*       DP    real  dimension(N-1)
*       RP    real  dimension(N-1)
*       GMAP  real  dimension(N)
*       TWOK  int   dimension(N)
      IXDP   = IXWORK
      IXRP   = IXDP   + N-1
      IXGMAP = IXRP   + N-1
      IXWORK = IXWORK + 3*N-2

      IYTWOK = IYWORK
      IYOMDP = IYTWOK + N
      IYOMRP = IYOMDP + N
      IYWORK = IYWORK + 3*N

*     Data per candidate (NCAND <= 2*ILEN), indexed from 0
*     Set by SELSHF:
*       CLOC (int)   ew index we shift to, sign gives direction
*       COFF (real)  absolute offset from the ew bound
*     Set after evaluation:
*       CTWI (int)   optimal twist index
*       CXIZ (int)   zero inertia
*       CEVL (real)  evaluation
*     General:
*       CFDG (int)   counter for number of small fudges to avoid
*                    ew-underflow
*       CGRT (int)   set to one for a retry with the perturbed rep
*                    to get a wider cluster or even a gap
      IYCLOC = IYWORK
      IYCTWI = IYWORK +   MAXNCAND
      IYCXIZ = IYWORK + 2*MAXNCAND
      IYCFDG = IYWORK + 3*MAXNCAND
      IYCGRT = IYWORK + 4*MAXNCAND
      IYWORK = IYWORK + 5*MAXNCAND

      IXCTAU = IXWORK
      IXCEVL = IXCTAU + MAXNCAND
      IXWORK = IXCEVL + MAXNCAND

*     Blowup of father bounds to get consistent zero inertias.
*     Indexed (2*ICBEG-1:2*ICEND), conformingly with EWL.
      IXBLUP = IXWORK - (2*ICBEG-1)
      IXWORK = IXWORK + 2*ICLEN
*
*     .. Handle Workspace Query ..
*
*      Additional ws for called routines:
*       2N real  for DLAXRS
*       5N real N+2 int  for _SELTW
*       2N int   for _SELSHF  [can use IWORK]
*
      WSREQR = IXWORK-1 + 5*N
      WSREQI = IYWORK-1 + N+2
      IF( LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )THEN
         IF( LIWORK.EQ.-1 )  LIWORK = WSREQI
         IF( LRWORK.EQ.-1 )  LRWORK = WSREQR
         RETURN
      ENDIF
      IF( LIWORK.LT.WSREQI .OR. LRWORK.LT.WSREQR )THEN
         INFO = -1
         RETURN
      ENDIF

*     -------------------------------------------------------------

      EPS     = DLAMCH('Epsilon')
      PREC    = DLAMCH('Precision')
      JXG     = (0)
      JYOMGA  = (4)
      PIVBASE = FREPR(4*N+3)
      CLB     = FEWL_LU(2*ICBEG-1)
      CUB     = FEWL_LU(2*ICEND)
      RCLWID  = (CUB - CLB) / MAX(ABS(CLB),ABS(CUB))

*     Setup tolerance for how wide the eigenvalues of the child have
*     to be so that we accept it (only used if FORCEGAP = .TRUE.).
*     [ Note: don't make this too small, otherwise the refinement
*       for revealing a gap becomes too expensive ]
*     We take the candidate immediately if its ews have relative width
*     exceeding WNTRW1, if it only exceeds WNTRW0 we may keep the
*     candidate as fallback.
      WNTRW = MAX( GAPTOL**2, MIN( GAPTOL, RCLWID / GAPTOL ) )

*
*     .. Select which random perturbations for the shift to use ..
*
*      Set SHPOFF=0 to deactivate, SHPOFF=1 otherwise
      SHPOFF = 0
      IF( (CUB-CLB) .LE. SQRT(PREC)*MAX(ABS(CLB),ABS(CUB)) )THEN
         SHPOFF = 1
      ENDIF
*
*     .. Prepare the candidates ..
*

      CALL DLAXRF_SELSHF( N, FREPR, FREPI,
     $       ICBEG, ICEND, FLGAP, FUGAP, FEWL_AE, FEWL_LU, FRGINFO,
     $       TAUBAR, GAPTOL, MAXCPO, MAXCPI, MAXNCAND, MAXNB,
     $       NCAND, IWORK(IYCLOC), RWORK(IXCTAU), NBATCH, ABEND,
     $       IWORK
     $     )
*     Note: can use full IWORK, as not yet in use

      IF( NCAND .EQ. 0 )THEN
         INFO = 1
         RETURN
      ENDIF
*
*     ===================================================================
*     =                          MAIN LOOP                              =
*     ===================================================================
*
*     Candidates go through several phases:
*       active if LOC != 0
*       evaluated if TWI != 0
*     Init all candidates as actuve but not evaluated and unfudged.
      DO IXF77A = IYCTWI, IYCTWI+MAXNCAND-1
         IWORK(IXF77A) = 0
      ENDDO
      DO IXF77A = IYCFDG, IYCFDG+MAXNCAND-1
         IWORK(IXF77A) = 0
      ENDDO
      DO IXF77A = IYCGRT, IYCGRT+MAXNCAND-1
         IWORK(IXF77A) = 0
      ENDDO
*
*     The blown up father bounds are zero as long as they are not set.
      DO IXF77A = IXBLUP + 2*ICBEG-1, IXBLUP + 2*ICEND
         RWORK(IXF77A) = ZERO
      ENDDO
*
*     BBEG:BEND hold the offsets (0-based) of the candidates in the
*     current batch.
      BBEG = 0
      BEND = ABEND(1)-1
      IBATCH = 1
      LASTHOPE = .FALSE.
      SUCCESS  = .FALSE.

      candloop: DO

         IF( IBATCH .GT. NBATCH+1 )THEN
            EXIT candloop
         ENDIF

*        ----------------------------------------------------------------
*                          Select Next Candidate
*        ----------------------------------------------------------------
*        We evaluate the candidates in the order they are given.
*        If for the current batch all are evaluated we proceed with the
*        best one.
         BBESTEVAL = -ONE
         NIBWEVAL = 0
         NIBACTIVE = 0
         DO I = BBEG, BEND
            IF( IWORK(IYCLOC + I).NE.0 )THEN
*              candidate is active ..
               NIBACTIVE = NIBACTIVE + 1
               IF( IWORK(IYCTWI + I).NE.0 )THEN
*                 .. and has been evaluated
                  NIBWEVAL = NIBWEVAL + 1
                  IF( NIBWEVAL .EQ. 1 )THEN
                     BBESTEVAL = RWORK(IXCEVL + I)
                  ELSE
                     BBESTEVAL = MIN( BBESTEVAL, RWORK(IXCEVL + I) )
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         ICAND = -1
         DO I = BBEG, BEND
            IF( IWORK(IYCLOC + I).NE.0 )THEN
*              candidate is active ..
               IF( IWORK(IYCTWI + I).EQ.0 )THEN
*                 .. but has not yet been evaluated
                  ICAND = I
                  EXIT
               ENDIF
*              .. and has already been evaluated ..
               IF( RWORK(IXCEVL + I).LE.EVALOK .OR. LASTHOPE )THEN
*                 .. and applies to current phase
                  IF( ICAND.EQ.-1 .OR.
     $                RWORK(IXCEVL + I).LT.RWORK(IXCEVL + ICAND) )
     $            THEN
                     ICAND = I
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         IF( ICAND .EQ. -1 )THEN
*           no candidate found in current batch, proceed to the next
            IBATCH = IBATCH+1
            IF( IBATCH .LE. NBATCH )THEN
               BBEG = ABEND(IBATCH-1)
               BEND = ABEND(IBATCH)-1
            ELSEIF( IBATCH .EQ. NBATCH+1)THEN
               BBEG = 0
               BEND = NCAND-1
C Take fallback if we have one
               LASTHOPE = .TRUE.
            ELSE
               BBEG = -1
               BEND = -1
            ENDIF
            CYCLE candloop
         ENDIF

         TWIST = 0
         EVAL  = -ONE

*        **************************************************************
*        *  Now a series of tests are performed to see if this is an  *
*        *  acceptable candidate. To keep the loop structure simple,  *
*        *  as soon as a test fails, we modify the candidate record   *
*        *  (eg to deactivate or signal as last resort) and cycle     *
*        *  back to the beginning of the loop.                        *
*        **************************************************************
         LOC = ABS( IWORK(IYCLOC + ICAND) )
         IF( LOC .GT. ICEND )THEN
*           Special case: Shift not coupled to an eigenvalue
            DIR = 0
            LOC = 0
         ELSE
            DIR = SIGN( 1, IWORK(IYCLOC + ICAND) )
         ENDIF
         TAU    = RWORK(IXCTAU + ICAND)
         NFUDGE = IWORK(IYCFDG + ICAND)
         SHPSEL = N+1 + DIR*N*SHPOFF
         HADEVAL    = ( IWORK(IYCTWI + ICAND) .NE. 0 )
         LIBNONEVAL = ( .NOT.HADEVAL .AND. NIBWEVAL.EQ.NIBACTIVE-1 )
         ISGAPRETRY = ( IWORK(IYCGRT + ICAND) .GT. 0 )
         DOKILL     = .FALSE.

*        Special cases:
*         LIBNONEVAL indicates that all other active candidates in the
*         batch have an evaluation already, only this candidate did not.
*         BBESTEVAL is the best evaluation of a candidate that we looked
*         at this round; this equals EVAL if HADEVAL is true.



*        If this is the first shift to this location, we may allow
*        backing off to get the zero inertia consistent.
*        Otherwise, if previous inflating of the bound came to close
*        to this shift, we kill it.
         INDFBB = 0
         IF( DIR .NE. 0 )THEN
            DOKILL = .FALSE.
            IF( DIR.EQ.+1 )THEN
               INDFBB = IXBLUP + 2*LOC
               IF( RWORK(INDFBB).NE.ZERO .AND. .NOT.HADEVAL .AND.
     $             NFUDGE.EQ.0 .AND. .NOT.ISGAPRETRY )
     $         THEN
                  DOKILL = TAU.LE.RWORK(INDFBB)
               ENDIF
            ELSE
               INDFBB = IXBLUP + 2*LOC-1
               IF( RWORK(INDFBB).NE.ZERO .AND. .NOT.HADEVAL .AND.
     $             NFUDGE.EQ.0 .AND. .NOT.ISGAPRETRY )
     $         THEN
                  DOKILL = TAU.GE.RWORK(INDFBB)
               ENDIF
            ENDIF
            IF( DOKILL )THEN
*              too close, deactivate
               IWORK(IYCLOC + ICAND) = 0
               CYCLE candloop
            ENDIF
         ENDIF
*
*        ----------------------------------------------------------------
*                            Perform the Shift
*        ----------------------------------------------------------------
*
*        Retry until the inertia is consistent.
         NFBTRF = 0
         retryloop: DO

            IF( HADEVAL )THEN
               KL = IWORK(IYCTWI + ICAND)
               KU = KL
            ELSE
               KL = 1
               KU = N
            ENDIF

            IF( .NOT. ISGAPRETRY )THEN
               CALL DLAXRS(
     $             N, KL, KU, FREPR, FREPI, TAU, SHFMOD(SHPSEL),
     $             RWORK(IXDP), IWORK(IYOMDP),
     $             RWORK(IXRP), IWORK(IYOMRP),
     $             RWORK(IXGMAP-1+KL), IWORK(IYTWOK-1+KL), RWORK(IXWORK)
     $         )
            ELSE
*              Same call as above, just with FREPR_GRT instead of FREPR
               CALL DLAXRS(
     $             N, KL, KU, FREPR_GRT, FREPI, TAU, SHFMOD(SHPSEL),
     $             RWORK(IXDP), IWORK(IYOMDP),
     $             RWORK(IXRP), IWORK(IYOMRP),
     $             RWORK(IXGMAP-1+KL), IWORK(IYTWOK-1+KL), RWORK(IXWORK)
     $         )
            ENDIF

*           Select twist
            IF( HADEVAL )THEN
               TWIST  = IWORK(IYCTWI + ICAND)
               XIZERO = IWORK(IYCXIZ + ICAND)
               EVAL   = RWORK(IXCEVL + ICAND)
            ELSE

*              Need full data here, KL=1, KU=N
               CALL DLAXRF_SELTW( N, KMINUS,
     $                RWORK(IXDP), IWORK(IYOMDP),
     $                RWORK(IXRP), IWORK(IYOMRP),
     $                E, RWORK(IXGMAP), IWORK(IYTWOK),
     $                ENV, SPDIAM, MIN(FLGAP,FUGAP),
     $                DIR, LOC, TAU, WANTINEROK.AND.(DIR.NE.0),
     $                TWIST, EVAL, XIZERO,
     $                RWORK(IXWORK), IWORK(IYWORK)
     $              )

               IF( TWIST .EQ. -1 )THEN
                  IWORK(IYCLOC + ICAND) = 0
                  CYCLE candloop

               ENDIF
            ENDIF

*           Determine if the inertia is consistent
            IF( TWIST .EQ. 0 )THEN
               INEROK = .FALSE.
            ELSE
               INEROK = (  DIR.EQ.0  .OR.
     $                    (DIR.EQ.-1 .AND. XIZERO.LE.(2*LOC-1)) .OR.
     $                    (DIR.EQ.+1 .AND. XIZERO.GE.(2*LOC-1)) )

            ENDIF

            IF( .NOT.INEROK .AND. WANTINEROK )THEN
*              Initiate a retry if possible, or kill the candidate.
*              An fbtrf-retry is only considered if the transformed
*              father bound is not yet set. This implies that the
*              candidate has no evaluation yet.
               IF( RWORK(INDFBB).EQ.ZERO .AND. NFBTRF.LE.MAXNFBTRF )THEN
                  NFBTRF = NFBTRF+1
                  TAU = TAU + DIR*ABS(TAU)*(FBTRFGRAN**NFBTRF)*N*EPS
C                  TAU = TAU + DIR*ABS(TAU)*((2**NFBTRF)*FBTRFGRAN*N*EPS)
                  CYCLE retryloop
               ELSE
*                 No further retry allowed. Kill the candidate.
                  IWORK(IYCLOC + ICAND) = 0
                  CYCLE candloop
               ENDIF
            ENDIF
*
            EXIT
         ENDDO retryloop

*        Now that the inertia is ok, update blown up father bounds
         IF( DIR .NE. 0 )THEN
            RWORK(INDFBB) = TAU
         ENDIF

         IWORK(IYCTWI + ICAND) = TWIST
         IWORK(IYCXIZ + ICAND) = XIZERO
         RWORK(IXCEVL + ICAND) = EVAL
         RWORK(IXCTAU + ICAND) = TAU

*        Look at evaluation.
*        We only proceed with this candidate if
*          (a) it is a priori ok   (EVAL <= 1)
*        or
*          (b) it was the best in the batch (HADEVAL) and either
*              is acceptable (EVAL <= EVALOK) or we are desperate
*              (LASTHOPE)
*        or
*          (c) it was the only one without eval in the batch, but
*              after evaluating it we can see now that it is the
*              best
*              [special case to avoid unnessesary recomputation, if
*               we put it back now it is selected immediately again]
*
         CONT_A = (EVAL .LE. ONE)
         CONT_B = (HADEVAL .AND. (EVAL.LE.EVALOK .OR. LASTHOPE))
         CONT_C = (LIBNONEVAL .AND. EVAL.LE.EVALOK .AND.
     $              (BBESTEVAL.EQ.-ONE .OR.
     $               EVAL.LE.MIN(EVALOK,BBESTEVAL)))

         IF( .NOT. (CONT_A .OR. CONT_B .OR. CONT_C) )THEN
            CYCLE candloop
         ENDIF
*
*        ----------------------------------------------------------------
*                            Init cached quantities
*        ----------------------------------------------------------------
*
         DO I = 1, TWIST-1
            SREPR(JXG + I)    = RWORK(IXDP-1 + I)
            SREPI(JYOMGA + I) = IWORK(IYOMDP-1 + I)
         ENDDO
*
         SREPR(JXG + TWIST)    = RWORK(IXGMAP-1 + TWIST)
         IF( TWIST .EQ. 1 )THEN
            SREPI(JYOMGA + TWIST) = IWORK(IYOMRP-1 + TWIST)
         ELSEIF( TWIST .EQ. N) THEN
            SREPI(JYOMGA + TWIST) = IWORK(IYOMDP-1 + TWIST)
         ELSE
            SREPI(JYOMGA + TWIST) = 0
         ENDIF
*
         DO I = TWIST+1, N
            SREPR(JXG + I)    = RWORK(IXRP-2 + I)
            SREPI(JYOMGA + I) = IWORK(IYOMRP-1 + I)
         ENDDO

         CALL DLAXRR( N, TWIST, FREPI(1), E,
     $                PIVBASE, SREPR, SREPI )
*
*        ----------------------------------------------------------------
*                            Check Outer Bounds
*        ----------------------------------------------------------------
*
         BISTRYCNT = XNUMFN
         CALL DLAXRF_COB( N,
     $             FLGAP, FEWL_LU(2*ICBEG-1), FEWL_LU(2*ICEND), FUGAP,
     $             SREPR, SREPI, TAU, GAPTOL, ICBEG, ICEND, XIZERO,
     $             SLGAP, SUGAP, SEWL_AE, SEWL_LU, COBSTA
     $           )

         IF( ABS(COBSTA).EQ.1 .OR. ABS(COBSTA).EQ.2 )THEN
*           Eigenvalue-"Undeflow". If possible, apply small fudge
*           to the shift and try again.
            IF( DIR.NE.0 .AND.
     $          FUDGEFAC.GT.0 .AND. NFUDGE.LT.MAXFUDGE )
     $      THEN
*              Reinit as new candidate, only increase its fudge
*              counter.
               TAU = TAU + ABS(TAU)*(DIR*FUDGEFAC*EPS)
               IWORK(IYCTWI + ICAND) = 0
               IWORK(IYCFDG + ICAND) = NFUDGE + 1
*              Normally we do expect that a small fudge won't change
*              the evaluation that much. Nevertheless, treating the
*              candidate again as new provides full flexibility and
*              avoids that the fudged shift becomes worse than
*              another candidate but is still preferred.
               XNBIS_WASTED = XNBIS_WASTED + (XNUMFN - BISTRYCNT)
               CYCLE candloop
            ENDIF
         ENDIF
         IF( COBSTA .NE. 0 )THEN
               XNBIS_WASTED = XNBIS_WASTED + (XNUMFN - BISTRYCNT)
            IWORK(IYCLOC + ICAND) = 0
            CYCLE candloop
         ENDIF
*
*        ----------------------------------------------------------------
*                     Init Inner Bounds from Father
*        ----------------------------------------------------------------
         CALL DLAXRF_IIB(
     $          N, SREPR, SREPI, TAU, ICBEG, ICEND,
     $          FEWL_AE, FEWL_LU,
     $          SLGAP, SUGAP, SEWL_AE, SEWL_LU
     $       )

C         IF( IIBSTA .NE. 0 )THEN
C            Previously we discarded candidates where for some subclusters
C            of the father no corresponding bound for the son could be
C            found. This was switched off (its main use was for checking
C            the shift-relation in the coupled case, also it required
C            using the cluster internal gaps, and causes problems with
C            parallelization).
C
C*           Failed to obtain inner bounds, which means something is
C*           not right with the shift relation.
C@LOGGING on
C            WRITE(FIDRRF,*)' IIB failed, killing candidate',ICAND
C@LOGGING !
C@STATS on
C            XNBIS_WASTED = XNBIS_WASTED + (XNUMFN - BISTRYCNT)
C@STATS !
C            IWORK(IYCLOC + ICAND) = 0
C            CYCLE candloop
C         ENDIF

*        ----------------------------------------------------------------
*              [Optional]   See if the child has a gap
*        ----------------------------------------------------------------
         IF( FORCEGAP )THEN
            DO
               XA = SEWL_LU(2*ICBEG-1)
               XB = SEWL_LU(2*ICBEG)
               XC = SEWL_LU(2*ICEND-1)
               XD = SEWL_LU(2*ICEND)
*              Note that we may have XA=XC, XB=XD

               IF( (XC-XB) .GE. WNTRW*MAX(ABS(XB),ABS(XC)) )THEN
*                 Jep, wide enough, take it
                  EXIT
               ELSEIF( (XD-XA) .LT. WNTRW*MAX(ABS(XA),ABS(XD)) )THEN
*                 No this child's ews are not wide enough and will
*                 never be.
                  IF( DOGRT .AND. .NOT.ISGAPRETRY )THEN
*                    Initiate a gap retry
                     IWORK(IYCGRT + ICAND) = 1
                     IWORK(IYCTWI + ICAND) = 0
                     IWORK(IYCFDG + ICAND) = 0
                  ELSE
*                    Discard the candidate
                     IWORK(IYCLOC + ICAND) = 0
                  ENDIF
                  XNBIS_WASTED = XNBIS_WASTED + (XNUMFN - BISTRYCNT)
                  CYCLE candloop
               ELSE
*                 We cannot decide yet, bisect the larger interval.
                  IF( (XB-XA) .GE. (XD-XC) )THEN
                     I = ICBEG
                     J = SEWL_AE(2*I)
                     MID = HALF*(XA + XB)
                  ELSE
                     J = ICEND
                     I = SEWL_AE(2*J-1)
                     MID = HALF*(XC + XD)
                  ENDIF
                  XI = DLAXRN( N, SREPR, SREPI, MID )
                  XNBIS_CLASS = XNBIS_CLASS + 1
                  CALL DLAXRL_REFINE(
     $                    ICBEG, ICEND, SLGAP, SUGAP, SEWL_AE, SEWL_LU,
     $                    I, J, MID, XI
     $                 )
               ENDIF
            ENDDO
         ENDIF

*        ----------------------------------------------------------------

*        At this point the candidate passed all tests.
*        The son's data (SEWL,...) and TAU are set already.

         SUCCESS = .TRUE.
         EXIT
      ENDDO candloop

      INFO = 0
      IF( .NOT. SUCCESS )THEN
         INFO = 2
      ENDIF


C@LOGGING on
C 99         FORMAT(1X,I5,1X,@(FL))
C            IF( N .LE. MAXLOGDIM )THEN
C               WRITE(FIDRRF,*) ' SUCCESS'
C               WRITE(FIDRRF,*) ' Bounds after COB:'
C               WRITE(FIDRRF,*) '  left gap  = ', SLGAP
C               I = ICBEG
C               DO
C                  J = SEWL_AE(2*I)
C                  WRITE(FIDRRF,FMT=99)  I, SEWL_LU(2*I-1)
C                  WRITE(FIDRRF,FMT=99)  J, SEWL_LU(2*J)
C                  IF( J .EQ. ICEND ) EXIT
C                  I = J+1
C                  WRITE(FIDRRF,*)
C     $              '  [abs gap =',SEWL_LU(2*I-1)-SEWL_LU(2*J)
C               ENDDO
C               WRITE(FIDRRF,*) '  right gap = ', SUGAP
C            ENDIF
C@LOGGING !
      END SUBROUTINE DLAXRF
*
************************************************************************
