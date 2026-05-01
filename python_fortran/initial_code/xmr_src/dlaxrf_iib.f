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
*
*  Purpose
*  =======
*
*    Use the father's inner bounds in FEWL to initialize bounds for
*    the son.
*    For each father bound of an interval, a couple of relative
*    relaxations are tried and the shift is applied to get a bound
*    for the son.
*
*  ======================================================================
*
*     .. Declarations ..
*
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
      INTERFACE
      FUNCTION DLAXRN(N, REPR, REPI, TAU)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN
      INTEGER,          INTENT(IN)  ::  N
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU
      END FUNCTION DLAXRN
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0

      INTEGER, PARAMETER  ::  NRELAX = 3

*
*     .. Parameters ..
*

*         Intervals will be refined via bisection up to a relative
*         accuracy of BISACCFAC*N*EPS
      INTEGER, PARAMETER  ::  BISACCFAC = 4

*     Father and son bounds are relaxed before checking them by
*         2**LOGFBNDRELAXFAC * (N*EPS)
*
      INTEGER, PARAMETER  ::  LOGFBNDRELAXFAC = 7

*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, ARXFAC(NRELAX), FBOUND, SBOUND
      INTEGER           I, J, IRXFAC, XI, ISON, JSON
      LOGICAL           GOTBND

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

*    Temporaries and Configuration
*        If true, calls to internal routines dlaxr[n,t,s] are not
*        counted. Higher level routines like DLAXRB are not affected.
*        We need this to hide bisections steps done within aux
*        routines like xmr_estrc.
         LOGICAL XSTEALTHMODE

      SAVE /XMRSTATS/
      INTEGER XNN0
*
*  ===== Executable Statements ==========================================
*
      XNN0 = XNUMFN
      EPS  = DLAMCH('Epsilon')
*
      ARXFAC(1) = 2*EPS
      ARXFAC(2) = BISACCFAC * (N * EPS)
      ARXFAC(3) = 2**LOGFBNDRELAXFAC * (N * EPS)

*     Loop over inner gaps
      J = ICBEG-1
      DO
         J = FEWL_AE(2*(J+1))
         IF( J .EQ. ICEND )  EXIT
         I = J+1
*
*        Upper bound for J
*
         GOTBND = .FALSE.
         IRXFAC = 1
         DO
            IF( GOTBND .OR. IRXFAC.GT.NRELAX )  EXIT
*
            FBOUND = FEWL_LU(2*J)
            FBOUND = FBOUND + ARXFAC(IRXFAC)*ABS(FBOUND)
            SBOUND = FBOUND - TAU
            SBOUND = SBOUND + ARXFAC(IRXFAC)*ABS(SBOUND)
            ISON   = SEWL_AE(2*J-1)
            JSON   = SEWL_AE(2*J)
*
            IF( SEWL_LU(2*J-1).LT.SBOUND .AND. SBOUND.LT.SEWL_LU(2*J) )
     $      THEN
               XI = DLAXRN( N, REPR, REPI, SBOUND )
               CALL DLAXRL_REFINE(
     $                ICBEG, ICEND, SLGAP, SUGAP, SEWL_AE, SEWL_LU,
     $                ISON, JSON, SBOUND, XI
     $              )
            ENDIF
*
            GOTBND = ( SEWL_LU(2*J) .LE. SBOUND )
            IRXFAC = IRXFAC + 1
         ENDDO
*
*        Lower bound for I
*
         GOTBND = .FALSE.
         IRXFAC = 1
         DO
            IF( GOTBND .OR. IRXFAC.GT.NRELAX )  EXIT
*
            FBOUND = FEWL_LU(2*I-1)
            FBOUND = FBOUND - ARXFAC(IRXFAC)*ABS(FBOUND)
            SBOUND = FBOUND - TAU
            SBOUND = SBOUND - ARXFAC(IRXFAC)*ABS(SBOUND)
            ISON   = SEWL_AE(2*I-1)
            JSON   = SEWL_AE(2*I)
*
            IF( SEWL_LU(2*I-1).LT.SBOUND .AND. SBOUND.LT.SEWL_LU(2*I) )
     $      THEN
               XI = DLAXRN( N, REPR, REPI, SBOUND )
               CALL DLAXRL_REFINE(
     $                ICBEG, ICEND, SLGAP, SUGAP, SEWL_AE, SEWL_LU,
     $                ISON, JSON, SBOUND, XI
     $              )
            ENDIF
*
            GOTBND = ( SBOUND .LE. SEWL_LU(2*I-1) )
            IRXFAC = IRXFAC + 1
         ENDDO
      ENDDO
      XNBIS_IIB = XNBIS_IIB + (XNUMFN - XNN0)

*
      END SUBROUTINE DLAXRF_IIB
*
*************************************************************************
