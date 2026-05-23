      FUNCTION DLAXRN(N, REPR, REPI, TAU)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN
      INTEGER,          INTENT(IN)  ::  N
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU
*
*  Purpose
*  =======
*
*    Perform a sturm count by doing a (suitably twisted) bidiagonal
*    factorization of T - TAU and then count the negative pivots.
*    However, for utmost information the function does not return the
*    negcount, but
*     xi := 2*negc + issing
*    where issing=1 if T-TAU is singular, and 0 otherwise. Thus it pro-
*    vides full information about the inertia.
*
*    Consequences / Usage
*    --------------------
*         xi = 2k-1 (odd)  means we are on ew k
*         xi = 2k   (even) means we are in interval (k,k+1)
*
*    TAU is
*      upper bound for ew k         <=>  xi >= 2*k-1
*      strict upper bound for ew k  <=>  xi >  2*k-1
*
*      lower bound for ew k         <=>  xi <= 2*k-1
*      strict lower bound for ew k  <=>  xi <  2*k-1
*
*      in (ew k, ew k+1)  <=>  xi = 2k
*      in [ew k, ew k+1)  <=>  xi = 2k-1 or xi = 2k
*      in (ew k, ew k+1]  <=>  xi = 2k   or xi = 2k+1
*      in [ew k, ew k+1]  <=>  2k-1 <= xi <= 2k+1
*
*    For given inertia xi and sample lambda,
*      the nearest eigenvalue smaller than (not equal to) lambda is
*        xi / 2,
*      the nearest eigenvalue larger than (not equal to) lambda is
*        xi / 2 + MOD(xi,2) + 1  =  (xi+1) / 2 + 1
*
*  ======================================================================
*
*
*     .. Parameters ..
*

*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRN_STAT(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, NCOUNT, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  NCOUNT
      DOUBLE PRECISION, INTENT(OUT)  ::  S
      END SUBROUTINE DLAXRN_STAT
      END INTERFACE
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  PIVBASE, PIVMIN, SN, SP, GMA
      INTEGER           NEGN, NEGP, NEG, IBB, I, K, NB
      INTEGER           IXG, IXBDET, IXNGN, IXGNSQ, IYLBBK, IYOMGA


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
*
*  ----- Executable Statements ------------------------------------------
*
      IF( .NOT. XSTEALTHMODE )  XNUMFN = XNUMFN + 1
*
*     .. Decode Representation Data ..
*
      IXG    = (0)
      IXBDET = (N)
      IXNGN  = (2*N+1)
      IXGNSQ = (3*N+2)
      IYLBBK = (N+6)
      IYOMGA = (4)

      PIVBASE = REPR(4*N+3)
      K       = REPI(2)
      NB      = REPI(3)

*
*     .. Special treatment for zero inertia ..
*
      IF( TAU .EQ. ZERO )THEN
*        NOTE: The source may have blocks, so here we need to
*              count them correctly (!)
         NEG = 0
*
         IBB = 0
         I   = 1
         DO
            DO
               IF( I .GE. REPI(IYLBBK + IBB) )  EXIT
               IF( REPR(IXG+I) .LT. ZERO )  NEG = NEG+1
               I = I+1
            ENDDO
            IF( I .GE. K )  EXIT
            IBB = IBB+1
            NEG = NEG+1
            I   =  I+2
         ENDDO
*
         IBB = NB
         I   = N
         DO
            DO
               IF( I .LE. REPI(IYLBBK + IBB) )  EXIT
               IF( REPR(IXG+I) .LT. ZERO )  NEG = NEG+1
               I = I-1
            ENDDO
            IF( I .LE. K )  EXIT
            IBB = IBB-1
            NEG = NEG+1
            I   = I-2
         ENDDO
*
         DLAXRN = 2 * NEG
         IF( REPI(IYOMGA + K) .EQ. 0 )THEN
            IF( REPR(IXG+K) .LT. ZERO )  DLAXRN = DLAXRN + 2
            IF( REPR(IXG+K) .EQ. ZERO )  DLAXRN = DLAXRN + 1
         ENDIF
         RETURN
*
      ENDIF
*
*     General case, tau is nonzero
*
      PIVMIN = PIVBASE

      IBB = 1
      NEGN = 0
      SN  = ZERO
      IF( K.GT.1 )THEN
         CALL DLAXRN_STAT(
     $          N, K, +1,
     $          REPR(IXG+1), IBB, REPI(IYLBBK), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, NEGN, SN )
      ENDIF
*
      IBB = NB+1
      NEGP = 0
      SP  = ZERO
      IF( K.LT.N )THEN
         CALL DLAXRN_STAT(
     $          N, K, -1,
     $          REPR(IXG+1), IBB, REPI(IYLBBK), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, NEGP, SP )
      ENDIF
*
      GMA = REPR(IXG + K) + ((SN + SP) - TAU)
*
      NEG = NEGN + NEGP
      IF( GMA.LT.ZERO )  NEG = NEG + 1

      DLAXRN = 2 * NEG
      IF( GMA.EQ.ZERO )  DLAXRN = DLAXRN + 1
*
      END FUNCTION DLAXRN
*
************************************************************************
