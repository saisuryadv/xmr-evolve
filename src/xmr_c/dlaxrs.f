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
*
*  Purpose
*  =======
*
*     Factorize N G N' - TTAU = Np Gp Np', where both NGN and NpGpNp' are
*     twisted block factorizations, and TTAU is a diagonal matrix.
*     The union of data for all twists in IA:IE in the target is computed.
*   ! At the moment progressive par only implemented for non-blocked !
*
*     There are two reasons why some twist may not be ok: a block would
*     end there (we don't allow that for now), or we cannot compute the
*     twist element stably. 1 and N are always ok.
*
*     Post conditions
*     ---------------
*       if IA = 1 then TWISTOK(1) = 1
*       if IE = N then TWISTOK(N) = 1
*       OMGADP(1) = 0
*       OMGARP(N) = 0
*
*  ======================================================================
*
*     .. Parameters ..
*

*
*     PIVTAUFAC
*         Adjust PIVMIN to PIVTAUFAC*EPS*TAU, set to zero to deactivate
*         a shift-dependent pivmin.
*
      DOUBLE PRECISION , PARAMETER  ::  PIVTAUFAC = 0.0D0


*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRS_STAT(
     $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, SHFPRT, GPLUS, OPLUS, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N)
*
      INTEGER,          INTENT(OUT)  ::  OPLUS(N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N)
      END SUBROUTINE DLAXRS_STAT
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRS_PROG(
     $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, SHFPRT, OFF, GPLUS, OPLUS, P
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N)
*
      INTEGER,          INTENT(OUT)  ::  OPLUS(N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N)
      END SUBROUTINE DLAXRS_PROG
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, PIVBASE, PIVMIN, BIGG
      DOUBLE PRECISION  ABSTAU, TN, TP, GMA
      INTEGER           IXG, IXBDET, IXNGN, IXGNSQ, IYOMGA, IYLBBK
      INTEGER           K, NB, IXTN, IXTP, I, J1, J2, IBB


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
*  ----- Executable Statements -----------------------------------------
*
      IF( .NOT. XSTEALTHMODE )THEN
         IF( IA .NE. IE )THEN
            XNUMFS_2 = XNUMFS_2 + 1
         ELSE
            XNUMFS_K = XNUMFS_K + 1
         ENDIF
      ENDIF

*
*     .. Decode Representation Data ..
*
      IXG    = (0)
      IXBDET = (N)
      IXNGN  = (2*N+1)
      IXGNSQ = (3*N+2)
      IYOMGA = (4)
      IYLBBK = (N+6)

      PIVBASE = REPR(4*N+3)
      K       = REPI(2)
      NB      = REPI(3)

      EPS  = DLAMCH('Epsilon')
      BIGG = DLAMCH('Overflow')

      ABSTAU = ABS(TAU)

*     Restrict range for twists if we have blocks ..
*     (progressive only implemented for no blocks in source yet)
      J1 = IA
      J2 = IE
      IF( NB .GT. 0 )THEN
*        We cannot do a progressive part that would range over a
*        block in the source.
         DO I = 0, NB
            IBB = REPI(IYLBBK + I)
            IF( IBB .LT. K )THEN
               J1 = MAX( J1, MIN(K, IBB+2) )
            ENDIF
            IF( IBB .GT. K )THEN
               J2 = MIN( J2, MAX(K, IBB-2) )
            ENDIF
         ENDDO
      ENDIF

      PIVMIN = MAX( PIVBASE, (PIVTAUFAC * EPS) * ABSTAU )

      IXTN = 1
      IXTP = N+1

      CALL DLAXRS_STAT(
     $       N, MIN(K,J2), 1, N-1,
     $       REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $       REPR(IXNGN+1), REPR(IXBDET+1),
     $       PIVMIN, TAU, SHFPRT, DPLUS, OMGADP, RWORK(IXTN)
     $     )
      IF( K.LT.J2 )THEN
         TN = RWORK(IXTN-1 + K)
         CALL DLAXRS_PROG(
     $          N, K, J2, 1, N-1,
     $          REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, SHFPRT, TN, DPLUS, OMGADP, RWORK(IXTN)
     $        )
      ENDIF

      CALL DLAXRS_STAT(
     $       N, MAX(K,J1), 2, N,
     $       REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $       REPR(IXNGN+1), REPR(IXBDET+1),
     $       PIVMIN, TAU, SHFPRT, RPLUS, OMGARP, RWORK(IXTP)
     $     )
      IF( J1.LT.K )THEN
         TP = RWORK(IXTP-1 + K)
         CALL DLAXRS_PROG(
     $          N, K, J1, 2, N,
     $          REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, SHFPRT, TP,
     $          RPLUS, OMGARP, RWORK(IXTP)
     $        )
      ENDIF
*
*     Compute gammas
*
*     Notes:
*     (1)  If one of the parts places a block end at i, 1 < i < n, then
*          we do not want to twist there, as that would again expose
*          the elg that was hidden by the block.
*      ==> it may well be that no twist is ok
*     (2)  For twisting at i where a block ends in the source, the
*          same special handling as in dlaxrt does apply.
*
      DO IXF77A = IA, IE
         TWISTOK(IXF77A) = 0
      ENDDO

      DO I = J1, J2
         GMA = BIGG

         IF( (I.NE.N .AND. OMGADP(I).NE.0) .OR.
     $       (I.NE.1 .AND. OMGARP(I).NE.0) )
     $   THEN
            TWISTOK(I) = 0
         ELSEIF( REPI(IYOMGA + I).EQ.0 .OR. I.EQ.1 .OR. I.EQ.N )THEN
*           neither source nor target have a block here
            TN = RWORK(IXTN-1 + I)
            TP = RWORK(IXTP-1 + I)
            IF( I.NE.K )THEN
               GMA = (TN + TP) - TAU*SHFPRT(I)
            ELSE
               GMA = REPR(IXG + K) + ((TN + TP) - TAU*SHFPRT(I))
            ENDIF
            TWISTOK(I) = 1
         ELSE
*           block ends in the source
C            TODO, see dlaxrt
            TWISTOK(I) = 0
         ENDIF
         IF( TWISTOK(I) .EQ. 1 )THEN
            GAMMAP(I) = GMA
         ENDIF
      ENDDO

C      DO I = IA,IE
C         IF( (I.LT.J1) .OR. (I.GT.J2) .OR.
C     $       (I.NE.N .AND. OMGADP(I).NE.0) .OR.
C     $       (I.NE.1 .AND. OMGARP(I).NE.0) )
C     $   THEN
C            TWISTOK(I) = 0
C         ELSE
C*           normal case (no block ending)
C            TN = RWORK(IXTN-1 + I)
C            TP = RWORK(IXTP-1 + I)
C            GMA = (TN + TP) - TAU*SHFPRT(I)
C            IF( I.EQ.K )THEN
C               GMA = GMA + REPR(IXG + K)
C            ENDIF
C            GAMMAP(I)  = GMA
C            TWISTOK(I) = 1
C         ENDIF
C      ENDDO

      END SUBROUTINE DLAXRS
*
************************************************************************
