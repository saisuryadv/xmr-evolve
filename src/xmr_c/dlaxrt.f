      SUBROUTINE DLAXRT(
     $             N, J1, J2, REPR, REPI,
     $             TAU, DPLUS, RPLUS, GAMMAP, TWISTOK, RWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J1, J2
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  TAU
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
*
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK(2*N)
*
      INTEGER,          INTENT(OUT)  ::  TWISTOK(J1:J2)
      DOUBLE PRECISION, INTENT(OUT)  ::  DPLUS(1:N-1), RPLUS(2:N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GAMMAP(J1:J2)
*
*  Purpose
*  =======
*
*     Factorize twisted blocked to twisted non-blocked, for all twists
*     in J1:J2.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRT_STAT(
     $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, GPLUS, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N)
      END SUBROUTINE DLAXRT_STAT
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRT_PROG(
     $             N, J1, J2, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, OFF, GPLUS, P
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J1, J2, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU, OFF
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), P(N)
      END SUBROUTINE DLAXRT_PROG
      END INTERFACE

      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION ZERO, ONE, FIVE, TEN
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, FIVE=5.0D0, TEN=10.0D0)

*
*     .. Parameters ..
*

*
*     PIVTAUFAC -- see laxrs  (might make sense to use different values)
*
      DOUBLE PRECISION , PARAMETER  ::  PIVTAUFAC = 0.0D0

*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, PIVBASE, PIVMIN, BIGG
      DOUBLE PRECISION  TN, TP, GMA
      DOUBLE PRECISION  TMPA, TMPB, TMPF, TMPG1, TMPG2, TMPH
      DOUBLE PRECISION  PMTAU, SMTAU
      INTEGER           I, IXTN, IXTP, TWOK
      INTEGER           K, IXG, IXBDET, IXNGN, IXGNSQ, IYOMGA
      LOGICAL           DOGMASTD

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
*  ----- Executable Statements -----------------------------------------
*
      IF( .NOT. XSTEALTHMODE )  XNUMFT = XNUMFT + 1
*
*     .. Decode Representation Data ..
*
      IXG    = (0)
      IXBDET = (N)
      IXNGN  = (2*N+1)
      IXGNSQ = (3*N+2)
      IYOMGA = (4)

      PIVBASE = REPR(4*N+3)
      K       = REPI(2)

      EPS  = DLAMCH('Epsilon')
      BIGG = DLAMCH('Overflow')


      PIVMIN = MAX( PIVBASE, (PIVTAUFAC * EPS) * ABS(TAU) )
      IXTN = 1
      IXTP = N+1

      CALL DLAXRT_STAT(
     $       N, MIN(K,J2), 1, N-1,
     $       REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $       REPR(IXNGN+1), REPR(IXBDET+1),
     $       PIVMIN, TAU, DPLUS, RWORK(IXTN)
     $     )
      IF( K.LT.J2 )THEN
         TN = RWORK(IXTN-1 + K)
         CALL DLAXRT_PROG(
     $          N, K, J2, 1, N-1,
     $          REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, TN, DPLUS, RWORK(IXTN)
     $        )
      ENDIF

      CALL DLAXRT_STAT(
     $       N, MAX(K,J1), 2, N,
     $       REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $       REPR(IXNGN+1), REPR(IXBDET+1),
     $       PIVMIN, TAU, RPLUS, RWORK(IXTP)
     $     )
      IF( J1.LT.K )THEN
         TP = RWORK(IXTP-1 + K)
         CALL DLAXRT_PROG(
     $          N, K, J1, 2, N,
     $          REPR(IXG+1), REPI(IYOMGA), REPR(IXGNSQ+1),
     $          REPR(IXNGN+1), REPR(IXBDET+1),
     $          PIVMIN, TAU, TP, RPLUS, RWORK(IXTP)
     $        )
      ENDIF
*
*     Compute gammas
*

      IF( J1.LE.K .AND. K.LE.J2 )THEN
*        k=1 or k=n also handled by this
         TN = RWORK(IXTN-1 + K)
         TP = RWORK(IXTP-1 + K)
         GAMMAP(K)  = (REPR(IXG + K) + (TN + TP)) - TAU
         TWISTOK(K) = 1
      ENDIF

*     Notes:
*     (1)  The following requires safe pivmin atm, otherwise
*          use safe mults
*     (2)  We access the next GAMMAP(i+1) as P before overwriting
*          it, so
*          DO NOT CHANGE THE LOOP ORDERS

*
*     Upper Part (analogous to below)
*
      DO I = J1, MIN(J2,K-1)
*        we set gma just to soothe compilers
         GMA  = BIGG
         TWOK = 0

         DOGMASTD = (REPI(IYOMGA + I).EQ.0 .OR. I.EQ.1 .OR. I.EQ.N)
         IF( .NOT. DOGMASTD )THEN
            DOGMASTD = ( REPR(IXG + I-1).EQ.ZERO )
         ENDIF

         IF( DOGMASTD )THEN
*           normal case (no inner block ending) or
*           block with zero d
            TN = RWORK(IXTN-1 + I)
            TP = RWORK(IXTP-1 + I)
            GMA = (TN + TP) - TAU
            TWOK = 1
         ELSE
*           twist where the source has a block-end

            IF( I+1.EQ.K )THEN
               PMTAU = RWORK(IXTP-1 + I+1) + (REPR(IXG + K) - TAU)
            ELSE
               PMTAU = RWORK(IXTP-1 + I+1) - TAU
            ENDIF
            SMTAU = RWORK(IXTN-1 + I-1) - TAU

            TMPA = ( REPR(IXBDET + I-1) * PMTAU ) /
     $             ( REPR(IXG + I-1) * RPLUS(I+1) )
            TMPB = ( REPR(IXGNSQ + I-1) * SMTAU ) /
     $             ( REPR(IXG + I-1) * DPLUS(I-1) )
C            TMPG1 = REPR(IXNGN + I) / RPLUS(I+1)
C            TMPG2 = REPR(IXG + I-1) / DPLUS(I-1)
            IF( (ABS(TMPA)+ABS(TMPB) .LE. TEN*ABS(TMPA+TMPB)) .OR.
     $           MIN(ABS(TMPA),ABS(TMPB)) .LE. FIVE*ABS(TAU) )
     $      THEN
*              method 1
               GMA = (TMPA + TMPB) - TAU
               TWOK = 1
C something is wrong with method 2, but turned out to be not absolutely
C required anyway
C            ELSEIF( ABS(TMPG1)+ABS(TMPG2) .LE. FIVE*(TMPG1-TMPG2) )
C     $      THEN
C*              method 2
C               TMPF  = REPR(IXG + I) * (PMTAU / RPLUS(I+1))
C               TMPH  = REPR(IXGNSQ + I-1) / REPR(IXG + I-1)
C               GMA = ( TMPF + TMPH*(TMPG1-TMPG2) ) - TAU
C               TWOK = 1
            ENDIF
            IF( TWOK .EQ. 0 )THEN
               TMPG1 = RWORK(IXTP-1 + I) - TAU
*              now TMPG1 = RP(i) as computed by progressive part
               TMPG2 = REPR(IXGNSQ + I-1) / DPLUS(I-1)
               TMPF  = TMPG1 - TMPG2
               IF( ABS(TMPG1)+ABS(TMPG2) .LE. TEN*ABS(TMPF) ) THEN
                  GMA = TMPF
                  TWOK = 1
               ENDIF
            ENDIF
         ENDIF
         TWISTOK(I) = TWOK
         IF( TWOK .EQ. 1 )THEN
            GAMMAP(I) = GMA
         ENDIF
      ENDDO

*
*     Lower Part (analogous to above)
*
      DO I = MAX(J1,K+1), J2
*        set gma just to soothe compilers
         GMA = BIGG
         TWOK = 0

         DOGMASTD = (REPI(IYOMGA + I).EQ.0 .OR. I.EQ.1 .OR. I.EQ.N)
         IF( .NOT. DOGMASTD )THEN
            DOGMASTD = ( REPR(IXG + I+1).EQ.ZERO )
         ENDIF

         IF( DOGMASTD )THEN
*           normal case (no block ending) or
*           block with zero r
            TN = RWORK(IXTN-1 + I)
            TP = RWORK(IXTP-1 + I)
            GMA = (TN + TP) - TAU
            TWOK = 1
         ELSE
*           twist where the source has a block-end

            IF( I-1.EQ.K )THEN
               PMTAU = RWORK(IXTP-1 + I-1) + (REPR(IXG + K) - TAU)
            ELSE
               PMTAU = RWORK(IXTP-1 + I-1) - TAU
            ENDIF
            SMTAU = RWORK(IXTN-1 + I+1) - TAU

            TMPA = ( REPR(IXBDET + I+1) * PMTAU ) /
     $             ( REPR(IXG + I+1) * DPLUS(I-1) )
            TMPB = ( REPR(IXGNSQ + I+1) * SMTAU ) /
     $             ( REPR(IXG + I+1) * RPLUS(I+1) )
C            TMPG1 = REPR(IXNGN + I) / DPLUS(I-1)
C            TMPG2 = REPR(IXG + I+1) / RPLUS(I+1)
            IF( (ABS(TMPA)+ABS(TMPB) .LE. TEN*ABS(TMPA+TMPB)) .OR.
     $           MIN(ABS(TMPA),ABS(TMPB)) .LE. FIVE*ABS(TAU) )
     $      THEN
*              method 1
               GMA = (TMPA + TMPB) - TAU
               TWOK = 1
C see above
C            ELSEIF( ABS(TMPG1)+ABS(TMPG2) .LE. FIVE*(TMPG1-TMPG2) )
C     $      THEN
C*              method 2
C               TMPF  = REPR(IXG + I) * (PMTAU / DPLUS(I-1))
C               TMPH  = REPR(IXGNSQ + I+1) / REPR(IXG + I+1)
C               GMA = ( TMPF + TMPH*(TMPG1-TMPG2) ) - TAU
C               TWOK = 1
            ENDIF
            IF( TWOK .EQ. 0 )THEN
               TMPG1 = RWORK(IXTP-1 + I) - TAU
*              now TMPG1 = RP(i) as computed by progressive part
               TMPG2 = REPR(IXGNSQ + I+1) / RPLUS(I+1)
               TMPF  = TMPG1 - TMPG2
               IF( ABS(TMPG1)+ABS(TMPG2) .LE. TEN*ABS(TMPF) ) THEN
                  GMA = TMPF
                  TWOK = 1
               ENDIF
            ENDIF
         ENDIF
         TWISTOK(I) = TWOK
         IF( TWOK .EQ. 1 )THEN
            GAMMAP(I) = GMA
         ENDIF
      ENDDO

*
      END SUBROUTINE DLAXRT
*
************************************************************************
