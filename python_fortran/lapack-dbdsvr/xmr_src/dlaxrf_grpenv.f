      SUBROUTINE DLAXRF_GRPENV(
     $            N, E, REPR, REPI, I, J,
     $            LGAP, LB, UB, UGAP,
     $            ENV,
     $            RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, I, J
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  LGAP, LB, UB, UGAP
      DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*N )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 6*N )
*
      DOUBLE PRECISION, INTENT(OUT)  ::  ENV(N)
*
*
*  Purpose
*  =======
*
*     Compute envelope for (sub-)cluster I:J.
*     It is *not* checked beforehand if it makes sense to try for an
*     envelope here, that should have been done already.
*     We use the Parlett/Voemel envelope strategy, or a plain
*     twisted facto for a singleton subcluster.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
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
      END SUBROUTINE DLAXRT
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRG0(
     $             N, E, REPR, REPI, CUTTOL, Z
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)   ::  N
      INTEGER,          INTENT(IN)   ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)   ::  REPR(4*N+3)
      DOUBLE PRECISION, INTENT(IN)   ::  E(N-1)
      DOUBLE PRECISION, INTENT(IN)   ::  CUTTOL
*
      DOUBLE PRECISION, INTENT(OUT)  ::  Z(N)
      END SUBROUTINE DLAXRG0
      SUBROUTINE DLAXRG(
     $             N, K, D, R, E, GAMMA, CUTTOL,
     $             Z, NORMZ, RESID, RQCORR
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K
      DOUBLE PRECISION, INTENT(IN)  ::  GAMMA, CUTTOL
      DOUBLE PRECISION, INTENT(IN)  ::  D(1:K-1), R(K+1:N), E(1:N-1)
*
      DOUBLE PRECISION, INTENT(OUT)  ::  NORMZ, RESID, RQCORR
      DOUBLE PRECISION, INTENT(OUT)  ::  Z(N)
      END SUBROUTINE DLAXRG
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH


*
*     .. Constants ..
*
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER       ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )

*
*     .. Parameters ..
*

*         Intervals will be refined via bisection up to a relative
*         accuracy of BISACCFAC*N*EPS
      INTEGER, PARAMETER  ::  BISACCFAC = 4

      INTEGER, PARAMETER  ::  TAUFAC = 2

      INTEGER, PARAMETER  ::  MAXRELCOND = 10



*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, BISACC, ABSGMA, MINGAP, RTMP
      DOUBLE PRECISION  NORMZ, RESID, RQCORR, SINBND
      DOUBLE PRECISION  MU, DELTA, BCKOFF, INVRHO, SL, SU, XENV
      DOUBLE PRECISION  INVMINRGAP, ENVFORCE
      INTEGER           JXTOKL, JXTOKU, JXWORK
      INTEGER           IXDP, IXRP, IXGMAL, IXGMAU, IXWORK
      INTEGER           TWIST, K


*     [XMRSTATS COMMON block removed for thread safety]
*
*  ===== Executable Statements ==========================================
*
      EPS = DLAMCH('Precision')
      BISACC = BISACCFAC * N * EPS

      MINGAP = MIN( LGAP, UGAP )


*
*     Integer Workspace
*
      JXTOKL = 1
      JXTOKU = JXTOKL + N
      JXWORK = JXTOKU + N
*       ---------
*              2N

*
*     Real Workspace
*
      IXDP   = 1
      IXRP   = IXDP   + N
      IXGMAL = IXRP   + N
      IXGMAU = IXGMAL + N
      IXWORK = IXGMAU + N
*       ---------
*              4N
*            + 2N for DLAXRT
*       =========
*              6N

*     --------------------
*      Build the envelope
*     --------------------

      IF( I .EQ. J )THEN
*        ----------------------
*         Singleton Subcluster
*        ----------------------
         CALL DLAXRT(
     $          N, 1, N, REPR, REPI, HALF*(LB+UB),
     $          RWORK(IXDP), RWORK(IXRP),
     $          RWORK(IXGMAL), IWORK, RWORK(IXWORK)
     $        )

*        get min gamma
         IWORK(N+1) = -1
         TWIST = 1
         DO
            IF( IWORK(TWIST) .NE. 0 )  EXIT
            TWIST = TWIST+1
         ENDDO
         ABSGMA = ABS(RWORK(IXGMAL-1 + TWIST))
         DO K = TWIST+1,N
            IF( IWORK(K) .NE. 0 )THEN
               RTMP = ABS(RWORK(IXGMAL-1 + K))
               IF( RTMP.LT.ABSGMA )THEN
                  ABSGMA = RTMP
                  TWIST  = K
               ENDIF
            ENDIF
         ENDDO

*        Compute eigenvector approximation
         CALL DLAXRG(
     $          N, TWIST, RWORK(IXDP), RWORK(IXRP-2 + TWIST+1), E,
     $          RWORK(IXGMAL-1 + TWIST), ZERO, RWORK(IXWORK),
     $          NORMZ, RESID, RQCORR
     $        )
         SINBND = RESID / MINGAP

*        Extract envelope
         DO K = 1, N
            ENV(K) = MIN( ONE, ABS(RWORK(IXWORK-1+K)) + SINBND )
         ENDDO

C@LOGGING on
C         WRITE(FIDRRF,*) '   -- Singleton --'
C         WRITE(FIDRRF,*) '     absgma =',ABSGMA
C         WRITE(FIDRRF,*) '     twist  =',TWIST
C         WRITE(FIDRRF,*) '     sinbnd =',SINBND
C@LOGGING !
      ELSE
*        ----------
*         Subgroup
*        ----------
*        We use one parameter BCKOFF in place of TAU*DELTA in the
*        Parlett/Voemel paper. To avoid zero gammas we backoff at
*        least a bit.
         MU    = HALF*(LB + UB)
         DELTA = HALF*(UB - LB)
         BCKOFF = MAX( TAUFAC*DELTA,
     $                 DELTA + BISACC*MAX(ABS(LB),ABS(UB)) )
         INVRHO = ONE / ( ((MINGAP+DELTA)/BCKOFF)**2 - ONE )
         SL    = MU - BCKOFF
         SU    = MU + BCKOFF
*        Compute twisted factorizations at both ends.
*        Only the gammas are needed.
         CALL DLAXRT(
     $          N, 1, N, REPR, REPI, SL, RWORK(IXDP), RWORK(IXRP),
     $          RWORK(IXGMAL), IWORK(JXTOKL), RWORK(IXWORK)
     $        )

         CALL DLAXRT(
     $          N, 1, N, REPR, REPI, SU, RWORK(IXDP), RWORK(IXRP),
     $          RWORK(IXGMAU), IWORK(JXTOKU), RWORK(IXWORK)
     $        )

         DO K = 1, N
            IF( IWORK(JXTOKL-1 + K).EQ.0     .OR.
     $          IWORK(JXTOKU-1 + K).EQ.0     .OR.
     $          RWORK(IXGMAL-1 + K).EQ.ZERO  .OR.
     $          RWORK(IXGMAU-1 + K).EQ.ZERO
     $        )
     $      THEN
               XENV = ONE
            ELSE
*              The smallest weight inside the cluster is at least
*              2 / BCKOFF.
               RTMP = (HALF * BCKOFF)
               RTMP = ABS(  RTMP / RWORK(IXGMAL-1 + K)
     $                    + RTMP / RWORK(IXGMAU-1 + K) )
     $                + (N - (J-I+1)) * INVRHO
               XENV = MIN( ONE, SQRT(RTMP) )
            ENDIF
            ENV(K) = XENV
         ENDDO

      ENDIF
*
*     Safeguard by forcing a minimal envelope entry. The rationale
*     here is that for a small change in the representation, a
*     componentwise change of
*       n eps relcond / relgap
*     has to be anticipated in each vector, so the envelope entries
*     really cannot fall below that.
*
      INVMINRGAP = MAX( ABS(LB) / LGAP, ABS(UB) / UGAP )
      ENVFORCE   = MIN( ONE, (MAXRELCOND*N*EPS)*INVMINRGAP )
      DO K = 1, N
         ENV(K) = MAX( ENVFORCE, ENV(K) )
      ENDDO

C@LOGGING on
C      ENVMIN = ONE
C      ENVMAX = ZERO
C      ENVAVG = ZERO
C      DO K = 1, N
C         ENVMIN = MIN(ENVMIN, ENV(K))
C         ENVMAX = MAX(ENVMAX, ENV(K))
C         ENVAVG = ENVAVG + ENV(K)
C      ENDDO
C      ENVAVG = ENVAVG / N
C      WRITE(FIDRRF,*) '   -- Summary --'
C      WRITE(FIDRRF,*) '     envforce = ',ENVFORCE
C      WRITE(FIDRRF,*) '     min = ',ENVMIN
C      WRITE(FIDRRF,*) '     avg = ',ENVAVG
C      WRITE(FIDRRF,*) '     max = ',ENVMAX
C      WRITE(FIDRRF,*) ' '
C@LOGGING !

*

      IF( I .EQ. J )THEN
      ELSE
      ENDIF

      END SUBROUTINE DLAXRF_GRPENV
*
************************************************************************
