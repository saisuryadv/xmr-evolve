      SUBROUTINE WSREQ_XRX(N,REQR,REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)     ::  N
      INTEGER, INTENT(INOUT)  ::  REQR, REQI
      REQR = 7*N
      REQI = N+1
      END SUBROUTINE WSREQ_XRX
*
************************************************************************
*
      SUBROUTINE DLAXRX(
     $             N, E, REPR, REPI,
     $             INDEX, MINGAP, LAMBDA,
     $             Z, ISUPPZ, RWORK, LRWORK, IWORK, LIWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, INDEX
      INTEGER,          INTENT(IN)     ::  LRWORK, LIWORK
      INTEGER,          INTENT(IN)     ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)     ::  MINGAP
      DOUBLE PRECISION, INTENT(IN)     ::  E(N-1), REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK ), LAMBDA
*
      INTEGER,          INTENT(OUT)    ::  ISUPPZ(2)
      DOUBLE PRECISION, INTENT(OUT)    ::  Z(N)
*
*  Purpose
*  =======
*
*     Compute an eigenvector for a singleton via RQI, backed up by
*     (improved) bisection.
*
*     Note:  We do not use initial bounds for the eigenvalue, as
*            these would be based on normal bisection. Trusting them
*            can lose to much accuracy.
*            This means we implicitly assume LAMBDA to have some
*            initial accuracy, and in particular that it is placed far
*            enough from other ews so that the RQI can converge.
*            If that is not the case, the routine will detect this
*            and produce a runtime abort (in debug mode only, this
*            should become an error code).
*
*  ======================================================================
*
*
*     .. Parameters ..
*

*
*    ACCTOLFAC
*        Accept eigenpair as soon as change in eigenvalue is less than
*        ACCTOLFAC*Eps*|lambda|
*
      DOUBLE PRECISION, PARAMETER  ::  ACCTOLFAC = 2.0D0
*
*    RESTOLFAC
*        Accept eigenpair as soon as resiudal falls below
*        RESTOL=RESTOLFAC*n*Eps*mingap.
*
      DOUBLE PRECISION, PARAMETER  ::  RESTOLFAC = 2.0D0
*
*     CUTTOLFAC
*        Allowed change to the residual through cutting the support is
*        CUTTOLFAC*RESTOL.
*
      DOUBLE PRECISION, PARAMETER  ::  CUTTOLFAC = 0.5D0
*
*     PBBLEN ('bisection burst length')
*        Number of steps in a bisection burst, >= 1
      INTEGER, PARAMETER  ::  PBBLEN = 4
*
*     RESDECFAC ('residual (norm) decrease factor')
*        Factor in (0,1) by how much we would like the residual norm to
*        decrease per iteration at least.
*
      DOUBLE PRECISION, PARAMETER  ::  RESDECFAC = 0.7D0

*
*     .. Declarations ..
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

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

      DOUBLE PRECISION, PARAMETER  ::    ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  ATHIRD = 5.0D0 / 16.0D0
      DOUBLE PRECISION, PARAMETER  ::    HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::     ONE = 1.0D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  TMP, RTMPA, RTMPB, ZETA, MINGMA, ABSGMA, DIR
      DOUBLE PRECISION  LB, UB
      DOUBLE PRECISION  NORMZ, RESID, RQCORR
      DOUBLE PRECISION  ACCTOL, RESTOL, CUTTOL, EPS
      DOUBLE PRECISION  OPTLDA, OPTRES, BCKLDA, XLAMBDA

      INTEGER           IXGMA, IXDP, IXRP, IXBUF(2), IXWRK, IYTWOK
      INTEGER           I, TWIST, BISCOUNT, STATUS, BUFSEL, NTWBAD
      INTEGER           ITER, NBIS, NRQI
      INTEGER           IXF77A

      LOGICAL           HAVEACC, HAVELB, HAVEUB, HAVEOPT



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
      EPS = DLAMCH('Epsilon')

      IXGMA    = 1                !  GMA(1:N)
      IXDP     = IXGMA    + N     !  DPLUS(1:N-1)
      IXRP     = IXDP     + N-1   !  RPLUS(2:N)
      IXBUF(1) = IXRP     + N-1
      IXBUF(2) = IXBUF(1) + N
      IXWRK    = IXBUF(2) + N
*     Further workspace:  2N for DLAXRT


*     We hold two buffers to temporarily store the computed vectors,
*     indicated by IXBUFA and IXBUFB. The index IXACTBUF points to the
*     one we are currently using. Once a vector was found (HAVEOPT=
*     .TRUE.) the other buffer will hold the so-far best one.
      BUFSEL = 1

      IYTWOK = 1

      ACCTOL = ACCTOLFAC * EPS
      RESTOL = (RESTOLFAC * (N * EPS)) * MINGAP
      CUTTOL = CUTTOLFAC * RESTOL

      HAVELB = .FALSE.
      HAVEUB = .FALSE.
      HAVEOPT = .FALSE.
      OPTRES = -ONE
*     Just to soothe compilers, these settings will never be used.
      LB = ZERO
      UB = ZERO
      OPTLDA = ZERO

      BISCOUNT = 0
      HAVEACC = .FALSE.

      ITER = 0
      NBIS = 0
      NRQI = 0


*     ========================
*      Special case: Lambda=0
*     ========================
      IF( LAMBDA.EQ.ZERO )THEN
         CALL DLAXRG0( N, E, REPR, REPI, CUTTOL, Z )
         OPTLDA = LAMBDA
         OPTRES = ZERO
         GOTO 100
      ENDIF

*     ===============
*        Main Loop
*     ===============
      DO
*
*        LAMBDA is set, either LB and UB are still undefined, or LAMBDA
*        is consistent with them.
*
         CALL DLAXRT(
     $          N, 1, N, REPR, REPI, LAMBDA,
     $          RWORK(IXDP), RWORK(IXRP), RWORK(IXGMA), IWORK(IYTWOK),
     $          RWORK(IXWRK)
     $        )
         ITER = ITER+1
*
*        Find minimal gamma
*
         NTWBAD = 0
         TWIST  = 0
         ABSGMA = -ONE
         DO I = 1,N
            IF( IWORK(IYTWOK-1 + I) .NE. 0 )THEN
               TMP = ABS(RWORK(IXGMA-1 + I))
               IF( TWIST.EQ.0 .OR. TMP.LT.ABSGMA )THEN
                  ABSGMA = TMP
                  TWIST  = I
               ENDIF
            ELSE
               NTWBAD = NTWBAD + 1
            ENDIF
         ENDDO
         MINGMA = RWORK(IXGMA-1 + TWIST)
*
*        Update guidance interval
*
         IF( MINGMA .EQ. ZERO )THEN
            LB = LAMBDA
            HAVELB = .TRUE.
            UB = LAMBDA
            HAVEUB = .TRUE.
         ELSEIF( MINGMA.LT.ZERO )THEN
            UB = LAMBDA
            HAVEUB = .TRUE.
         ELSE
            LB = LAMBDA
            HAVELB = .TRUE.
         ENDIF

         IF( HAVELB .AND. HAVEUB )THEN
            HAVEACC = HAVEACC .OR. ((UB-LB) .LE. ACCTOL*ABS(LAMBDA))
         ENDIF


*
*        Cancel a bisection request if full accuracy is reached or if
*        we can see that this is a good one.
*
         IF( (BISCOUNT .GT. 0) .AND.
     $       (HAVEACC .OR. ABSGMA.LE.5*RESTOL .OR.
     $        (HAVEOPT .AND. ABSGMA.LE.OPTRES)) )
     $   THEN
            BISCOUNT = 0
         ENDIF

         IF( BISCOUNT .NE. 0 )THEN
*
*           -----------
*            Bisection
*           -----------
*
            NBIS = NBIS+1
            LAMBDA   = HALF*(LB + UB)
            BISCOUNT = BISCOUNT - 1

         ELSE
*
*           ---------------
*            Full RQI step
*           ---------------
*
            CALL DLAXRG(
     $             N, TWIST,
     $             RWORK(IXDP), RWORK(IXRP-2 + TWIST+1), E, MINGMA,
     $             CUTTOL, RWORK(IXBUF(BUFSEL)), NORMZ, RESID, RQCORR
     $           )
            NRQI = NRQI+1


*           Note: The current optimum is needed to evaluate the quality
*           of this step. This is the reason why we record the optimum
*           only at the end of the loop body below, and the loop exit is
*           there as well.
            BCKLDA  = LAMBDA
            XLAMBDA = LAMBDA + RQCORR
            STATUS = 0

            IF( RESID .LE. RESTOL )THEN
*
*              converged: residual tolerance is reached
*
               STATUS = 1
C@CHECKS on
C               IF( INDEX .EQ. 150 )THEN
C                  CALL FABORT('','')
C               ENDIF
C@CHECKS !


            ELSEIF( XLAMBDA.EQ.LAMBDA .OR.
     $              ABS(RQCORR).LE.ACCTOL*ABS(LAMBDA) )
     $      THEN
*
*              converged: rqc tolerance reached
*
               STATUS = 2

            ELSEIF( HAVEACC )THEN
*
*              converged: guidance interval has full acc
*                (although this rqc might jump out of it)
*
               STATUS = 3

            ENDIF


*           If STATUS=0 still, we are not yet converged. The remainder of
*           the loop body centers around setting the next iterate LAMBDA.
*           Once we have one, STATUS is set to sth < 0.

            IF( STATUS.EQ.0 .AND. HAVELB .AND. HAVEUB )THEN
*
*              Check if we stay in the guidance interval, if not
*              switch to bisection
*
               IF( XLAMBDA.LE.LB .OR. XLAMBDA.GE.UB )THEN
                  BISCOUNT = PBBLEN
                  STATUS = -1
               ENDIF
            ENDIF

            IF( STATUS.EQ.0 )THEN
*
*              Check if the residual norm is sufficiently better than
*              the current optimum. If so, then continue normally.
*
               IF( .NOT.HAVEOPT .OR. RESID.LE.RESDECFAC*OPTRES )THEN
                  LAMBDA = XLAMBDA
                  STATUS = -2
               ENDIF
            ENDIF

            IF( STATUS.EQ.0 )THEN
*
*              The residual is not good enough, try alternatives to avoid
*              slow convergence.
*
               DIR  = SIGN(ONE,RQCORR)
               ZETA = LAMBDA + DIR*MAX(3*ABS(RQCORR),RESID)

               IF( .NOT.HAVELB .OR. .NOT.HAVEUB )THEN
                  LAMBDA = ZETA
               ELSE
                  IF( DIR.EQ.-ONE )THEN
                     RTMPA = HALF*(LB+UB)
                     RTMPB = UB
                  ELSE
                     RTMPA = LB
                     RTMPB = HALF*(LB+UB)
                  ENDIF

                  IF( RTMPA.LE.ZETA .AND. ZETA.LE.RTMPB )THEN
                     LAMBDA = ZETA
                  ELSE
                     RTMPA = (ONE-ATHIRD)*LB + ATHIRD*UB
                     RTMPB = ATHIRD*LB + (ONE-ATHIRD)*UB
                     IF( RTMPA.LE.XLAMBDA .AND. XLAMBDA.LE.RTMPB )THEN
                        LAMBDA = XLAMBDA
                     ELSE
                        BISCOUNT = PBBLEN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
*
*           Record optimal shift
*
            IF( .NOT.HAVEOPT .OR. RESID.LT.OPTRES )THEN
               OPTLDA  = BCKLDA
               OPTRES  = RESID
               HAVEOPT = .TRUE.
*              Keep this one as optimum, switch buffers
               BUFSEL = 3 - BUFSEL
            ENDIF
*
*           ===========
*            LOOP EXIT
*           ===========
*
            IF( STATUS.GT.0 )  EXIT
*
*           Setup a bisection iterate
*
            IF( BISCOUNT .NE. 0 )THEN
               LAMBDA = HALF*(LB + UB)
            ENDIF
*
         ENDIF
*
      ENDDO
*
*     The loop was only exited after a vector was computed into the active
*     buffer. However, the last computed one need not be the optimum if
*     the iteration stopped due to maximum accuracy being reached. The
*     vector with optimal residual is kept in the backup buffer.
*
      LAMBDA = OPTLDA

      BUFSEL = 3-BUFSEL
      DO IXF77A = 1, N
         Z(IXF77A) = RWORK( IXBUF(BUFSEL) - 1 + IXF77A )
      ENDDO
*
*     Determine the support
*
 100  I = 1
      DO
         IF( Z(I) .NE. ZERO ) EXIT
         I = I+1
      ENDDO
      ISUPPZ(1) = I
      I = N
      DO
         IF( Z(I) .NE. ZERO ) EXIT
         I = I-1
      ENDDO
      ISUPPZ(2) = I

      XNRQI = XNRQI + NRQI
      XNRQIBIS = XNRQIBIS + NBIS
      XMAXNRQI = MAX( XMAXNRQI, NRQI )
      XMAXNRQIBIS = MAX( XMAXNRQIBIS, NBIS )
      END SUBROUTINE DLAXRX
*
************************************************************************
