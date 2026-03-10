      SUBROUTINE DLAXRF_ENV(
     $            N, E, REPR, REPI, DEPTH, ICBEG, ICEND,
     $            LGAP, UGAP, EWL_AE, EWL_LU, RGINFO,
     $            MODE, ENV, MINENV,
     $            RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, DEPTH, ICBEG, ICEND, MODE
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      INTEGER,          INTENT(IN)  ::  RGINFO(ICBEG-1:ICEND)
      DOUBLE PRECISION, INTENT(IN)  ::  E(N-1), REPR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( 2*N )
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 8*N )
*
      DOUBLE PRECISION, INTENT(OUT)    ::  MINENV
      DOUBLE PRECISION, INTENT(INOUT)  ::  ENV(N)
*
*
*  Purpose
*  =======
*
*     Determine an envelope for the given cluster ICBEG:ICEND.
*
*     Modes:
*      1 - do full env (outer and inner, if possible)
*      2 - do only outer
*      3 - improve outer, meaning ENV is set to outer envelope already,
*          but now inner envs shall improve it
*
*  Notes
*  =====
*
*     We use ideas from Parlett&voemel (submatrix-method); check
*     their paper about envelope localization.
*
*     For even trying to do an envelope we require that a (sub-)cluster
*     I:J fulfills
*     (1)  |J-I+1| < N/3
*     (2)  width < mingap / ENVGAPFAC
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
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
      END SUBROUTINE DLAXRF_GRPENV
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO  = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF  = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE   = 1.0D0

      INTEGER, PARAMETER  ::  MODE_FULL = 1
      INTEGER, PARAMETER  ::  MODE_ONLYOUTER = 2
      INTEGER, PARAMETER  ::  MODE_IMPROVEOUTER = 3

      INTEGER, PARAMETER  ::  GI_NOFULL  = -2
      INTEGER, PARAMETER  ::  GI_UNKNOWN = -1
      INTEGER, PARAMETER  ::  GI_NOGAP   =  0
      INTEGER, PARAMETER  ::  GI_INTGAP  =  1
      INTEGER, PARAMETER  ::  GI_FULLGAP =  2

*
*     .. Parameters ..
*


*         Cluster boundaries will be refined to relative accuracy
*         RACLBFAC*N*EPS
      INTEGER, PARAMETER  ::  RACLBFAC = 4

      INTEGER, PARAMETER  ::  MAXRELCOND = 10

      INTEGER, PARAMETER  ::  ENVGAPFAC = 100
      LOGICAL, PARAMETER  ::  USEENV = .TRUE.


*
*     .. Locals ..
*
      DOUBLE PRECISION  LB, UB, GAPL, GAPU, MIENV
      INTEGER           I, J, K, NGRPS, IXIENV, IXGENV, IXWORK
      INTEGER IXF77A, IXF77B, IXF77C
      LOGICAL           GRPSOK, TOOMANY, TOOWIDE

*
*  ===== Executable Statements ==========================================
*

      IF( .NOT.USEENV )THEN
         DO IXF77A = 1, N
            ENV(IXF77A) = ONE
         ENDDO
         MINENV   = ONE
         RETURN
      ENDIF

*
*     Real Workspace
*
      IXIENV = 1
      IXGENV = N+1
      IXWORK = 2*N+1
*       --------
*             2N
*           + 6N for GRPENV

*     Integer Workspace
*             2N for GRPENV


*     ------------------------
*      Compute outer envelope
*     ------------------------
      IF( MODE.EQ.MODE_FULL .OR. MODE.EQ.MODE_ONLYOUTER )THEN
         LB = EWL_LU(2*ICBEG-1)
         UB = EWL_LU(2*ICEND)
         TOOMANY = ( (3*(ICEND-ICBEG+1)) .GT. N )
         TOOWIDE = ( ((UB-LB)*ENVGAPFAC).GT.MIN(LGAP,UGAP) )
         IF( TOOMANY .OR. TOOWIDE )THEN
            DO IXF77A = 1, N
               ENV(IXF77A) = ONE
            ENDDO
            MINENV   = ONE
         ELSE
            CALL DLAXRF_GRPENV(
     $             N, E, REPR, REPI, ICBEG, ICEND,
     $             LGAP, LB, UB, UGAP, ENV,
     $             RWORK(IXWORK), IWORK
     $          )
            MINENV = ONE
            DO K = 1, N
               MINENV = MIN( MINENV, ENV(K) )
            ENDDO
         ENDIF
      ENDIF


      IF( MODE .EQ. MODE_ONLYOUTER )THEN
         RETURN
      ENDIF

*     -----------------------------------
*      Compute envelope from subclusters
*     -----------------------------------

*     First pass
*     Determine number of subgroups and check that all satisfy the
*     necessary criteria to make an envelope useful
      NGRPS = 0
      GRPSOK = .TRUE.
      I = ICBEG
      DO
         J = EWL_AE(2*I)
         DO
            IF( RGINFO(J) .GT. 0 )  EXIT
            J = EWL_AE(2*(J+1))
         ENDDO
         NGRPS = NGRPS + 1
*
         LB = EWL_LU(2*I-1)
         UB = EWL_LU(2*J)
         IF( I .EQ. ICBEG )THEN
            GAPL = LGAP
         ELSE
            GAPL = LB - EWL_LU(2*I-2)
         ENDIF
         IF( J .EQ. ICEND )THEN
            GAPU = UGAP
         ELSE
            GAPU = EWL_LU(2*J+1) - UB
         ENDIF
         TOOMANY = ( (3*(J-I+1)) .GT. N )
         TOOWIDE = ( ((UB-LB)*ENVGAPFAC).GT.MIN(GAPL,GAPU) )
         GRPSOK = GRPSOK .AND. .NOT.(TOOMANY.OR.TOOWIDE)
*
         IF( J.EQ.ICEND .OR. .NOT.GRPSOK )  EXIT
         I = J+1
      ENDDO

      IF( .NOT. GRPSOK )THEN
         RETURN
      ENDIF
      IF( NGRPS .EQ. 1 )THEN
         RETURN
      ENDIF

*     Second pass
*     Accumulate envelopes for the subgroups
      DO IXF77A = IXIENV, IXIENV+N-1
         RWORK(IXF77A) = ZERO
      ENDDO
      MIENV = ZERO
*
      I = ICBEG
      DO
         J = EWL_AE(2*I)
         DO
            IF( RGINFO(J) .GT. 0 )  EXIT
            J = EWL_AE(2*(J+1))
         ENDDO

         LB = EWL_LU(2*I-1)
         UB = EWL_LU(2*J)
         IF( I .EQ. ICBEG )THEN
            GAPL = LGAP
         ELSE
            GAPL = LB - EWL_LU(2*I-2)
         ENDIF
         IF( J .EQ. ICEND )THEN
            GAPU = UGAP
         ELSE
            GAPU = EWL_LU(2*J+1) - UB
         ENDIF

         CALL DLAXRF_GRPENV(
     $          N, E, REPR, REPI, I, J,
     $          GAPL, LB, UB, GAPU,
     $          RWORK(IXGENV),
     $          RWORK(IXWORK), IWORK
     $        )
*
*        accumulate
*
         MIENV = ZERO
         DO K = 0, N-1
            RWORK(IXIENV + K) = RWORK(IXIENV + K) + RWORK(IXGENV + K)**2
            MIENV = MIN( MIENV, RWORK(IXIENV + K) )
         ENDDO

         IF( MIENV.GE.ONE .OR. J.EQ.ICEND )  EXIT
         I = J+1
      ENDDO

      IF( MIENV .GE. ONE )THEN
         RETURN
      ENDIF

*     at this point, the inner envelope holds the squares of the
*     entries of a valid envelope approximation
      MINENV = ONE
      DO K = 1, N
         ENV(K) = MIN( ENV(K), SQRT(RWORK(IXIENV-1 + K)) )
         MINENV = MIN( MINENV, ENV(K) )
      ENDDO




      END SUBROUTINE DLAXRF_ENV
*
************************************************************************
