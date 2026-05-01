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
*
*  Purpose
*  =======
*
*     Selects a twist with optimal element growth and relative condition
*     number. For this twist, the characteristics of the resulting rep
*     are determined and returned in ELG, RCOND, NCD and combined as
*     EVAL.
*
*     If all entries in TWISTOK are 0, the routine sets K=-1 and returns.
*
*  Arguments
*  =========
*
*  KMINUS  (input)  INTEGER
*          Twist in the source.
*
*  TWISTOK (input)  INTEGER array, dimension (N)
*          Flags where a twist is allowed, as delivered by DLAXRS.
*
*  TAU     (input)  DOUBLE PRECISION
*  LOC     (input)  INTEGER, in {1,N}
*  DIR     (input)  INTEGER, in {-1,+1}
*          Specify the eigenvalue we are close to, on which side, and
*          the shift itself.
*
*  WINOK   (input)  LOGICAL
*          Set to true if the resulting shift shall have a consistent
*          inertia wrt LOC and DIR. The routine will try to get a
*          twist with consistent inertia anyway, but the impact of
*          WINOK is that it causes only twists with a consistent
*          inertia to be considered in the first place.
*
*  GAMMA   (input)  DOUBLE PRECISION array, dimension (N)
*          In concert with looking for consistent inertias, the routine
*          may set very small gammas (wrt |tau|) to zero.
*
*  K       (output)  INTEGER
*          The selected twist.
*          Set to <= 0 if no suitable twist was found:
*           0 if no twist had a consistent inertia when WINOK=.TRUE.
*          -1 if no twist was possible
*
*  EVAL    (output)  DOUBLE PRECISION, >= 0
*          Evaluation, combined from the normalised single features.
*          Is <= 1 iff the rep passes the a priori criteria.
*          Only set if K!=-1.
*
*          The detailed features of the representation that constitute
*          eval are returned in the first entries of RWORK.
*          At the moment these are:
*            RWORK(1) = element growth without envelope / spdiam
*            RWORK(2) = element growth with envelope / spdiam
*            RWORK(3) = -1 [reserved for relcond without envelope]
*            RWORK(4) = relcond with envelope
*            RWORK(5) = -1 [reserved for ncd]
*            RWORK(6) = -1 [reserved for ncd]
*
*  ======================================================================
*
*     .. Constants ..
*
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
*
*     .. Parameters ..
*

      INTEGER, PARAMETER  ::  MAXRELCOND = 10

      INTEGER, PARAMETER  ::  MAXGROWTH = 8

*     Restrict possible choice of twist for node reps:
*      0 - no restriction
*      1 - always twist at 1
*      2 - always twist at n
      INTEGER, PARAMETER  ::  TWISTRESTR = 0

*
*     .. Declarations ..
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

      INTERFACE
      SUBROUTINE DLAXRF_SELTW_PART(
     $             N, IA, IE, DIR, G, GN, OMEGA, TWISTOK, S,
     $             ARCSQ, ABETA, AGMAX, AGSMAX, ABRCFAC
     $           )
      IMPLICIT NONE
*
      INTEGER, INTENT(IN)  ::  N, IA, IE, DIR
      INTEGER, INTENT(IN)  ::
     $   TWISTOK( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   OMEGA(   MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
      DOUBLE PRECISION, INTENT(IN) ::
     $   G(  MIN(IA,IE) : MAX(IA,IE) ),
     $   GN( MIN(IA,IE) : MAX(IA,IE) ),
     $   S( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
*
      DOUBLE PRECISION, INTENT(INOUT)  ::
     $   ARCSQ ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   ABETA ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   AGMAX ( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   AGSMAX( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) ),
     $   ABRCFAC( MIN(IA,IE+DIR) : MAX(IA,IE+DIR) )
      END SUBROUTINE DLAXRF_SELTW_PART
      END INTERFACE

*
*     .. Local Variables ..
*
      DOUBLE PRECISION  EPS, GMATHRESH
      DOUBLE PRECISION  BELG, BELG_NOENV, BRCOND, DELG, DRCOND
      DOUBLE PRECISION  FEATS(6), ACTEVAL
      DOUBLE PRECISION  GROWTHBOUND, INVGRBND, INVMAXRC, FAILXRRR
      INTEGER           IXRCSQ, IXBETA, IXGMAX, IXGSMAX, IXMXBRC, IYINER
      INTEGER           I, FOK, LOK, XI, NINOK, TWINOK
      LOGICAL           ACTINEROK, INEROK
*
*  ===== Executable Statements ==========================================
*
*     Determine range of possible twists
*
      FOK = 1
      DO
         IF( TWISTOK(FOK).NE.0 ) EXIT
         FOK = FOK+1
         IF( FOK.EQ.N+1 ) EXIT
      ENDDO
      LOK = N
      DO
         IF( TWISTOK(LOK).NE.0 ) EXIT
         LOK = LOK-1
         IF( LOK.EQ.0 ) EXIT
      ENDDO

*
*     Apply twist restriction
*
      IF( TWISTRESTR .EQ. 1 )  LOK = 1
      IF( TWISTRESTR .EQ. 2 )  FOK = N

*
*     Return if no twist possible
*
      IF( FOK .GT. LOK )THEN
         K = -1
         RETURN
      ENDIF

*     ----------------------------
*      Compute per-twist Inertias
*     ----------------------------
      IYINER = 1
      IWORK( IYINER+0 : IYINER+N+1 ) = 0

      EPS = DLAMCH('Epsilon')
      GMATHRESH = EPS*ABS(TAU)

*     We compute in two sweeps the inertias for the parts above
*     and below the twist candidate in IWORK(IYINER + [FOK:LOK]).
*     Entries where a block ends are not set, except if the block
*     is at the end.
      XI = 0
      I = 1
      DO
         IWORK(IYINER + I) = XI
         IF( I .GE. LOK ) EXIT
         IF( OMEGAD(I+1) .NE. 0 )THEN
            XI = XI + 2
            I = I+1
         ELSEIF( D(I) .LT. ZERO )THEN
            XI = XI + 2
         ENDIF
         I = I+1
      ENDDO
      IF( LOK.EQ.N .AND. OMEGAD(N).NE.0 )THEN
         IWORK(IYINER + N) = XI
      ENDIF
*
      XI = 0
      I = N
      DO
         IWORK(IYINER + I) = IWORK(IYINER + I) + XI
         IF( I .LE. FOK ) EXIT
         IF( OMEGAR(I-1) .NE. 0 )THEN
            XI = XI + 2
            I = I-1
         ELSEIF( R(I) .LT. ZERO )THEN
            XI = XI + 2
         ENDIF
         I = I-1
      ENDDO
      IF( FOK.EQ.1 .AND. OMEGAR(1).NE.0 )THEN
         IWORK(IYINER + 1) = XI
      ENDIF
*
*     Note: The two special case handlers after the loops above already
*     set correct inertias for blocks at end, so we can ignore indices
*     where blocks end now.
      TWINOK = -1
      NINOK = 0
      DO I = FOK, LOK
         IF( OMEGAD(I).EQ.0 .AND. OMEGAR(I).EQ.0 )THEN
            IF( ABS(GAMMA(I)) .LE. GMATHRESH )THEN
*              Set tiny gammas to zero (dlaxrs may or may not have done
*              that already).
*              Note that the inertia does not yet incorporate gamma.
               IF( DIR.EQ.0 .OR. IWORK(IYINER + I).EQ.2*LOC-2 )THEN
*                 We will land exactly on ew LOC
                  GAMMA(I) = ZERO
               ELSE
*                 May already be to far inside, try to rescue what
*                 we can
                  GAMMA(I) = -DIR * GMATHRESH
               ENDIF
            ENDIF
            IF( GAMMA(I) .LT. ZERO )THEN
               IWORK(IYINER + I) = IWORK(IYINER + I) + 2
            ELSEIF( GAMMA(I) .EQ. ZERO )THEN
               IWORK(IYINER + I) = IWORK(IYINER + I) + 1
            ENDIF
         ENDIF
         IF( TWISTOK(I) .NE. 0 )THEN
            XI = IWORK(IYINER + I)
            IF( (DIR.EQ.-1 .AND. XI.LE.(2*LOC-1)) .OR.
     $          (DIR.EQ.+1 .AND. XI.GE.(2*LOC-1)) )
     $      THEN
               NINOK = NINOK + 1
               IF( TWINOK.EQ.-1 .OR.
     $             ABS(I-KMINUS).LT.ABS(TWINOK-KMINUS) )
     $         THEN
                  TWINOK = I
               ENDIF
            ENDIF
         ENDIF
      ENDDO
*     Now, if TWINOK != -1 then it is set to a twist with consistent
*     inertia that is as close as possible to the source twist.

      IF( WINOK .AND. NINOK.EQ.0 )THEN
*        None of the twists give a consistent inertia
         K = 0
         RETURN
      ENDIF

*     ------------------------------------
*      Special case: Take guaranteed RRRs
*     ------------------------------------
      IF( TWINOK.NE.-1 .AND.
     $    ( (LOC.EQ.1 .AND. DIR.EQ.-1) .OR.
     $      (LOC.EQ.N .AND. DIR.EQ.+1) )
     $  )THEN
*        Rationale for eval: If definite then entries are all bounded
*        by the norm.
         K = TWINOK
         EVAL = ONE / MIN( MAXGROWTH, MAXRELCOND )
         XIZERO = IWORK(IYINER + TWINOK)
         FEATS(1) = ONE / SPDIAM
         FEATS(2) = FEATS(1)
         FEATS(4) = ONE
         RETURN
      ENDIF

*
*     ----------
*      Evaluate
*     ----------
*
      IXRCSQ  = 0
      IXBETA  = N
      IXGMAX  = 2*N
      IXGSMAX = 3*N
      IXMXBRC = 4*N
*
      RWORK(1:5*N) = ZERO
*
*     ---------------
*      Top-to-Bottom
*     ---------------
*
      IF( LOK .GT. 1 )THEN
         CALL DLAXRF_SELTW_PART(
     $          N, 1, LOK-1, +1, D, E, OMEGAD, TWISTOK, ENV,
     $          RWORK(IXRCSQ+1), RWORK(IXBETA+1),
     $          RWORK(IXGMAX+1), RWORK(IXGSMAX+1),
     $          RWORK(IXMXBRC+1)
     $        )
      ENDIF
*
*     ---------------
*      Bottom-to-Top
*     ---------------
*
      IF( FOK .LT. N )THEN
         CALL DLAXRF_SELTW_PART(
     $          N, N, FOK+1, -1,
     $          R(FOK+1), E(FOK), OMEGAR(FOK), TWISTOK(FOK), ENV(FOK),
     $          RWORK(IXRCSQ + FOK), RWORK(IXBETA + FOK),
     $          RWORK(IXGMAX + FOK), RWORK(IXGSMAX + FOK),
     $          RWORK(IXMXBRC + FOK)
     $        )
      ENDIF
*
*     ---------
*      Combine
*     ---------
*
      GROWTHBOUND = MAXGROWTH * SPDIAM
      INVGRBND = ONE / GROWTHBOUND
      INVMAXRC = ONE / DBLE( MAXRELCOND )
      FAILXRRR = ((N-1)*MINGAP) / (SPDIAM*SQRT(EPS))

      FEATS = -ONE
      K = -1
      EVAL = -ONE
      XIZERO = -1
      INEROK = .FALSE.
      DO I = FOK, LOK
         IF( TWISTOK(I) .NE. 0 )THEN
*           Compute rc and elg that would result from choosing this
*           twist.
            BRCOND = RWORK(IXRCSQ + I)
     $               + ABS( ENV(I)*(RWORK(IXBETA + I) + ENV(I)) )

            BRCOND = MAX( BRCOND * RWORK(IXMXBRC + I), ONE )
            DRCOND = BRCOND * INVMAXRC

            BELG_NOENV = MAX( RWORK(IXGMAX + I), ABS( GAMMA(I) ) )
            BELG = MAX( RWORK(IXGSMAX + I), ABS( ENV(I)*GAMMA(I) ) )

C            IF( GROWTHBOUND.LT.BELG .AND. BELG.LT.FAILXRRR )
C     $      THEN
               DELG = BELG
C            ELSE
C               DELG = BELG_NOENV
C            ENDIF
            DELG = DELG * INVGRBND

            ACTEVAL = DRCOND + DELG

            XI = IWORK(IYINER + I)
            ACTINEROK = ( (DIR.EQ.-1 .AND. XI.LE.(2*LOC-1)) .OR.
     $                    (DIR.EQ.+1 .AND. XI.GE.(2*LOC-1)) )

            IF( K.EQ.-1 .OR.
     $          (ACTEVAL.LT.EVAL .AND.
     $           (.NOT.INEROK .OR. ACTINEROK .OR. 2*ACTEVAL.LT.EVAL)
     $        ) )
     $      THEN
               K = I
               EVAL = ACTEVAL
               XIZERO = XI
               INEROK = ACTINEROK
               FEATS(1) = BELG_NOENV / SPDIAM
               FEATS(2) = BELG / SPDIAM
               FEATS(4) = BRCOND
            ENDIF
         ENDIF
      ENDDO

      IF( K .NE. -1 )THEN
         RWORK(1:6) = FEATS(1:6)
      ENDIF
      END SUBROUTINE DLAXRf_SELTW
*
************************************************************************
