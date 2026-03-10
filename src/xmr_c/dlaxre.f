      SUBROUTINE DLAXRE(
     $             N, D, E, GL, GU, ABSERR, GAPTOL,
     $             REPR, REPI, TAU,
     $             EWL_AE, EWL_LU, EMODE,
     $             RWORK,
     $             INFO
     $           )
      IMPLICIT NONE
*
      CHARACTER,        INTENT(IN)     ::  EMODE
      INTEGER,          INTENT(IN)     ::  N
      DOUBLE PRECISION, INTENT(IN)     ::  GL, GU, ABSERR, GAPTOL
      DOUBLE PRECISION, INTENT(IN)     ::  D(N)
*
      DOUBLE PRECISION, INTENT(INOUT)  ::  E( N )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( 6*N )
*
      INTEGER,          INTENT(OUT)    ::  INFO
      INTEGER,          INTENT(OUT)    ::  REPI( (6+N+N/2) )
      INTEGER,          INTENT(OUT)    ::  EWL_AE(1:2*N)
      DOUBLE PRECISION, INTENT(OUT)    ::  TAU
      DOUBLE PRECISION, INTENT(OUT)    ::  REPR( (4*N+3) )
      DOUBLE PRECISION, INTENT(OUT)    ::  EWL_LU(1:2*N)
*
*
*  Purpose
*  =======
*
*    Build root representation and init ew-list.
*
*      Given a symmetric tridiagonal matrix T by its entries in (D,E),
*    we determine a representation REP of a matrix M with
*        (T+err)  -  TAU  =  M,
*    where TAU is a shift and err is a matrix with norm bounded by
*    ABSERR (used for matrices with a nearly constant diagonal, to get
*    them into GK form and then use blocks).
*
*    The specific nature of M is kept hidden, it may be positive or
*    negative definite, or some other kind of good representation (RRR).
*
*    It is assumed that the input tridiagonal T is (numerically)
*    irreducible and properly scaled, that is, one block produced
*    by DLAXRA.
*
*    Upon exit, REPR, REPI and E combined define the representation.
*    (E may be modified, since the representation data does not include
*    the offdiagonal entries itself)
*
*    Once the representation has been found, the ew-list is initialized.
*    How this is done is determined by EMODE:
*     EMODE = 'd' or 'D'    Init ews based on ews that were computed
*                           to full accuracy by dqds.
*                           NOTE: Complexity O(n^2)
*
*     EMODE = 'o' or 'O'    Use Gershgorin Discs for initial setup, in
*                           particular to obtain outer bounds. Maybe
*                           get some more info, but not exceeding linear
*                           complexity.
*                           NOTE: Complexity O(n)
*
*     EMODE = 'g' or 'G'    Use full information from the Gershgorin
*                           Discs.
*                           NOTE: Complexity O(n logn)
*
*     EMODE = 'b' or 'B'    Start with Gershgorin discs (g), then refine
*                           all to full precision using bisection.
*                           NOTE: Complexity O(n^2)
*
* !! At the moment only modes 'd' and 'o' are implemented.
*
*  Arguments
*  =========


*
*  ======================================================================
*
*     .. Declarations ..
*
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
      SUBROUTINE DLAXRE_INITEWLDQDS(
     $             N, REPR, REPI, IL, IU,
     $             QDVALS, WIL, WIU, EWL_AE, EWL_LU
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)   ::  N, IL, IU, WIL, WIU
      INTEGER,          INTENT(IN)   ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)   ::  QDVALS(IL:IU)
      DOUBLE PRECISION, INTENT(IN)   ::  REPR(4*N+3)
*
      INTEGER,          INTENT(OUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*IL-1:2*IU)
*
      END SUBROUTINE DLAXRE_INITEWLDQDS
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRL_UPDATE(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, LAMBDA, XI
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  IL, IU, XI
      DOUBLE PRECISION, INTENT(IN)  ::  LAMBDA
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      END SUBROUTINE DLAXRL_UPDATE
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

      EXTERNAL          DLARNV, DLASQ1


*
*     .. Parameters ..
*



*     How many tries are allowed to find a root representation
*     This is just a stopgap for situations where something is very wrong,
*     normally 2 or at most 3 tries should always suffice.
      INTEGER, PARAMETER  ::  NTRYMAX = 10

*     By how much ULP to perturb the primary data of the root rep.
*     Set to 0 for no perturbation.
      INTEGER, PARAMETER  ::  NULPROOTPERT = 8

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0

*
*     .. Locals ..
*
      DOUBLE PRECISION  SFMIN, EPS, PREC, PIVBASE, EMAX, DMIN, DMAX
      DOUBLE PRECISION  VL, VU, AUXL, AUXU, DL, DU, RTMP, FAC
      DOUBLE PRECISION  VSTART, OFFSET, AUX, DP, SIGND
      DOUBLE PRECISION  GMA
      DOUBLE PRECISION  MGL, MGU, DIAG, EOFF

      INTEGER           JXG, JYOMGA, JXNGN
      INTEGER           I, J, TWIST, TYPE, TRY
      INTEGER           NEGL, NEGU, IXDP, IXRP
      INTEGER           ISEED(4), ITMP
      INTEGER           SQDIM, IXSQD, IXSQE, IXSQW, IINFO

      LOGICAL           GOTREP, DISDEF, RISDEF

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



      EPS   = DLAMCH('Epsilon')
      PREC  = DLAMCH('Precision')
      SFMIN = DLAMCH('Safe Minimum')

      JXG    = (0)
      JYOMGA = (4)
      JXNGN  = (2*N+1)

*     -------------
*      Set PIVBASE
*     -------------
      EMAX = E(1)
      DMIN = D(1)
      DMAX = D(1)
      DO I = 2, N
         DMIN = MIN( DMIN, D(I) )
         DMAX = MAX( DMAX, D(I) )
         EMAX = MAX( EMAX, E(I) )
      ENDDO

      PIVBASE = 4 * SQRT( SFMIN ) * MAX( EMAX, ONE )


      GOTREP = .FALSE.

      IF( .FALSE. )THEN
C        Support for GK-type matrices deactivated.
C        Reason: One should call a specific SVD routine in this case,
*         see below. We ran into problems for GRO-type synthetics,
*         mainly because full blocks in the root restrict possible
*         choices for blocks in the children too much, due to
*         overlap sequences.
C      IF( (DMAX-DMIN) .LE. 2*ABSERR )THEN
*        -----------------------------------
*         Matrices with a constant diagonal
*        -----------------------------------
*        Using an absolute perturbation bounded by ABSERR gives a
*        constant diagonal. We can then shift it away to obtain a
*        nicely blockable matrix in GK-form, possibly with an odd
*        dimension.

*        A bidiagonal matrix can have tiny singular values even
*        if the entries itself are not small, even below PIVBASE.
*        This is a problem since the shifting algorithms (DLAXR[N,T,S])
*        require a shift > PIVBASE/EPS, otherwise they cannot guarantee
*        relative accuracy. Right now we discard the option to interpret
*        the root as GK-matrix if that happens.

         VL = MAX( ABS(GL), ABS(GU), GU-GL ) * (EPS**2)
         NEGL = 0
         AUXL = ZERO
         I = 1
         DO
            DL = ( D(I) - VL ) - AUXL
            IF( ABS(DL).LT.PIVBASE )  DL = -PIVBASE
            IF( DL.LT.ZERO )  NEGL = NEGL+1
            IF( I .EQ. N )THEN
               EXIT
            ENDIF
            RTMP = E(I)**2
            AUXL = RTMP / DL
            I = I+1
         ENDDO

*        Only proceed if there are no ews between 0 and VL
         IF( NEGL .EQ. (N+1)/2 )THEN

            TAU = HALF * (DMIN + DMAX)

            IF( E(1) .LE. E(N-1) )THEN
               ITMP  = 0
               TWIST = N
            ELSE
               ITMP  = 1 - MOD(N,2)
               TWIST = 1
            ENDIF
            DO I = 1, N
               REPR(JXG    + I) = ZERO
               REPI(JYOMGA + I) = ITMP
               ITMP = 1-ITMP
            ENDDO

            TYPE   = (1)
            GOTREP = .TRUE.

         ENDIF

*        Note
*        ====
*          The above support fot this special case is rudimentary at
*          best. One should then really better handle the problem just
*          as a bSVD routine would to it:
*          (1) For tiny singular values, use a couple of sweeps of the
*              implicit zero-shift QR to get them away.
*          (2) Compute at most one half of the eigenpairs for GK, then
*              negate half of the vector entries to get the rest.
*
      ENDIF

      IF( .NOT. GOTREP )THEN
*        --------------------------------------
*         Treat as standard Tridiagonal Matrix
*        --------------------------------------


*        Determine preferred side for shifting
         VL = GL + (GU-GL)/4
         VU = GU - (GU-GL)/4
         NEGL = 0
         NEGU = 0
         AUXL = ZERO
         AUXU = ZERO
         I = 1
         DO
            DL = ( D(I) - VL ) - AUXL
            DU = ( D(I) - VU ) - AUXU
            IF( ABS(DL).LT.PIVBASE )  DL = -PIVBASE
            IF( ABS(DU).LT.PIVBASE )  DU = -PIVBASE
            IF( DL.LT.ZERO )  NEGL = NEGL+1
            IF( DU.LT.ZERO )  NEGU = NEGU+1
            IF( I .EQ. N )THEN
               EXIT
            ENDIF
            RTMP = E(I)**2
            AUXL = RTMP / DL
            AUXU = RTMP / DU
            I = I+1
         ENDDO

         IF( NEGL .GE. N-NEGU )THEN
            VSTART = GL
            SIGND  = +ONE
         ELSE
            VSTART = GU
            SIGND  = -ONE
         ENDIF

         OFFSET = (2 * N * EPS) * MAX( ABS(VSTART), GU - GL )


*        Setup workspace for factorization data
*         RWORK(1:N)    = DP(1:N)
*         RWORK(N+1:2N) = RP(1:N)
         IXDP  = 0
         IXRP  = N

         GOTREP = .FALSE.
         TRY = 0
         DO
            TAU = VSTART - TRY*SIGND*OFFSET

*           Top to Bottom
            DISDEF = .TRUE.
            AUX = ZERO
            I = 1
            DO
               DP = (D(I) - TAU) - AUX
               IF( ABS(DP) .LT. PIVBASE )  DP = SIGND*PIVBASE
               RWORK(IXDP + I) = DP
               IF( I.EQ.N .OR. .NOT.DISDEF )THEN
                  EXIT
               ENDIF
               DISDEF = DISDEF .AND. ( SIGN(ONE,DP) .EQ. SIGND )
               AUX = E(I)**2 / DP
               I = I+1
            ENDDO
            IF( I.EQ.N )THEN
               IF( ABS(RWORK(IXDP + I)) .LE.
     $             MAX(N*EPS*ABS(TAU),ABSERR) )
     $         THEN
                  RWORK(IXDP + I) = ZERO
               ENDIF
            ENDIF

*           Bottom to Top
            RISDEF = .TRUE.
            AUX = ZERO
            I = N
            DO
               DP = (D(I) - TAU) - AUX
               IF( ABS(DP) .LT. PIVBASE )  DP = SIGND*PIVBASE
               RWORK(IXRP + I) = DP
               IF( I.EQ.1 .OR. .NOT.RISDEF )THEN
                  EXIT
               ENDIF
               RISDEF = RISDEF .AND. ( SIGN(ONE,DP) .EQ. SIGND )
               I = I-1
               AUX = E(I)**2 / DP
            ENDDO
            IF( I.EQ.1 )THEN
               IF( ABS(RWORK(IXRP + 1)) .LE.
     $             MAX(N*EPS*ABS(TAU),ABSERR) )
     $         THEN
                  RWORK(IXRP + 1) = ZERO
               ENDIF
            ENDIF



*           We want that both directions are semi-definite, otherwise
*           the results may differ when called with the flipped
*           matrix.
            IF( DISDEF .AND. RISDEF )THEN
*              Take min gamma
               IF( ABS(RWORK(IXDP+1)) .LT. ABS(RWORK(IXRP+N)) )THEN
                  GMA = RWORK(IXRP+1)
                  TWIST = 1
               ELSE
                  GMA = RWORK(IXDP+N)
                  TWIST = N
               ENDIF

*              Set repdata if successful
               DO I = 1, TWIST-1
                  REPR( JXG + I ) = RWORK( IXDP + I )
               ENDDO
               REPR( JXG + TWIST ) = GMA
               DO I = TWIST+1, N
                  REPR( JXG + I ) = RWORK( IXRP + I )
               ENDDO
               DO IXF77A = JYOMGA+1, JYOMGA+N
                  REPI(IXF77A) = 0
               ENDDO
               TYPE = (0)
               GOTREP = .TRUE.
            ENDIF

            TRY = TRY + 1
            IF( TRY.EQ.NTRYMAX .OR. GOTREP )THEN
               EXIT
            ENDIF
         ENDDO

      ENDIF

      IF( .NOT. GOTREP )THEN
         INFO = 1
         RETURN
      ENDIF

*     --------------
*      Perturb Root
*     --------------
*     RWORK is not in use right now, so we can employ it freely.

      IF( NULPROOTPERT .GT. 0 )THEN
         ISEED(1) = 377
         ISEED(2) = 610
         ISEED(3) = 987
         ISEED(4) = 1597
         FAC   = NULPROOTPERT * PREC
         CALL DLARNV( 2, ISEED, 2*N-1, RWORK )
         DO I = 1, N-1
            E(I) = E(I) * (ONE + RWORK(I)*FAC)
         ENDDO
         I = 1
         J = N
         DO
            IF( I .GT. N )  EXIT
            REPR(JXG+I) = REPR(JXG+I) * (ONE + RWORK(J)*FAC)
            I = I+1
            J = J+1
         ENDDO
      ENDIF

      CALL DLAXRR( N, TWIST, TYPE, E, PIVBASE, REPR, REPI )

*     ===================================================================
*     =                        Init EW-List                             =
*     ===================================================================
*     RWORK is not in use
      TWIST = REPI(2)
      TYPE  = REPI(1)

      IF( EMODE.EQ.'d' .OR. EMODE.EQ.'D' )THEN
*        Use DQDS

         IXSQD = 1
         IXSQE = N+1
         IXSQW = IXSQE + N

         IF( TYPE .EQ. (1) )THEN
*           We can ignore if the matrix was flipped, just take the
*           bidiagonal for the GK-form from top to bottom.
            SQDIM = (N+1) / 2
            I = 1
            J = 1
            DO
               RWORK(IXSQD-1 + I) = E(J)
               J = J+1
               IF( I .EQ. SQDIM )  EXIT
               RWORK(IXSQE-1 + I) = E(J)
               J = J+1
               I = I+1
            ENDDO
            CALL DLASQ1(
     $             SQDIM, RWORK(IXSQD), RWORK(IXSQE), RWORK(IXSQW),
     $             IINFO
     $           )
            IF( IINFO .NE. 0 )THEN
               INFO = 2
               RETURN
            ENDIF
*           The singular values of the bidiagonal factor are in SQD,
*           but in descending order. We have to mirror them to obtain
*           eigenvalues of M.
            I = 1
            J = N
            DO
               IF( I .GE. J )  EXIT
               RWORK(IXSQD-1 + J) =   RWORK(IXSQD-1 + I)
               RWORK(IXSQD-1 + I) = - RWORK(IXSQD-1 + I)
               J = J-1
               I = I+1
            ENDDO

         ELSE
            IF( TWIST.EQ.1 )THEN
               SIGND = SIGN(ONE,REPR(JXG+N))
               I = N
               DO
                  RWORK(IXSQD-1 + I) = SQRT(SIGND*REPR(JXG + I))
                  IF( I .EQ. 1 )  EXIT
                  RWORK(IXSQE-1 + I-1) = SQRT(SIGND*REPR(JXNGN + I))
                  I = I-1
               ENDDO
            ELSE
               SIGND = SIGN(ONE,REPR(JXG+1))
               I = 1
               DO
                  RWORK(IXSQD-1 + I) = SQRT(SIGND*REPR(JXG + I))
                  IF( I .EQ. N )  EXIT
                  RWORK(IXSQE-1 + I) = SQRT(SIGND*REPR(JXNGN + I))
                  I = I+1
               ENDDO
            ENDIF

            CALL DLASQ1(
     $             N, RWORK(IXSQD), RWORK(IXSQE), RWORK(IXSQW), IINFO
     $           )
            IF( IINFO .NE. 0 )THEN
               INFO = 2
               RETURN
            ENDIF

*           Now SQD holds the singular values of the bidiagonal factor
*           in descending order. We have to square them, adjust to a
*           negative definite M (SIGND=-1) and invert their order.
            IF( SIGND .EQ. ONE )THEN
*              swap order
               I = 1
               J = N
               DO
                  IF( I .GE. J )  EXIT
                  RTMP = RWORK(IXSQD-1 + I)
                  RWORK(IXSQD-1 + I) = RWORK(IXSQD-1 + J)
                  RWORK(IXSQD-1 + J) = RTMP
                  I = I+1
                  J = J-1
               ENDDO
            ENDIF
*           square and sign
            DO I = 0, N-1
               RWORK(IXSQD + I) = SIGND * RWORK(IXSQD + I)**2
            ENDDO
         ENDIF
*
*
         CALL DLAXRE_INITEWLDQDS(
     $          N, REPR, REPI, 1, N,
     $          RWORK(IXSQD), 1, N, EWL_AE, EWL_LU
     $        )

      ELSEIF( EMODE.EQ.'o' .OR. EMODE.EQ.'O' )THEN
*        Use Gersgorin Discs for outer bounds only.

*        Recompute them for the RRR, instead of shifting GL,GU from
*        the original T.
         EOFF = E(1)
         DIAG = REPR(JXG+1)
         IF( 1 .EQ. TWIST )THEN
            DIAG = DIAG + REPR(JXNGN + 2)
         ENDIF
         MGL = DIAG - EOFF
         MGU = DIAG + EOFF
         DO I = 2, N
            EOFF = E(I-1) + E(I)
            DIAG = REPR(JXG + I)
            IF( I .LE. TWIST )  DIAG = DIAG + REPR(JXNGN + I-1)
            IF( I .GE. TWIST )  DIAG = DIAG + REPR(JXNGN + I+1)
            MGL = MIN( MGL, DIAG - EOFF )
            MGU = MAX( MGU, DIAG + EOFF )
         ENDDO

         DO I = 1, N
            EWL_AE(2*I-1) = 1
            EWL_AE(2*I)   = N
            EWL_LU(2*I-1) = MGL
            EWL_LU(2*I)   = MGU
         ENDDO

      ELSE
         INFO = -12
         RETURN
      ENDIF

*     Always include information from the zero inertia in the ewlist
      VL = ONE
      VU = ONE
      CALL DLAXRL_UPDATE(
     $       1, N, VL, VU, EWL_AE, EWL_LU,
     $       ZERO, DLAXRN(N,REPR,REPI,ZERO)
     $     )
      XNBIS_INIT = XNBIS_INIT + 1

      INFO = 0
      END SUBROUTINE DLAXRE
*
************************************************************************
