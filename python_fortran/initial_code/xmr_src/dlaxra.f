      SUBROUTINE DLAXRA(
     $             N, D, E, S, SCALE, ASPTOL, NBLCKS, ABINDS
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N
*
      DOUBLE PRECISION, INTENT(INOUT)  ::  D(N), E(N)
*
      INTEGER,          INTENT(OUT)    ::  NBLCKS
      INTEGER,          INTENT(OUT)    ::  S(N), ABINDS(2*N)
      DOUBLE PRECISION, INTENT(OUT)    ::  ASPTOL, SCALE
*
*
*  Purpose
*  =======
*
*    For a symmetric tridiagonal matrix T of dimension N given by its
*    entries in D and E, apply some preparatory steps:
*
*    (1)  Find a diagonal sign matrix S = diag(+- 1) such that
*         the offdiagonal entries of T := S T S are nonnegative.
*
*    (2)  Scale the matrix entries by SCALE into proper numerical range.
*         NOTE: Here it would be tempting to scale each block
*               individually after step (3), but then DLAXRI cannot
*               determine index ranges anymore (would have to carry
*               the per-block scalings all the way through).
*
*    (3)  Split the matrix into irreducible blocks.
*         To this end, determine a suitable tolerance ASPTOL and set
*         all ofdiagonal entries below that threshold to zero.
*         The resulting block structure is returned in NLBKCS and ABINDS.
*
*    Upon exit, these transformations are reflected in the entries of
*    D and E.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The dimension of the matrix T.
*
*  D       (input/output) DOUBLE PRECISION array, dimension ( N )
*          Upon entry, the diagonal entries of the matrix T.
*          Upon exit, the diagonal entries of the transformed matrix.
*
*  E       (input) DOUBLE PRECISION array, dimension ( N )
*          Upon entry, the offdiagonal entries of the matrix T.
*          Upon exit, the offdiagonal entries of the transformed matrix.
*          These will be zero where a block ends, that is, at
*          ABINDS(2),...,ABINDS(2*NBLCKS), and > ASPTOL elsewhere.
*
*  S       (output)  INTEGER array, dimension ( N )
*          The diagonal scaling matrix S, entries are +-1.
*
*  SCALE   (output)  DOUBLE PRECISION
*          The scaling factor used in step (2).
*
*  ASPTOL  (output)  DOUBLE PRECISION
*          The tolerance used for splitting in step (3). This will be
*          somewhere in the area of n*eps*norm, where norm is the
*          2-norm of the transformed matrix after steps (1) and (2).
*
*  NBLCKS  (output) INTEGER
*          How many irreducible blocks the matrix T was split into,
*          NBLCKS >= 1.
*
*  ABINDS  (output) INTEGER array, dimension ( 2*NBLCKS )
*          Indices where the blocks start and end.
*          For 1 <= i <= NBLCKS, the i'th irreducible block is the
*          principal submatrix ABLCKS(2*i-1):ABLCKS(2*i).
*          Purely for increased ease-of-use there is some redundancy
*          here since always
*            ABINDS(1)=1, ABINDS(2*i)+1=ABINDS(2*i+1),
*          and
*            ABINDS(2*NBLCKS)=N.
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL          DSCAL

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH


*
*     .. Parameters ..
*
      INTEGER, PARAMETER  ::  TOLFAC = 3
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0

*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, OFLOW
      DOUBLE PRECISION  GL, GU, OFF, ALPHA
      DOUBLE PRECISION  RTMP, EMIN, EMAX, TMINNZ, TMAXNZ

      INTEGER           BBEG, BEND, BLEN, IBLCK, I
*
*  ===== Executable Statements ==========================================
*
      EPS   = DLAMCH('Epsilon')
      OFLOW = DLAMCH('Overflow')

*     ----------------------------
*      Build S to make E positive
*     ----------------------------
      S(1) = +1
      DO I = 1, N-1
         IF( S(I)*E(I) .LT. ZERO )THEN
            S(I+1) = -1
         ELSE
            S(I+1) = +1
         ENDIF
         E(I) = ABS(E(I))
      ENDDO

*     -----------------------------------
*      Compute union of Gershgorin Discs
*     -----------------------------------
*
      GL = D(1) - E(1)
      GU = D(1) + E(1)
      DO I = 2, N-1
         OFF = E(I-1) + E(I)
         GL  = MIN( GL, D(I) - OFF )
         GU  = MAX( GU, D(I) + OFF )
      ENDDO
      GL = MIN( GL, D(N) - E(N-1) )
      GU = MAX( GU, D(N) + E(N-1) )
*
*     -----------------------------
*      Set tolerance for splitting
*     -----------------------------
*      Note: Even if every offdiagonal element is <= neps and split away,
*            the norm of the additive error matrix is still bounded by
*            2neps. Hence, a tolerance of c*neps is fine.

      ASPTOL = (TOLFAC * N * EPS) * MAX( ABS(GL), ABS(GU), GU-GL )
*
*     -------
*      Do it
*     -------
*
      E(N) = ZERO
      NBLCKS = 1
      ABINDS(1) = 1
      I = 1
      DO
         IF( E(I) .LE. ASPTOL )THEN
            E(I) = ZERO
            ABINDS(2*NBLCKS) = I
            IF( I .EQ. N )  EXIT
            NBLCKS = NBLCKS + 1
            ABINDS(2*NBLCKS-1) = I+1
         ENDIF
         I = I+1
      ENDDO


*     ------------------
*           Scaling
*     ------------------

*     Determine min/max nonzero entries
      TMINNZ = E(1)
      TMAXNZ = TMINNZ
      DO I = 1, N-1
         TMINNZ = MIN( TMINNZ, E(I) )
         TMAXNZ = MAX( TMAXNZ, E(I) )
      ENDDO
      DO I = 1, N
         IF( D(I) .NE. 0 )THEN
            RTMP = ABS(D(I))
            TMINNZ = MIN( TMINNZ, RTMP )
            TMAXNZ = MAX( TMAXNZ, RTMP )
         ENDIF
      ENDDO

*     Since the matrix is split already, we have the magnitudes of all
*     entries within [ EPS*TMAX, TMAX ], except possibly for smaller
*     entries on the diagonal. There is no reason not to move this
*     window such that it contains one.
      SCALE = ONE
      IF( TMAXNZ .GT. OFLOW/4 )THEN
         SCALE = ONE / TMAXNZ
      ELSEIF( TMINNZ.GT.ONE .OR. TMAXNZ.LT.ONE )THEN
         SCALE = 2 / (TMINNZ + TMAXNZ)
      ENDIF

      IF( SCALE .NE. ONE )THEN
         CALL DSCAL( N,   SCALE, D(1), 1 )
         CALL DSCAL( N-1, SCALE, E(1), 1 )
      ENDIF


      END SUBROUTINE DLAXRA
*
************************************************************************
