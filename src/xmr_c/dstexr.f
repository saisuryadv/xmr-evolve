      SUBROUTINE DSTEXR(
     $             N, D, E, WIL, WIU,
     $             W, Z, LDZ, ISUPPZ,
     $             RWORK, LRWORK, IWORK, LIWORK,
     $             INFO
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, WIL, WIU, LDZ
      DOUBLE PRECISION, INTENT(IN)     ::  D(N), E(N-1)
*
      INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK
      INTEGER,          INTENT(INOUT)  ::  IWORK( * )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( * )
*
      INTEGER,          INTENT(OUT)    ::  INFO
      INTEGER,          INTENT(OUT)    ::  ISUPPZ( 2*WIL-1 : 2*WIU )
      DOUBLE PRECISION, INTENT(OUT)    ::  W( WIL : WIU )
      DOUBLE PRECISION, INTENT(OUT)    ::  Z( LDZ, WIL : WIU )
*
*
*  Purpose
*  =======
*
*    DSTEXR computes selected eigenpairs with indices wil:wiu for a symmetric
*    tridiagonal matrix T. The indexing of eigenpairs is with respect to an
*    ascending order of the eigenvalues, and the delivered eigenvalues will be
*    ordered ascendingly.
*
*    The computed eigenpairs will be consistent in the following sense: If the
*    routine is called twice with index sets wil1:wiu1 and wil2:wiu2, and those
*    index sets are non-overlapping, say wiu1 < wil2, then the results will obey
*    (1)  The computed eigenvalues from the first call are all <= the computed
*         eigenvalues from the second call.
*    (2)  The computed vectors from both calls will be numerically orthogonal
*          to each other.
*
*    Esp. (2) is the key to an easy 'naive' parallelization approach, where
*    each processor computes one about equally sized chunk of the desired
*    eigenpairs separately from all other processors.
*    The benefit of (1) is that the computed eigenpairs won't need to be
*    sorted across processors afterwards, that is, the whole computation can
*    be performed without *any* communication.
*
*  Arguments
*  =========
*
*  D       (input) DOUBLE PRECISION ARRAY, dimension (N).
*          The diagonal entries of T.
*
*  E       (input) DOUBLE PRECISION ARRAY, dimension (N-1).
*          The offdiagonal entries of T.
*
*  WIL     (input) INTEGER
*  WIU     (input) INTEGER
*          The range of wanted eigenpairs, 1 <= WIL <= WIU <= N.
*
*  W       (output) DOUBLE PRECISION array, dimension ( WIL:WIU )
*          The computed eigenvalues WIL:WIU.
*          The i'th entry is only set if the eigenpair could be computed,
*          as indicated by ISUPPZ.
*
*  Z       (output) DOUBLE PRECISION array, dimension ( LDZ, WIL:WIU )
*          The computed orthonormal eigenvectors WIL:WIU.
*          The i'th column is only set if the eigenpair could be
*          computed, as indicated by ISUPPZ.
*
*  LDZ     (input) INTEGER
*          Leading dimension of array Z.
*
*  ISUPPZ  (output) INTEGER array, dimension (2*WIL-1:2*WIU)
*          If the i'th vector was succsessfully computed (cf INFO), where
*          WIL <= i <= WIU, then the vector Z(i) has nonzero entries
*          only within ISUPPZ(2*i-1):ISUPPZ(2*i), which is a subset
*          of 1:N.
*          Failure to compute the vector is indicated by
*          ISUPPZ(2*i-1) = 0, ISUPPZ(2*i) < 0, the precise value of
*          the latter may give information about the kind of failure:
*          -1:  depth limit exceeded
*          -2:  no rep found in DLAXRF
*          -3:  could not find a root rep
*
*  RWORK   (workspace) DOUBLE PRECISION, dimension (LRWORK)
*
*  LRWORK  (input/output) INTEGER
*          Dimension of the array RWORK. Set to -1 for a workspace
*          query, then the routine only computes the required real
*          workspace, sets LRWORK to this value and returns.
*            For a workspace query only N needs to be set, no other
*          argument is referenced. You can do a real and integer ws
*          query at once.
*
*  IWORK   (workspace) INTEGER, dimension (LIWORK)
*
*  LIWORK  (input/output) INTEGER
*          Dimension of the array IWORK. Set to -1 for a workspace
*          query, then the routine only computes the required integer
*          workspace, sets LIWORK to this value and returns.
*            For a workspace query only N needs to be set, no other
*          argument is referenced. You can do a real and integer ws
*          query at once.
*
*  INFO    (output) INTEGER
*          Indicate succes of the computation. For INFO=0 or INFO=1, the
*          result data structures are set, for other values this need not
*          be the case.
*          < 0 : Some argument had an illegal value
*            0 : Everything ok, all desired eigenpairs were computed.
*            1 : At least some eigenpairs could not be computed, refer
*                to ISUPPZ to see which.
*          > 1 : Some other error.
*
*  ======================================================================
Cextract -b parameters.inc.f procname=stexr
*
*     .. Parameters ..
*
*     For blocks of size <= QRDIM, the QR-Algorithm is invoked.
      INTEGER,          PARAMETER  ::  QRDIM = 0

      DOUBLE PRECISION, PARAMETER  ::  GAPTOL = 0.001D0

*     Activate the use of dqds to compute root-eigenvalues. Regardless
*     of this setting DQDS is only empleyed if all ews belonging to
*     a block are wanted.
      LOGICAL,          PARAMETER  ::  USEDQD = .TRUE.
*
*     .. Declarations ..
*
      INTERFACE
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
      END SUBROUTINE DLAXRA
      END INTERFACE
      INTERFACE
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
      END SUBROUTINE DLAXRE
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRI(
     $             N, D, E, WIL, WIU, NBLCKS, ABINDS,
     $             ABGERS, ABWIND, ABVSEP,
     $             RWORK, IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, WIL, WIU, NBLCKS
      INTEGER,          INTENT(IN)     ::  ABINDS(2*NBLCKS)
      DOUBLE PRECISION, INTENT(IN)     ::  D(N), E(N)
*
      INTEGER,          INTENT(INOUT)  ::  IWORK( 5 * NBLCKS )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( N )
*
      INTEGER,          INTENT(OUT)    ::  ABWIND(2*NBLCKS)
      DOUBLE PRECISION, INTENT(OUT)    ::  ABGERS(2*NBLCKS)
      DOUBLE PRECISION, INTENT(OUT)    ::  ABVSEP(2*NBLCKS)
*
      END SUBROUTINE DLAXRI
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRO( N, M, W, Z, LDZ, ISUPPZ, REVORD )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)     ::  N, M, LDZ
*
      INTEGER,          INTENT(INOUT)  ::  ISUPPZ(2*M)
      DOUBLE PRECISION, INTENT(INOUT)  ::  W(M), Z(LDZ,M)
*
      INTEGER,          INTENT(OUT)    ::  REVORD(M)
*
      END SUBROUTINE DLAXRO
      END INTERFACE
      INTERFACE
      SUBROUTINE WSREQ_XRV(N, REQR, REQI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)   ::  N
      INTEGER, INTENT(OUT)  ::  REQR, REQI
      END SUBROUTINE WSREQ_XRV
*
************************************************************************
*
      SUBROUTINE DLAXRV(
     $             N, E, ROOTR, ROOTI, EWL_AE, EWL_LU,
     $             WIL, WIU, SPDIAM, GAPTOL,
     $             W, Z, LDZ, ISUPPZ,
     $             RWORK, LRWORK, IWORK, LIWORK, INFO
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, WIL, WIU, LDZ
      INTEGER,          INTENT(IN)  ::  ROOTI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  SPDIAM, GAPTOL
      DOUBLE PRECISION, INTENT(IN)  ::  E( N-1 ), ROOTR(4*N+3)
*
      INTEGER,          INTENT(INOUT)  ::  LRWORK, LIWORK
      INTEGER,          INTENT(INOUT)  ::  EWL_AE( 1:2*N )
      INTEGER,          INTENT(INOUT)  ::  IWORK( LIWORK )
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU( 1:2*N )
      DOUBLE PRECISION, INTENT(INOUT)  ::  RWORK( LRWORK )
*
      INTEGER,          INTENT(OUT)  ::  INFO
      INTEGER,          INTENT(OUT)  ::  ISUPPZ( 2*WIL-1 : 2*WIU )
      DOUBLE PRECISION, INTENT(OUT)  ::  Z( LDZ, WIL : WIU )
      DOUBLE PRECISION, INTENT(OUT)  ::  W( WIL : WIU )
      END SUBROUTINE DLAXRV
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH

      EXTERNAL DSECND
      DOUBLE PRECISION DSECND

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  ONE  = 1.0D0


*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, ASPTOL
      DOUBLE PRECISION  BGERL, BGERU, BSDIAM, BSEPL, BSEPU
      DOUBLE PRECISION  SCALE, INVSC, BSHIFT

      INTEGER           WSREQR, WSREQI
      INTEGER           NBLCKS, IB, BBEG, BEND, BLEN
      INTEGER           IYWORK, IYS, IYBLCK, IYWIND, IYREP, IYEWAE
      INTEGER           IXWORK, IXD, IXE, IXGERS, IXVSEP
      INTEGER           IXREP, IXEWLU, WSI_XRV, WSR_XRV
      INTEGER           BWIL, BWIU, ZA, ZE, LRWREM, LIWREM
      INTEGER           I, J, M, IINFO, NCOMP

      CHARACTER         EMODE

      LOGICAL           GOTALL


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
      DOUBLE PRECISION  TSTAMP0, TSTAMP1
*
*  ===== Executable Statements ==========================================
*
*     ===================================================================
*     =                     Check input parameters                      =
*     ===================================================================


*     ----------------------
*      Workspace Management
*     ----------------------

      IXWORK = 1
      IYWORK = 1

*     Copies of the matrix data to modify
      IXD    = IXWORK
      IXE    = IXWORK + N
      IXWORK = IXWORK + 2*N

*     Outer scaling by DLAXRA to get offdiagonals positive
      IYS    = IYWORK
      IYWORK = IYWORK + N

*     Begin and end of blocks in the splitted matrix
*     ( we assume here N as upper bound for number of blocks, however
*       then we would not call the kernel, so this could be optimized )
      IYBLCK = IYWORK
      IYWORK = IYWORK + 2*N

*     Index ranges of desired eigenvalues for each block
      IYWIND = IYWORK
      IYWORK = IYWORK + 2*N

*     Union of Gershgorin Discs and Value Separators for each block
      IXGERS = IXWORK
      IXVSEP = IXWORK + 2*N
      IXWORK = IXWORK + 4*N

*     EW-List for the current block
      IYEWAE = IYWORK
      IXEWLU = IXWORK
      IYWORK = IYWORK + 2*N
      IXWORK = IXWORK + 2*N

*     Root Representation for the current block
      IYREP  = IYWORK
      IXREP  = IXWORK
      IYWORK = IYWORK + (6+N+N/2)
      IXWORK = IXWORK + (4*N+3)

      LRWREM = LRWORK - IXWORK + 1
      LIWREM = LIWORK - IYWORK + 1
*
*     --------------------------
*      Handle workspace queries
*     --------------------------
*
      CALL WSREQ_XRV( N, WSR_XRV, WSI_XRV )

      WSREQR = IXWORK-1 + WSR_XRV
      WSREQI = IYWORK-1 + WSI_XRV

      IF( LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )THEN
         IF( LRWORK.EQ.-1 )  LRWORK = WSREQR
         IF( LIWORK.EQ.-1 )  LIWORK = WSREQI
         RETURN
      ENDIF

      TSTAMP0 = DSECND()

      EPS = DLAMCH('Epsilon')

*     Copy D and E so that we can modify them
      DO I = 1, N
         RWORK(IXD + I - 1) = D(I)
      ENDDO
      DO I = 1, N-1
         RWORK(IXE + I - 1) = E(I)
      ENDDO

*     ===================================================================
*     =                          Splitting                              =
*     ===================================================================
*
*     Split the matrix into irreducible blocks, record where they begin
*     and end, and scale them individually into proper numerical range.
*

      CALL DLAXRA(
     $       N, RWORK(IXD), RWORK(IXE),
     $       IWORK(IYS), SCALE, ASPTOL, NBLCKS, IWORK(IYBLCK)
     $     )
      INVSC = ONE / SCALE

*     ========================
*      Determine Index ranges
*     ========================

      CALL DLAXRI(
     $       N, RWORK(IXD), RWORK(IXE), WIL, WIU, NBLCKS, IWORK(IYBLCK),
     $       RWORK(IXGERS), IWORK(IYWIND), RWORK(IXVSEP),
     $       RWORK(IXWORK), IWORK(IYWORK)
     $     )

*     ===================================================================
*     =                       Build the Vectors                         =
*     ===================================================================
      TSTAMP1 = DSECND()
      XTIME1  = XTIME1 + (TSTAMP1 - TSTAMP0)
      TSTAMP0 = TSTAMP1
      XNBLCKS = XNBLCKS + NBLCKS

      DO J = 1, WIU - WIL + 1
         DO I = 1, LDZ
            Z(I,J) = ZERO
         ENDDO
      ENDDO

      DO I = 1, WIU - WIL + 1
         W(I) = ZERO
      ENDDO

      M = WIU - WIL + 1
      GOTALL = .TRUE.
      ZE = WIL-1
      DO IB = 1, NBLCKS
         BBEG   = IWORK(IYBLCK-1 + 2*IB-1)
         BEND   = IWORK(IYBLCK-1 + 2*IB)
         BLEN   = BEND - BBEG + 1
         BWIL   = IWORK(IYWIND-1 + 2*IB-1)
         BWIU   = IWORK(IYWIND-1 + 2*IB)
         BGERL  = RWORK(IXGERS-1 + 2*IB-1)
         BGERU  = RWORK(IXGERS-1 + 2*IB)
         BSDIAM = BGERU - BGERL
         BSEPL  = RWORK(IXVSEP-1 + 2*IB-1)
         BSEPU  = RWORK(IXVSEP-1 + 2*IB)

         IF( BWIL .GT. BWIU )THEN
            CYCLE
         ENDIF

         ZA = ZE + 1
         ZE = ZA + (BWIU - BWIL)

         IF( BLEN .EQ. 1 )THEN

            Z( BBEG, ZA )  = ONE
            W( ZA )        = RWORK(IXD-1 + BBEG)
            ISUPPZ(2*ZA-1) = 1
            ISUPPZ(2*ZA)   = 1

C These special handlers are not implemented yet, since esp for debugging
C and testing we want to call the kernel on small matrices.
C         ELSEIF( BLEN .EQ. 2 )THEN
C
C            call special 2x2 handler
C
C         ELSEIF( BLEN .LE. QRDIM )THEN
C
C            call QR
C
         ELSE
*           Do MRRR
*           Build the representation and init eigenvalues

*           The decision to compute the eigenvalues with dqds or not
*           may not depend on the local index range of wanted eigen-
*           pairs, otherwise we lose constistency. The sole exception
*           is if all ews in the block are wanted (even if for the full
*           matrix only a subset is desired).
            IF( BWIL.EQ.1 .AND. BWIU.EQ.BLEN .AND. USEDQD )THEN
               EMODE = 'd'
            ELSE
               EMODE = 'o'
            ENDIF

            CALL DLAXRE(
     $             BLEN, RWORK(IXD-1 + BBEG), RWORK(IXE-1 + BBEG),
     $             BGERL, BGERU, ASPTOL, GAPTOL,
     $             RWORK(IXREP), IWORK(IYREP), BSHIFT,
     $             IWORK(IYEWAE), RWORK(IXEWLU), EMODE,
     $             RWORK(IXWORK), IINFO
     $           )
            IF( IINFO .NE. 0 )THEN
               DO IXF77A = 2*ZA-1, 2*ZE-1, 2
                  ISUPPZ(IXF77A) = 0
               ENDDO
               DO IXF77A = 2*ZA, 2*ZE, 2
                  ISUPPZ(IXF77A) = -3
               ENDDO
               CYCLE
            ENDIF


*           Compute the vectors
            CALL DLAXRV(
     $             BLEN, RWORK(IXE-1 + BBEG),
     $             RWORK(IXREP), IWORK(IYREP),
     $             IWORK(IYEWAE), RWORK(IXEWLU), BWIL, BWIU,
     $             BSDIAM, GAPTOL,
     $             W(ZA), Z(BBEG,ZA), LDZ, ISUPPZ(2*ZA-1),
     $             RWORK(IXWORK), LRWREM, IWORK(IYWORK), LIWREM,
     $             IINFO
     $           )
            GOTALL = (IINFO .EQ. 0)
            IF( IINFO.NE.0 .AND. IINFO.NE.1 )THEN
               INFO = 10
               RETURN
            ENDIF

*           Undo the shift from DLAXRE
            DO J = ZA, ZE
               W(J) = W(J) + BSHIFT
            ENDDO
         ENDIF
*
*        -----------------------------------
*         Adjust results to original matrix
*        -----------------------------------
*
*        Eigenvalues:
*        - Cap eigenvalues by the separators
*        - Undo the scaling from DLAXRA
*
*        Eigenvectors:
*        - Undo the outer scaling by signature matrix in DLAXRA
*        - Adjust support to full matrix


         DO J = ZA, ZE
            IF( W(J) .LT. BSEPL )THEN
               W(J) = BSEPL
            ELSEIF( W(J) .GT. BSEPU )THEN
               W(J) = BSEPU
            ENDIF
            W(J) = W(J) * INVSC
*
            IF( ISUPPZ(2*J-1) .GT. 0 )THEN
               ISUPPZ(2*J-1) = ISUPPZ(2*J-1) + BBEG-1
               ISUPPZ(2*J)   = ISUPPZ(2*J)   + BBEG-1
            ENDIF
         ENDDO

         DO I = BBEG, BEND
            IF( IWORK(IYS-1 + I) .NE. +1 )THEN
*              negate row of Z
               DO J = ZA, ZE
                  Z(I,J) = -Z(I,J)
               ENDDO
            ENDIF
         ENDDO

      ENDDO

      TSTAMP1 = DSECND()
      XTIME2  = XTIME2 + (TSTAMP1 - TSTAMP0)
      TSTAMP0 = TSTAMP1

*     ===================================================================
*                Sort the eigenpairs into ascending order
*     ===================================================================

      CALL DLAXRO( N, M, W, Z, LDZ, ISUPPZ, IWORK )


      INFO = 0
      IF( .NOT. GOTALL )  INFO = 1
*
      TSTAMP1 = DSECND()
      XTIME3  = XTIME3 + (TSTAMP1 - TSTAMP0)
      TSTAMP0 = TSTAMP1
      END SUBROUTINE DSTEXR
*
*************************************************************************
