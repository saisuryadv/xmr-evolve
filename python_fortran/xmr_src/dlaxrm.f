      SUBROUTINE  DLAXRM( N, REPR, REPI, M, ATAU, AXI )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, M
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), ATAU(M)
*
      INTEGER,          INTENT(OUT)  ::  AXI(M)
*
*  Purpose
*  =======
*
*    Compute multiple sturm counts with respect to the same
*    representation. The routine computes for each shift tau
*    in ATAU the inertia of T - tau into AXI---see dlaxrn for
*    information on how the inertias are to be interpreted.
*
*    Pre: No zero shift.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRM_STAT64(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(64)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(64)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(64)
      END SUBROUTINE DLAXRM_STAT64
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRM_STAT32(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(32)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(32)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(32)
      END SUBROUTINE DLAXRM_STAT32
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRM_STAT16(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(16)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(16)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(16)
      END SUBROUTINE DLAXRM_STAT16
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRM_STAT8(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(8)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(8)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(8)
      END SUBROUTINE DLAXRM_STAT8
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRM_STAT4(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(4)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(4)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(4)
      END SUBROUTINE DLAXRM_STAT4
      END INTERFACE
      INTERFACE
      SUBROUTINE DLAXRM_STAT2(
     $             N, K, DIR, G, IBB, LBBEGK, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, ANEGC, AAUX
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, K, DIR
      INTEGER,          INTENT(IN)  ::  LBBEGK(*)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU(2)
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
*
      INTEGER,          INTENT(INOUT)  ::  IBB
*
      INTEGER,          INTENT(OUT)  ::  ANEGC(2)
      DOUBLE PRECISION, INTENT(OUT)  ::  AAUX(2)
      END SUBROUTINE DLAXRM_STAT2
      END INTERFACE
      INTERFACE
      FUNCTION DLAXRN(N, REPR, REPI, TAU)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN
      INTEGER,          INTENT(IN)  ::  N
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU
      END FUNCTION DLAXRN
      END INTERFACE


Cextract -b parameters.inc.f procname=laxrm
*
*     The maximal number of negcounts that are to be performed at
*     once. Controls which of the routines dlaxrm_statXX are allowed.
*     As such, sensible values are those for which a corresponding
*     routine exists in the first place (1,2,4,8,16,32,64 atm).
*     Must be >= 1.
*
      INTEGER, PARAMETER  ::  MAXPARNEG = 2

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0
*
*     .. Local Variables ..
*
      DOUBLE PRECISION  PIVBASE, PIVMIN, GMAREP, GMA
      INTEGER           J, I, IBB, K, NB, XI
      INTEGER           IXG, IXBDET, IXNGN, IXGNSQ, JYLBBK

      DOUBLE PRECISION  AUXN( 64 ), AUXP( 64 )
      INTEGER           NEGN( 64 ), NEGP( 64 )


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
*     -- Executable Statements -----------------------------------------
*

*     Note: Here we use absolute positions to the beginning of data.
      IXG    = 1 + (0)
      IXBDET = 1 + (N)
      IXNGN  = 1 + (2*N+1)
      IXGNSQ = 1 + (3*N+2)

      JYLBBK = (N+6)

      PIVBASE = REPR(4*N+3)
      K       = REPI(2)
      NB      = REPI(3)

      GMAREP = REPR(IXG-1 + K)

      PIVMIN = PIVBASE


      J = 1

      IF( MAXPARNEG.GE.64 )THEN
         DO
            IF( J + 63 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT64(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT64(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 64
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF

      IF( MAXPARNEG.GE.32 )THEN
         DO
            IF( J + 31 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT32(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT32(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 32
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF

      IF( MAXPARNEG.GE.16 )THEN
         DO
            IF( J + 15 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT16(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT16(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 16
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF

      IF( MAXPARNEG.GE.8 )THEN
         DO
            IF( J + 7 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT8(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT8(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 8
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF

      IF( MAXPARNEG.GE.4 )THEN
         DO
            IF( J + 3 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT4(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT4(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 4
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF


      IF( MAXPARNEG.GE.2 )THEN
         DO
            IF( J + 1 .GT. M )  EXIT

            IBB = 1
            CALL DLAXRM_STAT2(
     $              N, K, +1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGN, AUXN
     $           )

            IBB = NB+1
            CALL DLAXRM_STAT2(
     $              N, K, -1, REPR(IXG), IBB, REPI(JYLBBK),
     $              REPR(IXGNSQ), REPR(IXNGN), REPR(IXBDET),
     $              PIVMIN, ATAU(J), NEGP, AUXP
     $           )

            DO I = 1, 2
               GMA = GMAREP + ((AUXN(I) + AUXP(I)) - ATAU(J))
               XI  = 2*(NEGN(I) + NEGP(I))
               IF( GMA.LT.ZERO )THEN
                  XI = XI + 2
               ELSEIF( GMA.EQ.ZERO )THEN
                  XI = XI + 1
               ENDIF
               AXI(J) = XI
               J = J+1
            ENDDO
         ENDDO
      ENDIF


      IF( .NOT. XSTEALTHMODE )  XNUMFN = XNUMFN + J-1

*     Do all remaining with the standard single negcount
      DO
         IF( J .GT. M )  EXIT

         AXI(J) = DLAXRN(N, REPR, REPI, ATAU(J))
         J = J+1
      ENDDO

      END SUBROUTINE DLAXRM
*
************************************************************************


