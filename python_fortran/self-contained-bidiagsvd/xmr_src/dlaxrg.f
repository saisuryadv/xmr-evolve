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
*
*  Purpose
*  =======
*
*     Computes the FP-vector for a zero eigenvalue, normed to unity.
*
*     Concerning the block structure, there should be no block ending
*     at k, regardless of k=1, k=n or in-between. This could be fixed,
*     but atm we only compute vectors with blocks for a singular matrix,
*     and those should have a single pivot being zero.
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL DSCAL

      EXTERNAL         DLAMCH
      DOUBLE PRECISION DLAMCH

*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::   ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::   ONE  = 1.0D0
*
*     .. Locals ..
*
      DOUBLE PRECISION  UTHRESH, L, OLDL, U, OLDU, AZ, PREVZ
      DOUBLE PRECISION  NORMSQ, INVNRM
      INTEGER           JXG, JXBDET, JYOMGA, K
      INTEGER           I, NZFIX, IACTZ

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
      IF( .NOT. XSTEALTHMODE )THEN
         XNUMGV0 = XNUMGV0 + 1
      ENDIF

      JXG    = (0)
      JXBDET = (N)
      JYOMGA = (4)
      K      = REPI(2)


      UTHRESH = DLAMCH('Underflow')

      Z(K) = ONE
      NORMSQ = ONE
      NZFIX = 0

*
*     Solve from twist on down to 1
*
      PREVZ = ONE

      OLDL = ZERO
      I    = K-1
      DO
         IF( I.EQ.0 )THEN
            EXIT
         ENDIF

* assume we have PREVZ such that Z(i+1) was computed as -OLDL*PREVZ
* if i+2 ends a block then PREVZ=Z(i+3), otherwise PREVZ=Z(i+2)
CCC
         IF( REPI(JYOMGA + I).NE.0 )THEN
*           I ends a block, implies i < k here --> L(i) = E(i)/invD(i)
            L = (E(I) * REPR(JXG + I-1)) / REPR(JXBDET + I-1)
            IACTZ = I+1
         ELSEIF( REPI(JYOMGA + I+1).NE.0 )THEN
*           I starts a block, implies i < k-1 here
            L = - (E(I+1) * E(I)) / REPR(JXBDET + I)
            IACTZ = I+2
         ELSE
            L = E(I) / REPR(JXG + I)
            IACTZ = I+1
         ENDIF
*
*        Compute next vector entry, handle possible creeping underflow
*
         IF( ABS(Z(IACTZ)).LE.UTHRESH )THEN
*           We should have I < K-1.
*           We also have |OLDL| < 1, |PREVZ| < 1, thus if there is
*           underflow anywhere in the following product, the exact
*           result would underflow as well.
            Z(IACTZ) = ZERO
            Z(I)   = (L * OLDL) * PREVZ
            NZFIX  = NZFIX+1
         ELSE
            Z(I)   = - L * Z(IACTZ)
         ENDIF
*
*        support cut
*
         AZ = ABS(Z(I)) + ABS(Z(I+1))
         IF( AZ.LE.UTHRESH .OR. AZ*E(I).LE.CUTTOL )THEN
            Z(1:I) = ZERO
            EXIT
         ELSE
            NORMSQ = NORMSQ + Z(I)**2
         ENDIF

         PREVZ = Z(IACTZ)
         OLDL  = L

         I = I-1
      ENDDO
      IF( ABS(Z(1)).LE.UTHRESH ) Z(1) = ZERO
*
*     Solve from twist on down to n
*     (completely analogous to upper part)
*
      PREVZ = ONE
      OLDU = ZERO
      I    = K+1
      DO
         IF( I.EQ.N+1 )THEN
            EXIT
         ENDIF

         IF( REPI(JYOMGA + I).NE.0 )THEN
*           I ends a block, implies i > k here --> U(i) = E(i-1)/invR(i)
            U = (E(I-1) * REPR(JXG + I+1)) / REPR(JXBDET + I+1)
            IACTZ = I-1
         ELSEIF( REPI(JYOMGA + I-1).NE.0 )THEN
*           I starts a block, implies i > k+1 here
            U = - (E(I-2) * E(I-1)) / REPR(JXBDET + I)
            IACTZ = I-2
         ELSE
            U = E(I-1) / REPR(JXG + I)
            IACTZ = I-1
         ENDIF

         IF( ABS(Z(IACTZ)).LE.UTHRESH )THEN
            Z(IACTZ) = ZERO
            Z(I)   = (U * OLDU) * PREVZ
            NZFIX  = NZFIX+1
         ELSE
            Z(I)   = - U * Z(IACTZ)
         ENDIF

         AZ = ABS(Z(I)) + ABS(Z(I-1))
         IF( AZ.LE.UTHRESH .OR. AZ*E(I-1).LE.CUTTOL )THEN
            Z(I:N) = ZERO
            EXIT
         ELSE
            NORMSQ = NORMSQ + Z(I)**2
         ENDIF

         PREVZ = Z(IACTZ)
         OLDU  = U

         I    = I+1
      ENDDO
      IF( ABS(Z(N)).LE.UTHRESH ) Z(N) = ZERO
*
*     Determine results
*
      INVNRM = ONE / SQRT(NORMSQ)
*
*     Normalize Z
*
      CALL DSCAL( N, INVNRM, Z, 1 )

      END SUBROUTINE DLAXRG0
*
************************************************************************
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
*
*  Purpose
*  =======
*
*     Computes the FP-vector, normed to unity.
*
*     We require that the factorization exists, in particular the
*     elements L(i) and U(i) should be finite and nonzero.
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL DSCAL

      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH

      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

*
*     .. Locals ..
*
      DOUBLE PRECISION UTHRESH, L, OLDL, U, OLDU, AZ, NORMSQ, INVNRM
      INTEGER I, NZFIX, ISA, ISE

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
      IF( .NOT. XSTEALTHMODE )THEN
         XNUMGV = XNUMGV + 1
      ENDIF

      UTHRESH = DLAMCH('Underflow')

      Z(K) = ONE
      NORMSQ = ONE
      NZFIX = 0
      ISA = 1
      ISE = N

*
*     Solve from twist on down to 1
*
      OLDL = ZERO
      I    = K-1
      DO
         IF( I.EQ.0 )THEN
            EXIT
         ENDIF

         L = E(I) / D(I)
*
*        Compute next vector entry, handle possible creeping underflow
*
         IF( ABS(Z(I+1)).LE.UTHRESH )THEN
*           We should have I < K-1.
*           We also have |OLDL| < 1, |Z(i+2)| < 1, thus if there is
*           underflow anywhere in the following product, the exact
*           result would underflow as well.
            Z(I+1) = ZERO
            Z(I)   = (L * OLDL) * Z(I+2)
            NZFIX  = NZFIX+1
         ELSE
            Z(I)   = - L * Z(I+1)
         ENDIF
*
*        support cut
*
         AZ = ABS(Z(I)) + ABS(Z(I+1))
         IF( AZ.LE.UTHRESH .OR. AZ*E(I).LE.CUTTOL )THEN
            Z(1:I) = ZERO
            ISA = I+1
            EXIT
         ELSE
            NORMSQ = NORMSQ + Z(I)**2
         ENDIF

         OLDL = L
         I = I-1
      ENDDO
*
*     Solve from twist on down to n
*     (completely analogous to upper part)
*
      OLDU = ZERO
      I    = K+1
      DO
         IF( I.EQ.N+1 )THEN
            EXIT
         ENDIF

         U = E(I-1) / R(I)
         IF( ABS(Z(I-1)).LE.UTHRESH )THEN
            Z(I-1) = ZERO
            Z(I)   = (U * OLDU) * Z(I-2)
            NZFIX  = NZFIX+1
         ELSE
            Z(I)   = - U * Z(I-1)
         ENDIF

         AZ = ABS(Z(I)) + ABS(Z(I-1))
         IF( AZ.LE.UTHRESH .OR. AZ*E(I-1).LE.CUTTOL )THEN
            Z(I:N) = ZERO
            ISE = I-1
            EXIT
         ELSE
            NORMSQ = NORMSQ + Z(I)**2
         ENDIF

         OLDU = U
         I    = I+1
      ENDDO
*
*     Determine results
*
      NORMZ  = SQRT(NORMSQ)
      INVNRM = ONE / NORMZ

      RESID  = ABS(GAMMA) * INVNRM
      RQCORR = GAMMA / NORMSQ
*
*     Normalize Z
*
      CALL DSCAL( N, INVNRM, Z, 1 )

      END SUBROUTINE DLAXRG
*
************************************************************************
