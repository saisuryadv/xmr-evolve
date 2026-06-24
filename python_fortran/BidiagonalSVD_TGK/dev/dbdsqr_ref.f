      PROGRAM DBDSQR_REF
*
*     Reference DBDSQR runner -- reads an (upper) bidiagonal matrix from
*     an STCollection-format .dat file (line 1 = N, lines 2..N+1 = IDX
*     D(i) E(i)) and computes only its singular values via DBDSQR with
*     NCVT = NRU = NCC = 0.  Used as the timing-and-accuracy reference
*     for the self-contained orchestrator path (eval_dbdsvdmr3.py).
*
*     Output (parseable):
*       INFO=...  T_SEC=<min-of-N seconds>  NREPS=<reps actually run>
*       <i>  <sigma_i>     (i = 1..n, descending order)
*
*     Benchmark protocol:
*       - DBDSQR overwrites D and E, so save them and restore each rep.
*       - One untimed warmup call (J.I.T. paging / library load).
*       - Adaptive loop: stop when either (a) the very first rep took
*         >= 0.5 s, (b) total bench time >= 0.2 s and NREPS >= 3, or
*         (c) NREPS >= 30.  Headline T_SEC = MIN across reps -- the
*         standard "best representative single sample" measure that is
*         least polluted by OS interrupts / scheduling jitter.
*       - Cache flush: between each timed call, walk a 32 MiB scratch
*         buffer (write+read every cache line) to evict the bidiagonal
*         input and DBDSQR's workspace from L1/L2/L3.  The clock starts
*         only AFTER the flush, so reported time is SVD work only.
*       - SYSTEM_CLOCK is called with INTEGER(KIND=8) ticks for the
*         highest resolution gfortran provides (typically nanoseconds
*         under glibc's CLOCK_MONOTONIC).
*
*     Recommended invocation: taskset -c 0 dbdsqr_ref input.dat
*     (so the OS doesn't migrate the process between cores mid-call).
*
      IMPLICIT NONE
      INTEGER, PARAMETER :: JUNK_SIZE = 4*1024*1024
      INTEGER N, INFO, I, IDX, IOS, NREPS, J
      INTEGER*8 T1, T2, TR
      DOUBLE PRECISION DD, EE, T_MIN, T_TOTAL, T_THIS, JSUM
      DOUBLE PRECISION, ALLOCATABLE :: D(:), E(:), DS(:), ES(:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), JUNK(:)
      DOUBLE PRECISION DUMMY(1,1)
      CHARACTER*4096 FNAME
      EXTERNAL DBDSQR
      INTRINSIC MIN, MAX, DBLE
*
      CALL GET_COMMAND_ARGUMENT( 1, FNAME )
      OPEN( UNIT=10, FILE=FNAME, STATUS='OLD', IOSTAT=IOS )
      IF( IOS.NE.0 ) THEN
         WRITE(*,*) 'cannot open ', TRIM(FNAME)
         STOP 1
      END IF
      READ(10,*) N
      IF( N.LT.1 ) THEN
         WRITE(*,*) 'bad N=', N
         STOP 1
      END IF
      ALLOCATE( D(N), E(N), DS(N), ES(N), WORK(4*N) )
      DO I = 1, N
         READ(10,*) IDX, DD, EE
         DS(I) = DD
         ES(I) = EE
      END DO
      CLOSE(10)
      ES(N) = 0.0D0
*
*     Allocate junk array to flush L1/L2/L3 between reps.
      ALLOCATE( JUNK(JUNK_SIZE) )
      DO I = 1, JUNK_SIZE
         JUNK(I) = DBLE(I)
      END DO
*
*     Warmup (untimed).
      D = DS
      E = ES
      CALL DBDSQR( 'U', N, 0, 0, 0, D, E, DUMMY, 1, DUMMY, 1,
     $             DUMMY, 1, WORK, INFO )
*
*     Adaptive benchmark loop.
      NREPS = 0
      T_MIN = 1.0D30
      T_TOTAL = 0.0D0
      JSUM = 0.0D0
   10 CONTINUE
         D = DS
         E = ES
*        Cache flush: write+read every cache line in 32 MiB scratch.
         DO J = 1, JUNK_SIZE
            JUNK(J) = JUNK(J) * 1.0000001D0 + 1.0D-30
            JSUM = JSUM + JUNK(J)
         END DO
         CALL SYSTEM_CLOCK( T1, TR )
         CALL DBDSQR( 'U', N, 0, 0, 0, D, E, DUMMY, 1, DUMMY, 1,
     $                DUMMY, 1, WORK, INFO )
         CALL SYSTEM_CLOCK( T2 )
         T_THIS = DBLE( T2 - T1 ) / DBLE( TR )
         IF( T_THIS.LT.T_MIN ) T_MIN = T_THIS
         T_TOTAL = T_TOTAL + T_THIS
         NREPS = NREPS + 1
*        Stop conditions
         IF( NREPS.EQ.1 .AND. T_THIS.GE.0.5D0 ) GOTO 20
         IF( NREPS.GE.3 .AND. T_TOTAL.GE.0.2D0 ) GOTO 20
         IF( NREPS.GE.30 ) GOTO 20
         GOTO 10
   20 CONTINUE
*
*     Anti-DCE: print something only on impossible value of JSUM.
      IF( JSUM.EQ.1.2345D-300 ) WRITE(*,*) 'never', JUNK(1)
*
      WRITE(*,900) INFO, T_MIN, NREPS
  900 FORMAT('INFO=',I6,' T_SEC=',ES17.10,' NREPS=',I0,
     $       ' CACHE_FLUSH=ON')
      DO I = 1, N
         WRITE(*,910) I, D(I)
      END DO
  910 FORMAT(I8,1X,ES23.16)
      DEALLOCATE( D, E, DS, ES, WORK, JUNK )
      STOP
      END
