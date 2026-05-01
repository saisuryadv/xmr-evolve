      SUBROUTINE DLAXRF_COB(
     $             N, FLGAP, FLB, FUB, FUGAP, REPR, REPI, TAU, GAPTOL,
     $             ICBEG, ICEND, XIZERO, SLGAP, SUGAP, EWL_AE, EWL_LU,
     $             STATUS
     $         )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND, XIZERO
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  FLGAP, FLB, FUB, FUGAP
      DOUBLE PRECISION, INTENT(IN)  ::  TAU, GAPTOL
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
*
      INTEGER,          INTENT(OUT)  ::  STATUS
      INTEGER,          INTENT(OUT)  ::  EWL_AE(2*ICBEG-1:2*ICEND)
      DOUBLE PRECISION, INTENT(OUT)  ::  SLGAP, SUGAP
      DOUBLE PRECISION, INTENT(OUT)  ::  EWL_LU(2*ICBEG-1:2*ICEND)
*
*  Purpose
*  =======
*
*    Check Outer Bounds: Determine initial outer bounds, verify that
*    there is no eigenvalue underflow, and initialize a sons ewi-list.
*
*    Failure of some checks is indicated by STATUS != 0, the value then
*    indicates the kind of problem:
*        -1|+1  ew underflow on negative|positive side
*        -2|+2  two very small (about eps*tau) nonzero ews
*        -3|+3  could not verify left|right near bound
*        -5|+5  could not verify left|right far bound (the gaps)
*
*    If STATUS=0, all test are passed. Then the sons ew-list will be
*    initialized properly.
*
*    Independent of STATUS, XIZERO is always set to the sons zero
*    inertia.
*
*  ======================================================================
*
*     .. Declarations ..
*
      INTERFACE
      SUBROUTINE DLAXRL_RESET(
     $             IL, IU, LGAP, UGAP, EWL_AE, EWL_LU, L, XIL, U, XIU
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  IL, IU, XIL, XIU
      DOUBLE PRECISION, INTENT(IN)  ::  L, U
*
      INTEGER,          INTENT(INOUT)  ::  EWL_AE(2*IL-1:2*IU)
      DOUBLE PRECISION, INTENT(INOUT)  ::  LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT)  ::  EWL_LU(2*IL-1:2*IU)
      END SUBROUTINE DLAXRL_RESET
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
      INTERFACE
      FUNCTION DLAXRN(N, REPR, REPI, TAU)
      IMPLICIT NONE
      INTEGER  ::  DLAXRN
      INTEGER,          INTENT(IN)  ::  N
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3), TAU
      END FUNCTION DLAXRN
      END INTERFACE

      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO = 0.0D0

*     The maximal number of bounds we generate to init the sons ew-list.
      INTEGER, PARAMETER  ::  MAXNBOUNDS = 15
*
*     .. Parameters ..
*

*     Father and son bounds are relaxed before checking them by
*         2**LOGFBNDRELAXFAC * (N*EPS)
*
      INTEGER, PARAMETER  ::  LOGFBNDRELAXFAC = 7

*
*     .. Locals ..
*
      DOUBLE PRECISION  EPS, PIVBASE, UTOL, L, LX, WRELGAP
      DOUBLE PRECISION  FRELAX, FBOUND, SBOUND, FBX, SBX1, SBX2
      DOUBLE PRECISION  FRELTIGHT
      DOUBLE PRECISION  ABNDS( MAXNBOUNDS )

      INTEGER           TYPE, XI, INEG, IPOS
      INTEGER           I, ISL, ISU
      INTEGER           AINRS( MAXNBOUNDS ), NBNDS

      LOGICAL  OK, GOTBND

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
*  ===== Executable Statements ==========================================
*

      STATUS = 0
      NBNDS  = 0

      EPS  = DLAMCH('Epsilon')
      PIVBASE = REPR(4*N+3)

      UTOL = PIVBASE * (N / EPS)
      L    = EPS*ABS(TAU)
      LX   = 2*L

      FRELAX  = 2**LOGFBNDRELAXFAC * (N * EPS)
      WRELGAP = GAPTOL
      FRELTIGHT = 4 * (N * EPS)
*
*     Determine the type of this shift based on the zero inertia.
*
      NBNDS = NBNDS + 1
      ABNDS(NBNDS) = ZERO
      AINRS(NBNDS) = XIZERO

*     set IPOS as index of first positive ew
*     set INEG as index of first negative ew
      IPOS = (XIZERO+1) / 2 + 1
      INEG = XIZERO / 2

      IF( XIZERO .LE. 2*ICBEG-1 )THEN
         TYPE = +1
      ELSEIF( XIZERO .GE. 2*ICEND-1 )THEN
         TYPE = -1
      ELSE
         TYPE = 0
      ENDIF

*
*     -------------------------------
*      Check positive very small ews
*     -------------------------------
*
      IF( TYPE.NE.-1 .AND. STATUS.EQ.0 )THEN
         XI = DLAXRN( N, REPR, REPI, L )
         NBNDS = NBNDS + 1
         ABNDS(NBNDS) = L
         AINRS(NBNDS) = XI

         IF( .NOT.(XI .LE. 2*IPOS-1) )THEN
*
*           L is not a lower bound for the first positive ew
*
            XI = DLAXRN( N, REPR, REPI, UTOL )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = UTOL
            AINRS(NBNDS) = XI
            IF( .NOT.(XI.EQ.2*IPOS-2 .OR. XI.EQ.2*IPOS-1) )THEN
*              underflow on positive side
               STATUS = +1
            ELSEIF( IPOS .LT. ICEND )THEN
               XI = DLAXRN( N, REPR, REPI, LX )
               NBNDS = NBNDS + 1
               ABNDS(NBNDS) = LX
               AINRS(NBNDS) = XI
               IF( .NOT.(XI .LE. 2*(IPOS+1)-1) )THEN
*                 we more than one very small ew between L and LX
                  STATUS = +2
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
*     -------------------------------
*      Check negative very small ews
*     -------------------------------
*
      IF( TYPE.NE.+1 .AND. STATUS.EQ.0 )THEN
         XI = DLAXRN( N, REPR, REPI, -L )
         NBNDS = NBNDS + 1
         ABNDS(NBNDS) = -L
         AINRS(NBNDS) = XI
         IF( .NOT.(XI .GE. 2*INEG-1) )THEN
*
*          -L is not an upper bound for first negative ew
*
            XI = DLAXRN( N, REPR, REPI, -UTOL )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = -UTOL
            AINRS(NBNDS) = XI
            IF( .NOT.(XI.EQ.2*INEG-1 .OR. XI.EQ.2*INEG) )THEN
*              underflow on negative side
               STATUS = -1
            ELSEIF( INEG .GT. ICBEG )THEN
               XI = DLAXRN( N, REPR, REPI, -LX )
               NBNDS = NBNDS + 1
               ABNDS(NBNDS) = -LX
               AINRS(NBNDS) = XI
               IF( .NOT.(XI .GE. 2*(INEG-1)-1) )THEN
*                 more than one very small ew between -LX and -L
                  STATUS = -2
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
*     ------------------
*      Check left bound
*     ------------------
*
      IF( STATUS .EQ. 0 )THEN
         GOTBND = .FALSE.

         FBOUND = FLB - ABS(FLB)*FRELTIGHT
         IF( TYPE.EQ.+1 .AND. TAU.GE.FBOUND )THEN
            GOTBND = .TRUE.
         ELSE
            SBOUND = (FBOUND-TAU) - ABS(FBOUND-TAU)*FRELTIGHT
            XI = DLAXRN( N, REPR, REPI, SBOUND )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBOUND
            AINRS(NBNDS) = XI
            IF( XI.EQ.2*ICBEG-2 .OR. XI.EQ.2*ICBEG-1 )THEN
               GOTBND = .TRUE.
            ENDIF
         ENDIF

         FBOUND = FLB - ABS(FLB)*FRELAX
         IF( .NOT.GOTBND .AND. .NOT.(TYPE.EQ.+1 .AND. TAU.GE.FBOUND) )
     $   THEN
            SBOUND = (FBOUND-TAU) - ABS(FBOUND-TAU)*FRELAX
            XI = DLAXRN( N, REPR, REPI, SBOUND )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBOUND
            AINRS(NBNDS) = XI
            IF( .NOT.(XI.EQ.2*ICBEG-2 .OR. XI.EQ.2*ICBEG-1) )THEN
*              Cannot verify left inner far bound. If this shift is
*              to the left we can live with that, since 0 is a bound.
               IF( TYPE .NE. +1 )THEN
                  STATUS = -3
               ENDIF
            ENDIF
         ENDIF

C         FBOUND = FLB - ABS(FLB)*FRELAX
C
C         IF( .NOT.(TYPE.EQ.+1 .AND. TAU.GE.FBOUND) )THEN
C
C            SBOUND = (FBOUND-TAU) - ABS(FBOUND-TAU)*FRELAX
C            XI = DLAXRN( N, REPR, REPI, SBOUND )
C            NBNDS = NBNDS + 1
C            ABNDS(NBNDS) = SBOUND
C            AINRS(NBNDS) = XI
C
C            IF( .NOT.(XI.EQ.2*ICBEG-2 .OR. XI.EQ.2*ICBEG-1) )THEN
C*              cannot verify left inner far bound
C               STATUS = -3
C            ENDIF
C         ENDIF
      ENDIF
*
*     -------------------
*      Check right bound
*     -------------------
*
      IF( STATUS .EQ. 0 )THEN
         GOTBND = .FALSE.

         FBOUND = FUB + ABS(FUB)*FRELTIGHT
         IF( TYPE.EQ.-1 .AND. TAU.LE.FBOUND )THEN
            GOTBND = .TRUE.
         ELSE
            SBOUND = (FBOUND-TAU) + ABS(FBOUND-TAU)*FRELTIGHT
            XI = DLAXRN( N, REPR, REPI, SBOUND )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBOUND
            AINRS(NBNDS) = XI

            IF( XI.EQ.2*ICEND-1 .OR. XI.EQ.2*ICEND )THEN
               GOTBND = .TRUE.
            ENDIF
         ENDIF

         FBOUND = FUB + ABS(FUB)*FRELAX
         IF( .NOT.GOTBND .AND. .NOT.(TYPE.EQ.-1 .AND. TAU.LE.FBOUND) )
     $   THEN
            SBOUND = (FBOUND-TAU) + ABS(FBOUND-TAU)*FRELAX
            XI = DLAXRN( N, REPR, REPI, SBOUND )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBOUND
            AINRS(NBNDS) = XI

            IF( .NOT.(XI.EQ.2*ICEND-1 .OR. XI.EQ.2*ICEND) )THEN
*              Cannot verifY right inner bound. If this shift is
*              to the right we can live with that, since 0 is a bound.
               IF( TYPE .NE. +1 )THEN
                  STATUS = +3
               ENDIF
            ENDIF
         ENDIF

C         FBOUND = FUB + ABS(FUB)*FRELAX
C         IF( .NOT.(TYPE.EQ.-1 .AND. TAU.LE.FBOUND) )THEN
C
C            SBOUND = (FBOUND-TAU) + ABS(FBOUND-TAU)*FRELAX
C            XI = DLAXRN( N, REPR, REPI, SBOUND )
C            NBNDS = NBNDS + 1
C            ABNDS(NBNDS) = SBOUND
C            AINRS(NBNDS) = XI
C
C            IF( .NOT.(XI.EQ.2*ICEND-1 .OR. XI.EQ.2*ICEND) )THEN
C*              cannot verifY right inner far bound
C               STATUS = +3
C            ENDIF
C         ENDIF
      ENDIF
*
*     -----------------------
*      Check left outer gap
*     -----------------------
*
      IF( STATUS.EQ.0 .AND. TYPE.NE.+1 )THEN
         FBOUND = FLB - ABS(FLB)*FRELAX
         FBX    = (FLB-FLGAP) + ABS(FLB-FLGAP)*FRELAX
         SBOUND = (FBOUND-TAU) - ABS(FBOUND-TAU)*FRELAX
         SBX1   = MAX( 2*SBOUND, (FBX-TAU) + ABS(FBX-TAU)*FRELAX )
         SBX2   = SBOUND - WRELGAP*ABS(SBOUND)
         OK = .FALSE.

*        1st try: go wide
         IF( SBX1 .LT. SBX2 )THEN
            XI = DLAXRN( N, REPR, REPI, SBX1 )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBX1
            AINRS(NBNDS) = XI

            OK = ( XI.EQ.2*(ICBEG-1)-1 .OR. XI.EQ.2*(ICBEG-1) )
         ENDIF

*        2nd try: go tight
         IF( .NOT. OK )THEN
            XI = DLAXRN( N, REPR, REPI, SBX2 )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBX2
            AINRS(NBNDS) = XI

            OK = ( XI.EQ.2*(ICBEG-1)-1 .OR. XI.EQ.2*(ICBEG-1) )
         ENDIF

         IF( .NOT. OK )THEN
*           Could not establish sons left outer gap
            STATUS = -5
         ENDIF
      ENDIF
*
*     -----------------------
*      Check right outer gap
*     -----------------------
*
      IF( STATUS.EQ.0 .AND. TYPE.NE.-1 )THEN
         FBOUND = FUB + ABS(FUB)*FRELAX
         FBX    = (FUB+FUGAP) - ABS(FUB+FUGAP)*FRELAX
         SBOUND = (FBOUND-TAU) + ABS(FBOUND-TAU)*FRELAX
         SBX1   = MIN( 2*SBOUND, (FBX-TAU) - ABS(FBX-TAU)*FRELAX )
         SBX2   = SBOUND + WRELGAP*ABS(SBOUND)
         OK = .FALSE.

*        1st try: go wide
         IF( SBX1 .GT. SBX2 )THEN
            XI = DLAXRN( N, REPR, REPI, SBX1 )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBX1
            AINRS(NBNDS) = XI

            OK = ( XI.EQ.2*ICEND .OR. XI.EQ.2*ICEND+1 )
         ENDIF

*        2nd try: go tight
         IF( .NOT. OK )THEN
            XI = DLAXRN( N, REPR, REPI, SBX2 )
            NBNDS = NBNDS + 1
            ABNDS(NBNDS) = SBX2
            AINRS(NBNDS) = XI

            OK = ( XI.EQ.2*ICEND .OR. XI.EQ.2*ICEND+1 )
         ENDIF

         IF( .NOT. OK )THEN
*           Could not establish sons right outer gap
            STATUS = +5
         ENDIF
      ENDIF
C     Don't count xizero
      XNBIS_COB = XNBIS_COB + NBNDS-1
*
*     ------------------------
*      Init the sons ewi-list
*     ------------------------
*
      IF( STATUS .EQ. 0 )THEN
         ISL = 0
         ISU = 0
         DO I = 1, NBNDS
            XI = AINRS(I)
            IF( XI .LE. 2*ICBEG-1 )  ISL = I
            IF( XI .GE. 2*ICEND-1 )  ISU = I
         ENDDO

         CALL DLAXRL_RESET(
     $          ICBEG, ICEND, SLGAP, SUGAP, EWL_AE, EWL_LU,
     $          ABNDS(ISL), AINRS(ISL), ABNDS(ISU), AINRS(ISU)
     $        )

         DO I = 1, NBNDS
            IF( I.NE.ISL .AND. I.NE.ISU )THEN
               CALL DLAXRL_UPDATE(
     $                ICBEG, ICEND, SLGAP, SUGAP, EWL_AE, EWL_LU,
     $                ABNDS(I), AINRS(I)
     $              )
            ENDIF
         ENDDO

*        If we have a zero bound, the outer gaps were not checked since
*        we know there is a relative gap, the outer gaps for the son
*        will then be zero. To avoid zero gaps, we just copy the one
*        from the father. The gap are not used anyway, since dlaxrx has
*        special handling for zero eigenvalues.
         IF( EWL_LU(2*ICBEG-1).EQ.ZERO )THEN
            SLGAP = FLGAP
         ENDIF
         IF( EWL_LU(2*ICEND).EQ.ZERO )THEN
            SUGAP = FUGAP
         ENDIF

      ENDIF


      END SUBROUTINE DLAXRF_COB
*
*************************************************************************
