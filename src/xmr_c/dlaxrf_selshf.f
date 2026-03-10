      SUBROUTINE DLAXRF_SELSHF(
     $             N, REPR, REPI,
     $             ICBEG, ICEND, LGAP, UGAP, EWL_AE, EWL_LU, RGINFO,
     $             TAUBAR, GAPTOL, MAXCPO, MAXCPI, MAXNC, MAXNB,
     $             NCAND, ACLOC, ACTAU, NBATCH, ABEND,
     $             IWORK
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, ICBEG, ICEND
      INTEGER,          INTENT(IN)  ::  MAXCPO, MAXCPI, MAXNC, MAXNB
      INTEGER,          INTENT(IN)  ::  RGINFO(ICBEG-1:ICEND)
      INTEGER,          INTENT(IN)  ::  REPI(6+N+N/2)
      DOUBLE PRECISION, INTENT(IN)  ::  REPR(4*N+3)
      DOUBLE PRECISION, INTENT(IN)  ::  TAUBAR, GAPTOL
*
      INTEGER,          INTENT(INOUT) :: EWL_AE(2*ICBEG-1:2*ICEND)
      INTEGER,          INTENT(INOUT) :: IWORK( 2*(ICEND-ICBEG+1) )
      DOUBLE PRECISION, INTENT(INOUT) :: LGAP, UGAP
      DOUBLE PRECISION, INTENT(INOUT) :: EWL_LU(2*ICBEG-1:2*ICEND)
*
      INTEGER,          INTENT(OUT) :: NCAND, NBATCH
      INTEGER,          INTENT(OUT) :: ACLOC( MAXNC )
      INTEGER,          INTENT(OUT) :: ABEND( MAXNB )
      DOUBLE PRECISION, INTENT(OUT) :: ACTAU( MAXNC )
*
*  Purpose
*  =======
*
*     Produces shift candidates to try around a cluster.
*
*     For candidate 1 <= i <= NCAND, ACTAU(i) is the shift (absolute),
*     and the location is encoded in ACLOC(i):
*        ABS( ACLOC(i) )  is the eigenvalue we are shifting to,
*        SIGN( ACLOC(i) ) is the direction of the offset/gap, that is,
*                         +1 if shifting to the right, -1 if left.
*     As example, ACLOC(1)=-5 means the first shift goes to the left
*     of eigenvalue 5. This implies that all shifts to this location
*     are <= the lower bound of ew 5.
*       A special kind of candidates are those which are placed blindly
*     and not inside a gap, thos do not have a location but instead
*     ACLOC(i) is set to ICEND+1.
*
*     The shifts are grouped into up to NBATCH <= 2*MAXCPL batches
*                   1, ... , ABEND(1),
*          ABEND(1)+1, ... , ABEND(2)
*                      ...
*          ABEND(NBATCH-1), ... , ABEND(NBATCH),
*     where ABEND(NBATCH)=NCAND
*
*     Out:
*       NCAND=0,NBATCH=0 is possible.
*       Batches contain at least one candidate, and each candidate is
*       contained in exactly one batch (see above).
*
*  Arguments
*  =========
*
*  MAXCPO  (input) INTEGER
*          Maximal number of candidates for both outside locations.
*          Set to 0 to deactivate outside shifts.
*
*  MAXCPI  (input) INTEGER
*          Maximal number of candidates for inside locations.
*          Set to 0 to deactivate inside shifts.
*
*  MAXNC   (input) INTEGER
*          Maximal number of candidates.
*          This is not a parameter, but must be set to
*            MAXNC  :=  2*MAXCPO + 2*(ICEND-ICBEG)*MAXCPI + 1.
*          It is checked at runtime that the correct value was set.
*          The purpose of this argument is just to make declaring the
*          array in the interface and at the caller easier.
*
*  MAXNB   (input) INTEGER
*          MAximal number of batches.
*          This is not a parameter, but must be set to
*            MAXNB  :=  3 * MAX(MAXCPO,MAXCPI)
*          It is checked at runtime that the correct value was set.
*          The purpose of this argument is just to make declaring the
*          array in the interface and at the caller easier.
*
*
*  ======================================================================
*
*     .. Declarations ..
*
      EXTERNAL          DLAMCH
      DOUBLE PRECISION  DLAMCH
*
*     .. Constants ..
*
      DOUBLE PRECISION, PARAMETER  ::  ZERO  = 0.0D0
      DOUBLE PRECISION, PARAMETER  ::  QUART = 0.25D0
      DOUBLE PRECISION, PARAMETER  ::  HALF  = 0.5D0
      DOUBLE PRECISION, PARAMETER  ::  ONE   = 1.0D0
      DOUBLE PRECISION, PARAMETER  ::  FOUR  = 4.0D0
*
*     .. Parameters ..
*
*
*     .. Parameters ..
*


*         Cluster boundaries will be refined to relative accuracy
*         RACLBFAC*N*EPS
      INTEGER, PARAMETER  ::  RACLBFAC = 4



*     Maximal number of inside shifts for primary batches
      INTEGER, PARAMETER  ::  MAXINPRIM = 2


*     For an inside shift to become primary, we require that the
*     smaller group of eigenvalues left/right has at least
*         MAX( 2, MININPRIMPART * (#ews in cluster) )
*     elements.
      DOUBLE PRECISION, PARAMETER  ::  MININPRIMPART = 0.4D0

*     Activate the last resort option to just throw a shift somewhere in the
*     middle of a tight cluster.
      LOGICAL, PARAMETER  ::  DOLASTRESORTTHROWIN = .TRUE.

*
*     .. Locals ..
*
      DOUBLE PRECISION  LB, UB, EWSIGN, EPS, PREC, RACLB
      DOUBLE PRECISION  SRNGLB, SRNGUB, GAPL, GAPU
      DOUBLE PRECISION  MINGAP, WIDTH, AVGAP, INGAP, BOUND, DELTA
      DOUBLE PRECISION  OFF0, OFF, OFFMAX

      INTEGER           MAXNLOC, MAXNBATCH, MAXCPLIN
      INTEGER           ICLEN, NLOCS, PREFER, NPRIM, NSCND, NTRNY, NPART
      INTEGER           I, J, I0, J0, M, KL, KU, ISLOT, IBATCH
      INTEGER           DIR, INDEX, KTRY, ILOC, FLAG, NINGAP
      INTEGER           BSLOT, BBTCH, TRYMAX, NLOCPH, NLOUT
      INTEGER           ILBOUT, ILBIN, NLBOUT, NLBIN, BLOC
      INTEGER           ITMP77
      INTEGER IXF77A, IXF77B, IXF77C

      LOGICAL           IOK, JOK, GOINSIDE, GOOUTSIDE
*
*  === Executable Statements ============================================
*
      ICLEN = ICEND - ICBEG + 1
      MAXNLOC = 2*ICLEN
      EPS = DLAMCH('Epsilon')
      PREC = DLAMCH('Precision')
      LB = EWL_LU(2*ICBEG - 1)
      UB = EWL_LU(2*ICEND)
      EWSIGN = SIGN(ONE,LB)
      RACLB = RACLBFAC * N * EPS


      GOINSIDE  = (MAXCPI .GT. 0)
      GOOUTSIDE = (MAXCPO .GT. 0)

      I0 = -1
      J0 = -1

      SRNGLB = LB - MIN( QUART*LGAP, GAPTOL*ABS(LB) )
      SRNGUB = UB + MIN( QUART*UGAP, GAPTOL*ABS(UB) )

      IF( TAUBAR.GT.ZERO .AND. LB.GT.ZERO )THEN
         SRNGUB = MIN(SRNGUB, TAUBAR + 2*LB )
      ELSEIF( TAUBAR.LT.ZERO .AND. UB.LT.ZERO )THEN
         SRNGLB = MAX(SRNGLB, TAUBAR + 2*UB )
      ENDIF
*
*     ===================================================================
*                         Select Shift Locations
*     ===================================================================
*
*     More precisely, we order them with the most promising first.
*     The locations are stored in IWORK.

*     Set PREFER TO +1 to emphasize shifts in the right half,
*     and to -1 for left.
*     TODO: prefer side "closer to min/max ew" (Parlett). Can
*           detect using taubar and spdiam
*     At the moment: prefer near side
      IF( LB.GT.ZERO )THEN
         PREFER = -1
      ELSE
         PREFER = +1
      ENDIF

      NINGAP = 0
      IF( GOINSIDE )THEN
         M = (ICBEG+ICEND) / 2
         I0 = EWL_AE( 2*(M + MOD(ICBEG+ICEND,2)) - 1 )
         J0 = EWL_AE( 2*M )
*        For gaps <= I0 we normally shift only to left ends of subgroups
*        except as a last resort, similar for J0.

*        Loop once over all inner locations and determine where we want
*        them:
*          1 - eligible even for primary
*          2 - eligible for secondary
*          3 - only last resort
*        We store these flags temporarily in ACLOC(1:2*ILEN), where
*          ACLOC(2*(I-ICBEG) + 1)  is for shift left of I
*          ACLOC(2*(I-ICBEG) + 2)  is for shift right of I.
*        These flags are later set to zero as soon as the location was
*        chosen. The outer ones are always used, we set their flags to
*        0 as stoppers.
         ACLOC(1) = 0
         I = ICBEG
         DO
            DO
               J = EWL_AE(2*I)
               I = J+1
               IF( RGINFO(J) .GT. 0 ) EXIT
            ENDDO
            IF( J .EQ. ICEND )THEN
               EXIT
            ENDIF
*           Inner gap ]J,I[
            NINGAP = NINGAP + 1

            NPART = MIN( I-ICBEG, ICEND-J )
            IF( NPART.LT.2 )THEN
               ACLOC(2*(I-ICBEG) + 1) = 3
               ACLOC(2*(J-ICBEG) + 2) = 3
            ELSE
               IF( DBLE(NPART).GE.MININPRIMPART*ICLEN )THEN
                  ACLOC(2*(I-ICBEG) + 1) = 1
                  ACLOC(2*(J-ICBEG) + 2) = 1
               ELSE
                  ACLOC(2*(I-ICBEG) + 1) = 2
                  ACLOC(2*(J-ICBEG) + 2) = 2
               ENDIF
               IF( I .GT. I0 ) ACLOC(2*(I-ICBEG) + 1) = 3
               IF( J .LT. J0 ) ACLOC(2*(J-ICBEG) + 2) = 3
            ENDIF
         ENDDO
         ACLOC(2*ICLEN) = 0
      ENDIF

      NLOCS = 0

*     First pass: Inside shifts for primary batches.
      IF( GOINSIDE .AND. NINGAP.GT.0 )THEN
         I = I0
         J = J0
         DO
            IF( (I.EQ.ICBEG .AND. J.EQ.ICEND) .OR.
     $           NLOCS.EQ.MAXINPRIM )
     $      THEN
               EXIT
            ENDIF
            DO
               IF( RGINFO(J) .GT. 0 ) EXIT
               J = EWL_AE(2*(J+1))
            ENDDO
            DO
               IF( RGINFO(I-1) .GT. 0 ) EXIT
               I = EWL_AE(2*(I-1)-1)
            ENDDO
*           Now there are gaps left of I and right of J

            IOK = ( ACLOC(2*(I-ICBEG) + 1) .EQ. 1 )
            JOK = ( ACLOC(2*(J-ICBEG) + 2) .EQ. 1 )

            IF( IOK .AND. (.NOT.JOK .OR. PREFER.NE.+1) )THEN
               NLOCS = NLOCS + 1
               IWORK(NLOCS) = -I
               ACLOC(2*(I-ICBEG) + 1) = 0
               I = I-1
            ELSEIF( JOK )THEN
               NLOCS = NLOCS + 1
               IWORK(NLOCS) = +J
               ACLOC(2*(J-ICBEG) + 2) = 0
               J = J+1
            ENDIF
            IF( .NOT.IOK .AND. I.GT.ICBEG )  I = I-1
            IF( .NOT.JOK .AND. J.LT.ICEND )  J = J+1
         ENDDO
      ENDIF

*     Outside locations
      NLOUT = 0
      IF( GOOUTSIDE )THEN
         IF( PREFER.EQ.-1 )THEN
            IWORK(NLOCS+1) = -ICBEG
            IWORK(NLOCS+2) = +ICEND
         ELSE
            IWORK(NLOCS+1) = +ICEND
            IWORK(NLOCS+2) = -ICBEG
         ENDIF
         NLOCS = NLOCS + 2
         NLOUT = NLOUT + 2
      ENDIF

      NPRIM = NLOCS

*     Secondary & Ternary: Remaining inside shifts, first those on the
*     right gap sides, then all others
      NSCND = 0
      NTRNY = 0
      IF( GOINSIDE .AND. NINGAP.GT.0 )THEN
         DO FLAG = 2,3
            I = I0
            J = J0
            DO
               IF( I.EQ.ICBEG .AND. J.EQ.ICEND )THEN
                  EXIT
               ENDIF
               DO
                  IF( RGINFO(J) .GT. 0 ) EXIT
                  J = EWL_AE(2*(J+1))
               ENDDO
               DO
                  IF( RGINFO(I-1) .GT. 0 ) EXIT
                  I = EWL_AE(2*(I-1)-1)
               ENDDO
*              Now there are gaps left of I and right of J

               IOK = ( ACLOC(2*(I-ICBEG) + 1) .EQ. FLAG )
               JOK = ( ACLOC(2*(J-ICBEG) + 2) .EQ. FLAG )

               IF( IOK .AND. (.NOT.JOK .OR. PREFER.NE.+1) )THEN
                  NLOCS = NLOCS + 1
                  IWORK(NLOCS) = -I
                  ACLOC(2*(I-ICBEG) + 1) = 0
                  I = I-1
               ELSEIF( JOK )THEN
                  NLOCS = NLOCS + 1
                  IWORK(NLOCS) = +J
                  ACLOC(2*(J-ICBEG) + 2) = 0
                  J = J+1
               ENDIF
               IF( .NOT.IOK .AND. I.GT.ICBEG )  I = I-1
               IF( .NOT.JOK .AND. J.LT.ICEND )  J = J+1
            ENDDO
            IF( FLAG.EQ.2 )THEN
               NSCND = NLOCS - NPRIM
            ELSE
               NTRNY = NLOCS - NSCND - NPRIM
            ENDIF
         ENDDO
      ENDIF
*
*     ===================================================================
*                         Place the Shift Candidates
*     ===================================================================
*
      DO ITMP77 = 1, MAXNB
         ABEND(ITMP77) = 0
      ENDDO
      DO ITMP77 = 1, MAXNC
         ACLOC(ITMP77) = 0
      ENDDO

*     We fill the candidate data fields ACLOC and ACTAU by scanning the
*     locations once from left to right. This leaves some entries unset
*     (ACLOC(i)=0). The fields will be compressed later.

*     At the same time we count in the entries of ABEND the number of
*     elements per batch.

*     The first slot of candidates within the current phase (primary,
*     secondary, ternary).
      BSLOT = 1
*     Number of the first batch within the current phase
      BBTCH = 1
*     Number of locations in the current phase
      NLOCPH = NPRIM

*     First location in current phase
      BLOC = 1

*     The current phase has NLBOUT outside and NLBIN inside locations.
      NLBOUT = NLOUT
      NLBIN  = NPRIM - NLBOUT

      ILBOUT = 0
      ILBIN  = 0

c      write(*,*) '================================================'
c      write(*,*) 'primary: nlocph = ',NLOCPH

      NCAND = 0
      DO ILOC = 1, NLOCS

         IF( ILOC .EQ. NPRIM+1 )THEN
*           We now move to secondary locations.
            BSLOT = BSLOT + NLOUT*MAXCPO + (NPRIM-NLOUT)*MAXCPI
            BBTCH = BBTCH + MAX( MAXCPO, MAXCPI )
            NLOCPH = NSCND
            BLOC = ILOC
c            write(*,*)
c     $         'secondary: bslot=',bslot,'bbtch=',bbtch,
c     $         'nlocph=',nlocph,'bloc=',bloc
            ILBOUT = 0
            ILBIN  = 0
            NLBOUT = 0
            NLBIN  = NSCND
         ENDIF
*        Note that we might have no secondary or ternary locations
*        at all.
         IF( ILOC .EQ. NPRIM+NSCND+1 )THEN
*           We now move to ternary locations.
            BSLOT = BSLOT + NSCND*MAXCPI
            BBTCH = BBTCH + MAXCPI
            NLOCPH = NTRNY
            BLOC = ILOC
c            write(*,*)
c     $         'ternary: bslot=',bslot,'bbtch=',bbtch,
c     $         'nlocph=',nlocph,'bloc=',bloc
            ILBOUT = 0
            ILBIN  = 0
            NLBOUT = 0
            NLBIN  = NTRNY
         ENDIF

         INDEX = ABS( IWORK(ILOC) )
         DIR   = SIGN(1, IWORK(ILOC))

*        We might be interested in just a chunk of the cluster
         IF( DIR .EQ. +1 )THEN
            KL = ICBEG
            KU = INDEX
            BOUND = EWL_LU(2*INDEX)
         ELSE
            KL = INDEX
            KU = ICEND
            BOUND = EWL_LU(2*INDEX - 1)
         ENDIF

*        Determine minimal absolute gap separating the chunk
         IF( KL .EQ. ICBEG )THEN
            GAPL = LGAP
         ELSE
            GAPL = EWL_LU(2*KL-1) - EWL_LU(2*KL-2)
         ENDIF
         IF( KU .EQ. ICEND )THEN
            GAPU = UGAP
         ELSE
            GAPU = EWL_LU(2*KL+1) - EWL_LU(2*KL)
         ENDIF
         MINGAP = MIN( GAPL, GAPU )

*        See if we have an internal gap inside the chunk, on the
*        other side of the outermost interval
         INGAP = ZERO
         IF( DIR .EQ. +1 )THEN
            I = EWL_AE(2*INDEX - 1)
            IF( I .GT. KL )THEN
               INGAP = EWL_LU(2*I-1) - EWL_LU(2*I-2)
            ENDIF
         ELSE
            J = EWL_AE(2*INDEX)
            IF( J .LT. KU )THEN
               INGAP = EWL_LU(2*J+1) - EWL_LU(2*J)
            ENDIF
         ENDIF

         WIDTH = EWL_LU(2*KU) - EWL_LU(2*KL-1)
         AVGAP = ZERO
         IF( KL .NE. KU )  AVGAP = WIDTH / (KU - KL)

*        Set maximal allowed offset
         IF( DIR .EQ. +1 )THEN
            OFFMAX = SRNGUB - BOUND
         ELSE
            OFFMAX = BOUND - SRNGLB
         ENDIF
         OFFMAX = MIN( OFFMAX, QUART*MINGAP, QUART*ABS(BOUND) )
*        now offmax is negative if bound is outside SRNG

*        Set maximal number of tries at this location, depending on
*        if we are inside or not.
         IF( IWORK(ILOC).EQ.-ICBEG .OR. IWORK(ILOC).EQ.+ICEND )THEN
            ILBOUT = ILBOUT + 1
            TRYMAX = MAXCPO
         ELSE
            ILBIN  = ILBIN + 1
            TRYMAX = MAXCPI
         ENDIF
*        The current location is the ILBOUT'th outer location and the
*        ILBIN'th inner location within the current phase.

*        Place the shifts at this location
         KTRY = 0
         OFF0 = ABS(BOUND)*(8*EPS)
         DELTA = MAX( OFF0, MAX( AVGAP, INGAP ) / 2**TRYMAX )
         IF( OFF0 .GT. OFFMAX )  OFF0 = ZERO
         DO
            OFF = OFF0 + (2**KTRY - 1)*DELTA
            IF( KTRY.GE.TRYMAX .OR. OFF.GT.OFFMAX )THEN
               EXIT
            ENDIF
            IBATCH = BBTCH + KTRY
C            ISLOT  = BSLOT + KTRY*NLOCPH + ILOC-BLOC

            ISLOT = BSLOT-1
     $            + MIN(KTRY,MAXCPO)*NLBOUT
     $            + MIN(KTRY,MAXCPI)*NLBIN
            IF( KTRY+1 .LE. MAXCPO )  ISLOT = ISLOT + ILBOUT
            IF( KTRY+1 .LE. MAXCPI )  ISLOT = ISLOT + ILBIN


            ACTAU(ISLOT) = BOUND + DIR*OFF
            ACLOC(ISLOT) = IWORK(ILOC)
c            write(*,*) 'batch',IBATCH,'slot',ISLOT,'loc',IWORK(ILOC)
C     $        'bound',bound,'dir',dir,'off',off,'tau',actau(islot)
            KTRY  = KTRY + 1
            NCAND = NCAND + 1
            ABEND(IBATCH) = ABEND(IBATCH) + 1
         ENDDO
      ENDDO
*
*     .. Decode Batch Boundaries ..
*
      NBATCH = 0
      J = 0
      IBATCH = 1
      DO IBATCH = 1, MAXNB
         M = ABEND(IBATCH)
*        Now M is the number of candidates in batch IBATCH
         IF( M .GT. 0 )THEN
            NBATCH = NBATCH + 1
            J = J + M
            ABEND(NBATCH) = J
         ENDIF
      ENDDO
      DO IXF77A = NBATCH+1, MAXNB
         ABEND(IXF77A) = 0
      ENDDO

*
*     Special Backup: If there are no inside gaps (tight cluster)
*     we just throw a last resort shift somewhere inside in the
*     ternary batch.
*     These shifts could be regarded as outer ones, since they do
*     only depend on the outer bounds of the eigenvalues. However,
*     we think they should only be used as last resort, so they are
*     only placed if no inside shifts would be possible and in the last batch
      IF( DOLASTRESORTTHROWIN .AND. NINGAP.EQ.0 .AND. NCAND.LT.MAXNC
     $    .AND.
     $    (UB-LB) .LE. SQRT(PREC)*MAX(ABS(LB),ABS(UB))
     $)THEN
         NCAND = NCAND + 1
*        mark as special
         ACLOC(NCAND)  = ICEND+1
*        just take the midpoint
         ACTAU(NCAND)  = HALF * (LB + UB)
         ABEND(NBATCH) = ABEND(NBATCH) + 1
      ENDIF

*
*     .. Compress the candidate fields ..
*
      I = 1
      DO
         IF( I .GT. NCAND ) EXIT
         IF( ACLOC(I) .EQ. 0 )THEN
            J = I+1
            DO
               IF( ACLOC(J) .NE. 0 )  EXIT
               J = J+1
            ENDDO
            ACLOC(I) = ACLOC(J)
            ACTAU(I) = ACTAU(J)
            ACLOC(J) = 0
         ENDIF
         I = I+1
      ENDDO


      END SUBROUTINE DLAXRF_SELSHF
*
************************************************************************
