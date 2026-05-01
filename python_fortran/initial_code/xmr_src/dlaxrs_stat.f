      SUBROUTINE DLAXRS_STAT(
     $             N, J, I0, I1, G, OMEGA, GNSQ, NGN, BDET,
     $             PIVMIN, TAU, SHFPRT, GPLUS, OPLUS, S
     $           )
      IMPLICIT NONE
*
      INTEGER,          INTENT(IN)  ::  N, J, I0, I1
      INTEGER,          INTENT(IN)  ::  OMEGA(0:N+1)
      DOUBLE PRECISION, INTENT(IN)  ::  PIVMIN, TAU
      DOUBLE PRECISION, INTENT(IN)  ::  G(N), GNSQ(N), NGN(N), BDET(N)
      DOUBLE PRECISION, INTENT(IN)  ::  SHFPRT(N)
*
      INTEGER,          INTENT(OUT)  ::  OPLUS(N)
      DOUBLE PRECISION, INTENT(OUT)  ::  GPLUS(I0:I1), S(N)
*
*  Purpose
*  =======
*
*     Computes stationary factorization from 1 to J or from N down to J.
*
*     I0=1 and I1=N-1
*       assumes that the twist of the source is >= J
*       produces
*         GPLUS(1:J-1), OPLUS(1:J) and S(1),...,S(J),
*         only S(1)=0, OPLUS(1)=0 if J=1
*
*     I0=2 and I1=N
*       assumes that the twist of the source is <= J
*       produces
*         GPLUS(J+1:N), OPLUS(J:N) and S(J),...,S(N),
*         only S(N)=0, OPLUS(J)=0 if J=N
*
*     Note that the decision about whether a block ends at J has to be
*     made here. But since J must assumed to be a twist index, the
*     routine will allow OPLUS(J)=1 only if J in {1,N}.
*
*  ======================================================================
*
Cextract -b parameters.inc.f procname=laxrs_stat
*
*     KBLOCK --  Base constant for block-determinant criterion
*
      DOUBLE PRECISION KBLOCK
      PARAMETER (KBLOCK = 1.0D0 / 8.0D0)
*
*     KBLMOD --  Adjusted KBLOCK for practical tests
*
      DOUBLE PRECISION KBLMOD
      PARAMETER (KBLMOD = KBLOCK * 0.999D0)
*
*     PK1, PK2
*
      DOUBLE PRECISION PK1, PK2
      PARAMETER (PK1 = KBLMOD, PK2 = KBLMOD / 3.01D0)
*
*     PRBRK, PRCRE, PROSQ
*
      DOUBLE PRECISION PRBRK, PRCRE, PROSQ
      PARAMETER (PRBRK = 5.0D0, PRCRE = 5.0D0, PROSQ = 0.25D0)
*
*     PMAXOSEQLEN,  >= 0
*       Set to 0 to only allow clean creation or those that are
*       guaranteed to be stable (small d, block at end).
      INTEGER, PARAMETER  ::  PMAXOSEQLEN = 1

*     USEBLOCKS
*       Set to false to deactivate creation of blocks completely.
      LOGICAL, PARAMETER  ::  USEBLOCKS = .TRUE.
*
*     .. Declarations ..
*
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

*
*     .. Locals ..
*
      DOUBLE PRECISION  ABSTAU, OLDSMT, SMT, AUX, CPLUS, BDETP, X
      INTEGER           I, DIR, IPREV, INEXT, BCOUNT
      LOGICAL           BRES
*
*     -- Executable Statements -----------------------------------------
*
      ABSTAU = ABS(TAU)

      IF( I0.EQ.1 .AND. I1.EQ.N-1 )THEN
         DIR = +1
         I   = 1
      ELSEIF( I0.EQ.2 .AND. I1.EQ.N )THEN
         DIR = -1
         I   = N
      ELSE
         DIR = 0
         I   = -1
      ENDIF
      INEXT = I+DIR
*
*     Main Loop
*
      AUX  = ZERO
      S(I) = ZERO
      OPLUS(I) = 0
      DO
         IF( I.EQ.J )THEN
            EXIT
         ENDIF

         SMT = AUX - TAU*SHFPRT(I)
         GPLUS(I) = G(I) + SMT
         IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = ZERO
         OPLUS(INEXT) = 0

         IF( OMEGA(INEXT).NE.0 )THEN
*           Check if the block should be kept
            CPLUS = G(INEXT) - TAU*SHFPRT(INEXT)
            BRES = ( ABS(GPLUS(I)*CPLUS) .LT. PK1*GNSQ(I) )
            IF( BRES .AND. (INEXT.NE.J .OR. J.EQ.N .OR. J.EQ.1) )THEN
*
*              -- Keep a block --
*
               AUX = ZERO
               OPLUS(INEXT) = 1
               IF( INEXT.NE.J )THEN
                  S(INEXT) = ZERO
                  GPLUS(INEXT) = CPLUS
*
                  BDETP = GPLUS(I)*CPLUS - GNSQ(I)
                  IF( G(I).EQ.ZERO )THEN
                     AUX = - (GNSQ(INEXT) / BDETP) * GPLUS(I)
                  ELSE
*                    recall that NGN(i) is not defined YET
C                    X = NGN(I)*SMT - GPLUS(I)*(TAU*SHFPRT(INEXT))
                     X = (GNSQ(I) / G(I)) * SMT
     $                   - GPLUS(I)*(TAU*SHFPRT(INEXT))
                     AUX = NGN(INEXT) * (X / BDETP)
                  ENDIF
                  I     = INEXT
                  INEXT = INEXT+DIR
                  OPLUS(INEXT) = 0
               ENDIF
            ELSE
*
*              -- Initiate breaking the block --
*
               AUX = - GNSQ(I) / GPLUS(I)
            ENDIF
         ELSE
*           check if a new block should be created
*           Have to avoid creating a block that would end at an inner
*           twist in the source.
            IF( .NOT.USEBLOCKS .OR.
     $          (INEXT.EQ.J .AND. J.NE.1 .AND. J.NE.N) )
     $      THEN
               BRES = .FALSE.
            ELSEIF( INEXT.EQ.J .OR. OMEGA(INEXT+DIR).EQ.0 )THEN
*              laxer criterion, includes creating block at end
               CPLUS = G(INEXT) + (NGN(I) - TAU*SHFPRT(INEXT))
               BRES = ( ABS(GPLUS(I)*CPLUS) .LT. PK1*GNSQ(I) )

            ELSE
               BRES = ( GPLUS(I).EQ.ZERO .OR. PMAXOSEQLEN.GT.0 )
     $          .AND. ( (ABS(GPLUS(I))
     $                    * MAX(ABS(TAU*SHFPRT(INEXT)),ABS(G(INEXT))))
     $                  .LT. PK2*GNSQ(I) )
     $          .AND. (ABS(GPLUS(I)) .LT. PK2*ABS(G(I)))
            ENDIF
*
            IF( BRES )THEN
*
*              -- Initiate creation of a block --
*
               OPLUS(INEXT) = 1
               AUX = NGN(I)
            ELSE
*
*              -- Standard dstqds --
*               Note: may have to avoid breakdown, eg if useblocks=false
               IF( GPLUS(I).EQ.ZERO )  GPLUS(I) = -PIVMIN
               AUX = NGN(I) * (SMT / GPLUS(I))
            ENDIF
         ENDIF
         OLDSMT = SMT
         S(INEXT) = AUX

         IPREV = I
         I = INEXT
         INEXT = INEXT+DIR

*        =======================================================
*        ==  Block structure change and overlap-control loop  ==
*        =======================================================
         BCOUNT = 0
         IF( OPLUS(I).NE.0 .AND. OMEGA(INEXT).NE.0 )  BCOUNT = 1

         DO
            IF( I.EQ.J .OR. ((OMEGA(I).EQ.0).EQV.(OPLUS(I).EQ.0)) )THEN
               EXIT
            ENDIF
            SMT = AUX - TAU*SHFPRT(I)
            GPLUS(I) = G(I) + SMT
            IF( ABS(GPLUS(I)).LT.PIVMIN )  GPLUS(I) = ZERO
            OPLUS(INEXT) = 0

*           reset overlap-counter if perturbations can be attributed
*           to the shift
            IF( ABS(AUX) .LT. PROSQ*ABSTAU )THEN
               BCOUNT = 0
            ENDIF

            IF( OPLUS(I).NE.0 )THEN
               CPLUS = GPLUS(I)
               BDETP = GPLUS(IPREV)*CPLUS - GNSQ(IPREV)
               IF( OMEGA(INEXT).NE.0 )THEN
*                 = - ngnp(i)
                  AUX = - (GNSQ(I) / BDETP) * GPLUS(IPREV)
               ELSE
*
*                 -- End by create or clean create--
*
                  IF( GPLUS(IPREV).EQ.ZERO )THEN
                     AUX = NGN(I)
                  ELSE
                     IF( OMEGA(IPREV).NE.0 .OR.
     $                   ABS(AUX).LE.PRCRE*ABSTAU .OR.
     $                   SIGN(ONE,GPLUS(IPREV)).NE.SIGN(ONE,G(IPREV)) )
     $               THEN
                        X = GPLUS(IPREV)*SMT - GNSQ(IPREV)
                     ELSE
                        X = AUX*OLDSMT - GPLUS(IPREV)*TAU*SHFPRT(I)
                     ENDIF
                     AUX = NGN(I) * (X / BDETP)
                  ENDIF
               ENDIF
            ELSE
*              check if a new block (with overlap) may be created
               IF( .NOT.USEBLOCKS .OR.
     $             (INEXT.EQ.J .AND. J.NE.1 .AND. J.NE.N) )
     $         THEN
                  BRES = .FALSE.
               ELSE
                  BRES = ( INEXT.EQ.J .OR. BCOUNT.LT.PMAXOSEQLEN .OR.
     $                     ABS(GPLUS(I))*GNSQ(INEXT)
     $                      .LT. (ONE - KBLMOD) * (PROSQ*ABSTAU) *
     $                           GNSQ(I) )
                  BRES = BRES .AND.
     $                   ( ABS(GPLUS(I)) *
     $                     MAX(ABS(TAU*SHFPRT(INEXT)),ABS(G(INEXT)))
     $                     .LT. PK2*GNSQ(I) ) .AND.
     $                   ( ABS(GPLUS(I)*G(IPREV))
     $                     .LT. PK2*ABS(BDET(IPREV)) )
               ENDIF
               IF( BRES )THEN
*
*                 -- Create block, continuing the sequence --
*
                  OPLUS(INEXT) = 1
                  AUX = NGN(I)
                  BCOUNT = BCOUNT+1
               ELSE
*
*                 -- End by break or clean break --
*
                  IF( GPLUS(I).EQ.ZERO )  GPLUS(I) = -PIVMIN
*
                  IF( G(IPREV).EQ.ZERO )THEN
                     AUX = - GNSQ(I) / GPLUS(I)
                  ELSE
                     IF( OPLUS(IPREV).NE.0 .OR.
     $                   ABS(AUX).LE.PRBRK*ABS(TAU*SHFPRT(I)) .OR.
     $                   SIGN(ONE,G(IPREV)).NE.SIGN(ONE,GPLUS(IPREV)) )
     $               THEN
                        X = G(IPREV)*SMT + GNSQ(IPREV)
                     ELSE
                        X = -AUX*OLDSMT - G(IPREV)*TAU*SHFPRT(I)
                     ENDIF
                     AUX = (GNSQ(I) * X) / (BDET(IPREV) * GPLUS(I))
                  ENDIF
               ENDIF
            ENDIF
            S(INEXT) = AUX
            OLDSMT   = SMT
            IPREV = I
            I     = INEXT
            INEXT = INEXT + DIR
         ENDDO
      ENDDO
      END SUBROUTINE DLAXRS_STAT
*
************************************************************************
