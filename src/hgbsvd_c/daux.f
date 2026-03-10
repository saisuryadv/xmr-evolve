c#define VERBOSE
c#define INFO_VERBOSE
c#define DT
c#define VS
c#define PR
c#define HL
c#define EVAL
c#define HLVAL
c#define HLVEC
c#define BLAS_BASED
c#define NaN_Check
c#define Assign_Check
c
c       ftncheck:                               okay
c       1111:                                   okay
c       invalid input token:                    okay
c       Praeprozessor:                          BLAS_BASED,NaN_Check,Assign_Check
c       Kommentare:                             okay
c       Semantik:                               okay
c       Liste subroutinen, ueberfluesige:       okay
c       Auskommentierte Anweisungen:            okay
c       NaN-Proof:                              okay
c       workarounds:                            okay
c       optimiert:                              okay
c       info:                                   okay
c       init:                                   okay
c       Dokumentation purpose / parameter
c
*******************************************************************************
*
*       daux.f : auxiliary-routines for dtest.f
*
*       This file contains a set of useful tools for the driver dtest :
*
*       Date:   Sep 10, 1999
*       Author: Benedikt Grosser
*
*       dcheck_by_max: ( show )
*       dcheck_by_norm:         ( nicht benutzt )
*       dminiorthcheck:         ( nicht benutzt )
*       dmview:                 ( nicht benutzt )
*       dvview:                 ( nicht benutzt )
*       ivview:                 ( nicht benutzt )
*       dlv:                    ( nicht benutzt )
*       dlv1:                   ( nicht benutzt )
*       dlv2:                   ( nicht benutzt )
*       dlv2rd:                 ( nicht benutzt )
*
*******************************************************************************
*
        SUBROUTINE DCHECK_BY_MAX(
     >    M,N,A,SIGMA,DDIM,B,EDIM,
     >    U,MMAX,V,NMAX,
     >    UPDU,UPDV,TRANSU,TRANSV,
     >    RES,ORTHU,ORTHV,
     >    MAXU,MAXV,MAXR,SIGMX,
     >    MAXUC,MAXVC,MAXRC,
     >    WORK,LWORK,WORK1,LWORK1,INFO)
*
        IMPLICIT NONE
*
*       purpose
*
*       Check accuracy:
*       Deviation from orthogonality and residual.
*
*       orthu = norm(UT*U-I,1)  maxu = max(max(UT*U-I)) maxuc = maxu/(m*eps)
*       orthv = norm(VT*V-I,1)  maxu = max(max(VT*V-I)) maxvc = maxv/(n*eps)
*       res = norm(U*Sigma*VT-B,1) or norm(UT*B*V-Sigma,1)
*       maxr = max(max(U*Sigma*VT-B)) or max(max(UT*B*V-Sigma))
*       maxrc = maxc/(n*eps*sigmx)
*
*       ..
*       .. parameters ..
*       ..
        INTEGER M,N,DDIM,EDIM,MMAX,NMAX,LWORK,LWORK1,INFO
        LOGICAL UPDU,UPDV,TRANSU,TRANSV
        DOUBLE PRECISION RES,ORTHU,ORTHV
        DOUBLE PRECISION MAXU,MAXV,MAXR,SIGMX
        DOUBLE PRECISION MAXUC,MAXVC,MAXRC
        DOUBLE PRECISION SIGMA(DDIM),A(DDIM),B(EDIM)
        DOUBLE PRECISION U(MMAX,MMAX), V(NMAX,NMAX)
        DOUBLE PRECISION WORK(LWORK),WORK1(LWORK1)
*
        DOUBLE PRECISION ZERO, ONE
        PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
*       ..
*       .. external functions ..
*       ..
        EXTERNAL DLANGE, DLAMCH, DDOT
        DOUBLE PRECISION DLANGE, DLAMCH, DDOT
*       ..
*       .. local variables ..
*       ..
        DOUBLE PRECISION DUMMY, MACHEPS, DN, TMP, TMP1, FACT
        INTEGER I,J,WIND
        INTEGER K
*       ..
*       .. executable statements ..
*       ..
        INFO = 0
        DUMMY = ZERO
        MACHEPS = DLAMCH('P')
        DN=DBLE(N)
        FACT = DN*MACHEPS
*
*       dummy statments to make ftnchek happy
*
        TMP1 = ZERO
        IF (TRANSU) TMP1 = ONE
        IF (TRANSV) TMP1 = ONE
        DUMMY = TMP1
*
*       determine orthu
*
        IF (.NOT.UPDU) THEN
          ORTHU = -ONE
        ELSE
          CALL DLASET('F',M,M,ZERO,ONE,WORK,M)
          CALL DGEMM('T','N',M,M,M,ONE,U,MMAX,U,MMAX,-ONE,WORK,M)
          ORTHU = DLANGE('1',M,M,WORK,M,DUMMY)
          MAXU = ZERO
          MAXUC = ZERO
          DO I = 1,M
            DO J = 1,M
              WIND = (I-1)*M+J
              TMP = ABS(WORK(WIND))
              IF (TMP.GT.MAXU) THEN
                MAXU = TMP
              ENDIF
              IF (TMP.GT.FACT) THEN
                WORK(WIND) = ABS(WORK(WIND))/FACT
                MAXUC = MAX(MAXU,ABS(WORK(WIND)))
              ELSE
                WORK(WIND) = ZERO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*
*       determine orthv
*
        IF (.NOT.UPDV) THEN
          ORTHV = -ONE
        ELSE
          CALL DLASET('F',N,N,ZERO,ONE,WORK,N)
          CALL DGEMM('T','N',N,N,N,ONE,V,NMAX,V,NMAX,-ONE,WORK,N)
          ORTHV = DLANGE('1',N,N,WORK,N,DUMMY)
          MAXV = ZERO
          MAXVC = ZERO
          DO I = 1,N
            DO J = 1,N
              WIND = (I-1)*N+J
              TMP = ABS(WORK(WIND))
              IF (TMP.GT.MAXV) THEN
                MAXV = TMP
              ENDIF
              IF (TMP.GT.FACT) THEN
                WORK(WIND) = TMP/FACT
                MAXVC = MAX(MAXV,ABS(WORK(WIND)))
              ELSE
                WORK(WIND) = ZERO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*
*       determine res
*
        IF (.NOT.(UPDU.AND.UPDV)) THEN
          RES = -ONE
        ELSE
        DO K = 1,N
          DO I = 1,N
*
*           compute B*V(:,i)
*
            WIND = (I-1)*N
            DO J = 1,N-1
              WORK(J+WIND) = A(J)*V(J,I)+B(J)*V(J+1,I)
            ENDDO
            WORK(N+WIND) = A(N)*V(N,I)
            WORK1(K+WIND) = DDOT(N,U(1,K),1,WORK(1+WIND),1)
          ENDDO
        ENDDO
        DO K = 1,N
          WORK1(K+(K-1)*N) = WORK1(K+(K-1)*N) - SIGMA(K)
        ENDDO
          RES = DLANGE('1',M,N,WORK1,M,DUMMY)
          MAXR = ZERO
          MAXRC = ZERO
          DO I = 1,N
            DO J = 1,M
              WIND = (I-1)*M+J
              IF (ABS(WORK1(WIND)).GT.MAXR) THEN
                MAXR = ABS(WORK1(WIND))
              ENDIF
              IF (ABS(WORK1(WIND)).GT.FACT*SIGMX) THEN
                WORK1(WIND) = ABS(WORK1(WIND))/(FACT*SIGMX)
                MAXRC = MAX(MAXR,ABS(WORK1(WIND)))
                WORK1(WIND) = ABS(WORK1(WIND))*(FACT*SIGMX)
              ELSE
                WORK1(WIND) = ZERO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*
        RETURN
        END
*
*******************************************************************************
*
        SUBROUTINE DCHECK_RESULTS(
     >    M,N,A,SIGMA,DDIM,B,EDIM,
     >    U,MMAX,V,NMAX,
     >    UPDU,UPDV,TRANSU,TRANSV,
     >    RES,ORTHU,ORTHV,
     >    MAXU,MAXV,MAXR,SIGMX,
     >    MAXUC,MAXVC,MAXRC,
     >    WORK,LWORK,WORK1,LWORK1,INFO)
*
        IMPLICIT NONE
*
*       purpose
*
*       Check accuracy:
*       Deviation from orthogonality and residual.
*
*       orthu = norm(UT*U-I,1)  maxu = max(max(UT*U-I)) maxuc = maxu/(m*eps)
*       orthv = norm(VT*V-I,1)  maxu = max(max(VT*V-I)) maxvc = maxv/(n*eps)
*       res = norm(U*Sigma*VT-B,1) or norm(UT*B*V-Sigma,1)
*       maxr = max(max(U*Sigma*VT-B)) or max(max(UT*B*V-Sigma))
*       maxrc = maxc/(n*eps*sigmx)
*
*       ..
*       .. parameters ..
*       ..
        INTEGER M,N,DDIM,EDIM,MMAX,NMAX,LWORK,LWORK1,INFO
        LOGICAL UPDU,UPDV,TRANSU,TRANSV
        DOUBLE PRECISION RES,ORTHU,ORTHV
        DOUBLE PRECISION MAXU,MAXV,MAXR,SIGMX
        DOUBLE PRECISION MAXUC,MAXVC,MAXRC
        DOUBLE PRECISION SIGMA(DDIM),A(DDIM),B(EDIM)
        DOUBLE PRECISION U(MMAX,MMAX), V(NMAX,NMAX)
        DOUBLE PRECISION WORK(LWORK),WORK1(LWORK1)
*
        DOUBLE PRECISION ZERO, ONE
        PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
*       ..
*       .. external functions ..
*       ..
        EXTERNAL DLANGE, DLAMCH, DDOT, DNRM2
        DOUBLE PRECISION DLANGE, DLAMCH, DDOT, DNRM2
*       ..
*       .. local variables ..
*       ..
        DOUBLE PRECISION DUMMY, MACHEPS, DN, TMP, TMP1, FACT
        INTEGER I,J,WIND
        INTEGER K
*       ..
*       .. executable statements ..
*       ..
        INFO = 0
        DUMMY = ZERO
        MACHEPS = DLAMCH('P')
        DN=DBLE(N)
        FACT = DN*MACHEPS
*
*       dummy statments to make ftnchek happy
*
        TMP1 = ZERO
        IF (TRANSU) TMP1 = ONE
        IF (TRANSV) TMP1 = ONE
        DUMMY = TMP1
*
*       determine orthu
*
        IF (.NOT.UPDU) THEN
          ORTHU = -ONE
        ELSE
          CALL DLASET('F',M,M,ZERO,ONE,WORK,M)
          CALL DGEMM('T','N',M,M,M,ONE,U,MMAX,U,MMAX,-ONE,WORK,M)
          ORTHU = DLANGE('1',M,M,WORK,M,DUMMY)
          MAXU = ZERO
          MAXUC = ZERO
          DO I = 1,M
            DO J = 1,M
              WIND = (I-1)*M+J
              TMP = ABS(WORK(WIND))
              IF (TMP.GT.MAXU) THEN
                MAXU = TMP
              ENDIF
              IF (TMP.GT.FACT) THEN
	write(*,*) 'orthu ',i,j,tmp
                WORK(WIND) = ABS(WORK(WIND))/FACT
                MAXUC = MAX(MAXU,ABS(WORK(WIND)))
              ELSE
                WORK(WIND) = ZERO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*
*       determine orthv
*
        IF (.NOT.UPDV) THEN
          ORTHV = -ONE
        ELSE
          CALL DLASET('F',N,N,ZERO,ONE,WORK,N)
          CALL DGEMM('T','N',N,N,N,ONE,V,NMAX,V,NMAX,-ONE,WORK,N)
          ORTHV = DLANGE('1',N,N,WORK,N,DUMMY)
          MAXV = ZERO
          MAXVC = ZERO
          DO I = 1,N
            DO J = 1,N
              WIND = (I-1)*N+J
              TMP = ABS(WORK(WIND))
              IF (TMP.GT.MAXV) THEN
                MAXV = TMP
              ENDIF
              IF (TMP.GT.FACT) THEN
c	write(*,*) 'orthv ',i,j,tmp
                WORK(WIND) = TMP/FACT
                MAXVC = MAX(MAXV,ABS(WORK(WIND)))
              ELSE
                WORK(WIND) = ZERO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*
*       determine res
*
        IF (.NOT.(UPDU.AND.UPDV)) THEN
          RES = -ONE
        ELSE
          RES = ZERO
          MAXR = ZERO
          DO I = 1,N
*
*           compute B*V(:,i)
*
            DO J = 1,N-1
              WORK(J) = A(J)*V(J,I)+B(J)*V(J+1,I)
            ENDDO
            WORK(N) = A(N)*V(N,I)
c	do k = 1,n
c	  write(*,*) k,work(k),u(k,i)
c	enddo
c	read(*,*)
	if (1.eq.0) then
	do j = 1,n
	    res = zero
	    do k = 1,n
	      res = res + work(k)*U(k,j)
	    enddo
	if (i.eq.j) res = res - sigma(i)
	if (abs(res).gt.fact*sigma(n)) then
	write(*,*) 'res ',i,j,res,sigma(i)
	endif
	enddo
	else
            DO J = 1,N
              WORK(J) =  WORK(J) - SIGMA(I)*U(J,I)
              MAXR = MAX( MAXR, ABS(WORK(J)) )
c	if (ABS(WORK(J)).gt.fact) then
c	  write(*,*) 'res ',i,j,ABS(WORK(J))
c	endif
            ENDDO
	endif
            WORK1(K) = DNRM2(N,WORK,1)
            RES = MAX(RES,ABS(WORK1(K)))
c       write(*,*) i,res,abs(work1(k))
          ENDDO
        ENDIF
*
        RETURN
        END
*
*******************************************************************************
*
        SUBROUTINE DMVIEW( A,M,N,LLD,MB,NB,NAME )
*
        IMPLICIT NONE
*
        DOUBLE PRECISION A( * )
        CHARACTER*5 NAME
        INTEGER M,N,LLD,MB,NB, I,J,NUMCOL
*
*
*******
*
        WRITE(*,*)
        WRITE(*,*) '+++++ Matrix ',NAME,' ++++++++++++++++++++++++'
        WRITE(*,*)
        WRITE(*,*) '=================================='
        DO I = 1,N,NB
          NUMCOL = NB
          IF ( (I+NB) .GT. N ) NUMCOL = N-I+1
          WRITE(*,*)
          WRITE(*,*) ' Spalte ',I,' bis ', I+NUMCOL-1
          WRITE(*,*) '=================================='
          DO J = 1,M
            IF ( NUMCOL .EQ. 1 )
     >      THEN
              WRITE(*,101) A( LLD*(I-1) + J )
            ENDIF
            IF ( NUMCOL .EQ. 2 )
     >      THEN
              WRITE(*,102) A( LLD*(I-1) + J ), A( LLD*I + J )
            ENDIF
            IF ( NUMCOL .EQ. 3 )
     >      THEN
              WRITE(*,103) A( LLD*(I-1) + J ), A( LLD*I + J ),
     >                     A( LLD*(I+1) + J )
            ENDIF
            IF ( NUMCOL .EQ. 4 )
     >      THEN
              WRITE(*,104) A( LLD*(I-1) + J ), A( LLD*I + J ),
     >                     A( LLD*(I+1) + J ), A( LLD*(I+2) + J )
            ENDIF
            IF ( NUMCOL .EQ. 5 )
     >      THEN
              WRITE(*,105) A( LLD*(I-1) + J ), A( LLD*I + J ),
     >                     A( LLD*(I+1) + J ), A( LLD*(I+2) + J ),
     >                     A( LLD*(I+3) + J )
            ENDIF
            IF ( NUMCOL .EQ. 6 )
     >      THEN
              WRITE(*,106) A( LLD*(I-1) + J ), A( LLD*I + J ),
     >                     A( LLD*(I+1) + J ), A( LLD*(I+2) + J ),
     >                     A( LLD*(I+3) + J ), A( LLD*(I+4) + J )
            ENDIF
            IF ( MOD( J,MB ) .EQ. 0 ) WRITE(*,*)
          ENDDO
        ENDDO
        WRITE(*,*)
        WRITE(*,*) '+++++ Ende ',NAME,' +++++++++++++++++++++++++'
        WRITE(*,*)
 101    FORMAT( D13.4 )
 102    FORMAT( D13.4 ,D13.4 )
 103    FORMAT( D13.4 ,D13.4 ,D13.4 )
 104    FORMAT( D13.4 ,D13.4 ,D13.4 ,D13.4 )
 105    FORMAT( D13.4 ,D13.4 ,D13.4 ,D13.4 ,D13.4 )
 106    FORMAT( D13.4 ,D13.4 ,D13.4 ,D13.4 ,D13.4 ,D13.4 )
        RETURN
        END
*
*******************************************************************************
