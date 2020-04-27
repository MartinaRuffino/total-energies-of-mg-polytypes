      SUBROUTINE LEQR2F(A,N,NDIM,IPIV,WORK)
C  GAUSS FACT. OF (N X N) MATRIX A. IPIV STORES PIVOTS FOR LEQR2S
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION A(NDIM,N),WORK(N),IPIV(N)
      DO 10 IR=1,N
      DO 11 IT=IR,N
        SUM=-A(IT,IR)
        DO 12 I=1,IR-1
  12    SUM=SUM+A(IT,I)*A(I,IR)
  11    WORK(IT)=-SUM
      IRMAX=IR
      R=WORK(IR)*WORK(IR)
      DO 13 IT=IR+1,N
        RR=WORK(IT)*WORK(IT)
        IF(R.GE.RR) GOTO 13
          IRMAX=IT
          R=RR
  13    CONTINUE
      DO 14 I=1,N
        CC=A(IR,I)
        A(IR,I)=A(IRMAX,I)
  14    A(IRMAX,I)=CC
      A(IR,IR)=WORK(IRMAX)
      WORK(IRMAX)=WORK(IR)
      IPIV(IR)=IRMAX
      DO 15 IT=IR+1,N
  15  A(IT,IR)=(1.D0/A(IR,IR))*WORK(IT)
      DO 16 IT=IR+1,N
        SUM=A(IR,IT)
        DO 17 I=1,IR-1
  17    SUM=SUM-A(IR,I)*A(I,IT)
  16    A(IR,IT)=SUM
  10  CONTINUE
      RETURN
      END
      SUBROUTINE LEQR2S(A,N,NDIM,B,NB,IPIV)
C  BACK SUBSTITUTION AFTER LEQR2F.
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION A(NDIM,N),B(NDIM,NB),IPIV(N)
      DO 20 IB=1,NB
      DO 10 IR=1,N
      IRMAX=IPIV(IR)
      SUM=B(IRMAX,IB)
      B(IRMAX,IB)=B(IR,IB)
      DO 11 I=1,IR-1
  11  SUM=SUM-A(IR,I)*B(I,IB)
  10  B(IR,IB)=SUM
      DO 12 IR=N,1,-1
      SUM=B(IR,IB)
      DO 13 I=IR+1,N
  13  SUM=SUM-A(IR,I)*B(I,IB)
  12  B(IR,IB)=SUM/A(IR,IR)
  20  CONTINUE
      RETURN
      END
