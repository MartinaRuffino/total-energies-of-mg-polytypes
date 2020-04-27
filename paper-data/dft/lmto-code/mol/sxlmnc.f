      SUBROUTINE SXLMNC(C,LMX)
C  NORMALIZATION CONSTANTS FOR THE XLM. USE WITH SXLM.
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION C(1)
      TPI=8.D0*DATAN(1.D0)
      FPI=2.D0*TPI
      Y0=1.D0/DSQRT(FPI)
      C(1)=Y0
      DO 2 L=1,LMX
      LP1=L+1
      TLP1=L+LP1
      LM0=(L*LP1)/2+1
      C(LM0)=DSQRT(TLP1/FPI)
      DO 2 M=1,L
      N2=LP1-M
      N1=N2+1
      N3=L+M
      FN2=DFLOAT(N2)
      DO 1 I=N1,N3
  1   FN2=FN2*DFLOAT(I)
  2   C(LM0+M)=DSQRT(TLP1/(FN2*TPI))
      RETURN
      END