      SUBROUTINE HYFSTF(ADEC,A1,D0,DIST,WIST,F,NDIST,NALF)
C  SET FUNCTIONS FOR THE NDIST DISTANCES
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (o)
      DIMENSION DIST(NDIST),F(NDIST,NALF),WIST(NDIST)
      CALL GETPR(IPR)
      IF(IPR.GE.20) WRITE(6,200) NALF,ADEC,A1,D0,NDIST
  200 FORMAT(/' HYFSTF:  NALF=',I2,'   ADEC=',F7.3,'   A1=',F7.3,
     .   '   D0=',F7.3,'   NDIST=',I3)
      DO 20 IDIST=1,NDIST
      X=IDIST/(NDIST+0.D0)
  20  DIST(IDIST)=D0-DLOG(X)/ADEC
C -------- MAKE F AND WIST -------------
      IF(IPR.GE.60) WRITE(6,331)
      DO 10 IDIST=1,NDIST
      D=DIST(IDIST)
      X=DEXP(A1*(D0-D))
      WIST(IDIST)=1.D0
      XM=1.D0
      DO 11 IALF=1,NALF
      XM=XM*X
  11  F(IDIST,IALF)=XM
      IF(IPR.GE.60) WRITE(6,330) IDIST,D,X,(F(IDIST,IA),IA=1,NALF)
  330 FORMAT(I5,2F11.5,4F12.6:/(27X,4F12.6))
  331 FORMAT(' IDIST    DIST',8X,'X',8X,'F(1..NALF)')
  10  CONTINUE
      RETURN
      END