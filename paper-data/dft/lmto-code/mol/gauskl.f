      SUBROUTINE GAUSKL(R,E,A,G,LMAX,KMAX,N0)
C  MAKE RADIAL PARTS OF G_KL, DIVIDED BY R**L.
C  G_KL IS:   (LAPLACE-OP)**K YL(-GRAD) G0
      PARAMETER( LX=16 )
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION G(0:N0,0:1)
      PI=4D0*DATAN(1.D0)
      TA2=2.D0*A*A
      IF(KMAX.LT.0.OR.LMAX.LT.0) RETURN
      G0=DEXP(-A*A*R*R) * (A*A/PI)**1.5D0 * DEXP(E/(4D0*A*A))
C ---------- DO EXPLICITLY FOR K=0,1 ------
      DO 6 L=0,LMAX
  6   G(0,L)=TA2**L * G0
      IF(KMAX.EQ.0) RETURN
      DO 7 L=0,LMAX
  7   G(1,L)=TA2**(L+1)*(2D0*A*A*R*R-3-2*L) * G0
C ---------- RECURSION FOR HIGHER K ----------
      DO 8 K=2,KMAX
      DO 8 L=0,LMAX
      X=2*(K-1)*(2*K+2*L-1)
      Y=4*K+2*L-1
  8   G(K,L)=TA2*(TA2*R*R*G(K-1,L)-Y*G(K-1,L)-X*TA2*G(K-2,L))
      RETURN
      END
