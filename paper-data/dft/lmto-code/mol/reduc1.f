      SUBROUTINE REDUC1(NAB,N,A,B,DL)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION A(NAB,1),B(NAB,1),DL(1)
      DO 1 I=1,N
      IM1=I-1
      DO 1 J=1,N
      X=B(I,J)
      IF(IM1.EQ.0) GOTO 20
      DO 2 IK=1,IM1
      K=I-IK
  2   X=X-B(I,K)*B(J,K)
  20  IF(I.NE.J) GOTO 3
      IF(X.LE.0.D0) GOTO 100
      Y=DSQRT(X)
      DL(I)=Y
      GOTO 1
  3   B(J,I)=X/Y
  1   CONTINUE
      DO 4 I=1,N
      IM1=I-1
      Y=DL(I)
      DO 4 J=1,N
      X=A(I,J)
      IF(IM1.EQ.0) GOTO 4
      DO 5 IK=1,IM1
      K=I-IK
  5   X=X-B(I,K)*A(J,K)
  4   A(J,I)=X/Y
      DO 6 J=1,N
      JM1=J-1
      DO 6 I=J,N
      IM1=I-1
      X=A(I,J)
      IF(IM1.LT.J) GOTO 9
      DO 7 IK=J,IM1
      K=IM1+J-IK
  7   X=X-A(K,J)*B(I,K)
  9   IF(JM1.EQ.0) GOTO 6
      DO 8 IK=1,JM1
      K=J-IK
  8   X=X-A(J,K)*B(I,K)
  6   A(I,J)=X/DL(I)
C ..... ADDED: FILL UPPER TRIANGULAR ALSO
C|    DO 25 I=1,N
C|    DO 25 J=1,I
C|25  A(J,I)=A(I,J)
      RETURN
  100 WRITE(6,101)
  101 FORMAT(' *** REDUC1: OVERLAP MATRIX NOT POSITIVE DEFINITE')
      STOP
      END
C --------------------------------
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      DO 100 I = 1, N
      DO 100 J = 1, I
  100 Z(I,J) = A(I,J)
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.D0
         SCALE = 0.D0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(Z(I,K))
         IF (SCALE .NE. 0.D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
         F = Z(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.D0
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.D0
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
         DO 280 K = 1, L
  280    Z(I,K) = SCALE * Z(I,K)
  290    D(I) = H
  300 CONTINUE
  320 D(1) = 0.D0
      E(1) = 0.D0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.D0) GO TO 380
         DO 360 J = 1, L
            G = 0.D0
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
  380    D(I) = Z(I,I)
         Z(I,I) = 1.D0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.D0
            Z(J,I) = 0.D0
  400    CONTINUE
  500 CONTINUE
      RETURN
      END
C --------------------------------
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      REAL*8 D(N),E(N),Z(NM,N),MACHEP
      DATA MACHEP/3.D-14/
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      E(N) = 0.D0
      DO 240 L = 1, N
         J = 0
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (DABS(E(M)) .LE. MACHEP * (DABS(D(M)) + DABS(D(M+1))))
     .         GO TO 120
  110    CONTINUE
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = DSQRT(G*G+1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF ( DABS(F) .LT. DABS(G)) GO TO 150
            C = G / F
            R =  DSQRT(C*C+1.0D0)
            E(I+1) = F * R
            S = 1.0D0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = DSQRT(S*S+1.0D0)
            E(I+1) = G * R
            C = 1.0D0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
  200    CONTINUE
         D(L) = D(L) - P
         E(L) = G
C        E(M) = DBLE(0.000)
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END
C -------------------------
      SUBROUTINE REBAKA(NDIM,N,M1,M2,B,DL,Z)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION B(NDIM,1),DL(1),Z(NDIM,1)
      DO 1 J=M1,M2
      DO 1 I=1,N
      IN=N+1-I
      X=Z(IN,J)
      IN1=IN+1
      IF(IN.EQ.N) GOTO 1
      DO 2 K=IN1,N
  2   X=X-B(K,IN)*Z(K,J)
  1   Z(IN,J)=X/DL(IN)
      RETURN
      END
