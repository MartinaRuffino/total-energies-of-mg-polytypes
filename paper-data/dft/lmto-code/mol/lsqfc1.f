C  PUTS CONSTRAINTS INTO A LEAST-SQUARE FIT (SEE BOOK I P65)
C
C   S         CHOLESKY FACTORIZED OVERLAP MATRIX FOR THE N BASIS FCTS
C   P(N,NC)   SCALPRODS OF BASISFCTS WITH THE NC CONSTRAINT FCTS
C   B(N,NB)   SCALPRODS OF THE NB FCTS TO BE FITTED WITH BASISFCTS,
C             FIT COEFFICIENTS ON OUTPUT
C   Q(NC,NB)  PRESCRIBED VALUES FOR SP OF FIT WITH CONSTRAINT FCTS
C   A         WORKSPACE OF DIMENSION (NC*(NC+1)/2
C   W         WORKSPACE OF DIMENSION N
C
C  LSQFC1(S,P,N,NC,A,W)        SETS UP P,A (NEEDED BY LSQFC2)
C  LSQFC2(B,Q,S,P,N,NB,NC,A)   DOES THE FIT FOR INPUT DATA B,Q
C ---------------------------------------------
      SUBROUTINE LSQFC1(S,P,N,NC,A,W)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION S(1),P(N,NC),W(N),A(1)
      IF(NC.LE.0) STOP' *** LSQFC1 CALLED WITH NC LE 0'
C OVERWRITE P WITH SINV*P AND MAKE A=PT*SINV*P
      IA=0
      DO 10 IC=1,NC
      DO 11 I=1,N
  11  W(I)=P(I,IC)
      CALL CHLR2S(S,P(1,IC),N,1)
      DO 10 JC=1,IC
      IA=IA+1
      SUM=0.D0
      DO 12 I=1,N
  12  SUM=SUM+P(I,JC)*W(I)
  10  A(IA)=SUM
C CHOLESKY FACTORIZATION OF A
      CALL CHLR2F(A,W,NC,NDEF)
      RETURN
      END
C ---------------------------------------
      SUBROUTINE LSQFC2(B,Q,S,P,N,NB,NC,A)
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION S(1),B(N,NB),P(N,NC),Q(NC,NB),A(1)
      IF(NC.LE.0) STOP '*** LSQFC2 CALLED WITH NC LE 0'
C OVERWRITE Q WITH Q-PT*SINV*B
      DO 15 IRHS=1,NB
      DO 15 I=1,N
      DO 15 IC=1,NC
  15  Q(IC,IRHS)=Q(IC,IRHS)-P(I,IC)*B(I,IRHS)
C MAKE LAMBDA, THEN C = SINV*B + SINV*P*LAMBDA
      CALL CHLR2S(S,B,N,NB)
      CALL CHLR2S(A,Q,NC,NB)
      DO 20 IRHS=1,NB
      DO 20 IC=1,NC
      DO 20 I=1,N
  20  B(I,IRHS)=B(I,IRHS)+P(I,IC)*Q(IC,IRHS)
      RETURN
      END