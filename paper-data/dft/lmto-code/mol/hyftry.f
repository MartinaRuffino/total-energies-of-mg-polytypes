      SUBROUTINE HYFTRY(R1,LP1,EP1,NX1,LX1,EX1,R2,LP2,EP2,NX2,LX2,EX2,
     .   S1,KP1,FP1,MX1,KX1,FX1,S2,KP2,FP2,MX2,KX2,FX2,
     .   JXI1,JXI2,LOK,RDIFF)
C  COMPARES PARAMETERS FOR WHICH THE TCF IS NEEDED WITH
C  SPECIFICATIONS OF A TABLE. RETURNS LOK=1 IF TABLE IS OK.
C  THEN RETURNS POINTERS FOR THE XI-FCTS AND RADII DIFFERENCE.
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (o)
      DIMENSION LX1(1),EX1(1),LX2(1),EX2(1),
     .   KX1(10),FX1(10),KX2(10),FX2(10),JXI1(10),JXI2(10)
      CALL GETPR(IPR)
      FUZZ=1.D-10
      LOK=0
      rdiff=0
C ------- COMPARE RADII AND PHI'S -------------
      IF(S1-R1.GT.FUZZ.OR.S2-R2.GT.FUZZ) RETURN
      IF(DABS(EP1-FP1).GT.FUZZ) RETURN
      IF(DABS(EP2-FP2).GT.FUZZ) RETURN
      IF(KP1.LT.LP1.OR.KP2.LT.LP2) RETURN
C ------- LOOK FOR TABLE-XI'S IN THOSE OF CALC -------
      DO 21 I=1,MX1
      J0=0
      DO 22 J=1,NX1
  22  IF(LX1(J).GE.KX1(I).AND.DABS(EX1(J)-FX1(I)).LT.FUZZ) J0=J
      IF(J0.EQ.0) RETURN
  21  JXI1(I)=J0
      DO 31 I=1,MX2
      J0=0
      DO 32 J=1,NX2
  32  IF(LX2(J).GE.KX2(I).AND.DABS(EX2(J)-FX2(I)).LT.FUZZ) J0=J
      IF(J0.EQ.0) RETURN
  31  JXI2(I)=J0
C ------- HERE IF THE TABLE CAN BE USED ------------
      RDIFF=DMAX1(R1-S1,R2-S2)
      LOK=1
      RETURN
      END