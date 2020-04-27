      SUBROUTINE AUGSKJ(E1,E2,R,LMAX,SAB,SKK,SKJ,SJK,SJJ)
C  INTEGRALS TO AUGMENT FUNCTIONS OF ENERGY E1,E2
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(O)
      DIMENSION SAB(4,1),SKK(1),SKJ(1),SJK(1),SJJ(1),
     .   AK1(10),AK2(10),DK1(10),DK2(10),
     .   AJ1(10),AJ2(10),DJ1(10),DJ2(10)

      CALL RADKJ(E1,R,LMAX,AK1,AJ1,DK1,DJ1,0)
      CALL RADKJ(E2,R,LMAX,AK2,AJ2,DK2,DJ2,0)
      DO 10 I=1,LMAX+1
      S00=SAB(1,I)
      S01=SAB(2,I)
      S10=SAB(3,I)
      S11=SAB(4,I)
      SKK(I)= AK1(I)*S00*AK2(I) + AK1(I)*S01*DK2(I)
     .      + DK1(I)*S10*AK2(I) + DK1(I)*S11*DK2(I)

      SKJ(I)= AK1(I)*S00*AJ2(I) + AK1(I)*S01*DJ2(I)
     .      + DK1(I)*S10*AJ2(I) + DK1(I)*S11*DJ2(I)

      SJK(I)= AJ1(I)*S00*AK2(I) + AJ1(I)*S01*DK2(I)
     .      + DJ1(I)*S10*AK2(I) + DJ1(I)*S11*DK2(I)

      SJJ(I)= AJ1(I)*S00*AJ2(I) + AJ1(I)*S01*DJ2(I)
     .      + DJ1(I)*S10*AJ2(I) + DJ1(I)*S11*DJ2(I)

c      WRITE(6,300) I-1,SKK(I),SKJ(I),SJK(I),SJJ(I)
c  300 FORMAT(I4,'  KK,KJ,JK,JJ',4F12.5)

   10 CONTINUE

      END
