      SUBROUTINE SVSVEC(E,LMXA,RMAX,IPCLAS,NBAS,GH,GJ,NLA)
C  MAKES VECTORS OF HANKEL AND BESSEL VALUES AND SLOPES INTO GH,GJ
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION LMXA(1),IPCLAS(1),RMAX(1),PHI(0:10),PSI(0:10),
     .   GH(NLA,2),GJ(NLA,2)
C ----------- START LOOP OVER ATOMS -----------
      ILA=0
      DO 10 IB=1,NBAS
      IC=IPCLAS(IB)
      LMAX=LMXA(IC)
      IF(LMAX.GT.9) CALL RX('SVSVEC:  INCREASE DIMENSIONS''')
      R=RMAX(IC)
      ER2=E*R*R
      CALL BESSL(ER2,LMAX+1,PHI,PSI)
      RL=1.D0/R
      DO 11 L=0,LMAX
      RL=RL*R
      VJ=RL*PHI(L)
      SJ=RL*(L*PHI(L)-ER2*PHI(L+1))/R
      VH=PSI(L)/(RL*R)
      SH=(L*PSI(L)-PSI(L+1))/(RL*R*R)
      DO 11 M=1,2*L+1
      ILA=ILA+1
      GH(ILA,1)=VH
      GH(ILA,2)=SH
      GJ(ILA,1)=VJ
  11  GJ(ILA,2)=SJ
  10  CONTINUE

      RETURN
      END