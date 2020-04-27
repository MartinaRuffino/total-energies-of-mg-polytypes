      SUBROUTINE AMSTRP(NPHI,LPHI,NEL,EL,LMXA,POS,N0,NBAS,IPS,ALAT,
     .   CG,JCG,INDXCG,CY,NLA,NHS,B,BP)
C  STRUCTURE CONSTANTS AND ENERGY DERIVATIVES
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (O)
      DIMENSION NPHI(1),LPHI(N0,1),EL(1),POS(3,1),DR(3),LMXA(1),
     .  CG(1),JCG(1),INDXCG(1),CY(1),IPS(1),B(NLA,NHS),BP(NLA,NHS)
      CALL GETPR(IPR)
      CALL DPZERO(B,     NLA*NHS)
      CALL DPZERO(BP,    NLA*NHS)

C ------ LOOP OVER LMTO ENERGIES -----
      I1=1
      DO 10 IE=1,NEL
      E=EL(IE)
C ------ LOOP OVER CENTERS OF FUNCTIONS -----
      DO 10 IB=1,NBAS
      IS=IPS(IB)
      NLMB=(LPHI(IE,IS)+1)**2
      IF(NLMB.LE.0) GOTO 10
C --- LOOP OVER SITES WHERE EXPANDED -------
      J1=1
      DO 11 JB=1,NBAS
      JS=IPS(JB)
      NLMA=(LMXA(JS)+1)**2
      IF(JB.NE.IB) THEN
        DO 14 M=1,3
   14   DR(M)=ALAT*(POS(M,IB)-POS(M,JB))
        CALL MSTRXP(E,DR,B(J1,I1),BP(J1,I1),NLMA,NLMB,NLA,
     .      CG,INDXCG,JCG,CY)
        ENDIF
   11 J1=J1+NLMA

   10 I1=I1+NLMB

c      do 20 i=1,nhs
c   20 write(6,200) i,(b(j,i),j=1,nla)
c  200 format(i4,9f8.5/(4x,9f8.5))



      END