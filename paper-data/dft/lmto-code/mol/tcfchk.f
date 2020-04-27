      SUBROUTINE TCFCHK(MMAX,LX1,NX1,LX2,NX2,XP,ZP,WP,NP,
     .  nx,lx,XIE,NLM,lmax,PHI1,PHI2,ri1,ri2,w1,w2,
     .  R1,R2,ZC1,ZC2,LP1,EP1,LP2,EP2,LSYM,B,NDIM,NRHS,
     .  NDIMX,NRHSX,ERR1,ERR2,prd,BOT,TOP,ERRS)
C- Generate rms and max errors of the TCF.
C  Each wave function is normalized to have an RMS average=1 at its RMT.
      PARAMETER( N0=10, NXMX=4 )
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (o)
      DIMENSION NDIM(0:20),XP(1),ZP(1),WP(1),IRHS(0:20),
     .   XI1(0:N0,0:N0,NXMX),XIM(200,0:20),
     .   PH1(0:N0,0:N0),PH2(0:N0,0:N0),TOP(NRHSX,0:MMAX),CX(200),
     .   EE1(1),EE2(1),B(NDIMX,NRHSX,0:MMAX),BOT(NRHSX,0:MMAX),
     .   ERR1(NRHSX,0:MMAX),ERR2(NRHSX,0:MMAX),NRHS(0:20),
     .   FAC1(0:N0),FAC2(0:N0),ERRS(0:LP1,0:LP2),prd(NRHSX,0:MMAX),
     .  lx(1),xie(np,nlm,1),phi1(np,1),LX1(6),phi2(np,1),LX2(6)

      CHARACTER*1 CL(0:8)
      DATA CL /'S','P','D','F','G','5','6','7','8'/

      CALL GETPR(IPR)
      IF(IPR.GE.20) WRITE(6,345)
  345 format(/' --------- tcfchk: -----------------')
      ERRALF=1.5D0
      EE1(1)=EP1
      EE2(1)=EP2
      CALL SXLMNC(CX,10)
      IF(MAX0(NX1,NX2).GT.NXMX) CALL RX('TCFINT: INCREASE LMXX''')
      ERR0=0.02D0
      DO 1 M=0,MMAX
      DO 1 I=1,NRHS(M)
      TOP(I,M)=-1000.D0
      BOT(I,M)=1000.D0
      prd(I,M) =0.d0
      ERR1(I,M)=0.D0
  1   ERR2(I,M)=0.D0
      ERRVOL=0.D0
C ------ GET AVG VALUES ON SPHERE TO SCALE ERRORS -------
      CALL HRMSAV(R1,LP1,EE1,FAC1)
      CALL HRMSAV(R2,LP2,EE2,FAC2)
      IF(IPR.GE.40) WRITE(6,330) 1,(FAC1(L),L=0,LP1)
      IF(IPR.GE.40) WRITE(6,330) 2,(FAC2(L),L=0,LP2)
  330 FORMAT(' RMS AVG OF PHI AT R',I1,': ',6F9.4)
C ------ START LOOP OVER BISINT POINTS ------------------
      DO 10 IP=1,NP
      X=XP(IP)
      Z=ZP(IP)
      D1=DSQRT(X*X+(Z-ZC1)**2)
      D2=DSQRT(X*X+(Z-ZC2)**2)

C ------ FOR EACH M, COLLECT FCTS IN THAT M-SPACE ------
      DO 30 M=0,MMAX
      lm = ((2*lmax-m+1)*m)/2 + 1
      II=0
      DO 31 IE=1,NX
      DO 31 L=M,LX(IE)
      II=II+1
  31  XIM(II,M)=xie(ip,lm+l,ie)
  30  IRHS(M)=0
C ------ ADD TO ERR INTEGRALS ------------------
      DO 20 L1=0,LP1
      DO 20 M1=0,L1
      LTOP=LP2
      IF(LSYM.EQ.1) LTOP=L1
      DO 21 L2=0,LTOP
      MTOP=L2
      IF(LSYM.EQ.1.AND.L2.EQ.L1) MTOP=M1
      DO 21 M2=0,MTOP
      XXX=phi1(ip,(l1*(l1+1))/2+m1+1)*phi2(ip,(l2*(l2+1))/2+m2+1)
      MP=M1+M2
      MM=IABS(M1-M2)
      IF(MP.LE.MMAX) THEN
        IRHS(MP)=IRHS(MP)+1
        JRHS=IRHS(MP)
        FIT=0.D0
        DO 25 I=1,NDIM(MP)
  25    FIT=FIT+XIM(I,MP)*B(I,JRHS,MP)
        ERR1(JRHS,MP)=DMAX1(ERR1(JRHS,MP),DABS(XXX-FIT))
        ERR2(JRHS,MP)=ERR2(JRHS,MP)+WP(IP)*(XXX-FIT)**2
        prd(jrhs,mp) =prd(jrhs,mp) +WP(IP)*XXX**2
        BOT(JRHS,MP)=DMIN1(BOT(JRHS,MP),XXX)
        TOP(JRHS,MP)=DMAX1(TOP(JRHS,MP),XXX)
        ENDIF
      IF(MM.LE.MMAX.AND.MM.NE.MP) THEN
        IRHS(MM)=IRHS(MM)+1
        JRHS=IRHS(MM)
        FIT=0.D0
        DO 26 I=1,NDIM(MM)
  26    FIT=FIT+XIM(I,MM)*B(I,JRHS,MM)
        ERR1(JRHS,MM)=DMAX1(ERR1(JRHS,MM),DABS(XXX-FIT))
        ERR2(JRHS,MM)=ERR2(JRHS,MM)+WP(IP)*(XXX-FIT)**2
        prd(jrhs,mm) =prd(jrhs,mm) +WP(IP)*XXX**2
        BOT(JRHS,MM)=DMIN1(BOT(JRHS,MM),XXX)
        TOP(JRHS,MM)=DMAX1(TOP(JRHS,MM),XXX)
        ENDIF
  21  CONTINUE
  20  CONTINUE

  10  CONTINUE
C ------ COLLECT ERRORS BY L1,L2 AND OUTPUT BIG TABLE -----
      IF(IPR.GE.40) WRITE(6,331) ERRVOL,ERRALF
  331 FORMAT(' RMS-ERR:  VOL=',F8.2,'   FROM RADII *',F6.2)
      YY=100D0
      IF(IPR.GE.40) WRITE(6,992)
  992   FORMAT(/'  CASE   L1 M1   L2 M2',
     .    '    MAX ABS VAL     M    RMS-ERR  MAX-ERR (%)')
      DO 45 L1=0,LP1
      DO 45 L2=0,LP2
      XI1(L1,L2,1)=0.D0
      PH1(L1,L2)=0.D0
  45  PH2(L1,L2)=0.D0
      DO 29 M=0,MMAX
  29  IRHS(M)=0
      DO 40 L1=0,LP1
      DO 40 M1=0,L1
      LTOP=LP2
      IF(LSYM.EQ.1) LTOP=L1
      DO 41 L2=0,LTOP
      MTOP=L2
      IF(LSYM.EQ.1.AND.L2.EQ.L1) MTOP=M1
      DO 41 M2=0,MTOP
      MP=M1+M2
      MM=IABS(M1-M2)
      FAC=1.D0/(FAC1(L1)*FAC2(L2))
C ---------- PART FOR M3=M1+M2 --------------------
      IF(MP.LE.MMAX) THEN
        IRHS(MP)=IRHS(MP)+1
        J=IRHS(MP)
        ERR1(J,MP)=ERR1(J,MP)*fac
        ERR2(J,MP)=DSQRT(ERR2(J,MP)/prd(j,mp))
        PMAX=DMAX1(TOP(J,MP),-BOT(J,MP))*FAC
        PH1(L1,L2)=DMAX1(PH1(L1,L2),ERR1(J,MP))
        PH2(L1,L2)=DMAX1(PH2(L1,L2),ERR2(J,MP))
        XI1(L1,L2,1)=DMAX1(XI1(L1,L2,1),PMAX)
        ERR=ERR2(J,MP)
        IF(IPR.GE.50.OR.IPR.GE.40.AND.ERR.GT.ERR0)
     .    WRITE(6,991) CL(L1),CL(L2),L1,M1,L2,M2,PMAX,
     .                  MP,ERR2(J,MP)*YY,ERR1(J,MP)*YY
  991   FORMAT(1X,A1,' * ',A1,2X,2I3,2X,2I3,F13.4,I8,1X,2F9.2)
        ENDIF
C ---------- PART FOR M3=IABS(M1-M2) --------------
      IF(MM.LE.MMAX.AND.MM.NE.MP) THEN
        IRHS(MM)=IRHS(MM)+1
        J=IRHS(MM)
        ERR1(J,MM)=ERR1(J,MM)*FAC
        ERR2(J,MM)=DSQRT(ERR2(J,MM)/prd(j,mm))
        PH1(L1,L2)=DMAX1(PH1(L1,L2),ERR1(J,MM))
        PH2(L1,L2)=DMAX1(PH2(L1,L2),ERR2(J,MM))
        IF(IPR.GE.50.OR.IPR.GE.40.AND.ERR.GT.ERR0) WRITE(6,993)
     .     MM,ERR2(J,MM)*YY,ERR1(J,MM)*YY
  993   FORMAT(37X,I6,1X,2F9.2)
        ENDIF
  41  CONTINUE
  40  CONTINUE
C ------ PUT INTO ERRS, OUTPUT SMALL TABLE ----------
      IF(IPR.GE.20) WRITE(6,671) ZC2-ZC1
CL      WRITE(71,671) ZC2-ZC1
      DO 65 L1=0,LP1
      LTOP=LP2
      IF(LSYM.EQ.1) LTOP=L1
      DO 65 L2=0,LTOP
      PX=XI1(L1,L2,1)
      IF(IPR.GE.20)
     . WRITE(6,670) CL(L1),CL(L2),PX,YY*PH2(L1,L2),YY*PH1(L1,L2)
CL      WRITE(71,670) CL(L1),CL(L2),PX,YY*PH2(L1,L2),YY*PH1(L1,L2)
      ERRS(L1,L2)=PH2(L1,L2)*100D0
      IF(LSYM.EQ.1) ERRS(L2,L1)=ERRS(L1,L2)
  65  CONTINUE
  670 FORMAT(3X,A1,' * ',A1,F13.4,1X,2F11.2)
  671 FORMAT(/' TCFCHK:   D=',F10.5/'    CASE',
     .   '    MAX ABS VAL    RMS-ERR   MAX-ERR (%)')

      END
