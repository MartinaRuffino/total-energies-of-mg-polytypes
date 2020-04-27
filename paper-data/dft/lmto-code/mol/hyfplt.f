      SUBROUTINE HYFPLT(LX1,EX1,NX1,LX2,EX2,NX2,DR,R1,R2,
     .    RSM1,RSM2,EP1,EP2,C1,C2,NF1,NF2,NLM1,NLM2)
C  ADD TOGETHER ROTATED TC-FIT AT POINTS, COMPARE TO PRODUCT
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (o)
      PARAMETER(  NX0=25, NY0=25, NPMX=3)
      DIMENSION RMAX(2),BOT(10),ILM0(NPMX),JLM0(NPMX),DR(3),
     .   TAU(3,2),TOP(10),XV(3),YV(3),ZV(3),CPAR(3,2)
      DIMENSION EX1(1),EX2(1),EP1(1),EP2(1),CY(100),PH1(100),PH2(100),
     .   P(3),P1(3),P2(3),LX1(1),LX2(1),
     .   FAC1(0:10),FAC2(0:10),FAC(NPMX),XI1(500),XI2(500),
     .   C1(NF1,NLM1,NLM2),C2(NF2,NLM2,NLM1)
      REAL PROD(NX0,NY0,NPMX),FIT(NX0,NY0,NPMX),XZSCM,YSZCM
      CALL GETPR(IPR)
      CALL SYLMNC(CY,8)
      ASM1=1.D0/RSM1
      ASM2=1.D0/RSM2
      LP1=LL(NLM1)
      LP2=LL(NLM2)
      DO 1 M=1,3
      TAU(M,1)=-DR(M)*0.5D0
  1   TAU(M,2)= DR(M)*0.5D0
      RMAX(1)=R1
      RMAX(2)=R2
      WRITE(6,338) TAU
  338 FORMAT(' TAU1=',3F9.4,'   TAU2=',3F9.4)
C ------- GET MAX VALUES ON SPHERE TO NORMALIZE PHI1,PHI2 -----
      CALL HRMSAV(R1,LP1,EP1,FAC1)
      CALL HRMSAV(R2,LP2,EP2,FAC2)
      IF(IPR.GE.30) WRITE(6,330) 1,(FAC1(L),L=0,LP1)
      IF(IPR.GE.30) WRITE(6,330) 2,(FAC2(L),L=0,LP2)
  330 FORMAT(' RMS AVG PHI AT R',I1,': ',6F9.4)
C ------- DEFINE PLANE ---------------
      WRITE(6,*) 'ENTER XV,YV,X1,X2,Y1,Y2,Z0'
      READ (5,*)  XV,YV,X1,X2,Y1,Y2,Z0
      XNRM=DSQRT(SP3(XV,XV))
      DO 35 M=1,3
  35  XV(M)=XV(M)/XNRM
      XX=SP3(XV,YV)
      DO 36 M=1,3
  36  YV(M)=YV(M)-XX*XV(M)
      YNRM=DSQRT(SP3(YV,YV))
      DO 37 M=1,3
  37  YV(M)=YV(M)/YNRM
      CALL CROSS(XV,YV,ZV)
      WRITE(6,882) XV,YV,ZV
  882 FORMAT(' XV=',3F7.3,'   YV=',3F7.3,'   ZV=',3F7.3)
      DO 38 IB=1,2
      CPAR(1,IB)=SP3(XV,TAU(1,IB))
      CPAR(2,IB)=SP3(YV,TAU(1,IB))
      ZDIST=SP3(ZV,TAU(1,IB))-Z0
      CPAR(3,IB)=.01
      XXX=RMAX(IB)**2-ZDIST**2
      IF(XXX.GT.0.D0) CPAR(3,IB)=DSQRT(XXX)
  38  WRITE(6,997) IB,CPAR(1,IB),CPAR(2,IB),CPAR(3,IB)
  997 FORMAT(' IB=',I5,'   CPAR=',3F9.5)
C ------- SELECT PICTURES ----------
      WRITE(6,*) 'ENTER NX,NY,NPICT, PAIRS ILM,JLM'
      READ (5,*)  nx,ny,NPICT,(ILM0(II),JLM0(II),II=1,NPICT)
      DO 47 II=1,NPICT
      L1=LL(ILM0(II))
      L2=LL(JLM0(II))
      FAC(II)=1.D0/(FAC1(L1)*FAC2(L2))
      BOT(II)=1.D10
  47  TOP(II)=-1.D10
C ------- LOOP OVER POINTS ----------
      DO 50 IX=1,NX
      DO 50 IY=1,NY
      xxx=x1
      yyy=y1
      if (nx > 1) XXX=X1+(IX-1)*(X2-X1)/(NX-1.D0)
      if (ny > 1) YYY=Y1+(IY-1)*(Y2-Y1)/(NY-1.D0)
      DO 55 M=1,3
      P(M)=XXX*XV(M)+YYY*YV(M)+Z0*ZV(M)
      P1(M)=P(M)-TAU(M,1)
  55  P2(M)=P(M)-TAU(M,2)
C ------- ADD TOGETHER FIT -----------
      CALL SOLHSM(P1,EP1(1),ASM1,LP1,PH1,CY)
      CALL SOLHSM(P2,EP2(1),ASM2,LP2,PH2,CY)
      CALL HYFXIV(LX1,EX1,NX1,P1,ASM1,CY,XI1)
      CALL HYFXIV(LX2,EX2,NX2,P2,ASM2,CY,XI2)
      DO 15 II=1,NPICT
        ILM1=ILM0(II)
        ILM2=JLM0(II)
        SUM=0.D0
        DO 11 J=1,NF1
  11    SUM=SUM+C1(J,ILM1,ILM2)*XI1(J)
        DO 12 J=1,NF2
  12    SUM=SUM+C2(J,ILM2,ILM1)*XI2(J)

C        print 345, p, ph1(ilm1)*ph2(ilm2)*fac(ii), sum*fac(ii)
C  345 format(3f8.4,3f11.5)

        FIT(IX,IY,II)=SUM*FAC(II)
        PRD=PH1(ILM1)*PH2(ILM2)*FAC(II)
        TOP(II)=DMAX1(TOP(II),PRD)
        BOT(II)=DMIN1(BOT(II),PRD)
  15    PROD(IX,IY,II)=PRD
  50  CONTINUE
      WRITE(6,561) 'BOT',(BOT(II),II=1,NPICT)
      WRITE(6,561) 'TOP',(TOP(II),II=1,NPICT)
  561 FORMAT(1X,A3,5F12.6)
C ------- DO THE PLOT --------------------------
C*    CALL BEGPLT(25.,18.,'H')
      YSZCM=4.8
      XSZCM=YSZCM*(X2-X1)/(Y2-Y1)
      IF(XSZCM.GT.24.) THEN
        YSZCM=YSZCM*24./XSZCM
        XSZCM=24.
        ENDIF
      WRITE(6,772) XSZCM,YSZCM
  772 FORMAT(' XSZCM=',F12.6,'   YSZCM=',F12.6)
C*    CALL POSICM(0.1,XSZCM,17.9-YSZCM,17.9)
C|    CALL LSIZE(.6,.6)
      DO 80 II=1,NPICT
      H1=BOT(II)
      H2=TOP(II)
C|    DH=(H2-H1)/9.99999999D0
      DH=(H2-H1)/19.99
        H1=-3.1D0
        H2=+3.1D0
        DH=0.199999
        IF(DMAX1(TOP(II),-BOT(II)).LE.1.D0) THEN
          H1=-1.05D0
          H2= 1.05D0
          DH=0.099999D0
          ENDIF
        IF(DMAX1(TOP(II),-BOT(II)).LE.0.2D0) THEN
          H1=-0.22D0
          H2= 0.22D0
          DH=0.03999999D0
          ENDIF
C|    WRITE(6,740) BOT(II),TOP(II)
C|740 FORMAT(' BOT=',F12.6,'   TOP=',F12.6,'    ENTER H1,H2,DH')
C|    READ (5,*)  H1,H2,DH
C*    IF(II.GT.1) CALL PSHIFT(0.0,-1.1)
C*    CALL AMVE(0.,-.08)
C*    CALL LABEL(' ILM=<I2>  JLM=<I2>$',ILM0(II),JLM0(II))
C*    CALL LABEL('  H1=<F7.4>  H2=<F7.4>  DH=<F7.4>$',SNGL(H1),
C*   .   SNGL(H2),SNGL(DH))
      print *, 'h=',h1,h2,dh
      open(81,file='prod')
      open(82,file='fit')
      write (81,*) nx,ny
      write (82,*) nx,ny
      do  77 i=1,nx
        write(81,810) (prod(i,j,ii), j=1,ny)
        write(82,810) (fit(i,j,ii), j=1,ny)
  810   format(5f15.7)
   77 continue
C      CALL PL2CP(PROD(1,1,II),NX,NY,X1,X2,Y1,Y2,0D0,0D0,H1,H2,DH,CPAR)
C      CALL PL2CP(FIT(1,1,II),NX,NY,X1,X2,Y1,Y2,1.1D0,0D0,H1,H2,DH,CPAR)
  80  CONTINUE
C*    CALL MVECM(0.,17.5)
C|    CALL LABEL('  NEXI=<I2>  EXI=$',NEXI)
C|    DO 46 IE=1,NEXI
C|46  CALL LABEL('<F7.2>  $',EXI(IE))
C|    CALL LABEL('   D=<F6.3><F6.3><F6.3>$',DR(1),DR(2),DR(3))
C|    CALL LABEL('   Z0=<F6.3>$',Z0)
C*    CALL ENDPLT
      STOP
      END
      DOUBLEPRECISION FUNCTION SP3(A,B)
      REAL*8 A(3),B(3),ADOTB
      SP3=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
      SUBROUTINE SOLHSM(R,E,A,LMAX,HL,CY)
C  REAL SMOOTHED SOLID HANKEL FUNCTIONS, NEGATIVE ENERGIES
      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
      DIMENSION R(3),HL(1),BL(1),CY(1),PSI(0:20)
      CALL SYLM(R,HL,LMAX,R2)
      R1=DSQRT(R(1)**2+R(2)**2+R(3)**2)

      if (e < -1d-6) call hansmr(r1,e,a,psi,lmax)
      if (dabs(e) < 1d-6) call ropgau(1/a,lmax,1,r1,psi,000)
      if (e > 1d-6) call ropgau(e,lmax,1,r1,psi,000)
      ILM=0
      DO 10 L=0,LMAX
      NM=2*L+1
      DO 10 M=1,NM
      ILM=ILM+1
  10  HL(ILM)=PSI(L)*CY(ILM)*HL(ILM)
      RETURN
      END
