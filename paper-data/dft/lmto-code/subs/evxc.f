C#ifdefC VOSKO
C      SUBROUTINE EVXC(RHO,RHOSP,EC,UC1)
CC- Makes XC energy density and potential
CC  ---------------------------------------------------
CCi Inputs:
CCi   RHO,RHOSP
CCi     RHO   is the total electron density
CCi     RHOSP is the electron density for one spin
CCo Outputs:
CCo   EC, UC1 : energy density and potential
CCo    Vosko-Ceperley-Alder
CC  ---------------------------------------------------
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/.0621814d0,-.10498d0,3.72744d0,
C     .  12.9352d0,6.1519908d0,1.2117833d0,1.1435257d0,-.031167608d0/
C      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/.0310907d0,-.32500d0,7.06042d0,
C     .  18.0578d0,4.7309269d0,2.9847935d0,2.7100059d0,-.1446006d0/
C
C      ec = 0d0
C      uc1 = 0d0
C      th = 1d0/3d0
C      th4 = 4d0/3d0
C      if( rho < 1d-18 ) return
C      rs = ( rho * 4.1887902d0 )**(-1d0/3d0)
C      s = 2d0 * rhosp / rho - 1d0
C      x = dsqrt(rs)
C      xpx = x*x+bp*x + cp
C      xfx = x*x+bf*x + cf
C      s4 = s**4 - 1d0
C      FS = ((1d0+S)**th4 + (1d0-S)**th4 - 2d0) / (2d0**th4 - 2d0)
C      BETA = 1d0/(1.74208d0 + 3.182d0*X+.09873d0*X*X+.18268d0*X**3)
C      DFS = th4*((1d0+S)**th-(1d0-S)**th)/(2d0**th4-2d0)
C      DBETA = -(.27402d0*X+.09873d0+1.591d0/X)*BETA**2
C      ATNP = dATAN(QP/(2d0*X+BP))
C      ATNF = dATAN(QF/(2d0*X+BF))
C      ECP = AP*(dLOG(X*X/XPX)+CP1*ATNP
C     . -CP3*(dLOG((X-XP0)**2/XPX)+CP2*ATNP))
C      ECF = AF*(dLOG(X*X/XFX)+CF1*ATNF
C     . -CF3*(dLOG((X-XF0)**2/XFX)+CF2*ATNF))
C      EC = ECP+FS*(ECF-ECP)*(1d0+S4*BETA)
C      TP1 = (X*X+BP*X)/XPX
C      TF1 = (X*X+BF*X)/XFX
C      UCP = ECP-AP/3d0*(1d0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))
C      UCF = ECF-AF/3d0*(1d0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))
C      UC0 = UCP+(UCF-UCP)*FS
C      UC10 = UC0-(ECF-ECP)*(S-1d0)*DFS
C      DUC = (UCF-UCP)*BETA*S4*FS
C     . +(ECF-ECP)*(-RS/3d0)*DBETA*S4*FS
C      DUC1 = DUC-(ECF-ECP)*BETA*(S-1)*(4d0*S**3*FS+S4*DFS)
C      UC1 = UC10+DUC1
Cc     WRITE(*,101)RS,S,EC,UC1
C  101 FORMAT(F6.1,F5.1,7F9.2)
C      UC1 = UC1  -1.2217736d0  * (1d0 + s)**th/ RS
C      EC  = EC -0.9163306d0/RS -0.2381735d0/RS * FS
Cc     WRITE(*,101)RS,S,EC,UC1
C      RETURN
C      END
C#else
      SUBROUTINE EVXC(RHO,RHOSP,EXC,VXC)
C- Makes XC energy density and potential
C  ---------------------------------------------------
Ci Inputs:
Ci   RHO,RHOSP
Ci     RHO   is the total electron density
Ci     RHOSP is the electron density for one spin
Co Outputs:
Co   EXC,VXC: energy density and potential (v = d (rho eps) /rho)
Cr Remarks
Cr   Von Barth, Hedin should be used with gradient corrections
C  ---------------------------------------------------
C ----------------------------------------------------------------------
      IMPLICIT NONE
C ... Passed parameters
      DOUBLE PRECISION RHO,RHOSP,EXC,VXC
C ... Local parameters
      DOUBLE PRECISION AF,AP,CF,CMUXP,CP,EFMEP,EPSCF,EPSCP,EPSX0,EPSX1,
     .  F,FOFXI,FOTHRD,FPI,FPIRO,GAMMA,MUCF,MUCP,NUCE,ONTHRD,RS,TAUCE,X,
     .  XI
      DATA ONTHRD,FOTHRD,FPI/.333333333D0,1.333333333D0,12.566370614D0/
      DATA GAMMA,CMUXP/5.1297628D0,1.2217741D0/
      DATA EPSX0,EPSX1/.9163305866D0,.23817361D0/
C ------ V.BARTH, HEDIN PARAMETRIZATION ------------
C      DATA CP,CF/-.0504D0,-.0254D0/,AP,AF/30D0,75d0/
C ------- TAKEN FROM ASW ---------------
      DATA CP,CF/-.045D0,-.0225D0/,AP,AF/21D0,52.9167D0/
C --------------------------------------
      F(X) = (1D0 + X*X*X)*DLOG(1D0 + 1D0/X) + .5D0*X - X*X - ONTHRD
      FOFXI(X) = (X**FOTHRD +(1-X)**FOTHRD - .793700526D0)/.206299474D0
      FPIRO = FPI*RHO
c      FPIROS = FPI*RHOSP
      VXC = 0D0
      EXC = 0D0
      IF (FPIRO .LE. 1D-18) RETURN
      RS = (3D0/FPIRO)**ONTHRD
      XI = RHOSP/RHO
      X = RS/AP
      MUCP = CP*DLOG(1D0 + 1D0/X)
      EFMEP = -CP*F(X)
      EPSCP = -EFMEP
      X = RS/AF
      MUCF = CF*DLOG(1D0 + 1D0/X)
      EPSCF = CF*F(X)
      EFMEP = EFMEP + EPSCF
      NUCE = GAMMA*EFMEP
      TAUCE = MUCF - MUCP - FOTHRD*EFMEP
      VXC = (-CMUXP/RS + NUCE)*((2D0*XI)**ONTHRD) - NUCE + MUCP
     .   +TAUCE*FOFXI(XI)
      EXC = -EPSX0/RS + EPSCP + FOFXI(XI)*(EPSCF - EPSCP - EPSX1/RS)
      END
C#endif
C#ifdefC PERD91
C      subroutine nlocxc(rho1,rho2,drho,dabgrh,ddrho,drho1,drho2,
C     .  dabgrh1,dabgrh2,ddrho1,ddrho2,nsp,vxc,exc)
CC- Perdew-Wang gradient correction, version CGA91
CC---------------------------------------------------------------
CCi Inputs:
CCi   rho1,rho2,drho,dabgrh,ddrho,drho1,drho2,dabgrh1,dabgrh2,
CCi   ddrho1,ddrho2,nsp
CCi   charge, gradient and laplacian for both spins.
CCo Outputs:
CCo   vxc, exc
CCr Remarks:
CCr   exc is the energy density times the charge.  Perdew's remarks:
CC   Generates the generalized gradient approximation for
CC   the exchange-correlation energy and potential (j.p. perdew
CC   and y. wang, 1991).  They are keyed to the "programmable
CC   expressions" in the notes of 1/23/91.
CC----------------------------------------------------------------
C      implicit none
CC Passed parameters
C      integer nsp
C      double precision rho1,rho2,drho,dabgrh,ddrho,drho1,drho2,
C     . dabgrh1,dabgrh2, ddrho1,ddrho2,vxc(2),exc
CC Local parameters
C      double precision pi,thrd,thrd2,d,zet,conrs,rs,fk,sk,g,y,t,uu,
C     . vv,ww,dvuaq,dvcup,h,haq,dvcdn,dvdaq
C      double precision rho, exl, ex1, ex2, vx1, vx2, polar, vx3
C      data pi,thrd,thrd2/3.14159265d0,0.333333333333d0,0.666666666667d0/
C
C      rho = rho1 + rho2
C      if (rho < 1d-18) then
C        vxc(1) = 0d0
C        vxc(2) = 0d0
C        exc = 0d0
C        return
C      endif
C      if (nsp == 1) then
C        call xperdew(rho,drho,ddrho,dabgrh,exl,vx1)
CC      ec = c*(drho*drho/(rho**fothr))*dexp(-phi)
CC      exc = ex+ec
CC      vcc = dexp(-phi)*drho**2/rho**fothr*(phi**2-phi-1d0)*dcdn
CC      vc1 = -2d0*c*dexp(-phi)*rho**(-oneth)*(ddrho*
CC     .      (1d0-phi/2d0)/rho - (twoth-phi*11d0/6d0+
CC     .      phi**2*7d0/12d0)*drho*drho/(rho*rho) + phi*
CC     .      (phi-3d0)*sig*dabgrh/(2d0*rho))
C        vxc(1) = vx1
C      else
C        call xperdew(rho,drho,ddrho,dabgrh,exl,vx3)
C        call xperdew(2*rho1,2*drho1,2*ddrho1,2*dabgrh1,ex1,vx1)
C        call xperdew(2*rho2,2*drho2,2*ddrho2,2*dabgrh2,ex2,vx2)
C        exl = (ex1+ex2)/2d0
CC      polar = (rho1-rho2)/rho
CC      if (polar > 1d0-1d-9) polar = 1d0-1d-9
CC      dd = (((1d0+polar)**fivth + (1d0-polar)**fivth)/2d0)**.5d0
CC      ec = c/dd*exp(-phi)*drho**2/rho**fothr
CC      exc = ex+ec
CC      vc1 = -c/rho**oneth*dexp(-phi)/dd*((2d0-phi)*ddrho/rho-
CC     .      (fothr-phi*11d0/3d0+phi**2*7d0/6d0)*drho*
CC     .      drho/(rho*rho) + phi*(phi -3d0)*sig*dabgrh/rho
CC     .      -5d0*rho**oneth/(6d0*dd*dd)*(rho1**
CC     .      twoth - rho2**twoth)/(rho**4d0)*drho*(2d0**twoth*
CC     .      (1d0-phi)*rho2*drho - 2d0**twoth*(2d0-phi)*
CC     .      rho*drho2))
CC      vc2 = -c/rho**oneth*dexp(-phi)/dd*((2d0-phi)*ddrho/rho-
CC     .      (fothr - phi*11d0/3d0 + phi**2*7d0/6d0)*drho*
CC     .      drho/(rho*rho) + phi*(phi -3d0)*sig*dabgrh/rho
CC     .      -5d0*rho**oneth/(6d0*dd*dd)*(rho2**
CC     .      twoth - rho1**twoth)/(rho**4d0)*drho*(2d0**twoth*
CC     .      (1d0-phi)*rho1*drho - 2d0**twoth*(2d0-phi)*
CC     .      rho*drho1))
CC      vcc = dexp(-phi)*drho**2/rho**fothr*(phi**2-phi-1d0)*dcdn/dd
C        vxc(1) = vx1
C        vxc(2) = vx2
C      endif
C
CC --- Nonlocal Correlation ---
C      D = rho1 + rho2
C      ZET = (rho1-rho2)/D
Cc     if (dabs(zet-1d0) < .000001) zet = sign(zet,1)*0.999999d0
C      if (d < 1d-18) return
C      CONRS = (3D0/(4D0*PI))**THRD
C      RS = CONRS/D**THRD
C      FK = 1.91915829D0/RS
C      SK = DSQRT(4D0*FK/PI)
C      G = ((1D0+ZET)**THRD2 + (1D0-ZET)**THRD2)/2D0
C      Y = 2D0*SK*G
C      T = dabs(drho)/(D*Y)
C      uu = drho*dabgrh/D/D/Y/Y/Y
C      vv = ddrho/D/Y/Y
C      ww = drho*(drho1-drho2)/D/Y/Y
Cc      print *,"on my way to corgga",rs,zet
C      CALL CORGGA(RS,ZET,T,UU,VV,WW,H,DVCUP,DVCDN,fk,sk,g)
Cc Change to Ry for LMTO package(CA: 23.04.91)
C      haq = -d*h*2d0
C      dvuaq = -dvcup*2d0
C      dvdaq = -dvcdn*2d0
Cc      print *,"after corgga; exl, vxcu, vxcd:",exl,vxc(1),vxc(2)
Cc      print *,"after corgga; h, dvcup, dvcdn:",haq,dvuaq,dvdaq
Cc     WT = 4.*pi*r*r*0.002
Cc     EXC =  D*H*WT  / 2D0
C      EXC = exl + 2d0*d*h
C      VXC(1) = vxc(1) + 2d0*dvcup
C      VXC(2) = vxc(2) + 2d0*dvcdn
Cc      EXC = (exl*0.5d0 - D*H)*2d0
Cc      VXC(1) = (vxc(1)*0.5d0 - DVCUP)*2d0
Cc      VXC(2) = (vxc(2)*0.5d0 - DVCDN)*2d0
Cc      print *, "ex, vxc1, vxc2", ex, vxc(1), vxc(2)
C      RETURN
C      END
C      SUBROUTINE CORGGA(RS,ZET,T,UU,VV,WW,H,DVCUP,DVCDN,fk,sk,g)
Cc  GGA91 CORRELATION
Cc  INPUT RS: SEITZ RADIUS
Cc  INPUT ZET: RELATIVE SPIN POLARIZATION
Cc  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
Cc  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
Cc  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
Cc  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
Cc  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
Cc  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
C      IMPLICIT double precision(A-H,O-Z)
Cc      COMMON/GAS/FK,SK,G,EC,ECRS,ECZET
C      DATA XNU,CC0,CX,ALF/15.75592D0,0.004235D0,-0.001667212D0,0.09D0/
C      DATA C1,C2,C3,C4/0.002568D0,0.023266D0,7.389D-6,8.723D0/
C      DATA C5,C6,A4/0.472D0,7.389D-2,100D0/
C      DATA THRDM,THRD2/-0.333333333333D0,0.666666666667D0/
C      BET = XNU*CC0
C      DELT = 2D0*ALF/BET
C      G3 = G**3
C      G4 = G3*G
C      call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
C      PON = -DELT*EC/(G3*BET)
C      B = DELT/(DEXP(PON)-1D0)
C      B2 = B*B
C      T2 = T*T
C      T4 = T2*T2
C      T6 = T4*T2
C      RS2 = RS*RS
C      RS3 = RS2*RS
C      Q4 = 1D0+B*T2
C      Q5 = 1D0 + B*T2+B2*T4
C      Q6 = C1+C2*RS + C3*RS2
C      Q7 = 1D0 + C4*RS + C5*RS2 + C6*RS3
C      CC = -CX + Q6/Q7
C      R0 = (SK/FK)**2
C      R1 = A4*R0*G4
C      COEFF = CC - CC0 - 3D0*CX/7D0
C      R2 = XNU*COEFF*G3
C      R3 = DEXP(-R1*T2)
C      H0 = G3*(BET/DELT)*DLOG(1D0+DELT*Q4*T2/Q5)
C      H1 = R3*R2*T2
C      H = H0+H1
Cc  LOCAL CORRELATION OPTION:
Cc     H = 0D0
Cc  ENERGY DONE. NOW THE POTENTIAL:
C      CCRS = (C2+2d0*C3*RS)/Q7 - Q6*(C4+2d0*C5*RS+3d0*C6*RS2)/Q7**2
C      RSTHRD = RS/3D0
C      R4 = RSTHRD*CCRS/COEFF
C      GZ = ((1D0+ZET)**THRDM - (1D0-ZET)**THRDM)/3D0
C      FAC = DELT/B + 1D0
C      BG = -3D0*B2*EC*FAC/(BET*G4)
C      BEC = B2*FAC/(BET*G3)
C      Q8 = Q5*Q5 + DELT*Q4*Q5*T2
C      Q9 = 1D0 + 2D0*B*T2
C      H0B = -BET*G3*B*T6*(2D0+B*T2)/Q8
C      H0RS = -RSTHRD*H0B*BEC*ECRS
C      FACT0 = 2D0*DELT-6D0*B
C      FACT1 = Q5*Q9 + Q4*Q9*Q9
C      H0BT = 2D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
C      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
C      H0Z = 3D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
C      H0T = 2d0*BET*G3*Q9/Q8
C      H0ZT = 3D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
C      FACT2 = Q4*Q5 + B*T2*(Q4*Q9+Q5)
C      FACT3 = 2D0*B*Q5*Q9 + DELT*FACT2
C      H0TT = 4D0*BET*G3*T*(2D0*B/Q8-(Q9*FACT3/Q8)/Q8)
C      H1RS = R3*R2*T2*(-R4+R1*T2/3D0)
C      FACT4 = 2D0 - R1*T2
C      H1RST = R3*R2*T2*(2D0*R4*(1D0-R1*T2) - THRD2*R1*T2*FACT4)
C      H1Z = GZ*R3*R2*T2*(3D0-4D0*R1*T2)/G
C      H1T = 2D0*R3*R2*(1D0-R1*T2)
C      H1ZT = 2D0*GZ*R3*R2*(3D0-11D0*R1*T2 + 4D0*R1*R1*T4)/G
C      H1TT = 4D0*R3*R2*R1*T*(-2D0+R1*T2)
C      HRS = H0RS+H1RS
C      HRST = H0RST+H1RST
C      HT = H0T+H1T
C      HTT = H0TT+H1TT
C      HZ = H0Z+H1Z
C      HZT = H0ZT+H1ZT
C      COMM = H + HRS+HRST + T2*HT/6D0 + 7D0*T2*T*HTT/6D0
C      PREF = HZ - GZ*T2*HT/G
C      FACT5 = GZ*(2D0*HT + T*HTT)/G
C      COMM = COMM - PREF*ZET - UU*HTT - VV*HT - WW*(HZT-FACT5)
C      DVCUP = COMM + PREF
C      DVCDN = COMM - PREF
Cc      print *,"h,dvcup,dvcdn:",h,dvcup,dvcdn
Cc  LOCAL CORRELATION OPTION:
Cc     DVCUP = 0D0
Cc     DVCDN = 0D0
C      RETURN
C      END
C      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
Cc  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
Cc  INPUT: SEITZ RADIUS(RS), RELATIVE SPIN POLARIZATION(ZET)
Cc  OUTPUT: CORRELATION ENERGY PER ELECTRON(EC), UP- AND DOWN-SPIN
Cc     POTENTIALS(VCUP,VCDN), DERIVATIVES OF EC WRT RS(ECRS) & ZET(ECZET)
Cc  OUTPUT: CORRELATION CONTRIBUTION(ALFC) TO THE SPIN STIFFNESS
C      IMPLICIT double precision(A-H,O-Z)
C      DATA GAM,FZZ/0.5198421D0,1.709921D0/
C      DATA THRD,THRD4/0.333333333333D0,1.333333333333D0/
C      F = ((1D0+ZET)**THRD4 + (1D0-ZET)**THRD4-2D0)/GAM
C      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
C     .    0.49294D0,1D0,RS,EU,EURS)
C      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
C     .    0.62517D0,1D0,RS,EP,EPRS)
C      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
C     .    0.49671D0,1D0,RS,ALFM,ALFRSM)
Cc  ALFM IS MINUS THE SPIN STIFFNESS ALFC
C      ALFC = -ALFM
C      Z4 = ZET**4
C      EC = EU*(1D0-F*Z4) + EP*F*Z4 - ALFM*F*(1D0-Z4)/FZZ
Cc  ENERGY DONE. NOW THE POTENTIAL:
C      ECRS = EURS*(1D0-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1D0-Z4)/FZZ
C      FZ = THRD4*((1D0+ZET)**THRD - (1D0-ZET)**THRD)/GAM
C      ECZET = 4D0*(ZET**3)*F*(EP-EU+ALFM/FZZ) + FZ*(Z4*EP - Z4*EU
C     .        -(1D0-Z4)*ALFM/FZZ)
C      COMM = EC -RS*ECRS/3D0 - ZET*ECZET
C      VCUP = COMM + ECZET
C      VCDN = COMM - ECZET
C      RETURN
C      END
C      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
Cc  CALLED BY SUBROUTINE CORLSD
C      implicit double precision (a-h,o-z)
C      P1 = P + 1D0
C      Q0 = -2D0*A*(1D0+A1*RS)
C      RS12 = DSQRT(RS)
C      RS32 = RS12**3
C      RSP = RS**P
C      Q1 = 2D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
C      Q2 = DLOG(1D0+1D0/Q1)
C      GG = Q0*Q2
C      Q3 = A*(B1/RS12+2D0*B2+3D0*B3*RS12+2D0*B4*P1*RSP)
C      GGRS = -2D0*A*A1*Q2 - Q0*Q3/(Q1**2+Q1)
C      RETURN
C      END
C      subroutine perdew(s,ff,f1,f2)
C      implicit none
C      double precision s,ff,f1,f2,ag1,ag2,g1,g1p,g1pp,g2,g2p,g2pp,div
C      double precision c1,f21,f22,beta,eta,gamma,alpha,delta
C      parameter (alpha=0.19645d0, beta=7.7956d0, gamma=0.2743d0,
C     .           delta=0.1508d0, eta=0.004d0)
C      ag1 = dsqrt((beta*s)**2+1d0)
C      ag2 = delta*dexp(-1d2*s*s)
C      g1 = alpha*dlog(beta*s+ag1)
C      g1p = alpha*beta/ag1
C      g1pp = -alpha*(beta**3)*s/(ag1**3)
C      g2 = gamma-ag2
C      g2p = 2d2*s*ag2
C      g2pp = 2d2*ag2*(1d0-2d2*s*s)
C      div = (1d0+(g1+eta*s*s*s)*s)
C      ff = (1d0+(g1+g2*s)*s)/div
C      f1 = ((((((eta*g2p)*s+eta*(g1p-2d0*g2))*s-3d0*eta*g1)*s+
C     .     (g1*g2p-4d0*eta-g1p*g2))*s+(g2p+g1*g2))*s+2d0*g2)*s / div**2
C      c1 = -2d0*((4d0*eta*s*s+g1p)*s+g1)
C      f21 = f1/s
C      f22 = ((((eta*g2pp*s+eta*(3d0*g2p+g1pp))*s+(eta*(g1p-
C     .     8d0*g2)))*s+(g1*g2pp-g1pp*g2-9d0*eta*g1))*s+
C     .     (3d0*g1*g2p-g1p*g2+g2pp-8d0*eta))*s+g1*g2+3d0*g2p
C      f2 = (f21*c1+f22/div)/div
C      end
C      subroutine xperdew(rho,drho,ddrho,dabgrh,ex,vx1)
C      implicit none
C      double precision rho,drho,ddrho,dabgrh,ex,vx1,s,t,u,ff,f1,f2,rs
C      double precision pi,ax,akf,fothr,oneth
C      pi = 4d0*datan(1d0)
C      fothr = 4d0/3d0
C      oneth = 1d0/3d0
C      rs = (4d0*pi*rho/3d0)**(-oneth)
C      akf = 1d0/(0.521062d0*rs)
C      ax = -1.4771175d0
C      s = dabs(drho)/(2d0*akf*rho)
C      t = ddrho/rho/(2d0*akf)**2
C      u = drho*dabgrh/rho**2/(2*akf)**3
C      call perdew(s,ff,f1,f2)
C      ex = ax*(rho)**fothr*(ff-1)
C      vx1 = ax*rho**(oneth)*(fothr*(ff-1) - t/s*f1 - (u-fothr*s**3)*f2)
C      end
C#else
      subroutine nlocxc(rhop,rhom,grho,abgrh,ggrho,grhop,grhom,
     .  abgrhp,abgrhm,ggrhop,ggrhom,nsp,vxc,exc)
C- Langreth-Mehl-Hu correction to the exc and vxc
C ---------------------------------------------------------------
Ci Inputs:
Ci   rhop,rhom (rho+,rho-); grho,grhop,grhom (grad rho)
Ci   ggrhop,ggrhom (Laplacians rho)
Ci   abgrh,abgrhp,abgrhm (grad abs grad rho)
Co Outputs:
Co   vxc, exc
Cr Remarks:
Cr   exc is the energy density time the charge
Cr   factor f is 1/6.  For pure gradient expansion put f=0 in data below
Cr   cutoff is a convergence factor which eliminates the blowup
Cr   of vxc for small rho
Cr   The Langreth Mehl form does not use abgrhp,abgrhm
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nsp
      double precision rhop,rhom,grho,grhop,grhom,ggrho,ggrhop,ggrhom,
     .  abgrh,abgrhp,abgrhm,vxc(2),exc
C Local parameters
      double precision pi,fivth,fothr,twoth,oneth,sevni,rho,
     .                 bigf,aa,dd,polar,f,sig,cutoff,cutof1,
     .                 cutof2,hh
      parameter (f=0.15d0,hh=1d-3,fivth=5d0/3d0,fothr=4d0/3d0,
     .  twoth=2d0/3d0, oneth=1d0/3d0, sevni=7d0/9d0)

      pi = 4d0*datan(1d0)
      vxc(1) = 0d0
      vxc(2) = 0d0
      exc = 0d0
      rho = rhop+rhom
      if (rho < 1d-18) return
      bigf = (9d0*pi)**(1d0/6d0)*f*dabs(grho)/(rho**(7d0/6d0))
      aa = (pi/(8d0*(3d0*pi*pi)**fothr))
      sig = dsign(1d0,grho)

      if (nsp == 1) then
      cutoff = dexp(-hh*(grho*grho/(rho*rho))*(rho/2)**(-twoth))
      exc = aa*(grho*grho/(rho**fothr))*(2d0*dexp(-bigf)-sevni)
c    .    *cutoff
      vxc(1) = 2d0*aa*rho**(-oneth)*(sevni*(ggrho/rho-twoth*
     .      grho*grho/(rho*rho))-2d0*dexp(-bigf)*(ggrho*
     .      (1d0-bigf/2d0)/rho-(twoth-bigf*11d0/6d0+
     .      bigf*bigf*7d0/12d0)*grho*grho/(rho*rho)+bigf*
     .      (bigf-3d0)*sig*abgrh/(2d0*rho)))*cutoff
      else
      polar = (rhop-rhom)/rho
      dd = (((1d0+polar)**fivth+(1d0-polar)**fivth)/2d0)**.5d0
      exc = aa*(-sevni/(2d0**oneth)*(grhop*grhop/(rhop**
     .   fothr)+grhom*grhom/(rhom**fothr))+2d0/dd*dexp(-bigf)*
     .   grho*grho/(rho**fothr))
      cutof1 = dexp(-hh*(grhop*grhop/(rhop*rhop))*rhop**(-twoth))
      vxc(1) = aa/rho**oneth*(-sevni/(2d0**oneth)*(fothr*
     .      grhop*grhop/(rhop*rhop)-2d0*ggrhop/rhop)*(rho/rhop)
     .      **oneth-2d0*dexp(-bigf)/dd*((2d0-bigf)*ggrho/rho-
     .      (fothr-bigf*11d0/3d0+bigf*bigf*7d0/6d0)*grho*
     .      grho/(rho*rho)+bigf*(bigf-3d0)*sig*abgrh/rho
     .      -5d0*rho**oneth/(6d0*dd*dd)*(rhop**
     .      twoth-rhom**twoth)/(rho**4d0)*grho*(2d0**twoth*
     .      (1d0-bigf)*rhom*grho-2d0**twoth*(2d0-bigf)*
     .      rho*grhom)))*cutof1
      cutof2 = dexp(-hh*(grhom*grhom/(rhom*rhom))*rhom**(-twoth))
      vxc(2) = aa/rho**oneth*(-sevni/(2d0**oneth)*(fothr*
     .      grhom*grhom/(rhom*rhom)-2d0*ggrhom/rhom)*(rho/rhom)
     .      **oneth-2d0*dexp(-bigf)/dd*((2d0-bigf)*ggrho/rho-
     .      (fothr-bigf*11d0/3d0+bigf*bigf*7d0/6d0)*grho*
     .      grho/(rho*rho)+bigf*(bigf-3d0)*sig*abgrh/rho
     .      -5d0*rho**oneth/(6d0*dd*dd)*(rhom**
     .      twoth-rhop**twoth)/(rho**4d0)*grho*(2d0**twoth*
     .      (1d0-bigf)*rhop*grho-2d0**twoth*(2d0-bigf)*
     .      rho*grhop)))*cutof2
      endif
      end
C#endif
      subroutine radgra(a,b,nr,r,v,gradv)
C- radial gradient
C ------------------------------------------------------------------
Ci Inputs:
Ci    a,b,nr,r(nr),v(nr)
Co Outputs:
Co    gradv(nr)
Cr Remarks:
Cr    makes the derivative of the function v defined in a mesh
Cr    of this kind:
Cr                r(i)=B (exp A(i-1)-1)
C ------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nr
      double precision a,b,r(nr),v(nr),gradv(nr)
C Local parameters
      integer nm2,i

C Forward diffs for first and second point (Handbook,25.3.9 with 25.1.1)
      gradv(1) = ((6d0*v(2)+20d0/3d0*v(4)+1.2d0*v(6))
     . -(2.45d0*v(1)+7.5d0*v(3)+3.75d0*v(5)+1d0/6d0*v(7)))/a
      gradv(2) = ((6d0*v(3)+20d0/3d0*v(5)+1.2d0*v(7))
     . -(2.45d0*v(2)+7.5d0*v(4)+3.75d0*v(6)+1d0/6d0*v(8)))/a

C Five points formula  (25.3.6)
      nm2 = nr-2
      do  i = 3, nm2
        gradv(i) = ((v(i-2)+8*v(i+1)) - (8*v(i-1)+v(i+2)))/12/a
      enddo

C Five points formula  (25.3.6)
      gradv(nr-1) = (-1d0/12d0*v(nr-4)+0.5d0*v(nr-3)-1.5d0*v(nr-2)
     .           +5d0/6d0*v(nr-1)+0.25d0*v(nr))/a
      gradv(nr) = (0.25d0*v(nr-4)-4d0/3d0*v(nr-3)+3d0*v(nr-2)
     .         -4d0*v(nr-1)+25d0/12d0*v(nr))/a
C Three points formula  (25.3.4)
C     gradv(nr-1)=(v(nr)-v(nr-2))/2d0/a
C     gradv(nr)=(v(nr-2)/2d0-2d0*v(nr-1)+1.5d0*v(nr))/a

      do  i = 1, nr
        gradv(i) = gradv(i)/(r(i)+b)
      enddo
      end
