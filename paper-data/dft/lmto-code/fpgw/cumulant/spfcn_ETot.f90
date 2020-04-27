! Calculate spectral function using the cumulant method given the self-energy
! and various quasiparticle properties.
! Written by Joshua Kas
! 10/05/2016
! Modified 11/23/2016 - Fixed segfault due to redefinition of nomega.
PROGRAM spfcn2
  USE parmod
  IMPLICIT NONE
  REAL(8) xk, xklast, tmp, tmp2, Ek, gamk, kFermi, norml, normg, norm, normqp, Mu, wp,slope,intercept, QPCorr, gamk_user, &
       & EFermi, deltagF, deltalF, domega, deltat, gam, gam0,gamGW, tmax, ETot, ETot2, ETotAH, kocc, EInt(2), EIntAH(2), EIntKin(2), E0, EKin, EKinGW, Eexch, Rs, &
       & omegaMax, omegaMin, ESpec, totocc, occint(2), ExInt(2), Ex, tol, ETotRPA, Rea, Ima, ETotGW, EIntGW(2), ntot, ntotGW, ImSigF, betap, betapp, ntot0, E0Test, MuGW, &
       & alpha, HFMult, Sigc, wmin, wmax, w0, dw0, dxHigh, dxLow, kT, deltac
  REAL(8),ALLOCATABLE :: t(:), betal(:), betag(:), betag2(:), omega(:), omega0(:), ReSig(:), ImSig(:), ReSigInt(:), ImSigInt(:), omegaOut(:), ERPA(:), kocc0(:), E0T(:), &
       &  dos(:), dosOcc(:), dosGW(:), dosqp(:), ampband(:), gamband(:), scale2(:), wtk(:)
  COMPLEX*16,ALLOCATABLE :: Cl(:), Cg(:), Cg2(:), CInt(:), CInt2(:),Gg(:), Gl(:), Gtg(:), Gtl(:), qpl(:), G_GW(:), G_SB(:), G0(:), G0t(:), EkSCF(:), deltal(:), deltag(:), deltaTot(:), koccArr(:), koccGW(:)
  INTEGER ik, nk, iomega,  nt, it, iSCF, nscf, kMax, iFermi, iFermiGW, itmp, ik2, ik3, nomega, nomega_fine, nomega0, iprint, ibroad, ioerr, nLow, nHigh, nband, iband
  REAL(8), parameter :: pi  = 3.1415926535897932384626433832795d0, ryd  = 13.605698, hart = 2.d0*ryd
  INTEGER, PARAMETER :: nomegamax = 1000000
  COMPLEX*16, PARAMETER :: coni = (0.d0,1.d0)
  COMPLEX*16 HFExc, al, ag, atot, Sk, ck, Ak, SkAk, Zk, ZkGW, ZkSC, DelStatic
  LOGICAL CTilde, GRet, Formal_a, Formal_del, MixedForm, PrintPercent, ImAxis, UserMu, G0Test, CalcSpfcn, AddHF, SelfConsistentZ, CalcExtInf, &
       & GWCorr
  CALL par_begin
  OPEN(UNIT=14,FILE='spfcn.in',STATUS='OLD')
  ! ibraod = 1, lorenzian broadening
  ! ibroad = 2, gaussian broadenin
  ibroad = 2
  Rea = 0.d0
  nband=1
  alpha = 1.d0
  GRet = .TRUE.
  Formal_a = .FALSE.
  Formal_del = .FALSE.
  MixedForm = .FALSE.
  ImAxis = .FALSE.
  SelfConsistentZ = .FALSE.

  !PRINT*, 'Enter number of time points: '
  READ(14,*) nomega
  !PRINT*, this_process, 'nomega'
  READ(14,*) nt
  !PRINT*, this_process, 'nt', nt
  !PRINT*, 'Enter broadening: '
  READ(14,*) gam0
  !PRINT*, this_process, 'gam0'
  !PRINT*, 'Enter Rs'
  !READ(14,*) Rs
  !PRINT*, this_process, 'Rs'
  !PRINT*, 'Enter wMin, wMax:'
  !READ(14,*) omegaMin, omegaMax
  omegaMin = -1.d10
  omegaMax = 1.d10
  !PRINT*, this_process, 'omegamin'
  !PRINT*, 'Use CTilde method? (.TRUE./.FALSE.)'
  !READ(14,*) CTilde
  CTilde = .FALSE.
  !PRINT*, this_process, 'ctilde'
  !PRINT*, 'Enter number of k points between 0 and kFermi, kMax:'
   READ(14,*) nk ! number of spectral functions, i.e., kpoints*bands
  !READ(14,*) nk, kMax
  !PRINT*, this_process, 'nk'
  !READ(14,*) CalcSpfcn
  CalcSpfcn = .TRUE.
  !PRINT*, this_process, 'CalcSpfcn'
  !READ(14,*) UserMu
  READ(14,*) Mu
  Mu = Mu/hart
  READ(14,*) kT
  ! GRet
  !READ(14,*) GRet
  GRet = .TRUE.
  !IF(.NOT.GRet) Formal_del = .TRUE.
  !READ(14,*) alpha
  alpha = 1.d0
  !READ(14,*) SelfConsistentZ
  SelfConsistentZ = .FALSE.
  READ(14,*) CalcExtInf
  READ(14,*) GRet
  READ(14,*) GWCorr
  READ(14,*) Formal_del
  READ(14,*) Formal_a
  IF(Formal_del) THEN
     READ(14,*) gamk_user
     gamk_user = gamk_user/27.2d0
  ELSE
     gamk_user = -1.d0
  END IF

  hfmult = 1.d0
  !read(14,*,iostat=ioerr) addhf
  addhf = .true.
  !if(ioerr == 0) then
  if(.not.addhf) hfmult = 0.d0
  !end if
  !if(.not.gret) then
  !   formal_a = .false.
  !   formal_del = .false.
  !   mixedform = .false.
  !end if
  !nk = nk*kmax
  tol = 0.d0

  allocate(cl(2*nt),cg(2*nt), cg2(2*nt),cint(max(nomegamax,2*nt)),cint2(max(nomegamax,2*nt)),t(2*nt),gg(nomegamax), qpl(nomegamax), g_gw(nomegamax), g_sb(nomegamax),scale2(nomegamax),wtk(nk))
  allocate(gl(nomegamax),gtg(2*nt),gtl(2*nt),g0(nomegamax),g0t(2*nt),ekscf(nk),deltal(nk),deltag(nk),deltatot(nk),erpa(nk),e0t(nk),koccarr(nk),koccgw(nk),kocc0(nk))
  allocate(omega(nomegamax), omega0(nomegamax), omegaout(3*(nomegamax)), betal(nomegamax), betag(nomegamax), betag2(nomegamax), resig(nomegamax),imsig(nomegamax),resigint(nomegamax),imsigint(nomegamax),dos(nomegamax),dosqp(nomegamax),dosocc(nomegamax),dosgw(nomegamax))
  allocate(ampband(nband),gamband(nband))
  ampband(:) = 1.d0
  !ampband(1) = 17.27d0
  !ampband(2) = 12.32d0
  !ampband(3) = 12.29d0
  !ampband(4) = 12.29d0
  !gamband(:) = 0.01d0
  !gamband(1) = 1.0d0
  !gamband(2) = 0.8d0
  !gamband(3) = 2.0d0
  !gamband(4) = 2.0d0

  open(unit=15,file='selfenergy.dat',status='old')
  open(unit=18,file='ek.dat',status='old')
  open(unit=19,file='QPCorr.dat',status='old')
  open(unit=20,file='wtk.dat',status='old')
  open(unit=21,file='deltac.dat',status='replace')

  !open(unit=20,file='wp.dat', status='old')
  if(master) then
     open(unit=26,file='nkrpa.dat', status='replace')
     !open(unit=27,file='etot.dat', status='replace')
     !open(unit=25,file='etot_k.dat', status='replace')
     !open(unit=20,file='wp.dat', status='old')
     !open(unit=18,file='se.dat', status='replace')
     open(unit=29,file='dos.dat', status='replace')
     if(calcspfcn) then
        !open(unit=21,file='spfcngw.dat',status='replace')
        !open(unit=28,file='ak.dat', status='replace')
        open(unit=16,file='beta.dat',status='replace')
        open(unit=17,file='spfcn.dat',status='replace')
        !open(unit=23,file='zk.dat', status='replace')
        !open(unit=24,file='sigonshell.dat', status='replace')
     else
        !open(unit=21,file='spfcngw.dat',status='old')
        !open(unit=28,file='ak.dat', status='old')
        open(unit=16,file='beta.dat',status='old')
        open(unit=17,file='spfcn.dat',status='old')
        !open(unit=23,file='zk.dat', status='old')
        !open(unit=24,file='sigonshell.dat', status='replace')
     end if
  end if
  !read(20,*) wp
  wp = 1.d0
  gam0 = gam0/hart
  ! find quantities at kfermi
  etot = 0.d0
  erpa = 0.d0
  etotah = 0.d0
  totocc = 0.d0
  occint = 0.d0
  espec = 0.d0
  ekin = 0.d0
  ekingw = 0.d0
  e0 = 0.d0
  dos=0.d0
  dosocc=0.d0
  dosgw=0.d0
  deltag = 0.d0

  ! for each k, make spectral function
  ! start k-loop
  itmp = 1
  xklast = 0.d0
  !print*, 'nt, nomega', nt, nomega
  do ik = 1, nk
     if(this_process > nk) then
        print*, 'more processors than k-points. reduce to ', nk
        call par_stop
     end if

     if(.not.(calcspfcn).and.(.not.master)) goto 15

     ! read imaginary self-energy im[sigma(k, ek+omega)] from file
     betal = 0.d0
     betag = 0.d0
     betag2 = 0.d0
     read(18,*) ek
     xk = SQRT(2.d0*ek)
     read(19,*) QPCorr
     read(20,*) wtk(ik)
     PRINT*, 'nomega = ', nomega
     do iomega = 1, nomega
        betag(iomega) = 0.d0
        betal(iomega) = 0.d0
        !read(15,*) omega(iomega), resig(iomega), imsig(iomega), ek, qpcorr, scale2(iomega), wtk(ik)
        read(15,*) omega(iomega), resig(iomega), imsig(iomega), scale2(iomega)
        if(.not.calcextinf) scale2(iomega) = 1.d0
        !if(.not.usermu) mu = mugw
        mugw = mu
        !if(((omega(iomega) - efermi)*imsig(iomega)) > 0.d0) then
        !   imsig(iomega) = 0.d0
        !end if
        !xk = xk + kfermi
        efermi = mugw
        !ck = xk
        !ekscf(ik) = xk**2/2.d0
        !ek = ekscf(ik)
        ekscf(ik) = ek

        omega(iomega) = omega(iomega) - dble(ek)

        ! get betag, betag2
        betag(iomega) = abs(imsig(iomega))/pi*scale2(iomega)
        betal(iomega) = betag(iomega)
        betag2(iomega) = betag(iomega)
        if(.not.gret) then
           if((ek < efermi).and.(omega(iomega) > (efermi-dble(ek)))) then
              betag(iomega) = 0.d0
           elseif((ek > efermi).and.(omega(iomega) < (efermi-dble(ek)))) then
              betag(iomega) = 0.d0
           end if
           if((ek < efermi).and.(omega(iomega) <= (efermi-dble(ek)))) then
              betal(iomega) = 0.d0
           elseif((ek > efermi).and.(omega(iomega) >= (efermi-dble(ek)))) then
              betal(iomega) = 0.d0
           end if
        end if

        !if(.not.gret) then
           ! use time-ordered single branch method
        !   if((omega(iomega) + ek) < efermi) then
        !      betag(iomega) = 0.d0
        !   else
        !      betal(iomega) = 0.d0
        !   end if
        !   if(xk < kfermi) then
        !      cint(iomega) = betag(iomega)/omega(iomega)
        !   else
        !      cint(iomega) = betal(iomega)/omega(iomega)
        !   end if
        !else
        !   cint(iomega) = betag(iomega)/omega(iomega)
        !end if
     end do
     ! get deltac
     !call terp(omega,resig,nomega,3,efermi-dble(ek),deltac)

     kfermi = 0.d0
     atot = -rea - coni*ima
     if(xk < kfermi) atot = conjg(atot) ! don't know why i need this, but it seems i do.
     zk = exp(-atot)
     zksc = 1.d0
     if(selfconsistentz) then
        zksc = solvezsc(atot)
        betag(:) = betag(:)*dble(zksc)
        betal(:) = betal(:)*dble(zksc)
        cint(:) = cint(:)*dble(zksc)
     else
        zksc = dble(zk)
     end if

     !print*, 'line 158'
     if(ik /= (itmp+this_process)) CYCLE ! mpi - each processor calculates one
     !print*, 'line 160'

     itmp = itmp + numprocs
     if(master) then
        print*
        print*, '########################################'
        print*, 'calculation # ', ik, "of", nk, 'process:', this_process
        print*, '########################################'
        print*
     end if
     domega = omega(2) - omega(1)
     deltal(ik) = 0.d0
     deltag(ik) = 0.d0
     if(xk < kfermi) then
        deltag(ik) = cintegrate(cint,nomega,domega)
     else
        deltal(ik) = cintegrate(cint,nomega,domega)
     end if

     !betap = (betag(nomega+2) - betag(nomega))/domega/2.d0
     !betapp = (betag(nomega+3) - 2.d0*betag(nomega+1) + betag(nomega-1))/4.d0/domega**2
     ! find im part of qp energy
     call terp(omega,betag,nomega,1,0.d0,gamk)
     call terp(omega,betag,nomega,1,0.00333333333d0,betap)
     call terp(omega,betag,nomega,1,-0.00333333333d0,tmp)
     !betap = (betap-gamk)/0.00333333333d0
     betap = 0.5d0*(betap-tmp)/0.00333333333d0
     DO iomega = 1, nomega
        !IF(.NOT.GRet) THEN
        ! Use time-ordered single branch method
        !   IF((omega(iomega) + Ek).LT.EFermi) THEN
        !      betag(iomega) = 0.d0
        !   ELSE
        !      betal(iomega) = 0.d0
        !   END IF
        !END IF
        ! Use separation of the quasiparticle shift
        IF(Formal_del) THEN
           betag(iomega) = betag(iomega) - gamk
           betal(iomega) = betal(iomega) - gamk
        END IF
        IF(gamk_user.GT.0.d0) gamk = gamk_user
        ! Use separation of the phase and amplitude reduction (alpha) - This doesn't work very well.
        IF(Formal_a) THEN
          betag(iomega) = betag(iomega) - omega(iomega)*betap
          betal(iomega) = betal(iomega) - omega(iomega)*betap
        END IF
     END DO
     ! make a new fine grid such that w = 0 is halfway through.
     domega = ABS(gam0)/5.d0
     nomega_fine = INT((omega(nomega)-omega(1))/domega)
     nomega_fine = nomega_fine/2
     nomega_fine = nomega_fine*2
     IF(nomega_fine.LT.nomega) nomega_fine = nomega
     domega = (omega(nomega) - omega(1))/nomega_fine
     IF(nomega_fine.GT.nomegamax) THEN
        nomega_fine = nomegamax
        PRINT*, 'nomegamax too small.'
        PRINT*, 'Please recompile with larger nomegamax.'
        CALL par_stop
     END IF

     DO iomega = 1, nomega_fine
        omega0(iomega) = domega*iomega
        omega0(iomega) = omega0(iomega) - domega*DBLE(nomega_fine-1)/2.d0
        IF((omega0(iomega).LT.omega(1)).OR.(omega0(iomega).GT.omega(nomega))) THEN
           betag2(iomega) = 0.d0
        ELSE
           CALL terp(omega,betag,nomega,1,omega0(iomega),betag2(iomega))
        END IF
        IF(.FALSE.) THEN !DEBUG
           !betag2(iomega) = 0.5d0
           !IF(omega0(iomega).GT.0.5d0) betag2(iomega) = 0.d0
           !IF(omega0(iomega).LT.-0.2d0) betag2(iomega) = 0.d0
           IF(ABS(omega0(iomega)).LT.0.1) betag2(iomega) = 0.d0
        END IF
     END DO
     betag(1:nomega_fine) = betag2(1:nomega_fine)
     DO iomega = 1, nomega_fine
        IF((omega0(iomega).LT.omega(1)).OR.(omega0(iomega).GT.omega(nomega))) THEN
           betag2(iomega) = 0.d0
        ELSE
           CALL terp(omega,betal,nomega,1,omega0(iomega),betag2(iomega))
        END IF
        IF(.FALSE.) THEN !DEBUG
           !betag2(iomega) = 0.5d0
           !IF(omega0(iomega).GT.0.5d0) betag2(iomega) = 0.d0
           !IF(omega0(iomega).LT.-0.2d0) betag2(iomega) = 0.d0
           IF(ABS(omega0(iomega)).LT.0.1) betag2(iomega) = 0.d0
        END IF
     END DO
     betal(1:nomega_fine) = betag2(1:nomega_fine)
     omega = omega0
     gam = gam0

     ! Form Cl and Cg by integrating beta
     Cg = 0.d0
     Cl = 0.d0


     ! Integrate to find deltaTot = correlation part of quasiparticle self-energy
     deltaTot(ik) = 0.d0
     deltal(ik) = 0.d0
     DO iomega = 1, nomega_fine - 1
        IF(omega(iomega+1).EQ.0.d0) THEN
           IF(iomega.LT.nomega_fine-1) THEN
              slope = (betag(iomega+2)-betag(iomega))/(omega(iomega+2)-omega(iomega))
              intercept = betag(iomega) - slope*omega(iomega)
              deltaTot(ik) = deltaTot(ik) + slope*(omega(iomega+2) - omega(iomega)) + intercept*LOG(ABS(omega(iomega+2)/omega(iomega)))
              slope = (betal(iomega+2)-betal(iomega))/(omega(iomega+2)-omega(iomega))
              intercept = betal(iomega) - slope*omega(iomega)
              deltal(ik) = deltal(ik) + slope*(omega(iomega+2) - omega(iomega)) + intercept*LOG(ABS(omega(iomega+2)/omega(iomega)))
           END IF
        ELSE IF(omega(iomega).NE.0.d0) THEN
           slope = (betag(iomega+1)-betag(iomega))/(omega(iomega+1)-omega(iomega))
           intercept = betag(iomega) - slope*omega(iomega)
           deltaTot(ik) = deltaTot(ik) + slope*(omega(iomega+1) - omega(iomega)) + intercept*LOG(ABS(omega(iomega+1)/omega(iomega)))
           slope = (betal(iomega+1)-betal(iomega))/(omega(iomega+1)-omega(iomega))
           intercept = betal(iomega) - slope*omega(iomega)
           deltal(ik) = deltal(ik) + slope*(omega(iomega+1) - omega(iomega)) + intercept*LOG(ABS(omega(iomega+1)/omega(iomega)))
        END IF
        CInt(iomega) = betag(iomega)/omega(iomega)
        WRITE(16,*) omega(iomega), betag(iomega), betal(iomega)
        !PRINT*, 'iomega, deltaTot', iomega, deltaTot(ik), omega(iomega)
     END DO
     WRITE(16,*)
     !deltaTot(ik) = CIntegrate(CInt,nomega,domega)
     !IF(GWCorr) THEN
     !   Sigc = DBLE(deltaTot(ik))
     !ELSE
     !   Sigc = 0.d0
     !END IF
     WRITE(21,*) ik, DBLE(deltaTot(ik)), DBLE(deltal(ik))
     deltac = DBLE(deltaTot(ik))
     IF(.NOT.GWCorr) THEN
        IF(GRet) THEN
           deltaTot(ik) = 0.d0
        ELSE
           deltaTot(ik) = -deltal(ik)
        END IF
     END IF
!     CALL terp(omega,betag,nomega,3,0.d0,tmp)
!     deltaTot(ik) = deltaTot(ik) + coni*tmp*pi
     !PRINT*, 'ik, deltaTot = ', ik, deltaTot(ik)
     !IF(ABS(DIMAG(deltaTot(ik))).LT.gam) gam0 = gam !- ABS(DIMAG(deltaTot(ik)))
     !gam0 = gamband(mod(ik-1,nband)+1)
     !PRINT*, 'gamband, iband =', gam0, mod(ik-1,nband)+1, ik
     ! Set tmax based on extra broadening parameter, or total broadening for this k.
     tmax = 10.d0/(gam0) !+ ABS(DIMAG(deltaTot(ik))))
     deltat = tmax/DBLE(2*nt-1)

     ZkGW = 1.d0/(1.d0+atot)

     ! Get aless and agreater
     ag = 0.d0
     al = 0.d0
     DO iomega = 2, nomega
        IF((omega(iomega)+DBLE(EKSCF(ik))).GT.EFermi) ag = ag + 0.5d0*(ABS(ImSig(iomega))/omega(iomega)**2 + ABS(ImSig(iomega-1))/(omega(iomega-1)**2))*domega/pi
        IF((omega(iomega)+DBLE(EKSCF(ik))).LT.EFermi) al = al + 0.5d0*(ABS(ImSig(iomega))/omega(iomega)**2 + ABS(ImSig(iomega-1))/(omega(iomega-1)**2))*domega/pi
     END DO
     !IF(.NOT.GRet) THEN
     !   IF(Ek.LT.Mu) betag = betal
     !END IF
     IF(.NOT.CalcSpfcn) GOTO 15
     IF(master) PRINT*, "Forming C(t)."
     t = 0.d0
     CInt(1) = 0.d0
     DO it = 1, 2*nt
        t(it) = (0.5d0 + DBLE(it-1))*deltat
        IF(master) THEN
           IF(MOD(it*20,2*nt).EQ.0) THEN
              WRITE(*,'(i3,1a)') it*100/(2*nt),'%'
           END IF
           !IF(MOD(DBLE(it*100)/DBLE(2*nt),50.d0).EQ.0.d0) PRINT*, (it*100)/nt/2, "%"
        END IF
        DO iomega = 1, nomega_fine
           IF(t(it).EQ.0.d0) THEN
              CInt(iomega) = 0.d0
           ELSEIF(ABS(omega(iomega)*t(it)).LT.0.001d0) THEN
              IF(Formal_a.AND.Formal_del) THEN
                 CInt(iomega) = betapp*EXP(-coni*omega(iomega)*t(it))/2.d0
              ELSEIF(Formal_del) THEN

                 CInt(iomega) = -coni*betap*t(it)
                 !CInt(iomega) = -betag(iomega)*t(it)**2/2.d0
              ELSE
                 CInt(iomega) = -betag(iomega)*t(it)**2/2.d0
              END IF
              CInt2(iomega) = -betag2(iomega)*t(it)**2/2.d0
           ELSE
              IF(Formal_a.AND.Formal_del) THEN
                 CInt(iomega) = betag(iomega) * EXP(-coni*omega(iomega) * t(it))/(omega(iomega)-coni*tol)**2
              ELSEIF(Formal_del) THEN
                 CInt(iomega) = betag(iomega) * (EXP(-coni*omega(iomega) * t(it)) + coni*omega(iomega)*t(it) - 1.d0)/(omega(iomega)-coni*tol)**2
              ELSE
                 CInt(iomega) = betag(iomega) * (EXP(-coni*omega(iomega) * t(it)) + coni*omega(iomega)*t(it) - 1.d0)/(omega(iomega)-coni*tol)**2
              END IF
              CInt2(iomega) = betag2(iomega) * (EXP(-coni*omega(iomega) * t(it)) + coni*omega(iomega)*t(it) - 1.d0)/(omega(iomega)-coni*tol)**2
           END IF
        END DO
        IF(Formal_a.AND.Formal_del) THEN
           Cg(it) = CIntegrate(CInt,nomega_fine,domega) - atot
        ELSEIF(Formal_del) THEN
           Cg(it) = CIntegrate(CInt,nomega_fine,domega)
        ELSE
           Cg(it) = CIntegrate(CInt,nomega_fine,domega)
        END IF
        Cg2(it) = CIntegrate(CInt2,nomega_fine,domega)
     END DO

     ! Create omega0 grid.
     wmin = omega(1) + Ek
     wmax = omega(nomega_fine) + Ek
     w0 = DBLE(Ek)
     dw0=gam0/10.d0
     nLow=INT(nomega_fine*ABS(wmin-DBLE(Ek))/(wmax-wmin))
     nHigh=INT(nomega_fine*ABS(wmax-DBLE(Ek))/(wmax-wmin))
     dxHigh=LOG((ABS(wmax-DBLE(Ek))+dw0)/dw0)/DBLE(nHigh-1)
     dxLow=LOG((ABS(wmin-DBLE(Ek))+dw0)/dw0)/DBLE(nLow-1)
     !WRITE(99,*) '# wmin = ', wmin
     !WRITE(99,*) '# wmax = ', wmax
     !WRITE(99,*) '# w0 = ', w0
     !WRITE(99,*) '# dw0 = ', dw0
     !WRITE(99,*) '# nLow = ', nLow
     !WRITE(99,*) '# nHigh = ', nHigh
     !WRITE(99,*) '# dxHigh = ', dxHigh
     !WRITE(99,*) '# dxLow = ', dxLow
     DO iomega = nlow, 1, -1
        ! Grid that goes from wmin to Ek with a final step of dw0=gam0/10 and using nLow points
        omega0(nlow-iomega+1) = w0 - dw0*EXP(dxLow*(iomega-1))
     END DO
     DO iomega = nlow + 1, nomega_fine
        ! Grid that goes from Ek to wmax with an initial step of dw0=gam0/10 and using nHigh points
        omega0(iomega) = w0 + dw0*EXP(dxHigh*(iomega-nlow-1))
     END DO

     ! Interpolate self-energy for this k point.
     omega=omega+Ek
     DO iomega = 1, nomega
        CALL terp(omega,ReSig,nomega,1,omega0(iomega),ReSigInt(iomega))
        CALL terp(omega,ImSig,nomega,1,omega0(iomega),ImSigInt(iomega))
     END DO
     !DO iomega = 1, nomega
     !   WRITE(99,*) omega0(iomega), 1.d0/((omega0(iomega)-DBLE(Ek))**2 + 1.d0)
     !END DO

     ! Fourier transform G(t) to get G(omega), A(omega)
     IF(master) PRINT*, 'Fourier transforming G(t)'
     iprint = 0
     PrintPercent = .TRUE.
     DO iomega = 1, nomega_fine
        IF(master) THEN
           iprint = iprint + 1
           IF(MOD(iomega*20,nomega_fine).EQ.0) THEN
              WRITE(*,'(i3,1a)') iomega*100/nomega_fine,'%'
           END IF
           IF(iprint*10.GT.nomega_fine) THEN
              !PRINT*, iomega*100/(nomega), "%"
              iprint = 0
           END IF
        END IF
        !IF(ImAxis) THEN
        !   omega0(iomega) = Mu
        !   gam0 = gam*iomega
        !ELSE
        !   omega0(iomega) = omega(iomega) + DBLE(Ek)
        !END IF
        CInt = 0.d0
        DO it = 1, 2*nt
           Gtg(it) = 0.d0
           Gtl(it) = 0.d0
           CInt(it) = 0.d0
           G0t(it) = -coni*EXP(coni*(omega0(iomega)-Ek)*t(it)-(gam0*t(it))**2)
           !IF(SelfConsistentZ) DelStatic = DelStatic - (1.d0 - ZkSC)*(deltaTot(ik) - HFExc(ck,kFermi**2/2.d0,kFermi))
           IF(Formal_del) THEN
              !Gtg(it) = -coni*EXP(coni*(omega0(iomega) - DBLE(EkSCF(ik)))*t(it) - (gam0*Abs(t(it)))**ibroad + Cg(it))
              ! Debug
              !Ek=0.d0
              !gamk=0.d0
              Gtg(it) = -coni*EXP(coni*(omega0(iomega) - Ek - DBLE(deltaTot(ik)) - QPCorr)*t(it) - Abs(gamk*t(it)) - (gam0*Abs(t(it)))**ibroad + Cg(it))
              Gtl(it) = -coni*EXP(coni*(omega0(iomega) - Ek)*t(it) - Abs(gamk*t(it)) - (gam0*Abs(t(it)))**ibroad + Cg2(it))
           ELSE
              Gtg(it) = -coni*EXP(coni*(omega0(iomega) - Ek - DBLE(deltaTot(ik)) - QPCorr)*t(it) - (gam0*Abs(t(it)))**ibroad + Cg(it))
              ! DEBUG below
              !Gtg(it) = -coni*EXP(coni*(omega0(iomega))*t(it) - (gam0*Abs(t(it)))**ibroad + Cg(it))
              !Gtg(it) = -coni*EXP(coni*(omega0(iomega) - Ek)*t(it) - (gam0*Abs(t(it)))**ibroad)
              !IF(.NOT.GRet) THEN
              !   IF(xk.LT.kFermi) THEN
              !      Gtg(it) = Gtg(it)*EXP(coni*deltag(ik)*t(it))
              !   ELSE
              !      Gtg(it) = Gtg(it)*EXP(coni*deltal(ik)*t(it))
              !   END IF
              !END IF
           END IF

           CInt(it)  = -coni*EXP(coni*(omega0(iomega) - Ek + deltac - QPCorr)*t(it) - ABS(DIMAG(deltaTot(ik)))*Abs(t(it)) - (gam0*ABS(t(it)))**ibroad)
           ! Debug below
           !CInt(it)  = -coni*EXP(coni*(omega0(iomega)+DBLE(deltaTot(ik)))*t(it) - ABS(DIMAG(deltaTot(ik)))*Abs(t(it)) - (gam0*ABS(t(it)))**ibroad)
           !IF((DBLE(it-1-nt)+0.5d0).GT.0.d0) THEN
           !   CInt(it)  = -coni*EXP(coni*(omega0(iomega) - DBLE(EkSCF(ik) + deltaTot(ik)))*(DBLE(it-1-nt)+0.5d0)*deltat - (gam0+ABS(DIMAG(deltaTot(ik))))*Abs((DBLE(it-1-nt)+0.5d0)*deltat))
           !END IF
        END DO
        Gl(iomega) = -SUM(Gtl)*deltat !CIntegrate(Gtl,2*nt,deltat)*SIGN(1.d0,kFermi-xk) + 0.25d0*(3.d0*Gtl(1)+Gtg(2))*deltat
        Gg(iomega) = -SUM(Gtg)*deltat !CIntegrate(Gtg,2*nt,deltat)*SIGN(1.d0,kFermi-xk) + 0.25d0*(3.d0*Gtg(1)+Gtg(2))*deltat
        qpl(iomega)= -SUM(CInt)*deltat !CIntegrate(CInt(nt-3),nt,deltat)
        G0(iomega) = -SUM(G0t)*deltat
        !IF(Formal_del.AND.(ABS(omega(iomega)/wp).GT.4.d0)) Gg(iomega) = Gl(iomega)
     END DO

     ! Form GW Green's function.
     gamGW = gam0 !1.35d0
     !gamGW = 0.d0
     !IF(ABS(DIMAG(deltatot(ik)).GT.gam0)) gamGW = 0.d0
     DO iomega = 1, nomega
        G_GW(iomega) = 1.d0/(omega0(iomega) - xk - (ReSigInt(iomega)-gamk) - coni*(ABS(ImSigInt(iomega)) + gamGW))
     END DO

15   CONTINUE
     iFermi = 0.d0
     kocc = 0.d0
     IF((.NOT.master).AND.CalcSpfcn) THEN
        ! Send data to master
        PRINT*, 'Process:', this_process, 'sending data for ik = ', ik
        CALL par_send_int_scalar(ik,1,0,this_process)
        CALL par_send_dbl(omega0,nomega_fine,0,this_process)
        CALL par_send_dbl(xk,1,0,this_process)
        CALL par_send_dbl(deltac,1,0,this_process)
        CALL par_send_dc(EkSCF(ik),1,0,this_process)
        CALL par_send_dbl(wtk(ik),1,0,this_process)
        CALL par_send_dc(Gg,nomega_fine,0,this_process)
        CALL par_send_dc(Gl,nomega_fine,0,this_process)
        CALL par_send_dc(G0,nomega_fine,0,this_process)
        CALL par_send_dc(G_GW,nomega_fine,0,this_process)
        CALL par_send_dc(Zk,1,0,this_process)
        CALL par_send_dc(ZkGW,1,0,this_process)
        CALL par_send_dc(ZkSC,1,0,this_process)
        CALL par_send_dc(qpl,nomega_fine,0,this_process)
        CALL par_send_dbl(ReSig,nomega,0,this_process)
        CALL par_send_dbl(ImSig,nomega,0,this_process)
        CALL par_send_dc(atot,1,0,this_process)
        CALL par_send_dc(deltag(ik),1,0,this_process)
        CALL par_send_dc(deltal(ik),1,0,this_process)
        CALL par_send_dc(al,1,0,this_process)
        CALL par_send_dc(ag,1,0,this_process)
     ELSE
        ! Now loop all processors, receive data, and do integrals for nk, Etot
        ! etc.
        IF(ik.EQ.1) THEN
           EIntAH(:) = 0.d0
           EIntGW(:) = 0.d0
           ETotAH = 0.d0
           ETotGW = 0.d0
           Ex = 0.d0
           E0 = 0.d0
           ntot = 0.d0
        END IF

        !IF(ik.EQ.1) PRINT*, 'Calculating total energies and occupations.'
        DO ik2 = ik, MIN(ik + numprocs - 1, nk)
           iband = MOD(ik2-1,nband) + 1
           IF(MOD(ik2,10).EQ.0) PRINT '(A,I3,A)', '   ', INT(DBLE(ik2)/DBLE(nk)*100.d0), '%'
           IF((ik2.NE.ik).AND.CalcSpfcn) THEN
              CALL par_recv_int_scalar(ik3,1,ik2-ik,ik2-ik)
              IF(ik2.NE.ik3) THEN
                 PRINT*, 'k-points not equal. Stopping.'
                 CALL par_stop
              END IF
              CALL par_recv_dbl(omega0,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dbl(xk,1,ik2-ik,ik2-ik)
              CALL par_recv_dbl(deltac,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(EkSCF(ik2),1,ik2-ik,ik2-ik)
              CALL par_recv_dbl(wtk(ik2),1,ik2-ik,ik2-ik)
              CALL par_recv_dc(Gg,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dc(Gl,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dc(G0,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dc(G_GW,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dc(Zk,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(ZkGW,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(ZkSC,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(qpl,nomega_fine,ik2-ik,ik2-ik)
              CALL par_recv_dbl(ReSig,nomega,ik2-ik,ik2-ik)
              CALL par_recv_dbl(ImSig,nomega,ik2-ik,ik2-ik)
              CALL par_recv_dc(atot,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(deltag(ik2),1,ik2-ik,ik2-ik)
              CALL par_recv_dc(deltal(ik2),1,ik2-ik,ik2-ik)
              CALL par_recv_dc(al,1,ik2-ik,ik2-ik)
              CALL par_recv_dc(ag,1,ik2-ik,ik2-ik)
              PRINT*, 'Data received from process', ik2-ik
           END IF
           IF(.NOT.CalcSpfcn) THEN
              ! Read data from files
              DO iomega = 1, nomega
                 ! spfcn.dat
                 READ(17,*) omega0(iomega), xk, tmp, tmp, tmp2
                 Gg(iomega) = (tmp + coni*tmp2*pi)/wp
                 omega0(iomega) = omega0(iomega)*wp
                 ! spfcnGW.dat
                 !READ(21,*) omega(iomega), xk, tmp, tmp, tmp2
                 G_GW(iomega) = (tmp + coni*tmp2*pi)/wp
                 omega(iomega) = omega(iomega)*wp
              END DO
              ! ak.dat
              !READ(28,*) xk, kFermi, al, ag
              ! Zk.dat
              !READ(23,*) xk, kFermi, Zk, ZkSC, ZkGW
              EFermi = kFermi**2/2.d0
              !WRITE(77, '(20e20.10)') xk, kFermi, DBLE(deltag(ik2)), DIMAG(deltag(ik2)), DBLE(deltal(ik2)), DIMAG(deltal(ik2))
           ELSE
              ! Write to ak.dat
              !IF(xk.LT.kFermi) THEN
                 !WRITE(28,'(20e20.10)') xk, kFermi, DBLE(al), DIMAG(al), DBLE(ag), DIMAG(ag), DBLE(EXP(-ag)), DBLE(atot), DIMAG(atot)
              !ELSE
                 !WRITE(28,'(20e20.10)') xk, kFermi, DBLE(al), DIMAG(al), DBLE(ag), DIMAG(ag), 1.d0 - DBLE(EXP(-al)), DBLE(atot), DIMAG(atot)
              !END IF
              ! Write to deltak.dat
              !WRITE(77, '(20e20.10)') xk, kFermi, DBLE(deltag(ik2)), DIMAG(deltag(ik2)), DBLE(deltal(ik2)), DIMAG(deltal(ik2))
              ! Write to Zk.dat
              !WRITE(23,'(20e20.10)') xk, kFermi, DBLE(Zk), DIMAG(Zk), DBLE(ZkSC), DIMAG(ZkSC), DBLE(ZkGW), DIMAG(ZkGW), DBLE(EXP(-atot*ZkSC)),DIMAG(EXP(-atot*ZkSC))

              ! Write to spfcn.dat, spfcnGW.dat
              WRITE(17,*) '# omega, Ek, A_RC, A_QP'
              ZkSC = 1.d0
              DO iomega = 1, nomega_fine
                 IF(ImAxis) gam0=gam*iomega
                 !IF((omega(iomega)/wp.GT.omegaMin).AND.(omega(iomega)/wp.LT.omegaMax)) THEN
                    WRITE(17,'(20e20.10)') (omega0(iomega)), DBLE(EkSCF(ik2)), DIMAG(Gg(iomega))/pi, (DIMAG(ZkSC*qpl(iomega))/pi)
                    !WRITE(21,'(20e20.10)') omega0(iomega)/wp, xk, DBLE(EkSCF(ik2))/wp, DBLE(G_GW(iomega)*wp), DIMAG(G_GW(iomega)*wp)/pi, ReSig(iomega), ImSig(iomega)
                 !END IF
              END DO
              WRITE(17,*)
              !WRITE(21,*)
           END IF

           !DO iomega = 1, nomega
           !   WRITE(18,'(20e20.10)') omega0(iomega), DBLE(1.d0/Gg(iomega)-HFMult*DelStatic+omega0(iomega)), DIMAG(- 1.d0/Gg(iomega))
           !END DO
           !WRITE(18,*)
           ! Integrate total energies and occupations
           koccArr(ik2) = 0.d0
           koccGW(ik2)  = 0.d0
           kocc0(ik2)   = 0.d0
           ERPA(ik2)    = 0.d0
           E0T(ik2)     = 0.d0
           norm         = 0.d0
           EIntGW(1)    = 0.d0
           EIntAH(1)    = 0.d0
           !domega = omega0(2) - omega0(1)
           DO iomega = 2, nomega_fine
              domega = omega0(iomega) - omega0(iomega-1)
              norm = norm + 0.5d0*(ABS(DIMAG(Gg(iomega-1) + Gg(iomega))))*domega/pi
              IF(omega0(iomega).GT.omegaMin) THEN
                 IF(GRet) THEN
                    IF(omega0(iomega).LT.Mu) THEN
                       ERPA(ik2) = ERPA(ik2) + 0.5d0*DIMAG(Gg(iomega-1)*omega0(iomega-1) + Gg(iomega)*omega0(iomega))*domega/pi
                       koccArr(ik2) = koccArr(ik2) + 0.5d0*(ABS(DIMAG(Gg(iomega-1) + Gg(iomega))))*domega/pi
                       iFermi = iomega
                    END IF
                 ELSE
                    IF(xk.LE.kFermi) THEN
                       ERPA(ik2) = ERPA(ik2) + 0.5d0*DIMAG(Gg(iomega-1)*omega0(iomega-1) + Gg(iomega)*omega0(iomega))*domega/pi
                       koccArr(ik2) = koccArr(ik2) + 0.5d0*(ABS(DIMAG(Gg(iomega-1) + Gg(iomega))))*domega/pi
                       iFermi = iomega
                    ELSE
                       ERPA(ik2) = 0.d0
                    END IF
                 END IF
                 IF(omega0(iomega).LT.MuGW) THEN
                    EIntGW(1) = EIntGW(1) + 0.5d0*DIMAG(G_GW(iomega-1)*(omega0(iomega-1)-EFermi+MuGW) + G_GW(iomega)*(omega0(iomega)-EFermi+MuGW))*domega/pi
                    koccGW(ik2) = koccGW(ik2) + 0.5d0*(ABS(DIMAG(G_GW(iomega-1) + G_GW(iomega))))*domega/pi
                    E0T(ik2) = E0T(ik2) + 0.5d0*DIMAG(G0(iomega-1)*omega0(iomega-1) + G0(iomega)*omega0(iomega))*domega/pi
                    kocc0(ik2) = kocc0(ik2) + 0.5d0*(ABS(DIMAG(G0(iomega-1) + G0(iomega))))*domega/pi
                    iFermiGW = iomega
                 END IF
              END IF
           END DO

           ! Fix end points of integration.
           IF(GRet) THEN
              IF(omega0(1).LT.Mu.AND.iFermi.GT.0.AND.iFermi.LT.nomega_fine) THEN
                 CALL terpc(omega0,Gg,nomega_fine,1,Mu,CInt(1))
                 koccArr(ik2) = koccArr(ik2) + 0.5d0*(ABS(DIMAG(CInt(1))) + ABS(DIMAG(Gg(iFermi))))*(Mu-omega0(iFermi))/pi
                 ERPA(ik2) = ERPA(ik2) + 0.5d0*(ABS(DIMAG(CInt(1)))*Mu + ABS(DIMAG(Gg(iFermi)))*omega0(iFermi))*(Mu-omega0(iFermi))/pi
              END IF
           END IF
           IF(omega0(1).LT.EFermi) THEN
              CALL terpc(omega0,G_GW,nomega_fine,1,EFermi,CInt(1))
              koccGW(ik2) = koccGW(ik2) + 0.5d0*(ABS(DIMAG(CInt(1))) + ABS(DIMAG(G_GW(iFermiGW))))*(EFermi-omega0(iFermiGW))/pi
              EIntGW(1) = EIntGW(1) + 0.5d0*(ABS(DIMAG(CInt(1)))*MuGW + ABS(DIMAG(G_GW(iFermiGW)))*(omega0(iFermiGW)-EFermi+MuGW))*(EFermi-omega0(iFermiGW))/pi
              CALL terpc(omega0,G0,nomega_fine,1,EFermi,CInt(1))
              kocc0(ik2) = kocc0(ik2) + 0.5d0*(ABS(DIMAG(CInt(1))) + ABS(DIMAG(G0(iFermiGW))))*(EFermi-omega0(iFermiGW))/pi
              E0T(ik2) = E0T(ik2) + 0.5d0*(ABS(DIMAG(CInt(1)))*EFermi + ABS(DIMAG(G0(iFermiGW)))*omega0(iFermiGW))*(EFermi-omega0(iFermiGW))/pi
           END IF
           EIntGW(1) = EIntGW(1) + xk**2/2.d0*koccGW(ik2)
           E0T(ik2) = E0T(ik2) + xk**2/2.d0*kocc0(ik2)

           IF(xk.LT.kFermi) THEN
              E0 = E0 + 0.25d0*(xk**4+xklast**4)*(xk - xklast)
              CInt(1) = xklast
              CInt(2) = xk
              !Ex = Ex + 0.5d0*DBLE(HFMult*HFExc(CInt(1),EFermi,kFermi)*xklast**2 + HFMult*HFExc(CInt(2),EFermi,kFermi)*xk**2)*(xk-xklast)
           END IF
           IF(xk.LT.kFermi) THEN
              CInt(1) = xk
              !EIntAH(1) = xk**2 + DBLE(HFMult*HFExc(CInt(1),EFermi,kFermi) - deltag(ik2))
           ELSE
              EIntAH(1) = 0.d0
           END IF
           ETotAH = ETotAH + 0.5d0*(EIntAH(1)*xk**2 + EIntAH(2)*xklast**2)*(xk - xklast)
           ETotGW = ETotGW + 0.5d0*(EIntGW(1)*xk**2 + EIntGW(2)*xklast**2)*(xk - xklast)
           EIntAH(2) = EIntAH(1)
           EIntGW(2) = EIntGW(1)
           IF(ik2.EQ.1) THEN
              ntot = ntot + koccArr(ik2)*xk**3/3.d0
              EKin = EKin + koccArr(ik2)*xk**5/10.d0
              ntotGW = ntotGW + koccGW(ik2)*xk**3/3.d0
              EKinGW = EKinGW + koccGW(ik2)*xk**5/10.d0
           ELSE
              ntot = ntot + 0.5d0*(koccArr(ik2)*xk**2 + koccArr(ik2-1)*xklast**2)*(xk-xklast)
              ntotGW = ntotGW + 0.5d0*(koccGW(ik2)*xk**2 + koccGW(ik2-1)*xklast**2)*(xk-xklast)
              EKin = EKin + 0.5d0*(koccArr(ik2)*xk**4/2.d0 + koccArr(ik2-1)*xklast**4/2.d0)*(xk-xklast)
              EKinGW = EKinGW + 0.5d0*(koccGW(ik2)*xk**4/2.d0 + koccArr(ik2-1)*xklast**4/2.d0)*(xk-xklast)
           END IF

           !IF(xk/kFermi.GT.1.8d0) THEN
              ! use high k limit of occupation.
              !koccArr(ik2) = DBLE(al)
              !ERPA(ik2) = deltal(ik2) + DBLE(al)*xk**2/2.d0
           !END IF

           !IF(.FALSE.) THEN
              ! Fix endpoint by interpolating onto fermi level
           !   CALL terpc(omega0,Gg,nomega,1,Mu,CInt(1))
           !   ERPA(ik2) = (ERPA(ik2) +  0.5d0*DIMAG(Gg(iFermi)*omega0(iFermi) + CInt(1)*Mu)*(Mu - omega0(iFermi)))*xk**2 + xk**4/2.d0
              !koccArr(ik2) = koccArr(ik2) +  0.5d0*ABS(DIMAG(Gg(iFermi) + CInt(1)))*(Mu-omega0(iFermi))
           !END IF


           ! Write to nkRPA.dat
           WRITE(26,'(1e20.10,I4,20e20.10)') DBLE(EkSCF(ik2)), iband, MuGW, DBLE(koccArr(ik2)), DBLE(koccGW(ik2)), kocc0(ik2), norm

           ERPA(ik2) = ERPA(ik2) + koccArr(ik2)*xk**2/2.d0
           IF(ik2.GT.1) THEN
              ETotRPA = ETotRPA + 0.5d0*(ERPA(ik2-1)*xklast**2 + ERPA(ik2)*xk**2)*(xk-xklast)
              E0Test = E0Test + 0.5d0*(E0T(ik2-1)*xklast**2 + E0T(ik2)*xk**2)*(xk-xklast)
           END IF

           ! Write to ETot_k.dat
           !WRITE(25,'(20e20.10)') xk, kFermi, Mu, ERPA(ik2)*3.d0/kFermi**3/2.d0, EIntGW(1), EIntAH(1), ETotRPA*3.d0/kfermi**3/2.d0, ETotGW*3.d0/kFErmi**3/2.d0, ETotAH*3.d0/kFermi**3/2.d0

           ! Get dos(iomega) for this k (sum over bands)
           domega = 5.d0/(nomega_fine/2.d0-1.d0)
           IF(wtk(ik2).LT.0.d0) wtk(ik2) = xk**2*(xk-xklast)
           PRINT*, wtk(ik2), xk, xklast
           DO iomega = 1, nomega_fine
              tmp = (iomega-nomega_fine/2)*domega
              IF((tmp.GT.omega0(1)).AND.(tmp.LT.omega0(nomega_fine))) THEN
                 CALL terpc(omega0,Gg,nomega_fine,3,tmp,CInt(1))
                 dos(iomega) = dos(iomega) + DIMAG(CInt(1))/pi*wtk(ik2)
                 CALL terpc(omega0,qpl,nomega_fine,3,tmp,CInt(1))
                 dosqp(iomega) = dosqp(iomega) + DIMAG(CInt(1))/pi*wtk(ik2)
                 IF(DBLE(EkSCF(ik2)).LE.EFermi) dosOcc(iomega) = dosOcc(iomega) + DIMAG(CInt(1))/pi*wtk(ik2)
                 CALL terpc(omega0,G_GW,nomega_fine,1,tmp,CInt(1))
                 dosGW(iomega) = dosGW(iomega) + DIMAG(CInt(1))/pi*ampband(iband)
              END IF
           END DO

           ! If done with all k, print dos
           IF(ik2.EQ.nk) THEN
              tmp=0.d0
              DO iomega = 1, nomega_fine
                 IF(iomega.GT.1) tmp = tmp + 0.5d0*domega*(dos(iomega) + dos(iomega-1))
                 WRITE(29,'(20e20.10)') domega*(iomega-nomega_fine/2.d0), dos(iomega), dosqp(iomega), dosGW(iomega), tmp, DBLE(EkSCF(ik2))
              END DO
              WRITE(29,*)
              dos = 0.d0
              dosGW = 0.d0
              dosOcc = 0.d0
           END IF
           xklast = xk
        END DO
     END IF
  END DO
  IF(master) THEN
     !dos = dos*3.d0/kFermi**3
     norm = 0.d0
     DO iomega = 1, nomega_fine
        tmp = (iomega-nomega+1)*domega
        IF(iomega.EQ.1) THEN
           norm = norm + dos(iomega)*0.5d0
        ELSE
           norm = norm + 0.5d0*(dos(iomega-1)+dos(iomega))*domega
        END IF
        !WRITE(29,*) tmp, dos(iomega), norm
     END DO
     !PRINT*, 'Total Energy: ', ETotRPA*3.d0/kFermi**3/2.d0
     !PRINT*, 'Total Energy (AH): ', ETotAH*3.d0/kFermi**3/2.d0
     !PRINT*, 'E0: ', E0/(kFermi**3/3.d0)/2.d0
     !PRINT*, 'E0Test: ', E0Test/(kFermi**3/3.d0)/2.d0
     !PRINT*, 'E0 Exact: ', 3.d0*kFermi**2/10.d0
     !PRINT*, 'Kinetic Energy (GC): ', EKin*3.d0/kFermi**3/2.d0
     !PRINT*, 'Kinetic Energy (GW): ', EKinGW*3.d0/kFermi**3/2.d0
     !PRINT*, 'Exchange Energy (Calc.): ', Ex*3.d0/kFermi**3/2.d0
     !PRINT*, 'Exchange Energy (Form.): ', Eexch
     !PRINT*, 'Total Occupation/Expected Occupation (GC, GW): ', ntot/(kFermi**3/3.d0), ntotGW/(kFermi**3/3.d0)
     !PRINT*, '1/2 Mu: ', 0.5d0*Mu
     !WRITE(27,'(A)'), 'EGC, EGW, ESB, Ex, E0, totocc/(kFermi**3/3.d0)'
     !WRITE(27,'(20e20.10)'), ETotRPA/(kFermi**3/3.d0)/2.d0, ETotGW/(kFermi**3/3.d0)/2.d0, ETotAH/(kFermi**3/3.d0)/2.d0, Ex/(kFermi**3/3.d0)/2.d0, E0/(kFermi**3/3.d0)/2.d0, ntot/(kFermi**3/3.d0)
  END IF
  CALL par_barrier
  IF(master) PRINT*, 'Finished with spectral function calculation.'
  CALL par_end

            CONTAINS
              COMPLEX*16 FUNCTION CIntegrate(CIntegrand,n,step)
                COMPLEX*16 CIntegrand(n)
                REAL(8) step
                INTEGER n, i
                IF(n.LT.3) THEN
                   PRINT*, 'Error in CIntegrate: n must be 3 or larger.'
                   STOP
                END IF
                CIntegrate = CIntegrand(1) + 2.d0*CIntegrand(2) + 2.d0*CIntegrand(n-1) + CIntegrand(n)
                DO i = 3, n-2
                   CIntegrate = CIntegrate + 4.d0*CIntegrand(i)
                END DO
                CIntegrate = 0.25d0*step*CIntegrate
                RETURN
              END FUNCTION CIntegrate

              ! COMPLEX*16 FUNCTION CIntegrate2(beta,omega,n,t)
              !   REAL(8), INTENT(IN) :: beta(n), omega(n), t
              !   REAL(8) a1, b1
              !   INTEGER i
              !   IF(t.EQ.0) THEN
              !      CIntegrate2 = 0.d0
              !      RETURN
              !   END IF
              !   DO i = 1, n-1


              COMPLEX*16 FUNCTION SolveZSC(a0)
              ! Solve equation Z = Exp(-a0*Z) for Z
                COMPLEX*16, INTENT(IN) :: a0
                COMPLEX*16 Z(2)
                REAL(8) tol, a0tmp, alpha
                INTEGER nIter
                REAL(8), PARAMETER :: tol2 = 1.d-5
                INTEGER, PARAMETER :: MaxIter = 1000
                tol = 1.d-6
                alpha = 0.1d0
                ! Start with Z = Exp(-a0)
                Z(1) = Exp(-a0)
10              CONTINUE
                DO nIter = 1, MaxIter
                  Z(2) = Exp(-a0*DBLE(Z(1)))
                  IF(ABS(Z(2)-Z(1)).LT.tol) EXIT
                  Z(1) = Z(1)*(1-alpha) + Z(2)*alpha
                END DO
                IF(nIter.EQ.MaxIter) PRINT*, "MaxIter reached in SolveZSC"

                ! Check that solution has been found.
                !IF(ABS(Z(2) - Exp(-Z(2)*a0)).LT.tol2) THEN
                ! If Re[Z] > 1, take Re[Z] = 1
                IF(DBLE(Z(2)).GT.1.d0) THEN
                  SolveZSC = 1.d0
                ELSE
                  SolveZSC = Z(2)
                END IF
                !ELSE
                !  tol = tol/10.d0
                !  GOTO 10
                !END IF
                RETURN
              END FUNCTION SolveZSC

END PROGRAM spfcn2
