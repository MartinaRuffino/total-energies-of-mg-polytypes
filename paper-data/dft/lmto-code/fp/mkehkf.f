      subroutine mkehkf(mode,s_ham,sev,valmom,sumtv,etot)
C- Make Harris energy or Kohn-Sham energy
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  eterms
Co     Stored:     eterms
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1 Harris Foulkes energy
Ci         :2 Kohn-Sham energy
Ci   sev   :band sum
Ci   valmom:magnetic moment (not used unless mode=1)
Cio Inputs/Outputs
Cio  sumtv :valence kinetic energy :
Cio        :output for mode=1
Cio        :input for mode=2
Co   sham->eterms
Ci         :... the following quantities are input
Ci         :(3)  utot   = total electrostatic energy
Ci         :(4)  valves = valence rho * estat potential (not used)
Ci         :(5)  cpnves = core+nuc * estat potential (not used)
Ci         :(6)  rhoexc = rho * exc
Ci         :(7)  rhovxc = rho * vxc (not used)
Ci         :(8)  sumec  = sum of foca=0 core eigenvalues
Ci         :(9)  sumtc  = sum of core kinetic energies
Ci         :(10) rhcvef1= int (rhoc*vef1) for sites w/ lfoca=0
Ci         :(11) valvef = smrhov * vsm + sum_ib valvef_ib
Ci                        valvef_ib = rhov * vtrue - smrhov * vsm)_ib
Ci         :(12) sumt0  = sum of frozen core kinetic energies
Co         :(13) dq1      (not used here)
Co         :(14) dq2      (not used here)
Ci         :(15) amom   = system magnetic moment (mode=2 only)
Ci         :              If mode=1, input moment is supplied by valmom
Co         :(16) sumev    (not used here) sphere sum-of-eigenvalues
Co         :(17) rinvxt   (not used here)
Co         :(18) rouvxt   (not used here)
Ci         :(19) rhosig = d.c. from self-energy sigma-vxc (eks only)
Ci                        Not used if rhosig=NULLI
Ci         :(20) Bext.M = d.c. from external field.
Ci         :(21) Beff.M = d.c. from effective field
Ci         :              (Fixed spin moment method)
Co         :... the following quantities are stored:
Co         :(1 )  eh    = etot = Harris energy (mode=1)
Co         :(2 )  eks   = etot = Hohenberg-Kohn-Sham energy (mode=2)
Co         :(13)  sumev = sum of eigenvalues
Co         :(15)  amom  = magnetic moment from bands (mode=1)
Co Outputs
Co   etot  :Harris energy (mode = 1)
Co         :Hohenberg-Kohn-Sham energy (mode = 2)
Cr Remarks
Cr   Information related to total energy is printed out, depending on mode
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   27 Jan 11 Change sign of B in double-counting, Bext.M and Beff.M
Cu   19 Jul 10 Correct double counting in fixed-spin-moment method
Cu   01 Apr 10 Double counting term from external B field for ehar
Cu   02 Jan 06 stores or prints magnetic moments
Cu   11 Jan 05 double-counting term rho*sig subtracted from ehks.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode
      double precision sev,etot,sumtv,valmom
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Local parameters
      integer stdo,stdl,ipr,ipl,nsp,nglob
      double precision eh,eks,rhoexc,sumec,utot,valvef,rhvext,rhcvef1,
     .  ttcor,sumt0,sumev,sumtc,ekin,rhosig,amom,bfam,bfsmom
C     double precision cpnves,rhovxc,sumtc,valves
C     This ordering must match sham->eterms; see uham
      integer, parameter :: NULLI=-99999
      double precision eterms(22)
      equivalence (eterms(1),eh)
      equivalence (eterms(2),eks)
      equivalence (eterms(3),utot)
C     equivalence (eterms(4),valves)
C     equivalence (eterms(5),cpnves)
      equivalence (eterms(6),rhoexc)
C     equivalence (eterms(7),rhovxc)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(10),rhcvef1)
      equivalence (eterms(11),valvef)
      equivalence (eterms(12),sumt0)
      equivalence (eterms(13),sumev)
      equivalence (eterms(15),amom)
      equivalence (eterms(17),rhvext)
      equivalence (eterms(19),rhosig)
      equivalence (eterms(20),bfam)
      equivalence (eterms(21),bfsmom)

C     stdo = lgunit(1)
C     stdl = lgunit(2)
      nsp  = nglob('nsp')
      stdo = nglob('stdo')
      stdl = nglob('stdl')
      call getpr(ipr)
      ipl = 2; if (stdl < 0) ipl = 0
      eterms = s_ham%eterms
      sumev = sev

      if (mode == 1) then
        sumtv = sev - valvef - rhvext
        ttcor = sumec - rhcvef1 + sumt0
        ekin  = sumtv + ttcor
        etot  = ekin + utot + rhoexc
C   ... double counting from bext * local moment
        if (bfam /= NULLI) then
          etot = etot - bfam
        endif
        if (bfsmom /= NULLI) then
          etot = etot - bfsmom
        endif
        eh    = etot
        amom  = valmom
      elseif (mode == 2) then
C   ... double counting from sigma
        if (rhosig /= NULLI) then
          sumtv = sumtv - rhosig
        endif
        if (bfsmom /= NULLI) then
          sumtv = sumtv - bfsmom
        endif
        ekin = sumtv + sumtc
        etot = ekin + utot + rhoexc
        eks  = etot
      endif

      s_ham%eterms = eterms
      if (mode == 0) return

C      print *, sumev(1,1)-sev
C      print *, 'xxsev=',sev
C      print *, 'xxvalvef*-1=',-valvef
C      print *, 'xxttcor=',ttcor
C      print *, 'xxutot=',utot
C      print *, 'xxrhoexc=',rhoexc

C --- Printout ---
      if (mode == 1) then
C        if (rhvext /= 0) then
C        call info2(30,1,-1,' Contribution to Harris energy from external potential:'//
C     .    '%N val*vext=%;13,6D  cor*vext=%;14,6D',rhvext,rhcvef1)
C      endif
      if (ipr >= 30) then
        write(stdo,660) sev,valvef,sumtv,sumec,rhcvef1,ttcor,rhoexc,utot,eh
  660 format(/' Harris energy:'
     .       /' sumev= ',f15.6,'  val*vef=',f15.6,'   sumtv=',f15.6
     .       /' sumec= ',f15.6,'  cor*vef=',f15.6,'   ttcor=',f15.6
     .       /' rhoeps=',f15.6,'     utot=',f15.6,'    ehar=',f15.6)
      elseif (ipr >= 10) then
        call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d  sumev=%,6;6d'//
     .    '  ehf=%,6;6d',' ',80,stdo,ekin,utot+rhoexc,sev,eh)
      endif
      if (ipr > 1 .and. ipl > 1) write (stdl,720) utot,ekin,sev,etot
  720 format('fp Har U',f15.6,'  T',f15.6,'  sev',f14.6,'  EH ',f14.6)

      elseif (mode == 2) then
      if (ipr >= 30) then
        if (rhosig /= NULLI .and. rhosig /= 0) then
          write (stdo,411) sumtv,sumtc,ekin,rhoexc,utot,rhosig,eks
  411     format(/' Kohn-Sham energy:'
     .           /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6
     .           /' rhoep= ',f15.6,'   utot=',f17.6,' rhosig=',f16.6
     .           /'  ehks=', f16.6)
        else
          write (stdo,410) sumtv,sumtc,ekin,rhoexc,utot,eks
  410     format(/' Kohn-Sham energy:'
     .           /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6
     .           /' rhoep= ',f15.6,'   utot=',f17.6,'   ehks=',f16.6)
        endif
        if (nsp == 2) write (stdo,412) valmom,amom
  412   format(' mag. mom=',f13.6,'  (bands)',f16.6,'  (output rho)')

      elseif (ipr >= 10) then
        call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d  ehf=%,6;6d'//
     .    '  ehks=%,6;6d',' ',80,stdo,ekin,utot+rhoexc,eh,eks)
      endif
      if (ipr > 1 .and. ipl > 1) write (stdl,721) utot,ekin,rhoexc,eks
  721 format('fp KS U',f16.7,' T',f16.7,' Exc',f15.7,' EKS',f15.7)

      endif

C      ifet=1456
C      open(ifet,file='ETOTLDA')
C      write(ifet,"(d23.16,a)")  eks,    ' ! EKS  (Ry)'
C      write(ifet,"(d23.16,a)")  rhoexc, ' ! \int rho exc '
C      write(ifet,"(d23.16,a)")  utot,   ' ! U '
C      write(ifet,"(d23.16,a)")  ekin,   ' ! T '
C      close(ifet)

      end
