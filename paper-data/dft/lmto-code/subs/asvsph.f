      subroutine asvsph(s_ctrl,s_lat,s_spec,s_ham,s_pot,vrl,mode,ehterm,lhave)
C- Make ASA sphere potential, potential parameters
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nbasp nclass nspec zbak nl nspin lncol loptc
Ci                 nclasp ldlm nccomp lpgf idcc nrc ncomp ics initc
Ci                 lgen3 lsx smalit lmet npadl npadr ipcp ips ipc rmax
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed: lasa lves ipc lbxc dclabl rmax initc nrc lncol lrel
Cio                lcd lscr lham ics ncomp ipcp idcc ips
Cio    Passed to:  atscpp asamad asavqm asaqmp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  vol avw alat plat as tol nkdmx nkqmx platl platr
Ci                 qlat awald nkd nkq nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:vol pos cg jcg indxcg qlv dlv symgr ag
Cio    Passed to:  asamad asavqm asaqmp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  eref idmod lmxa z a nr coreh coreq name lmxl lmxf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  atscpp gtpcor asamad getq dlmq asavqm
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  nmto kmto neula eterms bdots etrms eula hrs iaxs
Ci                 iprmb lmxa magf nprs offH qsig
Co     Stored:     eterms ehk thrpv seref amgm bdots etrms eula hrs
Co                 iaxs iprmb lmxa magf nprs offH qsig
Co     Allocated:  *
Cio    Elts passed:eula etrms
Cio    Passed to:  bcast_strx iinit mpibc1
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  gibbs bxc qt ves vconst qpp vmtz0 aamom cp ddpf
Ci                 dddpf ddpfr dlmwt dpf dpfr gma gmar grrme mad mxy
Ci                 palp papg pf pfnc pfr pmpol pnu pp ppn pprel pti qc
Ci                 qcorr qnu rhat rnew rhos rhrmx sop thetcl vdif vintr
Ci                 vrmax vshft smpot smrho smrout
Co     Stored:     vmtz vconst aamom bxc cp ddpf dddpf ddpfr dlmwt dpf
Co                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Co                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Co                 rhat rnew rhos rhrmx sop thetcl vdif ves vintr vrmax
Co                 vshft smpot smrho smrout
Co     Allocated:  *
Cio    Elts passed: thetcl pnu qnu ves bxc pp pprel qt rhrmx vrmax sop
Cio                grrme pmpol vintr vdif gibbs mad qc qpp
Cio    Passed to:  asamad asavqm asaqmp bcast_strx iinit mpibc1
Ci Inputs
Ci   vrl   :(pgf) difference in potential between left- and right- leads.
Ci         :      See lmasa-gf.f
Ci   mode  : controls program flow of atom program
Ci         :1s digit = imake : tells sphere program what to generate
Ci         :   Note: nonzero 100s digit may modify this digit.
Ci         :   Also see 100s digit for definition of "given potential"
Ci         :0  Read double counting terms from atom files.
Ci         :   No potential or pot pars calculated; no atom file written
Ci         :   This is a 'cheap' version of imake=4, where
Ci         :   instead of computing the d.c. terms from the given potential,
Ci         :   they are copied from atom files.
Ci         :1  Make self-consistent potential from given pnu,qnu
Ci         :   No potential parameters calculated.
Ci         :2  Like imake=0, but make ppars from given potential
Ci         :   Passed pnu,qnu written to disk.
Ci         :   NB: potential, ppars need not be related to pnu,qnu
Ci         :3  Make self-consistent potential from given pnu,qnu
Ci         :   and potential parameters from resulting potential
Ci         :4  Make sphere double-counting terms from given potential
Ci         :   and supplied moments P,Q.  No internal self-consistency
Ci         :   in the potential, nor are potential parameters generated;
Ci         :   ves = estat potential for supplied Q is not retained.
Ci         :   no atom file written.  Using this mode to make terms
Ci         :   in Kohn-Sham energy, e.g. v_in, rho_out.
Ci         :10s digit: deals with source of potential
Ci         :   In each case potential may or may not be available.
Ci         :0  potential read from disk
Ci         :1  potential copied from s_pot, never read from disk.
Ci         :2  potential copied from s_pot if available, otherwise from disk.
Ci         :   In modes 1 and 2, 1s digit mode
Ci         :100s digit: deals with generation of ppars from potential
Ci         :0  Standard operation: potential may generated internally by
Ci         :   sphere program (see 1s digit), and becomes the "given potential"
Ci         :1  Potential not generated self-consistently by sphere program:
Ci         :   potential parameters are made if given potential is available
Ci         :   Note: 1s modes 1 and 3 are nonsensical in this case.
Ci         :   If 1s mode= 1 -> aborts; if 1s digit 3, switched to mode 2
Ci         :1000s digit: deals with renormalization of qnu from potential
Ci         :   Potential at rmt
Co Outputs
Co   ehterm:(1) ehkk (no mad)
Co          (2) sphere band sum (ves=0 at rmt)
Co          (3) sum q_i v(r)
Co          (4) emad
Co          Note: ehterm is obsolete and will be phased out in favor of
Co          sham->eterms
Co   sham->eterms integrals for the total energy are accumulated
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves --- not used by ASA
Co         :(5)  cpnves --- not used by ASA
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc
Co         :(8)  sumec  = sum-of-core eigenvalues
Co         :(9)  sumtc  --- not used by ASA
Co         :(10) xcore  = rhoc * total potential
Co         :(11) valvef = rhov * total potential
Co         :(12) sumt0  --- not touched here
Co         :(13) dq1    --- not touched here
Co         :(14) dq2    --- not touched here
Co         :(15) amom   =  total magnetic moment
Co         :(16) sumev  =  sum-of-eigenvalues from moments
Co         :(17) rinvxt --- not touched here
Co         :(18) rouvxt --- not touched here
Co         :(19) bmval  = M<B> : magnetic contribution to valvef (asadc)
Co   lhave :1s digit
Co          0 insufficient input available to generate potentials
Co            sought
Co          1 All potentials sought were calculated or read from disk
Co         10s digit
Co          0 insufficient input available to generate ehkk,ehterm
Co          1 ehkk, ehterm have been generated
Co   See Remarks
Cl Local variables
Cl    imake :1s digit mode
Cl    ehkk  :ASA total energy as sum-of-sphere energies
Cl          :This is the HK total energy in versions 6.12 and earlier
Cb Bugs
Cb   s_pot%grrme should be merged with s_optic%rgrad
Cr Remarks
Cu Updates
Cu   28 Oct 17 Extended functionality of mode (formerly imake)
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   28 May 08 (Kirill) extensions for disordered local moments
Cu   21 Jul 07 (pgf) vne->vrl (for inequivalent left- and right- end layeers)
Cu   10 Feb 04 (S.Faleev) vrl added to argument list; passed
Cu              to asamad for non-equilibrium mode
Cu    1 Apr 04 Sphere program for orbital-dependent XC field
Cu   26 Apr 03 Added MPI calls
Cu   18 Mar 03 (A Chantis) relativistic potential parameters.
Cu   10 Mar 03 First cut at making proper Kohn-Sham energy.
Cu   17 Feb 03 Added double-counting terms to ehk for applied field
Cu   15 Feb 03 B-field requires SO parameters
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lhave,mode
      double precision ehterm(5),vrl
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      logical lehk,havedc,havedci
      integer ic,nclass,nclasp,nclspp,nbas,nl,nsp,i,nbasp,lpgf,imake,is,
     .  lncol,loptc,lnsph,isw,lgunit,lves,k,nrclas,nspec,neul
      integer mpipid,procid
      double precision thrpv,amgm,ehkk,emad,trumad,vmtz(2),zbak(2),
     .  wk(5),amag(3),dsqrt,bscal,facso,xx
      real(8), pointer :: p_ves(:)
C ... For atscpp
      character clabl*8,outs*256
      double precision amgmat,amomnc(3),ehkat,sevat,thrpva,vol,avw,seref
C ... For third-generation LMTO
      integer nmto
      double precision kmto(20)
C ... Dynamically allocated arrays
      real(8), allocatable :: eulat(:),wkl(:)
C ... For DLM
      integer ldlm,nangl,icp
      integer nclspd,nclsppd
      double precision wt,costh,ehbak
C     This ordering must match sham->eterms; see uham
      double precision eterms(22),eterma(22)

      procedure(integer) wordsw
      procedure(logical) cmdopt

      real(8), target :: tmp(2)
      procedure(real(8)) :: ddot
      procedure(integer) :: iclbsj
      procedure (logical) :: bittst

      imake = mod(mode,10)
      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      nclass = s_ctrl%nclass
      nspec = s_ctrl%nspec
      zbak = s_ctrl%zbak
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lncol = s_ctrl%lncol
      loptc = s_ctrl%loptc
      nclasp = s_ctrl%nclasp
      nmto = s_ham%nmto
      kmto(1:6) = s_ham%kmto
      neul = s_ham%neula
      vol = s_lat%vol
      avw = s_lat%avw
      ldlm = s_ctrl%ldlm
      nangl = s_ctrl%nccomp
      lpgf = s_ctrl%lpgf(1)
      lnsph = isw(IAND(s_ctrl%lasa,32) /= 0)
      lves = 2*IAND(s_ctrl%lves,1)
      if (lpgf /= 0) then
        lves = lves+100
      endif
      ehkk = 0
      amgm = 0
      thrpv = 0
      seref = 0
      nclspp = 2*nclasp-nclass
      nclspd = nclasp + nangl
      nclsppd = nclspp + nangl
      call dpzero(ehterm,4)
      call dpzero(amag,3)
      eterms = s_ham%eterms
C     asvsph only modifies eterms(3:11)
      call dvset(eterms,3,11,0d0)
      call dvset(eterms,15,16,0d0)

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid == 0) then

      call togprt()
      allocate(wkl(nclspp)); call dpzero(wkl,nclspp)
      lehk = .true.
      lhave = 1
      havedc = .true.

      if (.not. associated(s_pot%grrme)) s_pot%grrme => tmp
      if (.not. associated(s_pot%pmpol)) s_pot%pmpol => tmp
      if (.not. associated(s_pot%vintr)) s_pot%vintr => tmp

      if (cmdopt('--zerq',6,0,outs)) then
        call asazerq(outs(7:),1,s_ctrl,s_lat,s_spec,s_pot)
      endif

      do  ic = 1, nclspd
        havedci = .true.
        icp = ic
        if (ic > nclasp) icp = s_ctrl%idcc(ic)
        nrclas = s_ctrl%nrc(icp)
C   ... Classes to skip
        i = iclbsj(icp,s_ctrl%ipc,-nbasp,1)
        if (i < 1) cycle
C   ... Skip parent DLM classes
        if (ldlm > 0 .and. ic <= nclasp) then
          if (s_ctrl%ncomp(ic) >= 2) cycle
        endif
C   ... Retrieve weights for DLM classes
        wt = 1d0
        costh = 1d0
        if (ic > nclasp) then
          wt = s_pot%gibbs(ic)
          costh = dcos(s_pot%thetcl(ic-nclasp,1))
        endif
C   ... What to make
        if (lpgf == 2 .and. i <= nbas) cycle
        is = s_ctrl%ics(ic)
C       call pshpr(80)
        call dpzero(eterma,22)
        if (neul > 0) then
          allocate(eulat(neul*3))
          call dmscop(eulat,1,s_ham%eula,nbas,i,i,1,neul*3,1,1,1d0)
        else
          allocate(eulat(1))
        endif
        if (mod(s_ctrl%lbxc,2) == 1) then
          bscal = 1d0+s_pot%shfac(1,ic)
        elseif (mod(s_ctrl%lbxc,4) == 2) then
          bscal = dsqrt(1d0+s_pot%shfac(1,ic)**2)
        else
          bscal = 1d0
        endif
C        if (bscal /= 1d0) write (6,901) ic, bscal
C  901   format(' Bxc for ic = ',i3,' scaled by ',f9.5)
C        facso = 1 + s_pot%shfac(2,ic)
        facso = 1  ! May 2015 full SOC written to disk. scale them when they are used.

        call atscpp(s_ctrl,s_spec,s_pot,is,ic,s_ctrl%dclabl,s_ctrl%rmax,
     .    mode,nl,nsp,s_ctrl%initc(ic),zbak(1)/vol,avw,s_pot%pnu,
     .    s_pot%qnu,s_pot%ves,wkl,eulat,bscal,facso,neul,s_pot%bxc,0d0,
     .    s_pot%pp,s_pot%pprel,ehkat,sevat,s_pot%qt,amgmat,amomnc,
     .    s_pot%rhrmx,s_pot%vrmax,thrpva,s_pot%sop,s_pot%grrme,
     .    s_pot%pmpol,s_pot%vintr,clabl,eterma)
        if (neul > 0) then
          call dpsadd(amag,amomnc,3,1,1,dble(nrclas)*wt)
        endif
        deallocate(eulat)

C       If estat energy changed, double-counting terms are available
        havedci = eterma(3) /= 0 ! utot is nonzero
        if (mode == 3 .and. mod(s_ctrl%initc(ic),2) == 1 .and. s_spec(is)%z == 0) havedci = .true. ! ES may have no charge
        havedc = havedc .and. havedci
        if (havedci) then
          eterma(4) = s_pot%ves(ic) * s_pot%qt(ic) / 2
          call dpadd(eterms,eterma,3,11,dble(nrclas)*wt)
          eterma(15) = eterma(15) * costh
          call dpadd(eterms,eterma,15,16,dble(nrclas)*wt)
          call dpadd(eterms,eterma,19,19,dble(nrclas)*wt)
          if (ldlm > 0 .and. associated(s_ham%etrms)) call pketerms(imake,eterma,s_ham%etrms,ic)
        elseif (s_spec(is)%z /= 0) then
          call dvset(eterms,3,11,0d0)
        endif
C       call poppr
        seref = seref + s_spec(is)%eref*nrclas*wt
        if (imake == 4) cycle

C   ... Check that potential and other pars available if to continue
        k = s_ctrl%initc(ic)
        if (s_ctrl%lgen3 == 0) then
        if (mod(k/2,2) /= 1
     .    .or.(bittst(lncol,4).or.bittst(lncol,8)).and.(mod(k/4,2)) /= 1
     .    .or. loptc > 0 .and. (mod(k/32,2)) /= 1
     .    .or. lnsph /= 0 .and. (mod(k/16,2)) /= 1) then
          call awrit0(' LM (warning) class '//clabl//
     .     '%a missing pot, s-orbit or optical parms',' ',-80,lgunit(1))
          lhave = 0
          cycle
        endif
        endif

C   ... Only accumulate into ehkk if ehkat is set
        if (havedc .and. lehk) then
          ehkat = ehkat - s_spec(is)%eref
          ehkk = ehkk + ehkat*nrclas*wt
        else
          lehk = .false.
        endif
C  ...  Correct ehterm(2) if wt is not 1 (needed for DLM)
        if (wt /= 1d0) ehbak = ehterm(2)
        call atsev(ic,s_ctrl%nrc,s_pot%ves,nl,nsp,s_pot%qnu,sevat*wt,ehterm(2))
        if (wt /= 1d0) ehterm(2) = ehbak + wt*(ehterm(2)-ehbak)
        amgm = amgm + amgmat*nrclas*wt
        thrpv = thrpv + thrpva*nrclas*wt
      enddo

      if (associated(s_pot%grrme, tmp)) s_pot%grrme => null()
      if (associated(s_pot%pmpol, tmp)) s_pot%pmpol => null()
      if (associated(s_pot%vintr, tmp)) s_pot%vintr => null()

      deallocate(wkl)
      call togprt()

C --- Make ves, emad, vmtz for new sphere Q, potential ---
C     imake=4: don't update pot->ves; use temporary array
      if (imake == 4) then
        allocate(p_ves(nclsppd))
        call dpzero(p_ves,nclsppd)
      else
        p_ves => s_pot%ves
      endif
      call asamad(s_ctrl,s_pot,s_lat,s_spec,lves,s_pot%pnu,
     .  s_pot%qnu,vrl,p_ves,emad,trumad,vmtz,s_ham%etrms)
C ... V[q] needed for total energy if a potential shift
      call dpzero(s_pot%vdif,nclsppd)
      if (mod(lves,100) == 2) then
        call pshpr(1)
        allocate(wkl(22*nclsppd))
        call asamad(s_ctrl,s_pot,s_lat,s_spec,0,s_pot%pnu,
     .    s_pot%qnu,vrl,s_pot%vdif,wk(1),wk(2),wk(3),wkl)
        deallocate(wkl)
        call daxpy(nclsppd,-1d0,p_ves,1,s_pot%vdif,1)
        call poppr
      endif
      if (imake == 4) deallocate(p_ves)

C ... Add Madelung contribution to total electrostatic energy
      eterms(3) = eterms(3) + emad
C ... Total magnetization, noncollinear case
      if (neul > 0 .and. ldlm == 0) then
        amgm = dsqrt(ddot(3,amag,1,amag,1))
        eterms(15) = amgm
      endif
C ... Repack double-counting terms
      s_ham%eterms = eterms

C --- Make NMTO potential parameters ---
CC     if (lgors('ctrl lgen3',sctrl)) then
C      if (s_ctrl%lgen3 /= 0) then
C        call ppmesh(nl,nmto,kmto,s_ctrl%dclabl,
C     .    vmtz,s_ctrl%ics,w(oppn))
C      endif

C --- Save energy for HF (to recover if information lost) ---
C terms:
C (1) ehkk (no Mad) (2) sumeV(VH=0)
C (3) sum Q_i V(R) (4) emad
C (5) d.c. terms from applied fields.
      ehterm(4) = emad
      if (lehk) then
        ehterm(1) = ehkk
        ehkk = ehkk + emad - ehterm(5)
      endif

      thrpv = thrpv + trumad
      s_pot%vmtz = vmtz(1)
      s_ham%ehk = ehkk
      s_ham%thrpv = thrpv
      s_ham%seref = seref
      s_ham%amgm = amgm
      if (lehk) lhave = lhave+10

      endif

C --- MPI broadcast everything passed out of asvsph ---
      call mpibc1(s_ctrl%initc,nclspd,2,.false.,'asvsph','initc')
      call bcast_strx(2**4,xx,xx,s_ham,xx,xx,xx,xx,xx,xx,0,0)
      call mpibc1(s_pot%rhrmx,nclspd,4,.false.,'asvsph','rhrmx')
      call mpibc1(s_pot%ves,nclsppd,4,.false.,'asvsph','ves')
      call mpibc1(s_pot%vdif,nclsppd,4,.false.,'asvsph','vdif')
      call mpibc1(s_pot%vrmax,2*nclsppd,4,.false.,'asvsph','vrmax')
      if (associated(s_pot%vintr)) then
        i = nclspd*(nl*nsp)**2
        call mpibc1(s_pot%vintr,i,4,.false.,'asvsph','vintr')
      endif
      call mpibc1(s_pot%pnu,nl*nsp*nclspd,4,.false.,'asvsph','pnu')
      call mpibc1(s_pot%qnu,3*nl*nsp*nclspd,4,.false.,'asvsph','qnu')
      call mpibc1(s_pot%qt,nclspd,4,.false.,'asvsph','qt')
      call mpibc1(s_pot%pp,6*nl*nsp*nclspd,4,.false.,'asvsph','pp')
C     if (lrel == 2) then
      if (associated(s_pot%pprel)) then
        call mpibc1(s_pot%pprel,5*8*nl*nl*max(nclsppd,nspec),4,.false.,'asvsph','pprel')
      endif

      if (s_ctrl%lncol > 0) then
        i = size(s_pot%sop)
        call mpibc1(s_pot%sop,i,4,.false.,'asvsph','sop')
      endif

      if (s_ctrl%loptc > 0) then
         i = size(s_pot%grrme)
        call mpibc1(s_pot%grrme,i,4,.false.,'asvsph','grrme')
      endif
      if (associated(s_pot%pmpol)) then
       i = (2*nl-1)*nl**2*3*nsp*nclspd
       call mpibc1(s_pot%pmpol,i,4,.false.,'asvsph','pmpol')
      endif
C     call mpibc1(spot,nint(spot(1)),4,.false.,'asvsph','spot')
      call bcast_strx(2**5,xx,xx,xx,s_pot,xx,xx,xx,xx,xx,0,0)
C     Remaining arguments passed to asvsph
      call mpibc1(ehterm,5,4,.false.,'asvsph','ehterm')
      call mpibc1(lhave,1,2,.false.,'asvsph','lhave')

      end
      subroutine pketerms(imake,eterms,etrms,ic)
C- Package eterms into etrms (needed for DLM)
      double precision eterms(22),etrms(22,ic)
      integer ic,imake
      integer i
      do  i = 3, 11
        etrms(i,ic) = eterms(i)
      enddo
      if (imake /= 4) then
        etrms(16,ic) = eterms(16)
      endif
      end

      subroutine asazerq(opts,mode,s_ctrl,s_lat,s_spec,s_pot)
C- Adds or scales charges in a list of spheres to make system neutral
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   opts  : string containing switches
Ci   mode  : 0 => always act
Ci         : 1 => act only if opts contains 'qin'
Ci         : 2 => act only if opts contains 'qout'
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   15 Nov 17
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character(len=*) :: opts
      integer mode
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
!      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      integer, allocatable :: lmx(:),iclst(:),iwkv(:)
      real(8), allocatable :: dq(:),z(:)
C ... Local parameters
      character dc*1
      integer i,j,n,ic,icp,ipass,j1,j2,nclasp,nclspd,nrclas,nangl,nl,nsp,nlstc,iv(10),iv0,ib
      double precision vol,vsph,qtot,xv(10),fracq,qsph,qnow,qlat(3,3),h(3),dpos(3),avfracq,qin,qout
      real(8), parameter :: qmin = .0001d0
      procedure(logical) :: a2bin
      procedure(integer) wordsw,mkilsd
      procedure(real(8)) :: dsum,ddot

      if (s_ctrl%lpgf(1) == 2) return ! Not ready for pgf end layers
C      if (opts == ' ') then
C        call info0(20,1,0,' ASAZERQ: empty options ... nothing done')
C        return
C      endif
      dc = opts(1:1)

      if (wordsw(opts,dc,'qin','',j1) /= 0) then
        if (mod(mode,10) /= 1) return
      elseif (wordsw(opts,dc,'qout','',j1) /= 0) then
        if (mod(mode,10) /= 2) return
      endif

C     Get the total charge
      nclasp = s_ctrl%nclass  ! For pgf=1, scale charges in active region only
      nangl = s_ctrl%nccomp
      nclspd = nclasp + nangl
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      allocate(z(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'z',1,xv,z)
      allocate(lmx(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'lmxa',1,lmx,xv)
      allocate(dq(nclspd))
      call getq(nsp,nl,lmx,nclspd,z,s_pot%pnu,s_pot%qnu,
     .  s_ctrl%ics,s_spec,s_pot%qc,s_pot%qt,dq)

      qtot = 0
      do  ic = 1, nclspd
        icp = ic
        if (ic > nclasp) icp = s_ctrl%idcc(ic)
        nrclas = s_ctrl%nrc(icp)
        qtot = qtot + nrclas*s_pot%qt(ic)
      enddo
      call info5(25,1,0,' asazerq : system charge is now %d',qtot,2,3,4,5)

C --- Get list of classes
      nlstc = 0
      if (wordsw(opts,dc,'es','',j1) /= 0) then
        do  ipass = 1, 2
          nlstc = 0
          do  ic = 1, nclspd
            icp = ic
            if (ic > nclasp) icp = s_ctrl%idcc(ic)
            if (z(icp) > 0.5d0) cycle
            nlstc = nlstc+1
            if (ipass == 1) cycle
            iclst(nlstc) = ic
          enddo
          if (ipass == 1) allocate(iclst(nlstc))
        enddo
      elseif (wordsw(opts,dc,'expr:','',j1) /= 0) then
        call nwordg(opts,0,dc//' ',1,j1,j2)

C       This should be bundled into a separate routine. See stackel
        call numsyv(iv0)
        allocate(iwkv(s_ctrl%nbasp)); call ivset(iwkv,1,s_ctrl%nbasp,s_ctrl%nbasp+1)
        do  ib = 1, s_ctrl%nbasp
C         call shosyv(0,0,0,6)

C         projection of bas along qlat(i)
          call dinv33(s_lat%plat,1,qlat,vol)
          do  i = 1, 3
            h(i) = ddot(3,s_lat%pos(1,ib),1,qlat(1,i),1)
            dpos(i) = ddot(3,s_lat%pos(1,ib),1,s_lat%plat(1,i),1)
          enddo

          call lodsyv('ib',1,dble(ib),n)
          i = s_ctrl%ips(ib)
          call lodsyv('is',1,dble(i),n)
          call lodsyv('z',1,s_spec(i)%z,n)
          i = s_ctrl%ipc(ib)
          call lodsyv('ic',1,dble(i),n)
          call lodsyv('x1',1,s_lat%pos(1,ib),n)
          call lodsyv('x2',1,s_lat%pos(2,ib),n)
          call lodsyv('x3',1,s_lat%pos(3,ib),n)
          call lodsyv('h1',1,h(1),n)
          call lodsyv('h2',1,h(2),n)
          call lodsyv('h3',1,h(3),n)
C         call lodsyv('h',1,h(3),n)
          call lodsyv('p1',1,dpos(1),n)
          call lodsyv('p2',1,dpos(2),n)
          call lodsyv('p3',1,dpos(3),n)

C         call shosyv(0,0,0,6)
          j = 0
          if (.not. a2bin(opts(j1:j2)//' ',xv,4,0,' ',j,-1))
     .      call rx('asazerq failed to parse expr: '//opts(j1:j2))
          if (xv(1) /= 0) then
            ic = s_ctrl%ipc(ib)
            iwkv(ib) = ic
          endif
          call clrsyv(iv0)
        enddo
        allocate(iclst(s_ctrl%nbasp))
        call ivheap(1,s_ctrl%nbasp,iwkv,iclst,0d0,100)
        nlstc = 0
        do  ib = 1, s_ctrl%nbasp
          if (iwkv(ib) > s_ctrl%nbasp) cycle
          if (nlstc > 0) then
            if (iwkv(ib) == iclst(nlstc)) cycle
          endif
          nlstc = nlstc+1
          iclst(nlstc) = iwkv(ib)
        enddo
      elseif (wordsw(opts,dc,'ic=','',j1) /= 0) then
        call nwordg(opts,0,dc//' ',1,j1,j2)
        nlstc = mkilsd(opts(j1:j2),-1,iv)
        allocate(iclst(nlstc))
        if (nlstc <= 0) call rxs(' bad or null site list : ',opts(j1:j2))
        call mkilssr(11,opts(j1:j2),nlstc,iclst,[1,nclspd])
        if (maxval(iclst(1:nlstc)) > nclspd) call rxs(' bad or null site list : ',opts(j1:j2))
        if (nlstc <= 0 .or. nlstc > nclspd) goto 98
C     elseif (wordsw(opts,dc,'all','',j1) /= 0) then
      else
        nlstc = nclspd
        allocate(iclst(nlstc))
        forall (i = 1:nlstc) iclst(i) = i
      endif
      if (nlstc == 0) then
        call rx(' asazerq : no class list specified')
        return
      endif

C --- Assign weights for each class ---
      vol = 0; avfracq = 1; qin = 0 ; qout = 0
      do  ipass = 1, 2
C       First pass: compute total volume
        do  ic = 1, nlstc
          icp = iclst(ic)
          nrclas = s_ctrl%nrc(icp)
          vsph = 4.1887902d0*s_ctrl%rmax(icp)**3
          if (wordsw(opts,dc,'add','',j1) == 0) then  ! Check if charge is too small to scale
            qnow = dsum(lmx(ic)+1,s_pot%qnu(1,icp),3)
            if (nsp == 2) qnow = qnow + dsum(lmx(ic)+1,s_pot%qnu(1+3*nl,icp),3)
            if (abs(qnow) < qmin) then
              call info2(50,0,0,'%11fskipping class %i qnow = %,6;6d',ic,qnow)
              cycle
            endif
          endif
          if (ipass == 1) then
            vol = vol + nrclas*vsph
          else
            fracq = vsph/vol
            qsph = fracq*qtot
            if (wordsw(opts,dc,'add','',j1) /= 0) then
              s_pot%qnu(1,icp) = -qsph
              if (nsp == 2) call rx('asazerq: adding charge not implemented for spin pol case yet, sorry')
              call info5(45,0,0,'%11fadding charge %,6;6d to s channel, class %i',-qsph,ic,3,4,5)
            else
              qnow = dsum(lmx(ic)+1,s_pot%qnu(1,icp),3)
              if (nsp == 2) qnow = qnow + dsum(lmx(ic)+1,s_pot%qnu(1+3*nl,icp),3)
              if (abs(qnow) < qmin) then
                call rx('asazerq: attempt to scale sphere with too few electrons')
              endif
              qin = qin + qnow
              fracq = (qnow-qsph)/qnow
              avfracq = avfracq * fracq
              call dscal(lmx(ic)+1,fracq,s_pot%qnu(1,icp),3)
              xv(1) = dsum(lmx(ic)+1,s_pot%qnu(1,icp),3)
              if (nsp == 2) then
                call dscal(lmx(ic)+1,fracq,s_pot%qnu(1+3*nl,icp),3)
                xv(1) = xv(1) + dsum(lmx(ic)+1,s_pot%qnu(1+3*nl,icp),3)
              endif
              qout = qout + xv(1)
              call info5(45,0,0,'%11fscaling charge, class %i, factor %,6;6d  q = %,6;6d -> %,6;6d',ic,fracq,qnow,xv,5)
            endif
          endif
        enddo
      enddo

      if (.not. wordsw(opts,dc,'add','',j1) /= 0) then
        fracq = avfracq**(1d0/nlstc)
        call info5(25,0,0,'%11faverage scale factor  %,6;6d  qin = %,6;6d  qout %,6;6d',fracq,qin,qout,4,5)
      endif

      return
   98 continue
      call rx('asazerq : failed to parse '//trim(opts))

      end
