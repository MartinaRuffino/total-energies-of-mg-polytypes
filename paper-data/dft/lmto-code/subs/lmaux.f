      subroutine lmaux(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_pot,s_str,s_spec,s_site,s_strn,slabl,mode)
C- Routines for auxilliary programs
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Stored:     ef def w dos idtet ipq pdos qp star wtkp wtkb swtk
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 aioxtn
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin lncol modep lpgf npadl
Ci                 npadr nclasp omax1 omax2 wsrmax ips nbasp ipc cllst
Ci                 clp clssl group initc ics idcc ipcp mxcst nrc ncomp
Ci                 pgfsl pgfvl pgord pgplp spid dclabl pos rmax lmet
Ci                 zbak ldlm nccomp lrs nesabc rmines rmaxes sclwsr
Co     Stored:     cllst clp clssl group initc ics idcc ipc ipcp ips
Co                 mxcst nrc ncomp pgfsl pgfvl pgord pgplp spid dclabl
Co                 pos rmax nbas nbasp lrs
Co     Allocated:  *
Cio    Elts passed: rmax ics lrs pgfsl ips dclabl ipc nrc ncomp lasa
Cio                ipcp idcc
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 supot asamad
Cio                asavqm asaqmp aioxtn findes
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula qss
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:eula
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat nkd nsgrp nabc ag bgv cg cy dlv gv gvq
Ci                 indxcg ips0 istab jcg kv igv igv2 kv2 pos qlv symgr
Ci                 s_sym vol awald nkq as tol nkdmx nkqmx platl platr
Ci                 qlat plat2
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr s_sym
Co     Allocated:  pos
Cio    Elts passed:pos dlv symgr ag qlv vol cg jcg indxcg plat0 qlat
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 supot asamad
Cio                asavqm asaqmp plana aioxtn findes
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  ves aamom bxc cp ddpf dddpf ddpfr dlmwt dmatk dpf
Ci                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Ci                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Ci                 rhat rnew rhos rhrmx sop shfac thetcl vdif vintr
Ci                 vrmax vshft smpot smrho smrout vconst vmtz0
Co     Stored:     ves aamom bxc cp ddpf dddpf ddpfr dlmwt dmatk dpf
Co                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Co                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Co                 rhat rnew rhos rhrmx sop shfac thetcl vdif vintr
Co                 vrmax vshft smpot smrho smrout vconst
Co     Allocated:  mad
Cio    Elts passed:pnu qnu ves aamom pp mad vrmax gibbs qc qpp pmpol
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 supot asamad
Cio                asavqm asaqmp
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  mxnbr rmax
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z rmt rham hcr lmxa name a nr lmxl kmxt p pz lfoca
Ci                 qc idmod rsma lmxb kmxv rsmv rfoca ctail etail stc
Ci                 nxi exi chfa rsmfa rhoc coreh coreq lmxf
Co     Stored:     rmt name lmxa a nr z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa rhoc
Co     Allocated:  rhoc
Cio    Elts passed:mxcst rhoc
Cio    Passed to:  asars iorsa bcast_strx iosits asamad getq gtpcor dlmq
Cio                asavqm plana shopol aioxtn
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos spec class pnu force vel pz v0 v1 bxc cpawt omg
Ci                 omgn domg dmat gc gcu gcorr sfvrtx j0 pdos rho1 rho2
Ci                 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl
Ci                 eqkkl sighh sighk sigkk tauhh tauhk taukk pihh pihk
Ci                 pikk sohh sohk sokk sighhx sighkx sigkkx tauhhx
Ci                 tauhkx taukkx pihhx pihkx pikkx thet clabel eula pl
Ci                 relax vshft ndelta delta
Co     Stored:     pnu pos pos0 force vel spec clabel pz bxc cpawt omg
Co                 omgn domg dmat gc gcu gcorr sfvrtx j0 pdos rho1 rho2
Co                 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl
Co                 eqkkl sighh sighk sigkk tauhh tauhk taukk pihh pihk
Co                 pikk sohh sohk sokk sighhx sighkx sigkkx tauhhx
Co                 tauhkx taukkx pihhx pihkx pikkx thet v0 v1 eula pl
Co                 relax vshft
Co     Allocated:  v0 v1
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  asars asars1 iorsa bcast_strx shoshl iopos iosits
Cio                plana aioxtn aioxt1
Ci Inputs
Ci   prgnam:name of main program
Ci   slabl :vector of species labels
Ci   mode  :a compound of bits, which are independent of each other
Ci         :  2**0 Show neighbors (standard mode for lmchk)
Ci         :  2**1 Plane analysis
Ci         :  2**2 Generate input to xbs program
Ci         :  2**3 Shift moments, pp's to new linearization energy
Ci         :  2**4 Interpolate core to another mesh
Ci         :  2**5 Display poles of potential functions
Ci         :  2**6 Import data from other formats
Ci         :  2**7 Find empty spheres (--findes) or resize spheres (--omax)
Cs Command-line switches
Cs   --wpos=fn   : (lmchk) write site positions to file
Cs   --Cion      : Generate the nonanalytic dipole contribution to interatomic force constants
Cs   --findes    : (lmchk) empty sphere finder
Cs   --euler     : (lmchk) print inner products between Euler angles (nc case)
Cs   --getwsr    : (lmchk) Invoke algorithm to get MT radius
Cs   --angles    : (lmchk) print angles between atom triplets
Cs   --shell     : (lmchk) Print out neighbor tables.  Several formats are available
Cs   --mino      : (lmchk) Move atoms to minimize overlaps (useful for E.S.)
Cs   --terse     : (lmchk) terse printout of sphere overlaps
Cs   --basis=    : (lmchk) check whether s_lat%pos differs from site file by translation
Cs   -bs=        : (lmxbs) see www.questaal.org/docs/input/commandline/#switches-for-lmxbs
Cs   -dup=       : (lmxbs) see www.questaal.org/docs/input/commandline/#switches-for-lmxbs
Cs   -sp         : (lmxbs) see www.questaal.org/docs/input/commandline/#switches-for-lmxbs
Cs   -spec       : (lmxbs) see www.questaal.org/docs/input/commandline/#switches-for-lmxbs
Cs   -ss=        : (lmxbs) see www.questaal.org/docs/input/commandline/#switches-for-lmxbs
Cs   --mlog      : (MPI) write MPI commands to log file
Cs   --rs=       : Controls I/O with rst file; see Command-line-options.html
Cs   --shorten   : Suppress all shortening of basis vectors
Cs   --syml      : Write symmetry lines file
Cs   --wsite     : write site file
Cs   --wsitex    : write site file
Cs   --rs        : (lmimp) Controls I/O with rst file; see Command-line-options.html
Cs   -enu=       : Not documented
Cs   --basp      : Not documented
Cs   --dumprho   : Not documented
Cu Updates
Cu   26 Feb 18 New --syml switch
Cu   16 Apr 16 New --Cion switch
Cu   19 Jun 15 --wsite option can accept ~short argument
Cu   23 Apr 15 New --wsitep option for lmaux (write a modified POSCAR file)
Cu   31 Jan 14 New --wsite option for lmaux
Cu   08 Jul 13 Replace all f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   01 Sep 11 Begin migration to f90 structures
Cu   10 Feb 10 angtab re-implemented and updated
Cu   12 Aug 08 (L. Ke) empty sphere finder
Cu   04 Nov 04 Upgrade of rsta editor
Cu   26 Jan 03 Call to angtab changed
Cu   17 May 02 Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02 Added option (--getwsr) to find MT radii
Cu   01 Mar 02 Updated Import data mode
Cu   05 Oct 01 Adapted mode 2**3 to work with lm v6.11
Cu   24 Nov 97 changed ovmin to run quickly
C ----------------------------------------------------------------------
      use structures
      implicit none
      integer mode
C ... Passed parameters
      character*(*) prgnam*8
      character*8 slabl(*)
C ... For structures
!       include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_mix)::   s_mix
      type(str_pot)::   s_pot
      type(str_str)::   s_str
C     type(str_tb)::    s_tb
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_strn) :: s_strn(*)

C ... Dynamically allocated arrays
      integer, allocatable :: iax(:)
      integer,pointer :: ips(:)
      integer, allocatable :: it(:)
C      real(8),pointer :: rham(:),rmax(:),hcr(:),symgr(:)
      real(8),allocatable :: zc(:),rmts(:)
      integer, allocatable :: ics(:)
      integer, allocatable :: ips2(:)
      integer, allocatable :: iwk(:)
      integer, allocatable :: lmxa(:),lmxb(:)
      integer, allocatable :: lock(:)
      integer, allocatable :: lockc(:)
      integer, allocatable :: ntab(:)
      real(8), allocatable :: hcr(:)
      real(8), allocatable :: pos2(:)
      real(8), allocatable :: rham(:)
      real(8), allocatable :: rmax(:)
      real(8), allocatable :: rmtc(:)
      real(8), allocatable :: rmx(:)
      real(8), allocatable :: symgr(:)
      real(8), allocatable :: zs(:)
      real(8), allocatable :: zstar(:,:,:)
      real(8), allocatable :: mion(:,:,:)
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: cij(:,:,:,:)
      character(len=8),allocatable :: clabel(:)
C ... Local parameters
      character outs*256,fnam*8,out2*256
      logical swtmp
      logical :: lmino = .false. ! If T, shift ES to minimize overlap after they are found
      logical :: liosits = .false. ! Switch, used when writing site file
      logical, parameter :: T=.true., F=.false.
      integer, parameter :: ngmx=48,NULLI=-99999, niax=10
      integer cmplat,i,ic,ifi,ip,is,j,j1,j2,k,lncol,lpbc,lpgf,m,mxcsiz,mxnbr,mxesspec,
     .  nbas,nbasp,nbaspp,nbas0,nc,nclasp,nclass,nclspp,neul,ngrp,nkd,nl,npadl,npadr,mxaddes,mxclas,nclass0,
     .  nsgrp,nsp,nspec,nttab,stdo
      integer irs(5),modep(3),iv(1000),ishores(3)
      double precision xx,alat,facrmx,facrng,avw,enu,rmaxs,wsrmax,emad,trumad,vne,vrl,sclwse
      double precision xv(10),qss(4),qlat(9),plat(3,3),vmtz(2),omax1(3),omax2(3),omaxe(3),
     .  etrms(22),eps00(3,3),vconst(3)
      character dc*1
      procedure(logical) :: cmdopt
      procedure(integer) :: iprint,getdig,lgunit,fopn,fopna,fxst,
     .  iositp,iosits,iosite,cmdoptswx,wordsw,a2vec,parg,str_pack

C ... External calls
      external aioxtn,amagnc,angtab,asamad,asars,awrit2,dcopy,dgemm,
     .         dinv33,dpcopy,dpscop,fclose,fclr,fexit,findes,icopy,
     .         info0,info5,iopos,ioxbs,makrm0,ovmin,plana,poppr,
     .         pshpr,ptr_ctrl,rxs,sclwsr,shoang,shopol,shorps,shoshl,
     .         sitepack,skipbl,spec2class,supot,symlat,to


      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      lncol = s_ctrl%lncol
      modep = s_ctrl%modep
      lpgf = s_ctrl%lpgf(1)
C     lpbc = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc = 0
      if (lpgf > 0 .and. lpgf < 10) lpbc = 11
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = s_ctrl%nclasp
      neul = s_ham%neula
      qss = s_ham%qss
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      nkd = s_lat%nkd
      mxnbr = s_str%mxnbr
      rmaxs = s_str%rmax
      allocate(zs(nspec),rmts(nspec))
      do  is = 1, nspec
        zs(is) = s_spec(is)%z
        rmts(is) = s_spec(is)%rmt
      enddo
      allocate(zc(nclasp))
      allocate(ics(nclasp),lmxb(nclasp),lmxa(nclasp),rham(nclasp),rmax(nclasp),hcr(nclasp))
      call dcopy(nclasp,s_ctrl%rmax,1,rmax,1)
      call icopy(nclasp,s_ctrl%ics,1,ics,1)
      do  ic = 1, nclasp
        is = ics(ic)
        zc(ic) = zs(is)
        rham(ic) = s_spec(is)%rham
        hcr(ic) = s_spec(is)%hcr(1)
        lmxa(ic) = s_spec(is)%lmxa
        lmxb(ic) = s_spec(is)%lmxb
      enddo

      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      stdo = lgunit(1)

      if (cmdopt('--shorten',9,0,outs)) then
        if (.not. cmdopt('--shorten=no',12,0,outs)) then
        call shorps(nbasp,plat,modep,s_lat%pos,s_lat%pos)
        do i = 1, nbasp
          s_site(i)%pos(:) = s_lat%pos(:,i)
        enddo
      endif
      endif

C ... Build symmetry lines file
      if (cmdopt('--syml',6,0,outs)) then
        i = str_pack('syml',-2,s_strn,out2)
        if (outs(7:7) == ' ') then
          outs(7:) = out2
        endif
        call info0(10,1,0,' MKSYML: generating symmetry line file ...')
        call mksyml(outs(7:),plat)
        call rx0(' wrote symmetry lines to syml file')
      endif

C ... Ionic contribution to force constant matrix
      if (cmdopt('--Cion',6,0,outs)) then
        allocate(cij(3,3,nbas,nbas),zstar(3,3,nbas),mion(3,3,nbas),
     .            it(max(nbas,3)),wk(max(nbas,3)))
        call dpzero(zstar,size(zstar)); call dpzero(cij,2*size(cij))
        call dpzero(mion,size(mion))
        qss(1:3) = (/.1d0,.02d0,.03d0/)
        ip = 3
        if (cmdopt('-q=',ip,0,outs)) then
          call skipbl(outs,len(outs),ip)
          i = parg(' ',4,outs,ip,len(outs),', ',2,3,it,qss)
          if (i /= 3) call rx('failed to parse '//trim(outs))
        endif
        eps00 = 0; eps00(1,1) = 1d0; eps00(2,2) = 1d0; eps00(3,3) = 1d0
        ip = 3
        print *,' eps00', eps00
        if (cmdopt('-z=',ip,0,outs)) then
          call skipbl(outs,len(outs),ip)
          i = parg(' ',4,outs,ip,len(outs),', ',2,nbas,it,wk)
          if (i /= nbas) call rx('failed to parse '//trim(outs))
          do  i = 1, nbas
            zstar(1,1,i) = wk(i); zstar(2,2,i) = wk(i); zstar(3,3,i) = wk(i)
          enddo
        endif
        ip=3
        if (cmdopt('-M=',ip,0,outs)) then
          call skipbl(outs,len(outs),ip)
          i = parg(' ',4,outs,ip,len(outs),', ',2,nbas,it,wk)
          if (i /= nbas) call rx('failed to parse '//trim(outs))
          do  i = 1, nbas
            mion(1,1,i) = 1d0/sqrt(wk(i)*1822.9)
            mion(2,2,i) = 1d0/sqrt(wk(i)*1822.9)
            mion(3,3,i) = 1d0/sqrt(wk(i)*1822.9)
          enddo
        endif

        ip = 4
        call epsdipole(ip,nbas,s_lat,eps00,zstar,mion,qss)
C       call strxion(nbas,s_lat,eps00,qss,cij)
        call rx0('done')
      endif

C ... Read from restart file
      if (cmdopt('--rs=',5,0,outs)) then
        irs(1) = IAND(s_ctrl%lrs,7)
        if (irs(1) > 0) then
          ifi = fopna('rsta',-1,0)
          call asars(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .      s_pot%pnu,s_pot%qnu,.false.,ifi)
          call fclr('rsta',ifi)
C         call shoctl(s_ctrl,s_spec,s_pot,F,stdo)
C         call rx('done')
        endif
      endif

C ... Override given omax1 if command-line argument is present
      omax1 = s_ctrl%omax1
      omax2 = s_ctrl%omax2
      wsrmax = s_ctrl%wsrmax
      if (cmdopt('--omax=',7,0,outs)) then
        i = 7
        is = a2vec(outs,len_trim(outs),i,4,', ',2,3,3,iv,omax1)
        if (is < 0) call rx('suctrl: failed to parse '//trim(outs))
        if (is < 2) omax1(2) = omax1(1)
        if (is < 3) omax1(3) = omax1(2)
      endif
      s_ctrl%omax1 = omax1

C --- Neighbor tables and sphere overlaps ---
      if (getdig(mode,0,2) /= 0) then
      if (rmaxs <= 0d0) then
        rmaxs = 2.7d0*avw
        call info5(30,0,0,'%1fUse default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat',
     .    rmaxs,rmaxs/avw,rmaxs/alat,0,0)
      endif

C ... Get neighbor table iax for each atom in the cluster
      if (lpbc == 0) then  ! 3 dimensions
        i = 3
        j = -1
      elseif (lpbc == 1 .or. lpbc == 11) then
        i = 2
        j = 1
      else
        call rx('pairs:: not implemented for lpbc>1')
      endif
      mxcsiz = s_str%mxnbr
      call pshpr(iprint()-20)
      allocate(ntab(nbasp+1))
      call pairs(nbas,nbasp,alat,plat,[rmaxs/2],s_lat%pos,
     .  [-1],i,j,s_ctrl%pgfsl,nttab,ntab,iax,mxcsiz)
      call poppr

C --- Print out a few lattice lattice translation vectors ---
      j = 6
      if (cmdopt('--slat',j-1,0,outs)) then
      if (iprint() >= 10) then
        call info0(10,1,0,' LMCHK:  print multiples of plat%N'//
     .    '  i1  i2  i3%7fx%11fy%11fz%11flen')
        do  i = -2, 2
        do  j = -2, 2
        do  k = -2, 2
        xx = 0
        do  m = 1, 3
          xv(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
          xx = xx + xv(m)**2
        enddo
        xx = dsqrt(xx)
        print 368, i,j,k, xv(1), xv(2), xv(3), xx
  368   format(3i4, 3f12.7, 1x, f12.5)
        enddo
        enddo
        enddo
      endif
      endif

C --- Find sphere overlaps or resize sphere ---
      j = 9
      if (cmdopt('--getwsr',j-1,0,outs)) then
        call info(10,1,0,' ... Make sphere radii',0,0)
C        xx = dglob('lrel',1d0,1)
C        xx = dglob('nsp',1d0,1)
C       Initial estimate for sphere radii: overlapping atom potentials
        allocate(lock(nspec))
        do  i = 1, nspec
          lock(i) = IAND(s_spec(i)%mxcst,2)
        enddo
        if (lpbc == 0) then
          i = 3
        elseif (lpbc == 1 .or. lpbc == 11) then
          i = 2
        else
          call rx('LMAUX: not implemented for lpbc>1')
        endif
        call makrm0(101,nspec,nbas,alat,plat,s_lat%pos,slabl,
     .    s_ctrl%ips,modep,lock,zs,rmts)

C   ... Scale sphere radii satisfying constraints: scale rmts(1:nspec)
        wsrmax = s_ctrl%wsrmax
        call sclwsr(20,nbas,nbasp,nspec,alat,plat,s_lat%pos,s_ctrl%ips,
     .    modep,slabl,zs,lock,1d0,wsrmax,omax1,omax2,rmts)
        call spec2class(s_spec,nspec,-1,'-rmt',1,xx,rmts)
C        do  is = 1, nspec
C          s_spec(is)%rmt = rmts(is)
C        enddo
        nclspp = max(2*nclasp-nclass,nspec)
        deallocate(rmax)
        allocate(rmax(nclspp))
        call spec2class(s_spec,nclasp,ics,'rmt',1,xx,rmax)
C        do  ic = 1, nclspp
C          is = ics(ic)
C          rmax(ic) = rmts(is)
C        enddo

      endif

C --- Show neighbors by shell ---
      outs = ' '; j = 8
      if (.not. associated(s_ham%eula)) allocate(s_ham%eula(1,1))
      if (cmdopt('--shell',j-1,0,outs)) then
        call shoshl(outs(j:),nbasp,s_site,s_lat%pos,alat,plat,mxnbr,zc,
     .    slabl,s_ctrl%dclabl,s_ctrl%ips,s_ctrl%ipc,s_pot%ves,
     .    s_ham%eula,nclass)
      endif

C --- Show angles between neighbors ---
      j = 9
      if (cmdopt('--angles',j-1,0,outs)) then
       call shoang(outs(j:),nbas,s_lat%pos,plat,mxnbr,slabl,s_ctrl%ips,zs)
      endif

C --- Check whether basis s_lat%pos differs from site file by translation ---
      j = 9
      if (cmdopt('--basis=',j-1,0,outs)) then
        fnam = outs(j:)
        call info(20,1,0,' checking whether basis equivalent to file '
     .    //fnam//'...',0,0)
        j = iosits(8070,3d0,0,fnam,ifi,slabl,alat,plat,nbas,nspec,
     .    s_spec,s_site)
        allocate(pos2(3*nbas))
        allocate(ips2(nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,pos2)
        call sitepack(s_site,1,nbas,'spec',1,ips2,xx)
        allocate(symgr(9*ngmx))
        call symlat(plat,ngrp,symgr,j)
        j = cmplat(nbas,plat,plat,ngrp,symgr,s_ctrl%ips,s_lat%pos,
     .    ips2,pos2)
        deallocate(pos2,ips2)
        call fexit(j,1,' Exit %i lmchk --basis= ...',j)
      endif

C ... Write positions in Cartesian coordinates and as multiples plat
      if (iprint() >= 50) then
      write(stdo,357)
  357 format(/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
C     qlat = (plat+)^-1
      call dinv33(plat,1,qlat,xx)
      do  i = 1, nbas
        call dpscop(s_lat%pos,xv,3,3*i-2,1,1d0)
C       posp+ = (plat)^-1 pos+
        call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
        ip = s_ctrl%ips(i)
        print 345, i, slabl(ip), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345   format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
      enddo
      endif

C --- Print overlaps, optionally minimize wrt spec'd sites ---
      outs = ' '
      i = 6
      swtmp = cmdopt('-mino',5,0,outs)
      swtmp = cmdopt('--mino',6,0,outs)
      if (swtmp) i = 7
      j = 1
      if (iprint() < 30) j = 0
      call ovmin(outs(i:),nbas,nbasp,alat,plat,rmax,hcr,
     .  s_ctrl%clabl,s_ctrl%ipc,modep,zc,ntab,iax,s_lat%pos,j)
      forall (i = 1:nbasp) s_site(i)%pos(:) = s_lat%pos(:,i)
C ... Write positions to file
      if (cmdopt('--wpos=',7,0,outs)) call iopos(T,0,outs(8:),nbasp,s_lat%pos,s_site)
      deallocate(ntab)

C ... Inner products between Euler angles
      j = 8
      if (mod(lncol,2) == 1 .and. cmdopt('--euler',j-1,0,outs)) then
        call amagnc(nbas,nl,s_ctrl%ipc,xx,1,s_pot%qnu,s_ham%eula,
     .    neul,1,xv,s_pot%aamom,xx)
        call dinv33(plat,1,qlat,xx)
        print '(1x)'
        call angtab(outs(j:),nbas,s_lat%pos,alat,rmax,qss,qlat,
     .    s_lat%dlv,nkd,s_ctrl%ipc,slabl,s_ctrl%ips,zs,neul,s_ham%eula,
     .    nl,s_pot%qnu)
      endif

      nclass = nspec
      if (associated(s_ctrl%clabl)) deallocate(s_ctrl%clabl)
      allocate(s_ctrl%clabl(nclass))
      forall (i = 1:nclass) s_ctrl%clabl(i) = slabl(i)

      liosits = .true.  ! If write site file, use iosits
      goto 99  ! Cover possible --wsite

      endif ! Neighbor tables and sphere overlaps

C --- Plane analysis branch ---
      if (getdig(mode,1,2) /= 0) then
        nbas = s_ctrl%nbas
        nbasp = s_ctrl%nbasp
        nl = s_ctrl%nl
        call supot(1,s_ctrl,s_lat,s_pot)
        vne = s_bz%semsh(8); vrl = vne
        if (fxst('vshft') /= 0) then
C         if (procid == master) then
          ifi = fopn('vshft')
          ip = 1+4 ; if (lpgf /= 0) ip = 2+4
          call iovshf(nbasp,ip,'read potential shifts from file vshft',
     .      xx,xv,vconst,s_pot%vshft,ifi)
C         endif
C         call mpibc1(s_pot%vshft,nbas,4,.false.,'lmasa','vshft')
        endif
        call asamad(s_ctrl,s_pot,s_lat,s_spec,0,
     .    s_pot%pnu,s_pot%qnu,vrl,s_pot%ves,emad,trumad,vmtz,etrms)
        call plana(npadl,npadr,nbaspp,slabl,
     .    s_lat,s_spec,s_site,s_pot%ves,s_pot%vshft,vconst,vrl,s_pot%pnu,s_pot%qnu)
        return
      endif

C --- Generate input file to xbs program ---
      if (getdig(mode,2,2) /= 0) then
      ifi = fopn('XBS')
      facrmx = .5d0
      ip = 4
      if (cmdopt('-bs=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrmx)
      endif
      facrng = 1d0
      ip = 4
      if (cmdopt('-ss=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrng)
      endif
C ... Copy wsr*facrmx into rmax, wsr*facrng into rham (if nonzero)
      allocate(rmx(nclass))
      call dpcopy(rmax,rmx,1,nclass,facrmx)
C ... Copy wsr*facrng into rham, if zero
      do  i = 1, nclass
        is = ics(i)
        if (rham(i) == 0 .or. rham(i) == NULLI)
     .    rham(i) = facrng*rmts(is)
      enddo
      if (iprint() >= 20) then
        call awrit2('%N ball size = %d * sphere size;  '//
     .    'def stick length = %d * sum sphere sizes',
     .    ' ',80,stdo,facrmx,facrng)
      endif
      if (cmdopt('-spec',5,0,outs) .or. cmdopt('--spec',6,0,outs)) then
        nc = nspec
        ips => s_ctrl%ips
      else
        nc = nclass
        ips => s_ctrl%ipc
      endif
      call ioxbs(ifi,nbas,nc,alat,plat,rham,rmx,s_ctrl%dclabl,ips,zc,s_lat%pos)
      call fclose(ifi)
      deallocate(rmx)
      endif

C --- Shift pp's (and optionally) moments by enu ---
C     pp's are remade using the potential if available.
C     use -enu=val to shift all enu's to val.
C     Use -mom if to save shifted moments.  Potential NOT remade.
      if (getdig(mode,3,2) /= 0) then
      call rx('lmaux mode 3 needs to be updated')
      ip = 5
      if (cmdopt('-enu=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,enu)
        if (i == -1) call rxs('LMSHF: failed to parse ',outs)
      else
        call rx('LMSHF: missing argument -enu=val')
      endif
      endif

C --- Interpolate core to another mesh ---
      if (getdig(mode,4,2) /= 0) then
        call rx('patch clabl for call to coritp')
C       call coritp(nclass,nsp,s_ctrl%dclabl,w(onrmsh),w(oamsh),rmax)
      endif

C --- Display poles of potential functions ---
      if (getdig(mode,5,2) /= 0) then
        call rx('lmaux update call to iostr')
        call shopol(nl,nclass,nsp,s_spec,ics,s_pot%pp)
      endif

C --- Import data in other formats ---
      if (getdig(mode,6,2) /= 0) then
        call rx('update interface to import')
        call aioxtn(s_ctrl,s_spec,s_site,s_lat,s_bz,slabl,s_pot%pnu,
     .    s_pot%qnu)

C       Output to restart file
        if (cmdopt('--rs=',5,0,outs)) then
          irs(2) = IAND(s_ctrl%lrs,8+16)/8
          if (irs(2) > 0) then
            ifi = fopna('rsta',-1,0)
            call asars(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .        s_pot%pnu,s_pot%qnu,.false.,-ifi)
            call fclr('rsta',ifi)
          endif
        endif
      endif

C --- Empty sphere finder ---
      if (getdig(mode,7,2) /= 0 .and. cmdopt('--findes',8,0,outs)) then
        mxesspec = NULLI; mxaddes = 0
        nsgrp = s_lat%nsgrp
        allocate(lock(nspec)); call iinit(lock,nspec)
        forall (i = 1:nspec) lock(i) = IAND(s_spec(i)%mxcst,2)
        dc = outs(9:9)
        k = wordsw(outs,dc,'rmin=','',j1) + 5
        if (k > 5) then
          if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,iv,s_ctrl%rmines) < 1)
     .      call rx('lmchk failed to parse '//trim(outs))
        endif
        k = wordsw(outs,dc,'nspmx=','',j1) + 6
        if (k > 6) then
          if (a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,1,iv,mxesspec) < 1)
     .      call rx('lmchk failed to parse '//trim(outs))
        endif
        mxaddes = 0
        k = wordsw(outs,dc,'nesmx=','',j1) + 6
        if (k > 6) then
          if (a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,1,iv,mxaddes) < 1)
     .      call rx('lmchk failed to parse '//trim(outs))
        endif
        k = wordsw(outs,dc,'wsrmx=','',j1) + 6
        if (k > 6) then
          if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,iv,wsrmax) < 1)
     .      call rx('lmchk failed to parse '//trim(outs))
        endif
        lmino = wordsw(outs,dc,'mino','',j1) > 0
        i = cmdoptswx('--findes','1spec','') + 5
        if (i > 5) mxesspec = 1
        if (wordsw(outs,dc,'spec','',j1) > 0) then
          s_ctrl%nrc = 0
          nclass = nspec
          do  i = 1, nbas
            is = s_ctrl%ips(i)
            s_ctrl%ipc(i) = is
            s_site(i)%clabel = slabl(is)
            s_ctrl%nrc(is) = s_ctrl%nrc(is) + 1
          enddo
          do  i = 1, nspec
            s_ctrl%clabl(i) = slabl(i)
            s_ctrl%spid(i) = slabl(i)
            zc(i) = zs(i)
          enddo
        endif

        ishores = 0
        k = wordsw(outs,dc,'shorten=','',j1) + 8
        if (k > 8) then
          is = a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,3,iv,ishores)
          if (is /= 3) call rx('suctrl: failed to parse '//trim(outs))
        endif

        omaxe = NULLI
        k = wordsw(outs,dc,'omax=','',j1) + 5
        if (k > 5) then
          is = a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,3,iv,omaxe)
          if (is < 0) call rx('suctrl: failed to parse '//trim(outs))
          if (is < 2) omaxe(2) = omaxe(1)
          if (is < 3) omaxe(3) = omaxe(2)
        endif
        sclwse = NULLI
        k = wordsw(outs,dc,'sclwsr=','',j1) + 7
        if (k > 7) then
          if (a2vec(outs,len_trim(outs),k,4,', '//dc,3,2,1,iv,sclwse) < 1)
     .      call rx('lmchk failed to parse '//trim(outs))
        endif
        k = wordsw(outs,dc,'lock=','',j1) + 5
        if (k > 5) then
          is = a2vec(outs,len_trim(outs),k,2,', '//dc,3,2,nspec,iv,lock)
          if (is < 0) call rx('suctrl: failed to parse '//trim(outs))
          do  i = 1, is
            if (lock(i) == 1) lock(i) = 2
          enddo
        endif

C       Retain for backward compatibility .. same functionality as --findes~nesmx
        i = 9
        if (cmdopt('--nescut=',i,0,outs)) then
          is = a2vec(outs,len(outs),i,2,', ',2,2,1,iv,mxaddes)
          if (is /= 1) call rx('findes failed to parse '//trim(outs))
        endif

        if (allocated(zc)) deallocate(zc)
        mxclas = 5*nbas
        allocate(zc(mxclas),rmtc(mxclas))
        allocate(lockc(mxclas)); call iinit(lockc,mxclas)
        call ptr_ctrl(s_ctrl,2,'nrc',mxclas,0,0,xx)
        if (nclass >= mxclas) call rx('lmaux: increase mxclas')
        call spec2c(nspec,nclass,ics,rmts,rmtc,zs,zc,lock,lockc)
        allocate(iwk(size(lockc))); call icopy(size(lockc),lockc,1,iwk,1)
        allocate(clabel(nclass))
        forall (ic = 1:nclass) clabel(ic) = s_ctrl%clabl(ic)
        deallocate(s_ctrl%clabl); allocate(s_ctrl%clabl(mxclas))
        forall (ic = 1:nclass) s_ctrl%clabl(ic) = clabel(ic)

        if (omaxe(1) /= NULLI) then
          call info2(2,1,0,' Temporary scaling of spheres while finding ES: use omax =%s,%3:1;1d%%',100*omaxe,2)
          call sclwsr(100,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,modep,s_ctrl%clabl,zc,
     .      lockc,1d0,wsrmax,omaxe,omax2,rmtc)
          s_ctrl%omax1 = omaxe
          call icopy(size(lockc),iwk,1,lockc,1)
!         call ovlchk(nbas,nbas,s_lat%pos,alat,rmtc,rmtc,s_ctrl%clabl,s_ctrl%ipc,modep,plat,xv(1),xv(2))
        endif

        nbas0 = nbas
        nclass0 = nclass
        call findes(1,mxaddes,mxesspec,s_ctrl,s_lat,alat,nbas,nclass,nl,
     .    s_ctrl%nrc,mxclas,nsgrp,plat,s_lat%symgr,s_lat%ag,ishores,lockc,rmtc,zc)
        call icopy(size(lockc),iwk,1,lockc,1)

!       if (lesend .and. zc(nclass0) <= 0.5d0) nclass0 = nclass0-1 ! Maybe merge last species with new ES
        if (mxesspec /= NULLI .and. nclass-nclass0 > mxesspec) then
          do  i = nbas0+1, nbas
            if (s_ctrl%ipc(i) > nclass0+mxesspec) then
              s_ctrl%ipc(i) = nclass0+mxesspec
            endif
          enddo
          rmtc(nclass0+mxesspec) = min(rmtc(nclass0+mxesspec),rmtc(nclass)) ! choose the smallest size
          nclass = nclass0 + mxesspec
        endif

        forall (i = 1:nclass) s_site(i)%clabel = s_ctrl%clabl(i)

        if (omaxe(1) /= NULLI) then
          s_ctrl%omax1 = omax1
          if (sclwse /= NULLI) s_ctrl%sclwsr = sclwse
          call info2(2,0,0,' rescale spheres use omax =%s,%3:1;1d%%  sclwsr = %d',100*omax1,s_ctrl%sclwsr)
          k = s_ctrl%sclwsr/10
          call sclwsr(10*k+1,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,modep,s_ctrl%clabl,zc,
     .      lockc,1d0,wsrmax,omax1,omax2,rmtc)
          call icopy(size(lockc),iwk,1,lockc,1)
        endif

        if (lmino) then
          mxcsiz = 0
          call pshpr(iprint()-20)
          allocate(ntab(nbas+1))
          call pairs(nbas,nbas,alat,plat,[5d0],s_lat%pos,[-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
          call poppr

          call ovmin('~z',nbas,nbas,alat,plat,rmtc,rmtc,s_ctrl%clabl,s_ctrl%ipc,modep,zc,ntab,iax,s_lat%pos,1)
        endif

C   ... Unshear site positions and write basis to file poses
        call info0(20,1,0,' Writing new basis to file poses ...')
        call prpos(nbas,nclass,nclass0,s_ctrl%clabl,s_ctrl%ipc,
     .    s_lat%plat0,s_lat%qlat,s_lat%pos,lmxb,rmtc,zc)

        deallocate(lockc,zc,rmtc)
        allocate(wk(3*nbas),ips2(3*nbas))
        call dvset(wk,1,size(wk),0d0)
        call ivset(ips2,1,size(ips2),0)
        liosits = .false.        ! If write site file, use iosite

      elseif (cmdopt('--omax=',7,0,outs)) then
        wsrmax = s_ctrl%wsrmax
        mxclas = nclass
        allocate(lock(nspec))
        forall (i = 1:nspec) lock(i) = IAND(s_spec(i)%mxcst,2)
        deallocate(zc); allocate(zc(mxclas),rmtc(mxclas))
        allocate(lockc(mxclas)); call iinit(lockc,mxclas)
        call ptr_ctrl(s_ctrl,2,'nrc',mxclas,0,0,xx)
        call spec2c(nspec,nclass,ics,rmts,rmtc,zs,zc,lock,lockc)

        call info2(2,1,0,' Rescale sphere radii with omax =%s,%3:1;1d%%',100*omax1,2)
        call sclwsr(0,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,modep,s_ctrl%clabl,zc,
     .    lockc,1d0,wsrmax,omax1,omax2,rmtc)
!      call ovlchk(nbas,nbas,s_lat%pos,alat,rmtc,rmtc,s_ctrl%clabl,s_ctrl%ipc,modep,plat,xv(1),xv(2))

      endif

C --- Write site file ---
   99 continue
      if (cmdopt('--wsite',7,0,out2)) then
      outs = ' '
      if (out2(1:8) == '--wsitep') then
        if (iositp(1,s_lat,s_spec,s_site,j) < 0) call rx('failed to write POSCAR')
        call rx0('done writing file POSCAR')
      elseif (out2(1:8) == '--wsitex') then
        j = 1000*(1+2+4+8+32) + 1 + 10
        j2 = 9
      elseif (out2(1:7) == '--wsite') then
        j = 1000*(1+2+4+8+32) + 1
        j2 = 8
      endif
      dc = out2(j2:j2)
      if (dc == '=') then
        outs = out2(j2+1:)
      elseif (dc /= ' ') then
C   ... Return here to resume parsing for arguments
   10   continue
        j2 = j2+1
        if (out2(j2:j2) == dc) goto 10
        j1 = min(len(out2),j2)
        call nwordg(out2,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (out2(j1:j1+4) == 'short')  then
            j = j - 1000*iand(j/1000,32)
          elseif (out2(j1:j1+2) == 'fn=')  then
            outs = out2(j1+3:j2)
          else
            goto 11
          endif
          goto 10
   11     continue
          call rxs('lmaux: failed to parse --wsite option', out2)
        endif
      endif
      if (outs == ' ') outs = 'site'
      if (liosits) then
        if (iosits(j,3d0,0,outs,ifi,s_ctrl%clabl,alat,plat,nbas,nclass,
     .    s_spec,s_site) < 0) call rx('failed to write site file')
      else
        ic = iosite(j,3d0,0,trim(outs),ifi,s_ctrl%clabl,alat,s_lat%plat0,nbas,
     .    nclass,s_lat%pos,wk,wk,wk,s_ctrl%ipc,ips2,ips2)
      endif

      endif

      end

       subroutine spec2c(nspec,nclass,ics,rmts,rmtc,z,zc,lock,lockc)
C- Copy species data to class data
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspec
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   rmtc  :rmt by species
Ci   z     :z by species
Co Outputs
Co   rmtc  :rmt by class
Co   zc    :Z by class
co   lockc :lock by class
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Aug 08
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nspec,nclass,ics(nclass),lock(nspec),lockc(nclass)
      double precision rmts(nspec),rmtc(nclass),z(nspec),zc(nclass)
C ... Local parameters
      integer j,k

      do  k = 1, nclass
        j = ics(k)
        rmtc(k) = rmts(j)
        zc(k) = z(j)
        lockc(k) = lock(j)
C       if (iprint() > 60) write(*,310) k,rmtc(k)
      enddo

C 310 format(1x,'class ',I3,T15,'rmt = ',f10.7)

      end

