C#define INTEL_IFORT
      subroutine st2ary(s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,
     .  s_spec,s_site,s_str,s_move,
     .  nel,el,elink,
     .  alfsi,lswtch,nit,
     .  smear,dqval,etol,qrtol,seref,nfile,ntx,ndx,
     .  alat,plat,nft1,nft2,nft3,dx,dy,dz,
     .  atid,spid,z,amass,rmt,rsm,rsmfa,rint,rcut,rham,nphi,lphi,tphi,
     .  colxbs,radxbs,
     .  nxi,lxi,exi,lmxl,lmxa,pin,idmod,piz,idmoz,lsc,nbas,nbasp,nspec,
     .  ipclas,nr,arad,qin,
     .  pos,mpole,dpole,n0,nsx,nbx,
     .  as,ewtol,
     .  ntcpm,tcfpm,
     .  lrx,ldyn,tstep,tskip,mdtime,temp,taup,itmxrl,step,
     .  xtol,gtol,lrs,ifrlx,vsav,
     .  nvarms)
C- Convert data in lm structures to primitive array data for lmce lmcfit
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters

C ... Dimensioning parameters
      integer n0,nsx,nbx,nsmx
      parameter( nsmx=200 )

C ... Options, structure, basis
      integer nbas,nbasp,nspec,nsp,nel,lsc
      double precision plat(3,3),alat,el(10),alfsi(2),elink

C ... For iterations
      integer nmix,nit,nvarms
      double precision etol,qrtol(2),smear,dqval

C ... Switches
      integer lswtch(20)
      logical tpan

C ... Ewald sums
      double precision ewtol,as,alat0
      integer nkdmx,nkrmx

C ... Species parameters
      integer nphi(nsx),lphi(n0,nsx),
     .  nxi(nsx),lxi(n0,nsx),lmxa(nsx),lmxl(nsx),
     .  nr(nsx),idmod(n0,nsx),idmoz(n0,nsx)
      double precision eref(nsmx),seref,tphi(n0,5,nsx),arad(nsx),
     .  rmt(nsx),rsm(nsx),rsmfa(nsx),rint(nsx),rcut(nsx),rham(nsx),
     .  z(nsx),exi(n0,nsx),amass(nsx),
     .  vsav(3,nbx),colxbs(3,nsx),radxbs(nsx)
      double precision pin(n0,2,nsx),piz(n0,2,nsx),qin(n0,2,nsx)
      character*8 spid(nsx),atid(nbx),label

C ... Site parameters
      double precision pos(3,1),mpole(1),dpole(3,1)
      integer oipc,oips,ipclas(1)

C ... CD integration
      integer nabc(3),nft1,nft2,nft3
      double precision dabc(3),dx,dy,dz

C ... Tables, two-center fit
      integer nbisi(3),nalf,ncupl,ndust
      double precision adec,wztcf
      integer nfile,ntx,ndx,ntcpm
      double precision tcfpm(10,1)

C PBE kappa :
      double precision uk
      data uk /0.804d0/

C ... Dynamics and statics
      integer lrx,ldyn,lrs(5),ifrlx(3,*),itmxrl,tskip
      double precision tstep,mdtime,step,xtol,gtol,temp,taup,mdprm(7)

C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
      type(str_mix)::   s_mix
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_str)::   s_str
      type(str_move)::  s_move

C ... Local parameters
      integer i1mach,i,j,ib,is,iprint,lstart,ls,lxcf,lxcg,lmxb,lmxf
C#ifndefC LINUX_PGI | INTEL_IFORT
C      integer setsp,setrel,setxcf,setxcg
C#else
      integer lrel,lsp,lxc,lxg,lxcfun,lgcfun
      common /setrel/ lrel
      common /setsp/  lsp
      common /setxc/  lxc,lxg
C#endif
      logical iojob,cmdopt,a2bin
      double precision vtol(5),zbak(2),dasum
      character strn*80,s1*40,s2*80

C ... Heap, static
      integer w(1)
      common /w/ w

C --- Defaults ---
      data eref /nsmx*0d0/

      if (nsmx < nsx) call rx('st2ary: nsx gt nsmx')

      call dpzero(tphi,n0*5*nsx)
      call iinit(nr,nsx)
      call dpzero(arad,nsx)
      call iinit(lswtch,20)
C#ifndefC LINUX_PGI | INTEL_IFORT
C      j = setxcf(1)
C      j = setxcg(0)
C#else
      lxc = 1
      lxg = 0
C#endif

C ... species colours and radii for xbs
      call dcopy(3*nsx,-1d0,0,colxbs,1)
      call dcopy(nsx,-1d0,0,radxbs,1)

      nfile = 2
      ntx = 40
      ndx = 100000
      nsp = 1
      nel = 10
      elink = 0
      alfsi(1) = 0
      alfsi(2) = 0
      lrx = 0
      ldyn = 0
      tstep = .5d0
      tskip = 1
      mdtime = 1d6
      itmxrl = 1
      step = .05
      xtol = 1d-3
      gtol = 1d-3
      call icopy(3*nbx,1,0,ifrlx,1)
      qrtol(1) = 1d-5
      qrtol(2) = 1d-5
C      do  5  i = 1, 3*nbx
C    5 ifrlx(i,1) = 1
C lrs(1): no iter between rho dump; (2) not used
C    (3): use positions read in this file; (4) ditto, pnu (5) ditto v0
      lrs(1) = 1
      lrs(2) = 0
      lrs(3) = 0
      lrs(4) = 0
      lrs(5) = 0
      lsc = 0
      dqval = 0d0
C --- TCF parameters
C --- NBISI:
      tcfpm(3,1) = 8
      tcfpm(4,1) = 24
      tcfpm(5,1) = 24
C --- NALF:
      tcfpm(6,1) = 10
C --- ADEC:
      tcfpm(7,1) = 1
C --- WZTCF:
      tcfpm(8,1) = 0
C --- NCUPL:
      tcfpm(9,1) = 8
C --- NDUST:
      tcfpm(10,1) = 2
      seref = 0d0
      nmix = 0
      etol = .1d-3
      as = 2
      ewtol = 1d-6

C --- Initialise velocites to zero
      call dpzero(vsav,3*nbx)

C --- set point multipole moments to zero ---
      call dpzero(mpole,nbx)
      call dpzero(dpole,3*nbx)

C --- Disable two-panel
      tpan = .false.
C     if (lgors('ctrl lham,1',sctrl)) then
      if (IAND(s_ctrl%lham,1) /= 0) then
        call rx(' MC no longer implements two panels')
      endif
      call dcopy(n0*2*nsx,0d0,0,piz,1)

C --- ldyn and lrx
      mdprm = s_ctrl%mdprm
      itmxrl = s_ctrl%nitmv
      i = nint(mdprm(1))
      if (i > 0 .and. i < 4) then
        ldyn = i
        tstep = mdprm(2)
        temp = mdprm(3)
        taup = mdprm(4)
        mdtime = mdprm(5)
      endif
      if (i > 3) lrx = i
      call rxx(lrx > 0 .and. ldyn > 0,
     .  'ST2ARY: set either DYN or RELAX to zero')

C --- Integer switches (lswtch). Only 2,3,5 survive to version 2
C     if (lgors('ctrl lrel,-1',sctrl)) lswtch(2) = 1
      if (mod(s_ctrl%lrel,10) /= 0) lswtch(2) = 1
C     Freeze all cores, phi, phidot, vint'l
C     if (lgors('ctrl lbas,16',sctrl)) lswtch(3) = 1
      if (IAND(s_ctrl%lbas,16) /= 0) lswtch(3) = 1
C     XCQS: do XC by FFT
C     if (lgors('ctrl lcd,64',sctrl)) lswtch(5) = 1
      if (IAND(s_ctrl%lcd,64) /= 0) lswtch(5) = 1

C --- TCF coefficients: NBISI, NALF, NCUPL, NDUST, ADEC, WZTCF
      nbisi = s_str%nbisi
      nalf = s_str%nalf
      ncupl = s_str%ncupl
      ndust = s_str%ndust
      adec = s_str%adec
      wztcf = s_str%wztcf
      ntcpm = 1
      tcfpm(3,ntcpm)  = nbisi(1)
      tcfpm(4,ntcpm)  = nbisi(2)
      tcfpm(5,ntcpm)  = nbisi(3)
      tcfpm(6,ntcpm)  = nalf
      tcfpm(7,ntcpm)  = adec
      tcfpm(8,ntcpm)  = wztcf
      tcfpm(9,ntcpm)  = ncupl
      tcfpm(10,ntcpm) = ndust

C --- XC functionals
      lxcfun = mod(s_ctrl%lxcf,100)
      lgcfun = s_ctrl%lxcf/100
      if (lgcfun == 2) lxcfun = 3
      if (lgcfun == 3) lxcfun = 3
      if (lgcfun == 4) lxcfun = 3
      if (lxcfun == 3) call setuk(uk)
C#ifndefC LINUX_PGI | INTEL_IFORT
C      j = setxcg(lgcfun)
C      j = setxcf(lxcfun)
C#else
      lxg = lgcfun
      lxc = lxcfun
C#endif

C --- Unpack structures
      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      zbak = s_ctrl%zbak
      nbasp = s_ctrl%nbasp
      dqval = zbak(1)
      nit = s_ctrl%maxit
      vtol(1:3) = s_ctrl%tol(1:3)
      etol = vtol(1)
      qrtol(1) = vtol(2)
      qrtol(2) = vtol(3)
      alat = s_lat%alat
      plat = s_lat%plat
      nabc = s_lat%nabc
      el(1:6) = s_ham%kmto
      nel = s_ham%nmto
      alfsi = s_ham%alfsi
      dabc = s_ham%dabc
      smear = s_bz%w
C --- Species specific
      do  is = 1, nspec
        z(is) = s_spec(is)%z
        rmt(is) = s_spec(is)%rmt
        arad(is) = s_spec(is)%a
        nr(is) = s_spec(is)%nr
        rsm(is) = s_spec(is)%rg
        rint(is) = s_spec(is)%rint
        rcut(is) = s_spec(is)%rcut
        rham(is) = s_spec(is)%rham
        spid(is) = s_spec(is)%name
        pin(1:n0,1:nsp,is) = s_spec(is)%p(1:n0,1:nsp)
        qin(1:n0,1:nsp,is) = s_spec(is)%q(1:n0,1:nsp)
        eref(is) = s_spec(is)%eref
        lmxf = s_spec(is)%lxi
        exi(1:n0,is) = s_spec(is)%exi
        colxbs(1:3,is) = s_spec(is)%colxbs
        radxbs(is) = s_spec(is)%radxbs
        call lx2vec(lmxf,0,nxi(is),lxi(1,is))
        lmxb = s_spec(is)%lmxpb
        lmxa(is) = s_spec(is)%lmxa
        lmxl(is) = s_spec(is)%lmxl
        idmod(1:n0,is) = s_spec(is)%idmod
        do  i = 1, nel
          lphi(i,is) = -1
        enddo
        call lx2vec(lmxb,0,nphi(is),lphi(1,is))
        amass(is) = s_spec(is)%mass
        rsmfa(is) = s_spec(is)%rsmfa
      enddo
C --- Site specific
      seref = 0
      do ib = 1, nbasp
        pos(1:3,ib) = s_site(ib)%pos
        vsav(1:3,ib) = s_site(ib)%vel
        mpole(ib) = s_site(ib)%mpole
        dpole(1:3,ib) = s_site(ib)%dpole
        is = s_site(ib)%spec
        if (ib <= nbas) then
          atid(ib) = spid(is)
          ipclas(ib) = is
          seref = seref + eref(is)
        else
          atid(ib) = 'PM'
        endif
      enddo
      call numsyv(nvarms)

C --- Check if no max(lmxb)<0 for el, to reduce nel ---
      do i = nel, 1, -1
        j = i
        do is = 1, nspec
          if (lphi(i,is) >= 0) goto 14
        enddo
      enddo
   14 continue
      if (nel /= j) then
        if (iprint() >= 30) print 345, nel, j
  345   format(/' mcio: reducing nel from',i3,' to',i3)
        do  is = 1, nspec
          nphi(is) = j
        enddo
        nel = j
      endif

C --- Check whether velocites have been read from ctrl
      if (dasum(3*nbasp,vsav,1) > 1d-9) lrs(5) = 1

C --- Check if any of ebas equal, or equal to elink ---
      do  15  i = 1, nel
      do  16  j = i+1, nel
        if (dabs(el(i)-el(j)) < .2d0) call awrit2(' mcio (warning)'//
     .    '  el(%i) and el(%i) differ by less than 0.2eV',' ',80,
     .    i1mach(2),i,j)
        if (dabs(el(i)-el(j)) < 1d-3) call rx('mcio: el too close')
   16 continue
   15 continue

C --- Map LSWTCH and other vars ---
C#ifndefC LINUX_PGI | INTEL_IFORT
C      j = setsp(nsp-1)
C#else
      lsp = nsp-1
C#endif

      dx = dabc(1)
      dy = dabc(2)
      dz = dabc(3)
      nft1 = nabc(1)
      nft2 = nabc(2)
      nft3 = nabc(3)
      if (.not. tpan) lsc = 0
C#ifndefC LINUX_PGI | INTEL_IFORT
C      j = setrel(lswtch(2))
C#else
      lrel = lswtch(2)
C#endif
      call setcc(lrel)
      call setigv(1,iprint())

C ... Printout ...
      if (iprint() >= 20) then
        ls = len(strn)
        strn = ' '
        call awrit3(' st2ary:  nbas %i  nsp %i  verb %i',strn,ls,0,nbas,
     .    nsp,iprint())
        if (lswtch(2) == 0) call awrit0(strn//'%a  non-rel',strn,ls,0)
        if (lswtch(2) == 1) call awrit0(strn//'%a  scal-rel',strn,ls,0)
        if (lswtch(6) == 1) call awrit0(strn//'%a  linked-basis',
     .    strn,ls,0)
        if (lswtch(3) == 1) call awrit0(strn//'%a  frozen core',strn,
     .    ls,0)
        if (lsc == 1) call awrit0(strn//'%a  two-panel',strn,ls,0)
        call awrit0('%a',strn,-ls,-i1mach(2))
        strn = ' '
        call prxcf(s1,s2,0)
        call awrit0(' xc '//s1,strn,ls,0)
        if (s2 /= ' ') call awrit0(strn//'%a  gc '//s2,strn,ls,0)
        call awrit0('%a',strn,-ls,-i1mach(2))
        if (lxcg() == 1 .and. lxcf() /= 2)
     .    call rx('st2ary: incompatible gradient and local functionals')
      endif
      end
