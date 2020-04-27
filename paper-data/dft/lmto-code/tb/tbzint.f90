      subroutine tbzint(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_pot,s_str,s_spec,s_site,s_tb,s_strn,ltbe)
!C- TBE self-consistency loop
!C ----------------------------------------------------------------------
!Cio Structures
!Cio  s_bz   :struct for the Brillouin Zone; see structures.h
!Ci     Elts read:  lmull zval nkp efmax nkabc ntet lmet ndos dosw n
!Ci                 range w nevmx
!Co     Stored:     zval
!Co     Allocated:  *
!Cio    Elts passed:qp wtkp idtet
!Cio    Passed to:  getzv bndtb
!Cio  s_ctrl :struct for program flow parameters; see structures.h
!Ci     Elts read:  nbas nclass nl nspec nspin maxit lasa lncol lpgf
!Ci                 sdmod lstonr tol nitmv mdprm ltb zbak defm lfrce ldos
!Ci                 npadl npadr nclasp ipc
!Co     Stored:     rmax
!Co     Allocated:  rmax
!Cio    Elts passed:ldos quit lstonr nvario lasa lqp ics ipc dclabl nrc
!Cio                ips rmax initc
!Cio    Passed to:  subasi getzv rlxstp tbham bndtb secmtb tbtote relax
!Cio                tbfitmrq getfitptb mrqmintb bndtbf tbham bndtb tbtote
!Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
!Ci     Elts read:  *
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  subasi
!Cio  s_lat  :struct containing lattice information; see structures.h
!Ci     Elts read:  alat plat avw vol nsgrp awald nkd nkq gam as tol
!Ci                 rpad nkdmx nkqmx
!Co     Stored:     alat plat vol
!Co     Allocated:  *
!Cio    Elts passed:dlv qlv pos symgr istab
!Cio    Passed to:  tbham bndtb tbtote tbfitmrq mrqmintb bndtbf
!C
!Cio  s_mix :contains parameters for charge mixing; see structures.h
!Ci     Elts read:  *
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  *
!Cio  s_pot  :struct for information about the potential; see structures.h
!Ci     Elts read:  *
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:qnu pnu qc qt
!Cio    Passed to:  *
!Cio  s_str  :struct for parameters for screened strux; see structures.h
!Ci     Elts read:  mxnbr rmax
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  *
!Cio  s_spec :struct for species-specific data; see structures.h
!Ci     Elts read:  qpol vso rham iq1 iq2 stni eref z lmxa lmxl idxdn
!Ci                 mass idu uh jh coreh coreq lmxb ncomp rmt radxbs
!Ci                 colxbs name
!Co     Stored:     idxdn
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  getujh subasi getq gtpcor makidx tbxbs bndtb relax
!Cio  s_site :struct for site-specific data; see structures.h
!Ci     Elts read:  delta clabel relax spec
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  iopos rlxstp tbesel tbtote relax
!Cio  s_tb
!Ci     Elts read:  fmode nbfit ebfit alam alsc
!Co     Stored:     *
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  tbfitmrq mrqmintb
!Cio  s_strn :struct for global strings; see structures.h
!Ci     Elts read:  *
!Co     Stored:     strn
!Co     Allocated:  *
!Cio    Elts passed:*
!Cio    Passed to:  str_pack
!Ci Inputs
!Ci   prgnam:name of main program
!Ci   ltbe        :1 usual band pass
!Ci               :2 generate energy bands
!Ci               :4 band fitting mode
!Co Outputs
!Cs Command-line switches
!Cs   --DC        :
!Cs   --LUMO      :
!Cs   --RRR       :
!Cs   --X         :
!Cs   --allvecs   :
!Cs   --band      : Tabulate energy bands; see doc/Command-line-options.html
!Cs   --bangs     :
!Cs   --d=        :
!Cs   --diagn     :
!Cs   --eispack   :
!Cs   --flush     :
!Cs   --fmax=     :
!Cs   --grfac=    :
!Cs   --invbl     : Not documented
!Cs   --lc        :
!Cs   --md=       :
!Cs   --mlog      : (MPI) write MPI commands to log file
!Cs   --mv        :
!Cs   --mv=       :
!Cs   --mxq       :
!Cs   --mxq2      :
!Cs   --mxr       :
!Cs   --nochkme   :
!Cs   --nodefer   :
!Cs   --nodr      :
!Cs   --nomixkill :
!Cs   --point     :
!Cs   --rdham     :
!Cs   --sfly      :
!Cs   --st        :
!Cs   --wforce=   :
!Cs   --wpos=     :
!Cs   --wvecs     :
!Cs   --xbs       :
!Cs   --xtoll=    :
!Cs   --xyz       :
!Cs   --xyz=      :
!Cs   -akm        :
!Cs   -bmax=      :
!Cs   -cont       :
!Cs   -dumph      :
!Cs   -ef=        : Overwrite Fermi level; use with --band
!Cs   -efield=    :
!Cs   -excite=    :
!Cl Local variables
!Cl   switches in ltb:
!Cl   bit  decimal  token       switch
!Cl    0      1      OVLP       non orthogonal TB
!Cl    1      2      CRYSF      crystal field (empirical)
!Cl    2      4      OVCF        - ditto - with overlap
!Cl    3      8      ADDES      add e*S to H
!Cl    4     16      FORCES     calculate force
!Cl    5     32      FIJ        force on atom i due to atom j
!Cl    6     64      DONLY      dos only
!Cl    7    128      3PV        calculate pressure
!Cl    8    256      EVDISC     keep e'vecs on disc for BZ integration
!Cl    9    512      PAIR       pair potential only
!Cl   10   1024      TRH        calculate Tr[rho][H]
!Cl   11   2048      RHO        calculate local charges
!Cl   12   4096      SCALE      scale cut offs and rmaxh with alat
!Cl   13   8192      TBU        TB+U
!Cl   14  16384      NOUAVG     Old TB-L don't average U
!Cl   15  32768      UL         New TB-L using delta-n and average U (PRL)
!Cl   16  65536      IODEL      Read delta-H or delta-rho from disc
!Cl   17 131072      GAMMA      gamma-point only
!Cl   18 262144      MOL        do a molecule, or cluster
!Cr Remarks
!Cr   Aug 06 (ATP) TB+U: nsp is already used for S-O (JK) so
!Cr                use nspc for this as these are coupled spins;
!Cr                then nsp=2 will denote either TB+U or TB-L
!Cr                with spin polarisation. We will also have
!Cr                nsp1=2 if nsp==2 | nspc==2 for dimensioning
!Cu Updates
!Cu   12 Nov 12 migrated structures to f90 pointers
!Cu   11 Apr 11 (SL)  Species-dependent ME and more cutoff options
!Cu   10 Oct 09 (DP)  Colour weights for bands
!Cu   04 Jun 08 (ATP) Handles molecules as well as solids
!Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
!Cu   15 Feb 02 (ATP) Added MPI parallelization
!C ----------------------------------------------------------------------
      use tbprl
      use mpi
      use driver

!       use tbelstat
      use structures
      use mod_ctx

      implicit none
!C ... Passed parameters:
      character*(*) prgnam*8
      integer ltbe
!C ... For structures

      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_mix)::   s_mix
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_tb)::    s_tb
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_strn) :: s_strn(*)
!C ... Local variables
      character(len=128) :: outs
      character(len=1) :: ch
      character(len=20), parameter :: zvfrom(3) = &
     & ['ZVAL token          ', &
     &  'START category      ', &
     &  'Atomic number       ']
      logical bittst,cmdopt,lstart,lscale
      integer parg,ip,ii,jj,ldip
      logical, parameter :: t = .true., f = .false.
!       integer  oz
      integer fopn,i,j,i1,i2,ifi,ipr,iprint,i1mach,lasa,              &
     &  lpgf,lgunit,lncol,nkd,nkq,lstnr(3),maxit,nbas,nclasp,   &
     &  nclspp,nclass,nsgrp,nl,nspin,nlspc,nspc,nsp,isp,nsp1,nspec,   &
     &  sdmod,fopna,fopnn,inxsh(16),izvfr,nkp,nwin,os2
      integer, parameter :: niax = 10
!C     MPI
      logical mlog
!C ... For the Gaunt coefficients
!c lmax = 6
!c     integer indxcg(1300),jcg(6500)
!c     double precision cg(6500)
!c lmax = 10, see subs/setcg.f; Need 2*lmax(Q_L) + 1 (for forces) = 9
      integer indxcg(7400),jcg(62200)
      double precision cg(62200)
      double precision cy(289),gaunt(9,9,25),ak(9,9,9,9,3)
      double precision d1mach
!C ... Total energy and potential.  For ASA Harris functional:
      double precision thrpv,amgm,efermi(2),emad,zval
!C ... Spin dynamics
!C     double precision sdprm(6)
!C ... for PGF
      integer npadl,npadr,nbasp
!       integer opgfsl

!C     integer nzp,npl
!C     integer ozp,owz,ogll,olgii,oogll,opgplp
!C     double precision semsh(10)
!C ... tight-binding specific.  Equivalence order needed for isave
      integer ltb,mxnbr,pdim
      integer, parameter :: memx=500, NULLI = -99999
      double precision alat,rmaxh,emg,energy(7),sumev,entrpy,etot,efree, &
     &                 alpha,erep,vol,avw,ecorr,ppdip(3),fmax,tpvq
      equivalence (amgm,energy(1)), (sumev,energy(2)), (erep,energy(3)), &
     &            (etot,energy(4)), (emad,energy(5)),                    &
     &            (thrpv,energy(6)),(entrpy,energy(7))
!C ... for relaxation or dynamics
      logical md,defer,ipi
      logical :: xyzfrz(3)
      integer mdsw,nf,icom,nitrlx,natrlx,nvar,itrlx0,itrlx
      character*3 ensbl
      double precision mdprm(7),tstep,temp,pext,accrlx,taup,taub,ekin,tkin, &
               & cons,time0,time,mdtime,zeta,logs,eps,veps,zacc,zsqacc,alat0,dist0
      equivalence (tstep,mdprm(2)), (temp,mdprm(3)), (taup,mdprm(4)), (mdtime,mdprm(5)), (taub,mdprm(6)), (pext,mdprm(7))
!C ... For strux
      integer nsites,mxcsiz,nttab
      integer lmxst,nkdmx,nkqmx,nlmq,nlmq1
      integer(8) :: memstr
      double precision ekap(20),plat(3,3),qlat(3,3),as,ewtol,rpad,gam(4),gx,gy,gz,gt
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
!C ... For O-N
      integer mordrn
      logical lls
!C ... for charge mixing
      logical parmxp,uconv,lc
      double precision beta,betav(2),tjmax,ddot,cnvg,cnvm,rmsdel,rms2,dum(3),dum2,wt(3),wc
      integer nelts,neltsm,neltst,mmix,nitmax,it,it0,nkill,broy,magmix
!C ... for bzmaps,dos
      logical bzmp
      integer :: npln,nwmx, nqmx
      integer mull !,onq,onw,ovx,ovy,oxx,oyy,ozz,ozos
      double precision wtbzm,ckbas,cksumf
!C ... for TBBND
      integer ldim,idim,nlb,nlbmx,ifib,nq,iq,nevmx,nev,nfbn(2),fopno
      parameter (nlbmx=2000)
      double precision q(3),xx,efmax
      integer onesp 
      character*120 strn,plbopt
!C ... for the tight-binding hamiltonian
      logical rl, mixrho, mixh, mixQ, mixQ2
      logical ovlp,cryf,ocryf,addsll,force,pv,pair,rdham,point,          &
     &  mol,ul,tbu,sc,lso,tbe,tbbnd,tbfit,iodel,ltmp,iostr,fitpar,hdiffs,trh
      integer n,m,nterm,nset,npair,nlme,nlmesp,str_pack !,memode
      integer ntermx,memodx
      parameter (ntermx=9, memodx=7)
      integer nfilin
!       odrosl, ov0,otabme,odecay,oitab,oidec,otabcf,otabov,otbocf,odeccf,oew
!       odecov,odcocf,oitbcf,oitbov,oitocf,oidcf,oidov,oidocf,ocut,ocutov,
!       onpm,oiam,ohso,odh,odhcf,odov,odovcf,,oamassorhon,oqp,onpr,oifit,otmp,oemag,


      !i-PI helper vectors
      integer natoms
      real(8) stress(3,3), cell(3,3)
      real(8), allocatable :: atoms(:,:), forces(:,:)

      logical :: lroot, longrange
      type(tbc_t) :: tbc

      integer :: err, pid, nproc, kblocks, d2nproc, nprow, npcol, xl
      integer, allocatable, dimension(:,:) :: idxdn,iq1,iq2,npr

      integer, allocatable, dimension(:) :: ics, iprmb, idu, itab, idec, &
     & memod, itbcf, itbov, itocf, idcf, idov, idocf, ppmod, poly,cutmd, &
     & npm, iam, itmp, lmx, lmxl, pgfsl, ntab, ifit, ifbls, iax, indrx,  &
     & molidc

      real(8), allocatable, dimension(:,:) :: delta,qlv2,dlv2,qpol,vso,hdelta

      real(8), allocatable, dimension(:) :: stni, eref, esite, fnou,    &
     & fnom, fr, fe, rho, mmom,  wk, dq, qmpol, delh, delL, uh, jh,     &
     & rhoc, eband, tabme, decay, decov, v0, tabcf, tabov, tbocf, deccf,&
     & dcocf, cut, cutov, cutpp, rtab, rham, pot0, h0, ov, dh, dhcf, &
     & dov, dovcf, amass, dros, rhon, qp, vel, rho0, rhoit, mmom0, rhl0,&
     & rhlm0, awk, awk1, awk2, tj, tj1, tj2, rhlit, qit, mit, e, z,emag,&
     & p, w, zos, evl

      real(8), allocatable, target :: h(:)

      real(8), allocatable :: strx(:,:,:,:), dstrx(:,:,:,:)
!       type(parray_r8), allocatable ::  strxp(:,:), dstrxp(:,:)
      integer :: struxsize, dstrxsize
      integer, allocatable :: struxidx(:,:), dstrxidx(:,:)
      real(8), allocatable :: struxd(:), dstrxd(:)
      logical, parameter :: cmpstrx = .false.
      type(str_str0) :: strh
      interface
         subroutine safe_alloc_i(a,s,z)
            implicit none
            integer, intent(inout), allocatable :: a(:)
            integer, intent(in) :: s
            logical, intent(in) :: z
         end subroutine safe_alloc_i

         subroutine safe_alloc_r8(a,s,z)
            implicit none
            real(8), intent(inout), allocatable :: a(:)
            integer, intent(in) :: s
            logical, intent(in) :: z
         end subroutine safe_alloc_r8
      end interface

      complex(8), allocatable :: hso(:)

!     atom parallel index, see tbprl.F90 for details
      integer, allocatable :: pidx(:)

      integer, parameter :: ntm(0:memodx) = [1,1,1,1,9,5,2,4]
!c        memode 0 1 2 3 4 5 6 7 - - -

!       data oew /1/ onq /1/ onw /1/ ovx /1/ ovy /1/ oxx /1/ oyy /1/ ozz /1/ ozos /1/ oindex /1/ oemag /1/
!      ocut /1/ ocutov /1/
!      oitab /1/ oidec /1/ odeccf /1/ otabcf /1/ otabov /1/
!      &  otbocf /1/ odcocf /1/ oitbov /1/ oitocf /1/
!      &  odecov /1/ oitbcf /1/ oidcf /1/ oidocf /1/ oidov /1/

!       data xyzfrz /F,F,F/ npln /0/ nwmx /0/ nqmx /0/
!C     data autime /0.048377d0/
!       data zvfrom /'ZVAL token','START category','Atomic number'/


      call tcn('tbzint')

      call tcn('brlxlp')

      ipi = .false.
      ifi = 0
      npln = 0
      nwmx = 0
      nqmx = 0
      xyzfrz = [.false.,.false.,.false.]

!       Preliminary for io purposes. Properly done after nkp is known
      call mpi_comm_size(mpi_comm_world, nproc, err)
      call mpi_comm_rank(mpi_comm_world, pid  , err)
      lroot = pid == 0

      fitpar = .false.

      mlog = cmdopt('--mlog',6,0,outs)

      tbe   = ltbe == 1
      tbbnd = ltbe == 2
      tbfit = ltbe == 4

      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nspin = s_ctrl%nspin
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = s_ctrl%nclasp
      maxit = s_ctrl%maxit
      lasa = s_ctrl%lasa
      lncol = s_ctrl%lncol
      lpgf = s_ctrl%lpgf(1)
      sdmod = s_ctrl%sdmod
!C     sdprm(1:5) = s_ctrl%sdprm
      lstnr = s_ctrl%lstonr
      cnvg = s_ctrl%tol(1)
      cnvm = s_ctrl%tol(3)
      mull = s_bz%lmull
      allocate(ics(nclasp))
      call icopy(nclasp,s_ctrl%ics,1,ics,1)

      natoms = nbas
      allocate (atoms(3,natoms), forces(3,natoms))
!C     allocate(idmod(nl,nclasp))
!C     call spec2class(s_spec,nclasp,s_ctrl%ics,'idmod',nl,idmod,dum)

      allocate(stni(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'stni',1,stni,stni)

      alat = s_lat%alat
      plat = s_lat%plat
      avw = s_lat%avw
      vol = s_lat%vol
      nsgrp = s_lat%nsgrp
      nitrlx = s_ctrl%nitmv
      mdprm = s_ctrl%mdprm
      ltb = s_ctrl%ltb
      mxnbr = s_str%mxnbr
      rmaxh = s_str%rmax
!       awald = s_lat%awald
      nkd = s_lat%nkd
      nkq = s_lat%nkq
      gam = s_lat%gam
      as = s_lat%as
      ewtol = s_lat%tol
      rpad = s_lat%rpad
      nkdmx = s_lat%nkdmx
      nkqmx = s_lat%nkqmx
      call dpzero(betav,2)


!C ... Pre-initialize some variables that might not be initialised but
!C     referenced (depending on switches)
      nvar = 0
!C ... same for w-pointers
!       opot0 = 1
!       oov = 1
!       odelL = 1
!       ojh = 1

!C --- make larger arrays for lattice vectors for the case of NPT ---
      allocate(dlv2(3,nkdmx), qlv2(3,nkqmx))
      call dcopy(3*nkd,s_lat%dlv,1,dlv2,1)
      call dcopy(3*nkq,s_lat%qlv,1,qlv2,1)


!C ... switches
      ovlp   = bittst(ltb,1)
      cryf   = bittst(ltb,2)
      ocryf  = bittst(ltb,4)
      addsll = bittst(ltb,8) .and. ovlp
      force  = bittst(ltb,16)
!C     stress = bittst(ltb,32)
      pv     = bittst(ltb,128)
      pair   = bittst(ltb,512)
      point  = bittst(ltb,2**9) .or. cmdopt('--point',7,0,strn)
!c     if (bittst(ltb,2**9)) then
!c       point = .true.
!c     else
!c       point = cmdopt('--point',7,0,strn)
!c     endif
      trh    = bittst(ltb,2**10)
      lscale = bittst(ltb,2**12)
      TBU    = bittst(ltb,2**13)
      UL     = bittst(ltb,2**15)
      MOL    = bittst(ltb,2**18)
      iodel  = bittst(ltb,2**16)
      sc = UL .or. TBU
      lso    = bittst(lncol,4)
!C     bzmp = lgors('ctrl ldos,8',sctrl)
      bzmp = iand(s_ctrl%ldos,8)/=0
!C ... dipole correction to Ewald
!C...  under development. Hardwire ldip=3 for now
!c     if (cmdopt('--sphere',8,0,strn)) then
        ldip = 3
        call info0(30,1,0,' TBZINT: use spherical dipole correction')
!c       call info0(30,0,0,' Warning: spherical dipole correction is'//
!c    .     ' incomplete')
!c     elseif (cmdopt('--slab',6,0,strn)) then
!c       ldip = 2
!c       call info0(30,1,0,' TBZINT: use slab dipole correction')
!c       call rx0(' tbzint: slab dipole correction is not implemented')
!c     else
!c       ldip = 0
!c     endif
!c     if(ldip /= 0 .and. MOL) call
!c    .  rx0(' tbzint: dipole correction only applies to periodic BC')

!c ... scale or do not scale cutoff distances and RMAXH with alat
!C     (passed to tbham, tbtote, shotbm and used in the construction
!C      of the neighbour table)

!C     Pre-initialize some variables whose initialisation depends on switches
!C --- defer force calculation until self consistent ---
      defer = T
      if (cmdopt('--nodefer',9,0,outs)) defer = F

!C --- read hamiltonian and overlap as strux from disc ---
      rdham = cmdopt('--rdham',7,0,outs)
      if (rdham) sc = F

      lls = cmdopt('--linscl',8,0,outs)

!C --- retain nspin for dimensioning certain arrays, eg qnu ---
      nsp1 = nspin
      if (TBU .or. UL) then
        if (nspin == 1 .and. TBU) call rx('TBZINT: for TB+U restart with NSPIN=2')
        nsp = nspin
        nspc = 1
      else
        nspc = nspin
        nsp = 1
      endif
!C ... Modes of mixing
!C 1. Mix density matrix
      mixrho = F
!C 2. Mix multipole moments and magnetic moments separately: beta, betav
      mixQ2 = cmdopt('--mxq2',6,0,outs)
!C 3. Mix multipole moments and magnetic moments togther: beta
      mixQ = cmdopt('--mxq',5,0,outs) .or. mixQ2
!C 4. Mix hamiltonian (default); cannot use with overlap
      mixh = T
      if (cmdopt('--mxr',5,0,outs)) then
        mixrho = T
        mixh = F
      endif
      lc = F
      if (pair .or. ovlp .or. mixQ .or. mixQ2) then
        call rxx(bittst(ltb,2**14),' TBZINT: need ave. U for Q mixing')
        lc = cmdopt('--lc',4,0,outs)
        mixQ = T
        mixh = F
      endif

      hdiffs = force .or. pv
      ckbas = cksumf(s_lat%pos,3*nbas)
      mdsw = nint(mdprm(1))
      md = mdsw > 0 .and. mdsw < 5
      zval = s_bz%zval
      izvfr = 1
      nlspc = nl*nspc*nclass
      nbasp = nbas + npadl + npadr
      nclspp = 2*nclasp-nclass
!C ... Hardwire no order-N  for now
      mordrn = 0
!C     call defrr(ovdif,-nclspp)
!C ... eref, lmx,z by class
!C     call redfrr(oeref,2*nclasp)
      allocate(eref(2*nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'eref',1,[0],eref)
      call dvset(eref,1+nclasp,2*nclasp,0d0)
!       call defrr(oz,nclasp)
      allocate(z(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'z',1,[0],z)
      allocate(lmx(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxa',1,lmx,[0])
      allocate(lmxl(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxl',1,lmxl,[0])
      allocate(idxdn(nl,nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'idxdn',nl,idxdn,[0])


      ! init i-pi
      if (cmdopt('--ipi-host=',11,0,outs)) then
           ipi = .true.
           md = .true.
           call driver_init(lroot)
      endif

!C ... Spin-orbit coupling parameters stashed in alpha
!C     call defdr(odel,-nl*nsp1*nbasp)
      allocate(delta(nl*nsp1,nbasp)); delta = 0.0_8
!C     call spackv(10,'site delta',ssite,1,nbasp,w(odel))
      do  i = 1, nbasp
        delta(1:nl*nsp1,i) = s_site(i)%delta(1:nl*nsp1)
      enddo

!C ... TB polarisation parameters, potentials and Hubbard U's and J's
      if (sc) then

!C ...   Adjust lmxl and set nlmq, nlmq1
!C       call spec2class(s_spec,nclasp,s_ctrl%ics,0,10,dum,w(oqpol))
        allocate(qpol(10,nclasp))
        do  i = 1, nclasp
          qpol(:,i) = s_spec(ics(i))%qpol
        enddo

        call sulmxl(force,pv,point,nl,nclasp,qpol,lmxl,nlmq,nlmq1)

!         call defdr(oqmpol,nlmq*nbas)
!         call defdr(odelh,-nl**4*nbas*nsp)
!         call defdr(odelL,-nl**4*nbas*nsp)
        allocate(qmpol(nlmq*nbas), delh(nl**4*nbas*nsp), delL(nl**4*nbas*nsp), uh(4*nbas), jh(4*nbas) )
        allocate(idu(4*nbas))
        delh = 0.0_8; delL = 0.0_8; uh = 0.0_8; jh = 0.0_8; idu = 0
        if (ovlp) then
          allocate(pot0(nbas)); pot0 = 0.0_8
!         else
!           opot0 = 1
        endif

!         call defdr(ouh,-4*nbas)
!         call defdr(ojh,-4*nbas)
!         call defi(oidu,-4*nbas)
        if (nsp == 2) then
          call getujh(s_spec,nl,nbas,s_ctrl%ipc,s_ctrl%dclabl,idu,uh,jh)
        endif
      else
!C ... Default dimensions for multipole arrays
        nlmq = 9
        nlmq1 = 16
      endif

!C ... tbe-specific initialization
!C ... for alpha, see Phys Rev B, 53, 15381 (1996)
      alpha = 0d0
      entrpy = 0d0
      tpvq = 0
      if (bzmp) then
        call rxx(bzmp,'TBZINT: not ready for bzmaps')
        call rx('bzmp branch not checked')
        ifi = fopno('BZPL')
        rewind ifi
        nkp = 0
   93   read (ifi,*,err=96,end=96) nwin,npln
        do i = 1, nwin
          read (ifi,*,err=96,end=96) xx,xx
        end do
        do i = 1, npln
          read (ifi,*,err=96,end=96) xx,xx,xx,xx,xx,n,xx,xx,xx,xx,xx,m,xx
          nkp = nkp + n*m
        end do
        goto 93
   96   call rxx(nkp == 0,prgnam//': empty or badly formed BZPL file')
      else
        nkp = s_bz%nkp
      endif

!       call init_ctx(tbc % c3d, 3, dims_order = [2,0,1], root = 0, comm = mpi_comm_world)
! Figure out dimensions for the 3d process mesh such that priority is given to the k-point parallelization since it is more efficient.
      if (.not. cmdopt('--kblocks=',10,0,outs)) then
         kblocks = min(nkp,nproc)
         do while (mod(nproc, kblocks) /= 0)
            kblocks = kblocks - 1
         end do
      else
         read(outs(11:len(outs)), *) kblocks
      end if
      d2nproc = nproc/kblocks

      if (.not. cmdopt('--nprows=',9,0,outs)) then
         nprow = nint(sqrt(real(d2nproc)))
         do while (nprow > 1 .and. (d2nproc/nprow)*nprow /= d2nproc)
            nprow = nprow - 1
         end do
      else
         read(outs(10:len(outs)), *) nprow
      end if

      npcol = d2nproc/nprow

      call init_ctx(tbc % c3d, 3, sizes = [kblocks,nprow,npcol], reorder = .false., root = 0, comm = mpi_comm_world)

      tbc % c1d = sub_cart_ctx(tbc % c3d, [.true. , .false., .false.])
      tbc % c2d = sub_cart_ctx(tbc % c3d, [.false., .true. , .true. ])


      tbc % desc(2) = tbc % c2d % comm
      call blacs_gridinit(tbc % desc(2), 'r', nprow, npcol)
      if (nprow*npcol < 1 .and. tbc % c2d % sz > 1) then
         if (nkp < 2*nproc) then
            call rx('TBZINT: nproc > nkp, consider enabling SCALAPACK or launch fewer processes')
         end if
      end if

      tbc % sl = nprow*npcol > 1
      if (tbc%sl) then
         call blacs_gridinfo(tbc % desc(2), nprow, npcol, i, j)
         if (      (nprow /= tbc%c2d%szs(0)) .or. (npcol /= tbc%c2d%szs(1)) &
            & .or. (    i /= tbc%c2d%crd(0)) .or. (    j /= tbc%c2d%crd(1))) then
            print *, 'BLACS: size, coords: ', nprow, npcol, i, j
            print *, 'MPICR: size, coords: ', tbc%c2d%szs(0:1), tbc%c2d%crd(0:1)
            print *, "TBZINT: Mismatch between MPI 2D Cartesian topology and BLACS's context"
            call rx('Likely an inconvenient number of processes was launched, choose lesser or equal to nkp multiple of squares')
         end if
      end if

      lroot  = tbc % c3d % lrt
      pid    = tbc % c3d % id
      nproc  = tbc % c3d % sz

      allocate(tbc % kmap(0:kblocks), pidx(nkp))
      pidx = [(i, i = 1,  nkp)]
      call vbdist(nkp, pidx, kblocks, tbc % kmap)
      deallocate(pidx)

      flush(6); call mpi_barrier(mpi_comm_world, err)
      if (lroot) print '(x,a,3(x,i0))', 'k-point : process map : ', tbc % kmap

      allocate(tbc % kcount(0:kblocks-1))
      tbc % kcount = tbc % kmap(1:kblocks) - tbc % kmap(0:kblocks-1)

      allocate(tbc % d2amap(0 : tbc % c2d % sz), tbc % d2count(0 : tbc % c2d % sz - 1))

      if (.not. (tbe .or. tbbnd)) goto 499

!       nband = nbas*nl**2
!       call defdr(oeband,nband*nsp1*nkp)
!       allocate(eband(nband*nsp1*nkp)) ! now it is only passed to tbfitmrq and allocated right before it
      call subasi(s_ctrl,s_spec,s_ham)

!C This is now done in v7input/rdctrl2
!CC ... Read positions from file
!C      if (cmdopt('--rpos=',7,0,outs)) then
!C        if (lroot) then
!C          call iopos(F,0,outs(8:),nbasp,s_lat%pos,s_site)
!C        endif
!C        call mpibc1(s_lat%pos,nbasp,4,mlog,'tbzint','pos')
!C      endif

!C ... Write positions to file, nit=0
      if (nitrlx == 0) then
        if (cmdopt('--wpos=',7,0,outs) .and. lroot) call iopos(T,1,outs(8:),nbasp,s_lat%pos,s_site)
        call rx(prgnam//': no iterations')
      endif

!C ... Make only DOS weights; not total energy, E_F, etc.
!C     if (lgors('ctrl ldos,64',sctrl)) goto 499
      if (iand(s_ctrl%ldos,64)/=0) goto 499

!C --- Get zval. Priority is: 1. ZVAL from CTRL; 2. from moms; 3. from Z
      if (zval == 0 .or. zval == NULLI) then
        izvfr = 2
        call getmom(nsp1,nl,s_pot%qnu,nclass,s_ctrl%nrc,idxdn,zval)
        if (zval == 0) then
          izvfr = 3
          call pshpr(0)
          allocate(dq(nclasp))
          call getq(nsp1,nl,lmx,nclasp,z,s_pot%pnu,s_pot%qnu,s_ctrl%ics,s_spec,s_pot%qc,s_pot%qt,dq)
          call getzv(nclasp,s_ctrl%nrc,z,s_pot%qc,s_bz,s_ctrl,zval)
          deallocate(dq)
          call poppr
        endif
      endif
      if (iprint() >= 30) then
        call awrit1('%N TBZINT: zval=%d electrons from '//zvfrom(izvfr),' ',128,i1mach(2),zval)
      endif

!C --- Set up relaxation parameters ---
      allocate(itmp(6*nbas))
      call rlxstp(s_ctrl,s_site,natrlx,nvar,itmp,xyzfrz,pdim)
      icom = 0
      if (nvar /= 0 ) then
        allocate(indrx(2*natrlx))
        indrx(1:2*natrlx) = itmp(1:2*natrlx)
!         call defdr(ow,nvar*nvar)
!         call defdr(op,pdim)
        allocate(w(nvar*nvar))
        allocate(p(pdim))
      endif
      deallocate(itmp)

!C --- Set up force/charge symmetrization ---
  499 continue

!C ... A few sanity checks for tbfit
      if (tbfit) then
!c       if (memode >= 10)
!c    .    call rx(prgnam//': canonical ham not allowed for fit')
        call rxx(mordrn /= 0,prgnam//': order N not implemented for fit')
        call rxx(lpgf /= 0,prgnam//': PGF not implemented for fit')
        call rxx(UL,prgnam//': UL=T not implemented for fit')
        call rxx(nitrlx>1,prgnam//': cannot do relaxation with fit')
      endif

!C --- Read tight-binding matrix elements ---
      if (lroot) nfilin = fopna('CTRL',-1,1)

!C     An upper bound to the number of coefficients to be read in
      nterm = ntermx
!C     Note: ME are spin dependent in noncollinear case:
!C     Then groups of three blocks (++, --, +-)
      nlmesp = nlme(nl)*(2*nspc-1)
      nset = memx
!       call defdr(otabme,-nterm*nlmesp*3*nset)
      allocate(tabme(nterm*nlmesp*3*nset)); tabme = 0.0_8
      allocate(memod(nset)); memod = 0
      allocate(v0(9)); v0 = 0.0_8
      allocate(decay(nlmesp*3*nclass**2)); decay = 0.0_8
!C --- First pass to get dimensions ---
!      call rdtbh(nfilin,nclass,s_ctrl%dclabl,z,nlmesp,ltb,nterm,1,F,rmaxh,0,0, &
!       & tabme,tabcf,tabov,tbocf,decay,deccf,decov,dcocf,cut,cutov,itab,            &
!       & itbcf,itbov,itocf,idec,idcf,idov,idocf,v0,nset,memod,ppmod,poly,cutmd,cutpp)
      call rdtbh(nfilin,nclass,s_ctrl%dclabl,z,nlmesp,ltb,nterm,1,F,rmaxh,0,0, &
       & tabme,v0,v0,v0,decay,v0,v0,v0,v0,v0,v0,            &
       & v0,v0,v0,v0,v0,v0,v0,v0,nset,memod,v0,v0,v0,v0)
!C      call rdtbh(nfilin,nclass,s_ctrl%dclabl,z,nlmesp,ltb,nterm,1,F,0,
!C     &  0,tabme,tabcf,tabov,tbocf,decay,deccf,
!C     &  decov,dcocf,cut,cutov,itab,itbcf,
!C     &  itbov,itocf,idec,idcf,idov,idocf,v0,
!C     &  nset,memode)

!C...  find the min dimension for parameter matrices
      nterm = 0
!       memode = memod
      do i = 1, nset
        nterm = max(ntm(memod(i)),nterm)
      enddo
!       print *, nterm
!       stop
      if (nterm > ntermx) call rxi(' tbzint: nterm = %i exceeds ntermx. Increase ntermx or check the input',nterm)

      deallocate(tabme, memod, v0, decay)
      allocate(memod(nset)); memod = 0
      allocate(tabme(nterm*nlmesp*nset)); tabme = 0.0_8
      allocate(decay(nlmesp*nset)); decay = 0.0_8
      allocate(cut(2*nlmesp*nset)); cut = 0.0_8

      allocate(itab(1),idec(1),tabcf(1),deccf(1),itbcf(1),idcf(1))
      allocate(tabov(1),decov(1),cutov(1),itbov(1),idov(1),tbocf(1),dcocf(1),itocf(1),idocf(1))
      if (tbfit) then
        deallocate(itab); allocate(itab(nterm*nlmesp*nset)); itab = 0
        deallocate(idec); allocate(idec(nlmesp*nset)); idec = 0
      endif
      if (cryf) then
        deallocate(tabcf); allocate(tabcf(nterm*nlmesp*nset)); tabcf = 0.0_8
        deallocate(deccf); allocate(deccf(nlmesp*nset));       deccf = 0.0_8
        if (tbfit) then
          deallocate(itbcf); allocate(itbcf(nterm*nlmesp*nset)); itbcf = 0
          deallocate(idcf);  allocate(idcf(nlmesp*nset));        idcf  = 0
        endif
      endif
      if (ovlp) then
        deallocate(tabov); allocate(tabov(nterm*nlmesp*nset)); tabov = 0.0_8
        deallocate(decov); allocate(decov(nlmesp*nset));       decov = 0.0_8
        deallocate(cutov); allocate(cutov(2*nlmesp*nset));     cutov = 0.0_8
        if (tbfit) then
          deallocate(itbov); allocate(itbov(nterm*nlmesp*nset)); itbov = 0
          deallocate(idov);  allocate(idov(nlmesp*nset));        idov = 0
        endif
      endif
      if (ocryf) then
        deallocate(tbocf); allocate(tbocf(nterm*nlmesp*nset)); tbocf = 0.0_8
        deallocate(dcocf); allocate(dcocf(nlmesp*nset));       dcocf = 0.0_8
        if (tbfit) then
          deallocate(itocf); allocate(itocf(nterm*nlmesp*nset)); itocf = 0
          deallocate(idocf); allocate(idocf(nlmesp*nset));       idocf = 0
        endif
      endif
      allocate(V0(9*nset)); v0 = 0.0_8
      allocate(ppmod(nset)); ppmod = 0
      allocate(poly(nset));  poly = 0
      allocate(cutmd(nset)); cutmd = 0
      allocate(cutpp(nset*2)); cutpp = 0.0_8
      allocate(npm(2*nclass))
      allocate(itmp(4*3*nclass**2))
!C --- Second pass to get remaining parameters and pointer arrays ---
      npair = nset
      call rdtbh(nfilin,nclass,s_ctrl%dclabl,z,nlmesp,ltb,nterm,2, &
     &  tbfit,rmaxh,itmp,npm,tabme,tabcf,tabov,tbocf,decay,deccf,decov, &
     &  dcocf,cut,cutov,itab,itbcf,itbov,itocf,idec,idcf,idov,idocf,v0, &
     &  npair,memod,ppmod,poly,cutmd,cutpp)
!C      call rdtbh(nfilin,nclass,s_ctrl%dclabl,z,nlmesp,ltb,nterm,2,tbfit,
!C     &  iam,npm,tabme,tabcf,tabov,tbocf,
!C     &  decay,deccf,decov,dcocf,cut,cutov,
!C     &  itab,itbcf,itbov,itocf,idec,idcf,
!C     &  idov,idocf,v0,npair,memode)
!       deallocate(iam)
      allocate(iam(3*npair))
      iam = itmp(1:3*npair)
      deallocate(itmp)
      call fclose(nfilin)

      if (.not. cmdopt('--nochkme',9,0,outs)) then
!C ...   consistency check for ME and PP parameters
        call chkme(0,0,nl,nlmesp,nspc,ltb,nterm,tabme,decay,cut,v0,nset,memod,ppmod,poly,cutmd,cutpp)
        if (ovlp) call chkme(1,1,nl,nlmesp,nspc,ltb,nterm,tabov,decov,cutov,v0,nset,memod,ppmod,poly,cutmd,cutpp)
      endif

!C --- Printout ---
      if (iprint() > 30) then
        if (tbe .or. tbbnd) then
          i = 1
!c         if (memode == 6) i = 0
        else
          i = 0
        endif
        call shotbm(lscale,alat,tabme,decay,v0,nset,npair,nterm,nl,nspc, &
            & nlmesp,memod,ppmod,poly,cutmd,cutpp,cut,nclass,s_ctrl%dclabl,iam,i,'Hamiltonian')
        if (cryf) call shotbm(lscale,alat,tabcf,deccf,v0,nset,npair,nterm,nl,nspc,nlmesp, &
            & memod,ppmod,poly,cutmd,cutpp,cut,nclass,s_ctrl%dclabl,iam,0,'Hamiltonian crystal field')
        if (ovlp) call shotbm(lscale,alat,tabov,decov,v0,nset,npair,nterm,nl,nspc,nlmesp, &
            & memod,ppmod,poly,cutmd,cutpp,cutov,nclass,s_ctrl%dclabl,iam,0,'overlap')
        if (ocryf) call shotbm(lscale,alat,tbocf,dcocf,v0,nset,npair,nterm,nl,nspc,nlmesp,memod,ppmod,poly,cutmd, &
            & cutpp,cutov,nclass,s_ctrl%dclabl,iam,0,'overlap crystal field')
      endif

!C --- Set up for spin-orbit matrix elements ---
      if (lso) then
!C       call spec2class(s_spec,nclasp,s_ctrl%ics,0,nl,dum,w(ovso))
        allocate(vso(nl,nclasp))
        do  i = 1, nclasp
          vso(1:nl,i) = s_spec(ics(i))%vso(1:nl)
        enddo

!C       call shstru('spec',sspec,1,nspec)
!C       call prmx('vso',vso,nl,nl,nclasp)
        allocate(hso(4*nl**4)); hso = 0.0_8
        call skhso(nl,hso)
!       else
! C       call defrr(ovso,1)
!         call defrr(ohso,1)
      endif

!C --- Make mxnbr ---
      if (rmaxh > d1mach(3)) then
        if (mxnbr == 0) then
!c         mxnbr = 4*nint((rmaxh*alat/avw)**3*nbasp)
          if (lscale) then
            mxnbr = max(4*nint((rmaxh*alat/avw)**3*nbasp),nbas*nbas)
          else
            mxnbr = max(4*nint((rmaxh/avw)**3*nbasp),nbas*nbas)
          endif
        else
          mxnbr = mxnbr*nbas
        endif
      endif

!C --- Get LL' on-site hamiltonian melts from disk
      if (iodel .and. sc) then
        if (lroot) then
          call iodelL(ltb,1,nl,nbas,nsp,delta,delL)
        endif
        call mpibc1(delL,nl**4*nsp*nbas,4,mlog,'tbzint','delL')
      endif

      nelts = nl*nspc
      if (nspc == 2 .and. nsp == 2) call rx(' UL=T not set up for spin orbit coupling yet')

!C --- Set up for tight-binding MD ---
      if (nvar == 0 .and. .not. md) nitrlx = 1
      itrlx0 = 0
!       call defdr(ovel,-3*nbas)
      allocate(vel(3*nbas)); vel = 0.0_8
      if (md .and. .not. ipi) then
        if (lroot) ifi = fopn('STRT')
!C --- choose MD ensemble ---
        alat0 = alat
        ensbl = 'NVE'
        if (mdsw == 2) then
          ensbl = 'NVT'
        endif
        if (mdsw == 3) then
          call rxx(.not.pv,'TBZINT: set 3PV=T for NPT')
          ensbl = 'NPT'
        endif
        nitrlx = mdtime / tstep + 1
        time0 = 0d0
        eps = 0d0
!         call defdr(oamass,-nspec)
        allocate(amass(nspec)); amass = 0.0_8
!C       call spackv(10,'spec mass',sspec,1,nspec,amass)
        call spec2class(s_spec,nspec,-1,'mass',1,dum,amass)
        if (ddot(nspec,amass,1,amass,1) == 0) call rx('TBZINT: masses must be all nonzero')
!c ---   Is it a new run or continuation from the previous run? ---
!c ...   lstart = .true. if new MD run or restart file not found
        lstart = T
        if (.not. cmdopt('--st',4,0,outs)) then
!c ...   Attempt to read restart file
          call iostrt(nbas,nf,itrlx0,s_lat%pos,vel,eps,zeta,zacc,zsqacc,veps,time0,ifi,err)
          if (err == 0) then
            lstart = F
          else
            if (iprint() >= 10) print *,' TBZINT *warning* could not read strt file starting new MD'
          endif
        endif
        if (lstart) then
          call initv(ensbl,nbas,nf,tstep,s_ctrl%ips,amass,temp,pext,taup,taub, &
                     & min(nspec,nclasp),1,' ',s_ctrl%dclabl,zeta,zacc,zsqacc,veps,vel)
          call zercmv(nbas,vel,amass,s_ctrl%ips)
        else
          if (ensbl == 'NPT') then
!C ---       remake vectors in scaled alat ---
            alat = alat0 * exp(eps)
            lmxst = 6
            call pshpr(0)
            call lattc(as,ewtol,rpad,alat,alat,plat,gx,gy,gz,gt,plat,qlat, &
                        & lmxst,vol,s_lat%awald,dlv2,nkd,qlv2,nkq,nkdmx,nkqmx)
            call poppr
            s_lat%alat = alat
            s_lat%plat = plat
            s_lat%vol = vol
          endif
        endif
        time = time0
        if (ensbl == 'NVE') then
          logs = 0d0
        else
          logs = zacc * tstep
        endif
        mdtime = time0 + mdtime
!C ... end of MD set up
      endif
!C --- Make Gaunt (CG) coefficients and allocate strux ---
      if (sc) then
        call sylmnc(cy,16)
        call scg(9,cg,indxcg,jcg)
        call makcg9(indxcg,jcg,cg,gaunt,ak)
        if (.not. cmdopt('--sfly',6,0,outs)) then
!       Set process:atom-range mapping
          if (.not. allocated(tbc%esamap)) then
             allocate(tbc%esamap(0:nproc), pidx(nbas))

             xl = 1
             if (nlmq==nlmq1) xl = 2
             pidx(1) = (lmxl(s_ctrl%ipc(1))+xl)**2
             do i = 2, nbas
               pidx(i) = pidx(i-1) + (lmxl(s_ctrl%ipc(i))+xl)**2
             end do
             call vbdist(nbas, pidx, nproc, tbc%esamap)
             deallocate(pidx)
             allocate(tbc%escount(0:nproc-1))
             tbc%escount = tbc%esamap(1:nproc) - tbc%esamap(0:nproc-1)
          end if
!         transposed for more efficient, contiguous, access. same applies to dstrx
!         (L,L',R,R') -> (L',L,R',R)
!          if (.not. cmpstrx) then
!            allocate (strx(nlmq,nlmq1,nbas,tbc%esamap(pid+1)-tbc%esamap(pid)),stat=err)
!
!            if (err /= 0) call rx0(' tbzint: failed to allocate strx')
!            if (pv) then
!             allocate (dstrx(nlmq,nlmq,nbas,tbc%esamap(pid+1)-tbc%esamap(pid)),stat=err)
!              if (err/=0)call rx0(' tbzint: failed to allocate dstrx')
!            else
!              allocate (dstrx(1,1,1,1))
!            endif
! !C...      Memory estimation
!             memstr = nlmq1
!             if (pv) memstr = memstr + nlmq
!             memstr = memstr * nlmq * nbas*(tbc%esamap(pid+1)-tbc%esamap(pid))
!             memstr = (memstr * 8) / 2**20
!             if (memstr > 24) call info2(10,1,1,' TBZINT: allocating  %i MiB memory for strux',memstr,0)
!          else
            allocate(struxidx(nbas, tbc%escount(pid)))
            if (pv) then
               allocate(dstrxidx(nbas,tbc%escount(pid)))
            else
               allocate(dstrxidx(1,1))
            end if
            call mkstrxidx(tbc, pv, force, nbas, s_ctrl%ipc, lmxl, struxsize, dstrxsize, struxidx, dstrxidx)
            memstr = nbas*tbc%escount(pid)
            if (pv) memstr = memstr*2
            memstr = (memstr*4 + (struxsize+dstrxsize)*8)/2**20
            if (memstr > 24) call info2(10,1,1,' TBZINT: allocating  %i MiB memory for strux&dstrx + indices at root',memstr,0)
            allocate(struxd(struxsize))
            allocate(dstrxd(dstrxsize))
!          end if
         else
!c ... Need to allocate strux arrays anyway for passing into subroutines
          allocate (strx(1,1,1,1))
          allocate (dstrx(1,1,1,1))
        endif
      endif
!C --- set up arrays for forces, density matrix ---
!       call defdr(oesite,-nbas*nsp1)
!       call defdr(ofnou,-4*nbas)
!       call defdr(ofnom,-4*nbas)
!       call defdr(of,-3*nbas)
!       call defdr(ofe,-3*nbas)
!       call defdr(orho,-6*nbas)
!       call defdr(ommom,-nbas)
!       call defdr(orholm,-nl**2*2*nbas)

      allocate(iprmb(nbas*nl**2), esite(nbas*nsp1), fnou(4*nbas), fnom(4*nbas), &
                         & fr(3*nbas), fe(3*nbas), rho(nl*nsp1*nbas), mmom(nbas))

      esite= 0.0_8
      fnou = 0.0_8
      fnom = 0.0_8
      fr   = 0.0_8
      fe   = 0.0_8
      rho  = 0.0_8
      mmom = 0.0_8
      iprmb = 9999999

      call pshprt(0)
      call makidx(nl,1,1,nbas,0,s_spec,s_ctrl%ips,-1,iprmb,inxsh)
      call popprt
      ldim = inxsh(1)
      if (sc) then
!         call defdr(orhoc,-nl**4*nbas)
        allocate(rhoc(nl**4*nbas))
        rhoc = 0.0_8
!         odrosl = 1
        if (ovlp) then
          allocate(dros(3*nbas**2))
          dros = 0.0_8
        end if
      endif


      tbc % realsz = [ldim, ldim]

      tbc % blcksz = [32, 32]
      if (.not. tbc%sl) tbc % blcksz = [ldim, ldim]

      tbc % blocks = (tbc % realsz - 1) / tbc % blcksz + 1
      tbc % lcblks = (tbc % blocks - 1) / tbc % c2d % szs(0:1) + 1
      tbc % loclsz = tbc % lcblks * tbc % blcksz
      tbc % globsz = tbc % loclsz * tbc % c2d % szs(0:1)

      tbc % desc(1) = 1
      tbc % desc(3:4) = tbc % realsz
      tbc % desc(5:6) = tbc % blcksz
      tbc % desc(7:8) = 0
      tbc % desc(9) = tbc % loclsz(1)




!       orhon = 1
!       if (TBU) call defdr(orhon,-nl**2*nl**2*nsp*nbas)
      if (TBU) then
         allocate(rhon(nl**2*nl**2*nsp*nbas))
         rhon = 0.0_8
      endif
!C .. set up work arrays for mixing
      if (sc) then
        i1 = str_pack('mix',-2,s_strn,strn)
        tjmax = 10d0
        broy = 0
        dum(1) = 0
        dum(2) = 0
        dum2   = 0
        rmsdel = 1d0
        nkill = 0
        mmix = 0
        wt = 0

        if (lroot) call pshpr(80)

        if (.not. parmxp(1,strn,len_trim(strn),broy,mmix,wt,beta,dum2, &
     &      dum2,dum2,dum2,outs,wc,nkill,betav,rmsdel)) call rx(prgnam//': parse in parmxp failed')

        if (lroot) call poppr
!C .. input Mulliken charges (l and lm resolved) ..
!         call defdr(ommom0,-nbas)
!         call defdr(orhl0,-nl*2*nbas)
!         call defdr(orhlm0,-nl**2*2*nbas)
        allocate(mmom0(        nbas)); mmom0 = 0.0_8
        allocate(rhl0 (nl   *2*nbas)); rhl0  = 0.0_8
        allocate(rhlm0(nl**2*2*nbas)); rhlm0 = 0.0_8

        neltst = nl**4*nsp*nbas

        if (mixrho) then
          if (ovlp .or. nsp == 2) then
!             call defdr(orhlit,-nl*nsp*nbas*(mmix+2)*2)
            allocate(rhlit(nl*nsp*nbas*(mmix+2)*2)); rhlit = 0.0_8
            neltst = neltst + nl*nsp*nbas
          endif
!           call defdr(orho0,-nl**4*nsp*nbas)
!           call defdr(orhoit,-nl**4*nsp*nbas*(mmix+2)*2)
          allocate(rho0 (nl**4*nsp*nbas)); rho0 = 0.0_8
          allocate(rhoit(nl**4*nsp*nbas*(mmix+2)*2)); rhoit = 0.0_8
          call mkrho0(1,nl,nbas,nclass,s_ctrl%ipc,s_ctrl%dclabl,nsp,s_pot%qnu,rho0,rhl0,rhlm0,mmom0)
        elseif (mixh) then
!           call defdr(odelta,-nl**4*nsp*nbas*(mmix+2)*2)
          allocate(hdelta(nl**4*nsp*(mmix+2)*2,nbas)); hdelta = 0.0_8
          if (ovlp) then
            neltst = neltst + (nbas*(nbas-1))/2
          endif
        elseif (mixQ) then
          if (pair) then
!             call dcopy(nl*nsp1*nbas,0d0,0,rho,1)
            rho = 0.0_8
          else
            if (mixQ2) then
              nelts  = nlmq * nbas
              neltsm = nbas
              neltst = (nlmq + (nsp-1))*nbas
!               call defdr(oqit,-nelts*(mmix+2)*2)
              allocate(qit(nelts*(mmix+2)*2)); qit = 0.0_8
!               omit = 1
              if (nsp == 2) then
!                 call defdr(omit,-neltsm*(mmix+2)*2)
                  allocate(mit(neltsm*(mmix+2)*2)); mit = 0.0_8
              endif
!               call defdr(oawk1,-nelts*(mmix+2)*2)
!               call defdr(oawk2,-neltsm*(mmix+2)*2)
!               call defdr(otj1,mmix)
!               call defdr(otj2,mmix)
              allocate(awk1(nelts *(mmix+2)*2)); awk1 = 0.0_8
              allocate(awk2(neltsm*(mmix+2)*2)); awk2 = 0.0_8
              allocate(tj1(mmix))
              allocate(tj2(mmix))
            else
              neltst = (nlmq + (nsp-1))*nbas
!               call defdr(oqit,-neltst*(mmix+2)*2)
!               call defdr(omit,-nbas*(mmix+2)*2)
!               call defdr(oawk,-neltst*(mmix+2)*2)
!               call defdr(otj,mmix)
              allocate(qit(neltst*(mmix+2)*2)); qit = 0.0_8
              allocate(mit(nbas  *(mmix+2)*2)); mit = 0.0_8
              allocate(awk(neltst*(mmix+2)*2)); awk = 0.0_8
              allocate(tj(mmix))
            endif
            call mkrho0(0,nl,nbas,nclass,s_ctrl%ipc,s_ctrl%dclabl,nsp,s_pot%qnu,dum,rho,rhlm0,mmom0)
            if (lc) nelts = nbas
          endif
        else
          call rx(' TBZINT: No mixing option set')
        endif
        if (.not. mixQ) then
!           call defdr(oawk,-neltst*(mmix+2)*2)
!           call defdr(otj,mmix)
          allocate(awk(neltst*(mmix+2)*2)); awk = 0.0_8
          allocate(tj(mmix))
        endif
      endif


      call tcx('brlxlp')

!C ----------------- Relaxation or MD loop ---------------------
      itrlx = itrlx0+1
      do  2 while (itrlx <= itrlx0+nitrlx)

        if (ipi) then

           forces = reshape(fr,(/3,natoms/))
           stress = 0.0d0
           stress(1,1) =  thrpv/3
           stress(2,2) =  thrpv/3
           stress(3,3) =  thrpv/3
           call driver_step(natoms, atoms, cell, efree, forces, stress, lroot)
           s_lat%plat = cell/alat
           s_lat%pos = atoms/alat

           do  j = 1, natoms
              s_site(j)%pos(1:3) = atoms(1:3,j)
           enddo
          !C --- Lattice setup ---
          call setcg(s_lat,8,12)
          call lattic(s_lat,s_ctrl,s_site)
          hasdata = .true.
          itrlx = itrlx - 1  ! scales the iterator back, so as to hijack this into an infinite loop
        endif !i-PI loop!

        call tcn('bmkstrx')
        call parmx0(1,1,0d0)
        if (.not. iodel .and. sc) then
          call dcopy(nl**4*nbas*nsp,0d0,0,delL,1)
        endif

!C --- Make tight-binding Hamiltonian ---
!         call defi(onpr,2*nbas)
        allocate(npr(0:1,nbas))
!C --- Embedded cluster scheme (order N) ---
        if (mordrn == 1) then
          call rx('no order-N scheme implemented')
        endif
!C   --- Get neighbor table iax for each atom in the cluster ---
        if (lpgf /= 0) then
          i = 2
          j = 1
        else
          i = 3
          j = -1
        endif
!C   ... Make nttab,ontab,oiax
!         call defdr(opgfsl, nbas)
        allocate(pgfsl(nbas))

        if (rdham) then
          ckbas = cksumf(s_lat%pos,3*nbas)
          call pshpr(100)
          call free_str0(strh)
          ltmp = iostr(8,'HAM',nl,nbas,1,ekap,0,ckbas,-1,nttab,strh)
          nsites = strh%n(nbas+1)
          if (allocated(h)) deallocate(h); allocate(h(nsites*nl**4))
          call dcopy(nsites*nl**4, strh % s, 1, h, 1)
          if (allocated(h0)) deallocate(h0); allocate(h0(nsites*nl**4))
          call dcopy(nsites*nl**4, strh % s, 1, h0, 1)
          call free_str0(strh)
          if (ovlp) then
            call free_str0(strh)
            ltmp = iostr(8,'OVL',nl,nbas,1,ekap,0,ckbas,-1,nttab,strh)
            if (allocated(ov)) deallocate(ov);allocate(ov(nsites*nl**4))
            call dcopy(nsites*nl**4, strh%s, 1, ov, 1)
            call free_str0(strh)
          endif
          call poppr
        else
          if (rmaxh > d1mach(3)) then
            mxcsiz = mxnbr
            if (MOL) then
!               call defdr(ontab, nbas+1)
!               call defi(oiax,  300*nbas*niax)
              allocate(ntab(nbas+1))
!               call safe_alloc_i(ntab, nbas+1, F)
              allocate(rtab(300*nbas*3))
              allocate(rham(nspec))

              if (lscale) then
                call dcopy(nspec,rmaxh*alat/2d0,0,rham,1)
              else
                call dcopy(nspec,rmaxh/2d0,0,rham,1)
              endif
              call hpairm(nbas,s_ctrl%ips,alat,plat,s_lat%pos,rham,nsites,ntab,iax,rtab,mxcsiz,1,niax)
              deallocate(rtab, rham)
!               call rlse(oiax)
!               call defi(oiax, niax*nsites)
            else
!               call defi(ontab, nbasp+1)
              allocate(ntab(nbasp+1))
              xx = alat
              if (lscale) xx = 1
              call pairs(nbas,nbasp,xx,plat,[rmaxh/2],s_lat%pos,[-1],i,j,pgfsl,nsites,ntab,iax,mxcsiz)
            endif
          else
            if (MOL) call rx(' set RMAXH > 0 in CTRL')
            mxcsiz = 0
! !C           call defrr(ormax, nspec)
!             call ptr_ctrl(s_ctrl,1,'rmax',nspec,0,0,dum)
!             do  is = 1, nspec
!               rmax = s_spec(is)%rham
!               s_ctrl%rmax(is) = rmax/2
! !C             call dvset(s_ctrl%rmax,is,is,rmax/2d0)
!             enddo
            if (associated(s_ctrl%rmax)) deallocate(s_ctrl%rmax)
            allocate(s_ctrl%rmax(nspec))
            s_ctrl%rmax(1:nspec) = s_spec(1:nspec)%rham*0.5_8
!             call defi(ontab, nbasp+1)
            allocate(ntab(nbasp+1))
            xx = alat
            if (lscale) xx = 1
            call pairs(nbas,nbasp,xx,plat,s_ctrl%rmax,s_lat%pos,s_ctrl%ips,i,j,pgfsl,nsites,ntab,iax,mxcsiz)
          endif

          if (cmdopt('--xbs',5,0,outs) .and. itrlx == 1) then
            call tbxbs(s_spec,nbas,nspec,alat,s_lat%pos,s_ctrl%ipc,s_ctrl%dclabl)
            if (cmdopt('--bangs',7,0,outs)) then
              dist0 = 5
              if (cmdopt('--d=',4,0,outs)) then
                ip = 4
                call skipbl(outs,len(outs),ip)
                ii = parg(' ',4,outs,ip,len(outs),' ',1,1,jj,dist0)
                call rxx(ii /= 1,' TBZINT: error parsing --d=')
              endif
              call bangs(dist0,nbas,s_lat%pos,alat,s_ctrl%dclabl,s_ctrl%ipc)
            endif
          endif
        endif
!C   ... Patch iax(6) for the padded basis
        if (nbasp > nbas) call pairp6(nbas,npadl,npadr,iax,ntab)
        if (mordrn == 1) then
          call rx('pairec commented out')
!C          mxnbr = 2*rmaxs**3*nbasp
!C          call redfi(oiax, niax*mxnbr)
!C          call pairec(nbas,nbase,0,s_ctrl%ipc,alat/alat,plate,w(obas),
!C     &      w(orham),nttab,ntab,iax,w(owk),mxcsiz)
!C          call redfi(oiax, niax*nttab)
        endif
        call mkiaxd(nsites,lmx,s_ctrl%ipc,iax)
        call safe_alloc_r8(h, nl**4*nsites*min(nspc**2*nsp,4),T)

        if (sc .or. trh) call safe_alloc_r8(h0, nl**4*nsites*min(nspc**2*nsp,4),T)
        if (ovlp) call safe_alloc_r8(ov, nl**4*nsites*nspc**2, T)

        if (fitpar) then
          call safe_alloc_r8(dh, nvar*nl**4*nsites*nspc**2, T)
        elseif (hdiffs) then
          call safe_alloc_r8(dh, 4*nl**4*nsites*nspc**2, T)
          if (cryf) call safe_alloc_r8(dhcf, 4*nl**4*nsites*nspc**2, T)
          if (ovlp) then
            call safe_alloc_r8(dov, 4*nl**4*nsites*nspc**2, T)
            if(ocryf) call safe_alloc_r8(dovcf,4*nl**4*nsites*nspc**2,T)
          endif
        endif

!C   ... tb programs not converted to ntab from npr yet ...
        call npr2tb(1,nbasp,npr,ntab)
        call mk2didc(tbc, nbas, iax, niax, nsites, npr, s_ctrl%ipc, idxdn, nl)

        longrange = cmdopt('-lres',5,0,outs)
        if (longrange) then
            if (.not. allocated(molidc)) then
                allocate(molidc(nbas))
!                 Later it may be beneficial to initiale the list here and only update it later.
            end if
            molidc = 0
            call findmolidc(nbas, iax, niax, nsites, npr, molidc)
        end if

        call tcx('bmkstrx')
!C   --- Structure constants for Ewald sums ---
        if (sc .and. .not. cmdopt('--sfly',6,0,outs)) then
!          if (.not. cmpstrx) then
!           call mkstrx(tbc,ldip,nbas,nlmq1,nlmq,s_lat%pos,s_ctrl%ipc,nclasp,lmxl, &
!                      & s_lat%awald,alat,vol,dlv2,nkd,qlv2,nkq,indxcg,jcg,cg,cy,pv,MOL,strx,dstrx,plat)
!          else
          call mkstrxd(s_ctrl,s_ctrl%ipc,s_lat,tbc,nlmq, nlmq1,lmxl,ldip,dlv2,nkd,qlv2,nkq, &
                     & indxcg,jcg,cg,cy,struxd,dstrxd, struxidx, dstrxidx, struxsize, dstrxsize)
!          end if
        endif
        call tcn('bsc')

!C   --- Slater-Koster Hamiltonian ---
        if (.not. rdham) then
          call tcn('bldham')
!           call tbham(nsp,nspc,s_ctrl,s_lat,lscale,fitpar,0,dum,dum,
!      &         nlmesp,memod,decay,deccf,decov,
!      &         dcocf,poly,cutmd,cut,cutov,
!      &         iam,npm,nterm,nset,tabme,tabcf,tabov,
!      &         tbocf,nsites,npr,oiax,h,h0,dh,ov,dov,dhcf,dovcf)
! !
          if (.not. allocated(ov)) allocate(ov(1),dov(1))
          if (.not. allocated(dh)) allocate(dh(1))
          if (.not. allocated(dov)) allocate(dov(1))
          call bldham(nsp,nspc,s_ctrl,s_lat,tbc,lscale,F,0,dum,dum,nlmesp,memod,decay,decov,poly, &
                     & cutmd,cut,cutov,iam,npm,nterm,nset,tabme,tabov,nsites,npr,iax,lmx,h,dh,ov,dov)
          if (sc .or. trh) call dcopy(nl**4*nsites*min(nspc**2*nsp,4),h,1,h0,1)

!C#ifdefC TBPG
!C        call swapad(nl**2,nbas,npadl,npadr,ntab,iax,h)
!C#endif
          if (.not. allocated(delL)) allocate(delL(1))
          if (.not. allocated(h0)) allocate(h0(1))
          call tbdiag(ltb,nbasp,nl,nsp,nspc,nsp1,s_ctrl%ipc,nsites,npr, &
                        & s_ctrl%initc,s_pot%qnu,nl*nsp1,delta,delL,h,h0,ov)
          call tcx('bldham')
        endif

!C   --- Add to Hamiltonian MEs: hLL' -> hLL' + ebarLL' * sLL' ---
        if (addsll) then
          call tcn('addes')
          allocate(wk(nl*nbas))
          if (tbe) then
            call pshpr(iprint())
          else
            call pshpr(0)
          endif
          call addes(nsites,nl,nspc,nbas,ltb,npr,iax,ov,dov,wk,h,dh)
          call poppr
          deallocate(wk)
          call tcx('addes')
        endif

!C   --- Write Hamiltonian to STR file ---
        if (lpgf /= 0 .or. cmdopt('-dumph',6,0,outs) .and. lroot) then
          call npr2tb(0,nbasp,npr,ntab)
          allocate(strh % a(nl*nl*nbasp)); strh % a = 0.0_8
          strh % s => h

          ltmp = iostr(1,'STR',nl,nbasp,1,ekap,0,ckbas,-1,nsites,strh)
          nsites = strh % n(nbasp+1)
          call fclose(fopn('STR'))
          call npr2tb(1,nbasp,npr,ntab)
!           call rlse(oalpha)
          call free_str0(strh)
        endif

!C   --- Printout ---
        if (iprint() >= 50 .and. .not. pair .and. lroot) then
          call npr2tb(0,nbasp,npr,ntab)
!           call defdr(oalpha,-nbasp*nl**2)
          allocate(strh%a(nbasp*nl**2)); strh%a(:) = 0.0_8
!C     ... Printout Hamiltonian
          print *
          print *,'Real-space TB hamiltonian :'

          if (TBU) then
            do j = 1, 2
              print *
              if (j == 1) print *,'TB+U spin up'
              if (j == 2) print *,'TB+U spin down'
              os2 = 1 + (j-1)*nl**4*nsites
              call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,h(os2),1,1,0d0,0,1d0)
!C     ...  Printout Overlap
              if (ovlp .and. j == 1) then
                print *
                print *
                print *, 'Real-space Overlap matrix'
                os2 = 1 + (j-1)*nl**4*nsites
                call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,ov(os2),1,1,0d0,0,1d0)
              endif
            enddo
          else
            do  j = 1, nspc**2
              if (nspc == 2) print *
              if (j == 1 .and. nspc == 2) print *,'Spin: Up-Up, w/o spin-orbit'
              if (j == 2) print *,'Spin: Up-Down, w/o spin-orbit'
              if (j == 3) print *,'Spin: Down-Up, w/o spin-orbit'
              if (j == 4) print *,'Spin: Down-Down, w/o spin-orbit'
              os2 = 1 + (j-1)*nl**4*nsites
              call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,h(os2),1,1,0d0,0,1d0)
!C     ...  Printout Overlap
              if (ovlp) then
                print *, 'overlap matrix'
                os2 = 1 + (j-1)*nl**4*nsites
                call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,ov(os2),1,1,0d0,0,1d0)
              endif
            enddo
          endif
          if (iprint() > 50 .and. iand(ltb,16+128) /= 0) then
!C       ... Printout Hamiltonian derivatives
            print *
            print *, 'Real-space gradH :'
            do i = 1, 4
              print *
              if (i == 1) print *, 'x-component:'
              if (i == 2) print *, 'y-component:'
              if (i == 3) print *, 'z-component:'
              if (i == 4) print *, 'radial component:'
              do j = 1, nspc**2
                if (nspc == 2) print *
                if (j == 1 .and. nspc == 2) print *,'Spin: Up-Up'
                if (j == 2) print *,'Spin: Up-Down'
                if (j == 3) print *,'Spin: Down-Up'
                if (j == 4) print *,'Spin: Down-Down'
                os2 = 1 + (nspc**2*(i-1)+(j-1))*nl**4*nsites
                call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,dh(os2),1,1,0d0,0,1d0)
!C           ... Printout Crystal Field derivatives
                if (cryf) then
                  print *,'Crystal Field gradH'
                  os2 = 1+(nspc**2*(i-1)+(j-1))*nl**4*nsites
                  call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,dhcf(os2),1,1,0d0,0,1d0)
                endif
!C           ... Printout Overlap derivatives
                if (ovlp) then
                  print *,'gradO'
                  os2 = 1+(nspc**2*(i-1)+(j-1))*nl**4*nsites
                  call shostr(nl**2,nsites,nbasp,plat,s_lat%pos,0,strh%a,iax,ntab,dov(os2),1,1,0d0,0,1d0)
!C             ... Printout Overlap Crystal Field derivatives
                  if (ocryf) then
                    print *,'Crystal Field gradO'
                    os2 = 1+(nspc**2*(i-1)+(j-1))*nl**4*nsites
                    call shostr(nl**2,nsites,nbas,plat,s_lat%pos,0,strh%a,iax,ntab,dovcf(os2),1,1,0d0,0,1d0)
                  endif
                endif
              end do
            end do
          endif
          deallocate(strh%a)
          call npr2tb(1,nbasp,npr,ntab)
        endif
        if (.not. sc .and. (tbbnd .or. tbfit)) goto 3

!C#ifdefC TBPG
!C          if (lpgf /= 0) then
!CC       ... Setup energy mesh ...
!C            nzp = nint(semsh(1))
!C            call defdc(ozp,nzp)
!C            call defdc(owz,nzp)
!C            call emesh(semsh,w(ozp),w(owz))
!CC     ... ASA planar GF ...
!C            call defi(oogll,npl+2)
!C            call defi(olgii,npl+2)
!C             call defdr(opp,-6*nlspc)
!C          call rx('patch call to pgfasa')
!C             call pgfasa(plat,s_lat%pos,alat,nl,nbas,nspc,nclass,w(oiclsp)
!C     .         ,s_ctrl%nrc,clabl,wsr,rmaxs,lmx,nband,switch,w(opp)
!C     .         ,w(ovshfp),s_pot%qnu,totmom,qspirl,w(oeula),nkxyz(1)
!C     .         ,nkxyz(2),nl,norder,width,range,drange,npts,nkp,w(oqp)
!C     .         ,w(owght),semsh,w(ozp),w(owz),eband,idxdn,efermi
!C     .         ,vmtz,elin,z,s_pot%qc,ntet,w(ontet),sumev,w(opgslp)
!C     .         ,w(opgplp),gfopts,.true.,w(olgii),w(oogll))
!C            call rx0(prgnam//' done')
!C          else
!C            call rx('tbpg: pgf not set')
!C          endif
!C#endif
!C        if (lgors('ctrl quit,2',sctrl))
!C     .    call rx0(prgnam//': Q=ATOM encountered')
        if (iand(s_ctrl%quit,2)/=0) call rx0(prgnam//': Q=ATOM encountered')
!C#ifdefC STONER
!C        if (lstnr(1)/=0) then
!C          if (nsp == 2) call rx(prgnam//': for STONER=T set NSPIN=1')
!C          call defdr(ozos,iabs(npts)*nlspc)
!C          call defi(oindex,nclass)
!C          call defdr(onbar,mnpts)
!C          call defdr(oewk,mnpts)
!C          call defdr(omwk,mnpts)
!C          call defdr(ommom,-nbas)
!C          call defdr(oemag,-nbas)
!C        endif
!C#endif
        sumev = 0d0
        thrpv = 0d0
        emad = 0d0
        emg = 0d0
!C --- Vanderbilt order N ---
        if (mordrn == 2) then
          if (iand(ltb,1+2+2**7+2**13)/= 0 .or. nspc==2) &
               & call fexit2(-1,111,' Exit -1 : VDB not implemented for ltb=%i, nspc=%i',ltb,nspc)
          if (lstnr(1)/=0 .or. iand(ltb,1024+2048+4096)/=0) &
               & call fexit2(-1,111,' Exit -1 : VDB incompatible with ltb=%i, stoner=%l',ltb,lstnr)
!C#ifdefC VDB
!C          rmx(3) = rmaxr + rmaxh
!C          call defi(onp3,-3*nbas)
!C          call defi(otnp3,-3*nbas)
!C          call defi(oinorb,-nbas)
!C          stop 'convert iclass to ipc, vdbptr'
!C          call vdbptr(nbas,nspc,s_ctrl%ipc,nclass,nl,nhs,nrs,
!C     .      idxdn,w(oinorb),alat,plat,s_lat%pos,
!C     .      rmx,npr,w(onp3),w(otnp3),tnpr,
!C     .      oiax3,oihx,oirx,oirhx,ohs,odhs,odms,odm)
!C          stop 'convert iclass to ipc, shorth,shordh,shortr'
!C          call shorth(nbas,tnpr,w(otnp3),w(onp3),
!C     .      w(oiax3),nl,s_ctrl%ipc,w(oinorb),
!C     .      h,oh,ohs)
!C          call shordh(nbas,tnpr,w(otnp3),w(onp3),
!C     .      w(oiax3),nl,s_ctrl%ipc,w(oinorb),nhs,
!C     .      w(odh),odh,odhs)
!C          call shortr(nbas,tnpr,w(otnp3),w(onp3),
!C     .      w(oiax3),nl,s_ctrl%ipc,w(oinorb),
!C     .      w(odm),odm,odms)
!CC ... Can't release oh here anymore ...
!CC          call rlse(oh)
!C          call cgmin(nbas,nl,tnpr,nrs,
!C     .             oinorb,onp3,otnp3,
!C     .             oiax3,oihx,oirx,oirhx,ohs,odms,cgmtol,sumev)
!C          call hfvdb(nbas,nl,tnpr,nhs,
!C     .      oinorb,onp3,otnp3,
!C     .      oiax3,oihx,oirx,odms,odhs,of)
!C#else
          call rxx(.true.,prgnam//': Vanderbilt not implemented')
!C#endif VDB
        endif
!C --- Self-consistency loop ---
        nitmax = 1
        it0 = 1
        uconv = .not. sc
        if (sc) then
          nitmax = maxit
          if (pair) nitmax = 1
          if (nitmax == 1) uconv = T
          if (mixQ) it0 = 0
        endif
        call tcx('bsc')
  700   do  1  it = it0, nitmax
!C --- device to supress force calculation until self consistent ---
          if (defer .and. sc .and. .not. uconv) then
            if (iprint() >= 20) print *,'TBZINT: defer force calculation until self consistent'
            pv = F
            force = F
          endif
          if ((.not. sc .and. pair) .or. mordrn > 0) then
            entrpy = 0d0
            call dcopy(3*nbas,0d0,0,fr,1)
            goto 500
          endif
          if (it > 0) then
            if (.not. lls) then
              if (.not. allocated(dhcf)) allocate(dhcf(1))
              if (.not. allocated(vso)) allocate(vso(1,1),hso(1))
              if (.not. allocated(dovcf)) allocate(dovcf(1))
              if (.not. allocated(pot0)) allocate(pot0(1))
              if (.not. allocated(rhon)) allocate(rhon(1))
              if (.not. allocated(zos)) allocate(zos(1))
              if (.not. allocated(dros)) allocate(dros(1))
              call bndtb(s_ctrl,s_bz,s_lat,s_spec,tbc,nbas,nl,nspc,    &
     &          nsp,nsp1,lmx,idxdn,nclass,s_ctrl%ips,s_ctrl%ipc,       &
     &          s_ctrl%nrc,pv,force,zval,mull,                   &
     &          wtbzm,nsites,iax,npr,xyzfrz,h,h0,dh,dhcf,vso,hso,      &
     &          ov,dov,dovcf,pot0,efermi,sumev,entrpy,fr,thrpv,        &
     &          esite,rho,rhoc,rhon,zos,dros,ldim,iprmb)
            else
!-- Here goes all the alternative stuff (vanderbuilt & similar polynomials, ls++) at first to be easier and safer to experiment.
!   The communicators's setup should be simpler, the 2D comms are still present but their structure is not reaally used.
!   Shall be properly integrated once this starts being an issue.
              call bndtb_alt(s_ctrl,s_bz,s_lat,s_spec,tbc,nbas,nl,nspc,&
     &          nsp,nsp1,lmx,idxdn,nclass,s_ctrl%ips,s_ctrl%ipc,       &
     &          s_ctrl%nrc,pv,force,zval,mull,                   &
     &          wtbzm,nsites,iax,npr,xyzfrz,h,h0,dh,dhcf,vso,hso,      &
     &          ov,dov,dovcf,pot0,efermi,sumev,entrpy,fr,thrpv,        &
     &          esite,rho,rhoc,rhon,zos,dros,ldim,iprmb)
            end if
            if (bittst(ltb,2**6) .or. cmdopt('--dos',5,0,outs)) goto 3
          endif

!C     ... Electrostatic terms
          if (sc) then
            if (it > 1) then
!C             call upacks('strn mix',i1,i2)
              i1 = str_pack('mix',-2,s_strn,strn)
              dum = 0
              if (lroot) call pshpr(80)
              if (.not. parmxp(it,strn,len_trim(strn),broy,mmix,wt, &
     &          beta,dum2,dum2,dum2,dum2,outs,wc,nkill,betav,rmsdel))   &
     &          call rx(prgnam//': parse in parmxp failed')
              if (lroot) call poppr
            endif
            if (nkill == -999) then
              rmsdel = 1d-16
              goto 500
            endif
            if (nitmax > 1 .and. mixrho) then
              call tcn('mixing')
!C         ... mix density matrix (and Mulliken charges) now
              call rhomix(ovlp,nl,nsp,nbas,s_ctrl%dclabl,s_ctrl%ipc, &
     &          it,nitmax,cnvg,mmix,neltst,beta,tjmax,tj,rhl0, &
     &          rho,rho0,rhoc,rhlit,rhoit,awk,rmsdel)
              call tcx('mixing')
            endif
            call tbmpol(nbas,nsp,nl,nlmq,s_pot%qnu,s_ctrl%ipc,lmxl,gaunt,qpol,rho,rhoc,qmpol,mmom)
            if (nitmax > 1 .and. mixQ) then
!C         ... mix multipole moments now
              call tcn('mixing')
              if (mixQ2) then
                call qmix2(lc,nbas,nsp,nlmq,ltb,s_ctrl%dclabl,        &
     &            s_ctrl%ipc,it,nitmax,cnvg,cnvm,mmix,nkill,nelts,beta, &
     &            betav,tjmax,tj1,tj2,qmpol,qit,           &
     &            mmom,mmom0,mit,awk1,awk2,rmsdel,wc,broy,wt,magmix,rms2)
              else
                call qmix(lc,nbas,nsp,nlmq,ltb,s_ctrl%dclabl,         &
     &            s_ctrl%ipc,it,nitmax,cnvg,mmix,nkill,neltst,beta,   &
     &            tjmax,tj,qmpol,qit,mmom,mmom0,           &
     &            mit,awk,rmsdel,wc,broy,wt)
              endif
              call tcx('mixing')
            endif
!             if (.not. cmpstrx) then
!             call tbesel(tbc,mixQ,force,pv,s_site,ltb,nbas,nl,nlmq1,nlmq, &
!      &        nsp,nclass,lmxl,s_ctrl%ipc,s_ctrl%dclabl,idxdn,          &
!      &        s_lat%symgr,nsgrp,s_lat%istab,s_lat%pos,alat,              &
!      &        qpol,indxcg,jcg,cg,gaunt,ak,strx,dstrx,rho,mmom,           &
!      &        dros,ldim,rhoc,rhon,s_pot%qnu,stni,idu,                   &
!      &        uh,jh,qmpol,delh,pot0,ecorr,fe,fnou,fnom,ppdip,tpvq)
!             else
            if (.not. allocated(rhon)) allocate(rhon(1))
            call tbeseld(tbc, mixQ,force,pv,s_site,ltb,nbas,nl,nlmq1,nlmq,  &
     &        nsp,nclass,lmxl,s_ctrl%ipc,s_ctrl%dclabl,idxdn,        &
     &        s_lat%symgr,nsgrp,s_lat%istab,s_lat%pos,alat,            &
     &        qpol,indxcg,jcg,cg,gaunt,ak,struxd,dstrxd, struxidx, dstrxidx, &
     &        rho,mmom,dros,ldim,rhoc,rhon,s_pot%qnu,stni,idu,               &
     &        uh,jh,qmpol,delh,pot0,ecorr,fe,fnou,fnom,ppdip,tpvq)
!             endif
            if (pair) then
              entrpy = 0d0
!               call dcopy(3*nbas,0d0,0,fr,1)
              fr = 0.0_8
              goto 500
            endif


            if (nitmax > 1) then
              call tcn('mixing')
              if (mixrho .or. mixQ) then
                call tbadh2(nl,nbas,nsp,nsites,npr,delh,h0,h)
              else
!C           ... mix hamiltonian matrix elements now
                call dhmix(neltst,nl,nsp,nbas,idxdn,s_ctrl%dclabl,s_ctrl%ipc, &
                & it,nitmax,cnvg,mmix,nkill,beta,tjmax,tj,delh,delL,hdelta,awk,rmsdel)
                call tbadh1(nl,nbas,nsp,nsites,mmix,npr,hdelta,h0,h)
              endif
              call tcx('mixing')
            endif

            if (nitmax <= 1) then
!             This is just to get around the structuring of the time table.
              call tcn('mixing')
              call tcx('mixing')
            end if

          endif

!C --- Generalised Stoner model ---
!C#ifdefC STONER
!CC         if ((lgors('ctrl lstonr,1',sctrl))) then
!C          if ((iand(s_ctrl%lstonr,1))) then
!C            call dpzero(mmom,nbas)
!C            call dpzero(w(oemag),nbas)
!c             stop 'iclass in stoner'
!C            call stoner(nl,nclass,clabl,idxdn,efermi,
!C     .        drange,iabs(npts),w(ozos),w(oindex),stni,w(oammx),
!C     .        mnpts,switch(44),w(onbar),w(oewk),w(omwk),mmom,
!C     .        w(oemag))
!C          endif
!C#else
!C     if (lgors('ctrl lstonr,1',sctrl))
      if (iand(s_ctrl%lstonr(1),1)/=0) call rx('LM: recompile with ccomp -dSTONER ...')
!C#endif

  500 continue
!      call mpi_allreduce(mpi_in_place,rmsdel,1,mpi_double,mpi_max,mpi_comm_world,err)
!C --- TB total energy ---
      if (it > 0 .or. pair) then
         if (.not. allocated(emag)) allocate(emag(1))
         call tbtote(s_ctrl,s_lat,s_site,nl,nsp1,idxdn,nclass, &
     &        s_ctrl%nrc,s_ctrl%dclabl,pv,force,nbas,iax,npr, &
     &        iam,npm,lscale,v0,ppmod,poly,cutmd,cutpp,s_pot%qnu,sumev, &
     &        alpha,entrpy,emad,ecorr,rho,stni,jh,mmom,emag,esite,  &
     &        ppdip,thrpv,tpvq,fr,fe,fmax,fnou,fnom,eref,erep,etot,     &
     &        efree,emg,amgm,vol)
      endif

      if (pair) goto 600
      if (it > 1 .and. sc) then
        if (wt(3) == 0) then
          if (rmsdel < cnvg .or. it == maxit) then
            if (lroot) then
              call iodelL(ltb,0,nl,nbas,nsp,delta,delL)
            endif
            if (rmsdel > cnvg .and. .not. uconv) then
              if (iprint() > 10) call awrit3(                         &
     &            ' (warning) cnvg=%g, rms=%g, after %i iterations.',' '&
     &            ,80,i1mach(2),cnvg,rmsdel,it)
            endif
            uconv = .true.
!C --- device to supress force calculation until self consistent ---
            if (defer) then
              if ((bittst(ltb,2**7) .and. .not. pv) .or. &
     &            (bittst(ltb,2**4) .and. .not. force)) then
                it0 = 1
                nitmax = 1
                pv     = bittst(ltb,2**7)
                force  = bittst(ltb,2**4)
                goto 700
              else
                goto 600
              endif
            else
              goto 600
            endif
          endif
        elseif ((wt(3) == 1) .and. (broy == 1)) then
          if (cnvm == 0) then
            cnvm = cnvg
            if (iprint() > 30) then
              print *, "TBZINT: cnvm not set, defaulting to cnvg"
            endif
          endif
          call rxx(.not.(mixQ2),' TBZINT: need mxq2 to clock-mix')
          if (magmix == 1 .or. it == maxit) then
            if (rmsdel < cnvm .or. it == maxit) then
              if (lroot) then
                call iodelL(ltb,0,nl,nbas,nsp,delta,delL)
              endif
              if (rms2 > cnvm .and. .not. uconv .or. magmix == 0) then
                if (iprint() > 10) call awrit3(                         &
     &            ' (warning) cnvm=%g, rms=%g, after %i iterations.',' '&
     &            ,80,i1mach(2),cnvm,rms2,it)
              endif
              uconv = .true.
!C --- device to supress force calculation until self consistent ---
              if (defer) then
                if ((bittst(ltb,2**7) .and. .not. pv) .or. &
     &              (bittst(ltb,2**4) .and. .not. force)) then
                  it0 = 1
                  nitmax = 1
                  pv     = bittst(ltb,2**7)
                  force  = bittst(ltb,2**4)
                  goto 700
                else
                  goto 600
                endif
              else
                goto 600
              endif
            endif
            rms2 = rmsdel
            magmix = 0
          endif
        endif
      endif
!C   --- end self consistency loop ---
    1 continue
  600 continue
      if (tbbnd) goto 3

!C   --- Molecular statics ---
        if (nvar /= 0) then
          call query('step ',4,tstep)
          call query('accrlx ',4,accrlx)
          call relax(prgnam,s_ctrl,s_site,s_spec,itrlx,indrx,natrlx,nvar,fr,p,dum,nelts,delta,s_lat%pos,icom)
          if (lroot) then
            if (cmdopt('--mv',4,0,outs)) then
              ifi = fopn('MV')
              if (cmdopt('--st',4,0,outs) .and. itrlx == 1) then
                if (iprint() > 10) print*, 'Starting new mv file ..'
              else
                if (iprint() > 10) print*, 'Appending to mv file ..'
                call poseof(ifi)
              endif
              call awrit2('frame etot= %;4d Ry fmax=%;4da.u.',' ',128,ifi,etot,fmax)
              call iomv(nbas,s_lat%pos,alat,ifi)
              call fclose(ifi)
            endif
            if (cmdopt('--xyz',5,0,outs) .and. lroot) then
              ifi = fopn('XYZ')
              if (cmdopt('--st',4,0,outs) .and. itrlx == 1) then
                if (iprint() > 10) print*, 'Starting new xyz file ..'
              else
                if (iprint() > 10) print*, 'Appending to xyz file ..'
                call poseof(ifi)
              endif
              call awrit3('%i PLAT: %9:2;5d ALAT: %;5d',' ',144,ifi,nbas,plat,alat)
              call awrit2('etot= %;4d Ry fmax=%;4da.u.',' ',128,ifi,etot,fmax)
              call ioxyz(nbas,s_ctrl%ips,s_lat%pos,alat,vel,fr,0d0,nlmq,qmpol,s_ctrl%dclabl,sc,.true.,ifi)
              call fclose(ifi)
            endif
          endif
!C   --- Molecular dynamics ---
        elseif (md .and. .not. ipi) then
          call verlet(ensbl,itrlx,lstart,nbas,plat,nf,amass,tstep,   &
     &      fr,s_ctrl%ips,s_ctrl%dclabl,temp,efree,thrpv,vol,pext, &
     &      sc,nlmq,qmpol,taup,taub,mdtime,time,alat,s_lat%pos,      &
     &      vel,zeta,eps,veps,logs,zacc,zsqacc,ekin,tkin,cons)
          if (ensbl == 'NPT') then
            alat = alat0 * exp(eps)
            lmxst = 6
            call pshpr(0)
            call lattc(as,ewtol,rpad,alat,alat,plat,gx,gy,gz,gt,plat, &
                        & qlat,lmxst,vol,s_lat%awald,dlv2,nkd,qlv2,nkq,nkdmx,nkqmx)
            call poppr
            if (iprint() >= 10) call awrit3(' TBZINT scaling, eps=%d, alat=%d vol=%d',' ',128,i1mach(2),eps,alat,vol)
            s_lat%alat = alat
            s_lat%plat = plat
            s_lat%vol = vol
            do  j = 1, nbas
             s_site(j)%pos(1:3) = atoms(1:3,j)
            enddo
          endif

          if (lroot) ifi = fopn('STRT')

          call iostrt(nbas,nf,itrlx,s_lat%pos,vel,eps,zeta,zacc,zsqacc,veps,time,-ifi,err)
          if (lroot) call fclose(ifi)
        endif

!C   ... Write positions to file
        if (cmdopt('--wpos=',7,0,outs) .and. lroot) call iopos(T,1,outs(8:),nbasp,s_lat%pos,s_site)

!C   ... Write forces to file
        if (cmdopt('--wforce=',9,0,outs) .and. lroot) call iopos(T,0,outs(10:),nbasp,fr,s_site)

!C --- Write variables to SAVE file ---
        if (lroot) then
          ifi = fopn('SAVE')
          print *
          outs = 'mmom,sev,erep,etot,emad,3pv,TS'
          ch = 't'
          if (sc .and. uconv) ch = 'c'
          if (icom == 1) ch = 'r'
          if (md)          ch = 'm'
          call iosave(ch,outs,energy,-ifi,s_ctrl%nvario)
          call iosave(ch,outs,energy,-lgunit(1),s_ctrl%nvario)
          call fclose(ifi)
        endif
        if (icom == 1) then
          call info0(20,1,1,' '//prgnam//'%a:  Relaxation complete')
          if (MOL) then
            if (lroot) then
              dist0 = 5
              if (cmdopt('--d=',4,0,outs)) then
                ip = 4
                call skipbl(outs,len(outs),ip)
                ii = parg(' ',4,outs,ip,len(outs),' ',1,1,jj,dist0)
                call rxx(ii /= 1,' TBZINT: error parsing --d=')
              endif
              call bangs(dist0,nbas,s_lat%pos,alat,s_ctrl%dclabl,s_ctrl%ipc)
            endif
          endif
!           call rlse(onpr)
!           call rlse(opgfsl)
          deallocate(pgfsl)
          deallocate(ntab)
          deallocate(npr)

          goto 3
        endif

!         call rlse(opgfsl)
!         call rlse(onpr)

        deallocate(ntab)
        deallocate(pgfsl)
        deallocate(npr)

!C --- End of relaxation/MD loop ---
        itrlx = itrlx + 1
    2 continue
    3 continue

!       if (allocated(strx)) deallocate(strx)
!       if (allocated(dstrx)) deallocate(dstrx)
!       if (allocated(dlv2)) deallocate(dlv2)
!       if (allocated(qlv2)) deallocate(qlv2)
      if (tbe) then
        call tcx('tbzint')
        return
      endif
!C --- End of TBE branch ---

!C --- TBBND and TBFIT branch ---
      if (.not. (tbbnd .or. tbfit)) return
!C --- Get efmax, nevmx
!       if (tbbnd) then
!         nevmx = 0
!         efmax = s_bz%efmax
! !C ---   Get Fermi energy ---
!         nfilem = fopn('BAND')
!         efermi(1) = 999
!         i = getef(nfilem,111,efermi(1))
!         call info2(21,0,0,'%N Found fermi energy : %d',efermi(1),0)
!         call fclose(nfilem)
!       endif

!C --- Set up permutation matrix for Bloch ---
      if (iprint() >= 30) print *
      call makidx(nl,1,1,nbas,0,s_spec,s_ctrl%ips,-1,iprmb,inxsh)
!c...deb
!C ... Print orbital positions in hamiltonian, resolved by l
!c     if (iprint() > 50) then
!c       call showbs(ssite,sspec,nkaph,0,iprmb,ldham)
!c     endif
!c...deb
      ldim = inxsh(1)*nspc
      idim = inxsh(2)*nspc - ldim
      call rxx(idim /= 0,prgnam//': bad orbital switches')
      nlb  = ldim
      call rxx(nlb > nlbmx,'increase nlbmx in main')

      if (tbbnd) then
!         if (fxst('SYML') /= 1) then
!           print *, ' '
!           print *,' TBZINT can''t find symmetry line file!'
!           print *, ' '
!           return
!         endif
!         ifi =  fopno('SYML')
        if (cmdopt('-ef',4,1,strn)) then
!           if (.not. a2bin(strn,efermi(1),4,0,' ',4,-1)) call rxs2('TBZINF: failed to parse "',strn(1:30),'%a"')
            read(strn,*) efermi(1)
        endif
        ifib = fopnn('BNDS')
        write(ifib,'(i5,f10.5,i5)') nlb,efermi(1),0
        call fclose(ifib)
      endif
      if (tbfit) then
!C       Tony used nsp=1 for nspc=2.  Won't work until this is fixed.
        if (nspc == 2  .or. nsp == 2) call rx('need to restore nsp=2 to SO coupled case')
!C       fixed cutoffs not implemented in tbfit
        if (.not. lscale) &
     &    call rx('TBZINT: Fixed cutoff option is not implemented in the tbfit branch. Use scaled cutoffs instead.')
        ifib = fopno('BNDS')
        call pshpr(1)
!C       call suqlsr(1,ifib,nlb,.false.,.false.,nq,dum,dum)
        call suqlsr(1,ifib,max(nspc,nsp),nlb,i2,ldim,1,ldim,.false.,.false.,nq,efermi(2),dum,dum)
        call poppr
        call info2(10,1,0,' tbe --fit: fit to data in bnds file (contains %i kp, %i bands)',nq,nlb)
!         call defdr(oqp, 3*nq)
!         call defdr(oeband, -ldim*max(nspc,nsp)*nq)
        allocate(qp(3*nq))
        allocate(eband(ldim*max(nspc,nsp)*nq)); eband = 0.0_8
        call suqlsr(6,ifib,max(nspc,nsp),nlb,i2,ldim,1,ldim,.false.,.false.,nq,efermi(2),qp,eband)
        allocate(iq1(nl,nclasp))
        allocate(iq2(nl,nclasp))
        do  i = 1, nclasp
          iq1(1:nl,i) = s_spec(ics(i)) % iq1(1:nl)
          iq2(1:nl,i) = s_spec(ics(i)) % iq2(1:nl)
        enddo

!         call defi(oifit, 2*nq)
        allocate(ifit(2*nq))
!C       Until it is checked ...
        rl = .false.
        if (.not. allocated(hso)) allocate(hso(1),vso(1,1))
        if (.not. allocated(pot0)) allocate(pot0(1))
        call tbfitmrq( s_ctrl,s_lat,s_tb,nterm,nlmesp,nset,nl,nsp,nspc,nclass, &
                     & npair,nq,ldim,nelts,memod,mxnbr,idim,s_ctrl%dclabl,itab,idec,  &
                     & itbcf,idcf,itbov,idov,itocf,idocf,iq1,iq2,idxdn,cut,cutov,iam,  &
                     & npm,lmx,iprmb,rmaxh,eband,delta,qp,hso,v0,s_pot%pnu,ifit,   &
                     & tabme,decay,tabcf,deccf,tabov,decov,tbocf,dcocf,s_pot%qnu,vso,pot0,rl)
!         call rlse(oifit)
        deallocate(ifit)
        return
      endif
      if (tbbnd) then
        allocate(e(ldim*ldim*2),evl(ldim))
        nevmx = ldim
        efmax = s_bz%efmax
! !C ---   Get Fermi energy ---
!         nfilem = fopn('BAND')
!         efermi(1) = 999
!         i = getef(nfilem,111,efermi(1))
!         call info2(21,0,0,'%N Found fermi energy : %d',efermi(1),0)
!         call fclose(nfilem)
!C   ... Override Fermi level with command-line value
        if (cmdopt('--band',6,0,strn)) plbopt = strn(7:)
!         call defi(oifbls, -ldim*2*2)
        allocate(ifbls(ldim*2*2)); ifbls = 0
        nq = 0
        do ! loop over blocks of q-points
!C --- Get efmax, nevmx
          i = nsp1
          if (lso) i = 1
          call suqlst(s_lat,plbopt,0,nlb,efermi(1),i,dum,nfbn,ifbls,nq,q,onesp)
          if (nq <= 0) call rx0('TBZINT: finished bands.')
          do  iq = 1, nq
              call suqlst(s_lat,plbopt,0,nlb,efermi(1),i,dum,nfbn,ifbls,nq,q,onesp)
            do  isp = 1, nsp
!C --- Generate the eigenvalues into evl ---
              ipr = 0
              if (iprint() > 40) ipr = 3
              evl = 0.0_8
              call secmtb(s_ctrl,tbc,plat,nbas,nl,nspc,nsp,isp,lmx,s_ctrl%ipc, &
                        & iprmb,ldim,nevmx,efmax,iq,nq,q,nsites,iax,    &
                        & npr,h,vso,hso,ov,pot0,sc,ipr,tbbnd,nev,e,evl)
              call suqlsw(nlb,isp,nsp1,evl)
              if (nfbn(1) /= 0) then
                call ztoyy(e,nlb,nlb,nlb,nlb,0,1)
                call suqlse(nlb,isp,nsp1,nlb,1,nfbn,ifbls,nlb,e,evl)
                if (nfbn(2) /= 0) call suqlse(nlb,isp,nsp1,nlb,2,nfbn,ifbls,nlb,e,evl)
              endif
            enddo
          enddo
        enddo
        deallocate(evl,e)
      endif

      call del_cart_ctx(tbc % c3d)
      call del_cart_ctx(tbc % c2d)
      call del_cart_ctx(tbc % c1d)
      call tcx('tbzint')
      end subroutine tbzint


      subroutine safe_alloc_i(a,s,z)
         implicit none
         integer, intent(inout), allocatable :: a(:)
         integer, intent(in) :: s
         logical, intent(in) :: z

         if (.not. allocated(a)) then
            allocate(a(s));
         else
            if (size(a) < s) then
               deallocate(a)
               allocate(a(s))
            end if
         end if
         if (z) a = 0
      end subroutine safe_alloc_i

      subroutine safe_alloc_r8(a,s,z)
         implicit none
         real(8), intent(inout), allocatable :: a(:)
         integer, intent(in) :: s
         logical, intent(in) :: z

         if (.not. allocated(a)) then
            allocate(a(s));
         else
            if (size(a) < s) then
               deallocate(a)
               allocate(a(s))
            end if
         end if
         if (z) a = 0.0_8
      end subroutine safe_alloc_r8

      subroutine free_str0(s)
         use structures, only : str_str0
         implicit none
         type(str_str0), intent(inout) :: s
         if (associated(s%i )) deallocate(s%i )
         if (associated(s%n )) deallocate(s%n )
         if (associated(s%a )) deallocate(s%a )
         if (associated(s%s )) deallocate(s%s )
         if (associated(s%ng)) deallocate(s%ng)
      end subroutine free_str0
