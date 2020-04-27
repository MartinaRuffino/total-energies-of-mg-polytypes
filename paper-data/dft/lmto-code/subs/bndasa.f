      subroutine bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,
     .  s_pot,s_str,s_optic,vorb,dmatu,efermi,nit,ftmod,nvfit,nfcof,
     .  ivcof,gband,eband,nbmax,nevmx,qnu,sumev,rhos,amag,aamom)
C- k-integration of bands over BZ for moments.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin lncol lsx lham lasa lgen3
Ci                 loptc nbasp lrel nclasp ipc
Co     Stored:     lham lasa
Co     Allocated:  *
Cio    Elts passed: lscr lasa ldos lham ipc dclabl ics lstonr ips rmax
Cio                lgen3 lqp ncomp idcc nrc
Cio    Passed to:  nmpot secmat secmtn secmt2 pp2alp asaddq asaopm
Cio                optint optin2 optinq intdrv
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat nsgrp nkd nkq awald vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:vol pos dlv qlv cg jcg indxcg symgr
Cio    Passed to:  secmat secmtn sstrxq asaddq asaopm
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  idu lmxb lmxa a nr z rmt hcr name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  getidu sumlst nmpot secmat secmtn mullmf iorbtm
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class v0
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  getidu sumlst nmpot lmfitmr6 mullmf
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet lmet n range w nevmx efmax lio ndos
Ci                 dosw zval fsmom qp lshft
Co     Stored:     numq qp ef
Co     Allocated:  qp wtkb swtk
Cio    Elts passed:qp wtkp swtk wtkb idtet ipq n w def
Cio    Passed to:  subzi optint optin2 optinq intdrv
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  qss neula eterms nbf nmto nlibu lmaxu udiag lsig
Ci                 ldham kmto ndhrs
Co     Stored:     lham eterms
Co     Allocated:  *
Cio    Elts passed:eula bdots iprmb magf nprs iaxs hrs
Cio    Passed to:  nmpot secmat secmtn sstrxq asaddq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nrhos vmtz ves
Co     Stored:     *
Co     Allocated:  ppn
Cio    Elts passed:qpp pp sop bxc ppn pti pprel grrme
Cio    Passed to:  nmpot secmat secmtn secmt2 asaddq asaopm
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:alph nkaps iax kaps s sdot nitab adot
Cio    Passed to:  secmat secmtn secmt2 pp2alp
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window nchi2 axes esciss lpart iq
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  asaopm optint optin2 optinq intdrv
Ci Inputs
Ci   vorb  :(LDA+U) orbital dependent potential matrices
Ci   dmatu :(LDA+U) density matrix for LDA+U
Ci   nit   :iteration number (used for printout only)
Ci   ftmod :mode for fitting bands; see doc/levenberg-marquardt.html
Ci   nvfit :Used with ftmod
Ci         :nvfit>0 flags that gband is to be calculated.
Ci         :see nvar in lmfitmr.f
Ci   nfcof :Used with ftmode.  Number of coefficients to vary
Ci   ivcof :Used with ftmode.  Index to coefficients to vary
Ci   nbmax :leading dimension of eband
Co Outputs
Co   gband :(Used when nvfit>0)
Co         :gradient of band energies wrt particular variables.
Co   efermi:Fermi energy
Co   eband :energy bands; alias eb (sec*.f)
Co   nevmx :largest number of eigenvectors found
Co   qnu   :energy-weighted moments of the sphere charges
Co   sumev :sum of eigenvalues
Co   rhos  :spin density-matrix
Co   amag  :system magnetic moment
Co   aamom :local magnetic moments
Co   dmatu :site density-matrix corresponding to LDA+U blocks
Cs Command-line switches
Cs   --band      : Tabulate energy bands; see doc/Command-line-options.html
Cs   --bonly     : Use in conjuction with --lxb switch
Cs   --chklim    : Check band limits in optics calculations
Cs   --cv:       : Calculate electronic specific heat, eV
Cs   --cvK:      : Calculate electronic specific heat, Ry
Cs   --invbl     : Not documented
Cs   --jdosw     : Channels for optical properties; See doc/optics.html
Cs   --jdosw2    : Channels for optical properties; See doc/optics.html
Cs   --lxb       : Check band limits in optics calculations
Cs   --mull      : Mulliken analysis; see doc/Command-line-options.html
Cs   --oldbz     : Not documented
Cs   --onesp     : Generate bands for one spin; use with --band
Cs   --opt:read  : Read optical matrix elements; See doc/optics.html
Cs   --opt:write : Write optical matrix elements; See doc/optics.html
Cs   --pdos      : Partial DOS; see doc/Command-line-options.html
Cs   -ef=        : Overwrite Fermi level; use with --band
Cl Local variables
Cl   nspc: number of coupled spins (1 unless noncoll and/or S-O coupling).
Cl   nspx: number of independent spin channels containing eigenvalues
Cl         and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl   nev:  number of eigenvalues (and eigenvectors) calculated
Cl   nchan:number of l (or lm) channels summed over all classes,
Cl         e.g. nlo*nclass.
Cl   liomom:0 suppress writing dos weights to moments file
Cl         :1 write dos weights to moments file to be used in the
Cl            generation of output density
Cl         :2 write partial dos weights to moments file, for lmdos
Cl         :3 write Mulliken weights to moments file, for lmdos
Cl         :4 Same as 2, but parallel mode
Cl   lwtkb :0 weights are neither required a priori nor available
Cl         :1 weights are required a priori, and are available
Cl         :  (having read from disk or from a prior band pass)
Cl         :-1 weights are required a priori, but were not read
Cl   lnsph=1 make coffs for nonspherical density
Cl   ldos =a combination of the following integers
Cl         1 make dos
Cl         2 generate weights for partial dos
Cl         4 generate weights for m-decompos'n of pdos
Cl   lrho=T make the on-site spin-density matrix, needed for mag. forces.
Cl   ldens=T makes the density matrix (under construction)
Cl   lsx1   if 1, a response function will be computed
Cl          from eigenvectors created by bndasa; bndasa must generate
Cl          eigenvectors appropriate for this purpose.
Cl   nfilem is the logical units for moments file.
Cr Remarks
Cm  Memory:
Cm   Without exploiting hermiticity, and without folding down, the
Cm   minimum storage requirement is three complex matrices (H, O and Z)
Cm   dimension ndim^2. The last of these is not needed until the work-
Cm   space for the real-space structure constants has been released.
Cm   If ndim is separated into lower and intermediate sets
Cm   (ndim = ldim + idim) then if the i-waves remain in the basis the
Cm   requirement is 3ldim^2 + 6(ldim*idim) + 3idim^2. If the i-waves are
Cm   folded down the requirement becomes 3ldim^2 + 2(ldim*idim) + idim^2
Cm   but these must be declared simultaneously with the work-space for
Cm   the real-space structure constants.
Cf  Files:
Cf    For each qp, bndasa writes eband, and then accmom (and also accsm
Cf    if spin-pol bands coupled), provided metal=T and lwtkb=F.
Cf    If metal=F, writes doswt instead.
Cf    See iomoms for structure of the moments file.
Cb Bugs
Cb   lwtkb should be patterned after bndfp.
Cb   ldos should be decoupled from lwtkb (see makwts).  Requires caller to
Cb   allocate appropriate memory for accwt.
Cu Updates
Cu   09 Aug 13 Option to calculate bands, sumev by tetrahedron,
Cu             weights for density by sampling with lmet modes 1,2,3
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 Apr 11 Improved treatment of orbital-dependent noncollinearity
Cu   07 Apr 11 Passes work arrays from secmat to asaddq in memory,
Cu             to replace dist read
Cu   01 Mar 11 --pdos works in noncollinear case
Cu   07 Jan 11 fixed-spin-moment compatible w/ fp; proper total energy
Cu   17 Nov 10 L-M fitting passes through ftmod
Cu   15 Oct 10 Read strux through rdstrx
Cu   07 Aug 10 qpp moments with with MPIK
Cu   01 Jul 10 gradient of bands works properly for spin-pol case
Cu   05 Mar 10 Implemented Mulliken analysis and site-resolved partial dos
Cu   11 Sep 09 Reworked optics
Cu   11 Aug 09 Changes in energy bands wrt ASA parameters, first cut
Cu   26 Oct 08 (ATP) first attempt at inverse Bloch transform for band fitting
Cu   04 Jun 08 Orbital moments now generated when lmet=2 or lmet=3
Cu   08 Nov 07 (J. Xu) New LDA+U
Cu   08 Jun 06 optional weights for band plotting
Cu   13 Jul 05 (wrl) Nonlinear optics re-instated
Cu   27 Apr 04 Optics branch works with MPI
Cu   25 Apr 04 Reads efermi when lmet=2 and qp weights available.
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   20 Feb 04 (U.J.Rho) Optics code can to on-the-fly sampling integration
Cu   23 Sep 03 SX patterned after GW.  sigm(rs) obtained
Cu             by call to rdsigh.
Cu   06 Sep 03 (jek) extended optics code to metals
Cu   10 Mar 03 Poke double-counting from appl. field into eterms
Cu    1 Mar 03 Revised calculation of output magnetization
Cu             in noncollinear case; new double counting
Cu   27 Jan 03 In noncollinear mode, always accumulates density
Cu             matrix and prints output magnetization from it.
Cu    9 Jan 03 (K. Belashchenko) parallelize bndasa over k-points
Cu   30 May 02 Allow reading SX sigma in alpha repsn
Cu   18 Jan 02 Redesigned to handle accumulation of m-resolved weights
Cu             New command-line argument --pdos; see sumlst for syntax
Cu   03 Feb 01 patched code for m-resolved moments.
Cu   28 Apr 98 optical matrix elements and dielectric eps added
Cu    4 May 98 calls atfold so high not automatically folded to interm.
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
      integer nit,nvfit,nbmax,mpsord,ftmod,nfcof,ivcof(nfcof),nevmx
      double precision efermi,sumev,amag(3),qnu(*),
     .  rhos(2,*),aamom(*),eband(nbmax,*),gband(nvfit,nbmax,*)
      double complex vorb(*),dmatu(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_bz)::    s_bz
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_optic):: s_optic
C ... Dynamically allocated arrays
      integer,allocatable :: lmxa(:)
C ... Local parameters
      character strn*256,plbopt*256,sjdosw(2)*120
      logical metal,lrhos,ldens,savvec,T,F
      logical :: spintexture = .false. ! .true. when color weights made by spin texture
      integer idalloc,allocvb,fopna,fopno,iobzwt,iomoms,iopq,iprint
      integer ltet,ldham(16),lmdim,lwtkb,nlo,nrhos,stdo,bitand,fopn,i,
     .  idim,ifi,ikp,iprt,iscr,iskp,isp,isw,izval,j,k,kkp,lasa,lasas,
     .  lbzio,ldim,ldimx,ldos,lgen3,lgunit,lham,lhams,lidim,lihdim,
     .  liomom,lmet,lmfitmr5,lncol,lnsph,loptic,lrout,lrsig,lswtk,lsx,
     .  lsx1,lteto,moddos,n1,n2,n3,nbas,nbf,nbmaxx,nclass,ndos,ndos0,
     .  neul,nev,nevsav,nevx,nfbn(4),nfilem,nfstg,njdosw,nkabc(3),nkp,nl,
     .  nlspc,nmom,nmto,npol,nppn,nqpp,nschan,nsp,nspc,nspec,nspx,ntet,
     .  plbnd,stdl,onesp
      integer nflink,NULLI
      parameter (NULLI=-99999)
      integer,allocatable:: ivfit(:,:)
C     Used for site-resolved partial dos (--pdos) or mulliken (--mull)
      integer nsitmx,lmxch,nchan,nchmx
      parameter (nsitmx=1024)
      integer lsite(nsitmx)
      parameter (nppn=12)
      double precision dum(20),avw,alat,plat(3,3),qspirl(4),fsmom(3),
     .  efloc(2),efmax,swidth,srnge,dosw(2),qp(3),zval,ef0,swtkb,
     .  evtop,ecbot,eterms(22),vmtz(2)
      real(8), target :: xx(10)
      real(8), pointer :: qpsave(:,:)
      parameter(T=.true., F=.false.)
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (ldimx,ldham(5))
      equivalence (nspx,ldham(8))
      character*80 outs

      integer, allocatable :: idu(:)
      integer, allocatable :: chan(:)
      real(8), allocatable :: accwt(:)
      real(8), allocatable :: wk(:)
      real(8), allocatable :: strx(:),hrs(:),ors(:)
      real(8), allocatable :: bxcnw(:),doswq(:),ppsav(:)
      complex(8), allocatable :: accsm(:)
      complex(8), allocatable :: zll(:)
      complex(8), allocatable :: s(:)
      complex(8), allocatable :: h(:)
      complex(8), allocatable :: h2(:)
      complex(8), allocatable :: s2(:)

      type(dp_wk_vec) :: s_wk(3)

      procedure(logical) cmdopt,bittst,a2bin

C ... For optics
      integer nfilo,nfiup,nemlo,nemup,ocrng(2),unrng(2),nfilm,nempm,
     .  ne(3),npopt
      double precision optrng(4),dwopt(3)
      real(8),allocatable :: jdosw(:,:,:,:),wk3(:),doswtq(:,:,:),dos(:)
      real(8),allocatable :: orbtm(:),dnpp(:),optmt(:),velmt(:),eps(:)
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      integer,allocatable :: iwk1(:),iwk2(:),ifbls(:)

C ... For LDA+U
      integer nlibu,lmaxu,ludiag

C#ifdefC STONER
C      integer l,nchan,odoswt,oeband
C#endif
C ... For the density matrix
      integer nsite
C ... For symmetry operations
      integer nsgrp
C     integer nstar,oqstar,ormat,igstar,ig
C     double precision qlat(3,3)
C     bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) == bit)
C ... External calls
      external amagnc,asaddq,asaopm,awrit0,awrit1,awrit2,awrit3,awrit7,
     .  bdotsr,bzints,bzwtsf,dcopy,defdc,defdr,defrr,dosio, dpdump,
     .  dpscop,dpzero,dscal,dvset,fclose,fclr,fexit2,info, info0,iomomn,
     .  iomomq,iorbtm,rdstrx,isanrg,lsets,makdos,maknos,moment,
     .  nmpot,optint,pack1,pack2,poppr,pshpr,query,rlse,rsmsym,
     .  rx,rx0,rx1,secmat,setpr,shostr,subzi,sumlst,suqlst,suqlsw,
     .  togpr,upack,upack1,upack2, xxxdif
      integer procid,master,numprocs

      integer ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision buf
C      integer length(0:MAX_PROCS-1)
C      integer offset(0:MAX_PROCS-1)
      double precision sttime,entime
      logical mlog
      integer i1,i2,iq,isqp

      integer, dimension(:),allocatable :: kpproc,offset,length
      real(8),allocatable :: evll(:,:),bufv(:)

C --- Unpack parameters and parameter setup ---
      call tcn('bndasa')
      master = 0
      call mpi_comm_rank( mpi_comm_world, procid, ierr )
      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
      stdo = lgunit(1)
C     stde = lgunit(1)
      stdl = lgunit(2)
      liomom = 1
      if (numprocs > 1) then
        liomom = 0
        mlog = cmdopt('--mlog',6,0,strn)
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',i)
        namelen(procid) = i-1
      endif
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
      lmet = s_bz%lmet
      mpsord = s_bz%n
      srnge = s_bz%range
      swidth = s_bz%w
      nevmx = s_bz%nevmx
      efmax = s_bz%efmax
      lbzio = s_bz%lio
      ndos0 = s_bz%ndos
      dosw = s_bz%dosw
      zval = s_bz%zval
      fsmom(1:2) = s_bz%fsmom
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      qspirl = s_ham%qss
      neul = s_ham%neula
      eterms = s_ham%eterms
      nbf = s_ham%nbf
      nmto = s_ham%nmto
      nlibu = s_ham%nlibu
      lmaxu = s_ham%lmaxu
      ludiag = s_ham%udiag
      call dvset(eterms,17,18,0d0)
      nrhos = s_pot%nrhos
      lncol = s_ctrl%lncol
      lsx = s_ctrl%lsx
      lrsig= s_ham%lsig
      iscr = mod(s_ctrl%lscr,10)
      lsx1 = mod(lsx,2)
      if (mod(iscr,2) == 1) lsx1 = 1
C     Hang onto lham, lasa to restore after exiting bndasa
      lhams = s_ctrl%lham
      lasas = s_ctrl%lasa
      lgen3 = s_ctrl%lgen3
      loptic = s_ctrl%loptc
      lnsph = isw(IAND(s_ctrl%lasa,32) /= 0)
      vmtz = s_pot%vmtz
      ldham = s_ham%ldham
      nfbn = 0
      call dpzero(s_bz%sopertp,12)
      if (lnsph /= 0 .or. nlibu > 0) then
        i = nl**2
        nqpp = (i*(i+1))/2
        call zinit(s_pot%qpp,nqpp*4*nsp*nbas)
      else
        nqpp = 0
      endif
      if (allocated(lmxa)) deallocate(lmxa)
      allocate(lmxa(nclass))
      call spec2class(s_spec,nclass,s_ctrl%ics,'lmxa',1,lmxa,xx)

      if (lmet > 3) call rx1('BNDASA: METAL=%i not implemented',lmet)
      ef0 = 0
C     Flag indicating Fermi level is not known
      efermi = -99
      nullify(qpsave,s_wk(1)%p,s_wk(2)%p,s_wk(3)%p)
C     Debugging
C     call shoctl(s_ctrl,s_spec,s_pot,101,stdo); stop

CC ... pp's in alpha representation, consistent with strux
C      call pp2alp(0,s_ctrl,s_str,lham,s_pot%pp)

C ... Get nldau for ldau
      if (nlibu /= 0) then
        allocate(idu(4*nbas))
        call getidu(nbas,s_spec,s_site,idu)
      else
        allocate(idu(1))
      endif

C ... Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn)) then
        plbnd = 1
        plbopt = strn(7:)
        nevmx = -1
        loptic = 0
        s_bz%numq = 1
        allocate(ifbls(ldimx*2*2)); call iinit(ifbls,ldimx*2*2)
        lmet = 0
        lsx1 = 0
      else
        allocate(ifbls(1))
        plbnd = 0
      endif
C ... Band fitting mode (set plbnd=2 as switch for later)
      if (nvfit > 0) then
        if (plbnd == 1)
     .    call rx('BNDASA: cannot use fitting mode with --band')
        plbnd = 2
        allocate(ivfit(nfcof,2))
      endif

      nspc = 1
      if (bitand(lncol,1+2+4+8) /= 0) nspc = 2
C     lrhos  = bittst(lncol,16) .and. nspc == 2 .and. nevmx >= 0
      lrhos  = nspc == 2 .and. nevmx >= 0
      if (cmdopt('--pdos',6,0,outs)) lrhos = .false.
      metal  = lmet /= 0
      ltet = mod(s_ctrl%lmet/2,2) + 10*mod(s_ctrl%lmet/16,2)
      ldos   = IAND(s_ctrl%ldos,7)
      nsite = 0
C ... Partial DOS: get sites, number of channels
      if (cmdopt('--mull',6,0,strn).or.cmdopt('--pdos',6,0,strn)) then
        nsgrp = s_lat%nsgrp
        efmax = 1d3
        nevmx = ldimx
        nchmx = min(1024,nbas*nl**2)
        allocate(chan(nchmx))
        j = 0
        if (cmdopt('--mull',6,0,strn)) j = 1
        i = 0
        call sumlst(j,nchmx,nbas,nsgrp,s_site,s_spec,strn(7:),moddos,
     .    nsite,lsite,lmxch,nchan,lmdim,chan,i)
C       call redfi(ochan,lmdim*nsite)
        liomom = 3
        if (cmdopt('--pdos',6,0,strn)) then
          if (cmdopt('--mull',6,0,strn))
     .    call rx('--pdos and --mull not allowed in conjunction')
          if (ldos < 4 .and. moddos == 2 .or. moddos == 5)
     .      ldos = ldos+4
          liomom = 2
          if (numprocs > 1) liomom = 4
        endif
        nmom = nspc
C       Only can use pdos with metal=2 => lwtkb>0 and nschan=4
        call iomomn(.true.,2,.false.,1,nspc,nmom,1,nfstg)
C       nschan should be 1 in collinear case, 4 in noncollinear case
        nschan = mod(nfstg/10,10)
        if (liomom == 4) then
          allocate(doswtq(nchan,ldimx,max(nschan,nsp)*nkp))
        endif
        if (procid == master) then
          nfilem = fopna('moms',-1,4)
          i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .      ldim,nfstg,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
        endif
      else
        allocate(chan(1))
      endif
      if (nbf == nl*nl .or. neul == nl*nl) then
        if (ldos < 4) ldos = ldos+4
        call info0(30,1,0,
     .    ' bndasa: m-resolved noncollinearity: resolve DOS by m')
      endif
      nlo = nl
      if (ldos >= 4) nlo = nl*nl
      if (liomom /= 1 .and. lmet == 1) call rx('BNDASA: METAL=1 '
     .  //'incompatible with suppressing I/O to moms file')
      if (nsite /= 0 .and. lmet == 1 .and. iprint() >= 10) then
        call info(10,1,0,' ***WARNING*** bndasa: partial DOS mode not'//
     .    ' compatible with metal=1:',0,0)
        call info(10,0,0,' Program will not accumulate output moments.',
     .    0,0)
      endif
      if (nsgrp /= 1 .and. nlo /= nl .and. iprint() >= 10) then
        call info(10,1,0,' ***WARNING*** bndasa: ASA does not '//
     .    'symmetrize m-resolved moments.',0,0)
        call info(10,0,0,' Suppress symops for m-resolved DOS.',0,0)
      endif

      ldens  = bittst(lbzio,8)
      onesp  = isw(cmdopt('--onesp',7,0,outs))
      ndos   = iabs(ndos0)
      nlspc = nl * nsp * nclass
      nevsav = -1
      lrout  = 1
      if (nevmx < 0) lrout = 0
      nbmaxx = nbmax
      if (nspc == 2) nbmaxx = nbmax*nsp

C --- Setup pertaining to screened exchange ---
      if (lrsig+lsx1 /= 0) then
C   ... Turn on Gamma rep, ortho z; KKR rho; no downfolding
        if (lsx1 == 1) then
          s_ctrl%lham = IOR(s_ctrl%lham,128+32) -
     .                  IAND(s_ctrl%lham,4)
          s_ctrl%lasa = s_ctrl%lasa - IAND(s_ctrl%lasa,256)

          nevmx = nbmax
          efmax = 999
        endif
C   ... SX -> Gamma rep
        if (lrsig /= 0) then
          s_ctrl%lham = IOR(s_ctrl%lham,128)
        endif
      endif

C     Unpack (possibly altered) lham, lasa
      lham = s_ctrl%lham
      lasa = s_ctrl%lasa
      s_ham%lham = lham

C --- Printout at entry ---
      if (nit > 0) then
        call info8(30,1,0,
     .  ' --- BNDASA : iteration      ('//
     .    '%?#n#sigma, ##'//
     .    '%?#n#two-c-H, ##'//
     .    '%?#n#ortho-z, ##'//
     .    '%?#n==128#H-gamma, ##%-1j'//
     .    '%?#n==640#H-gamma-bar, ##'//
     .    '%?#n#KKR-qout, ##'//
     .    '%?#n#NMTO=%-1j%i, ##'//
     .    '%?#n#noncol, ##'//
     .    '%b%?#%c==,#%b)#%b# ---'//
     .    '%24p%i%f%c%o%a',
     .    lrsig,bitand(lham,1),bitand(lham,32),bitand(lham,128)+bitand(lasa,512),
     .    256-bitand(lasa,256),nmto,lncol,nit)
      else
        call info8(30,1,0,
     .  ' --- BNDASA : band pass ('//
     .    '%?#n#sigma, ##'//
     .    '%?#n#two-c-H, ##'//
     .    '%?#n#ortho-z, ##'//
     .    '%?#n#H-gamma, ##'//
     .    '%?#n#KKR-qout, ##'//
     .    '%?#n#NMTO=%-1j%i, ##'//
     .    '%?#n#noncol, ##'//
     .    '%b%b%?#p>24#)## ---',
     .    lrsig,bitand(lham,1),bitand(lham,32),bitand(lham,128),
     .    256-bitand(lasa,256),nmto,lncol,0)
      endif

C ... Starting Euler angles
      if (bittst(lncol,1)) then
        if (iprt(2) /= iprt(1)) then
          call pshpr(iprt(2))
        else
          call pshpr(iprint()-0)
        endif
C       call setpr(50)
        call amagnc(nbas,nl,s_ctrl%ipc,rhos,nrhos,qnu,s_ham%eula,neul,
     .    1,amag,aamom,xx)
        if (bittst(lncol,8)) then
          call bdotsr(nbas,nl,s_ctrl%ipc,rhos,nrhos,qnu,
     .      s_ham%bdots,lihdim,s_ham%iprmb,1,dum)
C         ehterm(5) = dum(1)
          eterms(17) = dum(1)
        endif
        if (neul > 1 .or. bittst(lncol,8)) then
          if (.not. bittst(lncol,8)) nbf = 99
          if (.not. associated(s_ham%magf)) s_ham%magf => xx
          call bsrhos(nbas,nl,s_ctrl%ipc,rhos,nrhos,qnu,s_pot%pp,
     .      s_pot%sop,s_ham%eula,neul,s_pot%bxc,s_ham%magf,nbf,lihdim,
     .      s_ham%iprmb,1,dum)
          eterms(17) = dum(2)
        endif

        call poppr
      endif

C ... Sanity checks
      if (loptic /= 0) then
        if (nevmx < ldim) then
          call info2(2,1,0,' BNDASA (warning) optics calculation with'//
     .      ' nevmx (%i) < hamiltonian dimension (%i)',nevmx,ldim)
        endif
        if (nevmx < 0) call rx('BNDASA: optics sought but NEVMX<0')
      endif

C ... Allocate space for density matrix.  For now, same range as s
      if (ldens) then
        call rx('bndasa: density matrix is not yet built')
      endif

C --- Setup for optics ---
      lteto = 0
      if (loptic /= 0) then
        lteto = s_optic%ltet
        ocrng = s_optic%ocrng(1:2)
        unrng = s_optic%unrng(1:2)
        dwopt = s_optic%dw
        optrng = s_optic%window
        i = 1
        if (dwopt(2) /= 0) i = 11
        call dosmsh(-i,optrng,dwopt,0,npopt,xx)

        if (metal) then
          if (nfilo <= 0) nfilo = 1
          if (nfiup <= 0) nfiup = ldimx
          if (nemlo <= 0) nemlo = 1
          if (nemup <= 0) nemup = ldimx
        else
          if (nfilo <= 0) nfilo = 1
          if (nfiup <= 0) nfiup = (nspc*nint(zval)+1)/2
          if (nemlo <= 0) nemlo = nfiup+1
          if (nemup <= 0) nemup = ldimx
        endif
        nemup = min(nemup,ldimx)
        nfiup = min(nfiup,ldimx)
        if (nfiup > nemup) then
          call info2(30,0,0,' (warning, optics) last filled '//
     .      'state (#%i) above last empty state (#%i)',nfiup,nemup)
          nfiup = nemup
        endif
        if (nemlo < nfilo) then
          call info2(30,0,0,' (warning, optics) first empty '//
     .      'state (#%i) below first filled state (#%i)',nemlo,nfilo)
          nemlo = nfilo
        endif
        s_optic%ocrng(1:2) = ocrng; s_optic%ocrng(3:4) = ocrng
        s_optic%unrng(1:2) = unrng; s_optic%unrng(3:4) = unrng
        nfilm = nfiup-nfilo+1
        nempm = nemup-nemlo+1
        s_optic%nfilm = nfilm
        s_optic%nempm = nempm

C       Allocate memory to hold optical matrix elements
        if (loptic /= 0 .and. lteto == -1) then ! On the fly sampling optics
          i = 3*nfilm*nspx
          allocate(optmt(i*nempm)); call dpzero(optmt,i*nempm)
          allocate(velmt(i)); call dpzero(velmt,i)
          i = npopt*(1+nspx*3); allocate(eps(i)); call dpzero(eps,i)
        elseif (loptic > 0) then
          i = 3*nfilm*nkp*nspx
          k = idalloc('optmt',allocvb()+2,2*i,nempm)
          allocate(velmt(i)); call dpzero(velmt,i)
          allocate(optmt(i*nempm)); call dpzero(optmt,i*nempm)
        else
          allocate(optmt(1))
        endif
        call ivset(nfbn,1,4,0)
      endif

C ... Parse switch '--jdosw': get list of projections
      sjdosw = ' '
      njdosw = 0
      if (loptic < 0 .and. cmdopt('--jdosw',7,0,strn)) then
        if (strn(1:8) == '--jdosw2')
     .    call rx('BNDASA: --jdosw must precede --jdosw2')
        sjdosw(1) = strn(8:)
        call getjdosw(1,8,strn,nfbn(1),nfbn(3),i,i)
        allocate(iwk1(maxval(nfbn)))
        allocate(iwk2(maxval(nfbn)))
        call getjdosw(3,8,strn,nfbn(1),nfbn(3),iwk1,iwk2)
        call ilst2a(iwk1,nfbn(1),outs)
        call ilst2a(iwk2,nfbn(3),strn)
        call info0(20,1,0,' BNDASA:  Mulliken projection of JDOS ')
        call info2(20,0,0,'%10fgroup:  '//
     .    '%?#n#%-1j%i channel(s): '//outs//'%a (occ)  ##' //
     .    '%?#n#%-1j%i channel(s): '//strn//'%a (unocc)  ##',
     .    nfbn(1),nfbn(3))
        if (loptic == -5 .or. loptic == -6) then
          if (nfbn(3) /= 0)
     .      call rxi(' unocc weights not allowed for loptic =',loptic)
        endif
        deallocate(iwk1,iwk2)
        allocate(jdosw(1,nfilm+nempm,nkp,nspx))
        njdosw = 1
      else
        allocate(jdosw(1,1,1,1))
      endif
      if (loptic < 0 .and. cmdopt('--jdosw2',8,0,strn)) then
        sjdosw(2) = strn(9:)
        call getjdosw(1,9,strn,nfbn(2),nfbn(4),i,i)
        allocate(iwk1(maxval(nfbn)))
        allocate(iwk2(maxval(nfbn)))
        call getjdosw(3,9,strn,nfbn(2),nfbn(4),iwk1,iwk2)
        call ilst2a(iwk1,nfbn(2),outs)
        call ilst2a(iwk2,nfbn(4),strn)
        call info2(20,0,0,'%10fgrp 2:  '//
     .    '%?#n#%-1j%i channel(s): '//outs//'%a (occ)  ##' //
     .    '%?#n#%-1j%i channel(s): '//strn//'%a (unocc)  ##',
     .    nfbn(2),nfbn(4))
        if (loptic == -5 .or. loptic == -6) then
          if (nfbn(4) /= 0)
     .      call rxi(' unocc weights not allowed for loptic =',loptic)
        endif
        deallocate(iwk1,iwk2)
        deallocate(jdosw)
        allocate(jdosw(2,nfilm+nempm,nkp,nspx))
        njdosw = 2
      endif
      call dpzero(jdosw,size(jdosw))

C --- Setup for BZ integration ---
      if (plbnd == 0) then
C       call pshpr(1)
        call subzi(s_bz,lmet,ltet,lrout > 0,1,ldim,nsp,nspc,nkp,
     .    zval,nevmx,lwtkb,efermi,lswtk,ef0)
C       call poppr
C       bndasa takes eband as a passed argument
      else
        nkp = 0
        ldos = 0
        loptic = 0
        lwtkb = 0
        nevmx = -1
        if (plbnd == 2) nevmx = ldimx
        nfstg = 1
        call ptr_bz(s_bz,8+1,'wtkb',1,0,xx)
      endif

C --- NMTO setup ---
      if (lgen3 /= 0) then
        if (nevmx >= 0 .and. iprint() >= 30) then
          call awrit0(' (warning) no output density yet in 3rd gen',
     .      outs,80,stdo)
        endif
        nevmx = -1
        nfstg = 1

C   ... NMTO Potential parameters
        call ptr_pot(s_pot,1,'ppn',nppn*lihdim*nmto*nsp,0,xx)
        call togpr
        call nmpot(40,s_ctrl,s_spec,s_site,s_ham,s_pot,lidim,
     .    lihdim,s_ham%iprmb,s_ctrl%dclabl,s_pot%ppn)
        call togpr
C       s_pot%oppn = oppn
      endif

C --- More k-point independent local arrays and initialization ---
      if (lrhos) call dpzero(rhos,2*3*nrhos*2*2*nclass)
      call dvset(eband,1,nbmax*nsp*nkp,9999d0)
      if (loptic /= 0 .and. nevmx < nemup) then
        if (iprint() >= 10) call awrit2(' BNDASA (warning): '//
     .    'nevmx reset from %i to %i',' ',80,stdo,nevmx,nemup)
        nevmx = nemup
      endif
      if (lrout /= 0) call dpzero(qnu,3*nlspc)
      if (nspc == 2) then
        allocate(orbtm(nlo*2*nclass)); call dpzero(orbtm,nlo*2*nclass)
      else
        allocate(orbtm(1))
      endif
      savvec = lsx1 == 1 .or. bittst(lham,64) .and. nevmx >= 0

C --- Start loop over k points; also, re-entry for second band pass ---
      nevsav = nevmx
   99 continue

      allocate(zll(ldimx**2))
      nevmx = nevsav
      if (lwtkb == -1) nevmx = -1
C     Information to write to moments file
      nmom = nspc
      if (lwtkb == 0 .and. metal) nmom = nspc+2
      call iomomn(metal,ldos,lrhos,nevmx,nspc,nmom,lwtkb,nfstg)

      if (lwtkb == -1)
     .  call info(20,0,0,' Start first of two band passes ...',0,0)

C      if (lwtkb == 1) then
C        call setpr(110)
C      endif

C ... Case bands at usual qp list
      if (plbnd == 0) then

C       Moments file setup
        if (liomom == 1) then
          nfilem = fopna('MOMS',-1,4)
          rewind nfilem
          write (nfilem) nl, nsp, nspc, nkp, ldim, nfstg
        elseif (liomom == 0 .or. liomom == 3  .or. liomom == 4) then
          nfilem = 0
        endif
C       Eigenvector file setup
        if (savvec) then
          i = fopna('EVEC',-1,4)
          rewind (i)
          write(i) nkp,nsp
        endif

      if (numprocs > nkp) call rxi('MPIK job cannot allocate more processors than nkp =',nkp)

C ... Case generate bands at fit points
      elseif (plbnd == 2) then

        if (procid == master) then
        ifi = fopno('refbnds')
        call suqlsr(21,ifi,max(nspc,nsp),i,j,1,1,1,.false.,.false.,
     .    nkp,xx,xx,xx)
        endif
        call mpibc1(nkp,1,2,.false.,'bndasa','nkp')
        if (nkp <= 0) call rx0('bndasa: no bands to fit')

C       Make backup for s_bz%qp, and allocate new one (nkp may change)
        qpsave => s_bz%qp; nullify(s_bz%qp)
        call ptr_bz(s_bz,1,'qp',3,nkp,xx)
C        print *, 'procid',procid,size(qpsave),size(s_bz%qp)
C        call rx('done')

        if (procid == master) then
        call suqlsr(22,ifi,max(nspc,nsp),i,j,1,1,1,.false.,.false.,nkp,
     .    xx,s_bz%qp,xx)
        endif
        call mpibc1(s_bz%qp,3*nkp,4,.false.,'bndasa','qp')
        call dpzero(gband,nvfit*nbmax*nsp*nkp)

C ... Case only bands at supplied qp generated: setup
      else

C   ... Try and read Fermi level from file
        if (procid == master) then
        ifi = fopna('wkp',-1,4)
        i = iobzwt(1,i,i,i,ef0,xx,ifi)
        call fclr('wkp',ifi)
        endif
        call mpibc1(ef0,1,4,.false.,'bndasa','ef0')
C   ... Override Fermi level with command-line value
        i = 4
        if (cmdopt('-ef=',i,0,strn)) then
          if (.not. a2bin(strn,ef0,4,0,' ',i,-1)) call
     .      rxs2('BNDASA: failed to parse "',strn(1:30),'%a"')
        endif
        i = nspx
        if (onesp /= 0) i = 1
        iopq = 0
C       suqlst in MPIK mode: return cumulative number of k-points
        if (numprocs > 1) iopq = 2
C       In parallel mode, suqlst call only serves to generate nkp
        if (procid == master) then
          call suqlstst(s_lat,plbopt,iopq,ldimx,ef0,i,xx,nfbn,ifbls,nkp,qp,onesp,spintexture)
C          if (nfbn(1) < 0) then
C            spintexture = .true.
C            nfbn(1) = -nfbn(1)
C          endif
        endif
        call mpibc1(nkp,1,2,.false.,'bndasa','nkp')
        call mpibc1(nfbn,4,2,.false.,'bndasa','nfbn')
        call mpibc1(onesp,1,2,.false.,'bndfp','onesp')
        if (nkp <= 0) call rx0('bndasa: finished bands')
C MPIK: Setup to assemble all k-points into single list with qp array ---
        if (numprocs > 1) then
        if (nfbn(1) > 0 .or. nfbn(2) > 0) then
          call rx('Cannot use color weights with MPIK')
        endif
C       Return which mode suqlst uses
        call suqlsm(i)
C       Re-allocate qp and evl arrays
        call ptr_bz(s_bz,1,'qp',3,nkp,xx)
        allocate(evll(ldimx,nsp*nkp))
        call info2(20,1,1,
     .    ' bndasa:  MPIK band plotting mode %i:  %i q-points',i,nkp)
C       Loop through all qp; accumulate vector of qp.
C       Use i2 in place of nkp to preserve nkp
        if (procid == master) then
C       iq = running index to big qp list, i1 = index to current line
          iq = 0
          i2 = 0
  199     continue
          i = 1
          call pshpr(1)
          iopq = 1
          call suqlst(s_lat,plbopt,iopq,ldimx,ef0,i,xx,nfbn,ifbls,i2,qp,onesp)
          call poppr
          if (i2 > 0) then
            call pshpr(1)
            do  i1 = 1, i2
              iq = iq+1
              iopq = 1
              call suqlst(s_lat,plbopt,iopq,ldimx,ef0,i,xx,nfbn,ifbls,i2,qp,onesp)
              call dpscop(qp,s_bz%qp,3,1,3*iq-2,1d0)
            enddo
            call poppr
            call suqlsm(i)
            if (i /= 2 .and. i /= 3) goto 199
          endif
        endif
        call mpibc1(s_bz%qp,3*nkp,4,.false.,'bndasa','qp')
        if (numprocs > nkp) call rxi('MPIK job cannot allocate more processors than nkp =',nkp)
        endif
      endif
      call pshpr(iprint())
      nevx = 0
      izval = zval/2
      evtop = -9999
      ecbot = -evtop

C ... Loop over k-points

      if (.not. allocated(strx)) then
        allocate(strx(1),hrs(1),ors(1))
      endif
      if (.not. allocated(kpproc)) allocate (kpproc(0:numprocs))
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_KLOOP,procid,"k-loop")
C#endif
      call info0(30,1,0, ' ... Start MPI k-loop')

      call dstrbp(nkp,numprocs,1,kpproc(0))
      sttime = MPI_WTIME()

      else
      kpproc(0:1) = [1,nkp+1]

      if (cmdopt('--invbl',7,0,outs)) then
        call rx('bndasa 941 update iostr --invbl')
      endif

C ... Skip to optics, reading data from file
      if (cmdopt('--opt:read',10,0,strn)) then
        goto 499
      endif

      if (fsmom(1) /= 0 .and. fsmom(1) /= NULLI .and. lmet < 2) then
        call rxi('FSMOM requires BZ_METAL>1, but METAL=',lmet)
      endif
      end if

      do  ikp = kpproc(procid), kpproc(procid+1)-1
        if (numprocs > 1) then
C       print *, 'start ikp',procid,ikp
        if (ikp == kpproc(procid)) then
          if (mlog) then
            call gettime(datim)
            call awrit4(' bndasa '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting k-points %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
          endif
        endif
        endif

      iskp = nsp*(ikp-1)
      call setpr(iprt(1))
      if (numprocs > 1) then
      if (ikp > 1 .and. mod(ikp,100) /= 1) call setpr(iprt(1)-26)
      isqp = nspx*(ikp-1)
      else
      if (ikp > 1 .and. mod(ikp,100) /= 1) call setpr(iprt(1)-6)
      endif

C ... Get qp either from qp list or read from file list
      if (plbnd == 2 .or. plbnd == 0) then
        call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
      else
        iskp = 0

        if (numprocs > 1) then
          call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
C         qp = [.1d0,.2d0,.3d0]; print *, 'set qp',qp
        else
          i = nspx
          if (onesp /= 0) i = 1
          call suqlst(s_lat,plbopt,0,ldimx,ef0,i,xx,nfbn,ifbls,nkp,qp,onesp)
          if (nfbn(1) > 0) then
            nevmx = ldimx
            efmax = 9999
          endif
        endif
      endif

C ... debugging to retrieve info about star of q
C      call gtstr(111,ikp,qp,s_lat%symgr,s_lat%ag,igstar,s_bz%wtkp,nl,plat,
C     .  nbas,s_lat%pos,nstar,w(oqstar),w(ormat),s_lat%istab)
C      cycle

      do  isp = 1, nsp
C   ... No second spin for coupled-spin case
        if (nspc == 2 .and. isp == 2) cycle
C   ... Skip over second spin in collinear antiferromagnetic case
        if (onesp /= 0 .and. isp == 2) then
          call dcopy(ldim,eband(1,1+iskp),1,eband(1,2+iskp),1)
          cycle
        endif

C   --- Eigenvalues and vectors of ASA hamiltonian ---
C       Allocate memory for temporary files
        idim = lidim - ldim
        if (nevmx > 0) then
          i = ldim
          if (bittst(lncol,2)) i = i*2
          allocate(s_wk(1)%p(max(i*idim*2,1)))
          allocate(s_wk(3)%p(max(i*ldim*2,1)))
        endif
        if (lncol /= 0.and.idim > 0) then
          allocate(s_wk(2)%p(2*ldimx**2))
        endif
        allocate(s(ldimx**2))
        call secmat(s_ctrl,s_spec,s_ham,s_pot,s_lat,s_str,ikp,
     .    nkp,qp,s_bz%wtkp,isp,lrsig,nevmx,efmax,nlibu,lmaxu,idu,
     .    ludiag,vorb,nev,zll,s,s_wk,eband(1,isp+iskp),strx,hrs,ors)
        deallocate(s)
        if (lncol /= 0.and.idim > 0) deallocate(s_wk(2)%p)
        evtop = max(evtop,eband(izval,isp+iskp))
        ecbot = min(ecbot,eband(izval+1,isp+iskp))

        if (cmdopt('--invbl',7,0,outs)) cycle

C       Shift bands by beff, if present
        if (fsmom(2) /= 0 .and. nspc /= 2) then
           call bzbfde(lswtk == 1,1,1,nbmax,ldimx,(3-2*isp)*fsmom(2),
     .      s_bz%swtk,eband(1,isp+iskp))
          if (iprint() >= 30) then
          call info5(0,0,0,' after Beff  kpt %i of %i, k=%3:2,5;5d',
     .      ikp,nkp,qp,0,0)
          j = min(9,ldimx)
          if (iprint() >= 35) j = ldimx
          write(stdo,'(255(9f8.4:/))') (eband(i,isp+iskp), i=1,j)
          endif
        endif

        if (plbnd == 0) then

C          if (check /= 0) then
C          call fixef0(zval,isp,nspc,nevx,nbmax,eband,dosw,ef0)
C          if (ef0+.1d0 > eband(max(nev,1),isp) .and. nevmx >= 0 .and.
C     .        lmet /= 0) then
C            call awrit3('%N evl(nev=%i)=%;3d but '//
C     .        'ef0=%;3d ... restart with larger efmax or nevmx',
C     .        ' ',80,stdo,nev,eband(max(nev,1),isp),ef0)
C            call rx('bndfp')
C          endif
C          endif
C
C          if (ef0+.1d0 > evl(max(nev,1),isp) .and. nevmx >= 0 .and.
C     .        lmet /= 0) then
C            call awrit3('%N evl(nev=%i)=%;3d but '//
C     .        'ef0=%;3d ... restart with larger efmax or nevmx',
C     .        ' ',80,stdo,nev,evl(max(nev,1),isp),ef0)
C            call rx('bndfp')
C          endif

C   --- Calculate gradient in band energy wrt parameters in perturbation theory
        elseif (plbnd == 2) then

C     ... Vector of tb screening parameters to transform pps
          i = 6*nl*nsp*nclass
          allocate(ppsav(i))
          call dcopy(i,s_pot%pp,1,ppsav,1) ! keep original

          allocate(h(ldimx**2),s(ldimx**2),h2(ldimx**2),s2(ldimx**2))
          if (lncol /= 0.and.idim > 0) then
            allocate(s_wk(2)%p(2*ldimx**2))
          endif
          do  i = 1, nvfit
            call lmfitmr6(10,s_site,nl,nbas,nsp,s_ham%iprmb,
     .        ldham,i,nfcof,ivcof,nflink,ivfit)

C           Bug if different members of linked list have diff spins
C           if (lmfitmr5(i,nvfit,isp,nspc,ivfit0) == 0) cycle
            if (lmfitmr5(1,nfcof,isp,nspc,ivfit) == 0) cycle

C      ...  H, S for parameter + eps
C           Require pp's in orthogonal rep'n for num. diff.
            if (ivfit(1,1) <= 6) then
             call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp)
            endif
C           call prmx('pp',s_pot%pp,6,6,nl*nsp*nclass)
            call lmfitmr2(10*ftmod+2,nflink,xx,nfcof,nl,nsp,ivfit,
     .        s_pot%pp)
C           call prmx('pp',s_pot%pp,6,6,nl*nsp*nclass)
C           Restore pp's to tb repsn
            if (ivfit(1,1) <= 6) then
              call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp)
            endif
            call secmat(s_ctrl,s_spec,s_ham,s_pot,s_lat,s_str,
     .        ikp,nkp,qp,s_bz%wtkp,isp,lrsig,-2,efmax,nlibu,lmaxu,idu,
     .        ludiag,vorb,nev,h,s,s_wk,eband(1,isp+iskp),strx,hrs,ors)
C           Restore parameters to unperturbed values
            if (ivfit(1,1) <= 6) then
              call dcopy(6*nl*nsp*nclass,ppsav,1,s_pot%pp,1)
            endif

C      ...  H, S for parameter decremented by eps
C           Require pp's in orthogonal rep'n for num. diff.
            if (ivfit(1,1) <= 6) then
              call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp)
            endif
C           call prmx('pp',s_pot%pp,6,6,nl*nsp*nclass)
            call lmfitmr2(10*ftmod+4,nflink,xx,nfcof,nl,nsp,ivfit,
     .        s_pot%pp)
C           call prmx('pp',s_pot%pp,6,6,nl*nsp*nclass)
C           Restore pp's to tb repsn
            if (ivfit(1,1) <= 6) then
             call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp)
            endif
            call secmat(s_ctrl,s_spec,s_ham,s_pot,s_lat,s_str,
     .        ikp,nkp,qp,s_bz%wtkp,isp,lrsig,-2,efmax,nlibu,lmaxu,idu,
     .        ludiag,vorb,nev,h2,s2,s_wk,eband(1,isp+iskp),strx,hrs,
     .        ors)
C           Restore parameters to unperturbed values
            if (ivfit(1,1) <= 6) then
              call dcopy(6*nl*nsp*nclass,ppsav,1,s_pot%pp,1)
            endif

C       ... Gradient of band energy
C           Given e,z = evals and evecs of (h,s):
C           The change in e is given by 1st order perturbation theory as
C             sum_kl (z+)_ik dh_kl z_li - (z+)_ik ds_kl z_li e_i

            call daxpy(2*ldimx**2,-1d0,h2,1,h,1)
            call dscal(2*ldimx**2,.5d0,h,1)
C           call yprm('dH',12,h,ldimx*ldimx,ldimx,ldimx,ldimx)
            call ygemm('C','N',ldimx,ldimx,ldimx,1d0,zll,ldimx**2,
     .        ldimx,h,ldimx**2,ldimx,0d0,h2,ldimx**2,ldimx)
C           call yprm('z+ dH',2,h2,ldimx*ldimx,ldimx,ldimx,ldimx)
            call ygemm('N','N',ldimx,ldimx,ldimx,1d0,h2,ldimx**2,
     .        ldimx,zll,ldimx**2,ldimx,0d0,h,ldimx**2,ldimx)
C           call yprm('z+ dH z',2,h,ldimx*ldimx,ldimx,ldimx,ldimx)

            call daxpy(2*ldimx**2,-1d0,s2,1,s,1)
            call dscal(2*ldimx**2,.5d0,s,1)
C           call yprm('dS',12,s,ldimx*ldimx,ldimx,ldimx,ldimx)
            call ygemm('C','N',ldimx,ldimx,ldimx,1d0,zll,ldimx**2,
     .        ldimx,s,ldimx**2,ldimx,0d0,s2,ldimx**2,ldimx)
C           call yprm('z+ dS',2,s2,ldimx*ldimx,ldimx,ldimx,ldimx)
            call ygemm('N','N',ldimx,ldimx,ldimx,1d0,s2,ldimx**2,
     .        ldimx,zll,ldimx**2,ldimx,0d0,s,ldimx**2,ldimx)
C           call yprm('z+ dS z',2,s,ldimx*ldimx,ldimx,ldimx,ldimx)

            call lmfitmr3(i,ldim,h,s,eband(1,isp+iskp),nvfit,
     .        gband(1,1,isp+iskp))

          enddo
          if (lncol /= 0.and.idim > 0) deallocate(s_wk(2)%p)
          deallocate(h,s,h2,s2,ppsav)
          if (nevmx > 0) deallocate(s_wk(1)%p,s_wk(3)%p)
          cycle  ! to loop over spins

        else

         if (numprocs > 1) then
C          isp should always be 1 in noncollinear case; isqp = nspx*(ikp-1)
           call dpscop(eband(1,isp),evll(1,isp+isqp),ldimx,1,1,1d0)
C          call dpscop(eband(1,isp),evll,ldimx,1,1+ldimx*(isp-1+isqp),1d0)
         else
          call suqlsw(ldimx,isp,nspx,eband(1,isp+iskp))
          if (nfbn(1) /= 0) then
            if (nmto == 0) call ztoyy(zll,ldimx,ldimx,ldimx,ldimx,0,1)
            do  j = 1, 3
              k = j ; if (spintexture) k = 0
              if (nfbn(j) > 0) call suqlse(0,ldimx,isp,nspx,ldimx,k,nfbn,ifbls,ldimx,1,zll,zll)
            enddo
          endif
          endif
          cycle  ! to loop over spins
        endif

C       Sanity check
        if (loptic /= 0 .and. lwtkb >= 0 .and. nev < nemup) then
          call rx('BNDASA: too few eigenvectors generated'//
     .      '...increase NEVMX or EFMAX in BZ')
        endif

C   ... Write energy bands to moms file and jump to end of k-loop if no evecs
        nevx = max(nevx,nev)
        if (nevmx <= 0 .and. liomom == 1) then
          write (nfilem) 0, ldimx
          call dpdump(eband(1,isp+iskp),ldimx,-nfilem)
          cycle ! to loop over spins
        endif

C   --- Mulliken analysis and partial DOS ---
        if (lwtkb /= -1 .and. procid == master) then
        if (cmdopt('--mull',6,0,strn)) then
          call rxx(nev /= ldimx,
     .      'Mulliken analysis requires all eigenvectors')
          i = nchan*nevmx*nspc
          allocate(dos(i)); call dpzero(dos,i)
          if (nmto == 0) call ztoyy(zll,ldimx,ldimx,ldimx,ldimx,0,1)
          call mullmf(nbas,s_site,s_spec,s_ham%iprmb,zll,ldim,nspc,ikp,
     .      isp,moddos,nsite,lsite,nchan,chan,lmdim,ldim,dos)
          if (nmto == 0) call ztoyy(zll,ldimx,ldimx,ldimx,ldimx,1,0)
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
          j = fopna('moms',-1,4)
          i = iomoms(-j,nl,nsp,nspc,nkp,ldim,i,nspc,1,1,ldimx,
     .      nevmx,nchan,nchan,nev,eband(1,isp+iskp),0d0,dos,
     .      0d0,0d0,0d0)
          deallocate(dos)
        endif
        endif

C   --- Contribution to output density and related quantities from qp --
C   ... Save eigenvectors
        if (savvec) then
          i = fopna('EVEC',-1,4)
C         Record: eb header.  Read 1st eb with mc -r:br:s=1
          write(i) ldimx,1,1
C         Record: eband
          call dpdump(eband(1,isp+iskp),ldimx,-i)
C         Record: evec header.  Read 1st evec with mc -r:br:s=3
          write(i) ldimx,nev,2
C         Record: evecs
          call dpdump(zll,ldimx**2*2,-i)
        endif

        if (lrout /= 0) then
        if (loptic > 0) then
          allocate(dnpp(nl**2*nbas*nspc*4*nev))
        else
          allocate(dnpp(1))
        endif
        if (nsite > 0)  then
          allocate(doswq(nchan*nev*8))
        else
          allocate(doswq(1))
        endif
        call asaddq(s_ctrl,s_lat,s_ham,s_pot,ikp,nkp,isp,nev,
     .    nbmax,s_wk,eband(1,isp+iskp),zll,nlibu,lmaxu,idu,vorb,
     .    nfilem,lwtkb,s_bz%wtkb,s_bz%wtkp,zval,metal,loptic > 0,nqpp,
     .    ldos,lrhos,moddos,nchan,chan,lmdim,nsite,lsite,0,qnu,
     .    s_pot%qpp,dnpp,rhos,orbtm,doswq)
        if (nsite > 0 .and. nev > 0 .and. liomom == 4) then
          nschan = mod(nfstg/10,10)
          j = isp+iskp
          if (nspc == 2) j = (ikp-1)*nschan + 1
          call dcopy(nchan*nev*nschan,doswq,1,doswtq(1,1,j),1)
        endif
        deallocate(doswq)
        endif

C       call yprm('zll',2,zll,ldimx*ldimx,ldimx,ldimx,ldimx)

C   ... Optical matrix elements
        if (loptic /= 0 .and. lwtkb /= -1) then
          if (nemup > nev) call fexit2(-1,111,' Exit -1 BNDASA: '
     .      //'optics needs %i evecs but calculated only %i',nemup,nev)
        endif
        if (loptic > 0 .and. lwtkb /= -1) then
          kkp = ikp
C         On-the-fly sampling
          if (lteto == -1) then
            kkp = 1
            if (isp == 1) then
              call dpzero(optmt,3*nfilm*nempm*nspx)
              call dpzero(velmt,3*nfilm*nspx)
            endif
          endif
          call asaopm(s_ctrl,s_pot,s_lat,s_optic,nl,isp,nsp,nspc,nclass,
     .      nbas,lmxa,ldimx,nev,qp,kkp,nkp,eband(1,isp+iskp),nbmaxx,
     .      dnpp,nfilm,nempm,optmt,velmt)
        endif
        if (allocated(dnpp)) deallocate(dnpp)

C       Free temporary arrays generated by secmat
        if (nevmx > 0) deallocate(s_wk(1)%p,s_wk(3)%p)

C   ... Color weights to project jdos
        if (loptic < 0 .and. lwtkb /= -1 .and. maxval(nfbn) > 0) then
          if (nmto == 0) call ztoyy(zll,ldimx,ldimx,ldimx,ldimx,0,1)
          allocate(wk3(ldimx))
          do  i = 1, njdosw
            allocate(iwk1(maxval(nfbn)))
            allocate(iwk2(maxval(nfbn)))
            call getjdosw(3,1,sjdosw(i),nfbn(i),nfbn(i+2),iwk1,iwk2)
            if (nfbn(i) /= 0) then    ! initial states, or just states
C             call info2(0,0,0,' weight list %n:1i',nfbn(i),iwk1)
              call suqlsz(ldimx,1,nfbn(i),iwk1,1,zll,wk3)
              call dcopy(nfilm,wk3(nfilo),1,jdosw(i,1,ikp,isp),njdosw)
            else
              jdosw(i,1:nfilm,ikp,isp) = 1
            endif
            if (nfbn(i+2) /= 0) then  ! final states, if they exist
C             call info2(0,0,0,' weight list %n:1i',nfbn(i+2),iwk2)
              call suqlsz(ldimx,1,nfbn(i+2),iwk2,1,zll,wk3)
              call dcopy(nempm,wk3(nemlo),1,jdosw(i,1+nfilm,ikp,isp),
     .          njdosw)
            else
              jdosw(i,1+nfilm:nfilm+nempm,ikp,isp) = 1
            endif
            deallocate(iwk1,iwk2)
          enddo
          deallocate(wk3)
        endif

C   ... On-the-fly sampling integration of optics
        if (loptic /= 0 .and. lwtkb /= -1 .and.
     .      lteto == -1 .and. nevmx > 0) then
          npol = 3
          if (loptic < 0) npol = 1
          call optint(s_optic,s_bz,lteto,eband,nbmaxx,isp,nspx,
     .      nspc,ikp,efermi,s_bz%idtet,s_bz%wtkp,s_lat%vol,nfilm,nempm,
     .      npopt,npol,optmt,eps)
        endif

C ...
      enddo ! Loop over spins
      enddo ! k-point loop

      if (nvfit > 0) deallocate(ivfit)

      if (cmdopt('--invbl',7,0,outs)) then
        call rx('bndasa 1353 update iostr --invbl')
C        ckbas = cksumf(s_lat%pos,3*nbas)
C        call rx('bndasa 1354 update iostr')
C        call iostrx(8,'STR',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,os)
C        nsite = s_str%ntab(nbas+1)
C        mxorb = nglob('mxorb')
CC ... symmetrise real space sqrdel S sqrdel
C        call rsmsym(0,plat,mxorb,s_ham%iprmb,ldim,nbas,s_lat%pos,nl,nsp,
C     .              1,nsite,w(ontab),w(oiax),s_lat%symgr,s_lat%istab,
C     .              nsgrp,0,nl*nl,strx,w(ostrxs))
CC ... symmetrise real space H
C        call rsmsym(0,plat,mxorb,s_ham%iprmb,ldim,nbas,s_lat%pos,nl,nsp,
C     .              1,nsite,w(ontab),w(oiax),s_lat%symgr,s_lat%istab,
C     .              nsgrp,0,nl*nl,hrs,w(ohrss))
CC ... symmetrise real space O
C        call rsmsym(0,plat,mxorb,s_ham%iprmb,ldim,nbas,s_lat%pos,nl,nsp,
C     .              1,nsite,w(ontab),w(oiax),s_lat%symgr,s_lat%istab,
C     .              nsgrp,0,nl*nl,ors,w(oorss))
C        call defdr(odon, nbas*nl**4)
C        call wrirsh(T,T,bittst(lham,128),nl,nsp,nclass,nbas,plat,
C     .    s_lat%pos,alat,s_ctrl%dclabl,s_ctrl%ipc,s_pot%pp,w(oalph),
C     .    nsite,w(oiax),w(ostrxs),w(ohrss),w(oorss),w(oh1rs),
C     .    w(oh2rs),w(oo2rs),w(odon))
C        call rx('bndasa 1376 update iostr')
C        call iostrx(1,'HAM',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,ohrss)
C        call iostrx(1,'HAM1',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,oh1rs)
C        call iostrx(1,'HAM2',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,oh2rs)
C        call iostrx(1,'OVL',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,oorss)
C        call iostrx(1,'OVL2',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
C     .             oiax,ontab,oo2rs)
C        call fexit(0,9,'Real space hamiltonian written to disc',0)
      endif
      call poppr
      deallocate(zll)

C ... Gather the results related to bands from all processes
      if (numprocs > 1) then
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',(entime-sttime),0)

C ... Band plotting mode: dump bands to disk
      if (plbnd == 1) then
        call info0(20,1,0,' Writing bands to bands file (MPIK) ...')
        call mpibc4(evll,evll,kpproc,ldim*nsp,4)
C       iq = running index to big qp list, i1 = index to current line
        if (procid == master) then
          iq = 0
          i2 = 0
  299     continue
          i = nspx
          if (onesp /= 0) i = 1
          iopq = 0
C         Sets up suqlst for writing
          call suqlst(s_lat,plbopt,iopq,ldimx,ef0,i,xx,nfbn,ifbls,i2,qp,onesp)
          if (i2 <= 0) call rx0('bndasa')
          do  i1 = 1, i2
            iq = iq+1
            isqp = nspx*(iq-1)
            call suqlst(s_lat,plbopt,1,ldimx,ef0,i,xx,nfbn,ifbls,i2,qp,onesp)
            call dpscop(qp,s_bz%qp,3,1,3*iq-2,1d0)
            do  isp = 1, nspx
              call dpscop(evll(1,isp+isqp),eband(1,isp),ldimx,1,1,1d0)
C              call dpscop(w(oevl),eband(1,isp),ldimx,
C     .          1+ldimx*(isp-1+isqp),1,1d0)
C              if (mod(i1,10) /= 1) call pshpr(iprint()-6)
              call info5(30,0,0,' bndasa:  kpt %i of %i, k=%3:2,5;5d',i1,nkp,qp,0,0)
              if (mod(i1,10) /= 1) call poppr
              call suqlsw(ldimx,isp,i,eband(1,isp))
            enddo
          enddo
          if (i /= 3) goto 299
        endif
        call rx0('done')
      endif

      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0

      if (procid == master) then
        call info0(20,0,-1,' Sharing data between processes...')
        sttime = MPI_WTIME()
      endif
      call MPI_ALLREDUCE(evtop,buf,1,mpi_real8,MPI_MAX,MPI_COMM_WORLD,ierr)
      evtop = buf
      call MPI_ALLREDUCE(ecbot,buf,1,mpi_real8,MPI_MIN,MPI_COMM_WORLD,ierr)
      ecbot = buf

C ... Collect eband from various processors
      call mpibc4(eband,eband,kpproc,nbmax*nsp,4)
      if (plbnd == 2) then
        call mpibc4(gband,gband,kpproc,nvfit*nbmax*nsp,4)
      endif

C     print 974, procid, eband(1,1:nkp)
C 974 format(i4/(6f12.6))
      endif

C ... Case generating bands: find next block of qp
      if (plbnd == 2) then
C       Restore original s_bz%qp
        if (associated(qpsave)) then
          deallocate(s_bz%qp); s_bz%qp => qpsave
        endif
        deallocate(strx,hrs,ors,chan,idu,ifbls,orbtm)
        if (allocated(eps)) deallocate(eps)
        if (nevmx > 0) then
        endif
        call tcx('bndasa')
        return
      endif

C --- Gather the results related to c.d. from all processes ---
      if (numprocs > 1) then
C ... Optical matrix elements
      if (lwtkb /= -1 .and. lrout > 0 .and. mod(loptic,10) > 0 .and. lteto >= 0) then
        call rx('bug in generating optmt with MPIK')
C       if (procid == master) print *, kpproc
C       call snot(procid,optmt,3*nfilm*nempm*nspx,nkp)
        call mpibc4(optmt,optmt,kpproc,3*nfilm*nempm*nspx,4)
      endif

C ... Collect doswtq from various processors; store on disk as 1 file
      if (lwtkb /= -1 .and. lrout > 0 .and. liomom == 4) then
        j = nchan*ldimx*max(nschan,nsp)
        call mpibc4(doswtq,doswtq,kpproc,j,4)
        if (procid == master) then
C       print 388, doswtq(nchan,ldimx,:)
C 388   format(f15.10)
        nfilem = fopna('MOMS',-1,4)
        nschan = mod(nfstg/10,10)
        do  ikp = 1, nkp
          iskp = nsp*(ikp-1)
          do  isp = 1, nsp
C           No second spin for coupled-spin case
            if (nspc == 2 .and. isp == 2) cycle
            j = isp+iskp
            if (nspc == 2) j = (ikp-1)*nschan + 1
            j = iomoms(-nfilem,nl,nsp,nspc,2,ldim,nfstg,nschan,1,1,
     .        ldimx,ldimx,nchan,0,ldimx,eband(1,isp+iskp),xx,doswtq(1,1,j),xx,xx,xx)
          enddo
        enddo
      endif
      endif

C     Add the qnu contributed by various processors
      if (lwtkb /= -1 .and. lrout > 0) then
        call mpibc2(qnu,3*nlspc,4,3,.false.,'bndasa','qnu')
      endif

C     Add the rhos contributed by various processors
      if (lwtkb /= -1 .and. lrout > 0 .and. lrhos) then
        call mpibc2(rhos,2*3*nrhos*2*2*nclass,4,3,.false.,'bndasa','rhos')
      endif

      if (lwtkb /= -1 .and. lrout > 0 .and. nqpp > 0) then
        call mpibc2(s_pot%qpp,2*nqpp*4*nsp*nbas,4,3,.false.,'bndasa','qpp')
      endif

C     Add the response function contributed by various processors
      if (lwtkb /= -1 .and. lrout > 0 .and. lteto == -1) then
        i = (1+nspx*3)*npopt
        allocate(bufv(i))
        call MPI_ALLREDUCE(eps,bufv,i,mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(i,bufv,1,eps,1)
        deallocate(bufv)
      endif

      deallocate(kpproc, stat=ierr)
      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)

      if (procid == master) then
        entime = MPI_WTIME()
        call info2(20,0,0,' MPI broadcast took %;1d sec',(entime-sttime),0)
      endif

C     Printout of eigenvalues
      call info0(30,0,0,' Printout of eigenvalues ...')
      call pshpr(iprint())
      do  ikp = 1, nkp
        call dpscop(s_bz%qp,qp,3,3*ikp-2,1,1d0)
        if (ikp > 1 .and. mod(ikp,100) /= 1) call setpr(iprt(1)-6)
        iskp = nsp*(ikp-1)
C       call setpr(iprt(1))
        do  isp = 1, nsp
C     ... No second spin for coupled-spin case
          if (nspc /= 2 .or. isp /= 2) then
            j = min(9,ldimx)
            if (iprint() >= 35) j = ldimx
            if (iprint() >= 30) then
              call awrit3(' kpt %i of %i, k=%3:2,5;5d',
     .          ' ',80,stdo,ikp,nkp,qp)
              write(stdo,'(255(9f8.4:/))') (eband(i,isp+iskp), i=1,j)
            endif
            if (procid == master) then
              call awrit3(' kpt %i of %i, k=%3:2,5;5d',
     .          ' ',80,stdl,ikp,nkp,qp)
              write(stdl,'(255(9f8.4:/))') (eband(i,isp+iskp), i=1,j)
            endif
          endif
        enddo
      enddo
      call poppr
      endif

C --- Extract dmatu ---
      if (nlibu > 0) then
        call asadmu(nbas,nsp,lmaxu,nl,idu,nqpp,s_pot%qpp,ldim,s_ham%iprmb,dmatu)
      endif

C ... Finish on-the-fly sampling integration of optics
      if (mod(loptic,10) == 1 .and.
     .    lteto == -1 .and. lwtkb >= 0 .and. procid == master) then
        npol = 3; if (loptic < 0) npol = 1
        dwopt = s_optic%dw; optrng = s_optic%window
        call dosmsh(-1,optrng,dwopt,0,ne(1),xx)
        allocate(wk3(ne(1)))
        call dosmsh(1,optrng,dwopt,ne,i,wk3); i = 13
        call optin2(s_optic,s_bz,i,nspx,nspc,s_lat%vol,npopt,
     .    npol,.false.,0,ne,wk3,xx,eps)
        deallocate(wk3)
        deallocate(eps)
      endif

C ... Case generating bands: find next block of qp
      if (plbnd /= 0) goto 99

C --- Interpolate density to Fermi energy (lmet=4) ---
C     sev = sumev(1,1)
      if (lmet == 4) then
        call rx('bndasa not ready for lmet=4')
      endif

C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet /= 4 .or. ltet > 0) then
        call bzwtsf(nbmax,ldim,nsp,nspc,n1,n2,n3,nkp,ntet,s_bz%idtet,
     .    zval,fsmom,isw(metal),ltet,mpsord,ndos,swidth,srnge,s_bz%wtkp,
     .    eband,efmax,-2,xx,efloc,sumev,xx,s_bz%wtkb,xx,xx(3),lwtkb)
C       call prmx('wtkb',s_bz%wtkb,ldim,ldim,nsp*nkp)
C       Double counting terms for global bfield
        eterms(21) = 0
        if (fsmom(2) /= 0) eterms(21) = fsmom(3)
        efermi = efloc(1)

        if (lmet /= 4) then
          ef0 = efermi
          s_bz%ef = ef0
        endif
        if (lwtkb == -1 .and. lrout > 0) then
          call info(20,0,0,' Start second band pass ...',0,0)
          lwtkb = 1
          goto 99
        endif
      endif

C ... Save k-point weights and Fermi level in wkp file
      if (lmet == 0) then
C        if (iprint() >= 20)
C     .  call awrit3(' Highest occ. level = %,5;5d '//
C     .  ' Lowest unocc. = %,5;5d  diff =  %,5;5d',' ',80,stdo,
C     .  evtop,ecbot,ecbot-evtop)
      endif
      if (procid == master) then
        ifi = fopna('wkp',-1,4)
        j = 0
        if (lmet == 0) j = 1
        i = iobzwt(j,ldimx,nkp,nspx,efermi,s_bz%wtkb,-ifi)
        call fclr('wkp',ifi)
      endif

C --- Generate density matrix ---
      if (ldens) then
C        call rx('bndasa: density matrix is not yet built')
C        complex(8), allocatable :: dmats(:)
C         allocate(dmats(nsite*nl**4)); call dpzero(dmats,2*nsite*nl**4)
C        call rx('bndasa: fix call to rsmsym')
CC        call rsmsym(1,plat,nbas,s_lat%pos,nl,nsp,1,nsite,w(ontab),w(oiax),
CC     .    s_lat%symgr,s_lat%ag,nsgrp,0,1,w(odmat),dmats)
C        call query('V>=40 to display rs s',-1,0)
C        call shostr(nl,nsite,nbas,1,plat,s_lat%pos,11000,w(oalph),
C     .    w(oiax),w(ontab),dmats,1,dum,1,1d0)
C        call rx0('done showing dmat')
      endif

C --- e'vector decomp'n and moments from supplied weights ---
      if (nevmx <= 0 .or. .not. metal .or. lwtkb /= 0) goto 1002
      if (efermi > efmax) call fexit2(-1,111,' Exit -1 BNDASA: '//
     .  'E_f = %1;6d exceeds efmax = %1;6d',efermi,efmax)
C ... Reread header in moments file, resetting internal variables
      if (liomom /= 1)
     .  call rx('BNDASA: no moments file ... can''t make output moms')
      call iomomq(nfilem,2,nl,nsp,nspc,nkp,ldim,0,i,nbmax,nlo*nclass,
     .  nrhos*nclass,nev,xx,xx,xx,xx,efermi,xx)
      i = nlo*nclass*nevx*mod(nfstg/10,10)
C     this should never happen
      if (mod(nfstg/10,10) /= (nspc+2)*nspc)
     .  call rx('BNDASA: moments file mismatch ... bug in bndasa')
      allocate(accwt(i)); call dpzero(accwt,i)
      i = 1
      if (lrhos) then
        i = nrhos*nclass*nevx*mod(nfstg/10,10)
        call dpzero(rhos,2*3*nrhos*2*2*nclass)
      endif
      if (allocated(accsm)) deallocate(accsm)
      allocate(accsm(i)); call dpzero(accsm,2*i)
      allocate(wk(nevx))
      swtkb = 0
      nschan = mod(nfstg/10,10)
      do  ikp = 1, nkp
      do  isp = 1, nsp
        if (onesp /= 0 .and. isp == 2) cycle
        if (nspc == 2 .and. isp == 2) cycle
        i = iomoms(nfilem,dum,nsp,nspc,1,ldim,nfstg,nschan,1,1,nevx,
     .    nevx,nlo*nclass,nrhos*nclass,nev,wk,accwt,accwt,
     .    accsm,dum,dum)
        i = isw(lrhos)
        call moment(i,nl,nlo,nrhos,isp,nsp,nspc,nkp,ldim,nevx,nev,ikp,
     .    s_bz%wtkb,nclass,accwt,accsm,qnu,rhos,orbtm,swtkb)
      enddo
      enddo
      deallocate(wk,accwt)
      if (lrhos) deallocate(accsm)
      if (dabs(zval-swtkb/nspc) > 1d-6) then
        call awrit3(' (warning) sum-of-weights (%;6d) does not'//
     .    ' match valence charge (%;6d)'//
     .    '%?#n<0#%N Consider increasing NEVMX##',
     .    ' ',160,stdo,swtkb,zval,nevmx-ldimx)
      endif
      if (onesp /= 0) call dscal(3*nl*nsp*nclass,2d0,qnu,1)
 1002 continue
      if (liomom /= 0 .and. procid == master) then
        nfilem = fopna('MOMS',-1,4)
        write (nfilem) efermi, vmtz(1)
        call fclose(nfilem)
      endif
      if (nevmx <= 0 .or. .not. metal .or. nspc /= 2) then
      elseif (cmdopt('--pdos',6,0,strn)) then
      else
        call mpibc2(orbtm,nlo*2*nclass,4,3,.false.,'bndasa','orbtm')
        call iorbtm(s_spec,0,0,s_ctrl%ics,nl,nlo,nclass,nsp,orbtm,s_bz%sopertp,xx)
C        if (n1*n2*n3 /= nkp)
C     .    call info0(20,0,0,' BNDASA warning! use --nosym to calculate orbital moments')
        if (n1*n2*n3 /= nkp)
     .    call logwarn(2,'%N BNDASA warning! use --nosym to calculate orbital moments')
      endif

C --- Restore potential parameters to orthogonal rep'n ---
      if (nmto == 0) call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp)

C --- Linear optics ---
  499 continue  ! Entry point for --opt:read
      if (cmdopt('--opt:read',10,0,strn) .or.
     .    cmdopt('--opt:write',11,0,strn)) then
        if (procid == master) then
          ifi = fopna('optdata',-1,4)
          rewind ifi
          if (cmdopt('--opt:read',10,0,strn)) then
            call info0(10,1,0,' reading optics data ...')
            read(ifi) efermi, i,j,k
            call sanrg(.true.,i,ldim,ldim,'bndasa:','ham dimsn')
            call sanrg(.true.,j,nsp,nsp,'bndasa:','nsp')
            call sanrg(.true.,k,nkp,nkp,'bndasa:','nsp')
            call dpdump(eband,ldim*nsp*nkp,ifi)
            read(ifi) i,j
            if (loptic < 0) then
              call sanrg(.true.,i,njdosw,njdosw,'bndasa:','njdosw')
              call sanrg(.true.,j,nfilm+nempm,nfilm+nempm,'bndasa:','nfilm+nempm')
              call dpdump(jdosw,njdosw*(nfilm+nempm)*nspx*nkp,ifi)
            elseif (loptic > 0) then
              call sanrg(.true.,i,nfilm,nfilm,'bndasa:','nfilm')
              call sanrg(.true.,j,nempm,nempm,'bndasa:','nempm')
              call dpdump(optmt,3*nfilm*nempm*nspx*nkp,ifi)
            endif
          else
            write(ifi) efermi, ldim,nsp,nkp
            call dpdump(eband,ldim*nsp*nkp,-ifi)
            if (loptic < 0) then
              write(ifi) njdosw,nfilm+nempm
              call dpdump(jdosw,njdosw*(nfilm+nempm)*nspx*nkp,-ifi)
            elseif (loptic > 0) then
              write(ifi) nfilm,nempm,nspx,nkp
              call dpdump(optmt,3*nfilm*nempm*nspx*nkp,-ifi)
            endif
          endif
        endif
        if (cmdopt('--opt:write',11,0,strn)) then
          call rx0('written optics data')
        endif
      endif

      if (nevmx > 0 .and. mod(iabs(loptic),10) /= 0 .and. lteto /= -1 .and.
     .    procid == master) then
        npol = 3
        if (loptic < 0) npol = 1
        if (loptic < 0 .and. njdosw > 0) npol = njdosw+1
        i = npopt*(1+nspx*npol)
        allocate(eps(i)); call dpzero(eps,i)
        if (lteto == 3) then
          call optinq(s_ctrl,s_optic,s_bz,eband,nbmaxx,nspx,nspc,efermi,
     .      s_bz%ipq,alat,plat,nfilm,nempm,npopt,npol,njdosw,optmt,jdosw,eps)
        else
          call optint(s_optic,s_bz,lteto,eband,nbmaxx,1,nspx,nspc,-1,efermi,
     .      s_bz%idtet,s_bz%wtkp,s_lat%vol,nfilm,nempm,npopt,npol,optmt,eps)
        endif
        deallocate(eps)
        if (cmdopt('--opt:read',10,0,strn)) call rx0('finished optics')
      endif

      if (bitand(loptic/10,2) /= 0) then
        call info0(30,0,0,' BNDASA: nonlinear optics')
        call intdrv(s_ctrl,s_optic,s_bz,eband,nbmax,nsp,nspc,efermi,
     .    s_bz%idtet,s_lat%vol,nfilm,nempm)
      endif

C --- Make local DOS for Stoner magnetism ---
C#ifdefC STONER
C       if (IAND(s_ctrl%lstonr,1) /= 0) then
C         stop 'upack zos,stoner index'
C        if (nlo /= nl)
C     .    call rx('BNDASA: (class,l,m) moments not allowed with STONER')
C        call rxx (nl < 3,'BNDASA: nl must be > 2 for Stoner')
C        if (iprint() >= 30) call awrit3(
C     .   ' BNDASA: Make number-of-states functions for Stoner model'//
C     .   ', %i bins in (%d,%d)',' ',100,stdo,ndos-1,dosw(1),
C     .   dosw(2))
C        i = 0
C        do  14  ic = 0, nclass-1
C        do  14  l = 0, nl-1
C        if (idxdn(l+nl*ic) /= 3) i = i+1
C   14   if (l == 2) index(ic) = i-1
CC       write (*,157) (ic, index(ic), ic = 0, nclass-1)
CC 157   format ('ic=',i3,' index=',i3)
C        if (liomom /= 0) call rx('BNDASA: no moments file')
C        nfilem = fopna('MOMS',-1,4)
C        allocate(wk(ndos))
C        nchan = nlo*nclass
C        if (ltet /= 1) then
C          real(8), allocatable :: eband(:)
C          real(8), allocatable :: doswt(:)
C          allocate(eband(ldim*nsp*nkp))
C          allocate(doswt(nchan*ldim*nsp*nkp)); call dpzero(doswt,nchan*ldim*nsp*nkp)
C          nfstg = 11
C          call iomomq(nfilem,22,nl,nsp,nspc,nkp,ldim,nfstg,i,ldim,
C     .      nchan,nchan,nev,eband,xx,doswt,xx,efermi,xx)
C          call rxx(i /= nkp,'BNDASA:  moments file missing qpts')
C          call dostet(ldim,nsp,nsp,nev,nchan,n1,n2,n3,ntet,
C     .      s_bz%idtet,eband,doswt,ndos,dosw(1),
C     .      dosw(2),.true.,wk,zos)
C        else
C          allocate(eband(ldim))
C          allocate(doswt(nchan*ldim))
C          nfstg = 0
C          call iomomq(nfilem,12,nl,nsp,nspc,nkp,ldim,nfstg,i,ldim,
C     .      nchan,nchan,nev,xx,xx,xx,xx,efermi,xx)
C          call rxx(i /= nkp,'BNDASA:  moments file missing qpts')
C          call dosspl(nfilem,ldim,nsp,1,nchan,mpsord,swidth,nkp,
C     .      s_bz%wtkp,eband,doswt,ndos,dosw(1),dosw(2),.true.,
C     .      wk,zos)
C        endif
C        deallocate(wk,eband,doswt)
CC   ... Scale zos to get per spin
C        call dscal(ndos*nchan,0.5d0,zos,1)
C      endif
C#else
      if (IAND(s_ctrl%lstonr(1),1) /= 0)
     .  call rx('BNDASA: recompile with ccomp -dSTONER ...')
C#endif

C --- Generate DOS on disk ---
      if (IAND(s_ctrl%ldos,1) /= 0 .and. procid == master) then
        allocate(dos(3*ndos))
        if (iprint() >= 30) call awrit1('%x%N ... Generating %?#n<0#'
     .    //'integrated#total# DOS',' ',80,stdo,ndos0)
        if (ltet > 0) then
          call bzints(n1,n2,n3,eband,dum,nkp,ldimx,nbmaxx,nspx,
     .      dosw(1),dosw(2),dos,ndos,efermi,1,ntet,s_bz%idtet,
     .      dum,xx)
          if (ndos0 > 0)
     .      call xxxdif(dosw(1),dosw(2),ndos,nspx,0,dos)
C         del = 0d0
        else
          if (mpsord >= 100) mpsord = mod(mpsord,100)
          if (ndos0 > 0)
     .      call makdos(nkp,ldimx,nbmaxx,nspx,s_bz%wtkp,eband,mpsord,
     .      swidth,-srnge,dosw(1),dosw(2),ndos,dos)
          if (ndos0 < 0)
     .      call maknos(nkp,ldimx,nbmaxx,nspx,s_bz%wtkp,eband,mpsord,
     .      swidth,-srnge,dosw(1),dosw(2),ndos,dos)
C         del = mpsord+swidth
        endif
        if (nspc == 2) call dscal(ndos,.5d0,dos,1)
C        call dosio(dos,ndos,nspx,ndos,1,dosw(1),dosw(2),nspx,
C     .    efermi,del,1,-fopn('DOS'))
        call iodos(3,-fopn('DOS'),dos,ndos,nspx,ndos,1,dosw(1),
     .    dosw(2),nspx,efermi,1)
        call fclose(fopn('DOS'))
        deallocate(dos)
      endif

C --- Magnetic moments corresponding to (output) density matrix ---
      if (lrhos) then
        if (iprt(2) /= iprt(1)) then
          call pshpr(iprt(2))
        else
          call pshpr(iprint()-0)
        endif

        allocate(bxcnw(3*nclass)); call dpzero(bxcnw,3*nclass)
        call amagnc(nbas,nl,s_ctrl%ipc,rhos,nrhos,qnu,s_ham%eula,neul,
     .    2,amag,aamom,bxcnw)
        if (bittst(lncol,8)) then
          call bdotsr(nbas,nl,s_ctrl%ipc,rhos,nrhos,xx,
     .      s_ham%bdots,lihdim,s_ham%iprmb,0,dum)
C         ehterm(6) = dum(1)
          eterms(18) = dum(1)
        endif
        if (neul > 1 .or. bittst(lncol,8)) then
          if (.not. bittst(lncol,8)) nbf = 99
          if (.not. associated(s_ham%magf)) s_ham%magf => xx
          call bsrhos(nbas,nl,s_ctrl%ipc,rhos,nrhos,xx,s_pot%pp,
     .      s_pot%sop,s_ham%eula,neul,s_pot%bxc,s_ham%magf,nbf,lihdim,
     .      s_ham%iprmb,0,dum)
          eterms(18) = dum(2)
        endif
        call poppr

        call amagn2(nl,nclass,nbas,s_ctrl%ipc,s_ham%eula,neul,
     .    bxcnw,qnu,nrhos,rhos)

C       call dcopy(3*nclass,bxcnw,1,s_pot%bxc,1)
        deallocate(bxcnw)

      endif

C --- Make Magnetization of r
C      call asvsph2(
C     .  s_ctrl,s_lat,s_spec,s_ham,s_pot,
C     .  rhos,max(neul,nl),1,0,1)
C ... To check quantities integrated throughout the sphere:
c      call asvsphint(rhos,
c     .  max(neul,nl),1,0,1)

C     call info0(30,0,0,' ')

C ... Repack d.c. terms from applied field
      s_ham%eterms = eterms

C ... Restore struc variables we changed
      s_ctrl%lham = lhams
      s_ctrl%lasa = lasas
      s_ham%lham = lhams

      deallocate(strx,hrs,ors,chan,idu,ifbls,orbtm,jdosw)
      i = idalloc('optmt',allocvb()+3,0,0)
      i = idalloc('wtkb',allocvb()+3,0,0)
      i = idalloc('swtk',allocvb()+3,0,0)

      call tcx('bndasa')
      end

C      subroutine gethkp(nsp,isp,nkp,ikp,ldim,hk,hkp)
C
C      integer nsp,isp,nkp,ikp,i,j
C      double precision hk(ldim,ldim*2)
C      double precision hkp(ldim*2,ldim)
C
C      call dpzero(hkp,ldim*ldim*2)
C      do i = 1, ldim
C        do j = 1, ldim
C
C          hkp(i*2-1,j) = hk(i,j)
C          hkp(i*2,j) = hk(i,j+ldim)
C        end do
C      end do
C
C      end
C      subroutine snot(procid,optmt,n,nkp)
C      integer procid
C      double precision optmt(n,nkp)
C
C      print *, procid, optmt(1,:)
C      end

