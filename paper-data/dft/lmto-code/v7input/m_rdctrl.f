      module m_rdctrl
C- Module to read input for LMTO package
C ----------------------------------------------------------------------
Cr This is the module that parses an input file for data; it is used
Cr by almost all the main programs in this package.
Cr
Cr *Routine readctrl(prgn,..) is the main entry point.
Cr  What tokens are sought depends on the argument 'prgn', which
Cr  specifies which main program is calling m_rdctrl.
Cr
Cr *Rules governing whether tag is to be parsed and in what manner:
Cr  Routine toksw_init (in this module) generates a list of tokens to be
Cr  parsed. There is a separate list for each main program.
Cr
Cr  Routines tksw and tksw2 (module m_toksw) use this list to return
Cr  a switch indicating whether and how the tag is to be parsed.
Cr  A tag may be required, optional, excluded, or ignored
Cr  and information is given whether it is enabled fro express mode.
Cr
Cr  Routine 'tkadd' (in m_toksw) provides the linkage between toksw_init
Cr  module m_toksw.
Cr
Cr *How tags are parsed:
Cr   The family of routines defined by interface gtv, i.e.
Cr     gtv_r8,gtv_r8v,gtv_i4,gtv_i4v,gtv_lg,gtv_char,gtv_none
Cr   (module _gtv) do the parsing for content of tag.
Cr   gtv requires, at a minimum, a tag name you supply and instructions
Cr   governing how it is to be parsed (this is returned by tksw2).
Cr
Cr   gtv has a great deal of flexibility in whether or how the contents
Cr   following a token or tag are to be parsed.
Cr   See gtv documentation, particularly after gtv_entrance().
Cr
Cr  Routine readctrlpq reads class-specific data, for ASA-based programs.
Cr
Cr *To add a new tag for input, select a category and token name.
Cr    1: include it in the 'tag' list for each main program that calls it;
Cr       see toksw_init below.
Cr    2. add a call to gtv somewhere in the appropriate location in the
Cr       category you have chosen.  subroutine readctrl supplies examples.
Cr    3. If you want to enable this tag for 'express' mode, append
Cr       a '%' or '&' to the tag.
Cr
Cr *Not yet implemented:
Cr  We could have express mode not skip the tree strructure, and use
Cr  the last part of the tag only.
Cr
Cu Updates
Cu   20 Apr 16 New express mode
C ----------------------------------------------------------------------
      implicit none
C ... Fixed parameters
      integer,parameter:: NULLI = -99999
      real(8),parameter::    NULLR = -99999
      logical,parameter:: T=.true.,F=.false.
      integer,parameter :: n0=10
      integer,parameter :: nkap0=4
      integer(2):: nono

C ... IO
      integer:: io_show=0,io_help=0,nvario=0,maxmem

C ... HEADER, SYMGRP
      character(256):: header,symg=' '

C ... HAM
      logical :: frzwf=.false.,ham_ewald=.false.
      integer:: ctrl_lfrce=0,lxcf=2,gga=0,ftmesh(3),nmto=0,lrsig=0,nsp=1,ham_qasa,
     .  lrel=1,lso=0,ham_udiag=0,autob=0,pnudef=1,ham_rdvext=0
C     integer:: lmxb_cut(2)=nulli
      logical:: ltbe  ! set to T if TBE program (used for defaults)
      integer:: lfp   ! set nonzero if FP program (used for defaults)
      real(8):: lat_gmax,tolft,dqval,kmto(10),rsrnge,vmtz,elind=-0.7d0
      real(8):: alfsi=0d0,dabc(3)=nullr,rsstol,v0beta=1d0
      real(8):: pmin(n0),pmax(n0),basopt(15)
C     sigp holds parameters for approximating self-energy sigma
C       arg 1: mode : specifies how to set its diagonal part
C              for states above the high-energy cutoff
C              0 constrain sigii to be > asig+bsig*e
C              1 constrain sigii to be = asig+bsig*e
C              2 constrain sigii to be > asig and < bsig
C              3 constraint same as mode 1.
C                Mode 3 differs in that the least-squares fit to
C                sigii (for informational purposes only, to help
C                estimate asig and bsig) is done for states between
C                efit and nmax or emax
C              4 fix sigii to eseavr given in sigm file.  No fitting
C       arg 2: nmin : sigma for states 1..nmin are approximated by sigii
C       arg 3: emin : (used only if nmin<0)
C                   : sigma for levels e<emin are approximated by sigii
C       arg 4: nmax : sigma for levels i>nmax are approximated by
C                     sigii AND constrained according to mode
C       arg 5: emax : (used only if nmax<=0)
C                   : sigma for levels e<emax are approximated by
C                     sigii AND constrained according to mode
C       arg 6: asig : constraint used to approximate
C                     sigii = asig + E * bsig  or
C                     asig < sigii < bsig
C       arg 7: bsig : constraint used to approximate
C                     sigii = asig + E * bsig  or
C                     asig < sigii < bsig
      real(8):: sigp(10)=0,sigp_emin,sigp_emax,sigp_a,sigp_b,sigp_efit
      integer:: sigp_mode=3,sigp_nmin=0,sigp_nmax=0
      equivalence (sigp_emin,sigp(3)),(sigp_emax,sigp(5)),
     .             (sigp_a,sigp(6)),(sigp_b,sigp(7)),(sigp_efit,sigp(8))

C ... OPTIONS
      logical ::
     .  lcd1,lcd2,lcd4=F,lcd8=F,lcd64,lasa4=T,lasa8=F,lasa32=F,lasa64=F,
     .  lasa512,lasa1024,lham1=F,lham4=F,lham8=F,lham16=F,lham32=F,lham64=F,lham128=F,
     .  lham256=F,lham512=F,lves=F,lcd16=F,lfp2=F,lvshft=F
C     Initialize, since read is only executed when NSPIN=2
      logical:: lncol1=F,lncol2=F,lncol4=F,lncol8=F,lncol16=F,
     .  lncol32=F,lncol64=F,lncol128=F,lBz=F,lsos=F,ctrl_lbxc4=F,ctrl_lbxc16=F
      character(256):: sxopt=' '
      integer:: lasa3=1,lsx=0,lscr=0,lrquad=0,smalit(2),nesabc(3),
     .  lstonr(3)=0,quit=0,nl=0,lham3=0,lekkl=1,ctrl_lbxc1=0,lrsa=0
      real(8):: asa_elin,rmines,rmaxes,ham_qss(4)

C ... STRUC
C     Initializing nbas and alat flags no nbas or alat has been input
      real(8):: lat_slat(9),dlat,alat=NULLR,plat(9),lat_gam(4),dalat
      real(8):: vol,avw,lat_dist(9)
      integer:: nbas=NULLI,nbasp=NULLI,nsite=NULLI,nspec=NULLI,nclass
      integer:: lat_ldist=0

C ... SPEC
      real(8):: omax1(3),omax2(3),wsrmax,sclwsr=0
      integer,parameter ::mxspec=256
      character*8, target ::  slabl(mxspec), ch*3
      integer:: lmxax,nkaph
      logical,allocatable:: mxcst1(:),mxcst2(:),mxcst4(:)
      integer,allocatable:: idxdn(:,:), grp(:),grp2(:),
     .  idu(:,:),lmxb(:),lmxa(:),idmod(:,:),iq1(:,:),iq2(:,:),
     .  kmxt(:),kmxv(:),lfoca(:),lmxl(:),lxi(:),nxi(:),nr(:),lmxpb(:)
      real(8):: hsfitd(2),hfacd(3)=0
      real(8),allocatable:: rsmh(:,:),rsmh2(:,:),eh(:,:),eh2(:,:),
     .  hcr(:,:),rsminl(:,:),rs3(:),rham(:),alpha(:,:),ehvl(:,:),
     .  dv(:), uh(:,:),jh(:,:),
     .  qpol(:,:),stni(:),tbvso(:,:),
     .  pnu(:,:,:),qnu(:,:,:),
     .  coreq(:,:),mass(:),colxbs(:,:),radxbs(:),hsclsm(:),hsfitp(:,:),
     .  rg(:),rsma(:),rfoca(:),
     .  rsmfa(:),rcfa(:,:),
     .  exi(:,:),rint(:),rcut(:),rmt(:),pz(:,:,:),
     .  amom(:,:),spec_a(:),z(:),eref(:),rsmv(:)
      character*(8),allocatable::pb1(:),pb2(:), coreh(:)
      integer,allocatable:: nthet(:),nang(:,:)
      integer,allocatable:: ncpa(:),nbeff(:)
      integer,allocatable:: iscpa(:,:)
      real(8),allocatable   :: xcpa(:,:),beff(:,:),hfacs(:,:)

C ... SITE
      character*8:: alabl
      real(8),allocatable   :: pos(:,:),vel(:,:),eula(:,:),vshft(:)
      integer,allocatable:: ips(:),ipl(:),plv(:),irlx(:,:),ndelta(:)
      real(8),allocatable   :: delta(:,:),mpole(:),dpole(:,:),rmaxs(:)

C ... Iterations
      character(128) :: iter_mix=' ',iter_amix=' '
      real(8):: smix(37)=0,ctrl_tol(5)=1d-4
      integer:: iter_maxit=1

C ... BZ
      integer:: bz_nabc(3)=NULLI,bz_lshft(3)=0,
     .  bz_lmet,bz_n=0,bz_nevmx,bz_lmull=0,ctrl_ldos=0,bz_ndos,ctrl_lmet2
      real(8):: bz_w,bz_ef,bz_def,bz_efmax,bz_zval,bz_fsmom(2)=NULLR,
     .  bz_semsh(10),zbak(2)=0d0,bz_dosw(2),bz_lcond(4)=0,bz_range=5d0,bz_zinsul
      integer bz_lio48
      logical:: bz_lio1=.false.,bz_lio2=.false.
     .  ,ctrl_lmet4=.true.,ctrl_lmet8=.true.,ctrl_lqp1=.false.,ctrl_lqp2=.false.
      character(256) :: strn_syml=' '

C ... Fit
      integer:: fit_mode=NULLI,fit_nbfit(2)=0,fit_nbfitf,
     .  fit_shft(2)=0,fit_ndos=NULLI
      real(8):: fit_ebfit(2)=0d0,fit_rmfit(2)=0d0,fit_alam=NULLR,
     .  fit_alsc=NULLR,fit_wt(3)=[0.0d0,0.0d0,0.1d0],fit_gausw(2)

C ... TB
      logical:: ltb1   =.false., ltb2   =.false., ltb4   =.false.
     .         ,ltb8   =.false., ltb16  =.false., ltb32  =.false.
     .         ,ltb64  =.false., ltb128 =.false., ltb256 =.false.
     .         ,ltb512 =.false., ltb1024=.false., ltb2048=.false.
     .         ,ltb4096=.false., ltb213 =.false., ltb214 =.false.
     .         ,ltb215 =.false., ltb216 =.false., ltb217 =.false.
     .         ,ltb218 =.false.

C ... Ewald
      real(8):: lat_as,lat_tol,lat_rpad=0
      integer:: lat_nkdmx

C ... STR
      real(8):: str_rmax=nullr,str_rmaxg=nullr,str_kaps(6)=NULLR,
     .  str_drwats=NULLR,iinv_tol,str_delrx,str_tolg,str_rsminl(n0)=NULLR
      integer:: str_mode=nulli,str_nkaps=1,str_mxnbr,str_lmaxw=NULLI
      integer:: iinv_nit,iinv_ncut,tcf_nalf,str_ivl=0,str_mnn
      integer:: str_mem=1
      logical:: str_lshow1,str_lequiv1

C ... TCF
      integer tcf_nbisi(3),tcf_ncupl,tcf_ndust
      real(8):: tcf_adec,tcf_wztcf

C ... CGF
      integer:: ctrl_lcgf = 0
      real(8)   :: ctrl_tdlm = 0
      integer:: ctrl_ldlm = 0

C ... PGF
      integer:: ctrl_lpgf(2)=0
      character(256):: gfopt=' '
      real(8):: platl(3)=0d0,platr(3)=0d0

C ... DYN
C   structure of mdprm:
C   arg 1: 0 no relaxation or dynamics
C          1 new dynamics 2  restart dynamics
C          4 relax with conjugate gradients
C          5 relax with variable metric
C          6 relax with Broyden
C   arg 2: statics: switch
C          1 read hessian matrix
C          dynamics:
C            number of iterations between printouts.
C   arg 3: (stat) relaxation x-tolerance
C          (dyn)  temperature
C   arg 4: (stat) relaxation g-tolerance
C          (dyn)  time step
C   arg 5: (stat) step length
C          (dyn)  relaxation time
C   arg 6: (stat) Remove hessian after this many steps
C          (dyn)  --
C   Structure of sdprm :parameters governing dynamics (see magtrq.f)
C   arg 1: scale factor amplifying magnetic forces
C   arg 2: time step tau for Landau dynamics
C   arg 3: reference energy etot0 for Nose thermostat
C          (not used unless fscal is zero)
C   arg 4: maximum allowed change in angle (not used unless nonneg)
C   arg 5: etol: set tau=0 this iter if etot-ehf>etol

C  prmint  Parameters for numerical integration
C          For Bulirsch-Stoer integration:
C  arg  1:   mode: (1 for BS-integration)
C  arg  2:   rst:  1 for start new integration, 0 to continue
C  arg  3:   ts0:  minimum time step size
C  arg  4:   tol:  tolerance in integration errors
C  arg  5:   mx:   order of rational function extrapolation
C  arg  6:   mi:   number of midpoint rules
C  arg  7-17 nseq: sequence of no. midpoint divisions
C  arg  18:        offset to bs workspace
      integer:: sdmod=-1
      real(8):: mdprm(7)=0,lat_defm(6),sdprm(5)
      real(8):: prmint(20)=(/1d0,1d0,0d0,2d-4,7d0,11d0,2d0,4d0,6d0,8d0,
     .  12d0,16d0,24d0,32d0,48d0,64d0,96d0,0d0,0d0,0d0/),
     .  prmint_ts0,prmint_tol,gd_ct(3),
     .  move_kt,move_ts,move_tsequ,move_tstot=NULLR
      logical:: prmint_new,lbsprm=F,lsdyn=F
      integer:: nitmv,prmint_mi,prmint_mx,prmint_nseq(11),
     .  gd_modt(3),gd_nmodt=1
      character(1024):: mmham=' '
      real(8), parameter:: fs = 20.67098d0, degK = 6.3333d-6   ! defaults for MD
      real(8), parameter:: kboltzeV=8.6173303d-5  ! Boltzmann contant in eV
      real(8), parameter:: RyeV=13.60569193d0 ! 1 Ry, in eV

C ... DMFT
      integer ncix   ! Total number of correlated subblocks among all sites
      integer nicix  ! Number of inequivalent correlated subblocks
      integer dmft_ncatom                 ! Equivalent to s_dmft%ncatom
      integer dmft_maxl                   ! Maximum L among all DMFT blocks, for dimensioning
      integer dmft_maxdim                 ! 2*dmft_maxl+1
      integer dmft_nzsig                  ! Total number of nonzero matrix elements for 1 spin, all sites
      logical:: dmft_imw = .true.         ! DMFT on imaginary axis
!     integer,allocatable:: cixmap(:,:)   ! cixmap(icix,1) points to a cix block this inequivalent block belongs to
      integer,allocatable:: dmft_l(:)     ! (1:nicix) Same as s_dmft_l
      integer,allocatable:: dmft_qsplit(:)! (1:nicix) Same as s_dmft_qsplit
      integer,allocatable:: dmft_umode(:) ! (1:nicix) Same as s_dmft_umode
      integer,allocatable:: dmft_ib(:)    ! (1:ncix) Same as s_dmft_ib
      integer,allocatable:: dmft_icix(:)  ! (1:ncix) Same as s_dmft_icix
      integer,allocatable:: dmft_sigind(:,:,:,:) ! (1:nicix) Will become s_dmft_sigind
C     integer,allocatable:: dmft_ndim(:)   ! number of channels in each independent cix block
      integer,allocatable:: dmft_ndsigi(:) ! accumulated number of channels in each independent cix block

!      character(8), allocatable :: cixlabl(:)

      integer :: dmft_nabc(3)=NULLI,dmft_nlohi(2)=NULLI,dmft_proj=NULLI,dmft_knorm,dmft_nomg
      integer :: dmft_nomf,dmft_nomb
      real(8):: dmft_wlohi(2)=NULLR,dmft_beta,dmft_broad
C ... GW
      integer:: gw_mksig=0,gw_nabc(3)=NULLI,gw_lshft(3),
     .  gw_nband=-NULLI,gw_nime=6,gw_code=0,gw_usebsew=0
      real(8):: gw_gcutb=NULLR,gw_gcutx=NULLR,gw_ecuts=NULLR,
     .  gw_delre(2)=NULLR,gw_deltax=NULLR,gw_deltaw=NULLR,
     .  gw_gsmear=NULLR,gw_pbtol(8)=NULLR,gw_qoffp=NULLR,autoqcut=0,gw_ecutpb=NULLR

C ... PLANE
      real(8):: slat_plat2(9)

C ... ORDERN
C      integer lordn
C      real(8):: plate(9)    !slat, 'lat plate'
C      real(8):: efre,ordn_rmaxg

C ... OPTICS
      logical:: optic_alltrans=.false.
      integer:: optic_ne,optic_ocrng(2),optic_unrng(2),ctrl_loptc=0,
     .  optic_nchi2=0,optic_axes(18),optic_ltet,optic_mode1,optic_iq(3),
     .  optic_mefac=0,optic_ffmt=0,optic_kee,optic_nmp
      real(8):: optic_window(3)=0d0,optic_esciss,optic_dw(3)=0d0,optic_imref(2)=NULLR,
     .  optic_esmr,optic_kt,optic_w

C ... APW basis
      integer:: pwmode=0,npwpad,ovncut=0
      real(8)::    pwemax,pwemin,oveps=0d0

      contains

      subroutine readctrl(prgnam,vstrn,ivn)
C- Reads data parsed from ctrl file into m_rdctrl
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   vstrn :version
Ci   ivn   :
Co Outputs
Cl Local variables
Cr Remarks
Cu Updates
Cu   26 Mar 19 Start on reading complete DMFT input from ctrl file
Cu   29 Jun 17 EH,EH2 read 1:lmxb data regardless of how many read into RSMH
Cu   25 Jun 17 Change default for GW_ECUTS to depend on HAM_SIGP_EMAX; GW_PBTOL now a vector
Cu   20 May 17 Change default rsrnge to scale with lattice volume^1/3
Cu   25 Oct 16 Change default values for ITER_CONV from 1d-4->1d-5 and ITER_CONVC from 1d-4->3d-5
Cu   12 Oct 16 New OPTICS_FERMI
Cu   20 Apr 16 New express mode
Cu             The following tag default values were modified:
Cu             IO_VERBOS=35    (was 30)
Cu             HAM_ELIND=-0.7  (was 0)
Cu             HAM_SIGP_MODE=4 (was 3)
Cu             HAM_FORCES=1, or 0 if SO is turned on (was 0)
Cu   15 Oct 15 New tags for DMFT
Cu   03 Aug 15 New ability to read NSPEC and species labels from site file
Cu             Extended switches to PFLOAT
Cu   15 Dec 14 New OPTICS_IMREF OPTICS_KT
Cu   27 Aug 14 New IO_MAXMEM, IO_SHOMEM
Cu   10 Mar 14 Allow args 2,3 of BZ_NKABC to be negative
Cu   30 Sep 13 New SPEC_HFAC_V
Cu   09 Sep 13 Handles missing SPEC_ATOM_LMX better
Cu   09 Aug 13 New SPEC_NANG, and SO=11 option
Cu   17 Jun 13 New SPEC_HFAC
Cu   13 May 13 New definition of ivn
Cu   17 Apr 13 New switch BZ_DMATK
Cu   28 Jan 13 New switch OPTIONS_NMCORE OPTIONS_SHORBZ
Cu   31 Aug 12 New input GW_BZJOB
Cu   25 Apr 12 Parameters for chemical CPA
Cu   22 Apr 12 New OPTIONS_RQUAD
Cu   22 Feb 12 Parameters for CPA, disordered local moments
Cu   16 Apr 10 New ITER_FIT_W
Cu   01 Apr 10 Bfield can have just z component; enable lmf to read it
Cu   27 Mar 10 New parameters for Levenberg-Marquardt fitting
Cu   02 Feb 10 command-line verbosity takes precedence for entire stack
Cu   09 Sep 09 New OTPICS inputs
Cu   11 Aug 09 Inputs for Levenberg-Marquardt fitting of ASA bands
Cu   07 Jul 09 New HAM_AUTOBAS_ELOC, GW_CODE
Cu   10 May 09 New input for automatic parameter generation
Cu   24 Apr 09 New str->delrx,mnnbr,tolg and site->rmaxs
Cu   01 May 08 Added DLM vars ldlm,tdlm; arrays opnud,oqnud,oppd (Kirill)
Cu   09 Aug 07 First created
C ----------------------------------------------------------------------
      use m_toksw
      use m_gtv
      implicit none
C ... Passed parameters
      character(len=*) :: prgnam
      character(6):: vstrn(2)
      integer ivn(3,2)
C ... Local allocatable arrays
      character*8, pointer :: slabll(:)
      integer, allocatable :: ipsl(:)
C ... Local variables
      integer :: nout,i0,i,j,jp,k,nn,ivec(49)
      integer :: io_iactive,io_tim(2),verbos(5),iapflt(2)=[1,1],llxcf(3)
      character(256)::  a,outs
      logical :: debug = .false.      ! Set to .true. to turn on debugging printout
      logical :: ldefaultl            ! Set to .true. if lmx not read explicitly
      logical :: linitsite = .false.  ! Set to .true. after site positions are read
      logical:: lread_site = .false.  ! set to .true. if basis info read from site file
      integer,parameter :: kmxvdefault=15
      integer,parameter :: master=0
      integer :: procid,nproc
      integer :: lmxbj=0,lmxaj=0,nlbj,nlbji,nlaj,nlaji
      integer :: mpipid,vsn
      double precision vers,xv(2*n0)
      character(256*64) :: bigstr=' '
      integer :: it(8)
      logical,parameter:: T=.true.,F=.false.,mlog=.false.

      integer :: lp1,lpzi,nphimx,nlmax,nrmx=5001
      real(8)::  xx,samix(37)
      integer :: izerv(n0)=(/(0,i=1,n0)/)
      real(8)::  zerov(n0)=(/(0d0,i=1,n0)/)

      real(8):: tbalph(6,6)
      integer :: ii,sw,stdo,stdl,stde
      character(128) :: nm
      character(len=10) :: prgn(2) !,swk(0:3)

      real(8):: nullrv(256),qlat(3,3)
      integer :: nulliv(256),jj(2),nkapsi
      integer(2) :: nono
      logical:: ltmp

      procedure(logical) cmdopt,parmxp,isanrg
      procedure(integer) iprint,lgunit,isw,iosite,fopna,a2vec,getdig,dmft_equiv_cix,wordsw
      procedure(real(8)) dasum,dglob,avwsr

      data tbalph/.214286D0,.000001D0,.000001D0,.000001D0,.001D-3,0d0,
     .            .287234D0,.025822D0,.000001D0,.000001D0,.001D-3,0d0,
     .            .348485D0,.053030D0,.010714D0,.000001D0,.001D-3,0d0,
     .            .385057D0,.073209D0,.022481D0,.006069D0,.001D-3,0d0,
     .            .404761D0,.087927D0,.032886D0,.012257D0,4.2735D-3,0d0,
     .            .404761D0,.087927D0,.032886D0,.012257D0,4.2735D-3,.001D-3/

      pmin = 0
      dmft_beta = 0
      io_iactive = 0
      vmtz = 0
      ncix = 0

C ... For MPI
      procid = mpipid(1)
      nproc  = mpipid(0)
C     pi = 4*datan(1d0)
C     mlog = cmdopt('--mlog',6,0,a)
      ltmp = .false.   ! So the compiler doesn't complain
      llxcf(1) = 2
      a = ''; express_cat = ''
      slabll => slabl  ! Unless read by iosite, slabll will correspond to slabl
      ctrl_tol(4) = 0

      prgn(1) = prgnam; prgn(2) = prgnam
      if (prgnam == 'LMFDMFT' .or. prgnam == 'LMFGWD' .or. prgnam == 'LMFGW' .or.
     .    prgnam == 'LMFGWS' .or. prgnam == 'LMFA') prgn(2) = 'LMF'
      if (prgnam == 'LMPG') prgn(2) = 'LMGF'
      if (prgnam == 'LMSTR' .or. prgnam == 'LMCHK' .or. prgnam == 'LMXBS' .or. prgnam == 'LMPLAN' .or.
     .    prgnam == 'LMCTL' .or. prgnam == 'LMSCELL' .or. prgnam == 'TBE' .or. prgnam == 'MMAG' .or.
     .    prgnam == 'LMDOS' .or. prgnam == 'LMMC') prgn(2) = 'BASE'

C     Used to set tbe-dependent or f-dependent defaults
      ltbe = prgnam == 'TBE'
      lfp = 0; if (prgn(2) == 'LMF' .or. prgnam=='LMMC') lfp = 1
C     Flags representing how cd is represented: default values
      lcd4 = (prgn(2) == 'LMF' .or. prgnam=='LMMC')
      lcd8 = (prgnam == 'LMMC')

C --- Initialize ---
      nullrv = nullr
      nulliv  =nulli
      bz_lmet = 1; if (lfp /= 0) bz_lmet = 5
      debug = cmdopt('--debug',6,0,a)
      call toksw_init(debug)
      stdo = lgunit(1); stdl = lgunit(2); stde = lgunit(1)
C     call gtv_setst(debug,stdo,stdl,stde)
      ltmp = cmdopt('--show',6,0,a)
      if (ltmp .and. scan(a(7:7),' p=')>0) io_show = 1
      if (ltmp .and. a(7:7)=='=') then
        i = 7
        i = a2vec(a,len(a),i,2,' ',1,1,1,it,io_show)
      endif

      if (procid /= master) io_show = 0
      if (cmdopt('--input',7,0,a)) io_help = 1
      if (cmdopt('--tags',6,0,a)) io_help = 2
      call gtv_setio(debug,io_show,io_help)
      if (io_help /= 0) then
        write(stdo,332)
  332   format(/' Tag',t25,'Input   cast  (size,min)'/
     .    ' ------------------------------------------')
      elseif (io_show /= 0) then
        write(stdo,333)
  333   format(/' Tag',t25,'Input   cast  (size,min,read,def)     result')
      endif

C --- Initial IO ---
C      nm='IO_SHOW'; call gtv(trim(nm),tksw2(prgn,nm),io_show,def_i4=
C     .  0,note='Echo data as it is read from input file')
C      if (cmdopt('--show',6,0,a)) io_show = 1
C      if (procid /= master) io_show = 0
C      if (cmdopt('--show=',7,0,a)) then
C        i = 7
C        i = a2vec(a,len(a),i,2,' ',1,1,1,it,io_show)
C      endif
C      call gtv_setio(debug,io_show,io_help) ! In case io_show changed

C      nm='IO_HELP'; call gtv(trim(nm),tksw2(prgn,nm),io_help,
C     .  def_i4=0, note='Show what input would be sought, '//
C     .  'without attempting to read data')
      if (cmdopt('--input',7,0,a)) io_help = 1  !optio=0 in old code
      if (io_help == 1) io_show = 1
      call gtv_setio(debug,io_show,io_help) ! In case io_help changed

C --- Express ---
   48 continue
      if (tksw2(prgn,'EXPRESS') >= 2) goto 49
      if (express == 0) then
        nm='EXPRESS'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,Texist=ltmp,
     .    note='$1 Express input for commonly used tags')
      else
        ltmp = .true.
      endif
      if (ltmp) then
        express_cat = 'EXPRESS'
        if (express == 0) then
          express = -1
C         call gtv_express(express)
        else
          express = 1
C         call gtv_express(express)
        endif
      endif
C ... End of LMF tokens
   49 continue

C --- Version ---
      if (express /= -1) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Version control ---')
      do  j = 1, 2
      vsn = ivn(1,j)
      if (vsn /= 0) then
        nm = 'VERS_'//vstrn(j)
        outs = 'Input style'; if (j == 2) outs = 'Program'
        sw = tksw2(prgn,nm)
        if (sw < 2) then
          call gtv(trim(nm),sw,vers,note=trim(outs)//
     .    ' version check: integer part must be at least as large as:')
          if (io_help >= 1) then
            call info2(2,0,0,'   '//trim(outs)//' version: %i',vsn,0)
          else
            if (io_show /= 0) call info2(2,0,0,'  Expected '//
     .        trim(outs)//' version: %i ... found in file: %d',vsn,vers)
            if (vers < vsn) call fexit(-1,1,
     .        '  Exit -1 RDCTRL: ctrl file requires '//trim(nm)//' >= %i',vsn)
            if (int(vers) > vsn) call fexit(-1,1,
     .        '  Exit -1 RDCTRL: ctrl file newer than v%i',vsn)
          endif
        endif
      endif
      enddo
      endif

C --- CMD ---
C      nm='CMD'; call gtv(trim(nm),tksw2(prgn,nm), bigstr, note=
C     .  'Contents are appended to command-line arguments',nout=nout)
C      if (nout == 1) then
C        call acmdop(bigstr,len_trim(bigstr),0)
C      endif

C --- CONST ---
      nm='CONST'; if (io_show+io_help/=0 .and. tkswp(prgn,express,nm) < 2)
     .  call info0(2,1,0,' --- Variable assignments ---')
      call numsyv(nvario)
      call gtv(trim(nm),tksw2(prgn,nm),bigstr,nout=nout,note=
     .  'Declare variables for use in expressions.'//
     . '%N%3fYou may also declare variables on the command-line:  -vnam=#')
      if (nout == 1) then
        i = 0
        call parsyv(bigstr,len_trim(bigstr),1999,0,i)
        if (io_show /= 0 .and. iprint() >= 40) call shosyv(0,0,0,stdo)
      endif

C --- IO ---
      nm='IO_SHOW'; if (io_show+io_help/=0 .and. tkswp(prgn,express,nm) < 2)
     .  call info0(2,1,0,' --- I/O management ---')
      if (io_show == 0 .and. express < 1) then
        call gtv(trim(nm),tksw2(prgn,nm),io_show,def_i4=
     .    0,note='Echo data as it is read from input file')
        call gtv_setio(debug,io_show,io_help) ! In case io_show changed
      endif
      nm='IO_VERBOS'; call gtv(trim(nm),tksw2(prgn,nm),verbos,
     .  note='Verbosity stack for printout.'//
     .  '%N   May also be set from the command-line: --pr#1[,#2]',
     .  def_i4v=(/35/),nout=i0)
C ... Override verbos w/ -pr commmand-line arg
      if ((cmdopt('--pr',4,0,a) .or. cmdopt('-pr',3,0,a))
     .    .and. .not. cmdopt('--pra',5,0,a)) then
        i = 4
        if (cmdopt('-pr',3,0,a)) i = 3
        i = a2vec(a,len(a),i,2,', ',2,2,5,it,verbos)
        if (i <= 0) call rxs('error parsing switch ',a)
C       i0 = max(i0,i)
        i0 = i
      endif
C     Copy verbosities to print stack
      if (i0 >= 1) then
        do  i = 0, 5
          call sprt(i,verbos(min(i+1,i0)))
        enddo
      endif
      if (procid /= master) then
        call pshpr(1)
        call pshpr(1)
        do  i = 1, 4
          call sprt(i,0)
        enddo
      endif

      nm='IO_IACTIV'; call gtv(trim(nm),tksw2(prgn,nm),
     .  io_iactive, note='Turn on interactive mode.'//
     .  '%N   May also be controlled from the command-line:'//
     .  '  --iactiv  or  --iactiv=no',
     .  def_i4=0)
      if( cmdopt('--no-iactiv',11,0,a)) io_iactive = 0
      if( cmdopt('--iactiv',8,0,a))     io_iactive = 1
      if( cmdopt('--iactiv=no',11,0,a)) io_iactive = 0
      if (io_iactive > 1) io_iactive = 1
      call initqu(io_iactive)

      nm='IO_TIM'; call gtv(trim(nm),tksw2(prgn,nm),io_tim,
     .  note='Turns CPU timing log.  Value sets tree depth.'//
     .  '%N   Optional 2nd arg prints CPU times on the fly.'//
     .  '%N   May also be controlled from the command-line:  --time=#1[,#2]',
     .  def_i4v=(/0,0/),nmin=2,nout=i0)
C     Override with '--time=' commmand-line arg
      if ( cmdopt('--time=',7,0,a) ) then
        i = 7
        i = a2vec(a,len(a),i,2,', ',2,2,2,it,io_tim)
        if (i < 0) call rxs('error parsing switch ',a)
        i0 = max(i0,i)
      endif
      if ( i0 >=1 ) call tcinit(io_tim(2),io_tim(1))

      nm='IO_MAXMEM'; call gtv(trim(nm),tksw2(prgn,nm),maxmem,def_i4=0,
     .  note='impose upper bound to dynamic memory allocation (0 => no bound)')
      i = 0
      nm='IO_SHOMEM'; call gtv(trim(nm),tksw2(prgn,nm),i,def_i4=0,note=
     .  'display dynamic memory allocation for arrays exceeding SHOMEM MB')
      if (cmdopt('--shomem',8,0,a)) i = 1
      if (i > 0) then
        call set_allocvb(1); call idasiz(i)
      endif

C --- Header ---
      nm='HEADER'; call gtv(trim(nm),tksw2(prgn,nm),header,note=
     .  '$1Contents displayed at beginning of program execution',nout=nout)
      if (nout > 0) then
        if (stdl>0 .and. iprint()>0) call headl2(prgnam,0,stdl)
        outs  = 'HEADER '//trim(header)
        if (procid == master) then
          if (io_show /= 0) write(stdo,'(1x)')
          write(stdo,'(1x,a)') trim(outs)
C          if (io_show /= 0) write(stdo,'(1x)')
C          if (stdl>0 .and. iprint()>0) write(stdl,'(1x,a)') trim(outs)
        endif
      endif

C --- Struc ---
      sw = tkswp(prgn,express,'STRUC'); if (sw >= 2) goto 59
      if (io_show+io_help/=0) call info0(2,1,0,
     .  ' --- Parameters for crystal structure ---')

C     gtv call is useful for printout, and also checks whether category is present
C      call gtv_none(trim('STRUC'),sw,nono,Texist=ltmp,note=
C     .  'Parameters for crystal structure')

C
      nm='STRUC_NSPEC'; call gtv(trim(nm),tksw2(prgn,nm),nspec,
     .  note='Number of species to read from SPEC category.'//
     .  '%N%3fIf not present, NSPEC will be obtained by:'//
     .  '%N%5f- identifying species from the SITE file if it is read, or by'//
     .  '%N%5f- counting entries in the SPEC category' )

      nm='STRUC_FILE'; call gtv(trim(nm),tksw2(prgn,nm),outs,nmin=10,
     .  nout=nout, note=
     .  '$$e(also alias SITE_FILE; see below)%N%3f$$e'//
     .  'Name of site file containing basis and lattice information.'//
     .  '%N%3fRead NBAS, PLAT, and optionally ALAT from site file, '//
     .  'if specified.'//
     .  '$$s%N%3fOtherwise, they are read from the ctrl file.')
      if (nout == 1) then
        xx = nullr
C       Create a separate species table since it might be reordered
        if (nspec == NULLI) allocate(slabll(mxspec))
        if (procid == master) then
          i0 = 3000
          if (nspec == NULLI) then  ! Read species labels
            i0 = 131000  ! Load species labels
            nspec = 0
          endif
C         Does nothing except return nbas
          j = iosite(256000+135000,3d0,0,trim(outs),i,slabll,alat,plat,nbas,
     .               nspec,xx,xx,xx,xx,xx,xx,xx)
          if (nbas <= 0) call rx('failed to read nbas from site file "'//trim(outs)//'"')
          j = iosite(i0,3d0,0,trim(outs),i,slabll,alat,plat,nbas,
     .               nspec,xx,xx,xx,xx,xx,xx,xx)
          if (io_show>0) call info2(2,0,0,
     .      ' ... site file contains %i species and %i atoms',nspec,nbas)
        endif
        j = nspec*len(slabll); call mpibc1(j,1,2,mlog,'','')
        call mpibcc(slabll,j,mlog,' ',' ')
        call mpibc1(nbas,1,2,mlog,'readctrl','nbas')
        call mpibc1(nspec,1,2,mlog,'readctrl','nspec')
        call mpibc1(alat,1,4,mlog,'readctrl','alat')
        call mpibc1(plat,9,4,mlog,'readctrl','plat')
        call rxx(nbas==0.or.nspec==0,'no sites read from file '//trim(outs))

C       Exclude further attempt to read these data
        if (nspec > 0) call tkexclude('STRUC_NSPEC')
        call tkexclude('STRUC_NBAS')
        call tkexclude('STRUC_NSPEC')
        if (alat /= NULLR) call tkexclude('STRUC_ALAT')
        call tkexclude('STRUC_PLAT')
      endif

C     call shosyv(0,0,0,stdo)
      if (alat == NULLR) then
        nm='STRUC_ALAT'; call gtv(trim(nm),tksw2(prgn,nm),
     .    alat, note= 'Scaling of lattice vectors, in a.u.')
      endif
      if (nbas == NULLI) then
        nm='STRUC_NBAS'; call gtv(trim(nm),tksw2(prgn,nm),
     .    nbas,note='Size of basis')
        nm='STRUC_PLAT'; call gtv(trim(nm),tksw2(prgn,nm),
     .    plat, nmin=9, nout=nout, note=
     .    'Primitive lattice vectors, in units of alat')
      endif
      if (io_help == 0 .and. express >= 0) then
        avw = avwsr(plat,alat,vol,nbas)
        call mkqlat(plat,qlat,xx)
      endif

C ... Count number of species in SPEC category
      if (io_help == 0 .and. nspec == NULLI .and. express >= 0) then
      nm='SPEC_ATOM'; sw = tksw2(prgn,nm)
      if (sw < 2) then
        j = 0; nspec = 0
        do  while (nspec <= 0)
          j = j+1; jj= (/1,j/)
          if (.not. debug) call pshpr(1)
          call gtv(trim(nm),0,nono,Texist=ltmp,cindx=jj)
          if (.not. debug) call poppr
          if (.not. ltmp) nspec = j-1
        enddo
        if (io_show>0) call info2(2,0,0,' ... found %i species in SPEC category',nspec,0)
      endif
      endif

C     Extra site positions for e.g. point multipoles
      nm='STRUC_NBASP'; call gtv(trim(nm),tksw2(prgn,nm),nbasp,
     .  def_i4=nbas,note='nbas + no. of point multipoles')
      if (nbasp == NULLI) nbasp = nbas
      nsite = nbasp

      nm='STRUC_SLAT'; call gtv(trim(nm),tksw2(prgn,nm),lat_slat,
     .  nmin=9,note='Supercell lattice vectors')
      nm='STRUC_DALAT'; call gtv(trim(nm),tksw2(prgn,nm),
     .  dalat, def_r8=0d0, note='added to alat after input is read')
      i = 3
      if (lfp == 1 .or. prgnam=='LMFA') i=4
      if (prgnam=='LMMC') i=3
      nm='STRUC_NL'; call gtv(trim(nm),tksw2(prgn,nm),nl,def_i4=i,
     .  note='global default lmax+1 for basis and augmentation')

C     Assume no subsequent tokens in STRUC belong to 'express' mode
      if (express < 0) goto 59

C ... Lattice distortion or rotation
      sw = tksw2(prgn,'STRUC_SHEAR'); lat_gam(1) = NULLR
      if (sw < 2) then
        j = 0
        nm='STRUC_SHEAR'; call gtv(trim(nm),0,lat_gam,nmin=4,or=T,
     .    note='Volume-conserving shear of PLAT (1=ideal)')
        if (lat_gam(1) /= NULLR .and. io_help == 0) then
          j = 0
          goto 880
        endif
        nm='STRUC_ROT'; call gtv(trim(nm),0,outs,nmin=10,nout=nout,or=T,
     .    note='Rotate PLAT.  Use standard Questaal conventions for ROT.')
        if (nout /= 0 .and. io_help == 0) then
          call a2rotm(outs,.false.,iprint()-10,lat_dist)
          j = 2
          goto 880
        endif
        nm='STRUC_DEFGRD'; call gtv(trim(nm),0,lat_dist,nmin=9, or=T,
     .    note='General shear of PLAT (3,3) matrix')
        if (lat_dist(1) /= NULLR .and. io_help == 0) then
          j = 2
          goto 880
        endif
        nm='STRUC_STRAIN'; call gtv(trim(nm),sw,lat_dist(1:6),nmin=6,
     .    note='strain PLAT.  Supply 6 Voigt strain matrix elements')
        if (lat_dist(1) /= NULLR .or. io_help /= 0) then
          nm='STRUC_ALPHA'; call gtv(trim(nm),1,lat_dist(7),nout=nout,
     .      note=
     .      'Amplitude of (Voigt) strain.  Only read if STRAIN input')
        if (nout /= 0) j = 3
        endif
  880   continue
        lat_ldist = j
      endif
      if (lat_gam(1) == NULLR) lat_gam=(/0d0,0d0,1d0,1d0/)

C ... End of struc tokens
   59 continue

C --- Options ---
C     Skip if category is ignored, or if pass=-1 and no tags to read
      if (tkswp(prgn,express,'OPTIONS') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for Options ---')

      nm='OPTIONS_HF'; call gtv(trim(nm),tksw2(prgn,nm),lcd2,def_lg=F,
     .  note='T for non-self-consistent Harris')
      i=0; nm='OPTIONS_SHARM'; call gtv(trim(nm),tksw2(prgn,nm),i,def_i4=0,
     .  note='Set to T, or 1 or 2 to use true spherical harmonics'//
     .  '%N%3f1 Spherical harmonics, m ordered -l:l'//
     .  '%N%3f2 Spherical harmonics, m ordered l:-l')
      if (i == 2) call setmorderdefault(0)
      lham256 = i > 0

      nm='OPTIONS_FRZ'; call gtv(trim(nm),tksw2(prgn,nm),lcd1,def_lg=F,
     .  note= 'Freeze core')
      nm='OPTIONS_NMCORE'; call gtv(trim(nm),tksw2(prgn,nm),lcd16,def_lg=F,
     .  note= 'Nonmagnetic core')
      nm='OPTIONS_WRONSK'; call gtv(trim(nm),tksw2(prgn,nm),ltmp,def_lg=F,
     .  note='Enforce Wronskian condition on partial waves')
      xv(1) = dglob('wronsk',dble(isw(ltmp)),1)

      nm='OPTIONS_RQUAD'; call gtv(trim(nm),tksw2(prgn,nm),lrquad,
     .  def_i4=0,note='Quadrature for radial mesh'//
     .  '%N%3f0 3-point (Simpson''s rule)'//
     .  '%N%3f1 5-point (Boole''s rule; NOT implemented)'//
     .  '%N%3f2 7-point')

      nm='OPTIONS_XCQS'; call gtv(trim(nm),tksw2(prgn,nm),lcd64,
     .  note='XC potential in reciprocal space')

      nm='OPTIONS_WRONSK'; call gtv(trim(nm),tksw2(prgn,nm),lcd64,
     .  note='Enforce Wronskian condition on partial waves')

C      nm='OPTIONS_TPAN'; call gtv(trim(nm),tksw2(prgn,nm),lham1,
C     .  def_lg=F,note='two-panel')
      nm='OPTIONS_SAVVEC'; call gtv(trim(nm),tksw2(prgn,nm),lham64,
     .  def_lg=F,note='Save eigenvectors on disk')

      nm='OPTIONS_RMINES'; call gtv(trim(nm),tksw2(prgn,nm),rmines,
     .  def_r8=1d0,note='Minimum MT radius when finding new ES')
      nm='OPTIONS_RMAXES'; call gtv(trim(nm),tksw2(prgn,nm),rmaxes,
     .  def_r8=2d0,note='Maximum MT radius when finding new ES')
      nm='OPTIONS_NESABC'; call gtv(trim(nm),tksw2(prgn,nm),nesabc,
     .  nmin=3,def_i4v=(/100,100,100/),
     .  note='No. divisions when searching for empty spheres')

C      nm='OPTIONS_STONER'; call gtv(trim(nm),tksw2(prgn,nm),
C     .  lstonr, note='Generalised Stoner rigid band calculation '//
C     .  '%N%3fSecond argument is number of points; third for graphical output')

      nm='OPTIONS_SHORBZ'; call gtv(trim(nm),tksw2(prgn,nm),lfp2,
     .  def_lg=T, note= 'globally shorten q vectors')

      outs = ' '
      nm='OPTIONS_Q'; call gtv(trim(nm),tksw2(prgn,nm),outs,nmin=10,
     .  note='Use Q=show, Q=atom, or Q=band to quit after '//
     .  'input, sphere calc or band pass.'//
     .  '%N%3fCommand-line `--quit=string'' overrides input file')
      call locase(outs)
      if (outs=='show' .or. cmdopt('--quit=show',11,0,a)) then
        quit=1
      elseif (outs=='atom' .or. cmdopt('--quit=atom',11,0,a)) then; quit=2
      elseif (outs=='band' .or. cmdopt('--quit=band',11,0,a)) then; quit=4
      elseif (outs=='ham' .or. cmdopt('--quit=ham',10,0,a)) then; quit=8
      elseif (outs=='dos' .or. cmdopt('--quit=dos',10,0,a)) then; quit=16
      elseif (outs=='rho' .or. cmdopt('--quit=rho',10,0,a)) then; quit=32
      elseif (outs=='pot' .or. cmdopt('--quit=pot',10,0,a)) then; quit=64
      elseif (outs=='bc' .or. cmdopt('--quit=bc',9,0,a)) then; quit=128
      elseif (cmdopt('--quit=',7,0,a)) then
        call rx('--quit=xx not recognized : xx must be one of SHOW, ATOM, HAM, POT, BAND, DOS, BC, or RHO')
      elseif (outs==' ') then
      else
        call rx('OPTIONS_Q= must contain one of SHOW, ATOM, HAM, POT, BAND, DOS, or RHO')
      endif
      nm='OPTIONS_SCR'; call gtv(trim(nm),tksw2(prgn,nm),lscr,
     .  def_i4=0,note='Use scr to accelerate convergence:'//
     .  '%N%3f0 do nothing'//
     .  '%N%3f1 Make ASA static response function (see documentation)'//
     .  '%N%3f2 Use response to screen output q and ves'//
     .  '%N%3f4 Use model response to screen output q'//
     .  '%N%3f6 Use response to screen output ves only'//
     .  '%N%5fAdd 1 to combine mode 1 with another mode'//
     .  '%N%5fAdd 10*k to compute intra-site contribution to vbare'//
     .  ' each kth iteration'//
     .  '%N%5fAdd 100*k to compute response function'//
     .  ' on every kth iteration')

C     ASA-specific
      sw = tksw2(prgn,'OPTIONS_ASA')
      if (io_show+io_help/=0 .and. sw<2)
     .  call info0(2,1,0,' ... The following are ASA-specific')
      if (sw < 2) then
C     This call does nothing except for printout
      nm='OPTIONS_ASA'; call gtv_none(trim(nm),sw,nono,Texist=ltmp,
     .  note= 'Contains the following ASA-specific tokens:')
      nm='OPTIONS_ASA_ADNF'; call gtv(trim(nm),tksw2(prgn,nm),lham4,
     .  def_lg=F,note= 'Turn on automatic downfolding')
      nm='OPTIONS_ASA_NSPH'; call gtv(trim(nm),tksw2(prgn,nm),lasa32,
     .  note='generate multipole moments in output density')
      nm='OPTIONS_ASA_TWOC'; call gtv(trim(nm),tksw2(prgn,nm),lham3,
     .  def_i4=0,note='Two-center ASA hamiltonian.'//
     .  '  Use TWOC=3 for 2C+pert corr')
      j = 0
      nm='OPTIONS_ASA_GAMMA'; call gtv(trim(nm),tksw2(prgn,nm),j,
     .  def_i4=0,note=
     .  '>0: gamma representation'//
     .  '%N%4f1: gamma is spin dependent for spin polarized cases'//
     .  '%N%4f2: gamma = gamma-bar = (gamma(1)+gamma(nsp))/2'//
     .  '%N%4f   (some noncollinear algorithms automatically convert gamma to gamma-bar)'//
     .  '%N%4f4: convert to S^gamma by rotation of S^alp instead of energy scaling (lmgf,lmpg)'//
     .  '%N%4f5: same as 4, but convert to S^gamma-bar')
      call sanrg(io_help==0,j,0,5,' rdctrl:','GAMMA')
      if (j == 3) call rxi('m_rdctrl: OPTIONS_ASA_GAMMA has illegal value',j)
      lham128 = j > 0
      lasa512 = j == 2 .or. j == 5
      lasa1024 = j == 4 .or. j == 5
      nm='OPTIONS_ASA_CCOR'; call gtv(trim(nm),tksw2(prgn,nm),lasa4,
     .  def_lg=T,note= 'Turn on combined correction')
      nm='OPTIONS_ASA_ELIN'; call gtv(trim(nm),tksw2(prgn,nm),asa_elin,
     .  def_r8=0d0,note=
     .  'Energy to linearize for CCOR (2nd gen ASA 2C hamiltonian)')
      nm='OPTIONS_ASA_NEWREP'; call gtv(trim(nm),tksw2(prgn,nm),lham8,
     .  def_lg=F,note=
     .  'Interactively transform representation, (2nd gen ASA)')
      nm='OPTIONS_ASA_NOHYB'; call gtv(trim(nm),tksw2(prgn,nm),lham16,
     .  def_lg=F,note='Turns off hybridisation ')
      nm='OPTIONS_ASA_MTCOR'; call gtv(trim(nm),tksw2(prgn,nm),lasa64,
     .  def_lg=F,note='Turns on Ewald MT correction')
      nm='OPTIONS_ASA_QMT'; call gtv(trim(nm),tksw2(prgn,nm),
     .  zbak(2),def_r8=0d0,note=
     .  'Override standard background charge for'//
     .  ' Ewald MT correction%N%3fInput only meaningful if MTCOR=T')
      endif
      endif ! OPTIONS

C --- Hamiltonian parameters ---
      if (tkswp(prgn,express,'HAM') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for hamiltonian ---')

      sw = tkswp(prgn,express,'HAM_GMAX')
      if (sw<2) then
        nm='HAM_GMAX'; call gtv_r8(trim(nm),sw,lat_gmax,nmin=1,
     .    nout=nout,note='Energy cutoff for plane-wave mesh (Ry)',or=T)
        if (nout > 0) then
          sw = 2
        else
          lat_gmax = 0
        endif
        if (sw/=2) then
          nm='HAM_GMAXDA'; call gtv(trim(nm),sw,lat_gmax,nmin=1,nout=nout,or=T,
     .      note='Energy cutoff for plane-wave mesh, scaled by ALAT/(ALAT+DALAT)')
          if (nout > 0) then
            sw = 2
            lat_gmax = lat_gmax * alat/(alat+dalat)
          else
            lat_gmax = 0
          endif
        endif
        nm='HAM_FTMESH'; call gtv(trim(nm),sw,ftmesh,nout=nout,note=
     .    'No. divisions for plane-wave mesh '//
     .    'along each of 3 lattice vectors.'//
     .    '%N%3fSupply one number for all vectors or a separate '//
     .    'number for each vector.')
        if (nout > 0) then
          call tkexclude('HAM_GMAX') ! Also exclude gmax
        else
          ftmesh = 0  ! Override NULLI (express mode)
        endif
        call fill3in(nout,ftmesh)
      endif
      nm='HAM_TOL'; call gtv(trim(nm),tksw2(prgn,nm),tolft,
     .  def_r8=1d-6, note='w.f. tolerance for FT mesh')
C   Parameters basopt(*) used generation of LMTO basis parms
C    1 = autob
C        1s   digit 1 or 3 (lmfa) Autogenerate RSMH,EH
C                   2 or 4 (lmfa) Autogenerate RSMH,EH, RSMH2,EH2
C                   1 or 3 (lmf)  Read RSMH,EH,RSMH from basp file
C                   2 or 4 (lmf)  READ RSMH,EH, RSMH2,EH2 from basp file
C        10s  digit 1 (lmfa) Find and estimate PZ from free atom wf
C                     (lmf)  Read P from basp file
C        100s digit 1 (lmfa) Estimate P from free atom wf
C                     (lmf)  Read P from basp file
C        (used in autoset basis)
C    2 = global cutoff to lmxb, 1st kappa
C        (no longer used)
C    3 = elcut Set local orbital PZ when E>autoecut
C        (used in autogenerating PZ)
C    4 = rsmmx: maximum rsm, in units of rmt
C        (used in autogenerating parameters)
C    5 = ehmx : maximum eh, Ry
C        (used in autogenerating parameters)
C    6 = esprd: default spread in EH,EH2
C        (used in autogenerating EH,EH2)
C    7 = modeV: specifies a mode to modify potential
C        (used in autogenerating EH,EH2)
C    8 = vbar : parameter in V modification
C        (used in autogenerating EH,EH2)
C    9 = autoqcut Set local orbital PZ when q(r>rmt)
C        (used in autogenerating PZ)
C   10 = set defaults for GW calculation
      if (express <= 0) then ! initialize these quantities only once
        basopt = 0; basopt(4)=2d0/3
        basopt(5)=NULLR; basopt(6)=0.8d0; basopt(10)=0d0
        nsp = 1
      endif
      if (prgnam == 'MMAG') then
        nsp = 2
        lncol1 = T
      endif
      nm='HAM_NSPIN'; call gtv(trim(nm),tksw2(prgn,nm),nsp,
     .  def_i4=1,note='Set to 2 for spin polarized calculations')
      call sanrg(io_help==0,nsp,1,2,'rdctrl','nsp')

      nm='HAM_REL'; call gtv(trim(nm),tksw2(prgn,nm),lrel,def_i4=1,nout=nout,
     .  note='0 for nonrelativistic Schrodinger equation'//
     .  '%N%3f1 for scalar relativistic Schrodinger equation'//
     .  '%N%3f2 for Dirac equation'//
     .  '%N%3f  10s digit 1: compute core density with full Dirac equation'//
     .  '%N%3f  10s digit 2: Like 1, but neglect coupling (1,2) pairs in 4-vector')

      call sanrg(io_help==0.and.nout >= 0,mod(lrel,10),0,2,' rdctrl:','HAM_REL')
C     Fully relativistic => spin-orbit coupling
      if (mod(lrel,10)==2) lncol4 = T

C      if (io_help /= 0 .and. tkswp(prgn,express,'HAM_SO') < 2)
C     .    call info0(2,1,0,' - To read the magnetic '//
C     .    'parameters below, HAM_NSPIN must be 2')
      if (nsp==2 .or. io_help/=0) then
      nm='HAM_SO'; call gtv_i4(trim(nm),tksw2(prgn,nm),
     .    lso,def_i4=0,note=
     .  'Spin-orbit coupling (for REL=1).'//
     .  ' You must also set NSPIN=2.'//
     .  '%N   Automatically sets HAM_NCOL to T.'//
     .  '%N%3f0 : no SO coupling'//
     .  '%N%3f1 : Add L.S to hamiltonian'//
     .  '%N%3f2 : Add Lz.Sz only to hamiltonian'//
     .  '%N%3f3 : Like 2, but also include <L.S-LzSz> by perturbation,'
     .           //' maintaining independent spins'//
     .  '%N%3f4 : Compute entire <L.S> by perturbation'//
     . '%N%2f11 : Same as 1, but additionally decompose SO by site')
      if (lso == 11) then
        lsos = T
        lso = 1
      endif
      call sanrg(io_help==0,lso,0,4,' rdctrl:','SO')
      if (lso == 1) lncol4 = T
      if (lso == 2) lncol32= T
      if (lso == 3) lncol64= T
      if (lso == 4) lncol128= T

C      nm='HAM_SHFAC'; call gtv(trim(nm),tksw2(prgn,nm),lham512,
C     .  def_lg=F, note=
C     .  'Read shfac file, if it exists')

      nm='HAM_NONCOL'; call gtv(trim(nm),tksw2(prgn,nm),lncol1,
     .  def_lg=F, note='Noncollinear magnetism')
      nm='HAM_SAXIS'; call gtv(trim(nm),tksw2(prgn,nm),lrsa,def_i4=0,
     .  note='Governs how noncollinear density is treated'//
     .  '%N%3f0 : Collinear rho; output noncoll rho only for Etot'//
     .  '%N%3f1 : Rigid spin approximation: input local rho '//
     .           'points along fixed given axis'//
     .  '%N%3f2 : Rigid spin approximation: local rho '//
     .           'rigid along average magnetization'//
     .  '%N%3f3 : Fully noncollinear density'//
     . '%N%2f10 : Add 10 to require site local axis '//
     .           ' parallel to |M| (default is parallel to M or -M)'//
     .'%N%1f100 : Add 100 to spin-average core density')
      nm='HAM_SS'; call gtv(trim(nm),tksw2(prgn,nm),ham_qss,nmin=4,
     .  nout=nout,note='Magnetic spin spiral, direction vector and angle')
      lncol2 = (nout == 4)
C     If spin-orbit or SS, also turn on noncollinear
      if (lncol2 .or. lncol4) lncol1 = T
      if (lncol1 .or. io_help/=0) then
C      nm='HAM_SPINTRQ'; call gtv(trim(nm),tksw2(prgn,nm),lncol16,
C     .  def_lg=F, note='Calculate spin torques (requires NONCOL=t)')
      endif
      nm='HAM_BFIELD'; sw = tksw2(prgn,nm); call gtv(trim(nm),sw,j,
     .  def_i4=0,nout=nout,note=
     .  'Apply site-dependent magnetic field (for REL=1). '//
     .  'Fields read from file bfield.'//
     .  '%N%3f0 : no external field coupling'//
     .  '%N%3f1 : Add external B field to hamiltonian (file bfield)'//
     .  '%N%3f2 : Add Bz.Sz only to hamiltonian')
      if (sw < 2 .and. nout > 0) then
        call sanrg(io_help==0.and.express >= 0,j,0,2,' rdctrl:','BFIELD')
        if (j == 1) lncol8 = T
        if (j == 2) lBz = T
      endif
      endif
      nm='HAM_BXCSCAL'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lbxc1,
     .  def_i4=0,note='Scale magnetic part of LDA XC field'//
     .  ' (Compare to HAM_BFIELD).'//
     .  '%N%3fScalings are site-dependent; coefficients read '//
     .  'from first column of file "shfac".'//
     .  '%N%3f0 : do not apply any scaling'//
     .  '%N%3f1 : scale Bxc by 1+shfac(1,ib).'//
     .  '%N%3f2 : scale Bxc by [1+shfac(1,ib)^2]^1/2'//
     .  ' (used for DLM constraining fields)'//
     .  '%N%7fib is site index (lmf) or class index (lm, lmgf)'//
     .  '%N%3fAdd 10 to scale interstitial part also (lmf only)')
      if (ctrl_lbxc1 > 0) lham512 = .true.

C      nm='HAM_SOSCAL'; call gtv(trim(nm),tksw2(prgn,nm),j,
C     .  def_i4=0,nout=nout,note='Scale L.S coupling.'//
C     .  '%N%3fScalings are site-dependent; coefficients read '//
C     .  'from second column of file "shfac".'//
C     .  '%N%3f0 : do not apply any scaling'//
C     .  '%N%3f1 : scale L.S by 1+shfac(2,ib).'//
C     .  '%N%7fib is site index (lmf) or class index (lm, lmgf)')
C      if (j /= 0 .and. nout >= 0) lham512 = .true.

      nm='HAM_BXC0'; call gtv(trim(nm),tksw2(prgn,nm),
     .  ctrl_lbxc16,def_lg=F,note='Render LDA Bxc=0')

      nm='HAM_FRZWF'; call gtv(trim(nm),tksw2(prgn,nm),frzwf,def_lg=F,
     . note='Set to freeze augmentation wave functions for all species')
      nm='HAM_V0BETA'; call gtv(trim(nm),tksw2(prgn,nm),v0beta,
     .  def_r8=1d0, note='mixing beta for admixing sphere potential V0 defining phi,phidot')
      nm='HAM_FORCES'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lfrce,
     .  def_i4=0,note=
     .  'Controls the ansatz for density shift in force calculation.'//
     .  '%N%3f-1 no force%3f0 no shift'//
     .  '%N%3f 1 free-atom shift  12 screened core+nucleus')
      nm='HAM_ELIND'; call gtv(trim(nm),tksw2(prgn,nm),elind,
     .  def_r8=elind,note='Lindhard energy for model screening')
      nm='HAM_XCFUN'; call gtv(trim(nm),tksw2(prgn,nm),llxcf,
     .  def_i4v=(/2,0,0/),
     .  note='Specifies the local exchange correlation functional.'//
C     .  '%N%3fUse a functional from libxc if arg #1=0.'//
     .  '%N%3f0,#2,#3 : Use libxc exchange functional #2'//
     .  ' and correlation functional #3.'//
     .  '%N%5farg #1 is nonzero, code uses an internally coded functional:'//
     .  '%N%3f1 Ceperly-Alder'//
     .  '%N%3f2 Barth-Hedin (ASW fit)'//
     .  '%N%3f3 PW91 (same as PBE) %N%3f4 PBE (same as PW91)')
      nm='HAM_GGA'; call gtv(trim(nm),tksw2(prgn,nm),gga,
     .  nout=nout, def_i4=0, note='Specifies GGA functional'//
     .  ' (used only if XCFUN arg #1 is > 0).'//
     .  '%N%3f0 LSDA%N%3f1 Langreth-Mehl%N%3f2 PW91'//
     .  '%N%3f3 PBE%N%3f4 PBE with Becke exchange')
      if (llxcf(1) == 0) then  ! a libxc functional
        call rxx(llxcf(2)==0,'m_rdctrl nonsensical XC functional')
        lxcf = llxcf(2)*2**16 + llxcf(3)
      else
        lxcf = llxcf(1) + 100*gga
      endif
      nm='HAM_AUTOBAS_GW'; call gtv(trim(nm),tksw2(prgn,nm),j,
     .  nout=nout,def_i4=0,note='If >0, tailor default '//
     .  'inputs for a GW calculation')
      if (nout > 0) then
        iapflt = 1
        if (j > 0) then
          iapflt(1) = 2
          basopt(10) = 10
        endif
      endif
      if (prgnam=='LMFA') then
        nm='HAM_AUTOBAS_LMTO';call gtv(trim(nm),tksw2(prgn,nm),i,
     .    nout=nout,def_i4=0,note='Controls lmfa''s autogeneration of '/
     .    /'LMTO basis parameters (RSMH,EH,RSMH2,EH2)'//
     .    '%N%3f0 standard minimal basis, (same as 3)'//
     .    '%N%3f1 hyperminimal basis (use in conjunction w/APWs)'//
     .    '%N%3f2 hyperminimal basis + 1 higher l'//
     .    '%N%3f3 standard minimal basis, typically spd, f for Z>Kr'//
     .    '%N%3f4 standard basis, similar to (3), f for Z>Ar'//
     .    '%N%3f5 Large basis')
        if (nout > 0) then
          call sanrg(io_help==0,i,0,5,'rdctrl','HAM_AUTOBAS_LMTO')
          basopt(10) = basopt(10) + (getdig(i,0,10)-getdig(nint(basopt(10)),0,10))
        endif
        nm='HAM_AUTOBAS_MTO';call gtv(trim(nm),tksw2(prgn,nm),i,
     .  nout=nout,def_i4=0,note='Controls lmfa''s autogeneration of '//
     .    'LMTO basis parameters (RSMH,EH,RSMH2,EH2)'//
     .    '%N%3f0 do not autogenerate basis parameters'//
     .    '%N%3f1 or 3 1-kappa parameters with Z-dependent LMX'//
     .    '%N%3f2 or 4 2-kappa parameters with Z-dependent LMX'//
     .    '%N%3f10s digit'//
     .    '%N%3f0 classical method'//
     .    '%N%3f1 Jackson''s classical turning point method')
        if (nout > 0) then
          basopt(10) = basopt(10) + 100*(getdig(i,1,10)-getdig(nint(basopt(10)),2,10))
          autob = autob + getdig(i,0,10)-getdig(autob,0,10)
        endif
        nm='HAM_AUTOBAS_PNU'; call gtv(trim(nm),tksw2(prgn,nm),i,
     .  nout=nout,def_i4=0,note='Controls lmfa''s autogeneration of '//
     .    'log derivative parameters P'//
     .    '%N%3f0 do not generate P'//
     .    '%N%3f1 Find P for l<lmxb from free atom wave function')
        if (nout > 0) then
          call sanrg(io_help==0,i,0,1,'rdctrl','HAM_AUTOBAS_PNU')
          autob = autob + 100*(getdig(i,0,10)-getdig(autob,2,10))
        endif
        nm='HAM_AUTOBAS_LOC'; call gtv(trim(nm),tksw2(prgn,nm),i,
     .  nout=nout,def_i4=0,note='Controls lmfa''s autogeneration '//
     .    'of local orbital parameters PZ'//
     .    '%N%3f0 do not autogenerate PZ'//
     .    '%N%3f1 or 2   autogenerate PZ'//
     .    '%N%3f1 On input, nonzero values from ctrl file '//
     .    'take precedence over basis file'//
     .    '%N%3f2 On input, contents of ctrl file are ignored')
        if (nout > 0 .and. i > 0) then
C         autob = autob + 20
          call sanrg(io_help==0,i,0,2,'rdctrl','HAM_AUTOBAS_LOC')
          autob = autob + 10*(getdig(i,0,10)-getdig(autob,1,10)) ! This one allows ctrl to override default
        endif
        nm='HAM_AUTOBAS_RSMMX'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(4),nout=nout,def_r8=2d0/3,
     .    note='Upper bound when autogenerating LMTO smoothing radius RSMH'//
     .    ' (units of RSM)')
CJJ
        nm='HAM_AUTOBAS_RSMN'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(2),nout=nout,def_r8=0.9d0,
     .    note='Lower bound when autogenerating LMTO smoothing radius RSMH'//
     .    ' (atomic units)')
        nm='HAM_AUTOBAS_EHMX'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(5),nout=nout,def_r8=NULLR,
     .    note='Upper bound when autogenerating LMTO Hankel energy EH')
        nm='HAM_AUTOBAS_ELOC'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(3),nout=nout,def_r8=-2d0,
     .    note='Local orbitals specified from atomic orbitals '//
     .    'where E is more shallow than ELOC')
        nm='HAM_AUTOBAS_QLOC'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(9),nout=nout,def_r8=5d-3,
     .    note='Local orbitals specified from atomic orbitals '//
     .    'where q(r>rmt) exceeds QLOC')
C        nm='HAM_AUTOBAS_LMXB'; call gtv(trim(nm),tksw2(prgn,nm),lmxb_cut,
C     .    nout=nout,note=
C     .    'Global maximum l when autogenerating LMTO basis parms'//
C     .    '%N%3f1st argument for 1st kappa, 2nd arg for 2nd kappa')
        nm='HAM_AUTOBAS_PFLOAT'; call gtv(trim(nm),tksw2(prgn,nm),pnudef,
     .  nout=nout,def_i4=iapflt(1),note=
     .    'Default starting values of P'//
     .    '%N%5f0 Use version 6 defaults'//
     .    '%N%5f1 Use defaults tailored for LDA'//
     .    '%N%5f2 Use defaults tailored for GW')
        if (nout > 0) then
          call sanrg(io_help==0,mod(pnudef,10),0,2,'rdctrl',
     .      '1st argument (1s digit) HAM_AUTOBAS_PFLOAT')
        endif
C       extra parameters for c.t.p. basis setup scheme
        nm='HAM_AUTOBAS_EIN'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(11),nout=nout,def_r8=0.d0,
     .    note='Lower energy bound for setting up basis using V(r)'//
     .    ' (Ryd)')
        nm='HAM_AUTOBAS_EOUT'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(12),nout=nout,def_r8=0.d0,
     .    note='Upper energy bound for setting up basis using V(r)'//
     .    ' (Ryd)')
        nm='HAM_AUTOBAS_EH'
        call gtv(trim(nm),tksw2(prgn,nm),basopt(13),nout=nout,def_r8=0.d0,
     .    note='Hankel energy used in V(r) basis scheme (Ryd)')
      else
        nm='HAM_AUTOBAS_MTO';call gtv(trim(nm),tksw2(prgn,nm),i,
     .    Texist=ltmp,nout=nout,def_i4=0,note='Autoset basis:  '//
     .    'controls what part of MTO basis is read from file basp'//
     .    '%N%3f0       No parameters are read from basp'//
     .    '%N%3f1 or 3: 1-kappa parameters may be specified by basp'//
     .    '%N%3f2 or 4: 2-kappa parameters may be specified by basp'//
     .    '%N%3f1 or 2: Nonzero parameters from ctrl file take '//
     .    'precedence over basp input'//
     .    '%N%3f3 or 4: Nonzero parameters from basp file take '//
     .    'precedence over ctrl file')
        if (nout > 0) then
          autob = autob + getdig(i,0,10)-getdig(autob,0,10)
        endif
        nm='HAM_AUTOBAS_PNU'; call gtv(trim(nm),tksw2(prgn,nm),i,
     .    nout=nout,def_i4=0,note='Read log derivative '//
     .    'parameters P from basis file basp'//
     .    '%N%3f0 Do not read parameters P from basp'//
     .    '%N%3f1 Read parameters P from basp, if they are present')
        if (nout > 0) then
          call sanrg(io_help==0,i,0,1,'rdctrl','HAM_AUTOBAS_PNU')
          autob = autob + 100*(getdig(i,0,10)-getdig(autob,2,10))
        endif
        nm='HAM_AUTOBAS_LOC'; call gtv(trim(nm),tksw2(prgn,nm),i,
     .    nout=nout,def_i4=0,note='Read autogenerated '//
     .    'local orbital parameters PZ from basis file'//
     .    '%N%3f1 or 2 read parameters PZ'//
     .    '%N%3f1 Nonzero values from ctrl file '//
     .    'take precedence over basis file input'//
     .    '%N%3f2 Contents of ctrl file are ignored')
        if (nout == 1) then
          call sanrg(io_help==0,i,0,2,'rdctrl','HAM_AUTOBAS_LOC')
          autob = autob + 10*(getdig(i,0,10)-getdig(autob,1,10))
        endif

C       if (prgnam=='LMF'. or. prgnam=='LMFA') then
        nm='HAM_AUTOBAS_PFLOAT';call gtv(trim(nm),tksw2(prgn,nm),
     .    ivec(1:2),nout=nout,def_i4v=iapflt,note=
     .    '&Governs the MT log derivative parameters P'//
     .    '%N%3f1s digit, 1st element, affects how P is floated'//
     .    '%N%5f0 Use version 6 bounds on P and initial defaults'//
     .    '%N%5f1 Use bounds on P and initial defaults designed for LDA'//
     .    '%N%5f2 Use bounds on P and initial defaults designed for GW'//
     .    '%N%3f10s digit, 1st element, controls hows PZ is floated'//
     .    '%N%5f0 PZ (low) is allowed to float and is spin-dependent'//
     .    '%N%5f1 PZ (low) is allowed to float and is spin-averaged'//
     .    '%N%5f2 PZ is frozen'//
     .    '%N%5f4 add to maintain pre-7.14 compatibility in semicore envelope parms (nsp=2)'//
     .    '%N%3f2nd element affects how band CG is determined, used to float P'//
     .    '%N%5f0 Band CG found by traditional method'//
     .    '%N%5f1 Band CG found from true energy moment of density')
        if (nout >= 1) then
          pnudef = ivec(1)
          call sanrg(io_help==0,mod(pnudef,10),0,2,'rdctrl',
     .      '1st argument (1s digit) HAM_AUTOBAS_PFLOAT')
C          call sanrg(io_help==0,mod(pnudef/10,10),0,2,'rdctrl',
C     .      '1st argument (10s digit) HAM_AUTOBAS_PFLOAT')
          lekkl = ivec(2)
          call sanrg(io_help==0,lekkl,0,1,'rdctrl','2nd argument HAM_AUTOBAS_PFLOAT')
        endif


C        nm='HAM_AUTOBAS_APW';
C        call gtv(trim(nm),tksw2(prgn,nm),autob,Texist=ltmp,
C     .  nout=nout,def_i4=0,note='Autoset basis:  '//
C     .  'controls size of APW basis'//
C     .  '%N%4f0 null   APW basis (default PWMODE->0)'//
C     .  '%N%4f1 small  APW basis (default PWEMAX->3)'//
C     .  '%N%4f2 medium APW basis (default PWEMAX->5)'//
C     .  '%N%4f3 large  APW basis (default PWEMAX->8)'//
C     .  '%N%4f4 larger APW basis (default PWEMAX->11)'//
C     .  '%N%3f10s digit refers to APW basis type '//
C     .  '%N%4f0 fixed APW basis (default PWMODE->1)'//
C     .  '%N%4f1 q-dependent APW basis (default PWMODE->11)')
      endif
      if (io_help == 0) then
        call sanrg(io_help==0,mod(autob,10),0,4,'rdctrl','HAM_AUTOBAS_MTO')
      endif
      basopt(1) = autob
C     basopt(2) = lmxb_cut(1)
C     basopt(3) = lmxb_cut(2)

C ... Switch to turn on APW basis ... set i = default value
      nm='HAM_PWMODE'; sw=tksw2(prgn,nm); call gtv(trim(nm),sw,pwmode,def_i4=0,note=
     .  '&Controls APW addition to LMTO basis'//
     .  '%N%3f1s digit:'//
     .  '%N%6f0: LMTO basis only'//
     .  '%N%6f1: Mixed LMTO+PW'//
     .  '%N%6f2: PW basis only'//
     .  '%N%3f10s digit:'//
     .  '%N%6f0: PW basis fixed'//
     .  '%N%6f1: PW basis q-dependent')

      call gtv('HAM_PWEMIN',sw,pwemin,def_r8=0d0,nout=nout,
     .  note='Include APWs with energy E > PWEMIN (Ry)')

C ... APW basis cutoff ... set i = default value
      xx = 3
      if (ltmp) then
        j = mod(autob/10,10)
        if (j == 0) xx = 0
        if (j == 1) xx = 3
        if (j == 2) xx = 5
        if (j == 3) xx = 8
        if (j == 4) xx =11
      endif
      call gtv('HAM_PWEMAX',sw,pwemax,def_r8=xx,nout=nout,
     .  note='&Include APWs with energy E < PWEMAX (Ry)')
      call gtv('HAM_NPWPAD',sw,npwpad,def_i4=-1,
     .  note='Overrides default padding of variable basis dimension')

      nm='HAM_EBAS'; call gtv(trim(nm),tksw2(prgn,nm),
     .  kmto,nout=nmto,note='LMTO envelope kinetic energies')
      nmto = max(nmto,0)
      nm='HAM_DABC' ; call gtv(trim(nm),tksw2(prgn,nm),
     .  dabc,nout=nout,
     .  note='Spacings for real-space interstital mesh')
      if (nout == 1) dabc(2) = dabc(1)
      if (nout <  3) dabc(3) = dabc(2)
!!    Merge with ?
      nm='HAM_DQVAL'; call gtv(trim(nm),tksw2(prgn,nm),dqval,note='Total charge')
      nm='HAM_RDSIG'; sw=tksw2(prgn,nm); call gtv(trim(nm),sw,lrsig,def_i4=0,note=
     .  'Controls how self-energy is added to '//
     .  'local exchange correlation functional:'//
     .  '%N%3f 1s digit:'//
     .  '%N%3f   0 do not read Sigma'//
     .  '%N%3f   1 add sigma to potential'//
     .  '%N%3f   Add 2 to symmetrize Sigma(T)'//
     .  '%N%3f   Add 4 to include Re(Sigma(T)) only'//
     .  '%N%3f 10s digit:'//
     .  '%N%3f   0 simple interpolation'//
     .  '%N%3f   1 approx high and low sigma by diagonal (see sigp)'//
     .  '%N%3f   3 interpolate sigma by LDA evecs'//
     .  '%N%3f     (In this mode use 100s digit for'//
     .  ' # interpolation points)'//
     .  '%N%3f Add 10000 to indicate sigma has no symops'//
     .  '%N%3f Add 20000 to use minimum neighbor table'//
     .  '%N%3f Add 40000 to allow file qp mismatch')

      if (io_help==0) then
        xx = 20*avwsr(plat,1d0,xv,nbas)
      else
        xx = NULLR
      endif
      call gtv('HAM_RSRNGE',sw,rsrnge,def_r8=xx,
     .  note='Maximum range in connecting vectors for r.s. sigma (units of alat)')
      call gtv('HAM_RSSTOL',sw,rsstol,def_r8=5d-6,
     .  note='Max tolerance in Bloch sum error for r.s. sigma ')


      nm='HAM_NMTO'; call gtv(trim(nm),tksw2(prgn,nm),nmto,def_i4=0,Texist=ltmp,
     .  note='Order of polynomial approximation for NMTO hamiltonian')
      call sanrg(io_help==0.and.nmto>0,nmto,2,5,'rdctrl','nmto')
      if (nmto > 0 .or. (io_help /= 0 .and. ltmp)) then
!       Warning: kmto also used by lmmc; meaning is somewhat different
        nm='HAM_KMTO'; call gtv(trim(nm),tksw2(prgn,nm),
     .    kmto(1:max(1,nmto)),nmin=nmto,
     .    note='Corresponding NMTO kinetic energies'//
     .    '%N%3fRead NMTO values, or skip if NMTO=0' )
      endif
      nm='HAM_EWALD'; call gtv(trim(nm),tksw2(prgn,nm),ham_ewald,
     .  def_lg=.false.,note='Make strux by Ewald summation')
      nm='HAM_VMTZ'; call gtv(trim(nm),tksw2(prgn,nm),vmtz,def_r8=0d0,
     .  note='Muffin-tin zero defining wave functions')
      nm='HAM_QASA'; call gtv(trim(nm),tksw2(prgn,nm),ham_qasa,
     .  def_i4=3,note='0 => Methfessel conventions for 2nd gen'//
     .  ' ASA moments Q0..2'//
     .  '%N%3f1 => Q2 = coff to phidot**2 - p phi**2 in sphere'//
     .  '%N%3f2 => Q1,Q2 accumulated as coffs to <phi*phidot> '//
     .                 'and <phidot**2>'//
     .  '%N%3f3 => 1+2 (Stuttgart conventions)'//
     .  '%N%3f4 => Add 4 to compute phi,phidot integrating only outward')

      nm='HAM_PMIN'; call gtv(trim(nm),tksw2(prgn,nm),pmin,
     .  def_r8v=zerov,nout=nout,note=
     .  'Global minimum in fractional part of P-functions, all species.'//
     .  '%N%3fEnter positive values for l=0..nl:'//
     .  '%N%3f0: no minimum constraint'//
     .  '%N%3f#: with #<1, floor of fractional P is #'//
     .  '%N%3f1: use free-electron value as minimum'//
     .  '%N%3fEnter one negative value, applies to all l'//
     .  '%N%2f-1: lower bound P set to free electron value'//
     .  '%N%2f-2: lower bound on P designed for LDA (like AUTOBAS_PFLOAT=1)'//
     .  '%N%2f-3: lower bound on P designed for GW  (like AUTOBAS_PFLOAT=2)')

      if (pmin(1) < 0) then
        i = -pmin(1)
        if (-i /= pmin(1)) call rx('pmin must be a negative integer')
        call sanrg(io_help==0,-i,-3,-1,'rdctrl','pmin')
        j = 1+mod(i+1,3)
        call defpq(20+j,0d0,n0-1,1,pmin,0)
        forall (j=1:n0) pmin(j) = mod(pmin(j),1d0)
      endif

      nm='HAM_PMAX'; call gtv(trim(nm),tksw2(prgn,nm),
     .  pmax, def_r8v=zerov, nout=nout, note=
     .  'Global maximum in fractional part of P-functions.'//
     .  '%N%3fEnter values for l=0..nl:'//
     .  '%N%3f0: no maximum constraint'//
     .  '%N%3f#: with #<1, ceiling of fractional P is #')

      if (prgnam == 'LM' .or. prgnam == 'LMGF') then
        a = '%N%3f1: Read shift in ASA PP''s'
      endif
      nm='HAM_RDVEXT'; call gtv(trim(nm),tksw2(prgn,nm),
     .  ham_rdvext, def_i4=0, note=
     .  'Read parameters for adding external potential, file vext'//
     .  '%N%3f0: nothing read'//trim(a))

      nm='HAM_ALFSI'; call gtv(trim(nm),tksw2(prgn,nm),alfsi,def_r8=0d0,note=
     .  'Coefficient to artificial addition of overlap to hamiltonian')
      nm='HAM_OVEPS'; call gtv(trim(nm),tksw2(prgn,nm),
     .  oveps, def_r8=0d0, nout=nout, note=
     .  'Diagonalize hamiltonian in reduced hilbert space,'//
     .  '%N%3fdiscarding part with evals of overlap < OVEPS')
      nm='HAM_OVNCUT'; call gtv(trim(nm),tksw2(prgn,nm),
     .  ovncut, def_i4=0, nout=nout, note=
     .  'Diagonalize hamiltonian in reduced hilbert space,'//
     .  '%N%3fdiscarding space belonging to lowest OVNCUT evals'//
     .  ' of overlap.  Supersedes OVEPS')

      nm='HAM_SX'; call gtv(trim(nm),tksw2(prgn,nm),lsx,
     .  def_i4=0,note='Screened exchange:'//
     .  '%N%3f 0 do nothing'//
     .  '%N%3f 1 Calculate SX Sigma'//
     .  '%N%3f11 Calculate SX Sigma w/ local W')
      nm='HAM_SXOPTS'; call gtv(trim(nm),tksw2(prgn,nm),sxopt,nmin=10,
     .  note='Options for screened exchange,.e.g.  SXOPTS=rssig;nit=3')

      sw = tksw2(prgn,'HAM_SIGP')
C     This call does nothing except for printout
      if (sw < 2) then
      nm='HAM_SIGP'; call gtv(trim(nm),sw,nono,Texist=ltmp,
     .  note= 'Parameters for replacing high (and possibly low) '//
     .  'subblocks of sigma'//
     .  '%N%3fwith a linear function of LDA eigenvalues (diagonal)')
      nm='HAM_SIGP_MODE'; call gtv(trim(nm),sw,sigp_mode,def_i4=4,
     .  nout=nout,
     .  note= 'Specifies linear fitting function'//
     .  '%N%3f 0 sigii > a + b*e'//
     .  '%N%3f 1 sigii = a + b*e'//
     .  '%N%3f 2 sigii > a  and  sigii < b'//
     .  '%N%3f 3 sigii = a + b*e'//
     .  '%N%3f 4 sigii supplied by sigm file (constant)')
      nm='HAM_SIGP_A'; call gtv(trim(nm),sw,sigp_a,def_r8=0.02d0,
     .  note='a parameter in fit')
      nm='HAM_SIGP_B'; call gtv(trim(nm),sw,sigp_b,def_r8=0.06d0,
     .  note='b parameter in fit')
      nm='HAM_SIGP_NMIN'; call gtv(trim(nm),sw,sigp_nmin,def_i4=0,
     . note='Replace subblock 1..nmin with fit')
      nm='HAM_SIGP_EMIN'; call gtv(trim(nm),sw,sigp_emin,def_r8=0d0,
     . note='Replace subblock elda<emin with fit (not used if nmin>=0)')
      nm='HAM_SIGP_NMAX'; call gtv(trim(nm),sw,sigp_nmax,def_i4=0,
     . note='Replace subblock nmax+1.. with fit')
      nm='HAM_SIGP_EMAX'; call gtv(trim(nm),sw,sigp_emax,def_r8=2d0,
     . note='Replace subblock elda>emax with fit (not used if nmax>0)')
      nm='HAM_SIGP_EFIT'; call gtv(trim(nm),sw,sigp_efit,def_r8=0d0,
     .  note='Energy cutoff for fitting sigii (not used in any calculation')
      endif

      nm='HAM_UDIAG'; call gtv(trim(nm),tksw2(prgn,nm),ham_udiag,
     .  def_i4=0,note='nonzero => diagonal-only LDA+U')

C      nm='HAM_SITED'
C      sw=tksw2(prgn,nm);call gtv(trim(nm),sw,ham_blockd,nmin=10,
C     .  note='Render blocks of the Hamiltonian diagonal.  Syntax:'//
C     .  '%N%3fBLOCKD=integer-list[;E=#]'//
C     .  '%N%3fOptional E=# sets the diagonal part of H to #'//
C     .  '%N%3fSee doc/Integer-list-syntax.html for the syntax of '//
C     .  'integer-lists.')

      endif ! HAM

C --- Symmetry group ---
      if (io_show+io_help/=0 .and. tkswp(prgn,express,'SYMGRP') < 2)
     .  call info0(2,1,0,' --- Symmetry group operations ---')
      nm='SYMGRP'; call gtv(trim(nm),tksw2(prgn,nm),symg,
     .  note='String containing generators for symmetry group'//
     .  '%N%3fInclude "SOC=1" for reduced magnetic symmetry'//
     .  '%N%3fInclude "NOSPIN" to exclude spinor rotation when symmetrizing over k'//
     .  '%N%3fInclude "NOORB" to exclude orbital rotation when performing spinor rotations')

C --- Crystal Green's function ---
      if (tksw2(prgn,'GF') >= 2) goto 269
      if (io_show+io_help/=0) call info0(2,1,
     .  0,' --- Parameters for GF ---')
      nm='GF_MODE'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lcgf,
     .  def_i4=0,note=' 0: do nothing'//
     . '%N%3f 1: self-consistent cycle'//
     .'%N%3f10: Transverse exchange interactions J(q), MST'//
     .'%N%3f11: Read J(q) from disk and print derivative properties'//
C    .'%N%3f12: generate G(R,R'')'//
     .'%N%3f14: Longitudinal exchange interactions J(q), MST'//
     .'%N%3f20: Transverse chi+- from ASA GF'//
     .'%N%3f21: Read chi from disk and print derivative properties'//
     .'%N%3f24: Transverse chi++,chi-- from ASA GF')

      nm='GF_GFOPTS'; call gtv(trim(nm),tksw2(prgn,nm),gfopt,nmin=10,
     .  note='Switches governing execution of GF')

      nm='GF_DLM'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_ldlm,def_i4=0,
     .  note='Switch for chemical CPA and DLM'//
     .  '%N%3f12: normal CPA/DLM calculation: '//
     .  'charge and omega both iterated'//
     .  '%N%3f32: only omega iterated to '//
     .  'prescribed tolerance for each z'//
     .  '%N%2f112:  not documented.')
      if (ctrl_ldlm /= 0) lham512 = .true.

C      nm='GF_BXY'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lbxc,def_i4=0,
C     .  note='(DLM) Constraining fields are converged if set to 1')
      nm='GF_BXY'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lbxc4,def_lg=F,
     .note='(DLM) If T, a site-dependent constraining field is added'//
     .  '%N%3fto properly align magnetic moments.'//
     .  '%N%3fSee HAM_BXCSCAL and documentation in tokens.html')
      if (ctrl_lbxc4) ctrl_lbxc1 = 2
      nm='GF_TEMP'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tdlm,
     .  def_r8=0d0,note='(DLM) temperature parameter in p(theta)')

  269 continue

C --- Planar Green's function ---
      if (tkswp(prgn,express,'PGF') >= 2) goto 219
      if (io_show+io_help/=0) call info0(2,1,
     .  0,' --- Parameters for PGF ---')
C     nm='PGF'; call gtv(trim(nm),tksw2(prgn,nm),Texist=ltmp,note=
C    .  'Parameters for planar Green''s function')
      ctrl_lpgf(1) = 0
      nm='PGF_MODE'; call gtv(trim(nm),tksw2(prgn,nm),
     .  ctrl_lpgf(1),note='0: do nothing'//
     .'%N%3f1: diagonal layer GF'//
     .'%N%3f2: left- and right- bulk GF'//
     .'%N%3f3: find k(E) for left bulk'//
     .'%N%3f4: find k(E) for right bulk'//
     .'%N%3f5: Calculate conductance'//
     .'%N%3f7: Calculate reflectance'//
     .'%N%3f8: Calculate spin torques')
      if (io_help == 0) then
        if (ctrl_lpgf(1) == NULLI) ctrl_lpgf(1) = 0
        if (ctrl_lpgf(1)==0) goto 219
      else
      if (io_show+io_help/=0) call info0(2,0,0,
     .    ' ... The following are read if MODE>0')
      endif
      nm='PGF_SPARSE'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lpgf(2),
     .  def_i4=0,note='0: Calculate GF layer by layer'//
     .  '%N%3f1: Calculate GF with LU decomposition')
      nm='PGF_PLATL'; call gtv(trim(nm),tksw2(prgn,nm),
     .  platl, nmin=3,note='third lattice vector of left bulk')
      nm='PGF_PLATR'; call gtv(trim(nm),tksw2(prgn,nm),
     .  platr, nmin=3,note='third lattice vector of right bulk')
      gfopt = ' '
      nm='PGF_GFOPTS'; call gtv(trim(nm),tksw2(prgn,nm),gfopt,nmin=10,
     .  note='Switches governing execution of PGF')
  219 continue

C --- Species ---
      if (tkswp(prgn,express,'SPEC') >= 2) goto 79
      if (io_help == 1) nspec = 1
      if (io_show+io_help/=0) call info0(2,1,0,
     .  ' --- Parameters for species data ---')
C     SPEC call is only helpful for printout
C      if (io_help /= 0 .or. io_show>0) write(stdo,'(1x)')
C      nm='SPEC'; call gtv(trim(nm),tksw2(prgn,nm),Texist=ltmp,note=
C     .  'Parameters for species data')
C      if (.not. ltmp) goto 79
      if (io_help /= 0) call info0(2,0,0,' ... The next four tokens '//
     .  'apply to the automatic sphere resizer')
      nm='SPEC_SCLWSR'; call gtv(trim(nm),tksw2(prgn,nm),
     .  sclwsr, def_r8=0d0, note=
     .  'Scales sphere radii, trying to reach volume = '//
     .  'SCLWSR * cell volume'//
     .  '%N%3fSCLWSR=0 turns off this option.'//
     .  '%N%3fAdd  10  to initially scale non-ES first;'//
     .  '%N%3f or  20  to scale ES independently.')
C     if (sclwsr /= 0) then
      xv(1:3) = (/.16d0,.18d0,.20d0/)
      if (lfp /= 0) xv(1:3) = 0d0
      nm='SPEC_OMAX1'; call gtv(trim(nm),tksw2(prgn,nm),omax1,
     .  def_r8v=xv(1:3),note=
     .  'Limits max sphere overlaps when adjusting MT radii')
C     xv(1:3) = (/.40d0,.45d0,.50d0/)
      xv(1:3) = (/.40d0,1d0,1d0/)
      if (lfp /= 0) xv(1:3) = 0d0
      nm='SPEC_OMAX2'; call gtv(trim(nm),tksw2(prgn,nm),
     .  omax2, def_r8v=xv(1:3),note=
     .  'Sphere overlap constraints of second type',nout=nout)
      nm='SPEC_WSRMAX'; call gtv(trim(nm),tksw2(prgn,nm),wsrmax,
     .  def_r8=0d0,note=
     .  'If WSRMAX is nonzero, no sphere radius may exceed its value')
C     endif

      nm='SPEC_HSFITK'; call gtv(trim(nm),tksw2(prgn,nm),hsfitd,
     .  def_r8v=(/1.5d0,1d0/),note='HSFITK*HCR = radius'//
     .  ' beyond which K.E. should approach asymptotic value.'//
     .  '%N%16fIndividual species not explicitly set with '//
     .  'SPEC_ATOM_HSFITK'//
     .  '%N%16ftake result of SPEC_HSFITK as a default')

C       For printout only
        nm='SPEC_HFAC'; call gtv(trim(nm),tksw2(prgn,nm),
     .  alabl,nmin=0,Texist=ltmp,note='Supplies default scale factors'//
     .  ' that can parameterize the LDA hamiltonian.'//
     .  '%N%3fThese values only supply input to file shfac '//
     .        'when it does not yet exist.'//
     .  '%N%3fIf shfac already exists these parameters are not used.'//
     .  '%N%3fFor species-specific factors, edit shfac or '//
     .  'use SPEC_ATOM_HFAC[B,L,V].')
        if (express <= 0) hfacd = 1 ! initialize once only
        nm='SPEC_HFAC_B'; call gtv(trim(nm),tksw2(prgn,nm),
     .    hfacd(1),def_r8=1d0,note='Scale BXC by B.'//
     .    '  To invoke this scaling, also turn on HAM_BXCSCAL.'//
     .    '%N%3fUse SPEC_ATOM_HFACB for species-specific scalings.')
C        nm='SPEC_HFAC_L'; call gtv(trim(nm),tksw2(prgn,nm),
C     .    hfacd(2),def_r8=1d0,note='L.S scaling.'//
C     .    '  To invoke this scaling, also turn on HAM_SOSCAL.'//
C     .    '%N%3fUse SPEC_ATOM_HFACL for species-specific scalings.')
        nm='SPEC_HFAC_D'; call gtv(trim(nm),tksw2(prgn,nm),
     .    hfacd(3),def_r8=1d0,note='Scale channel bandwidth by D')
        nm='SPEC_HFAC_V'; call gtv(trim(nm),tksw2(prgn,nm),
     .    lvshft,def_lg=F,note='potential shifts at sites.')

      if (io_help >= 1) then
        write(*,382)
  382   format(/' - ',
     .    'The following tokens are read for each species. ',
     .    'Data sandwiched'/3x,'between successive occurences of ',
     .    'token ATOM apply to one species.')
        nspec = 1
      endif

      if (nspec <= 0) goto 79

      if (.not. allocated(pnu)) then
      allocate(pnu(n0,nsp,nspec),qnu(n0,nsp,nspec),
     .  pz(n0,nsp,nspec),amom(n0,nspec),idmod(n0,nspec),
     .  rsmh(n0,nspec),eh(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec),
     .  pb1(nspec),pb2(nspec),lmxpb(nspec),
     .  ehvl(n0,nspec),
     .  qpol(n0,nspec),stni(nspec),tbvso(4,nspec),
     .  iq1(n0,nspec),iq2(n0,nspec),
     .  rg(nspec),rsma(nspec),rfoca(nspec),rsmfa(nspec),rcfa(2,nspec),
     .  rs3(nspec),rham(nspec),rmt(nspec),rsmv(nspec),
     .  nxi(nspec),exi(n0,nspec),rint(nspec),rcut(nspec),
     .  spec_a(nspec),z(nspec),nr(nspec),mass(nspec),eref(nspec),
     .  coreh(nspec),coreq(2,nspec),
     .  colxbs(3,nspec),radxbs(nspec),hsclsm(nspec),hsfitp(2,nspec),
     .  idxdn(n0,nspec),
     .  hcr(n0,nspec),rsminl(n0,nspec),alpha(n0,nspec),
     .  idu(4,nspec),uh(4,nspec),jh(4,nspec),
     .  dv(nspec),grp(nspec),grp2(nspec),
     .  mxcst1(nspec),mxcst2(nspec),mxcst4(nspec),
     .  kmxt(nspec),kmxv(nspec),
     .  lfoca(nspec),lmxl(nspec),lxi(nspec),lmxa(nspec),lmxb(nspec),
     .  nthet(nspec),nang(2,nspec),ncpa(nspec),nbeff(nspec),
     .  iscpa(n0,nspec),xcpa(n0,nspec),beff(n0,nspec),hfacs(3,nspec))
!     allocate(idm(4,nspec))
      endif

      lmxa = 0
c      jc = 0
c      if (optio == 2) nkaph = nglob('nkaph')
C     Loop over species
      nkaph = 1
      lpzi = 1
      qpol = NULLR
      rsmh = 0d0
      rsmh2 = 0d0
      eh = NULLR
      eh2 = NULLR
      hcr = NULLR
      idmod = NULLI
      nthet = NULLI
      nang = NULLI
      ehvl = NULLR
      dv = 0
      rsminl = NULLR
      j = 0; jp = 0
C     do  j = 1, nspec
      do while (j < nspec)  ! Until parse and retain nspec
        jp = jp+1  ! jp = index to current species in file
        j = j+1    ! j  = index to current species to keep
        colxbs(:,j) = NULLR; radxbs(j) = NULLR
        rcfa(:,j) = NULLR; rfoca(j) = 0d0; rg(j) = 0d0
        rham(j) = NULLR; rsma(j) = 0d0; rsmfa(j) = 0d0
        spec_a(j) = NULLR; nr(j) = NULLI
        exi(:,j) = NULLR
        coreh(j) = ' '; coreq(:,j) = NULLR
        eref(j) = NULLR
        mass(j) = NULLR

        if (io_help /= 0) then
          write(stdo,'(1x)')
        elseif (io_help == 0 .and. io_show>0) then
          call info(1,0,0,' ... Species %i',j,0)
        endif

        jj = (/1,jp/)

        nm='SPEC_ATOM'; call gtv(trim(nm),tksw2(prgn,nm),slabl(j),
     .    nmin=10,cindx=jj,note='Species label')

        if (.not. associated(slabll,slabl)) then
          call tokmat(slabl(j),slabll,nspec,8,' ',i0,nn,.false.)
          i0 = i0+1
          if (i0 <= 0) then
            call info0(10,0,0,' ... Discard unused species '//slabl(j))
            j = j-1             ! Discard species
            cycle
          endif
        endif

        nm='SPEC_ATOM_Z'; call gtv(trim(nm),tksw2(prgn,nm),z(j),
     .    cindx=jj,note='Atomic number')
        sw = tksw2(prgn,'SPEC_ATOM_R')
        if (sw < 2) then
          nout = 0
          nm='SPEC_ATOM_R'; call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .      nout=nout,note= 'Augmentation sphere radius rmax',or=T)
          if (nout /= 1) then
          nm='SPEC_ATOM_R/W';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .      nout=nout,note='rmax relative to average WS radius',or=T)
          if (nout == 1) rmt(j) =rmt(j)*avw
          if (nout /= 1) then
          nm='SPEC_ATOM_R/A';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .      nout=nout,note='rmax relative to lattice constant')
          if (nout == 1) rmt(j) =rmt(j)*alat
          endif
          endif
        endif
        if (ltbe) rmt(j) = 1d0

C   ... Radial mesh parameters: determine default value of a
        i0 = NULLI; xx = NULLR
        if (io_help == 0) then
          call pshpr(1)
          call rmesh(z(j),rmt(j),mod(lrel,10),.false.,nrmx,xx,i0)
          call poppr
          if (xx == .03d0) xx = .025d0
        endif
        nm='SPEC_ATOM_A'; call gtv(trim(nm),tksw2(prgn,nm),spec_a(j),
     .    def_r8=xx,cindx=jj,nout=nout,
     .    note='&Radial mesh point spacing parameter')
C       Determine default NR
        if (tksw2(prgn,'SPEC_ATOM_NR') /= 2) then
          i0 = 0
          call pshpr(1)
          call rmesh(z(j),rmt(j),mod(lrel,10),.false.,nrmx,spec_a(j),i0)
          call poppr
        endif
        nm='SPEC_ATOM_NR'; call gtv(trim(nm),tksw2(prgn,nm),nr(j),
     .    def_i4=i0,cindx=jj, note='Number of radial mesh points')
        if (nr(j) == 0) nr(j) = i0
C       Case lmxb will be set by program: lmxa be at least that large
C        i = nl-1
C        if (prgnam=='LMFA') then
C        if (mod(autob,10) /= 0 .or. cmdopt('--basp',6,0,nm)) then
C          call fadflb(int(basopt(10)),99,z(j),ivec,ivec(3))
C          i = max(ivec(1),ivec(2),nl-1)
C        endif
C        endif
        nm='SPEC_ATOM_LMX'; call gtv(trim(nm),tksw2(prgn,nm),lmxb(j),
     .    def_i4=max(nl-1,NULLI),Texist=ltmp,cindx=jj,
     .    note='&l-cutoff for basis')
        ldefaultl = .not. ltmp   ! Not read; default was used

C       Running account of maximum lmxb
        lmxbj = NULLI
        if (io_help == 0) then
        if (prgnam == 'LMMC') then
          call lx2vec(lmxb(j),0,nn,ivec)
          do  i = 1, nn
            lmxbj = max(lmxbj,ivec(i))
          enddo
        else
          lmxbj = lmxb(j)
        endif
        endif
C       nlbj = number of elements associated with lmxb
C              0 => no elements
C       nlbji: ditto, but used to specify number of default values
        nlbj = 1+lmxbj
        nlbji = nlbj
        if (io_help >= 1) then
          nlbji = NULLI
          nlbj = 1
        elseif (lmxbj == NULLI) then
          nlbji = NULLI
          nlbj = 0
        endif

C   ... Basis set for lmf
        xx = 0; if (prgnam=='LMFA') xx = -1; call dvset(xv,1,n0,xx)
        nm='SPEC_ATOM_RSMH'; call gtv(trim(nm),tksw2(prgn,nm),
     .    rsmh(1:nlbj,j),cindx=jj,nout=nout,nmin=nlbj,Texist=ltmp,
     .    def_r8v=xv,note='Smoothing radii for basis')
        ldefaultl = ldefaultl .and. ltmp  ! If RSMH was read but not lmx
                                          ! lmx is specified by it
C        nn = NULLI; if (nout > 0) nn = nout
CC       Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
C        do  i = 1, nlbj
C          if (rsmh(i,j) <= 0) cycle
C          nn = i
C        enddo

        nm='SPEC_ATOM_EH'; call gtv(trim(nm),tksw2(prgn,nm),eh(1:nlbj,j),
     .    nout=nout,cindx=jj,note='Kinetic energies for basis')
        if (nout > 0) then
          forall (i=nout+1:nlbj) eh(i,j) = eh(nout,j) ! Fill in missing EH with last entry
        endif

        xx = 0; if (prgnam=='LMFA') xx = -1; call dvset(xv,1,n0,xx)
        nm='SPEC_ATOM_RSMH2'; call gtv(trim(nm),tksw2(prgn,nm),
     .    rsmh2(1:nlbj,j),nmin=nlbj,def_r8v=xv,cindx=jj,nout=nout,
     .    Texist=ltmp,note='Basis smoothing radii, second group')
        if (ltmp) then
C         nn = NULLI
          sw = tksw2(prgn,nm)
          if (nout>0) then
            nkaph=2; sw = 1
C           nn = nout
C           Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
C            do  i = 1, nlbj
C              if (rsmh2(i,j) == 0) cycle
C              nn = i
C            enddo
          endif
          nm='SPEC_ATOM_EH2'; call gtv(trim(nm),sw,eh2(1:nlbj,j),
     .      nout=nout,cindx=jj,note='Basis kinetic energies, second group')
          if (nout > 0) then
            forall (i=nout+1:nlbj) eh2(i,j) = eh2(nout,j) ! Fill in missing EH2 with last entry
          endif
        endif

        nm='SPEC_ATOM_EHVL'; call gtv(trim(nm),tksw2(prgn,nm),
     .    ehvl(1:nlbj,j),cindx=jj,def_r8v=(/(-0.5d0,i=1,n0)/),
     .    nout=nout,nmin=nlbj, note='val-lap fit energies')

C   ... If ldefaultl wasn't specified but RSMH was
        if ((prgnam=='LMFA' .or. prgnam=='LMF') .and. ldefaultl .and. io_help==0) then
          nn = 0
          do  i = nlbj, 1, -1
            nn = i
            if (rsmh(i,j) > 0 .or. rsmh2(i,j) > 0) exit
          enddo
          if (io_show /= 0 .and. lmxb(j)+1 /= nn) then
            call info2(2,0,0,' ... Reducing LMX to %i',nn-1,0)
          endif
          lmxb(j) = min(lmxb(j),nn-1)
          lmxbj = lmxb(j)
        endif
C      lmxbx = max(lmxbx,lmxb(j))

C   ... Determine lmxa: floating orbitals have no lmxa
        nm='SPEC_ATOM_LMXA'; sw = tksw2(prgn,nm)
        if (rmt(j) == 0 .and. io_help /= 1) then
          lmxa(j) = -1
        elseif (sw >= 2) then   !lmxa not read: look for substitute
          if (lfp > 0) then     !Possibly replace with a formula
                                !that depends on rmt
            lmxa(j) = 4
          elseif (tksw2(prgn,'SPEC_ATOM_LMX') /= 2) then
            lmxa(j) = lmxb(j)
          elseif (tksw2(prgn,'STRUC_NL') /= 2) then
            lmxa(j) = nl-1
          else
            lmxa(j) = 0
          endif
        else
C         Possibly replace with a formula that depends on rmt
          call gtv(trim(nm),sw,lmxa(j),
     .      def_i4=max(nl-1,lmxbj,NULLI),cindx=jj,Texist=ltmp,note=
     .      '&l-cutoff for augmentation')
C         lmxb may not exceed lmxa
          if (io_help == 0 .and. lmxbj > lmxa(j))
     .      call rx2('species '//trim(slabl(j))//' : LMX=%i '//
     .      'exceeds LMXA=%i.  Revise input so that LMX<=LMXA.',
     .      lmxbj,lmxa(j))
        endif
        lmxaj = lmxa(j)
C       nlaj = number of elements associated with lmxa
C              0 => no elements
C       nlaji: ditto, but used to specify number of default values
        nlaj = 1+lmxaj
        nlaji = nlaj
        if (io_help >= 1) then
          nlaji = NULLI
          nlaj = 1
        elseif (lmxaj == NULLI) then
          nlaji = NULLI
          nlaj = 0
        endif
C        if (lmxaj == NULLI) nlaj = NULLI
C        if (lmxaj == NULLI) nlaji = NULLI

        kmxv(j) = 15
        nm='SPEC_ATOM_KMXV';call gtv(trim(nm),tksw2(prgn,nm),kmxv(j),
     .    def_i4=15,cindx=jj,note='cutoff to expand smoothed potential')

C   ... Parameters that depend on the existence of an augmentation sphere
C       lmxl = l-cutoff for numerical rep'sn of density in sphere
C              in TBE, l-cutoff for point multipoles and Madelung potential
C       lfoca = mode for treating core
C       kmxt = kmax for expansion of envelope wf tails
C       kmxv = cutoff to expand smoothed potential
        kmxv(j) = kmxvdefault
C       Cannot set default here: need set after rescaling of rmt
C       rsmv(j) = rmt(j)*.5d0   !Not input
        rsmv(j) = 0d0           !Not input
        kmxt(j) = -1            !If sought, default will be set below
        lfoca(j) = 0            !If sought, default will be reset below
        if (.not. ltbe)         !Use lmxaj in case not sought (ASA:mpol)
     .    lmxl(j) =  lmxaj      !In TBE, a default value of 2*(nl-1) is set in the gtv call
        lmxpb(j) = NULLI        !If sought, default will be set below
        lxi(j) = NULLI          !If sought, default will be set below
        nxi(j) = NULLI          !If sought, default will be set below
        pb1(j) = ' '
        pb2(j) = ' '
        call dpzero(pnu(1,1,j),n0*nsp)
        call dpzero(pz(1,1,j),n0*nsp)
        rint(j) = NULLI         !No default for now
        rcut(j) = NULLI         !No default for now
        rs3(j) = NULLI          !If sought, default will be set below
        idxdn(:,j) = 1
        pnu(1,1,j) = NULLI      !If sought, default will be set below
        qnu(1,1,j) = NULLI      !If sought, default will be set below
        mxcst1(j) = .false.
        mxcst2(j) = .false.
        mxcst4(j) = .false.
        idu(:,j) = 0
!       idm(:,j) = 0

C       Downfolding switches: if auto DNF turned on, default is zero.
        if (lham4) then
          nn = 0
C       Help mode: If ADNF COULD have been turned on, default unknown
        elseif (tksw2(prgn,'OPTIONS_ADNF')/= 2.and. io_help==1) then
          nn = NULLI
C       All other cases: default is 1
        else
          nn = 1
        endif
        nm='SPEC_ATOM_IDXDN'; call gtv(trim(nm),tksw2(prgn,nm),
     .    idxdn(1:nlbj,j),nmin=nlbj,
     .    def_i4v=(/(nn,i=1,nlbj)/),cindx=jj,
     .    note='downfolding index: 0, auto; 1, no dnf; 2, fold down;'//
     .    ' 3, neglect')

        if (ltbe) then
          nm='SPEC_ATOM_LMXL'; call gtv(trim(nm),tksw2(prgn,nm),lmxl(j),
     .        cindx=jj,def_i4=2*(nl-1),note=
     .        '&lmax up to which Q_L and V_L are constructed')
        endif

        if (nlaj /= 0) then
          if (.not. ltbe) then
            nm='SPEC_ATOM_LMXL'; call gtv(trim(nm),tksw2(prgn,nm),lmxl(j),
     .        cindx=jj,def_i4=lmxaj,note=
     .        '&lmax for which to accumulate rho,V in sphere')
          endif

C   ... Set up default P,Q in absence of explicit specification
        call dpzero(pnu(1,1,j),nsp*n0)
        call dpzero(qnu(1,1,j),nsp*n0)
C       This branch sets defaults beforehand ... too complicated
C        if (io_help == 0) then
CC         Default P,Q, nsp=1
C          if (debug) call pshpr(100)
C          call defpq(pnudef,z(j),lmxaj,1,pnu(1,1,j),qnu(1,1,j))
C          if (debug) call poppr
C        endif
C        nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw2(prgn,nm),
C     .    pnu(1:nlaj,1,j),def_r8v=pnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting log der. parameters for each l')
C        nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw2(prgn,nm),
C     .    qnu(1:nlaj,1,j),def_r8v=qnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting sphere charges for each l channel')
C        call snit
        nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw2(prgn,nm),
     .    pnu(1:nlaj,1,j),cindx=jj,note=
     .    'Starting log der. parameters for each l')
        nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw2(prgn,nm),
     .    qnu(1:nlaj,1,j),cindx=jj,note=
     .    'Starting sphere charges for each l channel')
C       Reset default P,Q in absence of explicit specification
        if (io_help == 0) then
          if (io_show /= 0) call pshpr(50)
          call defpq(mod(pnudef,10),z(j),lmxaj,1,pnu(1,1,j),qnu(1,1,j))
          if (io_show /= 0) call poppr
        endif
        if (nsp == 2) call dcopy(n0,pnu(1,1,j),1,pnu(1,2,j),1)
        if (nsp == 2 .or. io_help == 1) then
          nm='SPEC_ATOM_MMOM'; call gtv(trim(nm),tksw2(prgn,nm),
     .    qnu(1:nlaj,2,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .    'Starting mag. moms for each l channel')
        endif

        nm='SPEC_ATOM_PZ'; call gtv(trim(nm),tksw2(prgn,nm),
     .    pz(1:nlaj,1,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .    'Starting semicore log der. parameters'//
     .    '%N%10fAdd 10 to attach Hankel tail'//
     .    '%N%10fAdd 20 to include perturbatively',nout=nout)
        if (nout>0) then
          if (dasum(nlaj,pz(1,1,j),1) /= 0) lpzi = max(lpzi,2)
        endif

        i0 = 1
        if (z(j) <= 8) i0 = 0
        nm='SPEC_ATOM_LFOCA';call gtv(trim(nm),tksw2(prgn,nm),lfoca(j),
     .    def_i4=i0,cindx=jj,note='&FOCA switch 0, 1 or 2 (see docs)')
        i = 3
        if (pwmode /= 0) then
          i = max(i,int(pwemax+1))
        endif
        nm='SPEC_ATOM_KMXA'; call gtv(trim(nm),tksw2(prgn,nm),kmxt(j),def_i4=i,
     .    cindx=jj,note='&k-cutoff for projection of wave functions in sphere.')
        nm='SPEC_ATOM_KMXV'; call gtv(trim(nm),tksw2(prgn,nm),kmxv(j),def_i4=15,
     .    cindx=jj,note='&k-cutoff for overlap free atom density')
C       Cannot set default here: need set after rescaling of rmt
        nm='SPEC_ATOM_RSMA'; call gtv(trim(nm),tksw2(prgn,nm),rsma(j),
     .    def_r8=0d0,cindx=jj,note=
     .    '&Smoothing for projection of wave functions in sphere.'//
     .    '%N%3finput<0 => choose default * -input')
C       Cannot set default here: need set after rescaling of rmt
        nm='SPEC_ATOM_RSMG'; call gtv(trim(nm),tksw2(prgn,nm),rg(j),
     .    def_r8=0d0,cindx=jj,note=
     .    '&Smoothing for projection of charge in sphere.'//
     .    '%N%3finput<0 => choose default * -input')
C       Cannot set default here: need set after rescaling of rmt
        nm='SPEC_ATOM_RFOCA'; call gtv(trim(nm),tksw2(prgn,nm),rfoca(j),
     .    def_r8=0d0,cindx=jj,note=
     .   '&Smoothing for core tail.  input<0 => choose default * -input')
C       Cannot set default here: need set after rescaling of rmt
        nm='SPEC_ATOM_RSMFA'; call gtv(trim(nm),tksw2(prgn,nm),rsmfa(j),
     .    def_r8=0d0,cindx=jj,note=
     .   '&Smoothing for free atom.  input<0 => choose default * -input')
        nm='SPEC_ATOM_RCFA'; call gtv(trim(nm),tksw2(prgn,nm),rcfa(1:2,j),
     .    def_r8v=zerov,nmin=2,cindx=jj,note=
     .    'Cutoff radius for renormalization of free atom density.'//
     .    '%N%3fOptional 2nd argument = width'//
     .    '%N%3fRCFA<0 => renormalize potential instead of density')

C       Negative radii: convert to actual numbers
        if (rg(j) < 0) rg(j)    = -rg(j)*0.25d0*rmt(j)
        if (rsma(j) < 0) rsma(j)  = -rsma(j)*0.4d0*rmt(j)
        if (rfoca(j) < 0) rfoca(j) = -rfoca(j)*0.4d0*rmt(j)
        if (rsmfa(j) < 0) rsmfa(j) = -rsmfa(j)*0.5d0*rmt(j)

        nm='SPEC_ATOM_RS3'; call gtv(trim(nm),tksw2(prgn,nm),rs3(j),def_r8=1d0,
     .    cindx=jj, note='Lower bound to smoothing radius for local orbital')

C       ASA tight-binding alpha parameters
        call dpzero(xv,n0)
        lp1 = nlbj
        nm='SPEC_ATOM_ALPHA'; sw=tksw2(prgn,nm)
C       Set default values; poke into xv
        if (sw /= 2) then
C       Help mode, default -> unknown
        if (io_help == 1) then
          xv(1) = NULLI
        else
C         Determine effective 1+lmxb: neglected orbitals => reduce lmxb
          if (tksw2(prgn,'SPEC_ATOM_IDXDN') /= 2) then
   75       if (lp1 > 1 .and. idxdn(lp1,j) >= 4) then
              lp1 = lp1-1
              goto 75
            endif
          endif
          call dcopy(6,tbalph(1,lp1),1,xv,1)
        endif
        endif
        i0 = lp1
        if (lmxbj == NULLI) i0 = NULLI
C       Copy all 5 screening parameters, to avoid division zero
        call dcopy(6,tbalph(1,lp1),1,alpha(1,j),1)
        call gtv(trim(nm),sw,alpha(1:lp1,j),cindx=jj,nmin=i0,
     .   def_r8v=xv,note='Screening parameters for structure constants')

        nm='SPEC_ATOM_IDMOD'; call gtv(trim(nm),tksw2(prgn,nm),
     .    idmod(1:nlaj,j),nmin=nlaji,def_i4v=(/(0,i=1,n0)/),
     .    cindx=jj,note=
     .    'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu')

        nm='SPEC_ATOM_DV'; call gtv(trim(nm),tksw2(prgn,nm),dv(j),
     .    def_r8=0d0,cindx=jj,note='Artificial constant potential '//
     .    'shift added to spheres belonging to this species')

        nm='SPEC_ATOM_MIX'; call gtv(trim(nm),tksw2(prgn,nm),mxcst1(j),
     .    def_lg=F,cindx=jj,note='Set to suppress '//
     .    'self-consistency of classes in this spec')

        nm='SPEC_ATOM_CSTRMX'; call gtv(trim(nm),tksw2(prgn,nm),
     .    mxcst2(j),cindx=jj,def_lg=F,note='Set to exclude this'//
     .    ' species when automatically resizing sphere radii'//
     .    ' (SCLWSR>0)')
        if (sclwsr == 0) mxcst2(j) = F

        else
          call dcopy(5,tbalph(1,lmxbj+1),1,alpha(1,j),1)
        endif    ! end of input dependent on presence of aug sphere.

        nm='SPEC_ATOM_RHAM'; call gtv(trim(nm),tksw2(prgn,nm),rham(j),
     .    cindx=jj,note='Radius delimiting range of r.s. hamiltonian')

        nm='SPEC_ATOM_GROUP'; call gtv(trim(nm),tksw2(prgn,nm),grp(j),def_i4=0,cindx=jj)
        nm='SPEC_ATOM_GRP2'; call gtv(trim(nm),tksw2(prgn,nm),grp2(j),def_i4=0,note=
     .    'The l=0 density of species with a common positive value of GRP2 are symmetrized,'//
     .    '%N%3findependent of symmetry operations.'//
     .    '%N%3fGRP2<0: species with common |GRP2| are symmetrized '//
     .    '(but spins flipped for GRP2<0)',cindx=jj)
        nm='SPEC_ATOM_PBAS'; call gtv(trim(nm),tksw2(prgn,nm),pb1(j),
     .    nmin=10,cindx=jj,note='product basis for GW')
        nm='SPEC_ATOM_PBAS2'; call gtv(trim(nm),tksw2(prgn,nm),
     .    pb2(j),cindx=jj,note='second product basis for GW')
        nm='SPEC_ATOM_LMXPB'; call gtv(trim(nm),tksw2(prgn,nm),lmxpb(j),
     .    def_i4=4,cindx=jj,note='l-cutoff for product basis')

        sw = tksw2(prgn,'SPEC_ATOM_IDU')
        if (io_help > 0 .and. sw < 2) then
          call info0(2,0,0,' ... The next tokens apply to LDA+U or LDA+DMFT')
        endif
        nm='SPEC_ATOM_IDU'; call gtv(trim(nm),sw,idu(:,j),cindx=jj,def_i4v=(/(0,i=1,4)/),
     .    note='LDA+U mode:  0 not treated by LDA+U, 1 AMF, 2 FLL, 3 mixed '//
     .    '%N%3f10''s digit 0 => apply U to phi only'//
     .    '%N%3f10''s digit 1 => apply U to phi, phidot, LO')

C   ... to be removed in future
        nm='SPEC_ATOM_IDD'; call gtv(trim(nm),tksw2(prgn,'SPEC_ATOM_IDD'),ivec(1:4),
     .    cindx=jj,def_i4v=(/(0,i=1,4)/),Texist=ltmp,
     .    note='DMFT mode:  0 not treated by DMFT, 1 AMF, 2 FLL, 4 constrained occupancy')
        if (ltmp .or. io_help /= 0) then
          do  i = 1, 4
            if (ivec(i) /= 0 .and. idu(i,j) /= 0 .and. io_help == 0)
     .        call rx('m_rdctrl: cannot apply LDA+U and LDA+DMFT to the same orbital')
            idu(i,j) = idu(i,j) + 100*ivec(i)
          enddo
C          nm='SPEC_ATOM_IDM'; call gtv(trim(nm),tksw2(prgn,nm),idm(:,j),cindx=jj,def_i4v=(/(0,i=1,4)/),
C     .      note='DMFT mode:  One number for each l. Each number consists of a sequence of%N%3f'//
C     .      '(2l+1) digits prescribing m-dependent treatment of correlation, e.g. 11212')
        endif

C       Note sw is sw for SPEC_ATOM_IDU
        nm='SPEC_ATOM_UH'; call gtv(trim(nm),sw,uh(:,j),cindx=jj,
     .    def_r8v=zerov,note='Hubbard U for LDA+U or DMFT')
        nm='SPEC_ATOM_JH'; call gtv(trim(nm),sw,jh(:,j),cindx=jj,
     .    def_r8v=zerov,note='Exchange parameter J for LDA+U or DMFT')

C   ... lm-specific
        nm='SPEC_ATOM_FRZC'; call gtv(trim(nm),tksw2(prgn,nm),
     .    iq1(1:nlbj,j),cindx=jj,nmin=nlbji,def_i4v=izerv,
     .    note='1: freeze ASA C parameters (band fitting mode)')
        nm='SPEC_ATOM_FRZD'; call gtv(trim(nm),tksw2(prgn,nm),
     .    iq2(1:nlbj,j),cindx=jj,nmin=nlbji,def_i4v=izerv,
     .    note='1: freeze ASA Delta parameters (band fitting mode)')

C   ... DLM-specific
        nm='SPEC_ATOM_NTHET';call gtv(trim(nm),tksw2(prgn,nm),nthet(j),
     .    def_i4=0,cindx=jj,note='number of theta angles for DLM')

        nm='SPEC_ATOM_NANG';call gtv(trim(nm),tksw2(prgn,nm),nang(1:2,j),
     .    def_i4v=izerv,cindx=jj,note='number of solid angles for (relativistic) DLM')

        nm='SPEC_ATOM_BEFF'; call gtv(trim(nm),tksw2(prgn,nm),beff(:,j),
     .    def_r8v=zerov,cindx=jj,nout=nout,note='Effective field Legendre moments')
        nbeff(j) = nout

C   ... CPA components
        nm='SPEC_ATOM_CPA'; call gtv(trim(nm),tksw2(prgn,nm),iscpa(:,j),
     .    cindx=jj,nout=nout,note='Species indices for components to CPA.%N%3f'//
     .     'Number of of components fixed by number of indices entered')
        ncpa(j) = nout

C   ... CPA concentrations
        nm='SPEC_ATOM_C'; call gtv(trim(nm),tksw2(prgn,nm),xcpa(:,j),
     .    cindx=jj,nout=nout,note='Concentrations c for CPA components.%N%3f'//
     .    'Enter same number of elements as in SPEC_ATOM_CPA;%N%3f'//
     .    'also concentrations should sum to 1.')
        if (io_help == 0) then
          if (nout /= ncpa(j)) then
            call rxs('unequal number of elements, SPEC_ATOM_CPA and'//
     .        ' SPEC_ATOM_C for atom ',slabl(j))
          endif
          if (nout > 0) then
          if (abs(sum(xcpa(1:nout,j))-1) > 1d-8) then
            call rxs('species concentrations SPEC_ATOM_C do not sum '//
     .        ' to 1 for atom ',slabl(j))
          endif
          endif
        endif

C   ... tbe specific
        nm='SPEC_ATOM_QPOL'; call gtv(trim(nm),tksw2(prgn,nm),
     .    qpol(:,j),def_r8v=zerov,cindx=jj,note=
     .    'up to ten polarisability parameters')
        nm='SPEC_ATOM_I'; call gtv(trim(nm),tksw2(prgn,nm),stni(j),
     .    def_r8=0d0,cindx=jj,note='TB Stoner I (l=2)')
        nm='SPEC_ATOM_FRZQ1'; call gtv(trim(nm),tksw2(prgn,nm),
     .    iq1(1:nlbj*nsp,j),nmin=nlbji*nsp,def_i4v=izerv,
     .    note='(--fit mode) T: freeze site energies')
        nm='SPEC_ATOM_FRZVSO'; call gtv(trim(nm),tksw2(prgn,nm),
     .    iq2(2:nlbj,j),nmin=max(nlbji-1,NULLI),def_i4v=izerv,
     .    note='(--fit mode) T: freeze spin-orbit parameters')
        if (lncol4 .and. lmxaj >= 1 .or. io_help /= 0) then
          tbvso(:,j)=0
          nm='SPEC_ATOM_VSO'
          call gtv(trim(nm),tksw2(prgn,nm),tbvso(2:nlbj,j),
     .      def_r8v=zerov,cindx=jj,note='TB spin-orbit parameters')
        endif

        nm='SPEC_ATOM_LMXB'; call gtv(trim(nm),tksw2(prgn,nm),lmxpb(j),
     .    cindx=jj,note='compound l-cutoffs for basis (molecules)')
        nm='SPEC_ATOM_LXI'; call gtv(trim(nm),tksw2(prgn,nm),lxi(j),
     .    cindx=jj,note='l-cutoffs for interstitial product basis')
        nm='SPEC_ATOM_EXI'; call gtv(trim(nm),tksw2(prgn,nm),exi(:,j),
     .    cindx=jj,note='Hankel energies for interstitial product basis'
     .    ,nout=nout)
        if (nout >= 0) nxi(j) = nout
        nm='SPEC_ATOM_RINT'; call gtv(trim(nm),tksw2(prgn,nm),rint(j),
     .    cindx=jj,note='range of interstitial product basis')
        nm='SPEC_ATOM_RCUT'; call gtv(trim(nm),tksw2(prgn,nm),rcut(j),
     .    cindx=jj,note='radius over which true,smooth Hankels differ')

C       Some sanity checks
        if (io_help == 0 .and. lmxaj >= 0) then
          if (tksw2(prgn,'SPEC_ATOM_LFOCA') /= 2) call sanrg(io_help==0,lfoca(j),0,2,'rdctrl','lfoca')
          if (tksw2(prgn,'SPEC_ATOM_LMXL') /= 2) then
            if (.not. ltbe) then
              call sanrg(io_help==0,lmxl(j),min(0,lmxaj),max(0,lmxaj),'rdctrl','lmxl')
            else
              call sanrg(io_help==0,lmxl(j),0,2*(nl-1),'rdctrl','lmxl')
            endif
          endif
          if (tksw2(prgn,'SPEC_ATOM_KMXA') < 2) ltmp=isanrg(kmxt(j),2,25,' rdctrl (warning):','kmxa',F)
        endif

        coreh(j) = ' '
        nm='SPEC_ATOM_C-HOLE'; call gtv(trim(nm),tksw2(prgn,nm),coreh(j),nmin=10,
     .    cindx=jj,note='Channel for core hole')
        nm='SPEC_ATOM_C-HQ'; call gtv(trim(nm),tksw2(prgn,nm),coreq(:,j),
     .    def_r8v=(/-1d0,0d0/),cindx=jj,nmin=2,note=
     .    'Charge in core hole.  '//
     .    'Optional 2nd entry is moment of core hole:'//
     .    '%N%5fQ(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2'//
     .    '%N%5fQ(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2')

        nm='SPEC_ATOM_EREF'; call gtv(trim(nm),tksw2(prgn,nm),eref(j),def_r8=0d0,
     .    cindx=jj,note='Reference energy subtracted from total energy')

        nm='SPEC_ATOM_AMASS'; call gtv(trim(nm),tksw2(prgn,nm),mass(j),
     .    cindx=jj,note='Nuclear mass in a.u. (for dynamics)')

        nm='SPEC_ATOM_COLOUR'; call gtv(trim(nm),tksw2(prgn,nm),
     .    colxbs(:,j),def_r8v=zerov,cindx=jj,note=
     .    'Colour for xbs') !,-3,4,i1,TF(37))
        xx = rmt(j); if (io_help >= 1) xx = NULLI
        nm='SPEC_ATOM_RADIUS'; call gtv(trim(nm),tksw2(prgn,nm),
     .    radxbs(j),def_r8=xx,cindx=jj,note='Radius for xbs')

        if ((tksw2(prgn,'SPEC_ATOM_HCR') < 2 .or. tksw2(prgn,'SPEC_ATOM_FRZWF') < 2) .and. io_help > 0) then
          call info0(2,1,0,' ... The next tokens affect the shape of basis set:')
        endif

        nm='SPEC_ATOM_FRZWF'; call gtv(trim(nm),tksw2(prgn,nm),mxcst4(j),cindx=jj,
     .    def_lg=F,note='Set to freeze augmentation wave functions for this species')

        if (nlaj /= 0) then
C         Hard core radius defining NMTO, new FP
          nm='SPEC_ATOM_HCR'; call gtv(trim(nm),tksw2(prgn,nm),
     .      hcr(1:nlbj,j),nmin=-NULLI,cindx=jj,nout=nout,or=T,
     .      note='Hard sphere radii for structure constants')
          if (nout == 0) then   !nout=-1 if sw=2; otherwise nout=0 unless data was read
            xv(1) = .7d0        ! Intel compiler bug: must pass array to gtv
            nm='SPEC_ATOM_HCR/R'; call gtv(trim(nm),tksw2(prgn,nm),
     .        hcr(1:nlbj,j),nmin=-NULLI,cindx=jj,nout=nout,
     .        def_r8v=xv(1:1),note='Same as HCR, but in units of R')
            call dscal(nout,rmt(j),hcr(1,j),1)
          endif
          nm='SPEC_ATOM_RSMINL'; call gtv(trim(nm),tksw2(prgn,nm),
     .      rsminl(1:nlaj,j),nmin=-NULLI,cindx=jj,nout=nout,
     .      note='Lower bound to smoothing radius for this species (optionally l-dependent)')
          if (nout>0 .and. nout<nlbj .and. io_help==0) rsminl(nout+1:nlbj,j) = rsminl(nout,j)
        elseif (rmt(j) == 0) then
          hcr(1,j) = 0
        endif

        nm='SPEC_ATOM_HSFITK'; call gtv(trim(nm),tksw2(prgn,nm),hsfitp(:,j),def_r8v=hsfitd,
     .    cindx=jj,note='HSFITK*HCR = radius beyond which K.E. should approach asymptotic value'//
     .    '%N   2nd argument weights how much K.E. at RMT is added to minimization')

C        nm='SPEC_ATOM_HFAC'; call gtv(trim(nm),tksw2(prgn,nm),
C     .    alabl,nmin=0,Texist=ltmp,cindx=jj,note=
C     .    'Species-specific parameters overriding global defaults '//
C     .    'in SPEC_ATOM_HFAC.  Note that these data are used only '//
C     .    'to supply values for file shfac')
        hfacs(:,j) = 1
        nm='SPEC_ATOM_HFACB'; call gtv(trim(nm),tksw2(prgn,nm),
     .    hfacs(1,j),cindx=jj,def_r8=hfacd(1),
     .    note='Species-specific parameter overriding global SPEC_HFAC_B (which see)')
C        nm='SPEC_ATOM_HFACL'; call gtv(trim(nm),tksw2(prgn,nm),
C     .    hfacs(2,j),cindx=jj,def_r8=hfacd(2),
C     .    note='Species-specific parameter overriding global SPEC_HFAC_L (which see)')
        nm='SPEC_ATOM_HFACD'; call gtv(trim(nm),tksw2(prgn,nm),
     .    hfacs(3,j),cindx=jj,def_r8=hfacd(3),
     .    note='Scale channel bandwidth by D; see SPEC_HFAC')

      enddo                     ! Loop over species

C ... Cleanup after looping over species data
      if (io_help==0) then
C       Maximum L-cutoff
        call imxmn(nspec,lmxa,1,i0,lmxax)
        nlmax = (lmxax+1)**2
C       Set global variables nkaph and mxorb
        nkaph = nkaph+lpzi-1
C       print *, nkaph,nlmax,nkaph*nlmax
        xx = dglob('nkaph',dble(nkaph),1)
        xx = dglob('nlmax',dble(nlmax),1)
        xx = dglob('mxorb',dble(nkaph)*nlmax,1)
        nphimx = 2; if (lpzi > 1) nphimx = 3
        xx = dglob('nphimx',dble(nphimx),1)
      endif

C     Case site data read before species: reorder ips
      ltmp = .false.
      if (.not. associated(slabll,slabl)) then
        if (lread_site) then  ! Sites were loaded first: reorder ips
          allocate(ipsl(nsite)); call icopy(nsite,ips,1,ipsl,1)
          do  i = 1, nsite
            alabl = slabll(ipsl(i))
            call tokmat(alabl,slabl,nspec,8,' ',i0,nn,.true.)
            ips(i) = i0+1
            if (slabl(ips(i)) /= slabll(ipsl(i))) call rx('m_rdctrl: bug in species table')
          enddo
          deallocate(ipsl)
          ltmp = .true.
        endif
        deallocate(slabll)
        slabll => slabl
      endif

   79 continue

C --- Site ---
      nm='SITE'; sw = tkswp(prgn,express,nm); if (sw >= 2) goto 89

C     gtv call is only helpful for printout
C     call gtv(trim(nm),sw,Texist=ltmp,note='Site-specific parameters')

      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for site data ---')

      if (io_help >= 1 .and. express /= -1) then
        nbas  = 1
        nsite = 2
      endif
      if (nbas < 0 .or. lread_site) goto 89 ! STRUC data not yet given or site data already read

C     Initialize site data
      if (.not. linitsite) then
        allocate(pos(3,nsite),vel(3,nsite),eula(3,nsite),vshft(nsite),ips(nsite),
     .    ipl(nsite),plv(nsite),irlx(3,nsite),rmaxs(nsite),mpole(nsite),dpole(3,nsite))
        if (ltbe) allocate(delta(n0,nsite),ndelta(nsite))

C       Default values
        vel  = 0d0; eula = 0d0; vshft = 0d0; ipl  = 0; plv  = 0; irlx = 0
C       Should always be set, unless never read.  But see ips below
        ips  = NULLI; pos  = NULLR; mpole = 0d0; dpole = 0d0; rmaxs = NULLR

        linitsite = .true.
      endif

      nm='SITE_FILE'; call gtv(trim(nm),tksw2(prgn,nm),outs,nmin=10,
     .  nout=nout, note=
     .  '$$e(also alias STRUC_FILE; see above)%N%3f$$e'//
     .  'Name of site file containing basis and related information.'//
     .  '%N%3fIf token is given, read POS from site file,'//
     .  '%N%3fand VEC, EULA, VSHFT, PL, RLX, if they are present.')

      if (nout == 1) then
C       Read pos, and possibly vel, eula,vshft,ips,ipl,irlx
        if (procid == master) then
C         Get FILE lio to see what it contains
          j = iosite(16000,3d0,0,trim(outs),i,slabll,alat,plat,nbas,
     .               nspec,xx,xx,xx,xx,xx,xx,xx)
          ii = 8000 + 140
C         File contains vel,eula,PL,rlx,vshft
          if (mod(j/32,2) == 1) then
            ii = ii + 32000
C         Even if not, let iosite assign default values
          else
            ii = ii + 32000
          endif
! This specialisation assumes the TB site files interlace the onsite delta.
! Though the space for delta is too small (6) instead of the usual 16 tradi-
! tionally written in file delta.ext. In addition to making handling more
! difficult, the interlacing is also quite an eyesore. DMT
!           if (ltbe) ii = ii + 64000
C
C         Note : ips may need to be modified since SPEC may order species differently
C         See cleanup at the end of SPECIES category
          j = iosite(ii,3d0,0,trim(outs),i,slabll,alat,plat,nbas,
     .      nspec,pos,vel,eula,vshft,ips,ipl,irlx)
        endif
        call mpibc1(pos,3*nsite,4,mlog,'readctrl','pos')
        call mpibc1(vel,3*nsite,4,mlog,'readctrl','vel')
        call mpibc1(eula,3*nsite,4,mlog,'readctrl','eula')
        call mpibc1(vshft,nsite,4,mlog,'readctrl','vshft')
        call mpibc1(ipl,nsite,2,mlog,'readctrl','ipl')
        call mpibc1(irlx,3*nsite,2,mlog,'readctrl','irlx')
C       Must be be updated if species are read later; see associated(slabll,slabl) in SPEC
        call mpibc1(ips,nsite,2,mlog,'readctrl','ips')
        lread_site = .true.

        goto 89
      endif
      if (express == -1) goto 89  ! Assume no tags below will go participate in express read
      if (io_help >= 1 .and. express >= 0) then
        if (iprint() > 0) write(*,383)
  383   format(/' - The following tokens are input for each site. Data sandwiched'/
     .    3x,'between successive occurences of token ATOM apply to one site.'/
     .    3x,'Tags are not parsed if site data is read from a SITE file.')
        nbas = 1
        nsite = 2
      endif

C ... Site data, one pass for each atom
      do  j = 1, nbas

        if (io_help /= 0) then
          write(stdo,'(1x)')
        elseif (io_help == 0 .and. io_show>0) then
          call info(1,0,0,' ... Site %i',j,0)
        endif

        jj=(/1,j/)
        nm='SITE_ATOM'; call gtv(trim(nm),tksw2(prgn,nm),alabl,nmin=10,
     .    cindx=jj,note='Species label')
        if (io_help == 0) then
          do  i = 1, nspec
            if (trim(alabl) == trim(slabl(i)) ) then
              ips(j) = i
              goto 881
            endif
          enddo
          call rxs('Category SITE referred to'//
     .      ' nonexistent species: ',alabl)
        endif
  881   continue

C       call snit
C
C   ... Site positions
        sw = tksw2(prgn,'SITE_ATOM_XPOS')
        nm='SITE_ATOM_POS'; call gtv(trim(nm),tksw2(prgn,nm),pos(:,j),nout=nout,
     .    cindx=jj,note='Atom coordinates, in units of alat',or=(sw < 2))
        if (nout == 0 .or. tksw2(prgn,'SITE_ATOM_POS') >= 2) then !nout=-1 if sw=2; otherwise nout=0 unless data was read
          nm='SITE_ATOM_XPOS'; call gtv(trim(nm),tksw2(prgn,nm),pos(:,j),
     .      cindx=jj,note='Atom coordinates, as (fractional) '//
     .      'multiples of the lattice vectors')
          call dcopy(3,pos(1,j),1,xv,1)
          call dmpy(plat,3,1,xv,3,1,pos(1,j),3,1,3,1,3)
        endif
        nm='SITE_ATOM_DPOS'; call gtv(trim(nm),tksw2(prgn,nm),
     .    xv(1:3),def_r8v=zerov,nout=nout,cindx=jj,note=
     .    'Shift in atom coordinates added to pos')
        if (nout == 3) call daxpy(3,1d0,xv,1,pos(1,j),1)

        nm='SITE_ATOM_V0'; call gtv(trim(nm),tksw2(prgn,nm),vel(:,j),
     .    def_r8v=zerov,cindx=jj,note='Initial velocity for molecular dynamics')

        nm='SITE_ATOM_RELAX'; call gtv(trim(nm),tksw2(prgn,nm),irlx(:,j),
     .    def_i4v=(/(1,i=1,n0)/),cindx=jj,
     .    note='relax site positions (lattice dynamics) or Euler angles (spin dynamics)')

        nm='SITE_ATOM_VSHFT'; call gtv(trim(nm),tksw2(prgn,nm),vshft(j),
     .    def_r8=0d0,cindx=jj,note='Constant potential shift for this site')

        nm='SITE_ATOM_RMAXS'; call gtv(trim(nm),tksw2(prgn,nm),rmaxs(j),
     .    cindx=jj,note='Site-dependent radial cutoff for strux, in a.u.')

        nm='SITE_ATOM_ROT'; call gtv(trim(nm),tksw2(prgn,nm),outs,nmin=
     .    10,cindx=jj,nout=nout,note='Rotation of spin quantization axis at this site')
        if (io_help == 0 .and. nout == 1) then
          call numsyv(i0)
          call lodsyv('ib',0,dble(j),ii)
          call lodsyv('x', 0,pos(1,j),ii)
          call lodsyv('y', 0,pos(2,j),ii)
          call lodsyv('z', 0,pos(3,j),ii)
          call a2rotm(outs,F,iprint()-10,xv)
          call clrsyv(i0)
          call rm2eua(xv,eula(1,j),eula(2,j),eula(3,j))
          if (io_show>0) call info5(0,0,0,
     .      '%15pROT= : Euler alpha = %1;6d  beta = %1;6d'//
     .      '  gamma = %1;6d',eula(1,j),eula(2,j),eula(3,j),0,0)
        endif

        nm='SITE_ATOM_PL'; call gtv(trim(nm),tksw2(prgn,nm),ipl(j),
     .    def_i4=0,cindx=jj,note='Assign principal layer number to this site')

        nm='SITE_ATOM_PLV'; call gtv(trim(nm),tksw2(prgn,nm),plv(j),
     .    def_i4=0,cindx=jj,note='Assign PL potential index to this site')

        nm='SITE_ATOM_DELTA'; sw = tksw2(prgn,nm)
        if (sw < 2 .and. ltbe) then
          if (io_help >= 1) then
            nlbji = NULLI
            nlbj = 1
          elseif (lmxb(i) == NULLI) then
            nlbji = NULLI
            nlbj = 0
          else
            i = ips(j)
            nlbj = 1+lmxb(i)
            nlbji = nlbj
          endif
          ndelta(j)=nlbj*nsp
          call gtv(trim(nm),sw,delta(1:nlbj*nsp,j),nmin=nlbji*nsp,
     .      cindx=jj,def_r8v=zerov,note='Vector of on-site energy shifts')
        endif

C        nm='SITE_ATOM_SID'; call gtv(trim(nm),tksw2(prgn,nm),sid(j),
C     .    def_i4=0,cindx=jj,note=
C     .    'Site ID, used in the merging of coincident clusters')
C
      enddo

C ... Input for point multipoles
      sw = tksw2(prgn,'SITE_PM')
      if (sw < 2) then
      if (io_help == 1) nbasp=nbas+1
      if (nbasp > nbas) then
        if (io_show>0) write
     .    (stdo,'(/'' The following input is for point multipoles:'')')
        do j = nbas+1, nbasp

          jj=(/1,j-nbas/)
C          nm='SITE_PM_POS'; call gtv(trim(nm),tksw2(prgn,nm),pos(:,j),
C     .      nout=nout,cindx=jj,note='PM coordinates, in units of alat')

          nm='SITE_PM_POS'; call gtv(trim(nm),tksw2(prgn,nm),pos(:,j),
     .      cindx=jj,note='PM coordinates, in units of alat')

          nm='SITE_PM_Q'; call gtv(trim(nm),tksw2(prgn,nm),mpole(j),
     .      def_r8=0d0,cindx=jj,note='Point monopole moment')

          nm='SITE_PM_P'; call gtv(trim(nm),tksw2(prgn,nm),dpole(:,j),
     .      def_r8v=zerov,cindx=jj,note='Point dipole moment')
        enddo
      endif
      endif

   89 continue

C --- Structure constants ---
      if (tkswp(prgn,express,'STR') < 2) then
      if (io_show+io_help/=0 .and. tksw2(prgn,'STR') < 2) call info0(2,1,0,
     .    ' --- Parameters for structure constants ---')
      nm='STR_RMAXS'; call gtv(trim(nm),tksw2(prgn,nm),str_rmax,nout=nout,
     .  note='Radial cutoff for strux, in a.u.',or=T)
      if (nout == 0) then !nout=-1 if sw=2; otherwise nout=0 unless data was read
        nm='STR_RMAX'; call gtv(trim(nm),tksw2(prgn,nm),str_rmax,nout=nout,def_r8=0d0,
     .    note='Radial cutoff for strux, in units of avw',or=T)
        str_rmax = str_rmax*avw
      endif
      if (str_rmax == 0d0) then
        nm='STR_NEIGHB'; call gtv(trim(nm),tksw2(prgn,nm),str_mnn,Texist=ltmp,def_i4=30,
     .    note='Minimum number of neighbors in cluster')
        if (.not. ltmp) str_mnn=-30 ! sign flags nothing was read
      endif

C     Mode for structure constants
      nm='STR_ENV_MODE'; call gtv(trim(nm),tksw2(prgn,nm),str_mode,def_i4=0,
     .  note='Type of envelope functions:'//
     .  '%N%3f 0: 2nd generation E=0, screening defined by alpha'//
     .  '%N%7fUse 100 for Sdot fixed by old Kanpur conventions'//
     .  '%N%3f    For all other modes, screening defined by HCR.'//
     .  '%N%3f 1: unit-zero val Hankel-based SSWs, NEL energies'//
     .  '%N%3f 2: NMTOs, NEL energies (OKA conventions)'//
     .  '%N%3f11: unit-zero val-slo Hankel SSWs, 2 (linked) energies'//
     .  '%N%3f13: Similar to 11, but also Gaussian val-lap strux')

      nkapsi = 1
      if (str_mode == 2) nkapsi=2
      if (io_help >= 1) nkapsi = NULLI
      if (str_mode == 0) then
        str_nkaps = 1
      else
        nm='STR_ENV_NEL'; call gtv(trim(nm),tksw2(prgn,nm),str_nkaps,
     .    def_i4=nkapsi,note='Number of SSSW or NMTO energies')
        if (str_mode == 2) call sanrg(io_help==0,str_nkaps,2,5,'rdctrl','NEL')
        if (str_mode == 1) call sanrg(io_help==0,str_nkaps,1,1,'rdctrl','NEL')
        if (str_mode == 11 .or. str_mode == 13)
     .    call sanrg(io_help==0,str_nkaps,2,2,'rdctrl','NEL')
        nkapsi = str_nkaps
        if (io_help >= 1) then
          str_nkaps = 1
          nkapsi = NULLI
        endif
        nm='STR_ENV_EL';call gtv(trim(nm),tksw2(prgn,nm),
     .    str_kaps(1:str_nkaps),nmin=nkapsi,def_r8v=(/0d0,-1d0,2.3d0/),
     .    note='SSSW or NMTO energies, atomic units')
        nm='STR_ENV_RSMINL';call gtv(trim(nm),tksw2(prgn,nm),str_rsminl,nmin=-NULLI,
     .    nout=nout,note='Lower bound to smoothing radius (optionally l-dependent)')
      endif

      nm='STR_MXNBR'; call gtv(trim(nm),tksw2(prgn,nm),str_mxnbr,
     .  def_i4=0,note='Max number of nbrs (for dimensioning arrays)')
      nm='STR_SHOW'; call gtv(trim(nm),tksw2(prgn,nm),str_lshow1,
     .  def_lg=F,note='Show strux after generating them')
      nm='STR_EQUIV'; call gtv(trim(nm),tksw2(prgn,nm),str_lequiv1,
     .  def_lg=F,note='Look for equivalent strux')

      nm='STR_LMAXW'; call gtv(trim(nm),tksw2(prgn,nm),str_lmaxw,def_i4=-1,
     .  note='l-cutoff for Watson sphere, used to help localize strux')
      nm='STR_DELRW'; call gtv(trim(nm),tksw2(prgn,nm),str_drwats,
     .  def_r8=0.1d0,note='Padding beyond cluster for Watson sphere')

      nm='STR_DELRX'; call gtv(trim(nm),tksw2(prgn,nm),str_delrx,def_r8=3d0,
     .  note='Range of screened function beyond last site in cluster')

      nm='STR_TOLG'; call gtv(trim(nm),tksw2(prgn,nm),str_tolg,def_r8=1d-6,
     .  note='Tolerance in l=0 gaussians, which determines their range')

      if (io_help /= 0 .and. tksw2(prgn,'STR_IINV_NIT') < 2) call info0(2,1,0,
     .  ' * IINV parameters govern iterative solutions to screened strux')
      nm='STR_IINV_NIT'; call gtv(trim(nm),tksw2(prgn,nm),iinv_nit,def_i4=0,
     .  note='Number of iterations')
      nm='STR_IINV_NCUT'; call gtv(trim(nm),tksw2(prgn,nm),iinv_ncut,def_i4=0,
     .  note='Number of sites for inner block')
      nm='STR_IINV_TOL'; call gtv(trim(nm),tksw2(prgn,nm),iinv_tol,def_r8=0d0,
     .  note='Tolerance in errors')

      nm='STR_RVL/R'; call gtv(trim(nm),tksw2(prgn,nm),str_rmaxg,def_r8=0.7d0,
     .  note='Radial cutoff for val-lap, in units of RMAX')

      nm='STR_VLFUN'; call gtv(trim(nm),tksw2(prgn,nm),str_ivl,def_i4=0,
     .  note='Functions for val-lap basis'//
     .  '%N%10f0 G0 + G1'//
     .  '%N%10f1 G0 + Hsm'//
     .  '%N%10f2 G0 + Hsm-dot')

      endif ! STR

C --- Two-center fit (molecules code) ---
      if (tksw2(prgn,'TCF') < 2) then
      if (io_help >= 1) call info(1,1,0,
     .  ' * The following tokens are used for the 2-center fit',0,0)
      nm='TCF_NBISI'; call gtv(trim(nm),tksw2(prgn,nm),tcf_nbisi,
     .  note='Mesh parameters for TCF ')
      nm='TCF_NALF'; call gtv(trim(nm),tksw2(prgn,nm),tcf_nalf,
     .  note='Polynomial order for d-dependence of TCF fit coffs ')
      nm='TCF_NCUPL'; call gtv(trim(nm),tksw2(prgn,nm),tcf_ncupl,
     .  def_i4=8,note='Reduced polynomial order in TCF ')
      nm='TCF_NDUST'; call gtv(trim(nm),tksw2(prgn,nm),tcf_ndust,
     .  def_i4=2,note='No. iterations to improve on TCF ')
      nm='TCF_ADEC'; call gtv(trim(nm),tksw2(prgn,nm),tcf_adec,
     .  def_r8=1d0,note='Log spacing in distance for TCF ')
      nm='TCF_WZTCF'; call gtv(trim(nm),tksw2(prgn,nm),tcf_wztcf,
     .  def_r8=0d0,note='Weighting of z-points in TCF ')
      endif ! TCF

C --- Brillouin Zone ---
      if (tkswp(prgn,express,'BZ') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for Brillouin zone integration ---')

      nm='BZ_GETQP'; call gtv(trim(nm),tksw2(prgn,nm),bz_lio1,
     .  def_lg=F,note='Read qp from disk',or=T)
      nm='BZ_NKABC'; sw=tksw2(prgn,nm); if (bz_lio1) sw = 2
      call gtv(trim(nm),sw,bz_nabc,nout=nout,
     .  note='No. qp along each of 3 lattice vectors.'//
     .  '%N%3fSupply one number for all vectors or a separate '//
     .  'number for each vector.')
      if (nout > 0) then
      call fill3in(nout,bz_nabc)
      if (io_help == 0) call fill3inl(bz_nabc,plat,io_show /= 0,'NKABC')
      endif

      nm='BZ_PUTQP'; call gtv(trim(nm),tksw2(prgn,nm),bz_lio2,
     .  def_lg=F,note='Write qp to disk')

      nm='BZ_BZJOB';call gtv(trim(nm),tksw2(prgn,nm),bz_lshft,nout=nout,
     .  def_i4v=izerv(1:1),note=
     .  '0 centers BZ mesh at origin, 1 centers off origin'//
     .  '%N%3fSupply one number for all vectors or a separate '//
     .  'number for each vector.')
      call fill3in(nout,bz_lshft)

      nm='BZ_METAL'; call gtv(trim(nm),tksw2(prgn,nm),bz_lmet,def_i4=bz_lmet,
     .  note='Controls management of k-point integration weights'//
     .  '%N%3f0: assume insulator'//
     .  '%N%3f1: save condensed evecs on disk (ASA only)'//
     .  '%N%3f2: read kp weights from disk (from prior iteration)'//
     .  '%N%3f3: Always make 2 band passes'//
     .  '%N%3f4: Density by sampling, 3-point interpolation'//
     .  '%N%3f5: keep evecs in memory (lmf only)')

      nm='BZ_TETRA'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lmet2,def_i4=1,note=
     .  '1: Tetrahedron integration%N%2f'//
     . '11: Sampling integration for density, '//
     .      'but Ef and sumev by tetrahedron')
      if (nout > 0 .and. ctrl_lmet2 /= 0 .and.
     .  ctrl_lmet2 /= 1 .and. ctrl_lmet2 /= 11) then
        call rx('rdctrl: BZ_TETRA must be at 0,1, or 11')
      endif

      nm='BZ_N'; call gtv(trim(nm),tksw2(prgn,nm),bz_n,def_i4=0,note=
     .  'N>0: Polynomial order for Methfessel-Paxton sampling%N%3f'//
     .  'N=0: Conventional Gaussian sampling%N%3f'//
     .  'N<0: Broadening by Fermi-Dirac distribution%N%3f'//
     .  'To be used in conjunction with W, next')
      nm='BZ_W'; call gtv(trim(nm),tksw2(prgn,nm),bz_w,def_r8=5d-3,note=
     .  'N>=0: Line broadening for sampling integration%N%3f'//
     .  'N<0 : Temperature for Fermi distribution (Ry)')

      nm='BZ_EF0'; call gtv(trim(nm),tksw2(prgn,nm),bz_ef,
     .  def_r8=0d0,note='Initial guess at Fermi energy')
      nm='BZ_DELEF'; call gtv(trim(nm),tksw2(prgn,nm),bz_def,
     .  def_r8=0.05d0,note='Initial uncertainty in Fermi energy')
      nm='BZ_ZBAK'; call gtv(trim(nm),tksw2(prgn,nm),zbak(1),
     .  def_r8=0d0,note='Homogeneous background charge')
      if (prgnam=='LMF') then
        nm='BZ_SAVDOS'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_ldos,def_i4=0,
     .    note=' Write DOS to directly disk (NPTS and DOS also needed)')
        call sanrg(io_help==0.and.nout > 0,ctrl_ldos,0,1,'rdctrl','SAVDOS')
      else
        nm='BZ_SAVDOS'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_ldos,
     .    def_i4=0,note='Choose some combination of the following:'
     .    //'%N%3f1 Write DOS to directly disk (NPTS and DOS also needed)'
     .    //'%N%3f2 Write weights for partial DOS'
     .    //'%N%3f4 Same as (2), but weights m-resolved')
      endif

      nm='BZ_DOS'; call gtv(trim(nm),tksw2(prgn,nm),bz_dosw,
     .  def_r8v=(/-1d0,0d0/),note='Energy window over which DOS accumulated')
      nm='BZ_NPTS'; call gtv(trim(nm),tksw2(prgn,nm),bz_ndos,def_i4=1001,
     .  note='No. DOS points (sampling integration, and lmdos)')

      xx = 2d0; if (ltbe) xx = 5d0
      nm='BZ_EFMAX'; call gtv(trim(nm),tksw2(prgn,nm),bz_efmax,
     .  def_r8=xx,note='Find evecs up to efmax')
      nm='BZ_NEVMX'; call gtv(trim(nm),tksw2(prgn,nm),bz_nevmx,
     .  def_i4=0,note='Find at most nevmx eigenvectors'//
     .  '%N%3fIf NEVMX=0, program uses internal default'//
     .  '%N%3fIf NEVMX<0, no eigenvectors are generated')
      bz_zval = NULLR
      nm='BZ_ZVAL'; call gtv(trim(nm),tksw2(prgn,nm),bz_zval,
     .  note='Number of electrons to accumulate in BZ integration')
      nm='BZ_NOINV'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lqp1,def_lg=F,
     .  note='Suppress automatic inclusion of inversion symmetry for BZ')
      nm='BZ_FSMOM'; call gtv(trim(nm),tksw2(prgn,nm),bz_fsmom,
     .  def_r8v=(/NULLR,0d0/),note='Fixed-spin moment (fixed-spin moment method)'//
     .  '%N%3fOptional 2nd arg supplies initial Beff')
      nm='BZ_DMATK'; call gtv(trim(nm),tksw2(prgn,nm),bz_lio48,def_i4=0,
     .  note='Calculate density matrix n'//
     .  '%N%3f1: Rotate n(k) to orbital basis'//
     .  '%N%3f2: Rotate (n(k)-1/2) to orbital basis')
      call sanrg(io_help==0.and.nout > 0,bz_lio48,0,2,'rdctrl','BZ_DMATK')
      nm='BZ_INVIT'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_lqp2,
     .  def_lg=T,note='Use inverse iteration for diagonalization')
      i0 = 0; if (ltbe) i0 = -1
      nm='BZ_MULL'; call gtv(trim(nm),tksw2(prgn,nm),bz_lmull,
     .  def_i4=i0,note='Mulliken population analysis ')
      nm='BZ_EMESH'; call gtv(trim(nm),tksw2(prgn,nm),bz_semsh,nmin=10,
     .  def_r8v=(/10d0,0d0,-1d0,0d0,0.01d0,0d0,0d0,0d0,1d-10,0d0/),
     .  note='Mesh for energy contour integration'
     .  //'%N%3fentry 1   2    3    4    5   6    7    8    9    10'
     .  //'%N%9fnz mode ebot etop  ----  depend on mode ---- '
     .  //'(see documentation)')

      nm='BZ_COND_MODE'; call gtv(trim(nm),tksw2(prgn,nm),i0,nout=nout,
     .  def_i4=0,note='Mode for type conductivity calculation'//
     .  '%N%3f0: Calculate DOS D(E) instead of conductivity'//
     .  '%N%3f1: ballistic conductivity <v/2.DIR> (must also supply DIR)'//
     .  '%N%3f2: <v_i . v_j> (diffusive conductivity)'//
     .  '%N%3f3: <|v|> absolute value of velocity'//
     .  '%N%6fIn this mode one or two elts of DIR should be 1, '//
     .  'marking i or i,j')
      bz_lcond(1) = i0
      nm='BZ_COND_DIR'
      if (io_help /= 0 .or. i0 > 0 .and. i0 < 3) then
        sw = tksw2(prgn,nm)
      else
        sw = 2
      endif
      call gtv(trim(nm),sw,bz_lcond(2:),
     .  note='Conductivity direction vector DIR, used for MODE>0')

      nm='BZ_ZINSUL'; call gtv(trim(nm),tksw2(prgn,nm),bz_zinsul,def_r8=0d0,
     .  note='Deviation from charge neutrality for GF fermi level finder (insulator mode)')

      nm='BZ_SYML'; sw=tksw2(prgn,nm); call gtv(trim(nm),sw,strn_syml,nmin=10,
     .  nout=nout,note='String to generate symmetry lines')

      endif                     !BZ

C --- Ewald sums ---
      if (tkswp(prgn,express,'EWALD') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for Ewald sums ---')

      nm='EWALD_AS'; call gtv(trim(nm),tksw2(prgn,nm),lat_as,def_r8=2d0,
     .  note='Ewald smoothing parameter')
      nm='EWALD_TOL'; call gtv(trim(nm),tksw2(prgn,nm),lat_tol,def_r8=1d-8,
     .  note='Ewald tolerance')
      nm='EWALD_NKDMX'; call gtv(trim(nm),tksw2(prgn,nm),lat_nkdmx,def_i4=800,
     .  note='Ewald tolerance')
      nm='EWALD_RPAD'; call gtv(trim(nm),tksw2(prgn,nm),lat_rpad,def_r8=0d0,
     .  note='Scale rcutoff by rpad when lattice vectors padded in oblong geometries')
      endif

C --- Iterations (formerly MIX) ---
      smix(2) = NULLI           ! Not set
      if (tkswp(prgn,express,'ITER') >= 2) goto 99
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for iterations ---')
C     Default values for smix (array has same same structure as lstra smix)
      smix = 0
      smix(2) = 1               ! beta
      smix(4) = 1               ! bv
C     smix(5) =  0              ! elind
      call s8tor8('mixm',smix(6)) ! file name
      smix(7) =  0              ! nkill
      smix(8) =  0              ! lxpot
      smix(9) = -1              ! mmix
      smix(10) = 0              ! mode (0=Anderson)
C     smix(11) = 0              ! model = previous mixing mode
C     smix(12) = 0              ! n     = Number of iterations for this species
C     smix(13) = 1              ! nitu  = max number of LDA+U iterations
C     smix(14) = 1              ! nmix  = actual number mixed
      smix(15) = 8              ! nsave = # iter to save on disk
      call s8tor8(' ',smix(16)) ! expression for rmscst
      smix(32) = 0              ! tolu
      smix(33) = 1              ! umix (mixing parm for LDA+U)
      smix(34) = 1              ! w(1)
      smix(35) = 1              ! w(2)
      smix(37) = -1             ! wc
      smalit = NULLI
      call dcopy(37,smix,1,samix,1)
      call s8tor8('mixa',samix(6)) ! file name

      nm='ITER'; sw = tkswp(prgn,express,nm); if (sw >= 2) goto 99

      nm='ITER_NIT'; call gtv(trim(nm),tksw2(prgn,nm),iter_maxit,def_i4=1,
     .  note='maximum number of iterations in self-consistency cycle')

      nm='ITER_NRMIX'; call gtv(trim(nm),tksw2(prgn,nm),smalit,def_i4v=(/80,2/),note=
     .  'Sphere program, (1) max iter; (2) no. prior iter for Anderson mixing ')

      nm='ITER_MIX'; sw=tksw2(prgn,nm); call gtv(trim(nm),sw,iter_mix,nmin=10,
     .  nout=nout,note='Mixing rules for charge mixing.  Syntax:')

      if (io_help == 1 .and. tkswp(prgn,express,nm)<2) print 345
  345 format(3x,'A[nmix][,b=beta][,b2=b2][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill][,elind=#][;...]  or'/
     .3x,'B[nmix][,b=beta][,b2=b2][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn][,elind=#][,k=nkill]'/
     .3x,'Some codes do not use all parameters; see documentation.')

C ... Sanity check and printout
      if (nout==1 .and. io_help==0 .and. sw< 2) then
        call r8tos8(smix(6),alabl)
        if (io_show==0) call pshpr(1)
        xx = 0
        if (.not. parmxp(-1,iter_mix,len_trim(iter_mix),
     .    int(smix(10)),int(smix(14)),smix(34),smix(2),xx,smix(5),xx,
     .    xx,alabl,smix(37),smix(7),smix(4),xx)) then
          if (io_show==0) call poppr
          call rx( 'm_rdctrl: parse in parmxp failed')
        endif
        if (io_show==0) call poppr
      endif

      nm='ITER_AMIX'; sw=tksw2(prgn,nm); call gtv(trim(nm),sw,iter_amix,nmin=10,
     .  nout=nout,note='Mixing rules for Euler angles.  Syntax as in MIX')
C ... Sanity check and printout
      if (nout==1 .and. io_help==0 .and. sw< 2) then
        call r8tos8(samix(6),alabl)
        if (io_show==0) call pshpr(1)
        if (.not. parmxp(-1,iter_amix,len_trim(iter_amix),
     .    int(samix(10)),int(samix(14)),samix(34),samix(2),xx,samix(5),
     .    xx,xx,alabl,samix(37),samix(7),samix(4),xx)) then
          if (io_show==0) call poppr
          call rx('m_rdctrl: parse in parmxp failed')
        endif
        if (io_show==0) call poppr
      endif

      if (prgnam /= 'LMMC') then
        nm='ITER_CONV'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(3),def_r8=1d-5,
     .    note='Tolerance in energy change from prior iteration for self-consistency')
        nm='ITER_CONVC'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(1),def_r8=3d-5,
     .    note='Tolerance in output-input charge for self-consistency')
        nm='ITER_DHFTOL'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(4),def_r8=0d0,
     .    note='Tolerance in correction to Harris forces.  0=>not used')
      else
C   ... CONV and CONVC are packed in the opposite order in LMMC
        nm='ITER_CONV'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(1),def_r8=1d-4,
     .    note='Tolerance in energy change from prior iteration for self-consistency')
        nm='ITER_CONVC'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(3),def_r8=1d-4,
     .    note='Tolerance in output-input charge for self-consistency')
        ctrl_tol(2) = ctrl_tol(3)
c       nm='ITER_QTOLSP'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(2),
c    .    def_r8=1d-4,note='Tolerance in sphere charge'//
c    .    ' for self-consistency')
c       nm='ITER_QTOLI'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_tol(3),
c    .    def_r8=1d-4,note='Tolerance in interstitial charge'//
c    .    ' for self-consistency')
      endif

      i0 = 0
      nm='ITER_XIPMX'; call gtv(trim(nm),tksw2(prgn,nm),i0,def_i4=0,
     .  nout=nout,note='Mix potential independently of charge:'//
     .  '%N%3fXIPMX=1: mix vin and v(qmix)'//
     .  '%N%3fXIPMX=2: mix vin and v(qout)')
      if (nout > 0) then
        call sanrg(io_help==0,i0,0,3,'rdctrl','XIPMX')
        if (i0 > 0) lves = .true.
      endif
      smix(8) = i0

      nm='ITER_UMIX'; call gtv(trim(nm),tksw2(prgn,nm),smix(33),
     .  def_r8=1d0,note='Mixing parameter for densmat in LDA+U')
      nm='ITER_TOLU'; call gtv(trim(nm),tksw2(prgn,nm),smix(32),
     .  def_r8=0d0,note='Tolerance for densmat in LDA+U')
      nm='ITER_NITU'; call gtv(trim(nm),tksw2(prgn,nm),i0,
     .  def_i4=0,note='Max number of LDA+U iterations of densmat')
      smix(13) = i0

C     sw = tkswp(prgn,express,'ITER_FIT') ! not suitable for 'express'
      sw = tkswp(prgn,express,'ITER_FIT')
      if (sw < 2) then
      if (io_help >= 1) call info0(2,1,0,
     .  ' - The FIT tokens are used in fitting hamiltonian parameters to bands')
      if (prgnam == 'LM') then
        nm='ITER_FIT_MODE'; call gtv(trim(nm),tksw2(prgn,nm),fit_mode,
     .    note='0: No fitting'//
     .    '%N%3f1s digit governs parameters to fit:'//
     .    '%N%5f1: fit ASA C and Delta'//
     .    '%N%3f10s digit governs data to fit:'//
     .    '%N%5f0: fit to bands'//
     .    '%N%5f1: fit to site partial DOS (see documentation)')
      else
        nm='ITER_FIT_MODE'; call gtv(trim(nm),tksw2(prgn,nm),fit_mode,
     .    note='0: fit parms (use --fit)'//
     .    '%N%3f1: 1 fit in range of RFIT (not implemented)')
      endif
      if (fit_mode>=0 .or. io_help >= 1) then
      nm='ITER_FIT_NBFIT'; call gtv(trim(nm),tksw2(prgn,nm),fit_nbfit,
     .    def_i4v=(/1,0/),nmin=2,note=
     .    'Indices to lowest and highest bands to include in fit'//
     .    ' (10s dig mode=0)'//
     .    '%N%3fMissing second argument => highest band in basis')
      call sanrg(io_help==0.and.fit_nbfit(1) <= 0,fit_nbfit,1,1,'rdctrl','1st arg ITER_FIT_NBFIT')
      i = fit_nbfit(1); if (io_help >= 1) i = NULLI
      nm='ITER_FIT_NBFITF'; call gtv(trim(nm),tksw2(prgn,nm),fit_nbfitf,
     .  def_i4=i,note='Index to first reference band to be fit'//
     .    ' (10s dig mode=0)')
      if (tksw2(prgn,nm) /= 2) then
        if (fit_nbfitf == NULLI) fit_nbfitf = fit_nbfit(1)
        call sanrg(io_help==0.and.fit_nbfitf <= 0,fit_nbfitf,1,1,'rdctrl','ITER_FIT_NBFITF')
      endif
      nm='ITER_FIT_NDOS'; call gtv(trim(nm),tksw2(prgn,nm),fit_ndos,
     .  def_i4=NULLI,note='Number of energy points to fit DOS'//
     .    ' (10s dig mode=1)')
      nm='ITER_FIT_WG'; call gtv(trim(nm),tksw2(prgn,nm),fit_gausw,
     .  nout=nout,note='Broaden calculated DOS with Gaussian, width W'//
     .  '%N%3fOptional 2nd arg specifies different width for '//
     .  'reference DOS')
      if (nout == 0) call dpzero(fit_gausw,2)
      if (nout == 1) fit_gausw(2) = fit_gausw(1)
      nm='ITER_FIT_EBFIT'; call gtv(trim(nm),tksw2(prgn,nm),fit_ebfit,
     .  nmin=2,note='Energy range to fit eigenvalues or DOS')
      if (io_help >= 1 .or. fit_mode == 1) then
      nm='ITER_FIT_RFIT'; call gtv(trim(nm),tksw2(prgn,nm),fit_rmfit,
     .  nmin=2,note=
     .    '(MODE=1) range for fitting Slater-Koster parameters')
      endif
      nm='ITER_FIT_WT'; call gtv(trim(nm),tksw2(prgn,nm),fit_wt,
     .  def_r8v=(/0d0,0d0,.1d0/),note=
     .  'Fitting weights wt_kn of fit bands E_kn or DOS to file data'//
     .  '%N%3fArg #1 defines the functional form of the weights.'//
     .  '%N%3fArg 2 is an energy shift, labeled a below'//
     .  '%N%3fArg 3 is an energy scale, labeled b below'//
     .  '%N%3fform:'//
     .  '%N%3f0: unit weights'//
     .  '%N%3f1: Fermi function, wt = {1+exp[-(E_kn-(a+Ef))/b]}^-1'//
     .  '%N%3f2: Fermi-like gaussian cutoff, '//
     .  'wt = {1+exp[-((E_kn-(a+Ef))/b)^2]}^-1'//
     .  '%N%3f3: Exponential weights, wt = exp[-|(E_kn-(a+Ef))/b|]'//
     .  '%N%3f4: Gaussian weights, wt = exp[-((E_kn-(a+Ef))/b)^2]'//
     .  ' ')
      if (fit_wt(1)-nint(fit_wt(1)) /= 0) call rx1(
     .  'Illegal ITER_FIT_WT(1)=%d (must be integer)',fit_wt)
      nm='ITER_FIT_SHFT'; call gtv(trim(nm),tksw2(prgn,nm),fit_shft(1),
     .  def_i4=0,note='To set constant potential shift to align fit '//
     .  'with reference bands or DOS'//
     .  '%N%3f0: no shift is added'//
     .  '%N%3f1: Align Ef with reference bands or DOS'//
     .  '%N%3f2: Not implemented')

      nm='ITER_FIT_LAM'; call gtv(trim(nm),tksw2(prgn,nm),fit_alam,
     .  def_r8=1d0,note='Initial value of Levenberg-Marquardt lambda '//
     .  '(cf Numerical Recipes)')
      nm='ITER_FIT_SCL'; call gtv(trim(nm),tksw2(prgn,nm),fit_alsc,
     .  def_r8=5d0,note='Scale factor for Levenberg-Marquardt lambda '//
     .  '(cf Numerical Recipes)')

      endif
      endif
   99 continue

C --- TB ---
      if (tksw2(prgn,'TB') /= 2) then
      if (io_show+io_help/=0)
     .  call info0(2,1,0,' --- Tight binding parameters ---')
      nm='TB_OVLP'; call gtv(trim(nm),tksw2(prgn,nm),ltb1,
     .  def_lg=F,note='Non orthogonal tight-binding ')
      nm='TB_CRYSF'; call gtv(trim(nm),tksw2(prgn,nm),ltb2,
     .  def_lg=F,note='Crystal field terms in hamiltonian ')
      nm='TB_OVCF'; call gtv(trim(nm),tksw2(prgn,nm),ltb4,
     .  def_lg=F,note='Crystal field terms in overlap ')
      nm='TB_ADDES'; call gtv(trim(nm),tksw2(prgn,nm),ltb8,
     .  def_lg=F,note='Add ebarLL'' * sLL'' to hLL''')
      nm='TB_RMAXH'; call gtv(trim(nm),tksw2(prgn,nm),str_rmax,
     .  def_r8=0d0,note='Hamiltonian cut-off length in units of a')
      nm='TB_FORCES'; call gtv(trim(nm),tksw2(prgn,nm),ltb16,
     .  def_lg=T,note='Calculate forces ')
      nm='TB_FIJ'; call gtv(trim(nm),tksw2(prgn,nm),ltb32,
     .  def_lg=F,note='To get forces by atom pair, e.g., for stresses ')
      nm='TB_3PV'; call gtv(trim(nm),tksw2(prgn,nm),ltb128,
     .  def_lg=F,note='Calculate pressure')
      nm='TB_EVDISC'; call gtv(trim(nm),tksw2(prgn,nm),ltb256,
     .  def_lg=T,note='Can be F for insulators or to save space for metals ')
      nm='TB_PAIR'; call gtv(trim(nm),tksw2(prgn,nm),ltb512,
     .  def_lg=F,note='Pair potential only')
      nm='TB_TRH'; call gtv(trim(nm),tksw2(prgn,nm),ltb1024,
     .  def_lg=F,note='Calculate local projection of band energies')
      nm='TB_RHO'; call gtv(trim(nm),tksw2(prgn,nm),ltb2048,
     .  def_lg=F,note='Calculate local projection of charges')
      nm='TB_SCALE'; call gtv(trim(nm),tksw2(prgn,nm),ltb4096,
     .  def_lg=F,note='Scale cut offs and rmaxh with alat')
      nm='TB_UL'; call gtv(trim(nm),tksw2(prgn,nm),ltb215,
     .  def_lg=F,note='For electrostatics with L>=0')
      nm='TB_TBU'; call gtv(trim(nm),tksw2(prgn,nm),ltb213,
     .  def_lg=F,note='Tight-binding+U')
      nm='TB_GAMMA'; call gtv(trim(nm),tksw2(prgn,nm),ltb217,
     .  def_lg=F,note='Do gamma-point only')
      nm='TB_MOL'; call gtv(trim(nm),tksw2(prgn,nm),ltb218,
     .  def_lg=F,note='molecule: no PBCs')
      if (io_help >= 1) call info0(2,1,0,
     .  ' * if U1=T, UL=T or TBU=T, the next two tokens are used:')
      nm='TB_NOUAVG'; call gtv(trim(nm),tksw2(prgn,nm),ltb214,
     .  def_lg=F,note='Use individual U_l from Q=; don''t average them')
      nm='TB_IODEL'; call gtv(trim(nm),tksw2(prgn,nm),ltb216,
     .  def_lg=F,note='Attempt to read increments from disk')

      endif                     !TB parameters

C --- BANDS (do not maintain for now) ---
C      if (tksw2(prgn,'BANDS') >= 2) goto 259
C
C      nm='BANDS_MODE'; call gtv(trim(nm),tksw2(prgn,nm),j,def_i4=0,note=
C     .        '0: do nothing'//
C     .   '%N%3f1: write file for symmetry line mode'//
C     .   '%N%3f2: write file with list of qp'//
C     .   '%N%3f3: write file for 2D mesh mode')
C      if (j == 0 .and. io_help == 0) goto 259
C
C      nm='BANDS_FN'; call gtv(trim(nm),tksw2(prgn,nm),a,
C     .  note='k-points input file name',nmin=10,nout=nout)
C      if (nout == 0) a = 'qp'
C
C      if (j == 1 .or. io_help >= 1) then
C        nm='BANDS_SYML'; call gtv(trim(nm),tksw2(prgn,nm),bigstr,
C     .    note='Contents of qpts symmetry line file',nout=nout)
C      endif
C      if (j == 2 .or. io_help >= 1) then
C        nm='BANDS_QP'; call gtv(trim(nm),tksw2(prgn,nm),bigstr,
C     .    note='Contents of qpts file',nout=nout)
C      endif
C      if (j == 3 .or. io_help >= 1) then
C        nm='BANDS_MESH'; call gtv(trim(nm),tksw2(prgn,nm),bigstr,
C     .    note='Contents of qpts mesh file',nout=nout)
C      endif
C
C      if (nout == 1) then
C        j = fopna(a,-1,0)
C        i0 = 1
C        call spchar(1,ch(1:1))
C        call spchar(2,ch(2:2))
C        ch(3:3) = ' '
C        do while (.true.)
C          i = index(bigstr(i0:),ch(1:1))
C          if (i > 0) then
C            write(j,"(a)") trim(bigstr(i0:i0+i-2))
C            i0 = i0+i+verify(bigstr(i0+i:),ch(2:3))-1
C          else
C            write(j,"(a)") trim(bigstr(i0:))
C            exit
C          endif
C        enddo
C        call info0(10,0,0,' m_rctrl: wrote syml file to file '//
C     .    trim(a))
C      endif
C  259 continue                  ! end of bands

C --- Optics ---
      if (tkswp(prgn,express,'OPTICS') >= 2) goto 279
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for optics ---')

      nm='OPTICS_MODE'; call gtv(trim(nm),tksw2(prgn,nm),ctrl_loptc,def_i4=0,note=
     .        '0: do nothing;'//
     .   '%N%3f1: generate linear eps_2'//
     .   '%N%3f2: generate linear eps_2, spin 2'//
     .   '%N%3f   Add 20 for Second Harmonic Generation'//
     .   '%N%3f8: eps_2 for nonequilibrium absorption.  KT, IMREF are required inputs.'//
     .   '%N%3f9: eps_2 for photoluminescence.  KT, IMREF are required inputs.'//
     .   '%N%2f-1: generate joint DOS'//
     .   '%N%2f-2: generate joint DOS, spin 2 (LTET=3)'//
     .   '%N%2f-3: generate up-down joint DOS (LTET=3)'//
     .   '%N%2f-4: generate down-up joint DOS (LTET=3)'//
     .   '%N%2f-5: generate spin-up single DOS (LTET=3)'//
     .   '%N%2f-6: generate spin-dn single DOS (LTET=3)'//
     .   '%N%2f-8: Joint DOS for nonequilibrium absorption.  KT, IMREF are required inputs.'//
     .   '%N%2f-9: Joint DOS for photoluminescence.  KT, IMREF are required inputs.')

      if (io_show+io_help/=0) call info0(2,1,0,' ...  The following are read only if MODE is nonzero:')
      if (ctrl_loptc == 0 .and. io_help == 0) goto 279 !don't skip help mode

      nm='OPTICS_LTET'; call gtv(trim(nm),tksw2(prgn,nm),optic_ltet,def_i4=3,
     .  note='Integration mode:'//
     .   '%N%3f0: Methfessel-Paxton sampling'//
     .   '%N%3f1: standard tetrahedron package'//
     .   '%N%3f2: same as 1'//
     .   '%N%3f3: enhanced tetrahedron integration'//
     .   '%N%2f-1: on-the-fly MP sampling (conserves memory)')

      nm='OPTICS_PART'; call gtv(trim(nm),tksw2(prgn,nm),optic_mode1,
     .  def_i4=0,note='Decomposition of epsilon or jdos'//
     .  '%N%3f0:  No decomposition'//
     .  '%N%3f1:  Resolve eps or dos into (occ,unocc) contributions'//
     .  '%N%3f2:  Resolve eps or dos by k'//
     .  '%N%3f3:  Both 1 and 2'//
     .  '%N%3f    Add 10 to save popt as a binary file')

      nm='OPTICS_WINDOW'; call gtv(trim(nm),tksw2(prgn,nm),optic_window(1:2),
     .  def_r8v=(/0d0,1d0/),note='Energy window over which to calc. eps')
      nm='OPTICS_NPTS'; call gtv(trim(nm),tksw2(prgn,nm),optic_ne,nout=nout,
     .  note='Number of energy points',or=T)
      if (nout > 0) then       !nout>0 if NPTS was read
        if (optic_ne < 2 .and. io_help == 0) call rx('rdctrl: OPTICS_NPTS must be at least 2')
        optic_dw(1) = (optic_window(2)-optic_window(1))/(optic_ne-1)
      else
        nm='OPTICS_DW'; call gtv(trim(nm),tksw2(prgn,nm),optic_dw(1:2),note=
     .    'Energy mesh spacing.'//
     .    '%N%3fOne argument supplied => uniform spacing DW in WINDOW'//
     .    '%N%3fOptional 2nd argument OMG => space points '//
     .    'quadratically as:'//
     .    '%N%5fE(i) = WINDOW(1) + DW*(i-1) + DW**2/OMG/2*(i-1)**2')
        if (optic_dw(2) == NULLI) optic_dw(2) = 0
        if (optic_dw(1) <= 0 .and. io_help == 0) then
          call rx1('READCTRL: mesh spacing OPTICS_DE=%d '//
     .      '(must be positive)',optic_dw(1))
        endif
      endif

      nm='OPTICS_FILBND'; call gtv(trim(nm),tksw2(prgn,nm),optic_ocrng,
     .  def_i4v=izerv,note='Occ. energy bands from which to calc. eps.  Range of states for single DOS.')
      nm='OPTICS_EMPBND'; call gtv(trim(nm),tksw2(prgn,nm),optic_unrng,
     .  def_i4v=izerv,note='Unocc energy bands from which to calc. eps.  Not relevant for single DOS.')
      nm='OPTICS_NMP'; call gtv(trim(nm),tksw2(prgn,nm),optic_Nmp,def_i4=NULLI,
     .  note='If present, supersedes BZ_N in the optics calculation')
      nm='OPTICS_W'; call gtv(trim(nm),tksw2(prgn,nm),optic_w,def_r8=NULLR,
     .  note='If present, supersedes BZ_W in the optics calculation')
      nm='OPTICS_FERMI'; call gtv(trim(nm),tksw2(prgn,nm),optic_imref(1),def_r8=NULLR,
     .  note='If present, supersede calculated Fermi level with this value in the optics calculation')

      optic_axes = NULLI
      nm='OPTICS_CHI2'; call gtv(trim(nm),tksw2(prgn,nm),nono,
     .  Texist=ltmp,note='Parameters second harmonic generation.'//
     .    '  If present then:')
      if (ltmp) then   ! Only read if OPTICS_CHI2 present or help
        nm='OPTICS_CHI2_NCHI2'; call gtv(trim(nm),tksw2(prgn,nm),
     .    optic_nchi2,def_i4=0,note=
     .    'number of direction vectors for which to calc. chi2')
        nm = 'OPTICS_CHI2_AXES'
        nn = 3*optic_nchi2
        if (io_help == 1 .and. sw /= 2) nn = NULLI
        call gtv(trim(nm),tksw2(prgn,nm),optic_axes(1:nn),nmin=nn,
     .    note='a,b,c direction vectors for each of the nchi2 sets')
      endif
      nm='OPTICS_ESCISS'; call gtv(trim(nm),tksw2(prgn,nm),optic_esciss,def_r8=0d0,
     .  note='Scissors operator (energy added to unoccupied bands)')

      nm='OPTICS_MEFAC'; call gtv(trim(nm),tksw2(prgn,nm),optic_mefac,def_i4=0,
     .  note='If 1, include dSig/dk nonlocal sigm'//
     .  '%N   If 2, approximate correction to Im eps from nonlocal sigm using LDA eigenvalues')

      nm='OPTICS_FFMT'; call gtv(trim(nm),tksw2(prgn,nm),optic_ffmt,def_i4=0,
     .  note='Governs formatting of optics file'//
     .  '%N%3f0 floating format'//
     .  '%N%3f1 exponential format')

      nm='OPTICS_IQ'; call gtv(trim(nm),tksw2(prgn,nm),optic_iq,def_i4v=(/0,0,0/),
     .  note='q vector for JDOS(q), in multiples of qlat/NKABC')

      nm='OPTICS_ESMR'; call gtv(trim(nm),tksw2(prgn,nm),optic_esmr,def_r8=5d-2,
     .  note='Energy smearing width for determining (occ,unocc) window.'//
     .  '%N%3fStates are excluded for which occ<EF-ESMR or unocc>EF+ESMR.')
      nm='OPTICS_ALLTRANS'; call gtv(trim(nm),tksw2(prgn,nm),optic_alltrans,
     .  def_lg=F, note= 'Do not limit allowed transitions to (occ<EF-ESMR and unocc>EF+ESMR)')
      nm='OPTICS_FERMI'; call gtv(trim(nm),tksw2(prgn,nm),optic_imref(1),def_r8=NULLR,
     .  note='Supersede calculated Fermi level with this value for calculating optics')

      if (io_help==0 .and. abs(ctrl_loptc) /= 8 .and. abs(ctrl_loptc) /= 9) goto 279
      if (io_help/=0) call info0(2,1,0,' ...  The following apply to MODE=8 or 9:')

      optic_imref = NULLR
      nm='OPTICS_IMREF'; call gtv(trim(nm),tksw2(prgn,nm),optic_imref(1:2),
     .  def_r8v=(/NULLR,NULLR/),
     .  nmin=1, note='Quasi-Fermi levels for occ and unocc states.'//
     .  '%N%3fIf you use IMREF= or IMREF=-99999, crystal Ef will be substituted for both.')
      if (optic_imref(1) /= NULLR .and. optic_imref(2) == NULLR) optic_imref(2)=optic_imref(1)
      nm='OPTICS_KT'; call gtv(trim(nm),tksw2(prgn,nm),optic_kt,
     .  note='Temperature for Fermi functions (Ry)')

  279 continue                  ! end of optics

C --- Dynamics ---
      if (tkswp(prgn,express,'DYN') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for dynamics and statics ---')

      nm='DYN_NIT'; call gtv(trim(nm),tksw2(prgn,nm),nitmv,def_i4=1,
     .  note='maximum number of relaxation steps (statics)'//
     .  ' or time steps (dynamics)')

      nm='DYN_SSTAT'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,Texist=ltmp,
     .  note='Parameters for spin statics:')
      if (io_help /= 0 .or. ltmp) then
        nm='DYN_SSTAT_MODE'; call gtv(trim(nm),tksw2(prgn,nm),sdmod,def_i4=-1,note=
     .    '%b<0:  No spin statics or dynamics'//
     .    '%N%3f0:  output Euler angles generated from density matrix'//
     .    '%N%3f1:  relax along force'//
     .    '%N%3fAdd 10   to mix angles independently of P,Q'//
     .    '%N%3fAdd 1000 to suppress self-const in P,Q')
        nm='DYN_SSTAT_SCALE'; call gtv(trim(nm),tksw2(prgn,nm),sdprm(1),
     .    def_r8=1d0,note='Scale factor amplifying magnetic forces')
        if (io_help /= 0 .or. mod(sdmod,10) == 2) then
        nm='DYN_SSTAT_TAU'; call gtv(trim(nm),tksw2(prgn,nm),sdprm(2),
     .    note='Time step (only read if Landau-Gilbert dynamics specified)')
        endif
        nm='DYN_SSTAT_MAXT'; call gtv(trim(nm),tksw2(prgn,nm),sdprm(4),
     .    def_r8=0d0,note='Maximum allowed change in angle')
        nm='DYN_SSTAT_ETOL'; call gtv(trim(nm),tksw2(prgn,nm),sdprm(5),
     .    def_r8=0d0,note='Set tau=0 this iter if etot-ehf>ETOL')
      endif
      lncol16 = sdmod >= 0

      nm='DYN_MSTAT'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,Texist=ltmp,
     .  note='Parameters for molecular statics')
      mdprm(1) = 0
      if (io_help /= 0 .or. ltmp) then
      nm='DYN_MSTAT_MODE'; call gtv(trim(nm),tksw2(prgn,nm),i0,
     .  def_i4=0,note=
     .  '0: no relaxation  '//
     .  '4: conjugate gradients  '//
     .  '5: Fletcher-Powell  '//
     .  '6: Broyden')
      mdprm(1) = i0
      nm='DYN_MSTAT_HESS'; call gtv(trim(nm),tksw2(prgn,nm),ltmp,
     .  def_lg=T,note='Read hessian matrix')
      mdprm(2) = isw(ltmp)
      nm='DYN_MSTAT_XTOL'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(3),
     .  def_r8=1d-3,note=
     .  'Convergence criterion in displacements'//
     .  '%N%3fXTOL>0: use length;  <0: use max val;  =0: do not use')
      nm='DYN_MSTAT_GTOL'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(4),
     .  def_r8=0d0,note=
     .  'Convergence criterion in gradients'//
     .  '%N%3fGTOL>0: use length;  <0: use max val;  =0: do not use')
      nm='DYN_MSTAT_STEP'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(5),
     .  def_r8=0.015d0,note='Initial (and maximum) step length')
      nm='DYN_MSTAT_NKILL'; call gtv(trim(nm),tksw2(prgn,nm),i0,
     .  def_i4=0,note='Remove hessian after NKILL iter')
      mdprm(6) = i0
      nm='DYN_MSTAT_PDEF'; call gtv(trim(nm),tksw2(prgn,nm),lat_defm,
     .  def_r8v=zerov,note='Lattice deformation modes')
      endif

C --- if DYN_MSTAT_MODE < 4 then attempt to read MD ---
      if (mdprm(1) < 4) then
      nm='DYN_MD'; sw=tksw2(prgn,nm)
      call gtv_none(trim(nm),sw,nono,Texist=ltmp,note='Parameters for molecular dynamics')
      if (io_help /= 0 .or. ltmp) then
      nm='DYN_MD_MODE'; call gtv(trim(nm),tksw2(prgn,nm),i0,def_i4=0,note=
     .  '0: no MD '//
     .  '1: NVE '//
     .  '2: NVT '//
     .  '3: NPT')
      mdprm(1) = i0
      nm='DYN_MD_TSTEP'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(2),
     .  def_r8=fs,note='Time step (a.u.)')
      nm='DYN_MD_TEMP'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(3),
     .  def_r8=300d0*degK,note='Temperature (a.u.)')
      nm='DYN_MD_TAUP'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(4),
     .  def_r8=10d0*fs,note='Thermostat relaxation time (a.u.)')
      nm='DYN_MD_TIME'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(5),
     .  def_r8=1d6*fs,note='Total MD time (a.u.)')
      nm='DYN_MD_TAUB'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(6),
     .  def_r8=100d0*fs,note='Barostat relaxation time (a.u.)')
      nm='DYN_MD_P'; call gtv(trim(nm),tksw2(prgn,nm),mdprm(7),
     .  def_r8=0d0,note='External pressure')
      if (tksw2(prgn,'DYN_MD') /= 2) call info2(10,0,0,
     .  ' NB: 1 deg.K = %e a.u.; 1 fs = %d a.u.',degK,fs)
      endif
      if (sw /= 2 .and. mdprm(1) > 3)
     .  call rx('M_RDCTRL: illegal "MODE" parameter in token "DYN_MD"')
      endif

C Content of mdprm(i)
C   i = 1: statics/dynamics switch
C          0 do neither
C          1-3 MD (1: NVE, 2: NVT; 3: NPT)
C          >3 statics (4: conjugate gradients, 5: Fletcher-Powell, 6: Broyden)
C          >100 statics with lattice deformation
C   i = 2: (stat) 1 read hessian matrix
C          (dyn)  MD timestep
C   i = 3: (stat) relaxation x-tolerance
C          (dyn)  temperature
C   i = 4: (stat) relaxation g-tolerance
C          (dyn)  thermostat relaxation time
C   i = 5: (stat) step length
C          (dyn)  total MD time
C   i = 6: (stat) Remove hessian after this many steps
C          (dyn)  barostat relaxation time

      nm='DYN_SDYN'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,Texist=lsdyn,
     .  note='Parameters for spin dynamics')
      if (io_help /= 0 .or. lsdyn) then
      nm='DYN_SDYN_KT'; call gtv(trim(nm),tksw2(prgn,nm),move_kt,
     .  note='Temperature, in a.u.')
      nm='DYN_SDYN_TS'; call gtv(trim(nm),tksw2(prgn,nm),move_ts,
     .  note='Time step in a.u.')
      nm='DYN_SDYN_TEQU'; call gtv(trim(nm),tksw2(prgn,nm),move_tsequ,
     .  def_r8=0d0,note='Equilibration time, a.u.')
      nm='DYN_SDYN_TTOT'; call gtv(trim(nm),tksw2(prgn,nm),move_tstot,
     .  note='Duration of total simulation, in a.u.')
      endif

      nm='DYN_MMHAM'; call gtv(trim(nm),tksw2(prgn,nm),mmham,
     .  note='Rules for micromagnetics hamiltonian')

C ... Global daemons
      nm='DYN_GD'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,Texist=ltmp,
     .  note='Parameters for Global daemons in theromstats')
C    if (io_help /= 0 .or. ltmp) then
      nm='DYN_GD_NTHERM'; call gtv(trim(nm),tksw2(prgn,nm),gd_nmodt,
     .  def_i4=1,note='Number of GD thermostats')
      nn = gd_nmodt
      if (io_help /= 0 .or. tksw2(prgn,nm) >= 2) then
        nn = NULLI
      else
        call sanrg(io_help==0,nn,0,3,'rdctrl','GD_NTHERM')
      endif
      nm='DYN_GD_MODET'; call gtv(trim(nm),tksw2(prgn,nm),gd_modt(1:max(1,nn)),
     .  nmin=nn,def_i4v=(/31,32,33/),note='Thermostat mode(s)')
      nm='DYN_GD_CT'; call gtv(trim(nm),tksw2(prgn,nm),gd_ct(1:max(1,nn)),nmin=nn,
     .  def_r8v=(/1d0,1d0,1d0/),note='Corresponding thermostat coefficient(s)')
C     endif

C ... Bulirsch-Stoer integration parmaters
      nm='DYN_BSINT'; call gtv_none(trim(nm),tksw2(prgn,nm),nono,
     .  Texist=lbsprm,note='Parameters for Bulirsch-Stoer integration')
      if (io_help /= 0 .or. lbsprm) then
      nm='DYN_BSINT_NEW'; call gtv(trim(nm),tksw2(prgn,nm),prmint_new,
     .  def_lg=T,note='Start new SD run, Bulirsch-Stoer integration')
      nm='DYN_BSINT_TOL'; call gtv(trim(nm),tksw2(prgn,nm),prmint_tol,
     .  note='Tolerance in numerical integration')
      nm='DYN_BSINT_TS0'; call gtv(trim(nm),tksw2(prgn,nm),prmint_ts0,
     .  def_r8=0d0,note='Minimum time step in units of TS')
      nm='DYN_BSINT_MX'; call gtv(trim(nm),tksw2(prgn,nm),prmint_mx,
     . def_i4=7,note='Maximum order of rational function extrapolation')
      nm='DYN_BSINT_MI'; call gtv(trim(nm),tksw2(prgn,nm),prmint_mi,
     .  def_i4=11,note='Maximum number of midpoint rules to invoke')
      nn = prmint_mi
      if (io_help /= 0 .or. tksw2(prgn,nm) >=  2) then
        nn = NULLI
      else
        call sanrg(io_help==0,nn,2,11,'rdctrl','BSINT_MI')
      endif
      nm='DYN_BSINT_NSEQ'; call gtv(trim(nm),tksw2(prgn,nm),
     .  prmint_nseq(1:max(1,nn)),
     .  nmin=nn,def_i4v=(/2,4,6,8,12,16,24,32,48,64,96/),
     .  note='Sequence of number of midpoint divisions')
      endif
cxxx int
C      sdyn_sdmod= sdyn(2)       !xxxxx
C      sdyn_fscal= sdyn(3)
C      sdyn_tau = sdyn(4)
C      sdyn_NoseEtot0= sdyn(5)
C      sdyn_maxtheta= sdyn(6)
C      sdyn_etol= sdyn(7)
      endif ! DYN

C --- GW ---
      if (tkswp(prgn,express,'GW') < 2) then
      if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for GW ---')

      nm='GW_NKABC'; call gtv(trim(nm),tksw2(prgn,nm),gw_nabc,nout=nout,
     .  note='No. qp along each of 3 lattice vectors.'//
     .  '%N%3fSupply one number for all vectors or a separate '//
     .  'number for each vector.')
      call fill3in(nout,gw_nabc)
      if (io_help == 0) call fill3inl(gw_nabc,plat,io_show /= 0,'NKABC')
      if (gw_nabc(1) == NULLI .or. gw_nabc(1) == 0) gw_nabc = bz_nabc

      nm='GW_BZJOB';call gtv(trim(nm),tksw2(prgn,nm),gw_lshft,nout=nout,def_i4v=izerv(1:1),note=
     .  '0 centers BZ mesh at origin, 1 centers off origin'//
     .  '%N%3fSupply one number for all vectors or a separate '//
     .  'number for each vector.')
      call fill3in(nout,gw_lshft)

      if (io_show+io_help/=0 .and. tksw2(prgn,'GW_CODE')<2) call info0(2,1,0,' ... Parameters for GW driver')

      nm='GW_CODE'; call gtv(trim(nm),tksw2(prgn,nm),gw_code,
     .  def_i4=2,note='Driver for one of these GW codes:'//
     .  '%N%3f0 fpgw (v033a5)'//
     .  '%N%3f1 spex (v02.04)'//
     .  '%N%3f2 fpgw (Sep12)'//
     . '%N%2f10 DMFT')

      if (io_show+io_help/=0 .and. tksw2(prgn,'GW_MKSIG')<2)
     .  call info0(2,1,0,' ... Parameters for setup job=-1')

      nm='GW_MKSIG'; call gtv(trim(nm),tksw2(prgn,nm),gw_mksig,def_i4=3,note=
     .  'Make full sigma(k) for self-consistent calculations:'//
     .  '%N%3f0 do not create self-energy matrix sigma'//
     .  '%N%3f1 Make sigma_nn''(E_f)'//
     .  '%N%3f2 Make sigma_nn''(E_n/2+E_n''/2)'//
     .  '%N%3f3 Make sigma_nn''(E_n)/2 + sigma_nn''(E_n'')/2')

C      nm='GW_EVECL'; call gtv(trim(nm),tksw2(prgn,nm),gw_evecl,
C     .  def_i4=0,note='If nonzero, write LMTO-LDA evecs to disk')

C      nm='GW_ECUTPB'; call gtv(trim(nm),tksw2(prgn,nm),gw_ecutpb,
C     .  def_r8v=[-6d0,0d0],note=
C     .  'Include core states in product basis with energy > ECUTPB'//
C     .  '%N%3fOptional second element : Include core states charge > ECUTPB(2)')
      nm='GW_ECUTPB'; call gtv(trim(nm),tksw2(prgn,nm),gw_ecutpb,
     .  def_r8=-6.5d0,note='Include core states in product basis with energy > ECUTPB')
      nm='GW_GCUTB'; call gtv(trim(nm),tksw2(prgn,nm),gw_gcutb,
     . def_r8=2.7d0,note='G-vector cutoff for basis envelope functions')
      nm='GW_GCUTX'; call gtv(trim(nm),tksw2(prgn,nm),gw_gcutx,
     .  def_r8=2.2d0,note='G-vector cutoff for interstitial part of response function')
      nm='GW_ECUTS'; call gtv(trim(nm),tksw2(prgn,nm),gw_ecuts,def_r8=sigp_emax+.5d0,
     .  note='(self-consistency) max energy for which to calc. sigma')
      nm='GW_NBAND'; call gtv(trim(nm),tksw2(prgn,nm),gw_nband,def_i4=-NULLI,
     .  note='No. bands to include in calc. response function')
      nm='GW_NIME'; call gtv(trim(nm),tksw2(prgn,nm),gw_nime,def_i4=6,
     .  note='Number of imaginary energy points for energy contour')
      nm='GW_DELRE'; call gtv(trim(nm),tksw2(prgn,nm),gw_delre,def_r8v=(/.01d0,.04d0/),note=
     .  'Energy mesh spacing on real axis.'//
     .  '%N%3fSecond element is coefficient to nonlinear scaling')
      nm='GW_DELTA'; call gtv(trim(nm),tksw2(prgn,nm),gw_deltax,def_r8=-1d-4,note=
     .  'Imaginary energy added to denominator in generation of Pi.'//
     .  '%N%3fNot used in Faleev''s fast integration.'//
     .  '%N%3fDelta should be negative.')
      nm='GW_DELTAW'; call gtv(trim(nm),tksw2(prgn,nm),gw_deltaw,def_r8=0.02d0,note=
     .  'Width for finite difference in diff. of sigma for Z factor')
      nm='GW_GSMEAR'; call gtv(trim(nm),tksw2(prgn,nm),gw_gsmear,def_r8=3d-3,note=
     .  'Broadening in pole of GF for computation of sigma')
      nm='GW_PBTOL'; call gtv(trim(nm),tksw2(prgn,nm),gw_pbtol,
     .  def_r8v=[3d-4,3d-4,1d-3],
     .  note='retain prod-basis with overlap eigenvalues>PBTOL')
      nm='GW_QOFFP'; call gtv(trim(nm),tksw2(prgn,nm),gw_qoffp,def_r8=1d0,note=
     .  'k-point offset parameter for special BZ integration for xi')
      nm='GW_USEBSEW'; call gtv(trim(nm),tksw2(prgn,nm),gw_usebsew,def_i4=0,
     .  note='include BSE contributions to W')

      endif ! GW

C --- DMFT ---
      if (tksw2(prgn,'DMFT') < 2) then
        if (io_show+io_help/=0) call info0(2,1,0,' --- Parameters for DMFT ---')

        nm='DMFT_NKABC'; call gtv(trim(nm),tksw2(prgn,nm),dmft_nabc,nout=nout,Texist=ltmp,
     .    note='No. qp along each of 3 lattice vectors for DMFT driver.%N%3f'//
     .    'If not present, substitute BZ_NKABC.')
        if (dmft_nabc(1) == NULLI .or. dmft_nabc(1) == 0) dmft_nabc = bz_nabc
        call fill3in(nout,dmft_nabc)
        if (io_help == 0) call fill3inl(dmft_nabc,plat,io_show /= 0,'NKABC')

        nm = 'DMFT_NLOHI'; sw = tksw2(prgn,nm)
        if (sw < 2) then
          nout = 0
          call gtv(trim(nm),sw,dmft_nlohi,nout=nout,
     .      note='first, last eval to include in projector',or=T)
          if (nout /= 2) then
            nm = 'DMFT_WLOHI'; sw = tksw2(prgn,nm)
            call gtv(trim(nm),sw,dmft_wlohi,nout=nout,note='lower, upper'//
     .        ' bound to frequency to include in projector, relative to Ef')
          endif
        endif

        nm = 'DMFT_PROJ'; call gtv(trim(nm),tksw2(prgn,nm),dmft_proj,def_i4=2,note='DMFT projector type : '//
     .    '%N%3f1 eq (17) in PhysRevB.81.195107'//
     .    '%N%3f2 eq (22) in PhysRevB.81.195107')


        nm = 'DMFT_KNORM'; call gtv(trim(nm),tksw2(prgn,nm),dmft_knorm,
     .    def_i4=0,note='How local projectors are normalized:'//
     .    '%N%3f0 k-independent normalization'//
     .       '%N%3f1 k-dependent normalization')
        
        nm = 'DMFT_NOMF'; call gtv(trim(nm),tksw2(prgn,nm),dmft_nomf,
     .    def_i4=0,note='number of fermionic frequencies for the susceptibility')
        nm = 'DMFT_NOMB'; call gtv(trim(nm),tksw2(prgn,nm),dmft_nomb,
     .    def_i4=0,note='number of fermionic frequencies for the susceptibility')

C       This is not read for now
        nm = 'DMFT_IMW'; call gtv(trim(nm),tksw2(prgn,nm),dmft_imw,def_lg=T,
     .    note='DMFT on imaginary axis')

        nm = 'DMFT_BROAD'; call gtv(trim(nm),tksw2(prgn,nm),dmft_broad,def_r8=1/400d0,
     .    nout=nout,note='Broadening of sigma in eV (when computed on real axis only)')

        nm = 'DMFT_BETAR'; call gtv(trim(nm),tksw2(prgn,nm),dmft_beta,def_r8=38.61d0,or=T,
     .    nout=nout,note='Inverse temperature, in Ry^-1')
        if (nout /= 0) dmft_beta = dmft_beta * RyeV
        if (nout == 0) then     !nout=-1 if sw=2; otherwise nout=0 unless prior tag was read
          nm = 'DMFT_BETAK'; call gtv(trim(nm),tksw2(prgn,nm),dmft_beta,
     .      def_r8=0.00332715d0,or=T,nout=nout,note='Inverse temperature, in Kelvin^-1')
          if (nout /= 0) dmft_beta = dmft_beta * RyeV/kboltzeV
        endif
        if (nout == 0) then     !nout=-1 if sw=2; otherwise nout=0 unless prior tag was read
          nm = 'DMFT_BETA'; call gtv(trim(nm),tksw2(prgn,nm),dmft_beta,
     .      nout=nout,note='Inverse temperature, in eV^-1')
        endif

        nm = 'DMFT_NOMEGA'; call gtv(trim(nm),tksw2(prgn,nm),dmft_nomg,
     .    def_i4=2000,note='Number of points on frequency mesh')

C   ... Count number of DMFT correlated blocks, parameters for dimensioning
        nm='DMFT_BLOCK'; sw = tkswp(prgn,express,nm); if (sw >= 2) goto 109
        dmft_maxl = 0; ncix = 1; nicix = 1
        if (io_help == 0) then
          dmft_maxl = -1  ! Maximum l in any block
          nicix = 0       ! total number of independent sites
          ncix = 0        ! total number of sites participating in DMFT
          call pshpr(1)
          j = 0; ltmp = .true.
          do  while (ltmp)
            j = j+1; jj= (/1,j/)
            nm='DMFT_BLOCK'; call gtv(trim(nm),0,nono,Texist=ltmp,cindx=jj)
            if (.not. ltmp) exit
            nm='DMFT_BLOCK_L'; call gtv(trim(nm),tksw2(prgn,nm),i,cindx=jj) ! largest l
            dmft_maxl = max(dmft_maxl,i)
            nm='DMFT_BLOCK_SITES'; call gtv(nm,tksw2(prgn,nm),ivec,nout=nout,cindx=jj) ! equivalent blocks
            nicix = nicix+1
            ncix = ncix+nout
          enddo
          call poppr
        endif
        if (ncix <= 0) goto 109
        dmft_maxdim = 2*dmft_maxl+1

        if (io_show+io_help/=0 .and. ncix>0) call info0(2,1,0,
     .    ' --- Parameters for DMFT correlated subblocks ---')
        if (io_help == 0 .and. io_show>0) then
          call info(1,0,0,' ... Found %i cix blocks, %i independent',ncix,nicix)
        endif
        if (io_help >= 1) then
          write(*,309)
  309     format(/' - ',
     .      'The following tokens are read for each correlated subblock. ',
     .      'Data sandwiched'/3x,'between successive occurences of ',
     .      'token BLOCK apply to one DMFT correlated block.')
        endif

        if (.not. allocated(dmft_l)) then
          allocate(dmft_ib(ncix),dmft_icix(ncix),ipsl(ncix*dmft_maxdim**2))
          allocate(dmft_sigind(dmft_maxdim,dmft_maxdim,nicix,nsp))
          allocate(dmft_l(nicix),dmft_umode(nicix),dmft_qsplit(nicix),dmft_ndsigi(0:nicix))
C         allocate(dmft_ndim(nicix))
          call iinit(dmft_sigind,size(dmft_sigind))
        endif

C   ... DMFT correlated block data, one pass for each block
        ncix = 0                ! running index to current cix block
        jp = 0                  ! At end of loop, highest value sigind from prior block
        do  nn = 1, nicix       ! nicix should be correct number of independent cix blocks

          dmft_ndsigi(nn-1) = jp
          jj= (/1,nn/)
          if (io_help /= 0) then
            write(stdo,'(1x)')
          elseif (io_help == 0 .and. io_show>0) then
            call info(1,0,0,' ... Reading subblock %i for independent cix block %i',ncix+1,nn)
          endif
          nm='DMFT_BLOCK'; call gtv(trim(nm),0,nono,Texist=ltmp,cindx=jj)
          if (.not. ltmp) call rxi('readctrl failed to read DMFT_BLOCK for block',nn)

          nm='DMFT_BLOCK_L'; call gtv(trim(nm),tksw2(prgn,nm),dmft_l(nn),cindx=jj)   ! block l
          nm='DMFT_BLOCK_QSPLIT'; call gtv(trim(nm),tksw2(prgn,nm),dmft_qsplit(nn),def_i4=2,cindx=jj)   ! Haule's qsplit
          nm='DMFT_BLOCK_SITES'; call gtv(nm,tksw2(prgn,nm),ivec,nout=k,cindx=jj) ! sites with this block
          do  i = 1, k          ! k = number of sites
            ncix = ncix+1
            sw = isign(1,ivec(i)); ivec(i) = iabs(ivec(i))
            dmft_ib(ncix) = ivec(i)
            dmft_icix(ncix) = nn*sw
          enddo
          i0 = (2*dmft_l(nn)+1); if (io_help == 1) i0 = NULLI
          nm='DMFT_BLOCK_SIDXD'; call gtv(nm,tksw2(prgn,nm),ivec(1:max(i0,1)),
     .      nmin=i0,nout=nout,cindx=jj,or=T,note='(diagonal Sigma only) list components of Sigma to calculate')
          if (nout /= 0) then
            call icopy(nout,ivec,1,dmft_sigind(1,1,nn,1),dmft_maxdim+1)
          else                  !nout=-1 if sw=2; otherwise nout=0 unless prior tag was read
            i = i0**2; if (io_help == 1) i = NULLI
            nm='DMFT_BLOCK_SIDXM'; call gtv(nm,tksw2(prgn,nm),ivec(1:i),nmin=i,nout=nout,cindx=jj,
     .        note='(matrix sigma) list components of Sigma to calculate')
            ii = 0
            do  i = 1, i0
            do  j = 1, i0
              ii = ii+1
              dmft_sigind(i,j,nn,1) = ivec(ii)
            enddo
            enddo
          endif
C         Sanity check: indices must be contiguous, and contiguous with prior block
          if (io_help == 0) then
            ii = 0
            do  i = 1, i0
            do  j = 1, i0
              ii = ii+1
              ivec(ii) = dmft_sigind(i,j,nn,1)
            enddo
            enddo
            call ivheap(1,i0**2,ivec,ipsl,1)
            if (ivec(ipsl(i0**2))<=0) call rxi('readctrl no nonvanishing Sigma index for block',nn)
            ii = 0              ! Accumulate total number of nonvanish elements this block
            do  i = 1, i0**2
              if (ivec(ipsl(i)) == 0) cycle
C             print *, nn, ivec(ipsl(i)), jp
              if (ivec(ipsl(i)) /= jp .and. ivec(ipsl(i)) /= jp+1)
     .        call rx1('sigind element %i not contiguous',ivec(ipsl(i)))
              jp = ivec(ipsl(i))
              ii = ii+1
            enddo
C           Accumulate number of nonzero elements for all sites from this independent cix block
            dmft_nzsig = dmft_nzsig + k*ii
          endif
C         nm='DMFT_BLOCK_UMODE'; call gtv(nm,tksw2(prgn,nm),ivec(1:i),nmin=i,nout=nout,cindx=jj,
          nm='DMFT_BLOCK_UMODE'
C         call gtv(trim(nm),tksw2(prgn,nm),dmft_umode(nn),nmin=0,nout=nout,Texist=ltmp,cindx=jj,note='xx')
          call gtv(trim(nm),tksw2(prgn,nm),dmft_umode(nn),nmin=0,nout=nout,cindx=jj,or=T,
     .      note='Specify by integer the meaning of DMFT hubbard U for this block as follows'//
     .      '%N%3f1s digit specifies meaning of U'//
     .      '%N%5f0: u(1) = Hubbard U, u(2) = Hubbard J'//
     .      '%N%5f1: u(1) = F0, u(2) = F2, u(3) = F4'//
     .      '%N%5f2: u(1) = screening length, Yukawa potential'//
     .      '%N%3f10s digit specifies approximation for U'//
     .      '%N%5f0: density-density'//
     .      '%N%5f1: full matrix U'//
     .      '%N%3f100s digit specifies frequency dependence of U'//
     .      '%N%5f0: U is static'//
     .      '%N%5f1: U is dynamic')
          if (nout == 0) then   ! No integer read; try reading character string
            dmft_umode(nn) = 10
            call gtv(trim(nm),tksw2(prgn,nm),outs(2:),nmin=10,cindx=jj,
     .        nout=nout,
     .        note='Specify by a string information as follows.'//
     .        '%N%3fChoose one string from each line and separate arguments by "~".'//
     .        '%N%3fuj | slater | yukawa'//
     .        '%N%3fdensity | full')
            if (nout > 0) then
              outs(1:1) = '~'
              ii = dmft_umode(nn)
              if (wordsw(outs,'~ ','uj','~ ',i) > 0)      ii = ii +   0 -   1*(getdig(ii,0,10))
              if (wordsw(outs,'~ ','slater','~ ',i) > 0)  ii = ii +   1 -   1*(getdig(ii,0,10))
              if (wordsw(outs,'~ ','yukawa','~ ',i) > 0)  ii = ii +   2 -   1*(getdig(ii,0,10))
              if (wordsw(outs,'~ ','density','~ ',i) > 0) ii = ii +   0 -  10*(getdig(ii,1,10))
              if (wordsw(outs,'~ ','full','~ ',i) > 0)    ii = ii +  10 -  10*(getdig(ii,1,10))
              if (wordsw(outs,'~ ','static','~ ',i) > 0)  ii = ii +   0 - 100*(getdig(ii,2,10))
              if (wordsw(outs,'~ ','dynamic','~ ',i) > 0) ii = ii + 100 - 100*(getdig(ii,2,10))
              dmft_umode(nn) = ii
            endif
          endif

C         dmft_ndim(nn) = jp - dmft_ndsigi(nn-1)
        enddo
        dmft_ndsigi(nicix) = jp

C       Count dmft_ncatom and perform sanity check: all cix blocks must have distinct (ib,l) pairs
        call ivheap(1,ncix,dmft_ib,ipsl,1)
        dmft_ncatom = 1
        do  i = 2, ncix
C          print *, i, ipsl(i-1), ipsl(i), dmft_ib(ipsl(i-1)), dmft_ib(ipsl(i))
C          print *, dmft_icix(ipsl(i-1)), dmft_icix(ipsl(i)),
C     .             dmft_l(dmft_icix(ipsl(i-1))), dmft_l(dmft_icix(ipsl(i)))
          if (dmft_ib(ipsl(i)) == dmft_ib(ipsl(i-1))) then
            if (dmft_l(dmft_icix(ipsl(i-1))) /= dmft_l(dmft_icix(ipsl(i)))) cycle
            call rx2('duplicate (ib,l) in cix blocks %i and %i',ipsl(i-1),ipsl(i))
          endif
          dmft_ncatom = dmft_ncatom + 1
        enddo
        deallocate(ipsl)

  109   continue

C       Debugging check that all processes parse input
C        if (sw < 2) then
CC        print *, 'procid,ncix', procid, ncix
CC        print *, 'procid,nicix', procid, nicix
CC        print *, 'procid,dmft_qsplit', procid, dmft_qsplit
CC        print *, 'procid,dmft_sigind', procid, dmft_sigind
C        call rx0('done')
C        endif

      endif

C --- Plane ---
      if (io_show+io_help/=0 .and. tksw2(prgn,'PLANE')<2)
     .  call info0(2,1,0,' --- Parameters defining a plane ---')
      nm='PLANE_NORMAL'; call gtv(trim(nm),tksw2(prgn,nm),
     .  slat_plat2(7:9),def_r8v=zerov,note='Plane normal')
      nm='PLANE_X'; call gtv(trim(nm),tksw2(prgn,nm),
     .  slat_plat2(1:3),def_r8v=zerov,note='Plane X-axis')
      nm='PLANE_NORMAL'; call gtv(trim(nm),tksw2(prgn,nm),
     .  slat_plat2(4:6),def_r8v=zerov,note='Plane Y-axis')

      if (express == -1) goto 48

      if (.not. associated(slabll,slabl)) deallocate(slabll)

      end subroutine readctrl

      subroutine readctrlpq(prgnam,nclasp,nl,nsp,pnu,qnu,pp,z,ves,initc,ics,clabl)
C- Read class data from ctrl file
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam
Ci   nclasp
Ci   clabl :class name
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Oct 07
C ----------------------------------------------------------------------
      use m_toksw
      use m_gtv
      implicit none
C ... Passed parameters
      character*(*):: prgnam
      integer:: nclasp,nl,nsp
C ... Class-related arrays
      integer:: initc(nclasp),ics(nclasp)
      real(8):: pnu(nl,nsp,nclasp),qnu(3,nl,nsp,nclasp),
     .  pp(6,nl,nsp,nclasp),ves(nclasp),z(nclasp)
      character*8 clabl(nclasp)
C ... Local parameters
      character(128) :: nm
      character(256) :: outs
      logical:: ltmp,rdves,cmdopt
      integer:: stdo,lgunit,i,j,ic,lp1,jj(2),nout,nlaj,nlaji,iv(100),ii,n,js
      real(8):: xv(100),xx,pi
      character(len=10) :: prgn(2)

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)

      prgn(1) = prgnam; prgn(2) = prgnam
      if (prgnam == 'LMPG') prgn(2) = 'LMGF'

C --- Start ---
      if (tksw2(prgn,'START') >= 2) goto 159
      if (io_show+io_help/=0) call
     .  info0(2,1,0,' --- Read potential data (ASA) ---')

      nm='START_BEGMOM'; call gtv(trim(nm),tksw2(prgn,nm),lasa3,def_i4=1,note=
     .  '0 Start from potential parameters (on file)'//
     .  '%N%3f1 Create potential and pot pars from P,Q'//
     .  '%N%3f2 Create pot pars from file potential'//
     .  '%N%3f3 Create potential from restart file')

      nm='START_FREE'; call gtv(trim(nm),tksw2(prgn,nm),lasa8,
     .  def_lg=F,note=
     .  '(Free atom mode): make self-consistent density for free atoms')

      nm='START_CNTROL'; call gtv(trim(nm),tksw2(prgn,nm),ltmp,def_lg=F,
     .  note='if F, the remainder of this category is ignored')
      if (io_help==0 .and. .not. ltmp) goto 159

      rdves = .false. ! gtv does not set anything when tksw == 2
      nm='START_RDVES'; call gtv(trim(nm),tksw2(prgn,nm),rdves,
     .  def_lg=F,note='Read Ves(RMT) from this category')
C     Subtract current ves from pp now; add updated ves later
      if (rdves) call shftpp(nclasp,nl*nsp,pp,ves,ves,T,F)

C ... Read data for each occurence of ATOM=
C     j = index to jth occurence; ic = class index
      outs = ' '
      j = 1
      do  while (j > 0)
        if (io_help /= 0) then
          write(stdo,"(/' ... For each class, read the following:')")
        elseif (io_show>0) then
          write(stdo,"(1x)")
C         call info(1,0,0,' ... Site %i',j,0)
        endif

        jj=(/1,j/)
        nm='START_ATOM'; call gtv(trim(nm),tksw2(prgn,nm),alabl,nmin=10,
     .    cindx=jj,nout=nout,note='Class label')
        if (nout == 0 .and. io_help == 0) then
          j = 0
          cycle
        endif

        if (io_help == 0) then
          ic = 0
          do  i = 1, nclasp
            if (alabl == clabl(i)) then
              ic = i
              call awrit0('%a'//alabl//'%a,',outs,len(outs),0)
              goto 881
            endif
          enddo
          call info0(10,0,0,' rdctrl: (warning) class "'
     .      //alabl//'%a" not in class table')
          j = j+1; cycle
  881     continue
          nlaj = 1+lmxa(ic)
          nlaji = nlaj
        else                    !help mode; cycle only once
          nlaj = 1
          nlaji = NULLI
          ic = 1
          j = -1
        endif

        nm='START_ATOM_P'; call gtv(trim(nm),tksw2(prgn,nm),
     .    xv(1:nlaj*nsp),nmin=nlaji*nsp,nout=nout,cindx=jj,note=
     .    'Methfessel''s Potential functions')
C       Copy to pnu, spins separately in case nlja<nl
        if (nout > 0) then
          initc(ic) = 2*(initc(ic)/2) + 1
          call dpzero(pnu(1,1,ic),nl*nsp)
          call dcopy(nlaj,xv,1,pnu(1,1,ic),1)
          if (nsp == 2)
     .    call dcopy(nlaj,xv(1+nlaj),1,pnu(1,2,ic),1)
C         Fill in values for l>lmxa
          if (.not. cmdopt('--keepnu',8,0,nm)) then
          call config(pnu(1,1,ic),-1,z(ic),iv,n)
          do  i = 1, nsp
          do  lp1 = 1, nl
            ii = iv(lp1)
            xx = pnu(lp1,i,ic)
C           If P=0, use free-electron value
            if (xx == 0) then
              xx = ii + 0.5d0 - datan(dble(lp1-1))/pi
C           If fractional part only, add integer part
            elseif (xx < 1) then
              xx = xx + ii
            endif
            pnu(lp1,i,ic) = xx
          enddo
          enddo
          endif
        endif

        nm='START_ATOM_Q'; call gtv(trim(nm),tksw2(prgn,nm),
     .    xv(1:3*nlaj*nsp),nmin=3*nlaji*nsp,nout=nout,cindx=jj,
     .    note='Energy moments of the spherical density')
C       Copy to qnu, spins separately in case nlja<nl
        if (nout > 0) then
          call dpzero(qnu(1,1,1,ic),3*nl*nsp)
          call dcopy(3*nlaj,xv,1,qnu(1,1,1,ic),1)
          if (nsp == 2)
     .    call dcopy(3*nlaj,xv(1+3*nlaj),1,qnu(1,1,2,ic),1)
        endif

        nm='START_ATOM_ENU'; call gtv(trim(nm),tksw2(prgn,nm),
     .    xv(1:nlaj*nsp),nmin=nlaji*nsp,nout=nout,cindx=jj,note=
     .    'Linearization energies')
        if (nout > 0) then
          call dcopy(nlaj,xv,1,pp(1,1,1,ic),6)
          if (nsp == 2)
     .    call dcopy(nlaj,xv(1+nlaj),1,pp(1,1,2,ic),6)
        endif

        if (rdves .or. io_help /= 0) then
          if (io_help >= 1 .and. tksw2(prgn,'START_RDVES') < 2)
     .      call info0(2,0,0,
     .  ' ... The following tokens are sought if RDVES=T')
          nm='START_ATOM_V'; call gtv(trim(nm),tksw2(prgn,nm),
     .      xx,nmin=1,nout=nout,cindx=jj,note=
     .      'ASA potential at MT boundary')
          if (nout > 0) ves(ic) = xx
          nm='START_ATOM_DV'; call gtv(trim(nm),tksw2(prgn,nm),
     .      xx,nmin=1,nout=nout,cindx=jj,note=
     .      'Add to ASA potential at MT boundary')
          if (nout > 0) ves(ic) = ves(ic) + xx
        endif

        j = j+1
      enddo
      if (rdves) call shftpp(nclasp,nl*nsp,pp,ves,ves,F,T)

C ... Printout
      if (io_show /= 0 .and. io_help == 0) then
        call awrit0('%23o%0p%N rdctrl: read P,Q for %a%b%N P,Q,ppars:',
     .    outs,80,-stdo)
        do   j = 1, nclasp
          js = ics(j)
          nlaj = 1+lmxa(j)
          if (mod(initc(j),2) == 1) then
            print 335, 'class ', clabl(j), (pnu(lp1,1,j),lp1=1,nlaj)
            if (nsp == 2)
     .      print 335,' spin ', '2       ',(pnu(lp1,2,j),lp1=1,nlaj)
            print 335,'      ', 'Q       ',(qnu(1,lp1,1,j),lp1=1,nlaj)
            if (nsp == 2)
     .      print 335,' spin ', '2       ',(qnu(1,lp1,2,j),lp1=1,nlaj)
  335       format(4x,a,a,7f8.4)
          else
            print 335, 'class ', clabl(j) // '  missing pnu'
          endif
          if (mod(initc(j)/2,2) == 1) then
            print 336,'  enu       ',      (pp(1,lp1,1,j), lp1=1,nlaj)
            if (nsp == 2)
     .      print 335, ' spin ', '2       ',(pp(1,lp1,2,j),lp1=1,nlaj)
            print 336,'    c       ',      (pp(2,lp1,1,j), lp1=1,nlaj)
            if (nsp == 2)
     .      print 335, ' spin ', '2       ',(pp(2,lp1,2,j),lp1=1,nlaj)
            print 336,'srdel       ',      (pp(3,lp1,1,j), lp1=1,nlaj)
            if (nsp == 2)
     .      print 335, ' spin ', '2       ',(pp(3,lp1,2,j),lp1=1,nlaj)
  336       format(6x,a,7f8.4)
          else
            print 335, 'class ', clabl(j) // '  missing ppars'
          endif
          print *, ' '
        enddo
      endif

  159 continue                  ! End of START

      end subroutine readctrlpq

      end module

      subroutine fill3in(nin,res)
C- Fills res(2) or res(2:3) if res(1) or res(1:2) are given
      implicit none
      integer nin,res(3)
      if (nin==2) then
        res(3) = res(2)
      elseif (nin==1) then
        res(2:3) = res(1)
      endif
      end

      subroutine fill3inl(res,plat,lprt,labl)
C- Overwrite res(3)<0 and/or res(2)<0 with best approx to uniform spacing
C ----------------------------------------------------------------------
Ci Inputs
Ci   res
Ci   labl  :(printout only)
Ci   qlat  :set of 3 vectors
Cu Updates
Cu   10 Mar 14  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer res(3)
      logical lprt
      character labl*(*)
      double precision plat(3,3)
C ... Local parameters
      integer,parameter:: NULLI = -99999
      logical lnew
      integer i,res0(3)
      double precision qlat(3,3),xx,l123(3)
      procedure(real(8)) :: dlength

C     No information to determine anythng
      if (res(1) == 0 .or. res(1) == NULLI .or. all(plat == NULLI)) return

      call mkqlat(plat,qlat,xx)
      lnew = .false.; res0 = res

C     Length of each r.l.v.
      do  i  = 1, 3
        l123(i) = dlength(3,qlat(1,i),1)
      enddo

C     Infer each length from the volume of the microcell
      if (res(1) < 0) then
C       estimate for number of points
        xx = (-res(1)/(l123(1)*l123(2)*l123(3)))**(1d0/3)
        do  i  = 1, 3
          res(i) = max(nint(xx*l123(i)),1)
        enddo
        lnew = .true.
C     Infer lengths 2,3 as proportions of 1st length
      else
        do  i  = 2, 3
          if (res(i) < 0) then
            res(i) = max(nint((res(1)*l123(i))/l123(1)),1)
            lnew = .true.
          endif
        enddo
      endif

      if (lnew .and. lprt) then
        call info2(30,0,0,' set '//labl//' =%3:1i',res,0)
        if (res(1)*res(2)*res(3) <= 0) call rx('improper k-mesh')
      endif

      end

      subroutine toksw_init(debug)
C- Assembles list of program-specific tags (full name tokens)
C ----------------------------------------------------------------------
Ci Inputs
Ci   debug
Co Outputs
Co   swtok :assembles lists of tokens (tags) for each main program
Cr Remarks
Cr   Generates tags for each program. Routines tksw2 and tksw describe how
Cr   tags are parsed, and how to make a tag required, optional, ignored, or excluded.
Cr
Cr   For some main programs toksw_init returns all tags the main program will parse.
Cr     For example:  LM, LMF, LMGF
Cr   For others, the list modifies or augments a reference set of tags; for example
Cr     LMDOS tags add to or subtract from those in 'BASE'
Cr     LMFA  tags add to or subtract from those in 'LMF'
Cr
Cr   A tag, e.g. HAM_NSPIN, may have characters appended .e.g. HAM_NSPIN~
Cr    char   meaning
Cr     ~     tag is optional calling program doesn't require contents or
Cr           can use defaults
Cr     !     tag is excluded.  This is similar to the tag missing all together
Cr           but it differs in that further searching for the tag is precluded
Cr           (e.g. lmfa will includes tags in lmf unless they are excluded)
Cr     %     tag may be parsed in express mode, from express category
Cr           Input is required
Cr     &     tag may be parsed in express mode, from express category
Cr           Input is optional
Cr     #     tag may be parsed in express mode, from express category
Cr           Input is optional in express, required in the regular category
Cu Updates
Cu   19 Apr 16 Redesigned for express mode
Cu   18 Oct 14 Redesigned to enable optional inclusion of reference tags.
C ----------------------------------------------------------------------
      use m_toksw
      implicit none
C ... Passed parameters
      logical:: debug
C ... Local parameters
      character*(200):: io_sw, struc_sw, asa_start
      character*(220):: asa_sw
      character*(50) sscale
      character*(180) bz_sw
      character*(200):: dyn_mstat, dyn_md

C ... Store commonly used groups of tags into arrays
      call clear_swtok(debug)
      io_sw = ' EXPRESS~ HEADER~ CONST~ IO_SHOW~ IO_HELP~ IO_VERBOS~ '//
     .  'IO_IACTIV~ IO_TIM~ OPTIONS_Q~ IO_MAXMEM~ IO_SHOMEM~'
      struc_sw = ' STRUC% STRUC_NSPEC~ STRUC_NCLASS STRUC_FILE& STRUC_NBAS  '//
     .  'STRUC_PLAT STRUC_ALAT STRUC_DALAT~ STRUC_NL~ STRUC_SHEAR~ '//
     .  'STRUC_ROT~ STRUC_DEFGRD~ STRUC_STRAIN~ STRUC_ALPHA~ STRUC_LOADV~ '
      asa_sw = ' OPTIONS~ OPTIONS_ASA OPTIONS_ASA_ADNF~ '
     .  //'OPTIONS_ASA_NSPH~ OPTIONS_ASA_TWOC~ OPTIONS_ASA_CCOR~ '
     .  //'OPTIONS_ASA_GAMMA~ OPTIONS_ASA_ELIN~ OPTIONS_ASA_NEWREP~ '
     .  //'OPTIONS_ASA_NOHYB~ OPTIONS_ASA_MTCOR~ OPTIONS_ASA_QMT~ '
     .  //'HAM_FRZWF~'
      sscale = ' SPEC_SCLWSR~ SPEC_WSRMAX~ SPEC_OMAX1~ SPEC_OMAX2~ '
      asa_start = ' START~ START_BEGMOM~ START_CNTROL~ START_RDVES~ '//
     .  'START_ATOM~ START_ATOM_P~ START_ATOM_Q~ START_ATOM_ENU~ START_ATOM_V~ START_ATOM_DV~'
      dyn_mstat = ' DYN_MSTAT~ DYN_MSTAT_MODE DYN_MSTAT_HESS~ DYN_MSTAT_XTOL~'//
     .  ' DYN_MSTAT_GTOL~ DYN_MSTAT_STEP~ DYN_MSTAT_NKILL~ DYN_MSTAT_PDEF~ '
      dyn_md = ' DYN_MD~ DYN_MD_MODE DYN_MD_TSTEP~ DYN_MD_TEMP~ '//
     .  'DYN_MD_TAUP~ DYN_MD_TIME DYN_MD_TAUB~ DYN_MD_P~'
      bz_sw = ' BZ BZ_NKABC# BZ_BZJOB~ BZ_SAVDOS~ BZ_DOS~ BZ_NPTS~ BZ_DELEF~ BZ_EF0~ BZ_EFMAX~ BZ_GETQP~ '//
     .        ' BZ_PUTQP~ BZ_INVIT~ BZ_NOINV~ BZ_METAL& BZ_N~ BZ_W~ BZ_TETRA~ BZ_NEVMX~ BZ_ZBAK~ BZ_ZVAL~'

C --- BASE tags ---
      call nswadd()
      call tkadd(" BASE::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))

      call tkadd(" EWALD~ EWALD_AS~ EWALD_NKDMX~ EWALD_RPAD~ EWALD_TOL~")

      call tkadd(" HAM~ HAM_NSPIN&")

      call tkadd(" PGF~ PGF_GFOPTS~ PGF_MODE~ PGF_PLATL PGF_PLATR PGF_SPARSE~")

      call tkadd(" SPEC SPEC_SCLWSR~ SPEC_WSRMAX~ SPEC_OMAX1~ SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM  SPEC_ATOM_CSTRMX~ SPEC_ATOM_Z SPEC_ATOM_R SPEC_ATOM_R/A SPEC_ATOM_R/W")

      call tkadd(" SITE SITE_FILE&")
      call tkadd(" SITE_ATOM SITE_ATOM_POS SITE_ATOM_DPOS~ SITE_ATOM_XPOS")

      call tkadd(" SYMGRP~")

C --- LM tags ---
      call nswadd()
      call tkadd(" LM::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))

      call tkadd(" OPTIONS_FRZ~ OPTIONS_HF~ OPTIONS_SAVVEC~ OPTIONS_SCR~ OPTIONS_SHARM~ OPTIONS_RQUAD~")

      call tkadd(trim(asa_sw))

      call tkadd(" SPEC")
      call tkadd(sscale)
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_ALPHA~")
      call tkadd(" SPEC_ATOM_C-HOLE~ SPEC_ATOM_C-HQ~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_GROUP~ SPEC_ATOM_GRP2~")
      call tkadd(" SPEC_ATOM_HCR~ SPEC_ATOM_HCR/R~")
C     call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_L~ SPEC_HFAC_D~ SPEC_HFAC_V~")
      call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_D~ SPEC_HFAC_V~")
      call tkadd(" SPEC_ATOM_HFACB~ SPEC_ATOM_HFACL~ SPEC_ATOM_HFACD~")
      call tkadd(" SPEC_ATOM_IDXDN~ SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_IDU~ SPEC_ATOM_UH~ SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_FRZC~ SPEC_ATOM_FRZD~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_MIX~")
      call tkadd(" SPEC_ATOM_Z SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_DV~")
      call tkadd(" SPEC_ATOM_R SPEC_ATOM_R/A SPEC_ATOM_R/W")

      call tkadd(" SITE SITE_FILE&")
      call tkadd(" SITE_ATOM SITE_ATOM_POS SITE_ATOM_DPOS~ SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_RELAX~")
      call tkadd(" SITE_ATOM_ROT~")

      call tkadd(" SYMGRP~")

      call tkadd(trim(bz_sw))
      call tkadd(" BZ_FSMOM~ BZ_MULL~ BZ_DMATK~")

      call tkadd(" EWALD~ EWALD_AS~ EWALD_NKDMX~ EWALD_RPAD~ EWALD_TOL~")

      call tkadd(" ITER~ ITER_AMIX~ ITER_CONV& ITER_CONVC& ITER_MIX& ITER_NIT&")
      call tkadd(" ITER_NITU~ ITER_TOLU~ ITER_UMIX~")
      call tkadd(" ITER_NRMIX~")
      call tkadd(" ITER_XIPMX~")
      call tkadd(" ITER_FIT~ ITER_FIT_MODE~ ITER_FIT_NDOS~ ITER_FIT_WG~")
      call tkadd(" ITER_FIT_EBFIT~ ITER_FIT_NBFIT~ ITER_FIT_NBFITF~")
      call tkadd(" ITER_FIT_SHFT~ ITER_FIT_WT~")
      call tkadd(" ITER_FIT_LAM~ ITER_FIT_SCL~")

      call tkadd(asa_start)
      call tkadd(" START_FREE~")

      call tkadd(" HAM~")
C     call tkadd(" HAM_NSPIN& HAM_BXCSCAL~ HAM_SOSCAL~  HAM_BFIELD~")
      call tkadd(" HAM_NSPIN& HAM_BXCSCAL~ HAM_BFIELD~")
      call tkadd(" HAM_ELIND~")
      call tkadd(" HAM_EWALD~")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FTMESH~")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")
C     call tkadd(" HAM_GMAX")
      call tkadd(" HAM_KMTO~ HAM_NMTO~")
      call tkadd(" HAM_SPINTRQ~")
      call tkadd(" HAM_PMIN~ HAM_PMAX~")
      call tkadd(" HAM_REL~ HAM_SO& HAM_NONCOL~ HAM_SS~")
      call tkadd(" HAM_QASA~")
      call tkadd(" HAM_RDVEXT~")
      call tkadd(" HAM_RDSIG~ HAM_RSRNGE~ HAM_RSSTOL~")
      call tkadd(" HAM_TOL~")
      call tkadd(" HAM_VMTZ~")
      call tkadd(" HAM_SX~ HAM_SXOPTS~")
      call tkadd(" HAM_UDIAG~")

      call tkadd(" MAP~")

      call tkadd(" OPTICS~ OPTICS_MODE~ OPTICS_MEFAC~ OPTICS_FFMT~ OPTICS_IQ~")
      call tkadd(" OPTICS_EMPBND~ OPTICS_FILBND~ OPTICS_ESMR~ OPTICS_ALLTRANS~ OPTICS_NMP~ OPTICS_W~")
      call tkadd(" OPTICS_ESCISS~ OPTICS_NPTS OPTICS_DW OPTICS_WINDOW")
      call tkadd(" OPTICS_IMREF OPTICS_FERMI~ OPTICS_KT")
      call tkadd(" OPTICS_LTET~ OPTICS_PART~")
      call tkadd(" OPTICS_CHI2~ OPTICS_CHI2_NCHI2 OPTICS_CHI2_AXES")

      call tkadd(" DYN~")
      call tkadd(" DYN_SSTAT~ DYN_SSTAT_MODE~ DYN_SSTAT_SCALE~ DYN_SSTAT_TAU DYN_SSTAT_MAXT~ DYN_SSTAT_ETOL~")

C --- LMF tags ---
      call nswadd()
      call tkadd(" LMF::" )
      call tkadd(" VERS_LM VERS_FP")
C     call tkadd(" EXPRESS~")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))

      call tkadd(" OPTIONS~ OPTIONS_FRZ~ OPTIONS_HF~ OPTIONS_RQUAD~ OPTIONS_SHORBZ~ OPTIONS_WRONSK~")

      call tkadd(" SPEC")
      call tkadd(sscale)
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~ SPEC_ATOM_Z")
      call tkadd(" SPEC_ATOM_AMASS~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_GRP2~")
      call tkadd(" SPEC_ATOM_FRZWF~")
C     call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_L~ SPEC_HFAC_V~")
      call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_V~")
      call tkadd(" SPEC_ATOM_HFACB~ SPEC_ATOM_HFACL~")
      call tkadd(" SPEC_ATOM_IDMOD~ SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_IDU~ SPEC_ATOM_UH~ SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_LMX~ SPEC_ATOM_LMXA~ SPEC_ATOM_LMXL~ SPEC_ATOM_KMXA~ SPEC_ATOM_KMXV~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_PZ~")
      call tkadd(" SPEC_ATOM_R SPEC_ATOM_R/A SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_LFOCA~ SPEC_ATOM_RFOCA~")
      call tkadd(" SPEC_ATOM_RSMA~ SPEC_ATOM_RSMG~ SPEC_ATOM_RSMFA~ SPEC_ATOM_RS3~")
      call tkadd(" SPEC_ATOM_RSMH~ SPEC_ATOM_RSMH2~ SPEC_ATOM_EH~ SPEC_ATOM_EH2~")
      call tkadd(" SPEC_ATOM_C-HOLE~ SPEC_ATOM_C-HQ~")

C ... for LMFRS
C      call tkadd(" SPEC_ATOM_HCR~ SPEC_ATOM_HCR/R~")
C      call tkadd(" STR~")
C      call tkadd(" STR_ENV_MODE~")
C      call tkadd(" STR_ENV_NEL~")
C      call tkadd(" STR_ENV_EL")
C      call tkadd(" STR_RMAX~ STR_RMAXS~")
C      call tkadd(" STR_RVL/R~ STR_VLFUN~")
C      call tkadd(" STR_NEIGHB~ STR_DELRX~ STR_TOLG~")

      call tkadd(" SITE")
      call tkadd(" SITE_FILE&")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
      call tkadd(" SITE_ATOM_POS SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_RELAX~ SITE_ATOM_RMAXS~")

      call tkadd(" SYMGRP~")

      call tkadd(trim(bz_sw))
      call tkadd(" BZ_FSMOM~ BZ_MULL~ BZ_DMATK~")

      call tkadd(" EWALD~ EWALD_AS~ EWALD_NKDMX~ EWALD_RPAD~ EWALD_TOL~")

      call tkadd(" ITER~ ITER_AMIX~ ITER_CONV& ITER_CONVC& ITER_DHFTOL& ITER_MIX& ITER_NIT&")
      call tkadd(" ITER_NITU~ ITER_TOLU~ ITER_UMIX~")
      call tkadd(" ITER_NRMIX~")

      call tkadd(" HAM")
      call tkadd(" HAM_NSPIN&")
      call tkadd(" HAM_ELIND~")
      call tkadd(" HAM_EWALD~")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FRZWF~ HAM_V0BETA~")
      call tkadd(" HAM_GMAX# HAM_FTMESH#")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")
      call tkadd(" HAM_PMAX~ HAM_PMIN~ HAM_PWMODE~ HAM_OVEPS~ HAM_OVNCUT~")
C     call tkadd(" HAM_PMAX~ HAM_PMIN~ HAM_PWMODE~ HAM_OVEPS~ HAM_OVNCUT~ HAM_ALFSI~")
      call tkadd(" HAM_AUTOBAS_PNU& HAM_AUTOBAS_LOC& HAM_AUTOBAS_MTO& HAM_AUTOBAS_PFLOAT& HAM_AUTOBAS_GW&")
C     call tkadd(" HAM_REL~ HAM_SO& HAM_NONCOL~ HAM_SAXIS~ HAM_BFIELD~ HAM_BXCSCAL~ HAM_SOSCAL~ HAM_BXC0~")
      call tkadd(" HAM_REL~ HAM_SO& HAM_NONCOL~ HAM_SAXIS~ HAM_BFIELD~ HAM_BXCSCAL~ HAM_BXC0~")
      call tkadd(" HAM_RDSIG~ HAM_RSRNGE~ HAM_RSSTOL~ HAM_SIGP~")
      call tkadd(" HAM_TOL~")
C     call tkadd(" HAM_BLOCKD~")
C     call tkadd(" HAM_SHFAC~")
      call tkadd(" OPTICS~ OPTICS_MODE~ OPTICS_MEFAC~ OPTICS_FFMT~ OPTICS_IQ~")
      call tkadd(" OPTICS_EMPBND~ OPTICS_FILBND~ OPTICS_ESMR~ OPTICS_ALLTRANS~ OPTICS_NMP~ OPTICS_W~")
      call tkadd(" OPTICS_ESCISS~ OPTICS_NPTS OPTICS_DW OPTICS_WINDOW")
      call tkadd(" OPTICS_IMREF OPTICS_FERMI~ OPTICS_KT")
      call tkadd(" OPTICS_LTET~ OPTICS_PART~")

      call tkadd(" DYN~")
      call tkadd(dyn_mstat)
      call tkadd(dyn_md)
      call tkadd(" DYN_NIT~")

      call tkadd(" GW~ GW_NKABC~ GW_BZJOB~")

C     call tkadd(" DMFT~ ")

C --- LMFA tags (relative to lmf) ---
      call nswadd()
      call tkadd(" LMFA::" )

      call tkadd(" VERS_FP! BZ! OPTICS! DYN! GW! SYMGRP! EWALD")

      call tkadd(" OPTIONS_SHORBZ! OPTIONS_NMCORE~ OPTIONS_WRONSK~")

      call tkadd(" SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_RCFA~ SPEC_ATOM_GRP2! SPEC_ATOM_FRZWF! SPEC_ATOM_KMXA! SPEC_ATOM_KMXV!")

      call tkadd(" SITE_ATOM_RMAXS!")

C     call tkadd(" HAM_SO! HAM_NONCOL! HAM_SAXIS! HAM_SOSCAL! HAM_GMAX! HAM_FTMESH! HAM_TOL!")
      call tkadd(" HAM_SO! HAM_NONCOL! HAM_SAXIS! HAM_GMAX! HAM_FTMESH! HAM_TOL~")
C      call tkadd(" HAM_FRZWF! HAM_FORCES! HAM_ELIND! HAM_PWMODE! HAM_OVEPS! HAM_OVNCUT! HAM_ALFSI!")
      call tkadd(" HAM_FRZWF! HAM_FORCES! HAM_ELIND! HAM_PWMODE! HAM_OVEPS! HAM_OVNCUT!")
      call tkadd(" HAM_RDSIG! HAM_SIGP! HAM_EWALD!")
      call tkadd(" HAM_AUTOBAS_QLOC& HAM_AUTOBAS_ELOC& HAM_AUTOBAS_LMXB& HAM_AUTOBAS_LMTO&")
      call tkadd(" HAM_AUTOBAS_RSMMX& HAM_AUTOBAS_EHMX&")
!     For Jerome's new basis construction scheme
      call tkadd(" HAM_AUTOBAS_EIN& HAM_AUTOBAS_EOUT& HAM_AUTOBAS_EH&")

      call tkadd(" ITER_AMIX! ITER_CONV! ITER_CONVC! ITER_MIX! ITER_NIT!")
      call tkadd(" ITER_NITU! ITER_TOLU! ITER_UMIX!")

C --- LMFGWD tags (relative to lmf) ---
      call nswadd()
      call tkadd(" LMFGWD::" )

      call tkadd(" BZ_NKABC~") ! Possibly remove the whole cetegory?

      call tkadd(" OPTICS! ITER! DYN!")
      call tkadd(" OPTIONS_HF! OPTIONS_NMCORE~")
      call tkadd(" SPEC_ATOM_GRP2!")
      call tkadd(" SPEC_ATOM_LMXPB~ SPEC_ATOM_PBAS~ SPEC_ATOM_PBAS2~")

      call tkadd(" GW~")
      call tkadd(" GW_CODE~")
      call tkadd(" GW_DELRE~ GW_DELTA~ GW_DELTAW~")
      call tkadd(" GW_ECUTS~ GW_GCUTB~ GW_GCUTX~ GW_GSMEAR~ GW_ECUTPB~")
      call tkadd(" GW_MKSIG~")
      call tkadd(" GW_NBAND~")
      call tkadd(" GW_NIME~")
      call tkadd(" GW_NKABC~ GW_BZJOB~")
      call tkadd(" GW_PBTOL~")
      call tkadd(" GW_QOFFP~")
      call tkadd(" GW_USEBSEW~")

      call tkadd(" SITE_ATOM_RELAX! SITE_ATOM_RMAXS!")

C --- LMFGW tags (relative to lmf) ---
      call nswadd()
      call tkadd(" LMFGW::" )

      call tkadd(" OPTICS! ITER! DYN!")
      call tkadd(" OPTIONS_HF! OPTIONS_NMCORE~")
      call tkadd(" SPEC_ATOM_GRP2!")
      call tkadd(" SPEC_ATOM_LMXPB~ SPEC_ATOM_PBAS~ SPEC_ATOM_PBAS2~")

      call tkadd(" GW~")
      call tkadd(" GW_CODE~")
      call tkadd(" GW_DELRE~ GW_DELTA~ GW_DELTAW~")
      call tkadd(" GW_ECUTS~ GW_GCUTB~ GW_GCUTX~ GW_GSMEAR~ GW_ECUTPB~")
      call tkadd(" GW_MKSIG~")
      call tkadd(" GW_NBAND~")
      call tkadd(" GW_NIME~")
      call tkadd(" GW_NKABC~ GW_BZJOB~")
      call tkadd(" GW_PBTOL~")
      call tkadd(" GW_QOFFP~")
      call tkadd(" GW_USEBSEW~")

      call tkadd(" SITE_ATOM_RELAX! SITE_ATOM_RMAXS!")

C --- LMFGWS tags (relative to lmf) ---
      call nswadd()
      call tkadd(" LMFGWS::" )

      call tkadd(" BZ_NKABC~") ! Possibly remove the whole cetegory?

      call tkadd(" OPTICS! ITER! DYN! OPTIONS_HF!")
      call tkadd(" SPEC_ATOM_GRP2!")
      call tkadd(" SITE_ATOM_RELAX! SITE_ATOM_RMAXS!")

C --- LMFDMFT tags (relative to lmf) ---
      call nswadd()
      call tkadd(" LMFDMFT::" )
C     call tkadd(" DMFT DMFT_NKABC~ DMFT_PROJ~")
      call tkadd(" BZ_ZINSUL~")
      call tkadd(" DMFT DMFT_NKABC~ DMFT_NLOHI~ DMFT_WLOHI~ DMFT_PROJ~")
      call tkadd(" DMFT_BETAR DMFT_BETAK DMFT_BETA DMFT_NOMEGA DMFT_KNORM~ DMFT_BROAD~")
      call tkadd(" DMFT_NOMF~ DMFT_NOMB~")

      call tkadd(" DMFT_BLOCK~ DMFT_BLOCK_L DMFT_BLOCK_QSPLIT~")
      call tkadd(" DMFT_BLOCK_SITES DMFT_BLOCK_SIDXD~ DMFT_BLOCK_SIDXM~ DMFT_BLOCK_UMODE~")

C     To be removed
      call tkadd(" SPEC_ATOM_IDD~")
C     call tkadd(" SPEC_ATOM_IDD~ SPEC_ATOM_IDM~")

C --- LMMC tags (relative to base) ---
      call nswadd()
      call tkadd(" LMMC::" )
      call tkadd(" VERS_MOL")

      call tkadd(" STRUC_NBASP~")

      call tkadd(" OPTIONS OPTIONS_XCQS OPTIONS_RQUAD~")

      call tkadd(" SYMGRP! PGF! EWALD!")

      call tkadd(" SPEC_ATOM_CSTRMX!")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_IDU~ SPEC_ATOM_UH~ SPEC_ATOM_JH~")

      call tkadd(" SITE_ATOM_XPOS!")

      call tkadd(" SPEC_ATOM_AMASS~")
      call tkadd(" SPEC_ATOM_COLOUR~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_LMX~ SPEC_ATOM_LMXA~ SPEC_ATOM_LMXB SPEC_ATOM_LMXL~")
      call tkadd(" SPEC_ATOM_LXI SPEC_ATOM_EXI")
      call tkadd(" SPEC_ATOM_RADIUS~ SPEC_ATOM_RCUT SPEC_ATOM_RHAM SPEC_ATOM_RINT SPEC_ATOM_RSMG")

      call tkadd(" SITE_ATOM_RELAX~")
      call tkadd(" SITE_ATOM_V0~")
      call tkadd(" SITE_PM~ SITE_PM_POS SITE_PM_Q SITE_PM_P~")

      call tkadd(" BZ~ BZ_W~")

      call tkadd(" TCF TCF_ADEC TCF_NALF TCF_NBISI TCF_NCUPL~ TCF_NDUST~ TCF_WZTCF~")

      call tkadd(" ITER ITER_CONV& ITER_CONVC& ITER_MIX& ITER_NIT# ITER_NRMIX~")

      call tkadd(" HAM_GMAX#")
      call tkadd(" HAM_ALFSI")
      call tkadd(" HAM_DABC~")
      call tkadd(" HAM_DQVAL~")
      call tkadd(" HAM_EBAS")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FRZWF~")
C     call tkadd(" HAM_FTMESH")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")

      call tkadd(" DYN~")
      call tkadd(dyn_mstat)
      call tkadd(dyn_md)
      call tkadd(" DYN_NIT~")

C --- LMGF tags ---
      call nswadd()
      call tkadd(" LMGF::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))

      call tkadd(" OPTIONS~ OPTIONS_FRZ~ OPTIONS_HF~ OPTIONS_SAVVEC~ OPTIONS_SCR~ OPTIONS_SHARM~ OPTIONS_RQUAD~")
      call tkadd(trim(asa_sw))

      call tkadd(" SPEC")
      call tkadd(sscale)
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_ALPHA~")
      call tkadd(" SPEC_ATOM_C-HOLE~ SPEC_ATOM_C-HQ~")
      call tkadd(" SPEC_ATOM_C~ SPEC_ATOM_CPA~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_DV~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_GROUP~ SPEC_ATOM_GRP2~")
      call tkadd(" SPEC_ATOM_HCR~ SPEC_ATOM_HCR/R~")
C     call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_L~ SPEC_HFAC_D~")
      call tkadd(" SPEC_HFAC~ SPEC_HFAC_B~ SPEC_HFAC_D~")
      call tkadd(" SPEC_ATOM_HFACB~ SPEC_ATOM_HFACL~ SPEC_ATOM_HFACD~ SPEC_ATOM_BEFF~")
      call tkadd(" SPEC_ATOM_IDMOD~ SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_IDU~ SPEC_ATOM_UH~ SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_MIX~")
      call tkadd(" SPEC_ATOM_NTHET~ SPEC_ATOM_NANG~")
      call tkadd(" SPEC_ATOM_Z SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_R SPEC_ATOM_R/A SPEC_ATOM_R/W")

      call tkadd(" SITE")
      call tkadd(" SITE_FILE&")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~ SITE_ATOM_POS SITE_ATOM_XPOS SITE_ATOM_RELAX~ SITE_ATOM_ROT~ SITE_ATOM_SID~ SITE_ATOM_VSHFT~")

      call tkadd(" SYMGRP~")

      call tkadd(trim(bz_sw))
      call tkadd(" BZ_EMESH")

      call tkadd(" EWALD~ EWALD_AS~ EWALD_NKDMX~ EWALD_RPAD~ EWALD_TOL~")

      call tkadd(" ITER~ ITER_AMIX~ ITER_CONV& ITER_CONVC& ITER_MIX& ITER_XIPMX~ ITER_NIT& ITER_NRMIX~")
      call tkadd(" ITER_NITU~ ITER_TOLU~ ITER_UMIX~")

      call tkadd(asa_start)

      call tkadd(" HAM~")
      call tkadd(" HAM_UDIAG~")
C     call tkadd(" HAM_NSPIN& HAM_BXCSCAL~ HAM_SOSCAL~ HAM_BFIELD~")
      call tkadd(" HAM_NSPIN& HAM_BXCSCAL~ HAM_BFIELD~")
      call tkadd(" HAM_ELIND~")
      call tkadd(" HAM_EWALD~")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FTMESH~")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")
      call tkadd(" HAM_SPINTRQ~")
      call tkadd(" HAM_PMIN~ HAM_PMAX~")
      call tkadd(" HAM_REL~ HAM_SO& HAM_NONCOL~ HAM_SS~")
      call tkadd(" HAM_QASA~")
      call tkadd(" HAM_RDSIG~ HAM_RSRNGE~ HAM_RSSTOL~")
      call tkadd(" HAM_TOL~")
      call tkadd(" HAM_VMTZ~")
      call tkadd(" HAM_SX~ HAM_SXOPTS~")
      call tkadd(" HAM_RDVEXT~")

      call tkadd(" GF~ GF_GFOPTS~ GF_MODE GF_DLM~ GF_BXY~ GF_TEMP~")

      call tkadd(" DYN~ DYN_SSTAT~ DYN_SSTAT_MODE~ DYN_SSTAT_SCALE~ DYN_SSTAT_TAU DYN_SSTAT_MAXT~ DYN_SSTAT_ETOL~")

C --- LMPG tags (relative to lmgf) ---
      call nswadd()
      call tkadd(" LMPG::" )

      call tkadd(" GF!")

      call tkadd(" PGF PGF_MODE PGF_GFOPTS~ PGF_PLATL PGF_PLATR PGF_SPARSE~")

      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_PLV~")

      call tkadd(" HAM_SX! HAM_SXOPTS! HAM_RDSIG! HAM_RSRNGE! HAM_RSSTOL!")

C --- LMDOS tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMDOS::" )

      call tkadd(" PGF!")

      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_Z SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~")

      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_LMXA~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_ALPHA~") ! Is this needed?
      call tkadd(" SPEC_ATOM_DV~")  ! Is this needed?
      call tkadd(" SPEC_ATOM_MIX~") ! Is this needed?
      call tkadd(" SPEC_ATOM_EREF~") ! Is this needed?
      call tkadd(" SPEC_ATOM_GROUP~") ! Is this needed?
      call tkadd(" SPEC_ATOM_GRP2~") ! Is this needed?

      call tkadd(" SITE_ATOM_PL~")

      call tkadd(trim(bz_sw))
      call tkadd(" BZ_COND_MODE~ BZ_COND_DIR")
      call tkadd(" BZ_MULL~") ! needed?
      call tkadd(" BZ_FSMOM~") ! needed?

      call tkadd(" HAM_SO&")

C --- LMCHK tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMCHK:: VERS_LM")

      call tkadd(" BANDS~ BANDS_MODE~ BANDS_SYML~ BANDS_MESH~ BANDS_QP~ BANDS_FN~")

      call tkadd(" HAM_SO& HAM_NONCOL~ HAM_SS~")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")

      call tkadd(" OPTIONS~ OPTIONS_NESABC~ OPTIONS_RMAXES~ OPTIONS_RMINES~")

      call tkadd(" SPEC_ATOM_HCR~ SPEC_ATOM_HCR/R~")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~ SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_DV~")

      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_ROT~")

      call tkadd(" STR~ STR_MXNBR~ STR_RMAX~ STR_RMAXS~ STR_EQUIV~ STR_SHOW~")

      call tkadd(" BZ BZ_SYML~")

C --- LMXBS tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMXBS::" )

      call tkadd(" HAM! OPTIONS!")
      call tkadd(" SPEC_ATOM_RHAM~ SPEC_ATOM_Z~ SPEC_ATOM_LMX~")
      call tkadd(" SITE_ATOM_PL~")
      call tkadd(" STR~ STR_RMAX~ STR_RMAXS~")

C --- LMPLAN tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMPLAN::" )
!     call tkadd(" SPEC_ATOM_RHAM~")

      call tkadd(asa_start)

!     call tkadd(" HAM_REL~ HAM_SO& HAM_NONCOL~ HAM_SS~")
      call tkadd(" HAM_SO& HAM_NONCOL~ HAM_SS~")

      call tkadd(" BZ BZ_BZJOB~ BZ_COND~ BZ_DELEF~ BZ_DOS~ BZ_EF0~ BZ_FSMOM~ BZ_GETQP~ BZ_ZVAL~")
      call tkadd(" BZ_METAL& BZ_N~ BZ_NEVMX~ BZ_NKABC# BZ_NOINV~ BZ_NPTS~ BZ_TETRA~ BZ_W~ BZ_ZBAK~")
      call tkadd(" BZ_EMESH~") ! needed to pick up vne

      call tkadd(" PLANE~ PLANE_NORMAL~ PLANE_X~ PLANE_Y~")

      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~ SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_DV~")

      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_ROT~")

C --- LMCTL tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMCTL::" )
      call tkadd(" VERS_ASA")

      call tkadd(" HAM_SO& HAM_NONCOL~ HAM_SS~")

      call tkadd(" OPTIONS! BANDS! EWALD! STR! HAM_XCFUN!")

      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~ SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_DV~")
      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_ROT~")

      call tkadd(" SPEC_SCLWSR! SPEC_OMAX1! SPEC_OMAX2! SPEC_WSRMAX! SPEC_ATOM_CSTRMX!")

C --- LMSCELL tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMSCELL::" )

      call tkadd(" HAM_SO& HAM_NONCOL~")

      call tkadd(" SYMGRP! PGF! EWALD!")

      call tkadd(" DYN~ DYN_SDYN~ DYN_SDYN_KT DYN_SDYN_TS DYN_SDYN_TTOT DYN_SDYN_TEQU~")

      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_NR~ SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_P~ SPEC_ATOM_Q~ SPEC_ATOM_MMOM~ SPEC_ATOM_DV~")
      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_ROT~")

      call tkadd(" SPEC_SCLWSR! SPEC_OMAX1! SPEC_OMAX2! SPEC_WSRMAX!")
      call tkadd(" SPEC_ATOM_IDMOD~ SPEC_ATOM_DV! SPEC_ATOM_CSTRMX! SPEC_ATOM_HCR! SPEC_ATOM_HCR/R!")

      call tkadd(" PGF~ PGF_PLATL PGF_PLATR")

      call tkadd(" STRUC_SLAT~")

      call tkadd(" SITE_ATOM_RELAX~")

C     call tkadd(" EWALD~ EWALD_NKDMX~")

C --- LMSTR tags (relative to BASE) ---
      call nswadd()
      call tkadd(" LMSTR::" )
      call tkadd(" VERS_ASA")

      call tkadd(" HAM_KMTO~ HAM_NMTO~ HAM_EWALD~")
      call tkadd(" HAM_GGA~ HAM_XCFUN&")

      call tkadd(" SPEC_ATOM_Z~")
      call tkadd(" SPEC_ATOM_ALPHA~")
      call tkadd(" SPEC_ATOM_HCR~ SPEC_ATOM_HCR/R~ SPEC_ATOM_RSMINL~")
      call tkadd(" SPEC_ATOM_HSFITK~ SPEC_HSFITK~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_LMX~ SPEC_ATOM_KMXA~ SPEC_ATOM_RSMA~")
      call tkadd(" SPEC_ATOM_RSMH~")
      call tkadd(" SPEC_ATOM_EHVL~")

      call tkadd(" SITE_ATOM_PL~ SITE_ATOM_PLV~")

      call tkadd(" STR~")
      call tkadd(" STR_RMAX~ STR_RMAXS~")
      call tkadd(" STR_MXNBR~")
      call tkadd(" STR_IINV~ STR_IINV_NCUT~ STR_IINV_NIT~ STR_IINV_TOL~")
      call tkadd(" STR_LMAXW~ STR_DELRW~")
      call tkadd(" STR_ENV_MODE~ STR_ENV_NEL~ STR_ENV_EL~ STR_ENV_RSMINL~")
      call tkadd(" STR_RVL/R~ STR_VLFUN~")
      call tkadd(" STR_EQUIV~")
      call tkadd(" STR_SHOW~")

C --- TBE tags (relative to BASE) ---
      call nswadd()
      call tkadd(" TBE::" )
      call tkadd(" VERS_TB")

      call tkadd(" PGF!")

      call tkadd(" OPTIONS~ OPTIONS_SHARM~ OPTIONS_Q~")

      call tkadd(" HAM_SO&")

      call tkadd(" SPEC_ATOM_Z~ SPEC_ATOM_R!")
      call tkadd(" SPEC_ATOM_AMASS~ SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_ISOTOPE")
      call tkadd(" SPEC_ATOM_IDU~ SPEC_ATOM_I~ SPEC_ATOM_UH~ SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_FRZQ1~ SPEC_ATOM_FRZVSO~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_LMX~ SPEC_ATOM_LMXA~ SPEC_ATOM_LMXL~")
      call tkadd(" SPEC_ATOM_QPOL~")
      call tkadd(" SPEC_ATOM_RHAM~")
      call tkadd(" SPEC_ATOM_VSO~")
      call tkadd(" SPEC_ATOM_COLOUR~ SPEC_ATOM_RADIUS~")

      call tkadd(" SITE_ATOM_DELTA~ SITE_ATOM_RELAX~")

      call tkadd(" STR~ STR_MXNBR~ STR_SHOW~")

      call tkadd(trim(bz_sw))
      call tkadd(" BZ_FSMOM~")
      call tkadd(" BZ_MULL~")

      call tkadd(" ITER~")
      call tkadd(" ITER_MIX& ITER_XIPMX~")
      call tkadd(" ITER_CONV& ITER_CONVC& ITER_NIT&")
      call tkadd(" ITER_FIT~ ITER_FIT_MODE~ ITER_FIT_NBFIT ITER_FIT_EBFIT ITER_FIT_RFIT ITER_FIT_LAM~ ITER_FIT_SCL~")
C     call tkadd(" ITER_FIT_SHFT~ ITER_FIT_WT~")

      call tkadd(" TB")
      call tkadd(" TB_3PV~")
      call tkadd(" TB_ADDES~")
      call tkadd(" TB_CRYSF~")
      call tkadd(" TB_EVDISC~")
      call tkadd(" TB_FIJ~")
      call tkadd(" TB_FORCES~")
      call tkadd(" TB_GAMMA~")
      call tkadd(" TB_IODEL~")
      call tkadd(" TB_MOL~")
      call tkadd(" TB_NOUAVG~")
      call tkadd(" TB_OVCF~ TB_OVLP~")
      call tkadd(" TB_PAIR~")
      call tkadd(" TB_RHO~")
      call tkadd(" TB_RMAXH~")
      call tkadd(" TB_TBU~")
      call tkadd(" TB_TRH~")
      call tkadd(" TB_SCALE~")
      call tkadd(" TB_UL~")

      call tkadd(" DYN~")
      call tkadd(" DYN_KT~")
      call tkadd(" DYN_MODET~")
      call tkadd(dyn_mstat)
      call tkadd(dyn_md)
      call tkadd(" DYN_NIT~")
      call tkadd(" DYN_TEQU~")
      call tkadd(" DYN_TS~")
      call tkadd(" DYN_TTOT~")

      call tkadd(" START START_CNTROL~ START_ATOM~ START_ATOM_P~ START_ATOM_Q~")


C --- LMMAG tags (relative to BASE) ---
      call nswadd()
      call tkadd(" MMAG::" )
      call tkadd(" VERS_MM")

      call tkadd(" OPTIONS! SYMGRP! PGF!")
      call tkadd(" HAM_NSPIN! HAM_SO& HAM_SS~ HAM_BFIELD~")

      call tkadd(" DYN")
      call tkadd(" DYN_BSINT~ DYN_BSINT_MX~ DYN_BSINT_MI~ DYN_BSINT_NEW~ DYN_BSINT_NSEQ~ DYN_BSINT_TOL DYN_BSINT_TS0~")
      call tkadd(" DYN_GD~ DYN_GD_CT~ DYN_GD_MODET~ DYN_GD_NTHERM~")
      call tkadd(" DYN_SDYN DYN_SDYN_KT DYN_SDYN_TS DYN_SDYN_TTOT DYN_SDYN_TEQU~")
      call tkadd(" DYN_MD~ DYN_MMHAM")

      call tkadd(" SPEC_SCLWSR! SPEC_OMAX1! SPEC_OMAX2! SPEC_WSRMAX! SPEC_ATOM_CSTRMX!")
      call tkadd(" SPEC_ATOM_RHAM SPEC_ATOM_Z~ SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_R!")

      call tkadd(" SITE_ATOM_RELAX~ SITE_ATOM_ROT~")

      end subroutine toksw_init

C      integer function dmft_equiv_cix(ncix,dmft_ib,dmft_l,ib,l)
CC- Check to see if (ib,l) pair occurs in any of ncix blocks
C      implicit none
CC ... Passed parameters
C      integer ncix,ib,l
C      integer :: dmft_ib(ncix)    ! (1:ncix) Equivalent to s_dmft_ib
C      integer :: dmft_l(ncix)     ! (1:ncix) Not equivalent to s_dmft_l
CC ... Local parameters
C      integer i,j
C      procedure(integer) :: iinear
C
C      i = 0; dmft_equiv_cix = 0
C      do  while (i<ncix)
C        j = iinear(ncix-i,dmft_ib(1+i),ib,1)
C        i = i+j
C        if (dmft_ib(i) /= ib) cycle
C        if (dmft_l(i) /= l) cycle
C        dmft_equiv_cix = i
C        return
C      enddo
C
C      end
