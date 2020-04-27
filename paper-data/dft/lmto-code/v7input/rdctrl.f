      subroutine rdctrl(recrd,recln,nrecs,prgnam,vrsion,ivn,passn,
     .  slabl_,s_bz,s_ctrl,s_ham,s_pot,s_lat,s_mix,s_spec,s_site,
     .  nbmx,s_str,s_move,s_tb,s_optic,s_gw,s_dmft,s_strn)
C- Main input for LMTO programs
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  fsmom
Co     Stored:     def dosw ef efmax fsmom egap lio lcond lmet lmull
Co                 lopt lshft n ndos nevmx nkabc nkp range semsh w zval
Co     Allocated:  *
Cio    Elts passed:lmet lio lmull
Cio    Passed to:  *
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nspec nspin nl nclasp dclabl ics lasa
Ci                 lncol lham lgen3 lrsa lcd lrel lmet lxcf lpgf sdmod
Ci                 mdprm ltb lstonr nsite nbasp lordn pgfsl pgfvl
Co     Stored:     lasa defm elin lbas lcd nccomp maxmem lcgf ldos lfp
Co                 lfrce lgen3 ldlm lbxc tdlm lham lmet lncol loptc
Co                 lordn lrsa lpgf lqp lrel lrs lscr lrquad lstonr lsx
Co                 ltb lves lxcf maxit mdprm modep nbas nbasp nesabc
Co                 nitmv nl nmap nsite nspec nspin nvario omax1 omax2
Co                 quit plbnd rmaxes rmines sclwsr sdmod sdprm sdxsi
Co                 smalit tol wsrmax zbak pfloat pos ncl pgfsl pgfvl
Co                 ips npl pgplp npadl npadr
Co     Allocated:  *
Cio    Elts passed: initc ics lcd nl lham lmet lbas lqp ldos maxit
Cio                lncol lfrce lxcf lrs loptc lrel tol lasa lcgf lpgf
Cio                lsx lves lgen3 sdprm sdxsi ips clssl cllst clp pgfsl
Cio                pgfvl pgplp
Cio    Passed to:  susite
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     lxcf alfsi dabc elind nmto kmto nqsig ndhrs lncol
Co                 lham lgen3 lrsa nkaph pmax pmin lsig sigp pnudef
Co                 rdvext qss rsrnge rsstol udiag pwmode npwpad pwemin
Co                 pwemax basopt oveps ovncut nlibu lmaxu seref neula
Co                 nbf
Co     Allocated:  eula magf
Cio    Elts passed:eula magf
Cio    Passed to:  susite
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     vmtz0 vconst
Co     Allocated:  vshft aamom
Cio    Elts passed:pnu qnu pp ves vshft
Cio    Passed to:  susite
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  lsym pos alat plat platl platr plate
Co     Stored:     alat as avw nkdmx nkqmx gam gmax nabc ldist dist
Co                 plat plat2 rpad slat tol tolft vol platl platr lsym
Co                 pos
Co     Allocated:  pos igv2
Cio    Elts passed:lsym avw pos
Cio    Passed to:  ioxsf susite clsset
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  *
Co     Stored:     b bl bv elind fn kill lxpot mmix mode model nitu
Co                 nmix nsave tolu umix w wc
Co     Allocated:  *
Cio    Elts passed:lxpot
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa z name lmxb idxdn orbp p pz norp ntorb rhoc eref
Co     Stored:     eh3 etf lmxf norp vmtz pb1 pb2 coreh name a nr alpha
Co                 coreq lxi nxi exi group grp2 hcr hsfitp idmod
Co                 ehvl idxdn idu idm jh uh kmxt kmxv lfoca rsmv lmxa
Co                 lmxb lmxl lmxpb nthet nang ncpa nbeff beff iscpa
Co                 xcpa mass mxcst orbp p pz q z colxbs radxbs rcfa
Co                 rcut rfoca rg rham rmt rs3 rsma rsmfa iq1 iq2 stni
Co                 vso qpol dv eref rint shfac ntorb rhoc
Co     Allocated:  *
Cio    Elts passed:pz name idu
Cio    Passed to:  ioorbp bcast_strx uspecb suidx suldau ioxsf ioxsf2
Cio                susite pgfset
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec clabel pos pl plv bfield eula
Co     Stored:     spec class dpole mpole clabel pl plv pos vel vshft
Co                 relax eula rmaxs ndelta delta ncomp bfield
Co     Allocated:  *
Cio    Elts passed:pl plv
Cio    Passed to:  iopos suldau ioxsf ioxsf2 susite siteperm pvsub2
C
Ci   nbmx
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  rmax
Co     Stored:     adec amode drwats iinv lmem lequiv loka lmaxw lshow
Co                 mxnbr nalf nbisi ncupl ndust wztcf nkaps rmax rfit
Co                 kaps rmaxg ivl delrx mnnbr tolg
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     gyro nmodt modt ct prmint kt ts tsequ tstot
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_tb
Ci     Elts read:  *
Co     Stored:     alam alsc shft wg wt fmode ebfit ndos nbfit nbfitf
Co                 rmfit
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  *
Co     Stored:     dw window esmr ocrng unrng esciss ltet iq loptic
Co                 lpart mefac nchi2 axes
Co     Allocated:  *
Cio    Elts passed:cls
Cio    Passed to:  *
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  *
Co     Stored:     mksig gcutb gcutx qoffp nband gsmear code lgw nkabc
Co                 lshft delre ecuts nime deltax deltaw pbtol pb1 pb2
Co                 nkabcd nlohi wlohi iproj
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  str_pack
Ci Inputs
Ci   recrd (recln*nrecs) : preprocessed input
Ci   prgnam:name of main program
Ci   vrsion:string specifying expected program version
Ci   ivn describes version:
Ci     ivn(:,1) corresponds to input file system
Ci     ivn(:,2) may have a unique program version, but usually identical to ivn(:,1)
Ci     ivn(1,:) main version number.
Ci              Input file across different versions may have some incompatibilities.
Ci     ivn(2,:) minor version number.
Ci     ivn(3,1) revision within minor version.
Ci   passn :specifies what to read from input file
Ci         :1, read main contents of input file
Ci         :2, reads class-specific info, e.g. ASA moments
Co Outputs
Co   Input file is read and data is packed into these structures:
Co   slabl_:vector of species labels
Cg Global variables
Cg   The following global variables are set by rdctrl and may be accessed by
Cg   any routine via function call 'dglob' (for double) or 'nglob' (for int)
Cg   avw   :global length scale, usu. the average Wigner-Seitz radius,
Cg         :used in various places to set a length scale for a range,
Cg         :sometimes in generating structure constants, etc.
Cg   lrel  :specifies type of Schrodinger equation
Cg         :0 nonrelativistic Schrodinger equation
Cg         :1 scalar relativistic Schrodinger equation
Cg         :2 Dirac equation
Cg   lxcf  :specifies type of XC potential.  1s digit specifies local XC:
Cg         :1 for Ceperly-Alder
Cg         :2 for Barth-Hedin (ASW fit)
Cg         :3 for PW91
Cg         :4 for PBE
Cg         :100s digit specifies type of gradient correction
Cg         :0 no gradient correction
Cg         :1 Langreth-Mehl
Cg         :2 PW91
Cg         :3 PBE
Cg         :4 PBE with Becke exchange
Cg   mxorb :nkaph * (maximum number of lm channels in any sphere)
Cg         :Used for dimensioning the indexing arrays involved in
Cg         :assembling the hamiltonian;
Cg   nbas  :number of atoms in the basis
Cg   nbasp :number of atoms in the padded basis
Cg         :(when extensions are needed, e.g. in layer GF code)
Cg   nkape :NOT USED The maximum number of envelope functions centered at
Cg         :particular R and l channel
Cg         :NB: nkape is not used now.
Cg   nkaph :The maximum number of radial functions centered at
Cg         :particular R and l channel used in the lmto basis.
Cg   nl    :1+Maximum l-cutoff for augmentation
Cg   npl   :(not set by rdctrl) number of principal layers (layer geometries)
Cg   nkaph :The maximum number of "principal quantum" numbers centered
Cg         :at a particular R and l channel --- energies for one Rl
Cg         :at which augmentation (phi-phidot) functions are made.
Cg   nsp   :1 if not spin-polarized; otherwise 2
Cg   nspec :number of species
Cg   stde  :standard error file
Cg   stdl  :standard log file
Cg   stdo  :standard output file
Cr Remarks
Cr rdctrl does:
Cr  1. read input data specified by tokens
Cr  2. If passn is 2, read class parameters from START
Cu Updates
Cu   26 Mar 19 Start on reading complete DMFT input from ctrl file
Cu   06 Apr 18 (Jerome Jackson) Added GW variable GW_USEBSEW (lmfgwd)
Cu   19 Mar 18 New --revec switch
Cu   20 Jul 17 Allow --rs~switches~...
Cu   05 Aug 15 Header information modified in prepration for v7.12
Cu   18 Nov 14 Use 10s digit of 2nd argument --rs=#1,#2,... for conditional write of rs file
Cu   19 Aug 13 Parse sections of symg; set corresponding bits in structures
Cu   09 Aug 13 New s_spec%nang, 16s bit of s_ctrl%lmet, site-resolved SO
Cu   17 Jun 13 New s_spec%shfac
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   28 Jan 13 New switch OPTIONS_NMCORE
Cu   25 Apr 12 (Belashchenko) New inputs for CPA (lmgf)
Cu   01 Sep 11 Begin migration to f90 structures
Cu   01 Apr 10 Bfield can have just z component
Cu   27 Mar 10 New parameters for Levenberg-Marquardt fitting
Cu   15 Jan 10 (Walter Lambrecht) writes crystal data in xsf format
Cu   09 Sep 09 New OTPICS inputs
Cu   11 Aug 09 Inputs for Levenberg-Marquardt fitting of ASA bands
Cu   10 May 09 New input for automatic parameter generation
Cu   28 Apr 09 Make --rdbasp work
Cu   24 Apr 09 New str->delrx,mnnbr,tolg and site->rmaxs
Cu   01 May 08 Added DLM vars ldlm,tdlm; arrays opnud,oqnud,oppd (Kirill)
Cu   19 Sep 07 (TK+MvS) Adapted from rdctrl, 1st cut at new input
C ----------------------------------------------------------------------
      use m_rdctrl
      use m_gtv
      use structures
      implicit none
C     include "mpif.h"
C ... Passed parameters
      integer(4):: recln,nrecs
      character slabl_(1)*8
      character*(120) recrd
      integer passn,nbmx
      character(len=*) :: prgnam
c      character toksw(0:30)*(*), vrsion*6
C     character  vrsion*6
      character(6):: vrsion(2)
      integer ivn(3,2)
C ... For structures
!      include 'structures.h'
      type(str_bz):: s_bz
      type(str_ctrl):: s_ctrl
      type(str_lat):: s_lat
      type(str_ham):: s_ham
      type(str_pot):: s_pot
      type(str_mix):: s_mix
      type(str_move):: s_move
      type(str_str):: s_str
      type(str_tb):: s_tb
      type(str_optic):: s_optic
      type(str_gw):: s_gw
      type(str_dmft):: s_dmft
      type(str_strn) :: s_strn(*)
      type(str_spec):: s_spec(*)
      type(str_site):: s_site(*)
C ... Dynamically allocated local arrays
      real(8),allocatable:: pp(:,:,:,:),ves(:),zc(:),pnuc(:,:,:),qnuc(:,:,:,:)
      integer, pointer :: wk(:)
      character(len=8),allocatable::clabl(:)
C ... Local parameters
      integer(4) :: dvec1(3)=1, dvec2(3)=0
      integer, parameter :: lBextz=256,lSOsite=512
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      real(8), parameter :: ry2eV = 13.60569193d0
      character fileid*64,outs*256,dc*1
      character strn*(recln)
      integer procid,nproc,master
      logical ltmp,asa,mlog,lwarn
      double precision xx(n0*2)
      integer i,ifi,irs(5),is,ix(n0*nkap0),j,j1,j2,k,k1,k2,l,
     .  lasa,lbas,lcd,lfrzw,lham,lmet,lncol,lordn,lqp,lrhop,lrs,lstsym,lsx1,ltb,
     .  n,nat,nclasp,nlibu,nlmax,nphimx,nspc,scrwid,stde,stdl,stdo

C ... basis
      double precision orbp(n0,2,nkap0)
      procedure(logical) :: cmdopt,bittst,ioorbp,cmdstr,isanrg
      procedure(integer) :: a2vec,bitand,fopna,fxst,getdig,iprint,iprt,lgunit,nglob,
     .  wordsw,cmdoptswx,isw,mpipid,nglobxst,str_pack
      procedure(real(8)) :: dglob

C --- Start of execution ---
      procid = mpipid(1)
      nproc  = mpipid(0)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      scrwid = 100
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = stdo

C --- Initialize gtv; copy recrd to rcd ---
      call gtv_setst(stdo,stdl,stde)
      call gtv_setrcd(recrd,nrecs,recln)

      if (passn == 2) then
        nsite = s_ctrl%nbas
        nclass = s_ctrl%nclass
        nspec = s_ctrl%nspec
        nsp = s_ctrl%nspin
        nl = s_ctrl%nl
        nbas = nsite
        nclasp = s_ctrl%nclasp
        allocate(clabl(nclasp))
        allocate(lmxa(nclasp),pnuc(nl,nsp,nclasp),qnuc(3,nl,nsp,nclasp),
     .    pp(6,nl,nsp,nclasp),ves(nclasp),zc(nclasp))
        call dpcopy(s_pot%pnu,pnuc,1,nl*nsp*nclasp,1d0)
        call dpcopy(s_pot%qnu,qnuc,1,3*nl*nsp*nclasp,1d0)
        call dcopy(6*nl*nsp*nclasp,s_pot%pp,1,pp,1)
        call dpcopy(s_pot%ves,ves,1,nclasp,1d0)
        do  j = 1, nclasp
          xx = s_ctrl%dclabl(j)
          call r8tos8(xx,clabl(j))
          is = s_ctrl%ics(j)
          lmxa(j) = s_spec(is)%lmxa
          zc(j) = s_spec(is)%z
        enddo

        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
     .    ves,s_ctrl%initc,s_ctrl%ics,clabl)
        call dpcopy(pnuc,s_pot%pnu,1,nl*nsp*nclasp,1d0)
        call dpcopy(qnuc,s_pot%qnu,1,3*nl*nsp*nclasp,1d0)
        call dpcopy(pp,s_pot%pp,1,6*nl*nsp*nclasp,1d0)
        call dpcopy(ves,s_pot%ves,1,nclasp,1d0)
        deallocate(lmxa,pnuc,qnuc,pp,ves,zc,clabl)
        lasa = s_ctrl%lasa
        lasa = lasa + lasa3 + 8*isw(lasa8)
        s_ctrl%lasa = lasa
C       ldlm = s_ctrl%ldlm
C       tdlm = s_ctrl%tdlm
        return
      else
        call setmorderdefault(1)
      endif

C --- Read input parameters from contents of rcd ---
C     readctrl reads input file data
      call readctrl(prgnam,vrsion(1),ivn)
      if (io_help > 0) then
C        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
C     .    ves,s_ctrl%initc,s_ctrl%ics,clabl)
        call readctrlpq(prgnam,nclasp,nl,nsp,xx,xx,xx,xx,xx,irs,irs,slabl_)
        call cexit(0,1)
      endif
C     Allocate s_site, s_spec, etc
      if (nbasp > nbmx) call rxi('rdctrl: increase nbmx in main to at least',nbasp)
C     Initially allocate arrays with nbmx, since new sites may be added
      call ptr_lat(s_lat,1,'pos',3,nbmx,0,0)
C     call ptr_ctrl(s_ctrl,1,'ips',nbmx,0,0,0)
C     s_ctrl%pos => s_lat%pos
      if (nglobxst('wronsk') < 0) xx(1) = dglob('wronsk',0d0,1)

C     Set switches depending type of program
      lbas = 0
      lcd = 0
C     For now, LMF => fp; no screening; nfp-style hamiltonian
C     cd represented in plane waves
      if (mod(lfp,2) /= 0) then
        lbas = 3
      endif
      if (trim(prgnam) == 'LMMC') then
        lbas = 1
      endif
      asa = .false.
      if (prgnam == 'LM' .or. prgnam == 'LMGF' .or.
     .    prgnam == 'LMPG' .or. prgnam == 'LMCTL') then
        asa = .true.
      endif

C ... Optionally read positions from pos file
      if (cmdopt('--rpos=',7,0,fileid)) then
        call iopos(.false.,-1,fileid(8:),nbasp,pos,s_site)
      endif

C ... Other overrides
      if (cmdopt('--noinv',7,0,fileid)) ctrl_lqp1 = .true.

C ------------------- Copy to structures ----------------------
C --- Copy input to s_bz ---
      s_bz%def = bz_def
      s_bz%dosw = bz_dosw
      s_bz%ef = bz_ef
      s_bz%efmax = bz_efmax
      s_bz%fsmom = bz_fsmom
      s_bz%egap = 0

      if (cmdopt('--getqp',7,0,strn)) bz_lio1 = .true.
      if (cmdopt('--putqp',7,0,strn)) bz_lio2 = .true.
      bz_lio48 = min(bz_lio48,2)
      if (cmdopt('--dmatk',7,0,strn) ) bz_lio48 = 1
      if (cmdopt('--dmatk2',8,0,strn) ) bz_lio48 = 2
      i = 1*isw(bz_lio1)+2*isw(bz_lio2)+4*bz_lio48
      s_bz%lio = i

      s_bz%lcond = bz_lcond
      s_bz%lmet = bz_lmet
      s_bz%lmull = bz_lmull
      s_bz%lopt = 0

      s_bz%lshft = bz_lshft
      s_bz%n = bz_n
      s_bz%ndos = bz_ndos
      s_bz%nevmx = bz_nevmx
      s_bz%nkabc = bz_nabc
      s_bz%nkp = 0

      s_bz%range = bz_range
      s_bz%semsh = bz_semsh
      s_bz%w = bz_w
      s_bz%zval = bz_zval
      s_bz%zinsul = bz_zinsul

C     lasa: 1 Make V from P,Q  2 Make pp  4 ccor  8 free atm
C          16 map  32 nonspherical mpol moms 64 MT corr
C         128 interpretation of sphere Q2; see newrho.f
C         256 how ASA Q1,Q2 are accumulated; see makwts.f
C         512 (spin pol) alpha repsn = (gamma(1) + gamma(nsp))/2
      lasa=4*isw(lasa4)+32*isw(lasa32)+64*isw(lasa64)+128*mod(ham_qasa,4)+512*isw(lasa512)+1024*isw(lasa1024)
     .    +4096*mod(ham_qasa/4,2)
      if (.not. asa) lasa=0
C     lbas: 1 Hamiltonian has no screening transformation
C           2 Hamiltonian is nfp style
C          16 freeze phi,phidot for all species
      j = lbas + 16*isw(frzwf)
C     lcd: 1 freeze core
C          2 non-self-consistent Harris
C          4 represent full potential density on a uniform mesh
C          8 represent full potential density via TCF
C         16 When calculating core, calculate w/ spin-averaged potential
C         32 unused
C         64 (molecules) XC potential by FFT
      k = 1*isw(lcd1)+2*isw(lcd2)+4*isw(lcd4)+8*isw(lcd8)+64*isw(lcd64)
      k = k + 16*isw(lcd16)

      s_ctrl%defm = lat_defm
      s_ctrl%elin = asa_elin
      s_ctrl%lasa = lasa
      s_ctrl%lbas = j
      s_ctrl%lcd = k
      s_ctrl%nccomp = 0
      s_ctrl%maxmem = maxmem

      k = isw(nmto>1)
      if (nmto>1) k = k+2*isw(ham_ewald)
      s_ctrl%lcgf = ctrl_lcgf
      s_ctrl%ldos = ctrl_ldos
      s_ctrl%lfp = lfp + 2*isw(lfp2)
      s_ctrl%lfrce = ctrl_lfrce
      s_ctrl%lgen3 = k
      s_ctrl%ldlm = ctrl_ldlm
C     lbxc  0 no modification of LDA XC potential
C           1 Scale field as 1+bxc(ib), ib=site index
C           2 Scale field as 1+bxc^2(ib), ib=site index
C             bit lbxc can be 1 or 2 but not 3.
C           4 impose constraining fields for vector DLM  (rdctrl)
C             This automatically sets 1+2's bits lbxc to 2.
C           8 (lmf) scale Bxc for smooth density as well as local density
!r         16 (lmf) spin-average potential so LDA part of Bxc=0
      s_ctrl%lbxc = mod(ctrl_lbxc1,4) + 4*isw(ctrl_lbxc4)
     .            + 8*isw(ctrl_lbxc1 >= 10) + 16*isw(ctrl_lbxc16)
      s_ctrl%tdlm = ctrl_tdlm

C     lham  1 (ASA) 2-center
C           1 (molecules) two-panel
C           2 (ASA) 2-c + pert. corr
C           4 (ASA) auto-down-fold
C           8 (ASA) change rep interactively
C          16 (ASA) suppress d hybridization
C          32 (ASA) preserve ortho. evecs
C          64 (ASA) save evecs to disk
C         128 (ASA) gamma-rep
C         256       use true spherical harmonics
C         512       Read shfac file, if it exists
      lham = 1*isw(lham1)+4*isw(lham4)+8*isw(lham8)+16*isw(lham16)+
     .      32*isw(lham32)+64*isw(lham64)+128*isw(lham128)+
     .     256*isw(lham256) + 512*isw(lham512) + lham3
C     ctrl_lmet    1 metal  2 tetrahedron  16 sampling but use tetra to print Ef
Cr                 4 (PGF) V-shift1 is zero
Cr                 8 (PGF) V-shift2 is zero
      lmet = isw(bz_lmet/=0) + 2*mod(ctrl_lmet2,2) +
     .     4*isw(ctrl_lmet4)   + 8*isw(ctrl_lmet8) +
     .    16*mod(ctrl_lmet2/10,2)
C     lncol 1 noncollinear magnetism
C           2 spin spirals
C           4 spin-orbit coupling
C           8 External magnetic field
C          16 mag. forces
C          32 spin-orbit coupling, LzSz only
C          64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
C     If SS, spin-orbit or Bfield, also turn on noncollinear
      k = 2*isw(lncol2)+lSOf*isw(lncol4)+8*isw(lncol8)+16*isw(lncol16)
C     lncol1 = k /= 0
      if (k /= 0) lncol1=T
      if (k == 0) lrsa = 0
      lncol = 1*isw(lncol1) + 2*isw(lncol2) + 4*isw(lncol4) +
     .        8*isw(lncol8) +16*isw(lncol16)+lSzLz*isw(lncol32)+
     .        lSzLzp*isw(lncol64)+lSzLzp0*isw(lncol128) +
     .        lBextz*isw(lBz)+lSOsite*isw(lsos)
      lordn = 0

      s_ctrl%lham = lham
      s_ctrl%lmet = lmet
      s_ctrl%lncol = lncol
      s_ctrl%loptc = ctrl_loptc
      s_ctrl%lordn = lordn
      s_ctrl%lrsa = lrsa
C     The GF codes require CPA to be turned on for SO to work
      if (prgnam == 'LMGF' .or. prgnam == 'LMPG') then
        if ((lncol4 .or. mod(lrel,10)==2) .and. s_ctrl%ldlm==0) s_ctrl%ldlm=12
      endif

C     lrs  switches concerning restart mode.
C         1 Read from restart file
C         2 Read from restart file, ascii mode
C         4 Read from restart file, invoke smshft
C         8 Write new density to restart file
C        16 Write new density to restart file, ascii format
C        32 read site positions from input file
C        64 read starting fermi level from input file
C       128 read starting pnu level from input file
C       256 rotate local density after reading
      call ivset(irs,1,5,0)
      irs(1) = 1
      irs(2) = 1
C     if (cmdopt('--rs=',5,0,strn) .or. cmdopt('--rs',4,0,strn)) then ! Bug: may pick up another switch, e.g. --rsedit
      call readrssw(strn,i,irs)
      ltmp = isanrg(mod(irs(1),10),0,3,'RDCTRL','rs(1), 1s digit ',.true.)
      irs(1) = mod(mod(irs(1),10),4) + 4*getdig(irs(1),1,10) + 8*getdig(irs(1),1,100)
      lrs = 1*mod(irs(1),8) + 8*mod(irs(2),10) + 32*irs(3) + 64*irs(4)+128*irs(5) + 256*mod(irs(1)/8,2)
      if (mod(irs(2)/10,10) == 1) lrs = lrs + 512
C     lqp 1 do not add inversion 2 inverse iteration
      lqp = 1*isw(ctrl_lqp1)+2*isw(ctrl_lqp2)
C     lscr 0 do nothing
C          1 Make P0(0)
C          2 Screen output q and ves
C          3 Screen output ves only
C            Add 10*k to compute intra-site contribution to
C            vbare each kth iteration
C            Add 100*k to compute response function only
C            each kth iteration
C          4 Use model response to screen output q
C            Add 1 to combine mode 1 with another mode
C            Add 10*k to compute intra-site contribution to
C            vbare each kth iteration
C            Add 100*k to compute response function only each kth iteration
C     .  ctrl_lpgf,lqp,lrel,lrs,lscr)
      s_ctrl%lpgf = ctrl_lpgf
      s_ctrl%lqp = lqp
      s_ctrl%lrel = lrel
      s_ctrl%lrs = lrs
      s_ctrl%lscr = lscr
      s_ctrl%lrquad = lrquad
      xx(1) = dglob('lrquad',dble(lrquad),1)

C      ltb switches for empirical tight-binding
C         1 overlap        2 crystal-field     4 ovlp+CF
C         8 add ebarLL    16 forces           32 fij
C        64 not used     128 pressure        256 evdisc
C       512 pair pot    1024 TrH & local E  2048 local rho
C      2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
C      2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
      i = 1*isw(ltb1)+2*isw(ltb2)+4*isw(ltb4)+8*isw(ltb8)+16*isw(ltb16)
     .  +32*isw(ltb32)+64*isw(ltb64)+128*isw(ltb128)+256*isw(ltb256)
     .  +512*isw(ltb512)+1024*isw(ltb1024)+2048*isw(ltb2048)
     .  +4096*isw(ltb4096)+2**13*isw(ltb213)+2**14*isw(ltb214)
     .  +2**15*isw(ltb215)+2**16*isw(ltb216)+2**17*isw(ltb217)
     .  +2**18*isw(ltb218)

      s_ctrl%lstonr = lstonr
C     s_ctrl%lstr = 0
      s_ctrl%lsx = lsx
      s_ctrl%ltb = i
      s_ctrl%lves = isw(lves) + 2*isw(lvshft)

C     lxcf = parameter defining XC functional (rdctrl)
C     See structures.h for definition.
      s_ctrl%lxcf = lxcf
      s_ham%lxcf = lxcf

C     Set modep
      ix(1:3) = 2
      if (prgnam == 'LMPG') ix(3) = 0
      if (prgnam == 'LMMC') ix(1:3) = 0
      s_ctrl%maxit = iter_maxit
      s_ctrl%mdprm = mdprm
      s_ctrl%modep = ix(1:3)
      s_ctrl%nbas = nbas

C     Reset nl
      if (nl /= lmxax+1 .and. io_help == 0) then
        call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .    nl,lmxax+1)
        nl = lmxax+1
      endif

      s_ctrl%nbasp = nbasp
      s_ctrl%nesabc = nesabc
      s_ctrl%nitmv = nitmv
      s_ctrl%nl = nl
      s_ctrl%nmap = 0
      s_ctrl%nsite = nsite
      s_ctrl%nspec = nspec
      s_ctrl%nspin = nsp
      s_ctrl%nvario = nvario
      s_ctrl%omax1 = omax1
      s_ctrl%omax2 = omax2
      s_ctrl%quit = quit
      s_ctrl%plbnd = 0
      s_ctrl%rmaxes = rmaxes
      s_ctrl%rmines = rmines
      s_ctrl%sclwsr = sclwsr
      s_ctrl%sdmod = sdmod
      s_ctrl%sdprm = sdprm
      s_ctrl%sdxsi = NULLI
      s_ctrl%smalit = smalit
      s_ctrl%tol = ctrl_tol
      s_ctrl%wsrmax = wsrmax
      s_ctrl%zbak = zbak
      s_ctrl%lekkl = lekkl

C ... Rotation of lattice via command-line ... supersedes lat_dist
      i = 6
      if (cmdopt('--rot=',i,0,strn)) then
        call a2rotm(strn(i+1:),.false.,iprint()-10,lat_dist)
        lat_ldist = 2
      endif

      if (dalat == NULLR) dalat=0
      s_lat%alat = alat+dalat
      s_lat%as = lat_as
      s_lat%avw = avw
      s_lat%nkdmx = lat_nkdmx
      s_lat%nkqmx = lat_nkdmx
      s_lat%gam = lat_gam
      s_lat%gmax = lat_gmax
      s_lat%nabc = ftmesh
      s_lat%ldist = lat_ldist
      s_lat%dist = reshape(lat_dist(:9),(/3,3/))
      s_lat%plat = reshape(plat, (/3,3/))
      s_lat%plat2 = reshape(slat_plat2, (/3,3/))
      s_lat%rpad = lat_rpad
      s_lat%slat = reshape(lat_slat, (/3,3/))
      s_lat%tol = lat_tol
      s_lat%tolft = tolft
      s_lat%vol = vol
      call dcopy(6,plat,1,xx,1)
      call dcopy(3,platl,1,xx(7),1)
      s_lat%platl = reshape(xx(1:9), (/3,3/))
      call dcopy(3,platr,1,xx(7),1)
      s_lat%platr = reshape(xx(1:9), (/3,3/))

      s_ham%alfsi = alfsi
      s_ham%dabc = dabc
      s_ham%elind = elind
      s_ham%nmto = nmto
      s_ham%kmto = kmto(1:6)
      s_ham%nqsig = 0
      s_ham%ndhrs = 0
C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3,lrsa
      s_ham%lncol = s_ctrl%lncol
      s_ham%lham  = s_ctrl%lham
      s_ham%lgen3 = s_ctrl%lgen3
      s_ham%lrsa = s_ctrl%lrsa
C     Transfer integer parts of sigp
      sigp(1) = sigp_mode
      sigp(2) = sigp_nmin
      sigp(4) = sigp_nmax

      s_ham%nkaph = nkaph
      s_ham%pmax = pmax
      s_ham%pmin = pmin
      s_ham%lsig = lrsig
      s_ham%sigp = sigp
      s_ham%pnudef = pnudef
      s_ham%rdvext = ham_rdvext
      s_ham%qss = ham_qss
      s_ham%rsrnge = rsrnge
      s_ham%rsstol = rsstol
      s_ham%udiag = ham_udiag

C     Parameters for APW
      s_ham%pwmode = pwmode
      s_ham%npwpad = npwpad
      s_ham%pwemin = pwemin
      s_ham%pwemax = pwemax
      s_ham%basopt = basopt
      s_ham%oveps = oveps
      s_ham%ovncut = ovncut

      s_pot%vmtz0 = vmtz
      s_pot%vconst = 0
      s_pot%v0beta = v0beta   !mixing beta for admixing sphere potential V0 defining phi,phidot

C --- Copy input to s_mix ---
      s_mix%b  = smix(2)      !mixing beta
      s_mix%bl = 0            !mixing beta, previous iteration
      s_mix%bv = smix(4)      !extra potential mixing
      s_mix%elind = smix(5)   !Lindhard energy for model screening
      call r8tos8(smix(6),alabl)
      s_mix%fn = alabl
      s_mix%kill = smix(7)   !kill the mixing file after k iterations
      s_mix%lxpot = smix(8)   !for decoupling potential and charge
      s_mix%mmix = smix(9)    !max number prior iter to mix
      s_mix%mode  = smix(10)  !1 for Anderson,  2 for Broyden
      s_mix%model = smix(11)  !previous mixing mode
C     s_mix%n     = smix(12)
      s_mix%nitu  = smix(13)  !max number of LDA+U iterations
      s_mix%nmix  = smix(14)  !actual number prior iter mixed
      s_mix%nsave = smix(15)  !# prior iter to save on disk
      s_mix%tolu = smix(32)   !tolerance for LDA+U
      s_mix%umix = smix(33)   !mixing parm for LDA+U
      s_mix%w(1:3) = smix(34:36) !linear mixing weights
      s_mix%wc = smix(37)

C --- Copy input to s_move ---
      if (lbsprm) then          !Load Bulirsch-Stoer parameters into structure
        prmint(2) = isw(prmint_new)
        prmint(3) = prmint_ts0
        prmint(4) = prmint_tol
        prmint(5) = prmint_mx
        prmint(6) = prmint_mi
        prmint(7:6+prmint_mi) = prmint_nseq(1:prmint_mi)
      endif
      s_move%gyro = 2d0
      s_move%nmodt = gd_nmodt
      s_move%modt = gd_modt
      s_move%ct = gd_ct
      s_move%prmint = prmint
      s_move%kt = move_kt
      s_move%ts = move_ts
      s_move%tsequ = move_tsequ
      s_move%tstot = move_tstot

C --- Copy input to s_optic ---
      s_optic%Nmp = optic_Nmp; s_optic%w = optic_w
      s_optic%dw = optic_dw
      s_optic%window = NULLI
      s_optic%window(1:3) = optic_window
      s_optic%esmr = optic_esmr
      s_optic%ocrng = 0; s_optic%unrng = 0
      s_optic%ocrng(1:2) = optic_ocrng
      s_optic%unrng(1:2) = optic_unrng
      s_optic%esciss = optic_esciss
      s_optic%ltet = optic_ltet
      s_optic%iq = optic_iq
      s_optic%loptic = ctrl_loptc
      s_optic%lpart = optic_mode1
      s_optic%mefac = optic_mefac
      s_optic%ffmt = optic_ffmt
      s_optic%nchi2 = optic_nchi2
      s_optic%imref = optic_imref
      s_optic%kt = optic_kt
      s_optic%axes = reshape(optic_axes, [3,6])
      s_optic%alltrans = isw(optic_alltrans)

C --- Copy input to s_gw ---
      if (cmdopt('--gwcode=0',10,0,strn)) gw_code = 0
      if (cmdopt('--gwcode=1',10,0,strn)) gw_code = 1
      if (cmdopt('--gwcode=2',10,0,strn)) gw_code = 2
      if (cmdopt('--dmft',6,0,strn)) gw_code = 10
C      s_gw%lsigint = 0
      s_gw%mksig = gw_mksig
      s_gw%gcutb = gw_gcutb
      s_gw%gcutx = gw_gcutx
      s_gw%ecutpb = gw_ecutpb
      s_gw%qoffp = gw_qoffp
      s_gw%nband = gw_nband
      s_gw%gsmear = gw_gsmear
      s_gw%code = gw_code
C      s_gw%lsigint = gw_evecl
      s_gw%lgw = 1
      s_gw%nkabc = gw_nabc
      s_gw%lshft = gw_lshft
      s_gw%delre = gw_delre
      s_gw%ecuts = gw_ecuts
      s_gw%nime = gw_nime
      s_gw%deltax = gw_deltax
      s_gw%deltaw = gw_deltaw
      s_gw%pbtol = gw_pbtol
      s_gw%pb1 = '111'
      s_gw%pb2 = '1111'
      s_gw%usebsew = gw_usebsew

C ... Load DMFT parameters
      s_dmft%nkabc = dmft_nabc
      s_dmft%nlohi = dmft_nlohi
      s_dmft%wlohi = dmft_wlohi
      s_dmft%iproj = dmft_proj
      s_dmft%knorm = dmft_knorm
      s_dmft%beta = dmft_beta
      s_dmft%nomg = dmft_nomg
      s_dmft%nomf = dmft_nomf
      s_dmft%nomb = dmft_nomb
      if (cmdopt('--revec',7,0,outs)) then
        s_ctrl%lwsig = 8
        dc = outs(8:8)
        call info0(20,0,0,' rdctrl: read eigenvalues, eigenvectors from evec file')
        if (wordsw(outs,dc,'gw','',j1) > 0) then
          s_bz%nkabc = s_gw%nkabc
          s_bz%lshft = s_gw%lshft
          call info(10,0,0,'%9fusing BZ mesh from GW : nk = %s,%3i',s_bz%nkabc,2)
        endif
      endif

C --- Copy input to s_str ---
C     loka conventions for old 2nd gen and NMTOs
      i = 0
      if (mod(str_mode,100) == 0 .or. str_mode == 2) i = 1
C     Pack iinv parameters
      call dpzero(xx,5)
      xx(1) = iinv_nit
      xx(2) = iinv_ncut
      xx(3) = iinv_tol
      s_str%adec = tcf_adec
      s_str%amode = str_mode
      s_str%drwats = str_drwats
      s_str%iinv = xx(1:3)
      s_str%lmem = str_mem
      s_str%lequiv = isw(str_lequiv1)
      s_str%loka = i
      s_str%lmaxw = str_lmaxw
      s_str%lshow = isw(str_lshow1)
      s_str%mxnbr = str_mxnbr
      s_str%nalf = tcf_nalf
      s_str%nbisi = tcf_nbisi
      s_str%ncupl = tcf_ncupl
      s_str%ndust = tcf_ndust
      s_str%wztcf = tcf_wztcf
      s_str%nkaps = str_nkaps
      s_str%rmax = str_rmax
      s_str%rfit = 0.8d0
      s_str%kaps = str_kaps
      s_str%rmaxg = str_rmaxg
      s_str%ivl = str_ivl
      s_str%delrx = str_delrx
      s_str%mnnbr = str_mnn
      s_str%tolg = str_tolg

C --- Copy input to s_tb ---
      s_tb%alam = fit_alam
      s_tb%alsc = fit_alsc
      s_tb%shft = fit_shft
      s_tb%wg = fit_gausw
      s_tb%wt = fit_wt
      s_tb%fmode = fit_mode
      s_tb%ebfit = fit_ebfit
      s_tb%ndos = fit_ndos
      s_tb%nbfit = fit_nbfit
      s_tb%nbfitf = fit_nbfitf
      s_tb%rmfit = fit_rmfit

C --- Copy input to s_spec ---
      do  j = 1, nspec

C   ... Default values in the absence of declaration
        s_spec(j)%eh3  = -0.5d0
        s_spec(j)%etf  = -1d0
        s_spec(j)%lmxf = 2*nl-2
        s_spec(j)%norp = 2
        s_spec(j)%vmtz = -0.5d0

C   ... Pack species structure
        slabl_(j) = slabl(j)
        s_spec(j)%pb1 = pb1(j)
        s_spec(j)%pb2 = pb2(j)
        s_spec(j)%coreh = coreh(j)
        s_spec(j)%name = slabl(j)
        s_spec(j)%a = spec_a(j)
        s_spec(j)%nr = nr(j)
        s_spec(j)%alpha(1:n0) = alpha(1:n0,j)
        s_spec(j)%coreq(1:2) = coreq(1:2,j)
        s_spec(j)%lxi = lxi(j)
        s_spec(j)%nxi = nxi(j)
        s_spec(j)%exi(1:n0) = exi(1:n0,j)
        s_spec(j)%group = grp(j)
        s_spec(j)%grp2 = grp2(j)
        s_spec(j)%hcr(1:n0) = hcr(1:n0,j)
        s_spec(j)%rsminl = str_rsminl
        if (rsminl(1,j) /= NULLR) then
          do  i = 1, n0
            if (rsminl(i,j) /= NULLR) s_spec(j)%rsminl(i) = rsminl(i,j)
          enddo
        endif
        s_spec(j)%hsfitp(1:2) = hsfitp(1:2,j)
        s_spec(j)%idmod(1:n0) = idmod(1:n0,j)
        s_spec(j)%ehvl(1:n0) = ehvl(1:n0,j)
C       Set idxdn
        call ivset(ix,1,n0*nkap0,1)
        call icopy(1+lmxb(j),idxdn(1,j),1,ix,1)
        s_spec(j)%idxdn(1:n0,1:nkap0) = 1
        s_spec(j)%idxdn(1:1+lmxb(j),1) = idxdn(:,j)
        s_spec(j)%idu(1:4) = idu(1:4,j)
!       s_spec(j)%idm(1:4) = idm(1:4,j)
        s_spec(j)%jh(1:4) = jh(1:4,j)
        s_spec(j)%uh(1:4) = uh(1:4,j)
        s_spec(j)%kmxt = kmxt(j)
        s_spec(j)%kmxv = kmxv(j)
        s_spec(j)%lfoca = lfoca(j)
        s_spec(j)%rsmv = rsmv(j)
        s_spec(j)%lmxa = lmxa(j)
        s_spec(j)%lmxb = lmxb(j)
        s_spec(j)%lmxl = lmxl(j)
        s_spec(j)%lmxpb = lmxpb(j)
        if (nthet(j) == NULLI) nthet(j) = 0
        if (ctrl_ldlm == 0) nthet(j) = 0
        s_spec(j)%nthet = nthet(j)
        s_spec(j)%nang(1:2) = nang(1:2,j)
        s_spec(j)%ncpa  = ncpa(j)
        s_spec(j)%nbeff = nbeff(j)
        if (nbeff(j) /= 0) s_spec(j)%beff(1:n0) = beff(1:n0,j)
        if (ncpa(j) /= 0) then
          s_spec(j)%iscpa(1:n0) = iscpa(1:n0,j)
          s_spec(j)%xcpa(1:n0) = xcpa(1:n0,j)
        endif
C       pack mxcst(j)
        i = 1*isw(mxcst1(j))+2*isw(mxcst2(j))+4*isw(mxcst4(j))
        s_spec(j)%mxcst = i
C       pack orbp(j)
        call dpzero(orbp,n0*2*nkap0)
        call dcopy(n0,rsmh(1,j),1,orbp(1,1,1),1)
        call dcopy(n0,eh(1,j),1,orbp(1,2,1),1)
        call dcopy(n0,rsmh2(1,j),1,orbp(1,1,2),1)
        call dcopy(n0,eh2(1,j),1,orbp(1,2,2),1)
        s_spec(j)%mass = mass(j)
        s_spec(j)%orbp(1:n0,1:2,1:nkap0) = orbp(1:n0,1:2,1:nkap0)

C       Pack P,Q,PZ for both spins
        call dpzero(orbp,n0*6)
        call dcopy(n0*nsp,pnu(1,1,j),1,orbp(1,1,1),1)
        call dcopy(n0*nsp,pz(1,1,j),1,orbp(1,1,2),1)
        call dcopy(n0*nsp,qnu(1,1,j),1,orbp(1,1,3),1)
C     .    orbp(1,1,2),orbp(1,1,3),z(j))
        s_spec(j)%p(1:n0,1:nsp) = pnu(1:n0,1:nsp,j)
        s_spec(j)%pz(1:n0,1:nsp) = pz(1:n0,1:nsp,j)
        s_spec(j)%q(1:n0,1:nsp) = qnu(1:n0,1:nsp,j)
        s_spec(j)%z = z(j)
C     .    colxbs(1,j),radxbs(j),0,0)
        s_spec(j)%colxbs(1:3) = colxbs(1:3,j)
        s_spec(j)%radxbs = radxbs(j)
C     .    rcfa(1,j),rcut(j),rfoca(j),rg(j))
        s_spec(j)%rcfa(1:2) = rcfa(1:2,j)
        s_spec(j)%rcut = rcut(j)
        s_spec(j)%rfoca = rfoca(j)
        s_spec(j)%rg = rg(j)
        s_spec(j)%rham = rham(j)
C     .    rmt(j),rs3(j),rsma(j),rsmfa(j))
        s_spec(j)%rmt = rmt(j)
        s_spec(j)%rs3 = rs3(j)
        s_spec(j)%rsma = rsma(j)
        s_spec(j)%rsmfa = rsmfa(j)
        if (ltbe) then
C     .      iq1(1,j),iq2(1,j),stni(j),tbvso(1,j))
          s_spec(j)%iq1(1:n0) = iq1(1:n0,j)
          s_spec(j)%iq2(1:n0) = iq2(1:n0,j)
          s_spec(j)%stni = stni(j)
          s_spec(j)%vso(1:4) = tbvso(1:4,j)
          s_spec(j)%qpol(1:n0) = qpol(1:n0,j)
          if (ltb217 .or. ltb218) then
           s_bz%nkabc = dvec1
           s_bz%lshft = dvec2
          endif
        endif
        if (prgnam == 'LM') then
          s_spec(j)%iq1(1:n0) = iq1(1:n0,j)
          s_spec(j)%iq2(1:n0) = iq2(1:n0,j)
        endif
C     .    dv(j),eref(j),rham(j),rint(j))
        s_spec(j)%dv = dv(j)
        s_spec(j)%eref = eref(j)
        s_spec(j)%rham = rham(j)
        s_spec(j)%rint = rint(j)
        s_spec(j)%shfac(:) = hfacs(:,j)-1

      enddo

C --- Copy input to s_site ---
      do  j = 1, nsite

C   ... Pack s_site
        s_site(j)%spec = ips(j)
        s_site(j)%class = ips(j)
        s_site(j)%dpole(1:3) = dpole(1:3,j)
        s_site(j)%mpole = mpole(j)
        if (ips(j) > 0) s_site(j)%clabel = slabl(ips(j))
        s_site(j)%pl = ipl(j)
        s_site(j)%plv = plv(j)
        s_site(j)%pos(1:3) = pos(1:3,j)
        s_site(j)%vel(1:3) = vel(1:3,j)
        s_site(j)%vshft = vshft(j)
        s_site(j)%relax(1:3) = irlx(1:3,j)
        s_site(j)%eula(1:3) = eula(1:3,j)
        s_site(j)%rmaxs = rmaxs(j)
        if (ltbe) then
          s_site(j)%ndelta = ndelta(j)
          s_site(j)%delta(1:6) = delta(1:6,j)
        endif
        s_site(j)%ncomp = 1
      enddo

C --- Copy input to s_dmft ---
      if (ncix>0) then
        s_dmft%ncix = ncix       ! Not clear if it is ever used
        s_dmft%nicix = nicix     ! number of independent blocks
        s_dmft%nzsig = dmft_nzsig! Total number of nonzero matrix elements, for one spin now.  See readindmfl.
        s_dmft%lsigim = dmft_imw ! Whether real or imaginary axis
        s_dmft%ncatom = dmft_ncatom ! Number of atoms with correlated blocks
        allocate(s_dmft%ib(ncix),s_dmft%icix(ncix),s_dmft%nzsigi(0:ncix),s_dmft%ndim(nicix),s_dmft%umode(nicix))
        call icopy(ncix,dmft_ib,1,s_dmft%ib,1)
        call icopy(ncix,dmft_icix,1,s_dmft%icix,1)
        call icopy(ncix,dmft_umode,1,s_dmft%umode,1)
        allocate(s_dmft%l(nicix),s_dmft%qsplit(nicix))
        call icopy(nicix+1,dmft_ndsigi,1,s_dmft%nzsigi,1)
        call icopy(nicix,dmft_l,1,s_dmft%l,1)
        call icopy(nicix,dmft_qsplit,1,s_dmft%qsplit,1)
        do  j = 1, nicix
          s_dmft%ndim(j) = 2*s_dmft%l(j)+1
        enddo
        allocate(s_dmft%sigind(dmft_maxdim,dmft_maxdim,nicix,nsp))
        call icopy(size(dmft_sigind),dmft_sigind,1,s_dmft%sigind,1)
        deallocate(dmft_ib,dmft_icix,dmft_sigind,dmft_l,dmft_umode,dmft_qsplit,dmft_ndsigi)
        s_dmft%gammac = dmft_broad/ry2eV
      endif

C ... Copy string outputs to sstrn : amix, gfopt, jobid, mix, mmham, sxopt, symg
      j = str_pack('amix',1,s_strn,iter_amix)
      j = len_trim(iter_amix)          ! Euler angle mixing amix
      j = str_pack('mix',1,s_strn,iter_mix)
      j = str_pack('mmham',1,s_strn,mmham)
      j = str_pack('jobid',1,s_strn,header)
      j = str_pack('gfopt',1,s_strn,gfopt)
      j = str_pack('sxopt',1,s_strn,sxopt)
      j = str_pack('syml',1,s_strn,strn_syml)

C     j = str_pack('blockh',1,s_strn,ham_blockd)

C ... Translate special tokens in SYMOPS to ctrl variables
      call words(symg,n)
      j1 = 1
      lrhop = 0
      do  i = n, 1, -1
        call word(symg,i,j1,j2)
C       Translate GRP2 to s_ctrl%lcd and purge from symg
        if (symg(j1:j1+3) == 'GRP2' .or.
     .      symg(j1:j1+3) == 'grp2') then
          k = 1
          if (symg(j1+4:j1+4) == '=') then
            j = j1+4
            j = a2vec(symg,len(symg),j,2,' ',1,1,1,ix,k)
          endif
          ltmp = isanrg(k,0,3,'RDCTRL','GRP2',.true.)
          s_ctrl%lcd = s_ctrl%lcd +
     .      iand(32*k,32+64) - iand(s_ctrl%lcd,32+64)
          symg(j1:j2) = ' '
C       Translate RHOPOS to s_ctrl%lcd and purge from symg
        else if (symg(j1:j1+5) == 'RHOPOS' .or.
     .           symg(j1:j1+5) == 'rhopos') then
          lrhop = 1
          if (symg(j1+4:j1+6) == '=') then
            j = j1+6
            j = a2vec(symg,len(symg),j,2,' ',1,1,1,ix,lrhop)
          endif
          ltmp = isanrg(lrhop,0,3,'RDCTRL','RHOPOS',.true.)
          symg(j1:j2) = ' '
C       Translate SOC to s_lat%lsym and purge from symg
        else if (symg(j1:j1+2) == 'SOC' .or.
     .           symg(j1:j1+2) == 'soc') then
          if (iand(lncol,5+lSzLz+lSzLzp+lSzLzp0) /= 0) then  ! Symm reduction from SO coupling
          k = 1
          if (symg(j1+3:j1+3) == '=') then
            j = j1+3
            j = a2vec(symg,len(symg),j,2,' ',1,1,1,ix,k)
          endif
          ltmp = isanrg(k,0,3,'RDCTRL','SOC',.true.)
          s_lat%lsym = s_lat%lsym + k - iand(s_lat%lsym,3)
          endif
          symg(j1:j2) = ' '
C       Translate NOSPIN to bit 2^2 in s_lat%lsym and purge from symg
        else if (symg(j1:j1+5) == 'NOSPIN' .or.
     .           symg(j1:j1+5) == 'nospin') then
          k = 1
          if (symg(j1+6:j1+6) == '=') then
            j = j1+6
            j = a2vec(symg,len(symg),j,2,' ',1,1,1,ix,k)
          endif
          ltmp = isanrg(k,0,1,'RDCTRL','NOSPIN',.true.)
          s_lat%lsym = s_lat%lsym + 4*k - iand(s_lat%lsym,4)
          symg(j1:j2) = ' '
C       Translate NOORB to bit 2^3 in s_lat%lsym and purge from symg
        else if (symg(j1:j1+4) == 'NOORB' .or.
     .           symg(j1:j1+4) == 'noorb') then
          k = 1
          if (symg(j1+5:j1+5) == '=') then
            j = j1+5
            j = a2vec(symg,len(symg),j,2,' ',1,1,1,ix,k)
          endif
          ltmp = isanrg(k,0,1,'RDCTRL','NOORB',.true.)
          s_lat%lsym = s_lat%lsym + 8*k - iand(s_lat%lsym,8)
          symg(j1:j2) = ' '
        endif

      enddo
      if (cmdopt('--rhopos',8,0,strn)) lrhop = 1
      if (cmdopt('--soc',5,0,strn)) s_lat%lsym = 1
      s_ctrl%lcd = s_ctrl%lcd +
     .  iand(128*lrhop,128+256) - iand(s_ctrl%lcd,128+256)

C ... Suppress symmetry operations for special circumstances
      lstsym = 0
      if (mod(s_lat%lsym,4) == 0) then  !This is a bit dangerous
        if (lncol /= 0 .and. lncol /= 256) lstsym=1 !lstsym=1: noncollinear case
                                                        !      =2: turn off symops
      endif
C     Switches that automatically turn of all symops
      if ((mdprm(1) >= 1 .and. mdprm(1) <= 3) .or.
     .  cmdopt('--nosym',7,0,strn)) then
C    .  cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn)) then
        symg = 'e'
        lstsym = 2               !lstsym=2: turn off symops
      endif
C     Switches that turn off automatic finder, incl. inversion
      if (lstsym /= 0) then
        i = 1
        do while (i /= 0)
          i = index(symg,'find')
          if (i /= 0) then
            symg(i:i+3) = 'e'
          endif
        enddo
        if (symg == ' ') symg = 'e'  ! suppress ops if none explicit
        lqp = lqp-bitand(lqp,1)+1
      endif
      j = str_pack('symg',1,s_strn,symg)
C      j = len_trim(symg)               ! Symmetry group symg
C      if (j > 0) then
C        call lstra('strn symg',i,o,k)  ! o = index symg has (see ustrn),
C        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
C        sstrn(i:i+j-1) = symg
C      endif
C ... End of copy

C ... Read the basis from the atm file
      if (prgnam=='LMF' .or. prgnam=='LMFGWD' .or. prgnam=='LMFDMFT' .or.
     .    prgnam=='LMFA'.or. prgnam=='LMFGW' .or. prgnam=='LMFGWS') then
      j = nint(basopt(1))
      if (cmdopt('--rdbasp',8,0,strn) .or. j /= 0) then
        fileid = 'basp'
        if (strn(9:12) == ':fn=') then
          fileid = strn(13:13+63) !tk to avoid complaint; replace fileid = strn(13:)
        else
        endif

C       Number of envelope functions before modifications
        k1 = nglob('nkaph')
        k2 = k1
        nphimx = nglob('nphimx')
        if (procid == master) then
        call strip(fileid,i,j)
        call info0(20,1,0,' rdctrl: reading basis parameters from file '
     .    //fileid(1:j))
        i = fxst(fileid(1:j))
        if (i /= 1 .and. prgnam/='LMFA') then
          call rx('file does not exist ... aborting')
        elseif (i /= 1 .and. prgnam=='LMFA') then
          call info0(20,0,0,'%9ffile does not exist ... skipping read')
        endif
        if (i == 1) then
        ifi = fopna(fileid(1:j),-1,1)
        rewind ifi

C       basopt(1):
C        1s   digit 1 or 3 (lmfa) Autogenerate RSMH,EH
C                   2 or 4 (lmfa) Autogenerate RSMH,EH, RSMH2,EH2
C                   1 or 3 (lmf)  Read RSMH,EH,RSMH from basp file
C                   2 or 4 (lmf)  READ RSMH,EH, RSMH2,EH2 from basp file
C        10s  digit 1 (lmfa) Find and estimate PZ from free atom wf
C                     (lmf)  Read P from basp file
C        100s digit 1 (lmfa) Estimate P from free atom wf
C                     (lmf)  Read P from basp file
        if (prgnam=='LMFA') then
          i = 111 + 1000*mod(nint(basopt(1)),10)
        else
          i = 111 + 1000*mod(nint(basopt(1)),10)
        endif
        j = mod(nint(basopt(1))/100,10)
        i = i + 100000*j
        j = mod(nint(basopt(1))/10,10)
        i = i + 10000*j
        if (.not. ioorbp(i,2,1,nspec,s_spec,k2,ifi)) then
          if (prgnam/='LMFA') then
            call rxs2('rdctrl: failed to find BASIS: token in file "',trim(fileid),'"')
          else
            call info0(20,0,0,
     .        '    ...  failed to find BASIS: token in file ... skipping input')
          endif
        endif
C       Update nphimx if LO added
        nphimx = max(nglob('nphimx'),2+mod(k2/10,10))

        k2 = mod(k2,10) + mod(k2/10,10) ! max number of functions of a one l entering into h
        call fclr(' ',ifi)
        endif
        endif
        call mpibc1(k2,1,2,0,'rdctrl','k2')
        call mpibc1(nphimx,1,2,0,'rdctrl','nphimx')
        xx = dglob('nphimx',dble(nphimx),1)
        call bcast_strx(1,xx,xx,xx,xx,xx,xx,s_spec,xx,xx,nspec,0)

C       File read cause number of envelope functions to change
C       If so, reset nkaph
        if (k2 /= k1) then
          xx(1) = dglob('nkaph',dble(k2),1)
          call uspecb(0,-1,s_spec,1,nspec,xx,xx,xx,xx)
          nlmax = (s_ctrl%nl)**2
          xx(1) = dglob('mxorb',dble(k2)*nlmax,1)
          call info2(20,0,0,'%9freset nkaph from %i to %i',k1,k2)
        endif
      endif
      endif

C     Add dalat to alat ... already done
C     s_lat%alat = s_lat%alat+dalat
C     Replicate ctrl->lncol in ham->lncol
      s_ham%lncol = lncol
      lsx1 = mod(lsx,2)
      if (lasa32) s_mix%nsave = 3
C ... Use true spherical harmonics
C     if (bittst(lncol,4) .or. ctrl_loptc /= 0 .or. mod(lrel,10) == 2) then
C      print *, '!! suppress switch to sharm in rdctrl'
      if (bittst(lncol,4) .or. mod(lrel,10) == 2) then
        s_ctrl%lham = IOR(s_ctrl%lham,256)
      endif

C ... Dirac equation requires spin polarization
      if (nsp == 1 .and. mod(s_ctrl%lrel,10) == 2) then
        call rx('rdccat: Dirac equation requires NSPIN=2')
      endif

C ... Suppress inversion when noncollinear magnetism, SX, NLO
      if (lncol+lsx1 /= 0 .or. ctrl_loptc >= 10 .or.
     .   (mod(lscr,10) == 1 .and. prgnam == 'LM'))
     .  lqp = lqp-bitand(lqp,1)+1
      s_ctrl%lqp = lqp
C ... Special pgf initialization
      if (ctrl_lpgf(1) /= 0) then
        s_ctrl%lmet = s_ctrl%lmet - IAND(s_ctrl%lmet,2)
      endif

C ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
      j = IAND(s_ctrl%lham,4)/4
      j = 2*(1-j)
C     No screening => no downfolding; also lmxb<l<=lmxa => 'high'
C     Probably ought to have lmxb<l<=lmxa => 'high' always
      if (IAND(s_ctrl%lbas,1)/=0) j = 3
C     nfp-style basis:
      if (IAND(s_ctrl%lbas,2)/=0) j = j+10
      call suidx(nglob('nkaph'),j,nspec,s_spec)
      call uspecb(0,0,s_spec,1,nspec,xx,xx,xx,i)
C     lpz contains:
C     :1s   digit: maximum number of local orbitals of the 1st
C           type on a particular site
C     :10s  digit: maximum number of local orbitals of the 2nd
C           type on a particular site
C     :100s digit: maximum number of local orbitals of the 3rd
C           type on a particular site
C     :1000s digit: maximum number of local orbitals of any
C           type on a particular site
      xx(1) = dglob('lpz',dble(i/10),1)

C ... Set some global variables
      xx(1) = dglob('nspec',dble(nspec),1)
      xx(1) = dglob('nbas',dble(nbas),1)
      xx(1) = dglob('nbasp',dble(nbasp),1)
      xx(1) = dglob('nsp',dble(nsp),1)
      xx(1) = dglob('nl',dble(nl),1)
      xx(1) = dglob('avw',s_lat%avw,1)
      xx(1) = dglob('lrel',dble(lrel),1)
      lxcf  = s_ctrl%lxcf
      xx(1) = dglob('lxcf',dble(lxcf),1)
      xx(1) = dglob('stdo',dble(stdo),1)
      xx(1) = dglob('stdl',dble(stdl),1)
      xx(1) = dglob('stde',dble(stde),1)
      nspc = 1
      if (bitand(lncol,1+2+4+8) /= 0) nspc = 2
      xx(1) = dglob('nspc',dble(nspc),1)
      xx(1) = dglob('nspc3',dble(nspc),1)
      if (lncol64 .or. lncol128) xx(1) = dglob('nspc3',2d0,1)

C     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid == master) then
      nat = nbas
      do  i = 1, nbas
        j = s_site(i)%spec
        l = s_spec(j)%lmxa
        if (l == -1) nat = nat-1
      enddo
      endif
      call mpibc1(nat,1,2,0,'rdctrl','nat')
      xx(1) = dglob('nat',dble(nat),1)

C ... Set modep
      ix(1) = 2
      ix(2) = 2
      ix(3) = 2
      if (s_ctrl%lpgf(1)/=0) ix(3) = 0
      s_ctrl%modep = ix(1:3)

C ... Count LDA+U blocks (printout only)
      allocate(wk(nbas)); call iinit(wk,nbas)
      call pshpr(1)
      call suldau(nbas,s_spec,s_site,nlibu,k,wk)
      s_ham%nlibu = nlibu
      s_ham%lmaxu = k
      call poppr
      deallocate(wk)

C     Free arrays used to read input
      deallocate(pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .  lmxpb,qpol,stni,tbvso,iq1,iq2,rg,rsma,rfoca,rsmfa,rcfa,nxi,ehvl,hsclsm,
     .  exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,rsminl,
     .  rmt,alpha,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv,
     .  lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv,hsfitp,
     .  nthet,ncpa,iscpa,xcpa,beff,nbeff,nang,hfacs)
      deallocate(pos,vel,eula,vshft,ips,ipl,plv,irlx,mpole,dpole)
      if (ltbe) deallocate(delta,ndelta)

C --- Structural data output in xsf format ---
      if (cmdopt('--wxsf',6,0,strn)) then
        call ioxsf(1,s_lat,s_site,s_spec)
        call rx0('written structural data to file')
      endif

C --- Printout ---
      if (iprint() >= 20 .and. procid == master) then
        do  k = 1, 2

        strn = '  '//prgnam
        if (prgnam == 'LMMC') then
          call awrit3('%N %a:%12pnbas = %i%?#n#%-1j(+%i)##  nspec = %i',
     .      strn,scrwid,0,nbas,nsite-nbas,nspec)
        else
          call awrit5('%N %a:%12pnbas = %i%?#n#%-1j(+%i)##%?#n#+%i#%j#  nspec = %i',
     .      strn,scrwid,0,nat,nbas-nat,nbas-nsite,nsite,nspec)
        endif
        if (ivn(1,1) == 0) then
          call awrit5('%a  verb %i%?#n#,%i#%j#%?#n#,%i#%j#',strn,scrwid,-lgunit(k),iprint(),
     .      iabs(iprint()-iprt(1))+iabs(iprt(2)-iprt(1)),iprt(1),iprt(2)-iprt(1),iprt(2))
        else
!           outs = trim(strn) // '  vn'
!           call lmvsn(ivn,outs)
!           if (ivn(1,1)/=ivn(1,2) .or. ivn(2,1)/=ivn(2,2)) then
!             strn = trim(outs) // '(' // trim(prgnam)
!             call lmvsn(ivn(1,2),strn); i = len_trim(strn) + 1; strn(i:i) = ')'
!           else
!             strn = trim(outs)
!           endif
          call awrit5('%a  verb %i%?#n#,%i#%j#%?#n#,%i#%j',
     .      strn,scrwid,-lgunit(k),iprint(),iabs(iprint()-iprt(1))+iabs(iprt(2)-iprt(1)),
     .      iprt(1),iprt(2)-iprt(1),iprt(2))
        endif
        call lmorder(0,i,irs,irs); xx(1) = dglob('morder',dble(i),1)
        i = IAND(s_ctrl%lham,256)/256*(1+10*i) ! for printout
        lfrzw = isw(IAND(s_ctrl%lbas,16)/=0)

        call awrit8(' special:%10p'//
     .    '%?;n; forces,;;'//
     .    '%?;(n%10)==2; Dirac-eqn,;;%-1j'//
     .    '%?;(n/10)>=1; Dirac-core,;;'//
     .    '%?;(n==1); sharm(m=l:-l),;;%-1j'//
     .    '%?;(n==11); sharm(m=-l:l),;;%-1j'//
     .    '%?;(n==21); sharm(old),;;'//
     .    '%?;n==1; eps^-1,;;%-1j%?;(n>=2); scr-rho-out,;;'//
     .    '%-1j%?;(n>=4);%b(model eps),;;'//
     .    '%?;(n==1|n==2); LM fit to bands,;;%-1j'//
     .    '%?;(n==11|n==12); LM fit to DOS,;;'//
C     .    '%?;n; Order-N:?,;;%-1j'//
C     .    '%?;n==1;%2bEmbedded-Cluster,;;%-1j'//
C     .    '%?;n==2;%2bVanderbuilt,;;'//
     .    '%?;n; APW basis,;;%-1j%?;n>=10;%b(q),;;'//
     .    '%?;n; oveps,;;'//
     .    '%?;n; oveps:n=%-1j%i,;;',
     .    strn,scrwid,0,ctrl_lfrce,lrel,i,
     .    mod(lscr,10),fit_mode,
     .    pwmode,isw(oveps/=0d0),ovncut)
        call awrit1('%a%?;(n>1);%-1j  MPI=%i,;;',
     .    strn,scrwid,0,nproc)
        if (prgnam == 'LMGF') then
          call awrit2('%a%?;n;  CPA,;;'//
     .      '%?;n==10;  transverse J,;;%-1j'//
     .      '%?;n==11;  display J,;;',
     .      strn,scrwid,0,ctrl_ldlm,ctrl_lcgf)
        endif
        call awrit0('%a%b %b',strn,scrwid,0)
        if (strn /= ' special') call awrit0(strn,' ',-80,lgunit(k))

        ltmp = prgnam=='LM'  .or. prgnam=='LMFGWS' .or.
     .         prgnam=='LMF' .or. prgnam=='LMFGWD' .or. prgnam=='LMFGW' .or.
     .         prgnam=='LMFDMFT' .or. prgnam=='LMMC' .or.
     .         prgnam=='LMGF' .or. prgnam=='LMPG' .or.
     .         prgnam=='LMFA' .or. prgnam=='LMCHK'
        if (ltmp) then
          call awrit1('%x pot:%10p'//'%?;n==0; non-rel,;.;',strn,
     .      scrwid,0,mod(lrel,10))

          if (lxcf > 1000) then
            call awrit7('%a%?;%c==,;;%b%a%b;'//
     .        '%?;n>1; spin-pol,;;'//
     .        '%?;n; LDA+U,;;'//
     .        '%?;n; frozen-wave-functions,;;'//
     .        '%?;n==1; nsph-mpol,;;'//
     .        '%?;n; read Sigma,;;'//
     .        '%?;n; make SX,;;'//
     .        '%?;n;; v8 Bessels;'//
     .        '%b %b',
     .        strn,scrwid,0,nsp,nlibu,lfrzw,isw(lasa32),lrsig,lsx,cmdoptswx('--newvesmt','',''))
            call awrit0(strn,' ',-80,lgunit(k))
            call xc_libxc_name(lxcf,nsp,1,outs)
            call awrit0('%x libxc%5f'//outs,strn,scrwid,0)
            call xc_libxc_name(lxcf,nsp,2,outs)
            call awrit0('%a + '//outs,strn,scrwid,0)
            call xc_libxc_type(lxcf,0,ix)
            if (ix(1)>2 .or. ix(2)>2) then
              call awrit0(strn,' ',-80,lgunit(k)); call rx('functional not available, sorry')
            endif
            lwarn = ix(1)/=ix(2)  .and. ix(2)>0 ! functionals from a different family
            call xc_libxc_type(lxcf,nsp,ix)
            lwarn = lwarn
     .         .or. ix(1)/=0 .and. ix(2)/=0 .and. ix(1)/=2 ! Neither contains exchange
     .         .or. ix(1)/=1 .and. ix(2)/=1 .and. ix(1)/=2 ! Neither contains correlation
     .         .or. ix(1)==ix(2) .and. (ix(1)==0 .or. ix(1)==1) ! Duplicate of exchange or correlation
            if (lwarn) call awrit0('%a (warning: mismatch)',strn,len(strn),0)
          else
            call awrit8('%a%?;%c==,;;%b%a%b;'//
     .        '%?;n>1; spin-pol,;;'//
     .        '%?;n; LDA+U,;;'//
     .        '%?;n==1; XC:CA,;;%-1j'//
     .        '%?;n==2; XC:BH,;;%-1j'//
     .        '%?;n==3; XC:PW91,;;%-1j'//
     .        '%?;n==4; XC:PBE,;;'//
     .        '%?;n==1;%b+LMH(gga),;;%-1j'//
     .        '%?;n==2;%b+PW91(gga),;;%-1j'//
     .        '%?;n==3;%b+PBE(gga),;;%-1j'//
     .        '%?;n==4;%b+Becke(gga),;;'//
     .        '%?;n==1; nsph-mpol,;;'//
     .        '%?;n; read Sigma,;;'//
     .        '%?;n; make SX,;;'//
     .        '%?;(n>0); v8-Bessel,;;'//
     .        '%b %b',
     .        strn,len(strn),0,nsp,nlibu,mod(lxcf,100),lxcf/100,isw(lasa32),lrsig,lsx,cmdoptswx('--newvesmt','',''))
          endif
          call awrit0(strn,' ',-len(strn),lgunit(k))

          i = basopt(1); strn = ' '
          if (prgnam=='LMFA') then
            call awrit8(' autogen:%10p'//
     .        '%?;n; mto basis(%-1j%i),;;'//
     .        '%?;n; pz(%-1j%i),;;'//
     .        '%?;n; pnu(%-1j%i),;;'//
     .        '%?;n;%b   Autoread: pz(%-1j%i),;;'//
     .        '%b %b',strn,scrwid,0,mod(i,10),
     .        mod(i/10,10),mod(i/100,10),mod(i/10,10),0,0,0,0)

          elseif (prgnam=='LMF' .or. prgnam=='LMFDMFT' .or.
     .            prgnam=='LMFGWS' .or. prgnam=='LMFGWD' .or. prgnam=='LMFGW' .or.
     .            prgnam=='LMFA' .or. prgnam=='LMCHK') then

            if (lfrzw == 1) then
C              call info2(10,0,0,' float:%11p'//
C     .        'fixed basis'//
C     .        '%?;(flor(n/10)%10>=4);, old extended Pz;;',pnudef,0)

              call awrit8(' float:%11pfixed basis'//
     .        '%?;(flor(n/10)%10>=4);, old extended Pz;;',
     .          strn,scrwid,0,
     .          pnudef,0,0,0,0,0,0,0)

            else
C              call info2(10,0,0,' float:%11pfloat P '//
C     .        '%?;n%10==0;v6;;%-1j%?;n%10==1;LDA;;%-1j%?;n%10==2;GW;;%-1j-style'//
C     .        '%?;(flor(n/10)%4==1);, spin-av Pz;;%-1j'//
C     .        '%?;(flor(n/10)%4==2);, freeze Pz;;%-1j'//
C     .        '%?;(flor(n/10)%10>=4);, old extended Pz;;'//
C     .        '%?;n;;, v6-ebar;',pnudef,lekkl)

              call awrit8(' float:%11pfloat P '//
     .        '%?;n%10==0;v6;;%-1j%?;n%10==1;LDA;;%-1j%?;n%10==2;GW;;%-1j-style'//
     .        '%?;(flor(n/10)%4==1);, spin-av Pz;;%-1j'//
     .        '%?;(flor(n/10)%4==2);, freeze Pz;;%-1j'//
     .        '%?;(flor(n/10)%10>=4);, old extended Pz;;'//
     .        '%?;n;;, v6-ebar;',
     .          strn,scrwid,0,
     .          pnudef,lekkl,0,0,0,0,0,0)

            endif
            call awrit0(strn,' ',-recln,lgunit(k))

            call awrit8(' autoread:%10p'//
     .        '%?;n; mto basis(%-1j%i),;;'//
     .        '%?;n; pz(%-1j%i),;;'//
     .        '%?;n; pnu(%-1j%i),;;'//
     .        '%?;%c==,;%b;;',strn,scrwid,0,
     .        mod(i,10),mod(i/10,10),mod(i/100,10),0,0,0,0,0)
          endif
          i = len(trim(strn))
          if (i > 0) then
            if (strn(i:i) == ':') call awrit8('%11pnone',strn,scrwid,0,0,0,0,0,0,0,0,0)
          endif
          call awrit0(strn,' ',-recln,lgunit(k))
        endif

        if (lxcf < 1000) then
        ltmp = isanrg(mod(lxcf,100),0,4,prgnam,'XC functional',.true.)
        ltmp = isanrg(lxcf/100,0,4,prgnam,'GGA functional',.true.)
        if (mod(lxcf,100) == 2 .and. lxcf/100 > 1 .or.
     .      mod(lxcf,100) /= 2 .and. lxcf/100 == 1) call info0
     .    (10,0,0,'%10f(warning) mixing incompatible functionals')
        endif

        if (bitand(lncol,1+2+4+8+lSzLz+lSzLzp+lSzLzp0+lBextz)/=0) then
          i = s_ctrl%sdmod
          if (.not. bittst(lncol,16)) i = -1
          call awrit8(' noncoll: '//
     .    '%?;n; Non-coll,;;'//
     .    '%?;n; L.S coupling,;;%-1j%?;(n==128);%b (pert),;;'//
     .    '%?;n; LzSz,;;%-1j%?;(n>32);%b+L.S(pert),;;'//
     .    '%?;(n==8); B-field,;;%-1j%?;(n==128); BzSigz,;;'//
     .    '%?;n; spin-spiral,;;'//
     .    '%?;n; mag-forces:;;'//
     .    '%?;(n%10)<2&(n%10)>=0; relax,;;'//
     .    '%?;(n%10)>1; spin-dynamics,;;'//
     .    '%a%b %a',strn,scrwid,lgunit(k),bitand(lncol,1),
     .      bitand(lncol,4+lSzLzp0),bitand(lncol,lSzLz+lSzLzp),
     .      bitand(lncol,8+lBextz),bitand(lncol,2),bitand(lncol,16),i,i)
        endif

        if (asa) then
          lham = s_ctrl%lham
          lcd = s_ctrl%lcd
          call awrit8('%x asa:%10p'//
     .      '%?;n<>4; no-ccor,;;'//
     .      '%?;n==128; gam-rep,;;%-1j'//
     .      '%?;n==128+512; gamma-bar-rep,;;'//
     .      '%?;n; two-c-H,;;'//
     .      '%?;n;%b + pert-ev,;;'//
     .      '%?;n; map,;;'//
     .      '%?;n; indep-vmix,;;'//
     .      '%?;n; mt-corr,;;'//
     .      '%?;n; frozen core,;;'//
     .      '%a%b ',strn,scrwid,0,
     .      bitand(lasa,4),
     .      bitand(lham,128)+bitand(lasa,512),
     .      bitand(lham,3),bitand(lham,2),
     .      bitand(lasa,16),
     .      IAND(s_mix%lxpot,3),
     .      bitand(lasa,64),bitand(lcd,1))
          call awrit8('%a%?;(n);:%10p;,;'//
     .      '%?;(n==0); Q2:Power-moments,;;'//
     .      '%?;(n==1); phidot:outward-intgr,;;'//
     .      '%a%b ',strn,scrwid,0,
     .      isw(strn==' asa'),
     .      mod(ham_qasa,4),
     .      mod(ham_qasa/4,2),
     .      4,5,6,7,8)
          if (strn /= ' asa') call awrit0(strn,' ',-scrwid,lgunit(k))
        endif

        if (prgnam(1:2) == 'TB') then
          xx(1:7) = s_ctrl%mdprm
          ltb = s_ctrl%ltb
          ix(1) = nint(xx(1))
          ix(2) = nint(xx(2))
          xx = s_str%rmax
          call awrit1('%x TB: %11prmaxh = %d,',strn,scrwid,0,xx)
          call awrit6('%a'//
     .      '%?;n; s-c multipoles,;;'//
     .      '%?;n; s-c multipoles: MRS theory,;;'//
     .      '%?;n; read del,;;'//
     .      '%?;n; n-orthog TB,;;'//
     .      '%?;n==4; m-stat: Conj. grad.;;%-1j'//
     .      '%?;n==5; m-stat: F-P;;%-1j'//
     .      '%?;n==6; m-stat: Broy;;%-1j'//
     .      '%?;n==1; MD (NVE),;;%-1j'//
     .      '%?;n==2; MD (NVT),;;%-1j'//
     .      '%?;n==3; MD (NPT),;;'//
     .      '%?;n; pair-only,;;',
     .      strn,scrwid,0,bitand(ltb,2**15),bitand(ltb,2**12),
     .      bitand(ltb,2**16),bitand(ltb,1),ix(1),bitand(ltb,512))
          if (ix(1) >= 4) then
            call awrit2('%a'//
     .      '%?;n; rlx-vol,;;'//
     .      '%?;n; i/o-hess,;;',
     .      strn,scrwid,0,bitand(ix(2),1),bitand(ix(2),2))
          endif
          call awrit5('%a'//
     .      '%?;n; trh,;;'//
     .      '%?;n; rho,;;'//
     .      '%?;n; spin-orb,;;'//
     .      '%?;n; crysf,;;'//
     .      '%?;n;%b+ovlp-cf,;;'//
     .      '%a%b %a',strn,-scrwid,-lgunit(k),
     .      bitand(ltb,1024),bitand(ltb,2048),
     .      bitand(lncol,4),bitand(ltb,2),bitand(ltb,4))
        endif

        if (prgnam(1:3) == 'LMF') then
          call awrit1(' special:%10p'//
     .      '%?;n==1; core-level-optics,;;'//
     .      '%b %b',strn,scrwid,0,s_optic%cls)
          if (strn /= ' special:') call awrit0(strn,' ',-80,lgunit(k))
        endif

        if (ctrl_loptc /= 0 .and. k == 1) then
          i = 1
          if (optic_dw(2) /= 0) i = 11
          call dosmsh(-i,optic_window,optic_dw,0,j,xx)
          call info8(20,0,0,' optics:  '//
     .      '%?#n>0# Im(eps(w)),##%-1j'//
     .      '%?#n>=10# NLO,##%-1j'//
     .      '%?#(n==8)# noneq absorption,##%-1j'//
     .      '%?#(n==9)# noneq emission,##%-1j'//
     .      '%?#(n<0&n>-5)# JDOS,##%-1j'//
     .      '%?#(n==-5|n==-6)# DOS,##'//
     .      '%?#(n==1|n==3)# band-decomp,##%-1j'//
     .      '%?#(n==2|n==3)# k-decomp,##'//
     .      ' %i pts in [%d,%d]'//
     .      '%?#n>10# : dw=%;3g, omg=%;3g##',
     .      ctrl_loptc,optic_mode1,j,optic_window(1),optic_window(2),
     .      i,optic_dw(1),optic_dw(2))
          if ((ctrl_loptc == -5 .or. ctrl_loptc == -6) .and.
     .        optic_dw(2) /= 0) call rxi('RDCTRL: quadratic '//
     .      'energy mesh not allowed with OPTICS_MODE =',ctrl_loptc)
          if (abs(s_optic%loptic)==8 .or. abs(s_optic%loptic)==9) then
            if (s_optic%ltet/=0) call rxi('rdctrl : tetrahedron not allowed with optics mode',s_optic%loptic)
            if (s_optic%imref(1) > s_optic%imref(2)) call info0(2,0,0,'%11fWarning! imref(1) > imref(2)')
          endif
        endif

c       ltmp = prgnam=='LM'   .or. prgnam=='LMMC'  .or.
        ltmp = prgnam=='LM'   .or. prgnam=='LMFGWS' .or.
     .         prgnam=='LMF'  .or. prgnam=='LMFGWD' .or.
     .         prgnam=='LMFGW' .or. prgnam=='LMFDMFT' .or.
     .         prgnam=='LMGF' .or. prgnam=='LMPG' .or.
     .         prgnam=='LMDOS'.or. prgnam=='TBE'

        if (ltmp) call awrit7(' bz:%10p'//
     .    '%?;n; metal,; nonmetal,;'//
     .    '%-1j%?;n>1;%b(%-1j%i),;;'//
     .    '%?;n; tetra,;;'//
     .    '%?;n; get-qp,;;'//
     .    '%?;n; invit,;;'//
     .    '%?;n; dens-mat,;;'//
     .    '%?;(n>0); %-1jmull=%i,;;'//
     .    '%?;n; fixed-spin-mom,;;%b ',
     .    strn,scrwid,lgunit(k),
     .    s_bz%lmet,bitand(lmet,2),
     .    IAND(s_bz%lio,1),IAND(s_ctrl%lqp,2),IAND(s_bz%lio,4+8),
     .    s_bz%lmull,isw(s_bz%fsmom/=NULLR))

        if (prgnam == 'LMPG') then
          call awrit2(' pgf:%10p'//
     .    '%?;(n==1); self-consistency;;%-1j'//
     .    '%?;(n==2); self-consistent leads;;%-1j'//
     .    '%?;(n==3); k(E), left lead;;%-1j'//
     .    '%?;(n==4); k(E), right lead;;%-1j'//
     .    '%?;(n==5); transmittance;;%-1j'//
     .    '%?;(n==7); reflectance;;'//
     .    '  via %?;(n==0);embedding;LU decmpsn;',
     .      strn,scrwid,lgunit(k),
     .      s_ctrl%lpgf,
     .      mod(s_ctrl%lpgf(2),10))
        endif

        if (ncix>0) call awrit7(' dmft:%11p'//
     .    '%i blocks '//
     .    ' (%i independent)',
     .    strn,scrwid,lgunit(k),
     .    ncix,nicix,3,4,5,6,7)
        ltmp = prgnam=='LM' .or. prgnam=='LMGF' .or. prgnam=='LMPG'
        enddo
      endif

C --- Sanity checks and other initialization ---
      if (IAND(s_ctrl%ldos,8)/=0 .and. IAND(s_ctrl%ldos,4+2)/=0)
     .  call rx('inconsistent DOS options')
      call rxx(s_ctrl%lstonr/=0 .and. nsp==2,
     .  'Stoner model not compatible with nsp=2')
      call rxx(s_ctrl%lstonr/=0 .and.IAND(s_ctrl%lham,4)/=0,
     .  'Stoner not compatible with ADNF')
      if ((lncol64 .or. lncol128) .and. (lncol1.or.lncol2.or.lncol4))
     .  call rx('cannot use SO=3 with noncollinear magnetism')

C --- Summarize DMFT blocks ---
      if (ncix>0 .and. iprint() > 41) then
        call info(1,0,0,' DMFT: %,5i cix blocks, %i inequivalent',ncix,nicix)
        k = min(4,ncix)
        call arrprt(' cix site  l icix |','%,4i%,4i%,4i%,4i','Iiii',
     .    ncix,0,k,0,'  | ',i,dmft_ib,dmft_l,dmft_icix,i,i,i,i)
        if (maxval(dmft_ib)>nbas) call rx(trim(prgnam)//': cix block exceeds nbas')
        call rx0('done with dmft printout')
      endif

C --- Check and order principal layers, and sites by PL ---
      call susite(s_ctrl,s_ham,s_pot,s_lat,s_spec,s_site)
      nsite = s_ctrl%nsite
      nbasp = s_ctrl%nbasp
      xx(1) = dglob('nbasp',dble(nbasp),1)
C     call ptr_lat(s_lat,2,'pos',3,nsite,0,0)
C     Copy equivalent data
      s_ctrl%pos => s_lat%pos
      s_ham%lham = s_ctrl%lham

      if (procid == master) then
      if (iprint() >= 20) then
        if (lstsym == 1) then
          i = str_pack('symg',-2,s_strn,strn)
          write(stdo,357) trim(strn)
  357     format(/' Automatic symmetry finder turned off.  Use: ',a)
        elseif (lstsym == 2) then
          write(stdo,358)
  358     format(/' Symmetry operations suppressed')
        endif
      endif
      endif
C     Broadcast spec structure to fix floating orbitals case
C      call mpibc1(w(osspec),nspec*nint(dval(w(osspec),1)),4,0,
C     .  'rdctrl','sspec')
      call bcast_strx(1,xx,xx,xx,xx,xx,xx,s_spec,xx,xx,nspec,0)

C --- Debugging printout ---
      if (io_help == 0 .and. io_show > 1) then
        print "(/' ---------- contents of s_strn ------------')"
        i = str_pack('amix',-2,s_strn,strn)
        print *, 'amix:   ', trim(strn)
C        i = str_pack('blockh',-2,s_strn,strn)
C        print *, 'blockh: ', trim(strn)
        i = str_pack('mix',-2,s_strn,strn)
        print *, 'mix:    ', trim(strn)
        i = str_pack('gfopt',-2,s_strn,strn)
        print *, 'gfopt:  ', trim(strn)
        i = str_pack('mmham',-2,s_strn,strn)
        print *, 'mmham:  ', trim(strn)
        i = str_pack('symg',-2,s_strn,strn)
        print *, 'symg:   ', trim(strn)
        i = str_pack('sxopt',-2,s_strn,strn)
        print *, 'sxopt:  ', trim(strn)

        call info8(0,1,0,' --- s_ctrl structure ---%N '//
     .    'lham %i  lmet %i  maxit %i  lncol %i  '//
     .    'lfrce %i  lxcf %i  lrs %i  lqp %i',
     .    s_ctrl%lham,
     .    s_ctrl%lmet,
     .    s_ctrl%maxit,
     .    s_ctrl%lncol,
     .    s_ctrl%lfrce,
     .    s_ctrl%lxcf,
     .    s_ctrl%lrs,
     .    s_ctrl%lqp)
        call info8(0,0,0,' '//
     .    'ldos %i  loptc %i  lbas %i  lcd %i  lrel %i  tol%3:1g',
     .    s_ctrl%ldos,
     .    s_ctrl%loptc,
     .    s_ctrl%lbas,
     .    s_ctrl%lcd,
     .    s_ctrl%lrel,
     .    s_ctrl%tol,0,0)
        if (asa) then
          call info8(0,0,0,' '//
     .      'lasa %i  lcgf %i  lpgf%2:i  lsx %i  '//
     .      'lves %i  lgen3 %i',
     .      s_ctrl%lasa,s_ctrl%lcgf,s_ctrl%lpgf,s_ctrl%lsx,
     .      IAND(s_ctrl%lves,1),s_ctrl%lgen3,0,0)
          if (s_ctrl%lncol /= 0) then
            call info8(0,0,0,' '//
     .        'sdprm%5:1g  sdxsi%3:1g',
     .        s_ctrl%sdprm,s_ctrl%sdxsi,0,0,0,0,0,0)
          endif
        endif
         call rx0('done show')
      endif
C     call rx0('done')

C ... Artificial declarations avoid problems with DEC compiler
      call ptr_lat(s_lat,1,'igv2',1,1,0,xx)
      call ptr_ham(s_ham,1,'qsig',1,1,xx)

      end subroutine rdctrl
