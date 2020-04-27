C#define DEBUG
      subroutine makegw(s_ctrl,s_site,s_spec,s_lat,s_ham,s_pot,s_bz,s_gw,nbas,ppn)
C- Driver to set up GW
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfp lwsig
Co     Stored:     lfp
Co     Allocated:  *
Cio    Elts passed:lcd
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos pnu pz class clabel sighh sighk sigkk pihh
Ci                 pihk pikk sighhx sighkx sigkkx pihhx pihkx pikkx
Ci                 sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:v0 pihhx tauhh pikkx taukk pihkx tauhk
Cio    Passed to:  chkgwin sugwin wlattc mk_hamindex hambl augmbl
Cio                bstrux smhsbl hsibq hsubblock hsibq2 hsibq4 hambls
Cio                makusq pusq1 gwcphi pwmat pwmat2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa pz p idxdn a nr z rmt name pb1 pb2 lmxb orbp
Ci                 kmxt rsma ngcut lmxl
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:name
Cio    Passed to:  chkgwin sugwin wlattc mk_hamindex uspecb hambl
Cio                augmbl bstrux smhsbl hsibq tbhsi hsubblock hsibq2
Cio                hsibq4 hambls makusq pusq1 gwcphi pwmat pwmat2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nabc gmax npgrp nsgrp vol igv2 kv2 kv
Ci                 igv awald tol nkd nkq ng napw
Co     Stored:     napw igv2 kv2 igv
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: s_sym symgr nsgrp igv2 vol kv2 qlat igv gv ag pos
Cio                plat cg indxcg jcg cy qlv dlv kv
Cio    Passed to:  chkgwin sugwin sugvec mk_hamindex hambl augmbl
Cio                bstrux hxpbl ghibl hklbl gklbl hxpgbl ghigbl hklgbl
Cio                smhsbl hhibl phhibl hsmbl hsibq hambls makusq pusq1
Cio                pwmat pwmat2
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lsig nqsig oveps ovncut nlmto ldham pwmode pwemin
Ci                 pwemax ndham ndhrs nprs sigp eseavr
Co     Stored:     *
Co     Allocated:  qsig
Cio    Elts passed:qsig lncol iprmb eseavr lrsa iaxs hrs
Cio    Passed to:  sugwin mk_hamindex hambl augmbl hambls sopert3 makusq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vesrmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qval smpoth smvcnst lcplxp smpot
Cio    Passed to:  *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc lshft
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  chkgwin sugwin
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  gcutb code nkabc mksig lgw nband qoffp gcutx ecuts
Ci                 nime delre deltax deltaw pbtol gsmear
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nkabc lshft
Cio    Passed to:  chkgwin sugwin
Ci Inputs
Ci   sblockh: used to zero out M.E. of orbital subblock with neighbors
Ci   nbas  :size of basis
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   ppn   :NMTO pot pars; see potpus.f
Co Outputs
Co   Files written according to jobgw
Co   The following shows what files are written for jobgw=1
Co   and the records stored in each file.
Co .. gwcode=0,2:
Co   gw1:  evals of Hlda+sigma-vxc
Co        *for each q-point and spin:
Co         q, evl(i:ndimh)
Co   gw2:  evals of Hlda+sigma-vxc-vxc
Co        *for each q-point and spin
Co         q, evl(i:ndimh)
Co   gwb:  Information about eigenfunctions, matrix elements
Co         nat,nsp,ndima,ndham,alat,qlat,ef0,nqbz,plat,nqnum,nqi
Co         lmxa(1:nat), bas(1:3,1:nat)
Co         ngpmx,npmbx  -- largest number of G vectors for psi, vcoul
Co        *for each q-point and spin:
Co           q, ndimh
Co           evl, cphi
Co           ngp, ngc
Co           ngvecp, ngvecc, pwz  (G vectors for psi,vcou; PW expansion of z)
Co   gwa:  site data.
Co        *for each site:
Co           z, nr, a, b, rmt, lmaxa, nsp, ncore
Co           konfig(0:lmaxa) : note
Co           rofi
Co          *for each l, spin isp
Co             l, isp
Co             radial w.f. gval: phi
Co             radial w.f. gval: phidot
Co             radial w.f. gval: phiz    written if konfig(l)>10
Co             *for each l, spin, konf
Co                icore, l, isp, konf, ecore(icore)+vshft
Co                gcore(1:nr,1,icore)
Co   evec: eigenvectors.
Co         ndham, nsp, nnn, nqnum
Co        *for each q-point and spin:
Co           q, evec(1:ndimh,1:ndimh)
Co   vxc:  matrix elements of XC potential
Co         ndham, nsp, nnn
Co        *for each q-point and spin:
Co           q, vxc
Cs Command-line switches
Cs   --gradpb     : Also make gradient product basis
Cs   --checkrpb   : Check completeness of product basis
Cs   --checkrpb=2 : Detailed completeness check
Cs   --stop=prbas : Stop after constructing product basis
Cl Local variables
Cl   compat:compatiblitity modes to sync with old GW code
Cl         :2^0 = 1 => read from QGpsi
Cl         :2^1 = 1 => (not used now) reserve for writing files gwb,gwa
Cl         :2^2 = 1 => Read product basis info from GWinput
Cl         :2^3 = 1 => Read qp calling rdQIBZ
Cl         :2^4 = 1 => Write matrix elements of radial product basis functions to files PPBRD_V2_xx
Cl         :2^5 = 1 => Write Coulomb integrals to files Vcoud.xx
Cl   lpdiag:0 use standard diagonalization (zhev)
Cl         :1 use parallel diagonalization (pzhev)
Cl         :2 diagonalization done internally (hambls)
Cl   lsig  :switch to create files vxc and evec for making sigma
Cl   lwvxc :1 or 2, write to evec and vxc files for irreducible qp on regular mesh
Cl         :2 write evec file for all irreducible qp
Cl   lsapw :T, APW's indices must be modified because q vector
Cl             has been shortened (see comments before call to pwmat)
Cl   nnn   :number of qp in full BZ
Cl   nqnum :number of qp in full BZ * number of 'special gamma points'
Cl         :or generally the number of qp at which eigenfunctions calc.
Cl   nqibz :number of qp in irreducible BZ on regular mesh
Cl   nqi   :total number of qp in irreducible BZ, including offset gamma mesh
Cl   ngp   :no. G vectors for eigenfunction expansion (depends on q-pt)
Cl   ngc   :no. G vectors for coulomb interaction (depends on q-pt)
Cl   ispc  :2 when working on (2,2) block of noncollinear hamiltonian;
Cl         :otherwise 1
Cl   ipb   :index to true basis (excluding floating orbitals)
Cl         :given site index including those orbitals
Cl   ndima :number of augmentation channels
Ci   ndham :dimensioning parameter, at least as large as largest
Ci         :hamiltonian dimension
Ci   iprocq:procid for q where file data is generated
Ci         :In the MPIK case each processor writes a portion of evec and vxc
Ci         :to a procid-specific file (name has _procid appended)
Ci tprodbas:Switch to generate product basis
Ci         :1s digit radial product basis
Ci         :1  use valence partial waves only
Ci         :>1 Product basis for merged valence and core partial waves
Ci         :2  Modify noccc, noccv according to Kotani's prescription.  Rethink for future
Ci         :4  noccc, noccv taken from GWinput
Ci         :10s digit Matrix elements of partial waves and radial product basis
Ci         :1  make orthonormalized product basis and matrix elements with B
Ci         :2  ditto, orthonormalized product basis compatible with old GW code
Ci         :100s digit Full product basis; check completeness
Ci         :1  make full product basis
Ci         :2  check completeness
Cb Bugs
Cr Remarks
Cu Updates
Cu   24 Sep 18 Makes gradient product basis, radial part
Cu   07 Aug 18 Start on makegw, adapting from sugw
C ----------------------------------------------------------------------
      use structures
!     use bzdata
      implicit none
C ... Passed parameters
      integer nbas,n0,nkap0,nppn
      parameter (n0=10, nppn=12, nkap0=4)
      integer lh(n0)
C     real(8) sab(nab,n0,*)
      real(8) rsml(n0),ehl(n0),ppn(nppn,n0,nbas)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
      type(str_gw)::    s_gw
C     type(str_symops),allocatable::  s_sym(:)
C ... Dynamically allocated local arrays
      integer, allocatable :: kpproc(:) ! q loop blocking indices for multithreaded case
      integer, allocatable :: iprocq(:) ! Written to file mpiqlst
      integer,allocatable :: nlindx(:,:,:),konfa(:,:),lcuta(:)
      real(8),allocatable:: ww(:),rofi(:),rwgt(:),gcore(:,:,:),gval(:,:,:,:,:)
      real(8),pointer:: gvala(:,:,:,:,:),gtota(:,:,:,:,:),gtoto(:,:,:,:,:),gvalo(:,:,:,:,:)
      real(8),pointer:: gcora(:,:,:,:,:)
      real(8),allocatable:: evl(:,:),cphin(:,:,:),qpgw(:,:)
      complex(8),allocatable:: ham(:,:,:),ovl(:,:,:),evec(:,:),vxc(:,:,:),vxcsig(:,:,:,:),
     .  ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:),ppovld(:),testc(:,:),
     .  testcd(:),cphi(:,:,:),cphi_p(:,:,:),geig(:,:,:),geig_p(:,:,:),aus(:)
      integer,allocatable, target :: ngvecc(:,:,:)

C ... Local parameters
      logical :: bandmode=.false.,endofline=.false.
      integer :: tprodbas=0, lgrad=0,lwvxc=0
      logical rank0  ! Equal to  procid == master
      logical lsapw,passl
      logical keephi ! If T, stores partial waves into large array with all sites
      integer ierr,im2,konf,mnevl,nline,nqibz,gwcode,i,i1,i2,iat,ib,ibr,ic,icore,
     .  ifeigen,ifi,ifiqg,ifiqgc,ifqeigen,ifsyml,iix,iline,im1,imx,ipqn,ipr,iq,
     .  iqibz,irr,is,isp,ispc,j,jb,jfi,job,jobb,jsp,k,k1,k2,k3,l,lchk,ldim,lmaxa,
     .  lmxax,lpdiag,lrsig,lso,lwsig,mx,mxint,n1,n2,n3,napws,nat,nclasx,
     .  ncore,ndham,ndima,ndimh,ndimhs,ndimhx,ndimhxs,nev,nevl,ngc,ngp,ngp_p,
     .  ngpmx,nkp,nlinemax,nlmax,nlmto,nmb,nmcore,nmx,nn1,nn2,nnn,npgrp,nphimx,
     .  ndlmto,npmbx,npqn,nqbz,nqi,nqibze,nqnum,nqnumx,nqtot,nr,nrmx,nsgrp,
     .  nsp,nspc,nspc3,ntmb,ovncut,stdo,compat
      integer procid, master, mpipid, nproc ! For MPIK
      integer iflband(2),inn(3),ipb(nbas),konfig(0:n0),ngabc(3),ivec(10),nlnmaug(nbas,2)
      integer, parameter :: LW9=9,LW19=19
      integer n1q,n2q,n3q,nkgw(3)
      equivalence (n1q,nkgw(1)),(n2q,nkgw(2)),(n3q,nkgw(3))
      double precision dnn(3)
C     integer lshft(3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      real(8) q(3),QpGcut_psi,QpGcut_cou,dum,qpscale,xx(5),gmax,gcutb,
     .  pnu(n0,2),pnz(n0,2),ecore(50),a,z,rmt(nbas),b,vshft,alat,alfa,ef0,
     .  plat(3,3),qlat(3,3),qp(3),qpos,qx(3),q_p(3),epsovl,qb(3,3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
C     integer ndhamx,nspx,nnnx
      character strn*120,ext*20,dc*1
      integer,allocatable:: nevls(:,:)
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      procedure(logical) :: cmdopt,latvec
      procedure(integer) :: bitand,iprint,maxnphi,parg2,idalloc,allocvb,nglob,lgunit
      procedure(integer) :: fopna,fopnx,fopng,iclbas,isw,wordsw,a2vec,iofa
      procedure(real(8)) :: dlength,dval
C ... For PW basis
      integer,parameter:: NULLI = -99999
      real(8),parameter:: NULLR =-99999
      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax,eomin,ebot,etop
      complex(8), parameter :: zer=(0d0,0d0),one=(1d0,0d0)

C ... for reading self-energy
      integer nqsig

C ... for reading BZDATA (job=6,7)
      type str_bzdata

      integer    ::   nqbz     ! Number of k-points in the full BZ
      integer    ::   nqibz    ! Number of k-points in the irreducible part
      integer    ::   nqbzw    ! Dimensions qbzw,ib1bz.  Should be (n1+1)*(n2+1)*(n3+1)
      integer    ::   ntetf    ! Number of tetrahedra in the full BZ
      integer    ::   nteti    ! Number of tetrahedra in the irreducible part
      integer    ::   ngrp     ! Number of symmetry operations
      integer    ::   nqibz_r  !
      integer    ::   n123(3)  ! Number of k-divisions along the reciprocal lattice vectors
      integer,pointer::
     .  idtetf(:,:),           ! corners of tetrahedra for full BZ
     .  ib1bz(:),              ! maps k-points in the padded qp list (qbzw) to the original list
     .  idteti(:,:),           ! corners of tetrahedra for irreducible BZ
     .  nstar(:),              ! star of k
     .  irk(:,:),              ! irreducible
     .  nstbz(:)               !
      real(8),pointer::
     .  qbz(:,:),wbz(:),       ! k-points in the full Brillouin zone, and weights
     .  qibz(:,:),wibz(:),     ! k-points in the irreducible Brillouin zone, and weights
     .  qbzw(:,:),             ! k-points on the regular mesh, padded to repeat boundary points
     .  qibz_r(:,:)
      end type str_bzdata
      type(str_bzdata)::  s_bzdata

C ... External calls
      external atwf2l,awrit1,awrit2,awrit3,awrit8,bndconn_v2,chcase,
     .         daxpy,dgemm,dinv33,dpcopy,dpzero,dscal,dvset,fclose,fclr,
     .         fexit2,fftz30,getpr,gintsl,gvlst2,gwcphi,hambls,
     .         iinit,info0,info2,info5,iosigh,makusq,mk_hamindex,
     .         mshnrl,orbl,phmbls,poppr,prtev,pshpr,pvsug1,
     .         pwmat2,radmsh,radwgt,rx0,rxi,shorps,sitepack,sopert3,
     .         spex_read_k,sugvec,tcn,tcx,uspecb,wlattc,wsymops,zcopy,
     .         zgemm,zhevo,zinit,zqinvb,zslabl

C ... for band plotting
      real(8),allocatable :: ovvx(:,:,:)
      real(8) ::ovvpp
      real(8), allocatable :: double_dummy(:)

C ... for testing product basis
      logical :: lwriteVcoulTK = .false.
      logical :: lwrhoc = .false. ! If .true., write core densities (job 0)
      integer nl,nqc,ncoremx,maxcor(2),ndrphi,nradpb(3),ndrpbg(3),lcutmx,nlcutmx,jat,js,nlmi,nlmj
      integer iqfbz,nnegevl,ndrpb,mxnpba,ntpb,ntpba(0:nbas),ntgpba(0:nbas,2)
      equivalence (ncoremx,maxcor(1))
      integer, allocatable :: ncwf(:,:,:),ncwf2(:,:,:),nrpb(:,:),nrgpb(:,:,:)
      integer, pointer :: npmb(:),nrphi(:,:),nrphiv(:,:),nrphic(:,:),
     .                    nocc(:,:,:),noccv(:,:,:),noccc(:,:,:),
     .                    nunoccv(:,:,:),nunoccc(:,:,:),nunocc(:,:,:)
      real(8), allocatable :: rprodb(:,:),rgprodb(:,:)
      real(8), allocatable :: tolbas(:),rojb(:,:,:,:),sgbb(:,:,:,:),sgbc(:,:,:),
     .                        qc(:,:),qibze(:,:),wibze(:)
      real(8), allocatable :: rprbme(:,:,:,:,:,:),prbasme(:,:,:,:)
      complex(8), allocatable :: rojp(:,:,:),sgpb(:,:,:,:),fouvb(:,:,:,:),
     .                           strx(:,:,:,:),vcoul(:,:)
!     complex(8), allocatable :: sgpp(:,:,:)
      real(8) :: Ecoul,p(3),minnegevl
      real(8), parameter :: tolq = 1d-5, tolevl = 1d-10
      complex(8) :: phasep,img=(0d0,1d0)
      real(8) :: pi=4d0*datan(1d0)

      integer idxdn(n0,nkap0)
      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
      data lsym /'s','p','d','f','g','5','6','7','8','9'/
      data lorb /'p','d','l'/
      data dig /'1','2','3','4','5','6','7','8','9'/

C
CC      integer jfi
CC      integer,allocatable::  ngvecp_(:,:),ngvecc_(:,:)
CC      real(8),allocatable:: evl_(:,:)
CC      complex(8),allocatable:: evec_(:,:),cphi_(:,:,:),pwz_(:,:)
C
C --- Setup ---
      call tcn('makegw')
      call getpr(ipr)
      nproc = mpipid(0)
      procid = mpipid(1)
      master = 0
      rank0 = procid == master  ! true if master process
      gcutb = s_gw%gcutb
      gwcode = s_gw%code
      call sanrg(.true.,gwcode,0,2,'sugw:','gwcode')
      if (.not. associated(s_ham%qsig)) call ptr_ham(s_ham,1,'qsig',1,1,xx)

      xx = 0
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      ngabc = s_lat%nabc
      gmax = s_lat%gmax
      etop = NULLI
      ebot = -etop
      nmcore = 1000*(bitand(s_ctrl%lcd,16)/16)
      npgrp = s_lat%npgrp
      nsgrp = s_lat%nsgrp
      call fftz30(n1,n2,n3,k1,k2,k3)
      stdo = lgunit(1)
      lchk = 1
      nsp  = nglob('nsp')
      nspc = nglob('nspc')    ! 2 for noncollinear case
      nspc3 = nglob('nspc3')  ! 2 if spins coupled perturbatively (lso=3)
      lso =   isw(IAND(s_ham%lncol,4)/=0)
     .    + 2*isw(IAND(s_ham%lncol,lSzLz)/=0)
     .    + 3*isw(IAND(s_ham%lncol,lSzLzp)/=0)
     .    + 4*isw(IAND(s_ham%lncol,lSzLzp0)/=0)
      lrsig= s_ham%lsig
      lwsig = 0
C     lsig = 0
      nqsig = s_ham%nqsig
      compat = 2**0 + 2**2 + 2**3 + 2**4 + 2**5 ! Switches for compatibility with old GW
C ... For reduced hilbert space
      epsovl = s_ham%oveps
      ovncut = s_ham%ovncut
      mnevl = 9999 !for printout, smallest ham dimension
      nlmto = s_ham%nlmto
      tprodbas = 124

C     Number of partial waves for a given l.  For now:
      nphimx = 3
C     Do not shorten vectors
      if (s_ctrl%lfp/2 /= 0) s_ctrl%lfp = s_ctrl%lfp - 2

C ... Count number of atoms : exclude floating orbitals; also get lmxax, maxcor, konfa
      lmxax = -1
      nat = 0
      j = -1
      do  ib = 1, nbas
        is = s_site(ib)%spec
        ipb(ib) = 0
        lmaxa = s_spec(is)%lmxa
        if (lmaxa == -1) cycle
        if (lmaxa > -1) nat = nat + 1
        ipb(ib) = nat
C       ipbi(nat) = ib
        if (j == -1) j = lmaxa
        if (lmaxa /= j) call rx('sugw: lmfgwd requires a fixed lmxa for all sites, sorry')
        lmxax = max(lmxax,lmaxa)
      enddo
C     Get mxcor and konfa
      allocate(konfa(0:lmxax,nat))
      call GWdimpar(1,s_site,s_spec,nbas,lmxax+1,nsp,maxcor,konfa)

      napw = 0
      ldim = s_ham%ldham(1)
      pwmode = s_ham%pwmode; pwemin = s_ham%pwemin; pwemax = s_ham%pwemax; ndham = s_ham%ndham
      allocate(evl(ndham,nsp))
      if (nlmto/=ndham) call rx('makegw: LMTO basis only for now')

C ... Generate ndima,nlindx
      allocate(nlindx(3,0:lmxax,nat))
      nlindx = -1
      ndima = 0
      do  ipqn = 1, nphimx   ! loop over phi, phidot, LO
        do  ib = 1, nbas
          is = s_site(ib)%spec
          lmaxa = s_spec(is)%lmxa
          pnz = s_spec(is)%pz
          lmxax = max(lmxax,lmaxa)
          if (lmaxa > -1) then
            do  l = 0, lmaxa
              npqn = 2
              if (pnz(l+1,1) /= 0) npqn = 3
              if (ipqn <= npqn) then
                nlindx(ipqn,l,ipb(ib)) = ndima
                ndima = ndima + (2*l+1)
              endif
            enddo
          endif
        enddo
      enddo

C --- GW setup files, fpgw format, job 1 ---
C      if (gwcode == 0 .or. gwcode == 2) then

C ... Read QGpsi and QGcou
      if (mod(compat,2) > 0) then
!!       call mk_hamindex(s_ham,s_lat,s_site,s_spec) ! SYMOPS is rewritten
        ifiqg  = fopnx('QGpsi',2,4+1,-1); rewind ifiqg
        read(ifiqg) nqnum, ngpmx, QpGcut_psi, nqbz, nqi,imx,nqibz
      endif

C ... Determine nphimx and nrmx (needed for product basis manipulations)
C     Populate lmxa(:), bas(:,:) in case write header to file gwb
!!     if (mod(compat/2,2) > 0) allocate(lmxa(nat),bas(3,nat))
      iat = 0; nphimx = 0; nrmx = 0; a = 0 ; nr = 0 ; z = 0
      do  i = 1, nbas
        is = s_site(i)%spec
        pnu = s_site(i)%pnu
        pnz = s_site(i)%pz
C       a = s_spec(is)%a
C       nr = s_spec(is)%nr
C       z = s_spec(is)%z
        rmt(i) = s_spec(is)%rmt ! Not needed; avoid compiler complaint
        lmaxa = s_spec(is)%lmxa
        if (lmaxa < 0) cycle
        iat = iat + 1
        if (iat > nat) call rx('bug in sugw')
        call atwf(nmcore,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(i),z,
     .    xx,i1,ncore,konfig,ecore,xx,xx)
        nphimx = max(nphimx,i1)
        nrmx = max(nrmx,s_spec(is)%nr)
!!       if (mod(compat/2,2) == 0) cycle
!!       bas(1:3,iat) = s_site(i)%pos
!!       lmxa(iat) = lmaxa
      enddo

C ... Write header to file gwb
!!     if (mod(compat/2,2) > 0) then
!!       ifi = fopna('gwb',-1,4); rewind ifi
!!       if (rank0) then
!!         ef0 = 1d99 !dummy
!!         write(ifi) nat,nsp,ndima,ndham,alat,qlat,ef0,nqbz,plat,nqnum,nqi !takao 2012 Sep
!!         write(ifi) lmxa(1:nat), bas(1:3,1:nat)
!!       endif
!!       deallocate(lmxa,bas)
!!     endif

C ... Partial waves in augmentation spheres (gwa)
      call info0(30,1,1,' ... Generate radial valence and core partial waves')
!!     if (mod(compat/2,2) > 0 .and. rank0) then
!!       ifi = fopna('gwa',-1,4); rewind ifi
!!     endif
      allocate(gvala(nrmx,0:lmxax,nphimx,nsp,nat))
      allocate(gcora(nrmx,0:lmxax,ncoremx,nsp,nat))
      call dpzero(gvala,size(gvala))
      call dpzero(gcora,size(gcora))
      do  ib = 1, nbas
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        pnz = s_site(ib)%pz
C       p_v0 => s_site(ib)%v0
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        z = s_spec(is)%z
        rmt(ib) = s_spec(is)%rmt
        lmaxa = s_spec(is)%lmxa
        if (lmaxa < 0) cycle
        call atwf(nmcore,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,
     .    s_site(ib)%v0,i1,ncore,konfig,ecore,xx,xx)
        allocate(rofi(nr),rwgt(nr),gcore(nr,2,ncore))
        allocate(gval(nr,2,0:lmaxa,nphimx,nsp))
        call dpzero(gval,nr*2*(1+lmaxa)*nphimx*nsp)
C       Create augmented wave functions for this atom
        call uspecb(0,4,s_spec,is,is,lh,rsml,ehl,i)
        call atwf(nmcore+03,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,
     .    s_site(ib)%v0,nphimx,ncore,konfig,ecore,gcore,gval)
C       Header data for this atom
        b = rmt(ib)/(dexp(a*nr-a)-1d0)
        call radmsh(rmt(ib),a,nr,rofi)
        call radwgt(0,rmt(ib),a,nr,rwgt)

!!  ... Write to disk
!!       if (mod(compat/2,2) > 0 .and. rank0) then
!!         write(ifi) z, nr, a, b, rmt(ib), lmaxa, nsp, ncore
!!         write(ifi) konfig(0:lmaxa)
!!         write(ifi) rofi
!!
C!!        Write orthonormalized valence wave functions for this atom
!!         do  l = 0, lmaxa
!!           do  i = 1, nsp
!!             write(ifi) l,i
!!             write(ifi) gval(1:nr,1,l,1,i)
!!             write(ifi) gval(1:nr,1,l,2,i)
!!             if (konfig(l) >= 10) write(ifi) gval(1:nr,1,l,3,i)
!!           enddo
!!         enddo
C
C!!        Core wave functions for this atom
!!         icore = 0
!!         vshft = 0
!!         do  l = 0, lmaxa
!!           do  isp = 1, nsp
!!             do  konf = l+1, mod(konfig(l),10)-1
!!               icore = icore+1
!!               write(ifi) icore, l, isp, konf, ecore(icore)+vshft
!!               write(ifi) gcore(1:nr,1,icore)
C!!              print *, icore,gcore(nr,1,icore)
!!             enddo
!!           enddo
!!         enddo
!!       endif ! End of write for rank=0

C   ... Keep large component of g
        icore = 0
        do  l = 0, lmaxa
          do  isp = 1, nsp
            do  i = 1, nr
              gvala(i,l,1,isp,ipb(ib)) = gval(i,1,l,1,isp)
              gvala(i,l,2,isp,ipb(ib)) = gval(i,1,l,2,isp)
              if (konfig(l) >= 10) gvala(i,l,3,isp,ipb(ib)) = gval(i,1,l,3,isp)
            enddo
            do  konf = l+1, mod(konfig(l),10)-1
              icore = icore+1
              forall (i=1:nr) gcora(i,l,konf-l,isp,ipb(ib)) = gcore(i,1,icore)
            enddo
          enddo
        enddo
        deallocate(rofi,rwgt,gcore,gval)
      enddo ! Loop over ib
      if (mod(compat/2,2) > 0 .and. rank0) call fclr('gwa',ifi)

C --- Make product basis and coulomb matrix ---
      nnegevl = 0; minnegevl = 0
      nl = lmxax+1
      Ecoul = -1d-4
      lcutmx = 2*(nl-1)

      allocate (nrphiv(0:nl-1,nat),nrphic(0:nl-1,nat),
     .          noccv(0:nl-1,nphimx,nat),noccc(0:nl-1,ncoremx,nat),
     .          nunoccv(0:nl-1,nphimx,nat),nunoccc(0:nl-1,ncoremx,nat),
     .          ncwf(0:nl-1,ncoremx,nat),ncwf2(0:nl-1,ncoremx,nat))
      allocate (tolbas(0:lcutmx),lcuta(nat))
      if (mod(compat/4,2) > 0) then
        call GWinputprodbas(nl,nat,nphimx,ncoremx,lcuta,
     .    nrphiv,nrphic,noccv,nunoccv,noccc,nunoccc,ncwf,ncwf2,tolbas)
      else
        call rx('product basis spec must be read from GWinput for now')
      endif

      if (mod(tprodbas,10) == 1) then
        nocc => noccv
        nunocc => nunoccv
        nrphi => nrphiv
        ndrphi = maxval(nrphiv)
        gtota => gvala
        call rx('sugw need check tprodbas == 1')
      elseif (mod(tprodbas,10) == 2 .or. mod(tprodbas,10) == 4) then
        if (mod(tprodbas,10) == 2) then
C         Modify noccc, noccv according to Takao's prescription.
C         Rethink for future
C         noccc = mirror of nwsc; nunocc = 0;  nunoccv <- (or(nunoccv,occv)
          forall (l=0:nl-1, ic=1:ncoremx, ib=1:nat) noccc(l,ic,ib) = 1-ncwf2(l,ic,ib)
          nunoccc = 0
          forall (l=0:nl-1, ic=1:nphimx, ib=1:nat)  nunoccv(l,ic,ib) = min(nunoccv(l,ic,ib)+noccv(l,ic,ib),1)
        endif
        ndrphi = maxnphi(nrphiv,nrphic,nl,nat) ! Max number of (core+valence) waves in any l
        allocate(nocc(0:nl-1,ndrphi,nat),nunocc(0:nl-1,ndrphi,nat),nrphi(0:nl-1,nat))
        allocate(gtota(nrmx,0:nl-1,ndrphi,nsp,nat))
        nocc = NULLI; nunocc = NULLI
        call dpzero(gtota,size(gtota))
C       Make nocc,nunocc and merge core + valence states into gtota
        call mergephivc(11,nl,nsp,nat,nphimx,ncoremx,ndrphi,nrmx,
     .    nrphiv,nrphic,noccv,nunoccv,noccc,nunoccc,nocc,nunocc,nrphi,
     .    gvala,gcora,gtota)
      endif

C ... Make the product basis B for valence states
      allocate(nrpb(0:lcutmx,nat+1)); call iinit(nrpb,size(nrpb))
      allocate(nrgpb(0:lcutmx,nat+1,2)); call iinit(nrgpb,size(nrgpb))
C     First call only counts number of radial product basis functions for dimensioning (nradpb)
      call pshpr(1)
      i = 100                   ! nradpb for product basis only
      if (cmdopt('--gradpb',8,0,strn)) then
        i = 500                 ! product basis + gradient product basis
      endif
      nradpb = 0
      call rgprodbasa(i,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,nocc,nunocc,tolbas,
     .  nrmx,gtota,nradpb,ndrpbg,nrpb,ntpba,xx,nrgpb,ntgpba,xx)
      call poppr
      mxnpba = maxval(ntpba(1:nat)); ntpb = sum(ntpba(1:nat))
      allocate(rprodb(nrmx,max(1,nradpb(1))),rgprodb(nrmx,max(1,nradpb(2)+nradpb(3))))
C     Second call makes the product basis and possibly its gradient
      i = i + 20
      if (cmdopt('--checkrpb',10,0,strn)) i = i+2000 ! Check completeness
      if (cmdopt('--checkrpb=2',12,0,strn)) i = i+2000 ! Verbose output
      call rgprodbasa(i,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,lcuta,nocc,nunocc,tolbas,
     .  nrmx,gtota,nradpb,ndrpbg,nrpb,ntpba,rprodb,nrgpb,ntgpba,rgprodb)
      ndrpb = ndrpbg(1)
      ntpba(0) = mxnpba         ! prodbas may use ntpba(0), depending on 10s digit mode
C     Debugging printout
      call info0(1,1,0,' debugging printout')
      call info2(1,0,0,'%n:1i',lcutmx*nat+1,nrpb)
      call info2(1,0,1,'%n:1i',nat+1,ntpba)

C ... Make orthogonalized valence partial waves
      if (mod(tprodbas/10,10) > 0) then
        ivec(1:5) = shape(gvala)
        allocate(gvalo(ivec(1),0:ivec(2)-1,ivec(3),ivec(4),ivec(5)))
        i = 1
        if (mod(tprodbas/10,10) > 1) i = 21 ! Kotani conventions to check code
        call gworthphi(i,s_site,s_spec,nl,nsp,nbas,nat,nphimx,nrmx,nrphiv,gvala,gvalo)

        if (mod(tprodbas,10) == 1) then
          gtoto => gvalo
        else
C         Merge core + orthogonalized valence states into gtoto
          allocate(gtoto(nrmx,0:nl-1,ndrphi,nsp,nat))
          call dpzero(gtoto,size(gtoto))
          call mergephivc(10,nl,nsp,nat,nphimx,ncoremx,ndrphi,nrmx,
     .      nrphiv,nrphic,noccv,nunoccv,noccc,nunoccc,nocc,nunocc,nrphi,
     .      gvalo,gcora,gtoto)
        endif

        call nrpb2nbl(11,nl,nat,nrpb,offl,i,j) ! offl and j are dummies; return mxrbls in i
        allocate(rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),i))
        call dpzero(rprbme,size(rprbme))
        call prodbasmea(0,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,nrphi,nrmx,nrpb,rprodb,gtoto,rprbme)
        if (mod(compat/2**4,2) > 0) then
          i = 0; call iopbme(12,1,nat,nl,ndrphi,nrpb,i,rprbme)
        endif

C   ... Make prbasme
        if (mod(tprodbas/100,10) > 0) then
          call mtodim(3,s_site,s_spec,nl,nbas,nat,nrphiv,nrphic,nlnmaug)
          ndlmto = maxval(nlnmaug(1:nat,1))
          is = idalloc('prbasme',allocvb()+2,ndlmto*ndlmto,mxnpba*nat)
          allocate(prbasme(ndlmto,ndlmto,mxnpba,nat))  ! Old GW convention
          call dpzero(prbasme,size(prbasme))
          i = min(9,nsgrp) ! For debugging
          call prodbasa(3,s_lat,s_site,s_spec,i,nl,nsp,nbas,nat,ndrphi,nrphiv,nrphic,
     .      nrpb,ntpba,ndlmto,rprbme,prbasme)
          call prodbaschka(3,s_lat,s_site,s_spec,i,nl,nsp,nbas,nat,ndrphi,nrphiv,nrphic,
     .      ndrpb,nrpb,ntpba,ndlmto,rprbme,prbasme)
        endif
      endif

      if (cmdopt('--stop=prbas',12,0,strn)) call rx0('stop = prbas')

      lgrad = 1
      allocate(rojb(ndrpb,0:lcutmx,0:lgrad,nbas),sgbb(ndrpb,ndrpb,0:lcutmx,nbas),sgbc(ndrpb,0:1,nbas))
      call onsitecoul(s_lat,s_site,s_spec,3+4,lgrad,nbas,lcutmx,ndrpb,nrpb,Ecoul,
     .  nrmx,rprodb,xx,xx,xx,xx,rojb,sgbb,sgbc,xx,xx,xx,xx)

      call info2(1,1,0,' sumcheck sgbc %;12,12F', sum(sgbc),2)

C#ifdef DEBUG
        call info2(1,1,0,' sumcheck rojb %;12,12F  sgbb %;12,12F', sum(rojb),sum(sgbb))
C#endif

C ... For now, read qp from file QIBZ and Q0P and G vectors from file QGcou
      if (mod(compat/8,2) > 0) then
        i = 0; call rdQIBZ(3,nqibz,i,xx,xx)
        nqibze = nqibz + i
        allocate(qibze(3,nqibze),wibze(nqibze))
        call rdQIBZ(30,nqibz,i,qibze,wibze)

C   ... For now, read G vectors from file QGcou
        call rdQGcou(1,nqc,npmbx,QpGcut_cou,xx,xx,xx)
        allocate(ngvecc(3,npmbx,nqc)) ! called ngvecc in old GW code
        allocate(npmb(nqc))           ! called ngc in old GW code
        allocate(qc(3,nqc))
        call rdQGcou(10,nqc,npmbx,QpGcut_cou,qc,npmb,ngvecc)
        if (npmbx < maxval(npmb)) call rx('dimensioning error in file QGcou')
      else
        call rx('qp and G vectors must be read from QIBZ, Q0P, QGCou for now')
      endif

      call info2(10,1,0,' sugw: make matrices and strux for vcoul, E=%g',Ecoul,2)
      s_lat%lmxst = 2*lcutmx
      s_lat%nkdmx = 2*s_lat%nkdmx
      s_lat%nkqmx = 2*s_lat%nkqmx
      call setcg(s_lat,lcutmx,2*lcutmx)
      call lattic(s_lat,s_ctrl,s_site)

C ... For each qp in extended BZ + offset gamma, do
C     grep sumcheck out | grep rojp  |  awk '{print $(NF-1),$NF}' > 2
      nlcutmx = (lcutmx+1)**2
      do  iq = 1, nqibze

        call findqbz(nqc,qc,plat,tolq,qibze(1,iq),iqfbz)
        if (iqfbz < 0) call rx1('sugw cannot find qp=%s,%3d among qp read from QGcou',qibze(1,iq))
        s_lat%igvcc => ngvecc(:,:,iqfbz)
        allocate(rojp(npmb(iqfbz),nlcutmx,nbas))
        allocate(sgpb(npmb(iqfbz),ndrpb,nlcutmx,nbas))
        allocate(fouvb(npmb(iqfbz),ndrpb,nlcutmx,nbas))
C   ... Coulomb matrix
        nmb = npmb(iqfbz)+ntpb
        ntmb = npmbx+ntpb; allocate(vcoul(ntmb,ntmb))
        is = idalloc('sgpb',allocvb()+2,npmb(iqfbz)*ndrpb,nlcutmx*nbas*2)
        is = idalloc('vcoul',allocvb()+2,ntmb,ntmb*2)
        call dpzero(sgpb,2*size(sgpb)); call dpzero(fouvb,2*size(fouvb)); call dpzero(vcoul,2*size(vcoul))
        call onsitecoul(s_lat,s_site,s_spec,30,lgrad,nbas,lcutmx,ndrpb,nrpb,Ecoul,nrmx,rprodb,
     .    qibze(1,iq),ntpb,npmb(iqfbz),ntmb,xx,xx,xx,rojp,sgpb,fouvb,vcoul)

C   ... Structure constants
        i = idalloc('strx',allocvb()+2,nlcutmx*nbas,nlcutmx*nbas)
        allocate(strx(nlcutmx,nat,nlcutmx,nat)); call dpzero(strx,2*size(strx))
        iat = 0
        do  ib = 1, nbas
          is = s_site(ib)%spec
          if (s_spec(is)%lmxa < 0) cycle
          iat = iat + 1
          jat = 0
          do  jb = 1, nbas
            js = s_site(jb)%spec
            if (s_spec(js)%lmxa < 0) cycle
            jat = jat + 1
            p = s_site(jb)%pos - s_site(ib)%pos
            phasep = exp(img*2*pi*sum(qibze(1,iq)*p))
C           For now, strux made to lcutmx for all sites
            nlmi = (lcutmx+1)**2; nlmj = (lcutmx+1)**2
            call strxq(1,Ecoul,qibze(1,iq),p,nlmi,nlmj,nlcutmx*nat,
     .        alat,s_lat%vol,s_lat%awald,s_lat%nkd,s_lat%nkq,s_lat%dlv,s_lat%qlv,
     .        s_lat%cg,s_lat%indxcg,s_lat%jcg,strx(1,iat,1,jat),xx)
C            call zprm0('(1p9e22.12)')
C            call zprm('strx',2,strx(1,iat,1,jat),nlcutmx*nat,nlcutmx,nlcutmx)
          enddo ! jb
        enddo ! ib
        call dscal(2*size(strx),4*pi,strx,1)
C#ifdef DEBUG
        call info2(1,0,0,' sumcheck strx/1e6 iq=%i %2;18,6D',iq,sum(strx)/1d6)
C#endif

        call vcoulq(s_lat,lgrad,qibze(1,iq),ntpb,npmb(iqfbz),ntmb,
     .    nbas,nat,lcutmx,nrpb,ndrpb,qlat,strx,
     .    rojp,rojb,sgbb,sgpb,fouvb,Ecoul,vcoul)

        i = idalloc('strx',allocvb()+4,nlcutmx*nbas,nlcutmx*nbas)
        i = idalloc('sgpb',allocvb()+4,npmb(iqfbz)*ndrpb,nlcutmx*nbas*2)
        deallocate(rojp,sgpb,fouvb,strx)

C   ... Diagonalize
        call tcn('diagonalize')

C       Overlap matrix
        allocate(ppovl(ntmb,ntmb)); call dpzero(ppovl,2*size(ppovl))
        call pwmat(s_lat,s_site,s_spec,1,nbas,xx,xx,
     .    xx,qibze(1,iq),ntmb,npmb(iqfbz),s_lat%igvcc,xx,xx,ppovl(1+ntpb,1+ntpb),xx)
C#ifdef DEBUG
        call info2(1,0,0,' sumcheck ppovl iq=%i %2;18,6D',iq,sum(ppovl))
C#ifdefC DEBUG2
C          call zprm0('(1p9e22.12)')
C          call zprm('vcoul',2,vcoul,ntmb,nmb,nmb)
C          call zprm('ppovl',2,ppovl,ntmb,nmb,nmb)
C#endif
C#endif

C       Overlap matrix for PB part
        forall (i = 1:ntpb) ppovl(i,i) = 1

C       Eigenvectors, eigenvalues of vcoul
        if (allocated(evl)) deallocate(evl)
        allocate(evec(nmb,nmb),ww(max(11,2*nmb)*nmb),evl(nmb,1))
        call dscal(2*size(vcoul),-1d0,vcoul,1)  ! So that evals ordered large to small
        call zhevx(nmb,ntmb,vcoul,ppovl,1,.true.,nmb,9d99,nev,
     .    ww,.false.,evl,nmb,evec)
        call dscal(nmb,-1d0,evl,1)
C       call prmx('evals of coulomb matrix',evl,nmb,nmb,1)

C   ... Reset negative eigenvalues; collect list of evals for printout
        if (allocated(nevls)) deallocate(nevls); allocate(nevls(nmb,1))
        j = 0; k = 0
        do  i = 1, nev
          if (.not. (i > 10 .and. i < nev-5 .and. evl(i,1) > 0)) then
            k = k+1; nevls(k,1) = i; ww(k) = evl(i,1)
          endif
          if (evl(i,1) < tolevl) then
            minnegevl = min(minnegevl,evl(i,1))
            j = j+1; nnegevl = nnegevl+1; evl(i,1) = tolevl
          endif
        enddo

C   ... Printout
        call info5(20,1,0,
     .    ' evals of vcoul for iq = %i of %i  q =%3;9,5D  rank %i%?#n#%-1j  (%i neg evl)##',
     .    iq,nqibze,qibze(1,iq),nmb,j)
        if (iprint() > 45) then
          call arrprt(' state  evl','%,5i  %;7,7F','id',k,0,4,0,'  | ',nevls,ww,3,4,5,6,7,8)
        endif

        call tcx('diagonalize')

C   ... Write to disk, Kotani style
        if (mod(compat/2**5,2) > 0 .and. rank0) then ! Write Vcoud Kotani style
          strn = 'Vcoud.'
          i = 6
          call bin2a('i5',0,0,iq,2,0,20,strn,i)
          ifi = fopng(trim(strn),-1,4)
          call info0(10,0,0,' Write coulomb integrals to file '//trim(strn))
          write(ifi) nmb
          write(ifi) qibze(:,iq)
          write(ifi) evl(1:nmb,1)
          write(ifi) evec
          call fclr(' ',ifi)
        endif

        i = idalloc('vcoul',allocvb()+4,ntmb,ntmb*2)
        deallocate(ppovl,vcoul,ww,evec,evl)

      enddo

      call rx0('finished making bare coulomb integrals')

      end


