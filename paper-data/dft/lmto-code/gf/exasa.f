      subroutine exasa(s_bz,s_ctrl,s_ham,s_lat,s_pot,s_spec,s_site,s_str,s_strn,slabl)
C- Exchange interactions from ASA Green's function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc lshft semsh ef def w dos idtet ipq pdos qp
Ci                 star wtkp wtkb swtk
Co     Stored:     semsh ef def w dos idtet ipq pdos qp star wtkp wtkb
Co                 swtk
Co     Allocated:  qp
Cio    Elts passed:star qp wtkp ipq
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin lcgf lncol lham npl
Ci                 nccomp ldlm tdlm nbasp lgen3 lfp lpgf lrel lasa ldomg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrs lasa lgen3
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 subasi suham
Cio                dlminit suhamz mkpotf makpfz mkpotf mkcpa mkgint
Cio                gfibz asajft hamfbz shoctl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lsig ldham rdvext udiag pwmode pwemin pwemax npwpad
Ci                 lncol neula qss nbf hord nlibu lgen3 bandw ndhrs offH
Co     Stored:     lsig ndham ndofH ldham lmxax npwmin npwpad hord
Co     Allocated:  offH iprmb bdots etrms
Cio    Elts passed:iprmb offH eula magf bdots entropy nprs iaxs hrs
Cio    Passed to:  subasi suham dlminit suhamz mkpotf makpfz mkpotf
Cio                mkcpa mkgint gfibz gf1kp gfg2g asajft gfp0f1 pdglrs
Cio                gfdpp gfg0g
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat qlat nsgrp nabc ag bgv cg cy dlv gv
Ci                 gvq indxcg ips0 istab jcg kv igv igv2 kv2 pos qlv
Ci                 symgr ng tolft gmax
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:symgr istab pos ag gv igv igv2
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 suham sugvec
Cio                sugvec0 mkcpa asajft
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  dlmwt qt ves aamom bxc cp ddpf dddpf ddpfr dpf dpfr
Ci                 gibbs gma gmar grrme mad mxy palp papg pf pfnc pfr
Ci                 pmpol pnu pp ppn pprel pti qc qcorr qnu qpp rhat
Ci                 rnew rhos rhrmx sop thetcl vdif vintr vrmax vshft
Ci                 smpot smrho smrout vmtz
Co     Stored:     gibbs ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf
Co                 dpfr gma gmar grrme mad mxy palp papg pf pfnc pfr
Co                 pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt rhat
Co                 rnew rhos rhrmx sop thetcl vdif vintr vrmax vshft
Co                 smpot smrho smrout nlma nlml
Co     Allocated:  pti cp pf dpf ddpf dddpf pfr dpfr ddpfr gmar papg
Co                 palp gma
Cio    Elts passed: pnu qnu qc qt pp bxc dlmwt vshft pf vrmax gibbs pti
Cio                thetcl pprel dpf ddpf dddpf palp gma pfr ddpfr dpfr
Cio                papg gmar sop cp mxy qcorr ves
Cio    Passed to:  asars iorsa bcast_strx iinit mpibc1 suham dlminit
Cio                suhamz mkpotf makpfz mkpotf mkcpa mkgint gfibz
Cio                gf1kp gfg2g asajft gfp0f1 pdglrs gfdpp gfg0g pasajq
Cio                shoctl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z lmxa name a nr rmt lmxl kmxt p pz lfoca qc idmod
Ci                 rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc coreh coreq idxdn ncomp ngcut orbp
Ci                 iq1 iq2 idu uh jh ncpa nthet iscpa beff hcr
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa
Co                 rhoc idxdn ngcut orbp
Co     Allocated:  rhoc
Cio    Elts passed:rhoc xcpa beff
Cio    Passed to:  asars iorsa bcast_strx getq gtpcor subasi suham
Cio                atfold makidx nscpa showbs sugcut uspecb lmfitmr1
Cio                suldau sudmtu praldm rotycs symdmu getidu dlminit
Cio                dlmwgts cpaidx suhamz mkpotf makpfz mkpotf vorbydl
Cio                asajft shoctl dlmq
C
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  class cpawt dlmcl pnu pos spec force vel pz v0 v1
Ci                 bxc omg omgn domg gc gcu gcorr sfvrtx j0 pdos rho1
Ci                 rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl
Ci                 eqhkl eqkkl sighh sighk sigkk tauhh tauhk taukk pihh
Ci                 pihk pikk sighhx sighkx sigkkx tauhhx tauhkx taukkx
Ci                 pihhx pihkx pikkx thet clabel ncomp norb
Co     Stored:     pnu pos pos0 force vel spec clabel pz bxc cpawt omg
Co                 omgn domg gc gcu gcorr sfvrtx j0 pdos rho1 rho2 rhoc
Co                 rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl
Co                 sighh sighk sigkk tauhh tauhk taukk pihh pihk pikk
Co                 sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Co                 pihkx pikkx thet v0 v1 norb ncomp
Co     Allocated:  v0 v1 cpawt thet j0 bxc omg omgn domg gc gcorr gcu
Co                 sfvrtx
Cio    Elts passed: rhoc rho1 rho2 thet cpawt bxc gcu domg gc gcorr j0
Cio                sfvrtx norb ncomp dlmcl
Cio    Passed to:  asars asars1 iorsa bcast_strx suham setnorb showbs
Cio                pvioeu lmfitmr1 suldau sudmtu praldm rotycs symdmu
Cio                getidu dlminit pkgdlm rdomg cpaidx suhamz mkpotf
Cio                makpfz mkpotf vorbydl mkcpa mkgint asajft pasajq
Cio                pasajf pvgfed
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:s iax alph kaps
Cio    Passed to:  suham rdstrx pp2alp
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  suham str_pack
Ci Inputs
Ci   slabl :vector of species labels
Co Outputs
Co   Generates exchange integrals and writes to file 'jr'
Cs Command-line switches
Cs   --2xrpa    : Double k-mesh when calculate RPA approx to Tc
Cs   --Jq       : Write out Jq
Cs   --amom     :
Cs   --jqraw    :
Cs   --mem      : Use disk to conserve memory
Cs   --nosym    : Not documented
Cs   --nosymdm  : Not documented
Cs   --rs=      : Controls I/O with rst file; see Command-line-options.html
Cs   --sites    : Restrict exchange calc. to selected sites; see doc/gf.html
Cs   --swftitrp :
Cs   --tcqres   :
Cs   --tolJq=   : Sets tolerance for allowed deviation from symmetry in Jr
Cs   --tolimw=  :
Cs   --wmfj     :
Cs   --wsw      :
Cs   -ef=       : Overwrite Fermi level; use with --band
Cs    -long     :
Cs   --long     : Write spin wave energies with extra digits
Cs    -qp=fn    :
Cs   --qp=fn    : Read qp for SW energies from file fn
Cs    -rcut=#   :
Cs   --rcut=#   : Truncate J(R-R') to radius #
Cs   --wrsj     : Write r.s. exchange interactions to disk, for use by other programs
Cs   --vshft    : Read and apply constant potential shifts from vshft file
Cl Local variables
Cl   leres :1 => write exchange integrals (actually integrands)
Cl         :     for every energy iz>=leres
Cl   lwj00 :(switches for writing to rsj file)
Cl         :1 write J00 into rsj file
Cl         :2 write moments into rsj file
Cl         :  (default is to write moment of unit amplitude)
Cl         :4 scale rsj by sign of moments;
Cl         :  also write abs(moments) into rsj file
Cl         :1 may be taken in combination with 2 or 4
Cl   lcgf  :See Remarks
Cl  mxcomp :max number of components associated with any site
Cr Remarks
Cr   Branches of execution.  Branches are mutually exclusive.
Cr   ctrl->lcgf governs which branch is taken.
Cr    10: Transverse exchange interactions J(q) from MST, LW approx
Cr    11: Read J(q) from disk and generate derivative properties
Cr    12: Make G(R,R') (not implemented)
Cr    13: Transverse exchange interactions J(q) from MST (not implemented)
Cr    14: Longitudinal exchange interactions J(q) from MST, LW approx
Cr    20: Transverse chi+- from ASA GF
Cr    21: Read chi(q) from disk and print derivative properties
Cr    24: chi++ from ASA GF
Cr    25: chi-- from ASA GF
Cr    26: transverse susceptibility from G*G
Cr    27: attempt at non-collinear exchange
Cu Updates
Cu   25 Jun 16 (Vishina) extended to relativistic case, working version
Cu   02 Sep 15 (Vishina) extended to relativistic case, first attempt
Cu   08 May 13 Eliminate s_array
Cu   22 Feb 13 DLM case: checks representation dependence of omega file
Cu   17 Jan 13 First attempt at extending exchange calculations to CPA
Cu   31 Dec 12 transverse J can be calculated from given g, rather
Cu             than G (followed by scaling to bare g repsn)
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   25 Oct 11 Started migration to f90 structures
Cu   30 Jun 11 Some improvements to the ferrimagnetic case
Cu   15 Oct 10 Read strux through rdstrx
Cu   09 Aug 10 Proper implementation of MF Tc for many atoms/cell.
Cu   02 Aug 10 Tyablikov formula for Tc, 1 atom/cell
Cu    2 Jul 10 Read pp shifts from vext file
Cu   24 Jun 09 Exclude translation vectors sans full star in printout
Cu   13 Jan 09 First cut at longitudinal response function
Cu   16 Nov 07 New LDA+U (diagonal-only)
Cu   20 Feb 07 Spin waves for 2 atom AFM case
Cu   02 Aug 05 --jq and --jqraw switches
Cu   29 May 03 Correct implementation of --wrsj:scl=...
Cu   12 Apr 03 Better treatment of rcut
Cu   22 Nov 02 Added ability to write r.s. H to file ifi
Cu   21 Dec 01 Potential functions generated by new mkpotf.f
Cu   05 Jun 01 Added generation of L-resolved J_0
Cu   27 Feb 01 New contour where J(q) may be written for each energy
Cu   15 Mar 00 Extended to include iwaves in generation of gf for J
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*8 slabl(*)
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_str)::   s_str
      type(str_strn)::  s_strn(*)
C ... Local structures
      type(str_cpasite),pointer:: s_cpasitei(:),s_cpasitej(:)
C ... Dynamically allocated arrays
      complex(8), pointer :: bgv(:),cJ(:),csym(:),dmatu(:),
     .  gRR(:,:,:,:,:),gii(:),gij(:),gji(:),h(:),hb(:),vorb(:),
     .  wz(:),zp(:),JqRRz(:,:,:)
C     complex(8),pointer :: h2(:),sil(:,:),sii(:,:),gwk(:,:)
      real(8), pointer :: Jq2(:),JqRR(:,:,:,:,:,:),Jr(:),Om(:),Omb(:),
     .  aamom(:),agT(:),amom(:),dos(:,:),dvldu(:),ez(:),gv(:),
     .  sclj(:),secmf(:,:),wk(:),z(:),Jsum(:,:,:),JsCPA(:,:),qtotr(:,:),dqt(:,:)
!A    real(8), pointer :: Jqtest(:,:,:,:,:)
      complex(8),pointer :: sumRL(:,:,:,:,:,:),sumR(:,:,:,:,:,:)
!A    complex(8),pointer :: sumRLt(:,:,:,:),sumRt(:,:,:,:)
      integer, pointer :: idu(:),ilist(:),ips0(:),jlist(:),kv(:),
     .  lldau(:),lmx(:),sid(:),nrc(:),ibcomp(:,:),ipdc(:)
C ... Local parameters
      character outs*512,dc*1,out2*120,jfilnm*72,strn*160
C     character gfopts*120
      double precision avw,alat,plat(3,3),kboltz,chitot
      logical T,F,llshft(3),lso
      integer nkap0,n0H,LGAM,xycomp
      parameter (LGAM=128,nkap0=4,n0H=5)
      integer i,i1,ib,ib2,
     .  idim,ifi,ifirsj,ifib,ifj,ifj2,iq,isp,j,jb,jb2,k,leres,
     .  ld2,ldim,ldimx,lcgf,lham,lidim,lrel,lihdim,lwrite,
     .  ldham(16),lncol,lrsig,nbas,nbcmpi,nbcmpj,nblk,nsitei,nsitej,
     .  nclass,nclasp,nclspd,nbmx,npl,ng,ngdim,ngmx,ni,nj,nl,
     .  nlb,nlspc,nlspcd,nlst,nsp,nspec,nptab,offi,offj,rdm,stdo,pvgfe7,
     .  j1,j2,parg,iorsj,irs(5),lwj00,mode,irep,nbig,
     .  nsgrp,isw,nomg,str_pack,icomp,icmp,iib,jcmp,jib,mxcomp,jcomp,
     .  offai,nlma,nangl,nspc,lidimx,lhdimx,nlmu
      integer mpipid,procid,master
      integer nkxyz(3),is(3),ifac(3),lshft(3),nk1,nk2,nk3,nb1,nb2,nb3,
     .  nq,nzp,iz,nkfbz,nkp,k1,k2,k3,offr,offc
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (ldimx,ldham(5))
      equivalence (nk1,nkxyz(1)), (nk2,nkxyz(2)), (nk3,nkxyz(3))
      logical bittst,cmdopt,lrdqp,ltmp,a2bin,lgetTc,lweiss
C     integer oclp,oclssl,oiaxb
C     kboltz is in units eV K^-1
      integer, parameter :: NULLI=-99999
      parameter (T=.true.,F=.false.,nbmx=1024,kboltz=8.617d-5)
      double precision xx,rcutj,sclrsj,dq,tolrsj,tolJq,fac,facJ0
      double precision dval,pi,tdlm,vconst
      double precision J0(3),J0i(3),efermi(2),gmax,q(3),q1(3),
     .  q2(3),qb(3,3),qbb(3,3),qlat(3,3),rb(3,3),rbb(3,3),tau(3),wint(2),wq(2),wzi(2),zpi(2)
      double precision chioo,chiss,chiso,lttau,deltaG,chiii,DeltaC,avtau,j0cpa,lttau2,tau0,j00
      equivalence (nspc,ldham(4)),(lidimx,ldham(6)),(lhdimx,ldham(7))
C     double precision vconst(3)
C ... For DLM
      integer nccomp,ldlm
C ... For LDA+U
      integer nlibu,lmaxu,ludiag
C     integer nlibu,lmaxu
C     double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     integer odvldu
      procedure(integer) :: bitand,nglob,iprint,lmfitmr1,iprt,fopn,fopna,fopnn,fopno,fopnx,fxst,fhndl

      integer MNC,MNB,MNT,MLB,MPL,MZP,MCD,lio,iogfrs
      parameter (MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64)
C#ifdefC CHKGFJ
CC ... Check free-electron GF
C      integer ocg
C      double precision zkin(2)
C#endif

C --- Setup ---
      stdo = nglob('stdo')
      nbas = s_ctrl%nbas
C     nbasp = nbas
      nclass = s_ctrl%nclass
      nangl = s_ctrl%nccomp
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      lcgf = s_ctrl%lcgf
      lncol = s_ctrl%lncol
      lham = s_ctrl%lham
      npl = s_ctrl%npl
      nccomp = s_ctrl%nccomp
      nclasp = s_ctrl%nclasp
      nclspd = nclasp + nccomp
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      nsgrp = s_lat%nsgrp
      nkp = s_bz%nkp
      nkxyz = s_bz%nkabc
      lshft = s_bz%lshft
      ldlm = s_ctrl%ldlm
      allocate(ipdc(nbas))
      call sitepack(s_site,1,nbas,'dlmcl',1,ipdc,xx)
      lrel = mod(s_ctrl%lrel,10)
      lso = bittst(lncol,4)

      procid = mpipid(1)
      master = 0

      nlspc = nl*nsp*max(nclasp,nspec)
      nlspcd = nl*nsp*max(nclspd,nspec)

      lrsig = s_ham%lsig
      if (lrsig /= 0 .and. fxst('sigm') /= 1) then
        call info0(2,1,0,' EXASA (warning): no sigm file found ... DFT calculation only')
        s_ham%lsig = 0
      endif
      lrsig = s_ham%lsig
      if (lrsig /= 0) then
        call rx('exasa not ready for sigm')
      endif

C ... Input from restart file
      if (cmdopt('--rs=',5,0,outs)) then
      irs(1) = IAND(s_ctrl%lrs,7)
      if (irs(1) > 0) then
        ifi = fopna('rsta',-1,0)
        call asars(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .    s_pot%pnu,s_pot%qnu,.false.,ifi)
C       call clsprp(1,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz)
        call fclr('rsta',ifi)
C       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
C       call rx('done')
      endif
      endif

C     200's bit of lbloch makes -(S-P)
C      lbloch = 200
CC     Use spherical harmonics when making Bloch sum
C      if (bittst(lham,256)) lbloch = 1000
      if (nbas > nbmx) call rx('EXASA: increase nbmx')
C ... Class-dependent arrays
      allocate(z(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'z',1,xx,z)
      allocate(lmx(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'lmxa',1,lmx,xx)
      pi = 4*datan(1d0)

      call info5(20,0,0,'%x%N EXASA: '//
     .  '%?#(n==128)#(gamma repsn) ##%-1j'//
     .  '%?#(n==640)#(gamma-bar repsn) ##'//
     .  '%?#n#%b%b through S^gam) ##'//
     .  '%?;n==10;make transverse exchange J (LW);;%-1j'//
     .  '%?;n==11;print properties calculated from J;;%-1j'//
     .  '%?;n==14;make longitudinal exchange J (LW);;%-1j'//
     .  '%?;n==20;make transverse susceptibility chi;;%-1j'//
     .  '%?;n==21;print properties calculated from chi;;%-1j'//
     .  '%?;n==24;make longitudinal susceptibility chi++,chi--;;%-1j'//
     .  '%?;n==12;make G(R,R'');;',bitand(lham,128)+bitand(s_ctrl%lasa,512),bitand(s_ctrl%lasa,1024),lcgf,0,0)

      call fftz30(nk1,nk2,nk3,k1,k2,k3)
      if (nk1 /= k1 .or. nk2 /= k2 .or. nk3 /= k3)
     .  call rx('exasa: not ready for FFT w/ dimensions ne no div.')

C --- Check that relevant parameters are available ---
      if (lcgf /= 11) then
      do  i  = 1, nclasp+nangl
        if (ldlm > 0 .and. i <= nclasp) then ! skip parent classes
          if (s_ctrl%ncomp(i) >= 2) cycle
        endif
        if (iand(s_ctrl%initc(i),2) == 0) call rx('exasa: missing potential parameters')
        if (lrel == 1 .and. lso) then
          if (iand(s_ctrl%initc(i),4) == 0) call rx('exasa: missing spin orbit coupling parameters')
        endif
      enddo
      endif

C --- Read socscl file for spin-orbit scaling parameters ---
      if (lcgf /= 11) then
      call ptr_pot(s_pot,8+1,'socscl',nclasp+nangl,0,xx)
      if (lso) then
        call dvset(s_pot%socscl,1,size(s_pot%socscl),1d0)
        if (procid == master) then
        if (fxst('socscl') == 1) then
          ifi = fopn('socscl')
          rewind ifi
          call ioextf(2,'SO scaling from socscl',nclasp+nangl,1,s_pot%socscl,s_pot%socscl,1,1,ifi)
          call fclr('socscl',ifi)
          write(stdo,901) ! (s_pot%socscl(i),i=1,nclasp+nangl)
  901     format(' SOC coefficients scaled by: ':/16F6.2)
          call arrprt(' Class scale','%,5i%;7,3D','Id',nclasp+nangl,stdo,7,0,' | ',
     .      xx,s_pot%socscl,xx,xx,xx,xx,xx,xx)
        endif
        endif
        call mpibc1(s_pot%socscl,nclasp+nangl,4,.false.,'gfasa','socscl')
      endif
      endif

C ... Get valence magnetic moments
      allocate(amom(nclspd)); call dpzero(amom,nclspd)
      allocate(aamom(nbas)); call dpzero(aamom,nbas)
      if (.not. cmdopt('--amom',6,0,outs)) then
        call getq(nsp,nl,lmx,nclspd,z,s_pot%pnu,s_pot%qnu,
     .    s_ctrl%ics,s_spec,s_pot%qc,s_pot%qt,amom)
      else
        call getqsw(outs,nclspd,nbas,ipdc,amom,aamom)
      endif
      call cpvprm(0,1,nbas,ipdc,amom,aamom)

C ... Other initialization
C      nspc = 1
C      if (lncol /= 0) nspc = 2
C ... Get v_mtz ... should we shift pp's by del-ves?
C      call pshpr(1)
C      call supot(1,s_ctrl,s_lat,s_pot)
C      call asamad(sspec,s_ctrl,s_pot,s_lat,s_spec,
C     .  0,s_pot%pnu,s_pot%qnu,0d0,s_pot%ves,emad,trumad,vmtz,etrms)
C      call poppr

C ... Read Ef & v_mtz from disk
C      stop 'this should be superseded by call to getef'
C      if (lcgf == 10 .or. lcgf == 11 .or. lcgf == 13) then
C      nfilem = fopn('MOMS')
C      if (iprint() >= 30) print *
C      efermi(1) = 0
C      call iomomq(nfilem,31,nl,nsp,nspc,i,ldim,0,j,ldim,
C     .  nchan,nchan,nevmx,xx,xx,xx,xx,efermi,vmtz)
C      call fclose(nfilem)
C      call awrit2(' ef= %,6d   vmtz= %,6d',outs,80,stdo,
C     .  efermi,vmtz)
C      eferm0 = efermi(1)
C      endif

C ... Command-line override of fermi energy
      efermi(1) = s_bz%semsh(4)
      i = 4
      if (cmdopt('-ef=',i,0,outs)) then
        if (.not. a2bin(outs,efermi,4,0,' ',i,-1)) call
     .    rxs2('EXASA: failed to parse "',outs(1:30),'%a"')
        xx = efermi(1) - s_bz%semsh(4)
        call info2(10,1,1,' Override file fermi level, use ef= %,6d, '//
     .    'delef= %,6d',efermi,xx)
        i = mod(nint(s_bz%semsh(2)),100)
        if (i /= 2) then
          s_bz%semsh(3) = s_bz%semsh(3) + xx
          s_bz%semsh(4) = s_bz%semsh(4) + xx
        endif
      endif

C ... If pp's were generated by lmgf, vshft is built into the pp's.
C     If pp's were generated by lm it is not.
C     In the latter case this block should be called.
      if (cmdopt('--vshft',7,0,strn)) then
        if (procid == master) then
          if (fxst('vshft') == 0) call rx('no file vshft')
          ifi = fopn('vshft')
          xx = NULLI
          call iovshf(nbas,1,'read file vshft: ef=%d  vconst=%d',xx,xx,vconst,s_pot%vshft,ifi)
C         s_pot%vconst = vconst  ! Does this do anything?
          call fclr('VSHFT',-1)
          allocate(wk(nclspd))
          wk = vconst
C         call shoctl(s_ctrl,s_spec,s_pot,101,stdo)
          call shftppr(lrel,nclspd,nl,nsp,s_pot%pp,s_pot%pprel,wk,wk,F,T)
          deallocate(wk)
C         call shoctl(s_ctrl,s_spec,s_pot,101,stdo)
        endif

        if (cmdopt('--fileef',8,0,outs) .and. xx /= NULLI) then
          xx = xx - s_bz%semsh(4)
          call info2(10,1,1,' Use file ef= %,6d, delef= %,6d',xx+s_bz%semsh(4),xx)
          s_bz%semsh(3) = s_bz%semsh(3) + xx
          s_bz%semsh(4) = s_bz%semsh(4) + xx
        endif

      endif

C --- Set up hamiltonian; get alpha, salpa ---
      call subasi(s_ctrl,s_spec,s_ham)
      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,
     .  s_strn)
      ldham = s_ham%ldham
C     offH = s_ham%ooffH

C --- Remake empirical shifts in PPARS ---
      if (s_ham%rdvext /= 0) then
        if (nccomp > 0) call rx('fitting not ready for CPA')
C       No shift unless file present
        if (fxst('vext') == 1 .and. procid == master) then
          ifi = fopn('vext')
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
          call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp)
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
          i = s_ham%rdvext
          i = lmfitmr1(10*i+4,s_ctrl,s_site,s_spec,nl,nbas,nsp,nclasp,
     .      z,s_ham%iprmb,ldham,0,0,xx,s_pot%pp,xx)
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
C         The exasa code requires pp's in alpha repsn
          call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp)
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
C         stop
        else
          call info0(30,1,0,' EXASA (warning) missing file '//
     .      '"vext" required for HAM_RDVEXT ... skipping read')
        endif
C       call mpibc1(sw,1,1,.false.,'exasa','sw_vext')
C       if (sw) call mpibc1(s_pot%pp,6*nlspc,4,mlog,'exasa','shifted pp')
        call mpibc1(s_pot%pp,6*nlspcd,4,.false.,'exasa','shifted pp')
      endif

      idim = lidim - ldim
      ld2 = ldimx**2
      nptab = s_str%npr(nbas+1)

C ... Unpack nkxyz
      k = iabs(s_bz%star(1))
      call cmpi123(0,nkxyz(1),nkxyz(2),nkxyz(3),k)
      nkfbz    = nkxyz(1)*nkxyz(2)*nkxyz(3)
      do  8  i = 1, 3
    8 llshft(i) = lshft(i) /= 0
C ... Make is,ifac,qb,qlat,wght1
      call pshpr(1)
      call bzmsh0(plat,llshft,0,nkxyz(1),nkxyz(2),nkxyz(3),is,ifac,rb,
     .  qb)
      call poppr
C     call dinv33(plat,1,qlat,wght1)

C --- LDA+U initialization ---
      allocate(dvldu(lihdim*nsp)); call dpzero(dvldu,lihdim*nsp)
      allocate(lldau(nbas)); call iinit(lldau,nbas)
C     Check for LDA+U ... return nlibu > 0 if any U blocks.
      call suldau(nbas,s_spec,s_site,nlibu,lmaxu,lldau)
      ludiag = s_ham%udiag
      if (nlibu > 0) then
        i = nsp*nlibu*(lmaxu*2+1)**2
        allocate(vorb(i)); call dpzero(vorb,2*i)
        allocate(dmatu(i)); call dpzero(dmatu,2*i)
        allocate(idu(4*nbas))
C       Initialize density matrix for LDA+U
        i = 0
        if (bittst(lham,256)) i = 1
        call sudmtu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,i,lldau,
     .    nsgrp,s_lat%symgr,s_lat%istab,dmatu,vorb)
C       Diagonal part of LDA+U potential
        call getidu(nbas,s_spec,s_site,idu)
        call u2pph(1,nbas,lmaxu,nsp,nlibu,idu,vorb,
     .    s_ham%iprmb,ldim,lihdim,1d0,dvldu)

C        if (nlibu > 0 .and. bitand(lham,128) == 0)
C     .    call rx('LDA+U requires gamma repsn')
        if (ludiag == 0) call info0(20,1,1,' exasa (warning):  '//
     .    'exchange implemented for diagonal LDA+U only.'//
     .    '%N Restart with HAM UDIAG=1.')

      else
        allocate(idu(1))
        allocate(vorb(1))
        allocate(dmatu(1))
c       nullify(vorb,dmatu,idu)
      endif

C --- Generate energy mesh ---
      nzp = nint(s_bz%semsh(1))
      allocate(zp(nzp),wz(nzp))
      call emesh(s_bz%semsh,zp,wz)
      if (s_bz%semsh(2) /= 2) then
        leres = nzp
      else
        leres = 1
        if (iprint() >= 10) then
          write(stdo,369)
  369     format(9x,'Energy contour mode 2 selected.'/
     .           8x,'*Note : This mode generates the energy derivative',
     .              ' of J, dJ/dE'/
     .           9x,'All quantities labeled as J below actually '/
     .           9x,'correspond to energy derivatives, dJ/dE,'/
     .           9x,'and are written for each point on the contour.')
        endif
      endif

C --- Read shfac file for scaling parameters, incl XC field ---
      if (procid == master) then
        if (fxst('shfac') == 1) then
          ifi = fopn('shfac')
          rewind ifi
          call info0(20,0,0,' GFASA:  reading initial shfac file')
C         call iomagf(2,nclass+nangl,1,s_pot%shfac,s_pot%shfac,1,ifi)
          call ioextf(2,'shfac',nclass+nangl,1,s_pot%shfac,s_pot%shfac,3,1,ifi)
          call fclr('shfac',ifi)
        else
          ifi = 0
        endif
        call mpibc1(ifi,1,2,.false.,'lmasa','ifi')
        if (ifi /= 0) call mpibc1(s_pot%shfac,
     .    (nclass+nangl)*3,4,.false.,'lmasa-gf','shfac')
      endif

C --- DLM setup ---
      if (ldlm /= 0) then
        call dlminit(s_ctrl,s_ham,s_pot,s_spec,s_site,s_ctrl%ics,nbas,
     .    nzp)
        tdlm = s_ctrl%tdlm
        lweiss = mod(mod(ldlm/10,10),2) == 1
        if (procid == master) then
          if (lweiss) then
            if (tdlm > 0) then
             call awrit1('%N DLM: p(x) = -%,3d cos(t)',' ',80,stdo,tdlm)
            else
             call awrit1('%N DLM: p(x) = %,3d cos(t)',' ',80,stdo,-tdlm)
            endif
          elseif (tdlm /= 0d0) then
            call awrit1('%N DLM: Spin T = %,0d K',' ',80,stdo,tdlm)
          else
            call awrit0('%N DLM: Infinite spin temperature',' ',80,stdo)
          endif
        endif

C   ... Read/assign Gibbs weights for DLM
C       ogibbs = odlmwt
        s_pot%gibbs => s_pot%dlmwt
        call pkgdlm(nbas,s_site,s_pot%dlmwt,s_pot%bxc)
C   ... Assign class groups for DLM if necessary
C       call dlmgrp2(nclasp,s_ctrl%ncomp,s_ctrl%idcc,w(ogrp2))

C   ... Read interactor from disk; require repsn match
        i = 1
        if (bittst(s_ctrl%lham,LGAM)) then
          i = 2
          if (bittst(s_ctrl%lasa,512) .or. mod(s_ctrl%lrel,10) == 2) i=3
        endif
C        if (i == 2) call rx('DLM requires alpha or spin-avgd gamma, sorry (GAMMA=0 or 2)')
        call initomg(s_ctrl,s_site,nbas,nspc,nzp)
        call rdomg(100*(nspc-1)+10*i,s_site,nbas,nzp)
c       call rdomg(10*i,s_site,nbas,nzp)
C       call mkcpa(s_ctrl,s_site,s_pot,s_ham,s_lat,F,iz)
      endif
      if (lcgf == 26) then
        chiss = 0
        chiso = 0
        chioo = 0
        chitot = 0
        deltaG = 0
        chiii = 0
      endif

C -------------- Calculate exchange J(q) ---------------
      if (lcgf == 10 .or. lcgf == 13 .or. lcgf == 14 .or.
     .    lcgf == 20 .or. lcgf == 24 .or. lcgf == 26 .or. lcgf == 27) then
      if ((lcgf < 20) .or. (lcgf == 26) .or. (lcgf == 27)) then
        ifj = fopna('JR',-1,0)
      else
        ifj = fopna('chiq',-1,0)
        if (lcgf == 24) ifj2 = fopna('chiq2',-1,0)
      endif
      rewind ifj
      call dpzero(s_pot%qnu,3*nlspcd)
      do  i = 1, nclspd
        call qnu2qnur(-1,nl,nsp,q2,q2,s_pot%qnur(1,i))
      enddo

C     Calculate exchange J_ij for all j even if i is restricted
      allocate(s_cpasitej(1)); allocate(jlist(1)); jlist(1)=0
      call cpaidx(101,s_spec,s_site,s_cpasitej,nbas,0,0,s_ctrl%idcc,
     .  nsitej,nbcmpj,mxcomp,xx)
      deallocate(s_cpasitej); allocate(s_cpasitej(nsitej))
      allocate(ibcomp(4,nbcmpj))
      call cpaidx(131,s_spec,s_site,s_cpasitej,nbas,0,0,s_ctrl%idcc,
     .  nsitej,nbcmpj,mxcomp,ibcomp)

C ... Pick up list of sites (components in the CPA case) to calculate.
C     Syntax --sites[:pair]:site-list
      i = 1
      outs = ' '
      if (cmdopt('--sites',7,0,outs)) then
        allocate(ilist(nbas)); nlst=0; call iinit(ilist,nbas)
        i = 8
        dc = outs(i:i)
        k = i
        j = k+1
C       If option 'pair' specified, associate jlist with ilist
C       jlist used in this branch only when writing J to disk
        call nwordg(outs,0,dc//' ',1,j,k)
        if (outs(j:k) == 'pair') then
          deallocate(jlist); jlist => ilist
          i = k+1
        endif
        call baslst(0,11,outs(i:),j,ipdc,nbas,slabl,z,0,' ',xx,
     .    nlst,ilist)
        if (nlst == 0) then
          ilist(1) = 0
          allocate(s_cpasitei(nbas))
        else
          allocate(s_cpasitei(nlst))
        endif
        call cpaidx(111,s_spec,s_site,s_cpasitei,nbas,nlst,ilist,
     .    s_ctrl%idcc,nsitei,nbcmpi,mxcomp,xx)
C     Case no site list specified: use all sites
      else
        ilist => jlist
        nsitei = nsitej
        nbcmpi = nbcmpj
        s_cpasitei => s_cpasitej
      endif
C     mxcomp = max(mxcomp,1)

      if (ldlm /= 0) then
        call info2(20,1,0,' Site list has %i site(s), %i components',
     .    nsitei,nbcmpi)
      endif

      nlmu = nbas*nl*nl*2
      xycomp = 1
      if (lcgf == 27) xycomp = 3

!        allocate(JqRR(nkfbz,nbcmpi,nbcmpj,1,3,3))
!        call dpzero(JqRR,nkfbz*nbcmpi*nbcmpj*3*3)
      if (lcgf < 20 .or. lcgf == 26 .or. lcgf == 27) then
        allocate(JqRR(nkfbz,xycomp,xycomp,nbcmpi,nbcmpj,1))
        call dpzero(JqRR,nkfbz*nbcmpi*nbcmpj*xycomp*xycomp)
!A      allocate(Jqtest(2,nkfbz,nbcmpi,nbcmpj,1))
!A      call dpzero(Jqtest,2*nkfbz*nbcmpi*nbcmpj)
      else
        allocate(JqRR(nkfbz,1,1,nbcmpi,nbcmpj,2))
        call dpzero(JqRR,nkfbz*nbcmpi*nbcmpj*2)
!A      allocate(Jqtest(2,nkfbz,nbcmpi,nbcmpj,2))
!A      call dpzero(Jqtest,2*nkfbz*nbcmpi*nbcmpj*2)
      endif
      allocate(sumRL(11,xycomp,xycomp,0:mxcomp,0:mxcomp,lidim))
      call dpzero(sumRL,22*(1+mxcomp)**2*lidim*xycomp*xycomp)
      allocate(sumR(11,xycomp,xycomp,0:mxcomp,0:mxcomp,nbas))
      call dpzero(sumR,22*(1+mxcomp)**2*nbas*xycomp*xycomp)

!A    allocate(sumRLt(3,0:mxcomp,0:mxcomp,lidim))
!A    call dpzero(sumRLt,6*(1+mxcomp)**2*lidim)
!A    allocate(sumRt(3,0:mxcomp,0:mxcomp,nbas))
!A    call dpzero(sumRt,6*(1+mxcomp)**2*nbas)

C --- Count classes in given site list ---
      allocate(nrc(nclasp)); call iinit(nrc,nclasp)
      do  iib = 1, nsitei
        ib = s_cpasitei(iib)%ib
C       if (s_cpasitei(iib)%ncomp == 0) then
        i = s_site(ib)%class
        nrc(i) = nrc(i)+1
C       endif
      enddo

C --- For each energy, do ---
      do  iz = 1, nzp

        if (leres == 1) then
           call dpzero(JqRR,nkfbz*nbcmpi*nbcmpj*xycomp*xycomp)
        endif
!A      if (leres == 1) call dpzero(Jqtest,2*nkfbz*nbcmpi*nbcmpj)

C   ... Hamiltonian setup for this energy
        call dpscop(zp,zpi,2,2*iz-1,1,1d0)
        call dpscop(wz,wzi,2,2*iz-1,1,1d0)
        call suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,dvldu,s_pot%vshft,zpi)
        if (nlibu > 0) then
          call vorbydl(-2,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,s_pot%pp,vorb)
        endif
        if (ldlm /= 0) then  ! Only relevant for CPA sites
          call mkcpa(s_ctrl,s_site,s_pot,s_ham,s_lat,F,iz)
C         Always use 100s digit mode=1 since s_site%omgn not available
C         1000s digit makes J0
          call mkgint(2113,s_ctrl,s_site,s_pot,s_ham,iz,wzi,0d0,iz == nzp)
        endif

C   ... Make and save on disk for all points in the irr. BZ for this energy one of:
C          Generate              irep       purpose
C       g^ij (alp repsn)          0        exchange  J requires g^bet, bet spin independent
C       g^ij (gamma-bar rep)      0        exchange  J requires g^bet, bet spin independent
C       G^ij                      1        exchange in gamma -- asajft scales G->g(bare)
C       G^ij                      1        susceptibility; see gfp0f1
C       ... these remarks need updating
C       since G->g, AKA T matrix, merely restores g (pvgfe2)
C       if (lcgf < 20 .and. bittst(lham,128) .and. .not. bittst(s_ctrl%lasa,512)) irep = 1
C       g -> G if susceptibility or if gamma and not gamma-bar
        irep = 0
        if (lcgf >= 20 .or. bitand(lham,128)+bitand(s_ctrl%lasa,512) == 128 .and. lcgf /= 27) irep = 1
        if (lcgf == 26) irep = 1
        if (lcgf == 27) irep = 0
        do  isp = 1, nsp
          if (nspc == 2 .and. isp == 2) cycle
          mode = 100*irep + 22                  ! make g, write to disk.  irep=1 => make G
          if (isp == 1) mode = 100*irep + 42  ! also rewind file before starting
          idu = 0
          call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,mode,s_lat%pos,
     .      nlibu,lmaxu,idu,ludiag,vorb,nkp,s_bz%qp,s_bz%wtkp,nk1,nk2,nk3,
     .      s_bz%ipq,plat,s_str%s,s_str%iax,nptab,s_lat%istab,
     .      s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,qb,xx,xx)
        enddo
        if (nlibu > 0) then
          call vorbydl(-1,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,s_pot%pp,vorb)
        endif

C   ... Convert local g to G by energy scaling
C        if (irep == 0) then
C          call pshpr(iprint()+10)
C          i = 0
C          if (lso) i = 2
C          if (lso .and. bittst(s_ctrl%lbxc,4)) i = 12
C          if (lrel == 2) i = 1
C          if (i /= 0) then
C            call gfg2gr(i,s_site,nbas)
C          endif
C          call poppr
C        endif

C   --- For each (site,component), do ---
        avtau = 0
        tau0 = 0
        j0cpa = 0 ; j00 = 0
        icmp = 0
        do  iib = 1, nsitei
        ib = s_cpasitei(iib)%ib
        do  icomp = 0, s_cpasitei(iib)%ncomp
          if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle
          icmp = icmp+1

C          if (iz == nzp) then
C             print*,'ferm'
C          endif
C     ... Accumulate magnetic exchange interaction J(q) for this component
          call asajft(irep,s_ctrl,s_spec,s_site,s_ham,s_pot,s_lat,
     .      s_cpasitej,nsitej,nbcmpj,ib,icomp,icmp,nbcmpi,mxcomp,xycomp,dvldu,
     .      s_ham%offH,s_ham%iprmb,nbas,nrc,s_lat%pos,nkp,s_bz%qp,iz,nzp,
     .      zpi,wzi,nk1,nk2,nk3,s_bz%ipq,s_bz%star,ifac,qb,JqRR,sumRL,
     .      sumR,s_pot%qnu)

C     ... Last energy level in loop: I/O various quantities
          if (iz >= leres) then
C         if (iz >= leres .or. .true.) then

C       ... Print out sum rule quantities
            if (lcgf == 10) then
              offai = s_ham%offH(4,ib)
              nlma = s_ham%offH(4,ib+1) - s_ham%offH(4,ib)

              fac = 1000
              do  jcomp = 1, mxcomp
                sumR(1:11,1,1,icomp,0,ib) = sumR(1:11,1,1,icomp,0,ib) + sumR(1:11,1,1,icomp,jcomp,ib)
!A              sumRt(1:3,icomp,0,ib) = sumRt(1:3,icomp,0,ib) + sumRt(1:3,icomp,jcomp,ib)
              enddo
              do  jcomp = 0, s_cpasitei(iib)%ncomp*0
                k = 20; if (jcomp > 0) k = 40
                call info8(k,max(1-jcomp,0),0,
     .            ' Sum rule: ib %,2i '//
     .            '%?;n;%-1j ic %,2i ;;'//
     .            '%?;n;%-1j jc %,2i ;;'//
     .            '  pTpT = %,4;9D,%,4;9D; '//
     .            '-p(T+ - T-) = %,4;9D,%,4;9D',ib,icomp,jcomp,
     .            fac*dble(sumR(1,1,1,icomp,jcomp,ib)),
     .            fac*dimag(sumR(1,1,1,icomp,jcomp,ib)),
     .           -fac*dble(sumR(2,1,1,icomp,jcomp,ib)),
     .           -fac*dimag(sumR(2,1,1,icomp,jcomp,ib)),0)
                call info8(k,max(1-jcomp,0),0,
     .            ' Sum rule: ib %,2i '//
     .            '%?;n;%-1j ic %,2i ;;'//
     .            '%?;n;%-1j jc %,2i ;;'//
     .            '  lifetime (mRy) = %,4;9D,%,4;9D; '//
     .            'N/A = %,4;9D,%,4;9D',ib,icomp,jcomp,
     .            fac*dble(sumR(7,1,1,icomp,jcomp,ib)),
     .            fac*dimag(sumR(7,1,1,icomp,jcomp,ib)),
     .           0,
     .           0,0)
                 if (dble(sumR(7,1,1,icomp,jcomp,ib)) /= 0) then
                 lttau2 = 4.83776884e-17/dble(sumR(7,1,1,icomp,jcomp,ib)) ! add
!                 lttau2 = -4.83776884e-17/dimag(sumR(7,1,1,icomp,jcomp,ib))
                 print*,'Lifetime = ',lttau2,'Re Sum7',dble(sumR(7,1,1,icomp,jcomp,ib))!'sec =', lttau2*1e15,'fs' ! add
!                 if (s_cpasitei(iib)%ncomp == 0) then
!                   lttau = 4.83776884e-17/dble(sumR(7,1,1,icomp,jcomp,ib))
!                   print*,'Lifetime = ',lttau,'sec =', lttau*1e15,'fs'
!                 endif
                 end if
                 if (s_cpasitei(iib)%ncomp /= 0) then
!                    print*,'Lifetime (pTpT) = ',-4.83776884e-17/s_site(ib)%tau(icomp,1),'sec ='
                    print*,'J_0 (asajft,vert) = ',-fac*dimag(sumR(2,1,1,icomp,jcomp,ib))
     .                        -fac*dimag(sumR(3,1,1,icomp,jcomp,ib))
                    print*,'J_0 (pTpT) = ',fac*sum(s_site(ib)%j0(icomp,1:nlma))
                    j0cpa = j0cpa - fac*(dimag(sumR(2,1,1,icomp,jcomp,ib))
     .                        +dimag(sumR(3,1,1,icomp,jcomp,ib)))*s_site(ib)%cpawt(icomp)
                    avtau = avtau + dble(sumR(7,1,1,icomp,jcomp,ib))*s_site(ib)%cpawt(icomp)
!                    avtau = avtau + abs(dimag(sumR(7,1,1,icomp,jcomp,ib)))*s_site(ib)%cpawt(icomp)
!                    tau0 = tau0 + s_site(ib)%tau(icomp,1)*s_site(ib)%cpawt(icomp)
                    j00 = j00 + sum(s_site(ib)%j0(icomp,1:nlma))*s_site(ib)%cpawt(icomp)
                 endif

                 if (s_cpasitei(iib)%ncomp /= 0 .and. icomp == s_cpasitei(iib)%ncomp) then
                   lttau = abs(4.83776884e-17/avtau)
                   print*,'Final:'
                   print*,'Re Sum 7',avtau
                   print*,'Lifetime (asajft) = ',lttau,'sec =', lttau*1e15,'fs'
!                   print*,'Lifetime (pTpT) = ',-4.83776884e-17/tau0,'sec ='
                   print*,'J0 (asajft) = ',j0cpa, 'mRy'
                   print*,'J0 (pTpT) = ',j00*1000, 'mRy'
                 endif


!A               call info5(k,max(1-jcomp,0),0,
!A    .            ' Sum rule: ib %,2i '//
!A    .            '%?;n;%-1j ic %,2i ;;'//
!A    .            '%?;n;%-1j jc %,2i ;;'//
!A    .            '  lifetime = %,4;9D,%,4;9D; '//
!A    .            'N/A = %,4;9D,%,4;9D',ib,icomp,jcomp,
!A    .            fac*dble(sumR(7,icomp,jcomp,ib)),
!A    .            fac*dimag(sumR(7,icomp,jcomp,ib)))

!Ax                call info8(k,max(1-jcomp,0),0,
!Ax     .            ' O-O,S-O ib %,2i '//
!Ax     .            '%?;n;%-1j ic %,2i ;;'//
!Ax     .            '%?;n;%-1j jc %,2i ;;'//
!Ax     .            '  MTMT = %,4;9D,%,4;9D; '//
!Ax     .            ' MTPT = %,4;9D,%,4;9D',ib,icomp,jcomp,
!Ax     .            fac*dble(sumR(4,icomp,jcomp,ib)),
!Ax     .            fac*dimag(sumR(4,icomp,jcomp,ib)),
!Ax     .           fac*dble(sumR(5,icomp,jcomp,ib)),
!Ax     .           fac*dimag(sumR(5,icomp,jcomp,ib)),0)

!A              if (nspc == 2) then
!A              Sr(1) = fac*dble(sumRt(1,icomp,jcomp,ib))
!A              Sr(2) = fac*dimag(sumRt(1,icomp,jcomp,ib))
!A              Sr(3) = -fac*dble(sumRt(2,icomp,jcomp,ib))
!A              print*,'pTpT1:',Sr(1),'pTpT2',Sr(2),'p(T-T):',Sr(3)
!A              endif

              enddo

              if (mxcomp == 0) then
              if (abs(sumR(1,1,1,icomp,0,ib)+sumR(2,1,1,icomp,0,ib)) > 1d-6)
     .          call info(20,0,0,' warning : sum rule mismatch',0,0)
              endif
            endif

C       ... Print out sum rule quantities
            if (lcgf == 27) then
              offai = s_ham%offH(4,ib)
              nlma = s_ham%offH(4,ib+1) - s_ham%offH(4,ib)

              fac = 1000
              do  jcomp = 1, mxcomp
                sumR(1:11,:,:,icomp,0,ib) = sumR(1:11,:,:,icomp,0,ib) + sumR(1:11,:,:,icomp,jcomp,ib)
!A              sumRt(1:3,icomp,0,ib) = sumRt(1:3,icomp,0,ib) + sumRt(1:3,icomp,jcomp,ib)
              enddo
              do  jcomp = 0, s_cpasitei(iib)%ncomp*0
                k = 20; if (jcomp > 0) k = 40
!                call info8(k,max(1-jcomp,0),0,
!     .            'zz: Sum rule: ib %,2i '//
!     .            '%?;n;%-1j ic %,2i ;;'//
!     .            '%?;n;%-1j jc %,2i ;;'//
!     .            '  pTpT = %,4;9D,%,4;9D; '//
!     .            '-p(T+ - T-) = %,4;9D,%,4;9D',ib,icomp,jcomp,
!     .            fac*dble(sumR(1,3,3,icomp,jcomp,ib)),
!     .            fac*dimag(sumR(1,3,3,icomp,jcomp,ib)),
!     .           -fac*dble(sumR(2,3,3,icomp,jcomp,ib)),
!     .           -fac*dimag(sumR(2,3,3,icomp,jcomp,ib)),0)
!                call info8(k,max(1-jcomp,0),0,
!     .            'zz: Sum rule: ib %,2i '//
!     .            '%?;n;%-1j ic %,2i ;;'//
!     .            '%?;n;%-1j jc %,2i ;;'//
!     .            '  lifetime (mRy) = %,4;9D,%,4;9D; '//
!     .            'N/A = %,4;9D,%,4;9D',ib,icomp,jcomp,
!     .            fac*dble(sumR(7,3,3,icomp,jcomp,ib)),
!     .            fac*dimag(sumR(7,3,3,icomp,jcomp,ib)),
!     .           0,
!     .           0,0)

                 lttau = 4.83776884e-17/dble(sumR(7,3,3,icomp,jcomp,ib))
!                 print*,'Lifetime = ',lttau,'sec =', lttau*1e15,'fs'
                 print*,'x','y','z'
                 print*,'x',fac*(-dimag(sumR(2,1,1,icomp,jcomp,ib))-dimag(sumR(3,1,1,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,1,2,icomp,jcomp,ib))-dimag(sumR(3,1,2,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,1,3,icomp,jcomp,ib))-dimag(sumR(3,1,3,icomp,jcomp,ib)))
                 print*,'y',fac*(-dimag(sumR(2,2,1,icomp,jcomp,ib))-dimag(sumR(3,2,1,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,2,2,icomp,jcomp,ib))-dimag(sumR(3,2,2,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,2,3,icomp,jcomp,ib))-dimag(sumR(3,2,3,icomp,jcomp,ib)))
                 print*,'z',fac*(-dimag(sumR(2,3,1,icomp,jcomp,ib))-dimag(sumR(3,3,1,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,3,2,icomp,jcomp,ib))-dimag(sumR(3,3,2,icomp,jcomp,ib))),
     .        fac*(-dimag(sumR(2,3,3,icomp,jcomp,ib))-dimag(sumR(3,3,3,icomp,jcomp,ib)))

              enddo

!              if (mxcomp == 0) then
!              if (abs(sumR(1,3,3,icomp,0,ib)+sumR(2,3,3,icomp,0,ib)) > 1d-6)
!     .          call info(20,0,0,' warning : sum rule mismatch',0,0)
!              endif
            endif


C       ... Print out chi
            if (lcgf == 26) then
              offai = s_ham%offH(4,ib)
              nlma = s_ham%offH(4,ib+1) - s_ham%offH(4,ib)
              fac = 1
              do  jcomp = 1, mxcomp
                sumR(1:11,1,1,icomp,0,ib) = sumR(1:11,1,1,icomp,0,ib) + sumR(1:11,1,1,icomp,jcomp,ib)
!A              sumRt(1:3,icomp,0,ib) = sumRt(1:3,icomp,0,ib) + sumRt(1:3,icomp,jcomp,ib)
              enddo
              do  jcomp = 0, s_cpasitei(iib)%ncomp*0
                if (s_cpasitei(iib)%ncomp /= 0) fac = s_site(ib)%cpawt(icomp)
                chiss = chiss + dimag(sumR(1,1,1,icomp,jcomp,ib))
                chiso = chiso + dimag(sumR(3,1,1,icomp,jcomp,ib))
                chioo = chioo + dimag(sumR(4,1,1,icomp,jcomp,ib))
C        ... Gilbert damping
                chitot = chitot + fac*dimag(sumR(1,1,1,icomp,jcomp,ib))
!                chitot = (dabs(dimag(sumR(1,1,1,icomp,jcomp,ib)))+dabs(dimag(sumR(3,1,1,icomp,jcomp,ib)))+
!     .                        dabs(dimag(sumR(4,1,1,icomp,jcomp,ib))))*2.374e-6
!                alpha = sqrt(2d0)*1.7e-4*dabs(amom(ib))/8d0/0.25d0/274d0/274d0*(1d0 + 1d0/chitot)
                deltaG = deltaG + dimag(sumR(7,1,1,icomp,jcomp,ib))
                chiii = chiii + dimag(sumR(2,1,1,icomp,jcomp,ib))
!                print*,'Gilbert damping = ',alpha,'; amom = ',amom(ib)
                k = 20; if (jcomp > 0) k = 40
                call info8(k,max(1-jcomp,0),0,
     .            ' Chi: ib %,2i '//
     .            '%?;n;%-1j ic %,2i ;;'//
     .            '%?;n;%-1j jc %,2i ;;'//
     .            '  G_uu*G_dd (SS) = %,5;12D,%,5;12D; '//
     .            ' G<LS>G<LS> (SO) = %,5;12D,%,5;12D',ib,icomp,jcomp,
!     .            fac*dble(sumR(1,icomp,jcomp,ib)),
!     .            fac*dimag(sumR(1,icomp,jcomp,ib)),
!     .            fac*dble(sumR(2,icomp,jcomp,ib)),
!     .            fac*dimag(sumR(2,icomp,jcomp,ib)),0)
     .            dble(sumR(1,1,1,icomp,jcomp,ib)),
     .            dimag(sumR(1,1,1,icomp,jcomp,ib)),
     .            dble(sumR(3,1,1,icomp,jcomp,ib)),
     .            dimag(sumR(3,1,1,icomp,jcomp,ib)),0)

                call info8(k,max(1-jcomp,0),0,
     .            ' Chi: ib %,2i '//
     .            '%?;n;%-1j ic %,2i ;;'//
     .            '%?;n;%-1j jc %,2i ;;'//
     .            '  mi*G*mj*G (OO) = %,5;12D,%,5;12D; '//
     .            ' S*G*L*G = %,5;12D,%,5;12D',ib,icomp,jcomp,
     .            dble(sumR(4,1,1,icomp,jcomp,ib)),
     .            dimag(sumR(4,1,1,icomp,jcomp,ib)),
     .            dble(sumR(5,1,1,icomp,jcomp,ib)),!0,
     .            dimag(sumR(5,1,1,icomp,jcomp,ib)),0)! 0,0)
                call info8(k,max(1-jcomp,0),0,
     .            ' Chi: ib %,2i '//
     .            '%?;n;%-1j ic %,2i ;;'//
     .            '%?;n;%-1j jc %,2i ;;'//
     .            ' G<LS>G<LS> B||z %,5;12D,%,5;12D; '//
     .            ' G<LS>G<LS> B||x %,5;12D,%,5;12D',ib,icomp,jcomp,
     .            dble(sumR(3,1,1,icomp,jcomp,ib)),
     .            dimag(sumR(3,1,1,icomp,jcomp,ib)),
     .            dble(sumR(6,1,1,icomp,jcomp,ib)),!0,
     .            dimag(sumR(6,1,1,icomp,jcomp,ib)),0)! 0,0)

                 DeltaC = s_pot%pp((2*nl-1)*6+2,i)-s_pot%pp((nl-1)*6+2,i)
!                 print*,'C*chi_ss*C =', DeltaC*dimag(sumR(1,1,1,icomp,jcomp,ib))*DeltaC*1000,'mRy'
!                 print*,'C*chi_ii*C =', DeltaC*dimag(sumR(2,1,1,icomp,jcomp,ib))*DeltaC*1000,'mRy'
!                 print*,'C*DeltaG', DeltaC*dimag(sumR(7,1,1,icomp,jcomp,ib))*1000,'mRy'
                 print*,'DeltaG', dimag(sumR(7,1,1,icomp,jcomp,ib)),'Ry'

              enddo
              if (mxcomp == 0) then
              if (abs(sumR(1,1,1,icomp,0,ib)+sumR(2,1,1,icomp,0,ib)) > 1d-6)
     .          call info(20,0,0,' warning : sum rule mismatch',0,0)
              endif
              if (lcgf == 26 .and. iib == nsitei) then
                print*,'Spin and orbital susceptibilities'
                print*,'Chi_ss (GG) = ',chiss,'Ry^-1 = ',chiss*2.374,'10^-6 emu/G*mol'
                print*,'Chi_so (GLzGLz) = ',chiso,'Ry^-1 = ',chiso*2.374,'10^-6 emu/G*mol'
                print*,'Chi_oo (GmGm) = ',chioo,'Ry^-1 = ',chioo*2.374,'10^-6 emu/G*mol'
                print*,'Chi_tot = ',chitot,'Ry^-1 = ',chitot*2.374,'10^-6 emu/G*mol'
!                print*,'Chi_tot = ',dabs(chioo)+dabs(chiss)+dabs(chiso),'Ry^-1',
!     .                        (dabs(chioo)+dabs(chiss)+dabs(chiso))*2.374,'10^-6 emu/G*mol'
!                print*,'G_up - G_down = ',deltaG,'1/sqrt(Ry)'
!                print*,'chi_ii = ',chiii,'Ry^-1'

!                      do  i = 1, nclspd
!                do lorb = 0, (nl-1)
!                         print*,'c_up-c_down for l=',lorb,' i = ',i,'is',s_pot%pp((nl+lorb)*6+2,i)-s_pot%pp(lorb*6+2,i)
!                enddo

!                      enddo

              endif

            endif


C       ... Write to disk
            jcmp = 0
            do  jib = 1, nsitej
            jb = s_cpasitej(jib)%ib
            do  jcomp = 0, s_cpasitej(jib)%ncomp
              if (jcomp == 0 .and. s_cpasitej(jib)%ncomp > 0) cycle
              jcmp = jcmp+1

              if ( lcgf < 20 .or. lcgf == 26 .or. lcgf == 27) then
                i = pvgfe7(JqRR(1,xycomp,xycomp,icmp,jcmp,1),111,ifj,
     .            ib,icomp,jb,jcomp,nk1,nk2,nk3)

C              else
CC               allocate(Jq(nkxyz(1),nkxyz(2),nkxyz(3)))
CC               offi = nkfbz * (ib-1 + nbas*(jb-1))
CC               call dpscop(JqRR,Jq2,nkfbz,1+offi,1,1d0)
CC               call dpscop(JqRR,Jq,nkfbz,1+offi,1,1d0)
C
C                offi = nkfbz * (jb-1 + nbas*(ib-1))
C                call dpscop(JqRR,Jq2,nkfbz,1+offi,1,1d0)
CC               call dpscop(JqRR,Jq,nkfbz,1+offi,1,1d0)
C                sumR(1,1,jb) = dval(Jq2,1)
C                i = pvgfe7(Jq2,21,ifj,ib,jb,nk1,nk2,nk3)
C                offi = nkfbz * (jb-1 + nbas*(ib-1) + nbas**2)
C                call dpscop(JqRR,Jq2,nkfbz,1+offi,1,1d0)
C                sumR(2,1,jb) = dval(Jq2,1)
C                i = pvgfe7(Jq2,21,ifj2,ib,jb,nk1,nk2,nk3)
CC               deallocate(Jq)
              endif
            enddo
            enddo
C           deallocate(Jq2)

            if (icomp > 0) then
              sumR(1:11,:,:,0,0,ib) = sumR(1:11,:,:,0,0,ib) +
     .          sumR(1:11,:,:,icomp,0,ib)*s_site(ib)%cpawt(icomp)
            endif

            if (icmp == nbcmpi .and. leres == 1 .and. lcgf < 20) then
              call pvgfed(s_ham%offH,s_cpasitei,nsitei,s_site,mxcomp,
     .          sumRL,sumR)
              call dpzero(sumRL,22*(1+mxcomp)**2*lidim*xycomp*xycomp) !11*2
            elseif (icmp == nbcmpi .and. leres == 1 .and. lcgf == 26) then
              call pvgfed(s_ham%offH,s_cpasitei,nsitei,s_site,mxcomp,
     .          sumRL,sumR)
              call dpzero(sumRL,22*(1+mxcomp)**2*lidim)
            elseif (icmp == nbcmpi .and. leres == 1 .and. lcgf == 27) then
              call pvgfed(s_ham%offH,s_cpasitei,nsitei,s_site,mxcomp,
     .          sumRL(:,xycomp,xycomp,:,:,:),sumR(:,xycomp,xycomp,:,:,:))
              call dpzero(sumRL,22*(1+mxcomp)**2*lidim*xycomp*xycomp)
            endif
          endif

        enddo  ! Loop over components
        enddo  ! Loop over sites
        if (mod(iz,5) == 1 .or. iz == nzp) call info5(30,0,0,
     .    '%?;n==1;%N;;%-1j done iz=%i  z = (%d,%d)',iz,zpi,zpi(2),0,0)
        call clhamz(s_pot%pf)
        if (fhndl('GFQP') >= 0) call dfclos(fopn('GFQP'))
      enddo ! Loop over energy
      if (leres == 1) return

C ... Write susceptibility to disk
      if ( lcgf >= 20 .and. lcgf /= 26 .and. lcgf /= 27) then
        rewind ifj
        rewind ifj2

C       Write to disk susceptibility for each (ib,icomp) (jb,jcomp) pair
        icmp = 0
        do  iib = 1, nsitei
        ib = s_cpasitei(iib)%ib
        do  icomp = 0, s_cpasitei(iib)%ncomp
          if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle
          icmp = icmp+1

          do  jcomp = 1, mxcomp
            call rx('exasa update sum rule for susceptibility')
            sumR(1:11,1,1,icomp,0,ib) = sumR(1:11,1,1,icomp,0,ib) +
     .                             sumR(1:11,1,1,icomp,jcomp,ib)
          enddo

          jcmp = 0
          do  jib = 1, nsitej
          jb = s_cpasitej(jib)%ib
          do  jcomp = 0, s_cpasitej(jib)%ncomp
            if (jcomp == 0 .and. s_cpasitej(jib)%ncomp > 0) cycle
            jcmp = jcmp+1

            sumR(1,1,1,1,icmp,jcmp) = sumR(1,1,1,1,icmp,jcmp) + JqRR(1,1,1,jb,ib,1)
            i = pvgfe7(JqRR(1,1,1,jcmp,icmp,1),21,ifj,ib,icomp,jb,jcomp,nk1,nk2,nk3)
            sumR(2,1,1,1,icmp,jcmp) = sumR(2,1,1,1,icmp,jcmp) + JqRR(1,1,1,jb,ib,2)
            i = pvgfe7(JqRR(1,1,1,jcmp,icmp,2),21,ifj2,ib,icomp,jb,jcomp,nk1,nk2,nk3)

          enddo
          enddo
        enddo
        enddo
      endif

C ... Printout
      if ( lcgf < 20 .or. lcgf == 26 .or. lcgf == 27) then
!        if (lcgf == 27) print*,'xx,yy,zz:'

        do i = 1, xycomp
         if (lcgf == 27) then
         if (i == 1) then
          print*,'J_xx:'
         elseif (i == 2) then
          print*,'J_yy:'
         else
          print*,'J_zz:'
         endif
         endif
        call pvgfed(s_ham%offH,s_cpasitei,nsitei,s_site,mxcomp,
     .    sumRL(:,i,i,:,:,:),sumR(:,i,i,:,:,:))
        enddo
        call info(20,1,0,' Moments and charges for ef=%d:',efermi,0)
        if (nsitei == nbas)
     .    call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
        call getq(nsp,nl,lmx,nclspd,z,s_pot%pnu,s_pot%qnu,
     .    s_ctrl%ics,s_spec,s_pot%qc,s_pot%qt,amom)
        if (ldlm > 0) then
          call dlmq(nclasp,s_spec,s_ctrl%ics,s_pot%qt,
     .      s_pot%vrmax,s_ctrl%ncomp,s_pot%gibbs)
        endif
        outs = ' Class   dq     amom'
        if (lrel == 2) then
          if (lrel == 2) call qrel2z12(nl,nclspd,s_pot%qnur)
          write(stdo,"(/' Relativistic sphere charges')")
          allocate(qtotr(2,nclspd),dqt(nclspd,2))
          call qtotrel(nl,nclspd,s_pot%qnur,qtotr)
          do  i = 1, nclspd
            dqt(i,1) = qtotr(1,i) - (z(i)-s_pot%qc(i))
            dqt(i,2) = qtotr(2,i)
          enddo
          s_pot%qt(1:nclspd) = dqt(1:nclspd,1)
          call arrprt(outs,'%,4i%,4;9D%,4;8D','Idd',nclspd,0,3,
     .      0,'  | ',xx,dqt,dqt(1,2),xx,xx,xx,xx,xx)
        else
          write(stdo,"(1x)")
          call arrprt(outs,'%,4i%,4;9D%,4;8D','Idd',nclspd,0,3,
     .      0,'  | ',xx,s_pot%qt,amom,xx,xx,xx,xx,xx)
        endif
C       Compute the total charge
        dq = 0
        do  i = 1, nsitei
          ib = s_cpasitei(i)%ib
          j = s_site(ib)%class
          dq = dq + s_pot%qt(j) ! *s_ctrl%nrc(j)
        enddo
        call info(20,0,0,' Sum of sphere charges = %,6;6d',dq,0)


C     chi++ branch
      elseif (lcgf == 24) then
        call info0(10,1,0,' Resolve sum chi++(0), chi--(0) by site : '//
     .    'sum_j chi_ij(q=0)')

        do  jcomp = 1, mxcomp
        write(stdo,654) ' spin 1', sumR(1,xycomp,xycomp,1,1:nbcmpi,jcomp)
        write(stdo,654) ' total:', sum(sumR(1,xycomp,xycomp,1,1:nbcmpi,jcomp))
        write(stdo,654) ' spin 2', sumR(2,xycomp,xycomp,1,1:nbcmpi,jcomp)
        write(stdo,654) ' total:', sum(sumR(2,xycomp,xycomp,1,1:nbcmpi,jcomp))
  654   format(a, 9f12.6)
        enddo

      endif
      endif   ! Generation of J or chi

C -------------- Properties of J(q) ----------------
      if (lcgf == 11) then

      outs = ' Class   dq     amom'
      call arrprt(outs,'%,4i%:-3,4;4d%:-2,4;4d','Idd',nclasp,0,3,
     .  0,'  | ',xx,s_pot%qt,amom,xx,xx,xx,xx,xx)

      lwj00 = 0
      rcutJ = 5d0
      if (cmdopt('-rcut=',6,0,outs).or.cmdopt('--rcut=',7,0,outs)) then
        j = 6
        if (cmdopt('--rcut=',7,0,outs)) j = 7
        if (.not. a2bin(outs,rcutJ,4,0,' ',j,-1)) call
     .    rxs2('EXASA: failed to parse "',outs(1:30),'%a"')
        rcutJ = -rcutJ
      endif
      call info2(20,1,0,' ... Display properties of Jij (GF mode 11):'//
     .  '  rcut = %d*alat',-rcutJ,0)

      ifj = fopna('JR',-1,0)
      rewind ifj

C ... Pick up list of sites to calculate.
C     Syntax --sites[:pair]:site-list
      allocate(ilist(nbas)); call iinit(ilist,nbas)
      i = 1
      outs = ' '
      if (cmdopt('--sites',7,0,outs)) then
        i = 8
        dc = outs(i:i)
        k = i
        j = k+1
        call nwordg(outs,0,dc//' ',1,j,k)
C       If option 'pair' specified, associate jlist with ilist
        if (outs(j:k) == 'pair') then
          jlist => ilist
          i = k+1
        endif
      else
        jlist => ilist
      endif
      call baslst(0,11,outs(i:),j,ipdc,nbas,slabl,z,0,' ',xx,
     .  nlst,ilist)
      allocate(s_cpasitei(1))
      call cpaidx(101,s_spec,s_site,s_cpasitei,nbas,nlst,ilist,
     .  s_ctrl%idcc,nsitei,nbcmpi,mxcomp,xx)
      deallocate(s_cpasitei); allocate(s_cpasitei(nsitei))
      allocate(ibcomp(4,nbcmpi))
      call cpaidx(131,s_spec,s_site,s_cpasitei,nbas,nlst,ilist,
     .  s_ctrl%idcc,nsitei,nbcmpi,mxcomp,ibcomp)
C     (i,j) linked => link other j-associated quantities to i
      if (associated(ilist,jlist)) then
        nsitej = nsitei
        nbcmpj = nbcmpi
        s_cpasitej => s_cpasitei
        lgetTc = .true.
C     Otherwise the j index ranges over all sites
      else
        allocate(s_cpasitej(1)); allocate(jlist(1)); jlist(1)=0
        call cpaidx(101,s_spec,s_site,s_cpasitej,nbas,0,0,s_ctrl%idcc,
     .    nsitej,nbcmpj,mxcomp,xx)
        deallocate(s_cpasitej); allocate(s_cpasitej(nsitej))
        allocate(ibcomp(4,nbcmpj))
        call cpaidx(131,s_spec,s_site,s_cpasitej,nbas,0,0,s_ctrl%idcc,
     .    nsitej,nbcmpj,mxcomp,ibcomp)
        lgetTc = .false.
      endif

      allocate(sclj(nbcmpi)); call dvset(sclj,1,nbcmpi,1d0)
      allocate(Jsum(3,nbcmpi,0:nbcmpj),JsCPA(3,nbas))
      call dpzero(Jsum,3*nbcmpi*(1+nbcmpj)); call dpzero(JsCPA,3*nbas)

C ... Write Jij to file.
C     Syntax: --wrsj[:j00][:amom][:fn=name][:sscl][:scl=#,#,..][:tol=#]
      ifirsj = 0
      tolrsj = 1d-4
      if (cmdopt('--wrsj',6,0,outs) .or.
     .    cmdopt('-wrsj',5,0,outs)) then

        out2 = outs(7:)
        if (outs(1:5) == '-wrsj') out2 = outs(6:)
        jfilnm = 'rsj'

        dc = out2(1:1)
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 0
   50     continue
          j2 = j2+1
          if (out2(j2:j2) == dc) goto 50
          j1 = min(len(out2),j2)
          call nwordg(out2,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
            elseif (out2(j1:j1+2) == 'j00')  then
              lwj00 = 1 + 2*(lwj00/2)
            elseif (out2(j1:j1+3) == 'sscl')  then
              lwj00 = 4 + mod(lwj00,4)
            elseif (out2(j1:j1+3) == 'amom')  then
              if (lwj00 < 4) lwj00 = 2 + mod(lwj00,2)
            elseif (out2(j1:j1+4) == 'gscl=')  then
              if (ilist(1) == 0 .or. jlist(1) == 0) then
                print *, '*** use --wrsj:scl= '//
     .            ' only with --ssite:pair option'
                goto 59
              endif
              j = 0
              i = parg('gscl=',4,out2(j1:),j,j2-j1+1,
     .          ','//dc//' ',2,1,wq,xx)
              call dvset(sclj,1,nbas,xx)
C             call prmx('sclj',sclj,nbas,nbas,1)
            elseif (out2(j1:j1+3) == 'scl=')  then
              if (ilist(1) == 0 .or. jlist(1) == 0) then
                print *, '*** use --wrsj:scl= '//
     .            ' only with --ssite:pair option'
                goto 59
              endif
              j = 0
              allocate(wk(nbas))
              i = parg('scl=',4,out2(j1:),j,j2-j1+1,
     .          ','//dc//' ',2,nlst,wk,sclj)
              deallocate(wk)
              if (i < nlst) then
                if (i > 0) print *, '*** --wrsj:scl= option '//
     .            ' requires scale for',nlst,' sites'
                goto 59
              endif
              if (i <= 0) goto 59
C             call prmx('sclj',sclj,nlst,nlst,1)
            elseif (out2(j1:j1+3) == 'tol=')  then
              j = 0
              i = parg('tol=',4,out2(j1:),j,len(out2(j1:)),
     .          dc//' ',1,1,k,tolrsj)
              if (i <= 0) goto 59
            elseif (out2(j1:j1+2) == 'fn=')  then
              jfilnm = out2(j1+5:j2)
            else
              goto 59
            endif
            goto 50
   59       call rxs2('exasa: failed to parse --wrsj option: "',
     .        out2(1:j2),'%a ..."')
          endif
        endif

C   ... Write header to rsj; create array for effective moments
        call info0(20,0,0,' ... writing real-space Jij')
        ifirsj = fopna(jfilnm,-1,0)
        allocate(wk(nbas))
C       Write moments into rsj file:
C       Write either 1, amom, |amom|, depending on lwj00
        if (lwj00 > 1) then
          do  ib = 1, nbas
            xx = aamom(ib)
            if (lwj00 >= 4) xx = abs(xx)
            wk(ib) = xx
          enddo
        else
          call dvset(wk,1,nbas,1d0)
        endif
        i = iorsj(0,-ifirsj,nbas,alat,plat,wk,xx,xx,xx,xx,xx)
        deallocate(wk)
      endif

C      call awrit2('%N mag moments: %n:1,3;3d',' ',
C     .    100,6,nclasp,amom)

      allocate(JqRRz(nkfbz,nbcmpi,nbcmpj)) ! jqRR is complex in this branch
      allocate(Jr(nkfbz))
      allocate(Jq2(2*nkfbz))
      allocate(h(nkfbz))
C ... Setup for secular matrix to generate TC in MFA
      allocate(secmf(nbcmpi,nbcmpi),ez(nbcmpi))
      call dpzero(secmf,nbcmpi*nbcmpi)

C --- For each pair of sites, do ---
C ... Re-entry point for exchange parameters at multiple energies
      iz = 0
   30 continue
      iz = iz+1
      if (leres == 1) call dpscop(zp,zpi,2,2*iz-1,1,1d0)

      call dpzero(Jsum,3*nbcmpi*(1+nbcmpj))
      icmp = 0
      do  iib = 1, nsitei
      ib = s_cpasitei(iib)%ib
      do  icomp = 0, s_cpasitei(iib)%ncomp
      if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle
      icmp = icmp+1

        jcmp = 0
        do  jib = 1, nsitej
        jb = s_cpasitej(jib)%ib
        do  jcomp = 0, s_cpasitej(jib)%ncomp
        if (jcomp == 0 .and. s_cpasitej(jib)%ncomp > 0) cycle
        jcmp = jcmp+1

C       print *, 'doing ib,icomp,jb,jcomp',ib,icomp,jb,jcomp,jcmp

C   ... Position file pointer to first pair
        if (icmp == 1 .and. jcmp == 1) then
          i = pvgfe7(Jr,2,ifj,ib,0,jb,0,nk1,nk2,nk3)
          if (i == -2)
     .      call rx('exasa: improper or missing file jr')
          if (leres > 1) call info2(30,0,0,
     .      ' exasa: located (%i,%i) pair in file jr',ib,jb)
        endif

C   ... Shift tau = pos(jb)-pos(ib)
        call dpscop(s_lat%pos,tau,3,3*jb-2,1,1d0)
        call dpsadd(tau,s_lat%pos,3,1,3*ib-2,-1d0)

C   ... Read in real-space J for this ib,jb pair; make Jq; save into jqRRz
        i = pvgfe7(Jr,0,ifj,ib,icomp,jb,jcomp,nk1,nk2,nk3)
        if (i < 0) then
          call info8(10,0,0,' No data for ib=%i%?#n# icomp=%i#%j#, '//
     .      'jb=%i%?#n# jcomp=%i## ... skipping',
     .      ib,s_cpasitei(iib)%ncomp,icomp,
     .      jb,s_cpasitej(jib)%ncomp,jcomp,0,0)
          cycle
        endif

C       Write J(q)
C        rewind ifj
C        i = pvgfe7(Jr,11,ifj,ib,icomp,jb,jcomp,nk1,nk2,nk3)
C        call rx('done')

C   ... Scale exchanges for this pair by factor
        sclrsj = dsqrt(sclj(icmp)*sclj(jcmp))
C       Scale rsj by sign of moments
        if (lwj00 >= 4) then
          sclrsj = sclrsj*
     .      dsign(1d0,dval(aamom,ib)) * dsign(1d0,dval(aamom,jb))
        endif
        call dscal(nkfbz,sclrsj,Jr,1)

C   ... Make Jq (for this ib,jb pair) and JqRR (all ib,jb pairs)
        call dpzero(h,nkfbz*2)
        call dcopy(nkfbz,Jr,1,h,2)  ! Jr is real; h is complex
C       call zprm3('J(r)',0,h,nk1,nk2,nk3)
        call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,1)
C       call zprm3('J(q), FFT',0,h,nk1,nk2,nk3)
C       call pvgfeg(plat,s_bz%ipq,nk1,nk2,nk3,h)
C       call pvgfeb(plat,Jr,tau,nk1,nk2,nk3,h)
C       call zprm3('J(q) by brute force',0,h,nk1,nk2,nk3)
        call dpcopy(h,JqRRz(1,icmp,jcmp),1,2*nkfbz,1d0)
        call dpscop(h,Jq2,2*nkfbz,1,1,1d0)
C       call pvgfe9(plat,h,tau,nk1,nk2,nk3,Jq2)
C        call dpscop(Jq2,JqRRz,nkfbz,1,
C       .  1+nkfbz*(ib-1+nbas*(jb-1)),1d0)

C   ... Debugging: convolve Jq with a Gaussian
C        call mkqlat(plat,qlat,xx)
C        fac = (3d0/xx/16/datan(1d0))**(1d0/3d0) ! this is WSR in Q-space
C        xx = fac * .04d0*2
C        call dpscop(Jq2,h,2*nkfbz,1,1,1d0)
C        call pvgfem(1,plat,nk1,nk2,nk3,xx,-1,h,1)
C        call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,-1)
C        call dcopy(nkfbz,h,2,Jr,1)
C        call zprm3('J(r) convolved',0,h,nk1,nk2,nk3)

        ltmp = T
        if (nbas /= 1) then
          outs = ' '
          call awrit2('Display properties for ib=%i, jb=%i?',
     .      outs,80,0,ib,jb)
          call query(outs,0,ltmp)
        endif
        if (.not. ltmp) then
          call awrit2(' ... skipping ib=%i, jb=%i pair',
     .      outs,80,stdo,ib,jb)
          cycle
        endif

C   --- J summed over all connecting vectors for this ib,jb pair ---
C       J0(1) = sum_T J_i,j+T = J_ij(q=0)
C       J0(2) = J_i,i
        J0(1) = JqRRz(1,icmp,jcmp)
        call dpzero(h,nkfbz*2)
        J0(2) = dval(Jr,1)
        if (ib /= jb) J0(2) = 0
        if (lgetTc) then
          secmf(icmp,jcmp) = J0(1) - J0(2)
        endif

        if (s_cpasitei(iib)%ncomp == 0) then
          call info5(20,1,0,' Sum J_ij=%;9,6D  J_00=%;9,6D  J_0=%;9,6D'
     .      //' for ib=%i, jb=%i',J0(1),J0(2),J0(1)-J0(2),ib,jb)
        else
          call info8(20,1,0,' Sum J_ij=%;9,6D  J_00=%;9,6D  J_0=%;9,6D'
     .      //' for ib=%i icomp=%i, jb=%i jcomp=%i',
     .      J0(1),J0(2),J0(1)-J0(2),ib,icomp,jb,jcomp,0)
        endif

C       Jsum(1,ib) = sum_j sum_T J_ij+T = sum_j J_ij(q=0)
C       Jsum(2,ib) = J_ii
        Jsum(1,icmp,jcmp) = Jsum(1,icmp,jcmp) + J0(1)
        if (ib == jb) Jsum(2,icmp,jcmp) = J0(2)

C   --- display J(R) ---
C   ... Optionally suppress symmetry operations
        if (cmdopt('--nosym',7,0,outs)) nsgrp = 1

C   ... Shift tau = pos(jb)-pos(ib)
        call dpscop(s_lat%pos,tau,3,3*jb-2,1,1d0)
        call dpsadd(tau,s_lat%pos,3,1,3*ib-2,-1d0)

C   ... Cutoff in real space; make list of lattice vectors within cutoff
        call pshpr(iprint()-30)
        call gvctof(0,1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx)
        if (cmdopt('--allg',6,0,outs)) then
          gmax = 9999
          ngmx = nk1*nk2*nk3
        endif
        call poppr
        allocate(gv(ngmx*3))
        allocate(kv(ngmx*3))
C   ... Make list for a dimensionless lattice (alat=1)
        call ogvlst(1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx,ng,gv,kv)
        if (ng /= ngmx) call rx('exasa: bug in gvlist')
        allocate(ips0(ng))
        allocate(bgv(ng))
        allocate(agT(3*nsgrp))
        call dpcopy(s_lat%ag,agT,1,3*nsgrp,1d0)
        call shfspc(s_lat%symgr,s_lat%ag,nsgrp,s_lat%pos,jb,plat,agT)
        call srvsym(nsgrp,s_lat%symgr,agT,tau,plat,ng,gv,
     .    -gmax/2/pi,ips0,bgv)
        deallocate(agT)
C   ... Symmetrize J
        call dpzero(h,nkfbz*2)
        call dcopy(nkfbz,Jr,1,h,2)
        allocate(cJ(ng))
        allocate(csym(ng))
        call gvgetf(ng,1,kv,nk1,nk2,nk3,h,cJ)
C   ... debugging
C        print *, 'debugging .. muck up cj'
C        call dpzero(cJ, 2*ng)
C        call dvset(w(ocj),3,3,1d0/0.09375d0*.75d0)
        call gvsym(ng,gv,ips0,bgv,cJ,csym)
C       print *, 'xx',ib,jb,icmp,jcmp,sclrsj,leres
        call pvgfe1(ng,ips0,gv,cJ,csym,tau,rcutj,Jsum(1,icmp,jcmp),ib,
     .    jb,ifirsj,1d0,tolrsj,lwj00)
        deallocate(ips0,bgv)
C   ... Print out J(q)
        if (cmdopt('--Jq',4,0,outs)) then
          call gvputf(ng,1,kv,nk1,nk2,nk3,csym,h)
          if (cmdopt('--jqraw',7,0,outs)) then
            call dpzero(h,nkfbz*2)
            call dcopy(nkfbz,Jr,1,h,2)
          endif
          call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,1)
          call pvgfeg(plat,s_bz%ipq,nk1,nk2,nk3,h)
          call rx('finish Jq branch')
        endif

C   ... If to truncate J_R, ditto for JqRR
        if (cmdopt('-rcut=',6,0,outs).or.cmdopt('--rcut=',7,0,outs)) then
          call gvputf(ng,1,kv,nk1,nk2,nk3,csym,h)
          call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,1)
C         call prm3('J(q) before truncation',0,Jq2,nk1,nk2,nk3)
          call dpcopy(h,JqRRz(1,icmp,jcmp),1,2*nkfbz,1d0)
C         call dpscop(h,Jq2,2*nkfbz,1,1,1d0)
C         just to keep Jq used later.  But for now only valid for ib=jb
C          call pvgfe9(plat,h,tau,nk1,nk2,nk3,Jq2)
C          call dpscop(Jq2,JqRRz,nkfbz,1,
C       .    1+nkfbz*(ib-1+nbas*(jb-1)),1d0)
        endif

C       Note: Jq is only valid for one atom/cell.  Need to fix
        if (nbas /= 1) cycle

CC   --- Get J''(q=0) ---
CC       Scaling to convert from Ry-dimensionless to ev-AA^2
C        facAA = 13.6d3 * (alat*.529d0)**2
CC       cJr = R**2 * cJ
C        ocJr = ocJ
C        call pvgfea(1d0,ng,gv,csym,0d0,w(ocJr))
CC       h = R**2 * J on a uniform mesh
C        call gvputf(ng,1,kv,nk1,nk2,nk3,w(ocJr),h)
CC       call zprm3('R**2 J(R)',0,h,nk1,nk2,nk3)
C        call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,1)
CC       call zprm3('J''''(q)',0,h,nk1,nk2,nk3)
CC ...   Factor of 6 to convert nabla to coff of q**2
C        D1 = -dval(h,1)/6
C        D2 = D1
C        call awrit0('%N ... Extract stiffness from FT of R**2 J(R)',
C       .  ' ',80,stdo)
C        call pvgfe6(1,1,1,amom,D2)
C        call awrit3(' coff to q**2 expansion of J: %,6;6d'//
C       .            '   of w: %;6d Ry-au**2 = %;0d meV-A**2 ',
C       .            ' ',80,stdo,D1,D2*alat**2,D2*facAA)
C        deallocate(gv)
C
CC ---   Make J(q) on a doubled mesh ---
C        call zcopy(nkfbz,Jq2,1,h,1)
CC       call dpzero(h,2*nkfbz)
CC       call dcopy(nkfbz,Jq2,1,h,2)
C
CC        nbig = nkfbz*8
CC        allocate(hb(nbig))
CC        call chgmsh(3,qlat,1,nk1,nk2,nk3,h,nb1,nb2,nb3,hb)
CCC       call zprm3('J(q), fine mesh',0,hb,nb1,nb2,nb3)
CC        call dcopy(nbig,hb,2,w(oJqb),1)
CCC       call prm3('J(q), fine mesh',0,w(oJqb),nb1,nb2,nb3)
CC        deallocate(hb)
C
CC ...   Revert to old mesh
CC        oJqb = oJq2
C         nb1 = nk1
C         nb2 = nk2
C         nb3 = nk3
C         nbig = nkfbz
C
CC        call pshpr(1)
CC        call bzmsh0(plat,llshft,0,nb1,nb2,nb3,is,ifac,rbb,qbb)
CC        call poppr
C
C
CC ---   Lim J(q), q->0 --
C        call zcopy(nkfbz,Jq2,1,h,1)
CC       call dpzero(h,2*nkfbz)
CC       call dcopy(nkfbz,Jq2,1,h,2)
CC       call zprm3('J(q)',0,h,nk1,nk2,nk3)
C        call pvgfe3(plat,nk1,nk2,nk3,h,D1,D2)
C        D1 = D1/(2*pi)**2
C        D2 = D1
C        call awrit0('%N ... Extract stiffness from finite diff of J(q)',
C       .  ' ',80,stdo)
C        call pvgfe6(1,1,1,amom,D2)
C        call awrit3(' coff to q**2 expansion of J: %,6;6d'//
C       .            '   of w: %;6d Ry-au**2 = %;0d meV-A**2 ',
C       .            ' ',80,stdo,D1,D2*alat**2,D2*facAA)
CC ...   Repeat for finer mesh
CC        call dcopy(nkfbz,w(oJqb),1,hb,2)
CC        call pvgfe3(plat,nb1,nb2,nb3,hb,D1,D2)
CC        D1 = D1/(2*pi)**2
CC        D2 = D1
CC        call pvgfe6(1,1,1,amom,D2)
CC        call awrit2(' coff to q**2 expansion of J: %,6;6d'//
CC       .            '   of w: %;6d (2x mesh)',
CC       .            ' ',80,stdo,D1,D2)
CC        deallocate(hb)

        enddo  ! Loop over components jb
        enddo  ! Loop over sites jb
      enddo  ! Loop over components ib
      enddo  ! Loop over sites ib

C --- Printout of total J ---
      if (.true.) then
      if (lcgf == 26) then
      call info2(20,1,0,' Sum chi_ij for all pairs connected to'//
     .    ' listed sites (Ry)%?#n==1#, z=%2:1,6;6d##',leres,zpi)
      else
      call info2(20,1,0,' Sum J_ij for all pairs connected to'//
     .    ' listed sites (mRy)%?#n==1#, z=%2:1,6;6d##',leres,zpi)
      endif
      if (mxcomp > 1) then
        write(stdo,1) ' comp  '
      else
        write(stdo,1) '       '
      endif
      if (lcgf == 26) then
    2 format(' site',a,'Sum chi_ij',7x,'chi_00',4x,'chi_0(r<rcut)',6x,'chi_0',
     .  4x,'2/3 chi_0 (K)')
      else
    1 format(' site',a,'Sum J_ij',7x,'J_00',4x,'J_0(r<rcut)',6x,'J_0',
     .  4x,'2/3 J_0 (K)')
      endif

      J0(1) = 0
      J0(2) = 0
      J0(3) = 0
      facJ0 = 2d0/3d0 * 13.6d0 / kboltz
      if (lcgf == 26) then
       fac = 1
      else
      fac = 1000
      endif

      icmp = 0
      do  iib = 1, nsitei
      ib = s_cpasitei(iib)%ib
      do  icomp = 0, s_cpasitei(iib)%ncomp
      if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) cycle
      icmp = icmp+1

        do  i = 1, 3
          Jsum(i,icmp,0) = 0
        enddo
        do  jcmp = 1, nbcmpj
          if (ibcomp(3,jcmp) /= 0) then
            jb = ibcomp(1,jcmp)
            jcomp = ibcomp(4,jcmp)
            xx = s_site(jb)%cpawt(jcomp)
          else
            xx = 1
          endif
          do  i = 1, 3
            Jsum(i,icmp,0) = Jsum(i,icmp,0) + Jsum(i,icmp,jcmp)*xx
          enddo
        enddo

        if (s_cpasitei(iib)%ncomp /= 0) then
          xx = s_site(ib)%cpawt(icomp)
        else
          xx = 1
        endif
        JsCPA(1:3,ib) = JsCPA(1:3,ib) + Jsum(1:3,icmp,0)*xx

      enddo
      enddo

      icmp = 0
      do  iib = 1, nsitei
      ib = s_cpasitei(iib)%ib
      do  icomp = 0, s_cpasitei(iib)%ncomp
        if (icomp == 0 .and. s_cpasitei(iib)%ncomp > 0) then
          J0i(1:3) = JsCPA(1:3,ib)
        else
          icmp = icmp+1
          J0i(1:3) = Jsum(1:3,icmp,0)
        endif

        if (J0i(1) /= 0 .and. J0i(2) /= 0) then
          i = 10; if (icomp > 0) i = 30
          call info8(i,0,0,'%,4i%?#n#%-1j%,4i#    #'//
     .      '%,3;12D%,3;12D%,3;12D%,3;12D%,1;10D',ib,icomp,
     .        fac*J0i(1),fac*J0i(2),
     .        fac*(J0i(3)-J0i(2)),
     .        fac*(J0i(1)-J0i(2)),
     .        facJ0*(J0i(1)-J0i(2)),0)
        endif
        if (icomp == 0) then
          J0(1) = J0(1) + J0i(1)/nsitei
          J0(2) = J0(2) + J0i(2)/nsitei
          J0(3) = J0(3) + J0i(3)/nsitei
        endif
      enddo
      enddo

        write(stdo,827) fac*J0(1),fac*J0(2),
     .    fac*(J0(3)-J0(2)),
     .    fac*(J0(1)-J0(2)),
     .    facJ0 * (J0(1)-J0(2)),
     .    facJ0 * (J0(3)-J0(2))
  827   format(' -------'/' avg',4x,4f12.3,3x,f7.1,' (',f7.1,' r<rcut)')

      endif

      if (ldlm > 0) call rx0(
     .  'remainder of exchange analysis not implemented for CPA')

C     Output J(q=0) as input for mean-field calculation
      lwrite = 0
      if (cmdopt('--wmfj',6,0,outs)) then
        lwrite = 1
C       If option 'mom' specified, add 10 to lwrite to scale J w/ mom
        i = 7
        dc = outs(i:i)
        k = i
        j = k+1
        call nwordg(outs,0,dc//' ',1,j,k)
        if (outs(j:k) == 'mom') lwrite = lwrite+10
        call info0(20,0,0,' ... writing file jmf for mean-field hamiltonian')
        print *, 'warning: check ilist in pvgfee call'
        ifi = fopna('Jmf',-1,0)
        call pvgfee(nbas,ilist,jlist,nk1,nk2,nk3,lwrite,ipdc,amom,Jsum,JqRRz,ifi)
        call fclose(ifi)
        call rx0('finished writing Jmf')
      endif

C ... Mean-field estimate of Tc
      if (lgetTc) then
        call pvgfej(nbcmpi,secmf,ez)
        call arrprt(' Site  ez','%,4i%:-3,4;4d','id',nbcmpi,0,4,
     .    0,' | ',s_cpasitei(:)%ib,ez,xx,xx,xx,xx,xx,xx)
      endif

C ... Close of loop over energy points
      if (leres == 1 .and. iz < nzp) goto 30

C      if ((nlst == 1 .or. nbas == 1) .and. leres > 1) then
      if (leres > 1) then

C --- Symmetrize JqRR ---
      tolJq = 1e-12
      j = 8
      if (cmdopt('--tolJq=',j,0,outs)) then
      if (.not. a2bin(outs,tolJq,4,0,' ',j,-1)) call
     .  rxs2('EXASA: failed to parse "',outs(1:30),'%a"')
      endif
      call pvgfeo(nk1,nk2,nk3,s_cpasitei,nsitei,tolJq,nbcmpi,JqRRz)

C --- Tc, Tyablikov (RPA) ---
      if (lgetTc) then
        j = str_pack('mix',-2,s_strn,strn)
        if (cmdopt('--2xrpa',7,0,outs)) then
          nbig = nkfbz*8
          allocate(hb(nbig*nbcmpi**2))
          call pvgfen(1,nk1,nk2,nk3,nbcmpi,plat,JqRRz,nb1,nb2,nb3,hb)
          call pvgfek(nb1,nb2,nb3,strn,s_cpasitei,nsitei,nbcmpi,hb,ez)
          deallocate(hb)
        else
          call pvgfek(nk1,nk2,nk3,strn,s_cpasitei,nsitei,nbcmpi,JqRRz,ez)
      endif
      endif

C --- Spin wave spectra ---
      allocate(Om(nkfbz*nbcmpi)); call dpzero(Om,nkfbz*nbcmpi)
      call pvgfe8(nk1,nk2,nk3,plat,ipdc,amom,
     .  s_cpasitei,nsitei,nbcmpi,JqRRz,Om,nomg)

C --- Spin stiffness ---
C ... Cutoff in real space; make list of lattice vectors within cutoff
      if (nomg > 0) then
        call dpzero(tau,3)
        call pshpr(iprint()-30)
        call gvctof(0,1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx)
        if (nsgrp == 1) then
          gmax = 9999
          ngmx = nk1*nk2*nk3
        endif
        call poppr
        allocate(gv(ngmx*3))
        allocate(kv(ngmx*3))
C   ... Make list for a dimensionless lattice (alat=1)
        call ogvlst(1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx,ng,gv,kv)
        if (ng /= ngmx) call rx('exasa: bug in gvlist')
        allocate(Jq2(2*nkfbz))
        call pvgfef(nk1,nk2,nk3,nbas,alat,Om,ng,gv,kv,Jq2)
        deallocate(Jq2,gv,kv)
      endif

C --- Make omega(q) on a doubled mesh ---
      if (nomg > 0) then
        allocate(h(nkfbz))
        call dpzero(h,2*nkfbz)
        call dcopy(nkfbz,Om,1,h,2)
        nbig = nkfbz*8
        allocate(Omb(nbig))
        allocate(hb(nbig))
        call pvgfen(1,nk1,nk2,nk3,1,plat,h,nb1,nb2,nb3,hb)
C       call chgmsh(3,qlat,1,nk1,nk2,nk3,nk1,nk2,nk3,h,nb1,nb2,nb3,
C    .    2*nk1,2*nk2,2*nk3,hb)
C       call zprm3('J(q), fine mesh',0,hb,nb1,nb2,nb3)
        call dcopy(nbig,hb,2,Omb,1)
C       call prm3('omega(q), fine mesh',0,Omb,nb1,nb2,nb3)
        deallocate(hb)
      endif

C --- Tyablikov formula for Tc, one atom/cell ---
      if (nbas == 1) then

C       See Kudrnovsky PRB 64, 174402, Eq. 10
        fac = 2d0/3d0 * 13.6d0 * dabs(dval(amom,1)/4)
        outs = ' '
        ltmp = cmdopt('--tcqres',8,0,outs)
        call pvgfe4(plat,nk1,nk2,nk3,Om,fac,outs,wint)

        call info2(10,1,0,'  MF vs RPA estimate for Tc'//
     .    '%?;n;.  Tc k-resolved (file tcqres);;',
     .    isw(ltmp),0)
        write(stdo,567)
  567   format('  mesh  2/3 1/N sum_q [J(q)-J(0)]',5x,
     .    '2/3 [1/N sum_q [J(q)-J(0)]^-1]^-1'/
     .    15x,'meV     K',21x,'meV     K')

C        wint(1) = 2d0/3d0 * 13.6d0 * wint(1) * dval(amom,1)/4
C        wint(2) = 2d0/3d0 * 13.6d0 * wint(2) * dval(amom,1)/4
        write(*,568) '  ',
     .    wint(2)*1d3,wint(2)/kboltz,wint(1)*1d3,wint(1)/kboltz
  568   format(3x,a,6x,f8.3,f8.1,14x,f8.3,f8.1)

        call pvgfe4(plat,nb1,nb2,nb3,Omb,fac,outs,wint)
C        wint(1) = 2d0/3d0 * 13.6d0 * wint(1) * dval(amom,1)/4
C        wint(2) = 2d0/3d0 * 13.6d0 * wint(2) * dval(amom,1)/4
        write(*,568) '2x',
     .    wint(2)*1d3,wint(2)/kboltz,wint(1)*1d3,wint(1)/kboltz

C        call info2(10,1,0,' volq [integral 1/w(q) d^q]^-1 : '//
C     .    '%,3;3d (meV) = %,1;1d (K)',wint,wint(1)/kboltz)
C        call info2(10,0,0,' volq [integral 1/w(q) d^q]^-1 : '//
C     .    '%,3;3d (meV) = %,1;1d (K) (2x mesh)',wint,wint(1)/kboltz)
      endif

C --- omega(q) along selected symmetry lines ---
      if (cmdopt('--swftitrp',10,0,outs)) then
        call info0(10,1,0,' ... Make w(q) along lines in syml file'//
     .    ' (FT interp)')
      else
        call info0(10,1,0,' ... Make w(q) along lines in syml file'
     .    //' (mRy)')
      endif
      if (nomg < 1) then
        call info0(20,0,0,
     .    '     Positive SW freq not found at every kp ... '//
     .    'skipping SW spectra')
        goto 94
      endif

C     call prm3('J(q)-J(0)',0,Omb,nb1,nb2,nb3)
C     call prm3('w(q)',0,Om,nk1,nk2,nk3)

      call dpzero(h,2*nkfbz)
      call dcopy(nkfbz,Om,1,h,2)
C     call zprm3('omega(q)',0,h,nk1,nk2,nk3)
      call fftz3(h,nk1,nk2,nk3,nk1,nk2,nk3,1,0,-1)
C      call zprm3('omega(T)',0,h,nk1,nk2,nk3)
C     call dcopy(nkfbz,h,2,Om,1)

      nlb = 1
      lrdqp = cmdopt('-qp=',4,0,outs)
      if (lrdqp .or. cmdopt('--qp=',5,0,outs)) then
        call word(outs,1,j,k)
        if (lrdqp) then
          ifi = fopn(outs(5:k))
        else
          ifi = fopn(outs(6:k))
        endif
        lrdqp = T
      else
        lrdqp = F
        if (fxst('syml') /= 1) then
          call info0(20,0,0,
     .    ' ... no syml file ... skipping SW spectra')
          goto 94
        endif
        ifi =  fopno('SYML')

      endif
      ifib = fopnn('BNDS')
C      if (lrdqp .or. cmdopt('-con',4,0,outs)) then
C        rewind ifib
C      else
      write(ifib,336) nlb,0d0,0
C      endif
  336 format(i5,f10.5,i5)
  355 format(3f10.5/(10f8.4))

C ... Make qbb
      call pshpr(1)
      call bzmsh0(plat,llshft,0,nb1,nb2,nb3,is,ifac,rbb,qbb)
      call poppr

C ... Loop through all lines of k-points or read qp from file
      rewind ifi
      wint = 0
      wq(2) = 0
   91 continue
      if (lrdqp) then
        j = 3
        nq = 0
        if (rdm(ifi,10000,0,' ',xx,j,nq) /= 1) call
     .    rx('exasa:  bad qp file ' // outs(5:k))
        call ptr_bz(s_bz,1,'qp',3,nq,xx)
        rewind ifi
        j = rdm(ifi,10000,3*nq,' ',s_bz%qp,j,nq)
        call awrit1('%N exasa: read qp from file '//outs(5:k)//
     .    ', nq=%i',' ',80,stdo,nq)
        if (j /= 1) call rx('exasa: failed to read qp')
        call awrit2('%% rows %i cols %i',' ',80,ifib,nq,nlb+3)
      elseif (cmdopt('-con',4,0,outs)) then
        call rx('exasa not ready for contour plots')
      else
        read(ifi,*,end=92,err=92) nq,q1,q2
        write(*,785) nq,q1,q2
  785   format(' nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
      endif
      if (lrdqp .or. cmdopt('-con',4,0,outs)) then
      else
        write(ifib,337) nq
  337   format(2i5)
      endif
      call pshpr(iprint())
      do  iq = 1, nq
        call setpr(iprt(1))
        xx = dble(iq-1)/dble(nq-1)
        q(1) = xx*q2(1) + (1-xx)*q1(1)
        q(2) = xx*q2(2) + (1-xx)*q1(2)
        q(3) = xx*q2(3) + (1-xx)*q1(3)
        if (lrdqp) then
          call dpscop(s_bz%qp,q,3,iq*3-2,1,1d0)
        endif

        if (cmdopt('--swftitrp',10,0,outs)) then
C         This one can make wiggly spectra if not enough qp
          call pvgfeh(plat,h,q,nk1,nk2,nk3,wq)
        else
C         This one linearly interpolates spectra
C         call qpitrp(nk1,nk2,nk3,ifac,qb,q,Om,wq)
          call qpitrp(nb1,nb2,nb3,ifac,qbb,q,Omb,wq)
        endif

        wint = max(wint,dabs(wq(2)))
        k = nlst
        if (k == 0) k = nbas
        k = 1
        call dscal(k,1000d0,wq,1)

        if (iq >= 2) call setpr(iprt(1)-6)
        if (lrdqp) then
          if (cmdopt('-long',5,0,outs) .or.
     .      cmdopt('--long',6,0,outs)) then
            write(ifib,356) q, wq(1)
  356       format(3f10.6/(5f15.10))
          else
            write(ifib,357) q, wq(1)
  357       format(3f10.6/(8f10.6))
          endif
        else
          write(ifib,355) q, wq(1)
        endif
      enddo
      call poppr

      if (lrdqp) then
        return
      endif
      goto 91
   92 continue
      write(ifib,337) 0
C      if (wint > 1d-6) then
C        call info2(20,0,0,' ... finished bands:  max(Im omega) = %,3;3g'
C     .    ,wint,0)
C      endif
   94 continue

      endif
      endif

C -------------- Make G(R,R') for clusters ----------------
      if (lcgf == 12) then
      nblk = min(nbas,10)

C ... Potential functions P^alpha
      allocate(wk(nl*nsp*nbas))
      call pptrns(0,nl,ipdc,nclass,nsp,s_str%alph,
     .  size(ipdc),s_pot%pp,wk)
      deallocate(wk)

      i1 = 1
C#ifdefC DEBUG
C      j = 4
C      if (cmdopt('-i1=',j,0,outs)) then
C        if (.not. a2bin(outs,i1,2,0,' ',j,-1)) call
C     .    rxs2('EXASA: failed to parse "',outs(1:30),'%a"')
C        print *, 'use iz',i1
C      endif
C#endif
      if (nlibu > 0) then
        call vorbydl(-2,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,
     .    s_pot%pp,vorb)
      endif

C     Make G(E)
      if (lcgf == 12) then
        zpi(1) = efermi(1)
        zpi(2) = 1d-6
        call dpscop(zp,zpi,2,1,1,1d0)
        call dpscop(zp,zpi,2,2*nzp-1,1,1d0)
        call info2(0,1,-1,
     .    ' Energy at which to make r.s. GRR (file grr)? (def=%g,%g)? ',
     .    zpi(1),zpi(2))
        read(*,*) zpi

        call suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,dvldu,
     .    s_pot%vshft,zpi)

        allocate(dos(lidim,nsp))
        allocate(gRR(lidim,nspc,lidim,nspc,2))
        call dpzero(dos,lidim*nsp)
        do  isp = 1, nsp
          mode = 132
          if (isp == 1) mode = 152

          call dpzero(gRR,2*(lidim*nspc)**2)
          call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,mode,s_lat%pos,
     .      nlibu,lmaxu,idu,ludiag,vorb,nkp,s_bz%qp,s_bz%wtkp,
     .      nk1,nk2,nk3,s_bz%ipq,plat,s_str%s,s_str%iax,nptab,
     .      s_lat%istab,s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,qb,
     .      gRR,xx)

          i = lidim*nspc
          ifi = fopna('grr',-1,0)
          if (isp == 1) rewind ifi
          call ywrm(0,' ',2,ifi,'(9f12.6)',gRR,i**2,i,i,i)
C         call ztoyy(gRR,i,i,i,i,0,1)
          call daxpy(i,-1/pi,gRR(1,1,1,1,2),i+1,dos(1,isp),1)
        enddo

        call info0(0,1,0,' -1/pi Im G_ii resolved by channel:')
        do  isp = 1, nsp
          call info2(0,1,0,' spin %i',isp,0)
          print 543, dos(:,isp)
  543     format(9f9.4)
          print 544, sum(dos(1:i,isp))
  544     format(' sum: ', f12.4)
        enddo

        deallocate(dos,gRR)
        call fclose(ifi)
        call rx0('done')
        return

      endif

C ... Create cluster table for this lattice
      allocate(sid(nbas))
      call rx('exasa: update sid,ftiaxg.  ngdim not defined')
C     call spackv(10,'site sid',ssite,1,nbas,sid)
C      call ftiaxg(0,nk1,nk2,nk3,plat,s_lat%pos,nbas,s_ham%offH,sid,
C     .  ncl,oclp,oiaxb,oclssl,ngdim)
      call rx('ngdim is not defined')
      allocate(gii(ngdim*nsp))

      do  iz = i1, nzp

C   ... Hamiltonian setup for this energy
        call dpscop(zp,zpi,2,2*iz-1,1,1d0)
        call dpscop(wz,wzi,2,2*iz-1,1,1d0)
        call suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,dvldu,
     .    s_pot%vshft,zpi)

C   ... Make and save on disk G^ij in the irr. BZ for this energy
        isp = 1
        mode = 142
C#ifdefC CHKGRR
C        mode = 042
C#endif
        call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,mode,s_lat%pos,
     .    nlibu,lmaxu,idu,ludiag,vorb,nkp,s_bz%qp,s_bz%wtkp,
     .    nk1,nk2,nk3,s_bz%ipq,plat,s_str%s,s_str%iax,nptab,
     .    s_lat%istab,s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,qb,
     .    xx,xx)

C#ifdefC CHKGFJ
CC     Potential functions to scale to g^bare to g^alpha
C      call rx('chkgfj:use new mkpotf')
CC      call mkptfp(nl,nbas,nsp,lrel,ipdc,lihdim,s_ham%iprmb,130+order3,
CC     .  s_pot%pp,xx,xx,xx,dvldu,vshft,zpi,w(odpfb))
CC      call mkptfp(nl,nbas,nsp,lrel,ipdc,lihdim,s_ham%iprmb,120+order3,
CC     .  s_pot%pp,xx,xx,xx,dvldu,vshft,zpi,w(oddpfb))
CC     Potential functions to scale g^bare to G for free electrons
CC      call mkpotj(nl,nbas,nsp,ipdc,lmx,s_ctrl%rmax,avw,s_ham%iprmb,
CC     .  130,s_pot%pp,zpi,(0d0,0d0),w(odpfb))
CC      call mkpotj(nl,nbas,nsp,ipdc,lmx,s_ctrl%rmax,avw,s_ham%iprmb,
CC     .  120,s_pot%pp,zpi,(0d0,0d0),w(oddpfb))
CC     testing: scale g^bare to g^alpha; turn off g->G in gffbz
C      print *, 'testing: make g^alpha; turn off g->G in gffbz'
C      call mkptfp(s_site,s_pot,nl,nbas,nsp,lrel,.false.,ipc,ipdc,lihdim,pfdim,
C     .    s_ham%iprmb,40+order3,s_pot%pp,s_pot%pprel,xx,xx,xx,dvldu,vshft,zpi,w(odpfb))
C      call mkptfp(s_site,s_pot,nl,nbas,nsp,lrel,.false.,ipc,ipdc,lihdim,pfdim,
C     .  s_ham%iprmb,50+order3,s_pot%pp,s_pot%pprel,xx,xx,xx,dvldu,vshft,zpi,w(oddpfb))
CC     testing: scale g^bare to g^alpha; turn off g->G in gffbz
C      call gfzkin(1,zpi,efree,zkin)
C      call mkpotj(nl,nbas,nsp,ipdc,lmx,s_ctrl%rmax,avw,s_ham%iprmb
C     .  ,40,s_pot%pp,zkin,(0d0,0d0),w(odpfb))
C      call mkpotj(nl,nbas,nsp,ipdc,lmx,s_ctrl%rmax,avw,s_ham%iprmb
C     .  ,50,s_pot%pp,zkin,(0d0,0d0),w(oddpfb))
C
C#endif

C#ifdefC CHKGRR
C      ifi = fopnx('gfqp',100,16+8+4+0,-1)
C      rewind ifi
C      allocate(gij(k1*k2*k3*lidim*lidim))
C      call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,ib2,1441,
C     .  s_lat%pos,nkp,s_bz%qp,zpi,nk1,nk2,nk3,k1,k2,k3,s_bz%ipq,plat,
C     .  s_lat%istab,s_lat%symgr,s_lat%ag,s_bz%star,ifac,qb,offr,offc,
C     .  0,gij,xx,xx)
C
CC ... Potential functions P^alpha & GF for reference potential
C      zpi(1) = zpi(1) + .1d0
C      zpi(2) = zpi(2) + .1d0
C       call suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,s_pot%vshft,zpi)
C      mode = 042
C        call defcc(ogref,k1*k2*k3*lidim*lidim)
C        call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,mode,s_lat%pos,nkp,
C     .    s_bz%qp,w(owtkp),nk1,nk2,nk3,s_bz%ipq,plat,s_str%s,s_str%iax,nptab,
C     .    s_lat%istab,s_lat%symgr,s_lat%ag,s_bz%star,ifac,qb,xx,xx,0,F,xx))
C        call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,ib2,1441,
C     .    s_lat%pos,nkp,s_bz%qp,zpi,nk1,nk2,nk3,k1,k2,k3,s_bz%ipq,plat,
C     .    s_lat%istab,s_lat%symgr,s_lat%ag,s_bz%star,ifac,qb,offr,offc,
C     .    0,w(ogref),xx,xx)
C      opfb = s_pot%opf
C      odpfb = s_pot%odpf
C      oddpfb = s_pot%oddpf
C      zpi(1) = zpi(1) - .1d0
C      zpi(2) = zpi(2) - .1d0
C
CC ... Compare to dyson equation
C      call defcc(oh,k1*k2*k3*lidim*lidim)
C      call gfqdys(nk1,nk2,nk3,lidim,nl**2*nbas,s_pot%pf,w(opfb),
C     .  gij,w(ogref),h,w(oh2))
C
C      goto 40
C#endif

      ifi = fopnx('gfqp',100,16+8+4+0,-1)
      rewind ifi
      call pshpr(1)
      lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MPL*0+MCD*1+MZP*0)
      if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zpi,q,plat,xx,xx,xx,xx,0,0,
     .  0,xx) /= 0) call rx('exasa: failed to read file header')
      call poppr

C     for now, instead of clusters.  Also, should loop over spin
      jb = 1
      jb2 = nbas
C      offr = (isp-1)*nlma*ldimi
C      offc = (isp-1)*ldimi*nlma
      call offsHp(s_ham%iprmb,jb,jb2,1,0,lidim,offj,nj)
      do  ib = 1, nbas, 1

        ib2 = ib
        call offsHp(s_ham%iprmb,ib,ib2,1,0,lidim,offi,ni)

        allocate(gij(k1*k2*k3*ni*lidim*nsp))
        allocate(gji(k1*k2*k3*lidim*ni*nsp))
C       ib2 = min(ib-1+nblk,nbas)
C       jb2 = min(jb-1+nblk,nbas)
C       offr = w(offH+nkap0*n0H*(ib-1))
C       ni   = w(offH+nkap0*n0H*ib2) - offr
C       offc = w(offH+nkap0*n0H*(jb-1))
C       nj   = w(offH+nkap0*n0H*jb2) - offc

        offr = 0
        offc = 0

C   ... Expand G to full BZ for subblock connected to site ib
        call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,ib2,440,
     .    s_lat%pos,nkp,s_bz%qp,zpi,nk1,nk2,nk3,k1,k2,k3,s_bz%ipq,plat,
     .    s_lat%istab,s_lat%symgr,s_lat%ag,s_bz%star,ifac,qb,offr,offc,
     .    0,gij,gji,xx)

C#ifdefC CHKGFJ
CC   ... Shift tau = pos(jb)-pos(ib)
C        call dpscop(s_lat%pos,tau,3,3*jb-2,1,1d0)
C        call dpsadd(tau,s_lat%pos,3,1,3*ib-2,-1d0)
C        call gvctof(0,1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx)
C        allocate(gv(ngmx*3))
C        allocate(kv(ngmx*3))
CC   ... Make list for a dimensionless lattice (alat=1)
C        call ogvlst(1d0,qlat,tau,nk1,nk2,nk3,gmax,ngmx,ng,gv,kv)
C
C          call rxx(nblk /= 1,'not ready for nblk>1')
C          rmaxi = dval(s_ctrl%rmax,ipdc(ib))
C          rmaxj = dval(s_ctrl%rmax,ipdc(jb))
C          call chkgfj(nk1,nk2,nk3,nsp,ib,jb,ni,nj,offi,offj,gij,
C     .    alat,avw,rmaxi,rmaxj,zkin,gv,kv,ng,tau,nl**2*nbas,
C     .      w(odpfb),w(oddpfb),s_lat%indxcg,s_lat%jcg,w(ocg),s_lat%cy)
C          goto 42
C#endif

        call rx('update call to gfq2rs')
C        call gfq2rs(nk1,nk2,nk3,nsp/nsp,ib,ib2,jb,jb2,s_ham%offH,ni,nj,
C     .    gij,plat,s_lat%pos,ncl,w(oclp),w(oiaxb),w(oclssl),gii)
      enddo

      ifi = fopna('pemb',-1,4)
      i = nl**2*nbas*nsp
      call ywrm(1,' ',3,ifi,' ',s_pot%pf,i,i,i,1)
      call fclose(ifi)

      i = (32+16+8+4+1)*10000 + 1
      if (iz /= 1) i = 5
      print *, 'doctor i since iz ne 1'
      i = (32+16+8+4+1)*10000 + 1
      call dpzero(q,3)
      call rx('update call to iogfrs')
C      j = iogfrs(i,1,0,'gemb',ifi,1,ncl,nbas,0,zpi,q,plat,s_lat%pos,
C     .  w(oclp),w(oclssl),w(oiaxb),0,0,3,gii)
      call fclose(ifi)
      stop
      enddo
      if (nlibu > 0) then
        call vorbydl(-1,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,
     .    s_pot%pp,vorb)
      endif

      endif

      return

C --- Error exit ---
C  99 call rx('EXASA: failed to read evecs')

      end

      subroutine getqsw(strn,nclasp,nbas,ipc,amom,aamom)
C- Read sphere charges from command-line switch
      implicit none
C ... Passed parameters
      integer nclasp,nbas,ipc(*)
      double precision amom(nclasp),aamom(nbas)
      character strn*(*)
      integer h(max(nclasp,nbas))
C ... Local parameters
      integer i,j,a2vec
      character dc*1
      logical ltmp

      i = 7
      dc = strn(i:i)
      ltmp = dc == 's'        ! moments by site, not spec
      if (ltmp) then
        i = i+1
        j = a2vec(strn,len(strn),i,4,dc//', ',3,3,nbas,h,aamom)
      else
        j = a2vec(strn,len(strn),i,4,dc//', ',3,3,nclasp,h,amom)
      endif
      if (j >= 0 .and. ltmp) then
        call info5(10,0,0,'%8fRead %i moment(s) by site from switch :'//
     .    ' %n:1d',j,j,aamom,0,0)
        call cpvprm(1,1,j,ipc,aamom,amom)
      elseif (j >= 0) then
        call info5(10,0,0,'%8fRead %i moment(s) by class from switch :'
     .    //'%n:1d',j,j,amom,0,0)
        call cpvprm(0,1,nbas,ipc,amom,aamom)
      else
        call rxs('exasa: failed to parse switch',strn)
      endif


      end
