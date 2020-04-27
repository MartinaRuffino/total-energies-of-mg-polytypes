C#define SX
C#define v2
      subroutine gfasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,
     .  s_str,s_strn,lldau,vorb,dmatu,vconst,vshft,zp,wz,pgplp,efermi,
     .  moddos,sumev,ehterm,amag,aamom,sxad,rms)
C- Crystal Green's function, ASA
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspin lcgf lrel lncol lsx lordn lham
Ci                 lasa ldlm nccomp nbasp nspec ldomg lgen3 nclasp
Co     Stored:     lncol ldomg
Co     Allocated:  *
Cio    Elts passed: lncol lscr ldos ics ncomp ipc lrel nrc rmax lbxc
Cio                lham lasa dclabl initc
Cio    Passed to:  suhamz mkpotf makpfz mkcpa gfibz gfomg gfp0ft
Cio                hamfbz gfenint mixomg mkgint shoctl gfgw gfsigma
Cio                magcpa
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nsgrp awald nkd nkq alat nkdmx nkqmx as tol vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos istab symgr ag dlv qlv
Cio    Passed to:  mkcpa gfp0ft gfgw gfsigma
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z idu a nr rmt lmxa hcr name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  getidu suhamz mkpotf makpfz vorbydl shoctl iorbtm
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp dlmcl spec class pnu clabel v0 domg omg
Ci                 bxc omgn alr thet cpawt
Co     Stored:     *
Co     Allocated:  pdos dmat dpfr ddpfr pfr sfvrtx
Cio    Elts passed: j0 alr amr amrt pdos dmat thet dlmcl cpawt gc gcorr
Cio                norb ddpfr dpfr pfr omgn domg gcu sfvrtx ncomp omg
Cio    Passed to:  getidu suhamz mkpotf makpfz mkptfpd mkfrpf mksopf
Cio                mkcpa vorbydl gfomg mixomg womg mkgint gfg2gr cpasf
Cio                packdmat magcpa cpomg
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc lshft semsh dosw ndos lmet qp zval
Co     Stored:     *
Co     Allocated:  qp
Cio    Elts passed:qp wtkp ipq star sopertp
Cio    Passed to:  gfgw gfsigma
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  nlibu lmaxu udiag lsig eterms ldham neula lncol hord
Ci                 lgen3 bandw ndhrs qss offH
Co     Stored:     lncol eterms
Co     Allocated:  *
Cio    Elts passed: iprmb eula bdots offH etrms lncol neula qss nprs
Cio                iaxs hrs
Cio    Passed to:  suhamz mkpotf makpfz mkcpa gfibz gf1kp gfg2g gfomg
Cio                gfp0ft gfp0f1 pdglrs gfdpp gfg0g gfenint mkgint
Cio                gfg2gr gfgw gfsigma
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves pfr pf palp dpfr ddpfr dpf ddpf
Co     Stored:     socscl
Co     Allocated:  socscl qcorr cp pf dpf ddpf dddpf pfr dpfr ddpfr
Co                 gmar papg palp gma
Cio    Elts passed: gibbs qc socscl mxy pp qcorr sop pf vintr pprel dpf
Cio                ddpf dddpf palp gma pfr ddpfr dpfr papg gmar cp pnu
Cio                qnu ves qnur gfr
Cio    Passed to:  suhamz mkpotf makpfz mkcpa gfibz gf1kp gfg2g gfomg
Cio                gfp0ft gfp0f1 pdglrs gfdpp gfg0g gfenint mkgint
Cio                shoctl gfgw gfsigma
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:s iax
Cio    Passed to:  *
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  str_pack
Ci Inputs
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci   vorb  :orbital dependent potential matrices
Ci   dmatu :density matrix for LDA+U
Ci   zp    :points for complex energy contour
Ci   wz    :weights for complex energy integration
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          This is passed for consistency with pgfasa.
Ci   efermi:Fermi energy
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   qnur  :relativistic energy-weighted moments of the sphere charges
Ci   sxad
Cio Inputs/Outputs
Cio  vconst:global constant potential shift added to estat potential; it
Cio        :should be chosen so that the system is neutral for the given
Cio        :Fermi level (input, and fixed quantity).  vconst is not
Cio        :known in advance, but must be determined by gfasa.  Local
Cio        :variable vbar accumulates the correction to vconst until
Cio        :gfasa is finished.
Cio        :Input vconst is estimate for this shift;
Cio        :Output vconst is value that makes the system neutral.
Cio  vshft :array of site potential shifts that may be added.
Cio         See iovshf.f for array structure.
Co   vshft :(if output moments generated) adjusted for charge neutrality
Co Outputs:
Co   moddos: flags marking what quantities were made
Co   sumev :sum of single-particle energies
Co   ehterm:double-counting correction added to ehterm
Co   rhos  :spin density-matrix
Co   amag  :net magnetization
Co   aamom :local magnetic moments
Co   rms   :(CPA only) largest RMS mismatch for Omega among all CPA sites
Cs Command-line switches
Cs   --mem  : Use disk to conserve memory
Cs   --onesp: Generate bands for one spin (not implemented)
Cl Local variables
Cl   delef :  cumulative energy shift for lpade=2, in search for Fermi level
Cl  lagain :0 => quantities from diagaonal G calculated from true GF
Cl         :1 => quantities from diagaonal G estimated by Pade
Cl         :<0 => same as 0, but internal iteration
Cl   lpade :0 Diagonal part of GF is not kept
Cl         :1 Used to interpolate gf to a elliptical contour shifted
Cl         :  by a constant energy from the one actually used for
Cl         :  computation of the GF.  The interpolation is computed
Cl         :  using a Pade approximation.  In this mode, the diagonal
Cl         :  part of the gf are retained for later interpolation to a
Cl         :  shifted contour to the estimated Fermi energy, once the
Cl         :  deviation from charge neutrality is known.
Cl         :2 Used in conjunction with Simpson integration along the
Cl         :  real axis when searching for the Fermi level.  In this
Cl         :  mode, the diagonal part of the energy-dependent gf is
Cl         :  retained as in lpade=1; here, however, the points
Cl         :  themselves are updated as new points are generated by
Cl         :  explicit computation of the GF.  The Pade interpolation is
Cl         :  used to estimate the midpoint for the Simpson rule.
Cl         :  See Remarks above line ... elseif (lpade == 2) then
Cl   leres :  set to 1 if "energy resolution mode", intended for
Cl         :  generation of energy-resolved quantities.  Energy
Cl         :  points are on a uniform mesh close to real axis
Cl         :  In this mode, gfasa stops after quantities are generated
Cl   lrect :a parameter with the trapezoidal integration along the real
Cl         :  axis (see lpade=2 above):
Cl         :0 gf has only been created for one point (diagonal gf
Cl         :  only known for 1 point).
Cl         :2 diagonal gf has been saved for a prior point
Cl   gfopts:string containing switches for GF program; see Remarks
Cl   fgemb :string contain file info for embedded GF
Cl   vshfn :estimate of potential shift relative to prior iteration
Cl         :that makes the system charge neutral for the supplied Fermi level.
Cl   vhldp :a work array used by gfzerq to bracket potential shifts
Cl         :for charge-neutrality point, in a Pade approximation to G
Cl   vsint :a work array to bracket potential shifts over internal loops
Cl         :where g is explicitly calculated.  It plays a role similar
Cl         :to the Pade, except that points kept are real, not estimates.
Cl   vbar  :cumulative change potential shift for charge neutrality.
Cl         :At end of gfasa cycle, vconst is incremented by vbar
Cl   vbardc:same as vbar, but only shift that contributes to a shift in
Cl         :double-counting.  Some Fermi-level finders work by shifting
Cl         :the potential in the search for c.n, and require an extra
Cl         :double-counting term.  Fermi-level finders that vary the fermi
Cl         :level and add a potential shift later do not.
Cl   dosn  :keeps track of dos,idos for on-the-fly trapezoidal
Cl         :integration when integrating along the real axis.
Cl   dosi  :integrated dos and moments of dos
Cl          1: dos at this e : -wtkp/pi Im G
Cl          2: nos at this e : -wtkp/pi Im (wz*G)
Cl          3: 1st energy mom : -wtkp/pi Im (wz*G z)
Cl          4: 2nd energy mom : -wtkp/pi Im (wz*G z**2)
Cl          5: projection of dos to fermi level (not calculated here)
Cl   isx    1: generate screened exchange
Cl ninternal: number of iterations where the true G has been integrated over frequency
Cl          :Used in to estimate the true Fermi level.
Cr
Cr Remarks
Cr Branches of execution.  Branches are mutually exclusive.
Cr   'ctrl lcgf' governs which branch is taken.
Cr     1: Generate diagonal GF and output density
Cr Options, which may apply to some or all of the major branches.
Cr Those options inapplicable to a particular branch are ignored
Cr Options are passed in a string form, delimited by ';', in gfopts.  They are:
Cr   p1,p3, or pz  order of potential function
Cr   idos        integrated DOS
Cr   noidos      suppress integrated DOS
Cr   dmat        make the density-matrix
Cr   sdmat       make the site-diagonal density-matrix, write to file dmat.ext
Cr   pdos        accumulate partial DOS
Cr   emom,noemom accumulate output moments, or suppress accumulation
Cr   ibloch      inverse Bloch transform
Cr   frzvc       Suppress saving potential shift used to determine charge neutrality
Cr   shftef      Find charge neutrality point by shifting the Fermi level,
Cr               rather than adding a constant potential shift
Cr   padtol=#    tolerance for maximum potential shift permissible by Pade interpolation
Cr   qtolp=#     A second tolerance for deviation from charge neutrality.  
Cr               If tag is missing or #=0, this switch has no effect.
Cr               Otherwise, when Deltaq < #, code exits search for neutral point.
Cr               See Questaal web page for more details.
Cr  The following are specific to the CPA:
Cr   omgtol=#    Tolerance in the Omega potential, CPA self-consistency    
Cr   omgmix=#    How much of prior iterations to mix, CPA self-consistency 
Cr   nitmax=#    Maximum number of iterations for CPA self-consistency     
Cr   lotf        Learn on-the-fly                                           
Cr   specfun     Make spectral function                                     
Cr   specfr      Make spectral function resolved by site                    
Cr   specfrl     Make spectral function resolved by site and l              
Cr   dmsv        Record density matrix to file                              
Cr   dz=#        Shift Omega potential by dz                                
Cr   sfrot=#     Rotate the GF in spin space by this angle around the y axis
Cb Bugs
Cb   Order-N not working
Cu Updates
Cu   17 Jul 18 Implemented qtolp
Cu   18 Jun 18 Synchronize with updated spherical harmonics
Cu   25 Mar 16 (MvS) Restructured FR construction of ppars and scaling
Cu   25 Mar 16 (MvS) For scaling g->G use gfg2gnc
Cu   09 Apr 15 (Vishina, MvS) more changes preparing for Dirac case
Cu   08 Jun 14 (Belashchenko) modifications for the relativistic case
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   12 Jan 12 (Belashchenko) rewrote DLM case
Cu   23 Nov 10 (Larson) added DLM case
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Oct 10 Read strux through rdstrx
Cu   08 Nov 07 (J. Xu) LDA+U implementation, first cut
Cu   17 Mar 05 selfconsistent sigma (T.Sandu)
Cu   05 Mar 05 New frzvc and updated lpade=2 mode
Cu   15 Dec 04 first draft selconsistent sigma (T. Sandu)
Cu   20 Oct 04 Adapted to spin-polarized screened exchange(T. Sandu)
Cu   29 Sep 04 Fix some problems with c.n. pot shifts
Cu   27 Sep 04 implemented screened exchange(T. Sandu)
Cu   17 Sep 04 Adaptations to non-equilibrium mode for pgf
Cu   18 Jun 04 (A Chantis) orbital moments for fully relativistic case
Cu   02 Mar 04 first attempt to implement screened exchange(T. Sandu)
Cu   06 Jun 03 introduced semsh(2)=20: integration proceeds as in
Cu             semsh(2)=10 but no energy-integrated quantities are
Cu             calulated. Does no charge neutrality
Cu   04 Jun 03 Allow energy mesh mode 2 to write partial dos
Cu   21 Mar 03 redesigned potential shifts; new argument list.
Cu   10 Mar 03 Poke double-counting from appl. field into eterms
Cu   24 May 02 First cut at making q-dependent response function
Cu   21 May 02 deletes file gfqp on exit, if it exists.
Cu   26 Feb 02 small changes to adapt to changes in layer code
Cu             vshft is now distributed over sites internally.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer pgplp(6,0:1)
      double precision amag(3),aamom(*),vshft(*),
     .  efermi,zp(2,1),wz(2,1),sumev,ehterm(4),vconst,rms
C ... For LDA+U
      integer lldau(*)
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
C     type(str_ordn)::  s_ordn
      type(str_strn) :: s_strn(*)
C ... Dynamically allocated local arrays
      integer, pointer :: gflst(:),glstp(:)
      integer, allocatable :: idu(:),ncomp(:),dlmcl(:)
C     integer, allocatable :: lmx(:)

      real(8), allocatable :: z(:),pdos(:),vshfn(:),vhldp(:),vsint(:),dvldu(:),qnusav(:)
      real(8), allocatable :: potp(:),pph1(:),wk(:),oo(:),orbtm(:)  !,qlm(:,:,:)
C     real(8), allocatable :: sali(:)
      real(8), pointer :: rhos(:)

      complex(8), pointer :: gd(:,:,:,:),gds(:,:,:,:)
      complex(8), allocatable :: dgdk(:),dgd(:),gii(:,:),ghh(:,:),znp(:)
      complex(8), allocatable :: wscrl(:),wscr(:),wscrk(:)
      complex(8), allocatable :: v0(:),v01(:),v02(:),p0(:),madm(:)
      complex(8), allocatable :: sqrtdpp(:),gfdmtx(:)
      complex(8), allocatable :: dmat(:),gwk(:,:)
C ... Local parameters
      logical F,T,lccpa,lemom,lnos,lcon,
     .  lrhos,ltmp,lpdos,lfrzvc,lso,ldmsit,ldmsv,efshft
      integer hdim,i,ibzsw,ierr,ipr,iprt,iscr,iscr1,isp,isx,isx1,iz1,
     .  izp,lasa,lcgf,ldim,ldmat,ldos,leres,lham,lhdim,lhdimx,lidim,lidimx,lagain,
     .  lmet,lncol,lordn,lpade,lrel,lrsig,lx2,lswd,moddos,mxorb,nRlc,nbas,nbl,
     .  nclasp,nclass,ndos,neul,ninternal,nk1,nk2,nk3,nkp,nl,nl2,nlo,nlspc,
     .  nlspcr,npados,npdos,nptab,nrhos,nsp,nspc,nzp,offd,offpd,scrwid,stdl,stdo
      integer is(3),ifac(3),nkabc(3),ldham(16),lshft(3)
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lhdim,ldham(3))
      equivalence (nspc,ldham(4)),(lidimx,ldham(6)),(lhdimx,ldham(7))
      double precision xx,rmsdv,rmsdq,semsh(10),vbar,growdv,
     .  vbardc,zpi(2,3),wzi(2,3),dosi(5,2),dosn(3,2),eterms(22),delef
      double precision plat(3,3),padtol,dosw(2),zval,qtol,qtolp,qtot,qp2(3),dum(20),dosir(5,2)
      real(8), target :: xxp(1,1)
      character*120 gfopts,fgemb,strn,plbopt
      procedure(integer) :: bitand,fhndl,fopn,fopna,str_pack,iprint,nglob,nnrl,partok,isw
      procedure(logical) cmdopt,bittst

C ... For finding charge neutrality point
      integer nznp,nzi,lrect
      double precision de,qtotl,dqsimp,defzpn
C ... For symmetry operations and rotations of G
      logical llshft(3)
      integer nsgrp
      double precision rb(3,3),qb(3,3)
      parameter(T=.true.,F=.false.,scrwid=80)
      character*120 outs
      equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
C ... For Order-N
C      integer lemb,ngdim
C      double precision efree
C      debugging for relativistic case
C      complex(8),allocatable:: pfall(:,:,:,:,:)
C      double complex zxx
C      integer n1,l,m,n4,n5,n6,k,lmr,opfall,izz,nzz
C      common /snott/ izz,nzz,opfall
C ....new parameters for path integration of P0
C     double precision eminim,emaxim,zp1(2,1),wz1(2,1)
C     integer p,p1,nzpr1,nzpr2
C....starting definitions for screened exchange sigma
      logical aintra
      double precision qp(3),awald,alat,vol
      integer iq,nkd,nkq,lpdim,lp2
C     integer ov0,ov01,ov02,op0,owscr,owscrl,owscrk
      double precision awald0,tol
      integer nkdmx,nkqmx
      integer pas
c      parameter (T=.true., F=.false.)
C ... For CPA
      integer ldlm,ldlm1,nangl,nitmax,moddos2,moddos3,icomp
      integer nclassd,norb,nitomg,ib
      real(8) omgtol,omgmix,dz
      logical lmixom,lotf,lweiss,lgcorr,omonly,ldomg,lsev
      character *12 dosfil
      character *4 aib
C ... For spectral functions
      real(8), allocatable :: spf(:,:,:),spfr(:,:,:,:),spfrl(:,:,:,:,:)
      real(8), allocatable :: sfint(:,:)
      real(8) :: sfth
      real(8) :: efermsf = 0d0
      integer onesp,iopq,i1,i2,nfbn(4),modsf
      logical lspec,lspecr,lspecrl

C ....new parms selconsistency sigma
      logical sxad
      double precision tolsig
C ... For LDA+U
      integer nlibu,lmaxu,ludiag
c     logical mlog
      integer mpipid,procid,master,numprocs
C ... For full rel
      integer nclspd, ic

      call tcn('gfasa')
      procid = mpipid(1)
      master = 0
      numprocs = mpipid(0)
C     call wkprnt(1)
      nullify(gd,gds)
C     The DEC alpha compiler sometimes require passed arrays to be allocated
      if (.not. allocated(pdos )) allocate(pdos (1))
      if (.not. allocated(orbtm)) allocate(orbtm(1))
      if (.not. allocated(dmat )) allocate(dmat (1))
      if (.not. associated(gd  )) allocate(gd(1,1,1,1)) ! hopefully not a leak..

C --- Unpack parameters, set defaults ---
      nkp = s_bz%nkp
      nkabc = s_bz%nkabc
      lshft = s_bz%lshft
      semsh = s_bz%semsh
      dosw = s_bz%dosw
      ndos = s_bz%ndos
      lmet = s_bz%lmet
      nlibu = s_ham%nlibu
      lmaxu = s_ham%lmaxu
      ludiag = s_ham%udiag
      if (semsh(2) /= 2) then
        leres = 0
      else
        leres = 1
        ndos = semsh(1)
        dosw(1) = semsh(3)
        dosw(2) = semsh(4)
      endif
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      growdv = 1                ! May be set larger if vshft is systematically underestimated

C...for screened exchange sigma......................................
      awald = s_lat%awald
      nkd = s_lat%nkd
      nkq = s_lat%nkq
c      print *,'oqlv','awald',awald
c      call prmx('qlv',s_lat%qlv,3,3,nkd)
c
c      stop
      alat = s_lat%alat
C     oaamom = s_pot%oaamom
      nkdmx = s_lat%nkdmx
      nkqmx = s_lat%nkqmx
c      print *, 'nkdmx',nkdmx,'nkqmx',nkqmx,'nkd',nkd,'nkq',nkq
c      call prmx('nkdmx',nkdmx,1,1,1)
c
c      stop
      awald = s_lat%awald
      awald0 = s_lat%as
      tol = s_lat%tol
      vol = s_lat%vol
C.....end new unpacking for sigma

      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lcgf = s_ctrl%lcgf
      lrel = mod(s_ctrl%lrel,10)
C ... If Dirac, unset the spin-orbit flag
!      if (mod(s_ctrl%lrel,10) == 2 .and. bittst(s_ctrl%lncol,4)) then
!        s_ctrl%lncol = s_ctrl%lncol - 4
!        s_ham%lncol = s_ctrl%lncol
!      endif
      lncol = s_ctrl%lncol
      lso = bittst(lncol,4)
      if (lso .and. lrel == 2)
     .  call rx('lso and lrel=2 in gfasa: should never happen')

      isx = s_ctrl%lsx
      lordn = s_ctrl%lordn
      lrsig = s_ham%lsig
      sxad = lrsig /= 0
      iscr = mod(s_ctrl%lscr,10)
      isx1 = mod(isx,10)
      iscr1 = mod(mod(iscr,10),2)
      if (iscr1 == 4) iscr1 = 0
      aintra = F
      if (mod(isx,2) == 1 .or. mod(iscr1,2) == 1) aintra = T
      lham = s_ctrl%lham
      lasa = s_ctrl%lasa
      ldos   = IAND(s_ctrl%ldos,7)
      nlo = nl
      if (ldos >= 4) nlo = nl*nl
C      if (lordn == 1) then
C        osordn = s_ctrl%osordn
C      endif
C      gfopts = ' '
C      if (i2 >= i1) gfopts = sstrn(i1:i2)
      i = str_pack('gfopt',-2,s_strn,gfopts)
C      fgemb = ' '
C      if (i2 >= i1) fgemb = sstrn(i1:i2)
      i = str_pack('gemb',-2,s_strn,fgemb)

C     vbar = change in vconst to ensure charge neutrality
C     At the close of gfasa, vconst is incremented by vbar
      vbar = 0
      vbardc = 0

      padtol = .03d0; qtolp = 0; qtol = 0; lrect = 0
      omgtol = .5d0; omgmix = .5d0
      nitmax = 20
      s_ctrl%ldomg = 0; lotf = F; lspec = F; lspecr = F; lspecrl = F; ldmsit = F; ldomg = F
      dz = 1d-5
      sfth = 0
      efermi = semsh(4)
      stdo = nglob('stdo')
      stdl = nglob('stdl')

C     Debugging
C     call shoctl(s_ctrl,s_spec,s_pot,101,stdo); stop

C --- Global setup: array dimensions ---
      eterms = s_ham%eterms
      call dvset(eterms,17,18,0d0)

      ldham = s_ham%ldham
C     offH = s_ham%ooffH
      nzp = nint(semsh(1))
C     lrhos  = bittst(lncol,16)
      lrhos  = nspc == 2

      ldlm = s_ctrl%ldlm
      nangl = s_ctrl%nccomp
      nclasp = s_ctrl%nclasp
      if (ldlm /= 0 .and. nsp /= 2)
     .  call rx('Bug: CPA needs NSPIN=2')
      if (ldlm == 0) nangl = 0
      nclassd = nclass + nangl
      nclspd = nclasp + nangl

      if (ldlm /= 0) then
        if (isx /= 0)   call rx('GFASA: isx and CPA')
        if (iscr1 == 1)  call rx('GFASA: iscr and CPA')
        if (lordn /= 0) call rx('GFASA: lordn and CPA')
        if (sxad)         call rx('GFASA: sxad and CPA')
        if (procid == master .and. lrhos)
     .    print *, 'WARNING (gfasa): lrhos and CPA'
      endif
C      if (procid == master .and. bittst(s_ctrl%lbxc,4))
C     .  print *, 'NOTE: Constraining fields fixed as of 07/21/14'
      ldlm1  = mod(ldlm/10,10)
      lsev = mod(ldlm/100,10) == 1
      if (lsev .and. lrel == 2) call rx('lsev and rel=2')
      lweiss = mod(ldlm1,2) == 1
      omonly = ldlm1 >= 2

      nlspc = nl * nsp * nclassd
      nlspcr = 8*nl*nl*max(nclassd,s_ctrl%nspec)
      lx2 = lidimx**2
      nl2 = nl*nl
      mxorb = nglob('mxorb')
      hdim = lhdim-lidim
      rhos => s_pot%rhos; nrhos = s_pot%nrhos

      if (nlibu > 0 .and. ludiag==0 .and. bitand(lham,128) == 0)
     .  call rx('Off-diagonal LDA+U requires gamma repsn')

      if (nlibu > 0 .and. ludiag == 0 .and. ldlm > 0)
     .  call rx('Off-diagonal LDA+U does not allow CPA')

C ... Class-based arrays
      if (allocated(z)) deallocate(z); allocate(z(nclassd))
      call spec2class(s_spec,nclassd,s_ctrl%ics,'z',1,xx,z)
      if (ldlm /= 0) then
        call cpaz(nclass,z,s_ctrl%ncomp,s_pot%gibbs)
        call cpaz(nclass,s_pot%qc,s_ctrl%ncomp,s_pot%gibbs)
      endif

C --- Unpack GF options ---
C ... Default nos,emom ON for contour integral
      i = mod(nint(semsh(2)),100)
      lcon = i == 10 .or. i == 11 .or. i == 20
      lnos = lcon .or. i == 2
      lemom = lnos
      ldmat = 0; lpade = 0; lgcorr = F; lfrzvc = .false.; efshft = .false.
      if ((lnos .or. lcon) .and. leres == 0 .and. lmet /= 0 .and. isx1 == 0) then
        lpade = 1
        i = mod(nint(semsh(2))/100,10)
        if (i == 3) lpade = 2   ! vertical contour starting on the real axis
        if (i == 4) lpade = -2  ! ?? 10s digit mode = 4 is not documented
        semsh(8) = 6
      endif

      call partk0(0,len(gfopts),1,-1,0,len(gfopts),-1,31,F)
      i = partok(gfopts,'emom',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) lemom = T
      i = partok(gfopts,'noemom',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) lemom = F
      i = partok(gfopts,'idos',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) lnos = T
      i = partok(gfopts,'sdmat',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) ldmat = 1
      i = partok(gfopts,'dmat',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) ldmat = 2
      i = partok(gfopts,'noidos',' ;',ltmp,' ',0,0,0,0)
      if (ltmp) lnos = F
      if (ltmp) lemom = F
      i = partok(gfopts,'pdos',' ;',lpdos,' ',0,0,0,0)
      i = partok(gfopts,'dmsv',' ;',ldmsv,' ',0,0,0,0)
      i = partok(gfopts,'shftef',' ;',efshft,' ',0,0,0,0) ! Shift Fermi level for c.n., no vshft
      i = partok(gfopts,'frzvc',' ;',lfrzvc,' ',0,0,0,0)
      i = partok(gfopts,'ibloch',' ;',ltmp,' ',0,0,0,0)
      i = partok(gfopts,'padtol=','= ;',padtol,' ',-1,4,0,0)
      i = partok(gfopts,'qtolp=','= ;',qtolp,' ',-1,4,0,0)
!      i = partok(gfopts,'qtol=','= ;',qtol,' ',-1,4,0,0) Not ready yet
      i = partok(gfopts,'omgtol=','= ;',omgtol,' ',-1,4,0,0)
      i = partok(gfopts,'omgmix=','= ;',omgmix,' ',-1,4,0,0)
      i = partok(gfopts,'nitmax=','= ;',nitmax,' ',-1,2,0,0)
      i = partok(gfopts,'lotf',' ;',ltmp,' ',0,0,0,0)
      if (ltmp .and. ldlm /= 0) lotf = T
      i = partok(gfopts,'specfun',' ;',lspec,' ',0,0,0,0)
      i = partok(gfopts,'specfr',' ;',lspecr,' ',0,0,0,0)
      i = partok(gfopts,'specfrl',' ;',lspecrl,' ',0,0,0,0)
      i = partok(gfopts,'dz=','= ;',dz,' ',-1,4,0,0)
      i = partok(gfopts,'sfrot=','= ;',sfth,' ',-1,4,0,0)

c     if (dz /= 0d0 .and. .not. lsev) then
c       ldomg = T
c       omonly = T
c       s_ctrl%ldomg = 1
c     endif

      if (lspecrl) lspecr = T
      if (lspecr) lspec = T
      if (lspec) lotf = F  ! use Omega from disk
      if (lspec .and. lsev) call rx('Spectral function and DLM=100')

C ... Reentry point for the shifted energy loop
    6 continue
      if (omonly .or. lspec) then
        lpade = 0
        lemom = F
        lnos = F
        lpdos = F
      endif

      if (lspec) omonly = F
c     if (lsev) then
c       lpade = 0
c       lpdos = F
c       lotf = F
c     endif

C     Density matrix is required in LDA+U
      if (nlibu > 0) ldmat = 1
      if (ldlm /= 0 .and. nlibu > 0)
     .  print *,'WARNING: only IDU=4 or 5 are allowed with CPA'

C     Always make site density matrices in the spin-coupled case
      if (ldlm /= 0 .and. nspc == 2) ldmat = 1
c     if (ldlm == 0 .and. lspec) call rx('specfun needs CPA')

C ... If Pade, Omega(z), and Omega(z+dz) are already done
      if (lgcorr) then
        lemom = F
        lpade = 0
        lnos = F
        lpdos = F
        ldmat = 0
        lrhos = F
        lotf = F
      endif
C ... What of output rho or dos to get (see pgidos)
      moddos = 0
      if (lnos) moddos = 2                     ! Compute dosi
C ... No Pade approximant for this mesh
      if (lpdos .and. ldos >= 4 .and. procid == master)
     .  print *,'GFASA: make m-resolved DOS'
      if (lpdos .and. mod(nint(semsh(2)),100) < 10)
     .  moddos = moddos + 10
      if (lpade > 0)   moddos = moddos + 20   ! accumulate diagonal part of g
      if (ldos >= 4)    moddos = moddos + 40   ! decompose dos into both l and m
      if (lemom) moddos = moddos + 100         ! compute ASA moments of the DOS
      if (lgcorr) moddos = moddos + 200        ! compute Lloyd corrections (gfidos2 only)
      if (ldmat /= 0) moddos = moddos + 1000*ldmat ! Site-diagonal density-matrix
      if (lrhos) moddos = moddos + 4000        ! compute spin density matrix
C     Switches for gfidos2
      moddos2 = 0
      if (ldlm /= 0) then
        moddos2 = moddos
        if (.not. lso) moddos2 = moddos2 + 20000 ! gfidos2 processes CPA sites only
        if (lso .and. lemom) moddos2 = moddos2 + 400 ! gfidos2 calculates orbital moments
        if (mod(moddos2,10) == 2) moddos2 = moddos2 - 1 ! Add to DOS rather than zero it out first
        if (lrel == 2 .and. lemom) moddos2 = moddos2 + 100000
      endif
C ... Make bare ASA response function, q=0
      if (allocated(dgd)) deallocate(dgd,dgdk)
      if (iscr1 /= 0) then
        if (hdim /= 0) call rx('SCR not implemented with HIGH '//
     .    'channels ... replace IDXDN=3 with IDXDN=2')
        nRlc = nnrl(1,1,nbas,s_ham%iprmb,lidim)
        allocate(dgd(nRlc*nRlc*nsp*nspc)) !; call dpzero(dgd,2*nRlc*nRlc*nsp*nspc)
      else
        allocate(dgd(1))
      endif
C ... Make bare ASA response function, all irreducible q
      if (isx1 /= 0) then
        if (hdim /= 0) call rx('SX not implemented with HIGH '//
     .    'channels ... replace IDXDN=3 with IDXDN=2')
        nRlc = nnrl(1,1,nbas,s_ham%iprmb,lidim)
        allocate(dgdk(nRlc*nRlc*nsp*nspc*nkp)) !; call dpzero(dgdk,2*nRlc*nRlc*nsp*nspc*nkp)
      else
        allocate(dgdk(1))
      endif

C --- Printout at entry ---
      if (iprint() >= 30 .and. procid == master) then
        if (ldomg) then
          call awrit1('%x%N --- GFASA: Omega at shifted z, dz=%;7d',
     .      outs,100,stdo,dz)
        elseif (lgcorr) then
          call awrit0('%x%N --- GFASA: Correction from dO/dz',outs,100,
     .      stdo)
        elseif (lspec) then
          call awrit0('%x%N --- GFASA: Spectral function',outs,100,stdo)
        else
          call awrit4('%x%N --- GFASA ('//
     .      '%?#n#gamma,##'//
     .      '%?#n#%b-bar,##'//
     .      '%?#n#%b through S^gam,##'//
     .      '%?#n#frzvc,##'//
     .      '%b%-1j%?#(p>12)#)##',
     .      outs,100,0,bitand(lham,128),bitand(lasa,512),bitand(lasa,1024),isw(lfrzvc))
          call awrit8('%a : make'//
     .      '%?;n; qnu,;;'//
     .      '%?;n; SX-sigma,;;'//
     .      '%?;n; P0(0),;;'//
     .      '%?;n; spin-density matrix,;;'//
     .      '%?;n==1; site-density-matrix,;;%-1j'//
     .      '%?;n==2; density-matrix,;;'//
     .      '%?;n; idos,;;'//
     .      '%?;n; pdos,;;'//
     .      '%?;n; omega only,;;%a%b ',outs,100,-stdo,
     .      isw(lemom),isx1,iscr1,lrhos,ldmat,isw(lnos),isw(lpdos),
     .      isw(omonly))
        endif
      endif
      call info5(30,0,0,' Energy integration : pade=%i%?#n==1#  qtolp=%g#%j#',lpade,isw(qtolp>0),qtolp,4,5)

      if (ldomg .or. lgcorr) then
        lagain = 0
        lmixom = F
        goto 101
      endif

C ... Starting Euler angles
      if (bittst(lncol,1)) then
        neul = s_ham%neula
        if (iprt(2) /= iprt(1)) then
          call pshpr(iprt(2))
        else
          call pshpr(iprint()-0)
        endif
C       call setpr(50)
        call amagnc(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,s_ham%eula,neul,
     .    1,amag,aamom,xx)
        if (bittst(lncol,8)) then
          call bdotsr(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,s_ham%bdots,
     .      lhdim,s_ham%iprmb,1,dum)
C         ehterm(5) = dum(1)
          eterms(17) = dum(1)
        endif
C       External field not implemented in lmgf
C        if (neul > 1 .or. bittst(lncol,8)) then
C          if (.not. bittst(lncol,8)) nbf = 99
C          call bsrhos(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,s_pot%pp,
C     .      s_pot%sop,s_ham%eula,neul,w(obxc),w(omagf),nbf,
C     .      lhdim,s_ham%iprmb,1,dum)
C          eterms(17) = dum(2)
C        endif

        call poppr
      endif

C ... File read neighbor table and real-space strux
      nptab = s_str%npr(nbas+1)

C ... Make is,ifac,qb,qlat
      do  8  i = 1, 3
    8 llshft(i) = lshft(i) /= 0
      call pshpr(1)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr

C ... gflst for crystal GF.
C     gflst is needed when the GF is partitioned into subblocks.
C     It is not meaningful here, but is retained for compatibility w/
C     related programs that require it, e.g. PGF, order-N approaches
      allocate(gflst(3)); call ivset(gflst,1,3,1)
      glstp => gflst
C ... Arrays to adjust potential shifts, for charge neutrality
C     vshfn = current potential shift adjustment:
C     This quantity is estimated

      if (allocated(vshfn)) deallocate(vshfn); allocate(vshfn(glstp(1)))
      if (allocated(vhldp)) deallocate(vhldp); allocate(vhldp(12*gflst(1)))
      if (allocated(vsint)) deallocate(vsint); allocate(vsint(14))
C ... Entire crystal treated as one block
      nbl = 1

C ... Allocate and fill angle and weight arrays for CPA
      if (allocated(ncomp)) deallocate(ncomp)
      allocate(ncomp(nbas)); call iinit(ncomp,nbas)
      if (allocated(dlmcl)) deallocate(dlmcl)
      allocate(dlmcl(nbas))
      call sitepack(s_site,1,nbas,'ncomp',1,ncomp,xx)
      call sitepack(s_site,1,nbas,'dlmcl',1,dlmcl,xx)

C ... Initialize exchange constants J0 for all CPA sites
      do  ib = 1, nbas
        if (ncomp(ib) < 2) cycle
        call dvset(s_site(ib)%j0,1,ncomp(ib),0d0)
      enddo

C ... Initialize exchange constants tau for all CPA sites
      do  ib = 1, nbas
        if (ncomp(ib) < 2) cycle
        call dvset(s_site(ib)%tau,1,ncomp(ib),0d0)
      enddo

C ... Allocate room for partial dos and diagonal GF
      npados = lhdimx
      npdos  = nlo*nbas
      if (lpdos .and. mod(mod(moddos/10,10),2) == 1) then
        deallocate(pdos); allocate(pdos(nzp*npdos*nsp))
        call dpzero(pdos,nzp*npdos*nsp)
        do  ib = 1, nbas
c         call rx('check gfasa for partial dos')
C         if (ncomp(ib) < 2) cycle
          call ptr_site(s_site,8+1,'pdos',ib,nzp,nlo*ncomp(ib)*nsp,xx)
        enddo
      else
        forall (ib = 1:nbas) s_site(ib)%pdos => xxp
      endif
      if (lpade > 0) then
        deallocate(gd); allocate(gd(nzp,nspc,nsp,npados))
      endif

C --- Re-entry for internal self-consistency ---
      lagain = 0; lmixom = F; ninternal = 0; call dpzero(vsint,14)
  101 continue

      if (.not. allocated(dvldu)) allocate(dvldu(lhdim*nsp))
      call dpzero(dvldu,lhdim*nsp)
      if (nlibu /= 0) then
        if (.not. allocated(idu)) allocate(idu(4*nbas))
        call getidu(nbas,s_spec,s_site,idu)
        call u2pph(1,nbas,lmaxu,nsp,nlibu,idu,vorb,
     .    s_ham%iprmb,ldim,lhdim,1d0,dvldu)
      else
        if (.not. allocated(idu)) allocate(idu(1))
      endif

C --- Initialize arrays in preparation for energy, k loop ---
      if (lnos) call dpzero(dosi,5*2)
      if (lnos) call dpzero(dosir,5*2)
      call dpzero(dosn,6)

      if (lemom) then
        call dpzero(s_pot%qnu,3*nlspc)
        if (lrel == 2) then
          allocate(qnusav(nlspcr))
!          if (.not. allocated(qlm)) allocate(qlm(0:2,nl*nl,nsp*nsp))
!          call dpzero(qlm,size(qlm))
          call dcopy(nlspcr,s_pot%qnur(4,1),4,qnusav,1)
          call dpzero(s_pot%qnur,4*nlspcr)
          call dcopy(nlspcr,qnusav,1,s_pot%qnur(4,1),4) ! Preserve linearization energy
          deallocate(qnusav)
          allocate(qnusav(3*nlspc))
        endif
        if (ldlm /= 0) call dpzero(s_pot%mxy,nlspc)
      endif
      if (lgcorr) then
        call ptr_pot(s_pot,8+1,'qcorr',2*nlspc,0,xx)
      else
        call ptr_pot(s_pot,8+1,'qcorr',1,0,xx)
C        allocate(s_pot%qcorr(1))
      endif
C     if (lgcorr) call dpzero(s_pot%qcorr,2*nlspc)
      if (lpade > 0) call dpzero(gd,nzp*npados*nsp*nspc*2)
      if (lrhos) call dpzero(rhos,size(rhos))
      if (iscr1 == 1) call dpzero(dgd,nRlc*nRlc*nsp*nspc*2)
      if (isx1 == 1) call dpzero(dgdk,nRlc*nRlc*nsp*nspc*nkp*2)
      call dpzero(vshfn,gflst); call dpzero(vhldp,12*gflst)

C      print *, '!! read potential functions from file'
C      i = fopna('pfrel',-1,0)
C      allocate(pfall(nzp,lhdim,4,2,2))
C      pfall = 0
C      do  izp = 1, nzp
C      k=0
C      do  j = 1, 999999
C        read(i,*) n1,l,m,n4,n5,n6
C        if (l==0 .and. m==0 .and. n4==0 .and.
C     .    n5==1 .and. n6==1) k=k+1
C        if (k == 5) then
C          backspace i
C          goto 876
C        endif
C        read(i,*) zxx
CC       down-down block
C        lmr = nl*nl*(n1-1) + (l+1)**2 - (2*l+1) + (l+m) + 1
C        if (n5==1 .and. n6==1) pfall(izp,lmr,k,1,1) = zxx
C        if (n5==2 .and. n6==2) pfall(izp,lmr,k,2,2) = zxx
C        if (n5==1 .and. n6==2) pfall(izp,lmr,k,2,1) = zxx
C        if (n5==2 .and. n6==1) pfall(izp,lmr,k,1,2) = zxx
C      enddo
C  876 continue
C      enddo
C      call defcc(opfall,nbas*nzp*lhdim*4*2*2)
C      call zcopy(nbas*nzp*lhdim*4*2*2,pfall,1,w(opfall),1)
C      nzz = nzp
C      print *, '!!'
C      deallocate(pfall)

C --- Initializations for spectral function, including new k-point set
      if (lspec) then
        do  ib = 1, nbas
          norb = s_site(ib)%norb
          allocate(s_site(ib)%alr(norb,nspc,norb,nsp))
          allocate(s_site(ib)%amr(norb,nspc,norb,nsp))
          allocate(s_site(ib)%amrt(norb,nspc,norb,nsp))
        enddo
        if (lspecrl) then
          modsf = 2
        elseif (lspecr) then
          modsf = 1
        else
          modsf = 0
        endif
        if (cmdopt('--band',6,0,strn)) then
          plbopt = strn(7:)
          call info0(30,1,0,' ... use bands from file for specfun')
          nkp = 0 ; i = 1 ; nfbn = 0
          if (procid == master) then
            iopq = 2
            call suqlst(s_lat,plbopt,iopq,1,efermsf,i,xx,nfbn,i,nkp,qp,onesp)
          endif
          call mpibc1(nkp,1,2,.false.,'gfasa','nkp')
          if (nkp <= 0) call rx0('gfasa: no k-points specified')
          call ptr_bz(s_bz,1,'qp',3,nkp,xx)
C    ...  Generate new q-points for spectral function calculation
C    ...  Use i2 in place of nkp to preserve nkp
          if (procid == master) then
            call pshpr(1) ; i = 1 ; iq = 0 ; i2 = 0 ; iopq = 1
            do
              call suqlst(s_lat,plbopt,iopq,1,efermsf,i,xx,nfbn,i,i2,qp,onesp)
              if (i2 <= 0) exit
              do  i1 = 1, i2
                iq = iq+1
                call suqlst(s_lat,plbopt,iopq,0,efermsf,i,xx,nfbn,i,i2,qp,onesp)
                call dpscop(qp,s_bz%qp,3,1,3*iq-2,1d0)
              enddo
            enddo
            call poppr
          endif
          call mpibc1(s_bz%qp,3*nkp,4,.false.,'gfasa','qp')
          if (procid == master .and. iprint() > 50) then
            call info2(50,1,1,
     .        ' Calculating spectral function for %i q-points:',nkp,i)
            do  iq = 1, nkp
              write(6,902) iq,s_bz%qp(1:3,iq)
 902          format(1x,i4,3f8.4)
            enddo
          endif
        else
          modsf = modsf + 10
          if (procid == master)
     .      call info0(0,1,0,' ... use BZ mesh for specfun')
        endif
        if (lcon) then
          if (procid == master)
     .      call info0(0,0,0,'     - integrating SF over contour')
          modsf = modsf + 100
          allocate(sfint(2,nkp))
          sfint = 0
        else
          allocate(sfint(1,1))
        endif
        allocate(spf(2,nkp,nzp))
        if (lspecr) then
          allocate(spfr(2,nbas,nkp,nzp))
        else
          allocate(spfr(1,1,1,nzp))
        endif
        if (lspecrl) then
          allocate(spfrl(2,nl*nl,nbas,nkp,nzp))
        else
          allocate(spfrl(1,1,1,1,nzp))
        endif
      endif

C --- Setup for loop over energies, spin; k is embedded within gidos ---
C     Allocate space for orbital moment in nspc=2 case
      if (allocated(orbtm)) deallocate(orbtm);
      if (nspc == 2) then
        allocate(orbtm(nlo*2*nclassd)) ! DEC compiler requires
        call dpzero(orbtm,nlo*2*nclassd)                  ! orbtm always allocated
      else
        allocate(orbtm(1))
      endif
      if (lordn == 1) then
        call rx('gfordn not ready')
C        if (ldmat /= 0) call rx('dmat not ready order-N')
C        call suordn(s_ctrl,s_lat,s_spec,s_site,lordn,s_str%iax,
C     .    s_str%npr,ngdim)
C        allocate(gii(ngdim))
C        call dpzero(gii,ngdim*2)
CC        call gfrcls(1,1,ncl,avw,alat,plat,s_lat%pos,zeff,s_lat%cy,s_lat%cg,
CC     .    s_lat%indxcg,s_lat%jcg,rmax,s_ctrl%ips,iaxg,clp,w(ocllst),offcH,
CC     .    ndofH,offH,dpfj,ddpfj,gii)
C        call sugemb(s_ctrl,s_site,s_ordn,fgemb,nbas,plat,s_ham%offH,
C     .    s_lat%pos,gii)
      elseif (lordn == 2) then
        call rx('not ready for lordn=2')
      else
C   ... Density matrix is just a local array except for CPA sites
        i = 1
        if (ldmat == 1) i = nl**4*nbas*nsp*nspc
        if (ldmat == 2) i = lidim**2*nsp*nspc*nspc
        if (nspc == 2) i = i * 4
        if (allocated(dmat)) deallocate(dmat)
        allocate(dmat(i)); dmat = 0
        if (ldmat == 1) then
          do  ib = 1, nbas
            norb = s_site(ib)%norb
            i = 4*norb**2
            if (nspc == 2) i = i * 4
            call ptr_site(s_site,8+1,'dmat',ib,i,ncomp(ib),xx)
          enddo
        endif
        allocate(gii(lx2,3-nspc))
        i = 1 ; if (hdim > 0) i = hdim*nspc*mxorb
        allocate(ghh(i,3-nspc))
      endif

      iz1 = 1
C#ifdefC DEBUG
C      j = 4
C      if (cmdopt('-iz1=',j,0,outs)) then
C        if (.not. a2bin(outs,iz1,2,0,' ',j,-1)) call
C     .    rxs2('GFASA: failed to parse "',outs(1:30),'%a"')
C        print *, 'use iz',iz1
C      endif
CC      call zprm('zp',2,zp,nzp,nzp,1)
CC      call zprm('wz',2,wz,nzp,nzp,1)
C#endif

C --- Re-entry point for finding neutral point, lpade>1 ---
      nzi = nzp
      delef = 0
   10 continue

C --- Start loop over energies ---
      if (leres == 1 .and. (.not.lspec) .and. procid == master) then
        write(stdo,531)
  531   format(t6,'Re z',t17,'Im z',t24,'spin',t35,'dos',t49,'idos')
      endif

C --- Start of loop to generate SX potential ---
      if (isx1 /= 0)  then
        call info0(30,1,0,' ... start of loop to make SX sigma')
        i = lidim*lidim*nsp*nspc*nkp
        allocate(gfdmtx(i)); call dpzero(gfdmtx,2*i)
        print *, 'vbare depends spin for now'
      endif
      do  izp = iz1, nzi

      if (lrect == 0) then
        zpi(1,1) = zp(1,izp)
        zpi(2,1) = zp(2,izp)
        wzi(1,1) = wz(1,izp)
        wzi(2,1) = wz(2,izp)
      else
        call zmscop(0,1,1,nznp,nznp,izp-1,0,0,0,znp,zpi)
      endif
      if (ldomg) zpi(1,1) = zpi(1,1) + dz

C      print *, zpi, ' . zp=?'
C      read(*,*) zpi(1,1),zpi(2,1)

C ... Hamiltonian setup for this energy
      call suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,dvldu,vshft,zpi)
      if (lgcorr) goto 14
      nitomg = 0
C ... Re-entry for CPA self-consistency
  50  continue
      i = 1
      if (ldmat == 1) i = nl**4*nbas*nsp*nspc
      if (ldmat == 2) i = lidim**2*nsp*nspc*nspc
      if (leres == 1) then
        if (lemom) then
          call dpzero(s_pot%qnu,3*nlspc)
          if (lrel == 2) call dpzero(s_pot%qnur,3*nlspcr)
        endif
        if (lrhos) call dpzero(rhos,size(rhos))
        dmat = 0
      endif

C ... Make the coherent potential
      if (ldlm /= 0 .or. lspec) call mkcpa(s_ctrl,s_site,s_pot,s_ham,s_lat,lspec,izp)

C     sqrt(P_alpha-dot)
      if (isx1 /= 0 .and. sxad) then
        if (.not. allocated(sqrtdpp)) allocate(sqrtdpp(lidim*2))
        call dpzero(sqrtdpp,2*lidim*2)
      endif

      do  isp = 1, nsp
      if (nspc == 2 .and. isp == 2) cycle
      if (cmdopt('--onesp',7,0,outs).and.isp==2) call rx('not ready')

C --- Make (on-site) G_00 for this energy, integrated over BZ ---
      if (lordn == 1) then
        call rx('gfordn not ready')
CC       call defdr(osali,nl**4*s_str%npr(nbas+1))
C        allocate(sali(nl**4*s_str%npr(nbas+1)))
C        call gforni(11,nbas,s_ham%offH,s_str%npr,s_str%iax,s_pot%palp,
C     .    nl2,s_str%s,sali)
C
C        efree = s_ordn%efre
C        allocate(lmx(nclass))
C        call spec2class(s_spec,nclass,s_ctrl%ics,'lmxa',1,lmx,xx)
C        call sugfrs(s_ctrl,s_lat,zpi,s_ham%offH,lmx,
C     .    s_ham%iprmb,s_pot%pp,efree,s_str%iax,s_str%npr,s_str%s,
C     .    sali,s_pot%palp,lidimx,fgemb,gii,lemb)
C        call gfordn(s_ctrl,s_lat,s_ordn,zpi,s_ham%offH,lmx,
C     .    s_ham%iprmb,s_pot%pp,efree,s_str%iax,s_str%npr,s_str%s,
C     .    sali,s_pot%palp,lidimx,lemb,gii)
C        deallocate(lmx,sali)

      else

C   ... Have gfibz make GF with only diagonal symmetrized
        ibzsw = 11
C   ... Have gfibz make fully symmetrized GF
        if (ldmat == 2) ibzsw = 12
        if (ldlm /= 0) ibzsw = 12
C       print *, 'ibzsw=12'; ibzsw = 12
C   ... Have gfibz make and save on disk g for irreducible points
C       In this case, gii is not accumulated in gfibz, but later
        if (iscr1 == 1 .or. isx1 == 1) then
          ibzsw = 42
        else
          call dpzero(gii(1,isp),lx2*2)
          call dpzero(ghh(1,isp),hdim*nspc*mxorb*2)
        endif
        if (lspec) then
          ibzsw = 42
          if (nspc == 1) ibzsw = ibzsw + 1000 ! record gfqp1 and gfqp2 to disk
        endif
        ibzsw = ibzsw + 100000*(iand(s_lat%lsym,4)/4)
C        call tgf1kp(sspec,s_str%s,s_str%iax,
C     .    nptab,s_bz%qp)

C        if (cmdopt('-si+l',5,0,outs)) ldim = lidim
C        call troth(nl,nbas,isp,nsp,s_lat%pos,s_ham%offH,s_ham%iprmb,s_lat%istab,
C     .    nk1,nk2,nk3,s_bz%ipq,nkp,s_bz%qp,s_bz%wtkp,ldim,lidim-ldim,plat,
C     .    s_str%s,s_str%iax,nptab,s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,
C     .    qb)

        if (nlibu > 0) then
          call vorbydl(-2,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,s_pot%pp,vorb)
        endif
        if (isx1 /= 0 .and. sxad) then
C         make sqrt(delta-alpha)
          allocate(pph1(5*lhdim*nsp),potp(6*nlspc))
          call getnewpp(nl,nsp,nclass,s_pot%pp,potp)
c...for debugging
cc         call defdr(oo1,nlspc)
cc         call pptrns(1,nl,s_ctrl%ipc,nclass,nsp,xx,nbas,potp,
cc     .       w(oo1))
cc         call rlse(oo1)
          call makpph(nl,nsp,nbas,lhdim,s_ctrl%ipc,s_ham%iprmb,
     .      potp,pph1)
          call mksqdela(nsp,lidim,pph1,sqrtdpp)
          call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,ibzsw,s_lat%pos,
     .      nlibu,lmaxu,idu,ludiag,vorb,nkp,s_bz%qp,s_bz%wtkp,
     .      nk1,nk2,nk3,s_bz%ipq,plat,s_str%s,s_str%iax,nptab,
     .      s_lat%istab,s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,qb,
     .      gii(1,isp),ghh(1,isp))
          deallocate(pph1,potp)
C         call rlse(opph1)
        else
          call gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,ibzsw,s_lat%pos,
     .      nlibu,lmaxu,idu,ludiag,vorb,nkp,s_bz%qp,s_bz%wtkp,
     .      nk1,nk2,nk3,s_bz%ipq,plat,s_str%s,s_str%iax,nptab,
     .      s_lat%istab,s_lat%symgr,s_lat%ag,nsgrp,s_bz%star,ifac,qb,
     .      gii(1,isp),ghh(1,isp))
C   ...   Make omega for CPA sites
c         call gfgii(s_ctrl,s_site,s_ham,gii(1,isp),isp)
C          if (ldlm /= 0 .and. .not. lspec)
c    .      call gfomg(s_ctrl,s_site,s_pot,s_ham,izp,isp)
          if (omonly) cycle
        endif
        if (nlibu > 0) then
          call vorbydl(-1,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,s_pot%pp,vorb)
        endif

C   ... Accumulate linear response perturbation, and gii
        if (iscr1 == 1 .or. isx1 /= 0) then
          if (ldmat /= 0)
     .      call rx('density-matrix not implemented in this branch')
          call dpzero(gii(1,isp),lx2*2)
          i = 0
          if (isx1 /= 0) i = 101
          call gfp0ft(s_ctrl,s_ham,s_pot,s_lat,iscr1,i,lemom,
     .      isp,s_ham%offH,s_ham%iprmb,nbas,0,s_lat%pos,nkp,s_bz%qp,zpi,
     .      wzi,nk1,nk2,nk3,s_bz%ipq,s_bz%star,ifac,qb,nRlc,dgd,
     .      dgdk,xx,xx,gii(1,isp))
          if ((isx1 /= 0) .and. (lagain == 0)) then
            call gfenint(s_ctrl,s_ham,s_pot,nbas,s_ham%offH,nzi,zpi,
     .        wzi,izp,lidim,hdim*nspc*mxorb,nspc,nsp,isp,nkp,
     .        s_bz%qp,plat,gfdmtx)
          endif
        endif

      endif

C     call yprm('gii',2,gii,lidimx*lidimx,lidimx,lidimx,lidimx)

      enddo ! Loop over spins

      if (ldlm /= 0) then
        allocate(gwk(lx2,3-nspc))
        call zcopy(lx2*(3-nspc),gii,1,gwk,1)
        do  isp = 1, 3-nspc
          call ztoyy(gwk(:,isp),lidimx,lidimx,lidimx,lidimx,0,1)
        enddo
        call gfgii(s_ctrl,s_site,s_ham,gwk)
        deallocate(gwk)
        if (.not. lspec) call gfomg(s_ctrl,s_site,s_pot,s_ham,izp)
      endif

C --- On-the-fly mixing of Omega for CPA sites
      if (lotf) then
        call mixomg(s_ctrl,s_site,nspc,nzp,izp,omgmix,rms)
        nitomg = nitomg + 1
        i = 45
        if (mod(izp-1,10) == 0) i = 40
        lccpa = rms > omgtol .and. nitomg < nitmax
        call info5(i,0,0,' CPA iter %i iz %i:  RMS domega %;3g'//
     .    '  tol %;3g  repeat=%l',nitomg,izp,rms,omgtol,lccpa)
C  ...  If not converged to tolerance, repeat CPA pass
        if (lccpa) goto 50
      endif
      if (omonly) goto 13

C --- Make conditionally averaged local g for CPA
   14 if (ldlm /= 0) then
        i = 1000 + 104*isw(lgcorr) + 1
C   ... In fully relativistic case calculate unscaled g instead of G
c       Kirill uncomment this line:
c       if (lrel == 2 .or. lso) i = i + 1
        if (lrel == 2 .or. lso) i = i + 1 - 1000
C       call zprm('gii before mkgint',2,s_site(1)%gii,s_site(1)%norb*2,s_site(1)%norb*2,s_site(1)%norb*2)
        call mkgint(i,s_ctrl,s_site,s_pot,s_ham,izp,wzi,dz,0)
      endif

C --- Extract diagonal part of g for Pade, other properties ---
C     Accumulate some or all of: dosi,pdos,gd,s_pot%qnu,dmat,rhos
      do  isp = 1, nsp
        if (nspc == 2 .and. isp == 2) cycle

C   ... Convert g to proper Green's function
        call pshpr(iprint()+10)
C       Onsite g taken from gfg2gr in the SO and lrel=2 cases
c       Scalar relativistic onsite g's in s_site%gc, for CPA and for spectral functions
        i = 0
        if (lso) i = 2
        if (lso .and. bittst(s_ctrl%lbxc,4)) i = 12
        if (lrel == 2) i = 1
        if (i /= 0) then
          ldmsit = T            ! dmat is stored in s_site for all sites
          call gfg2gr(i,s_site,nbas)
        endif

C   ... Scale global g -> G.  Meaningful for non-CPA blocks only.
C       call yprm('gii before scaling g->G',2,gii,(lidimx)**2,lidimx,lidimx,lidimx)
C       if (nsp == 2 .and. lidim == lhdim) then
C       Alternate: only nc,CPA cases through gfg2gnc. gfg2g can manage all others.
        if (nspc == 2 .and. s_ctrl%nccomp /= 0 .and. lidim == lhdim) then
C         print *, 'scaling through gfg2gnc'
          call sugfg2gnc(lso,lrel,0,1,s_ctrl%nccomp,.true.,i2)
          call gfg2gnc(i2,1,s_pot,s_site,isp,nsp,nspc,s_ham%offH,s_ham%iprmb,
     .                 lidim,lidim,0,(/0,0/),gii(1,isp),nbas)
        else
C         Parts of g connected with CPA sites is returned null, if P is null for those sites
C         print *, 'scaling through gfg2g'
          call gfg2g(s_ham,s_pot,100,lrel,isp,1,nbas,1,nbas,lidim,lidim,gii(1,isp),ghh(1,isp),ierr)
        endif
        call poppr
C       print *, 'gii after scaling g->G',gii(1,:)
C       call yprm('gii after scaling g->G',2,gii(1,isp),(lidimx)**2,lidimx,lidimx,lidimx)

C  ...  Make DOS from G
C       Note: non-CPA sites (add 10000 to moddos to exclude CPA sites)
C       Note: in noncollinear case  dosi accumulates sum of local spins, not global spins.
C             This does not affect the total charge density.
        lswd = IAND(s_ctrl%lham,256)/256
        if (.not. lgcorr .and. .not. lso) then
          offd = 0 ; offpd = 0; i = 0
C         To make the spin density, with global quantization axis, uncomment the following line
C         and the second call to gifdos below.
C         if (nspc == 2) i = mod(moddos,10)  ! Separate calc. of dosi if nspc=2
          call gfidos(nl,nbas,moddos+10000-i,lswd,1,nbas,s_ham%iprmb,0,lidim,0,hdim,lidim,
     .      gii(1,isp),hdim,ghh(1,isp),1d0,zpi,wzi,isp,nsp,nspc,s_pot%pp,dvldu,vshft,
     .      s_ctrl%ipc,s_ctrl%nrc,ncomp,nrhos,izp,nzi,nzi,1,1,1,offd,offpd,
     .      dosi,pdos,gd,s_pot%qnu,dum,dum,orbtm,dmat,rhos)
C          Remake dosi with global quantization axis
C          if (i /= 0) then
CC           Rotate g to global quantization axis
C            call rotspn(30100,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,neul,
C     .        s_ham%qss(4),dum,dum,1,lidim,lidim,lidim,lidim,dum,gii)
CC           call yprm('gii in global axis',2,gii,(lidim*nspc)**2,lidim*nspc,lidim*nspc,lidim*nspc)
CC           Make dosi using global quantization axis
C            if (.not. associated(gdnc,gd) .and. lpade > 0) i = i + 20
C            offpd = 0
C            call gfidos(nl,nbas,10000+i,0,1,nbas,s_ham%iprmb,0,lidim,0,hdim,lidim,
C     .        gii,hdim,ghh,1d0,zpi,wzi,1,nsp,nspc,s_pot%pp,dvldu,vshft,
C     .        s_ctrl%ipc,s_ctrl%nrc,ncomp,nrhos,izp,nzi,nzi,1,1,1,offd,offpd,
C     .        dosi,dum,gdnc,dum,dum,dum,dum,dum,dum)
C          endif
        endif

C  ...  CPA sites (site-diagonal Green's function taken from s_site(ib)%gc)
        if (ldlm /= 0) then
          offd = 0 ; offpd = 0
          do  ib = 1, nbas
            call gfidos2(nl,nbas,moddos2,lswd,ib,ib,s_ham%iprmb,0,lidim,0,
     .        hdim,lidim,hdim,ghh,1d0,zpi,wzi,isp,nsp,nspc,s_pot%pp,
     .        dvldu,vshft,s_ctrl%ipc,s_ctrl%nrc,izp,nzi,nzi,1,1,1,
     .        offd,offpd,dosi,pdos,s_site(ib)%pdos,gd,s_pot%qnu,xx,xx,
     .        orbtm,s_site(ib)%dmat,rhos,ncomp(ib),s_site(ib)%thet,
     .        s_site(ib)%dlmcl,s_site(ib)%cpawt,s_site(ib)%gc,
     .        s_site(ib)%gcorr,s_pot%mxy,s_pot%qcorr,s_site(ib)%norb,
     .        s_pot%sop)
          enddo
        endif
      enddo
      if (lrel == 2) then
        offd = 0 ; offpd = 0
        do  ib = 1, nbas
          ic = s_ctrl%ipc(ib)
          call gfidosr(nl,nbas,102,0,ib,s_ham%iprmb,0,lidim,0,
     .      hdim,lidim,hdim,ghh,1d0,zpi,wzi,isp,nsp,nspc,s_pot%pprel,
     .      dvldu,vshft,s_ctrl%ipc,s_ctrl%nrc,izp,nzi,nzi,1,1,1,
     .      offd,offpd,dosir,pdos,s_site(ib)%pdos,gd,s_pot%qnur(1,ic),xx,
     .      orbtm,s_site(ib)%dmat,rhos,ncomp(ib),s_site(ib)%thet,
     .      s_site(ib)%dlmcl,s_site(ib)%cpawt,s_site(ib)%gc,
     .      s_site(ib)%gcorr,s_pot%mxy,s_pot%qcorr,s_site(ib)%norb,
     .      s_pot%GFr(1,ic))

C          call gfidosr(101,s_site,nl,ic,1,s_ctrl%nrc(ic),zpi,wzi,s_pot%pprel,vshft,s_ctrl%ipc,
C     .      s_ctrl%nrc,s_pot%gfr(1,ic),izp,nzi,nzi,1,dosir,s_pot%qnur(1,ic),q12)
C         print*,'ic = ',ic,'nrc = ',s_ctrl%nrc(ic),'GFR = ',s_pot%gfr(1:10,ic)
C         print*,'ic = ',ic,'nrc = ',s_ctrl%nrc(ic),'qnur(1:10,ic) = ',s_pot%qnur(1:3,ic)
!         print*,'dosir(1:5,1) = ',dosir(1:5,1),'gfr(1:5,ic) = ',s_pot%gfr(1:5,ic)
!         print*,'q12 = ',q12
        enddo
      endif

C --- Printout dos, nos, trapezoidal rule, on the fly
C     dosn(1..3) = (z,dos,idos) at prior energy
      if (leres == 1 .and. .not. lspec .and. iprint() >= 20) then
      do  isp = 1, nsp
        if (izp == iz1) then
          call dpzero(dosn(1,isp),3)
        else
          dosn(3,isp) = dosn(3,isp) + (zpi(1,1)-dosn(1,isp)) * (dosi(1,isp)+dosn(2,isp))/2
        endif
        if (isp < nsp/nspc)  then
          write(stdo,532) zpi(1,1),zpi(2,1), '1', dosi(1,1), dosn(3,1)
  532     format(2f11.6,3x,a1,2f14.5)
        elseif (isp == 2) then
          write(stdo,532) zpi(1,1),zpi(2,1), '2', dosi(1,2), dosn(3,2)
        endif
        if (isp == nsp/nspc) write(stdo,532) zpi(1,1),zpi(2,1), 't', dosi(1,1)+dosi(1,2), dosn(3,1)+dosn(3,2)
        dosn(1,isp) = zpi(1,1)
        dosn(2,isp) = dosi(1,isp)
      enddo
      endif

C ... End loop over energy
   13 call clhamz(s_pot%pf)
      if (iprint() > 30 .and. mod(izp,10) == 1 .and. nzp > 20)
     .    call awrit3('%?;n==1;%N;;%-1j done iz=%i  z = (%d,%d)',' ',
     .    80,stdo,izp,zpi,zpi(2,1))

C ... CPA spectral function
      if (lspec) then
        call cpasf(modsf,nspc,s_site,zpi,wz(1,izp),nl,nbas,lidim,s_ham%iprmb,
     .    nkp,s_bz%qp,sfth,spf(1,1,izp),sfint,spfr(1,1,1,izp),spfrl(1,1,1,1,izp))
      endif

      enddo ! Loop over energy

C --- Write spectral function to disk
      if (lspec) then
        if (procid == master) then
          if (.not. lcon) then
            call iospf(modsf,nl,efermsf,nzp,zp,nkp,s_bz%qp,spf,nbas,spfr,spfrl)
          else
            call iospf2(nkp,s_bz%qp,sfint)
            deallocate(sfint)
          endif
        endif

        deallocate(spf,spfr,spfrl)
        if (nspc == 2) then
          call dfclos(fopn('gfqp'))
        else
          call dfclos(fopn('gfqp1'))
          call dfclos(fopn('gfqp2'))
        endif
        ierr = mpipid(2)
        goto 100
      endif

C --- Pack the Hermitian part of dmat into s_site
      if (ldmat == 1 .and. nspc == 2) call packdmat(ldmsit,nbas,nl,dmat,s_site)

C --- Extract output density for LDA+U blocks
      if (nlibu > 0) call gfdmatu(nbas,nsp,nspc,lmaxu,nl,nlibu,idu,dmat,dmatu)

      if (lemom .and. iprint() > 45) call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
C ... Record the density matrices to file(s)
      if (ldmat /= 0 .and. ldmsv .and. procid == master) then
c     if (ldmat /= 0 .and. procid == master) then
        i = fopna('dmat',-1,0)
        call dmatio(ldmat,lidim,nl2,nbas,nsp,nspc,dmat,dmat,-i)
        call fclr('dmat',i)
        if (ldmat == 1) then
          call pshpr(1)
          do  ib = 1, nbas
            do icomp = 1, ncomp(ib)
              i = 0
              call bin2a('',0,0,1000*ib+icomp,2,0,4,aib,i)
              dosfil = 'dmat'//aib(1:i)
              i = fopna(dosfil,-1,0)
              call dmatio(1,1,s_site(ib)%norb,1,nsp,nspc,
     .          s_site(ib)%dmat(:,icomp),xx,-i)
              call fclr(dosfil,i)
            enddo
          enddo
          call poppr
        endif
      endif

c     if (lpade <= 1) deallocate(dmat,dmatb,gii,ghh)
      if (lpade <= 1) deallocate(dmat,gii,ghh)

      if (omonly) goto 110
      if (leres == 1) goto 100
      if (semsh(2) == 20) goto 90

C --- Save P0(q=0) ---
      if (iscr1 == 1 .and. lrect == 0) then
C   ... Convert RPA response fuction = 1/i delta G
        call zscal(nRlc*nRlc*nsp*nspc,(0d0,-1d0),dgd,1)
C   ... File I/O
        call gfp0io(000,'PSTA',nRlc,nsp,nspc,1,dgd)
      endif

C  --- Make screened-exchange Sigma ---
C#ifdef SX
      if (isx1 /= 0 .and. lrect == 0) then

C   ... Save P0(q)
        call zscal(nRlc*nRlc*nsp*nspc*nkp,(0d0,-1d0),dgdk,1)
        call p0herm(nRlc,nsp,nkp,dgdk)
        call p0c(nRlc,nsp,nkp,dgdk)
C   ... File I/O
        call gfp0io(001,'P0',nRlc,nsp,nspc,nkp,dgdk)

C  ... Big memory allocation (needs be reorganized)
        lpdim = nbas
        if (aintra) lpdim = nRlc
        lp2 = lpdim**2
        allocate(madm(nbas*nbas),v0(lp2),v01(lp2),v02(lp2),p0(lp2*2))
        allocate(wscrl(lp2),wscr(lp2),wscrk(lp2*nkp))
        call dpzero(wscrl,2*lp2); call dpzero(wscr,2*lp2)
        call dpzero(wscrk,2*lp2*nkp); call dpzero(p0,2*lp2*2)

C  ... Make wstat(q) by FFT
        do  iq =1, nkp

          call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C     ... q-dependent bare Coulomb interaction
          call madmtq(0,qp,nbas,s_lat%pos,awald,alat,vol,
     .      s_lat%dlv,nkd,s_lat%qlv,nkq,madm,xx,xx)
C         call yprm('madmq(q)',2,madm,nbas*nbas,nbas,nbas,nbas)
          call vbare(0,nsp,aintra,s_pot%vintr,madm,xx,nbas,nl,lpdim,
     .      s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,v0,xx)
c        call yprm('vbare(q)',2,v0,lpdim*lpdim,lpdim,lpdim,lpdim)
C#ifdef v2 | v3
          if (iq == 1) call dcopy(lp2*2,v0,1,v01,1)
          if (iq == 2) call dcopy(lp2*2,v0,1,v02,1)
C#endif

          call p0k2p0(aintra,nRlc,lpdim,nsp,iq,nkp,dgdk,p0)
C         call yprm('P0(q)',2,p0,lpdim*lpdim,lpdim,lpdim,lpdim)

C     ... Screened potential for this q
          call wstat(lpdim,p0,v0,wscr)
C         call yprm('wstat(q)',2,wscr,lpdim*lpdim,lpdim,lpdim,lpdim)
          call dpscop(wscr,wscrk,lp2*2,1,1+(iq-1)*lp2*2,1d0)

        enddo

C --- w(q = lim_x->0 q1*x) ---
      call dpscop(s_bz%qp,qp2,3,4,1,1d0)
C#ifdef v2 | v3
      call wsloc(wscrk,lpdim,v01,v02,alat**3/vol,
     .  s_bz%qp,nkp,nkabc,s_bz%wtkp,wscrl)
C     call yprm('wstl(q)',2,wscrl,lpdim*lpdim,lpdim,lpdim,lpdim)
      call wsloca(nsp,s_ham%iprmb,lidim,nbas,nl,s_ctrl%ipc,s_pot%qnu,lpdim,wscrl)
C     call yprm('wstla(q)',2,wscrl,lpdim*lpdim,lpdim,lpdim,lpdim)
C      Difference wstat - wloca
C#endif

C#ifdefC v1
CC --- Local part of screened potential  ---
C      do  40  iq = 1, nkp
C        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
CC   ... Remake q-dependent bare Coulomb interaction
C        call madmtq(0,qp,nbas,w(obas),awald,alat,vol,
C     .       s_lat%dlv,nkd,s_lat%qlv,nkr,madm,xx,xx)
C        call vbare(0,nsp,aintra,vintra,madm,xx,nbas,nl,lpdim,
C     .    s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,v0,xx)
CC   ... Screened response for this q, local part of P0 ---
C        call wstat(lpdim,w(op0l),v0,wscr)
C        call dpscop(wscr,wscrl,lp2*2,1,1+(iq-1)*lp2*2,1d0)
C   40 continue
C
CC --- w(q = lim_eps->0 q1*eps) ---
C      call madmtq(1,qp2,nbas,w(obas),awald,alat,vol,
C     .  s_lat%dlv,nkd,s_lat%qlv,nkr,w(omad0),w(omad1),v2)
C      call vbare(1,nsp,aintra,vintra,w(omad0),w(omad1),nbas,nl,lpdim,
C     .  s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,w(ov00),v01)
Cc ... Nonlocal, static response P0 for q1
C      call asxp0(lsharm,aintra,1,isp,plat,nl,nsp,nbas,s_ham%offH,
C     .  s_ham%iprmb,nkp,llshft,s_bz%ipq,s_bz%star,ldim,lpdim,eband,nbmx,
C     .  w(oevec),zdel,w(obas),s_lat%istab,s_lat%symgr,s_lat%ag,efermi,w(oh),w(oh2)
C     .  ,w(oh3),p0)
CC ... local part of q->0 limit of P0
C      call p0q0(qp2,lpdim,p0,w(op1),w(op00),w(op01),w(op02))
CC ... q->0 limit of w
C      call wq0(lpdim,w(op00),w(op01),w(op02),w(ov00),v01,v2,
C     .     w(ow0),w(ow1),w(ow2))
CC ... q->0 limit of w from local part of P
C      call p0q0(qp2,lpdim,w(op0l),w(op0l),w(op00),w(op01),w(op02))
C      call wq0(lpdim,w(op00),w(op01),w(op02),w(ov00),v01,v2,
C     .     w(ow0m),w(ow1m),w(ow2m))
CC --- Difference w(P0) and w(local-P0) ---
C      call wdiff(lpdim,nkp,wscrk,wscrl,w(ow0),w(ow1),w(ow2),
C     .           w(ow0m),w(ow1m),w(ow2m))
C#endif

       call gfwdiff(lpdim,nkp,wscrk,wscrl)

C       rewrite
C       sets the shift to the BZ mesh to calc. sigma
       do  82  i = 1, 3
   82  llshft(i) = lshft(i) /= 0
C      Sets the size of the slice handled by gfgw
       pas = 1

C  ... Calculate screened exchange sigma and save slices on disk
       call gfgw(aintra,s_ctrl,s_ham,s_pot,s_lat,s_bz,llshft,
     .   s_ham%offH,s_ham%iprmb,nbas,lpdim,s_lat%pos,nkp,s_bz%qp,nk1,
     .   nk2,nk3,s_bz%ipq,s_bz%star,ifac,qb,pas,wscrk,wscrl)

C      call rlse(owscrk)
       deallocate(madm,v0,v01,v02,p0,wscrl,wscr,wscrk)

C  ... ensembles back sigma
       call gfsigma(s_ctrl,s_ham,s_pot,s_lat,s_bz,s_ham%offH,
     .   nbas,nkp,s_bz%qp,nk1,nk2,nk3,pas,s_str%s,s_str%iax,nptab,sxad,
     .   tolsig,1)

C ... it enables the next iteration for sigma
       print*,'gfasa tolsig=',tolsig

      endif
C#endif
C ... end of screened exchange

C --- Integrated properties: NOS and band-structure energy ---
      lagain = 0
      if (lnos .and. leres == 0) then
        call info2(25,1,0,' GFASA:  integrated quantities to efermi = %;6d',efermi+delef,0)
      endif
C     Hang on to qnu to estimate Pade correction to qnur from (s_pot%qnu - qnusav)
      if (lrel == 2) call dcopy(3*nlspc,s_pot%qnu,1,qnusav,1)

C ... Accumulate N(vbar) from passes involving actual G
      if (lnos .and. leres == 0) then
        izp = vsint(13)
        xx = vbar
        call pshpr(iprint()-10)
        call rfalsi(xx,dosi(2,1)+dosi(2,2)-s_bz%zval,5d-8,1d-7,5d-8,.1d0,30,vsint(1),izp)
        call poppr
        vsint(13) = izp
        vsint(14) = xx
      endif

C     Include contribution from <dot|dot>_12
C     This way allows nonorthogonal (dot)_1 (dot)_2 to be included
C      if (lrel == 2) then
C        dosi(2,1) = dosi(2,1) + q12/2
C        dosi(2,2) = dosi(2,2) + q12/2
C      endif
C     Alternatively, remove qnur(2,:,:,1,2) and qnur(2,:,:,2,1)
C     so that sphere part and band part connnect only through qnur
      if (lrel == 2) call qrel2z12(nl,nclspd,s_pot%qnur)

C ... Re-entry point for Pade approximation
   16 continue
      if (lnos .and. leres == 0) then

C     call pshpr(60)
      if (lagain == 1) call pshpr(iprint()-10)
      ipr = iprint()
      zval = s_bz%zval
      call pgqtot(nbl,nsp,glstp,pgplp,dosi,s_ctrl%ipc,z,s_pot%qc,zval,qtot,sumev)

      xx = (dosi(1,1)+dosi(1,nsp))/2
      if (xx < 0) call info2(25,0,0,'%9p(warning) DOS(ef) = %;6d < 0 ',xx,0)

C --- Float potentials to ensure charge neutrality ---
C ... Est. pot. shift from dos, integrated dos at prior points
      call info0(25,0,0,' ')
      if (lmet == 0 .and. qtolp == 0) then ! Do not attempt to update potential shift
        xx = 0
        call info0(21,1,0,' gfasa : met=0 and qtolp=0 ... skip search for neutrality point')
      elseif (abs(qtot)<qtolp) then ! charge tolerance reached; no further interpolation
        call info2(21,1,0,
     .    ' gfasa : deviation from c.n. < qtolp (%;3g) ... end %?;n;Pade interpolation;search for neutrality point;',
     .    qtolp,lagain)
      else

C   ... Estimate correction to integrated dos by Pade approximant
        if (lpade == 1) then

C         Update vshfp(1); vhold(1) keeps input vshfp(1)
          call gfzerq(2,nbl,nsp,pgplp,s_ctrl%ipc,z,s_pot%qc,
     .      glstp,efermi,dosi,xx,0,zval,vshfn,vhldp,rmsdv,rmsdq)
          vbar = vbar + vshfn(1) - vhldp(1)
          vbardc = vbar
          if (vbar > 1d0) call rx0('gfasa - runaway vbar')

          allocate(wk(glstp(1)))
          call dcopy(glstp,vhldp,12,wk,1)
C         ncomp,dlmcl needed to update s_pot%qnu's for CPA classes
          call gippad(000,moddos,pgplp,glstp,nbl,nzp,nzp,zp,wz,nl,nsp,
     .      nspc,s_ham%iprmb,lidim,s_pot%pp,vshft,s_ctrl%ipc,s_ctrl%nrc,
     .      nrhos,ncomp,dlmcl,nzp,nzp,zp,gd,efermi,0d0,xx,wk,
     .      vshfn,dosi,s_pot%qnu,rhos)
          call pgevdc(0,nbl,nsp,dosi,wk(1),vshfn(1))
C         dosi(3,1) = dosi(3,1) + dosi(2,1)*de
C         dosi(3,2) = dosi(3,2) + dosi(2,2)*de
          deallocate(wk)

C   ... Numerically integrate along real axis, using the trapezoidal (or
C       Simpson rule; see below), between the trial fermi level and the
C       true fermi level corresponding to the charge neutrality condition.
C       semsh(7) supplies step size  de  for numerical integration.
C       Integration proceeds along a mesh of points spaced by  de  until
C       the charge neutrality condition is bracketed.  The fermi level
C       is fixed by linearly interpolating between these points.  For the
C       (small) imaginary part needed to avoid poles, we take Im zp(nzp).
C
C       For the trapezoidal rule (or Simpson rule; see below), G at the
C       lower and upper bounds are required.  For each interval excepting
C       the first, the upper bound of the prior interval serves as the
C       lower bound, and the contribution from each lower bound is
C       obtained by calling gippad.  For the first interval, we do not
C       have G at the lower bound, but only G at energy close to to it,
C       namely zp(nzp), does not exactly coincide with the lower bound.
C       Thus we actually integrate from zp(nzp) to zp(nzp)+de,
C       zp(nzp)+2de, etc.  This introduces a small error (which disappears
C       when the trial fermi level approaches the proper one), but avoids
C       pitfalls such as generating slight deviations from charge
C       neutrality.
C
C       Simpson rule integration proceeds as in the trapezoidal case,
C       except that additional midpoint values are estimated by Pade
C       interpolation.
        elseif (lpade == 2) then

C         If first interval, setup
          if (lrect == 0) then
C           Allocate new gd, copy last nzp points of gd
            nznp = min(nint(semsh(8)),nzp)
            gds => gd
            i = npados*nspc*nsp
            allocate(gd(nzp,nspc,nsp,npados))
            call zmscop(0,nznp,i,nzp,nznp,nzp-nznp,0,0,0,gds,gd)
C           Copy corresponding energies
            if (allocated(znp)) deallocate(znp); allocate(znp(nznp))
            call zmscop(0,nznp,1,nzp,nznp,nzp-nznp,0,0,0,zp,znp)

C           call zprm('gd',2,gds,nzp,nzp,4)
C           call zprm('new gd',2,gd,nznp,nznp-1,4)

            lrect = 2
            de = semsh(7) * dsign(1d0,-qtot)
C           zpi(1,1) = efermi
            zpi(1,1) = zp(1,nzp)
C           The 'fudge' in the energy window; see comments, above.
            defzpn = efermi - zp(1,nzp)
            zpi(2,1) = zp(2,nzp)
            zpi(2,2) = zp(2,nzp)
            zpi(2,3) = zp(2,nzp)
            wzi(2,1) = 0
            wzi(2,2) = 0
            wzi(2,3) = 0

C           Save charge at this energy to know when integral bracketed
            qtotl = qtot

C           Turn off extras for now
            iscr1 = 0
            ldmat = 0

C         Correction to trapezoidal rule using Pade to est midpoint
C         All integrated quantities are adjusted.
          elseif (lrect == 2) then
            zpi(1,3) = (zpi(1,1)+zpi(1,2))/2
            zpi(2,3) = zpi(2,1)
            wzi(1,1) = de/6 - de/2
            wzi(1,2) = de/6 - de/2
            wzi(1,3) = 4*de/6
            call info0(30,1,0,
     .        ' gfasa:  Simpson correction to trapezoidal rule ...')
            call gippad(002,moddos,pgplp,glstp,nbl,3,3,zpi,wzi,nl,nsp,
     .        nspc,s_ham%iprmb,lidim,s_pot%pp,vshft,s_ctrl%ipc,
     .        s_ctrl%nrc,nrhos,ncomp,dlmcl,nznp,nznp,znp,gd,
     .        efermi,0d0,xx,vshfn,vshfn,dosi,s_pot%qnu,rhos)
            
            call dvset(dosi,5,5,0d0)
            call dvset(dosi,5*nsp,5*nsp,0d0)
            dqsimp = (dosi(2,1)+dosi(2,nsp))/(3-nsp) - (zval+qtot)
            qtot = qtot + dqsimp
            call info5(30,0,0,'%9fde = %;4d  dq(Simp) = %;6d  '//
     .        'qtot(Simp) = %;6d ',de,dqsimp,zval+qtot,0,0)
          endif ! lrect

C     ... If Fermi level is bracketed ...
          if (qtot*qtotl <= 0) then
            call info2(30,1,0,' Fix Ef (bracketed between %;4g and %;4g) '
     .        //'by linear interpolation',
     .        defzpn+min(zpi(1,1),zpi(1,2)),defzpn+max(zpi(1,1),zpi(1,2)))

C           The prior call to gfzerq (i.e. before c.n. point was
C           bracketed, see below) assigned pot shift = de.  The call below
C           corresponds to pot shift = 0, and the call determines
C           the linearly interpolated estimate for the c.n. shift
            call dpzero(vshfn,gflst)

C           Estimate shift to charge-neutrality with linear interpolation
            call gfzerq(2,nbl,nsp,pgplp,s_ctrl%ipc,z,s_pot%qc,
     .        glstp,zpi,dosi,xx,0,zval,vshfn,vhldp,rmsdv,rmsdq)

C           For consistency, we require linear interpolation in dos, s_pot%qnu
C           between last two energy points.  Do this by integrating
C           backwards from de to current point, scaling weights.
            de = - vshfn(1)
C           print *, 'muck up de'
C           de = .02d0
            wzi(1,1) = de/2
            wzi(1,2) = de/2
            nzi = 2
            if (lrect == 2) then
              wzi(1,1) = de/6
              wzi(1,2) = de/6
              wzi(1,3) = 4*de/6
              nzi = 3
            endif

            call gippad(002,moddos,pgplp,glstp,nbl,nzi,nzi,zpi,wzi,nl,
     .        nsp,nspc,s_ham%iprmb,lidim,s_pot%pp,vshft,s_ctrl%ipc,
     .        s_ctrl%nrc,nrhos,ncomp,dlmcl,nznp,nznp,znp,gd,
     .        efermi,0d0,xx,vshfn,vshfn,dosi,s_pot%qnu,rhos)
C          print *, dosi(1:3,1)
C          call shoctl(s_ctrl,s_spec,s_pot,0,stdo)

C           Cumulative shift in Ef for charge neutrality is
C              zpi(1,1) - zp(1,nzp) - current vshfn(1)
C           To keep Ef fixed, transfer -shift to pot; store in vshfn
C           de = vshfn(1) + efermi - zpi(1,1)
            de = zpi(1,1) - zp(1,nzp) + de
            call dvset(vshfn,1,gflst,-de)
            vbar = -de

            call info2(20,1,0,' GFASA:  cum. Ef shift from rect.'
     .        //' int. = %;6d: transfer to pot shift',de,0)

            zval = s_bz%zval
            call pgqtot(nbl,nsp,glstp,pgplp,dosi,s_ctrl%ipc,z,
     .        s_pot%qc,zval,qtot,sumev)

            deallocate(dmat,gii,ghh,gd)
            gd => gds
            lrect = 0
C     ... If Fermi level is NOT bracketed ...
C         Last calculated point becomes 1st point for next interval
C         in trapezoidal rule.  Both first and last point are given
C         weight = de/2 and added to integrated properties.
C         This branch performs two tasks:
C         1. it adds the contribution to the first point (when g for the
C            second point is calculated, gfidos will add the contribution
C            from that point).
C         2. it shuffles the energy mesh around to enable quadratic
C            corrections using Simpson's rule together with Pade
          else

            if (iprint() >= 25) call awrit2('%N gfasa:  '//
     .        'add trapezoidal integration from %,4;6g to %,4;6g',
     .        ' ',80,stdo,zpi(1,1)+defzpn,zpi(1,1)+de+defzpn)

C           Initialization in anticipation of linear interpolation step
C           once Fermi level is bracketed
C           Set this point as first of 2 bracketing points for rfalsi
            call dpzero(vshfn,gflst)
            call dvset(vshfn,1,gflst,de)
            call dpzero(vhldp,12*gflst)
            call gfzerq(2,nbl,nsp,pgplp,s_ctrl%ipc,z,s_pot%qc,
     .        glstp,zpi(1,1)+de+defzpn,dosi,xx,0,zval,vshfn,vhldp,
     .        rmsdv,rmsdq)

C           Trapezoidal weights
            wzi(1,1) = de/2
            wzi(1,2) = de/2

C           Add contribution from last energy point, trapezoidal rule
            call gippad(002,moddos,pgplp,glstp,nbl,1,1,zpi,wzi,nl,nsp,
     .        nspc,s_ham%iprmb,lidim,s_pot%pp,vshft,s_ctrl%ipc,
     .        s_ctrl%nrc,nrhos,ncomp,dlmcl,nznp,nznp,znp,gd,efermi,0d0,xx,xx,
     .        xx,dosi,s_pot%qnu,rhos)
            call dvset(dosi,5,5,0d0)
            call dvset(dosi,5*nsp,5*nsp,0d0)

C           Shuffle gd,znp down one register to make room for next energy
            i = npados*nspc*nsp
            call zmscop(0,nznp-1,i,nznp,nznp,1,0,0,0,gd,gd)
            allocate(wk(2*i)); call dpzero(wk,2*i)
            call zmscop(0,1,i,1,nznp,0,0,nznp-1,0,wk,gd)
            deallocate(wk)
C           call zprm('new gd',2,gd,nznp,nznp,i)
            call zmscop(0,nznp-1,1,nznp,nznp,1,0,0,0,znp,znp)
C           Add de to last energy
            zpi(1,2) = zpi(1,1)
            zpi(1,1) = zpi(1,1) + de
            call zmscop(0,1,1,nznp,nznp,0,0,nznp-1,0,zpi,znp)
C           call zprm('new gd',2,gd,nznp,nznp,4)
C           call zprm('new znp',2,znp,nznp,nznp,1)

            delef = delef + de

C           Get diagonal part of gd at a new point
            iz1 = nznp
            nzi = nznp
            goto 10

          endif                 ! if (qtot*qtotl <= 0)

        elseif (lpade > 0) then
          call rxi('gfasa: bad value for lpade:',lpade)
        endif                   ! differing Pade treatments
      endif                     ! if (abs(qtot)<qtolp)

C ... Restore reduced verbosity
      if (lagain == 1) call poppr

C ... If estimate correction to integrated quantities with Pade
      if (abs(qtot)<qtolp) then ! charge tolerance reached; no further interpolation
        lagain = 0
      elseif (lpade==1) then
        lagain = isw(rmsdv > 1d-7 .and. lpade>0)
        if (lagain == 1) goto 16
        call info0(25,1,0,' Corrections to integrated quantities estimated by Pade interpolation')
        zval = s_bz%zval
        call pgqtot(nbl,nsp,glstp,pgplp,dosi,s_ctrl%ipc,z,s_pot%qc,zval,qtot,sumev)
      endif

C ... Update vshift, decide whether dq is large enough to recalculate G
      if (lpade > 0) then
        xx = vshfn(1)
        call advshp(nbl,0,pgplp,glstp,vshfn,vshft)
      else
        vshfn(1) = 0
      endif

C ... Mix Omega only if the potential shift was within padtol
C     If Omega is not converged to omgtol, redo contour integral
C     similar to the case when padtol is exceeded
      if (ldlm /= 0 .and. .not. lgcorr .and. .not. lspec .and. .not. (lmet == 0 .and. qtolp == 0)) then
        lmixom = abs(vshfn(1)) < padtol .or. lpade == 0
        if (lmixom) then
          call mixomg(s_ctrl,s_site,nspc,nzp,0,omgmix,rms)
          if (rms > omgtol) lagain = 1
        elseif (nangl > 0) then
          call info0(20,0,0,' Shift > padtol, Omega not mixed')
        endif
      endif

      if (procid == master) then
        call awrit3(' gfasa  vshft: %;4g   cumulative: %;4g %?;n; ... remake g;;',' ',
     .    scrwid,stdl,vshfn(1),vbar,isw(lagain))
        if (iprint() >= 20) then
          call awrit2('%N gfasa:  potential shift this iter = '//
     .      '%;6d.  Cumulative shift = %;6d',' ',scrwid,stdo,vshfn(1),vbar)
          if (lagain==1 .and. ipr >= 10) then
            if (lmixom) then
              call awrit2('%9f'//
     .          'Omega RMS (%;2g) exceeded tolerance (%;2g) '//
     .          '... redo GF pass',' ',scrwid,stdo,rms,omgtol)
            else
              call awrit2('%9f'//
     .          'vshift (%;2g) exceeded Pade tolerance (%;2g) '//
     .          '... redo GF pass',' ',scrwid,stdo,vshfn(1),padtol)
            endif
          endif
        endif
      endif

C --- Repeat internal self-consistency loop ---
      if (abs(qtot)<qtol) then
        lagain = 0
      else
        lagain = -isw(abs(vshfn(1)) > padtol)
        if (vsint(13) >= -4 .and. vsint(13) <= -2) then
          call info8(30,0,0,
     .      ' Interpolate shift for true G: v1=%1;4e N1=%1;4e  v2=%1;4e N2=%1;4e'//
     .      '%?;n<-3;  bracket=F;;%-1j'//
     .      '%?;n==-2;  lin=T;;%-1j'//
     .      '%?;n==-3;  qua=T;;'//
     .      '  est=%1;6g',
     .      vsint(2),vsint(5),vsint(1),vsint(4),nint(vsint(13)),vsint(14),0,0)
          if (vsint(13) == -4) growdv = growdv * 1.1d0
          if (vsint(13) /= -4) growdv = 1
        endif
        if (lagain == -1) then
          call info2(20,0,1,' shift (%;2g) exceeded'//
     .      ' Pade tolerance (%;2g) ... redo GF pass',vshfn(1),padtol)
        elseif (lpade > 0) then
          call info5(20,1,1,' shift (%;2g) within Pade'//
     .      ' tolerance (%;2g, tolq=%;2g) ... accept Pade correction to G',vshfn(1),padtol,qtol,4,5)
        endif
      endif

!     if (lagain /0 0 .or. (qtol>0 .and. qtot>qtol)) ! Need estimate of vbar for qtol to work
      if (lagain /= 0) then
        xx = vbar*growdv
        if (vsint(13) >= -4 .and. vsint(13) <= -2) then
          if (vsint(13) == -3 .or. vsint(13) == -2) xx = vsint(14)
          call awrit2(' vbar %;6g (dos) %;6g (rfalsi).  Modify vbar',outs,len(outs),0,vbar,xx)
          call query(trim(outs),4,xx)
        endif
        if (xx /= vbar) then
          call advshp(nbl,0,pgplp,glstp,xx-vbar,vshft)
          vbar = xx
        endif
        if (lrel == 2) deallocate(qnusav)
C       goto 101
      endif

      ltmp = lagain /= 0
      call query('redo gf pass',0,ltmp)
      if (ltmp) goto 101

C ... End of integrated-properties branch fixing charge neutrality
C     endif
      endif

C --- Estimate correction to qnur from Pade estimate to change in qnu ---
      if (lrel == 2) then
        if (nclassd /= nclass) call rx('Pade correction to qnur not implemented for CPA')
        call gippadr(nl,1,nclass,s_pot%qnu,qnusav,s_pot%qnur)
      endif

C     debugging
C      call shoctl(s_ctrl,s_spec,s_pot,1000,stdo)
CC      print *, 'from qnur'
CC      call prmx('1851 qnur',s_pot%qnur(:,1),4,4,size(s_pot%qnur(:,1))/4)
C      s_pot%qnu = 0
C      call qnu2qnur(41,nl,nsp,s_pot%qnu,qlm,s_pot%qnur)
C      call shoctl(s_ctrl,s_spec,s_pot,1000,stdo)
C      stop

C --- Magnetic moments corresponding to output density matrix ---
      if (bittst(lncol,1)) then
        call pshpr(iprint()-0)
C       if (bittst(lncol,16)) then
        call amagnc(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,s_ham%eula,neul,0,
     .    amag,aamom,xx)
c       if (ldlm /= 0) call magcpa(s_ctrl,s_site,nl,s_pot%sop,amag)
        if (lrel == 2) print *, '** magcpa needs to update amag for lrel = 2'
        if (ldmat == 1 .and. lrel /= 2) call magcpa(s_ctrl,s_site,nl,s_pot%sop,amag)
C       endif
        if (lrel == 2 .or. lso) call iorbtm(s_spec,0,0,s_ctrl%ics,nl,nlo,nclassd,nsp,orbtm,s_bz%sopertp,xx)
C       External field not implemented
C        if (bittst(lncol,8)) then
C          call bdotsr(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,
C     .      s_ham%bdots,lihdim,s_ham%iprmb,0,dum)
CC         ehterm(6) = dum(1)
C          eterms(18) = dum(1)
C        endif
C        if (bittst(lncol,8)) then
C          call bdotsr(nbas,nl,s_ctrl%ipc,rhos,nrhos,s_pot%qnu,
C     .      s_ham%bdots,lhdim,s_ham%iprmb,0,dum)
C          eterms(18) = dum(1)
C        endif
       call poppr
      endif

C --- Partial DOS by Pade approximant ---
 90   continue
      if (lpdos .and. mod(mod(moddos/10,10),2) == 0) then
        call info0(2,0,0,' gfasa (warning) p DOS for unshifted v')
        deallocate(pdos); allocate(pdos(ndos*npdos*nsp))
        call dpzero(pdos,ndos*npdos*nsp)
        call pgpdos(ndos,dosw,pgplp,glstp,nbl,nl,nlo,nsp,nspc,nzp,zp,gd,pdos)
      elseif (lpdos) then
        ndos  = nzp
        call dpcopy(semsh(3),dosw,1,2,1d0)
      endif

C --- Dump partial DOS, if generated ---
  100 continue
      if (lpdos .and. procid == master) then
        call info0(1,0,0,' ... saving partial DOS')
        call iodos(3,-fopn('DOS'),pdos,ndos,npdos*nsp,ndos,npdos,
     .    dosw(1),dosw(2),nsp,efermi,1)
        call fclose(fopn('DOS'))
        do  ib = 1, nbas
          if (ncomp(ib) < 2) cycle
          i = 0
          call bin2a('',0,0,ib,2,0,4,aib,i)
          dosfil = 'DOS'//aib(1:i)
          call iodos(3,-fopn(dosfil),s_site(ib)%pdos,ndos,nlo*ncomp(ib)*nsp,
     .      ndos,nlo*ncomp(ib),dosw(1),dosw(2),nsp,efermi,1)
          call fclose(fopn(dosfil))
        enddo
      endif

C --- If only Omega was made, record it to files ---
  110 if (omonly .or. (ldlm /= 0 .and. lpdos)) then
        call mixomg(s_ctrl,s_site,nspc,nzp,0,0d0,rms)
        if (procid == master .and. .not. ldomg)
     .    call awrit3('%x%N Omega %?;n;converged to;;%?;n; RMS:;; %;2g',
     .    ' ',scrwid,stdo,isw(rms<=omgtol),isw(rms>omgtol),rms)
      endif

      if ((omonly.and..not.ldomg) .or. leres==1 .or. semsh(2)==20) call rx0('gfasa')

C --- If not done yet, repeat the loop with shifted energies (for CPA) ---
      if (lsev .and. .not. lgcorr) then
C ...   Omega(z+dz) is not yet ready
        print *,'not done yet, repeat the loop with shifted energies'
        if (.not. ldomg) then
          ldomg = T
          s_ctrl%ldomg = 1
          omonly = T
C ...     Copy Omega to shifted Omega
          print *,'Copy Omega to shifted Omega'
          call cpomg(nbas,nspc,nzp,s_site)
          moddos3 = moddos
C ...   Both Omega(z) and Omega(z+dz) are ready
        else
          ldomg = F
          s_ctrl%ldomg = 0
          omonly = F
          lgcorr = T
          moddos = moddos3
        endif
        goto 6
      endif

      if (lsev .and. lgcorr) moddos = moddos3

      if (ldlm /= 0) call dlmsumev(nclass,nangl,nl,s_pot%qnu,lgcorr,
     .  s_pot%qcorr,s_pot%pp,efermi+delef,s_ham%etrms)

C --- Restore potential parameters to orthogonal rep'n ---
      allocate(oo(nlspc))
      call pptrns(1,nl,s_ctrl%ipc,nclassd,nsp,xx,size(s_ctrl%ipc),s_pot%pp,oo)
      deallocate(oo)
      if (lrel == 2) deallocate(qnusav)

C ... Add double-counting correction to eterms,ehterm
      call info(21,0,0,' gfasa: vbar = %,6d adds d.c. term '//
     .  'vbar * zval = %,7;7d',vbardc,vbardc*zval)
      ehterm(3)  = ehterm(3)  + vbardc * zval
      eterms(17) = eterms(17) + vbardc * zval
      eterms(18) = eterms(18) + vbardc * zval

C     Update constant potential shift and undo vshft
      if (efshft) then  ! Modify energy mesh instead of updating vconst
        efermi = efermi - vbar
        s_bz%semsh(3) = s_bz%semsh(3) - vbar
        s_bz%semsh(4) = s_bz%semsh(4) - vbar
        call emesh(s_bz%semsh,zp,wz)
        call advshp(nbl,0,pgplp,glstp,-vbar,vshft)
C       debugging
C        pause
C        call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp) ! pp's in tb rep
C        vbar = 0; goto 101
      else
        call advshp(nbl,1,pgplp,glstp,vbar,vshft)
      endif
      if (lfrzvc) then
        call info2(31,1,0,' gfasa (exit): suppressed vconst '//
     .    'update from %,6d to %,6d',vconst,vconst+vbar)
      elseif (efshft) then
        call info2(31,1,0,' gfasa (exit): Fermi energy updated from'//
     .    ' %,6d to %,6d',efermi+vbar,efermi)
      else
        call info2(31,1,0,' gfasa (exit): vconst updated from'//
     .    ' %,6d to %,6d',vconst,vconst+vbar)
        vconst = vconst + vbar
      endif
      vbar = 0

C     Remove gfqp file
      if (fhndl('GFQP') >= 0) call dfclos(fopn('GFQP'))
      if (mod(moddos,10) == 2) moddos = moddos-1

C ... Repack d.c. terms from potential shifts or applied field
      s_ham%eterms = eterms

      call tcx('gfasa')
      if (mod(moddos,10) == 2) moddos = moddos-1

      end
