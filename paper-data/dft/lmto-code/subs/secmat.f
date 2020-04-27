C#define BLAS3
      subroutine secmat(s_ctrl,s_spec,s_ham,s_pot,s_lat,s_str,ikp,nkp,
     .  qp,wtkp,isp,lrsig,nevmx,efmax,nlibu,lmaxu,idu,ludiag,vorb,
     .  nev,z,s,s_wk,eb,strx,rhrs,rors)
C- Set up Hamiltonian and Overlap and diagonalize secular matrix.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nl nspin nclass lgen3 lncol lham loptc lrel
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ips ipc rmax lgen3 lqp lasa lham
Cio    Passed to:  secmtn secmt2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxb hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  secmtn
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula qss ldham ndhrs nmto kmto
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula bdots nprs iaxs hrs
Cio    Passed to:  secmtn sstrxq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ppn sop pp pti pprel
Cio    Passed to:  secmtn secmt2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat alat nkd nkq awald vol avw
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos dlv qlv cg jcg indxcg
Cio    Passed to:  secmtn sstrxq
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:alph nkaps iax kaps s sdot nitab adot
Cio    Passed to:  secmtn secmt2
Ci Inputs
Ci   ikp   :k-point label
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci          Needed mainly for nonlocal exchange potential
Ci   qp    :k-point
Ci   wtkp  :k-point weights, including spin degeneracy (bzmesh.f)
Ci         :(only used with --invbl to make rsrs,rors)
Ci   isp   :current spin channel (1 or 2)
Ci   lrsig :0 no nonlocal (screened-exchange) potential
Ci         :1 or 2 add screened-exchange potential sigm
Ci         :       as bloch sum of s_ham%hrs
Ci         :-1     add screened-exchange potential sigm
Ci         :       by reading sigma from disk file 'sigm'
Ci   nevmx :largest number of eigenvectors to find
Ci         :If nevmx=-2, secmat does not return z or eb; instead
Ci         :H is returned in z, s is returned
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   nlibu :(LDA+U) number of U blocks
Ci   lmaxu :(LDA+U) dimensioning parameter for U matrix
Ci   idu   :(LDA+U) idu(l+1)=1 => this l has a nonlocal U matrix
Ci   ludiag:T assume LDA+U potential is diagonal in m
Ci   vorb  :(LDA+U) orbital dependent potential matrices
Co Outputs
Co   Case nevmx is not -2:
Co   nev   :actual number of eigenvectors generated
Co   z     :eigenvectors
Co   s     :Work array for overlap matrix (DESTROYED on exit)
Co   eb    :energy bands; alias eband (untouched if nevmx=
Co   Case nevmx = -2:
Co   z     :hamiltonian h
Co   s     :Overlap matrix s
Co   s_wk  :work arrays kept in memory for use by asaddq
Co         :Substitutes for data kept in file 'tmp'
Co   Case --invbl
Co   strx,rsrs,rors : inverse Bloch transformed S, H and O
Cs Command-line switches
Cs   --invbl : Not documented
Cr Remarks
Cu Updates
Cu   17 Dec 15 Repackaged SO setup to use also with GF code
Cu   21 Sep 13 New vshft added to pph
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   07 Apr 11 Copies work arrays formerly written to tmp file
Cu             to memory
Cu   15 Oct 10 Read strux through rdstrx
Cu   11 Aug 09 Handle nevmx=-2.
Cu      Transform 2nd gen pp's to alpha if given in gamma
Cu   13 Nov 08 (ATP) Inverse Bloch transform
Cu   15 Nov 07 Improved LDA+U; works in noncollinear case
Cu   08 Nov 07 (J. Xu) LDA+U, first cut
Cu   18 Jun 04 (A Chantis) fully relativistic hamiltonian
Cu   10 Oct 03 SX sigma can be spin polarized
Cu   23 Sep 03 SX patterned after GW.  sigm(rs) stored in sham->hrs
Cu   14 Feb 03 Added applied magnetic field
Cu   03 Oct 01 bug fix, 2nd gen case when gam-rep & neglected orbitals
Cu   30 Aug 00 Added NMTO ASA Hamiltonian
Cu   14 Sep 99 rewrite argument list using structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nevmx,nev,ikp,nkp,lrsig
      double precision eb(*),qp(3),wtkp(nkp),z(*),s(*),efmax
C     For LDA+U
      integer nlibu,lmaxu,idu(4,*),ludiag
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_str)::   s_str
C ... Local parameters
      type(dp_wk_vec) :: s_wk(3)
      integer lgen3,ldham(16),nbas,nclass,nl,nsp,neul,ldim,
     .  lidim,lihdim,fopnT,nfilet,ndim,nglob
      double precision dval
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      double precision qss(4),vmtz,cksumf,ckbas,plat(3,3),avw,dglob
      real(8), allocatable :: pph(:),oppov(:)
C     real(8), pointer:: alpha(:)
C     For self-energy sigma
      integer ndhrs
C     For inverse Bloch transform
      double precision strx(*),rhrs(*),rors(*)

      complex(8), target :: ctmp(1)
      integer   , target :: itmp(1)
      real(8)   , target :: rtmp(2)

      call tcn('secmat')

C --- Setup ---
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      nclass = s_ctrl%nclass
      lgen3 = s_ctrl%lgen3
      neul = s_ham%neula
      qss = s_ham%qss
      ldham = s_ham%ldham
      plat = s_lat%plat
      ckbas = cksumf(s_lat%pos,3*nbas)
      avw = dglob('avw',0d0,0)
      ndim = nbas * nglob('mxorb')
C     Local copy of alpha, so that it can be altered
C      allocate(alpha(size(s_str%alph)))
C      call dcopy(size(s_str%alph),s_str%alph,1,alpha,1)

      if (lgen3 /= 0) then
        call secmtn(s_ctrl,s_spec,s_lat,s_ham,s_pot,s_str,nl,nsp,nbas,
     .    s_ctrl%ips,s_ham%iprmb,qss,s_ham%eula,neul,ikp,nkp,ldim,
     .    lidim,lihdim,qp,s_pot%ppn,s_pot%sop,isp,nevmx,efmax,nev,z,eb)
      else
C       nfilet = fopn('TMP')
        nfilet = fopnT('TMP',-1,4,0)
        rewind nfilet
        vmtz = s_pot%vmtz
C       If pp's in gamma repsn, transform them to alpha (check 1st chan)
        if (dabs(dval(s_pot%pp,6)-dval(s_pot%pp,5)) < 1d-10+9999)then
          allocate(oppov(nl*nsp*nclass))
          call pptrns(0,nl,s_ctrl%ipc,nclass,nsp,s_str%alph,
     .      size(s_ctrl%ipc),s_pot%pp,oppov)
          deallocate(oppov)
        endif
C       Vector of pp, ordered by site
        allocate(pph(5*lihdim*nsp))
        call makpph(nl,nsp,nbas,lihdim,s_ctrl%ipc,s_ham%iprmb,s_pot%pp,pph)
        call u2pph(0,nbas,lmaxu,nsp,nlibu,idu,vorb,s_ham%iprmb,ldim,lihdim,1d0,pph)
        if (IAND(s_ctrl%lves,2) /= 0) then
          if (mod(s_ctrl%lrel,10) == 2) call rx('potential shifts not implemented for Dirac case')
          call shftpph(nl,nsp,nbas,lihdim,s_ctrl%ipc,s_ham%iprmb,
     .      s_pot%vshft,s_pot%vshft,.false.,.true.,pph)
        endif

        ndhrs = s_ham%ndhrs
C       Fast version works, but comment out for generality; use secmt2
C        if (IAND(s_ctrl%lham,3) /= 0) then
C          if (lrsig /= 0) call rx('update lrsig for call to 2C')
C          call rx('update lrsig for call to 2C --- no nfilet')
C          call secm2c(s_ctrl,ckbas,plat,nl,nsp,nbas,nclass,s_ctrl%ipc,
C     .      s_ham%iprmb,qss,s_ham%eula,neul,ikp,nkp,ldim,lidim,lihdim,ndim,
C     .      qp,s_pot%pp,s_pot%sop,vmtz,s_ctrl%rmax,avw,isp,s_pot%pti,pph,
C     .      nfilet,nevmx,efmax,nev,lrsig,z,eb)
C        else

        if (.not. associated(s_ham%bdots)) s_ham%bdots => ctmp
        if (.not. associated(s_ham%nprs )) s_ham%nprs  => itmp
        if (.not. associated(s_ham%iaxs )) s_ham%iaxs  => itmp
        if (.not. associated(s_ham%hrs  )) s_ham%hrs   => rtmp

        call secmt2(s_ctrl,s_pot,s_str,plat,nl,nsp,nbas,nclass,
     .    s_ctrl%ipc,nlibu,lmaxu,idu,ludiag,vorb,s_ham%iprmb,qss,
     .    s_ham%eula,neul,s_ham%bdots,ikp,wtkp,nkp,ldim,lidim,lihdim,
     .    ndim,qp,s_pot%pp,s_pot%sop,vmtz,s_ctrl%rmax,avw,isp,
     .    s_pot%pti,pph,nfilet,nevmx,efmax,nev,lrsig,ndhrs,
     .    s_ham%nprs,s_ham%iaxs,s_ham%hrs,s_wk,z,s,eb,strx,rhrs,rors)
        deallocate(pph)
        endif
C      endif
      call tcx('secmat')
      end
      subroutine secmt2(s_ctrl,s_pot,s_str,plat,nl,nsp,nbas,
     .  nclass,ipc,nlibu,lmaxu,idu,ludiag,vorb,indxsh,qss,eula,neul,
     .  bdots,ikp,wtkp,nkp,ldim,lidim,lihdim,ndim,qp,pp,sop,vmtz,wsr,
     .  avw,isp,pti,pph,nfilet,nevmx,efmax,nev,lrsig,ndsig,ntabs,iaxs,
     .  sigrs,s_wk,z,s,eb,strx,rhrs,rors)
C- Hamiltonian and Overlap, 2nd generation lmto
C ----------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: lham lncol loptc
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   nclass:number of inequivalent classes
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   qss   :Parameters defining spin-spiral
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   bdots :(magnetic field . Pauli spin matrices), downfolding order
Ci   ikp   :k-point label
Ci   wtkp  :k-point weights, including spin degeneracy (bzmesh.f)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci          Needed mainly for nonlocal exchange potential
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   ndim  :nl*nl*nbas
Ci   qp    :k-point
Ci   pp    :potential parameters (atomsr.f)
Ci   sop   :spin-orbit parameters (atomsr.f)
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci   wsr   :Wigner-Seitz radius, in a.u. (input; alias rmax)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   isp   :current spin channel (1 or 2)
Ci   pti   :inverse potential functions needed for downfolding
Ci   pph   :potential parameters in downfolding order (makpph.f), alpha rep'n
Ci   nfilet:logical units for temporary file.  File has structure:
Ci         sil (idim>0 and nevmx>0)
Ci         i-wave contr. to h (idim>0 and nspc=2)
Ci         sll, beta repsn (nevmx>0)
Ci   nevmx :max. no. evecs to generate.
Ci         -1 suppresses generation of z
Ci         -2 Do not diagonalize, but return overlap in z,
Ci            allocate oc for hamiltonian and place there
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   s     :work array of same dimension as hamiltonian
Co Outputs:
Ci   ccd:  diagonal matrices for 1- 2- & 3-centre CCOR integrals
Co   eigenvalues and eigenvectors are returned in eb, z
Co   nev:  number of evecs generated.
Co   eb    :energy bands
Co   z     :eigenvectors.
Co         :z is used as a work array, whether or not evecs generated
Cr Remarks
Cr  *Orbitals may be folded down as described in the Kanpur notes or
Cr   they may be discarded before diagonalising. Such orbitals are
Cr   indicated by IDXDN(i) = 2 or 3 respectively. If the automatic
Cr   downfolding switch is set, orbitals are assigned to a set using
Cr   a set of rules determined in ATFOLD.
Cr
Cr   Folding down is about the inverse, unscreened potential function
Cr   at e_nu, and this is the only choice allowed (the choice uncouples
Cr   the moments of the charge density of the lower and intermediate
Cr   sets, and is needed in order to ensure a representation-independent
Cr   construction of the i-wave eigenvectors, 3-centre integrals
Cr   and beta^dot).
Cr
Cr   Downfolding automatically turns on the combined correction.
Cr
Cr   bittst(lham,8) can be used to transform
Cr   structure constants to an arbitrary representation.
Cr
Cr   Hybridisation is turned off when bittst(lham,16) is set
Cr
Cr   Dimensions of pph,eb,z are doubled when spins are coupled.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   14 Feb 03 Added applied magnetic field
Cu   01 Jun 01 bug fix, iwave avg (noncol case) when lihdim<ndim
Cu   14 Sep 99 renamed from original secmat.
Cu   28 Apr 98 uses spherical harmonics when lsph set.
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nbas,nclass,neul,isp,nfilet,ldim,lidim,lihdim,ndim,
     .  nevmx,nev,ipc(nbas),indxsh(ndim),ikp,nkp,lrsig,ludiag
      double precision vmtz,avw,qp(3),efmax,
     .  qss(4),plat(3,3),eb(ldim*2),pp(6,nl,nsp,*),z(*),s(ldim,*),
     .  pti(ndim,nsp),eula(nbas,neul,3),wsr(*),pph(5,lihdim,nsp),
     .  sop(0:nl-1,nsp,nsp,9,*),wtkp(nkp)
C     For LDA+U
      integer nlibu,lmaxu,idu(4,nbas)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     For inverse Bloch transform
      double precision wtqp,strx(*),rhrs(*),rors(*)
C     double complex bdots(2,2,lihdim)
      double precision bdots(2,2,2,lihdim)
C     For SX
      integer ndsig,ntabs(nbas+1),iaxs(*)
      double precision sigrs(2,ndsig,ndsig,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_str)::   s_str
C ... Dynamically allocated local arrays
      integer, allocatable :: iaxb(:)
      real(8), allocatable::  alpha(:),adot(:),acopy(:)
      real(8), allocatable :: wk(:),wk2(:),wk3(:)
      real(8), allocatable :: newa(:),xsidt(:),soph(:),sod(:)
      real(8), allocatable :: ppov(:),gma(:),diawk(:),ccd(:),xsi(:)
      complex(8), allocatable :: h(:),sii(:),sil(:),A(:),ohp1(:)
C ... Local parameters
      type(dp_wk_vec) :: s_wk(3)
      logical ndnfld,ccor,lss,lnc,lso,lbf,lx,F,T,bittst
      logical TMPFILE
      parameter (F=.false., T=.true., TMPFILE=.false.)
      integer bit,i,ifi,i1mach,i2,idim,ii,iprint,j,l2,lbloch,ld2,ldimx,
     .  lham,li,linv,lncol,loptic,lov,n,nl2,nlspc,nsite,ncsw
C     For inverse Bloch transform
      integer nsmax
      parameter (nsmax=5000)
      logical invb,cmdopt
      character*20 outs
      double precision qpq(3),xx
C     For self-energy sigma
      integer hreal,nttabs,mxorb,lrel
      procedure(integer) :: nglob,fopna,lgunit

      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) == bit)

C     print *, '!!'
C     qp(1) = .1d0
C     qp(2) = .2d0
C     qp(3) = .3d0

C --- Setup ---
      invb = cmdopt('--invbl',7,0,outs)
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      loptic = s_ctrl%loptc
      nl2 = nl**2
      idim = lidim - ldim
C     hdim = lihdim - lidim
      nlspc = nl * nsp * nclass
C     Local copy of alpha,adot, so they can be altered
      if (allocated(alpha)) deallocate(alpha)
      allocate(alpha(size(s_str%alph)))
      call dcopy(size(s_str%alph),s_str%alph,1,alpha,1)
      if (allocated(adot)) deallocate(adot)
      allocate(adot(size(s_str%adot)))
      call dcopy(size(s_str%adot),s_str%adot,1,adot,1)

C     Noncollinear switches
      lnc = lncol /= 0
      lss = bittst(lncol,2)
      lso = bittst(lncol,4)
      lbf = bittst(lncol,8)
C     In noncollinear case use average of up+down pot. functions.
      if (lnc) then
        call saiwav(ldim,lidim,1,ndim,pti)
        call saiwav(ldim,lidim,5,lihdim,pph)
      endif

C     Possibly rotate to spherical harmonics when making Bloch sum
      lbloch = 0
      if (bittst(lham,256)) lbloch = 1000

C     Combined correction required if downfolding
      ndnfld = idim == 0
      ccor = .not. ndnfld .or. IAND(s_ctrl%lasa,4) /= 0
      if (.not. ccor  .and.  iprint() >= 30  .and.  ikp == 1)
     . print *, 'SECMAT : Combined Correction switched off'

C     Diagonalize by inverse iteration, or not
      linv = 0
      if ((nevmx > 0 .or. nevmx == -2) .and.
     .    IAND(s_ctrl%lqp,2) /= 0) linv = 1

C ... Sanity checks
      if (.not. ndnfld .and. bittst(lham,128))
     .  call rx('SECMAT: no downfolding in gamma rep')

C     Some dimensioning parameters and memory allocation
      ldimx = ldim
      ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      l2 = ldimx**2
      i2 = idim**2
      li = ldim * idim
      call zinit(z,l2)
      if (ccor) then
        allocate(ccd(3*lihdim)); call dpzero(ccd,3*lihdim)
      else
        allocate(ccd(1))
      endif
      allocate(wk(2*ndim))
      if ( .not.  ndnfld ) then
        allocate(h(l2)); call dpzero(h,2*l2)
      endif
      if (lss) then
        i = 2
      else
        i = 1
      endif

C     Arrays for downfolding
      if (allocated(sii)) deallocate(sii)
      if (allocated(sil)) deallocate(sil)
      if (allocated(A  )) deallocate(A  )

      allocate(sii(max(i2*i,1)))
      allocate(sil(max(li*i,1)))
      allocate(A(max(li*i,1)))

C --- Get screened strux from disc and Bloch-transform them ---
      if (invb) allocate(iaxb(10*nsmax))
      nsite = s_str%npr(nbas+1)
      call dcopy(3,qp,1,qpq,1)
C ... make copy of iax for invbl
      if (invb) then
        call rxx(nsite > nsmax,' Increase nsmax in SECMT2')
        call icopy(10*nsite,s_str%iax,1,iaxb,1)
      endif
      if (lss) then
        qpq(1) = qp(1) + qss(1)/2
        qpq(2) = qp(2) + qss(2)/2
        qpq(3) = qp(3) + qss(3)/2
        call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .    s_str%s,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,s,sil,sii)
C       call yprm('Sll(2)',12,s,ldim*ldim,ldim,ldim,ldim)
        call pvsec2(2*ld2,s)
        call pvsec2(2*li, sil)
        call pvsec2(2*i2, sii)
        qpq(1) = qp(1) - qss(1)/2
        qpq(2) = qp(2) - qss(2)/2
        qpq(3) = qp(3) - qss(3)/2
      endif

      call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .  s_str%s,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,s,sil,sii)

C ... debugging ... put sk into z and exit
C      print *, 'debugging ... put zk into z'
C      call dcopy(ld2*2,s,1,z,1)
C      return

C --- Transform potential parameters to gamma representation ---
      if (bittst(lham,128)) then
        allocate(gma(ndim))
C   ... Put pp in gamma; make vector of gamma's
        allocate(ppov(nlspc))
        i = -3
        if (isp == 2) i = -4
        if (lnc) i = -5
        if (IAND(s_ctrl%lasa,512) /= 0) i = -5
        call pptrns(i,nl,ipc,nclass,nsp,gma,nbas,pp,ppov)
C        call prmx('gamma',gma,ndim,ndim,1)
C        call prmx('alpha',alpha,ndim,ndim,1)
        deallocate(ppov)
        call makpph(nl,nsp,nbas,lihdim,ipc,indxsh,pp,pph)
        call u2pph(0,nbas,lmaxu,nsp,nlibu,idu,vorb,indxsh,ldim,
     .    lihdim,1d0,pph)
C       Put gamma-alpha in gam, gamma in alpha
        call daxpy(ndim,-1d0,alpha,1,gma,1)
        call daxpy(ndim, 1d0,gma,1,alpha,1)
C       call prmx('gamma-alpha',gma,ndim,ndim,1)
C --- Or else transform potential parameters to alpha rep'n ---
      else
        allocate(ppov(nlspc))
        call pptrns(0,nl,ipc,nclass,nsp,alpha,nbas,pp,ppov)
        deallocate(ppov)
      endif

C --- Change representation interactively ---
C#ifdefC IREP
C      if (bittst(lham,8)) then
C        call rxx(bittst(lham,128),'already changed to gamma rep')
C        call rxx(.not. ndnfld,'NEW_REP: unset downfolding switches')
C        if (ikp == 1) then
C            print *, 'SECMAT : Transforming to new representation'
C        endif
C        allocate(newa(ldim))
C        allocate(ppov(nlspc))
C        call newalp(ldim,nl,nbas,ipc,nsp,isp,pp,newa)
C        call pptrns(0,nl,ipc,nclass,nsp,newa,nbas,pp,ppov)
C        deallocate(ppov)
C      endif
C#else
      if (bittst(lham,8)) call rx('uncomment IREP in secmat')
C#endif

C --- S^alpha -> S^beta ---
      if (.not. ndnfld) then

C   ... Make a copy of screening alphas, transform i-waves to beta=pti
        allocate(acopy(ndim))
        call dcopy(ndim,alpha,1,acopy,1)
        call makbet(ndim,ldim,idim,indxsh,pti(1,isp),wk,alpha)

C   ... Transform S^alpha to S^beta for downfolding the i-waves
        call shffle(.true.,ndim,indxsh,acopy,wk)
        if (lss) then
C         Transform s(2)=s(q+qss/2) first
          call pvsec2(2*ld2,s)
          call pvsec2(2*li, sil)
          call pvsec2(2*i2, sii)
          call pvsec2(2*li, A)
          allocate(wk2(ndim))
          call dcopy(ndim,acopy,1,wk2,1)
          call trS(ldim,idim,wk2,pti(1,isp),s,sil,sii,A,wk)
C         call yprm('Sil2 (beta)',02,sil,idim*ldim,idim,idim,ldim)
          deallocate(wk2)
C         Restore s(1)=sdot(q-qss/2) to first array, and transform
          call pvsec2(2*ld2,s)
          call pvsec2(2*li, sil)
          call pvsec2(2*i2, sii)
          call pvsec2(2*li, A)
        endif
        call trS(ldim,idim,acopy,pti(1,isp),s,sil,sii,A,wk)
        if (iprint() >= 110) then
          call yprm('Sll (beta)',02,s,ldim*ldim,ldim,ldim,ldim)
          call yprm('Sil (beta)',02,sil,idim*ldim,idim,idim,ldim)
        endif
      endif

C --- Save S_il in temporary file, if needed ---
      if (nevmx > 0) then
        i = ldim
        if (lss) i = ldim*2
        if (TMPFILE) then
          if (idim > 0) call dpdump(sil,idim*i*2,-nfilet)
        endif
        if (idim > 0) call dcopy(idim*i*2,sil,1,s_wk(1)%p,1)
      endif

C --- Make 3-centre hamiltonian integrals over the i-waves ---
      if (idim > 0) then
        call i3cntr(ldim,idim,lnc,sil,pph(1,1,isp),vmtz,h,wk)
        if (lnc) then
C         SS: 3-center terms for s(q+qss)
          if (lss) then
            call pvsec2(2*li, sil)
            call pvsec2(2*ld2, h)
            call i3cntr(ldim,idim,lnc,sil,pph(1,1,isp),vmtz,h,wk)
            call pvsec2(2*li, sil)
            call pvsec2(2*ld2, h)
          endif
          allocate(wk2(ldim*2))
          call dcopy(ldim,pph(3,1,2),5,wk2,1)
          call dpscop(wk2,wk2,ldim,1,1+ldim,1d0)
          call dcopy(ldim,pph(3,1,1),5,wk2,1)
          ncsw = 2000
          if (lss) ncsw = ncsw + 20000
          call rotspn(ncsw,1,nbas,nbas,nl,indxsh,eula,neul,qss(4),
     .      wk2,wk2,ldim,ldim,ldim,ldim,ldim,h,z)
          if (TMPFILE) then
            call dpdump(z,l2*2,-nfilet)
          endif
          call dcopy(l2*2,z,1,s_wk(2)%p,1)
          deallocate(wk2)
          if (iprint() >= 110)
     .      call yprm('iwaves H',12,z,ldimx**2,ldimx,ldimx,ldimx)
        endif
      endif

C --- Make S-dot(alpha) for the combined correction; store in z ---
      if (ccor) then
        nsite = s_str%npr(nbas+1)
        if (lss) then
          qpq(1) = qp(1) + qss(1)/2
          qpq(2) = qp(2) + qss(2)/2
          qpq(3) = qp(3) + qss(3)/2
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,
     .      z(1+2*ld2),sil,sii)
          call pvsec2(2*li, sil)
          call pvsec2(2*i2, sii)
C         call yprm('Sd(2)',12,z(1+2*ld2),ldim*ldim,ldim,ldim,ldim)
          qpq(1) = qp(1) - qss(1)/2
          qpq(2) = qp(2) - qss(2)/2
          qpq(3) = qp(3) - qss(3)/2
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z,
     .      sil,sii)
        else
          if (iprint() >= 100) print *, '... bloch sum sdot'
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z,sil,sii)
        endif
      endif

C --- Transform sll, sdot to gamma representation ---
      if (bittst(lham,128)) then
        allocate(wk3(2*ld2))
        allocate(xsidt(ndim)); call dpzero(xsidt,ndim)

C   ... SS case: s(q+qss/2) also rotated to gamma
        if (lss) then
C         Use copies of gamma-alpha, xsidot since mksbet overwrites them
          call dpcopy(gma,wk,1,ndim,1d0)
          allocate(wk2(ndim))
          call dpcopy(xsidt,wk2,1,ndim,1d0)
C         swap 1st and 2nd B.T. sll, to avoid more workspace
          call pvsec2(2*ld2,s)
          if (ndim /= ldim) call rx('secmat: need shuffle wk')
          call mksbet(ccor,ldim,wk,wk2,s,z(1+2*ld2),wk3)
          deallocate(wk2)
          if (iprint() >= 110) then
            call yprm('Sll(2) (gamma)',02,s,ldim*ldim,ldim,ldim,ldim)
            if (ccor) call yprm('Sdt(2) (gamma)',02,z(1+2*ld2),
     .        ldim*ldim,ldim,ldim,ldim)
          endif
          call pvsec2(2*ld2,s)
        endif

        call shffle(.true.,ndim,indxsh,gma,wk)
        call shffle(.true.,ndim,indxsh,xsidt,wk)
        call mksbet(ccor,ldim,gma,xsidt,s,z,wk3)
        call shffle(.false.,ndim,indxsh,gma,wk)
        call shffle(.false.,ndim,indxsh,xsidt,wk)

        if (iprint() >= 110) then
          call yprm('Sll (gamma)',02,s,ldim*ldim,ldim,ldim,ldim)
          if (ccor) call yprm('Sdt (gam)',02,z,ldim*ldim,ldim,ldim,ldim)
        endif
        deallocate(wk3,xsidt,gma)
      endif

C --- Change representation interactively ---
C#ifdefC IREP
C      if (bittst(lham,8)) then
C        allocate(xsi(ldim))
C        call dmadd(newa,1,1,1d0,alpha,1,1,-1d0,xsi,1,1,
C     .             ldim,1)
C        allocate(xsidt(ldim)); call dpzero(xsidt,ldim)
C        allocate(wk3(2*ld2))
C        call mksbet(ccor,ldim,xsi,xsidt,s,z,wk3)
C        deallocate(wk3)
C        if (iprint() >= 110 .or. F) then
C          call yprm('Sll (new rep)',02,s,ldim*ldim,ldim,ldim,ldim)
C          if (ccor) call yprm('Sdt',02,z,ldim*ldim,ldim,ldim,ldim)
C        endif
C        deallocate(xsi)
C        call dcopy(ldim,newa,1,alpha,1)
C      endif
C#endif

C ... Save S_ll in temporary file
      if (TMPFILE) then
        if (nevmx > 0) call dpdump(s,ldim*i*2,-nfilet)
      endif
      if (nevmx > 0) call dcopy(ldim*i*2,s,1,s_wk(3)%p,1)

C --- Rotate S-dot(alpha) to S-dot(beta) for downfolding ---
      if (ccor) then
C   ... Diagonal matrices for 1- 2- and 3-centre CCOR integrals
        call makdia(nl,nbas,lihdim,indxsh,lidim,ipc,wsr,avw,alpha,adot,ccd)
        call dscal(3*lihdim,avw**2,ccd,1)

        if (.not. ndnfld) then
C     ... w^-2(beta^dot - alpha^dot) for the i-waves
          call shffle(.true.,ndim,indxsh,adot,wk)
          call mkbdot(ldim,idim,lihdim,lnc,ccd,pph(1,1,isp),avw,adot)

C     ... Transform w^-2 S^alpha^dot to the beta representation
          if (lss) then
C           Transform sdot(2)=s(q+qss/2) first
            call pvsec2(2*li, sil)
            call pvsec2(2*i2, sii)
            call pvsec2(2*li, A)
            allocate(wk2(lidim))
            call dpcopy(adot,wk2,1,lidim,1d0)
            call trSdot(ldim,idim,acopy,wk2,A,
     .        z(1+2*ld2),sil,sii)
            deallocate(wk2)
C           Restore sdot(1)=sdot(q-qss/2) to first array
            call pvsec2(2*li, sil)
            call pvsec2(2*i2, sii)
            call pvsec2(2*li, A)
          endif
          call trSdot(ldim,idim,acopy,adot,A,z,sil,sii)
          deallocate(sii,sil,A)
          deallocate(acopy)
        endif

C   ... Remove w^-2 scaling from Sdot
        call dscal(2*ld2,-avw**2,z,1)
        if (lss) call dscal(2*ld2,-avw**2,z(1+2*ld2),1)
      endif

C ... If h not yet allocated, do so now
      if (.not. allocated(h)) then
        allocate(h(l2)); call dpzero(h,2*l2)
      endif

C --- Noncollinear and/or S-O three-center hamiltonian ---
      if (nsp == 2 .and. lnc) then
C       if (lso .or. lbf) then
        if (lbf) then
          allocate(soph(32*ndim)); call dpzero(soph,32*ndim)
          call mksoph(nl,nsp,nbas,lihdim,ipc,indxsh,s_pot%socscl,sop,soph)
C          allocate(sod(ldim*12))
C          call mksod(ldim,lihdim,ndim,indxsh,pph,soph,sod)
C         call prmx('sod0',sod,ldim,ldim,12)
        endif
        if (lso) then
          allocate(sod(ldim*12))
          call mksodn(nbas,nl,nsp,ldim,ipc,indxsh,pp,sop,s_pot%socscl,sod)
C         call prmx('sod',sod,ldim,ldim,12)
        else
          allocate(sod(1))
        endif

C#ifdefC NC
C        lrel = mod(s_ctrl%lrel,10)
C        if (lrel == 2) then
C          if ( .not. bittst(lham,128)) call rx('for fully relativistic set gamma=T')
C          if (IAND(s_ctrl%lham,3) /= 0) then
C            call hmfr2c(nbas,nl,indxsh,qss,eula,neul,ldim,lihdim,ipc,nsp,s_pot%pprel,s,h,z)
C          else
C            call hmfr3c(nbas,nl,indxsh,qss,eula,neul,ldim,lihdim,ipc,nsp,s_pot%pprel,pp,sop,s,h,z)
C          endif
C
C        else
C#endif
          if (lss.and.lso) call rx('S-O coupling incompatible with SS')
          call hmltnc(0,nbas,nl,indxsh,qss,eula,neul,pph,sod,T,
     .      ccor,lss,lnc,lso,ccd,vmtz,ldim,lihdim,s,h,z,wk)
C     ... Add (1+oh)+ [V(LDA+U)-diagonal part] (1+oh) to the LDA hamiltonian
          if (nlibu > 0 .and. ludiag == 0) then
          allocate(ohp1(l2)); call dpzero(ohp1,2*l2)
          call hmltnc(4,nbas,nl,indxsh,qss,eula,neul,pph,sod,T,
     .      ccor,lss,lnc,lso,ccd,vmtz,ldim,lihdim,s,
     .      ohp1,z,wk)
          call asldau(nbas,ldim,lmaxu,isp,nsp,2,nlibu,idu,
     .      indxsh,ohp1,vorb,h)
          deallocate(ohp1,sod)
          endif
C#ifdefC NC
C        endif
C#endif
          if (allocated(sod)) deallocate(sod)
C        call yprm('H after hmltnc',12,h,ldimx*ldimx,ldimx,ldimx,
C     .    ldimx)

C   ... External magnetic field
        if (lbf) then
          call ztoyy(s,ldimx,ldimx,ldimx,ldimx,0,1)
          call ztoyy(h,ldimx,ldimx,ldimx,ldimx,0,1)
C         call zprm('h',2,h,ldimx,ldimx,ldimx)
          call hmladb(ldim,lihdim,ndim,indxsh,pph,bdots,soph,s,h)
          call ztoyy(h,ldimx,ldimx,ldimx,ldimx,1,0)
        endif

C        call yprm('H after hmltnc',12,h,ldimx*ldimx,ldimx,ldimx,
C     .    ldimx)
        if (lbf) deallocate(soph)

C   ... Add iwaves to hamiltonian
        if (lnc .and. idim > 0) then
          if (TMPFILE) then
            rewind nfilet
            if (nevmx > 0) call dpdump(s,1,nfilet)
            call dpdump(s,l2*2,nfilet)
          endif
          call dcopy(l2*2,s_wk(2)%p,1,s,1)
          call daxpy(l2*2,1d0,s,1,h,1)
        endif

C --- Collinear three-center hamiltonian ---
      else
        if (nlibu > 0) then
        endif
        call hmltns(0,ccor,ccd,vmtz,ldim,lihdim,pph(1,1,isp),s,h,z,wk)
        if (bittst(lham,16)) call remhyb(ldim,nl,nbas,h,z)
C   ... Add LDA+U Hamiltonian
        if (nlibu > 0 .and. ludiag == 0) then
c         call yprm('test H',12,h,ldimx*ldimx,ldimx,ldimx,ldimx)
C         This one merely adds V(LDA+U)  to H
C          call asaddu(3,nbas,ldim,lmaxu,isp,nsp,1,nlibu,idu,
C     .      indxsh,vorb,h,xx)
C         This one adds (1+oh)+ V(LDA+U) (1+oh) to the LDA hamiltonian
          allocate(ohp1(l2)); call dpzero(ohp1,2*l2)
          call hmltns(4,ccor,ccd,vmtz,ldim,lihdim,pph(1,1,isp),
     .      s,ohp1,xx,wk)
          call asldau(nbas,ldim,lmaxu,isp,nsp,1,nlibu,idu,
     .      indxsh,ohp1,vorb,h)
          deallocate(ohp1)
c         call ztoyy(h,ldim,ldim,ldim,ldim,1,0)
        endif
C       call yprm('after LDA+U',12,h,ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

      if (iprint() >= 110 .or. F) then
        call yprm('3-C H',12,h,ldimx*ldimx,ldimx,ldimx,ldimx)
        call yprm('3-C O',12,z,ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

C --- Add sigma potential ---
C ... For now, assume sigm not spin-pol; use j instead of ldimx
      j = ldimx
      j = ldim
      if (lrsig >= 1) then
C       Call to bloch : 104110=perm orb, transpose, no add, c*16
        hreal = 0
        if (lrsig >= 4) hreal = 1
        i = 100000 + 0000 + 40*(1-hreal)
        nttabs = ntabs(nbas+1)
        mxorb = nglob('mxorb')
        call bloch(lbloch+i,qp,nl,plat,mxorb,indxsh,1,nttabs,
     .    iaxs,sigrs,ndsig,isp,nsp,j,j,0,j,0,j,0,s,xx,xx)
      elseif (lrsig == -1) then
        ifi = fopna('SIGM',-1,4)
        read(ifi)
        call dpdump(s,ld2*2,ifi)
      endif
      if (lrsig /= 0) then
C       To guarantee Hermitian
        call dosymm(2,s,j,j)
        if (iprint() >= 110) then
          call yprm('SX sigma',2,s,j*j,j,j,j)
        endif
C   ... For now, since sigma not sp
        if (lso) then
          call pvsec1(ldim,s,h)
        else
          call daxpy(l2*2,-1d0,s,1,h,1)
        endif
        if (iprint() >= 110)
     .    call yprm('3-C H + sx',2,h,ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

C --- Return with hamiltonian in z if nevmx is -2 ---
      if (nevmx == -2) then
        call dcopy(2*l2,z,1,s,1)
        call dcopy(l2*2,h,1,z,1)
        deallocate(ccd,wk,h)
        return
      endif

C --- Inverse Bloch transform of S and H ---
      if (invb) then
        wtqp = wtkp(ikp)
        mxorb = nglob('mxorb')
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,iaxb,strx,
     .              strx,strx,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,s)
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,iaxb,rhrs,
     .              rhrs,rhrs,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,h)
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,iaxb,rors,
     .              rors,rors,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,z)
C        call invbl(T,T,nbas,nl,nsp,ldim,plat,nsite,iaxb,
C     .             indxsh,ikp,nkp,qp,s,h,z,strx,rhrs,rors)
        deallocate(iaxb)
        goto 999
      endif

C --- Eigenvalues and eigenvectors of H ---
C#ifdef BLAS3
      lx = .true.
C#elseC
C      lx = .false.
C#endif
      if (linv /= 0) then
        allocate(diawk(ldimx*11))
      else
        allocate(diawk(ldimx*5))
      endif
      call dcopy(2*l2,z,1,s,1)
      lov = 1
      if (bittst(lham,32)) lov = 2
C#ifdefC SIGMO
C      if (lrsig /= 0) lov = 3
C#endif
      call diagno(ldimx,h,s,diawk,lx,lov,linv,nevmx,efmax,nev,z,eb)
      deallocate(diawk)

C --- Printout ---
      if (iprint() >= 30) then
        j = min(9,ldimx)
        if (iprint() >= 35) j = ldimx
C#ifdefC LINUX_PGI
C        do  ii = 1, 1
C#else
        do ii = 1, 2
C#endif
        call awrit3(' SECMAT:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
        write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        enddo
        if (iprint() >= 36 .and. nev > 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,i1mach(2),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif
      if (iprint() >= 110 .and. nev > 0 .and. lov /= 3) then
        outs = 'evec'
        if (lov == 2) outs = 'ortho evec'
        call yprm(outs,2,z,ldimx*ldimx,ldimx,ldimx,nev)
        call yprm('eval',1,eb,ldimx*1,ldimx,nev,1)
        call query('V<110 to skip matrix printing',-1,0)
      endif

  999 continue
      deallocate(ccd,wk,h)

      end
      subroutine pvsec1(ldim,sigm,h)
      implicit none
      integer ldim,i,j
      double precision sigm(ldim,ldim,2),h(ldim,2,ldim,2,2)

      do j = 1, ldim
        do i = 1, ldim
          h(i,1,j,1,1) = h(i,1,j,1,1) - sigm(i,j,1)
          h(i,2,j,2,1) = h(i,2,j,2,1) - sigm(i,j,1)
          h(i,1,j,1,2) = h(i,1,j,1,2) - sigm(i,j,2)
          h(i,2,j,2,2) = h(i,2,j,2,2) - sigm(i,j,2)
        enddo
      enddo

      end
      subroutine pvsec2(ld,s)
      implicit none
      integer ld
      double precision s(ld,2)

      if (ld == 0) return
      call dswap(ld,s,1,s(1,2),1)
      end

