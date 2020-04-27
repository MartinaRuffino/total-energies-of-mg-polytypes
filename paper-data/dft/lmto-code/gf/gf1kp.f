      subroutine gf1kp(s_site,s_ham,s_pot,mode,lbloch,lrel,nl,nbas,isp,nsp,
     .  nspc,nlibu,lmaxu,vorb,idu,ludiag,lrsig,qp,plat,salp,iax,nptab,
     .  ldg,ldh,gf,ghh)
C- Green's function at one k-point
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham bandw lncol ndhrs qss offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nprs iprmb iaxs hrs eula neula
Cio    Passed to:  gfg2g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  gma gmar palp pf pfr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gma cp papg gmar palp pf dpfr ddpfr dpf ddpf
Cio    Passed to:  gfg2g
Cio  s_site :struct for site-specific data; see structures.h
Ci Inputs
Ci   mode: :1s digit:
Ci           0 require only diagonal part of g
Ci           1 require only only site diagonal g
Ci          >1 require entire g
Ci         10s digit:
Ci           0 return unscaled g
Ci           1 scale g to make G, calling gfg2g
Ci           2 transform g^alp to g^gam
Ci           3 combination of 1+2
Ci           4 transform g to gamma-representation by rotating S^alp -> S^gam
Ci             gma is supplied in s_pot%gma
Ci           5 combination of 1+4
Ci        100s digit
Ci           1 1 generate gf treating iwaves nonperturbatively.
Ci             Except for the permuted ordering, gf is equivalent
Ci             to treating iwaves as lwaves.
Ci       1000s digit
Ci           1 CPA mode : case use matrix cp instead of vector palp
Ci           2 CPA mode : like 1 but no CPA components
Ci  lbloch:1s digit pertains to storage of Bloch summed hamiltonian
Ci          0: s is stored in unpacked form
Ci          1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci          0: gf has real, imaginary separated
Ci             gf = gf(ldg,ldg2,2), with gf(*,*,1..2) = real..imag
Ci          1: gf is returned complex*16 format:
Ci             gf = gf(2,ldg,ldg2), with gf(1..2,*,*) = real..imag
Ci          2: gf has real, imaginary separated by columns
Ci             gf = gf(ldg,2,ldg2), with gf(*,1..2,*) = real..imag
Ci        100s digit
Ci           1 if to add hamiltonian to what is is passed in array gf
Ci           2 Make g (otherwise make -g)
Ci           3 combination of 1+2
Ci      1000s digit 1 if to convert s to spherical harmonics
Ci   lrel  :0 for non-relativistic
Ci          1 for scalar relativistic
Ci          2 for Dirac equation
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   qp    :k-point
Ci   plat  :primitive lattice vectors, in units of alat
Ci   salp  :real-space structure constant or hamiltonian matrix
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nptab :total number of pairs in neighbor and iax (pairc.f)
Ci   ldg   :formal leading dimension of gf; see 10s digit lbloch
Ci   ldh   :formal leading dimension of ghh; see 10s digit lbloch
Co Outputs
Co   gf    :gf for this qp for lower+intermediate blocks, in local spin quantization axis
Cl Local variables
Cl   lsgam :0 S = S^alp
Cl         :1 S = S^gambar (not spin dependent) with S^gambar in (1,1) block
Cl         :2 S = S^gamma, gamma is spin dependent
Cl   lglobz:F => P,S in local quantization axis; P is diagonal in the absence of SO coupling
Cl         :T => P,S in global quantization axis; S is diagonal in any beta repsn
Cr Remarks
Cr   gf1kp returns g in local quantization axis, though at some intermediate points
Cr   internally it may be in the global axis, depending on lglobz.
Cu Updates
Cu   15 Jul 18 Add Spin rotations in the CPA case
Cu   18 Jun 18 Synchronize with updated spherical harmonics
Cu   25 Mar 16 Restructured scaling of g to g^gam or G in noncollinear case
Cu   10 Jan 16 New ability to make (P^bet-S^bet)^1 for bet /= alp by rotating S^alp
Cu   08 Jun 14 (Belashchenko) First cut at spin-coupled CPA case
Cu   05 Jan 13 Patched gamma-repsn for CPA case
Cu   18 Dec 12 Completed migration to f90 pointers
Cu   11 Feb 12 (Belashchenko) Updates for DLM
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Nov 10 enable possiblity for LAPACK matrix inversion
Cu   08 Nov 07 (J. Xu) LDA+U implementation, first cut
Cu   15 Dec 04 added sigma part to GF for sx-sigma selconsistency (T. Sandu)
Cu   18 Mar 03 (A Chantis) first cut at relativistic case
Cu             Altered argument list.
Cu   01 Nov 02 bug fix for noncollinear, gamma-repsn
Cu    9 Mar 00 redesigned the iwaves downfolding, added h-waves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,nptab,isp,nsp,nspc,ldg,ldh,mode,niax,lrel,lrsig
      parameter (niax=10)
      integer iax(niax,nptab)
      double precision plat(3,3),salp(*),qp(3),ghh(*)
C     complex(8), target :: gf(ldg,nspc,ldg,nspc) ! kcplx=1
      real(8), target :: gf(ldg,nspc,ldg,nspc,2)  ! kcplx=0
C     For LDA+U
      integer n0
      parameter (n0=10)
      integer nlibu,lmaxu,idu(4,nbas),ludiag
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C ... For structures
!      include 'structures.h'
      type(str_ham) ::   s_ham
      type(str_pot) ::   s_pot
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      real(8), pointer :: wk(:,:),shl(:)
      complex(8), pointer:: pf(:,:)
      complex(8), allocatable :: pfa(:,:),gmanc(:,:,:),pfnc(:,:,:)
C     complex(8), pointer :: su(:,:,:,:) ! kcplx=1
      real(8), pointer :: su(:,:,:,:,:)  ! kcplx=0
C ... Local parameters
      integer PRTG,hdim,i,icpa,idim,ierr,kcplx,kl,lblch2,lbloch,ld11,ld21,ldg2,
     .  ldgx,ldham(16),ldim,ldimx,ldrg,ldrgn,ldrh,ldrhn,ldrs,lgen3,lhdim,lidim,
     .  lidimx,lncol,lsgam,lsph,modg2g,morder,mxorb,ncsw
      integer offp,ogi,osi,ogh
C     cyinv specifies method used to invert matrix
      character cyinv*1,strn*40
      double precision xx,qss(4)
      parameter (PRTG=100)
      integer pfdim
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lhdim,ldham(3))
      equivalence                 (ldimx,ldham(5)),(lidimx,ldham(6))
      equivalence (pfdim,ldham(9))
      logical lss,lso,ldiwav,lglobz,ltmp
      procedure(logical) :: bittst,cmdopt
      procedure(integer) :: bitand,isw,iprint,nglob,rotspd
      procedure(real(8)) :: dlength
C     For screened exchange
      double precision,allocatable:: sigma(:,:,:)
      integer lblchs
      integer hreal,ndhrs,nptabs

C     call pshpr(110)
      call tcn('gf1kp')
      kcplx = mod(lbloch/10,10)
      mxorb = nglob('mxorb')
      call cplxdm(kcplx,ldg*nspc,ldg*nspc,ld11,ld21,ldrgn,ogi)
      call cplxdm(kcplx,ldh*nspc,mxorb,ld11,ld21,ldrhn,ogh)
      ldgx = ldg*nspc
      ldrg = ldrgn/nspc
      ldrh = ldrhn/nspc
      lsgam = 0  ! Becomes nonzero if S is rotated to gamma
      modg2g = 100*mod(mode/10,10)
      lsph = mod(lbloch/1000,10)
      call lmorder(0,morder,[0],[0])

      call sanrg(.true.,modg2g/100,0,5,'gf1kp:','10s digit mode')
      ldiwav = mod(mode/100,10) /= 0
      cyinv = ' '
C     cyinv = '0'

      ldg2 = ldg
      lgen3 = s_ham%lgen3
      ldham = s_ham%ldham
      kl = s_ham%bandw

      lncol   = s_ham%lncol
      idim = lidim - ldim
      hdim = lhdim - lidim
      lglobz = .true.
      lso = bittst(lncol,4)
      icpa = mod(mode/1000,10)  ! >0 => CPA mode, =2 => CPA mode but no CPA elements

C     debugging : mode = 0 for SO, 10 for FP
C      if (iprint() >= 55) then
C      call pshpr(100)
C       print *, 'P'; call pokepf(10*morder,1,pfdim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,s_pot%pfr,gf)
C       gf = 0 ; call gfradp(10*morder,ldg,s_pot%pfr,gf)
C       print *, 'Pdot'; call pokepf(10*morder,1,pfdim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,s_pot%dpfr,gf)
C       print *, 'Pddot'; call pokepf(10*morder,1,pfdim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,s_pot%ddpfr,gf)
C       print *, 'Pa/Pg'; call pokepf(10*morder,1,pfdim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,s_pot%papg,gf)
C       print *, '(g-a)*Pa/Pg'; call pokepf(10*morder,1,pfdim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,s_pot%gmar,gf)
C      call poppr
C      stop
C      endif

C     In CPA gamma-repsn, transform to gamma repsn through rotation of S
C     User must choose GAMMA=4 or 5 now
C      if (icpa /= 0 .and. mod(modg2g/100,10) == 2) then
C        modg2g = modg2g + 200
C      endif

c     if (lncol /= 0 .and. icpa /= 0)
c    .  print *, 'WARNING (gf1kp): lncol and icpa'

C --- Create -1 * Bloch-sum strux of l+i,l+i block ---
C     100s digit should be 2 to generate -S
      lblch2 = 60 000 + mod(lbloch,10000)
      if (lncol /= 0) then
C       allocate(su(ldg,1,ldg,1))  ! kcplx=1
        allocate(su(ldg,1,ldg,1,2))  ! kcplx=0
      else
C       sigma(q)&S(q)-->S(q)-(1/sqrt(P_alpha-dot))*sigma_alpa*(1/sqrt(P_alpha-dot))
        if (lrsig /= 0) then
          allocate(sigma(ldg,ldg,2))
          hreal = 0
          lblchs = 100000 + 0000 + 40*(1-hreal)
          ndhrs = s_ham%ndhrs
          call getntabidx(nbas,s_ham%nprs,nptabs)
          call bloch(lblchs,qp,nl,plat,mxorb,s_ham%iprmb,1,nptabs,
     .      s_ham%iaxs,s_ham%hrs,ndhrs,isp,nsp,ldim,ldim,0,ldim,0,
     .      ldim,0,sigma,xx,xx)
        endif
        su => gf
      endif
      call bloch(lblch2,qp,nl,plat,nl**2,s_ham%iprmb,1,nptab,iax,
     .    salp,nl**2,1,1,lidim,lidim,0,ldg,idim,ldg2,kl,su,xx,xx)
C     call yprmi('-S(q) gf1kp q=%s,(%3;6d)',qp,0,2,su,ldg*ldg,ldg,lidim,lidim)

C     Rotate -S^alpha to -S^gamma = [-S^alpha^-1 + gamma - alpha]^-1
      if (modg2g/100 >= 4) then
        if (ldim /= lidim) call rx('gamma=4,5 not implemented for i-waves, sorry')
C       call prmx('gamma-alpha',s_pot%gma,lhdim,lidim,nsp)

        allocate(wk(lidim,lidim+1))
C       Noncollinear case when gma^up ne gma^dn
        if (lncol /= 0) then

C       should be re-designed, with lrel=2 and noncollinear in a single branch
        ltmp = lrel == 2 .or. lso
        if  (ltmp) then
          if (dlength(size(s_ham%eula),s_ham%eula,1) > 1d-32)
     .      call rx('gf1kp: gamma mode not implemented with local rotations')
        else
          ltmp = dlength(lidim,s_pot%gma(:,1)-s_pot%gma(:,2),1) > 1d-8
        endif
C       In this block noncollinear -S^gam placed directly into gf.  This block requires kcplx=0.
        if (ltmp) then
          if (kcplx /= 0) call rx('modg2g=4 needs kcplx=0')
          if (bittst(lncol,2)) call rx('modg2g=4 not implemented with SS')
C         Park noncollinear -S into 2x2 spinor gf
C          call dpzero(gf,ldgx**2*2)
C          call cplxdm(kcplx,ldg,ldg,ld11,ld21,ldrs,osi)
C          call ymscop(0,lidim,lidim,ldrs,ldrgn,0,0,0,0,su,osi,gf,ogi)
C          call ymscop(0,lidim,lidim,ldrs,ldrgn,0,0,ldrs,ldg,su,osi,gf,ogi)
C         call yprm('-S in (1,1) block',2,gf,ldgx*ldgx,ldgx,lidim,lidim)
C         call yprm('-S in (2,2) block',2,gf(1,2,1,2,1),ldgx*ldgx,ldgx,lidim,lidim)
          ncsw = 0; if (bittst(lncol,2)) ncsw = ncsw + 20000
          i = rotspd(1)     ! Forward rotation
          call rotspn(ncsw+100*i,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
     .      qss(4),xx,xx,ldg,lidim,lidim,ldg,ldg,su,gf)
C         call yprm('S after noncollinear rotation',2,gf,ogi,ldgx,lidimx,lidimx)
          deallocate(wk); allocate(wk(lidimx,lidimx+1))
          if (lidimx /= ldgx) call rx('gf1kp: dimension mismatch')
          call yqinv('N',gf,ldgx*ldgx,ldgx,2,ldgx,wk,ldgx,ierr) ! -S^alpha^-1
          if (ierr /= 0) call rxi('gf1kp: S is singular, ierr=',ierr)
C         call yprm('-Salp^-1',2,gf,ogi,ldgx,lidimx,lidimx)
          allocate(gmanc(ldg,2,2)) !! Make it complex, merge with rel case
          if (associated(s_pot%gmar)) then ! P^gam couples spins
            call zcopy(size(gmanc),s_pot%gmar,1,gmanc,1)
          elseif (associated(s_pot%gma)) then ! Assume P made as collinear pfun
            call dpzero(gmanc,2*size(gmanc))
            call dcopy(lidim,s_pot%gma(1,1),1,gmanc(1,1,1),2)
            call dcopy(lidim,s_pot%gma(1,2),1,gmanc(1,2,2),2)
          else
            call rx('no gamma-alpha available')
          endif
          call pokepf(1+10*morder,kcplx,lidim,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,gmanc,gf)
          deallocate(gmanc)
C         call yprm('-Salp^-1 + gamma-alpha',2,gf,ldgx*ldgx,ldgx,ldgx,ldgx)
          call yqinv('N',gf,ldgx*ldgx,ldgx,2,ldgx,wk,ldgx,ierr) ! -S^alpha^-1 + gamma - alpha]^-1
C         call yprm('-Sgam',2,gf,ldgx*ldgx,ldgx,ldgx,ldgx)

C         else  ! This branch works in special cases but comment out for now
CC        mch salp -i gma -coll 1 -v2dia -+ -i
C         call yqinv('N',gf,ldgx*ldgx,ldgx,2,lidim,wk,lidim,ierr) ! -S^alpha^-1
CC        call yprm('-Salp^-1',2,gf,ldgx*ldgx,ldgx,lidim,lidim)
C         call subdiag(kcplx,ldgx,lidim,lhdim,0,1,s_pot%gma,gf) ! -S^alpha^-1 + gamma(1) - alpha]
CC        call yprm('-Salp^-1 + gamma-alpha',2,gf,ldgx*ldgx,ldgx,lidim,lidim)
C         call yqinv('N',gf,ldgx*ldgx,ldgx,2,lidim,wk,lidim,ierr) ! -S^alpha^-1 + gamma - alpha]^-1
CC        call yprm('-Sgam (spin 1)',2,gf,ldgx*ldgx,ldgx,lidim,lidim)
C
CC        mch salp -i gma -coll 2 -v2dia -+ -i
C         call yqinv('N',gf(1,2,1,2,1),ldgx*ldgx,ldgx,2,lidim,wk,lidim,ierr) ! -S^alpha^-1
CC        call yprm('-Salp^-1',2,gf(1,2,1,2,1),ldgx*ldgx,ldgx,lidim,lidim)
C         call subdiag(kcplx,ldgx,lidim,lhdim,0,2,s_pot%gma,gf(1,2,1,2,1)) ! -S^alpha^-1 + gamma(2) - alpha]
CC        call yprm('-Salp^-1 + gamma-alpha',2,gf(1,2,1,2,1),ldgx*ldgx,ldgx,lidim,lidim)
C         call yqinv('N',gf(1,2,1,2,1),ldgx*ldgx,ldgx,2,lidim,wk,lidim,ierr) ! -S^alpha^-1 + gamma - alpha]^-1
CC        call yprm('-Sgam (spin 2)',2,gf(1,2,1,2,1),ldgx*ldgx,ldgx,lidim,lidim)
C         endif

          lglobz = .false.      ! S has been rotated to local quantization axis
          lsgam = 2             ! Flags S = S^gamma, gamma is spin dependent
C         call yprm('-Sgam^-1 (nc)',2,gf,ldgx*ldgx,ldgx,ldgx,ldgx)

          goto 10               ! Inversion accomplished ... skip collinear inversion
        endif                   ! Noncollinear S^gamma
        endif

C   ... Turn S -> S^gamma, gamma not spin dependent with S^gamma in (1,1) block
C       mch salp -i gma -coll 1 -v2dia -+ -i
C       call yprm('-Salp',2,su,ldg*ldg,ldg,lidim,lidim)
        call yqinv('N',su,ldg*ldg,ldg,2,lidim,wk,lidim,ierr) ! -S^alpha^-1
C       call yprm('-Salp^-1',2,su,ldg*ldg,ldg,lidim,lidim)
        call subdiag(kcplx,ldg,lidim,lhdim,0,isp,s_pot%gma,su) ! -S^alpha^-1 + gamma - alpha]
C       call yprm('-Salp^-1 + gamma-alpha',2,su,ldg*ldg,ldg,lidim,lidim)
        call yqinv('N',su,ldg*ldg,ldg,2,lidim,wk,lidim,ierr) ! -S^alpha^-1 + gamma - alpha]^-1
C       call yprm('-Sgam',2,su,ldg*ldg,ldg,lidim,lidim)
        lsgam = 1               ! Flags S = S^gamma, gamma not spin dependent with S^gamma in (1,1) block


   10   continue                ! Re-entry for noncollinear S^alpha->S^gamma
        deallocate(wk)
        modg2g = modg2g - 400   ! Strip from modg2g switches handling S^alpha->S^gamma
      endif

C --- 3nd generation NMTO --
      if (lgen3 /= 0) then
        call rx('gf1kp not ready for 3rd generation')

C --- 2nd generation ASA : Make (P-S)^-1 in gf ---
      else
C       Spin average i+h -wave potential functions
        allocate(pfa(pfdim,nsp))
C        if (lrel == 2 .and. lhdim /= ldim) then
C          call gfradpl(kcplx,ldim,idim,ldg,w(opfr),ogi,gf)
C          stop 'not ready'
C        endif
C  ...  NB: In the CPA case gf1kps only affects normal sites
C       and pfa is used only for h-waves (l+i are already in s_pot%cp)
        pf => s_pot%palp; if (lsgam /= 0) pf => s_pot%pf
C        call zprm('palp',2,s_pot%palp,lhdim,lhdim,2)
C        call zprm('pf',2,s_pot%pf,lhdim,lhdim,2)
        if (icpa == 0 .or. icpa == 2) then
          call gf1kps(nsp,ldim,pfdim,isw(.not.ldiwav),pf,pfa) ! pfa in hamiltonian order
        elseif (lrel /= 2 .and. .not. lso) then
          call gf1kps(nsp,lidim,pfdim,isw(.not.ldiwav),pf,pfa) ! pfa in hamiltonian order
        endif
        if (lncol /= 0) then
          lss = bittst(lncol,2)
          call rxx(lss,'gf1kp not set up for spin spiral')
          qss = s_ham%qss
          ncsw = 0
          if (lss) ncsw = ncsw + 20000
          if (lsgam == 2) ncsw = ncsw + 30000
          if (lrel == 2 .or. lso) then
            if (lsph /= 1) call rx('gf1kp: fully relativistic GF requires spherical harmonics')
            if (lhdim /= ldim) call rx('downfolding not ready lrel=2 or SO')
          endif
          if (icpa == 1 .and. mod(modg2g/100,10) == 2)
     .       call rx('gf1kp: noncollinear CPA not ready for gamma-rep')

C     ... Invert noncollinear P-S, brute force
          if (lss .or. lso .or. lrel == 2 .or. icpa /= 0 .or. .not. lglobz) then
C           Rotate spinor part of -S, unless rotated -S is already given
            if (lglobz) then
              i = rotspd(1)     ! Forward rotation
              call rotspn(ncsw+100*i,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
     .          qss(4),xx,xx,ldg,lidim,lidim,ldg,ldg,su,gf)
C             Rotate orbital part of S.  For now, do only in the presence of SO
              ltmp = icpa == 2 .and. (lrel == 2 .or. lso) .and. .not. cmdopt('--spinoronly',12,0,strn)
              if (ltmp) then
                i = 1000 + 100*kcplx + 10*lsph + rotspd(1) ! Forward rotation, orbital part
C               call yprmi('S (rel) q=%s,(%3;6d)',qp,i,2,gf,ogi,ldgx,lidimx,lidimx)
                call rotheu(i,s_ham,nl,1,nbas,nbas,0,lidim,ldg,ldg,gf)
C               call yprmi('S-orb-rot (rel) q=%s,(%3;6d) i=%i',qp,i,2,gf,ogi,ldgx,lidimx,lidimx)
              endif
              lglobz = .false.
            endif
C           call yprm('S after noncollinear rotation',2,gf,ogi,ldgx,lidimx,lidimx)

C       ... Add potential parameters or coherent potential
            if (icpa == 0) then
              if (lso) call rx('gf1kp: Turn on CPA to use SOC')
C             Add P to -S
              offp = 0
              call sblhm1(lbloch,nl**2,1,nbas,s_ham%iprmb,pfa,
     .          offp,ldim,idim*0,ldg*2,idim,ldg2*2,kl,gf,xx,xx)
              offp = lhdim
              call sblhm1(lbloch,nl**2,1,nbas,s_ham%iprmb,pfa,offp,
     .          ldim,idim*0,ldg*2,idim,ldg2*2,kl,gf(1,2,1,2,1),xx,xx)
C           Not true CPA: P still a vector
C           else if (icpa == 2 .and. (lrel == 2 .or. lso) .and. lsgam == 0) then ! also ok since s_pot%cp available
            else if (icpa == 2 .and. (lrel == 2 .or. lso)) then
C             In future this can be merged with branch using s_pot%cp
C             But s_pot%cp is generated in gamma-repsn (see mkcpa)
C             This works if GAMMA=4 or 5 (flagged by lsgam>0), but not if GAMMA=1,2
C             In such a case s_pot%cp should be in alpha.
C             Didn't make the change in case this messes up other uses of s_pot%cp
              pf => s_pot%palp; if (lsgam == 2) pf => s_pot%pfr

C               call zprm('pf-rel',2,pf,lhdim,lhdim,4)
C               call zprm('palp',2,s_pot%palp,lhdim,lhdim,4)
C               call zprm('sr(dotpf)-rel',2,s_pot%dpfr,lhdim,lhdim,4)
C               call zprm('ddotpf-rel',2,s_pot%ddpfr,lhdim,lhdim,4)
C               call gtpffr(w(opfall),nzp,lhdim,izz,3,pf,ldg)
C               call zprm('pf-rel',2,s_pot%palp,lhdim,lhdim,4)
C               call zprm('pf-srel',2,s_pot%palp,lhdim,lhdim,2)
C               call yprm('starting S',2,gf,ogi,ldg,lidim,lidim)
C                debugging
C                call rotpfr(10*mod(lbloch/1000,10),s_ham%eula(2,1:3),nl,
C     .            s_ham%iprmb,0,lidim,nl*nl*(2-1),pfdim,pf,nl**2,1,sloc,
C     .            s_ham%offh)
C                stop

              allocate(pfnc(ldg*nspc,ldg*nspc,1))
              if (lhdim /= ldim) then
                call rx('check remake gfradpl with updated m ordering')
C               pfnc = 0
C               call gfradpl(kcplx+10*morder,ldim,idim,ldg,pf,ogi,gf)
              else
C               call yprmi('-S (rel) q=%s,(%3;6d)',qp,0,2,gf,ogi,ldgx,lidimx,lidimx)
C               Either line should be equivalent ... the latter works with downfolding
C               pfnc = 0
C               call gfradp(i,ldg,pf,gf)  ! old --- should work
                call pokepf(10*morder,kcplx,ldg,ldg,0,lidim,s_ham%iprmb,0,1,nbas,nspc,nsp,pf,pfnc)
              endif
C             Rotate P to global spin quantization axis.  Not needed because S already rotated
C              i = 0000 + 100*kcplx + 10*lsph + rotspd(0) ! Backwards rotation
C              call yprmi('P-unrot (rel) q=%s,(%3;6d)',qp,0,2,pfnc,ogi,ldgx,lidimx,lidimx)
C              call rotheu(i,s_ham,nl,1,nbas,nbas,0,lidim,ldg,ldg,pfnc)
C              call yprmi('P-rot (rel) q=%s,(%3;6d) i=%i',qp,i,2,pfnc,ogi,ldgx,lidimx,lidimx)
C             Rotate S to local spin quantization axis ... already done above
C             i = 0000 + 100*kcplx + 10*lsph + rotspd(1) ! Forward rotation
C             i = i + 2000 ! To skip orbital part ...
C             call yprmi('S (rel) q=%s,(%3;6d)',qp,i,2,gf,ogi,ldgx,lidimx,lidimx)
C             call rotheu(i,s_ham,nl,1,nbas,nbas,0,lidim,ldg,ldg,gf)
C             call yprmi('S-rot (rel) q=%s,(%3;6d) i=%i',qp,i,2,gf,ogi,ldgx,lidimx,lidimx)
C             Add P to -S
              call daxpy(2*size(pfnc),1d0,pfnc,1,gf,1)
C             call yprmi('P-S (loc z,rel) q=%s,(%3;6d)',qp,0,2,gf,ogi,ldgx,lidimx,lidimx)
              deallocate(pfnc)
            else ! CPA branch
              i = mod(modg2g/100,10)
              if ((i==1 .or. i==2) .and. lsgam==0) call rx('gf1kp: s_pot%cp in gamma but require alpha')
C             call yprm('-S before adding s_pot%cp',2,gf,ogi,ldg*nspc,lidimx,lidimx)
              call ztoyy(gf,ldimx,ldimx,ldimx,ldimx,kcplx,1)
              call daxpy(2*(ldimx)**2,1d0,s_pot%cp,1,gf,1)
              call ztoyy(gf,ldimx,ldimx,ldimx,ldimx,1,kcplx)
C             call yprm('P-S with P from s_pot%cp',2,gf,ogi,ldg*nspc,lidimx,lidimx)
            endif

C     ... Fast inversion exploits fact that Pfnc is diagonal in RL
          else
            allocate(pfnc(ldg,nspc,nspc))
            call dpzero(pfnc,2*ldg*nspc**2)
            call dpscop(pfa,pfnc,ldgx,1,1,1d0)
            offp = lhdim
            call dpscop(pfa,pfnc,ldgx,1+2*offp,1+2*3*ldg,1d0)
C           call zprm('pfnc',2,pfnc,ldg,ldg,4)
C           Debugging : show rotspn works for kcplx=0
C           call ztoyy(pfnc,ldg,4,ldg,4,1,0)
C           call yprm('pfnc(0)',2,pfnc,ldg*4,ldg,ldg,4)
C           Rotate to global axis (backward rotation)
            i = rotspd(0)
            call rotspn(230010+100*i,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,xx,
     .        xx,xx,ldg,lidim,lidim,ldg,1,xx,pfnc)
C           call zprm('pfnc',2,pfnc,ldg,ldg,4)
            call cplxdm(kcplx,ldg,ldg,ld11,ld21,ldrs,osi)
C       ... Poke -S into (1,1) and (2,2) blocks of g ... skip if -S is already in g
            if (lsgam /= 2) then
              call dpzero(gf,(ldgx)**2*2)
              call ymscop(0,lidim,lidim,ldrs,ldrgn,0,0,0,0,su,osi,gf,ogi)
              call ymscop(0,lidim,lidim,ldrs,ldrgn,0,0,ldrs,ldg,su,osi,gf,ogi)
            endif

C       ... Make (Pfnc - S)^-1, where Pfnc is in global axis
            call gf1kpn(mod(mode,10),kcplx,nbas,s_ham%iprmb,ldim,lidim,ldg,ldrs,osi,ogi,pfnc,su,gf)
C           call yprm('gnc after ll block, fixed quantization axis',2,gf,ogi,ldgx,lidimx,lidimx)
C       ... Rotate ll block to local z quantization axis
C            call rotspn(30000,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
C     .        qss(4),xx,xx,ldg,ldim,ldim,ldg,ldg,su,gf)
CC            call yprm('gnc (local quant axis || z)',2,gf,ogi,ldgx,lidimx,lidimx)
C            call rotspn(30100,1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
C     .        qss(4),xx,xx,ldg,ldim,ldim,ldg,ldg,su,gf)
C            call yprm('gnc (fixed z again)',2,gf,ogi,ldgx,lidimx,lidimx)
C           print *, 'return with inv P-S'; return
C           g is already directly obtained; skip past inversion
            deallocate(pfnc)
            offp = 0
            goto 100
          endif
        else ! Collinear case
          offp = (isp-1)*lhdim
C     ... Generate (P-S) for l+i blocks
          if (icpa == 0 .or. icpa == 2) then
            call sblhm1(lbloch,nl**2,1,nbas,s_ham%iprmb,pfa,offp,
     .      lidim,idim*0,ldg,idim,ldg2,kl,gf,xx,xx)
          else
C           CPA case: cp is a matrix; convert S to gamma here if req'd ... already accomplished now
C            if (mod(modg2g/100,10) == 2) then
C              allocate(wk(lidimx,lidimx+1))
C              call yqinv('N',gf,ogi,ldrgn,2,lidimx,wk,lidimx,ierr)
C              call subdiag(kcplx,ldg,lidim,lhdim,0,isp,s_pot%gma,gf)
C              call yqinv('N',gf,ogi,ldrgn,2,lidimx,wk,lidimx,ierr)
C              deallocate(wk)
C              modg2g = modg2g - 200
C            endif
c           call sblhm2(lbloch,nl**2,1,nbas,s_ham%iprmb,s_pot%cp,isp,
c    .        lidim,ldg,ldg2,kl,gf)
            call ztoyy(gf,lidim,lidim,lidim,lidim,kcplx,1)
            call daxpy(2*lidim**2,1d0,s_pot%cp(1+lidim*lidim*(isp-1)),1,gf,1)
            call ztoyy(gf,lidim,lidim,lidim,lidim,1,kcplx)
C           call yprm('P-S',2,gf,ogi,ldg,lidim,lidim)
          endif
          if (iprint()>=110) call yprm('P-S',2,gf,ogi,ldg,lidim,lidim)

C     ... Correction to (P-S)_ll from i-waves
          if (.not. ldiwav) call gf1kpl(kcplx,nbas,s_ham%iprmb,ldim,idim,1d0,ldrgn,ogi,gf)
          if (iprint() >= 110 .and. idim > 0) call yprm('P-S + iwaves',2,gf,ogi,ldg,lidim,lidim)
        endif

C   ... g_ll = (P-S)_ll^-1
        if (ldiwav) ldimx = lidimx
C       call yprm('P-S',2,gf,ogi,ldg*nspc,ldimx,ldimx)
        call tcn('inversion')
        allocate(wk(ldimx,ldimx+1))
C       check that call is compatible with kcplx=0 and 2
        call yqinv(cyinv,gf,ogi,ldrgn,2,ldimx,wk,ldimx,ierr)
        call tcx('inversion')
        if (ierr /= 0) call rxi('gf1kp: g is singular, ierr=',ierr)
C       call yprmi('inv(P-S) isp=%i q=%s,(%3;6d)',isp,qp,2,gf,ogi,ldg*nspc,ldimx,ldimx)
        deallocate(wk)

C       Re-entry point : g_ll has already been calculated
  100   continue

C   --- g for i-waves ---
        if (idim > 0 .and. .not. ldiwav) then
          allocate(wk(2*ldim,idim))
C         call yprm('g before iwaves complete',2,gf,ogi,ldg*nspc,lidimx,lidimx)
          call gf1kpi(mod(mode,10),nbas,s_ham%iprmb,kcplx,ldim,idim,wk,wk,ldrg,nspc,ogi,gf)
C         call yprm('g after iwaves complete',2,gf,ogi,ldg*nspc,lidimx,lidimx)
          deallocate(wk)
        endif
C       call yprm('g_l+i,l+i',2,gf,ogi,ldg*nspc,lidimx,lidimx)

C   --- g for h-waves ---
        if (lhdim > lidim) then
          hdim = lhdim - lidim
C     ... Bloch sum s_hl
          allocate(shl(hdim*ldim*2))
          lblch2 = 50000 + mod(lbloch,10000)
          call bloch(lblch2,qp,nl,plat,nl**2,s_ham%iprmb,1,nptab,iax,
     .      salp,nl**2,1,1,lidim,ldim,hdim,0,hdim,ldim,kl,xx,shl,xx)
          call cplxdm(kcplx,hdim,ldim,ld11,ld21,ldrs,osi)
          allocate(wk(2*hdim,ldim))
C         Potential function offset starts at h-waves
          offp = offp + lidim
C     ... Make diagonal ghh
          call gf1kph(kcplx,nbas,nspc,s_ham%iprmb,ldim,lidim,hdim,
     .      pfa,offp,wk,ldrg,ogi,gf,ldrs,osi,shl,ldrh,ogh,ghh)
C         call yprm('ghh',2,ghh,ogh,ldrh,hdim,nspc*mxorb)
          deallocate(wk,shl)
        endif

        if (lncol /= 0 .and. lglobz) then
C     ... Rotate ll block to local z quantization axis (forward rotation)
          i = 0000 + 100*kcplx + 10*lsph + rotspd(1) ! Forward rotation
C         if (.true. .or. cmdopt('--spinoronly',12,0,strn)) i = i + 2000 ! Skip orbital part
          if (cmdopt('--spinoronly',12,0,strn)) i = i + 2000 ! Skip orbital part
          call rotheu(i,s_ham,nl,1,nbas,nbas,0,lidim,ldg,ldg,gf)
C         call yprm('gnc (local z)',2,gf,ogi,ldg*nspc,lidimx,lidimx)
        endif
        deallocate(pfa)
      endif ! 2nd gen ASA
      if (lncol /= 0) deallocate(su)
C ... End of block generating unscaled gf

C --- Add off-diagonal U contribution to LDA+U potential ---
      if (nlibu > 0 .and. ludiag == 0) then
        if (icpa == 1) call rx('update gf1kp for LDA+U, CPA case')
        call gfg2g(s_ham,s_pot,200,lrel,isp,1,nbas,1,nbas,lidim,lidim,gf,ghh,ierr)
        call tcn('inversion')
        allocate(wk(ldimx,ldimx+1))
        call yqinv(cyinv,gf,ogi,ldrgn,2,ldimx,wk,ldimx,ierr)
        call tcx('inversion')
C       call yprm('g^-1 before',2,gf,lidimx**2,lidimx,lidimx,lidimx)
        call asaddu(5,nbas,lidim,lmaxu,isp,nsp,nspc,nlibu,idu,
     .    s_ham%iprmb,vorb,gf,xx)
C       if (isp == 2) then;   print *, isp
C       call yprm('g^-1 after',2,gf,lidimx**2,lidimx,lidimx,lidimx)
C       endif
        call tcn('inversion')
        call yqinv(cyinv,gf,ogi,ldrgn,2,ldimx,wk,ldimx,ierr)
        deallocate(wk)
        call tcx('inversion')
        if (modg2g >= 200) modg2g = modg2g - 200
      endif

C     debugging checks for noncollinear cases.  Both gfg2g and gfg2gnc should handle these cases
C     fept tests downfolding (set up with gf/test/test.pgf 7)
C     run fept -vnk=8 -vnc=1 -vntht=6 -vnit=1 --pr31,30 --quit=rho -vidxp=2 -vgamma=0 --quit=rho
C     mnpt tests SO and fully relativistic cases (set up with cp mnpt-files/*.mnpt .)
C     run mnpt -vgamrep=1 -vrel=2 -vgfmod=26 -vso=1
C     call yprm('g before scaling',2,gf,ldgx*ldgx,ldgx,ldgx,ldgx)
C     print *, '!! gf1kp scale'; modg2g = modg2g + 100

C --- Convert g to g^gam by and/or to G by energy scaling, non-CPA case ---
C     if (modg2g/100 /= 0 .and. icpa /= 1 .and. .false.) then
      if (modg2g/100 /= 0 .and. icpa /= 1) then
        call gfg2g(s_ham,s_pot,modg2g,lrel,isp,1,nbas,1,nbas,lidim,lidim,gf,ghh,ierr)
C       call zprm('(P-S)^-1, gamma from gfg2g',2,gf,ldgx,ldgx,ldgx)

C --- Convert g to gamma, and/or to G by scaling, using general noncollinear gfg2gnc ---
C     This branch should correctly handle all magnetic cases implemented so far that with hdim=0
C     But we use gfg2g where applicable because it is more efficient.
      elseif ((nsp == 2 .or. icpa /= 0) .and. modg2g/100 /= 0 .and. hdim == 0) then
        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)

C       g^alp -> g^gam
        if (modg2g/100 >= 2) then
C         if (icpa == 1) call rx('gf1kp:  g^alp -> g^gam by scaling not available with CPA')
          call sugfg2gnc(lso,lrel,1,0,mod(icpa,2),.true.,i)
          call gfg2gnc(i,1,s_pot,s_site,isp,nsp,nspc,s_ham%offH,s_ham%iprmb,lidim,lidim,0,(/0,0/),gf,nbas)
C         call zprm('(P-S)^-1, gamma',2,gf,ldgx,ldgx,ldgx)
        endif
C       g -> G
        if (mod(modg2g/100,2) /= 0) then
C         if (icpa == 1) call rx('gf1kp:  gfg2gnc not equipped to scale g^alp -> G in CPA case')
          i = 1; if (lrel == 2 .or. lso) i = 2; if (icpa == 2) i = 5; if (icpa == 2 .and. lrel == 2) i = 6
          if (lrel == 2) i = 10000 + i
          call rx('gf1kp fix')
          call sugfg2gnc(lso,lrel,1,1,mod(icpa,2),.true.,i)
          call gfg2gnc(i,1,s_pot,s_site,isp,nsp,nspc,s_ham%offH,s_ham%iprmb,lidim,lidim,0,(/0,0/),gf,nbas)
C         call zprm('G',2,gf,ldgx,ldgx,ldgx)
        endif

        call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)

C --- Convert g to g^gam by and/or to G by energy scaling, Standard collinear case ---
      elseif (modg2g/100 /= 0) then
        call gfg2g(s_ham,s_pot,modg2g,lrel,isp,1,nbas,1,nbas,lidim,lidim,gf,ghh,ierr)
      endif

C --- Add sigma to g in gamma-representation ---
      if ((modg2g /= 0) .and. lrsig /= 0) then
        if (nspc == 2) call rx('gf1kp not ready for SX, nspc = 2')
        allocate(wk(ldimx,ldimx+1))
        call tcn('inversion')
        call yqinv(cyinv,gf,ogi,ldrgn,2,ldimx,wk,ldimx,ierr)
        call tcx('inversion')
        if (ierr /= 0) call rxi('gf1kp:sig1 is singular, ierr=',ierr)
c        call yprm('sigm-ad',2,sigma,lidim*lidim,lidim,lidim,lidim)
        call sigmadd(lidim,sigma,gf)
        call tcn('inversion')
        call yqinv(cyinv,gf,ogi,ldrgn,2,ldimx,wk,ldimx,ierr)
        call tcx('inversion')
        if (ierr /= 0) call rxi('gf1kp:sig2 is singular, ierr=',ierr)
        deallocate(wk,sigma)
      endif

      if (iprint() >= PRTG/1) then
        call yprm('g',2,gf,ogi,ldg*nspc,lidimx,lidimx)
        if (hdim/=0) call yprm('ghh',2,ghh,ogh,ldrhn,hdim*nspc,mxorb)
      endif
      call tcx('gf1kp')

C     call poppr
      end

      subroutine gf1kpn(mode,kcplx,nbas,indxsh,ldim,lidim,ldg,ldrg,osi,ogi,pfnc,wk,gf)
C- Noncollinear GF when off-diagonal spinor parts are diagonal in basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode: :0 require only diagonal parts of g
Ci          1 require only only site diagonal g
Ci         >1 require entire g
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
Ci          1: complex*16: gf = gf(ldg,ldg)
Ci          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
Ci   ldim  :dimension of hamiltonian matrix, lower block (makidx.f)
Ci   lidim :dimension of hamiltonian matrix (makidx.f)
Ci   ldg   :formal leading dimension of gf,wk; see kcplx
Ci   ldrg  :true leading dimension of gf,wk; see kcplx
Ci   ldrg  :leading dimension of real part of s; generate with
Ci          call cplxdm(kcplx,ldg,ldg,ld11,ld21,ldrg,osi)
Ci   osi   :offset to imaginary part of s; generate with
Ci          call cplxdm(kcplx,ldg,ldg,ld11,ld21,ldrg,osi)
Ci   ogi   :offset to imaginary part of g
Ci   pfnc  :spinor diagonal in basis, nondiagonal in spin space
Ci         :It is assumed that the spin-up and spin-down potential
Ci         :functions for the iwaves are the same
Ci   wk    :work array of dimension ldg*ldg*2
Cio Inputs/Outputs
Cio  gf    :On input, g = h, zero in off-diagonal spin blocks
Cio         On output, g = (pfnc + h)^-1
Cl Local variables
Cl   ldgn  :leading dimension of gf, treated as a 2d array
Cl          with spin and orbital indices combined
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
      integer nbas,indxsh(*)
      integer mode,kcplx,ldim,lidim,ldg,ldrg,osi,ogi
      double precision gf(ldrg,2,ldg,2),wk(ldrg,ldg)
      double precision pfnc(2,ldg,2,2)
C Local variables
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer ldgn,ierr,i,j,ib,i1,idim
      integer nlm,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb,offi
      double complex xx

      call tcn('gf1kpn')

      ldgn = 2*ldrg
      idim = lidim - ldim

C --- Make contribution to ll block of gf from intermediate waves ---
C     Do only once, since contribution is equal for both spin channels
C     Make correction in 21 block: copy 11 to 21 block
      if (idim /= 0) then
      call ymscop(0,lidim,lidim,ldrg*2,ldrg*2,0,0,0,0,gf,ogi,gf(1,2,1,1),ogi)
      do  1  i = ldim+1, lidim
        gf(i,2,i,1)= gf(i,2,i,1) + pfnc(1,i,1,1)
        gf(i+ogi,2,i,1)= gf(i+ogi,2,i,1) + pfnc(2,i,1,1)
    1 continue
C ... i-waves correction to gll.  Store in g(sig=2,sig'=1)
      call gf1kpl(kcplx,nbas,indxsh,ldim,idim,0d0,ldgn,ogi,gf(1,2,1,1))
C ... Add i-waves correction to gll in 11 and 22 blocks
      do  i = 1, 2
        call ymscop(1,ldim,ldim,ldrg*2,ldrg*2,0,0,0,0,gf(1,2,1,1),ogi,gf(1,i,1,i),ogi)
      enddo
C ... Copy g_li and gl+i,i from 21 block into 11 block; needed by gf1kpi
      call ymscop(0,lidim,idim,ldrg*2,ldrg*2,0,0,0,0,gf(1,2,1+ldim,1),ogi,gf(1,1,1+ldim,1),ogi)
      endif

C --- Add diagonal to g ---
C     call zprm('pfnc',2,pfnc,ldg,ldg,4)
      call daxpy(ldim,1d0,pfnc,2,gf,ldgn+1)
      call daxpy(ldim,1d0,pfnc(2,1,1,1),2,gf(1+ogi,1,1,1),ldgn+1)
      call daxpy(ldim,1d0,pfnc(1,1,2,2),2,gf(1,2,1,2),ldgn+1)
      call daxpy(ldim,1d0,pfnc(2,1,2,2),2,gf(1+ogi,2,1,2),ldgn+1)
C     call yprm('P-S',kcplx+2,gf,ogi,ldg*2,lidim*2,lidim*2)

C ... a11^-1
      call yyqinv('N',gf,gf(1+ogi,1,1,1),ldgn,2,ldim,wk,ldim,ierr)
      if (ierr /= 0) call rx('GF1KPN: matrix singular')
C     call yprm('s(1,1)^-1',kcplx+2,gf,ogi,ldg*2,ldim,ldim)

C --- g22 = (a22 - a21 a11^-1 a12)^-1 ---
      do  10  j = 1, ldim
      do  10  i = 1, ldim
        xx = dcmplx(pfnc(1,i,2,1),pfnc(2,i,2,1)) *
     .       dcmplx(gf(i,1,j,1),gf(i+ogi,1,j,1)) *
     .       dcmplx(pfnc(1,j,1,2),pfnc(2,j,1,2))
        gf(i,2,j,2) = gf(i,2,j,2) - dble(xx)
        gf(i+ogi,2,j,2) = gf(i+ogi,2,j,2) - dimag(xx)
   10 continue
      call yqinv('N',gf(1,2,1,2),ogi,ldgn,2,ldim,wk,ldim,ierr)
      if (ierr /= 0) call rx('GF1KPN: matrix singular')
C     call yprm('g22',kcplx+2,gf(1,2,1,2),ogi,ldg*2,ldim,ldim)

C --- g12 = -a11^-1 a12 g22 ---
      do  20  j = 1, ldim
      do  20  i = 1, ldim
        xx = dcmplx(gf(i,1,j,1),gf(i+ogi,1,j,1)) *
     .       dcmplx(pfnc(1,j,1,2),pfnc(2,j,1,2))
        wk(i,j)     = - dble(xx)
        wk(i+osi,j) = - dimag(xx)
   20 continue
      call ygemm('N','N',ldim,ldim,ldim,1d0,wk,osi,ldg,gf(1,2,1,2),ogi,
     .  ldg*2,0d0,gf(1,1,1,2),ogi,ldg*2)
C     call yprm('g12',kcplx+2,gf(1,1,1,2),ogi,ldg*2,ldim,ldim)

C ... wk = a21 a11^-1
      do  30  j = 1, ldim
      do  30  i = 1, ldim
        xx = dcmplx(pfnc(1,i,2,1),pfnc(2,i,2,1)) *
     .       dcmplx(gf(i,1,j,1),gf(i+ogi,1,j,1))
        wk(i,j)     = - dble(xx)
        wk(i+osi,j) = - dimag(xx)
   30 continue

C --- Make all or part of g11, g21, depending on mode --
      if (mode == 0 .or. mode == 1) then
        do  40  ib = 1, nbas

        call orbl(ib,0,ldim,indxsh,norb,ltab,ktab,offi,offl,nlm)
        if (nlm == 0) goto 40
        i1 = offl(1)+1

C   ... g11 = a11^-1 - g12 a21 a11^-1 for site ib
        call ygemm('N','N',nlm,nlm,ldim,1d0,gf(i1,1,1,2),ogi,
     .    ldg*2,wk(1,i1),osi,ldg,1d0,gf(i1,1,i1,1),ogi,ldg*2)

C   ... g21 = - g22 a21 a11^-1
        call ygemm('N','N',nlm,nlm,ldim,1d0,gf(i1,2,1,2),ogi,
     .    ldg*2,wk(1,i1),osi,ldg,0d0,gf(i1,2,i1,1),ogi,ldg*2)

   40   continue

      else
C   ... g11 = a11^-1 - g12 a21 a11^-1
        call ygemm('N','N',ldim,ldim,ldim,1d0,gf(1,1,1,2),ogi,ldg*2,
     .    wk,osi,ldg,1d0,gf(1,1,1,1),ogi,ldg*2)

C   ... g21 = - g22 a21 a11^-1
        call ygemm('N','N',ldim,ldim,ldim,1d0,gf(1,2,1,2),ogi,ldg*2,
     .    wk,osi,ldg,0d0,gf(1,2,1,1),ogi,ldg*2)
      endif

C      call yprm('g11',kcplx+2,gf(1,1,1,1),ogi,ldg*2,ldim,ldim)
C      call yprm('g21',kcplx+2,gf(1,2,1,1),ogi,ldg*2,ldim,ldim)

C     call yprm('g before rot',kcplx+2,gf,ogi,ldg*2,lidim*2,lidim*2)
      call tcx('gf1kpn')

      end

      subroutine gf1kpl(kcplx,nbas,indxsh,ldim,idim,fac,ldrg,ogi,gf)
C- Make contribution to gf_ll from intermediate waves, collinear case
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: gf = gf(ldrg,ldrg,1..2)
Ci          2: real, imaginary in columns : gf = gf(ldrg/2,1..2,ldrg/2)
Ci   nbas  :size of basis
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   fac   :replace ll block with iwaves correction + fac * ll block
Ci   ldrg  :leading dimension of gf
Ci   ogi   :offset to imaginary part of gf
Cio Inputs/Outputs
Cio                        ( P-Sll     -Sli)
Cio  gf    :On input, gf = (               )
Cio                        (  -Sil    P-Sii)
Cio
Cio         with Sil = Sli+
Cio
Cio                        ( P-Sll - dP   Sli diag(P-Sii)^-1 )
Cio        :On output, g = (                                 )
Cio                        (  -Sil        diag(P-Sii)^-1     )
Cio
Cio         where  diag(P-Sii)^-1 is an approximate inverse of (P-Sii)
Cio         and    dP= Sli diag(P-Sii)^-1 Sil
Cr Remarks
Cr                  ( P-Sll     -Sli)
Cr   Let  g   = inv (               )
Cr                  (  -Sil    P-Sii)
Cr
Cr   then gll = (P-Sll - (-Sli) (P-Sii)^-1 (-Sil))^-1 .
Cr   This routine subtracts  (-Sli) (P-Sii)^-1 (-Sil) ,
Cr   approximating P-Sii^-1 with a block diagonal inverse.
Cr   It is well justified for higher waves, i.e. when P>>S
Cu Updates
Cu   8 Mar 00 use block-diagonal g
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,nbas,indxsh(*),ldim,idim,ldrg,ogi
      double precision gf(ldrg,ldrg),fac
C ... Local parameters
      integer nlmx,n0,nkap0
      parameter (nlmx=81,n0=10,nkap0=4)
      integer i,j,ibas,lidim,nlm,norb,offi,ofa,
     .  ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      double complex pi(nlmx),pij

      if (idim == 0) return
      lidim = ldim+idim

      call tcn('gf1kpl')
      if (kcplx /= 0 .and. kcplx /= 2) call rx('gf1kpl needs kcplx=0 or 2')

C      call yprm0('(1p,9e18.10)')
C      call yprm('gii',2,gf,ogi,lidim,lidim,lidim)

C --- Overwrite diagonal blocks of g_ii with approximate inverse: ---
C     g_ij approx [diag (P-S)_i] g_ij [diag (P-S)_j]^-1
      do  ibas = 1, nbas
      call orbl(ibas,ldim,lidim,indxsh,norb,ltab,ktab,offi,offl,nlm)
      if (nlm == 0) cycle
      if (nlm > nlmx) call rx('gf1kpl: increase nlmx')
      ofa = offl(1)             ! Offset to intermediate block
      do  i = 1, nlm
        pi(i) = 1/dcmplx(gf(i+ofa,i+ofa),gf(i+ofa+ogi,i+ofa))
      enddo
      do  i = 1, nlm
      do  j = 1, nlm
        pij = pi(i)*dcmplx(gf(i+ofa,j+ofa),gf(i+ofa+ogi,j+ofa))*pi(j)
        if (i /= j) pij = -pij
        gf(i+ofa,j+ofa) = dble(pij)
        gf(i+ofa+ogi,j+ofa) = dimag(pij)
      enddo
      enddo

C ... Zero out gii except elements on block diagonal
C      do  32  j = ofa+1, ofa+nlm
C        do  33  i = ldim+1, ofa
C          gf(i,j) = 0
C          gf(i+ogi,j) = 0
C   33   continue
C        do  34  i = ofa+nlm+1, lidim
C          gf(i,j) = 0
C          gf(i+ogi,j) = 0
C   34   continue
C   32 continue
      enddo
C     call yprm('gii',2,gf,ogi,lidim,lidim,lidim)

C --- Overwrite g_li with  - (-Sil+) [diag (P-S)_ii]^-1  ---
      do  ibas = 1, nbas
        call orbl(ibas,ldim,lidim,indxsh,norb,ltab,ktab,offi,offl,nlm)
        if (nlm == 0) cycle
        ofa = offl(1)
        call ygemm('C','N',ldim,nlm,nlm,-1d0,gf(1+ofa,1),ogi,ldrg,
     .    gf(1+ofa,1+ofa),ogi,ldrg,0d0,gf(1,1+ofa),ogi,ldrg)
      enddo

C --- gf_ll -=  (-Sli) [diag (P-S)_ii]^-1 (-Sil) ---
      call ygemm('N','N',ldim,ldim,idim,1d0,gf(1,ldim+1),ogi,
     .  ldrg,gf(1+ldim,1),ogi,ldrg,fac,gf,ogi,ldrg)

C     call yprm('gii',2,gf(1,1),ogi,ldrg,lidim,lidim)
      call tcx('gf1kpl')

      end

      subroutine gf1kpi(mode,nbas,indxsh,kcplx,ldim,idim,wk,wk2,ldrg,nspc,ogi,gf)
C- gf for intermediate waves, g_ii
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode: :0 require only diagonal parts of g
Ci          1 require only only site diagonal g
Ci         >1 require entire g
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: gf = gf(ldrg,ldrg,1..2)
Ci          2: real, imaginary in columns : gf = gf(ldrg/2,1..2,ldrg/2)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   wk    :double complex work array of dimension ldim*idim
Ci   wk2   :double complex work array of same dimension as wk
Ci         :and may point to the same address space as wk
Ci   ldrg  :leading dimension of gf
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ogi   :offset to imaginary part of gf
Cio Inputs/Outputs
Cio        :               (  gll   -Sli      )
Cio  gf    :On input, gf = (                  )
Cio        :               ( -Sil  (P-Sii)^-1 )
Cio        :
Cio        :In the noncollinear case, gf is a 2x2 supermatrix of the above:
Cio        :               ( g11   g12 )
Cio        :           gf= (           )
Cio        :               ( g21   g22 )
Cio        :
Cio        :  The ii block of each
Cio        :
Cio        :               ( gll   gli )
Cio        :On output, gf= (           )
Cio        :               ( gil   gii )
Cio        :
Cio        : gli = -gll Sli (P-Sii)^-1
Cio        : gil = -(P-Sii)^-1 Sil gll
Cio        : gii = (P-Sii)^-1 + (P-Sii)^-1 Sil gll Sli (P-Sii)^-1
Cio        : wk2 = (P-Sii)^-1 Sil gll
Cio
Cio        : Noncollinear case: Make 11 and 22 blocks only
Cio        : Since the spin off-diagonal blocks of Sil, Sii are zero,
Cio        : gii_11 = (P-Sii)^-1 + (P-Sii)^-1 Sil gll_11 Sli (P-Sii)^-1
Cio        : gii_22 = (P-Sii)^-1 + (P-Sii)^-1 Sil gll_22 Sli (P-Sii)^-1
Cr Remarks
Cr   This routine must be called using a global quantization axis.
Cr
Cr  *1 It is assumed that P_ii for intermediate waves is independent of
Cr   spin, so P_ii is spin diagonal.  S is intrinsically spin diagonal.
Cr  *2 It is assumed that inv(P-S)_ii is site diagonal, and only the
Cr   site-diagonal part of g_ii is calculated. In any case, inv(P-S)_ii
Cr   is calculated approximately; see gf1kpl.
Cr  *3 In the noncollinear case,
Cr   gf1kpi assumes sil for up-up and down-down spins are the same
Cu Updates
Cu   8 Mar 00 use block-diagonal g
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,indxsh(*),kcplx,ldim,idim,ldrg,nspc,ogi
      double precision gf(ldrg,nspc,ldrg,nspc),wk(ldim,2,idim),wk2(idim,2,ldim)
C ... Local parameters
      integer n0,nlmx,nkap0
      parameter (nlmx=81,n0=10,nkap0=4)
      integer i,j,ldgn,is,lidim,ibas,nlm,norb,offi,ofa,ofb,
     .  ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)

C      call xx(ldrg,0,gf)
C      call xx(ldrg,1,gf)

      if (idim == 0) return
      ldgn = ldrg*nspc
      lidim = ldim+idim
      call rxx(kcplx /= 0.and.kcplx /= 2,'gf1kpi needs kcplx=0 or 2')
      call tcn('gf1kpi')

C --- Overwrite gil with -(P-Sii)^-1 (-Sil) ---
C     In the noncollinear case, this quantity is written to the 11 spin block.
C     It is the same for the 11- and 22- blocks and zero for the spin off-diagonal blocks.
      call ymscop(0,idim,ldim,ldgn,idim*2,0,0,0,0,gf(1+ldim,1,1,1),ogi,wk2,idim) ! wk2 <- -Sil
      do  ibas = 1, nbas  ! See Remarks *1
        call orbl(ibas,ldim,lidim,indxsh,norb,ltab,ktab,offi,offl,nlm)
        if (nlm == 0) cycle
        if (nlm > nlmx) call rx('gf1kpi: increase nlmx')
        ofa = offl(1)           ! Offset to intermediate block
        ofb = offl(1)-ldim      ! Is normally if i block follows l block
        call ygemm('N','N',nlm,ldim,nlm,-1d0,gf(1+ofa,1,1+ofa,1),ogi,
     .    ldgn,wk2(1+ofb,1,1),idim,idim*2,0d0,gf(1+ofa,1,1,1),ogi,ldgn)
      enddo
C ... Noncollinear case: copy inv(P-Sii)_11 to the 22 block
      if (nspc == 2) call ymscop(0,idim,idim,ldgn,ldgn,0,0,0,0,
     .  gf(1+ldim,1,1+ldim,1),ogi,gf(1+ldim,2,1+ldim,2),ogi)

C --- Make gii (noncollinear case: gii_11 and gii_22) ---
C ... wk <- (Sli (P-Sii)^-1)
      call ymscop(0,ldim,idim,ldgn,ldim*2,0,0,0,0,gf(1,1,1+ldim,1),ogi,wk,ldim)
      do  is = 1, nspc

C   ... gli(is,is)  <-  gll(is,is) (Sli (P-Sii)^-1)
        call ygemm('N','N',ldim,idim,ldim,1d0,gf(1,is,1,is),ogi,ldgn,
     .    wk,ldim,ldim*2,0d0,gf(1,is,1+ldim,is),ogi,ldgn)

C   ... Add to diag. gii :  -(P-Sii)^-1 (-Sil) gll Sli (P-Sii)^-1
        do  ibas = 1, nbas
          call orbl(ibas,ldim,lidim,indxsh,norb,ltab,ktab,offi,offl,nlm)
          if (nlm == 0) cycle   !
          if (nlm > nlmx) call rx('gf1kpi: increase nlmx')
          ofa = offl(1)           ! Offset to intermediate block
          call ygemm('N','N',nlm,nlm,ldim,1d0,gf(1+ofa,1,1,1),ogi,ldgn,
     .      gf(1,is,1+ofa,is),ogi,ldgn,1d0,gf(1+ofa,is,1+ofa,is),ogi,ldgn)
C     ... Zero out gii except elements on block diagonal
          do  j = ofa+1, ofa+nlm
            do  i = ldim+1, ofa
              gf(i,is,j,is) = 0
              gf(i+ogi,is,j,is) = 0
            enddo
            do  i = ofa+nlm+1, lidim
              gf(i,is,j,is) = 0
              gf(i+ogi,is,j,is) = 0
            enddo
C       ... Zero out spin-off-diagonal parts of gii
            if (is == 2) then
              do  i = ldim+1, lidim
                gf(i,1,j,2) = 0
                gf(i+ogi,1,j,2) = 0
                gf(i,2,j,1) = 0
                gf(i+ogi,2,j,1) = 0
              enddo
            endif
          enddo
        enddo
        enddo

C --- gil <- - (P-Sii)^-1 Sil gll ---
      if (mode > 1) then
        call ymscop(0,idim,ldim,ldgn,idim*2,0,0,0,0,gf(1+ldim,1,1,1),ogi,wk2,idim)
        do  is = 1, nspc
          call ygemm('N','N',idim,ldim,ldim,1d0,wk2,idim,idim*2,
     .      gf(1,is,1,is),ogi,ldgn,0d0,gf(1+ldim,is,1,is),ogi,ldgn)
        enddo
      endif

      call tcx('gf1kpi')

      end

      subroutine gf1kph(kcplx,nbas,nspc,indxsh,ldim,lidim,hdim,pfun,
     .  offp,wk,ldrg,ogi,gf,ldrs,osi,shl,ldrh,ogh,ghh)
C- gf for higher waves
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: gf = gf(ldrg,ldrg,1..2)
Ci          2: real, imaginary in columns : gf = gf(ldrg/2,1..2,ldrg/2)
Ci   nbas  :size of basis
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lidim :number of lower+intermediate orbitals
Ci   hdim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   pfun  :vector of potential functions
Ci         :It is assumed that the spin-up and spin-down potential
Ci         :functions for the hwaves are the same in noncollinear case
Ci   offp  :offset to 1st entry in pfun for this higher block
Ci   ldrg  :leading dimension of gf
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ogi   :offset to imaginary part of gf
Ci   gf    :gll (lower block of g)
Ci   ldrs  :leading dimension of shl
Ci   osi   :offset to imaginary part of shl
Ci   shl   :strux connecting higher and lower blocks
Ci   ldrh  :leading dimension of ghh
Ci   ogh   :offset to imaginary part of ghh
Co Outputs
Co   ghh   :higher block of diagonal part of g
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,nbas,indxsh(*),ldim,lidim,hdim,ldrg,nspc,ogi,offp,ldrs,osi,ldrh,ogh
      double precision gf(ldrg,nspc,ldrg,nspc),shl(ldrs,1),ghh(ldrh,nspc,1),wk(hdim,2,ldim)
      double complex pfun(*)
C ... Local parameters
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer i,j,is,lhdim,nlm,ibas,iorb,l,ofa,ofb,nlmi,ldhn
      integer ltab(n0*nkap0),ktab(n0*nkap0),norb,offi,offl(n0*nkap0)
      double complex psgsp,pfi,pfj

C     if (iprint() >= 110)
C     call yprm('shl',2,shl,osi,ldrs,hdim,ldim)
C     if (iprint() >= 110)
C     call yprm('gll',2,gf,ogi,ldrg,ldim,ldim)
C     if (iprint() >= 110)
C     call zprm('pfun',2,pfun,hdim,lidim+hdim,nspc**2)

      if (hdim == 0) return
      lhdim = lidim + hdim
      ldhn = nspc*ldrh

      call rxx(kcplx /= 0.and.kcplx /= 2,'gf1kph needs kcplx=0 or 2')
      call tcn('gf1kph')

      do  10  is = 1, nspc

C   --- wk = -shl gll(is,is) ---
        call ygemm('N','N',hdim,ldim,ldim,1d0,shl,osi,ldrs,
     .    gf(1,is,1,is),ogi,ldrg*nspc,0d0,wk,hdim,hdim*2)
C       call yprm('-shl gll',4,wk,ogi,hdim,hdim,ldim)

C   --- Add to diag. gii :  P^-1 Shl gll (Shl+) P^-1 ---
        do  20  ibas = 1, nbas
        call orbl(ibas,lidim,lhdim,indxsh,norb,ltab,ktab,offi,offl,nlm)
        if (nlm == 0) goto 20
        do  22  iorb = 1, norb
          l = ltab(iorb)
          ofa = offl(iorb)-lidim
          ofb = offl(iorb)-offl(1)
          nlmi = 2*l + 1
          call ygemm('N','C',nlmi,nlmi,ldim,1d0,wk(1+ofa,1,1),hdim,hdim*
     .      2,shl(1+ofa,1),osi,ldrs,0d0,ghh(1+ofa,is,1+ofb),ogh,ldhn)
          do  24  i = 1, nlmi
          do  24  j = 1, nlmi
            psgsp = dcmplx(ghh(i+ofa,is,j+ofb),ghh(i+ofa+ogh,is,j+ofb))
            pfi = pfun(i+ofa+offp)
            pfj = pfun(j+ofa+offp)
            psgsp = psgsp/pfi/pfj
            if (i == j) psgsp = psgsp + 1/pfi
            ghh(i+ofa,is,j+ofb) = dble(psgsp)
            ghh(i+ofa+ogh,is,j+ofb) = dimag(psgsp)
   24     continue
C       call yprm('ghh',kcplx+2,ghh(1+ofa,is,1),ogh,hdim,nlmi,nlmi)
   22   continue
   20   continue
   10 continue
      call tcx('gf1kph')

C     mxorb = nglob('mxorb')
C     call yprm('ghh',2,ghh,ogh,ldhn,hdim*nspc,mxorb)

      end
      subroutine gf1kps(nsp,ldim,lhdim,lbar,pf,pfa)
C- Spin-average potential functions
      implicit none
      integer nsp,ldim,lhdim,lbar
      double complex pf(lhdim,nsp),pfa(lhdim,nsp)
      integer i

      call dcopy(nsp*2*lhdim,pf,1,pfa,1)
      if (nsp == 2) then
        do  10  i = ldim+1, lhdim
          if (lbar == 0) then
            pfa(i,1)   = pf(i,1)
            pfa(i,nsp) = pf(i,nsp)
          else
            pfa(i,1)   = (pf(i,1)+pf(i,nsp))/2
            pfa(i,nsp) = (pf(i,1)+pf(i,nsp))/2
          endif
   10   continue
      endif

      end
