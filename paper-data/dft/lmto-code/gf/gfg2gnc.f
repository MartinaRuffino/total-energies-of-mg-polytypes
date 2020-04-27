      subroutine gfg2gnc(mode,nk,s_pot,s_site,isp,nsp,nspc,offH,iprmb,ldi,ldj,ib0,ijcomp,gij,nbas)
C- Convert g_ij to G_ij by energy scaling, general noncollinear case
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dpfr ddpfr
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit controls what potential functions are used.
Ci         :0-4 read potential functions in vector form, 5-6 in matrix form.  Key+ is below.
Ci             rank  rep  comp  make  scale  g               Add to diagonal
Ci         :0   v2    lms  CPA  g->G  sqrt(s_site(ib)%dpfr)  s_site(ib)%ddpfr
Ci         :1   v2    lms  1    g->G  sqrt(s_pot%dpf)        s_pot%ddpf
Ci         :2   v4    lms  1    g->G  s_pot%dpfr             s_pot%ddpfr
Ci         :3   v2    lms  1    ggam  s_pot%(palp/pgam)      s_pot%(palp/pgam)*s_pot%gma
Ci         :4   v4    lms  1    ggam  s_pot%papg             s_pot%gmar
Ci         :5   m4    lms  CPA  g->G  s_site(ib)%dpfr        s_site(ib)%ddpfr
Ci         :6   m4    kmu  CPA  g->G  s_site(ib)%dpfr        s_site(ib)%ddpfr
Ci         +Key
Ci           v2   potential functions are spin diagonal and have nsp spin indices, stored as vectors
Ci           v4   potential functions have 4 spin components, stored as vectors
Ci           m4   potential functions have 4 spin components, stored matrix form (1:norb,1:2,1:norb,1:2)
Ci           lms  potential functions are in lms representation
Ci           kmu  potential functions are in kappa-mu representation
Ci           CPA  means all components for CPA treatment are available (for CPA)
Ci           g->G used to scale g to G
Ci           ggam used to scale g to g^gam
Ci         :10s digit
Ci         :0 Scale lower block of g only
Ci         :1 Scale lower and intermediate blocks of g
Ci         :100s digit
Ci         :0 forward sense: make g->G
Ci         :1 inverse sense: make G->g
Ci         :1000s digit
Ci         :0 work with a full matrix G (g)
Ci         :1 gij(ldi=norb,2,ldj,2) rows 1:norb belong to site ib0
Ci         :2 gij(ldi,2,ldj=norb,2) cols 1:norb belong to site ib0
Ci         :10000s digit (applies to 1s digit modes 2-4)
Ci         :0 SR convention for m ordering (used with 1s digit mode=1)
Ci         :1 FR convention for m ordering (used with 1s digit mode=1)
Ci         :100000s digit kcplx (see ztoyy.f)
Ci         :0 Input gij has real, imaginary separated: gij = gij(ldi,nspc,ldj,nspc,1:2)
Ci         :1 Input gij in complex*16 format.  Internally gfg2gnc converts to this format
Ci         :2 Input gij has real, imaginary separated by columns: gij = gij(ldi,nspc,1:2,ldj,nspc)
Ci         :On output, gij is converted to the same format is input
Ci   nbas  :number of sites
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ldi   :leading dimension of gij
Ci   ldj   :second dimension of gij
Ci   ib0   :site index, to be used when 1000s digit mode is 1 or 2
Ci   gij   :g_ij if 100s digit eq 0, G_ij if 1
Co Outputs
Co   gij   :gij is scaled from g to G (G >> g if 100s digit eq 1)
Cr Remarks
Cb Bugs
Cb   Downfolding not implemented when P given in matrix form (modes 2 and 3)
Cb   For it to work, pgfg2gnc needs to permute rows and columns of P
Cu Updates
Cu   25 Jun 16 (MvS) Bundle loop over kp to speed energy scaling
Cu   04 Apr 16 (MvS) Extended to collinear, spin polarized case
Cu   25 Mar 16 (MvS) Unify different sources of P including collinear P. New ability to scale g->g^gam
Cu    3 Mar 16 (Vishina) added G->g and ability to scale one site only
Cu   21 Jan 16 (MvS) first created
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nk,nbas,ldi,ldj,ib0,ijcomp(2),isp,nsp,nspc
      integer, parameter :: nkap0=4, n0H=5
      integer offH(n0H,nkap0,nbas+1),iprmb(*)
      complex(8), target :: gij(nk,ldi,nspc,ldj,nspc)
C ... For structures
!      include 'structures.h'
      type(str_pot) ::  s_pot
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8),allocatable, dimension(:,:,:,:) :: dpfr, wk(:,:,:,:,:)
      complex(8),pointer :: pfv(:,:),pfm(:,:),gijl(:,:,:,:,:)
C ... Local parameters
      integer mod0,mod1,mod2,mod3,mod4,kcplx
      integer ik,ib,idx,norb,n2,mxorb,is,js,offg,offL,offR,norbi,offp,isw,idnf,lidim,lmr,icomp
      procedure(integer) :: nglob

C     print *, 'enter gfg2gnc',mode; call rx('test')
      if (nk <= 0) return
      call tcn('gfg2gnc')

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      mod3 = mod(mode/1000,10)
      mod4 = mod(mode/10000,10)
      kcplx = mod(mode/100000,10)

      call ztoyy(gij,nk*ldi*nspc,nk*ldi*nspc,ldj*nspc,ldj*nspc,kcplx,1)
      if (nk > 1) then
        allocate(gijl(1,ldi,nspc,ldj,nspc))
      else
        gijl => gij
      endif

      call sanrg(.true.,mod1,0,1,'gfg2gnc:','mod1')
      if (.not. (nspc == 2 .or. mod0 == 0 .or. mod0 == 1 .or. mod0 == 3))
     .  call rx('gfg2gnc : inconsistent mode')
      mxorb = nglob('mxorb')     ! dimensions dpfr and sufficient space for pfv (my
      allocate(wk(ldi,nspc,ldj,nspc,nk))           ! Intermediate sqrt(pdot) g
      allocate(dpfr(mxorb,nspc,mxorb,2))            ! pf in matrix form, for all forms
      allocate(pfv(mxorb,4))                        ! large enough to hold vector pf at any site
      idnf = min(offH(2,1,nbas+1)-offH(2,1,1),1)    ! 0 if no downfolding, 1 if there is downfolding
      lidim = offH(4,1,nbas+1)                      ! number of lower+intermediate orbitals in iprmb array
      if (mod0 >= 5 .and. idnf /= 0) call rx('downfolding in gfg2gnc not implemented for CPA style pf')
      dpfr = 0
      wk = 0

C --- Make G <- G - pdotdot/pdot ---
      if (mod2 == 1) then      ! for G -> g
      do  ib = 1, nbas
        if (mod3 /= 0 .and. ib /= ib0) cycle
        icomp = 1
        if (s_site(ib)%ncomp > 1 .and. (mod0 == 0 .or. mod0 >= 5)) icomp = ijcomp(1)
        if (icomp < 1 .or.  icomp > s_site(ib)%ncomp) cycle
        norb = s_site(ib)%norb ; n2 = nspc*norb
        lmr = mxorb*(ib-1)
        pfm => s_site(ib)%ddpfr
        call pgfg2gnc2(100*idnf+10+mod0,s_pot,s_site,ib,icomp,norb,iprmb(1+lmr),pfv) ! Set up pfv
        isw = 10000*mod4 + 100*idnf + mod0
C       ? Why is 1:norb*nspc*norb*nsp necessary?
        call pgfg2gnc(isw,isp,nsp,nspc,norb,mxorb,pfv,pfm(1:norb*nspc*norb*nsp,icomp),dpfr)

C   ... Subtract pdotdot/pdot from g
        offg = offH(1,1,ib)     ! Offset to lower block of g for this site
        offp = 0                ! Offset to lower block for local dpfr
        do  idx = 1, 1+mod1
          norbi = offH(idx,1,ib+1)-offH(idx,1,ib)  ! Size of downfolding block
          if (norbi == 0) cycle
          offL = offg; offR = offg; if (mod3 == 1) offL = 0; if (mod3 == 2) offR = 0
          do  ik = 1, nk
          do  is = 1, nspc
          do  js = 1, nspc
            gij(ik,offL+1:offL+norbi,is,offR+1:offR+norbi,js) =
     .      gij(ik,offL+1:offL+norbi,is,offR+1:offR+norbi,js) - dpfr(offp+1:norbi,is,offp+1:norbi,js)
          enddo
          enddo
          enddo
          offg = offH(2,1,ib) + offH(1,1,nbas+1)  ! Offset to intermediate block of g for this site
          offp = offp + norbi   ! Offset to intermediate block for local dpfr
        enddo                   ! downfolding blocks
      enddo                     ! sites
      endif

C --- Make wk <- sqrt(pdot) g ---
      do  ib = 1, nbas
C       if (s_site(ib)%ncomp > 1) cycle ! CPA case => g corresponds to g for first component
        if (mod3 == 1 .and. ib /= ib0) cycle
        icomp = 1
        if (s_site(ib)%ncomp > 1 .and. (mod0 == 0 .or. mod0 >= 5)) icomp = ijcomp(1)
        if (icomp < 1 .or.  icomp > s_site(ib)%ncomp) cycle
        norb = s_site(ib)%norb ; n2 = 2*norb
        lmr = mxorb*(ib-1)
        pfm => s_site(ib)%dpfr
        call pgfg2gnc2(100*idnf+0+mod0,s_pot,s_site,ib,icomp,norb,iprmb(1+lmr),pfv) ! Set up pfv
        isw = 10000*mod4 + 1000*mod2 + 100*idnf + mod0
C       ? Why is 1:norb*nspc*norb*nsp necessary?
        call pgfg2gnc(isw,isp,nsp,nspc,norb,mxorb,pfv,pfm(1:norb*nspc*norb*nsp,icomp),dpfr)
C       call zprm('sqrt(pdot)',2,dpfr,mxorb*nspc,norb*nspc,norb*nspc)  ! printout only works if mxorb = norb

C       It would be more efficient to copy transpose of ib subblock block only
        do  ik = 1, nk
        if (nk > 1) then
          call zcopy(ldi*nspc*ldj*nspc,gij(ik,1,1,1,1),nk,gijl,1)
        endif

        offg = offH(1,1,ib)     ! Offset to lower block of g for this site
        offp = 0                ! Offset to lower block for local dpfr
        do  idx = 1, 1+mod1
          norbi = offH(idx,1,ib+1)-offH(idx,1,ib)  ! Size of downfolding block
          if (norbi == 0) cycle
          offL = offg; offR = offg; if (mod3 == 1) offL = 0; if (mod3 == 2) offR = 0
          do  is = 1, nspc
          do  js = 1, nspc
            call zgemm('T','N',norbi,ldj,norbi,(1d0,0d0),dpfr(1+offp,1,1+offp,is),mxorb*nspc,
     .        gijl(1,offL+1,1,1,js),ldi*nspc,(0d0,0d0),wk(offL+1,is,1,js,ik),ldi*nspc)
            if (nspc == 2) then
            call zgemm('T','N',norbi,ldj,norbi,(1d0,0d0),dpfr(1+offp,2,1+offp,is),mxorb*nspc,
     .          gijl(1,offL+1,2,1,js),ldi*nspc,(1d0,0d0),wk(offL+1,is,1,js,ik),ldi*nspc)
            endif
          enddo
          enddo
          offg = offH(2,1,ib) + offH(1,1,nbas+1) ! Offset to intermediate block of g for this site
          offp = offp + norbi   ! Offset to intermediate block for local dpfr
        enddo                   ! downfolding blocks
      enddo                     ! k-points
      enddo                     ! sites
C     call zprm('sqrt(pdot) g ',2,wk,ldi*nspc,ldi*nspc,ldj*nspc)

C --- Make g <- sqrt(pdot) g sqrt(pdot) ---
      do  ib = 1, nbas
C       if (s_site(ib)%ncomp > 1) cycle ! CPA case => g corresponds to g for first component
        if (mod3 == 2 .and. ib /= ib0) cycle
        icomp = 1
        if (s_site(ib)%ncomp > 1 .and. (mod0 == 0 .or. mod0 >= 5)) icomp = ijcomp(2)
        if (icomp < 1 .or.  icomp > s_site(ib)%ncomp) cycle
        norb = s_site(ib)%norb ; n2 = 2*norb
        lmr = mxorb*(ib-1)
        pfm => s_site(ib)%dpfr
        call pgfg2gnc2(100*idnf+0+mod0,s_pot,s_site,ib,icomp,norb,iprmb(1+lmr),pfv) ! Set up pfv
        isw = 10000*mod4 + 1000*mod2 + 100*idnf + mod0
C       ? Why is 1:norb*nspc*norb*nsp necessary?
        call pgfg2gnc(isw,isp,nsp,nspc,norb,mxorb,pfv,pfm(1:norb*nspc*norb*nsp,icomp),dpfr)
C       call zprm('sqrt(pdot)',2,dpfr,mxorb*nspc,norb*nspc,norb*nspc)

        do  ik = 1, nk          ! k-points

C   ... g <- wk sqrt(pdot)
        offg = offH(1,1,ib)     ! Offset to lower block of g for this site
        offp = 0                ! Offset to lower block for local dpfr
        do  idx = 1, 1+mod1
          norbi = offH(idx,1,ib+1)-offH(idx,1,ib)  ! Size of downfolding block
          if (norbi == 0) cycle
          offL = offg; offR = offg; if (mod3 == 1) offL = 0; if (mod3 == 2) offR = 0
          do  is = 1, nspc
          do  js = 1, nspc
            call zgemm('N','N',ldi,norbi,norbi,(1d0,0d0),wk(1,is,offR+1,1,ik),ldi*nspc,
     .        dpfr(1+offp,1,1+offp,js),mxorb*nspc,(0d0,0d0),gijl(1,1,is,offR+1,js),ldi*nspc)
            if (nspc == 2) then
            call zgemm('N','N',ldi,norbi,norbi,(1d0,0d0),wk(1,is,offR+1,2,ik),ldi*nspc,
     .          dpfr(1+offp,2,1+offp,js),mxorb*nspc,(1d0,0d0),gijl(1,1,is,offR+1,js),ldi*nspc)
            endif
            if (nk > 1) then
              call zcopy(ldi*norbi,gijl(1,1,is,offR+1,js),1,gij(ik,1,is,offR+1,js),nk)
            endif
          enddo
          enddo
          offg = offH(2,1,ib) + offH(1,1,nbas+1)  ! Offset to intermediate block of g for this site
          offp = offp + norbi   ! Offset to intermediate block for local dpfr
        enddo                   ! downfolding blocks

      enddo                     ! k-points

      enddo                     ! sites
C     call zprm('sqrt(pdot) g sqrt(pdot)',2,gijl,ldi*nspc,ldi*nspc,ldj*nspc)

C --- Make G <- pdotdot/pdot + sqrt(pdot) g sqrt(pdot) ---
      if (mod2 == 0) then      ! g ->> G
      do  ib = 1, nbas
        if (mod3 /= 0 .and. ib /= ib0) cycle
        icomp = 1
        if (s_site(ib)%ncomp > 1 .and. (mod0 == 0 .or. mod0 >= 5)) icomp = ijcomp(1)
        if (icomp < 1 .or.  icomp > s_site(ib)%ncomp) cycle
        norb = s_site(ib)%norb ; n2 = nspc*norb
        lmr = mxorb*(ib-1)
        pfm => s_site(ib)%ddpfr
        call pgfg2gnc2(100*idnf+10+mod0,s_pot,s_site,ib,icomp,norb,iprmb(1+lmr),pfv) ! Set up pfv
        isw = 10000*mod4 + 100*idnf + mod0
C       ? Why is 1:norb*nspc*norb*nsp necessary?
        call pgfg2gnc(isw,isp,nsp,nspc,norb,mxorb,pfv,pfm(1:norb*nspc*norb*nsp,icomp),dpfr)

C   ... Add pdotdot/pdot to g
        offg = offH(1,1,ib)     ! Offset to lower block of g for this site
        offp = 0                ! Offset to lower block for local dpfr
        do  idx = 1, 1+mod1
          norbi = offH(idx,1,ib+1)-offH(idx,1,ib)  ! Size of downfolding block
          if (norbi == 0) cycle
          offL = offg; offR = offg; if (mod3 == 1) offL = 0; if (mod3 == 2) offR = 0
          do  ik = 1, nk
          do  is = 1, nspc
          do  js = 1, nspc
            gij(ik,offL+1:offL+norbi,is,offR+1:offR+norbi,js) =
     .      gij(ik,offL+1:offL+norbi,is,offR+1:offR+norbi,js) + dpfr(offp+1:norbi,is,offp+1:norbi,js)
          enddo
          enddo
          enddo
          offg = offH(2,1,ib) + offH(1,1,nbas+1)  ! Offset to intermediate block of g for this site
          offp = offp + norbi   ! Offset to intermediate block for local dpfr
        enddo                   ! downfolding blocks
      enddo                     ! sites
      endif

      deallocate(pfv,dpfr,wk)
      if (nk > 1) deallocate(gijl)

C     call zprm('pdotdot/pdot + sqrt(pdot) g sqrt(pdot)',2,gij,ldi*nspc,ldi*nspc,ldj*nspc)

      call ztoyy(gij,nk*ldi*nspc,nk*ldi*nspc,ldj*nspc,ldj*nspc,1,kcplx)

      call tcx('gfg2gnc')

      end
      subroutine sugfg2gnc(lso,lrel,kcplx,lgtoG,nccomp,lint,mode)
C- Determine mode for gfg2gnc given specified input conditions
C ----------------------------------------------------------------------
Ci Inputs
Ci   lso   :T => SO coupling included perturbatively
Ci   lrel  :0 for non-relativistic
Ci          1 for scalar relativistic
Ci          2 for Dirac equation
Ci   kcplx :0 Im part of g follows real:  gij(ldi,nspc,ldj,nspc,1:2)
Ci         :1 g is in complex*16 format
Ci         :2 Im part of g follows real after first column:  gij(ldi,nspc,1:2,ldj,nspc)
Ci   lgtoG :0 scale g^alp to g^gam
Ci         :1 scale g to G
Ci         :2 scale G to g
Ci   nccomp:number of CPA components.  (0 => no CPA components; see s_ctrl%nccomp)
Ci   lint  :T if to include intermediate orbitals in the scaling
Co Outputs
Co   mode  :mode with which to call gfg2gnc
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   26 Mar 16  First created
C ----------------------------------------------------------------------
      implicit none
      logical lso,lint
      integer lrel,nccomp,mode,lgtoG,kcplx

C ... Non CPA case, g->G
      if (nccomp == 0 .and. lgtoG /= 0) then
        mode = 1; if (lso) mode = 2; if (lrel == 2) mode = 10002
C       Matrix form can also be used in SO and lrel cases
C       But not required since data should be stored in both s_pot%dpfr and s_site%dpfr
C       if (lso) mode = 5; if (lrel == 2) mode = 6

        if (lgtoG == 2) mode = mode+100  ! Reverse sense of scaling

C ... Non CPA case, galp -> ggam
      elseif (nccomp == 0 .and. lgtoG == 0) then
        mode = 3; if (lso) mode = 4; if (lrel == 2) mode = 10004

C ... CPA case, g -> G
      elseif (lgtoG > 0) then
        mode = 0; if (lso) mode = 5; if (lrel == 2) mode = 6

C ... CPA case, g -> ggam
      else
        call rx('sugfg2gnc: no implementation for CPA scaling galp -> ggam')
      endif

      if (lint) mode = mode+10
      mode = mode + 100000*kcplx

      end
      subroutine pgfg2gnc(mode,isp,nsp,nspc,norb,ldpf,pfv,pfm,pfout)
C- Kernel that copies potential function for one site to matrix form
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :
Ci         :0,1,3  input is pfv, in collinear form (two spin channels)
Ci         :2,4    input is pfv, in noncollinear form (four spin channels).  For SO or relativistic cases.
Ci         :5      input is pfm, in matrix form
Ci         :6      input is pfm, in matrix form
Ci         :7 (noncollinear) adds s_pot%gmar to diagonal and scales gij with s_pot%papg
Ci         :10s digit no longer used
Ci         :100s digit : ordering of orbitals (NOT IMPLEMENTED)
Ci         :0 return P in orbital order (iprmb is not used)
Ci         :1 return P in local downfolding order : iprmb order
Ci         :1000s digit : option to return inverse of P in pfout
Ci         :0 return P
Ci         :1 return P^-1
Ci   mode  :1s digit controls how s_site%pfr and s_site%dpfr are stored.
Ci         :0 in vector form, collinear spins i.e. pfr(1:norb,1:2)
Ci         :1 in vector form for spinors, i.e. pfr(1:norb,1:2,1:2)  (not implemented)
Ci         :2 in matrix form, lms repsn  pfr(1:norb,1:2,1:norb,1:2)
Ci         :3 in matrix form, kappa-mu repsn  pfr(1:norb,1:2,1:norb,1:2)
Ci   norb  :number of orbitals belonging to this site
Ci   ldpf  :dimensions pfout
Ci   pfv   :potential functions, in vector form
Ci   pfm   :potential functions, in matrix form
Co Outputs
Ci  pfout  :potential functions, in matrix form, lms repsn, downfolding order
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cb Bugs
Cb   Downfolding not implemented when P given in matrix form (modes 2 and 3)
Cb   For it to work, pgfg2gnc needs to permute rows and columns of P
Cb   iprml needs to be passed to this routine
Cu Updates
Cu   21 Jan 16 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,norb,ldpf,isp,nsp,nspc
      complex(8) :: pfv(norb,nsp),pfm(norb,nspc,norb,nsp),pfout(ldpf,nspc,ldpf,nsp)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wrk(:)
C ... Local parameters
      integer, parameter :: ncomp=1
      integer i,j,is,js,mod0,mod2,mod3,mod4,iprint !,iprml(norb)
!     complex(8) :: pfloc(norb,2,norb,2)

      mod0 = mod(mode,10)
C     mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      mod3 = mod(mode/1000,10)
      mod4 = mod(mode/10000,10)

      if (nspc == 2 .or. nsp == 1)
     .  call sanrg(.true.,isp,1,1,'pgfg2gnc:','isp')

C ... Collinear vector form
      if (mod0 <= 1 .or. mod0 == 3) then
C       Zero out relevant subblock of pfout
        do  i = 1, norb
          do  j = 1, norb
            pfout(i,:,j,:) = 0
          enddo
        enddo
        do  i = 1, norb
          pfout(i,1,i,1) = pfv(i,isp)
          if (nspc == 2) pfout(i,2,i,2) = pfv(i,2)
        enddo

C ... Noncollinear vector form
      elseif (mod0 == 2 .or. mod0 == 4) then
        call rx('gfg2gnc: check pokepf')
        call pokepf(10*mod4,1,norb,ldpf,0,norb,[0],0,1,1,nspc,nsp,pfv,pfout)

C ... Matrix form
      elseif (mod0 == 5 .or. mod0 == 6) then
        do  is = 1, 2
        do  js = 1, 2
          call zmcpy('N',pfm(1,is,1,js),norb*2,1,pfout(1,is,1,js),ldpf*2,1,norb,norb)
        enddo
        enddo
      else
        call rxi('gfg2nc : 1s digit mod0 not recognized : ',mod0)
      endif

C ... kappa-mu to lms
      if (mod0 == 6) then
        if (norb /= ldpf) call rx('mstokm not ready for generalized dimensioning')
C       call zprm('pkmu',2,pfout,ldpf*2,norb*2,norb*2)
        call mstokm(1,0,ncomp,norb,pfout,pfm,i)
C       call zprm('plms',2,pfout,ldpf*2,norb*2,norb*2)
      endif

C ... Permute to local downfolding order
      if (mod2 /= 0 .and. mod0 >= 5) then
        call rx('tidy up the permutation branch.')
C        call zprm('pfout',2,pfout,ldpf*2,norb*2,norb*2)
C        pfloc = pfout
C        do  i = 1, norb
C          do  j = 1, norb
C            pfout(i,:,j,:) = pfloc(iprml(i),:,iprml(j),:)
C          enddo
C        enddo
C       call zprm('pfout',2,pfout,ldpf*2,norb*2,norb*2)
      endif

C ... Return inverse of P
      if (mod3 == 1) then
        allocate(wrk(66*norb*nspc))
        call zqinv('n',pfout,ldpf*nspc,-66*norb*nspc,norb*nspc,wrk,norb*nspc,i)
        if (i /= 0) call rx1('gfg2gnc: matrix inversion failed, ierr=%i',i)
        deallocate(wrk)
      endif

      if (iprint() >= 100) call zprm('plms',2,pfout,ldpf*nspc,norb*nspc,norb*nspc)

      end

      subroutine pgfg2gnc2(mode,s_pot,s_site,ib,icomp,norb,iprmb,pfv)
C- Extracts potential function for one site, when given in vector form
C ----------------------------------------------------------------------
Cio Structures
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  dpf ddpf dpfr ddpfr palp pf gma papg gmar
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pf
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  dpfr ddpfr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :Controls what potential functions are returned
Ci         :Potential functions are returned in local downfolding order.
Ci         :Number of spin channels returned is inferred from size of structure used
Ci         :1s-10s digit  pfun returned for orbitals : connected to site ib
Ci         :0   sqrt(s_site(ib)%dpfr(:,icomp))  for component icomp
Ci         :10  s_site(ib)%ddpfr(:)
Ci         :1   sqrt(s_pot%dpf(:))
Ci         :11  s_pot%ddpf(:)
Ci         :2   s_pot%dpfr(:)
Ci         :12  s_pot%ddpfr(:)
Ci         :3   s_pot%palp(:)/s_pot%pgam(:)
Ci         :13  same as 3, but scaled by s_pot%gma
Ci         :4   s_pot%papg(:)
Ci         :14  s_pot%gmar(:)
Ci         :100s digit : ordering of orbitals
Ci         :0 return P in orbital order
Ci         :1 return P in local downfolding order : order of increasing iprmb(1:norb)
Ci   ib    :site index
Ci   norb  :number of orbitals belonging to this site
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co Outputs
Ci   pfv   :potential functions for one site, vector form, orbital or site-downfolding order
Cr Remarks
Cr   This routine bundles various sources of potential functions, returning
Cr   a vector form of pfun in a uniform style, suitable for generic processing
Cb Bugs
Cb   Should loop over mxorb for s_pot arrays, not norb!
Cu Updates
Cu   25 Mar 16 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,norb,ib,icomp,iprmb(*)
      complex(8) :: pfv(norb,*)
C ... For structures
!      include 'structures.h'
      type(str_pot) ::  s_pot
      type(str_site)::  s_site(*)
C ... Local parameters
      integer i,pfdim(2),nspin,mod01,mod2,iprml(norb)
      complex(8),pointer :: pfp(:,:),pfl(:,:)
      equivalence(nspin,pfdim(2))

      mod01 = mod(mode,100)
      mod2 = mod(mode/100,10)

C ... pfun drawing from a vector in either s_pot or s_site
      select case (mod01)
        case (0);  pfp => s_site(ib)%dpfr; goto 10
        case (10); pfp => s_site(ib)%ddpfr; goto 10 ! already ordered by l
        case (1);  pfp => s_pot%dpf
        case (11); pfp => s_pot%ddpf
        case (2);  pfp => s_pot%dpfr
        case (12); pfp => s_pot%ddpfr
        case (3, 13)
          pfdim(1:2) = shape(s_pot%pf)
          allocate(pfp(pfdim(1),pfdim(2)))
          pfp(:,:) = s_pot%palp(:,:)/s_pot%pf(:,:)
          if (mod01 == 13) pfp(:,:) = pfp(:,:)*s_pot%gma(:,:)
C          call zprm('palp',2,s_pot%palp(:,:),pfdim(1),pfdim(1),pfdim(2))
C          call zprm('pgam',2,s_pot%pf(:,:),pfdim(1),pfdim(1),pfdim(2))
C          call prmx('gma',s_pot%gma,pfdim(1),pfdim(1),pfdim(2))
        case (4);  pfp => s_pot%papg
        case (14); pfp => s_pot%gmar
        case default; return
      end select

C     Unwind downfolding order: re-order by l
      pfdim(1:2) = shape(pfp)
      allocate(pfl(norb,nspin))
      do  i = 1, norb
        if (mod01 <= 1) then
          pfl(i,1:nspin) = sqrt(pfp(iprmb(i),1:nspin))
        else
          pfl(i,1:nspin) = pfp(iprmb(i),1:nspin)
        endif
      enddo

      if (mod01 == 3 .or. mod01 == 13) deallocate(pfp)
      goto 20

C ... pfp -> pfl to (uncompressed orbital)*(spin dimension)
   10 continue
      pfdim(1:2) = shape(pfp); nspin = pfdim(1)/norb
      if (norb*nspin /= pfdim(1)) call rx('pgfg2gnc2 has dimensioning problem')
      allocate(pfl(norb,nspin))
      if (mod01 <= 1) then
        pfl = sqrt(reshape(pfp(:,icomp),(/norb,nspin/)))
      else
        pfl = reshape(pfp(:,icomp),(/norb,nspin/))
      endif

C ... Reorder orbitals in site-local downfolding order with intermediate block following lower block
   20 continue
      if (mod2 == 0) then
        do  i = 1, norb
          iprml(i) = i
        enddo
      else
        call ivheap(1,norb,iprmb,iprml,101) ! downfolding order of the orbitals within a site
      endif

      do  i = 1, norb
        pfv(i,1:nspin) = pfl(iprml(i),1:nspin)
      enddo

      deallocate(pfl)

      end
