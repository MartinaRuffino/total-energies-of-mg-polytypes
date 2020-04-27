C#define BLAS3
      subroutine secm2c(s_ctrl,s_str,ckbas,plat,nl,nsp,nbas,nclass,ipc,
     .  indxsh,qss,eula,neul,ikp,nkp,ldim,lidim,lihdim,ndim,qp,pp,sop,
     .  vmtz,wsr,avw,isp,pti,pph,nfilet,nevmx,efmax,nev,addsx,z,eb)
C- Set up and diagonalize ASA Hamiltonian, 2-center approximation.
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat,nl,nsp,nbas,lmx,pp,vmtz,indxsh
Ci   eula: Euler angles for spin rotations
Ci   ldim: dimension of l - wave block of hamiltonian matrix
Ci   nevmx:-2 Do not diagonalize, but copy h into z
Co Outputs:
Co   eigenvalues and eigenvectors are returned in eb, z
Co   pph : vector of pot. par's. in alpha rep'n; pp are returned
Co   in alpha representation.
Cr Remarks
Cr   Folding down not allowed here.
Cr   z must be available as a work space of size ldim**2*2, even if
Cr   evecs are not calculated.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      use structures
      implicit none
      integer nl,nsp,nbas,nclass,neul,isp,nfilet,ldim,lidim,lihdim,ndim,
     .  nevmx,nev,ipc(nbas),indxsh(ndim),ikp,nkp,addsx
      double precision ckbas,plat(3,3),eb(ldim*2),pp(6,nl,nsp,0:*),vmtz,
     .  avw,qp(3),z(ldim,ldim*2),pti(ndim,nsp),eula(nbas,neul,3),wsr(*),
     .  pph(5,lihdim,nsp),efmax,qss(4),sop(0:nl-1,nsp,nsp,9,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_str)::  s_str
C ... Dynamically allocated local arrays
      real(8), allocatable:: alpha(:),adot(:)
      real(8), allocatable :: ov(:)
      complex(8), allocatable :: h(:)
      complex(8), allocatable :: sii(:)
      complex(8), allocatable :: sil(:)
      complex(8), allocatable :: a(:)
      complex(8), allocatable :: sd(:)
      real(8), allocatable :: wk(:)
      real(8), allocatable :: gma(:)
      real(8), allocatable :: wk2(:)
      real(8), allocatable :: wk3(:)
      real(8), allocatable :: xsidt(:)
      real(8), allocatable :: ccd(:)
      real(8), allocatable :: diawk(:)
C local variables
      integer i,j,ipr,ii,lgunit,nsite,i2,li,lncol,lham,ld2,ld22,n,
     .  bit,nlspc,nl2,i1mach,linv,ldimx,idim,lbloch,isw
      logical ndnfld,ccor,lx,lss,lnc,lso,bittst,ltmp,iostr
      double precision kap2(20),qpq(3),ddot,ercc,elin,xx
C     double precision Ze(2),Gij(2,ldim,ldim)
      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) == bit)

C --- Setup ---
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      elin = s_ctrl%elin
      call getpr(ipr)
      nl2 = nl**2
      idim = lidim - ldim
      if (idim /= 0) call fexit(-1,111,'  Exit -1 SECM2C: '//
     .  '%i intermediate waves sought but none allowed',idim)
      nlspc = nl * nsp * nclass
C     Local copy of alpha,adot, so they can be altered
      allocate(alpha(size(s_str%alph)))
      call dcopy(size(s_str%alph),s_str%alph,1,alpha,1)
      allocate(adot(size(s_str%adot)))
      call dcopy(size(s_str%adot),s_str%adot,1,adot,1)

C     Noncollinear switches
      lnc = lncol /= 0
      lss = bittst(lncol,2)
      lso = bittst(lncol,4)

C     Possibly rotate to spherical harmonics when making Bloch sum
      lbloch = 0
      if (bittst(lham,256)) lbloch = 1000

C     Combined correction required if downfolding
      ndnfld = idim == 0
C     ccor = .not. ndnfld .or. lgors('ctrl lasa,4',sctrl)
      ccor = .not. ndnfld .or. IAND(s_ctrl%lasa,4) /= 0
C     if (.not. ccor  .and.  ipr >= 30  .and.  ikp == 1)
C    . print *, 'SECM2C : Combined Correction switched off'
C     if (.not. ndnfld .and. bittst(lham,128))
C    .  call rx('secmat: no downfolding in gamma rep')

C     Diagonalize by inverse iteration, or not
      linv = 0
C     if (nevmx > 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
      if (nevmx > 0 .and. IAND(s_ctrl%lqp,2) /= 0) linv = 1

C ... Sanity checks
      if (.not. ndnfld .and. bittst(lham,128))
     .  call rx('SECMAT: no downfolding in gamma rep')
      if (lnc .and. .not. ndnfld)
     .  call rx('noncollinear magnetism not implemented with dnfolding')
      call sanrg(.true.,addsx,0,0,'secm2c:','addsx')
      call sanrg(.true.,isw(lso),0,0,'secm2c:','spin-orbit')


C     Some dimensioning parameters and memory allocation
      ldimx = ldim
      ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      ld22 = 2*ld2
      if (ndnfld) li = 1
      if (ndnfld) i2 = 1
      call zinit(z,ldimx**2)

      allocate(sii(i2),sil(li),a(li),wk(2*lihdim))
      allocate(h(ldimx**2)); call dpzero(h,2*ldimx**2)
      call dcopy(3,qp,1,qpq,1)

C --- Get screened sdot from disc and Bloch-transform ---
      allocate(sd(1))
      if (ccor) then
        if (lss) then
          deallocate(sd)
          allocate(sd(ld22))
        endif
        if (lnc .and. .not. lss)  then
          deallocate(sd)
          allocate(sd(ld2))
        endif
        nsite = s_str%npr(nbas+1)
        if (lss) then
          qpq(1) = qp(1) + qss(1)/2
          qpq(2) = qp(2) + qss(2)/2
          qpq(3) = qp(3) + qss(3)/2
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,sd,xx,xx)
          if (ipr >= 110)
     .    call yprm('Sd(2)',12,sd,ldim*ldim,ldim,ldim,ldim)
          call pvsec2(ld22,sd)
          qpq(1) = qp(1) - qss(1)/2
          qpq(2) = qp(2) - qss(2)/2
          qpq(3) = qp(3) - qss(3)/2
        endif
        if (lss .or. lnc) then
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,sd,xx,xx)
          if (ipr >= 110)
     .    call yprm('Sd(1)',12,sd,ldim,ldim,ldim,ldim)
        else
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .      s_str%sdot,nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,z,xx,xx)
          call dscal(ld22,-avw**2,z,1)
        endif
        allocate(ccd(3*lihdim))
      else
        allocate(ccd(1))
      endif

C --- Get screened strux from disc and Bloch-transform them ---
      nsite = s_str%npr(nbas+1)
      if (lss) then
        qpq(1) = qp(1) + qss(1)/2
        qpq(2) = qp(2) + qss(2)/2
        qpq(3) = qp(3) + qss(3)/2
        call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .    s_str%s,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z(1+ld22,1),
     .    xx,xx)
        qpq(1) = qp(1) - qss(1)/2
        qpq(2) = qp(2) - qss(2)/2
        qpq(3) = qp(3) - qss(3)/2
        if (ipr >= 110)
     .    call yprm('Sll(2)',12,z(1+ld22,1),ldim*ldim,ldim,ldim,ldim)
      endif
      call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,s_str%iax,
     .  s_str%s,nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,h,xx,xx)

C --- Transform to gamma representation ---
      if (bittst(lham,128)) then
C   ... s,sdot are correctly rotated but ccor should be done in pert th
        if (ccor) call rx('secm2c: ccor not working for Gamma repsn')
C   ... Put pp in gamma; make vector of gamma's
        allocate(gma(ldim))
        allocate(ov(nlspc))
        i = -3
        if (isp == 2) i = -4
        if (lnc) i = -5
        call pptrns(i,nl,ipc,nclass,nsp,gma,nbas,pp,ov)
        deallocate(ov)
        call makpph(nl,nsp,nbas,lihdim,ipc,indxsh,pp,pph)
C       call prm('gamma',0,gma,lihdim,lihdim,1)
C       call prm('alpha',0,alpha,lihdim,lihdim,1)
C   ... put gamma-alpha in gam, gamma in alpha
        call daxpy(ldim,-1d0,alpha,1,gma,1)
        call daxpy(ldim, 1d0,gma,1,alpha,1)
C   ... Transform sll = h and sdot into gamma
        allocate(wk3(ld22))
        allocate(xsidt(ldim)); call dpzero(xsidt,ldim)
C   ... Make adot so that 3C contribution to cc is exactly zero
*       call prmx('gamma-alpha',gma,ldim,ldim,1)
*       call prmx('adot',adot,ldim,ldim,1)
        if (lss) then
C     ... Use copies of gamma-alpha, xsidot since mksbet overwrites them
          call dpcopy(gma,wk,1,ldim,1d0)
          allocate(wk2(ldim))
          call dpcopy(xsidt,wk2,1,ldim,1d0)
C     ... swap 1st and 2nd B.T. sd, to avoid more workspace
          if (ccor) call pvsec2(ld22,sd)
          call mksbet(ccor,ldim,wk,wk2,z(1+ld22,1),sd,wk3)
          if (ipr >= 110) then
            if (ccor) call yprm('Sd(2) (gam)',12,sd,ldim*ldim,ldim,
     .        ldim,ldim)
            call yprm('Sll(2)(gam)',12,z(1+ld22,1),ldim*ldim,ldim,ldim,
     .        ldim)
          endif
          if (ccor) call pvsec2(ld22,sd)
          deallocate(wk2)
        endif
        call mksbet(ccor,ldim,gma,xsidt,h,sd,wk3)
        if (ipr >= 110) then
          if (ccor) call yprm('Sd (gam)',12,sd,ldim*ldim,ldim,ldim,
     .      ldim)
          call yprm('Sll (gamma)',12,h,ldim*ldim,ldim,ldim,ldim)
        endif
        deallocate(gma,wk3,xsidt)
C --- Or else transform potential parameters to alpha rep'n ---
      else
        allocate(ov(nlspc))
        call pptrns(0,nl,ipc,nclass,nsp,alpha,nbas,pp,ov)
        deallocate(ov)
      endif

C ... Copy h into z for noncoll case
      if (lnc .or. lss) then
        call dpcopy(h,z,1,ld22,1d0)
      endif

C ... Scale sdot.  All scalings should be unified
      if (ccor .and. size(sd) > 1)
     .  call dscal((ldimx/ldim)*ld22,-avw**2,sd,1)

C --- Diagonal matrices for combined correction---
      if (ccor) then
        call makdia(nl,nbas,lihdim,indxsh,lidim,ipc,wsr,avw,alpha,
     .    adot,ccd)
        call dscal(3*lihdim,avw**2,ccd,1)
        ercc = (ddot(3*ldim,ccd,1,ccd,1) -
     .    ddot(2*ldim,ccd,1,ccd,1))/ldim
        if (ercc > 1d-10 .and. ipr >= 10 .and. ikp == 1)
     .    call awrit1(' SECM2C (warning): 3-C CC should be zero '//
     .    'but has rms %1,5;5d',' ',80,lgunit(1),dsqrt(ercc))
      endif

C --- Non-collinear two-center hamiltonian ---
      if (nsp == 2 .and. (lss .or. lnc)) then
        allocate(wk2(ldim*2))
        call hml2nc(nbas,nl,indxsh,qss,eula,neul,pph,
     .    ccor,lss,lnc,ccd,wk2,vmtz,elin,ldim,z,sd,h)
        deallocate(wk2)
C --- Collinear two-center hamiltonian ---
      else
C    ... ASA 2-center hamiltonian + ccor + terms of order o*enu ---
        if (.true. .or. ccor) then
          call hmlt2c(ccor,ccd,vmtz,elin,ldim,lihdim,pph(1,1,isp),
     .      h,z,wk)
        else
C    ... Simple ASA 2-center hamiltonian; uses no work space
         call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph(1,1,isp),h,h)
         call daxpy(ldim,1d0,pph(2,1,isp),5,h,ldim+1)
        endif
      endif
      if (ipr >= 110)
     .  call yprm('H',12,h,ldimx*ldimx,ldimx,ldimx,ldimx)

C --- Return with hamiltonian in z if nevmx is -2 ---
      if (nevmx == -2) then
        call dcopy(ldimx**2*2,h,1,z,1)
        return
      endif

C --- Eigenvalues and eigenvectors of 2C Hamiltonian ---
C#ifdef BLAS3
      lx = .true.
C#elseC
C      lx = .false.
C#endif
      if (linv /= 0) then
        allocate(diawk(ldimx*11))
      else
C ...   Use 5*ldim for parallel implementations ...
        allocate(diawk(ldimx*5))
      endif
C      call prmx('h',h,ldim,ldim,ldim)
      call diagno(ldimx,h,xx,diawk,lx,.false.,linv,nevmx,efmax,nev,z,eb)
      deallocate(diawk)

C     call prmx('evl',eb,ldim,ldim,1)

C  999 print *, 'zp= ?'
C      read (*,*) ze
C      call ev2cgf(isp,ldim,ldim,ldim,nev,eb,z,ze,Gij)
C      goto 999

      if (ipr >= 30) then
        j = min(9,ldimx)
        if (ipr >= 35) j = ldimx
C#ifdefC LINUX_PGI
C        do  18  ii = 1, 1
C#else
        do  18  ii = 1, 2
C#endif
        call awrit3(' SECM2C:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
   18   write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (ipr >= 36 .and. nev > 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,i1mach(2),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif
      if (ipr >= 110)
     .  call yprm('evec',2,z,ldimx*ldimx,ldimx,ldimx,ldimx)

      deallocate(sii,sil,a,wk,h,sd,ccd)
      if (ipr >= 110) call query('V<110 to skip matrix printing',-1,0)
      end
C      subroutine ev2cgf(isp,ni,nj,ldim,nbmx,eval,evc,ze,Gij)
CC- Make GF from eigenvectors
C      implicit none
C      integer isp,ldim,nbmx,ni,nj
C      double precision evc(ldim,ldim,2),eval(nbmx,isp)
C      double complex Gij(ni,nj),Ze
C
C      integer i,j,nu
C      double complex wk
C
C      do  20  i = 1, ni
C      do  20  j = 1, nj
C        wk = 0
C        do  22  nu = 1, ldim
C          wk = wk + (dcmplx(evc(i,nu,1),evc(i,nu,2)) *
C     .               dcmplx(evc(j,nu,1),-evc(j,nu,2))) /
C     .              (Ze-eval(nu,isp))
C   22   continue
C        Gij(i,j) = wk
C
C   20 continue
C
C      call zprm('gf',2,Gij,ldim,ldim,ldim)
C
C      end
C
