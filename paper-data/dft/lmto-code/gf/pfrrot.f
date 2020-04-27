      subroutine rotpfr(mode,eula,nl,iprmb,ldimp,ldima,offp,pfdim,pfr,lds,kcplx,ps)
C- Expand, rotate (vector) relativistic Pot. functions for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         : 0 normal rotation
Ci         : 1 reverse the sense of rotation
Ci         :10s digit
Ci         : 0 real harmonics
Ci         : 1 spherical harmonics
Ci         :100s digit
Ci         : 0 initialize ps to zero
Ci         : 1 add to ps instead of initializing ps to zero
Ci   eula  :Euler angles alpha,beta,gamma
Ci   nl    :(maximum l) + 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci         :if iprmb(1) = 0, normal ordering is used
Ci   ldimp :(used only with downfolding)
Ci         :orbitals with index outside (ldimp+1,ldima) are not copied to ps
Ci   ldima :see ldimp
Ci   offp  :offset in pfr to start of potential fns for this site
Ci         :offp = 0 for first site, nl**2 for 2nd site, etc.
Ci   pfdim :leading dimension of pfr
Ci   pfr   :Relativistic potential functions
Ci         :pfr is diagonal in L, but has 2x2 matrix structure in spin
Ci         :pfr are ordered in one of the following ways:
Ci         :Case iprmb(1)=0 : lm order
Ci         :Case iprmb(1)>0 : downfolding order
Ci   lds   :leading dimension of ps.  lds must be at least nl**2.
Ci   kcplx :s is returne in the following format:
Ci          0: double complex with imaginary following real
Ci          1: complex*16
Ci          2: double complex with imaginary following real in columns
Co Outputs
Co    ps   :pfr (diagonal in L) is expanded to array ps(L,is,L',js)
Co         :and rotated through an angle specified by Euler angles.
Co         :Array ps is always returned in orbital order
Cr Remarks
Cu Updates
Cu   18 Jun 18 Supersede call to pfr2block with call to new pokepf
Cu   13 Aug 13 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,ldimp,ldima,offp,pfdim,lds,kcplx
      integer iprmb(*)
      double precision eula(3)
      double complex pfr(pfdim,2,2),ps(lds,2,lds,2)
C ... Local parameters
      integer mod2,morder,nsp,nspc

      nsp = 2; nspc = 2
      call lmorder(0,morder,mod2,mod2)
      mod2 = mod(mode/100,10)

      call rx('rotpfr: check rotation convention in call to pokepf')
C     call pfr2block(mod2,nl,iprmb,ldimp,ldima,offp,pfdim,pfr,lds,1,ps)
      call pokepf(mod2+10*morder,1,pfdim,lds,0,ldima,iprmb,offp,0,0,nspc,nsp,pfr(1+offp,1,1),ps)
c     call ztoyy(ps,18,18,18,18,1,2)
      call rot_LandS(mode,eula,nl,lds,kcplx,ps)

C     Debugging
C     if (kcplx == 1) then
C     call zprm('ps after spinor + orbital rotation',2,
C    .  ps,lds*2,lds*2,lds*2)
C     endif

      end

      subroutine rot_LandS(mode,eula,nl,lds,kcplx,s)
C- Rotate one (L,1:2,L',1:2) block of a spinor, both orbital and spin parts
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         : 0 normal rotation
Ci         : 1 reverse the sense of rotation
Ci         :10s digit
Ci         : 0 real harmonics
Ci         : 1 spherical harmonics
Ci         :100s digit
Ci         : 0 normal rotation
Ci         : 1 skip spinor part
Ci         : 2 skip orbital part
Ci   eula  :Euler angles alpha,beta,gamma
Ci   nl    :(maximum l) + 1
Ci   lds   :leading dimension of ps.  lds must be at least nl**2.
Ci   kcplx :s is returne in the following format:
Ci          0: double complex with imaginary following real
Ci          1: complex*16
Ci          2: double complex with imaginary following real in columns
Co Outputs
Co   s     : U+ s  U   if mode=0
Co         : U  s  U+  if mode=1
Co         : where U is the product of spinor, orbital rotations
Co         : s is returned in orbital order
Cl Local variables
Cl         :
Cr Remarks
Cr   Check of consistent definitions of spinor and orbital parts.
Cr   rotation z:pi/2 should rotate both spinor and orbital parts passively, i.e.
Cr   px -> -py, and py -> px,  sigx -> -sigy, and sigy -> sigx.
Cr   It does this correctly for the spinor part, because:
Cr     rotspu adopts an "active" convention but rotsp1 performs the inverse
Cr     operation U+ S U, making the net rotation passive.
Cr   It does this correctly for the orbital part (standard definition for orbital part)
Cu Updates
Cu   13 Aug 13 first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,lds,kcplx
      double precision eula(3)
      double complex s(lds,2,lds,2)
C ... Local parameters
      integer i,nl2,modspn,is,js,lsph,isw,optrot,roth,istab,ldsx
      integer n0H,nkap0
      parameter (nkap0=4,n0H=5)
      integer iprmb(nl*nl),offH(n0H,nkap0,2)
      double precision g(3,3),eulal(3)
      double complex u(2,2,nl*nl),xx,swk(lds,2,lds,2)
      double complex sdum(lds,lds)
      real(8),parameter::  pos(3)=0d0, ag(3)=0d0, q0(3)=0d0
      integer kcplxi
C     double precision rmat(nl**2,nl**2)

      if (kcplx /= 1) call rx('rot_LandS configured for kcplx=1 only')
      istab = 1
      ldsx = lds*2
      lsph = isw(mod(mode/10,10) /= 0) ! 0 => real harmonics, 1 => spherical harmonics
      nl2 = nl*nl
      forall (i = 1:nl2) iprmb(i) = i
C     Reverse sense of rotation: swap alpha with gamma and scale all angles by -1
      if (mod(mode,10) /= 0) then
        eulal(1) = eula(3)
        eulal(2) = eula(2)
        eulal(3) = eula(1)
        call dscal(3,-1d0,eulal,1)
      else
        eulal = eula
      endif

C ... Spinor part
      if (mod(mod(mode/100,10),2) == 0) then ! spin rotation

C     Spinor rotation matrix, U+ used in rotsp1
      call rotspu(0,1,1,1,nl,eulal,1,u) ! Spinor rotation matrices

C     Debugging: Copy spinor rotation into (LS,L'S') matrix
C      print *, 'euler', sngl(eulal)
C      s = 0
C      do  i = 1, nl2
C        do  is = 1, 2
C        do  js = 1, 2
C          s(i,is,i,js) = u(is,js,i)
C        enddo
C        enddo
C      enddo
C      call zprm('U+',2,s,ldsx,ldsx,ldsx)

C ... Spinor rotation: 2x2 rotation for each (L,L')
      modspn = 10000 + 10 ! c*16 array, input in spinor form
C     call zprm('S before spinor rotation',2,s,ldsx,ldsx,ldsx)
      call rotsp1(modspn,1,1,1,nl,iprmb,u,xx,xx,1,nl2,nl2,lds,
     .  lds,xx,xx,xx,s,s,s)
C     call zprm('S after spinor rotation',2,s,ldsx,ldsx,ldsx)
      endif

C ... Orbital part
      if (mod(mode/100,10) < 2) then ! spin rotation

C ... Orbital rotation LxL rotation for each (S,S')
C     setup to call roth
      offH(:,1,1) = 0
      offH(:,1,2) = nl2 ! 'last+1' site
      offH(2:5,1,1) = nl2 ! No intermediate or higher blocks
      call eua2rm(eulal(1),eulal(2),eulal(3),g) ! Rotation matrix from Euler angles

C     Debugging: copy matrix rotating Ylm into (LS,L'S') matrix (real harmonics)
C      s = 0
C      do  is = 1, nl2
C      do  js = 1, nl2
C        s(is,1,js,1) = rmat(is,js)
C        s(is,2,js,2) = rmat(is,js)
C      enddo
C      enddo
C      call zprm('rl',2,s,ldsx,ldsx,ldsx)

C     roth works internally with real,imaginary separated.
C     Could speed things up a little if we do the splitting beforehand.
C     Doesn't work ... not sure why
C     kcplxi = 2; call ztoyy(s,ldsx,ldsx,ldsx,ldsx,1,kcplxi)
      kcplxi = 1
      optrot = 100*kcplxi + 20*lsph + 1 ! , possibly spherical harmonics
      do  is = 1, 2
        do  js = 1, 2
          sdum(:,:) = s(:,is,:,js)
          if (roth(optrot,nl,1,pos,xx,offH,iprmb,istab,g,ag,q0,lds,
     .      lds,swk,sdum) < 0) call rx('bug in roth')
          s(:,is,:,js) = sdum(:,:)
        enddo
      enddo
C     call zprm('S after spinor + orbital rotation',2,s,ldsx,ldsx,ldsx)
C ... Put s in expected format
C     call ztoyy(s,ldsx,ldsx,ldsx,ldsx,1,kcplx)

      endif

      end

      subroutine pfr2block(mode,nl,iprmb,ldimp,ldima,offp,pfdim,pfr,lds,kcplx,ps)
C- Expand (vector) relativistic Pot. function to array for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 initialize ps to zero before adding P's
Ci         :1 add P's to ps passed as input
Ci         :10s digit
Ci         : 0 poke into matrix for one site
Ci         : 1 poke into a large matrix in downfolding order (needs iprmb)
Ci         :100s digit
Ci         : 0 use pfr stored in the format of mkfrpf.f (2 x 2)
Ci         : 1 use
Ci   nl    :(maximum l) + 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci         :if iprmb(1) = 0, normal ordering is used
Ci   ldimp :(used only with downfolding)
Ci         :orbitals with index outside (ldimp+1,ldima) are not copied to ps
Ci   ldima :see ldimp
Ci   offp  :offset in pfr to start of potential fns for this site
Ci         :offp = 0 for first site, nl**2 for 2nd site, etc.
Ci   pfdim :leading dimension of pfr
Ci   pfr   :Relativistic potential functions (mkfrpf.f) in lms repsn
Ci         :pfr is diagonal in L, but has 2x2 matrix structure in spin
Ci         :pfr are ordered in one of the following ways:
Ci         :Case iprmb(1)=0 : lm order
Ci         :Case iprmb(1)>0 : downfolding order
Ci   lds   :leading dimension of ps.  lds must be at least nl**2.
Co Outputs
Co    ps   :vector pfr is expanded to array ps(L,is,L',js)
Co         :Array ps is always returned in orbital order
Cb Bugs
Cb   This routine should be merged with pokepf
Cu Updates
Cu   05 Jun 18 Redesigned with consistent index ordering
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,ldimp,ldima,offp,pfdim,lds,kcplx
      integer iprmb(*)
      double complex pfr(pfdim,2,2),ps(lds,2,lds,2)
C ... Local parameters
      logical ldnfld
      integer mod0,mod1,i,l,ilm,m,lmi,morder,im
C ... External calls
      external dpzero,ztoyy

      call rxx(lds<nl**2,'pfr2block: bad dimension, lds')

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      morder = mod(mode/100,10)
      im = -1; if (morder==1 .or. morder==2) im = 1

      if (mod0 == 0) call dpzero(ps,lds*2*lds*2*2)
      ldnfld = iprmb(1) /= 0

      if (mod1 == 1 .and. .not. ldnfld)
     .  call rx('pfr2block: mode1=1 and iprmb=0')

C ... Copy pfr to ps, spin diagonal parts for the lower+intermediate block
      lmi = offp
      ilm = 0
      do  l = 0, nl-1
      do  m = -l, l
        lmi = lmi+1
        ilm = ilm+1
        if (ldnfld) then
          i = iprmb(lmi)
        else
          i = lmi
        endif
        if (mod1 == 1) ilm = i
        if (i > ldimp .and. i <= ldima .or. .not. ldnfld) then
C         print *, l,m,ilm,i,cmplx(pfr(i,1,2))
          ps(ilm,1,ilm,1) = ps(ilm,1,ilm,1) + pfr(i,1,1) ! (1,1) block
          ps(ilm,2,ilm,2) = ps(ilm,2,ilm,2) + pfr(i,2,2) ! (2,2) block
          if (ilm+im < 1 .or. ilm+im > lds) cycle
          ps(ilm+im,2,ilm,1) = ps(ilm+im,2,ilm,1) + pfr(i,2,1)    ! (2,1) block
          ps(ilm,1,ilm+im,2) = ps(ilm,1,ilm+im,2) + pfr(i+im,1,2) ! (1,2) block
        endif
      enddo
      enddo
C     call zprm('P in s format',2,ps,lds*2,lds*2,lds*2)

C     stop

C ... Put ps in expected format
      if (kcplx /= 1) then
        i = lds*2
        call ztoyy(ps,i,i,i,i,1,kcplx)
      endif

      end
C      subroutine fmain
CC- Tests conventions for rot_LandS
C      implicit none
C      integer nl,lds,i,j,is,js
C      double complex srm1,sig(2,2),sigy(2,2),s(4,2,4,2),sr(4,2,4,2)
C      double precision g(3,3),eula(3),pi,xx
C      procedure(real(8)) :: dglob
C
C      pi = 4*datan(1d0)
C
C      srm1 = dcmplx(0d0,1d0)
C      sigy = 0 ; sigy(1,2) = -srm1 ; sigy(2,1) = srm1
C
C      print *, '+pi/2 rotation around z, spinor part only : rotate sigy * S00, show sigy->sigx'
C
C      call a2rotm('z:pi/2',.false.,40,g)
C      call rm2eua(g,eula(1),eula(2),eula(3))
C      call info5(1,0,0,' Euler angles from g: eula(1)%;12,8D   eula(2)%;12,8D  eula(3)%;12,8D',eula(1),eula(2),eula(3),4,5)
C
C      nl = 1; lds = nl*nl
C      xx = dglob('nl',dble(nl),1)
C      xx = dglob('nkaph',1d0,1)
C      sig = sigy
C      call rot_LandS(200,eula,nl,lds,1,sig)
C
C  333 format(2f12.6,4x,2f12.6)
C      call info0(1,0,0,'%12fsigy%23frotated')
C      do j = 1, 2
C        print 333, (dble(sigy(j,i)), i=1,2),  (dble(sig(j,i)), i=1,2)
C      enddo
C      print *
C      do j = 1, 2
C        print 333, (dimag(sigy(j,i)), i=1,2),  (dimag(sig(j,i)), i=1,2)
C      enddo
C
C      nl = 2; lds = nl*nl
C      xx = dglob('nl',dble(nl),1)
C      call info0(1,1,-1,' pi/2 orbital-only rotation around z, s+p orbitals.  S = S_yz...')
C      call info0(1,0,0,' show that S_yz -> +S_xz')
C      s = 0
C      s(3,1,2,1) = 1
C      s(3,2,2,2) = 1
C      sr = s
C      call rot_LandS(100,eula,nl,lds,1,sr)
C
C  334 format(4f12.6,2x,4f12.6)
C      call info0(1,1,0,' Initial s (real part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dble(s(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
CC      print *
CC      do js = 1, 2
CC      do j = 1, 4
CC        print 334, ((dimag(s(j,js,i,is)), i=1,4), is=1,2)
CC      enddo
CC      enddo
C
C      call info0(1,0,0,' Rotated s (real part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dble(sr(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
C
C      call info0(1,1,-1,' Spin+orbital rogation of sigy*S_00 + 1*S_yz')
C      s = 0
C      s(1,:,1,:) = sigy
C      s(3,1,2,1) = 1
C      s(3,2,2,2) = 1
C      sr = s
C      call rot_LandS(000,eula,nl,lds,1,sr)
C
C      call info0(1,0,0,' Initial s (real part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dble(s(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
C      call info0(1,0,0,' Initial s (Imaginary part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dimag(s(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
C
C      call info0(1,0,0,' Rotated s (real part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dble(sr(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
C      call info0(1,0,0,' Rotated s (Imaginary part)')
C      do js = 1, 2
C      do j = 1, 4
C        print 334, ((dimag(sr(j,js,i,is)), i=1,4), is=1,2)
C      enddo
C      print *
C      enddo
C
C      end
