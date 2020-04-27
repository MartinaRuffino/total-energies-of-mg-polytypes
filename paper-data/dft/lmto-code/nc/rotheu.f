      subroutine rotheu(opt,s_ham,nl,ib1,ib2,nbas,ldmpa,ldima,ldha,ldhb,h)
C- Rotate subblock of h by a given set of site-dependent Euler angles
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         : 0 normal rotation
Ci         : 1 reverse the sense of rotation
Ci         :10s digit
Ci         : 0 real harmonics
Ci         : 1 spherical harmonics
Ci         :100s digit distinguishes how complex arithmetic is handled
Ci         :   0: h,h0 have real, imaginary separated
Ci         :      h = h(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci         :   1: h,h0 are in complex*16 format (see Bugs)
Ci         :      h = h(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci         :   2: h,h0 have real, imaginary separated by columns
Ci         :      h = h(ldha,2,ldhb), with h(*,1..2,*) = real..imag
Ci         :1000s digit if true, suppresses orbital rotation
Ci         :   0: both orbital and spinor rotation
Ci         :   1: skip spinor part
Ci         :   2: skip orbital part
Ci   nl    :(global maximum l) + 1
Ci   ib1
Ci   ib2
Ci   ldmpa :offset to first orbital in downfolding block of h to rotate
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h
Cio Inputs/Outputs
Cio   h    :On input, unrotated hamiltonian.  On output, rotated hamiltonian
Cl Local variables
Cl   rmat  :rotation matrix.
Cl         :rmat is real for real harmonics, complex for spherical harmonics
Cl   ldhr  :leading dimension of h
Cl   ofhi  :offset to imaginary part of h
Cr Remarks
Cr   Clarification of the sense of rotation. The following g corresponds to
Cr   a 90 degree counterclockwise rotation of the coordinate system
Cr   generated by rotation: z:pi/2  (right-hand rule)
Cr            ( 0   1   0)
Cr        g = (-1   0   0)
Cr            ( 0   0   1)
Cr   The x-axis xhat gets mapped into yhat, the y-axis mapped into -xhat.  Alternatively,
Cr   a vector initially parallel to xhat (yhat) becomes parallel to -yhat (+xhat).
Cr   This is the "intrinsic" sense defined in Wikipedia (en.wikipedia.org/wiki/Euler_angles)
Cr   Euler angles are alpha=0, beta=pi/4, gamma=0
Cr
Cr   This g generates the following pp block of the rotation of real harmonics:
Cr            (0    0   -1)
Cr     R_pp = (0    1    0)
Cr            (1    0    0)
Cr   (Note Y_1-1,Y_10,Y_11 correspond to y,z,x respectively.)
Cr   R_pp maps  px -> -py, and py -> px.    In other words,  a vector of
Cr   coefficients to a linear combination of Y_1-1,Y_10,Y_11 is transformed as
Cr
Cr            (C_1-1)   (-C_11)
Cr       R_pp (C_10 ) = (C_10)
Cr            (C_11 )   (C_1-1)
Cr
Cr   The p block of structure matrix S consisting only of S_yz -> S_xz
Cr    (0    0    0)    (0    0    0)          (0    0    0)    ( 0    0    0)
Cr    (1    0    0) -> (0    0    1)  whereas (0    0    1) -> (-1    0    0)
Cr    (0    0    0)    (0    0    0)          (0    0    0)    ( 0    0    0)
Cr
Cr   Note consistent definitions of spinor and orbital parts.  rotation z:pi/2
Cr   should rotate px -> -py, and py -> px,  and rotate sigx -> -sigy, and sigy -> sigx.
Cb Bugs
Cb  *kcplx=1 not tested
Cu Updates
Cu   28 Apr 17 First created.  A generalization of rot_LandS to hamiltonian blocks
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, parameter :: nkap0=4,n0H=5,nspc=2
      integer opt,ldha,ldhb,nl,ib1,ib2,nbas,ldmpa,ldima
      real(8) :: h(ldha,nspc,ldhb,nspc,2)
C ... For structures
!     include 'structures.h'
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      real(8),allocatable :: eulal(:,:),hwk(:,:,:,:,:)
C ... Local parameters
      integer ib,is,js,kcplx,ldh1,ldh2,ldhax,ldhbx,
     .  ldhr,ldimb,ldmpb,lsph,nl2,ofhi,ncsw,i1,i2,i3,optrot,morder
      double precision pi,q0(3)
      double precision wk(nl**4),rmat(nl**4*2),g(3,3)
      procedure(integer) :: isw,iprint
      procedure(real(8)) :: ddot,dasum

      call tcn('rotheu')

      q0 = 0
      ldhax = ldha*nspc
      ldhbx = ldhb*nspc
      kcplx = mod(opt/100,10)
      nl2 = nl*nl
      pi = 4d0*datan(1d0)
      lsph = isw(mod(opt/10,10)/=0) ! 0 => real harmonics, 1 => spherical harmonics
      call cplxdm(kcplx,ldhax,ldhbx,ldh1,ldh2,ldhr,ofhi) ! ldhr,ofhi for h
      call lmorder(0,morder,[0],[0])

C      if (iprint() >= 50) then
C        call yprm('starting h',kcplx+2,h,ofhi,ldhax,ldhax,ldhbx)
C      endif

      if (kcplx == 1) then
        call rx('rotheu: check kcplx=1')
        kcplx = 2
        call ztoyy(h,ldhax,ldhbx,ldhax,ldhbx,1,2)
      endif

C --- Reverse sense of rotation: swap alpha with gamma and scale all angles by -1 ---
      allocate(eulal(ib2,3))
      do  ib = ib1, ib2
        if (mod(opt,10) /= 0) then
          eulal(ib,1) = -s_ham%eula(ib,3)
          eulal(ib,2) = -s_ham%eula(ib,2)
          eulal(ib,3) = -s_ham%eula(ib,1)
        else
          eulal(ib,1) = s_ham%eula(ib,1)
          eulal(ib,2) = s_ham%eula(ib,2)
          eulal(ib,3) = s_ham%eula(ib,3)
        endif
      enddo

C --- Spinor rotation ---
      if (mod(opt/1000,10) /= 1) then  ! Execute branch only if include spinor part

        ncsw = 30000 + 000 + 10*kcplx
        if (ldhb /= ldha) call rx('implement spinor rotation for off-diagonal h')
        call rotspn(ncsw,ib1,ib2,ib2,nl,s_ham%iprmb,eulal,s_ham%neula,q0,
     .    q0,q0,0,ldha,ldhb,ldha,ldhb,q0,h) ! q0 is used as a dummy
        q0 = 0

C       if (iprint() >= 50) then
C        call yprm('h after spinor rotation',kcplx+2,h,ofhi,ldhax,ldhax,ldhbx)
C       endif
      endif  ! spinor rotation


C --- Orbital rotation ---
      if (mod(opt/1000,10) == 2) goto 999  ! Execute branch only if include orbital part

      allocate(hwk(ldha,2,ldhb,2,2))
      call dcopy(size(h),h,1,hwk,1)

C --- For each spinor component, do ---
      do  is = 1, nspc
      do  js = 1, nspc

        if (kcplx == 0) then
          i1 = is
          i2 = js
          i3 = 1
        else ! kcplx = 2
          i1 = 1
          i2 = is
          i3 = js
        endif

C   --- Rotate h into hwk = R h (rotate row index one site at a time) ---
        do  ib = ib1, ib2

C     ... Set up orbital rotation matrix R = rmat for this g
          if (dasum(3,eulal(ib,1),ib2) == 0) cycle
          call eua2rm(eulal(ib,1),eulal(ib,2),eulal(ib,3),g) ! Rotation matrix from Euler angles
          if (lsph == 0) then
            call ylmrtg(nl2,g,rmat)
C            call yprmi('g, ib=%i',ib,0,1,g,0,3,3,3)
C            call yprmi('rmat, ib=%i',ib,0,1,rmat,0,nl2,nl2,nl2)
          else
            call ylmrtg(nl2,g,wk)
C           prothl and prothr expect rmat in kcplx=0 mode
            call s2sph(0+100*morder,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C            call yprmi('g, ib=%i',ib,0,1,g,0,3,3,3)
C            call yprm('rmat',2,rmat,nl2*nl2,nl2,nl2,nl2)
          endif

C     ... Rotate rows of h into hwk
C         hwk <- R h
          optrot = 0; if (lsph /= 0) optrot = 20
          call prothl(optrot,nl2,nbas,wk,ldmpa,ldima,s_ham%iprmb,[0],g,q0,q0,
     .      rmat,1,ldhb,ib,ib,0,ldhr,[-1],0,hwk(1,i1,1,i2,i3),ofhi,ldhr,ldhbx,
     .      h(1,i1,1,i2,i3),ofhi)
C         mc rmat15 out -sub 1,9,1,nc -x -show out.copt -sub 1,9,1,nc -- -abs -max:g -px
C          call yprmi('H, is=%i js=%i',is,js,kcplx+2,h(1,i1,1,i2,i3),ofhi,ldhax,ldha,ldhb)
C         call yprmi('R H, is,js=%2:1i after ib=%i',[is,js],ib,kcplx+2,hwk(1,i1,1,i2,i3),ofhi,ldhax,ldha,ldhb)
        enddo                   ! Loop over ib

C   --- Rotate hwk = R h(orig) into h = R h(orig) R^-1 (Rotate column index one site at a time) ---
        ldmpb = ldmpa; ldimb = ldima
        do  ib = ib1, ib2

C     ... Set up orbital rotation matrix R = rmat for this g
          if (dasum(3,eulal(ib,1),ib2) == 0) cycle
          call eua2rm(eulal(ib,1),eulal(ib,2),eulal(ib,3),g) ! Rotation matrix from Euler angles
          if (lsph == 0) then
            call ylmrtg(nl2,g,rmat)
C            call yprmi('g, ib=%i',ib,0,1,g,0,3,3,3)
C            call yprmi('rmat, ib=%i',ib,0,1,rmat,0,nl2,nl2,nl2)
          else
            call ylmrtg(nl2,g,wk)
C           prothl and prothr expect rmat in kcplx=0 mode
            call s2sph(0+100*morder,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C           call yprm('rmat',2,rmat,nl2*nl2,nl2,nl2,nl2)
          endif

C     ... Rotate columns of hwk into h
C         h <- R h R+
          call prothr(optrot,nl2,nbas,wk,ldmpb,ldimb,s_ham%iprmb,[0],g,q0,q0,
     .      rmat,1,ldhb,ib,ib,0,ldhr,[-1],0,hwk(1,i1,1,i2,i3),ofhi,ldhr,ldhbx,
     .      h(1,i1,1,i2,i3),ofhi)
C          call yprmi('R H, is=%i js=%i',is,js,kcplx+2,hwk(1,i1,1,i2,i3),ofhi,ldhax,ldha,ldhb)
C          call yprmi('R H R+, is,js=%2:1i after ib=%i',[is,js],ib,kcplx+2,h(1,i1,1,i2,i3),ofhi,ldhax,ldha,ldhb)
        enddo

      enddo ! Loop over spinor index
      enddo ! Loop over spinor index

C      if (iprint() >= 50) then
C        call yprm('h after spinor + orbital rotation',kcplx+2,h,ofhi,ldhax,ldhax,ldhbx)
C      endif

  999 continue
      call tcx('rotheu')

      end
C      subroutine snot(h1,h2,ldha,ofhi)
C      implicit none
C      integer ofhi,ldha
C      double precision h1(ldha,2,2,ldha,2),h2(ldha,2,ldha,2,2)
C
C      print *, h1(2,1,1,2,1),h1(2,1,1,2,2)
C      print *, h1(2,2,1,2,1),h1(2,2,1,2,2)
C
C      print *, h1(2+ldha*4,1,1,1,1),h1(2+ldha*ldha*4+ldha*4,1,1,1,1)
C      print *, h1(2+ldha+ldha*4,1,1,1,1),h1(2+ldha+ldha*ldha*4+ldha*4,1,1,1,1)
C
C      print *, 0,ldha*ldha*4
C      print *, ldha,ldha+ldha*ldha*4
C
CC      print *, h2(2,1,2,1,1),h2(2,1,2,2,1)
CC      print *, h2(2,2,2,1,1),h2(2,2,2,2,1)
CC
CC      print *, h2(2+ldha*2,1,1,1,1),h2(2+ldha*ldha*2+ldha*2,1,1,1,1)
CC      print *, h2(2+ldha+ldha*2,1,1,1,1),h2(2+ldha+ldha*ldha*2+ldha*2,1,1,1,1)
CC
CC      print *, 0,ldha*ldha*2
CC      print *, ldha,ldha+ldha*ldha*2
C
CC      stop
C
C      end