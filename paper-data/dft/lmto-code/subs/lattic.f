      subroutine lattic(s_lat,s_ctrl,s_site)
C- Sets up the real and reciprocal space lattice vectors
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat as tol rpad nkdmx nkqmx gam plat platl platr
Ci                 ldist dist pos
Co     Stored:     vol plat0 plat qlat platl platr awald nkd nkq dist
Co     Allocated:  dlv qlv
Cio    Elts passed:pos
Cio    Passed to:  *
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas npadl npadr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     pos
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Cr Remarks
Cr    For historical reasons, lattice distortions may be EITHER
Cr    defined through gam (special-purpose volume conserving shear) OR
Cr    by one of the ldist modes:
Cr    ldist: 1: defgrd holds rot about spec'd angle
Cr           2, lattice deformed with a general linear transformation
Cr           3, lattice deformed by a shear.
Cu Updates
Cu  17 Jun 13 Replace f77 pointers with f90 ones
Cu  08 May 13 Complete migration to f90 structures; eliminate s_array
Cu  01 Sep 11 Begin migration to f90 structures
Cu   2 Mar 04 Pass rpad to lattc
Cu   5 Jun 01 (ATP) Now calls lattc after lattice transformation
Cu  19 Apr 00 Fixed rotations; new argument list
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      real(8), allocatable :: dlv(:)
      real(8), allocatable :: qlv(:)
C ... Local parameters
      integer ldist,lmxst,nkd,nkdmx,nkq,nkqmx,nbas
      integer i,nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .  xx,xx1,xx2,dotprd,pi,rpad,
     .  platl(3,3),platr(3,3),plat0(3,3),plat(3,3),qlat(3,3),dist(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
C ... External calls
      external daxpy,dcopy,defps2,defrr,lattc,lattdf,pack1,pack5,rdistn,
     .         redfrr,spackv,upack,upack1,upack2

C     call info(30,1,0,' Real and recip space lattices:',0,0)
C     .  alat,plat0)
      alat = s_lat%alat
      awald0 = s_lat%as
      tol = s_lat%tol
      rpad = s_lat%rpad
      nkdmx = s_lat%nkdmx
      nkqmx = s_lat%nkqmx
      gam = s_lat%gam
      alat = s_lat%alat
      plat0 = s_lat%plat
      nbas = s_ctrl%nbas
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nbaspp = nbas + 2*(npadl+npadr)

      if (nbaspp > nbas) then
        call dcopy(9,plat0,1,plat,1)
        platl = s_lat%platl
        platr = s_lat%platr
        call dcopy(6,plat0,1,platl,1)
        call dcopy(6,plat0,1,platr,1)
        call daxpy(3,2d0,platl(1,3),1,plat0(1,3),1)
        call daxpy(3,2d0,platr(1,3),1,plat0(1,3),1)

        pi = 4*datan(1d0)
        xx1 = 180/pi*dotprd(3,platl(1,3),1,plat(1,3),1)
        xx2 = 180/pi*dotprd(3,platr(1,3),1,plat(1,3),1)

        call info8(30,1,0,
     .    ' LATTIC:  Padding Plat(3) with end principal layers: '//
     .    '%N%3;11,6D Plat(3) as input'//
     .    '%N%3;11,6D PlatL:  angle (deg) with Plat(3) = %d'//
     .    '%N%3;11,6D PlatR:  angle (deg) with Plat(3) = %d'//
     .    '%N%3;11,6D Plat(3) including double padding',
     .    plat(1,3),platl(1,3),xx1,platr(1,3),xx2,plat0(1,3),0,0)

      endif

C ... Apply specified linear transformation of lattice and basis vectors
      ldist = s_lat%ldist
      dist = s_lat%dist
C     call prmx('pos from lat',s_lat%pos,3,3,nbaspp)

      if (abs(gt-1) > 1d-10) then
        call rdistn(s_lat%pos,s_lat%pos,nbaspp,gx,gy,gz,gt)
      elseif (ldist /= 0) then
        call lattdf(ldist,dist,plat0,nbaspp,s_lat%pos,0,[0d0],[0d0])
        if (nbaspp > nbas) then
          call lattdf(ldist,dist,platl,0,xx,0,[0d0],[0d0])
          call lattdf(ldist,dist,platr,0,xx,0,[0d0],[0d0])
        endif
      else
        call dpzero(dist,9)
        dist(1,1) = 1
        dist(2,2) = 1
        dist(3,3) = 1
      endif

      allocate(dlv(3*nkdmx),qlv(3*nkqmx))
      lmxst = 6
      if (s_lat%lmxst > 0) lmxst = s_lat%lmxst
      call lattc(awald0,tol,rpad,alat,alat,plat0,gx,gy,gz,gt,plat,qlat,
     .   lmxst,vol,awald,dlv,nkd,qlv,nkq,nkdmx,nkqmx)

      if (nbaspp > nbas) then
        if (cos(pi*xx1/180) < 0) call rx('lattic:  '//
     .    'angle betw/ PlatL and Plat(3) > 90 degrees ... fix PlatL')
        if (cos(pi*xx2/180) < 0) call rx('lattic:  '//
     .    'angle betw/ PlatR and Plat(3) > 90 degrees ... fix PlatR')
      endif

      s_lat%vol = vol
      s_lat%plat0 = plat0
      s_lat%plat = plat
      s_lat%qlat = qlat
      s_lat%platl = platl
      s_lat%platr = platr
      s_lat%awald = awald
      s_lat%nkd = nkd
      s_lat%nkq = nkq
      s_lat%dist = dist

C     call dcopy(3*nbaspp,s_lat%pos,1,s_lat%pos,1)
      s_lat%vol = vol
      s_lat%plat0 = plat0
      s_lat%plat = plat
      s_lat%qlat = qlat
      s_lat%platl = platl
      s_lat%platr = platr
      s_lat%awald = awald
      s_lat%nkd = nkd
      s_lat%nkq = nkq
      s_lat%dist = dist

      call ptr_lat(s_lat,4+1,'dlv',3,nkd,0,dlv)
      call ptr_lat(s_lat,4+1,'qlv',3,nkq,0,qlv)
      deallocate(dlv,qlv)

      do i = 1, nbaspp
        s_site(i)%pos(:) = s_lat%pos(:,i)
      enddo

      end
