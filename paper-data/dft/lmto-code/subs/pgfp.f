      subroutine pgfp(nbas,npadl,npadr,n1,n2,lint,v,vp,iv,ivp)
C- Pads a matrix for basis from bulk layers on LHS and RHS.
C ----------------------------------------------------------------
Ci   nbas  :size of basis
Ci   npadl :number of sites in left padding layer
Ci   npadr :number of sites in right padding layer
Ci   n1    :array v (or iv) is dimensioned v(n1,nbas,n2)
Ci   n2    :array v (or iv) is dimensioned v(n1,nbas,n2)
Ci   lint  :F, matrix is double v; T, matrix is integer iv.
Co vp (or ivp) padded matrix, dimensioned vp(n1,nbasp,n2)
Cr vp(*,nbas+1..nbas+PL(0),*) copied from left PL
Cr vp(*,nbas+PL(0)+1..nbas+PL(1)+PL(*),*) copied from right PL
C ----------------------------------------------------------------
      implicit none
      logical lint
      integer nbas,npadl,npadr,n1,n2
      integer iv(n1,nbas,n2),ivp(n1*(nbas+npadl+npadr)*n2)
      integer nbasp,ib,jb,i2,ioff
      double precision v(n1,nbas,n2),vp(*)

      nbasp = nbas + npadl + npadr
      do  i2 = 1, n2
        do  ib = 1, nbasp
          jb = ib
          if (jb > nbas+npadl) then
            jb = jb - npadl - npadr
          elseif (jb > nbas) then
            jb = jb - nbas
          endif
          ioff = nbasp*n1*(i2-1) + n1*(ib-1) + 1
          if (lint) then
            call icopy(n1,iv(1,jb,i2),1,ivp(ioff),1)
          else
            call dcopy(n1, v(1,jb,i2),1, vp(ioff),1)
          endif
        enddo
      enddo

c      call snit(nbas,nbasp,lint,n1,n2,v,vp,iv,ivp)
      end
      subroutine pgfpp(nbas,npadl,npadr,n1,n2,lint,lalloc,v,vpp)
C- Doubly pads a matrix for basis from bulk layers on LHS and RHS.
C ----------------------------------------------------------------------
Ci   nbas  :size of basis
Ci   npadl :number of sites in left padding layer
Ci   npadr :number of sites in right padding layer
Ci   n1    :array v is dimensioned v(n1,nbas,n2)
Ci   n2    :array v is dimensioned v(n1,nbas,n2)
Ci   lint  :F, matrix is double v; T, matrix is integer iv.
Ci   lalloc:true, if pgfpp is to allocate destination array vpp
Ci   v     :input matrix to pad
Co Outputs
Co   vpp   :padded matrix
Cu Updates
Cu   05 Jul 13 partially replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lint,lalloc
      integer nbas,npadl,npadr,n1,n2,v(*),vpp(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:)
C ... Local parameters
      integer nbasp,nbaspp
C heap:
      integer w(1)
      common /w/ w

      nbasp = nbas + npadl + npadr
      if (nbasp == nbas) return
      nbaspp = (2*nbasp-nbas)
      if (lalloc .and. lint) call defi(vpp(1),  n1*nbaspp*n2)
      if (lalloc .and. .not. lint) call defdr(vpp(1),  n1*nbaspp*n2)
      allocate(wk(n1*nbasp*n2))
      call pgfp(nbas,npadl,npadr,n1,n2,lint,v,wk,v,wk)
      if (lalloc) then
        call pgfp(nbasp,npadl,npadr,n1,n2,lint,wk,
     .    w(vpp(1)),wk,w(vpp(1)))
      else
        call pgfp(nbasp,npadl,npadr,n1,n2,lint,wk,vpp,wk,vpp)
      endif
      deallocate(wk)
      end
C      subroutine pgfpad(nbas,nbasp,npadl,npadr,nl,nkap,plat,platl,platr,
C     .  mode,iclass,bas,pgfsl,alpha,adot,oiclpp,obaspp,opgslp,oalpp,
C     .  oadotp)
CC- Create doubly padded alpha,adot,iclasss,pgfsl,bas for layer geometry
CC     mode  1s digit 2  doubly pad iclass
CC          10s digit 2  doubly pad bas
CC         100s digit 2  doubly pad pgfsl
CC        1000s digit 2  doubly pad alpha, adot
C      implicit none
C      integer nbas,nbasp,npadr,iclass(1),pgfsl(1),oalpp,oadotp,oiclpp,
C     .  opgslp,obaspp
C      double precision alpha(*),adot(*),plat(9),platl(9),platr(9),
C     .  bas(3,nbas)
CC heap:
C      integer w(1)
CC Local variable
C      integer nl,nkap,npadl,nbaspp,mode,obasp
C      logical T,F
C      parameter (T=.true., F=.false.)
C      common /w/ w
C
CC --- iclass ---
C      if (mod(mode,10) == 2)
C     .  call pgfpp(nbas,npadl,npadr,1,1,T,T,iclass,oiclpp)
CC --- Doubly pad bas.  Looks funny, but it's ok ---
C      if (mod(mode/10,10) == 2)  then
C        nbaspp = 2*nbasp - nbas
C        call defdr(obaspp, 3*nbaspp)
C        call defdr(obasp, 3*nbasp)
C        call pgbasp(nbas,npadl,npadr,bas,plat,platl,platr,
C     .    w(obasp))
C        call dscal(9,2d0,platl,1)
C        call pgbasp(nbasp,npadl,npadr,w(obasp),plat,platl,platr,
C     .    w(obaspp))
C        call dscal(9,.5d0,platl,1)
C        call rlse(obasp)
C      endif
CC --- pgfsl ---
C      if (mod(mode/100,10) == 2)  then
C        call pgfpp(nbas,npadl,npadr,1,1,T,T,pgfsl,opgslp)
C        call decc(w(opgslp+nbas),npadl)
C        call incc(w(opgslp+nbas+npadl),npadl)
C        call decc(w(opgslp+nbasp),npadl)
C        call decc(w(opgslp+nbasp),npadl)
C        call incc(w(opgslp+nbasp+npadl),npadl)
C        call incc(w(opgslp+nbasp+npadl),npadl)
C      endif
CC --- alpha, adot ---
C      if (mod(mode/1000,10) == 2)  then
C        call pgfpp(nbas,npadl,npadr,nl**2,nkap**2,F,T,alpha,oalpp)
C        call pgfpp(nbas,npadl,npadr,nl**2,nkap**2,F,T,adot,oadotp)
C      endif
C      end


C  for testing ...
C      subroutine snit(nbas,nbasp,lint,n1,n2,v,vp,iv,ivp)
C      implicit none
C      logical lint
C      integer nbas,nbasp,n1,n2
C      integer iv(n1,nbas,n2),ivp(n1,nbasp,n2)
C      integer ib,jb,i2,i
C      double precision v(n1,nbas,n2),vp(n1,nbasp,n2)
C
CC      do  10  i2 = 1, n2
CC        do  12  ib = 1, nbasp
CC          jb = ib
CC          if (jb > nbas+npadl) then
CC            jb = jb - npadl - npadr
CC          elseif (jb > nbas) then
CC            jb = jb - nbas
CC          endif
CC          if (lint) then
CC            call icopy(n1,iv(1,jb,i2),1,ivp(1,ib,i2),1)
CC          else
CC            call dcopy(n1,v(1,jb,i2),1,vp(1,ib,i2),1)
CC          endif
CC   12   continue
CC   10 continue
C      do  16  i2 = 1, n2
C      do  16  ib = 1, nbasp
C   16 print 333, ib,i2, (vp(i,ib,i2), i= 1,n1)
CC     .    (vp(i,ib,i2) - v(i,ib,i2), i= 1,n1)
C  333 format(2i4,10f12.6)
CC 333 format(2i4,10i12)
C      pause
C      end
