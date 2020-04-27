      subroutine ropbes(r,e,lmax,y,h,xi,n,job)
C- Radial Bessel functions divided by r**l (vectorizes)
C ----------------------------------------------------------------
Ci Inputs
Ci   r    :vector of points (but see job)
Ci   e    :vector of energies (but see job)
Ci   h    : work array of length n
Ci   job: :1s digit
Ci        :0, r and e are both vectors
Ci        :1, r is a vector and e is a scalar
Ci        :2, r is a scalar and e is a vector
Ci        :10s digit
Ci        :Remake bessel for any y>100, using dbesnu
Co Outputs
Co   y    y(i) = e*r(i)**2
Co   xi   J(r,l)/r**l, according to standard definition
Cr Remarks
Cr   J(r,lmax,lmax-1)  are calculated by a power series.
Cr   J for lower l are calculated by downward recursion
Cu Updates
Cu   21 Aug 17  New 10s digit job
Cu   13 May 15  Changed tolerance to 1d-15
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,n,job
      double precision e(n),r(n),xi(n,0:lmax),h(n),y(n)
C ... Local parameters
      integer i,l,job0
      double precision xx,phi(0:lmax)

      if (lmax < 0) return
      job0 = mod(job,10)
      if (job0 == 0) then
        do  i = 1, n
          y(i) = r(i)**2*e(i)
        enddo
      elseif (job0 == 1) then
        do  i = 1, n
          y(i) = r(i)*r(i)*e(1)
        enddo
      elseif (job0 == 2) then
        xx = r(1)**2
        do  i = 1, n
          y(i) = xx*e(i)
        enddo
      else
        call rx('ropbes: bad job0')
      endif

C --- Power series expansion for lmax, lmax-1 ---
      call ropbs1(y,lmax,xi(1,lmax),h,n)
      if (lmax >= 1) call ropbs1(y,lmax-1,xi(1,lmax-1),h,n)

C --- Downward recursion ---
      do  l = lmax-2, 0, -1
        xx = 2*l+3
        do  i = 1, n
          xi(i,l) = xx*xi(i,l+1)-y(i)*xi(i,l+2)
        enddo
      enddo

C --- Replace result with dbesnu for large argument ---
      if (job < 10) return

      do  i = 1, n
        if (y(i) > 90) then
          call bessjs(y(i),lmax,11,phi,[xx])
          xi(i,0:lmax) = phi(0:lmax)
        endif
      enddo

      end
      subroutine ropbs1(y,l,xi,h,n)
C- Evaluates bessel function for one l using power series expansion
      implicit none
      integer n,l
      double precision xi(n),h(n),y(n)
      integer i,k
      double precision df,tol,top,xx

      tol = 1d-12
      df = 1d0
      do  k = 3, 2*l+1,2
        df = df*k
      enddo
      do  i = 1, n
        xi(i) = 1d0
        h(i) = 1d0
      enddo
      do  k = 1, 1000
        xx = -1d0/( 2*k*(2*k+2*l+1) )
        do  i = 1, n
          h(i) = xx*h(i)*y(i)
          xi(i) = xi(i)+h(i)
        enddo
        top = 0d0
        do  i = 1, n
          top = dmax1(top,dabs(h(i)))
        enddo
        if (top <= tol) goto 10
      enddo
      call rx('ropbes: power series failed to converge')
   10 continue
      xx = 1d0/df
      do  i = 1, n
        xi(i) = xi(i)*xx
      enddo
      end
C     testing ...
C      subroutine fmain
C      implicit none
C      integer nr,lmax
C      parameter (nr=8,lmax=6)
CC     parameter (nr=8*50000,lmax=6)
C      double precision xi(nr,0:lmax),ri(nr),e(nr),y(nr),h(nr),pi
C      double precision phi(0:lmax+1),psi(0:lmax+1),errmxa,errmxd
C      integer i,l
CC     integer ierr
CC     double precision x,xx,phi2(0:lmax+1),swalltime,delwc,e2
C
C      pi = 4*datan(1d0)
C      ri(1) = 0
C      ri(2) = 0.001d0
C      do  10  i = 3, nr
C   10 ri(i) = 0.5*dble(i)**2/nr-.5d0
C      e(1) = 0
C      e(2) = 0.001d0
C      do  i = 3, nr
C        e(i) = (-1)**i*dble(i)/4
C      enddo
C
C  333 format(' i',i3,12f13.6)
C  334 format(5x,12x,12f13.3)
C
CC      e2 = 2d0  ! dbesnu handles only positive arguments
CC      call dscal(nr,3/ri(nr),ri,1)
CC      call info2(2,0,0,' timing study ropbes job 1 ...  e=%d',e2,0)
CC      swalltime = delwc()
CC      call ropbes(ri,[e2],lmax,y,h,xi,nr,1)
CC      swalltime = delwc()
CC      print 128, swalltime
CC  128 format(' finished ropbes call ... elapsed wall time',F9.4)
CC      do  i = 1, nr
CC        x = sqrt(e2*ri(i)**2)
CC        call dbesnu(x, lmax+0.5d0, 11, phi2, xx, ierr)  ! Only for positive x
CC        do  l = 0, lmax
CC          phi2(l) = phi2(l)*sqrt(pi/2/x)/x**l
CC        enddo
CC      enddo
CC      swalltime = delwc()
CC      print 129, swalltime
CC  129 format(' finished dbesnu call ... elapsed wall time',F9.4)
CC      do  i = nr, nr
CC        print 333, 0, ri(i),(xi(i,l), l=0,lmax)
CC        print 333, 0, ri(i),(phi2(l), l=0,lmax)
CC        print 334, ((xi(i,l)-phi2(l))*1e15/phi2(l), l=0,lmax)
CC      enddo
CC      stop
C
C      call info2(2,0,0,' testing ropbes job 0...',0,0)
C      call ropbes(ri,e,lmax,y,h,xi,nr,0)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        call bessl(e(i)*ri(i)**2,lmax,phi,psi)
C        print 333, i, e(i),(xi(i,l), l=0,lmax)
C        print 333, i, e(i),(phi(l), l=0,lmax)
C        print 334, ((xi(i,l)-phi(l))*1e15/phi(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((xi(i,l)-phi(l))/phi(l)))
C        enddo
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing ropbes job 1, e=%d ...',e(nr),0)
C      call ropbes(ri,e(nr),lmax,y,h,xi,nr,1)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        call bessl(e(nr)*ri(i)**2,lmax,phi,psi)
C        print 333, i, ri(i),(xi(i,l), l=0,lmax)
C        print 333, i, ri(i),(phi(l), l=0,lmax)
C        print 334, ((xi(i,l)-phi(l))*1e15/phi(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((xi(i,l)-phi(l))/phi(l)))
C        enddo
CC       Check against dbesnu
CC        if (e(nr) > 0) then
CC          x = sqrt(e(nr)*ri(i)**2)
CC          call dbesnu(x, lmax+0.5d0, 11, phi2, xx, ierr)
CC          do  l = 0, lmax
CC            phi2(l) = phi2(l)*sqrt(pi/2/x)/x**l
CC          enddo
CC          print 333, i, ri(i),(phi2(l), l=0,lmax)
CC          print 334, ((xi(i,l)-phi2(l))*1e15/phi2(l), l=0,lmax)
CC        endif
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing ropbes job 2, r=%d ...',ri(nr),0)
C      call ropbes(ri(nr),e,lmax,y,h,xi,nr,2)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        call bessl(e(i)*ri(nr)**2,lmax,phi,psi)
C        print 333, i, e(i),(xi(i,l), l=0,lmax)
C        print 333, i, e(i),(phi(l), l=0,lmax)
C        print 334, ((xi(i,l)-phi(l))*1e15/phi(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((xi(i,l)-phi(l))/phi(l)))
C        enddo
C      enddo
C
CC      call info2(2,0,0,' testing ropbes job 0, large energies, no check for large y...',0,0)
CC      e = [4d0, 6d0, 8d0, 10d0, 12d0, 15d0, 16d0, 20d0]
CC      call ropbes(ri,e,lmax,y,h,xi,nr,0)
CC      errmxa = 0; errmxd = 0
CC      do  i = 1, nr
CC        call bessl(e(i)*ri(i)**2,lmax,phi,psi)
CC
CC        print 333, i, e(i),(xi(i,l), l=0,lmax)
CC        print 333, i, e(i),(phi(l), l=0,lmax)
CC        print 334, ((xi(i,l)-phi(l))*1e15/phi(l), l=0,lmax)
CC        do  l  = 0, lmax
CC          errmxa = max(errmxa,abs((xi(i,l)-phi(l))/phi(l)))
CC        enddo
CC      enddo
CC      call info2(2,0,1,' maximum deviation in Bessel = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing ropbes job 0, large energies, check for large y...',0,0)
C      e = [4d0, 6d0, 8d0, 10d0, 12d0, 15d0, 16d0, 20d0]
C      call ropbes(ri,e,lmax,y,h,xi,nr,10)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        call bessjs(e(i)*ri(i)**2,lmax,11,phi,psi)
C
C        print 333, i, e(i),(xi(i,l), l=0,lmax)
C        print 333, i, e(i)*ri(i)**2,(phi(l), l=0,lmax)
C        print 334, ((xi(i,l)-phi(l))*1e15/phi(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((xi(i,l)-phi(l))/phi(l)))
C        enddo
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g',errmxa,errmxd)
C
C      end
