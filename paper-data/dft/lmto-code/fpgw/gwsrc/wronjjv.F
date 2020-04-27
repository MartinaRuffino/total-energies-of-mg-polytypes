      subroutine wronjje(job,e1,e2,r,n1,n2,lmax,wjj)
C- Wronskians of a pair of spherical Bessel functions for a vector of energies
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit (passed to ropbes as 100s digit)
Ci         :0 always use besslr
Ci         :1 call besnu for y>=90
Ci         :2 always call dbesnu
Ci   e1    :e1 fpr Wronskian W{j(e1),j(e2)}
Ci   e2    :e2 for Wronskian W{j(e1),j(e2)}
Ci   r     :radius
Ci   n1    :number of energies e1
Ci   n2    :number of energies e1
Ci   lmax  :Wronskians for l=0:lmax
Co Outputs
Co   wjj   :Wronskian W{j(e1),j(e2)}
Cs Command-line switches
Cr Remarks
Cr  wjj is continuous as e2->e1, and is symmetric in e1,e2.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,lmax,job
      double precision e1(n1),e2(n2),r
      real(8) :: wjj(n1,n2,0:lmax)
C ... Local parameters
      integer i1,i2,l
      real(8) :: r2,r3,rk,rj,efac
      real(8) :: aj1(n1,0:lmax),aj2(n2,0:lmax),dj1(n1,0:lmax),dj2(n2,0:lmax)
      real(8) :: ajd(0:lmax),djd(0:lmax)
      real(8),parameter :: tol = 1d-8

      if (r == 0) then
        call dpzero(wjj,size(wjj))
        return
      endif

      call radjjv(e1,[r],n1,lmax,aj1,dj1,100*job+20)
      call radjjv(e2,[r],n2,lmax,aj2,dj2,100*job+20)
      r2 = r*r

C --- Loop over (e1,e2) pairs ---
      do  i1 = 1, n1
      do  i2 = 1, n2
C   ... Case e1=e2=0
        if (dabs(e1(i1)) <= tol .and. dabs(e2(i2)) <= tol) then
          r3 = r**3
          rk = -1d0
          rj = 1d0/r
          do  l = 0, lmax
            rk = rk*(2*l-1)/r
            rj = rj*r/(2*l+1)
            wjj(i1,i2,l) = -rj*rj*r3/(2*l+3)
          enddo
          cycle

C   ... Case e1 /= e2
        elseif (dabs(e1(i1)-e2(i2)) > tol) then
          efac = 1d0/(e2(i2)-e1(i1))
          do  l = 0, lmax
            wjj(i1,i2,l) = efac*r2*(aj1(i1,l)*dj2(i2,l)-dj1(i1,l)*aj2(i2,l))
          enddo

C   ... Case e1 == e2 but not zero
        else
          call radjjv(e1(i1),[r],1,lmax,ajd,djd,100*job+11)
          do  l = 0, lmax
            wjj(i1,i2,l) = r2*(aj1(i1,l)*djd(l)-dj1(i1,l)*ajd(l))
          enddo
        endif
      enddo
      enddo

      end

      subroutine wronjjr(job,e1,e2,r,n,lmax,wjj)
C- Wronskians of a pair of spherical Bessel functions for a vector of radii
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit (passed to ropbes as 100s digit)
Ci         :0 always use besslr
Ci         :1 call besnu for y>=90
Ci         :2 always call dbesnu
Ci   e1    :First energy for Wronskian
Ci   e2    :Second energy for Wronskian
Ci   r     :vector of radii
Ci   n     :Number of points
Ci   lmax  :maximum l
Co Outputs
Co   wjj   :W(j1,j2) (Wronskian) for l=0..lmax
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  wjj is continuous as e2->e1, and is symmetric in e1,e2.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,lmax,job
      double precision e1,e2,r(n)
      real(8) :: wjj(n,0:lmax)
C ... Local parameters
      integer i,l
      real(8) :: r3,rk,rj,efac
      real(8) :: aj1(n,0:lmax),aj2(n,0:lmax),dj1(n,0:lmax),dj2(n,0:lmax)
      real(8),parameter :: tol = 1d-8

C --- Case e1=e2=0 ---
      if (dabs(e1) <= tol .and. dabs(e2) <= tol) then
        do  i = 1, n
          if (r(i) == 0) then
            wjj(i,:) = 0
            cycle
          endif
          r3 = r(i)**3
          rk = -1d0
          rj = 1d0/r(i)
          do  l = 0, lmax
            rk = rk*(2*l-1)/r(i)
            rj = rj*r(i)/(2*l+1)
            wjj(i,l) = -rj*rj*r3/(2*l+3)
          enddo
        enddo

C --- Case e1 /= e2 ---
      elseif (dabs(e1-e2) > tol) then
        efac = 1d0/(e2-e1)
        call radjjv([e1],r,n,lmax,aj1,dj1,100*job+10)
        call radjjv([e2],r,n,lmax,aj2,dj2,100*job+10)
        do  l = 0, lmax
        do  i = 1, n
          wjj(i,l) = efac*r(i)*r(i)*(aj1(i,l)*dj2(i,l)-dj1(i,l)*aj2(i,l))
        enddo
        enddo

C --- Case e1 == e2 but not zero ---
      else
        call radjjv([e1],r,n,lmax,aj1,dj1,100*job+10)
        call radjjv([e1],r,n,lmax,aj2,dj2,100*job+11)
        do  l = 0, lmax
        do  i = 1, n
          wjj(i,l) = r(i)*r(i)*(aj1(i,l)*dj2(i,l)-dj1(i,l)*aj2(i,l))
        enddo
        enddo
      endif
      end

      subroutine radjjv(e,r,n,lmax,aj,dj,job)
C- Spherical Bessel functions and derivatives on a mesh of points
C ----------------------------------------------------------------
Ci Inputs
Ci   e    :Make j_l(e,r)
Ci   r    :list of points
Ci   n    :number of radial points
Ci   lmax :Make j_l for l=0..lmax
Ci   job  :1s digit
Ci        :0, makes values and slopes
Ci        :1, makes energy derivatives
Ci        :10's digit (passed to ropbes as 1s digit)
Ci        :0, r and e are both vectors
Ci        :1, r is a vector and e is a scalar
Ci        :2, r is a scalar and e is a vector
Ci        :100s digit (passed to ropbes as 10s digit)
Ci        :0 always use besslr
Ci        :1 call besnu for y>=90
Ci        :2 always call dbesnu
Co Outputs
Co   If 1s digit of job is 0:
Co   aj   :Bessel function j_l(r,e,l)
Co   dj   :radial derivative of aj
Co   If 1s digit of job is 1 aj,dj are the energy derivatives of the above.
Cr Remarks
Cr   This is a vectorized version of radkj, returning only the Bessel part.
Cr   j(r,lmax,lmax-1)  are calculated by a power series.
Cr   j for lower l are calculated by downward recursion.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,n,job
      double precision e(n),r(n),aj(n,0:lmax),dj(n,0:lmax)
C ... Local parameters
      integer i,l,job1,job12
      double precision h(n),er2(n),phi(n,0:lmax+2),rl,r2,php(0:lmax+1),ri

      job1 = mod(job/10,10)
      job12 = mod(job/10,100)

      if (mod(job,10) == 0) then
        call ropbes(r,e,lmax+1,er2,h,phi,n,job12)
        do  i = 1, n
          ri = r(i); if (job1 == 2) ri=r(1)
          if (ri == 0) then
            aj(i,0) = 1; aj(i,1:) = 0; dj(i,:) = 0
            cycle
          endif
          rl = 1
          do  l = 0, lmax
            aj(i,l) = phi(i,l)*rl
            dj(i,l) = (l*phi(i,l) - er2(i)*phi(i,l+1))*rl/ri
            rl = rl*ri
          enddo
        enddo
      else
        call ropbes(r,e,lmax+2,er2,h,phi,n,job12)
        do  i = 1, n
          ri = r(i); if (job1 == 2) ri=r(1)
          if (ri == 0) then
            aj(i,:) = 0; dj(i,:) = 0
            cycle
          endif
          r2 = ri**2
          do  l = 0, lmax+1
            php(l) = -0.5d0*r2*phi(i,l+1)
          enddo
          rl = 1
          do  l = 0, lmax
            aj(i,l) = php(l)*rl
            dj(i,l) = (l*php(l) - er2(i)*php(l+1) - r2*phi(i,l+1))*rl/ri
            rl = rl*ri
          enddo
        enddo
      endif

      end
C     testing ...
C      subroutine fmain
C      implicit none
C      integer nr1,nr,lmax
C      parameter (nr1=2,nr=8,lmax=5)
C      double precision ri(nr),e(nr),aj(nr,0:lmax),dj(nr,0:lmax),wjj2(nr1,nr,0:lmax),wjj(nr,0:lmax)
CC     double precision phi(0:lmax+1),psi(0:lmax+1)
C      double precision aki(0:lmax+1),dki(0:lmax+1),aji(0:lmax+1),dji(0:lmax+1),errmxa,errmxd
C      double precision fkk(0:lmax),fkj(0:lmax),fjk(0:lmax),fjj(0:lmax)
C      integer i,l,i1,j
C
C
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
C  333 format(' i',i3,10f12.6)
C  334 format(5x,12x,10f12.3)
C  335 format(5x,24x,10f12.3)
C
C      call info2(2,0,0,' testing radjjv job 10, e=%d ...',e(nr),0)
C      print *, '          r           l=0 ...'
C      call radjjv(e(nr),ri,nr,lmax,aj,dj,10)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, ri(i),(aj(i,l), l=0,lmax)
C        if (ri(i) == 0) cycle
C        call radkj(e(nr),ri(i),lmax,aki,aji,dki,dji,0)
C        print 333, i, ri(i),(aji(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((aj(i,l)-aji(l))/aji(l)))
C          errmxd = max(errmxd,abs((dj(i,l)-dji(l))/dji(l)))
C        enddo
C        print 334, ((aj(i,l)-aji(l))*1e15/aji(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g   in derivative = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing radjjv job 11, e=%d ...',e(nr),0)
C      call radjjv(e(nr),ri,nr,lmax,aj,dj,11)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, ri(i),(aj(i,l), l=0,lmax)
C        if (ri(i) == 0) cycle
C        call radkj(e(nr),ri(i),lmax,aki,aji,dki,dji,1)
C        print 333, i, ri(i),(aji(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((aj(i,l)-aji(l))/aji(l)))
C          errmxd = max(errmxd,abs((dj(i,l)-dji(l))/dji(l)))
C        enddo
C        print 334, ((aj(i,l)-aji(l))*1e15/aji(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g   in derivative = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing radjjv job 20, r=%d ...',ri(nr),0)
C      print *, '          e           l=0 ...'
C      call radjjv(e,ri(nr),nr,lmax,aj,dj,20)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, e(i),(aj(i,l), l=0,lmax)
C        if (e(i)*ri(i) == 0) cycle
C        call radkj(e(i),ri(nr),lmax,aki,aji,dki,dji,0)
C        print 333, i, e(i),(aji(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((aj(i,l)-aji(l))/aji(l)))
C          errmxd = max(errmxd,abs((dj(i,l)-dji(l))/dji(l)))
C        enddo
C        print 334, ((aj(i,l)-aji(l))*1e15/aji(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g   in derivative = %g',errmxa,errmxd)
C
C      call info2(2,0,0,' testing radjjv job 21, r=%d ...',ri(nr),0)
C      print *, '          e           l=0 ...'
C      call radjjv(e,ri(nr),nr,lmax,aj,dj,21)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, e(i),(aj(i,l), l=0,lmax)
C        if (e(i)*ri(i) == 0) cycle
C        call radkj(e(i),ri(nr),lmax,aki,aji,dki,dji,1)
C        print 333, i, e(i),(aji(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((aj(i,l)-aji(l))/aji(l)))
C          errmxd = max(errmxd,abs((dj(i,l)-dji(l))/dji(l)))
C        enddo
C        print 334, ((aj(i,l)-aji(l))*1e15/aji(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Bessel = %g   in derivative = %g',errmxa,errmxd)
C
C      print *, 'testing wronjjr e1=e2=0 ...'
C      print *, '          r           l=0 ...'
C      call wronjjr(0,0d0,0d0,ri,nr,lmax,wjj)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, ri(i),(wjj(i,l), l=0,lmax)
C        if (ri(i) == 0) cycle
C        call wronkj(jobBes,0d0,0d0,ri(i),lmax,fkk,fkj,fjk,fjj)
C        print 333, i, ri(i),(fjj(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((wjj(i,l)-fjj(l))/fjj(l)))
C        enddo
C        print 334, ((wjj(i,l)-fjj(l))*1e15/fjj(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Wronskian = %g',errmxa,0)
C
C      call info2(2,0,0,' testing wronjjr e1=e2 = %d ...',e(nr),0)
C      print *, '          r           l=0 ...'
C      call wronjjr(0,e(nr),e(nr),ri,nr,lmax,wjj)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, ri(i),(wjj(i,l), l=0,lmax)
C        if (ri(i) == 0) cycle
C        call wronkj(jobBes,e(nr),e(nr),ri(i),lmax,fkk,fkj,fjk,fjj)
C        print 333, i, ri(i),(fjj(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((wjj(i,l)-fjj(l))/fjj(l)))
C        enddo
C        print 334, ((wjj(i,l)-fjj(l))*1e15/fjj(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Wronskian = %g',errmxa,0)
C
C      call info2(2,0,0,' testing wronjjr e1=%d e2=%d ...',e(nr-2),e(nr))
C      print *, '          r           l=0 ...'
C      call wronjjr(0,e(nr-2),e(nr),ri,nr,lmax,wjj)
C      errmxa = 0; errmxd = 0
C      do  i = 1, nr
C        print 333, i, ri(i),(wjj(i,l), l=0,lmax)
C        if (ri(i) == 0) cycle
C        call wronkj(jobBes,e(nr-2),e(nr),ri(i),lmax,fkk,fkj,fjk,fjj)
C        print 333, i, ri(i),(fjj(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((wjj(i,l)-fjj(l))/fjj(l)))
C        enddo
C        print 334, ((wjj(i,l)-fjj(l))*1e15/fjj(l), l=0,lmax)
C      enddo
C      call info2(2,0,1,' maximum deviation in Wronskian = %g',errmxa,0)
C
C      call info2(2,0,0,' testing wronjje r=%d...',ri(nr),0)
C      print *, '          e1          e2          l=0 ...'
C      call wronjje(0,e(nr-1),e,ri(nr),nr1,nr,lmax,wjj2)
C      errmxa = 0; errmxd = 0
C      j = 0
C      do i1 = nr-1, nr
C      j = j+1
C      do  i = 1, nr
C        print 333, i, e(i1),e(i),(wjj2(j,i,l), l=0,lmax)
CC       if (e(i) == 0) cycle
C        call wronkj(jobBes,e(i1),e(i),ri(nr),lmax,fkk,fkj,fjk,fjj)
C        print 333, i, e(i1),e(i),(fjj(l), l=0,lmax)
C        do  l  = 0, lmax
C          errmxa = max(errmxa,abs((wjj2(j,i,l)-fjj(l))/fjj(l)))
C        enddo
C        print 335, ((wjj2(j,i,l)-fjj(l))*1e15/fjj(l), l=0,lmax)
C      enddo
C      enddo
C      call info2(2,0,1,' maximum deviation in Wronskian = %g',errmxa,0)
C
C
C      end
