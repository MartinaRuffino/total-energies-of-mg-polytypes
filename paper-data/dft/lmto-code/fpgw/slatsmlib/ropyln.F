      subroutine ropyln(n,x,y,z,lmax,nd,yl,rsq)
C- Normalized spheric harmonic polynomials (vectorizes).
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of points for which to calculate yl
Ci   x     :x component of Cartesian coordinate
Ci   y     :y component of Cartesian coordinate
Ci   z     :z component of Cartesian coordinate
Ci   lmax  :compute yl for l=0..lmax
Ci   nd    :leading dimension of yl; nd must be >= n
Co Outputs
Co   yl    :Ylm(i,ilm) are the (real) spherical harmonics
Co         :for l=0..lmax and point i.
Co   rsq   :rsq(i) square of length of point i
Cr Remarks
Cr   yl = real harmonics (see Takao's GW note) * r^l
Cr   Explicit forms for s,p:
Cr     l  m  L     Y_L
Cr     0  0  1   1/sqrt(4*pi)
Cr     1 -1  2   sqrt(3/4/pi) y
Cr     1 -1  3   sqrt(3/4/pi) z
Cr     1 -1  4   sqrt(3/4/pi) x
Cu Updates
Cu  24 Apr 09 (Lozovoi) replace w() with allocate
Cu  25 Jun 03 (Kino) initialize cx to zero
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nd,n,m,lmax
      double precision x(n),y(n),z(n),yl(nd,(lmax+1)**2),rsq(n)
C ... Local parameters
      integer l,kk
      double precision fpi,f2m,cx(3)
      double precision, allocatable:: cm(:),sm(:),q(:),h(:)

C     call tcn('ropyln')
      if (n > nd) call rx('ropyln: nd too small')
      fpi = 16*datan(1d0)
      allocate (cm(n),sm(n),q(n*2),h(n))
c     call defrr(ocm,  n)
c     call defrr(osm,  n)
c     call defrr(oq,   n*2)
c     call defrr(oh,   n)
      rsq(1:n) = x(1:n)*x(1:n)+y(1:n)*y(1:n)+z(1:n)*z(1:n)
      cx(:) = 0

C --- Loop over m: cx are normalizations for l, l-1, l-2 ---
      f2m = 1d0
      do  m = 0, lmax
        call ropcsm(m,n,x,y,h,cm,sm)
        if (m == 0) then
          cx(1) = dsqrt(1/fpi)
        else
          f2m = f2m*(2*m*(2*m-1))
          cx(1) = dsqrt(((2*m+1)*2)/(fpi*f2m))
        endif
        do l = m, lmax
          call ropqln(m,l,n,rsq,z,cx,q,kk)
          call ropynx(m,l,kk,n,q,cm,sm,nd,yl)
          cx(3) = cx(2)
          cx(2) = cx(1)
          cx(1) = cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
        enddo
      enddo
      deallocate(cm,sm,q,h)
C     call tcx('ropyln')
      end
C     Separate subroutines below to avoid problems with some
C     optimizing compilers.
C      subroutine ropqln(m,l,n,r2,z,cx,q,kk)
CC- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
Cc  Returns kk, which points to the current component of q.
C      implicit none
C      integer mm,n,i,l,m,kk,k2,k1
C      double precision q(n,2),r2(n),z(n),cx(3)
C      double precision a,b,xx,yy
C
CC --- Case l=m ---
C      if (l == m) then
C        a = 1d0
C        do  1  mm = 0, m-1
C    1   a = a*(2*mm+1)
C        kk = 1
C        a = a*cx(1)
C        do  2  i = 1, n
C    2   q(i,kk) = a
C        return
C      endif
C
CC --- Case l=m+1 ---
C      if (l == m+1) then
C        b = 1d0
C        do  3  mm = 0, m
C    3   b = b*(2*mm+1)
C        b = b*cx(1)
C        kk = 2
C        do  4  i = 1, n
C    4   q(i,kk) = b*z(i)
C        return
C      endif
C
CC --- Case l=m+2 and higher by recursion ---
C      if (l >= m+2) then
C        k2 = kk
C        k1 = kk+1
C        if (k1 == 3) k1 = 1
C        xx = -(l+m-1d0)/(l-m)*cx(1)/cx(3)
C        yy = (2*l-1d0)/(l-m)*cx(1)/cx(2)
C        do  6  i = 1, n
C    6   q(i,k1) = xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
C        kk = k1
C        return
C      endif
C      end
C      subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
C      implicit none
C      integer lav,n,nd,l,i,m,kk
C      double precision q(n,2),cm(n),sm(n),yl(nd,1)
C      lav = l*(l+1)+1
C      do  1  i = 1, n
C    1 yl(i,lav+m) = cm(i)*q(i,kk)
C      if (m == 0) return
C      do  2  i = 1, n
C    2 yl(i,lav-m) = sm(i)*q(i,kk)
C      end
C#ifdefC TEST
CC Test program to check ropyln
C      subroutine fmain
C      implicit none
C      integer nrx,nlmx,nr,lmax,nlm1,ir,ii,l,ilm,i1,i2,nsize
C      parameter (nrx=20,nlmx=49,nsize=100000)
C      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),
C     .  ylv(nrx,nlmx),ylok(nrx,nlmx),dr(3),tops,ylm(nlmx)
C
C      call sylmnc(cy,15)
C
C      lmax = 2
C   99 print *, 'lmax='
C      read(*,*) lmax
C
C      call makr(0d0,nr,x,y,z)
C
CC ... Make nonvectorized ylm's
C      nlm1 = (lmax+1)**2
C      do  ir = 1, nr
C        dr(1) = x(ir)
C        dr(2) = y(ir)
C        dr(3) = z(ir)
C        call sylm(dr,ylm,lmax,r2)
C        do  ilm = 1, nlm1
C          ylok(ir,ilm) = cy(ilm)*ylm(ilm)
C        enddo
C      enddo
CC     Test: Y_1-1 = sqrt(3/4/pi) y
CC     print *, y(1) * dsqrt(0.75/4/atan(1d0))
CC     print *, ylok(1,2)
C
C      call ropyln(nr,x,y,z,lmax,nrx,ylv,r2)
C
C      tops = 0
C      do  10  ir = 1, nr
C        do  12  l = 0, lmax
C        i1 = l*l+1
C        i2 = (l+1)**2
C   12   print 333, (ylok(ir,ii),ii=i1,i2)
C  333   format(9f8.5)
C        print *
C        do  14  l = 0, lmax
C        i1 = l*l+1
C        i2 = (l+1)**2
C        do  16  ii = i1, i2
C   16   tops = max(tops,dabs(ylok(ir,ii)-ylv(ir,ii)))
C   14   print 333, (ylok(ir,ii)-ylv(ir,ii),ii=i1,i2)
C        print *, '----------------'
C   10 continue
C
C      print 335, tops
C  335 format(' max error for ylm:',1pe12.2)
C      end
C      subroutine makr(rsm,nr,x,y,z)
C      implicit none
C      integer nr,i,ir
C      double precision rs,rsm,x(1),y(1),z(1)
C      real ran1
C      rs = rsm
C      if (rsm < 1d-9) rs = .5d0
C      call ran1in(1)
C      nr = 5
C      do  10  i = 1, nr
C        ir = i+1
C        x(i) = abs((ran1()-.5d0)*5*rs)
C        y(i) = (ran1()-.5d0)*5*rs
C        z(i) = (ran1()-.5d0)*5*rs
C   10 continue
C
C      end
C#endif
