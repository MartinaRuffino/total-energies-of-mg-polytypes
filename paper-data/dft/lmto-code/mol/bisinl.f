      subroutine bisinl(r1,r2,d,n1,n2,x,z,w,n,nmax,zc1,zc2)
C- Weights and points for bispherical integration
C  of functions with cylindrical symmetry.  Makes 2*n1*n2 points.
C  Adapted from bisint, but uses gaussian quadrature for intgration vars
C  Weights are for integral over three dimensions.
C  For smooth functions, choosing n2 approximately 2*n1+12 works well.
      implicit none
      integer n1,n2,nmax
      double precision x(1),z(1),w(1),d,r1,r2,zc1,zc2
      integer i,iprint,irep,j,m,n
      double precision a,ajac,amumn,amumx,bb,chi1,chi2,chi1mx,chi1mn,
     .  amu,r1p2,r1m2,hta,xxx,pi,rho,sg,
     .  wgt1,wgt2,x1(200),w1(200),x2(200),w2(200)

C --- Setup ---
      if (n1 > 200 .or. n2 > 200) call rx('bisinl: n1 or n2 gt 200')
      call mklegw(n1,x1,w1,iprint())
      call mklegw(n2,x2,w2,iprint())
      pi = 4d0*datan(1d0)
      r1p2 = (r1+r2)/d
      r1m2 = (r1-r2)/d
      a = dsqrt((1-r1p2*r1p2)*(1-r1m2*r1m2))*d/2d0
      chi1mn = -dsqrt(1d0+(a/r1)**2)
      chi1mx = dsqrt(1d0+(a/r2)**2)
      zc1 = r1*chi1mn
      zc2 = r2*chi1mx
c ... acosh(chi1mx)
      amumx = dlog(dabs(chi1mx)+dsqrt(chi1mx*chi1mx-1d0))
      amumn = dlog(dabs(chi1mn)+dsqrt(chi1mn*chi1mn-1d0))

C --- loop over two halves ---
      bb = amumx
      sg = 1d0
      m = 0
      do  50  irep = 1, 2
C ---   Loop over n1 ---
        do  11  j = 1, n1
          amu = bb/2*(x1(j)+1)
          wgt1= bb/2*w1(j)
c ... cosh(amu)
          chi1 = sg*(dexp(amu) + dexp(-amu))*0.5d0
C ---   Loop over n2 ---
          do  10  i = 1, n2
            hta = pi/2*(x2(i)+1)
            wgt2= pi/2*w2(i)
            chi2 = dcos(hta)
            xxx = a/(chi1-chi2)
            rho = dsqrt(1d0-chi2*chi2)*dabs(xxx)
            ajac = xxx*xxx*rho
            m = m+1
            if (m > nmax) call rx('bisinl:  more than nmax points')
            x(m) = rho
            z(m) = dsqrt(chi1*chi1-1d0)*xxx
            w(m) = ajac*(2d0*pi*wgt1*wgt2)
   10     continue
   11   continue
        sg = -1d0
        bb = amumn
   50 continue
      n = m

C --- Printout ---
      if (iprint() < 50) return
      print 745, n1,n2,n,zc1,zc2
  745 format(/' bisinl:  n1,n2=',2i3,'   n=',i6,'   zc1,zc2=',2f11.5)
      if (iprint() < 100) return
      do  60  i = 1, n
   60 write(*,210) x(i),z(i)
  210 format(4f20.15)
C     call rx('bisinl')
      end
C Tests bisinl, bisint
C      subroutine fmain
C      implicit none
C      integer nmax
C      parameter (nmax=100000)
C      integer n1,n2,n
C      double precision x(nmax), z(nmax), w(nmax), out(nmax),
C     .  resl,rest,zc1,zc2,d,r1,r2
C      common /static/ x,z,w,out
C
C      r1 = 1
C      r2 = 1
C      d = 2.1
C      n1 = 25
C      n2 = 25
CC     call pshprt(70)
C
C      print *, 'd='
C      read(*,*) d
C
C
C      print *, 'bisint n1,n2='
C      read(*,*) n1, n2
CC      call bisinl(r1,r2,d,n1,n2,x,z,w,n,nmax,zc1,zc2)
CC      n = n1*n2
CC      CALL XX1(N1,N2,X,Z)
CC      CALL FXX(X,Z,W,N,ZC1,ZC2,OUT)
CC      CALL XX2(N1,N2,OUT)
C
C      call bisint(r1,r2,d,n1,n2,x,z,w,n,nmax,zc1,zc2)
C      call fxx(x,z,w,n,zc1,zc2,out)
C      call dpdot(out,w,n,rest)
C
C      call bisinl(r1,r2,d,n1,n2,x,z,w,n,nmax,zc1,zc2)
C      call fxx(x,z,w,n,zc1,zc2,out)
C      call dpdot(out,w,n,resl)
C
C      print *, 'rest,resl=',n1,n2,rest,resl
C      end
C      subroutine xx1(n1,n2,x,z)
C      double precision x(1),z(1)
C
C      ip = 0
C      do  10  i1 = 1, n1
C        do  10  i2 = 1, n2
C          ip = ip+1
C          x(ip) = -3 + 6*dble(i1-1)/(n1-1)
C          z(ip) = -3 + 6*dble(i2-1)/(n2-1)
C   10 continue
C      end
C      subroutine xx2(n1,n2,out)
C      double precision out(1)
C      print *, n1,n2
C      do  10  i1 = 1, n1*n2
C   10 print *, out(i1)
C      stop
C      end
C      subroutine fxx(x,z,w,n,zc1,zc2,res)
C      implicit none
C      integer n,ip
C      double precision x(1), z(1), w(1),
C     .  res(1),zc1,zc2,z1,z2,xi(10),xx,r1,r2
C
C
C      res(1) = 0
C      do  10  ip = 1, n
C
C        z1 = z(ip)-zc1
C        z2 = z(ip)-zc2
C        r1 = dsqrt(z1**2+x(ip)**2)
C        r2 = dsqrt(z2**2+x(ip)**2)
C
CC        res = res + w(ip)*
CC     .    (dexp((-(x(ip)**2 + 1.1d0*(z(ip)+zc1)**2)/2)) +
CC     .     dexp((-(x(ip)**2 + 0.9d0*(z(ip)+zc2)**2)/2)) )
C        call hansmr(r1,-1d0,2d0,xi,1)
C        xx = xi(1)*x(ip)
C        call hansmr(r2,-1d0,2d0,xi,1)
C        xx = xx*xi(1)
C        res(ip) = xx
C   10 continue
C      end
