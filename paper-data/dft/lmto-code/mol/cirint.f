      subroutine cirint(r1,r2,n1,n2,zc,fac,x,z,w,n,nmax)
C-  Integration points and weights for circle r1 to r2, around zc
      implicit real*8 (a-h,p-z), integer(o)
      dimension x(1),z(1),w(1),x1(100),w1(100),x2(100),w2(100)
      pi=4.d0*datan(1.d0)
      n=0
      if(n1 <= 0.or.n2 <= 0) return
      if(n1 > 100.or.n2 > 100) call rx('cirint: n1 or n2 gt 200')
      call mklegw(n1,x1,w1,000)
      call mklegw(n2,x2,w2,000)
      a=(r2-r1)/2d0
      b=(r2+r1)/2d0
      wfac=2d0*pi*fac
      m=0
      do 10 j=1,n1
      rad=a*x1(j)+b
      wad=a*w1(j)
      do 10 i=1,n2
      hta=0.5d0*pi*(x2(i)+1d0)
      wta=0.5d0*pi*w2(i)
      m=m+1
      if(m > nmax) call rx('cirint: more than nmax points''')
      x(m)=rad*dsin(hta)
      z(m)=rad*dcos(hta)+zc
  10  w(m)=wfac*wad*wta*x(m)*rad
      n=m
      end
