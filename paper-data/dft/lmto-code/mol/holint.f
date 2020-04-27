      subroutine holint(r1,r2,d,n1,n2,x,z,w,n,nmax,zc1,zc2)
C-  Points and weights for 3d integration with holes (bispherical)
      implicit real*8 (a-h,p-z), integer(o)
      dimension x(9),z(9),w(9)
      fac=-1d0
      s1=0.8d0*r1
      s2=0.8d0*r2
      m1=7
      m2=21
      call bisinl(s1,s2,d,n1,n2,x,z,w,k0,nmax,zc1,zc2)
      n=k0
      call cirint(s1,r1,m1,m2,zc1,fac,x(n+1),z(n+1),w(n+1),k1,nmax-n)
      n=n+k1
      call cirint(s2,r2,m1,m2,zc2,fac,x(n+1),z(n+1),w(n+1),k2,nmax-n)
      n=n+k2
      if(iprint() >= 30) write(6,731) m1,m2,k0,k1,k2,n
  731 format(' holint:  m1,m2=',2i3,'   n0,n1,n2=',3i5,'   n=',i5)
      end
