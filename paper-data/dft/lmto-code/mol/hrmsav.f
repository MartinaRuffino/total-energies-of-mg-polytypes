      subroutine hrmsav(r,lmax,e,hrms)
c  rms average of hankels over sphere with radius r
      implicit real*8 (a-h,p-z), integer (o)
      dimension hrms(0:lmax),phi(0:20),psi(0:20)
c|    write(6,898) r,e,lmax
c|898 format(' r=',f12.6,'   e=',f12.6,'   lmax=',i5)
      pi=4d0*datan(1d0)
      cxx=1d0/dsqrt(4d0*pi)
      call bessl(e*r*r,lmax,phi,psi)
      rfac=1.d0
      do 20 l=0,lmax
      rfac=rfac*(1d0/r)
 20   hrms(l)=rfac*cxx*dabs(psi(l))
      return
      end
