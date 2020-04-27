      subroutine hyfxsc(nxi,lxi,exi,lmax,rmt,np,nlmx,iperm,xi)
C- Scale to normalize Hankels
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),exi(1),xi(np,nlmx,1),fac(0:20)
      do 10 ie=1,nxi
      lmx=lxi(ie)
      call hrmsav(rmt,lmx,exi(ie),fac)
C|    write(6,300) exi(ie),rmt,(fac(l),l=0,lmx)
C|300 format(' e',f9.4,'   r',f9.4,'   fac',7f9.4)
      do 11 l=0,lmx
      xx=1d0/fac(l)
      do 11 m=0,l
      lm=(l*(l+1))/2+m+1
      if(iperm /= 0) lm=((2*lmax-m+1)*m)/2 + 1 + l
      do 12 ip=1,np
  12  xi(ip,lm,ie)=xi(ip,lm,ie)*xx
  11  continue
  10  continue
      end
