      subroutine ropxfg(m,l,nxi,lxi,rhoi,cy,lmax,ioff,xi,q,h,f,g,n)
c  for making total density, adds to f (gets multiplied by cm)
c  and to g (gets multiplied by sm)
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),cy(1),rhoi(1),ioff(1),xi(n,0:lmax,nxi),
     .   h(n),f(n),g(n),q(n)
      lav=l*(l+1)+1
      do 20 ie=1,nxi
      lx=lxi(ie)
      if(l > lx) goto 20
      do 1 i=1,n
  1   h(i)=q(i)*xi(i,l,ie)
      joff=ioff(ie)
      cc=cy(lav+m)*rhoi(lav+m+joff)
      do 2 i=1,n
  2   f(i)=f(i)+cc*h(i)
      if(m == 0) goto 20
      cc=cy(lav-m)*rhoi(lav-m+joff)
      do 3 i=1,n
  3   g(i)=g(i)+cc*h(i)
  20  continue
      return
      end
