      subroutine ropxfv(m,l,nxi,lxi,vxci,cy,lmax,ioff,xi,cm,sm,v,
     .   q,h,f,g,n)
c  for some (l,m), adds to vxci
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),cy(1),ioff(1),xi(n,0:lmax,nxi),
     .   h(n),f(n),g(n),q(n),v(n),cm(n),sm(n),vxci(1)
      lav=l*(l+1)+1
c ------- first, make f=v*qml*cm, g=v*qlm*sm --------
      do 7 i=1,n
      xx=q(i)*v(i)
      f(i)=xx*cm(i)
  7   g(i)=xx*sm(i)
c ------ start looping over energies -----------
      do 20 ie=1,nxi
      lx=lxi(ie)
      if(l > lx) goto 20
      joff=ioff(ie)
      sum=0d0
      do 1 i=1,n
  1   sum=sum+f(i)*xi(i,l,ie)
      vxci(lav+m+joff)=vxci(lav+m+joff)+sum*cy(lav+m)
      if(m == 0) goto 20
      sum=0d0
      do 3 i=1,n
  3   sum=sum+g(i)*xi(i,l,ie)
      vxci(lav-m+joff)=vxci(lav-m+joff)+sum*cy(lav-m)
  20  continue
      return
      end
