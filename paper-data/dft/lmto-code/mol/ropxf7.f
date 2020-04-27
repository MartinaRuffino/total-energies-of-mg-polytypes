      subroutine ropxf7(m,l,nxi,lxi,vxci,cy,lmax,ioff,xi,cm,sm,v,
     .   q,h,f,g,n,dx,dy,dz,dvxcx,dvxcy,dvxcz)
c  for some (l,m), adds to vxci
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),cy(1),ioff(1),xi(n,0:lmax,nxi),
     .   h(n),f(n),g(n),q(n),v(n),cm(n),sm(n),vxci(1),
     .   dx(n),dy(n),dz(n),dvxcx(1),dvxcy(1),dvxcz(1)
      lav=l*(l+1)+1
c ------- first, make f=qml*cm, g=qlm*sm --------
      do 7 i=1,n
      f(i)=cm(i)*q(i)
  7   g(i)=sm(i)*q(i)
c ------ start looping over energies -----------
      do 20 ie=1,nxi
      lx=lxi(ie)
      if(l > lx) goto 20
      joff=ioff(ie)
      sum=0d0
      sux=0d0
      suy=0d0
      suz=0d0
      do 1 i=1,n
      xx=f(i)*xi(i,l,ie)
      sum=sum+xx*v(i)
      sux=sux+xx*dx(i)
      suy=suy+xx*dy(i)
  1   suz=suz+xx*dz(i)
      vxci(lav+m+joff)=vxci(lav+m+joff)+sum*cy(lav+m)
      dvxcx(lav+m+joff)=dvxcx(lav+m+joff)+sux*cy(lav+m)
      dvxcy(lav+m+joff)=dvxcy(lav+m+joff)+suy*cy(lav+m)
      dvxcz(lav+m+joff)=dvxcz(lav+m+joff)+suz*cy(lav+m)
c -------- repeat for other m --------------
      if(m == 0) goto 20
      sum=0d0
      sux=0d0
      suy=0d0
      suz=0d0
      do 3 i=1,n
      xx=g(i)*xi(i,l,ie)
      sum=sum+xx*v(i)
      sux=sux+xx*dx(i)
      suy=suy+xx*dy(i)
  3   suz=suz+xx*dz(i)
      vxci(lav-m+joff)=vxci(lav-m+joff)+sum*cy(lav-m)
      dvxcx(lav-m+joff)=dvxcx(lav-m+joff)+sux*cy(lav-m)
      dvxcy(lav-m+joff)=dvxcy(lav-m+joff)+suy*cy(lav-m)
      dvxcz(lav-m+joff)=dvxcz(lav-m+joff)+suz*cy(lav-m)
  20  continue
      return
      end
