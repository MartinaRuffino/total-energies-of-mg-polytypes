      subroutine ropxvf(m,l,nxi,lxi,vxci,cy,lmax,ioff,xi,cm,v,
     .   q,h,f,g,n,rhoi,job)
c  for some (l,m), adds to vxci and density of this site
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),cy(1),ioff(1),xi(n,0:lmax,nxi),
     .   h(n),f(n),g(n),q(n),v(n),cm(n),vxci(1),rhoi(1),rho(1)
      lav=l*(l+1)+1
c ------- first, make f=qml*cm*vxc ----------------
      do 7 i=1,n
      f(i)=cm(i)*q(i)*v(i)
  7   h(i)=0d0
c ------ start looping over energies -----------
      do 20 ie=1,nxi
      lx=lxi(ie)
      if(l > lx) goto 20
      joff=ioff(ie)
      if(job == 1) then
      sum=0d0
      rrr=rhoi(lav+m+joff)*cy(lav+m)
      do 1 i=1,n
      sum=sum+f(i)*xi(i,l,ie)
  1   h(i)=h(i)+rrr*xi(i,l,ie)
      vxci(lav+m+joff)=vxci(lav+m+joff)+sum*cy(lav+m)
      endif
c -------- repeat for other m --------------
      if(m > 0.and.job == 2) then
      sum=0d0
      rrr=rhoi(lav-m+joff)*cy(lav-m)
      do 3 i=1,n
      sum=sum+f(i)*xi(i,l,ie)
  3   h(i)=h(i)+rrr*xi(i,l,ie)
      vxci(lav-m+joff)=vxci(lav-m+joff)+sum*cy(lav-m)
      endif
  20  continue
c -------- add h*q into g ------------------
      do 40 i=1,n
  40  g(i)=g(i)+h(i)*q(i)

      return
      end
