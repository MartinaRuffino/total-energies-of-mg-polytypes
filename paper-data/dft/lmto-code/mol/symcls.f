      subroutine symcls(ng,g,ag,nbas,pos,ips,spid,nclas,ipclas)
      implicit real*8 (a-h,p-z), integer (o)
      parameter( nsv=100 )
      dimension qb(3,3),g(9,1),ag(3,1),ips(1),pos(3,1),isv(nsv),
     .   ipclas(nbas),d(3)
      character*8 spid(1)
      call getpr(ipr)
      do 1 ib=1,nbas
  1   ipclas(ib)=0
      nclas=0
c ---- find next unclassified atom ---------
      if(ipr >= 30) write(6,838)
  71  ib=0
      do 2 jb=nbas,1,-1
  2   if(ipclas(jb) == 0) ib=jb
      if(ib == 0) goto 77
c ---- apply group ops, check that atom of same species there ---
      nclas=nclas+1
      ipclas(ib)=nclas
      isv(1)=ib
      nrc=1
      do 10 ig=1,ng
      call gpfnd1(g(1,ig),ag(1,ig),pos(1,ib),nbas,pos,d,jb)
      lpr=0
      if(ipr >= 50.or.jb == 0) lpr=1
      if(lpr == 1) write(6,200) ib,ig,(pos(m,ib),m=1,3),d,jb
  200 format(' ib',i3,'  op',i3,3f9.4,'  -->',3f9.4,'  jb',i3)
      if(jb == 0) call rx('symcls: no atom found at g(pos)')
      if(ips(jb) /= ips(ib)) write(6,200) ib,ig,(pos(m,ib),m=1,3),d,jb
      if(ips(jb) /= ips(ib)) call rx('symcl: wrong species at g(pos)')
      if(ipclas(jb) == 0) nrc=nrc+1
      if(nrc > nsv) call rx('symcls: incr nsv')
      if(ipclas(jb) == 0) isv(nrc)=jb
  10  ipclas(jb)=nclas
      if(ipr >= 30) write(6,837) nclas,spid(ips(ib)),nrc,
     .   (isv(j),j=1,nrc)
  837 format(i4,4x,a4,i5,3x,15i4:(20x,15i4))
  838 format(/' class  spec  size    atoms ...')
      goto 71
  77  continue
      end
c --------- gpfndz --------------------------------
      subroutine gpfndz(g,ag,p,nbas,pos,d,jb)
c  finds atom jb at map of p under group operation (g,ag).
      implicit real*8 (a-h,p-z), integer(o)
      dimension g(3,3),ag(3),pos(3,1),d(3),p(3)
      do 4 m=1,3
      d(m)=ag(m)
      do 4 n=1,3
  4   d(m)=d(m)+g(m,n)*p(n)
      jb=0
      do 5 k=1,nbas
      ddd=(pos(1,k)-d(1))**2+(pos(2,k)-d(2))**2+(pos(3,k)-d(3))**2
  5   if(dsqrt(ddd) < 1d-5) jb=k
      end
