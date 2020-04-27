      subroutine iomc4(ipr,id,symgrp,nel,el,lswtch,nit1,nit2,
     .  beta1,beta2,tol,amp,nfile,ntx,ndx,alfa,plat,nft1,nft2,nft3,
     .  dx,dy,dz,nc0,atr,tolch1,tolch2,ifi)
c  read general data from control file
      implicit real*8 (a-h,p-z), integer (o)
      dimension el(1),lswtch(20),plat(9)
      character t*72,w*10,wl*10,ttt*60
      character*1 id(60,2),symgrp(60)
      wl=' '
      nid=1
      rewind ifi
  90  call inline(w,t,leof,ifi)
      if(leof == 1 .or. w == 'exit') return
      if(w == '"') w=wl
c -------- cntrl,elmto,xcmesh,tables,switch,amp -----------
      if(w == 'cntrl')  read(t,*,end=9) ipr,nit1,nit2,beta1,beta2,tol
      if(w == 'elmto')  read(t,*,end=9) nel,(el(i),i=1,nel)
      if(w == 'xcmesh') read(t,*,end=9) dx,dy,dz,nc0,atr,tolch1,tolch2
      if(w == 'ftmesh') read(t,*,end=9) plat,nft1,nft2,nft3
C      if(w == 'ftmesh'.and.mod(nft1,2) /= 0)
C     .  call rx('iomc: ftmesh expects nft1 even')
      if(w == 'tables') read(t,*,end=9) nfile,ntx,ndx,alfa
      if(w == 'switch') read(t,*,end=9) (lswtch(i),i=1,3)
      if(w == 'amp')    read(t,*,end=9) amp
c -------- id ------------------------------
      if(w == 'id') then
         read(t,*,end=9) ttt
         do 1 i=1,60
         if(nid == 1) id(i,2)=' '
  1      id(i,nid)=ttt(i:i)
         nid=2
         endif
c -------- symgrp --------------------------
      if(w == 'symgrp') then
         read(t,*,end=9) ttt
         do 2 i=1,60
  2      symgrp(i)=ttt(i:i)
         endif
c -----------------------------------------
      wl=w
      goto 90
  9   call rx('rdcfil: eol reached; 1st word='//w)

      end
