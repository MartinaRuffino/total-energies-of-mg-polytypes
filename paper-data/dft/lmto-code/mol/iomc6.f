      subroutine iomc6(ipr,id,tabs,symgrp,nel,el,lswtch,nit1,nit2,beta1,
     .  beta2,tol,amp,namp,nfile,ntx,ndx,alfa,alat,plat,nft1,nft2,nft3,
     .  dx,dy,dz,nit0,lrx,gamrx,ldyn,tau,fric,temp,lrs,ifi)
c  read general data from control file
      implicit real*8 (a-h,p-z), integer (o)
      dimension el(10),lswtch(20),plat(9),amp(15),lrs(4)
      character t*73,w*10,wl*10,ttt*60,tabs*80
      character*1 id(60,2),symgrp(60)
      wl=' '
      nid=1
      rewind ifi
  90  call inline(w,t,leof,ifi)
      t(73:73)=' '
      if(leof == 1 .or. w == 'exit') goto 88
      if(w == '"') w=wl
c -------- files -----------
      if(w == 'files') call strcop(tabs,t,60,'#',i)
c -------- cntrl,xcmesh,elmto,tables,switch,amp -----------
      if(w == 'cntrl')  read(t,*,end=9) ipr,nit1,nit2,beta1,beta2,tol
      if(w == 'elmto')  call rdcnt(10,nel,el,t)
      if(w == 'amp')    call rdcnt(15,namp,amp,t)
      if(w == 'xcmesh') read(t,*,end=9) dx,dy,dz
      if(w == 'ftmesh') read(t,*,end=9) alat,plat,nft1,nft2,nft3
C      if(w == 'ftmesh'.and.mod(nft1,2) /= 0)
C     .  call rx('iomc: ftmesh expects nft1 even')
      if(w == 'tables') read(t,*,end=9) nfile,ntx,ndx,alfa
      if(w == 'switch') read(t,*,end=9) (lswtch(i),i=1,9)
      if(w == 'move')   read(t,*,end=9) nit0,lrx,gamrx,ldyn,
     .   tau,fric,temp
      if(w == 'rsta')   read(t,*,end=9) (lrs(i),i=1,3)
c -------- id ------------------------------
      if(w == 'id') then
C         read(t,*,end=9) ttt
         call strcop(ttt,t,60,'#',i)
         do 1 i=1,60
         if(nid == 1) id(i,2)=' '
  1      id(i,nid)=ttt(i:i)
         nid=2
         endif
c -------- symgrp --------------------------
      if(w == 'symgrp') then
         read(t,'(a)',end=9) ttt
         do 2 i=1,60
  2      symgrp(i)=ttt(i:i)
         endif

      wl=w
      goto 90
c ---- post-procesing and return ---------
  88  continue
      if(namp == 0) then
        namp=1
        amp(1)=0d0
        endif
      return
c -------- error exit ----------
  9   call rx('iomc6: eol reached; 1st word='//w)

      end
c ------- sub rdcnt: read array and count -------
      subroutine rdcnt(na0,na,a,t)
      implicit real*8 (a-h,p-z), integer (o)
      dimension a(na0)
      character*73 t
      do 8 i=1,na0
  8   a(i)=1d21
      t(73:73)='/'
      read(t,*) a
      na=0
      do 7 i=1,na0
  7   if(a(i) < 1d20) na=i
      end
