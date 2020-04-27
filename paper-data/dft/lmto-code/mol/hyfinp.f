      subroutine hyfinp(nfile,tcspec,tcdata,ntmx,ndmx)
c  reads in tcf tables from files 61,62 ...
      implicit real*8 (a-h,p-z), integer (o)
      character*1 ch(0:1)
      character*8 ctype(0:3)
      dimension tcspec(100,ntmx),tcdata(ndmx)
      data ch /' ','s'/
      data ctype/ '2c-unsym','2c-sym  ','   ?    ','1-center'/
      call getpr(ipr)
      ntbl=0
      ndata=0
      nalfx=0
      do 10 ifile=1,nfile
      ifi=60+ifile
      rewind ifi
      do 11 irep=1,99
      call hyfin0(ncof,nalf,tcspec(1,ntbl+1),ier,ifi)
      if(ier /= 0) goto 92
      ntbl=ntbl+1
      nalfx=max0(nalfx,nalf)
      if(ntbl > ntmx) call rx('hyfinp: too many tables')
      tcspec(25,ntbl)=ndata+1
      call hyfin2(ncof,nalf,tcdata(ndata+1),ier,ifi)
      ndata=ndata+ncof*nalf
      if(ndata > ndmx) call rx('hyfinp: too much data')
      if(ier /= 0) call rx('error reading tcdata')
  11  continue
  92  continue
  10  continue
      i=min0(iabs(ntbl-1),1)
      j=min0(iabs(nfile-1),1)

      if(ipr >= 20) write(6,924) ntbl,ch(i),ndata/1000d0,nfile,ch(j)
  924 format(/' hyfinp:',i4,' table',a1,'  (',f6.1,
     .  ' kw  of data )  loaded from',i3,' file',a1)
CL      write(71,711) ntbl,nfile,ndata/1000d0,ndmx/1000d0,nalfx
  711 format(' tables',i3,'   files',i3,'   data',f6.1,'  of',
     .   f6.1,' kw    max nalf',i4)
      tcspec(26,1)=ntbl
c ----- output some specifications ----
      if(ipr < 30) return
      write(6,400)
      do 30 it=1,ntbl
      r1=tcspec( 1,it)
      r2=tcspec( 2,it)
      ep1=tcspec(16,it)
      ep2=tcspec(18,it)
      lp1=idnint( tcspec(15,it) )
      lp2=idnint( tcspec(17,it) )
      lsym=idnint( tcspec(13,it) )
      ncof=idnint( tcspec( 9,it) )
      nalf=idnint( tcspec( 7,it) )
      size=ncof*nalf*0.001d0
      write(6,330) it,r1,r2,ep1,ep2,lp1,lp2,nalf,size,ctype(lsym)
  330 format(i3,1x,2f9.4,2f8.2,1x,2i3,i5,f7.1,3x,a8)
  30  continue
  400 format(/' tbl',5x,'r1',7x,'r2',7x,'ep1',5x,'ep2',
     .   3x,'lphi   nalf  kw    type')


      return
      end
