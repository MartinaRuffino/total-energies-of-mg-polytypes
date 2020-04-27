      subroutine hyfin2(ncof,nalf,cof,ier,ifi)
c  reads tcf data from file ifi. ier=1 for eof.
      implicit real*8 (a-h,p-z), integer (o)
      dimension cof(ncof*nalf)
      call getpr(ipr)
      if(ipr >= 50) write(6,100) ifi,ncof*nalf*0.001d0
  100 format(' hyfin2 (file',i3,')  read',f5.1,' k-w tcf data')
      ier=0
      read(ifi,end=99) cof
      return
c ----- here if eof reached -----
  99  ier=1
      if(ipr >= 50) write(6,220) ifi
  220 format(' hyfin2:  eof reached on file',i3)
      return
      end
