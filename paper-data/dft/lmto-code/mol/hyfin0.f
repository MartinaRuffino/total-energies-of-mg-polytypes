      subroutine hyfin0(ncof,nalf,t,ier,ifi)
c  reads tcf specifications from file ifi. ier=1 for eof.
c  specification vector is returned in t
      implicit real*8 (a-h,p-z), integer (o)
      dimension lx1(8),ex1(8),lx2(8),ex2(8),t(100),err(20)
      call getpr(ipr)
c ------- read t from file -----------------------
      ier=0
      read(ifi,end=99) t
      if(ipr >= 50) write(6,100) ifi
  100 format(' hyfin0 (file',i3,')  read in tcf specs')
      if(ipr >= 70) write(6,200) t
  200 format(1x,10f8.2)
      nalf        =idnint( t( 7) )
      ncof        =idnint( t( 9) )
      if(ipr >= 60) write(6,886) nalf,ncof
  886 format(' nalf=',i5,'   ncof=',i5)
      return
c ----- here if eof reached -----
  99  ier=1
      if(ipr >= 50) write(6,220) ifi
  220 format(' hyfin0 (file',i3,')  end-of-file')
      return
      end
