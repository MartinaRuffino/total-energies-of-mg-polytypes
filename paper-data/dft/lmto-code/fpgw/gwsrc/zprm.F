      subroutine zprm(strn,s,ns,nr,nc)
      implicit none
      integer ns,nr,nc
      double precision s(2,ns,nc)
      character*(14) fmt, strn*(*), outs*80
      integer i,j,ifi

      fmt = '(5f20.14)'
      fmt = '(1p9e22.12)'
      ifi = 19
      open(ifi,file='out')
      rewind ifi
      write(ifi,234) nr, nc
  234 format('% rows',i6,' cols',i6,' complex')
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(1,i,j), j=1,nc)
      do  20  i = 1, nr
   20 write(ifi,fmt) (s(2,i,j), j=1,nc)
      close(ifi)
      outs = 'zprm: done writing data '//strn
      print *, outs
      read(*,'(a80)') outs
      if (outs == 'q') stop 'quit in zprm'
      end

      subroutine dprm(strn,s,ns,nr,nc)
      implicit none
      integer ns,nr,nc
      double precision s(ns,nc)
      character*(10) fmt, strn*(*), outs*80
      integer i,j,ifi

      fmt = '(5f18.13)'
      ifi = 19
      open(ifi,file='out')
      rewind ifi
      write(ifi,234) nr, nc
  234 format('% rows',i5,' cols',i5)
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j), j=1,nc)
      close(ifi)
      outs = 'zprm: done writing data '//strn
      print *, outs
      read(*,'(a80)') outs
      if (outs == 'q') stop 'quit in zprm'
      end
