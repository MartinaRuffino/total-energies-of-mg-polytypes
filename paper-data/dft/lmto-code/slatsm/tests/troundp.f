      subroutine fmain
      implicit none
      logical a2bin
      integer ndig,ip
      double precision x,v,tiny,tiny2,roundp,xr,xr7
      character *50 t

      ndig = 7
C     afmt = '%;7d'
      tiny=1d-4
      tiny2=1d-7

      call info2(0,0,0,' Most calculations use resolution: %g',tiny,0)

      print *, 'A series of 1-decimal numbers:'
      x = 1
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = .3d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -.6d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      print *, '2-decimal numbers:'
      x = 5.01
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)
      x = .25d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      print *, 'Simple fractions:'
      x = 1d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -1d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 2d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -2d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 5d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 8d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 10d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -3d0/7
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      print *, 'Numbers close to simple fractions:'
      x = .6667d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      print *, 'Numbers close to simple fractions:'
      x = -0.428571d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 2.0003d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 2.00003d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 2.000003d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = 2.0000003d0/3
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      print *, 'A fraction with tiny and ndig in sync, out of sync:'
      call snot
      x = -5d0/9.000001d0
      call bin2ax(x,ndig,tiny2,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny2)
      xr7 = xr
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)
      write(*,333), xr7, xr
  333 format(' rounded using 7-decimal resolution:',f15.10,
     .  ';  using 4-decimal:',f15.10)

      print *, 'Numbers involving sqrt(3):'
      x = 5*sqrt(3d0)/12
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -5*sqrt(3d0)/12
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = -sqrt(3d0)/2
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      x = sqrt(3d0)/6
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      call snot


      print *, 'Numbers close to fractions of sqrt(3):'
      x = -0.7216878d0
      call bin2ax(x,ndig,tiny,t)
      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
      xr = roundp(x,tiny)
      call awrit3('   x = %;8g %20pstrng : '//t//
     .  '%50p delta = %g %75px-xr = %g',' ',
     .  100,6,x,v-x,x-xr)

      return

  999 continue
      stop 'something wrong with bin2ax'

      end
      subroutine snot
C     print *, 'hi'
      end
