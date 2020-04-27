      double precision function roundp(x,tiny)
C- Rounds a number of order unity, close to a rational number
C  or some simple rational number * sqrt(3)
C ----------------------------------------------------------------------
Ci Inputs
Ci   x     :double precision number to be rounded
Ci   tiny  :resolution in number
Co Outputs
Co  roundp :a number within tiny of x, but close to a simple rational
Co         :number, or rational number * sqrt(3)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 May 12 Rounds numbers less than 'small' to zero
Cu   10 Jul 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndig
      double precision x,tiny
C ... Local parameters
      integer j,iz,id
      double precision y,d,big,small
      parameter (small=1.000000000000001d-6)

C     Default: original number
      roundp = x

C     Number close to zero
      if (tiny <= small .and. abs(x) <= tiny) then
        roundp = 0
        return
      endif

C     Do not massage large numbers or numbers less than tiny
      if (abs(x) < tiny .or. x > 5d0) return
C     Do not massage numbers fewer decimals than -log(tiny)
      ndig = nint(dlog10(tiny))
      big = 10d0**(-ndig)
      y = big*abs(x)
      if (abs(y-nint(y)) < tiny) then
C       print *, 'bye', x,y,abs(y)-abs(nint(y))
        return
      endif

      y = x*dsqrt(3d0)*4
      do  j = 1, 9
C ... If close to x/j, j=1..9
        if (abs(j/x-nint(j/x)) < tiny) then
          id = nint(j/x)
          roundp = dble(j)/dble(id)
          exit
C ... If close to simple fraction * sqrt(3)
        elseif (abs(y-nint(y)) < tiny .and.
     .      abs(y) > .5d0) then
          d = 12
          iz = nint(abs(y))
          id = 12
          if (iz == 2) id = 6
          if (iz == 4) id = 3
          if (iz == 6) id = 2
          if (iz == 12) id = 1
          y = y/(12/id)
          roundp = nint(y)*dsqrt(3d0)/dble(id)
          exit
        endif
      enddo

      end
      subroutine bin2ax(x,ndig,tiny,sout)
C- Converts binary number of order unity to ascii string,
C  rounding numbers close to simple fractions, or sqrt(3)*fractions
C ----------------------------------------------------------------------
Ci Inputs
Ci   x     :double precision number to be rounded
Ci   ndig  :number of digits to print out
Ci   tiny  :resolution in number.  Normally tiny = 1/10^ndig
Co Outputs
Co  sout   :a ascii representation of number within tiny of x, but close
Co         :to a simple rational number, or rational number * sqrt(3)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Jul 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndig
      double precision x,tiny
      character sout*(*)
C ... Local parameters
      integer j,iz,id
      double precision y,d
      character afmt*(20)

C     Initially write in standard format
      afmt = ' '
      call awrit1('%%;%id',afmt,len(afmt),0,ndig)
      sout = ' '
      call awrit1(afmt,sout,len(sout),0,x)

C     Do not massage small or large numbers
      if (abs(x) < tiny .or. x > 5d0) return
      y = 10*abs(x)
C     Do not massage numbers with only one decimal
      if (abs(y/nint(y)-1) < tiny) then
        return
      endif

      y = x*dsqrt(3d0)*4
      do  j = 1, 9
C ... If close to x/j, j=1..9
        if (abs(j/x-nint(j/x)) < tiny) then
          if (x > 0) then
            call awrit2('%x%i/%d',sout,len(sout),0,j,j/x)
            exit
          endif
          if (x < 0) then
            call awrit2('%x-%i/%d',sout,len(sout),0,j,-j/x)
            exit
          endif
C ... If close to simple fraction * sqrt(3)
        elseif (abs(y-nint(y)) < tiny .and.
     .      abs(y) > .5d0) then
          d = 12
          iz = nint(abs(y))
          id = 12
          if (iz == 2) id = 6
          if (iz == 4) id = 3
          if (iz == 6) id = 2
          if (iz == 12) id = 1
          y = y/(12/id)
          if (y > 0) then
            call awrit3('%x%?#n==1#%j#%d*#sqrt(3)/%i',
     .        sout,len(sout),0,nint(y),y,id)
          else
            call awrit3('%x-%?#n==1#%j#%d*#sqrt(3)/%i',
     .        sout,len(sout),0,nint(-y),-y,id)
          endif
          exit
        endif
      enddo
      end
C      subroutine fmain
C      implicit none
C      logical a2bin
C      integer ndig,ip
C      double precision x,v,tiny,tiny2,roundp,xr,xr7
C      character *50 t
C
C      ndig = 7
CC     afmt = '%;7d'
C      tiny=1d-4
C      tiny2=1d-7
C
C      call info2(0,0,0,' Most calculations use resolution: %g',tiny,0)
C
C      print *, 'A series of 1-decimal numbers:'
C      x = 1
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = .3d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -.6d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      print *, '2-decimal numbers:'
C      x = 5.01
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C      x = .25d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      print *, 'Simple fractions:'
C      x = 1d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -1d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 2d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -2d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 5d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 8d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 10d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -3d0/7
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      print *, 'Numbers close to simple fractions:'
C      x = .6667d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      print *, 'Numbers close to simple fractions:'
C      x = -0.428571d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 2.0003d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 2.00003d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 2.000003d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = 2.0000003d0/3
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      print *, 'A fraction with tiny and ndig in sync, out of sync:'
C      call snot
C      x = -5d0/9.000001d0
C      call bin2ax(x,ndig,tiny2,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny2)
C      xr7 = xr
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C      write(*,333), xr7, xr
C  333 format(' rounded using 7-decimal resolution:',f15.10,
C     .  ';  using 4-decimal:',f15.10)
C
C      print *, 'Numbers involving sqrt(3):'
C      x = 5*sqrt(3d0)/12
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -5*sqrt(3d0)/12
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = -sqrt(3d0)/2
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      x = sqrt(3d0)/6
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      call snot
C
C
C      print *, 'Numbers close to fractions of sqrt(3):'
C      x = -0.7216878d0
C      call bin2ax(x,ndig,tiny,t)
C      ip = 0; if (.not. a2bin(t,v,4,0,' ',ip,-1)) goto 999
C      xr = roundp(x,tiny)
C      call awrit3('   x = %;8g %20pstrng : '//t//
C     .  '%50p delta = %g %75px-xr = %g',' ',
C     .  100,6,x,v-x,x-xr)
C
C      return
C
C  999 continue
C      stop 'something wrong with bin2ax'
C
C      end
C      subroutine snot
CC     print *, 'hi'
C      end
