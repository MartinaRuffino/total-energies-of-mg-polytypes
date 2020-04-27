      subroutine headl3(name,vrsn,wksiz,ifi)
C-  Puts a heading line into file, unit ifi
      implicit none
      integer ifi,wksiz
      double precision vrsn
      character*8 name,outs*80

      outs = ' '
      if (vrsn == 0) call headl2(name,wksiz,ifi)

      if (wksiz == 0)
     .  call awrit1(' -----------------------  START '
     .  //trim(name)//'%a (%d)  -----------------------',
     .  outs,-80,-ifi,vrsn)
      if (wksiz /= 0)
     .  call awrit2(' -----------------------  START '
     .  //name//'%a (v%d, %iK)  -----------------------',
     .  outs,-80,-ifi,vrsn,(wksiz+499)/1000)
      end

      subroutine headl2(name,wksiz,ifi)
C- Puts a heading line into file, unit ifi
      implicit none
      integer ifi,wksiz
      character*8 name,outs*80
      character(len=16) :: gitref
! this is sily but surprisingly it was not a simple matter to repeat a character, maybe awrite can do it?
      character(len=64), parameter :: dashes = '----------------------------------------------------------------'
      integer :: n, nf, ne

      gitref = ''
      include 'gitref.h'

      if (wksiz == 0) then
        if (len_trim(gitref) > 0 .and. index(gitref, 'v') /= 1) gitref = '- ref: '//trim(gitref)
        n = 48 - len_trim(name) - len_trim(gitref)
        nf = n/2
        ne = nf; if (mod(n,2) == 1) ne = ne + 1

        write(ifi,'(x,a,2x,"START",x,a,x,a,x,a)') dashes(1:nf), trim(name), trim(gitref), dashes(1:ne)
      else
        outs = ' '
        call awrit1(' -----------------------  START '//trim(name)//'%a (%iK)  -----------------------',
     .  outs,-80,-ifi,(wksiz+499)/1000)
      end if
      end
      subroutine poseof(iunit)
C- Positions file handle at end-of-file
      implicit none
C Passed parameters
      integer iunit
C Local parameters
      integer i,nrec
C#ifdefC SVS | DEC | DECA
C      double precision dummy
C#endif

      nrec = 0
      rewind iunit
      do  10  i = 1, 100000000
C ... Avoid read bug ...
C#ifdefC SVS | DEC | DECA
C        read(iunit,100,end=90,err=91) dummy
C#endif
C ... Avoid AIX bug
C#ifdefC AIX
C        read(iunit,100,end=90,err=90)
C#else
        read(iunit,100,end=90,err=91)
C#endif
        nrec = i
  100   format(a1)
   10 continue
      write(*,200) iunit
  200 format(' POSEOF: no EOF found for file',i3)
      return
   90 continue
C ... Avoid AIX bug in xlf90 compiler
C#ifndef AIX-xlf90
      backspace iunit
C#endif
C ... If backspace doesn't work, do it this way ...
C     rewind iunit
C     do 11 i=1,nrec
C  11 read(iunit,100)
   91 continue
      end
