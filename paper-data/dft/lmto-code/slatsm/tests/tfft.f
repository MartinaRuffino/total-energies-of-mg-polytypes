      subroutine fmain
C Convolve brach options : --convolve~off=#~npad=#~reverse
C Test paddding in convolve branch:
C   Copy the following to in1.dat [makes x = i for (-171,171)]
C     % var npts=171 rs=npts/10
C     % repeat i= -npts:npts
C     {i}
C     % end
C   For test, make an array with padding and without.  Cull out padding and compare.
C   a.out --convolve~off=290~npad=37
C   mcx -vl=53 -vnpad=37 out.dat -real -inc 'i<=l|i>l+npad' > snot
C   a.out --convolve~off=290~npad=0
C   mcx snot out.dat -- -px
C   First point should be 119, last point should be 118, whether padded or not
Cu Updates
Cu   30 May 17 Added --convolve switch
C ----------------------------------------------------------------------
      implicit none
      integer iset,isig,i,j,jm,k,plan1(2),plan2(2)
      integer n1,n2,n3
      integer k1,k2,k3
      double complex simple(2,2,4,2),ans(2,2,4,2),ans2(2,2,4,2)
      logical lok
C     integer w(1000)
C     for 1D check
      integer npts,ifi,iv(10),ioff,npad,ireverse,cyclic
      real(8) :: pi = 4*datan(1d0),rs,xv(10)
      character strn*160, dc*1
      procedure(logical) :: a2bin,cmdopt
      procedure(integer) :: fxst,fopna,rdm,cmdoptsw,a2vec,iprint
      procedure(real(8)) :: dlength
      real(8), allocatable :: gr(:)
      complex(8), allocatable :: g1(:),g2(:),c(:),g1save(:),g2save(:)
      cyclic(i,j) = i-floor(dble(i)/j)*j  ! Cyclic modulus: returns n>=0 and n<j

C      npts = 4
C      do  i = -npts, npts
C        jm = npts+2-i; if (i == 1) jm = 1
CC       print *, i, jm, 1+cyclic(1-i,npts)
C        print *, i, cyclic(i,npts)
C      enddo
C      stop

      call suarray(simple,2,2,4)
      call dpzero(ans,2*size(ans))
      ans(1,1,1,1) = dcmplx(1d0,1d0)
      ans(1,1,1,2) = 4d0
      ans(1,1,2,2) = dcmplx(2d0,2d0)
      ans(1,1,4,2) = dcmplx(2d0,-2d0)
      ans(2,1,1,2) = 2d0
      ans(2,1,2,2) = 2d0
      ans(2,1,3,2) = 2d0
      ans(2,1,4,2) = 2d0
      ans(2,2,1,2) = 2d0
      ans(2,2,2,2) = dcmplx(0d0,2d0)
      ans(2,2,3,2) =-2d0
      ans(2,2,4,2) = dcmplx(0d0,-2d0)

C --- Special test: FFT of 1D convolution g1 and g2 if available on disk --
C     Read g1, g2 from disk files in1.dat and in2.dat.  Functions assumed to be real
      if (.not. cmdopt('--convolve',10,0,strn)) goto 99
      dc = strn(11:11); ioff = 0; npad = 0
      if (fxst('in1') /= 1 .or. fxst('in2') /= 1)
     .  call rx('--convolve requires files in1.dat and in2.dat')
      i = cmdoptsw('--convolve','off=','') + 4-1
      if (i > 4-1) then
      if (a2vec(strn,len_trim(strn),i,2,dc//' ',2,1,1,iv,ioff) < 0) goto 911
      endif
      i = cmdoptsw('--convolve','npad=','') + 5-1
      if (i > 5-1) then
        if (a2vec(strn,len_trim(strn),i,2,dc//' ',2,1,1,iv,npad) < 0) goto 911
      endif
      ireverse = 0 ; if (cmdoptsw('--convolve','reverse','') > 0) ireverse = 100
      ifi = fopna('in1',-1,0)
      rewind ifi; npts = 0 ; j = 0

CC     Read g1 from in1.dat
C      call info5(1,0,0,' call suconvolve1D npts=%i ioff=%i l=%i npad=%i',npts,ioff,npts-ioff,npad,5)
C      if (rdm(ifi,0,0,' ',xv,npts,j) < 0) call rx('tfft : file in1.dat not readable')
C      if (j > 1) call rx('tfft : file in1.dat must have only 1 column')
C      allocate(gr(npts*2),g1(npts*2),g1save(npts*2),g2(npts*2),g2save(npts*2))
C      call dpzero(gr,size(gr)); call dpzero(g1,2*size(g1)); call dpzero(g2,2*size(g2))
C      rewind ifi
C      if (rdm(ifi,0,npts*j,' ',gr,npts,j) < 0) call rx('tfft : file in1.dat not readable')
C      call dcopy(npts*2,g1,1,g1save,1)
CC     call prmx('in1.dat',gr,npts,npts,1)
C      call suconvolve1D(ireverse,npts,ioff,npts-ioff,npad,1,gr,g1)
C
CC     debugging check: confirm inverse map recovers original data
CC      gr = 0
CC      call suconvolve1D(ireverse+10,npts,ioff,npts-ioff,npad,1,gr,g1)
CC      call prmx('check that original data restored',gr,npts,npts,1)
C
CC     Read g2 from in2.dat
C      ifi = fopna('in2',-1,0)
C      rewind ifi; k = 0 ; j = 0
C      if (rdm(ifi,0,0,' ',xv,k,j) < 0) call rx('tfft : file in2.dat not readable')
C      if (j > 1) call rx('tfft : file in2.dat must have only 1 column')
C      if (k /= npts) call rx('tfft : file in1.dat and in2.dat must be the same size')
C      rewind ifi
C      if (rdm(ifi,0,npts*j,' ',gr,npts,j) < 0) call rx('tfft : file in2.dat not readable')
C      call suconvolve1D(ireverse,npts,ioff,npts-ioff,npad,1,gr,g2)
C      call info2(1,0,0,' convolve files in1.dat in2.dat, %i points',npts,2)

      call rx0('done')

C --- Default tests ---
   99 continue
      n1 = 2
      n2 = 2
      n3 = 4
      iset = 0
      isig = 1

      call fftz30(n1,n2,n3,k1,k2,k3)

      do  i = 1, 2
C     print *, i; call zprm3('c',0,simple(1,1,1,i),n1,n2,n3)
      call fftz3(simple(1,1,1,i),n1,n2,n3,k1,k2,k3,1,iset,isig)
      print '(a,i2)', ' tfft: check fft for 2x2x4 matrix', i
      CALL diff('compare fftz3 ', simple(1,1,1,i), ans(1,1,1,i),
     .  n1, n2, n3, lok )
      enddo
      call fftzv(simple,n1,n2,n3,1,2,00,isig)

      print '(/1x,a)', 'Checks of fftzv'
C      call suarray(simple,2,2,4)
C      call fftzv(simple,n1,n2,n3,1,2,03,isig)
C      do  i = 1, 2
C      print '(a,i2)', ' tfft: check fft for 2x2x4 matrix', i
C      CALL diff('compare fftz3 ', simple(1,1,1,i), ans(1,1,1,i),
C     .  n1, n2, n3, lok )
C      enddo

      do  k = 0, 2
      iset = 10*k
      print '(/1x,a,i2,''+1, '',i2,''+2, 4'')',
     .  'call fftzv, iset=', iset,iset
      call fftzv(simple,n1,n2,n3,1,2,iset+1,isig)
      call suarray(simple,2,2,4)
      call fftzv(simple,n1,n2,n3,1,2,iset+2,isig)
      do  i = 1, 2
      print '(a,i2)', ' tfft: check fft for 2x2x4 matrix', i
      CALL diff('compare fftz3 ', simple(1,1,1,i), ans(1,1,1,i),
     .  n1, n2, n3, lok )
      enddo
      call fftzv(simple,n1,n2,n3,1,2,4,isig)
      enddo

      print '(/1x,a)', 'call fftzv, iset = 03'
      call suarray(simple,2,2,4)
      call fftzv(simple,n1,n2,n3,1,2,03,isig)
      do  i = 1, 2
      print '(a,i2)', ' tfft: check fft for 2x2x4 matrix', i
      CALL diff('compare fftz3 ', simple(1,1,1,i), ans(1,1,1,i),
     .  n1, n2, n3, lok )
      enddo

      print '(/1x,a)', 'Demonstrate fftzv with mutiple plans'
      call fftzv(simple,n1,n2,n3,1,2,21,-isig)
      call fftzvp(0,plan2)
      call fftzv(simple,n1,n2,n3,1,2,21,isig)
      call fftzvp(0,plan1)

      print '(/1x,a)', 'Assign plan for reverse FT; do it'
      call fftzvp(1,plan2)
      call dcopy(2*2*4*2*2,ans,1,simple,1)
      call fftzv(simple,n1,n2,n3,1,2,22,-isig)
      call suarray(ans2,2,2,4)
      do  i = 1, 2
      print '(a,i2)', ' Check for inverse FT, number', i
      CALL diff('compare fftz3 ', simple(1,1,1,i), ans2(1,1,1,i),
     .  n1, n2, n3, lok )
      enddo

      print '(/1x,a)', 'Assign plan for forward FT; do it'
      call fftzvp(1,plan1)
      call fftzv(simple,n1,n2,n3,1,2,22,isig)
      do  i = 1, 2
      print '(a,i2)', ' Check for inverse FT, number', i
      CALL diff('compare fftz3 ', simple(1,1,1,i), ans(1,1,1,i),
     .  n1, n2, n3, lok )
      enddo

C     Make g1, g2 in-line .. some decaying complex functions
      npts = 343; npad=0
      allocate(g1(npts),g2(npts))
      rs = dble(npts)/10
      do  i = 1, npts
        g1(i) = exp(-dble(30*i)/npts)  * cdexp((0d0,1d0)*dble(i)/npts)
        g2(i) = 1/sqrt(pi)/rs/exp((i/rs)**2) * cdexp((0d0,1d0)*dble(i+4)/2/npts)
      enddo
      do  i = 0, -10, -1
        g1(npts+i) = exp(dble(30*i)/npts)
        g2(npts+i) = 1/sqrt(pi)/rs/exp((i/rs)**2)
      enddo

      if (cmdoptsw('--show','a','') == 0) call pshpr(50)

      if (iprint() >= 50) then
        call zprm('g1',2,g1,npts+npad,npts+npad,1)
        call zprm('g2',2,g2,npts,npts,1)
      endif

C     FFTs of g1,g2
      call fftz30(n1,n2,n3,k1,k2,k3)
C      call fftz3(g1,n1,n2,n3,k1,k2,k3,1,iset,-isig)
C      call zprm('fft g1',2,g1,npts,npts,1)
C      call fftz3(g1,n1,n2,n3,k1,k2,k3,1,iset,isig)
C      call zprm('g1 restored from FFT',2,g1,npts,npts,1)

C      call fftz3(g2,n1,n2,n3,k1,k2,k3,1,iset,-isig)
C      call zprm('fft g2',2,g2,npts,npts,1)

C      call fftz3(c,n1,n2,n3,k1,k2,k3,1,iset,isig)
C      call zprm('convolution(g1,g2)',2,c,npts,npts,1)

C --- Do various kinds of convolution ---
      call info0(1,1,0,' Compare some brute force 1D convolutions to equivalent computed by FFT')
      call info0(1,0,0,' Use --show to print arrays out.dat as they are generated')
      n1 = npts
      n2 = 1
      n3 = 1
      iset = 0
      isig = 1

      allocate(g1save(npts),g2save(npts))
      call dcopy(npts*2,g1,1,g1save,1)
      call dcopy(npts*2,g2,1,g2save,1)

C ... Convolution c(1+cyclic(-i)) = sum_j (g1(1+j)*g2(1+cyclic(i+j))) by brute force
C     Analog of h(-x) = int[ du f(u) g(u+x) ]                           H(k) G(-k)
      allocate(c(npts))
      call dcopy(npts*2,g1save,1,g1,1)
      do  i = 0, npts-1
        k = cyclic(-i,npts)
        c(1+k) = 0  ! c(k) holds h(-x); c(i) holds h(-x)
        do  j = 0, npts-1
          c(1+k) = c(1+k) + g1(1+j)*g2(1+cyclic(i+j,npts))/npts
        enddo
      enddo
      call fftz3c(g1,g2,n1,n2,n3,k1,k2,k3,0,-1)
      if (iprint() >= 50) then
        call zprm('convolution(g1,g2) by brute force',2,c,npts,npts,1)
        call zprm('convolution(g1,g2) by FFT',2,g1,npts,npts,1)
      endif
      g1 = g1-c
      xv(1) = dlength(npts*2,g1,1)
      call info2(1,0,0,' rms difference convolution h(-x) = int[ du f(u) g(u+x) ] : h(2)=%s,(%2;5,5g)  diff*1e16 = %;3d',
     .  c(2),1d16*xv(1))
      deallocate(c)

C ... Convolution c(1+i) = sum_j (g1(1+j)*g2(1+cyclic(-i+j))) by brute force
C     Analog of h(x) = int[ du f(u) g(u-x) ]    H(k) G(-k)
C     Note: because index to x is inverted, c calculated this way is identical to previous test
      allocate(c(npts))
      call dcopy(npts*2,g1save,1,g1,1)
      do  i = 0, npts-1
        c(1+i) = 0
        do  j = 0, npts-1
          c(1+i) = c(1+i) + g1(1+j)*g2(1+cyclic(-i+j,npts))/npts
        enddo
      enddo
      call fftz3c(g1,g2,n1,n2,n3,k1,k2,k3,0,-1)
      if (iprint() >= 50) then
        call zprm('convolution(g1,g2) by brute force',2,c,npts,npts,1)
        call zprm('convolution(g1,g2) by FFT',2,g1,npts,npts,1)
      endif
      g1 = g1-c
      xv(2) = dlength(npts*2,g1,1)
      call info2(1,0,0,' rms difference convolution h(x)  = int[ du f(u) g(u-x) ] : h(2)=%s,(%2;5,5g)  diff*1e16 = %;3d',
     .  c(2),1d16*xv(2))
      call info0(1,0,0,' [Note: previous two tests should be identical since the first array holds discretized h(-x)]')
      deallocate(c)

C ... Convolution c(1+i) = sum_j (g1(1+j)*g2(1+cyclic(i-j))) by brute force
C     Analog of h(x) = int[ du f(u) g(x-u) ]   ->  H(k) G(k)  [standard definiton of convolution]
      allocate(c(npts))
      call dcopy(npts*2,g1save,1,g1,1)
      do  i = 0, npts-1
        c(1+i) = 0
        do  j = 0, npts-1
          c(1+i) = c(1+i) + g1(1+j)*g2(1+cyclic(i-j,npts))/npts
        enddo
      enddo
      call fftz3c(g1,g2,n1,n2,n3,k1,k2,k3,100,-1)
      if (iprint() >= 50) then
        call zprm('convolution(g1,g2) by brute force',2,c,npts,npts,1)
        call zprm('convolution(g1,g2) by FFT',2,g1,npts,npts,1)
      endif
      g1 = g1-c
      xv(3) = dlength(npts*2,g1,1)
      call info2(1,0,0,' rms difference convolution h(x)  = int[ du f(u) g(x-u) ] : h(2)=%s,(%2;5,5g)  diff*1e16 = %;3d',
     .  c(2),xv(3)*1d16)
      deallocate(c)

      if (maxval(xv(1:3)) < 1d-16) then
        call info0(1,0,0,' Direct and FFT convolutions are equivalent')
      else
        call info0(1,0,0,' Direct and FFT convolutions are NOT equivalent')
      endif

      call cexit(0,1)

  911 call rx('usage : --convolve[~origin=#][~npad=#][~reverse]')

      end

      SUBROUTINE diff (strn, sc, pc, n1, n2, n3, lok )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER n1, n2, n3, i, j, k
      logical lok
      double complex sc(n1,n2,n3)
      double complex pc(n1,n2,n3)

      DO i = 1, n1
        DO j = 1, n2
          DO k = 1, n3
            IF ( ABS( sc(i,j,k) - pc(i,j,k) ) .GT. 1.0d-6 ) THEN
               WRITE (6,100) strn
               lok = .false.
               RETURN
               END IF
            END DO
          END DO
        END DO
      WRITE (6,110) strn
      lok = .true.
      RETURN
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results')
  110 FORMAT(1X,a,'... arrays have the same results')
      END

      subroutine zprm3(sfmt,msga,s,n1,n2,n3)
C     implicit none
      character *(*) sfmt
      integer n1,n2,n3,msga
      double precision s(2,n1,n2,n3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna

      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i complex',' ',80,ifi,n1*n2,n3)
      do  10  i = 1, n1
      do  10  j = 1, n2
   10 write(ifi,fmt) (s(1,i,j,k), k=1,n3)
      do  20  i = 1, n1
      do  20  j = 1, n2
   20 write(ifi,fmt) (s(2,i,j,k), k=1,n3)
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in zprm3')
      end

      subroutine suarray(c,n1,n2,n3)
      implicit none
      integer n1,n2,n3
      double complex c(n1,n2,n3,2)


      call dvset(c,1,32,.0625d0)
      call dpzero(c(1,1,1,2),32)
      c(1,1,1,2) = 1d0
      c(1,2,1,2) = 1d0
      c(1,1,2,2) = 1d0
      c(2,2,2,2) = 1d0

      end

      subroutine prmx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(14) fmt, fmt0, strn*(*), outs*80
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#else
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(i,j,1),j=1,nc)
      enddo
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)

C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#else
      outs = ' prm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in prmx'
C#else
      if (outs == 'q') call rx0('quit in prmx')
C#endif
      return

      entry prmx0(fmt0)
      fmt = fmt0
      end

      subroutine zprm(strn,icast,s,ns,nr,nc)
      implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(10) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      fmt = '(9f15.10)'
      fmt = '(9f18.11)'
C     fmt = '(5f20.15)'
      outs = ' '
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex'
      if (icast == 12) outs = ' herm'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#else
      call awrit2('%% rows %i cols %i'//trim(outs),' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(1,i,j),j=1,nc)
      enddo
      if (mod(icast,10) > 1) then
      write(ifi,*)
        do  i = 1, nr
          write (ifi,fmt) (s(2,i,j),j=1,nc)
        enddo
      endif
      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#else
      outs = ' zprm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in zprm'
C#else
      if (outs == 'q') call rx0('quit in zprm')
C#endif
      end

Cc     Copyright (c) 2003, 2007-8 Matteo Frigo
Cc     Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
Cc
Cc     This program is free software; you can redistribute it and/or modify
Cc     it under the terms of the GNU General Public License as published by
Cc     the Free Software Foundation; either version 2 of the License, or
Cc     (at your option) any later version.
Cc
Cc     This program is distributed in the hope that it will be useful,
Cc     but WITHOUT ANY WARRANTY; without even the implied warranty of
Cc     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Cc     GNU General Public License for more details.
Cc
Cc     You should have received a copy of the GNU General Public License
Cc     along with this program; if not, write to the Free Software
Cc     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
Cc
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc
Cc     This is an example implementation of Fortran wisdom export/import
Cc     to/from a Fortran unit (file), exploiting the generic
Cc     dfftw_export_wisdom/dfftw_import_wisdom functions.
Cc
Cc     We cannot compile this file into the FFTW library itself, lest all
Cc     FFTW-calling programs be required to link to the Fortran I/O
Cc     libraries.
Cc
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
Cc     Strictly speaking, the '$' format specifier, which allows us to
Cc     write a character without a trailing newline, is not standard F77.
Cc     However, it seems to be a nearly universal extension.
C      subroutine write_char(c, iunit)
C      character c
C      integer iunit
C      write(iunit,321) c
C 321  format(a,$)
C      end
C
C      subroutine export_wisdom_to_file(iunit)
C      integer iunit
C      external write_char
C      call dfftw_export_wisdom(write_char, iunit)
C      end
C
Cc     Fortran 77 does not have any portable way to read an arbitrary
Cc     file one character at a time.  The best alternative seems to be to
Cc     read a whole line into a buffer, since for fftw-exported wisdom we
Cc     can bound the line length.  (If the file contains longer lines,
Cc     then the lines will be truncated and the wisdom import should
Cc     simply fail.)  Ugh.
C      subroutine read_char(ic, iunit)
C      integer ic
C      integer iunit
C      character*256 buf
C      save buf
C      integer ibuf
C      data ibuf/257/
C      save ibuf
C      if (ibuf < 257) then
C         ic = ichar(buf(ibuf:ibuf))
C         ibuf = ibuf + 1
C         return
C      endif
C      read(iunit,123,end=666) buf
C      ic = ichar(buf(1:1))
C      ibuf = 2
C      return
C 666  ic = -1
C      ibuf = 257
C 123  format(a256)
C      end
C
C      subroutine import_wisdom_from_file(isuccess, iunit)
C      integer isuccess
C      integer iunit
C      external read_char
C      call dfftw_import_wisdom(isuccess, read_char, iunit)
C      end
