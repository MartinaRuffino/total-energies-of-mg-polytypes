      subroutine prm(sfmt,msga,s,ns,nr,nc)
C- writes matrix into file with message (for debugging)
      implicit none
      integer nr,nc,ns,ifi,msga
      double precision s(ns,nc,2)
      character*(10) fmt, sfmt*(*), outs*80
      integer i,j,fopna,i1mach
      fmt = '(9f22.17)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
      do  i = 1, nr
        write (ifi,fmt) (s(i,j,1),j=1,nc)
      enddo
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in prmx')
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
      subroutine prmi(strn,s,ns,nr,nc)
C- writes integer matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      integer s(ns,nc)
      character*(14) fmt, strn*(*), outs*80
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(20i6)'/
C      fmt = '(1p9e20.10)'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#else
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(i,j),j=1,nc)
      enddo
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
      end

      subroutine prmr(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(10) fmt, strn*(*)
      integer i,j,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      write(ifi,*) nr, nc
      do  i = 1, nr
        write (ifi,fmt) (s(i,j,1),j=1,nc)
      enddo
      write(ifi,*)
      call fclose(ifi)
      print *, 'prm: pausing after writing data ',strn
      pause
      end
      subroutine yprmi(fmt,i1,i2,icast,s,ofi,ns,nr,nc)
C- Prints complex matrix, with formatting st
C- Formats and prints a string, followed by the contents of a complex matrix
      implicit none
      character*(*) fmt
      integer i1,i2,icast,ofi,ns,nr,nc
      double precision s(ns,nc,2)
      character*(120) strn

      strn = ' '
      call awrit2(fmt,strn,len(strn),0,i1,i2)
      call yprm(strn,icast,s,ofi,ns,nr,nc)

      end

      subroutine yprm(strn,icast,s,ofi,ns,nr,nc)
C- Prints a string, followed by the contents of a complex matrix
C ofi used only for kcplx=0
C ns,nr,nc are formal dimensions, not real ones
Ci  icast:   0 integer
Ci           1 double precision
Ci           2 double complex with imaginary following real
Ci           3 double complex
Ci           4 double complex with imaginary following real in columns
      implicit none
      integer icast,ofi,ns,nr,nc,ifi
*      double precision s(ns,nsc,2)
      integer s(ns,nc)
      character*(20) fmt, fmt0*(*), outs*120, strn*(*)
      integer i,j,fopna,i1mach
      save fmt
C     data fmt /'(1p9e20.12)'/
      data fmt /'(9f15.10)'/
C     data fmt /'(%9;10,10d)'/

      ifi = fopna('out',29,0)

C     call snot(s,s,s,icast,ns,nr,nc)

      if (icast /= 0) then
        call ywrm(0,' ',icast,ifi,fmt,s,ofi,ns,nr,nc)
      else
        call awrit2('%% rows %i cols %i',' ',80,ifi,nr,nc)
        do  i = 1, nr
          write (ifi,'(22i7)') (s(i,j),j=1,nc)
        enddo
      endif

      call fclose(ifi)
      call awrit0(' prm: wrote '//trim(strn)//'. Continue?',' ',120,i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in prmx')
      return
      entry yprm0(fmt0)
      fmt = fmt0
      end

      subroutine zprm(strn,icast,s,ns,nr,nc)
C- Print complex matrix to file out
      implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(20) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      character fmt0*(*)
      save fmt
C     data fmt /'(9f15.10)'/
C      data fmt /'(9f18.11)'/
      data fmt /'(4f20.15)'/

      outs = ' '
      if (icast == 1)  outs = ' real'   ! print real part only
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

      return

      entry zprm0(fmt0)
      fmt = fmt0
      end

      subroutine zprm3(sfmt,msga,s,n1,n2,n3)
      implicit none
      character *(*) sfmt
      integer n1,n2,n3,msga
      double precision s(2,n1,n2,n3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna
      fmt = '(9f15.10)'
C     fmt = '(1p9e20.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i complex',' ',80,ifi,n1*n2,n3)
      do  i = 1, n1
        do  j = 1, n2
          write (ifi,fmt) (s(1,i,j,k),k=1,n3)
        enddo
      enddo
      do  i = 1, n1
        do  j = 1, n2
          write (ifi,fmt) (s(2,i,j,k),k=1,n3)
        enddo
      enddo
      call fclose(ifi)
      i = 1; if (sfmt(i:i) == '$') i = 2
      outs = ' '
      call awrit1(sfmt(i:),outs,80,0,msga)
      if (i == 2) then
        call info0(1,0,0,outs)
        return
      endif
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in zprm3')
      end
      subroutine prm3(sfmt,msga,s,n1,n2,n3)
      implicit none
      character *(*) sfmt
      integer n1,n2,n3,msga
      double precision s(n1,n2,n3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,n1*n2,n3)
      do  i = 1, n1
        do  j = 1, n2
          write (ifi,fmt) (s(i,j,k),k=1,n3)
        enddo
      enddo
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in zprm3')
      end

      subroutine prm3x(sfmt,msga,s,k1,k2,k3,n1,n2,n3)
      implicit none
      character *(*) sfmt
      integer k1,k2,k3,n1,n2,n3,msga
      double precision s(k1,k2,k3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,n1*n2,n3)
      do  i = 1, n1
        do  j = 1, n2
          write (ifi,fmt) (s(i,j,k),k=1,n3)
        enddo
      enddo
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs == 'q') call rx0('quit in prm3x')
      end

      subroutine rdm3(s,n1,n2,n3)
      implicit none
      integer n1,n2,n3
      double precision s(n1,n2,n3)
      integer rdm,fopna,ifi,i,j,k
      double precision xx

      ifi = fopna('in',29,0)
      if (rdm(ifi,0,0,' ',xx,n1*n2,n3) /= 1)
     .     call rx('rdm3:  file mismatch')
      do  i = 1, n1
        do  j = 1, n2
          read (ifi,*) (s(i,j,k),k=1,n3)
        enddo
      enddo

      call fclose(ifi)

      end

      subroutine zrdm3(ifi,s,n1,n2,n3,k1,k2,k3)
C- Reads from file ifi an array on a 3D mesh, written in 2D form
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   n1..n3:number of elements along each of 3 dimensions
Ci   k1..k3:dimenions s
Co Outputs
Co   s     :array s(1:n1,1:n2,1:n3)
Cl Local variables
Cl         :
Cr Remarks
Cr  Array is returned complex, but file data may be real
Cu Updates
Cu   18 Aug 17 Redsigned for general read
C ----------------------------------------------------------------------
      implicit none
      integer ifi,n1,n2,n3,k1,k2,k3,ierr
      double precision s(2,k1,k2,k3)
      integer i,j,k
      procedure(integer) :: rdm

      ierr = rdm(ifi,10,0,' ',s,n1*n2,n3)
      if (ierr < 0) call rx('zrdm3 failed to read file')

      do  i = 1, n1
      do  j = 1, n2
      read(ifi,*) (s(1,i,j,k), k=1,n3)
      enddo
      enddo

      if (ierr == 1) goto 99

      do  i = 1, n1
      do  j = 1, n2
      read(ifi,*,end=99,err=99) (s(2,i,j,k), k=1,n3)
      enddo
      enddo
      return

   99 continue
      call info0(10,0,0,'# zrdm3 (warning) reading real array')
      forall(i=1:n1,j=1:n2,k=1:n3) s(2,i,j,k) = 0d0
      end

      subroutine zprmx(strn,icast,s,ns,nr,nc)
      implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(11) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      fmt = '(9f15.10)'
      fmt = '(1p9e20.11)'
      fmt = '(4f20.15)'
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
