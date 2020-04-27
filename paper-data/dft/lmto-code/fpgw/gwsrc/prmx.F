C#define NOAWRITE
      subroutine prmx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(14) fmt, fmt0, strn*(*), outs*80
      integer i,j,fopna
      procedure(integer) :: iopen,iclose,i1mach
      save fmt
      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
C     ifi = fopna('out',29,0)
      ifi = iopen('out',1,-1,0)
      rewind ifi
C#ifdef NOAWRITE
      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#elseC
C      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(i,j,1),j=1,nc)
      enddo
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      ifi = iclose('out')

C#ifdef NOAWRITE
      outs = ' prm: wrote '//strn//' continue?'
      print *, outs
C#elseC
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdef NOAWRITE
      if (outs == 'q') stop 'quit in prmx'
C#elseC
C      if (outs == 'q') call rx0('quit in prmx')
C#endif
      return

      entry prmx0(fmt0)
      fmt = fmt0
      end

      subroutine prrmsh(strn,rofi,s,ns,nr,nc)
C- Print YL expansion to file
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn  :printout string
Ci   rofi  :radial mesh points
Ci   s     :YL expansion to print
Ci   ns    :leading dimension of s
Ci   nr    :number of radial mesh points
Ci   nc    :number of columns
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   The following integrates the L=0 functions
Cr      mc out -coll 1,2 -av:nr,1 rmax -int 0 rmax
Cr   The following integrates r^2* L=1 function
Cr      mc out -e2 x1 'x2*x1*x1' -av:nr,1 rmax -int 0 rmax
Cu Updates
Cu   11 Nov 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nc,ns,ifi
      double precision rofi(nr),s(ns,nc)
C ... Local parameters
      character*(20) fmt, strn*(*), outs*80
      integer i,j,i1
      procedure(integer) :: iopen,iclose,i1mach
C     fmt = '(f10.6,10f15.10)'
      fmt = '(f11.8,1p,10g18.10)'
C     ifi = fopna('out',-1,0)

      ifi = iopen('out',1,0,0)
      i1 = 1
      if (s(1,1) == 0) i1 = 2
C#ifdef NOAWRITE
      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr-i1+1,nc+1
C#elseC
C     call awrit2('%% rows %i cols %i',' ',80,ifi,nr-i1+1,nc+1)
C#endif

      do  i = i1, nr
        write(ifi,fmt) rofi(i), (s(i,j), j=1,nc)
      enddo
      ifi = iclose('out')

C#ifdef NOAWRITE
      outs = 'prrmsh: done writing to file out data '//trim(strn)//'.  Continue?'
      print *, outs
C#elseC
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdef NOAWRITE
      if (outs == 'q') call rx0('quit in prrmsh')
C#elseC
C      if (outs == 'q') call rx0('quit in prrmsh')
C#endif


      end
