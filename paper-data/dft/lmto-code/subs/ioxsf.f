      subroutine ioxsf(opt,s_lat,s_site,s_spec)
C- File I/O of structure in xsf format
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: ioxsf2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: z
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: ioxsf2
Ci Inputs
Ci   opt   :1s digit
Ci         :0 for read, 1 for write
Co Outputs
Co: nothing
Cr Remarks
Cr   writes the primitive vectors and positions in Angstrom units
Cr   on a file xsf.ext  in xsf format that can be read by xcrsyden among others
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   January 7 2010 (Walter Lambrecht) first created
C ---------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: z(:)
      real(8), allocatable :: pos(:)
C ... Local parameters
      double precision alat,plat(3,3)
      integer i1,i2,nbas,ifi,nglob
      integer fopn
      character*10 fn

      if (opt /= 1) call rx('ioxsf not ready for opt ne 1')
      nbas = nglob('nbas')
      alat = s_lat%alat
      plat = s_lat%plat
      allocate(z(nbas))
      allocate(pos(nbas*3))
      fn = 'xsf'
      ifi = fopn(fn)
      rewind ifi
      write(ifi,'("CRYSTAL")')
      write(ifi,'("PRIMVEC")')
      write(ifi,'(3f15.10)')
     .  ((plat(i1,i2)*alat*0.529177208d0,i1=1,3),i2=1,3)
      write(ifi,'("PRIMCOORD")')
      write(ifi,'(2i5)') nbas,1
      call ioxsf2(ifi,nbas,z,pos,s_spec,s_site,alat)
      call fclose(ifi)
      deallocate(z,pos)
      end

      subroutine ioxsf2(ifi,nbas,z,pos,s_spec,s_site,alat)
C- Kernel used by ioxsf: reads/writes atomic numbers and site positions
C---------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   nbas  :size of basis
Ci   alat  :length scale of lattice and basis vectors, a.u.
Cio Inputs/Outputs
Cio   z    :nuclear charge
Cio  pos   :basis vectors in units of alat
Cio  sspec :struct for species-specific information; see routine uspec
Cio    Elts read (written): z
Cio  ssite :struct for site-specific information; see routine usite
Cio    Elts read (written): spec pos
Cr Remarks
Cr   Site positions and atomic numbers are read from or written to disk,
Cr   in xsf style; see http://www.xcrysden.org/doc/XSF.html
Cr   Note: xsf format takes Cartesian positions in Angstrom.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,ifi
      double precision z(nbas),pos(3,nbas),alat
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer i,i2,is

      do  i = 1, nbas
        is = s_site(i)%spec
        z(i) = s_spec(is)%z
        pos(1:3,i) = s_site(i)%pos
        write(ifi,'(i4,2x,3f15.10)') int(z(i)),
     .    (pos(i2,i)*alat*0.529177208d0,i2=1,3)
      enddo
      end
