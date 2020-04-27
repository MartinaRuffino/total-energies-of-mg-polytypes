      integer function iositp(lio,s_lat,s_spec,s_site,ifi)
C- File I/O of site data into s_spec and s_site strux
C ----------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec relax pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   lio  :controls file handling and header I/O
Ci         1s digit: file handling
Ci           0 for read
Ci           1 for write
Co Outputs
Co   ifi  :file logical unit for read/write
Cl Local variables
Cr Remarks
Cr  *File I/O in simple POSCAR format
Cu Updates
Cu   22 Apr 15  First created
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lio,ifi
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      integer,allocatable :: nrspec(:)
      character, allocatable:: slabl(:)*8
C ... Local parameters
      integer nspec,nbas,i1,i2,lio0,ib,is,i,k
      procedure(integer) :: nglob,fopng
      double precision alat,plat(3,3),plati(3,3),posx(3)
      character fnam*256
      integer, parameter :: procid=0, master=0
C ... External calls
      external dgemm,dinv33,rx,word

C ... Setup
      nspec = nglob('nspec')
      nbas = nglob('nbas')
      lio0   = mod(lio,10)
      allocate(nrspec(nspec),slabl(nspec))

C --- File write ---
      if (mod(lio0,2) == 1 .and. procid == master) then

        fnam = 'POSCAR'
        call word(fnam,1,i1,i2)
        ifi = fopng(fnam(i1:i2),-1,0)

C       Convert alat to AA from a.u.
        alat = s_lat%alat
        alat = alat * 0.529177d0
        write(ifi,"(f19.14)") alat

        plat = s_lat%plat
        write(ifi,"(1x,3f22.16)") plat ! POSCAR writes plat_ij in site file order

        do  i = 1, nspec
          slabl(i) = s_spec(i)%name
        enddo
        write(ifi,"(3x,100a5)") slabl(1:nspec)
        nrspec = 0
        do  ib = 1, nbas
          is = s_site(ib)%spec
          if (is > nspec) call rx('iositp: species out of range')
          nrspec(is) = nrspec(is) + 1
        enddo
        write(ifi,"(100i6)") nrspec(1:nspec)

        call dinv33(plat,0,plati,posx)
        write (ifi,"('Direct')")
        do  is = 1, nspec
          do  ib = 1, nbas
            if (is == s_site(ib)%spec) then
C           Forward: pos+ = plat posp+
C           call dgemm('N','N',3,1,3,1d0,plat,3,plx,3,0d0,posl,3)
C           Reverse:  posp+ = (plat)^-1 pos+
            call dgemm('N','N',3,1,3,1d0,plati,3,s_site(ib)%pos,3,0d0,posx,3)
            write(ifi,"(3f20.16,3x,3L2)") posx, (s_site(ib)%relax(k) /= 0, k=1,3)
           endif
          enddo
        enddo
        iositp = 0

C --- File read ---
      else
        call rx('reading data from poscar style is not implemented')
      endif

      end
