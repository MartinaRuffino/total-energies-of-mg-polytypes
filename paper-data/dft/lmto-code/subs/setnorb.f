      subroutine setnorb(s_site,offH,nbas)
      use structures
      implicit none
C ... Passed variables
      integer nkap0,n0H,nbas
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib

      do  ib = 1, nbas
        s_site(ib)%norb = offH(5,1,ib+1)-offH(5,1,ib)
      enddo

      end
