      subroutine getidu(nbas,s_spec,s_site,idu)
C- Get idu for all sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  idu
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Co Outputs
Co   idu   :idu(l+1,ib)=1 => this l and site has a nonlocal U matrix
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   Mar 29, 07 (Jialei) first created
C------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,idu(4,nbas)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer is,ib

      call iinit(idu,4*nbas)
      do  ib = 1, nbas
        is = s_site(ib)%spec
        idu(1:4,ib) = s_spec(is)%idu
      enddo

      end



