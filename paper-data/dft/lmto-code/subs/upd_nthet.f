      subroutine upd_nthet(nspec,s_spec)
C- Update nthet in s_spec when nang is specified
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  ncpa nthet iscpa
Co     Stored:     ncomp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nspec :number of species
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nspec
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer is,nang(2),ntot,maxang
      parameter (maxang=10000)
      double precision x(3,maxang),w(maxang)

      do is = 1, nspec
        nang(1:2) = s_spec(is)%nang(1:2)
        if (nang(1) == 0) cycle
        call fpiint(nang(1),nang(2),ntot,x,w)
        if (ntot > maxang) call rx('Increase maxang in upd_nthet')
        s_spec(is)%nthet = ntot
      enddo

      end

