      subroutine setofl(mode,s_site,s_spec,nbas,nvl,iiv0)
C- Sets vector of offsets for lmxl
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :not used now
Ci   nbas  :size of basis
Co Outputs
Co    nvl  :dimension of array containing charge channels
Co    iiv0 :offsets to start of channels for sites 1..nbas
Cr Remarks
Cr   Useful for, e.g. multiple threaded code.
Cu Updates
Cu   1 May 00 Adapted from nfp setofl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nvl,iiv0(nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,iv0,lmxl,nlm

      if (mode /= 0) call rx('setofl: bad mode')
      iv0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
      enddo
      nvl = iv0
      end
