      subroutine dlmentrp(s_site,nbas,s)
C- Entropy as -<log(P)>
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp cpawt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Co Outputs
Co   s     :entropy
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Apr 12
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas
      double precision s
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer i,ib,nthet
      double precision wt

      s = 0
      do  ib = 1, nbas
        nthet = s_site(ib)%ncomp
        if (nthet < 2) cycle
        do  i = 1, nthet
          wt = s_site(ib)%cpawt(i)
          s = s - wt * dlog(wt)
        enddo
      enddo

      end
