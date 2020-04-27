      subroutine GWdimpar(opt,s_site,s_spec,nbas,nl,nsp,maxcor,konfa)
C- Extract dimensioning parameters needed for GW
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :controls what parameters are returned
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   konfa :P.Q.N for partial waves
Co   maxcor:returned if 1s digit opt is nonzero
Co         :maxcor(1) = largest number of core radial functions in system
Co         :maxcor(2) = largest l for which a core state exists
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Dec 17 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,nbas,nl,nsp,maxcor(2),konfa(0:nl-1,nbas)
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer, parameter :: n0=10
      integer i,is,isp,l,iat
      real(8) pnu(n0,2),pnz(n0,2)

      if (mod(opt,10) == 1) then
        maxcor(1) = 0; maxcor(2) = -1
        iat = 0
        do  i = 1, nbas
          is = s_site(i)%spec
          if (s_spec(is)%lmxa < 0) cycle
          iat = iat+1
          pnu = s_site(i)%pnu
          pnz = s_site(i)%pz
          do  isp = 1, nsp
            do  l  = 0, s_spec(is)%lmxa
              konfa(l,iat) = pnu(l+1,isp)
              if (mod(pnz(l+1,isp),10d0) < pnu(l+1,isp) .and. pnz(l+1,isp) > 0)
     .          konfa(l,iat) = mod(pnz(l+1,isp),10d0)
              maxcor(1) = max(maxcor(1),konfa(l,iat)-l-1)
              if (konfa(l,iat)-l-1 > 0) maxcor(2) = max(maxcor(2),l)
            enddo
          enddo
        enddo
      endif

      end
