      subroutine suldau(nbas,s_spec,s_site,nlibu,lmaxu,lldau)
C- Finds lda+U sites and counts number of blocks
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa idu
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nbas  :size of basis
Co Outputs
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Co   nlibu :number of LDA+U blocks
Co   lmaxu :highest l for which a U block found, used as
Co         :dimensioning parameter for U matrix
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   27 Apr 05 (Lambrecht) first created
C------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nlibu,lmaxu,lldau(nbas),is,ib,l,lmxa,idu(4)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
C     integer stdo,nglob
C     stdo = nglob('stdo')

      nlibu = 0
      lmaxu = 0
      do  ib = 1, nbas
        lldau(ib) = 0
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        idu = mod(s_spec(is)%idu,100)
        do  l = 0, min(lmxa,3)
          if (idu(l+1) /= 0) then
            if (lldau(ib) == 0) lldau(ib) = nlibu+1
            nlibu = nlibu+1
            lmaxu = max(lmaxu,l)
          endif
        enddo
      enddo

      if (nlibu /= 0) then
        call info2(10,1,0,' suldau:  %i U block(s)  lmaxu = %i',nlibu,lmaxu)
      endif

      end
