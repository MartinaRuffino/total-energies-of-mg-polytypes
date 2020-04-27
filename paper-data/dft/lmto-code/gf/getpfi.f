      subroutine getpfi(s_site,nl2,nsp,ib,icomp,ncomp,iprmb,lidim,lhdim,
     .  ldpf,lpfi,pfun,pfi)
C- Returns potential functions for lower block, one CPA (site,component)
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   nl2   :spacing between offsets to pfun in successive sites:
Ci          offset to pfun for site ib is nl2*(ib-1).
Ci          See makidx.f
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ib    :potential function
Ci   icomp :CPA component. pfi returned for components icomp,ncomp
Ci   ncomp :number of components to copy; see Outputs
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lidim :number of lower+intermediate orbitals
Ci   lhdim :number of lower+intermediate+higher orbitals
Ci   ldpf  :dimensions pfun, aka pfdim
Ci   lpfi  :dimensions pfi
Ci   pfun  :vector of potential functions (mkptfp.f,mkptfpd.f)
Co Outputs
Co   pfi   :potential functions copied from pfun for orbitals in
Co         :lower+intermediate block, i.e. iprmb()<=lidim
Co          CPA site?   icomp      ncomp     pfi returned in :
Co            no          0          1       pfi(:,1,1:nsp)
Co            yes         0         >0       pfi(:,1:ncomp,1:nsp)
Co            yes        >0          1       pfi(:,1,1:nsp) for icomp
Co            Nonsensical inputs:
Co            any         any       <=0
Co            any   >s_site(ib)%ncomp
Co            no          any       <>1
Co            no          <>0       any
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   02 Apr 13  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl2,nsp,ib,lidim,lhdim,icomp,ncomp,ldpf,lpfi,iprmb(*)
      double complex pfun(ldpf,nsp),pfi(lpfi,ncomp,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer lm,offpi,ilm,idlmoff,jb,ip,isp,jcomp

      if (ncomp < 1) goto 999 ! nonsensical input
      if (icomp > s_site(ib)%ncomp) goto 999 ! component doesn't exist

C ... Offset to pfun for 1st CPA component, or icomp^th component
      if (s_site(ib)%ncomp > 1) then
        idlmoff = lhdim + (icomp-1)*nl2
        do  jb = 1, ib-1
          idlmoff = idlmoff + nl2*s_site(jb)%ncomp
        enddo
      else
        if (icomp /= 0 .or. ncomp /= 1) goto 999
      endif

C ... Extract pfi for orbitals for which iprmb < lidim
C     Possible bug: idlmoff in CPA for downfolded orbitals
      offpi = nl2*(ib-1)
      jcomp = 1; if (ncomp > 1) jcomp = icomp
      do
        ilm = 0
        do  lm = 1, nl2
          ip = iprmb(offpi+lm)
          if (ip > lidim) cycle
          if (s_site(ib)%ncomp > 1) then
            idlmoff = idlmoff + 1
            ip = idlmoff
          endif
          ilm = ilm+1
C         print *, ib, ilm, ip
          do  isp = 1, nsp
            pfi(ilm,jcomp,isp) = pfun(ip,isp)
          enddo
        enddo
        jcomp = jcomp + 1
        if (ncomp == 1 .or. jcomp > s_site(ib)%ncomp) exit
      enddo

      return

C --- Error exit --
  999 call rx('getpfi: bad input')

      end
