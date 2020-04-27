      subroutine pkgdlm(nbas,s_site,gibbs,shfac)
C- Package Gibbs weights for CPA in s_site
C ----------------------------------------------------------------------
Cu Updates
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   01 Nov 11 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas
      double precision gibbs(*),shfac(3,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ioff,ncomp(nbas),dlmcl(nbas)
      double precision xx

      call sitepack(s_site,1,nbas,'ncomp',1,ncomp,xx)
      call sitepack(s_site,1,nbas,'dlmcl',1,dlmcl,xx)
C     call sitepack(s_site,1,nbas,'ocpawt',1,ocpawt,xx)

      do  ib = 1, nbas
        if (ncomp(ib) < 2) cycle
        ioff = 3 * (dlmcl(ib) - 1) + 1
        call dpscop(gibbs,s_site(ib)%cpawt,ncomp(ib),dlmcl(ib),1,1d0)
        call dpscop(shfac,s_site(ib)%bxc,3*ncomp(ib),ioff,1,1d0)
      enddo

      end
