      subroutine symre(lov,trh,nl,nsp,nbas,ng,istab,rho,rholm,e)
C- Symmetrise e,rho,rholm
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lov, trh
Ci   e, rho, rholm (unsymmetrised, from tbfrce)
Co Outputs:
Co   e, rho, rholm (symmetrised)
Cr Remarks
Cr   Symmetrise rho, in the usual way.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nl,nsp,nbas,ng,istab(nbas,ng)
      double precision rho(0:nl-1,nsp,nbas),rwk(0:nl-1,nsp,nbas),
     .                 rholm(nl**2,nsp,nbas),rwklm(nl**2,nsp,nbas),
     .                 e(nbas,nsp),wk(nbas,nsp)
      logical trh,lov
C Local Variables
      integer isp,ib,jb,ig,l,ilm

      rwk = rho*(1d0/ng)
      rho = 0.0d0
      if (trh) then
        wk = e*(1d0/ng)
        e = 0.0d0
      endif
      if (lov) then
        rwklm = rholm*(1d0/ng)
        rholm = 0.0d0
      endif

      do   isp = 1, nsp
        do   ig = 1, ng
          do   ib = 1, nbas
            jb = istab(ib,ig)
            if (trh) e(ib,isp) = e(ib,isp) + wk(jb,isp)
            do  l = 0, nl-1
              rho(l,isp,ib) = rho(l,isp,ib) + rwk(l,isp,jb)
            enddo
            if (lov) then
              do  ilm = 1, nl*nl
                rholm(ilm,isp,ib) = rholm(ilm,isp,ib)
     .                            + rwklm(ilm,isp,jb)
              enddo
            endif
          enddo
        enddo
      enddo
      end
