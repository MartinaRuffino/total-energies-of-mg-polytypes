      subroutine tbaddh(nl,nbas,it,nitmax,npr,delta,hrs)
C-
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nl,nbas,it,nitmax,npr(0:1,nbas)
      double precision delta(nl**2,nl**2,nbas,nitmax+1,2),
     .                 hrs(nl**2,nl**2,1)

C Local Variables
      integer ib,j,ilm,ilmp,iprint,stdo,nglob

      stdo = nglob('stdo')
      if (iprint() > 50) write(stdo,1)
    1 format('TBADDH: site  ilm  H_in  Delta_last  Delta_now  H_out')
      do  ib = 1, nbas
        j = npr(1,ib) + 1
        do  ilm = 1, nl**2
          do  ilmp = 1, nl**2
            if (ilm == ilmp .and. iprint() > 50) write(stdo,2)
     .        ib,ilm,hrs(ilm,ilmp,j),delta(ilm,ilmp,ib,it,1),
     .        delta(ilm,ilmp,ib,it+1,1),hrs(ilm,ilmp,j)-
     .        delta(ilm,ilmp,ib,it,1)+delta(ilm,ilmp,ib,it+1,1)
            hrs(ilm,ilmp,j) = hrs(ilm,ilmp,j)
     .                      - delta(ilm,ilmp,ib,it,1)
     .                      + delta(ilm,ilmp,ib,it+1,1)
          enddo
        enddo
      enddo
    2 format (7x,2i4,4f10.4)
      end
