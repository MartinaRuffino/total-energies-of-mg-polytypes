      subroutine pnu2pl(pnu,lmaxa,n0,nsp,ib,pl)
C- copies pnu dimensioned (n0,nsp,ib) to pl dimensioned (9,2)
      implicit none
      integer n0,nsp,ib,lmaxa
      double precision pnu(n0,nsp,ib),pl(9,2)

      call dpzero(pl,9*nsp)
      call dpcopy(pnu(1,1,ib),pl,1,1+lmaxa,1d0)
      if (nsp == 2) call dpcopy(pnu(1,2,ib),pl(1,2),1,1+lmaxa,1d0)
      end
