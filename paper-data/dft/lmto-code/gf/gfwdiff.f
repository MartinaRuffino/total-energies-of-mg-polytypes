      subroutine gfwdiff(lpdim,nkp,wscrk,wscrl)
C- Difference w(P0) and w(local P0)
Cu Updates
Cu 08 Oct 04(T. Sandu) Redone for nsp=2(spin-polarized)
Cu-Made 10 May 04 (T.Sandu)
      implicit none
      integer lpdim,nkp,i,k
      double precision wscrk(lpdim,lpdim,2,nkp),
     . wscrl(lpdim,lpdim,2)
      print *,'gfwdiff at start: merge into wdiff?'
      do  21  k = 1, nkp
      do  20  i = 1, lpdim
      wscrk(i,i,1,k) = wscrk(i,i,1,k) - wscrl(i,i,1)
   20 wscrk(i,i,2,k) = wscrk(i,i,2,k) - wscrl(i,i,2)
c      print *,'iq',k
c      call yprm('wdiff-q',2,wscrk(1,1,1,k),lpdim*lpdim,lpdim,lpdim,lpdim)
 21   continue
      end

