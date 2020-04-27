      subroutine wdiff(lpdim,nkp,wscrk,wscrkm,w0,w1,w2,w0m,w1m,w2m)
C- Difference w(P0) and w(P0-local P0)
      implicit none
      integer lpdim,nkp,i,j,k
      double precision wscrk(lpdim,lpdim,*),wscrkm(lpdim,lpdim,*),
     .  w0(lpdim,lpdim), w1(lpdim,lpdim), w2(lpdim,lpdim),
     .  w0m(lpdim,lpdim),w1m(lpdim,lpdim),w2m(lpdim,lpdim)
      do  10  i = 1, lpdim
      do  10  j = 1, lpdim
        w0(j,i) = w0(j,i) - w0m(j,i)
        w1(j,i) = w1(j,i) - w1m(j,i)
        w2(j,i) = w2(j,i) - w2m(j,i)
   10 continue

      do  20  k = 1, 2*nkp
      do  20  i = 1, lpdim
      do  20  j = 1, lpdim
   20 wscrk(j,i,k) = wscrk(j,i,k) - wscrkm(j,i,k)

      end

