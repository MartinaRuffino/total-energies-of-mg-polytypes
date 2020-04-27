      subroutine gfsbar(nbas,nl,indxsh,lmin,lhdim,sbar)
C- Make average
      implicit none
      integer nbas,nl,lmin,lhdim,indxsh(*)
      double precision sbar(lhdim)
      integer n,lmr,ibas,l,ip,m
      double precision savgm,dsum

      lmr = 0
      do  3  ibas = 1, nbas
        do  2  l = 0, nl-1
          ip = indxsh(lmr+1)
          n = 2*l+1
          if (ip > lmin .and. ip <= lhdim) then
            savgm = dsum(n,sbar(ip),1) / n
            ip = ip+l
            do  1  m = -l, l
    1       sbar(ip+m) = savgm
          endif
          lmr = lmr + n
    2   continue
    3 continue

C      do  10  i = 1, lhdim
C        print *, i,sbar(i)
C   10 continue
      end
