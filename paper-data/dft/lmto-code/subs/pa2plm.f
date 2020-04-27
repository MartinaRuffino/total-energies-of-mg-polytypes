      subroutine pa2plm(pa,nbas,lpdim,ldim,indxsh,iopt,p)
C- Expand l- independent matrix into matrix l- or lm- components
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpdim  dimension of p
Ci   nbas
Ci   iopt   1s digit
Ci          0  one-dimensional matrix
Ci          1  two-dimensional matrix
Ci          10s digit
Ci          0  p has dimension nbas
Ci          1  p has dimension nbas*nl
Ci          2  p has dimension nbas*nlm
Ci          100s digit
Ci          0  copy pa(i,j) to each element il,jl in p(i,il,j,jl)
Ci          1  copy pa(i,j) to diagonal elements il,jl in p(i,il,j,jl)
Ci          1000s digit
Ci          1  do not initialize pa to zero before starting
Co Outputs
Co   p:  the l-indepent parts of pa copied into l- (lm-) dependent p
Cr Remarks
Cr   pa2plm performs the inverse function to plm2pa
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lpdim,nbas,ldim,indxsh(*),iopt
      double precision pa(nbas,nbas),p(lpdim,lpdim)
C ... Local parameters
      integer ib,ilm,jb,jlm,nli,nlj,norbi,norbj,ntorbi,ntorbj,nd,
     .  opt1,opt2,opt3,offj,offi
      integer, parameter :: n0=10, nkap0=4
      integer ltabi(n0*nkap0),ltabj(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
C     double precision scale

      nd = 1
      if (mod(iopt,10) == 1) nd = 2
      opt1 = mod(iopt/10,10)
      opt2 = mod(iopt/100,10)
      opt3 = mod(iopt/1000,10)

C ... Simple copy if no indices to contract over
      if (opt1 == 0) then
        call dcopy(nbas**nd,pa,1,p,1)
        return
      endif

      if (opt3 == 0) call dpzero(p, lpdim**nd)

      nlj = 0
      do  jb = 1, nbas
C     Uses norbj,ntorbj
      call orbl(jb,0,ldim,indxsh,norbj,ltabj,ktab,offj,offl,ntorbj)
      if (opt1 == 1) ntorbj = norbj
        do  jlm = 1, ntorbj
        nlj = nlj+1

        if (nd == 1) then
          p(nlj,1) = pa(jb,1)
        else
          nli = 0
          do  ib = 1, nbas
C         Uses norbi,ntorbi
          call orbl(ib,0,ldim,indxsh,norbi,ltabi,ktab,offi,offl,ntorbi)
          if (opt1 == 1) ntorbi = norbi
          do  ilm = 1, ntorbi
            nli = nli+1
            if (opt2 == 0 .or. ilm == jlm) then
              p(nli,nlj) = pa(ib,jb)
            else
              p(nli,nlj) = 0
            endif
          enddo
          enddo
        endif
        enddo
      enddo

      end
