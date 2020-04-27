      subroutine mstrug(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,s,gs)
C- Molecular structure constants and gradient
C      implicit none
      double precision s(ndim,nlmb),gs(ndim,nlmb,3),
     .  cy(1),cg(1),dr(3)

      integer indxcg(1),jcg(1),lx
      parameter (lx=14)
      double precision hl((lx+2)**2),efac(0:20),sig(0:20),
     .  ghl((lx+1)**2,3)
      integer icg,icg1,icg2,ii,indx,ipow,klm,l,lk,llm,lm,
     .  lmax,mlm,ndim,nlma,nlmb,nlm,ndh,ll
      double precision e,fpi,fac,sum,sumx,sumy,sumz

      if (nlma > ndim) call rx('mstrug: nlma gt ndim')
      if (nlma < 0 .or. nlmb < 0) return
C     if (lmax+1 > lx) call rx('mstrug: increase lx')

C --- Make solid Hankel functions HL ---
      lmax = ll(nlma)+ll(nlmb)
      ndh = (lx+1)**2
      call solhg(e,dr,lmax,ndh,hl,ghl,cy)

      fpi = 16*datan(1d0)
      nlm = (lmax+1)**2
      efac(0) = 1d0
      sig(0) = 1d0
      do  30  l = 1, lmax
      efac(l) = -e*efac(l-1)
   30 sig(l) = -sig(l-1)

C ---   Make strx ---
        do  40  mlm = 1, nlma
          lm = ll(mlm)
          do  40  klm = 1, nlmb
          lk = ll(klm)
          sum  = 0d0
          sumx = 0d0
          sumy = 0d0
          sumz = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  45  icg = icg1, icg2
            llm = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            fac = cg(icg)*efac(ipow)
            sum  = sum + fac*hl(llm)
            sumx = sumx+ fac*ghl(llm,1)
            sumy = sumy+ fac*ghl(llm,2)
            sumz = sumz+ fac*ghl(llm,3)
   45     continue
          s(mlm,klm) = fpi*sum*sig(lk)
          gs(mlm,klm,1) = fpi*sumx*sig(lk)
          gs(mlm,klm,2) = fpi*sumy*sig(lk)
          gs(mlm,klm,3) = fpi*sumz*sig(lk)
   40   continue
      end
