      subroutine mstrgg(e,dr,nlma,nlmb,nla,nlb,cg,indxcg,jcg,cy,gs,gd)
C- Gradient of mol strux and energy derivatives
C  Expansion theorem is:  H(k,r-dr) = sum_m S(m,k,dr)*J(m,r)
C  Alternatively, let dr=r2-r1 and H and J be centered at r2,r1.
C  Expansion theorem is:  H(k,r-r2) = sum_m S(m,k,r2-r1)*J(m,r-r1)
      implicit none
      integer nlma,nlmb,nla,nlb
      double precision gs(nla,nlb,3),gd(nla,nlb,3),cy(1),cg(1),dr(3)
      integer indxcg(1),jcg(1)
      double precision hl(144),hd(144),
     .  efac(0:20),sig(0:20),ghl(144,3),ghd(144,3)
      integer icg,icg1,icg2,ii,indx,ipow,klm,l,lk,llm,lm,
     .  lmax,m,mlm,nlm,ndh,ll,ixyz
      double precision e,fpi,fac,sum,sumg,sudg,ube

      if (nlma <= 0 .or. nlmb <= 0) return
      if (nlma > nla .or. nlmb > nlb) call rx('mstrgg: nlma gt ndim')

C --- Make solid Hankel functions HL ---
      lmax = ll(nlma)+ll(nlmb)
      ndh=144
      call solhpg(e,dr,lmax,ndh,hl,ghl,hd,ghd,cy)

      fpi = 16*datan(1d0)
      nlm = (lmax+1)**2
      efac(0) = 1d0
      sig(0) = 1d0
      do  30  l = 1, lmax
      efac(l) = -e*efac(l-1)
   30 sig(l) = -sig(l-1)

C ---   Make strx ---
        ube = 1d0/e
        do 20 ixyz = 1, 3
        do  40  mlm = 1, nlma
          lm = ll(mlm)
          do  40  klm = 1, nlmb
          lk = ll(klm)
          sumg = 0d0
          sudg = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  45  icg = icg1, icg2
            llm = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            fac = cg(icg)*efac(ipow)
            sumg = sumg+ fac*ghl(llm,ixyz)
            sudg = sudg+ fac*(ghd(llm,ixyz)+ ipow*ube*ghl(llm,ixyz))
   45     continue
          gs(mlm,klm,ixyz) = fpi*sumg*sig(lk)
          gd(mlm,klm,ixyz) = fpi*sudg*sig(lk)
   40   continue
   20 continue
      end
