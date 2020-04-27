      subroutine mstrpg(e,dr,nlma,nlmb,ndim,mdim,cg,indxcg,jcg,cy,
     .   s,sd,gs,gd)
C- Gradient of mol strux and energy derivatives
C      implicit none
      double precision s(ndim,nlmb),gs(ndim,mdim,3),
     .  cy(1),cg(1),dr(3),sd(ndim,nlmb),gd(ndim,mdim,3)

      integer indxcg(1),jcg(1)
      double precision hl(200),hd(200),
     .  efac(0:20),sig(0:20),ghl(200,3),ghd(200,3)
      integer icg,icg1,icg2,ii,indx,ipow,klm,l,lk,llm,lm,
     .  lmax,mlm,ndim,mdim,nlma,nlmb,nlm,ll,ndh
      double precision e,fpi,fac,sum,sud,ube,
     .  smgx,smgy,smgz,sdgx,sdgy,sdgz

      if (nlma > ndim) call rx('mstrpg: nlma gt ndim')
      if (nlma < 0 .or. nlmb < 0) return

C --- Make solid Hankel functions HL ---
      lmax = ll(nlma)+ll(nlmb)
      ndh=200
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
        do  40  mlm = 1, nlma
          lm = ll(mlm)
          do  40  klm = 1, nlmb
          lk = ll(klm)
          sum  = 0d0
          smgx = 0d0
          smgy = 0d0
          smgz = 0d0
          sud  = 0d0
          sdgx = 0d0
          sdgy = 0d0
          sdgz = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  45  icg = icg1, icg2
            llm = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            fac = cg(icg)*efac(ipow)
            sum  = sum + fac*hl(llm)
            smgx = smgx+ fac*ghl(llm,1)
            smgy = smgy+ fac*ghl(llm,2)
            smgz = smgz+ fac*ghl(llm,3)
            sud  = sud + fac*(hd(llm) + ipow*ube*hl(llm))
            sdgx = sdgx+ fac*(ghd(llm,1)+ ipow*ube*ghl(llm,1))
            sdgy = sdgy+ fac*(ghd(llm,2)+ ipow*ube*ghl(llm,2))
            sdgz = sdgz+ fac*(ghd(llm,3)+ ipow*ube*ghl(llm,3))
   45     continue
          s(mlm,klm) = fpi*sum*sig(lk)
          gs(mlm,klm,1) = fpi*smgx*sig(lk)
          gs(mlm,klm,2) = fpi*smgy*sig(lk)
          gs(mlm,klm,3) = fpi*smgz*sig(lk)
          sd(mlm,klm)= fpi*sud*sig(lk)
          gd(mlm,klm,1) = fpi*sdgx*sig(lk)
          gd(mlm,klm,2) = fpi*sdgy*sig(lk)
          gd(mlm,klm,3) = fpi*sdgz*sig(lk)
   40   continue
      end
