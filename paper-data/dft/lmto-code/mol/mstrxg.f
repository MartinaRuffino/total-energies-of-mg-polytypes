      subroutine mstrxg(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,
     .  job,s,sd,gs,gd)
C- Gradient of mol strux and energy derivatives
C  job=0: calculate s,gs only; otherwise calculate sd and gd also
      implicit none
      integer nlma,nlmb,ndim
      double precision    s(ndim,nlmb),gs(ndim,nlmb,3),
     .  cy(1),cg(1),dr(3),sd(ndim,nlmb),gd(ndim,nlmb,3)
      integer indxcg(1),jcg(1),job
      double precision hl(144),hd(144),psi(0:20),phi(0:20),
     .  efac(0:20),sig(0:20),ghl(144),ghd(144)
      integer icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,llm,lm,
     .  lmax,m,mlm,nlm,ixyz,ll
      double precision e,fpi,psidot,r2,fac,sum,sud,sumg,sudg,ube

      if (nlma > ndim) call rx('mstrxg: nlma gt ndim')
      if (nlma < 0 .or. nlmb < 0) return

C --- Make solid Hankel functions HL ---
      lmax = ll(nlma)+ll(nlmb)
      if (lmax+1 > 9) call rx('mstrxg: increase lmax')
      call sylm(dr,hl,lmax+2,r2)
      call bessl(e*r2,lmax+2,phi,psi)
      ilm = 0
      fac = dsqrt(r2)
      do  10  l = 0, lmax+1
        fac = fac/r2
        psidot = ((l+l+1)*psi(l) - psi(l+1))/(e+e)
        do  10  m = -l, l
        ilm = ilm+1
        hd(ilm) = fac*psidot*cy(ilm)*hl(ilm)
        hl(ilm) = fac*psi(l)*cy(ilm)*hl(ilm)
   10 continue

      fpi = 16*datan(1d0)
      nlm = (lmax+1)**2
      efac(0) = 1d0
      sig(0) = 1d0
      do  30  l = 1, lmax
      efac(l) = -e*efac(l-1)
   30 sig(l) = -sig(l-1)
C ---   Make gradient of HL ---
      do  20  ixyz = 1, 3
        ilm = ixyz
        if (ixyz == 1) ilm=4
        fac = -dsqrt(fpi/3)
        do  22  mlm =1, nlm
          lm=ll(mlm)
          ii = max0(mlm,ilm)
          ii = (ii*(ii-1))/2 + min(mlm,ilm)
          icg1 = indxcg(ii)
          icg2 = indxcg(ii+1)-1
          sum = 0d0
          sud = 0d0
          do  25  icg = icg1, icg2
            klm = jcg(icg)
            lk = ll(klm)
            ipow = (lm+1-lk)/2
            if (ipow == 0) then
              sum = sum + cg(icg)*hl(klm)
              sud = sud + cg(icg)*hd(klm)
            else
              sum = sum - cg(icg)*e*hl(klm)
              sud = sud - cg(icg)*(e*hd(klm)+hl(klm))
            endif
   25     continue
          ghl(mlm) = fac*sum
          ghd(mlm) = fac*sud
   22   continue

C ---   Make strx ---
        ube = 1d0/e
        do  40  mlm = 1, nlma
          lm = ll(mlm)
          do  40  klm = 1, nlmb
          lk = ll(klm)
          sum  = 0d0
          sumg = 0d0
          sud  = 0d0
          sudg = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  45  icg = icg1, icg2
            llm = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            fac = cg(icg)*efac(ipow)
            sum  = sum + fac*hl(llm)
            sumg = sumg+ fac*ghl(llm)
            sud  = sud + fac*(hd(llm) + ipow*ube*hl(llm))
            sudg = sudg+ fac*(ghd(llm)+ ipow*ube*ghl(llm))
   45     continue
          s(mlm,klm) = fpi*sum*sig(lk)
          gs(mlm,klm,ixyz) = fpi*sumg*sig(lk)
          if (job /= 0) then
            sd(mlm,klm)= fpi*sud*sig(lk)
            gd(mlm,klm,ixyz) = fpi*sudg*sig(lk)
          endif
   40   continue
   20 continue
      end
C#ifdefC TEST
CC Test program to check mstrxg
C      subroutine fmain
C      implicit none
CC ... For lmx=10
C      double precision cg(62152)
C      integer indxcg(7382),jcg(62152)
CC ... For lmx=8
CC      double precision cg(22662)
Cc      integer indxcg(3322),jcg(22662)
CC ... For lmx=6
Cc      double precision cg(6408)
Cc      integer indxcg(1226),jcg(6408)
C      double precision dr(3),e,hl(200),ghl(600),ghd(600),cy(16**2),
C     .  s(100**2),gs(100**2*3),gd(100**2*3),sp(100**2),sm(100**2),
C     .  sd(100**2),spd(100**2),smd(100**2),hd(200)
C      integer lmxa,lmxb,nlma,nlmb
C      common /static/ indxcg,cg,jcg,cy,s,sd,sp,sm,spd,smd,gs,gd
C
C      call scg(10,cg,indxcg,jcg)
C      call sylmnc(cy,16)
C
C   99 print *, 'lmxa,lmxb,e,dr='
C      read(*,*) lmxa,lmxb,e,dr
C
C      nlma=(lmxa+1)**2
C      nlmb=(lmxb+1)**2
Cc      call xxx(dr,e,lmxa,cy,cg,indxcg,jcg,hl,hd,ghl,ghd)
C      call yyy(dr,e,cy,cg,indxcg,jcg,nlma,nlmb,
C     .  s,sp,sm,sd,spd,smd,gs,gd)
Cc      goto 99
C
C      end
C      subroutine yyy(dr,e,cy,cg,indxcg,jcg,nlma,nlmb,
C     .  s,sp,sm,sd,spd,smd,gs,gd)
C      implicit none
C      integer indxcg(1),jcg(1),nlma,nlmb
C      double precision cg(1),cy(1),dr(3),e,dx,gd(nlma,nlmb,3),
C     .  s(nlma,nlmb),sp(nlma,nlmb),sm(nlma,nlmb),gs(nlma,nlmb,3),
C     .  sd(nlma,nlmb),spd(nlma,nlmb),smd(nlma,nlmb),xx,tops,topd
C      integer ndim,i1,i2,i
C
C      tops = 0d0
C      topd = 0d0
C      ndim = nlma
C      dx = 1d-6
C      do  10  i = 1, 3
C        dr(i)=dr(i)+dx
C        call mstrxg(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,1,sp,spd,gs,gd)
C        dr(i)=dr(i)-2*dx
C        call mstrxg(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,1,sm,smd,gs,gd)
C        dr(i)=dr(i)+dx
C        call mstrxg(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,1,s,sd,gs,gd)
C
C        do  12  i1 = 1, nlma
C   12   print 333, ((sp(i1,i2)-sm(i1,i2))/(2*dx),i2=1,nlmb)
C  333   format(9f8.5)
C        print *,
C        do  14  i1 = 1, nlma
C        do  15  i2 = 1, nlmb
C   15   tops = max(tops,dabs(gs(i1,i2,i)-(sp(i1,i2)-sm(i1,i2))/(2*dx)))
C   14   print 333, (gs(i1,i2,i)-(sp(i1,i2)-sm(i1,i2))/(2*dx),i2=1,nlmb)
C        print *, '--- For dot ---'
C        do  16  i1 = 1, nlma
C   16   print 333, ((spd(i1,i2)-smd(i1,i2))/(2*dx),i2=1,nlmb)
C        print *,
C        do  18  i1 = 1, nlma
C        do  17  i2 = 1, nlmb
C   17   topd =max(topd,dabs(gd(i1,i2,i)-(spd(i1,i2)-smd(i1,i2))/(2*dx)))
C   18   print 333, (gd(i1,i2,i)-(spd(i1,i2)-smd(i1,i2))/(2*dx),
C     .    i2=1,nlmb)
C        print *, '-----------'
C   10 continue
C
C      print 335, tops, topd
C  335 format(' max errors for grad s, grad sdot:',2f12.6)
C
C      end
C#endif
