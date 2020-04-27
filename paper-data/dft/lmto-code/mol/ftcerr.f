      subroutine ftcerr(d,mmax,wp,np,xp,zp,zc1,zc2,werr,wk,fit,
     .  nx,lx,xie,nlm,lmax,ph1,ph2,lp1,ep1,lp2,ep2,r1,r2,lsym,b,
     .  ndimx,nrhsx,err1,err2,errs)
C- Generate rms and max errors of the TCF.
C  MSM definition of errors.
      parameter (n0=10)
      implicit real*8 (a-h,p-z), integer (o)
      dimension irhs(0:20),lx(1),wp(1),wk(1),fit(1),
     .  xie(np,nlm,1),ph1(np,1),ph2(np,1),xp(1),zp(1),werr(1),
     .  arr1(0:n0,0:n0),arr2(0:n0,0:n0),cx(200),errs(0:lp1,0:lp2),
     .  b(ndimx,nrhsx,0:mmax),fac1(0:8),fac2(0:8),
     .  err1(nrhsx,0:mmax),err2(nrhsx,0:mmax)
      character*1 cl(0:8)
      data cl /'s','p','d','f','g','5','6','7','8'/
      call getpr(ipr)
      call sxlmnc(cx,10)
      if(ipr >= 30) write(6,455) d,ep1,ep2
  455 format(/' ftcerr:  d=',f9.4,'   ep1,ep2=',2f9.4)

C --- Make errvol and werr ----------------
      erralf=1.5d0
      errvol=0d0
      do 33 ip=1,np
      d1=dsqrt(xp(ip)**2+(zp(ip)-zc1)**2)
      d2=dsqrt(xp(ip)**2+(zp(ip)-zc2)**2)
      werr(ip)=0d00
      if(d1 < erralf*r1.or.d2 < erralf*r2) werr(ip)=2d0*wp(ip)
      if(d1 < r1.or.d2 < r2) werr(ip)=0d0
  33  errvol=errvol+werr(ip)
      if(ipr >= 35) write(6,331) errvol,erralf
  331 format(' rms-err:  vol=',f8.2,'   from radii *',f6.2)
      call hrmsav(r1,lp1,ep1,fac1)
      call hrmsav(r2,lp2,ep2,fac2)
      if(ipr >= 35) write(6,330) 1,(fac1(l),l=0,lp1)
      if(ipr >= 35) write(6,330) 2,(fac2(l),l=0,lp2)
  330 format(' rms avg of phi at r',i1,': ',6f9.4)

C --- For each (l1,m1,l2,m2) accumulate fit errors ---
      do 5 m=0,mmax
    5 irhs(m)=0
      do 20 l1=0,lp1
      do 20 m1=0,l1
      lmp1=(l1*(l1+1))/2+m1+1
      ltop=lp2
      if(lsym == 1)ltop=l1
      do 21 l2=0,ltop
      mtop=l2
      if (lsym == 1.and.l2 == l1) mtop=m1
      do 21 m2=0,mtop
      lmp2=(l2*(l2+1))/2+m2+1
      fac=1d0/(fac1(l1)*fac2(l2))
C ... Hold wave function product in work array
      do 22 ip=1,np
   22 wk(ip)=ph1(ip,lmp1)*ph2(ip,lmp2)
      mp=m1+m2
      mm=iabs(m1-m2)
C --- Fit and error for m1+m2 ---
      if(mp <= mmax) then
        irhs(mp)=irhs(mp)+1
        jrhs=irhs(mp)
        lm=((2*lmax-mp+1)*mp)/2 + 1
        call dpzero(fit,np)
        i=0
        do 25 ie=1,nx
        do 25 l=mp,lx(ie)
        i=i+1
        do 25 ip=1,np
   25   fit(ip)=fit(ip) + xie(ip,lm+l,ie)*b(i,jrhs,mp)
        err1(jrhs,mp)=0
        err2(jrhs,mp)=0
        do 26 ip=1,np
          diff=wk(ip) - fit(ip)
          if(dabs(werr(ip)) < 1d-10) diff=0d0
          err1(jrhs,mp)=dmax1(err1(jrhs,mp),dabs(diff)*fac)
   26     err2(jrhs,mp)=err2(jrhs,mp) + werr(ip)*diff**2*fac**2
      endif
C --- Fit and error for |m1-m2| ---
      if(mm <= mmax.and.mm /= mp) then
        irhs(mm)=irhs(mm)+1
        jrhs=irhs(mm)
        lm=((2*lmax-mm+1)*mm)/2 + 1
        call dpzero(fit,np)
        i=0
        do 27 ie=1,nx
        do 27 l=mm,lx(ie)
        i=i+1
        do 27 ip=1,np
   27   fit(ip)=fit(ip) + xie(ip,lm+l,ie)*b(i,jrhs,mm)
        err1(jrhs,mm)=0
        err2(jrhs,mm)=0
        do 28 ip=1,np
          diff=wk(ip) - fit(ip)
          if(dabs(werr(ip)) < 1d-10) diff=0d0
          err1(jrhs,mm)=dmax1(err1(jrhs,mm),dabs(diff)*fac)
   28     err2(jrhs,mm)=err2(jrhs,mm) + werr(ip)*diff**2*fac**2
      endif

   21 continue
   20 continue

C --- Collect errors by l1,l2 and output big table ---
      yy=100d0
      do  6  m=0, mmax
    6 irhs(m)=0
CL      if(ipr > 1) write(71,345) ep1,ep2,d,np
  345 format(/' ftcerr:  e1',f9.3,'  e2',f9.3,'  d',f9.4,
     .   '   np',i5)
      if (ipr >= 50) write(6,346)
  346 format(
     .  '  Prod   l1 m1   l2 m2    m    rms-err  max-err (%)')
      if (ipr >= 30.and.ipr < 50) write(6,347)
CL      if(ipr > 1) write(71,347)
  347 format(9x,' Prod     rms-err  max-err')
      do  31  l1=0, lp1
      do  31  l2=0, lp2
      arr1(l1,l2)=0d0
   31 arr2(l1,l2)=0d0
      do  35  l1=0, lp1
      do  35  m1=0, l1
      ltop=lp2
      if (lsym == 1) ltop=l1
      do  36  l2=0, ltop
      mtop=l2
      if (lsym == 1.and.l2 == l1) mtop=m1
      do  36  m2=0, mtop
      mp=m1+m2
      mm=iabs(m1-m2)
C --- Part for m3=m1+m2 ---
      if (mp <= mmax) then
        irhs(mp)=irhs(mp)+1
        jrhs=irhs(mp)
        err1(jrhs,mp)=err1(jrhs,mp)
        err2(jrhs,mp)=dsqrt(err2(jrhs,mp)/errvol)
        arr1(l1,l2)=dmax1(arr1(l1,l2),err1(jrhs,mp))
        arr2(l1,l2)=dmax1(arr2(l1,l2),err2(jrhs,mp))
        if(ipr >= 50) print 991, cl(l1),cl(l2),l1,m1,l2,m2,
     .    mp,err2(jrhs,mp)*yy,err1(jrhs,mp)*yy
  991   format(1x,a1,' * ',a1,2x,2i3,2x,2i3,i5,1x,2f9.2)
      endif
C --- Part for m3=iabs(m1-m2) ---
      if (mm <= mmax .and. mm /= mp) then
        irhs(mm)=irhs(mm)+1
        jrhs=irhs(mm)
        err1(jrhs,mm)=err1(jrhs,mm)
        err2(jrhs,mm)=dsqrt(err2(jrhs,mm)/errvol)
        arr1(l1,l2)=dmax1(arr1(l1,l2),err1(jrhs,mm))
        arr2(l1,l2)=dmax1(arr2(l1,l2),err2(jrhs,mm))
        if(ipr >= 50) print 993, mm,err2(jrhs,mm)*yy,err1(jrhs,mm)*yy
  993   format(22x,i5,1x,2f9.2)
      endif
   36 continue
   35 continue
C --- Accumulate final errs, output small table ---
      if(ipr >= 50) write(6,'('' Total:'')')
      do 38 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 38 l2=0,ltop
        xx2=yy*arr2(l1,l2)
        xx1=yy*arr1(l1,l2)
      if (ipr >= 50) write(6,670) cl(l1),cl(l2),xx2,xx1
      if (ipr >= 30.and.ipr < 50) write(6,671) cl(l1),cl(l2),xx2,xx1
CL      if(ipr > 1) write(71,671) cl(l1),cl(l2),xx2,xx1
      errs(l1,l2)=arr2(l1,l2)*100d0
      if(lsym == 1) errs(l2,l1)=errs(l1,l2)
   38 continue
  670 format(1x,a1,' * ',a1,22x,2f9.2)
  671 format(10x,a1,' * ',a1,1x,2f9.2)

      end
