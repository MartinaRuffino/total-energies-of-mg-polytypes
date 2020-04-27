      subroutine tcferr(d,mmax,wp,np,wk,fit,
     .  nx,lx,xie,nlm,lmax,ph1,ph2,nlmp1,nlmp2,
     .  lp1,ep1,lp2,ep2,idxj,nph,lsymj,b,
     .  ndimx,nrhsx,err1,err2,prd,top,errs)
C- Generate rms and max errors of the TCF.
C  RMS errors printed as sqrt(intgrl(diff**2)/intgrl(ph1*ph2)**2)
C     but returned as intgrl(diff**2) for HYF.
C  Absolute errors relative to largest tabulated product of phi's
      implicit none
      integer lp1(1),lp2(1),lsymj(1),nlm,n0,ndimx,np,nrhsx,nx,
     .  mmax,irhs(0:20),lx(1),idxj(2,1),nlmp1,nlmp2,nph
      parameter (n0=10)
      double precision wp(1),wk(1),fit(1),ep1(1),ep2(1),
     .  xie(np,nlm,1),ph1(np,nlmp1,1),ph2(np,nlmp2,1),d,
     .  phi1(0:n0,0:n0),phi2(0:n0,0:n0),top(nrhsx,0:mmax),cx(200),
     .  b(ndimx,nrhsx,0:mmax),errs(1),
     .  err1(nrhsx,0:mmax),err2(nrhsx,0:mmax),prd(nrhsx,0:mmax)
      double precision diff,yy,xx1,xx2,err2x
      integer i,ie,ip,ipr,jrhs,l,l1,l2,lm,lmax,lmp1,lmp2,ltop,
     .  m,m1,m2,mm,mp,mtop,lp1j,lp2j,iph,offerr
      character*1 cl(0:8)
      data cl /'s','p','d','f','g','5','6','7','8'/

      call tcn('tcferr')
      call getpr(ipr)
      call sxlmnc(cx,10)

C --- For each (l1,m1,l2,m2,iph) accumulate fit errors ---
      do  5  m = 0, mmax
    5 irhs(m) = 0
      do  20  iph = 1, nph
      lp1j = lp1(idxj(1,iph))
      lp2j = lp2(idxj(2,iph))
      do  20  l1 = 0, lp1j
      do  20  m1 = 0, l1
      lmp1 = (l1*(l1+1))/2+m1+1
      ltop = lp2j
      if (lsymj(iph) == 1) ltop = l1
      do  21  l2 = 0, ltop
      mtop = l2
      if (lsymj(iph) == 1 .and. l2 == l1) mtop = m1
      do  21  m2 = 0, mtop
      lmp2 = (l2*(l2+1))/2+m2+1

C ... Hold wave function product in work array
      do  22  ip = 1, np
   22 wk(ip) = ph1(ip,lmp1,idxj(1,iph))*ph2(ip,lmp2,idxj(2,iph))

      mp = m1+m2
      mm = iabs(m1-m2)
C --- Fit and error for m1+m2 ---
      if (mp <= mmax) then
        irhs(mp) = irhs(mp)+1
        jrhs = irhs(mp)
        lm = ((2*lmax-mp+1)*mp)/2 + 1
        call dpzero(fit,np)
        i = 0
        do  25  ie = 1, nx
        do  25  l = mp, lx(ie)
        i = i+1
        do  25  ip = 1, np
   25   fit(ip) = fit(ip) + xie(ip,lm+l,ie)*b(i,jrhs,mp)
        err1(jrhs,mp) = 0
        err2(jrhs,mp) = 0
        prd(jrhs,mp)  = 0
        top(jrhs,mp) = -1000
        do  26  ip = 1, np
          diff = wk(ip) - fit(ip)
          err1(jrhs,mp) = dmax1(err1(jrhs,mp),dabs(diff))
          err2(jrhs,mp) = err2(jrhs,mp) + wp(ip)*diff**2
          top(jrhs,mp) = dmax1(top(jrhs,mp),dabs(wk(ip)))
          prd(jrhs,mp)  = prd(jrhs,mp)  + wp(ip)*wk(ip)**2
C          if (jrhs == 1 .and. mp == 0)
C     .      print 333, wk(ip), fit(ip), diff, wp(ip)
C  333     format(1p,4e17.8)
   26   continue
      endif
C --- Fit and error for |m1-m2| ---
      if (mm <= mmax .and. mm /= mp) then
        irhs(mm) = irhs(mm)+1
        jrhs = irhs(mm)
        lm = ((2*lmax-mm+1)*mm)/2 + 1
        call dpzero(fit,np)
        i = 0
        do  27  ie = 1, nx
        do  27  l = mm, lx(ie)
        i = i+1
        do  27  ip = 1, np
   27   fit(ip) = fit(ip) + xie(ip,lm+l,ie)*b(i,jrhs,mm)
        err1(jrhs,mm) = 0
        err2(jrhs,mm) = 0
        prd(jrhs,mm)  = 0
        top(jrhs,mm) = -1000
        do  28  ip = 1, np
          diff = wk(ip) - fit(ip)
          err1(jrhs,mm) = dmax1(err1(jrhs,mm),dabs(diff))
          err2(jrhs,mm) = err2(jrhs,mm) + wp(ip)*diff**2
          prd(jrhs,mm)  = prd(jrhs,mm)  + wp(ip)*wk(ip)**2
          top(jrhs,mm)  = dmax1(top(jrhs,mm),dabs(wk(ip)))
   28   continue
      endif

   21 continue
   20 continue

C --- Collect errors by l1,l2,iph and output big table ---
      yy = 100d0
      do  6  m = 0, mmax
    6 irhs(m) = 0
      offerr = 1
      do  30  iph = 1, nph
      lp1j = lp1(idxj(1,iph))
      lp2j = lp2(idxj(2,iph))
      if (ipr >= 40) write(6,*) ' '
      if (ipr >= 30) write(6,345) ep1(idxj(1,iph)),ep2(idxj(2,iph)),d
CL      write(71,345) ep1(idxj(1,iph)),ep1(idxj(2,iph)),d
  345 format(' tcferr: e1=',f7.4,'  e2=',f7.4,'  d=',f8.5,'  ----')
      if (ipr >= 40) write(6,346)
  346 format(
     .  '  Prod   l1 m1   l2 m2    m    rms-err  max-err (% of max)')
      if (ipr >= 30 .and. ipr < 40) write(6,347)
CL      write(71,347)
  347 format(' Prod     rms-err  max-err')
      do  31  l1 = 0, lp1j
      do  31  l2 = 0, lp2j
      phi1(l1,l2) = 0d0
   31 phi2(l1,l2) = 0d0
      do  35  l1 = 0, lp1j
      do  35  m1 = 0, l1
      ltop = lp2j
      if (lsymj(iph) == 1) ltop=l1
      do  36  l2 = 0, ltop
      mtop = l2
      if (lsymj(iph) == 1 .and. l2 == l1) mtop = m1
      do  36  m2 = 0, mtop
      mp = m1+m2
      mm = iabs(m1-m2)
C --- Part for m3 = m1+m2 ---
      if (mp <= mmax) then
        irhs(mp) = irhs(mp)+1
        jrhs = irhs(mp)
        err1(jrhs,mp) = err1(jrhs,mp)/top(jrhs,mp)
        err2x         = dsqrt(dmax1(err2(jrhs,mp)/prd(jrhs,mp),0d0))
        phi1(l1,l2) = dmax1(phi1(l1,l2),err1(jrhs,mp))
        phi2(l1,l2) = dmax1(phi2(l1,l2),err2x)
        if (ipr >= 40) print 991, cl(l1),cl(l2),l1,m1,l2,m2,
     .    mp,err2x*yy,err1(jrhs,mp)*yy,top(jrhs,mp)
  991   format(1x,a1,' * ',a1,2x,2i3,2x,2i3,i5,1x,2f9.2,'  (',f7.5,')')
      endif
C --- Part for m3 = iabs(m1-m2) ---
      if (mm <= mmax .and. mm /= mp) then
        irhs(mm) = irhs(mm)+1
        jrhs = irhs(mm)
        err1(jrhs,mm) = err1(jrhs,mm)/top(jrhs,mm)
        err2x         = dsqrt(dmax1(err2(jrhs,mm)/prd(jrhs,mm),0d0))
        phi1(l1,l2) = dmax1(phi1(l1,l2),err1(jrhs,mm))
        phi2(l1,l2) = dmax1(phi2(l1,l2),err2x)
        if (ipr >= 40) print 993, mm,err2x*yy,err1(jrhs,mm)*yy
  993   format(22x,i5,1x,2f9.2)
      endif
   36 continue
   35 continue
C --- Accumulate final errs, output small table ---
      if (ipr >= 40) write(6,'('' Total:'')')
      do  38  l1 = 0, lp1j
      ltop = lp2j
      if (lsymj(iph) == 1) ltop = l1
      do  38  l2 = 0, ltop
        xx2 = yy*phi2(l1,l2)
        xx1 = yy*phi1(l1,l2)
      if (ipr >= 40) write(6,670) cl(l1),cl(l2),xx2,xx1
      if (ipr >= 30 .and. ipr < 40) write(6,671) cl(l1),cl(l2),xx2,xx1
CL      write(71,671) cl(l1),cl(l2),xx2,xx1
      errs(l1+l2*(lp1j+1)+offerr) = xx2
      if (lsymj(iph) == 1) errs(l2+l1*(lp1j+1)+offerr) = xx2
   38 continue
  670 format(1x,a1,' * ',a1,22x,2f9.2)
  671 format(1x,a1,' * ',a1,1x,2f9.2)
      offerr = offerr + (lp1j+1)**2
   30 continue

      call tcx('tcferr')
      end
