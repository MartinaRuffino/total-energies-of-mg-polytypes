      subroutine tcfnrm(mmax,lmax,lx,nx,wk,np,nbisi,xie,nlm,ph1,ph2,
     .  lp1,lp2,nlmp1,nlmp2,idxj,nph,lsymj,s,b,ndimx,nsdmx,nrhsx)
C- Makes integrals for tcf (s=overlap and b=rhs)
      implicit none
      integer irhs(0:20),lx(1),lp1(1),lp2(1),lsymj(1),ndimx,np,nph,
     .  nrhsx,nsdmx,nx,nlm,mmax,nbisi,idxj(2,1),nlmp1,nlmp2
      double precision wk(1),s(nsdmx,0:mmax),ph1(np,nlmp1,1),sum1,
     .   xie(np,nlm,1),b(ndimx,nrhsx,0:mmax),ph2(np,nlmp2,1),sum2
      integer i,ie,ip,is,j,je,jrhs,l,l1,l2,lmp1,lmp2,lp1j,lp2j,
     .  li,lj,lm,lmax,ltop,m,m1,m2,mm,mp,mtop,iprint,iph

      call tcn('tcfnrm')

C ------------------ Normal matrix ----------------------------
      do  40  m = 0, mmax
        lm = ((2*lmax-m+1)*m)/2 + 1
        is = 0
C ...   do  41  j = 1, ndim(m)
        j = 0
        do  41  je = 1, nx
        do  41  lj = m, lx(je)
        j = j+1
C ...   do  42  i = 1, j
        i = 0
        do  42  ie = 1, nx
        do  42  li = m, lx(ie)
        i = i+1
        if (i > j) goto 42
        is = is+1

        sum1 = 0
        sum2 = 0
        do  10  ip = 1, nbisi
   10   sum1 = sum1 + xie(ip,lm+li,ie)*xie(ip,lm+lj,je)
        do  11  ip = nbisi+1, np
   11   sum2 = sum2 + xie(ip,lm+li,ie)*xie(ip,lm+lj,je)
        s(is,m) = sum1 - sum2
C       if (m == 0) print *, 'tcfnrm:', i,j,s(is,m)
   42   continue
   41   continue
   40 continue

C ------------------ Make rhs ----------------------------
      do  66  m = 0, mmax
   66 irhs(m) = 0

C ... For each (l1,m1,l2,m2) do ...
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
C --- rhs for m1+m2 ---
      if (mp <= mmax) then
        irhs(mp) = irhs(mp)+1
        jrhs = irhs(mp)
        lm = ((2*lmax-mp+1)*mp)/2 + 1
        i = 0
C ...   do  24  i = 1, ndim(mp)
        do  24  ie = 1, nx
        do  24  l = mp, lx(ie)
          i = i+1
          sum1 = 0
          sum2 = 0
          do  25  ip = 1, nbisi
   25     sum1 = sum1 + wk(ip)*xie(ip,lm+l,ie)
          do  26  ip = nbisi+1, np
   26     sum2 = sum2 + wk(ip)*xie(ip,lm+l,ie)
          b(i,jrhs,mp) = sum1 - sum2
C          if (mp == 0 .and. jrhs == 1)
C     .      print *, 'tcf b:', i,b(i,jrhs,mp)
   24   continue
      endif
C --- rhs for |m1-m2| ---
      if (mm <= mmax .and. mm /= mp) then
        irhs(mm) = irhs(mm)+1
        jrhs = irhs(mm)
        lm = ((2*lmax-mm+1)*mm)/2 + 1
        i = 0
C ...   do  27  i = 1, ndim(mm)
        do  27  ie = 1, nx
        do  27  l = mm, lx(ie)
          i = i+1
          sum1 = 0
          sum2 = 0
          do  28  ip = 1, nbisi
   28     sum1 = sum1 + wk(ip)*xie(ip,lm+l,ie)
          do  29  ip = nbisi+1, np
   29     sum2 = sum2 + wk(ip)*xie(ip,lm+l,ie)
          b(i,jrhs,mm) = sum1 - sum2
   27   continue
      endif

   21 continue
   20 continue

      call tcx('tcfnrm')
      end
