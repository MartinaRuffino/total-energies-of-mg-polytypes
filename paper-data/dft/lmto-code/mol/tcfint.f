      subroutine tcfint(mmax,lmax,lx,nx,wk,np,nbisi,lmxl,xie,nlm,cx,
     .  phi1,phi2,lp1,lp2,lsym,s,b,ndimx,nsdmx,nrhsx)
C- Makes integrals for tcf (s=overlap and b=rhs)
C  np,nbisi,lmxl: number of points, number of points in strict
C  interstitial, and lmxl.  Points between nbisi and np are included
C  in inner product for m>=lmxl.
      implicit none
      integer irhs(0:20),lx(1),lp1,lp2,lsym,ndimx,np,nrhsx,nsdmx,nx,
     .  nlm,mmax,nbisi,lmxl
      double precision wk(1),s(nsdmx,0:mmax),phi1(np,1),phi2(np,1),
     .  b(ndimx,nrhsx,0:mmax),xie(np,nlm,1),sum,cx(200)
      integer i,ie,ip,is,j,je,jrhs,l,l1,l2,lmp1,lmp2,
     .  li,lj,lm,lmax,ltop,m,m1,m2,mm,mp,mtop,iprint,np0

C ------------------ Normal matrix ----------------------------
      if (iprint() > 60) call tm('start lhs')
C*$* ASSERT CONCURRENT CALL
      do  40  m = 0, mmax
        np0 = np
        if (m <= lmxl) np0 = nbisi
        lm = ((2*lmax-m+1)*m)/2 + 1
        is = 0
C ...   do  41  i = 1, ndim(m)
        i = 0
        do  41  ie = 1, nx
        do  41  li = m, lx(ie)
        i = i+1

C ...   do  42  j = 1, i
        j = 0
        do  42  je = 1, nx
        do  42  lj = m, lx(je)
        j = j+1
        if (j > i) goto 42
        is = is+1

C        sum = 0d0
C        do  10  ip = 1, np0
C   10   sum = sum + xie(ip,lm+li,ie)*xie(ip,lm+lj,je)
C       s(is,m) = sum
        call dpdot(xie(1,lm+li,ie),xie(1,lm+lj,je),np0,s(is,m))
   42   continue
   41   continue
   40 continue


C ------------------ Make rhs ----------------------------
      if (iprint() > 60) call tm('start rhs')
      do  66  m = 0, mmax
   66 irhs(m) = 0

C ... For each (l1,m1,l2,m2) do ...
      do  20  l1 = 0, lp1
      do  20  m1 = 0, l1
      lmp1 = (l1*(l1+1))/2+m1+1
      ltop = lp2
      if (lsym == 1) ltop = l1
      do  21  l2 = 0, ltop
      mtop = l2
      if (lsym == 1 .and. l2 == l1) mtop = m1
      do  21  m2 = 0, mtop
      lmp2 = (l2*(l2+1))/2+m2+1

C ... Hold wave function product in work array
      do  22  ip = 1, np
   22 wk(ip) = phi1(ip,lmp1)*phi2(ip,lmp2)

      mp = m1+m2
      mm = iabs(m1-m2)
C --- rhs for m1+m2 ---
      if (mp <= mmax) then
        np0 = np
        if (mp <= lmxl) np0 = nbisi
        irhs(mp) = irhs(mp)+1
        jrhs = irhs(mp)
        lm = ((2*lmax-mp+1)*mp)/2 + 1
        i = 0
C*$* ASSERT CONCURRENT CALL
        do  25  ie = 1, nx
        do  25  l = mp, lx(ie)
          i = i+1
C          sum = 0d0
C          do  26  ip = 1, np0
C   26     sum = sum + wk(ip)*xie(ip,lm+l,ie)
C          b(i,jrhs,mp) = sum
          call dpdot(wk,xie(1,lm+l,ie),np0,b(i,jrhs,mp))
   25   continue
      endif
C --- rhs for |m1-m2| ---
      if (mm <= mmax .and. mm /= mp) then
        np0 = np
        if (mm <= lmxl) np0 = nbisi
        irhs(mm) = irhs(mm)+1
        jrhs = irhs(mm)
        lm = ((2*lmax-mm+1)*mm)/2 + 1
        i = 0
C*$* ASSERT CONCURRENT CALL
        do  27  ie = 1, nx
        do  27  l = mm, lx(ie)
          i = i+1
C          sum = 0
C          do  28  ip = 1, np0
C   28     sum = sum + wk(ip)*xie(ip,lm+l,ie)
C          b(i,jrhs,mm) = sum
          call dpdot(wk,xie(1,lm+l,ie),np0,b(i,jrhs,mm))
   27   continue
      endif

   21 continue
   20 continue

      if (iprint() > 60) call tm('done rhs')
      end
