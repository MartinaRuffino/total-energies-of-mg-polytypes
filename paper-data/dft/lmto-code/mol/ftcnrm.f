      subroutine ftcnrm(mmax,lmax,lx,nx,wk,wp,np,xie,nlm,cx,
     .  ph1,ph2,lp1,lp2,lsym,s,b,ndimx,nsdmx,nrhsx)
C- Adds into integrals for tcf (s=overlap and b=rhs)
      implicit real*8 (a-h,p-z), integer (o)
      dimension irhs(0:20),lx(1),wk(1),cx(200),s(nsdmx,0:mmax),
     .   ph1(np,1),ph2(np,1),b(ndimx,nrhsx,0:mmax),xie(np,nlm,1),
     .   wp(1)

C ------------------ Normal matrix ----------------------------
      if (iprint() > 60) call tm('start lhs')
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
        call wpdot(xie(1,lm+li,ie),xie(1,lm+lj,je),wp,np,sum)
        s(is,m)=s(is,m)+sum
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
      if (lsym == 1.and.l2 == l1) mtop = m1
      do  21  m2 = 0, mtop
      lmp2 = (l2*(l2+1))/2+m2+1

C ... Hold wave function product in work array
      do  22  ip = 1, np
   22 wk(ip) = ph1(ip,lmp1)*ph2(ip,lmp2)

      mp = m1+m2
      mm = iabs(m1-m2)
C --- rhs for m1+m2 ---
      if (mp <= mmax) then
        irhs(mp) = irhs(mp)+1
        jrhs = irhs(mp)
        lm = ((2*lmax-mp+1)*mp)/2 + 1
        i = 0
C ...   do  25  i = 1, ndim(mp)
        do  25  ie = 1, nx
        do  25  l = mp, lx(ie)
          i = i+1
          call wpdot(wk,xie(1,lm+l,ie),wp,np,sum)
  25      b(i,jrhs,mp)=b(i,jrhs,mp)+sum
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
          call wpdot(wk,xie(1,lm+l,ie),wp,np,sum)
  27      b(i,jrhs,mm)=b(i,jrhs,mm)+sum
      endif

   21 continue
   20 continue

      end

      subroutine wpdot(x,y,w,n,sum)
      implicit real*8 (a-h,p-z), integer (o)
      dimension x(1),y(1),w(1)
      sum=0d0
      do 1 i=1,n
  1   sum=sum+x(i)*y(i)*w(i)
      end
