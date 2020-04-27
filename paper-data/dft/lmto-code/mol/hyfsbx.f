      subroutine hyfsbx(ncupl,ndust,ndimbx,nrhsx,mmax,nsdmbx,nsdmx,
     .  ndimx,ndim,nrhs,b,s,nz,nalf,f,wz,wk,wk2,wk3,err,sm,sb,bb)
C- Generate normal matrix and RHS for hyf from same of TCF; solve
C  in successive subblocks
      implicit none
      integer ndimbx,nrhsx,mmax,nsdmbx,nsdmx,ndimx,nz,nalf,ncupl,
     .  ndim(0:mmax),nrhs(0:mmax),ierr,ipr
      double precision sb(nsdmbx),wk(ndimbx,1),wk2(ndimx,1),
     .  wk3(ndimx,1),fac,wgt,err(nrhsx,3)
C#ifndef CRAY-DP
      double precision
C#elseC
C      real*8
C#endif
     . bb(ndimx,nrhsx,0:mmax,nalf),
     . s(nsdmx,0:mmax,nz),b(ndimx,nrhsx,0:mmax,nz),sm(ndimx,ndimx),
     . wz(1),f(nz,nalf),v
      integer i,j,ialf,jalf,ja,ir,ia0,nd,m,iz,isb,nda,i1,
     .  iblk,nblk,ndust,n,idust

C     call pshpr(70)
C     ndust = 6

      call getpr(ipr)
      call tcn('hyfsbx')
C --- For each m and subblock set up and solve normal equations ---
      do  10  m = 0, mmax
      do  10  idust = 1, max(ndust,1)
      call dpzero(err,nrhs(m))
      call dpzero(err(1,2),nrhs(m))
      do  11  iblk = 1, nalf, ncupl
      nblk = min(ncupl+iblk-1,nalf)
      nd = ndim(m)
      nda = nd*(nblk-iblk+1)

C ... zero out initial hyf S and RHS for sum_iz
      j = (nda*(nda+1))/2
      do  12  i = 1, j
   12 sb(i) = 0d0
      j = ndimbx*nrhs(m)
      do  14  i = 1, j
   14 wk(i,1) = 0d0

C --- Add to normal matrix and rhs for this iz and subblock ---
      do  20  iz = nz, 1, -1
        wgt = wz(iz)

C   ... Add to HYF S
        call dpmucp(s(1,m,iz),1,nd,1,nd,sm,ndimx)
        do  30  jalf = iblk, nblk
          do  30  ialf = iblk, jalf
          fac = (wgt*f(iz,ialf))*f(iz,jalf)
          do  31  j = 1, nd
          ja = j + (jalf-iblk)*nd
          isb = (ja*(ja-1))/2 + (ialf-iblk)*nd
          i1 = nd
          if (ialf == jalf) i1 = j
          do  31  i = 1, i1
            if (i+isb == 99) call snitt
   31     sb(i+isb) = sb(i+isb) + fac*sm(i,j)
   30   continue

C   ... Add to HYF B
        do  40  ialf = iblk, nblk
        ia0 = (ialf-iblk)*nd
        fac = f(iz,ialf)*wgt
        do  40  ir = 1, nrhs(m)
        do  42  i  = 1, nd
        err(ir,2) = err(ir,2) + bb(i,ir,m,ialf)*fac*b(i,ir,m,iz)
   42   wk(i+ia0,ir) = wk(i+ia0,ir) + fac*b(i,ir,m,iz)
   40   continue

C ...   Subtract HYF S  * coffs for previously determined coffs
C       (iblk-1 blocks if first time through; otherwise all blocks)
        do  50  ialf = iblk, nblk
          ia0 = (ialf-iblk)*nd
          call dpzero(wk3,ndimx*nrhs(m))
          n = iblk-1
          if (idust > 1) n = nalf
          do  51  jalf = 1, n
          do  52  ir = 1, nrhs(m)
          fac = f(iz,jalf)
          do  52  i  = 1, nd
   52     wk2(i,ir) = fac*bb(i,ir,m,jalf)
          call dgemm('N','N',nd,nrhs(m),nd,1d0,sm,ndimx,wk2,ndimx,
     .      1d0,wk3,ndimx)
   51     continue
          do  54  ir = 1, nrhs(m)
            fac = f(iz,ialf)*wgt
            do  55  i  = 1, nd
   55       wk(i+ia0,ir) = wk(i+ia0,ir) - fac*wk3(i,ir)
C ...       For monitoring errors ...  accumulate coff * S * coff
            do  56  i  = 1, nd
   56       err(ir,1) = err(ir,1) + bb(i,ir,m,ialf)*fac*wk3(i,ir)
   54     continue
   50   continue
   20 continue

C     call wkchk(' after 20')

C --- Solve the normal equations ---
C#ifndef CRAY-DP
C      print *, 'iblk,idust,m=',iblk,idust,m
C      call dupack(sb,nda)
C      call prmx('sb',sb,nda,nda,nda)
C      call prmx('wk',wk,nda,nda,1)
C      call dpack(sb,nda)
      call dspfa(sb,nda,wk2,ierr)
C#elseC
C      call xspfa(sb,nda,wk2,ierr)
C#endif
      if (ierr /= 0)
     .  call rx('hyfsbb: matrix singular for m='//char(m+ichar('0')))
      do  60  ir = 1, nrhs(m)
C#ifndef CRAY-DP
        call dspsl(sb,nda,wk2,wk(1,ir))
C#elseC
C      call xspsl(sb,nda,wk2,wk(1,ir))
C#endif
        if (ipr >= 70) then
          if (ir == 1) then
            call dpdot(wk(1,ir),wk(1,ir),nda,v)
            v = dsqrt(v)/nda
            v = 10d0**nint(dlog10(v*10))
            v = max(v,1d-9)
          endif
          print 220, m,idust,ir,-nint(dlog10(v)), (wk(i,ir)/v, i=1,nda)
  220     format(' change in coffs for m=',i1,', idust=',i1,', col',i5,
     .      ': scaled by 10**',i2/(1x,8f10.5))
        endif

C ---   Sum this iteration into bb ---
        do  74  ialf = iblk, nblk
        ia0 = (ialf-iblk)*nd
        do  74  i = 1, nd
        bb(i,ir,m,ialf) = bb(i,ir,m,ialf) + wk(i+ia0,ir)
   74   continue

   60 continue
   11 continue

C ---   Printout ---
      if ((ipr > 40 .or. (idust >= ndust .and. ipr >= 30))
     .    .and. idust > 1) then
        print 344, m, idust
  344   format(' hyfsbx: convergence (*100) of each rhs for m =',i2,
     .    ', idust =',i2)
        print 345, (100*(err(ir,1)-err(ir,2))/
     .    max(err(ir,1),1d-16), ir = 1,nrhs(m))
  345   format(9f8.4:/(9f8.4))
      endif

C      do  84  ir = 1, nrhs(m)
C        do  82  ialf = 1, nalf
C        ia0 = (ialf-1)*nd
C        do  82  i = 1, nd
C   82   wk(i+ia0,ir) = bb(i,ir,m,ialf)
C   84 continue

C     call prmx('b this iter',wk,nd*nalf,nd*nalf,1)
   10 continue

      call tcx('hyfsbx')
      end
      subroutine snitt
C     print *, 'snitt'
      end
