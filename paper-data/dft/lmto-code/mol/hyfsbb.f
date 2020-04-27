      subroutine hyfsbb(ndimbx,nrhsx,mmax,nsdmbx,nsdmx,ndimx,
     .  ndim,nrhs,b,s,nz,nalf,f,wz,wk,kpvt,sm,sb,bb)
C- Generate normal matrix and RHS for hyf from same of TCF; solve
      implicit none
      integer ndimbx,nrhsx,mmax,nsdmbx,nsdmx,ndimx,nz,nalf,
     .  ndim(0:mmax),nrhs(0:mmax),ierr,ipr,kpvt(nalf*ndimx)
      double precision sb(nsdmbx),wk(ndimbx,1),fac,wgt
C#ifndef CRAY-DP
      double precision
C#elseC
C      real*8
C#endif
     . bb(ndimx,nrhsx,0:mmax,nalf),
     . s(nsdmx,0:mmax,nz),b(ndimx,nrhsx,0:mmax,nz),sm(ndimx,ndimx),
     . wz(1),f(nz,nalf)
      integer i,j,ialf,jalf,ia,ja,ir,ia0,ja0,id,nd,m,is,iz,isb,nda,i1


      call getpr(ipr)
      call tcn('hyfsbb')
C --- For each m set up and solve normal equations ---
      do  10  m = 0, mmax
      nd = ndim(m)
      nda= nd*nalf

      j = (nda*(nda+1))/2
      do  12  i = 1, j
   12 sb(i) = 0d0

      j = ndimbx*nrhs(m)
      do  14  i = 1, j
   14 wk(i,1) = 0d0

C --- Add to normal matrix and rhs for this point ---
      do  20  iz = nz, 1, -1
        wgt = wz(iz)

C   ... Make sb
C        do  30  jalf = 1, nalf
C          ja0 = (jalf-1)*nd
C          do  30  ialf = 1, jalf
C          ia0 = (ialf-1)*nd
C          fac = wgt*f(iz,ialf)*f(iz,jalf)
C          is = 0
C          do  31  j = 1, nd
C          ja = j+ja0
C          isb = (ja*(ja-1))/2+ia0
C          do  31  i = 1, j
C            is = is+1
C            sb(i+isb) = sb(i+isb) + fac*s(is,m,iz)
C   31     continue
CC ...     For off-diagonal blocks, copy s(j,i) into sb(ia,ja)
C          if (ialf /= jalf) then
C          is = 0
C          do  32  i = 1, nd
C          ia = i+ia0
C          do  33  j = 1, i-1
C            ja = j+ja0
C            isb = (ja*(ja-1))/2
C            is = is+1
C            sb(ia+isb) = sb(ia+isb) + fac*s(is,m,iz)
C   33     continue
C          is = is+1
C   32     continue
C          endif
C   30   continue

C   ... Make sb
        call dpmucp(s(1,m,iz),1,nd,1,nd,sm,ndimx)
        do  30  jalf = 1, nalf
          ja0 = (jalf-1)*nd
          do  30  ialf = 1, jalf
          ia0 = (ialf-1)*nd
          fac = (wgt*f(iz,ialf))*f(iz,jalf)
          do  31  j = 1, nd
          ja = j+ja0
          isb = (ja*(ja-1))/2 + ia0
          i1 = nd
          if (ialf == jalf) i1 = j
          do  31  i = 1, i1
   31     sb(i+isb) = sb(i+isb) + fac*sm(i,j)
   30   continue

C   ... Add to bb
        do  40  ialf = 1, nalf
        ia0 = (ialf-1)*nd
        fac = f(iz,ialf)*wgt
        do  40  ir = 1, nrhs(m)
        do  40  id = 1, nd
   40   wk(id+ia0,ir) = wk(id+ia0,ir) + fac*b(id,ir,m,iz)
   20 continue

C --- Solve the normal equations ---
C#ifndef CRAY-DP
      call dspfa(sb,nda,kpvt,ierr)
C#elseC
C      call xspfa(sb,nda,kpvt,ierr)
C#endif
      if (ierr /= 0)
     .  call rx('hyfsbb: matrix singular for m='//char(m+ichar('0')))
      do  50  ir = 1, nrhs(m)
C#ifndef CRAY-DP
        call dspsl(sb,nda,kpvt,wk(1,ir))
C#elseC
C      call xspsl(sb,nda,kpvt,wk(1,ir))
C#endif
        if(ipr >= 70) write(6,220) ir,(wk(i,ir),i=1,nda)
  220   format(' jrhs=',i5/(1x,8f10.5))

C ...   Copy into bb
        do  54  ialf = 1, nalf
        ia0 = (ialf-1)*nd
        do  54  id = 1, nd
        bb(id,ir,m,ialf) = wk(id+ia0,ir)
   54   continue
   50 continue
   10 continue

      call tcx('hyfsbb')
      end
