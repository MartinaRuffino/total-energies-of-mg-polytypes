      subroutine tcfdm(mmax,lx1,nx1,lx2,nx2,lpj1,lpj2,nph,lsymj,
     .  ndim,nrhs,nrhsj,ndimx,nsdmx,nrhsx,ncof,nxi)
C- Get dimensions of the m-spaces, as well as number of rhs
C  Outputs:
C    ndim, nrhs:  dimension of normal matrix, and No of RHS for each m
C          nrhsj: nrhs for each energy pair
C    ndimx,nrhsx: max ndim and nrhs
C    nsdmx:       max dimensioning needed for normal matrix
C    ncof:        total number of coefficients.
C      implicit none
      integer ndim(0:mmax,0:1),nrhs(0:20),nrhsj(0:mmax,1),lx1(1),lx2(1),
     .  ipr,l1,l2,lpj1(1),lpj2(1),lsymj(1),ltop,m,m1,m2,mm,mmax,mp,mtop,
     .  ncof,nd,ndimx,nf1,nf2,nrhst,nrhsx,nsdm,nsdmx,nx1,nx2,nxi,iph,
     .  nph,itop,ip1,ip2,i

      call getpr(ipr)
C --- Count number of rhs for each m ---
      do  10  m = 0, mmax
      nrhsj(m,1) = 0
   10 nrhs(m) = 0
      do  25  iph = 1, nph
        nf1 = 0
        do  20  l1 = 0, lpj1(iph)
        do  20  m1 = 0, l1
          nf1 = nf1+1
          nf2 = 0
          ltop = lpj2(iph)
          if (lsymj(iph) == 1) ltop = l1
          do  21  l2 = 0, ltop
            mtop = l2
            if (lsymj(iph) == 1 .and. l2 == l1) mtop = m1
            do  21  m2 = 0, mtop
            nf2 = nf2+1
            mp = m1+m2
            mm = iabs(m1-m2)
            if (mp <= mmax) nrhs(mp) = nrhs(mp)+1
            if (mm /= mp .and. mm <= mmax) nrhs(mm) = nrhs(mm)+1
   21     continue
   20   continue
        do  24  m = 0, mmax
   24   nrhsj(m,iph+1) = nrhs(m)
   25 continue

C --- Get nxi ---
      nxi = 0
      do  28  i = 1, nx1
   28 nxi = nxi + (lx1(i)+1)**2
      do  29  i = 1, nx2
   29 nxi = nxi + (lx2(i)+1)**2

C --- Get ndim, nsdm and printout ---
      ndimx = 0
      nsdmx = 0
      nrhsx = 0
      nrhst = 0
      ncof = 0
      if (ipr >= 30) write(6,942)
      do  30  m = 0, mmax
      nd = 0
      do  32  i = 1, nx1
        if (m <= lx1(i)) nd = nd + (lx1(i)-m+1)
        ndim(m,i) = nd
   32 continue
      do  33  i = 1, nx2
        if (m <= lx2(i)) nd = nd + (lx2(i)-m+1)
        ndim(m,i+nx1) = nd
   33 continue
      ndim(m,0) = nd
      nsdm = (nd*(nd+1))/2
      ndimx = max0(ndimx,ndim(m,0))
      nsdmx = max0(nsdmx,nsdm)
      nrhsx = max0(nrhsx,nrhs(m))
      nrhst = nrhst + nrhs(m)
      if (ipr >= 30) write(6,941) m,ndim(m,0),nsdm,nrhs(m)
  941 format(i5,3i10)
  942 format(/'    m       ndim     nsdm       nrhs')
   30 ncof = ncof + ndim(m,0)*nrhs(m)
      if (ipr >= 30) print 950, ndimx,nrhsx,nxi,nsdmx,nrhst,ncof
  950 format( ' max ndim:',i5,'      max nrhs:',i5,
     .  i12,' fit functions'
     .       /' max nsdm:',i5,'      tot nrhs:',i5,
     .  i12,' nonzero coeffs')
      end



