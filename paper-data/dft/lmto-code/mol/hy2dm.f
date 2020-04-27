      subroutine hy2dm(mmax,lp1,ep1,np1,lx1,ex1,nx1,lp2,ep2,np2,
     .  lx2,ex2,nx2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,el1,tl1,el2,tl2,
     .  ndim,nrhs,nrhsj,ndimx,nsdmx,nrhsx,ncof,nxi,idxj,lsymj,nph,
     .  nalf,ndimb,ndimbx,nsdmbx)
C- Dimensions for hyfmak or hy2mak.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mmax   :
Ci   lp1    :lmto basis l-cutoff, function at first center
Ci   ep1    :lmto basis energy, function at first center
Ci   np1    :number of radial wave functions at first center
Ci   lx1    :product basis fitting l-cutoffs, first center
Ci   ex1    :product basis fitting energies, first center
Ci   nx1    :number of product basis energies, first center
Ci   lp2    :lmto basis l-cutoff, function at second center
Ci   ep2    :lmto basis energy, function at second center
Ci   np2    :number of radial wave functions at second center
Ci   lx2    :product basis fitting l-cutoffs, second center
Ci   ex2    :product basis fitting energies, second center
Ci   nx2    :number of product basis energies, second center
Ci   rmt1   :augmentation radius, in a.u., first center
Ci   rmt2   :augmentation radius, in a.u., second center
Ci   rsm1   :smoothing radius for p.b., first center
Ci   rsm2   :smoothing radius for p.b., second center
Ci   rsmp1  :smoothing radius basis fcn, first center
Ci   rsmp2  :smoothing radius basis fcn, second center
Ci   el1    :linking energy, linked basis, first center
Ci   tl1    :linking k.e.,   linked basis, first center
Ci   el2    :linking energy, linked basis, second center
Ci   tl2    :linking k.e.,   linked basis, second center
Ci   nrhsj  :
Ci   idxj   :
Ci   nph    :
Ci   nalf   :
Co Outputs
Co  lsymj   :1 if symmetric case, 0 otherwise
Co  ... the following are made in tcfdm
Co  ndim    :
Co  nrhs    :
Co  ndimx   :
Co  nsdmx   :
Co  nrhsx   :
Co  ncof    :
Co  nxi     :
Co  The following are analogs of ndim,ndimx,nsdmx for hyf
Co   ndimb  :
Co   ndimbx :
Co   nsdmbx :
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Jun 04
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mmax,ndim(0:mmax,0:1),nrhs(0:20),nrhsj(0:mmax,1),nph,
     .  ndimb(0:20),lx1(1),lx2(1),lp1(1),lp2(1),nalf,ncof,ndimbx,
     .  ndimx,nrhsx,nsdmbx,nsdmx,nx1,nx2,nxi,np1,np2,idxj(2,20),
     .  lsymj(20)
      double precision ex1(1),ex2(1),rsm1,rsm2,ep1(1),ep2(1),
     .  el1,tl1(1),el2,tl2(1),rmt1,rmt2,rsmp1,rsmp2
C ... Local parameters
      character*60 s1l,s2l, se(2)
      integer lpj1(20),lpj2(20),lsym,m
      integer iprint,ip1,ip2,jp,i,nsdb
      double precision epj1(20),epj2(20),fuzz,tlj1(20),tlj2(20)

C --- Print out specifications for xi's ---
      if (iprint() >= 20) then
        print 991, mmax,rmt1,rmt2
  991   format(/' hyfdm:  mmax=',i1,'  rmt=',2f8.5)
        call xxhypr('xi',rsm1,rsm2,lx1,ex1,nx1,lx2,ex2,nx2,lsymj)
      endif

C --- Find all inequivalent cases ---
      s1l = ' ph1:   l='
      s2l = ' ph2:   l='
      se(1) = '   el='
      se(2) = '   el='
      fuzz = 1d-10
      lsym = 1
      if (dabs(rsm1-rsm2)+dabs(rsmp1-rsmp2)+dabs(rmt1-rmt2)+
     .    dabs(el1-el2) > fuzz .or. nx1 /= nx2) lsym = 0
      do  4  i = 1, min0(nx1,nx2)
    4 if (lx1(i) /= lx2(i) .or. dabs(ex1(i)-ex2(i)) > fuzz) lsym = 0
C ... If now lsym=1, symmetry ok except for lp,ep mismatch
      nph = 0
      do  10  ip1 = 1, np1
      do  10  ip2 = 1, np2
        do  12  jp = 1, nph
C ...   Use earlier table if one is equivalent
        if (lp1(ip1) == lpj1(jp) .and. lp2(ip2) == lpj2(jp) .and.
     .    (dabs(tl1(ip1)-tlj1(jp)) < fuzz) .and.
     .    (dabs(tl2(ip2)-tlj2(jp)) < fuzz) .and.
     .    (dabs(ep1(ip1)-epj1(jp)) < fuzz) .and.
     .    (dabs(ep2(ip2)-epj2(jp)) < fuzz) .and. lsym == 1) goto 14
        if (lp2(ip2) == lpj1(jp) .and. lp1(ip1) == lpj2(jp) .and.
     .    (dabs(tl1(ip1)-tlj2(jp)) < fuzz) .and.
     .    (dabs(tl2(ip2)-tlj1(jp)) < fuzz) .and.
     .    (dabs(ep2(ip2)-epj1(jp)) < fuzz) .and.
     .    (dabs(ep1(ip1)-epj2(jp)) < fuzz) .and. lsym == 1) goto 14
   12   continue
        nph = nph+1
        lpj1(nph) = lp1(ip1)
        lpj2(nph) = lp2(ip2)
        epj1(nph) = ep1(ip1)
        epj2(nph) = ep2(ip2)
        tlj1(nph) = tl1(ip1)
        tlj2(nph) = tl2(ip2)
        lsymj(nph) = lsym
        idxj(1,nph) = ip1
        idxj(2,nph) = ip2
        if (lpj1(nph) /= lpj2(nph) .or.
     .      dabs(epj1(nph)-epj2(nph)) > fuzz) lsymj(nph) = 0
   14   continue
   10 continue
      if (iprint() >= 20)
     .  call xxhypr('ph',rsmp1,rsmp2,lpj1,epj1,nph,lpj2,epj2,nph,lsymj)

C --- Get the dimension ---
      call tcfdm(mmax,lx1,nx1,lx2,nx2,lpj1,lpj2,nph,lsymj,
     .   ndim,nrhs,nrhsj,ndimx,nsdmx,nrhsx,ncof,nxi)
      nsdmbx = 0
      ndimbx = 0
      do  3  m = 0, mmax
      ndimb(m) = ndim(m,0)*nalf
      nsdb = (ndimb(m)*(ndimb(m)+1))/2
      ndimbx = max0(ndimbx,ndimb(m))
    3 nsdmbx = max0(nsdmbx,nsdb)
      if (iprint() >= 30) write(6,300) nsdmbx,(ndimb(m),m=0,mmax)
  300 format(' nsdmbx=',i6,'  ndimb=',12i5)
      end

      subroutine xxhypr(s0,rsm1,rsm2,l1,e1,n1,l2,e2,n2,lsymj)
C- Printout for l's and e's for a pair of Hankels
      implicit none

      character sl(2)*60, se(2)*60, s0*2
      integer l1,n1,l2,n2,ise(2),isl(2),i,lsymj(1)
      double precision rsm1,rsm2,e1,e2


      do  10  i = 1, 2
        sl(i) = s0 // char(ichar('1')+i-1) // ':   rsm='
        se(i) = '   e='
        isl(i) = 11
   10 continue
      if (n2 == 0) sl(1) = s0 // ':    rsm='
      call bin2a('(f10.5)',0,0,rsm1,4,0,60,sl(1),isl(1))
      call bin2a('(f10.5)',0,0,rsm2,4,0,60,sl(2),isl(2))
      isl(1) = max(17, 2+max(isl(1),isl(2)))
      isl(2) = isl(1)
      do  12  i = 1, 2
   12 call bin2a('l=',0,0,0,1,0,60,sl(i),isl(i))
      ise(1) = 6


      do  20  i = 1, max(n1,n2)
        isl(2) = isl(1)
        ise(2) = ise(1)
        if (i <= n1) then
          call bin2a('',1,0,l1,2,i-1,60,sl(1),isl(1))
          call bin2a('(f10.5)',1,0,e1,4,i-1,60,se(1),ise(1))
        endif
        if (i <= n2) then
          call bin2a('',1,0,l2,2,i-1,60,sl(2),isl(2))
          call bin2a('(f10.5)',1,0,e2,4,i-1,60,se(2),ise(2))
        endif
        if (lsymj(i) == 1) then
          ise(1) = ise(1)+1
          se(1)(ise(1):ise(1)) = '*'
        endif
        ise(1) = max(ise(1),ise(2))
        isl(1) = max(isl(1),isl(2))
   20 continue
      ise(1) = min(ise(1),60)

      print 333, sl(1)(1:isl(1)), se(1)(1:ise(1))
      if (n2 > 0) print 333, sl(2)(1:isl(1)), se(2)(1:ise(1))
  333 format(1x,a,a)

      end
