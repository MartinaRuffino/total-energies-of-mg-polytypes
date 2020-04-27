      subroutine on2gen(n0,lx,ex,nx,lp,ep,np,llink,el,tl,
     .  rmt,rsm,rsmp,ri1,ri2,nph,ifi)
C- Expansion of products of on-site linked, unsmoothed hankels.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n0    :dimensioning parameter
Ci   lx    :product basis fitting l-cutoffs
Ci   ex    :product basis fitting energies
Ci   nx    :number of product basis fitting energies
Ci   lp    :lmto basis l-cutoff
Ci   ep    :lmto basis energy
Ci   np    :number of lmto basis energies
Ci   llink :1 for linked basis, otherwise 0
Ci   el    :lmto energies
Ci   tl    :(linked basis) linking energy
Ci   rmt   :augmentation radius, in a.u.
Ci   rsm   :smoothing radius
Ci   rsmp  :smoothing radius for basis function
Ci   ri1   :not used in this routine
Ci   ri2   :not used in this routine
Ci   ifi   :file logical unit (for saving fit)
Co Outputs
Co   nph   :number of product functions fit by this routine
Cl Local variables
Cl         :
Cr Remarks
Cr
Cb Bugs
Cb   eliminate ri1,ri2 from call
Cu Updates
Cu   24 Jun 04
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n0,lx(1),ifi,lp(1),nx,np,llink,nph
      double precision ex(1),rsmp,tl(n0,np)
      double precision el,ep(1),ri1,ri2,rmt,rsm
C ... Local parameters
      double precision epj1(20),epj2(20),chadd(25,5)
      integer lpj1(20),lpj2(20),lsymj(20),ip1,ip2,iph,
     .  ipj1(20),ipj2(20),ipr,nelx,nhl

      call getpr(ipr)
      do  5  ip1 = 1, np*np
    5 lsymj(ip1) = 0
      if (ipr >= 20) print 991, rmt
  991 format(/' on2gen:  rmt=',f7.4)
      if (llink == 1) then
C        if (ipr >= 20) print 333, el, tl
C  333   format(' linked basis, elink=',f8.4,'  tlink=',4f8.4)
C        if (ipr >= 20) print 333, el
C  333   format(' linked basis, elink=',f8.4)
C       call msadd(n0,np,rsmp,rmt,llink,el,tl,lp,ep,chadd)
        call msetpl(n0,np,'    ',1,1,rmt,np,1,llink,
     .    el,tl,lp,ep,nelx,nhl,chadd)
      endif
      call xxhypr('xi',rsm,rsm,lx,ex,nx,lx,ex,0,lsymj)

C --- Generate inequivalent fits ---
      nph = 0
      do  10  ip1 = 1, np
      do  10  ip2 = ip1, np
        nph = nph+1
        lpj1(nph) = lp(ip1)
        lpj2(nph) = lp(ip2)
        epj1(nph) = ep(ip1)
        epj2(nph) = ep(ip2)
        ipj1(nph) = ip1
        ipj2(nph) = ip2
        lsymj(nph) = 0
        if (ip1 == ip2) lsymj(nph) = 1
   10 continue
      if (ipr >= 20)
     .  call xxhypr('ph',0d0,0d0,lpj1,epj1,nph,lpj2,epj2,nph,lsymj)

      do  20  iph = 1, nph
        call on2mak(llink,lx,ex,nx,
     .    lpj1(iph),epj1(iph),chadd(1,ipj1(iph)),el,
     .    lpj2(iph),epj2(iph),chadd(1,ipj2(iph)),el,
     .    rmt,rsm,rsmp,ri1,ri2,ifi)
   20 continue

      end
