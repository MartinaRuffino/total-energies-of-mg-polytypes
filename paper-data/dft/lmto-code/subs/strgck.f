      subroutine strgck(ivl,nds,nl,iat,ips,rsm,ehvl,alat,plat,pos,
     .  iax,ntab,ntabg,xmrp,sg,slj,slk)
C- Evaluate value-Laplacian functions at point xmrp without using
C  one-center expansion (for testing)
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies the functions used to built the value-Laplacian
Ci          set U_iRL (see Remarks)
Ci          ivl = 0: G0 & G1
Ci          ivl = 1: G0 & Hs
Ci          ivl = 2: G0 & Hs-dot
Ci   nds   :leading dimensions of sg
Ci   nl    :maximum l quantum number + 1, leading dimension of rsm,ehvl
Ci   iat   :center R of the cluster for which the value-Laplacian functions
Ci          are constructed
Ci   ips   :index to which species each site belongs
Ci   rsm   :smoothing radii for each species and angular momentum
Ci   ehvl  :Hankel energies for each species and angular momentum
Ci          not used if ivl = 0
Ci   plat,alat :primitive lattice translation vectors and scale
Ci   pos   :basis vectors
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Ci   ntabg :ntabg(ib) no. of neighbors in cluster centered at ib (pairg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   xmrp  :Cartesian coordinates of the point at which to evaluate the
Ci          value-Laplacian functions
Ci   sg    :screened strux in (nl**2,nl**2) blocks for all pairs
Ci          (1..ntabg(ib),ib=1..nbas)
Co Outputs
Co   slj(i,L) :approximate U_iRL, kappa- and L-resolved value-Laplacian
Co            functions at xmrp
Co   slk(i,L) :slk(i,L) = \lap slj(i,L)
Cl Local variables
Cl   nlmh  :Number of L channels at iat = R
Cl   nlmg  :Number of L channels at R" entering into screened function
Cl   gex   : F_0R"L" in expression (1), Remarks
Cl   hex   : F_1R"L" in expression (1), Remarks
Cr Remarks
Cr   slj are calculated by applying the value-Laplacian strux
Cr   to a set of actual functions evaluated at point r (= xmrp) rather
Cr   then to their expansion into polynomials
Cr           slj(i,L) = \sum_i"R"L" sg(L",L,i",i,R") * F_i"R"L"(r)    (1)
Cr   where the summation runs over all sites and L channels of a given
Cr   cluster, and F_jR'L'(r) = F_jL'(r-R') are Gaussians if j=0 and
Cr   G1/Hs/Hs-dot (depending on ivl) if j=1.
Cr
Cr   Since sg were built for functions expanded into polynomials,
Cr   slj and slk are not exactly U_iRL and \lap U_iRL, but tend
Cr   to them as the polynomial expansion becomes more and more accurate.
Cu Updates
Cu   15 May 08 Adapted from strck.f
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ivl,iat,nds,nl
      integer niax,nkap
      parameter (niax=10, nkap=2)
      integer iax(niax,*),ntab(*),ntabg(*),ips(*)
      double precision alat,plat(3,3),pos(3,*),rsm(nl,*),ehvl(nl,*)
      double precision xmrp(3),slj(nkap,nds),slk(nkap,nds),
     .  sg(nds,nds,nkap,nkap,*)
C Local parameters
      integer nclus,nlmh,nlmg,ll,lmaxg,ikap,pmax
      integer offR,isg,it,isj,ib,jb,ii,itt,ilmg,ilmh,il
      integer ip
      integer n0,pmaxx,nlm0
      parameter (n0=10, pmaxx=2, nlm0=(n0+1)**2)
      double precision ddot,xmrpp(3),xmrp0(3),pclus(3)
      double precision rsmc(0:n0),ehvlc(0:n0),xcheck,tol
      double precision gex(0:pmaxx,nlm0),ggrad(3,nlm0)
      double precision hex(0:1,nlm0)
C      double precision hs(0:n0),dhs(0:n0),ddhs(0:n0)
C      double precision hsp(0:n0),dhsp(0:n0),ddhsp(0:n0)
      data xmrp0/1d-5,2d-5,3d-5/,tol/1d-12/


C --- Offset to iax table for cluster connected to R
      offR  = ntab(iat)
      nclus = ntabg(iat)
C --- Offset to strux for this cluster
      isg = 0
      if (iat > 1) then
        do  it = 1, iat-1
          isg = isg + ntabg(it)
        enddo
      endif

      if (ivl == 0) then
        pmax = 2
      else
        pmax = 1
C       jhd = 10*(ivl-1)+2
      endif

C --- Checks ---
      if (nds > n0*n0)
     .  call rxi('strgck: nds is bigger than n0^2. nds = ',nds)
      if (nclus < 1)
     .  call rxi('strgck: empty cluster encountered. nclus = ',nclus)
      if (nclus > ntab(iat+1)-ntab(iat))
     .  call rxi('strgck: cluster is too big. nclus = ',nclus)

C ... sort, lmax, and coordinates/alat of the head of the cluster
      ib = iax(1,offR+1)
C     isi = ips(ib)
      nlmh = iax(9,offR+1)

      call dpzero(slj,nkap*nds)
      call dpzero(slk,nkap*nds)

C ... Loop over all atoms R" in the cluster
      do  it = 1, nclus
        itt = it + offR
        jb = iax(2,itt)
        isj = ips(jb)
C       Number of L channels at R" entering into screened function
        nlmg = iax(9,itt)
        call acoord(ib,it,alat,plat,pos,iax,ntab,pclus)

        lmaxg = ll(nlmg)
        do  il = 0, lmaxg
          rsmc(il) = rsm(il+1,isj)
        enddo
        if (ivl /= 0) then
          do  il = 0, lmaxg
            ehvlc(il) = ehvl(il+1,isj)
          enddo
        endif

c ... coordinates xmrp relative to atom it: xmrpp = x - R'
        do  ii = 1, 3
          xmrpp(ii) = xmrp(ii) - pclus(ii)
        enddo

c ... if xmrp accidentally hits a regular site, make a small offset
c     and print a warning
        xcheck = ddot(3,xmrpp,1,xmrpp,1)
        if (xcheck <= tol*tol) then
          do  ii = 1, 3
            xmrpp(ii) = xmrp0(ii)
          enddo
          call info5(10,1,0,
     .      ' strgck: Warning! |xmrpp(R,R'')|=%g for R=%i and R''=%i',
     .      dsqrt(xcheck),iat,it,0,0)
          call info2(10,0,1,
     .      '%11f Point is offset by vector xmrp0 = %3:1;3g',
     .      xmrp0,0)
        endif

c   ... make solid Gaussians at xmrpp
        call solgsg(xmrpp,rsmc,lmaxg,pmax,pmaxx,gex,ggrad)
        if (ivl == 0) then
          do  ip = 0, pmax-1
            do  ilmg = 1, nlmg
              hex(ip,ilmg) = gex(ip+1,ilmg)
            enddo
          enddo
c   ... or sm H, or sm hdot at xmrpp
        else
C         call solhsg(xmrpp,rsmc,ehvlc,lmaxg,ivl-1,fac,hex,ggrad)
          call sole0g(xmrpp,rsmc,ehvlc,lmaxg,ivl-1,1,0d0,hex,0d0)
C         call prmx('hex',hex,pmaxx+1,pmax+1,nlmg)
        endif

c ... combine with sg to evaluate the value-Laplacian set slj
c     and their Laplacians slk
        do  ilmg = 1, nlmg
          do  ikap = 1, nkap
            do  ilmh = 1, nlmh
              slj(ikap,ilmh) = slj(ikap,ilmh) +
     .                       gex(0,ilmg)*sg(ilmg,ilmh,1,ikap,isg+it) +
     .                       hex(0,ilmg)*sg(ilmg,ilmh,2,ikap,isg+it)
              slk(ikap,ilmh) = slk(ikap,ilmh) +
     .                       gex(1,ilmg)*sg(ilmg,ilmh,1,ikap,isg+it) +
     .                       hex(1,ilmg)*sg(ilmg,ilmh,2,ikap,isg+it)
            enddo
          enddo
        enddo

c ... end loop over cluster atoms
      enddo

      end
