      subroutine hhiml(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,
     .  ndim1,ndim2,cg,indxcg,jcg,cy,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1's digit (not implemented: always vectors)
Ci         :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
Ci   e1    :energies of smooth Hankels at p1 (l-dependent)
Ci   e2    :energies of smooth Hankels at p2 (l-dependent)
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordan coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   1. s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cr   Row L corresponds to p1 and col M corresponds to p2.
Cr   Strux s(L,M,k) is computed for dp=p1-p2
Cr   See JMP 39, 3393, Section 10
Cr
Cr   2. s(L,M,k) are also integrals of H_k1,L(r-p1) * H_k2,M(r-p2)
Cr   provided k = k1 + k2 (see MvS's notes on smooth hankels)
Cb Bugs
Cb   Not tested for kmax > 0 and l-dependent e/rsm
Cu Updates
Cu   06 Mar 09  Adapted from fp/hhibl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nlm1,nlm2,kmax,ndim1,ndim2
      double precision p1(3),p2(3),rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
      double precision s(ndim1,ndim2,0:kmax)
C ... C.G. coefficients
      integer jcg(*),indxcg(*)
      double precision cg(*),cy(*)
C ... Local parameters
      integer lmx1,lmx2,ll,l1,l2,lm11,lm21,lm12,lm22,l1t,l2t
      integer lmaxx,nlmx,il,mode0
      integer, parameter :: lmx0 = 12
      double precision dp(3),dp2,yl((lmx0+1)**2)

      if (nlm1 == 0 .or. nlm2 == 0) return
      mode0 = mod(mode,10)

      lmx1 = ll(nlm1)
      lmx2 = ll(nlm2)
      lmaxx = lmx1 + lmx2
      if (lmaxx > lmx0) call rxi('hhiml: increase lmx0 up to',lmaxx)

C ... Stop if positive energies encountered
      if (e1(0) > 0 .or. e2(0) > 0) call rx('hhiml: e > 0 not implemented')
      if (mode0 > 0) then
      do  il = 0, lmx1
        if (e1(il) > 0) call rx('hhiml: e > 0 not implemented')
      enddo
      do  il = 0, lmx2
        if (e2(il) > 0) call rx('hhiml: e > 0 not implemented')
      enddo
      endif

      dp(:) = p1(:) - p2(:)
C ... Make normalized sph. harmonics for connecting vector dp
      call sylm(dp,yl,lmaxx,dp2)
      nlmx = (lmaxx+1)**2
      yl(1:nlmx) = yl(1:nlmx)*cy(1:nlmx)

C ... Initialize s
      s(1:nlm1,1:nlm2,0:kmax) = 0d0

C ... Make integrals for each (l1,l2)
      l1t = -1
      do  l1 = 0, lmx1
        if (l1 <= l1t) cycle

        if (mode0 == 0) then
          l1t = lmx1
        else
          call gtbsl2(l1,lmx1,e1,rsm1,l1t)
        endif

        lm11 = l1**2+1
        lm12 = (l1t+1)**2
        l2t = -1
        do  l2 = 0, lmx2
          if (l2 <= l2t) cycle

          if (mode0 == 0) then
            l2t = lmx2
          else
            call gtbsl2(l2,lmx2,e2,rsm2,l2t)
          endif

          lm21 = l2**2+1
          lm22 = (l2t+1)**2
c         if (mode/10 == 1 .and. rsm1(l1)*rsm2(l2) == 0) cycle

C         If mode == 0, l1 and l2 are always zero
          call phhiml(dp2,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,
     .      lm21,lm22,kmax,ndim1,ndim2,yl,cg,indxcg,jcg,s)

        enddo
      enddo

      end
      
      subroutine phhiml(dp2,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax,
     .  ndim1,ndim2,yl,cg,indxcg,jcg,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   dp    :|p1-p2|^2
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of Hankels at p2
Ci   e1    :energy of smooth Hankels at p1
Ci   e2    :energy of smooth Hankels at p2
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   mlm1  :fill up s(L,M,k) for m1m1 =< L =< nlm1
Ci   mlm2  :fill up s(L,M,k) for m1m2 =< M =< nlm2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   yl    :Sph. harmonic polynomials for vector dp
Ci   cg    :Clebsch Gordan coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cu Updates
Cu   06 Mar 09  Adapted from fp/hhibl.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2
      double precision dp2,rsm1,rsm2,e1,e2
      double precision yl(*),s(ndim1,ndim2,0:kmax)
C ... C.G. coefficients
      integer jcg(*),indxcg(*)
      double precision cg(*)
C ... Local parameters
      logical le
      integer icg,icg1,icg2,ii,ilm,ilm1,ilm2,indx,ip,k,ik,
     .  ktop,l1,l2,ll,lm,lmax1,lmax2,lmaxx,nlmx,il
      integer, parameter :: lmx0 = 12, ktop0=10
      double precision fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx
      double precision xi1(0:lmx0,0:ktop0),xi2(0:lmx0,0:ktop0)
      double precision xi(-1:lmx0),fi(-1:lmx0)
      double precision gkl(0:ktop0,0:lmx0),beta1,beta2
      double precision, parameter :: tole = 1d-5, tolp = 1d-5

      fpi = 16d0*datan(1.d0)
      gam1 = 0.25d0*rsm1*rsm1
      gam2 = 0.25d0*rsm2*rsm2
      gamx = gam1+gam2
      rsmx = 2d0*dsqrt(gamx)
      lmax1 = ll(nlm1)
      lmax2 = ll(nlm2)
      lmaxx = lmax1+lmax2
      nlmx = (lmaxx+1)**2
      ktop = max0(lmax1,lmax2)+kmax
      le = dabs(e1-e2) > tole

      if (lmaxx > lmx0) call rxi('hhiml: increase lmx0 up to',lmaxx)
      if (ktop > ktop0) call rxi('hhiml: increase ktop0 up to',ktop)

C ... Make energy independent Gaussians
      if (ktop > 0) call radgkv(1,ktop0,dp2,rsmx,1,ktop-1,lmaxx,gkl)

C ... Set up functions for connecting vector dp = p1-p2
      if (le) then                                              ! case e1 /= e2
         fac1 = dexp(gam2*(e2-e1))/(e1-e2)
         fac2 = dexp(gam1*(e1-e2))/(e2-e1)
         beta1 = fpi*dexp(e1*gamx)
         beta2 = fpi*dexp(e2*gamx)
         call hansrz(-rsmx,0,lmaxx,e1,dp2,1,1,10,xi1(0,0),fi(0))
         call hansrz(-rsmx,0,lmaxx,e2,dp2,1,1,10,xi2(0,0),fi(0))
c ... Make Laplacians H_k+1,L = -e*H_kL - 4*pi*G_kL
         if (ktop > 0) then
           do  ik = 1, ktop
             xi1(0:lmaxx,ik) =
     .          -e1*xi1(0:lmaxx,ik-1) - beta1*gkl(ik-1,0:lmaxx)
             xi2(0:lmaxx,ik) =
     .          -e2*xi2(0:lmaxx,ik-1) - beta2*gkl(ik-1,0:lmaxx)
           enddo
         endif
         do  ik = 0, ktop
           xi1(0:lmaxx,ik) = fac1*xi1(0:lmaxx,ik) + fac2*xi2(0:lmaxx,ik)
         enddo
      else                                                     ! case e1 = e2
        e = .5d0*(e1+e2)
        if (dabs(e) > tole) then
          call hansrz(-rsmx,-1,lmaxx,e,dp2,1,1,10,xi,fi)
C ...     xi2 contains h(r)/r^l and xi1 contains w(r)/r^l
C           where w(r) = hdot(r) - gamma*h(r) [JMP](7.12)
C           and hdot_l(r)/r^l = 1/2 xi(l-1) [JMP](7.5)
          xi2(0:lmaxx,0) = xi(0:lmaxx)
          do  il = 0, lmaxx
            xi1(il,0) = 0.5d0*xi(il-1) - gamx*xi(il)
          enddo
c ...     make h_pl in xi2 and w_pl in xi1 as:
C           h_p+1,l = -e*h_pl - 4*pi*G_pL [JMP]:(6.34)
c           w_p+1,l = -e*w_pl - H_pL      [JMP]:(7.13)
          if (ktop > 0) then
            beta1 = fpi*dexp(e*gamx)
            do  ik = 1, ktop
              xi2(0:lmaxx,ik) =
     .           -e*xi2(0:lmaxx,ik-1) - beta1*gkl(ik-1,0:lmaxx)
              xi1(0:lmaxx,ik) =
     .           -e*xi1(0:lmaxx,ik-1) - xi2(0:lmaxx,ik-1)
            enddo
          endif
        else
          call rx('hhiml: integral diverges if e1 = e2 = 0')
        endif
      endif

C ... Combine with Clebsch-Gordan coefficients and multiply by Y_L
      do  ilm1 = mlm1, nlm1
        l1 = ll(ilm1)
        do  ilm2 = mlm2, nlm2
          l2 = ll(ilm2)
          ii = max0(ilm1,ilm2)
          indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            ilm = jcg(icg)
            lm = ll(ilm)
            k = (l1+l2-lm)/2
            fac = fpi*(-1d0)**l1*cg(icg)*yl(ilm)
            do  ip = 0, kmax
              s(ilm1,ilm2,ip) = s(ilm1,ilm2,ip) +
     .          fac*xi1(lm,k+ip)
            enddo
          enddo
        enddo
      enddo

      end
