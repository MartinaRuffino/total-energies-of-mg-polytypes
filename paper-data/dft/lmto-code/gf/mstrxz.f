      subroutine mstrxz(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,iop,s,sdot)
C- Structure constants to h, for complex energy
C ----------------------------------------------------------------
Ci Inputs:
Ci   e:   energy of Hankel function
Ci   dr:  vector R_source - R_expansion
Ci   nlma:Generate coefficients S_R'L',RL for L' < nlma
Ci   nlmb:Generate coefficients S_R'L',RL for L  < nlmb
Ci   ndim:row dimension of structure constant matrix s
Ci   cy  :normalization constants for spherical harmonics (sylmnc.f)
Ci   cg,indxcg,jcg:Clebsch Gordan coefficients, and indices (scg.f)
Ci   iop :1's digit
Ci        1: calculate s only
Ci        2: calculate sdot only
Ci        any other number: calculate both s and sdot
Ci       10's digit: nonzero => use OKA's conventions
Ci        1 strux to satisfy Andersen's definitions, with factor (-1)
Ci      100s digit:
Ci        1 strux to one-center expansions of Neumann functions.
Co Outputs
Co   s   :structure constant matrix S_R'L',RL
Co   sdot:Energy derivative of s
Cb Bugs
Cb   sdot has not been checked.
Cr Remarks
Cr
Cr  Expansion Theorem: H_RL(r) = H_L(r-R)
Cr
Cr    H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
Cr                                          l   (l+l'-l")/2
Cr    S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1) (-E)         H_L"(E,R-R')
Cr
Cr    H_L  =  (i k)^(l+1) h_l (kr) Y_L   ?
Cr    J_L  =  (i k)^(-l)  j_l (kr) Y_L   ?
Cr
Cr  where  h = exp(ikr)/r, j = sin(kr)/r, k = sqrt(e), Im k >=0.
Cr  Generate w/ H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)
Cr
Cr  These satisfy:
Cr    (1)  cases e /= 0 and e == 0 same equations
Cr    (3)  S_KL (dr)  =  S^*_LK(-dr)
Cr    (4)  no extra factors (2l-1)!! etc as in Andersen's definitions.
Cr
Cr  Indices for L: middle is m=0; (z axis for m=0)
Cr                 after is real part of spherical harmonics: m=1,2, ...
Cr                 before is imaginary part m=1,2, ...
Cr  For s,p,d the indices 1-9 correspond to:
Cr         1     2     3     4     5     6     7     8     9
Cr         s     y     z     x     xy    yz  3z^2-1  xz  x^2-y^2
Cr
Cr  Andersen's definitions:
Cr  HOKA_L(k^2,r)  =  H_L / (2l-1)!!
Cr  JOKA_L(k^2,r)  = -J_L * (2l-1)!!/2
Cr
Cr  Thus Andersen structure constants differ by the ratio
Cr    S_RL,R'L' (OKA) / S_RL,R'L' =  2/((2l-1)!!(2l'-1)!!)
Cr  NB: Andersen defines H = -J S; mstrxz returns -S (OKA)
Cr
Cr  Note: the e-> zero limit of the ss(sigma) sdot diverges as
Cr  1/sqrt(e).  This singularity is removed for case e=0, making
Cr  a discontinuity in Sdot at e=0.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ndim,nlma,nlmb,iop
      integer indxcg(1),jcg(*)
      double complex e,s(ndim,nlmb),sdot(ndim,nlmb)
      double precision cy(*),cg(*),dr(3)
C Local parameters
      logical loka
      integer lmxx
      parameter (lmxx=11)
      integer icg,icg1,icg2,ii,mlm,indx,ipow,klm,l,lk,llm,lm,lmax,
     .        m,ll,klm0,sig(0:lmxx),iiop
      double complex psiz(-1:lmxx),phiz(-1:lmxx),gdot,sum,sumd
      double complex hz((lmxx+1)**2),hdtz((lmxx+1)**2)
      double complex edot(0:lmxx),efac(0:lmxx)
      double precision fpi,rfac,r2,xx,fac2l(0:lmxx),yl((lmxx+1)**2)
      external ll,sylm

C ... Setup
      lmax = ll(nlma) + ll(nlmb)
      fac2l(0) = 1
      do  5  l = 1, lmax
    5 fac2l(l) = fac2l(l-1) * (2*l-1)
      iiop = mod(iop,10)
      loka = mod(iop/10,10) /= 0
      fpi = 16.d0*datan(1.d0)
      lmax = ll(nlma) + ll(nlmb)
      if (nlma > ndim)
     .  call rxi('MSTRXZ: nlma gt ndim: need ndim=',nlma)
      if (lmax > lmxx)
     .  call rxi('MSTRXZ: lmax too big, need lmax=',lmax)
      call sylm(dr,yl,lmax,r2)
      if (r2 < 1d-10) return
C ... Hankel and Bessel functions, standard definitions
      call besslz(e*r2,(iop/100)*10+0,-1,lmax,phiz,psiz)
      mlm = 0
      rfac = dsqrt(r2)
C --- Generate H_L(DR) = 1/r^(2l+1) Psiz_l  Y_lm ---
C     and psidot, defined so that Hdot_L = H_L * psidot(l)
      do  10  l = 0, lmax
        gdot = psiz(l-1)*r2/2
        rfac = rfac/r2
        do  12  m = -l, l
          mlm = mlm+1
          xx = rfac*cy(mlm)*yl(mlm)
          hz(mlm)   = psiz(l)*xx
          hdtz(mlm) = gdot * xx
   12   continue
   10 continue

C ... efac(l)=(-e)**l; sig(l)=(-)**l * sign of iop
      efac(0) = 1
      edot(0) = 0
      sig(0)  = 1
      do  1  l = 1, lmax
        efac(l) = -e*efac(l-1)
        edot(l) = -l*efac(l-1)
        sig(l)  = -sig(l-1)
    1 continue

      if (iiop /= 2) then

C --- Generate S_MK ---
      do  11  mlm = 1, nlma
        lm = ll(mlm)
        klm0 = mlm
C   ... Handles subblock nlmb < mlm <= nlma, if there is one
        if (mlm > nlmb) klm0 = 1
        do  14  klm = klm0, nlmb
          lk = ll(klm)
          sum = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do  15  icg = icg1, icg2
            llm  = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            sum  = sum + cg(icg)*efac(ipow)*hz(llm)
   15     continue
          s(mlm,klm) = fpi*sig(lk)*sum
          if (loka) s(mlm,klm) = s(mlm,klm)*2/(fac2l(lm)*fac2l(lk))
          if (klm <= nlma .and. mlm <= nlmb)
     .      s(klm,mlm) = s(mlm,klm)*sig(lk)*sig(lm)
   14 continue
   11 continue

      endif

      if (iiop /= 1) then
C --- Generate Sdot_MK ---
      do  31  mlm = 1, nlma
        lm = ll(mlm)
        klm0 = mlm
C   ... Handles subblock nlmb < mlm <= nlma, if there is one
        if (mlm > nlmb) klm0 = 1
        do  31  klm = klm0, nlmb
        lk = ll(klm)
        sumd = 0d0
        ii = max0(mlm,klm)
        indx = (ii*(ii-1))/2 + min0(mlm,klm)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1) - 1
        do  32  icg = icg1, icg2
          llm = jcg(icg)
          ipow = (lm+lk-ll(llm))/2
          sumd = sumd + cg(icg)*
     .                  (efac(ipow)*hdtz(llm)+edot(ipow)*hz(llm))
   32   continue
        sdot(mlm,klm) = fpi*sumd*sig(lk)
        if (loka) sdot(mlm,klm) = sdot(mlm,klm)*2/(fac2l(lm)*fac2l(lk))
        if (klm <= nlma .and. mlm <= nlmb)
     .  sdot(klm,mlm) = sdot(mlm,klm)*sig(lk)*sig(lm)
   31 continue
      endif

      end
