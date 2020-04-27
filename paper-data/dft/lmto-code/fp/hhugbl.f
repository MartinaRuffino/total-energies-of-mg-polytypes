      subroutine hhugbl(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,ndim1,
     .  ndim2,s_lat,wk,dwk,s,ds)
C- Estatic energy integrals between Bloch Hankels, and gradients.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  vol alat plat qlat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  hhigbl phhigb hklbl gklbl fklbl hsmbl
Ci Inputs
Ci   mode  :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radius of Hankels at p2 (l-dependent)
Ci   e1    :energy  of Hankels at p1 (l-dependent)
Ci   e2    :energy  of Hankels at p2 (l-dependent)
Ci   nlm1  :L-max for  Hankels at p1
Ci   nlm2  :L-max for  Hankels at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   slat  :struct containing information about the lattice
Ci   wk    :work space of same size as s
Ci   dwk   :work space of same size as ds
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
Cu   22 Apr 00 Adapted from nfp hhug_bl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlm1,nlm2,ndim1,ndim2
      double precision rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*),
     .  p1(3),p2(3)
      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3),
     .   wk(ndim1,ndim2),dwk(ndim1,ndim2,3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer kmax,kdim,i2,i1,lmxx,nlmx,l,ll
      parameter (lmxx=6,nlmx=(lmxx+1)**2)
      double precision add,pi,q(3),fpi,vol,gam1,gam2,xx1,xx2
      double precision zer(0:lmxx),bet1(nlmx),fac(nlmx)
      data q /0d0,0d0,0d0/
      data zer /lmxx*0d0,0d0/

      pi = 4d0*datan(1d0)
      vol = s_lat%vol
      kmax = 0
      kdim = 0
      if (nlm1 > nlmx) call rx('hhugbl: increase lmxx')

      call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,ndim1,
     .  ndim2,kdim,s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy,s_lat,s,ds)

      call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,zer,nlm1,nlm2,kmax,ndim1,
     .  ndim2,kdim,s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy,s_lat,wk,
     .  dwk)

      do  i2 = 1, nlm2
        if (mode == 0) then
          l = 0
        else
          l = ll(i2)
        endif
        bet1(i2) = dexp(e2(l)*rsm2(l)*rsm2(l)/4d0)
        fac(i2) = 8d0*pi/e2(l)
      enddo
      do  i2 = 1, nlm2
        do  i1 = 1, nlm1
          s(i1,i2)    = fac(i2)*(s(i1,i2)    - bet1(i2)*wk(i1,i2))
          ds(i1,i2,1) = fac(i2)*(ds(i1,i2,1) - bet1(i2)*dwk(i1,i2,1))
          ds(i1,i2,2) = fac(i2)*(ds(i1,i2,2) - bet1(i2)*dwk(i1,i2,2))
          ds(i1,i2,3) = fac(i2)*(ds(i1,i2,3) - bet1(i2)*dwk(i1,i2,3))
        enddo
      enddo

C ... Extra term for l1=l2=0
      fpi = 4d0*pi
      gam1 = 0.25d0*rsm1(0)*rsm1(0)
      gam2 = 0.25d0*rsm2(0)*rsm2(0)
      xx1 = fpi*dexp(gam1*e1(0))/(vol*e1(0))
      xx2 = fpi*dexp(gam2*e2(0))/(vol*e2(0))
      add = -2*xx1*xx2*vol/e2(0)
      s(1,1) = s(1,1) + add

      end
