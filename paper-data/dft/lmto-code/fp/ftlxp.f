      subroutine ftlxp(nbas,s_site,s_spec,alat,ng,gv,cv,k0,nlm0,fkl)
C- Pkl expansion around all sites of a function given as Fourier series.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl rsmv kmxv
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cv    :FT of function, as list of G-vector coefficients
Ci   k0    :leadingd dimension of fkl
Ci   nlm0  :second dimension of fkl
Co Outputs
Co   fkl   :coefficient to P_kL expansion of function
Cr Remarks
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   31 May 00 Adapted from nfp ftlxp.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k0,nbas,ng,nlm0
      double precision gv(ng,3),alat
      double complex cv(ng),fkl(0:k0,nlm0,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: yl(:)
      real(8), allocatable :: g2(:)
      complex(8), allocatable :: h(:)
C ... Local parameters
      integer ib,ilm,ipr,is,k,kmxv,l,lgunit,lmxl,ltop,m,nlm,ntop,stdo
      double precision a,fpi,pfac,pi,rsmv,tip,top,tpiba,df(0:20),
     .  fact(0:20),tau(3)
      double complex cfac

      stdo = lgunit(1)
      call getpr(ipr)
      call setfac(20,fact)
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      fpi = 4*pi
      tpiba = 2d0*pi/alat

C ... Set up spherical harmonics for max l needed
      ltop = -1
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        ltop = max0(ltop,lmxl)
      enddo
      ntop = (ltop+1)**2
      allocate(yl(ntop*ng),g2(ng),h(ng))
      call ropyln(ng,gv(1,1),gv(1,2),gv(1,3),ltop,ng,yl,g2)

C ... Scale g2 by (2pi/a)**2
      call dpcopy(g2,g2,1,ng,tpiba*tpiba)

C --- Start loop over atoms ---
      do  ib = 1, nbas
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        lmxl = s_spec(is)%lmxl
        rsmv = s_spec(is)%rsmv
        kmxv = s_spec(is)%kmxv
        if (lmxl == -1) goto 10

        nlm = (lmxl+1)**2
        if (nlm > nlm0) call rxi('ftlxp: nlm > nlm0, need',nlm)
        if (kmxv > k0)  call rxi('ftlxp: kmxv > k0, need',kmxv)
        if (ipr >= 60) write(stdo,200) ib,is,nlm,rsmv,lmxl,kmxv,tau
  200   format(' ib=',2i3,'  nlm=',i3,'  rsmv=',f6.3,
     .     '  lv,kv=',2i2,'  pos=',3f7.4/
     .     ' Expansion coefficients:')

        call ftlxp2(rsmv,tau,kmxv,k0,nlm,ng,gv,g2,yl,h,cv,
     .    fkl(0,1,ib))

C   ... Put in factors independent of G
        a = 1d0/rsmv
        cfac = (1d0,0d0)
        ilm = 0
        do  l = 0, lmxl
          do m = -l,l
            ilm = ilm+1
            do  k = 0, kmxv
              pfac = (4*a*a)**k * a**l * fact(k) * df(2*l+1)
              fkl(k,ilm,ib) = fkl(k,ilm,ib)*cfac*fpi/pfac
            enddo
          enddo
          cfac = cfac*dcmplx(0d0,tpiba)
        enddo

C   ... Printout
        if (ipr >= 60) then
          do  ilm = 1, nlm
            top = 0d0
            tip = 0d0
            do  k = 0, kmxv
              top = dmax1(top,dabs(dble(fkl(k,ilm,ib))))
              tip = dmax1(tip,dabs(dimag(fkl(k,ilm,ib))))
            enddo
            if (top > 1d-8)
     .         write(stdo,400) ilm,(dble(fkl(k,ilm,ib)),k = 0,kmxv)
            if (tip > 1d-8)
     .         write(stdo,401)     (dimag(fkl(k,ilm,ib)),k = 0,kmxv)
  400       format(i4,' RE',5f12.6:/(7x,5f12.6))
  401       format(4x,' IM',5f12.6:/(7x,5f12.6))
          enddo
        endif

   10   continue
      enddo

      deallocate(yl,g2,h)

      end

      subroutine ftlxp2(rsmv,tau,kmxv,k0,nlm,ng,gv,g2,yl,h,cv,fkl)
      implicit none
C ... Passed parameters
      integer k0,kmxv,ng,nlm
      double precision rsmv,yl(ng,*),tau(3),gv(ng,3),g2(ng)
      double complex h(ng),cv(ng),fkl(0:k0,nlm)
C ... Local parameters
      integer i,ilm,k
      double precision scalp,sum1,sum2,pi,gam

      pi = 4d0*datan(1d0)
      gam = 0.25d0*rsmv**2

C ... Init h with conjg(phase) * exponential factor * coeff
      do  i = 1, ng
        scalp = -2*pi*(tau(1)*gv(i,1)+tau(2)*gv(i,2)+tau(3)*gv(i,3))
        h(i) = dcmplx(dcos(scalp),-dsin(scalp))*dexp(-gam*g2(i))*cv(i)
      enddo

C ... Loop over k and ilm, do the scalar products with cv
      do  k = 0, kmxv
        do  ilm = 1, nlm
          sum1 = 0d0
          sum2 = 0d0
          do  i = 1, ng
            sum1 = sum1 + dble(h(i))*yl(i,ilm)
            sum2 = sum2 + dimag(h(i))*yl(i,ilm)
          enddo
          fkl(k,ilm) = dcmplx(sum1,sum2)
        enddo
C   ... Multiply h to get next k-th power of (-g2)
        do  i = 1, ng
          h(i) = -h(i)*g2(i)
        enddo
      enddo

      end
