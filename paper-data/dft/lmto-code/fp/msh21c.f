      subroutine msh21c(ib,s_site,s_spec,s_lat,fncv,ng,f1cc,f1cr,nrx)
C- One-center expansion at one site of a function tabulated on a mesh
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rmt lmxl a nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc kv gv ng
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   ib    :One-center exansion around site ib
Ci         :ib<0 => sign used as a flag to return f1 as complex f1cc
Ci   fn    :F.T. of function for which to make one-center expansion
Ci   nrx   :leading dimension of f1cr,f1cc; must be at least as large as nr
Co Outputs
Co   f1cr  : The one-center expansion of fn is tabulated on a radial mesh
Co         : f1cr is real
Co   f1cc  : same as f1cr, but f1cc is complex
Cr Remarks
Cr   A PW exp(i.q.r) has a one-center expansion at radius r
Cr      sum_L C_L Y_L(r) where C_L = 4 pi i^l j_l(|rq|) Y_L(q)
Cr   Adapted from mshvmt.f
Cr   This routine is not optimized.
Cu Updates
Cu   16 Jun 16 Redesigned
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ib,nrx,ng
      double precision f1cr(nrx,*)
      double complex fncv(ng),f1cc(nrx,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays and pointers
      real(8), pointer :: gv(:,:)
      real(8),allocatable:: yl(:,:),gv2(:,:),
     .  agv(:),cgp(:),sgp(:),phil(:,:),rofi(:),rwgt(:)
C ... Local parameters
      logical lcmplx
      integer i,ibl,ilm,intopt,ipr,ir,is,l,lmxl,m,n1,n2,n3,nlm,nr,ngabc(3)
      double precision a,alat,fac,pi,r,rmt,tpiba
      double precision tau(3),plat(3,3)
      double complex vvali,fprli
      procedure(integer) :: nglob
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))

      call getpr(ipr)
      ibl = iabs(ib)
      lcmplx = ib < 0
      pi = 4d0*datan(1d0)
      alat = s_lat%alat
      plat = s_lat%plat
      ngabc = s_lat%nabc
      is = s_site(ibl)%spec
      tau = s_site(ibl)%pos
      call dscal(3,alat,tau,1)
      rmt = s_spec(is)%rmt
      lmxl = s_spec(is)%lmxl
      a = s_spec(is)%a
      nr = s_spec(is)%nr
      rmt = s_spec(is)%rmt
      nlm = (lmxl+1)**2
      tpiba = 2*pi/alat
      gv => s_lat%gv
      ng = s_lat%ng
      allocate(gv2(ng,3),agv(ng),cgp(ng),sgp(ng),phil(ng,0:lmxl),yl(ng,nlm))
      allocate(rofi(nr),rwgt(nr))
      intopt = 10*nglob('lrquad')

C --- YL(G)*G**l, agv=|G| for each g ---
      call dpcopy(gv,gv2,1,3*ng,tpiba)
      call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxl,ng,yl,agv)
      do  i = 1, ng
        agv(i) = sqrt(agv(i))
      enddo
      call radmsh(rmt,a,nr,rofi)
      call radwgt(intopt,rmt,a,nr,rwgt)

C --- For each radius, do ---
      if (lcmplx) then
        call dpzero(f1cc,2*nrx*nlm)
      else
        call dpzero(f1cr,nrx*nlm)
      endif
      do  ir = 1, nr
        r = rofi(ir)

C   ... j_l(|r*q|)/r**l for each G and l=0..lmax
C       Does not evolve correctly in the correct large r limit
        call ropbes(agv,r**2,lmxl,cgp,sgp,phil,ng,1)
C   ... Phases calculated straightforwardly.  Fast enough not to matter.
        do  i = 1, ng
          fac = -(tau(1)*gv2(i,1)+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
          cgp(i) = dcos(fac)
          sgp(i) = dsin(fac)
        enddo

C   ... Sum_G 4*pi*(i*r)**l j_l(|r*G|)/(r*G)**l YL(G) G**l
        ilm = 0
        fprli = 4*pi
        do  l  = 0, lmxl
          do  m = -l, l
            ilm = ilm+1
            vvali = 0
            do  i = 1, ng
              vvali = vvali + (phil(i,l)*yl(i,ilm))*
     .                        (fncv(i)*dcmplx(cgp(i),-sgp(i)))
            enddo
            if (lcmplx) then
              f1cc(ir,ilm) = fprli*vvali
            else
              f1cr(ir,ilm) = fprli*vvali
            endif
          enddo
          fprli = fprli*(0d0,1d0)*r
        enddo
      enddo

C     call prrmsh('one-center expansion',rofi,f1cr,nrx,nr,nlm)

      deallocate(gv2,agv,cgp,sgp,phil,yl,rofi,rwgt)
      end
