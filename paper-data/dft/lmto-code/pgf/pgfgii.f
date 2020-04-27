      subroutine pgfgii(mode,s_site,s_ham,ipl,npl,pgplp,isp,ldg1,ldg2,gf)
C- Extract the site-diagonal blocks of gii for one layer, and record in s_site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gii
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   mode :1s digit pertains to storage
Ci        : 0  copy s_site(ib)%gii(1,isp) to gf
Ci        : 1  copy gf to s_site(ib)%gii(1,isp)
Ci        : 2  copy s_site(ib)%gc(1,isp) to gf
Ci        : 3  copy gf to s_site(ib)%gc(1,isp)
Ci        10s digit distinguishes how complex arithmetic is handled
Ci          0: gf has real, imaginary separated
Ci             gf = gf(ldg1,ldg2,2), with gf(*,*,1..2) = real..imag
Ci          1: gf is returned complex*16 format:
Ci             gf = gf(2,ldg1,ldg2), with gf(1..2,*,*) = real..imag
Ci          2: gf has real, imaginary separated by columns
Ci             gf = gf(ldg1,2,ldg2), with gf(*,1..2,*) = real..imag
Ci   ipl   :index to current principal layer
Ci          -1 => strux for left bulk PL; npl => strux for right bulk
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   isp   :current spin channel (1 or 2) --- affects which spin of s_site(ib)%gii(1,isp) is accessed
Ci   ldg1  :leading dimension of gf.  Noncollinear case gf=gf(ldg1,2,ldg2,2)
Ci   ldg2  :second  dimension of gf
Ci   gii   :q-summed g (g_RR' for R,R' inequivalent sites in the basis)
Ci         :Generated in gfibz
Cio Inputs/Outputs
Cio  gf    :site diagonal part copied into or out of s_site(ib)%gii(1,isp)
Cio        :for all sites in PL ipl
Cl Debugging
Cl   set outfile=coptcoll/out0.1 (for spin 1)
Cl   set ib = 1; mc -f9f12.6 $outfile -split g `mlist -sep:, 1:55:9` `mlist -sep:, 1:55:9` g$ib$ib  out.copt -- -px | head -3
Cl   Noncollinear case loop over ij blocks
Cl   set ib=1 outfile=coptnc/out0
Cl   foreach i ( 1 2 )
Cl   foreach j ( 1 2 )
Cl    echo block $i $j
Cl    mc -f9f12.6 $outfile -split gnc 1,nr/2+1,nr+1 1,nr/2+1,nr+1 gnc$i$j -split g `mlist -sep:, 1:55:9` `mlist -sep:, 1:55:9`
Cl   end
Cl   end
Cu Updates
Cu   18 Jun 14 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ipl,npl,isp,pgplp(6,-1:npl),ldg1,ldg2
      complex(8) gf(ldg1,*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iprmb(:) ! ,lo(:),hi(:)
      complex(8), pointer :: gsite(:,:)
C ... For structures
!       include 'structures.h'
      type(str_ham)::   s_ham
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ib1,ib2,nspc,lidim,mxorb,ldg1x,ldg2x,kcplx,off,ld,mod0
C     character*40 outs
C     procedure(logical) :: bittst
      procedure(integer) :: nglob

      call gtibpl(ipl,npl,pgplp,ib1,ib2)
      mxorb = nglob('mxorb'); nspc = nglob('nspc')
!       allocate(lo(ib1:ib2),hi(ib1:ib2))
C     lso   = bittst(s_ctrl%lncol,4)
C     lrel  = mod(s_ctrl%lrel,10)
      lidim = s_ham%ldham(2)
      ldg1x = ldg1*nspc
      ldg2x = ldg2*nspc
      mod0  = mod(mode,10)
      kcplx = mod(mode/10,10)
      call ztoyy(gf,ldg1x,ldg2x,ldg1x,ldg2x,kcplx,1)

C     call offsHpi(s_ham%iprmb,ib1,ib2,mxorb,0,lidim,lo,hi)
      call offsHp(s_ham%iprmb,ib1,ib2,0,0,lidim,off,ld)
      allocate(iprmb(mxorb))
      do  ib = ib1, ib2
        call icopy(mxorb,s_ham%iprmb(1+mxorb*(ib-1)),1,iprmb,1)
        iprmb(:) = iprmb(:) - off
C       print *, ib,ipl,minval(iprmb),maxval(iprmb),size(s_site(ib)%gii)
        gsite => s_site(ib)%gii; if (mod0 >= 2) gsite => s_site(ib)%gc
        call pokeg0(mod(mod0,2),s_site(ib)%norb,iprmb,ldg1,ldg2,lidim,nspc,nspc,gsite,gf)
C       if (isp == 1) cycle
C        call awrit2('site g(ib=%i) ipl=%i',outs,len(outs),0,ib,ipl)
C        call zprm(outs,2,s_site(ib)%gii(1,isp),
C     .    s_site(ib)%norb*nspc,s_site(ib)%norb*nspc,s_site(ib)%norb*nspc)
      enddo
!       deallocate(lo,hi)
      deallocate(iprmb)
      call ztoyy(gf,ldg1x,ldg2x,ldg1x,ldg2x,1,kcplx)

      end

