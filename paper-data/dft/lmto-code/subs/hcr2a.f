      subroutine hcr2a(s_spec)
C- Set hcr according to alpha, E=0
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa alpha hcr
Co     Stored:    hcr
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Co Outputs
Cl Local variables
Cl   lalp  :T, generate hcr corresponding to species alpha
Cl   kap2  :interstitial kinetic energy
Cr Remarks
Cr
Cu Updates
Cu   06 Sep 11 Started migration to f90 structures
Cu   08 Aug 06 Redesigned
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical lalp
      integer lmxx,nspec,is,lmxa,l,stdo,nl,ipr
      parameter(lmxx=20)
      double precision alpha(lmxx),adot(lmxx),hcr(lmxx),kap2
      double precision fi(-1:lmxx),gi(-1:lmxx),avw,dglob
      procedure(integer) :: nglob

      stdo = nglob('stdo')
      call getpr(ipr)
      nspec = nglob('nspec')
      nl = nglob('nl')
      avw = dglob('avw',0d0,0)
      lalp = .true.
      kap2 = 0

      if (ipr >= 50) then
        call info0(50,1,0,
     .    ' Hard core radii corresponding to 2nd generation alpha''s')
        call awrit1(
     .  ' spec   energy    hcr ...%npalpha/adot',
     .  ' ',120,stdo,19+10*nl)
      endif
      do  is = 1, nspec

        lmxa = s_spec(is)%lmxa
        if (lmxa >= lmxx) call rx('hcr2a: increase lmxx')
        if (lalp) then
          alpha(1:lmxa+1) = s_spec(is)%alpha(1:lmxa+1)
          call bessl2(0d0,-1,lmxa+1,fi(-1),gi(-1))
          do  l = 0, lmxa
            hcr(l+1) = avw*(alpha(l+1)/fi(l)*gi(l))**(1d0/dble(2*l+1))
          enddo
          s_spec(is)%hcr(1:lmxa+1) = hcr(1:lmxa+1)
        endif
        hcr(1:lmxa+1)  = s_spec(is)%hcr(1:lmxa+1)
C       print *, '!!';hcr = 1.61d0

C      This will make alphas for any energy at this hcr(l) (not used)
        do  l = 0, lmxa
          call bessl2(kap2*hcr(l+1)**2,-1,l+1,fi(-1),gi(-1))
          alpha(l+1) = fi(l)/gi(l)*(hcr(l+1)/avw)**(l+l+1)
          adot(l+1)  = -0.5d0*(hcr(l+1)/avw)**(l+l+3)*
     .          (fi(l+1)*gi(l)/(l+l+1)+fi(l)*gi(l-1)/(l+l-1))/
     .          (gi(l)*gi(l))
        enddo

C       ik = 1
c       call dscal(lmxa+1,1/avw,hcr,1)
        call info8(50,0,0,
     .    '%,4i %5p%;10,6D '//
     .    '%n;10,6D'//
     .    '%47p%n;10,6D',
     .    is,kap2,lmxa+1,hcr,lmxa+1,alpha,0,0)
        call info2(50,0,0,'%47p%n;10,6D',lmxa+1,adot)

      enddo
      end
