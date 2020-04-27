      subroutine rmeshprm(opt,s_spec,nrmx,is1,is2,intopt,a,b,nr,rmt,rofi,rwgt)
C- Retrieve radial mesh parameters for a ranged of species from structure
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   opt   :1s digit
Ci         :1 return a,b,rmt
Ci         :2 return return rofi
Ci         :4 return return rwgt
Ci         :Any combination is allowed
Ci         :10s digit
Ci         :Use input intopt
Ci         :Make intopt from intopt = 10*nglob('lrquad')
Ci  is1,is2:Range of species
Cio Inputs/Outputs
Cio intopt :Parameter controlling how quadrature weights are set; see radwgt
Cio        :Input if 10s digit 0; otherwise output
Co Outputs
Co  a      :Parameter defining radial mesh; see radmsh.f
Co         :Returned if 1s bit opt nonzero
Co  b      :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Co         :Returned if 1s bit opt nonzero
Co  rmt    :augmentation radius, in a.u., by species
Co         :Returned if 1s bit opt nonzero
Co  nr     :Number of radial mesh points
Co         :Returned if 1s bit opt nonzero
Co  rofi   :radial mesh points
Co         :Returned if 2s bit opt nonzero
Co  rwgt   :radial mesh weights
Co         :Returned if 4s bit opt nonzero
Cu Updates
Cu   27 Aug 18 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer opt,is1,is2,intopt,nrmx,nr(is1:is2)
      double precision a(is1:is2),b(is1:is2),rmt(is1:is2)
      double precision rofi(nrmx,is1:is2),rwgt(nrmx,is1:is2)
C ... Local parameters
      integer is,opt0
      procedure(integer) :: nglob

      opt0 = mod(opt,10)
      if (mod(opt/10,10) > 0) intopt = 10*nglob('lrquad')

      do  is = is1, is2
        if (s_spec(is)%lmxa < 0) cycle
        if (mod(opt0,2) > 0) then
          a(is) = s_spec(is)%a
          b(is) = s_spec(is)%rmt / (dexp(a(is)*s_spec(is)%nr-a(is))-1d0)
          rmt(is) = s_spec(is)%rmt
          nr(is) = s_spec(is)%nr
        endif

        if (mod(opt0,10) < 2) cycle

        if (s_spec(is)%nr > nrmx) call rx('mmeshprm: increase nrmx')
        call radmsh(s_spec(is)%rmt,s_spec(is)%a,s_spec(is)%nr,rofi)
        if (mod(opt0,10) < 4) cycle
        call radwgt(intopt,s_spec(is)%rmt,s_spec(is)%a,s_spec(is)%nr,rwgt)

      enddo
      end
