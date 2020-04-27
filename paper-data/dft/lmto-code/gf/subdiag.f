      subroutine subdiag(kcplx,ldg,lidim,lhdim,ioff,isp,gma,s)
C- Subtract gma from the diagonal of S
Ci Inputs
Ci   kcplx :storage mode for gf
Ci   ldg   :leading dimension of gf
Ci   lidim :number of lower+intermediate orbitals to add
Ci   lhdim :dimensions gma
Ci   ioff  :offset for gma
Ci   isp   :spin index to gma
Ci   gma   :array of (gamma-alpha)
Cio Inputs/Outputs
Cio  s     :structure constant matrix, alpha representation on input
Cio        :gamma representation on output
Cr Remarks
Cr   This routine is used by the DLM branch, for which the crystal
Cr   Green's function can not be converted to gamma-rep
Cu Updates
Cu   05 Jan 13 Use spin-dependent gamma-alpha
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,ldg,lidim,lhdim,ioff,isp
      double precision gma(lhdim,2),s(ldg,ldg,2)
C ... Local parameters
      integer i

      if (kcplx /= 0) call rx('sa2gam: unsupported complex mode')

      do  i = 1, lidim
        s(i,i,1) = s(i,i,1) + gma(i+ioff,isp)
      enddo
      end
