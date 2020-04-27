        subroutine compressproj(s_dmft,dmftu,nlo,nhi,ldcix,nsp,nzsig,strans)
C- Full projector, all
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_dmft
Ci     Elts read:  iasig
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   dmftu :factored form of projector
Ci   ndham :dimensioning parameter, largest hamiltonian dimension
Ci   nev   :Number of eigenvectors to be used in projector
Ci   ldcix :dimensions dmftu
Ci   nzsig :Total number of nonzero matrix elements in all cix and spin blocks
Co Outputs
Co   strans: strans(ichan,ivec,jvec) = projector onto (iv,jv) pair for element ichan
Co         : Each element ichan as a definite spin isp and cix block
Co           strans(ichan,ivec,jvec) = dmftu^+(jv,mj,isp,cix)*dmftu(iv,mi,isp,cix)
Co           where mi and mj point to a particular (row,column) index in block cix.
Cr Remarks
Cr  See Eq 19 in arXiv:0907.0195v2
Cr  Note: in the spin polarized case, strans is calculated for both spins.
Cr  For collinear spins only on or the other spin sector will be used.
Cu Updates
Cu   16 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu      Nov 14 First created by P Pisanti.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, intent(in)     :: nlo,nhi,ldcix,nzsig,nsp
      complex(8), intent(out) :: strans(nzsig,nlo:nhi,nlo:nhi)
      complex(8), intent(in)  :: dmftu(nlo:nhi,ldcix,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_dmft)::  s_dmft
C ... Local parameters
      integer :: iv,jv,cix,mi,mj,isp,ichan

      call dpzero(strans,2*size(strans))
      do  ichan = 1, nzsig
        mi  = s_dmft%iasig(1,ichan)
        mj  = s_dmft%iasig(2,ichan)
        isp = s_dmft%iasig(3,ichan)
        cix = s_dmft%iasig(5,ichan)
C       Projectors do not flip spins with AFM mapping, only potential
C       if (s_dmft%icix(cix) < 0 .and. nsp == 2) isp = nsp+1-isp
        do  iv = nlo, nhi
          do  jv = nlo, nhi
            strans(ichan,iv,jv) = strans(ichan,iv,jv) +
     .        dconjg(dmftu(jv,mj,isp,cix))*dmftu(iv,mi,isp,cix)
          enddo
        enddo
      enddo

C     print *, 'checksum strans',sum(strans)

      end subroutine compressproj
