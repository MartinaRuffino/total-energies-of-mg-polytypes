      subroutine renorm_proj(dmftu, Olapm, nindo, ndham, nlmax, nicix)
      IMPLICIT NONE
      integer, intent(in) :: nicix,nlmax,ndham
      integer :: iorb1,nind1,ind1
      real(8) :: olocef
      complex(8), intent(inout) :: dmftu(ndham,nlmax,nicix)
      integer, intent(in) :: nindo(nicix)!,dmft(:,:,:)
      complex(8), intent(in) :: Olapm(nlmax,nlmax,nicix,nicix)
       do iorb1 = 1, nicix
        nind1 = nindo(iorb1)
         do ind1 = 1, nind1
          olocef = 1/sqrt(real(Olapm(ind1,ind1,iorb1,iorb1)))
          dmftu(:,ind1,iorb1) = dmftu(:,ind1,iorb1) * olocef
         enddo
       enddo
      end subroutine renorm_proj
