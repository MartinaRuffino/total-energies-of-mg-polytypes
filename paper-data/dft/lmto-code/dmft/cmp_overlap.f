      subroutine cmp_overlap(dmftu,nindo,nicix,nlmax,ndham,nev,Olapmk)
      IMPLICIT NONE
      integer, intent(in) :: nicix,nlmax,ndham,nev
      complex(8) :: olp(nlmax,nlmax)
      complex(8), intent(in) :: dmftu(ndham,nlmax,nicix)
      integer, intent(in) :: nindo(nicix)!,dmft(:,:,:)
      complex(8), intent(out) :: Olapmk(nlmax,nlmax,nicix,nicix)
      integer :: iorb1,iorb2,nind1,nind2
      olp = 0
        do iorb1 = 1, nicix
         nind1 = nindo(iorb1)
         do iorb2 = 1, nicix
          nind2 = nindo(iorb2)
          call zgemm('C','N',nind1,nind2,nev,(1.d0,0.d0),dmftu(:,:,iorb1),ndham,dmftu(:,:,iorb2),ndham,(0.d0,0.d0),olp,nlmax)
!           do ind1 = 1, nind1
!            do ind2 = 1, nind2
              Olapmk(:,:,iorb1,iorb2) = olp!DIMITAR CHANGES(not conjg)
!            enddo
!           enddo
         enddo
        enddo

      end subroutine cmp_overlap
