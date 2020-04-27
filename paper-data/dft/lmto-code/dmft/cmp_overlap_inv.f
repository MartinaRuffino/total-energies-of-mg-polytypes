      subroutine cmp_overlap_inv(olpinv, Olapm, nindo, nlmax, nicix)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   olpinv
Ci   Olapm
Ci   nindo
Ci   nlmax
Ci   nicix
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   19 Oct 14
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: nicix,nlmax
!      complex(8), intent(inout) :: dmftu(ndham,nlmax,nicix)
      integer, intent(in) :: nindo(nicix)!,dmft(:,:,:)
      complex(8), intent(in) :: Olapm(nlmax,nlmax,nicix,nicix)
      complex(8), intent(out) :: olpinv(nlmax*nicix,nlmax*nicix)
C ... Local parameters
      integer :: iorb,nind,nind1,info,i,j,bigl
!     real(8) :: olocef
      complex(8) :: olp(nlmax,nlmax), olpint(nlmax,nlmax)

      bigl = nlmax*nicix  ! This variable needs a more appropriate name
      olpinv = 0
      nind = 1
      do iorb = 1, nicix
       olp = Olapm(:,:,iorb,iorb)
       nind1 = nindo(iorb)+1
       forall (i = nind1:nlmax) olp(i,i) = cmplx(1.0d0, 0.0d0)
       forall (i = 1:bigl) olpinv(i,i) = cmplx(1.0d0, 0.0d0)

!       call print_mz(olp,'olpprint',label='olp')
!       call  zheev_paolo('olpeigen',olp_cp,nlmax)

       call zpotrf ('u', nlmax, olp, nlmax, info); if (info /= 0) stop 'zpotrf info /= 0'
       call ztrtri('u', 'n', nlmax, olp, nlmax, info); if (info /= 0) stop 'ztrtri info /= 0'
       do j = 1, nlmax-1
        do i = j+1, nlmax
            olp(i,j) = 0
        end do
       end do
!       dmftuu = dmftu(:,:,1)
       call zgemm('N', 'C', nlmax, nlmax, nlmax, (1.d0,0.d0), olp, nlmax, olp, nlmax, (0.d0,0.d0), olpint, nlmax)
       olpinv(nind:,nind:) = olpint(:,:)
       nind = nind1
      enddo
!       call print_mz(olpint,'olpint',mode='write')
!       call print_mz(olpinv,'olpinv',mode='write')

!       do iorb1=1,nicix
!        nind1 = nindo(iorb1)
!         do ind1=1,nind1
!          olocef = 1/sqrt(real(Olapm(ind1,ind1,iorb1,iorb1)))
!          dmftu(:,ind1,iorb1) = dmftu(:,ind1,iorb1) * olocef
!         enddo
!       enddo
      end subroutine cmp_overlap_inv
