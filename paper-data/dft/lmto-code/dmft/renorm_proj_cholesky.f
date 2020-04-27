      subroutine renorm_proj_cholesky(dmftu, Olapm, nindo, ndham, nev, nlmax, nicix)
C- Renormalize projectors.  This code needs a bit of cleaning up
C     use mod_paolo
      IMPLICIT NONE
      integer, intent(in) :: nicix,nlmax,ndham,nev
      complex(8), intent(inout) :: dmftu(ndham,nlmax,nicix)
      integer, intent(in) :: nindo(nicix)!,dmft(:,:,:)
      complex(8), intent(in) :: Olapm(nlmax,nlmax,nicix,nicix)
      integer :: nind1,info,i,j
!     real(8) :: olocef
      complex(8) :: olp(nlmax,nlmax),dmftuu(ndham,nlmax), olp_cp(nlmax,nlmax)
      logical, parameter :: debug=.true.

      nind1 = nindo(nicix) + 1
      olp = Olapm(:,:,nicix,nicix)
      forall (i = nind1:nlmax) olp(i,i) = cmplx(1d0, 0d0)

      olp_cp = olp

      if (debug) then
        call zwrmdmft(0,'olp','olpprint',' ',olp,nlmax,nlmax,nlmax)
      endif

      call zpotrf ('u', nlmax, olp, nlmax, info); if (info /= 0) stop 'zpotrf info /= 0'
      call ztrtri('u', 'n', nlmax, olp, nlmax, info); if (info /= 0) stop 'ztrtri info /= 0'
      do j = 1, nlmax-1
       do i = j+1, nlmax
           olp(i,j) = 0
       end do
      end do
      dmftuu = dmftu(:,:,1)
      call zgemm('N', 'N', nev, nlmax, nlmax, (1.d0,0.d0), dmftuu, ndham, olp, nlmax, (0.d0,0.d0), dmftu(:,:,1), ndham)

!       do iorb1=1,nicix
!        nind1 = nindo(iorb1)
!         do ind1=1,nind1
!          olocef = 1/sqrt(real(Olapm(ind1,ind1,iorb1,iorb1)))
!          dmftu(:,ind1,iorb1) = dmftu(:,ind1,iorb1) * olocef
!         enddo
!       enddo
      end subroutine renorm_proj_cholesky
