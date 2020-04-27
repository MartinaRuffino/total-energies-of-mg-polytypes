       subroutine print_sigf0(nev,nsp,sigbar,iqfbz,uout)
        ! the output is used by the flag --makesigqp which
        ! automatically converts from eV to Ry.
        ! So print_sigf0 writes in eV.
        implicit none
        real, parameter :: ry2ev = 13.60569193d0
        integer,    intent(in)  :: uout
        integer,    intent(in)  :: nev, nsp, iqfbz
        complex(8), intent(in)  :: sigbar(nev,nev,nsp)
        complex(8) :: tmp
        integer :: isp,i,j

        write(uout,'(" ")')
        write(uout,'("# KPT = ",i9)') iqfbz
        do isp=1,nsp
         write(uout,'("# SPIN = ",i4)') isp
         do i=1,nev
          do j=1,nev
           tmp = 0.5*(sigbar(i,j,isp)+conjg(sigbar(j,i,isp)))     ! hermitianize
           write(uout,'(f14.8)',advance='no') dble(tmp)*ry2ev
          enddo
          write(uout,*)
         enddo
        enddo



       end subroutine print_sigf0
