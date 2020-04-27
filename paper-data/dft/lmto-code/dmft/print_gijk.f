      subroutine print_gijk(nomg,nsp,nbd_print,eps,mu,omega,uout)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This function writes G_ijk(iomega) on file gijk.ext
      ! Energy units are eV
       implicit none
       integer,    intent(in) :: nomg,nsp,nbd_print,uout
       complex(8), intent(in) :: eps(nbd_print,nomg,nsp)
       real(8),    intent(in) :: mu, omega(nomg)
       complex(8), parameter :: z1=cmplx(1d0,0d0)
       real(8), parameter :: ry2ev = 13.60569193d0
       complex(8) :: gijk
       integer    :: iomg,i,isp


       do iomg=1,nomg
        write(uout,'(f14.8)',advance='no') omega(iomg)*ry2eV
        do isp=1,nsp
         do i=1,nbd_print
          gijk = z1/(cmplx(mu,omega(iomg))-eps(i,iomg,isp)) ! 1/[ii*w+mu-eps(i,w,spin)]
          write(uout,'(2(xx,f14.8))',advance='no') dble(gijk)/ry2eV,aimag(gijk)/ry2eV
         enddo
        enddo
        write(uout,*)
       enddo
      end subroutine print_gijk
