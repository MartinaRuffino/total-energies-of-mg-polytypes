        subroutine testdenmat(nbd,M,A,w,isp,iq)
         implicit none
         integer, intent(in)    :: nbd,isp,iq
         complex(8), intent(in) :: M(nbd,nbd),w(nbd),A(nbd,nbd)
         complex(8) :: wmat(nbd,nbd),Minv(nbd,nbd),tmp(nbd,nbd),R(nbd,nbd)
         complex(8) :: cwork(4*nbd)
         integer :: i,j,ipiv(nbd),ierr,u=9999


         write(u,*) " "
         write(u,'(" IQ, isp = ",2i6)')iq,isp

         ! construct diagonal matrix wmat
         wmat=cmplx(0d0,0d0)
         do i=1,nbd
          wmat(i,i)=w(i)
         enddo

         ! invert M
         Minv=M
         call zgetrf(nbd,nbd,Minv,nbd,ipiv,ierr)
         call zgetri(nbd,Minv,nbd,ipiv,cwork,4*nbd,ierr)


         ! control printing
         write(u,'("----- W -----")')
         do i=1,10
          write(u,'(10(3x,2f12.6))') (A(i,j) , j=1,10)
         enddo

         write(u,'("----- R (first 10x10)-----")')
         do i=1,10
          write(u,'(10(3x,2f12.6))') (M(i,j) , j=1,10)
         enddo

         write(u,'("----- inv(R) (first 10x10)-----")')
         do i=1,10
          write(u,'(10(3x,2f12.6))') (Minv(i,j) , j=1,10)
         enddo

         write(u,'("----- R*inv(R) (first 10x10)-----")')
         tmp=cmplx(0d0,0d0)
         call zgemm('n','n',nbd,nbd,nbd,(1d0,0d0),M,nbd,Minv,nbd,(0d0,0d0),tmp,nbd)
         do i=1,10
          write(u,'(10(3x,2f12.6))') (tmp(i,j) , j=1,10)
         enddo

         write(u,'("----- W*R ----- ")')
         tmp=cmplx(0d0,0d0)
         call zgemm('n','n',nbd,nbd,nbd,(1d0,0d0),A,nbd,M,nbd,(0d0,0d0),tmp,nbd)
         do i=1,10
          write(u,'(10(3x,2f12.6))') (tmp(i,j) , j=1,10)
         enddo

         write(u,'("----- R*w ----- ")')
         tmp=cmplx(0d0,0d0)
         call zgemm('n','n',nbd,nbd,nbd,(1d0,0d0),M,nbd,wmat,nbd,(0d0,0d0),tmp,nbd)
         do i=1,10
          write(u,'(10(3x,2f12.6))') (tmp(i,j) , j=1,10)
         enddo

         write(u,'("----- w -----")')
         do i=1,10
          write(u,'(10(3x,2f12.6))') (wmat(i,j) , j=1,10)
         enddo

         write(u,'("----- R*w*inv(R) -----")')
         tmp=cmplx(0d0,0d0) ; R=cmplx(0d0,0d0)
         call zgemm('n','n',nbd,nbd,nbd,(1d0,0d0),M,nbd,wmat,nbd,(0d0,0d0),tmp,nbd)
         call zgemm('n','n',nbd,nbd,nbd,(1d0,0d0),tmp,nbd,Minv,nbd,(0d0,0d0),R,nbd)
         do i=1,10
          write(u,'(10(3x,2f12.6))') (R(i,j) , j=1,10)
         enddo


         write(u,'("-----------------------------------------------")')
        end subroutine testdenmat
