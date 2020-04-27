         subroutine diagdenmat(nbd,denmat,evr,evrinv,w)
          implicit none
          integer,    intent(in)    :: nbd
          complex*16, intent(inout) :: denmat(nbd,nbd)
          complex*16, intent(out)   :: evr(nbd,nbd),evrinv(nbd,nbd)
          complex*16, intent(out)   :: w(nbd)
          complex*16 :: evl(nbd,nbd)
          complex*16 :: cworkvec(8*nbd)
          real*8     :: rworkvec(8*nbd)
          complex*16 :: tmp1(nbd,nbd)
          integer    :: ierr,ipiv(nbd),i,j

C     ... Hermitianize denmat
!         tmp1=denmat
!         do i=1,nbd
!          do j=1,nbd
!           denmat(i,j) = 0.5*(tmp1(i,j)+conjg(tmp1(j,i)) )
!          enddo
!         enddo

C     ... keep only the diagonal of denmat
C     ... This is just a test: when denmat is diagonal R is the identity.
!         tmp1=cmplx(0d0,0d0)
!         do i=1,nbd
!          tmp1(i,i)=denmat(i,i)
!         enddo
!         denmat=tmp1


          tmp1=denmat
C     ... Get eigenvalues and (right) eigenvectors
          call zgeev("n","v",nbd,tmp1,nbd,w,evl,nbd,evr,nbd,
     .               cworkvec,8*nbd,rworkvec,ierr)


C     ... Invert eigenvector matrix
          evrinv=evr
          call zgetrf(nbd,nbd,evrinv,nbd,ipiv,ierr)
          call zgetri(nbd,evrinv,nbd,ipiv,cworkvec,8*nbd,ierr)

         end subroutine diagdenmat
