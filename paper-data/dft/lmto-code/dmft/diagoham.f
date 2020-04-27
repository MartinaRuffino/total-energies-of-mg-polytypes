         subroutine diagoham(matdim,ham,zek)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  IT diagonalizes the matrix HAM     !
         !  and keeps only right eigenvalues.  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          use mod_utils, only : mergesort
          implicit none
          integer,    intent(in)  :: matdim             ! matrix dimension
          complex(8), intent(in)  :: ham(matdim,matdim) ! matrix
          complex(8), intent(out) :: zek(matdim)        ! eigenvalues
          complex(8) :: evl(matdim,matdim), evr(matdim,matdim)              ! eigenvectors (not referenced)
          complex(8) :: cworkvec(8*matdim)
          real(8)    :: rworkvec(8*matdim)
          integer    :: ierr

          real(8) :: zek_r(matdim)
          integer :: sk(matdim)

          ! perform diagonalization using Lapack routines
          call zgeev("N","N",matdim,ham,matdim,zek,evl,matdim,evr,
     .               matdim,cworkvec,8*matdim,rworkvec,ierr)


          zek_r = real(zek)
          call mergesort(zek_r, sk, matdim)
          call reorder(sk, matdim, zek, 2, 2)

          end subroutine diagoham


         subroutine diagoham_full(matdim,ham,zek,evl,evr)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  IT diagonalizes the matrix HAM                    !
         !  keeps Right and Left eigenvectors and eigenvalues !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          use mod_utils, only : mergesort
          implicit none
          integer,    intent(in)  :: matdim             ! matrix dimension
          complex(8), intent(in)  :: ham(matdim,matdim) ! matrix
          complex(8), intent(out) :: zek(matdim)        ! eigenvalues
          complex(8), intent(out) :: evl(matdim,matdim), evr(matdim,matdim) ! left and right column eigenvectors
          complex(8) :: cworkvec(8*matdim)
          real(8)    :: rworkvec(8*matdim)
          integer    :: ierr

          real(8) :: zek_r(matdim)
          integer :: sk(matdim)

          ! perform diagonalization using Lapack routines
          call zgeev("V","V",matdim,ham,matdim,zek,evl,matdim,evr,
     .               matdim,cworkvec,8*matdim,rworkvec,ierr)

          zek_r = real(zek)
          call mergesort(zek_r, sk, matdim)
          call reorder(sk, matdim, zek, 2, 2)
          call reorder(sk, matdim, evl, 2*matdim, 2*matdim)
          call reorder(sk, matdim, evr, 2*matdim, 2*matdim)

         end subroutine diagoham_full

