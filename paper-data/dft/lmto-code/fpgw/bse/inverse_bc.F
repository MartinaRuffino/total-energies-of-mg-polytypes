      MODULE inverse_bc

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: inversemat_comp
      PUBLIC :: inversemat_real

      CONTAINS

      FUNCTION inversemat_comp( mat ) RESULT(invmat_comp)
C! works out the inverse of a matrix

      COMPLEX*16, DIMENSION(:,:), INTENT(in)           :: mat
      COMPLEX*16, DIMENSION(size(mat,1),size(mat,2)) :: invmat_comp

      COMPLEX*16, DIMENSION(size(mat,1))             :: work!work array for LAPACK
      INTEGER, DIMENSION(size(mat,1))                    :: ipiv   ! pivot indices
      INTEGER                                            :: n, info

C      ! External procedures defined in LAPACK
      EXTERNAL ZGETRF
      EXTERNAL ZGETRI

C      ! Store A in Ainv to prevent it from being overwritten by LAPACK

      invmat_comp = mat
      n = SIZE(mat,1)

C     ! ZGETRF computes an LU factorization of a general M-by-N matrix A
C     ! using partial pivoting with row interchanges.

      CALL ZGETRF(n, n, invmat_comp, n, ipiv, info)

      if (info /= 0) then
        stop 'Matrix is numerically singular!'
      end if

C     ! ZGETRI computes the inverse of a matrix using the LU factorization
C     ! computed by DGETRF.

      CALL ZGETRI(n, invmat_comp, n, ipiv, work, n, info)

      if (info /= 0) then
        stop 'Matrix inversion failed!'
      end if

      END FUNCTION inversemat_comp

      FUNCTION inversemat_real( mat ) RESULT(invmat_real)

      REAL(8), DIMENSION(:,:), INTENT(in)           :: mat
      REAL(8), DIMENSION(size(mat,1),size(mat,2)) :: invmat_real

      REAL(8), DIMENSION(size(mat,1))             :: work!work array for LAPACK
      INTEGER, DIMENSION(size(mat,1))                    :: ipiv   ! pivot indices
      INTEGER                                            :: n, info

C      ! External procedures defined in LAPACK
      EXTERNAL DGETRF
      EXTERNAL DGETRI

C      ! Store A in Ainv to prevent it from being overwritten by LAPACK

      invmat_real = mat
      n = SIZE(mat,1)

C     ! DGETRF computes an LU factorization of a general M-by-N matrix A
C     ! using partial pivoting with row interchanges.

      CALL DGETRF(n, n, invmat_real, n, ipiv, info)

      if (info /= 0) then
        stop 'Matrix is numerically singular!'
      end if

C     ! DGETRI computes the inverse of a matrix using the LU factorization
C     ! computed by DGETRF.

      CALL DGETRI(n, invmat_real, n, ipiv, work,n,info)!last n maybe bigger,50*n?

      if (info /= 0) then
        stop 'Matrix inversion failed!'
      end if

      END FUNCTION inversemat_real

      END
