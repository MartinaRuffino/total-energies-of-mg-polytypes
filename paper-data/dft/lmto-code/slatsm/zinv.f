      subroutine zinv(ndim, a, lda)
!- Return the inverse of a complex matrix
      implicit none
       integer, intent(in)       :: ndim, lda
       complex(8), intent(inout) :: a(lda,*)
       ! locals
       integer    :: info, lwork
       integer    :: ipiv(ndim)
       complex(8) :: work(ndim*64)
       lwork = ndim*64

       call zgetrf( ndim, ndim, a, lda, ipiv, info )
       if (info /= 0) call rxi('zgetri info=',info)
       call zgetri( ndim, a, lda, ipiv, work, lwork, info )
       if (info /= 0) call rxi('zgetri info=',info)

      end subroutine zinv
