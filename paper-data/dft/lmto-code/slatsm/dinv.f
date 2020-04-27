      integer function dinv(cs,n,lda,a)
C- Inversion of a double precision matrix
      implicit none
      character*1 cs
      integer n,lda
      double precision a(lda,lda)
      integer ldw,i
      real(8), allocatable :: w(:)

      ldw = n
      allocate(w(ldw*(n+1)))
      call dqinv(cs,a,lda,2,n,w,ldw,i)
      deallocate(w)
      dinv = i
      end
