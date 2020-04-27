
   subroutine eigapprox(n, h, l)

!    use mod_io

   implicit none

   integer, parameter :: dp = 8

   integer, intent(in) :: n
   real(dp), intent(inout) :: h(n,n), l(n)

   integer, save :: cnt = 0
   integer :: i, j, u
   real(dp) :: dhc(n), dc(n), c(n,n), dh(n,n), cdhc, alpha

   real(dp), save, allocatable :: h0(:,:), l0(:), c0(:,:)

   character(len=100) :: vector_fmt

   procedure(real(dp)) :: ddot



   alpha = 0.05_dp



   cnt = cnt + 1


!    u = 6453
!    write(vector_fmt, "('(i5,',i0,'(x,f12.6))')") n

!    write(u, '("approx in iter: ", i0)') cnt
!    write(u, '(a,/,6(6(x,f12.6),/))') 'h', h
!    write(u, '(a,/,1(6(x,f12.6),/))') 'l', l
! !    write(u, '(a,/,6(6(x,f12.6),/))') 'c', c

   dhc = 0.0_dp

!    if (cnt < 2 .or. cnt > 51 ) then !or. cnt > 14.or. cnt > 12 .or. mod(cnt,5) == 4 .and. cnt < 30 .or. (cnt > 30 )
   if (cnt < 2000) then
      print *, 'diag', cnt
      if (.not. allocated(h0)) allocate(h0(n,n))
      if (.not. allocated(l0)) allocate(l0(n))
      if (.not. allocated(c0)) allocate(c0(n,n))

      h0 = h
      call dyn_dsyevd(n, h, l)
      l0 = l
      c0 = h
   else if (cnt < 52) then
      print *, 'power', cnt
      c = c0

      call dsymm ('l', 'l', n, n, -alpha, h, n, c0, n, 1.0_dp, c, n)

      do i = 1, n
         call daxpy(n, alpha*l0(i), c0(:,i), 1, c(:,i), 1)
      end do

      do i = 1, n
         call dsymv ('u', n, 1.0_dp, h, n, c0(:,i), 1, 0.0_dp, dhc, 1)
         l(i) = ddot(n, c0(:,i), 1, dhc, 1)
      end do

      call mgschm(n, c)

      do i = 1, n
         c(:,i) = c(:,i)/sqrt(ddot(n, c(:,i), 1, c(:,i), 1))
      end do



      l0 = l
      h0 = h
      h = c
      c0 = c

   else
      print *, 'pert', cnt
      dh = h - h0
      c = c0

      do i = 1, n
         call dsymv ('u', n, 1.0_dp, dh, n, c0(:,i), 1, 0.0_dp, dhc, 1)
         do j = 1, n
            if (i == j) cycle
            cdhc = ddot(n, c0(:,j), 1, dhc, 1)/(l0(i) - l0(j))
            call daxpy(n, cdhc, c0(:,j), 1, c(:,i), 1)
         end do
      end do

      call mgschm(n, c)

      do i = 1, n
         c(:,i) = c(:,i)/sqrt(ddot(n, c(:,i), 1, c(:,i), 1))
      end do

      do i = 1, n
         call dsymv ('u', n, 1.0_dp, h, n, c0(:,i), 1, 0.0_dp, dhc, 1)
         l(i) = ddot(n, c0(:,i), 1, dhc, 1)
      end do

      l0 = l
      h0 = h
      h = c
      c0 = c

   end if


!    write(u, '("approx out iter: ", i0)') cnt
!    write(u, '(a,/,6(6(x,f12.6),/))') 'h', h
!    write(u, '(a,/,1(6(x,f12.6),/))') 'l', l

!    write(u, trim(vector_fmt)) cnt, l
!
!    u = 7000+cnt
! !    open(newunit=u, file='ovlp.pgm', action='write')
!    write(u, '(a,/,a,x,i4,/,i4,x,i4,/,i6)') 'P2', "# |C| it: ", cnt, n, n, 65535
!    do i = 1, n
!       do j = 1, n
!          write(u, '(x,i9)' ,advance='no') int(32768*sqrt(abs(ddot(n, c0(:,i), 1, c0(:,j), 1))))
!       end do
!       write(u,'(" ")')
!    end do
!    close(u)

   end subroutine eigapprox





   subroutine mgschm(n, v)

   implicit none

   integer, parameter :: dp = 8

   integer, intent(in) :: n
   real(dp), intent(inout) :: v(n,n)


   real(dp) :: q(n)
   integer :: i, j

   procedure(real(dp)) :: ddot


   do i = 1, n
      q = v(:,i)*(1.0_dp/ddot(n, v(:,i), 1, v(:,i), 1))
      do j = i + 1, n
         v(:,j) = v(:,j) - ddot(n, q, 1, v(:,j), 1)*q
      end do
   end do

   end subroutine mgschm
















   subroutine dyn_dsyevd(n, a, l)


   implicit none

   integer, parameter :: dp = 8

   integer, intent(in) :: n
   real(dp), intent(inout) :: a(n,n)
   real(dp), intent(out) :: l(n)

   integer :: info, lwork, liwork

   integer, allocatable :: iwork(:)
   real(dp), allocatable :: work(:)



   allocate(work(1), iwork(1))

   lwork = -1
   liwork = -1
   call dsyevd('v', 'u', n, a, n, l, work, lwork, iwork, liwork, info)

   lwork = nint(work(1))
   liwork = iwork(1)

   deallocate(iwork, work)
   allocate(work(lwork), iwork(liwork))

   call dsyevd('v', 'u', n, a, n, l, work, lwork, iwork, liwork, info)

   deallocate(iwork, work)

   if (info /= 0) stop 'there is an issue with dsyevd'

   end subroutine dyn_dsyevd








!    subroutine carpar(n,h,l)
!
!
!    implicit none
!
!    integer, parameter :: dp = 8
!
!    integer, intent(in) :: n
!    real(dp), intent(inout) :: h(n,n), l(n)
!
!
!    integer, save :: cnt = 0
!    integer :: i
!    real(dp) :: c(n,n), alpha, dhc(n)
! !    real(dp), save, allocatable :: h0(:,:), l0(:)
!    real(dp), save, allocatable :: c0(:,:)
!
!
!    procedure(real(dp)) :: ddot
!
!
!
!    cnt = cnt + 1
!
!
!    alpha = 0.01_dp
!
!
!    if (cnt < 2  ) then !or. cnt > 14.or. cnt > 12 .or. mod(cnt,5) == 4
!       print *, 'diag'
! !       if (.not. allocated(h0)) allocate(h0(n,n))
! !       if (.not. allocated(l0)) allocate(l0(n))
!       if (.not. allocated(c0)) allocate(c0(n,n))
!
! !       h0 = h
!       call dyn_dsyevd(n, h, l)
! !       l0 = l
!       c0 = h
!    else
! !       dh = h - h0
!
!       c = c0
!
!       call dgemm ('t', 'n', n, n, n, -alpha, h, n, c0, n, 1.0_dp, c, n)
!
!       call mgschm(n, c)
!
!       do i = 1, n
!          c(:,i) = c(:,i)/sqrt(ddot(n, c(:,i), 1, c(:,i), 1))
!       end do
!
!       do i = 1, n
!          call dsymv ('u', n, 1.0_dp, h, n, c0(:,i), 1, 0.0_dp, dhc, 1)
!          l(i) = ddot(n, c0(:,i), 1, dhc, 1)
!       end do
!
!
! !       h0 = h
!       h = c
!       c0 = c
! !       l0 = l
!
!    end if
!
!
!
!
!
!    end subroutine carpar
