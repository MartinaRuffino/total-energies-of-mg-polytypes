      integer(4) function icompvv2(iv1,iv2,n,comment)
c takao modified from icompvv by ferdi
c compare two integer vectors iv1,iv2
      implicit none
      integer(4):: i,n, iv1(n),iv2(n)
      character*(*) comment
      icompvv2 = 0  ! (ML) bug fix
      do i = 1,n
        if(iv1(i) /= iv2(i)) then
          icompvv2 = iv1(i) - iv2(i)
          write (*,'(a)') comment
          call rx( 'icompvv2 err: ')
        endif
      end do
      end
