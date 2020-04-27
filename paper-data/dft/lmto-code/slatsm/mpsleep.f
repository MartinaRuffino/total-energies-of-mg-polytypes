      subroutine mpsleep(strn,i)
      implicit none
      character strn*(*)
      integer i
      double precision sum,start
      integer j,k,stdo,kloop
      integer :: nloop = 0

      if (nloop == 0) then
        call mpsleep0(nloop)
      endif

      sum = 0; stdo = 6

      call cpudel(0,'',start)
      kloop = nloop/sqrt(dble(i))
      call info2(20,0,0, ' mpsleep  %i sec  '//trim(strn)//' ...',i,2)


      do  j = 1, i*kloop
      do  k = 1, i*kloop
        
        sum = sum+j-k

      enddo
      enddo

      call cpudel(0,'',sum)
      call info2(20,0,0, ' elapsed time = %,3;3d',sum,2)

      end

      subroutine mpsleep0(nloop)
      implicit none
      integer nloop
      integer j,k,stdo
      double precision sum

      sum = 0

      call cpudel(0,'',sum)

      nloop = 1000  ! Adjust to fit time measured by cupdel

      do  j = 1, nloop
      do  k = 1, nloop
        sum = sum+j-k
      enddo
      enddo

      call cpudel(0,'',sum)

      nloop = nloop / sqrt(sum)

      end
