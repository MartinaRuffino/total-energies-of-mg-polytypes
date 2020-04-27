       subroutine setnth(nth1,nth3,is,nclasp,nangl)
       integer i,is,nth3,nclasp,nangl,nth1(nclasp+nangl)
       nth1(is) = nth3
       do i = nclasp+1,nclasp+nangl
         nth1(i) = 0
C         write(*,*)'AA-nth1(',i,')',nth1(i)
       enddo
       end
