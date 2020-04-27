      subroutine sleep1(j)

      print *, 'starting sleep'
      do  i = 1, 3000000
        call sleep2(i,j,l)
      enddo
      print *, 'end sleep',l
      end

      subroutine sleep2(i,j,l)

      l = 0
      do  k = 1, 100*j
        l = l + j + (-1)**k
      enddo
      end
