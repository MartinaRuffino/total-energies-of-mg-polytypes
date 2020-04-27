      subroutine dshell(n,array)
      implicit none
      integer n
      double precision array(n)
      integer i,j,k,inc
      double precision v

C ... Get the largest increment
      if (n <= 1) return
      inc = 1
   10 continue
      inc = 3*inc+1
      if (inc < n) goto 10

C ... Loop over partial sorts
   12 continue
        inc = inc/3
C   ... Outer loop of straight insertion
        do  11  i = inc+1, n
          v = array(i)
          j = i
C     ... Inner loop of straight insertion
   20     continue
          if (array(j-inc) > v) then
            array(j) = array(j-inc)
            j = j-inc
            if (j <= inc) goto 21
            goto 20
          endif
   21     continue
          array(j) = v
   11   continue
      if (inc > 1) goto 12
      end

