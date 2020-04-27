      subroutine vprds(nr,n,wt,a,b,ab)
C- Make a set of vector products with spec'd weighting function
      implicit none
      integer nr,nlml,n,i,ir
      double precision wt(nr),a(nr,n),b(nr,n),ab(0:n),sum

      ab(0) = 0
      do  10  i = 1, n
        sum = 0d0
        do  12  ir = 1, nr
   12   sum = sum + a(ir,i)*b(ir,i)*wt(ir)
        ab(i) = sum
        ab(0) = ab(0) + sum
   10 continue
      end
