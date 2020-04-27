      subroutine bsortl(vec,nkd)
C- Sorts a table of vectors by increasing length
      implicit none
      integer nkd
      double precision vec(3,1)
      double precision dsqr,xx
      integer i,is,j
      logical lswap

      dsqr(i) = vec(1,i)**2 + vec(2,i)**2 + vec(3,i)**2

 10   continue
      lswap = .false.
      do  is = 2, nkd
        if (dsqr(is) < dsqr(is-1)) then
        lswap = .true.
          do  j = 1, 3
          xx = vec(j,is)
          vec(j,is) = vec(j,is-1)
          vec(j,is-1) = xx
          enddo
        endif
      enddo
      if (lswap) goto 10

c      print *, nkd, ' lattice vectors'
c      print 300, (dsqr(is), vec(1,is), vec(2,is), vec(3,is), is=1,nkd)
c  300 format(4f12.5)
c      stop
      end
