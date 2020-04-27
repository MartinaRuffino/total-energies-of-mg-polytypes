C#define BLAS
      subroutine qyyaxpy(n,ar,ai,xr,xi,incx,yr,yi,incy,sw)
C- Complex qaxpy, using real arithmetic
      implicit none
C Passed parameters
      logical sw
      integer n,incx,incy
      real*16 ar,ai,xr(1),xi(1),yr(1),yi(1)
C Local variables
      integer ix,iy,i

      if (n <= 0) return

C#ifdef BLAS
      call qaxpy(n, ar,xr,incx,yr,incy)
      if (sw) then
        call qaxpy(n,-ai,xi,incx,yr,incy)
        call qaxpy(n, ar,xi,incx,yi,incy)
        call qaxpy(n, ai,xr,incx,yi,incy)
      endif
C#elseifC APOLLO | HP
CC --- Do real * real ---
C      if (incx == 1 .and. incy == 1) then
C        if (ar /= 0) then
C          call vec_$dmult_add(yr,xr,n, ar,yr)
C          if (sw) call vec_$dmult_add(yi,xi,n, ar,yi)
C        endif
C        if (ai /= 0 .and. sw) then
C          call vec_$dmult_add(yi,xr,n, ai,yi)
C          call vec_$dmult_add(yr,xi,n,-ai,yr)
C        endif
C      else
C        if (ar /= 0) then
C          call vec_$dmult_add_i(yr,incy,xr,incx,n, ar,yr,incy)
C          if (sw) call vec_$dmult_add_i(yi,incy,xi,incx,n, ar,yi,incy)
C        endif
C        if (ai /= 0 .and. sw) then
C          call vec_$dmult_add_i(yi,incy,xr,incx,n, ai,yi,incy)
C          call vec_$dmult_add_i(yr,incy,xi,incx,n,-ai,yr,incy)
C        endif
C      endif
C#elseC
C      ix = 1
C      iy = 1
C      if (incx < 0) ix = (1-n)*incx + 1
C      if (incy < 0) iy = (1-n)*incy + 1
C      if (sw .and. ai /= 0) then
C        if (ar /= 0) then
CC     --- Case ar != 0 && ai != 0 ---
C          if (incx == 1 .and. incy == 1 .and. .false.) then
C            do  10  ix = 1, n
C              yr(ix) = yr(ix) + ar*xr(ix) - ai*xi(ix)
C              yi(ix) = yi(ix) + ar*xi(ix) + ai*xr(ix)
C   10       continue
C            return
C          else
C            do  11  i = 1, n
C              yr(iy) = yr(iy) + ar*xr(ix) - ai*xi(ix)
C              yi(iy) = yi(iy) + ar*xi(ix) + ai*xr(ix)
C              ix = ix + incx
C              iy = iy + incy
C   11       continue
C          endif
C        else
CC     --- Case ar == 0 && ai != 0 ---
C          if (incx == 1 .and. incy == 1 .and. .false.) then
C            do  20  ix = 1, n
C              yr(ix) = yr(ix) - ai*xi(ix)
C              yi(ix) = yi(ix) + ai*xr(ix)
C   20       continue
C            return
C          else
C            do  21  i = 1, n
C              yr(iy) = yr(iy) - ai*xi(ix)
C              yi(iy) = yi(iy) + ai*xr(ix)
C              ix = ix + incx
C              iy = iy + incy
C   21       continue
C          endif
C        endif
C      else
C        if (ar == 0) return
CC     --- Case ar != 0 && ai == 0 ---
C        if (sw) then
C          if (incx == 1 .and. incy == 1) then
C            do  30  ix = 1, n
C              yr(ix) = yr(ix) + ar*xr(ix)
C              yi(ix) = yi(ix) + ar*xi(ix)
C   30       continue
C            return
C          else
C            do  31  i = 1, n
C              yr(iy) = yr(iy) + ar*xr(ix)
C              yi(iy) = yi(iy) + ar*xi(ix)
C              ix = ix + incx
C              iy = iy + incy
C   31       continue
C          endif
C        else
C          do  40  i = 1, n
C            yr(iy) = yr(iy) + ar*xr(ix)
C            ix = ix + incx
C            iy = iy + incy
C   40     continue
C        endif
C      endif
C#endif

      end
