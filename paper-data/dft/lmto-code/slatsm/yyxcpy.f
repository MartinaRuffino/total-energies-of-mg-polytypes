C#define BLAS
      subroutine yyxcpy(n,ar,ai,xr,xi,incx,yr,yi,incy,sw)
C- Complex daxpy, using real arithmetic and complex conjugate of x
      implicit none
C Passed parameters
      logical sw
      integer n,incx,incy
      double precision ar,ai,xr(1),xi(1),yr(1),yi(1)
C Local variables
      integer ix,iy,i

      if (n <= 0) return

C#ifdef BLAS
      call daxpy(n, ar,xr,incx,yr,incy)
      if (sw) then
        call daxpy(n, ai,xi,incx,yr,incy)
        call daxpy(n,-ar,xi,incx,yi,incy)
        call daxpy(n, ai,xr,incx,yi,incy)
      endif
C#elseifC APOLLO | HP
CC --- Do real * real ---
C      if (incx == 1 .and. incy == 1) then
C        if (ar /= 0) then
C          call vec_$dmult_add(yr,xr,n, ar,yr)
C          if (sw) call vec_$dmult_add(yi,xi,n,-ar,yi)
C        endif
C        if (ai /= 0 .and. sw) then
C          call vec_$dmult_add(yi,xr,n, ai,yi)
C          call vec_$dmult_add(yr,xi,n, ai,yr)
C        endif
C      else
C        if (ar /= 0) then
C          call vec_$dmult_add_i(yr,incy,xr,incx,n, ar,yr,incy)
C          if (sw) call vec_$dmult_add_i(yi,incy,xi,incx,n,-ar,yi,incy)
C        endif
C        if (ai /= 0 .and. sw) then
C          call vec_$dmult_add_i(yi,incy,xr,incx,n, ai,yi,incy)
C          call vec_$dmult_add_i(yr,incy,xi,incx,n, ai,yr,incy)
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
C          do  10  i = 1, n
C            yr(iy) = yr(iy) + ar*xr(ix) + ai*xi(ix)
C            yi(iy) = yi(iy) - ar*xi(ix) + ai*xr(ix)
C            ix = ix + incx
C            iy = iy + incy
C   10     continue
C        else
CC     --- Case ar == 0 && ai != 0 ---
C          do  20  i = 1, n
C            yr(iy) = yr(iy) + ai*xi(ix)
C            yi(iy) = yi(iy) + ai*xr(ix)
C            ix = ix + incx
C            iy = iy + incy
C   20     continue
C        endif
C      else
C        if (ar == 0) return
CC     --- Case ar != 0 && ai == 0 ---
C        if (sw) then
C          do  30  i = 1, n
C            yr(iy) = yr(iy) + ar*xr(ix)
C            yi(iy) = yi(iy) - ar*xi(ix)
C            ix = ix + incx
C            iy = iy + incy
C   30     continue
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