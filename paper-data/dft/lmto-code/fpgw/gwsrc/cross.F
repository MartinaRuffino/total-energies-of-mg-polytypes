      subroutine fermids(k1,k2,k3,ds)
      real(8):: k1(3),k2(3),k3(3),ds(3,3),ee(3)
      call cross_x(k1-k3,k2-k3,ee)

      do i=1,3
         do j=1,3
            ds(i,j)=ds(i,j)+ee(i)*ee(j)/sqrt(sum(ee**2))
         enddo
      enddo

      end

c --- external product
      subroutine cross_x(a,b,c)
      implicit none
      real*8 a(3), b(3), c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

!> This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) 
!! the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      subroutine minv33tp(plat,qlat)
      implicit none
      real(8),intent(in)::  plat(3,3)
      real(8),intent(out):: qlat(3,3)
      real(8):: det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

!>- Inverts 3X3 matrix
      subroutine minv33(matrix,inverse)
Co Outputs
Co   inverse, as modified according to iopt
Co   det:      determinant
      implicit none
c      integer:: iopt=0
      real(8), intent(in) :: matrix(3,3)
      real(8), intent(out) :: inverse(3,3)
      real(8) :: det,ddot
      call cross_x(matrix(1,2),matrix(1,3),inverse     )
      call cross_x(matrix(1,3),matrix     ,inverse(1,2))
      call cross_x(matrix     ,matrix(1,2),inverse(1,3))
      det = ddot(3,matrix,1,inverse,1)
      if (abs(det) ==0d0) call rx( 'minv33: vanishing determinant')
      inverse = transpose(inverse)
      inverse = inverse/det
      end

C      subroutine dinv33x(plat,qlat)
CC- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
CCr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
C      implicit none
C      double precision plat(3,3),qlat(3,3),det
C      call cross_x(plat(1,2),plat(1,3), qlat(1,1))
C      call cross_x(plat(1,3),plat(1,1), qlat(1,2))
C      call cross_x(plat(1,1),plat(1,2), qlat(1,3))
C      det  = sum( plat(1:3,1)*qlat(1:3,1) )
C      qlat = qlat/det
C      end

C      subroutine dinv33y(plat,qlat,det)
CC- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
CCr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
C      implicit none
C      double precision plat(3,3),qlat(3,3),det
C      call cross_x(plat(1,2),plat(1,3), qlat     )
C      call cross_x(plat(1,3),plat     , qlat(1,2))
C      call cross_x(plat     ,plat(1,2), qlat(1,3))
C      det  = sum( plat(1:3,1)*qlat(1:3,1) )
C      qlat = qlat/det
C      end

Cc-Taken from Ferdi's GW  -----------------------------------------------------
C      subroutine dinv33(matrix,iopt,inverse,det)
CC- Inverts 3X3 matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   inverse: input matrix
CCi   iopt:  if 0, usual inverse
CCi             1, transpose of inverse
CCo Outputs
CCo   inverse, as modified according to iopt
CCo   det:      determinant
CC ----------------------------------------------------------------
C      implicit none
C      integer iopt,i,j
C      double precision matrix(3,3),inverse(3,3),det,ddot
CCstop2rx 2013.08.09 kino      if(iopt<0.or.iopt>1) stop 'dinv33:wrong iopt'
C      if(iopt<0.or.iopt>1) call rx( 'dinv33:wrong iopt')
C      call cross_x(matrix(1,2),matrix(1,3),inverse(1,1))
C      call cross_x(matrix(1,3),matrix(1,1),inverse(1,2))
C      call cross_x(matrix(1,1),matrix(1,2),inverse(1,3))
C      det = ddot(3,matrix(1,1),1,inverse(1,1),1)
C
CCstop2rx 2013.08.09 kino      if (abs(det) ==0d0) stop 'dinv33: vanishing determinant'
C      if (abs(det) ==0d0) call rx( 'dinv33: vanishing determinant')
C      if (iopt == 0) inverse = transpose(inverse)
C      inverse = inverse/det
Cc      double precision xx
Cc      if (iopt >= 2) det = det/(8*datan(1d0))
Cc      if (mod(iopt,2) == 0) then
Cc        do  10  i = 1, 3
Cc          do  10  j = i+1, 3
Cc          xx = inverse(i,j)
Cc          inverse(i,j) = inverse(j,i)
Cc          inverse(j,i) = xx
Cc   10   continue
Cc      endif
Cc     call dscal(9,1/det,inverse,1)
Cc      return
C      end
c --- cross
      subroutine cross(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

c --- tripl
      doubleprecision function tripl(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
c     dimension a(3),b(3),c(3)
c ... g77 needs this rewriting for opt > O1
      dimension a(0:2),b(0:2),c(0:2)
c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c    .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c ... g77 needs this rewriting
      tmp = 0.d0
      do i = 0,2
        j = mod(i + 1, 3)
        k = mod(i + 2, 3)
        tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
      enddo
      tripl = tmp
      end

Cc --- sprod
C      doubleprecision function sprod(a,b)
C      implicit real*8 (a-h,p-z), integer(o)
C      dimension a(3),b(3)
C      sprod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
C      return
C      end
