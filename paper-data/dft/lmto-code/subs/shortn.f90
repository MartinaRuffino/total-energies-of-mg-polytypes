      subroutine shortn(p,p1,dlat,nkd)
!C- Get p1 = shortest vector such that p1-p is a lattice vector.
!C ----------------------------------------------------------------
!Ci Inputs
!Ci   p     :vector to shorten
!Ci   dlat  :lattice vectors, sorted by increasing length
!Ci   nkd   :number of dlat
!Co Outputs
!Co   p1    :shortened vector
!Cr Remarks
!Cr   A slightly skewed norm is used to make result unique.
!Cr   The first vector in the list must be the zero vector.
!c ----------------------------------------------------------------
      implicit none
!C ... Passed parameters
      integer nkd
      real(8), intent(in) :: p(3),dlat(3,nkd)
      real(8), intent(out) :: p1(3)
!C ... Local parameters
      integer irep,k0,k
      double precision p2,critk0,dd,crit,v(3)
!       anrm2(x,y,z) = x*x*1.00001d0 + y*y*1.00002d0 + z*z*1.00003d0
!      .                -x*0.000004d0 - y*0.000003d0 - z*0.000002d0

      p1(1:3) = p(1:3)
      do  irep = 1, 20
        p2 = anrm2(p1)
        k0 = 1
        critk0 = p2
        do  k = 1, nkd
          dd = sum(dlat(1:3,k)*dlat(1:3,k))
!C         Only valid if dlat is sorted:
          if (dd > p2*4d0) exit
          v = p1 + dlat(1:3,k)
          crit = anrm2(v)
          if (crit < critk0) then
            k0 = k
            critk0 = crit
          endif
        enddo
        if (k0 == 1) return
        p1(1:3) = p1(1:3) + dlat(1:3,k0)
      enddo
      call rx('shortn: shortest vector not found')
      contains
         function anrm2(v)
            implicit none
            real(8), intent(in) :: v(3)
            real(8) :: anrm2

!             anrm2 = sum(v*v)
            anrm2 = v(1)*v(1)*1.00001d0 + v(2)*v(2)*1.00002d0          &
                  + v(3)*v(3)*1.00003d0                                &
                  - v(1)*0.000004d0 - v(2)*0.000003d0 - v(3)*0.000002d0

         end function anrm2
      end subroutine shortn


      subroutine simpleshortn(p,p1,dlat,nkd)
         implicit none
         integer, intent(in) :: nkd
         real(8), intent(in) :: p(3),dlat(3,nkd)
         real(8), intent(out) :: p1(3)

         real(8) :: v(3), w(3), v2, w2
         integer :: k

         w = p
         w2 = sum(p*p)
         do k = 1, nkd
            v = p - dlat(1:3,k)
            v2 = sum(v*v)
            if (v2 < w2) then
               w = v
               w2 = v2
!                exit
            end if
         end do

         p1 = w

      end subroutine simpleshortn

      subroutine directshortn(p,plat,ilat)
!       Shorten vector 'p' to cell 'plat' using precomputed ilat = plat^-1
!       'plat' vectors shall be stored in the columns.
         implicit none
         integer, parameter :: dp = 8
         real(dp), intent(in) :: plat(3,3), ilat(3,3)
         real(dp), intent(inout) :: p(3)
         interface
            pure elemental subroutine p_e_s_r8(a)
               real(8), intent(inout) :: a
            end subroutine p_e_s_r8
         end interface
         real(dp) :: r(3)

         procedure(p_e_s_r8) :: sinclmod1

         r = matmul(ilat, p)
         call sinclmod1(r)
         p = matmul(plat, r)

      end subroutine directshortn

      pure elemental subroutine sinclmod1(a)
!     Shorten the N-dimensional vector 'a' to a vector within an N-dimensional cell
!     described by the identity matrix I. Preserve direction on boundary cases (inclusive boundary).
          implicit none
          integer, parameter :: dp = 8
          real(dp), intent(inout) :: a
          a = a + sign(1.0_dp, a)*floor(0.5_dp - abs(a))
      end subroutine sinclmod1

      pure elemental function finclmod1(a) result (r)
!     Shorten the N-dimensional vector 'a' to a vector within an N-dimensional cell
!     described by the identity matrix I. Preserve direction on boundary cases (inclusive boundary).
          implicit none
          integer, parameter :: dp = 8
          real(dp), intent(in) :: a
          real(dp) :: r
          r = a + sign(1.0_dp, a)*floor(0.5_dp - abs(a))
      end function finclmod1

      pure elemental function finclmod(a,b) result (r)
!     Shorten the N-dimensional vector 'a' to a vector within an N-dimensional cell
!     described by the matrix I*b. Preserve direction on boundary cases (inclusive boundary).
!     The I cell routines '[s|f]inclmod1' are likely more useful.
         implicit none
         integer, parameter :: dp = 8
         real(dp), intent(in) :: a, b
         real(dp) :: r, ab
         ab = a/b
         r = a + sign(1.0_dp, ab)*floor(0.5_dp - abs(ab))*b
      end function finclmod

