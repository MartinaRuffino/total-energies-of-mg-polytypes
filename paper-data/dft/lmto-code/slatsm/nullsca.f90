! NULL BLACS & SCALAPACK replacement library
! It is incomplete. New interfaces to be added as need arises.
! Dimitar Pashov <d.pashov@gmail.com>

subroutine blacs_gridinit(context,order,np_row,np_col)
   implicit none
   integer context, np_row, np_col
   character*1 order
   np_row = -1; np_col = -1
end subroutine blacs_gridinit

subroutine blacs_gridinfo(context, np_row, np_col, my_row, my_col)
   implicit none
   integer context, np_row, np_col, my_row, my_col
   np_row = 1
   np_col = 1
   my_row = 0
   my_col = 0
end subroutine blacs_gridinfo

subroutine blacs_barrier(context, scope)
   implicit none
   integer context
   character*1 scope
end subroutine blacs_barrier

subroutine blacs_gridexit(context)
   implicit none
   integer context
end subroutine blacs_gridexit

subroutine blacs_exit(contin)
   implicit none
   integer contin
end subroutine blacs_exit


subroutine pdsyevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, &
                                    & iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, il, info, iu, iz, ja, jz, liwork, lwork, m, n, nz
   double precision   abstol, orfac, vl, vu
   integer            desca( * ), descz( * ), iclustr( * ), ifail( * ), iwork( * )
   double precision   a( * ), gap( * ), w( * ), work( * ), z( * )
end subroutine pdsyevx


subroutine pdsyevr( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work,&
                                                                              & lwork, iwork, liwork, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, il, info, iu, iz, ja, jz, liwork, lwork, m, n, nz
   double precision vl, vu
   integer            desca( * ), descz( * ), iwork( * )
   double precision   a( * ), w( * ), work( * ), z( * )
end subroutine pdsyevr


subroutine pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
   implicit none
   character          jobz, uplo
   integer            ia, info, iz, ja, jz, liwork, lwork, n
   integer            desca( * ), descz( * ), iwork( * )
   double precision   a( * ), w( * ), work( * ), z( * )
end subroutine pdsyevd


subroutine pdsygvx( ibtype, jobz, range, uplo, n, a, ia, ja, desca, b, ib, jb, descb, vl, vu, il, iu, abstol,&
                  & m, nz, w, orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, ib, ibtype, il, info, iu, iz, ja, jb, jz, liwork, lwork, m, n, nz
   double precision   abstol, orfac, vl, vu
   integer            desca( * ), descb( * ), descz( * ), iclustr( * ), ifail( * ), iwork( * )
   double precision   a( * ), b( * ), gap( * ), w( * ), work( * ), z( * )
end subroutine pdsygvx


subroutine pzheevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, &
                           & jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, il, info, iu, iz, ja, jz, liwork, lrwork, lwork, m, n, nz
   double precision   abstol, orfac, vl, vu
   integer            desca( * ), descz( * ), iclustr( * ), ifail( * ), iwork( * )
   double precision   gap( * ), rwork( * ), w( * )
   complex*16         a( * ), work( * ), z( * )
end subroutine pzheevx


subroutine pzheevr( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, m, nz, w, z, iz, jz, descz, work,&
                                                                  & lwork, rwork, lrwork, iwork, liwork, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, il, info, iu, iz, ja, jz, liwork, lrwork, lwork, m, n, nz
   double precision vl, vu
   integer            desca( * ), descz( * ), iwork( * )
   double precision   w( * ), rwork( * )
   complex*16         a( * ), work( * ), z( * )
end subroutine pzheevr


subroutine pzheevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork,&
                                                                                             & liwork, info )
   implicit none
   character          jobz, uplo
   integer            ia, info, iz, ja, jz, liwork, lrwork, lwork, n
   integer            desca( * ), descz( * ), iwork( * )
   double precision   rwork( * ), w( * )
   complex*16         a( * ), work( * ), z( * )
end subroutine pzheevd


subroutine pzhegvx( ibtype, jobz, range, uplo, n, a, ia, ja, desca, b, ib, jb, descb, vl, vu, il, iu, abstol,&
   & m, nz, w, orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info )
   implicit none
   character          jobz, range, uplo
   integer            ia, ib, ibtype, il, info, iu, iz, ja, jb, jz, liwork, lrwork, lwork, m, n, nz
   double precision   abstol, orfac, vl, vu
   integer            desca( * ), descb( * ), descz( * ), iclustr( * ), ifail( * ), iwork( * )
   double precision   gap( * ), rwork( * ), w( * )
   complex*16         a( * ), b( * ), work( * ), z( * )
end subroutine pzhegvx

! function pdlamch( ictxt, cmach )
!    implicit none
!    character          cmach
!    integer            ictxt
!    real(8) :: pdlamch
!    procedure(real(8)) dlamch
!    pdlamch = dlamch(cmach)
! end function pdlamch

