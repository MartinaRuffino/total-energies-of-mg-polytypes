  module fmagma
  implicit none

  interface
  subroutine mdsyevd(nrgpu, jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info) bind(c, name='magma_dsyevd_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsyevd

  subroutine mdsyevx(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info) &
                                                                                       bind(c, name='magma_dsyevdx_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsyevx

  subroutine mdsye2x(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info) &
                                                                                 bind(c, name='magma_dsyevdx_2stage_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsye2x

  subroutine mdsygvd(nrgpu, itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info) &
                                                                                       bind(c, name='magma_dsygvd_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsygvd

  subroutine mdsygvx(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, &
                                                                  iwork, liwork, info) bind(c, name='magma_dsygvdx_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsygvx

  subroutine mdsyg2x(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, &
                                                            iwork, liwork, info) bind(c, name='magma_dsygvdx_2stage_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mdsyg2x


  subroutine mzheevd(nrgpu, jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info) &
                                                                                          bind(c, name='magma_zheevd_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzheevd

  subroutine mzheevx(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork, &
                                                                     iwork, liwork, info) bind(c, name='magma_zheevdx_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzheevx








  subroutine mzhee2x(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork, &
                                                               iwork, liwork, info) bind(c, name='magma_zheevdx_2stage_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzhee2x

  subroutine mzhegvd(nrgpu, itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info) &
                                                                                          bind(c, name='magma_zhegvd_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzhegvd

  subroutine mzhegvx(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork, &
                                                                     iwork, liwork, info) bind(c, name='magma_zhegvdx_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzhegvx

  subroutine mzheg2x(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork,&
                                                               iwork, liwork, info) bind(c, name='magma_zhegvdx_2stage_m')
    use iso_c_binding, only : c_char, c_int, c_double
    implicit none
    integer(c_int), intent(in), value :: nrgpu
    integer(c_int), intent(in), value :: itype
    character(kind=c_char), intent(in), value :: jobz
    character(kind=c_char), intent(in), value :: range
    character(kind=c_char), intent(in), value :: uplo
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: a(*)
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(inout) :: b(*)
    integer(c_int), intent(in), value :: ldb
    real(c_double), intent(in), value :: vl
    real(c_double), intent(in), value :: vu
    integer(c_int), intent(in), value :: il
    integer(c_int), intent(in), value :: iu
    integer(c_int), intent(inout) :: m
    real(c_double), intent(inout) :: w(*)
    real(c_double), intent(inout) :: work(*)
    integer(c_int), intent(in), value :: lwork
    real(c_double), intent(inout) :: rwork(*)
    integer(c_int), intent(in), value :: lrwork
    integer(c_int), intent(inout) :: iwork(*)
    integer(c_int), intent(in), value :: liwork
    integer(c_int), intent(out) :: info
  end subroutine mzheg2x

  end interface
  end module fmagma



