
  subroutine magma_dsyevd_m(nrgpu, jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsyevd_m

  subroutine magma_dsyevdx_m(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsyevdx_m

  subroutine magma_dsyevdx_2stage_m(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, iwork, liwork, info)&
                                                                                                                          bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsyevdx_2stage_m

  subroutine magma_dsygvd_m(nrgpu, itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsygvd_m

  subroutine magma_dsygvdx_m(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, iwork, &
                                                                                                       liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsygvdx_m

  subroutine magma_dsygvdx_2stage_m(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, iwork,&
                                                                                                             liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_dsygvdx_2stage_m


  subroutine magma_zheevd_m(nrgpu, jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zheevd_m

  subroutine magma_zheevdx_m(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork, iwork, &
                                                                                                         liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zheevdx_m

  subroutine magma_zheevdx_2stage_m(nrgpu, jobz, range, uplo, n, a, lda, vl, vu, il, iu, m, w, work, lwork, rwork, lrwork, &
                                                                                                iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zheevdx_2stage_m

  subroutine magma_zhegvd_m(nrgpu, itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info) &
                                                                                                                       bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zhegvd_m

  subroutine magma_zhegvdx_m(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, rwork, &
                                                                                        lrwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zhegvdx_m

  subroutine magma_zhegvdx_2stage_m(nrgpu, itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, m, w, work, lwork, rwork,&
                                                                                              lrwork, iwork, liwork, info) bind(c)
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
    info = -1
    call rx("MAGMA is not enabled. Thou shalt not be here!")
  end subroutine magma_zhegvdx_2stage_m
