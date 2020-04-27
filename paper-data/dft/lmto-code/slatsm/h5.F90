#define H5FOO "h5foo © 2018 Dimitar Pashov"
#define H5FOO_VERSION 0.1

! There is already https://github.com/scivision/oo_hdf5_fortran, which is pretty decent.
! The following is a more explicit take on the idea and an hdf5 learning experience.

#define assertm(E,M) if(.not.(E))then;write(stderr,'(a,":",i0,1x,a)')__FILE__,__LINE__,M;call exit_p(-1);endif
#define assert(E) if(.not.(E))then;write(stderr,'(a,":",i0,1x,"assert(E) failed")')__FILE__,__LINE__;call exit_p(-1);endif

#define H5_CLASS_INIT assertm(.not.self%initialised,'double init');if(uninitialised)call h5_init();self%initialised=.true.

#ifdef DEBUG
#define debug_print write(stderr,*) __FILE__//":",__LINE__,
#else
#define debug_print !
#endif

#define herr_t  integer(c_int)
#define hid_t   integer(c_int64_t)
#define hsize_t integer(hsize_t_kind)
#define hbool_t logical(c_bool)
#define htri_t  integer(c_int)
#define ssize_t integer(c_size_t)
! #define write_c(a) write(size(a), c_loc(a))

    module h5

    use iso_c_binding
    use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
    use mpi, only : mpi_comm_null, mpi_info_null, mpi_info_create, mpi_info_free

    implicit none

    private

    integer, parameter :: hsize_t_kind = c_long_long


    hid_t, public, parameter :: h5i_invalid_hid     = -1
    hid_t, public, parameter :: h5p_default         = 0
    hid_t, public, parameter :: h5s_all             = 0
    hid_t, public, parameter :: h5s_max_rank        = 32
    integer(c_int), public, parameter :: h5f_acc_rdonly      = 0
    integer(c_int), public, parameter :: h5f_acc_rdwr        = 1
    integer(c_int), public, parameter :: h5f_acc_trunc       = 2
    integer(c_int), public, parameter :: h5f_acc_excl        = 4
!     integer(c_int), public, parameter :: h5f_acc_creat       = 16
    integer(c_int), public, parameter :: h5f_acc_swmr_write  = 32
    integer(c_int), public, parameter :: h5f_acc_swmr_read   = 64
!     integer(c_int), public, parameter :: h5f_acc_default     = 65535

    integer(c_int), public, parameter :: h5f_obj_file        = 1
    integer(c_int), public, parameter :: h5f_obj_dataset     = 2
    integer(c_int), public, parameter :: h5f_obj_group       = 4
    integer(c_int), public, parameter :: h5f_obj_datatype    = 8
    integer(c_int), public, parameter :: h5f_obj_attr        = 16
    integer(c_int), public, parameter :: h5f_obj_all         = &
            ior(h5f_obj_file,ior(h5f_obj_dataset,ior(h5f_obj_group,ior(h5f_obj_datatype,h5f_obj_attr))))
    integer(c_int), public, parameter :: h5f_obj_local       = 32

    enum, bind(c)
        enumerator :: h5s_no_class = -1   ! error
        enumerator :: h5s_scalar   =  0   ! scalar variable
        enumerator :: h5s_simple   =  1   ! simple dataspace
        enumerator :: h5s_null     =  2   ! null dataspace
    end enum
    integer, parameter :: h5s_class_t_kind = kind(h5s_null)
#define h5s_class_t integer(kind=h5s_class_t_kind)

    enum, bind(c)
        enumerator :: h5s_select_noop      = -1  ! error
        enumerator :: h5s_select_set       = 0   ! Select "set" operation
        enumerator :: h5s_select_or              ! Binary "or" operation for hyperslabs
        enumerator :: h5s_select_and             ! Binary "and" operation for hyperslabs
        enumerator :: h5s_select_xor             ! Binary "xor" operation for hyperslabs
        enumerator :: h5s_select_notb            ! Binary "not" operation for hyperslabs
        enumerator :: h5s_select_nota            ! Binary "not" operation for hyperslabs
        enumerator :: h5s_select_append          ! Append elements to end of point selection
        enumerator :: h5s_select_prepend         ! Prepend elements to beginning of point selection
        enumerator :: h5s_select_invalid         ! Invalid upper bound on selection operations
    end enum
    integer, parameter :: h5s_seloper_t_kind = kind(h5s_select_noop)
#define h5s_seloper_t integer(kind=h5s_seloper_t_kind)

    enum, bind(c)
        enumerator :: h5s_sel_error       = -1  !  Error
        enumerator :: h5s_sel_none        = 0   !  Nothing selected
        enumerator :: h5s_sel_points      = 1   !  Sequence of points selected
        enumerator :: h5s_sel_hyperslabs  = 2   !  "New-style" hyperslab selection defined
        enumerator :: h5s_sel_all         = 3   !  Entire extent selected
        enumerator :: h5s_sel_n                 ! THIS MUST BE LAST
    end enum
    integer, parameter :: h5s_sel_type_kind = kind(h5s_sel_error)
#define h5s_sel_type integer(kind=h5s_sel_type_kind)

    enum, bind(c)
        enumerator :: h5t_no_class  = -1  ! error
        enumerator :: h5t_integer   = 0   ! integer types
        enumerator :: h5t_float     = 1   ! floating-point types
        enumerator :: h5t_time      = 2   ! date and time types
        enumerator :: h5t_string    = 3   ! character string types
        enumerator :: h5t_bitfield  = 4   ! bit field types
        enumerator :: h5t_opaque    = 5   ! opaque types
        enumerator :: h5t_compound  = 6   ! compound types
        enumerator :: h5t_reference = 7   ! reference types
        enumerator :: h5t_enum      = 8   ! enumeration types
        enumerator :: h5t_vlen      = 9   ! Variable-Length types
        enumerator :: h5t_array     = 10  ! Array types
    end enum
    integer, parameter :: h5t_class_t_kind = kind(h5t_no_class)
#define h5t_class_t integer(kind=h5t_class_t_kind)

    enum, bind(c)
        enumerator :: h5f_scope_local  = 0 ! specified file handle only
        enumerator :: h5f_scope_global = 1 ! entire virtual file
    end enum
    integer, parameter :: h5f_scope_t_kind = kind(h5f_scope_global)
#define h5f_scope_t integer(kind=h5f_scope_t_kind)

    enum, bind(c)
        enumerator :: h5g_storage_type_unknown = -1 ! Unknown link storage type
        enumerator :: h5g_storage_type_symbol_table ! Links in group are stored with a "symbol table"
                                                    ! (this is sometimes called "old-style" groups)
        enumerator :: h5g_storage_type_compact      ! Links are stored in object header
        enumerator :: h5g_storage_type_dense        ! Links are stored in fractal heap & indexed with v2 B-tree
    end enum
    integer, parameter :: h5g_storage_type_t_kind = kind(h5g_storage_type_unknown)
#define h5g_storage_type_t integer(kind=h5g_storage_type_t_kind)

    enum, bind(c)
        enumerator :: h5fd_mpio_independent = 0 ! zero is the default
        enumerator :: h5fd_mpio_collective
    end enum
    integer, parameter :: h5fd_mpio_xfer_t_kind = kind(h5fd_mpio_independent)
#define h5fd_mpio_xfer_t integer(kind=h5fd_mpio_xfer_t_kind)

    type, bind(c) :: h5g_info_t
        h5g_storage_type_t :: storage_type ! Type of storage for links in group
        hsize_t :: nlinks               ! Number of links in group
        integer(c_int64_t) :: max_corder           ! Current max. creation order value for group
        hbool_t ::    mounted              ! Whether group has a file mounted on it
    end type


    hid_t, public :: h5t_native_int8   = h5t_no_class
    hid_t, public :: h5t_native_int16  = h5t_no_class
    hid_t, public :: h5t_native_int32  = h5t_no_class
    hid_t, public :: h5t_native_int64  = h5t_no_class
    hid_t, public :: h5t_native_float  = h5t_no_class
    hid_t, public :: h5t_native_double = h5t_no_class
    hid_t, public :: h5t_std_i8be      = h5t_no_class
    hid_t, public :: h5t_std_i8le      = h5t_no_class
    hid_t, public :: h5t_std_i16be     = h5t_no_class
    hid_t, public :: h5t_std_i16le     = h5t_no_class
    hid_t, public :: h5t_std_i32be     = h5t_no_class
    hid_t, public :: h5t_std_i32le     = h5t_no_class
    hid_t, public :: h5t_std_i64be     = h5t_no_class
    hid_t, public :: h5t_std_i64le     = h5t_no_class
    hid_t, public :: h5t_ieee_f32be    = h5t_no_class
    hid_t, public :: h5t_ieee_f32le    = h5t_no_class
    hid_t, public :: h5t_ieee_f64be    = h5t_no_class
    hid_t, public :: h5t_ieee_f64le    = h5t_no_class

! beware the following real and cmplx are made to resemble real(kind) and complex(kind)
! while the above hdf5 defined types are all based on length in bits, for example  h5t_native_int8 is integer(1) in fortran
    hid_t, public :: h5t_native_real4  = h5t_no_class
    hid_t, public :: h5t_native_real8  = h5t_no_class
    hid_t, public :: h5t_native_cmplx4 = h5t_no_class
    hid_t, public :: h5t_native_cmplx8 = h5t_no_class


    hid_t, public :: h5p_attribute_create = h5p_default
    hid_t, public :: h5p_dataset_access   = h5p_default
    hid_t, public :: h5p_dataset_create   = h5p_default
    hid_t, public :: h5p_dataset_xfer     = h5p_default
    hid_t, public :: h5p_datatype_access  = h5p_default
    hid_t, public :: h5p_datatype_create  = h5p_default
    hid_t, public :: h5p_file_access      = h5p_default
    hid_t, public :: h5p_file_create      = h5p_default
    hid_t, public :: h5p_file_mount       = h5p_default
    hid_t, public :: h5p_group_access     = h5p_default
    hid_t, public :: h5p_group_create     = h5p_default
    hid_t, public :: h5p_link_access      = h5p_default
    hid_t, public :: h5p_link_create      = h5p_default
    hid_t, public :: h5p_object_copy      = h5p_default
    hid_t, public :: h5p_object_create    = h5p_default
    hid_t, public :: h5p_string_create    = h5p_default


    interface
! void h5_ftypes_init(size_t n, hid_t *types)
        subroutine h5_ftypes_init(n, types) bind(c, name='h5_ftypes_init')
            import c_size_t, c_int64_t
            implicit none
            integer(c_size_t), intent(in), value :: n
            hid_t, intent(out) :: types(*)
        end subroutine h5_ftypes_init

! herr_t H5open(void)
        herr_t function h5open() bind(c, name='H5open')
            import c_int
            implicit none
        end function h5open


! hid_t H5Fcreate( const char *name, unsigned flags, hid_t fcpl_id, hid_t fapl_id )
        hid_t function h5fcreate(name, flags, fcpl_id, fapl_id) bind(c, name='H5Fcreate')
            import c_int, c_char, c_int64_t
            implicit none
            character(c_char), intent(in) :: name(*)
            integer(c_int), value :: flags
            hid_t, value :: fcpl_id
            hid_t, value :: fapl_id
        end function h5fcreate

! hid_t H5Fopen( const char *name, unsigned flags, hid_t fapl_id )
        hid_t function h5fopen(name, flags, fapl_id) bind(c, name='H5Fopen')
            import c_int, c_char, c_int64_t
            implicit none
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), value :: flags
            hid_t, value :: fapl_id
        end function h5fopen

! herr_t H5Fclose( hid_t file_id )
        herr_t function h5fclose(file_id) bind(c, name='H5Fclose')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: file_id
        end function h5fclose

! herr_t H5Fflush(hid_t object_id, H5F_scope_t scope )
        herr_t function h5fflush(object_id, scope) bind(c, name='H5Fflush')
            import c_int, c_int64_t, c_char, h5f_scope_t_kind
            implicit none
            hid_t, intent(in), value :: object_id
            h5f_scope_t, intent(in), value :: scope
        end function h5fflush

! ssize_t H5Fget_obj_count( hid_t file_id, unsigned int types )
        ssize_t function h5fget_obj_count(file_id, types) bind(c, name='H5Fget_obj_count')
            import c_int64_t, c_int, c_size_t
            implicit none
            hid_t, intent(in), value :: file_id
            integer(c_int), intent(in), value :: types
        end function h5fget_obj_count

! htri_t H5Fis_hdf5(const char *name )
        htri_t function h5fis_hdf5(name) bind(c, name='H5Fis_hdf5')
            import c_char, c_int
            implicit none
            character(kind=c_char), intent(in) :: name(*)
        end function h5fis_hdf5


! hid_t H5Screate( H5S_class_t type )
        hid_t function h5screate(typ) bind(c, name='H5Screate')
            import c_int64_t, h5s_class_t_kind
            implicit none
            h5s_class_t, intent(in), value :: typ
        end function h5screate

! hid_t H5Screate_simple( int rank, const hsize_t * current_dims, const hsize_t * maximum_dims )
        hid_t function h5screate_simple(rank, current_dims, maximum_dims) bind(c, name='H5Screate_simple')
            import c_int, c_int64_t, hsize_t_kind !, c_ptr
            implicit none
            integer(c_int), intent(in), value :: rank
            hsize_t, intent(in) :: current_dims(*)
            hsize_t, intent(in) :: maximum_dims(*)
!             type(c_ptr), value :: maximum_dims ! this is to enable passing c_null_ptr
!             type(*), intent(in) :: maximum_dims(*) ! does not work as it does not allow passing of c_null_ptr
        end function h5screate_simple

! herr_t H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count, const hsize_t *block )
        herr_t function h5sselect_hyperslab(space_id, op, start, stride, count, block) &
            bind(c, name='H5Sselect_hyperslab')
            import c_int, c_int64_t, hsize_t_kind, h5s_seloper_t_kind
            implicit none
            hid_t, intent(in), value :: space_id
            h5s_seloper_t, intent(in), value :: op
            hsize_t, intent(in) :: start (*)
            hsize_t, intent(in) :: stride(*)
            hsize_t, intent(in) :: count (*)
            hsize_t, intent(in) :: block (*)
        end function h5sselect_hyperslab

! herr_t H5Sclose( hid_t space_id )
        herr_t function h5sclose(space_id) bind(c, name='H5Sclose')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: space_id
        end function h5sclose

! H5S_class_t H5Sget_simple_extent_type( hid_t space_id )
        h5s_class_t function h5sget_simple_extent_type(space_id) &
            bind(c, name='H5Sget_simple_extent_type')
            import c_int64_t, h5s_class_t_kind
            implicit none
            hid_t, intent(in), value :: space_id
        end function h5sget_simple_extent_type

! herr_t H5Sset_extent_simple( hid_t space_id, int rank, const hsize_t *current_size, const hsize_t *maximum_size )
        herr_t function h5sset_extent_simple(space_id, rank, current_size, maximum_size) &
            bind(c, name='H5Sset_extent_simple')
            import c_int, c_int64_t, hsize_t_kind
            implicit none
            hid_t, intent(in), value :: space_id
            integer(c_int), intent(in), value :: rank
            hsize_t, intent(in) :: current_size(*)
            hsize_t, intent(in) :: maximum_size(*)
        end function h5sset_extent_simple

! int H5Sget_simple_extent_ndims( hid_t space_id )
        integer(c_int) function h5sget_simple_extent_ndims(space_id) bind(c, name='H5Sget_simple_extent_ndims')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: space_id
        end function h5sget_simple_extent_ndims

! int H5Sget_simple_extent_dims(hid_t space_id, hsize_t *dims, hsize_t *maxdims )
        integer(c_int) function h5sget_simple_extent_dims(space_id, dims, maxdims) &
            bind(c, name='H5Sget_simple_extent_dims')
            import c_int, c_int64_t, hsize_t_kind
            implicit none
            hid_t, intent(in), value :: space_id
            hsize_t, intent(out) :: dims(*)
            hsize_t, intent(out) :: maxdims(*)
        end function h5sget_simple_extent_dims

! htri_t H5Sextent_equal( hid_t space1_id, hid_t space2_id )
        htri_t function h5sextent_equal(space1_id, space2_id) bind(c, name='H5Sextent_equal')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: space1_id
            hid_t, intent(in), value :: space2_id
        end function h5sextent_equal

! herr_t H5Sextent_copy(hid_t dest_space_id, hid_t source_space_id )
        herr_t function h5sextent_copy(dest_space_id, source_space_id) bind(c, name='H5Sextent_copy')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: dest_space_id
            hid_t, intent(in), value :: source_space_id
        end function h5sextent_copy

! herr_t H5Sget_regular_hyperslab( hid_t space_id, hsize_t start[], hsize_t stride[], hsize_t count[], hsize_t block[] )
        herr_t function h5sget_regular_hyperslab(space_id, start, stride, count, block) &
            bind(c, name='H5Sget_regular_hyperslab')
            import c_int64_t, c_int, hsize_t_kind
            implicit none
            hid_t, intent(in), value :: space_id
            hsize_t, intent(out) :: start (*)
            hsize_t, intent(out) :: stride(*)
            hsize_t, intent(out) :: count (*)
            hsize_t, intent(out) :: block (*)
        end function h5sget_regular_hyperslab

! htri_t H5Sis_regular_hyperslab( hid_t space_id )
        htri_t function h5sis_regular_hyperslab(space_id) bind(c, name='H5Sis_regular_hyperslab')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: space_id
        end function h5sis_regular_hyperslab

! H5S_sel_type H5Sget_select_type(hid_t space_id)
        h5s_sel_type function h5sget_select_type(space_id) bind(c, name='H5Sget_select_type')
            import c_int64_t, h5s_sel_type_kind
            implicit none
            hid_t, intent(in), value :: space_id
        end function h5sget_select_type


! hid_t H5Dcreate2( hid_t loc_id, const char *name, hid_t dtype_id, hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id)
        hid_t function h5dcreate2(loc_id, name, dtype_id, space_id, lcpl_id, dcpl_id, dapl_id) &
            bind(c, name='H5Dcreate2') ! the expansion of hid_t can push the line length beyond 132
            import c_int64_t, c_char
            implicit none
            hid_t, intent(in), value ::  loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: dtype_id
            hid_t, intent(in), value :: space_id
            hid_t, intent(in), value :: lcpl_id
            hid_t, intent(in), value :: dcpl_id
            hid_t, intent(in), value :: dapl_id
        end function h5dcreate2

! hid_t H5Dopen2( hid_t loc_id, const char *name, hid_t dapl_id )
        hid_t function h5dopen2(loc_id, name, dapl_id) bind(c, name='H5Dopen2')
            import c_int, c_char, c_int64_t
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: dapl_id
        end function h5dopen2

! herr_t H5Dread( hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void * buf )
        herr_t function h5dread(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf) &
            bind(c, name='H5Dread')
            import c_int, c_ptr, c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
            hid_t, intent(in), value :: mem_type_id
            hid_t, intent(in), value :: mem_space_id
            hid_t, intent(in), value :: file_space_id
            hid_t, intent(in), value :: xfer_plist_id
!             type(*), intent(inout) :: buf(..) ! a whole strct is passeed instead of the raw aray as in (*)
!             type(*), intent(inout) :: buf(*) ! ifort (up to v19.0.3) refuses to send (..) to this (gcc is ok)
            type(c_ptr), value :: buf
        end function h5dread

! herr_t H5Dwrite( hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void * buf )
        herr_t function h5dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf) &
            bind(c, name='H5Dwrite')
            import c_int, c_ptr, c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
            hid_t, intent(in), value :: mem_type_id
            hid_t, intent(in), value :: mem_space_id
            hid_t, intent(in), value :: file_space_id
            hid_t, intent(in), value :: xfer_plist_id
!             type(*), intent(in) :: buf(..) ! a whole strct is passeed instead of the raw aray as in (*)
!             type(*), intent(inout) :: buf(*) ! ifort (up to v19.0.3) refuses to send (..) to this (gcc is ok)
            type(c_ptr), intent(in), value :: buf
        end function h5dwrite

! herr_t H5Dclose(hid_t dataset_id )
        herr_t function h5dclose(dataset_id) bind(c, name='H5Dclose')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
        end function h5dclose

! herr_t H5Dflush(hid_t dataset_id )
        herr_t function h5dflush(dataset_id) bind(c, name='H5Dflush')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
        end function h5dflush

! hid_t H5Dget_space( hid_t dataset_id )
        hid_t function h5dget_space(dataset_id) bind(c, name='H5Dget_space')
            import c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
        end function h5dget_space

! hid_t H5Dget_type(hid_t dataset_id )
        hid_t function h5dget_type(dataset_id) bind(c, name='H5Dget_type')
            import c_int64_t
            implicit none
            hid_t, intent(in), value :: dataset_id
        end function h5dget_type


! hid_t H5Gcreate2( hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id )
        hid_t function h5gcreate2(loc_id, name, lcpl_id, gcpl_id, gapl_id) bind(c, name='H5Gcreate2')
            import c_int, c_int64_t, c_char
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: lcpl_id
            hid_t, intent(in), value :: gcpl_id
            hid_t, intent(in), value :: gapl_id
        end function h5gcreate2

! herr_t H5Gclose(hid_t group_id)
        herr_t function h5gclose(group_id) bind(c, name='H5Gclose')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: group_id
        end function h5gclose

! herr_t H5Gflush(hid_t group_id)
        herr_t function h5gflush(group_id) bind(c, name='H5Gflush')
            import c_int, c_int64_t
            implicit none
            hid_t, intent(in), value :: group_id
        end function h5gflush

! herr_t H5Gget_info( hid_t group_id, H5G_info_t *group_info )
!         herr_t function h5gget_info( hid_t group_id, H5G_info_t *group_info ) bind(c, name='H5Gget_info')
!         end function h5gget_info

! herr_t H5Gget_info_by_name( hid_t loc_id, const char *group_name, H5G_info_t *group_info, hid_t lapl_id )
        herr_t function h5gget_info_by_name(loc_id, group_name, group_info, lapl_id) &
            bind(c, name='H5Gget_info_by_name')
            import c_int, c_int64_t, c_char, h5g_info_t
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: group_name(*)
            type(h5g_info_t), intent(inout) :: group_info
            hid_t, intent(in), value :: lapl_id
        end function h5gget_info_by_name

! hid_t H5Gopen2( hid_t loc_id, const char * name, hid_t gapl_id )
        hid_t function h5gopen2(loc_id, name, gapl_id) bind(c, name='H5Gopen2')
            import c_int, c_int64_t, c_char
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: gapl_id
        end function h5gopen2


! hid_t H5Tcreate( H5T_class_t class, size_t size )
        hid_t function h5tcreate(class_t, size) bind(c, name='H5Tcreate')
            import h5t_class_t_kind, c_int64_t, c_size_t
            implicit none
            h5t_class_t, intent(in), value :: class_t
            integer(c_size_t), intent(in), value :: size
        end function h5tcreate

! herr_t H5Tclose( hid_t dtype_id )
        herr_t function h5tclose(dtype_id) bind(c, name='H5Tclose')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: dtype_id
        end function h5tclose

! hid_t H5Tarray_create2( hid_t base_type_id, unsigned rank, const hsize_t dims[/*rank*/] )
        hid_t function h5tarray_create2(base_type_id, rank, dims) bind(c, name='H5Tarray_create2')
            import c_int64_t, c_int, hsize_t_kind
            implicit none
            hid_t, intent(in), value :: base_type_id
            integer(c_int), intent(in), value :: rank
            hsize_t, intent(in) :: dims(*)
        end function h5tarray_create2

! hid_t H5Tcopy( hid_t dtype_id )
        hid_t function h5tcopy(dtype_id) bind(c, name='H5Tcopy')
            import c_int64_t
            implicit none
            hid_t, intent(in), value :: dtype_id
        end function h5tcopy

! herr_t H5Tinsert( hid_t dtype_id, const char * name, size_t offset, hid_t field_id )
        herr_t function h5tinsert(dtype_id, name, offset, field_id) bind(c, name='H5Tinsert')
            import c_int, c_int64_t, c_char, c_size_t
            implicit none
            hid_t, value :: dtype_id
            character(kind=c_char), intent(in) :: name(*)
            integer(c_size_t), value :: offset
            hid_t, intent(in), value :: field_id
        end function h5tinsert

! htri_t H5Tequal( hid_t dtype_id1, hid_t dtype_id2 )
        htri_t function h5tequal(dtype_id1, dtype_id2) bind(c, name='H5Tequal')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: dtype_id1
            hid_t, intent(in), value :: dtype_id2
        end function h5tequal

! herr_t H5Tcommit2( hid_t loc_id, const char *name, hid_t dtype_id, hid_t lcpl_id, hid_t tcpl_id, hid_t tapl_id )
        herr_t function h5tcommit2(loc_id, name, dtype_id, lcpl_id, tcpl_id, tapl_id) bind(c, name='H5Tcommit2')
            import c_int, c_int64_t, c_char
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, value :: dtype_id
            hid_t, intent(in), value :: lcpl_id
            hid_t, intent(in), value :: tcpl_id
            hid_t, intent(in), value :: tapl_id
        end function h5tcommit2

! herr_t H5Tflush(hid_t dtype_id)
        herr_t function h5tflush(dtype_id) bind(c, name='H5Tflush')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: dtype_id
        end function h5tflush

! hid_t H5Topen2( hid_t loc_id, const char * name, hid_t tapl_id )
        hid_t function h5topen2(loc_id, name, tapl_id) bind(c, name='H5Topen2')
            import c_int, c_int64_t, c_char
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: tapl_id
        end function h5topen2


! htri_t H5Oexists_by_name( hid_t loc_id, const char * name, hid_t lapl_id )
        htri_t function h5oexists_by_name(loc_id, name, lapl_id) bind(c, name='H5Oexists_by_name')
            import c_int64_t, c_char, c_int
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: lapl_id
        end function h5oexists_by_name


! htri_t H5Lexists( hid_t loc_id, const char *name, hid_t lapl_id )
        htri_t function h5lexists(loc_id, name, lapl_id) bind(c, name='H5Lexists')
            import c_int64_t, c_char, c_int
            implicit none
            hid_t, intent(in), value :: loc_id
            character(kind=c_char), intent(in) :: name(*)
            hid_t, intent(in), value :: lapl_id
        end function h5lexists


! hid_t H5Pcreate( hid_t cls_id )
        hid_t function h5pcreate(cls_id) bind(c, name='H5Pcreate')
            import c_int64_t
            implicit none
            hid_t, value :: cls_id
        end function h5pcreate

! herr_t H5Pclose(hid_t plist )
        herr_t function h5pclose(plist) bind(c, name='H5Pclose')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: plist
        end function h5pclose

! herr_t H5Pset_fapl_mpio_foo(hid_t fapl_id, MPI_Fint comm, MPI_Fint info)
        herr_t function h5pset_fapl_mpio_foo(fapl_id, comm, info) &
            bind(c, name='H5Pset_fapl_mpio_foo')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: fapl_id
            integer(c_int), intent(in), value :: comm
            integer(c_int), value :: info
        end function h5pset_fapl_mpio_foo

! herr_t H5Pset_dxpl_mpio( hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode )
        herr_t function h5pset_dxpl_mpio(dxpl_id, xfer_mode) bind(c, name='H5Pset_dxpl_mpio')
            import c_int64_t, c_int, h5fd_mpio_xfer_t_kind
            implicit none
            hid_t, intent(in), value :: dxpl_id
            h5fd_mpio_xfer_t, intent(in), value :: xfer_mode
        end function h5pset_dxpl_mpio


! htri_t H5Iis_valid( hid_t obj_id )
        htri_t function h5iis_valid(obj_id) bind(c, name='H5Iis_valid')
            import c_int64_t, c_int
            implicit none
            hid_t, intent(in), value :: obj_id
        end function h5iis_valid

    end interface

    interface
! #include <unistd.h>
! int access(const char *path, int amode);
!         integer(c_int) function access(const char *path, int amode) bind(c)
!         end function access
! let's just use f's inquire()

! void exit(int status);
        subroutine exit_p(rc) bind(c, name='exit')
            import c_int
            integer(c_int), intent(in), value :: rc
        end subroutine exit_p
! f's stop prints an ugly 'STOP' in gfortran and ifort

    end interface


    type, public :: h5space
        logical :: initialised = .false.
        hid_t :: id = h5i_invalid_hid
!         h5s_class_t :: typ = h5s_null
        h5s_class_t :: typ = h5s_simple
        integer(c_int) :: rank = 0
        hsize_t :: current_dims(h5s_max_rank) = 0
        hsize_t :: maximum_dims(h5s_max_rank) = 0

        contains

        procedure, private :: h5space_init_m
        procedure, private :: h5space_init_i
        generic :: init => h5space_init_m, h5space_init_i

        procedure, private :: h5space_resize_m
        procedure, private :: h5space_resize_s
        procedure, private :: h5space_resize_i
        generic :: resize => h5space_resize_m, h5space_resize_s, h5space_resize_i

        procedure, private :: h5space_select_hyperslab_m
        procedure, private :: h5space_select_hyperslab_s
        procedure, private :: h5space_select_hyperslab_i
        generic :: select_hyperslab => h5space_select_hyperslab_m, h5space_select_hyperslab_s, h5space_select_hyperslab_i

!         generic :: find => find_h5dtype_int1, find_h5dtype_int2, find_h5dtype_int4, find_h5dtype_int8, &
!              find_h5dtype_real4, find_h5dtype_real8, find_h5dtype_cmplx4, find_h5dtype_cmplx8
!         procedure, private :: find_h5dtype_int1, find_h5dtype_int2, find_h5dtype_int4, find_h5dtype_int8, &
!              find_h5dtype_real4, find_h5dtype_real8, find_h5dtype_cmplx4, find_h5dtype_cmplx8

        procedure, pass :: close => h5space_close
        final :: h5space_clear
    end type

    type, public, abstract :: h5location
        logical :: initialised = .false.
        hid_t :: id = h5i_invalid_hid
        hid_t :: typ = h5p_default
    end type

    type, public, extends(h5location) :: h5file
        character(len=:), allocatable :: name
        integer(c_int) :: flags
        hid_t :: fcpl_id = h5p_default
        hid_t :: fapl_id = h5p_default
        integer :: comm = mpi_comm_null

        contains
        procedure, pass :: init  => h5file_init
        procedure, pass :: read  => h5file_read
        procedure, pass :: write => h5file_write
        procedure, pass :: close => h5file_close
        final :: h5file_clear
    end type

    type, public, extends(h5location) :: h5dataset
!         hid_t :: id = h5i_invalid_hid
!         type(h5location), pointer :: location
!         type(h5location) :: location
!         class(*), pointer :: location
        character(len=:), allocatable :: name
        hid_t :: dtype = h5t_no_class
        type(h5space), pointer :: space => null()
        logical :: space_allocated = .false.

        hid_t :: lcpl_id = h5p_default
        hid_t :: dcpl_id = h5p_default
        hid_t :: dapl_id = h5p_default

        contains
        procedure, pass :: init  => h5dataset_init
        procedure, pass :: read  => h5dataset_read
        procedure, pass :: write => h5dataset_write
        procedure, pass :: close => h5dataset_close
        final :: h5dataset_clear
    end type

    type, public, extends(h5location) :: h5group
        character(len=:), allocatable :: name
        contains
        procedure, pass :: init => h5group_init
        procedure, pass :: close => h5group_close
        final :: h5group_clear
    end type

!! better integrate the various properties in the respective classes instead of making a p-class for every base class.
!! This does not follow the library as closely as other classes follow it but it makes more sense since there is not much to unite the different properties under one roof and the amount of work to mak.
!     type, public, abstract :: h5property
!         hid_t :: id = h5p_default
!     end type

    logical :: uninitialised = .true.

    character, parameter :: path_sep = '/'

    interface find_h5dtype
        module procedure find_h5dtype_int1, find_h5dtype_int2, find_h5dtype_int4, find_h5dtype_int8, &
            find_h5dtype_real4, find_h5dtype_real8, find_h5dtype_cmplx4, find_h5dtype_cmplx8
    end interface


!     public :: h5_init, c_loc
    public :: find_h5dtype, exit_p, stderr
    public :: h5_read
    public :: h5_write

    contains

    subroutine h5_init()
        integer(c_size_t), parameter :: n = 34
        herr_t :: err

        hid_t :: types(n)

        assertm(uninitialised, 'double init')

        call h5_ftypes_init(n, types)

        h5t_native_int8      = types( 1)
        h5t_native_int16     = types( 2)
        h5t_native_int32     = types( 3)
        h5t_native_int64     = types( 4)
        h5t_native_float     = types( 5)
        h5t_native_double    = types( 6)
        h5t_std_i8be         = types( 7)
        h5t_std_i8le         = types( 8)
        h5t_std_i16be        = types( 9)
        h5t_std_i16le        = types(10)
        h5t_std_i32be        = types(11)
        h5t_std_i32le        = types(12)
        h5t_std_i64be        = types(13)
        h5t_std_i64le        = types(14)
        h5t_ieee_f32be       = types(15)
        h5t_ieee_f32le       = types(16)
        h5t_ieee_f64be       = types(17)
        h5t_ieee_f64le       = types(18)
        h5p_attribute_create = types(19)
        h5p_dataset_access   = types(20)
        h5p_dataset_create   = types(21)
        h5p_dataset_xfer     = types(22)
        h5p_datatype_access  = types(23)
        h5p_datatype_create  = types(24)
        h5p_file_access      = types(25)
        h5p_file_create      = types(26)
        h5p_file_mount       = types(27)
        h5p_group_access     = types(28)
        h5p_group_create     = types(29)
        h5p_link_access      = types(30)
        h5p_link_create      = types(31)
        h5p_object_copy      = types(32)
        h5p_object_create    = types(33)
        h5p_string_create    = types(34)

        h5t_native_real4  = h5t_native_float
        h5t_native_real8  = h5t_native_double

#ifndef ARRAY_CMPLX
        h5t_native_cmplx4 = h5tcreate(h5t_compound, 2*4_c_size_t)
        err = h5tinsert(h5t_native_cmplx4, 'r'//c_null_char, 0_c_size_t, h5t_native_real4)
        err = h5tinsert(h5t_native_cmplx4, 'i'//c_null_char, 4_c_size_t, h5t_native_real4)

        h5t_native_cmplx8 = h5tcreate(h5t_compound, 2*8_c_size_t)
        err = h5tinsert(h5t_native_cmplx8, 'r'//c_null_char, 0_c_size_t, h5t_native_real8)
        err = h5tinsert(h5t_native_cmplx8, 'i'//c_null_char, 8_c_size_t, h5t_native_real8)
#else
        err = 0
        h5t_native_cmplx4 = h5tarray_create2(h5t_native_real4, 1, [2_hsize_t_kind])
        h5t_native_cmplx8 = h5tarray_create2(h5t_native_real8, 1, [2_hsize_t_kind])
#endif

! needs some mpi awareness
!         write(stdout, '(a)') H5FOO

        uninitialised = .false.
    end subroutine h5_init

#define find_h5dtype_m(name, ft, ht) hid_t function name(v)result(t);ft,intent(in)::v(..);if(uninitialised)call h5_init();t=ht;end function name

    find_h5dtype_m(find_h5dtype_int1, integer(1), h5t_native_int8)
    find_h5dtype_m(find_h5dtype_int2, integer(2), h5t_native_int16)
    find_h5dtype_m(find_h5dtype_int4, integer(4), h5t_native_int32)
    find_h5dtype_m(find_h5dtype_int8, integer(8), h5t_native_int64)
    find_h5dtype_m(find_h5dtype_real4, real(4), h5t_native_real4)
    find_h5dtype_m(find_h5dtype_real8, real(8), h5t_native_real8)
    find_h5dtype_m(find_h5dtype_cmplx4, complex(4), h5t_native_cmplx4)
    find_h5dtype_m(find_h5dtype_cmplx8, complex(8), h5t_native_cmplx8)


! when class(*), dimension(..) or (*) is implemented in compilers (gcc 8 and intel 19 do not have it), the preprocessor madness above shall be replaced with the following
!     hid_t function find_h5dtype(v) result (dtype)
!         class(*), intent(in) ::  v(..)
!
!         if (uninitialised) call h5_init()
!
!         select type (v)
! !         type is (real(16))
! !             dtype = h5t_no_class
!         type is (real(8))
!             dtype = h5t_native_real8
!         type is (real(4))
!             dtype = h5t_native_real4
!         type is (integer(8))
!             dtype = h5t_native_int64
!         type is (integer(4))
!             dtype = h5t_native_int32
!         type is (integer(2))
!             dtype = h5t_native_int16
!         type is (integer(1))
!             dtype = h5t_native_int8
! !         type is (complex(16))
!         type is (complex(8))
!             dtype = h5t_native_cmplx8
!         type is (complex(4))
!             dtype = h5t_native_cmplx4
! !         type is (logical)
!         class default
!             assertm(.false.,'type not found')
!         end select
!     end function find_h5dtype

    function find_h5shape(v) result(s)
        type(*) :: v(..)
        integer, allocatable :: s(:)
        if (rank(v) == 0) then
            s = [1]
        else
            s = shape(v)
        end if
    end function find_h5shape

!     subroutine h5_get_path(path, fname, dname, vname)
! maybe simply return index instead of copying all these strings or probably best to reuse the mamory from path?
!         character(len=*), intent(in) :: path
!         character(len=*), intent(out), optional :: fname
!         character(len=*), intent(out), optional :: dname
!         character(len=*), intent(out), optional :: vname
!     end subroutine h5_get_path

    subroutine h5_mkdir(location, path)
! eventhough it seems to work, according to the h5lexists documentation links shall be verified one at a time, from the root up and not the other way around!
        class(h5location), intent(in) :: location
        character(len=*), intent(in) :: path

        type(h5group), allocatable :: g(:)
        integer :: ng, np, nc, i, idx(0:len_trim(path))

! this splitting has to be exported to a separate utility
        nc = len_trim(path)
        assertm(nc > 0, 'empty path supplied')

        np = 0
        do i = 1, nc
            if (path(i:i) == path_sep) then
                idx(np) = i-1
                np = np + 1
            end if
        end do
        idx(np) = nc ! add '/' at the end to avoid having making a special case for the last element
        np = np + 1

        allocate(g(0:np-1))

        if (idx(0) > 0) then
            call g(0) % init(location, path(1:idx(0)))
        else
            call g(0) % init(location, '/')
        end if

        ng = 1
        do i = 1, np - 1
            if (idx(i) > idx(i-1) + 1) then ! guard against repeating path_sep
                call g(ng) % init(g(ng-1), path(idx(i-1)+2:idx(i)))
                ng = ng + 1
            end if
        end do

! not sure if the closing needs to be in reverse order
        do i = ng-1, 0, -1
            call g(i) % close()
        end do

        deallocate(g)
    end subroutine h5_mkdir

    subroutine h5_read(path, data, dtype, fdims, mdims, file_space, mem_space, comm, mpio_xfer)
        character(len=*), intent(in) :: path
        type(*), intent(inout) :: data(..)
        hid_t, intent(in) :: dtype
        integer, intent(in), optional :: fdims(:)
        integer, intent(in), optional :: mdims(:)
        type(h5space), intent(in), optional :: file_space ! optional if the dataset exists
        type(h5space), intent(in), optional :: mem_space
        integer, intent(in), optional :: comm
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        type(h5file) :: f
        integer :: idx

        idx = index(path, ':')
        assertm(idx > 1, 'filename part of argument path="'//trim(path)//'" not valid')

        call f % init(path(1:idx-1), flags = h5f_acc_rdonly, comm = comm)
        call f % read(path(idx+1:len_trim(path)), data, dtype, fdims, mdims, file_space, mem_space, mpio_xfer)
    end subroutine h5_read

! would have been nice if dimension(*) worked properly in generics or class(*), dimension(*) was implemented, till then dtype and file_space shall be mandatory
    subroutine h5_write(path, data, dtype, fdims, mdims, file_space, mem_space, comm, mpio_xfer)
        character(len=*), intent(in) :: path
        type(*), intent(in) :: data(..)
!         hid_t, intent(in), optional :: dtype ! when gcc get their act together and support CFI_Fortran_binding.h
        hid_t, intent(in) :: dtype
        integer, intent(in), optional :: fdims(:)
        integer, intent(in), optional :: mdims(:)
        type(h5space), intent(in), optional :: file_space ! optional if the dataset exists
        type(h5space), intent(in), optional :: mem_space
        integer, intent(in), optional :: comm
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        type(h5file) :: f
        integer :: idx

        idx = index(path, ':')
        assertm(idx > 1, 'filename part of argument path="'//trim(path)//'" not valid')

        call f % init(path(1:idx-1), comm = comm)
        call f % write(path(idx+1:len_trim(path)), data, dtype, fdims, mdims, file_space, mem_space, mpio_xfer)
    end subroutine h5_write


!%%%%%%%%%%%%%
    subroutine h5file_init(self, name, flags, comm)
        class(h5file), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, intent(in), optional :: flags
        integer, intent(in), optional :: comm
        herr_t :: err
        htri_t :: stat
        logical :: exists
        integer :: info

        H5_CLASS_INIT

        self % typ = h5f_obj_file
        self % name = trim(name)//c_null_char

        if (present(comm)) self % comm = comm
        if (self % comm /= mpi_comm_null) then
            self % fapl_id = h5pcreate(h5p_file_access)
            call mpi_info_create(info, err)
            err = h5pset_fapl_mpio_foo(self % fapl_id, self % comm, info)
            assertm(err > -1, 'h5pset_fapl_mpio_foo failed')
            call mpi_info_free(info, err)
        end if

        inquire(file = self % name, exist = exists)
        if (present(flags)) then
            assertm(iand(h5f_acc_rdonly, flags) == 0 .or. exists, "file "//trim(name)//" not found")
        end if

        if (exists) then
            stat = h5fis_hdf5(self % name)
            if (stat == 0) then
                write(*,*) 'file '//self % name//' exists and is not of hdf5 type'
                return
            end if
            self % flags = h5f_acc_rdwr
            if (present(flags)) self % flags = flags
            self % id = h5fopen(self % name, self % flags, self % fapl_id)
        else
            self % flags = h5f_acc_trunc
            if (present(flags)) self % flags = flags
            self % id = h5fcreate(self % name, self % flags, self % fcpl_id, self % fapl_id)
        end if
        debug_print 'h5file_init', self % id, self % name
    end subroutine h5file_init

    subroutine h5file_read(self, path, data, dtype, fdims, mdims, file_space, mem_space, mpio_xfer)
        class(h5file), intent(inout) :: self
        character(len=*), intent(in) :: path
        type(*), intent(inout) :: data(..)
        hid_t, intent(in) :: dtype
!         hid_t, intent(in), optional :: dtype
        integer, intent(in), optional :: fdims(:)
        integer, intent(in), optional :: mdims(:)
        type(h5space), intent(in), optional :: file_space ! optional if the dataset exists
        type(h5space), intent(in), optional :: mem_space
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        type(h5dataset) :: d

        call d % init(self, path, data, dtype, fdims, file_space)
        call d % read(data, mdims, fdims, mem_space, file_space, mpio_xfer) ! file_space sent again because the selection is not copied on init if the dataset already exists in the file.
    end subroutine h5file_read

! would have been nice if dimension(*) worked properly in generics or class(*), dimension(*) was implemented, till then dtype and file_space shall be mandatory
    subroutine h5file_write(self, path, data, dtype, fdims, mdims, file_space, mem_space, mpio_xfer)
        class(h5file), intent(inout) :: self
        character(len=*), intent(in) :: path
        type(*), intent(in) :: data(..)
!         hid_t, intent(in), optional :: dtype ! hopefully soon
        hid_t, intent(in) :: dtype
        integer, intent(in), optional :: fdims(:) ! either this or the space types have to be passed
        integer, intent(in), optional :: mdims(:)
        type(h5space), intent(in), optional :: file_space ! optional if the dataset exists
        type(h5space), intent(in), optional :: mem_space
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        type(h5dataset) :: d
!         herr_t :: err

        call d % init(self, path, data, dtype, fdims, file_space)
        call d % write(data, mdims, fdims, mem_space, file_space, mpio_xfer) ! file_space sent again because the selection is not copied on init if the dataset already exists in the file.
!         err = h5fflush(self % id, h5f_scope_local)
    end subroutine h5file_write

    subroutine h5file_close(self)
        class(h5file), intent(inout) :: self
        herr_t :: err
        debug_print 'h5file_close', self % id, self % name
        if (self % fapl_id /= h5p_default .and. h5iis_valid(self % fapl_id) > 0) then
            err = h5pclose(self % fapl_id)
            assertm(err > -1, 'h5pclose failed')
        end if
        err = h5fclose(self%id)
        assertm(err > -1, 'h5fclose failed')
    end subroutine h5file_close

    subroutine h5file_clear(self)
        type(h5file), intent(inout) :: self
        debug_print 'h5file_clear', self % id, self % name
        if (h5iis_valid(self % id) <= 0) return
        if (h5fget_obj_count(self % id, h5f_obj_file) > 0) call self % close()
        self % id = h5i_invalid_hid
    end subroutine h5file_clear


!%%%%%%%%%%%%%
    subroutine h5space_init_m(self, current_dims, maximum_dims, start, strid, count, block)
        class(h5space), intent(inout) :: self
        integer, intent(in) :: current_dims(:)
        integer, intent(in), optional :: maximum_dims(size(current_dims))
        integer, intent(in), optional :: start(size(current_dims))
        integer, intent(in), optional :: strid(size(current_dims))
        integer, intent(in), optional :: count(size(current_dims))
        integer, intent(in), optional :: block(size(current_dims))

        H5_CLASS_INIT

        self % id = h5screate(self % typ)

        call self % resize(current_dims, maximum_dims)
        if (present(start) .or. present(strid) .or. present(count) .or. present(block)) &
            call self % select_hyperslab(start, strid, count, block)
    end subroutine h5space_init_m

    subroutine h5space_init_i(self, space_id)
        class(h5space), intent(inout) :: self
        hid_t, intent(in) :: space_id

        H5_CLASS_INIT

! only h5s_simple is supported for now
        assert(h5sget_simple_extent_type(space_id) == self % typ)

        self % id = space_id
        call self % resize(self % id)
    end subroutine h5space_init_i

    subroutine h5space_resize_s(self, space)
        class(h5space), intent(inout) :: self
        type(h5space), intent(in) :: space

        call self % resize(space % id)
    end subroutine h5space_resize_s

    subroutine h5space_resize_i(self, space_id)
        class(h5space), intent(inout) :: self
        hid_t, intent(in) :: space_id

        herr_t :: err

        self % rank = h5sget_simple_extent_dims(space_id, self % current_dims, self % maximum_dims)
        self % current_dims(1 : self % rank) = self % current_dims(self % rank : 1 : -1)
        self % maximum_dims(1 : self % rank) = self % maximum_dims(self % rank : 1 : -1)

        if (self % id /= space_id) then
            if (h5sextent_equal(self % id, space_id) < 1) then
                err = h5sextent_copy(self % id, space_id)
                assertm(err > -1, 'h5sextent_copy failed')
            end if
        end if
    end subroutine h5space_resize_i

    subroutine h5space_resize_m(self, current_dims, maximum_dims)
        class(h5space), intent(inout) :: self
        integer, intent(in) :: current_dims(:)
        integer, intent(in), optional :: maximum_dims(size(current_dims))

        hsize_t :: current_dims_c(size(current_dims))
        hsize_t :: maximum_dims_c(size(current_dims))
        herr_t :: err

        self % rank = size(current_dims)
        self % current_dims(1 : self % rank) = current_dims
        self % maximum_dims(1 : self % rank) = current_dims
        if (present(maximum_dims)) self % maximum_dims(1 : self % rank) = maximum_dims

        current_dims_c = self % current_dims(self % rank : 1 : -1)
        maximum_dims_c = self % maximum_dims(self % rank : 1 : -1)

        err = h5sset_extent_simple(self % id, self % rank, current_dims_c, maximum_dims_c)
        assertm(err > -1, 'h5sset_extent_simple')
    end subroutine h5space_resize_m

    subroutine h5space_select_hyperslab_s(self, space)
        class(h5space), intent(inout) :: self
        type(h5space), intent(in) :: space
        call self % select_hyperslab(space % id)
    end subroutine h5space_select_hyperslab_s

    subroutine h5space_select_hyperslab_i(self, space_id)
        class(h5space), intent(inout) :: self
        hid_t, intent(in) :: space_id

        hsize_t :: start_c(self % rank)
        hsize_t :: strid_c(self % rank)
        hsize_t :: count_c(self % rank)
        hsize_t :: block_c(self % rank)
        herr_t :: err

        if (h5sextent_equal(self % id, space_id) > 0) then
            if (h5sget_select_type(space_id) == h5s_sel_hyperslabs) then
                if (h5sis_regular_hyperslab(space_id) > 0) then
                    err = h5sget_regular_hyperslab(space_id, start_c, strid_c, count_c, block_c)
                    assertm(err > -1, 'h5sget_regular_hyperslab failed')
                    err = h5sselect_hyperslab(self % id, h5s_select_set, start_c, strid_c, count_c, block_c)
                    assertm(err > -1, 'h5sselect_hyperslab failed')
                end if
            end if
        end if
    end subroutine h5space_select_hyperslab_i

    subroutine h5space_select_hyperslab_m(self, start, strid, count, block, op)
        class(h5space), intent(inout) :: self
        h5s_seloper_t, intent(in), optional :: op
        integer, intent(in), optional :: start(self % rank)
        integer, intent(in), optional :: strid(self % rank)
        integer, intent(in), optional :: count(self % rank)
        integer, intent(in), optional :: block(self % rank)

        herr_t :: err

        h5s_seloper_t :: op_c
        hsize_t :: start_c(self % rank)
        hsize_t :: strid_c(self % rank)
        hsize_t :: count_c(self % rank)
        hsize_t :: block_c(self % rank)

        op_c = h5s_select_set
        start_c = 0
        strid_c = 1
        count_c = self % current_dims(self % rank : 1 : -1) ! would it not make more sense to swap the values of count and block?
        block_c = 1

        if (present(op)) op_c = op
        if (present(start)) start_c = start(self % rank : 1 : -1)
        if (present(strid)) strid_c = strid(self % rank : 1 : -1)
        if (present(count)) count_c = count(self % rank : 1 : -1)
        if (present(block)) block_c = block(self % rank : 1 : -1)

        err = h5sselect_hyperslab(self % id, op_c, start_c, strid_c, count_c, block_c)
        assertm(err > -1, 'h5sselect_hyperslab failed')
    end subroutine h5space_select_hyperslab_m

    subroutine h5space_close(self)
        class(h5space), intent(inout) :: self
        herr_t :: err
        debug_print 'h5space_close', self % id
        err = h5sclose(self % id)
        assertm(err > -1, 'h5sclose failed')
    end subroutine h5space_close

    subroutine h5space_clear(self)
        type(h5space), intent(inout) :: self
        debug_print 'h5space_clear', self % id
        if (h5iis_valid(self % id) <= 0) return
        call self % close()
        self % id = h5i_invalid_hid
    end subroutine h5space_clear


!%%%%%%%%%%%%%
    subroutine h5dataset_init(self, location, name, data, dtype, ddims, space, lcpl_id, dcpl_id, dapl_id)
        class(h5dataset), intent(inout) :: self
        class(h5location), intent(in) :: location
        character(len=*), intent(in) :: name
        type(*), intent(in), optional :: data(..) !only used for type inference in new dataset creation and if dtype not passed
        hid_t, intent(in), optional  :: dtype
        integer, intent(in), optional :: ddims(:) ! either this or the space types have to be passed
        type(h5space), intent(in), target, optional :: space
        hid_t, intent(in), optional :: lcpl_id
        hid_t, intent(in), optional :: dcpl_id
        hid_t, intent(in), optional :: dapl_id

        integer :: idx
        herr_t :: err

        H5_CLASS_INIT

        assertm(.not.(present(ddims).and.present(space)), 'both, ddims and space supplied.')

        self % typ = h5f_obj_dataset
        self % name = trim(name)//c_null_char

!         if (present(dtype))   self % dtype = dtype
!         if (present(space))   self % space => space
        if (present(lcpl_id)) self % lcpl_id = lcpl_id
        if (present(dcpl_id)) self % dcpl_id = dcpl_id
        if (present(dapl_id)) self % dapl_id = dapl_id

        idx = index(name, path_sep, back = .true.)
        if (idx > 1) call h5_mkdir(location, name(1:idx))
        err = h5lexists(location % id, self % name, h5p_default)
        assertm(err > -1, 'h5lexists failed')
        if (err > 0) then
            self % id = h5dopen2(location % id, self % name, self % dapl_id)
            self % dtype = h5dget_type(self % id)
            if (present(dtype)) then ! cannot bring the assert here because it is a macro starting with 'if'
                assert(h5tequal(self % dtype, dtype) > 0)
            end if

            allocate(self % space); self % space_allocated = .true.
            call self % space % init(h5dget_space(self % id))
            if (present(ddims)) call self % space % resize(ddims)
            if (present(space)) call self % space % resize(space)
            if (present(space)) call self % space % select_hyperslab(space)
        else
!             assertm(present(dtype).or.present(data), "h5dataset_init: at least of dtype or data needs to be supplied")
            assertm(present(dtype), "h5dataset_init: dtype must be supplied for new dataset creation")
! instantiate a copt of potentially immutable type untill i figure out how to distinguish them, h5tclose(immutable) is a fail
            if (present(dtype)) self % dtype = h5tcopy(dtype)
!             if (.not. present(dtype)) self % dtype = find_h5dtype(data)
            if (present(space)) then
                self % space => space; self % space_allocated = .false.
            else if (present(ddims)) then
                allocate(self % space); self % space_allocated = .true.
                call self % space % init(ddims)
            else if (present(data)) then
                allocate(self % space); self % space_allocated = .true.
                call self % space % init(find_h5shape(data))
            else
                assertm(.false., 'neither ddims nor space or data supplied.')
            end if
            self % id = h5dcreate2(location % id, self % name, self % dtype, self % space % id, &
                                                self % lcpl_id, self % dcpl_id, self % dapl_id)
        end if
        debug_print 'h5dataset_init', self % id, self % name
    end subroutine h5dataset_init

    subroutine h5dataset_read(self, buf, mdims, fdims, mem_space, file_space, mpio_xfer)
        class(h5dataset), intent(inout) :: self
        type(*), intent(inout), target :: buf(..)
        integer, intent(in), optional :: mdims(:)
        integer, intent(in), optional :: fdims(:)
        type(h5space), intent(in), optional :: mem_space
        type(h5space), intent(in), optional :: file_space
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        herr_t :: err
        hid_t :: mem_space_id
        hid_t :: file_space_id
        hid_t :: xfer_plist_id

        mem_space_id = h5s_all
        file_space_id = h5s_all
        xfer_plist_id = h5p_default
        if (present(mem_space)) mem_space_id = mem_space % id
        if (present(file_space)) file_space_id = file_space % id

        if (present(mpio_xfer)) then
! 'ifndef' inside 'if' to avoid warnings about unused dummy argument 'mpio_xfer'
#ifdef USE_MPI
            xfer_plist_id = h5pcreate(h5p_dataset_xfer)
            err = h5pset_dxpl_mpio(xfer_plist_id, mpio_xfer)
            assertm(err > -1, 'h5pset_dxpl_mpio failed')
#else
            continue
#endif
        end if

        err = h5dread(self % id, self % dtype, mem_space_id, file_space_id, xfer_plist_id, c_loc(buf))
        assertm(err > -1, 'h5dread failed')

        if (xfer_plist_id /= h5p_default) then
            err = h5pclose(xfer_plist_id)
            assertm(err > -1, 'h5pclose failed')
        end if
    end subroutine h5dataset_read

    subroutine h5dataset_write(self, buf, mdims, fdims, mem_space, file_space, mpio_xfer)
        class(h5dataset), intent(inout) :: self
        type(*), intent(in), target :: buf(..)
        integer, intent(in), optional :: mdims(:)
        integer, intent(in), optional :: fdims(:)
        type(h5space), intent(in), optional :: mem_space
        type(h5space), intent(in), optional :: file_space
        h5fd_mpio_xfer_t, intent(in), optional :: mpio_xfer

        herr_t :: err

        type(h5space) :: mdims_space
        type(h5space) :: fdims_space
        hid_t :: mem_space_id
        hid_t :: file_space_id

        hid_t :: xfer_plist_id

        mem_space_id = h5s_all
        file_space_id = h5s_all
        xfer_plist_id = h5p_default

        assertm(.not.(present(mdims).and.present(mem_space)), 'both, mdims and mem_space supplied.')
        assertm(.not.(present(fdims).and.present(file_space)), 'both, fdims and file_space supplied.')

        if (present(mem_space)) mem_space_id = mem_space % id
        if (present(file_space)) file_space_id = file_space % id
        if (present(mdims)) then
            call mdims_space % init(mdims)
            mem_space_id = mdims_space % id
        end if
        if (present(fdims)) then
            call fdims_space % init(fdims)
            file_space_id = fdims_space % id
        end if
        if (.not. (present(mdims) .or. present(mem_space) .or. present(fdims) .or. present(file_space))) then
            call mdims_space % init(find_h5shape(buf))
            mem_space_id = mdims_space % id
        end if

        if (present(mpio_xfer)) then
! 'ifdef' inside 'if' to avoid warnings about unused dummy argument 'mpio_xfer'
#ifdef USE_MPI
            xfer_plist_id = h5pcreate(h5p_dataset_xfer)
            err = h5pset_dxpl_mpio(xfer_plist_id, mpio_xfer)
            assertm(err > -1, 'h5pset_dxpl_mpio failed')
#else
            continue
#endif
        end if

! https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Write
        err = h5dwrite(self % id, self % dtype, mem_space_id, file_space_id, xfer_plist_id, c_loc(buf))
        assertm(err > -1, 'h5dwrite failed')

        if (xfer_plist_id /= h5p_default) then
            err = h5pclose(xfer_plist_id)
            assertm(err > -1, 'h5pclose failed')
        end if
    end subroutine h5dataset_write

    subroutine h5dataset_close(self)
        class(h5dataset), intent(inout) :: self
        herr_t :: err
!         integer :: dealloc_stat
        debug_print 'h5dataset_close', self % id, self % name
        if (h5iis_valid(self % dtype) > 0) err = h5tclose(self % dtype)
        assertm(err > -1, 'h5tclose failed')
!         deallocate(self % space, stat=dealloc_stat) ! this is pretty silly but i'm not going to carry extra variable jut to signify whether a pointer has been associated or allocated, 'allocated()' should simply be made to accept pointers....  scratch this, stat does not work with gfortran (8.2.1) so extra variable it is ...
        if (self % space_allocated) deallocate(self % space)
        self % space_allocated = .false.
        self % space => null()
        err = h5dclose(self % id)
        assertm(err > -1, 'h5dclose failed')
    end subroutine h5dataset_close

    subroutine h5dataset_clear(self)
        type(h5dataset), intent(inout) :: self
        debug_print 'h5dataset_clear', self % id, self % name
        if (h5iis_valid(self % id) <= 0) return
        call self % close()
        self % id = h5i_invalid_hid
    end subroutine h5dataset_clear


!%%%%%%%%%%%%%
    subroutine h5group_init(self, location, name)
! this does not work recursively because the way it is currently structured it will go into quadratic difficulty, maybe h5mkdir shall be removed and functionality brought in here to keep it linear.
        class(h5group), intent(inout) :: self
        class(h5location), intent(in) :: location
        character(len=*), intent(in) :: name

        herr_t :: err

        H5_CLASS_INIT

        assertm(len_trim(name) > 0, 'empty path supplied')
        self % name = trim(name)//c_null_char

        debug_print name
        err = h5lexists(location % id, self % name, h5p_default)
        assertm(err > -1, 'h5lexists failed')
        if (err > 0) then
            self % id = h5gopen2(location % id, self % name, h5p_default)
        else
!         self % id = h5gcreate2(location % id, trim(self % name)//c_null_char, lcpl_id, gcpl_id, gapl_id)
            self % id = h5gcreate2(location % id, self % name, h5p_default, h5p_default, h5p_default)
        end if
        debug_print 'h5group_init', self % id, self % name
    end subroutine h5group_init

    subroutine h5group_close(self)
        class(h5group), intent(inout) :: self
        herr_t :: err
        debug_print 'h5group_close', self % id, self % name
        err = h5gclose(self % id)
        assertm(err > -1, 'h5gclose failed')
    end subroutine h5group_close

    subroutine h5group_clear(self)
        type(h5group), intent(inout) :: self
        debug_print 'h5group_clear', self % id, self % name
        if (h5iis_valid(self % id) <= 0) return
        call self % close()
        self % id = h5i_invalid_hid
    end subroutine h5group_clear

    end module h5
