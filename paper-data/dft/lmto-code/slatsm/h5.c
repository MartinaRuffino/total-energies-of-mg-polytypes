// h5 foo © 2018 Dimitar Pashov

#include <hdf5.h>
#include <assert.h>
#include <string.h>

#if (defined USE_MPI && !defined H5_HAVE_PARALLEL)
#error MPI support requires MPI enabled HDF5.
#endif

void h5_ftypes_init(size_t n, hid_t *types) {

    H5open();

    hid_t ta[] = {
          H5T_NATIVE_INT8
        , H5T_NATIVE_INT16
        , H5T_NATIVE_INT32
        , H5T_NATIVE_INT64
        , H5T_NATIVE_FLOAT
        , H5T_NATIVE_DOUBLE
        , H5T_STD_I8BE
        , H5T_STD_I8LE
        , H5T_STD_I16BE
        , H5T_STD_I16LE
        , H5T_STD_I32BE
        , H5T_STD_I32LE
        , H5T_STD_I64BE
        , H5T_STD_I64LE
        , H5T_IEEE_F32BE
        , H5T_IEEE_F32LE
        , H5T_IEEE_F64BE
        , H5T_IEEE_F64LE
        , H5P_ATTRIBUTE_CREATE
        , H5P_DATASET_ACCESS
        , H5P_DATASET_CREATE
        , H5P_DATASET_XFER
        , H5P_DATATYPE_ACCESS
        , H5P_DATATYPE_CREATE
        , H5P_FILE_ACCESS
        , H5P_FILE_CREATE
        , H5P_FILE_MOUNT
        , H5P_GROUP_ACCESS
        , H5P_GROUP_CREATE
        , H5P_LINK_ACCESS
        , H5P_LINK_CREATE
        , H5P_OBJECT_COPY
        , H5P_OBJECT_CREATE
        , H5P_STRING_CREATE
    };

    assert(sizeof(ta) == n*sizeof(hid_t));
    memcpy(types, ta, n*sizeof(hid_t));
}

herr_t H5Pset_fapl_mpio_foo(hid_t fapl_id, int comm, int info) {
#if (defined USE_MPI && defined H5_HAVE_PARALLEL)
    return H5Pset_fapl_mpio(fapl_id, MPI_Comm_f2c(comm), MPI_Info_f2c(info));
#else
    return (herr_t)0x1;
#endif
}

#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8)
htri_t H5Sis_regular_hyperslab( hid_t space_id ) {return (htri_t)0x0;}
herr_t H5Sget_regular_hyperslab( hid_t space_id, hsize_t start[], hsize_t stride[], hsize_t count[], hsize_t block[] ) {return (herr_t)0x1;}
#endif
