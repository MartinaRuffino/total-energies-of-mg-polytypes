void cudaGetDevice (int *devID) {}
void cudaSetDevice (int devID) {}
void cudaGetDeviceCount ( int *total_ng) { *total_ng = 0;}

// int device_is_usable(int devid) {return 0;}
void fcuda_get_local_dev_props( int *devid, char *name, int *major,
                                int *minor, int *totalGlobalMem,
                                int *multiProcessorCount, int *computeMode) {}

