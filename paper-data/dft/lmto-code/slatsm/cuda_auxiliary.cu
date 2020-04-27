#include <cuda.h>
#include <cuda_runtime.h>

// extern "C"
// int device_is_usable(int devid) {
//    int accepted = 0;
//    cudaDeviceProp p;
//    cudaGetDeviceProperties(&p, devid);
//    accepted = ((p.major * 10 + p.minor > 20) &&
//            (p.totalGlobalMem > (((size_t)1)<<32)) && /* 1<<32 == 2**32 == 4GiB */
//            (p.multiProcessorCount >= 4 ));
//    return accepted;
// }

extern "C"
void fcuda_get_local_dev_props( int *devid, char *name, int *major, int *minor,
                                int *totalGlobalMem, int *multiProcessorCount, int *computeMode) {
   cudaGetDevice(devid);
   cudaDeviceProp p;
   cudaGetDeviceProperties(&p, *devid);
   strncpy(name, p.name, 16);
   *major = p.major;
   *minor = p.minor;
   *totalGlobalMem = (p.totalGlobalMem>>20);
   *multiProcessorCount = p.multiProcessorCount;
   *computeMode = p.computeMode;
}