#define _GNU_SOURCE
#include <dlfcn.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>

struct cudadevmap_t {
   char remapped;
   int map[64], rmap[64];
   int total, usable;
   cudaError_t (*cudaSetDevice)(int);
   cudaError_t (*cudaGetDevice)(int*);
   cudaError_t (*cudaGetDeviceCount)(int*);
   cudaError_t (*cudaGetDeviceProperties)(struct cudaDeviceProp*, int);
};

struct cudadevmap_t cudadevmap = {.remapped = 0, .usable = 0};


int device_is_usable(int device) {
   int accepted = 0;
   struct cudaDeviceProp p;
   cudadevmap.cudaGetDeviceProperties(&p, device);
   accepted = ((p.major * 10 + p.minor > 20) &&
           (p.totalGlobalMem > (((size_t)1)<<32)) && /* 1<<32 == 2**32 == 4GiB */
           (p.multiProcessorCount >= 4 ));
   return accepted;
}

void remapCudaDevs() {

   char *dlerr = dlerror(); free(dlerr);

   cudadevmap.cudaSetDevice = (cudaError_t (*)(int)) dlsym(RTLD_NEXT, "cudaSetDevice");
   dlerr = dlerror(); if (dlerr != 0) fprintf(stderr, "cudaSetDevice SHIM dlsym err: %s\n", dlerr); free(dlerr);

   cudadevmap.cudaGetDevice = (cudaError_t (*)(int*)) dlsym(RTLD_NEXT, "cudaGetDevice");
   dlerr = dlerror(); if (dlerr != 0) fprintf(stderr, "cudaGetDevice SHIM dlsym err: %s\n", dlerr); free(dlerr);

   cudadevmap.cudaGetDeviceCount = (cudaError_t (*)(int*)) dlsym(RTLD_NEXT, "cudaGetDeviceCount");
   dlerr = dlerror(); if (dlerr != 0) fprintf(stderr, "cudaSetDevice SHIM dlsym err: %s\n", dlerr); free(dlerr);

   cudadevmap.cudaGetDeviceProperties = (cudaError_t (*)(struct cudaDeviceProp*,int)) dlsym(RTLD_NEXT, "cudaGetDeviceProperties");
   dlerr = dlerror(); if (dlerr != 0) fprintf(stderr, "cudaSetDevice SHIM dlsym err: %s\n", dlerr); free(dlerr);

   cudadevmap.cudaGetDeviceCount (&cudadevmap.total);

   if (cudadevmap.total > 16) {
      fprintf(stderr, "cudaSetDevice SHIM increase cudadevmap_t.map size to %d\n", cudadevmap.total);
      exit(-1);
   }

   cudadevmap.usable = 0;
   for (int i = 0; i < cudadevmap.total; ++i) {
      if (device_is_usable(i) == 1) {
         cudadevmap.map[cudadevmap.usable] = i;
         cudadevmap.rmap[i] = cudadevmap.usable;
         ++cudadevmap.usable;
      }
   }

   cudadevmap.remapped = 1;
}

cudaError_t cudaSetDevice(int device) {
   if (!cudadevmap.remapped) remapCudaDevs();
   return cudadevmap.cudaSetDevice(cudadevmap.map[device]);
}

cudaError_t cudaGetDevice(int *device) {
   if (!cudadevmap.remapped) remapCudaDevs();
   cudaError_t e = cudadevmap.cudaGetDevice(device);
   *device = cudadevmap.rmap[*device];
   return e;
}


cudaError_t cudaGetDeviceCount(int* c) {
   if (!cudadevmap.remapped) remapCudaDevs();
   *c = cudadevmap.usable;
   return cudaSuccess;
}

cudaError_t cudaGetDeviceProperties(struct cudaDeviceProp *p, int device) {
   if (!cudadevmap.remapped) remapCudaDevs();
   return cudadevmap.cudaGetDeviceProperties(p, cudadevmap.map[device]);
}


