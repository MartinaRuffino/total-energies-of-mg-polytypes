#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <mkl_cblas.h>
#include <time.h>
#include <cublas.h>
// #include <iostream>

// nvcc -O3 -arch=sm_35

#define BLOCK_SIZE 16

#ifdef MPICUD
#include <mpi.h>
extern "C"
void pcudgemm(char transa, char transb, int m, int n, int k, double alpha,
            double *a, int lda, double *b, int ldb, double beta, double *c, int ldc,
            MPI_Comm comm);
#endif

// template <int blocksize>
__global__ void m_plus_i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("b: %03d t: %04d i: %06d\n", blockIdx.x, threadIdx.x, t );
   if (t < n) m[t*n+t] += 1.0;
}

// template <int blocksize>
__global__ void m_mins_i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("b: %03d t: %04d i: %06d\n", blockIdx.x, threadIdx.x, t );
   if (t < n) m[t*n+t] -= 1.0;
}

// template <int blocksize>
__global__ void m_plus_35i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("b: %03d t: %04d i: %06d\n", blockIdx.x, threadIdx.x, t );
   if (t < n) m[t*n+t] += 35.0;
}


// template <int blocksize>
__global__ void m_mins_35i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("b: %03d t: %04d i: %06d\n", blockIdx.x, threadIdx.x, t );
   if (t < n) m[t*n+t] -= 35.0;
}

// template <int blocksize>
__global__ void m_plus_2i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
   if (t < n) m[t*n+t] += 2.0;
}

// template <int blocksize>
__global__ void m_plus_15i(double *m, int n) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
   if (t < n) m[t*n+t] += 15.0;
}


// template <int blocksize>
__global__ void m_plus_ai(double *m, int n, double a) {
   int t = blockIdx.x * blockDim.x + threadIdx.x;
   if (t < n) m[t*n+t] += a;
}

// template <int blocksize>
__global__ void i_mins_m(double *m, int n) {
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   if (x < n && y < n) {
      m[y*n + x] = - m[y*n + x];
      if (x == y) m[y*n + x] += 1.0;
   }
}

// template <int blocksize>
__global__ void i16_mins_5m(double *m, int n) {
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   if (x < n && y < n) {
      m[y*n + x] = -5.0*m[y*n + x];
      if (x == y) m[y*n + x] += 16.0;
   }
}




// template <int blocksize>
__global__ void i3m_mins_7i(double *m, int n) {
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   if (x < n && y < n) {
      m[y*n + x] *=  3.0 ;
      if (x == y) m[y*n + x] -= 7.0;
   }
}



extern "C" void r8_tofile(double *a, int m, int n, char *fln );


// void r8_fromfile(double *a, int m, int n, char *fln ) {
//
//    FILE * fl = fopen(fln, "w");
//    fprintf(fl, "%5d %5d\n", m, n);
//    for (int i = 0; i < m; ++i) {
//       for (int j = 0; j < n; ++j) {
//          fprintf(fl, " %5.2f", ha[i*n+j]);
//       }
//       fprintf(fl, "\n");
//    }
//    fclose(fl);
// }


// void init_cuda_stuff() {
//    int devID = 0;
//    cudaError_t error;
//    cudaDeviceProp deviceProp;
//    error = cudaGetDevice(&devID); if (error != cudaSuccess) printf("cudaGetDevice returned error code %d, line(%d)\n", error, __LINE__);
//
//    error = cudaGetDeviceProperties(&deviceProp, devID);
//    if (error == cudaSuccess) printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);
//    else printf("cudaGetDeviceProperties returned error code %d, line(%d)\n", error, __LINE__);
//    if ( CUBLAS_STATUS_SUCCESS != cublasInit()) printf("cublas init err\n");
// }

void handlerr() {
   cublasStatus st = cublasGetError();
   if ( CUBLAS_STATUS_SUCCESS          != st ) printf("dev rt err          : %d\n",st);
   if ( CUBLAS_STATUS_NOT_INITIALIZED  == st ) printf("CUBLAS_STATUS_NOT_INITIALIZED  : %d\n",st);
   if ( CUBLAS_STATUS_ALLOC_FAILED     == st ) printf("CUBLAS_STATUS_ALLOC_FAILED     : %d\n",st);
   if ( CUBLAS_STATUS_ARCH_MISMATCH    == st ) printf("CUBLAS_STATUS_ARCH_MISMATCH    : %d\n",st);
   if ( CUBLAS_STATUS_EXECUTION_FAILED == st ) printf("CUBLAS_STATUS_EXECUTION_FAILED : %d\n",st);
   if ( CUBLAS_STATUS_INVALID_VALUE    == st ) printf("CUBLAS_STATUS_INVALID_VALUE    : %d\n",st);
   if ( CUBLAS_STATUS_MAPPING_ERROR    == st ) printf("CUBLAS_STATUS_MAPPING_ERROR    : %d\n",st);
   if ( CUBLAS_STATUS_INTERNAL_ERROR   == st ) printf("CUBLAS_STATUS_INTERNAL_ERROR   : %d\n",st);

}


void r8_tofile(double *a, int m, int n, char *fln ) {

   FILE * fl = fopen(fln, "w");
   fprintf(fl, "%5d %5d\n", m, n);
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         fprintf(fl, " %5.2f", a[i*n+j]);
      }
      fprintf(fl, "\n");
   }
   fclose(fl);
}


void leave() {
   if ( CUBLAS_STATUS_SUCCESS != cublasShutdown()) printf("cublas down err\n");
   cudaDeviceReset();
   exit(0);
}

extern "C" void matsign_itr_cu(int d, double *x, int maxit, double efilt) {
//       X must be rescaled by the conditioning number    X = A/||A|| (largest abs eigenvalue)
//       approximation to the sign function: sign(X) = X(X^2)^(-1/2)
//       Xn+1 = 1/2 Xn(3I-Xn^2)
//       Xinf -> sign(A)

//    clock_t t[6];
//    double td;
   bool last;

//    t[0] = clock();

   cudaError_t error;
//    init_cuda_stuff();

//    cudaDeviceSynchronize();


//    t[1] = clock();

   double ix2fnrm, ex2fnrm;
//       bool last;
   double *p, *r;
//       int sparsity[2];

   int const blocksize = BLOCK_SIZE*BLOCK_SIZE;

   double esq = sqrt(efilt);

   int d2 = d*d;

//       double *t  = malloc(2*d2*sizeof(double));
//       double *x2 = malloc(d2*sizeof(double));
   double  *t0, *t1, *x2;
   size_t basemem = d2*sizeof(double);
   error = cudaMalloc((void **) &t0, basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);
   error = cudaMalloc((void **) &t1, basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);
   error = cudaMalloc((void **) &x2, basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);

   int it = 0;
   int k;
   p = x;

//    int blocks = d/blocksize + ((d%blocksize > 0)?1:0);
   int blocks = (d + blocksize - 1) / blocksize;

   int ip, np;
   ip = 0;
   np = 1;
#ifdef MPICUD
   MPI_Comm_rank(MPI_COMM_WORLD, &ip);
   MPI_Comm_size(MPI_COMM_WORLD, &np);
#endif
//    np = 1;

   for (;;) {
      k = it % 2;

      if (np > 1) {
#ifdef MPICUD
         pcudgemm('t','n',d,d,d,-1.0,p,d,p,d,0.0,x2,d,MPI_COMM_WORLD); handlerr();
#endif
      } else {
         cublasDgemm('t','n',d,d,d,-1.0,p,d,p,d,0.0,x2,d); handlerr();
      }
//       cublasDsymm('l','u',d,d,-1.0,p,d,p,d,0.0,x2,d); handlerr();
//          where ( -efilt < x2 .and. x2 < efilt ) x2 = 0.0_dp
//         sparsity(1) = count(-efilt < x2 .and. x2 < efilt)
//          cudaDeviceSynchronize();

      ex2fnrm = esq * cublasDdot(d2, x2, 1, x2, 1); handlerr();
//          cudaDeviceSynchronize();
//          printf("ex2fnrm esq: %18.12f %18.12f\n", ex2fnrm, esq);
//          leave();
//          for (int i = 0; i < d; ++i) x2[i*d+i] = x2[i*d+i] + 1.0;
//       m_plus_i<blocksize><<<blocks, blocksize>>>(x2, d);
      m_plus_i<<<blocks, blocksize>>>(x2, d);
//          cudaDeviceSynchronize();

      ix2fnrm = cublasDdot(d2, x2, 1, x2, 1); handlerr();
      last = ((ix2fnrm < ex2fnrm) || (it >= maxit)) && (it != 0);
//          cudaDeviceSynchronize();

      if (ip==0) printf("p%d: %s%3d%s %18.12f %18.12f%s%s", ip,
               "sign itr:", it," ||I-X^2||F, efilt^0.5*||X^2||F:", ix2fnrm, ex2fnrm," last:", last?"T":"F");

//          for (int i = 0; i < d; ++i) x2[i*d+i] = x2[i*d+i] + 2.0;
//       m_plus_2i<blocksize><<<blocks, blocksize>>>(x2, d);
      m_plus_2i<<<blocks, blocksize>>>(x2, d);
//          cudaDeviceSynchronize();
      r = (k==0)?t0:t1; //t + k*d2;
      if (last) r = x;

      if (np>1) {
#ifdef MPICUD
         pcudgemm('t', 'n',d,d,d,0.5,p,d,x2,d,0.0,r,d,MPI_COMM_WORLD); handlerr();
#endif
      } else {
         cublasDgemm('t', 'n',d,d,d,0.5,p,d,x2,d,0.0,r,d); handlerr();
      }
//       cublasDsymm('l','u',d,d,0.5,p,d,x2,d,0.0,r,d); handlerr();
//          cudaDeviceSynchronize();
//          where ( -efilt < r .and. r < efilt ) x = 0.0_dp
//          sparsity(2) = count( -efilt < r .and. r < efilt)
//          if (aproc==0) write (*,'(" sparsity:", 2(x,f6.3), " %", 2(x,i0))') 100*real(sparsity,dp)/real(d2,dp), sparsity

      if (ip==0) printf("\n");

      if (last) break;

//       p = (k==0)?t0:t1; //t + k*d2;
      p = r;

      ++it;
   }

   cudaFree(t0);
   cudaFree(t1);
   cudaFree(x2);
//    cudaDeviceSynchronize();
//    t[4] = clock();
//    cudaDeviceSynchronize();
//    t[5] = clock();
//
//
//    cudaDeviceSynchronize();
//
//    t[2] = clock();
//
//    cudaDeviceSynchronize();
//
//    t[3] = clock();
//
//    td = t[1] - t[0]; printf("start time: %12.8fs\n", td/CLOCKS_PER_SEC);
//    td = t[2] - t[1]; printf("polynomial time: %12.8fs\n", td/CLOCKS_PER_SEC);
//    td = t[3] - t[2]; printf("stop time: %12.8fs\n", td/CLOCKS_PER_SEC);
//    td = t[5] - t[4]; printf("copy out time: %12.8fs\n", td/CLOCKS_PER_SEC);
}



extern "C" void matsign_itr4t_cu(int d, double *x, int maxit, double efilt) {
//       X must be rescaled by the conditioning number    X = A/||A|| (largest abs eigenvalue)
//       approximation to the sign function: sign(X) = X(X^2)^(-1/2)
//       Xn+1 = 1/8 Xn(15I + (3Xn^2-10)Xn^2)
//       Xinf -> sign(A)

   bool last;

   cudaError_t error;

   double ix2fnrm, ex2fnrm;
   double *p, *r, *t[4];



   double esq = sqrt(efilt);

   int d2 = d*d;

   size_t basemem = d2*sizeof(double);

   for (int i = 0; i < 4; ++i) {
      error = cudaMalloc((void **) &t[i], basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);
   }

   int const blocksize = BLOCK_SIZE*BLOCK_SIZE;
//    int blocks = d/blocksize + ((d%blocksize > 0)?1:0);
   int blocks = (d + blocksize - 1) / blocksize;

//    int blocks2d = d/BLOCK_SIZE + ((d%BLOCK_SIZE > 0)?1:0);
   int blocks2d = (d + BLOCK_SIZE - 1) / BLOCK_SIZE;
   dim3 bgrid(blocks2d, blocks2d);
   dim3 tgrid(BLOCK_SIZE, BLOCK_SIZE);


   int it = 0;
   int k0 = 0, k1 = 1, k2 = 2, k3 = 3;
   p = x;

   for (;;) {

      cublasDgemm('t','n',d,d,d,1.0,p,d,p,d,0.0,t[k0],d);           // x - > x^2 : t[0]
//          where ( -efilt < x2 .and. x2 < efilt ) x2 = 0.0_dp
//         sparsity(1) = count(-efilt < x2 .and. x2 < efilt)
//          cudaDeviceSynchronize();

      ex2fnrm = esq * cublasDdot(d2, t[k0], 1, t[k0], 1); handlerr();   // e_f^0.5*||x^2||

                                                                  // x^2 -> x^2 : t[1]
      error = cudaMemcpy(t[k1], t[k0], basemem, cudaMemcpyDeviceToDevice); if (error != cudaSuccess) printf("cudaMemcpy (d,d) returned error code %d, line(%d)\n", error, __LINE__);

//       m_mins_i<blocksize><<<blocks, blocksize>>>(t[k1], d);
      m_mins_i<<<blocks, blocksize>>>(t[k1], d);                    // x^2 -> x^2 - I : t[1]

      ix2fnrm = cublasDdot(d2, t[k1], 1, t[k1], 1); handlerr();      // ||x^2 - I||

      last = ((ix2fnrm < ex2fnrm) || (it >= maxit)) && (it != 0);

      printf("%s%3d%s %18.12f %18.12f%s%s",
         "sign itr:", it," ||X^2-I||F, efilt^0.5*||X^2||F:", ix2fnrm, ex2fnrm," last:", last?"T":"F");

//       i3m_mins_i<BLOCK_SIZE><<<bgrid, tgrid>>>(x, d);
      i3m_mins_7i<<<bgrid, tgrid>>>(t[k1], d);                            // 3*(x^2 - I) - 7*I : t[1]

      cublasDgemm('t','n',d,d,d,1.0,t[k0],d,t[k1],d,0.0,t[k2],d); handlerr(); // x^2*(3*x^2-10) : t[2]

//       m_plus_15i<blocksize><<<blocks, blocksize>>>(t[k2], d); handlerr();
      m_plus_15i<<<blocks, blocksize>>>(t[k2], d); handlerr();               // x^2*(3*x^2-10) + 15*I : t[2]

      r = (last) ? x : (it == 0 ? t[k3] : t[k0]);

      cublasDgemm('t','n',d,d,d,0.125,p,d,t[k2],d,0.0,r,d); handlerr();       // x*(x^2*(3*x^2-10) + 15*I) : t[3]


      printf("\n");

      if (last) break;

      p = r;

      k0 = it       % 4;
      k1 = (it + 1) % 4;
      k2 = (it + 2) % 4;
      k3 = (it + 3) % 4;

      ++it;
   }

   for (int i = 0; i < 4; ++i) cudaFree(t[i]);

//    leave();
}






extern "C" void matsign_itr6t_cu(int d, double *x, int maxit, double efilt) {
//       X must be rescaled by the conditioning number    X = A/||A|| (largest abs eigenvalue)
//       approximation to the sign function: sign(X) = X(X^2)^(-1/2)
//       Xn+1 = 1/8 Xn(15I + (3Xn^2-10)Xn^2)
//       Xinf -> sign(A)

   bool last;

   cudaError_t error;

   double ix2fnrm, ex2fnrm;
   double *p, *r, *t[4];



   double esq = sqrt(efilt);

   int d2 = d*d;

   size_t basemem = d2*sizeof(double);

   for (int i = 0; i < 4; ++i) {
      error = cudaMalloc((void **) &t[i], basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);
   }

   int const blocksize = BLOCK_SIZE*BLOCK_SIZE;
//    int blocks = d/blocksize + ((d%blocksize > 0)?1:0);
   int blocks = (d + blocksize - 1) / blocksize;

//    int blocks2d = d/BLOCK_SIZE + ((d%BLOCK_SIZE > 0)?1:0);
   int blocks2d = (d + BLOCK_SIZE - 1) / BLOCK_SIZE;
   dim3 bgrid(blocks2d, blocks2d);
   dim3 tgrid(BLOCK_SIZE, BLOCK_SIZE);


   int it = 0;
   int k = 3;
   p = x;

//    printf("0: %p ; 3: %p\n", t[0], t[3]);
   for (;;) {

      cublasDgemm('t','n',d,d,d,1.0,p,d,p,d,0.0,t[1],d);           // x - > x^2 : t[0]
//          where ( -efilt < x2 .and. x2 < efilt ) x2 = 0.0_dp
//         sparsity(1) = count(-efilt < x2 .and. x2 < efilt)
//          cudaDeviceSynchronize();
//       printf("begin: %p * p -> 1\n", p);

      ex2fnrm = esq * cublasDdot(d2, t[1], 1, t[1], 1); handlerr();   // e_f^0.5*||x^2||

                                                                  // x^2 -> x^2 : t[1]
      error = cudaMemcpy(t[2], t[1], basemem, cudaMemcpyDeviceToDevice); if (error != cudaSuccess) printf("cudaMemcpy (d,d) returned error code %d, line(%d)\n", error, __LINE__);
//       printf("1 -> 2\n");

//       m_mins_i<blocksize><<<blocks, blocksize>>>(t[k1], d);
      m_mins_i<<<blocks, blocksize>>>(t[2], d);                    // x^2 -> x^2 - I : t[2]

      ix2fnrm = cublasDdot(d2, t[2], 1, t[2], 1); handlerr();      // ||x^2 - I||

      last = ((ix2fnrm < ex2fnrm) || (it >= maxit)) && (it != 0);

      printf("%s%3d%s %18.12f %18.12f%s%s",
         "sign itr:", it," ||X^2-I||F, efilt^0.5*||X^2||F:", ix2fnrm, ex2fnrm," last:", last?"T":"F");

//       i3m_mins_i<BLOCK_SIZE><<<bgrid, tgrid>>>(x, d);
      i16_mins_5m<<<bgrid, tgrid>>>(t[2], d);                            // -5*(x^2 - I) + 26*I : t[2]

      cublasDgemm('t','n',d,d,d,1.0,t[1],d,t[2],d,0.0,t[k],d); handlerr(); // x^2*(-5*x^2 + 21) : t[k] : t[3] ili t[0]
//       printf("1 * 2 -> %d\n", k);

//       m_plus_15i<blocksize><<<blocks, blocksize>>>(t[k2], d); handlerr();
      m_mins_35i<<<blocks, blocksize>>>(t[k], d); handlerr();               // x^2*(-5*x^2 + 21) - 35 : t[k]

      cublasDgemm('t','n',d,d,d,1.0,t[1],d,t[k],d,0.0,t[2],d); handlerr(); // x^2*(x^2*(-5*x^2 + 21) - 35) : t[2]
//       printf("1 * %d -> 2\n", k);

      m_plus_35i<<<blocks, blocksize>>>(t[2], d); handlerr();               // 35 + x^2*(x^2*(-5*x^2 + 21) - 35) : t[2]

      r = (last) ? x : ( it == 0 ? t[0] : t[k]);

      cublasDgemm('t','n',d,d,d,0.0625,p,d,t[2],d,0.0,r,d); handlerr();       // x*(35 + x^2*(x^2*(-5*x^2 + 21) - 35)) : t[3] ili t[0] ili x

//       printf("%p * 2 -> %p\n", p, r);

      printf("\n");

      if (last) break;

      p = r;

      k = 3*((it + 1) % 2);
/*
      if ((it%2) != 0) {
         int k = k0;
         k0 = k3;
         k3 = k;
      }*/

      ++it;
   }

   for (int i = 0; i < 4; ++i) cudaFree(t[i]);

//    leave();
}













extern "C" void h2rho_host(int d, double *hx, double ef, int maxit, double efilt) {
   cudaError_t error;

   double *x;
   size_t basemem = d*d*sizeof(double);
   error = cudaMalloc((void **) &x , basemem); if (error != cudaSuccess) printf("cudaMalloc returned error code %d, line(%d)\n", error, __LINE__);
   error = cudaMemcpy(x, hx, basemem, cudaMemcpyHostToDevice); if (error != cudaSuccess) printf("cudaMemcpy (d_A,h_A) returned error code %d, line(%d)\n", error, __LINE__);

   int const blocksize = BLOCK_SIZE*BLOCK_SIZE;
//    int blocks = d/blocksize + ((d%blocksize > 0)?1:0);
   int blocks = (d + blocksize - 1) / blocksize;

//    m_plus_ai<blocksize><<<blocks, blocksize>>>(x, d, -ef);
   m_plus_ai<<<blocks, blocksize>>>(x, d, -ef);

   matsign_itr_cu(d, x, maxit, efilt);
//    matsign_itr4t_cu(d, x, maxit, efilt);
//    matsign_itr6t_cu(d, x, maxit, efilt);


//    blocks = d/BLOCK_SIZE + ((d%BLOCK_SIZE > 0)?1:0);
   blocks = (d + BLOCK_SIZE - 1)/BLOCK_SIZE;
   dim3 bgrid(blocks, blocks);
   dim3 tgrid(BLOCK_SIZE, BLOCK_SIZE);

//    i_mins_m<BLOCK_SIZE><<<bgrid, tgrid>>>(x, d);
   i_mins_m<<<bgrid, tgrid>>>(x, d);

   error = cudaMemcpy(hx, x, basemem, cudaMemcpyDeviceToHost); if (error != cudaSuccess) printf("cudaMemcpy (h_C,d_C) returned error code %d, line(%d)\n", error, __LINE__);
   cudaFree(x);
}

