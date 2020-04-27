#include <cmath>
#include <cstdio>
#include <stdlib.h>
// #include <mkl_cblas.h>
#include <ctime>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <iostream>
#include <cassert>

using real  = double;
using cmplx = cuDoubleComplex;


std::ostream& operator<<(std::ostream& o, cudaDeviceProp& p) {
    o << p.name
      << ", total mem: " << (p.totalGlobalMem>>20) << " MiB"
      << ", multiproc count: " << (p.multiProcessorCount)
      << ", compute mode: " << p.computeMode;
    return o;
}


/**
 * Since this is ad-hoc style thing just to see if it works ok, all that device/context handle stuff that has to be passed from above will be saved in static vars.
 */
namespace cuconf {
    static bool blasinited {false};
    static cublasHandle_t ctx;

//     device allocated
    static cmplx* one;
    static cmplx* nul;
    static cmplx* img;

    constexpr int nstm {16};
    static cudaStream_t st0;
    static cudaStream_t* stm;

    constexpr auto n_ = CUBLAS_OP_N;
    constexpr auto t_ = CUBLAS_OP_T;
    constexpr auto c_ = CUBLAS_OP_C;

    constexpr auto lo_ = CUBLAS_FILL_MODE_LOWER;
    constexpr auto up_ = CUBLAS_FILL_MODE_UPPER;

    constexpr auto ls_ = CUBLAS_SIDE_LEFT;
    constexpr auto rs_ = CUBLAS_SIDE_RIGHT;
};


extern "C" void init_cublas() {
    using namespace cuconf;

    cublasCreate(&ctx);
    blasinited = true;
    int devid;
    cudaDeviceProp p;
    cudaGetDevice(&devid);
    cudaGetDeviceProperties(&p, devid);
    std::cout << p << std::endl;

    auto one_l = make_cuDoubleComplex(1.0, 0.0);
    auto nul_l = make_cuDoubleComplex(0.0, 0.0);
    auto img_l = make_cuDoubleComplex(0.0, 1.0);

    cudaMalloc(&one, sizeof(cmplx));
    cudaMalloc(&nul, sizeof(cmplx));
    cudaMalloc(&img, sizeof(cmplx));

    cudaMemcpy(one, &one_l, sizeof(cmplx), cudaMemcpyHostToDevice);
    cudaMemcpy(nul, &nul_l, sizeof(cmplx), cudaMemcpyHostToDevice);
    cudaMemcpy(img, &img_l, sizeof(cmplx), cudaMemcpyHostToDevice);

//     Supposedly faster access to the scalars (one and null in our case).
    cublasSetPointerMode(ctx, CUBLAS_POINTER_MODE_DEVICE);

// Streams are necessary when the matrixes are not large enough to utilise the device.
// Streams work independently but do not overlap with the default system stream. For
// this reason computing shall not be send to st0.
    stm = new cudaStream_t[nstm];
    cublasGetStream(ctx, &st0);
    for (int i = 0; i < nstm; i++) cudaStreamCreate(&stm[i]);

//     need this to flush stdout when called from fortran.
    cudaDeviceSynchronize();
}

extern "C" void stop_cublas() {
    using namespace cuconf;
    cublasDestroy(ctx);
    blasinited = false;

    for (int i = 0; i < nstm; i++) cudaStreamDestroy(stm[i]);
    delete [] stm;

    cudaFree(one);
    cudaFree(nul);
    cudaFree(img);
}


/**
 * m: mixed basis size
 * n: single basis size
 * k: number of states to reduce
 * the above would naturaly be size_t in c/c++ but here they are int for compatibility with the f code.
 *   _h suffix refers to host
 * w: screened potential
 * z: mixed product basis elements: <ψψ|M>
 * s: the zwz^* product to return
 * same_z: if 0 copy z from host memory otherwise reuse the one already in device memory
 *         This had to be integer instead of bool because the fortran 'logical' does not map to any C type.
 */
template<typename cmplx>
void zwz(const int m, const int n, const int k,
         const cmplx *w_h, const int ldw,
         const cmplx *z_h, const int ldz,
               cmplx *s_h, const int lds,
         const int same_z
        ) {

    using namespace cuconf;

//    constexpr bool clockon {false};
    constexpr bool clockon {true};
    cudaError_t err;

    clock_t t[8];
    auto timepoint = [&](const int i){cudaDeviceSynchronize(); t[i] = clock();};


    static cmplx *w {nullptr};
    static cmplx *z {nullptr};
    static cmplx *s {nullptr};
    static cmplx *c {nullptr};

    if (clockon) timepoint(0);

    if (!same_z) {
        if (w != nullptr) cudaFree(w);
        if (z != nullptr) cudaFree(z);
        if (s != nullptr) cudaFree(s);
        if (c != nullptr) cudaFree(c);

//         printf("cudaMalloc returned error code %d, line(%d)\n", err, __LINE__);
        err = cudaMalloc(&w, m*m     *sizeof(cmplx)); assert(err == cudaSuccess);
        err = cudaMalloc(&z, m*n*k   *sizeof(cmplx)); assert(err == cudaSuccess);
        err = cudaMalloc(&s, n*n*k   *sizeof(cmplx)); assert(err == cudaSuccess);
        err = cudaMalloc(&c, m*n*nstm*sizeof(cmplx)); assert(err == cudaSuccess);

        cublasSetMatrix(m, n*k, sizeof(cmplx), z_h, ldz, z, m);
    }

    cublasSetMatrix(m, m, sizeof(cmplx), w_h, ldw, w, m);

    if (clockon) timepoint(1);

    for (int i = 0; i < k; i++) {
        cublasSetStream(ctx, stm[i%nstm]);
        cublasZhemm(ctx, ls_, lo_, m, n, one, w, m, &z[m*n*i], m, nul, &c[m*n*(i%nstm)], m);
        cublasZgemm(ctx, c_, n_, n, n, m, one, &z[m*n*i], m, &c[m*n*(i%nstm)], m, nul, &s[n*n*i], n);
    }

    cublasSetStream(ctx, st0);

    if (clockon) timepoint(2);

    cublasGetMatrix(n, n*k, sizeof(cmplx), s, n, s_h, lds);

    if (clockon) timepoint(3);

    if (clockon) {
        double dt;
        dt = t[1] - t[0]; printf("  cuzwz prep: %12.8e\n", dt/CLOCKS_PER_SEC);
        dt = t[2] - t[1]; printf("  cuzwz calc: %12.8e\n", dt/CLOCKS_PER_SEC);
        dt = t[3] - t[2]; printf("  cuzwz copo: %12.8e\n", dt/CLOCKS_PER_SEC);
        std::cout << std::endl;
//         dt = t[4] - t[3]; printf("  cuzwz deal: %12.8e\n", dt/CLOCKS_PER_SEC);
    }
}

extern "C"
void fcuzwz(const int m, const int n, const int k,
         const double *w, const int ldw,
         const double *z, const int ldz,
               double *s, const int lds,
         const int same_z ) {
    zwz(m,n,k,
        reinterpret_cast<const cmplx*>(w), ldw,
        reinterpret_cast<const cmplx*>(z), ldz,
        reinterpret_cast<      cmplx*>(s), lds, same_z);
}


// int main() {
//     init_cublas();
//     stop_cublas();
// }
