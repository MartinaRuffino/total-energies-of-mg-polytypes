#define CUDA

#ifdef CUDA
#include <cuda.h>
#endif
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>


#define LMAX  3
#define LMMAX 16
#define MAXBLOCKSZ 256


// void leave() {
//    if ( CUBLAS_STATUS_SUCCESS != cublasShutdown()) printf("cublas down err\n");
//    cudaDeviceReset();
//    exit(0);
// }



/*
__device__
void besslr(double y, int loka, int lmin, int lmax, double *fi, double *gi) {
//       double precision y,fi(lmin:lmax),gi(lmin:lmax), fac2l(-nlmax:nlmax*2+3)
      int i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2;
      int const nlmax=20;
      double dt,dt2,exppr,my,srmy,g1,t, dum[nlmax*4+2], fac2l[nlmax*3+4];
      bool const lhank = true;
      double const tol=1e-15;

      lmx = lmax > 2 ? lmax : 2;
      if (lmin > 0) printf(" BESSL : lmin gt 0\n");
      if (lmx > nlmax+nlmax) printf(" BESSL : lmax gt nlmax*2, lmax=%d",lmx);

//C --- A table of fac2l(l)=(2l-1)!!
//c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      fac2l(0) = 1d0
      do  l = 1, lmx+1
        fac2l(l) = fac2l(l-1)*(l+l-1)
      enddo
      do  l = -1, lmin,-1
        fac2l(l) = fac2l(l+1)/(l+l+1)
      enddo

// C --- Case akap=0 ---
      if (y == 0) then
        do  l = lmin, lmax
          fi(l) = 1/fac2l(l+1)
          gi(l) = fac2l(l)
        enddo
        goto 100
      endif
      my = -y

// C --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
      tlp1 = lmx+lmx+1
      dt = 1d0
      t = 1d0
      i = 0
      do  k = 1, 1000
        if (dabs(dt) < tol) goto 21
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
      enddo
      call rx('BESSLR: series not convergent')
   21 continue
      dum(1) = t/fac2l(lmx+1)

// C --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
      tlp1 =  tlp1-2
      dt = 1d0
      t = 1d0
      i = 0
      do  k = 1, 1000
        if (dabs(dt) < tol) goto 31
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
      enddo
      call rx('BESSLR: series not convergent')
   31 continue
      dum(2) = t/fac2l(lmx)

// C --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
      ll1 = lmx + lmx + 1
      ll2 = ll1 + 1
      nf = ll1
      do  k = 3, ll2
        nf = nf-2
        dum(k) = nf*dum(k-1) - y*dum(k-2)
      enddo

// C --- Get fi and gi from dum ---
      lmxp1 = lmx+1
      lmxp2 = lmx+2
      isn = (-1)**lmin
      do  k = lmin, lmax
        j1 = lmxp1-k
        j2 = lmxp2+k
        fi(k) = dum(j1)
// c   ... n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi(k) = dum(j2)*isn
        isn = -isn
      enddo

// C --- For E<0, use Hankel functions rather than Neumann functions ---
      if (lhank .and. y < 0d0) then
        srmy = dsqrt(-y)
        gi(0) = 1d0
        g1 = 1d0+srmy
        if (lmax >= 1) gi(1) = g1
        if (lmax >= 2) then
          tlp1 = 1
          do  l = 2, lmax
            tlp1 = tlp1+2
            gi(l) = tlp1*gi(l-1) - y*gi(l-2)
          enddo
        endif
        if (lmin <= -1) then
          gi(-1) = (gi(0) - g1)/y
          tlp1 = 1
          if (lmin <= -2) then
            do  l = -2, lmin,-1
              tlp1  = tlp1-2
              gi(l) = ((l+l+3)*gi(l+1) - gi(l+2))/y
            enddo
          endif
        endif
        exppr = 1d0/dexp(srmy)
        do  l = lmin, lmax
          gi(l) = gi(l)*exppr
        enddo
      endif

// C --- Scaling to Andersen's 2nd generation LMTO conventions ---
  100 continue
      if (loka == 1) then
        do  l = lmin, lmax
        fi(l) = fi(l)*fac2l(l)*0.5d0
        gi(l) = gi(l)/fac2l(l)
        enddo
      endif
}

*/


 


__device__
double sylm_cu(double r[3], double *yl, int lmx) {
//- Generate unnormalized spherical harmonic polynomials
// ----------------------------------------------------------------
//i Inputs
//i   r     :vector with 3 components
//i   lmax  :maximum l for which ylm is calculated
//o Outputs:
//o   ylm   :unnormalized spherical harmonic polynomials
//o   r2s   :length of dr**2
//r Remarks:
//r   polar axis along 001 axis. (adapted from ASW programs)
//r   use together with sylmnc:
//r   The true spherical harmonic is:  Y_L = ylm(lm,r) / r^l
//r   The first 9 values are for ylm/cy:
//r     l=0:                 1
//r     l=1:        y        z       x
//r     l=2:  6xy  3yz  3zz/2-rr/2  3xz  3(xx-yy)
//r   Factors cy are (F = 4*pi)
//r     l=0:              sqrt(1/4F)
//r     l=1:              sqrt(3/F)   sqrt(3/F)  sqrt(3/F)
//r     l=2:  sqrt(5/12F) sqrt(5/3F)  sqrt(5/F)  sqrt(5/3F) sqrt(5/12F)
//u Updates
//u   18 Sep 06 Set tolerance to smaller number for gen. Ylm
// ----------------------------------------------------------------
   int lav,lp1,lm1,n,nt;
   double r2,st,z2, r2s;
   double c[LMMAX], s[LMMAX], p[LMMAX][LMMAX];
   double &x = c[1];
   double &y = s[1];
   double &z = p[1][0];

   c[0] = 1.0;
   s[0] = 0.0;
   p[0][0] = 1.0;
   p[1][1] = 1.0;

   n = lmx+1;
   n *= n;
   yl[0] = 1.0;
   x = r[0];
   y = r[1];
   z = r[2];
   st = x*x + y*y;
   z2 = z*z;
   r2 = st+z2;
   r2s = r2;
   if (n < 2) return r2s;
   if (r2 <= 1e-28) {
      for (int i = 1; i < n; ++i )  yl[i] = 0.0;
      return r2s;
   }

   yl[1] = y;
   yl[2] = z;
   yl[3] = x;
   nt = 1;
   for ( int l = 2; l <= lmx; ++l) {
      lp1 = l+1;
      lm1 = l-1;
      lav = l*lp1;
      p[l][0] = ((2*l-1)*z*p[l-1][0] - lm1*r2*p[lm1-1][0]) / l;
      yl[lav] = p[l][0];
      nt = nt+2;
      p[l][l] = p[lm1][lm1]*nt;
      c[l] = x*c[lm1] - y*s[lm1];
      s[l] = x*s[lm1] + y*c[lm1];
      yl[lav+l] = p[l][l]*c[l];
      yl[lav-l] = p[l][l]*s[l];
      if (st > z2) {
         for (int m = 1; m <= lm1; ++m) {
            p[l][m] = ((lm1+m)*r2*p[lm1][m-1]-(lp1-m)*z*p[l][m-1])/st;
            yl[lav+m] = p[l][m]*c[m];
            yl[lav-m] = p[l][m]*s[m];
         }
      } else {
         for ( int lmm = 1; lmm <= lm1; ++lmm ) {
            int m = l-lmm;
            p[l][m] = ((l+m)*r2*p[lm1][m]-st*p[l][m+1])/(z*(l-m));
            yl[lav+m] = p[l][m]*c[m];
            yl[lav-m] = p[l][m]*s[m];
         }
      }
   }
   return r2s;
}

// // __global__
// void sylm_cu_glob(double r[3], double *yl, int lmx, double *r2s) {
//    *r2s = sylm_cu(r, yl, lmx);
// }



// #ifdef MAIN
// #ifdef CUDA
// void init_cuda_stuff() {
//    int devID = 0;
//    cudaError_t err;
//    cudaDeviceProp deviceProp;
//    err = cudaGetDevice(&devID); if (err != cudaSuccess) printf("cudaGetDevice: err: %d, line(%d)\n", err, __LINE__);
//
//    err = cudaGetDeviceProperties(&deviceProp, devID);
//    if (err == cudaSuccess) printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);
//    else printf("cudaGetDeviceProperties: err: %d, line(%d)\n", err, __LINE__);
// //    if ( CUBLAS_STATUS_SUCCESS != cublasInit()) printf("cublas init err\n");
// }
// #endif
//

// int main(void) {
//
//
// //    std::cout << "cuda" << std::endl;
//
//    double r[3] = {-0.079935200000 ,  -27.986444848657,     0.569759600000};
//
//    double ylref[16] = {
//       1.000000000000E+00,
//       -2.798644484866E+01,
//       5.697596000000E-01,
//       -7.993520000000E-02,
//       1.342261239760E+01,
//       -4.783663686718E+01,
//       -3.912991164498E+02,
//       -1.366315427338E-01,
//       -2.349704116892E+03,
//       3.287939586943E+05,
//       3.823831135305E+01,
//       3.282595804246E+04,
//       -6.692092017682E+02,
//       9.375787226657E+01,
//       -6.693832388794E+03,
//       2.817376350574E+03,
//    };
//
//    double r2sref = 783.572110904925;
//    double yl[16], r2s;
//
//    clock_t t[6];
//    double td;
//
// #ifndef CUDA
//    t[0] = clock();
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    sylm_h(r,yl,3,r2s);
//    t[1] = clock();
//    td = t[1] - t[0]; printf("sylm_h time: %12.8fs\n", td/CLOCKS_PER_SEC);
// #else
//    t[0] = clock();
//    init_cuda_stuff();
//    t[1] = clock();
//    sylm_hcu(r,yl,3,r2s);
//    t[2] = clock();
//    td = t[1] - t[0]; printf("init time: %12.8fs\n", td/CLOCKS_PER_SEC);
//    td = t[2] - t[1]; printf("sylm_hcu time: %12.8fs\n", td/CLOCKS_PER_SEC);
// #endif
//
//    printf("yl:\n");
//    for (int i=0; i<16; ++i) {
//       printf("%22.12f %22.12f %22.12f \n", yl[i], ylref[i], yl[i] - ylref[i]);
//    }
//
//    printf("r2s:\n %22.12f %22.12f %22.12f \n", r2s, r2sref, r2s - r2sref);
//
//
//
// //    leave();
//    return 0;
// }
// #endif

