#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>


using std::cout;
using std::endl;


// #define MAXBLOCKSZ 256

void leave();
// void leave() {
//    if ( CUBLAS_STATUS_SUCCESS != cublasShutdown()) printf("cublas down err\n");
//    cudaDeviceReset();
//    exit(0);
// }


__device__
int mone_on_l(int l) { return (1-2*(l&1));}


// __device__ void ropqln(int m, int l, double r2, double z, double cx[3], double q[2], int &kk) {
// // C- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
// // c  Returns kk, which points to the current component of q.
// // C  These subroutines are utility routines called by ropyln.f.
// // C  Kept separate from ropyln because some optimizing compilers have bugs
// // C  (e.g. intel ifort version 11).
// // C  These routines are the time-critical steps.
//
//    int mm,k2,k1;
//    double a,b,xx,yy;
//
//    if (l == m) {
//       a = 1.0;
//       for ( mm = 0; mm <= m-1; ++mm) a = a*(2*mm+1);
//       kk = 0;
//       a = a*cx[0];
//       q[kk] = a;
//    } else if (l == m+1) {
//       b = 1.0;
//       for ( mm = 0; mm <= m; ++mm) b = b*(2*mm+1);
//       b = b*cx[0];
//       kk = 1;
//       q[kk] = b*z;
//    } else if (l >= m+2) {
//       k2 = kk;
//       k1 = kk+1;
//       if (k1 == 2) k1 = 0;
//       xx = -(l+m-1.0)/(l-m)*cx[0]/cx[2];
//       yy = (2*l-1.0)/(l-m)*cx[0]/cx[1];
//       q[k1] = xx*r2*q[k1]+yy*z*q[k2];
//       kk = k1;
//    }
//
// }
//
// __device__
// void ropynx(int m, int l, int kk, double q[2], double cm, double sm, double *yl) {
//       int lav = l*(l+1);
//       yl[lav+m] = cm*q[kk];
//       if (m == 0) return;
//       yl[lav-m] = sm*q[kk];
// }

// __device__
// void ropcsm(int m, double x, double y, double w, double &cm, double &sm) {
// // C- Makes cm and sm. Must be called in sequence m=0,1,2...
//
//    if (m == 0) {
//       cm = 1.0;
//       sm = 0.0;
//    } else if (m == 1) {
//       cm = x;
//       sm = y;
//    } else if (m >= 2) {
//       w  = cm;
//       cm = x*cm - y*sm;
//       sm = y*w + x*sm;
//    }
// }

__device__
void ropyln(double x, double y, double z, int lmax, int *df, double *yl, double &rsq) {
// C- Normalized spheric harmonic polynomials (vectorizes).
// C ----------------------------------------------------------------------
// Ci Inputs
// Ci   n     :number of points for which to calculate yl
// Ci   x     :x component of Cartesian coordinate
// Ci   y     :y component of Cartesian coordinate
// Ci   z     :z component of Cartesian coordinate
// Ci   lmax  :compute yl for l=0..lmax
// Ci   nd    :leading dimension of yl; nd must be >= n
// Co Outputs
// Co   yl    :Ylm(i,ilm) are the (real) spherical harmonics
// Co         :for l=0..lmax and point i.
// Co   rsq   :rsq(i) square of length of point i
// Cr Remarks
// Cr   yl = real harmonics (see Takao's GW note) * r^l
// Cu Updates
// Cu  24 Apr 09 (Lozovoi) replace w() with allocate
// Cu  25 Jun 03 (Kino) initialize cx to zero
// C ----------------------------------------------------------------------
   int m,l,kk;
// c     integer i,ocm,osm,oq,oh
//       yl(nd,(lmax+1)**2)
   double cx[3];
   double fpi,f2m;

   fpi = 4.0*M_PI;
   double cm, sm, q[2], h;

   rsq = x*x + y*y + z*z;
   for (int i = 0; i < 3; ++i) cx[i] = 0.0;

// C --- Loop over m: cx are normalizations for l, l-1, l-2 ---
   f2m = 1.0;
   for (m = 0; m <= lmax; ++m) {
      if (m == 0) {
         cx[0] = sqrt(1/fpi);
         cm = 1.0;
         sm = 0.0;
      } else {
         if (m == 1) {
            cm = x;
            sm = y;
         } else if (m >= 2) {
            h  = cm;
            cm = x*cm - y*sm;
            sm = y*h + x*sm;
         }
         f2m = f2m*(2*m*(2*m-1));
         cx[0] = sqrt(((2*m+1)*2)/(fpi*f2m));
      }
      for ( l = m; l <= lmax; ++l) {
         if (l == m) {
            kk = 0;
            q[kk] = cx[0]*df[l];
         } else if (l == m+1) {
            kk = 1;
            q[kk] = z*cx[0]*df[l];
         } else if (l >= m+2) {
            int k2 = kk;
            int k1 = kk+1;
            if (k1 == 2) k1 = 0;
            double xx = -(l+m-1.0)/(l-m)*cx[0]/cx[2];
            double yy = (2*l-1.0)/(l-m)*cx[0]/cx[1];
            q[k1] = xx*rsq*q[k1]+yy*z*q[k2];
            kk = k1;
         }

         int lav = l*(l+1);
         yl[lav+m] = cm*q[kk];
         if (m != 0) yl[lav-m] = sm*q[kk];

         cx[2] = cx[1];
         cx[1] = cx[0];
         cx[0] = cx[0]*sqrt(double((l+1-m)*(2*l+3))/double((l+1+m)*(2*l+1)));
      }
   }
}



template<int llmaxl>
__device__
void soldhj(double r[3], int lmax, int df[llmaxl+1], double hl[(llmaxl+1)*(llmaxl+1)]){
// C- Real solid hankel and bessel functions.
// C ----------------------------------------------------------------
// Ci Inputs
// Ci   r     :radius (or radius/avw using OKA's conventions)
// Ci   e     :energy (or energy*avw**2 using OKA's conventions)
// Ci   loka  :conventions for Hankels, Bessels; see besslr
// Ci   lmax  :maximum l for a given site
// Co Outputs
// Co   HL,BL: Hankel and Bessel functions:  calculated to (lmax+1)**2
// Cr Remarks
// Cr   Generates bessel function * real spherical harmonic
// Cr   MSM's standard defs, notes IV-43.
// Cu Updates
// Cu   19 May 04 Changed loka from logical to integer
// C ----------------------------------------------------------------

   int tid = threadIdx.y*blockDim.x + threadIdx.x;
   int nths = blockDim.x*blockDim.y;

   int const llmax = (llmaxl+1)*(llmaxl+1);
   double dl[llmax];
   for (int i = 0; i < 100; ++i) dl[i] = 0.0;

   double r2;
// C     call sylm(r,hl,lmax,r2)
   ropyln(r[0],r[1],r[2],lmax,df,dl,r2);
   int ilm = 0;
   double rfac = sqrt(r2);
   if (r2 < 1e-10) r2 = 1.0;
   for ( int l = 0; l <= lmax; ++l) {
// C       rfac = 1/r**(2l+1), or 1/(r/w)**(2l+1) using OKA conventions
      rfac /= r2;
      for ( int m = -l; m <= l; ++m ) {
         dl[ilm] *= rfac*df[l];
         ++ilm;
      }
   }

   for (int i = tid; i < 100; i += nths) hl[i] = dl[i];

}


template<int llmaxl>
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
   double c[llmaxl+1], s[llmaxl+1], p[llmaxl+1][llmaxl+1];

   double &x = c[1];
   double &y = s[1];
   double &z = p[1][0];

   c[0] = 1.0;
   s[0] = 0.0;
   p[0][0] = 1.0;
   p[1][1] = 1.0;

   n = (lmx+1)*(lmx+1);
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


__device__
double inclmod1(double a) {
// ! Alternative implementation. Should give the same result
// !       s = sign(1.0_dp, a)
// !       r = -s*(modulo(-s*a+0.5_dp, 1.0_dp) - 0.5_dp)
// ! plot it if not clear

// alternative: a - rint(a). rint is implementation and settings dependent...

// this whole block finds the nearest integer with the bounday case of 0.5 rounded towards 0;
//    double s = copysign(1.0, a);
//    double t = a + 0.5*s;
//    double r = trunc(t);
//    if (t == r) r = r - s;

   double r;

   if (a < 0.0) {
      double t = a - 0.5;
      r = trunc(t);
      if (t == r) r = r + 1.0;
   } else {
      double t = a + 0.5;
      r = trunc(t);
      if (t == r) r = r - 1.0;
   }

   return a - r;
}

__device__
void ploughmans_3x1_mv(double *m, double u[3], double v[3]) {
//    follows conventional mm in F so the access looks improper here in C
   for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) v[j] += m[3*i+j]*u[i];
}

__device__
void shorten(double p[3], double *plat, double *ilat) {

   double r[3];
   for (int i = 0; i < 3; ++i) r[i] = 0.0;
   ploughmans_3x1_mv(ilat, p, r);

   for (int i = 0; i < 3; ++i) r[i] = inclmod1(r[i]);

   for (int i = 0; i < 3; ++i) p[i] = 0.0;
   ploughmans_3x1_mv(plat, r, p);
}

template<int llmax, int nths>
__device__
void hstr_cu( bool pv, double *strux, double *dstrx,
              int nlmi1, int nlmi, int nlmj,
              double *hl, double *cg, int *indxcg, int *jcg,
              bool mol, int ldip, double vol, int *df) {
//- Make structure constants from reduced strux at energy zero
// ----------------------------------------------------------------------
//i Inputs:
//i  MOL   : if T skip dipole correction (use non-periodic setup)
//i  pv    :if T calculate pressure
//i  ldip  : 3 include 'spherical' dipole correction to Ewald
//i        : 2 include 'slab' dipole correction to Ewald
//i        : any other number - skip the dipole correction
//i  nlmf,nlm :make B_LL' for L = 1, nlmf, L'= 1, nlm
//i  nlmq1 :leading dimensions of B, nlmq1=(ll(nlmq)+1)**2
//i  nlmq  :max L-cutoff for multipoles, leading dimension of B'
//i  hl    :radial Hankels for given point |tau|
//i  cg,indxcg,jcg : Clebsch Gordan coefficients in condensed form
//i  vol   :unit cell volume
//o Outputs: strx,drstrx
//o  strx(nlmq1,nlm)  :coefficients B_LL'(tau) (structure constants)
//o  dstrx(nlmq,nlm) :derivatives of structure constants (x dB_LL'/dx) at x = tau
//o         if pv = F dstrx is not touched
//r Remarks: B_LL'(tau) = 4\pi \sum_L" (-1)^l' H_l"(tau) Y_L"(tau)
//r         HY are reduced structure constants (see rcnsl0, soldh)
//r         If pv is set, returns tau*d/dtau B in drstrx
//u Updates
//u   28 Apr 10 (SL) Slab dipole correction to Ewald
// ----------------------------------------------------------------------

   int ilm, jlm, li, lj, ii, indx, icg1, icg2, llm, lp;

//    initializer not allowed for __shared__ variables...
   int const msm[8] = {2,3,1,4,5,8,6,7};

   int tid = threadIdx.y*blockDim.x + threadIdx.x;


   __syncthreads();

   __shared__ int ll[llmax];

   for (int lm = tid; lm < llmax; lm += nths) ll[lm] = int(sqrtf(lm));

   __syncthreads();

//    int bid = blockIdx.y * gridDim.x + blockIdx.x;
//    if (tid == 0 && bid == 0) {
//       printf("ll:");
//       for (int i = 0; i < llmax; ++i) printf(" %d",ll[i]);
//       printf("\n");
//    }

   int els = nlmj*nlmi1;

   for (int iel = tid; iel < els; iel += nths) {

      ilm = iel/nlmj;
//       jlm = els%nlm;
      jlm = iel - ilm*nlmj;

      li = ll[ilm];
      lj = ll[jlm];

      bool npv = pv && (ilm < nlmi);
      double struxel = 0.0;
      double dstrxel = 0.0;

   //    FILE *fl = fopen("hstrinternal","a");

   // --- add together Gaunt-sums ---
      ii = ilm > jlm ? ilm : jlm;
      indx = (ii*(ii+1))/2 + (ilm > jlm ? jlm : ilm);
      icg1 = indxcg[indx];
      icg2 = indxcg[indx+1]-1;
      for (int icg = icg1 - 1; icg < icg2; ++icg) {
         llm = jcg[icg]-1;
         lp = ll[llm];
         if (li + lj == lp) {
            double t = cg[icg]*hl[llm];
            struxel += t;
            if (npv) dstrxel -= (lp+1) * t;
         }
      }


   // --- the following includes extra p terms 'implicitly' ---
      if ((ilm == jlm) && (li == 1) && ((ldip == 2) || (ldip == 3)) && (!mol)) {
         struxel += 1.0/vol;
         if (pv) dstrxel -= 3.0/vol;
      }

      double fac = 4.0 * M_PI * mone_on_l(lj) * sqrt(double((2*lj+1)*(2*li+1)))/double(df[lj]*df[li]);

      if (ilm > 0 && ilm < 9) ilm = msm[ilm-1];
      if (jlm > 0 && jlm < 9) jlm = msm[jlm-1];

      strux[ilm*nlmj + jlm] = struxel *fac;
      if (npv) dstrx[ilm*nlmj + jlm] = dstrxel *fac;
   }
//    fclose(fl);
//    __syncthreads();
//    return 0;
}

template<int llmaxl>
__device__
void rcsl01_cu(double tau[3], double a, int lmax, double alat, double *rlat, int nkr, double vol, double *dl) {
//  k-space part of reduced structure constants for e=0 and q=0
   double r2;

   int ilm,nlm;
   double fpibv,gamma,scalp,tpiba,eiphi[2];
   double r[3], yl[(llmaxl+1)*(llmaxl+1)];
   double const tpi = 2*M_PI;

   int tid = threadIdx.y*blockDim.x + threadIdx.x;
   int tot = blockDim.x * blockDim.y;


   if (lmax > llmaxl) {
      printf("rcnsl0: increase lmaxx\n");
      return;
   }

   gamma = 0.25/(a*a);
   fpibv = 2.0*tpi/vol;
   tpiba = tpi/alat;
   nlm = (lmax+1)*(lmax+1);
   for (int ilm = 0; ilm < nlm; ++ilm) dl[ilm] = 0.0;
//    dl[0] = -fpibv*gamma;
   if (tid == 0) dl[0] = -fpibv*gamma;
//    for ( int ir = 1; ir < nkr;  ++ir ) {
   for ( int ir = 1 + tid; ir < nkr;  ir += tot ) {
      for (int i = 0; i < 3; ++i) r[i] = tpiba*rlat[3*ir + i];
//       printf(" rlat: %20.12f %20.12f %20.12f\n",rlat[3*ir + 0], rlat[3*ir + 1], rlat[3*ir + 2]);
//       printf("r: %20.12f %20.12f %20.12f\n",r[0],r[1],r[2]);
      scalp = 0.0; for (int i = 0; i < 3; ++i) scalp += r[i]*tau[i];
      scalp *= alat;
      eiphi[0] = cos(scalp);
      eiphi[1] = sin(scalp);
//       printf("eiphi: %20.12f %20.12f\n",eiphi[0],eiphi[1]);
      r2 = sylm_cu<llmaxl>(r,yl,lmax);
//       printf("  r2: %20.12f\n",r2);
      double yyy = fpibv*exp(-gamma*r2)/r2;
      ilm = 0;
      for ( int l = 0; l <= lmax; ++l ) {
         double yyye0 = yyy*eiphi[0];
         for ( int m = 1; m <= 2*l+1; ++m ) {
            dl[ilm] += yl[ilm]*yyye0;
//             printf("dl: %20.12f\n",dl[ilm]);
            ++ilm;
         }
// eiphi *= (0,-1):
         double xxx = eiphi[0];
         eiphi[0] = eiphi[1];
         eiphi[1] = -xxx;
      }
   }
}



template<int llmaxl>
__device__
void rcsl02_cu(double tau[3], double a, int lmax, double alat, double *dlat, int nkd, double *dl) {
//  real space summation
   double r2;
   int ilm,ir1;
   double a2,cc,gl,r1,ta2,dl0;
   double r[3],yl[(llmaxl+1)*(llmaxl+1)],chi;//chi[LMAX+1];
   double const twoinvsqpi = M_2_SQRTPI; //1.12837916709551257390_8 ! 2/sqrt(pi)

   if (lmax > llmaxl) {
      printf("rcnsl0: increase lmaxx\n");
      return;
   }

   ir1 = 1;
   if ((tau[0]*tau[0]+tau[1]*tau[1]+tau[2]*tau[2]) > 1e-6) ir1=0;
//       if (sum(tau*tau) > 1d-6) ir1=1

   int tid = threadIdx.y * blockDim.x + threadIdx.x;
   int tot = blockDim.x * blockDim.y;

   if (lmax > 0) {
      a2 = a*a;
      ta2 = 2*a2;
      cc = ta2*a*twoinvsqpi;
//       for ( int ir = ir1; ir < nkd; ++ir ) {
      for ( int ir = ir1 + tid; ir < nkd; ir += tot ) {
         for (int i = 0; i < 3; ++i) r[i] = alat*(tau[i]-dlat[3*ir + i]);

         r2 = sylm_cu<llmaxl>(r,yl,lmax);

         double invr2 = 1.0/r2;
         r1 = sqrt(r2);
         chi = erfc(a*r1)/r1;
         gl = -cc*exp(-a2*r2)/ta2;
         dl[0] += yl[0]*chi;
         ilm = 1;
         for ( int l = 1; l <= lmax; ++l ) {
            chi = ((2*l-1)*chi - gl)*invr2;
            gl *= ta2;
            for ( int m = 1; m <= 2*l+1; ++m ) {
               dl[ilm] += yl[ilm]*chi;
               ++ilm;
            }
         }
      }
// ...In case lmax = 0 do everything explicitly
   } else {
      dl0 = 0.0;
//       for ( int ir = ir1; ir < nkd; ++ir) {
      for ( int ir = ir1 + tid; ir < nkd; ir += tot ) {
         for (int i = 0; i < 3; ++i) r[i] = tau[i]-dlat[3*ir + i];
         r2 = 0.0; for (int i = 0; i < 3; ++i) r2 += r[i]*r[i];
         r1 = alat*sqrt(r2);
         dl0 += erfc(a*r1)/r1;
      }
      dl[0] += dl0;
   }

// --- add dl3 for diagonal sructure constants ------
//    if (ir1 == 1) dl[0] -= a*twoinvsqpi;
   if (ir1 == 1 && tid == 0) dl[0] -= a*twoinvsqpi;
}




// __global__
template<int llmaxl, int nths>
__device__
void rcnsl0_cu(double tau[3], double a, int lmax, double alat,
                                     double *rlat, int nkr,
                                     double *dlat, int nkd,
                                     double vol, double *cy, double *hl) {
//  reduced structure constants on lattice for e=0 and q=0.
//  result is periodic sum of (2*l-1)!!*ylm/r**(l+1). additive
//  constant for l=0 is chosen so that function averages to zero.


   int tid = threadIdx.y*blockDim.x + threadIdx.x;


   int maxnk = nkd > nkd ? nkd : nkr;
   int maxusedth = maxnk < nths ? maxnk : nths;

   int nlm = (lmax+1)*(lmax+1);

//  insanity....  there is no atomicAdd for type double!!!!
// this is a shoddy reduce of the local dl to the shared hl;
//    __shared__ double *dl;
//    if (tid == 0) dl = (double*)malloc(maxusedth*nlm*sizeof(double));
   __shared__ double dl[nths][(llmaxl+1)*(llmaxl+1)];

   __syncthreads();

   if (tid < maxusedth) {
      rcsl01_cu<llmaxl>(tau,a,lmax,alat,rlat,nkr,vol,dl[tid]);
//    printf("rc01: %22.12f %22.12f %22.12f\n", dl[0], dl[1], dl[2]);

      rcsl02_cu<llmaxl>(tau,a,lmax,alat,dlat,nkd,dl[tid]);
//    printf("rc02: %22.12f %22.12f %22.12f\n", dl[0], dl[1], dl[2]);

   }

   __syncthreads();

   for (int ilm = tid; ilm < nlm; ilm += nths) {
      hl[ilm] = 0.0;
      for (int th = 0; th < maxusedth; ++th) hl[ilm] += dl[th][ilm];
      hl[ilm] *= cy[ilm];
   }

   __syncthreads();

}



template<int lmaxl, int nths>
__global__ void
// void __launch_bounds__(MAXBLOCKSZ)
mkstrxd_cu(int ldip, int nbas, int nbas1, int ib0,
                           double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
                           double *dlat, int nkd, double *glat, int nkg,
                           int *indxcg, int *jcg, double *cg, double *cy,
                           bool pv, bool mol, bool force, double *strux, double *dstrx, int *sidx, int *didx,
                           double *plat, double *qlat) {

// strx(nlmq,nlmq1,nbas,nbas1),dstrx(nlmq,nlmq,nbas,nbas1)

//    printf("gs: (%3d,%3d), block: (%3d,%3d), thread (%2d,%2d) start\n",
//          gridDim.x, gridDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);

//    return;

   int ib = blockIdx.y;
   int jb = blockIdx.x;
   int tid = threadIdx.y*blockDim.x + threadIdx.x;

//    bool rth = blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0
//            && threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;

//    if (tid == 0) printf("block: %d %d %d speaking\n", blockIdx.x,blockIdx.y,blockIdx.z );
   int const llmax = (2*lmaxl+2)*(2*lmaxl+2);

   int ibr, li, li1, lj, nlmj, nlmi, nlmi1, lmxst; //lmxf,  lmax
   double tau[3];
   __shared__ double hl[llmax];
   __shared__ int df[2*lmaxl+2];
//    int const df[]= {1,1,3,15,105,945,10395,135135,2027025,34459425};

//       df(l) = (2l+1)!!
   if (tid == 0) {
      df[0] = 1;
      for (int l = 0; l < 2*lmaxl+1; ++l) df[l+1] = df[l]*(2*l+1);
   }
   __syncthreads();


   ibr = ib + ib0;
   li = lmxl[ipc[ibr]-1];
//c       nlmi = (li+1)**2
//C...  nlmq1 == nlmq if no force or pressure are required.
//C     Otherwise need to go up to lmax+1 in the first index of B_LL'

//    int li1 =  (nlmq1 > nlmq) ? li + 1 : li;



//C...  The lines below are because we shall be filling the jb > ib triangle
//C     using symmetry properties of B (see tbesel.f) AND because we need l+1
//C     for forces. No complications if forces are not required.
//    if (nlmq1 > nlmq) {
//       lmax = max(li,lj);
// //       lmax = li;
// //             lmxf = lmax+1;
//       nlx = (lmax+1)*(lmax+1);
//       nlf = (lmax+2)*(lmax+2);
//       lmxst = 2*lmax+1;
//
//
//    } else {
//       nlx = (lj+1)*(lj+1);
//       nlf = (li+1)*(li+1);
//       nlf = max(nlf,nlx);
//       lmxst = li + lj;
//    }

//    if (nlmq1 > nlmq) ++li;
//    lmxst = li + lj;
//    nlf = (li+1)*(li+1);
//    nlx = (lj+1)*(lj+1);

   nlmi  = (li+1)*(li+1);
   nlmi1 = nlmi;
   li1 = li;

   if (force || pv) {
      li1 = ++li;
      nlmi1  = (li1+1)*(li1+1);
   }

   lj = lmxl[ipc[jb]-1];
   nlmj = (lj+1)*(lj+1);

   lmxst = li1 + lj;

//          printf("block: (%3d,%3d), thread (%2d,%2d) working\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);

   for (int i = 0; i < 3; ++i) tau[i] = bas[3*jb+i] - bas[3*ibr+i];

   if (!mol) {
      shorten(tau,plat,qlat);
      rcnsl0_cu<(2*lmaxl+1), nths>(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl);
   } else {
      for (int i = 0; i < 3; ++i) tau[i] *= alat;
      soldhj<(2*lmaxl+1)>(tau,lmxst,df,hl);
   }

   __syncthreads();

   int si;
   int di = 0;

   if (tid < nlmj*nlmi1) {
      si = sidx[nbas*ib + jb]-1;
      if (pv) di = didx[nbas*ib + jb]-1;
      hstr_cu<llmax,nths>(pv, &strux[si], &dstrx[di], nlmi1, nlmi, nlmj, hl, cg, indxcg, jcg, mol, ldip, vol, &df[1]);
   }
}



extern "C"
void mkstrxd_hcu(int ldip, int nbas, int nbas1, int ib0, int lmaxl,
                  double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
                  double *dlat, int nkd, double *glat, int nkg,
                  int *indxcg, int *jcg, double *cg, double *cy,
                  int fpv, int fmol, int fforce, double *strx, double *dstrx,
                                    int *sidx, int *didx, int nstrx, int ndstrx,
                  double *plat, double *qlat) {

   int *d_ipc, *d_lmxl, *d_indxcg, *d_jcg, *d_sidx, *d_didx;
   double *d_bas, *d_dlat, *d_glat, *d_cg, *d_cy, *d_strx, *d_dstrx, *d_plat, *d_qlat;



   bool mol = fmol != 0;
   bool pv = fpv != 0;
   bool force = fforce != 0;

   cudaError_t err;

   size_t r8 = sizeof(double);
   size_t i4 = sizeof(int);

//    int nstrx = nlmq*nlmq1*nbas*nbas1;
//    int ndstrx = nlmq*nlmq*nbas*nbas1;
   int idxsz = nbas*nbas1;

   err = cudaMalloc((void **) &d_ipc   , nbas  *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_lmxl  , nclas *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_indxcg, 7400  *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_jcg   , 62200 *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_sidx  , idxsz *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_didx  , idxsz *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_bas   , 3*nbas*r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_dlat  , 3*nkd *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_glat  , 3*nkg *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_cg    , 62200 *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_cy    , 289   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_strx  , nstrx *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_dstrx , ndstrx*r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_plat  , 3*3   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMalloc((void **) &d_qlat  , 3*3   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);

   err = cudaMemcpy(d_ipc   , ipc   , nbas  *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_lmxl  , lmxl  , nclas *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_indxcg, indxcg, 7400  *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_jcg   , jcg   , 62200 *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_sidx  , sidx  , idxsz *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_didx  , didx  , idxsz *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_bas   , bas   , 3*nbas*r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_dlat  , dlat  , 3*nkd *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_glat  , glat  , 3*nkg *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_cg    , cg    , 62200 *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_cy    , cy    , 289   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_plat  , plat  , 3*3   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(d_qlat  , qlat  , 3*3   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);

   err = cudaMemset ( d_strx , 0, nstrx *r8); if (err != cudaSuccess) printf("cudaMemset: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemset ( d_dstrx, 0, ndstrx*r8); if (err != cudaSuccess) printf("cudaMemset: err: %d, line(%d)\n", err, __LINE__);



   dim3 blocks (nbas, nbas1);
//    dim3 threads(bsize);
//    dim3 threads(bsize);

//
//    int const bszs[8] = {4,9,16,36,64,81,144,256};
//    for (int i = 0; i < 8; ++i)
//       if (nlmq*nlmq1 == bszs[i])
//          mkstrx_cu <bszs[i]> <<<blocks, threads>>> (ldip, nbas, nbas1, ib0, nlmq1, nlmq,d_bas, d_ipc, nclas, d_lmxl, awld, alat, vol,
//                              d_dlat, nkd, d_glat, nkg,d_indxcg, d_jcg, d_cg, d_cy, pv, mol, d_strx, d_dstrx, d_plat, d_qlat);



#define mkstrxd_cup(lmaxl,bsize) \
   printf(" mkstrx_cu <%d,%d> <<<(%d,%d),(%d)>>>\n", lmaxl, bsize, blocks.x, blocks.y, bsize); \
   mkstrxd_cu <(lmaxl),(bsize)> <<<blocks, bsize>>> (ldip, nbas, nbas1, ib0,d_bas, d_ipc, nclas, d_lmxl, awld, alat, vol, \
            d_dlat, nkd, d_glat, nkg,d_indxcg, d_jcg, d_cg, d_cy, pv, mol, force, d_strx, d_dstrx, d_sidx, d_didx, d_plat, d_qlat)

// the abreviations are
//    lmaxl : the maximal small l reached by any atom.
//    lmmax : the maximal capital L reached
//    llmaxl: the maximal small l  for the expansion term (2*lmaxl + 1)
//    llmax : the maximal capital L=lm  for the expansion term (2*lmaxl + 2)**2

//    int lmaxl = lround(sqrt(nlmq1))-(nlmq1>nlmq?2:1);
   if (!mol) {
      int const bsize = 32;
      switch (lmaxl) {
         case 0: mkstrxd_cup(0,bsize); break;
         case 1: mkstrxd_cup(1,bsize); break;
         case 2: mkstrxd_cup(2,bsize); break;
         case 3: mkstrxd_cup(3,bsize); break;
         case 4: mkstrxd_cup(4,bsize); break;
         default:
            printf(" mkstrx_cu: lmaxl: %d, lmaxl>4 not supported!\n", lmaxl);
            leave();
      }
   } else {
      int const bsize = 4;
      switch (lmaxl) {
         case 0: mkstrxd_cup(0,bsize); break;
         case 1: mkstrxd_cup(1,bsize); break;
         case 2: mkstrxd_cup(2,bsize); break;
         case 3: mkstrxd_cup(3,bsize); break;
         case 4: mkstrxd_cup(4,bsize); break;
         default:
            printf(" mkstrx_cu: lmaxl: %d, lmaxl>4 not supported!\n", lmaxl);
            leave();
      }
   }

//    cudaDeviceSynchronize();

   err = cudaMemcpy(strx  , d_strx  , nstrx *r8, cudaMemcpyDeviceToHost); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
   err = cudaMemcpy(dstrx , d_dstrx , ndstrx*r8, cudaMemcpyDeviceToHost); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);

//    cudaDeviceSynchronize();

   err = cudaFree(d_ipc   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_lmxl  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_indxcg); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_jcg   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_sidx  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_didx  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_bas   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_dlat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_glat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_cg    ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_cy    ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_strx  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_dstrx ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_plat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
   err = cudaFree(d_qlat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);


//    FILE *fl = fopen("strx-cu", "w");
//    for ( int ib = 0; ib < nbas1; ++ib ) {
//       for ( int jb = 0; jb < nbas; ++jb ) {
//          fprintf(fl, "jb, ib: %4d %4d\n", jb, ib);
//          for (int j = 0; j < nlmq1 ; ++j) {
//             for (int i = 0; i < nlmq; ++i) {
//                int idx = (nbas*ib + jb)*nlmq*nlmq1 + nlmq*j+i;
//                fprintf(fl, " %16.8f", strx[idx]);
//             }
//             fprintf(fl, "\n");
//          }
//       }
//    }
//    fflush(fl);
//    fclose(fl);
//    leave();

//    leave();
}

































// template<int lmaxl, int nths>
// __global__ void
// // void __launch_bounds__(MAXBLOCKSZ)
// mkhl_ewld_cu(int ldip, int nbas, int nbas1, int ib0, int nlmq1, int nlmq,
//                            double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
//                            double *dlat, int nkd, double *glat, int nkg,
//                            int *indxcg, int *jcg, double *cg, double *cy,
//                            bool pv, bool mol, double *strux, double *dstrx,
//                            double *plat, double *qlat) {
//
// // strx(nlmq,nlmq1,nbas,nbas1),dstrx(nlmq,nlmq,nbas,nbas1)
//
// //    printf("gs: (%3d,%3d), block: (%3d,%3d), thread (%2d,%2d) start\n",
// //          gridDim.x, gridDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
// //    return;
//
//    int ib = blockIdx.y;
//    int jb = blockIdx.x;
//    int tid = threadIdx.y*blockDim.x + threadIdx.x;
//
// //    bool rth = blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0
// //            && threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;
//
// //    if (tid == 0) printf("block: %d %d %d speaking\n", blockIdx.x,blockIdx.y,blockIdx.z );
//    int const llmax = (2*lmaxl+2)*(2*lmaxl+2);
//
//    int ibr, li, lj, lmax, nlx, nlf, lmxst; //lmxf
//    double tau[3];
//    __shared__ double hl[llmax];
//    __shared__ int df[2*lmaxl+2];
// //    int const df[]= {1,1,3,15,105,945,10395,135135,2027025,34459425};
//
// //       df(l) = (2l+1)!!
//    if (tid == 0) {
//       df[0] = 1;
//       for (int l = 0; l < 2*lmaxl+1; ++l) df[l+1] = df[l]*(2*l+1);
//    }
//    __syncthreads();
//
//
//    ibr = ib + ib0;
//    li = lmxl[ipc[ibr]-1];
// //c       nlmi = (li+1)**2
// //C...  nlmq1 == nlmq if no force or pressure are required.
// //C     Otherwise need to go up to lmax+1 in the first index of B_LL'
//
// //    int li1 =  (nlmq1 > nlmq) ? li + 1 : li;
//
//    lj = lmxl[ipc[jb]-1];
//
// //C...  The lines below are because we shall be filling the jb > ib triangle
// //C     using symmetry properties of B (see tbesel.f) AND because we need l+1
// //C     for forces. No complications if forces are not required.
//    if (nlmq1 > nlmq) {
//       lmax = max(li,lj);
// //             lmxf = lmax+1;
//       nlx = (lmax+1)*(lmax+1);
//       nlf = (lmax+2)*(lmax+2);
//       lmxst = 2*lmax+1;
//    } else {
//       nlx = (lj+1)*(lj+1);
//       nlf = (li+1)*(li+1);
//       nlf = max(nlf,nlx);
//       lmxst = li + lj;
//    }
//
//
// //          printf("block: (%3d,%3d), thread (%2d,%2d) working\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
//    for (int i = 0; i < 3; ++i) tau[i] = bas[3*jb+i] - bas[3*ibr+i];
//
//    shorten(tau,plat,qlat);
//    rcnsl0_cu<(2*lmaxl+1), nths>(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl);
//
//    __syncthreads();
//
// }
//
//
//
//
//
//
// template<int lmaxl, int nths>
// __global__ void
// // void __launch_bounds__(MAXBLOCKSZ)
// mkhl_mol_cu(int ldip, int nbas, int nbas1, int ib0, int nlmq1, int nlmq,
//                            double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
//                            double *dlat, int nkd, double *glat, int nkg,
//                            int *indxcg, int *jcg, double *cg, double *cy,
//                            bool pv, bool mol, double *strux, double *dstrx,
//                            double *plat, double *qlat) {
//
// // strx(nlmq,nlmq1,nbas,nbas1),dstrx(nlmq,nlmq,nbas,nbas1)
//
// //    printf("gs: (%3d,%3d), block: (%3d,%3d), thread (%2d,%2d) start\n",
// //          gridDim.x, gridDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
// //    return;
//
//    int ib = blockIdx.y;
//    int jb = blockIdx.x;
//    int tid = threadIdx.y*blockDim.x + threadIdx.x;
//
// //    bool rth = blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0
// //            && threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;
//
// //    if (tid == 0) printf("block: %d %d %d speaking\n", blockIdx.x,blockIdx.y,blockIdx.z );
//    int const llmax = (2*lmaxl+2)*(2*lmaxl+2);
//
//    int ibr, li, lj, lmax, nlx, nlf, lmxst; //lmxf
//    double tau[3];
//    __shared__ double hl[llmax];
//    __shared__ int df[2*lmaxl+2];
// //    int const df[]= {1,1,3,15,105,945,10395,135135,2027025,34459425};
//
// //       df(l) = (2l+1)!!
//    if (tid == 0) {
//       df[0] = 1;
//       for (int l = 0; l < 2*lmaxl+1; ++l) df[l+1] = df[l]*(2*l+1);
//    }
//    __syncthreads();
//
//
//    ibr = ib + ib0;
//    li = lmxl[ipc[ibr]-1];
// //c       nlmi = (li+1)**2
// //C...  nlmq1 == nlmq if no force or pressure are required.
// //C     Otherwise need to go up to lmax+1 in the first index of B_LL'
//
// //    int li1 =  (nlmq1 > nlmq) ? li + 1 : li;
//
//    lj = lmxl[ipc[jb]-1];
//
// //C...  The lines below are because we shall be filling the jb > ib triangle
// //C     using symmetry properties of B (see tbesel.f) AND because we need l+1
// //C     for forces. No complications if forces are not required.
//    if (nlmq1 > nlmq) {
//       lmax = max(li,lj);
// //             lmxf = lmax+1;
//       nlx = (lmax+1)*(lmax+1);
//       nlf = (lmax+2)*(lmax+2);
//       lmxst = 2*lmax+1;
//    } else {
//       nlx = (lj+1)*(lj+1);
//       nlf = (li+1)*(li+1);
//       nlf = max(nlf,nlx);
//       lmxst = li + lj;
//    }
//
//
// //          printf("block: (%3d,%3d), thread (%2d,%2d) working\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
//    for (int i = 0; i < 3; ++i) tau[i] = bas[3*jb+i] - bas[3*ibr+i];
//
//    if (!mol) {
//       shorten(tau,plat,qlat);
//       rcnsl0_cu<(2*lmaxl+1), nths>(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl);
//    } else {
//       for (int i = 0; i < 3; ++i) tau[i] *= alat;
//       soldhj<(2*lmaxl+1)>(tau,lmxst,df,hl);
//    }
//
//    __syncthreads();
//
//    if (tid < nlf*nlx) {
//       hstr_cu<llmax,nths>(pv, &strux[(nbas*ib + jb)*nlmq*nlmq1], &dstrx[(nbas*ib + jb)*nlmq*nlmq],
//                                           nlf, nlx, nlmq1, nlmq, hl, cg, indxcg, jcg, mol, ldip, vol, &df[1]);
//    }
// }
//
//
//
// template<int lmaxl, int nths>
// __global__ void
// // void __launch_bounds__(MAXBLOCKSZ)
// mkstrhl_cu(int ldip, int nbas, int nbas1, int ib0, int nlmq1, int nlmq,
//                            double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
//                            double *dlat, int nkd, double *glat, int nkg,
//                            int *indxcg, int *jcg, double *cg, double *cy,
//                            bool pv, bool mol, double *strux, double *dstrx,
//                            double *plat, double *qlat) {
//
// // strx(nlmq,nlmq1,nbas,nbas1),dstrx(nlmq,nlmq,nbas,nbas1)
//
// //    printf("gs: (%3d,%3d), block: (%3d,%3d), thread (%2d,%2d) start\n",
// //          gridDim.x, gridDim.y, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
// //    return;
//
//    int ib = blockIdx.y;
//    int jb = blockIdx.x;
//    int tid = threadIdx.y*blockDim.x + threadIdx.x;
//
// //    bool rth = blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0
// //            && threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;
//
// //    if (tid == 0) printf("block: %d %d %d speaking\n", blockIdx.x,blockIdx.y,blockIdx.z );
//    int const llmax = (2*lmaxl+2)*(2*lmaxl+2);
//
//    int ibr, li, lj, lmax, nlx, nlf, lmxst; //lmxf
//    double tau[3];
//    __shared__ double hl[llmax];
//    __shared__ int df[2*lmaxl+2];
// //    int const df[]= {1,1,3,15,105,945,10395,135135,2027025,34459425};
//
// //       df(l) = (2l+1)!!
//    if (tid == 0) {
//       df[0] = 1;
//       for (int l = 0; l < 2*lmaxl+1; ++l) df[l+1] = df[l]*(2*l+1);
//    }
//    __syncthreads();
//
//
//    ibr = ib + ib0;
//    li = lmxl[ipc[ibr]-1];
// //c       nlmi = (li+1)**2
// //C...  nlmq1 == nlmq if no force or pressure are required.
// //C     Otherwise need to go up to lmax+1 in the first index of B_LL'
//
// //    int li1 =  (nlmq1 > nlmq) ? li + 1 : li;
//
//    lj = lmxl[ipc[jb]-1];
//
// //C...  The lines below are because we shall be filling the jb > ib triangle
// //C     using symmetry properties of B (see tbesel.f) AND because we need l+1
// //C     for forces. No complications if forces are not required.
//    if (nlmq1 > nlmq) {
//       lmax = max(li,lj);
// //             lmxf = lmax+1;
//       nlx = (lmax+1)*(lmax+1);
//       nlf = (lmax+2)*(lmax+2);
//       lmxst = 2*lmax+1;
//    } else {
//       nlx = (lj+1)*(lj+1);
//       nlf = (li+1)*(li+1);
//       nlf = max(nlf,nlx);
//       lmxst = li + lj;
//    }
//
//
// //          printf("block: (%3d,%3d), thread (%2d,%2d) working\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
//
//    for (int i = 0; i < 3; ++i) tau[i] = bas[3*jb+i] - bas[3*ibr+i];
//
//    if (!mol) {
//       shorten(tau,plat,qlat);
//       rcnsl0_cu<(2*lmaxl+1), nths>(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl);
//    } else {
//       for (int i = 0; i < 3; ++i) tau[i] *= alat;
//       soldhj<(2*lmaxl+1)>(tau,lmxst,df,hl);
//    }
//
//    __syncthreads();
//
//    if (tid < nlf*nlx) {
//       hstr_cu<llmax,nths>(pv, &strux[(nbas*ib + jb)*nlmq*nlmq1], &dstrx[(nbas*ib + jb)*nlmq*nlmq],
//                                           nlf, nlx, nlmq1, nlmq, hl, cg, indxcg, jcg, mol, ldip, vol, &df[1]);
//    }
// }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// extern "C"
// void mkstrx_hcu_composite(int ldip, int nbas, int nbas1, int ib0, int nlmq1, int nlmq,
//                   double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
//                   double *dlat, int nkd, double *glat, int nkg,
//                   int *indxcg, int *jcg, double *cg, double *cy,
//                   int fpv, int fmol, double *strx, double *dstrx,
//                   double *plat, double *qlat) {
//
//    int *d_ipc, *d_lmxl, *d_indxcg, *d_jcg;
//    double *d_bas, *d_dlat, *d_glat, *d_cg, *d_cy, *d_strx, *d_dstrx, *d_plat, *d_qlat;
//
//
//
//    bool mol = fmol == 1;
//    bool pv = fpv == 1;
//
//    cudaError_t err;
//
//    size_t r8 = sizeof(double);
//    size_t i4 = sizeof(int);
//
//    int nstrx = nlmq*nlmq1*nbas*nbas1;
//    int ndstrx = nlmq*nlmq*nbas*nbas1;
//
//    err = cudaMalloc((void **) &d_ipc   , nbas  *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_lmxl  , nclas *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_indxcg, 7400  *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_jcg   , 62200 *i4); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_bas   , 3*nbas*r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_dlat  , 3*nkd *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_glat  , 3*nkg *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_cg    , 62200 *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_cy    , 289   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_strx  , nstrx *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_dstrx , ndstrx*r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_plat  , 3*3   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMalloc((void **) &d_qlat  , 3*3   *r8); if (err != cudaSuccess) printf("cudaMalloc: err: %d, line(%d)\n", err, __LINE__);
//
//    err = cudaMemcpy(d_ipc   , ipc   , nbas  *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_lmxl  , lmxl  , nclas *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_indxcg, indxcg, 7400  *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_jcg   , jcg   , 62200 *i4, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_bas   , bas   , 3*nbas*r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_dlat  , dlat  , 3*nkd *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_glat  , glat  , 3*nkg *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_cg    , cg    , 62200 *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_cy    , cy    , 289   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_plat  , plat  , 3*3   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(d_qlat  , qlat  , 3*3   *r8, cudaMemcpyHostToDevice); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//
//    err = cudaMemset ( d_strx , 0, nstrx *r8); if (err != cudaSuccess) printf("cudaMemset: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemset ( d_dstrx, 0, ndstrx*r8); if (err != cudaSuccess) printf("cudaMemset: err: %d, line(%d)\n", err, __LINE__);
//
//
//
//    dim3 blocks (nbas, nbas1);
// //    dim3 threads(bsize);
// //    dim3 threads(bsize);
//
// //
// //    int const bszs[8] = {4,9,16,36,64,81,144,256};
// //    for (int i = 0; i < 8; ++i)
// //       if (nlmq*nlmq1 == bszs[i])
// //          mkstrx_cu <bszs[i]> <<<blocks, threads>>> (ldip, nbas, nbas1, ib0, nlmq1, nlmq,d_bas, d_ipc, nclas, d_lmxl, awld, alat, vol,
// //                              d_dlat, nkd, d_glat, nkg,d_indxcg, d_jcg, d_cg, d_cy, pv, mol, d_strx, d_dstrx, d_plat, d_qlat);
//
//
//
// #define mkstrx_cup(lmaxl,bsize) \
//    printf(" mkstrx_cu <%d,%d> <<<(%d,%d),(%d)>>>\n", lmaxl, bsize, blocks.x, blocks.y, bsize); \
//    mkstrx_cu <(lmaxl),(bsize)> <<<blocks, bsize>>> (ldip, nbas, nbas1, ib0, nlmq1, nlmq,d_bas, d_ipc, nclas, d_lmxl, awld, alat, vol, \
//                              d_dlat, nkd, d_glat, nkg,d_indxcg, d_jcg, d_cg, d_cy, pv, mol, d_strx, d_dstrx, d_plat, d_qlat)
//
// // the abreviations are
// //    lmaxl : the maximal small l reached by any atom.
// //    lmmax : the maximal capital L reached
// //    llmaxl: the maximal small l  for the expansion term (2*lmaxl + 1)
// //    llmax : the maximal capital L=lm  for the expansion term (2*lmaxl + 2)**2
//
//    int lmaxl = lround(sqrt(nlmq1))-(nlmq1>nlmq?2:1);
//    if (!mol) {
//       int const bsize = 32;
//       switch (lmaxl) {
//          case 0: mkstrx_cup(0,bsize); break;
//          case 1: mkstrx_cup(1,bsize); break;
//          case 2: mkstrx_cup(2,bsize); break;
//          case 3: mkstrx_cup(3,bsize); break;
//          case 4: mkstrx_cup(4,bsize); break;
//          default:
//             printf(" mkstrx_cu: lmaxl: %d, lmaxl>4 not supported!\n", lmaxl);
//             leave();
//       }
//    } else {
//       int const bsize = 4;
//       switch (lmaxl) {
//          case 0: mkstrx_cup(0,bsize); break;
//          case 1: mkstrx_cup(1,bsize); break;
//          case 2: mkstrx_cup(2,bsize); break;
//          case 3: mkstrx_cup(3,bsize); break;
//          case 4: mkstrx_cup(4,bsize); break;
//          default:
//             printf(" mkstrx_cu: lmaxl: %d, lmaxl>4 not supported!\n", lmaxl);
//             leave();
//       }
//    }
//
// //    cudaDeviceSynchronize();
//
//    err = cudaMemcpy(strx  , d_strx  , nstrx *r8, cudaMemcpyDeviceToHost); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaMemcpy(dstrx , d_dstrx , ndstrx*r8, cudaMemcpyDeviceToHost); if (err != cudaSuccess) printf("cudaMemcpy: err: %d, line(%d)\n", err, __LINE__);
//
// //    cudaDeviceSynchronize();
//
//    err = cudaFree(d_ipc   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_lmxl  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_indxcg); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_jcg   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_bas   ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_dlat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_glat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_cg    ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_cy    ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_strx  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_dstrx ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_plat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//    err = cudaFree(d_qlat  ); if (err != cudaSuccess) printf("cudaFree: err: %d, line(%d)\n", err, __LINE__);
//
//
// //    FILE *fl = fopen("strx-cu", "w");
// //    for ( int ib = 0; ib < nbas1; ++ib ) {
// //       for ( int jb = 0; jb < nbas; ++jb ) {
// //          fprintf(fl, "jb, ib: %4d %4d\n", jb, ib);
// //          for (int j = 0; j < nlmq1 ; ++j) {
// //             for (int i = 0; i < nlmq; ++i) {
// //                int idx = (nbas*ib + jb)*nlmq*nlmq1 + nlmq*j+i;
// //                fprintf(fl, " %16.8f", strx[idx]);
// //             }
// //             fprintf(fl, "\n");
// //          }
// //       }
// //    }
// //    fflush(fl);
// //    fclose(fl);
// //    leave();
//
// //    leave();
// }

