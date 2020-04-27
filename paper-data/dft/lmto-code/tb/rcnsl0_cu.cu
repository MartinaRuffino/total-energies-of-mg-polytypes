#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>


#define LMAX  3
#define LMMAX 16
#define MAXBLOCKSZ 256



extern __device__ double sylm_cu(double *r, double *yl, int lmx);

__device__
void rcsl01_cu(double tau[3], double a, int lmax, double alat, double *rlat, int nkr, double vol, double *dl) {
//  k-space part of reduced structure constants for e=0 and q=0
   double r2;

   int ilm,nlm;
   double fpibv,gamma,scalp,tpiba,eiphi[2];
   double r[3], yl[LMMAX];
   double const tpi = 2*M_PI;

   int tid = threadIdx.y*blockDim.x + threadIdx.x;
   int tot = blockDim.x * blockDim.y;


   if (lmax > LMAX) {
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
      r2 = sylm_cu(r,yl,lmax);
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




__device__
void rcsl02_cu(double tau[3], double a, int lmax, double alat, double *dlat, int nkd, double *dl) {
//  real space summation
   double r2;
   int ilm,ir1;
   double a2,cc,gl,r1,ta2,dl0;
   double r[3],yl[LMMAX],chi;//chi[LMAX+1];
   double const twoinvsqpi = M_2_SQRTPI; //1.12837916709551257390_8 ! 2/sqrt(pi)

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

         r2 = sylm_cu(r,yl,lmax);

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


