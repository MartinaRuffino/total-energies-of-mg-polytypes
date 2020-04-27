#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>

#define LMAX  3
#define LMMAX 16
#define MAXBLOCKSZ 256

__device__
void hstr_cu( bool pv, double &strux, double &dstrx,
               int nlmf, int nlm, int nlmq1, int nlmq,
               double *hl, double *cg, int *indxcg, int *jcg) {
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

   int mlm, klm, lm, lk, lmax, ii, indx, icg1, icg2, llm, lp;
   int const ll[LMMAX] = {
      0,
      1,1,1,
      2,2,2,2,2,
      3,3,3,3,3,3,3
   };

   mlm = threadIdx.y;
   klm = threadIdx.x;
   lm = ll[mlm];
   lk = ll[klm];


//    FILE *fl = fopen("hstrinternal","a");
   lmax = ll[nlmf-1] + ll[nlm-1];

   if (lmax > LMMAX) {
      printf(" change dimensions in hstr\n");
      return;
   }

// --- add together Gaunt-sums ---
   strux = 0.0;
   dstrx = 0.0;
   ii = mlm > klm ? mlm : klm;
   indx = (ii*(ii+1))/2 + (mlm > klm ? klm : mlm);
   icg1 = indxcg[indx];
   icg2 = indxcg[indx+1]-1;
   for (int icg = icg1 - 1; icg < icg2; ++icg) {
      llm = jcg[icg]-1;
      lp = ll[llm];
      if (lm + lk == lp) {
         double t = cg[icg]*hl[llm];
         strux += t;
         if (pv) dstrx -= (lp+1) * t;
      }
   }


//    fclose(fl);
//    __syncthreads();
//    return 0;
}
