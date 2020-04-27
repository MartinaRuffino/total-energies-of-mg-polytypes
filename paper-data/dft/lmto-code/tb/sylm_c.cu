#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>


using std::cout;
using std::endl;

void leave();


// __device__
int mone_on_l(int l) { return (1-2*(l&1));}


// __device__
void besslr(double y, int loka, int lmin, int lmax, double *fi, double *gi) {
//    fi(lmin:lmax),gi(lmin:lmax), fac2l(-nlmax:nlmax*2+3)
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
      fac2l[nlmax] = 1.0;
      for ( l = 1; l<= lmx+1; ++l)
        fac2l[l+nlmax] = fac2l[l-1+nlmax]*(l+l-1);

      for ( l = -1; l >= lmin; --l )
        fac2l[l+nlmax] = fac2l[l+1+nlmax]/(l+l+1);

// C --- Case akap=0 ---
      if (y == 0) {
        for (l = lmin; l <= lmax; ++l) {
          fi[l-lmin] = 1.0/fac2l[l+1+nlmax];
          gi[l-lmin] = fac2l[l+nlmax];
        }
        goto mtogen2;
      }
      my = -y;

// C --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
      tlp1 = lmx+lmx+1;
      dt = 1.0;
      t = 1.0;
      i = 0;
      for ( k = 1; k <= 1000; ++k) {
        if (fabs(dt) < tol) goto conv;
        i = i+2;
        dt2 = i+tlp1;
        dt = dt*my/(i*dt2);
        t = t+dt;
      }
      printf("BESSLR: series not convergent\n");
      return;
   conv:
      dum[0] = t/fac2l[lmx+1+nlmax];

// C --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
      tlp1 =  tlp1-2;
      dt = 1.0;
      t = 1.0;
      i = 0;
      for ( k = 1; k <= 1000; ++k) {
        if (fabs(dt) < tol) goto conv1;
        i = i+2;
        dt2 = i+tlp1;
        dt = dt*my/(i*dt2);
        t = t+dt;
      }
      printf("BESSLR: series not convergent\n");
      return;
   conv1:
      dum[1] = t/fac2l[lmx+nlmax];

// C --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
      ll1 = lmx + lmx + 1;
      ll2 = ll1 + 1;
      nf = ll1;
      for ( k = 3; k <= ll2; ++k) {
        nf = nf-2;
        dum[k-1] = nf*dum[k-1-1] - y*dum[k-2-1];
      }

// C --- Get fi and gi from dum ---
      lmxp1 = lmx+1;
      lmxp2 = lmx+2;
      isn = mone_on_l(lmin);
      for ( k = lmin; k <= lmax; ++k) {
        j1 = lmxp1-k;
        j2 = lmxp2+k;
        fi[k-lmin] = dum[j1-1];
// c   ... n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi[k-lmin] = dum[j2-1]*isn;
        isn = -isn;
      }

// C --- For E<0, use Hankel functions rather than Neumann functions ---
      if (lhank && y < 0.0) {
        srmy = sqrt(-y);
        gi[0-lmin] = 1.0;
        g1 = 1.0+srmy;
        if (lmax >= 1) gi[1-lmin] = g1;
        if (lmax >= 2) {
          tlp1 = 1;
          for ( l = 2; l <= lmax; ++l) {
            tlp1 = tlp1+2;
            gi[l-lmin] = tlp1*gi[l-1-lmin] - y*gi[l-2-lmin];
          }
        }
        if (lmin <= -1) {
          gi[-1-lmin] = (gi[0-lmin] - g1)/y;
          tlp1 = 1;
          if (lmin <= -2) {
            for ( l = -2; l >= lmin; --l) {
              tlp1  = tlp1-2;
              gi[l-lmin] = ((l+l+3)*gi[l+1-lmin] - gi[l+2-lmin])/y;
            }
          }
        }
        exppr = 1.0/exp(srmy);
        for ( l = lmin; l <= lmax; ++l )
          gi[l - lmin] = gi[l-lmin]*exppr;
      }

// C --- Scaling to Andersen's 2nd generation LMTO conventions ---
  mtogen2:
      if (loka == 1) {
        for ( l = lmin; l <= lmax; ++l ) {
         fi[l-lmin] = fi[l-lmin]*fac2l[l+nlmax]*0.5;
         gi[l-lmin] = gi[l-lmin]/fac2l[l+nlmax];
        }
      }
}

// __device__
void ropqln(int m, int l, int n, double *r2, double *z, double cx[3], double *q, int &kk) {
// C- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
// c  Returns kk, which points to the current component of q.
// C  These subroutines are utility routines called by ropyln.f.
// C  Kept separate from ropyln because some optimizing compilers have bugs
// C  (e.g. intel ifort version 11).
// C  These routines are the time-critical steps.
//    q(n,2) r2(n),z(n)
   int mm,i,k2,k1;
   double a,b,xx,yy;

   if (l == m) {
      a = 1.0;
      for ( mm = 0; mm <= m-1; ++mm) a = a*(2*mm+1);
      kk = 0;
      a = a*cx[0];
      for ( i = 0; i < n; ++i ) q[n*kk + i] = a;
   } else if (l == m+1) {
      b = 1.0;
      for ( mm = 0; mm <= m; ++mm) b = b*(2*mm+1);
      b = b*cx[0];
      kk = 1;
      for ( i = 0; i < n; ++i) q[n*kk + i] = b*z[i];
   } else if (l >= m+2) {
      k2 = kk;
      k1 = kk+1;
      if (k1 == 2) k1 = 0;
      xx = -(l+m-1.0)/(l-m)*cx[0]/cx[2];
      yy = (2*l-1.0)/(l-m)*cx[0]/cx[1];
      for ( i = 0; i< n; ++i) q[n*k1 + i] = xx*r2[i]*q[n*k1 + i]+yy*z[i]*q[n*k2 + i];
      kk = k1;
   }

//    ++kk;
}

// __device__
void ropynx(int m, int l, int kk, int n, double *q, double *cm, double *sm, int nd, double *yl) {
      int lav,i;
//       double q(n,2),cm(n),sm(n),yl(nd,1)
//       --kk;
      lav = l*(l+1);
      for ( i = 0; i < n; ++i) yl[nd*(lav+m)+i] = cm[i]*q[n*kk+i];
      if (m == 0) return;
      for ( i = 0; i < n; ++i) yl[nd*(lav-m)+i] = sm[i]*q[n*kk+i];
}

// __device__
void ropcsm(int m, int n, double *x, double *y, double *w, double *cm, double *sm) {
// C- Makes cm and sm. Must be called in sequence m=0,1,2...
   int i;
//    double x(n),y(n),w(n),cm(n),sm(n)

   if (m == 0) {
      for (i = 0; i < n; ++i) cm[i] = 1.0;
      for (i = 0; i < n; ++i) sm[i] = 0.0;
   } else if (m == 1) {
      for (i = 0; i < n; ++i) cm[i] = x[i];
      for (i = 0; i < n; ++i) sm[i] = y[i];
   } else if (m >= 2) {
      for (i = 0; i < n; ++i) w[i] = cm[i];
      for (i = 0; i < n; ++i) cm[i] = x[i]*cm[i] - y[i]*sm[i];
      for (i = 0; i < n; ++i) sm[i] = y[i]*w[i] + x[i]*sm[i];
   }
}

// __device__
void ropyln(int n, double *x, double *y, double *z, int lmax, int nd, double *yl, double *rsq) {
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

   if (n > nd) {
      printf("ropyln: nd too small\n");
      return;
   }
   fpi = 4.0*M_PI;
//       allocate (cm(n),sm(n),q(n*2),h(n))
   double *cm = (double*)malloc(n*sizeof(double));
   double *sm = (double*)malloc(n*sizeof(double));
   double *q  = (double*)malloc(2*n*sizeof(double));
   double *h  = (double*)malloc(n*sizeof(double));

   for (int i = 0; i < n; ++i) rsq[i] = x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
   for (int i = 0; i < 3; ++i) cx[i] = 0.0;

// C --- Loop over m: cx are normalizations for l, l-1, l-2 ---
   f2m = 1.0;
   for (m = 0; m <= lmax; ++m) {
      ropcsm(m,n,x,y,h,cm,sm);
//       printf("ropcsm: m,n,x,y,h,cm,sm: %d %d %16.8f %16.8f %16.8f %16.8f %16.8f\n",
//          m,n,x[0],y[0],h[0],cm[0],sm[0]);
      if (m == 0) {
         cx[0] = sqrt(1/fpi);
      } else {
         f2m = f2m*(2*m*(2*m-1));
         cx[0] = sqrt(((2*m+1)*2)/(fpi*f2m));
      }
      for ( l = m; l <= lmax; ++l) {
         ropqln(m,l,n,rsq,z,cx,q,kk);
         ropynx(m,l,kk,n,q,cm,sm,nd,yl);
         cx[2] = cx[1];
         cx[1] = cx[0];
         cx[0] = cx[0]*sqrt(double((l+1-m)*(2*l+3))/double((l+1+m)*(2*l+1)));
      }
   }
   free(cm);
   free(sm);
   free(q);
   free(h);
}

// __device__
void soldhj(double r[3], double e, int loka, int lmax, double *hl){
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
   int ilm,l,m;
   double rfac,r2;
   double phi[11],psi[11];

//    double dl[100];

// C     call sylm(r,hl,lmax,r2)
   ropyln(1,&r[0],&r[1],&r[2],lmax,1,hl,&r2);
//    printf("\nropyln: r: %16.8f %16.8f %16.8f %16.8f",r[0],r[1],r[2],r2);
//    for (int i = 0; i<100; ++i) printf("%s %16.8f", hl[i],((i%5)==0?"\n":""));

   besslr(e*r2,loka,0,lmax,phi,psi);
//    printf("\npsi:");
//    for (int i = 0; i<11; ++i) printf(" %16.8f", psi[i]);
//    printf("\n");

   ilm = 0;
   rfac = sqrt(r2);
   if (r2 < 1e-10) r2 = 1.0;
   for ( l = 0; l <= lmax; ++l) {
// C       rfac = 1/r**(2l+1), or 1/(r/w)**(2l+1) using OKA conventions
      rfac = rfac/r2;
      for ( m = -l; m <= l; ++m ) {
// C       xx = cy(ilm)*hl(ilm)
// C       bl(ilm) = phi(l)*xx
// C       hl(ilm) = (rfac*psi(l))*xx
//          bl[ilm] = phi[l]*hl[ilm]
         hl[ilm] = (rfac*psi[l])*hl[ilm];
         ++ilm;
      }
   }


//    int tid = threadIdx.y*blockDim.x + threadIdx.x;
//    int nths = blockDim.x*blockDim.y;

//    for (int i = tid; i < 100; i += nths) hl[i] = dl[i];

}


double inclmod1_c(double a) {
// ! Alternative implementation. Should give the same result
// !       s = sign(1.0_dp, a)
// !       r = -s*(modulo(-s*a+0.5_dp, 1.0_dp) - 0.5_dp)
// ! plot it if not clear

// alternative: a - rint(a). rint is implementation and settings dependent...

// this whole block finds the nearest integer with the bounday case of 0.5 rounded towards 0;
   double s = copysign(1.0, a);
   double t = a + 0.5*s;
   double r = trunc(t);
   if (t == r) r = r - s;

   return a - r;
}


void ploughmans_3x1_mv_c(double *m, double u[3], double v[3]) {
//    follows conventional mm in F so the access looks improper here in C
   for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) v[j] += m[3*i+j]*u[i];
}


void shorten_c(double p[3], double *plat, double *ilat) {

   double r[3];
   for (int i = 0; i < 3; ++i) r[i] = 0.0;
   ploughmans_3x1_mv_c(ilat, p, r);

   for (int i = 0; i < 3; ++i) r[i] = inclmod1_c(r[i]);

   for (int i = 0; i < 3; ++i) p[i] = 0.0;
   ploughmans_3x1_mv_c(plat, r, p);
}


// static clock_t sylm_c_tr = 0, sylm_c_tk = 0;

double sylm_c(double r[3], double *yl, int lmx) {
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

   int const lmaxx = 3;
   int const nxx = (lmaxx+1)*(lmaxx+1);

   int lav,lp1,lm1,n,nt;
   double r2,st,z2, r2s;
   double c[nxx], s[nxx], p[nxx][nxx];
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
            p[l][m] = (r2*(l+m)*p[lm1][m]-st*p[l][m+1])/(z*(l-m));
            yl[lav+m] = p[l][m]*c[m];
            yl[lav-m] = p[l][m]*s[m];
         }
      }
   }
   return r2s;
}


void rcsl01_c(double tau[3], double a, int lmax, double alat, double *rlat, int nkr, double vol, double *dl) {
//  k-space part of reduced structure constants for e=0 and q=0
   double r2;

   int ilm,nlm;
   int const lmaxx = 3;
   int const n = (lmaxx+1)*(lmaxx+1);
   double fpibv,gamma,scalp,tpiba,yyy,eiphi[2];
   double r[3], yl[n];
   double const tpi = 2*M_PI;


   if (lmax > lmaxx) {
      printf("rcnsl0: increase lmaxx\n");
      return;
   }

   gamma = 0.25/(a*a);
   fpibv = 2.0*tpi/vol;
   tpiba = tpi/alat;
   nlm = (lmax+1)*(lmax+1);
   for (int ilm = 0; ilm < nlm; ++ilm) dl[ilm] = 0.0;
   dl[0] = -fpibv*gamma;
   for ( int ir = 1; ir < nkr;  ++ir ) {
      for (int i = 0; i < 3; ++i) r[i] = tpiba*rlat[3*ir + i];
//       printf(" rlat: %20.12f %20.12f %20.12f\n",rlat[3*ir + 0], rlat[3*ir + 1], rlat[3*ir + 2]);
//       printf("r: %20.12f %20.12f %20.12f\n",r[0],r[1],r[2]);
      scalp = 0.0; for (int i = 0; i < 3; ++i) scalp += r[i]*tau[i];
      scalp *= alat;
      eiphi[0] = cos(scalp);
      eiphi[1] = sin(scalp);
//       printf("eiphi: %20.12f %20.12f\n",eiphi[0],eiphi[1]);
//       clock_t t1 = clock();
      r2 = sylm_c(r,yl,lmax);
//       sylm_c_tk += clock() - t1;
//       printf("  r2: %20.12f\n",r2);
      yyy = fpibv*exp(-gamma*r2)/r2;
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





void rcsl02_c(double tau[3], double a, int lmax, double alat, double *dlat, int nkd, double *dl) {
//  real space summation
   double r2;
   int ilm,ir1;
   int const lmaxx = 3;
   int const n = (lmaxx+1)*(lmaxx+1);
   double a2,cc,gl,r1,ta2,dl0;
   double r[3],yl[n],chi[lmaxx+1];
   double const twoinvsqpi = M_2_SQRTPI; //1.12837916709551257390_8 ! 2/sqrt(pi)

   ir1 = ((tau[0]*tau[0]+tau[1]*tau[1]+tau[2]*tau[2]) > 1e-6) ? 0 : 1;
//       if (sum(tau*tau) > 1d-6) ir1=1

   if (lmax > 0) {
      a2 = a*a;
      ta2 = 2*a2;
      cc = ta2*a*twoinvsqpi;
      for ( int ir = ir1; ir < nkd; ++ir ) {
         for (int i = 0; i < 3; ++i) r[i] = alat*(tau[i]-dlat[3*ir + i]);

//          clock_t t1 = clock();
         r2 = sylm_c(r,yl,lmax);
//          sylm_c_tr += clock() - t1;

         r1 = sqrt(r2);
         chi[0] = erfc(a*r1)/r1;
         gl = -cc*exp(-a2*r2)/ta2;
         for ( int l = 0; l < lmax; ++l ) {
            chi[l+1] = ((2*l+1)*chi[l] - gl)/r2;
            gl *= ta2;
         }
         ilm = 0;
         for ( int l = 0; l <= lmax; ++l ) {
            for ( int m = 1; m <= 2*l+1; ++m ) {
               dl[ilm] += yl[ilm]*chi[l];
               ++ilm;
            }
         }
      }
// ...In case lmax = 0 do everything explicitly
   } else {
      dl0 = 0.0;
      for ( int ir = ir1; ir < nkd; ++ir) {
         for (int i = 0; i < 3; ++i) r[i] = tau[i]-dlat[3*ir + i];
         r2 = 0.0; for (int i = 0; i < 3; ++i) r2 += r[i]*r[i];
         r1 = alat*sqrt(r2);
         dl0 += erfc(a*r1)/r1;
      }
      dl[0] += dl0;
   }

// --- add dl3 for diagonal sructure constants ------
   if (ir1 == 1) dl[0] -= a*twoinvsqpi;
}


static clock_t rcsl01_c_t = 0, rcsl02_c_t = 0, rcsl03_c_t = 0;

// __global__

void rcnsl0_c(double tau[3], double a, int lmax, double alat,
                                     double *rlat, int nkr,
                                     double *dlat, int nkd,
                                     double vol, double *cy, double *dl) {
//  reduced structure constants on lattice for e=0 and q=0.
//  result is periodic sum of (2*l-1)!!*ylm/r**(l+1). additive
//  constant for l=0 is chosen so that function averages to zero.


   clock_t t1 = clock();
   rcsl01_c(tau,a,lmax,alat,rlat,nkr,vol,dl);
//    printf("rc01: %22.12f %22.12f %22.12f\n", dl[0], dl[1], dl[2]);
   clock_t t2 = clock();
   rcsl01_c_t += t2 - t1;
   rcsl02_c(tau,a,lmax,alat,dlat,nkd,dl);
//    printf("rc02: %22.12f %22.12f %22.12f\n", dl[0], dl[1], dl[2]);
   clock_t t3 = clock();
   rcsl02_c_t += t3 - t2;

   int nlm = (lmax+1)*(lmax+1);
   for (int ilm = 0; ilm < nlm; ++ilm) dl[ilm] *= cy[ilm];

   rcsl03_c_t += clock() - t3;


}

double one_on_l(int l) {
   return double(1-2*(l&1));
}






void hstr_c(bool mol, bool pv, int ldip,
                        double *strx, double *drstrx,
                        int nlmf, int nlm, int nlmq1, int nlmq,
                        double *hl, double *cg, int *indxcg, int *jcg,
                        double vol) {
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

   int icg1,icg2,indx,lk,llm,lm,lp,lmax;
   int const lmxx = 16;
   double sig[lmxx+1],fpibv;
   double const fourpi = 4.0*M_PI;


   int const ll[100] = {
      0,
      1,1,1,
      2,2,2,2,2,
      3,3,3,3,3,3,3,
      4,4,4,4,4,4,4,4,4,
      5,5,5,5,5,5,5,5,5,5,5,
      6,6,6,6,6,6,6,6,6,6,6,6,6,
      7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
      8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
      9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9
   };

   int const msm[8] = {2,3,1,4,5,8,6,7};

//    df(l) = (2l+1)!!
   int const df[6] = {1,3,15,105,945,10395};

//    FILE *fl = fopen("hstrinternal","a");
   lmax = ll[nlmf-1] + ll[nlm-1];

   if (lmax > lmxx) {
      printf(" change dimensions in hstr lmxx,lmax: %d, %d\n", lmxx, lmax);
//       leave();
      exit(1);
   }
// --- (-1)^l ---
   sig[0] = 1.0;
   if (lmax > 0) {
      for (int l = 0; l < lmax; ++l) sig[l+1] = -sig[l];
   }
// --- add together Gaunt-sums ---
   for (int mlm = 0; mlm < nlmf; ++mlm) {
      lm = ll[mlm];
      for (int klm = 0; klm < nlm; ++klm) {
         lk = ll[klm];
         double cghl = 0.0;
         double cghlr = 0.0;
         int ii = mlm > klm ? mlm : klm;
         indx = (ii*(ii+1))/2 + (mlm > klm ? klm : mlm);
         icg1 = indxcg[indx];
         icg2 = indxcg[indx+1]-1;
         for (int icg = icg1 - 1; icg < icg2; ++icg) {
            llm = jcg[icg]-1;
            lp = ll[llm];
            if (lm + lk == lp) {
               double t = cg[icg]*hl[llm];
//     fprintf( fl,
//  "mlm,klm,lm,lk,icg,llm,cg[icg],hl[llm],t: %3d %3d %3d %3d %3d %3d %16.8f %16.8f %16.8f\n",
//   mlm,klm,lm,lk,icg,llm,cg[icg],hl[llm],t);
               cghl += t;
               if (pv) cghlr -= (lp+1) * t;
            }
         }
         double ff = fourpi * sig[lk] * (sqrt(double((2*lk+1)*(2*lm+1)))/double(df[lk]*df[lm]));

         int mslm = (mlm > 0 && mlm < 9) ? msm[mlm-1] : mlm;
         int mslk = (klm > 0 && klm < 9) ? msm[klm-1] : klm;

         strx[mslm * nlmq + mslk] = cghl * ff;
// need to cut off l+1 terms for B'. Quick fix for now
         if (pv && (mlm < nlmq)) drstrx[mslm*nlmq + mslk] = cghlr * ff;
      }
   }
// --- the following includes extra p terms 'implicitly' ---
   if ((nlm > 1) && ((ldip == 2) || (ldip == 3)) && (!mol)) {
      fpibv = fourpi/vol;
//           do  ilm = 2, min0(4,nlm,nlmf)
      for (int ilm = 1; ilm < 4; ++ilm) {
         int sl = ll[ilm];
         double ff = fpibv * sqrt(double((2*sl+1)*(2*sl+1)))/double(df[sl]*df[sl]);
         int msli = (ilm > 0 && ilm < 10)? msm[ilm-1] : ilm;
         strx[msli * nlmq + msli] -= ff;
         if (pv) drstrx[msli * nlmq + msli] += 3.0*ff;
      }

   }
//    fclose(fl);

//    __syncthreads();
//    return 0;
}


extern "C"
void mkstrx_c(int ldip, int nbas, int nbas1, int ib0, int nlmq1, int nlmq,
                  double *bas, int *ipc, int nclas, int *lmxl, double awld, double alat, double vol,
                  double *dlat, int nkd, double *glat, int nkg,
                  int *indxcg, int *jcg, double *cg, double *cy,
                  int fpv, int fmol, double *strx, double *dstrx,
                  double *plat, double *qlat) {

// strx(nlmq,nlmq1,nbas,nbas1),dstrx(nlmq,nlmq,nbas,nbas1)

//C- Make lookup table of strux and radial derivatives
//C ----------------------------------------------------------------------
//Ci Inputs:
//Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
//Ci        : 2 include 'slab' dipole correction to Ewald
//Ci        : any other number - skip the dipole correction
//Ci   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
//Ci   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
//Ci   nlmq : max L-cutoff for multipoles, leading dimensions of dstrx,
//Ci           second dimension of strx
//Ci   lmxl : l-cutoff for multipoles for each class of species
//Ci   pv   : true if radial derivatives required
//Ci   MOL  : true for molecule (cluster) branch
//Co Outputs:
//Co  strx  : coefficients B_LL'(R'-R) (structure constants)
//Co  dstrx : derivatives of structure constants (x dB_LL'(x)/dx) at x = |R'-R|
//Co          if pv = F dstrx is not touched (see hstra.f)
//Cr Remarks
//Cr   Instead of making B on each sc iteration, the program prepares the whole B
//Cr   in the beginning of the self-consistency cycle to be used as a lookup table.
//Cr   This results in faster calculation at an expence of memory.
//Cr
//Cr   Calling mkstrx is set as the default option. To call instead hstr on each
//Cr   iteration (and therefore to reduce memory but increase computational time),
//Cr   use switch --sfly
//Cr
//Cr   Efficiency issues: there are two symmetry properties of zero energy
//Cr   structure constants that are used in mkstrx:
//Cr     B_LL'(R'-R) = B_L'L(R-R')                (1)
//Cr     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
//Cr   same properties hold for the radial derivative of B.
//Cr
//Cr   According to (1), mkstrx calculates strx(:,:,ib,jb) only for the
//Cr   lower triangle {ib = 1:nbas, jb = 1:ib}
//Cr   property (2) is used in hstra in the same way, ie strx(ilm,jlm,:,:)
//Cr   is sought only for the triangle {ilm = 1:Lmax, jlm=1:ilm}.
//Cr   hstra should therefore be called with Lmax >= Lmax', otherwise it stops.
//Cr
//Cr   mkstrx fills only lower triangles of strx and dstrx in the above sense;
//Cr   properties (1) and (2) are used later, when strx and dstrx are copied into
//Cr   'local' arrays fstrx and fdstrx in tbesel.f We decided to keep it this way
//Cr   because in future we plan to reduce the size of strx and dstrx by cutting of
//Cr   the unused parts of arrays strx and dstrx.
//Cr
//Cr   hstra vs hstr: a call to hstr (a routine that actually makes strux) is
//Cr   replaced with a call to hstra, an 'optimised' version of hstr. The difference
//Cr   between hstr and hstra is that hstra computes only the lower triangle (see above)
//Cr   of strux but is restricted to Lmax >= Lmax'. hstr does not have any restrictions
//Cr   and is called if the --sfly switch is on.
//Cr
//Cb Bugs
//Cb   mkstrx is not written in the most efficient way and could be further refined wrt to
//Cb   amount of calculated elements of strux. At the moment mkstrx is a result of certain
//Cb   trade off between performance and clarity of the program.
//Cb
//Cu Updates
//Cu    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
//Cu                      jb<=ib only)
//Cu    19 Jan 2010 (ATP) first written
//C ----------------------------------------------------------------------



   bool mol = fmol == 1;
   bool pv = fpv == 1;

   clock_t t[16], td;

   for (int i = 0; i < 16; ++i) t[i] = 0;

   t[0] = clock();

//    mol = true;

//    memset( strx, 0, nlmq*nlmq1*nbas*nbas1);
//    memset(dstrx, 0, nlmq*nlmq *nbas*nbas1);

   int ibr, li, lj, lmax, nlx, nlf, lmxst; //lmxf
   double tau[3], hl[100];
//    double *blk = (double *) malloc(nlmq1*nlmq*sizeof(double));
//    double *dblk = (double *) malloc(nlmq*nlmq*sizeof(double));

   for ( int ib = 0; ib < nbas1; ++ib ) {
      ibr = ib + ib0;
      li = lmxl[ipc[ibr]-1];
//c       nlmi = (li+1)**2
//C...  nlmq1 == nlmq if no force or pressure are required.
//C     Otherwise need to go up to lmax+1 in the first index of B_LL'

//    int li1 =  (nlmq1 > nlmq) ? li + 1 : li;

      for ( int jb = 0; jb < nbas; ++jb ) {
         lj = lmxl[ipc[jb]-1];

//C...  The lines below are because we shall be filling the jb > ib triangle
//C     using symmetry properties of B (see tbesel.f) AND because we need l+1
//C     for forces. No complications if forces are not required.
         if (nlmq1 > nlmq) {
            lmax = max(li,lj);
//             lmxf = lmax+1;
            nlx = (lmax+1)*(lmax+1);
            nlf = (lmax+2)*(lmax+2);
            lmxst = 2*lmax+1;
         } else {
            nlx = (lj+1)*(lj+1);
            nlf = (li+1)*(li+1);
            nlf = max(nlf,nlx);
            lmxst = li + lj;
         }

         for (int i = 0; i < 3; ++i) tau[i] = bas[3*jb+i] - bas[3*ibr+i];

//          printf("MOL BRANCH comes: %s\n", (mol?"yes":"no"));

         if (!mol) {
            clock_t tl = clock();
            shorten_c(tau,plat,qlat);
            clock_t tl1 = clock();
            t[1] += tl1 - tl;
            rcnsl0_c(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl);
            t[2] += clock() - tl1;
         } else {
//             printf("MOL BRANCH\n");
            for (int i = 0; i < 3; ++i) tau[i] *= alat;
//             printf("\nc ib,jb: %4d %4d",ib,jb);
            soldhj(tau,0.0,0,lmxst,hl);
//             printf("\nhl:");
//             for (int i = 0; i<100; ++i) printf("%s %16.8f", hl[i],((i%5)==0?"\n":""));

         }

//          exit(0);

//          for (int j = 0; j < nlmq ; ++j)
//             for (int i = 0; i < nlmq1; ++i)
//                blk[nlmq1*j+i] = 0.0;
//
//          for (int j = 0; j < nlmq ; ++j)
//             for (int i = 0; i < nlmq; ++i)
//                dblk[nlmq1*j+i] = 0.0;

         clock_t tl = clock();
         hstr_c(mol, pv, ldip, &strx[(nbas*ib + jb)*nlmq*nlmq1],
               &dstrx[(nbas*ib + jb)*nlmq*nlmq],nlf,nlx,nlmq1,nlmq,hl,cg,indxcg,jcg,vol);
         t[3] += clock() - tl;
//          for (int j = 0; j < nlmq1 ; ++j)
//             for (int i = 0; i < nlmq; ++i)
//                strx[(nbas*ib + jb)*nlmq*nlmq1 + nlmq*j+i] = blk[nlmq1*i+j];
//
//          if (pv)
//             for (int j = 0; j < nlmq ; ++j)
//                for (int i = 0; i < nlmq; ++i)
//                   dstrx[(nbas*ib + jb)*nlmq*nlmq + nlmq*j+i] = dblk[nlmq*i+j];
//          if (jb > 2) exit(5);
      }
   }

   t[4] = clock();

   td = t[4] - t[0]; printf("mkstrx_c: %12.8fs\n", double(td)/double(CLOCKS_PER_SEC));
   printf("shorten_c: %12.8fs    \n", double(t[1])/double(CLOCKS_PER_SEC));
   printf("rcnsl0_c: %12.8fs     \n", double(t[2])/double(CLOCKS_PER_SEC));
   printf("   rcsl01_c: %12.8fs  \n", double(rcsl01_c_t)/double(CLOCKS_PER_SEC));
//    printf("        k Y: %12.8fs  \n", double(sylm_c_tk)/double(CLOCKS_PER_SEC));
   printf("   rcsl02_c: %12.8fs  \n", double(rcsl02_c_t)/double(CLOCKS_PER_SEC));
//    printf("        r Y: %12.8fs  \n", double(sylm_c_tr)/double(CLOCKS_PER_SEC));
   printf("   elmul: %12.8fs     \n", double(rcsl03_c_t)/double(CLOCKS_PER_SEC));
   printf("hstr_c: %12.8fs       \n", double(t[3])/double(CLOCKS_PER_SEC));


//    FILE *fl = fopen("strx-c", "w");
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
//    exit(0);

//    free(blk);
//    free(dblk);
}

