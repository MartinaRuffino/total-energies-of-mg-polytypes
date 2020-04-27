      subroutine besslz(y,iopt,lmin,lmax,fi,gi)
C- Radial part of Hankel and Bessel functions, complex energy
C ----------------------------------------------------------------------
Ci Inputs:
Ci   y     :e*r**2 = (kappa*r)**2 = x**2
Ci   lmin  :minimum l
Ci   lmax  :maximum l
Ci   iopt  :1s digit:
Ci          1 scale fi,gi to Andersen's definitions (see below)
Ci         :10s digit:
Ci          1 return -1* Neumann function in gi, rather than Hankel
Co Outputs:
Co   fi    :   j_l(x)/x^l
Co   gi    : i*h_l(x)*x^(l+1) if 10s digit of iopt = 0
Co   gi    :-1*n_l(x)*x^(l+1) if 10s digit of iopt = 1
Cr Remarks
Cr   i h_0 = exp(ikr)/ikr, j_0 = sin(kr)/kr, k = sqrt(e), Im k >= 0.
Cr
Cr   Special cases:
Cr     fi(y,l=0) = sin(kappa*r)/(kappa*r)
Cr     gi(y,l=0) = exp(i*kappa*r)
Cr     fi(y,l) -> 1/(2l+1)!!     for y->0
Cr     gi(y,l) -> (2l-1)!!       for y->0
Cr
Cr   For Andersen's definitions, fi and gi are scaled as follows:
Cr     fi(OKA)=  fi*(2l-1)!!/2
Cr     gi(OKA)=  gi/(2l-1)!!
Cr
Cr   Relations between Hankel h_l, Bessel j_l and Neumann functions n_l:
Cr   are defined differently from Gunnarsson (PRB 27, 7144, 1983):
Cr      i h_l =  n_l + i j_l  = i h^(2)_l <- Gunnarsson, eq. A8
Cr      i h_l = -n_l + i j_l  = i h^(1)_l <- here
Cr
Cr   Methfessel defines Hankels as follows:
Cr      h_0(kappa,r) = exp(i kappa r)/r  h_1 = h_0 (1 - i kappa r) / r^2
Cr      H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)
Cr
Cr   Method of computation:
Cr   For Im k >0, Hankel functions h_0 and h_1 = (1 - ikr) h_0 are
Cr   directly calculated and h_l by recursion for higher l.
Cr
Cr   Bessel functions are calculated for lmax and lmax-1 by a expansion
Cr   in powers of x^2=y:
Cr
Cr                      (-x^2/2)^k
Cr        fi =  Sum_k  --------------  = dum(lmx+1-l)
Cr                     k! (2l+1+2k)!!
Cr
Cr        gi = dum(lmx+2+l), y>0
Cr
Cr   The remaining j are obtained by backward recursion:
Cr
Cr      j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x)    l=lmx-2,lmx-3,..,-lmx-1
Cr
Cr  <==> dum(k) = (2*lmx+5-2*k)*dum(k-1)-y*dum(k-2)   k=3,4,...,2*lmx+2
Cr
Cr   and the Neumann function are given through:
Cr
Cr         n_l(x) = j_{-l-1}*(-1)^{l+1}
Cr
Cr   and finally the hankel functions by  i h_l(x) = -n_l(x) + i j_l(x)
Cr
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iopt
      integer lmin,lmax
      double complex y,fi(lmin:lmax),gi(lmin:lmax)
C Local variables:
      logical lneu
      integer i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,lgunit,nf,tlp1,ll1,
     .        ll2,nlmax
      parameter (nlmax=20)
      double complex my,t,dt,dum(nlmax*4+2),exppr,srmy,ikr
      double precision dt2,tol,fac2l(-nlmax:nlmax*2+3)
      parameter(tol=1.d-15)
C External calls:
      external lgunit
C Intrinsic functions:
      intrinsic dabs,dexp,dsqrt,max0

      lmx = max0(lmax,2)
      if (lmx > nlmax+nlmax)
     .  call rxi(' BESSLZ : lmax gt nlmax*2, lmax=',lmx)
      lneu = iopt/10 /= 0

C --- A table of fac2l(l)=(2l-1)!!
      fac2l(0) = 1.d0
      do  10  l = 1, lmx+1
   10 fac2l(l) = fac2l(l-1) * (l+l-1)
      do  11  l = -1, lmin, -1
   11 fac2l(l) = fac2l(l+1) / (l+l+1)

C --- Case kappa=0 ---
      if (y == (0d0,0d0)) then
        do  12  l = lmin, lmx
          fi(l) = 1/fac2l(l+1)
          gi(l) = fac2l(l)
   12   continue
        goto 100
      endif
      my = -y

C --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
      tlp1 = lmx+lmx+1
      dt = 1
      t = 1
      i = 0
      do  20  k = 1, 1000
        if (abs(dt) < tol) goto 21
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
   20 continue
      call rx('BESSLZ: series not convergent')
   21 continue
      dum(1) = t/fac2l(lmx+1)

C --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
      tlp1 =  tlp1-2
      dt = 1
      t = 1
      i = 0
      do  30  k = 1, 1000
        if (abs(dt) < tol) goto 31
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
   30 continue
      call rx('BESSLZ: series not convergent')
   31 continue
      dum(2) = t/fac2l(lmx)

C --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
      ll1 = lmx + lmx + 1
      ll2 = ll1 + 1
      nf = ll1
      do  40  k = 3, ll2
        nf = nf-2
        dum(k) = nf*dum(k-1) - y*dum(k-2)
   40 continue

C ... Choose root for which Im(kap)>0
      srmy = sqrt(y)
      if (dimag(srmy) < 0d0) srmy = -srmy
      ikr = (0d0,1d0)*srmy

C --- Get fi and gi from dum ---
      lmxp1 = lmx+1
      lmxp2 = lmx+2
      isn = (-1)**lmin
      dt = y**lmin
      do  50  k = lmin, lmx
        j1 = lmxp1-k
        j2 = lmxp2+k
        fi(k) = dum(j1)
C   ... n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi(k) = dum(j2)*isn
C   ... i h_l(x) = -n_l(x) + i j_l(x) => gi += ikr e^l fi
        if (.not. lneu) gi(k) = gi(k) + fi(k)*ikr*dt
        isn = -isn
        dt = dt*y
   50 continue

C --- Case Im(kap) > 0 ... should be more stable numerically  ---
      if (dble(ikr) < 0d0) then
        gi(0) = 1
        gi(1) = 1 - ikr
        if (lmx >= 2) then
          tlp1 = 1
          do  62  l = 2, lmx
            tlp1 = tlp1+2
            gi(l) = tlp1*gi(l-1) - y*gi(l-2)
   62     continue
        endif
        tlp1 = 3
        if (lmin <= -1) then
          do  64  l = -1, lmin, -1
            tlp1  = tlp1-2
            gi(l) = ((l+l+3)*gi(l+1) - gi(l+2))/y
   64     continue
        endif
        exppr = exp(ikr)
        do  66  l = lmin, lmx
          gi(l) = gi(l)*exppr
   66   continue
C   ... i h_l(x) = -n_l(x) + i j_l(x) => gi -= (ikr) e^l fi
        if (lneu) then
          isn = (-1)**lmin
          dt = y**lmin
          do  68  l = lmin, lmx
            gi(l) = gi(l) - fi(l)*ikr*dt
            dt = dt*y
   68     continue
        endif
      endif

C --- Scaling to Andersen's conventions ---
  100 continue
      if (mod(iopt,10) == 0) return
      do  70  l = lmin, lmx
        fi(l) = fi(l)*fac2l(l)*0.5d0
        gi(l) = gi(l)/fac2l(l)
   70 continue
      end
