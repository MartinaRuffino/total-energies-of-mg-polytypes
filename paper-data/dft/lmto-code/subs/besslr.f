      subroutine besslr(y,loka,lmin,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C ----------------------------------------------------------------------
Ci Inputs:
Ci   y     :y = e*r**2 = z**2 = (kappa*r)**2 = -(akap*r)**2, Im kappa>=0
Ci   loka  :1s digit
Ci         :0 Methfessel's conventions
Ci         :1 Andersen conventions from 2nd generation LMTO
Ci         :2 the following convention:
Ci         :  Anderson's convention for gi so that
Ci         :  Anderson's convention for fi * 2 * (2*l+1)
Ci         :  This definition has the property
Ci         :  phi -> 1 as y->0
Ci         :  gi  -> 1 as y->0
Ci         :10s digit
Ci         : 0 always use besslr
Ci         : 1 call besnu for y>=90
Ci         : 2 always call dbesnu
Ci   lmin  :minimum l
Ci   lmax  :maximum l
Co Outputs:
Co   fi    :proportional to (Bessel function) / r**l.  Constant
Co         :of proportionality depends on conventions; see Remarks.
Co         :fi is evaluated in a power series expansion.
Co         :fi(l=0) = sinh(x)/x, x = sqrt(-y) for y<=0
Co   gi    :proportional to :
Co         :Neumann function * r**(l+1) for y>0
Co         :Hankel function  * r**(l+1) for y<0.
Co         :See Remarks for constant of proportionality.
Co         :For y>0 gi is evaluated in a power series expansion.
Co         :gi(l=0) = cos(sqrt(y))
Co         :For y<0:
Co         :gi(l=0) = exp(-x) = exp(i kappa r)  where x=sqrt(-y)
Co         :gi(l=1) = (1+x)*exp(-x) = (1-i kappa r)gi(l=0)
Co         :gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
Cr Remarks:
Cr  *Bessel and Hankel functions from fi and gi.
Cr   Let j, n, h be spherical Bessel, Neumann and Hankel functions.
Cr   We use Im kappa > 0.
Cr   Conventions for Hankel functions vary somewhat and are defined below.
Cr   See radhjz for generation of j,n,h following different conventions.
Cr
Cr  *Jackson, and Morse and Feshback conventions: z = kappa * r
Cr     h = j + i*n,  h = Hankel function of the first kind h^1, and
Cr     h^1,2(z) = j(z) +/- i n(z)
Cr     h^1_0(z) = exp(iz)/iz
Cr     h^1_1(z) = exp(iz)/iz * (-1 - i/z)
Cr   i h^1_l(z) = gi_l(z**2) / z^(l+1) (computed if y<0)
Cr              ->(2l-1)!! / z^(l+1) for z -> 0
Cr       j_0(z) = sin(z)/z
Cr       j_l(z) = fi_l(z**2) * z^l
Cr       n_0(z) = -cos(z)/z
Cr     Limiting cases for z = kappa*r -> 0:
Cr       n_l(z) ->  -(2l-1)!!/z^(l+1) ( 1 - z^2/(2(2l-1) + ...)
Cr       j_l(z) ->   z^l / (2l+1)!!   ( 1 - z^2/(2(2l+3) + ...)
Cr       h^1_l(z) -> -i (2l-1)!!/z^(l+1)
Cr     Limiting cases for z = kappa*r -> infty:
Cr       h^1_l(z) ->  exp(i(z-l*pi/2))/iz
Cr       n_l(z)   -> -cos(z-l*pi/2)/z
Cr       j_l(z)   ->  sin(z-l*pi/2)/z
Cr
Cr  *Gunnarsson's conventions (PRB 27, 7144, 1983).
Cr   Somewhat confusing.  Gunnarsson chooses  k = -kappa, Im k <= 0.
Cr   Apparently his definitions correspond to -1 * standard h^1.
Cr     h_l = j_l(k r) - i n_l(k r)     see eqn prior to eq. A8
Cr   Taking his A15 as a definition deriving from Andersen, we have
Cr   i h_l = -gi(z**2) / (kr)^(l+1)
Cr         -> - (2l-1)!! / (kr)^(l+1) for k->0  A10 (sign error fixed)
Cr     h_0 = -exp(-ikr)/ikr
Cr     h_1 = h_0 (1 - ikr) / kr
Cr     j_0 = sin(kr)/kr
Cr     Functions of odd order are imaginary for e<0
Cr     h_l(kr) = h^1_l(-kr)
Cr             = i gi_l((kr)**2) / (kr)^(l+1)
Cr     j_l(kr) = sin(kr)/kr
Cr     Limiting cases for energy->0:
Cr     h_l -> i (2l-1)!!/(kr)^(l+1)  j_l -> (kr)^l / (2l+1)!!
Cr     Wronskian: {h,j} = -1 / i kap
Cr
Cr  *Methfessel's conventions (label as an,ak,aj) have the properties:
Cr     (1)  functions ak and aj, and an are real for e<=0
Cr     (2)  cases e /= 0 and e == 0 same equations
Cr     Define akap = -i kappa = sqrt(-e):  akap is real and >0 for e<0.
Cr       ak_0 = exp(-akap r)/r;  ak_1 = (1 + akap*r) ak_0 / r
Cr     Relation to standard conventions:
Cr       ak_l = -i kappa^(l+1) h_l = -i (-i akap)^(l+1) h_l
Cr            = gi_l / r^(l+1)
Cr       aj_l = 1/kappa^l j_l      = 1/(-i akap)^l j_l
Cr            = fi_l r^l
Cr     Define Neumann function as:
Cr       an = ah - i kap e^l aj = ah + akap e^l aj
Cr     Solid Hankels, Bessels are defined as (CHECK)
Cr       H_L = Y_L(-grad) ak(r);   J_L = E^(-l) Y_L (-grad) aj(r)
Cr     Limiting cases for energy->0:
Cr       ak -> (2l-1)!!/r^(l+1)      aj -> r^l/(2l+1)!!
Cr
Cr  *Andersen conventions (label as K,J)
Cr       K_l = -1/(2l-1)!! (kappa avw)^(l+1) i h_l(kr)
Cr           =  1/(2l-1)!! (avw)^(l+1) ak_l
Cr           =  gi(OKA) / (r/w)^(l+1)
Cr       J_l = 1/2*(2l-1)!! (kappa avw)^(-l) j_l(kr)
Cr           = 1/2*(2l-1)!! (avw)^(-l) aj
Cr           = fi(OKA) * (r/w)^l
Cr     Define Neumann function as:
Cr       N_l = K_l - i kap e^l 2/((2l-1)!!)^2 J_l
Cr     avw is some arbitrary length scale, e.g. average WSR.
Cr     By setting loka=1 the fi and gi are rescaled as follows:
Cr       fi -> fi(OKA) = fi * (2l-1)!!/2
Cr       gi -> gi(OKA) = gi / (2l-1)!!
Cr     Limiting cases for energy->0
Cr       H_l -> (r/w)^(-l-1)         J_l -> (r/w)^l/(2(2l+1))
Cr     NB: in exact MTO paper, Andersen makes definitions j~ and n~:
Cr       n_l~(kappa r) =  (kappa r)^l+1/(2l-1)!! n_l(kappa r)
Cr                     -> 1 - (kappa r)^2/2(2l+3) + ... as kappa r -> 0
Cr       j_l~(kappa r) =  (2l+1)!!/(kappa r)^l j_l(kappa r)
Cr                     -> 1 - (kappa r)^2/2(2l+3) + ... as kappa r -> 0
Cr
Cr  *loka=2 conventions
Cr       Kbar_l = (2l-1)!! * h_l(akap r) * i*(akap)^(l+1) [Same as OKA]
Cr       Jbar_l = (2l+1)!! * j_l(akap r) / (akap)**l  [This is JOKA * 2 * (2l+1)]
Cr
Cr  *Generation of fi and gi:  fi (and gi for y>0) are calculated for
Cr   lmax and lmax-1 by an expansion in powers of x^2=y:
Cr
Cr                     (-x^2/2)^k
Cr       fi =  Sum_k  --------------  = dum(lmx+1-l)
Cr                    k! (2l+1+2k)!!
Cr
Cr       gi = dum(lmx+2+l), y>0
Cr
Cr     For y<0:
Cr       gi(l=0) = exp(-sqrt(-y))
Cr       gi(l=1) = (1d0+sqrt(-y))*exp(-sqrt(-y))
Cr       gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
Cr
Cr     The remaining functions are obtained by the recursion formula:
Cr
Cr       j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x)   l=lmx-2,lmx-3,..,-lmx-1
Cr
Cr     => dum(k) = (2*lmx+5-2*k)*dum(k-1)-y*dum(k-2)   k=3,4,...,2*lmx+2
Cr
Cr     and the Neumann function are given through:
Cr
Cr           n_l(x) = j_{-l-1}*(-1)^{l+1}
Cr
Cr  *Special cases:
Cr     fi(y,l) -> 1/(2l+1)!!     for y->0
Cr     gi(y,l) -> (2l-1)!!       for y->0
Cr     fi(y,l=0) = sinh(x)/x     for y<=0    x=sqrt(-e)*r, akap=sqrt(e) Im(akap)>=0
Cr     gi(y,l=0) = exp(i*akap*r) for y<=0    akap=sqrt(e) Im(akap)>=0
Cr               = exp(-x)       for y<=0    x=sqrt(-e)*r, akap=sqrt(e) Im(akap)>=0
Cr     gi(y,l=0) = cos(akap*r)   for y>0     akap=sqrt(e) Im(akap)=0
Cr   As mentioned above for OKA=1
Cr     fi(OKA) = fi * (2l-1)!!/2  and
Cr     gi(OKA) = gi / (2l-1)!!    so
Cr     J(E,r) = fi(OKA) * (r/w)^l     -> (r/w)^l/(2(2l+1)) for e->0
Cr     K(E,r) = gi(OKA) / (r/w)^(l+1) -> (r/w)^(-l-1) for e->0
Cr   And for OKA = 2
Cr     fi(OKA) = fi * (2l+1)!!    and
Cr     gi(OKA) = gi / (2l-1)!!    so
Cr     Jbar(E,r) -> (r/w)^l for e->0
Cr     Kbar(E,r) -> (r/w)^(l+1) -> for e->0
Cb Bugs
Cb   Program never checked for lmax < 0
Cb   For y > 100 the internal algorithm is numerically unstable !!!
Cu Updates
Cu   03 Jan 18 Added loka=2 convention
Cu   22 Aug 17 Optionally calls dbesnu for high accuracy.  Especially important for y>100.
Cu   23 Jul 08 bug fix: besslr doesn't make fi/gi(lmax+1) when lmax=0
Cu   19 May 04 Changed loka from logical to integer
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer loka
      integer lmin,lmax
      double precision y,fi(lmin:lmax),gi(lmin:lmax)
C Local variables:
      integer i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2,opt1,lmin1
      integer, parameter :: nlmax=20
      double precision dt,dt2,exppr,my,srmy,g1,t,
     .                 dum(nlmax*4+2),fac2l(-nlmax:nlmax*2+3)
C     double precision x,fack(0:2*nlmax+2),ak(0:2*nlmax+2)

      double precision fil(0:1),gil(0:1)
      real(8), parameter :: tol=1.d-15
      logical, parameter :: lhank = .true.
C Intrinsic functions:
      intrinsic dabs,dexp,dsqrt,max0

      if (lmin > 0) call rx('BESSLR : lmin gt 0')
      if (lmax < 0) call rx('BESSLR : lmax lt 0')

      opt1 = mod(loka/10,10)
      lmx = max0(lmax,2)
      if (lmx > nlmax+nlmax) call rxi(' BESSL : lmax gt nlmax*2, lmax=',lmx)

C --- A table of fac2l(l)=(2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      fac2l(0) = 1d0
      do  l = 1, lmx+1
        fac2l(l) = fac2l(l-1)*(l+l-1)
      enddo
      do  l = -1, lmin,-1
        fac2l(l) = fac2l(l+1)/(l+l+1)
      enddo

C --- Case call bessjs ---
      if (y >= 90 .and. opt1 == 1 .or. opt1 == 2) then
        call bessjs(y,lmax,13,fi(0),gi(0))
        if (lmin >= 0) goto 100
        lmin1 = -1
        if (lmax < 1) then
          call bessjs(y,1,13,fil,gil)
          fi(-1) =  fil(0) - y*fil(1)
          gi(-1) = (gil(0) - gil(1))/y
          lmin1 = -2
        endif
C       Warning ! if lmin is very negative and y large, recursion can be unstable
        do  l = lmin1, lmin, -1
          dt = 2*l+3
          fi(l) =  dt*fi(l+1) - y*fi(l+2)
          gi(l) = (dt*gi(l+1) - gi(l+2))/y
        enddo
        goto 100
      endif

C --- Case akap=0 ---
      if (y == 0) then
        do  l = lmin, lmax
          fi(l) = 1/fac2l(l+1)
          gi(l) = fac2l(l)
        enddo
        goto 100
      endif
      my = -y

C --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
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

C --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
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

C --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
      ll1 = lmx + lmx + 1
      ll2 = ll1 + 1
      nf = ll1
      do  k = 3, ll2
        nf = nf-2
        dum(k) = nf*dum(k-1) - y*dum(k-2)
      enddo

C --- Get fi and gi from dum ---
      lmxp1 = lmx+1
      lmxp2 = lmx+2
      isn = (-1)**lmin
      do  k = lmin, lmax
        j1 = lmxp1-k
        j2 = lmxp2+k
        fi(k) = dum(j1)
c   ... n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi(k) = dum(j2)*isn
        isn = -isn
      enddo

C --- Bessel function for large argument ... not finished.  See dlmf.nist.gov/10.49 ---
C      x = dsqrt(y)
C      fack(0) = 1
C      do  l = 1, 2*lmx+2
C        fack(l) = fack(l-1)*l
C      enddo
C      dt = 1
C      do  k = 0, lmx
C        ak(k) = fack(lmx+k)/fack(lmx-k)/fack(k)/dt
C        dt = dt*2
C      enddo
C      t = 0
C      srmy = 1
C      dt2 = 1
C      do  k = 0, (lmx+1)/2
C        dt = srmy*ak(2*k)/dt2
C        dt2 = dt2*y
C        t = t+dt
C        srmy = -srmy
C      enddo


C --- For E<0, use Hankel functions rather than Neumann functions ---
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
          if (lmin <= -2) then
            do  l = -2, lmin,-1
              gi(l) = ((l+l+3)*gi(l+1) - gi(l+2))/y
            enddo
          endif
        endif
        exppr = 1d0/dexp(srmy)
        forall (l=lmin:lmax) gi(l) = gi(l)*exppr
      endif

  100 continue
C --- Scaling to Andersen's 2nd generation LMTO conventions ---
      if (mod(loka,10) == 1) then
        do  l = lmin, lmax
          fi(l) = fi(l)*fac2l(l)*0.5d0
          gi(l) = gi(l)/fac2l(l)
        enddo
      endif

C --- Scaling to conventions mode 2
      if (mod(loka,10) == 2) then
        do  l = lmin, lmax
          fi(l) = fi(l)*fac2l(l)*(2*l+1)
          gi(l) = gi(l)/fac2l(l)
        enddo
      endif

      end
      subroutine bessl2(y,lmin,lmax,fi,gi)
C- Radial part of Bessel functions, Andersen's definitions
C  See besslr for definitions.
C  For y=0, fi(l) = 0.5/(2l+1), gi(l) = 1
      implicit none
C Passed variables:
      integer lmin,lmax
      double precision y,fi(lmin:lmax),gi(lmin:lmax)

      call besslr(y,1,lmin,lmax,fi,gi)
      end
      subroutine besslm(y,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C  See besslr for definitions.
C  For y=0, fi(l) = (2l-1)!!, gi(l) = 1/(2l+1)!!
      implicit none
C Passed variables:
      integer lmax
      double precision y,fi(0:lmax),gi(0:lmax)

      call besslr(y,0,0,lmax,fi,gi)
      end
      subroutine bessl(y,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C  See besslr for definitions.
C  For y=0, fi(l) = (2l-1)!!, gi(l) = 1/(2l+1)!!
      implicit none
C Passed variables:
      integer lmax
      double precision y,fi(0:lmax),gi(0:lmax)

      call besslr(y,0,0,lmax,fi,gi)
      end
      subroutine bessjs(y,lmax,opt,jr,hr)
C- Returns spherical Bessel and/or Hankel functions, calling dbesnu
C ----------------------------------------------------------------------
Ci Inputs
Ci   y     :y = x**2, where x is argument to Bessel
Ci         :y > 0 => return spherical Bessel/Neumann function
Ci         :y < 0 => return modified spherical Bessel/Hankel function
Ci   lmax  :maximum l for a given site
Ci   opt   :1s digit
Ci         : 0 return nothing in jr
Ci         : 1 return spherical Bessel in jr if y>0
Ci         :   return spherical Bessel of Im x in jr if y<0
Ci         : 2 return spherical Hankel in hr if y>0
Ci         :   return spherical Hankel of Im x in hr if y<0
Ci         : combination of 1+2 is allowed
Ci         :10s digit
Ci         : 0 return spherical Bessel function in jr and/or spherical Hankel in hr
Ci         :   depending on 1s digit opt
Ci         : 1 scale jr by 1/x**l and/or hr by x**(l+1)
Ci         :100s digit
Ci         : 0 always call dbesnu
Ci         : 1 call besslr for y<90
Ci         : 2 always call besslr
Co Outputs
Ci   jr    : spherical Bessel and/or functions, possibly scaled by r^l
Cr Remarks
Cr   This routine and besslr generate the same results to a relative
Cr   precision of better 10^-14 for y<90.  For y>90 besslr becomes unstable
Cu Updates
Cu   21 Aug 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,opt
      double precision y,jr(0:lmax),hr(0:lmax)
C ... Local parameters
      logical limg
      integer ierr,i,l,opt0,opt1,opt2
      integer, parameter :: ltop = 50
      real(8), parameter :: pi = 4*datan(1d0)
      double precision x,xx(1),dfac,fi(0:ltop),gi(0:ltop)
      integer k2
      double precision twolpk,tlp1,my
      real(8), parameter :: tol=1d-15

      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt2 = mod(opt/100,10)

      if (y == 0) then
        if (mod(opt0,2) == 1) then
          jr(0) = 1
          dfac = 1
          do  l  = 0, lmax
            dfac = dfac*(2*l+1)
            jr(l) = 1/dfac
          enddo
          if (opt1 == 0) jr(1:lmax) = 0
        endif
        if (opt0 >= 2) then
          hr(0) = 1
          dfac = 1
          do  l  = 0, lmax
            hr(l) = dfac
            dfac = dfac*(2*l+1)
          enddo
          if (opt1 == 0) call rx('bessls: dbesnu invalid argument for hankel, x=0')
        endif
        return

      elseif (y < 90 .and. opt2 == 1 .or. opt2 == 2) then
        if (opt0 > 2) then
          call besslr(y,0,0,lmax,jr,hr)
        else
          if (lmax > ltop) call rx('bessjs: increase ltop')
          call besslr(y,0,0,lmax,fi,gi)
          if (mod(opt0,2) == 1) jr(0:lmax) = fi(0:lmax)
          if (opt0 >= 2) hr(0:lmax) = gi(0:lmax)
        endif
        if (opt1 == 0) then
          call rx('bessls: add scaling to branch alling besslr')
        endif
        return
      endif

      limg = y<0
      x = sqrt(abs(y))
      l = 11 ; if (limg) l = 13
      ierr = 0
      if (mod(opt0,2) == 1) then
        call dbesnu(x, lmax+0.5d0, l, jr, xx, ierr)
C       dbesnu has trouble for small x and high l.  Use polynomial expansion
C       print *, '!!'; ierr = 1 ! For debugging
        if (ierr > 0 .and. x < 1) then
          my = -y
          tlp1 = 2*lmax+1
          ierr = 0
          jr(ierr:lmax) = 1
          fi(0:lmax+1-ierr) = 1
          do  k2 = 2, 1000, 2
            do  i = 0, lmax-ierr
              twolpk = k2+tlp1
              fi(i) = fi(i)*my/(k2*(twolpk-2*i))
              jr(lmax-i) = jr(lmax-i)+fi(i)
            enddo
            if (dabs(fi(lmax)) < tol) exit
          enddo
          dfac = 1
          do  i = 0, lmax
            dfac = dfac*(2*i+1)
C           if (i < ierr) cycle
C            if (i == 8) then
C              print *, 'hi'
C            endif
            jr(i) = jr(i)/dfac*x**i
          enddo
        else
          jr(0:lmax) = sqrt(pi/2/x) * jr(0:lmax)
        endif
        if (ierr /= 0) call rx('bessls: dbesnu returned with error')
      endif
      l = l+1
      if (opt0 >= 2) then
        call dbesnu(x, lmax+0.5d0, l, hr, xx, ierr)
        if (ierr /= 0) call rx('bessls: dbesnu returned with error')
        dfac = -sqrt(pi/2/x); if (limg) dfac = sqrt(2/pi/x)
        hr(0:lmax) = dfac * hr(0:lmax)
      endif

      if (opt1 /= 0) then
        if (mod(opt0,2) == 1) then
          forall (l=1:lmax) jr(l) = jr(l)/x**l
        endif
        if (opt0 >= 2) then
          forall (l=0:lmax) hr(l) = hr(l)*x**(l+1)
        endif
      endif

      end

C      subroutine fmain
C      implicit none
C      integer lmin,lmax,il
C      double precision fi(-3:10),gi(-3:10),y
C
C      fi = -99d0
C      gi = -99d0
C      lmin = -1
C      lmax =  2
C      y = 1.1d0
C
C      call besslr(y,0,lmin,lmax,fi(lmin),gi(lmin))
C      print *, sngl(fi(lmin:lmax+1))
C      print *, sngl(gi(lmin:lmax+1))
C
C      end
