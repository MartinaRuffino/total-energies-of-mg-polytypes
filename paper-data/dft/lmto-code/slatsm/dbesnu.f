      subroutine dbesnu(x, nu, opt, fr, fi, ierr)
C- Bessel functions of fractional order, real argument
C ----------------------------------------------------------------------
Ci Inputs
Ci   x     : evaluate Bessel at x
Ci   nu    : Bessel order
Ci   opt   :1s digit specifies which kind of Bessel function
Ci         : 0  do nothing
Ci         : 1  Bessel function J_nu
Ci         : 2  Neumann function Y_nu, aka N_nu
Ci         : 3  Modified Bessel of 1st kind I_nu(x) = i^nu J_nu(i x)
Ci         : 4  Modified Bessel of 2nd kind K_nu(x) (see Remarks)
Ci         :10s digit specifies which orders to return
Ci         : 0  returns Bessel_nu(x)
Ci         : 1  return Bessel_(i+alpha)(x), i=0...int(nu)
Ci         :100s digit specifies whether to make function or derivative
Ci         : 0  returns Bessel_nu
Ci         : 1  returns derivative of Bessel_nu wrt x
Ci         : 2  returns 2nd derivative of Bessel_nu wrt x
Co Outputs
Co   fr    : real part of Bessel function at x
Co   fi    : imaginary part of Bessel function at x
Co         : (not used now)
Co   ierr  : 0 function call was successful
Co         :-1 function call was not successful
Co         :>0 only bessels of order
Co              0+(fractional part nu) ... ierr+(fractional part nu)
Co             could be calculated.  Bessel(nu) is not returned
Cl Local variables
Cl         :
Cr Remarks
Cr   Hankel functions of 1st and 2nd kind are related to J and Y as:
Cr     H_nu^(1) (x) =  J_nu(x) + i Y_nu(x)
Cr                  = [J_(-nu)(x) - e^(i pi nu) J_nu(x)]/[i sin(pi nu)]
Cr     H_nu^(2) (x) =  J_nu(x) - i Y_nu(x)
Cr   K_nu is related to H (and also J) as
Cr     K_nu (x) = pi/2 [I_(-nu)(x) - I_nu(x)] / sin(pi nu)
Cr              = pi/2 i^(nu+1) H_nu^(1) (i x)
Cr   Spherical Bessel, Neumann, and Hankel functions:
Cr     j_n(x) = sqrt(pi/2/x) J_(n+1/2) (x)
Cr     n_n(x) = sqrt(pi/2/x) Y_(n+1/2) (x)
Cr     h_n^(1) (x) =  j_n(x) + i y_n(x)
Cr     h_n^(2) (x) =  j_n(x) - i y_n(x)
Cu Updates
Cu   19 May 11  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer  opt, ierr
      double precision nu, x, fr(*), fi(*)
C ... Local parameters
      integer nb,opt0,opt1,opt2,nextra,i
      double precision alpha,b(1000),fac,fac2

C     Number of bessels to calculate
      nb = int(nu)+1
C     Fractional order
      alpha = nu+1 - nb
      if (alpha >= 1) then
        alpha = alpha - 1
        nb = nb+1
      endif
      if (nb > 1000) call rx('dbesnu: increase size of b')
      ierr = 0
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt2 = mod(opt/100,10)

      nextra = 0
      if (opt2 /= 0) nextra=1
C ... Bessel J
      if (opt0 == 1) then
        call djbesl(x, alpha, nb+nextra, b, ierr)
        fac = -1
        fac2 = 1
C ... Neumann N or Y
      elseif (opt0 == 2) then
        call dybesl(x, alpha, nb+nextra, b, ierr)
        fac = -1
        fac2 = 1
C ... Modified Bessel I, or Bessel of imaginary argument
      elseif (opt0 == 3) then
        call dibesl(x, alpha, nb+nextra, 1, b, ierr)
        fac = 1
        fac2 = -1
C ... Hankel H, or Neumann of imaginary argument
      elseif (opt0 == 4) then
        call dkbesl(x, alpha, nb+nextra, 1, b, ierr)
        fac = -1
        fac2 = -1
      else
        call rx('dbesnu: bad opt')
      endif

      if (ierr < 0) return
      if (ierr >= nb+nextra) ierr = 0

C     Derivative from recurrence relation
      if (opt1 > 0) then
        if (opt2 == 0) then
          call dcopy(nb,b,1,fr,1)
        else
          do  i = 1, nb
            fr(i) = b(i)*(i-1+alpha)/x + fac*b(i+1)
          enddo
C         Warning ... never checked
          if (opt2 == 2) then
          do  i = 1, nb
            fr(i) = b(i)*((nu/x)**2-fac2) - fr(i)/x
          enddo
          endif
        endif
      else
        if (opt2 == 0) then
          fr(1) = b(nb)
        else
C         1st derivataive
          fr(1) = b(nb)*nu/x + fac*b(nb+1)
          if (opt2 == 2) then
            fr(1) = b(nb)*((nu/x)**2-fac2) - fr(1)/x
          endif
        endif
      endif

      end
C     Tests or illustrates various branches of dbesnu.
C     link with: lk dbesnu.o ~/lm/subs/subs.a
C     To compare Bessels of 1/2 order to besslr:
C       a.out --sph
C     ... the following check functions for wide range of x:
C       a.out --sph~chkbes
C       a.out --sph~chkhan
C       a.out --sph~chkbes~negx
C       a.out --sph~chkhan~negx
C       a.out --sph
C     To tabulate Bessels for fixed x vs nu
C       a.out --jnu
C     To tabulate radial derivative of Bessels for fixed x vs nu
C       a.out --djnu
C     To tabulate Bessels for fixed nu vs x
C       a.out --jx
C     To tabulate radial derivative of Bessels for fixed nu vs x
C       a.out --djx
C     To tabulate 2nd radial derivative of Bessels for fixed nu vs x
C       a.out --ddjx
C     To evaluate kinetic energy in MT potential at a point
C       a.out --ke
C     To check Wronskian at x
C       a.out --wronsk
C     Solutions to SE in an exponential-MT potential
C       a.out --psi1mt

C      subroutine fmain
C      implicit none
C      logical lneg
C      integer ierr,i,j,ir,nx,mode,fast
C      double precision nu,x,xx,pi,psip,psipp
C      real(8), target :: bess(10),hank(10),bes3(10),han3(10)
C      double precision a,b,s,potmt,v,e,z,coff(2,4)
CC     double precision psi(1000,2)
C      character outs*120, dc*1
C      real(8), pointer :: f1(:),f3(:)
C      double precision dnu,xn,fn,xtol,ftol,dxmn,wk(12),oldnu,kap,y
C      procedure(logical) cmdopt
C      procedure(integer) wordsw,isw
C      data xtol/1D-12/,ftol/1D-12/,dxmn/1D-13/
C
C      pi = 4*datan(1d0)
C
C      if (cmdopt('--sph',5,0,outs)) then
C      nu = 0.5d0
C      dc = outs(6:6)
CC ... Compare dbesnu and besslr tabulated on a uniform mesh
C      mode = 0
C      if (wordsw(outs,dc,'chkbes',dc//' ',j) > 0) mode = 1
C      if (wordsw(outs,dc,'chkhan',dc//' ',j) > 0) mode = 2
C      if (mode > 0) then
C
C      lneg = wordsw(outs,dc,'negx',dc//' ',j) > 0
C      fast = wordsw(outs,dc,'fast',dc//' ',j)
C      if (fast > 0) fast = 1
C
C      call info2(1,0,0,
C     .    '# ... testing %?#(n==1)#Bessel#Hankel# functions of half integer order, '//
C     .    '%?#(n==1)#negative#positive# argument x^2',
C     .    mode,isw(lneg))
C
C      do  j = -10, -1, 2
C        x = 10d0**j; if (j == -10) x = 0
C
C
C        y = x*x; if (lneg) y = -y
C        call bessjs(y,4-1,100*fast+10+mode,bess,hank)
C
CC       Substitute a direct call to dbesnu for debugging purposes
CC        if (x > 0) then
CC          if (lneg) then
CC            call dbesnu(x, 4+nu, 13, bess, xx, ierr)
CC            call dbesnu(x, 4+nu, 14, hank, xx, ierr)
CC            bess(1:4) = sqrt(pi/2/x) * bess(1:4)
CC            hank(1:4) = sqrt(2/pi/x) * hank(1:4)
CC          else
CC            call dbesnu(x, 4+nu, 11, bess, xx, ierr)
CC            call dbesnu(x, 4+nu, 12, hank, xx, ierr)
CC            bess(1:4) =  sqrt(pi/2/x) * bess(1:4)
CC            hank(1:4) = -sqrt(pi/2/x) * hank(1:4)
CC          endif
CC          if (ierr /= 0) stop 'oops'
CC          do  i = 2, 4
CC            bess(i) = bess(i)/x**(i-1)
CC          enddo
CC          do  i = 1, 4
CC            hank(i) = hank(i)*x**i
CC          enddo
CC        endif
C
C        call besslr(y,0,0,4-1,bes3,han3)
CC       call ropbes(x,1d0,4-1,xn,wk,bes3,1,1)
C
C        f1 => bess; f3 => bes3
C        if (mode == 2) then
C          f1 => hank; f3 => han3
C        endif
C        print 335, x, f1(1:4), (f1(1:4)-f3(1:4))/f1(1:4)
C  335   format(5f16.10,2x,1p,5e12.2)
C      enddo
C
C      x = 1
C      do while (x < 10d0)
C        y = x*x; if (lneg) y = -y
C
C        call bessjs(y,4-1,100*fast+10+mode,bess,hank)
C        call besslr(y,0,0,4-1,bes3,han3)
C        f1 => bess; f3 => bes3
C        if (mode == 2) then
C          f1 => hank; f3 => han3
C        endif
C        print 335, x, f1(1:4), (f1(1:4)-f3(1:4))/f1(1:4)
C        x = x + 0.5d0
C      enddo
C
C      x = 10
C      do while (x < 100d0)
C        y = x*x; if (lneg) y = -y
C        call bessjs(y,4-1,100*fast+10+mode,bess,hank)
C        call besslr(y,0,0,4-1,bes3,han3)
C        f1 => bess; f3 => bes3
C        if (mode == 2) then
C          f1 => hank; f3 => han3
C        endif
C        print 335, x, f1(1:4), (f1(1:4)-f3(1:4))/f1(1:4)
C        x = x + 5d0
C      enddo
C
C      f1(1:4) = [1d0,2d0,3d0,4d0]
C      f3(1:4) = [2d0,4d0,6d0,8d0]
C
C
C      call cexit(0,1)
C      endif
C
C      print *, '... testing bessel functions of half integer order'
C
CC     Show connection to spherical bessel functions
C      x = 2.7d0
CC     print *, '!!' ; x = 1
C      call dbesnu(x, 4+nu, 11, bess, xx, ierr)
C      if (ierr /= 0) stop 'oops'
C      bess(1:4) = sqrt(pi/2/x) * bess(1:4)
C      print
C     .  "(/' 1st 4 spherical Bessel functions from dbesnu, x=',f10.6)",x
C      print 333, ' ', bess(1:4)
C  333 format(a,10f15.10)
C      call besslr(x*x,0,0,4,bes3,han3)
C      print *, '1st 4 Bessel functions from besslr'
C      do  i = 2, 4
C        bes3(i) = bes3(i)*x**(i-1)
C      enddo
C      do  i = 1, 4
C        han3(i) = han3(i)/x**(i)
C      enddo
C      print 333, ' ', bes3(1:4)
CC     print *, 'difference'
C      print 334, 'd', bess(1:4)-bes3(1:4)
C  334 format(a,1p,10e15.5)
C
C      print 301, 'spherical Neumann', x
C      call dbesnu(x, 4+nu, 12, bess, xx, ierr)
C      bess(1:4) = -sqrt(pi/2/x) * bess(1:4)
C      print 333, ' ', bess(1:4)
C      print *, '1st 4 Neumann functions from besslr'
C      print 333, ' ', han3(1:4)
CC     print *, 'difference'
C      print 334, 'd', bess(1:4)-han3(1:4)
C
C      call dbesnu(x, 4+nu, 13, bess, xx, ierr)
C      if (ierr /= 0) stop 'oops'
C      bess(1:4) = sqrt(pi/2/x) * bess(1:4)
C      print 301, 'modified spherical Bessel', x
C  301 format (/' 1st 4 ',a,' functions',
C     .  ' from dbesnu, x=',f10.6)
C      print 333, ' ', bess(1:4)
C      call besslr(-x*x,0,0,4,bes3,han3)
C      print *, '1st 4 Bessel functions from besslr, Im x'
C      do  i = 2, 4
C        bes3(i) = bes3(i)*x**(i-1)
C      enddo
C      do  i = 1, 4
C        han3(i) = han3(i)/x**(i)
C      enddo
C      print 333, ' ', bes3(1:4)
CC     print *, 'difference'
C      print 334, 'd', bess(1:4)-bes3(1:4)
C
C      print 301, 'spherical Hankel', x
C      call dbesnu(x, 4+nu, 14, bess, xx, ierr)
C      if (ierr /= 0) stop 'oops'
C      bess(1:4) = sqrt(2/pi/x) * bess(1:4)
C
C      call besslr(-x*x,0,0,4,bes3,han3)
C      print 333, ' ', bess(1:4)
C      print *, '1st 4 Hankel functions from besslr, Im x'
C      do  i = 1, 4
C        han3(i) = han3(i)/x**(i)
C      enddo
C      print 333, ' ', han3(1:4)
CC     print *, 'difference'
C      print 334, 'd', bess(1:4)-han3(1:4)
C
C      else if (cmdopt('--jnu',5,0,outs)) then
C      x = 3d0
C      x = 10d0
C      print 302, x
C  302 format(' Plot bessel functions at x = ',f8.4,
C     .  ' as function of nu (file 80)')
C      print *, 'The following',
C     .  ' shows transition from oscillatory'
C      print *, 'to exponential behavior near x=nu',
C     .  ' (compare x=3 and x=20)'
C      print *,
C     .  "  fplot -frme:ly 0,1,0,1 -ord 'abs(y)' -colsy 2:5 fort.80"
C
C      do  nu = 0.0d0, 2*x, .1d0
CC      do  nu = 0.0d0, 20, .1d0
C        call dbesnu(x, nu, 1, bess(1), xx, ierr)
CC        if (ierr /= 0) then
CC          print *, int(nu+1),int(nu)+1
CC        call dbesnu(x, nu, 1, bess(1), xx, ierr)
CC        endif
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 2, bess(2), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 3, bess(3), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 4, bess(4), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C
CC       For plotting cut at maximum 1000
C        do  i = 1, 4
C          bess(i) = dmin1(dmax1(bess(i),-1d3),1d3)
C        enddo
C
C        write(80,401) nu, bess(1:4)
C  401   format(f8.4, 1p, 4e18.10)
C      enddo
C
CC fplg -frme 0,1,0,.7 -x 0,14 -y -.4,.4 -colsy 2 -lt 1,col=1,0,0 djnu10 -lt 1,col=0,1,0 djnu5 -lt 1,col=0,0,1 djnu3 -lbl 0,.4 rc "&{dJ}/&{dz}" -xl '~{n}'
CC fplg -frme 0,1,0,.7 -x 0,14 -y -1,1   -colsy 3 -lt 1,col=1,0,0 djnu10 -lt 1,col=0,1,0 -colsy 3 djnu5 -lt 1,col=0,0,1  -colsy 3 djnu3 -lbl 0,.8 rc "&{dN}/&{dz}" -xl '~{n}'
C      else if (cmdopt('--djnu',6,0,outs)) then
C      x = 10d0
CC      x = 20d0
C      x = 10d0
C
C      print 303, x
C  303 format(' Plot derivative of bessel functions at x = ',f8.4,
C     .  ' as function of nu (file 80)')
C      print *, 'The following',
C     .  ' shows transition from oscillatory'
C      print *, 'to exponential behavior near x=nu',
C     .  ' (compare x=3 and x=20)'
C      print *,
C     .  "  fplot -frme:ly 0,1,0,1 -ord 'abs(y)' -colsy 2:5 fort.80"
C
CC      do  nu = 0.0d0, 2*x, .1d0
C       do  nu = 0.0d0, 20, .1d0
C        call dbesnu(x, nu, 101, bess(1), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 102, bess(2), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 103, bess(3), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 104, bess(4), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C
CC       For plotting cut at maximum 1000
C        do  i = 1, 4
C          bess(i) = dmin1(dmax1(bess(i),-1d3),1d3)
C        enddo
C
C        write(80,401) nu, bess(1:4)
C      enddo
C
C      elseif (cmdopt('--jx',4,0,outs)) then
C      nu = 2.7d0
C      print 304, nu
C  304 format(' Plot bessel functions at nu = ',f8.4,
C     .  ' as function of x (file 80)')
C      print *,
C     .  "  fplot -frme:ly 0,1,0,1 -ord 'abs(y)' -colsy 2:5 fort.80"
C
CC      do  nu = 0.0d0, 2*x, .1d0
C      do  x = 0.1d0, 20, .1d0
C        call dbesnu(x, nu, 1, bess(1), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 2, bess(2), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 3, bess(3), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 4, bess(4), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C
CC       For plotting cut at maximum 1000
C        do  i = 1, 4
C          bess(i) = dmin1(dmax1(bess(i),-1d3),1d3)
C        enddo
C
C        write(80,401) x, bess(1:4)
C      enddo
C
C      else if (cmdopt('--djx',5,0,outs)) then
C      nu = 2.7d0
C      print 305, nu
C  305 format(' Plot bessel function derivative for nu = ',f8.4,
C     .  ' as function of x (file 80)')
C      print *,
C     .  "  fplot -frme:ly 0,1,0,1 -ord 'abs(y)' -colsy 2:5 fort.80"
C
CC      do  nu = 0.0d0, 2*x, .1d0
C      do  x = 0.1d0, 20, .1d0
C        call dbesnu(x, nu, 101, bess(1), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 102, bess(2), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 103, bess(3), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 104, bess(4), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C
CC       For plotting cut at maximum 1000
C        do  i = 1, 4
C          bess(i) = dmin1(dmax1(bess(i),-1d3),1d3)
C        enddo
C
C        write(80,401) x, bess(1:4)
C      enddo
C
C
C      else if (cmdopt('--ddjx',6,0,outs)) then
C      nu = 2.7d0
C      print 205, nu
C  205 format(' Plot bessel function 2nd derivative for nu = ',f8.4,
C     .  ' as function of x (file 80)')
C      print *,
C     .  "  fplot -frme:ly 0,1,0,1 -ord 'abs(y)' -colsy 2:5 fort.80"
C
CC      do  x = 1.2d0, 3d0, 3d0
C       do  x = 0.1d0, 20, .1d0
C        call dbesnu(x, nu, 201, bess(1), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 202, bess(2), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 203, bess(3), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C        call dbesnu(x, nu, 204, bess(4), xx, ierr)
C        if (ierr /= 0) call rx('oops')
C
CC       For plotting cut at maximum 1000
C        do  i = 1, 4
C          bess(i) = dmin1(dmax1(bess(i),-1d3),1d3)
C        enddo
C
C        write(80,401) x, bess(1:4)
C      enddo
C
C      else if (cmdopt('--wronsk',8,0,outs)) then
C        nu = 2.7d0
C        x = 1.9d0
CC       Wronskian J Y' - J' Y
C        call dbesnu(x, nu, 001, bess(1), xx, ierr)
C        call dbesnu(x, nu, 101, bess(2), xx, ierr)
C        call dbesnu(x, nu, 002, bess(3), xx, ierr)
C        call dbesnu(x, nu, 102, bess(4), xx, ierr)
C        print 306, x,nu,
C     .    " J Y' - J' Y =",
C     .    (bess(1)*bess(4) - bess(2)*bess(3)),
C     .    "      2/pi/x =",
C     .    2/pi/x,
C     .    (bess(1)*bess(4) - bess(2)*bess(3)) - 2/pi/x
C  306   format(/' Wronskian at x =',f8.4,' nu = ',f8.4/
C     .    a,f15.10/
C     .    a,f15.10/
C     .    "  difference =",1p,e18.10)
C
CC       Wronskian I K' - K' I
C        call dbesnu(x, nu, 003, bess(1), xx, ierr)
C        call dbesnu(x, nu, 103, bess(2), xx, ierr)
C        call dbesnu(x, nu, 004, bess(3), xx, ierr)
C        call dbesnu(x, nu, 104, bess(4), xx, ierr)
C        print 306, x,nu,
C     .    " I K' - K' I =",
C     .    (bess(1)*bess(4) - bess(2)*bess(3)),
C     .    "        -1/x =",
C     .    -1/x,
C     .    (bess(1)*bess(4) - bess(2)*bess(3)) - (-1/x)
C
C
C      else if (cmdopt('--ke',4,0,outs)) then
C        x = 2.7d0
C        s = 2.7d0
C        b = 1.5d0
C        a = 2.5d0
C        e = -1.09565656838237
C
C        kap = sqrt(-e)
C        nu = kap*a
C        z = a*b*exp(-x/a)
C
C        xx = potmt(x,a,b,s)
C        print *,'v',xx
C        print *, 'e-v',(e-xx)
C
C        z = a*b*exp(-(x+1d-4)/a)
C        call dbesnu(z, nu, 001, bess(2), xx, ierr)
C        psip = bess(2)
C        z = a*b*exp(-(x-1d-4)/a)
C        call dbesnu(z, nu, 001, bess(3), xx, ierr)
C        psip = (bess(2)-bess(3))/2d-4
C        z = a*b*exp(-x/a)
C        call dbesnu(z, nu, 001, bess(1), xx, ierr)
C        psipp = (bess(2)+bess(3)-2*bess(1))/1d-4**2
C        call dbesnu(z, nu, 101, bess(2), xx, ierr)
C        print *, psip,-z/a*bess(2),psip-(-z/a*bess(2))
C        call dbesnu(z, nu, 201, bess(3), xx, ierr)
C        print *, psipp,  (z/a)**2*bess(3) - psip/a,
C     .    psipp - ((z/a)**2*bess(3) - psip/a)
C        psipp = (z/a)**2*bess(3) - psip/a
C
C        xx = psipp/bess(1)
C        print *, 'psipp/psi',xx
C
CC     Eigenvalues in an exponential potential
C      else if (cmdopt('--evexp',7,0,outs)) then
C        b = 1.5d0
C        a = 2.5d0
C        dnu = .1d0
C
C        print 299, b, a
C  299   format('  Eigenvalues of exponential potential.',
C     .    '  b = ',f8.4,'  a = ',f8.4)
C
C        do  nu = 0, 2*b*b, dnu
C          z = a*b
CC         zp = -z/a
C          call dbesnu(z, nu, 1, bess(1), xx, ierr)
C          if (ierr /= 0) call rx('oops')
C          call dbesnu(z, nu, 101, bess(2), xx, ierr)
C          if (ierr /= 0) call rx('oops')
C
CC        print *, nu, bess(1), bess(2)
C
C          if (nu > 0) then
C
CC         Look for a root in J
C          if (bess(1)*bess(3) < 0) then
C            ir = 0
C            xn = oldnu
C            fn = bess(3)
C            do  j = 1, 10
C              call pshpr(1)
C              call rfalsi(xn,fn,xtol,ftol,dxmn,dnu,0,wk,ir)
C              call poppr
C              if (ir == 0) exit
C              call dbesnu(z, xn, 1, fn, xx, ierr)
C              if (ierr /= 0) call rx('oops')
C            enddo
C            if (ir /= 0) call rx('failed to locate root')
C            kap = xn/a
C            e = -kap**2
C            print 307, 'jnu',xn,kap,e
C  307       format('#  found root in ',a,' nu = ',f20.15,
C     .        ' kappa =', f20.15,' e =', f20.15)
C          endif
C
CC         Look for a root in J'
C          if (bess(2)*bess(4) < 0) then
C            ir = 0
C            xn = oldnu
C            fn = bess(4)
C            do  j = 1, 10
C              call pshpr(1)
C              call rfalsi(xn,fn,xtol,ftol,dxmn,dnu,0,wk,ir)
C              call poppr
C              if (ir == 0) exit
C              call dbesnu(z, xn, 101, fn, xx, ierr)
C              if (ierr /= 0) call rx('oops')
C            enddo
C            if (ir /= 0) call rx('failed to locate root')
C            kap = xn/a
C            e = -kap**2
C            print 307, 'jnup',xn,kap,e
C
C          endif
C
C          endif
C
C          bess(3) = bess(1)
C          bess(4) = bess(2)
C          oldnu = nu
C
C        enddo
C
C
CC     SE in an exponential-MT potential (see notes)
C      else if (cmdopt('--psi1mt',7,0,outs)) then
C        b = 1.5d0
C        a = 2.5d0
C        s = 2.7
C
C        print 298, b, a, s, potmt(s,a,b,s)
C  298   format('  w.f. in exponential MT potential.',
C     .    '  b =',f7.4,'  a =',f7.4,'  s =',f7.4,'  V(s) =',f8.4)
C
C        print *, ' Write potential to file 80'
C        do  x = -3*s, 3*s, 0.01d0
C          v = potmt(x,a,b,s)
C          write(80,"(f8.4,f15.10)") x, v
C        enddo
C
CC       Find exact eigenvalue
C        e = -1d0
C        call pshpr(1)
C        call psi1mt(0,coff,x,a,b,s,e,xx)
C        call poppr
C        fn = coff(2,1)
C        ir = 0
C        do  j = 1, 10
C          xx = e
C          call pshpr(1)
C          call rfalsi(e,fn,xtol,ftol,dxmn,.1d0,0,wk,ir)
C          call poppr
C          if (ir == -4) then
C            e = 1.1d0*e-.1*xx
C          endif
C          if (ir == 0) exit
C          call pshpr(1)
C          call psi1mt(0,coff,x,a,b,s,e,xx)
C          call poppr
C          fn = coff(2,1)
C        enddo
C        if (ir /= 0) call rx('failed to locate root')
C        kap = sqrt(-e)
C        xx = sqrt(-e + potmt(s,a,b,s))
C        print 309, e,kap,xx
C  309   format('  Exact eigenvalue found at e=',f15.10,
C     .         ' kap =',f15.10,'  kapp =',f15.10)
C
CC  gen k.e. with:
CC  mc -f9f18.10 fort.81 -p -p -diff -diff -tog -de -e1 x2 -ccat >dat
CC  fplg -colsy 3 dat fort.80
CC  Confirm SE satisfied k.e. with:
CC  mc dat -p -coll 3 fort.80 -coll 2  -- -ccat -e2 x1 x4 >dat2
C        nx = 0
CC       e = -0.141236055252858d0
C        e = -1.083656568382370d0     ! eigenvalue of exponential V
CC       e = -1.0846999443d0          ! eigenvalue of MT V
C        print 310, e
C  310   format('  Write psi to file 81 for energy',f15.10)
C        call psi1mt(0,coff,x,a,b,s,e,xx)
C        do  x = -3*s, 3*s, 0.01d0
C          nx = nx+1
CC         psi(nx,1) = x
C          call psi1mt(1,coff,x,a,b,s,e,xx)
CC         psi(nx,2) = xx
C          if (abs(xx) < 10) then
C            write(81,"(f8.4,f15.10)") x, xx
C          endif
C        enddo
C
C      else
C        print *, 'usage: a.out [option]'
C        print *, 'Options:'
C        print *, '  --sph      Relation between MSM sph bessels and std'
C        print *, '  --jnu      Bessel function vs nu, fixed x'
C        print *, '  --djnu     d(Bessel)/dx vs nu, fixed x'
C        print *, '  --jx       Bessel function vs x, fixed nu'
C        print *, '  --djx      1st derivative of Bessel'
C        print *, '  --ddjx     2nd derivative of Bessel'
C        print *, '  --wronsk   Wronskian betw/ two types of Bessels'
C        print *, '  --ke       check k.e. of Bessel(z(x))'
C        print *, '  --evexp    eigenvalues of an exponential potential'
C        print *, '  --psi1mt   w.f. in 1 MT potential'
C      endif
C
C      end
C      double precision function potmt(x,a,b,s)
CC- MT Potential
C      double precision x,a,b,s
C
C      if (abs(x) > s) then
C        potmt = -b**2*exp(-2*s/a)
C      else
C        potmt = -b**2*exp(-2*abs(x)/a)
C      endif
C      end
C
C      subroutine psi1mt(mode,coff,x,a,b,s,e,psi)
CC- Wave function in single MT Potential
CC     mode=0 set up coff
C      implicit none
C      integer mode
C      double precision a,b,s,e,x,psi,coff(2,4)
C      double precision z,f1l,df1l,f2l,df2l
C      integer ierr,iprint
C      double precision kap,kapp,nu,f1,df1,f2,df2,wronsk,xx,zp
C      double precision env,envp
C      double precision pi,potmt
C      wronsk(f1,df1,f2,df2) = f1*df2 - df1*f2
C
C      kap = sqrt(-e)
CC     kapp = sqrt(-e - b**2*exp(-2*s/a))
C      kapp = sqrt(-e + potmt(s,a,b,s))
C      nu = kap*a
C
CC     Get coefficients to match psi at each junction
C      if (mode == 0) then
C
C        pi = 4*datan(1d0)
C
C        coff(1,4) = 0
C        coff(2,4) = 1
C        env = exp(-kapp*s)
C        envp = -kapp*env
C        z = a*b*exp(-s/a)
C        zp = -z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C
CC        env = f1l
CC        envp = df1l
C
C        coff(1,3) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,3) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
CC       Match at MT center
C        z = a*b
C        zp = -z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C        env = coff(1,3)*f1l + coff(2,3)*f2l
C        envp = coff(1,3)*df1l + coff(2,3)*df2l
C        df1l = -df1l
C        df2l = -df2l
C
C        coff(1,2) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,2) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
C        env = coff(1,2)*f1l + coff(2,2)*f2l
C        envp = coff(1,2)*df1l + coff(2,2)*df2l
C
CC       Match at left boundary
C        z = a*b*exp(-s/a)
C        zp = z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C        env = coff(1,2)*f1l + coff(2,2)*f2l
C        envp = coff(1,2)*df1l + coff(2,2)*df2l
C        f1l = exp(kapp*(-s))
C        df1l = kapp * f1l
C        f2l = exp(-kapp*(-s))
C        df2l = -kapp * f2l
C
C        coff(1,1) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,1) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
C
C        if (iprint() > 0) then
C          print 101, kap, kapp,
C     .      coff(:,1),coff(:,2),coff(:,3),coff(:,4)
C        endif
C  101   format('  kap = ',f15.10,'   kapp = ',f15.10/
C     .    '  c(I)  ', 2f15.10/
C     .    '  c(II) ', 2f15.10/
C     .    '  c(III)', 2f15.10/
C     .    '  c(IV) ',  2f15.10)
C        return
C
C      endif
C
C      if (.false.) then
C      elseif (x <= -s) then
C        f1l = exp(kapp*x)
C        f2l = exp(-kapp*x)
C        psi = coff(1,1)*f1l + coff(2,1)*f2l
C      elseif (x < 0) then
C        z = a*b*exp(-abs(x)/a)
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        psi = coff(1,2)*f1l + coff(2,2)*f2l
C      elseif (x <= s) then
C        z = a*b*exp(-x/a)
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        psi = coff(1,3)*f1l + coff(2,3)*f2l
C
CC        xx = potmt(x,a,b,s)
CC        print *, xx,(e-xx),(e-xx)*psi
C
CCC       First derivatives
CC        env = psi
CCC       First  derivative:
CC        call dbesnu(z, nu, 101, df1l, xx, ierr)
CC        call dbesnu(z, nu, 102, df2l, xx, ierr)
CC        psi = -z/a*(coff(1,3)*df1l + coff(2,3)*df2l)
CCC       Log derivative
CC        psi = psi/env
CCC       Kinetic energy:
CC        call dbesnu(z, nu, 101, df1l, xx, ierr)
CC        call dbesnu(z, nu, 102, df2l, xx, ierr)
CC        psi = z/a/a*(coff(1,3)*df1l + coff(2,3)*df2l)
CC        call dbesnu(z, nu, 201, df1l, xx, ierr)
CC        call dbesnu(z, nu, 202, df2l, xx, ierr)
CC        psi = (z/a)**2*(coff(1,3)*df1l + coff(2,3)*df2l) + psi
CC        print *, psi,psi/env
CC        stop
C
C      else
C        f1l = exp(kapp*x)
C        f2l = exp(-kapp*x)
C        psi = coff(1,4)*f1l + coff(2,4)*f2l
CC       Log derivative
CC       psi = (kapp*coff(1,4)*f1l - kapp*coff(2,4)*f2l)/psi
CC       First derivative:
CC       psi = kapp*coff(1,4)*f1l - kapp*coff(2,4)*f2l
CC       Kinetic energy:
CC       psi = kapp*kapp*(coff(1,4)*f1l + coff(2,4)*f2l)
C
C      endif
C
C      end
C
C      subroutine psi2mt(mode,coff,x,a,b,s,e,psi)
CC- Wave function in a double MT Potential
CC     mode=0 set up coff
C      implicit none
C      integer mode
C      double precision a,b,s,e,x,psi,coff(2,4)
C      double precision z,f1l,df1l,f2l,df2l
C      integer ierr,iprint
C      double precision kap,kapp,nu,f1,df1,f2,df2,wronsk,xx,zp
C      double precision env,envp
C      double precision pi,potmt
C      wronsk(f1,df1,f2,df2) = f1*df2 - df1*f2
C
C      kap = sqrt(-e)
CC     kapp = sqrt(-e - b**2*exp(-2*s/a))
C      kapp = sqrt(-e + potmt(s,a,b,s))
C      nu = kap*a
C
CC     Get coefficients to match psi at each junction
C      if (mode == 0) then
C
C        pi = 4*datan(1d0)
C
C        coff(1,4) = 0
C        coff(2,4) = 1
C        env = exp(-kapp*s)
C        envp = -kapp*env
C        z = a*b*exp(-s/a)
C        zp = -z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C
CC        env = f1l
CC        envp = df1l
C
C        coff(1,3) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,3) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
CC       Match at MT center
C        z = a*b
C        zp = -z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C        env = coff(1,3)*f1l + coff(2,3)*f2l
C        envp = coff(1,3)*df1l + coff(2,3)*df2l
C        df1l = -df1l
C        df2l = -df2l
C
C        coff(1,2) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,2) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
C        env = coff(1,2)*f1l + coff(2,2)*f2l
C        envp = coff(1,2)*df1l + coff(2,2)*df2l
C
CC       Match at left boundary
C        z = a*b*exp(-s/a)
C        zp = z/a
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 101, df1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        call dbesnu(z, nu, 102, df2l, xx, ierr)
C        df1l = zp*df1l
C        df2l = zp*df2l
C        env = coff(1,2)*f1l + coff(2,2)*f2l
C        envp = coff(1,2)*df1l + coff(2,2)*df2l
C        f1l = exp(kapp*(-s))
C        df1l = kapp * f1l
C        f2l = exp(-kapp*(-s))
C        df2l = -kapp * f2l
C
C        coff(1,1) = wronsk(env,envp,f2l,df2l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C        coff(2,1) =-wronsk(env,envp,f1l,df1l)/
C     .              wronsk(f1l,df1l,f2l,df2l)
C
C
C        if (iprint() > 0) then
C          print 101, kap, kapp,
C     .      coff(:,1),coff(:,2),coff(:,3),coff(:,4)
C        endif
C  101   format('  kap = ',f15.10,'   kapp = ',f15.10/
C     .    '  c(I)  ', 2f15.10/
C     .    '  c(II) ', 2f15.10/
C     .    '  c(III)', 2f15.10/
C     .    '  c(IV) ',  2f15.10)
C        return
C
C      endif
C
C      if (.false.) then
C      elseif (x <= -s) then
C        f1l = exp(kapp*x)
C        f2l = exp(-kapp*x)
C        psi = coff(1,1)*f1l + coff(2,1)*f2l
C      elseif (x < 0) then
C        z = a*b*exp(-abs(x)/a)
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        psi = coff(1,2)*f1l + coff(2,2)*f2l
C      elseif (x <= s) then
C        z = a*b*exp(-x/a)
C        call dbesnu(z, nu, 001, f1l, xx, ierr)
C        call dbesnu(z, nu, 002, f2l, xx, ierr)
C        psi = coff(1,3)*f1l + coff(2,3)*f2l
C
CC        xx = potmt(x,a,b,s)
CC        print *, xx,(e-xx),(e-xx)*psi
C
CCC       First derivatives
CC        env = psi
CCC       First  derivative:
CC        call dbesnu(z, nu, 101, df1l, xx, ierr)
CC        call dbesnu(z, nu, 102, df2l, xx, ierr)
CC        psi = -z/a*(coff(1,3)*df1l + coff(2,3)*df2l)
CCC       Log derivative
CC        psi = psi/env
CCC       Kinetic energy:
CC        call dbesnu(z, nu, 101, df1l, xx, ierr)
CC        call dbesnu(z, nu, 102, df2l, xx, ierr)
CC        psi = z/a/a*(coff(1,3)*df1l + coff(2,3)*df2l)
CC        call dbesnu(z, nu, 201, df1l, xx, ierr)
CC        call dbesnu(z, nu, 202, df2l, xx, ierr)
CC        psi = (z/a)**2*(coff(1,3)*df1l + coff(2,3)*df2l) + psi
CC        print *, psi,psi/env
CC        stop
C
C      else
C        f1l = exp(kapp*x)
C        f2l = exp(-kapp*x)
C        psi = coff(1,4)*f1l + coff(2,4)*f2l
CC       Log derivative
CC       psi = (kapp*coff(1,4)*f1l - kapp*coff(2,4)*f2l)/psi
CC       First derivative:
CC       psi = kapp*coff(1,4)*f1l - kapp*coff(2,4)*f2l
CC       Kinetic energy:
CC       psi = kapp*kapp*(coff(1,4)*f1l + coff(2,4)*f2l)
C
C      endif
C
C      end

