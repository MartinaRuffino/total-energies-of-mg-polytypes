      subroutine radmwt(opt,rmax,a,nr,rofi,wt)
C- Makes mesh and weights for numerical integration on shifted log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 for uniform weight, i.e. int dr f(r)
Ci         :1 for r^2 weight,     i.e. int dr r^2 f(r)
Ci         :10s digit
Ci         :0 for 3-point quadrature (Simpson's rule)
Ci         :1 for 5-point quadrature (NOT IMPLEMENTED)
Ci         :2 for 7-point quadrature
Ci         :See Remarks for special treatment of endpoints
Ci         :100s digit
Ci         :0 return both rofi and wt
Ci         :1 return both wt only; rofi is not touched
Ci         :1000s digit
Ci         :0 shifted log mesh
Ci         :1 standard log mesh
Ci   rmax  :augmentation radius, in a.u.
Ci   a     :Mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Co Outputs
Co   rofi  :radial mesh points for shifted log mesh
Co         :Note: rofi is not touched if 100s digit opt is nonzero
Co   wt    :mesh weights; see Remarks
Cl Local variables
Cl         :
Cr Remarks
Cr    opt     integral    rule
Cr     0       f(r)       3-point quadrature (Simpson's rule)
Cr     1       f(r)*r^2   3-point quadrature (Simpson's rule)
Cr    20       f(r)       7-point quadrature
Cr    21       f(r)*r^2   7-point quadrature
Cr  Thus sum_i wt_i approx rmax (opt=0), or rmax^3/3 (opt=1)
Cr
Cr  Shifted log mesh:
Cr    r(i) = b[exp(a*(i-1))-1]
Cr    b is fixed so r(nr) = rmax.  Note also r(1)=0
Cr  Let I = int_0^rmax f[r(i)] dr
Cr    I = int_1^nr f[r(i)] dr/di di = int_1^nr f[r(i)] (r(i)+b)a di
Cr  Replace exact integral with quadrature:
Cr    I = sum_i=1:nr f[r(i)] w(i) with  w(i) = w0(i) * [r(i)+b] * a
Cr    w0(i) are weights for Newton-Cotes (uniformly spaced) quadrature
Cr    Note: if sum_i g(i) w0(i) is exact for g(i) = i^p (p is some power)
Cr    then: for f = i^p/(r(i)+b), I=a*[nr^(p+1)-1]/(p+1) is integrated exactly
Cr    3-point (Simpson) rule: exact for p=0,1,2
Cr    7-point rule: exact for p=0..7
Cr    I is NOT exact for f=constant.
Cr
Cr  Endpoints:
Cr  The 3-point (or 7-point) rule is used except at the beginning of
Cr  the mesh.  Special treatment is given there, to synchronize the
Cr  given mesh to the number of points the quadrature requires:
Cr    3-point rule requires 2n+1 points (odd number)
Cr    7-point rule rule requires 6n+1 points
Cr
Cr  The starting point is adjusted to meet this requirement, and points
Cr  up to the adjusted starting point are treated as follows:
Cr  For the 3-point rule:
Cr      (nr%2)      initial points          3-point mesh begins at
Cr         1        no special treatment    1
Cr         0        4-point                 4
Cr  For the 7-point rule:
Cr     (nr-1%6)+1   initial points          7-point mesh begins at
Cr         1        no special treatment    1
Cr         2        4-point+5-point         8
Cr         3        3-point                 3
Cr         4        4-point                 4
Cr         5        5-point                 5
Cr         6        4-point+3-point         6
Cr    where
Cr    4-point is "Simpson's 3/8 rule", wt=(1,3,3,1)/8*3
Cr    5-point is "Boole's rule", wt=(7,32,12,32,7)/90*4
Cu Updates
Cu   04 Apr 12 Repackaged, and extended for 7-point quadrature
Cu   22 Feb 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,opt
      double precision rmax,a,wt(nr)
      real(8),target :: rofi(nr)
C ... Local parameters
      integer ir,i07,opt0,opt1,opt2,opt3,i1
      double precision fac,b,wtl(0:6),xx
      real(8),target :: rofil(nr)
      real(8),pointer:: rloc(:)
      real(8), parameter :: wt5(0:4) = [7d0,32d0,12d0,32d0,7d0]
      real(8), parameter :: wt7(0:6) = [41d0,216d0,27d0,272d0,27d0,216d0,41d0]

      if (nr < 3) call rx('radmwt: nr is < 3')
      call dpzero(wt,nr)
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      if (opt1 /= 0 .and. opt1 /= 2)
     .  call rxi('radmwt: not implemented, 10s digit opt=',opt1)
      opt2 = mod(opt/100,10)
      opt3 = mod(opt/1000,10)
C      call isanrg(opt3,0,1,'radmwt: ','mesh type',.true.)

      i07 = mod(nr-1,6) + 1     ! First point for 7-point rule
      if (opt1 == 0) i07 = 0

C ... Case nr is even: first points with 4-point rule
      i1 = 1
      if (mod(nr,2) == 0) then
        wt(1) = 1d0/8d0*3
        wt(2) = 3d0/8d0*3
        wt(3) = 3d0/8d0*3
        wt(4) = 1d0/8d0*3
        i1 = i1+3
C       print *, '4 point'
      endif

C ... Next points with 3-point rule
      if (i07 == 3 .or. i07 == 6) then
        wt(i1+0) = wt(i1+0) + 1d0/3d0
        wt(i1+1) = wt(i1+1) + 4d0/3d0
        wt(i1+2) = wt(i1+2) + 1d0/3d0
        i1 = i1+2
C       print *, '3 point'
      endif

C ... Next points with 4-point rule ... no such cases
C     if (i07 == 3 .or. i07 == 6) then
C     endif

C ... Next points with 5-point rule
      if (i07 == 2 .or. i07 == 5) then
        wt(i1+0) = wt(i1+0) + wt5(0)/90*4
        wt(i1+1) = wt(i1+1) + wt5(1)/90*4
        wt(i1+2) = wt(i1+2) + wt5(2)/90*4
        wt(i1+3) = wt(i1+3) + wt5(3)/90*4
        wt(i1+4) = wt(i1+4) + wt5(4)/90*4
        i1 = i1+4
C       print *, '5 point'
      endif

C ... Simpson rule (3 point rule) for i1:nr
      if (opt1 == 0) then
        if (i1 < nr) then
          wt(i1) = wt(i1) + 1d0/3d0
          wt(nr) = wt(nr) + 1d0/3d0
        endif
        do  ir = i1+2, nr-2, 2
          wt(ir) = 2d0/3d0 + wt(ir)
        enddo
        do  ir = i1+1, nr, 2
          wt(ir) = 4d0/3d0 + wt(ir)
        enddo
C       print *, 'simpson'
C ... 7 point rule i1:nr
      else
        call dpcopy(wt7,wtl,1,7,6d0/840d0)
        do  ir = i1, nr-6, 6
          wt(ir+0) = wt(ir+0) + wtl(0)
          wt(ir+1) = wt(ir+1) + wtl(1)
          wt(ir+2) = wt(ir+2) + wtl(2)
          wt(ir+3) = wt(ir+3) + wtl(3)
          wt(ir+4) = wt(ir+4) + wtl(4)
          wt(ir+5) = wt(ir+5) + wtl(5)
          wt(ir+6) = wt(ir+6) + wtl(6)
        enddo
C       print *, '7 point'
      endif

C ... Generate radial mesh: either local array, or passed one
      rloc => rofi
      if (opt2 /= 0) rloc => rofil
      if (opt3 == 0) then
        call radmsh(rmax,a,nr,rloc)
      else
        xx = rmax
        fac = exp(-a)
        do  ir = nr, 1, -1
          rloc(ir) = xx
          xx = xx*fac
        enddo
      endif

C ... Jacobian for shifted log mesh
      if (opt3 == 0) then
        b = rmax / (dexp(a*nr-a)-1d0)
        do  ir = 1, nr
          wt(ir) = wt(ir) * a*(rloc(ir)+b)
          if (opt0 == 1) wt(ir) = wt(ir)*rloc(ir)**2
        enddo
      else
        do  ir = 1, nr
          wt(ir) = wt(ir) * a*rloc(ir)
          if (opt0 == 1) wt(ir) = wt(ir)*rloc(ir)**2
        enddo
      endif
      end

      subroutine radmwtp(opt,rofi,a,i1,nr,wt)
C- Mesh and weights for integration on portion of shifted log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 for uniform weight, i.e. int dr f(r)
Ci         :1 for r^2 weight,     i.e. int dr r^2 f(r)
Ci         :10s digit
Ci         :0 for 3-point quadrature (Simpson's rule)
Ci         :1 for 5-point quadrature (NOT IMPLEMENTED)
Ci         :2 for 7-point quadrature
Ci         :See Remarks for special treatment of endpoints
Ci   rofi  :radial mesh points for shifted log mesh; see radmsh
Ci   a     :Mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   i1    :integration weights for region rofi(i1:nr)
Ci   nr    :number of radial mesh points
Co Outputs
Co   wt    :mesh weights; see Remarks
Cl Local variables
Cl         :
Cr Remarks
Cr   This routine is similar to radmwt except that it is designed
Cr   to integrate function on a shifted log mesh, starting at
Cr   some intermediate point rofi(i1).
Cr
Cr  Endpoints:
Cr  The 3-point (or 7-point) rule is used except at the end of
Cr  the mesh.  Special treatment is given there, to synchronize the
Cr  given mesh to the number of points the quadrature requires:
Cr    3-point rule requires 2n+1 points (odd number)
Cr    7-point rule rule requires 6n+1 points
Cr
Cr  The end point is adjusted to meet this requirement, and points
Cr  between this point and the true endpoint are treated as follows:
Cr  For the 3-point rule:
Cr      (nr-i1%2)   last points             3-point mesh ends at
Cr         0        no special treatment    nr
Cr         1        4-point                 nr-2
Cr  For the 7-point rule:
Cr     (nr-i1%6)+1  last points             7-point mesh ends at
Cr         1        no special treatment    nr-0
Cr         2        4-point+5-point         nr-7
Cr         3        3-point                 nr-2
Cr         4        4-point                 nr-3
Cr         5        5-point                 nr-4
Cr         6        4-point+3-point         nr-5
Cr    where
Cr    4-point is "Simpson's 3/8 rule", wt=(1,3,3,1)/8*3
Cr    5-point is "Boole's rule", wt=(7,32,12,32,7)/90*4
Cu Updates
Cu   04 Apr 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer i1,nr,opt
      double precision a,rofi(nr),wt(nr)
C ... Local parameters
      integer ir,opt0,opt1,opt2,i07,ie
      double precision b,wtl(0:6)
      real(8), parameter :: wt5(0:4) = [7d0,32d0,12d0,32d0,7d0]
      real(8), parameter :: wt7(0:6) = [41d0,216d0,27d0,272d0,27d0,216d0,41d0]

      if (nr-i1+1 < 3) call rx('radmwtp: nr is < 3')
      call dpzero(wt,nr)
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      if (opt1 /= 0 .and. opt1 /= 2)
     .  call rxi('radmwtp: not implemented, 10s digit opt=',opt1)
      opt2 = mod(opt/100,10)

      i07 = mod(nr-i1,6) + 1    ! First point for 7-point rule
      if (opt1 == 0) i07 = 0
      ie = nr                   ! Current index to last point for which wt needed

C ... Last points with 3-point rule
      if (i07 == 3 .or. i07 == 6) then
        wt(ie-0) = wt(ie-0) + 1d0/3d0
        wt(ie-1) = wt(ie-1) + 4d0/3d0
        wt(ie-2) = wt(ie-2) + 1d0/3d0
        ie = ie-2
C       print *, '3 point'
      endif

C ... Case nr is even: last points with 4-point rule
      if (mod(nr-i1,2) == 1) then
        wt(ie-0) = wt(ie-0) + 1d0/8d0*3
        wt(ie-1) = wt(ie-1) + 3d0/8d0*3
        wt(ie-2) = wt(ie-2) + 3d0/8d0*3
        wt(ie-3) = wt(ie-3) + 1d0/8d0*3
        ie = ie-3
C       print *, '4 point'
      endif

C ... Last points with 4-point rule ... no more cases
C     if (i07 == 3 .or. i07 == 6) then
C     endif

C ... Last points with 5-point rule
      if (i07 == 2 .or. i07 == 5) then
        wt(ie-0) = wt(ie-0) + wt5(0)/90*4
        wt(ie-1) = wt(ie-1) + wt5(1)/90*4
        wt(ie-2) = wt(ie-2) + wt5(2)/90*4
        wt(ie-3) = wt(ie-3) + wt5(3)/90*4
        wt(ie-4) = wt(ie-4) + wt5(4)/90*4
        ie = ie-4
C       print *, '5 point'
      endif

C ... Simpson rule (3 point rule) for i1:ie
      if (opt1 == 0) then
        if (i1 < ie) then
          wt(i1) = wt(i1) + 1d0/3d0
          wt(ie) = wt(ie) + 1d0/3d0
        endif
        do  ir = i1+2, ie-2, 2
          wt(ir) = 2d0/3d0 + wt(ir)
        enddo
        do  ir = i1+1, ie, 2
          wt(ir) = 4d0/3d0 + wt(ir)
        enddo
C       print *, 'simpson'
C ... 7 point rule i1:ie
      else
        call dpcopy(wt7,wtl,1,7,6d0/840d0)
        do  ir = i1, ie-6, 6
          wt(ir+0) = wt(ir+0) + wtl(0)
          wt(ir+1) = wt(ir+1) + wtl(1)
          wt(ir+2) = wt(ir+2) + wtl(2)
          wt(ir+3) = wt(ir+3) + wtl(3)
          wt(ir+4) = wt(ir+4) + wtl(4)
          wt(ir+5) = wt(ir+5) + wtl(5)
          wt(ir+6) = wt(ir+6) + wtl(6)
        enddo
C       print *, '7 point'
      endif

C ... Jacobian for shifted log mesh
      b = rofi(nr) / (dexp(a*nr-a)-1d0)
      do  ir = i1, nr
        wt(ir) = wt(ir) * a*(rofi(ir)+b)
        if (opt0 == 1) wt(ir) = wt(ir)*rofi(ir)**2
      enddo

      end

      subroutine radprm(nr,rofi,a,b)
C- Deduces shifted mesh parameters from given radial mesh
C Note: code not fully debugged!
      implicit none
      integer nr
      double precision a,b,rofi(nr)
      double precision rmax,aold
      integer ir0,imid
      logical have0

      have0 = rofi(1) == 0
      ir0 = 1
      if (.not. have0) ir0 = 0
      rmax = rofi(nr)

C     First guess for a,b
      a = dlog(rofi(nr)/rofi(nr-1))
      imid = nr/2
C     if (nr > 100) imid = nr-20
      b = rofi(imid)/(exp(a*(imid-ir0))-1)

C     print *, 1.6926912750585881d-05 * (exp(a*(nr-ir0))-1)

C     Iterate until converged
C     b from 2nd point, assuming a is correct
   10 continue
      b = rofi(imid)/(exp(a*(imid-ir0))-1)
      aold = a
      a = dlog(rmax/b+1) / (nr-ir0)

C      print *, a-aold
C     check for convergence
      if (abs(a-aold) > 1d-8) goto 10

C      print *, a,b
C      do  ir = 2, nr+1-ir0
C       print 333, ir, rofi(ir+ir0-1), b*(exp((ir+ir0-1)*a)-1),
C     .    rofi(ir+ir0-1)-b*(exp((ir+ir0-1)*a)-1)
C  333  format(i4,2f15.10,1pe10.3)
C      enddo
      end

      subroutine radsum(nrx,nr,nlml,nsp,wt,rho,sum)
C- Numerical integration of a function on a shifted log mesh
Cu   19 Jun 00 added extra arguments
      implicit none
      integer nrx,nr,nlml,nsp
      double precision wt(nr),rho(nrx,nlml,nsp),sum,ddot

      sum = ddot(nr,wt,1,rho,1)
      if (nsp == 2) sum = sum + ddot(nr,wt,1,rho(1,1,2),1)
      end

      subroutine radext(mode,nr,nrx,fac,a,rmax,nrbig,rbig,rofi,rwgt)
C- Find radius, mesh suitable for extending orbitals outside MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0  nrbig,rbig are input; do not make them
Ci         :1  set rbig = smaller of   rofi(nrx)  and  fac*rmax
Ci         :10s digit
Ci         :if nonzero, make rofi and rwgt
Ci   nr    :number of radial mesh points on regular mesh
Ci   nrx   :maximum allowed number of radial mesh points
Ci   fac   :approximate factor to scale rmax, rbig ~ fac*rmax
Ci         :NB: true factor is constrained because rbig must
Ci         :conform to radial mesh specified by (rmax,a,nr)
Ci   a     :mesh points are given by
Ci         :rofi(i) = rmax [e^(a(i-1))-1] / [e^(a(nr-1))-1]
Ci   rmax  :augmentation radius, in a.u.,
Co Outputs
Cio  nrbig :number of points on extended mesh.
Cio        :NB: nrbig is input if 1s digit mode=0
Cio        :In the latter case, nrbig must be consistent with the mesh
Cio        :points specified by (a,nr,rmax) and also rbig.
Cio  rbig  :sphere radius of extended mesh
Cio        :NB: rbig is input if 1s digit mode=0
Co   rofi  :(10s digit mode > 0)
Co         :radial mesh points: rofi(1..nrbig) will be generated
Co         :rofi(nrbig) is rmax for extended mesh
Co   rwgt  :(10s digit mode > 0)
Co         :radial mesh weights: rwgt(1..nrbig) will be generated
Co         :rwgt is actually designed for two integration radii:
Co         :int(0,rmax) = I(1..nr) and int(rmax,rbig) = I(nr..nrbig).
Co         :Integral int(1..nrbig) must be done in two steps, by summing
Co         :I(1..nr) and I(nr..nrbig)
Cr Remarks
Cr
Cu Updates
Cu   24 Sep 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nr,nrx,nrbig
      double precision rmax,fac,rbig,a,rofi(*),rwgt(*)
C ... Local parameters
      integer idn

      if (mod(mode,10) == 1) then
      rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
C     If rbig>fac*rmax, estimate from exp((nrbig-nr)a) = fac
      if (rbig > fac*rmax) then
        idn = dlog(fac)/a
        if (mod(idn,2) == 1) idn = idn-1
        nrbig = min(nr+idn,nrx)
        rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
      endif
      endif

C --- Points and weights on extended mesh ---
      if (mod(mode/10,10) /= 0) then
        call radmsh(rbig,a,nrbig,rofi)
        call radwgt(0,rbig,a,nrbig,rwgt)
        if (nr < nrbig) rwgt(nr) = rwgt(nr)/2
      endif
      end
C      subroutine fmain
CC- Tests radmsh
CC  fcd radmsh.f; lk radmsh.o prmx.o prrmsh.o
C      implicit none
C      integer nr,ir,k,p,opt
C      double precision a
C
C                                          ! err for:
C                                          ! p=10  int(1)
C
CC     Integrates ir**p/(r+b) exactly for p<=3
CC     parameter(nr=4,opt=0,a=.1d0)        ! 7e-3  3e-6
CC     parameter(nr=3,opt=0,a=.1d0)        ! 2e-3  2e-6
CC     parameter(nr=7,opt=0,a=.1d0)        ! 1e-3  1e-6
CC     parameter(nr=20,opt=0,a=.1d0)       ! 5e-5  1e-6
CC     parameter(nr=21,opt=0,a=.1d0)       ! 4e-5  1e-6
C
CC     Integrates ir**p/(r+b) exactly for p<=7
CC     parameter(nr=7,opt=20,a=.1d0/2)     ! 3e-5  1e-13
CC     parameter(nr=13,opt=20,a=.1d0)      ! 1e-6  3e-11
C
CC     Integrates ir**p/(r+b) exactly for p<=5
CC     parameter(nr=5,opt=20,a=.1d0)       ! 1e-3  5e-9
CC     parameter(nr=11,opt=20,a=.1d0)      ! 3e-6  2e-9
CC     parameter(nr=11,opt=20,a=.1d0/2)    ! 2e-6  3e-11
C
CC     Integrates ir**p/(r+b) exactly for p<=3
CC     parameter(nr=6,opt=20,a=.1d0)       ! 2e-3  2e-6
CC     parameter(nr=6,opt=20,a=.1d0/2)     ! 8e-4  1e-7
CC     parameter(nr=8,opt=20,a=.1d0/2)     ! 8e-5  8e-8
C
CC     parameter(nr=9,opt=20,a=.1d0/2)     ! 6e-6  2e-8
CC     parameter(nr=10,opt=20,a=.1d0/2)    ! 3e-6  6e-8
C      parameter(nr=12,opt=20,a=.1d0/2)    ! 2e-6  6e-8
CC     parameter(nr=14,opt=20,a=.1d0/2)    ! 6e-7  3e-8
CC     parameter(nr=15,opt=20,a=.1d0/2)    ! 2e-7  9e-9
CC     parameter(nr=16,opt=20,a=.1d0/2)    ! 1e-7  3e-8
CC     parameter(nr=4*6+1,opt=20,a=.1d0/2) ! 5e-9  1e-13
CC     parameter(nr=6*6+1,opt=20,a=.1d0/2) ! 3e-10 1e-13
CC     parameter(nr=6*6+3,opt=20,a=.1d0/2) ! 2e-10 2e-9
C
CC     parameter(nr=60*6+6,opt=1020,a=.1d0/2) ! 2e-10 2e-9
C
CC     parameter(nr=211,opt=21,a=.1d0/4)
C
C      integer i0
C      parameter (i0=0)
C      double precision rmax,b,fn(i0+nr,12),res(12),Iexact,xx
C      double precision rofi(i0+nr),wt(i0+nr) !,wt2(nr)
C
C      rmax = 2.5d0
C      b = rmax / (dexp(a*(i0+nr)-a)-1d0)
C      if (mod(opt/1000,10) == 1) b = 0
C      call info5(0,1,1,' Numerical integration on shifted'//
C     .  ' log mesh:  opt=%i  nr=%i  rmax=%d  a=%d  b=%;10g',
C     .  opt,i0+nr,rmax,a,b)
C
C      if (i0 == 0) then
C        call radmwt(opt,rmax,a,nr,rofi,wt)
C      else
C        call radmsh(rmax,a,i0+nr,rofi)
C        call radmwtp(opt,rofi,a,i0+1,i0+nr,wt)
C      endif
C
CC      call radwgt(intopt,rmax,a,i0+nr,wt2)
CC       do  ir = 1, nr
CC         print 333, ir, rofi(ir),wt(ir),wt2(ir),wt(ir)-wt2(ir)
CC  333    format(i3,f15.10,2f20.15,1pe10.2)
CC       enddo
C
CC     call prrmsh('wt',rofi,wt,nr,nr,1)
C
C      call dpzero(res,12)
C      do  ir = 1, i0+nr
C        fn(ir,1) = wt(ir)
CC       Integral of ir**power/(r+b) for power=0..10
C        res(1) = res(1) + wt(ir)
C        do  k = 2, 12
C          fn(ir,k) = 1/(rofi(ir)+b)*dble(ir)**(k-2)
C          res(k) = res(k) + wt(ir)*fn(ir,k)/dble(i0+nr)**(k-2)
C        enddo
C      enddo
C      xx = rmax - rofi(i0+1)
C      call info5(0,0,1,' I for f = const:%30p %;20,15D  exact:'//
C     .  '%;20,15D  diff %,2;g',res(1),xx,res(1)-xx,0,0)
C      do  k = 2, 12
C        p = k-2
C        Iexact = a*(dble(i0+nr)**(p+1)-(i0+1d0)**(p+1))/
C     .           (p+1)/dble(i0+nr)**p
CC       Iexact = a*(dble(i0+nr)**(p+1)-3**(p+1))/(p+1)/dble(i0+nr)**p
CC       Iexact = a*(3**(p+1)-1**(p+1))/(p+1)/dble(i0+nr)**p
C        call info5(0,0,0,'   for f = ir^%i/(r+b):%30p %;20,15D  exact:'
C     .    //'%;20,15D  diff %,2;2g',k-2,res(k),Iexact,res(k)-Iexact,0)
C      enddo
C
CC     call prrmsh('f',rofi,fn,i0+nr,i0+nr,10)
C
C      end
