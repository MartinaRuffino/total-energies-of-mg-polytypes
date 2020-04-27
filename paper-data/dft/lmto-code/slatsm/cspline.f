C#define RX
      subroutine csplint(xa,ya,y2a,n,x,y,yp)
C- Returns a cubic-spline interpolated value at given point
C ----------------------------------------------------------------------
Ci Inputs
Ci   xa    :table of abscissas, n points; must be ordered
Ci   ya    :table of ordinates, n points
Ci   y2a   :second derivatives at tabulated points
Ci   n     :number of points
Ci   x     :point at which to interpolate function
Co Outputs
Co   y     :interpolated function at x
Co   yp    :dy/dx at x
Cr Remarks
Cr   To generate y2a, call cspline
Cu Updates
Cu   18 Sep 09 Adapted from Numerical Recipes Chapter 3
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      double precision x,y,yp,xa(n),y2a(n),ya(n)
C ... Local parameters
      integer khi,klo
      double precision a,b,h

C     Find points in table that bracket x
      call huntx(xa,n,x,-1,klo)
      if (klo == 0) klo = 1
      if (klo == n) klo = n-1
      khi = klo+1
C     Use bisection instead
C      klo = 1
C      khi = n
C    1 if (khi-klo > 1) then
C        k = (khi+klo)/2
C        if (xa(k) > x) then
C          khi = k
C        else
C          klo = k
C        endif
C      goto 1
C      endif
C     klo and khi now bracket the input x.

      h = xa(khi)-xa(klo)
C#ifdef RX
      if (h == 0d0) call rx('csplint: two points same in abscissa')
C#elseC
C      if (h == 0d0) stop 'csplint: two points same in abscissa'
C#endif
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) +
     .  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
      yp = (ya(khi)-ya(klo))/h +
     .  ((3*b**2-1)*y2a(khi)-(3*a**2-1)*y2a(klo))*h/6d0

      end

      subroutine csintg(xa,ya,y2a,n,xmin,xmax,res)
C- Returns integral of a cubic spline between two endpoints
C ----------------------------------------------------------------------
Ci Inputs
Ci   xa    :table of abscissas, n points; must be ordered
Ci   ya    :table of ordinates, n points
Ci   y2a   :second derivatives at tabulated points; see cspline below
Ci   n     :number of points
Ci   xlo   :lower bound of integration
Ci   xhi   :upper bound of integration
Co Outputs
Co   y     :interpolated function at x
Co   yp    :dy/dx at x
Cr Remarks
Cr   To generate y2a, call cspline
Cr   A spline has this property for 2nd derivative y2(x)
Cr   In the interval between points i and i+1,
Cr     y2(x) = a_i(x) ypp(i) + b_i(x) ypp(i+1)
Cr   where
Cr     a_i = [x_i+1-x]/h  and  b_i = [x-x_i]/h  and   h=[x_i+1-x_i]
Cr   Therefore the (constant) third derivative is
Cr     y3 = [ypp(i+1)-ypp(i)]/h
Cu Updates
Cu   27 Feb 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      double precision xa(n),y2a(n),ya(n),xmin,xmax,res
C ... Local parameters
      integer khi,klo
      double precision a,b,h,xlo,xhi,y,yp,ypp,yppp,dx

C     Order limits of integrand
      xlo = min(xmin,xmax)
      xhi = max(xmin,xmax)
      res = 0
      if (xlo == xhi) return

C     Find points in table that bracket lower bound
      call huntx(xa,n,xlo,-1,klo)
      if (klo == 0) klo = 1
      if (xlo == xa(min(klo+1,n))) klo = klo+1
      if (klo == n) klo = n-1

C ... Integrate by panels, current panel from xlo to xlo + deltax
C     xlo,klo are a running lower bound for current panel
C       Cases:        deltax
C     xlo<xa(klo)    xa(klo)-xlo
C     xlo>xa(khi)    xhi-xlo
C     otherwise      xa(khi)-xlo
C     Entry point for new panel
   10 continue
      khi = klo+1
      h = xa(khi)-xa(klo)
      a = (xa(khi)-xlo)/h
      b = (xlo-xa(klo))/h
C#ifdef RX
      if (h == 0d0) call rx('csplint: two points same in abscissa')
C#elseC
C      if (h == 0d0) stop 'csplint: two points same in abscissa'
C#endif

C ... First 3 derivatives at xlo
C     xlo at xa(klo): a=1,b=0
      if (xlo == xa(klo)) then
        y = ya(klo)
        yp = (ya(khi)-ya(klo))/h -
     .       (y2a(khi)+2*y2a(klo))*h/6d0
C     xlo at xa(khi): a=0,b=1
      else if (xlo == xa(khi)) then
        y = ya(khi)
        yp = (ya(khi)-ya(klo))/h +
     .    (2*y2a(khi)+y2a(klo))*h/6d0
C     xlo at some general point
      else
        call csplint(xa,ya,y2a,n,xlo,y,yp)
      endif
      ypp = a*y2a(klo) + b*y2a(khi)
      yppp = (y2a(khi) - y2a(klo))/h

C ... Integrate from xlo to xlo+dx
C       Cases:       upper             update for next panel
C     xlo<xa(klo)    xa(klo)           xlo -> xa(klo)
C     xlo>xa(khi)    xhi               xlo -> xhi
C     otherwise      min(xa(khi),xhi)  xlo->min(xa(khi),xhi), klo->klo+1
      if (xlo < xa(klo)) then
        dx = min(xa(klo),xhi)-xlo
        xlo = min(xa(klo),xhi)
      elseif (xlo >= xa(khi)) then
        dx = xhi-xlo
        xlo = xhi
      else
        dx = min(xa(khi),xhi)-xlo
        xlo = min(xa(khi),xhi)
        klo = min(klo+1,n-1)
      endif
C     Integral of [y + yp dx + 1/2 ypp dx^2 + 1/6 yppp dx^3] is
C                  y dx + 1/2 yp dx^2 + 1/6 ypp dx^3 + 1/24 yppp dx^4
      res = res + y*dx + yp*dx**2/2 + ypp*dx**3/6 + yppp*dx**4/24

      if (xlo < xhi) goto 10

C     Cleanup: if limits of integrand reverse, reverse sign of integral
      if (xmin > xmax) res = -res

      end

      subroutine cspline(x,y,n,yp1,ypn,y2)
C- Return second derivatives of an interpolating spline, from tabulated function
C ----------------------------------------------------------------------
Ci Inputs
Ci   x     :table of abscissas, n points.  x must be ordered
Ci   y     :table of ordinates, n points
Ci   n     :number of points
Ci   yp1   :first derivative of the interpolating function at x(1)
Ci   ypn   :first derivative of the interpolating function at x(n)
Co Outputs
Co   y2    :second derivatives of the interpolating function
Co         :at tabulated points x(1:n)
Cl Local variables
Cl         :
Cr Remarks
Cr   If yp1 and/or ypn are equal to 10^30 or larger, a natural spline
Cr   is generated (zero second derivative on the respective boundary)
Cu Updates
Cu   18 Sep 09 Adapted from Numerical Recipes Chapter 3
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      double precision yp1,ypn,x(n),y(n),y2(n)
C ... Local parameters
      integer i,k
      double precision p,qn,sig,un,u(n)

C     Case natural spline for lower boundary
      if (yp1 >= 0.99d30) then
        y2(1) = 0d0
        u(1) = 0d0
C     Case lower boundary has a specified first derivative
      else
        y2(1) = -0.5d0
        u(1) = (3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do  i = 2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2d0
        y2(i) = (sig-1d0)/p
        u(i) = (6d0*((y(i+1)-y(i))/
     .    (x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/
     .    (x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
C     Case natural spline for upper boundary
      if (ypn > 0.99d30) then
        qn = 0d0
        un = 0d0
C     Case upper boundary has a specified first derivative
      else
        qn = 0.5d0
        un = (3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1d0)
      do  k = n-1, 1, -1
        y2(k) = y2(k)*y2(k+1) + u(k)
      enddo
      end
C      program xspline
CC     driver to test routine spline
C      implicit none
C      integer n
C      double precision pi
C      parameter(n=20)
C      integer i
C      double precision yp1,ypn,x(n),y(n),y2(n)
C
C      write(*,*) 'Second derivatives for sin(x) from 0 to pi'
CC     generate array for interpolation
C      pi = 4*datan(1d0)
C      do  i = 1, 20
C        x(i) = i*pi/N
C        y(i) = dsin(x(i))
C      enddo
C
CC ... Calculate 2nd derivative
C      yp1 = dcos(x(1))
C      ypn = dcos(x(N))
CC     ypn = 1d30
C      call cspline(x,y,n,yp1,ypn,y2)
C
CC     test result
C      write(*,'(t19,a,t35,a,t50,a)') 'spline','actual','diff'
C      write(*,'(t6,a,t17,a,t33,a)') 'angle','2nd deriv','2nd deriv'
C      do  i = 1, n
C        write(*,'(1x,f8.2,2f16.6,1pe15.2)')
C     .    x(i),y2(i),-dsin(x(i)),y2(i)--dsin(x(i))
C      enddo
C      end
C      PROGRAM xsplint
CC     driver for routine splint, which calls spline
C      implicit none
C      INTEGER NP
C      double precision PI
C      PARAMETER(NP=10)
C      INTEGER i,nfunc
C      double precision x,yp1,ypn,xa(NP),ya(NP),y2(NP)
C      double precision f,fp,y,yp
C
C      pi = 4*datan(1d0)
C      do 14 nfunc=1,2
C        if (nfunc == 1) then
C          write(*,*) 'Sine function from 0 to pi'
C          do 11 i=1,np
C            xa(i)=i*pi/np
C            ya(i)=dsin(xa(i))
C   11     continue
C          yp1=dcos(xa(1))
C          ypn=dcos(xa(np))
C        else if (nfunc == 2) then
C          write(*,*) 'Exponential function from 0 to 1'
C          do 12 i=1,np
C            xa(i)=1.0d0*i/np
C            ya(i)=dexp(xa(i))
C   12     continue
C          yp1=dexp(xa(1))
C          ypn=dexp(xa(np))
C        else
C          stop
C        endif
C        call cspline(xa,ya,np,yp1,ypn,y2)
C        write(*,'(1x,t10,a1,t20,a4,t28,a13,t43,a,8x,a,3x,a,4x,a)')
C     .    'x','f(x)','interpolation','diff',
C     .    'f''(x)','interpolation','diff'
C        do 13 i=1,10
C          if (nfunc == 1) then
C            x=(-0.05d0+i/10.0d0)*pi
C            f=dsin(x)
C            fp=dcos(x)
C          else if (nfunc == 2) then
C            x=-0.05d0+i/10.0d0
C            f=dexp(x)
C            fp = f
C          endif
C          call csplint(xa,ya,y2,np,x,y,yp)
C          write(*,'(1x,3f12.6,1pe12.2,0p,2f12.6,1pe12.2)')
C     .      x,f,y,f-y,fp,yp,fp-yp
C   13   continue
C        write(*,*) '-----------------------------------'
CC       write(*,*) 'Press RETURN'
CC       read(*,*)
C   14 continue
C      end
C      subroutine rx(string)
C      character *(*) string
C
C      print *, string
C      stop
C      end
