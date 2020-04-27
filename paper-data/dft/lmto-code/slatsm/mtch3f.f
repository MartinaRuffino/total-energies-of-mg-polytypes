      subroutine mtch2f(mode,lmax,v0,v1,v2,s0,s1,s2,v,s,c1,c2)
C- Matches sum of 2 functions to two conditions, for each l
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0  Total function is v0 + c1*v1 + c2*v2
Ci         :1  Total function is c1*v1 + c2*v2
Ci         :   In this mode v0 and s0 are not used
Ci   v0    :   value of function 0, for 1st condition
Ci   v1    :   value of function 1, for 1st condition
Ci   v2    :   value of function 2, for 1st condition
Ci   s0    :   value of function 0, for 2nd condition
Ci   s1    :   value of function 1, for 2nd condition
Ci   s2    :   value of function 2, for 2nd condition
Ci   v     :   1st condition function must satisfy
Ci   s     :   2nd condition function must satisfy
Co Outputs
Co   c1    :   coefficient to 1st function
Co   c2    :   coefficient to 2nd function
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 May 11 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmax
      double precision v0(0:lmax),v1(0:lmax),v2(0:lmax),
     .                 s0(0:lmax),s1(0:lmax),s2(0:lmax),
     .                 v(0:lmax),s(0:lmax),c1(0:lmax),c2(0:lmax)
C ... Local parameters
      integer l
      double precision det

      if (mod(mode,10) == 0) then
        do l = 0, lmax
          det = v1(l)*s2(l) - s1(l)*v2(l)
          c1(l) = (s2(l)*(v(l)-v0(l)) - v2(l)*(s(l)-s0(l)))/det
          c2(l) =-(s1(l)*(v(l)-v0(l)) - v1(l)*(s(l)-s0(l)))/det
        enddo
      else
        do l = 0, lmax
          det = v1(l)*s2(l) - s1(l)*v2(l)
          c1(l) = (s2(l)*v(l) - v2(l)*s(l))/det
          c2(l) =-(s1(l)*v(l) - v1(l)*s(l))/det
        enddo
      endif

      end

      subroutine mtch3f(mode,lmax,v0,v1,v2,v3,s0,s1,s2,s3,t0,t1,t2,t3,
     .  v,s,t,nc,c123)
C- Matches sum of three functions to three conditions, for each l
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0  Total function is f0 + c1*f1 + c2*f2 + c3*f3
Ci         :1  Total function is c1*v1 + c2*v2
Ci         :10s digit
Ci         :0  c123 saved as c123(1:3,0:lmax)
Ci         :1  c123 saved as c123(0:lmax,1:3)
Ci         :   In this mode v0 and s0 are not used
Ci   v0    :   value of function 0, for 1st condition (used if mode=1)
Ci   v1    :   value of function 1, for 1st condition
Ci   v2    :   value of function 2, for 1st condition
Ci   v3    :   value of function 3, for 1st condition
Ci   s0    :   value of function 0, for 2nd condition (used if mode=1)
Ci   s1    :   value of function 1, for 2nd condition
Ci   s2    :   value of function 2, for 2nd condition
Ci   s3    :   value of function 3, for 2nd condition
Ci   t0    :   value of function 0, for 3rd condition (used if mode=1)
Ci   t1    :   value of function 1, for 3rd condition
Ci   t2    :   value of function 2, for 3rd condition
Ci   t3    :   value of function 3, for 3rd condition
Ci   v     :   1st condition function must satisfy
Ci   s     :   2nd condition function must satisfy
Ci   t     :   3rd condition function must satisfy
Co Outputs
Co   c123  :   coefficients to functions 1,2,3
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 May 11 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmax,nc
      double precision v0(0:lmax),v1(0:lmax),v2(0:lmax),v3(0:lmax),
     .                 s0(0:lmax),s1(0:lmax),s2(0:lmax),s3(0:lmax),
     .                 t0(0:lmax),t1(0:lmax),t2(0:lmax),t3(0:lmax),
     .                 v(0:lmax),s(0:lmax),t(0:lmax),c123(nc,3)
C ... Local parameters
      integer l,mode0,mode1
      double precision norm(3,3),rhs(3),inorm(3,3),det,res(3)

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      do l = 0, lmax
        norm(1,1) = v1(l)
        norm(1,2) = v2(l)
        norm(1,3) = v3(l)
        norm(2,1) = s1(l)
        norm(2,2) = s2(l)
        norm(2,3) = s3(l)
        norm(3,1) = t1(l)
        norm(3,2) = t2(l)
        norm(3,3) = t3(l)
C       call prmx('norm',norm,3,3,3)
        call dinv33(norm,0,inorm,det)
        if (mode0 == 0) then
          rhs(1) = v(l)-v0(l)
          rhs(2) = s(l)-s0(l)
          rhs(3) = t(l)-t0(l)
C         call prmx('rhs',rhs,3,3,1)
        else
          rhs(1) = v(l)
          rhs(2) = s(l)
          rhs(3) = t(l)
        endif
        if (mode1 == 0) then
          call dgemm('N','N',3,1,3,1d0,inorm,3,rhs,3,0d0,c123(1,l+1),nc)
        else
          call dgemm('N','N',3,1,3,1d0,inorm,3,rhs,3,0d0,res,3)
          c123(l+1,1) = res(1)
          c123(l+1,2) = res(2)
          c123(l+1,3) = res(3)
        endif
      enddo

      end
C      subroutine fmain
CC     test: fit log(x)  to x/2 + c1*cos(x) + c2*sin(x) + c3*exp(-x/2) l = 0
CC               cos(x)  to x/2 + c1*log(x) + c2*sin(x) + c3*exp(-x/2) l = 1
CC               sin(x)  to x/2 + c1*cos(x) + c2*log(x) + c3*exp(-x/2) l = 2
CC            exp(-x/2)  to x/2 + c1*cos(x) + c2*sin(x) + c3*log(x)    l = 3
CC     Fit conditions:
CC     Use value at x=1, slope at x=1, value at x=2
C      implicit none
C      integer lmax,l,i,nc
C      parameter (lmax=3,nc=7)
C      double precision v0(0:lmax),v1(0:lmax),v2(0:lmax),v3(0:lmax),
C     .                 s0(0:lmax),s1(0:lmax),s2(0:lmax),s3(0:lmax),
C     .                 t0(0:lmax),t1(0:lmax),t2(0:lmax),t3(0:lmax),
C     .                 v(0:lmax),s(0:lmax),t(0:lmax),c123(nc,lmax+1)
C      double precision x,xc,xs,xe,xl,x1,x2
C
C      x1 = 1
C      x2 = 2
C
C      do l = 0, lmax
C        v0(l) = x1/2
C        s0(l) = 1d0/2
C        t0(l) = x2/2
C      enddo
C      xc = dcos(x1)
C      xs = dsin(x1)
C      xe = dexp(-x1/2)
C      xl = dlog(x1)
C      v1(0) = xc
C      v2(0) = xs
C      v3(0) = xe
C      v(0)  = xl
C      v1(1) = xl
C      v2(1) = xs
C      v3(1) = xe
C      v(1)  = xc
C      v1(2) = xc
C      v2(2) = xl
C      v3(2) = xe
C      v(2)  = xs
C      v1(3) = xc
C      v2(3) = xs
C      v3(3) = xl
C      v(3)  = xe
C
C      xc = -dsin(x1)
C      xs = dcos(x1)
C      xe = -dexp(-x1/2)/2
C      xl = 1/x1
C      s1(0) = xc
C      s2(0) = xs
C      s3(0) = xe
C      s(0)  = xl
C      s1(1) = xl
C      s2(1) = xs
C      s3(1) = xe
C      s(1)  = xc
C      s1(2) = xc
C      s2(2) = xl
C      s3(2) = xe
C      s(2)  = xs
C      s1(3) = xc
C      s2(3) = xs
C      s3(3) = xl
C      s(3)  = xe
C
C      xc = dcos(x2)
C      xs = dsin(x2)
C      xe = dexp(-x2/2)
C      xl = dlog(x2)
C      t1(0) = xc
C      t2(0) = xs
C      t3(0) = xe
C      t(0)  = xl
C      t1(1) = xl
C      t2(1) = xs
C      t3(1) = xe
C      t(1)  = xc
C      t1(2) = xc
C      t2(2) = xl
C      t3(2) = xe
C      t(2)  = xs
C      t1(3) = xc
C      t2(3) = xs
C      t3(3) = xl
C      t(3)  = xe
C
C      call mtch3f(0,3,v0,v1,v2,v3,s0,s1,s2,s3,t0,t1,t2,t3,
C     .  v,s,t,nc,c123)
C      print *, 'writing x/2 + c1 f1 + c2 f2 + c3 f3 to file fort.88'
C      x = .5d0
C      do i = 1, 25
C        xc = dcos(x)
C        xs = dsin(x)
C        xe = dexp(-x/2)
C        xl = dlog(x)
CC        res(i,1) = x
CC        res(i,2) = xl
CC        res(i,3) = xc
CC        res(i,4) = xs
CC        res(i,5) = xe
C
C        write(88,"(f8.4,8f12.6)") x,
C     .    xl, x/2 + c123(1,1)*xc + c123(2,1)*xs + c123(3,1)*xe,
C     .    xc, x/2 + c123(1,2)*xl + c123(2,2)*xs + c123(3,2)*xe,
C     .    xs, x/2 + c123(1,3)*xc + c123(2,3)*xl + c123(3,3)*xe,
C     .    xe, x/2 + c123(1,4)*xc + c123(2,4)*xs + c123(3,4)*xl
C        x = x+.1d0
C      enddo
C
C      call mtch3f(1,3,v0,v1,v2,v3,s0,s1,s2,s3,t0,t1,t2,t3,
C     .  v,s,t,nc,c123)
CC     call prmx('c123',c123,nc,3,lmax+1)
C      print *, 'writing c1 f1 + c2 f2 + c3 f3 to file fort.89'
C      x = .5d0
C      do i = 1, 25
C        xc = dcos(x)
C        xs = dsin(x)
C        xe = dexp(-x/2)
C        xl = dlog(x)
CC        res(i,1) = x
CC        res(i,2) = xl
CC        res(i,3) = xc
CC        res(i,4) = xs
CC        res(i,5) = xe
C
C        write(89,"(f8.4,8f12.6)") x,
C     .    xl, 0*x/2 + c123(1,1)*xc + c123(2,1)*xs + c123(3,1)*xe,
C     .    xc, 0*x/2 + c123(1,2)*xl + c123(2,2)*xs + c123(3,2)*xe,
C     .    xs, 0*x/2 + c123(1,3)*xc + c123(2,3)*xl + c123(3,3)*xe,
C     .    xe, 0*x/2 + c123(1,4)*xc + c123(2,4)*xs + c123(3,4)*xl
C        x = x+.1d0
C      enddo
C
C      call mtch3f(11,3,v0,v1,v2,v3,s0,s1,s2,s3,t0,t1,t2,t3,
C     .  v,s,t,nc,c123)
CC     call prmx('c123',c123,nc,lmax+1,3)
C      print *, 'redo fit but generate c in transpose order'
C      print *, 'write to fort.90; should match fort.89'
C      x = .5d0
C      do i = 1, 25
C        xc = dcos(x)
C        xs = dsin(x)
C        xe = dexp(-x/2)
C        xl = dlog(x)
C        write(90,"(f8.4,8f12.6)") x,
C     .    xl, 0*x/2 + c123(1,1)*xc + c123(1,2)*xs + c123(1,3)*xe,
C     .    xc, 0*x/2 + c123(2,1)*xl + c123(2,2)*xs + c123(2,3)*xe,
C     .    xs, 0*x/2 + c123(3,1)*xc + c123(3,2)*xl + c123(3,3)*xe,
C     .    xe, 0*x/2 + c123(4,1)*xc + c123(4,2)*xs + c123(4,3)*xl
C        x = x+.1d0
C      enddo
C
C      end
