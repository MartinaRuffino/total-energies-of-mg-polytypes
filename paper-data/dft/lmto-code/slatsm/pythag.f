      double precision function pythag(a,b)
      implicit none
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u,d1mach
      p = dmax1(dabs(a),dabs(b))
      if (p == 0.0d0) go to 20
      if (p <= d1mach(1)) then
        p = 0d0
        goto 20
      endif

      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t == 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C      subroutine fmain
C      implicit none
C      double precision f,pythag,a,b,d1mach
C      ;
C
C      a = 0d0
C      b = d1mach(1)
CC      b = 1.4804388770234955d-11/2.225073858507201D-8*b
C
C      f = pythag(0d0,b)
C      print *, f
C      end
