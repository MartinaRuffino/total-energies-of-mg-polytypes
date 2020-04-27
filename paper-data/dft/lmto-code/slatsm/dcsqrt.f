      double complex function dcsqrt (z)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
      implicit none
      double complex z

      double precision cdabs,dsqrt,xtmp,ytmp,x,y,r
c
      x = dble  (z)
      y = dimag (z)
      r = cdabs (z)
c
      if (r == 0.d0) dcsqrt = (0.0d0, 0.0d0)
      if (r == 0.d0) return
c
      xtmp = dsqrt ((r+dabs(x))*0.5d0)
      ytmp = y*0.5d0/xtmp
c
      if (x >= 0.d0) dcsqrt = dcmplx (xtmp, ytmp)
      if (y == 0.d0) y = 1.d0
      if (x < 0.d0) dcsqrt = dcmplx (dabs(ytmp), dsign(xtmp,y))
c
      return
      end
