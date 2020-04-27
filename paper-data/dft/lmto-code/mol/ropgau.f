      subroutine ropgau(rsm,lmax,n,r,xi,job)
C- Vector xi of energy-independent gaussians.
c  Set job=1 to multiply xi_l by r**l (true radial part)
      implicit none
      integer lmax,n,job,ir,l
      double precision rsm,r(n),xi(n,0:lmax)
      double precision c,r2,a22,pi,a2

      pi = 4*datan(1d0)
      a2 = 1/(rsm**2)
      c = (a2/pi)**1.5d0
      if(lmax < 0) return
C --- Do l=0 explicitly ---
      do  10  ir = 1, n
   10 xi(ir,0) = c*dexp(-a2*r(ir)**2)
      if(lmax <= 0) return
C --- Radial part/r**l if job == 0, else true radial part ---
      a22 = 2*a2
      if (job == 0) then
        do  20  l = 1, lmax
        do  20  ir = 1, n
   20   xi(ir,l) = xi(ir,l-1)*a22
      else
        do  30  l = 1, lmax
        do  30  ir = 1, n
   30   xi(ir,l) = xi(ir,l-1)*(a22*r(ir))
      endif

      end
