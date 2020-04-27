      subroutine z0head(nrad,rad,nlml,lmax,rmt,rsm,xi,vhal,zet0s,sg0)
c  On-site head terms of zet0s (integrals of gaussian times i-pot).
c  ves must be: potential on mesh with weights multiplied in.
      implicit real*8 (a-h,p-z), integer (o)
      dimension vhal(nrad,nlml),zet0s(1),xi(nrad,0:1),rad(1),
     .   sg0(1),xi0(0:3)
      pi = 4d0*datan(1d0)
      lmaxl = ll(nlml)
      lx = min0(lmax,lmaxl)
      nlmx = (lx+1)**2
C --- Tabulate gaussians, l=0 explicitly, the rest by recursion ---
      a2 = 1/rsm**2
      a22 = 2*a2
      do  10  ir = 1, nrad
   10 xi(ir,0) = (a2/pi)**1.5d0*dexp(-a2*rad(ir)**2)
      do  12  l = 1, lx
      do  12  ir = 1, nrad
   12 xi(ir,l) = xi(ir,l-1)*a22*rad(ir)
C --- Integral of gaussians * ves ---
      ilm = 0
      do  20  l = 0, lx
      do  20  m = -l, l
        ilm = ilm+1
        do  21  ir = 1, nrad
   21   zet0s(ilm) = zet0s(ilm) + xi(ir,l)*vhal(ir,ilm)
   20 continue

C --- Add to integral of gaussians over interstitial, sg0 ----
      call hansmr(rmt,0d0,1d0/rsm,xi0,1)
      s00=(1d0-rmt**3*xi0(1))/dsqrt(4d0*pi)
      sg0(1)=sg0(1)+s00

      end
