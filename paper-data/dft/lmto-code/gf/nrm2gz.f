      subroutine nrm2gz(a,b,kap2,g,l,liophi,nr,hcr,slo,val,rmax,rofi,
     .  phi0,phio)
C- Normalizes wave function so that phi^0(a)=1
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   kap2  :kinetic energy e-vmtz (complex)
Ci   l     :angular momentum
Ci  liophi:T tabulate phi0 on the radial mesh from hcr to rmax
Ci   nr    :number of mesh points
Ci   hcr   :hard-core screening radius
Ci   slo   :slope of phi=g/r at rmax; see Remarks.
Ci   val   :value of phi=g/r at rmax; see Remarks.
Ci   rmax  :Wigner-Seitz sphere radius, in atomic units
Cio Inputs/Outputs:
Cio  g     :On input, g is phi times r (complex)
Cio        :On output, g is renormalized
Co Outputs:
Co   phio  :value of back extrapolated wavefunction at hcr
Co   phi0  :phi^0 on the radial mesh from hcr to rmax
Cr  Remarks
Cr    This is a complex analog of norm2g, which see.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer l,nr
      logical liophi
      double precision a,b,hcr,rmax,rofi(nr)
      double complex kap2,g(nr),slo,val,phi0(*),phio
C Local variables:
      integer i,imin,iprint,lmxx
      parameter(lmxx=10)
      double precision sigl,ri
      double complex er2,fi(0:lmxx+1),gi(0:lmxx+1),wn,wj,dphi,dn,dj,
     .                 phi,phioi
C External calls:
      external besslz,dscal

      er2 = kap2*rmax*rmax
C     call bessl2(er2,0,l+1,fi(0),gi(0))
      call besslz(er2,1,0,l+1,fi(0),gi(0))
C     phi,dphi are value, slope at rmax
      phi  = val/rmax
      dphi = rmax*slo/val - 1d0
C     dj,dn are the logarithmic derivatives of bessl and hankels
      dn   = l - gi(l+1)/gi(l)*(l+l+1)
      dj   = l - fi(l+1)/fi(l)/(l+l+1)*er2
C     wj,wn are the amounts of bessl and hankel making up phi0
      wj   = (dphi-dj)*fi(l)*phi
      wn   = (dphi-dn)*gi(l)*phi
      sigl = hcr/rmax
      er2  = kap2*hcr*hcr
C     call bessl2(er2,0,l+1,fi(0),gi(0))
      call besslz(er2,1,0,l+1,fi(0),gi(0))
      phio = 2d0*(wn*fi(l)*sigl**l - wj*gi(l)*sigl**(-l-1))
      if (iprint() >= 120) print 301, phi,phio
      call zscal(2*nr,1d0/phio,g,1)
      if (iprint() >= 120) print 300, kap2,val,slo,dphi,g(10),g(100)

C --- Generate phi0 on a radial mesh between hcr and rmax ---
      if (liophi) then
         phioi = 1d0/phio
         imin = int(1d0+dlog(hcr/b+1d0)/a) - 1
         do  10  i = 1, imin-1
   10    phi0(i) = 0d0
         do  20  i = imin, nr
           ri = rofi(i)
           er2 = kap2*ri*ri
C          call bessl2(er2,0,l+1,fi(0),gi(0))
           call besslz(er2,1,0,l+1,fi(0),gi(0))
           phi0(i) = (2d0*(wn*fi(l)*(ri/rmax)**l -
     .                     wj*gi(l)*(ri/rmax)**(-l-1)))*phioi
   20     continue
      endif

  300 format(' NRM2GZ: kap2, val slo dphi g1 g2',12f10.5)
  301 format(' NRM2GZ: phi, phio',4f10.5)
      end
