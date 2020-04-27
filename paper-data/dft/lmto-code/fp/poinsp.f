      subroutine poinsp(z,vval,nlm,a,b,v,rofi,rho,rho0,nr,rhoves,rhves1,
     .  vnucl,vsum)
C- Solves Poisson equation inside sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   vval  :boundary conditions: potential for channel ilm = vval(ilm)
Ci         :Note : vval(1) includes nuclear contribution -2z/y0/rmax
Ci         :See remarks on v below
Ci   nlm   :L-cutoff for density
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rho   :density, represented as sum_ilm rho(r,ilm) Y_L(ilm)
Ci   rho0  :work array, holding spherical charge density times 4*pi*r*r
Ci   nr    :number of radial mesh points
Co Outputs
Co   rofi  :radial mesh points
Co   v     :v(ir,ilm) * YL(ilm) = electrostatic potential,
Co         :excluding nuclear contribution -2z/y0/rmax
Co         :Note vval(1) does include nuclear contribution
Co         :Thus v(nr,1) = vval(1) + 2*z/y0/rmax
Co   rhoves:electrostatic energy
Co   rhves1:contribution to electrostatic energy from nonspherical rho
Co   vnucl :potential at the nucleus
Co   vsum  :integral over that potential which is zero at rmax.
Cr Remarks
Cr   Poisson's equation is d2u/dr2 = u*l*(l+1)/r**2 - 8pi*rho/r
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   1 May 00 Adapted from nfp poinsp.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlm,nr
      double precision z,a,b,rhoves,rhves1,vnucl,vsum
      double precision rho(nr,nlm),vval(nlm),v(nr,nlm),rofi(nr),
     .  vhrho(2),rho0(nr)
C ... Local parameters
      integer ilm,intopt,ir,l,ll,lp1,nglob,nsp
      double precision a2b4,alfa,atepi,df,drdi,ea,f,fllp1,fpi,g,r,rmax,
     .  rpb,srdrdi,srfpi,sum,vhom,vnow,x,y0,y2,y3,y4,wt(nr)

      nsp = 1
      fpi = 16d0*datan(1d0)
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi
      atepi = 2d0*fpi
      ea = dexp(a)
      a2b4 = a*a/4d0
      rpb = b
      do  ir = 1, nr
        rofi(ir) = rpb-b
        rpb = rpb*ea
      enddo
      rmax = rofi(nr)
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rofi(nr),a,nr,wt)

C --- Call poiss0 for l=0 to get good pot for small r ---
      do  ir = 1, nr
        rho0(ir) = rho(ir,1)*srfpi
      enddo
      call poiss0(z,a,b,rofi,rho0,nr,vval(1)*y0,v,vhrho,vsum,nsp)
      rhoves = vhrho(1)
      vnucl = v(1,1)
      do  ir = 1, nr
        v(ir,1) = v(ir,1)*srfpi
      enddo

C --- For each ilm>1, do ---
      rhves1 = 0d0
      do  ilm = 2, nlm
      l = ll(ilm)
      lp1 = l+1
      fllp1 = l*(l+1d0)

C ... Numerov for inhomogeneous solution
      v(1,ilm) = 0d0
      df = 0d0
      do  ir = 2, 3
        r = rofi(ir)
        drdi = a*(r+b)
        srdrdi = dsqrt(drdi)
        g = (r**lp1)/srdrdi
        v(ir,ilm) = r**l
        x = fllp1*drdi*drdi/(r*r) + a2b4
        f = g*(1d0-x/12d0)
        if (ir == 2) y2 = -atepi*rho(2,ilm)*drdi*srdrdi/r
        if (ir == 3) y3 = -atepi*rho(3,ilm)*drdi*srdrdi/r
        df = f-df
      enddo
      do  ir = 4, nr
        r = rofi(ir)
        drdi = a*(r+b)
        srdrdi = dsqrt(drdi)
        y4 = -atepi*drdi*srdrdi*rho(ir,ilm)/r
        df = df+g*x+(y4+10d0*y3+y2)/12d0
        f = f+df
        x = fllp1*drdi*drdi/(r*r) + a2b4
        g = f/(1d0-x/12d0)
        v(ir,ilm) = g*srdrdi/r
        y2 = y3
        y3 = y4
      enddo

C ... Add homogeneous solution alpha*r^l s.t. v(nr,ilm) = vval(ilm)
      vnow = v(nr,ilm)
      vhom = rmax**l
      alfa = (vval(ilm)-vnow)/vhom
      sum = 0d0
      do  ir = 2, nr
        r = rofi(ir)
        v(ir,ilm) = v(ir,ilm) + alfa*(r**l)
        sum = sum + wt(ir)*rho(ir,ilm)*v(ir,ilm)
      enddo
      rhves1 = rhves1 + sum
      v(1,ilm) = 0d0
      enddo
      rhoves = rhoves + rhves1
      end
