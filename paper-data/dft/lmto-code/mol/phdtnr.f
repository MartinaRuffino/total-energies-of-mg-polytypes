      subroutine phdtnr(z,l,v,e,a,b,rofi,nr,g,val,slo,gp,gpp,
     .   phi,dphi,phip,dphip,p)
c  make phidot,phidotdot to given e for non-relativistic case.
c  g,val,slo should be wvfct and bcs for normalized solution on input
      implicit real*8 (a-h,p-z), integer(o)
      dimension v(1),rofi(1),g(1),gp(1),gpp(1)
      rmax=rofi(nr)
      call rsqnri(e,l,gp,g,v,z,nodp,vlp,slp,qp,a,b,rofi,nr)
      call gint(g,gp,a,b,nr,sgp)
      call gint(gp,gp,a,b,nr,p)
      do 58 ir=1,nr
  58  gp(ir)=gp(ir)-sgp*g(ir)
      p=p-sgp*sgp
      call rsqnri(e,l,gpp,gp,v,z,npp,vpp,spp,qp,a,b,rofi,nr)
      call gint(g,gpp,a,b,nr,sum)
      xxx=-p-2.d0*sum
      do 50 ir=1,nr
  50  gpp(ir)=2.d0*gpp(ir)+xxx*g(ir)
      srdrdi=dsqrt(a*(rmax+b))
      phi=val/rmax
      dphi=slo/rmax-val/rmax/rmax
      phip=vlp/rmax-sgp*phi
      dphip=(slp-vlp/rmax)/rmax-sgp*dphi
      return
      end
c#
      subroutine gint(g1,g2,a,b,nr,sum)
c integrate product of two fcts, both given as g=u/sqrt(drdi)
      implicit real*8 (a-h,p-z), integer(o)
      dimension g1(nr),g2(nr)
      ea2=dexp(a+a)
      ea4=ea2*ea2
      sum=0.d0
      drdi2=b*b*a*a*ea2
      do 10 ir=2,nr-1,2
      sum=sum+g1(ir)*g2(ir)*drdi2
  10  drdi2=drdi2*ea4
      sum=sum+sum
      drdi2=b*b*a*a*ea4
      do 11 ir=3,nr-2,2
      sum=sum+g1(ir)*g2(ir)*drdi2
  11  drdi2=drdi2*ea4
      rmpb=b*dexp(a*nr-a)
      sum=(sum+sum+g1(1)*g2(1)*(a*b)**2+g1(nr)*g2(nr)*(a*rmpb)**2)/3.d0
      return
      end
c#  rsqnri
      subroutine rsqnri(e,l,g,g0,v,z,nodes,val,slo,q,a,b,rofi,nr)
c  integrates outward, making a special solution of
c  the inhomogeneous schroedinger equation (hl-e)g=g0.
c  g,b are both defined as:  r*(real radial function)/srdrdi.
      implicit real*8 (a-h,p-z), integer(o)
      dimension g(nr),g0(nr),v(nr),rofi(nr),cc(20)
      ncc=15
      a2b4=a*a/4.d0
c ---- power series near zero ---
      lp1=l+1
      r2=rofi(2)
      r3=rofi(3)
      x2=g0(2)*dsqrt(a*(r2+b))/r2**lp1
      x3=g0(3)*dsqrt(a*(r3+b))/r3**lp1
      blp1=-(r2*x3-r3*x2)/(r2-r3)
      blp20=-(x3-x2)/(r3-r2)
       blp2=-z*blp1/lp1
      fllp1=l*lp1
      zz=2.d0*z
      do 10 icc=1,lp1
  10  cc(icc)=0.d0
      vme=v(1)-e
      cc(lp1)=1.d0
      cc(l+2)=-z*cc(lp1)/lp1
      cc(l+3)=(-zz*cc(l+2)+vme*cc(l+1)+blp1)/(4*l+6)
      cc(l+4)=(-zz*cc(l+3)+vme*cc(l+2)+blp2)/(6*l+12)
      do 11 icc=l+5,ncc
  11  cc(icc)=(-zz*cc(icc-1)+vme*cc(icc-2))/(icc*(icc-1)-fllp1)
      g(1)=0.d0
      q=0.d0
      nodes=0
      slo=dsqrt(a*b)*2*a*b*cc(2)
      do 12 ir=2,3
        r=rofi(ir)
        sum=0.d0
        do 13 icc=1,ncc
  13      sum=r*sum+cc(ncc-icc+1)
 12     g(ir)=r*sum/dsqrt(a*(r+b))
c --------------------------
      xb12=1.d0/12.d0
      x=0.d0
      do 16 ir=2,3
      r=rofi(ir)
      s=(a*(r+b))**2
      wgt=2*(mod(ir+1,2)+1)
      t=((fllp1/r-zz)/r+v(ir)-e)*s+a2b4
      slo=slo+wgt*(g(ir)*t-g0(ir)*s)
      q=q+wgt*g(ir)*g(ir)*s
      y=g(ir)*(1.d0-t/12.d0)
      if(ir == 2) b2=-g0(2)*s
      if(ir == 3) b3=-g0(3)*s
  16  x=y-x
      ir=3
  17  ir=ir+1
      r=rofi(ir)
      s=(a*(r+b))**2
      b4=-g0(ir)*s
      wgt=2*(mod(ir+1,2)+1)
      x=x+t*g(ir-1)+(b4+10.d0*b3+b2)*xb12
      y=y+x
      t=((fllp1/r-zz)/r+v(ir)-e)*s+a2b4
      g(ir)=y/(1.d0-t*xb12)
      if(g(ir)*g(ir-1) < 0.d0) nodes=nodes+1
      slo=slo+wgt*(g(ir)*t+b4)
      q=q+wgt*g(ir)*g(ir)*s
      b2=b3
      b3=b4
      if(ir < nr) goto 17
      q=(q-g(nr)*g(nr)*s)/3.d0
      slo=(slo-(g(nr)*t+b4))/3.d0+dsqrt(a*b)*cc(1)
      slo=(.5d0*a*g(nr)+slo)/dsqrt(a*(r+b))
      val=g(nr)*dsqrt(a*(r+b))
      return
      end
