      subroutine hanszd(mode,r,e,rsml,lmax,hs,dhs,ddhs,hsp,dhsp,ddhsp)
C- Value and some derivatives of smoothed radial Hankel (e < 0)
C  or Neumann functions (e > 0).
C ---------------------------------------------------------------
Ci Inputs
Ci   mode  : tells hansmd what derivatives to make.
Ci         : 1s digit concerns 2nd radial derivative
Ci         : 0 make neither 1st or 2nd radial derivative.
Ci         : >0 make 1st and second radial derivative:
Ci         : 1 ddhs = radial part of Laplacian, 1/r d^2 (r*h_l) / dr^2
Ci         : 2 ddhs = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l
Ci         :   NB This is laplacian of 3-dimensional hs_l YL since
Ci         :      lap(hs_l YL) = lap(hs_l) YL + hs_l lap(YL)
Ci         :                   = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l
Ci         : 3 ddhs = d^2 (h_l) / dr^2
Ci         : 10s digit concerns energy derivative
Ci         : 0 make none of hsp,dhsp,ddhsp
Ci         : >0 make hsp
Ci         : >0 AND 1s digit >0 make all of hsp,dhsp,ddhsp
Ci   r     : radius
Ci   e     : hankel/neumann energy
Ci   rsml  : vector of l-dependent smoothing radii of Hankels
Ci         : EITHER must be specified for lmin..lmax
Ci         : OR     rsml(0) < 0. Implies rsml(l) = -rsml(0) for all l
Ci   lmax  : make function values for l between (0:lmax)
Co  Outputs:
Co   hs    : function values of the radial sm. hankel (e < 0)
Co         : or neumann (e > 0) function h(e,r,0:lmax)
Co         : A solid hankel/neumann is H=h*Y_L, where Y_L are the
Co         : spherical harmonics for unit radius (no scaling by r**l)
Co   dhs   : (mode>0) radial derivative of hs
Co   ddhs  : (mode>0) radial part of Laplacian of hs, i.e. 1/r d^2 (r h) /dr^2
Co         : OR some other second derivative, depending on 1s digit mode
Co   hsp   : (10s digit mode>0) energy derivative of hs
Co   dhsp  : (10s digit mode>0) mixed energy + radial derivative of hs
Co   ddhsp : (10s digit mode>0) radial part of nabla^2 hsp
Cr Remarks
Cr  Note connection with hansmr and hansmz, which make xi(l) = h(l)/r^l
Cr
Cr  See J. Math. Phys. 39, 3393 (1998).
Cr    For radial derivative, see JMP 39, 3393, Eq. 4.7
Cr      h'  = l/r h_l - h_l+1
Cr    Second radial derivative:
Cr      h'' = l(l-1)/r^2 xi_l - (2l+1)/r xi_l+1 + xi_l+2
Cr      1/r d^2/dr^2 (r*h_l) = l(l+1)/r^2 h_l - (2l+3)/r h_l+1 + h_l+2
Cr    Energy derivative:  see JMP 39, 3393, Eq. 7.5
Cr      hp_l = r/2 h_l-1  Special case l=0: hp_0 = 1/2 h_-1 ?
Cr    Mixed energy + radial derivative:
Cr      hp'_l = h(l-1)*l/2 - h(l)*r/2
Cr    Laplacian hp'' = lap(hp_l YL) = lap(hp_l) YL + hp_l lap(YL)
Cr      hp'' = -(2l+3)/2 h_l + h_l+1*r/2
Cu Updates
Cu   03 May 11 Handles mode 11
Cu   29 Feb 08 (S. Lozovoi) adapted from hansmd.f to handle e > 0
C ---------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmax
      double precision r,e,rsml(0:lmax)
      double precision hs(0:lmax),dhs(0:lmax),ddhs(0:lmax)
      double precision hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax)
C ... Local parameters
      logical lrsm0
      integer l,mode0,mode1,n0
      double precision tol
      parameter (n0=12, tol=1d-15)
      double precision xi(-1:n0),fi(-1:n0),rsm,rsx

C     Negative 1st smoothing radius => constant rsm
      lrsm0 = rsml(0) <= 0
C     if (rsml < 0) call rx('hanszd: negative rsml')
      if (lmax+2 > n0) call rx('hanszd: increase parameter n0')
      if (lmax < 0) return

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

      rsm = dabs(rsml(0))
      rsx = -1d3
      do  l = lmax, 0, -1
        if (.not. lrsm0) rsm = rsml(l)
        if (dabs(rsm-rsx) > tol) then
          call hansrz(-rsm,-1,l+2,e,r*r,1,1,11,xi,fi) ! True xi,fi (no scaling)
          rsx = rsm
        endif

        hs(l) = xi(l)
        if (mode0 /= 0) then
          dhs(l)  = xi(l)*l/r - xi(l+1)
          if (mode0 == 1)
     .      ddhs(l) = xi(l)*l*(l+1)/r**2 - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 == 2)
     .      ddhs(l) =                    - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 == 3)
     .      ddhs(l) = xi(l)*l*(l-1)/r**2 - (2*l+1)/r*xi(l+1) + xi(l+2)
        endif
        if (mode1 /= 0) then
          if (e == 0) call rx('hanszd: attempt to take energy derivative with e=0')
          hsp(l)   = xi(l-1)*r/2
          if (mode0 /= 0) then
            dhsp(l)  = (xi(l-1)*l - xi(l)*r)/2
            ddhsp(l) = - (2*l+3)*xi(l)/2 + xi(l+1)*r/2
            if (mode0 == 1) then
              ddhsp(l) = ddhsp(l) + l*(l+1)/r**2*hsp(l)
            elseif (mode0 == 3) then
              call rx(
     .        'hanszd: for energy derivatives mode0=3 not implemented')
            endif
          endif
        endif
      enddo

c The lines below are not needed since xi(-1) is now scaled inside hansrz
c     if (mode1 /= 0) then
c       hsp(0) = xi(-1)/2
c     endif

      end
C      subroutine fmain
C      implicit none
C      integer n0,kmax,lmax,l
C      parameter (n0=4,kmax=3,lmax=4)
C      double precision xi(0:n0),phi(0:n0),hsml0(0:lmax),rsml(0:lmax)
C      double precision hsm(0:lmax),dhsm(0:lmax),ddhsm(0:lmax)
C      double precision hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax),g0l(0:lmax)
C      double precision r,rsm,pi,a,g0,laphs0(0:lmax),eh,xx
C
C      xx = -1
C
C      r = 1.1d0; rsm = 1.198d0
CC     r = 1; rsm=1
C      r=2.3d0 ; rsm=2*r/3; eh=-.5d0
CC      r=2.3d0 ; rsm=.5d0; eh=-.5d0*1d-10
CC      r=2.3d0 ; rsm=2*r/3; eh=-.5d0*1d-10
CC      r=1.92d0 ; rsm=2*r/3; eh=-1d-10
C
CC      r = 2*r
C      do  l = 0, lmax
C        rsml(l) = rsm*(1-l/20d0)
C      enddo
C      pi = 4d0*datan(1d0)
C      a = 1d0/rsm
C      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0
C      print *, 'r=',sngl(r),'   rsm=',sngl(rsm),'   eh=',sngl(eh)
C
CC      print *, '!! scale rsm by .85 to check l-dep matches f orbital'
CC      rsm = rsm * .85d0
C
C      print *, 'compare hanr, hansmz, hanszd, true radial parts of H'
CC     Use as reference standard sm hankel * 1/r^l
CC     hankel * 1/r^l
C      call hanr(r**2,0,lmax,1,1,eh,xi)
C      print 333, 'hanr','h ',(xi(l)*r**(l), l=0,lmax)
C      call hansmz(r,eh,-rsm,xi,phi,lmax)
C      print 333, 'hansmz','hs',(xi(l)*r**(l), l=0,lmax)
C      call hanszd(0,r,eh,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      print 333, 'hanszd','hs',hsm
C  333 format(1x,a6,2x,a3,5f16.10)
C  334 format(1x,a11,5f16.10)
C      do  l = lmax, 0, -1
C      call hanszd(0,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,xx,xx,xx)
C      enddo
CC     print 333, 'l-dep','h',hsm
C      call hanszd(1,r,eh,rsml,lmax,hsm,dhsm,ddhsm,xx,xx,xx)
C      print 333, 'l-dep','h',hsm
C
C      print *
C      print *, 'dh=dhsm/dr, compared to num diff'
C      call hanszd(0,r+1d-4,eh,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      hsml0 = hsm(:)
C      call hanszd(0,r-1d-4,eh,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      hsml0 = (hsml0 - hsm(:))/2d-4
C      print 333, 'hanszd','num',hsml0
C      call hanszd(1,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,xx,xx,xx)
C      print 333, 'hanszd','dh',dhsm
C      do  l = lmax, 0, -1
C      call hanszd(1,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,xx,xx,xx)
C      enddo
CC     print 333, 'l-dep','dh',dhsm
C      call hanszd(1,r,eh,rsml,lmax,hsm,dhsm,ddhsm,xx,xx,xx)
C      print 333, 'l-dep','dh',dhsm
C
C      print *
C      print *, 'hp = dhsm/de, compared to num diff'
C      call hanszd(0,r,eh+1d-4,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      hsml0 = hsm(:)
C      call hanszd(0,r,eh-1d-4,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      hsml0 = (hsml0 - hsm(:))/2d-4
C      print 333, 'hanszd','num',hsml0
C      call hanszd(10,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'hanszd','hp',hsp
C      do  l = lmax, 0, -1
C      call hanszd(10,r,eh,-rsm,l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      enddo
CC     print 333, 'l-dep','hp',hsp
C      call hanszd(10,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','hp',hsp
C
CC     print *, 'scale rsm'; rsm = rsm * 0.80d0
C      print *
C      print *, 'dh=dhp/dr, compared to num diff'
C      call hanszd(10,r+1d-4,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      hsml0 = hsp(:)
C      call hanszd(10,r-1d-4,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      hsml0 = (hsml0 - hsp(:))/2d-4
C      print 333, 'hanszd','num',hsml0
C      call hanszd(11,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'hanszd','dh',dhsp
C      do  l = lmax, 0, -1
C      call hanszd(11,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      enddo
C      print 333, 'l-dep','dh',dhsp
C      call hanszd(11,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','dh',dhsp
C
C      print *
C      print *, 'nabla^2 hs = 1/r d^2(r*hs)/dr^2,',
C     .  ' compared to num diff'
C      call hanszd(0,r+1d-4,eh,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      laphs0 = (r+1d-4)*hsm(:)
C      call hanszd(0,r-1d-4,eh,-rsm,lmax,hsm,dhsm,xx,xx,xx,xx)
C      laphs0 = laphs0 + (r-1d-4)*hsm(:)
C      call hanszd(1,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,xx,xx,xx)
C      laphs0 = (laphs0 - 2*(r)*hsm(:))/(1d-4)**2/r
C      print 333, 'hanszd','num',laphs0
C      print 333, 'hanszd','ddh',ddhsm
C
C      print *
C      print *, '3D laplacian 1/r d^2(r*hs)/dr^2 - l(l+1)/r^2 hs,',' compared to num diff'
C      call hanszd(2,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,xx,xx,xx)
C      print 333, 'mode 2','num',(laphs0(l)-l*(l+1)/r**2*hsm(l),l=0,lmax)
C      print 333, 'mode 2','ddh',ddhsm
C      print *, '... the next two lines compare (nabla+e)h to -4*pi*g and should match'
C      call radgkl(r,rsm,0,lmax,0,g0l)
C      forall (l = 0:lmax) g0l(l) = g0l(l) * r**l * dexp(eh*rsm**2/4)
C      print 334, 'ddh+e*hs',ddhsm + eh*hsm
C      print 334, '-4pi*g0',-4*pi*g0l
C
C
CC     do  l = lmax, 0, -1
CC     call hanszd(2,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
CC     enddo
CC     print 333, 'l-dep','ddh',ddhsm
C      call hanszd(2,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','ddh',ddhsm
C
C
C      print *
C      print *, 'nabla^2 hp = 1/r d^2(r*hp)/dr^2, compared to num diff'
C      call hanszd(11,r+1d-4,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      laphs0 = (r+1d-4)*hsp(:)
C      call hanszd(11,r-1d-4,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      laphs0 = laphs0 + (r-1d-4)*hsp(:)
C      call hanszd(11,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      laphs0 = (laphs0 - 2*(r)*hsp(:))/(1d-4)**2/r
C      print 333, 'hanszd','num',laphs0
C      print 333, 'hanszd','ddh',ddhsp
C      call hanszd(12,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'mode 2','num',(laphs0(l)-l*(l+1)/r**2*hsp(l),l=0,lmax)
C      print 333, 'mode 2','ddh',ddhsp
C      do  l = lmax, 0, -1
C      call hanszd(12,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      enddo
CC     print 333, 'l-dep','ddh',ddhsp
C      call hanszd(12,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','ddh',ddhsp
C
C      print *
C      print *, 'nabla^2 hs / hs'
C      call hanszd(12,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'mode 2','ddh',(ddhsm(l)/hsm(l),l=0,lmax)
C
C      do  l = lmax, 0, -1
C      call hanszd(12,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      enddo
CC     print 333, 'l-dep','ddh',(ddhsm(l)/hsm(l),l=0,lmax)
C      call hanszd(12,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','ddh',(ddhsm(l)/hsm(l),l=0,lmax)
C
C      call hanszd(12,r,eh,-rsm,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print *, 'nabla^2 hp / hp'
C      print 333, 'mode 2','ddh',(ddhsp(l)/hsp(l),l=0,lmax)
C      do  l = lmax, 0, -1
C      call hanszd(12,r,eh,-rsml(l),l,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      enddo
CC     print 333, 'l-dep','ddh',(ddhsp(l)/hsp(l),l=0,lmax)
C      call hanszd(12,r,eh,rsml,lmax,hsm,dhsm,ddhsm,hsp,dhsp,ddhsp)
C      print 333, 'l-dep','ddh',(ddhsp(l)/hsp(l),l=0,lmax)
C
C
C      end
