      subroutine hsmml(p,rsml,e,lmax,hsm,hsmp)
C- Solid smooth hankel function and energy derivative
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsml  : vector of l-dependent smoothing radii of Hankels
Ci         : EITHER must be specified for lmin..lmax
Ci         : OR     rsml(0) < 0. Implies rsml(l) = -rsml(0) for all l
Ci   e     :energy of smoothed Hankel
Ci   lmax  :l-cutoff for hsm
Co Outputs
Co   hsm   :Bloch-summed smoothed Hankels
Co   hsmp  :Energy derivatives of hsm
Cl Local variables
Cl         :
Cr Remarks
Cr   Correctly handles the p=0 limit, since yl(ilm>1) = 0
Cu Updates
Cu   18 May 11 smoothing radius made l-dependent
Cu   24 Jan 07 New routine
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax
      double precision p(3),rsml(0:lmax),e,yl((lmax+1)**2)
      double precision hsm((lmax+1)**2),hsmp((lmax+1)**2)
C ... Local parameters
      integer ilm,l,m
      double precision r2,r,rl
      double precision wk(0:20),wkd(0:20),xx
C     double precision eps,wd1,wd2,wd
C ... External calls
      external hansmr,ropyln

C  ... Debugging check: xidot by brute force
C      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
C      r = dsqrt(r2)
C      call hansmr(r,e,1/rsml(0),wk(0),lmax)
C      eps = .001d0
C      call hansmr(dsqrt(r2),e+2*eps,1/rsm,wkd(1),0)
C      call hansmr(dsqrt(r2),e+1*eps,1/rsm,wkd(2),0)
C      call hansmr(dsqrt(r2),e-1*eps,1/rsm,wkd(3),0)
C      call hansmr(dsqrt(r2),e-2*eps,1/rsm,wkd(4),0)
C      wd1 = (wkd(1)-wkd(4))/4/eps
C      wd2 = (wkd(2)-wkd(3))/2/eps
C      wd  = (4*wd2-wd1)/3
C      wk(-1) = 2*wd
C      print '(4f20.15)', wk(-1:2)

C ... Make YL*r^l
      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
      r = dsqrt(r2)

C ... Radial part of smooth Hankel, and energy derivative
      call hanszd(10,r,e,rsml,lmax,wk,xx,xx,wkd,xx,xx)

C ... Convert to xi(l) = h(l)/r^l
      if (r /= 0) then
        rl = r
        do  l = 1, lmax
          wk(l) = wk(l)/rl
          wkd(l) = wkd(l)/rl
          rl = rl*r
        enddo
      endif
C     print '(4f20.15)', wk(-1:2)

C ... Solid Hankel is product of radial part and YL*r^l
      ilm = 0
      do  l = 0, lmax
        do  m = -l, l
          ilm  =  ilm+1
          hsm(ilm)  = yl(ilm)*wk(l)
          hsmp(ilm) = yl(ilm)*wkd(l)
        enddo
      enddo
      end

C     test rsm-dependent call
C      subroutine fmain
C      integer lmax
C      parameter (lmax=3)
C      double precision p(3),rsml(0:lmax),rsm1,eh
C      double precision hsm((lmax+1)**2),hsmp((lmax+1)**2)
C      double precision hsm0((lmax+1)**2),hsmp0((lmax+1)**2)
C      integer l,ilm
C
C      data p /-2.4854171306613027d0,
C     .        -2.2983342613226059d0,
C     .        -2.1112513919839087d0/
C
C      rsm1 = 2.3d0 * 2.0d0/3
C      do  l = 0, lmax
C        rsml(l) = rsm1*(1-1d0*l/20d0)
C      enddo
C      eh = -0.5d0
C
C      do  l = lmax, 0, -1
C        call hsmml(p,-rsml(l),eh,l,hsm0,hsmp0)
C      enddo
C      call hsmml(p,rsml,eh,lmax,hsm,hsmp)
C
C      do  ilm = 1, (lmax+1)**2
C        print 333, ilm, hsm(ilm), hsm(ilm)-hsm0(ilm),
C     .    hsmp(ilm), hsmp(ilm)-hsmp0(ilm)
C  333   format(i4,4f15.10)
C      enddo
C
C      end
