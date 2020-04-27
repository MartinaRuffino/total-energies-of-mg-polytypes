      subroutine hyfzw(alfaw,alfal,modef,nalf,a1,z0,nz,z,wz,f)
C- Functions, points, weights for the z integration of the TCF
C  The z-dependence of a fit to TC product p(z) is carried through
C  the fit functions f = c_i f_i(z), with
C  f_i(z) = exp(i*a1(z0-z)), i=1..nalf (modef=0)
C  Note: for now, no other form for f is permitted.
C  Coefficients c_i are generated in a least-squares sense over the
C  interval (z0,infty) wrt a weighting exp(-alfaw z), viz
C  min intgrl_(z0,infty) exp(-alfaw z) (sum (c_i f_i(z)) - p(z))**2
C  Exponenent alfaw may be of either sign.
C
C  Integrals are done numerically with Laguerre gaussian quadrature,
C  with nz points (except when nz=1, where the point z0 is used).
C  alfal is used to generate the z_i and wz_i for the numerical
C  integration, and should approximate the exponent in the asymtotic
C  dependence of functions integrated, i.e. the exponent to f_i f_j
C  exp(-alfaw*z).  For sufficiently large nz, integration is exact
C  independently of alfal; however a judicious choice of alfal
C  permits precise integrals with a minimal nz.
C      implicit none
      double precision a1,alfaw,alfal,z0,z(nz),f(nz,nalf),wz(nz),x,d,xm
      integer nz,ia,ialf,nalf,iz,ipr,modef

C --- Setup ---
      call getpr(ipr)
      if (ipr >= 20) print 200, nalf,a1,alfaw,alfal,z0,nz
  200 format(/' hyfzw:  nalf=',i1,'  a=',f6.3,'  alfaw=',f6.3,
     .  '  alfal=',f6.3,'  z0=',f6.3,'  nz=',i2)
      if (modef == 0) then
        if (alfal < 1d-6) alfal = 2*a1 + alfaw
      else
        call rx('hyfzw:  bad modef')
      endif
      if (alfal < .1d0) call rx('hyfzw:  alfal too small')

C --- z and wz with Laguerre gaussian quadrature (unless nz=1) ---
      z(1) = z0
      wz(1) = 1
      if (nz /= 1) call mklagw(nz,1,alfal,z0,z,wz,0)

C --- Make f, scale wz by exp(-alfaw z) ---
      if (ipr >= 60) write(6,331)
      do  10  iz = 1, nz
        d = z(iz)
        x = dexp(a1*(z0-d))
        wz(iz) = wz(iz)*dexp(-alfaw*d)
        xm = 1d0
        do  11  ialf = 1, nalf
        xm = xm*x
   11   f(iz,ialf) = xm

        if (ipr >= 60) write(6,330) iz,d,wz(iz),(f(iz,ia),ia=1,nalf)
  330   format(i5,2f11.6,4f12.6:/(27x,4f12.6))
  331   format('   iz      z',10x,'w',8x,'f(1..nalf)')
   10 continue
      end
