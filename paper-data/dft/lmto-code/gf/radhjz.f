      subroutine radhjz(e,r,lmax,ak,aj,dk,dj,job)
C- Spherical Hankels and Bessels, and their radial and energy derivatives
C ----------------------------------------------------------------
Ci Inputs
Ci    e : complex energy
Ci        NB: use e*avw**2 for Andersen convention; see Remarks
Ci    r : radius
Ci        NB: use r/avw for Andersen convention; see Remarks
Ci lmax : make functions from l=0..lmax
Ci  job :  1s digit
Ci         0, makes values and slopes
Ci         1, makes energy derivatives
Ci        10s digit:
Ci         0  use Gunnarsson's conventions; see Remarks
Ci         1  use Andersen's conventions; see Remarks
Ci       100s digit:
Ci         1  Use Neumann functions for the ak
Co Outputs
Co   If 1s digit of job is 0:
Co   ak   :value of Hankel function (Neuman function with 100s digit iopt)
Co   aj   :value of Bessel function
Co   dk   :radial derivative of ak
Co   dj   :radial derivative of aj
Co   Otherwise, ak,aj,dk,dj are the energy derivatives of the above.
Cr Remarks
Cr  Gunnarsson's conventions (PRB 27, 7144, 1983).
Cr     h_0 = -exp(ikr)/ikr, j_0 = sin(kr)/kr, k = sqrt(e), Im k >= 0.
Cr     Functions of odd order are imaginary for k->0.
Cr     Program returns ak = -i k^(l+1) h_l and aj = k^-l j_l.
Cr     These follow Methfessel's conventions for e->0, viz:
Cr    -i h_l -> (2l-1)!!/(kr)^(l+1)  j_l -> (kr)^l/(2l+1)!!
Cr       a_k -> (2l-1)!!/r^(l+1)     a_j -> r^l/(2l+1)!!
Cr  Andersen's conventions (See Gunnarsson paper, above)
Cr     K_l = -1/(2l-1)!! (k avw)^(l+1) i h_l(kr)
Cr     J_l = 1/2*(2l-1)!! (k avw)^(-l) j_l(kr)
Cr     Here avw is some arbitrary length scale, e.g. average WSR.
Cr     By passing scaled energies and radii, viz e*avw**2 and r/avw,
Cr     program returns ak =  K_l and aj = J_l.
Cr     Note that radial derivative is wrt r/avw.
Cr
Cr  The energy derivative of ak(l=0) diverges as 1/sqrt(e) for e->0
Cr  when the ak are Hankel functions.  Because Neumann functions
Cr  are well behaved, they are used for e=0, regardless of job.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer job,lmax
      double precision r
      double complex e,ak(0:lmax),aj(0:lmax),dk(0:lmax),dj(0:lmax)
C Local parameters
      logical loka
      integer l,lp1
      double precision rl,r2,fac2l(0:20)
      double complex er2,phi(-1:20),psi(-1:20),php(0:20),psp(0:20)

      loka = mod(job,100) >= 10
      r2 = r*r
      er2 = e*r2
      if (mod(job,10) == 0) then
C       Returns phi = j_l(kr)/(kr)**l, psi = i h_l(kr)*(kr)**(l+1)
C       with special scaling for Andersen conventions.
        call besslz(er2,job/10,0,lmax+1,phi(0),psi(0))

        rl = 1/r
        do  10  l = 0, lmax
          lp1 = l+1
          rl = rl*r
          ak(l) = psi(l)/(rl*r)
          aj(l) = phi(l)*rl
          if (loka) then
            dk(l) = (l*psi(l) - psi(l+1)*(l+lp1))/(rl*r2)
            dj(l) = (l*phi(l) - er2*phi(l+1)/(l+lp1))*rl/r
          else
            dk(l) = (l*psi(l) - psi(l+1))/(rl*r2)
            dj(l) = (l*phi(l) - er2*phi(l+1))*rl/r
          endif
   10   continue
      else
        call besslz(er2,(job/100)*10+0,-1,lmax+2,phi,psi)
        fac2l(0) = 1
        do  11  l = 0, lmax+1
          fac2l(l+1) = fac2l(l) * (l+l+1)
          php(l) = -0.5d0*r2*phi(l+1)
          psp(l) = psi(l-1)*r2/2
C         psp(l) = +0.5d0*((2*l+1)*psi(l) - psi(l+1))/e
   11   continue
        rl = 1/r
        do  12  l = 0, lmax
          lp1 = l+1
          rl = rl*r
          ak(l) = psp(l)/(rl*r)
          aj(l) = php(l)*rl
          dk(l) = (l*psp(l) - psp(l+1))/(rl*r*r)
          dj(l) = (l*php(l) - er2*php(l+1) - r*r*phi(l+1))*rl/r
          if (loka) then
            ak(l) = ak(l)/fac2l(l)
            aj(l) = aj(l)*fac2l(l)/2
            dk(l) = dk(l)/fac2l(l)
            dj(l) = dj(l)*fac2l(l)/2
          endif
   12  continue
      endif
      end
CC Tests rahjz
C      subroutine fmain
C      implicit none
C      double precision dr,e1(2),e2(2)
C      integer lmxa,i
CC heap:
C      integer w(100000)
C      common /w/ w
C
C      call finits(2,0,0,i)
C      dr = .7d0
C      e1(1) = -.45d0
C      e1(2) = .2d0
C      lmxa = 4
C   99 print *, 'lmax,z,r='
C      read(*,*) lmxa,e1,dr
C
C      call chkrhj(lmxa,e1,dr)
C      end
C      subroutine chkrhj(lmax,e1,r)
CC - Numerically check derivatives generated by radhjz
C      implicit none
C      integer lmax
C      double complex e1
C      double precision r
C
C      logical loka,cmdopt
C      integer l,ioka
C      character*72 strn
C      double precision fac2l(0:20),dx
C      double complex dr,
C     .  ak1(0:20),aj1(0:20),dk1(0:20),dj1(0:20),
C     .  ak2(0:20),aj2(0:20),dk2(0:20),dj2(0:20)
C
C      fac2l(0) = 1.d0
C      do  10  l = 1, 20
C   10 fac2l(l) = fac2l(l-1) * (l+l-1)
C      loka = .false.
C      if (cmdopt('-oka',4,0,strn)) loka = .true.
C      ioka = 0
C      if (loka) ioka = 10
C      if (cmdopt('-neu',4,0,strn)) then
C        ioka = ioka + 100
C      endif
C
C      write(*,345) lmax,e1,r,loka
C  345 format('  lmax=',i1,'  e=',2f8.5,'  r=',f8.5,'  loka=',L1)
C
C      print 346
C  346 format('   l',14x,'ak',12x,'           ',12x,
C     .          '          aj',12x,'           ')
C      call radhjz(e1,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      do  52  l = 0, lmax
C        print 335, l, ak1(l), aj1(l)
C  335   format(i4,6(2f12.6,24x))
C   52 continue
C
C      print *, '  ... Numerically check radial derivatives ...'
C      dr = (3d-5,5d-5)
C      dx = abs(dr)
C      call radhjz(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,ioka)
C      call radhjz(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,ioka)
C      do  20  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dx
C   20 dk2(l) = (ak2(l)-ak1(l))/dx
C      call radhjz(e1,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      print 341
C  341 format('   l',14x,'dk',12x,'        err',12x,
C     .          '          dj',12x,'        err')
C      do  22  l = 0, lmax
C        print 333, l, dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C  333   format(i4,12f12.6)
C   22 continue
C
C      print *, '  ... Numerically check energy derivatives ...'
C      call radhjz(e1-dr/2,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      call radhjz(e1+dr/2,r,lmax,ak2,aj2,dk1,dj1,ioka)
C      do  30  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dr
C   30 dk2(l) = (ak2(l)-ak1(l))/dr
C      call radhjz(e1,r,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      print 341
C      do  32  l = 0, lmax
C        print 333, l, ak1(l), ak1(l)-dk2(l), aj1(l), aj1(l)-dj2(l)
Cc        print 333, l, ak1(l), dk2(l), aj1(l), dj2(l)
C   32 continue
C
C      print *, '  ... Numerically check radial+energy derivatives ...'
C      call radhjz(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      call radhjz(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,ioka+1)
C      do  40  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dx
C   40 dk2(l) = (ak2(l)-ak1(l))/dx
C      call radhjz(e1,r,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      print 341
C      do  42  l = 0, lmax
C        print 333, l, dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C   42 continue
C
C      stop
C
C      end