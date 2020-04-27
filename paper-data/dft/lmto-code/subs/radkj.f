      subroutine radkj(e,r,lmax,ak,aj,dk,dj,job)
C- Radial parts of spherical hankels and bessels, with radial and energy derivatives
C ----------------------------------------------------------------
Ci Inputs
Ci   e     :energy (can be either sign)
Ci   r     :radius
Ci   lmax  :maximum l
Ci   job   :1s digit
Ci         : 0, makes values and slopes; 1, makes energy derivatives.
Ci         :10s digit: use Andersen's conventions
Ci         :0 Methfessel's conventions
Ci         :  Radial derivative dk(l) = l*ak(l)/r - ak(l+1)
Ci         :  Radial derivative dj(l) = l*aj(l)/r - e*aj(l+1)
Ci         :1 Andersen conventions from 2nd generation LMTO
Ci         :  Radial derivative dk(l) = l*ak(l)/r - ak(l+1)*(2*l+1)
Ci         :  Radial derivative dj(l) = l*aj(l)/r - e*aj(l+1)/(2*l+1)
Ci         :2 the following convention:
Ci         :  aj = (2l+1)!! * j_l(i sqrt(abs(E)) r) / (i sqrt(abs(E)))**l
Ci         :  ak = (2l-1)!! * h_l(i sqrt(abs(E)) r) * i*(i sqrt(abs(E)))**(l+1)
Ci         :  This definition has the property
Ci         :  ak -> r*(-l-1) as E->0 (OKA definition; see besslr)
Ci         :  aj -> r**l     as E->0 (OKA definition * (2(2l+1))
Ci         :  Radial derivative dk(l) = l*ak(l)/r - ak(l+1)*(2*l+1)
Ci         :  Radial derivative dj(l) = l*aj(l)/r - e*aj(l+1)/(2*l+3)
Ci         :100s digit
Ci         : 0 always use besslr
Ci         : 1 call besnu for y>=90
Ci         : 2 always call dbesnu
Co Outputs
Co   If 1s digit of job is 0:
Co   ak   :value of Hankel function (Neuman function for e>0)
Co   aj   :value of Bessel function
Co   dk   :radial derivative of ak
Co   dj   :radial derivative of aj
Co   Otherwise, ak,aj,dk,dj are the energy derivatives of the above.
Cr Remarks
Cr  For radial derivative of h, see JMP 39, 3393, Eq. 4.7
Cr
Cr  The e->0- limit of dk(l=0) diverges as 1/sqrt(e).
Cr  The e->0+ limit is well behaved, where the ak are Neumann functions.
Cu Updates
Cu   03 Jan 18  Added OKA=2 convention
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer job,lmax
      double precision e,r
      double precision ak(0:lmax),aj(0:lmax),dk(0:lmax),dj(0:lmax)
C Local parameters
      integer l,lp1,loka
      double precision er2,rl,r2,fac2l(0:20)
      double precision phi(-1:20),psi(-1:20),php(0:20),psp(0:20)

      loka = mod(job/10,100)

      r2 = r*r
      er2 = e*r2
      if (mod(job,10) == 0) then
        call besslr(er2,loka,0,lmax+1,phi(0),psi(0))

        rl = 1/r
        do  l = 0, lmax
          lp1 = l+1
          rl = rl*r
          ak(l) = psi(l)/(rl*r)
          aj(l) = phi(l)*rl
          if (loka == 1) then
            dk(l) = (l*psi(l) - psi(l+1)*(l+lp1))/(rl*r2)  ! l/r*H_l - H_l+1*(2l+1)
            dj(l) = (l*phi(l) - er2*phi(l+1)/(l+lp1))*rl/r ! l/r*H_l - e*J_l+1/(2l+1)
          elseif (loka == 2) then  ! Same as OKA for dk
            dk(l) = (l*psi(l) - psi(l+1)*(l+lp1))/(rl*r2)  ! l/r*H_l - H_l+1*(2l+1)
            dj(l) = (l*phi(l) - er2*phi(l+1)/(l+lp1+2))*rl/r ! l/r*H_l - e*J_l+1/(2l+3)
          else
            dk(l) = (l*psi(l) - psi(l+1))/(rl*r2)   ! l/r*H_l - H_l+1
            dj(l) = (l*phi(l) - er2*phi(l+1))*rl/r  ! l/r*J_l - e*J_l+1
          endif
        enddo
C       MSM definitions
C        do l = 0, lmax-1
C          print *, l, dk(l), l*ak(l)/r - ak(l+1) - dk(l)
C        enddo
C        print *
C        do l = 0, lmax-1
C          print *, l, dj(l), l*aj(l)/r - e*aj(l+1) - dj(l)
C        enddo

C       OKA=1 definitions
C        print *
C        do l = 0, lmax-1
C          print *, l, dk(l), l*ak(l)/r - ak(l+1)*(2*l+1) - dk(l)
C        enddo
C        print *
C        do l = 0, lmax-1
C          print 333, l, aj(l), dj(l), l*aj(l)/r - e*aj(l+1)/(2*l+1) - dj(l)
C  333     format(i4,5f15.10)
C        enddo

C       OKA=2 definitions
C        print *
C        do l = 0, lmax-1
C          print *, l, dk(l), l*ak(l)/r - ak(l+1)*(2*l+1) - dk(l)
C        enddo
C        print *
C        do l = 0, lmax-1
C          print 333, l, aj(l), dj(l), l*aj(l)/r - e*aj(l+1)/(2*l+3) - dj(l)
C  333     format(i4,5f15.10)
C        enddo

      else
        if (loka == 1) call rx('radkj OKA=1 not properly implemented for energy derivatives')
        if (loka == 2) call rx('radkj OKA=2 not implemented for energy derivatives')
        call besslr(er2,loka,-1,lmax+2,phi,psi)
        fac2l(0) = 1
        do  l = 0, lmax+1
          fac2l(l+1) = fac2l(l) * (l+l+1)
          php(l) = -0.5d0*r2*phi(l+1)
          psp(l) = psi(l-1)*r2/2
        enddo
        rl = 1/r
        do  l = 0, lmax
          lp1 = l+1
          rl = rl*r
          ak(l) = psp(l)/(rl*r)
          aj(l) = php(l)*rl
          dk(l) = (l*psp(l) - psp(l+1))/(rl*r*r)
          dj(l) = (l*php(l) - er2*php(l+1) - r*r*phi(l+1))*rl/r
          if (loka /= 0) then
            ak(l) = ak(l)/fac2l(l)
            aj(l) = aj(l)*fac2l(l)/2
            dk(l) = dk(l)/fac2l(l)
            dj(l) = dj(l)*fac2l(l)/2
          endif
        enddo
      endif
      end
CC Tests radkj
C      subroutine fmain
C      implicit none
C      double precision dr,e1(2),e2(2)
C      integer lmxa,i
C
C      call finits(2,0,0,i)
C      dr = .95d0
C      e1(1) = -.15d0
C      e1(2) = .2d0
C      lmxa = 4
C   99 print *, 'lmax,e1,e2,r='
C      read(*,*) lmxa,e1(1),e2(1),dr
C
C      call chkrkj(lmxa,e1,dr)
C      end
C      subroutine chkrkj(lmax,e1,r)
CC - Numerically check derivatives generated by radkj
C      implicit none
C      integer lmax
C      double precision e1(2), r
C      logical loka,cmdopt
C      integer l,ioka
C      double precision dx
C      character*32 strn
C
C      double precision fac2l(0:20),
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
C
C      write(*,345) lmax,e1(1),r,loka
C  345 format('  lmax=',i1,'  e=',f8.5,'  r=',f8.5,'  loka=',L1)
C
C      print *, '  ... Numerically check radial derivatives ...'
C      dx = 1d-4
C      call radkj(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,ioka)
C      call radkj(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,ioka)
C      do  l = 0, lmax
C        dj2(l) = (aj2(l)-aj1(l))/dx
C        dk2(l) = (ak2(l)-ak1(l))/dx
C      enddo
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      print *, '  l        ak        aj           dk        err          dj        err'
C      do  l = 0, lmax
C        print 333, l, ak1(l), aj1(l), dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C  333   format(i4,6f12.6)
C      enddo
C
C      print *, '  ... Numerically check radial derivatives, OKA=1 ...'
C      dx = 1d-4
C      call radkj(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,10)
C      call radkj(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,10)
C      do  l = 0, lmax
C        dj2(l) = (aj2(l)-aj1(l))/dx
C        dk2(l) = (ak2(l)-ak1(l))/dx
C      enddo
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,10)
C      print *, '  l        ak        aj           dk        err          dj        err'
C      do  l = 0, lmax
C        print 333, l, ak1(l), aj1(l), dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C      enddo
C      print *, '  ... Numerically check radial derivatives, OKA=2 ...'
C      dx = 1d-4
C      call radkj(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,20)
C      call radkj(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,20)
C      do  l = 0, lmax
C        dj2(l) = (aj2(l)-aj1(l))/dx
C        dk2(l) = (ak2(l)-ak1(l))/dx
C      enddo
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,20)
C      print *, '  l        ak        aj           dk        err          dj        err'
C      do  l = 0, lmax
C        print 333, l, ak1(l), aj1(l), dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C      enddo
C
C      print *, '  ... Numerically check energy derivatives ...'
C      dx = 1d-4
C      call radkj(e1(1)-dx/2,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      call radkj(e1(1)+dx/2,r,lmax,ak2,aj2,dk1,dj1,ioka)
C      do  30  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dx
C   30 dk2(l) = (ak2(l)-ak1(l))/dx
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      print *, '  l        dk        err          dj        err'
C      do  32  l = 0, lmax
C        print 333, l, ak1(l), ak1(l)-dk2(l), aj1(l), aj1(l)-dj2(l)
C   32 continue
C
C      print *, '  ... Numerically check radial+energy derivatives ...'
C      dx = 1d-4
C      call radkj(e1,r-dx/2,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      call radkj(e1,r+dx/2,lmax,ak2,aj2,dk1,dj1,ioka+1)
C      do  40  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dx
C   40 dk2(l) = (ak2(l)-ak1(l))/dx
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,ioka+1)
C      print *, '  l        dk        err          dj        err'
C      do  42  l = 0, lmax
C        print 333, l, dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C   42 continue
C
C      end
