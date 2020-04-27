      subroutine wronkj(job,e1,e2,r,lmax,fkk,fkj,fjk,fjj)
C- Wronskians for spherical hankels and bessels on one sphere.
C ----------------------------------------------------------------------
Ci Inputs
Ci   e1    :Wronskian for k or j (e1) with k or j (e2)
Ci         :where k = Hankel for e<0 and Neumann for e>0
Ci   e2    :Wronskian for k or j (e1) with k or j (e2)
Ci   r     :radius
Ci   lmax  :maximum l for a given site
Ci   job   :passed as 100s digit job to radkj
Co Outputs
Co   fkk   :W(k,k)
Co   fkj   :W(k,j)
Co   fjk   :W(j,k)
Co   fjj   :W(j,j)
Cr Remarks
Cr  fxy is continuous as e2->e1.
Cr  fkk,fjj are symmetric in e1,e2.  fkj(e1,e2)=fjk(e2,e1).
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmax,job
      double precision e1,e2,r
      real(8) :: fkk(*),fkj(*),fjk(*),fjj(*),ak1(100),aj1(100),
     .   ak2(100),aj2(100),dk2(100),dj2(100),dk1(100),dj1(100)
C ... Local parameters
      integer l,lp1
      real(8) :: r3,rk,rj,efac
      real(8),parameter :: tol=1d-8

C --- Case e1=e2=0 ---
      if (dabs(e1) <= tol .and. dabs(e2) <= tol) then
        r3 = r*r*r
        rk = -1d0
        rj = 1d0/r
        do  l = 0, lmax
          lp1 = l+1
          rk = rk*(2*l-1)/r
          rj = rj*r/(2*l+1)
          fkk(lp1) = rk*rk*r3/(2*l-1)
          fjj(lp1) = -rj*rj*r3/(2*l+3)
          fkj(lp1) = -0.5d0*rj*rk*r3
          fjk(lp1) = fkj(lp1)
        enddo

C --- Case e1 /= e2 ---
      elseif (dabs(e1-e2) > tol) then
        efac = 1d0/(e2-e1)
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,100*job+0)
        call radkj(e2,r,lmax,ak2,aj2,dk2,dj2,100*job+0)
        do  lp1 = 1, lmax+1
          fkk(lp1) = efac*r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1) = efac*r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1) = efac*r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))-efac
          fjk(lp1) = efac*r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))+efac
        enddo

C --- Case e1 == e2 but not zero ---
      else
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,100*job+0)
        call radkj(e1,r,lmax,ak2,aj2,dk2,dj2,100*job+1)
        do  lp1 = 1, lmax+1
          fkk(lp1) = r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1) = r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1) = r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))
          fjk(lp1) = r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))
        enddo
      endif

      end

      subroutine radkj(e,r,lmax,ak,aj,dk,dj,job)
C- Radial parts of spherical hankels and bessels, with radial and energy derivatives
C ----------------------------------------------------------------
Ci Inputs
Ci   e,r,lmax
Ci   job   :1s digit
Ci         : 0, makes values and slopes; 1, makes energy derivatives.
Ci         :10s digit: use Andersen's conventions
Ci         :0 Methfessel's conventions
Ci         :1 Andersen conventions from 2nd generation LMTO
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
Cr  The e->0- limit of dk(l=0) diverges as 1/sqrt(e).
Cr  The e->0+ limit is well behaved, where the ak are Neumann functions.
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
          if (loka /= 0) then
            dk(l) = (l*psi(l) - psi(l+1)*(l+lp1))/(rl*r2)
            dj(l) = (l*phi(l) - er2*phi(l+1)/(l+lp1))*rl/r
          else
            dk(l) = (l*psi(l) - psi(l+1))/(rl*r2)
            dj(l) = (l*phi(l) - er2*phi(l+1))*rl/r
          endif
        enddo
C        do l = 0, lmax-1
C          print *, l, dk(l), l*ak(l)/r - ak(l+1) - dk(l)
C        enddo
C        stop
      else
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
C      dr = .7d0
C      e1(1) = -.45d0
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
C      do  20  l = 0, lmax
C      dj2(l) = (aj2(l)-aj1(l))/dx
C   20 dk2(l) = (ak2(l)-ak1(l))/dx
C      call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,ioka)
C      print *, '  l        dk        err          dj        err'
C      do  22  l = 0, lmax
C        print 333, l, dk1(l), dk1(l)-dk2(l), dj1(l), dj1(l)-dj2(l)
C  333   format(i4,6f12.6)
C   22 continue
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
C      stop
C
C      end
