      subroutine lindhard(omega,q,kf,eps)
C- Return Lindhard function for given frequency, q, fermi level
C ----------------------------------------------------------------------
Ci Inputs
Ci   omega :frequency w
Ci   q     :wave number in atomic units
Ci   kf    :Fermi wave number in atomic units
Ci         :Note: Ef = kf^2 in Atomic Rydberg units
Co Outputs
Co   eps:  :Dielectic function eps(q,omega;kf)
Cr Remarks
Cr   The Lindhard function is
Cr     eps = 1 + (kf/q)^3/(\pi kf a_0) [..]  where
Cr     [..] = [2q/kf + f(z1) + f(z2)]
Cr     f(z) = (1 - (z^2/4)) \ln[(2+z+i\delta)/(2-z-i\delta)]
Cr     z1   = Q + W/Q
Cr     z2   = Q - W/Q
Cr     Q = q/kf and w = omega/ef where ef=hbar^2 kf^2/2m
Cr   Corresponds to Hedin's 5.33: Use \alpha r_s = 1/k_F (A&M 2.22)
Cr   except that the infinitessimal part is treated slightly differently
Cr   i\delta is be added if z is real:
Cr   Im ln[..] .. => +i\pi if z>0 and -i\pi if z<0
Cr
Cr   CASE omega->0 limit : z1=z2=z=q/kf
Cr   [..] = [2z + 2f(z)]
Cr
Cr   CASE q->0 limit and W>0 and Im W >= 0: use q->0 => z->infinity.  For large z,
Cr     f(z) = 32/15z^3 + 8/3z - z \pm i\pi(1-z^2/4) depending on the branch
Cr   If z=z1=Q+W/Q : choose principal branch assuming Im W is small
Cr     f(z1) = -W/Q + (-1+8/3W)Q - [8(-4+5W)/15W^3]Q^3 + i\pi(1-z1^2/4)
Cr   If z=z2=Q-W/Q : choose principal branch assuming Im W is small
Cr     f(z2) = +W/Q + (-1-8/3W)Q - [8(-4+5W)/15W^3]Q^3 - i\pi(1-z2^2/4)
Cr   Note that 1-z1^2/4 + 1-z2^2/4 = -W
Cr     [..] = 2Q + f(z1) + f(z2) = 2Q - 2Q - 16Q^3/3W^2 - i\pi W
Cr   If moreover W is REAL
Cr     [..] = 2Q + f(z1) + f(z2) = 2Q - 2Q - 16Q^3/3W^2 = -16Q^3/3W^2
Cr     1/(pi kf a_0 Q^3) Re [..]  = -(16/3)/(pi kf a_0)/W^2
Cr     eps(q->0,omega) = 1 - wp^2/omega^2 where
Cr     wp^2 = (16/3)/(pi kf a_0) Ef^2/hbar^2
Cr     Cast in traditional form for wp:
Cr     Use Ef^2/hbar^2 = hbar^2(kf^2/2m)^2, a0=hbar^2/me^2, kf^3 = 3 pi^2 n
Cr     wp^2 = (16/3)/(pi kf a_0) hbar^2 (kf^2/2m)^2
Cr          = (16/3)/(pi hbar^2/me^2) (3 pi^2 n) hbar^2 / 4m^2
Cr          = 4 pi n e^2/ m
Cr
Cr   Atomic Rydberg units:  hbar = 2m = e^2/2 = a0 = 1 so  Ef = kf^2
Cu Updates
Cu   26 Apr 14 First written.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8) :: q,kf
      complex(8) :: omega,eps
C ... Local parameters
      double precision pi,n,wp2,rbrak,qbyk,ef,delta
      complex(8) :: zlind,z,w,zbrak
      zlind(z) = (1-z**2/4) * log((2+z)/(2-z)+(0d0,1d0)*delta)

      pi = 4d0*datan(1d0)

      if (q == 0d0 .and. w == 0) then
        call rx('lindhard : illegal q = w = 0')
      elseif (q == 0d0) then
        n = kf**3/(3*pi**2)
        wp2 = 4*pi*n*4    !e^2/m = 4.  Rydberg units: e^2/m = 2 / (1/2) = 4
        wp2 = (16d0/3)/(pi*kf)*kf**4 ! Rydberg units: ef=kf**2
        eps = 1 - wp2/w**2
      elseif (w == 0d0) then  ! From Ashcroft and Mermin
        qbyk = q / (2*kf)
        rbrak = 2 + (1-qbyk**2)/qbyk*dlog(dabs((1+qbyk)/(1-qbyk)))
        eps = 1 + kf/(pi*q**2)*rbrak
      else
        qbyk = q/kf ; ef = kf**2  ! In Atomic Rydberg units
        W = omega/ef

        z = (qbyk + W/qbyk)
        delta = 0
        if (dimag(z) == 0) delta = dsign(1d-16,dble(z))

        delta = 0


        print *, 'z',z,delta
        print *, 'W',W
        print *, 'arg log',(2+z)/(2-z) + (0d0,1d0)*delta
        print *, log((2+z)/(2-z) + (0d0,1d0)*delta)
        print *, 'in line',(1-z**2/4)*log((2+z)/(2-z))
        print *, 'zlind',zlind(z)
        print *, 32d0/15/z**3 + 8/3/z - z
        print *, -W/qbyk+(-1+8/(3*W))*qbyk-8*(-4+5*W)*qbyk**3/(15*W**2)
        zbrak = zlind(z)

        print *, '---'
        z = (qbyk - W/qbyk)
        if (dimag(z) == 0) delta = dsign(1d-16,dble(z))

        delta = 0

        print *, 'z',z
        print *, 'W',W
        print *, 'arg log',(2+z)/(2-z)
        print *, log((2+z)/(2-z))
        print *, 'in line',(1-z**2/4)*log((2+z)/(2-z))
        print *, 'zlind',zlind(z)
        print *, 32d0/15/z**3 + 8/3/z - z
        print *, W/qbyk+(-1-8/(3*W))*qbyk-8*(-4+5*W)*qbyk**3/(15*W**2)
        zbrak = 2*qbyk + zbrak + zlind(z)

        print *, '---'
C        zbrak = 2*qbyk + zlind(qbyk+W/qbyk) + zlind(qbyk-W/qbyk)
        eps = 1 + 1/(pi*kf*qbyk**3)*zbrak

        print *, 'zbrak',zbrak
        print *, 'analytic',-16*qbyk**3/(3*W**2) - (0d0,1d0)*pi
        print *, eps
        wp2 = (16d0/3)/(pi*kf)*kf**4 ! Rydberg units: ef=kf**2
        print *, 1-wp2/omega**2 - (0d0,1d0)*pi*W/(pi*kf*qbyk**3)
        stop
      endif
      end

      subroutine fmain
      double precision q,kf
      double complex omega,eps

      q = .001d0 ; kf = .1d0
      print *, 'q,kf=?'; read(*,*) q,kf
   10 continue
      omega = (.01d0,0d-10)
      print *, 'omega=?'; read(*,*) omega
      call lindhard(omega,q,kf,eps)
      print 333, omega,q,eps
  333 format(2f16.6,1x,f12.6,1x,2f16.6)
      goto 10
      end
