      subroutine wrnhjz(e1,e2,r,lmax,avw,wkk,wkj,wjk,wjj,job)
C- Wronskians of hankels and bessels, complex energy
C ----------------------------------------------------------------
Ci Inputs
Ci   e1 :  Energy of first hankel and bessel functions
Ci   e2 :  Energy of second hankel and bessel functions
Ci   avw:  length scale in Andersen's definitions of H,J.
Ci   r  :  true radius, in a.u.
Ci   job:  1s digit
Ci         0  makes Wronskians of Hankels and Bessels.
Ci         1  makes Wronskians / (e1-e2) for making
Ci            interstitial integrals, analogous to radkj.
Ci        10s digit:
Ci         0  use Gunnarsson conventions; see Remarks
Ci         1  use Andersen conventions; see Remarks
Ci       100s digit:
Ci         0  K = Hankel function
Ci         1  K = Neumann function
Co Outputs
Co   If 1s digit of job is 0:
Co   wkk,wkj,wjk,wjj : W{K(e1),K(e2)}, etc, with possible
Co                     normalizations as described in Remarks.
Co   If 1s digit of job is 1:
Co   wkk,wkj,wjk,wjj : W{K(e1),K(e2)}/(e1-e2), etc
Cr Remarks
Cr  Program evaluates wronskians of functions ak and aj, which
Cr  conform to Methfessel's or Andersen's conventions; see radkjz.
Cr
Cr  For Gunnarsson's conventions :
Cr  Program uses ak = -i k^(l+1) h_l and aj = k^-l j_l.
Cr  Thus, the Wronskians conforming to the Gunnarsson convention are
Cr    W(h,h) = -k^(-2l-2) wkk
Cr    W(h,j) = i k^-1 wkj
Cr    W(j,h) = i k^-1 wjk
Cr    W(j,j) = (k1 k2)^l wjj
Cr
Cr  The e->0- limit of dk(l=0) diverges as 1/sqrt(e).
Cr  The e->0+ limit would be well behaved, if the ak were
Cr  Neumann functions instead of Hankel functions.  In that case,
Cr  the singularity is removed, thus following the Neumann convention.
Cb Bugs
Cb   wjj was checked for Andersen definitions, and wkk etc
Cb   all agree for real e with old FP makfkj.  These latter were
Cb   not checked for Andersen definitions.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer job,lmax
      double precision r,avw
      double complex e1,e2
      double complex wkk(0:lmax),wkj(0:lmax),wjk(0:lmax),wjj(0:lmax)
C Local parameters
      integer ioka,l
      double precision r2,esmall
      double complex ak1(0:20),aj1(0:20),dk1(0:20),dj1(0:20),efac
      double complex ak2(0:20),aj2(0:20),dk2(0:20),dj2(0:20),erfac
      parameter (esmall=1d-7)

      ioka = 10*mod(job/10,100)
      call radhjz(e1*avw**2,r/avw,lmax,ak1,aj1,dk1,dj1,ioka)
      if (e1 == e2) ioka = ioka+1
      call radhjz(e2*avw**2,r/avw,lmax,ak2,aj2,dk2,dj2,ioka)

C ... Extra factor avw because radkjz differentiated wrt (r/avw)
      r2 = r**2 / avw
      if (mod(job,10) == 0) then
        do  10  l = 0, lmax
          wkk(l) = r2*(ak1(l)*dk2(l) - dk1(l)*ak2(l))
          wjj(l) = r2*(aj1(l)*dj2(l) - dj1(l)*aj2(l))
          wkj(l) = r2*(ak1(l)*dj2(l) - dk1(l)*aj2(l))
          wjk(l) = r2*(aj1(l)*dk2(l) - dj1(l)*ak2(l))
   10   continue
      else
        if (abs(e1-e2) > esmall) then
          efac = 1/(e2-e1)
          erfac = efac*r2
        else
          efac = 0
          erfac = avw**2 * r2
        endif
        do  20  l = 0, lmax
          wkk(l) = erfac*(ak1(l)*dk2(l) - dk1(l)*ak2(l))
          wjj(l) = erfac*(aj1(l)*dj2(l) - dj1(l)*aj2(l))
          wkj(l) = erfac*(ak1(l)*dj2(l) - dk1(l)*aj2(l)) - efac
          wjk(l) = erfac*(aj1(l)*dk2(l) - dj1(l)*ak2(l)) + efac
   20   continue
      endif
      end
CC Tests wrnhjz
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
C      e2(1) = .45d0
C      e2(2) = .2d0
C      lmxa = 4
C   99 print *, 'lmax,e1,e2,r='
C      read(*,*) lmxa,e1,e2,dr
C
C      call chkwkj(lmxa,e1,e2,dr)
Cc     call chkrhj(lmxa,e1,dr)
C      end
C      subroutine chkwkj(lmax,e1,e2,r)
C      implicit none
C      integer lmax
C      double complex e1,e2
C      double precision r
C      logical loka,lgun,cmdopt,a2bin,lneu
C      integer l,ioka,j,n
C      parameter (n=20)
C      character*32 strn
C      double precision fac2l(0:20),dx,avw,xp(n),wp(n)
C      double complex dr,y,wkk(0:20),wkj(0:20),wjk(0:20),wjj(0:20)
C      double precision fkk(0:20),fkj(0:20),fjk(0:20),fjj(0:20)
C      double complex phi1(0:20),psi1(0:20),phi2(0:20),psi2(0:20),srmy,
C     .  ka,sumkk(0:20),sumkj(0:20),sumjk(0:20),sumjj(0:20),sre1,sre2
C
C      avw = 1
C      fac2l(0) = 1.d0
C      do  10  l = 1, 20
C   10 fac2l(l) = fac2l(l-1) * (l+l-1)
C      loka = .false.
C      lgun = .false.
C      if (cmdopt('-neu',4,0,strn)) lneu = .true.
C      if (cmdopt('-gun',4,0,strn)) lgun = .true.
C      if (cmdopt('-oka',4,0,strn)) loka = .true.
C      if (cmdopt('-oka:',5,0,strn)) then
C        j = 5
C        call rxx(.not.a2bin(strn,avw,4,0,' ',j,len(strn)),
C     .    'failed to parse '//strn)
C      endif
C      ioka = 0
C      if (loka) ioka = 10
C      if (lneu) ioka = ioka+100
C
C      write(*,345) lmax,e1,e2,r,loka,lgun,lneu
C  345 format(' lmax=',i1,'  e1=',2f8.5,'  e2=',2f8.5,'  r=',f8.5,
C     .  '  loka=',L1,'  lgun=',L1,'  neu=',L1)
C      print *, '... if e1 = e2, W{h,j} = 1 for loka=f'
C      if (loka) write(*,346) avw
C  346 format('  avw=',f8.5)
C
C      print 341
C  341 format(3x,'l',12x,'wkk',21x,'wkj',21x,'wjk',21x,'wjj')
C
C      call wrnhjz(e1,e2,r,lmax,avw,wkk,wkj,wjk,wjj,ioka)
C      do  22  l = 0, lmax
C        print 333, l, wkk(l), wkj(l), wjk(l), wjj(l)
C  333   format(i4,12f12.6)
C   22 continue
C
C      sre1 = sqrt(e1)
C      if (dimag(sre1) < 0d0) sre1 = -sre1
C      sre2 = sqrt(e2)
C      if (dimag(sre2) < 0d0) sre2 = -sre2
C
C      print *, ' '
C      print *, ' ... test integral J(E1) J2(E2)'
C      print '(''  l'',10x,''numerical'',17x,''-Wjj'',21x,''diff'')'
C      call gausq(n,0d0,r,xp,wp,2,0)
C      call dpzero(sumkk,2*lmax+2)
C      call dpzero(sumkj,2*lmax+2)
C      call dpzero(sumjk,2*lmax+2)
C      call dpzero(sumjj,2*lmax+2)
C      do  40  j = 1, n
C        y = e1*xp(j)**2
CC       Gunnarssonn convention
C        srmy = sqrt(y)
C        if (dimag(srmy) < 0d0) srmy = -srmy
CC       Methfessel convention
C        if (.not. lgun) srmy = xp(j)
CC       Andersen convention
C        if (loka) srmy = xp(j)/avw
C        call besslz(y,ioka/10,0,lmax,phi1(0),psi1(0))
CC       besslz returns x**l j_l, x**2 = e r*r
C        do  42  l = 0, lmax
C          phi1(l) = phi1(l) * srmy**l
C          psi1(l) = psi1(l) / srmy**(l+1)
CC          if (loka) then
CC            phi1(l) = phi1(l) / (sre1*avw)**l
CC            psi1(l) = psi1(l) * (sre1*avw)**(l+1)
CC          endif
C   42   continue
C        y = e2*xp(j)**2
C        call besslz(y,ioka/10,0,lmax,phi2(0),psi2(0))
CC       Gunnarssonn convention
C        srmy = sqrt(y)
C        if (dimag(srmy) < 0d0) srmy = -srmy
CC       Methfessel convention
C        if (.not. lgun) srmy = xp(j)
CC       Andersen convention
C        if (loka) srmy = xp(j)/avw
C        do  44  l = 0, lmax
C          phi2(l) = phi2(l) * srmy**l
C          psi2(l) = psi2(l) / srmy**(l+1)
CC          if (loka) then
CC            phi2(l) = phi2(l) / (sre2*avw)**l
CC            psi2(l) = psi2(l) * (sre2*avw)**(l+1)
CC          endif
C   44   continue
C        do  45  l = 0, lmax
C          sumjj(l) = sumjj(l) + phi1(l) * phi2(l) * wp(j)
C          sumkk(l) = sumkk(l) + psi1(l) * psi2(l) * wp(j)
C   45   continue
C   40 continue
C
C      call wrnhjz(e1,e2,r,lmax,avw,wkk,wkj,wjk,wjj,ioka+1)
C      do  46  l = 0, lmax
C        if (lgun) then
C          wjj(l) = (sre1*sre2)**l * wjj(l)
C        endif
C        print 333, l,sumjj(l),-wjj(l),sumjj(l)+wjj(l)
C   46 continue
C
C      print *, ' '
C      print *, ' ... Compare wkk etc against fkk from old fp'
C      if (dimag(e1) == 0 .and. dimag(e2) == 0 .and.
C     .    .not. loka .and. .not. lgun) then
C        call makfkj(e1,e2,r,lmax,fkk,fkj,fjk,fjj)
C      print 351
C  351 format(3x,'l',6x,'fkk',9x,'wkk',9x,'fkj',9x,'wkj',9x,'fjk',9x,
C     .  'wjk',9x,'fjj',9x,'wjj')
C
C        do  56  l = 0, lmax
C          print 333, l,fkk(l),dble(wkk(l)),fkj(l),dble(wkj(l)),
C     .      fjk(l),dble(wjk(l)),fjj(l),dble(wjj(l))
C   56   continue
C        print *, '     NB: agreement expected only if:'
C        print *, '     for e1<0, e2<0, no -neu  and',
C     .           ' for e1>0, e2>0, use -neu'
C      else
C        print *, '     skip comparison: ',
C     .    'need e real, and no oka, gun'
C      endif
C      stop
C      end
C
C      subroutine makfkj(e1,e2,r,lmax,fkk,fkj,fjk,fjj)
CC- Modified wronskians for hankels and bessels on one sphere.
CC ----------------------------------------------------------------
CCi Inputs
CCi   e1,e2,r,lmax
CCi   r<0 is treated as r>0 except that f's are true Wronskians
CCi   These do not have the property that fxy continuous as e2->e1.
CCo Outputs
CCo   fkk,fkj,fjk,fjj
CCr Remarks
CCr   fkk =  w(k1,k2) / (e2-e1)         e1 ne e2
CCr          w(k,kdot)                  e1 eq e2
CCr   fkj = (w(k1,j2) - 1)/ (e2-e1),    e1 ne e2
CCr          w(k,jdot)                  e1 eq e2
CCr   fjk = (w(j1,k2) + 1)/ (e2-e1),    e1 ne e2
CCr          w(jdot,k)                  e1 eq e2
CCr   fjj =  w(j1,k2) / (e2-e1),        e1 ne e2
CCr          w(j,jdot)                  e1 eq e2
CCr   fxy is continuous as e2 -> e1.
CCr   fkk,fjj are symmetric in e1,e2.  fkj(e1,e2)=fjk(e2,e1).
CCr   Not implemented for OKA's definitions!
CC ----------------------------------------------------------------
CC     implicit none
CC Passed parameters
C      integer lmax
C      double precision e1,e2,r
C      double precision fkk(0:1),fkj(0:1),fjk(0:1),fjj(0:1)
CC Local parameters
C      integer l,job
C      double precision efac,r3,rj,rk,erfac,esmall,rr
C      double precision ak1(0:12),aj1(0:12),ak2(0:12),aj2(0:12),
C     .                 dk2(0:12),dj2(0:12),dk1(0:12),dj1(0:12)
C      parameter (esmall=1d-7)
C
CC#ifdefC OKA
CC      stop 'MAKFKJ:  not implemented for OKA'
CC#endif
C
C      rr = dabs(r)
C
CC --- Special case e1 = e2 = 0 ---
C      if (dabs(e1) <= esmall .and. dabs(e2) <= esmall) then
C        r3 = rr**3
C        rk = -1d0
C        rj = 1d0/rr
C        do  20  l = 0, lmax
C          rk = rk*(2*l-1)/rr
C          rj = rj*rr/(2*l+1)
C          fkk(l) = rk*rk*r3/(2*l-1)
C          fjj(l) = -rj*rj*r3/(2*l+3)
C          fkj(l) = -0.5d0*rj*rk*r3
C          fjk(l) = fkj(l)
C   20   continue
C        return
C      endif
C      if (dabs(e1-e2) > esmall) then
CC ---   Case e1 /= e2 ---
C        efac = 1/(e2-e1)
C        erfac = efac*rr**2
C        if (r < 0) efac = 0
C        job = 0
C      else
CC ---   Case e1 == e2 but not zero ---
C        erfac = rr**2
C        efac = 0
C        job = 1
C      endif
C      call radkj(e1,rr,lmax,ak1,aj1,dk1,dj1,0)
C      call radkj(e2,rr,lmax,ak2,aj2,dk2,dj2,job)
C      do  10  l = 0, lmax
C        fkk(l) = erfac*(ak1(l)*dk2(l) - dk1(l)*ak2(l))
C        fjj(l) = erfac*(aj1(l)*dj2(l) - dj1(l)*aj2(l))
C        fkj(l) = erfac*(ak1(l)*dj2(l) - dk1(l)*aj2(l)) - efac
C        fjk(l) = erfac*(aj1(l)*dk2(l) - dj1(l)*ak2(l)) + efac
C   10 continue
C      end
C      subroutine radkj(e,r,lmax,ak,aj,dk,dj,job)
CC- Radial parts of spherical hankels and bessels.
CC ----------------------------------------------------------------
CCi Inputs
CCi   e,r,lmax
CCi   job:  0, makes values and slopes; 1, makes energy derivatives.
CCo Outputs
CCo   ak,aj,dk,dj for l=0..lmax (dk is partial ak / partial r)
CCr Remarks
CCr   Energy derivative does not work with OKA'S defs!
CC ----------------------------------------------------------------
C      implicit none
CC Passed parameters
C      integer Job,Lmax
C      double precision E,R
C      double precision ak(lmax),aj(lmax),dk(lmax),dj(lmax)
CC Local parameters
C      integer l,lp1,lmxx
C      parameter (lmxx=12)
C      double precision er2,rl
C      double precision phi(lmxx+2),psi(lmxx+2),php(lmxx+2),psp(lmxx+2)
C      external bessl
C
C      if (lmax > lmxx) call rx('radkj: increase lmxx')
C      er2 = e*r**2
C      if (job == 0) then
C        call besslm(er2,lmax+1,phi,psi)
C        rl = 1.d0/r
C        do  10  l = 0, lmax
C          lp1 = l+1
C          rl = rl*r
C          ak(lp1) = psi(lp1)/(rl*r)
C          aj(lp1) = phi(lp1)*rl
C          dk(lp1) = (l*psi(lp1) - psi(l+2))/(rl*r**2)
C          dj(lp1) = (l*phi(lp1) - er2*phi(l+2))*rl/r
C   10   continue
C      else
C        call besslm(er2,lmax+2,phi,psi)
C        do  11  lp1 = 1, lmax+2
C          php(lp1) = -0.5d0*r**2*phi(lp1+1)
C          psp(lp1) = +0.5d0*((2*lp1-1)*psi(lp1) - psi(lp1+1))/e
C   11   continue
C        rl = 1.d0/r
C        do  12  l = 0, lmax
C          lp1 = l+1
C          rl = rl*r
C          ak(lp1) = psp(lp1)/(rl*r)
C          aj(lp1) = php(lp1)*rl
C          dk(lp1) = (l*psp(lp1) - psp(l+2))/(rl*r*r)
C          dj(lp1) = (l*php(lp1) - er2*php(l+2) - r*r*phi(l+2))*rl/r
C   12  continue
C      endif
C      end
