      subroutine convolve(mode,gc,ef,e0,ep,sigma,n,x,f,fc)
C- Convolves a function with Gaussian or Lorentzian
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 convolve with Lorentzian
Ci         :2 convolve with Gaussian
Ci   gc    :Gamma_c =  Core width    (Lorentzian broadening)
Ci   ef    :E_F     =  Fermi energy  (Lorentzian broadening)
Ci   e0    :E_0     =  Bottom of CB  (Lorentzian broadening)
Ci   ep    :E_p     =  Plasma energy (Lorentzian broadening)
Ci   sigma :sigma   =  Instrumental Gaussian width (Gaussian broadening)
Ci   n     :number of points
Ci   x     :x(1:n) = energy (abscissa in general)
Ci   f     :f(1:n) = function
Co Outputs
Co   fc    :fc(1:n-1) convolved function
Cl Local variables
Cl   pf    :(for Lorentzian broadening)
C          :final state width, calculated relative to bottom of CB
Cl         :if ep != 0, pf = pi*pi*sqrt(3)*ep/(128*(ef-e0)*(ef-e0))
Cr Remarks
Cr
Cu Updates
Cu   16 Apr 10 Adapted from Paxton's convolute.c
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,n
      double precision gc,ef,e0,ep,sigma,x(n),f(n),fc(n)
C ... Local parameters
      integer i,j,bin1
      double precision pi,twopi,pf,sum,g,f1,f2,s2,xx

      if (mode == 0) return
      pi = 4*datan(1d0)
      twopi = 2*pi

      if (mode == 1) then
        bin1 = 1
        do while (x(bin1) < ef)
          bin1 = bin1+1
        end do
        if (ep /= 0d0) then
          pf = pi*pi*dsqrt(3d0)*ep/(128d0*(ef-e0)*(ef-e0))
        else
          pf = 0d0
        endif
        do  i = 1, n-1
          sum = 0
          do  j = bin1, n-1
            g = gc + pf * (x(j) - ef)*(x(j) - ef)
            f1 = g*f(j)/((x(i)-x(j))*(x(i)-x(j)) + 0.25d0*g*g)
            g = gc + pf * (x(j+1)-ef)*(x(j+1) - ef)
            f2 = g*f(j+1)/((x(i)-x(j+1))*(x(i)-x(j+1)) + 0.25d0*g*g)
            sum = sum + 0.5d0*(f1+f2)*(x(j+1)-x(j))
          end do
          fc(i) = sum/twopi
        end do
      elseif (mode == 2) then
        s2 = sigma*sigma
        do  i = 1, n-1
          sum = 0
          do  j = 1, n-1
            xx = x(i)-x(j)
            f1 = f(j) * dexp(-0.5d0*xx*xx/s2)
            xx = x(i)-x(j+1)
            f2 = f(j+1) * dexp(-0.5d0*xx*xx/s2)
            sum = sum + 0.5d0*(f1+f2)*(x(j+1)-x(j))
          end do
          fc(i) = sum/(sigma*dsqrt(twopi))
        end do
      else
        call rx('bad mode')
      endif
      end
C#ifdefC DEBUG
C      subroutine fmain
C      integer ifi,fopng,rdm,i
C      real(8),allocatable:: dat(:,:),x(:),f(:),fc(:)
C      double precision xx,gc,ef,ep,sigma
C
C      ifi = fopng('dosr.inp',-1,0)
C      nr = 0
C      nc = 0
C      xx = -1
C      i = rdm(ifi,0,0,' ',xx,nr,nc)
C      if (i /= 1)
C     .  call rx('file dosr.inp does not contain a real matrix')
C      if (nr < 2)
C     .  call rx('file dosr.inp has < 2 rows')
C      if (nc < 2)
C     .  call rx('file dosr.inp should have 2 columns')
C
C      allocate(dat(nr,nc),x(nr),f(nr),fc(nr))
C      rewind ifi
C      i = rdm(ifi,0,nr*nc,' ',dat,nr,nc)
C      x(:) = dat(:,1)
C      f(:) = dat(:,2)
C      gc = .01d0
C      ef = -.56d0
C      e0 = -.6d0
C      ep = 0.5d0*0
C      sigma = 0
CC     call convolve(1,gc,ef,e0,ep,sigma,nr,x,f,fc)
C
C      gc = 0
CC      ef = -.56d0
CC      e0 = -.6d0
CC      ep = 0.5d0*0
C      sigma = 0.02d0
C      call convolve(2,gc,ef,e0,ep,sigma,nr,x,f,fc)
C
C      dat(:,2) = fc(:)
C      ifi = fopng('out',-1,0)
C      call ywrm(0,' ',1,ifi,'(3f14.8)',dat,0,nr,nr-1,2)
C
C      end
