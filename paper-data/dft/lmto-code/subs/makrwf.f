      subroutine makrwf(mode,z,rmax,l,v,a,nr,rofi,pnu,mu,nptdif,g,gp,
     .  enu,phi,dphi,phip,dphip,p)
C- Radial wave functions and energy derivative
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit specifies boundary conditions
Ci         :0 input boundary conditions
Ci         :1 input energy
Ci         :10s digit
Ci         :1 set dphip to satisfy Wronskian condition, rather than compute numerically
Ci         :100s digit for phidot
Ci         :1, allow finite difference dele to shrink if slope
Ci         :   is different sign within finite difference.
Ci         :   Appropriate for deep states
Ci         :1000s digit for Dirac partial waves, if calculated
Ci         :1  Do not compute Dirac energy derivative (difficult for core-like states)
Ci         :
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   l     :l-quantum number
Ci   v     :spherical potential
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   mu    :If half-integer and |2*mu|<l+1,
Ci         :replace g with full Dirac rwf, alpha=1, no B-field
Co Outputs
Co   g     :r * wave function corresponding to b.c. pnu
Co   gp    :r * energy derivative of g; dimensioned 8*nr
Co   enu   :eigenvalue of g
Co   phi   :wave function at rmax, i.e. g/rmax
Co   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
Co   phip  :energy derivative of wave function at rmax
Co   dphip :energy derivative of slope of wave function at rmax
Co   p     :<gp**2> (potential parameter)
Cr Remarks
Cr   This routine makes r*phi and r*phidot, where phi and phidot
Cr   are true radial wave function and energy derivatives.
Cr   phi is normalized, and p = <phidot**2>
Cb Bugs
Cb   calls info ... should not do this since may be called in parallel!
Cu Updates
Cu   11 Jan 10 Patch phidx call for deep semicore states
Cu   16 Aug 04 10s digit to explicitly satisfy Wronskian
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   16 May 00 New routine
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: mode,l,nr,nptdif
      real(8), intent(in) :: a,rmax,z,mu,rofi(nr),v(nr),pnu
      real(8), intent(out) :: g(nr,2),gp(nr,2,4),phi,phip,dphi,dphip,p,enu
C ... Local parameters
      integer :: konf,nn,nre,modep,mode3,kc,imu,nsp,intopt,i
      real(8) :: b,dnu,eb1,eb2,pi,slo(5),sum,tol,val(5)
      real(8) :: gsmt(2,2,2),vsnr(2,2),gfn(2,2,2,2)
      real(8), allocatable :: psi(:,:,:,:),psd(:,:,:,:),riw(:,:),bf(:)
      procedure(integer) :: nglob,getdig

      pi = 4d0*datan(1d0)
      tol = 1d-12
C     for now, rmax must match to rofi(nr)
      call fsanrg(rmax,rofi(nr),rofi(nr),1d-8,'makrwf:','rmax',.true.)
      if (mod(mode,10) == 0) then
        b   = rmax/(exp(a*nr-a)-1d0)
        konf = mod(pnu,10d0)
        dnu = tan(pi*(.5d0-mod(pnu,10d0)))
        nn = konf-l-1
        val(1) = rmax
        slo(1) = dnu+1d0
        eb1 = -20d0
        eb2 =  20d0
        if (z == 1 .and. l > 2) eb2 = 100
        enu = -0.5d0
        call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre,kc)
        val(1) = val(1)/dsqrt(sum)  ! So that val = g(nr,1)
        slo(1) = slo(1)/dsqrt(sum)  ! and slo = dg/dr at nr

C        call phidot(z,l,v,enu,a,b,rofi,nr,g,val,slo,tol,nn,
C     .    gp,phi,dphi,phip,dphip,p)
        modep = 1
      else
        modep = 2
      endif

      if (mod(mode/100,10) == 1) modep = modep + 100
      modep = modep + 10*getdig(mode,1,10)

      call phidx(modep,z,l,v,0d0,0d0,rofi,nr,nptdif,tol,enu,val,slo,
     .  nn,g,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
C      dphip = (sloi(2)-phip)/rmax

      sum = abs(mu)
      if (pnu > 100 .and. sum <= l+0.5d0 .and. nint(2*sum) == 2*sum .and. mod(nint(2*sum),2) == 1)  then
        imu = mu + l + 1.5d0
        nsp = 1
        mode3 = mod(mode/1000,10)
        intopt = 10*nglob('lrquad')
        allocate(riw(nr,2),psi(nr,2,2,2),psd(nr,2,2,2),bf(nr))
        call dpzero(bf,nr)
        call dcopy(nr,rofi,1,riw,1)
        call radwgt(intopt,rofi(nr),a,nr,riw(1,2))

        call dpzero(gfn,size(gfn))
        gfn(1,1,1,1) = g(nr,1)    ! large component from scalar rel, alpha=1
        gfn(2,1,1,1) = g(nr,2)    ! small component from scalar rel, alpha=1
        gfn(1,2,1,2) = g(nr,1)    ! large component from scalar rel, alpha=2
        gfn(2,2,1,2) = g(nr,2)    ! small component from scalar rel, alpha=2

        gfn(1,1,2,1) = gp(nr,1,1) ! Corresponding components for dot
        gfn(2,1,2,1) = gp(nr,2,1)
        gfn(1,2,2,2) = gp(nr,1,1)
        gfn(2,2,2,2) = gp(nr,2,1)

        call info8(1,0,0,'  sr enu%;12,6D  phi, dphi%;12,6D%;12,6D   phip, dphip%;12,6D%;12,6D  p%;12,6D',
     .    enu,phi,dphi,phip,dphip,p,7,8)

        i = 0; if (pnu > 200) i = i+12
        if (mod(mode3/2,2) == 1 .or. .true.) i = i+100

        call rdeqmu(enu,gfn,i,z,v,riw,nr,nsp,kc,a,b,l,-1,imu,psi,psd,gsmt,g) ! g is dummy

        if (g(nr,1)*psi(nr,1,1,1) < 0) then
          call dscal(size(psi),-1d0,psi,1)
          call dscal(size(psd),-1d0,psd,1)
        endif

        call rdeqvs(1,l,imu,nr,nr,enu,rofi,z,v,bf,psi,gsmt,vsnr)
        phi  = vsnr(1,1)/rofi(nr)
        dphi = vsnr(1,2)/rofi(nr) - vsnr(1,1)/rofi(nr)**2
        if (mod(i/100,10) == 0) then
          call rdeqvs(1,l,imu,nr,nr,enu,rofi,z,v,bf,psd,gsmt,vsnr)
          phip  = vsnr(1,1)/rofi(nr)
          dphip = vsnr(1,2)/rofi(nr) - vsnr(1,1)/rofi(nr)**2
          call prodrw(nr,4,psd,psd,riw(1,2),p)
        endif

C       call info8(1,0,0,'  fr enu%;12,6D  phi, dphi%;12,6D%;12,6D   phip, dphip%;12,6D%;12,6D  p%;12,6D',
        call info8(1,0,0,'  fr enu%;12,6D  phi, dphi%;12,6D%;12,6D',
     .    enu,phi,dphi,phip,dphip,p,7,8)
C       call prrmsh('g',rofi,g,nr,nr,2)
C       call prrmsh('psi',rofi,psi,nr,nr,8)

        call dcopy(2*nr,psi,1,g,1)
        call dcopy(2*nr,psd,1,gp,1)
      endif

C     call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)

      end
