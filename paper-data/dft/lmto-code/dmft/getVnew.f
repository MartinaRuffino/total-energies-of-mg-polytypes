      subroutine getVnew(V,Nref,N,derN,min_req,rfwk,nit)
C- Controls iterative search to find chemical potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   Nref  :target number of electrons
Ci   N     :Number of electrons for this V
Ci   derN  :dN/dV = -DOS
Ci   rfwk  :work array for rfalsi
Cio Inputs/Outputs
Ci    V    :chemical potential is ef0 - V
Ci         :On input trial V that generated N
Ci         :On output estimate for new V
Cio   nit  :number of prior iterations; incremented by 1 on output
Co Outputs
Co  min_req:F when converged to tolerance
Cl Local variables
Co  grfac  :Growth factor extrapolation a bit farther than linear when ir is returned -4
Co         :Accelerates bracketing of root with minor loss of precision
Cr Remarks
Cr   This subroutine initializes quantities to be used  by rfalsi,
Cr   which actually performs the search of the root of N-Nef.
Cr   The simplest estimate for V is:
Cr   V = (N-Nref)/dos = -(N-Nref)/derN
Cu Updates
Cu   22 Aug 15 (L. Sponza) first created
C ----------------------------------------------------------------------
      implicit none
      real(8), intent(in)    :: Nref,N,derN
      real(8), intent(inout) :: V
      real(8), intent(inout) :: rfwk(12) ! working vector for rfalsi
      integer, intent(inout) :: nit
      logical, intent(out)   :: min_req
C ... Local parameters
      integer :: isw            ! rfalsi convergence criterion
      integer :: ir             ! switch for rfalsi
      real(8) :: dvcap=.2d0     ! cap on delta V estimated from DOS
      real(8) :: deltaN         ! number of excess electrons
      real(8) :: dxmx           ! max step for rfalsi
      real(8),parameter :: grfac = 1.1d0
      procedure(integer) :: iprint

      ir = nint(rfwk(12))
      deltaN = N-Nref
      if (ir == 0) then
        dxmx = 0
        if (derN /= 0) then
          dxmx = -deltaN/derN
        endif
        if (abs(dxmx) > dvcap) dxmx = sign(dvcap,dxmx)
      else
       dxmx = dvcap
      endif

      call pshpr(max(1,iprint()-20))
C     isw = 30 => tolerance in EITHER deltaN or change in V
C     isw =  0 => tolerance in deltaN only
      isw = 30
      call rfalsi(V,deltaN,5d-10,1d-7,5d-10,dxmx,isw,rfwk,ir)
      call poppr
      if (ir == -4 .and. abs(v-rfwk(1)) < dvcap/2) then
        v = rfwk(1) + (v-rfwk(1))*grfac
      endif

      if (ir /= -1) then
        call info8(30,0,0,' getVnew: v1=%1;4e N1=%1;4e  v2=%1;4e N2=%1;4e  bracket=%l  est=%1;6g',
     .    rfwk(2),rfwk(5),rfwk(1),rfwk(4),(rfwk(5)*rfwk(4)<0),V,0,0)
      endif
      if (ir == 0 .or. ir == 1) min_req = .false.
      rfwk(12) = ir
      nit = nit+1
      end subroutine getVnew
