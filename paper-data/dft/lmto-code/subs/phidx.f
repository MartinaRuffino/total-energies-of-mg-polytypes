      subroutine phidx(job,z,l,v,hcr,vmtz,rofi,nr,nptdif,tol,
     .  e,val,slo,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
C- Generate potential parameters for a prescribed energy or b.c.
C ----------------------------------------------------------------
Ci Inputs:
Ci   job:   1s digit specifies boundary conditions
Ci          0  boundary conditions specified val,slo,nn (see Remarks)
Ci          1  same as 0, but also g and e assumed to be generated
Ci          2  boundary conditions specified by energy e (see Remarks)
Ci          3  Same as 2, but radial integration only outward (see Remarks)
Ci             In this case val,slo,nn are output
Ci          4  Same as 3, but sign of g(nr) fixed to match sign of input val(1)
Ci          10s digit
Ci          1  set dphip to satisfy Wronskian condition
Ci          100s digit
Ci          1  allow finite difference dele to shrink if slope
Ci             is different sign within finite difference.
Ci             Appropriate for deep states
Ci   z:     nuclear charge
Ci   l:     l quantum number for this g
Ci   v:     spherical potential on shifted logarithmic mesh
Ci   vmtz   flat potential used to generate dla,dlap
Ci          Not used if hcr=0
Ci   hcr:   hard sphere radius.  If nonzero, dla,dlap are generated
Ci   rofi:  list of points
Ci   nr:    number of mesh points
Ci   nptdif:2 or 4 for 3- or 5- point differentiation
Ci         :You may also set nptdif=0.  Then quantities related to
Ci         :energy differences are not calculated (dlap,phip,dphip,p)
Ci   tol:   precision to which wave function is integrated
Cio Inputs/Outputs:
Cio         Which subset of these quantities are needed for input and
Cio         which quantities phidx alters depends on job; see Remarks.
Cio
Cio  e:     On input (job=1,2,3,4), energy eigenvalue
Cio         On output (job=0), energy eigenvalue
Cio  val:   On input (job=0,1), val(1) = value of g(r) at rmax
Cio         On output (job=0,2,3,4), val(1) = value of g(r) at rmax, g normalized
Cio         Also on output, val(2..1+nptdif) = energy derivatives of val
Cio  slo:   On input (job=0,1), slo(1) = dg(r)/dr at rmax
Cio         On output (job=0,2), slo(1) = dg(r)/dr at rmax, g normalized
Cio         Also on output, slo(2..1+nptdif) = energy derivatives of slo
Cio  nn:    On input (job=0,1), number of nodes
Cio         On output (job=2,3,4), number of nodes
Cio  g:     On input (job=1) r * radial part of partial wave
Cio         (assumed normalized so that int (g*g) dr = 1)
Cio         On output (job=0,2,3,4) normalized wave function times r
Co  Outputs:
Co   gp:    first nptdif energy derivatives to G
Co   phi:   true partial wave at rmax, i.e. g/rmax
Co   dphi  :slope of true partial wave at rmax, i.e. (d(g/r)/dr)_rmax
Co   phip:  energy derivative of true partial wave at rmax
Co   dphip: energy derivative of slope of true partial wave at rmax
Co   p:     <gp**2> (potential parameter)
Co   phia:  (hcr>0) value of phi at hcr boundary, i.e. g(hcr)/hcr
Co   phiap: (hcr>0) energy derivatives of phia
Co   dla:   (hcr>0) hcr * logarithmic derivative of phi0 at hcr boundary
Co                  where phi0 is back-extrapolated true partial wave
Co          (hcr=0) not calculated
Co   dlap:  (hcr>0) energy derivatives of dla
Co          (hcr=0) not calculated
Cr Remarks:
Cr   This version makes parameters related to true partial wave g(r)/r
Cr   defined by potential v(r).
Cr
Cr   Boundary conditions either by (val,slo) or be energy, according to 1s digit job
Cr     job=0   val,slo,nn are specified.  On output,
Cr             val,slo are renormalized so that val=g(nr), with <gg>=1
Cr             Here energy eigenvalue e is calculated
Cr             are assumed to correspond with g.)
Cr     job=1   Assumes that all quantities val,slo,nn,e,g have been
Cr             generated, and calculates none of them.
Cr     job=2   the energy eigenvalue e is specified.
Cr             val,slo, and nn are calculated.
Cr     job=3   the energy eigenvalue e is specified.
Cr             val,slo, and nn are calculated by outward radial integration only
Cl
Cl Local variables
Cl   g(r,1) : large component of g.  psi = (g(r)/r)*ylm
Cl
Cb Bugs
Cb   calls info ... should not do this since may be called in parallel!
Cb   Not really a bug, but phidx returns redundant information in
Cb   the following variables:
Cb      phi   = val(1)/rmax
Cb      dphi  = (slo(1) - phi)/rmax = d(g/r)/dr = [dg/dr - 1/r g] / r  at rmax
Cb      phip  = vali(1)/rmax
Cb      dphip = (sloi(1) - phip)/rmax
Cu Updates
Cu   23 Mar 17 New option to compute phi,phidot from outward radial integration only
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   11 Jan 10 Fix problem in num difference for very deep valence states
Cu   21 Jul 04 possibly set dphip to satisfy Wronskian condition
Cu   19 Dec 01 Return val,slo,phiap,dlap for nptdif energy derivatives
Cu    7 May 01 Allow nptdif=0
Cu   29 Nov 00 printout when phidx fails to converge, softer tole
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,l,nr,nn,nptdif
      double precision z,e,vmtz,hcr,val(*),slo(*),phi,dphi,phip,dphip,
     .  dla,dlap(*),p,phia,phiap(*),tol,v(nr),rofi(nr),
     .  g(2*nr),gp(2*nr,4)
C ... Local parameters
C     logical :: debug=.false.
      integer kc,nre,i,iprint,stdo,lgunit,job0,job1,job2
      double precision rmax,eb1,eb2,dele,ddde,sum,a,b,aold,dmach,tola,
     .  sloi(5),vali(5),phiai(4),dlai(4),ei(4),de1,de2,del1,del2,xx
      real(8) :: gsign  ! possibly scale g by sign (job=4)
      real(8), parameter :: tole=1d-10
C External calls
      external dfphi,gintsr,iprint,makdla,rseq,rsq1,rx,lgunit

C ... Setup: extract a and b
      job0 = mod(job,10)
      job1 = mod(job/10,10)
C     if (job1/=0) call rx('oops')
      job2 = mod(job/100,10)
      stdo = lgunit(1)
      rmax = rofi(nr)
      dele = .002d0
      a = log(rofi(nr)/rofi(nr-1))
      tola = 8*dmach(1)
      gsign = 1
      do   i = 1, 100
        aold = a
        a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
        if (i > 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
        if (abs(a-aold) <= tola) exit
        if (i < 100) cycle
        call rx('phidx failed to determine ''a'' parameter')
      enddo
      b = rmax/(dexp(a*(nr-1)) - 1d0)

C --- Find energy e corresponding to val,slo ---
      if (job0 == 0) then
        eb1 = -20d0
        eb2 =  20d0
        e = (eb1+eb2)/2
C       This generates g, normalized to <g g> = 1
        call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre,kc)
C       Scale val, slo to correspond to normalization <g g> = 1
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)

C --- Find val,slo corresponding to energy e ---
      elseif (job0 >= 2 .and. job0 <= 4) then
C       Initial estimate for val,slo.  Final result for job=3,4, except for normalization
        xx = val(1)
        call rsq1(0,e,l,z,v,nr,g,val,slo,nn,a,b,rofi,nr)
        if (job0 == 4 .and. xx*val(1) < 0) gsign = -1 ! Change sign of normalization
        if (job0 == 3 .or. job0 == 4) then ! Normalize
          call gintsr(g,g,a,nr,z,e,l,v,rofi,sum)
          call dscal(2*nr,gsign/dsqrt(sum),g,1)
        else
C       Adjust slope iteratively until ei(slope)-e is within tole
          ei(1) = e
          eb1 = e-.1d0
          eb2 = e+.1d0
          do  i = 1, 5+5
            if (i == 5+5) call rx('phidx failed to converge')
            call rseq(eb1,eb2,ei,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre,kc)
            if (iprint() > 0 .and. i > 5) then
              call awrit7(' PHIDX  Z=%d  l=%i  nod=%i  bc=%;4g %;4g'//
     .          '  e(bc)=%;4g  e(bc)-e=%;4g',' ',80,stdo,z,l,nn,val(1),
     .          slo(1),ei(1),ei(1)-e)
            endif
            if (abs(ei(1)-e) < tole) exit
            slo(1) = slo(1) + (ei(1)-e) * val(1) / g(nr)**2
          enddo
        endif
        val(1) = val(1)/dsqrt(sum) * gsign
        slo(1) = slo(1)/dsqrt(sum) * gsign
      elseif (job0 /= 1) then
        call rx('phidx: bad job')
      endif
      if (hcr /= 0) call makdla(e-vmtz,l,hcr,slo,val,rmax,phia,dla)

C --- Energy derivatives ---
      if (nptdif /= 0) then
        ddde = -rmax/g(nr)**2
   20   continue
        ei(1) = 1
        ei(2) = -1
        ei(3) = 1.5d0
        ei(4) = -1.5d0
C       Case finite diff should be small enough that slo preserves sign
C       For deep states, avoids change in node count
        if (job2 == 1) then
          do  i = 1, nptdif
            sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
C           Energy spread too wide: change of nodes.  Reduce de
            if (sloi(i)*slo(1) < 0 .and. dele > 1d-7) then
              dele = dele/4
              call info(50,0,0,'phidx l=%i: warning, scale dele to %g',l,dele)
C             if (dele < 1d-8) call rx('phidx dele too small.')
              goto 20
            endif
          enddo
        endif

C   ... Construct g at mesh of energies
        eb1 = e-.1d0
        eb2 = e+.1d0
        do  i = 1, nptdif
C         vali(i) = val(1)  ! This causes a problem for the intel ifort 15.0.2 compiler
          sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
          ei(i) = e + dele*ei(i)
          if (job0 == 3 .or. job0 == 4) then
            call rsq1(0,ei(i),l,z,v,nr,gp(1,i),vali(i),sloi(i),nn,a,b,rofi,nr)
            call gintsr(gp(1,i),gp(1,i),a,nr,z,e,l,v,rofi,sum)
            call dscal(2*nr,gsign/dsqrt(sum),gp(1,i),1)
            vali(i) = vali(i)/dsqrt(sum) * gsign
            sloi(i) = sloi(i)/dsqrt(sum) * gsign
          else
            call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
     .                sum,a,b,rofi,nr,nre,kc)
            vali(i) = val(1)/dsqrt(sum)
            sloi(i) = sloi(i)/dsqrt(sum)
          endif
          if (hcr /= 0) call makdla(ei(i)-vmtz,l,hcr,sloi(i),vali(i),
     .      rmax,phiai(i),dlai(i))

        enddo

C        if (debug) then
C          print *, 'before differentiation, l',l
C          call prrmsh('v',rofi,v,nr,nr,1)
C          call prrmsh('g(0)',rofi,g,nr,nr,2)
C          call prrmsh('g(ei)',rofi,gp,nr,nr,8)
C        endif

        de1  = (ei(1) - ei(2))/2
        del1 = (ei(1) + ei(2))/2 - e
        de2  = (ei(3) - ei(4))/2
        del2 = (ei(3) + ei(4))/2 - e
C       Energy derivatives of value and slope
        call dfphi(de1,del1,de2,del2,1,val,vali,nptdif == 4)
        call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif == 4)
C       Energy derivatives of dla
        if (hcr /= 0) then
          call dfphi(de1,del1,de2,del2,1,dla,dlai,nptdif == 4)
          call dfphi(de1,del1,de2,del2,1,phia,phiai,nptdif == 4)
            do  i = 1, nptdif
            dlap(i)  = dlai(i)
            phiap(i) = phiai(i)
            enddo
        endif
C       Energy derivatives of g
        call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif == 4)
C       p = integral <gp gp>
        call gintsr(gp,gp,a,nr,z,e,l,v,rofi,p)
      endif

C      if (debug) then
C        print *, 'after differentiation'
C        call prrmsh('gdot,gdotdot',rofi,gp,nr,nr,4)
C      endif

C     phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
      phi = val(1)/rmax
      dphi = (slo(1) - phi)/rmax  ! d(g/r)/dr = [dg/dr - 1/r g] / r
      if (nptdif /= 0) then
        phip = vali(1)/rmax
        dphip = (sloi(1) - phip)/rmax
      endif

C     Copy energy derivatives sloi to slo(2..)
      if (nptdif /= 0) then
        call dcopy(nptdif,sloi,1,slo(2),1)
        call dcopy(nptdif,vali,1,val(2),1)
      endif

C     Set dphip from Wronskian condition:
C     phi*dphip - dphi*phip = -1/rmax**2 =>
C     dphip = (dphi*phip - 1/rmax**2)/phi
      if (nptdif /= 0 .and. job1 /= 0) then
C        print *, dphip,(dphi*phip - 1/rmax**2)/phi,
C     . dphip - (dphi*phip - 1/rmax**2)/phi
        dphip = (dphi*phip - 1/rmax**2)/phi
      endif

      if (iprint() >= 111) write (stdo,1) phi,dphi,phip,dphip
    1 format(' PHIDOT:  phi,phip,phip,dphip=',4F12.6)

      end
