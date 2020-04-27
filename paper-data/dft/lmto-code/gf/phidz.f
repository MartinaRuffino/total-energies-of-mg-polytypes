      subroutine phidz(job,z,l,v,hcr,vmtz,rofi,nr,nptdif,tol,
     .  e,val,slo,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
C- Generate potential parameters for a prescribed energy
C ----------------------------------------------------------------
Ci Inputs:
Ci   job:   0, boundary conditions specified val,slo,nn (see Remarks)
Ci          1, same as 0, but also g and e assumed to be generated
Ci          2, boundary conditions specified by energy e (see Remarks)
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
Cio         which quantities phidz alters depends on job; see Remarks.
Cio
Cio  e:     On input (job=1,2), energy eigenvalue
Cio         On output (job=0), energy eigenvalue
Cio  val:   On input (job=0,1), val(1)=value of g(r)=r*phi(r) at rmax
Cio         On output (job=0,2), val(1)=value of normalized g(r) at rmax
Cio         Also on output, val(2..1+nptdif) = energy derivatives of val
Cio  slo:   On input (job=0,1), slo(1)=radial derivative of g(r) at rmax
Cio         On output (job=0,2), slo(1)=der. of normalized g(r) at rmax
Cio         Also on output, slo(2..1+nptdif) = energy derivatives of slo
Cio  nn:    On input (job=0,1), number of nodes
Cio         On output (job=2), number of nodes
Cio  g:     On input (job=1) Wave function times r
Cio         (assumed normalized so that int (g*g) dr = 1)
Cio         On output (job=0,2) normalized wave function times r
Co  Outputs:
Co   gp:    first nptdif energy derivatives to G
Co   phi:   wave function at rmax, i.e. g/rmax
Co   dphi:  slope of wave function at rmax, i.e. d(g/rmax)/dr
Co   phip:  energy derivative of wave function at rmax
Co   dphip: energy derivative of slope of wave function at rmax
Co   p:     <gp**2> (potential parameter)
Co   phia:  (hcr>0) value of phi0 at hcr boundary, where phi0 is
Co          back extrapolated function
Co   phiap: (hcr>0) energy derivatives of phia
Co   dla:   (hcr>0) hcr * logarithmic derivative of phi0 at hcr boundary
Co                  where phi0 is back-extrapolated wave function
Co          (hcr=0) not calculated
Co   dlap:  (hcr>0) energy derivatives of dla
Co          (hcr=0) not calculated
Cr Remarks:
Cr   phidz is an analog of phidx for complex energies, which see.
Cr
Cr  *Case job=0,1:
Cr   In this mode, the number of nodes nn is an input.  A 'node' is
Cr   defined when the real part of the wave function changes sign.
Cr   Because nodes are not necessarily meaningful for complex energies,
Cr   caller can pass nn<0 and no check is made for the number of nodes.
Cr   In that case caller is advised to call phidz with an energy
Cr   reasonably close to the one matching val,slo.
Cb Bugs
Cb   Not really a bug, but phidx returns redundant information in
Cb   the following variables:
Cb      phi   = val(1)/rmax
Cb      dphi  = (slo(1) - phi)/rmax
Cb      phip  = vali(1)/rmax
Cb      dphip = (sloi(1) - phip)/rmax
Cu Updates
Cu   19 Dec 01 Return val,slo,phiap,dlap for nptdif energy derivatives
Cu   10 Dec 01 First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,l,nr,nn,nptdif
      double precision z,vmtz,hcr,tol,v(nr),rofi(nr)
      double complex e,val(*),slo(*),phi,dphi,phip,dphip,dla,dlap(*),p,
     .  g(2*nr),gp(2*nr,4),phia,phiap(*)
C ... Local parameters
      logical lagain
      integer nre,i,j,iprint,stdo,lgunit,nnl,kc
      double precision rmax,eb1,eb2,dele,sum,a,b,aold,dmach,tola,tole,
     .  fac,diff1,diff2
      double complex ddde,sloi(4),vali(4),phiai(4),dlai(4),ei(4),de1,
     .  de2,del1,del2
      parameter (tole=1d-10)
C External calls
      external dfphiz,gintsr,iprint,makdlz,rseq,rx,lgunit

C ... Setup: extract a and b
      stdo = lgunit(1)
      rmax = rofi(nr)
      dele = .002d0
      a = log(rofi(nr)/rofi(nr-1))
      tola = 8*dmach(1)
      do   i = 1, 100
      aold = a
      a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
      if (i > 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
      if (abs(a-aold) <= tola) goto 1
      enddo
      call rx('phidz failed to determine ''a'' parameter')
    1 continue
      b = rmax/(dexp(a*(nr-1)) - 1d0)

C --- Find energy e corresponding to val,slo ---
      if (job == 0) then
        eb1 = -20d0
        eb2 =  20d0
        e = (eb1+eb2)/2
C       Generate g and e, for Re val, Re slo normalized to <g g> = 1
        if (nn >= 0)
     .    call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,
     .    nre,kc)
        call zseq(e,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
C       Scale val, slo to correspond to normalization <g g> = 1
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)

C --- Find val,slo corresponding to energy e ---
      elseif (job == 2) then
C       Initial estimate for val,slo
        call zsq1(e,l,z,v,nr,g,val,slo,nnl,a,b,rofi,nr)
        if (nn >= 0) nn = nnl
C       Adjust slope iteratively until ei(slope)-e is within tole
        ei(1) = e
C        eb1 = e-.1d0
C        eb2 = e+.1d0
        j = 5
        if (tol <= 1d-12) j = 7
        fac = 1d0
        diff1 = 0
C       Iterate to find energy with inward-outward matching of b.c.
C       lagain = T if delta e < 1/2 that of 2nd iter prior to current
        i = 0
    5   continue
        i = i+1
        diff2 = diff1
        diff1 = abs(ei(1)-e)
        call zseq(ei,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
        lagain = abs(ei(1)-e) < diff2/2 .or. i <= 2
        if (iprint() > 10 .and. i > j .and. .not. lagain) then
          call awrit8(' PHIDZ  i=%i  Z=%d  l=%i  nod=%i  bc=%2:1;4g'//
     .      ' %2:1;4g  e(bc)=%2:1;4g  diffe=%2:1;4g',' ',120,stdo,i,z,
     .      l,nn,val(1),slo(1),ei(1),(ei(1)-e))
C         fac = .6d0
        endif
        if (abs(ei(1)-e) < tole) goto 2
C       Estimate change in slope that corrects for energy mismatch
        slo(1) = slo(1) + (ei(1)-e) * val(1) / g(nr)**2 * fac
        if (lagain .or. i < j+j) goto 5
        call rx('phidz failed to converge')
    2   continue
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)
      elseif (job /= 1) then
        call rx('phidz: bad job')
      endif
      if (hcr /= 0) call makdlz(e-vmtz,l,hcr,slo,val,rmax,phia,dla)

      if (nptdif /= 0) then
      ddde = -rmax/g(nr)**2
      ei(1) = 1
      ei(2) = -1
      ei(3) = 1.5d0
      ei(4) = -1.5d0
      eb1 = e-.1d0
      eb2 = e+.1d0
      do  10  i = 1, nptdif
        sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
        ei(i) = e + dele*ei(i)
        call zseq(ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
     .            sum,a,b,rofi,nr,nre)
        vali(i) = val(1)/dsqrt(sum)
        sloi(i) = sloi(i)/dsqrt(sum)
        if (hcr /= 0) call makdlz(ei(i)-vmtz,l,hcr,sloi(i),vali(i),
     .    rmax,phiai(i),dlai(i))

   10 continue

      de1  = (ei(1) - ei(2))/2
      del1 = (ei(1) + ei(2))/2 - e
      de2  = (ei(3) - ei(4))/2
      del2 = (ei(3) + ei(4))/2 - e
C     Energy derivatives of value and slope
      call dfphiz(de1,del1,de2,del2,1,val,vali,nptdif == 4)
      call dfphiz(de1,del1,de2,del2,1,slo,sloi,nptdif == 4)
C     Energy derivatives of dla
      if (hcr /= 0) then
        call dfphiz(de1,del1,de2,del2,1,dla,dlai,nptdif == 4)
        call dfphiz(de1,del1,de2,del2,1,phia,phiai,nptdif == 4)
        do  12  i = 1, nptdif
          dlap(i)  = dlai(i)
          phiap(i) = phiai(i)
   12   continue
      endif
C     Energy derivatives of g
      call dfphiz(de1,del1,de2,del2,2*nr,g,gp,nptdif == 4)
C     p = integral <gp gp>
      call gintz(gp,gp,a,b,nr,z,e,l,v,rofi,p)
      endif

C     phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
      phi = val(1)/rmax
      dphi = (slo(1) - phi)/rmax
      if (nptdif /= 0) then
        phip = vali(1)/rmax
        dphip = (sloi(1) - phip)/rmax
      endif

C     Copy energy derivatives vali,sloi to val(2..),slo(2..)
      if (nptdif /= 0) then
        call zcopy(nptdif,sloi,1,slo(2),1)
        call zcopy(nptdif,vali,1,val(2),1)
      endif

      if (iprint() >= 111) write(stdo,749) phi,dphi,phip,dphip
  749 format(' PHIDZ:  phi,phip,phip,dphip=',8f12.6)

      end
