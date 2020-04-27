      subroutine emesh(semsh,z,w)
C- Energy contour for numerical integration in complex plane
C ----------------------------------------------------------------
Ci   semsh :array containing information about energy contour
Ci         :1   nzp   number of energy points on equilibrium contour
Ci         :2   mode  specifies contour; see zmesh.
Ci         :          Basic contours spec'd by 1s and 10s digit
Ci         :          0   Real part of energy z on a uniform mesh
Ci                        of nz points between eb and ef
Ci         :          1   Same as 0, but Im z opposite sign
Ci         :          2   Same as 0, but all weights w = 1
Ci         :          10  elliptical contour
Ci         :          11  elliptical contour, but Im z opposite sign
Ci         :          20  elliptical contour as in 10
Ci         :          30  vertical contour starting on the real axis
Ci         :          Higher digits for special modifications
Ci         :          100's digit=1,2 => points for nonequilibrium
Ci         :                             layer Green's function
Ci         :          100's digit=3,4 => mixed elliptical+real axis
Ci                                       to find fermi level
Ci         :          Add 100 for nonequil part using Im(z)=delne
Ci         :          Add 200 for nonequil part using Im(z)=del00
Ci         :          Add 300 for mixed elliptical contour +
Ci                            real axis to find fermi level
Ci         :          Add 1000 to set nonequil part only.
Ci         :3   emin  lower bound for energy contour
Ci         :4   emax  upper bound for energy contour
Ci         :5-6 ecc,eps  (1s digit modec>=10, elliptical mesh)
Ci                    eccentricity for ellipse (0 => circle, 1 => line)
Ci                    and bunching parameter for points near Ef
Ci         :5   delta (applies when Im z=const, i.e. when modec=0,1)
Ci                    Imaginary part of z
Ci         :5   alfa  (applies when when mode=30)
Ci                    exponent for Laguerre quadrature on vertical contour
Ci         :... The following apply to the layer Green's function,
Ci         :    nonequilibrium case (100s digit modec = 1,2)
Ci         :7   nzne  number of points on nonequil. contour
Ci         :8   vne   difference in fermi energies of right and
Ci                    left leads vne=ef(R)-ef(L); may be <0
Ci         :9   delne plays the role of delta for energy points on
Ci                    on nonequilibrium contour z=E+i*delne
Ci                    in evaluation of G^<
Ci         :10  del00 substitutes for delne when making the surface
Ci                    self-energy.
Ci         :... The following apply to 100s digit modec = 3,4
Ci         :7   delz  mesh spacing searching for charge. neutrality
Ci         :... The following applies to the layer Green's function,
Ci         :    when 1+10s digit mode = 0,1
Ci         :10  del00 Im z just for surface self-energy (i.e. delta
Ci                    specifically for calculating surface self-energy)
Co Outputs
Co   z     :(complex) energy
Co   w     :(complex) weights
Cu Updates
Cu   04 Mar 05 100s digit mode=3 or 4
Cu   10 Jan 04 (S.Faleev) Added 100 to semsh(2) for non-equilibrium mode
Cu   26 Feb 01 Added 1s digit modec = 2
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      double precision semsh(10)
      double precision z(2,*),w(2,*)
C ... Local parameters
      integer nzp,modec,modne,nzp0,nzne,iprint,stdo,nglob,neonly
      double precision emin,emax,ecc,eps,delta,vne,delne,del00,delz,alfa

C     call tcn('emesh')
      stdo = nglob('stdo')
      nzp   = semsh(1)
      nzp0  = nzp
      modec = mod(nint(semsh(2)),100)
      modne = mod(nint(semsh(2))/100,10)
      neonly= mod(nint(semsh(2))/1000,10)
      emin  = semsh(3)
      emax  = semsh(4)
      delta = semsh(5)
      ecc   = semsh(5)
      alfa  = semsh(5)
      eps   = semsh(6)
      del00 = semsh(10)
      nzne  = 0

C     modec = 300,400
      delz = 0
      if (modne == 3 .or. modne == 4) then
        delz = semsh(7)
        modne = 0
      endif

C     Non-equil. mode with nzp=0 used only for DOS from G^< (lmpg=1)
C     and a calculation using uniform mesh
      if (modne > 0 .and. nzp0 <= 0) then
        if (modec /= 0) call
     .    rx('emesh: non-equilib. mode with nzp=0 requires modec=0')
        nzp   = semsh(7)
        delta = semsh(9)
      endif

      if (neonly == 0) then
        if (iprint() >= 20) write(stdo,*) ' '
        call zmesh(modec,emin,emax,ecc,eps,delta,delz,nzp,z,w)
        if ((modec == 0 .or. modec == 1) .and. del00 /= 0) then
        call info2(20,0,0,'%16pUse Im(z)=%1;4g for self-energy',del00,0)
        endif
      endif

C     Non-equil. mode with nzp>0 used for lmpg=5 (nzne=0), and lmpg=1,
C     where nzne>0 used for self-consistent calculations and nzne=0
C     used for DOS by Im[G^r]
      if (modne /= 0 .and. nzp0 > 0) then
        call isanrg(modne,0,2,'emesh:','100s digit mode',.true.)
        nzne  = nint(semsh(7))
        if (nzne > 0) then
        vne   = semsh(8)
        delne = semsh(9)
        del00 = semsh(10)
        if (modne == 1) delta = delne
        if (modne == 2) delta = del00

        call info5(20,0,0,'%8fNonequil:  nzne=%i  Im(z)=%g (%g,SE)'
     .    //'  ef=%,5g  bias=%;4g',nzne,delne,del00,emax,vne)
C       call zmshne(emax,vne,delta,nzp,nzne,z,w)
C       if (nzne < 2) call rx('emesh: illegal mesh: nzne<2')

        call pshpr(iprint()-20)
        call zmesh(mod(modec,10),emax,emax+vne,ecc,eps,delta,delz,nzne,
     .    z(1,1+nzp),w(1,1+nzp))
        call poppr
        endif
      endif

C      do  modne = 1, nzp+nzne
C        print 333, z(1,modne),z(2,modne),w(1,modne),w(2,modne)
C  333   format(2f18.12,2x,2f18.12)
C      enddo
C      call zprm('z',2,z,nzp+nzne,nzp+nzne,1)
C      call zprm('w',2,w,nzp+nzne,nzp+nzne,1)
C      stop

C     call tcx('emesh')
      end
      subroutine zmesh(modec,eb,ef,ecc,eps,delta,delz,nz,z,w)
C- Makes a complex energy mesh for contour integration
C ----------------------------------------------------------------
Ci modec:  0 Real part of energy z: uniform mesh of nz points
Ci                                  between eb and ef
Ci           Imaginary part of z: constant shift i*delta
Ci           Weights w set for trapezoidal rule, i.e.
Ci           w = energy increment = (ef-eb)/(nz-1)
Ci         1 same as mode 0, except Im(z) = -i*delta
Ci         2 same path as mode 0, except that all weights w = 1
Ci        10 elliptical contour, top half of complex plane
Ci        11 elliptical contour, bottom half of complex plane
Ci        20 integration exactly as in 10 (this mode will not
Ci           calculate energy-integrated quantities, see gfasa.f)
Ci delta: imaginary part for uniform mesh (constant Im z)
Ci        (or alpha for modec=30, vertical contour)
Ci eb,ef: integration from eb to ef
Ci ecc:   eccentricity for ellipse (0 => circle, 1 => line)
Ci eps:   parameter controlling logarithmic bunching of points near Ef
Ci        eps=0 => no bunching
Ci nz:    number of points
Co Outputs
Co   z,w: complex points and weights for contour integration
Cr Remarks
Cr   Contour integration by Gaussian quadrature of an angle theta on
Cr   an elliptical path, for theta in [-pi,0].  When eps>0, there is
Cr   an additional mapping to a logarithmic mesh to bunch points near
Cr   the real axis.  Choose eps in [0,1], the closer to 1 the more
Cr   bunching:  theta(x)=(1-exp(x))/eps, x(theta)=log(1-theta*eps)
Cr   Points are ReZ = z0 + zr*cos(theta), Im Z = zr*(1-ecc)*sin(theta),
Cr   where z0=(ef+eb)/2 and zr=(ef-eb)/2.
Cr
Cr   For uniform mesh, weights are defined for trapezoidal integration
Cu Updates
Cu   26 Feb 01 Added 1s digit modec = 2
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nz,modec
      double precision z(2,nz),w(2,nz),wk(nz,2),delz
C     double complex z(nz),w(nz)
      double precision eb,ef,ecc,eps,delta
C ... Local parameters
      logical lellip,lwt,lsignz,ltmp,lvert
      integer i,nglob,iprint,isw,stdo,modloc
      double precision wgt,expx,theta,ct,st,pi,z0,zr,rm,de,del,eta,datan
      double precision alfa
      character*80 outs

      modloc = modec
      if (modloc == 20) modloc=10

      call isanrg(mod(modloc,10),0,2,'zmesh:','1s digit modec',.true.)
      lsignz = mod(modloc,10) == 1
      lellip = mod(modloc/10,10) == 1
      lvert = mod(modloc/10,10) == 3
      lwt = .false.
      alfa = ecc
      if (.not. lellip .and. .not. lvert) then
        lwt  = mod(modloc,10) == 2
      endif
      stdo = nglob('stdo')

C --- Setup and default settings ---
      pi = 4*datan(1d0)
      eta = 1-ecc
      if (ecc <= 0d0) eta = 1d0
      del = delta
*     if (delta == 0d0) del = 1d-2
      i = isw(delz /= 0)
      call awrit8('%x ZMESH: nz=%i (ellipse%?#n#+line##)  ecc,eps='//
     .  '%1;4,1d,%1;4d  e1,e2=%1;4,2d,%1;4,2d  %?#n#delz=%1,1;3g##',
     .  outs,80,0,nz,i,1-eta,eps,eb,ef,i,delz)
      if (.not. lellip .and. .not. lvert) then
        if (lwt) then
          de = 1
        else
          if (nz /= 1) then
            de = (ef-eb)/(nz-1)
          else
            de = ef-eb
          endif
        endif
        call awrit6('%x ZMESH: nz=%i (real axis: Im(z)=%1;4g'//
     .    '%?;n==1;(*-1);;)  e1,e2=%1;4g,%1;4g  wz=%1;4g',
     .  outs,80,0,nz,del,isw(lsignz),eb,ef,de)
      endif
      ltmp = nz <= 0 .or. nz <= 1 .and. (lellip .or. lvert .or. lwt)
C     ltmp = nz <= 0 .or. nz <= 1 .and. .not. (lellip .or. lwt)
      if (ltmp .or. iprint() >= 20) call awrit0('%a',outs,-80,-stdo)
      call rxx(ltmp,'... illegal mesh')

C --- Contour case ---
      if (lellip) then
C   ... Gaussian quadrature on interval [0,log(1-eps)] or [0,1] if eps=0
C   ... Use z, and w as a work space
c       call mklegw(nz,z,z(1+nz,1),0)
        call mklegw(nz,wk,wk(1,2),0)
C   ... Flip to  put z of highest last
        call dscal(nz,-1d0,wk,1)
C   ... Copy (pts,wgt) to (real, imag) of array w for temp storage
        do  5  i = 1, nz
        w(1,i) = wk(i,1)
    5   w(2,i) = wk(i,2)
        rm = 1
        if (eps /= 0) rm = dlog(1d0-eps)
        do  10  i = 1, nz
        w(1,i) = rm/2*w(1,i) + rm/2
   10   w(2,i) = w(2,i)*rm/2
C   ... Additional mapping to to logarithmic mesh in [0,1]
        if (eps /= 0) then
          do  20  i = 1, nz
          expx = dexp(w(1,i))
          w(1,i) = (1d0-expx)/eps
   20     w(2,i) = -expx*w(2,i)/eps
        endif
C   ... Map to points on an ellipse
        zr = (ef-eb)/2
        z0 = (ef+eb)/2
        do  30  i = 1, nz
          theta = -pi + pi*w(1,i)
          wgt = w(2,i)
          ct = dcos(theta)
          st = dsin(theta)
          z(1,i) = z0 + zr*ct
          z(2,i) = -zr*eta*st
          w(1,i) = -pi*zr*st*wgt
          w(2,i) = -pi*zr*eta*ct*wgt
          if (lsignz) then
            z(2,i) = -z(2,i)
            w(2,i) = -w(2,i)
          endif
   30   continue
      elseif (lvert) then
C   ... Laguerre quadrature on a vertical contour starting on the real axis---
        call mklagw(nz,1,alfa,0,wk,wk(1,2),0)
        do  i = 1, nz
          z(1,i) = ef
          z(2,i) = wk(i,1)
          w(1,i) = 0
          w(2,i) = - wk(i,2)
        enddo
C --- Uniformly spaced points shifted by delta off the real axis---
      else
        if (lsignz) del = -del
        i = 3
        if (lwt) i = 2
        call pshtzi(i,0,nz,eb,ef,del,z,w)

C        de = (ef-eb)/max(nz-1,1)
C        do  40  i = 1, nz
C          z(1,i) = eb + (i-1)*de
C          z(2,i) = del
C          if (lsignz) z(2,i) = -z(2,i)
C          if (lwt) then
C            w(1,i) = 1
C          else
C            w(1,i) = de
C          endif
C          w(2,i) = 0d0
C   40   continue
C        if (.not. lwt) then
C          w(1,1) = w(1,1)/2
C          w(1,nz) = w(1,nz)/2
C        endif
      endif

C --- Printout ---
      if (iprint() >= 50 .and. lellip .or. iprint() >= 80) then
        write(stdo,332)
  332   format(10x,'z',21x,'w')
        do  50  i = 1, nz
   50   write(stdo,333) z(1,i),z(2,i),w(1,i),w(2,i)
  333   format(2f10.6,2x,2f10.6)
      endif
      end
      subroutine shftzi(mode,semsh,izp,z,w)
C- Assign Im z for some energy points, depending on mode
C-----------------------------------------------------------------
Ci Inputs
Ci   mode  :1 contour for active layers
Ci         :2 contour for surface GF
Ci   semsh :array containing information about energy contour (cf emesh)
Ci   izp   :if  0, assign Im z for all appropriate points on contour
Ci         :if >0, assign Im z only for point izp, if this point is
Ci                appropriate (see Remarks)
Co  Outputs
Co   zp,w  :points, weights for nonequilibrium portion of energy contour
Co         :are adjusted for whether surface or active layer
Cr Remarks
Cr   This routine shifts delta = Im z depending on mode.
Cr   For the layer GF case, we want to have a separate Im z
Cr   for surface GF than for the active region, for stability.
Cr
Cr   Case 1: uniform mesh on real axis.
Cr           Active region,  Im z = delne = semsh(5)
Cr           Surface region, Im z = del00 = semsh(10)
Cr
Cr   Case 2: nonequilibrium mesh, which at present,
Cr   the nonequilibrium contour is on real axis.
Cr   Im(z) must be very small (like 1d-10) in order that GF<
Cr   equals imaginary part of retarded GF in equilibrium.
Cr           Active region,  Im z = delne = semsh(9)
Cr           Surface region, Im z = del00 = semsh(10)
Cr
Cr   Note: if del00, no shifts are made
Cr
Cr   By assigning a larger Im(z) for the surface GF, it
Cr   actually defines the imaginary part of GF in all layers
Cr   through the self-energy contribution to the GF
Cr   The imaginary part for the surface GF is del00 (see emesh.f)
Cr
Cr   For this contour, w is not affected
Cr
Cu Updates
Cu   23 Apr 05 Adapted for modec=0,1
Cu   15 Sep 04 (MvS) Adapted from original shftdel
Cu    8 Jun 04 (S.Faleev) First created
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,izp
      double precision semsh(10)
      double precision z(2,*),w(2,*)
C ... Local parameters
      integer nzne,nzp,modec,modne,nz1,nz2
      double precision delne,del00,del,xx

      nzp   = semsh(1)
      nzne  = nint(semsh(7))
      delne = semsh(9)
      del00 = semsh(10)
      if (del00 == 0) return
      modec = mod(nint(semsh(2)),100)
      modne = mod(nint(semsh(2))/100,10)

C     Lower bound =0 if main contour is uniform mesh
      if (modec == 0 .or. modec == 1) then
        nz1 = 0
      else
        nz1 = nzp
      endif

C     Upper bound depends on whether nonequil case, or not
      if (modne == 1) then
        nz2 = nzp+nzne
      else
        nz2 = nzp
      endif

      if (nz1 >= nz2) return

C     delta for central region
      if (mode == 1) then
        if (modec == 0) then ! uniform mesh with fixed Im z
          del = semsh(5)
        elseif (modec == 1) then
          del = -semsh(5)
        elseif (modne == 1) then ! Nonequil part of contour has fixed Im z
          del = semsh(9)
        else
          return
        endif
      elseif (mode == 2) then
        del = semsh(10)
        if (modec == 1) del = -del
      else
        call isanrg(mode,1,2,'shftzi:','mode',.true.)
      endif

C      call info5(0,0,0,
C     .  ' shftzi:  shifting  Im z  to %g, pts %i..%i',
C     .  del,nz1+1,nz2,0,0)

C     Shift just 1 point if izp>0
      if (izp > 0) then
        if (izp <= nz1 .or. izp > nz2) return
        nz1 = izp-1
        nz2 = izp
      endif

C     Do the shift
      xx = 0  ! To avoid run time check faults
      call pshtzi(0,nz1,nz2,xx,xx,del,z,[xx])

      end

      subroutine pshtzi(mode,nz1,nz2,e1,e2,zim,z,wz)
C- Mesh of uniformly spaced points shifted by zim off the real axis
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 assign imaginary part of z only
Ci         :1 assign z (both real and imaginary parts)
Ci         :2 assign z and weight wz=1
Ci         :3 assign z and weight wz according to trapezoidal rule
Ci   nz1   :index to point before first energy point
Ci   nz2   :index to last  energy point
Ci   e1    :first energy point (used only for mode>0)
Ci   e2    :last  energy point (used only for mode>0)
Ci   zim   :imaginary part of z
Co Outputs
Co   z     :complex energy mesh : z(nz1+1..nz2) assigned
Ci   wz    :(mode>0) weights for energy integration : w(nz1+1..nz2) assigned
Cr Remarks
Cr
Cu Updates
Cu   23 Apr 05 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nz1,nz2
      double precision e1,e2,z(2,nz2),wz(2,nz2),zim
C ... Local parameters
      integer iz
      double precision de

C     Spacing on real axis
      if (mode > 0) then
        de = (e2-e1)/max(nz2-nz1-1,1)
      endif

C     Assign uniform mesh
      do  iz = nz1+1, nz2
        z(2,iz) = zim
        if (mode > 0) then
          z(1,iz) = e1 + (iz-1)*de
          if (mode == 2) then
            wz(1,iz) = 1
            wz(2,iz) = 0d0
          elseif (mode == 3) then
            wz(1,iz) = de
            wz(2,iz) = 0d0
          endif
        endif
      enddo

      if (mode == 3 .and. nz1+1 < nz2) then ! Trapezoidal rule
        wz(1,nz1+1) = wz(1,nz1+1)/2
        wz(1,nz2)   = wz(1,nz2)/2
      endif

      end

C      subroutine zmshne(ef,vne,delne,nzp,nzne,z,w)
CC-Makes a uniform energy mesh for G^< calculation in non-equlibrium mode
CC-----------------------------------------------------------------------
CCi ef,vne: mesh from min(ef,ef+vne) to max(ef,ef+vne)
CCi delne : small shift to top half-plane of complex energy: z=E+i*delne
CCi nzp   : number of energy points on equilibrium contour
CCi nzne  : number of energy points on non-equilibrium contour
CCi z,w complex countors and weights for points 1:nzp
CCo Outputs
CCo   z   :(complex) energies for points 1:nzp+nzne
CCo   w   :(complex) weights  for points 1:nzp+nzne
CCu Updates
CCu 10 Jan 04 (S.Faleev) First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nzp,nzne
C      double precision z(2,nzp+nzne),w(2,nzp+nzne)
CC ... Local parameters
C      integer iprint,i,stdo,nglob
C      double precision ef,vne,delne,de
C
C      call info5(20,0,0,'%8fNonequil:  nzne=%i  Im(z)=%g'//
C     .  '  ef=%,4d  bias=%,4d',nzne,delne,ef,vne,0)
C
C      if (nzne <= 1) call rx('zmshne: illigal mesh nzne<2')
C      if (abs(vne) < 1d-5)
C     .      call rx('zmshne: too small vne<1d-5 while nzne>1')
C
C      de = vne/(nzne-1)
C      do  i = nzp+1, nzp+nzne
C        z(1,i) = ef+de*(i-nzp-1)
C        z(2,i) = delne
CC       w(1,i) is negative for vne<0, as it should be
C        w(1,i) = de
C        w(2,i) = 0d0
C      enddo
C      w(1,nzp+1)  = de/2d0
C      w(1,nzp+nzne) = de/2d0
C
C      if (iprint() >= 80) then
C        stdo = nglob('stdo')
C        do  i = nzp+1, nzp+nzne
C          write(stdo,333) z(1,i),z(2,i),w(1,i),w(2,i)
C  333     format(2f10.6,2x,2f10.6)
C        enddo
C      endif
C
C      end
C      subroutine xshftzi(mode,delne,del00,nzp,nzne,zp)
CC- Shift imaginary part of energy on non-equilibrium contour
CC-----------------------------------------------------------------
CCi Inputs
CCi   mode  :1 shift by del00
CCi         :2 shift by delne
CCi   nzp   :number of energy points on equilibrium contour
CCi   nzne  :number of points on nonequil. contour
CCio Inputs/Outputs
CCi   zp    :points for complex energy contour: imaginary part
CCi         :is shifted to del00 (mode=1) or delne (mode=2)
CCi         :for points nzp+1..nzp+nzne
CCr Remarks
CCr   delne should be small (like 1d-10) in order that GF<
CCr   equals imaginary part of retarded GF in equilibrium.
CCr   del00 actually defines the imaginary part of GF in all layers
CCr   through the self-energy term
CCu Updates
CCu    8 Jun 04 (S.Faleev) First created
CC ----------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nzp,nzne,mode
C      double precision zp(2,nzp+nzne), del00,delne
CC ... Local parameters
C      integer iz
C
C      if (nzne < 1) call rx('shftzi: nzne should be > 0')
C
CC     print *, '!! shftzi'
CC      open(unit=999,file='delne')
CC      read(999,*)delne,del00
CC      close(999)
C
C      do  iz = nzp+1, nzp+nzne
C        if (iz > nzp) then
C          if (mode == 1) then
CC           zp(2,iz) was delne, becomes del00
C            if (abs(zp(2,iz)-delne) > 1d-10)
C     .        call rx('shftzi: wrong delne')
C            zp(2,iz) = del00
C          elseif (mode == 2) then
CC           zp(2,iz) was del00, becomes delne
C            if (abs(zp(2,iz)-del00) > 1d-10)
C     .        call rx('shftzi: wrong del00')
C            zp(2,iz) = delne
C          else
C            call rx('shftzi: wrong mode')
C          endif
C        endif
C      enddo
C      end
C      subroutine fmain
CC Tests zmesh and Pade approximations to GF
CC There are three test cases: simple pole, simple s-band,
CC and semi-infinite s-band.  Tests for integrated properties
CC are written to stdout; also tests for DOS, which are calculated
CC directly and estimated by a Pade approximant.  These
CC are written to file dos.dat.
CC
CC To check stdout and dos with reference
CC echo / | a.out | tee out
CC diff out out.pade; mc -f'(7f12.6/12x,6f12.6)' dos.dat dos.pade -- -px| head
CC
CC To compare the exact DOS with the Pade approximant for these cases:
CC fplot -frme 0,1,1,1.5 -lt 2 -colsy 2 dos.dat -colsy 8 -lt 3,bold=6 dos.dat\
CC       -frme 0,1,.5,1 -lt 2 -colsy 4 dos.dat -colsy 10 -lt 3,bold=6 dos.dat\
CC       -frme 0,1,0,.5 -lt 2 -colsy 6 dos.dat -colsy 12 -lt 3,bold=6 dos.dat
C      implicit none
C      integer nzmx,nz,npad
C      parameter(nzmx=200)
C      double complex zm(nzmx),wz(nzmx),g(nzmx),sum,z3,pdot,dcsqrt
C      double complex zpad(nzmx),wpad(nzmx),gp(nzmx),gb(nzmx),gs(nzmx)
C      double precision eb,ef,ecc,eps,delta,result,pi,fac,xx
C      double precision dosb(nzmx),doss(nzmx),gfrb(nzmx),gfrs(nzmx),
C     .  dosp(nzmx),gfrp(nzmx),cpad(nzmx,nzmx+2)
C      double precision pole,psmall
C      integer i,ifi,jfi,fopng,sgn
C      pi = 4*datan(1d0)
C      eb = -2.1d0
C      ef = 0d0
C      nz = 20
C      ecc = .5d0
C      eps = 0d0
C      pole = -.1d0
C      delta =  .001d0
C      psmall = .02d0
C      sgn = -1
C      print *, 'input nz,ecc,eps,pole,psmall,delta'
C      read(*,*) nz,ecc,eps,pole,psmall,delta
C      npad = 200
C      i = (1+sgn)/2
C      call zmesh(i,1.5d0*eb,-1.5d0*eb,0d0,0d0,delta,npad,zpad,wpad)
C      call pshpr(80)
C      call zmesh(10+i,eb,ef,ecc,eps,delta,nz,zm,wz)
C      ifi = -1
C      ifi = fopng('zpoints.dat',ifi,0)
C      do  5  i = 1, nz
C      write(ifi,100) dble(zm(i)), dimag(zm(i))
C    5 continue
C  100 format(2e16.7)
C      call awrit6(' nz=%i  ecc=%d  eps=%d  pole=%d  psmall=%d'//
C     .  '  delta=%d',' ',80,6,nz,ecc,eps,pole,psmall,delta)
C
Cc --- test 1: simple pole 1/(z-a) with pole a inside integration interval
Cc     (Test is strictest when pole close to fermi energy)
Cc     pole = -0.01d0
C      sum = dcmplx(0d0,0d0)
C      xx = 0
C      do  10  i = 1, nz
C        g(i) = 1d0/(zm(i)-pole)
C        xx = max(xx,cdabs(g(i)))
C        sum = sum + g(i)*wz(i)
C   10 continue
C      result = sgn*dimag(sum)/pi
C      call awrit3(' Simple pole at e=%1;4d: gmax='//
C     .  '%1;3g  sum = %1;10d',' ',80,6,pole,xx,result)
C      call padcof(nz,zm,g,nz,cpad)
C      call pade(npad,zpad,nz,zm,nz,cpad,gp)
C
Cc --- test 2: contour integration of full-chain s-band gf x 2 for spin
Cc     integrated up to half band should give 1
C      sum = dcmplx(0d0,0d0)
C      xx = 0d0
C      do  20  i = 1, nz
C        fac = -dsign(1d0,dble(zm(i)))
C        z3 = zm(i) + zm(i)**3*psmall
C        pdot = 1 + 3*zm(i)*zm(i)*psmall
C        g(i) = -2d0*fac*pdot/dcsqrt(z3**2-4d0)
C        xx = max(xx,cdabs(g(i)))
C        sum = sum + g(i)*wz(i)
C   20 continue
C      result = sgn*dimag(sum)/pi
C      call awrit2(' s band chain: gmax='//
C     .  '%1;3g  sum = %1;10d',' ',80,6,xx,result)
C      call padcof(nz,zm,g,nz,cpad)
C      call pade(npad,zpad,nz,zm,nz,cpad,gb)
C
Cc --- test 3: integration of semi-infinite s-band gf x 2 for spin
C      sum = dcmplx(0d0,0d0)
C      xx = 0d0
C      do  30  i = 1, nz
C        fac = -dsign(1d0,dble(zm(i)))
C        z3 = zm(i) + zm(i)**3*psmall
C        pdot = 1 + 3*zm(i)*zm(i)*psmall
C        g(i) = pdot*(z3 + fac*dcsqrt(z3**2 - 4))
C        xx = max(xx,cdabs(g(i)))
C        sum = sum + g(i)*wz(i)
C   30 continue
C      result = sgn*dimag(sum)/pi
C      call awrit2(' semi-infinite s-band: gmax='//
C     .  '%1;3g  sum = %1;10d',' ',80,6,xx,result)
C      call padcof(nz,zm,g,nz,cpad)
C      call pade(npad,zpad,nz,zm,nz,cpad,gs)
C
CC --- DOS on a fine mesh ---
C      nz = npad
C      do  50  i = 1, nz
C        g(i) = 1d0/(zpad(i)-pole)
C        dosp(i) = sgn*dimag(g(i))/pi
C        gfrp(i) = dble(g(i))/pi
C        z3 = zpad(i) + zpad(i)**3*psmall
C        pdot = 1 + 3*zpad(i)*zpad(i)*psmall
C        fac = -dsign(1d0,dble(zpad(i)))
C        g(i) = -2*fac*pdot/dcsqrt(z3**2 - 4)
C        dosb(i) = sgn*dimag(g(i))/pi
C        gfrb(i) = dble(g(i))/pi
C        g(i) = pdot*(z3 + fac*dcsqrt(z3**2 - 4))
C        doss(i) = sgn*dimag(g(i))/pi
C        gfrs(i) = dble(g(i))/pi
C   50 continue
C      print *, 'dumping dos in file dos.dat ...'
C      jfi = -1
C      jfi = fopng('dos.dat',jfi,0)
C      rewind jfi
C      call awrit2('%% rows %i cols %i',' ',80,jfi,nz,7+6)
C      write(jfi,301)
C  301 format('#',5x,'z',9x,'dos(pole)  Re(G)',8x,
C     .  'dos(sband) Re(G)',8x,'dos(semi-inft) Re(G)/pade')
C      do  60  i = 1, nz
C   60 write(jfi,300) dble(zpad(i)),
C     .    dosp(i),gfrp(i),dosb(i),gfrb(i),doss(i),gfrs(i),
C     .    sgn*dimag(gp(i))/pi,dble(gp(i))/pi,
C     .    sgn*dimag(gb(i))/pi,dble(gb(i))/pi,
C     .    sgn*dimag(gs(i))/pi,dble(gs(i))/pi
C      call fclose(ifi)
C      call fclose(jfi)
C  300 format(7f12.6/12x,6f12.6)
C      call cexit(1,1)
C      end
