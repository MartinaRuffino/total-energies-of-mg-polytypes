      subroutine gippad(mode,modi,pgplp,glist,npl,nz,nzp,zp,wz,nl,nsp,
     .  nspc,indxsh,lidim,pp,vshft,ipc,nrclas,nrhos,nthet,dlmcl,
     .  nzg,nzpg,zpg,gd,ef0,vne,qne,vshfo,vshfn,dos,qnu,rhos)
C- Estimate some integrated quantities from Pade interpolation
C-----------------------------------------------------------------------
Ci Inputs
Ci   mode:  Controls what integration is performed (see Remarks)
Ci          1s digit
Ci          0  compute difference in integrated quantities
Ci             int dz (I(z-vshfn) - I(z-vshfo))
Ci          1  compute difference in integrated quantities
Ci             int dz (I(z+vshfn) - I(z+vshfo))
Ci          2  Generate integrated quantities without potential shifts
Ci         10s digit
Ci          0 dv is a constant;  1  dv depends on PL
Ci        100s digit
Ci          1 exclude PL -1, npl from shifts
Ci       1000s digit
Ci          1 extrapolate dos(ef) to real axis, estimated ef
Ci            (dos must be set)
Ci   modi:  1s digit
Ci          1 compute dos (see Outputs)
Ci        100s digit
Ci          1 compute emom (see Outputs)
Ci       1000s digit
Ci          4 compute spin density matrix rhos
Ci   pgplp :index and dimensioning information for each PL
Ci          1: number of basis atoms
Ci          2: index labelling kind of potential of this PL
Ci          3: source (column) dimension of GF for this PL
Ci          4: row (field) matrix dimension for this PL
Ci          5: matrix dimension for this PL including i-waves
Ci          6: offset to diagonal part of g
Ci   glist :list of PL over which to generate array
Ci   npl   :number of principal layers (pgfset.f)
Ci   nz    :number of energy points including points on
Ci          non-equilibrium contour (nz>nzp in non-equilibrium mode)
Ci   nzp   :number of energy points on equilibrium contour
Ci   zp    :points for complex energy contour.  g is interpolated
Ci          from the set of zpg points to these points, and the
Ci          selected quantities are integrated at these points
Ci   wz    :weights for complex energy integration
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   pp    :potential parameters, needed for the enu when calculating
Ci          the energy moments of the charge (atomsr.f)
Ci   vshft :array of site potential shifts, needed only to make qnu.
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   nth   :nth(i) = number of angles at the ith basis site
Ci         :nth(i) should be 1 except in the disordered local moments case
Ci   dlmcl :class index pointers for extra DLM classes
Ci   gd    :diagonal GF
Ci   nzg   :number of points for which diagonal GF has been calculated,
Ci          including points on non-equlibrium contour
Ci   nzpg  :number of equilibrium points for which diagonal GF has been
Ci          calculated (in non-equiibrium mode, when nz>nzp, it is
Ci          required that nz=nzg and nzp=nzpg)
Ci   zpg   :energy points for which diagonal GF has been calculated
Ci   ef0:   estimate for Fermi level is ef0+vshfn
Ci          (used only to estimate DOS(ef0+vshfn)
Ci   vne   :(non-equilibrium mode) difference in fermi energies of
Ci           right and left leads vne=ef(R)-ef(L); may be <0
Co   qne   :addition to charge by integration along nonequilib. contour
Ci   vshfo :initial potential; see Remarks
Ci   vshfn :final potential for which selected quantities are integrated
Cio Inputs/Outputs
Cio  dos   :PL-integrated dos and moments of dos
Cio         Input is dos at vshfo; output is dos at vshfn
Cio         Input dos is dos at vshfo
Cio          1: dos at this e : -1/pi sum_k Im G
Cio          2: nos at this e : -1/pi sum_k Im (wz*G)
Cio          3: 1st energy mom : -1/pi sum_k Im (wz*G (z-enu))
Cio          4: 2nd energy mom : -1/pi sum_k Im (wz*G (z-enu)**2)
Cio          5: projection of dos to fermi level (real axis)
Cio  qnu   : moments of DOS
Cio          Input is qnu at vshfo; output is qnu at vshfn
Cio  rhos   :spin density matrix
Cio          Input is rhos at vshfo; output is rhos at vshfn
Cl Local variables
Cl   offgi0:offset to first orbital of site ib1, in indxsh table
Cl   offgd :offset to gd array for current site
Cl   offgd0:offset to gd array for current PL
Cl   lmRL  :entry in indxsh for current orbital
Cl   lmgd  :current index to gd array relative to start of PL
Cr Remarks
Cr   gippad estimates a change in selected integrated properties from a
Cr   constant potential shift.  It uses a rational function (Pade)
Cr   approximation to interpolate the diagonal gii calculated at the
Cr   points {zpg} to points {zp}, and integrates the selected
Cr   properties using the interpolated values.
Cr   ... Case 1s digit mode=0,1
Cr     Integrated properties are ADDED at points {zp}+vshfn,
Cr     and SUBTRACTED at points {zp}-vshfo, thus estimating a
Cr     correction to the original integral owing to potential shifts
Cr     vshfn-vshfo
Cr   ... Case 1s digit mode=2
Cr     Integrated properties are simply evalated at points {zp};
Cr     vshfn and vshfo are not used.
Cu Updates
Cu   20 Dec 12 Pade extrap of Im G to Ef capped so DOS>0 for each channel
Cu   01 Nov 11 (Belashchenko) Parameters nthet,dlmcl used for DLM
Cu   27 Jan 04 (S.Faleev) non-equilibrium mode.  Changed argument list.
Cu   26 Feb 02 gippad works for a hamiltonian constructed from a
Cu             permuted list of sites.  Bug fix, 100s digit mode
Cu             Altered argument list.
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nz,nsp,nspc,npl,ipc(*),nrclas(*),nl,mode,modi,nzg,lidim
      integer glist(-1:*),pgplp(6,-1:npl),indxsh(*),nzp,nzpg,nrhos
      double precision zp(2,*),wz(2,*),gd(2,nzg,nsp,nspc,0:*),ef0,
     .  vshfo(2),vshfn(2),dos(5,nsp,*),qnu(3,nl,nsp,*),
     .  pp(6,nl,nsp,*),vshft(*),rhos(2,3,nrhos,2,2,*),
     .  zpg(2,nzg),vne,qne(nsp,*)
      integer nthet(*),dlmcl(*)
C ... Dynamically allocated arrays
      complex(8), allocatable :: cof(:,:),cof2(:,:)
C ... Local parameters
      logical ldosp,lmom,skippl
      integer i,ib,ib1,ib2,ic,ie,ipl,is,isp,jsp,kpl,kplv,l,
     .  ldmat,lmRL,lmgd,m,offgd,offgd0,offgi0,scrwid,ilm,k2
      integer ic0,ic1,ic2
      double precision bync,gpad(2),z(2),pi,zmenu(2),zf(2)
      double precision xx(2,0:2),fac,wt,wtz(2),vshfi
      parameter (pi=3.141592653589793d0,scrwid=80)
C ... Parameters for non-equilibrium mode (nz>nzp)
      integer modne,n1,n2,j,modlin
      double precision vshsm,vshmx,vshmd,f(2),g1,g2,gsave
      parameter (vshsm=0.03d0,vshmx=0.2d0,vshmd=0.1d0,modlin=1)
      procedure(integer) :: iprint,isw,lgunit

      integer mpipid,procid,master
c     logical mlog

      procid = mpipid(1)
      master = 0

C      call fopna('gd',-1,0)
C      call zprm('gd',2,gd(1,1,1,1,0),nzg,nzg,4)

      ldosp = mod(modi,10) == 1 .or. mod(modi,10) == 2
      lmom  = mod(modi/100,10) == 1
      ldmat = mod(modi/1000,10)
C      if (ldmat >= 4 .and. procid == master)
C     .  print *,'WARNING(gippad): sdmat with DLM'
      if (iprint() >= 20) then
        call awrit4(' gippad: integrate %i point%-1j%?#n>1#s## '//
     .    '%?#n<2#corr. from dv shift#by Pade interpolation#:'//
     .    '%?#n# idos,##'//
     .    '%?#n# qnu,##'//
     .    '%b',
     .    ' ',scrwid,lgunit(1),nz,mod(mode,10),isw(ldosp),isw(lmom))
        if (iprint() >= 40 .and. mod(mod(mode/10,10),2) == 1) then
        call awrit2(' vold =%n;10,6D',' ',scrwid,lgunit(1),glist,vshfo)
        call awrit2(' vnew =%n;10,6D',' ',scrwid,lgunit(1),glist,vshfn)
        endif
      endif

C     Sanity checks
      if ( nz > nzp .and. (nz /= nzg .or. nzp /= nzpg) )
     .  call rx('gippad: non-equilibrium mode and nz!=nzg or nzp!=nzpg')
      if ( nz > nzp .and. abs(ef0-zp(1,nzp+1)) > 1d-10)
     .  call rx('gippad: wrong ef(L) for non-equilibrium contour')
      if ( nz > nzp+1 .and. abs(ef0+vne-zp(1,nz)) > 1d-10)
     .  call rx('gippad: wrong ef(R) for non-equilibrium contour')
      if ( nz == nzp .and. nzg /= nzpg)
     .  call rx('gippad: equilibrium mode (nz=nzp) and nzg!=nzpg')

C ... Allocate memory for Pade coefficients
      allocate(cof(nzpg,nzpg+2))
      if (nspc == 2) allocate(cof2(nzpg,nzpg+2))
      do  isp = 1, nsp
C     jsp = 1 in collinear case, isp in the noncollinear; thus
C     (isp,jsp) are indices to the diagonal in the spinor part of gd
      jsp = min(isp,nspc)

C --- For each PL, do ---
      offgd0 = 0
      do  kpl = 1, glist(-1)
      ipl = glist(kpl-1)
C     Determine now whether to skip this PL.
C     But defer skipping until properly run over loop counters
C     to keep track of, e.g., offgd
      skippl = mod(mode/100,10) == 1 .and. (ipl < 0 .or. ipl == npl)
      call gtibpl(ipl,npl,pgplp,ib1,ib2)
      kplv = 1
      if (mod(mode/10,10) == 1) kplv = kpl
      if (ldosp .and. .not. skippl) then
        dos(1,isp,kpl) = 0
        dos(5,isp,kpl) = 0
      endif

C --- For each site in ib1..ib2, do ---
      offgd = offgd0
      offgi0 = nl**2*(ib1-1)
      lmRL = offgi0
C     call offsHp(indxsh,ib1,ib2,0,0,lidim,i,ldimi)
C     ldimi = ldimi+i
      do  ib = ib1, ib2
        ic0 = ipc(ib)
        bync = 1d0 / nrclas(ic0)
        lmgd = -1
        ilm = 0
        do  l = 0, nl-1
        do  m = -l, l

          ilm = ilm+1
          k2 = l+1
          if (nrhos == nl*nl) k2 = ilm

C         Skip orbitals outside size of g lower+interm. block
          lmRL = lmRL + 1
          if (indxsh(lmRL) > lidim) cycle
          lmgd = lmgd + 1
          if (skippl) cycle

C         call zprm('g',2,gd(1,1,isp,jsp,offgd+lmgd),nzg,nzg,1)
          call padcof(nzpg,zpg,gd(1,1,isp,jsp,offgd+lmgd),nzpg,cof)
          if (nspc == 2)
     .    call padcof(nzpg,zpg,gd(1,1,isp,3-isp,offgd+lmgd),nzpg,cof2)

C     ... Loop twice, adding at vshfn, subtracting at vshfo
          fac = -1
          if (mod(mode,10) == 1) fac = 1
          do  is = 1, 2

          fac = -fac
          if (mod(mode,10) /= 2) then
            if (is == 1) vshfi =  fac*vshfn(kplv)
            if (is == 2) vshfi = -fac*vshfo(kplv)
          else
            vshfi = 0
          endif

C     --- For each energy, do ---
          do  ie = 1, nzp

C     ... Pade approximation to g for z-vshfi
          z(1) = zp(1,ie)
          z(2) = zp(2,ie)
          wt = -1/pi * fac
          if (nsp == 1) wt = wt*2
          wtz(1) = wz(1,ie)*wt
          wtz(2) = wz(2,ie)*wt
          z(1) = z(1) - vshfi
          call pade(1,z,nzpg,zpg,nzpg,cof,gpad)
C     ... xx(0) = fac * (-1/pi w(z) * g)
          xx(1,0) = wtz(1)*gpad(1) - wtz(2)*gpad(2)
          xx(2,0) = wtz(1)*gpad(2) + wtz(2)*gpad(1)

C     ... Make PL-dos, add to PL-integrated dos and moments
          if (ldosp) then
            if (ie == nzp .and. is == 1) then
              dos(1,isp,kpl) = dos(1,isp,kpl) + wt * gpad(2)

C             Pade approximation to g for ef-vshfi
              zf(1) = ef0 - vshfi
              zf(2) = 0
              gsave = gpad(2)
              call pade(1,zf,nzpg,zpg,nzpg,cof,gpad)
              if (gpad(2) > 0) gpad(2) = gsave
              dos(5,isp,kpl) = dos(5,isp,kpl) + wt*gpad(2)
C             print *, 'ie,isp',ie,isp,sngl(wt),sngl(dos(5,isp,kpl))
            endif
            xx(1,1) = z(1)*xx(1,0) - z(2)*xx(2,0)
            xx(2,1) = z(2)*xx(1,0) + z(1)*xx(2,0)
            xx(1,2) = z(1)*xx(1,1) - z(2)*xx(2,1)
            xx(2,2) = z(2)*xx(1,1) + z(1)*xx(2,1)

            do  i = 0, 2
            dos(2+i,isp,kpl) = dos(2+i,isp,kpl) + xx(2,i)
C       ... Non-equilibrium branch
            if (ie == nzp .and. nz-nzp >= 2) then
C             CONSTANT interpolation of DOS using G^< for vshfi>vshsm
              modne = 0
C             LINEAR interpolation of DOS using G^< for vshfi<vshsm
C             DO NOT use linear interpolation until all checks passed
              if (abs(vshfi) <= vshsm) modne = modlin
              do j = 1,2
                n1 = (nzp+1)*(2-j) + nz*(j-1)
                n2 = (nzp+2)*(2-j) + (nz-1)*(j-1)
                g1 = -1/pi*gd(2,n1,isp,jsp,offgd+lmgd)
                g2 = -1/pi*gd(2,n2,isp,jsp,offgd+lmgd)
                if (abs(vshfi) > vshmx) then
C                 For large shifts neglect the addition to integrated
C                 dos from shift by non-equilibrium contour:
C                 interpolation is not accurate in that case
                  f(j) = 0d0
                else
                  call pgnei(modne,zp(1,n1),zp(1,n2),g1,g2,
     .              vshfi,0d0,i,f(j))
                  if (abs(vshfi) > vshmd)
     .              f(j)=f(j)*(vshmx-abs(vshfi))/(vshmx-vshmd)
                endif
              enddo
C             Addition to integrated dos from shift of non-equil. contour
              dos(2+i,isp,kpl) = dos(2+i,isp,kpl) + (f(1)-f(2))*fac
              if (i == 0)
     .          qne(isp,kpl) = qne(isp,kpl) + (f(1)-f(2))*fac
            endif
C       ... End of non-equilibrium branch
            enddo
          endif
C         if (is == 1 .and. ie == nz)
C    .    print 333,ie,is,ib,kpl,dos(1,1,kpl)
C 333     format(4i4,f15.10)

C     ... Make emom ... pot shift also shifts enu
          if (lmom) then
            if (nthet(ib) > 1) then
              ic1 = dlmcl(ib)
              ic2 = ic1 + nthet(ib) - 1
            else
              ic1 = ic0
              ic2 = ic0
            endif
            do  ic = ic1, ic2
            zmenu(1) = z(1) - (pp(1,l+1,isp,ic)+vshft(ib))
            zmenu(2) = z(2)
            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
            do  i = 0, 2
            qnu(i+1,l+1,isp,ic) = qnu(i+1,l+1,isp,ic) + xx(2,i)*bync
            if (ldmat >= 4) then
              rhos(1,i+1,k2,isp,isp,ic) = rhos(1,i+1,k2,isp,isp,ic) + xx(2,i)*bync
            endif
C       ... Non-equilibrium branch
            if (ie == nzp .and. nz-nzp >= 2) then
C             CONSTANT interpolation of DOS using G^< for vshfi>vshsm
              modne = 0
C             LINEAR interpolation of DOS using G^< for vshfi<vshsm
C             DO NOT use linear interpolation until all checks passed
              if (abs(vshfi) <= vshsm) modne = modlin
              do j = 1,2
                n1 = (nzp+1)*(2-j) + nz*(j-1)
                n2 = (nzp+2)*(2-j) + (nz-1)*(j-1)
                g1 = -1/pi*gd(2,n1,isp,jsp,offgd+lmgd)
                g2 = -1/pi*gd(2,n2,isp,jsp,offgd+lmgd)
                if (abs(vshfi) > vshmx) then
C                 For large shifts neglect the addition to integrated
C                 dos from shift by non-equilibrium contour:
C                 interpolation is not accurate in that case
                  f(j) = 0d0
                else
                  call pgnei(modne,zp(1,n1),zp(1,n2),g1,g2,
     .              vshfi,pp(1,l+1,isp,ic)+vshft(ib),i,f(j))
                  if (abs(vshfi) > vshmd)
     .              f(j) = f(j)*(vshmx-abs(vshfi))/(vshmx-vshmd)
                endif
              enddo
C             Addition to moms from shift of non-equil. contour
              qnu(i+1,l+1,isp,ic) = qnu(i+1,l+1,isp,ic) + (f(1)-f(2))*fac*bync
            endif
C       ... End of non-equilibrium branch
            enddo
C       ... End of loop over ic added to cover DLM sites
          enddo
          endif

C     ... Spin density matrix
          if (ldmat >= 4) then
c           if (nthet(ib) > 1)
c    .        call rx('gippad not ready for sdmat with DLM')
            if (nthet(ib) > 1) then
              ic1 = dlmcl(ib)
              ic2 = ic1 + nthet(ib) - 1
            else
              ic1 = ic0
              ic2 = ic0
            endif
            do  ic = ic1, ic2
            call pade(1,z,nzg,zpg,nzg,cof2,gpad)
            xx(1,0) = wtz(1)*gpad(1) - wtz(2)*gpad(2)
            xx(2,0) = wtz(1)*gpad(2) + wtz(2)*gpad(1)
            xx(1,1) = zmenu(1)*xx(1,0) - zmenu(2)*xx(2,0)
            xx(2,1) = zmenu(2)*xx(1,0) + zmenu(1)*xx(2,0)
            xx(1,2) = zmenu(1)*xx(1,1) - zmenu(2)*xx(2,1)
            xx(2,2) = zmenu(2)*xx(1,1) + zmenu(1)*xx(2,1)
            do  i = 1, 3
C             Im rhos = Re g;  Re rhos = -Im g
              rhos(1,i,k2,isp,3-isp,ic) = rhos(1,i,k2,isp,3-isp,ic) + xx(2,i-1)*bync
              rhos(2,i,k2,isp,3-isp,ic) = rhos(2,i,k2,isp,3-isp,ic) - xx(1,i-1)*bync
            enddo
          enddo
          endif

C       End of energy loop
        enddo
        if (mod(mode,10) == 2) exit
C       End of + vshfn - vshfo loop
        enddo ! loop over spin
        enddo ! loop over m
      enddo   ! loop over l
      offgd = offgd+lmgd+1
      enddo ! loop over site
      offgd0 = offgd
      enddo
      enddo

C     call awrit1('%10:1,d',' ',120,6,dos)
c     call yprm('dosi',1,dos,0,5*nsp,5*nsp,7)

      end
      subroutine pgevdc(mode,npl,nsp,dos,eold,enew)
C- Shift sum-of-eigenvalues from Pade approximant
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :0 eold,enew a constant, independent of PL
Ci         :1 eold,enew depend on PL
Ci   npl   :number of principal layers
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   de    :potential shift of Green'
Co Inputs/Outputs
Cio  dos   :PL-integrated dos and moments of dos
Cio          1: dos at this e : -1/pi sum_k Im G
Cio          2: nos at this e : -1/pi sum_k Im (wz*G)
Cio          3: 1st energy mom : -1/pi sum_k Im (wz*G (z-enu))
Cio          4: 2nd energy mom : -1/pi sum_k Im (wz*G (z-enu)**2)
Cio          5: projection of dos to fermi level (real axis)
Cio        :On output, 1st energy moment is shifted
Cr Remarks
Cr
Cu Updates
Cu   22 Mar 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsp,npl,mode
      double precision dos(5,nsp,npl),eold(npl),enew(npl)
C ... Local parameters
      integer kpl,isp
      double precision de

      do  kpl = 1, npl
        if (mode == 0) then
          de = enew(1)-eold(1)
        else
          de = enew(kpl)-eold(kpl)
        endif
        do  isp = 1, nsp
          dos(3,isp,kpl) = dos(3,isp,kpl) + dos(2,isp,kpl)*de
C         print *, kpl,dos(2,isp,kpl),dos(3,isp,kpl)
        enddo
      enddo

      end

      subroutine gippadr(nl,ic1,ic2,qnu,qnusav,qnur)
C- Convert Pade estimate for change in qnu to change in qnur
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   ic1
Ci   ic2
Ci   qnu   :Pade corrected energy-weighted moments of the sphere charges (lms)
Ci   qnusav:qnu before Pade correction
Cio Inputs/Outputs
Cio  qnur  :On input, Dirac moments before Pade correction
Cio  qnur  :On output, Pade estimate added to Dirac moments
Cr Remarks
Cr
Cu Updates
Cu   09 Apr 15 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ic1,ic2
      integer, parameter :: nsp=2
      double precision qnu(3,nl,nsp,ic2),qnusav(3,nl,nsp,ic2),qnur(4,nl,2*nl,2,2,ic2)
C ... Local parameters
      integer ic
      double precision qnuloc(3,nl,2),qnurloc(4,0:nl-1,2*nl,2,2)

      do  ic = ic1, ic2
        call dcopy(3*nl*2,qnu(1,1,1,ic),1,qnuloc,1)
        call daxpy(3*nl*2,-1d0,qnusav(1,1,1,ic),1,qnuloc,1)

        call pshpr(10)
        call dpzero(qnurloc,size(qnurloc))
        call qnu2qnur(2,nl,nsp,qnuloc,qnuloc,qnurloc) ! Map to lambda-mu, zeroing out off-digonal parts Q2
        call poppr
        call daxpy(4*8*nl*nl,1d0,qnurloc,1,qnur(1,1,1,1,1,ic),1) ! Add to qnur

      enddo

      end
