      subroutine smshft(s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot)
C- Estimate the smooth density for a shift in atomic positions.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos pos0
Co     Stored:     pos
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvsms1 rhgcmp rhogkl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z p pz lmxa name lmxl a nr rmt nxi exi chfa rsmfa
Ci                 coreh coreq rg lfoca rfoca qc ctail etail stc lmxb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvsms1 gtpcor rhgcmp corprm rhogkl
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nabc ng vol nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:kv gv ips0 bgv
Cio    Passed to:  pvsms1 rhgcmp symsmr
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfrce
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  elind
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  smrho
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rhat
Cio    Passed to:  *
Ci Inputs
Co Outputs
Co   smrho :a perturbation is added to smrho, depending on job
Cr Remarks
Cr   job describes which ansatz for charge shift is used for correction
Cr     <=0  do not calculate correction to force
Cr       1  shift in free-atom density
Cr       2  shift in core+nuclear density
Cr     +10  to screen the rigid shift by the Lindhard function
Cr   (job taken from ctrl->lfrce)
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   17 Sep 01 Adapted for local orbitals.  Altered argument list
Cu    3 Jul 00 Adapted from nfp smshft.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      complex(8), allocatable :: cgr(:)
      complex(8), allocatable :: cgs(:)
      complex(8), allocatable :: cwk(:)
      complex(8),pointer :: smrho(:,:)
C ... Local parameters
      integer i,ib,iprint,is,k1,k2,k3,kcor,lcor,lmxa,n0,n1,
     .  n2,n3,nbas,ng,ngabc(3),nglob,nsp
      integer kmax,job
      parameter (n0=10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,elind,pi,qc,qsc,qcor(2),plat(3,3),
     .  qlat(3,3),qv,qval,tpiba,vol,z,pnu(n0,2),pnz(n0,2)

      job = s_ctrl%lfrce
      if (job <= 0) return
      call tcn('smshft')

C --- Setup and printout ---
      nsp  = nglob('nsp')
      nbas = nglob('nbas')
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol
      call fftz30(n1,n2,n3,k1,k2,k3)
C     call zprm3('input smrho',0,smrho,k1,k2,k3)
      smrho => s_pot%smrho

C ... Hold on to original smrho (cgr)
C     call defcc(ocgr,  ng*nsp)
      allocate(cgr(ng*nsp))
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng,nsp,s_lat%kv,k1,k2,k3,smrho,cgr)

C --- Shift in unscreened density at the two positions ---
C      call defcc(ocgs, ng*nsp)
C      call defcc(ocwk, ng)
      allocate(cgs(ng*nsp),cwk(ng))
      kmax = 0
      call pvsms1(s_site,s_spec,s_lat,nbas,nsp,kmax,
     .  ng,s_lat%gv,s_pot%rhat,cwk,cgs,job)
      deallocate(cwk)

C ... Debugging: print unscreened shift in pseudo core density
C      call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cgs,smrho)
C      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('unscreened local density',0,smrho,k1,k2,k3)

C --- Screened shift ---
      if (job > 10) then
C       Compute elind if not given
        qval = 0d0
        do  ib = 1, nbas
          is = s_site(ib)%spec
          z = s_spec(is)%z
          pnu = s_spec(is)%p
          pnz = s_spec(is)%pz
          lmxa = s_spec(is)%lmxa
          if (lmxa == -1) cycle
          call gtpcor(s_spec,is,kcor,lcor,qcor)
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
          qval = qval+qv
        enddo
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
        elind = s_ham%elind
        if (elind < 0d0) elind = -(3*pi**2*qval/vol)**.66666d0*elind
        if (nsp == 2) call dsumdf(ng*2,1d0,cgs,0,1,cgs,ng*2,1)
C        call lindxx(122,n1,n2,n3,k1,k2,k3,ng,s_lat%kv,cgs,s_lat%gv,
C     .    tpiba,elind,w,w,w,w)
        call lindsc(2,ng,s_lat%gv,tpiba,elind,cgs)
C        call dscal(2*ng,.001d0,cgs,1)

        if (nsp == 2) call dsumdf(ng*2,.5d0,cgs,0,1,cgs,ng*2,1)
      endif

C ... Debugging: show delta smrho
      call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cgs,smrho)
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('screened delta smrho',0,smrho,k1,k2,k3)
C      print *, 'returning with delta smrho'
C      return

C --- Add shift to smrho, ensuring no shift in <rho> ---
      do  i = 1, nsp
        call dvset(cgs,1+ng*(i-1),1+ng*(i-1),0d0)
      enddo
      call dpadd(cgr,cgs,1,ng*2*nsp,1d0)
      call gvputf(ng,nsp,s_lat%kv,k1,k2,k3,cgr,smrho)
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)

C --- Symmetrize the shifted density ---
      call symsmr(s_lat,nsp,0,k1,k2,k3,smrho)

      if (iprint() > 100)
     .  call zprm3('shifted smrho',0,smrho,k1,k2,k3*nsp)

      deallocate(cgr,cgs)
      call tcx('smshft')

      end

      subroutine pvsms1(s_site,s_spec,s_lat,nbas,nsp,kmax,ng,gv,s_rhat,
     .  cwk,cg,job)
C- Shift smoothed density according to job
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos pos0 rho1 rho2 rhoc rho1x rho2x rhocx
Co     Stored:    pos
Cio    Passed to: rhgcmp rhogkl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa
Ci                coreh coreq rg lfoca rfoca qc ctail etail stc lmxb
Co     Stored:    *
Cio    Passed to: gtpcor rhgcmp corprm rhogkl
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat vol plat qlat nabc ng gv
Co     Stored:    *
Cio    Passed to: rhgcmp
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   kmax  :polynomial cutoff for augmentation expansion
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cwk   :work array of same dimension as cg
Ci   job :describes which ansatz for charge shift is used for correction
Ci         <=0  do not calculate correction to force
Ci           1  shift in free-atom density
Ci           2  shift in core+nuclear density
Ci         +10  to screen the rigid shift by the Lindhard function
Co Outputs
Co  cg   coefficients to FT of shifted density
Cl Local variables
Cl  qloc   :difference in true and smoothed local charge.
Cl         :If the charge is assembled from overlapping atom-centered
Cl         :densities, qloc is the difference between the smoothed
Cl         :and head densities.
Cr Remarks
Cr   Shift of the "free atom densities."
Cr   The table below shows the densities and corresponding charges, and
Cr   parameters that hold their representations (true and smoothed
Cr   approximate forms):
Cr      density     charge    reps'n     smooth -reps'n
Cr      rho(smH)     qfat    cofh,ceh    already smooth
Cr      rho1-rho2    qloc     rhoat      qkl
Cr      rhoc         qc
Cr   This routine constructs the following difference
Cr   rhat(final) - rhat(initial)  positions, in the smooth reps'n, where
Cr      rhat = rho(smH) + qg * g(r)
Cr   where
Cr       qg = qval+qsc-qfat-qloc
Cr
Cr   In the special case rho is assembled from a superposition of
Cr   free-atom densities, and rho1-rho2 = rhoval(free-atm)-rho(smH)
Cr   (see ovlcor.f).  Thus in this case:
Cr      rho(free-atom) = rho(smH) + rho1-rho2 + rhoc
Cr   with the corresponding integrated charges
Cr       qval=z-qc     = qfat     + qloc      - qc
Cr   Thus in this special case qg=0: the only shift comes from rho(smH).
Cr   Because the local density (which contains the remaining part of
Cr   the free-atom density) will automatically be shifted, it follows
Cr   that the shifted smooth density will correspond to the
Cr   smooth sum-of-FA densities constructed at the shifted positions.
Cr
Cr   In the general case, qg is not zero.  By shifting the a gaussian
Cr   along with the sm-Hankels, the integrated total density of charge
Cr   shifted (local density + mesh density) is neutral.
Cr
Cr   Improvements: if the tail density were also shifted inside each
Cr   augmentation sphere, the total density would correspond exactly
Cr   to the sum-of-FA densities at the shifted positions, when the
Cr   starting density is also a sum-of-FA densities.
Cr
Cr   Shift of the "core + valence densities."
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,ng,job
      double precision gv(ng,3)
      double complex cg(ng,nsp),cwk(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat) ::  s_lat
      type(str_rhat)::  s_rhat(*)
C ... Local parameters
      integer ib,i,is,iv0,kmax,lmxl,lmxa,nr,nrmx,n0,nlml,nxi,ie,ixi,ig,
     .  ipr,iprint,kcor,lcor,stdo,lgunit,intopt,nglob
      parameter (nrmx=5001, n0=10)
      double precision a,aa,alat,df(0:20),e,
     .  exi(n0),gam,hfc(n0,2),pi,pnew(3),pnu(n0,2),pnz(n0,2),pold(3),
     .  pp,qall,qc,qcor(2),qsc,qfat,qg,qloc,qval,rmt,
     .  rsmfa,rwgt(nrmx),scalp,sum,tpiba,v(3),v2,vol,volsp,y0,z
      character*35 strn,spid*8
      double complex phase

C ... Setup
      call tcn('pvsms1')
      call stdfac(20,df)
      alat = s_lat%alat
      vol = s_lat%vol
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      tpiba = 2*pi/alat
      ipr = iprint()
      stdo = lgunit(1)
      volsp = vol*nsp
      intopt = 10*nglob('lrquad')

C --- For each site, accumulate shift in density ---
      call dpzero(cg,2*ng*nsp)
      iv0 = 0
      strn = 'free atom densities'
      if (job == 11) strn = 'screened free atom densities'
      if (job == 12) strn = 'screened core+multipole densities'
      if (ipr >= 30) write (stdo,1) strn
    1 format(/' smshft:  add shifted ',a/'   site',16x,'old pos',22x,
     .  'new pos',14x,'shift')

      do  ib = 1, nbas
        is = s_site(ib)%spec
        spid = s_spec(is)%name
        lmxl = s_spec(is)%lmxl
        nlml = (lmxl+1)**2
        if (lmxl == -1) cycle

        is = s_site(ib)%spec
        pnew = s_site(ib)%pos
        pold = s_site(ib)%pos0
        pp = alat*dsqrt((pnew(1)-pold(1))**2 + (pnew(2)-pold(2))**2
     .                + (pnew(3)-pold(3))**2)
        if (ipr >= 30) write(stdo,2) ib,spid,pold,pnew,pp/alat
    2   format(i4,':',a,f8.5,2f9.5,2x,3f9.5,2x,f9.6)

C       Skip this site if shift is negligible
        if (pp > 1d-6) then

C   --- Shift in mesh density, job 1 ---
        if (mod(job,10) == 1) then
          z = s_spec(is)%z
          pnu = s_spec(is)%p
          pnz = s_spec(is)%pz
          lmxa = s_spec(is)%lmxa
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          nxi = s_spec(is)%nxi
          exi = s_spec(is)%exi
          hfc = s_spec(is)%chfa
          rsmfa = s_spec(is)%rsmfa
          call gtpcor(s_spec,is,kcor,lcor,qcor)
          if (nr > nrmx) call rx('dfrce: nr gt nrmx')
          call radwgt(intopt,rmt,a,nr,rwgt)
          call radsum(nr,nr,nlml,nsp,rwgt,s_rhat(ib)%rho1,qloc)
          call radsum(nr,nr,nlml,nsp,rwgt,s_rhat(ib)%rho2,sum)
          qloc = (qloc-sum)/y0
          qfat = 0d0
          do  i = 1, nsp
            do  ie = 1, nxi
              gam  = 0.25d0*rsmfa**2
              qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
              qfat = qfat + hfc(ie,i)*qall
            enddo
          enddo
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
C         Excess sphere charge.  See Remarks above.
          qg = qval+qsc-qfat-qloc

C     ... Shift in smoothed free atom density
          do  i = 1, nsp
          do  ixi = 1, nxi
          e = exi(ixi)
          do  ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + hfc(ixi,i)*aa*phase*y0/vol
          enddo
          enddo
          enddo

C     ... Add gaussian to conserve local charge; see Remarks
          do  i = 1, nsp
          do  ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + qg*phase*dexp(-gam*v2)/volsp
          enddo
          enddo

C   --- Shift in mesh density, job 12 ---
        elseif (job == 12) then

C     ... Core + valence at old position
          call dpzero(cwk,ng*2)
          s_site(ib)%pos = pold
          call rhgcmp(131,ib,ib,s_site,s_spec,s_lat,s_rhat,kmax,ng,cwk)
          call dscal(ng*2,-1d0,cwk,1)
C     ... Core + valence at new position
          s_site(ib)%pos = pnew
          call rhgcmp(131,ib,ib,s_site,s_spec,s_lat,s_rhat,kmax,ng,cwk)

C     ... Add to cg
          do  i = 1, nsp
            call daxpy(ng*2,1d0/nsp,cwk,1,cg(1,i),1)
          enddo

        else
          call rxi('smshft: bad job:',job)
        endif
        endif

        iv0 = iv0+nlml
      enddo

      call tcx('pvsms1')
      end

      subroutine pvsms2(s_site,s_spec,rotm,nbas,nsp)
C- Rotate local densities by specified rotation
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Co   orhoat:On output the different m-channels of rhoat(1) and rhoat(2)
Co         :are mixed by the rotation
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name nr lmxl
Ci     Stored:    *
Ci     Passed to: *
Ci   rotm  :3x3 cartesian rotation matrix
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Co Outputs
Cl Local variables
Cr Remarks
Cr   For a rotation matrix R, The density is stored in the 1-center form
Cr      rho_l(r) YL(rhat)
Cr   Given a rotation matrix R, this it transforms as
Cr      rho_l(r) YL(R rhat) = rho_l(r) rYL(rhat)
Cr   where rYL is made by ylmrtg
Cr
Cb Bugs
Cb   No ability is supplied when the Yl are true instead of real
Cb   spherical harmonics
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   21 Dec 04 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp
      double precision rotm(3,3)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,i,j,is,lmxl,nr,nlml,ipr,stdo,nlx,nl2,nglob
      parameter (nlx=9, nl2=nlx*nlx)
      double precision rYL(nl2,nl2),det33,det
      character spid*8

C ... Setup
      call getpr(ipr)
      stdo = nglob('stdo')

C ... Rotation matrix for real spherical harmonics
C     call prmx('pvsms2 rotm',rotm,3,3,3)
      call ylmrtg(nl2,rotm,rYL)
C     call prmx('rYL',rYL,nl2,nl2,nl2)

C --- For each site and l, rotate the m-components ---
      if (ipr >= 20) then
        call info0(20,0,0,' Rotate local densities using R=')
        write (stdo,1) ((rotm(i,j),j = 1,3),i = 1,3)
    1   format(3F11.6)
      endif
      det = det33(rotm)
      if (abs(abs(det)-1) > 1d-6)
     .  call info2(10,0,0,
     .  ' (warning) rotation is not unitary: determinant = %d',det,0)

      do  ib = 1, nbas
        is = s_site(ib)%spec
        spid = s_spec(is)%name
        nr = s_spec(is)%nr
        lmxl = s_spec(is)%lmxl
        if (lmxl == -1) cycle
        nlml = (lmxl+1)**2
        if (nlml > nl2) call rx('increase nl2 in pvsms2')

        call pvsms3(nr,nr,nlml,nsp,rYL,nl2,s_site(ib)%rho1)
        call pvsms3(nr,nr,nlml,nsp,rYL,nl2,s_site(ib)%rho2)

C       print "(i4,100f15.10)", ib,s_site(ib)%rho1(nr,:)

      enddo

      end

      subroutine pvsms3(nrx,nr,nlml,nsp,rYL,nl2,rho)
C- Rotation of an l-dependent density
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrx   :leading dimension of rho
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cl   rYL   :rotation matrix that rotates Y_lm
Ci   nl2   :leading dimension of rYL
Co Outputs
Co   rho   :On output the different m-channels of rho are
Co         :mixed by rYL
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   21 Dec 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrx,nr,nlml,nsp,nl2
      double precision rho(nrx,nlml,nsp)
      double precision rYL(nl2,nl2)
C ... Local parameters
      integer isp
C     integer l,lmax,ll,nlmi,offri
      double precision rwk(nrx,nlml)

C     lmax = ll(nlml)

C     call prmx('starting rho',rho,nrx,nr,nlml*nsp)

      if (nlml == 0) return
      do  isp = 1, nsp

        call dgemm('N','T',nr,nlml,nlml,1d0,rho(1,1,isp),nrx,
     .    rYL,nl2,0d0,rwk,nrx)
        call dcopy(nrx*nlml,rwk,1,rho(1,1,isp),1)

C        faster if done l-by-l
C        do  l = 0, lmax
C
C          nlmi = 2*l + 1
C          offri = l**2
C          print *, l, nlmi,offri
C
C          call dgemm
C        enddo
      enddo

C     call prmx('ending rho',rho,nrx,nr,nlml*nsp)

      end
