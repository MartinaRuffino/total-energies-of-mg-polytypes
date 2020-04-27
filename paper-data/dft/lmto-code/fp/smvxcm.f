      subroutine smvxcm(s_site,s_spec,s_lat,nbas,nsp,bscal,lfrce,k1,k2,k3,
     .  smrho,smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,focexc,focex,focec,focvxc,f)
C- XC potential for smooth mesh density
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  smcorm smvxc4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lfoca rfoca qc z ctail etail stc lmxb p pz rmt rg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  smcorm corprm smvxc4
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc ng vol alat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv kv cy
Cio    Passed to:  smcorm smvxc2 vxcnlm
Ci Inputs
Ci   nbas  :size of basis
Ci   bscal :scale vxc+ - vxc- by bscal (nsp=2 only)
Ci   lfrce :1, calculate contribution to forces
Ci   k1,k2,k3 dimensions smooth crystal densities, potentials on a mesh
Ci   smrho :smooth valence density on uniform mesh
Co Outputs
Co   smvxc :ex-corr  potential of smoothed density + core corrections
Co   smvx  :exchange potential of smoothed density + core corrections
Co   smvc  :correlation potential of smoothed density + core corrections
Co   smexc :ex-corr  energy density of smoothed density + core corrections
Co   smpot :smooth total potential; smvxc is added to smpot
Co   repsm :integrated exchange-correlation energy
Co   repsmx:integrated exchange energy
Co   repsmc:integrated correlation energy
Co   rmusm :int (smrho + smcor1) * vxc[rhosm+smcor1]
Co         :where smcor1 = portion of core treated directly
Co   rvmusm:int (smrho) * vxc[rhosm+smcor1]
Co   rvepsm:int (smrho) * exc[rhosm+smcor1]
Co   focexc:FOCA exchange-correlation energy:
Co         :int (smcor2) * vxc[rhosm+smcor1]
Co         :where smcor2 = portion of core treated perturbatively
Co   focex :exchange part of focexc
Co   focec :correlation part of focexc
Co   focvxc:integral of FOCA exchange-correlation potential:
Co         :int (smcor2) * (smrho) * (dvxc/drho)
Co   f     :contribution to forces from xc potential
Cl Local variables
Cl  lxcfun :switch defining xc functional; see lxcf in structures.h
Cr Remarks
Cr   smoothed core is partition into core1 + core2.  All atoms with
Cr   lfoc1=1 belong to core1; those with lfoc2=1 belong to core2.
Cr  *core1 is included directly into smrho; the nonlinear XC potential
Cr   is computed from vxc[smrho+smcor1].
Cr  *core2 is included perturbatively: its contribution to the vxc
Cr   is computed from the expansion
Cr     vxc[rho + smcor2] = vxc[rho] + smcor2 * (dvxc/drho)
Cr                       = vxc[rho] + dvxc
Cr   The perturbation correction to int (smrho * vxc) is then
Cr     focvxc = int smrho * smcor2 * (dvxc/drho)
Cr   If the perturbation approach is exact,
Cr     (focvxc+rvmusm) -> rvmusm when computed with smcor2=0
Cr   The corresponding XC energy density is
Cr     exc[rho + smcor2] = exc[rho] + smcor2 * (dexc/drho)
Cr                       = exc[rho] + smcor2 * (vxc-exc)/rho
Cr   The perturbation correction to the XC energy is then
Cr     int smcor2 * (vxc-exc) = focexc - int smcor2 exc[smrho+smcor1]
Cr
Cr   If vxc+ - vxc- is scaled by bscal, terms exc and integrate quantities
Cr   such as rhoeps, rhomu, lose their meaning. The potential is merely
Cr   scaled; no attempt is made to construct a functional corresponding
Cr   to the scaling.
Cu Updates
Cu   31 Jul 14 XC B-field can be scaled by bscal
Cu   09 Dec 13 First cut at using libxc functional for XC potential
Cu   03 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   21 Apr 09 Handles GGA functionals
Cu   02 Jul 05 skip sites for which cofh=0
Cu   25 Jun 04 return smexc,rvepsm
Cu   14 Jun 02 rhoex and rhoec (T. Miyake)
Cu    8 Feb 02 smvx and smvc (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp sxc_smooth.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,lfrce,k1,k2,k3
      double precision f(3,nbas),repsm(2),repsmx(2),repsmc(2),rmusm(2),
     .  rvmusm(2),rvepsm(2),focexc(2),focex(2),focec(2),focvxc(2),bscal
      double complex smrho(k1,k2,k3,2),smpot(k1,k2,k3,2),
     .  smvxc(k1,k2,k3,2),smvx(k1,k2,k3,2),smvc(k1,k2,k3,2),
     .  smexc(k1,k2,k3)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      complex(8), allocatable,target :: smrhol(:)
      complex(8), allocatable :: dxcv(:),cgh1(:),cgh2(:)
      complex(8), pointer :: smcor(:)
C ... Local parameters
      integer i,iprint,k123,lfoc1,lfoc2,lgunit,lxcfun,n1,
     .  n2,n3,ngabc(3),ng,nglob,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C     integer ocgh1,ocgh2,odxcv,osmcor,osmrho
      double precision vol,sum1,sum2,vxcavg(2),x1,x2,alat
      character*180 outs

      stdo = lgunit(1)
      lxcfun = nglob('lxcf')
      call tcn('smvxc')
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol
      alat = s_lat%alat
      vol = s_lat%vol

C ... Sum of foca hankel heads; break into direct and pert. parts
      allocate(smrhol(k1*k2*k3*nsp)); call dpzero(smrhol,2*k1*k2*k3*nsp)
      allocate(cgh1(ng),cgh2(ng))
      call smcorm(nbas,s_site,s_spec,s_lat,ng,s_lat%gv,cgh1,cgh2,
     .  lfoc1,lfoc2)
      if (lfoc2 /= 0 .and. mod(lxcfun/100,100) /= 0) call
     .  rx('smvxcm: perturbation treatment not implemented'//
     .  ' with this functional')

C ... smrhol = smrho + smoothed core from foca hankel heads
      k123 = 2*k1*k2*k3
      if (lfoc1 == 1) then
        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cgh1,smrhol)
        call fftz3(smrhol,n1,n2,n3,k1,k2,k3,1,0,1)
        if (nsp == 2) then
          call dscal(k123,.5d0,smrhol,1)
          call dpscop(smrhol,smrhol,k123,1,1+k123,1d0)
        endif
        call dpadd(smrhol,smrho,1,k123*nsp,1d0)
      else
        call dpcopy(smrho,smrhol,1,k123*nsp,1d0)
      endif

C ... Force density strictly positive definite
C      print *, 'density strictly pos def?'
C      call smrpos(smrhol,k1,k2,k3,n1,n2,n3)

      rvmusm(2) = 0
      focexc(2) = 0; focex(2) = 0; focec(2) = 0; focvxc(2) = 0

C --- Direct branch (lfoc2 == 0) ---
      if (lfoc2 == 0) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrhol,sum1,sum2)
        call smvxc2(0,s_lat,nsp,lxcfun,bscal,vol,n1,n2,n3,k1,k2,k3,
     .    smrhol,smvxc,smvx,smvc,smexc,x1,
     .    repsm,repsmx,repsmc,rmusm,vxcavg)
        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),
     .      smrho(1,1,1,i),rvmusm(i),x2)
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),
     .      smrho(1,1,1,i),rvepsm(i),x2)
        enddo
        focexc(1) = 0d0
        focex(1)  = 0d0
        focec(1)  = 0d0
        focvxc(1) = 0d0
      endif

C --- Perturbation branch (lfoc2 /= 0) ---
      if (lfoc2 /= 0) then
        allocate(dxcv(k1*k2*k3*nsp))
        call smvxc2(1,s_lat,nsp,lxcfun,bscal,vol,n1,n2,n3,k1,k2,k3,
     .   smrhol,smvxc,smvx,smvc,smexc,dxcv,
     .   repsm,repsmx,repsmc,rmusm,vxcavg)
        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),
     .      smrho(1,1,1,i),rvmusm(i),x2)
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),
     .      smrho(1,1,1,i),rvepsm(i),x2)
        enddo

C   ... Assemble core tails for linearized treatment on mesh
C       w(osmcor) = portion of core treated perturbatively
        smcor => smrhol
        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cgh2,smcor)
        call fftz3(smcor,n1,n2,n3,k1,k2,k3,1,0,1)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smcor,sum1,sum2)
        call dpzero(focexc,2)
        call dpzero(focex,2)
        call dpzero(focec,2)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvx(1,1,1,i),smcor,x1,x2)
          focex(i)  = x1/nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvc(1,1,1,i),smcor,x1,x2)
          focec(i)  = x1/nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),smcor,x1,
     .      x2)
          focexc(i) = x1/nsp
C         Add this term to focexc to make focexc=pert corr to rhov*exc
C          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),smcor,
C     .      x1,x2)
        enddo
C       Peturbation correction to smvxc
        call smvxc3(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,smcor,
     .    dxcv,smvxc,focvxc)
        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        if (iprint() >= 30) then
          outs = ' '
          call awrit8('%x   foca'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N   foca#%2j#'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',
     .      outs,len(outs),0,
     .      focexc(1)+focexc(2),nsp,focexc,focexc(2),
     .      focvxc(1)+focvxc(2),nsp,focvxc,focvxc(2))
          call awrit1('%a  charge  =%;12,6D',outs,len(outs),-stdo,sum1)
        endif
        deallocate(dxcv)
      endif

      deallocate(smrhol,cgh1,cgh2)

C --- Force from foca sm-head; cgh1 is workspace ---
      if (lfrce /= 0) then
      allocate(cgh1(ng*nsp))
      call dpzero(f,3*nbas)
      if (lfoc1 > 0 .or. lfoc2 > 0) then
        call fftz3(smvxc,n1,n2,n3,k1,k2,k3,nsp,0,-1)
        call gvgetf(ng,nsp,s_lat%kv,k1,k2,k3,smvxc,cgh1)
        call smvxc4(nbas,nsp,s_site,s_spec,
     .    alat,vol,s_lat%cy,ng,s_lat%gv,cgh1,f)
      endif
      deallocate(cgh1)
      endif

C     Debugging
C      print *, bscal
C      print *, dble(smvxc(1,1,1,1)+smvxc(1,1,1,2)),
C     .    dble(smvxc(1,1,1,1)-smvxc(1,1,1,2))

      call tcx('smvxc')

      end

      subroutine smvxc2(mode,s_lat,nsp,lxcfun,bscal,vol,n1,n2,n3,k1,
     .  k2,k3,smrho,smvxc,smvx,smvc,smexc,dsmvxc,rhoeps,rhoex,rhoec,
     .  rhomu,vxcavg)
C- Makes smooth local part of xc potential smvxc and optionally dsmvxc/drho
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not make dsmvxc/drho
Ci         :1 make dsmvxc/drho
Ci         :10s digit
Ci         : 0 calculated LDA for density as is.
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci   nsp   :number of spin channels
Ci   lxcfun:switch defining xc functional; see lxcf in structures.h
Ci   bscal :scale vxc+ - vxc- by bscal (nsp=2 only)
Ci   vol   :cell volume
Ci   slat  :struct containing information about the lattice
Ci   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :smooth density on uniform mesh
Cl Local variables
Cl   lxcfl :If nonzero, evaluate a local XC functional through vxcnsl
Cl         :If >1000, the libxc library is used.
Cl   lxcfnl:If nonzero, indicates a GGA from from the libxc library
Cl         :is to be evaluated through vxcnls
Cl   lxcg  :If nonzero, flags that this functional requires gradients
Cl         :If not a libxc functional, lxcg specifies which GGA
Co Outputs
Co   smvxc :xc potential of smoothed density (no core contr.)
Co   smvx  :exchange potential of smoothed density + core corrections
Co   smvc  :correlation potential of smoothed density + core corrections
Co   dsmvxc:dvxc/drho (mode=1)
Co   rhoeps:integrated exchange-correlation energy
Co   rhoex :integrated exchange energy
Co   rhoec :integrated correlation energy
Co   rhomu :integrated exchange-correlation potential
Co   vxcavg:average xc potential
Cr Remarks
Cr   For perturbation treatment, take numerical derivatives
Cr   df/dr = d/dr (vxc*r**alfa) instead of d/dr vxc because
Cr   f is nearly linear for alpha=2/3.
Cr
Cr   In the spin polarized case, the smooth core density is not
Cr   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
Cr   we use dm=0 and thus drho1 = drho2 = drho/2; thus
Cr     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
Cr
Cr   If vxc+ - vxc- is scaled by bscal, terms exc and integrate quantities
Cr   such as rhoeps, rhomu, lose their meaning. The potential is merely
Cr   scaled; no attempt is made to construct a functional corresponding
Cr   to the scaling.
Cu Updates
Cu   20 Nov 09 New 10s digit for mode
Cu   21 Apr 09 Handles GGA functionals
Cu   14 Jun 02 rhoex and rhoec (T. Miyake)
Cu    8 Feb 02 smvx and smvc (T. Miyake)
Cu   12 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp vxcd_smooth.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nsp,k1,k2,k3,n1,n2,n3,lxcfun
      double precision rhoeps(2),rhoex(2),rhoec(2),rhomu(2),
     .                 vxcavg(2),vol,bscal
      double complex smvxc(k1,k2,k3,2),smrho(k1,k2,k3,2),
     .               smvx(k1,k2,k3,2),smvc(k1,k2,k3,2),
     .               smexc(k1,k2,k3),dsmvxc(k1,k2,k3,2)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      integer i,i1,i2,i3,nx,stdo,iprint,n1x,nglob,mode1
      integer xctype(2),lxcfl,lxcg,lxcfnl
      parameter (n1x=512)
      double precision alfa,dfdr,dvdr,f,f1,f2,rrho,fac,dmach,rrmin
      double precision repnl(2),rmunl(2),vavgnl(2),rvx(2),rvc(2)
      double precision rho(n1x),rhos(n1x,2),
     .                 vxc2(n1x,2),vxc1(n1x,2),
     .                 vx1(n1x,2),vc1(n1x,2),
     .                 exc2(n1x),exc1(n1x),
     .                 ex1(n1x),ec1(n1x),xx(1)
      character*180 outs

C ... Setup
      if (n1 > n1x) call rxi('smvxc2: increase n1x, need',n1)
      stdo = nglob('stdo')
      if (lxcfun >= 1000) then   ! a libxc functional.  Check if need gradients
        call xc_libxc_type(lxcfun,0,xctype)
        if (xctype(1) > 2 .or. xctype(2) > 2)
     .    call rx('smvxcm: functional not implemented')
        if (xctype(1) == 2 .or. xctype(2) == 2) then ! This is GGA
          lxcfl = 0        ! No potential calculated by vxcnsl (LDA)
          lxcfnl = lxcfun  ! Potential calculated by vxcnls (GGA)
          lxcg = 1         ! Require gradients
        else ! This is LDA
          lxcfl = lxcfun ! Evaluate XC functional through vxcnsl (LDA)
          lxcfnl = 0     ! No potential calculated by vxcnsl (GGA)
          lxcg = 0       ! No gradients
        endif
      else
        lxcfl = mod(lxcfun,100) ! Evaluate XC functional through vxcnsl
        lxcfnl = 0              ! No libxc potential calculated by vxcnls
        lxcg = mod(lxcfun/100,100) ! nonlocal part of GGA by vxcnls
      endif

      alfa = 2d0/3d0  ! For numerical differentiation of rho, pert treatment
      fac = dmach(1)**(1d0/3d0) ! Size of finite difference
      rhoeps=0; rhoex=0; rhoec=0; rhomu=0; vxcavg=0; rvx=0; rvc=0
      repnl=0; rmunl=0; vavgnl=0; rrmin=0; nx=0
      mode1 = mod(mode/10,10)
      call dpzero(smvxc,2*k1*k2*k3*nsp)
      call dpzero(smvx,2*k1*k2*k3*nsp)
      call dpzero(smvc,2*k1*k2*k3*nsp)
      call dpzero(smexc,2*k1*k2*k3)

      if (lxcfl == 0) goto 100 ! No local functional; GGA only

C ... LD treatment: do by rows
      do  i3 = 1, n3
        do  i2 = 1, n2

C     ... rho = total rho for this vec; rhos = spin pol rho
          call dcopy(n1,smrho(1,i2,i3,1),2,rho,1)
          call dcopy(n1,smrho(1,i2,i3,1),2,rhos,1)
          if (nsp == 2) then
            call dcopy(n1,smrho(1,i2,i3,2),2,rhos(1,2),1)
            call daxpy(n1,1d0,rhos(1,2),1,rho,1)
          endif

C     ... Perturbation treatment: put df/dr into dsmvxc
          if (mod(mode,10) /= 0) then
            if (lxcfl < 1 .or. lxcfl > 4 .or. lxcg /= 0)
     .        call rx('smvxc3: perturbation treatment '//
     .        'not implemented with this functional')
C           Add rho*fac/2 into rhos+, rhos- and fac*rho into rho
            do  i = 1, nsp
              call daxpy(n1,fac/2,rho,1,rhos(1,i),1)
            enddo
            call dscal(n1,1+fac,rho,1)
C           Exchange potential at rho+drho
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,1,lxcfl,
     .          exc2,ex1,ec1,vxc2(1,i),vx1(1,i),vc1(1,i))
            enddo
            if (mode1 == 2) then
              do  i1 = 1, n1
                if (rhos(i1,1) < 0 .or. rhos(i1,nsp) < 0) then
                  vxc2(i1,1) = 0
                  vxc2(i1,nsp) = 0
                endif
              enddo
            endif
C           Restore rho,rhos; also add -drho to rho and -drho/2 to rhos
            call dscal(n1,1/(1+fac),rho,1)
            do  i = 1, nsp
              call daxpy(n1,-fac,rho,1,rhos(1,i),1)
            enddo
            call dscal(n1,(1-fac),rho,1)
C           Exchange potential at rho-drho
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,1,lxcfl,
     .          exc1,ex1,ec1,vxc1(1,i),vx1(1,i),vc1(1,i))
            enddo
            if (mode1 == 2) then
              do  i1 = 1, n1
                if (rhos(i1,1) < 0 .or. rhos(i1,nsp) < 0 .or.
     .              vxc2(i1,1) == 0 .or. vxc2(i1,nsp) == 0) then
                  vxc1(i1,1) = 0
                  vxc1(i1,nsp) = 0
                  vxc2(i1,1) = 0
                  vxc2(i1,nsp) = 0
                endif
              enddo
            endif

C           Restore rho,rhos
            call dscal(n1,1/(1-fac),rho,1)
            do  i = 1, nsp
              call daxpy(n1,fac/2,rho,1,rhos(1,i),1)
            enddo

            do  i = 1, nsp
            do  i1 = 1, n1
              if (rho(i1) > 0) then
                f1 = vxc1(i1,i)*(rho(i1)*(1-fac))**alfa
                f2 = vxc2(i1,i)*(rho(i1)*(1+fac))**alfa
                dfdr = (f2-f1)/(2d0*fac*rho(i1))
                vxc2(i1,i) = dfdr
              else
                vxc2(i1,i) = 0
              endif
            enddo
            enddo
          endif

C     ... Exchange into smvxc
          if (lxcfl > 1000) then ! A local libxc functional
            call xc_libxc(n1,nsp,-lxcfl,rhos(1,1),rhos(1,nsp),
     .        xx,xx,xx,xx,xx,xx,xx,xx,xx,xx,ex1,ec1,vx1,vc1,
     .        vx1(1,nsp),vc1(1,nsp),exc1,vxc1,vxc1(1,nsp),xx)
          elseif (lxcfl == 3 .or. lxcfl == 4) then
            call evxcp(rhos(1,1),rhos(1,nsp),n1,nsp,lxcfl,ex1,ec1,exc1,
     .        vx1(1,1),vx1(1,2),vc1(1,1),vc1(1,2),vxc1(1,1),vxc1(1,2))
          else
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,nsp,lxcfl,exc1,ex1,ec1,
     .          vxc1(1,i),vx1(1,i),vc1(1,i))
            enddo
C          print *, sngl(exc1(n1)),sngl(vxc1(n1,1)),sngl(vxc1(n1,2))
          endif
          if (mode1 == 2) then   ! Treat negative densities
            do  i = 1, nsp
            do  i1 = 1, n1
              if (rho(i1) <= 0 .or. rhos(i1,i) <= 0) then
                vxc1(1,i) = 0
              endif
            enddo
            enddo
          endif

C         call prmx('vxc',vxc1(1,1),n1,n1,1)
C         call prmx('vxc',vxc1(1,2),n1,n1,1)
C         call prmx('exc',exc1,n1,n1,1)
          do  i = 1, nsp
            call dcopy(n1,vxc1(1,i),1,smvxc(1,i2,i3,i),2)
            call dcopy(n1,vx1(1,i),1,smvx(1,i2,i3,i),2)
            call dcopy(n1,vc1(1,i),1,smvc(1,i2,i3,i),2)
            call dcopy(n1,exc1,1,smexc(1,i2,i3),2)
          enddo

C     ... Perturbation dv/dr into dsmvxc
          if (mod(mode,10) /= 0) then
            do  i = 1, nsp
              do  i1 = 1, n1
                rrho = rho(i1)
                if (rrho > 0) then
                  f = vxc1(i1,i) * rrho**alfa
                  dvdr = (vxc2(i1,i) - alfa*f/rrho) / rrho**alfa
                  dsmvxc(i1,i2,i3,i) = dvdr
                else
                  dsmvxc(i1,i2,i3,i) = 0
                endif
              enddo
            enddo
          endif

C     ... Add to integrals
          do  i = 1, nsp
            do  i1 = 1, n1
              rrho = rhos(i1,i)
              rrmin = min(rrho,rrmin)
              if (rrho < 0d0) nx = nx+1
              rhomu(i)  = rhomu(i)  + rrho*vxc1(i1,i)
              rhoeps(i) = rhoeps(i) + rrho*exc1(i1)
              rhoex(i)  = rhoex(i)  + rrho*ex1(i1)
              rhoec(i)  = rhoec(i)  + rrho*ec1(i1)
              vxcavg(i) = vxcavg(i) + vxc1(i1,i)
            enddo
          enddo
        enddo
      enddo

      f = vol/(n1*n2*n3)
      do  i = 1, nsp
        rhoeps(i) = rhoeps(i)*f
        rhoex(i) = rhoex(i)*f
        rhoec(i) = rhoec(i)*f
        rhomu(i) = rhomu(i)*f
        vxcavg(i) = vxcavg(i)/(n1*n2*n3)
      enddo

C      do  i = 1, nsp
C        call zprm3('LDA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
C        call zprm3('dsmvxc (isp=%i)',i,dsmvxc(1,1,1,i),n1,n2,n3)
C      enddo

C ... Gradient corrections to potential, energy
  100 continue
C     If lxcfnl>0, returns full exchange, correlation potentials, and integrals
C                  repnl,rmunl,vavgnl,smvx,smvc,smvxc are complete.
      if (lxcg /= 0) then
        call vxcnlm(lxcfnl,lxcg,nsp,k1,k2,k3,s_lat,smrho,
     .    rhoex,rhoec,repnl,rvx,rvc,rmunl,vavgnl,smvx,smvc,smvxc,smexc)
      endif
C      do  i = 1, nsp
C        call zprm3('LDA+GGA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
C      enddo

C     call setpr(30)

      if (nsp == 2 .and. bscal /= 1) then
        i = 2*k1*k2*k3  ! potential arrays are complex
        call bxcscale(i,nsp,bscal,i,smvxc)
        call bxcscale(i,nsp,bscal,i,smvx)
        call bxcscale(i,nsp,bscal,i,smvc)
        if (mod(mode,10) /= 0) then
          call bxcscale(i,nsp,bscal,i,dsmvxc)
        endif
C       call bxcscale(1,nsp,bscal,1,rhomu)  ! Not meaningful
        call bxcscale(1,nsp,bscal,1,vxcavg)
        if (lxcg /= 0) then
          call bxcscale(1,nsp,bscal,1,vavgnl)
        endif
      endif

C ... LDA, GGA Printout
      if (nx > 0) call info5(20,0,0,' smvxcm (warning) mesh density '
     . //'negative at %i point%?#n>1#s##:  rhomin=%;3g',nx,nx,rrmin,0,0)

      if (iprint() >= 30 .and. lxcg /= 0 .and. lxcfnl == 0) then
        call awrit8('%x sm LDA'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm LDA'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .    0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
     .    nsp,rhomu,rhomu(2))
        call awrit5('%a%?#n==2#%N%3f#  #avg LDA vxc ='//
     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .    -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))

        call awrit8('%x sm GGA'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm GGA'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .    0,repnl(1)+repnl(2),nsp,repnl,repnl(2),rmunl(1)+rmunl(2),
     .    nsp,rmunl,rmunl(2))
        call awrit5('%a%?#n==2#%N%3f#  #avg GGA vxc ='//
     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .    -stdo,nsp,(vavgnl(1)+vavgnl(nsp))/2,nsp,vavgnl,vavgnl(2))
      endif

      if (lxcg /= 0) then
        do  i = 1, nsp
          rhoeps(i) = rhoeps(i) + repnl(i)
          rhomu(i) = rhomu(i) + rmunl(i)
          vxcavg(i) = vxcavg(i) + vavgnl(i)
        enddo
      endif

C ... Printout, total potential
      if (iprint() >= 30) then
        call awrit8('%x smooth'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%7f#%2j#'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .    0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
     .    nsp,rhomu,rhomu(2))
        call awrit5('%a%?#n==2#%N%7f#  #'//
     .    'avg vxc =%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .    -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))
      endif

      end

      subroutine smvxc3(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,smcor,dsmvxc,
     .  smvxc,rmuxcc)
C- Smooth core density times dvxc/drho
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :smooth density on n1,n2,n3 mesh
Ci   smcor :smooth core density on n1,n2,n3 mesh
Ci   dsmvxc:dvxc/drho on n1,n2,n3 mesh mesh
Co Outputs
Co   smvxc :(dvxc/drho * smcor) = pert. correction to expansion
Co         : vxc[rho + rhoc] = vxc[rho] + rhoc * dvxc/drho
Co   rmuxcc:integral smrho * (dvxc/drho * smcor)
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3,n1,n2,n3
      double precision rmuxcc(nsp),vol
      double complex smvxc(k1,k2,k3,nsp),smcor(k1,k2,k3),
     .  dsmvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer i,i1,i2,i3
      double complex cadd,csum(2)

      rmuxcc(2) = 0
      do  i = 1, nsp
        csum(i) = 0d0
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              cadd = dsmvxc(i1,i2,i3,i)*smcor(i1,i2,i3)
              smvxc(i1,i2,i3,i) = cadd
              csum(i) = csum(i) + smrho(i1,i2,i3,i)*cadd
            enddo
          enddo
        enddo
        csum(i) = csum(i)*vol/(n1*n2*n3)
        rmuxcc(i) = dble(csum(i))
      enddo

C     write(stdo,862) csum
C 862 format(' csum=',2f14.8)

      end
      subroutine smvxc4(nbas,nsp,s_site,s_spec,alat,vol,cy,ng,gv,cvxc,f)
C- For foca, adds force from shift of smH-head against Vxc.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: corprm
Ci   cy    :Normalization constants for spherical harmonics
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cvxc  :Fourier transform of smooth vxc potential.
Co Outputs
Co   f     :force from shift of smH-head against Vxc added to f.
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   02 Jul 05  skip sites for which cofh=0
Cu    1 May 00  Adapted from nfp smc_force.f
C ----------------------------------------------------------------------
      use structures
      use mpi

      implicit none
C ... Passed parameters
      integer nbas,nsp,ng
      double precision gv(ng,3),alat,vol,cy(*),f(3,nbas)
      double complex cvxc(ng,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer stdo,k0,nlmx,nglob,kmax,ib,is,lfoc,i,kb,iprint
      double precision tau(3),v(3),pi,tpiba,qcorg,qcorh,qsc,cofg,
     .  cofh,ceh,rfoc,z,sum1,sum2,sum3,xx
      parameter (k0=3, nlmx=9)
      double complex gkl(0:k0,nlmx),ccc,cvxci

      integer :: iproc, nproc, comm, err, fl_mt, idx(nbas)
      integer, allocatable :: pmap(:), smap(:)
      real(8) :: fl(4,nbas), flav(4)

      call tcn('smvxc4')

      stdo = nglob('stdo')
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat
      kmax = 0

C --- Loop over sites ---
      if (iprint() >= 50) write(stdo,400)

      comm = mpi_comm_world
      call mpi_comm_size(comm, nproc, err)
      call mpi_comm_rank(comm, iproc, err)
      allocate(pmap(0:nproc))
      idx = [(ib, ib=1,nbas)]
      call vbdist(nbas, idx, nproc, pmap)

      flav = 0
      do ib = pmap(iproc)+1, pmap(iproc+1)
        fl(:,ib) = 0
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        if (lfoc > 0 .and. cofh /= 0) then
          do  i = 1, ng
            v = gv(i,1:3)
            call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy,gkl)
            ccc = cofh*gkl(0,1)/vol
            cvxci = 0.5d0 * (cvxc(i,1) + cvxc(i,nsp))
            xx = -dimag(dconjg(cvxci) * ccc)
            fl(1:3,ib) = fl(1:3,ib) + xx*gv(i,1:3)
          enddo
          fl(:,ib) = fl(:,ib)*vol*tpiba
          flav = flav + fl(:,ib)
        endif
      enddo

      if (nproc > 1) call mpi_allreduce(mpi_in_place, flav, 4, mpi_real8, mpi_sum, comm, err)

      flav = flav/nbas

      do ib = pmap(iproc)+1, pmap(iproc+1)
        fl(:,ib) = fl(:,ib) - flav
      end do

      if (nproc > 1) then
        allocate(smap(0:nproc-1))
        smap = pmap(1:nproc) - pmap(0:nproc-1)

        call mpi_type_contiguous(4, mpi_real8, fl_mt, err)
        call mpi_type_commit(fl_mt, err)
        call mpi_allgatherv(mpi_in_place, smap(iproc), fl_mt, fl, smap, pmap, fl_mt, comm, err)
        call mpi_type_free(fl_mt, err)
      end if

      do ib = 1, nbas
        f(1:3,ib) = f(1:3,ib) + fl(1:3,ib)
      end do

      if (iprint() >= 50)
     .  write(stdo,340) (ib,f(1,ib),f(2,ib),f(3,ib),ib = 1,nbas)
  340 format(i4,3f12.6)
  400 format(/' xc-force from foca:')

      call tcx('smvxc4')

      end
