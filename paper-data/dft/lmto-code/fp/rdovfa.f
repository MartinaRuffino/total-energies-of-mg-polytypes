      subroutine rdovfa(nbas,nspec,s_site,s_spec,s_lat,s_pot,s_ham,qbg)
C- Read and overlap free atom densities.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec rho1 rho2 rhoc pos
Co     Stored:     *
Co     Allocated:  v0 v1 rho1 rho2 rhoc
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  dfratm ovlpfa ovlocr adbkql
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lfoca rfoca lmxl qc rhoc coreh coreq
Ci                 kmxv rsmv ctail etail stc lmxb p pz rg
Co     Stored:     a nr qc nxi exi chfa rsmfa ctail etail stc rhoc
Co     Allocated:  rhoc
Cio    Elts passed:*
Cio    Passed to:  bcast_strx dfratm gtpcor ovlocr corprm adbkql
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat vol nabc ng qlat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv kv cg indxcg jcg cy qlv dlv
Cio    Passed to:  ovlpfa ovlocr hxpbl ghibl hklbl gklbl
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:smrho
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   qbg   :constant background charge
Co Outputs
Co   smrho :smoothed interstitial density
Co         :* for smrho = smoothed mesh density, smrho is complex and
Co         :  smrho = smrho(k1,k2,k3)
Cl Local variables
Cl   k1,k2,k3 dimensions of smrho for smoothed mesh density
Cl   ns4 = 1 if local density is not spin polarized
Cl         2 if local spin density is collinear along z
Cl         4 if local spin density is to be rotated;
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 May 07 package mpi-specific calls
Cu   02 Jan 06 generates core magnetic moment, checks against spec->qcor
Cu   01 Jul 05 Zero-radius empty spheres treated as local orbitals
Cu   12 Apr 03 (WRL) Added constant charge background
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   13 Jun 00 spin polarized
Cu   21 Apr 00 Adapted from nfp rdovfa.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nspec
      double precision qbg
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      real(8), allocatable  :: v0(:),v1(:)
      real(8), allocatable :: rwgt(:)
      complex(8), allocatable :: cv(:)
      type(dp_wk_vec), allocatable :: s_rhofa(:),s_rhoca(:),s_v0a(:)
C ... Local parameters
      integer procid, master, mpipid
      integer nrmx, n0
      parameter ( nrmx=5001, n0=10 )
      integer nxi(nspec)
      double precision rsmfa(nspec),pnu(n0,2),exi(n0,nspec),
     .  hfc(n0,2,nspec),hfct(n0,2,nspec),plat(3,3),qcor(2),dlength
      character*8 spid(nspec),spidr
      integer fopna,fxst,i,i1,ib,ifi,intopt,iofa,ipr,iprint,is,k1,k2,k3,
     .  kcor,lcor,lfoc,lgunit,lmxl,n1,n2,n3,nch,ng,ngabc(3),nglob,nlml,
     .  nr,nr0,nsp,ns4,nspc,stdl,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision a,a0,alat,ccof,ceh,corm,ctot,dq,fac,qc,rfoc,rmt,
     .  rmt0,slmom(0:3),smom(0:3),sqloc,stc,sum,sum1,sum2,vol,xx,z,z0,
     .  ztot
      character msg*23, strn*120
      logical mlog,cmdopt,lfail

      call tcn('rdovfa')
      ipr   = iprint()
      stdo  = nglob('stdo')
      stdl  = lgunit(2)
      nsp   = nglob('nsp')
      nspc = nglob('nspc')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10)==0) ns4 = nsp ! See "Local variables"

      msg   = '         File mismatch:'
      intopt = 10*nglob('lrquad')
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

      call info0(10,1,0,' rdovfa: read and overlap free-atom densities'
     .  //' (mesh density) ...')
      if (mod(s_ham%lrsa,10) /= 0) call info2(10,0,0,
     .  '%9fRotate atomic densities by Euler angles'//
     .  '%?#n>0#;  M parallel to |M|##',mod(s_ham%lrsa/10,10),0)
      if (s_ham%lrsa >= 100) call info0(10,0,0,
     .  '%9fspin average core densities')

      alat = s_lat%alat
      plat = s_lat%plat
      vol = s_lat%vol
      ngabc = s_lat%nabc
      call fftz30(n1,n2,n3,k1,k2,k3)
      call dpzero(hfc,n0*2*nspec)
      call dpzero(pnu,n0*2)

C  --- Read free-atom density for all species ---
C      if (ipr >= 30) write(stdo,*) ' '
      if (procid == master) then
        if (fxst('atm') /= 1) then
          call info0(10,0,0,' rdovfa: '//
     .      'atm file missing or unreadable ... aborting')
          lfail = .true.
        else
          ifi = fopna('atm',-1,0)
          rewind ifi
          lfail = .false.
        endif
      endif
      call mpibc1(lfail,1,1,mlog,'rdovfa','read error')
      if (lfail) call rx('missing atm file')

      allocate(s_rhofa(nspec),s_rhoca(nspec),s_v0a(nspec))
      do  is = 1, nspec
        allocate(s_rhofa(is)%p(nrmx*nsp))
        call dpzero(s_rhofa(is)%p,nrmx*nsp)
        allocate(s_rhoca(is)%p(nrmx*nsp))
        call dpzero(s_rhoca(is)%p,nrmx*nsp)
        allocate(s_v0a(is)%p(nrmx*nsp))
        call dpzero(s_v0a(is)%p,nrmx*nsp)

        spid(is) = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        z = s_spec(is)%z
        lfoc = s_spec(is)%lfoca
        rfoc = s_spec(is)%rfoca
        lfail = .false.
        if (procid == master) then
        if (z == 0 .and. rmt == 0) then
          nxi(is) = 0
          call dpzero(exi(1,is),n0)
          call dpzero(hfc(1,1,is),2*n0)
          call dpzero(hfct(1,1,is),2*n0)
          rsmfa(is) = 0
          z0=0
          rmt0=0
          a0=0
          nr0=0
          qc=0
          ccof=0
          ceh=0
          stc=0
        else
          lfail = iofa(63-2,spidr,n0,nxi(is),exi(1,is),hfc(1,1,is),
     .      hfct(1,1,is),rsmfa(is),z0,rmt0,a0,nr0,nsp,qc,ccof,
     .      ceh,stc,s_rhofa(is)%p,s_rhoca(is)%p,s_v0a(is)%p,ifi) < 0
        endif
        endif
        call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
        call mpibc1(qc,1,4,mlog,'rdovfa','qc')
        call mpibc1(ccof,1,4,mlog,'rdovfa','ccof')
        call mpibc1(ceh,1,4,mlog,'rdovfa','ceh')
        call mpibc1(stc,1,4,mlog,'rdovfa','stc')
        call mpibc1(lfail,1,1,mlog,'rdovfa','read error')
        if (lfail) call rxs('missing species data, species ',spid(is))

C   ... Broadcast file data
        call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
        call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
        call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
        call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
        call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
        call mpibc1(a0,1,4,mlog,'rdovfa','a0')
        if (nr0 > 0) then
        call mpibc1(s_rhofa(is)%p,nr0*nsp,4,mlog,'rdovfa','rhofa')
        call mpibc1(s_rhoca(is)%p,nr0*nsp,4,mlog,'rdovfa','rhoca')
        call mpibc1(s_v0a(is)%p,nr0*nsp,4,mlog,'rdovfa','v0a')
        endif
        i = mpipid(2)

C ...   Defaults
        if (procid == master) then
          call strip(spid(is),i1,nch)
          if (ipr >= 30 .and. rmt0 /= 0)
     .      write(stdo,400) spid(is)(1:nch),spidr,rmt0,nr0,a0
  400     format(' rdovfa: expected ',a,',',T27,' read ',a,
     .           ' with rmt=',f8.4,'  mesh',i6,f7.3)
        endif
        if (nr <= 0)   nr = nr0
        if (a <= 1d-6) a = a0
        if (z == 0 .and. rmt == 0) then
          a = 0
          nr = 0
        endif

C ...   Sanity checks
        if (procid == master) then
          call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
          call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
          call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
          call sanrg(.true.,nr0,nr,nr,msg,'nr')
        endif

C ...   Pack into sspec
        s_spec(is)%a = a
        s_spec(is)%nr = nr
        s_spec(is)%qc = qc
        s_spec(is)%nxi = nxi(is)
        s_spec(is)%exi = exi(:,is)
        s_spec(is)%chfa = hfc(:,:,is)
        s_spec(is)%rsmfa = rsmfa(is)
        s_spec(is)%ctail = ccof
        s_spec(is)%etail = ceh
        s_spec(is)%stc = stc
        if (nr > 0) then
          call ptr_spec(s_spec,4+1,'rhoc',is,nr*nsp,0,s_rhoca(is)%p)
        endif
      enddo
C     Wait for all proccesses to synchronize
      i = mpipid(2)
C     Re-broadcast entire species structure, and arrays used below
      call bcast_strx(1,xx,xx,xx,xx,xx,xx,s_spec,xx,xx,nspec,0)

      if (procid == master) then
        call fclose(ifi)
      endif

C --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
      ztot = 0d0
      ctot = 0d0
      corm = 0d0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        lmxl = s_spec(is)%lmxl
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        lfoc = s_spec(is)%lfoca

        nlml = (lmxl+1)**2
        call dfratm(s_site,s_spec,1,ib,ib,xx)
        allocate(v0(nr*nsp),v1(nr*nsp))

C       Core magnetic moment (possible if magnetized core hole)
        if (nsp == 2 .and. lmxl > -1) then
          allocate(rwgt(nr))
          call radwgt(intopt,rmt,a,nr,rwgt)
          call radsum(nr,nr,1,nsp,rwgt,s_rhoca(is)%p,sum)
          call radsum(nr,nr,1,1,rwgt,s_rhoca(is)%p,sum1)
          sum2 = sum - sum1
          call gtpcor(s_spec,is,kcor,lcor,qcor)

          if (dabs(qcor(2)-(sum1-sum2)) > 0.01d0) then
            call info5(10,0,0,' (warning) core moment mismatch spec %i:'
     .        //'  input file=%;6d  atom file=%;6d',
     .        is,qcor(2),sum1-sum2,0,0)
          endif
          corm = corm + qcor(2)
          deallocate(rwgt)
        endif

        if (lmxl > -1) then
          call dpcopy(s_v0a(is)%p,  v0,  1,nr*nsp,1d0)
          call dpcopy(s_v0a(is)%p,  v1,  1,nr*nsp,1d0)
          call dpcopy(s_rhoca(is)%p,s_site(ib)%rhoc,1,nr*nsp,1d0)

          if (lfoc == 0) then
            allocate(rwgt(nr))
            call radwgt(intopt,rmt,a,nr,rwgt)
            call radsum(nr,nr,1,nsp,rwgt,s_site(ib)%rhoc,sum)
            fac = 1d0
            if(dabs(sum) > 1d-7) fac = qc/sum
            if (ipr >= 40) write(stdo,787) is,qc,sum,fac
  787       format(' scale foca=0 core species',i2,': qc,sum,scale=',
     .        3f12.6,f12.6)
            call dpcopy(s_site(ib)%rhoc,s_site(ib)%rhoc,1,nr*nsp,fac)
            deallocate(rwgt)
          endif
        endif

        call ptr_site(s_site,4+1,'v0',ib,nr,nsp,v0)
        call ptr_site(s_site,4+1,'v1',ib,nr,nsp,v1)
        ztot = ztot+z
        ctot = ctot+qc

C       end loop over sites
        deallocate(v0,v1)
      enddo

C --- Overlap smooth hankels to get smooth interstitial density ---
      ng = s_lat%ng
      allocate(cv(ng*ns4))
      call ovlpfa(s_site,s_lat,s_ham,1,nbas,nxi,n0,exi,hfc,rsmfa,ng,ng,
     .  s_lat%gv,cv)
      call gvputf(ng,ns4,s_lat%kv,k1,k2,k3,cv,s_pot%smrho)
      deallocate(cv)

C ... FFT to real-space mesh
      call fftz3(s_pot%smrho,n1,n2,n3,k1,k2,k3,ns4,0,1)
      if (ns4==4) then ! Transpose to packed format
        i = k1*k2*k3
        call rotspv(251,xx,i,i,i,xx,s_pot%smrho,xx,s_pot%smrho)
      endif

C ... Add compensating uniform electron density to compensate background
      call addbkgsm(s_pot%smrho,k1,k2,k3,nsp,qbg,vol,-1d0)
C ... integrate for mesh density and magnetic moment
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,s_pot%smrho,sum1,sum2)
      smom = 0
      if (nsp == 2) then  ! Collinear spin polarized case
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho,smom,xx)
        smom(3) = 2*smom(0) - sum1
        smom(0) = smom(3)
      endif
      if (ns4 == 4) then  ! Noncollinear case
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho(1,3),smom(1),xx)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho(1,4),smom(2),xx)
        smom(1) = 2*smom(1)
        smom(2) =-2*smom(2)
        smom(0) = dlength(3,smom(1),1)
      endif

C --- Set up local densities using rmt from atm file ---
      call ovlocr(nbas,s_site,s_spec,s_lat,s_ham,n0,nxi,exi,
     .            hfc,rsmfa,s_rhofa,sqloc,slmom)
C --- Add compensating uniform electron density to compensate background
      call adbkql(nbas,nsp,qbg,vol,-1d0,s_spec,s_site)
      if (abs(qbg) /= 0) call info(10,0,0,' Uniform '//
     .  'density added to neutralize background, q=%;6,6d',qbg,0)

      do  is = 1, nspec
        deallocate(s_rhofa(is)%p)
        deallocate(s_rhoca(is)%p)
        deallocate(s_v0a(is)%p)
      enddo
      deallocate(s_rhofa,s_rhoca,s_v0a)

C --- Print charges ---
      dq = sum1+sqloc+ctot-ztot+qbg
      if (nsp == 1 .and. ipr >= 10) then
        write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
  895   format(/' Smooth charge on mesh:    ',f16.6
     .     /    ' Sum of local charges:     ',f16.6
     .     /    ' Total valence charge:     ',f16.6
     .     /    ' Sum of core charges:      ',f16.6
     .     /    ' Sum of nuclear charges:   ',f16.6
     .     /    ' Homogeneous background:   ',f16.6
     .     /    ' Deviation from neutrality:',f16.6)
        write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
  710   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .    '   bg',f10.6,'  dQ',f10.6)
      elseif (ipr >= 10) then
        write(stdo,896) sum1,smom(0),sqloc,slmom(0)
        if (ns4 == 4) write(stdo,898) smom(1:3), slmom(1:3)
        write(stdo,897)
     .    sum1+sqloc,smom(0)+slmom(0),ctot,corm,-ztot,qbg,dq
  896   format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6/
     .          ' Sum of local charges:     ',f16.6,4x,'moments',f11.6)
  897   format( ' Total valence charge:     ',f16.6,4x,'moment', f12.6,
     .     /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6,
     .     /    ' Sum of nuclear charges:   ',f16.6
     .     /    ' Homogeneous background:   ',f16.6
     .     /    ' Deviation from neutrality:',f16.6)
        write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom(0)+slmom(0)
  898   format(' Mesh magnetization (xyz):',6x,3f11.6/
     .         ' Local magnetization (xyz):',5x,3f11.6)
  711   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .    '   bg',f11.6,' mm',f11.6)
      endif

      if (dabs(dq) > 1d-4 .and. ipr > 1)
     .  call awrit1(' rdovfa (warning) overlapped'
     .  //' density not neutral'//', dq= %d',' ',80,stdo,dq)

      call tcx('rdovfa')
      end
