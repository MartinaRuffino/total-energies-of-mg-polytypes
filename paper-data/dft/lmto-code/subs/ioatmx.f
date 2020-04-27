      subroutine ioatmx(opt,ib1,ib2,nspec,s_site,s_spec,s_ham,ifi)
C- Read and overlap free atom densities.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec rho1
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rhoc v0
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr rmt lmxl z qc lfoca
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed: name nxi exi chfa rsmfa z rmt a nr qc ctail etail stc
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrsa
Cio    Passed to:  *
Ci Inputs
Ci   opt   :For file read:
Ci         :0 : nothing --- don't read
Ci         :>0 : read for FP
Ci         :<0 : read for ASA
Ci         :FP case:
Ci         :1  : read s_spec(is)%z,s_spec(is)%rmt,s_spec(is)%a,s_spec(is)%nr,s_spec(is)%qc
Ci         :   : If they are not read, they much match passed data
Ci         :2  : read rhoc into s_site(ib)%rhoc and core K.E. into s_spec(is)%stc
Ci         :4  : read interstitial interstitial fitting data
Ci               s_spec(is)%nxi,s_spec(is)%exi,s_spec(is)%chfa,s_spec(is)%rsmfa,s_spec(is)%ctail,s_spec(is)%etail
Ci         :8  : read spherical part of density, store in s_site(ib)%rho1
Ci         :16 : read v0, store in s_site(ib)%v0
Ci         :switches may be taken in combination
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   ifi   :|ifi|=file logical unit, with ifi>0 for read, ifi<0 for write
Cl Local variables
Cl   ns4 = 1 if local density is not spin polarized
Cl         2 if local spin density is collinear along z
Cl         4 if local spin density is to be rotated;
Cr Remarks
Cu Updates
Cu   23 Dec 16 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,ib1,ib2,nspec,ifi
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      real(8), allocatable  :: rhov(:,:),rhoc(:,:),v0(:,:)
!      real(8), allocatable :: rwgt(:)
C ... Local parameters
      integer procid, master, mpipid
      integer nrmx, n0
      parameter ( nrmx=5001, n0=10 )
      double precision a,qc,rmt,z,srfpi,y0,rsmfa,ccof,ceh,stc
      double precision pnu(n0,2),exi(n0),hfc(n0,2),hfct(n0,2)
      character*8 spidr
      integer ib,intopt,iofa,ipr,is,lfoc,nlml,nr,nsp,ns4,nspc,
     .  stdl,stdo,jfi,isp,nxi
      character msg*23, strn*120
      logical mlog,cmdopt,lfail
      procedure(integer) :: nglob,iprint

      ipr   = iprint()
      stdo  = nglob('stdo')
      stdl  = nglob('stdl')
      nsp   = nglob('nsp')
      nspc = nglob('nspc')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"
      srfpi = dsqrt(4*4d0*datan(1d0))
      y0    = 1d0/srfpi
      intopt = 10*nglob('lrquad')

      msg   = '         File mismatch:'
      intopt = 10*nglob('lrquad')
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

C      call info0(10,1,0,' rdovfa: read and overlap free-atom densities'
C     .  //' (mesh density) ...')
C      if (mod(s_ham%lrsa,10) /= 0) call info2(10,0,0,
C     .  '%9fRotate atomic densities by Euler angles'//
C     .  '%?#n>0#;  M parallel to |M|##',mod(s_ham%lrsa/10,10),0)
C      if (s_ham%lrsa >= 100) call info0(10,0,0,
C     .  '%9fspin average core densities')

      call dpzero(hfc,n0*2*nspec)
      call dpzero(pnu,n0*2)

C --- Input ---
      if (ifi > 0) then

        if (opt == 0) return  ! Nothing to do
        if (opt < 0) call rx('ioatmx not ready for asa')

        do  ib = 1, ib2
          is = s_site(ib)%spec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          nlml = (s_spec(is)%lmxl+1)**2
          z = s_spec(is)%z
          qc = s_spec(is)%qc
          lfoc = s_spec(is)%lfoca

          allocate(rhov(nr,nr*nspc),rhoc(nr,nr*nspc),v0(nr,nr*nspc))
          lfail = (iofa(63-2,spidr,n0,nxi,exi,hfc,hfct,rsmfa,z,rmt,a,nr,qc,ccof,
     .      ceh,stc,rhov,rhoc,v0,ifi) < 0)

          call mpibc1(nr,1,2,mlog,'ioatmx','nr')
          call mpibc1(qc,1,4,mlog,'ioatmx','qc')
          call mpibc1(ccof,1,4,mlog,'ioatmx','ccof')
          call mpibc1(ceh,1,4,mlog,'ioatmx','ceh')
          call mpibc1(stc,1,4,mlog,'ioatmx','stc')
          call mpibc1(lfail,1,1,mlog,'ioatmx','read error')
          if (lfail) call rxi('missing atom data, site ',ib)

          if (ib >= ib1) then

C         Broadcast file data
          call mpibc1(nxi,1,2,mlog,'ioatmx','nxi')
          call mpibc1(exi,nxi,4,mlog,'ioatmx','exi')
          call mpibc1(hfc,nsp*n0,4,mlog,'ioatmx','hfc')
          call mpibc1(hfct,nsp*n0,4,mlog,'ioatmx','hfct')
          call mpibc1(rsmfa,1,4,mlog,'ioatmx','rsmfa')
          call mpibc1(a,1,4,mlog,'ioatmx','a')
          if (nr > 0) then
          call mpibc1(rhov,nr*nsp*nspc,4,mlog,'ioatmx','rhov')
          call mpibc1(rhoc,nr*nsp*nspc,4,mlog,'ioatmx','rhoc')
          call mpibc1(v0,nr*nsp*nspc,4,mlog,'ioatmx','v0')
          endif

C         Check, and optionally copy, sphere parameters
          if (mod(opt,2) == 1) then
            s_spec(is)%z = z
            s_spec(is)%rmt = rmt
            s_spec(is)%a = a
            s_spec(is)%nr = nr
            s_spec(is)%qc = qc
          endif
          call fsanrg(z,s_spec(is)%z,s_spec(is)%z,0d-9,msg,'z',.true.)
          call fsanrg(rmt,s_spec(is)%rmt,s_spec(is)%rmt,1d-6,msg,'rmt',.true.)
          call fsanrg(a,s_spec(is)%a,s_spec(is)%a,0d0,msg,'a',.true.)
          call sanrg(.true.,nr,s_spec(is)%nr,s_spec(is)%nr,msg,'nr')

C         Read rhoc into s_site(ib)%rhoc and core K.E. into s_spec(is)%stc
          if (mod(opt/2,2) == 1) then
            call dcopy(nr*nsp*nspc,rhoc,1,s_site(ib)%rhoc,1)
            s_spec(is)%stc = stc
          endif

C         Read interstitial interstitial fitting data
          if (mod(opt/4,2) == 1) then
            s_spec(is)%nxi = nxi
            s_spec(is)%exi = exi
            s_spec(is)%chfa = hfc
            s_spec(is)%rsmfa = rsmfa
            s_spec(is)%ctail = ccof
            s_spec(is)%etail = ceh
          endif

C         Read spherical part of density
          if (mod(opt/8,2) == 1) then
            do  isp = 1, nsp*nspc
              s_site(ib)%rho1(:,1+(isp-1)*nlml) = rhov(:,isp)*y0
            enddo
          endif

C         Read v0
          if (mod(opt/16,2) == 1) then
            call dcopy(nr*nsp*nspc,v0,1,s_site(ib)%v0,1)
          endif

          endif                   ! Read for this ib
          deallocate(rhov,rhoc,v0)

        enddo

CC  --- Read free-atom density for all species ---
CC      if (ipr >= 30) write(stdo,*) ' '
C      if (procid == master) then
C        if (fxst('atm') /= 1) then
C          call info0(10,0,0,' ioatmx: '//
C     .      'atm file missing or unreadable ... aborting')
C          lfail = .true.
C        else
C          ifi = fopna('atm',-1,0)
C          rewind ifi
C          lfail = .false.
C        endif
C      endif
C      call mpibc1(lfail,1,1,mlog,'ioatmx','read error')
C      if (lfail) call rx('missing atm file')
C
C      allocate(s_rhofa(nspec),s_rhoca(nspec),s_v0a(nspec))
C      do  is = 1, nspec
C        allocate(s_rhofa(is)%p(nrmx*nsp))
C        call dpzero(s_rhofa(is)%p,nrmx*nsp)
C        allocate(s_rhoca(is)%p(nrmx*nsp))
C        call dpzero(s_rhoca(is)%p,nrmx*nsp)
C        allocate(s_v0a(is)%p(nrmx*nsp))
C        call dpzero(s_v0a(is)%p,nrmx*nsp)
C
C        spid(is) = s_spec(is)%name
C        a = s_spec(is)%a
C        nr = s_spec(is)%nr
C        rmt = s_spec(is)%rmt
C        z = s_spec(is)%z
C        lfoc = s_spec(is)%lfoca
C        rfoc = s_spec(is)%rfoca
C        lfail = .false.
C        if (procid == master) then
C        if (z == 0 .and. rmt == 0) then
C          nxi(is) = 0
C          call dpzero(exi(1,is),n0)
C          call dpzero(hfc(1,1,is),2*n0)
C          call dpzero(hfct(1,1,is),2*n0)
C          rsmfa(is) = 0
C          z0=0
C          rmt0=0
C          a0=0
C          nr0=0
C          qc=0
C          ccof=0
C          ceh=0
C          stc=0
C        else
C          lfail =
C     .     (iofa(63-2,spidr,n0,nxi(is),exi(1,is),hfc(1,1,is),
C     .      hfct(1,1,is),rsmfa(is),z0,rmt0,a0,nr0,qc,ccof,
C     .      ceh,stc,s_rhofa(is)%p,s_rhoca(is)%p,s_v0a(is)%p,ifi)
C     . < 0)
C        endif
C        endif
C        call mpibc1(nr0,1,2,mlog,'ioatmx','nr0')
C        call mpibc1(qc,1,4,mlog,'ioatmx','qc')
C        call mpibc1(ccof,1,4,mlog,'ioatmx','ccof')
C        call mpibc1(ceh,1,4,mlog,'ioatmx','ceh')
C        call mpibc1(stc,1,4,mlog,'ioatmx','stc')
C        call mpibc1(lfail,1,1,mlog,'ioatmx','read error')
C        if (lfail) call rxs('missing species data, species ',spid(is))
C
CC   ... Broadcast file data
C        call mpibc1(nxi(is),1,2,mlog,'ioatmx','nxi')
C        call mpibc1(exi(1,is),nxi(is),4,mlog,'ioatmx','exi')
C        call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'ioatmx','hfc')
C        call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'ioatmx','hfct')
C        call mpibc1(rsmfa(is),1,4,mlog,'ioatmx','rsmfa')
C        call mpibc1(a0,1,4,mlog,'ioatmx','a0')
C        if (nr0 > 0) then
C        call mpibc1(s_rhofa(is)%p,nr0*nsp,4,mlog,'ioatmx','rhofa')
C        call mpibc1(s_rhoca(is)%p,nr0*nsp,4,mlog,'ioatmx','rhoca')
C        call mpibc1(s_v0a(is)%p,nr0*nsp,4,mlog,'ioatmx','v0a')
C        endif
C        i = mpipid(2)
C
CC ...   Defaults
C        if (procid == master) then
C          call strip(spid(is),i1,nch)
C          if (ipr >= 30 .and. rmt0 /= 0)
C     .      write(stdo,400) spid(is)(1:nch),spidr,rmt0,nr0,a0
C  400     format(' ioatmx: expected ',a,',',T27,' read ',a,
C     .           ' with rmt=',f8.4,'  mesh',i6,f7.3)
C        endif
C        if (nr <= 0)   nr = nr0
C        if (a <= 1d-6) a = a0
C        if (z == 0 .and. rmt == 0) then
C          a = 0
C          nr = 0
C        endif
C
CC ...   Sanity checks
C        if (procid == master) then
C          call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
C          call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
C          call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
C          call sanrg(.true.,nr0,nr,nr,msg,'nr')
C        endif
C
CC ...   Pack into sspec
C        s_spec(is)%a = a
C        s_spec(is)%nr = nr
C        s_spec(is)%qc = qc
C        s_spec(is)%nxi = nxi(is)
C        s_spec(is)%exi = exi(:,is)
C        s_spec(is)%chfa = hfc(:,:,is)
C        s_spec(is)%rsmfa = rsmfa(is)
C        s_spec(is)%ctail = ccof
C        s_spec(is)%etail = ceh
C        s_spec(is)%stc = stc
C        if (nr > 0) then
C          call ptr_spec(s_spec,4+1,'rhoc',is,nr*nsp,0,s_rhoca(is)%p)
C        endif
C      enddo
CC     Wait for all proccesses to synchronize
C      i = mpipid(2)
CC     Re-broadcast entire species structure, and arrays used below
C      call bcast_strx(1,xx,xx,xx,xx,xx,xx,s_spec,xx,xx,nspec,0)
C
C      if (procid == master) then
C        call fclose(ifi)
C      endif
C
CC --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
C      ztot = 0d0
C      ctot = 0d0
C      corm = 0d0
C      do  ib = 1, nbas
C        is = s_site(ib)%spec
C        a = s_spec(is)%a
C        nr = s_spec(is)%nr
C        rmt = s_spec(is)%rmt
C        lmxl = s_spec(is)%lmxl
C        z = s_spec(is)%z
C        qc = s_spec(is)%qc
C        lfoc = s_spec(is)%lfoca
C
C        nlml = (lmxl+1)**2
C        call dfratm(s_site,s_spec,1,ib,ib,xx)
C        allocate(v0(nr*nsp),v1(nr*nsp))
C
CC       Core magnetic moment (possible if magnetized core hole)
C        if (nsp == 2 .and. lmxl > -1) then
C          allocate(rwgt(nr))
C          call radwgt(intopt,rmt,a,nr,rwgt)
C          call radsum(nr,nr,1,nsp,rwgt,s_rhoca(is)%p,sum)
C          call radsum(nr,nr,1,1,rwgt,s_rhoca(is)%p,sum1)
C          sum2 = sum - sum1
C          call gtpcor(s_spec,is,kcor,lcor,qcor)
C
C          if (dabs(qcor(2)-(sum1-sum2)) > 0.01d0) then
C            call info5(10,0,0,' (warning) core moment mismatch spec %i:'
C     .        //'  input file=%;6d  atom file=%;6d',
C     .        is,qcor(2),sum1-sum2,0,0)
C          endif
C          corm = corm + qcor(2)
C          deallocate(rwgt)
C        endif
C
C        if (lmxl > -1) then
C          call dpcopy(s_v0a(is)%p,  v0,  1,nr*nsp,1d0)
C          call dpcopy(s_v0a(is)%p,  v1,  1,nr*nsp,1d0)
C          call dpcopy(s_rhoca(is)%p,s_site(ib)%rhoc,1,nr*nsp,1d0)
C
C          if (lfoc == 0) then
C            allocate(rwgt(nr))
C            call radwgt(intopt,rmt,a,nr,rwgt)
C            call radsum(nr,nr,1,nsp,rwgt,s_site(ib)%rhoc,sum)
C            fac = 1d0
C            if(dabs(sum) > 1d-7) fac = qc/sum
C            if (ipr >= 40) write(stdo,787) is,qc,sum,fac
C  787       format(' scale foca=0 core species',i2,': qc,sum,scale=',
C     .        3f12.6,f12.6)
C            call dpcopy(s_site(ib)%rhoc,s_site(ib)%rhoc,1,nr*nsp,fac)
C            deallocate(rwgt)
C          endif
C        endif
C
C        call ptr_site(s_site,4+1,'v0',ib,nr,nsp,v0)
C        call ptr_site(s_site,4+1,'v1',ib,nr,nsp,v1)
C        ztot = ztot+z
C        ctot = ctot+qc
C
CC       end loop over sites
C        deallocate(v0,v1)
C      enddo
C
CC --- Overlap smooth hankels to get smooth interstitial density ---
C      ng = s_lat%ng
C      allocate(cv(ng*ns4))
C      call ovlpfa(s_site,s_lat,s_ham,1,nbas,nxi,n0,exi,hfc,rsmfa,ng,ng,
C     .  s_lat%gv,cv)
C      call gvputf(ng,ns4,s_lat%kv,k1,k2,k3,cv,s_pot%smrho)
C      deallocate(cv)
C
CC ... FFT to real-space mesh
C      call fftz3(s_pot%smrho,n1,n2,n3,k1,k2,k3,ns4,0,1)
C      if (ns4 == 4) then ! Transpose to packed format
C        i = k1*k2*k3
C        call rotspv(251,xx,i,i,i,xx,s_pot%smrho,xx,s_pot%smrho)
C      endif
C
CC ... Add compensating uniform electron density to compensate background
C      call addbkgsm(s_pot%smrho,k1,k2,k3,nsp,qbg,vol,-1d0)
CC ... integrate for mesh density and magnetic moment
C      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,s_pot%smrho,sum1,sum2)
C      smom = 0
C      if (nsp == 2) then  ! Collinear spin polarized case
C        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho,smom,xx)
C        smom(3) = 2*smom(0) - sum1
C        smom(0) = smom(3)
C      endif
C      if (ns4 == 4) then  ! Noncollinear case
C        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho(1,3),smom(1),xx)
C        call mshint(vol,1,n1,n2,n3,k1,k2,k3,s_pot%smrho(1,4),smom(2),xx)
C        smom(1) = 2*smom(1)
C        smom(2) =-2*smom(2)
C        smom(0) = dlength(3,smom(1),1)
C      endif
C
CC --- Set up local densities using rmt from atm file ---
C      call ovlocr(nbas,s_site,s_spec,s_lat,s_ham,n0,nxi,exi,
C     .            hfc,rsmfa,s_rhofa,sqloc,slmom)
CC --- Add compensating uniform electron density to compensate background
C      call adbkql(nbas,nsp,qbg,vol,-1d0,s_spec,s_site)
C      if (abs(qbg) /= 0) call info(10,0,0,' Uniform '//
C     .  'density added to neutralize background, q=%;6,6d',qbg,0)
C
C      do  is = 1, nspec
C        deallocate(s_rhofa(is)%p)
C        deallocate(s_rhoca(is)%p)
C        deallocate(s_v0a(is)%p)
C      enddo
C      deallocate(s_rhofa,s_rhoca,s_v0a)
C
CC --- Print charges ---
C      dq = sum1+sqloc+ctot-ztot+qbg
C      if (nsp == 1 .and. ipr >= 10) then
C        write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
C  895   format(/' Smooth charge on mesh:    ',f16.6
C     .     /    ' Sum of local charges:     ',f16.6
C     .     /    ' Total valence charge:     ',f16.6
C     .     /    ' Sum of core charges:      ',f16.6
C     .     /    ' Sum of nuclear charges:   ',f16.6
C     .     /    ' Homogeneous background:   ',f16.6
C     .     /    ' Deviation from neutrality:',f16.6)
C        write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
C  710   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
C     .    '   bg',f10.6,'  dQ',f10.6)
C      elseif (ipr >= 10) then
C        write(stdo,896) sum1,smom(0),sqloc,slmom(0)
C        if (ns4 == 4) write(stdo,898) smom(1:3), slmom(1:3)
C        write(stdo,897)
C     .    sum1+sqloc,smom(0)+slmom(0),ctot,corm,-ztot,qbg,dq
C  896   format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6/
C     .          ' Sum of local charges:     ',f16.6,4x,'moments',f11.6)
C  897   format( ' Total valence charge:     ',f16.6,4x,'moment', f12.6,
C     .     /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6,
C     .     /    ' Sum of nuclear charges:   ',f16.6
C     .     /    ' Homogeneous background:   ',f16.6
C     .     /    ' Deviation from neutrality:',f16.6)
C        write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom(0)+slmom(0)
C  898   format(' Mesh magnetization (xyz):',6x,3f11.6/
C     .         ' Local magnetization (xyz):',5x,3f11.6)
C  711   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
C     .    '   bg',f11.6,' mm',f11.6)
C      endif
C
C      if (dabs(dq) > 1d-4 .and. ipr > 1)
C     .  call awrit1(' ioatmx (warning) overlapped'
C     .  //' density not neutral'//', dq= %d',' ',80,stdo,dq)

      endif

C --- Output ---
      if (ifi < 0)  then

        call info0(30,0,0,' ioatmx: saving atomic data by site')
        jfi = -ifi

        call dpzero(hfct,size(hfct))
        do  ib = ib1, ib2
          is = s_site(ib)%spec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          nlml = (s_spec(is)%lmxl+1)**2
          z = s_spec(is)%z
          qc = s_spec(is)%qc
          lfoc = s_spec(is)%lfoca

          allocate(rhov(nr,4))
          do  isp = 1, nsp*nspc
            rhov(:,isp) = s_site(ib)%rho1(:,1+(isp-1)*nlml)/y0
          enddo
          lfail =
     .     iofa(63-2,s_spec(is)%name,n0,s_spec(is)%nxi,s_spec(is)%exi,s_spec(is)%chfa,hfct,s_spec(is)%rsmfa,
     .      s_spec(is)%z,s_spec(is)%rmt,s_spec(is)%a,s_spec(is)%nr,nsp,s_spec(is)%qc,s_spec(is)%ctail,
     .      s_spec(is)%etail,s_spec(is)%stc,rhov,s_site(ib)%rhoc,s_site(ib)%v0,ifi) < 0
          if (lfail)  call rx('failed to write atm data')
          deallocate(rhov)

        enddo
        call rx0('atomic data saved by site')

      endif
      end
