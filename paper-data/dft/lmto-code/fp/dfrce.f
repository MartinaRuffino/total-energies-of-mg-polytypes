      subroutine dfrce(s_site,s_spec,s_lat,s_ctrl,k1,k2,k3,nvl,
     .  s_rhat,s_rnew,elind,qmom,smrho,smrout,dfh)
C- Correction to force theorem, Harris functional
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvdf4 pvdf2 smvxcm smcorm smvxc4 rhomom pvdf1
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl rg lfoca rfoca qc z ctail etail stc lmxb p pz
Ci                 rmt a nr lmxa nxi exi chfa rsmfa coreh coreq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvdf4 corprm pvdf2 smvxcm smcorm smvxc4 rhomom pvdf1
Cio                gtpcor
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc ng vol alat plat qlat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv kv cy
Cio    Passed to:  pvdf4 pvdf2 smvxcm smcorm smvxc2 vxcnlm pvdf1
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfrce
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_rnew
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rhomom
Ci Inputs
Ci   k1..3 :dimensions smrho
Ci   nvl   :sum of local (lmxl+1)**2, lmxl = density l-cutoff
Ci   orhoat:vector of offsets containing site density
Ci   s_rnew:pointer to local densities
Ci   elind :Lindhard parameter, used for Lindhard screening
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   smrho :smooth density on uniform mesh
Ci   smrho :smooth (input) density that generated the hamiltonian
Cio  smrout:smooth (output) density that the hamiltonian generated
Co Outputs
Co   dfh   :correction to the HF force
Cl Local variables
Cl    job  :describes which ansatz for charge shift is used for correction
Cl         :<=0  do not calculate correction to force
Cl         :  1  shift in free-atom density
Cl         :  2  shift in core+nuclear density
Cl         :+10  to screen the rigid shift by the Lindhard function
Cr Remarks
Cu Updates
Cu   31 Jul 17 zeros out correction for points when rho+ or rho- is small
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   18 Dec 03 adapted to modified smvxc
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   17 Sep 01 Adapted for local orbitals
Cu   21 Jun 00 spin polarized
Cu   18 Jun 98 adapted from nfp dfrce.f
Cu   16 Jun 98 MvS parallelized for SGI
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPE
C      include "mpef.h"
C#endif
      integer procid, master, numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      logical mlog,cmdopt
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
      integer nvl,k1,k2,k3
      double precision dfh(3,*),qmom(*)
      double complex smrho(k1,k2,k3,*),smrout(k1,k2,k3,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ctrl)::  s_ctrl
      type(str_rhat):: s_rhat(*),s_rnew(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iv(:)
      real(8), allocatable :: cs(:),sn(:),qmout(:),yl(:),g(:),g2(:)
      complex(8), allocatable :: dvxc(:),smro(:),vxcp(:),vxcm(:)
      complex(8), allocatable :: wk1(:),wk2(:),wk3(:)
      complex(8), allocatable :: ceps(:),cnomi(:),cvin(:),cdvx(:)
      complex(8), allocatable :: cdn(:)
      complex(8), allocatable :: cdn0(:)
      complex(8), allocatable :: cdv(:)
C      real(8), allocatable :: cs(:)
C      real(8), allocatable :: sn(:)
C Local variables
      integer nbas,job,n1,n2,n3,ng,iprint,nsp,nglob,ib,is,lmxl,iv0,nlm,
     .  ip,i,ngabc(3),ltop,nspec,nlmtop,nn,lgunit,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C      integer ocdv,ocdn,ocdn0,ocs,osn
      double precision vol,plat(3,3),qlat(3,3),alat,vsum,pi,tpiba,elind,
     .  fes1(3),fes2(3),fxc(3),c,avgdf(3)
      integer npmx
      parameter (npmx=32)
C     integer oicdn(npmx),oicdn0(npmx),oicdv(npmx),oics(npmx),oisn(npmx)
C      type(dc_wk_vec) :: s_cdn(npmx),s_cdn0(npmx),s_cdv(npmx)
C      type(dp_wk_vec) :: s_cs(npmx),s_sn(npmx)

      integer, dimension(:), allocatable :: bproc, iiv0
      integer pid,jb

C     double precision dqsmo,dval
      character*40 strn
C ... Heap
      integer w(1)
      common /w/ w

      job = s_ctrl%lfrce
      if (job <= 0) return
      call tcn('dfrce')


      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs > 1) then
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      if (mlog) then
        do  pid = 0, numprocs-1
          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,MPI_COMM_WORLD,ierr)
        enddo
      endif
      endif

C --- Setup ---
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      c = 1000
      stdo = lgunit(1)
      nsp  = nglob('nsp')
      nbas = nglob('nbas')
      nspec= nglob('nspec')
      nn   = k1*k2*k3

C ... Arrays needed for pvdf1
      allocate(ceps(ng),cnomi(ng),cvin(ng),cdvx(ng*nsp))

C ... Set up for vectorized Y_lm and gaussians
      ltop = 0
      do  is = 1, nspec
        lmxl = s_spec(is)%lmxl
        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2
      allocate(yl(ng*nlmtop),g2(ng),g(ng*3))
      call suylg(ltop,alat,ng,s_lat%gv,g,g2,yl)
      deallocate(g)
      allocate(iv(ng*3))
      call suphs0(plat,ng,s_lat%gv,iv)

C --- Make ves(rhoin,q) ---
      allocate(smro(nn),cs(ng),sn(ng))
      call dpcopy(smrho,smro,1,2*nn,1d0)
      if (nsp == 2) call dpadd(smro,smrho(1,1,1,2),1,2*nn,1d0)
      call fftz3(smro,n1,n2,n3,k1,k2,k3,1,0,-1)
      call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smro,cvin)
      call pvdf4(s_site,s_spec,s_lat,qmom,ng,
     .  g2,yl,cs,sn,iv,qlat,cvin)
C     call gvputf(ng,1,s_lat%kv,k1,k2,k3,cvin,smro)
C     call zprm3('ves(input rho)',0,smro,k1,k2,k3)
      deallocate(smro,cs,sn)

C --- Make dVxc(in)/dn ---
      allocate(dvxc(nn*nsp),smro(nn*nsp))
      allocate(vxcp(nn*nsp),vxcm(nn*nsp))
      allocate(wk1(nn*nsp),wk2(nn*nsp),wk3(nn*nsp))
      call dpcopy(smrho,smro,1,2*nn*nsp,1d0)
      call pvdf2(nbas,nsp,s_site,s_spec,s_lat,n1,n2,n3,k1,k2,k3,
     .  smro,vxcp,vxcm,wk1,wk2,wk3,dvxc)
      deallocate(vxcp,vxcm,wk1,wk2,wk3)

C --- cdvx = FFT ((n0_out-n0_in) dVxc/dn) ---
C     Use total n0_out-n0_in but keep vxc spin polarized
      call dpzero(smro,2*nn)
      do  i = 1, nsp
        call dpadd(smro,smrout(1,1,1,i),1,2*nn,1d0)
        call dpadd(smro,smrho(1,1,1,i),1,2*nn,-1d0)
      enddo
C     call zprm3('drho',0,smro,k1,k2,k3)
      call pvdf3(n1,n2,n3,k1,k2,k3,nsp,smro,dvxc)
      call fftz3(dvxc,n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng,nsp,s_lat%kv,k1,k2,k3,dvxc,cdvx)

C --- Cnomi = (n0_out(q) - n0_in(q)) ---
      call fftz3(smro,n1,n2,n3,k1,k2,k3,1,0,-1)
c     call zprm3('FFT smrout-smrin',smro,k1,k2,k3*nsp)
      call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smro,cnomi)
C     dqsmo = vol*dval(cnomi,1)

C ... Debugging slot smrho(out) for out-in
C      print *, '*** debugging ... subs smrout for out-in'
C      call dpcopy(smrout,smro,1,2*nn,1d0)
C      call fftz3(smro,n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smro,cnomi)
C      call zprm3('rho-out(q)',smro,k1,k2,k3)

      deallocate(dvxc,smro)

C --- Multipole moments of the output density ---
      allocate(qmout(nvl))
      call pshpr(1)
      call rhomom(0,1,nbas,s_site,s_spec,s_rnew,qmout,vsum)
      call poppr
      call dpadd(qmout,qmom,1,nvl,-1d0)

C --- Lindhard dielectric function ---
      if (job > 10) then
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
C        call lindxx(123,n1,n2,n3,k1,k2,k3,ng,s_lat%kv,ceps,s_lat%gv,
C     .    tpiba,elind,w,w,w,w)
        call lindsc(3,ng,s_lat%gv,tpiba,elind,ceps)
      endif

C --- For each site, get correction to force ---
C     if (numprocs > 1) then
      if (iprint() >= 30) then
        strn = 'shift in free-atom density'
        if (job == 11) strn = 'screened shift in free-atom density'
        if (job == 12) strn = 'screened shift in core+nuclear density'
        write(stdo,201) strn
      endif
C     endif
  201 format(/' Harris correction to forces: ',a/
     .  '  ib',9x,'delta-n dVes',13x,'delta-n dVxc',15x,'total')
C  201   format(/' Harris correction to forces:'/
C     .    '  ib',11x,'dn0 dVes',15x,'dnloc dVes',15x,
C     .    'dn0 dVxc',15x,'total')

C ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
      if (numprocs > 1) then
      iv0 = 0
      allocate(iiv0(1:nbas), stat=ierr)
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
      enddo
      endif

C      mp = 1
C      do  ip = 1, min(nbas,mp)
C        allocate(s_cdn(ip)%p(ng))
C        allocate(s_cdn0(ip)%p(ng*nsp))
C        allocate(s_cdv(ip)%p(ng))
C        allocate(s_cs(ip)%p(ng))
C        allocate(s_sn(ip)%p(ng))
C      enddo

C ... Estimate shift in density for each site, and corresponding force
      ip = 1
      iv0 = 0
      allocate(cdn(ng),cdn0(ng*nsp),cdv(ng))
      allocate(cs(ng),sn(ng))
      allocate (bproc(0:numprocs))
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_DFRCE,procid,"dfrce")
C#endif
      call dstrbp(nbas,numprocs,1,bproc(0))
      else
      bproc(0:1) = [1,nbas+1]
      end if
      do  ib = bproc(procid), bproc(procid+1)-1
        if (numprocs > 1) then
        if (mlog .and. ib == bproc(procid)) then
          call gettime(datim)
          call awrit4(' dfrce '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting atoms %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        iv0 = iiv0(ib)
        endif
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        if (lmxl == -1) cycle

        nlm = (lmxl+1)**2
C        ocs = oics(ip)
C        osn = oisn(ip)
C        ocdn = oicdn(ip)
C        ocdn0 = oicdn0(ip)
C        ocdv = oicdv(ip)

        call pvdf1(job,s_site,s_spec,s_lat,nsp,ib,iv0,qmom,qmout,
     .    ng,s_lat%gv,g2,yl,cs,sn,iv,qlat,cnomi,ceps,cdn0,cdn,cdv,
     .    cdvx,cvin,s_rhat(ib),fes1,fes2,fxc)

C        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cdn,smpot)
C        call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
C        call zprm3('dn',smpot,k1,k2,k3)
C
        do  i = 1, 3
          dfh(i,ib) = -(fes1(i) + fes2(i) + fxc(i))
        enddo

        if (numprocs == 1) then
C        if (iprint() >= 30)
C     .    write(stdo,200) ib,(c*(fes1(m)+fes2(m)),m=1,3),
C     .    (c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
C  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)
        call info5(30,0,0,'%,4i%3;8,2D %3;8,2D %3;8,2D',
     .    ib,c*(fes1(1:3)+fes2(1:3)),c*fxc(1:3),c*dfh(1:3,ib),0)
        endif

        iv0 = iv0+nlm
      enddo
      deallocate(cdn,cdn0,cdv)
      deallocate(cs,sn)
      deallocate(qmout,yl,g2,iv)

      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_DFRCE,procid,"dfrce")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_BCAST,procid,"broadcast")
C#endif
      do  pid = 0, numprocs-1
        ib = bproc(pid)
        jb = bproc(pid+1) - ib
        call MPI_BCAST(dfh(1,ib),3*jb,mpi_real8,pid,
     .                 MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit6(' dfrce '//datim//' Process %i of %i on '
     .          //shortname(procid)(1:namelen(procid))//
     .          ' bcast dfh(1-3,%i-%i) %i d.p. numbers'//
     .          ' from process %i on '
     .          //shortname(pid)(1:namelen(pid)),' ',
     .          256,lgunit(3),procid,numprocs,bproc(pid),bproc(pid+1)-1,
     .          3*jb,pid)
        endif
      enddo
      deallocate (bproc, stat=ierr)
      deallocate (iiv0, stat=ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BCAST,procid,"broadcast")
C#endif
      endif
      call dpzero(avgdf,3)
      do  ib = 1, nbas
        do  i = 1, 3
          avgdf(i) = avgdf(i)+dfh(i,ib)/nbas
        enddo
      enddo

C ... MPI printout
      if (numprocs > 1) then
      if (iprint() >= 30) then
        strn = 'shift in free-atom density'
        if (job == 11) strn = 'screened shift in free-atom density'
        if (job == 12) strn = 'screened shift in core+nuclear density'
        write(stdo,201) strn
        do  ib = 1, nbas
          write (stdo,2) ib,(c*dfh(i,ib),i = 1,3)
    2     format(i4,50x,3F8.2)
        enddo
      endif
      endif

C ... Shift all forces to make avg correction zero
      do  ib = 1, nbas
        do  i = 1, 3
          dfh(i,ib) = dfh(i,ib) - avgdf(i)
        enddo
      enddo
C      if (iprint() >= 30) write (stdo,3) (c*avgdf(m),m = 1,3)
C    3 format(' shift forces to make zero average correction:',8x,3F8.2)
      call info2(30,0,0,' shift forces to make zero average correction:%8f%3;8,2D',
     .  c*avgdf,0)
      deallocate(ceps,cnomi,cvin,cdvx)
C      do  ip = 1, min(nbas,mp)
C        deallocate(s_cdn(ip)%p)
C        deallocate(s_cdn0(ip)%p)
C        deallocate(s_cdv(ip)%p)
C        deallocate(s_cs(ip)%p)
C        deallocate(s_sn(ip)%p)
C      enddo

      call tcx('dfrce')
      end

      subroutine pvdf1(job,s_site,s_spec,s_lat,nsp,ib,iv0,qmom,qmout,ng,
     .  gv,g2,yl,cs,sn,iv,qlat,cnomin,ceps,cdn0,cdn,cdv,cdvxc,cvin,
     .  s_rhat,fes1,fes2,fxc)
C- Estimate shift in local density for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: z p pz lmxa a nr rmt lmxl nxi exi chfa rsmfa rg coreh
Ci                coreq lfoca rfoca qc ctail etail stc orhoc lmxb
Ci     Stored:    *
Ci     Passed to: gtpcor corprm
Ci   s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: nabc alat vol
Ci     Stored:    *
Ci     Passed to: *
Ci   ng,gv
Ci   job: 1  shift in free-atom density
Ci        2  shift in core+nuclear density
Ci      +10  to screen the rigid shift by the response function
Ci   ib      which site is being shifted
Ci   iv0     offset to qmom
Ci   qmom,qmout moments of input and output densities
Ci   cnomin  difference betw. smoothed output and input density n0
Ci   cvin    electrostatic potential of input density Ves[n0~_in]
Ci   ceps    response function
Ci   cdvxc   dVxc/dn (nout-nin)
Co Outputs
Co   cdn0:   Job 1:  shift in valence part of the free atom density
Co           Job 12: shift in atom density (1/eps - 1)
Co   cdn:    Job 1:  dn^(u) where dn is the unscreened shift in
Co           in the free-atom density.
Co           Job 12: dn^(u) 1/eps where dn is unscreened shift in
Co           the charge density.  Local density approximated
Co   NB:     In all cases, the local part of density is approximated
Co           by a gaussian of the equivalent multipole moment.
Co   cdv:    shift in the electrostatic potential
Co   fes1,fes2,fxc
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ng,nsp,iv0,ib,job,iv(ng,3)
      double precision qmom(*),qmout(*),
     .  gv(ng,3),tau(3),fes1(3),fes2(3),fxc(3),g2(ng),yl(ng,1),
     .  cs(ng),sn(ng),qlat(3,3)
      double complex cdn0(ng,nsp),cdn(ng),cdv(ng),ceps(ng),
     .  cnomin(ng),cdvxc(ng,nsp),cvin(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_rhat):: s_rhat
C ... Local parameters
      integer ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jv0,jb,nbas,nglob,js,ll,n0,nrmx,intopt
      parameter (nlmx=81, nrmx=5001, n0=10)
      integer lmxa,nr,nxi,ie,ixi,job0,kcor,lcor,lfoc,i,
     .  ngabc(3),n1,n2,n3,nlml
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision pi,alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .  vol,y0,z,v(3),df(0:20),feso(3),qcor(2),gpot0(nlmx,3),fesdn(3),
     .  fesgg(3),pnu(n0,2),pnz(n0,2),a,rmt,qloc,rsmfa,exi(n0),hfc(n0,2),
     .  qfat,gam,qall,qc,qval,qg,e,aa,q0(3),sum
      double precision rwgt(nrmx),cc,gamf,cfoc,cvol
C     parameter (k0=3)
C     double complex gkl(0:k0,nlmx)
      double complex tpia,cxx,phase,gc0,xc0,cof(nlmx)
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf1')
      ngabc = s_lat%nabc
      alat = s_lat%alat
      vol = s_lat%vol
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      nbas = nglob('nbas')
C     tpiba = 2*pi/alat
      tpia = 2*pi*dcmplx(0d0,-1d0)/alat
      job0 = mod(job,10)

      call dpzero(cdn,2*ng)
      call dpzero(cdn0,2*ng*nsp)
      cdv(1) = 0d0
      call dpzero(fes1,3)
      call dpzero(fxc,3)
      call dpzero(fesdn,3)
      call dpzero(gpot0,nlmx*3)
      is = s_site(ib)%spec
      tau = s_site(ib)%pos
      call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)

C --- Unscreened rigid charge density shift, job 1, in cdn0 ---
      if (job0 == 1) then
C        is = s_site(ib)%spec
C        tau = s_site(ib)%pos
        z = s_spec(is)%z
        pnu = s_spec(is)%p
        pnz = s_spec(is)%pz
        lmxa = s_spec(is)%lmxa
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        lmxl = s_spec(is)%lmxl
        nxi = s_spec(is)%nxi
        exi = s_spec(is)%exi
        hfc = s_spec(is)%chfa
        rsmfa = s_spec(is)%rsmfa
        call gtpcor(s_spec,is,kcor,lcor,qcor)
        if (nr > nrmx) call rx('dfrce: nr gt nrmx')
        intopt = 10*nglob('lrquad')
        call radwgt(intopt,rmt,a,nr,rwgt)
        nlml = (lmxl+1)**2
        call radsum(nr,nr,nlml,nsp,rwgt,s_rhat%rho1,qloc)
        call radsum(nr,nr,nlml,nsp,rwgt,s_rhat%rho2,sum)
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
        qg = qval+qsc-qfat-qloc
C   ... Shift in free atom density
        do  i = 1, nsp
        do  ixi = 1, nxi
        e = exi(ixi)
C       do  14  ig = 1, ng
C         v(1) = gv(ig,1)*tpiba
C         v(2) = gv(ig,2)*tpiba
C         v(3) = gv(ig,3)*tpiba
C         v2 = v(1)**2+v(2)**2+v(3)**2
C         aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
C         scalp = -alat*(tau(1)*v(1)+tau(2)*v(2)+tau(3)*v(3))
C         phase = dcmplx(dcos(scalp),dsin(scalp))
C         cdn0(ig) = cdn0(ig) + hfc(ixi,i)*aa*phase*y0/vol
C     ... Vectorized version
          cc = -4d0*pi*hfc(ixi,i)*y0/vol
          do  ig = 1, ng
            aa = cc*dexp(gam*(e-g2(ig)))/(e-g2(ig))
            cdn0(ig,i) = cdn0(ig,i) + aa*dcmplx(cs(ig),sn(ig))
          enddo
        enddo
        enddo

C   ... Add gaussian to conserve local charge
C     ... Add gaussian to conserve local charge.  If density corresponds
C         to the free-atom density, qfat+qloc = qval+qsc; then qg=0

C       do  16  ig = 1, ng
C         v(1) = gv(ig,1)*tpiba
C         v(2) = gv(ig,2)*tpiba
C         v(3) = gv(ig,3)*tpiba
C         v2 = v(1)**2+v(2)**2+v(3)**2
C         scalp = -alat*(tau(1)*v(1)+tau(2)*v(2)+tau(3)*v(3))
C         phase = dcmplx(dcos(scalp),dsin(scalp))
C         cdn0(ig) = cdn0(ig) + qg*phase*dexp(-gam*v2)/vol
C  16   continue
C   ... Vectorized version
        cc = qg/vol/nsp
        do  i = 1, nsp
        do  ig = 1, ng
        cdn0(ig,i)=cdn0(ig,i)+cc*dcmplx(cs(ig),sn(ig))*dexp(-gam*g2(ig))
        enddo
        enddo
      endif

C --- Coefficients defining local valence + core density ---
      is = s_site(ib)%spec
      tau = s_site(ib)%pos
      lmxl = s_spec(is)%lmxl
      rg = s_spec(is)%rg
      call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
      nlm = (lmxl+1)**2
      if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
      ilm = 0
      cxx = dcmplx(0d0,1d0)
      do  l = 0, lmxl
        cxx = cxx*dcmplx(0d0,-1d0)
        do  m = -l, l
          ilm = ilm+1
          cof(ilm) = cxx*qmom(ilm+iv0)*4*pi/df(2*l+1)
        enddo
      enddo
C     cof(1) = cof(1) + 4*pi*y0*(qcorg+qsc-z)
      cof(1) = cof(1) + 4*pi*y0*(qcorg-z)

C --- Shift in n0, ves~ for list of G vectors ---
      gam = 0.25d0*rg*rg
      gamf = 0.25d0*rfoc*rfoc
      cfoc = -4d0*pi*y0*cofh/vol
      cvol = 1d0/vol
      do  ig = 2, ng

        v(1) = gv(ig,1)
        v(2) = gv(ig,2)
        v(3) = gv(ig,3)

C   ... Accumulate unscreened smoothed core+nuclear density

CC       Old, serial version
C        call gkl_ft(v,rg,0d0,tau,alat,kmax,nlm,k0,cy,gkl)
C        do  32  ilm = 1, nlm
C   32   cdn(ig) = cdn(ig) + cof(ilm)*gkl(0,ilm)/vol
C
CC       This part for (grad g) ves(in)
C        do  33  k = 1, 3
C        cxx = tpia*v(k)*cvin(ig)
C        do  33  ilm = 1, nlm
C   33   gpot0(ilm,k) = gpot0(ilm,k) + dconjg(cxx)*gkl(0,ilm)

C   ... Vectorized version (absorb -i**l)
        phase = dcmplx(cs(ig),sn(ig))
        gc0 = phase*dexp(-gam*g2(ig))*cvol
        xc0 = dcmplx(0d0,1d0)*dconjg(tpia*cvin(ig))*gc0*vol
        ilm = 0
        do  l = 0, lmxl
          xc0 = xc0*dcmplx(0d0,-1d0)
          do  m = -l, l
            ilm = ilm+1
            cdn(ig) = cdn(ig) + yl(ig,ilm)*cof(ilm)*gc0
            gpot0(ilm,1) = gpot0(ilm,1) + yl(ig,ilm)*gv(ig,1)*xc0
            gpot0(ilm,2) = gpot0(ilm,2) + yl(ig,ilm)*gv(ig,2)*xc0
            gpot0(ilm,3) = gpot0(ilm,3) + yl(ig,ilm)*gv(ig,3)*xc0
          enddo
        enddo

C   ... Accumulate unscreened foca density
        aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
        cdn(ig) = cdn(ig) + cfoc*aa*phase
C       A slow, unvectorized version
C       call hkl_ft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy, gkl)
C       cdn(ig) = cdn(ig) + cofh*gkl(0,1)/vol

C   ... Make the screened shift in input density n0~
C       Job 1: cdn0 = (valence part of) cdn^u ; cdn = cdn^u
        if (job0 == 1) then
          cdn(ig) = cdn(ig) + (cdn0(ig,1) + cdn0(ig,nsp))/(3-nsp)
          if (job > 10) cdn(ig) = cdn(ig) / ceps(ig)
C       Job 12: cdn0 = cdn^u (1/eps - 1); cdn = cdn^s = cdn^u / eps
        elseif (job == 12) then
          do  i = 1, nsp
            cdn0(ig,i) = cdn(ig) * (1/ceps(ig)-1) / nsp
          enddo
          cdn(ig) = cdn(ig) / ceps(ig)
        else
          call rxi('dfrce: nonsensical job',job)
        endif

C   ... Electrostatic potential shift = 1/eps dv [n0~]
C       g2 = tpiba*tpiba*(gv(ig,1)**2+gv(ig,2)**2+gv(ig,3)**2)
        cdv(ig)  = cdn(ig) * (8*pi/g2(ig))

C       fes1 = (n0_out - n0_in) d ves[n0~]
C       fxc  = dVxc/dn (nout-nin) d n0~
        do  k = 1, 3
          fes1(k) = fes1(k) + dconjg(cnomin(ig)) * tpia*v(k)*cdv(ig)
          do  i = 1, nsp
          fxc(k)  = fxc(k)  + dconjg(cdvxc(ig,i)) * tpia*v(k)*cdn0(ig,i)
          enddo
C         fesdn(k)= fesdn(k) + dconjg(cvin(ig))  * tpia*v(k)*cdn0(ig,i)
        enddo

      enddo

      do  k = 1, 3
C       fesdn(k)  = fesdn(k)*vol
        fxc(k)  = fxc(k)*vol
        fes1(k) = fes1(k)*vol
      enddo

C --- Integral of grad g (output-input local charge) ves~ ---
      call dpzero(fesgg,3)
      do  k = 1, 3
      do  ilm = 1, nlm
        l = ll(ilm)
        gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)

C       fesgg(k) = fesgg(k) + qmom(iv0+ilm)*gpot0(ilm,k)
        fesgg(k) = fesgg(k) + qmout(iv0+ilm)*gpot0(ilm,k)
      enddo
      enddo

C      print 339, 'n0(out-in) * g dves ',fes1
C      print 339, 'd(g) qmom(out-in) ves[n0~]',fesgg
C      print 339, 'n0~(out-in) * dvxc   ',fxc
C  339 format(a,6p,3f8.2)

C --- Integral of dves~ (output-input local charge) for all sites ---
      call dpzero(fes2,3)
      call dpzero(feso,3)
      jv0 = 0
      do  jb = 1, nbas
        js = s_site(jb)%spec
        tau = s_site(jb)%pos
        lmxl = s_spec(js)%lmxl
        rg = s_spec(js)%rg
        nlm = (lmxl+1)**2

C ... For this jb, mesh density for all G vectors
        if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
        call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
        call dpzero(gpot0,nlmx*3)
        gam = 0.25d0*rg*rg
        do  ig = 2, ng

C          v(1) = gv(ig,1)
C          v(2) = gv(ig,2)
C          v(3) = gv(ig,3)
C          call gkl_ft(v,rg,0d0,tau,alat,0,nlm,k0,cy,gkl)
C          do  55  k = 1, 3
C          cxx = tpia*v(k)*cdv(ig)
C          do  55  ilm = 1, nlm
C   55     gpot0(ilm,k) = gpot0(ilm,k) + dble(dconjg(cxx)*gkl(0,ilm))

C ... This is the vectorized version
          aa = dexp(-gam*g2(ig))
          gc0 = dcmplx(0d0,1d0)*aa*
     .          dconjg(tpia*cdv(ig))*dcmplx(cs(ig),sn(ig))
          ilm = 0
          do  l = 0, lmxl
            gc0 = gc0*dcmplx(0d0,-1d0)
            do  m = -l, l
              ilm = ilm+1
              gpot0(ilm,1) = gpot0(ilm,1)+dble(gc0)*yl(ig,ilm)*gv(ig,1)
              gpot0(ilm,2) = gpot0(ilm,2)+dble(gc0)*yl(ig,ilm)*gv(ig,2)
              gpot0(ilm,3) = gpot0(ilm,3)+dble(gc0)*yl(ig,ilm)*gv(ig,3)
            enddo
          enddo

C          print 357, ig, cxx, dconjg(cxx)*gkl(0,1)
C  357     format(i4,1p,4e18.8)

        enddo

C   ... Multiply factors into gpot0, accumulate force
        ilm = 0
        do  l = 0, lmxl
          do  m = -l, l
            ilm = ilm+1
            do  k = 1, 3
              gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)
              feso(k) = feso(k) + qmom(jv0+ilm)*gpot0(ilm,k)
              fes2(k) = fes2(k) + qmout(jv0+ilm)*gpot0(ilm,k)
            enddo
          enddo
        enddo

        jv0 = jv0+nlm
      enddo

C      print 339, 'qmom(in) dv         ',feso
C      print 339, 'qmom(out-in) dv     ',fes2

      call dpadd(fes2,fesgg,1,3,-1d0)
      call tcx('pvdf1')

      end

      subroutine pvdf2(nbas,nsp,s_site,s_spec,s_lat,
     .  n1,n2,n3,k1,k2,k3,smrho,vxcp,vxcm,wk1,wk2,wk3,dvxc)
C- Makes derivative of smoothed xc potential wrt density.
C ----------------------------------------------------------------------
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: smvxcm smcorm smvxc4
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lfoca rfoca qc z ctail etail stc orhoc lmxb p pz rmt
Ci                rg
Ci     Stored:    *
Ci     Passed to: smvxcm smcorm corprm smvxc4
Ci   s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: nabc ng gv kv vol alat ocy
Ci     Stored:    *
Ci     Passed to: smvxcm smcorm smvxc2 vxcnlm
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,n1,n2,n3,k1,k2,k3
      double complex vxcp(k1,k2,k3,nsp),vxcm(k1,k2,k3,nsp),
     .               dvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp),
     .               wk1(k1,k2,k3,nsp),wk2(k1,k2,k3,nsp),
     .               wk3(k1,k2,k3,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      logical,allocatable :: lnzero(:,:,:,:)
C ... Local parameters
      integer i1,i2,i3,i,nn
      double precision fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr,
     .  rmusm(2),rvmusm(2),rvepsm(2),repsm(2),repsmx(2),repsmc(2),
     .  fcexc0(2),fcex0(2),fcec0(2),fcvxc0(2),bscali
      real(8), parameter :: tol = 1d-12

      fac = dmach(1)**(1d0/3d0)
      alfa = 2d0/3d0
      nn = k1*k2*k3
      call pshpr(1)
      bscali = 1  ! Energy and potential not yet coupled when XC is scaled
      allocate(lnzero(k1,k2,k3,nsp))

C ... Add fac (rho+ + rho-)/2 into rho+, rho- for spin pol case,
C     Add fac * rho into rho if not spin polarized
      if (nsp == 1) then
        call dpcopy(smrho,smrho,1,nn*2,1d0+fac)
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) + rrho*fac/2
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) + rrho*fac/2
        enddo
        enddo
        enddo
      endif

C ... vxcp = vxc (smrho+drho)
      call dpzero(vxcp, nn*2*nsp)
      call dpzero(wk1, nn*2*nsp)
      call dpzero(wk2, nn*2*nsp)
      call dpzero(wk3, nn*2)
      call smvxcm(s_site,s_spec,s_lat,nbas,nsp,bscali,0,k1,k2,k3,smrho,
     .  vxcp,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... Replace fac*rho with -fac*rho
      if (nsp == 1) then
        call dpcopy(smrho,smrho,1,nn*2,(1d0-fac)/(1d0+fac))
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = (smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2))/(1d0+fac)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) - rrho*fac
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) - rrho*fac
        enddo
        enddo
        enddo
      endif

C ... vxcm = vxc (smrho-drho)
      call dpzero(vxcm, nn*2*nsp)
      call smvxcm(s_site,s_spec,s_lat,nbas,nsp,bscali,0,k1,k2,k3,smrho,
     .  vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... Restore rho+, rho-
      if (nsp == 1) then
        call dpcopy(smrho,smrho,1,nn*2,1/(1d0-fac))
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = (smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2))/(1d0-fac)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) + rrho*fac/2
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) + rrho*fac/2
        enddo
        enddo
        enddo
      endif

C ... Overwrite vxcp with df/drho
      do  i = 1, nsp

CML        do  i1 = 1, nn
CML          rrho = (smrho(i1,1,1,1)+smrho(i1,1,1,nsp))/(3-nsp)
CML          if (rrho > 0) then
CML            f1 = vxcm(i1,1,1,i)*(rrho*(1-fac))**alfa
CML            f2 = vxcp(i1,1,1,i)*(rrho*(1+fac))**alfa
CML            dfdr = (f2-f1)/(2d0*fac*rrho)
CML            vxcp(i1,1,1,i) = dfdr
CML          else
CML            vxcp(i1,1,1,i) = 0
CML          endif
CML        enddo

        do i3=1,k3
          do i2=1,k2
            do i1=1,k1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho > tol) then
                f1 = vxcm(i1,i2,i3,i)*(rrho*(1-fac))**alfa
                f2 = vxcp(i1,i2,i3,i)*(rrho*(1+fac))**alfa
                dfdr = (f2-f1)/(2d0*fac*rrho)
                lnzero(i1,i2,i3,i) = .true.
                vxcp(i1,i2,i3,i) = dfdr
              else
                lnzero(i1,i2,i3,i) = .false.
                vxcp(i1,i2,i3,i) = 0
              endif
            enddo
          enddo
        enddo

      enddo

C ... vxcm = vxc (smrho)
      call dpzero(vxcm, nn*2*nsp)
      call smvxcm(s_site,s_spec,s_lat,nbas,nsp,bscali,0,k1,k2,k3,smrho,
     .  vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... dvxc/drho into dvxc
      do  i = 1, nsp

CML        do  i1 = 1, nn
CML          rrho = (smrho(i1,1,1,1)+smrho(i1,1,1,nsp))/(3-nsp)
CML          if (rrho > 0) then
CML            f = vxcm(i1,1,1,i) * rrho**alfa
CML            dvdr = (vxcp(i1,1,1,i) - alfa*f/rrho) / rrho**alfa
CML            dvxc(i1,1,1,i) = dvdr
CML          else
CML            dvxc(i1,1,1,i) = 0
CML          endif
CML        enddo

        do  i3 = 1, k3
          do i2 = 1, k2
            do i1 = 1, k1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho > tol .and. lnzero(i1,i2,i3,i)) then
                f = vxcm(i1,i2,i3,i) * rrho**alfa
                dvdr = (vxcp(i1,i2,i3,i) - alfa*f/rrho) / rrho**alfa
                dvxc(i1,i2,i3,i) = dvdr
              else
                dvxc(i1,i2,i3,i) = 0
              endif
            enddo
          enddo
        enddo

      enddo

C     call zprm3('d vxc / dn',0,dvxc,k1,k2,k3*nsp)
      call poppr
      deallocate(lnzero)

      end

      subroutine pvdf3(n1,n2,n3,k1,k2,k3,nsp,deln0,dvxc)
C- Overwrites dvxc with (nout-nin)*dvxc
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,nsp
      double complex deln0(k1,k2,k3),dvxc(k1,k2,k3,nsp)
C ... Local parameters
      integer i1,i2,i3,i

      do  i = 1, nsp
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              dvxc(i1,i2,i3,i) = dvxc(i1,i2,i3,i)*deln0(i1,i2,i3)
            enddo
          enddo
        enddo
      enddo

C     call zprm3('dvxc/dn * (nout-nin)',0,dvxc,k1,k2,k3*nsp)

      end

      subroutine pvdf4(s_site,s_spec,s_lat,qmom,ng,g2,yl,cs,sn,iv,qlat,cv)
C- Makes smoothed ves from smoothed density and qmom, incl nuc. charge
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxl rg lfoca rfoca qc z ctail etail stc orhoc lmxb p
Ci                pz rmt
Ci     Stored:    *
Ci     Passed to: corprm
Ci   s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: nabc vol
Ci     Stored:    *
Ci     Passed to: *
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   ng    :number of G-vectors
Ci   g2    :square of G-vectors
Ci   yl    :spherical harmonics
Ci   cs    :vector of cosines for the ng vectors
Ci   sn    :vector of sines for the ng vectors
Ci   iv
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Co Outputs
Co   cv    :local gaussian density added to cv
Co         :estatatic potential make from density
Cr Remarks
Cr   Local charge consists of a sum of gaussians that compensate for
Cr   the difference in multipole moments of true and smooth local charge
Cr   and a contribution from the smooth core charge.
Cr     g(qmpol) + g(qcore-z) + h(ncore)
Cr
Cr   Adapted from vesgcm to make strictly FT ves(nloc)
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ng,iv(ng,3)
      double precision qmom(*),g2(ng),yl(ng,*),cs(ng),sn(ng),qlat(3,3)
      double complex cv(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer ig,ib,ilm,is,iv0,l,lmxl,m,nbas,nlm,nlmx,nglob,n1,n2,n3,
     .  ngabc(3),lfoc
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (nlmx=81)
      double precision tau(3),df(0:20),pi,y0,vol,rg,qcorg,qcorh,qsc,
     .  cofg,cofh,ceh,rfoc,z,q0(3),gam,gamf,cfoc,cvol,aa
      double complex cof(nlmx),cfac,phase
C      parameter (k0=3)
C      double precision gv(ng,3),v(3)
C      double complex gkl(0:k0,nlmx)
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf4')
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      nbas = nglob('nbas')
      ngabc = s_lat%nabc
      vol = s_lat%vol

C --- FT of gaussian density, all sites, for list of G vectors ---
      iv0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        lmxl = s_spec(is)%lmxl
        rg = s_spec(is)%rg
        if (lmxl == -1) cycle

        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
        nlm = (lmxl+1)**2
        if (nlm > nlmx) call rxi('pvdf4: increase nlmx to',nlm)
        ilm = 0
        cfac = dcmplx(0d0,1d0)
        do  l = 0, lmxl
          cfac = cfac*dcmplx(0d0,-1d0)
          do  m = -l, l
            ilm = ilm+1
            cof(ilm) = cfac*qmom(ilm+iv0)*4*pi/df(2*l+1)
          enddo
        enddo
        cof(1) = cof(1) + 4*pi*y0*(qcorg-z)

        gam = 0.25d0*rg*rg
        gamf = 0.25d0*rfoc*rfoc
        cfoc = -4d0*pi*y0*cofh/vol
        cvol = 1d0/vol
        do  ig = 1, ng
          phase = dcmplx(cs(ig),sn(ig))
          aa = dexp(-gam*g2(ig))*cvol
          do  ilm = 1, nlm
            cv(ig) = cv(ig)+aa*yl(ig,ilm)*cof(ilm)*phase
          enddo
C     ... Add foca hankel part
          aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
          cv(ig) = cv(ig) + cfoc*aa*phase
        enddo

        iv0 = iv0+nlm
      enddo

C --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
      cv(1) = (0d0,0d0)
      do  ig = 2, ng
        cv(ig) = (8*pi)*cv(ig)/g2(ig)
      enddo

      call tcx('pvdf4')
      end
