      subroutine ovlocr(nbas,s_site,s_spec,s_lat,s_ham,nxi0,nxi,exi,hfc,
     .  rsmfa,s_rhofa,sqloc,slmom)
C- Makes the site densities for overlapped free atoms.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl kmxv rsmv a nr rmt lfoca rfoca qc z ctail etail
Ci                 stc lmxb p pz rg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat awald tol vol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  hxpbl ghibl hklbl gklbl
Cio  s_rhofa : free-atom density, by species
Ci Inputs
Ci   nbas  :size of basis
Ci   nxi   :number of Hankels
Ci   nxi0  :leading dimension of hfc
Ci   exi   :smoothed Hankel energies; see Remarks
Ci   hfc   :coefficients to smoothed Hankels
Ci   rsmfa :Hankel smoothing radius
Co Outputs
Co   sqloc :sum of local charges (integral over rho1-rho2)
Co   slmom :sum of local magnetic moments
Co         :Collinear case:
Co         :slmom(0) = sum of local moments
Co         :Noncollinear case:
Co         :slmom(1:3) is a vector
Co         :slmom(0) = |slmom(1:3)|
Cl Local variables
Cl   lrots :Density is noncollinear for this site
Cl   ns4 = 1 if local density is not spin polarized
Cl         2 if local spin density is collinear along z
Cl         4 if local spin density is to be rotated;
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 May 07 parallelized (MPI)
Cu   01 Jul 05 Zero-radius empty spheres skip as having no local part
Cu   14 Jun 00 spin polarized
Cu   24 Apr 00 Adapted from nfp ovlocr.f
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nxi(*),nxi0
      double precision rsmfa(*),exi(nxi0,*),
     .  hfc(nxi0,2,1),sqloc,slmom(0:3)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(dp_wk_vec) :: s_rhofa(nbas)
C ... Local parameters
      logical lrots
      integer nrmx,nlmx,ktop0,stdo
      parameter(nrmx=5001,nlmx=64,ktop0=25)
      integer i,ib,intopt,ipr,iprint,is,jb,je,js,kmxv,lfoca,lgunit,lmxl,
     .  nglob,nlmh,nlml,nr,nsp,nspc,ns4
      double precision a,amom(0:3),ceh,cofg,cofh,eh,qcorg,qcorh,qcsm,
     .  qloc,qsc,rfoca,rmt,rsmh,rsmv,z,eula(3)
      double precision rofi(nrmx),rwgt(nrmx),p1(3),
     .  p2(3),b0(ktop0+1),q(3)
      double complex b((ktop0+1),nlmx) !,u(2,2)
      real(8),target :: acof((ktop0+1)*nlmx*4)
      real(8),pointer :: acofl(:)


      integer, dimension(:),allocatable :: ibproc
      integer ierr,pid
      integer procid, master, numprocs
      logical mlog,cmdopt
      character strn*120
C     double precision dval
C ... Heap
      data q /0d0,0d0,0d0/

      call tcn('ovlocr')
      ipr  = iprint()
      stdo = lgunit(1)
      nsp  = nglob('nsp')
      nspc = nglob('nspc')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"

      intopt = 10*nglob('lrquad')
      call mpi_comm_rank( mpi_comm_world, procid, ierr)
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

C --- Loop over sites where density is made ---
      sqloc = 0
      slmom = 0

      if (ns4 == 4) then
        allocate(acofl((ktop0+1)*nlmx*4))
      else
        acofl => acof
      endif
      allocate (ibproc(0:numprocs))
      if (numprocs > 1) then
      call info0(31,1,0,' ovlocr: make sphere densities from overlapping FA densities (parallel)')

      call pshpr(ipr-10)
      call dstrbp(nbas,numprocs,1,ibproc(0))
      call poppr
      ipr = 0
      call pshpr(10)
      else
      ibproc(0:1) = [1,nbas+1]
      endif
      do  ib = ibproc(procid), ibproc(procid+1)-1
        is = s_site(ib)%spec
        lrots = .false.
        p1 = s_site(ib)%pos
        lmxl = s_spec(is)%lmxl
        if (numprocs == 1 .and. lmxl >= 0) then
C        if (ipr >= 30) write (stdo,1) ib
C    1     format(/' Free atom and overlapped crystal site charges site',i5,':'/
C     .        '   ib    true(FA)    smooth(FA)  true(OV)',
C     .   '    smooth(OV)    local')
        call info2(30,1,0,' Free atom and overlapped crystal site charges  spec %i  site %i:%N'//
     .        '   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local',is,ib)
        endif

        kmxv = s_spec(is)%kmxv
        rsmv = s_spec(is)%rsmv
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoca,rfoca,z)
        qcsm = qcorg+qcorh
        nlml = (lmxl+1)**2
C       if (z == 0 .and. rmt == 0) goto 10
        if (lmxl == -1) cycle
        if (nr > nrmx) call rxi('ovlocr: increase nrmx, need',nr)
        if (kmxv > ktop0) call rxi('ovlocr: increase ktop0 need',kmxv)
        if (nlml > nlmx) call rxi('ovlocr: increase nlmx need',nlml)

        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)

C   ... Loop over other sites, add up tail expansion
        call dpzero(acof,(kmxv+1)*nlml*ns4)
        do  jb = 1, nbas
          if (ns4 == 4) then
            call dpzero(acofl,(kmxv+1)*nlml*ns4)
          endif
          js = s_site(jb)%spec
          p2 = s_site(jb)%pos
          do  je = 1, nxi(js)
            rsmh = rsmfa(js)
            eh   = exi(je,js)
            nlmh = 1
            call hxpbl(p2,p1,q,rsmh,rsmv,eh,kmxv,nlmh,nlml,kmxv,
     .        nlml,s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy,s_lat,b)
            call dpzero(b0,(kmxv+1)*nlmh)
            if (ib == jb) call hxpos(rsmh,rsmv,eh,kmxv,nlmh,kmxv,b0)
            do  i = 1, nsp
              call p1ovlc(kmxv,nlml,hfc(je,i,js),b,b0,
     .          acofl(1+(1+kmxv)*nlml*(i-1):))
            enddo
          enddo

c     ... rotate to quantization axis specified by euler angles
          if (ns4 == 4 .and. s_ham%neula /= 1) then
            call rx('ovlocr: eula should be l-independent')
          endif
c         if (ipr >= 40)
c     .      write(stdo,700) 'euler:',(s_ham%eula(ib,je),je=1,3)
c  700     format(2x,a7,6f11.3)
C         lrotj = ns4 == 4    ! analog of lrots for this jb
c           if (lrotj) then  ! skip this shortcut since maybe m<0 -> m>0
c             lrotj = dlength(3,s_ham%eula(jb,1),nbas) > 0
c           endif
          if (ns4 == 4) then
            i = (kmxv+1)*nlml
            call rotspv(644,s_ham%eula(jb,1:3),i,i,i,
     .        acofl,acofl,acofl,acofl)
C            call rotspu(0,1,1,1,1,s_ham%eula(jb,1),1,u)
C            call rotspv(44,u,i,i,i,acofl,acofl,acofl,acofl)
            call daxpy(i*4,1d0,acofl,1,acof,1)
            lrots = .true. ! spin density at ib is noncollinear
          endif
        enddo

        if (lrots) eula(1:3) = s_ham%eula(ib,1:3)
        if (s_ham%lrsa >= 100 .and. nsp == 2) then
          call splrho(10*(1+4),nsp,nr,1,z,z,s_site(ib)%rhoc)
          s_site(ib)%rhoc(:,2) = 0d0
          call splrho(10*(1+4)+1,nsp,nr,1,z,z,s_site(ib)%rhoc)
        endif
        call p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,lrots,eula,
     .    acof,rofi,rwgt,nxi0,nxi(is),exi(1,is),hfc(1,1,is),rsmfa(is),
     .    s_rhofa(is)%p,s_site(ib)%rhoc,lfoca,qcsm,qloc,amom,
     .    s_site(ib)%rho1,s_site(ib)%rho2)

        sqloc = sqloc + qloc
        slmom = slmom + amom
      enddo

C ... Combine sphere densities from separate threads
      if (numprocs > 1) then
        call poppr

        do  pid = 0, numprocs-1
          do  ib = ibproc(pid), ibproc(pid+1)-1
            is = s_site(ib)%spec
            nr = s_spec(is)%nr
            lmxl = s_spec(is)%lmxl
            if (lmxl<0) cycle
            nlml = (lmxl+1)**2
            call mpibc3(s_site(ib)%rho1,nr*nlml*nsp,4,pid,0,'ovlocr','rhoat(1)')
            call mpibc3(s_site(ib)%rho2,nr*nlml*nsp,4,pid,0,'ovlocr','rhoat(1)')
            call mpibc3(s_site(ib)%rhoc,nr*nsp,4,pid,0,'ovlocr','rhoat(1)')
C           call MPI_BCAST(s_site(ib)%rho1,nr*nlml*nsp,mpi_real8,pid,MPI_COMM_WORLD,ierr)
C           call MPI_BCAST(s_site(ib)%rho2,nr*nlml*nsp,mpi_real8,pid,MPI_COMM_WORLD,ierr)
C           call MPI_BCAST(s_site(ib)%rhoc,nr*nsp,mpi_real8,pid,MPI_COMM_WORLD,ierr)
          enddo
        enddo
C   ... Combine sum-of-sphere charges from separate threads
        call mpibc2(sqloc,1,4,3,mlog,'ovlocr','sqloc')
        call mpibc2(slmom,4,4,3,mlog,'ovlocr','slmom')

C       Debugging printout
C       print 345, procid, (dval(s_site(ib)%rho1,445), ib=1,20)
C       print 345, procid, (dval(s_site(ib)%rho2,445), ib=1,20)
C       print 345, procid, (dval(s_site(ib)%rhoc,445), ib=1,20)
C  345  format(i4,20f10.5)
C       print *, procid, iprint()

        deallocate(ibproc, stat=ierr)
C     call rx0('done')
      endif

      call tcx('ovlocr')
      end

      subroutine p1ovlc(kmxv,nlml,hfc,b,b0,a)
C- Adds contribution to P_kl expansion of density from one basis function
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlml  :density expanded to nlml
Ci   kmxv  :k-cutoff for P_kl expansion
Ci   hfc   :coefficient to basis function
Ci   b     :P_kl expansion of density from one basis function
Ci   b0    :P_kl expansion of on-site density
Co Outputs
Co   a     :cumulative P_kl expansion of density for this site
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlml,kmxv
      double precision a(0:kmxv,nlml),b0(0:kmxv,*),hfc
      double complex b(0:kmxv,nlml)
C ... Local parameters
      integer k,ilm

      do  k = 0, kmxv
        do  ilm = 1, nlml
          a(k,ilm) = a(k,ilm) + hfc*dble(b(k,ilm))
        enddo
        a(k,1) = a(k,1) - hfc*b0(k,1)
      enddo
      end

      subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,lrots,eula,acof,rofi,
     .  rwgt,nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,
     .  amom,rho1,rho2)
C- Assemble local density from P_kl expansion for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ib    :site for which to assemble local density
Ci   nsp   :number of spin channels
Ci   rsmv  :smoothing radius for P_kl expansion
Ci   kmxv  :k-cutoff for P_kl expansion
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for P_kl expansion
Ci   lrots :coefficients acof are noncollinear; see Remarks
Ci   eula  :Euler angles for spin rotation of head density (lrots=T)
Ci   acof  :coefficients to P_kl expansion
Ci   rofi  :radial mesh points for tabulation on a radial mesh
Ci   rwgt  :radial mesh weights for integration on a radial mesh
Ci   nxi0  :leading dimension of hfc
Ci   nxi   :number of smoothed Hankel energies in head expansion
Ci   exi   :smoothed Hankel energies in head expansion
Ci   hfc   :coefficients to Hankel energies in head expansion
Ci   rsmfa :Hankel smoothing radius in head expansion
Ci   rhofa :head free-atom density
Ci   rhoc  :core density --- used to integrate core charge
Ci   lfoca :switch specifying treatment of core density.
Ci          0 => val,slo = 0 at sphere boundary
Ci          1 => core tails included explicitly with valence
Ci          2 => tails included perturbatively
Ci   qcsm  :smoothed core density
Co Outputs
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Co   qloc  :sphere charge
Co   amom  :sphere magnetic moment (x,y,z components)
Cl Local variables
Cl   rhohd :on-site smoothed head density
Cr Remarks
Cr   In the noncollinear case acof has 4 spin components
Cr   stored in compressed format (see mode=4 in rotspv)
Cr   A noncollinear density will be stored in rho1,rho2,rhoc
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lrots
      integer nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca
C     Noncollinear case: rho1,rho2,rhoc have 4 components
      double precision qcsm,qloc,rsmv,rsmfa,amom(0:3),eula(3),
     .  rho1(nr,nlml,nsp*2),rofi(nr),rwgt(nr),exi(1),hfc(nxi0,nsp),
     .  rho2(nr,nlml,nsp*2),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp*2)
      real(8),target :: rhofa(nr,nsp)
C ... Local parameters
      integer kmx,lmx,stdo
      parameter (kmx=20, lmx=6)
      integer i,ie,ilm,ipr,iprint,k,l,lgunit,ll,lmax,lmxl,isp,nspc,isw
      double precision asm,gam,pi,qall,qexa,qin,qlc,qnum,qout,qsmo,qut,
     .  r,rl,rmt,srfpi,sum,sumfa,sumhd,sumsm,sumtr,y0,dlength,
     .  amomfa,                 ! magnetic moment of free atom
     .  amomsm,                 ! magnetic moment of smoothed free atom
     .  amomt(0:3),             ! true moment of overlapped atom
     .  amoms(0:3),             ! smoothed moment of overlapped atom
     .  pkl(0:kmx,0:lmx),xi(0:10),x0(0:2),ddot,
     .  rhohd(nr,4)   ! 2 if collinear magnetic; 4 if noncollinear
      double complex u(2,2)
      real(8),pointer :: rhofl(:,:)

      if (lrots .and. nsp /= 2) call rx('p2ovlc: nonsensical input')

      ipr   = iprint()
      stdo  = lgunit(1)
      pi    = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
      y0    = 1d0/srfpi
      lmxl  = ll(nlml)
      nspc  = 1 + isw(lrots)
      if (nspc == 2) then
        allocate(rhofl(nr,4)); call dpzero(rhofl,nr*4)
        call dcopy(2*nr,rhofa,1,rhofl,1)
      else
        rhofl => rhofa
      endif
      amom = 0

      if (lmxl > lmx) call rxi('ovlocr: increase lmx, need',lmxl)

C     do  ilm = 1, nlml
C       do  k = 0, kmxv
C         if (dabs(acof(k,ilm,1)) > 1d-6)
C    .      write(stdo,780) ilm,k,acof(k,ilm,1),acof(k,ilm,nsp)
C 780     format('ilm,k',2i5,2f14.8)
C       enddo
C     enddo

C --- Assemble smooth on-site head density in rhohd ---
      qnum = 0d0
      qexa = 0d0
      qsmo = 0d0
      qut = 0d0
      call dpzero(rhohd, nr*nsp*nspc)
      asm = 1d0/rsmfa
      lmax = 0
      do  ie = 1, nxi
        sum = 0d0
        do  i = 1, nr
          r = rofi(i)
          call hansmr(r,exi(ie),asm,xi,lmax)
          sum = sum + srfpi*rwgt(i)*xi(0)*r*r
          do  isp = 1, nsp
            rhohd(i,isp) = rhohd(i,isp) + srfpi*hfc(ie,isp)*xi(0)*r*r
          enddo
        enddo
        gam = 0.25d0*rsmfa**2
        qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
        rmt = rofi(nr)
        call hansmr(rmt,0d0,1/rsmfa,x0,1)
        call hansmr(rmt,exi(ie),1/rsmfa,xi,1)
        qout = srfpi/exi(ie)*(-dexp(rsmfa**2/4*exi(ie))
     .     - rmt**3*(xi(1)-dexp(rsmfa**2/4*exi(ie))*x0(1)))
        qin = qall-qout
        do  isp = 1, nsp
          qnum = qnum + hfc(ie,isp)*sum
          qexa = qexa + hfc(ie,isp)*qin
          qsmo = qsmo + hfc(ie,isp)*qall
          qut  = qut  + hfc(ie,isp)*qout
        enddo
      enddo

C     Magnetic moments before spin rotations
      amomfa = ddot(nr,rwgt,1,rhofa,1)-ddot(nr,rwgt,1,rhofa(1,nsp),1)
      amomsm = ddot(nr,rwgt,1,rhohd,1)-ddot(nr,rwgt,1,rhohd(1,nsp),1)

C     Rotate rhofa and rhohd
      if (lrots) then
        call rotspu(0,1,1,1,1,eula,1,u)
        call rotspv(44,u,nr,nr,nr,rhofl,rhofl,rhofl,rhofl)
        call rotspv(44,u,nr,nr,nr,rhohd,rhohd,rhohd,rhohd)
      endif

c|      write(stdo,917) qnum,qexa,qsmo,qut
c|  917 format('summed smooth charge:  num',f14.8,'   exact',f14.8
c|     .   /' total smooth q',f14.8,'  outside',f14.8)

C --- Assemble overlapped tail density in rho2 ---
      if (kmxv > kmx) call rx('ovlocr: increase kmx')
      call dpzero(rho2,  nr*nlml*nsp*nspc)
      do  i = 1, nr
        r = rofi(i)
        call radpkl(r,rsmv,kmxv,lmxl,kmx,pkl)
        do  isp = 1, nsp*nspc
        do  ilm = 1, nlml
          l = ll(ilm)
          rl = 0.d0
          if ( r > 0.d0 ) rl = r**l
          do  k = 0, kmxv
            rho2(i,ilm,isp) = rho2(i,ilm,isp) +
     .                        acof(k,ilm,isp)*pkl(k,l)*r*r*rl
          enddo
        enddo
        enddo
      enddo
C     Copy the tails to rho1
      call dpcopy(rho2,rho1,1,nr*nlml*nsp*nspc,1d0)

C ... Make the true density in rho1, smooth density in rho2
      do  isp = 1, nsp*nspc
        do  i = 1, nr
          rho1(i,1,isp) = rho1(i,1,isp) + y0*rhofl(i,isp)
          rho2(i,1,isp) = rho2(i,1,isp) + y0*rhohd(i,isp)
        enddo
      enddo

C ... Do some integrals
      sumfa = 0d0
      sumsm = 0d0
      sumhd = 0d0
      sumtr = 0d0
      qlc = 0d0
      do  isp = 1, nsp
        do  i = 1, nr
        sumfa = sumfa + rwgt(i)*rhofa(i,isp)
        sumhd = sumhd + rwgt(i)*rhohd(i,isp)
        sumtr = sumtr + rwgt(i)*rho1(i,1,isp)
        sumsm = sumsm + rwgt(i)*rho2(i,1,isp)
        qlc = qlc + rwgt(i)*rhoc(i,isp)
        enddo
      enddo
      sumsm = sumsm*srfpi
      sumtr = sumtr*srfpi
      qloc = sumtr-sumsm

C     Magnetization (-> z component if nspc=2)
      amomt(0) = srfpi*
     .  (ddot(nr,rwgt,1,rho1,1)-ddot(nr,rwgt,1,rho1(1,1,nsp),1))
      amoms(0) = srfpi*
     .  (ddot(nr,rwgt,1,rho2,1)-ddot(nr,rwgt,1,rho2(1,1,nsp),1))
      if (nspc == 2) then
        amomt(1) = 2*srfpi*ddot(nr,rwgt,1,rho1(1,1,3),1)
        amomt(2) =-2*srfpi*ddot(nr,rwgt,1,rho1(1,1,4),1)
        amoms(1) = 2*srfpi*ddot(nr,rwgt,1,rho2(1,1,3),1)
        amoms(2) =-2*srfpi*ddot(nr,rwgt,1,rho2(1,1,4),1)
        amomt(3) = amomt(0)
        amoms(3) = amoms(0)
        amom     = amomt - amoms
        amom(0) = dlength(3,amom(1),1)
        amomt(0) = dlength(3,amomt(1),1)
        amoms(0) = dlength(3,amoms(1),1)
      else
        amom(0) = amomt(0)-amoms(0)
      endif

      if (lfoca == 0) qloc = qloc + qlc - qcsm
      if (ipr >= 30) then
        write(stdo,810) ib,sumfa,sumhd,sumtr,sumsm,qloc
        if (nspc == 2) then
          write(stdo,812) '  Mx', amomt(1),amoms(1),
     .      'alpha =',eula(1)
          write(stdo,812) '  My', amomt(2),amoms(2),
     .      'beta  =',eula(2)
          write(stdo,812) '  Mz', amomt(3),amoms(3),
     .      'gamma =',eula(3)
        endif
        if (nsp == 2) then
          write(stdo,811) 'amom',
     .    amomfa,amomsm,amomt(0),amoms(0),amom(0)
        endif

      endif
  810 format(i5,6f12.6)
  811 format(1x,a,6f12.6)
  812 format(1x,a,24x,2f12.6,4x,a,f12.6)

      end
