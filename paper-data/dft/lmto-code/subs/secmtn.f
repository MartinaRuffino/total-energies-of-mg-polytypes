C#define BLAS3
      subroutine secmtn(s_ctrl,s_spec,s_lat,s_ham,s_pot,s_str,
     .  nl,nsp,nbas,ips,indxsh,qss,eula,neul,ikp,nkp,
     .  ldim,lidim,lihdim,qp,ppn,sop,isp,nevmx,efmax,nev,z,eb)
C- Hamiltonian and Overlap, NMTO
C ----------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: lncol
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:lgen3 lqp
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxb hcr
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: plat alat nkd nkq awald vol avw
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:pos dlv qlv cg jcg indxcg
Cio    Passed to: sstrxq
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: nmto kmto ldham
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:iprmb
Cio    Passed to: sstrxq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: vmtz
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read: npr
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:nkaps iax kaps s sdot nitab
Cio    Passed to: *
Ci Inputs:
Ci   nl,nsp,nbas,pp
Ci   ldim : dimension of l - wave block of hamiltonian matrix
Ci   lidim: ldim+dimension of i - wave block of hamiltonian matrix
Co   pph  : vector of pot. par's. in alpha rep'n (suham.f)
Ci   eula : Euler angles for spin rotations
Ci   nevmx: max. no. evecs to generate.
Ci          -1 suppresses generation of z
Ci          -2 Do not diagonalize, but return overlap in z,
Ci             allocate oc for hamiltonian and place there
Ci   z    :used as a work array, whether or not evecs generated
Co Outputs:
Ci   ccd:  diagonal matrices for 1- 2- & 3-centre CCOR integrals
Co   eigenvalues and eigenvectors are returned in eb, z
Co   nev:  number of evecs generated.
Cr Remarks
Cr   Downfolding automatically turns on the combined correction.
Cr   bittst(lham,8) can be used to transform
Cr   structure constants to an arbitrary representation.
Cr   Hybridisation is turned off when bittst(lham,16) is set (see remhyb.f).
Cr   Dimensions of pph,eb,z are doubled when spins are coupled.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   08 Feb 03 altered dmensions of sop
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nbas,neul,isp,ldim,lidim,lihdim,nevmx,nev,
     .  ips(nbas),ikp,nkp,indxsh(1)
      integer nppn
      parameter (nppn=12)
      double precision eb(ldim*2),ppn(nppn,lihdim,*),qp(3),qss(4),
     .  z(ldim,ldim*2),eula(nbas,neul,3),efmax,sop(0:nl-1,nsp,nsp,9,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
C ... Dynamically allocated arrays
      integer, allocatable :: lmx(:)
      real(8), allocatable :: tral(:),trad(:),hcr(:)
      complex(8), allocatable :: hk(:),ok(:),wk(:),sll(:),sdll(:),shl(:)
C ... Local parameters
      logical lewald,lnc,lx
      integer hdim,i,ii,ipr,j,l2,ldimx,lgunit,linv,lncol,nspec,
     .  nmto,n0,nsite,lbloch,nl2,ndoti
      parameter (n0=10)
      double precision vmtz,avw,dglob,kmto(n0),plat(3,3),xx
      character*80 outs
      procedure(integer) :: nglob

      lncol = s_ctrl%lncol
      nmto = s_ham%nmto
      kmto(1:6) = s_ham%kmto

      call getpr(ipr)
C     lewald = lgors('ctrl lgen3,2',sctrl)
      lewald = IAND(s_ctrl%lgen3,2) /= 0
      nspec = nglob('nspec')
      avw = dglob('avw',0d0,0)
      vmtz = s_pot%vmtz
      lnc = lncol /= 0
      nl2 = nl**2
      hdim = lihdim - lidim
C     idim = lidim - ldim
      ldimx = ldim
C     ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      l2 = ldimx**2
      lbloch = 10
      plat = s_lat%plat
C     ckbas = cksumf(s_lat%pos,3*nbas)

C     Space for nmto strux
      ii = lidim**2*nmto
      allocate(sll(ii),sdll(ii))
      if (hdim /= 0) then
        allocate(shl(hdim*lidim))
      else
        allocate(shl(1))
      endif

C --- Screened strux by Ewald summation ---
      if (lewald) then

C   ... Make the tral matrix on the fly
        ii = 4*nl**2*nbas*nmto
        allocate(tral(ii)); call dpzero(tral,ii)
        allocate(trad(ii)); call dpzero(trad,ii)
        allocate(lmx(nspec),hcr(nl*nspec))
        call spec2class(s_spec,nspec,0,'lmxb',1,lmx,xx)
        call spec2class(s_spec,nspec,0,'hcr',nl,xx,hcr)
        call dscal(nmto,avw**2,kmto,1)
        call dscal(nl*nspec,1/avw,hcr,1)
        call pshpr(1)
        if (ikp == 1 .and. isp == 1 .and. ipr >= 50) call setpr(50)
        call mktra2(1,1,nbas,ips,nl,lmx,avw,4,kmto,nmto,hcr,1,
     .    tral,trad,xx,xx,xx,xx)
        call poppr
        call dscal(nmto,1/avw**2,kmto,1)
        call dscal(nl*nspec,avw,hcr,1)
        deallocate(lmx,hcr)

        call sstrxq(1,s_lat,s_ham,qp,nmto,kmto,tral,trad,
     .    sll,sdll,shl,xx)

        deallocate(tral,trad)

C --- Screened strux by Bloch sum of R.S. strux ---
      else
        nsite = s_str%npr(nbas+1)
C       Energy-interpolate and Bloch sum strux
        ndoti = min(s_str%nkaps,2)
        call dscal(nmto,avw**2,kmto,1)
        call blochi(lbloch,qp,nl,plat,indxsh,1,nsite,s_str%iax,
     .    s_str%nkaps,ndoti,s_str%kaps,nmto,kmto,s_str%s,s_str%sdot,nl2,
     .    s_str%nitab,0,lidim,0,lidim,lidim,lidim,0,sll)
        call blochi(10000+lbloch,qp,nl,plat,indxsh,1,nsite,s_str%iax,
     .    s_str%nkaps,ndoti,s_str%kaps,nmto,kmto,s_str%s,s_str%sdot,nl2,
     .    s_str%nitab,0,lidim,0,lidim,lidim,lidim,0,sdll)
        call dscal(nmto,1/avw**2,kmto,1)
      endif

C --- Kink matrix from structure constants ---
      call dscal(2*lidim**2*nmto,avw**2,sdll,1)
      call sstr2k(ldim,lidim,lihdim,nmto,isp,ppn,sll,sdll)

C --- NMTO Hamiltonian and overlap ---
      if (allocated(hk)) deallocate(hk)
      allocate(hk(l2)); call dpzero(hk,2*l2)
      if (allocated(ok)) deallocate(ok)
      allocate(ok(l2)); call dpzero(ok,2*l2)
      allocate(wk(l2))
      call nmham(1,lidim,ldim,nmto,kmto,vmtz,sll,sdll,wk,hk,ok)
      deallocate(wk)

C --- Diagonalize ---
C#ifdef BLAS3
      lx = .true.
C#elseC
C      lx = .false.
C#endif
C     Diagonalize by inverse iteration, or not
      linv = 0
C     if (nevmx > 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
      if (nevmx > 0 .and. IAND(s_ctrl%lqp,2) /= 0) linv = 1
      allocate(wk(ldimx*11))
      call zhev(ldimx,hk,ok,.true.,lx,nevmx,efmax,nev,wk,linv,0,eb,z)
      deallocate(wk)

C --- Printout ---
      if (ipr >= 30) then
        j = min(9,ldimx)
        if (ipr >= 35) j = ldimx
C#ifdefC LINUX_PGI
C        do  18  ii = 1, 1
C#else
        do  18  ii = 1, 2
C#endif
        call awrit3(' SECMTN:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
        write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
   18   continue
        if (ipr >= 36 .and. nev > 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,lgunit(1),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif

      if (ipr >= 110) then
        outs = 'evec'
        call yprm(outs,2,z,ldimx*nev,ldimx,ldimx,nev)
        call zprm('eigenvectors',2,z,ldimx,ldimx,ldimx)
        call yprm('eval',1,eb,ldimx*1,ldimx,nev,1)
        call query('V<110 to skip matrix printing',-1,0)
      endif

      deallocate(sll,sdll,shl)
      end

      subroutine sstr2k(ldim,lidim,lihdim,nmto,isp,ppn,sll,sdll)
C- Overwrite screened structure constants with Kink matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   nmto :number of types of one l-quantum number in the basis
Ci   isp   :current spin channel (1 or 2)
Ci   ppn  :nmto potential parameters, in downfolding order
Cio Inputs/Outputs
Cio  sll   :On input, screened structure constants
Cio        :On output, K = aD - S
Cr Remarks
Cr  For downfolding, the KKR-equations have the simple form:
Cr
Cr    K^a *  u = 0   with  K^a = S^a + aD
Cr
Cr and where K^a is the kink matrix in the a-representation.
Cr Divide orbitals into lower and intermediate sets:
Cr Then u_i can be expressed by u_l:
Cr
Cr    u_i = - (K^a_ii)^-1 * K^a_il * u_l
Cr
Cr  and the equation for u_l is:
Cr
Cr    (K^a_ll - K^a_li * (K^a_ii)^-1 * K^a_il) u_l \equiv K^b_ll u_l = 0
Cr
Cr   with beta = inverse potential function P^0.
Cr
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldim,lidim,lihdim,nppn,nmto,isp
      parameter (nppn=12)
      double precision ppn(nppn,lihdim,nmto,isp)
      double complex sll(lidim,lidim,nmto),sdll(lidim,lidim,nmto)
C ... Local parameters
      integer inc,ik,i,j,idim,ierr,iprint
      double precision ad,add,phia,phida
      double complex alpha,beta
      character*40 outs
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: zwk(:)

C --- Downfolded kink matrix and its energy derivative ---
      idim = lidim-ldim
      if (idim /= 0) then

      do  ik = 1, nmto

        alpha = 1
        beta = 0

C   ... Scale the SSW ii block to match the back extrapolated phi(a)
        do  i = ldim+1, lidim
          do  j = ldim+1, lidim
            sdll(i,j,ik) = sll(i,j,ik) * ppn(6,j,ik,isp) +
     .                    sdll(i,j,ik) * ppn(5,j,ik,isp)
            sll(i,j,ik) =  sll(i,j,ik) * ppn(5,j,ik,isp)
          enddo
        enddo

C   ... Add diagonal D to make iwaves K, using SSW rescaled to phi0(a)
        do  i = ldim+1, lidim
             aD = ppn(3,i,ik,isp)
            aDd = ppn(4,i,ik,isp)
           phia = ppn(5,i,ik,isp)
          phida = ppn(6,i,ik,isp)

          sll(i,i,ik)  =  sll(i,i,ik) + aD*phia
          sdll(i,i,ik) = sdll(i,i,ik) + aDd*phia + aD*phida
        enddo

C   ... Kii -> Kii^-1; Overwrite Sil with Kii^-1 Sil = Kii^-1 Sli+
        allocate(wk(idim*max(idim+1,3)))
        i = ldim+1
        call zqinv('g',sll(i,i,ik),lidim,0,idim,wk,idim,ierr)
        if (ierr /= 0) call rx('sstr2k:  Kii is singular')
        deallocate(wk)
        call zgemm('N','C',idim,ldim,idim,alpha,sll(i,i,ik),lidim,
     .    sll(1,i,ik),lidim,beta,sll(i,1,ik),lidim)

C   ... Store in zwk : Sdil - Kdii Kii^-1 Sil
        allocate(zwk(idim*ldim))
        call zgemm('N','N',idim,ldim,idim,alpha,sdll(i,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,zwk,idim)
        call dscal(idim*ldim*2,-1d0,zwk,1)
        call zmsadd(idim,ldim,lidim,idim,0,0,0,0,alpha,alpha,
     .    sdll(i,1,ik),zwk)

C   ... Scale the SSW li block to match the back extrapolated phi(a)
C       Also scale by -1 to prepare for the final multiplication steps
        do  i = ldim+1, lidim
          phia  = ppn(5,i,ik,isp)
          phida = ppn(6,i,ik,isp)
          do  j = 1, ldim
            sdll(j,i,ik) =-sll(j,i,ik) * phida - sdll(j,i,ik) * phia
            sll(j,i,ik) = -sll(j,i,ik) * phia
          enddo
        enddo

C   ... Overwrite Sdil with : Kii^-1 (Sdil - Kdii Kii^-1 Sil)
        i = ldim+1; beta = 0
        call zgemm('N','N',idim,ldim,idim,alpha,sll(i,i,ik),lidim,
     .    zwk,idim,beta,sdll(i,1,ik),lidim)
        deallocate(zwk)

C   ... Kll(b) = Kll(a) - Kli * (Kii^-1 * Kil)
        i = ldim+1; beta = 1
C       call zprm('-Kli',2,sll(1,i,ik),lidim,ldim,idim)
C       call zprm('Kii^-1 Kil',2,sll(i,1,ik),lidim,idim,ldim)
        call zgemm('N','N',ldim,ldim,idim,alpha,sll(1,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,sll(1,1,ik),lidim)

C       call zprm('K^a - Kli Kii^-1 Kil',2,sll(1,1,ik),lidim,ldim,ldim)

C   ... Kdll(b) = Kdll(a) - Kdli * (Kii^-1 Kil)dot - Kdli * (Kii^-1 Kil)
        call zgemm('N','N',ldim,ldim,idim,alpha,sll(1,i,ik),lidim,
     .    sdll(i,1,ik),lidim,beta,sdll(1,1,ik),lidim)
        call zgemm('N','N',ldim,ldim,idim,alpha,sdll(1,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,sdll(1,1,ik),lidim)

C       call zprm('Kd^b',2,sdll(1,1,ik),lidim,ldim,ldim)

      enddo
      endif

C --- Kink matrix of ll block ---
      inc = 2*(1+lidim)
      do  ik = 1, nmto
        call daxpy(ldim,1d0,ppn(3,1,ik,isp),nppn,sll(1,1,ik),inc)
        call daxpy(ldim,1d0,ppn(4,1,ik,isp),nppn,sdll(1,1,ik),inc)

        if (iprint() > 110) then
          call awrit2('%xK^b%?#n#-dot##, ik=%i',outs,len(outs),0,0,ik)
          call zprm(outs,2,sll(1,1,ik),lidim,ldim,ldim)
          call awrit2('%xK^b%?#n#-dot##, ik=%i',outs,len(outs),0,1,ik)
          call zprm(outs,2,sdll(1,1,ik),lidim,ldim,ldim)
        endif

      enddo

      end
