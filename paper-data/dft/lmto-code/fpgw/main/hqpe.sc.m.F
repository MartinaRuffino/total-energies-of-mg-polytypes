      program hqpe_sc
C- make Sigma-Vxc
c------------------------------------------------------------------
c     calculates quasiparticle energies
c     E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
c     e(k,t) = LDA eigenvalue
c     Z      = [1 - dSEc(e(k,t))/dw]^(-1)
c     SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
c     SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
c     xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
c     SEx and xcLDA are in file SEX
c     SEc is in file SEC
c-----------------------------------------------------------------------
C sigm file is read, written, read, written, etc.
C Warning: this code needs cleaning up!
C File inputs:
C vxc.ext evec.ext GWinput BZDATA
C SEX2U SECU SEC2U SEXcoreU SEXcore2U SEXU XCU
C SEX2D SECD SEC2D SEXcoreD SEXcore2D SEXD XCD
C sigm from a prior iteration may be present --- mixed with current sigm
C ANFcond ldima
!     This is a developing version, not completely documented...
!     This uses amix in order to guess a better Sigma-Vxc from previous iterations.
!     * If ANFcond exists, it add dummy part (zero Sigma-Vxc ) for isp=2.
!     * iSigma_en==5 is for diagonal-only Sigma-Vxc in LDA basis set.
!     (Then you need evec0, which contains eigenvector of LDA).
!
Cl    ndham  :upper bound to the rank of the hamiltonian
Cl           :formerly ndimh and also nhq
Cl    ndsigm :upper bound to the number of states for which <i|sig|j> calculated (read for now from SEXU)
Cl           :formerly ntq
Cl    ndsigq :ndsigq(iq) = number of states for which <i|sig|j> calculated at iq
Cl    ndimh  :rank of the hamiltonian for this q
Cl           :formerly nz
Cl    ndimhq :ndimhq(iq) = hamiltonian rank for qp iq (read from vxc file)
Cl           :formerly nhqx
Cl    nmto   :rank of the LMTO part of hamiltonian (independent of q)
Cl    ndz    :rank of eigenvectors used to rotate hamiltonian to LMTO basis
Cl           :Takao's June 2014 code: this was always ndimh (formerly nz)
Cl    ndimsig:Dimensions sigm written to disk. nmto if mtosigmaonly()=T, ndham otherwise
Cl    nqsig  :number of qp at which sigm is calculated.
Cl           :Formerly nq
Cl    ... the following q vectors are confusing and should be cleaned up. For now:
Cl    qsig   :list of k-points for v_xc
Cl    qx     :list of k-points read from SEX file, and used for sigm file
Cl    qsigf  :list of k-points read from sigm file, and also used to write sigm file
Cu Updates
Cu   29 Apr  :Read mode=-1 => reads file vxc and dumps to file vxcsig in k-point order and format of file sigm
Cu           :Inputs required for this mode: sigm GWinput mpiqlst vxc.ext
      use keyvalue
      use m_read_bzdata, only: read_bzdata, nstar, nqibz2=>nqibz
      implicit none

      integer NULLI
      parameter (NULLI=-99999)

c     local data
      logical laf, lsigin
      integer ifianf,iopenxx,ifi,jin,nspin,nqsig,ndsigm,ndimsigin,
     .  nqx,lshft1,lshft2,lshft3,is,iSigma_en,izl,
     .  ndham,nsp,nspc,nqvxc,nspx,iq,nspinx,nsigmx,nqbz,nqibz,
     .  nqbzx,nqibzx,ip,i,isx,ix,icompvv2,itp,itpp,ikp,info,
     .  ildima,natom,ibas,ilmto,ibas1,levec,nstarsum,nevec,ndimx,neseavr
      integer,allocatable:: ndsigq(:)
      integer :: mode_sigfile = 0

      integer, allocatable :: iprocq(:) ! procid that generated data for given iq
      integer, allocatable :: ifevec__(:),ifvxc__(:),ifevecl__(:)  ! procid-dependent file pointers

      double precision deltaw,alat,ef,esmr,qxxx1,qxxx2,qxxx3,exxx,
     .  eshift0,eshift02,eshtlda,fac,eavr,eavr2,eseavr0,eseavr02,eee,beta
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      integer ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2),
     .  iftote(2),iftote2(2),ifsex2(2),ifsexcore2(2),ifsec2(2)
      integer(4),allocatable :: itxc(:),itc(:),itx(:)
      real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:),
     &     qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:),
     &     qx(:,:,:),eldax(:,:),rsec(:,:),csec(:,:) !sf..13May
     &     ,qsig(:,:,:)          !sf..3June
     &     ,qsigf(:,:,:)
      complex(8), allocatable :: sex2(:,:,:),sexcore2(:,:,:),sec2(:,:,:),se(:,:,:),
     .  work(:),evec_inv(:,:),evec_invt(:,:),ev_se_ev(:,:),v_xc(:,:,:,:)

      complex(8), target, allocatable :: evec(:,:,:,:),sigin(:,:,:,:)
      complex(8), pointer :: evecl(:,:,:,:),sigma_m(:,:,:,:)

C     complex(8), allocatable :: evec00(:,:,:,:) ,evec00inv(:,:)
C     complex(8), allocatable :: se_ev(:,:),sen(:,:),sen2(:,:),se_in(:,:),evec0(:,:)
      real(8) qx2(3),qtmp(3)
      integer, allocatable :: ipiv(:) !sf.end
      real(8) :: rydberg,hartree
      real(8) :: eseavrmean(2), eseavrmix(2), eseavr_file(2), threshold, eseadd
      integer(4) ::  n1,n2,n3, ifse_out
C     logical :: evec0ex=.false.

C     complex(8),allocatable:: sigma_m_out(:,:,:,:),pmat(:,:),pmatd(:,:)
      complex(8),allocatable:: sigmv(:,:,:)
      complex(8),allocatable:: zwk(:,:)

c----------------------------------------------
      integer(4),allocatable:: ibasf(:),ldima(:)
      real(8),allocatable :: pos(:,:)

      real(8):: anfvec(3)
      character(2):: soflag
      real(8) :: wex
      logical :: exonly,mtosigmaonly
      integer GWversion
      integer(4):: ret

      character(3):: iaaa
      real(8),allocatable:: eseavr(:,:,:),eseavr_in(:,:)
      real(8) :: xx(2), vxcscale

      integer,allocatable:: ndimhq(:,:),ndzq(:,:)
      integer:: ndimh, nmto,ndimsig,ndimsig2,ndz

      integer:: procid,nqcalc,nrank,ifvxc_,ifevec_,ifevecl_
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      integer, parameter :: master = 0

cccccccccccccccccccccccccccccccccccccccccccccccc
c     logical :: core3ptest=.false. ! This is for test for NiO_bigbas
c     ! Project out the 3p(LDA) part from Sigma-Vxc.
c     ! This requires evec0(LDA eigenfunctions)
cccccccccccccccccccccccccccccccccccccccccccccccc
c     if(core3ptest) print *,' core3ptest=',core3ptest

c --- Setup ---
      hartree=2d0*rydberg()
c     shift quasiparticle energies (eV)

      if (command_argument_count() == 0) then
        write (*,*) 'enter band index for Ef=0 [optional vxc_scale]'
        ext = ' '
        read (*,"(a100)",end=1100) ext ! ext used as char work array
 1100   continue
      else
        call get_command_argument(1, ext)
        if (command_argument_count() > 1) then
          call get_command_argument(2, ext(2+len_trim(ext):))
        endif
      endif

      jin = 0; vxcscale = 1
      call words(ext,ret)
      if (ret == 1) read(ext,*) jin
      if (ret >= 2) read(ext,*) jin, vxcscale

      call headver('hqpe_sc',jin)
      if (vxcscale /= 1) print "(' vxcscale =',f9.6)", vxcscale

      if (jin == -1) goto 301

      levec = IAND(mod(GWversion()/100,10),2)/2 ! levec=1 re-read evec from evecl for saving SE-vxc
      call getkeyvalue("GWinput","EXonly",wex,default=0d0,status=ret)
      if(wex==0d0) then
        exonly=.false.
      else
        exonly=.true.
        print *,' exonly=T wex=',wex
      endif

      inquire(file='ANFcond',EXIST=laf)
      if(laf) then
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*) soflag
        if(soflag=='SO') then
          print *,' ANFcond found. but SO case is not implimented. Not using ANF symm'
          laf=.false.
        else
          print *, 'Find ANFcond--- up only'
        endif
        close(ifianf)
      endif

      ifsex(1) = iopenxx ('SEXU')
      ifxc(1)  = iopenxx ('XCU')
      open(UNIT=3001, file='SEX2U',form='UNFORMATTED', status='OLD') !sf.beg
      ifsex2(1)=3001

      if(.not.exonly) then
        ifsec(1)     = iopenxx ('SECU')
        open(UNIT=3003,   file='SEC2U',form='UNFORMATTED', status='OLD')
        ifsec2(1)=3003         !sf.end
      endif
      ifsexcore(1) = iopenxx ('SEXcoreU')
      open(UNIT=3002,file='SEXcore2U',form='UNFORMATTED',status='OLD')
      ifsexcore2(1)=3002

      call readx (ifsex(1),50)
      read (ifsex(1),*) nspin,nqsig,ndsigm
      if(nspin == 2 .and. .not. laf) then
        ifsex(2)   = iopenxx ('SEXD')
        ifxc(2)    = iopenxx ('XCD')
        open(UNIT=3004, file='SEX2D',form='UNFORMATTED', status='OLD') !sf.beg
        ifsex2(2)=3004
        if(.not.exonly) then
          ifsec(2)   = iopenxx ('SECD')
          open(UNIT=3006, file='SEC2D',form='UNFORMATTED', status='OLD')
          ifsec2(2)=3006      !sf.end
        endif
        ifsexcore(2)   = iopenxx ('SEXcoreD')
        open(UNIT=3005,file='SEXcore2D',form='UNFORMATTED',status='OLD')
        ifsexcore2(2)=3005
      endif

      rewind (ifsex(1))

      ifqpe(1)   = iopenxx ('QPU')
      iftote(1)  = iopenxx ('TOTE.UP')
      iftote2(1) = iopenxx ('TOTE2.UP')

      allocate(eseavr(nqsig,nspin,2),eseavr_in(nqsig,nspin)) ! June 2009
      eseavr = 0

C ... Re-entry point for mode -1
  301 continue

C --- Read prior sigma into sigin, to pass to mixing routine ---
      INQUIRE (FILE='sigm', EXIST = lsigin)
      open(UNIT=3009, file='sigm',form='UNFORMATTED') !sf
      ifse_out=3009
      eseavr_file = 0
      if (lsigin) then
        rewind ifse_out
        print *, ' ... reading input sigma from file sigm'
C       read(ifse_out) nspin,ndimsigin,n1,n2,n3,nqx
        call iosigh(0,mode_sigfile,nspin,nspc,ndimh,ndimsigin,n1,n2,n3,nqx,is,
     .    lshft1,lshft2,lshft3,ifse_out,eseavr_file)

        if (jin == -1) nqsig = nqx

        if (nqx /= nqsig) then
          print 368, nqx,nqsig
 368      format (6x,' (warning) file mismatch : file nq=',i4,' but expected',i4)
          lsigin = .false.
        else

          if (jin == -1) then
!           allocate(eseavr(nqsig,nspin,2),eseavr_in(nqsig,nspin))
            allocate(eseavr_in(nqsig,nspin))
!           eseavr = NULLI
            izl = 0
          endif

          rewind ifse_out
          allocate(sigin(ndimsigin,ndimsigin,nqsig,nspin),qsigf(3,nqsig,nspin))
          call rwsigma(ifse_out,mode_sigfile,sigin,qsigf,nspin,ndimsigin,n1,n2,n3,nqsig,eseavr_in,eseavr_file)
          if (jin /= -1) deallocate(qsigf)  ! jin == -1 we won't re-read file, so keep qsigf

          print *, ' ... read successful: eseavr =', eseavr_file(1:nspin)
          do  is = 1, nspin
            if (eseavr_file(is) == NULLI) eseavr_file = 0
          enddo
        endif
      elseif (jin == -1) then
        call rx('no file sigm to read')
      endif
      rewind ifse_out

      if (nspin == 2 .and. jin /= -1) then
        ifqpe(2)   = iopenxx ('QPD')
        iftote(2)  = iopenxx ('TOTE.DN')
        iftote2(2) = iopenxx ('TOTE2.DN')
        !open(UNIT=3010, file='SED',form='UNFORMATTED')                !sf
        !ifse_out(2)=3010                                              !sf
      endif

      if(jin == -101) goto 9998

c-----------reading self-energy mode parameter from file 'GWIN_V2'
C     call readd_iSigma_en(ifinin,iSigma_en) !sf..21May02
      call getkeyvalue("GWinput","iSigMode",iSigma_en )
C      if(isigma_en==5) then     !.or.core3ptest) then
Cc     inquire(file='evec0',exist=evec0ex)
C        evec0ex=.false.    !true before 12Aug2006 ---> but it caused a problem maybe because of degeneracy.
C        if(evec0ex) then
C          ifevec0 =3018
C          open(ifevec0,file='evec0',form='UNFORMATTED',status='OLD')
C          ifevecchk=4019
C          if(isigma_en==5) open(ifevecchk,file='evecfix.chk')
C        endif
C      endif

C ... Read extension, number of threads, thread block array
      ifi=4001
      open(unit=ifi,file='mpiqlst',status='old')
      read(ifi,*) ext
      read(ifi,*) ext, ix, nqcalc, nrank
C     ix = fext(trim(ext))     ! Store file extension
      write(6,"(' hqpe_sc: files vxc',a,' and evec',a,' made for',i4,' qp by',i4,' thread(s)')")
     .  trim(ext),trim(ext),nqcalc,nrank
      allocate(iprocq(nqcalc))
      read(ifi,*) iprocq

C ... Open files generated by threads other than the master thread
      allocate(extp(0:nrank-1),ifvxc__(0:nrank-1),ifevec__(0:nrank-1),ifevecl__(0:nrank-1))
      extp(0) = trim(adjustl(ext))
      ix=5001
      do  procid = 0, nrank-1
        if (procid > 0) then  ! Make extension for MPI lmfgwd
          write(extn,"(i10)") procid
          write(extp(procid),"(a,i6)")
     .      trim(extp(0))//'_'//trim(adjustl(extn))
        endif
        ix=ix+1
        open(ix,file='vxc'//trim(extp(procid)),form='unformatted')
        ifvxc__(procid)=ix
        ix=ix+1
        open(ix,file='evec'//trim(extp(procid)),form='unformatted')
        ifevec__(procid)=ix
C        if (levec > 0) then
C          ix=ix+1
C          open(ix,file='evecl'//trim(extp(procid)),form='unformatted')
C          ifevecl__(procid)=ix
C          ifevecl_= ifevecl__(0)
C        endif
      enddo

      ifvxc_ = ifvxc__(0)
      ifevec_= ifevec__(0)
C     ndham probably shouldn't be read from vxc file, but for now ...
      read(ifvxc_) ndham,nsp,nqvxc
C     if(evec0ex) read(ifevec0)

C     evecs not needed when jin == -1
      if (jin /= -1) then

      rewind ifevec_
C     read(ifevec_) ndimx,nspx,nqx,iq,izl
      nspc = 1
C     ndimx should be the same as ndham.  izl=0 for now until issues with PMT are resolved.
      call iosigh(2,5,nspx,nspc,ndimx,nmto,n1,n2,n3,nqx,izl,lshft1,lshft2,lshft3,ifevec_,xx)
      izl = 0

      write(6,"(' sigma calculated with iSigMode =',i2,' for',i4,'  k-points and',i2,' spins')") iSigma_en,nqvxc,nsp
      write(6,"(' Maximum hamiltonian rank    = ',i4)") ndham
      write(6,"(' Maximum rank of <i|sigma|j> = ',i4)") ndsigm
      if (izl /= 0) write(6,"(' Use LMTO evecs to rotate sigm to orbital basis, dim=',i4)") izl

      !In PMT, hamiltonian rank ndimh is q-dependent. See lm*/gwd/sugw.F  june2009 takao
      if (nspin /= nsp)   call rx('nspin mismatch, vxc file')
      if (nspin /= nspx)  call rx('nspin mismatch, evec file')
      if (nqsig /= nqvxc) call rx('qp mismatch in vxc file')
      if (nqvxc /= nqx .or. ndimx < ndham)  call rx('mismatch in evec and vxc files')
      endif

      allocate(qsig(3,nqvxc,nspin))
      allocate(v_xc(ndham,ndham,nqvxc,nspin),evec(ndham,ndham,nqvxc,nspin))
C     if(evec0ex) allocate(evec00(ndham,ndham,nqvxc,nspin))
      allocate(ndimhq(nqvxc,nspin),ndzq(nqvxc,nspin))
      evecl => evec
      if (izl /= 0) allocate(evecl(izl,izl,nqvxc,nspin))

C ... Read eigenvectors and matrix elements of vxc
      write(6,"(' reading v_xc, evec ...'/'  iq  is  ndimh')")
      do  iq = 1, nqvxc
        do  is = 1, nspin
          ifvxc_  = ifvxc__ (iprocq(iq))
          ifevec_ = ifevec__(iprocq(iq))
          read(ifvxc_) ndimh,nmto ! ndimh is q-dependent hamiltonian rank (PMT)
                                  ! nmto  is rank of LMTO part of hamiltonian (independent of q)
          write(6,"(2i4,i6)") iq,is,ndimh
          ndimhq(iq,is) = ndimh   ! k-dependent dimension of hamiltonian
          read(ifvxc_) qsig(1:3,iq,is),v_xc(1:ndimh,1:ndimh,iq,is) ! qsig should be called qx here
          if (jin == -1) cycle
          if (izl /= 0) then
            call rx('hqpe_sc : izl > 0 branch needs to be revisited')
            read(ifevec_) qtmp
            read(ifevec_) ! eigenvalues
            read(ifevec_) evecl(:,:,iq,is)
            if (sum(abs( qsig(:,iq,is )-qtmp(:) )) > 1d-10)
     .        call rx('hqpe: mismatch qsig and q from evec file')
            ndzq(iq,is) = nmto
          else
            ndzq(iq,is) = ndimh
          endif
          read(ifevec_) qtmp
          read(ifevec_)         ! eigenvalues
          read(ifevec_) evec(:,:,iq,is)
          if (sum(abs( qsig(:,iq,is )-qtmp(:) )) > 1d-10)
     .      call rx('hqpe: mismatch qsig and q from evec file')
C         if(evec0ex) read(ifevec0) qtmp, evec00(1:ndimh,1:ndimh,iq,is)
        enddo
      enddo
      if (vxcscale /= 1d0) call dscal(2*size(v_xc),vxcscale,v_xc,1)

C     print *, '!!'; v_xc = 0

      if (mtosigmaonly()) then
        ndimsig = nmto
      else
        ndimsig = ndham
      endif
      if (lsigin .and. ndimsigin /= ndimsig) then
        deallocate(sigin)
        print *, '... input sigma dimension mismatch ... discarding'
        lsigin = .false.
      endif
      allocate(sigmv(ndimsig,ndimsig,nqsig))  ! should be merged with ev_se_ev
      if (jin /= -1 .and. izl /= 0 .and. ndimsig /= izl) call rx('mismatch sig, evec file')
C     call dpzero(sigmv,2*ndimsig*ndimsig*nqsig)
C ... Re-entry point for jin=-101
 9998 continue

      if (jin == -1) then
        sigma_m => sigin
        do  is = 1, nspin
          do  ip = 1, nqsig
C           Find v_xc(iq) corresponding to sigm(ip)
            iq = -1
            do ikp = 1, nqvxc
              if (sum ( (qsig(1:3,ikp,is)-qsigf(1:3,ip,1))**2 ) < 1d-10 ) then
                iq = ikp
                exit
              endif
            enddo
            if (iq == -1) call rx('hqpe.sc: no qp in vxc matches current qp')
            call dcopy(2*ndimh*ndimh,v_xc(1,1,iq,is),1,sigma_m(1,1,ip,is),1)
          enddo
        enddo
        eseavrmix = 0
        allocate(eseavr(nqsig,nspin,2)); eseavr = 0
        close(ifse_out)
        open(UNIT=ifse_out, file='vxcsig',form='UNFORMATTED')
        goto 3001
      endif

C --- For each spin, do ---
      do  1001  is = 1, nspin

        if(is==1) write(6,"(' reading XCU, SEXU, SEXcoreU, SECU data')")
        if(is==2) write(6,"(' reading XCD, SEXD, SEXcoreD, SECD data')")

        call readx(ifsex(is),50)
        read (ifsex(is),*) nspinx,nqx,nsigmx
        read (ifsex(is),*)
        read (ifsex(is),*) deltaw
        read (ifsex(is),*) alat
        read (ifsex(is),*) ef
        read (ifsex(is),*) esmr
        if (nspin /= nspinx)   call rx('nspin mismatch with SEX file')
        if (nqsig /= nqx)      call rx('qp mismatch with SEX file')
        if (ndsigm /= nsigmx)  call rx('ndsigm mismatch with SEX file')

        call readx(ifxc(is),50)
        read (ifxc(is),*) nspinx,nqx,nsigmx
        if (nspin /= nspinx)   call rx('nspin mismatch with XC file')
        if (nqsig /= nqx)      call rx('qp mismatch with XC file')
        if (ndsigm /= nsigmx)  call rx('ndsigm mismatch with XC file')

        read(ifsex2(is))  nspinx, nqx, nsigmx, nqbz,  nqibz,  n1,n2,n3 !sf.beg
        if (nspin /= nspinx)   call rx('nspin mismatch with SEX2 file')
        if (nqsig /= nqx)      call rx('qp mismatch with SEX2 file')
        if (ndsigm /= nsigmx)  call rx('ndsigm mismatch with SEX2 file')

        read(ifsexcore2(is)) nspinx,nqx,nsigmx,nqbzx,nqibzx
        if (nspin /= nspinx)   call rx('nspin mismatch with SEXcore2 file')
        if (nqsig /= nqx)      call rx('qp mismatch with SEXcore2 file')
        if (ndsigm /= nsigmx)  call rx('ndsigm mismatch with SEXcore2 file')
        if (nqbz /= nqbzx)    call rx('nqbz mismatch with SEXcore2 file')
        if (nqibz /= nqibzx)   call rx('nqibz mismatch with SEXcore2 file')

        if(.not.exonly)  then
          call readx(ifsec(is),50)
          read (ifsec(is),*) nspinx,nqx,nsigmx
          if (nspin /= nspinx) call rx('nspin mismatch with SEC file')
          if (nqsig /= nqx)    call rx('qp mismatch with SEC file')
          if (ndsigm /= nsigmx)call rx('ndsigm mismatch with SEC file')

          read(ifsec2(is))  nspinx, nqx, nsigmx, nqbzx, nqibzx,n1,n2,n3 !sf.end
          if (nspin /= nspinx) call rx('nspin mismatch with SEC2 file')
          if (nqsig /= nqx)    call rx('qp mismatch with SEC2 file')
          if (ndsigm /= nsigmx)call rx('ndsigm mismatch with SEC2 file')
          if (nqbz /= nqbzx)   call rx('nqbz mismatch with SEC2 file')
          if (nqibz /= nqibzx) call rx('nqibz mismatch with SEC2 file')
        endif

C   ... read XC data
        allocate( itxc(ndsigm),qxc(3,ndsigm,nqsig),eldaxc(ndsigm,nqsig),vxc(ndsigm,nqsig) )
        call readx(ifxc(is),50)
        read(ifxc(is),*)
        do iq = 1,nqsig
          do i = 1,ndsigm
            read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itxc(i),ix,ix, qxc(1:3,i,iq), eldaxc(i,iq), vxc(i,iq)
          enddo
        enddo

C   ... read SEX data
        allocate( itx(ndsigm), qx (3,ndsigm,nqsig),eldax (ndsigm,nqsig),sex(ndsigm,nqsig) )
        allocate( sex2(ndsigm,ndsigm,nqsig)) !sf..3June
        call readx(ifsex(is),50)
        read(ifsex(is),*)
        do iq = 1,nqsig
          do i  = 1,ndsigm
            read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itx(i),ix,ix, qx(1:3,i,iq), eldax(i,iq),
     &              sex(i,iq)
          enddo
          read(ifsex2(is)) isx,qx2,sex2(1:ndsigm,1:ndsigm,iq) !sf..3June
        enddo

C   ... read SEXcore data
        allocate( sexcore(ndsigm,nqsig),sexcore2(ndsigm,ndsigm,nqsig) ) !sf..3June
C       if(exonly) then         !zero for exonly case
C         print *,' set sexcore=0 for exonly case'
C         sexcore=0d0; sexcore2=0d0
C       else
        call readx   (ifsexcore(is),50)
        call readx   (ifsexcore(is),50)
        read(ifsexcore(is),*)
        do iq = 1,nqsig
          do i  = 1,ndsigm
            read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              ix,ix,ix, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,iq)
          enddo
          read(ifsexcore2(is)) isx,qx2,sexcore2(1:ndsigm,1:ndsigm,iq) !sf..3June
        enddo
C     endif

C   ... read SEC data
        allocate( itc(ndsigm), qc (3,ndsigm,nqsig),eldac (ndsigm,nqsig)
     &        ,rsec(ndsigm,nqsig),csec(ndsigm,nqsig),sec2(ndsigm,ndsigm,nqsig)) !sf..3June
        if(exonly) then        !zero for exonly case
          print *,' set sec=0 for exonly case'
          itc=0d0;qc=0d0;eldac=0d0;rsec=0d0;csec=0d0;sec2=0d0
        else
          call readx(ifsec(is),50)
          read(ifsec(is),*)
          do iq = 1,nqsig
            do i  = 1,ndsigm
              read(ifsec(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16,
     &                 3x,d24.16)") !sf13May2002     ,3x,d24.16)")
     &                 itc(i),ix,ix, qc(1:3,i,iq), eldac(i,iq),
     &                 rsec(i,iq),csec(i,iq) !sf..13May2002   ,zfac(i,iq)
            enddo
            read(ifsec2(is)) isx,qx2,sec2(1:ndsigm,1:ndsigm,iq) !sf..3June
          enddo
        endif

        if(.not.exonly) ix=icompvv2(itx, itc,ndsigm,' state-label x c')
        if(.not.exonly) ix=icompvv2(itx,itxc,ndsigm,' state-label x xc')

C   ... Calculate quasiparticle energies and write them to file QPE
        call qpe1_sc (ifqpe(is),iftote(is),iftote2(is),itx,qx, !itc,qc
     i        eldax,vxc,sex,sexcore,
     i        rsec,csec,jin,deltaw,alat,ef, !sf..13May |it was zfac after csec
     d        ndsigm,nqsig,is,
     o        eshift0,eshift02,eshtlda)
c     o              w(ieqp),w(iwqp))
        close(ifqpe(is))
        close(iftote(is))
        close(iftote2(is))
        if(jin==-101) cycle

C   --- For each q, assemble SE_ij.  ij are indices in the eigenfunction basis
        allocate(ipiv(ndham),work(ndham*ndham),evec_inv(ndham,ndham),evec_invt(ndham,ndham))
        allocate(se(ndham,ndham,nqsig),ev_se_ev(ndimsig,ndimsig))
        se = 0
C       if(evec0ex) allocate(evec00inv(ndham,ndham))
        do  iq = 1, nqsig
          do itp = 1, ndsigm
            do itpp = 1, ndsigm    !make Sigma hermitian
              if (.not.exonly) then
                se(itpp,itp,iq) = sex2(itpp,itp,iq)+sexcore2(itpp,itp,iq)
     &                   +.5d0*(sec2(itpp,itp,iq)+dconjg(sec2(itp,itpp,iq)))
              else
                se(itpp,itp,iq)= wex*sex2(itpp,itp,iq)+sexcore2(itpp,itp,iq)
              endif
cccccccccccccccccccccccc
c     c takao test ex only test
c     se(itpp,itp,iq)= 0.3*sex2(itpp,itp,iq)+ sexcore2(itpp,itp,iq)
ccccccccccccccc
c     c takao test2---reduce screen
c     se(itpp,itp,iq)=sex2(itpp,itp,iq)+sexcore2(itpp,itp,iq)+
c     &        .1d0*(sec2(itpp,itp,iq)+dconjg(sec2(itp,itpp,iq)))
cccccccccccccc
            enddo
          enddo
        enddo

!   ... Check that Sigma is hermitian
        do iq=1,nqsig
          do itp=1,ndsigm
            do itpp=1,ndsigm
              if(abs(se(itpp,itp,iq)-dconjg(se(itp,itpp,iq)) ) > 1d-8) then
                write(*,*)'diff=', se(itpp,itp,iq)-dconjg(se(itp,itpp,iq))
                write(*,*)'se   ', se(itpp,itp,iq), se(itp,itpp,iq)
                write(*,*)'sex2 ',sex2(itpp,itp,iq),sex2(itp,itpp,iq)
                write(*,*)'sexc2',sexcore2(itpp,itp,iq),sexcore2(itp,itpp,iq)
                call rx( "hqpe: Sigma_nn' is not hermitian")
              endif
            enddo
          enddo
        enddo

C   --- Make QSGW sigma-vxc in basis of eigenfunctions for all kp ---

C   ... Find k-index(vxc) associated with k-index(sigma)
C       ip is k-index to self-energy, iq = k-index to v_xc and evec
C       Note: qsig and qx should be renamed. Now qsig=q(v_xc), and qx=q(sigma)
        allocate(ndsigq(nqsig)); neseavr = 0
        do  2001 ip = 1, nqsig
C         find which iq = index to the v_xc array which corresponds to ip
          do ikp = 1, nqvxc
            if (sum ( (qsig(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < 1d-10 ) then
              iq=ikp      !qc(:,i,:) does not depend on band index i=1:ndsigm
              goto 100
            endif
          enddo               !ikp
          call rx('hqpe.sc: no qp in vxc matches current qp')
  100     continue
          ndimh = ndimhq(iq,is)      !june 2009 takao. we use ndimh instead of ndham

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     takao test only diagonal
c     call onlydiag(se(1:ndsigm,1:ndsigm,ip),ndsigm)
c     do itp=1,ndsigm
c     se(itp,itp,ip)=se(itp,itp,ip) - 0.5*vxc(itp,ip)/rydberg()
c     enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     ... Get ndsigq = number of states for which se is calculated
          ndsigq(ip)=ndsigm
          do  itp = ndsigm, 1, -1
            if (se(itp,itp,ip)/=0d0) then
              ndsigq(ip) = itp
              exit
            endif
          enddo

C          if(ip==1) then  ! Not actually used
C            ndsigqmin=ndsigq(ip)
C          else
C            if(ndsigq(ip)<ndsigqmin) ndsigqmin=ndsigq(ip)
C          endif

          fac = -0.5d0; if (exonly) fac = -wex
          if (isigma_en == 5) then
            do  itp = 1, ndsigq(ip)
              do  itpp = 1, ndsigq(ip)
                if (itp/=itpp) then
                  se(itp,itpp,ip) = 0d0
                else
                  se(itp,itpp,ip) = se(itp,itpp,ip)
     .              + fac * sum(dconjg(evec(1:ndimh,itp,iq,is))*
     .              matmul(v_xc(1:ndimh,1:ndimh,iq,is),evec(1:ndimh,itpp,iq,is)))
                endif
              enddo
            enddo

C         Add fac*(z+ v_xc z) into se.  Then se = sigm-vxc in eigenfunction basis
          else
            allocate(zwk(ndimh,ndsigq(ip)))
C           call zprm('evec',evec(1,1,iq,is),ndham,ndimh,ndsigq(ip))
C           call zprm('v_xc before rotation',v_xc(1,1,iq,is),ndham,ndimh,ndimh)
C           call zprm('se before vxc',se(1,1,ip),ndham,ndsigq(ip),ndsigq(ip))
            call zgemm('N','N',ndimh,ndsigq(ip),ndimh,one,v_xc(1,1,iq,is),ndham,
     .        evec(1,1,iq,is),ndham,zer,zwk,ndimh)
C            call zprm('zwk',zwk,ndimh,ndimh,ndsigq(ip))
            call zgemm('C','N',ndsigq(ip),ndsigq(ip),ndimh,dcmplx(fac),evec(1,1,iq,is),ndham,
     .        zwk,ndimh,one,se(1,1,ip),ndham)
C           call zprm('se after vxc subtracted',se(1,1,ip),ndham,ndsigq(ip),ndsigq(ip))
            deallocate(zwk)
          endif

ctakao 2009June. Find average diagonal of se above a threshold
c      eldaxc(ndsigm,nqsig)
!     Takao June 2014 New definition of average (eseavr)
!     e-weighted average
!     eavr  : average of eigenvalues within threshold (itp<=ndsigq)
!     eavr2  : square average of eigenvalues within threshold (itp<=ndsigq)
!     eseavr: average of se*eigenvalue
          eavr  = 0d0; eavr2  = 0d0; eseavr0= 0d0; eseavr02 =0d0
          threshold = 5d0  ! in eV. Include states at ef+threshold and above
          ix=0
          do  itp = 1, ndsigq(ip)
C           eee = eldax(itp,ip) - rydberg()*ef
            eee = eldax(itp,ip)  ! Bug fix MvS 13 Nov 2013. eldax(itp,ip) already relative to ef (in eV)
            if (eee > threshold) then
              if (ix == 0) ix = itp
              neseavr = neseavr+1
              eavr = eavr + eee
              eavr2   = eavr2   + eee**2
              eseavr0 = eseavr0 + eee*se(itp,itp,ip)
              eseavr02 = eseavr02 + eee**2*se(itp,itp,ip)
            endif
          enddo
          if (eavr == 0) then
             eavr=1d0; eseavr0=0d0
             eavr2=1d0; eseavr02=0d0
            call rx(' hqpe_sc (fatal): no evals within threshold')
          endif
          eseavr(ip,is,1) = eseavr0/eavr !energy-weighted average self-energy above threshold
          eseavr(ip,is,2) = eseavr02/eavr2 !new average
          write(6,"( ' qp',i5,' spin',i4,' sigma calculated from',i5,' to',i5,' of',i5,' eseavr(1:2) =',2f12.6)")
     .      ip,is,ix,ndsigq(ip),ndsigm,2*eseavr(ip,is,:)
          if (IAND(mod(GWversion()/100,10),1) == 1) eseavr(ip,is,1) = eseavr(ip,is,2)
C          do itp = ndsigq(ip)+1,ndsigm  ! Sep12 style: eseavr is k-dependent
C            se(itp,itp,ip) = eseavr(ip,is,1)
C          enddo
          eseavr(ip,is,:) = 2d0*eseavr(ip,is,:) !Hartree to  Ry
 2001   continue

C   --- k-integrated mean value of eseavr ---
        if (is == 1) then
          call read_BZDATA()
          write(*,"(' read from bzdata : nqibz2 nqibz nqsig',3i4)") nqibz2,nqibz,nqsig
          nstarsum= sum(nstar(:))
          if (nstarsum/= nqbz ) call rx(' nstarsum/= nqbz')
C         if_eseavr = iopenxx ('ESEAVR')
          eseavrmean = 0d0
        endif
        do ip=1,nqsig
          eseavrmean(is) = eseavrmean(is) + nstar(ip)*eseavr(ip,is,1)
        enddo
        eseavrmean(is) = eseavrmean(is)/nqbz
        call getkeyvalue("GWinput","AddToESEAVR",eseadd,default=0d0,status=ret)
        eseavrmean(is) =eseavrmean(is)+eseadd
C       write(if_eseavr,"(d23.15,i3,i8)") eseavrmean(is),is,ndsigqmin
        if (laf) then
          eseavrmean(2)= eseavrmean(1)
C         write(if_eseavr,"(d23.15,i3)") eseavrmean(2),2,ndsigqmin
          exit
        endif
C       write(if_eseavr,"(a)")' e-weighted average of self-energy for unoccupied bands. ispin'
C       close(if_eseavr)

C   --- Replace sigm with diagonal part for states above ndsigq ---
        write(6,"( ' k-averaged eseavr = ',1f10.6,' calculated from',f8.3,' states per qp on average')")
     .    eseavrmean(is),dble(neseavr)/nqsig
        do  2002  ip = 1, nqsig
          do  itp = ndsigq(ip)+1, ndimh
            se(itp,itp,ip) = eseavr(ip,is,1)/2 ! September 2012
            if (IAND(mod(GWversion()/100,10),1) == 1)
     .      se(itp,itp,ip) = eseavrmean(is)/2
          enddo
 2002   continue

!   --- Write header to sigm file ---
        if (is == 1) then
C         write(ifse_out) nspin,ndimsig,n1,n2,n3,nqsig,0,0,0
          call iosigh(0,0,nspin,1,ndimsig,ndimsig,n1,n2,n3,nqsig,nqsig,0,0,0,-ifse_out,eseavr)
        endif

C   ... Possibly read new set eigenvectors for rotating sigm
C        evecl => evec
C        if (levec > 0) then
C          ifevec_ = ifevecl__(0)
C          read(ifevec_) ndimx,nspx,nqx
C          do iq = 1, nqvxc
C            ifevec_ = ifevecl__(iprocq(iq))
C            read(ifevec_) qtmp(1:3),evecl(1:nmto,1:nmto,iq,is)
C            if (sum(abs( qsig(:,iq,is )-qtmp(:) )) > 1d-10)
C     .        call rx( 'hqpe: mismatch qsig and q from evecl file')
C          enddo
C          ndzq(:,is) = nmto
C        else
C          ndzq = ndimhq
C        endif

C   --- Project modified sigm onto LMTO part of orbital basis and write to disk ---
        do 2003 ip=1, nqsig
          do ikp=1, nqvxc
            if (sum ((qsig(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < 1d-10 ) then
              iq=ikp          !qc(:,i,:) does not depend on band index i=1:ndsigm
              goto 102
            endif
          enddo                 !ikp
          call rx( 'hqpe.sc: not find ikp 102')
  102     continue
          ndimh = ndimhq(iq,is) !full hamiltonian rank for this q
          ndz   = ndzq(iq,is)   !number of eigenvectors to make sigm for this q
          if(mtosigmaonly()) then ! rank of sigm to retain for this q
            ndimsig = nmto
          else
            ndimsig = ndz
          endif

!!        Make inverse evec_inv(n,i) of eigenvector matrix evecl(i,n)
!!        \psi_n = sum_i evecl(i,n)\phi_i, \psi is eigenfunction and \phi is basis function
          evec_inv(1:ndz,1:ndz) = evecl(1:ndz,1:ndz,iq,is) !evecl(i,n), i-basis, n-band argumants
          call zgetrf(ndz,ndz,evec_inv,ndham,ipiv,info)   !ndham is leading dimension of evec_inv
          if (info /= 0)  call rx('hqpe: zgetrf info/=0')
!         inverse matrix evec_inv(n,i), i-basis, n-band argumants
          call zgetri(ndz,evec_inv,ndham,ipiv,work,ndham**2,info)
          if (info /= 0)  call rx('hqpe: zgetri info/=0')
c         Takao: gfortran did not allow next line.  So I made a replacement.
c         evec_invt= dconjg(transpose(evec_inv))  !this did not work for gfortran
C         evec_invt= transpose(dconjg(evec_inv)) !Can fail for large systems
C         call zprm('evec',evecl,ndham,ndz,ndz)
C         call zprm('evec_inv',evec_inv,ndham,ndz,ndz)
          do  itp = 1, ndz
            do  itpp = 1, ndz
              evec_invt(itp,itpp) = dconjg(evec_inv(itpp,itp))
            enddo
          enddo

C          if(evec0ex.and.iSigma_en==5) then
C            call rx('evec0ex.and.iSigma_en==5 no longer supported')
C            evec00inv(1:ndz,1:ndz)=evec00(1:ndz,1:ndz,iq,is)
C            !evecl(i,n), i-basis, n - band argumants
C            call zgetrf(ndz,ndz,evec00inv,ndham,ipiv,info)
C            if (info /= 0)  call rx('hqpe: zgetrf 00 info/=0')
C            call zgetri(ndz,evec00inv,ndham,ipiv,work,ndham**2,info)
C            !inverse matrix evec_inv(n,i), i-basis, n - band argumants
C            if (info /= 0)  call rx('hqpe: zgetri 00 info/=0')
C            itppxs=0
C            do itp=1,ndz
C              sum0=0d0
C              sum1=0d0
C              do itpp=1,ndz
Cc     print *,evec_inv(itp,itpp)
C                aaa=abs(sum(evec_inv(itp,1:ndz)*evec00(1:ndz,itpp,iq,is)))**2
C                if(aaa>sum1) then
C                  sum1 =aaa
C                  itppx=itpp
C                endif
C                sum0=sum0+aaa
C              enddo
C              evec_inv(itp,1:ndz)    =  evec00inv(itppx,1:ndz)
C              evecl(1:ndz,itp,iq,is)=  evec00(1:ndz,itppx,iq,is)
C              itppxs = itppxs+itppx
C              write(ifevecchk,"(3i5,2d13.5)")itp,iq,is, sum0, sum1
Cc     if(abs(sum0-1d0)>0.1) stop "!!! sum0 is not unity ==="
C            enddo
C            if(ndz*(ndz+1)/2 /= itppxs) then
C              print *,'ndz_sum itppxs=', ndz*(ndz+1)/2,itppxs
C              stop "hqpe: ndz*(ndz+1)/2 /= itppxs"
C            endif
C          endif

C     ... Make ev_se_ev = sigm in orbital basis = [z+]^-1 sigm z^-1
C         contribution from the first ndz evecs
C!         Inline code
C          nevec = ndsigm         ! All states
C          nevec = ndsigq(ip)     ! Include states only to ndsigq
C!         subset of z^-1 that extends to ndsigq only
C          ev_se_ev(1:ndz,1:ndz) = matmul(evec_invt(1:ndz,1:nevec)
C     &    ,matmul(se(1:nevec,1:nevec,ip),evec_inv(1:nevec,1:ndz)))
CC          call zprm('ev_se_ev without high lying additions',ev_se_ev,ndz,ndz,ndz)
C!     ... for ndsigq+1..ndsigm sigm is replaced by a constant, diagonal array, rotated to LMTO
C          exxx = eseavrmean(is); if (IAND(mod(GWversion()/100,10),1) == 1) exxx = eseavr(ip,is,1)
C          do itp =1,ndz
C          do itpp=1,ndz
C            ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) +
C     &      sum(evec_invt(itp,ndsigq(ip)+1:ndz)*evec_inv(ndsigq(ip)+1:ndz,itpp))*exxx/2d0 ! in Hartree.
C          enddo
C          enddo
CC         call zprm('ev_se_ev',2*ev_se_ev,ndz,ndz,ndz)

C     ... Vectorized ev_se_ev = SE-Vxc in orbital basis = [z+]^-1 sigm z^-1
          allocate(zwk(ndz,ndimsig))
C         call zprm('se(eigenfunction basis, Ry)',2*se(:,:,ip),ndham,ndimh,ndimh)
C         call zprm('evec_inv',evec_inv,ndham,ndz,ndz)
          call zgemm('N','N',ndz,ndimsig,ndz,(1d0,0d0),
     .      se(1,1,ip),ndham,evec_inv,ndham,(0d0,0d0),zwk,ndz)
          call zgemm('N','N',ndimsig,ndimsig,ndz,(1d0,0d0),
     .      evec_invt,ndham,zwk,ndz,(0d0,0d0),ev_se_ev,ndimsig)
          deallocate(zwk)
C         call zprm('se(orbital basis, Ry)',2*ev_se_ev,ndimsig,ndimsig,ndimsig)

!     ... write SE_ij-Vxc_ij, orbital basis (ij are orbital basis indices)
          do  itp = 1, ndimsig
            do  itpp = 1, ndimsig
              if(abs(ev_se_ev(itpp,itp) - dconjg(ev_se_ev(itp,itpp))) > 1d-8 ) then
                print *,itp,itpp
                print *,ev_se_ev(itpp,itp),ev_se_ev(itp,itpp)
                call rx( 'hqpe: Sigma_ij is not hermitian')
              endif
              if(abs(v_xc(itp,itpp,iq,is)-dconjg(v_xc(itpp,itp,iq,is)))> 1d-8)
     &          call rx( 'hqpe: v_xc is not hermitian')
            enddo
          enddo

C         Written here, to be read again at the mixing stage.  This should be cleaned up!
          write(ifse_out) qsig(1:3,iq,is),eseavr(ip,is,1),ndsigq(ip)
          if (mtosigmaonly()) then
            ev_se_ev = 2*ev_se_ev !in Ry.
            write(ifse_out) ev_se_ev
          else
            sigmv(:,:,ip) = 1d20
            sigmv(1:ndimsig,1:ndimsig,ip) =  2d0*ev_se_ev(1:ndimsig,1:ndimsig) !in Ry.
            write(ifse_out) sigmv(:,:,ip)
          endif
 2003   continue
        deallocate(ndsigq)
        write(6,*)

!!     - ANFcond case. We will make this an independent program ---
        if(laf) then
          print *, ' Find ANFcond --- Call anfsig.'
          ildima = 211
          open(ildima,file='ldima')
          read(ildima,'(a3)') iaaa
          close(ildima)

C     ... June2007 for floating orbitals
!     ldima is generaged by lmfgw; it contains number of MTO including
!     floating orbital, and positions.
          open(ildima,file='ldima')
          if(iaaa=='***') then
            read(ildima,*) iaaa,natom
            allocate( ibasf(natom),ldima(natom),pos(1:3,natom) )
            do ibas = 1,natom
              read(ildima,*) ldima(ibas),pos(1:3,ibas)
              write(6,"('ldima pos=',i5,3f10.4)")ldima(ibas),pos(1:3,ibas)
            enddo
          else
            ilmto = 212
            open(ilmto,file='LMTO')
            do ix=1,12
              read(ilmto,*)
            enddo
            read(ilmto,*)natom
            allocate(ibasf(natom),ldima(natom)
     &              ,pos(1:3,natom))
            write(6,"(' natom=',i3)") natom
            read(ilmto,*)
            do ibas = 1,natom
              read(ilmto,*) pos(1:3,ibas)
              write(6,*) pos(1:3,ibas)
            enddo
            close(ilmto)

            do ibas=1,natom
              read(ildima,*) ldima(ibas)
              write(6,*) 'ldima=',ldima(ibas)
            enddo
          endif
          close(ildima)

          ifianf = 211
          open(ifianf,file='ANFcond')
          read(ifianf,*)
          read(ifianf,*)
          read(ifianf,*) anfvec(1:3)
          write(6,'(" ANFcond=",3f13.6)') anfvec
          do ibas=1,natom
            read(ifianf,*) ibas1, ibasf(ibas)
            if(ibas1/=ibas) call rx( ' readin ANFcond error:ibas1/=ibas')
            write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
          enddo
          close(ifianf)

          do ip = 1,nqsig
c     do ikp=1,nqbz
c     if (sum ( (qsig(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < 1d-10 ) then
c     iq=ikp !qx(:,i,:) does not depend on band index i=1:ndsigm
c     goto 101
c     endif
c     enddo !ikp
C           sigmv(:,:,ip)=1d20
            call anfsig( natom, ldima, pos, ndz,
     i              qx(1:3,1,ip),
     i              anfvec, ibasf, ! these are antiferro informations.
     i              sigmv(1:ndz,1:ndz,ip)) ! i/o
            write(ifse_out) qx(1:3,1,ip),eseavr(ip,is,1) !june2009 add eseavr
            write(ifse_out) sigmv(:,:,ip)
          enddo
        endif
        deallocate(sex2,sexcore2,sec2,se,ipiv,work)
        deallocate(evec_inv,evec_invt,ev_se_ev) !,se_ev
C       if(evec0ex) deallocate(evec00inv)
C     - end making SE_ij-VXC_ij cccccccccccccccccccccccccccccccccc  !sf..3June
        deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,eldac,
     &        sexcore ,rsec,csec, !zfac,  !sf..13May2002
     &        itx, qx ,eldax,sex)
        if (laf) exit
        if (jin > 0) jin = 999999
 1001 continue !isp


      if(jin==-101) goto 9999
      deallocate(v_xc,evec,sigmv)
      close(ifse_out)

C --- Mix sigma with prior iteration and rewrite ---
c     eseavr wasn't written to file, so read it to dummy.
      open(UNIT=ifse_out, file='sigm',form='UNFORMATTED')
C     Re-read sigm just created (should be cleaned up)
      rewind ifse_out
      allocate(sigma_m(ndimsig,ndimsig,nqsig,nspin),qsigf(3,nqsig,nspin))
      print *,"========= Mix Sigma with prior iterations read from mixsigma file ======="
      call rwsigma(ifse_out,mode_sigfile,sigma_m,qsigf,nspin,ndimsig,n1,n2,n3,nqsig,eseavr_in,xx)

C     MvS Feb 04 Allocate sigin if not already done so
      if (lsigin .and. ndimsigin /= ndimsig) then
        deallocate(sigin)
        print *, '... input sigma dimension mismatch ... discarding'
        lsigin = .false.
      endif
      if (.not. lsigin) then
        allocate(sigin(1,1,1,1))
        sigin = 0d0
      endif
      call mixsigma(sigma_m, lsigin, sigin, 2*ndimsig**2*nqsig*nspin, beta)
      eseavrmix = beta*eseavrmean + (1-beta)*eseavr_file

C ... Write sigm Binary format
 3001 continue  ! Re-entry point for jin == -1
      if (jin == -1) then
        mode_sigfile = 8
        write(6,"(/' === Write vxc to file vxcsig ===')")
      else
        mode_sigfile = 0
        write(6,"(/' === Write sigma-vxc to file sigm ===')")
      endif
      rewind ifse_out
      call rwsigma(-ifse_out,mode_sigfile,sigma_m,qsigf, nspin,ndimsig,n1,n2,n3,nqsig,eseavr,eseavrmix)
      close(ifse_out)

C     ... Write sigma ASCII format
C      open(UNIT=ifse_out, file='sigma')
C      call rwsigmaa (-ifse_out,sigma_m,qsigf, nspin,ndimsig,n1,n2,n3,nqsig,eseavr,
C     .  eseavrmix)
C      close(ifse_out)
 9999 call rx0(' OK! hqpe_sc ')
      end

c------------------------------------------------------
      subroutine testfff(a,b,nnn)
      real(8):: a(nnn),b(nnn)
      do i=1,nnn
        if(i/=nnn) then
          b(i)= +  (-1d0)**i*a(i+1) + i/10d0
        else
          b(i)=   a(i)**2 - 1d0
        endif
      enddo
      end

c------------------------------------------------------
      subroutine onlydiag(sigm,ndham)
      complex(8):: sigm(ndham,ndham)
      do ix=1,ndham
        do iy=1,ndham
          if(ix/=iy) sigm(ix,iy)=0d0
        enddo
      enddo
      end

C----------------------------------------------------------------------
      subroutine rwsigma(ifs,mode_sigfile,sigm,qsigf,nspin,ndsigm,n1,n2,n3,nq,eseavr,eseavrmean)
      implicit none
      integer ifs,mode_sigfile,n1,n2,n3,nq,nspin,ndsigm
      complex(8)::sigm(ndsigm,ndsigm,nq,nspin)
      real(8)::qsigf(3,nq,nspin),eseavr(nq,nspin),eseavrmean(nspin)
      integer :: ntm1=0, ntm2=0, ntm3=0, is,ip,nspc,ndham

C     For now
      nspc = 1
      ndham = ndsigm

      if(ifs>0) print *, " rwsigma: Reading sigm (Binary)"
      if(ifs<0) print *, " rwsigma  Writing sigm (Binary)"
      call iosigh(0,mode_sigfile,nspin,nspc,ndham,ndsigm,n1,n2,n3,nq,is,
     .  ntm1,ntm2,ntm3,ifs,eseavrmean)
      do is=1,nspin
        do ip=1,nq
          if(ifs>0) read (ifs)  qsigf(1:3,ip,is),eseavr(ip,is)
          if(ifs>0) read (ifs)  sigm(1:ndsigm,1:ndsigm,ip,is)
          if(ifs<0) write(-ifs) qsigf(1:3,ip,is),eseavr(ip,is)
          if(ifs<0) write(-ifs) sigm(1:ndsigm,1:ndsigm,ip,is)
          write(6,"('  === ',i5,i3,3f10.5)")ip,is,qsigf(1:3,ip,is)
        enddo
      enddo
      print *,' === rwsigma:  sum check of sigm=',sum(abs(sigm))
      print *
      end
c
      subroutine rwsigmaa(ifs,sigm,qsigf, nspin,ndsigm,n1,n2,n3,nq,eseavr,eseavrmean)
C Write without version in header, for compatibility with Takao
      complex(8)::sigm(ndsigm,ndsigm,nq,nspin)
      real(8)::qsigf(3,nq,nspin),eseavr(nq,nspin),eseavrmean(nspin)
      ntm1=0;ntm2=0;ntm3=0
      if(ifs>0) print *, " rwsigmaa: Reading sigma (ASCII)==="
      if(ifs<0) print *, " rwsigmaa: Writing sigma (ASCII)==="
      if(ifs>0)read ( ifs,*)      nspin,ndsigm,n1,n2,n3,nq,ntm1,ntm2,ntm3
      if(ifs<0)write(-ifs,"(9i10,2d24.16)")nspin,ndsigm,n1,n2,n3,nq,ntm1,ntm2,ntm3,eseavrmean
      do is=1,nspin
        do ip=1,nq
          if(ifs>0) read (ifs,*)  qsigf(1:3,ip,is),eseavr(ip,is)
          if(ifs<0) write(-ifs,"(3d24.16)") qsigf(1:3,ip,is),eseavr(ip,is)
          do j=1,ndsigm
            do i=1,ndsigm
              if(ifs>0) read (ifs,"(2d24.16)")  sigm(i,j,ip,is)
              if(ifs<0) write(-ifs,"(2d24.16,2i5)") sigm(i,j,ip,is),i,j
            enddo
          enddo
          write(6,"('  === ',i5,i3,3f10.5)")ip,is,qsigf(1:3,ip,is)
        enddo
      enddo
      print *,' === rwsigmaa: sum check of sigm=',sum(abs(sigm))
      print *
      end

C----------------------------------------------------------------------
      subroutine mixsigma(sss, lsigin, sigin, nda, beta)
c  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
C- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
C- Anderson mixing of a vector
Ci  mmix: number of iterates available to mix
Cio nmix: nmix > 0: number of iter to try and mix
Ci        nmix < 0: use mmix instead of nmix.
Co  nmix: (abs)  number of iter actually mixed.
Co        (sign) <0, intended that caller update nmix for next call.
Co  beta:  mixing beta
C  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
C             Add mixnit as parameter
      use keyvalue
      implicit none
      logical lsigin
      integer nda,mmix
      integer(4),parameter:: mxsav=3
      double precision rms2,tj(mxsav),beta
      integer imix,amix
      integer iprintxx,ifi,nitr,ndaf
      real(8)::sss(nda),sigin(nda)
      real(8):: tjmax

      real(8),allocatable::norm(:),a(:,:,:)
      integer(4),allocatable:: kpvt(:)

c      logical :: readgwinput,ret
      integer(4)::ret

c heap
c      integer(4),parameter::iwksize=10*1000*1000
c      integer(4):: w
c      common /w/ w(iwksize)
      character*8 fff
      logical fexist

      real(8):: acc
      integer(4):: ido

c      call wkinit(iwksize)
      iprintxx = 30

      beta=1d0
C      if(readgwinput()) then
      call getkeyvalue("GWinput","mixbeta",beta,default=1d0,status=ret)
C      endif
c      fff="MIX"
c      INQUIRE (FILE =fff, EXIST = fexist)
c      if(fexist) then
c        open(ifi,file=fff,form='formatted')
c        read(ifi,*) beta
c        close(ifi)
c      endif

      print '('' mixsigma: Anderson mixing sigma with mixing beta ='',f12.6)',beta

C' ... reads prior iteration INCLUDING starting sigma for current iteration
      allocate ( a(nda,0:mxsav+1,2) )
      a = 0
      fff="mixsigma"
      INQUIRE (FILE =fff, EXIST = fexist)
      if(fexist)      print *,'... reading file mixsigma'
      if(.not.fexist) print *,'... No file mixsigma'
      ifi=3021
      open(ifi,file=fff,form='unformatted')
      if(fexist) then
        read(ifi,err=903,end=903) nitr,ndaf
        if (ndaf /= nda) goto 903
C       print *,' read mixsigma: nitr=',nitr
        read(ifi,err=903,end=903) a
C       print *,' read end of mixsigma'
        goto 902
      endif
      goto 901

  903 continue
      print 368
  368 format(5x,'(warning) file mismatch ... mixing file not read')
  901 continue
      nitr = 0
C      goto 902
C      rewind(ifi)
C      print *,'--- New mixsigma  --- '
C      print *,'--- writing sigma to mixsigma  --- '
C      print *
C      a = 0d0
C      a(:,0,2) = sss
C      nitr = 1
C      write(ifi) nitr,nda
C      write(ifi) a
C      close(ifi)
C      return

 902  continue
c      beta     =  .5d0
c      if(nitr>10)  beta=1d0
C     beta=1d0

      a(:,0,1) = sss      !output
C     if input sigma available, use it instead of file a(:,0,2)
      if (lsigin) then
        print *,'... using input sigma read from sigm file'
        a(:,0,2) = sigin  !input
      endif

      print *,'sum sss=',sum(abs(sss))

C     Restrict maximum number of prior iterations
      call getkeyvalue("GWinput","mixpriorit",imix,default=9,status=ret)
      mmix = min(max(nitr-1,0),imix)

      if (mmix > mxsav) mmix = mxsav

c      nmix = mmix
c      print *,'mixsigma: 3333'
C ... default nmix
c      if (nmix < 0) nmix = 2
c      nmix = min(mmix,nmix)

c      if (beta < 0) nmix = -nmix
C ... imix is a local copy of nmix
c      imix = nmix

c      if (imix < 0) imix = mmix
C ... save PQ array
c      call defdr(oa,nda*(mxsav+2)*2)
c      call dcopy(nda*(mxsav+2)*2,a,1,w(oa),1)

c      print *,'mixsigma: goto daxpy'
C ... for iterations 1,2,... amix needs F-x
c    1 continue
c      do  10  jmix = 1, nmix
c   10 call daxpy(nda,-1d0,a(1,jmix,2),1,a(1,jmix,1),1)
ccc      call query('beta',4,beta)
c      jmix = min(mmix,iabs(imix))

C     print *,'mixsigma: goto amix'
c      call defdr(onorm,mxsav**2)
c      call defi (okpvt,mxsav)
      allocate(norm(mxsav**2),kpvt(mxsav))
C     this information already printed out by amix
C     print *,'mixing parameters for amix are fixed in mixsigma'
C     print *,'   beta       =', beta
C     print *,'   tjmax      =', tjmax
C     print *,'   mmix mxsav =', mmix,mxsav

C     iprintxx = 61

      call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
      if(acc/=0d0) then
        print *,' readin mixtj from GWinput: mixtj=',acc
        tjmax=abs(acc)+1d-3
        if(mmix==1) then
          tj(1)=1d0
        else
          tj(1)= acc
          tj(2)= 1-acc
          mmix=2
        endif
        ido=2
      else
        tjmax=5d0
        ido=0
      endif
      imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax,
c     .  w(onorm),w(okpvt),a,tj,rms2)
     .  norm,kpvt,a,tj,rms2)
      deallocate(norm, kpvt)
c      call rlse(onorm)

C ... Restore PQ array, updating new x
c      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
c      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
c ...
      sss = a(:,0,2)
      rewind(ifi)
      write(ifi) nitr+1,nda
      write(ifi) a
      close(ifi)
      end

c
C      subroutine readd_iSigma_en(ifinin,iSigma_en)
C      use keyvalue
C      integer(4):: iSigma_en
Cc      logical :: readgwinput
Cc      if(readgwinput()) then
C      call getkeyvalue("GWinput","iSigMode",iSigma_en )
Cc      else
Cc        ifinin=iopen('GWIN_V2',1,0,0)
Cc        read(ifinin,*)
Cc        read(ifinin,*)
Cc        read(ifinin,*) blank,blank,iSigma_en
Ccc        iSigma_en parameter which determines approximation for  self-energy,
Ccc        stored in GWIN_V2.
Ccc        iSigma_en==0 SE_nn'(ef)+imag intgrl:delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ccc        iSigma_en==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ccc        iSigma_en==2 SE_nn'((e_n+e_n')/2)
Ccc        iSigma_en==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
Ccc        iSigma_en==5 delta_nn' SE_nn(e_n)
Ccc        output file should contain hermitean part of SE for energies to be real
Ccc        (for example, hermitian conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
Cc        is=iclose('GWIN_V2')
Cc      endif
CC      write(6,"(' Sigma mode',i2)") iSigma_en
C      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine anfsig( natom,ldima,bas, nbloch,
     i  q,           ! for q+G
     i  anfvec,iaf,  ! these are antiferro informations.
     i  zxq)  ! i/o
C- antiferro part is added to sism
c We assume that the crystal has a magnetic symmetry described by (translataion + spin flip).
c The translation is specified by a vector,
c   AFvector = anfvec(1:3)*alat, which is the true real vector in Cartesian coodinate.
c   The anfvec(1:3) is given in ANFcond file.
c
c Each mixed basis is mapped to the other mixed basis.
c E.g. the product basis B({\bf r}-{\bf a}) is mapped to
c  B({\bf r}-{\bf a}-{\bf A}) = B({\bf r}-{\bf a}'-{\bf T}_0),
c  by the translation specified by AFvec={\bf A}.
c  Here {\bf T}_0 is some crystal tralslation vector.
c  In this code you see,
c      bas(1:3,ia1)+ anfvec  = bas(1:3,iaf(ia1)) + transaf(1:3,ia1)
c  ==   {\bf a}    +{\bf A}  = {\bf a}'          + {\bf T}_0
c
c ---- The corresponding atoms should have the same product basis.
c
      implicit none
      integer natom, nbloch
      integer ldima(natom), iaf(natom), iof(natom),
     &  ia1,i,im1,im2
      integer(4):: iam(nbloch),imf(nbloch)
      real(8) :: qt(natom),q(3),bas(3,natom),
     &         anfvec(3),transaf(3,natom)
      complex(8) :: zxq (nbloch,nbloch),fac(nbloch)
      complex(8):: imag=(0d0,1d0)
      complex(8),allocatable :: zxqw(:,:)
      real(8) :: pi=3.1415926535897932d0
cr
cr    True_q(1:3)     = 2*pi/alat * q(1:3)
cr  True G is given by
cr    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvecc(1:3,igp)) ,igp=1,ngp
c------------------------
c      imagtwopi = 2d0*(0d0,3.1415926535897932d0)
      write(6,"(' q=', 3d23.15)") q
      do ia1 = 1, natom
        transaf(1:3,ia1)= bas(1:3,ia1)+ anfvec - bas(1:3,iaf(ia1))
        write(6,"(' ia1 transaf=',i3, 3d23.15)") ia1,transaf(1:3,ia1)
      enddo
      iof(1) = 0
      do ia1 = 1, natom-1
        iof(ia1+1)= iof(ia1) + ldima(ia1)
        iam(iof(ia1)+1:iof(ia1+1)) = ia1
      enddo
      iam(iof(natom)+1:nbloch) = natom
      write(6,*) (ia1, ldima(ia1),ia1 = 1, natom)


c (we had a problem at June18 2007)
c  this sanitiy check should be removed in fufture...
      if( nbloch /= iof(natom) +ldima(natom) )
     &  call rx( ' anfsig: nbloch /= ... not yet for floting orbital case')

c phase shifts
      do ia1 = 1, natom
        qt(ia1) = 2d0*pi*sum(q*transaf(1:3,ia1))
c        write( 6, "(i3,2x,d23.15)") iaf(ia1),sum(q*transaf(1:3,ia1))
c        write( 6, "(3x,2x,d23.15)") 2d0*pi*sum(q*transaf(1:3,ia1))
      enddo
      write(6,*) ' anfsig:  nbloch=',nbloch
      do im1 = 1,nbloch
        ia1  = iam(im1)
        i    = im1 - iof(ia1)
        imf (im1) = i + iof(iaf(ia1))
        fac (im1) = exp( imag*qt(ia1) )
cccccccccccccccccccccccc
        write( 6, "(5i6,2d13.5)") ia1,iaf(ia1), i, im1, imf(im1)
     &  ,fac(im1)
cccccccccccccccccccccccc
      enddo
      allocate(zxqw(nbloch,nbloch))
      zxqw = zxq
      do im1 = 1,nbloch
        do im2 = 1,nbloch
          zxq(im1, im2)=
     &               + zxqw(imf(im1),imf(im2))*fac(im1)*dconjg(fac(im2))
        enddo
      enddo
      write(6,*) ' anfsig: end '
      return
      end
