      module m_ReadEfermi
      real(8),protected:: bandgap
      real(8)::ef
      contains

      subroutine readefermi()
      implicit none
      integer:: ifief,ifile_handle
      ifief=ifile_handle()
      open(ifief,file='EFERMI')
      read(ifief,*) ef,bandgap
      close(ifief)
      write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
      end subroutine
      end module m_ReadEfermi

      module m_genallcf_v3
C- get basic settings of crystal structure and nlm info
c----------------------------------------------------------------
c> structure
c     o                   plat,alat,natom,nclass,pos,
c     o                   ngrp, symgg,
c     o                   invg, ef,
c>> l,n and dimensions
c     o                   clabl, nspin,nl,nn,nnv,nnc,
c     o                   nindx, nindxv, nindxc, iclass,
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c
c     o                   z,
c>> l,n,m indices for Phi
c     o                   il, in, im,  ilnm,  nlnm,
c     o                   ilv,inv,imv, ilnmv, nlnmv,
c     o                   ilc,inc,imc, ilnmc, nlnmc,
c>> core
c     o                   ncwf, ecore, konf, icore, ncore,nctot,
c> frequency
c     o                   niw,diw,nw,dw,delta,deltaw,esmr, freq)
c             symgrp
c             ,nocc, nunocc, occv, unoccv, occc, unoccc
c----------------------------------------------------------------
      implicit none

      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot, niw,nw
      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
      character(120):: symgrp
      character(6),allocatable :: clabl(:)
      integer(4),allocatable:: iclass(:)
     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
     &    occv(:,:,:),unoccv(:,:,:)
     &   ,occc(:,:,:),unoccc(:,:,:),
     o    nocc(:,:,:),nunocc(:,:,:)
      real(8), allocatable::
     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:)   !w(igrp)
c-----------------------------------------------------------------------

      contains

      subroutine genallcf_v3(nwin,efin,incwfx)
C- Readin GWIN_V2 and LMTO(crystal) data and allocate all required.
C
Cr Return iclass=ibas.
Cr nwin,efin,incwfx, are used as switches.
Ci input: nwin,efin,incwfx,
Ci        GWIN_V2, LMTO
Co output: All the output are given in the declarations above.
Cr Remarks
Cr
C    For MPI, genallcf_v3 needs to be called on all ranks, and
C    will handle rank=0 reading and broadcasting internally.
C--------------------------------------------------------

      use mpi_mod, rank => mpi_rank
      implicit none
      integer(4)::iflmto,ifinin,nwin,incwfx,ifec,i,j,
     & lmx, lmx2,nlmto2,nprodxc,nlnaxc,nlnaxv,nprodx,ifi,ig,is
     & ,iopen,iclose,nprodxv,nlnax
     & ,noflmto,maxnn
      integer(4):: infwfx
      integer(4):: n1,n2,n3,imagw,lcutmx,n,ic
      logical :: nocore
      real(8)::efin
      real(8),allocatable::tolbas(:)
      character(120):: symgrpt
      real(8),   allocatable:: ecoret(:)
      integer(4),allocatable::ncwf2(:,:,:)

      print "(/' genallcf: read data from GWinput, LMTO, ECORE, EFERMI, NW')"

c      logical ::readgwinput
c----------------------------------------------
      if(rank == 0) then
        iflmto = iopen('LMTO',1,0,0)
        if (iflmto < 0) call rx( 'unit file for < 0')

        ifi = iflmto
        nw  = nwin
        ef  = efin
        read(ifi,*);  read(ifi,*)
        read(ifi,*)symgrpt        !SYMMETRY
        symgrp='  '//trim(adjustl(symgrpt))
        print *,'symgrp = ', trim(symgrp)

        read(ifi,*)
        read(ifi,*)
        read(ifi,*)
        read(ifi,*)alat              !lattice constant
        allocate(plat(3,3))          !primitive lattice vectors
        read(ifi,*)
        read(ifi,*)plat(1:3,1)
        read(ifi,*)plat(1:3,2)
        read(ifi,*)plat(1:3,3)
        read(ifi,*)
        read(ifi,*) natom            !Number of atoms
      endif
#ifdef USE_MPI
      call MPI_Bcast(symgrp, 120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(alat, 1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(plat, 9, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(alat, 1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(natom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif
      nclass = natom            !This code does not distinguish class yet.  Each site gets its own class
      allocate(pos(3,natom))    !positions of atoms

      if(rank == 0) then
        read(ifi,*)
        do n = 1,natom
          read(ifi,*) pos(1,n),pos(2,n),pos(3,n)
        end do
        read(ifi,*)
        read(ifi,*)
        read(ifi,*)
        read(ifi,*)nspin           !spin (1=paramagnetic  2=ferromagnetic)
        read(ifi,*)
        read(ifi,*)nl              !max. no. valence and core l
        read(ifi,*)
        read(ifi,*)nnv,nnc   !max. no. valence and core n
        write(*,321) nspin,nl,nnv,nnc
  321   format(' nspin=',i1,'  nl=',i1,' number of valence channels per l=',i1,' number of core states =',i3)

        if(nnv==1) nnv=2        ! for backward compatibility!takao apr 2002
                                ! nnv=2 corresponds to phi and phidot
                                ! nnv=3 corresponds to

        read(ifi,*)
        read(ifi,*)             !nrx is not readin
        read(ifi,*)
      endif
#ifdef USE_MPI
      call MPI_Bcast(pos,3*natom,MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nspin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nl,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nnv,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nnc,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      allocate(clabl(nclass),z(nclass),iclass(natom)) !class-label, atomic number, class index
      allocate(nindxv(nl,nclass), nindxc(nl,nclass),
     &        occv(nl,nnv,nclass),unoccv(nl,nnv,nclass),
     &        occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))

      allocate(tolbas(0:2*(nl-1)))

      if(rank == 0) then
        do ic = 1,nclass
          read(ifi,*) clabl(ic),z(ic) !,nrofi is not readin
        end do
        do n = 1,natom          !!We set nclass = natom through the GW calculations
          iclass(n)=n
        end do
      endif
#ifdef USE_MPI
      call MPI_Bcast(clabl, nclass*6, MPI_CHARACTER,  0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(z, nclass, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(iclass, natom, MPI_INTEGER,      0, MPI_COMM_WORLD, ierror)
#endif

      if(rank == 0) then
        call rgwinaf(ifi,ifinin,nl,nnv,nnc,nclass,
     o               n1,n2,n3,ef,  !BZ
     o               niw,diw,nw,dw,delta,deltaw,esmr,imagw,  ! frequencies
     o               tolbas,lcutmx,nindxv,nindxc,occv,unoccv,occc,unoccc, ! product basis
     o               ncwf,ncwf2) ! core
      endif
#ifdef USE_MPI
!      call MPI_Bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
!      call MPI_Bcast(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
!      call MPI_Bcast(n3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ef,  1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(niw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(diw, 1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nw,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(dw,  1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(delta,  1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(deltaw, 1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(esmr,   1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(imagw,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(tolbas, 2*(nl-1)+1,  MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(lcutmx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nindxv, nl*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nindxc, nl*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(occv,   nl*nnv*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(unoccv, nl*nnv*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(occc,   nl*nnc*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(unoccc, nl*nnc*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ncwf,   nl*nnc*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ncwf2,  nl*nnc*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      if(rank == 0) inquire(file='NoCore',exist=nocore)
#ifdef USE_MPI
      call MPI_Bcast(nocore,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
#endif

      if(nocore) then
        occc=0
        unoccc=0
        ncwf  =0
      elseif( incwfx==-1 ) then
        print *,' ### incwf=-1 Use ForSxc for core'
        ncwf = ncwf2  !call icopy(nl*nnc*nclass,w(incwf2),w(incwf))
      elseif( incwfx==-2 ) then
        print *,' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
        call notbit(nl*nnc*nclass, ncwf2)
        ncwf  = ncwf2
        occc= ncwf
        unoccc= 0
      elseif( incwfx==-3 ) then
        call ibiton(nclass,nl,nnc,nindxc, occc, ncwf)
        unoccc= 0
        print *,' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
      elseif( incwfx==-4 ) then
        print *,' ### incwf=-4  occ=0 and unocc=0 for all core '
        occc=0
        unoccc=0
        ncwf=0
      elseif(incwfx==0) then
        print *,' ### Use unocc occ ForX0 for core'
      else
        call rx( ' ### proper incwf is not given for genallcf2:rgwinf ')
      endif
      deallocate(ncwf2)
C... End of rgwinf section ----------------------------


c> dimensions and constants
      lmx        = 2*(nl-1)
      lmx2       = (lmx+1)**2
      nlmto      = noflmto(nindxv,iclass,nl,nclass,natom)
      nlmto2     = nlmto*nlmto
      nn         = maxnn (nindxv,nindxc,nl,nclass)

c>> combine nocc,nunocc,nindx
      allocate(nindx(nl,nclass))
      allocate(nocc(nl,nn,nclass),nunocc(nl,nn,nclass))
      call reindx  (occv,unoccv,nindxv,
     i              occc,unoccc,nindxc,
     d              nl,nn,nnv,nnc,nclass,
     o              nocc,nunocc,nindx)
      call maxdim  (occc,unoccc,nindxc,nl,nnc,nclass,
     o              nprodxc,nlnxc,nlnmxc,nlnaxc)
      call maxdim  (occv,unoccv,nindxv,nl,nnv,nclass,
     o              nprodxv,nlnxv,nlnmxv,nlnaxv)
      call maxdim  (nocc,nunocc,nindx,nl,nn,nclass,
     o              nprodx,nlnx,nlnmx,nlnax)

      if(rank == 0) print *,' nw from rgwinaf=',nw

      if(nw>0) then
        allocate(freq(nw))
        call genfreq (nw,dw,0.d0,freq)
      endif

      allocate(icore(nl**2*nnc,nclass),ncore(nclass))
      icore=9999999
      ncore=9999999
      call incor   (ncwf,nindxc,iclass,
     d              nl,nnc,nclass,natom,
     o              icore,ncore,nctot )

      allocate(konf(nl,nclass),ecore(nctot,2))
      konf=0

      if(rank == 0) then
        ifec = iopen ('ECORE',1,0,0)
        allocate(ecoret(nl*nnc*2*nclass))
        ecoret = 0d0
        call rcore2(ifec,nl,nnc,nclass, konf,ecoret)
        call gencor  (ecoret,ncwf,iclass,nl,nnc,natom,nclass,nctot, ecore)
        deallocate(ecoret)
      endif
#ifdef USE_MPI
      call MPI_Bcast(konf,  nl*nclass, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ecore, 2*nctot,   MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
#endif

c     indices for core and LMTO basis
      allocate(
     & il(nlnmx,nclass),
     & in(nlnmx,nclass),
     & im(nlnmx,nclass),
     & ilnm(nn*nl*nl*nclass),
     & ilv(nlnmxv*nclass),
     & inv(nlnmxv*nclass),
     & imv(nlnmxv*nclass),
     & ilnmv(nnv*nl*nl*nclass),
     & ilc(nlnmxc*nclass),
     & inc(nlnmxc*nclass),
     & imc(nlnmxc*nclass),
     & ilnmc(nnc*nl*nl*nclass)
     & )

      call idxlnmc ( nindxv,nindxc,
     d              nl,nn,nnv,nnc,nlnmx,nlnmxv,nlnmxc,nclass,
     o              il,in,im,ilnm,
     o              ilv,inv,imv,ilnmv,
     o              ilc,inc,imc,ilnmc)

      allocate(nlnmv(nclass),nlnmc(nclass),nlnm(nclass))
      call nolnma  (nindxv,nl,nclass,
     o              nlnmv )
      call nolnma  (nindxc,nl,nclass,
     o              nlnmc )
      call nolnma  (nindx,nl,nclass,
     o              nlnm )
      i=2 !see previous definition of symgrp
      if(symgrp(i+1:i+13)/= 'UseSYMOPSfile') then
        call rx( " Not: UseSYMOPSfile in LMTO file")
      endif
      if(rank == 0) then
        write(6,*) ' symgrp==UseSYMOPSfile'
        ifi = 6661
        open (ifi, file='SYMOPS')
        read(ifi,*) ngrp
      endif
#ifdef USE_MPI
      call MPI_Bcast(ngrp,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      allocate(symgg(3,3,ngrp))

      if(rank == 0) then
        do ig = 1,ngrp
          read(ifi,*)  ! translation is not read; calculated on the fly in mptauof
          do i=1,3
            read(ifi,"(3d24.16)") symgg(i,1:3,ig)
          enddo
        enddo
        close(ifi)
      endif
#ifdef USE_MPI
      call MPI_Bcast(symgg, 9*ngrp, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
#endif

c> index for inverse rotations
c      call defi    (iinvg,ngrp)
      allocate(invg(ngrp))
      call invgrp  (symgg,ngrp,
     o              invg)

      if(rank == 0) then
        is = iclose ('LMTO')
        is = iclose ('ECORE')
      endif

      call cputid(0)
      if(rank == 0) write (*,*) 'genallcf_v3: end '

      end subroutine
      end module


c-------------------------------------------------------------------
c taken from bas.f
      subroutine gencor (eco,ncwf,iclass,
     d                   nl,nnc,natom,nclass,nctot,
     o                   eclnm)

c 92.03.25
c generate core energies for each l,n,m

c eco(l,n,cl)  = core energies
c ncwf(l,n,cl) = 1 ==> allowed core state
c                0 ==> not allowed
c nnc          = max. no. core states/l
c nctot        = total no. core states (all atoms)

c eclnm({l,n,m,atom}) = core energies

      implicit real*8 (a-h,o-z)
      dimension eco(0:nl-1,nnc,2,nclass),ncwf(0:nl-1,nnc,nclass),
     i          iclass(natom)
      dimension eclnm(nctot,2)

      i          = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        do       l = 0,nl-1
          do       n = 1,nnc
            do       m = -l,l
              if (ncwf(l,n,ic) == 1) then
                i          = i + 1
                if (i > nctot) call rx( 'gencor: wrong nctot')
cccccccccccccccc
c      write(6,"(' eco=',4i4,d13.5)")i,l,n,ic,eco(l,n,1,ic)
cccccccccccccccc
                eclnm(i,1) = eco(l,n,1,ic)
                eclnm(i,2) = eco(l,n,2,ic)
              endif
            end do
          end do
        end do
      end do
      end

c--------------------------------------------------------------------
      subroutine rcore2 (ifec, ! minor changed version to rcore in rw.f
     d                  nl,nnc,nclass,
     o                  konfig,eco)

c 92.03.17
c read core eigenvalues both up and down

      implicit real*8 (a-h,o-z)
      character*32 char32
      dimension konfig(0:nl-1,nclass),eco(0:nl-1,nnc,2,nclass)
c loop over class
      do      ic = 1,nclass
        write(6,*) ' rcore2 : ic=',ic
        read (ifec,6000) blank
        read (ifec,6000) blank
        read (ifec,6000) blank
        read (ifec,*) z,ic1,nr,a,b,nsp
        read (ifec,6000) blank
        read (ifec,*) (konfig(l,ic),l=0,nl-1)
        read (ifec,6000) blank
c
        write(6,*) ' rcore : goto loop over l n'

c loop over l,n
        do       l = 0,nl-1
          ncore      = konfig(l,ic)-l-1
          if (ncore > nnc) call rx( 'rcore: wrong nnc')
          do       n = 1,ncore
            read (ifec,*) lt,nt,(eco(l,n,isp,ic),isp=1,nsp) !takao
            if(nsp==1) eco(l,n,2,ic) = eco(l,n,1,ic)        !
cccccccccccccccc
c      write(6,"(' eco=',3i4,d13.5)")l,n,ic,eco(l,n,isp,ic)
cccccccccccccccc
            if (lt /= l) call rx( 'rcore: wrong l')
            if (nt /= n) call rx( 'rcore: wrong n')
          end do
        end do

      end do
 6000 format(a)
 6300 format(1x,2i3,2d13.6)
 6500 format(1x,10i3)
      return
      end
