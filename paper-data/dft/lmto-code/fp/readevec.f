      subroutine readevec(mode,ndham,nsp,nspc,nkabc,nkp,lshft,nproc,evlq,zq)
C- Read eigenvalues and eigenvectors for all q from evec file
C ----------------------------------------------------------------------
Cio Structures
Ci Inputs
Co Outputs
Ci   sigdc:
Ci   mode  :1s digit (file READ) specifies how to handle non existent file
Ci         :0 Always return, whether file read is successful or not
Ci         :1 missing file: abort
Ci         :2 abort if passed ndham,nkabc,nkp,lshft do not match file
Ci         :3 combination of 1+2
Ci         :4 return dimensioning parameters ndham,nkabc,nkp,lshft
Ci         :5 combination of 1+4
Cs Command-line switches
Cl Local variables
Cr Remarks
Cr   Requires that evecs reside on disk in file 'evec'
Cr   Handles MPI case where evecs are distributed amone processor-dependent files evec.ext_procid
Cu Updates
Cu   28 Sep 17  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters:
      integer mode,ndham,nsp,nspc,lshft(3),nkabc(3),nkp,nproc
      real(8)    :: evlq(ndham,nsp,nkp)   ! Eigenvalues
      complex(8) :: zq(ndham*nspc,ndham*nspc,nsp,nkp) ! Eigenvectors

C ... Dynamically allocated local arrays
      integer, allocatable :: kpproc(:) ! For MPI distribution of processes
      complex(8), pointer :: z(:,:)     ! Eigenvectors

C ... Local variables
      integer, parameter :: LW5=5
      integer procid,master             ! For MPI
      integer ifiz,iq,isp,ndimh,ndimhx,ndhamx,mode0
      integer i1,i2,i3,i6,i7,nk123(3),lshf(3)
      double precision qp(3),eseavr(2) ! Constant term for high-lying QSGW sigma
      procedure(integer) :: fopna,mpipid

C ... Setup
      call tcn('readevec')
      procid = mpipid(1); master = 0
      mode0 = mod(mode,10)

C     Open evecs file. Header is read from master only
      ifiz = fopna('evec',-1,4); rewind ifiz
      if (procid == master) then
        call iosigh(2,LW5,i1,i3,i2,i2,nk123(1),nk123(2),nk123(3),i6,i7,lshf(1),lshf(2),lshf(3),ifiz,eseavr)
      endif
      call mpibc1(i1,1,2,.false.,'','')
      call mpibc1(i3,1,2,.false.,'','')
      call mpibc1(i2,1,2,.false.,'','')
      call mpibc1(nk123,3,2,.false.,'','')
      call mpibc1(i6,1,2,.false.,'','')
      call mpibc1(lshf,3,2,.false.,'','')
      if (mod(mode0/2,2) == 1) then
        call sanrg(.true.,i1,nsp,nsp,'readevec:','file''s nsp')
! until iosigh reads nspc        call sanrg(.true.,i3,nspc,nspc,'readevec:','file''s nspc')
        call sanrg(.true.,i2,ndham,ndham,'readevec:','file''s ndham')
        call sanrg(.true.,nk123(1),nkabc(1),nkabc(1),'readevec:','file''s nk1')
        call sanrg(.true.,nk123(2),nkabc(2),nkabc(2),'readevec:','file''s nk2')
        call sanrg(.true.,nk123(3),nkabc(3),nkabc(3),'readevec:','file''s nk3')
        call sanrg(.true.,i6,nkp,nkp,'readevec:','file''s nkp')
        call sanrg(.true.,lshf(1),lshft(1),lshft(1),'readevec:','file''s lshft1')
        call sanrg(.true.,lshf(2),lshft(2),lshft(2),'readevec:','file''s lshft2')
        call sanrg(.true.,lshf(3),lshft(3),lshft(3),'readevec:','file''s lshft3')
      elseif (mod(mode0/4,2) == 1) then
        nsp = i1
        nspc = i3
        ndham = i2
        nk123 = nkabc
        nkp = i6
        lshft = lshf
        return
      endif

C     Distribute k over processors
      allocate (kpproc(0:nproc))
      if (nproc > 1) then
        call info0(30,1,0, ' ... MPI read evals and evecs from file evec')
        if (nproc > nkp) call rxi('MPIK job cannot allocate more processors than nkp =',nkp)
        call dstrbp(nkp,nproc,1,kpproc(0))
      else
        kpproc(0:1) = [1,nkp+1]
      end if

C      ifiz = fopna('evec',-1,4); rewind ifiz
C      if (procid == master) then
C        call iosigh(3,LW5,nsp,nspc,ndham,nlmto,nk1,nk2,nk3,nkp,i7,lshft(1),lshft(2),lshft(3),ifiz,eseavr)
C      endif

C --- Read evl,z from processor-dependent file ---
      call dpzero(evlq,size(evlq))
      call dpzero(zq,2*size(zq))
      ndhamx = ndham*nspc
      do  iq = kpproc(procid), kpproc(procid+1)-1
        do  isp = 1, nsp
          read(ifiz) qp, ndimh
          ndimhx = ndimh * nspc
          allocate(z(ndimhx,ndimhx))
          call dpdump(evlq(1,isp,iq),ndhamx,ifiz)
          call dpdump(z,ndimhx**2*2,ifiz)
          zq(1:ndimhx,1:ndimhx,isp,iq) = z(:,:)
          deallocate(z)
          if (nspc == 2) exit
        enddo                   ! Loop over spins
      enddo                     ! Loop over irreducible qp
      deallocate (kpproc)

C --- Distribute each file input to all processors ---
      call mpibc2(evlq,size(evlq),4,3,.false.,'','')
      call mpibc2(zq,2*size(zq),4,3,.false.,'','')

      call tcx('readevec')
      end
