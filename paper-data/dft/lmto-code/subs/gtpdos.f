      subroutine gtpdss(mode,nfilem,s_bz,lidos,nchan,nsp,
     .  npts,emin,emax,wg,dos)
C- Partial density-of-states generator from weights file
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: nkabc nkp ntet n w
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:idtet wtkp
Cio    Passed to: *
Ci Inputs
Ci   mode  :0 if any of passed nchan,nsp,nkp are zero, use values from
Ci         :  file data.  Note: it is dangerous to pass nchan=0.
Ci         :  If you do, ensure dos(npts,nchan) is sufficiently large
Ci         :  to handle any nchan in the file.
Ci         :1 require match in nchan,nsp,nkp
Ci   nfilem:file logical unit containing dos weights
Ci   lidos :true => return integrated DOS
Ci   nchan :number of channels expected
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   npts  :number of points in DOS window
Ci   emin  :beginning of DOS energy window
Ci   emax  :end of DOS energy window
Ci   wg    :DOS is broadened by a gaussian of width wg
Co Outputs
Co   dos   :density of states for each of nchan channels
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Aug 10 works in parallel mode
Cu   16 Apr 10 Allow gaussian broadening of DOS
Cu   18 Mar 10 Adapted from asados
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lidos
      integer mode,nfilem,nchan,npts,nsp
      double precision emin,emax,wg,dos(npts,nchan*nsp)
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      real(8),allocatable:: emesh(:),fc(:)
C ... Local parameters
      logical ldum
      integer mpsord,nkabc(3),ntet,nkp,nchanl,nkpl,nspl,mode0,i,ich
      double precision swidth,def,xx
C ... MPI
      integer mpipid,procid,master

      procid = mpipid(1)
      master = 0
      mode0 = mod(mode,10)

C     Parameters for tetrahedron integration
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
C     Sampling parameters
      mpsord = s_bz%n
      swidth = s_bz%w

      if (mode0 == 0) then
C       nkabc, ntet are not used in this call
        call gtpdos(0,nfilem,ldum,nchanl,nkpl,nspl,nkabc,ntet,xx,
     .    xx,xx,xx,xx,xx,xx,xx)
        if (nkp == 0) nkp = nkpl
        if (nsp == 0) nsp = nspl
        if (nchan == 0) nchan = nchanl
      endif

      call gtpdos(mode-mode0+1,nfilem,lidos,nchan,nkp,nsp,nkabc,ntet,
     .  s_bz%idtet,mpsord,swidth,s_bz%wtkp,npts,emin,emax,dos)

      if (wg /= 0) then
        allocate(emesh(npts),fc(npts))
        def = (emax-emin)/(npts-1)
        do  i = 1, npts
          emesh(i) = emin + (i-1)*def
        enddo
        do  ich = 1, nchan*nsp
          call dcopy(npts,dos(1,ich),1,fc,1)
C         call prrmsh('before convolve',emesh,dos(1,ich),npts,npts,1)
          call convolve(2,0d0,0d0,0d0,0d0,wg,npts,emesh,fc,dos(1,ich))
C         call prrmsh('after convolve',emesh,dos(1,ich),npts,npts,1)
        enddo
        deallocate(emesh,fc)
      endif

      end
      subroutine gtpdos(mode,nfilem,lidos,nchan,nkp,nsp,
     .  nkabc,ntet,idtet,mpsord,swidth,wtkp,npts,emin,emax,dos)
C- Partial density-of-states generator from weights file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 read nchan,nsp,nkp and exit
Ci         :1 generate DOS
Ci   nfilem:file logical unit containing dos weights
Ci   lidos :true => return integrated DOS
Ci   nkabc :no. divisions for the 3 recip. latt. vectors
Ci   ntet  :number of tetrahedra.  ntet=0 => sampling integration
Ci   idtet :(from TETIRR) information for tetrahedron integration
Ci   mpsord:Polynomial order, Methfessel-Paxton sampling integration
Ci   swidth:gaussion broadening, Methfessel-Paxton sampling integration
Ci   wtkp  :k-point degeneracies for M-P sampling integration
Ci   npts  :number of points in DOS window
Ci   emin  :beginning of DOS energy window
Ci   emax  :end of DOS energy window
Cio Inputs/Outputs
Cio  If mode=0 the following are read from disk, otherwise they are input
Cio  nchan :number of channels
Cio        :match to file required if mode>0
Cio  nkp   :number of irreducible k-points (bzmesh.f)
Cio        :match to file required if mode>0
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio        :match to file required if mode>0
Co Outputs
Co    dos  :densities-of-states for each of nchan channels
Cl Local variables
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   15 Aug 10 works in parallel mode
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lidos
      integer mode,nfilem,nchan,nkp,nsp,nkabc(3),ntet,idtet(0:4,*),
     .  npts,mpsord
      double precision emin,emax,swidth,dos(npts,nchan*nsp),wtkp(nkp)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk(:)
      real(8), allocatable :: doswt(:)
C ... Local variables
      real(8),allocatable:: eband(:,:)
      double precision efermi,vmtz,xx(1)
      integer nlf,nspl,nspc,nfstg,nkpl,nband,ndum,j,nevmx,nbandx,nspx
C ... MPI
      integer mpipid,procid,master

      procid = mpipid(1)
      master = 0

C ... Get dimensions (nspl,nspc,nkpl,nband,nfstg); allocate memory for eband
      if (procid == master) then
        call iomomq(nfilem,0,nlf,nspl,nspc,nkpl,nband,nfstg,ndum,ndum,
     .    ndum,ndum,ndum,xx,xx,xx,xx,xx,xx)
      endif
      call mpibc1(nspl,1,2,.false.,'gtpdos','nspl')
      call mpibc1(nspc,1,2,.false.,'gtpdos','nspc')
      call mpibc1(nkpl,1,2,.false.,'gtpdos','nkpl')
      call mpibc1(nband,1,2,.false.,'gtpdos','nband')
C     call mpibc1(nfstg,1,2,.false.,'gtpdos','nfstg')
      nspx = nsp / nspc
      if (mode == 0) then
        if (procid == master) then
          read (nfilem) nchan
        endif
        call mpibc1(nchan,1,2,.false.,'gtpdos','nchan')
        nkp = nkpl
        nsp = nspl
        return
      endif
      nbandx = nband*nspc
      call sanrg(.true.,nkpl,nkp,nkp,'GTPDOS','file nkp')
      call sanrg(.true.,nspl,nsp,nsp,'GTPDOS','file nsp')

C ... not parallel for now
      if (procid == master) then
      allocate(eband(nband,nkp*nsp))
      j = nchan*nbandx*nsp*(nkp+1)
      allocate(doswt(j)); call dpzero(doswt,j)

C ... Read eband,doswt
      if (procid == master) then
        call iomomq(nfilem,32,nlf,nsp,nspc,nkp,nband,nfstg,j,nband,
     .    nchan,nchan,nevmx,eband,doswt,doswt,xx,efermi,vmtz)
      endif
C     call mpibc1(j,1,2,.false.,'gtpdos','nchan')
      if (j /= nkp) call rx('gtpdos: moments file missing qpts')

C --- make DOS or NOS ---
      allocate(wk(npts))
      if (ntet > 0) then
        call dostet(nbandx,nsp,nspx,nevmx,nchan,nkabc(1),
     .    nkabc(2),nkabc(3),ntet,idtet,eband,
     .    doswt,npts,emin,emax,lidos,wk,dos)

C ... dosspl reads bands and moms on the fly to save work space
      else
        call dosspl(nfilem,nbandx,nsp,nspc,nchan,mpsord,swidth,nkp,
     .    wtkp,eband,doswt,npts,emin,emax,lidos,wk,dos)
      endif
      deallocate(wk,doswt,eband)
      endif

      call mpibc1(dos,npts*nchan*nsp,4,.false.,'gtpdos','dos')

      end
