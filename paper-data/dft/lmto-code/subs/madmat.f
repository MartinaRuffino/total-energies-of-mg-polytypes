      subroutine madmat(nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg,dmad)
C- Coefficients to Madelung matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   awald :Ewald parameter, scales with the lattice as: as/(vol)**(1/3)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci         :Sign of volume is used as a flag: if <0 => serial mode
Ci   dlat  :direct lattice vectors
Ci   nkd   :number of dlat for Ewald sums
Ci   glat  :reciprocal lattice vectors
Ci   nkg   :number of reciprocal lattice vectors for Ewald sums
Co Outputs
Co   dmad  :Madelung matrix
Cr Remarks
Cr   The potential from a basis of more than one atom is obtained by
Cr   superposing the potential from each sublattice.  Matrix element
Cr   dmad(i,j) is the (1/2)potential at position tau(i) from unit
Cr   charges on sublattice j, compensated by a uniform background of
Cr   density 1/vol.  (if tau(i)=0, the potential from the singularity
Cr   0 is not included.)
Cr
Cr   Call lattc to generate inputs awald,vol,dlat,nkd,glat,nkg.
Cu Updates
Cu   17 Oct 17 madmat can operate in serial or MPI parallel mode
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nkd,nkg
      double precision awald,alat,vol
      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
      double precision dmad(nbas,nbas)
C ... Local parameters
      integer i,j,nproc,procid,master,inest,nnest,nproci
      double precision tau(3),truevol
      integer, allocatable :: kpproc(:),index(:,:)
      procedure(integer) :: iprint,nglob,mpipid,fopna

      call tcn('madmat')

      truevol = abs(vol)
      nproc = mpipid(0); procid = mpipid(1); master = 0

C --- Generate Madelung matrix ---
C ... Serial mode
      if (vol < 0) then
      do  i = 1, nbas
        do  j = i, nbas
          tau(1:3) = bas(1:3,j)-bas(1:3,i)
          call shortn(tau,tau,dlat,nkd)
          call strx00(tau,awald,alat,truevol,glat,nkg,dlat,nkd,dmad(i,j))
C         call strxq(0,0d0,[0d0,0d0,0d0],tau,1,1,1,alat,truevol,awald,nkd,nkg,dlat,glat,cg,indxcg,jcg,s,sd)
          dmad(j,i) = dmad(i,j)
        enddo
      enddo
      goto 99
      endif

C ... MPI mode
      allocate (kpproc(0:max(nproc,1))); kpproc = 0 ! kpproc=0 for any procid not used
      nnest = (nbas*(nbas+1))/2  ! Make big list of (i,j) elements.  Use dmad(i,j)=dmad(j,i)
      nproci = min(nproc,nnest)  ! Limit the number of processors to the number of elements
      allocate(index(nnest,2))   ! table of (i,j) indices for each memmber of list
      if (nproci > 1) then
        call info2(30,1,0,' MADMAT: Generate Madelung matrix distributed over %i processors ...',
     .    nproci,2)
      endif
      call dstrbpx(nnest,nproci,[-1d0],kpproc)

C     Create a vector of (i,j) pairs
      inest = 0
      do  i = 1, nbas
        do  j = i, nbas
          inest = inest+1
          index(inest,1) = i; index(inest,2) = j
        enddo
      enddo
      if (inest /= nnest) call rx('bug in setting up nested loop')

      call dpzero(dmad,size(dmad))  ! Results from all processes will be summed together
      if (kpproc(procid) > 0) then
      do  inest = kpproc(procid), kpproc(procid+1)-1  ! Loop over all (i,j) pairs
        i = index(inest,1); j = index(inest,2)
!       print 333, procid, inest, i, j
!  333 format(6i6)
        tau(1:3) = bas(1:3,j)-bas(1:3,i)
        call shortn(tau,tau,dlat,nkd)
        call strx00(tau,awald,alat,truevol,glat,nkg,dlat,nkd,dmad(i,j))
        dmad(j,i) = dmad(i,j)
      enddo
      endif
      deallocate(index)
      call mpibc2(dmad,nbas*nbas,4,3,.false.,'madmat','dmad')

C ... Re-entry point for serial mode
   99 continue
      call tcx('madmat')

C --- Write to file madmat if verbosity is high enough ---
      if (iprint() < 90) return
      call info0(10,1,0,' ... writing file madmat')
      if (procid == master) then
        i = fopna('madmat',-1,0)
        call ywrm(0,' ',1,i,'(10f14.8)',dmad,0,nbas,nbas,nbas)
      endif
!     call rx0('done')
      end

      subroutine smadmat(nbas,bas,e,awald,alat,vol,dlat,nkd,glat,nkg,
     .  cg,indxcg,jcg,dmad)
C- Coefficients to screened Madelung matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   awald :Ewald parameter, scales with the lattice as: as/(vol)**(1/3)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci         :Sign of volume is used as a flag: if <0 => serial mode
Ci   dlat  :direct lattice vectors
Ci   nkd   :number of dlat for Ewald sums
Ci   glat  :reciprocal lattice vectors
Ci   nkg   :number of reciprocal lattice vectors for Ewald sums
Co Outputs
Co   dmad  :Madelung matrix
Cr Remarks
Cr   The potential from a basis of more than one atom is obtained by
Cr   superposing the potential from each sublattice.  Matrix element
Cr   dmad(i,j) is the (1/2)potential at position tau(i) from unit
Cr   charges on sublattice j, compensated by a uniform background of
Cr   density 1/vol.  (if tau(i)=0, the potential from the singularity
Cr   0 is not included.)
Cr
Cr   Call lattc to generate inputs awald,vol,dlat,nkd,glat,nkg.
Cu Updates
Cu   17 Oct 17 madmat can operate in serial or MPI parallel mode
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nkd,nkg,indxcg(*),jcg(*)
      double precision e,awald,alat,vol,cg(*)
      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
      double precision dmad(nbas,nbas)
C ... Local parameters
      integer i,j,nproc,procid,master
      double precision tau(3),truevol
      integer, allocatable :: kpproc(:),index(:,:)
      procedure(integer) :: iprint,nglob,mpipid,fopna
      complex(8) :: zx(1)
!     real(8) :: xx = 0d0
      call tcn('madmat')

      truevol = abs(vol)
      nproc = mpipid(0); procid = mpipid(1); master = 0


C --- Generate Madelung matrix ---
C ... Serial mode ... no MPI yet.
      if (vol < 0 .or. .true.) then
      do  i = 1, nbas
        do  j = i, nbas
          tau(1:3) = bas(1:3,j)-bas(1:3,i)
          call shortn(tau,tau,dlat,nkd)
!         call strx00(tau,awald,alat,truevol,glat,nkg,dlat,nkd,dmad(i,j))
          call strxq(100,e,[0d0,0d0,0d0],tau,1,1,1,alat,truevol,awald,
     .      nkd,nkg,dlat,glat,cg,indxcg,jcg,zx,zx)
C         print *, dmad(i,j),dble(zx(1))
          dmad(i,j) = zx(1)
          dmad(j,i) = zx(1)
        enddo
      enddo
      goto 99
      endif

C ... Re-entry point for serial mode
   99 continue
      call tcx('madmat')

C --- Write to file madmat if verbosity is high enough ---
      if (iprint() < 90) return
      call info0(10,1,0,' ... writing file smadmat')
      if (procid == master) then
        i = fopna('smadmat',-1,0)
        call ywrm(0,' ',1,i,'(10f14.8)',dmad,0,nbas,nbas,nbas)
      endif
C     call rx0('done')
      end
