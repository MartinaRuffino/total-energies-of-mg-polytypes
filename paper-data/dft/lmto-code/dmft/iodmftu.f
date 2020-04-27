      subroutine iodmftu(s_dmft,rot,ndham,nlo,nhi,ldcix,nsp,ncix,iq,ig,qp,evlk,dmftu,ifi)
C- Reads dmft projectors U
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Cio  s_dmft
Ci     Elts read:  l
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:icix
Cio    Passed to:  *
Ci Inputs
Ci   ndham :dimensioning parameter, largest hamiltonian dimension
Ci   nlo   :dmftU is read for bands ilo:ihi
Ci   nhi   :dmftU is read for bands ilo:ihi
Ci   ldcix :Second dimension of dmftU
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ncix  :Number of correlated blocks
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   rot   : if T, the cubic harmonics of the correlated block cix are rotated in order to correspond to the cubic harmonic of icix(cix)
Co Outputs
Ci   iq    :iq for which projectors were made
Ci   ig    :index to which entry in star of k.  ig=1 => first entry in this star
Ci   qp    :k-point
Ci   evlk  :QP levels for this iq
Ci   dmftu :Projectors for this iq
Cs Command-line switches
Cl Local variables
Cr Remarks
Cu Updates
Cu   27 Nov 17
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, intent(in) :: ndham,nlo,nhi,ldcix,nsp,ncix,ifi
      integer, intent(out) :: iq,ig
      real(8), intent(out) :: qp(3),evlk(ndham,nsp)
      complex(8) :: dmftu(nlo:nhi,ldcix,nsp,ncix)

C ... For structures
!      include 'structures.h'
      type(str_dmft)::  s_dmft
C ... Dynamically allocated local arrays
      complex(8), allocatable :: dmftUi(:,:,:)
C ... Local parameters
      integer cix,i,icix,nl1,j,jsp,nevn
      integer procid,master
      double precision qpr(3)
!     for rotation
      logical, intent(in) :: rot
      complex(8),allocatable :: rmattmp(:,:),tmp(:,:)


      procedure(integer) :: mpipid

      procid = mpipid(1)
      master = 0

      if (ifi < 0) call rx('iodmftu : not ready for writing')

      call dpzero(dmftu,2*size(dmftu))
      do  cix = 1, ncix
        icix = iabs(s_dmft%icix(cix))
        nl1 = 2*s_dmft%l(icix)+1
        if (procid == master) then
          read(ifi) i,nevn,j,qpr,iq,ig,qp ! Read from proj file
        endif
C       print *, 'iodmftu read iq,ig',iq,ig
        call mpibc1(i,1,2,.false.,'','')
        call mpibc1(iq,1,2,.false.,'','')
        call mpibc1(nevn,1,2,.false.,'','')
        call mpibc1(j,1,2,.false.,'','')
        call mpibc1(ig,1,2,.false.,'','')
        if (i /= cix .or. j /= nl1*nsp) call rx('sudmft: file mismatch')
        allocate(dmftUi(nevn,nl1,nsp)) ! Allocate here to ensure contiguous spin 1,2 blocks
        if (procid == master) then
           call dpdump(evlk,ndham*nsp,ifi)
           call dpdump(dmftUi,2*nevn*nl1*nsp,ifi)
           if (rot) then
              do  jsp = 1, nsp
                 if( s_dmft%ig(cix) /= 1 ) then
                    if (.not. allocated(rmattmp)) then
                       allocate(rmattmp(ldcix,ldcix))
                       allocate(tmp(nevn,ldcix))
                    endif
                    forall(i=1:ldcix,j=1:ldcix) rmattmp(i,j) = (s_dmft % rmat(i,j,cix))
                    call zgemm('N','N',nevn,ldcix,ldcix,(1d0,0d0),
     .                   dmftUi(1,1,jsp),nevn,rmattmp,ldcix,(0d0,0d0),tmp,nevn)
                    dmftUi(:,:,jsp) = tmp(:,:)
                 endif
              enddo
           endif
        endif
        call mpibc1(evlk,ndham*nsp,4,.false.,'','')
        call mpibc1(dmftUi,2*nevn*nl1*nsp,4,.false.,'','')
        do  jsp = 1, nsp
           dmftu(nlo:nlo-1+nevn,1:nl1,jsp,cix) = dmftUi(1:nevn,1:nl1,jsp)
        enddo
        deallocate(dmftUi)
      enddo

      end
