      subroutine packdmat(ldmsit,nbas,nl,dmat,s_site)
C- Poke the density matrices for non-DLM sites into s_site
C- and replace dmat for all sites by its Hermitian part
C- This routine is only called for nspc=2
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dmat ncomp norb
Cio    Passed to:  *
Ci Inputs
Ci   ldmsit
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   dmat  :density matrix
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   08 Jun 14 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical ldmsit
      integer nbas,nl
      double complex dmat(nl**2,nl**2,nbas,2,2,4)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,norb
      complex(8), allocatable ::  dmatl(:,:,:,:,:)

C ... Pack dmat into s_site for non-CPA sites unless it is already there
C ... This part is only needed in the noncollinear case without SOC
      if (.not. ldmsit) then
        do ib = 1, nbas
          if (s_site(ib)%ncomp >= 2) cycle
          norb = s_site(ib)%norb
          allocate(dmatl(norb,norb,2,2,4))
          dmatl(1:norb,1:norb,:,:,:) = dmat(1:norb,1:norb,ib,:,:,:)
          call dpscop(dmatl,s_site(ib)%dmat,norb*norb*2*2*2*4,1,1,1d0)
          deallocate(dmatl)
        enddo
      endif

C ... Take the Hermitian part
      do ib = 1, nbas
        call dmahp(s_site(ib)%ncomp,s_site(ib)%norb,s_site(ib)%dmat)
      enddo
      end

      subroutine dmahp(ncomp,norb,dmat)
      implicit none
      integer ncomp,norb
      double complex dmat(norb,norb,2,2,4,ncomp)

      integer icomp,j
      complex(8) dml(norb,norb)

      do icomp = 1, ncomp
c       call yprm('dmat bef 11',3,dmat(1,1,1,1,icomp),0,norb,norb,norb)
c       call yprm('dmat bef 22',3,dmat(1,1,2,2,icomp),0,norb,norb,norb)
c       call yprm('dmat bef 12',3,dmat(1,1,1,2,icomp),0,norb,norb,norb)
c       call yprm('dmat bef 21',3,dmat(1,1,2,1,icomp),0,norb,norb,norb)
        do j = 1, 4
          dml = dconjg(transpose(dmat(:,:,1,1,j,icomp)))
          dmat(:,:,1,1,j,icomp) = (dmat(:,:,1,1,j,icomp) + dml)/2
          dml = dconjg(transpose(dmat(:,:,2,2,j,icomp)))
          dmat(:,:,2,2,j,icomp) = (dmat(:,:,2,2,j,icomp) + dml)/2
          dml = dconjg(transpose(dmat(:,:,1,2,j,icomp)))
          dmat(:,:,2,1,j,icomp) = (dmat(:,:,2,1,j,icomp) + dml)/2
          dmat(:,:,1,2,j,icomp) =
     .      dconjg(transpose(dmat(:,:,2,1,j,icomp)))
        enddo
c       call yprm('dmat aft 11',3,dmat(1,1,1,1,icomp),0,norb,norb,norb)
c       call yprm('dmat aft 22',3,dmat(1,1,2,2,icomp),0,norb,norb,norb)
c       call yprm('dmat aft 12',3,dmat(1,1,1,2,icomp),0,norb,norb,norb)
c       call yprm('dmat aft 21',3,dmat(1,1,2,1,icomp),0,norb,norb,norb)
      enddo

      end
