      subroutine symcp(optrot,cp,lidim,nspc,nbas,nl,offH,idxsh,nsgrp,
     .  istab,g,ag,nsp,kcplx)
C- Symmetrize the coherent potential or the Green's function
C ----------------------------------------------------------------------
Ci Inputs
Ci   optrot:rotation mode for roth
Ci   lidim :dimension of cp
Ci   nl    :maximum l
Ci   nbas  :number of basis sites
Ci   offH  :Hamiltonian offsets (for roth)
Ci   idxsh :Hamiltonian permutation table (for roth)
Ci   nsgrp :number of symmetry operations
Ci   istab :table of site permutations for group op (mksym.f,symtbl.f)
Ci   g     :rotation part of space group operations
Ci         :NB: under suitable conditions space group may be extended
Ci         :by the addition of the inversion operation; see mksym.f
Ci   ag    :translation part of space group
Ci   nsp   :number of spins to process
Ci   kcplx :complex storage mode for cp
Cio Inputs/ Outputs
Cio  cp    :if kcplx=1, coherent potential in double complex format
Cio        :if kcplx=0, Green's function with real/imag parts separated
Cr Remarks
Cr   This routine performs the symmetrization of a matrix over the
Cr   crystal symmetry group. cp is either the crystal coherent potential
Cr   (called from mkcpa) or the crystal GF (called from gfibz)
Cr   The declaration cp(lidim,lidim,*) is formally correct for cp
Cr   but not for GF. However, since for GF there is only one spin,
Cr   this declaration causes no trouble (indices are always one in this
Cr   case).
Cu Updates
Cu   09 May 15 Redesigned so noncollinear case uses proper magnetic symmetry
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   21 Jan 12 (Belashchenko) Now also treats the Green's function
Cu   12 Jan 12 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer optrot,lidim,nl,nbas,nsgrp,istab(nbas,*),offH(*),idxsh(*)
      integer nsp,nspc,kcplx
      double precision g(3,3,*),ag(3,*)
      double complex cp(lidim,nspc,lidim,nsp)
C ... Dynamically allocated local arrays
      complex(8), allocatable, dimension(:,:,:,:) :: wk,wk2
C ... Local parameters
      double complex alpha,beta
      double precision q(3),pos(3,nbas),xx
      integer ig,isp,l2,rothnc,lidimx

      alpha = dcmplx(1d0/nsgrp,0d0)
      beta = dcmplx(1d0,0d0)
      l2 = lidim*lidim
      call dpzero(q,3)
      call dpzero(pos,3*nbas)
      lidimx = lidim*nspc

      allocate(wk(lidim,nspc,lidim,nspc),wk2(lidim,nspc,lidim,nspc))

C ... Convert to double complex
      do isp = 1, nsp
        call ztoyy(cp(1,1,1,isp),lidimx,lidimx,lidimx,lidimx,kcplx,1)
        if (nspc == 2) exit
      enddo
C     call zprm('cp before symmetrization',2,cp,lidimx,lidimx,lidimx)

C ... Symmetrize using c*16 format
      do isp = 1, nsp
        call dpcopy(cp(1,1,1,isp),wk,1,2*size(wk),1d0/nsgrp)
        do  ig = 2, nsgrp
          if (rothnc(100+optrot,nl,nspc,nbas,pos,xx,offH,idxsh,istab(1,ig),
     .      g(1,1,ig),ag(1,ig),q,lidim,lidim,cp(1,1,1,isp),wk2) < 0)
     .      call rx('SYMCP: error in rothnc')
          call daxpy(2*size(wk),1d0/nsgrp,wk2,1,wk,1)
        enddo
        call dpcopy(wk,cp(1,1,1,isp),1,2*size(wk),1d0)
        if (nspc == 2) exit
      enddo
      deallocate(wk,wk2)
C     call zprm('cp after symmetrization',2,cp,lidimx,lidimx,lidimx)

C ... Restore kcplx format
      do isp = 1, nsp
        call ztoyy(cp(1,1,1,isp),lidimx,lidimx,lidimx,lidimx,1,kcplx)
        if (nspc == 2) exit
      enddo

      end
