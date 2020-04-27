      subroutine dlminit(s_ctrl,s_ham,s_pot,s_spec,s_site,ics,nbas,nzp)
C- Allocations and initializations for DLM
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nclass nclasp nccomp lrel ldlm
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  etrms
Cio    Elts passed:lncol entropy
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:thetcl dlmwt
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  ncomp name ncpa nthet iscpa nang beff
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:xcpa nthet beff
Cio    Passed to:  dlmwgts
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class
Co     Stored:     ncomp cpawt
Co     Allocated:  cpawt thet j0 tau bxc omg omgn domg gii gc gcu gcorr
Co                 dmat sfvrtx
Cio    Elts passed:thet
Cio    Passed to:  *
Ci Inputs
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   nbas  :size of basis
Ci   nzp   :number of energy points
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Aug 13 space for s_site%thet and s_pot%thetcl to hold two angles
Cu   18 Dec 12 Completed migration to f90 pointers
Cu   01 Nov 11 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nzp,ics(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer nclasp,nl,nccomp,ncomp,nelts,nt
      integer ib,ic,is,ioff,i,nspc,nglob
      double precision xx
      logical lso,bittst

      nl = s_ctrl%nl
C     nclass = s_ctrl%nclass
      nclasp = s_ctrl%nclasp
      nspc = nglob('nspc')      ! 2 for relativistic case
      nccomp = s_ctrl%nccomp
      lso = (mod(s_ctrl%lrel,10) == 2 .or. bittst(s_ham%lncol,4))
      if (nccomp > 0) call dlmwgts(nclasp,nccomp,s_spec,ics,lso,
     .  s_pot%thetcl,s_pot%dlmwt,s_ham%entropy)

      call ptr_ham(s_ham,8+1,'etrms',(nclasp+nccomp)*22,0,xx)

      nelts = (nl*nl)*(nl*nl)*2*nspc
      do  ib = 1, nbas
        is = s_site(ib)%spec
        ncomp = s_spec(is)%ncomp
        if (ncomp < 2) ncomp = 1
        s_site(ib)%ncomp = ncomp
        call ptr_site(s_site,8+1,'cpawt',ib,ncomp,0,xx)
        if (ncomp == 1) s_site(ib)%cpawt(1) = 1d0
        call ptr_site(s_site,8+1,'thet',ib,ncomp,nspc,xx)
        call ptr_site(s_site,8+1,'j0',ib,ncomp,nl*nl,xx)
        call ptr_site(s_site,8+1,'tau',ib,ncomp,nl*nl,xx)
        call ptr_site(s_site,8+1,'bxc',ib,3*ncomp,0,xx)
c       call ptr_site(s_site,8+1,'omg',ib,nelts,nzp,xx)
c       call ptr_site(s_site,8+1,'omgn',ib,nelts,nzp,xx)
c       if (mod(s_ctrl%ldlm/100,10) == 1) then
c         call ptr_site(s_site,8+1,'domg',ib,nelts,nzp,xx)
c       else
c         call ptr_site(s_site,8+1,'domg',ib,1,1,xx)
c       endif
C       Allocate twice as many elements for gc to include off-diagonal
        call ptr_site(s_site,8+1,'gii',ib,nl**4,2*nspc,xx)
        call ptr_site(s_site,8+1,'gc',ib,nl**4,ncomp*4,xx)
        call ptr_site(s_site,8+1,'gcu',ib,nl**4,ncomp*4,xx)   ! for lrel=2 and exasa
        call ptr_site(s_site,8+1,'gcorr',ib,nl**2,ncomp*2,xx)
C       Component-specific density matrices
        call ptr_site(s_site,8+1,'dmat',ib,1,1,xx)
C       Placeholders to avoid problems with some compilers
        call ptr_site(s_site,1,'sfvrtx',ib,1,1,xx)
        if (ncomp == 1) cycle

C ...   Determine the offset in class arrays for this ib
        ic = s_site(ib)%class
        ioff = 1
        if (ic == 1) goto 5
        do  i = 1, ic-1
          nt = s_spec(ics(i))%ncomp
          if (nt > 1) ioff = ioff + nt
        enddo
    5   continue
        call dpscop(s_pot%thetcl,s_site(ib)%thet,ncomp,ioff,1,1d0)
        if (nspc == 2) then
          call dpscop(s_pot%thetcl,s_site(ib)%thet,ncomp,
     .      ioff+nccomp,1+ncomp,1d0)
        endif
      enddo

      end

      subroutine initomg(s_ctrl,s_site,nbas,nspc,nzp)
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nspc,nzp
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,norb,nelts
      real(8) xx

      do  ib = 1, nbas
        norb = s_site(ib)%norb
        nelts = norb*nspc*norb*2
        call ptr_site(s_site,8+1,'omg',ib,nelts,nzp,xx)
        call ptr_site(s_site,8+1,'omgn',ib,nelts,nzp,xx)
        if (mod(s_ctrl%ldlm/100,10) == 1) then
          call ptr_site(s_site,8+1,'domg',ib,nelts,nzp,xx)
        else
          call ptr_site(s_site,8+1,'domg',ib,1,1,xx)
        endif
      enddo

      end
