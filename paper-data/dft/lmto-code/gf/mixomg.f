      subroutine mixomg(s_ctrl,s_site,nspc,nzp,izp,beta,rms)
C- Mixes CPA Omega-array
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp ldomg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp domg omg
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:omgn
Cio    Passed to:  womg
Ci Inputs
Ci   nspc  :number of coupled spins
Ci   nzp   :number of z-points
Ci   izp   :if 0, mix all nzp z-points and write omega to files
Ci         :if >0, mix one z-point izp, don't write to files
Ci   beta  :mixing parameter
Co Outputs
Co   omega :written to file
Co   rms   :largest RMS mismatch for Omega among all DLM sites
Cu Updates
Cu   19 Aug 13 (Belashchenko) Added relativistic case (2 x 2)
Cu   22 Feb 13 Writes screening representation to omega file
Cu   19 Jan 12 (Belashchenko) Complete rewrite
Cu   08 Dec 08 (P. Larson) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed parameters:
      integer nspc,nzp,izp
      double precision beta,rms
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
C Local parameters:
      logical bittst
      integer nbasp,norb,ib,nth,ioff,LGAM,irep,lrel
      integer mode,nelts,nglob,stdo
      double precision rms1
      complex(8), pointer :: p_omg(:,:)
      complex(8), pointer :: wk(:)
      parameter (LGAM=128)
      integer mpipid,procid,master

      procid = mpipid(1)
      master = 0
      stdo = nglob('stdo')
      nbasp = s_ctrl%nbasp
      lrel = mod(s_ctrl%lrel,10)

      mode = 0
      if (s_ctrl%ldomg == 1) mode = 1
      rms = 0d0
      do  ib = 1, nbasp
        norb = s_site(ib)%norb
        nth = s_site(ib)%ncomp
        if (mode == 1) then
          p_omg => s_site(ib)%domg
        else
          p_omg => s_site(ib)%omg
        endif
        if (nth < 2) cycle
        if (izp == 0) then
          nelts = norb*norb*nspc*2*nzp
          allocate(wk(nelts)); call dpzero(wk,2*nelts)
          call dpscop(p_omg,wk,2*nelts,1,1,-1d0)
          call dpsadd(wk,s_site(ib)%omgn,2*nelts,1,1,1d0)
          call dpdot(wk,wk,2*nelts,rms1)
          call dpsadd(p_omg,wk,2*nelts,1,1,beta)
          deallocate(wk)
          rms1 = dsqrt(rms1/nelts*norb)
          if (rms1 > rms) rms = rms1
          if (procid == master) then
            irep = 1
            if (bittst(s_ctrl%lham,LGAM)) then
              irep = 2
              if (bittst(s_ctrl%lasa,512) .or. mod(s_ctrl%lrel,10) == 2) irep=3
            endif
C           write(stdo,500) ib,rms1
            call info2(30,0,0,' Mixed Omega for site%,3i with RMS %,3;3g',ib,rms1)
            call womg(100*(nspc-1)+10*irep+mode,s_site,ib,nzp)
          endif
        else
          nelts = norb*norb*nspc*2
          allocate(wk(nelts))
          ioff = (izp - 1) * 2*nelts
          call dpscop(p_omg,wk,2*nelts,1+ioff,1,-1d0)
          call dpsadd(wk,s_site(ib)%omgn,2*nelts,1,1+ioff,1d0)
          call dpdot(wk,wk,2*nelts,rms1)
          call dpsadd(p_omg,wk,2*nelts,1+ioff,1,beta)
          deallocate(wk)
          rms1 = dsqrt(rms1/nelts*norb)
          if (rms1 > rms) rms = rms1
        endif
      enddo

C 500 format(' Mixed Omega for site',i3,' with RMS ',g10.3)
      end
