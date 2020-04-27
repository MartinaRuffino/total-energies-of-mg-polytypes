      integer function nnrl(mode,ib1,ib2,iprmb,ndim)
C- Returns number of RL channels between sites ib1, ib2
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0 number of RL channels
Ci          1 number of Rl channels (contracted over m)
Ci          2 number of sites containing a basis (contracted over l and m)
Ci         10s digit
Ci          0 count number of channels in ib1..ib2
Ci          1 find maximum number of channels for sites ib1..ib2
Ci   ib1   :starting site index
Ci   ib2   :ending site index
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndim  :cutoff for iprmb: orbitals for which iprmb<ndim are included
Co Outputs
Co   nnrl  :Number of RL or Rl channels, depending on mode.
Cr Remarks
Cu Updates
Cu   10 Jan 09 Handles 1s digit mode=2
Cu   29 Jun 00 Handles multiple-kappa case; also added 10s digit mode
C ----------------------------------------------------------------------
      implicit none
      integer mode,ib1,ib2,iprmb(*),ndim
C Local
      integer nglob,lmri,ib,nnrli,li,nl,nkaph,ik,mode0,nnrlx
      logical lmark
      save nl,nkaph
      data nl /-1/, nkaph /-1/

C     mxorb = nglob('mxorb')
      if (nl == -1) nl = nglob('nl')
      if (nkaph == -1) nkaph = nglob('nkaph')
      mode0 = mod(mode,10)

      nnrli = 0
      nnrlx = 0
      do  ib = ib1, ib2
        lmark = .false.
        lmri = nl*nl*nkaph*(ib-1)
        do  ik = 1, nkaph
        do  li = 0, nl-1
          lmri = lmri + 2*li+1
          if (iprmb(lmri) > ndim) cycle
          if (mode0 == 0) then
            nnrli = nnrli + 2*li+1
          elseif (mode0 == 1) then
            nnrli = nnrli + 1
          elseif (.not. lmark) then
            nnrli = nnrli + 1
            lmark = .true.
          endif
          enddo
        enddo
        if (mode >= 10) then
          nnrlx = max(nnrlx,nnrli)
          nnrli = 0
        endif
      enddo

      if (mode >= 10) then
        nnrl = nnrlx
      else
        nnrl = nnrli
      endif

      end
