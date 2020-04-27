      subroutine suhamz(s_ctrl,s_site,s_ham,s_pot,s_spec,vRLshf,vshft,zp)
C- Set up potential for Green's function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lrel nl nbasp nspin lham lasa ldlm nclass nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb class dlmcl pnu spec clabel v0
Co     Stored:     *
Co     Allocated:  pfr dpfr ddpfr
Cio    Elts passed:pfr dpfr ddpfr
Cio    Passed to:  mkpotf makpfz mkptfp mkfrpf mksopf gvbma
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lncol hord ldham nlibu udiag
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  mkpotf makpfz
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves
Co     Stored:     *
Co     Allocated:  cp pf dpf ddpf dddpf pfr dpfr ddpfr gmar papg palp
Co                 gma
Cio    Elts passed: pf dpf ddpf dddpf palp pp gma pfr ddpfr dpfr papg
Cio                gmar pprel sop socscl dlmwt
Cio    Passed to:  mkpotf makpfz mkptfp mksopf gvbma
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkpotf makpfz
Ci Inputs
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   zp    :complex energy
Cio Inputs/ Outputs
Cio  spot  : pot. information is packed in structure spot; see Remarks
Cr Remarks
Cr   This routine sets up hamiltonian- and energy-dependent potential
Cr   parameters needed to create the hamiltonian for a Green's function.
Cr
Cr  *2nd generation LMTO returns potential functions in the following cases.
Cr   What form beta it takes depends on which repsn:
Cr   repsn beta
Cr     alpha       screening corresponds to repsn in which strux S^alp are made
Cr     gamma       orthogonal repsn.  In relativistic case, must distinguish:
Cr     gam (SR)    In relativistic case, gamma has a spinor form.
Cr                 a screening beta may be defined from scalar gamma
Cr                 This form is used when scaling 1/(P-S)^gam from 1/(P-S)^alp
Cr     gam (FR)    The full spinor form is used when making 1/(P-S)^gam from S^gam
Cr     gamma-bar   spin-averaged gamma repsn.  Either SR or spinor forms are possible
Cr     beta=0      no longer used.
Cr
Cr  *suhamz returns the following depending on which beta repsn and which relativistic treatment
Cr   Case A refers to scalar rel, B to SR with SO or CPA added, or fully rel
Cr                               alpha repsn                         beta repsn, g^a->g^g            direct (P-S)^bet
Cr                          A           B                        A           B                     A           B
Cr                          ----------- --------------           ----------------------------    ----------------------------
Cr   P^bet                  s_pot%pf    Note1,Note2              s_pot%pf    Note1,Note2           s_pot%pf    Note1,Note2
Cr   P^bet-dot              s_pot%dpf   Note1,Note2              s_pot%dpf   Note1,Note2           s_pot%dpf   Note1,Note2
Cr   sqrt(P^bet-dot)        Note1       Note2                    Note1       Note2                 Note1       Note2
Cr   -1/2 P-dotdot/P-dot    Note1       Note2                    Note1       Note2                 Note1       Note2
Cr   energy deriv of ..     s_pot%dddpf Not implemented          s_pot%dddpf Not implemented       s_pot%dddpf Not implemented
Cr   P^alp/P^bet            Note1       s_pot%papg               Note1       s_pot%papg            Note1       s_pot%papg
Cr   P^alp                  --          s_site%pfr               s_pot%palp  s_pot%palp            --          --
Cr   bet-alp                --          --                       s_pot%gma   s_pot%gmar            s_pot%gma   s_pot%gmar
Cr   P^alp/P^bet (bet-alp)  --          --                       Note1       s_pot%gmar            Note1       s_pot%gmar
Cr
Cr   Note1: in SR case sqrt(P^bet-dot), P^alp/P^bet, P^alp/P^bet*(bet-alp) are scalars made
Cr          on-the-fly from P^alp, P^bet, P^bet-dot, (bet-alp).
Cr          in the spin-coupled case these quantites are spinors and are stored explicitly
Cr   Note2: s_pot%pf     holds P (case A) s_pot%pfr (case B) vector form, ham order, lms.  But not made in SR+SO case!
Cr          s_site%pfr   holds P (case A and B) matrix form, orbital order.  If lrel=2, s_site%dpfr is in kappa-mu.
Cr          Similarly:
Cr          s_pot%dpf    holds Pdot (case A) and sqrt(Pdot) (case B) ham order, lms form  But not made in SR+SO case!
Cr          s_site%dpfr  holds Pdot (case A) and sqrt(Pdot) (case B) orbital order.  If lrel=2, s_site%dpfr is in kappa-mu.
Cr          Similarly:
Cr          s_pot%ddpf   holds -1/2 P-dotdot/P-dot (case A) s_pot%ddpfr (case B) ham order, lms form  But not made SR+SO!
Cr          s_site%ddpfr holds -1/2 P-dotdot/P-dot (case A and B) orbital order.  If lrel=2, s_site%ddpfr is in kappa-mu.
Cr   Note3: s_pot%gma is real, or SR.  s_pot%gmar is complex, (2x2) spinor form in SR+SO, FR cases
Cr   Note4: To convert g to G, require :
Cr          sqrt(P-dot) (see Note2)
Cr          -1/2 P-dotdot/P-dot
Cr   Note5: To generate g^bet, either ... g^alp=1/(P-S)^alp scaled=> g^bet. Require
Cr            P^alp
Cr            P^alp/P^bet
Cr            P^alp/P^bet (bet-alp)
Cr            (Note P^bet is never needed)
Cr                                 or ... (P^bet-S^bet) is made directly. Require
Cr            P^bet
Cr            bet-alp
Cr            (Note P^alp is never needed)
Cr          Plans for Future:
Cr          (1) ALWAYS make sqrt(Pdot)
Cr          (2) s_pot%pf and s_pot%pfr etc shouldn't be distinct.  Even better:
Cr          (3) Should no longer require s_pot%pf or s_pot%dpf or s_pot%ddpf or there relativistic counterparts
Cr          To remove this redundancy, creae a mapping from s_site to vector form, ham order, lms repsn on-the-fly
Cr          (4) Should need EITHER P^bet or P^alp but not both.
Cr          Problem : mkcpa which copies P to array s_pot%cp, in SR, SO, and CPA cases, always copies P^bet.
Cr          Solution: copy P^alp when relevant (but much checking/fix further dependencies in s_pot%cp)
Cr          (5) Resolve whether s_site%pfr, s_site%dpfr, s_site%dpfr be consistently in kappa-mu or lms repsn?
Cu Updates
Cu  05 Jun 16 When making P^alp write P^alp to s_site(:)%pfra when P^bet written s_site(:)%pfr
Cu  06 Mar 16 Some redesign to partially unify the various forms P-functions take
Cu  10 Jan 16 Changes to enable making (P^bet-S^bet)^1 by rotating S^alp
Cu  11 Feb 12 (Belashchenko) Additions for DLM case
Cu  10 Nov 11 Begin migration to f90 structures
Cu  28 Nov 07 NC branch now uses avgd gamma rep only when specified
Cu  08 Nov 07 (J. Xu) Additions for LDA+U implementation
Cu  18 Mar 03 (A Chantis) relativistic potential parameters.
Cu  21 Dec 01 Potential functions generated by new mkpotf.f
Cu   9 Nov 99 added generation of dddpf for linear response.
Cu   2 Dec 99 pass through vshft explicitly
Cu  21 Dec 99 set up potential parameters for gamma-representation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision vRLshf(*),vshft(*),zp(2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated arrays
      integer,pointer:: ipc(:),ipdc(:)
C ... Local parameters
      logical bittst,lso
      integer LGAM,hord,iopt,lalloc,lasa,ldham(16),ldim,lham,ldlm,lidim,
     .  lihdim,pfdim,lrel,nbasp,ib,nl,nsp,nspc,ludiag,modbma,ncomp,norb,lncol
      integer opfl
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (pfdim,ldham(9))
      parameter (LGAM=128)
      equivalence (nspc,ldham(4))
      double precision zps(2),xx
      procedure(integer) :: parg,bitand
C     For LDA+U
      integer nlibu
      save zps,lalloc
      data lalloc /0/ zps /0d0,0d0/

C ... Early return when potential setup is already complete
      if (lalloc /= 0 .and. zps(1) == zp(1) .and. zps(2) == zp(2)) return

      zps(1) = zp(1)
      zps(2) = zp(2)
      lalloc = 1
      lrel = mod(s_ctrl%lrel,10)
      lncol = s_ham%lncol
      lso = bittst(lncol,4)

      nl = s_ctrl%nl
      nbasp = s_ctrl%nbasp
      nsp = s_ctrl%nspin
      lham = s_ctrl%lham
      lasa = s_ctrl%lasa
      ldlm = s_ctrl%ldlm
c     if (ldlm /= 0 .and. lrel == 2)
c    .  print *, 'WARNING (suhamz): lrel=2 and CPA'
      hord = s_ham%hord
      ldham = s_ham%ldham
      nlibu = s_ham%nlibu
      ludiag = s_ham%udiag

C ... Various initializations.  Layer code has no CPA for now.
C     if (ldlm /= 0) call ptr_pot(s_pot,8+1,'cp',lihdim*lihdim*nsp*nspc,0,xx)
      if (ldlm /= 0 .and. s_ctrl%lpgf(1) == 0) call ptr_pot(s_pot,8+1,'cp',lihdim*lihdim*nsp*nspc,0,xx)
      allocate(ipc(nbasp))
      call sitepack(s_site,1,nbasp,'class',1,ipc,xx)
      allocate(ipdc(nbasp))
      call sitepack(s_site,1,nbasp,'dlmcl',1,ipdc,xx)

      iopt = 0
      if (hord == 3) iopt = 1
      if (hord == 1) iopt = 2
      if (hord == 4) iopt = 3 + 3400000

      call ptr_pot(s_pot,8+1,'pf',pfdim,nsp,xx)
      call ptr_pot(s_pot,8+1,'dpf',pfdim,nsp,xx)
      call ptr_pot(s_pot,8+1,'ddpf',pfdim,nsp,xx)
      call ptr_pot(s_pot,8+1,'dddpf',pfdim,nsp,xx)
C     Potential functions for these and CPA branches also stored in s_site(:)%pfr and related
C     s_site(ib)%pfr etc stored as norb*2 x norb*2 matrix
      if (lrel == 2 .or. lso) then ! pot functions stored as 2x2 matrices
        call ptr_pot(s_pot,8+1,'pfr',pfdim,4,xx)
        call ptr_pot(s_pot,8+1,'dpfr',pfdim,4,xx)
        call ptr_pot(s_pot,8+1,'ddpfr',pfdim,4,xx)
        call ptr_pot(s_pot,8+1,'papg',pfdim,4,xx)
        call ptr_pot(s_pot,8+1,'gmar',pfdim,4,xx)
        do ib = 1, nbasp
          ncomp = s_site(ib)%ncomp ; norb = s_site(ib)%norb
          if (associated(s_site(ib)%pfra,s_site(ib)%pfr)) then
            nullify(s_site(ib)%pfra)
          endif
          call ptr_site(s_site,8+1,'pfr',ib,norb*norb*4,ncomp,xx)
          call ptr_site(s_site,8+1,'dpfr',ib,norb*norb*4,ncomp,xx)
          call ptr_site(s_site,8+1,'ddpfr',ib,norb*norb*4,ncomp,xx)
          if (bittst(lham,LGAM)) then
            call ptr_site(s_site,8+1,'pfra',ib,norb*norb*4,ncomp,xx)
          else
            s_site(ib)%pfra => s_site(ib)%pfr
          endif
        enddo
        s_ctrl%lasa = s_ctrl%lasa + 2048 - bitand(s_ctrl%lasa,2048)
      else
        do ib = 1, nbasp
          ncomp = s_site(ib)%ncomp ; norb = s_site(ib)%norb
          if (associated(s_site(ib)%pfra,s_site(ib)%pfr)) nullify(s_site(ib)%pfra)
          call ptr_site(s_site,8+1,'pfr',ib,norb*nsp,ncomp,xx)
          call ptr_site(s_site,8+1,'dpfr',ib,norb*nsp,ncomp,xx)
          call ptr_site(s_site,8+1,'ddpfr',ib,norb*nsp,ncomp,xx)
          if (bittst(lham,LGAM)) then
            call ptr_site(s_site,8+1,'pfra',ib,norb*nsp,ncomp,xx)
          else
            s_site(ib)%pfra => s_site(ib)%pfr
          endif
        enddo
      endif

C ... Potential functions, beta repsn where beta = some form of gamma
C     lasa,128  => some gamma repsn
C     lasa,512  => turn on spin-avg gamma
C     lasa,1024 => gamma via S^alp->S^gam
      if (bittst(lham,LGAM)) then
C       Make in beta repsn:  P, P-dot,  Pdotdot/Pdot, (Pdotdot/Pdot)dot, plus P^alp
        iopt = iopt + 10*(2**0 + 2**1 + 2**2 + 2**4 + 2**5) + 20000
C        if ((lrel /= 2 .and. lso) .and. iand(lasa,512+1024) == 0) then
C          call rx('suhamz: GAMMA must be one of 0,2,4,5 in spin-coupled case, sorry')
C        endif
        if (ldlm /= 0 .and. s_ctrl%nccomp /= 0 .and. iand(lasa,1024) == 0) then
          call rx('DLM requires OPTIONS_GAMMA= 0, 4 or 5, sorry')
        endif
        modbma = 2
        if (bittst(lasa,512)) then  ! Use spin-averaged gamma
          iopt = iopt + 10000       ! P for gamma representation averaged over spins
          modbma = modbma + 1
        endif
        if (nlibu > 0) iopt = iopt + 10000000 ! Add LDA+U to shift C

        if (lrel == 2 .or. lso) then
C         Switches to make the following [required if g made through scaling (P-S)^alp]
C         P, sqrt(P-dot), -1/2 P-dotdot/P-dot, P-alpha/P, P^alpha/P^gamma*(gamma-alpha)
          iopt = iopt + 10*(2**3 - 2**1 - 2**4) ! Add sqrt(Pdot), subtract Pdot, (Pdotdot/Pdot)dot
          call ptr_pot(s_pot,8+1,'palp',pfdim,4,xx)  ! Shouldn't be needed in future ... use s_site%pfra
C         Modify switches when g made directly through (P-S)^bet
          if (bittst(lasa,1024)) then
            iopt = iopt + 20000 ! make bet-alp instead of P^alp/P, P^alp/P*(bet-alp)
            iopt = iopt - 10*(2**5) ! Don't make P^alp
          endif
          call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,lrel,lso,pfdim,ipc,ipdc,vRLshf,vshft,zp,
     .      s_pot%pfr,s_pot%dpfr,s_pot%ddpfr,xx,s_pot%papg,s_pot%palp,s_pot%gmar,s_pot%gmar)
C       Scalar relativistic parameters w/out SO coupling
        else
C         Make in beta repsn P, P-dot,  Pdotdot/Pdot, (Pdotdot/Pdot)dot, plus P^alp, gma
C         In this simpler case, P^alp/P, P^alp/P*(bet-alp) can be made on the fly, given P^alp, gma
          call ptr_pot(s_pot,1,'palp',pfdim,nsp,xx)
          call ptr_pot(s_pot,1,'gma',lihdim,nsp,xx)
          call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,vRLshf,vshft,
     .      zp,s_pot%pf,s_pot%dpf,s_pot%ddpf,s_pot%dddpf,xx,s_pot%palp,[xx],[xx])
          call gvbma(s_site,s_pot,modbma,nl,nsp,nbasp,lihdim,s_ham%iprmb,s_pot%pp,s_pot%gma)
        endif

C ... Potential functions, alpha repsn, downfolding order
      else
C       Make P, P-dot, -1/2 P-dotdot/P-dot, (-1/2 P-dotdot/P-dot)-dot
C       The last is used for linear response
        iopt = iopt + 10*(2**0 + 2**1 + 2**2 + 2**4)
        if (nlibu > 0) then
          iopt = iopt + 10000000
          if (ludiag == 0) call rx('GF-LDA+U not implemented with UDIAG=0 and GAMMA=0')
        endif

        if (lrel == 2 .or. lso) then
C         Make P, sqrt(P-dot), -1/2 P-dotdot/P-dot
          iopt = iopt - 10*(2**1+2**4) + 10*(2**3)
          call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,lrel,lso,pfdim,ipc,ipdc,
     .      vRLshf,vshft,zp,s_pot%pfr,s_pot%dpfr,s_pot%ddpfr,xx,xx,xx,xx,xx)
        else
C         Make P, P-dot, -1/2 P-dotdot/P-dot, (-1/2 P-dotdot/P-dot)-dot
          call mkpotf(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,0,.false.,pfdim,ipc,ipdc,
     .      vRLshf,vshft,zp,s_pot%pf,s_pot%dpf,s_pot%ddpf,s_pot%dddpf,xx,xx,xx,xx)
        endif
        if (lrel == 2 .or. lso) then ! Should be just a soft link!
          call ptr_pot(s_pot,4+1,'palp',pfdim,4,s_pot%pfr)
        else
          call ptr_pot(s_pot,4+1,'palp',pfdim,nsp,s_pot%pf)
        endif
      endif

      deallocate(ipc,ipdc)
C     stop

C      call zprm('palp',2,s_pot%palp,lihdim,lihdim,nsp)
C      call zprm('pfr',2,s_pot%pfr,lihdim,lihdim,4)

      return

      entry clhamz(opfl)
      lalloc = 0
      end
