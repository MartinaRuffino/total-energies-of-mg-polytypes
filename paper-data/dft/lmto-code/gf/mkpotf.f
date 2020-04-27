      subroutine mkpotf(job,s_ctrl,s_spec,s_site,s_ham,s_pot,lrel,lso,
     .  pfdim,ipc,ipdc,vRLshf,vshft,zp,P0,P1,P2,P3,P4,P5,P6,P7)
C- Potential functions for a complex energy, for ASA GF, DLM version
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nbasp nspin nclass nspec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makpfz
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  a nr z rmt lmxa hcr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makpfz
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pnu spec class clabel v0 norb ncomp dlmcl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr dpfr ddpfr
Cio    Passed to:  makpfz mkptfp mkfrpf mksopf
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  makpfz
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pp pprel sop socscl dlmwt
Cio    Passed to:  makpfz mkptfp mksopf
Ci Inputs
Ci   job   :1s digit
Ci         : 0 second-order P, parameterized from pp
Ci         : 1 third-order P, parameterized from pp
Ci         : 2 first-order P, parameterized from pp
Ci         : 3 make P parameterized from direct integration of wave
Ci         :   equation at complex energy zp
Ci         :  (passed as 1s digit to mkptfp; also if 3, parameters
Ci             C,gam,del are made by integration of wf at energy zp).
Ci         :10s-1000s digit determines what to make:
Ci         : 2^0 P0 <- potential function P
Ci         : 2^1 P1 <- P-dot.  Note conflict 2^3: make either P1=P-dot or P1=sqrt(P-dot)
Ci         : 2^2 P2 <- -1/2 P-dotdot/P-dot
Ci         : 2^3 P1 <- sqrt(P-dot), choosing abs(sqrt(delta)).  See 2^1 option
Ci         : 2^4 P3 <- (-1/2 P-dotdot/P-dot)-dot
Ci         : 2^5 P5 <- potential function P^alp, stored in s_site(:)%pfra
Ci         :           (may be needed together with p^bet)
Ci         : 2^6 P4 <- P^alp/P^bet to convert g^alp -> g^bet, where bet = e.g. gamma
Ci         : 2^7 P5 <- P^alpha/P^gamma*(gamma-alpha) to convert g^alp -> g^bet
Ci         : Relativistic case, if gamma-bar repsn sought (10000s digit job = 3)
Ci         :     P4 <- P^alpha/P^gamma to scale 1/(P-S)^alp -> 1/(P-S)^gam
Ci         :     P6 <- (gamma-alpha) if work directly with S^gam (10000s digit>3)
Ci         :     P7 <- P^alpha/P^gamma (gamma-alpha) to scale 1/(P-S)^alp -> 1/(P-S)^gam
Ci         : Any combination of these may be taken.  Note, however,
Ci         : that different options write to P4 and P5.  Normally
Ci         : these kinds of functions are not both needed
Ci         :10000s digit Make P for particular screening transformation
Ci         : See Remarks.  Digit is passed as 100s digit to mkptfp
Ci         : 0 make P for beta=alpha
Ci         : 1 make P for beta=0 (bare representation)
Ci         : 2 make P for beta=gamma (gamma repsn).  Requires P^a/P^g and P^a/P^g (gam-alp)
Ci         : 3 Like 2, but beta=gamma-bar (gamma repsn averaged over spins)
Ci         : 4 Like 2, but S is scaled to S^gam, so P^alpha is not needed
Ci         : 5 Like 4, but beta=gamma-bar
Ci         :100000s digit applies when P computed from num. int. w.f.
Ci         :  (passed as 10s  digit to makpfz)
Ci         :0 Read potential from class file
Ci         :1 Use potential from site->ov0
Ci         :4 add 4 to 10s digit if pot->ves(ic) should be added
Ci         :  to spherical potential
Ci         :1000000s digit applies when P computed from num. int. w.f.
Ci         :  (passed as 100s  digit to makpfz)
Ci         :1 Use hcr = rmt, rather than spec->hcr
Ci         :  NB: 2nd gen LMTO uses this switch
Ci         :2 Envelope function has fixed k.e.=0 (as in 2nd gen LMTO)
Ci         :3 combination of 1+2
Ci         :10000000s digit applies when vRLshf is to be used
Ci         :0 vRLshf is not used
Ci         :1 RL-dependent shift of enu,C supplied by vRLshf
Ci         :2 Special mode for DLM (not used for now, reserved for later)
Ci   lrel  :0 for nonrelativistic
Ci         :1 for scalar relativistic
Ci         :2 for fully  relativistic
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ipdc  :DLM class index: first DLM class for site ib is ipdc(ib) (mksym.f)
Ci   pp    :potential parameters (atomsr.f)
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   vshft :array of site potential shifts
Ci   zp    :complex energy
Co Outputs
Co   These are returned in the beta repsn (beta=0,alp,gamma,or gamma-bar) in the following:
Co                               alpha reps        beta repsn, g^a->g^g     direct (P-S)^bet
Co                          SR   SR+SO,SR+CPA,FR   SR   SR+SO,SR+CPA,FR    SR   SR+SO,SR+CPA,FR
Co                          -------------------------------------------  ----------------------
Co   P^bet                  --                     P0   P0, s_site%pfr     P0   P0, s_site%pfr
Co   P^bet-dot              P1   P1, s_site%dpfr   P1   --                 P1   --
Co   sqrt(P-bet-dot)        --                     --   s_site%pfr         --   s_site%pfr
Co   -1/2 P-dotdot/P-dot    P2   P2, s_site%ddpfr  P2   s_site%ddpfr       P2   s_site%ddpfr
Co   energy deriv of ..     P3                     P3   --                 --   --
Co   P^alp/P^bet            --                     --   P4                 --   P4
Co   P^alp                  P0   P0, s_site%pfr    P5   s_site%pfra        P5
Co   bet-alp                --                     --   P6                 --   P6
Co   P^alp/P^bet (bet-alp)  --                     --   P7                 --   P7
Co
Co   Notes
Co    1.  Generic P0,P1,P2,... not meaningful in true CPA case with > 1 components / site
Co    2.  P0,P1,P2,... are stored in downfolding order, s_site%{pfr,dpfr,ddpfr} in orbital order
Co        P0    :potential function P
Co        P1    :P-dot or sqrt(P-dot)      ... for scaling g -> G
Co        P2    :-1/2 P-dotdot/P-dot       ... for scaling g -> G
Co        P3    :(-1/2 P-dotdot/P-dot)-dot ... for linear response
Co        P4    :P^alp/P^bet               ... for scaling g^alp -> g^gam
Co        P5    :P^alp                     ... when (P-S)^alp -> (P-S)^gam by energy scaling
Co        P6    :beta-alpha                ... to make S^gam
Co        P7    :P^alp/P^bet (bet-alp)     ... for scaling g^alp -> g^gam
Co    3.  right now P^alp and P^bet, but only one or the other is needed
Co        Note that the CPA s_pot%cp is always in bet.
Co    4.  s_site%pfr is in kappa-mu repsn in fully rel case, not in SO case for now
Co    5.  beta-alpha only returned in relativistic case, where it is a 2x2 spinor
Cr Remarks
Cr   To convert g to G, require :
Cr      sqrt(P-dot) (returned pdot in SR case, as sqrt(pdot) in REL=2 or SO cases)
Cr      -1/2 P-dotdot/P-dot
Cr   To scale g^alp => g^gam, if P-S is made in alp repsn, require:
Cr      P^alp
Cr      P^alp/P^bet
Cr      P^alp/P^bet (bet-alp)
Cr    ! P^bet is never needed
Cr   If (P^bet-S^bet) is made directly, require
Cr      P^bet
Cr      bet-alp
Cr    ! P^alp is never needed
Cu Updates
Cu   05 Jun 16 When making P^alp (1s digit = 2^5) write P to s_site(:)%pfra
Cu   06 Mar 16 Some redesign to improve consistency in various treatments
Cu             sqrt(P-dot) -> P1 (was P3) and P^bet/P^alp -> P4 (was P6)
Cu   10 Jan 16 Changes to enable making (P^bet-S^bet)^1 by rotating S^alp
Cu   00 xxx 13 (Kirill) mkpotf and mkpotf have been merged
Cu   15 May 08 (Kirill) First created from mkpotf to handle DLM case
Cu   16 Nov 07 Added vRLshf
Cu   09 Jun 04 (A Chantis) spin-orbit parameters for modification of enu
Cu             when P3 relativistic is constructed. Altered argument list.
Cu   18 Mar 03 (A Chantis) relativistic potential parameters.
Cu             Altered argument list.
Cu   21 Dec 01 First created
C     ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lso
      integer job,pfdim,ipc(*),ipdc(*),lrel
      double precision vRLshf(*),vshft(*)
      double complex zp,P0(*),P1(*),P2(*),P3(*),P4(*),P5(*),P6(*),P7(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer nl,nbasp,nsp,iopt,npfz,job1,iprint,lsgam
C     maps 10-1000's digit iopt into 10s digit
      procedure(logical) :: bittst
      integer ldim,lidim,lihdim
      parameter (npfz=8)
      double complex pfz(npfz,pfdim,2)

C     print *, '!! 162', job; call setpr(55)

      nl = s_ctrl%nl
      nbasp = s_ctrl%nbasp
      nsp = s_ctrl%nspin
      ldim = s_ham%ldham(1)
      lidim = s_ham%ldham(2)
      lihdim = s_ham%ldham(3)

      if (mod(job,10) == 3) then
        if (lrel == 2) call rx('mkpotf: not ready for lrel=2')
        iopt = (mod(job/100000,100))*10
        call pshpr(iprint()-20)
        call makpfz(iopt,s_ctrl,s_spec,s_site,s_ham,s_pot,lidim,lihdim,
     .    s_ham%iprmb,vshft,1,zp,pfz)
        call poppr
      endif

C     ...       order              repsn              effects of vRLshf
      iopt = mod(job,10) + 100*mod(job/10000,10) + 1000*mod(job/10000000,10)
      lsgam = 0   ! nonzero if to construct g = (P^gam - S^gam) instead of by energy scaling
      if (mod(job/10000,10) >= 4) then
        lsgam = mod(job/10000,10) - 3
        iopt = iopt - 200  ! gamma=4 -> gamma=1 and gamma=5 -> gamma=2
      endif

C --- Scaling parameters for gamma-repsn, relativistic or SO ppars
      if ((lrel == 2 .or. lso) .and. iopt/100 >= 2) then

C   ... P^alpha/P^beta, P^alpha/P^beta * (bet-alpha) to scale g^alp -> g^gam
        if (lsgam == 0) then
          call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,s_ham%iprmb,
     .      iopt+10040,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P4)

          call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,s_ham%iprmb,
     .      iopt+10050,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P7)
C   ... beta-alpha to convert S^alp -> S^gam
        else
          call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,s_ham%iprmb,
     .      iopt+10070,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P6)
        endif
      endif

C ... Potential function P
      job1 = mod(job/10,1000)
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,
     .  pfdim,s_ham%iprmb,iopt+0,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P0)

C ... P-dot
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,
     .  pfdim,s_ham%iprmb,iopt+10,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P1)

C ... -1/2 P-dotdot/P-dot
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,
     .  pfdim,s_ham%iprmb,iopt+20,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P2)

C ... Sqrt(P-dot)
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,
     .  pfdim,s_ham%iprmb,iopt+30,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P1)

C ... (-1/2 P-dotdot/P-dot)-dot
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,
     .  pfdim,s_ham%iprmb,iopt+60,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P3)

C ... P^alpha, stored in s_site(:)%pfra
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,
     .  s_ham%iprmb,iopt+0+10000-100*mod(iopt/100,10),s_pot%pp,s_pot%sop,s_pot%socscl,pfz,
     .  vRLshf,vshft,zp,P5)   ! Add 10000 to iopt to park P in s_site%pfra

C ... P-dot(bet)/P-dot(alpha)
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,
     .  s_ham%iprmb,iopt+40,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P4)

C ... beta-alpha
      job1 = job1/2
      if (mod(job1,2) /= 0) call mkptfp(s_site,s_pot,nl,nbasp,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,
     .  s_ham%iprmb,iopt+70,s_pot%pp,s_pot%sop,s_pot%socscl,pfz,vRLshf,vshft,zp,P6)

      end
