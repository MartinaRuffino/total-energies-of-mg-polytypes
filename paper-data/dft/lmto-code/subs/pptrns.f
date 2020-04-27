      subroutine pp2alp(mode,s_ctrl,s_str,lham,pp)
C- Transform 2nd gen potential parameters to alpha, read from disk
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl:struct contains pointers to various arrays; see structures.h
Ci     Elts read: nclasp ipc
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:ncomp ipc idcc
Cio    Passed to: *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:alph
Cio    Passed to: *
Ci Inputs
Ci   mode  :0  transform 2nd gen pp's to alpha as read from STR file
Ci         :1  transform pp's to gamma
Ci     10s :0  normal call
Ci         :1  DLM case (only needed if mode = 0)
Ci   lham  :ASA hamiltonian parameters (see structures.h)
Ci         :1's bit => 2-center
Ci         :2's bit => 2-c + pert. corr
Co Outputs
Co   pp are transformed to new representation (alph,c,sqdel,gam,palph)
Cl Local variables
Cr Remarks
Cr   The DLM mode is designed to take care of the fact that the extra
Cr   DLM classes do not appear in s_ctrl%ipc, but they are
Cr   mapped to the parent classes by idcc array. In this case a new
Cr   array icp is set up, which maps classes to lattice sites. This
Cr   array is passed on to pptrns, which has a special mode designed
Cr   for that purpose. Non-DLM calls to pp2alp and pptrns are thereby
Cr   preserved intact.
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   18 Jan 12 (Belashchenko) DLM mode added, 100's digit iopt for DLM
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Oct 10 Read strux through rdstrx
Cu   24 Aug 09 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,lham
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_str)::   s_str
C ... Local parameters
      double precision pp(6,*)
      integer nclasp,nbasp,nl,nsp,nglob,iopt,mode0
      logical bittst
C ... Dynamically allocated local arrays
      integer,target,allocatable :: icp(:)
      integer, pointer :: p_ipc(:)
      real(8),allocatable :: oold(:,:,:)
      real(8), pointer :: alph(:)
C ... For DLM
      integer idlm,nangl,isum,iclbas,ic,iclbsj

      nl = nglob('nl')
      nsp = nglob('nsp')

      if (mode == 11) call rx('pp2alp: overkill mode 11')
      idlm = mod(mode/10,10)
      mode0 = mod(mode,10)

C ... Get screening alpha's
      if (mode0 == 0) then
        iopt = 0
        if (bittst(lham,1) .or. bittst(lham,2)) iopt = 10
      else
        iopt = 1
      endif

      nclasp = s_ctrl%nclasp
      nbasp = nglob('nbasp')
      if (idlm /= 0) then
        nangl = isum(nclasp,s_ctrl%ncomp,1)
c       if (nangl < 2) call rx('pp2alp: DLM but nangl<2')
        allocate(icp(nclasp+nangl))
        do  ic = 1, nclasp
          icp(ic) = iclbsj(ic,s_ctrl%ipc,-nbasp,1)
        enddo
        do  ic = nclasp+1, nclasp+nangl
          icp(ic) = iclbas(s_ctrl%idcc(ic),s_ctrl%ipc,
     .      size(s_ctrl%ipc))
        enddo
        nclasp = nclasp + nangl
        iopt = iopt + 100
        p_ipc => icp
      else
        p_ipc => s_ctrl%ipc
      endif

C ... Do the transformation
      allocate(oold(nl,nsp,nclasp))
      if (associated(s_str%alph)) then
        alph => s_str%alph
      else
        if (mod(iabs(iopt),10) == 0) call rx('pp2alp has to s_str%alph')
        allocate(alph(1))
      endif

      call pptrns(iopt,nl,p_ipc,nclasp,nsp,alph,nbasp,pp,oold)
      deallocate(oold)
      if (.not. associated(s_str%alph)) deallocate(alph)

      end

      subroutine pptrns(iopt,nl,ipc,nclass,nsp,alpha,nbas,pp,oold)
C- Transform the set of potential parameters into another repsn
C ----------------------------------------------------------------
Ci Inputs
Ci   iopt: 1s digit
Ci         0: to passed alpha
Ci         1: to alpha = gamma
Ci         2: to alpha=0
Ci         3: to alpha = gamma(spin1)
Ci         4: to alpha = gamma(spin2)
Ci         5: to alpha = (gamma(spin1) + gamma(spin2)/2
Ci         10s digit
Ci         1: set p(gamma) to zero
Ci         100s digit
Ci         0: ipc contains mapping site -> class
Ci         1: ipc contains mapping class -> site (iopt > 0 in this case)
Ci         NB: sign of iopt<0 used to flag returning new alpha in alpha
Ci   nl,nclass,nsp,eny
Ci   alpha,ipc,nbas (needed only for iopt<=0)
Ci   alph,c,sqdel,gam,palph (old)
Ci   nbas  : size of ipc; also used if iopt<0
Co Outputs
Co   pp are transformed to new representation (alph,c,sqdel,gam,palph)
Co   oold -- overlap in old alpha representation
Co   alpha   (iopt<0)
Cr Remarks
Cr   pp(1) : enu
Cr   pp(2) : calpha
Cr   pp(3) : srdel = sqrt(delta) but with proper sign (that of phi-).
Cr   pp(4) : palpha
Cr   pp(5) : gamma, or Q in Varenna notes
Cr   pp(6) : alpha, or qbar in Varenna notes
Cr  Transformations use the following (see Varenna, p88)
Cr    (1) srdnew/srdold = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
Cr  where srdnew,srdold are sqrt(delta) for new and old representations
Cr    (2) 1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
Cr  (calculated by an external function);
Cr  change in palpha = change in oalpha**2 since
Cr    (3)  p = o**2 + p^gam
Cr    (4)  p^gam = <phidot | phidot> where phidot is in the gamma representation
Cr  Other relations can be verified.  Unscripted C and delta => gamma repsn
Cr    (5)  1 / o^a = delta^a / (alpha - gamma) - (C^a - enu)
Cr                 = delta   / (alpha - gamma) + (C - enu),
Cr                 where delta, C are in the gamma representation
Cr    (6)  delta^a = (alpha - gamma)^2 / delta / o^2
Cr    (7)  C^a-enu = (alpha - gamma) (C-enu) / delta / o
Cb Bugs
Cb   iopt=-1 ALWAYS returns alpha(spin 2), though alpha differs for
Cb   spins 1,2
Cu Updates
Cu   18 Apr 05 New 10's digit iopt
C ----------------------------------------------------------------
      implicit none
      integer iopt,nl,nclass,nsp,nbas,ipc(*)
      double precision alpha(0:nl**2-1,nbas)
      double precision pp(6,nl,nsp,nclass),oold(nl,nsp,nclass)
      integer isp,jc,il,iclbas,jb,m,ib,nl2
      double precision xx,enu,gamma,pgam
      double precision cold,srdold,alpold,cnew,srdnew,alpnew,pold,pnew
      double precision oalpha
      external oalpha,iclbas
      logical lrev

C     call yprm('enter pptrns',1,pp,0,6,6,nl*nsp*nclass)

      lrev = mod(iabs(iopt)/100,10) == 1

      do  jc = nclass, 1, -1
      do  isp = 1, nsp
      do  il = 1, nl

        if (mod(iabs(iopt),10)+1 == 2) then
          alpnew = pp(5,il,isp,jc)
        elseif (mod(iabs(iopt),10)+1 == 3) then
          alpnew = 0
        elseif (mod(iabs(iopt),10)+1 == 4) then
          alpnew = pp(5,il,1,jc)
        elseif (mod(iabs(iopt),10)+1 == 5) then
          alpnew = pp(5,il,nsp,jc)
        elseif (mod(iabs(iopt),10)+1 == 6) then
          alpnew = (pp(5,il,1,jc)+pp(5,il,nsp,jc))/2
        else
          if (lrev) jb = ipc(jc)
          if (.not.lrev) jb = iclbas(jc,ipc,nbas)
          if (jb == 0) cycle
          alpnew = alpha((il-1)**2,jb)
        endif
        if (iopt < 0) then
          jb = iclbas(jc,ipc,nbas)
          if (jb == 0) cycle
          do  m = 1, 2*il-1
            alpha((il-1)**2+m-1,jb) = alpnew
          enddo
        endif

C --- Calculate potential parameters in new representation from old ---
        enu = pp(1,il,isp,jc)
        cold = pp(2,il,isp,jc)
        gamma = pp(5,il,isp,jc)
        alpold = pp(6,il,isp,jc)
        pold = pp(4,il,isp,jc)
        srdold = pp(3,il,isp,jc)

C   ... delta=0 => no potential parameters for this l channel
        if (alpnew == alpold) cycle
        if (srdold == 0) cycle

        xx = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
        srdnew = srdold*xx
        cnew = enu + (cold-enu)*xx

        oold(il,isp,jc) = oalpha(enu,cold,srdold**2,alpold,gamma)
        pgam = pold - oold(il,isp,jc)**2
        if (mod(iabs(iopt),100) >= 10) pgam = 0
C        pnew = pold - oold(il,isp,jc)**2 +
C     .         oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2
        pnew = pgam + oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2

C        print *, 'Verify relation (5): calculate 1/o calling oalpha and by both Eqns 5'
C        print *, 1/oalpha(enu,cnew,srdnew**2,alpnew,gamma)
C        print *, srdnew**2/(alpnew-gamma) - (cnew-enu)
C        print *, srdold**2/(alpnew-gamma) + (cold-enu)
C
C        print *, 'Verify relation (6)'
C        print *, srdnew, (alpnew-gamma)/oalpha(enu,cnew,srdnew**2,alpnew,gamma)/srdold
C
C        print *, 'Verify relation (7)'
C        print *, cnew-enu, (cold-enu)*(alpnew-gamma)/oalpha(enu,cnew,srdnew**2,alpnew,gamma)/srdold**2

        pp(2,il,isp,jc) = cnew
        pp(3,il,isp,jc) = srdnew
        pp(6,il,isp,jc) = alpnew
        pp(4,il,isp,jc) = pnew
      enddo
      enddo
      enddo

C --- If alpha is returned, copy alpha to all ib ---
      if (iopt < 0) then
        nl2 = nl*nl
        do  ib = 1, nbas
          jb = iclbas(ipc(ib),ipc,nbas)
          if (jb == 0) cycle
          call dpscop(alpha,alpha,nl2,nl2*(jb-1)+1,nl2*(ib-1)+1,1d0)
        enddo
      endif

C     call yprm('exit pptrns',1,pp,0,6,6,nl*nsp*nclass)

      end
      double precision function oalpha(enu,c,delta,alpha,gamma)
C- Calculate overlap in alpha representation from pp's
C ----------------------------------------------------------------
Ci Inputs
Ci   enu,c,delta,alpha,gamma
Co Outputs
Co   oalpha
Cr Remarks
Cr   Varenna p 88, Eq 91 has:
Cr   1 / o^a =  (C^gam - Enu) - delta^gam / (gamma - alpha)
Cr   More generally, it can be written:
Cr   1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision enu,c,delta,alpha,gamma
C Local parameters
      double precision xx

      xx = (alpha-gamma)/delta
      oalpha = xx/(1 - xx*(c-enu))

      end
