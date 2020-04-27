      subroutine qmix(lc,nbas,nsp,nlmq,ltb,dclabl,ipc,it,itmax,cnvg,
     .                mmix,nkill,neltst,beta,tm,tj,qmpol,qits,
     .                mmom,mmom0,qdits,a,rms,wc,broy,wt)
C- Mixing multpole moments for TB-L (l-independent U)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lc: switch, if T mix only l=0 monopole moments
Ci   nbas,nsp,ltb,it,itmax,cnvg,mmix,nkill,beta,tm
Ci   qmpol: multipole moments from tbmpol
Ci   mmom: magnetic moments by site
Ci   qits, qdits, a: work arrays
Co Outputs:
Co   qmpol (and mmom) mixed, rms
Cr Remarks
Cr   If U is independent of l (but site dependent), i.e., Uav=T then
Cr   the most efficient is to mix the Q_RL and magnetic moments.
Cr
Cr   Mixing works like this:
Cr   At some iteration, i, we have a density (represented here simply
Cr   by multipoles, rather like the ASA) rho_i. This makes a hamiltonian
Cr   and the Schrodinger equation (bndtb, tbfrce) acts like a black box
Cr   to produce a function f(rho_i). The Anderson mixing consults all
Cr   previous rho_i and f(rho_i) and proposes a new density rho_i+1
Cr   which is fed back into the black box until rho_i+i is equal to
Cr   within an rms difference to rho_i.
Cr
Cr   qits is a work array which holds all previous rho_i, f(rho_i) and
Cr   rho_i+1, keeping just the previous mmix iterations. At each entry
Cr   to qmix these are rolled back in the array and iteration mmix+1 is
Cr   discarded so that after rolling the array is structured like this:
Cr
Cr                  i=0     i=1    2    3   ...
Cr   qits(..,i,1)  empty*   rho_i from previous iterations**
Cr   amix notation          x_0   x_1  x_2 ...
Cr
Cr    * mixed Q_RL will be put in here after amix
Cr   ** note i=1 is most recent, i=2 next most recent, etc.
Cr      if it=1 then for i=1 these are zero or Q_RL from disc (restart)
Cr
Cr   qits(..,i,2)  i=0     i=1    2    3   ...
Cr                 empty*   f(rho_i) from previous iterations**
Cr   amix notation f(x_0) f(x_1) f(x_2) ...
Cr
Cr    * the most recent rho_i from tbmpol is copied into here before
Cr      mixing
Cr   ** note i=1 is most recent, i=2 next most recent, etc.
Cr
Cr   amix requires as input f(x_i) and d(x_i) = f(x_i) - x_i
Cr   hence the call to daxpy in constructing the work array, a.
Cr
Cr   If nsp=2, then qits(1,.) holds just the spin up charge, the spin
Cr   down charge is maintained in qdits and tacked onto the
Cr   end of the mixing work array, "a".
Cr
Cb Bugs
Cb   No doubt keeping both qits and a is redundant, maybe someone clever
Cb   can save memory by getting rid of qits and building "a" directly
Cb   from qmpol.
Cu Updates
Cu   DMT cleanup a couple of things (Apr 2016)
Cu   ELS and ATP added Broyden mixing (Feb 2016)
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      logical lc
      integer nbas,nsp,nlmq,ltb,it,itmax,mmix,nkill,ipc(*),neltst,broy
      double precision qmpol(nlmq,nbas),qits(nlmq,nbas,0:mmix+1,2),
     .                 mmom(nbas),mmom0(nbas),qdits(nbas,0:mmix+1,2),
     .                 dclabl(1),a(neltst,0:mmix+1,2)
      double precision cnvg,beta,tm,tj(*),rms,wc,wt(3)
C Local Variables
      integer nelts,neltsm,nelts0,nmix,amix,npmix,i,ipr,iprint,ido,ib,
     .        ic,ierr,jump,i1mach
!       integer onorm,okpvt
      real(8), allocatable :: norm(:)
      integer, allocatable :: kpvt(:)
      double precision b,qmp(9),dqtot,d1mach,dabs,dsum,wctrue
      character clabl*8, outs*20
      logical IO, kill, bittst, cmdopt
C For MPI
      integer procid,master,mpipid
      logical mlog
C Local iteration count
      integer LOCIT
      save LOCIT

      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      IO = bittst(ltb,2**16)
      nelts = nlmq*nbas
      jump = 1
      if (lc) then
        nelts = nbas
        jump = nlmq
      endif
      neltsm = nbas
      nelts0 = nelts
      if (nsp == 2) nelts0 = nelts + neltsm
      call rxx(nelts0 /= neltst,' QMIX: bug, nelts0 /= neltst')

C --- "zero'th iteration" ---
      if (it == 0) then
        ierr = -1
        if (IO) then
          if (procid. eq. master) then
            call ioqm(nsp,nlmq*nbas,neltsm,qmpol,mmom,ierr)
          endif
          call mpibc1(ierr,1,2,mlog,'qmix','ierr')
          if (ierr /= -1) then
            call mpibc1(qmpol,nlmq*nbas,4,mlog,'qmix','qmpol')
            call mpibc1(mmom,neltsm,4,mlog,'qmix','mmom')
          endif
        endif
        if (ierr == -1) then
          call dcopy(nlmq*nbas,0d0,0,qmpol,1)
          if (nsp == 2) then
            call dcopy(nbas,mmom0,1,mmom,1)
          endif
        endif
        return
      endif

C --- kill "mix files" according to nkill or rms (see parmxp) ---
      kill = .true.
      if (cmdopt('--nomixkill',11,0,outs)) then
        kill = .false.
      endif

      if (it == 1) then
        LOCIT = 0
      endif
      if ( ( nkill < 0 .or. rms < 0d0 .or.
     .     ( nkill > 0 .and. mod(it,nkill) == 0 ) )
     .    .and. kill ) then
        LOCIT = 1
      else
        LOCIT = LOCIT + 1
      endif
      npmix = min(LOCIT-1,mmix)

C --- check for charge neutrality ---
      dqtot = dsum(nbas,qmpol,nlmq)
      if (dabs(dqtot) > d1mach(3)) then
        if (iprint() > 40 .or. (iprint() >= 30 .and.
     .      dabs(dqtot) > 1d-4)) then
          call awrit1(' QMIX: input qtot=%;2e',' ',120,i1mach(2),dqtot)
        endif
      endif

C --- Roll back previous iterations ---
      do  i = mmix, 0, -1
        call dcopy(nelts,qits(1,1,i,1),jump,qits(1,1,i+1,1),jump)
        call dcopy(nelts,qits(1,1,i,2),jump,qits(1,1,i+1,2),jump)
        if (nsp == 2) then
          call dcopy(neltsm,qdits(1,i,1),1,qdits(1,i+1,1),1)
          call dcopy(neltsm,qdits(1,i,2),1,qdits(1,i+1,2),1)
        endif
      enddo

C --- Copy new Q_RL and mmom for this iteration from input qmpol, mmom
      if (it == 1) then
        ierr = -1
        if (IO) then
          call pshprt(0)
          if (procid. eq. master) then
            call ioqm(nsp,nlmq*nbas,neltsm,qits(1,1,1,1),
     .                qdits(1,1,1),ierr)
          endif
          call mpibc1(ierr,1,2,mlog,'qmix','ierr')
          if (ierr /= -1) then
            call mpibc1(qits(1,1,1,1),nlmq*nbas,4,mlog,'qmix','qits')
            call mpibc1(qdits(1,1,1),neltsm,4,mlog,'qmix','qdits')
          endif
          call popprt
        endif
        if (ierr == -1) then
          call dcopy(nlmq*nbas,0d0,0,qits(1,1,1,1),1)
          if (nsp == 2) then
            call dcopy(neltsm,mmom0,1,qdits(1,1,1),1)
          endif
        endif
        if (nsp == 2) then
          call m2q(nbas,nlmq,qits(1,1,1,1),qdits(1,1,1))
        endif
      endif
      call dcopy(nlmq*nbas,qmpol,1,qits(1,1,0,2),1)
      if (nsp == 2) then
        call dcopy(neltsm,mmom,1,qdits(1,0,2),1)
        call m2q(nbas,nlmq,qits(1,1,0,2),qdits(1,0,2))
      endif

C --- Build work array for amix ---
      do  i = 0, npmix
        call dcopy(nelts,qits(1,1,i+1,1),jump,a(1,i,2),1)
        call dcopy(nelts,qits(1,1,i,2),jump,a(1,i,1),1)
        if (nsp == 2) then
          call dcopy(neltsm,qdits(1,i+1,1),1,a(1+nelts,i,2),1)
          call dcopy(neltsm,qdits(1,i,2),1,a(1+nelts,i,1),1)
        endif
        if (i /= 0 .and. broy == 0) then
          call daxpy(neltst,-1d0,a(1,i,2),1,a(1,i,1),1)
        endif
      enddo

C --- Mix; don't chatter about it ---
      if (broy == 0) then
        b = beta
        call pshprt(0)
        ipr = iprint()
        ido = 0
        allocate(norm(mmix*mmix))
        allocate(kpvt(mmix))
        nmix = amix(neltst,npmix,mmix,ido,b,ipr,tm,norm,kpvt,a,tj,rms)
        call popprt
        deallocate(kpvt,norm)
      else
        call qmixb(neltst,npmix,mmix,beta,wc,rms,a,wctrue,nmix)
      endif

C --- Get new Q_RL (and mmom) from work array ---
      call dcopy(nelts,a(1,0,2),1,qits(1,1,0,1),jump)
      call dcopy(nlmq*nbas,0d0,0,qmpol,1)
      call dcopy(nelts,a(1,0,2),1,qmpol,jump)
      if (nsp == 2) then
        call dcopy(neltsm,a(1+nelts,0,2),1,qdits(1,0,1),1)
        call dcopy(neltsm,a(1+nelts,0,2),1,mmom,1)
        call q2m(nbas,nlmq,qmpol,mmom)
      endif

C --- check for charge neutrality ---
      if (.not. cmdopt('--charged',9,0,outs)) then
        dqtot = dsum(nbas,qmpol,nlmq)
        if (dabs(dqtot) > d1mach(3)) then
          if (iprint() > 40 .or. (iprint() >= 30 .and.
     .      dabs(dqtot) > 1d-4)) then
            call awrit1(' QMIX: adding q=%;2e to conserve charge',' ',
     .        120,i1mach(2),-dqtot)
          endif
          dqtot = dqtot / nbas
          call daxpy(nbas,dqtot,-1d0,0,qmpol,nlmq)
        endif
      endif

C --- write moments to disc ---
      if (procid == master) then
        if (iprint() > 40) then
          print *, ' '
          print *,' QMIX: writing moments to disc..'
        endif
        ierr = 1
        call ioqm(nsp,nlmq*nbas,neltsm,qmpol,mmom,ierr)
      endif

C --- Printout ---
      if (iprint() < 10) return
      print 100
      if (broy == 1) then
        call awrit6(
     .  ' Broyden iteration %,3i. %i elements; mixed %i of %i, '
     .  //'rms diff: %g (tol: %g)',' ',90,i1mach(2),it,neltst,nmix,npmix,rms,cnvg)
      else
        call awrit6(
     .  ' Anderson iteration %,3i. %i elements; mixed %i of %i, beta=%d, '
     .  //'rms diff: %g',' ',90,i1mach(2),it,neltst,nmix,npmix,b,rms)
        if (nmix > 0) write (*,110) (tj(i),i=1,nmix)
      endif
      if (iprint() < 40) return
      do  ib = 1, nbas
        call dcopy(9,0d0,0,qmp,1)
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        call awrit1(' Atom %i '//clabl//'%cmultipole moments:',
     .        ' ',180,i1mach(2),ib)
        call dcopy(nlmq,qits(1,ib,1,1),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,1,1)
        call awrit3('        Q(in) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        call dcopy(nlmq,qits(1,ib,0,2),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,0,2)
        call awrit3('       Q(out) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        call dcopy(nlmq,qits(1,ib,0,1),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,0,1)
        call awrit3('     Q(mixed) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        if (nsp == 2) then
          call awrit3('    magnetic moment in %d, out %d, mixed %d',' ',
     .                120,i1mach(2),qits(1,ib,1,1)-qdits(ib,1,1),
     .                qits(1,ib,0,2)-qdits(ib,0,2),mmom(ib))
        endif
      enddo
  100 format(' QMIX mixing multipole moments:')
  110 format(' t_j :',10f8.4)
      end

      subroutine m2q(nbas,nlmq,qits,qdits)
C- Widget to exchange moments into up and down spin charges
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, nlmq, qits, qdits
Co Outputs:
Co   qits, qdits
Cr Remarks
Cr   This is really a bug fix. Previously we mixed qmpol and mmom; now
Cr   we want to mix up and down spin charges. This widget takes qits
Cr   and qdits which on input contain total charge and moments and
Cr   replace the relevant entries with up spin charge and down spin
Cr   charge.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nbas,nlmq
      double precision qits(nlmq,nbas),qdits(nbas)
C Local Variables
      integer ibas
      double precision q,m
      do  ibas = 1, nbas
        q = qits(1,ibas)
        m = qdits(ibas)
        qits(1,ibas) = 0.5d0*(q + m)
        qdits(ibas)  = 0.5d0*(q - m)
      enddo
      end

      subroutine q2m(nbas,nlmq,qits,qdits)
C- Widget to exchange up and down spin charges into moments
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, nlmq, qits, qdits
Co Outputs:
Co   qits, qdits
Cr Remarks
Cr   This is really a bug fix. Previously we mixed qmpol and mmom; now
Cr   we want to mix up and down spin charges. This widget takes qits
Cr   and qdits which on input contain up spin charge and down spin
Cr   charge and replace the relevant entries with total charge and
Cr   moments.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nbas, nlmq
      double precision qits(nlmq,nbas),qdits(nbas)
C Local Variables
      integer ibas
      double precision qu,qd
      do  ibas = 1, nbas
        qu = qits(1,ibas)
        qd = qdits(ibas)
        qits(1,ibas) = qu + qd
        qdits(ibas)  = qu - qd
      enddo
      end

      subroutine qmixb(neltst,npmix,mmix,beta,wc,rms,a,wctrue,jmix)
C- Broyden mixing of a vector, Duane Johnson's approach
C ------------------------------------------------------------------
Ci  npmix: number of iterates available to mix
Ci  a:    (*,i,1)  output values for prev. iteration i
Ci        (*,i,2)  input  values for prev. iteration i
Cio mmix: mmix > 0: number of iter to try and mix
Ci        mmix < 0: use npmix instead of mmix.
Co  npmix: (abs)  number of iter actually mixed.
Co        (sign) <0, intended that caller update npmix for next call.
Cr  Notations:
Cr  x^(m): input vector for iteration m
Cr  F^(m): difference between output and input vector in iteration m
Cu Updates
C ------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer neltst,npmix,mmix
      double precision beta,rms,wctrue,a(neltst,0:mmix+1,2)
C ... Dynamically allocated local arrays
      real(8), allocatable :: xmp1(:)
      real(8), allocatable :: dx(:)
      real(8), allocatable :: wk(:)
C ... Local parameters
      double precision ddot,dval,wc
      integer im,km,i,iprint,imix,stdo,nglob,broyj
      integer, intent(out) :: jmix
      real(8) :: dxnorm, neltstsq

C --- Setup ---

      allocate(xmp1(neltst))
      allocate(dx(neltst))
      neltstsq = sqrt(real(neltst,8))
C ... imix is a local copy of npmix
      imix = npmix
C --- Starting from iteration jmix, build the Jacobian matrix ---
   10 continue
      jmix = min(mmix,iabs(imix))
C     call defdr(owk,neltst*2*(jmix+2))
      allocate(wk(neltst*2*(jmix+2)))
      do  km = 1, jmix
C   ... this loops from most-distant to most-recent
        im = jmix-km+1
        call dcopy(neltst,a(1,im-1,1),1,dx,1)
        call daxpy(neltst,-1d0,a(1,im-1,2),1,dx,1)
        dxnorm = sqrt(ddot(neltst,dx,1,dx,1))
        rms = dxnorm/neltstsq
C ---   Determine wc_true if wc < 0 ---
        if (wc < 0) then
          wctrue = -0.01_8*wc/dxnorm
          wctrue = min(max(wctrue,1d0),1d4)
        else
          wctrue = wc
        endif
        if (km == 1) wctrue = .01d0

        i = iprint()
        if (km /= jmix) i = i-20
        i = broyj(neltst,a(1,im-1,2),dx,km,0,i,beta,0d0,0d0,0d0,wctrue,wk,neltst,xmp1)
      enddo

      deallocate(dx)

! C --- Check for interactive change of npmix ---
!       im = imix
! C      if (iprint() > 30) call query('redo, npmix=',2,imix)
!       if (iabs(imix) > mmix .and. imix /= im .and. iprint() > 30)
!      .  call awrit1(' (warning) only %i iter available',' ',80,mmix)
!       if (im /= imix) then
!         deallocate(wk)
!         goto 10
!       endif
      npmix = imix
C ... If no prior iter allowed, give up on npmix
C --- Save x^(m+2) into a(*,0,2) and exit ---
      if (npmix /= 0) call dcopy(neltst,xmp1,1,a(1,0,2),1)


      deallocate(xmp1)

      end
