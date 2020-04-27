      subroutine gfibz(s_site,s_ctrl,s_ham,s_pot,nbas,isp,nsp,mode,pos,
     .  nlibu,lmaxu,idu,ludiag,vorb,nkp,qp,wtkp,nk1,nk2,nk3,ipq,plat,
     .  salp,iax,nttab,istab,g,ag,nsgrp,igstar,ifac,qb,gii,ghh)
C- Green's function, integrated over BZ
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl ldlm lrel
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lsig ldham lgen3 bandw lncol ndhrs neula qss offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula neula qss offH nprs iaxs hrs
Cio    Passed to:  gf1kp gfg2g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  palp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:palp cp gma papg gmar pf dpfr ddpfr dpf ddpf
Cio    Passed to:  gf1kp gfg2g
Cio  s_site :struct for site-specific data; see structures.h
Ci Inputs
Ci   nbas  :number of atoms in the basis
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1 (input)
Ci   mode  :1s digit
Ci          0 make g only at irr qp (kamikaze approach).
Ci            Only class-averaged quantities are ok, such as qnu.
Ci          1 generate only site diagonal part of g
Ci          2 generate entire g
Ci          3 calculate g at every nk1*nk2*nk3 point;
Ci            do not use any information about symmetry operations.
Ci         :10s digit handles what is made
Ci           1 make gii = integral over bz of g(qp)
Ci           2 save g(qp) to disk
Ci           4 save g(qp) to disk, rewinding file before first point
Ci         :100s digit handles special scaling of gf
Ci           0 g is unscaled
Ci           1 g is scaled to G by calling gfg2g
Ci         :1000s digit
Ci          0 normal behavior
Ci          1 record separate files for spin 1 and 2 if saving g(qp) to disk
Ci         :10000 digit not used
Ci         :100000 digit suppress noncollinear spinor rotation
Ci   pos   :basis vectors
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :list of irreducible k-points
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   salp  :real space screened structure constants
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :rotation part of space group operations
Ci         :NB: under suitable conditions space group may be extended
Ci         :by the addition of the inversion operation; see mksym.f
Ci   ag    :translation part of space group
Ci   igstar:contains info needed to map an irreducible qp to its
Ci          original point (bzmesh should be called with igstar(0)=-2)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone
Co Outputs
Ci   gii   :q-summed g (g_RR' for R,R' inequivalent sites in the basis)
Cs Command-line switches
Cs   --mem : Use disk to conserve memory
Cl Local variables
Cl   optrot:0 rotate whole g  1 rotate only diagonal part of g
Cl          Add 20 for if rotating to spherical harmonics.
Cl   mod1kp:mode passed to gf1kp
Cr Remarks
Cu Updates
Cu   18 Jun 18 Synchronize with updated spherical harmonics
Cu   12 Aug 13 Noncollinear symmetrization includes spinor rotations
Cu   18 Dec 12 Completed migration to f90 pointers
Cu   10 Nov 11 Begin migration to f90 structures
Cu   08 Nov 07 (J. Xu) Added LDA+U
Cu   15 Dec 04 (T. Sandu) Read sigma to add to hamiltonian
Cu   18 Mar 03 (A Chantis) first cut at relativistic case.
Cu              Altered argument list.
Cu   05 Feb 01 roth was changed into a function call and is capable
Cu             of handling inversion extensions to the space group.
Cu   15 Mar 00 redesigned the iwaves downfolding, added h-waves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nsp,nttab,nkp,nk1,nk2,nk3,niax,istab(nbas,*),
     .  ipq(*),igstar(0:*),ifac(3),nsgrp
      parameter (niax=10)
      integer iax(niax,nttab)
C     for LDA+U
      integer nlibu,lmaxu,idu(4,*),ludiag
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     double precision gii(lidimx,lidimx,2)
      double precision plat(3,3),salp(*),gii(*),wtkp(nkp),qp(3,nkp)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3),ghh(*)
C ... Dynamically allocated local arrays
      complex(8), pointer :: grot(:,:)
      complex(8), allocatable :: ghhl(:),ghhl1(:,:),ghhl2(:,:)
      complex(8), allocatable, target :: sll(:,:)
      integer, allocatable :: kpproc(:)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_site)::  s_site(*)
C ... Local parameters
C     logical :: debug=.false.
      logical lsph,bittst,cmdopt,reads,lgrotu
      integer PRTG,hdim,ig,iq,iq1,isp,k,l2,jj1,jj2,jj3,lbloch,ldham(16),ldim,lham,
     .  lhdim,lidim,lidimx,lix,lrel,lix2,lncol,mod1kp,mode0,mode1,mode2,mode3,
     .  mxorb,nfilet,nglob,nkfbz,nl,nlx,nspc,optrot,optrth,hhx,lrsig,ldlm,i1,i2,i3
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lhdim,ldham(3))
      equivalence (lidimx,ldham(6)),(nspc,ldham(4))
      character strn*80
      parameter (PRTG=90,nlx=5)
      double precision q1(3),qk,beta(2),alpha(2)
      procedure(integer) :: rotspd,fopnT,iprint,roth,rothnc

C ... For file I/O of gf
      integer clp(9,2),fopnx,ifi,iogfrs,kcplx,lrd
      character*8 fnam
C     Bits for selecting which data to I/O, iogfrs
      integer ISH
      parameter (ISH=1)
      double precision zp0(2),xx(1)
      integer mpipid,procid,master,numprocs
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
      call tcn('gfibz')
      procid = mpipid(1)
      numprocs = mpipid(0)
      master = 0
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
      ldlm = s_ctrl%ldlm
      if (nl > nlx) call rxi('gfibz: increase nlx to',nl)
C     200's bit of lbloch tells gf1kp to make -(S-P)
      lbloch = 200
C     Use spherical harmonics when making Bloch sum
      if (bittst(lham,256)) then
        lsph = .true.
        lbloch = lbloch + 1000
      else
        lsph = .false.
      endif
      lrel = mod(s_ctrl%lrel,10)
      lrsig = s_ham%lsig
C     Hamiltonian dimensioning and downfolding parameters
C     offH = s_ham%ooffH
      ldham = s_ham%ldham
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)

C ... Set mode for gf1kp
C     10s digit controls repsn and whether to scale g->G
      mod1kp = mode0 + 10*mode2 ! scale g to G
      if (ldlm /= 0) then
        if (s_ctrl%nccomp == 0) then
          mod1kp = mod1kp + 2000 ! CPA mode but with no CPA components
        else
          if (lrel == 2) call rx('check CPA in relativistic case')
          mod1kp = mod1kp + 1000 ! CPA mode with CPA components
        endif
      endif
      if (bittst(lham,128) .and. bittst(s_ctrl%lasa,1024)) then
        mod1kp = mod1kp + 40    ! g^gam directly by making S^gam from S^alp
      else if (bittst(lham,128)) then
        mod1kp = mod1kp + 20    ! g^alp -> g^gam
      endif

      optrot = 0000
      optrth = 33005
      if (mode0 == 1) optrot = 0001
      if (lsph) optrot = optrot + 20
      if (lsph) optrth = optrth + 20
C     Rotate iwaves block
      if (lidim > ldim) optrot = optrot + 4000
      optrot = optrot + 100000*mod(mode/100000,10) ! suppress symmetrization of spinors, nc case
C     Switch should not be added. A bit confusing; but confirm by
C     comparing with and w/out symops in e.g. FR fept
C     This is because rotsp1 and roth both generates passive rotations.  See rotsp1.
C     if (mod(optrot/100000,10) == 0) optrot = optrot + 200000

      kcplx = 0
C ... Allocate memory for g at this qp
      l2  = lidim**2
      lix = lidim*nspc
      lix2 = lix**2
      hdim = lhdim-lidim
      mxorb = nglob('mxorb')
      hhx =  hdim*nspc*mxorb
C ... Open nfilet if to conserve memory
      nfilet = 0
C     if (cmdopt('--mem',5,0,strn)) nfilet = fopna('TMP',-1,4)
      if (cmdopt('--mem',5,0,strn)) nfilet = fopnT('TMP',-1,4,0)
      allocate(sll(lix,lix),ghhl((max(hhx,1))))
      if (nspc==2 .and. nsgrp>1) allocate(grot(lix,lix))
C     call dvset(sll,1,lix**2*2,-99d0)
      call dpzero(ghhl,max(2*hhx,1))
      nkfbz = nk1*nk2*nk3
C     Parameters needed to accumulate g(1kp) into gii
      alpha(1) = 1d0/nkfbz
      alpha(2) = 0
      beta(1) = 1
      beta(2) = 0
C     No gf to make; just exit
      if (mode1 == 0) then
        goto 999
C     Setup, file opening to save g on disk
      else if (mode1 >= 2) then
        call iinit(clp,9*2)
        clp(3,1) = lix
        clp(4,1) = lix
        clp(5,1) = lix
        fnam = 'gfqp'
        if (mode3 == 1) then
          fnam = 'gfqp1'
          if (isp == 2) fnam = 'gfqp2'
        endif
        ifi = fopnx(fnam,100,16+8+4+0,-1)
C   ... rewind file, write header
        if (mode1 >= 4) then
          rewind ifi
          call pshpr(1)
          zp0(1) = 0
          zp0(2) = 0
          lrd = iogfrs(3,0,0,fnam,ifi,1,1,nbas,0,zp0,qp,plat,xx,clp,
     .      xx,xx,0,0,0,xx)
          call poppr
        endif
      endif

C --- For each irreducible qp, make G at each qp in star of q ---
      if (mode0 < 3) then
      if (numprocs > 1) call info0(50,1,0, ' ... Start MPI k-loop')
      allocate (kpproc(0:numprocs))
      if (numprocs > 1) then
        call pshpr(1)
        call dstrbp(nkp,numprocs,1,kpproc(0))
        call poppr
      else
        kpproc(0) = 1 ; kpproc(1) = nkp+1
      endif
C       print *, '!!'; qp(1,1) = .1d0; qp(2,1) = .2d0; qp(3,1) = .3d0
C       print *, '!!'; qp(1,1) = -.1d0; qp(2,1) = -.2d0; qp(3,1) = .3d0
      do  iq = kpproc(procid), kpproc(procid+1)-1

C       print *, '!!'; if (iq == 1) cycle
C   ... sll <- g at this irreducible qp
C       call debugstop(0)
        call gf1kp(s_site,s_ham,s_pot,mod1kp,lbloch,lrel,nl,nbas,isp,nsp,nspc,
     .    nlibu,lmaxu,vorb,idu,ludiag,lrsig,qp(1,iq),plat,salp,iax,
     .    nttab,lidim,hdim,sll,ghhl)
        grot => sll             ! for first symop (ig=1), rotated g = unrotated g
        lgrotu = .false.        ! lgrotu=T when spinors are rotated to global axis
C        print "(' after gf1kp sum sll for iq',i4,2f18.12)", iq, sum(dcmplx(sll))
C       call yprm('g(q)',2,sll,lix*lix,lix,lix,lix)
C   ... Save g at this irreducible qp on disk
        if (mode1 >= 2) then
          call pshpr(iprint()-10)
          lrd = iogfrs(10000*ISH+7,1,0,' ',ifi,1,1,0,0,zp0,
     .      qp(1,iq),xx,xx,clp,xx,xx,0,0,kcplx,sll)
          call poppr
        endif

C   ... If gii is not sought, skip over accumulation into gii
        if (mod(mode1,2) == 0) then
          cycle

C   ... Add to gii, kamikaze style
        else if (mode0 == 0) then
          call daxpy(lix2*2,abs(wtkp(iq))/2,sll,1,gii,1)
          call daxpy(hhx*2,abs(wtkp(iq))/2,ghhl,1,ghh,1)

C   ... Symmetrize over all stars, add each to gii
        else
          iq1 = 0
          reads = .false.

          iq1 = 0
          do  while (.true.)  ! loop over star of q
            call iqstar(1,iq,nk1,nk2,nk3,ipq,ifac,qb,iq1,q1,[0])
            if (iq1 == 0) exit

C          do  i3 = 1, nk3
C          do  i2 = 1, nk2
C          do  i1 = 1, nk1
C
C            iq1 = iq1+1
CC       ... skip this qp unless it is related to iq
C            if (ipq(iq1) /= iq) cycle

C       ... Make g by rotation of g(iq): symop relating iq1 to iq
            ig = igstar(iq1)

C           print "(6i4)", iq,iq1,ig

C       ... q into which h or z is rotated
C            q1(1) = qk(1,i1,i2,i3)
C            q1(2) = qk(2,i1,i2,i3)
C            q1(3) = qk(3,i1,i2,i3)

C       ... Rotate g making all spinors parallel to z
C           Should orbital parts be rotated too !?
            if (nspc==2 .and. ig>1 .and. .not. lgrotu) then
              call rotspn(30000+100*rotspd(0),1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
     .          s_ham%qss(4),xx,xx,lidim,ldim,ldim,lidim,lidim,xx,sll)
              lgrotu = .true.
C             call yprm('gii(fixed z)',2,sll,lix2,lix,lix,lix)
            endif

C       ... Rotate g from qp(iq) to q1, add to gii
C           call yprmi('gii iq=%i',iq1,0,2,sll,lix2,lix,lix,lix)

            if (ig /= 1) then

C         ... Rotate ll block; put in grot
              if (nfilet == 0) then
              if (associated(grot,sll)) allocate(grot(lix,lix))
C             debug = ig == 2
C              if (debug) call yprmi('before rot from %s,q=(%3;6d) iq,ig=%2:1i',qp(1,iq),[iq,ig],2,sll,lix2,lix,lix,lix)
              if (rothnc(optrot,nl,nspc,nbas,pos,xx,s_ham%offH,s_ham%iprmb,
     .          istab(1,ig),g(1,1,ig),ag(1,ig),q1,lidim,lidim,sll,grot) < 0) goto 999
C              if (debug) call yprmi('after rot to %s,q=(%3;6d) iq,ig=%2:1i',q1,[iq,ig],2,grot,lix2,lix,lix,lix)
              else
                if (nspc == 2) then
                  call rx('file save branch not ready for nspc=2')
                endif
C                rewind nfilet
C                call dpdump(gii,2*l2,-nfilet)
C                call dpdump(sll,2*l2,-nfilet)
CC               Flag that sll was overwritten and must be re-read
C                reads = .true.
CC                if (is+js /= 2) then
CC                  call ymscop(0,lidim,lidim,lix,lidim,(is-1)*lidim,
CC     .              (js-1)*lidim,0,0,sll,lix2,gii,l2)
CC                  call dcopy(2*l2,gii,1,sll,1)
CC                endif
C                if (roth(optrot,nl,nbas,pos,xx,s_ham%offH,s_ham%iprmb,
C     .            istab(1,ig),g(1,1,ig),ag(1,ig),q1,lidim,lidim,gii,
C     .            sll) < 0) goto 999
C                rewind nfilet
C                call dpdump(gii,2*l2,nfilet)
              endif

C#ifdefC DEBUG
CC         ... debugging check
C              print *, qp(1,iq),qp(2,iq),qp(3,iq)
C              print *, q1
C              call yprm('g(qp)',2,sll,lix2,lix,lix,lix)
C              print *, i1,i2,i3,ig,grot(lix+1,lix)
C              call yprm('g(q1)',2,grot,lix2,lix,lix,lix)
C#endif

C       ... Restore to local quantization axis
            if (lgrotu) then
C             call yprm('rotated gii(fixed z)',2,grot,lix2,lix,lix,lix)
              call rotspn(30000+100*rotspd(1),1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,
     .          s_ham%neula,s_ham%qss(4),xx,xx,lidim,ldim,ldim,lidim,
     .          lidim,xx,grot)
C             call snot(ldim*nspc,i1,i2,i3,ig,iq,grot)
C             call yprm('gii(local z)',2,grot,lix2,lix,lix,lix)
            endif

C      ... Accumulate rotated g, into gii
            call ymsadd(lix,lix,lix,lix,0,0,0,0,alpha,beta,grot,lix2,gii,lix2)
            if (reads) call dpdump(sll,2*l2,nfilet)
            reads = .false.

C       ... Rotate hh block
            if (hdim /= 0) then
              allocate(ghhl1(hdim,mxorb),ghhl2(hdim,mxorb))
              do  k = 1, nspc  ! diagonal spin blocks only
              call dpzero(ghhl1,2*hdim*mxorb)
              call ymscop(0,hdim,mxorb,hdim*nspc,hdim,(k-1)*hdim,
     .          0,0,0,ghhl,hhx,ghhl1,hdim*mxorb)
              if (roth(optrth,nl,nbas,pos,1,s_ham%offH,s_ham%iprmb,
     .          istab(1,ig),g(1,1,ig),ag(1,ig),q1,hdim,mxorb,ghhl2,
     .          ghhl1) < 0) goto 999

C           ..Accumulate higher block (diagonal in spinors)
              call ymsadd(hdim,mxorb,hdim,hdim*nspc,0,0,(k-1)*hdim,
     .          0,alpha,beta,ghhl1,hdim*mxorb,ghh,hhx)
              deallocate(ghhl1,ghhl2)
              enddo
            endif

          else  ! if no rotation of g
C           call snot(ldim*nspc,i1,i2,i3,ig,iq,grot)
            call ymsadd(lix,lix,lix,lix,0,0,0,0,alpha,beta,grot,lix2,gii,lix2)
            call ymsadd(hdim*nspc,mxorb,hdim*nspc,hdim*nspc,0,0,0,0,
     .        alpha,beta,ghhl,hhx,ghh,hhx)
          endif ! of rotation of g

          enddo
C          enddo
C          enddo
          if (.not. associated(grot,sll)) deallocate(grot)


C         call yprm('after rot',2,gii,lix2,lix,lix,lix)
        endif

C#ifdefC DEBUG
C        print *, 'gfibz DEBUG done iq=',iq
C#endif
      enddo
      deallocate(kpproc)
      call mpibc2(gii,lix2,6,3,.false.,'gfibz','gii')
      if (numprocs > 1) call info0(50,0,0, ' ... Done MPI k-loop')
      if (hdim /= 0) call mpibc2(ghh,hhx,6,3,.false.,'gfibz','ghh')

C --- G calculated at every qp in the BZ ---
      else
        do  i3 = 1, nk3
        do  i2 = 1, nk2
        do  i1 = 1, nk1
C ...   This branch needs to be parallelized

          q1(1) = qk(1,i1,i2,i3)
          q1(2) = qk(2,i1,i2,i3)
          q1(3) = qk(3,i1,i2,i3)

          call gf1kp(s_site,s_ham,s_pot,mod1kp,lbloch,lrel,nl,nbas,isp,nsp,
     .      nspc,nlibu,lmaxu,vorb,idu,ludiag,lrsig,q1,plat,salp,iax,
     .      nttab,lidim,hdim,sll,ghhl)
C         call yprm('gf',2,sll,lidim*lidim,lidim,lidim,lidim))

C     ... Add to gii
          call daxpy(lix2*2,1d0/nkfbz,sll,1,gii,1)

        enddo
        enddo
        enddo
      endif

C#ifdefC DEBUG
C      call rx0( 'debugging check done')
C#endif

C --- Symmetrize the k-integrated Green's function gii
      if (s_ctrl%nccomp /= 0 .and. mod(mode1,2) /= 0 .and. lrel /= 2) then
        call symcp(optrot,gii,lidim,nspc,nbas,nl,s_ham%offH,s_ham%iprmb,
     .    nsgrp,istab,g,ag,nspc,kcplx)
C       call yprm('after sym',2,gii,lix2,lix,lix,lix)
      endif

      if (iprint() >= PRTG/1 .and. mod(mode1,2) /= 0) then
C        call yprm('sum_q gii',2,gii,lix2,lix,lix,lix)
C        call rotspn(30000+100*rotspd(0),1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,
C     .    s_ham%neula,s_ham%qss(4),xx,xx,lidim,lidim,lidim,lidim,lidim,xx,gii)
        call yprm('sum_q gii',2,gii,lix2,lix,lix,lix)
        if (hdim /= 0)
     .    call yprm('sum_q ghh',2,ghh,hhx,hdim*nspc,hdim*nspc,mxorb)
      endif

  999 continue
      if (allocated(sll)) deallocate(sll)
      if (allocated(ghhl)) deallocate(ghhl)
      call tcx('gfibz')
C      print *, 'exit gfibz sum gii = ',
C     .  sum(sngl(gii(1:lix2))), sum(sngl(gii(lix2+1:2*lix2)))

      end
C      subroutine snot(ldg,i1,i2,i3,ig,iq,g)
C      integer ldg,i1,i2,i3,ig,iq
C      double precision g(ldg,ldg)
C
C      print 345, i1,i2,i3,ig,iq,g(55,2)
C  345 format('g(55,2)',3i3,2i4,2f12.6)
C      end

C#ifdefC DEBUG
C      subroutine snot(opt,offH,indxsh,nbas,iq,iq1,ig,ldim,lidim,hdim,
C     .  nlmh,nspc,is,js,g1,g2,ghh,ghh2)
C      implicit none
C      integer n0H,nkap0
C      parameter (nkap0=4,n0H=5)
C      integer opt,nbas,ldim,lidim,iq,iq1,ig,is,js,nspc,offH(n0H,nkap0,1)
C      integer hdim,nlmh,indxsh(*)
C      double precision g1(lidim,lidim,2),g2(lidim,nspc,lidim,nspc,2)
C      double precision ghh(hdim,nlmh,2),ghh2(hdim,nspc,nlmh,2)
C      double precision tol,errmx,errmxl
C      integer i,j,lix,imax,jmax,i1,i2,j1,j2,opt0,off0,jb,iblb,n0,lhdim,
C     .  iprint,nkap0
C      parameter (n0=10,nkap0=4)
C      integer norb,nlmi,ixx,iorb,li,offj,
C     .  ltab(n0*nkap0),ktab(n0*nkap0),offlj(n0*nkap0)
C      double complex wk
C
C      opt0 = mod(opt,10)
C
CC      print 321, opt,iq,iq1,is,js
CC  321 format('check opt,iq,iq1,is,js',i6,5i4)
C      errmx = 0
C      errmxl = 0
C      i1 = 1
C      i2 = lidim
C      j1 = 1
C      j2 = lidim
C      do  10  j  = j1, j2
C
C      if (opt0 == 1) then
CC       Find basis atom corresponding to site
C        iblb = 1
C        off0 = 0
C        if (j > ldim) iblb = 2
C        if (j > ldim) off0 = ldim
C        do  12  jb = 1, nbas
C          if (j <= offH(iblb,1,jb)+off0) goto 12
C          if (j > offH(iblb,1,jb+1)+off0) goto 12
C          i1 = offH(iblb,1,jb)+1 + off0
C          i2 = offH(iblb,1,jb+1) + off0
C          goto 14
C   12   continue
C        goto 10
C   14   continue
C      endif
C
C      if (j == j1 .and. iprint() >= 60) print 357,'li',0,j1,j2,i1,i2
C  357 format(' check ',a,' jb=',i3,'  j1,j2=',2i4,'  i1,i2=',2i4)
C      do  11  i  = i1, i2
C        wk = dcmplx(g1(i,j,1)-g2(i,is,j,js,1),g1(i,j,2)-g2(i,is,j,js,2))
C        if (cdabs(wk) > errmx) then
C          imax = i
C          jmax = j
C          if (i <= ldim .and. j <= ldim) errmxl = cdabs(wk)
C          errmx = cdabs(wk)
C        endif
C   11 continue
C   10 continue
C
C      tol = 1d-7
CC     tol = 2d-5
C      if (errmxl > tol) then
C        print 334, iq,iq1,ig,imax,jmax,errmx,errmxl
C  334   format('mismatch iq=',2i5,' ig=',i2,' i,j(max)=',2i3,
C     .    '  max err=',1pe9.2,' in l block',1pe9.2)
CC        call yprm('g1',2,g1,lidim*lidim,lidim,lidim,lidim)
CC        lix = lidim*nspc
CC        call yprm('g2',2,g2(1,is,1,js,1),lix**2,lidim,lidim,lidim)
C      endif
C
C      if (hdim == 0 .or. is /= js) return
C
C      errmx = 0
C      errmxl = 0
C
C      lhdim = lidim+hdim
C      do  110  jb = 1, nbas
C        call orbl(jb,lidim,lhdim,indxsh,norb,ltab,ktab,ixx,offlj,nlmi)
C        i2 = 0
C        do  112  iorb = 1, norb
C          li = ltab(iorb)
C          nlmi = 2*li + 1
C          offj = offlj(iorb) - lidim
C          j1 = offj+1
C          j2 = offj+nlmi
C          i1 = i2+1
C          i2 = i2+nlmi
C          if (iprint() >= 60) print 357, 'hh',jb,j1,j2,i1,i2
C          do 120  j = j1, j2
C          do 120  i = i1, i2
C            wk = dcmplx(ghh(j,i,1)-ghh2(j,is,i,1),
C     .                  ghh(j,i,2)-ghh2(j,is,i,2))
C            if (cdabs(wk) > errmx) then
C              imax = i
C              jmax = j
C              errmx = cdabs(wk)
C            endif
C  120     continue
C  112   continue
C  110 continue
C
C      tol = 1d-7
CC     tol = 2d-5
C      if (errmx > tol) then
C        print 335, iq,iq1,ig,jmax,imax,errmx
C  335   format('mismatch iq=',2i5,' ig=',i2,' i,j(max)=',2i3,
C     .    '  max err in hh block',1pe9.2)
CC        call yprm('g1',2,g1,lidim*lidim,lidim,lidim,lidim)
CC        lix = lidim*nspc
CC        call yprm('g2',2,g2(1,is,1,js,1),lix**2,lidim,lidim,lidim)
C      endif
C
C      end
C
C#endif
