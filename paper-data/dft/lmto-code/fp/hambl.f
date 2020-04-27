      subroutine hambl(mode,s_site,s_spec,s_lat,s_ham,isp,q,
     .  k1,k2,k3,smpot,vconst,lcplxp,alfa,ndimh,napw,igapw,h,s,hso)
C- Make the LDA hamiltonian and overlap for one k-point
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  augmbl bstrux smhsbl hsibq hsubblock hsibq2 hsibq4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  augmbl bstrux uspecb smhsbl hsibq tbhsi hsubblock
Cio                hsibq2 hsibq4
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat qlat vol plat awald tol nkd nkq gmax nabc ng kv
Ci                 kv2 igv igv2
Co     Stored:     igv igv2
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:qlat plat cg indxcg jcg cy qlv dlv gv igv kv igv2
Cio    Passed to:  augmbl bstrux hxpbl ghibl hklbl gklbl hxpgbl ghigbl
Cio                hklgbl smhsbl hhibl phhibl hsmbl hsibq sugvec
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  pwmode
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit specifies LDA matrix elements
Ci         :  0 compute LDA hamiltonian and overlap
Ci         :  1 compute the overlap only
Ci         :10s digit
Ci         :  0 do not compute hso
Ci         :  1 compute hs0
Ci         :    Note: only a portion of hso is computed for a
Ci         :    particular isp.  The total hso is assembled
Ci         :    after isp loops from 1..2.  hso should not be
Ci         :    initialized between isp=1 and isp=2 loops.
Ci         :100s digit
Ci         :  0 augmentation matrices from sig,tau,ppi
Ci         :  1 augmentation matrices from sigx,taux,ppix
Ci   isp   :spin index
Ci   q     :Bloch vector (k-point)
Ci   k1,k2,k3 dimensions of smpot
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst:additional constant potential
Ci   lcplxp:0 if site potential s_site%ppi is real; 1 if ppi is complex
Ci   alfa  :add alfa * overlap to hamiltonian
Ci         :This is for stability in evals.  Preferably alfa=0
Ci   ndimh :dimension of hamiltonian and ovrlap h,s
Ci   napw  :number of APW's
Ci   igapw :APWs in units of reciprocal lattice vectors
Co Outputs
Co   h     :Hamiltonian matrix
Co   s     :overlap matrix
Co   hso   :spin off-diagonal block of spin-orbit hamiltonian
Cr Remarks
Cu Updates
Cu   08 Feb 13 Internally shortens q vector
Cu   10 Nov 11 Begin migration to f90 structures
Cu   04 Jul 08 (T. Kotani) New PW addition to basis
Cu   03 Feb 05 (A. Chantis) calculate hso
Cu    1 Sep 04 Adapted to handle complex ppi; S.O. put into ppi
Cu   25 Aug 04 modifications for extended local orbitals
Cu   15 Jul 04 (Chantis) Add Lz.Sz spin-orbit coupling
Cu   10 Jan 03 Remove addition from hambl.  See hambls.f
Cu   10 Jan 03 put sigma back into Bloch transform only
Cu   14 Aug 02 Added overlap-only option and option for orthog sigm
Cu   20 Jul 02 Can add Bloch transform of sigma matrix to ham
Cu   18 May 00 Adapted from nfp mk_hamiltonian.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isp,k1,k2,k3,ndimh,lcplxp,napw
      integer,target :: igapw(3,napw)
      double precision alfa
      double precision q(3),vconst
      double complex smpot(k1,k2,k3),h(ndimh,ndimh),s(ndimh,ndimh)
      double complex hso(ndimh,ndimh)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Local parameters
      character outs*256
      integer mode0,mode1,mode2,ipr
      double precision vavg,qs(3),tol
      parameter (tol=1d-8)
      complex(8), allocatable :: zwk(:,:,:)
      real(8), allocatable :: e0(:)
      integer,pointer :: igapwl(:,:)
      procedure(real(8)) :: dlength
      procedure(logical) :: cmdopt

      call tcn('hambl')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)

C ... Shorten q; shift APW G vectors to correspond
      igapwl => igapw
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        allocate(igapwl(3,napw))
        call shorigv(q,qs,s_lat%plat,napw,igapw,igapwl)
      endif

C ... Initialization of quantities to be computed
      if (mode0 == 0) call dpzero(h,ndimh*ndimh*2)
      if (mode0 <= 1) call dpzero(s,ndimh*ndimh*2)
      if (mode1 == 1 .and. isp == 1) call dpzero(hso,ndimh*ndimh*2)
C ... Skip LDA H, S
      if (mode0 > 1) goto 100

C --- Augmentation part of hamiltonian and overlap ---
      call augmbl(mode0+10*mode1+100*mode2,s_site,s_spec,s_lat,s_ham,
     .  isp,lcplxp,qs,0,ndimh,napw,igapwl,h,hso,s)
C      ipr = 0
C      if (ipr>0) call zprm('h after aug',2,h,ndimh,ndimh,ndimh)
C      if (ipr>0) call zprm('s after aug',2,s,ndimh,ndimh,ndimh)
C      if (ipr>0) h = 0
C      if (ipr>0) goto 200

C ... Optionally set average istl potential to zero
      vavg = 0
C     call pvhmb1(0,s_lat,k1,k2,k3,smpot,vavg)
C     call pvhmb1(1,s_lat,k1,k2,k3,smpot,-vavg)
C --- Smoothed hamiltonian and overlap, S and T done analytically ---
C     print *, '!!'; s = 0
      call smhsbl(mode0,s_site,s_spec,s_lat,vavg+vconst,qs,ndimh,
     .  s_ham%iprmb,napw,igapwl,h,s)
C     call zprm('s',12,s,ndimh,ndimh,ndimh)

      if (mode0 == 1) goto 100
C#ifdefC MPI
C      call hsibl(s_site,s_spec,s_lat,k1,k2,k3,smpot,
C     .  isp,qs,ndimh,s_ham%iprmb,napw,igapwl,h)
C#else
      call hsibq(s_site,s_spec,s_lat,k1,k2,k3,smpot,isp,qs,ndimh,
     .  s_ham%iprmb,napw,igapwl,h)
C#endif

C ... Restore istl potential
C      call pvhmb1(1,s_lat,k1,k2,k3,smpot,vavg)

C --- Alternatively, smoothed H and S done completely numerically ---
C      ngabc = s_lat%nabc
C      call hsmi_q(n1,n2,n3,k1,k2,k3,smpot,
C     .  vconst,qs,ndimh,h,s)

C ... End of assembling LDA H, S
  100 continue

C --- Fill second half of matrix ---
C 200 continue
      if (mode0 == 0) call z2herm('U',ndimh,ndimh,h)
      if (mode0 <= 1) call z2herm('U',ndimh,ndimh,s)

C --- Add alfa I to s and  alfa (s^-1 h + cc)/2 to h ---
CC... Add ds (sx+ds)^-1 hx or ds (sx)^-1 hx  ... latter is better for high evals but must invert negative matrix
CCmc -f4f18.14 -valfa=1e-4 -1:4 -salfa -w ds -a ds ds sx ds -s1 -+ -i -x -herm  hx -x -herm > deltah
CCmc -f4f18.10 sx ds  -+ hx deltah -+ -gevl  sx hx -gevl -ccat -e3 x1 x2 x1-x2
C      if (dabs(alfa) > 1d-10 .and. mode0 == 0)  then
CC        call zprm('s',2,zwk,ndimh,ndimh,ndimh)
CC        call zprm('h',2,zwk(1,1,2),ndimh,ndimh,ndimh)
C
CC Add alpha I to s and alpha (s + alpha I)^-1 h as perturbation ... problems
CC        forall (ipr = 1:ndimh) s(ipr,ipr) = s(ipr,ipr) + alfa
CC        allocate(zwk(ndimh,ndimh+1,2))
CC        call dcopy(2*ndimh**2,s,1,zwk,1)
CC        call dcopy(2*ndimh**2,h,1,zwk(1,1,2),1)
CC        call zpotrf('U',ndimh,zwk,ndimh,ipr) ! Cholesky decomposition
CC        if (ipr /= 0) call rxi('hambl: zpotrf failed to decompose overlap info=',ipr)
CC        call zpotrs('U',ndimh,ndimh,zwk,ndimh,zwk(1,1,2),ndimh,ipr)
CC Add alpha I to s and alpha s^-1 h as perturbation
C        allocate(zwk(ndimh,ndimh+1,3))
C        call dcopy(2*ndimh**2,s,1,zwk,1)
C        call dcopy(2*ndimh**2,h,1,zwk(1,1,2),1)
C        call zqinvb('lh',zwk,ndimh,ndimh,ndimh,zwk(1,1,3),ndimh,zwk(1,1,3),zwk(1,1,2),ndimh,ipr)
C        forall (ipr = 1:ndimh) s(ipr,ipr) = s(ipr,ipr) + alfa
CC       call zprm('s^-1 h',2,zwk(1,1,2),ndimh,ndimh,ndimh)
C        call z2herm('A',ndimh,ndimh,zwk(1,1,2))
C        call daxpy(2*ndimh**2,alfa,zwk(1,1,2),1,h,1) ! h += alfa [S^-1 h + cc]/2
CC       call zprm('(alfa s^-1 h + cc)/2',2,h,ndimh,ndimh,ndimh)
C        deallocate(zwk)
C      endif
C
CC      This does not work
CC      if (dabs(alfa) > 1d-10 .and. mode0 == 0)  then
CC        allocate(e0(ndimh))
CC        call zprm('unmodified s',12,s,ndimh,ndimh,ndimh)
CC        call zhevo(ndimh,ndimh,s,s,0,0d0,alfa,0,1,0,e0,e0,ndimh,s)
CC        call zprm('modified s',2,12,ndimh,ndimh,ndimh)
CC        deallocate(e0)
CC      endif

C --- Zero or modify subblock of h, s ---
      if (cmdopt('--zhblock',9,0,outs)) then
        call zhblock(outs(10:),ndimh,isp,h,s)
        call pshpr(1)
        if (mode1 == 1) call zhblock(outs(10:),ndimh,110,h,s)
        call poppr
      endif

C --- Make and display smallest eigenvalues of s ---
      if (cmdopt('--show-ovl-evl',14,0,outs)) then
        allocate(zwk(ndimh,ndimh,2),e0(ndimh))
        call dcopy(2*ndimh*ndimh,s,1,zwk,1)
        call zhevx(ndimh,ndimh,zwk,zwk,0,.true.,ndimh*0,9d99,ipr,zwk(1,1,2),.false.,e0,ndimh,zwk)
        call info5(10,1,1,' hambl: smallest evals of ovlp for q=%3:1;6,6d'//
     .    '%N sevl%n:1;4,4g',qs,min(ndimh,8),e0,4,5)
        deallocate(zwk,e0)
      endif

      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        deallocate(igapwl)
      endif

      call getpr(ipr)
      if (ipr >= 90) then
        call info(0,0,0,' h and s for q=%3:1;6,6d',qs,0)
        if (mode0 == 0) call zprm('h',12,h,ndimh,ndimh,ndimh)
        if (mode0 <= 1) call zprm('s',12,s,ndimh,ndimh,ndimh)
        if (mode1 == 1) call zprm('hso',2,hso,ndimh,ndimh,ndimh)
      endif
      call tcx('hambl')
      end

      subroutine pvhmb1(mode,s_lat,k1,k2,k3,smpot,vavg)
C- Add subtract a constant from smooth potential
C ----------------------------------------------------------------------
Ci Inputs
Ci    mode :0 do not adjust smpot; just return avg smpot
Ci         :1 add vavg to smpot
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc
Ci   k1,k2,k3 dimensions of smpot
Cio Inputs/Outputs
Cio  smpot :(input for mode=0, altered for mode=1)
Cio        :smooth potential on uniform mesh (mkpot.f)
Cio  vavg  :(output for mode=0, input for mode=1) average potential
Cr Remarks
Cr
Cu Updates
Cu   16 Aug 04 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,k1,k2,k3
      double complex smpot(k1,k2,k3)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Local parameters
      double precision vavg
      double complex xx
      integer ngabc(3),i

      if (mode == 0) then
        ngabc = s_lat%nabc
        call mshint(1d0,1,ngabc(1),ngabc(2),ngabc(3),
     .    k1,k2,k3,smpot,vavg,xx)
      elseif (mode == 1) then
        do i = 1, k1*k2*k3
          smpot(i,1,1) = smpot(i,1,1) + vavg
        enddo
      else
        call rx('pvhmb1: bad mode')
      endif
      end

C      subroutine pvhmb2(ndimh,h,s)
CC- test routine to doctor hamiltonian
C      implicit none
C      integer ndimh
C      double complex h(ndimh,ndimh),s(ndimh,ndimh)
C
C      integer low,high,i1,i2,k1,k2
C      return
C
CC ... Zeros out off-diagonal blocks between low+1, high
CC      low = 20
CC      high = 23
CC      print *, 'hambl zero out off-diagonal blocks', low+1,high
CC      do  i1 = 1, ndimh
CC        k1 = -1
CC        if (i1 > low) k1 = 0
CC        if (i1 > high) k1 = 1
CC        do  i2 = 1, ndimh
CC          k2 = -1
CC          if (i2 > low) k2 = 0
CC          if (i2 > high) k2 = 1
CC
CC          if (k1 == 0  .neqv.  k2 == 0) then
CC            h(i1,i2) = 0
CC            s(i1,i2) = 0
CC          endif
CC        enddo
CC      enddo
C
C      end
