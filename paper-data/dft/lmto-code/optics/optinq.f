      subroutine optinq(s_ctrl,s_optic,s_bz,eband,nbmax,nsp,nspc,efermi,
     .  ipq,alat,plat,nfilm,nempm,nptmx,npol,njdosw,optmt,jdosw,opt)
C- BZ integration of Im(eps) or joint DOS by tetrahedron method
C ----------------------------------------------------------------------
Cio Structures
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic lpart dw window ocrng unrng esciss iq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp lmet lshft ntet n range w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Ci Inputs:
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   efermi:Fermi level
Ci         :NB: efermi=-99 is a flag, indicating that the Fermi
Ci         :level is not set.
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nfilm,nempm: dimensions optmt.  Should be nfiup-nfilo+1,nemup-nemlo+1
Ci   nptmx :leading dimension of opt
Ci   npol  :number of polarizations
Ci         :Case njdosw=0 or loptic>0: npol should be 3
Ci         :Case njdosw>0 and loptic<0: npol should be 1+njdosw
Ci   optmt :<i|grad|j> connecting occ i with unocc j
Ci         :optmt is not used when making joint DOS (loptc=-1)
Ci         :only states i between nfilo and nfiup, and states j
Ci         :between nemlo and nemup are supplied
Co Outputs:
Co   popt  :opt, resolved in one of the following ways.
Co   opt   :Im(epsilon), unscaled by weight 16*pi**2/vol/e**2
Co         :opt is either :
Co         :(1) generated and the scaled epsilon written to file 'opt'
Co          (when information for all k-points is available; see ikp);
Co         :(2) partial contribution added to the unscaled epsilon
Co          (when information for 1 k-points is available; see ikp);
Co         :NB: output may be joint DOS instead of Im eps; see jdos
Co         :    below
Co         :NB: opt is NOT initialized by this routine if accumulating
Co         :    only 1 kp (ikp>0)
Co         :Caller should set opt=0 before first call to optinq.
Cl Local variables and switches
Cl  loptic  What optical property to calculate
Cl          Should be same as s_ctrl%loptc
Cl          1 generate Im(eps) (spin 1 only, if partial decompsn)
Cl          2 generate Im(eps) spin 2 (needed for partial decompsn)
Cl         -1 generate joint DOS
Cl         -2: generate joint DOS, spin 2 (use with LTET=3)
Cl         -3: generate up-down joint DOS (LTET=3)
Cl         -4: generate down-up joint DOS (LTET=3)
Cl         -5: generate spin-up single DOS (LTET=3)
Cl         -6: generate spin-dn single DOS (LTET=3)
Cl         Add 20 for Second Harmonic Generation
Cl   lpart :Switch controlling how to resolve JDOS or Im eps into parts
Cl         :1s digit
Cl         :0  No decomposition
Cl         :1  Resolve eps or dos into (occ,unocc) parts
Cl         :2  Resolve eps or dos by k
Cl         :3  Both 1 and 2
Cl         :10s digit
Cl         :0  save as ASCII file
Cl         :1  save as binary file
Cl   nfilo :first occupied band to include in Im eps or JDOS
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nfiup :last  occupied band to include in Im eps or JDOS
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nemlo :first unoccupied band
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   nemup :last  unoccupied band
Cl         :Also optmt = optmt(:,nfilm,nempm,nsp,iq)
Cl   liq   :True if to compute JDOS or Im eps with unocc shifted by q
Cb Bugs
Cb   Set of nfilo,nfiup,nemlo,nemup has duplicate nfilox,nfiupx,nemlox,nemupx
Cb   2nd set is used in tetwtq, but not clear why it should be distinct.
Cb   Possibly affects tetrahedron weights?
Cr Remarks
Cr   Uses tetwtq to makes joint density of states or Im(eps)
Cr   tetwtq calculates a (occ,unocc,k) resolved set of weights
Cr   for sequence of histograms 1,2,3,...
Cr     whw(i,ib,jb,k) =
Cr        \int_om_i^om_i+1 d\omega \times
Cr        \int d^3k f(e(k))(1-f(e(q+k)))\delta(omg - e(q+k) + e(k))
Cr
Cr   tetwtt uses these weights to assemble Im(eps) or JDOS.
Cr
Cu Updates
Cu   27 Aug 14 Better printout for memory management
Cu   21 Jun 14 s_optic%loptic=1 handles both spins unless partial decompsn
Cu   06 Jan 14 Implemented s_optic%loptic=2
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   09 Aug 12 jdosw is spin polarized
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Sep 09 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbmax,nsp,nspc,nfilm,nempm,nptmx,ipq(*),npol,njdosw
      double precision efermi,alat,plat(3,3)
      double precision opt(nptmx,npol,nsp),jdosw(max(njdosw,1),nfilm+nempm,*)
      double precision optmt(3,nfilm,nempm,nsp,*),
     .                 eband(nbmax,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical jdos,metal,ltop,liq,ltmp,memcheck
      integer i,ilo,ipr,iq1,iup,j,jpm,k,lgunit,isp,jsp,ipol,lmet,loptic,
     .  lpart,mxbnds,nband,nemlo,nemup,nfilo,nfiup,nkp,npts,ntet,nwhis,
     .  stdo,nfilmx,nempmx,nkpw,idum,nfilmw,nempw,nfilox,nfiupx,nemlox,
     .  nemupx,ibw,jbw,lpart0,lpart1,ksp,npass
      integer idalloc,allocvb
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3,lshft(3),i1,i2,i3,nqbz,nqbzw,iqshft(3)
      integer npm,nctot,ncc,nwgtx,jhwtot,nhwtot,iq,jj1,jj2,jj3,ifac(3)
      double precision emin,emax,pi,esciss,vol,de(3),saveef,fac,
     .  optrng(4),q(3),qlat(3,3),tiny,qb(3,3),xv(10),qoff(3)
      real demaxx
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      real(8),allocatable :: qbzw(:,:),ecore(:,:)
      real(8),allocatable,target :: efbz1(:,:,:,:),efbz2(:,:,:,:)
      real(8),pointer:: eocc(:,:,:,:),eunocc(:,:,:,:)
      real(8),allocatable :: emesh(:),frhis(:),popt(:,:,:,:,:)
      real(8),allocatable :: wtthis(:,:),wtphis(:,:,:,:,:),whw(:)
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      integer(4),allocatable::
     .  ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      integer,allocatable :: ib1bz(:),idtet(:,:),nwgt(:,:),nwgttt(:,:),
     .  nib(:,:,:),njb(:,:,:)
      integer :: noccxvx(2)=-9999
      logical,allocatable :: iwgt(:,:,:,:)
      character strn*32
      logical, parameter :: T=.true., F=.false.
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      double precision qk
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)
C ... External calls
      external cspline,csplint,dinv33,dosmsh,dpzero,dscal,getpr,
     .         hisrange,info2,info5,info8,optin2,optind,pair_index,
     .         pralloc,rxi,rxx,tcn,tcx,tetfbz,tetwtq,tetwtt

C ... Setup
      call tcn('optinq')
      call getpr(ipr)
      stdo = lgunit(1)
      loptic = s_optic%loptic
      if (loptic == 0) return
      lpart = s_optic%lpart
      de = s_optic%dw
      optrng = s_optic%window
      ocrng = s_optic%ocrng(1:2)
      unrng = s_optic%unrng(1:2)
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      lmet = s_bz%lmet
      lshft = s_bz%lshft
      esciss = s_optic%esciss
      iqshft = s_optic%iq
      jdos = loptic < 0
      liq = iqshft(1)+iqshft(2)+iqshft(3) /= 0
      ifac = 1
      metal = lmet /= 0
      call rxx(metal .and. dabs(esciss) > 1d-8,
     .  'OPTINQ: metal not compatible with scissors operator')
      pi = 4d0*datan(1d0)
      if (loptic >= -4 .and. loptic <= -2 .and. nsp == 1) call
     .  rx1('OPTINQ: loptic=%i but calculation is not spin pol',loptic)
      if (nfiup < nfilo) call rx('OPTINQ: FILBND has negative range')
      if (nemup < nemlo) call rx('OPTINQ: EMPBND has negative range')
C     A number around 4*10^-6 that has just one bit.
      tiny = 1d0/2**18
C ... Mesh setup ... get npts for printout
      i = 1; if (de(2) /= 0) i = 11
      call dosmsh(-i,optrng,de,nptmx,npts,xv)
      saveef = efermi
      lpart0 = mod(lpart,10)
      lpart1 = mod(lpart/10,10)

C       i = njdosw*(nfilm+nempm)
CC       print *, '!! reproduce bug'
CC      call dcopy(i*nkp,jdosw(1,1,1+nkp*(2-1)),1,jdosw,1)
C       call prmx('jdosw',jdosw,i,i,nkp*nsp)

C      Testing
C      print *, '!!'
C      efermi = 0
C      eband(1,:,1:nkp) = efermi - .1d0
C      eband(2,:,1:nkp) = efermi + .1d0
C      eband(3:,:,1:nkp) = efermi + 999
C      nfilo = 1; nfiup = 1
C      nemlo = 2; nemup = 2

      if (mod(loptic,4) > 0 .and. npol /= 3) then
        call rxi('optinq: illegal value, npol =',npol)
      elseif (njdosw > 0 .and. npol /= njdosw+1) then
        call rxi('optinq: illegal value, npol =',npol)
      endif
      call dpzero(opt,nptmx*npol*nsp)

      if (jdos .and. (loptic == -5 .or. loptic == -6)) then
        strn='DOS'; if (loptic == -6) strn='DOS (spin 2)'
        call info8(20,1,0,
     .    ' OPTINQ:  '//trim(strn)//' from occ bands (%i,%i)'//
     .    '%?#(n)#, %-1j%i Mulliken projection(s)##'//
     .    '%?#(n==1)#%N%10fResolve into band contributions##%-1j'//
     .    '%?#(n==2)#%N%10fResolve by k##%-1j'//
     .    '%?#(n==3)#%N%10fResolve into band contributions and by k##'//
     .    '%N%10f%i points in energy window (%,1d,%,1d)  Efermi=%;6d',
     .    nfilo,nfiup,njdosw,lpart0,npts,emin,emax,efermi)
      elseif (jdos) then
        call info5(20,1,0,
     .    ' OPTINQ:  JDOS '//
     .    '%?#(n==-1)#(++) ##%-1j'//
     .    '%?#(n==-2)#(--) ##%-1j'//
     .    '%?#(n==-3)#(+-) ##%-1j'//
     .    '%?#(n==-4)#(-+) ##'//
     .    'using occ bands (%i,%i) and unocc (%i,%i)',
     .    loptic,nfilo,nfiup,nemlo,nemup)
        call info5(20,0,0,'%10f%i points in energy window (%d,%d)'//
     .    '  Efermi=%;6d',npts,emin,emax,efermi,0)
        if (lpart0 /= 0) then
          call info2(20,0,0,
     .    '%?#(n==1)#%10fResolve into (occ,unocc) band pair '//
     .      'contributions##%-1j'//
     .    '%?#(n==2)#%10fResolve by k##%-1j'//
     .    '%?#(n==3)#%10fResolve into band contributions and by k##'//
     .    '%?#(n)#, %-1j%i Mulliken projection(s)##',
     .    lpart0,njdosw)
        endif
      else
        strn='Im eps(w)'; if (loptic == 2) strn='Im eps(w) (spin 2)'
        if (loptic == 1 .and. nsp == 2) strn='Im eps(w) (spin 1)'
        call info8(20,1,0,
     .    ' OPTINQ:  '//trim(strn)//' using occ bands (%i,%i) '//
     .    'and unocc (%i,%i)%N%10f%i points in energy window (%d,%d)'//
     .    '  Efermi=%;6d',
     .    nfilo,nfiup,nemlo,nemup,npts,emin,emax,efermi)
        if (lpart0 /= 0) then
          call info2(20,0,0,
     .    '%?#(n==1)#%10fResolve into (occ,unocc) band pair '//
     .      'contributions##%-1j'//
     .    '%?#(n==2)#%10fResolve by k##%-1j'//
     .    '%?#(n==3)#%10fResolve into band contributions and by k##',
     .    lpart0,0)
        endif
      endif

C ... Tetrahedron setup for tetwtq
      call dinv33(plat,1,qlat,vol)
      vol = dabs(vol)
      nqbz  = n1*n2*n3
      nqbzw = (n1+1)*(n2+1)*(n3+1)
      ntet = 6*nqbz
      allocate(ib1bz(nqbzw))
      allocate(qbzw(3,nqbzw))
      allocate(idtet(4,ntet))
C      call tetfbz(qlat,n1,n2,n3,lshft,idtet,qbzw,ib1bz)
C      print *, 'cksum',sum(idtet),sum(qbzw),sum(ib1bz)
      qoff(1:3) = dble(lshft(1:3))/2
      call tetfbz(qlat,n1,n2,n3,0,.false.,1,qoff,1d0,idtet,qbzw,ib1bz,ntet)
      qb(1:3,1) = qlat(1:3,1)/n1
      qb(1:3,2) = qlat(1:3,2)/n2
      qb(1:3,3) = qlat(1:3,3)/n3
      q(1) = qk(1,iqshft(1)+1,iqshft(2)+1,iqshft(3)+1)
      q(2) = qk(2,iqshft(1)+1,iqshft(2)+1,iqshft(3)+1)
      q(3) = qk(3,iqshft(1)+1,iqshft(2)+1,iqshft(3)+1)
      if (liq) then
        call info2(20,0,0,'%10fq for unocc states: '//
     .    'iq =%3:1i, q =%3:1,6;6d',iqshft,q)
      endif

C ... Setup for to determine range of energy bands to copy to full BZ
C     Note: tetwtq has different nfilo,nfiup than popt!
C     For tetwtq, tetwtt, use these four: nfilox,nfiupx,nemlox,nemupx
      if (loptic /= -5 .and. loptic /= -6) then
C       Maximum number of bands with energies < efermi + emax
        i1 = mxbnds(0,eband,nbmax,nbmax,nkp*nsp,efermi+emax,i)
C       Fewest number of bands with energies < efermi - emax
        i2 = mxbnds(1,eband,nbmax,nbmax,nkp*nsp,efermi-emax,i)
C       Reduce number of states to size of largest window that contributes
        nfiupx = min(nfiup,i1)
        nemupx = min(nemup,i1)
        nfilox = max(nfilo,i2)
        nemlox = max(nemlo,i2)
        if (nfilox > nfiupx) call rx('bug in optinq')
        if (nemlox > nemupx) call rx('bug in optinq')
C ... Single DOS case: efermi at emax
      else
        if (lpart0 == 1)
     .    call rxi('OPTINQ: lpart0=1 not ready with loptic=',loptic)
        efermi = emax + tiny
        isp = 1
        if (loptic == -6 .or. loptic == -2 .or. loptic == 2) isp = 2
        i1 = 0; i2 = nbmax
C       At most i1 bands below emax; emin starts at i2 and above
        do  k  = 1, nkp
          i1 = max(mxbnds(0,eband(1,isp,k),nbmax,nbmax,1,emax,i),i1)
C          i2 = mxbnds(1,eband(1,isp,k),nbmax,nbmax,1,emin,i)
C          print *, i2,eband(i2,isp,k),eband(i2+1,isp,k)
          i2 = min(mxbnds(1,eband(1,isp,k),nbmax,nbmax,1,emin,i),i2)
        enddo
C       Use lowest band for i2 (tetrahedra may use band)
        i2 = 1

C       Need all bands up to emax
        nfiupx = min(nfiup,i1)
C       Need no bands below emin
        nfilox = max(nfilo,i2)
        if (nfilox > nfiupx) call rx('bug in optinq')
C       Put one (constant) band artificially above highest band
C       no ... just copy occ
        nemlox = nfiupx
        nemupx = nemlox
      endif

C ... Make ilo,iup,nfilox,nfiupx,nemlox,nemupx to fit tetwq style
C     tetwtq requires same number of bands for eocc,eunocc
C     ilo,iup dimension nband
      ilo = min(nfilox,nemlox)
      iup = max(nfiupx,nemupx)
      nband = iup - ilo + 1
      nfilox = ilo
      nfiupx = iup
      nemlox = ilo
      nemupx = iup

      if (loptic == -5 .or. loptic == -6) then
        call info2(20,0,0,'%10fReduce bands to range (%i,%i)',ilo,iup-1)
      else
        call info2(20,0,0,'%10fReduce bands to range (%i,%i)',ilo,iup)
      endif

      if (loptic == 3) loptic = 1
      ksp = 1 ! spin channel where to write opt and popt
      npass = 1 ! Loop over only one spin channel
      if (loptic == 1 .and. lpart0 == 0 .and. nsp > nspc) npass = 2

C --- Re-entry point for 2nd spin s_optic%loptic == 3
   10 continue

C ... Copy subblock of energy bands to efbz1: use full BZ
C     call pralloc('OPTINQ','bands',nband*dble(n1*n2*n3),100,40,1)
      allocate(efbz1(nband,n1,n2,n3))
      eocc => efbz1
      eunocc => efbz1
      if (liq) then
        allocate(efbz2(nband,n1,n2,n3))
        eunocc => efbz2
      endif
      if (iabs(loptic) <= 2) then
        isp = iabs(loptic)
      elseif (loptic == -3) then
        if (.not. allocated(efbz2)) allocate(efbz2(nband,n1,n2,n3))
        isp = 1
        jsp = 2
        eunocc => efbz2
      elseif (loptic == -4) then
        if (.not. allocated(efbz2)) allocate(efbz2(nband,n1,n2,n3))
        isp = 2
        jsp = 1
        eocc => efbz2
        eunocc => efbz1
      elseif (loptic == -5 .or. loptic == -6) then
        if (liq)
     .    call rxi('OPTICS_IQ nonsensical with loptic=',loptic)
        allocate(efbz2(nband,n1,n2,n3))
        eunocc => efbz2
        isp = 1
        if (loptic == -6) isp = 2
      endif

      do  iq = 1, nkp
        iq1 = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          iq1 = iq1+1
C     ... Copy bands if this qp is symmetry-equivalent to iq
          if (ipq(iq1) == iq) then
            jj1 = i1; jj2 = i2; jj3 = i3
            if (liq) then
              jj1 = mod(i1+iqshft(1)-1,n1) + 1
              jj2 = mod(i2+iqshft(2)-1,n2) + 1
              jj3 = mod(i3+iqshft(3)-1,n3) + 1
            endif
C           print "(3i4,2x,3i4)", i1,i2,i3,jj1,jj2,jj3
            ltop = .false.
            do  i = ilo, iup
              if (iabs(loptic) <= 2) then
                eocc(i-ilo+1,i1,i2,i3) = eband(i,isp,iq)
                if (liq) then
                eunocc(i-ilo+1,jj1,jj2,jj3) = eband(i,isp,iq)
                endif
              elseif (loptic == -3 .or. loptic == -4) then
                eocc  (i-ilo+1,i1,i2,i3) = eband(i,isp,iq)
                eunocc(i-ilo+1,jj1,jj2,jj3) = eband(i,jsp,iq)
              elseif (loptic == -5 .or. loptic == -6) then
                eocc(i-ilo+1,i1,i2,i3) = eband(i,isp,iq)
                eunocc(i-ilo+1,i1,i2,i3) = emax + (emax-emin) + tiny
                if (i == ilo) eunocc(1,i1,i2,i3) = efermi+tiny
              endif
            enddo
          endif
        enddo
        enddo
        enddo
      enddo

C ... Set up dimensions to wtphis, for resolved DOS or eps
      if (lpart0 == 0) then
        nfilmx = 0
        nempmx = 0
      elseif (lpart0 == 1 .or. lpart0 == 3) then
        nfilmx = nfiupx-nfilox+1
        nempmx = nemupx-nemlox+1
      elseif (lpart0 == 2) then
        nfilmx = 1
        nempmx = 1
      endif

C ... Initialization : get demin,demax,iwgt,nwgt
      nqbz = n1*n2*n3
      npm = 1                   ! Assume time reversal symmetry for now
      nctot = 0                 ! No cores for now
      allocate(ecore(1,1))      ! dummy core array
      if (npm == 1) then
        ncc = 0
      else
        ncc = nctot
      endif

      i = nband+nctot; j = nband+ncc
      idum = idalloc('demin',allocvb()+2,i*j,nqbz*npm*2)
      allocate(demin(i,j,nqbz,npm),demax(i,j,nqbz,npm))
      allocate(iwgt(i,j,nqbz,npm))
      allocate(nwgt(nqbz,npm),nwgttt(nqbz,npm))

      i = 0
      if (loptic == -5 .or. loptic == -6) i = 10
      call tetwtq(i,npm,ncc,.false.,esciss,qlat,1d0,
     .  efermi,eocc,eunocc,nctot,ecore(1,1),
     .  ntet,nqbzw,nband,nqbz,idtet,qbzw,ib1bz,
     x  xv,[0],1,[0],1,[0],[0],[0],
     x  q,iq,1,1,nkp,
     o  iwgt,nwgt,demin,demax,
     x  [0d0])

C ... Sanity check, and global maximum demaxx of demax
C     Note!  demin,demax are in Hartree
      demaxx = -9999
      j = nband+ncc
      if (loptic == -5 .or. loptic == -6) j = 1
      do  i1 = 1, nband+nctot
        do  i2 = 1, j
          do  i3 = 1, nqbz
            if (iwgt(i1,i2,i3,1)) then
C             print *, i1,i2,i3, demin(i1,i2,i3,1), demax(i1,i2,i3,1)
              if (abs(demin(i1,i2,i3,1)) == 1d10 .or.
     .            abs(demax(i1,i2,i3,1)) == 1d10)
     .          call rx('bug in tetwtq, job=0')
              if (demax(i1,i2,i3,1) < demin(i1,i2,i3,1))
     .          call rx('bug in tetwtq, job=0')
              demaxx = max(demaxx,2*demax(i1,i2,i3,1))
            endif
          enddo
        enddo
      enddo

C ... Make indices ihw,jhw,nhw for tetwtq
      nwgtx = maxval(nwgt(1:nqbz,1:npm))
      call pralloc('OPTINQ','nib',2*nwgtx*dble(nqbz)*npm,100,30,0)
      idum = idalloc('nib',allocvb()+2,nwgtx*nqbz*npm,1)
      allocate(nib(nwgtx,nqbz,npm),njb(nwgtx,nqbz,npm))
C     nwgttt should be the same as nwgt
      do  jpm = 1, npm
        call pair_index(jpm,iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc,
     .    nwgtx,nib(1,1,jpm),njb(1,1,jpm),noccxvx(jpm),nwgttt(1,jpm))
      enddo
      deallocate(nwgttt)
C     call pralloc('OPTINQ','ihw',3*nwgtx*dble(nqbz)*npm,100,30,0)
      idum = idalloc('ihw',allocvb()+2,nwgtx*nqbz*3/2,npm)
      allocate(ihw(nwgtx,nqbz,npm),jhw(nwgtx,nqbz,npm),
     .         nhw(nwgtx,nqbz,npm))

C ... Energy histogram frhis. frhis must contain enough points to:
C     encompass largest demax
C     encompass emax
C     span the range(emin,emax)
      optrng(3) = max(emin,0d0)
      optrng(4) = max(dble(demaxx),emax-min(emin,0d0))
CC     Add a little margin following hx0fp0
C      if (de(2) /= 0) then
C        nwhis = int(de(2)/de(1)*
C     .    (dsqrt(1d0+2*(optrng(4)-optrng(3))/de(2))-1d0))
C     .    +1+3 !+3 for margin
C      endif

      i = 1; if (de(2) /= 0) i = 11
      call dosmsh(-i,optrng(3),de,0,nwhis,xv) ! Return nwhis
      allocate(frhis(nwhis+1))
      call dosmsh(100+i,optrng(3),de,nwhis+1,idum,frhis)
      if (nwhis < npts) call rx('optinq: frhis improperly dimensioned')
C      if (de(2) == 0) then
C        frhis(nwhis+1) = max(emin,0d0) + de(1)*nwhis
C      else
C        frhis(nwhis+1) = max(emin,0d0) +
C     .    de(1)*nwhis + de(1)**2/de(2)/2*nwhis**2
C      endif
C     frhis must be in Hartrees for hisrange and tetwtq
      call dscal(nwhis+1,.5d0,frhis,1)

C ... Range of bands for given window
      jhwtot = 1
      do  jpm = 1, npm
        do  k = 1, nqbz
          do  i = 1, nwgt(k,jpm)
            i1 = nib(i,k,jpm)
            i2 = njb(i,k,jpm)
            call hisrange(frhis,nwhis,
     .        demin(nib(i,k,jpm),njb(i,k,jpm),k,jpm),
     .        demax(nib(i,k,jpm),njb(i,k,jpm),k,jpm),
     .        ihw(i,k,jpm),nhw(i,k,jpm))
            jhw(i,k,jpm) = jhwtot
            jhwtot = jhwtot + nhw(i,k,jpm)
            if (nhw(i,k,jpm) < 0 .or. jhwtot < 0) then
              call info2(0,0,0,' OPTINQ:  aborting at k = %i of '//
     .          'nqfbz = %i',k,nqbz)
              call rxi('OPTINQ integer overflow: '//
     .          'too many bin allocations n =',jhwtot-nhw(i,k,jpm))
            endif
          enddo
        enddo
      enddo
      deallocate(demin,demax)
      idum = idalloc('demin',allocvb()+4,0,0)

C ... Allocate whw; make ibjb
      nhwtot = jhwtot-1
      i = nband+nctot; j = nband+ncc
      idum = idalloc('whw',allocvb()+2,nhwtot,1)
      idum = idalloc('ibjb',allocvb()+2,i*j/2,nqbz*npm)
      ltmp = memcheck('optinq','optics band weights',s_ctrl%maxmem,T)
      allocate(whw(nhwtot),ibjb(i,j,nqbz,npm))
      whw = 0d0
      ibjb = 0
      do  jpm = 1, npm
        do  k = 1, nqbz
          do  i = 1, nwgt(k,jpm)
            i1 = nib(i,k,jpm)
            i2 = njb(i,k,jpm)
            ibjb(i1,i2,k,jpm) = i
          enddo
        enddo
      enddo
      deallocate(nib,njb)
      idum = idalloc('nib',allocvb()+4,0,0)

C ... Tetrahedron weights, resolved by (ib,jb) pair and kpoint
      i = 1
      if (loptic == -5 .or. loptic == -6) i = 11
      call tetwtq(i,npm,ncc,.false.,esciss,qlat,1d0,
     .  efermi,eocc,eunocc,nctot,ecore,
     .  ntet,nqbzw,nband,nqbz,idtet,qbzw,ib1bz,
     .  frhis,nwhis,nwgtx,ibjb,nhwtot,ihw,nhw,jhw,
     x  q,iq,1,1,nkp,
     i  iwgt,nwgt,xv,xv,
     o  whw)
      deallocate(efbz1,ecore,iwgt,nwgt)
      if (npass == 1 .or. loptic == 2) deallocate(ib1bz,qbzw,idtet)

      allocate(wtthis(nwhis,npol))
C     Default number of (occ,unocc,k) channels for partial decomp'sn
      nfilmw = 1; nempw = 1; nkpw = 1
      if (lpart0 > 0) then
C       Number of (k) channels for partial decomposition
        if (lpart0 > 1) nkpw = nkp
        if (lpart0 == 1 .or. lpart0 == 3) then
          nfilmw = nfilm; nempw = nempm
        endif
        allocate(wtphis(nwhis,npol,nfilmw,nempw,nkpw))
        call dpzero(wtphis,nwhis*npol*nfilmw*nempw*nkpw)
      else
        nkpw = 0
        allocate(wtphis(1,1,1,1,1))
      endif

      i = 0
      if (loptic > 0) i = 1
      if (loptic < 0 .and. njdosw > 0) i = 2
      call tetwtt(lpart0*10+i,nfilox,nfiupx,nemlox,nemupx,npm,n1,n2,n3,
     .  nqbz,ipq,nwhis,nwgtx,ibjb,nhwtot,ihw,nhw,jhw,whw,nfilo,nfiup,
     .  nemlo,nemup,nfilmw,nempw,nsp,optmt(1,1,1,isp,1),njdosw,
     .  jdosw(1,1,1+nkp*(isp-1)),npol,nkp,wtthis,wtphis)

C ... Optics branch: double Im eps to include both spins if nsp=1
      fac = 1/pi
      if (loptic > 0 .and. nsp == 1) fac=fac*2
C     if (loptic > 0 .and. nsp == 1) xx = xx/2

C     Restore frhis to Ry scale
      call dscal(nwhis+1,2d0,frhis,1)

C     emesh : first npts of frhis, shifted by optrng(1)-optrng(3)
      allocate(emesh(npts))  ! frhis has at least npts+1 points
C     forall (i=1:npts) emesh(i) = optrng(1)-optrng(3) + (frhis(i)+frhis(i+1))/2d0
      forall (i=1:npts) emesh(i) = optrng(1)-optrng(3) + frhis(i)

C     Interpolate histogram to frhis mesh
      if (loptic == -5 .or. loptic == -6) then
C        do  i = 1, npts
C          opt(npts-i+1,1,ksp) = wtthis(i,1)/xx/pi
C        enddo
        do  ipol = 1, npol
          opt(1,ipol,ksp) = 0
          call optind(npts,fac,frhis,wtthis(1,ipol),wtthis(1,ipol))
          do  i = 1, npts-1
            opt(npts-i+1,ipol,ksp) = wtthis(i,ipol)
          enddo
        enddo
      else
        do  ipol = 1, npol
C         opt(1:npts,ipol,ksp) = wtthis(1:npts,ipol)/xx/pi
          call optind(npts,fac,frhis,wtthis(1,ipol),opt(1,ipol,ksp))
        enddo

      endif
      deallocate(ihw,nhw,jhw,ibjb,whw,wtthis)
      idum = idalloc('ihw',allocvb()+4,0,0)
      idum = idalloc('whw',allocvb()+4,0,0)
      idum = idalloc('ibjb',allocvb()+4,0,0)
C     idum = idalloc(' ',10,1,1)

C ... Interpolate partial decomposition
      if (lpart0 > 0) then
        if (nfilm /= nfiup-nfilo+1)
     .    call rx('improper dimensioning parameter nfilm')
        if (nempm /= nemup-nemlo+1)
     .    call rx('improper dimensioning parameter nempm')

        allocate(popt(nptmx,npol,nfilmw,nempw,nkpw))
        call dpzero(popt,nptmx*npol*nfilmw*nempw*nkpw)

        do  ipol = 1, npol
        do  ibw = 1, nfilmw
        do  jbw = 1, nempw
        do  iq1 = 1, nkpw
          if (loptic == -5 .or. loptic == -6) then ! Simple DOS
C            do  i = 1, npts
C              popt(npts-i+1,ipol,ibw,jbw,iq1) =
C     .          wtphis(i,ipol,ibw,jbw,iq1)/xx/pi
C            enddo
            popt(1,ipol,ibw,jbw,iq1) = 0
            call optind(npts-1,fac,frhis,wtphis(1,ipol,ibw,jbw,iq1),
     .        wtphis(1,ipol,ibw,jbw,iq1))
            do  i = 1, npts-1
             popt(npts-i+1,ipol,ibw,jbw,iq1) =wtphis(i,ipol,ibw,jbw,iq1)
            enddo

          else
C            popt(1:npts,ipol,ibw,jbw,iq1) =
C     .      wtphis(1:npts,ipol,ibw,jbw,iq1)/xx/pi
            call optind(npts,fac,frhis,wtphis(1,ipol,ibw,jbw,iq1),
     .        popt(1,ipol,ibw,jbw,iq1))
          endif
        enddo
        enddo
        enddo
        enddo
      endif
      if (allocated(wtphis)) deallocate(wtphis)
      if (allocated(popt) .and. ksp /= npass) deallocate(popt)
      if (allocated(frhis)) deallocate(frhis)

C ... Optics for 2nd spin if conditions apply
      if (loptic == 1 .and. npass == 2) then
        loptic = 2
        ksp = 2
        if (allocated(emesh)) deallocate(emesh)
        goto 10
      endif

      i = 13
      if (loptic < 0) i = 1
      if (loptic == -5 .or. loptic == -6) i = 21
      if (nkpw == 0) allocate(popt(1,1,1,1,1)) ! so compiler is happy
      call optin2(s_optic,s_bz,i,ksp,nspc,vol*alat**3,nptmx,npol,
     .  lpart1 /= 0,nfilmw*nempw*nkpw,npts,emesh,popt,opt)

      deallocate(emesh)

      efermi = saveef
      call tcx('optinq')
      end
      subroutine optind(npts,fac,frhis,hisn,hisd)
C- Generate NOS and DOS from histogram
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :Not implemented
Ci   npts  :number of DOS tabulation points (input; for sampling only)
Ci   fac   :scale hisd by constant
Ci   frhis :histogram energies.  Note that it must have one extra point
Ci   hisn  :histogram NOS
Co Outputs
Ci   hisd  :histogram DOS (may occupy same address as hisn)
Cl Local variables
Cl         :
Cr Remarks
Cr   This routine turns the histogram generated from tetwtq into an
Cr   integrated quantity, and differentiates it.
Cu Updates
Cu   18 Sep 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npts
      double precision fac,frhis(npts+1),hisn(npts),hisd(npts)
C ... Local parameters
      integer i
      double precision sum,wk(npts+1),y2(npts+1),x,y

C ... Accumulate histogram NOS to make NOS at frhis points
      sum = 0
      wk(1) = 0
      do  i = 1, npts
        sum = sum + hisn(i)
        wk(i+1) = sum
      enddo
      do  i = 2, npts+1
        wk(i) = wk(i)*fac
      enddo

C ... Get second derivatives for spline
      call cspline(frhis,wk,npts+1,1d99,1d99,y2)

C ... Make DOS on frhis mesh
      do  i = 1, npts
        x = frhis(i)
        call csplint(frhis,wk,y2,npts,x,y,hisd(i))
      enddo

      end
