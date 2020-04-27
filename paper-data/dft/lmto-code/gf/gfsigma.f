      subroutine gfsigma(s_ctrl,s_ham,s_pot,s_lat,s_bz,offH,
     .  nbas,nkp,qp,nk1,nk2,nk3,pas,salp,iax,nttab,sxad,tolsig,gamm)
C- Repacking of sigma in sigm-file
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl nspin lgen3 nclass
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham hord
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  lshft zval
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nbas  :size of basis
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   zpp    :complex energy
Ci   wzp    :weights for complex energy integration
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   pas   :number of sites considered in the slice of sigma.
Ci   sxad  :T if sigma already exists on disk
Co Outputs
Co  gfsigma :sigma written on disk.  sigm taken from file 'sgm'
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   30 Oct 04 (T. Sandu) Updated to handle spin-polarized systems
Cu   10 Apr 04 (T. Sandu) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed variables
      integer nbas,nkp,nk1,nk2,nk3
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas+1)
      double precision qp(3,nkp)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      complex(8), allocatable :: hc(:)
      real(8), allocatable :: gfsigmal(:)
C Local variables
      integer hord,ib,ib1,ib2,ldham(16),ldim,ldimx,lgen3,lham,lhdim,
     .  lhdimx,lidim,lidimx,lncol,nl,nsp,nspc,
     .  iblk,nsgrp,nlmaa,mxorb,nglob
c      integer oipc,onrc
      double precision plat(3,3)
      equivalence (ldim,ldham(1)),(ldimx,ldham(5)),(nspc,ldham(4))
      equivalence (lidim,ldham(2)),(lidimx,ldham(6))
      equivalence (lhdim,ldham(3)),(lhdimx,ldham(7))
C ... For file I/O of gf
      integer clp(9,2)
C ... For file I/O of gfwsq
c      integer ifis
      integer fisgm,iq1,iq,isp,isp1
c      integer kcplx
      integer fopna, ifisq,nlmac,nlmac1
      integer lshft(3),pas,pas1
      double precision zval
      character*8 fisnam
c      double complex wsqlm(lidim,lidim,nsp)
C      integer ISH
C      parameter (ISH=1)
      double precision qp0(3)
C...start new slfconsistency
      integer oh2,oh3,oh4,ig,ormat,ld2
C....new selfconsistency
      integer nsite
      integer osll,lbloch,bit
      integer idim,nl2,nclass,nlspc,li,i2,l2
      logical ltmp,bittst,iostr
      double precision vmtz,ckbas,cksumf,kap2(20)

      double precision g(3,3),ag(3)
      logical sxad
      integer ifisold,osigold,gamm,opotp,opotph,oo,niax,nttab
      double precision tolsig,qp1(3),xx,salp(*)
      parameter (niax=10)
      integer iax(niax,nttab)
      double precision qp2(3),qpoff(3)

      call tcn('gfsigma')

      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lgen3 = s_ctrl%lgen3
      nclass = s_ctrl%nclass
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      ldham = s_ham%ldham
      hord = s_ham%hord
      lshft = s_bz%lshft
      zval = s_bz%zval
C ... selfconsistent
      vmtz = s_pot%vmtz
      ckbas = cksumf(s_lat%pos,3*nbas)
c      kcplx = 0
      mxorb = nglob('mxorb')

      if (nspc == 2) call rx('gfsigma: not ready for non-collinear')

C     Possibly rotate to spherical harmonics when making Bloch sum
      lbloch = 60000
      if (bittst(lham,256)) lbloch = lbloch+1000

C     Some dimensioning parameters and memory allocation
C     List of sites

      ib1 = 1
      ib2 = nbas
      lidim = offH(4,1,nbas+1)
      call iinit(clp,9*2)
      clp(3,1) = lidim*nspc
      clp(4,1) = lidim*nspc
      clp(5,1) = lidim*nspc
      nl2 = nl**2
      idim = lidim-ldim
      nlspc = nl*nsp*nclass
      l2 = ldim**2
      i2 = idim**2
      li = ldim * idim

C   making potential parameters in alpha representation
      print *,'into gfsigma'

cccccc      call getnewpp(nl,nsp,nclass,s_pot%pp,w(opotp))
      print*,'aici 0'
cccccc      call makpph(nl,nsp,nbas,lhdim,s_ctrl%ipc,s_ham%iprmb,w(opotp),
cccccc     .    w(opotph))
c      call makpph(nl,nsp,nbas,lhdim,s_ctrl%ipc,s_ham%iprmb,s_pot%pp,
c     .    w(opotph))


C-...openig of the source file for sigma
      fisnam = 'sgm'
      fisgm = fopna(fisnam,-1,4)
C   ... rewind file, read header
      rewind fisgm

c        call yprm('str-q gf1kp',2,gf,lidim*lidim,lidim,lidim,lidim)
C.....sigma file
C ... Write sigma header
      ifisq = fopna('SIGM',-1,4)
      rewind ifisq
      call dpzero(qpoff,3)
      call iosigh(3,0,nsp,1,ldim,ldim,nk1,nk2,nk3,nkp,nkp,lshft(1),lshft(2),lshft(3),-ifisq,0d0)

      tolsig = 0d0

C ... If sxad, read head from either sigm-old (gam>0) or sigm1 (gam=0)
C      if (sxad .and. (gamm /= 0)) then
C      ifisold = fopna('SIGM-OLD',-1,4)
C      rewind ifisold
C      call iosigh(3,0,nsp,1,ldim,ldim,nk1,nk2,nk3,nkp,nkp,lshft(1),lshft(2),lshft(3),ifisold)
C
C      elseif (sxad .and. (gamm == 0)) then
C      ifisold = fopna('SIGM1',-1,4)
C      rewind ifisold
C      call iosigh(3,0,nsp,ldim,nk1,nk2,nk3,nkp,lshft(1),lshft(2),,lshft(3),ifisold)
C      tolsig = 1d0
C
C      endif

C ... Write sigma to disk, reading slices from file sgm
C     'sgm' writes records in order ((sig(i,iq), iq=1..nq), i=1,nbas)
C     Double k-loop needed to read (sig(i,iq), i=1,nbas) for 1 qp
      do  isp = 1, nsp
      do  iq = 1, nkp
      rewind fisgm
      pas1 = pas-1
C     Allocate memory for 1 k-point sigma new
      allocate(gfsigmal(lidim*lidim*2*nspc*nspc))
      call dpzero(gfsigmal,lidim*lidim*2*nspc*nspc)
C     Old sigma file

      do  isp1 = 1, nsp
      nlmac = 0
      iblk = 0
C --- For each site, extract ib-slice for iq and merge ---
      do  ib = ib1, ib2, pas

C --- sgm file for all qp in subblock connected to site ib ---
C ... Get gf for subblock connected to site ib

      if (iblk < ib) then
        iblk = min(ib+pas1,ib2)
        nlmaa = offH(4,1,iblk+1) - offH(4,1,ib)
        nlmac1 = nlmac
        nlmac = nlmac + nlmaa
        allocate(hc(nlmaa*lidim*nspc*nspc));
        call dpzero(hc,2*nlmaa*lidim*nspc*nspc)
      endif

C ... Assemble sigma for all sites for this qp
      do  iq1 = 1, nkp
        call dpdump(hc,nlmaa*lidim*nspc*2,fisgm)
C       call yprm('reading',2,hc,nlmaa*lidim,nlmaa,nlmaa,lidim)
        if (iq1 == iq .and. isp1 == isp) then
C         call yprm('bfr-sgfll-sig',2,hc,nlmaa*lidim,nlmaa,nlmaa,lidim)
          call gfsigfll(nlmac1,nlmaa,lidim,hc,gfsigmal)
        endif
      enddo

      if (ib == iblk) deallocate(hc)
      enddo
      enddo

C-.....writing on disk
c      print*,'before dump'
      call dpscop(qp,qp0,3,3*iq-2,1,1d0)

C....  making structure constant in k-space to transf. sigma to alpha
cccc      ltmp = iostrx(8,'STR',nl,nbas,1,kap2,0,ckbas,-1,nsite,oalph,
cccc     .  oiax,ontab,os)
cccc      nsite = ntab(nbas)
ccc      call zprm('str-rs',2,w(os),lidim,lidim,lidim)
cc      call bloch(lbloch,qp0,nl,plat,nl2,s_ham%iprmb,1,nsite,w(oiax),w(os),
cc     .  nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(osll),w(osil),w(osii))
ccc      call bloch(lbloch,qp0,nl,plat,nl2,s_ham%iprmb,1,nsite,w(oiax),w(os),
ccc     .  nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(osll),w,w)
      print*,'qp0',qp0
      print*,'nttab',nttab
      print*,'lbloch',lbloch
cccccc      call bloch(lbloch,qp0,nl,plat,nl2,s_ham%iprmb,1,nttab,iax,salp,
cccccc     .  nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(osll),w,w)
cc      print *,'lbloch',lbloch
c      call yprm('str-q',12,w(osll),ldim*ldim,ldim,ldim,ldim)

C Clean this up!

C ... sigm for 1 qp is written as 2 records: qp, then sigm (iosigh)
      write(ifisq) qp0

C ... Generated sigma = sigma(alpha), no prior sigma
      if ((gamm == 0) .and. (.not. sxad)) then
       print*,'last call .not. sxad & gfsigma gamm shd 0=',gamm
cc       call yprm('sig-q',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
       call dpdump(gfsigmal,lidim*lidim*2,-ifisq)

C ... Generated sigma = sigma(alpha), but prior sigma
      elseif ((gamm == 0) .and. sxad) then
       print*,' gfsigma sxad &  gamm shd 0=',gamm
       read(ifisold) qp2
       call dpdump(gfsigmal,lidim*lidim*2,ifisold)
cc       call yprm('sig-q',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
ccc       print *,'before sxsgmtrns'
cccccc       call sxsgmtrns(isp,nsp,gfsigmal,ldim,lhdim,w(opotph),w(osll))
ccc       print *,'after sxsgmtrns'
cccccc       call yprm('sig-2',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
c       call yprm('sig-q',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
       call dpdump(gfsigmal,lidim*lidim*2,-ifisq)

C ... Generated sigma = sigma(gamma)
      else

C       call yprm('sig-q',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
C       sigma(gam) -> sigm(alpha)
C       call sxsgmtrns(isp,nsp,gfsigmal,ldim,lhdim,w(opotph),w(osll))
C       call yprm('sig-2',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
        call dpdump(gfsigmal,lidim*lidim*2,-ifisq)

C       Compare new, old sigma
C       if (sxad) then
C         read(ifisold) qp1
C         call dpdump(w(osigold),lidim*lidim*2,ifisold)
Cc        call yprm('sig-new',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
Cc        call yprm('sig-old',2,w(osigold),lidim*lidim,lidim,lidim,lidim)
CC ...difference new_sigma - old_sigma
C         call delsigerr(gfsigmal,w(osigold),lidim,nkp,nsp,tolsig)
C         print*,'tolsig',tolsig
C         call rlse(osigold)
C       endif
      endif

      deallocate(gfsigmal)
      enddo
      enddo

C     End of file save

C      if ((gamm == 0) .and. (.not. sxad)) then
C        call rlse(ogfsigmal)
CCCC        call rlse(opotph)
C        call fclose(fisgm)
C        call fclose(ifisq)
C        print *, 'gfsigma at end gamma repres'
C
C        call tcx('gfsigma')
C        return
C      elseif ((gamm == 0) .and. (sxad)) then
C        call rlse(ogfsigmal)
CCCC        call rlse(opotph)
C        call fclose(fisgm)
C        call fclose(ifisq)
C       call fclose(ifisold)
C        print *, 'gfsigma at end gamma repres'
C
C        call tcx('gfsigma')
C        return
C      else
C
CC.....sigma file
CC ... read sigma header to write it into sigma old for next iteration
C       rewind ifisq
Cc       print*,'after rewinding sigm-ifisq'
C
C      call iosigh(3,0,nsp,1,ldim,ldim,nk1,nk2,nk3,nkp,nkp,lshft(1),lshft(2),lshft(3),ifisq)
C
Cc       print*,'after reading header sigm-ifisq'
C
C       if (sxad) then
CC ... Write old sigma file header
C        rewind ifisold
C
Cc        print*,'after rewinding sigm-old-ifisold'
C
C        call iosigh(3,0,nsp,1,ldim,ldim,nk1,nk2,nk3,nkp,nkp,
C     .  lshft(1),lshft(2),lshft(3),-ifisold)
C
Cc        print*,'after writing header sigm-ifisq'
C
C       else
CC ....create old sigma file
CC ... writes old sigma file header
C        ifisold = fopna('SIGM-OLD',-1,4)
Cc        print*,'after opening-creating sigm old-ifisold'
C        rewind ifisold
C        tolsig = 1d0
C
Cc        print*,'after rewinding sigm-old-ifisold'
C
C        call iosigh(3,0,nsp,1,ldim,ldim,nk1,nk2,nk3,nkp,nkp,
C     .   lshft(1),lshft(2),lshft(3),-ifisold)
Cc        print*,'after writing header sigm-old-ifisold'
C       endif
CC ... reads the new sigma from SIGM and writes it into SIGM-OLD
C      do  600 isp = 1, nsp
C      do  600 iq = 1, nkp
Cc       print*,'nkp iq',nkp,iq,'nsp isp',nsp,isp
Cc       print*,'reading/dumping '
C       read(ifisq) qp0
C       call dpdump(gfsigmal,lidim*lidim*2,ifisq)
Cc       call yprm('sgnew-rd',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
Cc       print*,'between reading writing'
C       write(ifisold) qp0
C       call dpdump(gfsigmal,lidim*lidim*2,-ifisold)
Cc       call yprm('sgolddmp',2,gfsigmal,lidim*lidim,lidim,lidim,lidim)
C 600  continue
Cc      print*,'gfsigma tolsig=',tolsig
C      endif
C      call rlse(ogfsigmal)
Ccc      call rlse(osigold)
CCCC      call rlse(opotph)


      call fclose(fisgm)
      call fclose(ifisq)
C      call fclose(ifisold)

      call tcx('gfsigma')

      end
       subroutine delsigerr(signew,sigold,lidim,nkp,nsp,tolsig)
C- Calculate delta sigma from one iteration to the old
      integer lidim,nkp,nsp
      double precision signew(lidim,lidim,2),sigold(lidim,lidim,2)
      double precision tolsig,diff,sumsig
      integer i,j,k
cc      tolsig = 0d0
      sumsig = 0d0
      diff = 0d0
      do  100  k = 1, 2
       do  100  j = 1, lidim
        do  100  i = 1, lidim
          diff = diff + (signew(i,j,k) - sigold(i,j,k))**2
          sumsig = sumsig + sigold(i,j,k)**2
cccc          tolsig = tolsig + diff*diff/(nkp*nsp)
cc          print*,'ijk',i,j,k
cc          print*,'diff=',diff
cc      call yprm('signew',2,signew,lidim*lidim,lidim,lidim,lidim)
cc      call yprm('sigold',2,sigold,lidim*lidim,lidim,lidim,lidim)
cc          print*,'tolsig=',tolsig
 100  continue
      tolsig = tolsig + diff/sumsig
cc      print*,'tolsig=',tolsig
c      call yprm('signew',2,signew,lidim*lidim,lidim,lidim,lidim)
c      call yprm('sigold',2,sigold,lidim*lidim,lidim,lidim,lidim)
      end


      subroutine gfsigfll(nlmac1,nlma,lidim,hr,sgm)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma,lidim  :dimensions of a slice of h
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   iq    :current k-point
Ci   hc    :h matrix in complex form
Co Outputs
Co    hr   :h matrix in real form
Cl Local variables
Cr Remarks
Cr Made April 15 2004
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlmac1,nlma,lidim
      double precision sgm(lidim,lidim,2)
C ... Local parameters
      integer ii,ij
      double precision hr(nlma,lidim,2)

      do  100  ij = 1, lidim
      do  100  ii = 1, nlma
      sgm(nlmac1+ii,ij,1) = hr(ii,ij,1)
      sgm(nlmac1+ii,ij,2) = hr(ii,ij,2)
c      print *,'nlmac1+ii',nlmac1+ii
  100 continue
c      call yprm('gfsigfll',2,sgm,lidim*lidim,lidim,lidim,lidim)
c      print *,'nlmac1',nlmac1

      end

      subroutine getnewpp(nl,nsp,nclass,ppold,ppnew)
      implicit none
      integer nl,nsp,nclass
      double precision ppold(6,nl,nsp,nclass),ppnew(6,nl,nsp,nclass)
      integer i,j,k,l
      do 1 i = 1, 6
      do 1 j = 1,nl
      do 1 k = 1,nsp
      do 1 l = 1, nclass

       ppnew(i,j,k,l) = ppold(i,j,k,l)

 1    continue

      end
