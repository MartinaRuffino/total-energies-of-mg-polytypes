      subroutine optint(s_optic,s_bz,lteto,eband,nbmax,ksp,nsp,
     .  nspc,ikp,efermi,idtet,wgts,vol,nfilm,nempm,nptmx,npol,optmt,opt)
C- BZ integration of Im(eps) or joint DOS
C ----------------------------------------------------------------------
Cio Structures
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic lpart dw window ocrng unrng esciss
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet lmet n range w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  optin2
Ci Inputs:
Ci   lteto :Integration mode
Ci         : 0 sampling
Ci         : 1 original tetrahedron integration (now merged with lteto=2)
Ci         : 2 tetrahedron, calling bzjnos
Ci         : 3 tetrahedron, calling tetwtq
Ci         :-1 on-the-fly sampling
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   ksp   :starting spin index passed to mkjdos.  Used for
Ci         :on-the-fly sampling integration (partial contribution),
Ci         :when only information for second-spin is available.
Ci         :See Remarks
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   ikp   :current k-point (sampling integration only)
Ci         :If ikp>0, partial contribution from a single
Ci         :k-point (ikp) and spin (ksp) accumulated (see Remarks)
Ci         :ikp<0 signifies that information from all k-points is
Ci         :available to be accumulated.
Ci   efermi:Fermi level
Ci         :NB: efermi=-99 is a flag, indicating that the Fermi
Ci         :level is not set.
Ci   idtet :(0,i) no. of tetrahedra of the i'th kind (tetirr)
Ci         :(1-4,i) identifies the i'th tetrahedron in terms of
Ci         :four irreducible k-points
Ci   vol   :volume
Ci   nfilm,nempm: dimensions optmt.  Should be nfiup-nfilo+1,nemup-nemlo+1
Ci   nptmx :max number of bins; dimension for following 3 work arrays
Ci   npol  :number of polarizations; needed for dimensioning.
Ci         :npol should be 1 for jdos, 3 for eps.
Ci   optmt :<i|grad|j> connecting occ i with unocc j
Ci         :optmt is not used when making joint DOS (loptc<0)
Co Outputs:
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
Co         :Caller should set opt=0 before first call to optint.
Cr Remarks
Cr   Adapted from bzints to make joint density of states or Im(eps)
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cr
Cr   (Sampling only) The partial contribution to Im(eps) may be made from
Cr   a single k-point and spin.  This routine accumulates contributions from
Cr   all k-points uness ikp>0, in which case it accumulates a
Cr   partial contribution from (ikp,ksp) only.
Cl Local variables and switches
Cl  loptic  What optical property to calculate
Cl          Should be same as s_ctrl%loptc
Cl          1 generate Im(eps)
Cl          2 same as 1 in this routine
Cl          8 Im eps for  nonequilibrium absorption : independent Ef for electrons and holes
Cl            Transition probability scaled by f(E-imrefp) * [1-f(E-imrefn)]
Cl          9 Im eps for  nonequilibrium emission : independent Ef for electrons and holes
Cl            Transition probability scaled by [1-f(E-imrefp)] * f(E-imrefn)
Cl         Add 20 for Second Harmonic Generation
Cl         -1 generate joint DOS
Cl         -2: generate joint DOS, spin 2 (use with LTET=3)
Cl         -3: generate up-down joint DOS (LTET=3)
Cl         -4: generate down-up joint DOS (LTET=3)
Cl         -5: generate spin-up single DOS (LTET=3)
Cl         -6: generate spin-dn single DOS (LTET=3)
Cl         -8 Joint DOS for  nonequilibrium absorption : independent Ef for electrons and holes
Cl            Transition probability scaled by f(E-imrefp) * [1-f(E-imrefn)]
Cl         -9 Joint DOS for  nonequilibrium emission : independent Ef for electrons and holes
Cl            Transition probability scaled by [1-f(E-imrefp)] * f(E-imrefn)
Cl   popt  :opt resolved by (unocc,occ) pairs.
Cu Updates
Cu   12 Oct 16 If s_optic%imref(1) is not NULL, use s_optic%imref(1) for the Fermi level
Cu   18 Dec 14 New modes loptic=8,9,-8,-9
Cu   06 Jan 14 Merged lteto=1 into lteto=2
Cu             lteto=1 had partial decompsn', which is now lost.
Cu             For this feature, use lteto=3 (built into optinq, not optint)
Cu   10 Nov 11 Begin migration to f90 structures
Cu   31 Dec 10 optmt from optdme is now calculated for all (nfilm,nempm)
Cu   11 Sep 09 Some reworking to accomodate future tetwtq
Cu   15 Jan 04 (Uk-Jin Roh) add partial sampling contribution to Im(eps)
Cu   20 Nov 02 (jek) extended to metals
Cu   02 Mar 01 Added scissors operator
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lteto
      integer nbmax,nsp,nspc,nfilm,nempm,nptmx,ksp,ikp,npol
      integer idtet(0:4,*)
      double precision efermi,vol
      double precision opt(nptmx,npol,nsp)
      double precision optmt(3,nfilm,nempm,nsp,*),wgts(*),
     .                 eband(nbmax,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Local parameters
      logical jdos,metal,lonekp
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet,ibf,ibm,k,ip,lmet
      integer fopn,lgunit,isw,loptic,lpart,lpart0,lpart1,ifil,ipr,stdo
      integer ib,jb,isp,ocrng(2),unrng(2),npts
      integer nkabc(3),n1,n2,n3
      integer mpsord,iq,ifmin,ifmax,iwmin,iwmax,npopt
      real(8),allocatable :: popt(:,:,:,:),emesh(:)
C     real(8),allocatable :: yy(:),optpar(:,:)
      double precision swidth,srnge,del
      double precision optrng(4),emin,emax,fac,pi,ds2,ds3,de(3),volwgt
      double precision esciss
      logical cmdopt,samp,bonly
      character*120 strn
      integer,parameter :: NULLI =-99999
      real(8),parameter :: NULLR =-99999
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
C     integer iq1,iq2,iq3,iq4,ipts,itet
C     double precision wt,eci(4),ecj(4),wtm
C ... External calls
      external awrit2,awrit3,awrit8,bzjnos,dosmsh,dpzero,dscal,fclose,getpr,
     .  info0,info2,info5,inta,mkjdos,optin2,rx2,rxx,tcn,tcx,xxxdif,ywrm

C ... Setup
      call tcn('optint')
      call getpr(ipr)
      stdo = lgunit(1)
C     loptic = s_ctrl%loptc
      loptic = s_optic%loptic
      lpart = s_optic%lpart
      de = s_optic%dw
      optrng = s_optic%window
      ocrng = s_optic%ocrng(1:2)
      unrng = s_optic%unrng(1:2)
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
      lmet = s_bz%lmet
      esciss = s_optic%esciss
      mpsord = s_optic%nmp; if (mpsord == NULLI) mpsord = s_bz%n
      swidth = s_optic%w; if (swidth == NULLI) swidth = s_bz%w
      srnge = s_bz%range
      if (iabs(loptic) == 8 .or. iabs(loptic) == 9) then
        if (s_optic%imref(1) == NULLR) s_optic%imref(1) = efermi
        if (s_optic%imref(2) == NULLR) s_optic%imref(2) = efermi
      endif
      if (iabs(loptic) == 1 .or. iabs(loptic) == 2) then
        if (s_optic%imref(1) /= NULLR) efermi = s_optic%imref(1)
      endif

C ... Mesh setup
      k = 1
      if (de(2) /= 0) k = 11
      if (k == 11) call rx('optint does not allow quadratic meshes')
!       call dosmsh(-k,optrng,de,0,npts,fac)
!       if (npts > nptmx)
!      .  call rxi('OPTINT: nptmx too small: need at least',npts)
!       allocate(emesh(npts))
!--- If the memory savings are indeed important then use additional
!    variable or add 'volatile' attribute to nptmx in dosmsh to notify the compiler of the aliasing. (dmt)
      allocate(emesh(nptmx))
      call dosmsh(k,optrng,de,nptmx,npts,emesh)

      jdos = loptic < 0
      mpsord = mod(mpsord,100)
      samp = lteto <= 0
      lonekp = ikp > 0
      metal = lmet /= 0
      npopt = 0
      lpart0 = mod(lpart,10)
      lpart1 = mod(lpart/10,10)
      if (lpart0 /= 0) then
        npopt = nfilm*nempm
      else
        allocate(popt(1,1,1,1))
      endif

C      call rxx(metal .and. dabs(esciss) > 1d-8,
C     .  'OPTINT: metal not compatible with scissors operator')
      call rxx(npts > nptmx,' OPTINT: npts gt nptmx')
      call rxx(npts < 1,' OPTINT: npts < 1')
      if (lonekp .and. lteto >= 0) call rx
     .  ('optint: on-the-fly integration but ltet0>=0')

C      print *, '!!'
C      efermi = 0
C      eband(1,:,1:nkp) = efermi - .1d0
C      eband(2,:,1:nkp) = efermi + .1d0
C      eband(3:,:,1:nkp) = efermi + 999
C      nfilo = 1; nfiup = 1
C      nemlo = 2; nemup = 2

C ... Skip this initialization step if on-the-fly integration
      if (.not. lonekp) then
        if (jdos) then
          if (npol /= 1) call rx('OPTINT: npol should be 1 for jdos')
        else
          if (npol /= 3) call rx('OPTINT: npol should be 3 for eps')
        endif
        pi = 4*datan(1d0)
        volwgt = dble(3-nsp)/(n1*n2*n3*6)
        if (jdos) volwgt = 1d0/(n1*n2*n3*6)
        ds2 = dsqrt(2d0)*1d-6
        ds3 = dsqrt(3d0)*1d-6
        fac = 32d0*pi**2/vol
        if (nspc == 2) fac = fac / 2d0
        call dpzero(opt,nptmx*npol*nsp)
      endif

      if ((ikp == 1 .and. ksp == 1 .or. .not. lonekp) .and. ipr >= 20) then
      if (jdos) then
        call awrit8('%N OPTINT:  JDOS using occ bands (%i,%i) '//
     .    'and unocc (%i,%i)%N%10f%i points in energy window (%d,%d)'//
     .    '%?;n;%N%10fResolve partial (occ,unocc) band contributions;;',
     .    ' ',180,stdo,nfilo,nfiup,nemlo,nemup,npts,emin,emax,
     .    lpart0)
      else
        call info8(20,1,0,' OPTINT:  eps2(w) using occ bands (%i,%i) '//
     .    'and unocc (%i,%i)%N%10f%i points in energy window (%d,%d)'//
     .    '%?;n;%N%10fResolve partial (occ,unocc) band contributions;;',
     .    nfilo,nfiup,nemlo,nemup,npts,emin,emax,lpart0)
        if (loptic == 8 .or. loptic == 9) call info8(20,0,0,
     .    '%10fSimulate %?;(n==8);Absorption;Emission; with '//
     .    '%?;n;Ef=%d%j;Ef(p,n)=(%d,%d); at kT=%d Ry',
     .    loptic,isw(s_optic%imref(1) == s_optic%imref(2)),
     .    s_optic%imref(1),s_optic%imref(2),s_optic%kt,npts,emax,lpart0)
      endif
      if (efermi == -99d0 .and. lonekp) then
        call rx('optint: Fermi level not set .. rerun with METAL=2 or METAL=3')
      elseif (efermi == -99d0) then
        call rx('OPTINT: Fermi level not set ... something wrong')
      endif
      if (iabs(loptic) /= 8 .and. iabs(loptic) /= 9) then
      call info5(20,0,0,'%10fEf=%d  Treating bands as %?;n;metal;insulator w/ Escissor=%d;',efermi,isw(metal),esciss,0,0)
      endif
      if (lteto > 0) call info2(10,0,0,'%10fIntegration with tetrahedron method, mode %i',lteto,0)
      if (samp) call info5(10,0,0,'%10fM-P sampling integration'//
     .  '%?;(n>0); (on-the-fly);;:  polynomial order=%i, width=%d',
     .  ikp,mpsord,swidth,0,0)
      endif

      if (.not. samp .and. ntet == 0)
     .  call rx('OPTINT: tetrahedra missing.  Restart with metal and TETRA=t')

C --- Write out file for LXB code ---
      if (cmdopt('--lxb',5,0,strn)) then
        if (lonekp) call
     .    rx('optint: --lxb not allowed with on-the-fly integration')
        ifil = fopn('LXB')
        rewind ifil
        write (ifil,*) nkp,nfiup,efermi
        bonly = cmdopt('--bonly',7,0,strn)
        do  50  isp = 1, nsp
          do  40  iq = 1, nkp
            ibf = 0
            do  30  ib = nfilo, nfiup
              ibf = ibf+1
              write (ifil,*) eband(ib,isp,iq)
              if (bonly) goto 30
              ibm = 0
              do  20  jb = nemlo, nemup
                ibm = ibm+1
                do  10  k = 1, 3
                  write (ifil,*) optmt(k,ibf,ibm,isp,iq)
   10           continue
   20         continue
   30       continue
   40     continue
   50   continue
        call fclose(ifil)
      endif

C --- Check band limits ---
      if (cmdopt('--chklim',8,0,strn)) then
        if (lonekp) call rx('optint: --chklim not allowed with on-the-fly integration')
        del = 0d0
        if (samp) del = 2.5d0*swidth
        ifmin =  10*nbmax
        iwmin =  10*nbmax
        ifmax = -10*nbmax
        iwmax = -10*nbmax
        do  150  isp = 1, nsp
          do  140  iq = 1, nkp
            do  60  ib = 1, nbmax
              if (eband(ib,isp,iq) > efermi+emin-del) goto 70
   60       continue
   70       ifmin = min0(ifmin,ib-1)
            do  80  ib = 1, nbmax
              if (eband(ib,isp,iq) > efermi+del) goto 90
   80       continue
   90       ifmax = max0(ifmax,ib)
            do  100  ib = 1, nbmax
              if (eband(ib,isp,iq) > efermi-emax-del) goto 110
  100       continue
  110       iwmin = min0(iwmin,ib-1)
            do  120  ib = 1, nbmax
              if (eband(ib,isp,iq) > efermi+emax+del) goto 130
  120       continue
  130       iwmax = max0(iwmax,ib)
  140     continue
  150   continue
        ifmin = max0(ifmin,1)
        iwmin = max0(iwmin,1)
        if (iwmin < nfilo .and. nfilo > 1)
     .    write(stdo,500) nfilo,iwmin
        if (ifmax > nfiup .and. nfiup < nbmax)
     .    write(stdo,510) nfiup,ifmax
        if (ifmin < nemlo .and. nemlo > 1)
     .    write(stdo,520) nemlo,ifmin
        if (iwmax > nemup .and. nemup < nbmax)
     .    write(stdo,530) nemup,iwmax
  500   format(' ***WARNING***  lower occupied band index too high:'/
     .    18x,'NFILO=',i4,' > ',i4)
  510   format(' ***WARNING***  upper occupied band index too low:'/
     .    18x,'NFIUP=',i4,' < ',i4)
  520   format(' ***WARNING***  lower unoccupied band index too high:'/
     .    18x,'NEMLO=',i4,' > ',i4)
  530   format(' ***WARNING***  upper unoccupied band index too low:'/
     .    18x,'NEMUP=',i4,' < ',i4)
        write(stdo,540) iwmin,ifmax,nfilo,nfiup,ifmin,iwmax,nemlo,nemup
  540   format(/' Check band index limits:'/
     .    10x,'Calculated',12x,'Input'/
     .    ' Occ:    (',i4,', ',i4,')',8x,'(',i4,', ',i4,')'/
     .    ' Unocc:  (',i4,', ',i4,')',8x,'(',i4,', ',i4,')')
      endif

C --- Use bjnos tetrahedron BZ integration ---
      if (lteto == 1 .or. lteto == 2) then
        if (lpart0 /= 0) call rx2(
     .  'OPTICS_PART=%i not implemented for LTET=%i',lpart0,lteto)
        call bzjnos(n1,n2,n3,eband,nkp,nbmax,nsp,npol,nfilo,nfiup,
     .    nemlo,nemup,emin,emax,esciss,jdos,optmt,opt,npts,
     .    efermi,ntet,idtet)
        if (jdos .and. nsp == 1) then
          call dscal(npts,0.5d0,opt,1)
        endif
        call xxxdif(emin,emax,npts,npol*nsp,0,opt)

C --- Use sampling instead of tetrahedron BZ integration ---
      elseif (lteto == 0 .and. ikp == -1) then
        if (lpart0 /= 0) call rx2(
     .  'OPTICS_PART=%i not implmented for OPTICS_LTET=%i',lpart0,lteto)

C   ... Nonequilibrium sampling
        if (iabs(loptic) == 8 .or. iabs(loptic) == 9) then
          call nesjdos(loptic,nkp,nbmax,nfilm,nempm,nsp,nfilo,nfiup,nemlo,nemup,
     .      wgts,eband,mpsord,swidth,-srnge,emin,emax,esciss,jdos,optmt,
     .      s_optic%imref,s_optic%kt,npts,opt)
        else
          call mkjdos(nkp,nbmax,nfilm,nempm,nsp,nfilo,nfiup,nemlo,nemup,
     .      wgts,eband,mpsord,swidth,-srnge,emin,emax,esciss,jdos,optmt,
     .      efermi,npts,opt)
        endif
        if (jdos .and. nsp == 1) then
          call dscal(npts,0.5d0,opt,1)
        endif

C --- Partial accumulation of opt for 1 kpt, spin ---
      elseif (lteto == -1) then
        if (lpart0 /= 0) call rx2(
     .  'OPTICS_PART=%i not implmented for OPTICS_LTET=%i',lpart0,lteto)
        if (jdos .and. nsp == 2)
     .    call rx('OPTINT: fix indexing opt(:,:,isp)->opt(:,isp,:)')
        call mkjdos(1,nbmax,nfilm,nempm,1,nfilo,nfiup,nemlo,nemup,
     .    wgts(ikp)/nsp,eband(1,ksp,ikp),mpsord,swidth,-srnge,emin,emax,
     .    esciss,jdos,optmt,efermi,npts,opt(1,1,ksp))
       call tcx('optint')
       return

C --- Tetrahedron integration as originally implemented ---
      else

        call rx('optint: ltet=1 no longer supported; use ltet=2')

C        if (lpart0 > 1) call rx2(
C     .  'OPTICS_PART=%i not implmented for OPTICS_LTET=%i',lpart0,lteto)
C
C        allocate(yy(npts),optpar(nptmx,3))
C        if (npopt > 0) then
C          allocate(popt(nptmx,npol,nfilm,nempm))
C          call dpzero(popt,nptmx*nfilm*nempm)
C        endif
C
CC   ... Loop over spins
C        do  isp = 1, nsp
CC   ... Loop over (occ,unocc) pairs
C        ibf = 0
C        do  ib = nfilo, nfiup
C          ibf = ibf+1
C
C          ibm = 0
C          do  jb = nemlo, nemup
C            ibm = ibm+1
C            if (jb <= ib) cycle
C
CC       ... Partial ib --> jb transition analysis
C            call dpzero(optpar, npol*nptmx)
C
CC       ... Loop over tetrahedra
C            do  itet = 1, ntet
C              iq1 = idtet(1,itet)
C              iq2 = idtet(2,itet)
C              iq3 = idtet(3,itet)
C              iq4 = idtet(4,itet)
C
CC         ... Energies at 4 corners of tetrahedron for ib
C              eci(1) = eband(ib,isp,iq1) + ds2*0d0
C              eci(2) = eband(ib,isp,iq2) + ds2*1d0
C              eci(3) = eband(ib,isp,iq3) + ds2*2d0
C              eci(4) = eband(ib,isp,iq4) + ds2*3d0
C              if (dmin1(eci(1),eci(2),eci(3),eci(4)) > efermi) cycle
CC         ... Energies at 4 corners of tetrahedron for jb
C              ecj(1) = eband(jb,isp,iq1) + ds3*0d0 + esciss
C              ecj(2) = eband(jb,isp,iq2) + ds3*1d0 + esciss
C              ecj(3) = eband(jb,isp,iq3) + ds3*2d0 + esciss
C              ecj(4) = eband(jb,isp,iq4) + ds3*3d0 + esciss
C              if (dmax1(ecj(1),ecj(2),ecj(3),ecj(4)) < efermi+emin)
C     .          cycle
C
C              do  k = 1, npol
C                wt = volwgt*idtet(0,itet)
C                if (.not. jdos) then
C                  wtm =optmt(k,ibf,ibm,isp,iq1)+optmt(k,ibf,ibm,isp,iq2)
C     .               + optmt(k,ibf,ibm,isp,iq3)+optmt(k,ibf,ibm,isp,iq4)
C                  wt = wt*wtm/4
C                endif
C                call opt0(eci,ecj,emin,emax,npts,yy,efermi,wt)
C                do  ipts = 1, npts
C                  optpar(ipts,k) = optpar(ipts,k) + yy(ipts)
C                enddo
C              enddo
CC       ... End of loops over energy, polarizations, tetrahedra
C            enddo
C            if (jdos) then
C              call daxpy(nptmx*npol,1d0,optpar,1,opt(1,isp,1),1)
C            else
C              call daxpy(nptmx*npol,1d0,optpar,1,opt(1,1,isp),1)
C            endif
C
CC       ... Copy band-to-band contributions
C            if (lpart0 /= 0) then
C              do  ip = 1, npol
C                call dcopy(npts,optpar(1,ip),1,popt(1,ip,ibf,ibm),1)
C              enddo
C            endif
C          enddo
C        enddo
C        enddo
C        deallocate(yy,optpar)
      endif
C ... End of loops over final band, initial band, spins

C --- Write popt,opt to disk ---
      ip = 13
      if (loptic < 0) ip = 1
      call optin2(s_optic,s_bz,ip,nsp,nspc,vol,nptmx,npol,
     .  lpart1 /= 0,npopt,npts,emesh,popt,opt)

      deallocate(emesh)
      if (allocated(popt)) deallocate(popt)
      call tcx('optint')

      end

      subroutine opt0(e0,h0,emin,emax,numb,y,efermi,vol)
C- Contribution to matrix elements from one tetrahedron
C ----------------------------------------------------------------------
Ci Inputs
Ci   e0    :Energies at corners of tetrahedra, occupied states
Ci   h0    :Energies at corners of tetrahedra, unoccupied states
Ci   emin  :beginning of energy window
Ci   emax  :end of energy window
Ci   numb  :number of points in window
Ci   efermi:Fermi energy
Ci   vol   :weight for this tetrahedron (e.g. for joint DOS, volume)
Co Outputs
Co   y     :contribution to matrix elements added to y(1..numb)
Cr Remarks
Cr   This routine neglects energy dependence of numerator.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer numb
      double precision e0(4),h0(4),emin,emax,efermi,vol,y(numb)
C ... Local parameters
      integer i,j,jj,ind(4)
      double precision e(4),h(4),om0(4),om(4),a0(3),b0(3)
      double precision de,omg,dom,dom1,dom2,dom3,f,aa,bb
      double precision gb,ga,omg3,om23,emdl,hmdl,tmprr

      de = (emax-emin)/(numb-1)

C ... om0 = energy difference (unocc - occ) for each corner
      do  11  i = 1, 4
   11 om0(i) = h0(i)-e0(i)

C ... om = om0, sorted in reverse order; ditto for e,h
      do  10  i = 1, 4
   10 ind(i) = 1
      do  1  i = 1, 3
      do  1  j = i+1, 4
      if (om0(i) > om0(j)) goto 2
      ind(i) = ind(i)+1
      goto 1
    2 ind(j) = ind(j)+1
    1 continue
      do  20  i = 1, 4
      jj = ind(i)
      om(jj) = om0(i)
      e(jj) = e0(i)
      h(jj) = h0(i)
   20 continue

C ... For each omega in energy window, do
      do  100  i = 1, numb
      omg = emin + de*(i-1)

C ... Case energy smaller than smallest point of tetrahedron
      if (omg <= om(4)) then
        y(i) = 0d0

C ... Case energy between points 3 and 4
      elseif (omg <= om(3)) then

        dom = omg - om(4)
        dom1 = om(1) - om(4)
        dom2 = om(2) - om(4)
        dom3 = om(3) - om(4)
        f = dom**2/dom1/dom2/dom3

        aa = efermi
        a0(1) = e(4) + (e(1) - e(4))*dom/dom1
        a0(2) = e(4) + (e(2) - e(4))*dom/dom2
        a0(3) = e(4) + (e(3) - e(4))*dom/dom3

        bb = efermi
        b0(1) = h(4) + (h(1) - h(4))*dom/dom1
        b0(2) = h(4) + (h(2) - h(4))*dom/dom2
        b0(3) = h(4) + (h(3) - h(4))*dom/dom3

        call inta(a0,aa,ga)
        call inta(b0,bb,gb)

        y(i) = -3*vol*f*(gb-ga)

C ... Case energy between points 2 and 3
      elseif (omg <= om(2)) then

        omg3 = omg - om(3)
        om23 = om(2) - om(3)
        emdl = e(3) + (e(2)-e(3))*omg3/om23
        hmdl = h(3) + (h(2)-h(3))*omg3/om23

C       <I>
        dom = omg - om(4)
        dom1 = om(1) - om(4)
        dom2 = om(2) - om(4)
        f = dom/dom1/dom2 * (om(2)-omg)/om23

        aa = efermi
        a0(1) = e(4) + (e(1)-e(4))*dom/dom1
        a0(2) = e(4) + (e(2)-e(4))*dom/dom2
        a0(3) = emdl

        bb = efermi
        b0(1) = h(4) + (h(1)-h(4))*dom/dom1
        b0(2) = h(4) + (h(2)-h(4))*dom/dom2
        b0(3) = hmdl

        call inta(a0,aa,ga)
        call inta(b0,bb,gb)
        tmprr = -f*(gb-ga)

C       <II>
        dom = om(1) - omg
        dom1 = om(1) - om(4)
        dom2 = om(1) - om(3)
        f = dom/dom2/dom1 * (omg-om(3))/om23

        aa = efermi
        a0(1) = emdl
        a0(2) = e(1) - (e(1)-e(3))*dom/dom2
        a0(3) = e(1) - (e(1)-e(4))*dom/dom1

        bb = efermi
        b0(1) = hmdl
        b0(2) = h(1) - (h(1)-h(3))*dom/dom2
        b0(3) = h(1) - (h(1)-h(4))*dom/dom1

        call inta(a0,aa,ga)
        call inta(b0,bb,gb)

        y(i) = 3*vol*(tmprr-f*(gb-ga))

C ... Case energy between points 1 and 2
      elseif (omg <= om(1)) then

        dom = om(1) - omg
        dom1 = om(1) - om(2)
        dom2 = om(1) - om(3)
        dom3 = om(1) - om(4)
        f = dom**2/dom1/dom2/dom3

        aa = efermi
        a0(1) = e(1) - (e(1)-e(2))*dom/dom1
        a0(2) = e(1) - (e(1)-e(3))*dom/dom2
        a0(3) = e(1) - (e(1)-e(4))*dom/dom3

        bb = efermi
        b0(1) = h(1) - (h(1)-h(2))*dom/dom1
        b0(2) = h(1) - (h(1)-h(3))*dom/dom2
        b0(3) = h(1) - (h(1)-h(4))*dom/dom3

        call inta(a0,aa,ga)
        call inta(b0,bb,gb)

        y(i) = -3*vol*f*(gb-ga)

C ... Case energy above point 1
      else
        y(i) = 0d0
      endif

  100 continue

      end


      subroutine inta(a0,aa,g)
c
      double precision a(3),a0(3)
      double precision aa
      double precision daa,dal,dam,das
      double precision g
      integer i,j,jj
      integer ind(3)

      do  10  i = 1, 3
   10 ind(i) = 1

      do  1  i = 1, 2
      do  1  j = i+1, 3
      if (a0(i) > a0(j)) goto 2
      ind(i) = ind(i)+1
      goto 1
    2 ind(j) = ind(j)+1
    1 continue

      do  20  i = 1, 3
      jj = ind(i)
   20 a(jj) = a0(i)

      if (aa > a(3)) goto 4
c     type*,'case 0'
      g = 0d0
      return

    4 if (aa > a(2)) goto 5
c     print *,'case 1'
      daa = aa-a(3)
      dam = a(2)-a(3)
      dal = a(1)-a(3)
      g = (daa/dam)*(daa/dal)
      return

    5 if (aa > a(1)) goto 6
c     print *,'case 2'
      daa = a(1)-aa
      dam = a(1)-a(2)
      das = a(1)-a(3)
      g = (1d0-(daa/dam)*(daa/das))
      return

    6 g=1d0
c     print *,'case 3'
      return

      end

      subroutine optin2(s_optic,s_bz,lfac,nsp,nspc,vol,nptmx,
     .  npol,lbin,npopt,npts,emesh,popt,opt)
C- File output for Im(eps) or joint DOS, generated by optint or optinq
C ----------------------------------------------------------------------
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic window ocrng unrng esciss ffmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet lmet n range w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs:
Ci   lfac  :scaling factor
Ci         :1s digit controls constant factor
Ci         :1 scale opt by 1
Ci         :2 scale opt by 2
Ci         :3 scale opt by 32d0*pi**2/vol  (Im eps)
Ci         :10s digit:
Ci         :1 scale opt by additional factor 1/E^2 (Im eps)
Ci         :2 print data including zero and negative frequencies
Ci   nsp   :number of spins if nspc=1, 1 if nspc=2
Ci   nspc  :2 for noncoll or S/O; otherwise 1
Ci   vol   :volume
Ci   nptmx :leading dimension of opt,popt
Ci   npol  :number of polarizations; needed for dimensioning.
Ci         :npol must be 1 for jdos, 3 for eps.
Ci   lbin  :0 save popt as ASCII file
Cl         :1 save popt as binary file
Ci   npopt :number of (ib,jb) pairs eps or jdos is resolved into
Ci         :If 0, no partial decomposition is written
Ci   npts  :number of energy points
Ci   emesh :energy mesh
Ci   popt  :Im(eps) or joint DOS, resolved by (occ,unocc) or k, or both
Ci   opt   :Im(eps) or joint DOS
Co Outputs:
Co   opt:  :written to disk
Co   popt  :written to disk if npopt>0
Cr Remarks
Cr   Sanity checks:
Cr   mc -qr popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -csum -ccat opt.ogan -e2 x1 x2+x3+x4 --
Cr   mc popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -coll 1:84:3 -csum -ccat opt.ogan -e2 x1 x2 --
Cr   or:
Cr   mc popt.ogan -split a 1,nr+1 1,2,nc+1 a11 a12 -csum -ccat jdos.ogan --
Cu Updates
Cu   27 Feb 15 s_optic%ffmt allows writing of opt file in other formats
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Sep 09 Rearrangement to accomodate tetwtq.  Now writes popt.
Cu   15 Jan 04 (Uk-Jin Roh) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lbin
      integer lfac
      integer nsp,nspc,npts,nptmx,npol,npopt
      double precision vol,emesh(npts)
      double precision opt(nptmx,npol*nsp),popt(nptmx,npol,npopt)
C ... For structures
!      include 'structures.h'
      type(str_optic):: s_optic
      type(str_bz)::    s_bz
C ... Dynamically allocated arrays
      real(8),allocatable :: poptk(:,:)
      real(8),allocatable :: pfile(:,:)
C ... Local parameters
      logical jdos
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet,ipts,i,j,ip,lmet,ic
      integer ifio,ifiop,fopn,fopna,loptic,nchan
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3
      integer mpsord
      double precision swidth,srnge
      double precision wt,ei,optrng(4),emin,emax,fac,pi
      double precision esciss
C     character*120 strn
      equivalence (emin,optrng(1)),(emax,optrng(2))
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))

C --- Setup ---
      loptic = s_optic%loptic
      optrng = s_optic%window
      ocrng = s_optic%ocrng(1:2)
      unrng = s_optic%unrng(1:2)

      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
      lmet = s_bz%lmet
      esciss = s_optic%esciss
      mpsord = s_bz%n
      srnge = s_bz%range
      swidth = s_bz%w
      mpsord = mod(mpsord,100)
      jdos = loptic < 0
      if (jdos) then
C       if (npol /= 1) call rx('OPTIN2: npol should be 1 for jdos')
        ifio = fopn('JDOS')
        rewind ifio
        call info0(30,0,0,' OPTINT:  writing jdos file ...')
      else
        if (npol /= 3) call rx('OPTIN2: npol should be 3 for eps')
        ifio = fopn('OPT')
        rewind ifio
        call info0(30,0,0,' OPTINT:  writing opt file ...')
      endif

C --- Write opt to disk ---
      pi = 4*datan(1d0)
      fac = 1
      if (mod(lfac,10) == 2) fac = 2
      if (mod(lfac,10) == 3) fac = 32d0*pi**2/vol
      if (nspc == 2) fac = fac / 2d0
      nchan = npol*nsp
      write(ifio,'(''% cols '',i2,''  spin'',i2)') 1+npol*nsp,nsp
      do  ipts = 1, npts
        ei = emesh(ipts)
        if (ei > 1d-5 .or. mod(lfac/10,10) >= 2) then
          wt = fac
          if (mod(lfac/10,10) == 1) wt = wt / ei**2
          if (s_optic%ffmt == 0) then
C            write (ifio,550) ei,(wt*opt(ipts,ip,1),ip=1,nchan)
C  550       format(7f13.6)
            call awrit3('%;13,6D%n;13,6D',' ',120,ifio,ei,nchan,wt*opt(ipts,1:nchan))
          elseif (s_optic%ffmt == 1) then
            write (ifio,551) ei,(wt*opt(ipts,ip),ip=1,nchan)
  551       format(1p7e14.6)
          endif
        endif
      enddo
      call fclose(ifio)

C --- Write popt to disk ---
      if (npopt /= 0) then
      i = npts
      if (emin < 1d-5 .and. mod(lfac/10,10) < 2) i = npts-1
      if (lbin) then
        call info2(30,0,0,'          writing binary file poptb:  '//
     .  '%i channels%?;(n>1);%-1j, %i polarizations;;',npopt,npol)
        ifiop = fopna('poptb',-1,4)
        rewind ifiop
        allocate(pfile(i,1+npol*npopt))
      else
        call info2(30,0,0,'          writing ascii file popt:  '//
     .  '%i channels%?;(n>1);%-1j, %i polarizations;;',npopt,npol)
        ifiop = fopn('popt')
        rewind ifiop
        call awrit2('%% rows %i cols %i',' ',80,ifiop,i,1+npol*npopt)
      endif
      i = 0
      allocate(poptk(npol,npopt))
      do  ipts = 1, npts
C       ei = emin + (ipts-1)*de
        ei = emesh(ipts)
        if (ei > 1d-5 .or. mod(lfac/10,10) >= 2) then
          wt = fac
          if (mod(lfac/10,10) == 1) wt = wt / ei**2
          if (lbin) then
            i = i+1
            pfile(i,1) = ei
            do  ip = 1, npol
            do  j = 1, npopt
              ic = 2 + (j-1) + npopt*(ip-1)
              pfile(i,ic) = wt*popt(ipts,ip,j)
            enddo
            enddo
          else
            poptk(:,:) = wt*popt(ipts,:,:)
            call zersmalln(npol*npopt,5d-7,poptk)
            write(ifiop,334) ei, ((poptk(ip,j),j=1,npopt),ip=1,npol)
  334       format(f12.6,6f13.6/(12x,6f13.6))
          endif
        endif
      enddo
      deallocate(poptk)

      if (lbin) then
        call ywrm(1,' ',1,ifiop,' ',pfile,0,i,i,1+npol*npopt)
        deallocate(pfile)
      endif

      call fclose(ifiop)
      endif

      end
