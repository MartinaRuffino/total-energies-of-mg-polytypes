      subroutine iomomx(ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,
     .  nband,nchan,nchan2,nchds,nevmx,nlistc,listc,eb,accwt,doswt,
     .  dosw2,efermi,vmtz)
C- Read selected channels of data from moments file
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file handle
Ci   mode  : 1's digit 0 read nl,nsp,nspc,nkp,ldim,nfstg
Ci         :           1 require a match in nl,nsp; read nspc,ldim,nkp
Ci         :           2 require a match in nl,nsp,nspc,nkp,ldim
Ci         :10's digit 0 exit after reading header info
Ci         :       1,2,3 read number of iq before EOF encountered
Ci         :         2,3 read all information sought by nfstg
Ci         :           3 read efermi,vmtz if available
Ci   nband :leading dimension of eb,dosw2,doswt,accts
Ci   nchds :number of DOS channels to read
Ci   nlistc:number of classes, species or sites (see Remarks)
Ci   listc :class, species or site list
Ci   nchan :number of channels (l+m+class) for accwt,doswt
Ci   nchan2:number of channels for dosw2
Cio Inputs/Outputs
Cio    ... The following are read for 1s digit mode=0; else they are input
Cio  nl    :(global maximum l) + 1
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Cio  nkp   :number of irreducible k-points (bzmesh.f)
Cio  ldim  :dimension of hamiltonian matrix (makidx.f)
Cio  nfstg :describes information contained in moments file (see iomoms).
Cio        :iomomq passes nfstg to iomoms; moments file (ifi) may contain
Cio        :more information than sought by nfstg; but if some information
Cio        :is missing and 10's digit of mode is 2 or 3, iomomq aborts.
Cio        :For 1s digit mode > 0, nfstg is an input.
Co Outputs:
Co   iq     :(1s digit mode>0) number of qp available in moments file.
Co   nevmx  :largest # evals encountered
Cf  Files:
Cf    Moments file has the following records
Cf    1.   nl nsp nspc nkp ldim nfstg
Cf    ... For each qpt (and spin), the following records:
Cf    2.   nchan  nev (if nfstg nonzero)
Cf         eband      (if 1s   digit of nfstg nonzero)
Cf         accwt      (if 10s  digit of nfstg 1)
Cf         doswt      (if 10s  digit of nfstg 2)
Cf         dosw2      (if 100s digit of nfstg 1)
Cf    3.   efermi, vmtz
Cb Bugs
Cb   dosw2 shouldn't be dimensioned nchds!
Cr Remarks
Cr   Same as iomomq, but reads a subset of DOS channels.
Cr   The dos channels are assumed be ordered by class, species, or site
Cr   (it doesn't matter for this routine's purposes) with nl channels
Cr   per class, species, or site.
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu    7 Apr 04 dosw2 can have different nchan than doswt
Cu   27 Apr 02 Bug fix
Cu   18 Jan 02 Cleaner conventions involving 10s digit nfstg.
Cu             Old accwt and doswt merged into one array.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,nband,nchan,nchan2,
     .  nchds,nevmx,nlistc,listc(nlistc)
      double precision eb(nband*nsp,1),doswt(nchds,nband*nspc,nspc,1),
     .  accwt(nchds,nband*nspc,nspc,3,1),dosw2(2*nchds*nband*nspc*3,1),
     .  efermi,vmtz
C ... Dynamically allocated local arrays
      real(8), allocatable :: ebl(:)
      real(8), allocatable :: accm(:)
      real(8), allocatable :: dosw(:)
C ... Local parameters
      double precision xx
      integer nlf,nspf,nspcf,nkpf,ldimf,nfstgf,iprint,i1mach,jq,iomoms,
     .  nev,nschan
      character outs*80

C --- Read header, checking req'd matches and copying the rest ---
      rewind ifi
      if (mod(mode,10) == 0) then
        read (ifi,err=999,end=999) nl, nsp, nspc, nkp, ldim, nfstg
        if (mode == 0) return
      else
        read (ifi,err=999,end=999) nlf, nspf, nspcf, nkpf, ldimf, nfstgf
        if (mod(mode,10) == 1) then
          nkp  = nkpf
          nspc = nspcf
          ldim = ldimf
        endif
      endif
      iq = 0
      rewind ifi
      if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,1,iq,1,nevmx,nevmx,
     .  nchan,nchan2,nevmx,eb,accwt,doswt,dosw2,efermi,vmtz) < 0)
     .  goto 999
      if (mod(mode/10,10) == 0) return

C --- Determine number of qp available ---
      iq = 0
      do  jq = 1, nkp*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,0,1,jq,1,nband*nspc,
     .    nband*nspc,nchan,nchan2,nev,eb,accwt,doswt,dosw2,efermi,vmtz)
     . < 0) goto 14
        iq = jq-1
      enddo
C ... If read all qp, and another record also present, let iq=nkp
      read (ifi,err=14,end=14) xx
      iq = iq+1
   14 continue
      iq = (nspc*iq)/nsp
      call awrit1(' IOMOMX: read %i qp',outs,80,0,iq)

C --- Read info spec'd by nfstg for each qp until error ---
      if (mod(mode/10,10) /= 1) then
      nfstgf = nfstg
      nschan = mod(nfstg/10,10)
      rewind ifi
      read (ifi,err=999,end=999) nlf
      nevmx = 0
      allocate(ebl(nband))
      allocate(accm(nchan*nband*nfstg))
      allocate(dosw(nchan*nband*nfstg))
      if (nfstg >= 100) call rxi('iomomx not ready for nfstg=',nfstg)
      do  jq = 1, iq*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,nschan,1,1,nband*nspc,
     .    nband*nspc,nchan,nchan2,nev,ebl,accm,dosw,dosw2,
     .    efermi,vmtz) < 0)goto 999
        call pviomx(nchan,nchds,nband*nspc,nl,nspc,listc,nlistc,nfstg,
     .    nschan,jq,ebl,dosw,eb,doswt)
        nevmx = max(nev,nevmx)
      enddo
      deallocate(ebl,accm,dosw)

        if (mod(mode/10,10) /= 2 .and. iq == nkp) then

C --- Read efermi, vmtz if sought and info available ---
      read (ifi,err=999,end=999) efermi,vmtz
      call awrit2('%a  efermi=%,6d  vmtz=%,6d',outs,80,0,efermi,vmtz)
      endif
      endif

C --- Exit ---
      if (iprint() < 30) return
      call awrit0('%a from moms file',outs,80,-i1mach(2))
      return

C --- Error exit ---
  999 continue
      if (iprint() > 10)
     .  print *, 'IOMOMX (warning): empty or bad moments file'
      iq = -1
      end

      subroutine pviomx(nchan,nchds,ndim,nl,nspc,listc,nlistc,nfstg,
     .  ldwt3,iq,eb0,dosw0,eb,doswt)
C- Copy a subset of channels
      implicit none
      integer nchan,nchds,ndim,nspc,nl,ldwt3,nfstg,iq,nlistc,
     .  listc(nlistc)
      double precision eb0(ndim),dosw0(nchan,ndim,nspc)
      double precision eb(ndim,iq),doswt(nchds,ndim,ldwt3,iq)
      integer il,ic,ichan,ichds,imax,i
      if (nspc /= 1) call rx('pviomx implemented only for nspc=1')

      il = 1
      ichan = -nl
      ichds = -nl
      imax = listc(nlistc)
      if (mod(nfstg,10) == 1) call dcopy(ndim,eb0,1,eb(1,iq),1)
      do  ic = 1, imax
        ichan = ichan + nl
        if (ic < listc(il)) cycle
        ichds = ichds + nl
        if (mod(nfstg/10,10) >= 1) then
          do  i = 1, nl
            call dcopy(ndim*nspc,dosw0(ichan+i,1,1),nchan,
     .        doswt(ichds+i,1,1,iq),nchds)
          enddo
        endif

        il = il+1
      enddo
      end
