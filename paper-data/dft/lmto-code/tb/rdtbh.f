      subroutine rdtbh(ifi,nclass,dclabl,z,nlmesp,ltb,nterm,ipass,
     .  fitpar,rmaxh,iam,npm,tabme,tabcf,tabov,tbocf,
     .  decay,deccf,decov,dcocf,cut,cutov,itab,itbcf,itbov,itocf,
     .  idec,idcf,idov,idocf,v0,k,memode,ppmode,poly,cutmod,cutpp)
C- Disk-read of tight-binding hamiltonian and pair potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   : file handle for CTRL
Ci   nlmesp: (number of ME/spin channel) * (# spin channels)
Ci           for which ME are input
Ci   nclass,nsp,ltb
Ci   nterm : max number of parameters for each matrix element,
Ci           leading dimension of itab,itbcf,itbov, and itocf
Ci   ipass : 1, read tabme; 2, generate iam, v0
Ci   fitpar: if true then fitting parameters, read freeze switches
Ci   rmaxh : Hamiltonian cutoff (used to set default PP and ME cutoffs)
Ci           rmaxh passed into rdtbh in the same units as PP and ME
Ci           cutoffs (ie either in Bohr or in units of alat)
Co Outputs
Co   iam(1..3,kk): are a group of three numbers that associate a pair
Co           with a set of matrix elements.  For a pair of classes
Co           iam(1,kk) and iam(2,kk), iam(3,kk) is an index specifying which
Co           rule in tabme (and possibly tabcf,tabov,tbocf) that defines the
Co           matrix elements between these classes.  It is possible that
Co           there is no rule connecting two particular classes, in which
Co           case the matrix element is taken to be zero.
Co           iam is ordered according to increasing class.
Co   npm   : is a table of offsets to iam so that every rule associated
Co           with a given class can be easily identified.  npm(0,i) is
Co           the number of entries in iam containing class i, and
Co           npm(1,i) is the offset to iam for the first entry.
Co   tabme : a set of parameters that correspond to coefficients of
Co           Slater-Koster, or Slater-Koster-like matrix elements.
Co           The meaning of the coefficients, depends on mode;
Co           see Remarks, and memode, below.  Which pair of atoms
Co           the coefficients apply is defined by a set of rules,
Co           as described in Remarks.
Co           tabme is a list of parameter, for as many rules
Co           are are read in.
Co           When ipass=1, the number of rules is returned (in k)
Co           and the corresponding coeficients are read in tabme
Co   tabcf : table of crystal field MEs [see D.J. Chadi in ``Atomistic
Co           Simulation of Materials Beyond Pair Potentials'', edited
Co           by V. Vitek and D. Srolovitz (Plenum, 1989), page 309].
Co           tabcf has the same structure as tabme
Co   tabov : table of overlap matrix elements, structured as tabme
Co   tbocf : table of crystal field MEs for overlap, structured as tabme
Co   decay,deccf,decov,dcocf: exponential or power decay parameters
Co           matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Co   cut,cutov: For each ME a pair of distances (r1,rc). For some forms
Co           of decay the ME is cut-off smoothly to zero between r1 and rc
Co           (see tbham, also makv0)
Co   itab,itbcf,itbov,itocf,idec,idcf,idov,idocf: switches for all the
Co           matrix element and decay parameters: if 0 vary param., if 1 fix;
Co           only read if fitpar=.true.
Co   v0    : parameters for pair potential
Co           Ordering: a_1 b_1 c_1 a_2 ... c_3
Co           V0 = \sum_i=1,3 a_i d^b_i exp(-c_i d) unless b_1 > 0 and then
Co           V0 = a_1 eps + a_2 eps^2 where eps = (d - c_1) / c_1
Co           If b_1 > 0 and c_1 < 0 then Goodwin-Skinner-Pettifor
Co           V0 = A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)]
Co             ordering: A (b_1 > 0) (c_1 < 0) m mc r0 rc
Co   k     : (pass 1) number of ME rules
Co           (pass 2) number of ME pairs
Co   memode: 0, fixed MEs
Co           1, Harrison universal MEs
Co           2, exponential decay
Co           3, power decay
Co           4, ME = \sum_i=1,3 a_i d^b_i exp(-c_i d), the ordering is:
Co              a_1 b_1 c_1 a_2 ... c_3 for ss-sigma, then sp-sigma, etc
Co           5, Goodwin-Skinner-Pettifor,
Co              v_ij (r0/d)^n exp[n (-{d/rc}^nc + {r0/rc}^nc)]
Co              ordering: v_ij n nc r0 rc for ss-sigma, etc
Co           6, (nl=1) Slater-Koster + Srini's extra term
Co              NB: NOT implemented for spin pol.
Co           7, a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc)
Co              ordering: a b c d0 for ss-sigma, etc
Co              nterm=1 for memode=0-3, nterm=9 for memode=4,
Co              nterm=5 for memode=5, nterm=2 for memode=6,
Co              nterm=4 for memode=7
Co              memode >= 10, use canonical TB Hamiltonian
Co   ppmode: type of repulsive pair potential
Co           10s digit:
Co           1, sum of exponents times power law decay
Co              V0(d) = \sum_i=1,3 a_i d^b_i exp(-c_i d)
Co           2, quadratic
Co              V0(d) = a_1 eps + a_2 eps^2 where eps = (d - c_1) / c_1
Co           3,  Goodwin-Skinner-Pettifor
Co              V0(d) = A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)]
Co           1s digit defines particular adjustments to above types (see makvpp.f)
Co   poly  : degree of polynomial Pn employed to smoothly cut off V0 to zero
Co           4  cutoff tail matches value, slope, and curvature of V0(r) at r=r1
Co              whereas only value and slope are zero at r=rc
Co           5  (default) same as above except that value, slope, and curvature
Co              of Pn(r1) should all turn to zero
Co   cutmod: cutoff mode for pair potentials, hopping integrals, and overlap matrix
Co           (if applicable)
Co           0  no cutoff
Co           1  augmentative cutoff: V0(r) is augmented with a polynomial P_n(r)
Co              at r = r1 up to the second derivatives, whereas P_n(rc),
Co              P'_n(rc) and P"_n(rc) (if poly = 5) are all zeros
Co           2  multiplicative cutoff (default mode): V0(r) is multiplied by
Co              a polynomial P_n(r): P_n(r1)=1 and P_n(rc)=P'_n(rc)=P'_n(r1)=P"_n(r1)=0.
Co              If poly=5, also P"_n(rc)=0
Co   cutpp : pair potential cutoffs cutpp(1,*) = r1 and cutpp(2,*) = rc, see above.
Co           (previously contained within v0)
Co   ppmode, poly, cutmod, and cutpp are all arrays defined for each pair of species
Cr Remarks
Cr *rdtbh parses category ME, looking for a set of rules defining the
Cr  matrix elements.  For each rule, rdtbh reads a vector of numbers;
Cr  the precise meaning of these numbers depends on memode.  But
Cr  typically, the numbers are coefficients to Slater-Koster matrix
Cr  elements, whose ordering is:
Cr     ss-sigma, sp-sigma, pp-sigma, pp-pi, sd-sigma, pd-sigma,
Cr     pd-pi, dd-sigma, dd-pi, dd-delta
Cr
Cr *For each rule, rdtbh reads in two lists of classes (precisely how
Cr  the lists are defined is described below), followed by a vector
Cr  of numbers that define the matrix elements (eg ss-sigma, etc), then
Cr  optionally followed by another vector of numbers for the overlap,
Cr  crystal-field terms, etc, and finally another vector of numbers
Cr  for the pairwise potential.
Cr
Cr *A rule takes the following schematic form:
Cr     i-class-list j-class-list | me ! V0
Cr  If crystal field terms and/or overlap the form is
Cr     i-class-list j-class-list | tb me & cf me @ ovl me % ovl cf ! V0
Cr                          or
Cr     i-class-list j-class-list | tb me @ ovl me & cf me ! V0
Cr                          or etc.
Cr
Cr  where me (or tb me, cf me, ovl me, V0) are a vector of numbers.  How
Cr  many numbers there are depends on nl and memode, except for the
Cr  pair potential vector, which is always of length 9.
Cr
Cr *In the spin-orbit case all the (+,+) MEs for a given atom pair are
Cr  listed first, followed by the (-,-) MEs and then the (+,-) MEs
Cr
Cr *Syntax of class-list (i- and j-class-list): This string specifies
Cr  which groups of classes are to be associated with this rule.
Cr  It can be a simple number, eg for memode=0 and nl=2, a rule can be:
Cr     1 2 | vsss vsps vpps vppp
Cr  More generally,  the i- and j- class-list can be a list of integers
Cr  such as 1,2,3 or 1:5; see mkilst for the syntax of an integer list.
Cr
Cr  It is also possible to define the list in one of two other alternate
Cr  styles.   The second style is to define the list according to
Cr  an expression, that is i- or j-class-list is some expression
Cr  the variables ic and z; thus class-list may look like ic<6&z==14.
Cr  Every class for which the expression evaluates to true belongs in
Cr  the list.  Thus the rule
Cr     ic<6 ic>=6 | vsss vsps vpps vppp
Cr  defines the rule that will connect every class up to the fifth
Cr  class with every class beyond the fifth class.
Cr
Cr  The third alternative is specifically for unix systems.  i- and j-
Cr  class-list are filenames with the usual unix wildcards, eg
Cr  a[1-6].  Each class has an associated file; rdtbh invokes a
Cr  shell 'ls class-list | grep classname'.  Any class which ls finds
Cr  is added to the list.  The following example
Cr     a[1-6] b* | vsss vsps vpps vppp
Cr  creates a rule connecting any class a1, a2, a3, ...  a6 to any class
Cr  whose name begins with 'b'
Cr
Cr *You tell rdtbh to use one of these styles at the start of
Cr  the ME category, directly after memode, eg
Cr    ME  3 CLSTYL=2
Cb Bugs
Cb    In some cases if the pair potential line (!) is missing then rdtbh
Cb    fails. Specifically, at line 296 iend is not properly shifted in
Cb    which case at line 337 the overlap parameters are read again
Cb    into a rule k+1 for which there is no space to store them.
Cb    The workaround is to always include a pair potential line even if
Cb    it is a series of zeros.
Cu Updates
Cu   19 Apr 11 (SL)  Redesigned due to species-dependent ME and new cutoff
Cu                   options
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
Cu   21 Jul 02 (MvS) Bug fix in partok call
Cu   15 Feb 02 (ATP) Added MPI parallelization
C ----------------------------------------------------------------------
      use iolib_cb
      implicit none
C Passed parameters
      integer ifi,nclass,nlmesp,nterm,ipass,k,ltb
      integer memode(k),ppmode(k),poly(k),cutmod(k)
      integer procid, master, mpipid
      double precision dclabl(nclass)
      integer iam(3,k),npm(0:1,k),itab(nterm,nlmesp,k),
     .  itbcf(nterm,nlmesp,k),itbov(nterm,nlmesp,k),
     .  itocf(nterm,nlmesp,k),idec(nlmesp,k),idcf(nlmesp,k),
     .  idov(nlmesp,k),idocf(nlmesp,k)
      double precision tabme(nterm,nlmesp,k),tabcf(nterm,nlmesp,k),
     .  tabov(nterm,nlmesp,k),tbocf(nterm,nlmesp,k),decay(nlmesp,k),
     .  deccf(nlmesp,k),decov(nlmesp,k),dcocf(nlmesp,k),rmaxh,
     .  cut(2,nlmesp,k),cutov(2,nlmesp,k),v0(9,k),cutpp(2,k),z(nclass)
      logical fitpar
! C For iolib
!       integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .        iend,ichoos,nchoos,optio
!       logical noerr
!       common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
!      .               iend,ichoos,nchoos,noerr,optio
C Local parameters
      integer memx,i,j,inext,owk,partok,iprint,is,js,is0,is1,nt,
     .  memod0,ival,ii,jj,jjj,kk,recln0,recln
      character(len=8) :: clabl, clabls(nclass), blabl(2)
      parameter (recln0=120)
      integer ij(2),clstyl,iv0,ic,fextg
      double precision decay0,v00(9)
      logical scat,lsequ,cryf,ovl,ocryf,a2bin,sw
      character*1 ext*40
c ... ME modes
      integer, parameter :: memodx = 7      ! max ME mode number
      character*24 mode(0:memodx)
      integer ntm(0:memodx),ntm0,subsix,nset
c        memode 0 1 2 3 4 5 6 7 - - -
      data ntm /1,1,1,1,9,5,2,4/
      parameter (memx=500)
      character*(recln0) strn
      character*(recln0*memx) mfile
C       character(len=1) :: mfile(memx*recln0)
      integer iwk(nterm*nlmesp)
      double precision wk(nterm*nlmesp)
      integer i1,j1,k1
C ... for rdfiln
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,ctlen=120)
      character vnam(mxlev)*16,ctbl(mxchr,2)*(ctlen),outs*72
      logical loop0(0:mxlev),mlog,cmdopt,bittst
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
      procedure(integer) :: chaidx
      procedure(logical) :: isnum, c_isnum
      integer, allocatable :: wka(:)
      data mode /'Fixed','Universal ME','Exp. decay','Power decay',
     .  'a d^b exp(-c d)', 'Goodwin-Skinner-Pettifor',
     .  'Fixed+extension','a d^-b/{1+exp[c(d-d0)]}'/

C --- Initialization ---
      i = recln(recln0)
      optio = 1
      recoff = 0
      decay0 = 0
      clstyl = 0
      ovl   = bittst(ltb,1)
      cryf  = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)
      if (ipass == 2) nset = k

C --- Read in ME category from CTRL file ---
      if (procid == master) then
        if (.not. scat(ifi,'ME ',' ',.true.))
     .  call rx('RDTBH: missing category ME')
        backspace ifi

        nrecs = 0
        do  i = 1, memx*recln0, recln0
          mfile(i:i+recln0-1) = 'z'
          call rdfiln(ifi,'#{}% c',mxlev,loop0,nlin,list,lstsiz,
     .                ilist,nlist,vnam,ctbl,mxchr,strn,
     &                mfile(i:i+recln0-1),recln0,nrecs)
          if (i > 1 .and. mfile(i:i) /= ' ') goto 20
        enddo
        call rx('RDTBH: increase memx')

   20   continue

        nrecs = nrecs - 1
      endif
      call mpibc1(nrecs,1,2,mlog,'rdtbh','nrecs')
      call mpibcc(mfile,memx*recln0,mlog,'rdtbh','ctrl-ME')
      call getcat(mfile,'ME ',' ',.true.)
c ... save original length of mfile in subsix
      subsix = subsiz

C --- Determine ME mode, clstyl ---
C ... attempt to read the global ME mode
      memod0 = 999
      j = partok(mfile,'MEMOD0=','=',memod0,' ',-1,2,0,0)
C ... MEMOD0 token doesn't exist, let's try for just a number (old style)
      if (j == 0)
     .j = partok(mfile,'ME ',' ',memod0,' ',-1,2,0,0)
      if (j /= 0 .and. ipass == 2)
     .  call info2(20,1,0,' RDTBH: found global matrix elements'//
     .    ' with mode %i: '//trim(mode(memod0)),memod0,0)

C ... Pick up optional clstyl
      is = iend
      call skipbl(mfile,memx*recln0,is)
      js = is
      is = is+1
      call skp2bl(mfile,memx*recln0,js)
      if (mfile(is:is+6) == 'CLSTYL=') then
        j = partok(mfile,'CLSTYL=','=',clstyl,' ',-1,2,0,1)
      endif

c     if (iprint() >= 20 .and. ipass == 2 .and. memod0 >= 0) print
c    .  '(/'' RDTBH:  matrix elements with mode: '', a)', mode(memod0)

C --- Attempt to read in global exponential or power decay parameter ---
      inext = iend
      j = partok(mfile,'DECAY0=','=',decay0,' ',-1,4,0,0)
c ... move the start position in mfile to after
      if (j /= 0) inext = max(iend,inext)

      if (ipass == 2) then
        call dcopy(k*nlmesp,decay0,0,decay,1)
        if (cryf)  call dcopy(k*nlmesp,decay0,0,deccf,1)
        if (ovl)   call dcopy(k*nlmesp,decay0,0,decov,1)
        if (ocryf) call dcopy(k*nlmesp,decay0,0,dcocf,1)
        if (fitpar) then
          call iinit(idec,k*nlmesp)
          if (cryf)  call iinit(idcf, k*nlmesp)
          if (ovl)   call iinit(idov, k*nlmesp)
          if (ocryf) call iinit(idocf,k*nlmesp)
        endif
      endif

      k = 0
      kk = 0
      if (ipass == 2) then
C ... set defaults
        poly(1:nset) = 5
        ppmode(1:nset) = 0
        cutmod(1:nset) = 0
C ...   defaults for [r1, rc] are rmaxh*[1, 1.1]
        cutpp(1,1:nset) = rmaxh
        cutpp(2,1:nset) = rmaxh*1.1d0
        cut(1,1:nlmesp,1:nset) = rmaxh
        cut(2,1:nlmesp,1:nset) = rmaxh*1.1d0
        if (ovl) then
          cutov(1,1:nlmesp,1:nset) = rmaxh
          cutov(2,1:nlmesp,1:nset) = rmaxh*1.1d0
        endif
      endif

      do ic = 1, nclass
         call r8tos8(dclabl(ic),clabls(ic))
      end do

C --- Loop over input lines ---
      do 60  i = 1, nrecs
        k = k+1
        is = inext
        jj = 0

C --- Read v0 parameters ---
        ii = 1
        if (ipass == 2) ii = k

        jj = partok(mfile,'!','!',v00,' ',-9,4,is,0)
!         print *, 'k,ii,v01ii', k,ii,v00, nrecs, jj
        if (jj /= 0) then
          v0(1:9,ii) = v00
          inext = iend
          ntm0 = 1
C ...     read ME modes as optional parameters between is and iend
          subsiz = iend-1
          memode(k) = memod0
          jjj = partok(mfile,'MEMODE=','=',memode(k),' ',-1,2,is,0)

          if (memode(k) /= 999) then
            ntm0 = ntm(memode(k))
!             print *, 'memode', k, memode(k), ntm0
          else
C ...     memode should be set at least in one way
            call rxi(' rdtbh: either MEMODE or MEMOD0 should be'//
     .               'specified explicitly for pair = %i',k)
          endif
          if (ipass == 2) then
C ...     read other switches and pp cutoffs as optional parameters
C         between is and iend
            jjj = partok(mfile,'PPMODE=','=',ppmode(ii),' ',-1,2,is,0)
            jjj = partok(mfile,'POLY=','=',poly(ii),' ',-1,2,is,0)
            jjj = partok(mfile,'CUTMOD=','=',cutmod(ii),' ',-1,2,is,0)
C ...       as soon as CUTMOD is on, CUTPP becomes a compulsory parameter
            if (cutmod(ii) /= 0)
     .      jjj = partok(mfile,'CUTPP=','=',cutpp(1,ii),' ',-2,4,is,1)
          endif
          subsiz = subsix
        endif

        if (ipass == 2) then
C ---   Read crystal field ME and decay parameters ---
          if (cryf) then
            jjj = partok(mfile,'&','&',wk,' ',-ntm0*nlmesp,4,is,0)
            if (jjj /= 0) then
C ...         need to fill up array tabcf explicitly as ntm0 does not
C             necessarily match its leading dimension
C             (same for some other arrays below)
              do k1 = 1, ntm0*nlmesp
                j1 = (k1-1)/ntm0 + 1
                i1 = k1 - (j1-1)*ntm0
                tabcf(i1,j1,k) = wk(k1)
              enddo
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',deccf(1,k),' ',
     .          -nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZ=',4,'=',is0)) then
c    .            jjj = partok(mfile,'FRZ=','=',itbcf(1,k),' ',
c    .            -ntm0*nlmesp,2,is1,0)
                  jjj = partok(mfile,'FRZ=','=',iwk,' ',
     .            -ntm0*nlmesp,2,is1,0)
                  if (jjj /= 0) then
                    do k1 = 1, ntm0*nlmesp
                      j1 = (k1-1)/ntm0 + 1
                      i1 = k1 - (j1-1)*ntm0
                      itbcf(i1,j1,k) = iwk(k1)
                    enddo
                  endif
                endif
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idcf(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif

C --- Read overlap ME and decay parameters ---
          if (ovl) then
c           jjj = partok(mfile,'@','@',tabov(1,k),' ',
c    .        -ntm0*nlmesp,4,is,0)
            jjj = partok(mfile,'@','@',wk,' ',-ntm0*nlmesp,4,is,0)
            if (jjj /= 0) then
              do k1 = 1, ntm0*nlmesp
                j1 = (k1-1)/ntm0 + 1
                i1 = k1 - (j1-1)*ntm0
                tabov(i1,j1,k) = wk(k1)
              enddo
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',decov(1,k),' ',
     .          -nlmesp,4,is1,0)
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'CUT=',4,'=',is0))
     .          jjj = partok(mfile,'CUT=','=',cutov(1,1,k),' ',
     .          -2*nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZ=',4,'=',is0)) then
c    .            jjj = partok(mfile,'FRZ=','=',itbov(1,k),' ',
c    .            -ntm0*nlmesp,2,is1,0)
                  jjj = partok(mfile,'FRZ=','=',iwk,' ',
     .            -ntm0*nlmesp,2,is1,0)
                  if (jjj /= 0) then
                    do k1 = 1, ntm0*nlmesp
                      j1 = (k1-1)/ntm0 + 1
                      i1 = k1 - (j1-1)*ntm0
                      itbov(i1,j1,k) = iwk(k1)
                    enddo
                  endif
                endif
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idov(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif

C --- Read overlap crystal field ME and decay parameters ---
          if (ocryf) then
c           jjj = partok(mfile,'%','%',tbocf(1,k),' ',
c           -ntm0*nlmesp,4,is,0)
            jjj = partok(mfile,'%','%',wk,' ',-ntm0*nlmesp,4,is,0)
            if (jjj /= 0) then
              do k1 = 1, ntm0*nlmesp
                j1 = (k1-1)/ntm0 + 1
                i1 = k1 - (j1-1)*ntm0
                tbocf(i1,j1,k) = wk(k1)
              enddo
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'DECAY=',6,'=',is0))
     .          jjj = partok(mfile,'DECAY=','=',dcocf(1,k),' ',
     .          -nlmesp,4,is1,0)
              if (fitpar) then
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZ=',4,'=',is0)) then
c    .            jjj = partok(mfile,'FRZ=','=',itocf(1,k),' ',
c    .            -ntm0*nlmesp,2,is1,0)
                  jjj = partok(mfile,'FRZ=','=',iwk,' ',
     .            -ntm0*nlmesp,2,is1,0)
                  if (jjj /= 0) then
                    do k1 = 1, ntm0*nlmesp
                      j1 = (k1-1)/ntm0 + 1
                      i1 = k1 - (j1-1)*ntm0
                      itocf(i1,j1,k) = iwk(k1)
                    enddo
                  endif
                endif
                nt = recln0*memx
                is1 = iend
                call skipbl(mfile,nt,is1)
                is0 = 1
                if (lsequ(mfile(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .            jjj = partok(mfile,'FRZDEC=','=',idocf(1,k),' ',
     .            -nlmesp,2,is1,0)
              endif
            endif
          endif
        endif

C --- Read tight-binding ME and decay parameters ---
        j = 0
        if (ipass == 1 .or. k <= nset) then ! to avoid going outside nset
          nt = 1
          ntm0 = ntm(memode(k))
          if (ipass == 2) then
            nt = ntm0*nlmesp
            if (memode(k) == 6) call rx('check branch')
C           if (memod0 == 6) j = 8
          endif
          j = partok(mfile,'|','|',wk,' ',-nt,4,is,0)
          if (j /= 0) then
            do k1 = 1, nt
              j1 = (k1-1)/ntm0 + 1
              i1 = k1 - (j1-1)*ntm0
c             i1 = mod(k1-1,ntm0) + 1
              tabme(i1,j1,k) = wk(k1)
            enddo
          endif
          if (j /= 0) then
            nt = recln0*memx
            is1 = iend
            call skipbl(mfile,nt,is1)
            is0 = 1
            if (lsequ(mfile(is1+1:is1+1),'DECAY=',6,'=',is0))
     .        jjj = partok(mfile,'DECAY=','=',decay(1,k),' ',
     .        -nlmesp,4,is1,0)
            if (ipass == 2) then
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'CUT=',4,'=',is0))
     .          jjj = partok(mfile,'CUT=','=',cut(1,1,k),' ',
     .          -2*nlmesp,4,is1,0)
            endif
            if (ipass == 2 .and. fitpar) then
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'FRZ=',4,'=',is0)) then
c    .          jjj = partok(mfile,'FRZ=','=',itab(1,k),' ',
c    .          -ntm0*nlmesp,2,is1,0)
                jjj = partok(mfile,'FRZ=','=',iwk,' ',
     .            -ntm0*nlmesp,2,is1,0)
                if (jjj /= 0) then
                  do k1 = 1, ntm0*nlmesp
                    j1 = (k1-1)/ntm0 + 1
                    i1 = k1 - (j1-1)*ntm0
                    itab(i1,j1,k) = iwk(k1)
                  enddo
                endif
              endif
              nt = recln0*memx
              is1 = iend
              call skipbl(mfile,nt,is1)
              is0 = 1
              if (lsequ(mfile(is1+1:is1+1),'FRZDEC=',7,'=',is0))
     .          jjj = partok(mfile,'FRZDEC=','=',idec(1,k),' ',
     .          -nlmesp,2,is1,0)
            endif
          endif
        end if

        if (j /= 0) inext = max(inext,iend)
        if (j == 0 .and. jj == 0) goto 70

        if (ipass == 2) then
C --- Assemble {ij} ---
          do ii = 1, 2
            call skipbl(mfile,memx*recln0,is)
            js = is
            call skp2bl(mfile,memx*recln0,js)
!               call mkilst(mfile(is+1:js+1),nij(ii),ij(1,ii))
                blabl(ii) = ' '
                blabl(ii)(1:js-is+1) = mfile(is+1:js)
                ij(ii) = chaidx(blabl(ii),clabls,nclass)
            is = js
          end do
C --- Generate iam ---
          if (      0 < ij(1) .and. ij(1) <= nclass
     &        .and. 0 < ij(2) .and. ij(2) <= nclass) then
              kk = kk+1
              iam(1,kk) = ij(1)
              iam(2,kk) = ij(2)
              iam(3,kk) = k
          else if (isnum(trim(blabl(1))).and.isnum(trim(blabl(2))))then
          ! let this be .and. instead of .or. just to be sure it is an old format file.

              read(blabl(1),*) ij(1)
              read(blabl(2),*) ij(2)

              if (any(ij > nclass) .or. any(ij < 1)) cycle

              print '(a)',
     &         " RDTBH: The pair specification format you are"//
     &         " trying to use is no longer supported!"
              print '(8x,a)', "Please replace the specie referring "
     &                        //"indices, by the respective labels."
              print'(8x,9a)','For example, in the current case, pair "',
     &            trim(blabl(1)),' ',trim(blabl(2)),'" should become "',
     &            trim(clabls(ij(1))),' ',trim(clabls(ij(2))),'" etc...'
              call rx("RDTBH: obsolete ME format detected")
          else
              ii = 1
              if (0 < ij(1) .and. ij(1) <= nclass) ii = 2
              if (procid == master) print '("No atoms of type ",a,
     &   " were found in the structure. Discarding pair ",a," ",a,".")',
     &          trim(blabl(ii)), trim(blabl(1)), trim(blabl(2))
          endif
        endif

   60 continue
      call rx('RDTBH:  increase memx')

   70 continue
      k = k-1

C --- Sort iam, generate npm ---
      if (ipass == 2) then
        k = kk
        allocate(wka(k))
        call ivshel(3,k,iam,wka,.false.)
        deallocate(wka)

        if (iprint() > 40) then
          print '(/'' RDTBH: '',i0, '' matrix element pairs found'')',k
          if (iprint() >= 50) print 333, ((iam(i,j), i=1,3), j=1,k)
  333     format(2i5,3x,i5)
        endif

        call iinit(npm,2*nclass)
        do  80  j = 1, k
          if (iam(1,j) > nclass .or. iam(2,j) > nclass .or.
     .        iam(1,j) < 1    .or. iam(2,j) < 1) then
            print *, ' j, iam(1,j), iam(2,j)=', j, iam(1,j), iam(2,j)
            call rx('RDTBH:  bad basis index in category ME')
          endif
          npm(0,iam(1,j)) = npm(0,iam(1,j))+1
   80   continue

        npm(1,1) = 0
        if (nclass >= 2) then
          do  i = 2, nclass
            npm(1,i) = npm(1,i-1) + npm(0,i-1)
          enddo
        endif

      endif

      end


      function c_isnum(c) result (r)
C- Figures whether 'c' is a digit.
C ----------------------------------------------------------------------
          implicit none
          character, intent(in) :: c
          logical :: r
          integer :: i
          i = iachar(c)
          r = (iachar('0') <= i) .and. (i <= iachar('9'))
      end function c_isnum


      function isnum(a)
C- Tells whether the string 'a' is composed of digits only.
C  A shoddy implementation
C ----------------------------------------------------------------------
          implicit none
          character(len=*), intent(in) :: a
          logical :: isnum
          procedure(logical) :: c_isnum

          integer :: i, n
          n = len(a)

          isnum = .true.

          do i = 1, n
              isnum = isnum .and. c_isnum(a(i:i))
          end do

      end function isnum

      function chaidx(e,a,n)
C- Find the position of the element "e" in array "a" searching linearly
C                                from left. (for character(len=*) type)
C ----------------------------------------------------------------------
Ci Inputs
Ci   e  : a string
Ci   a  : array of strings
Ci   n  : size of "a"
Co Outputs
Co   chaidx : the index of the element in the array which matches 'e'.
Co            Counting from 1. If not found return -1.
!DMT
          implicit none
          integer, intent(in) :: n
          character(len=*), intent(in) :: e, a(n)

          integer :: chaidx

          integer :: ne

          ne = len(e)

          chaidx = 1
          do while (chaidx <= n)
            if (a(chaidx)(1:ne) == e) exit ! this could have been put in the while but some compiler and valgrind may complain for the case chaidx==n
            chaidx = chaidx + 1
          end do
          if (chaidx == n + 1) chaidx = -1

      end function chaidx


