      subroutine iomagf(mode,nbas,nl,bfield,bxc,nbf,ifi)
C- I/O of external magnetic field, or exchange-correlation field
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  :0 I/O external magnetic field; may be orbital-dependent
Ci         :1 I/O exchange-correlation field: only site-dependent
Ci         :2 I/O scaling factors
Ci   nbas  :number of sites or classes
Ci   nl    :(global maximum l) + 1
Ci   ifi   :file logical unit, <0 for write, >0 for read
Co Inputs/Outputs
Cio  bfield:(mode 0) bfield vectors are read in or written out
Cio  bxc   :(mode 1,2) bxc vectors are read in or written out
Cio        :bxc is dimensioned (3,ib), where ib is site or class index.
Cio        :    (nbas can correspond either to sites or classes)
Cio        :The meaning of bxc depends on whether mode is 1 or 2.
Cio        :mode=1 bxc = vector field in the noncollinear magnetic code
Cio        :mode=2 bxc = scaling factors affecting the LDA hamiltonian
Cio        :       It is called shfac elsewhere; see structures.h
Cio        :       bxc(1,:) = scaling factor for the XC field
Cio        :       Its precise meaning depends on s_ctrl%lbxc.
Cio        :       s_ctrl%lbxc = 1 => scale field as 1+bxc(1,ib)
Cio        :       s_ctrl%lbxc = 2 => scale field as 1+bxc^2(1,ib)
Cio        :       bxc(2,:) = scale SO parameter lambda by 1+bxc(1,ib)
Cio        :       bxc(3,:) = scale channel bandwidth by 1+bxc(1,ib)
Cio  nbf  is read in or written out: one of 1, nl or nl**2 (mode=0 only)
Cr Remarks
Cr   Aborts on read when nbas does not match file
Cr   19 Nov 97 use rdfiln for reading first line
Cu Updates
Cu   07 Jun 12 New mode 2
Cu   12 Apr 11 Extended to read bxc
Cu   27 Mar 04 Extended to read bfield
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,nbf,nl,ifi
      double precision bfield(nbas,nbf,3),bxc(3,nbas)
C ... Local parameters
      integer i,j,k,ipr,nc,rdm
      logical parstr,a2bin
      character*80 ss2*4
      double precision tmp(10000)
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),ss*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
      equivalence (a,ss)
      data rdarg /'#{}% c'/

      if (mode == 1 .and. nbf /= 1) call rx('iomagf: illegal nbf')

C --- Write ---
      if (ifi < 0) then
C       rewind (-ifi)
        if (nbf /= 1 .and. nbf /= nl .and. nbf /= nl*nl)
     .    call rxi('iomagf: bad dim for bfield: nbfield=',nbf)
        nc = 3
        call awrit3('%% rows %i cols %i nbf %i',ss,80,0,nbas*nbf,nc,nbf)
        call awrit0('%a',ss,80,ifi)
        do  i = 1, nbas
          if (nbf > 1) call awrit1('# ib %,4i',' ',80,-ifi,i)
          if (mode == 0) then
            write(-ifi,'(3f16.12)') ((bfield(i,j,k),k=1,3),j=1,nbf)
          else
            write(-ifi,'(3f16.12)') (bxc(k,i),k=1,3)
          endif
        enddo
C --- Read ---
      else
        call getpr(ipr)
        rewind ifi
        nr = 0
        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,ss,recrd,recl,nr)
C   ... No error if file is empty, but don't bother reading
        if (nr == 0) then
          call info0(30,0,0,' IOMAGF:  empty file')
          return
        endif
        rewind ifi
C   ... Read nbf if it is there
        if (mode == 0) then
          if (ss(1:1) == '%') then
            nbf = 1
            i = 0
            if (parstr(ss,'nbf ',len(ss)-4,4,' ',i,j)) then
              if (.not. a2bin(ss,nbf,2,0,' ',j,len(ss)))
     .          call rx('iomagf: failed to parse for nbf')
            endif
          else
            nbf = 1
          endif
        endif

C   ... Expect nbas*nbf*3 rows in file ... skip reading if not
        if (mode == 0 .or. mode == 1) then
          if (nbf /= 1 .and. nbf /= nl .and. nbf /= nl*nl)
     .      call rxi('iomagf: illegal dimension: nbf =',nbf)
          ss2 = 'in'
          if (nbf == nl) ss2 = ' '
          if (nbf == nl*nl) ss2 = '%bm-'
          call info2(31,0,0,' IOMAGF: reading l-'//ss2//
     .      '%adependent magnetic field (file %i)',ifi,0)
          nc = 3
        else
          call info2(31,0,0,' IOMAGF: reading '//
     .      'constraining fields for CPA (file %i)',ifi,0)
          nc = 3
        endif

C   ... Read the data, 1 atom at a time ...
        if (mode == 0) call dpzero(bfield,nbas*nbf*nc)
        if (mode == 1 .or. mode == 2) call dpzero(bxc,nbas*nc)
        if (mode == 0) then
        do  i = 1, nbas
          call rxx(rdm(ifi,0,nbf*nc,' ',tmp,nbf,nc) /= 1,
     .      'IOMAGF:  bad file')
C         if (mode == 0)
          call dmcpy(tmp,nbf,1,bfield(i,1,1),nbas*nbf,nbas,nbf,nc)
C         if (mode == 1)
C    .      call dcopy(3,tmp,1,bxc(1,i),1)
          if (ipr >= 50) call pvioeu(11,0d0,bfield,nbas,nbf)
        enddo
      elseif (mode == 1 .or. mode == 2) then
          if (nbas*nbf*nc > 10000) call rx('IOMAGF: Increase tmp')
          k = rdm(ifi,100000,nbas*nbf*nc,' ',tmp,nbas*nbf,nc)
          call rxx(k /= 1,'IOMAGF:  bad file')
          do  i = 1, nbas
            call dcopy(nc,tmp(i),nbas,bxc(1,i),1)
          enddo
        endif
      endif
      end
