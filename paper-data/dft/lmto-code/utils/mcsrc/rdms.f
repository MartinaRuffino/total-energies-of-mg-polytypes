      integer function rdms(ifi,sw,mxelt,filel,isprs,s,nr,nc)
C- Matrix input, sparse format parsing file with rdfiln.
C ----------------------------------------------------------------
Ci  nr,nc: number of rows and columns (also output; see Remarks)
Ci  sw:    one's digit: 1 'quick read'--- read array with fortran read
Ci                      2 'binary read'--- read a binary file
Ci                      3 'binary read'--- nr,nc,ic not read
Ci                        (nr and nc must be passed; array assumed real)
Ci         ten's digit: 0, accept real matrix s(nr,nc) only
Ci                      1, accept real s(nr,nc) or cmplx s(nr,nc,2)
Ci                      2, accept real s(nr,nc) or cmplx s(2,nr,nc)
Ci                         NB: implemented only for binary read
Ci                      3, accept cmplx s(nr,nc,2) only (See remarks)
Ci                      4, accept cmplx s(2,nr,nc) only
Ci        1000's digit: 0, ignore label, if any exists
Ci                      1, return label in string filel
Ci       10000's digit: 0, matrix elements in normal order, ie
Ci                         s(1,1), s(1,2), ..., s(1,nc),
Ci                         s(2,1), s(2,2), ..., s(2,nc)...
Ci                      1, matrix elements in transpose order ie
Ci                         s(1,1), s(2,1), ..., s(nr,1),
Ci                         s(1,2), s(2,2), ..., s(nr,2)...
Ci                       NB: for this case, nr and nc must be spec'd
Ci                       explicitly, ie by priority 1 or 2 in Remarks
Ci                       NB: can also be set in the first line of the
Ci                       input; see Remarks
Ci  mxelt:number of nonzero elements
Ci        mxelt=0 => determine nr,nc,mxelt but make no attempt to read s
Co  s:    array for which
Co  rdm:  1 for real, 2 for complex, 3 for complex in complex*16 format,
Co        provided it evaluates all expressions with no errors and reads
Co        enough expressions to be consistent with nr and nc specified.
Co        If not, returns rdm = -1.
Cr Remarks
Cr  rdm attempts to read in matrix s by parsing file ifi, using routine
Cr  rdfiln for conditional reading of lines, looping syntax and
Cr  expression substitutions.  The number of rows nr and columns nc of s
Cr  need not be prescribed in advance.  rdm determines nr and nc
Cr  according to the following priority.
Cr  1. If on entry nc>0, nc always takes that value.  rdm will not set
Cr     nc.  If file also specifies nc (see below) not compatible
Cr     with passed nc, rdm fails and returns -1.  Similarly for nr.
Cr  ASCII read:
Cr  2. rdm reads the first line of the file.  If it has the syntax
Cr     % rows expr1 cols expr2
Cr     nr is assigned value expr1 and nc value expr2,
Cr     provided neither has been set in step 1
Cr     (NB: either rows declaration or cols declaration may be absent).
Cr  3. rdm parses the first line, and sets nc to the number of
Cr     expressions it finds on that line.  At this point, rdm has an
Cr     unambiguous value for nc (or else it has aborted)
Cr  4. If nr has not yet been assigned, the file is parsed until EOF is
Cr     reached or it encounters a line beginning with '%'.
Cr     nr is assigned the (total no expressions evaluated) / nc.
Cr     When nr is determined in this way, mxelt must be double the
Cr     number of elements in s, since rdm must read in s-transpose,
Cr     and copy s-transpose back to s.
Cr     NB: nr may not be determined this way using 'quick-read'
Cr  BINARY read:
Cr  2. rdm reads the first record in the file.  It must contain
Cr     nr nc ic
Cr     It MAY contain additionally
Cr     nr nc ic iswitch   <-- iswitch indicates next recrd is label
Cr     If nr is >0 on input, file nr must match passed nr; ditto for nc.
Cr     NB: this step is missing if 1's digit of sw is 3
Cr  Complex matrices: A file is unambiguously
Cr     specified to be real or complex if the first line begins
Cr     with % and contain 'real' or 'complex'  (eg, '% rows 7 complex').
Cr     Otherwise, it is assumed to real unless a complex matrix was
Cr     requested (10's digit of sw=3,4).  rdm will convert if the cast
Cr     sought conflicts with the cast specified.  The resultant cast
Cr     depends on the ten's digit of sw, as follows:
Cr             |   sw=0         sw=3,4       sw=1,2
Cr             |  (want R)     (want C)     (unspec)
Cr      -------|---------------------------------------
Cr     Have R  |      R          R->C          R
Cr     Have C  |   C->R             C          C
Cr     unspec  |      R             C          R
Cr  Tranpose of matrices: rdm will read the matrix elements in transpose
Cr    order if the appropriate switch is set (see sw, above) or if the
Cr    first line contains
Cr    % .... trans   ...
Cr
Cb  Bugs: no check to see whether matrix asserted to be symmetric
Cb        (hermitian) actually is.
Cu  Updates
Cu  17 Nov 98  rdm can read real sparse matrices
Cu  17 Mar 99  rdm can read integer array (converting it to real)
C ----------------------------------------------------------------
      implicit none
      integer ifi,mxelt,nr,nc,sw,isprs(mxelt,2)
      character*(*) filel
      double precision s(1)
C Local variables
      integer recl,nl,wantC,haveC,sw10,ncx,ioff
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=1500,ctlen=120)
      character recrd*(recl), a*(recl)
      character ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev),a2bin,qread,ltrnsi,ltrns,pvrdm1,cmdopt,
     .  lddump,lidump,lpr,lsprse
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev),i1mach,i,ip,k,nelt,ic,ir,retval,ii,ilbl,nrx,ipr,
     .  ndir,a2vec,ix(3),ilast
      double precision xwk(3)
      parameter (ndir=8)
      character vnam(mxlev)*16,dir(ndir)*7,llbl*1000,rdarg*6
      data rdarg /'#{}% c'/
      data dir /'rows','cols','real','symm','complex','herm',
     .  'trans','sparse'/

C --- Initialization ---
      call getpr(ipr)
      lpr = ipr > 0
      ltrnsi = .false.
      ltrns = mod(sw/10000,10) == 1
      lsprse = .false.
      qread = mod(sw,10) == 1
      sw10 = mod(sw/10,10)
C     wantC and haveC:  0, unspecified, 1 real, 2 complex
      wantC = 1
      if (sw10 == 1 .or. sw10 == 2) wantC = 0
      if (sw10 == 3 .or. sw10 == 4) wantC = 2
      haveC = 0
      nelt = 0
      ir = 0
      ic = 0
      ioff = 0
      rdms = -max(wantC,1)
      nl = 0
      if (mxelt > 0) then
        call dpzero(s,mxelt)
      endif
      retval = 0

C --- Determine nr,nc from % rows and % cols ---
      recrd = ' '
   21 call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
C ... Empty file kosher unless nr or nc nonzero
      if (nl <= 0 .and. nr == 0 .and. nc == 0) rdms = max(wantC,1)
      if (nl <= 0) return
C ... Skip blank record
      if (recrd == ' ') goto 21
C ... First record: check for % rows and % cols
      if (recrd(1:1) == '%') then
        ip = 1
   20   call skipbl(recrd,recl,ip)
        if (ip >= recl) goto 30
        k = ip-1
        call tokmat(recrd(ip+1:recl),dir,ndir,7,' ',i,ip,.false.)
        ip = ip+k
C       print *, 'after tokmat', i,ip
        if (i < 0) then
          call skp2bl(recrd,recl,ip)
          goto 20
        endif
C   ... Read matrix elements in transpose order?
        if (i == 6) then
          ltrns = .true.
          goto 20
C   ... Read sparse matrix?
        elseif (i == 7) then
          lsprse = .true.
          s(1) = 0
          goto 20
        endif
C   ... Matrix is specified real or complex
        if (i >= 2 .and. i <= 5) then
          haveC = 1
          if (i > 3) haveC = 2
          if (i == 3 .or. i == 5) retval = 10
          goto 20
        endif
        if (ip >= recl) goto 30
        if (.not. (a2bin(recrd,s,4,0,' ',ip,recl))) return
        if (i == 0 .and. nr /= 0 .and. nr /= nint(s(1))) goto 196
        if (i == 1 .and. nc /= 0 .and. nc /= nint(s(1))) goto 196
        if (i == 0 .and. nr == 0) nr = nint(s(1))
        if (i == 1 .and. nc == 0) nc = nint(s(1))
        goto 20
      else
C   ... First record not '%': skip reading this first record
        if (nr /= 0 .and. nc /= 0 .and. mxelt == 0) goto 98
        goto  31
      endif

C --- Read next record from file ---
   30 continue
      recrd = ' '
      if (.not. (qread .and. lsprse))
     .  call rdfiln(ifi,'#{}%',mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
   31 continue

      if (.not. lsprse)
     .  call rx('rdms only set up to read sparse matrix input')

      if (ltrns)
     .  call rx('rdms not set up to transpose sparse matrix')

      if ((nr == 0 .or. nc == 0) .and. lsprse)
     .  call rx('rdms not set up to read sparse matrix'//
     .  ' without explicit nr,nc')

C --- Before reading, decide whether file is real or complex ---
      if (nelt == 0) then
        ilast = 0
        if (haveC /= 0) then
          if (wantC == 0) wantC = haveC
        elseif (wantC == 1 .or. wantC == 2) then
          haveC = wantC
        else
          wantC = 1
          haveC = 1
        endif
        if (haveC == 1 .and. wantC == 2 .and. lpr) print '(a)',
     .    '#rdms (warning): sought complex matrix but file is real'
        if (haveC == 2 .and. wantC == 1 .and. lpr) print '(a)',
     .    '#rdms (warning): sought real matrix but file is complex'
        if (wantC == 0 .or. haveC == 0) stop 'bug in rdms'
      endif

C --- Parse this record for ir,ic,s(nelt) ---
C ... Exit if we are done reading this file
      if (recrd(1:1) == '%') goto 70
      if (nl <= 0) goto 70
      ip = 0
   61 call skipbl(recrd,recl,ip)
C   ... Pick up next record when this one parsed to the end
        if (ip >= recl) goto 30
C   ... Read into s(nelt) (or scratch if mxelt is zero)
        if (lsprse) then
          if (a2vec(recrd,recl,ip,4,' ',1,2,3,ix,xwk) /= 3) return
          ir = xwk(1)
          ic = xwk(2)
C         first two elements must be integers < nr,nc
          if (abs(ir-xwk(1))+abs(ic-xwk(2)) /= 0 .or.
     .        ir > nr .or. ic > nc .or. ic < ilast) goto 96
          nelt = nelt+1
          if (mxelt == 0) goto 30
C         New column
          if (ic > ilast) then
            do  62  i = ilast+1, ic
              isprs(i,1) = nelt
   62       continue
            ilast = ic
          endif
C         Poke row index
          isprs(nelt,2) = ir
          s(nelt) = xwk(3)
          goto 30
        endif
      goto 61

C --- Cleanup.  If nr is zero, set nr and copy s-transpose to s ---
   70 continue
      if (mxelt == 0) then
        mxelt = nelt
      else
        isprs(nc+1,1) = mxelt+1
      endif

C --- Normal exit ---
   98 rdms = 100 + max(wantC,1) + retval
      return

C --- Exit when inappropriate format ---
   96 continue
      print '(a)', '#rdms: matrix in inappropriate format ...'
      return

C --- Exit when not enough space to read s ---
   99 continue
      print '(a)', '#rdms: insufficent space ...'
      return

C --- Exit for incompatible input ---
  196 continue
      print '(a)', '#rdms: input dimensions incompatible with file'
      rdms = -1
      return

C --- Exit for binary read ---
  198 continue
      if (lpr) print '(a)', '#rdms: binary read failed...'
      rdms = -1
  199 continue
      end
      subroutine cprms(opt,filel,icast,ifi,fmt,nnz,isprs,s,nr,nc)
C- Writes sparse matrix to file ifi, in dense format
C ----------------------------------------------------------------
Ci Inputs:
Ci   opt:   1s digit:
Ci           0 writes in ascii mode
Ci           1 writes in binary mode
Ci  filel:   reference string
Ci  icast:   0 integer
Ci           1 real
Ci           2 complex with imaginary following real
Ci           3 double complex
Ci      Add 10 to indicate symmetric or hermitian
Ci    ifi:   file logical unit
Ci    fmt:   fortran format to write ascii string, e.g. '(5f12.6)'
Ci           You can use awrite format, e.g. '(%5,6g)'
Ci      s:   matrix to be printed out
Ci  ns,nr:   leading dimension of s and number of rows to write
Ci nsc,nc:   column dimension of s and number of columns to write
Cu  Updates
Cu  17 Mar 99  cprm can write integer array
C ----------------------------------------------------------------
      implicit none
      integer opt,icast,nr,nc,ifi,ns,nsc,nnz,isprs(nnz,2)
      character*(*) filel,fmt
      double precision s(nnz)
      character outs*80,lfmt*40,sout*500
      integer i,j,nw,iw(100),lf,ic,jc,k,ip,ir,imm,lbin
      logical a2bin,ltmp,gettok
      double precision sijc(2),sijr,siji
      lbin = mod(opt,10)

C --- Binary write ---
      if (lbin > 0) then
        call rxi('cprms: bad option',opt)
      endif

      if (mod(icast,10) == 0)
     .  call rx('cprms not ready for ascii write of integer')
      outs = ' '
      ic = mod(icast,100)
      if (ic == 1)  outs = ' real l="' //filel
      if (ic == 11) outs = ' symm l="' //filel
      if (ic == 2)  outs = ' complex l="' //filel
      if (ic == 12) outs = ' herm l="' //filel
      if (ic == 3)  outs = ' complex l="' //filel
      if (ic == 13) outs = ' herm l="' //filel
      call word(outs,2,i,j)
      if (outs(i:) == 'l="') then
        outs(i:) = ' '
      else
        call skpblb(outs,len(outs),i)
        outs(i+2:i+2) = '"'
      endif
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)

C --- Case awrite format ---
      if (fmt(2:2) == '%') then
        if (ic == 3 .or. ic == 13)
     .   call rx('cprm cannot handle %% fmt for complex*16')

C   ... copy to lfmt to get around SGI compiler bug
        lfmt = fmt
        lf = len(lfmt)
C   ... Find end of first expression
        j = 2
   30   ltmp = gettok(lfmt,j,' ',lf)
        if (ltmp) goto 30
        nw = 1
        if (j > 2) then
          i = 2
          call rxx(.not. a2bin(lfmt,nw,2,0,' ',i,j),
     .      'cprm failed to parse format %'//lfmt(2:j+2))
        endif
C   ... Loop through all rows to find col-dependent size
        nw = min(nw,100)
        call iinit(iw,nw)
        call skpblb(lfmt,lf,i)
        outs = '%'//lfmt(j+1:i)
        call rx('cprms not implemented this mode')
        do  32  i = 1, nr
C          call prmwid(s(i,1,1),nw,outs,iw,ns,nc)
C          if (mod(ic,10) > 1)
C     .    call prmwid(s(i,1,2),nw,outs,iw,ns,nc)
   32   continue
*       call awrit2('%n:1i',' ',80,6,nw,iw)

C   ... Write out matrix
        imm = 1
   44   continue
        do  40  ir = 1, nr
        do  40  ic = 1, nc, nw
          jc = min(ic+nw-1,nc)
          k = 0
          sout = ' '
          ip = 0
          do  42  i = ic, jc
            k = k+1
            call sij(ir,i,icast,isprs,s,sijc(1),sijc(2))
            call awrit1(outs,lfmt,len(lfmt),0,sijc(imm))
            call skpblb(lfmt,iw(k),j)
            call sij(ir,i,icast,isprs,s,sijc(1),sijc(2))
            call awrit2('%np'//outs,sout,len(sout),0,
     .        ip+iw(k)-j,sijc(imm))
            ip = ip + iw(k)
   42     continue
          call awrit0('%a',sout,-len(sout),-ifi)
   40   continue
        if (imm == 1 .and. mod(ic,10) > 1) then
          imm = 2
          write(ifi,'(1x)')
          goto 44
        endif
        return
      endif

      do  10  i = 1, nr
   10 write(ifi,fmt) (sijr(i,j,icast,nnz,isprs,s), j=1,nc)
      if (mod(ic,10) > 1) then
        write(ifi,'(1x)')
        do  20  i = 1, nr
   20   write(ifi,fmt) (siji(i,j,icast,nnz,isprs,s), j=1,nc)
      endif
      end
      double precision function sijr(ir,ic,icast,nnz,isprs,s)
C return real part of element s(ir,ic) from sparse matrix
      implicit none
      integer ir,ic,icast,nnz,isprs(nnz,2)
      double precision s(nnz)
      double precision sijc(2)
      call sij(ir,ic,icast,nnz,isprs,s,sijc)
      sijr = sijc(1)
      end
      double precision function siji(ir,ic,icast,nnz,isprs,s)
C return imaginary part of element s(ir,ic) from sparse matrix
      implicit none
      integer ir,ic,icast,nnz,isprs(nnz,2)
      double precision s(nnz)
      double precision sijc(2)
      call sij(ir,ic,icast,nnz,isprs,s,sijc)
      siji = sijc(2)
      end
      subroutine sij(ir,ic,icast,nnz,isprs,s,sijc)
C extract element s(ir,ic) from sparse matrix
      implicit none
      integer ir,ic,icast,nnz,isprs(nnz,2)
      double precision s(nnz),sijc(2)
      integer i

      if (mod(icast,10) == 2) call rx('sij not impl. for cmplx')
      sijc(1) = 0
      sijc(2) = 0
C     Try and find a row matching ir
      do  10  i = isprs(ic,1), isprs(ic+1,1)-1
        if (ir == isprs(i,2)) then
          sijc(1) = s(i)
          return
        endif
   10 continue

      end
