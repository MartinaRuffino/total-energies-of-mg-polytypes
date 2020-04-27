      integer function rdm(ifi,sw,mxelt,filel,s,nr,nc)
C- Matrix input, parsing file with rdfiln.
C ----------------------------------------------------------------
Ci  nr,nc: number of rows and columns (also output; see Remarks)
Ci  sw:    one's digit: 1 'quick read'--- read array with fortran read
Ci                      2 'binary read'--- read a binary file
Ci                      3 'binary read'--- nr,nc,ic not read
Ci                        (nr and nc must be passed; array assumed real)
Ci                      4 'binary read'--- same as 3 except that
Ci                        array is assumed to be c16 format
Ci                        (NOT IMPLEMENTED)
Ci         ten's digit: 0, accept real matrix s(nr,nc) only
Ci                      1, accept real s(nr,nc) or cmplx s(nr,nc,2)
Ci                      2, accept real s(nr,nc) or cmplx s(2,nr,nc)
Ci                         NB: implemented only for binary read
Ci                      3, accept cmplx s(nr,nc,2) only (See remarks)
Ci                      4, accept cmplx s(2,nr,nc) only
Ci        1000's digit: 0, ignore label if any exists
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
Ci       100000's digit: Governs what rdm does if number of elements
Ci                       not sufficient to fill array
Ci                      0, rdm returns > 0 accompanied by a warning
Ci                      1, rdm returns < 0
Ci  mxelt: maximum dimension of s (considered as a 1-dimensional array)
Ci         mxelt=0 => determine nr,nc but make no attempt to read s
Co  filel :text label describing file contents (see 1000s digit sw)
Co  s:    Data is read into this array
Co  rdm:  1 for real, 2 for complex, 3 for complex in complex*16 format,
Co        provided it evaluates all expressions with no errors and reads
Co        enough expressions to be consistent with nr and nc specified.
Co        If not, returns rdm = <0.
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
Cr     nr nc ic iswitch
Cr       where:
Cr     If nr is >0 on input, file nr must match passed nr; ditto for nc.
Cr     ic is the cast (0=int, 1=dble, 2=cmplx, split real, imag, 3=c*16)
Cr     iswitch indicates next recrd is label
Cr   NB: this step is missing if sw = 1's digit of ic is 3
Cr  Complex matrices: a file is unambiguously
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
Cu  07 Feb 13 When transposing array, allocates space internally
Cu  25 Sep 09 bug fix: rdm returns after reading nr,nc in binary mode
Cu            when mxelt=0
Cu  09 Nov 05 rdm can read complex arrays into c*16 format
Cu  16 jan 01 rdm can read complex arrays, sparse format
Cu  17 Mar 99 rdm can read integer array (converting it to real)
Cu  17 Nov 98 rdm can read real sparse matrices
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,mxelt,nr,nc,sw
      character*(*) filel
      double precision s(*)
C ... Local parameters
      integer recl,nl,wantC,haveC,sw10,ncx,ioff,nrl,ncl
      integer mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=2000,recl=1024,ctlen=120)
      character recrd*(recl), a*(recl)
      character ctbl(mxchr,2)*(ctlen)
      logical loop0(0:mxlev),a2bin,qread,ltrnsi,ltrns,pvrdm1,
     .  lddump,lidump,lpr,lsprse
C#ifdefC CCAR
C      logical cmdopt
C#endif
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev),i1mach,i,ip,k,nelt,ic,ir,retval,ii,ilbl,nrx,ipr,
     .  ndir,a2vec,ix(4),sw1
      double precision xwk(4),xx
      real(8), pointer:: swk(:)
      parameter (ndir=8)
      character vnam(mxlev)*16,dir(ndir)*7,llbl*1000,rdarg*7
      data rdarg /'#{}% ct'/
      data dir /'rows','cols','real','symm','complex','herm',
     .  'trans','sparse'/

C --- Initialization ---
      call getpr(ipr)
      lpr = ipr > 1
      ltrnsi = .false.
      ltrns = mod(sw/10000,10) == 1
      lsprse = .false.
      sw1 = mod(sw,10)
      qread = sw1 == 1
      sw10 = mod(sw/10,10)
      call dpzero(xwk,4)
C     wantC and haveC:  0, unspecified, 1 real, 2 complex
      wantC = 1
      if (sw10 == 1 .or. sw10 == 2) wantC = 0
      if (sw10 == 3 .or. sw10 == 4) wantC = 2
      haveC = 0
      nelt = 0
      ir = 0
      ic = 0
      ioff = 0
      rdm = -max(wantC,1) ! Until rdm succeeds, assume failure
      nl = 0
      if (mxelt > 0) call dpzero(s,mxelt)
      retval = 0
      nrl = nr ! Work with local copy so as not to overwrite
      ncl = nc ! nr,nc until program exits

C --- Binary read ---
      if (sw1 == 2 .or. sw1 == 3) then
        if (sw1 == 2) then
C         Read number of rows, columns, cast
          read(ifi,err=198,end=198) nrl,ncl,ic
          if (nr == 0) nr = nrl
          if (nc == 0) nc = ncl
          if (nrl /= nr .or. ncl /= nc) goto 196 ! incompatible dims
C     ... look for a label
          backspace ifi
          read(ifi,err=197,end=197) ii,ii,ii,ilbl
          if (ilbl == 0) goto 195 ! if no label
C     ... there is a label; read into llbl
          llbl = ' '
          read(ifi,err=198,end=198) llbl(1:min(len(llbl),ilbl))
          if (mod(sw/1000,10) /= 0) filel = llbl
          goto 195
C     ... recover from missing ilbl
  197     continue
          backspace ifi
          read(ifi,err=198,end=198) ii
C         read(ifi) s(1)
  195     continue
        else
          ic = 1
          if (nr <= 0 .or. nc <= 0) goto 199
        endif
        haveC = mod(ic,10)
        if (haveC == 3) haveC = 2
C       if (ic == 0 .and. wantC > 0) then
C       If mxelt0, just return nr,nc
        if (mxelt == 0) then
          rdm = ic
          return
        endif
        if (ic == 0) then
          if (mxelt < 2*nr*nc) goto 99
          if (.not. lidump(s(1+nr*nc),nr*nc,ifi)) goto 198
          call idscop(s(1+nr*nc),s,nr*nc,1,1)
          haveC = 1
          ic = 1
C       elseif (ic == 0) then
C          if (mxelt < nr*nc) goto 99
C         if (.not. lidump(s,nr*nc,ifi)) goto 198
C         rdm = ic
C         return
        else
          if (nr*nc*haveC > mxelt) goto 99
          if (.not. lddump(s,nr*nc*haveC,ifi)) goto 198
          nelt = nr*nc*haveC
        endif
        if (haveC == 1 .and. wantC == 2 .and. lpr) print '(a)',
     .    '#rdm (warning): sought complex matrix but file is real'
        if (haveC >= 2 .and. wantC == 1 .and. lpr) print '(a)',
     .    '#rdm (warning): sought real matrix but file is complex'
        rdm = ic
        ltrnsi = ltrns
        goto 199
      endif

      if (sw10 == 2)
     .  call rx('rdm not ready for complex*16 s')

C --- Determine nr,nc from % rows and % cols ---
      recrd = ' '
C#ifdefC CCAR
C      if (cmdopt('-C',2,0,a)) rdarg(1:1) = a(3:3)
C#endif
   21 call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
C ... Empty file kosher unless nr or nc nonzero
      if (nl <= 0 .and. nr == 0 .and. nc == 0) rdm = max(wantC,1)
      if (nl <= 0) return
C ... Skip blank record
      if (recrd == ' ') goto 21
C ... First record: check for % rows and % cols
      if (recrd(1:1) == '%') then
        ip = 1
   20   call skipbl(recrd,recl,ip)
        if (ip >= recl) goto 31
        k = ip-1
        call tokmat(recrd(ip+1:recl),dir,ndir,7,' ',i,ip,.false.)
        ip = ip+k
C       print *, 'after tokmat', i,ip
        if (i < 0) then  ! tag not recognized
          call skp2bl(recrd,recl,ip)
          goto 20
        endif

        if (i == 6) then ! "trans" => Read matrix elements in transpose order
          ltrns = .true.
          if (sw10 == 4) call rx('rdm not ready for complex*16 s')
          goto 20

        elseif (i == 7) then ! "sparse" => Read sparse matrix
          lsprse = .true.
          s(1) = 0
          if (sw10 == 4) call rx('rdm not ready for complex*16 s')
          goto 20
        endif

        if (i >= 2 .and. i <= 5) then
          haveC = 1 ! Default to real
          if (i > 3) haveC = 2 ! "complex" or "herm"
          if (i == 3 .or. i == 5) retval = 10 ! "symm" or "herm"
          goto 20
        endif

        if (ip >= recl) goto 31
        if (.not. (a2bin(recrd,xx,4,0,' ',ip,recl))) return

C       'rows'
        if (i == 0) then
          nrl = nint(xx)
          if (nr /= 0 .and. nr /= nrl) goto 196
          if (nr == 0) nr = nrl
        endif

C       'cols'
        if (i == 1) then
          ncl = nint(xx)
          if (nc /= 0 .and. nc /= ncl) goto 196
          if (nc == 0) nc = ncl
        endif
        goto 20

C ... First record not '%': jump to reading array body
      else
        if (nr /= 0 .and. nc /= 0 .and. mxelt == 0) goto 98
        nrl=nr; ncl=nc ! current line already contains the array body:
        goto  32       ! skip past reading of this record
      endif

C --- Entry point for reading array body ---
   31 continue
      if (nr /= 0 .and. nc /= 0 .and. mxelt == 0) goto 98

C --- Read next record from file ---
   30 continue
      recrd = ' '
      if (.not. (qread .and. lsprse))
     .  call rdfiln(ifi,'#{}%',mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nl)
   32 continue

C --- Before reading, decide whether file is real or complex ---
      if (nelt == 0) then
        if (haveC /= 0) then
          if (wantC == 0) wantC = haveC
        elseif (wantC == 1 .or. wantC == 2) then
          haveC = wantC
        else
          wantC = 1
          haveC = 1
        endif
        if (haveC == 1 .and. wantC == 2 .and. lpr) print '(a)',
     .    '#rdm (warning): sought complex matrix but file is real'
        if (haveC == 2 .and. wantC == 1 .and. lpr) print '(a)',
     .    '#rdm (warning): sought real matrix but file is complex'
        if (wantC == 0 .or. haveC == 0) stop 'bug in rdm'
      endif

C --- nr and nc must be spec'd if lsprse ---
      if ((nrl == 0 .or. ncl == 0) .and. lsprse)
     .  call rx('rdm not set up to read sparse matrix'//
     .  ' without explicit nr,nc')

C --- if nc zero, set nc to number of expressions in 1st recrd ---
      k = 0
      if ((ncl == 0 .and. .not. ltrns .or.
     .     nrl == 0 .and. ltrns .or. qread) .and. .not. lsprse) then
        if (ltrns) then ! Case ltrns: nc and nr exchange roles
          ip = nrl
          nrl = ncl
          ncl = ip
        endif

C       Entry point for next expression in this record
        ip = 0
   41   call skipbl(recrd,recl,ip) ! Next expression in record
        if (ip >= recl) then ! We are at end of record
          if (ncl == 0) then
            ncl = k
            if (ltrns) nr = ncl
            if (.not. ltrns) nc = ncl
          endif
          goto 42
        endif
C   ... Read into scratch unless qread (then load as if s a vector)
        i = 0
        if (qread .and. mxelt > 0) i=k
        if (i > mxelt) goto 99
        if (.not. a2bin(recrd,s,4,i,' ',ip,recl)) return
        k = k+1
        goto 41
   42 continue ! no more elements in this record

C ... Restore nr, nc to proper roles
      if (ltrns) then
        ip = nrl
        nrl = ncl
        ncl = ip
        endif
      endif

C ... nrx,ncx are  nr,nc  or if ltrns,  nc,nr
      nrx = nrl
      ncx = ncl
      if (ltrns) then
        nrx = ncl
        ncx = nrl
      endif

C ... Now ncx's value is known.  ncx=0 => empty file.
      if (ncx == 0 .or. (nrx /= 0 .and. mxelt == 0)) then
        rdm = max(wantC,1) + retval
        return
      endif

C --- Quick-read: read remaining elements with unformatted read ---
   55 continue
      if (qread) then
C   ... Quick-read option requires that nr,nc are known by now.
        if (nrl*ncl == 0) then
          print
     .      '(a)', '#rdm: quick-read set but nr or nc unspecified ...'
          return
        endif
        if (lsprse) then
          ltrns = .false.
          ltrnsi = .false.
          i = 3
          if (haveC == 2) i = 4
          if (.not. pvrdm1(ifi,xwk,i)) goto 70
          i = xwk(1)
          ii = xwk(2)
          if (abs(i-xwk(1))+abs(ii-xwk(2)) /= 0 .or.
     .      i > nrl .or. ii > ncl) return
          s(i+nrl*(ii-1)) = xwk(3)
          if (wantC == 2) s(i+nrl*(ii-1)+nrl*ncl) = xwk(4)
          goto 55
        else
          if (.not.pvrdm1(ifi,s(k+1),min(haveC,wantC)*nrl*ncl-k)) return
          ltrnsi = .not. ltrns
          nelt = nrl*ncl*min(haveC,wantC)
          goto 70
        endif
      endif

C --- Parse this record for array elements; load into s(ir,ic) ---
C     NB: if nr unknown, read in s-transpose, since nc is known.
C ... Exit if we are done reading this file
      if (recrd(1:1) == '%') goto 70
      if (nl <= 0) goto 70
      ip = 0
   61 call skipbl(recrd,recl,ip)
C   ... Pick up next record when this one parsed to the end
        if (ip >= recl) goto 30
C   ... Read into s(ir,ic) (or scratch if mxelt is zero)
        i = max(nrl,1)*ic + ir + ioff
C       s(2,i,j) --- complex --- format
        if (sw10 == 4) then
          i = 2*max(nrl,1)*ic + 2*ir + ioff
        endif
        if (mxelt == 0) i = 0
        if (i > mxelt) goto 99
        if (lsprse) then
          i = 3
          if (haveC == 2) i = 4
C         Abort if insufficient number of elements
          if (a2vec(recrd,recl,ip,4,' ',1,2,i,ix,xwk) /= i) return
          i = xwk(1)
          ii = xwk(2)
C         Abort if first two entries are not integers
          if (abs(i-xwk(1))+abs(ii-xwk(2)) /= 0 .or.
     .      i > nrl .or. ii > ncl) return
          s(i+nrl*(ii-1)) = xwk(3)
          if (wantC == 2) s(i+nrl*(ii-1)+nrl*ncl) = xwk(4)
          goto 30
        else
          if (.not. a2bin(recrd,s,4,i,' ',ip,recl)) return
        endif
C        print '('' nelt,i'',2i4,''  ir,ic,ioff='',3i4,f12.6)',
C     .    nelt,i,ir,ic,ioff,s(i+1)
        nelt = nelt+1
        ic = ic+1
        if (ltrns) then
          ic = ic-1
          ir = ir+1
        endif
        if (ncl > 0 .and. ltrns) then
          if (ir > nrl) call rx('bug in rdm')
          if (ir == nrl) then
            ic = ic+1
C       ... patch exit for now ... need to read to EOF too!
            if (ncl > 0 .and. ic == ncl) then
              if (min(wantC,haveC) == 1 .or. ioff /= 0) goto 70
              ioff = ncl*nrl
              ic = 0
C             s(2,i,j) --- complex --- format
              if (sw10 == 4) then
                ioff = 1
              endif
            endif
            ir = 0
          endif
        elseif (nrl > 0 .and. .not. ltrns) then
          if (ic > ncl) call rx('bug in rdm')
          if (ic == ncl) then
            ir = ir+1
C       ... patch exit for now ... need to read to EOF too!
            if (nrl > 0 .and. ir == nrl) then
              if (min(wantC,haveC) == 1 .or. ioff /= 0) goto 70
              ioff = nrl*ncl
              ir = 0
C             s(2,i,j) --- complex --- format
              if (sw10 == 4) then
                ioff = 1
              endif
            endif
            ic = 0
          endif
        endif
      goto 61

C --- Cleanup.  If nr is zero, set nr and copy s-transpose to s ---
   70 continue
      if (ncl == 0 .and. ltrns) then
        ncl = nelt/nrl/haveC
        nc = ncl
      elseif (nrl == 0 .and. .not. ltrns) then
        nrl = nelt/ncl/haveC
        nr = nrl
        ltrnsi = mxelt /= 0
      elseif (nelt < nrl*ncl*haveC .and. .not. (qread .or. lsprse)
     .    .and. lpr) then
        call awrit4('# rdm (warning) expected %i elements (%i rows,'//
     .    ' %i cols) but read %i',' ',80,i1mach(2),nrl*ncl*haveC,nrl,
     .    ncl,nelt)
      endif
      ncx = ncl
      if (ltrnsi .and. wantC == 2) ncx = 2*ncl
      if (ltrnsi) then
C       if (nelt*2 > mxelt .and. mxelt > 0) goto 99
        ioff = nr*ncx
        allocate(swk(nr*nc))
        call dcopy(nr*nc,s,1,swk,1)
        do  71  ir = 1, nr
        do  71  ic = 1, nc
   71   s(ir+nr*(ic-1)) = swk(ic+nc*(ir-1))
        if (ncx == 2*nc) then
          k = nr*nc
          call dcopy(k,s(1+k),1,swk,1)
          do  72  ir = 1, nr
          do  72  ic = 1, nc
   72     s(ir+nr*(ic-1)+k) = swk(ic+nc*(ir-1))
        endif
        deallocate(swk)
      endif

C --- Force matrix symmetric or hermitian if so stated ---
      if (retval == 10) call dosymm(max(wantC,1),s,nr,nc)

C --- Normal exit ---
   98 rdm = max(wantC,1) + retval
      if (mod(sw/100000,10) == 1 .and. nelt < nrl*ncl*haveC
     .    .and. .not. (qread .or. lsprse)) rdm = -1
      return

C --- Exit when not enough space to read s ---
   99 continue
      print '(a)', '#rdm: insufficent space ...'
      return

C --- Exit for incompatible input ---
  196 continue
      call info5(2,0,0,'# rdm: given dimensions (%i,%i) '//
     .  'incompatible with file (%i,%i)',nr,nc,nrl,ncl,0)
      return

C --- Exit for binary read ---
  198 continue
      if (lpr) print '(a)', '#rdm: binary read failed...'
      rdm = -1

C --- General exit ---
  199 continue
      if (nr == 0) nr = nrl
      if (nc == 0) nc = ncl
      end
      logical function pvrdm1(ifi,s,ns)
C- Kernel to read array s with unformatted read
      implicit none
      integer ifi,ns
      double precision s(ns)

      pvrdm1 = .false.
      if (ns > 0) read(ifi,*,err=99,end=99) s
      pvrdm1 = .true.
   99 continue
      end
      subroutine dosymm(icast,s,nr,nc)
C- Kernel to symmetrize array s
      integer nr,nc,nn
      double precision s(nr,nc,2),xx
      integer i,j,icast

      nn = min(nr,nc)
      do  10  i = 1, nn
      do  10  j = 1, i
        xx = (s(i,j,1) + s(j,i,1))/2
        s(i,j,1) = xx
        s(j,i,1) = xx
   10 continue

      if (icast == 2) then
        do  20  i = 1, nn
        do  20  j = 1, i
          xx = (s(i,j,2) - s(j,i,2))/2
          s(i,j,2) = xx
          s(j,i,2) = -xx
   20   continue
      endif

      end
