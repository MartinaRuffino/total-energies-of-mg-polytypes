C#define DEBUG
C#define TEST
C#define unix
C#define PATH
C#define FEXTENS
      integer function fopna(fnam,lunit,ista)
C- File opening
C ----------------------------------------------------------------
Ci Inputs
Ci    fnam:  file name
Ci   lunit:  file logical lunit.  If passed lunit is -1, program finds
Ci           the next available logical unit; see Remarks.
Ci   ista:   status switches governing mode of file opening; see Remarks
Co Outputs
Co   fopna returns logical lunit number for file name
Cr Remarks
Cr   fopnx, below, is an all-purpose file-handling routine that
Cr   keeps an internal table that maintains an association of filenames
Cr   with their logical unit and status.  It is intended
Cr   to facilitate reference to files without carrying around
Cr   the logical unit number or file status.  Several 'front end'
Cr   routines listed here call fopnx with special arguments.
Cr
Cr   fopng(fnam,lunit,ista):
Cr     opens file 'fnam' with attributes set by 'ista',
Cr     assigns logical unit 'lunit', and returns the logical unit.
Cr     If lunit<0, the logical unit will be assigned if there is no
Cr     record of this file name or logical unit in the internal table
Cr     maintained by fopnx.
Cr     In this case the file name and lunit are added to the table.
Cr     If fopnx has a record of the logical unit, that unit and
Cr     the corresponding filename are used for file handling.  Else:
Cr     If fopnx has a record of the file name, that name and
Cr     corresponding lunit are used for file handling.
Cr
Cr     ista is a composite of integers a, b, c, d, stored as digits
Cr     abcd, base 2.
Cr     Bits 01: 0, open the file as 'UNKNOWN'
Cr              1, open the file as 'OLD'
Cr              2, open the file as 'NEW'
Cr     Bit 2:   1, open file as unformatted.
Cr     Bit 3:   1, do not add this file to the file table
Cr     Thus ista=5 opens file unformatted, status='old'
Cr
Cr   fopna(fnam,lunit,ista):
Cr     is identical to fopng, except that 'fnam' a path is prepended
Cr     to, and extension is appended to 'fnam'.  'path' and 'ext'
Cr     are set and recovered with:
Cr
Cr   fext(ext)     sets extension
Cr   fpath(path)   sets path.  NB: path should end with '/'
Cr   fextg(ext)    recovers extension
Cr   fpathg(ext)   recovers path
Cr
Cr   fadd(fnam,lunit,ista) adds 'fnam' to the internal table.
Cr
Cr   fopn(fnam)  is equivalent to fopna(fnam,0,-1)
Cr
Cr   fopnn(fnam)   same as fopn, but open with status='NEW'
Cr   fopno(fnam)   same as fopn, but open with status='OLD'.
Cr
Cr   fxst(fnam)    inquires as to whether file exists
Cr
Cr   ffnam(fnam,fnam) returns full file name
Cr
Cr   fhndl(fnam) returns with logical lunit associated with name,
Cr     and -1 if none exists.
Cr
Cr   fclr(fnam,unit) clears the file name from the internal table,
Cr     and closes the file if open.
Cr
Cr   ftflsh(unit) flushes output buffer of logical unit.
Cr     This routine is nonstandard, and works only some some machines.
Cr     unit = stdout calls a C-language routine that flushes stdout
Cr     unit < 0 calls a C-language routine that flushes all buffers
Cr
Cr   fshow lists the internal table
Cu Updates
Cu   15 Feb 02 (ATP) added MPI log for lgunit(3)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) fnam
      integer ista,lunit
C Local parameters
      integer fext,fextg,fopn,fopno,fopnn,fhndl,fadd,fxst,fpath,fpathg,fopng,fopnx
      integer mode,bitand

C#ifndef AIX
      mode = 100 + 70*mod(ista/8,2)
      fopna = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
      return

      entry fopng(fnam,lunit,ista)
      mode =   2 + 70*mod(ista/8,2)
      fopng = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
      return

      entry fadd(fnam,lunit,ista)
      mode = 114 + 70*mod(ista/8,2)
      if (lunit < 0) mode = mode-4
      fadd = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
      return

      entry fopn(fnam)
      fopn = fopnx(fnam,100,0,-1)
      return

      entry fopno(fnam)
      fopno = fopnx(fnam,100,8+1,-1)
      return

      entry fopnn(fnam)
      fopnn = fopnx(fnam,100,8+2,-1)
      return

      entry fxst(fnam)
      fxst = fopnx(fnam,170,-1,-1)
      return

      entry fhndl(fnam)
      fhndl = fopnx(fnam,184,-1,-1)
      return

      entry fext(fnam)
      fext = fopnx(fnam,120,0,-1)
      return

      entry fpath(fnam)
      fpath = fopnx(fnam,140,0,-1)
      return

      entry fextg(fnam)
      fextg = fopnx(fnam,130,0,-1)
      return

      entry fpathg(fnam)
      fpathg = fopnx(fnam,150,0,-1)
      return
      end

      subroutine fclr(fnam,lunit)
      implicit none
C Passed parameters
      character*(*) fnam
      integer lunit,fopnx,i

      i = fopnx(fnam,111,0,lunit)
      return

      entry fshow
      i = fopnx(' ',190,0,-1)
      end
      subroutine ffnam(fnam,filnam)
C- Appends extension to fnam
      implicit none
      character*(*) fnam,filnam
      integer i,fopnx
      filnam = fnam
      i = fopnx(filnam,160,0,-1)
      end

C#elseC
C      mode = 100 + 70*mod(ista/8,2)
C      fopna = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
C      end
C      integer function fopng(fnam,lunit,ista)
CC      implicit none
C      character*(*) fnam
C      integer mode,lunit,fopnx,ista
C      mode =   2 + 70*mod(ista/8,2)
C      fopng = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
C      end
C      integer function fadd(fnam,lunit,ista)
CC      implicit none
C      character*(*) fnam
C      integer ista,lunit
C      integer fopnx,mode
C      mode = 114 + 70*mod(ista/8,2)
C      if (lunit < 0) mode = mode-4
C      fadd = fopnx(fnam,mode,mod(ista,8)+bitand(ista,32),lunit)
C      end
C      integer function fopn(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fopn = fopnx(fnam,100,0,-1)
C      end
C      integer function fopno(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fopno = fopnx(fnam,100,8+1,-1)
C      end
C      integer function fopnn(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fopnn = fopnx(fnam,100,8+2,-1)
C      end
C      integer function fxst(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fxst = fopnx(fnam,170,-1,-1)
C      end
C      integer function fhndl(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fhndl = fopnx(fnam,184,-1,-1)
C      end
C      integer function fext(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fext = fopnx(fnam,120,0,-1)
C      end
C      integer function fpath(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fpath = fopnx(fnam,140,0,-1)
C      end
C      integer function fextg(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fextg = fopnx(fnam,130,0,-1)
C      end
C      integer function fpathg(fnam)
CC      implicit none
C      character*(*) fnam
C      integer fopnx
C      fpathg = fopnx(fnam,150,0,-1)
C      end
C      subroutine fclr(fnam,lunit)
CC      implicit none
CC Passed parameters
C      character*(*) fnam
C      integer lunit,fopnx,i
C
C      i = fopnx(fnam,111,0,lunit)
C      return
C
C      entry fshow
C      i = fopnx(' ',190,0,-1)
C      end
C      subroutine ffnam(fnam,filnam)
CC- Appends extension to fnam
CC      implicit none
C      character*(*) fnam,filnam
C      integer i,fopnx
C      filnam = fnam
C      i = fopnx(filnam,160,0,-1)
C      end
C#endif
      integer function fopnT(fnam,lunit,ista,mode)
C- File opening, using a separate directory
C ----------------------------------------------------------------------
Ci Inputs
Ci   fnam  :file name WITHOUT path or extension
Ci   lunit :file logical lunit.  If passed lunit is -1, program finds
Ci         :the next available logical unit; see Remarks in fopna
Ci   ista  :status switches governing mode of file opening;
Ci         :see Remarks in fopna
Ci   mode  :selects which path is used for file opening
Ci         :0 fopnT opens file using tmpdir for path.
Ci         :  The extension used by fopna is appended
Ci         :1 fopnT opens file using tmpdir for path.
Ci         :  No extension is appended.
Ci         :..Add 2 to cause fopnT returns logical unit associated
Ci         :  with name if it exists in table, but do not touch file
Ci         :11 fopnT opens no files, but copies fnam into tmpdir
Ci         :  In this case, lunit and ista are not used.
Ci         :12 fopnT opens no files, but gets environment variable
Ci         :  TMPDIR and copies it into tmpdir
Ci         :  In this case, fnam, lunit and ista are not used.
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   30 Apr 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character fnam*(*)
      integer ista,lunit,mode
C ... Local parameters
      integer j1,j2,fextg,fopng,fopnx
      character tmpdir*100,pfnam*100,ext*100
      save tmpdir
      data tmpdir /' '/

      if (mode == 12) then
        call gtenv('TMPDIR',tmpdir)
      elseif (mode == 11) then
        tmpdir = fnam
      endif
      if (mode == 11 .or. mode == 12) then
        call word(tmpdir,1,j1,j2)
        if (j2 >= j1) then
          if (tmpdir(j2:j2) /= '/') tmpdir(j2+1:j2+1) = '/'
        endif
        fopnT = 0
        return
      endif

      call isanrg(mode,0,4,'fopnt','mode',.true.)

      if (mod(mode,2) == 0) then
        j1 = fextg(ext)
        call namcat(pfnam,fnam,ext,.false.)
        ext = pfnam
        call namcat(pfnam,tmpdir,ext,.true.)
      else
        call namcat(pfnam,tmpdir,fnam,.true.)
      endif

      if (mode <= 1) then
        fopnT = fopng(pfnam,lunit,ista)
      else
        fopnT = fopnx(pfnam,186,-1,-1)
      endif
      end
      integer function fopnx(pnam,mode,ista,lunit)
C- Kernel called by fopna, etc to open files
C ----------------------------------------------------------------
Ci Inputs
Ci   pnam:   file name, or portion of file name, depending on mode
Ci   mode:   1s digit (for 10s digit 0)
Ci           0  Find lunit, ista from name.
Ci              pnam is name sans path, extension
Ci              If name does not exist, add to table, assign lunit.
Ci              Subsequent file handling depends on ista.
Ci           1  Remove name from internal table.  Close file if open.
Ci              Error if name is missing from table.
Ci          ... Add 2 if pnam is already full name
Ci          ... Add 4 to require lunit match internal table
Ci                    returns -1 if lgunit<0 and pnam missing from table
Ci              Special meanings when used in conjunction with 10s digit 9
Ci              2 => Find name associated with lunit and return in pnam (full name)
Ci              4 => Find name associated with lunit and return in pnam (sans extension)
Ci          10s digit: for modifying or retrieving file's path and extension
Ci           ... 1..6 If resulting file name (after addition of path and
Ci                     extension) is new, it is added to internal table.
Ci           ... 1..6,8,9 suppress file handling
Ci           1  Suppress any file handling
Ci           2  copy pnam to extension
Ci           3  copy extension into pnam
Ci           4  copy pnam to path
Ci           5  copy path into pnam
Ci           6  copy full file name into pnam.
Ci           7  Do not add file name to internal table, but allow
Ci              removal of name from table
Ci           8  Do not add file name to internal table
Ci           9  AND 1s digit mode=0 : show internal tables; no file handling
Ci           9  AND 1s digit mode=4 : return full name associated with lunit
Ci          100s digit
Ci           1  convert pnam to lower case
Ci   ista:  -1  return 1 if file exists, 0 otherwise
Ci         >=0  sets flags file file status on opening.  This status
Ci              is kept and used along with the file table.  Internal
Ci              values are used except when missing or preempted:
Ci           ...The following bits work as independent switches:
Ci           ...bits 0,1:
Ci           0  open the file as 'UNKNOWN', if not open already
Ci           1  open the file as 'OLD'
Ci           2  open the file as 'NEW'
Ci           ...bit 2:
Ci           4  open the file as 'UNFORMATTED'
Ci           ...bits 3,4
Ci           8  Override internal ista for bits 0,1
Ci          16  Override internal ista for bit 2
Ci           ...bit 5
Ci          32  Open file with position='append'
Ci   lunit:  (lunit>0) if in table, uses lunit to match instead of name
Ci              UNLESS mode, 10s digit set.
Cb Bugs
Cb   This routine mixes up many different functions.
Cb   Operations handling internal tables of file names and logical units
Cb   should be separately packaged from those accessing (opening or
Cb   inquiring about) files.
Cr Remarks
Cr   mxnam is maximum number of files that may be open at one time.
Cr   All file closings should be made through entry fclose
Cu Updates
Cu   01 May 15 Added modes 92 and 94
Cr   17 Jun 97 internal ista bits 0,1 ALWAYS 0,
Cr             making default file opening 'unknown'
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ista,lunit,mode
      character*(*) pnam
      character*11 ftnfmt
C ... Local parameters
      logical bittst,isopen,ldum,mod100
      integer bit,bitand,getdig,i,i1,i2,iprint,locuni,lsta,matchn,
     .  matchu,mode1,mode10,n,IDBG,isw
      character ftnsta*7, outs*80
C Parameters concerning limits to string size in filenames:
C plen: path length; extlen: ext length; maxsiz: name (including ext)
      integer maxsiz,extlen,plen
C#ifdef unix
      parameter (extlen=20, plen=80, maxsiz=100)
C#ifdef DEBUG
      parameter (IDBG=10)
C#elseC
C      parameter (IDBG=110)
C#endif
C Save list of filenames
      integer mxnam
      parameter (mxnam=30)
      character*(maxsiz) tabnam(mxnam),pfnam,fnam
      character*(extlen) ext
      character*(plen) path
      integer tabsw(mxnam), tabuni(mxnam),nnam,lgunit
      save tabnam, tabsw, tabuni, nnam, ext, path
      data nnam /0/
C#ifdef FEXTENS
      data ext /'.dat'/
C#elseC
C      data ext /' '/
C#endif
      data path /' '/

      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) == bit)

      mode1   = getdig(mode,0,10)
      mode10  = getdig(mode,1,10)
      mod100 = getdig(mode,2,10) /= 0

C --- Cases handling partial file names ---
      fopnx = 0
C ... Change extension
      if (mode10 == 2) then
        call namcat(ext,pnam,' ',mod100)
        if (iprint() >= IDBG) call awrit0(
     .    ' FOPNX: set ext to "'//ext//'%a"',' ',80,lgunit(1))
        return
C ... Copy extension into pnam
      elseif (mode10 == 3) then
        call namcat(pnam,ext,' ',.false.)
        return
C ... Change internal path name
      elseif (mode10 == 4) then
        call word(pnam,1,i1,i2)
        if (i2 <= 0) then
          path = ' '
C        elseif (pnam(i2:i2) == '/') then
C          call namcat(path,pnam,' ',mod100)
C        else
C          call namcat(path,pnam,'/',mod100)
        else
          call namcat(path,pnam,' ',mod100)
        endif
        if (iprint() >= IDBG) call awrit0(
     .    ' FOPNX: set path to "'//path//'%a"',' ',80,lgunit(1))
        return
C ... Copy path into pnam
      elseif (mode10 == 5) then
        call namcat(pnam,path,' ',.false.)
        return
C ... Return name associated with lunit
      elseif (mode10 == 9 .and. (mode1 == 2 .or. mode1 == 4)) then
        pnam = ' '
        do  i = 1, nnam
          if (tabuni(i) /= lunit) cycle
          i1 = 0
          if (mode1 == 4) i1 = index(tabnam(i),ext,.true.)
          if (i1 > 1) then
            pnam = tabnam(i)(1:i1-1)
          else
            pnam = tabnam(i)
          endif
          fopnx = 0
          return
        enddo
        fopnx = -1
        return
C ... Show internal table
      elseif (mode10 == 9) then
        call awrit3('%N fopen:  %i files in table.'//
     .    '%?;n;  path='''//path//'%a'';;'//
     .    '%?;n;  ext='''//ext//'%a'';;',' ',-80,lgunit(1),nnam,
     .    isw(path /= ' '),isw(ext /= ' '))
        n = 0
        do  10  i = 1, nnam
          call word(tabnam(i),1,i1,i2)
          n = max(n,i2)
   10   continue
        n = n+2
        call awrit1(' file%npunit   mode  open',' ',80,lgunit(1),n)
        call pshpr(1)
        do  i = 1, nnam
          call awrit4(' '//tabnam(i)//
     .      '%np%,2i    %?;n==0;  fmt;unfmt;   %l',' ',80,lgunit(1),n,
     .      tabuni(i),bitand(tabsw(i),4),isopen(tabuni(i),.false.))
        enddo
        call poppr
        print *, ' '
        return

      endif

C --- Make full file name ---
      if (bittst(mode1,2)) then
        call namcat(pfnam,pnam,' ',mod100)
      else
        call namcat(pfnam,pnam,' ',mod100)
        call namcat(fnam,pfnam,ext,.false.)
        call namcat(pfnam,path,fnam,.false.)
      endif

C ... Return full file name in pnam
      if (mode10 == 6) then
        call namcat(pnam,pfnam,' ',.false.)
        return
      endif

C --- Match lunit with table.  matchu=-1 unless a match exists ---
      matchu = -1
      if (lunit > 0) then
        do  110  i = 1, nnam
          matchu = i
  110   if (tabuni(i) == lunit) goto 111
        matchu = -1
  111   continue
      endif

C --- If no match, get next available unit number => locuni ---
      locuni = lunit
      if (matchu == -1 .and. lunit < 0) then
        locuni = 7
   12   locuni = locuni+1
        if (isopen(locuni,.false.)) goto 12
        do  14  i = 1, nnam
   14   if (tabuni(i) == locuni) goto 12
      endif

C --- Match pfnam with table, unless not needed ---
      matchn = matchu
      if (matchu <= 0 .or. bittst(mode1,4)) then
        do  30  i = 1, nnam
        matchn = i
   30   if (tabnam(i) == pfnam) goto 31
        matchn = -1
   31   continue
      endif
      if (lunit < 0 .and. matchn > 0) then
        matchu = matchn
        locuni = tabuni(matchn)
      elseif (lunit < 0 .and. bittst(mode1,4)) then
C        call rxs('fopnx: no unit or name specified, file ',pnam)
        fopnx = -1
        return
      endif

      fopnx = locuni

C ... Force match between name and unit, or error
      call awrit1('%xFOPNX: file '//pfnam//'%a, unit %i, '//
     .  'conflicts with prior use',outs,80,0,locuni)
      if (matchu /= matchn .and. bittst(mode1,4)) then
        if (matchu > 0)
     .    call awrit0('%a%10b file '//tabnam(matchu),outs,80,0)
        call rx(outs)
      endif

C --- Add a file to internal table ---
      if (matchu < 0 .and. mode10 <= 6) then
        if (mod(mode1,2) == 1) then
          if (lunit > 0)
     .    call rxi('FOPNX: attempt to remove nonexistent unit',lunit)
          call rxs('FOPNX: attempt to remove nonexistent name ',pfnam)
        endif
        nnam = nnam+1
        if (nnam > mxnam) call
     .    rxi('FOPNX: too many file names.  max=',mxnam)
        tabnam(nnam) = pfnam
        tabuni(nnam) = locuni
        tabsw(nnam)  = mod(ista,8) - mod(ista,4)
        matchu = nnam
C --- Remove a file entry from the internal table ---
      elseif (mod(mode1,2) == 1 .and. mode10 <= 7) then
        i = tabuni(matchu)
        if (isopen(i,.false.)) call fclose(i)
        do  130  i = matchu+1, nnam
          tabnam(i-1) = tabnam(i)
          tabuni(i-1) = tabuni(i)
          tabsw(i-1)  = tabsw(i)
  130   continue
        nnam = nnam-1
        return
      endif

C#ifdef DEBUG
      if (iprint() >= IDBG) then
        call pshpr(1)
        call awrit3(' FOPEN: file '//pfnam//'%a unit %i  status %i'//
     .  ' open=%l',' ',80,lgunit(1),locuni,ista,isopen(locuni,.false.))
        call poppr
      endif
C#endif

C --- Cases ista=-2,-1 ---
      if (ista == -2 .or. mode10 /= 0 .and. mode10 /= 7) return
      if (ista == -1) then
        INQUIRE(FILE=pfnam,EXIST=ldum)
        fopnx = 0
        if (ldum) fopnx = 1
        return
      else
      endif

C ... Format status
      lsta = ista
      if (matchu > 0) then
        lsta = bitand(lsta,1023-4) + bitand(tabsw(matchu),4)
        if (bittst(ista,16)) then
          lsta = bitand(lsta,1023-4) + bitand(ista,4)
          tabsw(matchu) = bitand(tabsw(matchu),1023-4) + bitand(ista,4)
        endif
      endif
      ftnfmt = 'FORMATTED'
      if (bittst(lsta,4)) ftnfmt = 'UNFORMATTED'

C ... Handle fortran status (slow, but apparently stable)
      ftnsta = 'UNKNOWN'
      if (mod(lsta,4) == 1) ftnsta = 'OLD'
      if (mod(lsta,4) == 2) then
        ftnsta = 'NEW'
        call fclose(locuni)
        open(locuni,FILE=pfnam,FORM=ftnfmt,STATUS='UNKNOWN')
        close(unit=locuni,status='DELETE')
      endif

      if (.not. isopen(locuni,.true.)) then
        call word(pfnam,1,i1,i2)
        call chrpos(pfnam,' ',plen+maxsiz,i)
        if (bitand(ista,32) == 0) then
          if (iprint() >= IDBG) print '(a,a,a,a,a,a,a,i3)',' FOPEN: opening ',
     .    ftnfmt," file '",pfnam(1:i2),"', status=",ftnsta,', unit=',locuni
          open(locuni,FILE=pfnam(i1:i2),FORM=ftnfmt,STATUS=ftnsta)
          rewind locuni
        else
          if (iprint() >= IDBG) print '(a,a,a,a,a,a,a,a,i3)',' FOPEN: opening ',ftnfmt,
     .      " file '",pfnam(1:i2),"', status=",ftnsta,", position=append",', unit=',locuni
          open(locuni,FILE=pfnam(i1:i2),FORM=ftnfmt,STATUS=ftnsta,POSITION='append')
        endif
      endif

      end
      logical function isopen(unit,sw)
C- Returns whether unit aleady open; optionally adds unit to list if not
C ----------------------------------------------------------------
Ci Inputs
Ci   unit: if unit not among list of open files, add to list
Ci   sw:   if false, return isopen but do not add to list
Co Outputs
Co   isopen
Cr Remarks
Cr   maxfil is maximum number of files that may be open at one time
Cr   All file closings should be made through entry fclose
C ----------------------------------------------------------------
      implicit none
      integer unit
      logical sw
      integer i,iprint
      integer maxfil
      parameter (maxfil=15)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      if (iprint() > 110) print 20, unit, (unitab(i), i=0, nopen-1)
   20 format(/' ISOPEN: check logical unit',i3,
     .        ' among open units: ',15i3)

      isopen = .true.
      do  10  i = nopen-1, 0, -1
   10 if (unitab(i) == unit) return

      isopen = .false.
      if (.not. sw) return
      unitab(nopen) = unit
      nopen = nopen+1
      if (nopen > maxfil) call rx('ISOPEN: too many files')
      end
      subroutine fclose(unit)
C- Closes an open file, removing unit from the stack
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   none
Cr Remarks
Cr   use in conjunction with isopen
C ----------------------------------------------------------------
      implicit none
      integer unit
      integer i,iprint
      character*6 clstat
      integer maxfil
      parameter (maxfil=15)
      integer unitab(0:maxfil-1),nopen
      common /funits/ unitab,nopen

      clstat = 'KEEP'
      goto 10

   10 continue
      if (iprint() >= 100) print 20, unit
   20 format(' FCLOSE: closing',i3)

      do  30  i = nopen-1, 0, -1
   30 if (unitab(i) == unit) unitab(i) = 999
      call ishell(nopen,unitab)
      if (nopen == 0 .or. unitab(max(nopen-1,0)) /= 999) then
        if (iprint() >= 100)
     .  print *, 'FCLOSE: attempt to close unopened file'
      else
        close(unit=unit,status=clstat)
        nopen = nopen-1
      endif
      return

      entry dfclos(unit)
      clstat = 'DELETE'
      goto 10

C --- Close all open files ---
      entry closea
      do  50  i = 0, nopen-1
        if (iprint() >= 100) print *,'CLOSEA: closing file ',unitab(i)
        close(unit=unitab(i),status='KEEP')
        unitab(i) = 999
   50 continue
      nopen = 0
      return
      end
      integer function lgunit(i)
C- Returns stdout for i=1, log for i=2, mlog for i=3 (MPI logfile)
      implicit none
      integer i, fopn, i1mach, fhndl
      lgunit = i1mach(2)
      if (i == 1) return
      if (i == 2) then
        lgunit = -1
        if (fhndl('LOG') < 0) return
        lgunit = fopn('LOG')
      elseif (i == 3) then
        lgunit = -1
        if (fhndl('MLOG') < 0) return
        lgunit = fopn('MLOG')
      endif
      end
      subroutine namcat(dest,fnam1,fnam2,lcase)
C- Concatenates two names, optionally lowering the case
      implicit none
      logical lcase
      character*(*) fnam1,fnam2,dest
      integer iprint,i1,i2,j1,j2

      call word(fnam1,1,i1,i2)
      call word(fnam2,1,j1,j2)
      i2 = max(i1-1,i2)
      j2 = max(j1,j2)
      if (i2 >= i1) then
        dest = fnam1(i1:i2)//fnam2(j1:j2)
      else
        dest = fnam2(j1:j2)
      endif
      if (i2-i1+j2-j1+2 > len(dest) .and. iprint() > 0) then
        print *, 'namcat: string "',fnam1(i1:i2),fnam2(j1:j2),
     .    '" was truncated to "',dest,'"'
      endif

C#ifdef unix
      if (lcase) call locase(dest)
C#endif

      end
      subroutine ftflsh(ifi)
C- Flushes output buffer for file ifi
C  NB: this routine is nonstandard, and machine-dependent.
      implicit none
      integer ifi
      integer i1mach
C#ifdefC AIX
C      integer j1,j2
C      character fname*101,fmt*33
C      logical open
C#endif

      if (ifi == i1mach(2)) then
        call flushs(1)
      elseif (ifi < 0) then
        call flushs(-1)
C#ifdefC AIX
C      else
C        inquire (UNIT=ifi,OPENED=open,NAME=fname,FORM=fmt)
C        call word(fname,1,j1,j2)
CC        write (6,220) ifi,open,fname(j1:j2),fmt
CC  220   format(' unit',i3,' open=',l1,'  name=',a,'  fmt=',a)
C        if (open) then
C          close(ifi)
C          open(ifi, file=fname, form=fmt, status='old')
C        endif
C#endif
      endif
      end
C#ifdef TEST
      subroutine fmain
C      implicit none
      integer fopna,fext,fextg,fopn,fopno,fopnn,fhndl,fadd,fxst,fpath,
     .  fopng,fopnx,fopnT
      integer ifi,jfi,kfi,lgunit,i
      double precision s(2,4)
      character fnam*20,strn*120
*      call pshpr(110)

C     Check fopnT
C      i = fopnT(' ' ,-1,0,12)
C      i = fopnT('tmp' ,-1,0,0)
C      ifi = fopn('tmp')
C      call fshow
C      i = fhndl('tmp')
C      print *, i
C      call dfclos(fopn('TMP'))
C      i = fopnT('tmp' ,-1,0,2)
C      print *, i
C      call fshow
C      i = fhndl('tmp')
C      print *, i
C      call rx('done')

      s(1,1) = 11
      s(1,2) = 12
      s(2,1) = 21
      s(2,2) = 22

      fnam = 'xyz.123'
      ifi = -1
      print "(' ... test fopng(fnam=''xyz.123'',-1,4)')"

      ifi = fopng(fnam,-1,4)
      call ywrm(1,'testing fopng',1,ifi,' ',s,2,2,2,2)
      call fshow
      print "(' ... test fclose')"
      call fclose(ifi)
      call fshow

      ifi = 10
      print "(' ... test fopna(fnam=''pqr'',ifi,0), ifi=',i3)", ifi
      ifi = fopna('pqr' ,ifi,0)
      call fshow
      call dscal(4,2d0,s,1)
      call ywrm(0,' testing fopng',1,ifi,'(5f12.2)',s,2,2,2,2)

      ifi = fhndl('xyz')
      jfi = fhndl('pqr')
      print "(/' ... test fhndl(''xyz.123''),fhndl(''pqr.dat'') : ',2i4)", ifi,jfi

      print "(/' ... Close, reopen file ''pqr.dat'' with position=append, ifi='i4)", jfi
      call fclose(ifi)
      call fshow
      ifi = fopna('pqr' ,ifi,32)
      call fshow
      print "(' ... Append ''this is last line'' to file ''pqr.dat''')"
      write(ifi,'(a)') 'this is the last line'

      print "(' ... Close, reopen file ''pqr.dat'' with rewind, ifi='i4)", jfi
      call fclose(jfi)
      ifi = fopna('pqr' ,ifi,0)
      print "(/' ... Print out file ''pqr.dat''')"
      do i = 1,4
        read(ifi,'(a)') strn
        print '(a)', trim(strn)
      enddo

      print "(/' ... retest fopng(fnam=''xyz.123'',-1,4)')"
      ifi = fopng(fnam,-1,4)
      call fshow
      print "(' ... read array from binary file xyz.123, ifi=',i4)", ifi
      call rdm(ifi,0002,8,' ',s,2,2)
      call ywrm(0,' testing fopng',1,lgunit(1),'(5f12.2)',s,2,2,2,2)

      print "(/' ... test fopng(fnam=''xyz.123'',-1,4), file already open')"
      ifi = fopng(fnam,-1,4)
      call fshow

      print "(' ... test fadd(''abc'',ifi,0)')"
      ifi = 14
      ifi = fadd('abc',ifi,0)
      call fshow

      print "(' ... call fext(''.dat2'')')"
      ifi = fext('.dat2')
      call fshow

      print "(' ... test fadd(''abc'',ifi,0) with different extension')"
      ifi = fadd('abc',-1,0)
      call fshow
      print "(' ... test fopn(''abc'')')"
      ifi = fopn('abc')
      call ywrm(0,' testing fopn',1,ifi,'(5f12.2)',s,2,2,2,2)
      ifi = fopn('abc')
      call ywrm(0,' testing fopn',1,ifi,'(5f12.2)',s,2,2,2,2)
      call fshow
      jfi = fopna('abc',-1,0)
      if (ifi /= jfi) stop 'bug in fopna'

C      call awrit0('%N ... test ftflsh (stdout, abc.dat2)',' ',80,lgunit(1))
C      write(*,"(' ... testing ftflsh(6) ... ')",advance='no')
C      call ftflsh(6)
C      write(*,*) 'after ftflsh(6)'
C      write(*,*) '... calling ftflsh(-1)'
C      write(ifi,*) '... before ftflsh(-1)'
C      call ftflsh(-1)
C      write(*,*) '... calling ftflsh(ifi)'
C      write(ifi,*) '... before ftflsh(ifi)'
C      call ftflsh(ifi)
C      write(ifi,*) '... after ftflsh(ifi)'
C
C      print *, '... invoking system  call ''cat abc.dat2'' ...'
C      call ftflsh(6)
C      call fsystm('cat abc.dat2',i)
C      call awrit1('%N ... system  call ''cat abc.dat2'' returned... %i',' ',80,lgunit(1),i)

      call awrit0('%N ... restore ext to dat',' ',80,lgunit(1))
      ifi = fext('.dat')

      print "(/' ... test fclr(''abc'',-1) --- refer to file by name')"
      call fclr('abc',-1)
      call fshow

      print "(' ... test fclr('' '',10) -- refer to file ''pqr'' by logical unit')"
      call fclr(' ',10)
      call fshow

      print "(/' ... test fopna(''xyz'',-1,8) to show internal file table not affected')"
      ifi = fopna('xyz',-1,8)
      write (ifi,*) 'test'
      call fclose(ifi)
      call fshow

      print "(/' ... reopen with fopna(''xyz'',-1,0) to show modified internal file table')"
      ifi = fopna('xyz',-1,0)
      read(ifi,'(a10)') fnam
      print *, 'after opening, read line from file (should be ''test''): ',trim(fnam)
      ifi = fxst('xyz')
      jfi = fxst('zzz')
      kfi = fopnx('xyz.123',172,-1,-1)
      call awrit3(' ifi for xyz is %i, for zzz is %i, for xyz.123 %i',' ',80,lgunit(1),ifi,jfi,kfi)
      call fshow

      call awrit0('%N ... show fopn works without fadd',' ',80,lgunit(1))
      ifi = fopn('TMP')
      call fshow

      call awrit0('%N ... call fopnx to test locase, "NEW","OLD"',' ',80,lgunit(1))
      ifi = fopnx('TMP',2,2,-1)
      call fclose(ifi)
      call fshow
      ifi = fopnx('TMP',2,8+1,-1)
      call fshow

      call awrit0('%N ... use path as pre-extension',' ',80,lgunit(1))
      ifi = fext(' ')
      ifi = fpath('dat.')

      ifi = fopn('tmp')
      call fshow
      call fclr('TMP',-1)
      call fshow

      call awrit0('%N ... check opening as NEW',' ',80,lgunit(1))
      ifi = fopnx('STR',100,16+8+4+2,-1)
      call fshow
      call fclose(ifi)

      call awrit0('%N ... check fopn, UNFORMATTED',' ',80,lgunit(1))
      ifi = fopn('STR')
      call fshow

      call closea
      call fshow

      end
C#endif
