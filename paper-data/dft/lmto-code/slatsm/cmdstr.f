      logical function cmdstr(iarg,a)
C- Returns generalized command-line argument iarg
C ----------------------------------------------------------------------
Ci Inputs
Ci   iarg  :return generalized command-line argument i
Co Outputs
Co   a     :argument returned in string a
Co  cmdstr :returns false if iarg is outside the range of arguments
Cr Remarks
Cr   The usual command-line arguments may be enlarged; see cmdopt.f.
Cr   The functionality of cmdstr is similar to that of cmdopt, but
Cr   with a simpler interface.
Cu Updates
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------------
      implicit none
      character*(*) a
      integer iarg,nargf,n1,n2

      cmdstr = .true.
      n1 = nargf()
      if (iarg < n1) then
        call getarf(iarg,a)
      else
        call ncmdop(n2)
        if (iarg < n1+n2) then
          call gcmdop(iarg-n1+1,a)
        else
          cmdstr = .false.
        endif
      endif
      end
      integer function nxargc()
C- Number of generalized command-line arguments
      implicit none
      integer n,nargf

      call ncmdop(n)
      nxargc = nargf() + n
      end

      logical function cmdstrsyv(iarg,a)
C- Returns command-line argument, with variable subsitution
C ----------------------------------------------------------------------
Ci Inputs
Ci   iarg  :which command-line argument to return
Ci   a     :argument
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   14 Aug 17
C ----------------------------------------------------------------------
      implicit none
      character*(*) a
      logical ltmp
      integer iarg,recl
      character(len=8) :: cch = ' {}     '
      procedure(logical) :: cmdstr

      recl = len(a)
      ltmp = cmdstr(iarg,a)
      cmdstrsyv = ltmp
      if (.not. ltmp) return
      call substrsyv(a,recl,cch)
      end

      subroutine substrsyv(a,recl,cch)
C- Parses a string, making substitutions according to cch
      implicit none
      integer recl
      character(len=*) :: a,cch

      integer nrec,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=1,lstsiz=200,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),vnam(mxlev)*16
C     character(len=8) :: rdarg = ' {}     '
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),nlist(0:mxlev)

      nrec = 0
      call rdfiln(0,cch,mxlev,loop0,nlin,list,lstsiz,
     .  ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nrec)

      end
