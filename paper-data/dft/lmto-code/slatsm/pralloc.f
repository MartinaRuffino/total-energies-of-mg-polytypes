      subroutine pralloc(subnam,varnam,n,ncut,ipr,cast)
C- Print to standard out the size of a memory allocation
C ----------------------------------------------------------------------
Ci Inputs
Ci   subnam,
Ci   varnam: output is "SUBNAM: allocate # MB for VARNAM"
Ci   n     :number of elements to allocate (note: n is double)
Ci   ipr   :print if verbosity exceeds ipr, and if size exceeds ncut
Ci   ncut  :print if verbosity exceeds ipr, and if size exceeds ncut
Ci         :(ncut in MB)
Ci   cast  : 0 integer
Ci           1 double precision
Ci           2 double complex
Ci           3 double complex
Ci           4 double complex
Cu Updates
Cu   29 Sep 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer cast,ipr,ncut
      double precision n
      character subnam*(*),varnam*(*)
C ... Local parameters
      integer size

      size = 4
      if (cast > 0) size = 8
      if (cast > 1) size = 16

      if (n*size > 1d6*dble(ncut)) then
        call info2(ipr,0,0,' '//trim(subnam)//
     .    ':  allocate %;0d MB for '//trim(varnam),(n/1d6)*size,0)
      endif

      end
