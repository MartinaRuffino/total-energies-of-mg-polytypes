      subroutine logwarn(opt,strn)
C- Print out warning and increment warning log
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   opt   : <0 reset warnings counter to zero
Ci         :  0 return number of warnings logged in opt
Ci         :    Note: opt is modified on output
Ci         : >0 increment warnings counter.
Ci         :    strn is printed with verbosity check given by opt
Ci         :    if its length is nonzero.
Ci   strn
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr    It is an error to call logwarn with opt>=0 unless the global
Cr    variable 'warning' has been previously set.
Cu Updates
Cu   13 May 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt
      character *(*) strn
C ... Local parameters
      integer nwarn
      data nwarn /0/
      save nwarn

      if (opt < 0) then
        nwarn = 0
      elseif (opt > 0) then
        nwarn = nwarn + 1
        if (len(strn) > 0) call info0(opt,0,0,strn)
      elseif (opt == 0) then
        opt = nwarn
      endif

      end
