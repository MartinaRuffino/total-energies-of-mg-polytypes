      subroutine huntst(slist,n,strn,low)
C- Brackets a string within an ordered array of strings
C ----------------------------------------------------------------
Ci Inputs
Ci    slist,n: array of points and number
Ci    strn: value to bracket
Ci    low: initial guess for output low
Co Outputs
Co    low: slist(low) <= strn < slist(low+1)
Cr Remarks
Cr     Adapted from Numerical Recipes
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer n,low
      character*(*) slist(1),strn
C Local variables
      integer inc,jhi,jm
      logical ascnd

      ascnd = slist(n) > slist(1)
      if (low <= 0  .or. low > n) then
        low = 0
        jhi = n+1
        goto 3
      endif
      inc = 1
      if (strn >= slist(low) .eqv. ascnd) then
    1   jhi = low+inc
        if (jhi .gt . n)  then
          jhi = n+1
        else if (strn >= slist(jhi) .eqv. ascnd) then
          low = jhi
          inc = inc+inc
          goto 1
        endif
      else
        jhi = low
    2   low = jhi-inc
        if (low < 1) then
          low = 0
        else if (strn < slist(low) .eqv. ascnd) then
          jhi = low
          inc = inc+inc
          goto 2
        endif
      endif
    3 if (jhi-low == 1) then
        if (n <= 0) return
C ... Increment low if slist(low)<>strn and slist(low+1)=strn
        if (low < 1 .and. strn == slist(1)) low = 1
        if (low < 1 .or. low == n) return
        if (strn > slist(low) .and. strn == slist(low+1)) low=low+1
        return
      endif
      jm = (jhi+low)/2
      if (strn > slist(jm) .eqv. ascnd) then
        low = jm
      else
        jhi = jm
      endif
      goto 3
      end
