      subroutine hunti(xa,n,x,iprm,low)
C- Brackets a value within an ordered array of integer points
C ----------------------------------------------------------------
Ci Inputs
Ci   xa  :array of points
Ci   n   :size of array
Ci   x   :value to bracket
Ci   iprm:permutation table by which array xa is ordered
Ci        iprm(1) <= 0 => assumes iprm(i) = i; iprm not referenced
Ci   low : initial guess for output low
Co Outputs
Co   low : xa(low) < x <= xa(low+1)
Co       : when x cannot be bracketed,
Co       : low = 0 if x<=xa(1)
Co       : low = n if xa(n)<x
Co   ... if xa is ordered in descending order:
Co   low : xa(low) > x >= xa(low+1)
Co       : when x cannot be bracketed,
Co       : low = 0 if x>xa(1)
Co       : low = n if xa(n)>=x
Cu Updates
Cu   13 May 10 Handle special case xa(1) == xa(n)
Cu   13 Sep 01 Handle case n=0
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer n,low,iprm(n)
      integer xa(n),x
C Local variables
      integer inc,jhi,jm
      logical ascnd,liprm
      integer xn

      if (n == 0) then
        low = 0
        return
      endif
C ... Pathological case: all points equal (treat as ascending)
      if (xa(n) == xa(1)) then
        if (xa(1) >= x) then
          low = 0
        else
          low = n
        endif
        return
      endif

      liprm = iprm(1) > 0
      if (liprm) then
        ascnd = xa(iprm(n)) > xa(iprm(1))
      else
        ascnd = xa(n) > xa(1)
      endif
      if (low <= 0 .or. low > n) then
        low = 0
        jhi = n+1
        goto 3
      endif
      inc = 1
      if (liprm) then
        xn = xa(iprm(low))
      else
        xn = xa(low)
      endif
      if (x >= xn .eqv. ascnd) then
    1   jhi = low+inc
        if (jhi > n) then
          jhi = n+1
        else
          if (liprm) then
            xn = xa(iprm(jhi))
          else
            xn = xa(jhi)
          endif
          if (x >= xn .eqv. ascnd) then
            low = jhi
            inc = inc+inc
            goto 1
          endif
        endif
      else
        jhi = low
    2   low = jhi-inc
        if (low < 1) then
          low = 0
        else
          if (liprm) then
            xn = xa(iprm(low))
          else
            xn = xa(low)
          endif
          if (x < xn .eqv. ascnd) then
            jhi = low
            inc = inc+inc
            goto 2
          endif
        endif
      endif
    3 if (jhi-low == 1) then
C   ... Find the first of values equal to x
    4   continue
        if (low > 1) then
          if (liprm) then
            xn = xa(iprm(low-1))
          else
            xn = xa(low-1)
          endif
          if (xn == x) then
            low = low-1
            goto 4
          endif
        endif
C     ... if xa(low) = x, decrement low
        if (low >= 1) then
          if (liprm) then
            xn = xa(iprm(low))
          else
            xn = xa(low)
          endif
          if (xn == x) low = low-1
        endif
        return
      endif
      jm = (jhi+low)/2
      if (liprm) then
        xn = xa(iprm(jm))
      else
        xn = xa(jm)
      endif
      if (x > xn .eqv. ascnd) then
        low = jm
      else
        jhi = jm
      endif
      goto 3
      end
