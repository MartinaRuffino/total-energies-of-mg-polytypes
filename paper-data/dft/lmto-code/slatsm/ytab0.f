      recursive subroutine ytab0(opt,xn,yn,n,val,npoly,lrat,tol,xp,i1)
C- Find next zero in a tabulated function
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 look for roots for x>xp
Ci         :1 look for roots for x<xp
Ci         :2 look for roots for closest to xp
Ci   xn    :ordered list of points at which function tabulated
Ci   yn    :vector of function values at points xn
Ci   n     :number of points
Ci   val   :Use yn-val as actual function value
Ci   npoly :maximum order of polynomial
Ci   lrat  :T, use rational function interpolation
Ci         :F, use polynomial interpolation
Cio Inputs/Outputs
Cio  i1    :Index to point in xn bracketing root.  On input:
Cio        :i1<=0 => start new search starting at max(-i1,1)
Cio        :         New search considers initial yn=val as root.
Cio        :i1>=1 => continue search: skip over initial points until
Cio        :         yn is not val.
Cio        :i1 is incremented until i1=np or (yn(i1)-val)*(yn(i1+1)-val)<0
Cio        :i1=-1 => failed to find root
Cl Local variables
Cl   np    :number of points used for interpolation
Cio Inputs/Outputs
Cio  xp    :Input:  Initial value of x for which to seek root yn-val=0
Cio        :Output: estimate for root
Cr Remarks
Cu Updates
Cu   18 Aug 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lrat
      integer opt,n,npoly,i1
      double precision xn(n),yn(n),val,tol
C ... Local parameters
      double precision xtol,xp,dymx,yp,dy,wk(12),xl,xr
      integer np,ilow,il,ir,j

C ... End of list reached, w/out bracket
      if (opt == 2) then
C       Special cases
        i1 = n-1
        if (xp == xn(n) .and. yn(n) == val) return
        i1 = 1
        if (xp == xn(1) .and. yn(1) == val) return
        i1 = -1
        if (xp > xn(n) .or. xp < xn(1)) return
        if (n == 1) return
C       Bracket xp
        il = 0
        call huntx(xn,n,xp,0,il) ! point on left closest to xp
        if (il == 0 .or. il == n) return
C       Find first root to the right
        ir = -il
        xr = 0                  ! xr is output ... only purpose is to prevent compiler from complaining
        call ytab0(0,xn,yn,n,val,npoly,lrat,tol,xr,ir)
C       Find first root to the left
        il = -il
        call ytab0(1,xn,yn,n,val,npoly,lrat,tol,xl,il)
        if (il < 0 .and. ir < 0) return ! No root
        if (il == ir) then ! occurs when yn(il) is root
          if (il > 1 .and. yn(max(il-1,1)) == val) then
            if (abs(xn(il-1)-xp) < abs(xl-xp)) then
              il=il-1
              xl = xn(il)
            endif
          endif
          if (ir < n .and. yn(min(ir+1,n)) == val) then
            if (abs(xn(ir+1)-xp) < abs(xr-xp)) then
              ir=ir+1
              xr = xn(ir)
            endif
          endif
        endif
   10   if (il < 0) then  ! Only root to left
          i1 = ir
          xp = xr
          return
        endif
        if (ir < 0) then ! Only root to right
          i1 = max(il-1,1)
          xp = xl
          return
        endif
        if (xr-xp <= xp-xl) then ! Choose closest root
          il = -1
        else
          ir = -1
        endif
        goto 10
      elseif (opt == 0) then
        call brack0(0,yn,n,val,i1)
        if (i1 == n .and. yn(n) == val) then
          xp = xn(i1)
          return
        elseif (i1 >= n .or. xn(min(i1,n)) == xn(n)) then
          xp = xn(n)
          i1 = -1
          return
        endif
      elseif (opt == 1) then
        call brack0(1,yn,n,val,i1)
        if (i1 == 1 .and. yn(1) == val) then
          xp = xn(i1)
          return
        elseif (i1 <= 1 .or. xn(max(i1,1)) == xn(1)) then
          xp = xn(1)
          i1 = -1
          return
        endif
        i1 = i1-1
      else
        call rx('ytab0: improper value of argument opt')
      endif

C ... Setup for finding the root: get points nearest to midpoint
      np = min(npoly+1,n)
      ilow = min(max(1,i1-np/2),n-np+1)
      xp = xn(i1)/2+xn(i1+1)/2
      ilow = ilow-1
      do while (.true.)
        ilow = ilow+1
        if (ilow+np >= n) exit
        if (xp-xn(ilow) < xn(ilow+np)-xp) exit
      enddo
C      This does the same thing but above is more correct
C      np = min(npoly+1,n)
C      j = i1
C      call hunty(xn,n,np,xn(i1)/2+xn(i1+1)/2,j,ilow)

C ... Finding root with poly or rational function interpolation
      xtol = tol*(xn(n)-xn(i1))
      dymx = tol*max(maxval(yn),-minval(yn))
      ir = 0
      xp = xn(i1)
      do  j = 1, 100
C       call pshpr(1)
        if (lrat) then
          call ratint(xn(ilow),yn(ilow),np,xp,yp,dy)
        else
          call polinx(xn(ilow),yn(ilow),np,xp,dymx,yp,dy)
        endif
        yp = yp-val
        call rfalsi(xp,yp,xtol,0d0,xtol/10,xn(i1+1)-xn(i1),10,wk,ir)
C       call poppr
        if (ir == 0 .or. ir == 1 .or. ir == 3) exit
C       Root search failed
        if (j == 50) then
          i1 = -1
          return
        endif
      enddo
      if (ir /= 0 .and. ir /= 1 .and. ir /= 3) i1 = -1
      if (opt == 1 .and. i1 /= -1) i1=i1+1

      end

      integer function nbrack(yn,np,val)
C- Number of times a tabulated function brackets a value
C ----------------------------------------------------------------------
Ci Inputs
Ci   yn    :vector of points
Ci   np    :number of points
Ci   val   :Use yn-val as actual sequence
Co Outputs
Ci   nbrack:number of crossings
Cr Remarks
Cu Updates
Cu   18 Aug 12 First created
C ----------------------------------------------------------------------
      implicit none
      integer np
C ... Passed parameters
      double precision yn(np),val
C ... Local parameters
      integer i1

      nbrack = 0
      i1 = 0
      do  while (i1 <= np)
C       j = -i1; call brack0(0,yn,np,val,j); i1 = j
        call brack0(0,yn,np,val,i1)
        if (i1 < np .or. yn(i1) == val) nbrack = nbrack+1
        i1 = i1+1
        do  while (i1 <= np)
          if (yn(i1) /= val) exit
          i1 = i1+1
        enddo
      enddo
      end

      subroutine brack0(opt,yn,np,val,i1)
C- Bracket the next point where a tabulated sequence changes sign
C ----------------------------------------------------------------------
Ci Inputs
Co   opt   :0 search forward
Co         :1 search backward
Ci   yn    :vector of points
Ci   np    :number of points
Ci   val   :Use yn-val as actual sequence
Cio Inputs/Outputs
Cio  i1    :On input:
Cio        :|i1|>=np => nothing done
Cio        :i1<0  => start new search starting at max(-i1,1)
Cio        :         New search considers initial yn=val as root.
Cio        :i1=0  => Same as i1<0, but if opt=1, set initial i1=np
Cio        :i1>=1 => continue search: skip over initial points until
Cio        :         yn is not val.
Cio        :If opt=0,
Cio        :i1 is incremented until i1=np or (yn(i1)-val)*(yn(i1+1)-val)<0
Cio        :Otherwise
Cio        :i1 is decremented until i1=1 or (yn(i1)-val)*(yn(i1+1)-val)<0
Cr Remarks
Cr   To make search with i1>0 but "new search each time", do
Cr     j = -i1; call brack0(0,yn,np,val,j); i1 = j
Cu Updates
Cu   18 Aug 12 First created
C ----------------------------------------------------------------------
      implicit none
      integer opt,np,i1
C ... Passed parameters
      double precision yn(np),val
C ... Local parameters
      double precision tol,res,resold
      parameter (tol=0d0)
      integer k,i

      if (opt == 0 .and. iabs(i1) >= np) return
      if (opt /= 0 .and. iabs(i1) > np) return
C     Initial point
      if (i1 < 1) then
        if (opt == 0) then
          i1 = max(-i1,1)
        else
          i1 = min(-i1,np)
          if (i1 == 0) i1 = np
        endif
C     Skip initial points exactly matching val
      else
        if (opt == 0) then
          do  while (i1 < np)
            if (yn(i1) /= val) exit
            i1 = i1+1
          enddo
        else
          do  while (i1 > 1)
            if (yn(i1) /= val) exit
            i1 = i1-1
          enddo
        endif
      endif
      if (opt == 0 .and. i1 >= np) return
      if (opt == 1 .and. i1 <= 1) return
      res = yn(i1)
      if (opt == 0) then
        k = i1+1
        do  i = k, np
          resold = res
          res = yn(i)
          i1 = i-1
          if ((res-val)*(resold-val) <= 0) then
C         if (resold /= val .and. res == val) i1=i
            return
          endif
        enddo
        i1 = np
      else
        k = i1-1
        do  i = k, 1, -1
          resold = res
          res = yn(i)
          i1 = i+1
          if ((res-val)*(resold-val) <= 0) then
            return
          endif
        enddo
        i1 = 1
      endif

      end

C     Test brack0, ytab0
C      subroutine fmain
C      implicit none
C      integer np,i0,i1,nbrack,j
C      parameter (np=10)
C      double precision xn(np),yn(np),val,tol,xpr,xpp,xp
C      parameter (tol=1d-12)
C      data xn / 1d0,10d0,20d0,22d0,25d0,30d0,31d0,33d0,40d0,41d0/
CC     data yn / 1d0,2d0,5d0,3d0,-1d0,-2d0,2d0,4d0,1d0,2d0/
CC     data yn / 1d0,2d0,5d0,3d0,-1d0,-2d0,2d0,1d0,1d0,2d0/
C      data yn / 1d0,2d0,5d0,3d0,-1d0,-2d0,2d0,1d0,2d0,1d0/
CC     data yn / 1d0,2d0,5d0,3d0,-1d0,-2d0,2d0,1d0,1d0,1d0/
CC     data yn / 2d0,2d0,5d0,3d0,2d0,2d0,2d0,3d0,3d0,3d0/
C      val = 1
C
C      call info2(0,0,0,' xn = %10:2d',xn,0)
C      call info2(0,0,0,' yn = %10:2d',yn,0)
C
C      i1 = nbrack(yn,np,val)
C      call info2(0,0,0,' nbrack found %i zeros.',i1,0)
C
CC     goto 999
C
C      call info2(0,0,0,' bracket roots with brack0, forward ...',i1,0)
C      i1 = 0
C      do  while (i1 <= np)
CC       j = -i1; call brack0(0,yn,np,val,j); i1 = j
C        call brack0(0,yn,np,val,i1)
C        if (iabs(i1) < np) print "(i4,2f10.5)",
C     .    i1, yn(i1)-val,yn(i1+1)-val
C        if (i1 == np .and. yn(i1) == val) print "(i4,2f10.5)",
C     .    i1, yn(i1)-val
C        i1 = i1+1
C      enddo
C
C      i1 = 0
CC     print *, '!!';i1 = 7
C      do while (i1 < np)
C        i0 = i1
C        call pshpr(1)
C        call ytab0(0,xn,yn,np,val,6,.false.,tol,xpp,i0)
C        if (i0 < 0) exit
C        call ytab0(0,xn,yn,np,val,6,.true.,tol,xpr,i1)
C        if (i1 == -1) then
C          print *, 'ratint failed to find root'
C          i1 = i0
C        endif
C        call poppr
C        call info5(0,0,0,' i1=%i  yl=%d  yr=%d  '//
C     .    'xp(rat) = %d  xp(pol) = %d',
C     .    i1,yn(i1),yn(min(i1+1,np)),xpr,xpp)
C        i1 = i1+1
C      enddo
C
C      call info2(0,0,0,'  bracket roots with brack0, backward ...',i1,0)
C      i1 = -np
C      do  while (iabs(i1) >= 1)
C       call brack0(1,yn,np,val,i1)
C        if (iabs(i1) > 1) print "(i4,2f10.5)",
C     .   i1, yn(i1)-val,yn(i1-1)-val
C        if (i1 == 1 .and. yn(i1) == val) print "(i4,2f10.5)",
C     .    i1, yn(i1)-val
C        i1 = i1-1
C      enddo
C
C      i1 = -np
C      do while (iabs(i1) > 1)
C        i0 = i1
C        call pshpr(1)
C        call ytab0(1,xn,yn,np,val,6,.false.,tol,xpp,i0)
C        if (i0 < 0) exit
C        call ytab0(1,xn,yn,np,val,6,.true.,tol,xpr,i1)
C        if (i1 == -1) then
C          print *, 'ratint failed to find root'
C          i1 = i0
C        endif
C        call poppr
C        call info5(0,0,0,' i1=%i  yl=%d  yr=%d  '//
C     .    'xp(rat) = %d  xp(pol) = %d',
C     .    i1,yn(i1),yn(min(i1+1,np)),xpr,xpp)
C        i1 = i1-1
C      enddo
C
C  999 continue
C      xp = xn(1)-1
CC     xp = 39
C      do while (xp < xn(np))
C        xp = xp+1
C        xpp = xp
C        call pshpr(1)
C        call ytab0(2,xn,yn,np,val,6,.false.,tol,xpp,i0)
C        call poppr
C        print "(i4,2f10.5)", i0,xp,xpp
C      enddo
C      end
