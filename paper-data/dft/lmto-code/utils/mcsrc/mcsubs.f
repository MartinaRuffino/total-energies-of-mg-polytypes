C      subroutine cust(ns,osprs,os,nr,nc,icast)
CC- Customized function
C      implicit none
C      integer ns,osprs(ns),os(ns),nr(ns),nc(ns),icast(ns)
CC heap:
C      integer w(1)
C      common /w/ w
C
CC ... This makes t.k into exp(t.k)
C      call redfrr(os(ns),nr(ns)*nc(ns)*2)
C      call defrr(os(ns+1),1)
C      icast(ns) = 2
C      call cust2(w(os(ns)),nr(ns),nc(ns))
C      print '(''# -cust : make t.k into exp(2.pi.i t.k)'')'
C      end
C      subroutine cust2(s,nr,nc)
C      implicit none
C      integer nr,nc
C      double precision s(nr,nc,2)
C      double precision pi2,xx
C      integer i,j
C
C      pi2 = 8*datan(1d0)
C      do  10  i = 1, nr
C      do  10  j = 1, nc
C        xx = pi2 * s(i,j,1)
C        s(i,j,1) = dcos(xx)
C        s(i,j,2) = dsin(xx)
C   10 continue
C      end
      subroutine zcptoy(s1,nr,nc,s2)
C- Copy complex*16 w1 to (real,imag) s2
      implicit none
      integer nr,nc
      double precision s1(2,nr,nc),s2(nr,nc,2)
      integer ic,ir

      do  10  ic = 1, nc
      do  10  ir = 1, nr
        s2(ir,ic,1) = s1(1,ir,ic)
        s2(ir,ic,2) = s1(2,ir,ic)
   10 continue
      end
      subroutine pwrs(nr,nc,icast,pwr,sin,sout)
C- Raise matrix to a power
      implicit none
      integer nr,nc,icast
      double precision sin(nr,nc,*),sout(nr,nc,*),pwr
      double complex zz
      integer ir

      if (nc /= 1) call rx('pwrs not ready for nc>1')
C ... Real case
      if (mod(icast,10) == 1) then
        do  10  ir = 1, nr
   10   sout(ir,1,1) = sin(ir,1,1)**pwr
      else
        do  20  ir = 1, nr
          zz = dcmplx(sin(ir,1,1),sin(ir,1,2))**pwr
          sout(ir,1,1) = dble(zz)
          sout(ir,1,2) = dimag(zz)
   20   continue
      endif

      end

      subroutine mpy3(n1a,n2a,n3a,n1b,n2b,n3b,a,b,c)
C- Makes C(i,j,k) = sum_m A(i,j,m) B(i,m,k)
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1a  :1st dimension of A and C
Ci   n2a  :2nd dimension of A and C
Ci   n3a  :3rd dimension of A and 2nd dimension of B
Ci   n1b  :1st dimension of B
Ci   n2b  :2nd dimension of B (see n3a)
Ci   n3b  :3rd dimension of B and C
Ci   A     :left matrix
Ci   B     :right matrix
Co Outputs
Co   C     :result matrix
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Feb 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1a,n2a,n3a,n1b,n2b,n3b
      double precision A(n1a,n2a,n3a),B(n1b,n2b,n3b),C(n1a,n2a,n3b)
C ... Local parameters
      integer i,j,k,m

      if (n3a /= n2b) call rx('mpy3: need n3a=n2b')

      call dpzero(C,n1a*n2a*n3b)

      do  i = 1, n1a
      do  j = 1, n2a
      do  k = 1, n3b
      do  m = 1, n3a
        C(i,j,k) = C(i,j,k) + A(i,j,m)*B(i,m,k)
      enddo
      enddo
      enddo
      enddo

      end

      subroutine zmpy3(n1a,n2a,n3a,n1b,n2b,n3b,a,b,c)
C- Makes C(i,j,k) = sum_m A(i,j,m) B(i,m,k) (complex case)
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1a  :1st dimension of A and C
Ci   n2a  :2nd dimension of A and C
Ci   n3a  :3rd dimension of A and 2nd dimension of B
Ci   n1b  :1st dimension of B
Ci   n2b  :2nd dimension of B (see n3a)
Ci   n3b  :3rd dimension of B and C
Ci   A     :left matrix
Ci   B     :right matrix
Co Outputs
Co   C     :result matrix
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Feb 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1a,n2a,n3a,n1b,n2b,n3b
      double complex A(n1a,n2a,n3a),B(n1b,n2b,n3b),C(n1a,n2a,n3b)
C ... Local parameters
      integer i,j,k,m

      if (n3a /= n2b) call rx('mpy3: need n3a=n2b')

      call dpzero(C,n1a*n2a*n3b*2)

      do  i = 1, n1a
      do  j = 1, n2a
      do  k = 1, n3b
      do  m = 1, n3a
        C(i,j,k) = C(i,j,k) + A(i,j,m)*B(i,m,k)
      enddo
      enddo
      enddo
      enddo

      end

      subroutine intgx(intrps,x0,lsw,np,uplim,s,nr,nc,icast)
C- Replace second col by integral wrt first, in a polynomial approx
C ----------------------------------------------------------------------
Ci Inputs
Ci   intrps:string with integration options
Ci   x0    :lower limit for integration
Ci   lsw   :switches for integration (see itrpsw)
Ci         :lsw(1)=2 -> integration on a uniform mesh of points
Ci         :lsw(2)=polynomical order
Ci   np    :number of upper limits
Ci   s     :array holding abscissa (col 1) and ordinates (cols 2..)
Ci   nr    :number of rows
Ci   nc    :number of rows
Ci   icast :cast
Cio Inputs/Outputs
Ci  uplim :uplim(1..np,1) = list of upper limits for integration
Co        :uplim(1..np,2..nc) = values of integration
Cl Local variables
Cl  xnow  :last upper integration limit (either x0 or a mesh point)
Cl  is =  :index to current point in s
Cl  iup   :index to current upper limit working on
Cl  ls0   :flags whether xnow is a mesh point (old integ. limit)
Cl  ls1   :flags whether x1 is a mesh point (new integration limit)
Cl  summ  :integral from x0 to last mesh point
Cr Remarks
Cr   abscissa and uplim are both assumed to be ordered
Cu Updates
Cu   13 Dec 02
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character*(*) intrps
      integer nr,nc,icast,np,lsw(2)
      double precision s(nr,nc),uplim(np,nc),x0
C ... Local parameters
      logical ls0,ls1
      integer n,nn,i,low,hi,ic,is,iup
C     integer jx,ip
      double precision dymx,sumi(2),fuzz,fuzz0,xx,sum(nc),summ(nc),x1,xnow
      parameter (nn=50,fuzz0=1d-10)
C     double precision dy,xy(nn,2),wp(nn)

C ... setup and defaults
      if (nc < 2) return
      n = min(nr,lsw(2)+1)
      fuzz = fuzz0*(s(nr,1)-s(1,1))

C ... So far, trapezoidal rule on a mesh.  Only uplim(1) is used
      if (lsw(1) == 2) then
C  ... Get first, last point
        call huntx(s,nr,x0,0,low)
        low = max(1,low)
        hi = min(low+1,nr)
        if (abs(x0-s(low,1)) > abs(x0-s(hi,1))) low = low+1
        call huntx(s,nr,uplim,0,hi)
        hi = min(nr,hi)
        if (hi < nr) then
          if (abs(uplim(1,1)-s(hi,1)) > abs(uplim(1,1)-s(hi+1,1)))
     .      hi = hi+1
        endif
        np = hi-low+1
        call pvintg(low,nr,np,s,uplim)
        return
      endif

C ... Initial uplim = 0 and set dymx
      call dpzero(summ,nc)
      do  i = 1, np
        do  ic = 2, nc
          uplim(i,ic) = 0
        enddo
      enddo
      dymx = 0
      do  i = 1, nr
        do  ic = 2, nc
          dymx = max(dymx,s(i,ic))
        enddo
      enddo
      dymx = fuzz0*dymx

C ... Set up first point
      is = 0
      call huntx(s,nr,x0,0,is)
      xnow = x0
      is = max(1,is)
      ls0 = .false.
      if (abs(x0-s(is,1)) < fuzz) then
        ls0 = .true.
      elseif (is < nr .and. abs(x0-s(is+1,1)) < fuzz) then
        is = is+1
        ls0 = .true.
      endif
      iup = 1

C --- Integrate in interval xnow, x1 ---
C     Here, x1 is the minimum of the next mesh point or uplim
C     The possibilities:
C     xnow and x1 are both mesh points
C     xnow and x1 are both mesh points
   10 continue

        if (is > nr) then
          is = is-1
          x1 = uplim(iup,1)
        else
          x1 = min(s(is,1),uplim(iup,1))
        endif
        ls1 = x1 == s(is,1)
C   ... Case integrate from one mesh point to next one
        if (ls0 .and. ls1) then
          if (abs(xnow-x1) < fuzz) then
            xnow = x1
            is = is+1
            goto 10
          endif
          do  ic = 2, nc
            call politg(s,s(1,ic),nr,n,is-1,is,0,xx,sumi)
            summ(ic) = summ(ic) + sumi(2)
          enddo
          xnow = x1
C         Mesh point also happens to be a limit point
          if (abs(x1-uplim(iup,1)) < fuzz) then
            do  ic = 2, nc
              uplim(iup,ic) = summ(ic)
            enddo
            iup = iup+1
            if (iup > np) goto 20
          endif

C   ... Case either xnow or x1 is not a mesh point
C       Note: xnow must be either a mesh point or x0
        else
          call dpzero(sum,nc)
          call pvint2(s,nr,nc,n,xnow,x1,dymx,sum)

C         x1 is a limit point
          if (abs(x1-uplim(iup,1)) < fuzz) then
            do  ic = 2, nc
              uplim(iup,ic) = sum(ic) + summ(ic)
            enddo
            iup = iup+1
          endif

C         Case x1 is mesh point (then xnow is 1st point)
C         Put xnow on a mesh; accumulate summ to mesh
          if (ls1) then
            do  ic = 2, nc
              summ(ic) = summ(ic) + sum(ic)
            enddo
C           is = is+1
            xnow = x1
            ls0 = .true.
          endif
        endif

        if (iup <= np) goto 10

C       End of integration
   20   continue

C      if (is < 1) then
C        is = 0
C        iup = 0
C   10   continue
C        iup = iup+1
C        x1 = min(s(1,1),uplim(iup,1))
C        call pvint2(s,nr,nc,n,x0,x1,dymx,np,iup,uplim)
C        if (iup < np) then
C          if (uplim(iup+1,1) < s(1,1)) goto 10
C        endif
C      endif
C      is = max(1,is)


C ... dumb, but for now ...
C      do  10  ip = 1, np
C        call gausq(nn,x0,uplim(ip,1),xy,wp,0,000)
C        jx = 0
C        do  12  ic = 2, nc
C        do  20  i = 1, nn
C   20   call polint(s,s(1,ic),nr,n,xy(i,1),dymx,0,jx,xy(i,2),dy)
CC       call prm(0,' ',1,6,'(2f12.6)',0,xy,nn,2)
C        sum = 0
C        do  30  i  = 1, nn
C   30   sum = sum + wp(i)*xy(i,2)
C        uplim(ip,ic) = sum
C   12 continue
C   10 continue

      end
      subroutine pvint2(s,nr,nc,n,x0,x1,dymx,sum)
C- Evaluate integral of s(2..nc) between x0 and x1
C ----------------------------------------------------------------------
Ci Inputs
Ci   s     :array holding abscissa (col 1) and ordinates (cols 2..)
Ci   n     :number of points to include in interpolating polynomial
Ci   nr    :number of rows
Ci   nc    :number of rows
Ci   x0    :lower limit for integration
Ci   x1    :upper limit for integration
Ci   dymx  :tolerance in integration
Co Outputs
Co   uplim :result of integration added to uplim(ip,2..nc)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   13 Dec 02
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nc,n
      double precision s(nr,nc),sum(nc),x0,x1,dymx
C ... Local parameters
      integer nn,jx,i,ic
C     integer low,hi
      double precision dy,sumi
      parameter (nn=50)
      double precision xy(nn,2),wp(nn)

      call gausq(nn,x0,x1,xy,wp,0,000)
      jx = 0
      do  12  ic = 2, nc
        do  20  i = 1, nn
   20   call polint(s,s(1,ic),nr,n,xy(i,1),dymx,0,jx,xy(i,2),dy)
C       call prm(0,' ',1,6,'(2f12.6)',0,xy,nn,2)
        sumi = 0
        do  30  i  = 1, nn
   30   sumi = sumi + wp(i)*xy(i,2)
        sum(ic) = sumi
   12 continue

      end
      subroutine pvintg(low,nr,np,s,uplim)
C- Evaluate integral of s(2..nc) on mesh of points s(1) by trap. rule.
Co  uplim(1..np,1) = abscissa of integration
Co  uplim(1..np,2..nc) = values of integration
      implicit none
      integer low,nr,np,ip
      double precision s(nr,2),uplim(np,2),dx

      uplim(1,1) = s(low,1)
      uplim(1,2) = 0
      if (np < 2) return
      dx = s(2,1) - s(1,1)
      do  10  ip = 2, np
        uplim(ip,1) = s(ip+low-1,1)
        uplim(ip,2) = uplim(ip-1,2) + (s(ip+low-2,2)+s(ip+low-1,2))/2*dx
   10 continue
      end

      subroutine minv(lqinv,lreal,nr,s,iwk,wk,info)
C- inverse of a matrix
      implicit none
      logical lqinv,lreal
      integer nr,iwk(nr),info
      double precision s(nr,nr,2),wk(nr,2),det(4)
      integer ow,lwork
C heap:
      integer w(1)
      common /w/ w

      if (lreal) then
        if (lqinv) then
          call defrr(ow,nr*nr)
          call dqinv('n',s,nr,2,nr,w(ow),nr,info)
          call rlse(ow)
        else
C          call dgefa(s,nr,nr,iwk,info)
C          if (info /= 0) return
C          call dgedi(s,nr,nr,iwk,det,wk,1)
          call dgetrf(nr,nr,s,nr,iwk,info)
          if (info /= 0) return
          lwork = nr*64
          call defrr(ow,lwork)
          call dgetri(nr,s,nr,iwk,w(ow),lwork,info)
          call rlse(ow)
        endif
      else
        if (lqinv) then
          call defrr(ow,nr*(nr+1))
          call yyqinv('4',s,s(1,1,2),nr,2,nr,w(ow),nr,info)
          call rlse(ow)
        else
          call yygefa(s,s(1,1,2),nr,nr,iwk,info)
          if (info /= 0) call rx('minv: matrix singular')
          call yygedi(s,s(1,1,2),nr,nr,iwk,det,wk,wk(1,2),1)
        endif
      endif

      end
      subroutine s0mps1(s1,nr,nc,icast,sindx,nsi,lerr,smap)
C- Copy s1(s0(i)) into smap(i)
      implicit none
      integer nr,nc,icast,nsi
      logical lerr
      double precision s1(nr,nc,2),sindx(nsi),smap(nsi,nc,2)
      integer i,i1,j

      do  i = 1, nsi
        i1 = nint(sindx(i))
        call sanrg(lerr,i1,1,nr,'mc','index')
        if (i1 >= 1 .and. i1 <= nr) then
          do  j = 1, nc
            smap(i,j,1) = s1(i1,j,1)
          enddo
          if (mod(icast,10) == 2) then
            call rx('this branch of s0mps1 never checked')
            do  j = 1, nc
              smap(i,j,2) = s1(i1,j,2)
            enddo
          endif
        else
          do  j = 1, nc
            smap(i,j,1) = 0
          enddo
          if (mod(icast,10) == 2) then
            call rx('this branch of s0mps1 never checked')
            do  j = 1, nc
              smap(i,j,2) = 0
            enddo
          endif
        endif
      enddo
      end

      subroutine crss(s1,s2,n1,nr,nc,s12)
C- cross product of two vectors
C n1 = 1  => s12(i,1:3) = s1(1,1:3) x s2(i,1:3)
C n1 = nr => s12(i,1:3) = s1(i,1:3) x s2(i,1:3)
C Any other value of n1 is an error
      implicit none
      integer n1,nr,nc
      double precision s1(n1,nc),s2(nr,nc),s12(nr,nc)
      integer i

      call rxx(nc /= 3,'cross product needs vectors with nc=3')
      if (n1 /= 1) then
        call rxx(n1 /= nr,
     .    'cross product: first row dimension must be 1 or '//
     .    'match 2nd array')
      endif

      do  i = 1, nr
        if (n1 == 1) then
          s12(i,1) = s1(1,2)*s2(i,3) - s1(1,3)*s2(i,2)
          s12(i,2) = s1(1,3)*s2(i,1) - s1(1,1)*s2(i,3)
          s12(i,3) = s1(1,1)*s2(i,2) - s1(1,2)*s2(i,1)
        else
          s12(i,1) = s1(i,2)*s2(i,3) - s1(i,3)*s2(i,2)
          s12(i,2) = s1(i,3)*s2(i,1) - s1(i,1)*s2(i,3)
          s12(i,3) = s1(i,1)*s2(i,2) - s1(i,2)*s2(i,1)
        endif
      enddo
      end


      subroutine cmprd(s1,s2,nr,nc,scale,lreal,lsum,lcc,linv,s12)
C- divides or multiplies each element of b into each elemnt of a
      implicit none
      logical lreal,lsum,lcc,linv
      integer nr,nc,i,j
      double precision s1(nr,nc,2),s2(nr,nc,2),s12(nr,nc,2),scale,sr,si,
     .  c1,c2,fuzz,xxr
      double complex xx

      fuzz = 1d-10
      c1 = scale
      c2 = scale
      if (lcc) c2 = -scale

C --- Real matrix ---
      if (lreal) then
        if (.not. lsum) then
          if (linv) then
            do  38  i = 1, nr
            do  38  j = 1, nc
              if (dabs(s2(i,j,1)) < fuzz) then
                xxr = 0
              else
                xxr = s1(i,j,1)/s2(i,j,1)
              endif
              s12(i,j,1) = c1*xxr
   38       continue
          else
            do  40  i = 1, nr
            do  40  j = 1, nc
   40       s12(i,j,1) = c1*s1(i,j,1)*s2(i,j,1)
          endif
        else
          sr = 0
          do  42  i = 1, nr
          do  42  j = 1, nc
   42     sr = sr + c1*s1(i,j,1)*s2(i,j,1)
          s12(1,1,1) = sr
          nr = 1
          nc = 1
        endif
        return
      endif

C --- Complex matrix ---
      if (.not. lsum) then
        if (linv) then
          do  8  i = 1, nr
          do  8  j = 1, nc
            if (dabs(s2(i,j,1)) < fuzz) s2(i,j,1) = 0d0
            if (dabs(s2(i,j,2)) < fuzz) s2(i,j,2) = 0d0
            if (dcmplx(s2(i,j,1),s2(i,j,2)) == (0d0,0d0)) then
              xx = 0
            else
            xx = dcmplx(s1(i,j,1),s1(i,j,2))/dcmplx(s2(i,j,1),s2(i,j,2))
            endif
            s12(i,j,1) = c1*dble(xx)
            s12(i,j,2) = c2*dimag(xx)
    8     continue
        else
          do  10  i = 1, nr
          do  10  j = 1, nc
            s12(i,j,1) = c1*s1(i,j,1)*s2(i,j,1) - c2*s1(i,j,2)*s2(i,j,2)
            s12(i,j,2) = c1*s1(i,j,1)*s2(i,j,2) + c2*s1(i,j,2)*s2(i,j,1)
   10     continue
        endif
      else
        sr = 0
        si = 0
        do  20  i = 1, nr
        do  20  j = 1, nc
          sr = sr + c1*s1(i,j,1)*s2(i,j,1) - c2*s1(i,j,2)*s2(i,j,2)
          si = si + c1*s1(i,j,1)*s2(i,j,2) + c2*s1(i,j,2)*s2(i,j,1)
   20   continue
        s12(1,1,1) = sr
        s12(2,1,1) = si
        nr = 1
        nc = 1
      endif

      end
      subroutine dmprd(n,a,b)
C- multiplies each element of b into each element of a
      implicit none
      integer n,i
      double precision a(n),b(n)
      do  10  i = 1, n
   10 a(i) = a(i)*b(i)
      end
      subroutine zmprd(n,a,b)
C- multiplies each element of b into each element of a (complex)
      implicit none
      integer n,i
      double precision a(n,2),b(n,2),xr,xi
      do  10  i = 1, n
        xr = a(i,1)*b(i,1) - a(i,2)*b(i,2)
        xi = a(i,1)*b(i,2) + a(i,2)*b(i,1)
        a(i,1) = xr
        a(i,2) = xi
   10 continue
      end
      subroutine sumr(icast,src,ir1,ir2,nr,nc,res)
C- Sums rows of src into res
      implicit none
      integer nr,nc,i,icast,ir1,ir2
      double precision src(nr,nc,2),res(nc,2),dsum

      do  10  i = 1, nc
   10 res(i,1) = dsum(ir2-ir1+1,src(ir1,i,1),1)

      if (mod(icast,10) == 2) then
        do  20  i = 1, nc
   20   res(i,2) = dsum(ir2-ir1+1,src(ir1,i,2),1)
      endif
      end
      subroutine sumrx(icast,src,nrs,ncs,dest,nr,nc,num,mxnum)
C- Sums rows of src into dest.  Rows given by num (see sum columns or rows, mc.f)
      implicit none
      integer icast,nrs,ncs,nr,nc,mxnum
      double precision src(nrs,ncs,2),dest(nr,nc,2)
      integer num(mxnum,nr)
      integer ir,irs,k

      call dpzero(dest,nr*nc*mod(icast,10))
      do  ir = 1, nr
        do  irs = 1, mxnum
          k = num(irs,ir)
          if (k == 0 .or. k > nrs) cycle
          call daxpy(nc,1d0,src(k,1,1),nrs,dest(ir,1,1),nr)
          if (mod(icast,10) == 1) cycle
          call daxpy(nc,1d0,src(k,1,2),nrs,dest(ir,1,2),nr)
        enddo
      enddo

      end
      subroutine nams(ns,os,nr,nc,icast,os0,nr0,nc0,ic0)
C- Pushes top element below stack and assigns o0 to that element
      implicit none
      integer ns,os(20),nr(20),nc(20),icast(20),os0,nr0,nc0,ic0
      integer is

      if (ns == 0) return

C      call wkprnt(1)
C      call wkinfo
C      call wkprnt(0)
C      print *, '----'

C ... Push top matrix to bottom of stack (poor man's way) ...
      do  10  is = 0, ns-1
   10 call stktog(is,ns,os,nr,nc,icast)

C ... Assign bottom element of stack to '0' ...
      nr0 = nr(1)
      nc0 = nc(1)
      ic0 = icast(1)
      os0 = os(1)

C ... Purge stack of bottom element ...
      do  20  is = 1, ns-1
        nr(is) = nr(is+1)
        nc(is) = nc(is+1)
        icast(is) = icast(is+1)
        os(is) = os(is+1)
   20 continue
      call clears(ns,nr,nc,icast)
C      call wkprnt(1)
C      call wkinfo
C      call wkprnt(0)

      end

      subroutine sscast(ns,os,nr,nc,icast)
C- Force top two elements of stack to take the larger cast
      implicit none
      integer ns,os(20),nr(20),nc(20),icast(20)
      integer sz,i
      logical ltog
      sz(i) = mod(icast(i),10)

      if (ns < 2) return
      if (sz(ns) == sz(ns-1)) then
        return
      elseif (sz(ns) < sz(ns-1)) then
        ltog = .false.
      else
        call stktog(0,ns,os,nr,nc,icast)
        ltog = .true.
      endif
      call nwcast(ns,os,nr,nc,icast,max(sz(ns),sz(ns-1)))
      if (ltog) call stktog(0,ns,os,nr,nc,icast)
      end
      subroutine nwcast(ns,os,nr,nc,icast,nwc)
C- Changes cast of top- or second-level matrix
C (not checked for complex->real)
      implicit none
      integer ns,os(20),nr(20),nc(20),icast(20),nwc
      integer sz,i
C heap:
      integer w(1)
      common /w/ w

      sz(i) = mod(icast(i),10)

C     call wkprnt(1)
      i = nr(ns)*nc(ns)*mod(nwc,10)
      call rlse(os(ns))
      call defrr(os(ns),i)
C ... real -> complex
      if (sz(ns) < mod(nwc,10)) then
C       call prm(0,' ',2,6,'(5f12.6)',0,w(os(ns)),nr(ns),nc(ns))
        call defrr(os(ns+1),i/2)
        call dcopy(i/2,w(os(ns)),1,w(os(ns+1)),1)
        call dpzero(w(os(ns)),i)
        call dcopy(i/2,w(os(ns+1)),1,w(os(ns)),1)
C       call prm(0,' ',2,6,'(5f12.6)',0,w(os(ns)),nr(ns),nc(ns))
        call rlse(os(ns))
        call defrr(os(ns),i)
        icast(ns) = icast(ns) - sz(ns) + 2
C ... complex -> real
      elseif (sz(ns) > mod(nwc,10)) then
        icast(ns) = icast(ns) - sz(ns) + 1
      endif
C     call wkprnt(0)
      end
      integer function getev(n,h,o,icast,lov,sw,e,z,iwk)
C sw: 0, return evals, 1, return evals e and evecs z
      implicit none
      integer n,sw,icast,iwk(n)
      logical lov
      double precision h(n,n*2),o(n*n*2),e(n*2),z(n*n*2)
      double precision abnrm
C Local variables
C     integer ofv1,ofv2,ofv3
      integer i,nev,lovi,owk,ierr
C heap:
      integer w(1)
      common /w/ w

      getev = -1
      if (icast == 1)  then
        call rx('getev not ready for real')
        if (lov) return
      elseif (icast == 11) then
        call defdr(owk,n*11)
        i = n
        if (sw == 0) i = 0
        call dsev1(n,h,o,w(owk),0,.true.,lov,0,n,9d69,nev,z,e)
C   ... this one makes evecs by inverse iteration
C       call dsev1(n,h,o,w(owk),0,.true.,lov,1,i,9d69,nev,z,e)
        getev = 1
        call rlse(owk)
      elseif (icast == 2) then
        if (lov) return

        call ztoyy(h,n,n,n,n,0,1)
        call zgeevs('P','N','V','N',n,h,n,e,e,1,z,n,abnrm,ierr)
        call ztoyy(e,n,1,n,1,1,0)
        call ztoyy(z,n,n,n,n,1,0)
C        call defdr(ofv1,n)
C        call defdr(ofv2,n)
C        call defdr(ofv3,n)
C        call cg(n,n,h,h(n**2+1,1),e,e(n+1),sw,z,z(n**2+1),
C     .    w(ofv1),w(ofv2),w(ofv3),ierr)
C        call rlse(ofv1)
        if (ierr /= 0) call rx('failed to diagonalize matrix')
        getev = 2
        nev = n
C       nev = 30
C   ... Sort eigenvalues/vectors by increasing |e|
        do  50  i = 1, nev
          h(2*i-1,1) = e(i)
          h(2*i,1) = e(n+i)
   50   continue
        call dvshel(2,nev,h,iwk,11)
        call dvperm(1,nev,e,h,iwk,.true.)
        call dvperm(1,nev,e(n+1),h,iwk,.true.)
        call dvperm(n,nev,z,h,iwk,.true.)
        call dvperm(n,nev,z(n**2+1),h,iwk,.true.)
C        call dvheap(2,nev,h,iwk,0d0,11)
C        call dvprm(1,nev,e,h,iwk,.true.)
C        call dvprm(1,nev,e(n+1),h,iwk,.true.)
C        call dvprm(n,nev,z,h,iwk,.true.)
C        call dvprm(n,nev,z(n**2+1),h,iwk,.true.)


      elseif (icast == 12) then
        call defdr(owk,n*11)
        i = n
        if (sw == 0) i = -1
C       call prm(0,' ',2,6,'(5f12.6)',0,h,n,n)
C       call prm(0,' ',2,6,'(5f12.6)',0,o,n,n)
        lovi = 0
        if (lov) lovi = 1
        call diagno(n,h,o,w(owk),.false.,lovi,.false.,i,9d69,nev,z,e)
C       call diagno(n,h,o,w(owk),.true.,lovi,.false.,i,9d69,nev,z,e)
        getev = 2
        if (i == -1) getev = 1
        call rlse(owk)
      endif
      end
      subroutine clrop(opstk,nops)
C- Pop the last operator of the the operator stack
      implicit none
      integer nops,opstk(nops),i
      do  10  i = 1, nops-1
   10 opstk(i) = opstk(i+1)
      nops = nops-1
      end
      subroutine mtrns(icast,s,st,nr,nc)
C- Transpose of a matrix
      implicit none
      integer nr,nc,icast
      double precision s(nr,nc,2), st(nc,nr,2)
      integer i,j
      do  10  i = 1, nr
      do  10  j = 1, nc
   10 st(j,i,1) = s(i,j,1)

      if (mod(icast,10) == 2) then
        do  20  i = 1, nr
        do  20  j = 1, nc
   20   st(j,i,2) = s(i,j,2)
      endif

      end
      subroutine mapdat(nexpr,expr,sincl,iwk,nr,nc,dat,dat2)
C- Replace columns of dat with algebraic expressions of them
      implicit none
C Passed Parameters
      integer nexpr,nr,nc,iwk(nr)
      character*(*) expr(1), sincl
      double precision dat(nr,1),dat2(*)
C Local Variables
      integer ir,iv0,ival,i,j,ii,jr
      character*4 xn, outs*80
      logical a2bin,logi

      call numsyv(iv0)
      jr = 0
      do  20  ir = 1, nr
        call clrsyv(iv0)
        call addsyv('i',dble(ir),ival)
C ---   Load data table ---
        do  22  j = 1, nc
          ii = 1
          xn = 'x   '
          call bin2a(' ',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ir,j),ival)
   22   continue
C ---   Exclude points if not satisfy sincl ---
        logi = .true.
        if (sincl /= ' ') then
          j = 0
          if (.not. a2bin(sincl,logi,0,0,' ',j,-1))
     .      call rx('mapdat:  error parsing sincl')
        endif
C ---   Put expressions of these vars into dat ---
        do  24  i = 1, nexpr
          j = 0
          if (.not. a2bin(expr(i),dat(ir,i),4,0,' ',j,-1)) then
            outs = expr(i)
            call skpblb(expr(i),len(expr(i)),ii)
            call fexit(-1,1,'MAPDAT: failed to parse expr: '
     .        //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
          endif
   24   continue
        if (logi) then
          jr = jr+1
          iwk(jr) = ir
        endif
   20 continue
      if (nexpr /= 0) nc = nexpr
      if (jr /= nr) then
        call dcopy(nr*nc,dat,1,dat2,1)
        call rowpmt(nr,jr,nc,iwk,dat2,dat)
        nr = jr
      endif
      call clrsyv(iv0)

      end
      subroutine rowpmt(nf,nt,nc,ipmt,afrom,ato)
      implicit none
      integer nf,nt,nc,ipmt(nt),i
      double precision afrom(nf,nc), ato(nt,nc)

      do  10  i = 1, nt
   10 call dcopy(nc,afrom(ipmt(i),1),nf,ato(i,1),nt)
      end
      subroutine colswp(nr,nc,s2,icm,ncout,s)
      implicit none
      integer nr,nc,icm(1),ncout
      double precision s(nr,nc),s2(nr,nc)
      integer ic,ir,ix
      do  20  ic = 1, ncout
        ix = icm(ic)
        if (ix > nc) print *, 'MCAT:  col',ix,' exceeds max col=',nc
        if (ix > nc) stop 'MCAT:  bad column mapping'
        do  30  ir = 1, nr
   30   s(ir,ic) = s2(ir,ix)
   20 continue
      end
      subroutine msort(nr,nc,icast,iwk,wk,s,ssort)
C- sort matrix according to increasing wk.
      implicit none
      integer nr,nc,icast,iwk(nr)
      double precision wk(nr),s(nr,nc),ssort(nr,nc)
      integer ir,ic,ix

      if (mod(icast,10) /= 1) call rx('msort ready only for real')
      call dvshel(1,nr,wk,iwk,1)
      do  20  ir = 1, nr
        ix = iwk(ir)+1
        do  30  ic = 1, nc
   30   ssort(ir,ic) = s(ix,ic)
   20 continue
      end
      logical function rowswp(nr,nc,s2,irm,irmadd,irmscl,nrout,s)
C - row swapping
      implicit none
      integer nr,nc,irm(1),nrout,irmadd,irmscl
      double precision s(nr,nc),s2(nr,nc)
      integer ir,ic,ix
      rowswp = .false.
      do  20  ir = 1, nrout
        ix = irm(ir)*irmscl+irmadd
        if (ix > nr .or. ix < 1) return
        do  30  ic = 1, nc
   30   s(ix,ic) = s2(ir,ic)
   20 continue
      rowswp = .true.
      end

      subroutine rccatm(opnum,ns,icast,nr,nc,os)
C- Concatenate top two matrices on stack
C ----------------------------------------------------------------------
Ci Inputs
Ci   opnum :1s digit
Ci         :0 do nothing
Ci         :1 concatenate columns
Ci         :2 concatenate rows
Ci         :10s digit
Ci         :0 concatenate in normal order
Ci         :1 concatenate in reverse order
Ci   ns
Ci   icast
Ci   nr    :
Ci   nc    :
Ci   w(os) :Screened structure constant matrix (asastr.f)
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 Jun 05
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opnum,ns,os(ns+1),nr(ns+1),nc(ns+1),icast(ns+1)
C ... Local parameters
      integer sz
      integer i,i1mach
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      procedure(integer) :: getdig
      external awrit0,catm,defrr,fexit2,sscast
      sz(i) = mod(icast(i),10)

      if (mod(opnum,10) == 0) return

      if (sz(ns-1) /= sz(ns)) then
        call awrit0('#ccat (warning): cast mismatch ...',' ',80,i1mach(2))
        call sscast(ns,os,nr,nc,icast)
      endif
      if (mod(opnum,10) == 1) then
        nr(ns+1) = nr(ns-1)
        nc(ns+1) = nc(ns-1)+nc(ns)
        if (nr(ns-1) /= nr(ns))
     .    call fexit2(-1,1,'Exit -1 mc: ccat expected nr=%i '//
     .    'but found nr=%i',nr(ns-1),nr(ns))
      elseif (mod(opnum,10) == 2) then
        nr(ns+1) = nr(ns-1)+nr(ns)
        nc(ns+1) = nc(ns-1)
        if (nc(ns-1) /= nc(ns))
     .    call fexit2(-1,1,'Exit -1 mc: rcat expected nc=%i '//
     .    'but found nc=%i',nc(ns-1),nc(ns))
      endif
      i = min(sz(ns-1),sz(ns))
      icast(ns+1) = i
      if (icast(ns-1) == icast(ns)) icast(ns+1) = icast(ns)
      icast(ns+1) = icast(ns+1) - 10*mod(icast(ns+1)/10,10)
      call defrr(os(ns+1),nr(ns+1)*nc(ns+1)*i)

      if (mod(opnum/10,10) == 0) then
        call catm(icast(ns+1),
     .    w(os(ns-1)),nr(ns-1),nc(ns-1),
     .    w(os(ns)),nr(ns),nc(ns),
     .    w(os(ns+1)),nr(ns+1),nc(ns+1))
      else
        call catm(icast(ns+1),
     .    w(os(ns)),nr(ns),nc(ns),
     .    w(os(ns-1)),nr(ns-1),nc(ns-1),
     .    w(os(ns+1)),nr(ns+1),nc(ns+1))
      endif
      end

      subroutine catm(icast,s1,n1r,n1c,s2,n2r,n2c,s12,n12r,n12c)
C- Kernel to concatenate two matrices called by rccatm
      implicit none
      integer n1r,n1c,n2r,n2c,n12r,n12c,icast
      double precision s1(n1r,n1c),s2(n1r,n1c),s12(n12r,n12c)

      call dmcpy(s1,n1r,1,s12,n12r,1,n1r,n1c)
      call dmcpy(s2,n2r,1,s12(1+n12r-n2r,1+n12c-n2c),n12r,1,n2r,n2c)
      if (mod(icast,10) == 2) then
        call dmcpy(s1(1+n1r*n1c,1),n1r,1,s12(1+n12r*n12c,1),n12r,
     .    1,n1r,n1c)
        call dmcpy(s2(1+n2r*n2c,1),n2r,1,
     .    s12(1+n12r*n12c+n12r-n2r,1+n12c-n2c),n12r,1,n2r,n2c)
      endif

      end
      subroutine mergesn(n,ns,icast,nr,nc,os)
C- Collapse part of stack, shifting s(ns+1) to s(ns-n)
C ----------------------------------------------------------------------
Ci Inputs
Ci   opnum :0 do nothing
Ci         :1 concatenate columns
Ci         :2 concatenate rows
Ci   n     :Merge n with ns
Ci   ns    :
Ci   icast
Ci   nr    :number of radial mesh points
Ci   nc
Ci   w(os) :Screened structure constant matrix (asastr.f)
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   03 Jun 05
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,ns,os(ns+1),nr(ns+1),nc(ns+1),icast(ns+1)
C ... Local parameters
      integer sz
      integer i
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external awrit0,catm,defrr,fexit2,sscast
      sz(i) = mod(icast(i),10)

C ... Copies s(ns+1) to s(ns-n) and reduces stack size by n ...
      call rlse(os(ns-n))
      call defrr(os(ns-n),nr(ns+1)*nc(ns+1)*sz(ns+1))
      call dcopy(nr(ns+1)*nc(ns+1)*sz(ns+1),w(os(ns+1)),1,w(os(ns-n)),1)
      nr(ns-n) = nr(ns+1)
      nc(ns-n) = nc(ns+1)
      icast(ns-n) = icast(ns+1)
      do  i = 1, n
        call clears(ns,nr,nc,icast)
      enddo

      end
      subroutine v2dia(ic,n,svec,mode,sm)
C mode=0: copy svec into sm.  mode=1: copy diagonal sm into svec
      implicit none
      integer ic,n,mode
      double precision svec(n,2), sm(n,n,2)
      if (mode == 0) then
        call dcopy(n,svec,1,sm,1+n)
        if (mod(ic,10) == 2) call dcopy(n,svec(1,2),1,sm(1,1,2),1+n)
      else
        call dcopy(n,sm,1+n,svec,1)
        if (mod(ic,10) == 2) call dcopy(n,sm(1,1,2),1+n,svec(1,2),1)
      endif
      end

      subroutine prm(opt,filel,icast,ifi,fmt,osprs,s,nr,nc)
C- Writes complex matrix to file ifi, Real & Im separated.  File is not rewound before writing
C ----------------------------------------------------------------
Ci Inputs:
Ci   opt:   1s digit:
Ci           0 writes in ascii mode
Ci           1 writes in binary mode
Ci          10s digit:
Ci           1 omit header when writing
Ci  filel:   reference string
Ci  icast:   0 integer
Ci           1 double precision
Ci           2 double complex with imaginary following real
Ci           3 double complex
Ci           4 double complex with imaginary following real in columns
Ci      Add 10 to indicate symmetric or hermitian
Ci    ifi:   file logical unit
Ci    fmt:   fortran format to write ascii string, e.g. '(5f12.6)'
Ci           You can use awrite format, e.g. '(%5,6g)'
Ci      s:   matrix to be printed out
Ci    ofi:   offset to imaginary part of s (used only if icast=2)
Ci    ns:    leading dimension of s
Ci           If icast=2, the given s is dimensioned s(ns,nc) (real part)
Ci                       and the imaginary part is dimensioned the same
Ci                       but offset by ofi
Ci           If icast=3, the given s is dimensioned s(ns,2,nc)
Ci           If icast=4, the given s is dimensioned s(2,ns,nc)
Ci    nr:    Number of rows to write
Ci    nc:    Number of columns to write
Cr Remarks
Cr   Binary write first record is
Cr     nr  nc  cast  optional-string-length
Cr   If a optional-string-length > 0, second record contains string
Cr   Next record is entire array, written as nr*nc elements
Cu Updates
Cu   3 Sep 99  argument list changed: offset to Im s is passed.
Cu  17 Mar 99  ywrm can write integer array
Cu  19 Jul 99  implement icast=4 format for ascii write
C ----------------------------------------------------------------
      implicit none
      integer opt,icast,nr,nc,ifi,osprs(2)
      character*(*) filel,fmt
      double precision s(nr,nc,2)
      integer nnz
C heap:
      integer w(1)
      common /w/ w

      if (icast < 100) then
        call ywrm(opt,filel,icast,ifi,fmt,s,nr*nc,nr,nr,nc)
      else
        nnz = osprs(1)
        call cprms(opt,filel,icast,ifi,fmt,nnz,w(osprs(2)),s,nr,nc)
      endif
      end
      subroutine apprm(filel,icast,ifi,fmt,s,nr,nc)
C- writes complex matrix to file ifi in (amplitude,phase) format
      implicit none
      character*(*) filel
      integer icast,nr,nc
      double precision s(nr,nc,2)
      character*(*) fmt, outs*40
      integer i,j,ifi
      double precision xx(nc)

      outs = ' '
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex amp'
      if (icast == 12) outs = ' herm  amp'
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
      if (mod(icast,10) == 1) then
        do  10  i = 1, nr
   10   write(ifi,fmt) (s(i,j,1), j=1,nc)
        return
      endif

      do  20  i = 1, nr
   20 write(ifi,fmt) (dsqrt(s(i,j,1)**2+s(i,j,2)**2), j=1,nc)
      write(ifi,'(1x)')
      do  22  i = 1, nr
        do  23  j = 1, nc
          xx(j) = 0
          if (s(i,j,2) /= 0 .or. s(i,j,1) /= 0) then
            xx(j) = datan2(s(i,j,2),s(i,j,1))
          endif
   23   continue
   22 write(ifi,fmt) (xx(j), j=1,nc)
      end

      subroutine xprm(nprec,icast,ifi,osprs,s,nr,nc)
      implicit none
      integer nprec,nr,nc,ifi,icast,osprs(2)
      double precision s(nr,nc,2),cdabs2,fac
      integer i,j,nnz,jfi
      character outs*10
C heap:
      integer w(1)
      common /w/ w

      if (icast >= 100) then
        nnz = osprs(1)
        call xprmp(nprec,icast,ifi,nnz,w(osprs(2)),s,nr,nc)
        return
      endif
      jfi = iabs(ifi)
      fac = 1/10d0**nprec
      outs = ' '
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex'
      if (icast == 12) outs = ' herm'
      call awrit2('%% rows %i cols %i sparse'//outs,' ',80,jfi,nr,nc)
      if (mod(icast,10) == 1) then
        if (ifi > 0) then
        do  10  i = 1, nr
        do  10  j = 1, nc
   10   if (dabs(s(i,j,1)) > fac) call awrit3(' %i %i %1;15g',
     .      ' ',80,jfi,i,j,s(i,j,1))
        else
        do  11  j = 1, nc
        do  11  i = 1, nr
   11   if (dabs(s(i,j,1)) > fac) call awrit3(' %i %i %1;15g',
     .      ' ',80,jfi,i,j,s(i,j,1))
        endif
      else
        if (ifi > 0) then
        do  20  i = 1, nr
        do  20  j = 1, nc
   20   if (cdabs2(s(i,j,1),s(i,j,2)) > fac) call awrit4(
     .      ' %i %i   %1;15g %1;15g',' ',80,jfi,i,j,s(i,j,1),s(i,j,2))
        else
        do  21  i = 1, nr
        do  21  j = 1, nc
   21   if (cdabs2(s(i,j,1),s(i,j,2)) > fac) call awrit4(
     .      ' %i %i   %1;15g %1;15g',' ',80,jfi,i,j,s(i,j,1),s(i,j,2))
        endif
      endif

      end
      subroutine xprmp(nprec,icast,ifi,nnz,isprs,s,nr,nc)
      implicit none
      integer nprec,nr,nc,ifi,icast,nnz,isprs(nnz,2)
      double precision s(nnz),fac
      integer i,j,ir
      character outs*10
C heap:
      integer w(1)
      common /w/ w

      fac = 1/10d0**nprec
      outs = ' '
      if (icast < 100) then
        call xprm(nprec,icast,ifi,w,s,nr,nc)
        return
      endif
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex'
      if (icast == 12) outs = ' herm'
      call awrit2('%% rows %i cols %i sparse'//outs,' ',80,ifi,nr,nc)
      if (mod(icast,10) == 1) then
        do  10  j = 1, nc
        do  10  i = isprs(j,1), isprs(j+1,1)-1
          ir = isprs(i,2)
          if (dabs(s(i)) > fac) call awrit3(' %i %i %1;15g',
     .      ' ',80,ifi,ir,j,s(i))
   10 continue
      else
        call rx('not ready for complex')
      endif

      end

      subroutine mscop(opt,ns,os,nr,nc,icast,nblk,scale)
C- Copies or adds subblock scale*s(ns-1) to subblock of s(ns)
C  opt = 0 => copy opt = 1 => add
      implicit none
      integer opt,nblk(4),ns,os(*),nr(ns),nc(ns),icast(ns)
      double precision scale(2)
      integer i1,i2,j1,j2,lgunit
      double precision beta(2)
C heap:
      integer w(1)
      common /w/ w


      if (ns < 2) call rx('mscop: needs 2 matrices')
      i1  = max(nblk(1),1)
      i2 = max(min(nblk(2),nr(ns)),i1-1)
      j1  = max(nblk(3),1)
      j2 = max(min(nblk(4),nc(ns)),j1-1)
      i2 = min(i2,i1-1+nr(ns-1))
      j2 = min(j2,j1-1+nc(ns-1))

      if (i1 /= nblk(1) .or. i2 /= nblk(2) .or.
     .    j1 /= nblk(3) .or. j2 /= nblk(4)) call awrit4(
     .  '#mclip (warning): subblock out of range; use %i,%i,%i,%i',
     .  ' ',80,lgunit(1),i1,i2,j1,j2)
      if (i2 < i1 .or. j2 < j1) return
      call sscast(ns,os,nr,nc,icast)
      icast(ns) = mod(icast(ns),10)
      beta = 0
      if (opt == 1) beta = [1d0, 0d0]
      if (mod(icast(ns),10) == 2) then
        call ymsadd(i2-i1+1,j2-j1+1,nr(ns-1),nr(ns),0,0,i1-1,j1-1,
     .      scale,beta,w(os(ns-1)),nr(ns-1)*nc(ns-1),
     .      w(os(ns)),nr(ns)*nc(ns))
      else
        if (opt == 0) then
          call dmscop(w(os(ns)),nr(ns),w(os(ns-1)),nr(ns-1),
     .      1,i2-i1+1,1,j2-j1+1,i1,j1,scale)
        else
          call dmsadd(w(os(ns)),nr(ns),w(os(ns-1)),nr(ns-1),
     .      1,i2-i1+1,1,j2-j1+1,i1,j1,scale)
        endif

      endif
      end
      subroutine mclip(lsubc,val,s,icast,nr,nc,sclip,nclip)
C- Clips a matrix to range of nclip
C     lsubc
C       1   exclude subblock of s0
C       2   scales  subblock of s0
C       3   copies # to subblock of s0
C       4   extract subblock of s0
      implicit none
      integer lsubc,icast,nr,nc,nclip(4)
      double precision s(nr,nc,*), sclip(1),val(2),xx(2)
      integer i,j,ii,jj,nrnew,ncnew,i1,i2,j1,j2,nrc,lgunit
C      integer imax
C      imax = 0

      i1  = max(nclip(1),1)
      i2 = max(min(nclip(2),nr),i1-1)
      j1  = max(nclip(3),1)
      j2 = max(min(nclip(4),nc),j1-1)
      if (lsubc == 1) then
        if (nclip(1) == -1) return
        nrnew = nr-i2+i1-1
        ncnew = nc-j2+j1-1
        nrc = nrnew*ncnew
        ii = 0
        do  10  i = 1, nr
          if (i >= i1 .and. i <= i2) goto 10
          ii = ii+1
          jj = 0
          do  12  j = 1, nc
            if (j >= j1 .and. j <= j2) goto 12
            jj = jj+1
            sclip(ii+nrnew*(jj-1)) = s(i,j,1)
            if (icast == 2) sclip(ii+nrnew*(jj-1)+nrc) = s(i,j,2)
   12     continue
   10   continue
        nr = nrnew
        nc = ncnew
      elseif (lsubc == 2 .or. lsubc == 3) then
        if (i1 /= nclip(1) .or. i2 /= nclip(2) .or.
     .      j1 /= nclip(3) .or. j2 /= nclip(4))
     .    call awrit0('#mclip (warning): subblock out of range',' ',
     .    80,lgunit(1))
        xx(1) = val(1)
        xx(2) = val(2)
        do  30  i = i1, i2
          do  32  j = j1, j2
            if (lsubc == 2) then
              xx(1) = s(i,j,1)*val(1)
              if (icast == 2) then
                xx(1) = s(i,j,1)*val(1) - s(i,j,2)*val(2)
                xx(2) = s(i,j,2)*val(1) + s(i,j,1)*val(2)
              endif
            endif
            s(i,j,1) = xx(1)
            if (icast == 2) s(i,j,2) = xx(2)
   32     continue
   30   continue
      elseif (lsubc == 4) then
        nrnew = nclip(2) - nclip(1) + 1
        ncnew = nclip(4) - nclip(3) + 1
        nrc = nrnew*ncnew
        call dpzero(sclip,nrc*icast)
        ii = 0
        do  20  i = nclip(1), nclip(2)
          ii = ii+1
          jj = 0
          do  22  j = nclip(3), nclip(4)
          jj = jj+1
          if (i >= 1.and.i <= nr .and. j >= 1.and.j <= nc) then
C           imax = max(imax,ii+nrnew*(jj-1)+nrc)
            sclip(ii+nrnew*(jj-1)) = s(i,j,1)
            if (icast == 2) sclip(ii+nrnew*(jj-1)+nrc) = s(i,j,2)
          else
          endif
   22     continue
   20   continue
        nr = nrnew
        nc = ncnew
      else
        call rx('mclip: bad lsubc')
      endif
C     print *, 'imax=',imax
      end
      subroutine addnws(nrnew,ncnew,icnew,ns,os,nr,nc,icast)
C- Allocate new matrix at top of stack
      implicit none
      integer nrnew,ncnew,icnew,ns,os(20),nr(20),nc(20),icast(20)
C heap:
      integer w(1)
      common /w/ w
C     sz(i) = mod(icast(i),10)

      ns = ns+1
      nr(ns) = nrnew
      nc(ns) = ncnew
      icast(ns) = icnew
      call defrr(os(ns),nrnew*ncnew*mod(icnew,10))
      end

      subroutine clears(ns,nr,nc,icast)
C- Clear matrix at top of stack
      implicit none
      integer ns,nr(ns),nc(ns),icast(ns)

      ns = max(ns-1,0)
      nr(ns+1) = 0
      nc(ns+1) = 0
      icast(ns+1) = 0

      end

      subroutine stkpsh(is,ns,os,nr,nc,icast)
C- Push matrix (is) onto stack.
      implicit none
      integer is,ns,os(20),nr(20),nc(20),icast(20),sz
      integer i
C heap:
      integer w(1)
      common /w/ w
      sz(i) = mod(icast(i),10)

      call addnws(nr(is),nc(is),icast(is),ns,os,nr,nc,icast)
      call dcopy(nr(is)*nc(is)*sz(is),w(os(is)),1,w(os(ns)),1)
      end
      subroutine stktog(is,ns,os,nr,nc,icast)
C- Toggles two matrices at levels (is,is-1) on stack.
      implicit none
      integer is,ns,os(20),nr(20),nc(20),icast(20),sz,i,j,ks
C heap:
      integer w(1)
      common /w/ w
      sz(i) = mod(icast(i),10)

      ks = ns-is
      if (ks < 2) return
      i = nr(ks)*nc(ks)*sz(ks)
      call defrr(os(ns+1),i)
      call dcopy(i,w(os(ks)),1,w(os(ns+1)),1)
C ... Copy array ks-1 to work space
      j = nr(ks-1)*nc(ks-1)*sz(ks-1)
      call defrr(os(ns+2),j)
      call dcopy(j,w(os(ks-1)),1,w(os(ns+2)),1)
C ... Reassign work space for ks-1 and ks, and copy arrays back there
      call rlse(os(ks-1))
      call defrr(os(ks-1),i)
      call defrr(os(ks),j)
      call dcopy(i,w(os(ns+1)),1,w(os(ks-1)),1)
      call dcopy(j,w(os(ns+2)),1,w(os(ks)),1)
      i = nr(ks-1)
      j = nc(ks-1)
      nr(ks-1) = nr(ks)
      nc(ks-1) = nc(ks)
      nr(ks) = i
      nc(ks) = j
      i = icast(ks-1)
      icast(ks-1) = icast(ks)
      icast(ks) = i
C ... Resssign work space for ks+1...ns
      do  10  i = ks+1, ns
   10 call defrr(os(i),nr(i)*nc(i)*sz(i))

      end
      subroutine itrpsw(itrpst,lsw)
C- Parse switches for interpolating functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   itrpst:string containing switches. Can consist of the following,
Ci         :separated by commas:
Ci         : rat
Ci         : mesh
Ci         : ord=#
Co Outputs
Co   lsw   :lsw(1)= 0 (default)
Co         :        1 if 'rat' parsed
Co         :        2 if 'mesh' parsed
Co         :lsw(2)= 4 (default)
Co         :        # if ord=# parsed
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   13 Dec 02
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lsw(2)
      character*(*) itrpst
C ... Local parameters
      integer it(10),k,parg
      double precision res(10)

      k = 0
      lsw(1) = parg('rat',0,itrpst,k,len(itrpst),', ',2,0,it,res)
      k = 0
      if (parg('mesh',0,itrpst,k,len(itrpst),', ',2,0,it,res) == 1)
     .  lsw(1) = 2
      k = 0
C ... default to second-order
      lsw(2) = 4
      if (parg('ord=',2,itrpst,k,len(itrpst),', ',2,1,it,lsw(2)) < 0)
     .  call rxs('itrpsw: failed to parse',itrpst)
      end

      subroutine intrp(dat,nrd,ncd,icast,res,nres,lsw)
      implicit none
      integer nrd,ncd,icast,nres,lsw(2)
      double precision dat(nrd,ncd),res(nres,ncd)
      logical lrat
      integer i,ic,j1,jx,lgunit,ozp,ofp
      double precision dy,dymx
C heap:
      integer w(1)
      common /w/ w

      if (ncd < 2) return
      jx = 0
      dymx = 1d-10
      lrat = lsw(1) == 1
      if (lsw(2) <= 0) lsw(2) = min(10,nrd)

      if (mod(icast,10) == 1) then
      do  10  ic = 2, ncd
      do  10  i = 1, nres
        if (lrat) then
C         Problem: discontintities between points.
C         call hunty(dat,nrd,lsw(2),res(i,1),jx,j1)
C         Solution: pick midway between enclosing points,
C         so change of interpolation occurs at points.
          call hunty(dat,nrd,2,res(i,1),jx,j1)
          if (min(dat(j1,1),dat(j1+1,1)) < res(i,1) .and.
     .        max(dat(j1,1),dat(j1+1,1)) > res(i,1)) then
            call hunty(dat,nrd,lsw(2),dat(j1,1)/2+dat(j1+1,1)/2,jx,j1)
          else
            call hunty(dat,nrd,lsw(2),res(i,1),jx,j1)
          endif
          call ratint(dat(j1,1),dat(j1,ic),min(lsw(2),nrd-j1+1),
     .      res(i,1),res(i,ic),dy)
C          print *, i,j1,min(lsw(2),nrd-j1+1),res(i,1),res(i,2),dy
          print *, 'i,j1=',i,j1,dat(j1,1)
        else
          call polint(dat,dat(1,ic),nrd,lsw(2),res(i,1),dymx,0,jx,
     .      res(i,ic),dy)
C         print *,  jx,res(i,1),res(i,2),dy
        endif
   10 continue
      call awrit3('#intrp: interpolate to list of %i pts)'//
     .  ' order=%i rat=%l',' ',80,lgunit(1),nres,lsw(2),lrat)
      else
        call defdc(ozp, nrd)
        call defdc(ofp, nrd)
        if ((nres/2)*2 /= nres)
     .    call rx('intrp: read odd number of points for complex data')
C       separate real, imaginary parts of res
        call ztoy(res,nres/2*ncd,nres/2,1,0)
        nres = nres/2
        call pditrp(dat,nrd,ncd,icast,w(ozp),w(ofp),res,nres)
        call rlse(ozp)

      endif

      end

      subroutine pditrp(dat,nrd,ncd,icast,zp,fp,res,nres)
      implicit none
      integer nrd,ncd,icast,nres
      double precision dat(nrd,ncd),res(nres,ncd,2)
      integer i,ic,ocof
      double precision zp(2,nrd),fp(2,nrd),z(2),fpad(2)
C heap:
      integer w(1)
      common /w/ w

      call dcopy(nrd,dat,1,zp,2)
      call dcopy(nrd,dat(1,ncd+1),1,zp(2,1),2)
      if (mod(icast,10) /= 2) call rx('pditrp: bad cast')
      call defdc(ocof, nrd*(nrd+2))

      do  20  ic = 2, ncd
        call dcopy(nrd,dat(1,ic),1,fp,2)
        call dcopy(nrd,dat(1,ncd+ic),1,fp(2,1),2)
        call padcof(nrd,zp,fp,nrd,w(ocof))
c        call zprm('fp',2,fp,nrd,nrd,1)
      do  20  i = 1, nres
        z(1) = res(i,1,1)
        z(2) = res(i,1,2)
        call pade(1,z,nrd,zp,nrd,w(ocof),fpad)
        res(i,ic,1) = fpad(1)
        res(i,ic,2) = fpad(2)
   20 continue

      end

      logical function diffx(lsw,s,nr,nc,icast)
C- Replace second col by derivative wrt first, in a polynomical approx
      implicit none
      integer nr,nc,icast,lsw(2)
      double precision s(nr,nc)
      integer oy,ic,lerr
      logical lrat
      double precision tol
C heap:
      integer w(1)
      common /w/ w

C ... setup and defaults
      if (lsw(2) <= 0) lsw(2) = 4
      lsw(2) = min(nr,lsw(2))
      lrat = lsw(1) == 1
      tol = 0d0
      diffx = .false.
      if (nc < 2) return
      diffx = .true.

      call defdr(oy,nr)
      do  10  ic = 2, nc
        call dcopy(nr,s(1,ic),1,w(oy),1)
        call poldvm(s,w(oy),nr,lsw(2),lrat,tol,lerr,s(1,ic))
        diffx = diffx .and. lerr == 0
   10 continue
      call rlse(oy)

      end
      subroutine expand(nt,it,xxv,dat,nr,nc)
C- Expands a list of points into a table
      implicit none
      integer nt,it(nt),nr,nc,ir,i,ipr,ic,nc0,i1mach
      double precision dat(10),xxv(nt),dx,x0,x1

      call getpr(ipr)
      ir = 0
      i = 0

C --- Check if number of columns specified ---
      nc0 = 1
      if (it(1) == 3) then
        nc0 = iabs(nint(xxv(1)))
        i = i+1
      endif
C --- Expand xxv into table ----
   10 i = i+1
      ir = ir+1
      dat(ir) = xxv(i)
C ...   Generate points in range
      call rxx(it(i) == 3,'expand: ! allowed only as first element')
      if (it(i) == 2) then
        x0  = xxv(i)
        i = i+1
        x1  = xxv(i)
        dx = 1
        if (it(i) == 2) then
          i = i+1
          dx = xxv(i)
        endif
   20   continue
        if (x0*(1-1d-12)-x1 < 0) then
          ir = ir+1
          x0 = x0+dx
          dat(ir) = x0
          goto 20
        endif
        ir = ir-1
      endif
C --- Quit if last: rearrange points according to ir,nc ---
      if (it(i) == 4) then
        nr = ir
        if (ipr >= 40) call awrit1(' expand:  %i points generated',
     .    ' ',80,i1mach(2),nr)
        if (nc0 /= 1 .and. mod(nr,nc0) /= 0)
     .    call rx('expand:  nr not a multiple of nc')
        nr = nr/nc0
        nc = nc0
        if (it(1) /= 3 .or. nc0 == 0) then
          nc = 2
          call dpzero(dat(nr+1),nr)
          return
        elseif (xxv(1) > 0) then
          call dcopy(nr*nc,dat,1,dat(1+nr*nc),1)
          call xxpand(nr,nc,dat(1+nr*nc),dat)
        endif
        if (ipr >= 80) then
          print *, nr,nc
          do  12  ir = 1, nr
   12     print 333, (dat(ir+(ic-1)*nr), ic=1,nc)
  333     format(3f12.6)
        endif
        return
      endif
      goto 10
      end
      subroutine xxpand(nr,nc,d1,d2)
C- Swaps around rows and columns
      implicit none
      integer nr,nc,ir,ic
      double precision d1(nc,nr), d2(nr,nc)

      do  10  ir = 1, nr
      do  10  ic = 1, nc
   10 d2(ir,ic) = d1(ic,ir)
      end
      subroutine crowl(lrow,lpermc,nt,iprm,icast,dat,nr,nc,s,nsr,nsc)
C- Create a array from a list of rows and columns of an existing one
C ----------------------------------------------------------------------
Ci Inputs
Ci   lrow  :T copy rows in permuted order, else columns
Ci   lpermc:T source file dat contains a table of columns,
Ci          columns are permuted??
Ci   nt    :number of rows (columns) to (permute and) copy
Ci   iprm  :permutation table
Ci   icast :
Ci   dat   :source matrix
Ci   nr    :row dimension of dat
Ci   nc    :col dimension of dat
Ci   s     :destination matrix
Ci   nsr   :row dimension of s
Ci   nsc   :col dimension of s
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lrow,lpermc
      integer nt,icast,nr,nc,nsr,nsc,iprm(nt)
      double precision dat(nr,nc,*),s(nsr,nsc,*)
C ... Local parameters
      integer i,k,ic

      if (lrow) then
        do  i = 1, nt
          if (iprm(i) > nr) call rxi('crowl: no row entry ',iprm(i))
          call dcopy(nc,dat(iprm(i),1,1),nr,s(i,1,1),nsr)
          if (icast == 2)call dcopy(nc,dat(iprm(i),1,2),nr,s(i,1,2),nsr)
      enddo
      elseif (lpermc) then
        do  k = 1, nr
          do  i = 1, nt
            ic = nint(dat(k,iprm(i),1))
            s(k,i,1) = dat(k,ic,1)
            if (icast == 2) s(k,i,2) = dat(k,ic,2)
          enddo
        enddo
      else
        do  i = 1, nt
          if (iprm(i) > nc) call rxi('crowl: no col entry ',iprm(i))
          call dcopy(nr,dat(1,iprm(i),1),1,s(1,i,1),1)
          if (icast == 2) call dcopy(nr,dat(1,iprm(i),2),1,s(1,i,2),1)
      enddo
      endif

C     call prm(0,' ',icast,6,'(5f12.6)',0,s,nsr,nsc)

      end

      subroutine unx(parms,s,nr,nc,icast)
C- Uncross lines
      implicit none
      integer nr,nc,icast
      double precision parms(4),s(nr,nc,1)
      integer i1,i2,ir,ir1
      double precision dy,dynew

      i1 = nint(parms(1))
      i2 = nint(parms(2))
      ir1 = nint(parms(3))

      dy = (s(ir1,i1,1)-s(ir1,i2,1))**2
      do  10  ir = ir1, nr
        dynew = (s(ir,i1,1)-s(ir,i2,1))**2
        if (dynew <= dy) then
          dy = dynew
        else
          ir1 = ir
          goto 11
        endif
   10 continue

C ... Uncross from here
   11 continue
      call dswap(nr-ir1+1,s(ir1,i1,1),1,s(ir1,i2,1),1)
      end

      subroutine unx2(parms,s,nr,nc,icast)
C- Cross lines that are not crossed, but
      implicit none
      integer nr,nc,icast
      double precision parms(4),s(nr,nc,1)
      integer i1,i2,ir,ir1,irmin
      double precision dy,dynew,wk1(3,2),wk2(3,2),cof(3)
      double precision sqrcf1,sqrcf2

      i1 = nint(parms(1))
      i2 = nint(parms(2))
      ir1 = nint(parms(3))
      if (ir1 >= nr) return

C ... Initial point: if data drawing closer, set irmin=1
      irmin = 0
      dy = (s(ir1,i1,1)-s(ir1,i2,1))**2
      dynew = (s(ir1+1,i1,1)-s(ir1+1,i2,1))**2
      if (dynew < dy) irmin = ir1

C ... Find point of closest approach,
      dy = (s(ir1,i1,1)-s(ir1,i2,1))**2
      do  ir = ir1+1, nr-1
        dynew = (s(ir,i1,1)-s(ir,i2,1))**2
        if (dynew < dy) then
          irmin = ir
          dy = dynew
C         If the lines already cross, ignore this point
          if ((s(ir-1,i1,1)-s(ir-1,i2,1))*(s(ir,i1,1)-s(ir,i2,1)) < 0)
     .      then
            irmin = 0
          endif
          if ((s(ir+1,i1,1)-s(ir+1,i2,1))*(s(ir,i1,1)-s(ir,i2,1)) < 0)
     .      then
            irmin = 0
          endif

C       If no point found yet, reset dy until points begin to draw closer
        elseif (irmin == 0) then
          dy = dynew
        endif
      enddo

      if (irmin == 0) then
        call info0(0,0,0,'#mc: no uncrossed point found')
        return
      endif

C ... If point of closest approach, look for smoothest change
      call dcopy(3,s(irmin-1,1,1),1,wk1(1,1),1)
      call dcopy(3,s(irmin-1,1,1),1,wk2(1,1),1)
      call dcopy(3,s(irmin-1,i1,1),1,wk1(1,2),1)
      call dcopy(3,s(irmin-1,i2,1),1,wk2(1,2),1)
      call dswap(1,wk1(3,2),1,wk2(3,2),1)
C     Sum of squares curvatures, if start swap at irmin+1
      call polcof(0,0d0,wk1,wk1(1,2),3,cof)
      sqrcf1 = cof(3)**2
      call polcof(0,0d0,wk2,wk2(1,2),3,cof)
      sqrcf1 = sqrcf1 + cof(3)**2
C     Curvature if swap middle point
      call dswap(1,wk1(2,2),1,wk2(2,2),1)
      call polcof(0,0d0,wk1,wk1(1,2),3,cof)
      sqrcf2 = cof(3)**2
      call polcof(0,0d0,wk2,wk2(1,2),3,cof)
      sqrcf2 = sqrcf2 + cof(3)**2
C     If 1st ordering has lower curvature, swap after irmin
      if (sqrcf1 <= sqrcf2) irmin = irmin+1
      call info2(10,0,0,'#mc: swapping cols at point %i',irmin,0)
      call dswap(nr-irmin+1,s(irmin,i1,1),1,s(irmin,i2,1),1)
      end

      subroutine findx(mode,parms,expr,dat,nr,nc,s,nrs)
C- Find two rows that straddle expr=val; return x1 that interpolates
C  mode 0: return nrs
C  mode 1: return s and s(1:nrs)
      implicit none
      integer nr,nc,mode,nrs
      double precision parms(4),dat(nr,nc),val,s(nrs,2)
      character*(*) expr,xn*4
      integer ir,ir1,iv0,ival,j,ii,nbrak
      double precision res,resold
      logical a2bin
      character*(120) outs

      val = parms(1)
      ir1 = nint(parms(2))
      if (ir1 >= nr) return

      call numsyv(iv0)
      res = 0
      nbrak = 0
      do  ir = ir1, nr-1
        resold = res
C   ... Load data table
        call clrsyv(iv0)
        call addsyv('i',dble(ir),ival)
        do  j = 1, nc
          ii = 1
          xn = 'x   '
          call bin2a(' ',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ir,j),ival)
        enddo
        j = 0
C       call shosyv(0,0,0,6)
        if (.not. a2bin(expr,res,4,0,' ',j,-1)) then
          outs = expr
          call skpblb(expr,len(expr),ii)
          call fexit(-1,1,'MAPDAT: failed to parse expr: '
     .      //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
        endif
C       No first point yet
        if (ir == ir1) then
          if (res == val) then
            nbrak = 1
            if (mode == 1) then
              s(nbrak,1) = dat(ir1,1)
              s(nbrak,2) = val
            endif
          endif
          cycle
        endif
C       No action unless expression is bracketed
        if ((res-val)*(resold-val) > 0) cycle
        if ((res-val)*(resold-val) < 0 .or. res == val) then
          nbrak = nbrak+1
          if (mode == 1) then
          if (res == val) then
            s(nbrak,1) = dat(ir,1)
            s(nbrak,2) = val
          else
            s(nbrak,1) = (dat(ir,1)*(val-resold)-dat(ir-1,1)*(val-res))/
     .                   (res-resold)
            s(nbrak,2) = val
          endif
          endif
        endif
      enddo

      if (mode == 0) then
        nrs = nbrak
      endif
      call clrsyv(iv0)

      end

      subroutine smooth(parms,dat,nr,nc,s,nsr,nsc)
C- Create a array from a list of rows and columns of an existing one
Ci parms:  1: gaussian width, 2,3,4 : xmin, xmax, dx
C
      implicit none
      integer nr,nc,i,nsr,nsc,n1,n2,n
      double precision parms(4),dat(nr,nc,1),s(nsr,nsc,1),srpi,xn

      srpi = dsqrt(4*datan(1d0))

C ... List of points on x axis
      do  10  i = 1, nsr
        s(i,1,1) = parms(2) + (i-1)*parms(4)
        s(i,2,1) = 0
   10 continue

C ... For each point in dat, smear out with a gaussian
C     Normalization is width*sqrt(pi)
      do  20  i = 1, nr
C   ... Points at gaussian peak +/- 10*width
        xn = (dat(i,1,1)-parms(2))/parms(4)+1
        n1 = min(max(nint(xn-10*parms(1)/parms(4)),1),nsr)
        n2 = min(max(nint(xn+10*parms(1)/parms(4)),1),nsr)
        do  22  n = n1, n2
          xn = parms(2) + parms(4)*n - dat(i,1,1)
          s(n,2,1) = s(n,2,1) + exp(-(xn/parms(1))**2)/(srpi*parms(1))
C         print *, n,xn,exp(-(xn/parms(1))**2)/(srpi*parms(1))
   22   continue
   20 continue

C     call prm(0,' ',1,6,'(5f12.6)',0,s,nsr,nsc)
c     stop

      end
C      subroutine rwargs(sfargs,slbl)
CC- Parses switches and arguments for read/write
C      implicit none
C      character sfargs*256,slbl*(*)
C      character dc*2,ca*1
C      integer j1,j2
C
C      j1 = 1
C      dc = sfargs(j1:j1) // ' '
C      j1 = j1+1
CC ... Return here to resume parsing for arguments
C   10 continue
C      call nwordg(sfargs,0,dc,1,j1,j2)
C      j2 = j2+1
C
CC ... Parse special arguments
C      if (sfargs(j1:j2) == ' ')  return
C      if (sfargs(j1:j1+1) == 'l=') then
C        call rx('-r:l= not implemented')
C      elseif (sfargs(j1:j1) == 's') then
C        slbl = sfargs(j1:j2)
C      else
C        call rxs2('mc failed to parse arg to -r or -w "',
C     .    sfargs(j1:j2),'"')
C      endif
C
C      j1 = j2+1
C      goto 10
C
C      end
      subroutine abss(s,nr,nc,icast)
C- Replace each element in s with absolute value
      implicit none
      integer nr,nc,icast,ir
      double precision s(nr,nc,*)

      do  10  ir = 1, nr*nc
        if (mod(icast,10) == 2) then
          s(ir,1,1) = abs(dcmplx(s(ir,1,1),s(ir,1,2)))
        else
          s(ir,1,1) = abs(s(ir,1,1))
        endif
   10 continue
      icast = icast - mod(icast,10) + 1
      end
      subroutine tonint(s,nr,nc,icast,facint)
C- Replace each element in s with nearest value
      implicit none
      integer nr,nc,icast,ir
      double precision s(nr,nc,*),facint

      do  ir = 1, nr*nc
        if (mod(icast,10) == 2) then
          s(ir,1,1) = abs(dcmplx(s(ir,1,1),s(ir,1,2)))
        endif
        s(ir,1,1) = nint(facint*s(ir,1,1))/facint
      enddo
      icast = icast - mod(icast,10) + 1
      end
      subroutine maxs(mode,s,nr,nc,icast)
C- Replace s(*,1) with max value in each row, or s(1) with global max
C  1s  digit mode: 0 max value is computed for each row
C                  1 max value is global
C  10s digit mode: 0 s(*,1) is replaced with max value
C                  1 s(*,1) is replaced with index to max value
      implicit none
      integer mode,nr,nc,icast,ir,ic,mr,mc
      double precision s(nr,nc,*),smax,sij
      logical lglob

      lglob = mod(mode,10) == 1
      if (lglob .and. mode >= 10) call rx('incompatible switches')
      mr = 0
      smax = 0
      do  ir = 1, nr
        if (.not. lglob) mr = 0

        do  ic = 1, nc
          if (icast == 2) then
            sij = abs(dcmplx(s(ir,ic,1),s(ir,ic,2)))
          else
            sij = s(ir,ic,1)
          endif
          if (sij > smax .or. mr == 0) then
            mr = ir
            mc = ic
            smax = sij
          endif
        enddo

        if (lglob) cycle

        if (mode >= 10) then
          s(ir,1,1) = mc
          if (icast == 2) s(ir,2,1) = 0
        else
          s(ir,1,1) = s(mr,mc,1)
          if (icast == 2) s(ir,2,1) = s(mr,mc,2)
        endif
      enddo

      if (lglob) then
        s(1,1,1) = mr
        s(2,1,1) = mc
        s(3,1,1) = smax
      endif
      end
      subroutine mcrep(icast,nr,nrepr,nc,nrepc,s,srep)
C- Replicas of a matrix
      implicit none
      integer nr,nc,icast,nrepr,nrepc
      double precision s(nr,nc,2),srep(nr,nrepr,nc,nrepc,2)
      integer i,j,ir,jr

      do  10  jr = 1, nrepc
      do  10  j = 1, nc
      do  10  ir = 1, nrepr
      do  10  i = 1, nr
        srep(i,ir,j,jr,1) = s(i,j,1)
        if (mod(icast,10) == 2) srep(i,ir,j,jr,2) = s(i,j,2)
   10 continue

      end
      subroutine mcroll(icast,nr,nrollr,nc,nrollc,s,sroll)
C- Shifts rows, columns of a matrix with cyclic endpoints
      implicit none
      integer nr,nc,icast,nrollr,nrollc
      double precision s(nr,nc,2),sroll(nr,nc,2)
      integer i,j,ir,jr

      do  while (nrollr<0)
        nrollr = nrollr+nr
      enddo
      do  while (nrollc<0)
        nrollc = nrollc+nc
      enddo

      do  j = 1, nc
        jr = (mod(j-1+nrollc,nc)+1)
        do  i = 1, nr
          ir = (mod(i-1+nrollr,nr)+1)
C         print 333, i,ir, j,jr
C 333     format(2i4,2x,2i4)
          sroll(ir,jr,1) = s(i,j,1)
          if (mod(icast,10) == 2) sroll(ir,jr,2) = s(i,j,2)
        enddo
      enddo

      end

      subroutine mcmsg(vsn,wksize)
      implicit none
      integer wksize,iout,i1mach
      double precision vsn
      iout = i1mach(2)
    1 format(2x,a,':',T20,a)
    2 format(2x,a,': ',a)
    3 format(T20,a)
      call awrit2(' usage: mc [-switches] (array) [ops] ... '//
     .  '(vsn %d  wksize=%ik)',' ',80,iout,vsn,wksize/1000)
      print *, '       fnam is an array name or a disk file name'
      print *, '...Switches:'
      print 1, '--version','display the version and exit'
      print 1, '-nc=# (nr=#)','stipulate that next matrix read'
     .             //' has # cols (rows)'
      print 1, '-vvar=#','define variable ''var'', assign value'
      print 1, '-show','show array stack and any operations pending'
      print 2, '-w[:l=string][:nohead] filnam','write top array to file filnam'
      print 2, '-bw[:l=string] filnam','ditto, for binary file'
      print 1, '-wap','write top array as (amplitude,phase)'
      print 2, '-a[:nr=#|:nc=#] nam','assign top array to nam'
      print 3, 'Optional nr=# redefines nr,nc'
      print 1, '-ap nam','same as -a, but keep top array on stack'
      print 1, '-av[:ir,ic] nam','assigns scalar variable ''nam'' to'//
     .         ' element from top array on stack'
C     print 1, '-qr','(quick-read) read fnam with fortran read'
      print 1, '-r:switches','switches for file reading'
     .            //' separated by '':''.  switches are'
      print 3, 'qr        read with fortran read (fast, no algebra)'
      print 3, 's=#       skips # records before reading'
      print 3, 'open      leaves file open after reading'
      print 3, 'br        read from binary file'
      print 3, 'spc       store in sparse format'
      print 3, 'br,nr,nc: file missing 1st record and nr,nc'//
     .         ' are supplied'
C      print 1, '-br','(binary-read) read binary file fnam'
C      print 3,  'NB: -br:nr,nc skips 1st binary record'
      print 1, '-px[:nprec]','(-pxc) write in row (column) compressed '
     .  //'storage format'
      print 3, 'showing elements exceeding tolerance'
      print 1
      print *, '...Unops (operate on top matrix)'
      print 1, '-p','pushes top array onto stack'
      print 1, '-p+n (-p-n)','pushes nth array (from bottom) '
     .            //' onto stack'
      print 1, '-pop','pops stack'
      print *, '-rsum[:lst[;lst2;lst3...]]'
      print *, '-csum[:lst[;lst2;lst3...]]'
      print 3, 'combines (sums) list (or a sequences of lists separated by ;) of rows (columns).'
      print 3, 'The result array has as many rows as there are sequences of lists.'
      print 1, '-s#','scale matrix by # (may use -s#real,#imag)'
      print 1, '-shft=#','add constant # to matrix'
      print 1, '-sort','expr sorts rows matrix according to expr'
      print 1, '-i,-iq','inverts matrix'
      print 1, '-lu','replaces hermitian s0 by L where s0 = L L+'
      print 1, '-1:n','pushes onto stack unit matrix, dimension n'
      print 1, '-tp [nc~]list','generates matrix from "list"'
      print 1, '-evl  (evc)','replaces matrix by its eigenvalues'
     .            //' (eigenvectors)'
      print 1, '-t','transposes matrix'
      print 1, '-cc','take complex conjugate of matrix'
      print 1, '-herm','forces matrix hermitian (or symmetric)'
      print 1, '-real','take the real part of a complex matrix'
      print 1, '-v2dia','expands vector into diagonal matrix'
      print 1, '-rep:n1,n2','concatenates replicas of matrix to create'
     .            //' (nr*n1,nc*n2) array'
      print 1, '-roll:#1[,#2]',
     .           'cyclically shifts rows (and columns, respectively)'
     .           //' by #1 (and #2) elements.'
      print 1, '-pwr=#','raises (col-1) matrix to a power'
      print 1
      print 1, '-rot=strn','specify a rotation matrix from elementary rotations'
      print 1, '-roteu','Performs the inverse operation of -rot='
      print 3, 'It converts a rotation matrix resident on the stack to Euler angles.'
      print 3, '  mc -rot=x:pi/2         pushes a new 3x3 array (rotation) on the stack while'
      print 3, '  mc -rot=x:pi/2 -roteu  returns the Euler angles for that rotation'
      print 1, '-ylm~switches','Makes rotation matrix for spherical harmonics, spinors, or both'
      print 3, 'from top-level stack array (which must be a 3x3 rotation matrix)'

      print 1
      print 2, '...Row and column unops:'
      print 1, '-rowl list'
      print 1, '-coll list','create array from list'
     .            //' of rows (cols) of top array'
      print 1, '-rowl:mirr'
      print 1, '-coll:mirr','rearrange rows (cols) in reverse order'
      print 1, '-rowl:pf=fnam'
      print 1, '-coll:pf=fnam','same as -rowl (-coll) but'
     .            //' list read from file'
      print 1, '-rowl:ipf=fnam'
      print 1, '-coll:ipf=fnam','same as -rowl (-coll) but permutations'
     .            //' inverted.'
      print 1, '       Example','perm file contains 2 3 1 4.'
      print *, ' -rowl:pf=perm    returns array w/ rows 2,3,1,4'
      print *, ' -rowl:ipf=perm   returns array w/ rows 3,1,2,4'
      print 1, '-inc expr','include rows that satisfy expr<>0'
      print *, ' -sub  t,b,l,r | -sub t,l'
      print *, ' -subx t,b,l,r | -sub t,l'
      print 3, 'extracts or excludes a subblock of top matrix'
      print 3, 'In second form, bottom right corner = (nr,nc)'
      print 1, '-subs:# t,b,l,r','scales a subblock of top matrix by #'
      print 1, '-subv:# t,b,l,r','copies # to subblock of top matrix'
      print 2, '-split nam x1,x2..xn  y1,y2,..yn','split into subblocks'
      print *, ' -array[#1,#2] #,#,...'
      print 3, 'Push onto the stack a new real array of dimension (#1,#2).'
      print 3, '#,#,... must contain of #1 x #2 expressions, which specify elements'
      print 3, '(1,1), (1,2), ..., (1,#2), (2,1), ..., (#1,#2).'
      print 2, '-sarray[#1,#2]',' Assemble a superarray of #1 x #2 arrays preexisting on the stack.'
      print 3, 'The top rows of superarray are constructed from the first #1 x #2 arrays'
      print 3, 'while the bottom rows are constructed from the last #1 x #2 arrays.'
      print 2, '-e# expr1 ... expr#','Create new matrix of # columns '
     .            //'with values expr1 ... expr#'
      print 3,  'evaluated from matrix currently top of stack'

      print 1
      print *, '...The following unops treat data as discretized',
     .             ' continuous functions of col1:'
      print *, '   Syntax for optional [opts] : [mesh],[rat],[nord=#]'
      print 2, '-intrp[:opts] list','interpolates col 2 at list of'//
     .          ' points; see above for opts'
      print 2, '-int[:opts] x0 list','integrates col 2 '//
     .             'from x0 to a list of upper-bounds'
      print 1, '-diff[:opts]','differentiates cols 2..nc wrt col1'
      print 2, '-smo width,xmin,xmax,dx','smooths vector of '//
     .          'delta-functions with gaussians'
      print 1, '-abs takes the absolute value of each element'
      print 1, '-max[:i|g]','puts the largest element into the first'
     .         //' column'
      print 3, 'Optional g returns max value of the entire array'
      print 3, 'Optional i returns index to max value'
      print 2, '-unx[:i1] #1,#2','(uncross) exchanges points in '//
     .         'columns #1 and #2'
      print 3, 'after their point of closest approach'
      print 2, '-unx2[:i1] #1,#2','(cross) exchanges points in '//
     .         'columns #1 and #2'
      print 3, 'at their point of closest approach, '//
     .  'if they do not cross'
      print 2, '-at[:i1] val expr','find adjacent rows that bracket '//
     .  'expr=val.'
      print 3, 'Array contains linearly interpolated x1 and expr'

      print 1, '-nint[#1]','Replaces each element with nearest integer.'
      print 3, 'Optional #1: scale array by #1 before operation, '//
     .  'and by 1/#1 after.'
      print 3, 'Thus -nint:1000 rounds to nearest thousanth.'
      print 1
      print *, '...Binops (operate on top two matrices s0 and s1)'
      print 1, '-tog','toggles top arrays on stack'
      print 1, '-+','adds s1 and s0'
      print 1, '--','adds s1 and -s0'
      print 1, '-x','multiplies s1 and s0'
      print 1, '-xe','multiplies s1 and s0 element by element'
      print 1, '-x3','multiplies s1 and s0 as 3D arrays:'
      print 3, 's1=s1(n11,n21,n31); s0=s0(n10,n20,n30);'
      print 3, 'where n10=nr(0)/n20,n20=nc(1),n30=nc(0);'//
     .         '  n11=n10,n21=nr(1)/n11,n31=n20'
      print 3, 'result(i,j,k) = sum_m s1(i,j,m) s0(i,m,k)'//
     .         ' is condensed into 2D (nr(1),nc(0))'
      print 1, '-de','divides s1/s0 element by element'
      print 1, '-gevl','like -evl, but s1 nonorthogality.'
      print 1, '-gevc','like -evc, but s1 nonorthogality.'
      print 1, '-orthos','replaces s0 with s1^-1/2 s0 s1^-1/2'
      print 1, '-ccat','concatenates columns'
      print 1, '-rcat','concatenates rows'
      print 1, '-cross','cross product s1(1,1..3) x s0(:,1..3)  '//
     .  'or  s1(:,1..3) x s0(:,1..3)'
      print *, ' -subc[:#] t,b,l,r | -subc[:#] t,l'
      print *, ' -suba[:#] t,b,l,r | -suba[:#] t,l'
      print 3, 'Copies or adds s1 to s0.  If # is present, copy or add #*s1 to s0'
      print 3, 'In second form, bottom right corner = (nr,nc)'
      print 1, '-index','Uses s0 as an row index to s1.'
      print 3, 's0(i) is overwritten by s1(s0(i)).  s1 is preserved.'
      print 3, 'New s0 has row dimensions of s0 and '//
     .         'column dimensions of s1'
      print 1
      print *, ' Compare numerical values in two files'
      print '(xx,a,a)',
     .  '-cmpf[~ln=#[,#][~col=#[,#][~sep=c][~quiet][~verb][~long][~char|~tol=#]',
     .  '[incl][excl][~nchar|~nword|~ndiff|~max]~fn1=f1~fn2=f2'
      print 1
      print *, '... Loop through sequence of command-line arguments'
      print 1, '[ list  arg1  arg2 ... ] or [ var=list  arg1  arg2 ... ]'
      print *, '   arg1 arg2 ... are parsed for each element in the list'
      print *, '   Within the list an argument can be ?expr'
      print *, '   If expr is false, the next argument is passed over'

!       call cexit(-1,1)
      end
      subroutine prmxx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
      implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(14) fmt, fmt0, strn*(*), outs*80
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#else
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(i,j,1),j=1,nc)
      enddo
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)

C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#else
      outs = ' prm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in prmx'
C#else
      if (outs == 'q') call rx0('quit in prmx')
C#endif
      return

      entry prmx0(fmt0)
      fmt = fmt0
      end

      subroutine zprm(strn,icast,s,ns,nr,nc)
C- Print complex matrix to file out
      implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(20) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      character fmt0*(*)
      save fmt
C     data fmt /'(9f15.10)'/
      data fmt /'(9f18.11)'/
C     data fmt /'(5f20.15)'/

      outs = ' '
      if (icast == 1)  outs = ' real'
      if (icast == 11) outs = ' symm'
      if (icast == 2)  outs = ' complex'
      if (icast == 12) outs = ' herm'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#else
      call awrit2('%% rows %i cols %i'//trim(outs),' ',80,ifi,nr,nc)
C#endif
      do  i = 1, nr
        write (ifi,fmt) (s(1,i,j),j=1,nc)
      enddo
      if (mod(icast,10) > 1) then
      write(ifi,*)
      do  i = 1, nr
        write (ifi,fmt) (s(2,i,j),j=1,nc)
      enddo
      endif
      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = 'prm: done writing data '//strn//' continue?'
C      print *, outs
C#else
      outs = ' zprm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in zprm'
C#else
      if (outs == 'q') call rx0('quit in zprm')
C#endif
      return

      entry zprm0(fmt0)
      fmt = fmt0
      end

      subroutine rdperm(first,nr,nc,np,nt,iperm)
Ci   np    :maximum allowed dimension of perm
Co   nt    :number of rows (columns) to (permute and) copy
Co         :(may not exceed np)
Co   perm  :array of rows or columns to copy
      implicit none
      character first*256
      integer nr,nc,np,nt

      double precision, allocatable:: perm(:)
      integer, allocatable:: jperm(:)
      double precision xx
      integer iperm(np)
      integer ifi,fopng,i,rdm,ncr
      logical ireverse
      character rowcol*4

      ireverse = .false.
      if (first(1:9) == '-rowl:pf=') then
        ifi = fopng(first(10:),-1,1)
        rowcol = 'rows'
        ncr = nr
      elseif (first(1:10) == '-rowl:ipf=') then
        ifi = fopng(first(11:),-1,1)
        rowcol = 'rows'
        ireverse = .true.
        ncr = nr
      elseif (first(1:9) == '-coll:pf=') then
        ifi = fopng(first(10:),-1,1)
        rowcol = 'cols'
        ncr = nc
      elseif (first(1:10) == '-coll:ipf=') then
        ifi = fopng(first(11:),-1,1)
        rowcol = 'cols'
        ireverse = .true.
        ncr = nc
      else
        call rx('rdperm: failed to parse '//trim(first))
      endif

      nt = 0
      i = rdm(ifi,1000,0,' ',xx,nt,1)
      allocate(perm(nt))
      rewind ifi
      i = rdm(ifi,1000,nt,' ',perm,nt,1)
      call fclose(ifi)
      if (nt > np) call rx('too many '//rowcol)
      if (ireverse) then
        if (nt /= ncr) call rx
     .    ('rdperm: size of perm file must match number of '//rowcol)
        allocate(jperm(nt))
        do  i = 1, nt
          jperm(i) = nint(perm(i))
        enddo
        do  i = 1, nt
          iperm(i) = i
        enddo
        do  i = 1, nt
          if (jperm(i) < 1 .or. jperm(i) > nt) call rx
     .      ('rdperm: permutation element out of range')
          iperm(jperm(i)) = i
        enddo
      else
        do  i = 1, nt
          iperm(i) = nint(perm(i))
        enddo
        deallocate(perm)
      endif

C     call info2(0,0,0,'# %n:1i',nt,iperm)

      end

      subroutine makerots(mode,nl,eula,rmat)
C- 2x2 spinor rotation part of rotation of Ylm
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode : 1s digit
Ci        : 0 do nothing
Ci        : 1 return spinor part of rotation matrix only
Ci        : 2 return ylm part of rotation matrix only
Ci        : 3 combination of 1 and 2
Ci        : 10s digit
Ci        : 1 rotate cubic to spherical harmonics
Ci        : 2 rotate cubic to spherical harmonics, m ordered l:-l
Ci   nl   : lmax+1 for rmat
Ci   eula : Euler angles
Co Outputs:
Co   rmat : rotation matrix
Cr Remarks:
Cr   rmat is block diagonal in l.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl
      double precision eula(3)
      double complex rmat(nl**2,2,nl**2,2)
C ... Local parameters
      integer nlm,i,j,is,js,mod0
      double precision eulal(3),rotm(3,3)
      double precision roty(nl**2,nl**2)
      double complex u(2,2,nl*nl),zroty(nl**2,nl**2)

      call dpzero(rmat,2*size(rmat))
      mod0 = mod(mode,10)
      if (mod0 == 0) return

!     call rm2eua(rotm,eulal(1),eulal(2),eulal(3))

      eulal = eula
      if (mod(mod0,2) == 0) eulal = 0
!     call rotspu(0,1,1,1,nl,eulal,1,u) ! Spinor rotation matrices
      call rotspu(1,1,1,1,nl,eulal,1,u) ! Spinor rotation matrices, passive rotation

      nlm = nl**2
      if (mod(mod0/2,2) > 0) then
        call eua2rm(eula(1),eula(2),eula(3),rotm)
        call ylmrtg(nlm,rotm,roty)
C       call prmxx('r(ylm)',roty,nlm,nlm,nlm,nlm)
        i = mod(mode/10,10)
        if (i > 0) then
          i = 2-i
          call s2sph(1+100*i,nl,nl,roty,nlm,nlm,nlm,nlm,zroty)
        else
          zroty = roty
        endif
C       call zprm('zroty',2,zroty,nlm,nlm,nlm)
      else
        call dpzero(zroty,2*size(zroty))
        forall(i=1:nlm) zroty(i,i) = 1
      endif

C ... Combine spin and orbital parts
      do  is = 1, 2
      do  js = 1, 2
        do  i = 1, nlm
        do  j = 1, nlm
          rmat(i,is,j,js) = u(is,js,i) * zroty(i,j)
        enddo
        enddo
      enddo
      enddo

      end

      subroutine clrtri(icast,n,a)
C- Clear upper triangle of matrix a
      integer icast
      double precision a(n,n,2)

      do  i = 1, n
        forall(j=1:i-1) a(j,i,1) = 0
      enddo

      if (mod(icast,10) <= 1) return
      do  i = 1, n
        forall(j=1:i-1) a(j,i,2) = 0
      enddo

      end
      subroutine ipoly(n,xitrp,xn,wk)
C- Coefficients to lagrange interpolating polynomial with n points
      implicit none
      integer n
      real(8) :: xitrp,xn(n),wk(n)

      integer i,j

      do  j = 1, n
        wk(j) = 1
        do  i = 1, n
          if (j == i) cycle
          wk(j) = wk(j)*(xitrp-xn(i))/(xn(j)-xn(i))
        enddo
      enddo

      call dcopy(n,wk,1,xn,1)

      end
      subroutine ludcmp(a,n,icast)
      implicit none
      integer n,icast
      double complex a(n,n)

      integer info,i,j

      if (mod(icast/10,10) /= 1)
     .  call fexit2(-1,1,'ludcmp not ready for nonhermitian matrix; sorry',1,2)

      do  i = 1, n
      do  j = i+1, n
        a(i,j) = 0
      enddo
      enddo

C     call zprm('a',2,a,n,n,n)
      call zpotrf('L', n, a, n, INFO )
c     call zprm('a',2,a,n,n,n)

      end
