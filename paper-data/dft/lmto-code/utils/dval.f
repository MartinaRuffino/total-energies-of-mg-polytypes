C dval evaluates algebraic expressions, sending results to stdout
C Syntax:
C   dval [ -fformat | -a# aformat ]  [switches] [nam=expr | expr ...]
C The result of any expr is printed to stdout.  Thus
C    dval pi/2 pi/3 pi/4
C prints out
C    1.5708 1.0472 0.785398
C the string nam=expr evaluates expr and assigns to nam.  It does not
C print the result of expr, but 'nam' can be used in later expressions.
C Thus
C     dval 'x=atan(1)'
C prints out nothing, but
C    dval 'x=atan(1)' 'y=sin(-x*2)' x y y^2
C prints out
C    0.785398 -1 1
C see file a2bin.f for a table of operators.
C
C macros:
C   Macros may be defined, which act in a manner similar to the
C   operation of a function.  Example: compare analytic second derivative
C   of g(x)=e^(-x^2/rs^2).  Note x is parenthesized in the macro definition.
C      dval -m'g(x)=exp(-((x)^2/rs^2))' x=2 rs=2.2 \
C      '(g(x+1d-4)+g(x-1d-4)-2*g(x))/(1e-4)^2' '2/rs^2*(2*x^2/rs^2-1)*g(x)'
C   Evaluate smoothed exponential:
C   dval -m'env1(x)=exp(kap*(x))*erfc(kap*rs/2+(x)/rs)+exp(-kap*(x))*erfc(kap*rs/2-(x)/rs)' rs=1.1 kap=1.2 'env1(2)'
C   This generates 0.177273
C
C Formatting the output is accomplished with -fformat or -a# aformat.
C format looks like d#, e# or g#, which format the result in,
C respectively, a fortran-like F format, a fortran-like 1p,E format,
C or, either of the above, whichever is shorter.
C       Example           precision   printout
C    dval -fd3 pi/1000    absolute     0.003
C    dval -fd9 pi/1000    absolute     0.003141593
C    dval -fe9 pi/1000    relative     3.14159265e-3
C    dval -fg9 pi/1000    relative     0.00314159265
C Alternatively, you can use -a# aformat.  This version formats the
C string according to the awrite formatter (see awrite.f).  Example
C   dval -a "B(%1;3d)=%1;3d ep=%e epp=%1,3;3g"  .12 .023456 .034567 13
C prints out
C   B(0.12)=0.023 ep=3.4567e-2 epp=13.0
C
C Switches:
C The only switch implemented so far is the solution of simultaneous
C equations.  Syntax is -s # "z1=expr1 z2=expr2 .. zn=exprn", where
C expr1..n are expressions of x1..xn (and any other variables).
C dval attempts to find x1..xn that make z1..zn equal to zero.  Example:
C   dval -a2 'x1=%,9;9d x2=%,9;9d and z1,z2=%;10,10d %;10,10d'\
C   x1=0 x2=0 -s 2 "z1=x1+2*x2+6 z2=2*x1+1*x2*x2+2" x1 x2 z1 z2
C prints out:
C   x1=-2.516685226 x2=-1.741657387 and z1,z2=0.0000000000 0.0000000000
C Updates
C   19 Jul 2011  better treatment of -i switch
C ----------------------------------------------------------------------
      subroutine fmain
      implicit none
      double precision s(100),val(100),pi
      character first*1024,fmt*40,outs*256,awrits*256
      character*1 f2(20)
      logical cmdstr,lsequ,a2bin
      integer iwk(100),iarg,n,i,j,k,ix,ixa,naw,ifi,garg,ip,ips
      integer nsim,ssvdef,ntry
      parameter (ssvdef=256)
      character*(ssvdef) svdef
      equivalence (f2,first)
      double precision tol,xv(10),fv(30),dx,damp,ddot,err
      integer last,iprint,retval
      integer nxargc,nargf,iargsave(2)

      tol = 1d-10
      awrits = ' '
      outs = ' '
      ips = 0
      damp = 1
      nsim = 0
      ntry = 30
      dx = 1d-5
      pi = 4*datan(1d0)
      call addsyv('pi',pi,n)
      ix = 0
      ixa = 0
      naw = 0
      retval=0
      goto 10
   20 print 333
  333 format(
     .' usage: dval [-fformat -a# fmt -ifilename -s args]',
     .  '  expr | var=expr  ... '/
     .' -a fmt uses awrit# format to display output'/
     .' -mmacro(args)=expr',
     ."    e.g. dval -m'g(x)=exp(-((x)^2/rs^2))' rs=1.1 'g(2)'"/
     .' -s nvar[:damp:nit:dx] "expr ..." to solve simultaneous eqns:'/
     .'    x1 .. xn are independent var; z1 .. zn dependent var;'/
     .'    "expr ..." must set zn .. zn using x1 .. xn')
      call cexit(-1,1)

   10 continue
      fmt = 'g'
      iarg = 1
      if (.not. cmdstr(iarg,first)) goto 20
   15 continue
C     print *, 'iarg,nxargc()=', iarg, nxargc()
      if (iarg >= nxargc()) then
        if (iargsave(1) /= 0) then
          iarg = iargsave(1)+1
          iargsave(1) = 0
        else
          goto 80
        endif
      endif
      if (iarg == nargf() .and. iargsave(1) == 0) then
        goto 80
      endif
      if (.not. cmdstr(iarg,first)) goto 80

      if (.not. lsequ(first,'-',1,' ',n)) goto 30

      if (lsequ(first,'-f',2,' ',n)) then
        fmt = ' '
        call strcat(fmt,1,' ',f2(3),99,' ',n)
      elseif (first(1:2) == '-m') then
        j = 0
        call chrpos(first,'=',len(first)-1,j)
        if (j. lt. len(first)-1) first(j+1:j+1) = ' '
        call macset(first(3:),j)
        if (j. lt. 0) call rxs('dval: bad macro def: ',first)
      elseif (lsequ(first,'- ',2,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,first)) goto 20
        goto 30
      else if (first(1:2) == '-v') then
        j = 2
        call parsyv(first,len(first),999,0,j)
      else if (lsequ(first,'-i',2,' ',n)) then
        ifi = 10
        j = 0
        call chrpos(first,' ',100,j)
        open(ifi, file=first(3:j), status='OLD', err=20)
   24   read(ifi,'(a)',err=22, end=22) first
        iargsave(1) = iarg
        iargsave(2) = nxargc()
        rewind ifi
        do  j = 1, 9999
          first = ' '
          read(ifi,'(a)',err=23, end=23) first
          call acmdop(first,len(first),0)
C         call acmdop(first,len(first),1)
        enddo
   23   continue
        if (nxargc() > iargsave(2)) then
          iarg = iargsave(2)-1
C          do  i = 1, nxargc()
C            if (.not. cmdstr(i-1,first)) goto 20
C            print *, i, first(1:20)
C          enddo
        endif
        goto 24
   22   continue
      elseif (garg('-pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      elseif (garg('-a',iarg,2,' ',1,1,i,j,naw) /= 0) then
        if (naw == 0) naw = 8
        iarg = iarg+1
        if (.not. cmdstr(iarg,awrits)) goto 20
C      else if (lsequ(first,'-a',2,' ',n)) then
C        call bin2a(first(3:40),1,0,0,1,0,255,outs,ips)
C --- Simultaneous equations ---
      else if (lsequ(first,'-s',2,' ',n)) then
        ip = 0
C#ifndef DEC
        if (garg('',iarg+1,4,': ',2,4,i,iwk,xv) < 1) goto 20
C#elseC
C        if (garg(' ',iarg+1,4,': ',2,4,i,iwk,xv) < 1) goto 20
C#endif
        nsim = nint(xv(1))
        if (i >= 2) damp = xv(2)
        if (i >= 3) ntry = nint(xv(3))
        if (i >= 4) dx = xv(4)
        if (.not. cmdstr(iarg+2,svdef)) goto 20
        iarg = iarg + 3
        if (nsim > 0) then
          call rxx(nsim > 10,'dval: nsim gt 10')
          last = 0
          do  40  j = 1, ntry
            call linsv(nsim,xv,svdef,ssvdef,j,dx,s,fv)
            err = dsqrt(ddot(nsim,fv,1,fv,1)/nsim)
            call dgefa(s,nsim,nsim,fv(11),i)
            call rxx(i /= 0,'dval: normal matrix singular')
            first = '  x:'
            ip = 4
            do  44  i = 1, nsim
   44       call bin2a(' ',1,0,xv(i),4,0,255,first,ip)
            call bin2a('z:',2,0,0,1,0,255,first,ip)
            do  45  i = 1, nsim
   45       call bin2a(' ',1,0,fv(i),4,0,255,first,ip)
            call bin2a('dx:',2,0,0,1,0,255,first,ip)
            call dgesl(s,nsim,nsim,fv(11),fv,0)
            do  46  i = 1, nsim
   46       call bin2a(' ',1,0,-damp*fv(i),4,0,255,first,ip)
            if (iprint() >= 40) then
              call cwrite(outs,1,ips-1,1)
              print *, 'it ',j, first(1:ip)
            endif
            if (last == 1) goto 15
            err = dsqrt(ddot(nsim,fv,1,fv,1)/nsim)
            call daxpy(nsim,-damp,fv,1,xv,1)
            if (err < tol) last = 1
   40     continue
        endif
        goto 15
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15

C --- End of switches ---
   30 continue

C --- Evaluate expressions ---
   72 continue
      iarg = iarg+1
      j = 0
      call chrpos(first,'=',40,j)
C ... If an expression, print it out
      if (j >= 40) then
        j = 0
        if (.not. a2bin(first,val,4,0,' ',j,-1)) then
          retval = 1
          goto 80
        endif
        ix = ix+1
        if (awrits /= ' ') then
          ixa = ixa+1
          val(ixa+1) = val(1)
C          if (naw > 0 .and. naw == ixa) then
C            stop 'now'
C          endif
        else
          call bin2a(fmt,1,0,val,4,0,255,outs,ips)
        endif
      else
        n = j+1
        k = 0
        call parsyv(first,len(first),100,1,k)
      endif
      goto 15

   80 continue
      if (awrits /= ' ') then
C ... fix ... cludge for now
C want -a to mean -a8 (ok) -a#: only that many for awrite
C ixa should reset and dump if new -a encountered
        call dumpaw('%a'//awrits,val(2),outs(2:len(outs)),len(outs)-1)
        call skpblb(outs,len(outs),ips)
        ips = ips+3
      else
      endif
      call cwrite(outs,1,ips-1,1)
      if (retval == 1) call cexit(-1,1)

      end
      subroutine dumpaw(fmt,v,sout,mxlen)
      implicit none
      integer naw,mxlen
      double precision v(10)
C      character *(*), fmt,sout
      character fmt*80,sout*80

      call awrit8(fmt,sout,mxlen,0,v(1),v(2),v(3),v(4),
     .  v(5),v(6),v(7),v(8))

      end

      subroutine linsv(n,xv,svdef,ssvdef,iter,del,s,fv)
C- Linearize functions svdef around xv
      implicit none
      integer n,ssvdef,iter
      character*(*) svdef
      double precision xv(n),fv(n,3),s(n,n),del
      character*4 xn,zn
      integer ipr,ival,iv0,i,j,k,ii
      double precision val

      call numsyv(iv0)
      call getpr(ipr)

C ---   Load data table ---
      do  22  j = 1, n
        ii = 1
        xn = 'x   '
        call bin2a('',0,0,j,2,0,4,xn,ii)
        if (iter == 1) call getsyv(xn,xv(j),ival)
        call chsyv(xn,xv(j),ival)
        call rxx(ival == 0,'LINSV: variable '//xn//' not found')
C       call shosyv(0,0,0,6)
   22 continue

C --- Find values z_i and partial(z_i)/partial(x_j) ---
      do  44  j = 1, n

        ii = 1
        xn = 'x   '
        call bin2a('',0,0,j,2,0,4,xn,ii)

C --- Do for xv(j), xv+, xv-
        do  46  k = 1, 3

          if (k == 1) call chsyv(xn,xv(j)+del,ival)
          if (k == 2) call chsyv(xn,xv(j)-del,ival)
          if (k == 3) call chsyv(xn,xv(j)    ,ival)

          ii = 0
          call parsyv(svdef,ssvdef,999,1,ii)
C         call shosyv(0,0,0,6)
          do  45  i = 1, n
            ii = 1
            zn = 'z   '
            call bin2a('',0,0,i,2,0,4,zn,ii)
            call getsyv(zn,val,ival)
            if (ival == 0) call rx('LINSV: missing variable '//zn)
            fv(i,4-k) = val
c           print *, zn,val,ival
   45     continue
   46   continue
        do  47  i = 1, n
   47   s(i,j) = (fv(i,3) - fv(i,2))/(2*del)
   44 continue

      if (ipr >= 100) then
      do  54  i = 1, n
   54 print 345, (s(i,j),  j=1,n)
      print 345, (fv(j,1), j=1,n)
  345 format(5f15.10)
      pause
      endif

      end
