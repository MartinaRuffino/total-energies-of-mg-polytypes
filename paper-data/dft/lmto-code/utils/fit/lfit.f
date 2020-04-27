C Least-squares fit of an arbitrary linear combination of functions zn
C to columns of tabulated data.  By default, the last column of data
C represents the function being fit to (or the second-to-last column
C if the -w switch used).  (To change default, use "-y=expr".)
C Tabulated data is referred to as (x1,x2,...xn).  "x" is an alias
C for "x1" and "y" is an alias for the independent function.
C Example: lfit 2 "z1=x z2=x*exp(-x)" dat
C   fits a1*x + a2*x*exp(-x) to data in "dat."
C Example: lfit 6 "z1=1 z2=x z3=z2*x z4=z3*x z5=z4*x z6=z5*x" dat
C   fits a fifth-degree polynomial to last column of dat.
C Example: lfit 5 "z1=1 z2=x z3=x^2 z4=-y*x z5=-y*x^2" dat
C   approximately fits dat to (a1+a2*x+a3*x^2)/(1+a4*x+a5*x^2)
C   With the resulting 5 coffs, check fit with e.g.
C   nlfit a1= a2= a3= a4= a5=  \
C         "(a1+a2*x+a3*x^2)/(1+a4*x+a5*x^2)" a1,a2,a3,a4,a5 dat
C Example: lfit "-y=x2-.5" 4 "z1=x z2=x^2 z3=-x2*x z4=-x2*x^2" dat
C   approximately fits dat to (.5+a1*x+a2*x^2)/(1+a3*x+a4*x^2)
C Weights if used, go in last column; use w=1/y^2 to fit to rel err.
      subroutine fmain
      implicit none
      integer mxdim,mxbas,mxdat,ssvdef,ncf
      parameter (mxdim=40,mxbas=12,mxdat=2000,ssvdef=256)
      logical fail(3)
      real(8), allocatable :: dat(:),a(:),q(:),c(:)
      double precision errmx,w(mxdat),coffs(mxdim),res(mxdat),r(mxdim*mxdim),y(mxdat),xx,ddot,d1mach
      character*72 first
      character*(ssvdef) svdef,sydef,outs,out2
      character*1 f2(20),fmt*20
      logical cmdstr,lsequ,a2bin,rdwgt
      EXTERNAL cmdstr
      integer iarg,n,n1r,n1c,ifi,i1mach,i,garg,iprint,
     .  rdm,fopng,nargf,a2vec
      integer nfcn,n1,ipivot(mxdim),ifault
      integer it(50),nfmt
C For tabulating fit on a mesh ...
      integer nmsh,ix(10)
      double precision xxmin(mxbas),xxmax(mxbas),dx(mxbas),xv(10)
      logical ltab
      equivalence (f2,first)

      common /static/ coffs,w,res,r


C#ifdef DEBUG
      double precision snorm,ss,z
      integer m1,kz,nw,j,k,l,ndf,m1p1,n1p1,k2,m
      data m1 /0/, kz /0/, l /1/
C#endif

      allocate(dat((mxbas+1)*mxdat),a(mxdim*mxdat),q(mxdim*mxdat),c(4*(mxdat+mxdim)+2))

      rdwgt = .false.
      sydef = ' '
      ltab = .false.
      n1r = 0
      n1c = 0
      ifi = 10
      nfmt = 0
      goto 10
   20 print 333
  333 format(' usage: lfit [-switches] nfcn fcns fnam'/
     .  1x,'Least-squares fit of data in file "fnam" to ',
     .     'a linear combination'/1x,'of nfcn functions fcns.',
     .     '  The first columns contain the'/1x,'independent variables',
     .     ' x1,x2,..., and the last contains the'/
     .  1x,'dependent variable y.'/
     .  1x,'Switches:'/
     .  3x,'-pr# to set the verbosity'/
     .  3x,'-vnam=val sets variable "nam" to val'/
     .  3x,'-w uses last column as fit weights (the penultimate',
     .      ' column is used for y)'/
     .  3x,'-y defines the dependent variable as a function of x1..xn'/
     .  3x,'-g# specifies the precision coffs are printed out'/
     .  3x,'-t[n] x1,x2,dx [x1,x2,dx ..] tabulates fit data on a mesh'/
     .  3x,'-tc2 x1,x2,dx x1,x2,dx  x1,x2 data in fplot contour format'/
     .  1x, 'Examples:'/
     .  1x, 'Fit a fifth-degree polynomial to last column of "dat":'/
     .  3x,'lfit 6 "z1=1 z2=x z3=z2*x z4=z3*x z5=z4*x z6=z5*x" dat'/
     .  1x, 'Approximately fit dat to (a1+a2*x+a3*x^2)/(1+a4*x+a5*x^2)'/
     .  3x,'lfit 5 "z1=1 z2=x z3=x^2 z4=-y*x z5=-y*x^2" dat')
      call cexit(-1,1)

   10 continue
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (garg('-pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      else if (lsequ(first,'-v',2,' ',n)) then
        j = 2
        first = first // ' '
        call chrpos(first,'=',19,j)
        n = j+1
        if (.not. a2bin(first,xx,4,0,' ',n,-1)) goto 20
        call addsyv(first(3:j),xx,n)
      else if (lsequ(first,'-t',2,' ',n)) then
        nmsh = 1
        if (first == '-tc2') then
          nmsh = -2
        elseif (first(3:3) /= ' ') then
          i = 0
          if (.not. a2bin(f2(3),nmsh,2,0,' ',i,-1)) goto 20
        endif
        do  j = 1, iabs(nmsh)
          iarg = iarg+1
          if (.not. cmdstr(iarg,first))
     .      call rx('lfit: not enough arguments following -t ...')
          i = 0
          i = a2vec(first,len(first),i,4,', ',2,3,3,ix,xv)
          if (i /= 3) call rx('lfit failed to parse '//trim(first))
          xxmin(j) = xv(1)
          xxmax(j) = xv(2)
          dx(j)    = xv(3)
        enddo
        ltab = .true.
      else if (lsequ(first,'-g',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),nfmt,2,0,' ',i,-1)) goto 20
      elseif (garg('-nc=',iarg,2,' ',1,1,i,it,n1c) /= 0) then
        if (i < 0) goto 20
      else if (lsequ(first,'-y=',3,' ',n)) then
        if (.not. cmdstr(iarg,sydef)) goto 20
      else if (lsequ(first,'-w',2,' ',n)) then
        rdwgt = .true.
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15
   30 continue

C --- Input nfcn, svdef ---
      read(first,'(i4)',err=20) nfcn
      iarg=iarg+1
      if (.not. cmdstr(iarg,svdef)) goto 20
      iarg=iarg+1

C --- Read data from file ---
      if (.not. cmdstr(iarg,first) .and. nargf() < iarg) goto 20
      if (first == '.' .or. nargf() == iarg) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif
      i = rdm(ifi,0,mxbas*mxdat,' ',dat,n1r,n1c)
      if (i < 0)
     .  call fexit(-1,9,'lfit failed to parse file '//first,0)
      ncf = n1c
      if (rdwgt) ncf = n1c-1
      if (rdwgt .and. n1c < 3) stop 'lfit: missing weights col'
      if (n1c*n1r > mxbas*mxdat) stop 'lfit: input array too big'
      if (ifi /= i1mach(1)) call fclose(ifi)
      if (iprint() >= 5) then
        call awrit3('# lfit: read %i rows, %i cols.  Fit to col %i.',
     .    ' ',80,i1mach(2),n1r,n1c,ncf)
      endif

C --- Generate normal matrix ---
      call nrmgen(svdef,sydef,ssvdef,dat,n1r,ncf,nfcn,rdwgt,a,y,w)

C --- Least-squares fit ---
      call l2a(n1r,nfcn,0,1,a,y,w,0d0,n1r,nfcn,
     . n1, ipivot, coffs, res, r, q, c, ifault)

C --- Tabulate fit if switch set ---
      if (ltab) then
        call tfit(nmsh,xxmin,xxmax,dx,svdef,sydef,ssvdef,nfcn,coffs,res)
        stop
      endif

C --- Print computed results ---
      nw = i1mach(2)
      m = n1r
      l = 1
  210 continue
      print *
      call awrit4(' fault=%i  nfcn=%i  rank=%i  RMS error= %;3g',first,
     .  80,i1mach(2),ifault,nfcn,n1,dsqrt(ddot(n1r,res,1,res,1)/n1r))
C     write (nw,99890) ifault, nfcn, n1
      if (ifault >= 1 .and. ifault <= 4) go to 390
      if (n1 == 0) go to 390
      if (n1 < m1) go to 390
C     write (nw,99870)
C     write (nw,99860) (ipivot(j),j=1,n1)
      if (n1 == nfcn) go to 220
      n1p1 = n1 + 1
      write (nw,99850)
      write (nw,99860) (ipivot(j),j=n1p1,nfcn)
  220 ndf = m - n1 - kz
      do 240 k=1,l
        k2 = l + k
        if (c(k) < 0.0d0) go to 230
        fail(k) = .false.
c        write (nw,99830) c(k),c(k2)
        go to 240
  230   fail(k) = .true.
        c(k) = -c(k)
        write (nw,99820) k,c(k),c(k2)
  240 continue
      if (ifault == 7) go to 390
      do 350 k=1,l
        if (fail(k)) go to 350
        ss = 0.0d0
        if (m == m1) go to 310
        m1p1 = m1 + 1
        errmx = 0d0
        do  i=m1p1,m
          ss = ss + (res(i)*w(i))**2
          errmx = max(errmx,dabs(res(i)))
        enddo
  310   write (nw,99790)
        j = -dlog(d1mach(3))/dlog(10d0)
        if (nfmt /= 0) j = nfmt
        call awrit2('%x(i6,g%i.%i)',fmt,len(fmt),0,j+10,j)
        do 320 j=1,nfcn
  320   write (nw,fmt) j,coffs(j)
        if (iprint() >= 30) then
        write (nw,99780)
        do  i=1,m
          z = y(i) - res(i)
          write (nw,99770) i, y(i),z,res(i), (a(i+(j-1)*n1r), j=1,nfcn)
        enddo
        endif
        snorm = dsqrt(ss)
        write (nw,99760) ss,snorm,errmx
C       if ((mode == 2 .and. n1 < n) .or. (m == n1+kz)) go to 350
C       write (nw,99750) sd
  350 continue


      outs = ' '
      j = 10
      if (nfmt /= 0) j = nfmt
      call awrit1('%%;%ig',fmt,len(fmt),0,j)
      do  j = 1, nfcn
        call awrit2('%a a%i='//trim(fmt),outs,len(outs),0,j,coffs(j))
      enddo

C     outs = trim(out2) // ' ' // trim(svdef)
      out2 = svdef
      i = len(trim(out2)) + 1
      do  j = 1, nfcn
        call awrit1('%aa%i*%-1jz%i+',out2(i+1:),len(out2)-i,0,j)
      enddo
      i = len(trim(out2))
      out2(i:i) = ' '
      call info0(1,1,0,' To evaluate, use e.g.%N  dval x=.. \')

      print '(1x,a)', trim(outs)
      print '(2x,a)', trim(out2)


C      if (beta <= 0) then
C        call numsyv(nvar)
C        ip = -1
C        s = ' '
C        do  80  i = 4, nvar
C          call watsyv(s(ip+2:ssvdef),xx,i)
C          ip = awrite('%a=%1;nd',s,ssvdef,0,nfmt,xx,0,0,0,0,0,0)
C   80   continue
C        print 349, s(1:ip)
C  349   format(a)
C        s = ' '
C        xx = ddot(n1r,res,1,res,1)
C        do  82  j = 1, nfcn
C          call awrit0('%a d'//cnams(j),s,ssvdef,0)
C          stdj = dsqrt(xx*r(j+(j-1)*nfcn))
C          ip = awrite('%a=%1;6d',s,ssvdef,0,stdj,0,0,0,0,0,0,0)
C   82   continue
C        print 349, s(2:ip)
C        call fexit(0,0,' ',0)
C      endif



c
c print lower triangular portion of symmetric unscaled covariance
c matrix.
c
c      if (mode == 2 .and. n1 < n) go to 390
c      write (nw,99740)
c      if (mode == 2) go to 370
c      do 360 i=1,nfcn
c        write (nw,99910) (r(i,j),j=1,i)
c  360 continue
c      go to 390
c  370 do 380 i=1,nfcn
c        write (nw,99910) (qr(i,j),j=1,i)
c  380 continue
c  390 read (nr,99730) ifdone
c
c
  390 continue

c format statements.
99910 format (1x,8g15.8)
99890 format (' fault=',i2,'  nfcn=',i2,'  computed rank=',i3)
99870 format (50h columns of h = (sqrt(w))*a were selected by the p,
     * 37hivoting scheme in the following order)
99860 format (30i4)
99850 format (50h the following columns of h are linearly dependent,
     * 48h.  if mode 1, they did not enter the regression./6h if mo,
     * 24hde 2, they entered last./)
99830 format (' Converged in ',f4.0, ' iter(s) starting from',
     .  ' estimated ',g15.8,' digits')
99820 format (i7,9x,18hfailed to converge,5x,f4.0,10x,g15.8)
99810 format (26h solution for b-vector no.,i3)
99800 format (1h ,27x,18hstandard deviation/5x,1hj,4x,10hcoefficien,
     * 4ht(j),5x,17hof coefficient(j)/)
99790 format (/2x,'function',2x,'Coefficient')
99780 format (/2x,
     .  'i   Observed        Predicted       Residual    x ...')
99775 format (i6,g25.15)
99770 format (i5,2g16.8,g12.4,2x,15g12.6:/'   ...',7g12.6)
99760 format (30h sum of squared residuals    =,g15.8/1x,8hnorm of ,
     * 9hresiduals,11x,1h=,g15.8/' Max error',19x,'=',g15.8)
99750 format (30h residual standard deviation =,g15.8)
99740 format (27h unscaled covariance matrix/)
99730 format (i5)
99720 format(25h number of zero weights =,i3,5x,17hdeg. of freedom =,i3)

      end
      subroutine tfit(nmsh,x1,x2,dx,svdef,sydef,ssvdef,nfcn,coffs,res)
C- Tabulation of fit data
C  nmsh = -2 => contour format
      implicit none
      integer nfcn,ssvdef
      character*(*) svdef,sydef
      double precision x1(1),x2(1),dx(1),coffs(1)
      integer nmsh
C     integer nbas,npoly,ptab(1),nmsh
      double precision x(20),res,xx
      integer ixx(20)
      integer j,ifi,fopng

c --- Print out the number of points along each coordinate ---
      do  j = 1, iabs(nmsh)
        ixx(j) = 0
        do  xx = x1(j), x2(j), dx(j)
          ixx(j) = ixx(j)+1
        enddo
      enddo
      print 334, (ixx(j), j=1, iabs(nmsh)), iabs(nmsh)+1
  334 format('% rows ',10i4)
      ifi = -999
      if (nmsh == -2) then
        ifi = fopng('con.dat',-1,0)
        write(ifi,333), ixx(2), ixx(1)  ! Transpose is written!
  333   format('% rows ',i4,' cols ',i4)
      endif

      call tfitx(ifi,iabs(nmsh),iabs(nmsh),x,x1,x2,dx,svdef,sydef,
     .  ssvdef,nfcn,coffs,res)

      end
      recursive subroutine tfitx(ifi,nmsh,j,x,x1,x2,dx,svdef,sydef,
     .  ssvdef,nfcn,coffs,res)
      implicit none
      integer ifi,nfcn,ssvdef  ! ifi >= 0 => write res to file ifi
      character*(*) svdef,sydef
      double precision x(nmsh),x1(nmsh),x2(nmsh),dx(nmsh),coffs(1),res
      integer nmsh,j,k
      double precision xx,a(20),ddot

      do  xx = x1(j), x2(j), dx(j)
        x(j) = xx
        if (j == 1) then
C ...     y and w mean nothing in this context
          call nrmgen(svdef,sydef,ssvdef,x,1,nmsh,nfcn,
     .      .false.,a,res,res)
          res = ddot(nfcn,a,1,coffs,1)
          print 333, (x(k),k=1,nmsh), res
C          if (x(j)+dx(j) > x2(j)) print *, 'hi!'
  333     format(10f16.8)
          if (ifi > 0) write(ifi,*) res
        else
          call tfitx(ifi,nmsh,j-1,x,x1,x2,dx,svdef,sydef,ssvdef,nfcn,
     .      coffs,res)
        endif
      enddo
      end

      subroutine nrmgen(svdef,sydef,ssvdef,dat,np,ncf,nfcn,rdwgt,
     .  norm,y,w)
C- Normal matrix for linear least-squares fit
      implicit none
      logical rdwgt,a2bin
      integer np,nfcn,ssvdef,ncf,iprint
      character*(*) svdef,sydef
      double precision dat(np,1),norm(np,nfcn),y(1),w(np),val
      integer ival,ic1,iv0,ip,idim,j,ii
      character*4 xn,zn

      ic1 = ichar('1')
      call numsyv(iv0)

C -- for each row of data ... --
      do  20  ip = 1, np
        w(ip) = 1
        if (rdwgt) w(ip) = dat(ip,ncf+1)
        call clrsyv(iv0)
        call addsyv('i',dble(ip),ival)
        call addsyv('x',dat(ip,1),ival)
        y(ip) = dat(ip,ncf)
        call addsyv('y',y(ip),ival)
C ---   Load data table ---
        do  22  j = 1, ncf
          ii = 1
          xn = 'x   '
          call bin2a('',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ip,j),ival)
   22   continue
C ---   Change y= if sought ---
        if (sydef(1:1) /= ' ') then
          j = 3
          if (.not. a2bin(sydef,y(ip),4,0,' ',j,-1))
     .      call rx('NRMGEN: parse error reading y')
          call chsyv('y',y(ip),ival)
        endif
        j = 0
C ---   Obtain functions z as expressions of declared vars ---
        call parsyv(svdef,ssvdef,999,1,j)
c       call shosyv(0,0,0,6)
        do  25  idim = 1, nfcn
          zn(2:2) = char(ic1+idim-1)
          ii = 1
          zn = 'z   '
          call bin2a('',0,0,idim,2,0,4,zn,ii)
          call getsyv(zn,val,ival)
          if (ival == 0) call rx('NRMGEN: missing variable'//zn)
c         print *, zn,val,ival
          norm(ip,idim) = val
   25   continue
        if (iprint() >= 100) call shosyv(0,0,0,6)
        if (iprint() >= 100) pause
   20 continue

      call clrsyv(iv0)
      end

      subroutine prm(fmt,s,nr,nc)
      double precision s(nr,nc)
C#ifdef APOLLO_BUG
      character*(20) fmt
C#elseC
C      character*(*) fmt
C#endif
      print *, nr, nc
      do  10  i = 1, nr
   10 print fmt, (s(i,j), j=1,nc)
      end
