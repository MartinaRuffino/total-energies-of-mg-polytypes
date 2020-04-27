C Nonlinear least-squares fit by linearization about initial guess.
C Example: nlfit a1=1 a2=1.7 "a1*cos(a2*x1)" a1,a2 nlfit.dat
C   finds a1,a2 that minimize error in nlfit.dat
C Weights if used, go in last column; use w=1/y^2 to fit to rel err.
C Note: at Query prompt, beta can be changed interactively.
C       Setting beta=0 makes nlfit print current values and exit,
C       or if -tp is set, nlfit prints results of fit and exits.
      subroutine fmain
      implicit none
      integer mxdim,mxbas,mxdat,ssvdef,ncf,garg,ip,iter
C     parameter (mxdim=40,mxbas=12,mxdat=1000,ssvdef=512)
      parameter (mxdim=40,ssvdef=512)
      logical fail(3)
      real(8), dimension(:,:), allocatable :: dat,datt,a,q
      real(8), dimension(:), allocatable :: w,res,c,y,f0,fnew
      double precision coffs(mxdim),r(mxdim*mxdim),
     .  stdj,xx,av(20),del,beta,xxv(50),errmx,ddot,vmax,d1mach,xmax
      integer it(50),nrt,nct
      character*(ssvdef) f,sydef,first,tlist,sube,s,xn*4
      character f2(ssvdef)*1, cnams(20)*4,outs*128
      logical cmdstr,lsequ,rddim,a2bin,rdwgt
      EXTERNAL cmdstr
      integer iarg,n,n1r,n1c,ifi,i1mach,i,ich,ich0,nvar,
     .  awrite,rdm,fopng,nargf,nfmt,iv0
      integer nfcn,n1,ipivot(mxdim),ifault,iprint
      equivalence (f2,first)

      common /static/ coffs,r

C#ifdef DEBUG
      double precision snorm,sd,ss,z
      integer m1,kz,nw,mode,j,k,l,ndf,m1p1,n1p1,k2,m
      data m1 /0/, kz /0/, l /1/
C#endif

      nfmt = 6
      nrt = 0
      nct = 0
      call initqu(.true.)
      rdwgt = .false.
      sydef = ' '
      sube  = ' '
      n1r = 0
      n1c = 0
      ifi = 10
      beta = 1d0
      goto 10
   20 print 333
      print 334
  333 format(' usage: nlfit [-switches] [initial-vals ..]  ',
     .  'f  var1,var2,... fnam'/
     .  8x,'Least-squares fit of analytic function to data.'/
     .  8x,'The function may depend in a nonlinear way on coefficients;'
     .  /
     .  8x,'fitting proceeds by linearization of coefficients about',
     .     ' trial values.'//
     .  8x,'Data (both dependent and independent variables) are',
     .     ' contained in file fnam;'/
     .  8x,'the columns of fnam are labelled x1, x2, ...'/
     .  8x,'Each row in fnam holds a data point; the number of rows',
     .     ' in fnam are the number of'/
     .  8x,'points that will be included in the least-squares fit.'/
     .  /
     .  8x,'The analytic function is represented as an ascii string f,',
     .     ' which consists of an'/
     .  8x,'algebraic expression comprised',
     .     ' of coefficients and variables.'/
     .  8x,'Variables (which depend on file data) are labelled x#,'/
     .  8x,'e.g. x1, x2, where #=i corresponds to the ith column of',
     .     ' fnam.'/
     .  8x,'The least-squares procedure sums over all data points',
     .     ' (rows in fnam).'//
     .  8x,'Coefficients can be given any name not',
     .     ' beginning with x.'/
     .  8x,'  Example for f:  "a1+a3*cos(2*a2*x1)"'/
     .  8x,'Any number of coefficients may be used in specifying f:'/
     .  8x,'each coefficient must be assigned an initial value',
     .     ' (see below).'//
     .  8x,'Not all coefficients need be varied.'/
     .  8x,'The user specifies which coefficients are to be',
     .     ' included'/
     .  8x,'in the fit by listing them after f.  For example:'/
     .  8x,'    nlfit a0=1 a1=1 a2=1.7 "a0+a1*cos(a2*x1)" a1,a2 fnam'/
     .  8x,'fits f(a0,a1,a2;x)=a0+a1*cos(a2*x1) to',
     .     ' data in file fnam,'/
     .  8x,'varying a1 and a2 to minimize the difference between f',
     .     ' and the dependent variable.'/
     .  8x,'In this example,  a0=1 a1=1 a2=1.7  assign',
     .     ' initial values to coefficients;'/
     .  8x,'the dependent variable y is taken from the last column',
     .     ' of fnam'/
     .  8x,'(This is default for the dependent variable; it may be'/
     .  8x,' overridden, as described in Switches below)'
     .  )

  334 format(/
     .  8x,'Switches:'/
     .  10x,'-pr# set the verbosity (higher # => more printout)'/
     .  10x,'-b=beta mixes  (1-beta)*starting + beta*estimate'/
     .  10x,'        (on query, vary beta.  0 => exit, dumping vals'/
     .  10x,'-w uses the last column of fnam as fit weights'/
     .  10x,'-tp [nc]~xlist ... tabulates fit at table of points'/
     .  10x,'   nc=0 (default) maps points to 2 cols: xi=point, yi=0'/
     .  10x,'   nc=1  puts points in a single column'/
     .  10x,'   nc=n  orders points as dat(1,1:nc); dat(2,1:nc)...'/
     .  10x,'   nc=-n orders points as dat(1:nr,1); dat(1:nr,2)...'/
     .  10x,'-sub var=val (or "var1=val1 var2=val2 ...")'/
     .  10x,'   subexpressions evaluated prior to evaluation of  f'/
     .  10x,'   to simplify representation of  f'/
     .  10x,'-ifile reads some initial vals from file "file"'/
     .  10x,'-y= defines the dependent variable as an algebraic',
     .      ' function of x1..xn'/
     .  10x,'-g# specifies the precision coffs are printed out'/
     .  )
      call cexit(-1,1)

   10 continue
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (garg('-pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      elseif (garg('-b=',iarg,4,' ',1,1,i,j,beta) /= 0) then
      else if (first(1:4) == '-tp ') then
        if (garg(' ',iarg+1,4,',:~ ',4,50,i,it,xxv) < 1) goto 22
        iarg = iarg+1
        call expand(0,i,it,xxv,datt,nrt,nct)
        allocate(datt(nrt,2*nct))
        call expand(1,i,it,xxv,datt,nrt,nct)
        call pshpr(0)
      else if (lsequ(first,'-g',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),nfmt,2,0,' ',i,-1)) goto 20
      else if (lsequ(first,'-i',2,' ',n)) then
        ifi = 10
        j = 0
        call chrpos(first,' ',100,j)
        open(ifi, file=first(3:j), status='OLD', err=20)
   24   first = ' '
        read(ifi,'(a)',err=22, end=22) first
        i = 0
        call parsyv(first,len(first),99,1,i)
        goto 24
   22   continue
      elseif (garg('-nc=',iarg,2,' ',1,1,i,it,n1c) /= 0) then
        if (i < 0) goto 20
      else if (lsequ(first,'-y=',3,' ',n)) then
        if (.not. cmdstr(iarg,sydef)) goto 20
      else if (lsequ(first,'-sub ',5,' ',n)) then
        if (.not. cmdstr(iarg+1,sube)) goto 20
        iarg = iarg+1
      else if (lsequ(first,'-w',2,' ',n)) then
        rdwgt = .true.
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15
   30 continue

C --- Evaluate possible assignment; continue if one  ---
      iarg = iarg+1
      j = 0
      call chrpos(first,'=',40,j)
C ... If an expression, this is function to fit!
      if (j < 40) then
        n = j+1
        k = 0
        call parsyv(first,len(first),100,1,k)
        goto 15
      endif

C --- Input f, coffs ---
      f = first
      if (.not. cmdstr(iarg,first)) goto 20
      ich  = 0
      ich0 = 0
      nfcn = 0
   19 call chrps2(first,', ',2,ssvdef,ich,i)
      nfcn = nfcn+1
      cnams(nfcn) = first(ich0+1:ich)
      if (i == 1) then
        ich0 = ich+1
        ich = ich0
        if (ich < ssvdef) goto 19
      endif
      iarg=iarg+1

C --- Read data from file ---
      if (.not. cmdstr(iarg,first) .and. nargf() < iarg) goto 20
      if (first == '.' .or. nargf() == iarg) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif
C     Get dimensioning
      i = rdm(ifi,0,0,' ',xxv,n1r,n1c)
      mxdat = n1r
      mxbas = n1c
      allocate(dat(mxdat,mxbas+1),a(mxdat,mxdim),q(mxdat,mxdim))
      allocate(w(mxdat),res(mxdat),y(mxdat),f0(mxdat),fnew(mxdat),
     .  c(4*(mxdat+mxdim)+2))

      rewind ifi
      i = rdm(ifi,0,mxbas*mxdat,' ',dat,n1r,n1c)
      if (i < 0)
     .  call fexit(-1,9,'nlfit failed to parse file '//first,0)
      ncf = n1c
      if (rdwgt) ncf = n1c-1
      if (rdwgt .and. n1c < 3) stop 'nlfit: missing weights col'
      if (n1c*n1r > mxbas*mxdat) stop 'lfit: input array too big'
      if (ifi /= i1mach(1)) call fclose(ifi)
      if (iprint() >= 20) then
        call awrit3(' nlfit: read %i rows, %i cols.  Fit to col %i.',
     .    ' ',80,i1mach(2),n1r,n1c,ncf)
      endif

C --- Generate normal matrix ---
      iter = 1
      del = d1mach(3)**(1/3d0)
   40 continue
      call nrmgen(cnams,sube,f,ssvdef,sydef,dat,n1r,ncf,nfcn,iter,
     .  rdwgt,av,del,a,y,f0,w)

C --- Least-squares fit and standard deviations ---
      call l2a(n1r,nfcn,0,1,a,y,w,0d0,n1r,nfcn,
     . n1, ipivot, coffs, res, r, q, c, ifault)
C ... NB: l2a replaced w by its square root!
      do  42  i = 1, n1r
   42 res(i) = res(i)*w(i)
      xx = ddot(n1r,res,1,res,1)
      if (iprint() > 0) then
        print *
        call awrit5(' New coffs:  fault=%i  ncoff=%i  np=%i  rank=%i'//
     .    '  RMS= %1,3;3g(est)',first,80,i1mach(2),ifault,nfcn,n1r,n1,
     .    dsqrt(xx/n1r))
        call awrit1(' Coff      Old        Diff        Mixed       '//
     .    'Std dev   (beta=%d)',first,80,i1mach(2),beta)
      endif
      if (ifault == 11) then
        print *, ' (warning) rounding error => diagonal cov. elt < 0'
        ifault = 0
      endif
      if (ifault /= 0 .or. n1 /= nfcn) call cexit(-1,1)
      do  50  j = 1, nfcn
        stdj = dsqrt(xx*r(j+(j-1)*nfcn))
        if (iprint() > 0)
     .  print 320, cnams(j),av(j),coffs(j),av(j)+beta*coffs(j), stdj
  320   format(1x,a4,1p4e15.6)
c       outs = ' '//cnams(j)
C        call info5(1,0,0,'%a   %12;6g  %12;6g  %12;6g %12;6g',
C        call awrit4('%a %,6;6g',outs,len(outs),i1mach(2),
C        call awrit5('%a  %,6;6g   %,6;6g   %,6;6g   %,6;6g',
C        call awrit4('%a  %;12,6D   %,6;6g   %,6;6g   %,6;6g',
C     .    outs,80,-i1mach(2),
C     .    av(j),coffs(j),av(j)+beta*coffs(j), stdj)
C        call awrit4('%a',outs,len(outs),i1mach(2),
C     .    av(j),coffs(j),av(j)+beta*coffs(j), stdj, 0)


   50 continue

      call daxpy(nfcn,beta,coffs,1,av,1)

      if (iprint() >= 40) then
      call nrmgen(cnams,sube,f,ssvdef,sydef,dat,n1r,ncf,nfcn,iter,
     .  rdwgt,av,del,a,y,fnew,w)

      print 459
  459 format (/2x,
     .  'i   Observed        Starting       Predicted       Residual')
      errmx = 0d0
      vmax  = 0d0
      ss    = 0d0
      do  60  i = 1, n1r
        xx = dat(i,ncf)
        if (sydef /= ' ') then
          call numsyv(iv0)
          call clrsyv(iv0)
          call addsyv('i',dble(i),k)
          call addsyv('x',dat(i,1),k)
          do  124  j = 1, ncf
            k = 1
            xn = 'x   '
            call bin2a('',0,0,j,2,0,4,xn,k)
            call addsyv(xn,dat(i,j),k)
  124     continue
          j = 3
          if (.not. a2bin(sydef,xx,4,0,' ',j,-1)) then
            call rx('a2bin failed to parse y function')
          endif
          call clrsyv(iv0)
        endif
        vmax = max(vmax,dabs(xx))
        print 458, i, xx,f0(i),fnew(i),fnew(i)-xx
  458   format (i3,3g16.8,g12.4)
  340   continue
        ss = ss + (fnew(i)-xx)**2*w(i)
        errmx = max(errmx,dabs(fnew(i)-xx))
   60 continue
      snorm = dsqrt(ss/n1r)
      call awrit3(' true RMS err = %1,3;3g  Max err = %1,3;3g  RMS '//
     .  'err/Max val = %1,3;3g',first,80,i1mach(2),snorm,errmx,
     .  snorm/vmax)
      endif

      if (nrt /= 0 .and. nct /= 0) then
C   ... This call to nrmgen needed to produce fnew
        call nrmgen(cnams,sube,f,ssvdef,sydef,datt,nrt,nct,nfcn,iter,
     .    rdwgt,av,del,a,y,fnew,w)
        call awrit2('%% rows %i cols %i',' ',80,i1mach(2),nrt,nct+1)
        do  70  i = 1, nrt
        print 345, (datt(i,j), j=1,nct), fnew(i)
  345   format(10g14.6)
   70   continue
C        call awrit0('%%',' ',80,i1mach(2))
        if (beta <= 0) call fexit(0,0,' ',0)
      endif

C --- Allow beta to change; if beta=0, print variables and exit ---
      iter = iter+1
      call query('beta=',4,beta)
      if (beta <= 0) then
        call numsyv(nvar)
        ip = -1
        s = ' '
        do  80  i = 4, nvar
          call watsyv(s(ip+2:ssvdef),xx,i)
          ip = awrite('%a=%1;ng',s,ssvdef,0,nfmt,xx,0,0,0,0,0,0)
   80   continue
        print 349, s(1:ip)
  349   format(a)
        s = ' '
        xx = ddot(n1r,res,1,res,1)
        do  82  j = 1, nfcn
          call awrit0('%a d'//cnams(j),s,ssvdef,0)
          stdj = dsqrt(xx*r(j+(j-1)*nfcn))
          ip = awrite('%a=%1;6g',s,ssvdef,0,stdj,0,0,0,0,0,0,0)
   82   continue
        print 349, s(2:ip)
        call fexit(0,0,' ',0)
      endif
      goto 40
      end
      subroutine nrmgen(cnams,sube,f,ssvdef,sydef,dat,np,ncf,nfcn,iter,
     .  rdwgt,av,del,norm,y,y0,w)
C- Normal matrix by linearization of coffs around starting values
C ----------------------------------------------------------------------
Ci Inputs
Ci   cnams :names of parameters to be varied
Ci   sube  :subexpressions evaluated prior to evaluation of f
Ci   f     :algebraic, ascii representation of analytic function
Ci   ssvdef:not used
Ci   sydef :algebraic, ascii representation of the dependent variable
Ci         :as a function of x1..xn
Ci   dat   :data to be fit
Ci   np    :number of data points
Ci   ncf   :column in dat defining independent variable
Ci   nfcn  :number of parameters to fit
Ci   iter  :iteration number (iter=1 ...?)
Ci   rdwgt :T, then a points have a weighting factor
Ci   av    :current parameter values
Ci   del   :"small" parameter for numerical differentiation
Ci   w     :fitting weights
Co Outputs
Co   y     :difference between fit and data
Co   y0    :fit function values
Co   norm  :normal matrix
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   01 Aug 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer np,nfcn,ssvdef,iter,ncf
      double precision dat(np,1),norm(np,nfcn),y(np),y0(np),w(np),val,
     .  fv(20,3),d(20),del
      character*(*) f,sube,cnams(*),sydef
C ... Local parameters
      character*4 an,xn
      logical rdwgt,a2bin
      integer ii,ival,iv0,j,ip,k,iprint
      double precision av(20)

C --- Reset variables with new point for linearization ---
      do  22  j = 1, nfcn
        ii = 1
        an = cnams(j)
        if (iter == 1) call getsyv(an,av(j),ival) ! First pass, use table value
        call chsyv(an,av(j),ival)
        call rxx(ival == 0,'LINSV: variable '//an//' not found')
   22 continue
      if (iprint() >= 50) call shosyv(0,0,0,6)

      call numsyv(iv0)

C --- For each data point do ---
      do  45  ip = 1, np

        call clrsyv(iv0)

C ---   Load data variables from data array ---
        w(ip) = 1
        if (rdwgt) w(ip) = dat(ip,ncf+1)
        call addsyv('i',dble(ip),ival)
        call addsyv('x',dat(ip,1),ival)
        do  24  j = 1, ncf
          ii = 1
          xn = 'x   '
          call bin2a('',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ip,j),ival)
   24   continue
        y(ip) = dat(ip,ncf)
        if (sydef /= ' ') then
          ii = 3
          if (.not. a2bin(sydef,y(ip),4,0,' ',ii,-1)) then
            call rx('a2bin failed to parse y function')
          endif
        endif
        call addsyv('y',y(ip),ival)
        if (iprint() > 50) call shosyv(0,0,0,6)

C ---   Linearize around origin ---
        do  44  j = 1, nfcn
          ii = 1
          an = cnams(j)

C ...   Do for av(j), av+, av-
          do  46  k = 1, 3

            if (k == 1) call chsyv(an,av(j)+del,ival)
            if (k == 2) call chsyv(an,av(j)-del,ival)
            if (k == 3) call chsyv(an,av(j)    ,ival)

            ii = 0
            if (sube /= ' ')  then
              call parsyv(sube,len(sube),999,1,ii)
C             call shosyv(0,0,0,6)
              ii = 0
            endif
            if (.not. a2bin(f,val,4,0,' ',ii,-1)) then
              print *, f
              call shosyv(0,0,0,6)
              call rx('a2bin failed to parse function')
            endif
            fv(j,4-k) = val
   46     continue

          norm(ip,j) = (fv(j,3) - fv(j,2))/(2*del)
C         print *, ip,j, norm(ip,j), fv(j,1)
   44 continue
      y0(ip) = fv(1,1)
      y(ip) = y(ip) - fv(1,1)
   45 continue

C --- Printout ---
      if (iprint() >= 100) then
        print *, 'nrmgen: normal matrix'
        call prm('(2f12.5)',norm,np,nfcn)
        print *, 'nrmgen: rhs'
        call prm('(2f12.5)',y,np,1)
      endif
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
      subroutine expand(opt,nt,it,xxv,dat,nr,nc)
C- Expands a list of points into a table
C   it : indices to that flag meaning of arguments in xxv
Ci  it(i) = 1 => delimiter = ','
Ci  it(i) = 2 => delimiter = ':'
Ci  it(i) = 2 => delimiter = '~'
Ci  it(i) = 4 => delimiter = ' '
      implicit none
      integer opt,nt,it(nt),nr,nc,ir,i,ipr,ic,nc0,i1mach
      double precision dat(10),xxv(nt),dx,x0,x1
      double precision xnint

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
      if (opt /= 0) dat(ir) = xxv(i)
C ... Generate points in range
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
          if (opt /= 0) dat(ir) = x0
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
          if (opt /= 0) call dpzero(dat(nr+1),nr)
          return
        elseif (xxv(1) > 0 .and. opt /= 0) then
          call dcopy(nr*nc,dat,1,dat(1+nr*nc),1)
          call xxpand(nr,nc,dat(1+nr*nc),dat)
        endif
        if (ipr >= 80 .and. opt /= 0) then
          call awrit2('%% rows %i cols %i',' ',80,i1mach(2),nr,nc)
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
