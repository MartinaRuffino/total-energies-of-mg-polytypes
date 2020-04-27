C Nonlinear least-squares fit by linearization about initial guess.
C Example: nlfitlm a1=1 a2=1.7 'a1*cos(a2*x1)' a1,a2 nlfit.dat
C   finds a1,a2 that minimize error in nlfit.dat
C Weights if used, go in last column; use w=1/y^2 to fit to rel err.
C Note: at Query prompt, beta can be changed interactively.
C       Setting beta=0 makes nlfit print current values and exit,
C       or if -tp is set, nlfit prints results of fit and exits.
      subroutine fmain
      implicit none
      integer mxdim,mxbas,mxdat,ssvdef,ncf,garg,ip,iter
      parameter (mxdim=40,mxbas=12,mxdat=1000,ssvdef=512)
      double precision
     .  w(mxdat),coffs(mxdim),res(mxdat),r(mxdim*mxdim),q(mxdim*mxdat),
     .  stdj,xx,d1mach,alsc,alam,sig(mxdat),chi(2),chi2,wk(mxdim,5),
     .  del,beta,xxv(50),datt(mxbas*mxdat),errmx,ddot,vmax
      integer,allocatable:: ivar(:)
      real(8),allocatable :: alp(:,:),alp2(:,:),beta2(:),cov(:,:),a(:,:)
      real(8),allocatable :: dat(:,:),fy(:),fnew(:),y(:)
      real(8),allocatable :: av(:),avsav(:),av0(:)

      integer it(50),nrt,nct,j,k,ntab,mkdlst
      character*(ssvdef) f,sydef,ssdef,first,sube,s,s2,ptslst,savefile
      character f2(ssvdef)*1, cnams(20)*4
      logical cmdstr,lsequ,a2bin,rdwgt
      EXTERNAL cmdstr
      integer iarg,iargnam,iargfnam,n,n1r,n1c,ifi,i1mach,i,ich,ich0,
     .  nvar,awrite,rdm,fopng,nargf,nfmt,maxit
      integer nfcn,iprint
      equivalence (f2,first)
      integer namlen
      parameter (namlen=16)
      character an*(namlen)
      double precision rmsdel

      common /static/ coffs,datt,w,res,q,r

      nfmt = 6
      nrt = 0
      nct = 0
      call initqu(.true.)
      rdwgt = .false.
      sydef = ' '
      ssdef = ' '
      sube  = ' '
      ptslst = ' '
      savefile = ' '
      n1r = 0
      n1c = 0
      ifi = 10
      maxit = 0
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
     .  10x,'-pr#      set the verbosity (higher # => more printout)'/
     .  10x,'-y=expr   defines the dependent variable in terms of an',
     .      ' algebraic function of x1..xn'/
     .  10x,'-sig=expr uses the last column of fnam as fit weights'/
     .  10x,'-t=list   tabulates fit at a list of points,'/
     .  10x,'          varying x1 through the list.'/
     .  10x,'          (valid only for functions of a single variable)'/
     .  10x,'-sub var=val (or "var1=val1 var2=val2 ...")'/
     .  17x,'   subexpressions evaluated prior to evaluation of  f'/
     .  17x,'   to simplify representation of  f'/
     .  10x,'-ifile    reads some initial vals from file "file"'/
     .  10x,'-maxit=#  set lambda to 0 after maxit iterations'/
     .  10x,'-g#       specifies the precision coffs are printed out'/
     .  10x,'-save=fn  On exit (set lambda=0), append fit to data;',
     .  ' save in file fn.'/
     .  )
      call cexit(-1,1)

   10 continue
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (garg('-pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      elseif (garg('--pr',iarg,2,' ',1,1,i,j,k) /= 0) then
        call pshpr(k)
      elseif (garg('-maxit=',iarg,2,' ',1,1,i,j,maxit) /= 0) then
      else if (first(1:3) == '-t=') then
        ptslst = first(4:)
      else if (first(1:4) == '-tp ') then
        if (garg(' ',iarg+1,4,',:~ ',4,50,i,it,xxv) < 1) goto 22
        iarg = iarg+1
        call expand(i,it,xxv,datt,nrt,nct)
        call pshpr(0)
      else if (lsequ(first,'-g',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),nfmt,2,0,' ',i,-1)) goto 20

        i = 0
        if (.not. a2bin(f2(3),nfmt,2,0,' ',i,-1)) goto 20
      else if (first(1:6) == '-save=') then
        savefile = first(7:)
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
      else if (first(1:5) == '-sig=') then
        if (.not. cmdstr(iarg,ssdef)) goto 20
      else if (lsequ(first,'-sub ',5,' ',n)) then
        if (.not. cmdstr(iarg+1,sube)) goto 20
        iarg = iarg+1
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

C --- read string for f, names of coffs to vary ---
      f = first
      iargnam = iarg
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
      iargfnam = iarg
      if (.not. cmdstr(iarg,first) .and. nargf() < iarg) goto 20
      if (first == '.' .or. nargf() == iarg) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif

C     Number of variables declared
      call numsyv(nvar)

C --- Read data from file or stdin ---
      i = rdm(ifi,0,mxbas*mxdat,' ',datt,n1r,n1c)
      if (i < 0)
     .  call fexit(-1,9,'nlfit failed to parse file '//first,0)
      ncf = n1c
      if (rdwgt) ncf = n1c-1
      if (rdwgt .and. n1c < 3) stop 'nlfit: missing weights col'
      if (n1c*n1r > mxbas*mxdat) stop 'lfit: input array too big'
      if (ifi /= i1mach(1)) call fclose(ifi)
      if (iprint() >= 20) then
        call awrit3(' nlfitml: read %i rows, %i cols.  Fit to col %i.',
     .    ' ',80,i1mach(2),n1r,n1c,ncf)
      endif
      allocate(dat(n1r,n1c),y(n1r),fy(n1r),fnew(n1r))
      call dcopy(n1r*n1c,datt,1,dat,1)
      allocate(alp(nfcn,nfcn),cov(nfcn,nfcn),ivar(nfcn))
      allocate(alp2(nfcn,nfcn),beta2(nfcn))
      allocate(av(nfcn),avsav(nfcn),av0(nfcn))
      allocate(a(nfcn,n1r))

C --- Initial values for coefficients to vary ---
      do  j = 1, nfcn
        ivar(j) = j  ! all listed coefficients are to be varied
        an = cnams(j)
        call getsyv(an,av(j),i)
        call rxx(i == 0,'NLFITML: variable '//an//' not found')
        av0(j) = av(j)
        avsav(j) = av(j)
      enddo

C --- Start of iteration loop ---
      del = d1mach(3)**(1/3d0)
      alsc = 5
      alam = 1
      iter = 0
      call dpzero(fy,n1r)

C ... Entry for next iteration
   40 continue
C     Function values, chi for new trial cofs
      call gradc(cnams,sube,f,sydef,ssdef,dat,n1r,ncf,nfcn,
     .  av,del,y,fnew,a,sig)
C     This call only used to get chi2 for current coffs
      call mrqcof(n1r,nfcn,0,sig,y,fnew,a,chi2,xxv,xxv)
      rmsdel = dsqrt(chi2/n1r)

      call info5(10,1,0,' Coffs, iter %i:  ncoff=%i'//
     .  '  old chi^2=%1,3;3g  new chi^2=%1,3;3g  lambda=%;3g',
     .  iter,nfcn,chi,chi2,alam)
      call info0(10,0,0,' Coff%09pOld%24pChange%39pMixed%54pVariance')

      do  j = 1, nfcn
        stdj = dsqrt(cov(j,j))
        if (iprint() > 0)
     .  print 320, cnams(j),avsav(j),av(j)-avsav(j),av(j), stdj
  320   format(1x,a4,1p4e15.6)

      enddo

      if (iprint() >= 40) then
      print 459
  459 format (/2x,
     .  'i   Observed        Old fit        New fit         Residual')
      endif
      errmx = 0d0
      vmax  = 0d0
      do  i = 1, n1r
        vmax = max(vmax,dabs(y(i)-fnew(i)))
        if (iprint() >= 40) then
        print 458, i, y(i),fy(i),fnew(i),fnew(i)-y(i)
  458   format (i3,3g16.8,g12.4)
        endif
        errmx = max(errmx,dabs(fnew(i)-y(i)))
      enddo
      call info5(10,0,0,' RMS err = %1,3;3g  Max err = %1,3;3g  '//
     .  'RMS err/Max val = %1,3;3g',rmsdel,errmx,rmsdel/vmax,0,0)

C --- Allow user to change lambda; if lamda=0, print variables and exit ---
      if (rmsdel < 1d-10) alam = 0
      if (maxit > 0 .and. iter == maxit) alam = 0
      call query('lambda=',4,alam)
      if (alam <= 0) then

        alam = 0
        call dcopy(nfcn,av,1,avsav,1)
        call dcopy(n1r,fnew,1,fy,1)
        call mrqstp(n1r,nfcn,nfcn,ivar,y,fy,a,sig,wk,av,
     .    alsc,alp,alam,cov,chi,iter)
        call info0(10,1,0,' Coff%09pStart%24pNow%37pstd deviation')
        do  j = 1, nfcn
          stdj = dsqrt(chi(1)*cov(j,j))
          if (iprint() > 0)
     .      print 321, cnams(j),av0(j),av(j), stdj
  321     format(1x,a4,1p,3e15.6)
        enddo

        first = 'nlfitlm'
        if (maxit /= 0) then
          call awrit1('nlfitlm -maxit=%i',first,80,0,maxit)
        endif
        ip = len(trim(first))

        do  i = 4, nvar
          call watsyv(first(ip+2:ssvdef),xx,i)
          ip = awrite('%a=%1;ng',first,ssvdef,0,nfmt,xx,0,0,0,0,0,0)
        enddo
        s = trim(first) // " '" // trim(f) // "'"
        if (.not. cmdstr(iargnam,first)) goto 20
        if (.not. cmdstr(iargfnam,s2)) goto 20
        print 349, 'New command:',trim(s),trim(first),trim(s2)
  349   format(/a/a,1x,a,1x,a)

C        s = ' '
C        xx = ddot(n1r,res,1,res,1)
C        do  82  j = 1, nfcn
C          call awrit0('%a d'//cnams(j),s,ssvdef,0)
C          stdj = dsqrt(xx*r(j+(j-1)*nfcn))
C          ip = awrite('%a=%1;6d',s,ssvdef,0,stdj,0,0,0,0,0,0,0)
C   82   continue
C        print 349, 'Standard deviations:',s(2:ip)

C   ... Append fit to data; write to file
        if (savefile /= ' ') then
          call info0(10,1,0,'Writing data to file '//trim(savefile))
          deallocate(alp)
          allocate(alp(n1r,n1c+1))
          alp(:,1:n1c) = dat
          alp(:,n1c+1) = fnew
          ifi = fopng(savefile,-1,0)
          rewind ifi
          call ywrm(0,' ',1,ifi,'(1p20e16.8)',alp,0,n1r,n1r,n1c+1)
          call fclose(ifi)
        endif

C   ... Tabulate fit on mesh of points (function of 1 variable only)
        if (ptslst /= ' ') then
          ntab = mkdlst(ptslst,1d-8*maxval(abs(dat(:,1))),-1,xx)
          deallocate(dat,y,fnew)
          allocate(dat(ntab,1),y(ntab),fnew(ntab))
          ntab = mkdlst(ptslst,1d-8*maxval(abs(dat(:,1))),ntab,dat)
          call gradc(cnams,sube,f,' ',' ',dat,ntab,1,0,
     .      av,del,y,fnew,a,sig)

          write(*,543) ntab,2
  543     format('% rows',i4,' cols',i2)
          do  i = 1, ntab
            print 544, dat(i,1),fnew(i)
  544       format(1p,2e20.10)
          enddo
        endif

        call fexit(0,0,' ',0)

      endif

C ... New trial coefficients
      call dcopy(nfcn,av,1,avsav,1)
      call dcopy(n1r,fnew,1,fy,1)
      call mrqstp(n1r,nfcn,nfcn,ivar,y,fy,a,sig,wk,av,
     .  alsc,alp,alam,cov,chi,iter)

      goto 40
      end

      subroutine gradc(cnams,sube,f,sydef,ssdef,dat,np,ncf,nvar,
     .  av,del,depv,y0,grady,sig)
C- Function values and gradients wrt coefficients, for all data points
C ----------------------------------------------------------------------
Ci Inputs
Ci   cnams :names of parameters to be varied
Ci   sube  :subexpressions evaluated prior to evaluation of f
Ci   f     :algebraic, ascii representation of analytic function
Ci   sydef :algebraic, ascii representation of the dependent variable y
Ci         :as a function of x1..xn.  If missing, y is is assigned
Ci         :to the last column of dat
Ci   ssdef :algebraic, ascii representation of the standard deviations
Ci         :sig (used as fitting weights).  If missing, sig=1
Ci   dat   :data to be fit
Ci   np    :number of data points
Ci   ncf   :number of columns in data table
Ci   nvar  :number of coefficients to vary
Ci   av    :current parameter values (av initialized, first iteration)
Ci   del   :"small" parameter for numerical differentiation
Co Outputs
Co   depv  :values of the dependent variable to fit
Co   y0    :function values
Co   grady :gradient matrix, d(y0)/d(cof)
Co   sig   :standard deviations, if specified by ssdef
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   01 Aug 09  Adapted from nrmgen (very similar)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer np,nvar,ncf
      double precision dat(np,1),grady(nvar,np),depv(np),y0(np),
     .  av(nvar),fv(nvar,3),sig(np),del
      character*(*) f,sube,cnams(*),sydef,ssdef
C ... Local parameters
      character*4 an,xn
      logical a2bin
      integer ii,ival,iv0,j,ip,k,iprint
      double precision val

C --- Reset symbolic variables table with new coefficients ---
      do  j = 1, nvar
        ii = 1
        an = cnams(j)
        call chsyv(an,av(j),ival)
        call rxx(ival == 0,'LINSV: variable '//an//' not found')
      enddo
      if (iprint() >= 50) call shosyv(0,0,0,6)

      call numsyv(iv0)

C --- For each data point do ---
      do  ip = 1, np

        call clrsyv(iv0)

C   ... Load data variables from data array
        call addsyv('i',dble(ip),ival)
        call addsyv('x',dat(ip,1),ival)
        do  j = 1, ncf
          ii = 1
          xn = 'x   '
          call bin2a('',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ip,j),ival)
        enddo
        depv(ip) = dat(ip,ncf)
        if (sydef /= ' ') then
          ii = 3
          if (.not. a2bin(sydef,depv(ip),4,0,' ',ii,-1)) then
            call rx('a2bin failed to parse y function')
          endif
        endif
        call addsyv('y',depv(ip),ival)
        if (iprint() > 50) call shosyv(0,0,0,6)

C   ... Standard deviations
        sig(ip) = 1
        if (ssdef /= ' ') then
          ii = 5
          if (.not. a2bin(ssdef,sig(ip),4,0,' ',ii,-1)) then
            call rx('a2bin failed to parse sigma function')
          endif
        endif

C   ... Linearize around origin
        do  j = 1, nvar
          ii = 1
          an = cnams(j)

C     ... Do for av(j), av+, av-
          do  k = 1, 3

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
          enddo
          grady(j,ip) = (fv(j,3) - fv(j,2))/(2*del)
C         print *, ip,j, grady(j,ip), fv(j,1)
        enddo
        if (nvar == 0) then
C         call shosyv(0,0,0,6)
          ii = 0
          if (sube /= ' ')  then
            call parsyv(sube,len(sube),999,1,ii)
            ii = 0
          endif
          if (.not. a2bin(f,val,4,0,' ',ii,-1)) then
            print *, f
            call shosyv(0,0,0,6)
            call rx('a2bin failed to parse function')
          endif
          fv(1,1) = val
        endif
        y0(ip) = fv(1,1)
      enddo

C --- Printout ---
      if (iprint() >= 100) then
        print *, 'nrmgen: normal matrix'
        call prm('(2f12.5)',grady,nvar,np)
        print *, 'nrmgen: rhs'
        call prm('(2f12.5)',depv,np,1)
      endif
      call clrsyv(iv0)

      end

      subroutine prm(fmt,s,nr,nc)
      integer nr,nc,i,j
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
