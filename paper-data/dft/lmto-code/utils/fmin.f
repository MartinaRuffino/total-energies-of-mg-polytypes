C#define MULTIMF
Cr Minimization a function of several variables.  Program interacts
Cr with an external programs to generate fcn call in one of two
Cr ways.
Cr
Cr Mode 1: fmin and some function are called in alternate
Cr succession, until convergence is reached.  In this mode, fmin
Cr interacts with the function via a file 'fmin.val.'  At each step,
Cr fmin appends a string of variables at the values it wants next to
Cr the end of the file.  The function appends the function value
Cr corresponding to those variables to the end of the file.  fmin
Cr works from the beginning each time, playing itself out like the
Cr previous time, but advancing one more step and appending one more
Cr point to the file at each step.  To the file 'fmin.val' fmin
Cr appends the current values of variables for which it wants the
Cr function evaluated in lines looking like
Cr     val: var1=val1 var2=val2 var3=val3 ...
Cr If the function is capable of generating the gradient, it
Cr should do so by appending a line of the form
Cr     gra: grad1 grad2 grad3 ...
Cr to the end of the file.  grad1, grad2, ... are ascii representations
Cr of the gradient of f at the point (val1,val2,...).
Cr If your function make does make the gradent, invoke fmin with
Cr switch -eps=0.  If the function is not capable of generating
Cr the gradient, it should just append the function value to
Cr the end of the file.
Cr Example (mode 1): minimize x*x+y*y+(x+1)*z wrt x,y for z=7
Cr   fmin x=1 y=2 "z=3 z+=4" x,y
Cr   dval -fg10 `tail -1 fmin.val` 'x*x+2*y*y+(x+1)*z' >> fmin.val
Cr   The above lines are repeated until convergence
Cr
Cr Example (mode 1 with gradients): minimize x*x+2*y*y+(x+1)*z wrt x,y
Cr   fmin -eps=0 x=1 y=2 "z=3 z+=4" x,y
Cr   dval -a 'gra: %;10d %;10d   %;10d' `tail -1 fmin.val
Cr   | sed s/val://` '2*x+z' '4*y' 'x*x+2*y*y+(x+1)*z' >> fmin.val
Cr   The above lines are repeated until convergence
Cr
Cr Mode 2 is similar to Mode 1, except that fmin calls the function
Cr by creating a subprocess, which must return information in the last
Cr line of the output.  As before, the function reads the current
Cr values of variables from the last line of fmin.val.
Cr Example (mode 2): minimize x*x+2*y*y+(x+1)*z wrt x,y for z=7
Cr   fmin -sh='dval -fg10 `tail -1 fmin.val` "x*x+2*y*y+(x+1)*z" >> fmin.val'\
Cr         x=1 y=2 "z=3 z+=4" x,y
Cr Example (mode 2, with gradients): minimize this fcn wrt x,y for z=7
Cr   fmin -sh='dval -a "gra: %;10d %;10d %;10d" `tail -1 fmin.val|sed s/val://`
Cr   "x*x/12+2*x+z" "4*y" "x*x*x/36+x*x+2*y*y+(x+1)*z" >>fmin.val'
Cr   -eps=0 x=1 y=2 "z=3 z+=4" x,y
Cr Example illustrating non-positive Hessian:
Cr   fmin -pr50 -sw=3100 -gtol=1d-5 -sh='dval -a "gra: %;10d %;10d %;10d"
Cr   `tail -1 fmin.val|sed s/val://` "x*x/12-2*x-z/2" "1*y"
Cr   "x*x*x/36-x*x+y*y/2-(x/2+1)*z" >>fmin.val' -eps=0 x=1 y=2 "z=3 z+=4" x,y
Cr
Cl Local variables
Cl   lint : flags -inline set
Cu Updates
Cu  28 Oct 99 (ATP) added Brent; renamed some command switches
C ----------------------------------------------------
      subroutine fmain
      implicit none
      integer iarg,i,j,n,ifi,ific,ifip,garg,i1mach,ipr,iprint
      logical cmdstr,lsequ,a2bin,vtabio,lbrent,lkorr,lcont,lint
C     logical swtmp
C ... Local variables specific to fmin
      integer ndim,chrln,nfcall,fopng
      parameter (ndim=100,chrln=4096)
      character*(chrln) outs,cnams(ndim),first,shcmd,fnam,fnamc,
     .                  fnami,fnamp
      double precision pos(ndim),wk(0:26),
     .  dxmn,xtol,gtol,atol,rtol,fxx,x0,x00,dxmx,delp,eps
      equivalence (x0,wk)
C ... For old style, need wk; for gradzr ,embedded into p
      double precision p(ndim*(2*ndim+11)),hess(ndim,ndim)
      equivalence (p,pos)
C     double precision wk(ndim,9+2*ndim),p(0:6*ndim)
      integer iter,isw,isw2,getdig,maxit,ir,iter0,nargf,nfcn,ich,ich0
C     integer broyj
      external fxx
      double precision hessin(ndim*ndim),hmax(ndim)
C#ifdefC OLD
C      double precision g(ndim),xnew(ndim)
C      double precision diff,wc,beta
C#endif

C ... Local variables specific to praxis and korr
      external zfn,restri,gleich,tkontr,gaussn
      double precision zfn1,h,btol,d(ndim),scbd,val,val0,
     .                 yy(ndim),zz(ndim),q0(ndim),q1(ndim),
     .                 tvec(ndim),v(ndim**2)
      double precision xkor(ndim),epsilo(4),deltas,deltai,deltap,zstern,
     .                 zbest,tgrenz,cnst
      double precision d1mach
      logical lq1,lq2,lc,ilsw
      integer ktm,ipm,iscale
      integer m,ns,np,ny,ielter,nachko,irekom,konvkr,ifallk,kanal,
     .        icnst,ncall
      logical bkomma,bkorrl
      integer oskor,opkor,oykor
      integer :: err
      data ktm /0/ ipm /0/ iscale /1/
     .     lq1,lq2 /.false.,.false./ lc /.false./
      data bkomma /.true./ irekom /3/ bkorrl /.true./ tgrenz /1000d0/
     .     konvkr /1/  deltas,deltai /1d0,1d0/ deltap /0.087d0/
      common /cnstr/ icnst(ndim), cnst(2,ndim), m
      common /calls/ ncall

      integer wksize
      parameter (wksize=1000000)
      integer w(wksize)
      common /w/ w

      call pshpr(0)
      call wkinit(wksize)
      call poppr

C --- Setup ---
      lcont = .false.
      lbrent = .false.
      lkorr = .false.
      lint = .false.
      isw2 = 0
      eps = 1d-4
      dxmx = 1d0
      dxmn = 1d-3
      xtol = 1d-2
      atol = 1d-2
      rtol = 1d-2
      gtol = 1d-2
      ncall = 1
C number of parents
      ielter = 5
C number of offspring
      nachko = 50
C number of constraints
      m = 0
      isw = 41
      maxit = 20
      shcmd = ' '
      fnam = 'fmin.val'
      fnamc = 'fmin.cnst'
      fnamp = ' '
      fnami = ' '
      call dpzero(hmax,ndim)
      call dpzero(hessin,ndim**2)
      err = -1

      goto 10
   20 print 333
  333 format(
     .  ' usage: fmin [switches] var1=expr var2=expr ...  var1,var2,...'
     . /'     or fmin -lint [switches]'
     .//'        var1,var2,... are varied to find zero gradient'/
     .  '        caller reads last line of file fmin.val for requested',
     .           ' values'/8x,'of input variables, and appends to',
     .           ' this file functions or gradients'/8x,
     .           'corresponding to those values.'//
     .  '        Switches:'/
     .  '        -file=# parameter-vals filename'/
     .  '        -ifile=# initial value filename'/
     .  '        -pfile=# parameter-list filename'/
     .  '        -cfile=# constraints filename'/
     .  '        -broyden use Broyden minimization'/
     .  '        -fp      use Fletcher-Powell minimization'/
     .  '        -brent   use Brent minimization '/
     .  '        -genetic use multimembered evolution strategy'/
     .  '        -eps=#   sets the width for finite differences',
     .           ' (numerical derivatives)'/
     .  '        -eps=0   flags fmin that function will return',
     .           ' gradients'/
     .  '        -eps=#   ''machine'' precision for Brent (see praxis)'/
     .  '        -pr#     verbosity')
      print 334
  334 format(
     .  '        -parents=#  number of parents (for genetic)'/
     .  '        -offs=#  number of offspring (ditto)'/
     .  '        -cons=#  number of constraints (ditto)'/
     .  '        -conv=#  convergence criterion (1 or >=2N) (ditto)'/
     .  '        -recom=# recombination type (1-5) (ditto)'/
     .  '        -plus    use plus not comma selection (ditto)'/
     .  '        -atol=#  absolute tolerance in x (ditto)'/
     .  '        -rtol=#  relative tolerance in x (ditto)'/
     .  '        -dels=#  deltas (ditto)'/
     .  '        -deli=#  deltai (ditto)'/
     .  '        -delp=#  deltap (ditto)')
      print 335
  335 format(
     .  '        -xtol=#  allowed tolerance in x (not used in genetic)'/
     .  '        -gtol=#  allowed tolerance in function'/
     .  '         (tol for Brent)'/
     .  '        -dxmx=#  maximum step to take in x (h for Brent)'/
     .  '        -dxmn=#  minimum step to take in x'/
     .  '        -sw=#    switches passed to gradzr'/
     .  '        -inline  fmin will use function call compiled',
     .          ' with this executable'/
     .  '        -sh=command : fmin will invoke in a subshell',
     .           ' ''command'''/17x,'to append f or grad f'/
     .  '        -vname=val create variable name and assign val'/
     .  '                 (for vars not associated with minimization)'/
     .  '        -lq1 -lq2 logical switches for Brent (see praxis)'/
     .  '        -ipm=# -ktm=# -scbd=# for Brent (see praxis)'/
     .  '        -cont    (with -sh only) resume using initial data',
     .           ' from fmin.val')
      call cexit(err,1)

   10 continue
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first) .and. lint) goto 30
C     no arguments: print usage
      if (.not. cmdstr(iarg,first)) then
        err = 0
        goto 20
      end if
C     finished all options (-..) goto pick off initial values
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (lsequ(first,'-pr',3,' ',n)) then
        i = 0
        if (.not. a2bin(first(4:40),j,2,0,' ',i,-1)) goto 20
        call pshpr(j)
      else if (first == '-broyden') then
        isw2 = 2
      else if (first == '-fp') then
        isw2 = 1
      else if (first == '-brent') then
        lbrent = .true.
      else if (first == '-genetic') then
        lkorr = .true.
      else if (first == '-lq1') then
      else if (first == '-inline') then
        lint = .true.
        nfcn = -1
      else if (first == '-lq2') then
        lq2 = .true.
      else if (first == '-plus') then
        bkomma = .false.
      else if (garg('-xtol=',iarg,4,' ',1,1,i,j,xtol) /= 0) then
      else if (garg('-atol=',iarg,4,' ',1,1,i,j,atol) /= 0) then
      else if (garg('-rtol=',iarg,4,' ',1,1,i,j,rtol) /= 0) then
      else if (garg('-gtol=',iarg,4,' ',1,1,i,j,gtol) /= 0) then
      else if (garg('-dxmx=',iarg,4,' ',1,1,i,j,dxmx) /= 0) then
      else if (garg('-dxmn=',iarg,4,' ',1,1,i,j,dxmn) /= 0) then
      else if (garg('-eps=',iarg,4,' ',1,1,i,j,eps) /= 0) then
      else if (garg('-ktm=',iarg,2,' ',1,1,i,j,ktm) /= 0) then
      else if (garg('-ipm=',iarg,2,' ',1,1,i,j,ipm) /= 0) then
      else if (garg('-scbd=',iarg,2,' ',1,1,i,j,iscale) /= 0) then
      else if (garg('-parents=',iarg,2,' ',1,1,i,j,ielter) /= 0) then
      else if (garg('-offs=',iarg,2,' ',1,1,i,j,nachko) /= 0) then
      else if (garg('-cons=',iarg,2,' ',1,1,i,j,m) /= 0) then
      else if (garg('-conv=',iarg,2,' ',1,1,i,j,konvkr) /= 0) then
      else if (garg('-recom=',iarg,2,' ',1,1,i,j,irekom) /= 0) then
      else if (garg('-dels=',iarg,4,' ',1,1,i,j,deltas) /= 0) then
      else if (garg('-deli=',iarg,4,' ',1,1,i,j,deltai) /= 0) then
      else if (garg('-delp=',iarg,4,' ',1,1,i,j,deltap) /= 0) then
      else if (first(1:4) == '-sh=') then
        shcmd = first(5:)
      else if (first(1:6) == '-file=') then
        fnam = first(7:)
      else if (first(1:7) == '-ifile=') then
        fnami = first(8:)
      else if (first(1:7) == '-pfile=') then
        fnamp = first(8:)
      else if (first(1:7) == '-cfile=') then
        fnamc = first(8:)
      else if (first(1:5) == '-cont') then
        lcont = .true.
C ... Variables declaration
      else if (first(1:2) == '-v') then
        first = first(1:len(first)-1) // ' '
        call wordg(first,11,'a-zA-Z0-9_*/+^-',1,i,j)
        if (first(j+1:j+1) == '=') then
          i = 2
          call parsyv(first,len(first),999,0,i)
        else
          goto 20
        endif
      else if (garg('-sw=',iarg,2,' ',1,1,i,j,isw) /= 0) then
      else if (lsequ(first,'-max',4,' ',n)) then
        call rx('-max not implemented')
      else if (first(1:6) == '--help' .or. first(1:2) == '-h') then
        err = 0
        go to 20
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15

C --- End of switches: Evaluate starting expressions ---
C     now first holds the first parameter=value pair
   30 continue
C     In-line case: call setup and skip reading variables
      if (lint) then
        call fcall0(nfcn,pos)
        if (nfcn <= 0)
     .    call fexit(-1,9,'fcall0 did not set number of variables',0)
        goto 40
      endif
      ifi = fopng(fnam,-1,0)
      call getpr(ipr)
      iarg = iarg+1
      j = 0
      call chrpos(first,'=',40,j)
C     BUG: following means expecting at least one more arg. Not if -cfile
      if (iarg < nargf()) then
        n = 0
        call parsyv(first,chrln,999,0,n)
        if (.not. cmdstr(iarg,first)) goto 20
        goto 30
      endif
      if (isw2 /= 0) isw = isw + 100 * (isw2-getdig(isw,2,10))

C --- Load table of variables to minimize ---
      ich  = 0
      ich0 = 0
      nfcn = 0
      if (fnamp /= ' ') then
        call awrit0(' fmin reading parameters from '//fnamp,' ',
     .              128,i1mach(2))
        ifip = fopng(fnamp,-1,0)
        read (ifip, '(a)') first
        call fclose(ifip)
      endif
   19 call chrps2(first,', ',2,chrln,ich,i)
      nfcn = nfcn+1
      call rxx(nfcn > ndim,'FMIN: increase ndim')
      cnams(nfcn) = first(ich0+1:ich)
      if (i == 1) then
        ich0 = ich+1
        ich = ich0
        if (ich < chrln) goto 19
      endif
      call vload(cnams,nfcn,pos,.true.)

C --- Starting printout ---
   40 if (ipr >= 20) then
        if (lbrent) then
          call awrit7(' fmin (Brent):  %i coffs.  h=%d tol=%1;3g '/
     .    /'ipm=%i ktm=%i scbd=%i eps=%1;3g',' ',80,i1mach(2),nfcn
     .    ,dxmx,gtol,ipm,ktm,iscale,eps)
        elseif (lkorr) then
          call awrit6(' fmin (genetic):  %i coffs. '/
     .    /'%i parents %i offspring %i constraints, min/max step:%d %d',
     .    ' ',256,i1mach(2),nfcn,ielter,nachko,m,dxmn,dxmx)
          call awrit5('                           '/
     .      /'%?;n;comma;plus; version, '/
     .      /'recom type %i, abs tol=%1;3g, rel tol=%1;3g, konvkr=%i'
     .      ,' ',256,i1mach(2),bkomma,irekom,atol,rtol,konvkr)
          call awrit3('                           '/
     .      /'deltas=%d deltai=%d deltap=%d',' ',256,i1mach(2),
     .      deltas,deltai,deltap)
        else
          call awrit7(' fmin:  %i coffs  sw=%i  dxmn=%1;3g  dxmx=%1;3g'/
     .      /'  xtol=%1;3g  gtol=%1;3g  eps=%1;3g',' ',80,i1mach(2),nfcn
     .      ,isw,dxmn,dxmx,xtol,gtol,eps)
        endif
      endif
      if (eps >= dxmn .and. .not. (lbrent .or. lkorr))
     .  print '(a)', ' fmin (warning): dmxn < eps'

C --- Iterate to a minimum ---
      if (lbrent .or. lkorr) then
        if (lint) call fexit(-1,9,'Brent/korr not set up for inline',0)
        if (shcmd == '')
     .        call fexit(-1,9,'Need a shell command for Brent/korr',0)
C --- First call to function (create static data semantics!) ---
        val0 = zfn1(ifi,fnam,cnams,nfcn,pos,shcmd,lcont,ir)
      endif
      if (lbrent) then
        h = dxmx
        btol = gtol
        ilsw = .false.
        scbd = iscale
C --- Minimise ---
        call praxis(zfn,nfcn,pos,btol,h,eps,scbd,ktm,ilsw,ipm,lq1,
     .    lq2,lc,v,d,yy,zz,q0,q1,tvec,val,ir)
        if (ir == -1) call fexit(-1,9,'praxis returned ir=-1',0)
        call fexit(0,0,0,0)
      elseif (lkorr) then
C constraints file
        if (m > 0) then
          ific = fopng(fnamc,-1,0)
          print *, ' constraints, min < parameter < max :'
          do  i = 1, m
            read (ific, *,err=1,end=1) icnst(i), (cnst(j,i), j = 1, 2)
            call awrit4(' %:16i: %d < %d < %d',' ',128,i1mach(2),
     .        icnst(i),cnst(1,i),pos(i),cnst(2,i))
          enddo
          goto 2
    1     continue
          call rx(' FMIN cannot read from constraints file '//fnamc)
    2     continue
        endif
C number of distinct step sizes
        ns = nfcn
        if (bkorrl) then
          np = nfcn*(ns - 1) - ( (ns - 1)*ns ) / 2
          ny = (nfcn + ns + np + 1)*ielter*2
        else
          np = 1
          ny = (nfcn + ns + 1)*ielter*2
        endif
        call defdr(oskor,ns)
        call defdr(opkor,np)
        call defdr(oykor,ny)
        call dcopy(ns,dxmx,0,w(oskor),1)
        call dcopy(np,0.5d0,0,w(opkor),1)
        epsilo(1) = dxmn
        epsilo(2) = 1
        epsilo(3) = atol
        epsilo(4) = rtol
        call ran1in(-1)
        kanal = i1mach(2)
        call korr(ielter,bkomma,nachko,irekom,bkorrl,konvkr,ifallk,
     .            tgrenz,epsilo,deltas,deltai,deltap,nfcn,m,ns,np,ny,
     .            zstern,pos,zbest,xkor,w(oskor),w(opkor),w(oykor),
     .            zfn,restri,gaussn,gleich,tkontr,kanal)
        if (ifallk == -2) then
          call awrit3(' Constraints not satisfied on startup.'//
     .      ' Best estimated w/o satisfying constraints:%N'//
     .      ' Objective function: %d, parameter set: %n:1d',' ',
     .      256,i1mach(2),zstern,nfcn,pos)
        elseif (ifallk == -1) then
          call awrit0(' Constraints not satisfied on startup.'//
     .      ' Ran out of time looking for feasible parameter set.',' ',
     .      256,i1mach(2))
        elseif (ifallk == 0) then
          call awrit2(' Constraints not satisfied on startup.'//
     .      ' Feasible parameter set found: %n:1d%N'//
     .      ' Recommend restart with this set.',' ',256,i1mach(2),
     .      nfcn,pos)
        elseif (ifallk == 1) then
          call awrit3(' Ran out of time. Currently,%N'//
     .      ' objective function: %d, parameter set: %n:1d',' ',
     .      256,i1mach(2),zstern,nfcn,pos)
        elseif (ifallk == 2) then
          call awrit4(' KORR terminated normally, %i function calls.'//
     .      ' Objective function: %d, Parameter set: %n:1d',' ',
     .      256,i1mach(2),ncall,zstern,nfcn,pos)
        elseif (ifallk == 3) then
          call awrit0(' Ran out of time. No estimate available',' ',128,
     .      i1mach(2))
        else
          call awrit1(' Exeptional undocumented exit from KORR,'//
     .      ' ifallk = %i',' ',256,i1mach(2),ifallk)
        endif
        call fexit(0,0,0,0)
      endif

C     open(ifi,file=fnam)
      rewind ifi
      nfcall = 0
      ir = 0
      iter = 0
      x0 = 0
      iter0 = iter
   50 continue
      nfcall = nfcall+1
      if (ir == 0) call dpcopy(pos,p(1),1,nfcn,1d0)
C ... Try and get the next point from vtabio of fcall
      if (lint) then
        call fcall(nfcn,eps,p,p(nfcn+1))
      else
        if (.not. vtabio(lcont,shcmd,2d-10,eps,
     .      cnams,p(nfcn+1),nfcn,ifi,fnam)) then
          call rx('fmin: something went wrong in vtabio')
        endif
C        outs = ' '
C        call awrit1(' fmin call %i seeks grad vars',outs,
C     .    chrln,0,nfcall)
C        do  52  j = 1, nfcn
C   52   call awrit0('%a '//cnams(j),outs,chrln,0)
C        call awrit0('%a. Var table now:',outs,-len(outs),-i1mach(2))
C        outs = ' val:'
C        if (eps /= 0) outs = ' '
C        call avtab(outs)
C        call awrit0(outs,' ',-len(outs),i1mach(2))
C        if (ir == 0) then
C          close(ifi,status='DELETE')
C          open(ifi,file='fmin.val')
C          rewind ifi
C        endif
C        swtmp = vtabio(lcont,shcmd,2d-10,eps,cnams,
C     .    p(nfcn+1),nfcn,-ifi,fnam)
C        close(ifi)
C        call cexit(0,1)
      endif

C ... Next step in minimization
      x00 = x0
      call dpdot(p(2*nfcn+1),p(2*nfcn+1),nfcn,delp)
      delp = x0*dsqrt(delp)
      iter0 = iter
      call pshpr(iprint()+10)
      call gradzr(nfcn,p,hess,dxmn,dxmx,xtol,gtol,1.2d0,wk,isw,ir)
      if (ir >= 0) then
        isw2  = mod(isw/100,10)
        call awrit3(' fmin:  '//
     .    '%?#n==0#converged#aborting (ir=%-1j%i)# '//
     .    '%?#n==0#C.G.##%-1j'//
     .    '%?#n==1#F.P.##%-1j'//
     .    '%?#n==2#Broyden# minimization# in %i function calls',
     .    ' ',80,i1mach(2),ir,isw2,nfcall)

        if (lint) call fcallx(nfcn,eps,p,p(nfcn+1))

        if (ir == 0) call cexit(0,1)
        call cexit(-1,1)
      endif
      call poppr

C#ifdefC OLD
C      if (lbroy) then
CC#ifndef BROYJ
C        call brmin(nfcn,pos,p(nfcn+1),isw,iprint(),dxmx,xtol,
C     .    gtol,hmax,wk,diff,hessin,ir)
C        if (ir == 0) call cexit(0,1)
CC#else
CC        ir = ir+1
CC        beta = 1d0
CC        wc = 1d0
CC        if (ir == 1) wc = .01d0
CC        call dscal(nfcn,-1d0,p(nfcn+1),1)
CC        j = broyj(nfcn,pos,p(nfcn+1),ir,isw,iprint(),beta,dxmx,
CC     .    xtol,gtol,wc,wk,xnew)
CC        if (j == 0) call cexit(0,1)
CC        call dscal(nfcn,-1d0,p(nfcn+1),1)
CC        call dcopy(nfcn,xnew,1,pos,1)
CC#endif
C      else
C        if (mod(isw,10)/4 == 0) isw = isw+4
C        call frpmin(nfcn,fxx,p,g,dxmn,dxmx,xtol,gtol,x0,maxit,
C     .    isw,iter,ir)
C        if (ir == 0 .or. ir == -10) call cexit(0,1)
C        do  62  j = 1, nfcn
C   62   pos(j) = p(j) + x0*p(2*nfcn+j)
C      endif
C#endif
      if (.not. lint) then
        call vload(cnams,nfcn,pos,.false.)
      endif
C#ifdefC OLD
CC ... Check for convergence or new line minimization
C      if (iter0 /= iter .and. iter /= 1 .and. ipr > 0) then
C        call awrit2(' fmin: cg line %i  delta p=%1;3g',
C     .    ' ',80,i1mach(2),iter0,delp)
C      else
CC        call awrit5(' atmve: cg line %i e=%1;6,6d  x0,ftop,shmx='//
CC     .    '%1;4,4d %1;4,4d %1;4,4d',' ',80,i1mach(2),iter,etot,
CC     .    x00,ftop,shmx)
CC      endif
C#endif
      goto 50

      end
      logical function vtabio(lcont,shcmd,tol,eps,cnams,grad,nfcn,
     .  ifi,fnam)
C- Read variable table's value from file ifi, return true if matches
C ----------------------------------------------------------------------
Ci Inputs
Ci   lcont :Applies to file read only (ifi>0)
Ci         :if T, tell rdline to read data using rdline mode 2,
Ci         :i.e. to search file ifi for a line that matches a key
Ci         :which is constructed out of an ascii representation of the
Ci         :symbolic variables to be varied.
Ci         :(See ifi below for descr. of how key is read or written.)
Ci         :In this mode, the line following the key is assumed to
Ci         :hold value of the function at the variables' value as
Ci         :specified by the key.
Ci         :vtabio continues builds up the gradient from these values
Ci         :until the entire gradient is made, in which case vtabio
Ci         :returns T.  Otherwise vtabio exits.
Ci         :if F, tell rdline to read data using rdline mode 1
Ci         :It is an alternative mode which vtabio can communicate
Ci         :with the external program to make function values
Ci         :and/or gradients (via shell command)
Ci         :This mode also requires that shcmd is specified.
Ci   shcmd :(optional) containing shell command vtabio may invoke to
Ci         :generate the data it seeks.  Not used if lcont is .false.
Ci   tol   :Not used now
Ci   eps   :<> 0: take derivatives numerically, by repeated calls
Ci         :to generate function values
Ci         :0 : call that generates function can return analytic
Ci         :derivatives
Ci   cnams :names of the independent variables
Ci   nfcn  :number of independent variables
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci         :In case of file write, vtabio writes a line
Ci         :containing an ascii representation of variables
Ci         :of the form
Ci         : val: vnam1=val1 vnam2=val2 ...
Ci   fnam  :file name passed to rdline.  Used only when it obtains
Ci         :its data via a shell command
Co Outputs
Co   vtabio:returns T when it finds all the data it needs.
Ci   grad  :gradient of function wrt variables
Cl Local variables
Cl   s0    :an ascii representation of the variables table,
Cl         :used as a key for rdline
Cl   s     :string read by rdline
Cr Remarks
Cr   Data is read from file using function rdline, in one of two modes,
Cr   which see
Cu Updates
Cu   26 Sep 02
C ----------------------------------------------------------------------
      implicit none
      logical lcont
      integer ifi,nfcn
      character*(*) cnams(nfcn),shcmd*4096,fnam*256
      double precision tol,grad(nfcn),eps,val,vp,vm
      integer j,k,i1mach,rdline,iprint,isw
      character*500 strn,outs*512,s0*512,s*512

      s0 = 'val:'
      if (eps /= 0) s0 = ' '
      call avtab(s0)

C --- Write ---
      if (ifi < 0) then
        call poseof(-ifi)
        call awrit0(s0,' ',-len(trim(s0)),-ifi)
        return
      endif

C --- Read ---
      vtabio = .false.
      k = rdline(ifi,fnam,lcont,shcmd,eps /= 0,s0,s,val)
      call skpblb(s0,len(s0),j)
      if (k == 0) goto 25
      if (iprint() >= 30) call awrit3(' vtabio: seek '//s0(1:j+1)//
     .  ' ... found=%i%?#n# val=%,8d##',outs,len(outs),i1mach(2),
     .  k,isw(k > 0.and.eps /= 0),val)
      if (k == -1) goto 21
      if (eps /= 0) then
        do  20  j =  1, nfcn
          s0 = ' '
          call avtab(s0)
          strn = cnams(j)
          call awrit0('%a '//strn,s0,len(s0),0)
          call awrit1('%a+=%1;5g',s0,len(s0),0,eps)
          k = rdline(ifi,fnam,lcont,shcmd,eps /= 0,s0,s,vp)
          if (k == 0) goto 25
          if (k == -1) goto 21
          s0 = ' '
          call avtab(s0)
          strn = cnams(j)
          call awrit0('%a '//strn,s0,len(s0),0)
          call awrit1('%a-=%1;5g',s0,len(s0),0,eps)
          k = rdline(ifi,fnam,lcont,shcmd,eps /= 0,s0,s,vm)
          if (k == 0) goto 25
          if (k == -1) goto 21
          grad(j) = (vp-vm)/(2*eps)
   20   continue
        if (iprint() >= 20) print
     .    '('' vtabio: grad= '',20(6f10.6:/15x))',(grad(j), j=1,nfcn)
        vtabio = .true.
        return
C ...   Exit when a value is missing
   24   continue
        strn = ' vtabio: seek '//s0
        call awrit0(strn,' ',-len(strn),i1mach(2))
        call poseof(ifi)
        call awrit0(s0,' ',-len(trim(s0)),ifi)
        close(ifi)
        if (shcmd == ' ') call cexit(1,1)
      else
C       read(ifi,'(a)',end=25,err=25) strn
        if (s(1:4) /= 'gra:') goto 25
        read(s(5:len(s)),*) grad
        vtabio = .true.
        return
      endif

   25 continue
      call skpblb(s0,len(s0),j)
      if (iprint() >= 30) call awrit1(' vtabio: seek '//s0(1:j+1)//
     .  ' ... found=%i',outs,len(outs),i1mach(2),k)
      call poseof(ifi)
      call awrit0(s0,' ',-len(trim(s0)),ifi)
      call fclose(ifi)
      call cexit(1,1)

C ... Abnormal exit for value read
   21 continue
      strn = ' vtabio: missing value for:'//s0
      call awrit0(strn,' ',-len(strn),i1mach(2))
      call rx('vtabio')
C ... Abnormal exit for variables read
   22 continue
      strn = 'vtabio: expected '//s0
      call awrit0(strn,' ',-len(strn),i1mach(2))
      strn = '       but found '//s
      call awrit0(strn,' ',-len(strn),i1mach(2))
      call rx('vtabio')

      end
      subroutine vload(cnams,n,pos,sw)
C- Load into position array variables
      implicit none
      logical sw
      integer n
      character*(*) cnams(n), s*40
      double precision pos(n)
      integer i,j
C ... get variables from table
      if (sw) then
        do  10  i = 1, n
          call getsyv(cnams(i),pos(i),j)
          if (j == 0) then
            s = cnams(i)
            call rx('vload: variables table is missing variable '//s)
          endif
   10   continue
C ... load variables into table
      else
        do  20  i = 1, n
          call chsyv(cnams(i),pos(i),j)
          if (j == 0) then
            s = cnams(i)
            call rx('vload: variables table is missing variable '//s)
          endif
   20   continue
      endif
      end

      double precision function fxx()
      fxx = 0d0
      end
      subroutine avtab(s)
C- builds up an ascii representation of the variables table
      implicit none
      character*(*) s
      character*40 cnam
      integer i,j,nvar,len
      double precision xx

      call numsyv(nvar)
      do  80  i = 4, nvar
        cnam = ' '
        call watsyv(cnam,xx,i)
        call skpblb(cnam,len(cnam),j)
        call awrit1('%a '//cnam(1:j+1)//'=%1;10g',s,len(s),0,xx)
C        call awrit1('%a '//cnam(1:j+1)//'=%1;5g',s,len(s),0,xx)
   80 continue
      end
      integer function rdline(ifi,fnam,lcont,shcmd,lcnvt,s0,s,val)
C- Reads from file ifi looking for a line containing data sought
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit
Ci   fnam  :file name for shell command (see Remarks, mode 1)
Ci   lcont :T tell rdline to return line following match to key s0
Ci         :(see Remarks, mode 2)
Ci   shcmd :shell command that should append desired data to
Ci         :end of file with name 'fnam'.
Ci         :Used ONLY in mode 1, described in Remarks.
Ci   lcnvt :T, convert s to a floating point val
Ci   s0    :key to match against file (see Remarks, mode 2)
Co Outputs
Co   rdline:  0 if value missing
Co         :  1 if data found
Co         : -1 if abnormal exit
Co   s     : ascii form of data returned by rdline
Co   val   :if lcnvt=T, s is also converted to a floating-point
Co         :number which is returned in val
Cl Local variables
Cl         :
Cr Remarks
Cr   rdline reads one line of data from logical unit ifi, and returns
Cr   the result in string s
Cr   It does this in one of two ways:
Cr     1. rdline passes a command to a subshell, whose job
Cr        it is to append to file ifi the data is requires.
Cr        shcmd contains a string that is passed to a subshell,
Cr        This method is invoked only if shcmd is nonblank and
Cr        lcont is .false.
Cr     2. it reads the next line from file ifi.  If this line matches
Cr        the key s0 passed, rdline reads the line following
Cr        into string s.  If not, rdline returns -1.
Cu Updates
Cu   26 Sep 02
C ----------------------------------------------------------------------
      implicit none
      logical lcnvt,a2bin,lcont
      character*(*) s0,s,shcmd*4096,outs*1055,fnam*256
      integer i1mach,j,ifi,iprint
      double precision val

C ... Re-entry
   10 continue
      if (shcmd == ' ' .or. lcont) then

        rdline = 0
        read(ifi,'(a)',end=25,err=25) s
        if (s0 /= s) goto 25
        rdline = -1
        read(ifi,'(a)',end=25,err=25) s
        rdline = 1

      else

        call poseof(ifi)
        call awrit0(s0,' ',-len(trim(s0)),ifi)
        close(ifi)
        call fsystm(shcmd,j)
        outs = ' rdline:"'//shcmd
        if (iprint() > 40 .or. iprint() >= 30 .and. j /= 0)
     .  call awrit1('%a returned %i',outs,len(outs),-i1mach(2),j)
        open(ifi,file=fnam)
        call flushs(1)
        call poseof(ifi)
        backspace ifi
        rdline = -1
        read(ifi,'(a)',end=25,err=25) s
        rdline = 1
      endif

      if (rdline <= 0) return
      if (lcnvt) then
        j = 0
        call skipbl(s,len(s),j)
        if (a2bin(s,val,4,0,' ',j,-1)) return
        rdline = -1
      endif

   25 continue
      if (lcont .and. shcmd /= ' ' .and. rdline <= 0) then
        lcont = .false.
        goto 10
      endif

      end

      double precision function zfn1(ifi0,fnam0,cnams0,nfcn0,x,shcmd0,
     .                                 lcont0,ir)
C- Function to be called from praxis invoked by fmin
      implicit none
C Passed variables
      integer ir,nfcn0,ifi0
      character*(*) shcmd0*4096,cnams0(nfcn0),fnam0*256
      logical lcont0
C Local variables
      logical lcont
      integer nfcn,ifi,j,i1mach,iprint,ndim,ncall
      parameter (ndim=100)
      double precision zfn,val,x(1)
      character*4096 cnams(ndim),shcmd
      character*4096 outs
      integer iend,slen
      character*256 fnam
      save ifi,nfcn,shcmd,cnams,lcont,fnam
      common /calls/ ncall

      call rxx(nfcn0 > ndim,'Increase ndim in zfn1')
      lcont = lcont0
      ifi = ifi0
      fnam = fnam0
      nfcn = nfcn0
      shcmd = shcmd0
      do  1  j = 1, nfcn
        cnams(j) = ' '
        cnams(j) = cnams0(j)
    1 continue
      ir = 0
      zfn1 = 0d0
      open(ifi,file=fnam)
      rewind ifi
      return

      entry zfn(x,nfcn0,ir)
    2 continue
      if (lcont) then
        read (ifi,'(a)',end=10,err=10) outs
        read (ifi,*,end=10,err=10) val
        zfn = val
        return
   10   continue
        lcont = .false.
        goto 2
      else
        call flushs(1)
        call togprt
        call vload(cnams,nfcn,x,.false.)
        outs = ' '
        call avtab(outs)
        open(ifi,file=fnam)
        call poseof(ifi)
        call awrit0(outs,' ',-len(trim(trim(outs))),ifi)
        close(ifi)
        outs = ' '
        call awrit0(' zfn:'//shcmd,outs,len(trim(outs)),0)
        call fsystm(shcmd,j)
        ncall = ncall + 1
        if (iprint() >= 60)
     .    call awrit1('%a returned %i',outs,len(outs),-i1mach(2),j)
        open(ifi,file=fnam)
        call poseof(ifi)
        backspace ifi
        read(ifi,*,end=20,err=20) val
        close(ifi)
        zfn = val
        ir = 0
        call togprt
        return
   20   continue
        ir = -1
        call togprt
      endif
      end

C#ifdef MULTIMF
      subroutine fcall0(n,val)
C- This is the initialization step to set up for later fcall
C  fcall0 returns n and starting values val
      implicit none
      integer n
      double precision val(*)

C --- multivariable mean-field problem ---
C     Here, the val are the z-projection of the magnetic spins m
C     INPUTS
C     *read Jij in units of K from FILE 'j0'
C      For n spins, this file should have n lines
C      Each line has n+1 points:
C        spin amplitude, followed by n values for Jij
C     *optionally read in file magdata.
C      Last line is used to guess magnetizations
C     *vars set through command-line arguments
C      -vtemp=# -vh=#
C
C     OUTPUTS
C    *fmin will append result to file 'magdata'
C     Each row in magdata is:
C       T(K)  eint(eV)  <mag>  mag1 mag2 ...
C
C     EXAMPLE
C     fmin -fp -sw=51 -xtol=1d-6 -dxmx=.01 -eps=0 -inline -vtemp=298
C
C ... You can extract j0 or jij from lmgf output:
C     mc -f12f8.2 -vkboltz=8.617d-5 dat -s13.6d0/kboltz -rsum -s2/3
C     mc -f12f8.2 -vkboltz=8.617d-5 dat -csum dat -ccat -s13.6d0/kboltz
C Local variables
      integer n0
      parameter (n0=50)
      integer i,j,ifi,fopng,a2vec,nw,rdm,j1,j2,ix(n0),m
      double precision beta,h,spina(n0),Jij(n0,n0),kboltz,T,wk(n0*n0+n0)
      double precision sumi
      character*512 a
C     parameter (kboltz=8.61724d-5)
      parameter (kboltz=1)
      common /fcallv/ beta,h,spina,Jij

      call dvset(val,1,n0,.999999d0)

      T = 100
      call getsyv('temp',T,i)
      if (i <= 0) call rx('fcall0: need specify temp using -vtemp=#')
      h = 0
      call getsyv('h',h,i)

C     Input values for Jij
      ifi = fopng('j0',-1,0)
      rewind ifi
      n = 0
      m = 0
      i = rdm(ifi,10,n0*n0,a,wk,n,m)
      call dmscop(spina,n0,wk,n,1,n,1,1,1,1,1d0)
      call dmscop(Jij,n0,wk,n,1,n,2,n+1,1,1,1d0)
      if (i < 0 .or. n < 0)
     .  call rx('fcall0: failed to read data from file j0')
      if (n+1 /= m)
     .  call rx('fcall0: for now # cols must = # rows + 1')
      call info(0,1,0,' fcall0: read %i element(s) from file j0',n,0)
      call info(0,1,0,'%5ftemp = %;10,2D K  h = %d',T,h)
      call info(0,0,0,'     |s|        Jij ...',0,0)
      do  i = 1, n
        write(*,356) spina(i), (Jij(i,j), j=1, n)
  356   format(f10.2,50f10.2)
      enddo

C     Initial guesses for magnetization
      ifi = fopng('magdata',-1,0)
      call poseof(ifi)
      backspace ifi
      read(ifi,'(a512)',end=90,err=90) a
      call words(a,nw)
      if (nw /= n+3) then
        call info(0,0,0,' fcall0: (warning) file magdata '//
     .    'should have %i columns but has %i',n+3,nw)
      else
        call word(a,4,j1,j2)
        j2 = 0
        i = a2vec(a(j1:),len(a)-j1,j2,4,', ',2,3,n,ix,val)
        if (i /= n) call rx(
     .    'fcall0: failed to read vals from last line of file magdata')
      endif
      goto 91
C     Re-entry point when nothing on last line of file magdata
   90 continue
      rewind  ifi
      write(ifi,'(''% cols'',i4)') n+3
C     Set sign of initial values of m according to J
      do  i = 2, n
        sumi = 0
        do  j = 1, n
          sumi = sumi + Jij(i,j)
        enddo
        if (sumi < 0) val(i) = -val(i)
      enddo
C     Re-entry point when existing data in file magdata
   91 continue
      call fclose(ifi)
      call info(0,0,0,' starting m = %n:1;10,6D ',n,val)

C     Input J's are in units of temperature; scale to energy
      beta = 1/kboltz/T
      do  i = 1, n
        do  j = 1, n
          Jij(i,j) = Jij(i,j)*kboltz
        enddo
      enddo

      return
   99 continue
      call rx('no temperature in file magdata')
      end

      subroutine fcall(n,eps,val,grad)
C- This is the in-line function call.  Must return gradient
      implicit none
      integer n
      double precision eps
      double precision val(n),grad(n)
C --- multivariable mean-field problem ---
C Local variables
      integer n0,i,j,k
      double precision atanzm,kboltz
      parameter (n0=50)
C     parameter (kboltz=8.61724d-5)
      parameter (kboltz=1)

      double precision si(n),w(n),gfi,z(0:n),dwik,sumfi2,fi,swi
      double precision beta,h,spina(n0),Jij(n0,n0)
      common /fcallv/ beta,h,spina,Jij

C     recursive call to get all 2^n states for partition f
C     maybe later?
C      z = 0
C      call fcall2(n,1,0,si,z)

C     Weiss field at site i = sum_j J_ij m_j
C     w(i) = w_i = beta (h + sum J_ij m_j)
C     NB: m_j = val(j)
      do  i = 1, n
        w(i) = 0
        do  j = 1, n
          w(i) = w(i) + (Jij(i,j)+Jij(j,i))*val(j)
        enddo
        w(i) = beta*(h + w(i))
      enddo

C     Solve simultaneous equations for fi
C     As a minimization of sum fi**2
C     f = sum_i f_i^2;  f_i = m_i + 1/w_i - s_i/tanh(s_i w_i)
C     df/dm_k   = sum_i 2 * f_i * df_i/dm_k
C     Use d(1/tanh(x))/dx = -1/sinh^2(x);
C     df_i/dm_k = delta_ik - w'_ik/w_i^2 + s_i^2/sinh^2(s_i w_i) * w'_ik
C     w'_ik = dw_i / dm_k = beta * J_ik
      call dpzero(grad,n)
      sumfi2 = 0
      do  i = 1, n
        swi = spina(i)*w(i)
        fi = val(i) + 1/w(i) - spina(i)/tanh(swi)
        sumfi2 = sumfi2 + fi**2
C       diagonal part
        grad(i) = grad(i) + 2 * fi
        do  k = 1, n
          dwik = beta*(Jij(i,k)+Jij(k,i))
          grad(k) = grad(k) +
     .              2 * fi * dwik * ((spina(i)/sinh(swi))**2-1/w(i)**2)
        enddo
C       print *, 'fi=',fi
      enddo
C     print *, 'sumfi2=',sumfi2

C      i = 2
C      k = 1
C      swi = spina(i)*w(i)
C      fi = 1
C      if (i /= k) fi = 0
C      dwik = beta*(Jij(i,k)+Jij(k,i))
C      print *, 'df/dm=',fi + dwik * ((spina(i)/sinh(swi))**2-1/w(i)**2)
C      print *, dwik, (spina(i)/sinh(swi))**2, 1/w(i)**2

C     print *, val,grad
C      call info(0,0,0,' val = %n:1;6g ',n,val)
C      call info(0,0,0,' grad = %n:1;6g ',n,grad)

      end

C     recursive call to get all 2^n states for partition f
C     use later?
C      subroutine fcall2(n,k,kg,si,z)
C      implicit none
C      integer n,k,kg
C      double precision si(n),z(0:n)
C      integer i,pm
C
C      do  pm = 1, -1, -2
C        si(k) = pm
C        if (k < n) then
C          call fcall2(n,k+1,kg,si,z)
C        else
C          print 333, (int(si(i)), i=1, n)
C  333     format(20i2)
C        endif
C      enddo
C      end

      subroutine fcallx(n,eps,val,grad)
C- This is the cleanup step after fmin quits
      implicit none
      integer n
      double precision eps
      double precision val(n),grad(n)
C Local variables
      integer n0,i,j,ifi,fopng
      double precision atanzm,kboltz,fi,gfi,T,eint
      parameter (n0=50)
C     parameter (kboltz=8.61724d-5)
      parameter (kboltz=1)

      double precision beta,h,spina(n0),Jij(n0,n0),x(n0)
      common /fcallv/ beta,h,spina,Jij

      eint = 0
      do  i = 1, n
      do  j = 1, n
        eint = eint + val(i)*Jij(i,j)*val(j)
      enddo
      enddo
      eint = eint * 8.61724d-5

      T = 1/kboltz/beta
      write(*,333) (val(i), i=1,n)
      write(*,332) T, sum(val(1:n))/n, eint
  332 format(' Average Magnetization for T=',f12.6,':',f12.6,
     .  ';  eint =',f12.6,' eV')
  333 format(60f12.6)

      ifi = fopng('magdata',-1,0)
      call poseof(ifi)
      write(ifi,333) T, eint, sum(val(1:n))/n, (val(i), i=1,n)
      call fclose(ifi)

      end
C#elseifC SIMPLMF
C      subroutine fcall0(n,val)
CC- This is the initialization step to set up for later fcall
CC  fcall0 returns n and starting values val
C      implicit none
C      integer n
C      double precision val(*)
C
CC --- simple mean-field problem ---
CC     read 2/3 j0(1..n) in units of K from FILE 'j0'
CC     fmin will OUTPUT result in file 'magdata'
CC Local variables
C      integer n0
C      parameter (n0=50)
C      integer i,ifi,fopng,a2vec,nw,rdm,j1,j2,ix(n0)
C      double precision beta,h,J0(n0),kboltz,T
C      character*120 a
C      parameter (kboltz=8.61724d-5)
C      common /fcallv/ beta,h,J0
C
C      call dvset(val,1,n0,.999999d0)
C
C      call getsyv('temp',T,i)
C      h = 0
C      call getsyv('h',h,i)
C
C      if (i <= 0) call rx('fcall0: need specify temp using -vtemp=#')
C
CC     Input values for J0
C      ifi = fopng('j0',-1,0)
C      rewind ifi
C      n = 0
C      i = rdm(ifi,10,n0,a,J0,n,1)
C      if (i < 0 .or. n < 0)
C     .  call rx('fcall0: failed to read data from file j0')
C      call info(0,1,0,' fcall0: read %i element(s) from file j0',n,0)
C      call info(0,1,0,'%7ftemp = %;10,2D K  h = %d',T,h)
C      call info(0,0,0,'%5f2/3 j0 = %n:1;10,2D K',n,J0)
C
CC     Initial guesses for magnetization
C      ifi = fopng('magdata',-1,0)
C      call poseof(ifi)
C      backspace ifi
C      read(ifi,'(a120)',end=90,err=90) a
C      call words(a,nw)
C      if (nw /= n+2) then
C        call info(0,0,0,' fcall0: (warning) file magdata '//
C     .    'should have %i columns but has %i',n+2,nw)
C      else
C        call word(a,3,j1,j2)
C        j2 = 0
C        i = a2vec(a(j1:),len(a)-j1,j2,4,', ',2,3,n,ix,val)
C        if (i /= n) call rx(
C     .    'fcall0: failed to read vals from last line of file magdata')
C      endif
CC     re-entry point when no data in file magdata
C   90 continue
C      call fclose(ifi)
C
C      call info(0,0,0,' starting m = %n:1;10,6D ',n,val)
C
C      beta = 1/kboltz/T
C      n = 1
C      do  i = 1, n
C        J0(i) = 4d0/3*1.5d0*J0(i)*kboltz
C      enddo
C
C      return
C   99 continue
C      call rx('no temperature in file magdata')
C      end
C
C      subroutine fcall(n,eps,val,grad)
CC- This is the in-line function call.  Must return gradient
C      implicit none
C      integer n
C      double precision eps
C      double precision val(n),grad(n)
CC --- simple mean-field problem ---
CC Local variables
C      integer n0,i
C      double precision atanzm,kboltz
C      parameter (n0=50,kboltz=8.61724d-5)
C      double precision beta,h,J0(n0),x(n0),fi,gfi
C      common /fcallv/ beta,h,J0
C
C
C      do  i = 1, n
C        x(i) = beta*(h+J0(i)*val(i)/2)
C        fi = val(i)/2-0.5d0*tanh(x(i))
C        gfi = 0.5d0-0.25d0*beta*J0(i)/cosh(x(i))**2
C        grad(i) = 2*fi*gfi
C      enddo
C
C      end
C
C      subroutine fcallx(n,eps,val,grad)
CC- This is the cleanup step after fmin quits
C      implicit none
C      integer n
C      double precision eps
C      double precision val(n),grad(n)
CC Local variables
C      integer n0,i,ifi,fopng
C      double precision atanzm,kboltz
C      parameter (n0=50,kboltz=8.61724d-5)
C      double precision beta,h,J0(n0),x(n0),fi,gfi,T
C      common /fcallv/ beta,h,J0
C
C      T = 1/kboltz/beta
C      write(*,333) (val(i), i=1,n)
C      write(*,332) T, sum(val(1:n))/n
C  332 format(' Average Magnetization for T=',f12.6,':',f12.6)
C  333 format(20f12.6)
C
C      ifi = fopng('magdata',-1,0)
C      call poseof(ifi)
C      write(ifi,333) T, sum(val(1:n))/n, (val(i), i=1,n)
C      call fclose(ifi)
C
C      end
C#endif
      double precision function tkontr(d)
      implicit none
      double precision d, cpusec
      d = 0d0
C#ifdefC TIME_LIMIT
C      tkontr = cpusec()
C#else
      tkontr = 0d0
C#endif
      end

      double precision function gleich(d)
      implicit none
      real*8 d
      real*4 gleich4
      real ran1
      gleich4 = ran1()
      gleich = real(gleich4,8)
      end

      double precision function restri(j,n,x)
C- Constraints for the genetic algorithm
C ----------------------------------------------------------------------
Ci Inputs: j: constraint number
Ci         n: number of parameters
Ci         x: vector of parameters
Ci
Co Outputs:
Co   restri returns 0 if constraint j is satisfied, -1 otherwise
Cr Remarks
Cr  restri is the user provided function to korr.f (q.v.) which enquires
Cr  whether the parameters in vector x satisfy the j^th constraint.
Cr
Cr  As currently implemented in fmin, the only constraints allowed are
Cr  simple bounds on each parameter x(n). These are read by fmin from
Cr  the file whose name is entered at the -cfile=# command line argument
Cr  and whose default name is "fmin.cnst" These constraints are passed
Cr  into restri through the COMMON block cnst: m is the number of
Cr  constraints, icnst(j) labels the parameter bounded by the j^th
Cr  constraint and cnst(2,j) are the lower and upper bounds.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer j,n
      double precision x(n)
C Local Variables
      integer ndim, iprint, i1mach
      parameter (ndim=100)
      integer m,i
      integer icnst
      double precision cnst
      common /cnstr/ icnst(ndim), cnst(2,ndim), m
      restri = 0d0
      if (x(icnst(j)) < cnst(1,j)) restri = -1d0
      if (x(icnst(j)) > cnst(2,j)) restri = -1d0
      if (iprint() > 40 .or. (iprint() >= 30 .and. restri < 0))
     .  then
        call awrit5(' RESTRI: constraint %i: %d < %d < %d, restri=%d',
     .    ' ',256,i1mach(2),j,cnst(1,j),x(icnst(j)),cnst(2,j),restri)
      endif
      end
