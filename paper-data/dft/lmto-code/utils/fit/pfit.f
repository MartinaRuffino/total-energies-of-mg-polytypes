
      module pfit_const
      implicit none
      integer, parameter :: mxdim=128, mxbas=1024, mxdat=1000
      end module pfit_const

      module pfit_static
      use pfit_const
      implicit none
      integer :: ndim,nbas,ptab(mxdim*mxbas)
      real(8) :: coffs(mxdim)
      end module pfit_static

      module pfit_static2
      use pfit_const
      implicit none
      real(8) :: dat((mxbas+1)*mxdat), datt((mxbas+1)*mxdat),a(mxdim*mxdat),
     & res(mxdat),r(mxdim**2),q(mxdim*mxdat),wp(mxdat),xp(mxdat),yp(mxdat)
      end module pfit_static2




C --- Least-squares fit to a polynomial of several variables ---
C*pfit uses rdm to reads in an array from a disk file.
C Normally pfit assumes the last column is for data, and all
C but the last are the independent variables.
C ----------------------------------------------------------------
C#define unix
C#ifdef unix
      subroutine fmain
C#endif
      use pfit_static
      use pfit_static2

      implicit none
      logical fail(3)
      double precision
     .  xxv(50),
     .  c(4*(mxdat+mxdim)+2),resmx,vmax,ddot
      character*20 first,fmt
      character(len=80) :: pstrn = ' '
      character*1 f2(20)
      logical cmdstr,lsequ,a2bin
      EXTERNAL cmdstr
      procedure(integer) :: nargf
      integer iarg,n,n1r,n1c,ifi,i1mach,i,garg,fopng,nrt,nct,rdops,rdm
      integer pwr(mxbas),sump(mxdim),npoly,
     .  n1,ipivot(mxdim),ifault,nx,it(50)
C Mapping of abscissa and ordinate
      character*72 mlist
      integer ix(2),ssvdef,np
      double precision wk(mxdat)
      parameter (ssvdef=256)
      character*(ssvdef) sxdef,sydef,outs,outs2
C For tabulating fit on a mesh ...
      integer nmsh
      double precision xmin(mxbas),xmax(mxbas),dx(mxbas)
      logical ltab,lw
      equivalence (f2,first)


C#ifdef DEBUG
      double precision snorm,ss,z
      integer m1,nw,j,k,l,m1p1,n1p1,k2,m
      data m1 /0/, l /1/
C#endif

      n1r = 0
      n1c = 0
      rdops = 0
      mlist = ' '
      ix(1) = 1
      sxdef = ' '
      sydef = ' '
      ltab = .false.
      nrt = 0
      nct = 0
      lw = .false.
      ifi = 10
      goto 10
   20 print 336
      print 377
  336 format(
     .  'usage: pfit [-switches] npoly [file]'/
     .  '       npoly = polynomial order'/
     .  '       file  = file name, holding n+1+nw columns of data:'/
     .  '               or n+2 columns if weights used (see -w below)'/
     .  '               columns 1..n are the independent variables'/
     .  '               column  1+n is the dependent variable'/
     .  '               column  2+n is the fitting weight, if used.')
  377 format(
     . /'       Switches:'//
     .  '       -pp1p2...pn'/
     .  '              individually specify polynomial order of',
     .                 ' each independent variable'//
     .  '       -tp [nc]~xlist .. '/
     .  '              use specified table of points instead of file'//
     .  '       -t[n] x1 x2 dx'/
     .  '              tabulate fit on a mesh from x1..x2, spacing dx'//
     .  '       -nc=#  specify number of columns the file has'//
     .  '       -w     fit with weighting factor'//
     .  '       -g     use fortran "g" print format'//
     .  '       -pr    printout verbosity'//
     .  '       -ord expr  use "expr" for independent variable')
      call cexit(0,1)

   10 continue
      fmt = '(5f12.6)'
      iarg = 1
   15 continue
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (lsequ(first,'-pr',3,' ',n)) then
        i = 0
        if (.not. a2bin(f2(4),j,2,0,' ',i,-1)) goto 20
        call pshprt(j)
      else if (lsequ(first,'-w',2,' ',n)) then
        lw = .true.
      else if (lsequ(first,'-g',2,' ',n)) then
        fmt = '(5g14.6)'
      else if (first(1:4) == '-tp ') then
        if (garg(' ',iarg+1,4,',:~ ',4,50,i,it,xxv) < 1) goto 20
        iarg = iarg+1
        call dpzero(datt,(mxbas+1)*mxdat)
        call expand(i,it,xxv,datt,nrt,nct)
      else if (lsequ(first,'-t',2,' ',n)) then
        nmsh = 1
        if (f2(3) /= ' ') then
          i = 0
          if (.not. a2bin(f2(3),nmsh,2,0,' ',i,-1)) goto 20
        endif
        do  17  j = 1, nmsh
          i = 0
          if (.not. cmdstr(iarg+1,first) .or.
     .        .not. a2bin(first,xmin(j),4,0,' ',i,-1)) goto 20
          iarg = iarg+1
          i = 0
          if (.not. cmdstr(iarg+1,first) .or.
     .        .not. a2bin(first,xmax(j),4,0,' ',i,-1)) goto 20
          iarg = iarg+1
          i = 0
          if (.not. cmdstr(iarg+1,first) .or.
     .        .not. a2bin(first,dx(j),4,0,' ',i,-1)) goto 20
          iarg = iarg+1
   17   continue
        ltab = .true.
      else if (lsequ(first,'-map',4,' ',n)) then
        if (.not. cmdstr(iarg+1,mlist)) goto 20
        iarg = iarg+1
      else if (lsequ(first,'-ord',4,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,sydef)) goto 20
      else if (lsequ(first,'-p',2,' ',n)) then
        pstrn = first
      elseif (garg('-nc=',iarg,2,' ',1,1,i,it,n1c) /= 0) then
        if (i < 0) goto 20
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15
   30 continue

      read(first,'(i4)',err=20) npoly
      iarg=iarg+1

      if (.not. cmdstr(iarg,first) .and. nargf() < iarg) goto 20
      if (first == '.' .or. nargf() == iarg) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif
      i = rdm(ifi,rdops,mxbas*mxdat,' ',dat,n1r,n1c)
      if (i < 0)
     .  call fexit(-1,9,'pfit failed to parse file '//first,0)
      if (n1r > mxdat)
     .  call fexit(-1,1,'pfit: input data exceeds %i rows',mxdat)
      ix(2) = n1c
      if (lw) ix(2) = n1c-1
      if (n1c*n1r > mxbas*mxdat) stop 'pfit: nput array too big'
      if (ifi /= i1mach(1)) call fclose(ifi)
C ... done with input

C --- Extract maximum powers from pstrn ---
      nx = n1c-1
      if (lw) nx = nx-1
      if (nx <= 0) call rx('no independent variables')
      do  40  i = 1, nx
        pwr(i) = npoly
        if (pstrn(2+i:2+i) /= ' ') read(pstrn(2+i:2+i),'(i1)') pwr(i)
        if (pwr(i) > npoly) pwr(i) = npoly
   40 continue

C --- Generate powers table ---
      call mkptab(pwr,nx,npoly,mxdim,ptab,sump,ndim)

C --- Map the data (a used as a work array)
      np = n1r
      call mapdat(dat,np,n1c,ix,ix(2),sxdef,sydef,mlist,wk,xp,yp)

C --- Weights (last column for now) ---
      call dcopy(np,1d0,0,wp,1)
      if (lw) then
        call dcopy(np,dat(1+(n1c-1)*n1r),1,wp,1)
      endif

C --- Generate normal matrix ---
      call nrmgen(dat,n1r,np,nx,ptab,a,ndim)

C --- Least-squares fit ---
      call l2a(np,ndim,0,1,a,yp,wp,0d0,np,ndim,
     . n1, ipivot, coffs, res, r, q, c, ifault)

C --- Tabulate fit if switch set (comment out unless recursive calls)---
      if (ltab) call tfit(nmsh,xmin,xmax,dx,ndim,nx,npoly,ptab,coffs)

C --- Tabulate fit if switch set ---
      if (nrt /= 0 .and. nct /= 0) then
        call nrmgen(datt,nrt,nrt,nx,ptab,a,ndim)
        call awrit2('%% rows %i cols %i',' ',80,i1mach(2),nrt,nct+1)
        do  50  i = 1, nrt
          resmx = ddot(ndim,a(i),nrt,coffs,1)
          print 337, (datt(i+k*nrt),k=0,nct-1), resmx
  337     format(10f16.8)
   50   continue
        call cexit(0,1)
      endif

C --- Print computed results ---
      nw = i1mach(2)
      m = np
      nbas = nx
      l = 1
  210 continue
      call awrit3(' ndim = %i  rank = %i  fault = %i',' ',
     .  80,i1mach(2),ndim,n1,ifault)
      if (ifault /= 0) call fexit(-1,009,'fault encountered',0)
C      write (nw,99890) ifault, ndim, n1
      if (ifault >= 1 .and. ifault <= 4) go to 390
      if (n1 == 0) go to 390
      if (n1 < m1) go to 390
C      write (nw,99870)
C      write (nw,99860) (ipivot(j),j=1,n1)
      if (n1 == ndim) go to 220
      n1p1 = n1 + 1
      write (nw,99850)
      write (nw,99860) (ipivot(j),j=n1p1,ndim)
  220 continue
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
c
c compute sum of squared residuals, norm of residuals, residual
c standard deviation, standard deviations of coefficients, and
c predicted values.  print these quantities, together with
c coefficients, observed values and residuals.
c
        if (fail(k)) go to 350
        ss = 0.0d0
        if (m == m1) go to 310
        m1p1 = m1 + 1
        do 250 i=m1p1,m
          ss = ss + (res(i)*wp(i))**2
  250   continue

  310   write (nw,99790)
        do 320 j=1,ndim
          write (nw,99775) j,coffs(j), (ptab(i+(j-1)*nbas), i=1,nbas)
  320   continue

        if (nbas == 1) then
          outs = '((((((((((((((((((((((('
          outs2 = outs
          outs(ndim:) = '0 '
          outs2(ndim-1:) = '0 '
          do  j  = ndim, 1, -1
            if (coffs(j) > 0) call awrit0('%a+',outs,len(outs),0)
            if (coffs(j) > 0 .and. j > 1)
     .        call awrit0('%a+',outs2,len(outs2),0)
            if (j > 1) then
              call awrit1('%a%;15e)*x',outs,len(outs),0,coffs(j))
              if (j > 2) then
                call awrit1('%a%;15e)*x',outs2,len(outs2),0,
     .            (j-1)*coffs(j))
              else
                call awrit1('%a%;15e',outs2,len(outs2),0,coffs(j))
              endif
            else
              call awrit1('%a%;15e',outs,len(outs),0,coffs(j))
            endif
          enddo
          write (nw,987) 'polynomial'
  987     format(' Evaluate ',a,' with this string:')
          call awrit0('%a',outs,len(outs),-nw)
          write (nw,987) 'derivative'
          call awrit0('%a',outs2,len(outs2),-nw)
        endif

  330   write (nw,99780)
        resmx = 0
        vmax = 0
        do 340 i=1,m
          vmax = max(vmax,dabs(yp(i)))
          if (dabs(res(i)) > dabs(resmx)) resmx = res(i)
          z = yp(i) - res(i)
          forall (j = 1:nbas) xxv(j) = dat(i+(j-1)*n1r)
          call info8(1,0,0,'%,3i   %;9F%7f%;9F    %;10F %n:3;8F',i,yp(i),z,res(i),nbas,xxv,7,8)
C          write (nw,99770) i,yp(i),z,res(i),(dat(i+(j-1)*n1r), j=1,nbas)
  340   continue

        snorm = dsqrt(ss/m)
        call awrit3(' RMS err = %1,3;3g  Max err = %1,3;3g  Max err'//
     .   '/Max val = %1,3;3g',first,80,nw,snorm,resmx,resmx/vmax)

C       write (nw,99760) ss,snorm,resmx
C       if ((mode == 2 .and. n1 < n) .or. (m == n1+kz)) go to 350
C       write (nw,99750) sd
  350 continue
c
c print lower triangular portion of symmetric unscaled covariance
c matrix.
c
c      if (mode == 2 .and. n1 < n) go to 390
c      write (nw,99740)
c      if (mode == 2) go to 370
c      do 360 i=1,ndim
c        write (nw,99910) (r(i,j),j=1,i)
c  360 continue
c      go to 390
c  370 do 380 i=1,ndim
c        write (nw,99910) (qr(i,j),j=1,i)
c  380 continue
c  390 read (nr,99730) ifdone
c
c
  390 continue
      call evfit(ndim,nbas,npoly,pwr,ptab,coffs)


c format statements.
99910 format (1x,8g15.8)
99890 format (' fault=',i2,'  ndim=',i2,'  computed rank=',i3)
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
99790 format (/5x,'j',2x,'Coefficient(j)',5x,'Powers ...')
99780 format (/2x,
     .  'i   Observed        Predicted      Residual    x ...')
99775 format (i6,g17.8,3x,20i3)
C99770 format (i3,2g16.8,g10.4,2x,5g12.6)
99760 format (30h sum of squared residuals    =,g15.8/1x,8hnorm of ,
     * 9hresiduals,11x,1h=,g15.8/1x,13hmax residual=,16x,g15.8)
99750 format (30h residual standard deviation =,g15.8)
99740 format (27h unscaled covariance matrix/)
99730 format (i5)
99720 format(25h number of zero weights =,i3,5x,17hdeg. of freedom =,i3)

      end

      subroutine nrmgen(dat,n1r,ndat,nbas,ptab,norm,ndim)
      use pfit_const
C- Generate normal matrix for polynomial fit
      implicit none
      integer ndat,nbas,n1r,ndim,ptab(nbas,ndim)
      double precision dat(n1r,1), norm(ndat,ndim)
      double precision x(mxbas),evterm
      integer i,idim

      if (nbas > mxbas) call rx('nrmgen: nbas gt mxbas')

C --- for each data set do ... ---
      do  20  i = 1, ndat
        call dcopy(nbas,dat(i,1),n1r,x,1)
        do  25  idim = 1, ndim
   25   norm(i,idim) = evterm(nbas,ptab(:,idim),x)
   20 continue
      end
      double precision function evpoly(ndim,nbas,ptab,coffs,x)
C- Evaluates a polynomial of several variables at a point
      implicit none
      integer ndim,nbas,ptab(nbas,ndim)
      double precision x(1),coffs(1)
      integer idim
      double precision evterm

      evpoly = 0
      do  25  idim = 1, ndim
C       print *, x, coffs(idim), evterm(nbas,ptab(1,idim),x)
        evpoly = evpoly + coffs(idim)*evterm(nbas,ptab(1,idim),x)
   25 continue

      end

      double precision function evterm(nbas,ptab,x)
C- Evaluates term in a polynomial of several variables at a point
      implicit none
      integer nbas,ptab(nbas),i
      double precision x(nbas)

      evterm = 1
      do  10  i = 1, nbas
   10 if (ptab(i) /= 0) evterm = evterm * x(i)**ptab(i)
      end
      subroutine gradp(ndim,nbas,ptab,coffs,x,grad)
C- Evaluates gradient of a multivariate polynomial
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndim  :number of terms
Ci   nbas  :number of independent variables
Ci   ptab  :table of polynomial orders for each of ndim terms
Ci   coffs :coefficients to each of 1..ndim terms
Ci   x     :independent variables at which to evaluate derivative
Co Outputs
Co   grad  :first derivative of polynomial at x
C ----------------------------------------------------------------------
      implicit none
      integer nbas,ndim,ptab(nbas,ndim)
      double precision x(1),coffs(1),grad(nbas)
      double precision dvpoly
      integer igrad(nbas)
      integer i

      do  20  i = 1, nbas
        call iinit(igrad,nbas)
        igrad(i) = 1
        grad(i) = dvpoly(ndim,nbas,ptab,coffs,x,igrad)
   20 continue
C     call awrit2('grad: %n:1d',outs,80,-6,nbas,grad)

      end
      double precision function dvpoly(ndim,nbas,ptab,coffs,x,igrad)
C- Evaluates some derivative of a multivariate polynomial
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndim  :number of terms
Ci   nbas  :number of independent variables
Ci   ptab  :table of polynomial orders for each of ndim terms
Ci   coffs :coefficients to each of 1..ndim terms
Ci   x     :independent variables at which to evaluate derivative
Ci   igrad :order of derivative which to calculate:
Ci         :there is an igrad for each 1..nbas
ci         :igrad(i) = 0 for zeroth-derivative,
ci         :         = 1 for 1st derivative etc
Co Outputs
Co  dvpoly :igrad^th derivative of the polynomial
C ----------------------------------------------------------------------
      implicit none
      integer nbas,ndim,ptab(nbas,ndim),igrad(1)
      double precision x(1),coffs(1)
      integer idim
      double precision dvterm

      dvpoly = 0
      do  25  idim = 1, ndim
c        print *, coffs(idim), dvterm(nbas,ptab(1,idim),x,igrad)
        dvpoly = dvpoly + coffs(idim)*dvterm(nbas,ptab(1,idim),x,igrad)
   25 continue

      end
      double precision function dvterm(nbas,ptab,x,igrad)
C- Evaluates derivative of one term in a multivariate polynomial
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :number of independent variables
Ci   ptab  :table of polynomial orders for this term
Ci   x     :point at which to evaluate derivative
Ci   igrad :order of derivative which to calculate
Co Outputs
Co  dvterm :igrad^th derivative of this term, igrad for each 1..nbas
C ----------------------------------------------------------------------
      implicit none
      integer nbas,ptab(1),igrad(1)
      double precision x(3)
      double precision dfact
      integer i,j

      dvterm = 1
      do  10  i = 1, nbas
        j = ptab(i)-igrad(i)
        if (j < 0) then
          dvterm = 0
          return
        endif
        dvterm = dvterm * dfact(ptab(i)) / dfact(j)
        if (j > 0) dvterm = dvterm * x(i)**j
   10 continue
      end
      double precision function dfact(n)
C- Factorial of n
      implicit none
      integer n
      integer i

      dfact = max(n,1)
      do  10  i = n-1, 1, -1
   10 dfact = dfact * i
      end

      subroutine mkptab(pwr,nbas,npoly,mxdim,ptab,sump,ndim)
C- Make table of exponents for terms in polynomial
Co ptab, sump, ndim
      implicit none
      integer npoly, nbas, pwr(1), ptab(nbas,1), sump(1)
      integer ib,idim,ip,mxdim,ndim,newdim

      call iinit(ptab,nbas*mxdim)
      ndim = 1
      sump(1) = 0
      do  10  ib = 1, nbas
      do  10  idim = 1, ndim
        newdim = 0
        do  12  ip = 1, pwr(ib)
          if (sump(idim)+ip <= npoly) then
            newdim = newdim+1
            call icopy(nbas,ptab(1,idim),1,ptab(1,ndim+newdim),1)
            sump(ndim+newdim) = sump(idim)+ip
            ptab(ib,ndim+newdim) = ip
          endif
   12   continue
        ndim = ndim+newdim
   10 continue
      if (ndim > mxdim) then
        print *, 'mkptab (stop): mxdim,ndim=',mxdim,ndim
        stop
      endif
      end

      subroutine tfit(nmsh,x1,x2,dx,ndim,nbas,npoly,ptab,coffs)
      implicit none
      double precision x1(1),x2(1),dx(1),coffs(1)
      integer ndim,nbas,npoly,ptab(1),nmsh
      double precision x(20),evpoly,res,xx
      integer ixx(20)
      integer j

c --- Print out the number of points along each coordinate ---
      do  10  j = 1, nmsh
        ixx(j) = 0
        do  20  xx = x1(j), x2(j), dx(j)
   20   ixx(j) = ixx(j)+1
   10 continue
      print 334, (ixx(j), j=1, nmsh)
  334 format(10i4)

      call tfitx(nmsh,nmsh,x,x1,x2,dx,ndim,nbas,npoly,ptab,coffs,res)
      stop

C      if (nbas > 1) goto 20
C      do  10  x = x1, x2, dx
C   10 print 333, x, evpoly(ndim,nbas,ptab,coffs,x)
C  333 format(10f16.8)

      end
      recursive subroutine tfitx(nmsh,j,x,x1,x2,dx,ndim,nbas,npoly,ptab,
     .  coffs,res)
      implicit none
      double precision x(1),x1(1),x2(1),dx(1),coffs(1),res
      integer ndim,nbas,npoly,ptab(1),nmsh,j,k
      double precision evpoly,xx

      do  10  xx = x1(j), x2(j), dx(j)
        x(j) = xx
        if (j == 1) then
          res = evpoly(ndim,nbas,ptab,coffs,x)
          print 333, (x(k),k=1,nmsh), res
  333     format(10f16.8)
        else
          call tfitx(nmsh,j-1,x,x1,x2,dx,ndim,nbas,npoly,ptab,coffs,res)
        endif
   10 continue
      end
      subroutine evfit(ndim,nbas,npoly,pwr,ptab,coffs)
C- Find minimum of quadratic or higher order form
      implicit none
      integer ndim,nbas,npoly,ptab(nbas,ndim),pwr(nbas)
      double precision coffs(1)
      integer jmax
      double precision evpoly, dvpoly, xx
C For multivariate fit:
      integer nbmx,i1mach,a2vec,ip,i1,i2,info,iter,i,j
      parameter (nbmx=10)
      integer ix(10),kpvt(nbmx),ierr,isw,ir,nitm,nlinmx,nlin
      double precision xi(nbmx**2),b(nbmx),a(nbmx,nbmx),w(nbmx,3),
     .  pmin,func,res(0:20),p(0:6*nbmx),xtol,xtoll,grfac,gtol,
     .  dxmn,dxmx,f1d,wk(28)
      external f1d
      character fm*20,outs*80
      procedure(logical) :: a2bin
      logical lexit

      lexit = .false.
      call dpzero(a,nbmx*nbmx)
      call dpzero(b,nbmx)

C --- Minimize quadratic form for nbas=1 or 2 analytically ---
      if (npoly == 2) then
        do  10  i = 1, ndim

C ...   Determine whether term is constant, linear, bilinear
          i1 = 0
          i2 = 0
          do  12  j = 1, nbas
            if (ptab(j,i) /= 0) then
              if (i1 == 0) then
                i1 = j
                if (ptab(j,i) == 2) i2 = j
              else
                i2 = j
              endif
            endif
   12     continue
C         print *, i, i1,i2, coffs(i), (ptab(j,i), j=1,nbas)

C ...     Throw away constant term
          if (i1 == 0 .and. i2 == 0) then
C ...     Linear term becomes -rhs
          else if (i2 == 0) then
            b(i1) = -coffs(i)
C ...     Bilinear term becomes lhs
          else
            if (i1 == i2) then
              a(i1,i1) = 2*coffs(i)
            else
              a(i1,i2) = coffs(i)
              a(i2,i1) = coffs(i)
            endif
          endif
   10   continue
        call dsifa(a,nbmx,nbas,kpvt,info)
        if (info /= 0) call rx('evfit: dsifa cannot factor matrix')
        call dsisl(a,nbmx,nbas,kpvt,b)
        xx = evpoly(ndim,nbas,ptab,coffs,b)
        call awrit1(' %10zf''=0: f=%1;8g at x=',outs,80,0,xx)
        do  14  i = 1, nbas
   14   call awrit1('%a  %g',outs,80,0,b(i))
C        call awrit1(' %10zf''=0: f=%1;8g at',outs,80,0,xx)
C        do  22  i = 1, nbas
C          call awrit1('%%a x%i=%%g',fm,20,0,i)
C          call awrit1(fm,outs,80,0,b(i))
C   22   continue
        call awrit0('%a',outs,80,-i1mach(2))
C       if (nbas > 1) goto 60
      endif

C --- Input ---
   30   continue
        call cwrite(' x=? ',0,4,0)
        read(*,'(a80)',err=999,end=999) outs
        if (outs == '?') then
          print *, 'Enter one of the following:'
          print *, '  x   to evaluate f and derivatives'
          print *, '  /   same, at at current x'
          print *, '  princ for principal axes of local quadratic form'
          if (npoly > 2) then
            print *, '  min[:dxmx] to search for minimum'
            print *, '  min[:dxmx] x=#1,#2.. does the same,',
     .        'but assigns x before search',
     .        ' and exits if search is successful'
          endif
          print *, '  blank or q to quit'
          goto 30
        endif
C ...   Then input for special cases
        if (outs == ' ' .or. outs == 'q')     call cexit(0,1)
        if (outs == 'princ') goto 60
        if (outs(1:3) == 'min')   goto 70
C ...   Then input must be to enter x
        if (outs /= '/') then
          ip = 0
          i = a2vec(outs,len(outs),ip,4,', ',2,3,nbas,ix,b)
          if (i == -1) call fexit(-1,009,'Error parsing for x',0)
          if (i < 0) then
            outs = ' x='
            do  16  i = 1, nbas
   16       call awrit1('%a  %g',outs,80,0,b(i))
            call awrit0('%a',outs,80,-i1mach(2))
          endif
        endif

C --- Function and all derivatives for nbas=1 or 2 ---
      if (nbas <= 2) then
        pmin = func(b)
        if (nbas == 1) then
          call awrit1(' %10zf=%1;8g  f''=',outs,80,0,pmin)
          do  41  i = 1, ndim-1
            ix(1) = i
            xx = dvpoly(ndim,nbas,ptab,coffs,b,ix)
            call awrit1('%a  %g',outs,80,0,xx)
   41     continue
        call awrit1('%a',outs,80,-i1mach(2),i)
C ...   nbas = 1 ... print rows der wrt x, col der wrt y
        else
          print *, ' '
          do  42  i = 0, pwr(1)
          jmax = -1
          do  44  j = 0, pwr(2)
            if (i+j > npoly) goto 44
            jmax = jmax+1
            kpvt(1) = i
            kpvt(2) = j
            res(j) = dvpoly(ndim,nbas,ptab,coffs,b,kpvt)
   44     continue
          print 347, i, (res(j), j=0, jmax)
  347     format(i2,2x,5g15.7:/(' ...',5g15.7))
   42   continue
        endif
        goto 30
      endif

C --- Function values and first derivatives (default) ---
      xx = func(b)
      call awrit1(' %10zf=%1;8g f''=',outs,80,0,xx)
      do  50  i = 1, nbas
        call iinit(kpvt,nbas)
        kpvt(i) = 1
        res(i) = dvpoly(ndim,nbas,ptab,coffs,b,kpvt)
        call awrit1('%a %g',outs,80,0,res(i))
   50 continue
      call awrit1('%a',outs,80,-i1mach(2),res(1))
      goto 30

C ---  Determine principal axes (assuming at minimum) ---
   60 continue
      do  64  i = 1, nbas
      do  64  j = 1, nbas
        call iinit(kpvt,nbas)
        kpvt(i) = 1
        kpvt(j) = kpvt(j)+1
        a(i,j) = dvpoly(ndim,nbas,ptab,coffs,b,kpvt)
   64 continue
      call rs(nbmx,nbas,a,w,1,xi,w(1,2),w(1,3),ierr)
      if (ierr /= 0)
     .  call fexit(-1,009,'Cannot find principal axes',0)
      call dcopy(nbmx*nbmx,xi,1,a,1)
      call awrit0(' Principal axes centred at x=',outs,80,0)
      do  65  i = 1, nbas
   65 call awrit1('%a  %g',outs,80,0,b(i))
      call awrit0('%a',outs,80,-i1mach(2))
      print 321, 'eval', (w(i,1), i=1,nbas)
      do  66  i = 1, nbas
        print 321, '    ', (a(i,j), j=1,nbas)
  321   format(1x,a,5f12.6)
   66 continue
      goto 30

C --- Minimum in multivariate poly, order >2, using Powell's method ---
   70 continue
      call word(outs,2,i,j)
      if (j >= i .and. outs(i:i+1) == 'x=') then
        lexit = .true.
        ip = i+1
        i = a2vec(outs,len(outs),ip,4,', ',2,3,nbas,ix,b)
        if (i /= nbas) call fexit(-1,009,'Error parsing for x',0)
      endif

      call dpzero(xi,nbas**2)

      if (outs(4:4) == ':') then
        i = 0
        if (.not. a2bin(outs(5:),dxmx,4,0,' ',i,-1)) then
          print *, 'failed to parse',outs(4:)
          goto 30
        endif
      else
        dxmx = 0
        do  51  i = 1, nbas
          xi(i+nbas*(i-1)) = b(i)/100
          if (b(i) == 0) xi(i+nbas*(i-1)) = .01d0
          dxmx = max(dxmx,abs(b(i)/20))
   51   continue
        dxmx = max(dxmx,.1d0)
      endif
C     Initial function values
      call dcopy(nbas,b,1,p(1),1)

C      xtol = f1d(0d0,-1d0,ir)
C      p(0) = nbas
C      xtol = f1d(0d0,p,ir)

      nlin = 0
      nlinmx = nbas+1
      nitm = 1000
      xtol = 1d-6
      xtoll= 1d-6
      gtol = 1d-6
      dxmn = 1d-6

      grfac = 1.3d0
C      isw = 1011
C      iter = 0
C      call frpmin(nbas,f1d,p,g,dxmn,dxmx,xtol,gtol,x0,300,isw,
C     .  iter,ir)
      ir = 0
      isw = 111
C     Re-entry point for iteration
      iter = 0
   20 continue
      iter = iter+1
C     gradients at p
      call gradp(ndim,nbas,ptab,coffs,p(1),p(1+nbas))
      call gradzr(nbas,p(1),xi,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
      call dcopy(nbas,p(1),1,b,1)
      if (ir < 0) then

C        call awrit3('iter %i    x: %n:1d',outs,80,
C     .    -i1mach(2),iter,nbas,p(1))
C
C        call awrit3('iter %i grad: %n:1d',outs,80,
C     .    -i1mach(2),iter,nbas,p(1+nbas))

        if (iter > nitm) then
          print *, 'failed to find minimum after', nitm,' iterations'
          call awrit2(' x is now: %n:1g',outs,80,-i1mach(2),nbas,b)
          call awrit1(' gradzr returned ir = %i',' ',120,i1mach(2),ir)
          goto 30
        endif
        if (mod(isw/100,10) < 2 .and. ir == -1) then
          nlin = nlin+1
          if (nlin > nlinmx .and. nlinmx /= 0) then
            nlin = 0
            ir = 0
            call dpzero(xi,nbas**2)
          endif
        endif
        goto 20
      endif

      pmin = func(b)
      call awrit1(' gradzr returned ir = %i',' ',120,i1mach(2),ir)
      call awrit1('%%a  x=%%%i:1g',fm,20,0,nbas)
      call awrit3(' min=%1;8g found after %i iter:',outs,80,
     .  0,pmin,iter,b)
      call awrit1(fm,outs,80,-i1mach(2),b)
C      print 335, pmin, iter, (b(i), i=1,nbas)
C  335 format('min=',f15.7,' found after',i2,' iter: x=',5f11.5)
      if (lexit) goto 999
      goto 30

C --- EOF or ERR encountered
  999 continue
      call cexit(0,1)
      end
      double precision function func(x)
      use pfit_static
      implicit none
      double precision x(2)

      double precision evpoly
      integer iprint,i
      double precision xx


      if (iprint() > 110) then
        xx = evpoly(ndim,nbas,ptab,coffs,x)
        print 335, xx,(x(i),i=1,nbas)
  335 format('f=',f15.7,' x=',5f11.5)
      endif

      func = evpoly(ndim,nbas,ptab,coffs,x)
      end

      double precision function f1d(x,p,ir)
C- Generic function call for projection grad f in a specified direction
C  for subroutine frpmin.
Ci p(0): vector length; 1..n pos'n vec; 2n+1..3n direction vec
Co p(n+1..2n): grad f(posn-vec + x * direction-vec)
Co p(5n+1..6n): (posn-vec + x * direction-vec)
Co f1d: grad f(posn-vec + x * direction-vec) . (unit direction-vec)
Cr If p(0) is zero, f1d returns the total number of function calls.
      use pfit_static
      implicit none
      double precision x,p(0:1)
      integer j,ir,n,i1mach,ipr
      character *80 outs
      double precision dvpoly,ddot
      integer i
      integer igrad(mxbas)

      integer jj
      save jj
      data jj /0/


C ... Return total number of function calls if p(0) is 0
      if (nint(p(0)) == 0) then
        f1d = jj
        return
      endif
      if (nint(p(0)) < 0) then
        jj = 0
        return
      endif

      call getpr(ipr)
      jj = jj+1

C ... (posn-vec + x * direction-vec).  Here n should be same as nbas
      n = nint(p(0))
      do  10  j = 1, n
   10 p(5*n+j) = p(j) + x*p(2*n+j)

C ... Gradient (posn-vec + x * direction-vec) and inner product
      outs = ' '
      do  20  i = 1, nbas
        call iinit(igrad,nbas)
        igrad(i) = 1
        p(nbas+i) = dvpoly(ndim,nbas,ptab,coffs,p(5*n+1),igrad)
        if (ipr >= 80) call awrit1('%a %g',outs,80,0,p(nbas+i))
   20 continue
      if (ipr >= 80) call awrit1('%a',outs,80,-i1mach(2),p(nbas+1))
      f1d = ddot(n,p(nbas+1),1,p(2*nbas+1),1)
      end

      subroutine mapdat(dat,np,nc,ix,iy,sxdef,sydef,mlist,itp,xp,yp)
C- Maps data into (x,y) pair xp,yp
C  if nc=1, xp becomes 0,1,2,3...
      implicit none
C Passed Parameters
      integer np,nc,ix,iy,itp(18)
      character*(*) sxdef,sydef,mlist
      double precision dat(np,*),xp(np),yp(np)
C Local variables
      integer ip,ipr,iv0,ic1,ival,j,i1mach,np0
      character*2 xn
      logical a2bin

      call getpr(ipr)

C --- Case nc=1: Copy col 1 to col 2, make col1 = 0,1,2,3,... ---
      if (nc == 1) then
        do  ip = np, 1, -1
          dat(ip,2) = dat(ip,1)
        end do
        do ip = 1, np
          dat(ip,1) = ip-1
        end do
      endif

C --- Look for data point mapping
      if (mlist == ' ') then
        do  16  ip = 1, np
   16   itp(ip) = ip
        np0 = np
      else
        call mkilst(mlist,np0,itp)
      endif

C --- Copy dat, or subset, into (x,y) ---
      do  14  ip = 1, np0
        if (itp(ip) > np) call rx('mapdat: bad index')
        xp(ip) = dat(itp(ip),ix)
        yp(ip) = dat(itp(ip),iy)
   14 continue
      np = np0

C --- Check for variables mapping ---
      if (sxdef /= ' ' .or. sydef /= ' ') then
        xn = 'x1'
        ic1 = ichar('1')
        call numsyv(iv0)
        do  20  ip = 1, np
          call clrsyv(iv0)
          call addsyv('i',dble(ip),ival)
          call addsyv('x',dat(itp(ip),ix),ival)
          yp(ip) = dat(itp(ip),iy)
          call addsyv('y',yp(ip),ival)
C ---   Load data table ---
          do  22  j = 1, nc
            xn(2:2) = char(ic1+j-1)
            call addsyv(xn,dat(itp(ip),j),ival)
   22     continue
C ---   Change x if sought ---
          if (sxdef(1:1) /= ' ') then
            j = 0
            if (.not. a2bin(sxdef,xp(ip),4,0,' ',j,-1))
     .        call rx('mapdat:  error parsing x')
          endif
C ---   Change y if sought ---
          if (sydef(1:1) /= ' ') then
            j = 0
            if (.not. a2bin(sydef,yp(ip),4,0,' ',j,-1))
     .        call rx('mapdat:  error parsing y')
          endif

          if (ipr >= 100) print 345, ip, xp(ip), yp(ip)
  345     format(' mapdat:  ip=',i4,'  xp,yp=',2f15.6)
          if (ipr >= 100) call shosyv(0,0,0,6)
C         if (ipr >= 100) pause
   20   continue
        call clrsyv(iv0)
      endif

      if (ipr >= 60) then
        call awrit3(' MAPDAT:  data for %i points.  ix=%i  iy=%i',
     .    ' ',80,i1mach(2),np,ix,iy)
        print '(''   i          x(i)              y(i)'')'
        do  30  ip = 1, np
   30   print 347, ip, xp(ip), yp(ip)
  347   format(i4,2f18.6)
      endif

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
