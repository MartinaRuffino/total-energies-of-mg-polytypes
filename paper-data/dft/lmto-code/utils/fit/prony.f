C Least-squares fit, Prony's method (first point at origin)
C See Hildebrand, Introduction to Numerical Analysis, p379.
      subroutine fmain
      implicit none
      integer mxdim,mxbas,mxdat,ncf,nlsq
      parameter (mxdim=40,mxdat=1000)
      logical fail(3)
      double precision dat(mxdat),xp(mxdat),yp(mxdat),a(mxdim*mxdat),
     .  w(mxdat),coffs(mxdim),res(mxdat),r(mxdim*mxdim),q(mxdim*mxdat),
     .  c(4*(mxdat+mxdim)+2),y(mxdat),rhs(mxdat),ddot,xx(mxdat),
     .  delx,err,errmx,ymax
      complex*16 roots(mxdim),expn(mxdim)
      character*72 first,tlist,mlist
      character*1 f2(20)
      logical cmdstr,lsequ,a2bin,period,wline,rdm
      EXTERNAL cmdstr
      integer iarg,n,n1r,n1c,ifi,i,np,iprint,fopng,nargf,i1mach
      integer nfcn,n1,ipivot(mxdim),ifault,nf,mkdlst

C For output wline
      character cc*3, outs*256, sfmt*80, out2*256
      double precision aj,pj,wj,cj,sj
C For tabulating fit on a mesh ...
      integer nmsh
      double precision dx
      logical ltab
      equivalence (f2,first)
C Mapping of abscissa and ordinate
      integer ix(2),it(50),ssvdef,garg,isp(50)
      double precision wsp(50)
      parameter (ssvdef=256)
      character*(ssvdef) sxdef,sydef


      common /static/ coffs,dat,a,w,res,q,r

C#ifdef DEBUG
      double precision snorm,sd,sval,z
      integer m1,kz,nw,mode,j,k,l,ndf,m1p1,n1p1,k2,m
      data m1 /0/, kz /0/, l /1/
C#endif

      tlist = ' '
      mlist = ' '
      ix(1) = 1
      ix(2) = 2
      sxdef = ' '
      sydef = ' '
      period = .false.
      wline = .false.
      ltab = .false.
      n1r = 0
      n1c = 0
      ifi = 10
      goto 10
   20 continue
      print 321
  321 format(
     .' usage:',
     .     'prony [-pr# -p -w -t list -map list -ord "expr"] nexp file'/
     .  7x,'-pr# sets verbosity to #'/
     .  7x,'-p fits periodic functions'/
     .  7x,'-w writes dval and nlfit commands with calculated coffs'/
     .  7x,'-t list tabulates fit for a list of points'/
     .  7x,'-map list maps a subset of the data array'/
     .  7x,'-ord expr substitutes expression for the ordinate;'/
     .  7x,'     use expression of x1,x2,... and i (point no)')
      call cexit(-1,1)

   10 continue
      iarg = 1
   15 continue
      first = ' '
      if (.not. cmdstr(iarg,first)) goto 20
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (lsequ(first,'-pr',3,' ',n)) then
        i = 0
        if (.not. a2bin(f2(4),j,2,0,' ',i,-1)) goto 20
        call pshpr(j)
      else if (lsequ(first,'-v',2,' ',n)) then
        j = 2
        call chrpos(first,'=',19,j)
        n = j+1
        if (.not. a2bin(first,xx,4,0,' ',n,-1)) goto 20
        call addsyv(first(3:j),xx,n)
      else if (lsequ(first,'-map',4,' ',n)) then
        if (.not. cmdstr(iarg+1,mlist)) goto 20
        iarg = iarg+1
      else if (lsequ(first,'-t ',3,' ',n)) then
        if (.not. cmdstr(iarg+1,tlist)) goto 20
        iarg = iarg+1
        call pshpr(0)
      else if (lsequ(first,'-p',2,' ',n)) then
        period = .true.
      else if (lsequ(first,'-w',2,' ',n)) then
        wline = .true.
      else if (lsequ(first,'-col',4,' ',n)) then
        if (garg(' ',iarg+1,2,', ',2,2,i,it,ix) < 2) goto 20
        iarg = iarg+1
      else if (lsequ(first,'-ab',3,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,sxdef)) goto 20
      else if (lsequ(first,'-ord',4,' ',n)) then
        iarg = iarg+1
        if (.not. cmdstr(iarg,sydef)) goto 20
      elseif (garg('-nc=',iarg,2,' ',1,1,i,it,n1c) /= 0) then
        if (i < 0) goto 20
      else
        goto 20
      endif
      iarg = iarg+1
      goto 15
   30 continue

C --- Input nfcn, svdef ---
      read(first,'(i4)',err=20) nfcn
      iarg=iarg+1

C --- Read data from disk file ---
      if (.not. cmdstr(iarg,first)) goto 20
      if (first == '.' .or. (nargf() == iarg)) then
        ifi = i1mach(1)
      else
        ifi = fopng(first,-1,1)
      endif
      if (.not. rdm(ifi,0,mxdat,' ',dat,n1r,n1c))
     .  call fexit(-1,9,'pfit failed to parse file '//first,0)
      if (ifi /= i1mach(1)) close(ifi)
      ncf = n1c
C ... Map the data (a used as a work array)
      call mapdat(dat,n1r,n1c,ix,ix(2),sxdef,sydef,mlist,a,xp,yp)
      if (n1c*n1r > mxdat) stop 'prony: input array too big'
      if (period) then
        if (n1r < 3*nfcn) stop 'prony: must have np >= 3*nexp'
        nlsq = n1r - 2*nfcn
      else
        if (n1r < 2*nfcn) stop 'prony: must have np >= 2*nexp'
        nlsq = n1r - nfcn
      endif
C ... Warning if abscissa not evenly spaced
      dx = 1
      if (n1c > 1) then
        dx = xp(2)-xp(1)
        do  35  i = 2, n1r
          err = dx-(xp(i)-xp(i-1))
          if (dabs(err/dx) > 1d-6)
     .      call awrit4(
     .      ' PRONY: expected delta x between points %i-%i to be '//
     .      '%1;6d but found %1;6d',
     .      ' ',80,i1mach(2),i-1,i,dx,xp(i)-xp(i-1))
   35   continue
      endif

C --- Generate normal matrix for polynomial coefficients ---
      call nrmalf(n1r,nfcn,nlsq,period,yp,a,rhs)
      call dcopy(n1r,1d0,0,w,1)

C --- Least-squares fit for polynomical coefficents ---
      call l2a(nlsq,nfcn,0,1,a,rhs,w,0d0,nlsq,nfcn,
     . n1, ipivot, coffs, res, r, q, c, ifault)

      if (ifault /= 0) call fexit(-1,1,
     .'(/i3,'' PRONY: problem in normal eqns: fault='',f3.0,)',
     .  dble(ifault))

      call dpzero(y,nfcn+1)
      y(1) = 1
      if (period) then
        call daxpy(nfcn,-1d0,coffs,1,y(2),1)
      else
        call daxp(nfcn,-1d0,coffs(nfcn),-1,y(2),1)
      endif
C --- Solve the polynomial ---
      call psolv(nfcn,y,period,roots)
C --- Printout for the exponents ---
      if (iprint() > 30) then
      print 456, 'polynomial',nfcn,n1,dsqrt(ddot(nlsq,res,1,res,1)/nlsq)
  456 format(' Coffs to ',a,':  nfcn=',i2,'  rank=',i3,
     .  '  RMS error=',f12.6)
      print 457, (y(j),j=1,nfcn+1)
  457 format (1x,8g15.8)
      endif
      if (period .and. iprint() >= 10) then
        print *, ' cos(omega)              ...    omega'
        do  42  j = 1, nfcn
        print 457, roots(j), cdlog(roots(j)+cdsqrt(roots(j)**2-1))
   42 continue
      else if (iprint() >= 10) then
        print *, ' Roots to polynomial     ...    Exponents'
        do  40  j = 1, nfcn
          print 457, roots(j), cdlog(roots(j))
   40   continue
      endif
      if (dabs(dx-1) > 1d-8 .and. iprint() >= 10) then
        call awrit1(' Exponents scaled by 1/dx = 1/%1;6d',
     .    ' ',80,i1mach(2),dx)
        if (period) then
          do  44  j = 1, nfcn
   44     print 455, cdlog(roots(j)+cdsqrt(roots(j)**2-1))/dx
        else
          do  46  j = 1, nfcn
   46     print 455, cdlog(roots(j))/dx
        endif
  455   format(31x,8g15.8)
      endif

C --- Least-squares fit for coefficients to exponents ---
      call nrmgen(n1r,nfcn,dx,xp,roots,period,a)
      nf = nfcn
      if (period) nf = 2*nfcn
      call l2a(n1r,nf,0,1,a,yp,w,0d0,n1r,nf,
     . n1, ipivot, coffs, res, r, q, c, ifault)
      ifi = i1mach(2)
      if (period .and. iprint() > 30) then
       print 456, 'exponents',nfcn,n1,dsqrt(ddot(nlsq,res,1,res,1)/nlsq)
       call awrit1('%i',cc,1,0,nf)
       call awrit1('%'//cc(1:1)//':1,12;12g',' ',80,ifi,coffs)
      endif
      if (period  .and. iprint() >= 10) then
        call awrit1('%i',cc,1,0,nfcn)
        do  47  j = 1, nfcn
   47   xx(j) = dsqrt(coffs(j)**2+coffs(j+nfcn)**2)
        call awrit1(' amplitudes  :%'//cc(1:1)//':1,8;8g',' ',80,ifi,xx)
        do  48  j = 1, nfcn
   48   xx(j) = datan2(-coffs(j+nfcn),coffs(j))
        call awrit1(' phases (cos):%'//cc(1:1)//':1,8;8g',' ',80,ifi,xx)
        do  49  j = 1, nfcn
   49   xx(j) = datan2(coffs(j),coffs(j+nfcn))
        call awrit1(' phases (sin):%'//cc(1:1)//':1,8;8g',' ',80,ifi,xx)
      endif
      if (.not. period  .and. iprint() >= 10) then
       print 456, 'exponents',nfcn,n1,dsqrt(ddot(nlsq,res,1,res,1)/nlsq)
       do  j = 1, nfcn
         print 557, coffs(j)
  557    format(30x,1pe15.6)
       enddo
      endif

C --- Printout of fit compared with data ---
      if (iprint() > 10) then
        print 459
  459   format (/2x,
     .    'i    Abscissa   Observed        Predicted      Residual')
        errmx = 0d0
        ymax = 0
        do  50  i = 1, n1r
          ymax = max(ymax,dabs(yp(i)))
          print 458, i, xp(i), yp(i),yp(i)-res(i),res(i)
  458     format (i3,f12.6,2g16.8,g12.4,2x,5g12.6:/'   ...',7g12.6)
          errmx = max(errmx,dabs(res(i)))
   50   continue
        snorm = dsqrt(ddot(n1r,res,1,res,1)/n1r)
        sval  = dsqrt(ddot(n1r,yp,1,yp,1)/n1r)
        call awrit3(' RMS err=%1,3;3g  Max err=%1,3;3g  RMS err/RMS'//
     .    ' val=%1,3;3g',' ',80,i1mach(2),snorm,errmx,snorm/sval)
      endif

C --- Tabulate fit on a uniform mesh ---
      if (tlist /= ' ') then
        np = mkdlst(tlist,0d0,mxdat,xp)
        call nrmgen(np,nfcn,dx,xp,roots,period,a)
        call poppr
        call iotab(xp,np,nf,coffs,a)
        call cexit(0,1)
      endif

C --- Make nlfit and dval commands ---
      if (wline) then
        outs = 'dval x=$x'
        out2 = outs
        print *, ' '
        if (period) then
          do  510  j = 1, nfcn
            aj = dsqrt(coffs(j)**2+coffs(j+nfcn)**2)
            pj = datan2(-coffs(j+nfcn),coffs(j))
            wj = dimag(cdlog(roots(j)+cdsqrt(roots(j)**2-1))/dx)
            call awrit6('%a  a%i=%1;6g p%i=%1;6g w%i=%1;6g',outs,256,
     .        0,j,aj,j,pj,j,wj)
            call awrit6('%a  c%i=%1;6g s%i=%1;6g w%i=%1;6g',out2,256,
     .        0,j,coffs(j),j,coffs(j+nfcn),j,wj)
  510     continue
          cc = '%a '
          do  520  j = 1, nfcn
            call awrit3(cc//'"a%i*cos(w%i*x+p%i)"',outs,256,0,j,j,j)
            call awrit4(cc//'"c%i*cos(w%i*x)+s%i*sin(w%i*x)"',out2,
     .        256,0,j,j,j,j)
            cc = '%a+'
  520     continue
          call awrit0('%a',outs,256,-i1mach(2))
          call awrit0('%a',out2,256,-i1mach(2))

          call awrit1('nlfit -b=1%a',outs,256,0,j)
          call awrit1('nlfit -b=1%a',out2,256,0,j)
          cc = '%a '
          do  530  j = 1, nfcn
            call awrit3(cc//'a%i,w%i,p%i',outs,256,0,j,j,j)
            call awrit3(cc//'c%i,s%i,w%i',out2,256,0,j,j,j)
            cc = '%a,'
  530     continue
          call awrit1('%a '//first,outs,256,0,j)
          call awrit0('%a',outs,256,-i1mach(2))
          call awrit1('%a '//first,out2,256,0,j)
          call awrit0('%a',out2,256,-i1mach(2))
        else

          do  610  j = 1, nfcn
            aj = coffs(j)
            pj = datan2(-coffs(j+nfcn),coffs(j))
            wj = cdlog(roots(j))/dx
            call awrit4('%a  a%i=%1;6g w%i=%1;6g',outs,256,
     .        0,j,aj,j,wj)
  610     continue
          cc = '%a '
          do  620  j = 1, nfcn
            call awrit3(cc//'"a%i*exp(w%i*x)"',outs,256,0,j,j,j)
            cc = '%a+'
  620     continue
          call awrit0('%a',outs,256,-i1mach(2))

          call awrit1('nlfit -b=1%a',outs,256,0,j)
          cc = '%a '
          do  630  j = 1, nfcn
            call awrit2(cc//'a%i,w%i',outs,256,0,j,j)
            cc = '%a,'
  630     continue
          call awrit1('%a '//first,outs,256,0,j)
          call awrit0('%a',outs,256,-i1mach(2))


          call rx('-w not set up now')
        endif
      endif

      end
      subroutine iotab(xp,np,nf,coffs,norm)
      implicit none
      integer np,nf,i,iprint,i1mach
      double precision norm(np,nf),xp(np),coffs(nf),ddot

      if (iprint() < 30) return

C --- Printout of fit compared with data ---
C      print 459
C  459 format (/2x,'  Abscissa   Fn Value')
C     print '(i4,i2)', np, 2
      call awrit2('%% rows %i cols %i',' ',80,i1mach(2),np,2)
      do  10  i = 1, np
        print 458, xp(i), ddot(nf,coffs,1,norm(i,1),np)
  458   format (f12.6,2g16.8,g12.4,2x,5g12.6:/'   ...',7g12.6)
   10 continue

      end
      subroutine prm(ifi,fmt,s,nr,nc)
      integer nr,nc
      double precision s(nr,nc)
      character*(*) fmt
      integer i,j
      call awrit2('%% rows %i cols %i',' ',80,ifi,nr,nc)
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j), j=1,nc)
      end
      subroutine nrmalf(np,nfcn,nlsq,period,y,norm,rhs)
C- Normal matrix for polynomical coefficients
      implicit none
      integer np,nfcn,ip,nlsq,iprint,i1mach
      logical period
      double precision norm(nlsq,nfcn),y(np),rhs(nlsq)

      call dpzero(norm,nlsq*nfcn)
      if (period) then
       call dcopy(nlsq,y,1,rhs,1)
       call daxpy(nlsq,1d0,y(2*nfcn+1),1,rhs,1)
        do  10  ip = 1, nlsq
          call dcopy(nfcn,y(ip+1),1,norm(ip,1),nlsq)
          call daxp(nfcn,1d0,y(ip-1+2*nfcn),-1,norm(ip,1),nlsq)
   10   continue
      else
        call dcopy(nlsq,y(nfcn+1),1,rhs,1)
        do  20  ip = 1, nlsq
   20   call dcopy(nfcn,y(ip),1,norm(ip,1),nlsq)
      endif

C --- Printout ---
      if (iprint() >= 100) then
        print *, 'nrmalf: normal matrix'
        call prm(i1mach(2),'(5f12.5)',norm,nlsq,nfcn)
        print *, 'nrmalf: rhs'
        call prm(i1mach(2),'(5f12.5)',rhs,nlsq,1)
      endif

      end
      subroutine psolv(nord,alfa,period,roots)
      implicit none
      integer nord
      logical period
      double precision alfa(nord+1),aa(10)
      integer iter,ir,n,iprint,falsi
      double precision xl,xr,fl,fr,poly,x,fx
      external poly
      complex*16 roots(3)
      common /sharep/ aa,n

      if (period) then
        if (nord == 1) then
          roots(1) = -alfa(2)/alfa(1)
        elseif (nord == 2) then
          aa(1) = alfa(1)*2
          aa(2) = alfa(2)
          aa(3) = alfa(3) - alfa(1)
          roots(1) = -aa(2)/2/aa(1) +
     .      cdsqrt(dcmplx((aa(2)/2/aa(1))**2 - aa(3)/aa(1), 0d0))
          roots(2) = -aa(2)/2/aa(1) -
     .      cdsqrt(dcmplx((aa(2)/2/aa(1))**2 - aa(3)/aa(1), 0d0))
        elseif (nord == 3) then
          aa(1) = alfa(1)*4
          aa(2) = alfa(2)*2
          aa(3) = alfa(3) - 3*alfa(1)
          aa(4) = alfa(4) - alfa(2)

          if (iprint() >= 50) print 457, (aa(ir), ir=1,4)
  457     format (/' Coffs to polynomial mapped from Chebychev coffs:'/
     .      1x,8g15.8)

          n = nord
          xl = -5
          xr =  5
          call pshpr(iprint()-35)
          call mnbrak(xl,x,xr,fl,fx,fr,poly,30d0,0,100,iter,ir)
          call poppr
          if (iter >= 100) call rx('prony: cannot bracket real root')
          if (falsi(xl,xr,fl,fr,poly,aa,ir,0d0,100,x,fx) < 0)
     .      call rx('prony: falsi cannot find root')
          if (iprint() >= 60) print 345, x
  345     format(' psolv: falsi found root',f20.8)
          roots(1) = x
C ...     Extract coffs to 2nd order poly
          xl = -1
          xr = 1
          if (dabs(xl-x) < 1d-5 .or. dabs(xr-x) < 1d-5) then
            xl = xl*2
            xr = xr*2
          endif
          if (aa(4) == 0) then
            call rx('fix special case zero root')
          endif
          fx = poly(0d0,0,ir)/(0d0-x)
          fl = poly(xl,0,ir)/(xl-x)
          fr = poly(xr,0,ir)/(xr-x)

          aa(1) = (fr+fl-2*fx)/(xr-xl)
          aa(2) = (fr-fl)/(xr-xl)
          aa(3) = fx

          roots(2) = -aa(2)/2/aa(1) +
     .      cdsqrt(dcmplx((aa(2)/2/aa(1))**2 - aa(3)/aa(1), 0d0))
          roots(3) = -aa(2)/2/aa(1) -
     .      cdsqrt(dcmplx((aa(2)/2/aa(1))**2 - aa(3)/aa(1), 0d0))

        else
          call rx('prony not set up for this order root')
        endif
        return
      endif

      if (nord == 1) then
        roots(1) = -alfa(2)/alfa(1)
      elseif (nord == 2) then
        roots(1) = -alfa(2)/2/alfa(1) +
     .    cdsqrt(dcmplx((alfa(2)/2/alfa(1))**2 - alfa(3)/alfa(1), 0d0))
        roots(2) = -alfa(2)/2/alfa(1) -
     .    cdsqrt(dcmplx((alfa(2)/2/alfa(1))**2 - alfa(3)/alfa(1), 0d0))
      else
        stop 'not set up for this order root'
      endif
      end

      double precision function poly(z,k,m)
      implicit none
      integer n,m,k,i
      double precision aa(10),z,xx
      common /sharep/ aa,n

      m = 1
      xx = aa(1)
      do  10  i = 1, n
   10 xx = xx*z + aa(i+1)
      poly = xx
      end

      subroutine nrmgen(np,nfcn,dx,xp,roots,period,norm)
C- Normal matrix for coefficients, given roots

      implicit none
      integer np,nfcn,ip,jp,i,j,iprint,i1mach
      double precision norm(np,nfcn),xp(np),dx
      logical period
      complex*16 roots(1),w

      do  10  i = 1, nfcn
        if ( dimag(roots(i)/dble(roots(i))) > 1d-10)
     .    call rx('complex roots')
   10 continue

C --- Real roots ---
      if (period) then
        if (iprint() >= 30) print *,
     .    'nrmgen: functions cos(x(i)*omega_1..j), sin(x(i)*omega_1..j)'
        do  22  j = 1, nfcn
          w = cdlog(roots(j)+cdsqrt(roots(j)**2-1))/dx
          do  22  ip = 1, np
          jp = ip-1
          norm(ip,j)      = (cdexp(xp(ip)*w)+cdexp(xp(ip)*dconjg(w)))/2
          norm(ip,j+nfcn) = (cdexp(xp(ip)*w)-cdexp(xp(ip)*dconjg(w)))/
     .      (0d0,2d0)
   22   continue
        if (iprint() >= 50)call prm(i1mach(2),'(8f14.7)',norm,np,nfcn*2)
      else
        print *,'nrmgen: functions root_j**(x(i)/dx)'
        do  20  ip = 1, np
          if (xp(ip) /= 0) then
            do  21  j = 1, nfcn
   21       norm(ip,j) = roots(j)**(xp(ip)/dx)
          else
            do  24  j = 1, nfcn
   24       norm(ip,j) = 1d0
          endif
   20   continue
        if (iprint() >= 50) call prm(i1mach(2),'(8f14.7)',norm,np,nfcn)
      endif

      end
      subroutine daxp(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c     Taken from daxpy, start when incx or incy <0 is still at 1.
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,n
c
      ix = 1
      iy = 1
      do  10  i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end
      subroutine mapdat(dat,np,nc,ix,iy,sxdef,sydef,mlist,itp,xp,yp)
C- Maps data into (x,y) pair xp,yp
C  if nc=1, xp becomes 0,1,2,3...
      implicit none
C Passed Parameters
      integer np,nc,ix,iy,ixv(4),it(4),itp(18)
      character*(*) sxdef,sydef,mlist
      double precision dat(np,*),xp(np),yp(np)
C Local variables
      integer ip,ipr,iv0,ic1,ival,j,i1mach,np0
      character*2 xn
      logical a2bin

      call getpr(ipr)

C --- Case nc=1: Copy col 1 to col 2, make col1 = 0,1,2,3,... ---
      if (nc == 1) then
        do  10  ip = np, 1, -1
   10   dat(ip,2) = dat(ip,1)
        do  12  ip = 1, np
   12   dat(ip,1) = ip-1
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
          if (ipr >= 100) pause
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
