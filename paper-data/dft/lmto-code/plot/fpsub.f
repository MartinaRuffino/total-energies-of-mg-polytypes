       subroutine dupdat(nr,nc,nrd,ncd,ndup1,ndup2,nx1,nx2,dat,dupd)
C- Duplicates xy data (for contour plots)
C  nx1 : 1 if 1st and last rows are duplicates of each other
C  nx2 : 1 if 1st and last cols are duplicates of each other
      implicit none
      integer nr,nc,nrd,ncd,ndup1,ndup2,nx1,nx2
      double precision dat(nr,nc),dupd(nrd,ncd)
      integer idr,idc,iaddr,iaddc,ir,ic

      call dmscop(dupd,nrd,dat,nr,1,nr,1,nc,1,1,1d0)

C      call prmx('dat',dat,nr,nr,nc)

      iaddr = 0
      do  idr = 1, ndup1

C     print *, '!!', idr, iaddr
      do  ir = 1, nr
C       print *, ir, ir + iaddr

        iaddc = 0
        do  idc = 1, ndup2
C       print *, '!!', idc, iaddc
        do  ic = 1, nc
          dupd(ir+iaddr,ic+iaddc) = dat(ir,ic)
C         print *, ir, ir+iaddr, ic, ic+iaddc
        enddo
        iaddc = iaddc + nc-nx2
        enddo


      enddo
      iaddr = iaddr + nr-nx1
      enddo

C     call prmx('dup',dupd,nrd,nrd,ncd)

C     stop
      end

      subroutine plboun(xp,yp,nr,xmin,xmax,ymin,ymax)
C- Calculates plot bounds for data points
      implicit none
      integer nr
      double precision xp(nr), yp(nr)
      double precision xmin,xmax,ymin,ymax
      integer ir
      do  10  ir = 1, nr
        xmax = max(xmax,xp(ir))
        xmin = min(xmin,xp(ir))
        ymax = max(ymax,yp(ir))
        ymin = min(ymin,yp(ir))
   10 continue
      end
      subroutine expand(nt,it,xxv,dat,nr,nc)
C- Expands a list of points into a table
C  nt = -2 -> data array has no elements : garg called with "nc~" ... no array elements
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
      if (nt == -2) then  ! Data array has no elements
        nr = 0
        return
      endif
C --- Expand xxv into table ----
   10 i = i+1
      ir = ir+1
      dat(ir) = xxv(i)
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
        if (dx*(x0*(1-1d-12)-x1) < 0) then
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
C         print *, nr,nc
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
      subroutine mapdat(dat,nr,np,nc,ix,iy,ie,iw,iw2,iw3,sxdef,sydef,
     .  sfxdef,sincl,sexcl,nx,ny,pntlst,lsort,xxpl,logx,logy,itp,wk,
     .  xp,yp,ep,wp,wp2,wp3,fxp)
C- Maps data into (x,y) pair
C ----------------------------------------------------------------------
Ci   dat   :input data
Ci   nr    :number of rows in dat
Ci   nc    :number of columns in dat
Ci   ix    :column in dat to which assign (unmapped) abscissa
Ci   iy    :column in dat to which assign (unmapped) ordinate
Ci   ie    :column in dat to which assign (unmapped) error data
Ci   sxdef :string with expression redefining abscissa
Ci   sydef :string with expression redefining ordinate
Ci   sfxdef:function to map abscissa to function of absiccsa?
Ci   sincl :string specifying points to include in mapped list
Ci   sexcl :string specifying points to exclude in mapped list
Ci   nx    :T, normalize abscissa
Ci   ny    :T, normalize ordinate
Ci   pntlst:if not blank, a list of points specifying map
Ci         :Rules are described in mkilst.f
Ci   lsort :T, sort (mapped) data by increasing values of abscissa
Ci   xxpl  :parameters for interpolation of data to a mesh
Ci         :xxpl(1..3): xmin,xmax,dx;
Ci         :xxpl(4): <> if rational;  xxpl(5): order
Ci   logx  :if interpolating, inperpolate logx when logx=T
Ci   logy  :if interpolating, inperpolate logy when logy=T
Ci   itp   :work array of length at least np
Ci   wk    :work array of length at least np
Co Outputs
Co   xp    :mapped abscissa
Co   yp    :mapped ordinates
Co   ep    :mapped third column (error, weight, or z axis for 3D)
Co   wp    :mapped third column (second weight)
Co   wp2   :mapped fourth column (third weight)
Co   wp3   :mapped fifth column (third weight)
Ci   fxp   :mapped function of abscissa
Cr Remarks
Cr   if nc=1, first column becomes 1,2,3...
C    if xxpl nonzero, data is interpolated onto a mesh after mapping
Cu Updates
Cu   17 May 07 Interpolation for logx,logy done on log scale;
Cu             try to stabilize points used for interpolation
Cu   12 Jun 02 Replaced map list with string pntlst
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nr,np,nc,ix,iy,ie,iw,iw2,iw3,itp(np)
      logical nx,ny,lsort,logx,logy
      character*(*) sxdef,sydef,sfxdef,sincl,sexcl,pntlst
      double precision dat(nr,*),xp(np),yp(np),ep(np),wp(np),wp2(np),
     .  wp3(np),fxp(np),wk(1),xxpl(5)
C Local variables
      integer ip,ipr,iv0,ival,j,np0,i1mach
      character*4 xn, outs*80
      logical a2bin,lrat,log1,log2
      double precision xmin,xmax,ymin,ymax,d1mach,xx
      integer jx,nitrp,j1,jp,nxclud,ii
      double precision dy,dymx,x

      call getpr(ipr)
      nxclud = 0

C --- Case nc=1: Copy col 1 to col 2, make col1 = 1,2,3,... ---
      if (nc == 1) then
        do  10  ip = np, 1, -1
   10   dat(ip,2) = dat(ip,1)
        do  12  ip = 1, np
   12   dat(ip,1) = ip
      endif

C --- Make index for map, if a map exists, else index maps i->i ---
      if (pntlst == ' ') then
        do  16  ip = 1, np
   16   itp(ip) = ip
        np0 = np
      else
        call mkils0(pntlst,np0,itp)
        if (np0 > np)
     .    call rx('mapdat: number of points exceeds data size')
        call mkilst(pntlst,np0,itp)
      endif

C --- Using index, map dat into (x,y) ---
      ii = 0
      do  14  ip = 1, np0
        if (itp(ip) <= np) then
          ii = ii+1
          xp(ii) = dat(itp(ip),ix)
          yp(ii) = dat(itp(ip),iy)
          ep(ii) = dat(itp(ip),ie)
          wp(ii) = dat(itp(ip),iw)
          if (iw2 > 0) wp2(ii) = dat(itp(ip),iw2)
          if (iw3 > 0) wp3(ii) = dat(itp(ip),iw3)
        endif
   14 continue
      if (ii /= np0 .and. ipr > 30) then
        call awrit2(' mapdat (warning): only %i allowed points in map'//
     .    ' of %i  expected',' ',80,i1mach(2),ii,np0)
      endif
      np0 = ii
      np = np0

C --- Map expressions redefining ordinate or abscissa ---
      if (sxdef /= ' ' .or. sydef /= ' ' .or. sincl /= ' ' .or.
     .   sfxdef /= ' ') then
        call numsyv(iv0)
        jp = 0
        do  20  ip = 1, np
          call clrsyv(iv0)
          call addsyv('i',dble(ip),ival)
          call addsyv('x',dat(itp(ip),ix),ival)
          yp(ip) = dat(itp(ip),iy)
          ep(ip) = dat(itp(ip),ie)
          wp(ip) = dat(itp(ip),iw)
          call addsyv('y',yp(ip),ival)
          call addsyv('e',ep(ip),ival)
          call addsyv('w',wp(ip),ival)
C     --- Load data table ---
          do  22  j = 1, nc
            ii = 1
            xn = 'x   '
            call bin2a('',0,0,j,2,0,4,xn,ii)
            call addsyv(xn,dat(itp(ip),j),ival)
   22     continue
C     --- Change x if sought ---
          if (sxdef(1:1) /= ' ') then
            j = 0
            if (.not. a2bin(sxdef,xp(ip),4,0,' ',j,-1))
     .        call rx('mapdat:  error parsing x')
          endif
C     --- Change y if sought ---
          if (sydef(1:1) /= ' ') then
            j = 0
            if (.not. a2bin(sydef,yp(ip),4,0,' ',j,-1))
     .        call rx('mapdat:  error parsing y')
          endif
C     --- map abscissa to function of abscissa
          if (sfxdef(1:1) /= ' ') then
            j = 0
            if (.not. a2bin(sfxdef,fxp(ip),4,0,' ',j,-1))
     .        call rx('mapdat:  error parsing y')
          endif
C     --- Exclude points if not satisfied expressions sincl and sexcl ---
          log1 = .false.
          log2 = .false.
C ...     this resets the variable table for altered xp and yp
          if (sincl /= ' ' .or. sexcl /= ' ') then
            call lodsyv('x',1,xp(ip),ival)
            call lodsyv('y',1,yp(ip),ival)
            call lodsyv('e',1,ep(ip),ival)
            call lodsyv('w',1,wp(ip),ival)
          endif
          if (sincl /= ' ') then
            j = 0
C           call shosyv(0,0,0,i1mach(2))
            if (.not. a2bin(sincl,log1,0,0,' ',j,-1))
     .        call rx('mapdat:  error parsing sincl')
            log1 = .not. log1
          endif
          if (sexcl /= ' ') then
            j = 0
            if (.not. a2bin(sexcl,log2,0,0,' ',j,-1))
     .        call rx('mapdat:  error parsing sexcl')
          endif
          if (log1 .or. log2) then
            if (ipr >= 60) print 346, ip, xp(ip), yp(ip)
  346       format(' mapdat:  exclude ip=',i4,'  xp,yp=',2f15.6)
            nxclud = nxclud+1
          else
            jp = jp+1
            if (ipr >= 60) print 345, ip, xp(ip), yp(ip), ep(ip)
  345       format(' mapdat:  ip=',i4,'  xp,yp=',4f15.6)
            if (ipr >= 100) then
              call shosyv(0,0,0,i1mach(2))
C             pause
            endif
            xp(jp) = xp(ip)
            yp(jp) = yp(ip)
            ep(jp) = ep(ip)
            wp(jp) = wp(ip)
            if (sfxdef(1:1) /= ' ') fxp(jp) = fxp(ip)
          endif
   20   continue
        np = jp
        call clrsyv(iv0)
      endif

C --- Normalize abscissas to (0,1) if nx set ---
      if (nx) then
        xmin = d1mach(2)/4
        xmax = -xmin
        do  30  ip = 1, np
        xmax = max(xmax,xp(ip))
   30   xmin = min(xmin,xp(ip))
        if (ipr >= 40) print 333, xmin, xmax, 'abscissa'
  333   format(' mapdat: renormalize limits',2f12.5,' to (0,1)--',a)
        do  33  ip = 1, np
   33   xp(ip) = (xp(ip)-xmin)/(xmax-xmin)
      endif

C --- Normalize ordinates to (0,1) if ny set ---
      if (ny) then
        ymin = d1mach(2)/4
        ymax = -ymin
        do  40  ip = 1, np
        ymax = max(ymax,yp(ip))
   40   ymin = min(ymin,yp(ip))
        if (ipr >= 40) print 333, ymin, ymax, 'ordinate'
        do  43  ip = 1, np
   43   yp(ip) = (yp(ip)-ymin)/(ymax-ymin)
      endif

C --- Sort mapped data points by increasing abscissa ---
      if (lsort) then
        call dvshel(1,np,xp,itp,1)
        call dpcopy(xp,wk,1,np,1d0)
        do  51  ip = 1, np
   51   xp(ip) = wk(itp(ip)+1)
        call dpcopy(yp,wk,1,np,1d0)
        do  52  ip = 1, np
   52   yp(ip) = wk(itp(ip)+1)
        call dpcopy(ep,wk,1,np,1d0)
        do  53  ip = 1, np
   53   ep(ip) = wk(itp(ip)+1)
        call dpcopy(wp,wk,1,np,1d0)
        do  54  ip = 1, np
   54   wp(ip) = wk(itp(ip)+1)
        if (ipr >= 50)
     .    call awrit1(' mapdat: sorted %i points',' ',80,i1mach(2),np)
C        do 55 ip = 1, np
C   55   print 346, ip, xp(ip), yp(ip)
C  346   format(' mapda:  ip=',i4,'  xp,yp=',2f12.6)
      endif

      call awrit1(' mapdat: %i points',outs,80,0,np)
      if (nxclud /= 0)
     .  call awrit1('%a data, %i excluded',outs,80,0,nxclud)

C --- Map into smoothed data by polynomial interpolation ---
      if (xxpl(3) /= 0) then
        jx = 0
        dymx = 1d-10
        nitrp = xxpl(5)
        if (nitrp <= 0) nitrp = min(10,np)
        lrat = xxpl(4) /= 0
        np0 = 0
        if (logx) then
          do  ip = 1, np
            xp(ip) = dlog(xp(ip))
          enddo
        endif
        if (logy) then
          do  ip = 1, np
            yp(ip) = dlog(yp(ip))
          enddo
        endif
        do  60  x = xxpl(1), xxpl(2), xxpl(3)
          np0 = np0+1
          if (logx) then
            xx = dlog(x)
          else
            xx = x
          endif
          if (lrat) then
C           Problem: discontintities between points.
C           call hunty(xp,np,nitrp,xx,jx,j1)
C           Solution: pick midway between enclosing points,
C           so change of interpolation occurs at points.
            call hunty(xp,np,2,xx,jx,j1)
            if (min(xp(j1),xp(j1+1)) < xx .and.
     .          max(xp(j1),xp(j1+1)) > xx) then
              call hunty(xp,np,nitrp,xp(j1)/2+xp(j1+1)/2,jx,j1)
            else
              call hunty(xp,np,nitrp,xx,jx,j1)
            endif
C           call hunty(xp,np,nitrp,xx,jx,j1)
            call ratint(xp(j1),yp(j1),min(nitrp,np-j1+1),xx,wk(np0),dy)
C           print *, np0,j1,min(nitrp,np-j1+1),x,wk(np0),dy
          else
            call polint(xp,yp,np,nitrp,xx,dymx,0,jx,wk(np0),dy)
C           print *,  jx,x,wk(np0),dy
          endif
          if (logy) wk(np0) = exp(wk(np0))
   60   continue
        if (logx) then
          do  ip = 1, np
            xp(ip) = exp(xp(ip))
          enddo
        endif
        if (logy) then
          do  ip = 1, np
            yp(ip) = exp(yp(ip))
          enddo
        endif

        if (ipr >= 30) then
          call awrit4('%a: interpolate to mesh%3:1;4g (%i pts)'//
     .      ' order=%i rat=%l',outs,80,0,xxpl,np0,nitrp,lrat)
          call awrit0('%a',outs,80,-i1mach(2))
          if (ipr >= 50) print '(a)', ' points:'
        endif
        np0 = 0
        do  64  x = xxpl(1), xxpl(2), xxpl(3)
          np0 = np0+1
          xp(np0) = x
          yp(np0) = wk(np0)
          if (ipr >= 50) print '(2f12.6)', xp(np0), yp(np0)
   64   continue
        np = np0
      else
        if (ipr >= 40) call awrit1('%a%?;n>=50;:;;',
     .    outs,80,-i1mach(2),ipr)
        if (ipr >= 50) then
          do  68  ip = 1, np
   68     print 347, xp(ip), yp(ip)
  347     format(2f12.6)
        endif
      endif

      end
      subroutine ixpand(it,ixv,np,itp)
C- Expands a list of integer points into a table
      implicit none
      integer it(4),np,ip,i,ipr,itp(10),ixv(4),idx,ix0,ix1,i1mach
      character s*80, ch*3
C      data  ch /' ,:'/
      data  ch /',: '/

      call getpr(ipr)
      ip = 0
      i = 0
      s = ' ixpand:  "'

C --- Expand ixv into table ----
   10 i = i+1
      call awrit1('%a%i'//ch(it(i):it(i)),s,80,0,ixv(i))
      ip = ip+1
      itp(ip) = ixv(i)
C ...   Generate points in range
      if (it(i) == 2) then
        ix0  = ixv(i)
        i = i+1
        call awrit1('%a%i'//ch(it(i):it(i)),s,80,0,ixv(i))
        ix1  = ixv(i)
        idx = 1
        if (it(i) == 2) then
          i = i+1
          call awrit1('%a%i'//ch(it(i):it(i)),s,80,0,ixv(i))
          idx = ixv(i)
        endif
   20   continue
        if ((ix0-ix1)*isign(1,idx) <= 0) then
          ip = ip+1
          ix0 = ix0+idx
          itp(ip) = ix0
          goto 20
        endif
        ip = ip-1
      endif
C --- Quit if last: rearrange points according to ip,nc ---
      if (it(i) == 3) then
        np = ip
        call awrit1('%a" generated a map to %i points',s,80,0,np)
        if (ipr == 40) call awrit0('%a',s,-80,-i1mach(2))
        if (ipr > 40) then
          call awrit0('%a:',s,80,-i1mach(2))
          ip = 0
          s = ' '
          call bin2av(' ',0,1,0,itp,2,0,np-1,' ',len(s),s,ip)
          call awrit0('%a',s,80,-i1mach(2))
        endif
        return
      endif
      goto 10
      end
      subroutine map3d(dat,nr,np,nc,p3d,sxdef,sydef,sincl,sexcl,it,ixv,
     .  itp,wk,xp,yp,zp,xxmin,xxmax,yymin,yymax,zzmin,zzmax,
     .  xmin,xmax,ymin,ymax,zmin,zmax,frm3d,rotm)
C- Maps data into (x,y,z) triplet
C  itp and wk are work arrays of length at least np
C: pi/2,th,-pi/2 rotates by th around y
Co xmin,xmax,ymin,ymax,zmin,zmax,rotm
      implicit none
C Passed Parameters
      integer nr,np,nc,ix,iy,iz,ixv(4),it(4),itp(1)
      character*(*) sxdef,sydef,sincl,sexcl
      double precision dat(nr,1),xp(np),yp(np),zp(np),wk(1),p3d(6),
     .  frm3d(4)
      double precision xmin,xmax,ymin,ymax,zmin,zmax,rotm(3,3),
     .    xxmin,xxmax,yymin,yymax,zzmin,zzmax
C Local variables
      integer ip,ipr,iv0,ival,np0,i1mach
      character*4 xn, outs*80
      logical a2bin,log1,log2
      double precision d1mach
      integer jp,nxclud,ii,i,j,k
      double precision xprj,zprj,yprj,xr,yr,zr,
     .  abg(3),big

      call getpr(ipr)
      big = sqrt(d1mach(2))/100
      nxclud = 0
      ix = 1
      iy = 2
      iz = 3

C --- Make index for map, if a map exists, else index maps i->i ---
      if (it(1) == 0) then
        do  16  ip = 1, np
   16   itp(ip) = ip
        np0 = np
      else
        call ixpand(it,ixv,np0,itp)
      endif

C --- Using index, map dat into (x,y) ---
      ii = 0
      do  14  ip = 1, np0
        if (itp(ip) <= np) then
          ii = ii+1
          xp(ii) = dat(itp(ip),ix)
          yp(ii) = dat(itp(ip),iy)
          zp(ii) = dat(itp(ip),iz)
        endif
   14 continue
      if (ii /= np0 .and. ipr > 30) then
        call awrit2(' map3d (warning): only %i allowed points in map'//
     .    ' of %i  expected',' ',80,i1mach(2),ii,np0)
      endif
      np0 = ii
      np = np0

C --- Map expressions redefining ordinate or abscissa ---
      if (sxdef /= ' ' .or. sydef /= ' ' .or. sincl /= ' ') then
        call numsyv(iv0)
        jp = 0
        do  20  ip = 1, np
          call clrsyv(iv0)
          call addsyv('i',dble(ip),ival)
          call addsyv('x',dat(itp(ip),ix),ival)
          yp(ip) = dat(itp(ip),iy)
          call addsyv('y',yp(ip),ival)
C ---   Load data table ---
          do  22  j = 1, nc
            ii = 1
            xn = 'x   '
            call bin2a('',0,0,j,2,0,4,xn,ii)
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
C ---   Exclude points if not satisfied expressions sincl and sexcl ---
          log1 = .false.
          log2 = .false.
C ...     this resets the variable table for altered xp and yp
          if (sincl /= ' ' .or. sexcl /= ' ') then
            call lodsyv('x',1,xp(ip),ival)
            call lodsyv('y',1,yp(ip),ival)
          endif
          if (sincl /= ' ') then
            j = 0
C           call shosyv(0,0,0,i1mach(2))
            if (.not. a2bin(sincl,log1,0,0,' ',j,-1))
     .        call rx('mapdat:  error parsing sincl')
            log1 = .not. log1
          endif
          if (sexcl /= ' ') then
            j = 0
            if (.not. a2bin(sexcl,log2,0,0,' ',j,-1))
     .        call rx('mapdat:  error parsing sexcl')
          endif
          if (log1 .or. log2) then
            if (ipr >= 60) print 346, ip, xp(ip), yp(ip)
  346       format(' mapdat:  exclude ip=',i4,'  xp,yp=',2f15.6)
            nxclud = nxclud+1
          else
            jp = jp+1
            if (ipr >= 60) print 345, ip, xp(ip), yp(ip)
  345       format(' mapdat:  ip=',i4,'  xp,yp=',2f15.6)
            if (ipr >= 100) then
              call shosyv(0,0,0,i1mach(2))
C             pause
            endif
            xp(jp) = xp(ip)
            yp(jp) = yp(ip)
            zp(jp) = zp(ip)
          endif
   20   continue
        np = jp
        call clrsyv(iv0)
      endif

      call awrit1(' map3d: %i points',outs,80,0,np)
      if (nxclud /= 0)
     .  call awrit1('%a data, %i excluded',outs,80,0,nxclud)
      if (ipr >= 40) call awrit0('%a',outs,80,-i1mach(2))

C --- Rotate and shift data, get plot bounds ---
      if (ipr >= 10) then
        call rm2eua(rotm,abg(1),abg(2),abg(3))
        call awrit2(' map3d: shift y,x,z = (%3:1;5d )  Euler angles'//
     .    ' = (%3:1;5d)',' ',80,i1mach(2),p3d,abg)
        if (ipr >= 50) print 335, ((rotm(i,j),j=1,3),i=1,3)
  335   format(' Rotation matrix:'/(3f15.9))
      endif
      if (ipr >= 40) print 333
  333 format(6x,'xp',10x,'zp',10x,'yp',10x,'xrot',8x,'zrot',8x,'yrot')
      do  30  ip = 1, np
        xmax = max(xmax,xp(ip))
        xmin = min(xmin,xp(ip))
        ymax = max(ymax,yp(ip))
        ymin = min(ymin,yp(ip))
        zmax = max(zmax,zp(ip))
        zmin = min(zmin,zp(ip))
        if (xxmin /= big) xmin = xxmin
        if (xxmax /= big) xmax = xxmax
        if (yymin /= big) ymin = yymin
        if (yymax /= big) ymax = yymax
        if (zzmin /= big) zmin = zzmin
        if (zzmax /= big) zmax = zzmax
        call prj3dp(xp(ip),yp(ip),zp(ip),.true.,rotm,p3d,
     .    xr,yr,zr,xprj,yprj,zprj)
        if (ipr >= 40) print '(6f12.6)',
     .    xp(ip), zp(ip), yp(ip), xr, zr, yr
        frm3d(1) = min(frm3d(1),xprj)
        frm3d(2) = max(frm3d(2),xprj)
        frm3d(3) = min(frm3d(3),zprj)
        frm3d(4) = max(frm3d(4),zprj)
        xp(ip) = xr
        yp(ip) = yr
        zp(ip) = zr
   30 continue
C --- Use 2D projections for edges from xmin,xmax etc if all known ---
      if (xxmin /= big .and. xxmax /= big .and.
     .    yymin /= big .and. yymax /= big .and.
     .    zzmin /= big .and. zzmax /= big) then
        frm3d(2) = -big
        frm3d(4) = frm3d(2)
        frm3d(1) =  big
        frm3d(3) = frm3d(1)
        abg(1) = xxmin
        do  32  i = 1, 2
        abg(2) = yymin
        do  33  j = 1, 2
        abg(3) = zzmin
        do  34  k = 1, 2
          call prj3dp(abg(1),abg(2),abg(3),.true.,rotm,p3d,
     .    xr,yr,zr,xprj,yprj,zprj)
          frm3d(1) = min(frm3d(1),xprj)
          frm3d(2) = max(frm3d(2),xprj)
          frm3d(3) = min(frm3d(3),zprj)
          frm3d(4) = max(frm3d(4),zprj)
          abg(3) = zzmax
   34   continue
        abg(2) = yymax
   33   continue
        abg(1) = xxmax
   32   continue
      endif
      end
      subroutine frme3d(tx0,tsx,mtx,ntx,rmtx,xnum,xnfmt,xlabel,
     .                  ty0,tsy,mty,nty,rmty,ynum,ynfmt,ylabel,
     .                  tz0,tsz,mtz,ntz,rmtz,znum,znfmt,zlabel,
     .                  xmn,xmx,ymn,ymx,zmn,zmx,p3d,rotm,
     .                  title,lt,cs)
C- 3D frame around the plot and either tic marks or grid lines.
C ----------------------------------------------------------------
Ci Inputs
Ci   tx0:  Force a major tic through this point
Ci         tx0 < plul makes default
Ci   tsx:  spacing in user's units between tic marks along the x axis
Ci         tsx <= 0 makes default of approximately 1/10 grid width
Ci   mtx:  number of major tics per tic, x axis
Ci   ntx:  number of tics, x axis
Ci         ntx < 0 makes default to fill grid
Ci   rmtx: size of major tic, in proportion to plut-plub (.03 is good)
Ci   lt:   (1) line thickness, (2) whether 3 lines drawn or 8
Ci   Similarly for ty0,tsy,mty,nty,rmty; tz0,tsz,mtz,ntz,rmtz
Ci
Co Outputs
Co   Axes are drawn.
Cr Remarks
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer ntx,nty,ntz,mtx,mty,mtz,lt(2)
      double precision tx0,tsx,rmtx,ty0,tsy,rmty,tz0,tsz,rmtz,
     .  xmn,xmx,ymn,ymx,zmn,zmx,p3d(3),rotm(3,3),cs(0:3)
      character*20 xnfmt,ynfmt,znfmt
      character*60 xlabel,ylabel,zlabel,title
      integer xnum,ynum,znum

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .     plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .     plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .       cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .     xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .     linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .                plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .                plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .                     cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .                xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .                lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .                medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      integer i1mach,ipr
C     double precision sx,sy,sz,ticspc
      character*80 s
C      th(x,y,z) = datan2(dsqrt(x*x+z*z),y)
C      ph(x,z)   = datan2(z,x)
C      prjx(x,y,z) = dcos(datan2(z,x))*dsin(datan2(dsqrt(x*x+z*z),y))
C      prjy(x,y,z) = dsin(datan2(z,x))*dsin(datan2(dsqrt(x*x+z*z),y))

      call getpr(ipr)

C --- tic mark spacing ---
C      sx = tsx
C      sy = tsy
C      sz = tsz
C      if (sx <= 0) sx = ticspc(dabs(xmx-xmn)/10,1)
C      if (sy <= 0) sy = ticspc(dabs(ymx-ymn)/10,1)
C      if (sz <= 0) sz = ticspc(dabs(zmx-zmn)/10,1)

C --- Setup ---
      call awrit6(' FRME3D: x=%1;4g %1;4g y=%1;4g %1;4g z=%1;4g %1;4g',
     .  s,80,i1mach(2),xmn,xmx,ymn,ymx,zmn,zmx)
      if (medium == 3)
     .call awrit6('%%FRME3D: x=%1;4g %1;4g y=%1;4g %1;4g z=%1;4g %1;4g',
     .  s,80,stream,xmn,xmx,ymn,ymx,zmn,zmx)
      if (medium == 3) write (stream, 333)
  333 format('/max 0 def ')

      call plntys
      call plntyp(1,lt,0d0,0d0,0d0,0d0)

C --- Draw 3 lines from (xmn,ymn,zmn) ---
      if (lt(2) == 0) then
        call mve3d(xmn,ymn,zmn,rotm,p3d)
        call drw3d(xmx,ymn,zmn,rotm,p3d)
        call mve3d(xmn,ymn,zmn,rotm,p3d)
        call drw3d(xmn,ymx,zmn,rotm,p3d)
        call mve3d(xmn,ymn,zmn,rotm,p3d)
        call drw3d(xmn,ymn,zmx,rotm,p3d)
        call mve3d(xmn,ymn,zmn,rotm,p3d)
      else
C ... zmn face
        call mve3d(xmn,ymn,zmn,rotm,p3d)
        call drw3d(xmx,ymn,zmn,rotm,p3d)
        call drw3d(xmx,ymx,zmn,rotm,p3d)
        call drw3d(xmn,ymx,zmn,rotm,p3d)
        call drw3de(xmn,ymn,zmn,rotm,p3d,-2d0)
C ... zmx face
        call mve3d(xmn,ymn,zmx,rotm,p3d)
        call drw3d(xmx,ymn,zmx,rotm,p3d)
        call drw3d(xmx,ymx,zmx,rotm,p3d)
        call drw3d(xmn,ymx,zmx,rotm,p3d)
        call drw3de(xmn,ymn,zmx,rotm,p3d,-2d0)
C ... ymn face
        call mve3d(xmn,ymn,zmn,rotm,p3d)
        call drw3d(xmx,ymn,zmn,rotm,p3d)
        call drw3d(xmx,ymn,zmx,rotm,p3d)
        call drw3d(xmn,ymn,zmx,rotm,p3d)
        call drw3de(xmn,ymn,zmn,rotm,p3d,-2d0)
C ... ymx face
        call mve3d(xmn,ymx,zmn,rotm,p3d)
        call drw3d(xmx,ymx,zmn,rotm,p3d)
        call drw3d(xmx,ymx,zmx,rotm,p3d)
        call drw3d(xmn,ymx,zmx,rotm,p3d)
        call drw3de(xmn,ymx,zmn,rotm,p3d,-2d0)
      endif

      call plntyg

C     write (stream,'(''% we are done'')')

      if (medium /= 3) return
      end
      subroutine plsy3(xp,yp,zp,rotm,shft,type,syma,symcs)
C- Draw a symbol at a set of (x,y,z) points
C ----------------------------------------------------------------------
Ci Inputs
Ci   xp,yp,zp: coordinates and number
Ci   type: 1 for flat arrow, 2 for cone
Ci   arrow:   tip at (xp,yp,zp);
Ci            syma  (1,2,3): xp,yp,zp of tail
Ci                  (4,5): head length,width (fraction of arrow length)
Ci                  (6,7): phi,theta of (degrees)
Ci   cone:    tip at (xp,yp,zp);
Ci            syma  (1,2,3,4,5): ditto arrow
Ci                  (6,7): phi, angle between umbrella "wires" and
Ci                         theta, angle between cone and shaft (degrees)
Ci   syma: symbol parameters, and number of elts separating parms
Ci   symcs: symbol grey/color scale
Co Outputs
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer type
      double precision syma(7),symcs(0:3),xp,yp,zp,rotm(3,3),shft(3)

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .     plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .     plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .       cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .     xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .     linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .                plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .                plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .                     cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .                xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .                lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .                medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      character*7 synam(8),s*100
      integer is,ipr,nsyma(8),i1mach,awrite,i
      double precision tpi,tpix,rloc(3,3),vprj(3),dx,dy,dz,
     .  phi0,th0,d,ephi(2),dphi,chngy,xx,cs(0:3)

      data tpi /6.283185307179586d0/
      data synam /'farrow','cone',' ',' ',' ',' ',' ',' '/
      data nsyma /7,7,0,0,0,0,0,0/

C --- Setup ---
      if (type == 0) return
      if (type > 2) call fexit(-1,1,
     .  'Exit -1 plsy3: type %i not defined',type)
      call getpr(ipr)
C      is = 0
C      s = ' '
C      call plcrve(medium,pstat,-2d0,s,is)
C ... For white-out
      call dvset(cs,1,4,1d0)
      cs(0) = 3

      s = '%'
      call awrit6(' plsy3: '//synam(type)//'gs=%3:1;5d  attr=%'//
     .  'n:1;5d  x,y,z=%1;4g %1;4g %1;4g',
     .  s(2:),len(s)-1,0,symcs(1),nsyma(type),syma,xp,yp,zp)
      if (ipr >= 40) call awrit0('%a',s(2:),-len(s)-1,-i1mach(2))
      if (medium == 3) call awrit0('%a',s,-len(s),-stream)

      call plntys
      call plntyp(1,ltbld,0d0,0d0,0d0,0d0)

      call mve3d(xp,yp,zp,rotm,shft)
      goto (10,20), iabs(type)
      goto 5

C --- Symbol type flat arrow ---
   10 continue
C ... draw from tail to tip
      call mve3d(xp+syma(1),yp+syma(2),zp+syma(3),rotm,shft)
      call drw3d(xp,yp,zp,rotm,shft)

C ... Arrowhead oriented at angle relative to arrow shaft.  The latter
C     defines a local coordinate system with Euler angles (phi0,th0,0),
C     and tip at origin.  In that coordinate system, arrowhead is
C     along a line (phi,theta) = (syma(6),syma(7)) .
      tpix = tpi/360
      dx = syma(1)
      dy = syma(2)
      dz = syma(3)
      d = dsqrt(dx**2+dy**2+dz**2)
      phi0 = datan2(dy,dx)
      th0 = dacos(dz/d)
      call eua2rm(phi0,th0,0d0,rloc)
C ... Arrowhead corner, rotate from local coords w/ shaft parallel to z
      call arcdrw(d*syma(4),tpix*syma(6),0d0,tpix*syma(7),rloc,-1,
     .  xp,yp,zp,rotm,shft)
C ... Arrowhead middle, global coordinates
      call drw3d(xp+dx*syma(5),yp+dy*syma(5),zp+dz*syma(5),rotm,shft)
C ... Second arrowhead corner
      call arcdrw(d*syma(4),tpix*(syma(6)+180),0d0,tpix*syma(7),rloc,-1,
     .  xp,yp,zp,rotm,shft)
C ... Back to tip
      call drw3d(xp,yp,zp,rotm,shft)
      goto 5

C --- Symbol type cone ---
   20 continue
C ... draw from tail to tip
      call mve3d(xp+syma(1),yp+syma(2),zp+syma(3),rotm,shft)
      call drw3d(xp,yp,zp,rotm,shft)
C ... Cone oriented at angle relative to arrow shaft.  The latter
C     defines a local coordinate system with Euler angles (phi0,th0,0),
C     and tip at origin.  In that coordinate system, cone outline at
C     angle theta = syma(7).  Local phi determined by maximal projection
C     onto plane, or minimal projection onto y.
      tpix = tpi/360
      dx = syma(1)*syma(4)
      dy = syma(2)*syma(4)
      dz = syma(3)*syma(4)
      d = dsqrt(dx**2+dy**2+dz**2)
      if (d == 0) goto 5
      phi0 = 0
      if (dx /= 0 .or. dy /= 0) phi0 = datan2(dy,dx)
      th0 = dacos(dz/d)
      call eua2rm(phi0,th0,0d0,rloc)
C ... Extremal value of phi
      call dgemm('N','N',3,1,3,1d0,rotm,3,syma,3,0d0,vprj,3)
      call ephi3d(xp,yp,zp,rotm,shft,tpix*syma(7),rloc,ephi)
      xx = chngy(xp,yp,zp,syma(1),syma(2),syma(3),rotm,shft)
      if (ipr >= 50) call awrit2(' plsy3: yproj = %1;6d  '//
     .  'ephi =%2:1;6d',' ',80,i1mach(2),xx,ephi)
C ... Case tip pointing away from viewer:
      if (xx < 0) then
C   ... Draw entire arc and white out
        call arcdrw(d,0d0,tpi,tpix*syma(7),rloc,20,xp,yp,zp,rotm,shft)
        call drwre(0d0,0d0,cs)
C   ... Redraw tail-to-tip
        call mve3d(xp+syma(1),yp+syma(2),zp+syma(3),rotm,shft)
        call drw3d(xp,yp,zp,rotm,shft)
C   ... Draw tip-forward arc-tip, if it is visible
        i = 1
        if (ephi(1) == 0 .and. ephi(2) == 0) i = 0
C ... Case tip pointing toward viewer:
      else
C   ... Case partial arc showing ... draw the back side of the arc
        if (ephi(1) /= 0 .or. ephi(2) /= 0) then
          call arcdrw(d,ephi,ephi(2),tpix*syma(7),rloc,
     .      nint(tpi*dabs(ephi(1)-ephi(2))),xp,yp,zp,rotm,shft)
C     ... Set flag to draw tip-forward_arc-tip
          i = 1
C   ... Case whole arc showing ... draw and color entire arc
        else
          call arcdrw(d,0d0,tpi,tpix*syma(7),rloc,20,xp,yp,zp,rotm,shft)
          call drwre(0d0,0d0,symcs)
          ephi(2) = tpi
          i = 0
        endif
      endif

C ... Draw tip-forward_arc-tip, if flag set
      if (i /= 0) then
        if (ephi(2) > ephi(1)) then
          ephi(2) = ephi(2)-tpi
        else
          ephi(2) = ephi(2)+tpi
        endif
        call mve3d(xp,yp,zp,rotm,shft)
        call arcdrw(d,ephi,ephi(2),tpix*syma(7),rloc,
     .    -nint(tpi*dabs(ephi(1)-ephi(2))),xp,yp,zp,rotm,shft)
        call drw3d(xp,yp,zp,rotm,shft)
        call drwre(0d0,0d0,symcs)
      endif

C ... Draw lines from tip to skirt
      if (syma(6) /= 0) then
        xx = dsign(syma(6)*tpix,ephi(2)-ephi(1))
        do  24  dphi = ephi(1), ephi(2), xx
          call mve3d(xp,yp,zp,rotm,shft)
          call arcdrw(d,dphi,0d0,tpix*syma(7),rloc,-1,
     .      xp,yp,zp,rotm,shft)
   24   continue
      endif
      goto 5

C --- Cleanup for this point ---
    5 continue
      is = 0
      s = ' '
      if (medium == 3) then
        is = awrite(' closepath gsave %1;4d setgray fill grestore',
     .    s,len(s),0,symcs(1),0,0,0,0,0,0,0)
      endif
      call plcrve(medium,pstat,-2d0,s,is)
      if (is > 0) write(stream,*) s(1:is)
    2 continue
  430 format(' closepath gsave',f8.4,' setgray fill grestore')
      call plntyg
      return

      end
      double precision function chngy(x,y,z,dx,dy,dz,rotm,shf)
C- Change in yprx in movement from x to x+dx
      implicit none
      double precision x,y,z,dx,dy,dz,rotm(3,3),shf(3)
      double precision xr,yr,zr,xp,yp,zp,xp0,yp0,zp0

      call prj3dp(x,y,z,.true.,rotm,shf,xr,yr,zr,xp0,yp0,zp0)
      call prj3dp(x+dx,y+dy,z+dz,.true.,rotm,shf,xr,yr,zr,xp,yp,zp)
      chngy = yp - yp0
      end
      subroutine arcdrw(d,ph1,ph2,theta,rloc,np,xp,yp,zp,rotm,shft)
C- Draws an arc in local axis (d,ph1..ph2,theta) about xp,yp,zp
Ci np: 0: move pen to ph1; >0 draw arc, moving to starting point
Ci    -1: draw to ph1;    <-1 draw arc, drawing to starting point
      implicit none
      integer np
      double precision d,ph1,ph2,theta,xp,yp,zp,
     .  rloc(3,3),rotm(3,3),shft(3)
      integer mp
      double precision vloc(3),vglob(3),phi

      mp = 0
      phi = ph1
   10 vloc(1) = d*dcos(phi)*dsin(theta)
      vloc(2) = d*dsin(phi)*dsin(theta)
      vloc(3) = d*dcos(theta)
      call dgemm('T','N',3,1,3,1d0,rloc,3,vloc,3,0d0,vglob,3)
C ... First move
      if (mp == 0 .and. np >= 0) then
        call mve3d(xp+vglob(1),yp+vglob(2),zp+vglob(3),rotm,shft)
      else
        call drw3d(xp+vglob(1),yp+vglob(2),zp+vglob(3),rotm,shft)
      endif
      mp = mp+1
      if (mp >= iabs(np)) return
      phi = dble(mp)/dble(iabs(np)-1)*(ph2-ph1) + ph1
      goto 10

      end
      subroutine rotpos(d,phi,theta,rot,vrot)
C- Returns rotated position given spherical coordinates in local axis
      implicit none
      double precision d,phi,theta,rot(3,3),vrot(3)
      double precision vloc(3)
      vloc(1) = d*dcos(phi)*dsin(theta)
      vloc(2) = d*dsin(phi)*dsin(theta)
      vloc(3) = d*dcos(theta)
      call dgemm('T','N',3,1,3,1d0,rot,3,vloc,3,0d0,vrot,3)
      end
      subroutine ephi3d(xp,yp,zp,rotm,shf,theta,rloc,ephi)
C- Determines extremal angles on x-z plane wrt rotation about const. th.
      implicit none
      double precision xp,yp,zp,shf,rotm(3,3),rloc(3,3),theta,ephi(2)
      double precision r(3,3),pi,A,B
*     double precision C,rhs
      double precision vglob(3),dy,phibar,xpr0,ypr0,zpr0,
     .  xpr(-1:1),ypr(-1:1),zpr(-1:1),xr,yr,zr,ph,dph,phi
      integer i,j
      pi = 4*datan(1d0)
      call dgemm('N','T',3,3,3,1d0,rotm,3,rloc,3,0d0,r,3)
C     print 335, ((r(i,j),j=1,3),i=1,3)
C 335 format((3f15.9))

c ...  analytic formula, forgetting additional shift from perspective
C      A =  r(3,2)*r(1,3) - r(1,2)*r(3,3)
C      B = -r(3,1)*r(1,3) + r(1,1)*r(3,3)
C      C =  r(1,1)*r(3,2) - r(1,2)*r(3,1)
C      rhs = -C*dtan(theta)/dsqrt(A**2+B**2)
C      ephi(2) = 0
C      ephi(1) = 0
C      print *, 'rhs is=', rhs
C      if (dabs(rhs) > 1) return
C      ephi(1) = dasin(rhs) - datan2(A,B)
C      ephi(2) = pi - ephi(1) - 2*datan2(A,B)
C      phibar = (ephi(1)+ephi(2))/2
C      dy = (r(2,1)*dcos(phibar)    + r(2,2)*dsin(phibar)
C     .     -r(2,1)*dcos(phibar+pi) - r(2,2)*dsin(phibar+pi))*dsin(theta)
C      if (dy < 0) ephi(2) = ephi(2) - 2*pi
C      print *, 'ephi=',ephi(1), ephi(2), rhs, dy

C --- Determine extrema in projected rotation numerically ---
      call prj3dp(xp,yp,zp,.true.,rotm,shf,xr,yr,zr,xpr0,ypr0,zpr0)
      dph = pi/360
      j = 0
      ephi(1) = 0
      ephi(2) = 0
      do  10  ph = 0, 2*pi+1d-10, dph
        do  12  i = -1, 1
          phi = ph + i*dph
          if (ph == 0 .or. i == 1) then
            call rotpos(1d0,phi,theta,rloc,vglob)
            call prj3dp(xp+vglob(1),yp+vglob(2),zp+vglob(3),.true.,
     .        rotm,shf,xr,yr,zr,xpr(i),ypr(i),zpr(i))
            xpr(i) = xpr(i) - xpr0
            ypr(i) = ypr(i) - ypr0
            zpr(i) = zpr(i) - zpr0
          else
            xpr(i) = xpr(i+1)
            ypr(i) = ypr(i+1)
            zpr(i) = zpr(i+1)
          endif
   12   continue
        if (dabs(xpr(0)) > dabs(zpr(0))) then
          a = zpr(1)/xpr(1)-zpr(0)/xpr(0)
          b = zpr(0)/xpr(0)-zpr(-1)/xpr(-1)
        else
          a = xpr(1)/zpr(1)-xpr(0)/zpr(0)
          b = xpr(0)/zpr(0)-xpr(-1)/zpr(-1)
        endif
        if (a*b <= 0) then
C         print *, a,b
          phibar = ph - dph/2*(a+b)/(a-b)
          if (j >= 1 .and. dabs(phibar-2*pi-ephi(1)) < 1d-4)goto 10
          if (j >= 2) call fexit2(-1,1,' Exit -1 EPHI: '//
     .      'too many ephi:%2:1;6d %1;6d',ephi,ph-dph/2*(a+b)/(a-b))
          j = j+1
          ephi(j) = ph - dph/2*(a+b)/(a-b)
C         print *, ph - dph/2*(a+b)/(a-b),
C    .      xpr(-1)/zpr(-1), xpr(0)/zpr(0), xpr(1)/zpr(1)
        endif
   10 continue
      if (j == 1) call rx('ephi3d: odd number of extrema')
      if (dabs(ephi(1)-ephi(2)) < .01d0) then
        j = 0
        ephi(1) = 0
        ephi(2) = 0
      endif
      if (j == 2) then
        phibar = (ephi(1)+ephi(2))/2
        dy = (r(2,1)*dcos(phibar)    + r(2,2)*dsin(phibar)
     .    -r(2,1)*dcos(phibar+pi) - r(2,2)*dsin(phibar+pi))*dsin(theta)
        if (dy < 0) ephi(2) = ephi(2) - 2*pi
      endif
C     print *, 'ephi=',ephi
      end
      subroutine mve3d(x,y,z,rotm,shft)
C- Move point or draw to (x,y,z)
      implicit none
      double precision x,y,z,rotm(3,3),shft(3),cs(0:3)
      double precision xr,yr,zr,xprj,yprj,zprj
      integer iprint

      call prj3dp(x,y,z,.true.,rotm,shft,xr,yr,zr,xprj,yprj,zprj)
      if (iprint() >= 100) print 333, 'mve3d', x,y,z,xprj,zprj
  333 format(1x,a,5f12.6)
      call mve(xprj,zprj)
      return

      entry drw3d(x,y,z,rotm,shft)
      call prj3dp(x,y,z,.true.,rotm,shft,xr,yr,zr,xprj,yprj,zprj)
      if (iprint() >= 100) print 333, 'drw3d', x,y,z,xprj,zprj
      call drw(xprj,zprj)
      return

      entry drw3de(x,y,z,rotm,shft,cs)
      call prj3dp(x,y,z,.true.,rotm,shft,xr,yr,zr,xprj,yprj,zprj)
      if (iprint() >= 100) print 333, 'drw3de', x,y,z,xprj,zprj
      call drwe(xprj,zprj,cs)

      end
      subroutine prj3dp(x,y,z,sw,rotm,shft,xr,yr,zr,xprj,yprj,zprj)
C- Projection of 3d points onto plane
Ci x,y,z,rotm,shft
Ci sw: T, rotate x,y,z by rotm, shift y by shft before projection
Co xprj,yprj,zprj
      implicit none
      logical sw
      double precision x,y,z,rotm(3,3),shft(3),xprj,yprj,zprj,xr,yr,zr
      double precision datan2,dsin,dcos,pi,th,ph,p(3,2)

      pi = 4*datan(1d0)
      p(1,1) = x
      p(2,1) = y
      p(3,1) = z
      if (sw) then
        call dmpy(rotm,3,1,p,3,1,p(1,2),3,1,3,1,3)
        p(1,1) = p(1,2) + shft(2)
        p(2,1) = p(2,2) + shft(1)
        p(3,1) = p(3,2) + shft(3)
        xr = p(1,1)
        yr = p(2,1)
        zr = p(3,1)
      endif
      if (p(2,1) == 0d0) then
        th = pi/2
      else
        th = datan2(dsqrt(p(1,1)**2+p(3,1)**2),p(2,1))
      endif
      yprj = dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)
      if (th > pi/2) yprj = -yprj
      if (th == 0) then
        xprj = 0
        zprj = 0
        ph = 0
      else
        if (p(1,1) == 0d0 .and. p(3,1) == 0d0) then
          ph = 0
        else
          ph = datan2(p(3,1),p(1,1))
        endif
        xprj = dcos(ph)*dsin(th)
        zprj = dsin(ph)*dsin(th)
      endif
C     print '(6f12.6)', p(1,1),p(2,1),p(3,1),th,xprj,zprj
      end
      subroutine pls3da(xp,yp,zp,np,syp,ndsyp,n1syp,type,syma,symcs,
     .  bsrad,bslist,ntp,tp,symat,nsymc)
C- Accumulates points symbol printing at projection of (x,y,z)
C ----------------------------------------------------------------------
Ci Inputs
Ci   See plsym for xp,yp,np,type,syma,symcs; zp is normal projection
Ci   syp,n1syp:if type<0, symbol parameters are passed in
Ci                syp(ip,n1syp)      symbol type
Ci                syp(ip,n1syp+1..3) color
Ci                syp(ip,n1syp+4..)  attributes
Ci   ndsyp    :dimensions syp
Ci   type     :symbol type:
Ci            : type=-1 => draw 3D symbol of type syp(:,n1syp)
Ci            : type<0  => same as type=-1
Ci            : type>0  => draw 2D symbol of this type (plsym)
Ci   bsrad    :radius limiting which pairs are to be connected
Ci   bslist   :list of allowed symbols to connect each point to.
Ci   ntp      :leading dimension of tp
Co Outputs
Co   tp(ntp,*) points accumulated into this table.  Cols:
Co             1 xp  2 yp  3 zp  4 ptr-symat
Co   symat     symbol attributes accumulated into this table
Co             1..10  symbol-specific attributes
Co                    (Get passed as syma in plsy3)
Co             11     symbol type
Co             12..14 colors
Co             15     bold
Co             16     bsrad
Co   nsymc     accumulates total number of symbol classes
Cr Remarks
Cr   After all symbols accumulated, follow with a call to pls3dp
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer type,np,ntp,nsymc,ndsyp,n1syp
      double precision syma(2),symcs(0:3),bsrad,xp(np),yp(np),zp(np),
     .  tp(ntp,4),symat(16,5),syp(ndsyp,8)
      character*40 bslist(0:30)
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .     plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .     plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .       cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .     xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .     linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .                plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .                plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .                     cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .                xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .                lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .                medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      logical err
      integer ipt,ipr,i1mach,i1,ibold,iprint,nsyma(7),ich,i
      double precision dum,dasum
      character*100 s,synam(7)*7
      save ipt
      data synam /'x','square','diamond','+','polygon','circle','arrow'/
      data nsyma /1,1,1,1,2,1,5/
      data ipt /0/

C --- Accumulate points (and symbol table if type<0) ---
      if (type == 0) return
      if (ipt+np > ntp) call rx('pls3da: tp too small')
      call getpr(ipr)
      call dcopy(np,xp,1,tp(1+ipt,1),1)
      call dcopy(np,yp,1,tp(1+ipt,2),1)
      call dcopy(np,zp,1,tp(1+ipt,3),1)
      if (type > 0) then
        call dcopy(np,dble(nsymc+1),0,tp(1+ipt,4),1)
        nsymc = nsymc+1
      else
        call plntyr(i1,ibold,dum,dum,dum,dum,err)
        ich = 0
        call chrpos(bslist,';',40,ich)
        ich = ich+2
        do  10  i = 1, np
          nsymc = nsymc+1
          tp(i+ipt,4) = nsymc
          call dcopy(10,syp(i,n1syp+4),ndsyp,symat(1,nsymc),1)
C         Symbol type: (-) sign indicates 3D symbol
          symat(11,nsymc) = -syp(i,n1syp)
C         Symbol color
          call dcopy(3,symcs(1),1,symat(12,nsymc),1)
          if (dasum(3,syp(i,n1syp+1),ndsyp) > 0)
     .      call dcopy(3,syp(i,n1syp+1),ndsyp,symat(12,nsymc),1)
C Symbol bold,bsrad
          symat(15,nsymc) = ibold
          symat(16,nsymc) = bsrad
          bslist(nsymc) = ' '
          if (ich <= 40) bslist(nsymc) = bslist(0)(ich:40)
   10   continue
        if (iprint() >= 40)
     .    call awrit1(' pls3da: symbols read from data file, '//
     .    '%i points.',s,len(s),i1mach(2),np)
      endif
      ipt = ipt+np
      if (type < 0) return

C --- Accumulate symbol table for type>0 ---
      call dcopy(10,syma,1,symat(1,nsymc),1)
      symat(11,nsymc) = type
      symat(12,nsymc) = symcs(1)
      symat(13,nsymc) = symcs(2)
      symat(14,nsymc) = symcs(3)
      call plntyr(i1,ibold,dum,dum,dum,dum,err)
      symat(15,nsymc) = ibold
      symat(16,nsymc) = bsrad
      ich = 0
      call chrpos(bslist,';',40,ich)
      ich = ich+2
      bslist(nsymc) = ' '
      if (ich <= 40) bslist(nsymc) = bslist(0)(ich:40)
      if (iprint() >= 40)
     .  call awrit3(' pls3da: '//synam(type)//' %i points  '//
     .  'gs=%3:1;5d  attr=%'//char(ichar('0')+nsyma(type))//':1;5d  ',
     .  s,len(s),i1mach(2),np,symcs(1),syma)
      end

      subroutine pls3dp(ntp,tp,symat,bslist,tabpr,iwk,rotm,shf)
C- Draw symbols in order of descending distance (yprj)
C ----------------------------------------------------------------------
Ci Inputs
Ci   tp(ntp,*):points for which to draw symbols tp(*,4) holds
Ci             index to symbol information to use in symat array.
Ci   symat    :symbol attributes
Ci             1..10  symbol-specific attributes
Ci                    (passed to plsym or plsy3)
Ci             11     symbol type
Ci                    type > 0 => 2D symbol (drawn by plsym)
Ci                                See plsym for each type and attributes
Ci                    type < 0 => 3D symbol
Ci                                See plsy3 for each type and attributes
Ci             12..14 colors
Ci             15     bold
Ci             16     bsrad
Ci   bslist   :list of allowed symbols to connect each point to.
Ci   tabpr    :work array, length 3*(max number of pairs)
Ci   iwk      :integer work array, length 2*ntp
Ci   rotm,shf :rotation and translation matrix, needed for 3d symbols
Co Outputs
Co   symbols are drawn
Cr Remarks
Cr   tp,symat generated by pls3da
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical err
      integer ntp,iwk(ntp,2)
      character*40 bslist(0:30)
      double precision tp(ntp,4),symat(16,1),tabpr(3,1),shf(3),rotm(3,3)
C Local parameters
      integer itype,ibold,iprint,i1mach,ip,ipp,ips,i1,ipr,npr,jp,nbslst
      double precision dum1(9),dum2,dum3,dum4,xprj,yprj,zprj,bsrad,
     .  cs(0:7)
      character*100 s

      call getpr(ipr)
      cs(0) = 3
      cs(4) = -1
      cs(5:7) = 0
      if (ntp == 0) call rx('problem with ntp')
      call dvshel(1,ntp,tp(1,2),iwk,1)
      call plntyr(i1,ibold,dum1,dum2,dum3,dum4,err)
      if (iprint() >= 30) call awrit1(' pls3dp: generating '//
     .  'symbols for %i points',s,len(s),i1mach(2),ntp)
      call pshpr(iprint()-10)
        if (ipr >= 50) print 333
  333   format(' indx',4x,'xrot',8x,'zrot',8x,'yrot',8x,
     .    'xprj',8x,'zprj',8x,'yprj',4x,'npr')
      do  10  ip = ntp, 1, -1
        ipp = iwk(ip,1)+1
        ips = nint(tp(ipp,4))
        itype = nint(symat(11,ips))
        ibold = nint(symat(15,ips))
        bsrad = symat(16,ips)
        if (bslist(ips) == ' ') then
          nbslst = 0
          iwk(1,2) = -1
        else
          call mkilst(bslist(ips),nbslst,iwk(1,2))
        endif
C ...   Draw lines connecting all pairs ...
        call pairs(ipp,ntp,tp,nbslst,iwk(1,2),symat,npr,tabpr)
        if (npr > 0) then
          call plntys
          call plntyp(1,3,0d0,0d0,0d0,0d0)
          do  20  jp = 1, npr
            call prj3dp(tp(ipp,1),tp(ipp,2),tp(ipp,3),.false.,
     .        dum1,dum1,dum1,dum1,dum1,xprj,yprj,zprj)
            call mve(xprj,zprj)
            call prj3dp(tabpr(1,jp),tabpr(2,jp),tabpr(3,jp),.false.,
     .        dum1,dum1,dum1,dum1,dum1,xprj,yprj,zprj)
            call drw(xprj,zprj)
   20     continue
          call plntyg
        endif
C ...   Draw this symbol
        call plntyp(i1,ibold,dum1,dum2,dum3,dum4)
        call prj3dp(tp(ipp,1),tp(ipp,2),tp(ipp,3),.false.,
     .    dum1,dum1,dum1,dum1,dum1,xprj,yprj,zprj)
        if (ipr >= 50) print '(i3,6f12.6,i4)', ips,
     .    tp(ipp,1),tp(ipp,3),tp(ipp,2),xprj,zprj,yprj,npr
        if (yprj < 0) goto 10
        if (itype > 0)  then
          symat(1,ips) = symat(1,ips)/yprj
          call dcopy(3,symat(12,ips),1,cs(1),1)
          call plsym(xprj,zprj,0d0,1,itype,symat(1,ips),cs)
          symat(1,ips) = symat(1,ips)*yprj
        else
C ...     Undo rotation!
          dum1(1) = tp(ipp,1) - shf(2)
          dum1(2) = tp(ipp,2) - shf(1)
          dum1(3) = tp(ipp,3) - shf(3)
          call dgemm('T','N',3,1,3,1d0,rotm,3,dum1,3,0d0,dum1(4),3)
          call dcopy(3,symat(12,ips),1,cs(1),1)
          call plsy3(dum1(4),dum1(5),dum1(6),rotm,shf,
     .      -itype,symat(1,ips),cs)
        endif
   10 continue
      call poppr
      end

      subroutine pairs(ip,ntp,tp,nlist,list,symat,npr,rtab)
C- Accumulate table of all pairs for which yprj<0 and dr < ri+rj
C  and for which the fourth element of tp belongs to list
Co rtab: table of connecting points
      implicit none
      integer ip,ntp,npr,nlist,list(nlist)
      double precision tp(ntp,4),symat(16,1),rtab(3,1)
      integer is,js,k,jp,ipr
      double precision rr,r1,r2,xprj,yprj,zprj,dum1(1)
      call getpr(ipr)
      npr = 1
      is = nint(tp(ip,4))
      r1 = symat(16,is)
      do  20  jp = 1, ntp
        if (jp == ip) goto 20
        js = nint(tp(jp,4))
        r2 = symat(16,js)
C ...   Reject unless js in list
        if (list(1) /= 0) then
          do  22  k = 1, nlist
   22     if (js == list(k)) goto 24
          goto 20
   24     continue
        endif
C ...   Add to list if connecting vector < r1+r2
        rtab(1,npr) = tp(jp,1)-tp(ip,1)
        rtab(2,npr) = tp(jp,2)-tp(ip,2)
        rtab(3,npr) = tp(jp,3)-tp(ip,3)
        rr = dsqrt(rtab(1,npr)**2 + rtab(2,npr)**2 + rtab(3,npr)**2)
C        if (rr < r1+r2 .and. rtab(2,npr) <= 0) then
        if (rr < r1+r2) then
          call prj3dp(rtab(1,npr),rtab(2,npr),rtab(3,npr),.false.,
     .      dum1,dum1,dum1,dum1,dum1,xprj,yprj,zprj)
          if (yprj < 0) then
            if (ipr > 50)
     .      print '('' pairs: ip,jp ='',2i3,''  vec ='',6f12.6)',
     .      ip,jp,rtab(1,npr),rtab(2,npr),rtab(3,npr)
            rtab(1,npr) = rtab(1,npr)+tp(ip,1)
            rtab(2,npr) = rtab(2,npr)+tp(ip,2)
            rtab(3,npr) = rtab(3,npr)+tp(ip,3)
            npr = npr+1
          endif
        endif
   20 continue
      npr = npr-1
      end
      logical function setuu(xmin,xmax,ymin,ymax,big,logx,logy,pad,plu,
     .  uuytox)
C- Set user's units
Co plu
      implicit none
      logical logx,logy
      integer j
      double precision big,xmin,xmax,ymin,ymax,plu(6),uuytox,pad
      double precision xx,yy

      setuu = .false.
      if (xmax == -big .or. ymax == -big) return

C ... Widen by padding
      xx = xmax-xmin
      yy = ymax-ymin
      if (logx) xx = dlog(xmax) - dlog(xmin)
      if (logy) yy = dlog(ymax) - dlog(ymin)
      plu(1) = xmin-pad*xx/2
      plu(2) = xmax+pad*xx/2
      plu(3) = ymin-pad*yy/2
      plu(4) = ymax+pad*yy/2
      if (logx) then
        plu(1) = xmin/dexp(pad*xx/2)
        plu(2) = xmax*dexp(pad*xx/2)
      endif
      if (logy) then
        plu(3) = ymin/exp(pad*yy/2)
        plu(4) = ymax*exp(pad*yy/2)
      endif
      plu(5) = xx
      plu(6) = yy

      call lodsyv('xmin',1,plu(1),j)
      call lodsyv('xmax',1,plu(2),j)
      call lodsyv('ymin',1,plu(3),j)
      call lodsyv('ymax',1,plu(4),j)
C     call shosyv(0,0,0,6)

      setuu = .true.
      end
      subroutine err(test,vsn,strn)
C- Error test for fplot
      implicit none
      integer test
      double precision vsn
      character*(*) strn

      if (test <= 0) return
      if (len(strn) > 0 .and. strn /= ' ') print *, strn
      if (test <= 1) goto 100
      print 300
      if (test <= 2) then
        call info2(0,1,0,' fplot version %,2d',vsn,0)
        goto 100
      endif
    1 format(2x,a,:':',T20,a:/T20,a)
    3 format(T20,a:/T20,a)
    4 format(2x,a,':':/T20,a:':'/T20,a)
    5 format(T20,a,':',T32,a:/T32,a:/T32,a)
    6 format(2x,a,': ',a)

C      print 301
C      print 302
C      print 310
C      print 311
C      print 320
C      print 321
C      print 330
  300 format(' usage: fplot [-INIT-switches] [-FORMAT-switches]',
     .  ' [-DATA-switches] data-file ... '/
     .  8x,'fplot -h will print command line options')

      print 1,'INIT switches (must occur first):'

      print 1,'-h ','print a short description of switches and exit'

      print 1,'-rot ','rotates frame 90 degrees'

      print 1,'-pr# ','set the output verbosity to #'

      print 1,'-plm=#,#,#,# ','set the medium plot boundaries'

      print 1,'-shftm=#,# ','shift the medium plot boundaries.',
     .  'Note: top right corner at (612,792)'

      print 1,'-disp[:lus] ','displays picture in a subshell'
      print "(t20,'Optional landscape, upside-down, seascape')"

C     print 1,'-plm=l,b,r,t ','Sets absolute medium units'

      print 1,'-plaintext ','strings have no super or subscripts'

      print 1,'-f fnam','read remaining arguments from file fnam'

      print "(/'  FORMAT switches (govern frame layout and labels)')"

      print 4,
     .  '-frme[:lx|:ly|:lxy][:xor=#][:yab=#][:nofill][:col=#,#,#,][:theta=#][:font=#] l,r,b,t',
     .  'starts a new frame, in box l,r,b,t(GU)'
      print 5,':lx','| :ly | :lxy for log scales in x,y, or both'
      print 5,':xor=#','draw vertical axis at x=#'
      print 5,':yab=#','draw horizontal axis at y=#'
      print 5,':nofill','suppresses filling box with background before drawing frame'
      print 5,':col=#,#,#','specify background color when filling frame'
      print 5,':theta=#','angle between abscissa and ordinate (radians)'
      print 5,':font=#','font for frame numbering and labels'

      print 4,'-x x1,x2 (or -y y1,y2)','specifies x (y) plot boundaries'


      print 4,'-tmx (-tmy) spacing[:mt][,pos][;rmt][~rmnt][@mode]',
     .  'tic mark specification and attributes.'
      print 5,'spacing ','spacing between tics (see mode)'
      print 5,'mt ','no. tics per major tic'
      print 5,'pos ','position and size of major tic'
      print 5,'rmt ','size of major tic'
      print 5,'rmnt ','size of minor tic, relative to major tic'
      print 5,'mode ','1,2,3 for log;  mode=5 => caller specifies tics,'
     .  ,'eg 0@5,1,4,6 puts tic marks at 1,4,6'

      print 4,'-frmt  [col=#,#,#,][th=#1[,#2,#3]]','frame parameters'
      print 5,'th=#1[,#2,#3]] ','line specifications'
      print 3,'  #1 line thickness (default=3; -1 makes no frame)'
      print 3,'  #2 one of 0,1,2,3 to draw bottom and top axis, ',
     .             '     or one only, or neither'
      print 3,'  #3 one of 0,1,2,3 to draw left and right axis, ',
     .             '     or one only, or neither'
      print 1,'-p# ','pads frame boundary (as fraction of frame)'
      print 1,'-1p ','skip parsing of remaining args for first pass'
      print 1,'-ndpi=# ',
     .  'specifies resolution of medium (dots per inch)'
      print 1,'-nx (-ny) ','normalizes each column of abscissas and'//
     .             '(ordinates) to (0,1)'

      print 1,'-aspect # ',
     .  'resizes frame to make dy/dx(UU) / dy/dx(MU) = #'
      print 4,'-3d shy[,shx,shz] ',
     .  'for 3-d perspective drawing.'
      print 3,'Points are shifted by shy (shx,shz).'
      print 3,'Projection is x to right, z up, y normal to paper.'

      print 4,'-rotp rot1[,rot2,...] ','rotates points (3d only).'
      print 3,'Example: (0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2'
      print 3,'Alternatively: z:pi/4,y:pi/3,z:pi/2'

      print "(/'  LABELLING and NUMBERING switches')"
      print 4, "  Where strings are referenced below, parts "//
     .  "inside {} may be treated specially"
      print 4, "  ^{..}, _{..}, ~{..}, @{..}, &{..}, "//
     .  "are turned into"
      print 1, "  superscripts, subscripts, Greek, bold, italic,"//
     .  " respectively."
      print 1, '  Example: "&\{k}_\{~\{\{\136}}}" makes italic k, followed by subscript perp symbol'
      print 1

      print 4,'-lbl[um] x,y[:blk] cc[,rot=#] str [tex-string]',
     .'label and attributes.'
      print 3,'puts `str'' at (x,y),  centered at cc, in'//
     .             ' user (medium) units.'
      print 3,'cc: one of chars l,c,r followed by one of '//
     .             'u,c,d, eg `ld'''
      print 3,"OR: cc is `tx,' which flags ps file will be run through"
      print 3,'latex with string substitutions.  In which case:'
      print 3,'  `tex-string'' is the actual (TeX) string and'//
     .        ' is required;'
      print 3,'  `str'' is merely a tag needed for latex processing'
      print 4,'-lblx xlst y[:blk] cc str1 str2 ...'//
     .             '   or -lbly ylst x[:blk] cc str1 str2 ...',
     .  'puts labels (strings) at list of x-points (y-points)'

      print 1,'-xl ',
     .  'str (-yl str) (-tl str) for x (y) axis or title labels'
      print 1,'-font ','t# (h#) changes Times (Helvetica) font and size'
      print 4,'-k x,y[:len][,spacing][;style]','specifies key '//
     .  'position (and length, spacing or style)'
      print 6,'-fmtnx:string (-fmtny:) ','awrite format for axis label'
      print 6,'-xn:0 | -noxn (-yn:0 | -noyn) ','no x (y) axis numbering'
      print 1,'-xn:t (-yn:r) ',
     .  'place x (y) axis numbering on top (right)'

      print "(/'  DATA switches (govern how data are drawn as lines)')"
      print 4,'-lt n[:bold=#][:col=#,#,#][:colw=#,#,#]'//
     .        '[:fill=#][:brk=#][:la[,lb,lc,ld]]',
     .  'line type specification and attributes.'
      print 5,'n','line type (0=none, 1=solid, 2=dashed, 3=dotted)'
      print 5,'bold','line thickness (use 0-9; default=3)'
      print 5,'col=#,#,#','line color (rgb)'
      print 5,'colw=#,#,#','secondary color when weighted by point'
      print 5,'fill=#','1 color line  2 fill with color  3 both 1 and 2'
      print 5,'brk=1','starts new line if x_i > x_i-1'
      print 5,'la,lb','(dashed lines) lengths of dash and space'
      print 5,'lc,ld','(dashed lines) lengths of second dash and space (dot-dashed lines)'
      print 4,'-s S[~col=#,#,#][~clip][~bold=#][~fill=#]~sym1[,sym2 ..] for symbol S',
     .        'symbol specification and attributes.'
      print 5,'S is one of',
     .        'x square diamond + polygon circle arrow errbar timelin hist row wiggle'
      print 5,'         or','an index 1-12'
      print 5,'         or','-1 to read symbol parameters from data'//
     .  ' file.','In that case, columns must hold:',
     .  '4: symbol type  (1=>arrow 2=>cone)  5-7: color  8-*, symbol attributes'
      print 3,'fill=#, #=  0 do not fill symbol; 1 fill with gray; 2 fill with color'
      print 5,'bold=#',
     .  'bold for symbol contour. 0 => no contour (just fill symbol)'
      print 5,'sym1,sym2,..','attributes that alter size and shape of symbol'
      print 1,'-l[0] legend for key (option 0 suppresses blanking)'
      print 1,'-tp [nc~]xlist','generates a table of points:'
      print 3,'nc=0  (default) maps points to 2 cols: xi=point, yi=0'
      print 3,'nc=1  puts points in a single column'
      print 3,'nc=n  orders points as dat(1,1:nc); dat(2,1:nc)...'
      print 3,'nc=-n orders points as dat(1:nr,1); dat(1:nr,2)...'
      print 4,'-map [-i expr] list',
     .  'maps original data array to another set of points'//
     .  ' defined by "list"'
      print 3,'-i includes only points satisfying logical expr'
      print 1,'-itrp x1,x2,dx','maps data to interpolated polynomial'
      print 3,'Optional args x1,x2,dx,ifrat,nord also specifies '
      print 3,'rational function interpolation and polynomial order'
      print 1,'-ins[f] string','insert string (or file named by string)'
     .      //' directly into output file'
      print 1,'-sort','sorts data in ascending order (after mapping)'
      print 4,'-bs radius[;file-list]',
     .  '(3D only) connects points within radius and file-list'

      print 1, '-con[:opts] list','special mode for x-y contour plots.'
      print 3, 'Contours of constant f=# are interpolated'//
     .  ' for each # in ''list''.'
      print 3, 'File data must correspond to function values f(xi,yi)'
      print 3, 'for a uniform mesh of points (xi,yi)'
      print 3, 'Vector in a fixed COLUMN tabulates line parallel '//
     .  'to abscissa, f(xi,y=const)'
      print 3, 'Vector in a fixed ROW tabulates line parallel '//
     .  'to ordinate, f(x=const,yi).'
C      print 3, 'Abscissa, f(xi,y=const): row varies, col fixed'
C      print 3, 'Ordinate, f(x=const,yi): row fixed, col varies'
      print 3, 'Bottom and top (left and right) edges of frame'
      print 3, 'correspond to first and last columns (rows) of data.'
      print 3, 'Ordinate f(x=const,yi): row fixed, col varies'
      print 3, 'Options:'
      print 3, ':dup=#1[,#2] Duplicate row (col) data #1 (#2) times'
      print 3, ':fn=name dumps xy pen moves to file "name"'
      print 3, ':noclose suppresses closing contours on boundaries'

C      print 1, '-nc=# (nr=#)','stipulate that next matrix read'
C     .             //' has # cols (rows)'
      print 1, '-r:switches','switches for file reading'
     .            //' separated by ''~''.  Switches are:'
      print 3, 'nr=#      stipulate that next data read has nr rows'
      print 3, 'nc=#      stipulate that next data read has nc columns'
      print 3, 'qr        read with fortran read (fast, no algebra; nr,nc must be stipulated)'
      print 3, 'br        read from binary file'
      print 3, 's=#       skips # records before reading'
      print 3, 'open      leaves file open after reading'
C      print 3, 'spc       store in sparse format'
C      print 3, 'br,nr,nc: file missing 1st record and nr,nc'//
C     .         ' are supplied'
      print 1,'-qr','same function as -r:qr'
      print 1,'-br','same function as -r:br'
      print 1, '-nr=#','same function as -r:nr=#'
      print 1, '-nc=#','same function as -r:nc=#'
      print 1,'-ab  expr','uses expression "expr" for the abscissa.  '//
     .        'Expression'
      print 3,'can contain variables x# where #'//
     .        ' signifies value of column j'
      print 1,'-abf expr','maps tic marks on abscissa to expr'
      print 1,'-ord expr','substitutes expression for the ordinate'
      print 1,'-col cx,cy[,cw]',
     .        'plots cols cx and cy as abscissa and ordinate'
      print 1,'-colsy list','makes series of plots for columns in list'
      print 3,'If the data consists of a single column, it is copied to column 2'
      print 3,'Thus "fplot -colsy 2 ..." plots the column with row index as abscissa'
      print 1,'-colsw list','corresponding list of columns'//
     .        ' for (color) weights'
      print 1,'-ey #[,dx,shfy]','error bars in y using data from col #'
      if (test == 3) return ! just print usage and do not return error.
  100 call cexit(-1,1)
      end

      subroutine parsefont(fontstr,fontname,fontsz,lerr)
C- Parse fontstr to return Postscript font string and size
C lerr :  On input if .true. => abort with error if illegal font
C         otherwise if illegal font return with fontsz=-1
      implicit none
      integer fontsz
      character(len=*) :: fontstr,fontname
      logical lerr

      logical a2bin
      integer j,k
      integer, parameter :: nfont=5,strnsiz=19
      character(len=nfont), parameter :: fontlist='thibs'
      character(len=nfont*strnsiz), parameter ::
     .  fontnamtab =
     . '/Times-Roman       /Helvetica         /Times-Italic      /Times-Bold        /Symbol'

      k = 0
      call chrps2(fontstr,fontlist,nfont,1,k,j)
      if (j == 0) goto 999
      fontname = fontnamtab(1+(j-1)*19:j*19)
      k = 1
      if (.not. a2bin(fontstr,fontsz,2,0,' ',k,-1)) goto 999
      return

  999 continue
      if (lerr) call rx('parsefont : font not recognized : '//trim(fontstr))
      call info0(1,0,0,'warning!  font not recognized : '//trim(fontstr))
      fontsz = -1
      end

      subroutine frmeprms(mode,fillcs,xor,yab,frmefont)
C- Set/retrieve parameters for frame
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 0 do nothing
Ci         : 1 set current values to default settings
Ci         : 2 push current values to stack
Ci         : 3 set current values and stack to default settings
Ci         : 4 exchange current values with stack
Ci         : 5 Do (1), then (4)
Ci         : 6 copy stack into current values
Ci   fillcs
Ci   xor
Ci   yab
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   13 Oct 16
C ----------------------------------------------------------------------
      implicit none
      integer mode
      double precision fillcs(0:4),xor,yab
      character*(*) frmefont

      real(8), parameter :: NULLR = -99999
      real(8) :: stackp(10),stackn(10),stackd(10),frmefontd
      save stackp

C ... Default values
      stackd(1:9) = [103d0,1d0,1d0,1d0,1d0,   NULLR, NULLR, NULLR, NULLR]
      call s8tor8(' ',stackd(10))
      call s8tor8(frmefont,frmefontd)

C ... Push current values onto temporary stack (modes 2,4,5) ... this will become new stack
      if (mode == 2 .or. mode == 4) then
        stackn(1:5) = fillcs(0:4); stackn(8) = xor; stackn(9) = yab; stackn(10) = frmefontd
      endif
      if (mode == 3 .or. mode == 5) stackn = stackd  ! New stack for modes 3,5

C ... Pop stack or defaults into current values: use stackd as nondestructive staging for copy
      if (mode >= 4) stackd = stackp
      if (mode >= 3) then
        fillcs(0:4) = stackd(1:5); xor = stackd(8); yab = stackd(9); frmefontd = stackd(10)
      endif

C ... Copy temporary stack to stack
      if (mode >= 2 .and. mode < 6) then
        stackp = stackn
      endif
      call r8tos8(frmefontd,frmefont)

      end

      subroutine frmelabel(mt,ntx,nty,ntz,ffont,ltfrme,modxy,xnum,ynum,znum,
     .  ts,rmt,rmnt,tcp,xnfmt,ynfmt,znfmt,xlabel,ylabel,zlabel,title)
C- Sets default values for frame label
      implicit none

      integer mt(3),ntx,nty,ntz,ffont(2),ltfrme(3),modxy(2),xnum,ynum,znum
      double precision ts(50,3),rmt(3),rmnt(3),tcp(3)
      character*(*) xnfmt,ynfmt,znfmt,xlabel,ylabel,zlabel,title

      tcp(1) = -1d30; tcp(2) = tcp(1); tcp(3) = tcp(1)
      call dpzero(ts,150); ts(1,1:3) = [-1d0,-1d0,-1d0]
      mt(1:3) = -1
      rmt(1:3) = [.025d0,.025d0,.025d0]
      rmnt(1:3)= [.6d0,.6d0,.6d0]
      ltfrme(1) = 3; ltfrme(2) = 0; ltfrme(3) = 0
      modxy(1:2) = -1
      xnum = 1; ynum = 1; znum = 1
      xnfmt = '%;4d'; ynfmt = '%;4d'; znfmt = '%;4d'
      xlabel = ' '; ylabel = ' '; zlabel = ' '; title  = ' '
      ntx = -1; nty = -1; ntz = -1

      end

