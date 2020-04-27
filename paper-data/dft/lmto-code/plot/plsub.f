C#define INCH_UNITS
C --- Fortran primitives for plot procedures. ---
Cr This is a collection of routines for drawing lines and symbols.
Cr The hardware dependence is mostly concentrated into an initialization
Cr procedure (pltini) and one (pltam) that draws straight lines, to the
Cr extent it was feasible to do so.  The core routines are:
Cr   pltini(character device, integer opt)
Cr     defines and possibly initializes the current plot device
Cr   pltam(real x, real y, integer opt)
Cr     Moves the pen from the current position to point x-y on the
Cr     current device, in the units appropriate to the plot medium.
Cr     opt=1 plots with pen down; opt=0 plots with pen up.
Cr     All pen moves are done by calling pltam.
Cr     (These units are set in pltini).
Cr   pltdmp(io)
Cr     Dumps the current plot to the output device.  This is hardware-
Cr     dependent.  When it is applicable and io=1, pltdmp prompts for
Cr     a carriage return before executing.
Cr   plsym(xp,yp,zp,np,type,syma,symcs)
Cr     draws a symbol at a collection of points (xp,yp)
Cr
Cr Besides the (hardware-dependent) medium units, there are two sets of
Cr device-independent units.  By default, graphics units range from (0,0)
Cr to (1,1) for the largest square that fits within the plot medium.
Cr User's units define the user's values of x and y at the corners
Cr of the graph, and obviously will depend on the specific plot.
Cr There are several routines for defining the various units
Cr   pltstu(l,r,b,t):
Cr     defines the corners of a graph, in user's units
Cr   pltstp(l,r,b,t):
Cr     restricts (clips) a graph to a portion of the medium's surface
Cr     (in graphics units).  no check to see within (l,r,b,t)
Cr   pltstg(l,r,b,t):
Cr     Change the corners of the medium (in graphics units) from their
Cr     default values (0,0) and (1,1).  Not normally needed.
Cr   plthet(alpha):
Cr     Shears the orthogonal (x,y) axes by the linear transformation:
Cr       x->x
Cr       y->y*sin(alpha) + x*cos(alpha)
Cr      Default alpha is pi/2 (orthogonal x,y)
Cr The following procedures are for drawing straight lines.
C  Arguments (x,y) in user's units, except when stated:
Cr   mve(x,y)  moves the pen to (x,y) from current posn, pen up
Cr   drw(x,y)  draws the pen to (x,y) from current posn, pen down
Cr   mver(x,y) like mve, but motion relative to current x,y
Cr   drwr(x,y) like drw, but motion relative to current x,y
Cr   pposu(x,y) returns the current pen position in (x,y).
Cr   drwe(x,y,cs) draws the last point in a curve.  Area enclosed by
Cr             curve may be be filled in, according to cs.
Cr   drwre(x,y,cs) like drwe, but motion relative to current x,y
Cr   mveg(x,y) like mve, but x,y in GU
Cr   mvegr(x,y) like mveg, but motion relative to current x,y
Cr   pposg(x,y) like pposu, but in GU
Cr
Cr Other routines:
Cr   plcrv  draws a smooth curve for a collection of points.
Cr          For postscript, it writes the file directly instead of
Cr          calling the generic pltam, to preserve continuity of
Cr          dotted, dashed lines, and some other features.
Cr   plcrve closes a curve, possibly shading area enclosed by curve
Cr   frme,ticset,ticnxt draw a frame and tic marks.
Cr   pltsts sets current plot to linear, log scale
Cr   pstr   puts a string at current x,y
Cr   setfnt sets the font for drawing strings.
Cr   plntyp inputs variables to make dashed, broken or dotted lines.
Cr   pltcol inputs variables for color (not implemented)
Cr   plndpi sets the medium resolution (dots per inch)
Cr   key    draws a key
Cr   plstmv sets or reads the plot medium boundaries as symbolic
Cr          variables.  This permits user to redefine, e.g. the edges
Cr          of the page.
Cr
Cr Some Postscript-specific information and routines are:
Cr   psinit writes to the file stream some postscript macros
Cr   pslabl writes a string
Cr For postscript, the "edges" of the default to something smaller than
Cr a 8+1/2 x 11 inch page.  plstmv permits you to override the internal
Cr defaults (now set at (l,r,t,b) = (144,612-72,144,792-72) points).
Cr
Cr The initialization routine does the following:
Cr   1. sets the corners of the medium, in MU
Cr   2. defines scalxc,scalyc converting from MU into centimeters
Cr   3. defines xb,yb, the "line thickness" (for devices where
Cr      this is not adjustable, e.g. Hercules graphics card)
Cr   4. other device-specific initialization (eg calls psinit)
Cr   5. Defines the default graphics units and sets the clip region
Cr      to the largest available square.
Cr
Cv Common block variables:
Cv The variables corresponding to these units are:
Cv   plgl,r,b,t: graphics units demarcating value of x and y at four
Cv               corners of the plotting medium.  They default to (0,0)
Cv               and (1,1) but can be set explicitly.
Cv   plpl,r,b,t: values of x and y, in graphics units, at four corners
Cv               of a plot (useful when a plot is clipped to a portion
Cv               of the entire medium).
Cv   plul,r,b,t: user's units demarcating value of x and y at four
Cv               corners of a plot.
Cv   plml,r,b,t: units of the actual plot medium
Cv               (set automatically by pltini)
Cv   scalxg,scalyg,offxg,offyg (calculated in pltstg), make the
Cv               linear transformation between medium units
Cv               and graphics units:
Cv               xm = (offxg+xg)*scalxg, ym = (offyg+yg)*scalyg
Cv   scalxu,scalyu,offxu,offyu (calculated in pltstu), make the
Cv               linear transformation between medium units
Cv               and user units:
Cv               xm = (offxu+xu)*scalxu, ym = (offyu+yu)*scalyu
Cv   xcm,ycm:    current pen position in medium's units
Cv   xcu,ycu:    current pen position in user's units
Cv   xcc,ycc:    not used
Cv   strtch      the ratio of MU unit lengths in x,y that produces a
Cr               true square as seen by the viewer.  Unit lengths in x,y
Cr               may have differing true dimensions, as they probably
Cr               will when units refer to pixels.
Cv   stream:     file logical unit for plot output
Cr   xb,yb       the "line thickness" in MU (for devices where
Cr               this is not adjustable)
Cr   lintyp,ltbld define current
Cr               line type and boldness.  lintyp:
Cr               0 for none, 1 for solid, 2 for dashed, 3 for dotted
Cr   linl1,linl2,linl3,linl4, sequence of lengths for dashed lines
Cr               (dash, space, dash, space)
Cv   lnprd       scaling of line segments to medium units.  Defined
Cv               so that "typically sized" line segments have size of
Cv               order 1, independently of the medium.
Cv   noclip      true if no clipping has been set
Cv   medium      index defining current medium
Cv   logx,logy   true when user's units refer to a logarithmic scale
Cv   pstat       current pen status: 1 is down, 0 is up
Cv   ndpi,ndwrt  resolution, in dots per inch. number of digits
Cv               needed to reach medium's resolution
C ----------------------------------------------------------------
      subroutine pltini(device,opt,lrotp)
C- Initialization for plot
C ----------------------------------------------------------------
Ci Inputs
Ci   device: device name
Ci   opt:    see Remarks
Ci   lrot:   rotate figure
Co Outputs
Co   Initializes plot medium
Cr Remarks
Cr   pltini can be invoked using a number (opt) or (if 0) a device name.
Cv Variables
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) device
      integer opt
      logical lrotp

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
      integer nnames,idum,itx(2)
      procedure(logical) :: cmdopt
      procedure(integer) :: fopna,fopng,i1mach,isw,a2vec
      double precision transx,transy,scale,xx(2)
      character prmed*72
      parameter (nnames=3)
      character*10 names(nnames)
      data names /'ascii','hercules','ps'/

      lrot = lrotp

C --- If no number specified, seek one corresponding to name ---
      medium = opt
      if (medium == 0) then
   10   medium = medium+1
        if (names(medium) == device .or. medium >= nnames) goto 20
        goto 10
      endif

   20 goto (110,120,130), medium
      call rx('PLTINI: device not present')

C --- Ascii: write to file unit 'stream' ---
  110 continue
c      stream = i1mach(2)
      if (stream /= i1mach(2)) idum = fopna('wplt',stream,0)
C#ifdef INCH_UNITS
C --- Set medium's boundaries  ---
      lnprd = .1d0
      plml = 1.5d0
      plmr = 7.5d0
      plmb = 0d0
      plmt = 9d0
C --- Define size of cm units in terms of medium units ---
      scalxc = 1/2.54d0
      scalyc = 1/2.54d0
      strtch = dabs(scalxc/scalyc)
      offxc = 0d0
      offyc = 0d0
      xb = 1d0/150d0
      yb = xb
      goto 30
C#endif INCH_UNITS
C#ifdefC CM_UNITS
CC --- Set medium's boundaries  ---
C      lnprd = .27
C      plml = 0.
C      plmr = 24.
C      plmb = 0.
C      plmt = 18.
CC --- Define size of cm units in terms of medium units ---
C      scalxc = 1.
C      scalyc = 1.
C      strtch = abs(scalxc/scalyc)
C      offxc = 0.
C      offyc = 0.
C      goto 30
C#endif INCH_UNITS

C --- Hercules graphics card ---
  120 continue
      lnprd = 8
C --- Set medium's boundaries  ---
      plml = 0d0
      plmr = 719d0
      plmb = 317d0
      plmt = 0d0
C --- Define size of cm units in terms of medium units ---
C Calculate scaling by matching physical dimensions of the screen
C and putting (0,0) in cm units at (plml,plmb)
      scalxc = (plmr-plml)/19.5d0
      scalyc = (plmt-plmb)/14.0d0
      strtch = dabs(scalxc/scalyc)
      offxc = plml/scalxc - 0d0
      offyc = plmb/scalyc - 0d0
      goto 30

C --- PostScript ---
  130 continue
C     idum = fopna('ps',stream,0)
      idum = fopng('fplot.ps',stream,0)
C --- Set medium's boundaries  ---
      lnprd = 7.2d0
      plml = 144
      plmr = 612-72
      plmb = 144
      plmt = 792-72
      call plstmv(plml,plmr,plmb,plmt,' ',.true.)
      idum = 7
      if (cmdopt('-shftm=',idum,0,prmed)) then
        idum = a2vec(prmed,len(prmed),idum,4,', ',2,2,2,itx,xx)
        if (idum /= 2) call rx('failed to parse 2 arguments following -shftm')
        plml = plml + xx(1)
        plmr = plmr + xx(1)
        plmb = plmb + xx(2)
        plmt = plmt + xx(2)
      endif

      bbox(1) = 9999
      bbox(2) = -9999
      bbox(3) = 9999
      bbox(4) = -9999

C --- Define size of cm units in terms of medium units ---
      scalxc = 1d0/(2.54d0*72d0)
      scalyc = 1d0/(2.54d0*72d0)
      strtch = dabs(scalxc/scalyc)
      offxc = 0d0
      offyc = 0d0
      xb = 1d0/150d0
      yb = xb
      transx = 0
      transy = 0
      scale = 1
      call psinit(stream,isw(lrot),transx,transy,scale)
      goto 30

C --- Setup common to all devices ---
   30 continue
C      if (lrot) then
C        tmp  = plmt
C        plmt = plmr
C        plmr = plmb
C        plmb = plml
C        plml = tmp
C        tmp    = scalxc
C        scalxc = scalyc
C        scalyc = tmp
C        strtch = dabs(scalxc/scalyc)
C        tmp   = offxc
C        offxc = offyc
C        offyc = tmp
C        tmp = xb
C        xb = yb
C        yb = tmp
C      endif

C --- Define GU ---
      if (dabs((plmt-plmb)*strtch) > dabs(plmr-plml)) then
        call pltstg(0d0,1d0,0d0,dabs((plmt-plmb)*strtch/(plmr-plml)))
      else
        call pltstg(0d0,dabs((plmr-plml)/(plmt-plmb)/strtch),0d0,1d0)
      endif

C --- Clip to largest available centered square ---
      plpl = (plgr+plgl)/2 - .5d0
      plpr = (plgr+plgl)/2 + .5d0
      plpb = (plgt+plgb)/2 - .5d0
      plpt = (plgt+plgb)/2 + .5d0
      end
      subroutine plstmv(plml,plmr,plmb,plmt,s,lread)
C- Set or read plot medium boundaries from as symbolic variables
Cr Routine is initially called from main with lread=.false.
Cr If string s has a form =#1,#2,#3,#4
Cr then symbolic variables are set as
Cr plml = #1  plmb = #2   plmr = #3   plmt = #4
Cr Routine is later called from pltini with lread=.true.
Cr plml,plmb,plmr,plmt have been set to default values.
Cr Symbolic variables plml, plmb, plmr, plmt overwrite defaults
      implicit none
      logical lread
      character *(*) s
      double precision plml,plmr,plmb,plmt
      character*4 sp,st,pnam
      logical a2bin,lxst
      integer i,n,ich,ls,it,ich0,iprint,i1mach
      double precision xx(4)
      data sp /'lbrt'/ st /', '/

      if (lread) then
        lxst = .false.
        do  i = 1, 4
          call getsyv('plm'//sp(i:i),xx(i),n)
C         call shosyv(0,0,0,i1mach(2))
          if (n /= 0) then
            lxst = .true.
            if (i == 1) plml = xx(i)
            if (i == 2) plmb = xx(i)
            if (i == 3) plmr = xx(i)
            if (i == 4) plmt = xx(i)
          endif
        enddo
        if (lxst .and. iprint() >= 0) call awrit4(
     .      ' plstmv:  plml=%d  plmb=%d  plmr=%d  plmt=%d',
     .      ' ',80,i1mach(2),plml,plmb,plmr,plmt)
      else
        ls = len(s)
C        print *, s
        ich = 0
        call chrpos(s,'=',ls,ich)
        if (ich >= ls) goto 22
        ich = ich+1
        ich0 = ich
        do  i = 1, 4
          ich = ich0
          call chrps2(s,st,2,ls,ich,it)
C          print *, ich0, ich, s(ich0+1:ich+1)
          if (ich > ich0) then
            if (a2bin(s,xx(i),4,0,st(it:it),ich0,-1)) then
              pnam = 'plm'//sp(i:i)
              call lodsyv(pnam,1,xx(i),n)
C             call lodsyv('plm'//sp(i:i),1,xx(i),n)
              if (st(it:it) == ' ') goto 22
            endif
          else
            if (s(ich0:ich0) == ' ') goto 22
            ich0 = ich0+1
          endif
        enddo
   22   continue
      endif
      end
      subroutine pltstg(l,r,b,t)
C- set coordinates of graphics units
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t

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
      integer iprint,i1mach

      plgl = l
      plgr = r
      plgb = b
      plgt = t
      lnused = 0

C --- Set the scale and offsets ---
C Calculate scaling by taking the smaller of the x and y scales
C and offset by centering both graphics units in medium units.
      scalxg = dabs((plmr-plml)/(plgr-plgl)*strtch)
      scalyg = dabs((plmt-plmb)/(plgt-plgb))
      if (scalxg <= scalyg) then
        scalxg = dsign(scalxg,plmr-plml)
        scalyg = dsign(scalxg/strtch,plmt-plmb)
      else
        scalxg = dsign(scalyg*strtch,plmr-plml)
        scalyg = dsign(scalyg,plmt-plmb)
      endif

      offxg = ((plml+plmr)/scalxg - (plgl+plgr))/2
      offyg = ((plmb+plmt)/scalyg - (plgb+plgt))/2

      if (iprint() >= 50)
     .call awrit4(' pltstg: medium corners (GU): lb = (%d %d)'//
     .  '  rt = (%d %d)',' ',80,i1mach(2),plgl,plgb,plgr,plgt)
      if (iprint() > 100)
     .  call awrit4('         scal(x,y)g = %d %d  off(x,y)g = '//
     .  '(%d %d)',' ',80,i1mach(2),scalxg,scalyg,offxg,offyg)
      if (iprint() > 100)
     .  call awrit4(' pltstg: medium corners (MU): lb = (%d %d)'//
     .  '  rt = (%d %d)',' ',80,i1mach(2),plml,plmb,plmr,plmt)

      end
      subroutine pltstu(l,r,b,t,uuytox)
C- set coordinates of user's units
C ----------------------------------------------------------------
Ci Inputs
Ci   l,r,b,t: frame corners (UU)
Ci   uuytox: if zero, does nothing.  Otherwise relative lengths in UU
Ci     are linked to those in MU as:  dy/dx (UU) = uuytox * dy/dx (MU).
Ci     This is accomplished by resizing the current frame to the largest
Ci     size permitted by the above constraint.
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t,uuytox

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
      integer iprint,i1mach,is,scrwid
      double precision xx1,xx2
      logical error
      parameter (scrwid=100)
      character*(scrwid) s

      call lodsyv('ul',1,l,is)
      call lodsyv('ur',1,r,is)
      call lodsyv('ub',1,b,is)
      call lodsyv('ut',1,t,is)
      plul = l
      plur = r
      plub = b
      plut = t
      error = .true.
      if (logx) then
        if (l <= 0 .or. r <= 0) goto 99
        plul = dlog(l)
        plur = dlog(r)
      endif
      if (logy) then
        if (b <= 0 .or. t <= 0) goto 99
        plub = dlog(b)
        plut = dlog(t)
      endif
      error = .false.

C --- Set initial frame to max size if uuytox spec'd and noclip ---
      if (uuytox /= 0 .and. noclip) then
        if (iprint() >= 50) then
          call awrit1(' pltstu: resizing frame in x to satisfy '//
     .   'dy/dx(UU) = %1;4g*dy/dx(MU).',' ',80,i1mach(2),uuytox)
          print '(10x,''Frame initially set at ...'')'
        endif
        if (.not. noclip) call pltstp(plpl,plpr,plpb,plpt)
        if (noclip) call pltstp(plgl,plgr,plgb,plgt)
        if (iprint() >= 50) print '(10x,''and reset to ...'')'
      endif

C --- Set the scale and offsets in this frame ---
C Calculate by equating points at corners, eg
C scalxu (plul + offxu) = scalxg (plpl + offxg)
      scalxu = scalxg*(plpr-plpl)/(plur-plul)
      scalyu = scalyg*(plpt-plpb)/(plut-plub)
      offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (plul+plur)/2)
      offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (plut+plub)/2)
c      if (lrot) then
c        scalxu = scalxg*(plpr-plpl)/(t-b)
c        scalyu = scalyg*(plpt-plpb)/(l-r)
c        offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (t+b)/2)
c        offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (l+r)/2)
c      endif

C --- Clip frame to largest permitted by uuytox ---
      if (uuytox /= 0) then
        if (scalxu > scalyu*strtch*uuytox) then
          scalxu = scalyu*uuytox
          xx1 = (plur-plul)*scalxu/scalxg/2
          xx2 = (plpr+plpl)/2
          call pltstp(xx2-xx1,xx2+xx1,plpb,plpt)
        else
          scalyu = scalxu/uuytox
          xx1 = (plut-plub)*scalyu/scalyg/2
          xx2 = (plpt+plpb)/2
          call pltstp(plpl,plpr,xx2-xx1,xx2+xx1)
        endif
        scalxu = scalxg*(plpr-plpl)/(plur-plul)
        scalyu = scalyg*(plpt-plpb)/(plut-plub)
        offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (plul+plur)/2)
        offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (plut+plub)/2)
      endif

   99 continue
      if (iprint() >= 50 .or. error) then
        s = ' '
        is = 0
        call bin2a(' pltstu: frame corners (',1,0,0,1,0,scrwid,s,is)
        if (logx) call bin2a('logx',1,0,0,1,0,scrwid,s,is)
        if (logy) call bin2a('logy',1,0,0,1,0,scrwid,s,is)
        call bin2a('UU ): lb = (%1;4g %1;4g)',1,0,0,1,0,scrwid,s,is)
        call bin2a('  rt = (%1;4g %1;4g)',1,0,0,1,0,scrwid,s,is)
        call awrit4(s(1:is),' ',scrwid,i1mach(2),l,b,r,t)
        if (error) call rx(
     .  'pltstu (abort): log scale incompatible with negative limit')
        if (iprint() >= 100)
     .  call awrit4(' pltstu: scal(x,y)u = %1;4g %1;4g  off(x,y)u = '//
     .  '(%1;4g %1;4g)',' ',scrwid,i1mach(2),scalxu,scalyu,offxu,offyu)
      endif

      if (iprint() >= 100) call shosyv(0,0,0,i1mach(2))

      end
      subroutine pltstp(l,r,b,t)
C- Bound area of current plot, in GU; redefine UU at these corners
C ----------------------------------------------------------------
Ci Inputs
Ci   l,r,b,t: left, right, bottom, top
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t
      integer i

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
      integer iprint,i1mach

      call lodsyv('gl',1,l,i)
      call lodsyv('gr',1,r,i)
      call lodsyv('gb',1,b,i)
      call lodsyv('gt',1,t,i)
      call lodsyv('ml',1,(offxg+l)*scalxg,i)
      call lodsyv('mr',1,(offxg+r)*scalxg,i)
      call lodsyv('mb',1,(offyg+b)*scalyg,i)
      call lodsyv('mt',1,(offyg+t)*scalyg,i)
      plpl = l
      plpr = r
      plpb = b
      plpt = t
      noclip = .false.
      if (iprint() >= 30) then
        call awrit4(' pltstp:  frame corners (GU): lb = (%d %d)'//
     .  '  rt = (%d %d)',' ',80,i1mach(2),plpl,plpb,plpr,plpt)
        if (medium == 3)
     .  call awrit4('%% pltstp: frame corners (GU): lb= (%d %d)'//
     .  '  rt = (%d %d)',' ',80,stream,plpl,plpb,plpr,plpt)
      endif

C ... Update plot bounds in user units
C     call pltstu(plul,plur,plub,plut)

      end
      subroutine pltsts(plogx,plogy,io)
C- Set scale (linear or log)
C ----------------------------------------------------------------
Ci Inputs
Ci   plogx,plogy,io
Cr Remarks
Cr   Sets logx,logy from ilogx or ilogy (io=1=> set ilog* from log*)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer io
      logical plogx,plogy

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
      integer iprint,i1mach

      if (io == 0) then
        logx = plogx
        logy = plogy
        if (iprint() >= 30) then
          call awrit2(' pltsts:  set  logx=%l  logy=%l',
     .      ' ',80,i1mach(2),logx,logy)
          if (medium == 3) call awrit2('%% pltsts:  set  logx=%l'//
     .      '  logy=%l',' ',80,stream,logx,logy)
        endif
      else
        plogx = logx
        plogy = logy
      endif

      end
      subroutine plthet(alpha)
C- Set angle for nonorthogonal shear of Cartesian axes
C ----------------------------------------------------------------
Ci Inputs
Ci   alpha : shearing angle for linear transformation of (x,y)
Cr Remarks
Cr   By default, the axes are orthogonal; their angle is pi/2.
Cr   For other values of alpha, (x,y) are sheared by
Cr   the linear transformation:
Cr       x->x
Cr       y->y*sin(alpha) + x*cos(alpha)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision alpha

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


      cosa = cos(alpha)
      sina = sin(alpha)

      end
      subroutine pltam(x,y,opt,cs)
C- Plot absolute position, medium's units
C ----------------------------------------------------------------
Ci Inputs
Ci   x,y
Ci   opt: -1, return the current pen position
Ci         0, lift pen before move
Ci         1, drop pen before move
Ci         2, drop pen before move, lift after
Ci   cs:   passed to plcrve (for PostScript)
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y,cs(0:3)
      integer opt

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

C Local variables
      integer iprint,is
      double precision d,ds,xs,ys,cost,sint,d1mach
      logical tog
      character s*100

      is = 0
      s = ' '
      if (opt == -1) goto 5
      if (opt >= 1) goto 2
      goto 3

C --- Draw all line types, Postscript ---
    2 continue
      if (medium == 3) then
        if (pstat == 0) then
          call pbbox(xcm,ycm)
          call awrit4('newpath %1;nd %1;nd moveto',s,100,0,
     .      ndwrt,xcm,ndwrt,ycm)
C          call bin2a('newpath',0,0,0,1,0,80,s,is)
C          call bin2a('(f8.1)',1,0,xcm,4,0,100,s,is)
C          call bin2a('(f8.1)',1,0,ycm,4,0,100,s,is)
C          call bin2a('moveto',1,0,0,1,0,80,s,is)
        endif
C        call bin2a('(f8.1)',1,0,x,4,0,100,s,is)
C        call bin2a('(f8.1)',1,0,y,4,0,100,s,is)
C        call bin2a('lineto',1,0,0,1,0,80,s,is)
        call pbbox(x,y)
        call awrit4('%a %1;nd %1;nd lto',s,100,0,ndwrt,x,ndwrt,y)
        call skpblb(s,100,is)
        is = is+1
        pstat = 1
        goto 5
      endif
      pstat = 1

      goto (10,20) lintyp
      goto 3

C --- Line type 1 Draw, ASCII  ---
   10 continue
      if (medium == 1) then
        write(stream,110) x, y
        d = dsqrt((x-xcm)**2 + (y-ycm)**2) + d1mach(1)
        cost = (x-xcm)/d
        sint = (y-ycm)/d
        if (ltbld >= 1) then
          write(stream,310) xcm-xb*sint/1.75d0, ycm+yb*cost/1.75d0
          write(stream,110) x-xb*sint/1.75d0, y+yb*cost/1.75d0
          write(stream,310) xcm+xb*sint/1.75d0, ycm-yb*cost/1.75d0
          write(stream,110) x+xb*sint/1.75d0, y-yb*cost/1.75d0
          if (ltbld == 1) write(stream,310) x, y
        elseif (ltbld >= 2) then
          write(stream,310) xcm-xb*sint, ycm+yb*cost
          write(stream,110) x-xb*sint, y+yb*cost
          write(stream,310) xcm+xb*sint, ycm-yb*cost
          write(stream,110) x+xb*sint, y-yb*cost
          write(stream,310) x, y
        endif
  110 format('PA ',2f14.5)
      else
        call rx('PLTAM: unknown medium')
      endif
      goto 5

C --- Line type 2 draw, ASCII ---
   20 continue
      d = dsqrt((x-xcm)**2 + (y-ycm)**2) + d1mach(1)
      cost = (x-xcm)/d
      sint = (y-ycm)/d
      tog = .true.
C Move pen up when tog false, down when tog true
   25 if (tog) then
        ds = min(max(0d0,linl1-linl2-lnused),d)
      else
        ds = min(max(0d0,linl1-lnused),d)
      endif
      xs = ds*cost + xcm
      ys = ds*sint + ycm
C ... Pen down motion
      if (tog .and. ds > 0) then
        if (medium == 1) then
          write(stream,110) xs, ys
          if (ltbld >= 1) then
            write(stream,310) xcm-xb*sint/1.75d0, ycm+yb*cost/1.75d0
            write(stream,110) xs-xb*sint/1.75d0, ys+yb*cost/1.75d0
            write(stream,310) xcm+xb*sint/1.75d0, ycm-yb*cost/1.75d0
            write(stream,110) xs+xb*sint/1.75d0, ys-yb*cost/1.75d0
            write(stream,310) xs, ys
          elseif (ltbld >= 2) then
            write(stream,310) xcm-xb*sint, ycm+yb*cost
            write(stream,110) xs-xb*sint, ys+yb*cost
            write(stream,310) xcm+xb*sint, ycm-yb*cost
            write(stream,110) xs+xb*sint, ys-yb*cost
            write(stream,310) xs, ys
          endif
        else
          call rx('PLTAM: unknown medium')
        endif
C ... Pen up motion
      else
        if (medium == 1) then
          write(stream,310) xs, ys
        endif
      endif
      xcm = xs
      ycm = ys
      d = d - ds
      lnused = mod(lnused+ds,linl1)
      tog = .not. tog
      if (d > 0) goto 25
      goto 5

C --- Move (pen up) ---
    3 continue
C ... Close previously opened curve
      call plcrve(medium,pstat,[-2d0,0d0,0d0,0d0],s,is)
      lnused = 0
      if (medium == 1) then
        write(s,310) x, y
  310   format('MA ',2f14.5)
        is = 31
      elseif (medium == 3) then
      else
        call rx('PLTAM: unknown medium')
      endif
      goto 5

C --- Exit ---
    5 continue
      if (opt == -1) then
        x = xcm
        y = ycm
      else
        xcm = x
        ycm = y
        if (opt == 2) call plcrve(medium,pstat,cs,s,is)
        if (is > 0) write(stream,'(a)') s(1:is)
      endif

    6 continue
      if (iprint() >= 120)
     .  print 333, x, y, opt, lintyp,linl1,linl2,linl3,linl4
  333 format(' Pltam: x,y=', 2f10.3, '  opt=',i1,
     .  '  lintyp=',i2,4f7.3)

      end
      subroutine mve(x,y)
C- Plot absolute position, user's units
C ----------------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   Nonorthogonal shear transformation.  General linear transformation:
Cr
Cr      (xs)   (a  b) (x)      (x)
Cr      (  ) = (    ) ( ) =  R ( )
Cr      (ys)   (c  d) (y)      (y)
Cr
Cr   We seek special shear with the following properties:
Cr
Cr         (1)   (1)         |   (1) |             |   (0) |
Cr       R ( ) = ( )   and   | R ( ) | = 1   and   | R ( ) | = 1
Cr         (0)   (0)         |   (0) |             |   (1) |
Cr
Cr   Also the the angle between xs,ys should be cos(theta). Then
Cr      a=1, c=0, b=cos(theta), d=sin(theta),   or
Cr
Cr      (xs)   (1  cos(theta)) (x)
Cr      (  ) = (             ) ( )
Cr      (ys)   (0  sin(theta)) (y)
Cr
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y,cs0(0:3)
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

C Local variables
      integer io,iprint,i1mach
      double precision cs(0:3),xs,ys

      cs(0) = -2d0
      io = 0
      goto 10

      entry drwe(x,y,cs0)
      call dcopy(4,cs0,1,cs,1)
      io = 2
      if (mod(cs(0),100d0) == 10) cs(0) = 5
      if (mod(cs(0),100d0) == 11) cs(0) = 7
      goto 10

      entry drw(x,y)
      io = 1
   10 continue
      xcu = x
      ycu = y
      goto 30

      entry mver(x,y)
      io = 0
      goto 20

      entry drwr(x,y)
      io = 1
      goto 20

      entry drwre(x,y,cs0)
      call dcopy(4,cs0,1,cs,1)
      io = 2
      if (mod(cs(0),100d0) == 10) cs(0) = 5
      if (mod(cs(0),100d0) == 11) cs(0) = 7

   20 continue
      if (logx) then
        xcu = xcu*dexp(x)
      else
        xcu = xcu+x
      endif
      if (logy) then
        ycu = ycu*dexp(y)
      else
        ycu = ycu+y
      endif

   30 continue
      xs = xcu
      if (logx) then
        if (xcu <= 0)
     .    call rx1('x=%g encountered, but log x scale',xcu)
        xs = dlog(xcu)
      endif
      ys = ycu
      if (logy) then
        if (ycu <= 0)
     .    call rx1('y=%g encountered, but log y scale',ycu)
        ys = dlog(ycu)
      endif
C  ... Nonorthogonal transformation
       xs = (offxu+xs)*scalxu + cosa*(offyu+ys)*scalyu
       ys = sina*(offyu+ys)*scalyu

      if (iprint() >= 110) then
        call awrit5(' mve (opt=%i): %1;4g %1;4g(UU)  %d %d(MU)',
     .    ' ',80,i1mach(2),io,xcu,ycu,xs,ys)
      endif

      call pltam(xs,ys,io,cs)
C      if (lrot) then
C        call pltam(ys,xs,io,cs)
C      else
C        call pltam(xs,ys,io,cs)
C      endif
      return

      entry pposu(x,y)
C- Return, update current pen pos'n by scaling from the medium's units
      xcu = xcm/scalxu - offxu
      ycu = ycm/scalyu - offyu
      if (logx) xcu = dexp(xcu)
      if (logy) ycu = dexp(ycu)
      x = xcu
      y = ycu
      end
      subroutine mveg(x,y)
C- Plot absolute position, graphics units
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y
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

C Local variables
      integer io
      double precision xcg,ycg,cs0(4)

      io = 0
      goto 10

      entry drwg(x,y)
      io = 1
   10 continue
      xcg = x
      ycg = y
      goto 30

      entry mvegr(x,y)
      io = 0
      goto 20

      entry drwgr(x,y)
      io = 1
   20 continue
      call pposg(xcg,ycg)
      xcg = x+xcg
      ycg = y+ycg

   30 continue
      cs0(1) = -2d0
      call pltam((offxg+xcg)*scalxg,(offyg+ycg)*scalyg,io,cs0)
C the following call updates xcu,ycu
      call pposu(xcg,ycg)
      end
      subroutine pposg(x,y)
C- Return absolute position, graphics units
C Passed parameters
      double precision x,y
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
C- Return current pen position by scaling from the medium's units
      x = xcm/scalxg - offxg
      y = ycm/scalyg - offyg
      end
      subroutine cdrw(x,y,xold,yold,lclip,box,cs)
C- Draw with explicit clips of boundaries
C ----------------------------------------------------------------------
Ci Inputs
Ci   lx     :T if x is on a log scale
Ci   ly     :T if y is on a log scale
Ci   x,y    :draw line to this point
Ci   xold,yold :draw line from this point
Ci   box    :bounding box:
Ci          :xmin=box(1),ymin=box(2),xmax=box(3),ymax=box(4)
Ci   lclip  :T clip line segments outside bounding box
Ci   cs     :if cs(0) = 0, has no effect
Ci          :otherwise, substitute drwe for drw, passing cs
Co Outputs
Cr Remarks
Cb Bugs
Cb   not implemented for log scale
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lclip
      double precision x,y,xold,yold,box(4),cs(0:3)
C ... Local parameters
      logical le1,le2
      integer n
      double precision xe1,ye1,xe2,ye2

      if (lclip) then
        call gtxeff(xold,yold,x,y,box,xe1,ye1,xe2,ye2,le1,le2,n)
C       Both points within clip region
        if (le1 .and. le2) then
          if (cs(0) == 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Only first point within clip region
        elseif (le1) then
          call drw(xe1,ye1)
          call mve(x,y)
C       First point outside, second point on the boundary
        elseif (le2 .and. n == 0) then
          if (cs(0) == 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Only second point within clip region
        elseif (le2) then
          call mve(xe1,ye1)
          if (cs(0) == 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Neither point within clip region
        else
          if (n == 0) goto 99
          call mve(xe1,ye1)
          if (cs(0) == 0) then
            call drw(xe2,ye2)
          else
            call drwe(xe2,ye2,cs)
          endif
          call mve(x,y)
        endif

C     Neither point within clip region
      else
        if (cs(0) == 0) then
          call drw(x,y)
        else
          call drwe(x,y,cs)
        endif
      endif

C     Exit, setting xold,yold to current x,y
   99 continue
      xold = x
      yold = y
      end
      subroutine gtxeff(x1,y1,x2,y2,box,xe1,ye1,xe2,ye2,le1,le2,n)
C- Given a line segment, return a segment lying within a clipped region
C ----------------------------------------------------------------------
Ci Inputs
Ci   lx     :T if x is on a log scale
Ci   ly     :T if y is on a log scale
Ci   x1,y1  :first x-y pair
Ci   x2,y2  :second x-y pair
Ci   box    :bounding box:
Ci          :xmin=box(1),ymin=box(2),xmax=box(3),ymax=box(4)
Ci   xe1,ye1:first point of subsegment
Ci   xe2,ye2:second point of subsegment
Co Outputs
Co   le1    :T if xe1,ye1 and x1,y1 are the same point
Co   le2    :T if xe2,ye2 and x2,y2 are the same point
Co    n     :number of boundary points
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical le1,le2
      double precision box(4),x1,y1,x2,y2,xe1,ye1,xe2,ye2
C ... Local parameters
      logical inside,betw
      integer n,nb
      double precision x,y,fuzz,alpha,beta,xx
      parameter (fuzz=1d-7)
      inside(x,y) = x+fuzz >= box(1) .and. x-fuzz <= box(3) .and.
     .              y+fuzz >= box(2) .and. y-fuzz <= box(4)
      betw(x,alpha,beta) =
     .  x >= min(alpha,beta) .and. x <= max(alpha,beta)

C     Sanity check
      if (box(3) <= box(1)) call rx('gtxeff: bad box')
      if (box(4) <= box(2)) call rx('gtxeff: bad box')

      le1 = inside(x1,y1)
      le2 = inside(x2,y2)
      xe1 = x1
      xe2 = x2
      ye1 = y1
      ye2 = y2
      if (le1 .and. le2) then
        return
      elseif (le1 .or. le2) then
        call gtxff2(1,x1,x2,y1,y2,box,xe1,ye1,n,nb)
        if (n /= 1 .and. nb == 0) then
          call rx('bug in gtxeff')
        endif
      else
        call gtxff2(1,x1,x2,y1,y2,box,xe1,ye1,n,nb)
        call gtxff2(2,x1,x2,y1,y2,box,xe2,ye2,n,nb)
        if (n == 0) return
        if (n /= 2) call rx('bug in gtxeff')
        if (betw(xe2,x1,xe1)) then
          xx = xe2
          xe2 = xe1
          xe1 = xx
          xx = ye2
          ye2 = ye1
          ye1 = xx
        endif
      endif

      end
      subroutine gtxff2(n,x1,x2,y1,y2,box,xe,ye,k,kb)
C- Return nth occurrence of intersection of box w/ line seg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n      :nth occurence to seek.  Should be 1 or 2.
Ci   x1,y1  :first x-y pair of line segment
Ci   x2,y2  :second x-y pair of line  segment
Ci   box    :box
Co Outputs
Co   xe     :x coordinate of last occurence
Co   ye     :y coordinate of last occurence
Co   k      :last occurence found.  If successful, k=n. See Remarks
Co   kb     :number of points on the boundary. See Remarks
Cr Remarks
Cr   A point on the box boundary will be returned in (xe,ye),
Cr   but k will not reflect this.  Caller can check by evaluating
Cr   distance between (xe,ye) and box.
Cu Updates
Cu   Added argument kb
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,k,kb
      double precision x1,x2,y1,y2,xe,ye,box(4)
C ... Local parameters
      logical betw
      integer m
      double precision isec,xa,dy,sloi,fuzz,alpha,beta,x,xmin,xmax
      parameter (fuzz=1d-4)
      isec(sloi,xa,dy) = xa + dy*sloi
C     betw => x between alpha and beta and between xmin and xmax
      betw(x,alpha,beta,xmin,xmax) =
     .  x >= min(alpha,beta) .and. x <= max(alpha,beta) .and.
     .  x >= xmin .and. x <= xmax

      k = 0
      kb = 0

C     Intersection of two horizontals with line segment
      if (y2 /= y1) then
        do  m = 2, 4, 2
          ye = box(m)
          xe = isec((x2-x1)/(y2-y1),x1,ye-y1)
          if (betw(xe,x1,x2,box(1),box(3))) then
C         Exclude points coincident with (x1,y1) or (x2,y2)
          if (min(dabs(xe-x1)+dabs(ye-y1),dabs(xe-x2)+dabs(ye-y2))
     . > fuzz) then
            k = k+1
            if (k == n) return
          else
            kb = kb+1
          endif
          endif
        enddo
      endif

C     Intersection of two verticals with line segment
      if (x2 /= x1) then
        do  m = 1, 3, 2
          xe = box(m)
          ye = isec((y2-y1)/(x2-x1),y1,xe-x1)
          if (betw(ye,y1,y2,box(2),box(4))) then
C         Exclude points coincident with (x1,y1) or (x2,y2)
          if (min(dabs(xe-x1)+dabs(ye-y1),dabs(xe-x2)+dabs(ye-y2))
     . > fuzz) then
            k = k+1
            if (k == n) return
          else
            kb = kb+1
          endif
          endif
        enddo
      endif

C     call rx('bug in gtxeff')

      end
      subroutine plwstr(string,mode)
C- Write string directly to output stream
Ci mode=1: string is name of a file. Insert contents of the file
      implicit none
      character*(*) string
      integer mode
C Local
      integer ifi,fopng,i,iprint,i1mach
      character*1024 s
      logical rdstrn
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

      if (mode == 1) then
        ifi = fopng(string,-1,1)
        s(1:60) = string
        if (iprint() >= 40) call awrit0(' plwstr:  '//
     .    'inserting file "'//s(1:60)//'%a"',' ',80,i1mach(2))
        call awrit0('%N%% --- plwstr:  inserting file "'//
     .    s(1:60)//'%a"',' ',80,stream)
   10   if (rdstrn(ifi,s,len(s),.false.)) then
          call skpblb(s,len(s),i)
          write (stream,'(a)') s(1:i+1)
          goto 10
        endif
        call fclose(ifi)
      else
        s(1:60) = string
        if (iprint() >= 40) call awrit0(' plwstr:  '//
     .    'inserting string "'//s(1:60)//'%a"',' ',80,i1mach(2))
        call awrit0('%N%% --- plwstr:  inserting string "'//
     .    s(1:60)//'%a"',' ',80,stream)
        write (stream,'(a)') string
      endif

      end
      subroutine pltdmp(io)
C- Dumps plot file
C ----------------------------------------------------------------
Ci Inputs
Ci   io:  0, dump plot immediately; 1, wait for carriage return
Co Outputs
Cr Remarks
Cr  This is only meaningful if plot is dumped to device at one time.
Cr  Hercules card: restore screen to ASCII mode when plot is complete.
Cr  When applicable and io=1, prompt for a carriage return first.
Cr  PostScript file: write out showpage
C ----------------------------------------------------------------
      implicit none
      integer io
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
      if (medium == 1) then
      elseif (medium == 2) then
      elseif (medium == 3) then
        if (lrot) write(stream,031)
  031   format(' 0 612 translate -90 rotate')
        call awrit4('%%%%BoundingBox: %;1d %;1d %;1d %;1d',' ',
     .    80,stream,bbox(1),bbox(3),bbox(2),bbox(4))
        write (stream,10)
   10   format('showpage')
      endif
C     if (isopen(stream,.false.)) call fclose(stream)

      call fclose(stream)

      end
      subroutine pltcol(color)
C- Set pen color of current medium
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer color
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

      pcol = color
      end
      subroutine plndpi(n)
C- Set number of dots-per-inch
      implicit none
      integer n
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


      ndpi = n
      ndwrt = int(dlog10(ndpi*2.54d0*scalxc)) + 1

      end

      subroutine plntyp(type,bold,l1,l2,l3,l4)
C- Set line type
C ----------------------------------------------------------------
Ci Inputs
Ci   type: 0, move, pen up
Ci         1, draws simple lines
Ci         2, draws broken lines, pen down for length (l1-l2) followed
Ci            by pen up for length l2.
Cr   l1:   repeat length of the line type.  Scale is
Cr         hardware dependent but defined so that 1 corresponds
Cr         to a visually pleasing period.  Scale can be changed
Cr         in pltini.
Cr   l2:   length for which pen is up during a draw.
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer type,bold
      double precision l1,l2,l3,l4

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
      integer is
      character s*100
      integer lintys,ltblds
      double precision linl1s,linl2s,linl3s,linl4s
      save lintys,ltblds,linl1s,linl2s,linl3s,linl4s

C --- Scale linetype lengths from "sensible" scales (numbers of order 1) to
C     medium units ---
      lintyp = type
      ltbld  = bold
      linl1 = l1*lnprd
      linl2 = l2*lnprd
      linl3 = l3*lnprd
      linl4 = l4*lnprd
      lnused = 0

   10 if (medium /= 3) return
      s = ' '
      is = 0
      if (lintyp > 0) then
        call bin2a('d2',1,0,dble(ltbld)/2,4,0,100,s,is)
        call bin2a('setlinewidth 1 setlinejoin [',1,0,0,1,0,80,s,is)
        if (lintyp == 1) then
          is = is+1
        elseif (lintyp == 2) then
          call bin2a('d2',0,0,linl1,4,0,100,s,is)
          call bin2a('d2',1,0,linl2,4,0,100,s,is)
          if (linl3 /= 0) then
            call bin2a('d2',1,0,linl3,4,0,100,s,is)
            call bin2a('d2',1,0,linl4,4,0,100,s,is)
          endif
        elseif (lintyp == 3) then
          is = is-1
          call bin2a('1 setlinecap [0',0,0,0,1,0,80,s,is)
          call bin2a(' ',1,0,ltbld,2,0,100,s,is)
        endif
        call bin2a('] 0 setdash',0,0,0,1,0,80,s,is)
      endif
      if (is > 0) write (stream,*) s(1:is)
      return

      entry plntyr(type,bold,l1,l2,l3,l4,err)
C- Retrieves current line type
      err = lnused /= 0
      if (err) return
      type = lintyp
      bold = ltbld
      l1 = linl1/lnprd
      l2 = linl2/lnprd
      l3 = linl3/lnprd
      l4 = linl4/lnprd
      return

      entry plntys
C- Saves current line type
      lintys = lintyp
      ltblds  = ltbld
      linl1s = linl1
      linl2s = linl2
      linl3s = linl3
      linl4s = linl4
      return

      entry plntyg
C- Restores saved line type
      lintyp = lintys
      ltbld  = ltblds
      linl1 = linl1s
      linl2 = linl2s
      linl3 = linl3s
      linl4 = linl4s
      goto 10

      end
      subroutine ticset(llog,mode,mt0,bot,top,ts0,tx0,ch,lmode,
     .  tic,mt,mt1,itic,s,rmt,rmnt)
C- Set a tic pattern
C ----------------------------------------------------------------------
Ci Inputs
Ci   llog : T if log scale, F if linear scale
Ci   mode : default tic pattern.  Actual tic pattern is returned in
Ci        : lmode, which see for description.
Ci        : mode>=0 => lmode is assigned to mode
Ci        : mode<0 => ticset chooses a default lmode
Ci   mt0  : default number of tics per major tic.  Actual number is
Ci        : returned in mt, which see for description.
Ci        : mt0<0 => ticset chooses a value for mt.
Ci bot,top: axis bounds
Ci   ts0  : (mode=0 only) default tic spacing.  Uses only ts0(1).
Ci        : ts0(1)<0 => ticset chooses ts
Ci        : (mode>10) vector of positions at which to set tic marks
Ci   tx0  :force tic through tx0
Ci        :tx0 < bot => ticset chooses a value
Ci   ch   :a character of length 1, used for printout
Ci  rmt,rmnt :major and minor tic sizes (for printout only)
Co Outputs
Co  lmode :specifies tic pattern; see input 'mode' for description
Co        : mode=0: uniformly spaced tics; this is the
Co        :         usual mode for linear scale, but may be
Co        :         appropriate for log scale also.
Co        : ... modes 1..3 are appropriate for log scale:
Co        : mode=1: tics at 1,2,3,4,5,6,7,8,9 * 10^integer
Co        : mode=2: tics at 1,2,5 * 10^integer
Co        : mode=3: tics 1 * 10^integer
Co        :
Co        : mode>10: use supplied vector of tic marks, passed
Co        :        : in array ts0 (mode-10 points)
Co   tic  : position of 1st tic
Co   mt   : number of of tics between major tics, i.e. for tic mark i,
Co        : major tic at any (i+mt1 modulus mt) = 0
Co        : minor tic at any (i+mt1 modulus mt) > 0
Co   mt1  : offset for modulus of mt at 1st tic; see description for mt
Co   ts   : lmode=0, lmode>10; spacing between first and second tics
Co        : other values of lmode: not set (ts has no meaning)
Cr  Remarks
Cr   *tic marks are of two types, major and minor.
Cr    mt and mt1 specify of which type each tic mark belongs.
Cr   *The spacing between tic marks depends on lmode.
Cr    tic is the position of the first tic mark.
Cr    ts specifies spacing between 1st and 2nd mark (lmode=0,10)
Cr    For the log scale modes (lmode=1..3), see description of lmode
Cr    ts0 is a list of all marks (lmode=10)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical llog
      character ch*1, s*(*)
      double precision bot,top,ts0(*),tx0,ts,tic,rmt,rmnt
      integer mode,mt0,mt,mt1,lmode,itic
C ... Local parameters
      double precision botl,topl,l10,width,ticspc,xx,fuzfac
      parameter (fuzfac=1d-10)
      integer j,k

C --- Make lmode ---
      lmode = mode
      if (lmode < 0 .and. .not. llog) lmode = 0
      if (lmode < 0 .and. llog) then
        botl = dlog(bot)
        topl = dlog(top)
        l10  = dlog(10d0)
        width = dabs(topl-botl)/l10
        lmode = 3
        if (width < 3d0)  lmode = 2
        if (width < 1d0)  lmode = 1
        if (width < .5d0) lmode = 0
      endif

C --- Make mt ---
      mt = mt0
      if (mt == -1) then
        mt = 1
        if (lmode == 0) mt = 2
        if (lmode == 2) mt = 3
      endif

C     Shift bot,top by fuzz factor
      botl = min(top,bot)
      topl = max(top,bot)
      if (llog) then
C       fuzz = (topl/botl)*fuzfac
        botl = botl/(1+fuzfac)
        topl = topl*(1+fuzfac)
      else
C       fuzz = (topl - botl)*fuzfac
        botl = botl - fuzfac
        topl = topl + fuzfac
      endif

C --- ts, tic and mt1 for mode=0 ---
      if (lmode == 0) then
        ts = ts0(1)
        if (ts0(1) <= 0) ts = ticspc(dabs(top-bot)/10,1)
        if (top < bot) ts = -ts
        tic = tx0
        mt1 = 0
        if (tic < botl) then
          j = int(botl/ts)
          if (j*ts < botl) j = j+1
          tic = j*ts
          mt1 = mod(abs(j),mt)
        else
   11     continue
          if (tic-ts >= botl) then
            tic = tic-ts
            mt1 = mt1+1
            goto 11
          endif
          mt1 = mod(mt1,mt)
        endif
        itic = 1
C ...   For printout %10z => No leading 0, 1 absolute precision
        call awrit5('%10z%a'//ch//'t1=%1;4g ts'//ch//'=%1;4g mt'//ch//
     .    '=%i h=(%g,%g)',s,len(s),0,tic,ts,mt,rmt,rmnt*rmt)
        ts0(1) = ts
        return
      endif

C --- ts, tic and mt1 for mode>10 ---
      if (lmode > 10) then
        tic = ts0(1)
        if (mode == 11) then
          ts = 0
        else
          ts = ts0(2) - ts0(1)
        endif
        itic = 1
        k = 1
        if (tx0 < botl .or. tx0 > topl) goto 22
C ...   Determine mt1 by finding point closest to tx0
        xx = dabs(tx0-ts0(1))
        do  20  j = 1, mode-10
          if (dabs(tx0-ts0(j)) < xx) then
            k = j
            xx = dabs(tx0-ts0(j))
          endif
   20   continue
   22   continue
        mt1 = mod(k-1,mt)
        call awrit5('%10z%a  '//ch//'t1=%1;4g mode=%i(%i tics) mt'//
     .    ch//'=%i:%i',s,len(s),0,tic,5,lmode-10,mt,mt1)
        return
      endif

C --- ts, tic, mt1 and itic for mode=1,2,3 ---
      if (lmode >= 1 .and. lmode <= 3) then
        ts = 0
        tic = dmax1(tx0,botl)
        itic = 0
        call ticnxt(lmode,tic,topl,[0d0],itic)
        if (tic < botl) call ticnxt(lmode,tic,topl,[0d0],itic)
        mt1 = mod(mt-itic+1,mt)
        call awrit3('%10z '//ch//'t1=%1;4g mode=%i mt'//
     .    ch//'=%i',s,len(s),0,tic,lmode,mt)
        return
      endif

C      print *, 'ticset:', j,mt1
C      pause
      call rx('ticset: mode not defined')

      end
      subroutine ticnxt(mode,tic,top,ts0,inew)
C- Finds next tic, depending on mode (see ticset)
Ci mode; tic: last tic; ts: tic spacing (mode=0)
Co inew,tic: next tic, and index to next tic.
      implicit none
      integer mode
      double precision tic,top,ts0(1),d1mach
      integer is,i,inew,nx(3),iprint
      double precision tx(11,3),scale,ss,err
      save scale
      data tx /1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,
     .         1d0,2d0,5d0,10d0,20d0,0d0,0d0,0d0,0d0,0d0,0d0,
     .         1d0,10d0,100d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/,
     .  nx /10,4,2/

      if (mode == 0) then
        tic = tic + ts0(1)
        inew = inew+1
        if (tic > top) inew = 0
        goto 100
      endif

      if (mode > 10) then
        inew = inew+1
        if (inew > mode-10) then
          inew = 0
        else
          tic = ts0(inew)
        endif
        goto 100
      endif

C ... for the logarithmic modes, reset tic to closest table value
      if (inew == 0) then
        scale = nint(dlog10(dabs(tic)))
        if (10**scale > dabs(tic)) scale = scale-1
        ss = dabs(tic)/10**scale

        err = 9d9
        do  10  i = 1, nx(mode)
          if (dabs(ss/tx(i,mode)-1) > err) goto 10
          is = i
          err = dabs(ss/tx(i,mode)-1)
   10   continue
        inew = mod(is-1,nx(mode)+1)
      endif
      inew = inew+1
      if (inew > nx(mode)) then
        inew = 2
        scale = scale+1
      endif
      tic = tx(inew,mode)*10**scale
      if (tic > top*(1+4*d1mach(3))) inew = 0

  100 continue
      if (iprint() > 90) print 333, tic, mode, inew
  333 format(' ticnxt:  tic=',f9.6,' mode=',i2,' inew=',i2)

      end
      function ticspc(width,n)
C- Procedure determines tic spacing given a width, linear axes
C ----------------------------------------------------------------
Ci Inputs
Ci   width
Ci   n:  not used for now
Co Outputs
Co   ticspc
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision width,ticspc
      integer n
C Local parameters
      double precision scale,tmp
C This construction avoids taking integer part of negative number
      scale = dlog10(width)
      if (scale > 0) then
        scale = int(scale)
      else
        scale = -int(-scale) - 1
      endif
      scale = 10d0**(-scale)
C tmp is the width scaled to a number between 1 and 10
      tmp = width*scale
      if (tmp > 7.5d0) then
        tmp = 10
      elseif (tmp > 3.5d0) then
        tmp = 5
      elseif (tmp > 1.5d0) then
        tmp = 2.5d0
      else
        tmp = 1
      endif
      ticspc = int(tmp + .00001d0)/scale
      end
      subroutine frme(tx0,tsx,mtx,rmtx,rmntx,modx,xnum,xnfmt,xlabel,yab,
     .                ty0,tsy,mty,rmty,rmnty,mody,ynum,ynfmt,ylabel,xor,
     .                sfxdef,np,xp,fxdef,title,lt,frmefont,frmcs,fillcs)
C- Draw a frame around the plot and either tic marks or grid lines.
C ----------------------------------------------------------------
Ci Inputs
Ci   ... the following apply to the x-axis
Ci   tx0:  Force a major tic through this point
Ci         tx0 < plul makes default
Ci   tsx:  spacing in user's units between tic marks along the x axis
Ci         tsx <= 0 makes default of approximately 1/10 grid width
Ci   mtx:  number of major tics per tic, x axis
Ci   rmtx: size of major tic, in proportion to plut-plub (.03 is good)
Ci   rmntx:size of minor tic, in proportion to major (.6 is good).
Ci   modx :tic mode (see ticset)
Ci   xnum :1  to number axis on TOP
Ci        :0  to number axis on
Ci        :-1 to number axis on BOTTOM
Ci   xnfmt:awrite format string for x-axis numbering
Ci   xlabel:label for abscissa
Ci   ... similar variables for y-axis
Ci   --- these variables are a work in progress ---
Ci   sfxdef:a string for mapping absicssa tic marks to a function
Ci   np    :number of points describing function
Ci   xp    :x-points for function
Ci   fxdef :y-points for function
Ci   ---
Ci   title :
Ci   lt    :line thickness:
Ci          lt(1) is line thickness
Ci          lt(2) = 0 draws both top and bottom axis
Ci                  1 draws bottom only
Ci                  2 draws top only
Ci                  3 draws neither
Ci          lt(3) = 0 draws both left and right axes
Ci                  1 draws left only
Ci                  2 draws right only
Ci                  3 draws neither
Ci          NB: the frame is closed and filled only if lt(2)+lt(3)=0
Ci   frmefont : if not blank, specifies font for frame labelling (see routine parsefont)
Ci   frmcs :frame color specification.  frmcs(0) is not used.
Ci   fillcs:frame fill specification.
Ci         :fillcs(4)=1 => fill frame before drawing
Co Outputs
Cr Remarks
Cr   Axes are drawn.
Cu Updates
Cu   19 Feb 03 Added numbering to rhs to top
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mtx,mty,lt(3),modx,mody,np
      double precision tx0,tsx(1),rmtx,rmntx,ty0,tsy(1),rmty,rmnty,
     .  frmcs(0:3),fillcs(0:4),xp(np),fxdef(np),yab,xor
      character*(*) xlabel,ylabel,title,sfxdef,frmefont
      integer xnum,ynum

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
      integer NULLI
      parameter (NULLI=-99999)
      integer i,mt,it,ipr,lmode,mt1,side,jx,stdo,fontsz
      double precision xtic,ytic,ticlen,xl,yl,xxtic,dxtic,csloc(0:3)
      character*20 xnfmt,ynfmt,font*19
      character*120 s,s2,fmt
      procedure(integer) :: i1mach,lgunit

      double precision l,r,b,t,ll,rr,bb,tt,tolx,toly
      logical maktic,twolin

C --- Setup ---
      twolin = xor /= NULLI .or. yab /= NULLI
      call getpr(ipr)
      stdo = lgunit(1)
      l = plul
      r = plur
      if (scalxu < 0) then
        r = plul
        l = plur
      endif
      b = plub
      t = plut
      if (logx) then
        l = dexp(plul)
        r = dexp(plur)
        if (scalxu < 0) then
          r = dexp(plul)
          l = dexp(plur)
        endif
      endif
      ll = l
      rr = r
      if (xor /= NULLI .and. .not. logx) ll = xor
      if (xor /= NULLI .and. logx) ll = dexp(xor)
      if (logy) then
        b = dexp(plub)
        t = dexp(plut)
      endif
      bb = b
      tt = t
      if (yab /= NULLI .and. .not. logy) bb = yab
      if (yab /= NULLI .and. logy) bb = dexp(yab)
      tolx = (r-l)*1d-6
      toly = (t-b)*1d-6

      if (ipr >= 10) then
        fmt = ' FRME : axes at x=%1;4g %1;4g'
        if (abs(l) > 1e10) fmt(22:23) = '3e'
        if (abs(r) > 1e10) fmt(28:29) = '3e'
        call awrit2(trim(fmt),s,len(s),0,l,r)
        if (logx) call awrit1('%a (log)',s,len(s),0,l)
        call awrit2('%a  y=%1;4g %1;4g',s,len(s),0,b,t)
        if (logy) call awrit1('%a (log)',s,len(s),0,l)
        call awrit1('%a  bold=%i',s,len(s),0,lt)
        if (lt(2)+lt(3) > 0)
     .    call awrit4('%a: no axes for: '//
     .    '%?#n>0&n<>2# top,##%?#n>1# bottom,##'//
     .    '%?#n>0&n<>2# right,##%?#n>1# left,##%b',
     .    s,len(s),0,lt(2),lt(2),lt(3),lt(3))
      endif

C ... Set font size, informational printout
      if (medium == 3) then

C     Information about fill data, font
      s2 = ' '
      if (fillcs(4) /= 0)
     .  call awrit2('%a  fill cs=%d %s,col=%3;4d ',s2,len(s2),0,fillcs,fillcs(1))
      if (frmefont /= ' ')
     .  call awrit0('%a  font='//frmefont,s2,len(s2),0)

      call awrit0('%% ---'//s//'%a'//trim(s2),' ',len(s),stream)

C ... Set frame font
      if (frmefont /= ' ') then
        call parsefont(frmefont,font,fontsz,.true.) ! return font, fontsz
        call setfnt(font,fontsz)
      endif

      endif

      if (ipr >= 10) then
        call awrit0('%a'//s2,s,len(s),0)
        write(stdo,"(a)") trim(s)
      endif

      write (stream, 333)
  333 format('/max 0 def ')
      call plntys
      call plntyp(1,lt(1),0d0,0d0,0d0,0d0)

C ... Fill frame with background color
      csloc = frmcs
      if (fillcs(4) /= 0) then
        call mve (l, b)
        call drw (l, t)
        call drw (r, t)
        call drwe(r,b,fillcs)   ! fillcs(0) should be 103 for fill without lines
        csloc(0) = 0            ! No fill when frame is drawn
      endif

C ... Set frame color
      call awrit1('gsave %3:1;4d setrgbcolor',' ',80,stream,frmcs(1))

C --- Draw frame ---
      call mve (ll, b)
      if (twolin) then
        call drw (ll, t)
        call mve (l, bb)
        call drwe(r, bb, csloc)
      else
      if (lt(2) == 0 .or. lt(2) == 1) then
        call drw (r, b)
      else
        call mve (r, b)
      endif
      if (lt(3) == 0 .or. lt(3) == 2) then
        call drw (r, t)
      else
        call mve (r, t)
      endif
      if (lt(2) == 0 .or. lt(2) == 2) then
        call drw (ll, t)
      else
        call mve (ll, t)
      endif
      if (lt(2)+lt(3) == 0) then
        call drwe(ll, b, csloc)
      elseif (lt(3) == 0 .or. lt(3) == 1) then
        call drw (ll, b)
C       if (lt(2) == 0 .or. lt(2) == 1) call drw (r, b)
C      else
C        call mve(ll, b)
C        call drwre(0d0, 0d0, csloc)
      endif
      endif
*     write (stream,'(''% we are done'')')

c --- Draw and label tic marks, x axis ---
C ... tic marks for each side
      ytic = bb
      ticlen = rmtx*(plut-plub)
      if (logy) ticlen = rmtx*dlog(t/b)
      s = ' '
      do  side = 1, 2
      if (twolin .and. side == 2) cycle
C ... Tic mark setup
      s2 = ' '
      call ticset(logx,modx,mtx,l,r,tsx,tx0,'x',lmode,xtic,mt,mt1,it,s2,rmtx,rmntx)
      if (side == 1) s(2:) = s2
      maktic = (lt(2) == 0 .or. lt(2) == side) .and. lmode /= 10
      if (.not. maktic) goto 16
C ... for each tic do (mt1 is zero at major tics)
      do  i = 1, 999
C   ... for now
        xxtic = xtic
        if (sfxdef(1:1) /= ' ') then
          xxtic = xtic
          call polint(fxdef,xp,np,4,xxtic,0d0,0,jx,xtic,dxtic)
          if (ipr >= 50) then
            call awrit2(' plsym: abscissa %1:-1,4;4g mapped'//
     .        ' to %1:-1,4;4g',' ',80,stdo,xxtic,xtic)
          endif
        endif
        call mve(xtic, ytic)
        if (mt1 == 0) then
          if (twolin .and. ytic /= b) call mver(0d0, -ticlen/2)
          call drwre(0d0, ticlen, [-2d0,0d0,0d0,0d0])
          if (ytic == bb .and. xnum > 0) then
            if (medium == 3) then
              call mve(xtic,ytic)
              if (twolin .and. ytic /= b) call mver(0d0, -ticlen/2)
C             Omit label if the y axis passes through it
              if (twolin .and. abs(xtic-ll) < tolx) then
              else
              call pslabl('a',xnfmt,xxtic,' ','l',' ',1,0d0)
              endif
            endif
          endif
          if (ytic == t .and. xnum < 0) then
            if (medium == 3) then
              call mve(xtic,ytic)
              call pslabl('A',xnfmt,xxtic,' ','l',' ',1,0d0)
            endif
          endif
          mt1 = mt
        else
          if (twolin .and. ytic /= b) call mver(0d0, -rmntx*ticlen/2)
          call drwre(0d0, rmntx*ticlen, [-2d0,0d0,0d0,0d0])
        endif
   14   continue
        mt1 = mt1-1
        call ticnxt(lmode,xtic,r,tsx,it)
        if (it == 0) goto 16
      enddo
      call rx('FRME:  too many tics in x')
   16 continue
      ytic = t
      ticlen = -ticlen
      enddo

c --- Draw and label tic marks, y axis ---
C ... tic marks for each side
      xtic = ll
      ticlen = rmty*(plur-plul)
      if (logx) ticlen = rmty*dlog(r/l)
      do  side = 1, 2
      if (twolin .and. side == 2) cycle
C ... Tic mark setup
      s2 = ' '
      if (side == 2) s2 = s
      call ticset(logy,mody,mty,b,t,tsy,ty0,'y',lmode,ytic,mt,mt1,it,s2(3+len_trim(s2):),rmty,rmnty)
      maktic = (lt(3) == 0 .or. lt(3) == side) .and. lmode /= 10
      if (.not. maktic) goto 26
C ... for each tic do (mt1 is zero at major tics)
      do  i = 1, 999
        call mve(xtic, ytic)
        if (mt1 == 0) then
          if (twolin .and. xtic /= l) call mver(-ticlen/2, 0d0)
          call drwre(ticlen, 0d0, [-2d0,0d0,0d0,0d0])
          if (xtic == ll .and. ynum > 0) then
            if (medium == 3) then
              call mve(xtic,ytic)
              if (twolin .and. xtic /= l) call mver(-ticlen/2, 0d0)
              if (twolin .and. abs(ytic-bb) < toly) then
              else
              call pslabl('o',ynfmt,ytic,' ','l',' ',1,0d0)
              endif
            endif
          endif
          if (xtic == r .and. ynum < 0) then
            if (medium == 3) then
              call mve(xtic,ytic)
              call pslabl('O',ynfmt,ytic,' ','l',' ',1,0d0)
            endif
          endif
          mt1 = mt
        else
          if (twolin .and. xtic /= l) call mver(-rmnty*ticlen/2, 0d0)
          call drwre(rmnty*ticlen, 0d0, [-2d0,0d0,0d0,0d0])
        endif
   24   continue
        mt1 = mt1-1
        call ticnxt(lmode,ytic,t,tsy,it)
        if (it == 0) goto 26
      enddo
      call rx('FRME:  too many tics in y')
   26 continue
      xtic = r
      ticlen = -ticlen
      enddo
      if (ipr >= 30) call awrit0('%a',s2,-len(s),-i1mach(2))
      call plntyg
      if (medium /= 3) then
        call plntyg
        return
      endif
      s = s2
      s = '% --- '//s2(2:)
      call awrit0('%a',s,len(s),-stream)

C --- Label axes ---
      write (stream,*) '% FRME: Label axes ...'
C ... abcissa
      xl = 0.5d0*(plul + plur)
      yl = plub
      write (stream,130) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('0',' ',0d0,xlabel,'h','f',1,0d0)
      if (xlabel /= ' ') call pbbox((offxu+xl)*scalxu,(offyu+yl)*scalyu-3*fsiz)
C ... ordinate
      xl = plul
      yl = 0.5d0*(plut + plub)
      write (stream,230) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('0',' ',0d0,ylabel,'v','r',1,0d0)
      if (ylabel /= ' ') call pbbox((offxu+xl)*scalxu-3*fsiz,(offyu+yl)*scalyu)

C --- Add title ---
      write (stream,*) '% FRME: adding title ...'
      xl = plur
      yl = plut
      write (stream,330) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('4',' ',0d0,title,'h','r',1,0d0)
      if (title /= ' ')
     .  call pbbox((offxu+xl)*scalxu,(offyu+yl)*scalyu+fsiz)

      call awrit0('%% pop gsave for frame color%Ngrestore',' ',80,stream)
      call plntyg

  130 format(2f8.1,' moveto 0 h 3 mul neg rmoveto')
  230 format(2f8.1,' moveto h 2 mul max add neg 0 rmoveto')
  330 format(2f8.1,' moveto 0 h 0.5 mul rmoveto')
      end
      subroutine plcrv(xp,yp,wp,wp2,wp3,wp4,np,xmin,xmax,ymin,ymax,opts,
     .  cs,csw,iclos)
C- Draws or continues a smooth curve
C ----------------------------------------------------------------------
Ci Inputs:
Ci   xp,yp,np: set of x-y pairs, and number
Ci   wp,wp2 : color weights, used if csw is set
Ci   xmin..ymax: bounding box
Ci   opts :1s digit
Ci        :if nonzero, start a new curve whenever abscissa decreases
Ci        :10s digit
Ci        :if nonzero, plcrv clips curves to when they go outside
Ci         bounding box
Ci   iclos: 00: continues on existing curve, if one already open;
Ci              otherwise opens a new curve.
Ci          01: closes curve after drawing last point
Ci          10: closes any existing curve before starting new curve
Ci          11: both 10 and 01 above
Co Outputs:
Co   points smoothly connected on the plot stream
Cr Remarks
Cr   When a single color weight wp is present, color determined by
Cr       (1-wp)*cs + wp*csw
Cr   When two color weights are present, color determined by
Cr       (1-wp-wp2)*cs + wp*csw(1) + wp2*csw(2)
Cr   When three color weights are present, color determined by
Cr       (1-wp-wp2-wp3)*cs + wp*csw(1) + wp2*csw(2) + wp3*csw(3)
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer np,iclos,opts
      double precision xp(np),yp(np),wp(np),wp2(np),wp3(np),wp4(np),
     .  xmin,xmax,ymin,ymax,xi,yi,cs(0:3),csw(0:3,4)
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
C Local Variables
      logical lclip,lxbrk
      integer ip,ipr,getdig,i1,is,i1mach
      character s*120, s2*120, fmt*20
      double precision xsmin,ysmin,ysmax,xsmax,csl(0:3),xold,yold
      double precision box(4),cspt(0:3),w1,w2,w3,w4
      save xold,yold

      lxbrk = mod(opts,10) /= 0
      lclip = mod(opts/10,10) /= 0

      if (lintyp == 0) return
      call awrit2('%x(2f%i.%i,'' lto'')',fmt,20,0,ndwrt+8,ndwrt)
      call getpr(ipr)
C ... Local copy of cs which 10s digit of cs(0) is stripped
      call dcopy(4,cs,1,csl,1)
      if (cs(0) >= 10) csl(0) = mod(cs(0),10d0)
C ... Close curve if requested
      is = 0
      s = ' '
      if (getdig(iclos,1,10) == 1)
     & call plcrve(medium,pstat,[-2d0,0d0,0d0,0d0],s,is)
      if (is > 0) write(stream,*) s(1:is)
      is = 0
      s = ' '
      if (pstat == 0) then
        call bin2a('plcrv: new',1,0,0,1,0,120,s,is)
        i1 = 2
        xold = xp(1)
        yold = yp(1)
      else
        call bin2a('plcrv: continue',1,0,0,1,0,120,s,is)
        i1 = 1
        call pposu(xold,yold)
      endif
      call awrit3('%a curve, %i pts lt %i bold %i',s,120,0,np,lintyp,
     .  ltbld)
      if (cs(0) > -2d0) call awrit1('%a col=%3:1d',s,120,0,cs(1))
      if (csw(0,1) > 0d0)
     .  call awrit1('%a colw=%3:1d',s,120,0,csw(1,1))
      if (csw(0,2) > 0d0)
     .  call awrit1('%a colw2=%3:1d',s,120,0,csw(1,2))
      if (csw(0,3) > 0d0)
     .  call awrit1('%a colw3=%3:1d',s,120,0,csw(1,3))
      if (csw(0,4) > 0d0)
     .  call awrit1('%a colw4=%3:1d',s,120,0,csw(1,4))
      if (lintyp == 2 .and. linl3 == 0)
     .  call awrit2('%a len=(%d %d)',s,120,0,linl1,linl2)
      if (lintyp == 2 .and. linl3 /= 0) call
     .  awrit4('%a len=(%d %d %d %d)',s,120,0,linl1,linl2,linl3,linl4)
      if (ipr >= 40) call awrit1('%a',s,120,-i1mach(2),1)

C --------------- Postscript -----------------
      if (medium /= 3) goto 100
      s2 = '% ---' // s
      call awrit0('%a ---',s2,120,-stream)
C...  Clip beyond plot limits
      if (pstat == 0) then
        xsmin = (offxu+xmin)*scalxu
        ysmin = (offyu+ymin)*scalyu
        ysmax = (offyu+ymax)*scalyu
        xsmax = (offxu+xmax)*scalxu
        if (logx) then
          xsmin = (offxu+dlog(xmin))*scalxu
          xsmax = (offxu+dlog(xmax))*scalxu
        endif
        if (logy) then
          ysmin = (offyu+dlog(ymin))*scalyu
          ysmax = (offyu+dlog(ymax))*scalyu
        endif
        call pbbox(xsmax,ysmax)
        call pbbox(xsmin,ysmin)
        call awrit8('gsave newpath %1;nd %1;nd moveto %1;nd %1;nd '//
     .    'lineto',' ',120,stream,ndwrt,xsmin,ndwrt,ysmin,ndwrt,xsmin,
     .    ndwrt,ysmax)
        call awrit8('              %1;nd %1;nd lineto %1;nd %1;nd '//
     .    'lineto closepath clip',' ',120,stream,ndwrt,xsmax,
     .    ndwrt,ysmax,ndwrt,xsmax,ndwrt,ysmin)
        if (cs(0) >= 10) call awrit1('%3:1;4d setrgbcolor',' ',80,
     .    stream,cs(1))
        call mve(xp(1),yp(1))
      endif

      if (lclip) then
        box(1) = xmin
        box(2) = ymin
        box(3) = xmax
        box(4) = ymax
      endif

C --- Draw the curve ---
      if (lxbrk) then
        if (csw(1,1) /= 0) call rx('not ready for lxbrk and colw')
        do  1  ip = i1, np
        if (ip > i1) then
          if (xp(ip) > xp(ip-1)) then
            call drw(xp(ip),yp(ip))
          else
            if (ipr >= 4) call awrit1(
     .        '    ... break at point %i',' ',80,i1mach(2),ip)
            call mve(xp(ip),yp(ip))
          endif
        else
          call drw(xp(ip),yp(ip))
c         write(stream,fmt)  (offxu+xp(ip))*scalxu,(offyu+yp(ip))*scalyu
        endif
    1   continue
      else
        cspt(0) = 0
        if (csw(0,1) /= 0) then
          w1 = (wp(i1)+wp(max(i1-1,1)))/2
          w2 = 0
          w3 = 0
          w4 = 0
          if (csw(0,2) /= 0) then
            w2 = (wp2(i1)+wp2(max(i1-1,1)))/2
          endif
          if (csw(0,3) /= 0) then
            w3 = (wp3(i1)+wp3(max(i1-1,1)))/2
          endif
          if (csw(0,4) /= 0) then
            w4 = (wp4(i1)+wp4(max(i1-1,1)))/2
          endif
          call dpcopy(cs,cspt,1,4,1-w1-w2)
          call dpadd(cspt,csw(0,1),1,4,w1)
          call dpadd(cspt,csw(0,2),1,4,w2)
          call dpadd(cspt,csw(0,3),1,4,w3)
          call dpadd(cspt,csw(0,4),1,4,w4)
          cspt(0) = 15
        endif
        if (i1 <= np) call cdrw(xp(i1),yp(i1),xold,yold,lclip,box,cspt)
        if (csw(0,1) /= 0 .and. ipr > 90) print 421, xp(i1),yp(i1),cspt
  421   format(' drawing x,y',2f12.6,'  col=',4f8.3)
C        if (csw(0,1) /= 0) then
C          call mve(xp(i1),yp(i1))
C        endif
        if (logx .or. logy .or. lrot .or. lclip .or. csw(0,1) /= 0) then
          do  ip = i1+1, np-1
            if (csw(0,1) /= 0) then
              w1 = (wp(ip)+wp(max(ip-1,1)))/2
              w2 = 0
              if (csw(0,2) /= 0) then
                w2 = (wp2(ip)+wp2(max(ip-1,1)))/2
              endif
              if (csw(0,3) /= 0) then
                w3 = (wp3(ip)+wp3(max(ip-1,1)))/2
              endif
              if (csw(0,4) /= 0) then
                w4 = (wp4(ip)+wp4(max(ip-1,1)))/2
              endif
              call dpcopy(cs,cspt,1,4,1-w1-w2)
              call dpadd(cspt,csw(0,1),1,4,w1)
              call dpadd(cspt,csw(0,2),1,4,w2)
              call dpadd(cspt,csw(0,3),1,4,w3)
              call dpadd(cspt,csw(0,4),1,4,w4)
              cspt(0) = 15
            endif
            if (csw(0,1) /= 0 .and. ipr > 90)
     .        print 421, xp(ip),yp(ip),cspt
            call cdrw(xp(ip),yp(ip),xold,yold,lclip,box,cspt)
C            if (csw(0,1) /= 0) then
C              call mve(xp(ip),yp(ip))
C            endif
          enddo
        else
C     ... Above is more generic, but this way is faster
          do  3  ip = i1+1, np-1
    3     write(stream,fmt)  (offxu+xp(ip))*scalxu,(offyu+yp(ip))*scalyu
        endif
C       Last point
        if (i1+1 <= np) then
          if (csw(0,1) /= 0) then
            w1 = (wp(ip)+wp(max(ip-1,1)))/2
            w2 = 0
            if (csw(0,2) /= 0) then
              w2 = (wp2(ip)+wp2(max(ip-1,1)))/2
            endif
            if (csw(0,3) /= 0) then
              w3 = (wp3(ip)+wp3(max(ip-1,1)))/2
            endif
            if (csw(0,4) /= 0) then
              w4 = (wp4(ip)+wp4(max(ip-1,1)))/2
            endif
            call dpcopy(cs,cspt,1,4,1-w1-w2)
            call dpadd(cspt,csw(0,1),1,4,w1)
            call dpadd(cspt,csw(0,2),1,4,w2)
            call dpadd(cspt,csw(0,3),1,4,w3)
            call dpadd(cspt,csw(0,4),1,4,w4)
            cspt(0) = 15
          endif
          if (csw(0,1) /= 0 .and. ipr > 90)
     .      print 421, xp(np),yp(np),cspt
          call cdrw(xp(np),yp(np),xold,yold,lclip,box,cspt)
        endif
      endif
      goto 20

C --- Generic plcrv ---
  100 continue
C In this case, end does not matter
      xi = xp(1)
      yi = yp(1)
      if (ipr >= 40) call awrit2('    ... start line %1,6;6d %1,6;6d',
     .  ' ',80,i1mach(2),xi,yi)
      call mve(xi,yi)
C      if (xi >= xmin .and. xi <= xmax .and.
C     .    yi >= ymin .and. yi <= ymax) then
C        berr = .false.
C      else
C        berr = .true.
C      endif
      do  10  ip = 2, np
C cludge for now...
        xi = xp(ip)
        yi = yp(ip)
        if (xi >= xmin .and. xi <= xmax .and.
     .      yi >= ymin .and. yi <= ymax) then

          if (lxbrk) then
            if (xp(ip) > xp(ip-1)) then
              call drw(xp(ip),yp(ip))
            else
              if (ipr >= 40) call awrit1(
     .          '    ... break at point %i',' ',80,i1mach(2),ip)
              call mve(xp(ip),yp(ip))
            endif
          else
            call drw(xp(ip),yp(ip))
          endif
          if (ipr >= 50) print 345, 'draw within boundary',xi,yi
  345     format(' plcrv: ',a20,2f12.6)
        else
          call mve(xi,yi)
          if (ipr >= 50) print 345, 'mve outside boundary',xi,yi
        endif
   10 continue

C --- Close the curve ---
   20 continue
      if (mod(iclos,10) == 1) then
        is = 0
        s = ' '
        call plcrve(medium,pstat,csl,s,is)
        if (medium == 3) call bin2a('grestore',1,0,0,1,0,120,s,is)
        if (is > 0) write(stream,'(/a)') s(1:is)
      endif

      end
      subroutine plsym(xp,yp,zp,np,stype,syma,symcs)
C- Draw symbol at a collection of (x,y) pairs
C ----------------------------------------------------------------------
Ci Inputs
Ci   xp,yp (zp) x-y pairs where to draw symbols.
Ci   stype: compound integer of 1+10's+100 digit, 100's digit, and 1000's digit:
Ci          1+10's digit for symbol type: 1-12 for:
Ci                 x, square, diamond, +, polygon, circle, arrow,
Ci                 errbar, timelin, hist, row, wiggle
Ci          100's  digit = j>0: use zp in place of syma(j)
Ci         1000's  digit 1 => plot only symbols inside (plul,plur,plub,plut)
Ci   symcs:    color for symbol
Ci          (0) is used by plcrve and specifies fill and color
Ci              0 => no fill
Ci              1 => fill no color
Ci              2 => fill gray
Ci              3 => fill color
Ci              4 => no fill gray?
Ci              5 => no fill color?
Ci          .. Whether it is in color or B&W:
Ci             symcs(1:3)  specify the RGB fill color
Ci             symcs(4)    specifies the line boldness for symbol lines
Ci             symcs(4)=-1  => take bold value from -lt specification
Ci             symcs(4)=-2  => no strokes for symbol; only fill
Ci             symcs(4)>100 => Use symcs(5:7) for the RGB color of strokes
Ci             symcs(5:7)      Optional RGB color of the strokes
Ci   syma      Parameters defining attributes of symbol given type
Ci             Unless stated below, syma(1) is symbol-dependent size
Ci             Symbol-dependent attributes
Ci      stype    syma(i)
Ci     1  x        (1): width
Ci                 (2): height
Ci     2  square   (1): width
Ci                 (2): height
Ci     3  diamond  (1): width
Ci                 (2): height
Ci     4  +        (1): width
Ci                 (2): height
Ci     5  polygon  (1): size
Ci                 (2): number of sides
Ci                 (3): angle
Ci     6  circle   (1): size
Ci                 (2): ?
Ci     7  arrow    tip at (xp,yp), and syma:
Ci                 (1,2): tail coordinates relative to tip (GU)
Ci                 (3): head length (fraction of arrow length)
Ci                 (4): head angle rel to tail (degrees)
Ci                 (5): head length along axis (fraction of arrow length)
Ci                 (6): 0 tip  at x(ip),y(ip)
Ci                    : 1 tail at x(ip),y(ip)
Ci                    : 1/2 arrow centered at x(ip),y(ip)
Ci     8  errbar:  (1): width
Ci                 (2): scale factor for height (given by zp)
Ci     9  timelin: zp  bar length (UU)
Ci                 (1) endbar height (UU)
Ci                 (2) bar thickness
Ci                 (3) endbar thickness
Ci    10  hist:    (1) width of a histogram bar
Ci                 (2) bar thickness
Ci    11  row:     row index converted to number as : syma(1)*row + syma(2)
Ci                 (1) scale of row index
Ci                 (2) offset added to row index
Ci                 (3) centering ... should be 0..4
Ci                 (4) xshft
Ci                 (4) yshft
Ci    12 wiggle:   tip at (xp,yp), and syma:
Ci                 (1,2): tail coordinates relative to head (GU)
Ci                 (3): number of periods
Ci                 (4): excursion of wiggle about line
Ci                 (5): number of points in figure (0 => determined internally)
Cr Remarks
Cu Updates
Cu   29 Oct 16 (3.50) symcs(4)=-2 => no line around symbol; circle symbol extended to partial arc
Cu   25 Jan 15 Symbols have two color sets: one for boundary lines and one for fill
Cu   28 Aug 12 1000s digit
Cu   10 Mar 11 argument 6 for arrow
Cu   13 May 04 Extra option in errbar
Cu   05 Jun 02 Added histogram symbol
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer stype,np
      double precision syma(8),symcs(0:7),xp(np),yp(np),zp(np)

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
      double precision x,y,xold,yold
      integer nsym
      parameter (nsym=12)
      integer nb,ib,ipr,ip,is,nsyma(nsym),type,msyma,ttype,lclip
      double precision tpi,phi0,radx,rady,d,dx,dy,hsyma,xi,yi,dlength,ct,st,symcsl(0:3)
      character*80 s,slabl,nfmt*20,spos*2
      character*7 synam(nsym),nam*(*)
      data synam /'x','square','diamond','+','polygon','circle',
     .  'arrow','errbar','timelin','hist','row','wiggle'/
      data nsyma /1,1,1,1,2,1,5,2,3,2,1,5/
      data tpi /6.283185307179586d0/

      goto 1

      entry plsyn(nam,ttype)
C- return index to name of symbol
      call tokmat(nam,synam,nsym,len(nam),' ',ttype,is,.false.)
      ttype = ttype+1
      return

C --- Setup ---
    1 continue
      slabl = ' '
      nfmt = '%d'
      type = mod(stype,100)
      if (type == 0) return
      msyma = mod(stype/100,10)
      if (msyma > 0) hsyma = syma(msyma)
      lclip = mod(stype/1000,10)
      call getpr(ipr)
      is = 0
      s = ' '
      call plcrve(medium,pstat,[-2d0,0d0,0d0,0d0],s,is)  !  Close any prior curve
      if (is > 0) write(stream,*) s(1:is)
      call awrit5('plsym: '//synam(type)//' %i pts  gs=%1;5d  attr=%'
     .  //char(ichar('0')+nsyma(type))//':1;5d  x,y(1)=%1;4g %1;4g',
     .  s,80,0,np,symcs(1),syma,xp,yp)
      call skpblb(s,79,is)
      if (ipr >= 40) print *, s(1:is+1)
      if (medium == 3)
     .  write (stream,'(/a)') '% --- '//s(1:is+1)//' ---'

      call plntys
      if (symcs(4) >= 0) ltbld = nint(symcs(4))  ! bold locally specified
      if (ltbld >= 100) ltbld = ltbld - 100      ! 100s digit is a flag
      call plntyp(1,ltbld,0d0,0d0,0d0,0d0)
      if (symcs(4) >= 10 .and. medium == 3) then
        call awrit1(' gsave%3:1;4d setrgbcolor',s,80,stream,symcs(5))
      endif

C...  Clip beyond plot limits
C      if (pstat == 0) then
C        xsmin = (offxu+xmin)*scalxu
C        ysmin = (offyu+ymin)*scalyu
C        ysmax = (offyu+ymax)*scalyu
C        xsmax = (offxu+xmax)*scalxu
C        if (logx) then
C          xsmin = (offxu+dlog(xmin))*scalxu
C          xsmax = (offxu+dlog(xmax))*scalxu
C        endif
C        if (logy) then
C          ysmin = (offyu+dlog(ymin))*scalyu
C          ysmax = (offyu+dlog(ymax))*scalyu
C        endif
C        call pbbox(xsmax,ysmax)
C        call pbbox(xsmin,ysmin)
C        call awrit8('gsave newpath %1;nd %1;nd moveto %1;nd %1;nd '//
C     .    'lineto',' ',120,stream,ndwrt,xsmin,ndwrt,ysmin,ndwrt,xsmin,
C     .    ndwrt,ysmax)
C        call awrit8('              %1;nd %1;nd lineto %1;nd %1;nd '//
C     .    'lineto closepath clip',' ',120,stream,ndwrt,xsmax,
C     .    ndwrt,ysmax,ndwrt,xsmax,ndwrt,ysmin)
C      endif

C --- For each point, do ----
      do  2  ip = 1, np

        if (lclip /= 0) then
          if (xp(ip) < plul .or. xp(ip) > plur) cycle
          if (yp(ip) < plub .or. yp(ip) > plut) cycle
        endif

        if (msyma > 0) syma(msyma) = zp(ip)
        call mve(xp(ip),yp(ip))
        x = .015d0*syma(1)*scalxg/scalxu
        if (syma(2) /= 0 .and. iabs(type) <= 4) then
          y = .015d0*syma(2)*scalyg/scalyu
        else
          y = .015d0*syma(1)*scalyg/scalyu
        endif
        goto (10,20,30,40,50,60,70,80,90,100,110,120), iabs(type)
        goto 5

C --- Symbol type cross ---
   10   continue
        call drwr(x,y)
        call drwr(-2*x,-2*y)
        call drwr(x,y)
        call drwr(-x,y)
        call drwr(2*x,-2*y)
        call drwr(-x,y)
        goto 5

C --- Symbol type square (rectangle) ---
   20   continue
        call mver(x,y)
        call drwr(-2*x,0d0)
        call drwr(0d0,-2*y)
        call drwr(+2*x,0d0)
        call drwr(0d0,+2*y)
        goto 5

C --- Symbol type diamond ---
   30   continue
        call mver(0d0,y)
        call drwr(-x,-y)
        call drwr(+x,-y)
        call drwr(+x,+y)
        call drwr(-x,+y)
        goto 5

C --- Symbol type +  ---
   40   call drwr(0d0,y)
        call drwr(0d0,-2*y)
        call drwr(0d0,y)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        goto 5

C --- Symbol type poly: syma(1..3): radius, no. sides, angle 1st pt(rad)
   50   continue
        radx = dsqrt(2d0)*x
        rady = dsqrt(2d0)*y
        nb = nint(syma(2))
        phi0 = tpi/(2*nb) + syma(3)
        x = radx*dsin(phi0)
        y =-rady*dcos(phi0)
        call mver(x,y)
        do  55  ib = 1, nb
          xold = x
          yold = y
          x = radx*dsin(ib*tpi/nb+phi0)
          y =-rady*dcos(ib*tpi/nb+phi0)
          call drwr(x-xold,y-yold)
   55   continue
        goto 5

C --- Symbol type circle, radius a. PostScript only ---
   60   continue
        if (medium == 1) then
          syma(2) = 8
          goto 50
        endif
        if (medium /= 3) return
        call awrit5(' newpath%;8,2D%;8,2D%;8,2D %d %d arc',s,80,stream,
     .    xcm,ycm,syma(1)*6,syma(2),syma(3))
C        write (stream,530) xcm,ycm,syma(1)*6
        pstat = 1
        goto 5

C --- Symbol type arrow ---
   70   continue
C ...   Draw tail
        xi = xp(ip) - syma(6)*syma(1)*(plur-plul)
        yi = yp(ip) - syma(6)*syma(2)*(plut-plub)
        if (xi /= xp(ip) .or. yi /= yp(ip)) call mve(xi,yi)
        call mver(syma(1)*(plur-plul),syma(2)*(plut-plub))
        call drwe(xi,yi,[-2d0,0d0,0d0,0d0])
        if (syma(3) == 0) goto 5
C ...   Position of three arrow points, UU. phi0 is arrow direction.
        dx = syma(1)*(plur-plul)*scalxu
        dy = syma(2)*(plut-plub)*scalyu
        phi0 = datan2(dy,dx)
        d = dsqrt(dx**2+dy**2)
        radx = phi0 + tpi*syma(4)/360
        call drw(d*syma(3)*dcos(radx)/scalxu + xi,
     .           d*syma(3)*dsin(radx)/scalyu + yi)
        call drw(d*syma(5)*dcos(phi0)/scalxu + xi,
     .           d*syma(5)*dsin(phi0)/scalyu + yi)
        radx = phi0 - tpi*syma(4)/360
        call drw(d*syma(3)*dcos(radx)/scalxu + xi,
     .           d*syma(3)*dsin(radx)/scalyu + yi)
        call drw(xi,yi)
        goto 5

C --- Symbol type errbar ---
   80   continue
        d = yp(ip)+zp(ip)*(1+syma(2))
        if (logy) then
          if (d <= 0) d = yp(ip)
        endif
        call drw(xp(ip),d)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        d = yp(ip)-zp(ip)*(1-syma(2))
        if (logy) then
          d = yp(ip)**2/(yp(ip)+(zp(ip))*(1-syma(2)))
          if (d <= 0) d = yp(ip)
        endif
        call drw(xp(ip),d)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        goto 5

C --- Symbol type timelin ---
   90   continue
        call plntyp(1,nint(syma(3)),0d0,0d0,0d0,0d0)
        call drwr(0d0,syma(1)/2)
        call drwr(0d0,-syma(1))
        call drwre(0d0,syma(1)/2,symcs)
        call plntyp(1,nint(syma(2)),0d0,0d0,0d0,0d0)
        call drwre(zp(ip),0d0,symcs)
        call plntyp(1,nint(syma(3)),0d0,0d0,0d0,0d0)
        call drwr(0d0,syma(1)/2)
        call drwr(0d0,-syma(1))
        call drwre(0d0,syma(1)/2,symcs)
        goto 5

C --- Symbol type histogram ---
  100   continue
        call mve(xp(ip),syma(2))
C       call drwr(-syma(1)/2*(plur-plul),0)
        call drwr(-syma(1)/2,0d0)
        call drwr(0d0,yp(ip)-syma(2))
C       call drwr(syma(1)*(plur-plul),0)
        call drwr(syma(1),0d0)
        call drwr(0d0,syma(2)-yp(ip))
C       call drwr(-syma(1)/2*(plur-plul),0d0)
        call drwr(-syma(1)/2,0d0)
        goto 5

C --- Symbol is ascii reps'n of row index ---
  110   continue
        is = 0
        d = syma(1)*dble(ip) + syma(2)
        spos = 'tx'
        if (syma(3) >= 0 .and. syma(3) <= 4)
     .    call awrit1('%i',spos,2,0,nint(syma(3)))
        call mver(syma(4)*(plur-plul),syma(5)*(plut-plub))
        call pslabl(spos,nfmt,d,' ','l',' ',1,0d0)
        goto 5

C --- Symbol type wiggly line ---
  120   continue
C   ... Number of points to compose the wiggly line
        nb = nint(syma(5)); if (nb == 0) nb = nint(12*syma(3))
C   ... Draw the curve
        call mve(xp(ip),yp(ip))
        do  is = 1, nb

C         Straight line
C         xi = dble(is-1)/(nb-1)*syma(1)*(plur-plul)
C         yi = dble(is-1)/(nb-1)*syma(2)*(plut-plub)
C         print *, xi, yi

          dx = dble(is-1)/(nb-1) ! unrotated, unnormalized
          dy = syma(4)*(plur-plul)*dsin(tpi*syma(3)*dble(is-1)/(nb-1))
          ct = syma(1)*(plur-plul)
          st = syma(2)*(plut-plub)
          xi = dx*ct - dy*st
          yi = dx*st + dy*ct
C          print *, xi, yi
C          print *, dx, dy
C          print *
          call drw(xp(ip)+xi,yp(ip)+yi)
        enddo
C       call drwre(0d0,0d0)
        goto 5

C --- Cleanup for this point ---
    5   continue
        is = 0
        s = ' '
        symcsl = symcs(0:3)
        if (symcs(4) == -2) symcsl(0) = symcsl(0)+100
        call plcrve(medium,pstat,symcsl,s,is)
        if (is > 0) write(stream,*) s(1:is)
    2 continue
      if (symcs(4) >= 10 .and. medium == 3) then
        call awrit0(' grestore',s,80,stream)
      endif
      call plntyg
      if (msyma > 0) syma(msyma) = hsyma

  530 format(' newpath',3f8.2,' 0 360 arc')
      end
      subroutine plcrve(medium,pstat,cs,s,is)
C- Returns string signaling termination of a curve; sets pstat to 0
C ----------------------------------------------------------------------
Ci Inputs
Ci   medium
Ci   pstat  If pstat is zero, nothing is done
Ci   cs     color scale.  0 flags filling (see Outputs) 1..3 are colors
Co Outputs
Co   pstat is set to 0
Co   Appended to string s (Postscript mode):
Co   cs(0)<1:  s = 'stroke'                    (draw a line)
Co   cs(0)=1:  s = 'closepath stroke'
Co   cs(0)=2:  s = 'closepath gsave # setgray fill grestore stroke'
Co   cs(0)=3:  s = 'closepath gsave ### setrgbcolor fill grestore stroke'
Co   cs(0)=4:  s = '# setgray stroke'
Co   cs(0)=5:  s = '### setrgbcolor stroke'    (draw a colored line)
Co   cs(0)=6:  s = '# setgray closepath stroke'
Co   cs(0)=7:  s = '### setrgbcolor closepath stroke'
Co   cs(0)=8:  s = '### setrgbcolor'           (set up color)
Co   ... add 100 to cs0 => remove final stroke (fill without line)
Cu Updates
Cu   13 Feb 01 revised what plcrve does depending on cs(0)
C ----------------------------------------------------------------------
      double precision cs(0:3)
      integer medium,pstat,is,ls
      character s*(*)
      double precision cs0
      logical lstroke

      if (pstat /= 0 .and. medium == 3) then

C        if (cs(0) >= 10) then
C          print *, 'plcrve',cs(0)
C        endif
        cs0 = mod(cs(0),10d0)
        lstroke = cs(0) < 100
        if (cs0 >= 7) then
          ls = len(s)
          call awrit1(' %3:1;4d setrgbcolor closepath',
     .      s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 6) then
          ls = len(s)
          call awrit1(' %1;4d setgray closepath',s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 5) then
          ls = len(s)
          call awrit1(' %3:1;4d setrgbcolor ',s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 4) then
          ls = len(s)
          call awrit1(' %1;4d setgray ',s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 3) then
          ls = len(s)
          call awrit1(' closepath gsave%3:1;4d setrgbcolor '//
     .      'fill grestore',s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 2) then
          ls = len(s)
          call awrit1(' closepath gsave %1;4d setgray fill grestore',
     .      s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 >= 1) then
          call bin2a('closepath',1,0,0,1,0,len(s),s,is)
        endif
        if (lstroke) call bin2a('stroke',1,0,0,1,0,len(s),s,is)
      endif
      pstat = 0
      end
      subroutine key(style,xk,yk,xlen,xmin,xmax,ymin,ymax,
     .  symtyp,syma,symcs,crvcs,legend,blank)
C- Draw the key's legend at (xk,yk)
C ----------------------------------------------------------------------
Ci Inputs: xk,yk (user units), xlen: length of line to draw, syma,b,typ:
Ci         style: sets style for legend (See Remarks)
Ci                0 let this routine choose legend
Ci                1 or 2: use this style for legend
Ci         symbol, legend: legend
Cr Remarks
Cr   Draws legend containing symbol s. drawn in one of these styles:
Cr                 lt and s         lt only        s only
Cr     style=1:  "---- s ----"    "---------"    "   s    "
Cr     style=2:  "s"              "."            "s"
Cr     style=3:  "s---------s"    "---------"    "s      s"
Cr     style=4:  "s         s"    "         "    "s      s"
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      double precision xk,yk,xlen,syma(8),symcs(0:7),crvcs(0:3,2)
      integer style,symtyp,blank
      character*60 legend

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

C Local Parameters
      double precision x(2),y(2),xmax,xmin,ymax,ymin
      double precision xm,ym,csw(0:3,4)
      integer i1mach,isw,lstyle
      character*80 s

      lstyle = style
      if (style == 0) lstyle = 1
      if (style == 0 .and. lintyp == 0) lstyle = 2
      call dpzero(csw,4*4)

      s = ' '
      call awrit5(' key:  x,y=%1;4g,%1;4g(UU)  lt %i  sym %i  len=%1;4g'
     .  ,s,80,0,xk,yk,lintyp,symtyp,xlen)
      s = trim(s) // '  legend='//legend
C     call awrit1('%a  legend='//legend,s,80,0,xk)
      call awrit1('%a',s,80,-i1mach(2),xk)
      if (medium == 3) call awrit1(' %% KEY%a',s,80,-stream,xk)
      x(1) = xk
      x(2) = xk + xlen
      if (logx) x(2) = xk*dexp(xlen)
      y(1) = yk
      y(2) = yk
      if (lstyle == 2) x(2) = x(1)
      if (lstyle /= 4 .or. symtyp == 0) then
      call plcrv(x,y,x,x,x,x,2,xmin,xmax,ymin,ymax,0,crvcs,csw,11)
      endif
      if (lstyle < 3) then
        if (.not. logx) x(1) = (x(1)+x(2))/2
        if (      logx) x(1) = dsqrt(x(1)*x(2))
      endif
      call plsym(x,y,[0d0],1,symtyp,syma,symcs)
      if (lstyle == 3 .or. lstyle == 4) then
        call plsym(x(2),y(2),[0d0],1,symtyp,syma,symcs)
      endif
C      if (lintyp == 0) then
C        x(1) = x(1) + xlen/10
C        x(2) = x(2) - xlen/10
C        call plsym(x,y,0d0,2,symtyp,syma,symcs)
C        x(1) = x(1) - xlen/10
C        x(2) = x(2) + xlen/10
C      else
C        if (.not. logx) x(1) = (x(1)+x(2))/2
C        if (      logx) x(1) = dsqrt(x(1)*x(2))
C        call plsym(x,y,0d0,1,symtyp,syma,symcs)
C      endif
      if (medium == 3) then
        call mve(x(2),y(2))
        call pltam(xm,ym,-1,[1d0,0d0,0d0,0d0])
        call awrit3(' %1;1d %1;1d moveto h %?#n#.8#.5# mul h .5 mul neg'
     .    //' rmoveto',' ',120,stream,xm,ym,isw(x(2) == x(1)))
C       write (stream,030) xm,ym
        call pslabl('1',' ',0d0,legend,'h','f',blank,0d0)
      endif
  030 format(2f8.1,' moveto h .8 mul h .5 mul neg rmoveto')
      end

      subroutine pslabl(shift,nfmt,val,strn,hvl,fr,blank,rots)
C- Put a PostScript label or string at the current postion
C ----------------------------------------------------------------------
Ci Inputs:
Ci    shift: one of the following
Ci           'a','o' for abscissa, ordinate labelling,
Ci                   In this case, val is converted into a string
Ci           '0'-'4' for centring or placing in quadrants 1-4 about
Ci                   the current position
Ci           'tx'    for simple text 'show'
Ci           'cc'    or any of the nine pairs of characters
Ci                   lu lc ld  cu cc cd  ru rc rd
Ci                   1st character justfies left, center, or right
Ci                   2nd character justfies top, middle, or bottom
Ci    nfmt:  awrite format for val (see hvl='l')
Ci    val:   (hvl='l'): value to be converted into label
Ci    strn:  the string to print, according to hvl. A '\' in the string
Ci           forces a line break. Multiple lines are written forwards
Ci           starting at the current point or reversed ending at the
Ci           current point, according to fr
Ci    hvl:   if 'l', then convert val to a string, using awrite format.
Ci           if 'h' or 'v', label is horizontal or vertical
Ci    fr:    'f' or 'r' for whether multiple lines are written forward
Cr           or reversed
Ci    blank: 1 => blank out area to be written over
Ci    rots:  rotation of string.
Cr Remarks
Cr  pslabl writes a  PostScript segment assuming currentpoint has been
Cr  previously set. It does NOT take the current position from COMMON
Cb Bugs
Cb   This is a hoplessly confusing routine.  Plans are to move some
Cb   of its functionality to pstr.
Cu Updates
Cu   31 Jan 15 Enable sequence {\..} ... gets written as \..
Cu   30 Apr 07 bug fix: rots passed to fwrite
Cu   25 Apr 07 hvl='l' translates fortran aEn format to a.10^n
C ----------------------------------------------------------------------
      implicit none
C Passed Paramters
      character*1 shift*(*),hvl,fr,strn(0:*)
      character(len=*) :: nfmt ! this is slightly diff from the original but shall work fine
      integer blank
      double precision val,rots
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .     plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .     plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .     cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
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

C local parameters
      integer maxlen,maxlns,isw,i1,i2
      parameter (maxlns=10,maxlen=120)
      character*(maxlen) label,show*6,newl*9,s,ss,line(maxlns),strn2,bs*1
      integer lbrk(maxlns),nlines,ich,count,il,ip,ilen,l1,l2,l3,
     .        icop,llen(maxlns),iprint,i1mach,i,j
      logical cmdopt,flg,lshft

      call spchar(1,bs)

      show = 'show'//shift
      lshft = .false.
      if (len(shift)>1) then
        if (shift(1:2) == 'tx') lshft = .true.
      endif
      if (lshft) show = 'show'
      newl = 'newline'//hvl//fr
      if (hvl == 'l') then
        ip = 0
        count = 0
        label = ' '
C       if (ndec > 4) stop 'PSLABL: more than four decimal places !?'
C       call bin2a('(f10.4)',0,ndec,val,4,count,maxlen,label,ip)
C       call bin2a('d4:0',0,ndec,val,4,count,maxlen,label,ip)
        call bin2a0(0)
        call awrit1(nfmt,label,maxlen,0,val)

C  ...  Convert fortran 'e' format to 10^... and use fwrite
        call chrps2(label,'eEdD',4,len(label),ip,i2)
        if (i2 /= 0) then
          call word(label,1,i1,i2)
          ss = label
          ss = label(1:ip) // '^{.}10^{' // label(ip+2:i2) // '}'
          call word(ss,1,i1,i2)
C         Shorten 1.10^... to 10^...
          if (ss(i1:i1+4) == '1^{.}') i1=i1+5
          if (ss(i1:i1+5) == '-1^{.}') then
            i1=i1+5
            ss(i1:i1) = '-'
          endif
          write (stream,200)
  200     format('gsave /font /Symbol def fontsize font choosefont')
          if (shift == 'a'.or. shift == 'A') then
            write (stream,030) xcm,ycm
            call fwrite(stream,blank,'h',0d0,'showcD','newlinevf',
     .        ss(i1:i2),i2-i1+1,label)
          else if (shift == 'o' .or. shift == 'O') then
            write (stream,030) xcm,ycm,' h hshift mul neg 0 rmoveto'
            call fwrite(stream,blank,'h',0d0,'showlc','newlinevf',
     .        ss(i1:i2),i2-i1+1,label)
          endif
          write (stream,300)
  300     format ('grestore')
C  ...  Just use standard awrite
        else
        call skpblb(label,maxlen,ilen)
        ilen = ilen + 1
        if (iprint() > 50) call awrit2(
     .    ' pslabl: label "'//label(1:ilen)//'" at %1;2d,%1;2d (MU)',
     .    s,maxlen,i1mach(2),xcm,ycm)
        call pbbox(xcm,ycm)
        write (ss,031) xcm,ycm,label(1:ilen),show
        call awrit3('%xgsave blank-%a%?#n#on#off#%?#n# %d rotate#%j#',
     .    s,maxlen,0,isw(blank /= 0),isw(rots /= 0),rots)
C       s = 'gsave blank-on'
C       if (blank == 0) s = 'gsave blank-off'
        call awrit1('%a '//ss,s,maxlen,0,1)
        call awrit1('%a',s,maxlen,-stream,1)
        endif
      else
        if (hvl == 'v') then
          if (shift == '1') show = 'show2'
          if (shift == '2') show = 'show3'
          if (shift == '3') show = 'show4'
          if (shift == '4') show = 'show1'
        endif
        call skpblb(strn,maxlen,ilen)
        if (ilen < 0) return
        nlines = 1
        ich = 0
        do  1  il = 1, maxlns
C         call chrpos(strn,'\\',ilen,ich)
          call chrpos(strn,bs,ilen,ich)
          if (ich > 1) then
            if (strn(ich-1) == '{') cycle
          endif
          lbrk(il) = ich
          if (ich == ilen) goto 2
          ich = ich + 1
          nlines = nlines + 1
    1   continue
    2   continue
        ip = 0
        do  3  il = 1, nlines
          if (nlines == 1) then
            icop = ilen+1
          else
            if (il == 1) then
              icop = lbrk(1)
            elseif (il == nlines) then
              icop = lbrk(nlines) - lbrk(nlines-1)
            else
              icop = lbrk(il)-1 - lbrk(il-1)
            endif
          endif
!         call strcop(line(il),strn(ip),icop,bs,ich)
          j = 1; flg = .false.
          do i = 1, icop
            line(il)(j:j) = strn(ip)
            if (j > 1 .and. line(il)(max(j-1,1):j) == '{'//bs) then  ! Handle '{\' sequence
              line(il)(j-1:j-1) = bs
              flg = .true.
            elseif (flg .and. line(il)(j:j) == '}') then ! end of sequence
              flg = .false.
            else
              j = j+1
            endif
            ip = ip + 1
          enddo
          ip = ip+1
          llen(il) = j-1
          if (iprint() > 50) write (*,100) line(il)(1:llen(il))
    3   continue
        if (fr == 'f') then
          l1 = 1
          l2 = nlines
          l3 = 1
        else
          l1 = nlines
          l2 = 1
          l3 = -1
        endif
        do  4  il = l1, l2, l3
          s = 'gsave blank-'
          call awrit3('%a%?#n#on#off#%?#n# %d rotate#%j#',s,maxlen,0,
     .      blank /= 0,isw(rots /= 0),rots)
          if (cmdopt('-plaintext ',11,0,strn2) .or. lshft) then
            if (hvl == 'h') then
C             write (stream,130) line(il)(1:llen(il)),show,newl
              write (ss,131) line(il)(1:llen(il)),show,newl
            else
C             write (stream,230) line(il)(1:llen(il)),show,newl
              write (ss,231) line(il)(1:llen(il)),show,newl
            endif
            call awrit1('%a '//ss,s,maxlen,0,1)
            call awrit1('%a',s,maxlen,-stream,1)
          else
            call fwrite(stream,blank,hvl,rots,show,newl,
     .                  line(il)(1:llen(il)),llen(il),strn2)
          endif
    4   continue
      endif
  030 format(2f8.1,' moveto':/a)
  031 format(2f8.1,' moveto (',a,') ',a6,' grestore')
  130 format('gsave (',a,') ',a6,' grestore ',a9)
  131 format('(',a,') ',a6,' grestore ',a9)
  230 format('gsave 90 rotate (',a,') ',a6,' grestore ',a9)
  231 format('90 rotate (',a,') ',a6,' grestore ',a9)
  100 format(' PSLABL new line: ',a)
      end

      subroutine fwrite(stream,blank,hvl,rots,show,newl,line,ilen,ss)
C- Generalised write using FORTRAN
C ----------------------------------------------------------------------
Ci Inputs:
Ci   stream
Ci   blank: 1 => blank out area to be written on
Ci   hvl:   if 'v', rotate string 90 degrees; otherwise:
Ci   rots:  if nonzero, rotate string by angle
Ci   show:  same as shift in pslabl EXCEPT these options are not allowed
Ci          'a','o'
Ci          'tx'    for simple 'show'
Ci   newl:  controls justification in newlines and should be one of:
Ci          newlinehf, newlinevf, newlinehr, newlinevr
Ci   line:string to be written; ilen its length (including control chars)
Ci   ss: work string
Co Outputs:
Co
Cr Remarks
Cr   fplot will accept a string with control characters ^ _ ~ @ & to
Cr   make super- and sub-scripts,and  greek, bold or italic fonts.
Cr   The substring to be converted is enclosed in { .. }.
Cr   The following character pairs act as flages that cause fwrite to
Cr   change fonts for these strings:
Cr      superscript   subscript     Greek    bold  italic
Cr         ^{            _{           ~{      @{     &{
Cr   The conversion stops at when } is encountered
Cr   For example
Cr   ~{m}_{0}=4~{p}.10^{7} is equivalent to
Cr   $\mu$_0=4$\pi$.10$^7$ in TeX. The { .. } may contain any number
Cr   of characters and may also be nested, say, to make sub-sub-scripts.
Cu Updates
Cu   23 Oct 09 Added italics (& key)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ilen,maxlen,stream,blank
      parameter (maxlen=120)
      character*(maxlen) show*6,newl*9,hvl*1,ss
      character*(ilen) line
      double precision rots
C ... Local parameters
      integer i,tlen,iw
      character c*1,cp*1,cn*1,cw*(maxlen)
      character bkslsh*1

      call spchar(1,bkslsh)

C --- strip line of its control chars to determine true stringwidth ---
      tlen = 0
      do  i = 1, ilen
        c = line(i:i)
        if (c /= '^' .and. c /= '_' .and. c /= '~' .and.
     .      c /= '@' .and. c /= '&' .and.
     .      c /= '{' .and. c /= '}') then
          tlen = tlen + 1
          ss(tlen:tlen) = c
        endif
      enddo

C --- save the graphics state ---
      write (stream,100)

C --- inform PostScript of the true stringwidth ---
      write (stream,102) ss(1:tlen)

C --- Rotate text (do after previous write to avoid pdf bug) ---
      if (hvl == 'v') then
        write (stream,101) 90d0
      elseif (rots /= 0) then
        write (stream,101) rots
      endif

C --- align text according to right/left up/down rules ---
      if (show == 'showru' .or. show == 'show1') then
      elseif (show == 'showcu') then
        write (stream,104)
      elseif (show == 'showlu' .or. show == 'show4') then
        write (stream,105)
      elseif (show == 'showrc') then
        write (stream,106)
      elseif (show == 'showcc' .or. show == 'show0') then
        write (stream,107)
      elseif (show == 'showlc') then
        write (stream,108)
      elseif (show == 'showrd' .or. show == 'show2') then
        write (stream,109)
      elseif (show == 'showld' .or. show == 'show3') then
        write (stream,110)
      elseif (show == 'showcd') then
        write (stream,111)
      elseif (show == 'showcD') then
        write (stream,112)
      endif

C --- blank out area to be written on ---
      if (blank == 1) then
        write (stream,103)
      endif

C --- generalised write ---
      cw = ' '
      iw = 0
      write (stream,121)
      do i = 1, ilen
        c = line(i:i)
        cp = ' '
        if (i > 1) cp = line(i-1:i-1)
        cn = ' '
        if (i < ilen) cn = line(i+1:i+1)
        if ((c == '^' .and. cn == '{') .or.
     .      (c == '_' .and. cn == '{') .or.
     .      (c == '~' .and. cn == '{') .or.
     .      (c == '@' .and. cn == '{') .or.
     .      (c == '&' .and. cn == '{')) then
C         Flush internal buffer
          if (iw > 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
        elseif (c == '{' .and. cp == '^') then
          write (stream,122)
        elseif (c == '{' .and. cp == '_') then
          write (stream,123)
        elseif (c == '{' .and. cp == '~') then
          write (stream,124)
        elseif (c == '{' .and. cp == '@') then
          write (stream,127)
        elseif (c == '{' .and. cp == '&') then
          write (stream,126)
        elseif (c == '}') then
C         Flush internal buffer
          if (iw > 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          if (cn /= '}') then
            write (stream,125)
          endif
        elseif (c == '(') then
C         Flush internal buffer
          if (iw > 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          write (stream,129) bkslsh,'050'
        elseif (c == ')') then
C         Flush internal buffer
          if (iw > 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          write (stream,129) bkslsh,'051'
        else
          iw = iw + 1
          cw(iw:iw) = c
        endif
      enddo
C     Flush internal buffer
      if (iw > 0) then
        write (stream,130) cw(1:iw)
        iw = 0
      endif

      write (stream,200) newl

  100 format ('gsave')
  101 format (f6.1,' rotate')
  102 format ('/w (',a,') stringwidth pop def')
  103 format ('savecurrentpoint gsave w /x exch def'/
     .        'newpath xbak ybak h tails mul sub moveto'/
     .        'x 0 rlineto 0 h 1.1 mul h tails mul add rlineto x',
     .        ' neg 0 rlineto'/'closepath 1 setgray fill grestore')
  104 format ('w 2 div neg 0 rmoveto')
  105 format ('w neg 0 rmoveto')
  106 format ('0 h 2 div neg rmoveto')
  107 format ('w 2 div neg h 2 div neg rmoveto')
  108 format ('w neg h 2 div neg rmoveto')
  109 format ('0 h neg rmoveto')
  110 format ('w neg h neg rmoveto')
  111 format ('w 2 div neg h neg rmoveto')
  112 format ('w 2 div neg h 1.7 mul neg rmoveto')
  121 format ('/fntsiznow fontsize def /accshift 0 def /frac 1 def')
  122 format ('/fntsiznow fntsiznow smaller mul def'/
     .        'fntsiznow font choosefont'/
     .        '/accshift accshift h shiftu mul frac mul add def'/
     .        '0 h shiftu mul frac mul rmoveto'/
     .        '/frac frac 0.7 mul def')
  123 format ('/fntsiznow fntsiznow smaller mul def'/
     .        'fntsiznow font choosefont'/
     .        '/accshift accshift h shiftd mul frac mul sub def'/
     .        '0 h shiftd mul frac mul neg rmoveto'/
     .        '/frac frac 0.7 mul def')
  124 format ('fontsize /Symbol choosefont')
  126 format ('fontsize /Times-Italic choosefont')
  127 format ('fontsize /Times-Bold choosefont')
  125 format ('/fntsiznow fontsize def fontsize font choosefont'/
     .        '0 accshift neg rmoveto /accshift 0 def /frac 1 def')
  129 format ('(',a1,a,') show')
  130 format ('(',a,') show')

  200 format ('grestore ',a9)
      end

      subroutine pstr(x,y,rots,string,arg1,arg2,arg3,pos,units,optio)
C- Paint a string in the current font at coordinates x,y
C ----------------------------------------------------------------------
Ci Inputs:
Ci   optio:  a compound of digits specifying options
Ci           1s digit
Ci           1  blank-on
Ci          10s digit
Ci           0  raw string
Ci           1  string is a format entering into awrit2
Ci           2  string is a format ...?
Ci         100s digit
Ci           1  string a current coordinates; thus x,y are not used
Ci   x,y:    coordinates at which to paint string
Ci   rots:   rotation angle at which to write
Ci   units:  'u' for x,y in user units, 'm' for medium (cm) units
Ci   string: string to paint, or formatting command; see optio.
Ci   arg1,arg2: optional arguments entering into string format
Ci   pos:    compound of two characters, the first one of
Ci           l,c,r followed by one of u,c,d
Cb   Bugs
Cb     10's and 100's digit of optio never tested.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      double precision x,y,rots,arg1,arg2,arg3
      character pos*2,units*1,string*(*)
      integer optio
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
      integer iprint,i1mach,isw,blank,maxlen,ifmt  ! ,is,is1,is2
      parameter (maxlen=120)
      character*(maxlen) s, strn
      character*(2) lit
      double precision xx,yy

      if (medium /= 3) return

C --- Postscript label ---
      blank = mod(optio,10)
      ifmt  = mod(optio/10,10)
      lit = '{'; call spchar(1,lit(2:2)) ! Strip brackets in "{\..}"
      if (ifmt == 0) then
        strn = string
C      print *, index(string,lit),index(string,'}')
C      is1 = index(string,lit); is2 = index(string,'}')
C      if (ifmt == 0 .and. (is1 == 0 .or. is1 >= is2)) then
C        strn = string
C      elseif (ifmt == 0) then ! Strip {} from {\..}
C        is = 1
C        strn = ' '
C   10   continue ! Re-entry for another quote
C        if (is1 > is) strn = trim(strn)//string(is:is1-1) ! characters up to strip
C        strn = trim(strn)//string(is1+1:is2-1)  ! Add string contents up to is2
C        is = is2+1
C        print *, string(is:)
C        is1 = index(string(is:),lit)+is-1; is2 = index(string(is:),'}')+is-1
C        print *, string(is1:is2)
C        if (is1 > 0 .and. is1 < is2) goto 10
C

C        do while
C        if (string(1:1) == '{')
C          is = 1
C          is1 = 1
C          call nwordg(string,10,'{0123456789}',1,is1,is2)
C          print *, string(is1:is2)
C        else
C         strn = string
C        endif

        call awrit3(' pstr: '//pos//': '//'%?#n#rot=%d #%j#'//
     .    'blk=%l  "',s,80,0,isw(rots /= 0),rots,blank /= 0)
        s = trim(s) // strn

      elseif (ifmt == 1) then
        strn = ' '
        call awrit3(string,strn,maxlen,0,arg1,arg2,arg3)
        call awrit3(' pstr: '//pos//': '//'%?#n#rot=%d #%j#'//
     .    'blk=%l  "'//strn,s,80,0,isw(rots /= 0),rots,blank /= 0)
      else
        call rxi('pstr: bad formatting switch',ifmt)
      endif
      if (mod(optio/100,10) == 0) then
        call awrit0('%a" at',s,80,0)
        xx = x
        yy = y
        if (units == 'm' .and. scalxu*scalyu /= 0) then
          xx = x/scalxu - offxu
          yy = y/scalyu - offyu
          if (logx) xx = exp(xx)
          if (logy) yy = exp(yy)
          call awrit2('%a %1;4g,%1;4g(u)',s,80,0,xx,yy)
          xx = x
          yy = y
        elseif (units /= 'm') then
          xx = (offxu+x)*scalxu
          yy = (offyu+y)*scalyu
          if (logx) xx = (offxu+dlog(x))*scalxu
          if (logy) yy = (offyu+dlog(y))*scalyu
          call rxx(scalxu*scalyu == 0,'pstr: user''s units not set')
          call awrit2('%a %1;4g,%1;4g(u)',s,80,0,x,y)

        endif
        call awrit2('%a %1,1;1d,%1,1;1d(m)',s,80,0,xx,yy)
      endif
      if (iprint() >= 30) call awrit0('%a',s,len(s),-i1mach(2))
      if (medium == 3) call awrit0('%5o%0p%% ---%a',s,len(s),-stream)
      if (mod(optio/100,10) == 0) then
        write (stream,130) xx,yy
        call pbbox(xx,yy)
  130   format(2f8.1,' moveto')
      endif
      call pslabl(pos,' ',0d0,strn,'h','f',blank,rots)
      end
      subroutine pbbox(x,y)
C- Update bbox
      implicit none
C Passed Parameters
      double precision x,y
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

C      if (x < bbox(1) .or. x > bbox(2) .or.
C     .    y < bbox(3) .or. y > bbox(4)) then
C        print *, 'x,y=',x,y
C      endif
      bbox(1) = min(bbox(1),x)
      bbox(2) = max(bbox(2),x)
      bbox(3) = min(bbox(3),y)
      bbox(4) = max(bbox(4),y)


      end

      subroutine setfnt(font,fontsz)
C- Set font in PostScript protocol
C ----------------------------------------------------------------------
Ci Inputs:
Ci   font and font_size
C ----------------------------------------------------------------------
      implicit none

      character font*19
      integer fontsz,i
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

      fsiz = fontsz
      if (medium == 3) write (stream,130) trim(font),fontsz
      call lodsyv('ch',1,.6d0*fontsz,i)
  130 format(' % FONT: set font and fontsize ..'/
     .       '/font ',a,' def /fontsize ',i3,' def'/
     .       '/h fontsize 0.6 mul def % h is about char height'/
     .       'fontsize font choosefont')
      end
      subroutine psinit(stream,irot,transx,transy,scale)
C- Initialise PostScript file with defs and defaults
C ----------------------------------------------------------------------
Ci Inputs
Ci   stream:file logical unit
Ci   transx:shift of origin from bottom left corner of paper
Ci   transy:shift of origin from bottom left corner of paper
Ci   scale :shrinks or expands entire plot
Co Outputs
Co   definitions are written to file stream
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer stream,irot
      double precision transx,transy,scale

      write(stream,030)
C --- translation, scaling, rotating ---
      call awrit4('/inch{72 mul}def /cm{72 mul 2.54 div}def' //
     . ' %1;1d %1;1d translate  %1;1d %1;1d scale',
     .  ' ',90,stream,transx, transy, scale, scale)
      if (irot /= 0) write(stream,031)
C --- fontsize and definitions ---
      write (stream,230)
      write (stream,231)
      write (stream,240)
      write (stream,241)
      write (stream,242)
      write (stream,250)

  030 format('%!PS-Adobe-2.0 generated by plsub.f ...')
  031 format(' 90 rotate 0 -612 translate')
  130 format(2f8.1,' translate ',2f8.1,' scale')
  230 format(
     . '/tails 0.4 def % set to 0 for blanking chars without tails'/
     . '/lto{lineto}def % shorthand for lineto'/
     . '/blank-on{/b 1 def}def /blank-off{/b 0 def}def blank-on'/
     . '/smaller 0.67 def % sub/superscripts are this much smaller'/
     . '/shiftu  0.55 def % are shifted up by this much of the char'/
     . '/shiftd  0.33 def % are shifted down by this much of the char',
     . 'height'/
     . '/hshift  0.10 def % horizontal space, fraction of char ',
     . 'height'/'% push size then font onto stack to select font'/
     . '/choosefont {findfont exch scalefont setfont} def'/
     . '% paint object on top of stack at current position')
  231 format(
     . '/debug {dup /str 30 string def str cvs show} def'/
     . '% finding width of longest ordinate label (current width on ',
     . 'stack) ...'/'/max 0 def /mxw{/w exch def w max gt {/max w ',
     . 'def} if}def'/'% string centring, justification in the four ',
     . 'quadrants'/'% and absicca and ordinate label centring ...'/
     . '/showru{show}def'/'/showcu{dup stringwidth pop 2 div neg ',
     . '0 rmoveto show}def'/'/showlu{dup stringwidth pop neg 0 ',
     . 'rmoveto show}def'/'/showrc{0 h 2 div neg rmoveto show}def'/
     . '/showcc{dup stringwidth pop 2 div neg h 2 div neg rmoveto ',
     . 'show}def'/'/showlc{dup stringwidth pop neg h 2 div neg ',
     . 'rmoveto show}def'/'/showrd{0 h neg rmoveto show}def')
  240 format(
     . '/showld{dup stringwidth pop neg h neg rmoveto show}def'/
     . '/showcd{dup stringwidth pop 2 div neg h neg rmoveto show}def'/
     . '/showa{gsave fontsize /Symbol choosefont'/
     . '       dup stringwidth pop 2 div neg h 1.5 mul neg rmoveto ',
     . 'show grestore}def'/
     . '/showA{gsave fontsize /Symbol choosefont'/
     . '       dup stringwidth pop 2 div neg h 0.5 mul rmoveto ',
     . 'show grestore}def'/
     .  '/showo{gsave fontsize /Symbol choosefont'/
     . '       dup stringwidth pop dup mxw h 2 div add neg h 2 div ',
     . 'neg rmoveto show grestore}def'/
     .  '/showO{gsave fontsize /Symbol choosefont'/
     . '       dup h 2 div h 2 div neg rmoveto show grestore}def'/
     .  '% aliases for pslabl ...'/
     . '/show0{showcc}def /show1{show}def /show2{showrd}def ',
     . '/show3{showld}def'/'/show4{showlu}def')
  241 format('% newline macros for ',
     . 'horizontal and vertical labelling ...'/
     . '/newlinehf{currentpoint fontsize sub moveto}def'/
     . '/newlinevf{currentpoint exch fontsize add exch moveto}def'/
     . '/newlinehr{currentpoint fontsize add moveto}def'/
     . '/newlinevr{currentpoint exch fontsize sub exch moveto}def'/
     . '/savecurrentpoint {currentpoint /ybak exch def /xbak exch ',
     . 'def} def')
  242 format('% For other special symbols:'/
     . '/Times-Roman findfont'/
     . 'dup length dict begin'/
     . '   {1 index /FID ne {def} {pop pop} ifelse} forall'/
     . '   /Encoding ISOLatin1Encoding def'/
     . '   currentdict'/
     . 'end'/
     . '/Times-Roman-ISOLatin1 exch definefont pop')

  250 format('% end preamble ...'/)

      end


      block data plt
      implicit none
C --- Common variables for plot parameters ---
C plml,r,b,t: mediums's units for left,right,bottom and top coordinates
C             demarcating edge of plot medium.
C plul,r,b,t: user's units for left,right,bottom and top coordinates
C             demarcating edge of plot medium.
C plgl,r,b,t: user's units marking frame of interior of a current plot.
C xcm,ycm:    current pen position in medium's units.
C xcu,ycu:    current pen position in user's units.
C pcol:       current pen color
C stream:     stream where to send plot output (hardware dependent)
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

      data cosa /0d0/, sina /1d0/
      data plml /0d0/, plmr /1d0/, plmb /0d0/, plmt /1d0/
      data plgl /0d0/, plgr /1d0/, plgb /0d0/, plgt /1d0/
      data plpl /0d0/, plpr /1d0/, plpb /0d0/, plpt /1d0/
      data plul /0d0/, plur /1d0/, plub /0d0/, plut /1d0/
      data xcm  /0d0/, ycm  /0d0/, strtch /1d0/ ndpi /600/ ndwrt /1/
      data lintyp /1/, linl3 /0d0/, linl4 /0d0/ ltbld /2/
      data stream /11/, pstat /0/, noclip /.true./
      data logx /.false./, logy /.false./
      end
