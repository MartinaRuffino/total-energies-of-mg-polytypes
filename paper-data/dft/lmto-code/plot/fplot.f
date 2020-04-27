C#define unix
C Main for plotting file data.
C Updates
Cu  24 Jun 15 (3.50) New -shftm=#,#, circle symbol extended to partial arc
Cu  24 Jun 15 (3.49) New wiggly line as symbole
Cu  25 Jan 15 (3.48) Lines in symbols use default color from line type, or can be specified by coll
Cu                   Also enable special characters through \{\nnn}}},
Cu                   e.g. &\{k}_\{~\{\{\136}}} makes italic k, followed by subscript perp symbol
Cu  15 Mar 12 (3.45) Added row index as symbol
Cu  25 May 11 (3.44) Options -frme:xor:yab
Cu  01 Mar 11 (3.43) Four color weights for color weight plots
Cu  11 Mar 10 (3.42) Small changes to circumvent bugs in ps -> pdf conversion
Cu   9 Dec 09 (3.41) Added options to contour plots; contour plot bug fix
Cu  23 Oct 09 (3.40) Added italics: font -i# or & key, analog to ~ key in strings
Cu  22 Apr 09 (3.39) pldos -ref handles 2 color weights
Cu  08 Jun 08 (3.37) Continuously varying color lines, 3 colors
Cu   8 Jul 07 (3.36) Contour plot maker prints out area under contour
Cu  30 Jun 07 (3.35) Contour plots can replicate images
Cu  17 May 07 (3.34) Data Interpolation for logx,logy done on log scale;
Cu                   Try to stabilize points used for interpolation
Cu  25 Apr 07 (3.33) numbers written in fortran 'E' format are converted
Cu  02 Apr 07 (ATP)  Rewrote string write for postscript
Cu  23 Aug 06 (3.31) New LaTeX postprocessing of strings
Cu   6 Jun 06 (3.29) Continuously varying color lines, two colors
Cu   6 Jun 06 (3.28) Continuously varying color when drawing lines
Cu  19 Aug 04 (3.27) fixed -rot switch; also text after '#'
Cu                   is ignored for # anywhere in line
Cu  01 Jul 04        symbol can have its old bold
Cu  20 May 04        Extra file read options -r:
Cu  13 May 04        Extra option in errbar
Cu   4 Jan 04        can specify style in key
Cu  13 Jan 03        fplot loads xmin,xmax,ymin,ymax if available
Cu  12 Jun 02        data mapping now passed as string to mapdat
Cu  18 Apr 02        Added nr,nc to variables table when
Cu  26 Mar 02 (3.20) Added 'configure' script
Cu  22 Feb 01 (3.18) Some small bug fixes
Cu  31 Jan 01        revised -frme option; added :theta
Cu  19 Dec 00        -disp option
Cu  12 Jul 00        added clip to lt switches
Cu   3 May 00 (3.17) substituted calls to cmdstr for calls to nxtarg
Cu                   in preparation for macro generation.
Cu                   Requires SLATSM.36c.
Cu   6 Aug 98 (3.16) added -abf.  Requires SLATSM.35
Cu  31 Jul 97 (3.11) added new options to pstr.
C     Note for contour plots:
C     ROWS in data file  correspond to the y axis
C     COLUMNS in data correspond to the x axis
C#ifdef unix
      subroutine fmain
C#endif
      implicit none
      character*120 first,second,dspopt,rdcsw
      character*1 f2(40),dc
C For argument parsing and misc.
C iarg0 : index to argument starting current frame
C iargsh: index to argument starting new frame
      integer imx, NULLI
      parameter (imx=750, NULLI=-99999)
      integer garg,partok,fopng,fopnx,mkdlst,rdm,parg
      integer it(500),isp(imx),wsp(imx),iarg0,iargsh,n,i,j,
     .  k,j1,j2,ifi,i1mach,ipr,ip,ixv(imx),lxbrk,
     .  lclip,nxargs,rdops
      logical nxtarg,lsequ,a2bin,ltmp,setuu,ldisp,noclos
      double precision xxv(imx),xxpl(5),wk(imx),vsn,plu(6)
      logical flgnx,flgny,cmdopt,logx,logy,nlogx,nlogy
      equivalence (f2,first)
      double precision bound(6),xxmin,xxmax,yymin,yymax,zzmin,zzmax,
     .  xmin,xmax,ymin,ymax,zmin,zmax,pad,d1mach,big
      equivalence (bound(1),xxmin),(bound(2),xxmax),
     .            (bound(3),yymin),(bound(4),yymax),
     .            (bound(5),zzmin),(bound(6),zzmax)
C for work and data arrays
      integer nr,np,nc,ndpi,nr0,nc0
      real(8), allocatable, dimension(:) :: dat,xp,yp,zp,wp,wp2,wp3,fxp,rwk,wk2
C for frame
      double precision ts(50,3),rmt(3),rmnt(3),tcp(3),fcs(0:3),ffcs(0:4),theta,yab,xor
      integer mt(3),ntx,nty,ntz,fontsz,ltfrme(3),modxy(2),ffont(2)
      character font*19
      integer xnum,ynum,znum
C for 3D: parm3d(1..3) are rotations about z, x', z'
C         frm3d 1-4: projections of 3D coordinates onto plane
C                 5: total number of symbols
      integer nsymc !osym,osyma,oiwk,otabpr,
      real(8), allocatable, dimension(:) :: sym, rsyma
      integer, allocatable, dimension(:) :: iwk, tabpr
      double precision parm3d(6),frm3d(10),rotm(3,3),bsrad,pi
      parameter(pi=3.14159265358979324d0)
      logical flg3d
      character*40 bslist(0:imx)
C for a general label
      double precision rots
      character string*120, units*1
c for restricting portion of the medium
      double precision clip(4),clin(4),uuytox
C for key (xkey: x,y; length; spacing; if-blank)
      logical lkey
      double precision xkey(6)
      character legend*120
C for plot line type, symbols and attributes
C crvcs(0) is a flag indicating how line is to be colored and filled
C          For postscript (see plcrve):
C           1s digit
C         <=0  lines terminated with 'stroke'
C           1  lines terminated with 'closepath'
C           2  region filled with grey scale only
C           3  region filled with color
C           10s digit
C           0  lines drawn black
C           1  lines drawn colored
C           100s digit
C           1  closes curve without 'stroke' (for fill without border)
C Variables ending in 'cs' keeps information about line type and color.
C   var      contents
C   fcs      default color for frame
C   ffcs     frame fill color specification
C   symcs    symbol specification
C   crvcs    line specification
C __cs(0) is used by plcrve and specifies how a path is to be ended,
C         whether it is in color or B&W
C __cs(1..3) specify the RGB color
C __cs(4) specifies the symbol bold (used in cymcs only for now)
C         cs(4)<0 => take value from -lt specification
C ffcs(4) (frame fill) 1 => fill  0 => do not fill
C         cs(4)<0 => take value from -lt specification
      integer symtyp,lintyp,ltbld,lntypD,ltbldD
      double precision syma(10),symcs(0:7),crvcs(0:7),crvcsw(0:3,4),
     .  ltypn(10),ltypa,ltypb,ltypc,ltypd
      equivalence (ltypa,ltypn(1)),(ltypb,ltypn(2)),
     .            (ltypc,ltypn(3)),(ltypd,ltypn(4))
C Plot labels
      character*20 xnfmt,ynfmt,znfmt
      character*120 xlabel, ylabel, zlabel, title
C Mapping of abscissa, ordinate, and third column (usage depends on mode)
      integer ix(6)
C Special to plot
      integer nfile,nfrme
      logical pass1,flagf,lsort,fstplt,errbar,afcmd
      double precision xx,werrby(2)
      character prmed*72, ppstr*80, outs*80, dispfn*80
C For column mapping
      integer ssvdef,ncolsy
      parameter (ssvdef=256)
      character*(ssvdef) sxdef,sfxdef,sydef,scolsy,scolsw(4),
     .  sincl,sexcl,filel,strn,strn2,pntlst
C For latex postprocessing
      integer ntxstr,lbbox
      character*(ssvdef) texstr(2,100)
      character bkslsh*1
C For contour plotting
      logical flagcp
      integer ncval,ncval0,dupctr(2,2),ificon,cclose,ncvcol
      double precision cval(imx),lntypc(4,imx),contcol(3,imx),frmefontd
      character term*1,ct*1,frmefont*4
      external drwc
      integer, parameter :: nfont=5,strnsiz=13
      character(len=nfont), parameter :: fontlist='thibs'
C      character(len=nfont*strnsiz), parameter ::
C     .  fontnamtab =
C     . '/Times-Roman /Helvetica   /Times-Italic/Times-Bold  /Symbol      '
! C Heap allocation
!       integer wksize
!       parameter(wksize=140 000 000)
!       integer w(wksize)
!       common /w/ w
      procedure(integer) a2vec,iprint,nargf,wordsw

      vsn = 3.50d0
!       call pshpr(0)
!       call wkinit(wksize)
!       i = 0
!       call wkprnt(i)
!       call poppr
      call bin2a0(10)
      i = 0
      ltmp = a2bin('1 ',ix,2,0,' ',i,1)
C      dspopt = ' '
C     dispfn = 'ghostview -geometry +000+000'
C     dispfn = 'gs'
      dispfn = 'gv -geometry +000+000'
C ... Parse preliminary arguments; get index to starting argument
      iarg0 = 0
      if (nargf() == 1) call err(2,vsn,' ')
      if (cmdopt('-h',2,0,prmed) .or. cmdopt('--h',3,0,prmed) .or. cmdopt('--help',6,0,prmed)) then
        call info2(0,0,0,' fplot version %,2;2d',vsn,0)
        call err(3,vsn,' ')
        call cexit(0,1)
      endif
      if (cmdopt('--version',9,0,prmed) .or.
     .    cmdopt('-V ',3,0,prmed) .or.
     .    cmdopt('-v ',3,0,prmed)) then
        call info2(0,0,0,' fplot version %,2;2d',vsn,0)
        call cexit(0,1)
      endif

      prmed = ' ps'
      if (cmdopt('-plm=',4,0,prmed)) iarg0 = iarg0+1
      call plstmv(xx,xx,xx,xx,prmed,.false.)
      if (cmdopt('-shftm=',7,0,prmed)) iarg0 = iarg0+1
      if (cmdopt('-ascii',6,0,prmed)) iarg0 = iarg0+1
C ... Default verbosities for first, second pass
      call pshpr(40)
      call pshpr(30)
      if (cmdopt('-pr',3,0,first)) then
        n = 3
        if (a2bin(first,ipr,2,0,' ',n,-1)) then
          call pshpr(ipr)
          call pshpr(ipr)
        endif
      endif
      iargsh = iarg0
      call pltini(prmed(2:10),0,cmdopt('-rot ',5,0,strn))
      if (cmdopt('-db',3,0,first)) call setpr(121)
      call getpr(ipr)
      clip(1) = 0d0
      clip(2) = 1d0
      clip(3) = 0d0
      clip(4) = 1d0
      clin(1) = 0d0
      clin(2) = 1d0
      clin(3) = 0d0
      clin(4) = 1d0

C ... true until a frame defined or a set of points has been input
c     nofdef = .true.
C ... Global defaults
      pad = .1d0
      ndpi = 600
      font = '/Times-Roman'
      fontsz = 24
      call setfnt(font,fontsz)
      ltfrme(1) = 3
      ltfrme(2) = 0
      ltfrme(3) = 0
      big = sqrt(d1mach(2))/100
      lntypD = 1
      ltbldD = 3
      call pshpr(0)
      call pltsts(.false.,.false.,0)
      call poppr
      nlogx= .false.
      nlogy= .false.
      logx = .false.
      logy = .false.
      call plndpi(ndpi)
      ldisp = .false.
C     Angle between (x,y) axes
      theta = pi/2
      noclos = .false.
C     number of strings to convert to TeX
      ntxstr = 0
C     If to recalculate BoundingBox
      lbbox = 0
      nfrme = 0
      frmefont = ' '
      call frmeprms(3,ffcs,xor,yab,frmefont) ! Initialize frame parameters

C --- Start of new frame ---
    2 continue
      call dpzero(xkey,6)
      xkey(4) = -1
      uuytox = 0
      fcs(0:3) = [3d0,0d0,0d0,0d0]
      ncvcol = 0 ! Contours use default color

C      call dpzero(ts,150)
C      ts(1,1:3) = [-1d0,-1d0,-1d0]
C      tcp(1) = -1d30
C      tcp(2) = tcp(1)
C      tcp(3) = tcp(1)
C      rmt(1:3) = [.025d0,.025d0,.025d0]
C      rmnt(1:3)= [.6d0,.6d0,.6d0]
C      ltfrme(1) = 3
C      ltfrme(2) = 0
C      ltfrme(3) = 0
C      mt(1) = -1
C      mt(2) = -1
C      mt(3) = -1
C      ntx = -1
C      nty = -1
C      ntz = -1
C      modxy(1) = -1
C      modxy(2) = -1
C      xnfmt = '%;4d'
C      ynfmt = '%;4d'
C      znfmt = '%;4d'
C      xnum = 1
C      ynum = 1
C      znum = 1
C      xlabel = ' '
C      ylabel = ' '
C      zlabel = ' '
C      title  = ' '

      call frmelabel(mt,ntx,nty,ntz,ffont,ltfrme,modxy,xnum,ynum,znum,
     .  ts,rmt,rmnt,tcp,xnfmt,ynfmt,znfmt,xlabel,ylabel,zlabel,title)

      call dpzero(rotm,9)
      rotm(1,1) = 1
      rotm(2,2) = 1
      rotm(3,3) = 1

      lkey = .false.
      pass1 = .true.
      xmax = -big
      ymax = xmax
      zmax = xmax
      xmin =  big
      ymin = xmin
      zmin = xmin
      xxmin = big
      xxmax = big
      yymin = big
      yymax = big
      zzmin = big
      zzmax = big
      plu(2) = xmax
      plu(4) = ymax
      call dpzero(parm3d,6)
      call dpzero(frm3d,10)
      flg3d  = .false.
      frm3d(2) = -big
      frm3d(4) = frm3d(2)
      frm3d(1) =  big
      frm3d(3) = frm3d(1)
      flgnx = .false.
      flgny = .false.
      ificon = 0

C --- Start of second pass ---
    3 continue
      lintyp = lntypD
      ltbld  = ltbldD
      ltypa  = 2     ! size of dash in dashed line
      ltypb  = .5d0  ! size of space in dashed line
      ltypc  = 0     ! size of second dash
      ltypd  = ltypb ! size of second space
      nfile = 0
      nsymc = 0
      flagf = .true.
C     Set or reset current argument to iarg0+1
      k = nxargs(iarg0+1)
C     Argument doesn't exist ; we are done
      if (.not. nxtarg(first)) goto 199
      k = nxargs(nxargs(-1)-1)

C --- Start of new file ---
    5 continue

C     10 October 2016 revert to default line type after each data set
C      lintyp = lntypD
C      ltbld  = ltbldD
C      ltypa  = 2
C      ltypb  = .5d0
C      ltypc  = 0
C      ltypd  = ltypb

      rdcsw = ' '
      rdops = 0
      ncval = 0
      ncval0= 0
      ix(1) = 1
      ix(2) = 2
      ix(3) = 3
      ix(4) = 4
      ix(5) = 0
      ix(6) = 0
      symtyp = 0
      call dpzero(syma,10)
      syma(1) = 1
      symcs(0) = 3
      symcs(1) = 1
      symcs(2) = 1
      symcs(3) = 1
      symcs(4) = -1
      call dpzero(crvcs,8)
      crvcs(0) = -2
      call dpzero(crvcsw,4*4)
      call iinit(dupctr,4)
      sxdef = ' '
      sfxdef = ' '
      sydef = ' '
      sincl = ' '
      sexcl = ' '
      pntlst = ' '
      scolsy = ' '
      scolsw(1) = ' '
      scolsw(2) = ' '
      scolsw(3) = ' '
      scolsw(4) = ' '
      ncolsy = 0
      isp(1) = 0
      lsort  = .false.
      errbar = .false.
      nr = 0
      nc = 0
      flagcp = .false.
      legend = ' '
      lxbrk = 0
      lclip = 0
      xxpl(3) = 0
      xxpl(4) = 0
      xxpl(5) = 0

C --- Read next argument ---
   15 continue
      if (.not. nxtarg(first)) goto 99

      if (ipr > 50)
     .  call awrit0(' fplot: parsing '//first,outs,-80,i1mach(2))

C --- Case -frme: new frame ---
C     This switch serves a dual role:
C     1 as a terminator to a prior frame, if it exists
C     2 Defining parmeters for the next frame.
C     * Note also that the frame cannot be drawn until
C       subsequent switches such as -frmt or -x, or data files are parsed.
C       They can all affect how the frame is drawn.
C     Point(*) together with (1) and (2) complicates the algorithm,
C     because parameters for two frames must be kept at the same time.
C     This is accomplished with a "frame parameters" stack.
C     Routine frmeprms can (1) set default values, (2) push or pop current
C     values on the stack and (3) exchange current values with the stack.
C     nfrme keeps track of what is on the stack:
C       nfrme      this frme            next frme
C         0        current values         default
C         1        current values      current values
C        -1           stack            current values
      if (first(1:5) == '-frme') then
        call pshpr(0)
        call pltsts(logx,logy,0)
        call poppr
        yab = NULLI
        xor = NULLI
        frmefont = ' '
        ip = 0
        nlogx = .false.
        nlogy = .false.
        dc = first(6:6)
        if (dc /= ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 6-1
   10   continue
        j2 = j2+1
        if (first(j2:j2) == dc) goto 10
        j1 = min(len(first),j2)
        call nwordg(first,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (first(j1:j1+5) == 'theta=')  then
            i = j1+5
            j = a2vec(first,j2,i,4,dc//' ',2,2,1,ixv,theta)
            if (j /= 1) goto 21
          elseif (first(j1:j1+3) == 'col=')  then
            i = j1+3
            j = a2vec(first,j2,i,4,dc//', ',3,3,3,ixv,ffcs(1))
            if (j /= 0) then
              ffcs(0) = 103
              if (j < 0) goto 22
              if (j == 1) call dvset(ffcs(1),2,3,ffcs(1))
            endif
          elseif (first(j1:j1+3) == 'yab=')  then
            i = j1+3
            j = a2vec(first,j2,i,4,dc//' ',2,2,1,ixv,yab)
            if (j /= 1) goto 21
          elseif (first(j1:j1+4) == 'font=')  then
            i = j1
            call chrps2(first,dc//' ',2,len(first),i,j)
            frmefont = first(j1+5:i)
          elseif (first(j1:j1+3) == 'xor=')  then
            i = j1+3
            j = a2vec(first,j2,i,4,dc//' ',2,2,1,ixv,xor)
            if (j /= 1) goto 21
          elseif (first(j1:j2) == 'lx')  then
            nlogx = .true.
          elseif (first(j1:j2) == 'nofill')  then
            ffcs(4) = 0
          elseif (first(j1:j2) == 'fill')  then
            ffcs(4) = 1
          elseif (first(j1:j2) == 'ly')  then
            nlogy = .true.
          elseif (first(j1:j2) == 'lxy')  then
            nlogx = .true.
            nlogy = .true.
          else
            goto 21
          endif
          goto 10
        endif
        endif
        if (garg(' ',0,4,', ',2,4,i,it,clin) /= 4) goto 22
C       Set current cmd argument index to one after '-frme'
        iargsh = nxargs(-1)-1
!!      if (nfrme == 0) goto 99  ! nothing to do
        if (pass1) then     ! Retain new parms on 2nd pass only
          call s8tor8(frmefont,frmefontd)
          call frmeprms(6,ffcs,xor,yab,frmefont)
          nfrme = 1
        else
C         If new frme parameters are ahead of current frame, exchange with stack
C         if (nfrme == 1) then
            call frmeprms(4,ffcs,xor,yab,frmefont)
            nfrme = -1
C         endif
        endif
        goto 99
      endif

C --- Switches for this frame ---
C ... Avoid lots of tests for switches if this arg not a switch
      ip = 0
      if (.not. lsequ(first,'-',1,' ',n)) goto 30
      if (first(1:5) == '-rot ') then
      elseif (first(1:11) == '-plaintext ') then
      elseif (first(1:5) == '-disp') then
        ldisp = .true.
        outs = ' '
        if (first(1:7) == '-disp:p') then
          outs = ' -portrait'
        elseif (first(1:7) == '-disp:l') then
          outs = ' -landscape'
        elseif (first(1:7) == '-disp:u') then
          outs = ' -upsidedown'
        elseif (first(1:7) == '-disp:s') then
          outs = ' -seascape'
        elseif (first(1:11) == '-disp:exec=') then
          dispfn = first(12:)
        endif
        if (pass1) then
C          first = dspopt
C          call strip(first,i,j)
C          dspopt = first(i:max(j,1)) // outs
          first = dispfn
          call strip(first,i,j)
          dispfn = first(i:max(j,1)) // outs
        endif
      elseif (first(1:3) == '-pr') then
      else if (lsequ(first,'-f ',3,' ',n)) then
        if (.not. nxtarg(first)) goto 22
        if (pass1) then
          i = fopng(first,-1,0)
          rewind i
          if (.not. afcmd('fplot',i)) goto 22
          call fclose(i)
        endif
        if (iprint() >= 100) call acmdop(first,1,1)
      else if (first(1:3) == '-x ') then
        if (garg(' ',0,4,', ',2,2,i,it,xxmin) /= 2) goto 22
      else if (first(1:3) == '-y ') then
        if (garg(' ',0,4,', ',2,2,i,it,yymin) /= 2) goto 22
      else if (first(1:3) == '-z ') then
        if (garg(' ',0,4,', ',2,2,i,it,zzmin) /= 2) goto 22
      elseif (garg('-tm',-1,1,'xy',2,0,i,it,xxv) == 1) then
        k = it(1)
        if (garg(' ',0,4,':,;~@ ',6,50,i,it,xxv) < 1) goto 22
C ...   Spacing
        ts(1,k) = xxv(1)
        i = i-1
        j = 1
  55    continue
C ...   Number of tics per major tic
        if (i > 0 .and. it(j) == 1) then
          j = j+1
          i = i-1
          mt(k) = nint(xxv(j))
          goto 55
        endif
C ...   Position of major tic
        if (i > 0 .and. it(j) == 2) then
          j = j+1
          i = i-1
          tcp(k) = xxv(j)
          goto 55
        endif
C ...   Size of major tic
        if (i > 0 .and. it(j) == 3) then
          j = j+1
          i = i-1
          rmt(k) = xxv(j)
          goto 55
        endif
C ...   Size of major tic
        if (i > 0 .and. it(j) == 4) then
          j = j+1
          i = i-1
          rmnt(k) = xxv(j)
          goto 55
        endif
C ...   Mode of tic (log only)
        if (i > 0 .and. it(j) == 5) then
          j = j+1
          i = i-1
C          if (k == 1 .and. .not. logx .or.
C     .        k == 2 .and. .not. logy .or. k > 2) goto 22
          modxy(k) = nint(xxv(j))
          if  (modxy(k) == 5) then
            call dcopy(i,xxv(j+1),1,ts(1,k),1)
            modxy(k) = 10+i
          endif
        endif
      else if (first(1:5) == '-frmt') then
        if (.not. nxtarg(second)) goto 22
        k = 0
        k = parg('col=',4,second,k,len(second),', ',2,3,it,fcs(1))
        if (k < 0) goto 22
        fcs(0) = 10
        k = 0
        if (parg('th=',2,second,k,len(second),', ',2,3,it,ltfrme) < 0) goto 22
        if (ltfrme(1) == -1) flagf = .false.
      elseif (garg('-p',-1,4,' ',1,1,i,it,pad) /= 0) then
        if (i < 0) goto 21
      elseif (garg('-ndpi=',-1,2,' ',1,1,i,it,ndpi) /= 0) then
        if (i < 0) goto 21
        call plndpi(ndpi)
      elseif (garg('-lbl',-1,1,'umxy ',5,0,i,it,xxv) == 1) then
        units = 'u'
        if (it(1) == 2) units = 'm'
        xxv(3) = 1
C   ... Get vector of abscissas or ordinates to put labels
        if (it(1) == 3 .or. it(1) == 4) then
          if (.not. nxtarg(string)) goto 22
          j = mkdlst(string,-1024*d1mach(3),50,wk)
          if (garg(' ',0,4,': ',2,2,i,it(2),xxv(2)) < 1) goto 22
          if (it(1) == 4) xxv(1) = xxv(2)
        else
          if (garg(' ',0,4,',: ',3,3,i,it,xxv) < 2) goto 22
          j = 1
          wk(1) = xxv(1)
        endif
C   ... Pick up centering info
        if (.not. nxtarg(ppstr)) goto 22
        k = 0
        rots = 0
        if (parg('rot=',4,ppstr,k,len(ppstr),', ',2,3,it(2),rots) < 0) goto 22
        do  19  i = 1, j
          if (it(1) == 4) then
            xxv(2) = wk(i)
          else
            xxv(1) = wk(i)
          endif
          if (.not. nxtarg(string)) goto 22
          if (.not. pass1) call pstr(xxv(1),xxv(2),rots,string,
     .      0d0,0d0,0d0,ppstr,units,nint(xxv(3)))
          if (ppstr(1:2) == 'tx') then
            if (pass1) then
              ntxstr = ntxstr+1
              texstr(1,ntxstr) = string
              lbbox = max(lbbox,1)
            endif
            if (.not. nxtarg(string)) goto 22
            if (pass1) then
              texstr(2,ntxstr) = string
            endif
          endif
   19   continue
      else if (first(1:3) == '-xl') then
        if (.not. nxtarg(xlabel)) goto 22
      else if (first(1:3) == '-yl') then
        if (.not. nxtarg(ylabel)) goto 22
      else if (first(1:3) == '-tl') then
        if (.not. nxtarg(title)) goto 22
      else if (first(1:5) == '-font') then
        if (.not. nxtarg(second)) goto 22
        call parsefont(second,font,fontsz,.true.) ! return font, fontsz
        if (.not. pass1) call setfnt(font,fontsz)
      else if (first(1:2) == '-k') then
        lkey = .true.
        if (garg(' ',0,4,',:; ',4,5,k,it,xxv) < 2) goto 21
C        if (garg(' ',0,4,',: ',3,4,k,it,xxv) < 2) goto 21
        xkey(1) = xxv(1)
        xkey(2) = xxv(2)
C ...   Resolve the syntax [:len][,spacing][;style]
        i = 2
   17   i = i+1
        if (i <= k) then
          if (it(i-1) == 1) xkey(3) = xxv(i)
          if (it(i-1) == 2) xkey(4) = xxv(i)
          if (it(i-1) == 3) xkey(6) = xxv(i)
          goto 17
        endif
      else if (first(1:5)=='-noxn' .or. first(1:5)=='-xn:0') then
        xnum = 0
      else if (first(1:5) == '-xn:t') then
        xnum = -1
      else if (first(1:5) == '-xn:b') then
        xnum = 1
      else if (first(1:5)=='-noyn' .or. first(1:5)=='-yn:0') then
        ynum = 0
      else if (first(1:5) == '-yn:r') then
        ynum = -1
      else if (first(1:5) == '-yn:l') then
        ynum = 1
      elseif (first(1:7) == '-fmtnx:') then
        xnfmt = first(8:)
      elseif (first(1:7) == '-fmtny:') then
        ynfmt = first(8:)
      else if (first(1:7) == '-aspect ') then
        if (garg(' ',0,4,' ',1,1,i,it,uuytox) /= 1) goto 22
      elseif (garg('-3d',-1,2,' ',1,1,i,it,j) /= 0) then
        flg3d = .true.
        if (garg(' ',0,4,', ',2,3,i,it,parm3d) < 1) goto 21
C        call awrit1('%6:2;6d',parm3d,80,i1mach(2),parm3d)
      else if (first(1:5) == '-rotp ') then
        if (.not. flg3d) then
          print *, '-rotp implemented only for 3d plots ...'
          goto 199
        endif
        if (.not. nxtarg(strn)) goto 21
        call a2rotm(strn,.true.,ipr-20,rotm)
C --- File-specific switches ---
      else if (first(1:3) == '-r:') then
        rdcsw = first(3:)
      else if (first(1:3) == '-qr') then
        rdops = 10*(rdops/10) + 1
      else if (first(1:3) == '-br') then
        rdops = 10*(rdops/10) + 2
      else if (first(1:3) == '-tp') then
        k = garg(' ',0,4,',:~ ',4,50,i,it,xxv)
        if (k < 1 .and. k /= -2) goto 22 ! k=-2 -> read something like n~
        allocate(rwk(300000))
        call expand(i,it,xxv,rwk,nr,nc)
        allocate(dat(max(nr,1)*max(nc,2)))
        dat(1:max(nr,1)*max(nc,2)) = rwk(1:max(nr,1)*max(nc,2))
        deallocate(rwk)
!         call rlse(odat)
!         call redfrr(odat, nr*max(nc,2))
        goto 30
      else if (first(1:4) == '-ins') then
        if (.not. nxtarg(second)) goto 22
        if (.not. pass1) then
          call skpblb(second,len(second),i)
          j = 0
          if (first(5:5) == 'f') j = 1
          call plwstr(second(1:i+1),j)
        endif
      else if (first(1:5) == '-itrp') then
        if (garg(' ',0,4,', ',2,5,k,it,xxpl) < 3) goto 21
      else if (first(1:4) == '-inc') then
        if (.not. nxtarg(sincl)) goto 22
        pntlst = '1:nr'
      else if (first(1:4) == '-map') then
   32   continue
        if (garg('-',0,1,'ie ',3,0,i,it,xxv) == 1) then
          if (it(1) == 1) then
            if (.not. nxtarg(sincl)) goto 22
            goto 32
          endif
          if (it(1) == 2) then
            if (.not. nxtarg(sexcl)) goto 22
            goto 32
          endif
        endif
        k = nxargs(nxargs(-1)-1)
        if (.not. nxtarg(pntlst)) goto 22
      else if (first(1:5) == '-sort') then
        lsort =.true.
      else if (first(1:3) == '-lt') then
        if (.not. nxtarg(second)) goto 22
C   ... Pick up the line type: get terminating character
        call wordg(second,10,'0-9-',1,i,j)
        term = second(j+1:j+1)
        if (term == ' ') term = ','
        i = a2vec(second,len(second),ip,2,term//' ',2,2,1,it,lintyp)
        if (i < 1) goto 22
        call awrit1('%x fplot: set line type %i',strn,120,0,lintyp)
C   ... Return here to resume parsing for arguments
   33   continue
        if (it(i) == 1) then
          if (second(ip+1:ip+4) == 'col=') then
            i = parg('col=',4,second,ip,len(second),term//', ',3,3,it,crvcs(1))
            if (i < 0) goto 22
            crvcs(0) = 10
            if (i < 3) call dcopy(3-i,crvcs(1),1,crvcs(i+1),1)
            call awrit1('%a col=%3:1;4d',strn,120,0,crvcs(1))
          elseif (second(ip+1:ip+5) == 'colw=') then
            i = parg('colw=',4,second,ip,len(second),term//', ',3,3,it,
     .        crvcsw(1,1))
            if (i < 0) goto 22
            crvcsw(0,1) = 10
            if (i < 3) call dcopy(3-i,crvcsw(1,1),1,crvcsw(i+1,1),1)
            call awrit1('%a colw=%3:1;4d',strn,120,0,crvcsw(1,1))
          elseif (second(ip+1:ip+6) == 'colw2=') then
            i = parg('colw2=',4,second,ip,len(second),term//', ',3,3,it,
     .        crvcsw(1,2))
            if (i < 0) goto 22
            crvcsw(0,2) = 10
            if (i < 3) call dcopy(3-i,crvcsw(1,2),1,crvcsw(i+1,2),1)
            call awrit1('%a colw2=%3:1;4d',strn,120,0,crvcsw(1,2))
          elseif (second(ip+1:ip+6) == 'colw3=') then
            i = parg('colw3=',4,second,ip,len(second),term//', ',3,3,it,
     .        crvcsw(1,3))
            if (i < 0) goto 22
            crvcsw(0,3) = 10
            if (i < 3) call dcopy(3-i,crvcsw(1,3),1,crvcsw(i+1,3),1)
            call awrit1('%a colw3=%3:1;4d',strn,120,0,crvcsw(1,3))
          elseif (second(ip+1:ip+6) == 'colw4=') then
            i = parg('colw4=',4,second,ip,len(second),term//', ',3,3,it,
     .        crvcsw(1,4))
            if (i < 0) goto 22
            crvcsw(0,4) = 10
            if (i < 3) call dcopy(3-i,crvcsw(1,4),1,crvcsw(i+1,4),1)
            call awrit1('%a colw4=%3:1;4d',strn,120,0,crvcsw(1,4))
          elseif (second(ip+1:ip+5) == 'clip=') then
            i = parg('clip=',2,second,ip,len(second),term//' ',2,1,it,
     .        lclip)
            if (i < 0) goto 22
            call awrit1('%a clip %i',strn,120,0,lclip)
          elseif (second(ip+1:ip+4) == 'brk=') then
            i = parg('brk=',2,second,ip,len(second),term//' ',2,1,it,
     .        lxbrk)
            if (i < 0) goto 22
            call awrit1('%a brk %i',strn,120,0,lxbrk)
          elseif (second(ip+1:ip+5) == 'fill=') then
            i = parg('fill=',2,second,ip,len(second),term//' ',2,1,it,k)
            if (i < 0) goto 22
            call awrit1('%a fill %i',strn,120,0,k)
            if (k >= 3) crvcs(0) = 13
            if (k == 2) crvcs(0) = 3
            if (k == 1) crvcs(0) = 10
          elseif (second(ip+1:ip+5) == 'bold=') then
            i =parg('bold=',2,second,ip,len(second),term//' ',2,1,it,ltbld)
            if (i < 0) goto 22
            call awrit1('%a bold %i',strn,120,0,ltbld)
          else
            ltypb = .5d0; ltypc = 0; ltypd = ltypb
            i = a2vec(second,len(second),ip,4,term//', ',3,3,8,it,ltypn)
            if (i < 1) goto 22
            call awrit2('%a pat (%n:1;4d)',strn,120,0,i,ltypn)
            goto 34
          endif
          goto 33
        endif
   34   if (iprint() >= 30 .and. .not. pass1) then
          call awrit1('%a',strn,120,-i1mach(2),0)
        endif
      else if (first(1:3) == '-s ') then
        if (.not. nxtarg(second)) goto 22
        call wordg(second,10,'a-z+',1,i,j)
C   ... Get symbol type, if represented as a name
        symtyp = 0
        if (j > 0) then
          call plsyn(second(i:j),symtyp)
          if (symtyp < 1) goto 22
        endif
C   ... Otherwise, it must be a number
        if (symtyp <= 0) then
          call wordg(second,10,'0-9-',1,i,j)
          i = 0
          if (.not. a2bin(second(1:j),symtyp,2,0,' ',i,j-1)) goto 22
        endif
        if (symtyp < -1) goto 22
        if (symtyp == 6) then ! Defaults for syma(2,3)
          syma(2) = 0; syma(3) = 360
        endif
        term = second(j+1:j+1)
        symcs(5:7) = 0 ! Assign color for symbol lines
        if (crvcs(0) == 10) symcs(5:7) = crvcs(1:3)
C   ... Return here to resume parsing for arguments
        j = j+1
   51   continue
        i = j+1
        call nwordg(second,0,term//' ',1,i,j)
        if (j >= i) then
          if (second(i:i+4) == 'fill=') then
            j = j+1
            k = i-1
            k = parg('fill=',4,second,k,len(second),', '//term,2,1,it,symcs(0))
C           if (k < 0 .or. symcs(0) > 2 ) goto 22
C            symcs(0) = symcs(0)+1
            goto 51
          elseif (second(i:i+3) == 'clip') then
            symtyp = symtyp+1000
            j = j+1
            goto 51
          elseif (second(i:i+4) == 'bold=') then
            j = j+1
            k = i-1
            k = parg('bold=',4,second,k,len(second),', '//term,2,1,it,symcs(4))
            if (k < 0) goto 22
            if (symcs(4) < 0) symcs(4) = -2  ! suppress line for symbol
            goto 51
          elseif (second(i:i+4) == 'coll=') then
            j = j+1
            k = i-1
            k = parg('coll=',4,second,k,len(second),', '//term,2,3,it,
     .        symcs(5))
            if (k /= 0) then
              if (k < 0) goto 22
              if (k == 1) call dvset(symcs(5),2,3,symcs(5))
              goto 51
            endif
          elseif (second(i:i+3) == 'col=') then
            j = j+1
            k = i-1
            k = parg('col=',4,second,k,len(second),', '//term,2,3,it,
     .        symcs(1))
            if (k /= 0) then
              if (k < 0) goto 22
              if (k == 1) call dvset(symcs(1),2,3,symcs(1))
              goto 51
            endif
          endif
          k = i-1
          call dpzero(syma,10)
          if (symtyp == 6) then ! Defaults for syma(2,3)
            syma(2) = 0; syma(3) = 360
          endif
          i = parg(' ',4,second,k,len(second),', ',2,8,it,syma)
          if (i < 0) goto 22
        endif
      else if (first(1:4) == '-con') then
        ificon = 0
        cclose = 1
        j = 4
        dc = first(j+1:j+1)
C   ... Return here to resume parsing for arguments
        j = j+1
   47   continue
        i = j+1
        if (wordsw(first,dc,'dup=','',k) > 0) then
          k = k-1
          call rxx(a2vec(first,len_trim(first),k,2,', '//dc,3,2,2,it,dupctr)<1,
     .      'fplot failed to parse '//first)
        endif
        if (wordsw(first,dc,'noclos','',k) > 0) cclose = 0
        if (wordsw(first,dc,'fn=','',k) > 0) then
          call nwordg(first,0,dc//' ',1,k,i)
          ificon = fopng(first(k:i),-1,0)
          rewind ificon
        endif
        pad = 0
        flagcp = .true.
        if (.not. nxtarg(second)) goto 21
        j = mkdlst(second,-1024*d1mach(3),imx-ncval0,xxv)
        if (wordsw(first,dc,'col=','',k) > 0) then
          k = k-1
          i = a2vec(first,len_trim(first),k,4,', '//dc,3,2,3*j,it,contcol)
          if (mod(i,3) /= 0) call rx('-con~col must contain triplets of numbers')
          ncvcol = i/3
        endif
C        call snot

        do  i  = 1, j
          cval(ncval0+i) = xxv(i)
          lntypc(1,ncval0+i) = lintyp
          lntypc(2,ncval0+i) = ltbld
          lntypc(3,ncval0+i) = ltypa
          lntypc(4,ncval0+i) = ltypb
        enddo
        ncval = ncval0+j

C ... Printout of contours
        if (iprint() >= 40 .and. .not. pass1) then
          print 357
  357     format(' fplot:ic  contour     line type')
          do  i = ncval0+1, ncval
            strn = ' '
            write(strn(1:8),'(i8)') i
            ip = 8
            call bin2a(' ',2,0,cval(i),4,0,72,strn,ip)
            ip = 20
            do  j = 1, 4
              call bin2a(' ',2,0,lntypc(j,i),4,0,72,strn,ip)
            enddo
            print *, strn(1:ip)
          enddo
          ncval0 = ncval
        endif
      else if (first(1:3) == '-bs') then
        if (garg(' ',0,4,'; ',2,1,k,it,bsrad) < 1) goto 21
C       decrement argument index, to look at it again
        k = nxargs(nxargs(-1)-1)
        if (.not. nxtarg(bslist)) goto 21

      elseif (garg('-nc=',-1,2,' ',1,1,i,it,nc) /= 0) then
        if (i < 0) goto 21

      elseif (garg('-l',-1,1,'01 ',1,0,i,it,xxv) == 1) then
        xkey(5) = min(it(1)-1,1)
        if (.not. nxtarg(legend)) goto 199

      else if (first(1:3) == '-1p') then
        if (pass1) then
          call awrit0(' fplot: start second pass this frame',' ',80,i1mach(2))
          goto 99
        endif
      else if (first(1:3) == '-nx') then
        flgnx = .true.
      else if (first(1:3) == '-ny') then
        flgny = .true.
      else if (first(1:4) == '-ord') then
        if (.not. nxtarg(sydef)) goto 22
      else if (first(1:4) == '-abf') then
        if (.not. nxtarg(sfxdef)) goto 22
      else if (first(1:4) == '-ab ') then
        if (.not. nxtarg(sxdef)) goto 22
      else if (first(1:6) == '-colsy' .or. first(1:5) == '-coll') then
        if (.not. nxtarg(scolsy)) goto 22
      else if (first(1:7) == '-colsw4') then
        if (.not. nxtarg(scolsw(4))) goto 22
      else if (first(1:7) == '-colsw3') then
        if (.not. nxtarg(scolsw(3))) goto 22
      else if (first(1:7) == '-colsw2') then
        if (.not. nxtarg(scolsw(2))) goto 22
      else if (first(1:6) == '-colsw') then
        if (.not. nxtarg(scolsw)) goto 22
      else if (first(1:4) == '-col') then
        if (garg(' ',0,2,', ',2,3,i,it,ix) < 2) goto 22
      else if (first(1:3) == '-ey') then
        if (garg(' ',0,4,', ',2,3,i,it,xxv) < 1) goto 22
        ix(3) = xxv(1)
        werrby(1) = 1
        werrby(2) = 0
        if (i >= 2) werrby(1) = xxv(2)
        if (i >= 3) werrby(2) = xxv(3)
        errbar = .true.
      else if (first(1:2) == '-v') then
        j = 2
        call parsyv(first,len(first),999,0,j)
C        call shosyv(0,0,0,i1mach(2))
C        pause
      else if (first(1:2) == '-c') then ! char defs handled by rdfiln
      else
        goto 21
      endif
      goto 15
C --- End of switches ---

   30 continue
      nfile = nfile+1
      i = 2
      if (pass1) i = 1
      strn = ' '
      call awrit1(' fplot file "'//first,strn,80,0,1)
      call awrit1('%a"  pass %i',strn,80,0,i)
      if (errbar .and. .not. pass1)
     .  call awrit1('%a (error bars, col %i)',strn,80,0,ix(3))
      fstplt = nfile == 1

C --- Read data from disk, single column -> second col ---
C     Check read options from rdcsw
      if (rdcsw /= ' ') then
        ct = rdcsw(1:1)
        call partk0(0,len(rdcsw),1,-1,0,len(rdcsw),-1,31,.false.)
        j = partok(rdcsw,'qr',' '//ct,ltmp,' ',0,0,0,0)
        if (ltmp) rdops = 10*(rdops/10) + 1
        j = partok(rdcsw,'br',' '//ct,ltmp,' ',0,0,0,0)
        if (ltmp) then
          rdops = 10*(rdops/10) + 2
C          j = partok(rdcsw,'br,',', ,'//ct,it,' ',-2,2,0,0)
C          if (j == 2) then
C            nr(ns+1) = it(1)
C            nc(ns+1) = it(2)
C            rdops = 10*(rdops/10) + 3
C          endif
        endif
      endif

C ... Open the file
      if (nr == 0 .and. first /= '-tp') then
        if (first == '.') then
          ifi = i1mach(1)
        else
          i = 1
          if (mod(rdops,10) == 2 .or. mod(rdops,10) == 3) i = 5
          ifi = fopnx(first,72,-1,-1)
          call word(first,1,j1,j2)
          if (ifi == 0) then
            print *, 'Exit -1 fplot: attempt to open nonexistent file "'//first(j1:j2)//'"'
            call cexit(-1,1)
          endif
          ifi = fopng(first,-1,i)
        endif

!         call defask(i)
!         i = i/i1mach(18) - 1
!         call defdr(odat,i)
        if (allocated(dat)) deallocate(dat)
        i = 67108864 ! 512MB or 8B elements
        allocate (dat(i))
        filel = ' '
        if (rdcsw /= ' ') then
          ct = rdcsw(1:1)
          call partk0(0,len(rdcsw),1,-1,0,len(rdcsw),-1,31,.false.)
          j = partok(rdcsw,'nc=','= '//ct,k,' ',-1,2,0,0)
          if (j == 1) nc = k
          j = partok(rdcsw,'nr=','= '//ct,k,' ',-1,2,0,0)
          if (j == 1) nr = k
          j = partok(rdcsw,'s=','= '//ct,k,' ',-1,2,0,0)
          if (j == 1) then
            do  j = 1, k
              if (mod(rdops,10) >= 2) read(ifi)
              if (mod(rdops,10) < 2) read(ifi,*)
            enddo
          endif
          j = partok(rdcsw,'open',' '//ct,noclos,' ',0,0,0,0)
          ltmp = .false.
C          j = partok(rdcsw,'spc',' '//ct,ltmp,' ',0,0,0,0)
C          if (ltmp) icast(ns) = 100
          rdcsw = ' '
        endif
        allocate(rwk(size(dat)))
        i = rdm(ifi,1000+rdops,i,filel,rwk,nr,nc)
        strn2 = ' '
        if (filel /= ' ' .and. pass1)
     .    call awrit0(' fplot label '''//first//'%a'': "'//filel//
     .    '%a"',strn2,len(strn2),-i1mach(2))
        if (i < 0)
     .    call fexit(-1,9,'fplot failed to parse file '//first,0)
!         call rlse(odat)
!         call defdr(odat, nr*max(nc,2))
        if (allocated(dat))  deallocate(dat)
        allocate(dat(nr*max(nc,2)))
        dat(1:nr*max(nc,2)) = rwk(1:nr*max(nc,2))
        deallocate(rwk)
C       Contour plotting
        if (dupctr(1,1) /= 0 .or. dupctr(2,1) /= 0) then
          rewind ifi
          deallocate(dat)
          do  i = 1, 2
            dupctr(i,2) = 0
            if (dupctr(i,1) < 0) dupctr(i,2) = 1
          enddo
          dupctr(1,1) = max(iabs(dupctr(1,1)),1)
          dupctr(2,1) = max(iabs(dupctr(2,1)),1)
          allocate(dat(nr*dupctr(1,1)*nc*dupctr(2,1)))
          allocate(rwk(nr*dupctr(1,1)*nc*dupctr(2,1)))
          i = rdm(ifi,1000+rdops,nr*dupctr(1,1)*nc*dupctr(2,1),filel,rwk,nr,nc)
          if (i < 0)
     .      call fexit(-1,9,'fplot failed to parse file '//first,0)
          i = nr + (nr-dupctr(1,2))*(dupctr(1,1)-1)
          j = nc + (nc-dupctr(2,2))*(dupctr(2,1)-1)
          call dupdat(nr,nc,i,j,dupctr(1,1),dupctr(2,1),
     .      dupctr(1,2),dupctr(2,2),rwk,dat)
          deallocate(rwk)
          nr = i
          nc = j
        endif

        if (ifi /= i1mach(1) .and. .not. noclos) call fclr(first,ifi)
C       call fclose(ifi)
        noclos = .false.
      endif

C --- Allocate xp,yp,zp.  Note: for 2d plots, zp are for error bars ---
      allocate(xp(max(nr,5001))); call dpzero(xp,size(xp))
      allocate(yp(max(nr,5001))); call dpzero(yp,size(yp))
      allocate(zp(max(nr,5001))); call dpzero(zp,size(zp))
      allocate(wp(max(nr,5001))); call dpzero(wp,size(wp))
      if (scolsw(4) /= ' ') then
        allocate(wp3(max(nr,5001))); wp3 = 0
      else
        allocate(wp3(1))
      endif
      if (scolsw(3) /= ' ') then
        allocate(wp2(max(nr,5001))); wp2 = 0
      else
        allocate(wp2(1))
      endif
      allocate(fxp(max(nr,5001))); fxp = 0
      if (flagcp) then
        allocate(rwk(nr*nc*ncval))
      else
        allocate(rwk(max(nr,5001)))
      endif
      allocate(wk2(max(nr,5001)))

C --- Entry point for multiple column plots ---
   40 continue
      np = nr
      outs = strn

C     Apr 2002 Add nr and nc to symbolic variables table
      call numsyv(ip)
      xx = 0
      call getsyv('nr',xx,j)
      nr0 = xx
      call getsyv('nc',xx,j)
      nc0 = xx
      call lodsyv('nr',1,dble(nr),j)
      call lodsyv('nc',1,dble(nc),j)

C      if (pntlst /= ' ') then
C        ip = 0
C        nread = a2vec(pntlst,len(pntlst),ip,2,',: ',3,3,100,isp,wsp)
C        pntlst = ' '
C      endif
C       if (garg(' ',-1,2,',: ',3,50,i,isp,wsp) < 1) then
C          k = nxargs(nxargs(-1)-1)
C          goto 22
C        endif
C

      if (scolsy /= ' ') then
        call mkilst(scolsy,j,ixv)
        if (j > imx) call rxi('fplot: increase imx, need',j)
   41   ncolsy = ncolsy+1
        if (ncolsy > j) goto 42
        ix(2) = ixv(ncolsy)
        call awrit1('%a column %i',outs,80,0,ix(2))
        if (ix(2) > nc .and. nc > 1) then
          call awrit0('%a (missing)',outs,80,-i1mach(2))
          goto 41
        endif
        if (scolsw(1) /= ' ') then
          call mkilst(scolsw(1),k,ixv)
          if (k /= j) call rx('mismatch between colsy and colsw')
          ix(3) = ixv(ncolsy)
        endif
        if (scolsw(2) /= ' ') then
          call mkilst(scolsw(2),k,ixv)
          if (k /= j) call rx('mismatch between colsy and colsw2')
          ix(4) = ixv(ncolsy)
        endif
        if (scolsw(3) /= ' ') then
          call mkilst(scolsw(3),k,ixv)
          if (k /= j) call rx('mismatch between colsy and colsw3')
          ix(5) = ixv(ncolsy)
        endif
        if (scolsw(4) /= ' ') then
          call mkilst(scolsw(4),k,ixv)
          if (k /= j) call rx('mismatch between colsy and colsw4')
          ix(6) = ixv(ncolsy)
        endif
        fstplt = nfile == 1 .and. ncolsy == 1
      endif
      if (ipr >= 30) call awrit1('%a',outs,80,-i1mach(2),1)
      if (.not. flg3d) then
        if (yymin /= big) then
          call lodsyv('ymin',1,yymin,j)
        endif
        if (yymax /= big) then
          call lodsyv('ymax',1,yymax,j)
        endif
        if (xxmin /= big) then
          call lodsyv('xmin',1,xxmin,j)
        endif
        if (xxmax /= big) then
          call lodsyv('xmax',1,xxmax,j)
        endif
        call mapdat(dat,nr,np,nc,ix,ix(min(2,max(nc,2))),ix(min(3,nc)),
     .    ix(min(4,nc)),ix(5),ix(6),sxdef,sydef,sfxdef,sincl,sexcl,
     .    flgnx,flgny,pntlst,lsort,xxpl,logx,logy,rwk,wk2,xp,
     .    yp,zp,wp,wp2,wp3,fxp)
      else
        call map3d(dat,nr,np,nc,parm3d,sxdef,sydef,sincl,sexcl,
     .    isp,wsp,rwk,wk2,xp,yp,zp,
     .    xxmin,xxmax,yymin,yymax,zzmin,zzmax,
     .    xmin,xmax,ymin,ymax,zmin,zmax,frm3d,rotm)
        if (pass1 .and. symtyp /= 0) frm3d(5) = frm3d(5) + np
      endif

C     Apr 2002 restore variables table
      call lodsyv('nr',1,dble(nr0),j)
      call lodsyv('nc',1,dble(nc0),j)
      call clrsyv(ip)

C --- Draw frame and line ---
      if (pass1 .and. .not. flg3d) then
        if (.not. flagcp) then
          call plboun(xp,yp,np,xmin,xmax,ymin,ymax)
        endif
      elseif (.not. pass1 .and. .not. flg3d) then
        call plntyp(lintyp, ltbld, ltypa, ltypb, ltypc, ltypd)
C   ... For contour plots, set plot boundaries
        ltmp = .false.
        if (flagcp) then
          if (xmax == -big) then
            xmin = -1
            xmax =  1
            ltmp = .true.
          endif
          if (ymax == -big) then
            ymin = -1
            ymax =  1
            ltmp = .true.
          endif
        else
          if (plu(2) == -big) then
            ltmp = .true.
            plu(1) = xmin
            plu(2) = xmax
          endif
          if (plu(4) == -big) then
            ltmp = .true.
            plu(3) = ymin
            plu(4) = ymax
          endif
        endif
C   ... Case to set or seset user's units
        if (ltmp) then
          ltmp =
     .      setuu(xmin,xmax,ymin,ymax,big,logx,logy,pad,plu,uuytox)
          if (ltmp) then
            call pltstu(plu(1),plu(2),plu(3),plu(4),uuytox)
          else
            call rx('fplot: missing x and/or y units for frame')
          endif
        endif

        if (flagf .and. fstplt) then
C         If new frme parameters are ahead of current frame, exchange with stack
          if (nfrme == -1) call frmeprms(4,ffcs,xor,yab,frmefont)
          call frme(
     .    tcp(1),ts(1,1),mt(1),rmt(1),rmnt(1),modxy(1),xnum,xnfmt,xlabel,yab,
     .    tcp(2),ts(1,2),mt(2),rmt(2),rmnt(2),modxy(2),ynum,ynfmt,ylabel,xor,
     .    sfxdef,np,xp,fxp,title,ltfrme,frmefont,fcs,ffcs)
C         Restore possible next frame parameters; reset stack to defaults
          if (frmefont /= ' ') call setfnt(font,fontsz) ! Restore fonts
          call frmeprms(5,ffcs,xor,yab,frmefont); nfrme = 0
C          call frmelabel(mt,ntx,nty,ntz,fontsz,ltfrme,modxy,xnum,ynum,znum,
C     .      ts,rmt,rmnt,tcp,xnfmt,ynfmt,znfmt,xlabel,ylabel,zlabel,title)
        endif
        if (lkey .and. legend /= ' ') then
          if (xkey(3) == 0d0) xkey(3) = plu(6)/15
          if (xkey(4) < 0) xkey(4) = plu(5)/8
          call key(nint(xkey(6)),xkey(1),xkey(2),xkey(4),xmin,xmax,ymin,
     .      ymax,symtyp,syma,symcs,crvcs,legend,nint(xkey(5)))
          if (.not. logy) xkey(2) = xkey(2) - xkey(3)
          if (      logy) xkey(2) = xkey(2)/dexp(xkey(3))
        endif
        if (flagcp) then
          call drwcin(nr,plu(1),plu(2),nc,plu(3),plu(4),ificon,cclose)
          call gcontr(dat,nr,nr,nc,cval,ncval,d1mach(2),
     .      crvcs,ncvcol,contcol,rwk,drwc,lntypc)
        else
          call plcrv(xp,yp,zp,wp,wp2,wp3,np,
     .      xmin,xmax,ymin,ymax,lxbrk+10*lclip,crvcs,crvcsw,11)
          if (errbar) call plsym(xp,yp,zp,np,8,werrby,symcs)
          call plsym(xp,yp,zp,np,symtyp,syma,symcs)
        endif
      elseif (.not. pass1 .and. flg3d) then
        if (logx .or. logy)
     .    call rx('fplot:  3D plots not implemented with log scale')
        if (flagcp) call rx('fplot:  3D contour plots not implemented')
        call plntyp(lintyp, ltbld, ltypa, ltypb, ltypc, ltypd)
        if (fstplt)
     .    call pltstu(frm3d(1),frm3d(2),frm3d(3),frm3d(4),uuytox)
        if (flagf .and. fstplt)
     .  call frme3d(tcp(1),ts(1,1),mt(1),ntx,rmt(1),xnum,xnfmt,xlabel,
     .              tcp(2),ts(1,2),mt(2),nty,rmt(2),ynum,ynfmt,ylabel,
     .              tcp(3),ts(1,3),mt(3),ntz,rmt(3),znum,znfmt,zlabel,
     .              xmin,xmax,ymin,ymax,zmin,zmax,parm3d,rotm,
     .              title,ltfrme,ffcs)
C       call plcrv(w(oxp),w(oyp),np,frm3d(1),frm3d(2),frm3d(3),frm3d(4),
C    .    lxbrk+10*lclip,crvcs,11)
        call pls3da(xp,yp,zp,np,dat,np,4,symtyp,
     .    syma,symcs,bsrad,bslist,nint(frm3d(5)),sym,rsyma,nsymc)
      else
        if (symtyp >= 0) nsymc = nsymc+1
        if (symtyp < 0) nsymc = nsymc+np
      endif
C ... Look for next file, or next column if scolsy set
      if (scolsy /= ' ') goto 40
C ... Re-entry if scolsy exhausted
   42 deallocate(xp,yp,zp,wp,fxp,rwk,wk2,dat)
      if (allocated(wp3)) deallocate (wp3)
      if (allocated(wp2)) deallocate (wp2)
      goto 5

C --- End of first or second pass ---
   99 continue
      call togprt
      call getpr(ipr)
      if (flg3d) then
        if (pass1) allocate(sym(nint(frm3d(5)+1)*4),rsyma(nsymc*16))
        if (pass1) rsyma = 0
        if (.not. pass1 .and. nint(frm3d(5)) /= 0) then
          allocate(iwk(2*nint(frm3d(5)+1)), tabpr(3*nint(frm3d(5)+1)**2))
          call pls3dp(nint(frm3d(5)),sym,rsyma,bslist,tabpr,iwk,rotm,parm3d)
          deallocate(iwk,tabpr)
        endif
        if (.not. pass1) deallocate(sym,rsyma)
      endif

      pass1 = .not. pass1

C ... Close of first  pass ... get plot boundaries
      if (.not. pass1) then
        if (xxmin /= big) xmin = xxmin
        if (xxmax /= big) xmax = xxmax
        if (yymin /= big) ymin = yymin
        if (yymax /= big) ymax = yymax
        ltmp = setuu(xmin,xmax,ymin,ymax,big,logx,logy,pad,plu,uuytox)
        if (ltmp .and. .not. flg3d)
     .    call pltstu(plu(1),plu(2),plu(3),plu(4),uuytox)
      endif
      if (.not. pass1) goto 3

C ... We are finished unless -frme terminated last frame ---
      if (iargsh <= iarg0) goto 199
      print *, 'fplot: ------ starting new frame ------'
      clip(1) = clin(1)
      clip(2) = clin(2)
      clip(3) = clin(3)
      clip(4) = clin(4)
      logx = nlogx
      logy = nlogy
      call plthet(theta)
      call pltsts(logx,logy,0)
      call pltstp(clip(1),clip(2),clip(3),clip(4))
C     Reset current argument to first argument of new frame
      iarg0 = iargsh
      goto 2

C --- Normal exit ---
  199 continue
      call pltdmp(0)

C ... Below assumes a postscript file ...
      if (lbbox /= 0) then
        call spchar(1,bkslsh)
        strn = 'awk -vbbox="`echo quit|gs -sDEVICE=bbox -dNOPAUSE '//
     .    'fplot.ps 2>&1|grep BoundingBox:|head -1`" '//
     .    '''{print ; if (NR == 1) '//
     .    'printf "%s'//bkslsh//'n", bbox;}'' fplot.ps >ps.ps'
        call strip(strn,j1,j2)
        write(*,454)  strn(j1:j2)
  454   format(/' File fplot.ps created ...',
     .          ' remake Bounding Box (file ps.ps) with system call:':/
     .          1x,a)
        call fsystm(strn(j1:j2),j)
        if (j /= 0) call rxi('oops: shell returned value',j)
      endif
      if (ntxstr /= 0) then
        call psfrag(ntxstr,texstr)
        print 455
  455   format(/
     .    ' File psfrag.tex was created ... '/
     .    ' Do the following to make TeX substitutions: (file fplot.ps)'/
     .    ' latex psfrag.tex'/' dvips -E -f psfrag.dvi > fplot.ps'
     .    )
  555   format(a)
      endif
      if (ldisp) then

        if (ntxstr /= 0) then
        print 456
  456   format(/' ... invoking latex with system call')
        call fsystm('latex psfrag.tex',j)
        if (j /= 0) call rxi('oops: shell returned value',j)
        print 457
  457   format(/' ... invoking dvips with system call')
        call fsystm('dvips -E -f psfrag.dvi > fplot.ps',j)
        if (j /= 0) call rxi('oops: shell returned value',j)

        write(*,454) strn(j1:j2)
        call fsystm(strn(j1:j2),j)
        if (j /= 0) call rxi('oops: shell returned value',j)

        endif

C        call strip(dspopt,i,j)
C        outs = 'ghostview ' // dspopt(i:max(j,1)) // ' fplot.ps'
C        print *, outs
        call strip(dispfn,i,j)
        if (lbbox /= 0) then
          outs = dispfn(1:max(j,1)) // ' ps.ps'
        else
          outs = dispfn(1:max(j,1)) // ' fplot.ps'
        endif
        call strip(outs,j1,j2)
        print 458, outs(j1:j2)
  458   format(/' viewing file with system call:':/1x,a)
        call fsystm(outs,j)
      endif
      return

C --- Print out error ---
   22 if (.not. nxtarg(second)) goto 21
      call skpblb(first,len(first),ip)
      call err(1,vsn,'fplot: parse error at ' // first(1:ip+2)//second)
   21 call err(1,vsn,'fplot: parse error at ' // first)
      end
      subroutine psfrag(ntxstr,texstr)
C- For latex postprocessing
      implicit none
      integer ntxstr
      character texstr(2,100)*(*)
      integer ifi,fopng,i,i1,i2,j1,j2
      character bkslsh*1

C     Generic call for backslash character
      call spchar(1,bkslsh)

C ... Create head
      ifi = fopng('psfrag.tex',-1,0)
      rewind ifi
      write(ifi,333)
     .  bkslsh,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh
     .  ,bkslsh,bkslsh
     .  ,bkslsh,bkslsh
     .  ,bkslsh,bkslsh
     .  ,bkslsh,bkslsh
     .  ,bkslsh,bkslsh,bkslsh
     .  ,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh
     .  ,bkslsh,bkslsh,bkslsh,bkslsh,bkslsh
  333 format(
     . a1,'documentclass{slides}'/
     . a1,'usepackage{graphicx}'/
     . a1,'usepackage{psfrag}'/
     . a1,'usepackage{latexsym}'/
     . a1,'usepackage{amsmath}'/
     . a1,'pagestyle{empty}'/
     . a1,'begin{document}'//
     . a1,'fontfamily{ptm}',a1,'selectfont'/
     . '%',a1,'fontfamily{cmr}',a1,'selectfont'/
     . '%',a1,'fontfamily{phv}',a1,'selectfont'//
     . '%',a1,'fontseries{sb}',a1,'selectfont'//
     . '% Customize fonts ',
     . 'within text: e.g. {',a1,'Huge',a1,'textsl{',a1,'textsf{end}}}'/
     . '%',a1,'small'/'%',a1,'normalsize'/'%',a1,'large'/'%',a1,'Large'/
     .    ,a1,'LARGE            % A sensible global average'/
     . '%',a1,'huge'/'%',a1,'Huge'//
     . '%',a1,'textrm'/'%',a1,'textsf'/'%',a1,'textbf'/
     . '%',a1,'textsl'/'%',a1,'textit'/
     .  )
      write(ifi,433)
  433 format(/'% --- psfrag substitutions ---')
      do  i = 1, ntxstr
        call strip(texstr(1,i),i1,i2)
        call strip(texstr(2,i),j1,j2)
        write(ifi,334) bkslsh,texstr(1,i)(i1:i2),texstr(2,i)(j1:j2)
  334   format(a1,'psfrag{',a,'}{',a,'}')
      enddo
      write(ifi,335) bkslsh,bkslsh,bkslsh
  335 format(/
     .  a1,'includegraphics[width=1.0',a1,'textwidth,clip]{ps.ps}'//
     .  a1,'end{document}'
     .  )

      call fclose(ifi)

      end
