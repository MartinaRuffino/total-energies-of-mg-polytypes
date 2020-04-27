      subroutine fmain
C- Band plotter using bnds file in standard format for lmto programs.
C ----------------------------------------------------------------------
Cl Local variables
Cl   ifi bands file
Cl   ifm file logical unit for merge file
Cl   ifw file logical unit for new merge file to be created
Cl   spintexture
Cl       >0 => treat color weights as spin texture
Cl             bands file must have 4 color weights, corresponding to
Cl             charge, x-, y-, and z- spin components
Cl       1s digit
Cl       1 => renormalize (x,y,z) -> (x,y,z)/d where d = |(x,y,z)|
Cl       2 => renormalize (x,y,z) -> (x,y,z)/q where q = charge component
Cl       3 => do not renormalize (x,y,z)
Cl       10s digit 1
Cl       0 => col1 => |x|, col2 => |y|, col3 => |z|
Cl       1 => Two color weights for x : col1 for x>0 col2 for x<0
Cl       2 => Two color weights for y : col1 for y>0 col2 for y<0
Cl       4 => Two color weights for z : col1 for z>0 col2 for z<0
Cr Remarks
Cr   Spin texture
Cr   spectral mode examples
Cr     plbnds -sp~out=eps~atop=10 spq.fe
Cr     plbnds -sp~out=png~atop=10 spq.fe
Cr     plbnds -sp~fmt=2~escl=1/13.6 spq.fe
Cr     plbnds -sp~atop=10~window=-4,4 spq.fe
Cr     plbnds -sp~fs~abot=1~atop=10 spq.fese
Cr     plbnds -sp~atop=10~window=-4,4~drawqp:bold=4:col=0,2,.3 spq.fe
Cr     plbnds -sp~atop=17~window=-1.2,1~drawqp:ef=-.41 spq-qsgwdmft.sr2ruo4
Cu Updates
Cu   11 Jan 18 New -sp~drawqp
Cu   15 Nov 17 New -sp~writeqp
Cu   08 Jun 17 New -fplot~ib1 and -fplot~ib2
Cu   07 Jun 17 Spin texture
Cu   09 Apr 17 tidied up spectral mode
Cu   20 Nov 16 New spectral mode -sp
Cu   24 Aug 15 New -skip= and -wscal options
Cu   16 Mar 11 Fermi level in header may be algebraic expression
Cu   08 Jul 08 Read bnds file with q-dependent number of bands
Cu   07 Jun 08 printout 'mc' merge command for 3 color weights
Cu   01 May 07 plbnds used variables ymin,ymax when making fplot cmd
Cu   12 Jan 07 New switch -nocol to suppress using color weight info
Cu   06 Jun 06 can read files with color weights
C ----------------------------------------------------------------------
      implicit none
      integer nqmx,maxlin,louts,nouts
      parameter(nqmx=2000, maxlin=50, louts=180, nouts=80 )
      double precision xline(0:maxlin),xdum(2),ydum(2)
      double precision xq(0:nqmx),eq(nqmx)
      real(8),allocatable:: eb(:,:),ebw(:,:,:),eb2(:),ebw2(:,:)
      double precision dq,emin,ebot,ebotb,etopb,emax,htcm,
     .  etik,xmin,xmax,ql1,ql2,ql3,e,h,x,q1,q2,q3,eferm,efermm,efqp,widcm,
     .  etop,xq0,efnew,abot,atop,scl,escl,syma(4),clip(4),frmcs(0:3),fillcs(0:4),sgs(0:7),
     .  qm1,qm2,qm3,csw(16),xx(4),wscal(2),wshft(2),window(2)
      integer ib,jb,kb,lb,mb,nb,nbf,nbfil,nbm,iq,nq,nqsave,nlab,iq1,nq1,i,ifi,jfi,
     .  ifc,ifm,ifmc,ifw,iline,iline0,nline,nline0,i1mach,nql(maxlin),ib1,ib2,j,
     .  ltfrme(3),psym,iouts,iiouts,ipr,icolwt,ncolwt,ncolwtx,lcol4,ifmt,ifig
      integer ix(4)
      integer, parameter :: NULLI=-99999
      integer, allocatable :: iblst(:)
      logical ls1,ls2,lpad,lskip,lsp
      logical :: lwriteqp = .false.  ! T => write qp levels to bnds file (spectral mode)
      logical :: ldrawqp = .false.  ! T => draw qp levels on top of spectral function (spectral mode)
      logical :: lfs = .false.  ! T => spectral function in FS (contour) mode
      character dc*1,dcfp*1,title*120,xlabel*120,ylabel*120,font*19,ltype*120,
     .  ll*120,s*256,fpstrn*160,outs(nouts)*(louts),prmed*40,ss*80,smerge*256,
     .  xnfmt*20,ynfmt*20,datnam*10,fncol4*80,colab(4)*256
      character*10 lab(maxlin)
      integer fontsz,nplst,iplst(20),lgcmd,lenf,iarg,nblst,j1,j2,ltbold
      integer :: spintexture=0
      procedure(logical) :: cmdopt,a2bin,rdstrn
      procedure(integer) :: fopna,fopnx,fopng,awrite,iinear,a2vec,nargf,fxst,cmdoptswx,isw,wordsw
      procedure(real(8)) :: dlength

C --- Parse all arguments before last - : -- should be switches ---
      if (nargf() <= 1) goto 999
      call getarf(1,s)
      if (s(1:6) == '--help' .or. s(1:2) == '-h' .or. s(1:3) == '--h') goto 999
      ifm = 0
      ifw = 0
      lpad = .false.
      lsp = .false.
      smerge = ' '
      ltype = '-lt 1,bold=3,col=0,0,0'
      datnam = 'dat'
      ls1 = .false.; ls2 = .false.
      nblst = 0
      scl = 1
      escl = 1
      nplst = 0
      window = NULLI
      call dpzero(csw,16)

      do  iarg = 1, nargf()-2
        call getarf(iarg,s)

        if (s(1:1) /= '-') goto 999
        if (s(1:7) == '-merge=') smerge = s(8:)
        if (s(1:8) == '-mergep=') then
          smerge = s(9:)
          lpad = .true.
        elseif (s(1:6) == '--help' .or. s(1:2) == '-h') then
          goto 999
        elseif (s(1:6) == '-skip=')  then
          call mkils0(s(7:),nplst,iplst)
          if (nplst > 20) call rx('increase size of iplst')
          call mkilst(s(7:),nplst,iplst)
        elseif (s(1:4) == '-lt=')  then
          ltype = s
          ltype(1:4) = '-lt '
        elseif (s(1:5) == '-dat=')  then
          datnam = s(6:)
        elseif (s(1:6) == '-spin1') then
          ls1 = .true.
        elseif (s(1:6) == '-spin2') then
          ls2 = .true.
C       fplot switches
        elseif (s(1:6) == '-fplot')  then
          dcfp = s(7:7)
          fpstrn = s
          i = cmdoptswx('-fplot','spintexture','')
          if (i <= 0) then
            spintexture = 0
          else
            select case (s(i+11:i+11))
            case ('0'); spintexture = 3; i=i+12
            case ('q'); spintexture = 2; i=i+12
            case default; spintexture = 1; i=i+11
            end select
            ncolwtx = 0
            if (s(i:i) == '=') then
              if (index(s(i+1:i+2),'x') > 0) then
                ncolwtx = ncolwtx + 2; spintexture = spintexture + 10
              endif
              if (index(s(i+1:i+2),'y') > 0) then
                ncolwtx = ncolwtx + 2; spintexture = spintexture + 20
              endif
              if (index(s(i+1:i+2),'z') > 0) then
                ncolwtx = ncolwtx + 2; spintexture = spintexture + 40
              endif
            endif
            if (spintexture/10 == 0) ncolwtx = 3
          endif
!         if (cmdoptswx('-fplot','2sp','') > 0) then
          i = cmdoptswx('-fplot','lt=','')  ! Parse for line type
          if (i > 0) then
            call nwordg(s,0,dcfp//' ',1,i,j) ! i now terminates switch and must be preserved
            ltype(2:) = s(i:j)
            ltype(1:4) = '-lt '
          endif

          i = cmdoptswx('-fplot','ext=','')  ! Parse for line type
          if (i > 0) then
            i = i+4
            call nwordg(s,0,dcfp//' ',1,i,j) ! i now terminates switch and must be preserved
            datnam = s(i:j)
          endif

          call awrit0(' plbnds creating setup for fplot',s,len(s),0)
          if (spintexture > 0) then
            call awrit5('%a w/ spin texture, %i colors:'//
     .        '%?#n# (x,-x)##'//
     .        '%?#n# (y,-y)##'//
     .        '%?#n# (z,-z)##'//
     .        '%?#n# (x,y,z)##',
     .        s,len(s),0,ncolwtx,
     .        mod(mod(spintexture/10,10),2),
     .        mod(mod(spintexture/10,10),4)-mod(mod(spintexture/10,10),2),
     .        mod(mod(spintexture/10,10),8)-mod(mod(spintexture/10,10),4),
     .        spintexture)
          endif
          if (cmdoptswx('-fplot','ext=','') > 0)
     .    call awrit0('%a  extension='//trim(datnam),s,len(s),0)
          call info0(1,1,0,trim(s))

C       Special Spectral function mode
        elseif (s(1:3) == '-sp')  then
          ifmt = 1
          ifig = 1
          atop = 0
          abot = 0
          lsp = .true.
          allocate(iblst(1))
          dc = s(4:4)
          if (dc == ' ') cycle
C     ... Return here to resume parsing for arguments
          j2 = 4-1
   10     continue
          j2 = j2+1
          if (s(j2:j2) == dc) goto 10
          j1 = min(len(s),j2)
          call nwordg(s,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
            elseif (s(j1:j1+1) == 'fs')  then
              lfs = .true.
            elseif (s(j1:j1+3) == 'lst=')  then
              if (j1+4 > j2) call rx('plbnds: bad list, lst=..')
              call mkils0(s(j1+4:j2),nblst,iblst)
              if (nblst < 0) goto 997
              deallocate(iblst); allocate(iblst(nblst))
              call mkilst(s(j1+4:j2),nblst,iblst)
            elseif (s(j1:j1+3) == 'out=')  then
              ifig = 1
              if (s(j1+4:j1+4) == dc) goto 10
              i = j1+3
              if (a2vec(s,j2,i,2,dc//' ',2,2,1,ix,ifig) < 0) then
                ifig = -1
                if (s(j1:j1+6) == 'out=eps') ifig = 1
                if (s(j1:j1+6) == 'out=svg') ifig = 2
                if (s(j1:j1+6) == 'out=png') ifig = 3
              endif
              if (ifig < 0) goto 997
              call sanrg(.true.,ifig,1,3,'plbnds:','output format')
            elseif (s(j1:j1+6) == 'writeqp')  then
              lwriteqp = .true.
            elseif (s(j1:j1+3) == 'fmt=')  then
              ifmt = 1
              if (s(j1+4:j1+4) == dc) goto 10
              i = j1+3
              if (a2vec(s,j2,i,2,dc//' ',2,2,1,ix,ifmt) < 0) then
                ifmt = -1
                if (s(j1:j1+6) == 'fmt=gnu') ifmt = 1
                if (s(j1:j1+7) == 'fmt=lmgf') ifmt = 2
              endif
              if (ifmt < 0) goto 997
              call sanrg(.true.,ifmt,1,2,'plbnds:','fmt')

            elseif (s(j1:j1+4) == 'escl=')  then
              i = j1+4
              call rxx(a2vec(s,j2,i,4,dc//' ',2,2,1,ix,escl) < 0,
     .          'failed to parse '//trim(s(j1:)))
            elseif (s(j1:j1+4) == 'ascl=')  then
              i = j1+4
              call rxx(a2vec(s,j2,i,4,dc//' ',2,2,1,ix,scl) < 0,
     .          'failed to parse '//trim(s(j1:)))
            elseif (s(j1:j1+6) == 'window=') then
              i = j1+6
              call rxx(a2vec(s,j2,i,4,dc//', ',2,3,2,ix,window) /= 2,
     .          'failed to parse '//trim(s(j1:)))
            elseif (s(j1:j1+4) == 'atop=')  then
              i = j1+4
              call rxx(a2vec(s,j2,i,4,dc//' ',2,2,1,ix,atop) < 0,
     .          'failed to parse '//trim(s(j1:)))
            elseif (s(j1:j1+4) == 'abot=')  then
              i = j1+4
              call rxx(a2vec(s,j2,i,4,dc//' ',2,2,1,ix,abot) < 0,
     .          'failed to parse '//trim(s(j1:)))
            elseif (s(j1:j1+5) == 'drawqp')  then
              ldrawqp = .true.; ltbold = 3; efqp = 0
              dcfp = s(j1+6:j1+6)
              if (.not. (dcfp == dc .or. dcfp == ' ')) then
                if (wordsw(s(j1+6:),dcfp,'ef=','',i) > 0) then
                  i = j1+4+i
                  call rxx(a2vec(s,j2,i,4,dc//dcfp//' ',3,1,1,ix,efqp)<0,
     .              'failed to parse '//trim(s(j1:)))
                endif
                if (wordsw(s(j1+6:),dcfp,'bold=','',i) > 0) then
                  i = j1+4+i
                  call rxx(a2vec(s,j2,i,2,dc//dcfp//' ',3,1,1,ix,ltbold)<0,
     .              'failed to parse '//trim(s(j1:)))
                endif
                if (wordsw(s(j1+6:),dcfp,'col=','',i) > 0) then
                  i = j1+4+i
                  call rxx(a2vec(s,j2,i,4,dc//dcfp//', ',4,5,3,ix,csw(2))<0,
     .              'failed to parse '//trim(s(j1:)))
                  csw(1) = 1
                endif
              endif
            endif
            goto 10
          endif
        elseif (s(1:3) == '-ef')  then
        elseif (s(1:4) == '-scl') then
        elseif (s(1:5) == '-wscl') then
        elseif (s(1:5) == '-fast') then
        elseif (s(1:4) == '-lbl') then
        elseif (s(1:3) == '-tl') then
        elseif (s(1:6) == '-fplot') then
        elseif (s(1:6) == '-ascii') then
        elseif (s(1:5) == '-col4' .or. s(1:6) == '--col4') then
        elseif (s(1:6) == '-nocol' .or. s(1:7) == '--nocol') then
        else
          print *, 'plbnds: switch "',trim(s),'" not recognized'
        endif

      enddo

C ... Sanity checks
      if (ls1 .and. ls2) call rx('plbnds: only spin1 or spin2 allowed')
      if (spintexture /= 0 .and. ls2) call rx('plbnds: spintexture and -spin2 incompatible')
      if (spintexture /= 0 .and.  cmdoptswx('-col4','','') > 0)
     .  call rx('plbnds: spintexture and -col4 not allowed')

C --- Other defaults ---
      xnfmt = '%;4d'
      ynfmt = '%;4d'
      sgs(0) = 2
      sgs(1) = 0
      sgs(2) = 0
      sgs(3) = 0
      sgs(4) = -1
      sgs(5:7) = 0
      frmcs(0:3) = [3d0,0d0,0d0,0d0]
      fillcs(0:4) = [103d0,1d0,1d0,1d0,1d0]
      iouts = 1
      call finits(2,0,0,i)
      forall (i = 1:nouts) outs(i) = ' '

C  18 continue
C ... Last argument is file name : open file
C     First check whether name exists
      call getarf(nargf()-1,s)
      j = fopnx(s,172,-1,-1)
C     Not merely extension if s contains '.'
      if (index(s,'.') /= 0) j = 1
      if (j == 0 .and. lsp) then
        call rx('plbnds requires true file name in sp mode')
      elseif (j /= 0) then
        ifi = fopng(s,-1,1)
C     Else check whether extension exists
      else
        j = fxst('bnds')
        if (j == 0) then
          call strip(s,i,j)
          print 345, s(i:j), s(i:j)
  345     format(' plbnds: failed to open file "',
     .      a,'" or file "bnds.',a,'".')
          goto 999
        endif
        ifi = fopna('bnds',-1,0)
      endif

      title = ' '
      nlab = 0
      forall (i=1:maxlin) lab(i) = ' '
      if (prmed /= '-ascii' .and.
     .    (cmdopt('-lbl',4,0,s) .or. cmdopt('-tl',3,0,s))) then
        if (cmdopt('-tl=',4,0,s)) then
          title = s(5:)
        elseif (.not. cmdopt('-lbl=',5,0,s)) then
C         call cwrite(' Title ? ',0,28,0)
          print 432, ' Title ? '
  432     format(a,1x,$)
          read(*,'(a80)') title
        endif
        s = ' '
        if (cmdopt('-lbl=',5,0,s)) then
          call wrdsg(s(6:len(s)),0,',',nlab)
          if (s(6:len(s)) == ' ') nlab = 0
          do  46  j = 1, nlab
          call wordg(s(6:len(s)),0,',',j,ib1,ib2)
   46     lab(j) = s(5+ib1:5+ib2)
        elseif (cmdopt('-lbl',4,0,s)) then
C         call cwrite(' Symmetry labels ? ',0,28,0)
          print 432, ' Symmetry labels ? '
          read(*,'(a80)') s
          call words(s,nlab)
          do  45  j = 1, nlab
          call word(s,j,ib1,ib2)
   45     lab(j) = s(ib1:ib2)
        endif
      endif

C ... Spectral Function mode
      if (lsp) then
        ss = 'spf'
        j = index(s,'.')
        if (j > 0) ss = 'spf'//s(j:)
        if (trim(ss) == trim(s)) ss = 'spf2'//s(j:)
        jfi = fopng(trim(ss),-1,0)
        title = ' spectral mode, write data to file '//trim(ss)
        if (nblst > 0) call awrit2('%a bands%n:1i',title,len(title),0,nblst,iblst)
        call awrit2('%a, %?#n==1#gnuplot##%-1j%?#n==2#lmgf## format',title,len(title),0,ifmt,0)
        call info0(1,0,0,trim(title))
        call suse(lfs,nblst,iblst,escl,scl,Abot,Atop,window,ifmt,lwriteqp,ifig,ifi,jfi,
     .    ldrawqp,ltbold,csw,efqp,lab,nlab)
        call rx0('done')
      endif

C ... Open merge file, if specified, file unit ifm
C     First check whether name exists
      if (smerge /= ' ') then
        call wordg(smerge,0,', ',1,i,j)
        if (smerge(i:j) /= smerge) then
          ss = smerge(i:j)
          call wordg(smerge,0,', ',2,i,j)
          ifw = fopng(smerge(i:j),-1,2)
          smerge = ss
        endif
        if (smerge == s)
     .      call rx('merge and original file of same name')
        s = smerge
        j = fopnx(s,172,-1,-1)
        if (j /= 0) then
          ifm = fopng(s,-1,1)
C       Else check whether extension exists
        else
          s = 'bnds.'//smerge
          j = fopnx(s,172,-1,-1)
          if (j == 0) then
            call strip(smerge,i,j)
            print 345, smerge(i:j), smerge(i:j)
            goto 999
          endif
          ifm = fopng(s,-1,1)
        endif
      endif

C ... Open 3rd color file
      lcol4 = 0
      ifc = 0
      if (cmdopt('-col4:',6,0,s) .or. cmdopt('--col4:',7,0,s)) then
        fncol4 = s(7:)
        if (s(1:2) == '--') fncol4 = s(8:)
        lcol4 = 1
      endif
      if (lcol4 == 1) then
        if (ifm /= 0) call rx('plbnds: cannot use -merge and -col4 together')
        call wordg(fncol4,0,', ',1,i,j)
        ib = fopnx(fncol4(i:j),172,-1,-1)
C       Not merely extension if s contains '.'
        if (index(fncol4(i:j),'.') /= 0) ib = 1
        if (ib /= 0) then
          ifc = fopng(fncol4(i:j),-1,1)
C       Else check whether file+extension exists
        else
          call rxs('plbnds: failed to open col4 file ',fncol4)
        endif
C       See if merge col4 file
        call wordg(fncol4,0,', ',2,i,j)
        ifw = fopng(fncol4(i:j),-1,0)
        rewind ifw
      endif

C      j = fxst('bnds')
C      if (j == 0) goto 999
C      ifi = fopna('bnds',-1,0)
      htcm = 16d0
C     gtoxcm = 6d0
      xlabel = ' '
      prmed = ' ps'
      psym = 6
      scl = 1
      syma(1) = .2d0
      syma(2) = 0
      syma(3) = 360
      if (cmdopt('-ascii',6,0,prmed)) then
        print *, 'writing ascii'
        psym = 1
      endif
      call getpr(ipr)
      ylabel = 'Energy (Ry)'
      font = '/Times-Roman'
      fontsz = 24
      ltfrme(1:3) = [3,0,0]

C --- Band file header, number of bands, Fermi level, color weights ---
C     For backwards compatibility, read 3rd argument as real; if so, icolwt=0
C     smerge used as a temporary character variable
      rewind ifi
      if (.not. rdstrn(ifi,smerge,len(smerge),.false.)) goto 999
      call words(smerge,i)
      icolwt = 0; lcol4 = 0
      xx(3) = 0
      j = 0; j = a2vec(smerge,len(smerge),j,4,', ',2,-3,min(i,3),ix,xx)
      if (j /= min(i,3)) call rx('plbnds: failed to read header in bands file')
      nbf = xx(1); eferm = xx(2)
      if (nbf /= xx(1)) call rx('plbnds: 1st arg in band file must be integer (no. of bands)')
C     For color weights, icolwt must be a positive integer
      if (abs(nint(xx(3))-xx(3)) == 0) icolwt = nint(xx(3))
      if (icolwt < 0) icolwt = 0
      if (spintexture /= 0 .and. icolwt /= 4) call rx('plbnds must have 4 color weights for spin texture')
      call word(smerge,3,i,j); i = i-1; call chrpos(smerge,'.',j,i); if (i < j) icolwt = 0
      kb = 0
      do  ib = 1, icolwt
        kb = kb+1; call word(smerge,3+kb,i,j); i = i + index(smerge(i:j),'=')
        if (i > j) then
          kb = kb+1; call word(smerge,3+kb,i,j)
        endif
        if (i > j) call info2(1,0,0,' plbnds (warning) '//
     .    'could not find color weight header col %i',ib,2)
        colab(ib) = smerge(i:j)
        if (spintexture /= 0) exit
      enddo

      if (icolwt >= 2) then
        call info2(1,0,0,' plbnds : bands file contains '//
     .    '%?#n==2#two##%-1j%?#n==3#three##%-1j%?#n==4#four##'//
     .    ' sets of color weights',icolwt,0)
      elseif (icolwt /= 0) then
        call info0(1,0,0,' plbnds : bands file contains color weights')
      endif
      if (lcol4 /= 0) then
        call info0(1,0,0,' plbnds : merge third color weight')
      endif
      if (ifc /= 0) then
      if (.not. rdstrn(ifc,smerge,len(smerge),.false.)) goto 999
      call words(smerge,i)
      rewind ifc
      if (i < 3) then
        if (i /= nbf) call rx('col4 file has no color weights')
      else
        read(ifc,*,end=999,err=999) i,h,x
        if (i /= nbf) call rx('col4 file has wrong number of bands')
        if (h /= eferm) call info0(0,0,0,
     .    ' plbnds (warning): col4 file fermi level mismatch')
C       For color weights, icolwt must be a positive integer
        lcol4 = nint(x)
        if (abs(nint(x)-x) /= 0)
     .    call rx('col4 file has no color weights')
        if (lcol4 /= 1 .and. lcol4 /= 2) call rxi
     .    ('col4 file must have one or two color weights, but has',lcol4)
        kb = 0
        do  ib = 1, lcol4
          kb = kb+1; call word(smerge,3+kb,i,j); i = i + index(smerge(i:j),'=')
          if (i > j) then
            kb = kb+1; call word(smerge,3+kb,i,j)
          endif
          colab(icolwt+ib) = smerge(i:j)
        enddo
      endif
      endif
      ncolwt = icolwt + lcol4
      if (ncolwt > 4) call rx(' plbnds: too many color weights')
      if (cmdopt('-nocol',6,0,s) .or. cmdopt('--nocol',7,0,s)) then
        ncolwt = 0
        call info0(0,0,0,' plbnds : color weights suppressed')
      endif

      allocate(eb(nbf,nqmx),eb2(nbf))
      if (icolwt /= 0) then
        allocate(ebw(nbf,nqmx,max(ncolwt,1)),ebw2(nbf,max(ncolwt,1)))
      endif

      nb = nbf
      if (ifm /= 0) then
        read(ifm,*,end=999,err=999) nbm,efermm
        if (nbm /= nbf)
     .   call info2(0,0,0,' plbnds (warning) : original file has %i '//
     .    'bands while merged file has %i bands',nbf,nbm)
        nb = min(nbm,nbf)
        if (lpad) nb = max(nbm,nbf)
        if (eferm /= efermm)
     .   call info2(0,0,0,' plbnds (warning) : original file has ef='//
     .    '%d while merged file has ef=%d',eferm,efermm)
        if (ls1 .or. ls2)
     .    call rx('-merge not compatible with -spin1 or -spin2')
      endif
      if (ifw /= 0) then
        if (ncolwt == 0) then
          write(ifw,'(i5,f10.5,i5)') nb,eferm,0
        else
          write(ifw,876) nb,eferm,ncolwt,(trim(colab(ib)),ib=1,ncolwt)
  876     format(i5,f10.5,i5:,
     .      2x,'col= ',a:2x,'col2= ',a:2x,'col3= ',a:2x,'col4= ',a)
        endif
      endif

      nq = 0
      nline = 0    ! Actual number of lines to draw
      nline0 = 0   ! Number of lines contained in file
      xq(0) = 0d0
      ebot = 1d10
      etop = -1d10
      xq0 = 0
      xline(0) = 0
      efnew = eferm
      wscal = 1; wshft = 0
      if (cmdopt('-ef=',4,0,s)) then
        j = 4
        call rxx(.not.a2bin(s,efnew,4,0,' ',j,len(s)),'failed to parse '//s)
      endif
      if (cmdopt('-scl=',5,0,s)) then
        j = 5
        call rxx(.not.a2bin(s,scl,4,0,' ',j,len(s)),'failed to parse '//s)
      endif
      if (cmdopt('-wscl=',6,0,s)) then
        j = 6; j = a2vec(s,len(s),j,4,', ',2,-2,2,ix,wscal)
        if (j < 1) call rx('failed to parse '//s)
        if (j == 1) wscal(2) = wscal(1)
      endif
      if (cmdopt('-wshft=',7,0,s)) then
        j = 7; j = a2vec(s,len(s),j,4,', ',2,-2,2,ix,wshft)
        if (j < 1) call rx('failed to parse '//s)
      endif
      call awrit2(' plbnds: %i bands  fermi=%;6d',s,80,0,nb,eferm)
      if (scl /= 1) call awrit1('%a  scaled by %;6d',s,80,0,scl)
      if (efnew /= eferm) call awrit1('%a  shifted to %;6d',s,80,0,efnew)
      call awrit0('%a',s,80,-i1mach(2))

C --- For each panel, do ---
      ifmc = ifm ! ifmc is ifm or ifc, if either is nonzero
      if (ifc /= 0) ifmc = ifc
   91 continue
      ebotb = 1d10
      etopb = -1d10
      read(ifi,*) nq1
C     Check for line match with merge file, if one is to be read
      if (ifmc /= 0) then
        read(ifmc,*) j
        if (nq1 /= j)
     .   call info5(0,0,0,' plbnds (fatal) : mismatch panel %i: file'//
     .    ' has %i qp while merged file has %i qp',nline0+1,nq1,j,0,0)
        if (nq1 /= j) goto 995
        if (ifmc == ifm) nq1 = nq1*2
      endif
C     If new merge file is to be created, write number of points on line
      if (ifw /= 0) then
        write(ifw,'(i6)') nq1
      endif

C     Exit read or possibly set skip for this line
      lskip = .false.
      if (nq1 <= 0) goto 90   ! No more lines ... exit read
      if (nplst /= 0) then
        ib = iinear(nplst,nline0+1,iplst,1)
        lskip = iplst(ib) == nline0+1
      endif

      if ((ls1 .or. ls2) .and. mod(nq1,2) /= 0) then
        call setpr(30)
        call rx('plbnds: spin pol case needs even number of bands')
      endif
      if (ls1 .or. ls2 .or. ifm /= 0) nq1 = nq1/2
C     Read energy bands along this line from bands file
      nqsave = nq
      do  iq1 = 1, nq1
        nq = nq+1
        if (nq > nqmx) stop 'increase nqmx in plbnds'
        nbfil = 0
        call rdbnd(ifi,q1,q2,q3,nbf,nbfil,-NULLI,eb(1,nq))
C       read(ifi,500) q1,q2,q3,(eb(i,nq),i=1,nbf)
        if (icolwt /= 0) then
C        Read associated file weights
          do  j = 1, icolwt
            if (ncolwt > 0 .or. j == 1) then
              call rdbnd(ifi,qm1,qm2,qm3,nbf,nbfil,0,ebw(1,nq,j))
            else
              call rdbnd(ifi,qm1,qm2,qm3,nbf,nbfil,0,ebw2)
            endif
C           read(ifi,500) qm1,qm2,qm3,(ebw(i,nq,j),i=1,nbf)
            if (qm1 /= q1) goto 995
            if (qm2 /= q2) goto 995
            if (qm3 /= q3) goto 995
          enddo
        endif
C       Read corresponding energy bands from merge file
        if (ifm /= 0) then
          read(ifmc,*) qm1,qm2,qm3,(eb2(i),i=1,nbm)
          if (qm1 /= q1) goto 995
          if (qm2 /= q2) goto 995
          if (qm3 /= q3) goto 995
C       If only one spin is to be extracted, read the second spin
        elseif (ifc /= 0) then
          nbfil = 0
          call rdbnd(ifc,qm1,qm2,qm3,nbf,nbfil,-NULLI,eb2)
          if (qm1 /= q1) goto 995
          if (qm2 /= q2) goto 995
          if (qm3 /= q3) goto 995
          if (abs(sum(eb(:,nq)-eb2(:))) /= 0) call info2(0,0,0,
     .      ' plbnds (warning) : bands mismatch, iq=%i',iq,nline0)
C         Read associated file weights
          do  j = 1, lcol4
            call rdbnd(ifc,qm1,qm2,qm3,nbf,nbfil,0,ebw2(1,j))
            if (qm1 /= q1) goto 995
            if (qm2 /= q2) goto 995
            if (qm3 /= q3) goto 995
            ebw(:,nq,icolwt+j) = ebw2(:,j)
          enddo

        else
          if (ls1 .or. ls2) then
            nbfil = 0
            call rdbnd(ifi,q1,q2,q3,nb,nbfil,-NULLI,eb2)
C           read(ifi,500) q1,q2,q3,(eb2(i),i=1,nb)
            if (icolwt /= 0) then
              do  j = 1, icolwt
                i = j
                if (ncolwt == 0) i = 1
                call rdbnd(ifi,q1,q2,q3,nb,nbfil,0,ebw2(1,i))
              enddo
            endif
          endif
        endif
        if (lpad) then
          do  i= nbf+1,nb
            eb(i,nq) = -NULLI
          enddo
          do  i= nbm+1,nb
            eb2(i) = -NULLI
          enddo
        endif
C       If new merge file is to be created, write the bands
        if (ifw /= 0 .and. ifc /= 0) then
          write(ifw,500) q1,q2,q3,(eb(i,nq),i=1,nb)
          do  j = 1, ncolwt
            write(ifw,500) q1,q2,q3,(ebw(i,nq,j),i=1,nb)
          enddo
        elseif (ifw /= 0) then
          write(ifw,500) q1,q2,q3,(eb(i,nq),i=1,nb)
          write(ifw,500) q1,q2,q3,(eb2(i),i=1,nb)
        endif
C       If 2nd spin is to be extracted, poke bands into spin1
        if (ls2) then
          call dcopy(nb,eb2,1,eb(1,nq),1)
          if (icolwt /= 0 .and. ncolwt > 0) then
            do  j = 1, icolwt
              call dcopy(nb,ebw2(1,j),1,ebw(1,nq,j),1)
            enddo
          endif
        endif
  500   format(3f10.5/(10f8.4))
  501   format((10f8.4))

C   ... Spin texture: redefine color weights
        if (spintexture /= 0) then
          do  ib = 1, nb
C           Scaling
            select case (mod(spintexture,10))
              case (1); xx(4) = dlength(3,ebw(ib,nq,2:4),1)
              case (2); xx(4) = ebw(ib,nq,1)
              case default; xx(4) = 1
            end select
            if (xx(4) == 0) xx(4) = 1 ! If no weights in x,y,z
            xx(1:3) = ebw(ib,nq,2:4)/xx(4)
C           Alternate choice of weights
            i = 1
            if (mod(mod(spintexture/10,10),2) == 1) then
              ebw(ib,nq,i) = max(xx(1),0d0); ebw(ib,nq,i+1) = max(-xx(1),0d0); i = i+2
            endif
            if (mod(mod(spintexture/10,10),4) >= 2) then
              ebw(ib,nq,i) = max(xx(2),0d0); ebw(ib,nq,i+1) = max(-xx(2),0d0); i = i+2
            endif
            if (mod(mod(spintexture/10,10),8) >= 4) then
              ebw(ib,nq,i) = max(xx(3),0d0); ebw(ib,nq,i+1) = max(-xx(3),0d0); i = i+2
            endif
C           Default weights
            if (mod(spintexture/10,10) == 0) then
              ebw(ib,nq,1:3) = abs(xx(1:3))
            endif
          enddo
        endif

C       Shift and scale bands
        forall (i=1:nb) eb(i,nq) = eb(i,nq)*scl + efnew - eferm*scl
C       Find upper and lower bounds (global and for this panel)
        if (.not. lskip) then
          do  i = 1, nb
            ebot = dmin1(ebot,eb(i,nq))
            etop = dmax1(etop,eb(i,nq))
            ebotb = dmin1(ebotb,eb(i,nq))
            etopb = dmax1(etopb,eb(i,nq))
          enddo
        endif
        dq = 0
        if (iq1 > 1) dq = dsqrt((q1-ql1)**2+(q2-ql2)**2+(q3-ql3)**2)
        xq(nq) = xq(nq-1)+dq
        if (ifm /= 0) then
          xq(nq+1) = xq(nq)
          nq = nq+1
        endif
        if (ipr > 40)
     .  call awrit7(' panel %i iq=%i  dq=%;6d (%;6d %;6d %;6d) xq=%;6d',
     .   ' ',80,i1mach(2),nline0,nq,dq,(q1-ql1),(q2-ql2),(q3-ql3),xq(nq))
        ql1 = q1
        ql2 = q2
        ql3 = q3
      enddo ! End of read block for energy bands for this line
      if (ifm /= 0) nq1 = nq1*2

      nline0 = nline0+1
      if (lskip) then
        call info2(2,0,0,' panel %i  nq=%i  skipped ...',nline0,nq1)
        nq = nqsave
        goto 91  ! next panel
      endif
      nline = nline+1
      nql(nline) = nq1
      xline(nline) = xq(nq)
      call awrit5(' panel %i  nq=%i  ebot=%1;6d  etop=%1;6d'//
     .  '  delta q=%1;6d',' ',80,i1mach(2),nline0,nq1,ebotb,etopb,
     .  xq(nq)-xq0)
      xq0 = xq(nq)
      goto 91  ! next panel

C --- Clean up setup ---
   90 continue
C      nline0 = nline0-1
C      if (.not. lskip) nline = nline-1
      if (ifw /= 0) then
        call fclose(ifw)
        call rx0('merged file written')
      endif
      if (nq == 0) call rx0('no bands to draw')

C --- Get dimensions for plot ---
      xmin = xq(1)
      xmax = xq(nq)
      call awrit6(' nq=%i  npan=%i  emin=%1;6d  '//
     .  'ef=%1;6d  emax=%1;6d  sum dq=%1;6d',
     .  ' ',80,i1mach(2),nq,nline,ebot,efnew,etop,xmax-xmin)
      emin = ebot
      emax = etop
C     gtoxcm = 5
      htcm = 10
   93 print '('' emin, emax, width(cm), height(cm) ?'')'
      read(5,*)  emin,emax,widcm,htcm
      if (widcm+2d0 > 25d0) then
        call awrit1(' Total width (%1;1dcm) cannot exceed 23cm',
     .    ' ',80,i1mach(2),widcm)
        goto 93
      endif
      call dcopy(4,0d0,0,clip,1)
      clip(2) = widcm/14d0
      clip(4) = htcm/14d0

C --- Make data files in format suitable for fplot ---
      if (cmdopt('-fplot',6,0,s)) then
        lgcmd = fopng('plot.plbnds',-1,0)

C   ... Scale and/or shift
        j = cmdoptswx('-fplot','scl=','') + 3
        if (j > 3) then
          j = a2vec(fpstrn,len_trim(fpstrn),j,4,', '//dcfp,3,-2,2,ix,wscal)
          if (j < 1) call rxs('failed to parse string ',fpstrn)
          if (j == 1) wscal(2) = wscal(1)
        endif
        j = cmdoptswx('-fplot','shft=','') + 4
        if (j > 4) then
          j = a2vec(fpstrn,len_trim(fpstrn),j,4,', '//dcfp,3,-2,2,ix,wshft)
          if (j < 1) call rxs('failed to parse string ',fpstrn)
          if (j == 1) wshft(2) = 0
        endif

C   ... Assign line type to character variable
        if (ltype /= ' ') then
        write(lgcmd,'(''%char0 ltb="'',a,''"'')') trim(ltype(5:))
        ltype(5:) = '{ltb}'
        endif
        call awrit2('fplot%N%% var ymin=%;3d ymax=%;3d%N',outs(iouts),louts,0,emin,emax)
        nq1 = 0
        iline = 0
        do  iline0 = 1, nline0
          if (nplst /= 0) then  ! possibly skip this panel
            ib = iinear(nplst,iline0,iplst,1)
            if (iplst(ib) == iline0) cycle
          endif
          iline = iline+1
          call awrit4('%a -frme %;3d,%;3d,%;3d,%;3d -p0 -x 0,1 -y {ymin},{ymax} -noxn -1p',
     .      outs(iouts),louts,0,
     .      wshft(1)+wscal(1)*xline(iline-1)/(xmax-xmin),
     .      wshft(1)+wscal(1)*xline(iline)/(xmax-xmin),
     .      wshft(2),
     .      wshft(2)+wscal(2))
          j = cmdoptswx('-fplot','ts=','') + 2
          if (j > 2) then
            j = a2vec(fpstrn,len_trim(fpstrn),j,4,', '//dcfp,3,-2,2,ix,xx)
            if (j < 1) call rxs('failed to parse string ',fpstrn)
            call awrit3('%a -tmx %d%?#(n>1)#;%d##',outs(iouts),louts,0,xx,j,xx(2))
          endif
          if (title /= ' ') then
            call strip(outs(1),ib1,ib2)
            outs(1)(ib2+1:) = ' -tl "' // title
            call strip(outs(1),ib1,ib2)
            outs(1)(ib2+1:ib2+1) = '"'
            title = ' '
          endif
C         make ib1,ib2 = indices to lowest, highest band for this panel
          ib1 = nb+1
          ib2 = -1
!          else
          do  ib = 1, nb
            ebotb = 1d10
            etopb = -1d10
            do  iq = 1+nq1, nql(iline)+nq1
              ebotb = dmin1(ebotb,eb(ib,iq))
              etopb = dmax1(etopb,eb(ib,iq))
            enddo
            if (ebotb < emax .and. etopb > emin) then
              ib1 = min(ib1,ib)
              ib2 = max(ib2,ib)
            endif
          enddo
          j = cmdoptswx('-fplot','ib1=','') + 3
          if (j > 3) then
            if (a2vec(fpstrn,len_trim(fpstrn),j,2,dcfp//' ',2,3,1,ix,ib1) < 0)
     .        call rxs('failed to parse string ',fpstrn)
          endif
          j = cmdoptswx('-fplot','ib2=','') + 3
          if (j > 3) then
            if (a2vec(fpstrn,len_trim(fpstrn),j,2,dcfp//' ',2,3,1,ix,ib2) < 0)
     .        call rxs('failed to parse string ',fpstrn)
          endif

C         Open band file for this panel
          call awrit3(' write file bnd%i.'//datnam//'%a, bands %i - %i',
     .      ' ',80,i1mach(2),iline0,ib1,ib2)
          lenf = awrite('%xbnd%i.'//datnam//'%a',s,len(s),0,iline0,0,0,0,0,0,0,0)
          open(10, file=s(1:lenf), status='UNKNOWN', err=99)

C         Two distinct spin channels: write them both out
          if (spintexture == 0) ncolwtx = ncolwt  ! ncolwtx already set if spintexture
          if (cmdoptswx('-fplot','2sp','') > 0) then
            if (mod(nql(iline),2) /= 0)
     .        call rx('spin pol bands should have even nq')
            j = 2*(ib2-ib1+1)
            if (ncolwtx /= 0) j = j*(ncolwtx+1)
            call awrit2('%% rows %i cols %i',' ',80,10,nql(iline)/2,1+j)
            do  iq = 1+nq1, nql(iline)+nq1,2
              write(10,333)
     .          (xq(iq)-xq(1+nq1))/(xq(nql(iline)+nq1)-xq(1+nq1)),
     .          (eb(ib,iq),eb(ib,iq+1),ib=ib1,ib2)
              if (ncolwtx /= 0) then
                do  j = 1, ncolwtx
                  write(10,334) (ebw(ib,iq,j),ebw(ib,iq+1,j),ib=ib1,ib2)
                enddo
              endif
            enddo

C         Write out bands and optional weights
          else
            j = ib2-ib1+1
            if (ncolwtx /= 0) j = j*(ncolwtx+1)
            call awrit2('%% rows %i cols %i',' ',80,10,nql(iline),1+j)
            if (maxval(eb(ib1:ib2,1+nq1:nql(iline)+nq1)) >= 1000 .or.
     .          minval(eb(ib1:ib2,1+nq1:nql(iline)+nq1)) <= -100) then
              ifmt = 3
            else
              ifmt = 0
            endif
            do  iq = 1+nq1, nql(iline)+nq1
              if (ifmt == 0) then
              write(10,333)
     .          (xq(iq)-xq(1+nq1))/(xq(nql(iline)+nq1)-xq(1+nq1)),
     .          (eb(ib,iq),ib=ib1,ib2)
              else
              write(10,335)
     .          (xq(iq)-xq(1+nq1))/(xq(nql(iline)+nq1)-xq(1+nq1)),
     .          (eb(ib,iq),ib=ib1,ib2)
              endif
              if (ncolwtx /= 0) then
                do  j = 1, ncolwtx
                  if (ifmt == 0) then
                  write(10,334) (ebw(ib,iq,j),ib=ib1,ib2)
                  else
                  write(10,336) (ebw(ib,iq,j),ib=ib1,ib2)
                  endif
                enddo
              endif
            enddo
          endif
  333     format(f8.5,9f9.4:/(8x,9f9.4))
  334     format((8x,9f9.4))
  335     format(f8.5,9f10.4:/(8x,9f10.4))
  336     format((8x,9f10.4))
          close(10)
          nq1 = nq1 + nql(iline)
          if (iline > 1) call awrit1('%a -noyn',outs(iouts),louts,0,1)
          iouts = iiouts(outs,iouts,louts+64)
          if (iouts > nouts) call rx('plbnd output string too long')
          ib = 1+(ib2-ib1+1)

          if (cmdoptswx('-fplot','2sp','') > 0) ib = 1+2*(ib2-ib1+1)
          jb = 1+2*(ib2-ib1+1)
          if (cmdoptswx('-fplot','2sp','') > 0) jb = 1+4*(ib2-ib1+1)
          kb = 1+3*(ib2-ib1+1)
          lb = 1+4*(ib2-ib1+1)
          mb = 1+5*(ib2-ib1+1)
          if (cmdoptswx('-fplot','2sp','') > 0) kb = 1+6*(ib2-ib1+1)

          call awrit2('%a -lt 2,bold=3,col=0,0,0,1,.5,.3,.5 '//
     .      '-tp 2~0,%1;4d,1,%1;4d -lt 1',
     .      outs(iouts),louts,0,efnew,efnew)
          iouts = iiouts(outs,iouts,louts+90)
          if (iouts > nouts) call rx('plbnd output string too long')
          if (ncolwtx == 4) then
            call awrit8('%a -colsy 2:%i'//
     .        ' -colsw %i:%i -colsw2 %i:%i -colsw3 %i:%i',
     .        outs(iouts),louts,0,ib, ib+1,jb, jb+1,kb, kb+1,lb,0)
            call awrit8('%a -colsw4 %i:%i '//ltype//
     .        '%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,
     .        0, lb+1,mb, iline0,0,0,0,0,0)
          elseif (ncolwtx == 3) then
            call awrit8('%a -colsy 2:%i'//
     .        ' -colsw %i:%i -colsw2 %i:%i -colsw3 %i:%i '//ltype//
     .        '%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,
     .        0,ib,ib+1,jb,jb+1,kb,kb+1,lb,iline0)
          elseif (ncolwtx == 2) then
            call awrit6('%a -colsy 2:%i'//
     .        ' -colsw %i:%i -colsw2 %i:%i '//ltype//
     .        '%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,
     .        0,ib,ib+1,jb,jb+1,kb,iline0)
          elseif (ncolwtx /= 0) then
            call awrit4('%a -colsy 2:%i'//
     .        ' -colsw %i:%i '//ltype//'%a -qr bnd%i.'//datnam//'%a',
     .        outs(iouts),louts,0,ib,ib+1,jb,iline0)
          else
            call awrit2('%a -colsy 2:%i '//
     .        ltype//'%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,0,
     .        ib,iline0)
          endif

C     ... Symmetry label, if there is one
          if (iline <= nlab) then
            iouts = iiouts(outs,iouts,louts+30)
            call awrit3('%a -font %?#n#s#t#24'//
     .        ' -lblu %;2d,{ymin-(ymax-ymin)/60} cd '//
     .        lab(iline),outs(iouts),louts,0,isw(lab(iline) == 'G'),
     .        0d0,emin-(emax-emin)/60)
          endif
          iouts = iiouts(outs,iouts,louts+10)
          if (iouts > nouts) call rx('plbnd output string too long')
        enddo
C   ... Last symmetry label, if there is one
        if (nline+1 <= nlab) then
          iouts = iiouts(outs,iouts,louts+30)
          call awrit3('%a -font %?#n#s#t#24'//
     .      ' -lblu %;2d,{ymin-(ymax-ymin)/60} cd '//
     .      lab(nline+1),outs(iouts),louts,0,isw(lab(nline+1) == 'G'),
     .      1d0,emin-(emax-emin)/60)
        endif
        do  25  ib = 1, iouts
          call skpblb(outs(ib),louts,i)
C         call cwrite(outs(ib),0,i,1)
          write(lgcmd,'(a)') outs(ib)(1:i+1)
   25   continue
        call fclose(lgcmd)
        if (cmdoptswx('-fplot','sh','') > 0) then
          call info0(2,1,0,' plbnds: invoking system call: fplot -pr10 -f plot.plbnds')
          call fsystm('fplot -pr10 -f plot.plbnds',i)
        else
          print *, ' ... to make a postscript file, invoke:'
          print *, ' fplot -f plot.plbnds'
        endif
        call cexit(0,1)
   99   call rx(' file open error, file='//s(1:lenf))
      endif

C --- Plot frame ---
      if (nplst /= 0) call rx('fplot not setup for -lst=')
      call pltini(prmed(2:10),0,.false.)
      call pltstp(clip(1),clip(2),clip(3),clip(4))
      call pltstu(xmin,xmax,emin,emax,0d0)
      etik = 0.1d0
      if (dabs(emax-emin) > 5d0) etik = 0.5d0
      if (dabs(emax-emin) > 10d0) etik = 1d0
      if (dabs(emax-emin) > 20d0) etik = 2d0
      call setfnt(font,fontsz)
      xx(1) = NULLI
      call frme(-1d30,xmax-xmin,1,0d0,   .6d0,0,0,xnfmt,xlabel,xx,
     .          -1d30,etik,     2,.025d0,.6d0,0,1,ynfmt,ylabel,xx,
     .          ' ',xdum,xdum,xdum,title,ltfrme,' ',frmcs,fillcs)
!fgs,fgs)
      do  4  iline = 1, nline+1
        call mve(xline(iline),emax)
        call drw(xline(iline),emin)
    4 continue

C --- Put in vertical lines and symmetry labels ---
      do  5  i = 0, nline
        ll = ' '
        ll = lab(i+1)
C       print *, ll
        call setfnt(font,fontsz)
        if (ll == 'G') call setfnt('/Symbol            ',fontsz)
        if (i == 0) then
          x = xmin
        elseif (i == nline) then
          x = xmax
        else
          x = xline(i)
        endif
        h = fontsz / (72 * 2.54 * 2) * scl
        call pstr(x,emin-h,0d0,ll,0d0,0d0,0d0,'cc','u',1)
C       call mve(x,emin)
C       call pslabl('a',0,0d0,ll,'h','f',1,0d0)
    5 continue

C --- Plot points ---
      do  7  ib = 1, nb
        if (mod(ib,10) == 1)
     .    call awrit1(' Plotting band %i',' ',80,i1mach(2),ib)
        do  8  iq = 1, nq
    8   eq(iq) = eb(ib,iq)
        if (cmdopt('-fast',5,0,s)) then
          call plsym(xq(1),eq,0d0,nq,psym,syma,sgs)
        else
          do  6  iq = 1, nq
            e = eb(ib,iq)
            if (e >= emin .and. e <= emax) then
              call plsym(xq(iq),e,0d0,1,psym,syma,sgs)
            endif
    6     continue
        endif
    7 continue

C --- Draw in Fermi energy ---
      if (efnew >= emin .and. efnew <= emax) then
        xdum(1) = xmin
        xdum(2) = xmax
        ydum(1) = efnew
        ydum(2) = efnew
        call plntyp(2,3,2d0,1d0,0d0,0d0)
        call plcrv(xdum,ydum,xdum,xdum,xdum,xdum,2,xdum(1),xdum(2),
     .    ydum(1)-0.1d0,ydum(2)+0.1d0,0,-2d0,csw,11)
      endif
      call pltdmp(0)
      call cexit(0,1)

C --- Error exit ---
  997 print *, 'plbnds: failed to parse "',s(j1:j2),'"'
      print *

  999 continue
      print '(a,a)',' usage: plbnds [-switches] filename'
      print 399, 'Last argument can be a true file name or "ext" with true filename = bnds.ext'
      print 399, 'Switches:'
  399 format(6x,a:t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,
     .          a:t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19,a:/t19)
      print 399, '--help | -h','','print this message and exit'
      print 399, '-ef=#','shifts fermi energy to #'
      print 399, '-scl=#','scales bands by #'
      print 399, '-spin1 | -spin2','','plots only spin-1 bands or spin-2 bands'
      print 399, '-lbl=CC..','specify symmetry point labels'
      print 399, '-lbl','plbnds prompts for symmetry labels and title'
      print 399, '-tl=string','specify title'
      print 399, '-skip=lst','skip panels in lst, e.g. -skip=1,3'
      print 399, '--nocol | -nocol','','suppress color weights'
      print 399, '-col4:fn2,outfile','','merges fn2, bnds file with 3rd or 3rd+4th color',
     .           'Example:  plbnds -col4:bnds.pd2,bnds.dat bnds.pd'
      print 399, '-merge=file2[,file3]'
      print 399, '-mergep=file2[,file3]','','reads two input files (one for each spin in spin-pol case)',
     .  'Optional file3 causes plbnds to write merged file to file3',
     .  '-mergep pads channel with fewer bands'
      print 399, '-sp[~switches]','','Band plot of spectral function.',
     .  'Switches: lst=list  out=eps|png|svg  writeqp  drawqp[:options] escl=#  ascal=#  atop=#  window=#,#',
     .  'See web site for documentation.'
      print 399, '-fplot[~options]','','generates an fplot script and supportng data files.'
      print 399, '','Options:','   ~2sp          writes data for two spins',
     .                         '   ~sh           invoke fplot in a subshell to make fplot.ps',
     .                         '   ~lt=string    specify fplot line type',
     .                         '   ~ext=nam      specify extension of data files',
     .                         '   ~ib1=#        specify lowest band to include in plot',
     .                         '   ~ib2=#        specify highest band to include in plot',
     .                         '   ~scl=#1[,#2]  scales graphics window by #1 or (#1,#2)',
     .                         '   ~shft=#1[,#2] shifts graphics window by #1 or (#1,#2)',
     .                         '   ~ts=#1[,#2]   sets tic spacing (ts=.5 => one tic), optional height'
      print 399, '',' ',       '   ~spintexture  spin texture defines color weights.',
     .                         '                 Default normalization : (|x,y,z)|=1',
     .                         '   ~spintexture0 no renormalization of x,y,z',
     .                         '   ~spintextureq normalization : (|x,y,z)|/w0=1',
     .                         '                 Default: col1=|x|, col2=|y|, col3=|z|.',
     .                         '                 Append =x | =y | =z to use 2 color weights',
     .                         '                 for one Cartesian component',
     .                         '                 (1st weight for >0, 2nd for <0)',
     .                         '                 One or two components can be specified.'
      print 399, '-lt=strn | -dat=nam | -wscl=#1[,#2] | -wshft=#1[,#2] (use with fplot)'
      print 399, '','same as fplot~lt=strn or fplot~ext=nam or fplot~scl=#1[,#2] or fplot~shft=#1[,#2]'
      if (s(1:2) == '-h' .or. s(1:6) == '--help') call cexit(0,1)
      call cexit(-1,1)

  995 continue
      call pshpr(40)
      call rx('mismatch between original and merged files')

      end
      integer function iiouts(outs,iln,itol)
      implicit none
      integer iln,itol,louts,i
      character*(*) outs(1)
      louts = len(outs(1))
      call skpblb(outs(iln),louts,i)
      if (i > louts-itol) then
C       outs(iln)(i+2:i+2) = '\\'
        iln = iln+1
      endif
      iiouts = iln
      end
      subroutine rdbnd(ifi,q1,q2,q3,nbf,nbfil,ipad,eb)
C- Read bands for one qp from standard bands file
C  and pads bands from nbfil .. nbf
Ci nbfil number of bands to read
Ci       input nbfil 0 => attempt to read from file.
Ci       If missing, attempt to read nbf bands.
Co       On output, nbfil = number actually read
      implicit none
      integer ifi,nbf,ipad,nbfil
      double precision q1,q2,q3,eb(nbf)
      character strn*120
      integer nw,j,j1,j2,ix(4),a2vec,i
      double precision fac(4)
      logical rdstrn

C     Check whether format contains number of bands to read
      if (.not. rdstrn(ifi,strn,len(strn),.false.)) goto 99
      call words(strn,nw)
      if (nw < 3) goto 99
      if (nw == 3) then
        if (nbfil == 0) nbfil = nbf
        j = 0
        j = a2vec(strn,len(strn),j,4,', ',2,-3,3,ix,fac)
        if (j /= 3) goto 99
      elseif (nw == 4) then
        call word(strn,4,j1,j2)
        j = 0
        j = a2vec(strn,len(strn),j,4,', ',2,-3,4,ix,fac)
        if (j /= 4) goto 99
        if (nbfil == 0) then
          nbfil = nint(fac(4))
          if (nbfil /= nint(fac(4))) then
            call rx('plbnds:  expected integer for 4th '//
     .        'argument, line: '//trim(strn))
          endif
        endif
      endif
      q1 = fac(1)
      q2 = fac(2)
      q3 = fac(3)

      read(ifi,501) (eb(i),i=1,nbfil)
  501 format((10f8.4))
      do  i = nbfil+1, nbf
        eb(i) = ipad
      enddo
      return
   99 continue
      call rx('plbnds:  failed to read qp from file, line:  '//
     .  trim(strn))

      end

      subroutine suse(lfs,nblst,iblst,Escl,Ascl,Abot,Atop,window,ifmt,lwriteqp,ifig,ifi,jfi,
     .  ldrawqp,ltbold,csw,efqp,lab,nlab)
C- Set up spectral function plots with dynamical self-energy
      implicit none
      integer maxlin
      parameter( maxlin=50 )
C ... Passed parameters
      logical lfs,lwriteqp
      integer nblst,ifmt,ifig,ifi,jfi,iblst(nblst),nlab,ltbold
      double precision Escl,Ascl,Abot,Atop,window(2),efqp,csw(0:3,1)
      logical ldrawqp
      character*10 lab(maxlin)
C ... Dynamically allocated arrays
      integer, allocatable :: ifblst(:),nqlst(:)
      real(8),allocatable:: qp(:,:),agq(:,:),qcut(:),qxlen(:),sag(:,:)
      real(8),allocatable:: eigiq(:,:)
C     real(8),allocatable:: sigiq(:,:),sexiq(:,:),vxciq(:,:)
C ... Local parameters
      real(8), parameter:: tolq=1d-7
      character spfile*120,fignam*7,prfmt*40,bndfilen*120,itos*15,spext*2,labstr*512,labf*20
      integer fmode,i,ib,icut,iq,iw,kb,lascii,kfi
      integer bnd,nband,ncut,nfblst,nincl,nomg,npan,nq,nqp,nsp,nspse,iqlst,bndfile
      double precision ommin,ommax,xx,distprev,dq(3),qdist,dw0,omgi,smax,Amax,chempot
      integer, parameter :: NULLI=-99999
      procedure(integer) :: iinear,fopng,fopnx
      procedure(real(8)) :: dlength,dsum

      spext = ''
      if (ldrawqp) spext = ',\'

      if (ifmt == 1) then
        kfi = fopng('gnu.plt',-1,0)
      endif

      nomg = 0

C ... Read header
      lascii = 20
      i = 3000 + 200 + lascii + 4
      call ioseh(i,fmode,nsp,nspse,nband,nqp,nomg,ommin,ommax,chempot,nfblst,[0],npan,ifi)
      if (npan == -3) then
        if (lfs) then           ! lfs requires k-points supplied on a raster
          if (nomg /= 1) call rx('Fermi surface mode expects single frequency')
        else
          call rx('3D spectral functions not ready: input contains q-points on a raster')
        endif
      endif
C     Require A(omega).  OK if file contains extra sexiq or vxciq
      if (mod(fmode/100,10) /= 1) call rx('file does not appear to contain spectral function')
C     if (mod(mod(fmode/10,10),2) /= 1) call rx('file does not appear to contain sigiq')
      if (mod(fmode,10) /= 2) call rx('file is not written in a format usable by plbnds')

      allocate(ifblst(nfblst),nqlst(0:max(2,npan)))
      i = i + 4000
      rewind ifi
      call ioseh(i,fmode,nsp,nspse,nband,nqp,nomg,ommin,ommax,chempot,nfblst,ifblst,nqlst,ifi)
C     nqlst(0) = 0

      if (nqp < 2) call rx0('plbnds: no bands to draw')

C     Scale energies
      ommin = ommin*escl
      ommax = ommax*escl

C ... Read q, eqp, sigqp. sigiq,sexiq,vxciq, A are not read
C     allocate(qp(3,nqp),eigiq(nfblst,nqp),sigiq(nfblst,nqp),sexiq(nfblst,nqp),vxciq(nfblst,nqp))
      allocate(qp(3,nqp),eigiq(nfblst,nqp))
      call iose2(lascii,fmode,0,1,nfblst,nqp,nomg,qp,eigiq,[xx],[xx],[xx],[xx],ifi)

      if (lfs) goto 40

C --- Write eigiq to bnds file ---
      if (lwriteqp) then
        spfile = 'bnds.dat'
        call info0(20,1,0,' Writing qp levels to file '//trim(spfile)//' ...')
        kfi = fopng(trim(spfile),-1,0)
        write(kfi,335) nfblst,chempot,0
  335   format(i5,f10.5,i6)
        iqlst = 0; nq = 0
        prfmt = '(3f10.5/(10f8.4))'
        do  iq = 1, nqp
          if (iq == nq+1 .or. nq == 0) then  ! new panel
            iqlst = iqlst+1
            write(kfi,337) nqlst(iqlst) - nq
  337       format(2i5)
            nq = nqlst(iqlst)
          endif
          write(kfi,prfmt) qp(:,iq), (eigiq(i,iq),i=1,nfblst)
        enddo
        write(kfi,337) 0
        call rx0('finished writing to bands file')
      endif

C --- Find the path cuts in the set of k-points (heuristic)
C     Clean this up; all we need is qxlen (it was adapted from gf/specfun.f)
      allocate(qxlen(nqp))
      icut = 1 ; qxlen = 0 ; distprev = 0
      do  iq  = 2, nqp
        dq = qp(:,iq)-qp(:,iq-1) ; qdist = dlength(3,dq,1)
        if (qdist < 1d-6) then
          icut = 1 ; qxlen(iq) = qxlen(iq-1) + qdist
        elseif (qdist > distprev*1.01 .and. icut == 0) then
          qxlen(iq) = qxlen(iq-1) ! + 0.02
          icut = 1
        else
          qxlen(iq) = qxlen(iq-1) + qdist ; icut = 0
        endif
        distprev = qdist
      enddo

C     Write header to jfi
      allocate(qcut(0:nqlst(0)))
      qcut(0) = 0
      ncut = 0
      do  i  = 1, nqlst(0)
        ncut = ncut+1
        qcut(ncut) = qxlen(nqlst(i))
      enddo

C ... Read spectral function and add components into sag
      allocate(agq(nomg,nfblst),sag(nomg,nqp))
      call dpzero(sag,nomg*nqp)
C     print *,iq; call prmx('agq',agq,nomg,nomg,nfblst)
      nq = 0
      do  iq = 1, nqp
        call dfdump(agq,nomg*nfblst,ifi)
        if (iq > 1) then
          if (abs(qxlen(iq)-qxlen(iq-1)) < tolq) cycle
        endif
        nq = nq+1
        nincl = 0
        do  kb = 1, nfblst
          ib = ifblst(kb)
          if (nblst /= 0) then
            i = iinear(nblst,ib,iblst,1)
            if (iblst(i) /= ib) ib = 0
C            if (iinear(nblst,ifblst(kb),iblst,1)) then
C              print *, kb
C            endif
          endif
          if (ib == 0) cycle
          nincl = nincl+1
          call daxpy(nomg,1d0,agq(1,kb),1,sag(1,nq),1)
        enddo
      enddo

C ... Pick top of color scale
      amax = maxval(sag)*ascl
      if (atop == 0) then
        atop = amax/3
      endif


   40 continue

      if (lfs) then
        if (nomg > 1) call rx('fs mode should have nomg=1')
        allocate(agq(nomg,nfblst),sag(nomg,nqp))
        do  iq = 1, nqp
          call dfdump(agq,nomg*nfblst,ifi)
          sag(1,iq) = dsum(nfblst,agq,1)*Ascl
        enddo
        do  iq = 1, nqp
          if (sag(1,iq) < abot) sag(1,iq) = 0
        enddo

        i = fopnx(spfile,94,0,jfi)
        if (i /= 0) call rx('problem with logical unit')
        call ywrm(0,' ',1,jfi,'(9f16.10)',sag,1,nqlst(1),nqlst(1),nqlst(2))

C        iq = 0
C        do  i = 1, nqlst(1)
C          do  j = 1, nqlst(2)
C            iq = iq+1
C            write(jfi,902) dble(i),dble(j),sag(1,iq)*Ascl,sag(1,iq)*Ascl
C          enddo
C          write(jfi,'(1x)')
C
C        enddo

        call rx0('done')

C ... Setup for gnuplot mode : make gnu.plt
      elseif (ifmt == 1) then

C       Get name of spectral function file
        i = fopnx(spfile,94,0,jfi)
        if (i /= 0) call rx('problem with logical unit')

C   ... Save in gnu script the command used to generate it
        call get_command(labstr)
        write (kfi,"('# ',a)") trim(labstr)

        do  i  = 0, nqlst(0)
          call awrit2(
     .    'set arrow from  graph 0,first %,5;5d, graph 0 to graph 1,first %,5;5d,0 nohead front lt -1',
     .    ' ',120,kfi,qcut(i),qcut(i))
        enddo

        if (nqlst(0) == nlab - 1) then
          labstr = 'set ytics border ( '
          do i = 0, nlab-1
            write(itos,'(F0.5)') qcut(i)
            labf = lab(i+1)
            if (trim(labf) == 'G') labf = '{/Symbol G}'

            labstr = trim(labstr)//'"'//trim(labf)//'" '//trim(itos)
            itos = ''

            if (i /= nlab-1) then
              labstr = trim(labstr)//','
            else
              labstr = trim(labstr)//')'
            endif
          enddo
        endif

        if (nlab > 0 .and. nqlst(0) == nlab - 1) then
          call awrit0(trim(labstr),' ',120,kfi)
        else
          call awrit0('unset ytics',' ',120,kfi)
        endif

        write(kfi,100)
  100   format(
     .    'set view 180,90'/
     .    'set palette defined (0 "white", .1 "blue", 1 "red");'/
     .    '#set palette rgbformulae -10,-13,-26'/
     .    '#unset colorbox'/
     .    '#set colorbox vertical user origin .02, .1 size .04,.8'/
     .    'set yzeroax lt 2 lw 3'/
     .    '#set zr [0:]'//
     .    'set cbrange resto')
        call awrit1('set cbrange[0:%,5;5d]',' ',120,kfi,atop)
        if (ifig == 1) then
          call awrit0('set term pos enh col eps 24',' ',120,kfi)
          fignam = 'spf.ps'
        elseif (ifig == 2) then
          call awrit0('set term svg fsize 22',' ',120,kfi)
          fignam = 'spf.svg'
        elseif (ifig == 3) then
          call awrit0('set terminal png enhanced  nocrop  font "Vera" size 800,600',' ',120,kfi)
          fignam = 'spf.png'
        endif
        call awrit0('set out "'//trim(fignam)//'"',' ',120,kfi)
        if (window(1) == NULLI) then
          window(1) = ommin; window(2) = ommax
        endif
        if (escl /= 1) then
          call awrit3('sp [%;3g:%;3g][][] "'//trim(spfile)//'" u ($1*%;5d):2:3 notitle with pm3d at b'//spext,' ',160,kfi,
     &                window(1),window(2),1/escl)
        else
          call awrit2('sp [%;3g:%;3g][][] "'//trim(spfile)//'" u 1:2:3 notitle with pm3d at b'//spext,' ',160,kfi,
     &                window(1),window(2))
        endif
        if (ldrawqp) then
          prfmt = 'lc rgb "yellow"'  ! Default color
          if (csw(0,1) > 0) then
            prfmt = 'lc rgb "#      "'
            write(prfmt(10:11),'(z2)') nint(255*min(max(csw(1,1),0d0),1d0))
            write(prfmt(12:13),'(z2)') nint(255*min(max(csw(2,1),0d0),1d0))
            write(prfmt(14:15),'(z2)') nint(255*min(max(csw(3,1),0d0),1d0))
            if (prfmt(10:10) == ' ') prfmt(10:10) = '0'
            if (prfmt(12:12) == ' ') prfmt(10:10) = '0'
            if (prfmt(14:14) == ' ') prfmt(10:10) = '0'
          endif
C         Line thickness
          call awrit5('%a'//
     .      '%?#n>1#%-1j lw %i##',
     .      prfmt,len(prfmt),0,ltbold,2,3,4,5)

          do  bnd = 1, nband
            if (bnd == nband) spext = ''
            write(itos,'(I0)') bnd
            bndfilen = 'bnd'//trim(itos)//'.dat'
            if (escl /= 1) then
              call awrit2('"'//trim(bndfilen)//'" u ($2*%;5d):1:(%,5;5d) '//trim(prfmt)//' with lines notitle'//spext,' ',160,kfi,
     &                    1/escl,atop+1)
            else
              call awrit1('"'//trim(bndfilen)//'" u 2:1:(%,5;5d) '//trim(prfmt)//' with lines notitle'//spext,' ',160,kfi,atop+1)
            endif
            itos = ''
          enddo
        endif
        call awrit0('unset out',' ',120,kfi)

C       Additions to spf file
C       call awrit2('# %n;10,5D',' ',80,jfi,ncut+1,qcut)

C ... Setup for SpectralFunction.sh mode
      elseif (ifmt == 2) then
        call awrit2('# %n;10,5D',' ',80,jfi,ncut+1,qcut)
      else
        call rxi('format not recognized:',ifmt)
      endif

      smax = 0
      dw0 = (ommax-ommin)/(nomg-1)
      do  iw = 1, nomg
        omgi = ommin + dw0*(iw-1)
        nq = 0
        if (ifmt == 1) then
          write(jfi,"(1x)")
        endif
        do  iq = 1, nqp
          if (iq > 1) then
            if (abs(qxlen(iq)-qxlen(iq-1)) < tolq) cycle
          endif
          nq = nq+1
          smax = max(smax,sag(iw,nq)*Ascl)
          write(jfi,902) omgi,qxlen(iq),sag(iw,nq)*Ascl,sag(iw,nq)*Ascl
  902     format(2f10.5,50f14.6)
        enddo
      enddo

      if (ldrawqp) then
        do bnd = 1, nband
          write(itos,'(I0)') bnd
          bndfilen = 'bnd'//trim(itos)//'.dat'
          bndfile = fopng(bndfilen,-1,0)
          do iq = 1, nqp
            write(bndfile,902) qxlen(iq),(eigiq(bnd,iq)-efqp)*Ascl
          enddo
          itos = ''
          call fclr(bndfilen,bndfile)
        enddo
      endif

      call info8(1,0,0,' %i bands  window=%s,(%2;3d)  nw=%i  nq=%i of %i  Atop=%;3g  Amax*Ascl=%;3g  escl=%;3g',
     .  nincl,window,nomg,nq,nqp,atop,smax,escl)
      if (ldrawqp) call info2(1,0,0,' Drawing QP bands with ef shifted by %d',efqp,2)
! these are not error states so they should not return error codes (which trigger the error code backtrace system).
      if (ifmt == 1) then
        call info0(1,0,0,' Make figure with the command shown.  gnuplot creates file  '//trim(fignam))
        call info0(1,0,0,' gnuplot gnu.plt')
        call cexit(0,1)
      elseif (ifmt == 2) then
        call info0(1,0,0,' Make figure with the command shown.')
        call info2(1,0,0,' echo fe | SpectralFunction.sh -e -f 1 -w %;3g:%;3g',ommin/escl,ommax/escl)
        call cexit(0,1)
      endif
      end
