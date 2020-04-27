      subroutine fmain
C - Band plotter using DOS file in standard format for lmto programs.
C ----------------------------------------------------------------------
Cl Local variables
Cl   ifr   If 1, there is a reference dos to be drawn, with color weights
Cs Command-line switches
Cs   -dif      : differentiate DOS
Cs   -ef=#     : shift Fermi level to #
Cs   -shiftef= : shift Fermi level by #
Cs   -escl=    : scale energy axis by #
Cs   -esclxy=  : scale E by #
Cs   -ext      : scale E by # and DOS by 1/#
Cs   -fplot    : makes data files and script for fplot command.
Cs   -int      : integrate DOS
Cs   -lst=     : specify groups of partial DOS to sum
Cs   -lst2[=]  : same as -lst, but for second spin
Cs   -ref      : name reference for plotting single DOS with color weights
Cu Updates
Cu   10 Sep 17 Added seq tag in -lst
Cu   25 May 17 Added switches to -fplot, removed -inc, -long, -ext
Cu   16 Mar 11 Added -inc and -esclxy switches
Cu   25 Apr 09 -ref option handles 2 color weights
Cu   07 Jun 06 Added -ref options
Cu   04 Nov 02 -lst2 w/out arguments automatically takes members
Cu                   of -lst and increments by 1 (for spin pol)
Cu   13 Jun 02 Added -ef= and -escl= command-line arguments
Cr Remarks
Cr   The -ref option is designed to draw a single DOS-like object, specified
Cr   by fn in switch
Cr     -ref:fn=...
Cr   In this case the last argument, another DOS-style file, serves to
Cr   supply the color weights file for the figure being drawn.
C ----------------------------------------------------------------------
      implicit none
      integer nemx,ndinmx
      integer i,ii,j,k,j1,j2,id,idtop,ie,ifmt,is,ist,istart,is1,is2,i1mach,
     .  ifi,ifr,ifinc,md,nd,ne,ne1,ne2,nnum,num1,num2,nnumsav(2),nnseq,npseq,nrnde,
     .  nspin,npdif,lerr,iouts,iiouts,iarg,lgcmd,ner,mdi,nchan,iscale,iseq
      integer, parameter:: ndmx=25, maxnum=300, NULLI=-99999
      integer num(maxnum),numsav(maxnum),nseq(maxnum),fontsz,ix(maxnum),ltfrme(3),iv(10)
      double precision cmhgt,dcm,dssum,dbot,dtop,dtopmx,dtoplo,dxtik,dy,dytik,
     .  dytk,e1,e2,eferm,emax,emin,estart,top,xmax,ymax,ypos1,ypos2,emun,dosscl,
     .  tmy,emux,efurm,dul,frmcs(0:3),fillcs(0:4),de,eshift,escl,eminr,emaxr,
     .  csw(0:15)
      double precision dmx(ndmx+1),xd(2),yd(2),xx,xxv(10),wgt(maxnum)
      character font*19,dc*1,dcfp*1
      integer xnum,ynum
      logical lrat,l2chan,prlong,panel1,esclxy
      character*100 outs(99)*80,xnfmt*20,ynfmt*20,extnam*40
      character xlab*120,ylab*120,title*120,fn*120,ltdos*80
      character*512 strn,doslst,dosls2,fpstrn
      character*1 strn2(0:1)
      procedure(logical) :: cmdopt,a2bin,rdstrn
      procedure(integer) :: cmdoptswx,isw,fopng,fopnx,mkilsd,a2vec,awrite,nargf
      equivalence (strn,strn2)
C#ifdef F90
      real(8),allocatable:: e(:),dtot(:),din(:,:),d(:,:,:),dref(:)
C#elseC
C      parameter(nemx=5001,ndinmx=2000)
C      double precision e(nemx),dtot(nemx),dref(nemx)
C      double precision din(nemx,ndinmx),d(nemx,2,ndmx)
C#endif

C Heap allocation
      integer wksize
      parameter(wksize=100000)
C#ifndefC F90
C      common /static/ d,din,e,dtot
C#endif

      call pshpr(0)
      call wkinit(wksize)
      call poppr

C --- Setup and defaults ---
C     check for obsolete switches
      if (cmdopt('-long',5,0,strn)) call rx(
     .  'switch -long no longer supported. Use -fplot~long')
      if (cmdopt('-inc',4,0,strn)) call rx(
     .  'switch -inc no longer supported. Use -fplot~inc=fn')
      if (cmdopt('-ext',4,0,strn)) call rx(
     .  'switch -ext no longer supported. Use -fplot~ext')

      iseq = 0
      call dpzero(csw,16)
      l2chan = .false.
      xnfmt = '%;4d'
      ynfmt = '%;4d'
      frmcs(0:3) = [3d0,0d0,0d0,0d0]
      fillcs(0:4) = [103d0,1d0,1d0,1d0,1d0]
      iouts = 1
      forall (i = 1:99) outs(i) = ' '
      call finits(2,0,0,i)
      dtopmx = 9999d0
C ... multiple of dos to round to
      nrnde = 5
C ... trial tick spacings
      dytik = 10d0
      dxtik = .5d0
C ... page size in cm
      xmax = 18d0
      ymax = 25d0
C ... character fonts
      font = '/Times-Roman'
      fontsz = 24
C ... Frame line type
      ltfrme(1) = 3
      ltfrme(2) = 0
      ltfrme(3) = 0
C ... If switch present, pldos generates an fplot command and appropriate files
      dcfp = ' '
      if (cmdopt('-fplot',6,0,fpstrn)) then
        dcfp = fpstrn(7:7)   ! separator for fplot switches
        extnam = 'dat'; j1 = cmdoptswx('-fplot','ext=','')+4
        if (j1-4 > 0) then
          call nwordg(fpstrn,0,dcfp,1,j1,j2)
          extnam = fpstrn(j1:j2)
        endif
        ifi = fopng('dosp.'//extnam,10,0)
      endif
      prlong = cmdoptswx('-fplot','long','') > 0

C ---- Input from file ifi -----
C     if (fxst('dos') < 1) goto 199
C     ifi = fopna('dos',-1,0)
      do  iarg = 1, nargf()-1
        call getarf(iarg,strn)
        if (strn(1:1) /= '-') then
          ifi = fopng(strn,-1,1)
          goto 8
        endif
      enddo
      goto 199
    8 continue
      print *

C --- Case using reference dos.  Integer flags this option ---
      ifr = 0
      iscale = 0
      if (cmdopt('-ref',4,0,strn)) then
        iscale = 0
        nchan = 0
        dc = strn(5:5)
        if (dc == ' ') then
          call info0(0,0,0,' pldos needs file name to accompany '//strn//'%a')
          goto 199
        endif
C   ... Return here to resume parsing for arguments
        j2 = 5-1
   10   continue
        j2 = j2+1
        if (strn(j2:j2) == dc) goto 10
        j1 = min(len(strn),j2)
        call nwordg(strn,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (strn(j1:j1+2) == 'fn=')  then
            if (j1+3 > j2) call rx('pldos: bad file name')
            fn = strn(j1+3:j2)
            k = fopnx(strn(j1+3:j2),172,-1,-1)
            if (k == 0) then
              call rxs2('pldos: missing file "',strn(j1+3:j2),'"')
            endif
            ifr = fopng(strn(j1+3:j2),-1,1)
          elseif (strn(j1:j1+4) == 'chan=')  then
            if (j1+5 > j2) call rx('pldos: bad list, chan=..')
            call mkils0(strn(j1+5:j2),nchan,ix)
            if (k > 100) call rx('increase maxnum')
            call mkilst(strn(j1+5:j2),k,ix)
            if (k < 0) call rx('pldos: bad list, chan=..')
          elseif (strn(j1:j2) == 'scale')  then
            iscale = 1
          else
            call rxs('pldos: fail to recognize argument ',
     .        strn(j1:j2))
          endif
          goto 10
        endif
        if (ifr == 0) then
          call info0(0,0,0,
     .      ' pldos: -ref option requires a file be specified')
          goto 199
        endif
        if (nchan == 0) then
          nchan = 1
          ix(1) = 1
        endif
C#ifdef F90
        rewind ifr
        read(ifr,*) eminr,emaxr,ner,nd,nspin
        rewind ifr
        nemx = ner
        ndinmx = nd*nspin
        allocate(dref(nemx),din(nemx,ndinmx))
C#endif
C       Read reference dos, combine channels to make dref
        call dosio(din,nemx,ndinmx,ner,nd,eminr,emaxr,nspin,
     .    efurm,dul,ifmt,ifr)
C       Assemble reference DOS
        call dpzero(dref,ner)
        do  i = 1, nchan
          k = ix(i)
          if (k > nd*nspin) call rxi(
     .      'pldos: reference channel outside file''s range 1 ..',
     .      nd*nspin)
          call daxpy(ner,1d0,din(1,k),1,dref,1)
        enddo
C#ifdef F90
        deallocate(din)
C#endif
        call fclose(ifr)
        ifr = 1
      endif

C#ifdef F90
      rewind ifi
      read(ifi,*) emun,emux,ne,nd,nspin
      rewind ifi
      nemx = ne
      ndinmx = nd*nspin
      allocate(e(nemx),dtot(nemx),din(nemx,ndinmx),d(nemx,2,ndmx))
C#endif
C     Read DOS
      call dosio(din,nemx,ndinmx,ne,nd,emun,emux,nspin,efurm,dul,ifmt,ifi)
      call fclose(ifi)
      nd = nd*nspin

      if (ifr /= 0) then
        call info0(20,0,1,' pldos: '//
     .    ' verifying reference and dos files have equivalent mesh ...')
        call sanrg(.true.,ner,ne,ne,' ','ne in reference file')
        call fsanrg(eminr,emun,emun,0d-9,' ','emin in reference file',.true.)
        call fsanrg(emaxr,emux,emux,0d-9,' ','emax in reference file',.true.)
      endif

      ltdos = '1,bold=3,col=0,0,0'
      if (iscale /= 0) ltdos = '1,bold=3,col=1,0,0,colw=0,0,1'

      eshift = 0
      if (cmdopt('-ef=',4,0,strn)) then
        ie = 4
        if (.not. a2bin(strn,eshift,4,0,' ',ie,-1)) call rxs('fplot: failed to parse argument  ',strn)
        eshift = efurm-eshift
      elseif (cmdopt('-shiftef=',9,0,strn)) then
        ie = 9
        if (.not. a2bin(strn,eshift,4,0,' ',ie,-1)) call rxs('fplot: failed to parse argument  ',strn)
      endif
      emin = emun-eshift
      emax = emux-eshift
      eferm = efurm-eshift
      escl = 1
      dosscl = 1

      ie = 0
      if (cmdopt('-escl=',6,0,strn)) ie = 6
      if (cmdopt('-esclxy=',8,0,strn)) ie = 8
      dtoplo = 5d0
      esclxy = .false.
      if (ie == 6 .or. ie == 8) then
        k = ie
        if (.not. a2bin(strn,escl,4,0,' ',k,-1)) call rxs('fplot: failed to parse argument  ',strn)
        emin = emin*escl
        emax = emax*escl
        eferm = eferm*escl
        if (ie == 8) then
          esclxy = .true.
          dosscl = 1/escl
          dtoplo = dtoplo/escl
        endif
      endif

      call awrit5(' pldos:  %i dos read  ef=%d  emin=%d  emax=%d  ne=%i',' ',80,i1mach(2),nd,eferm,emin,emax,ne)
      if (ne > nemx) stop '**** DIMENSION NEMX TOO SMALL'
      if (nd > ndinmx) stop '**** DIMENSION NDINMX TOO SMALL'
      if (escl /= 1) then
        call info2(0,0,0,' scaling energy axis by %;3d%?#n#, dos axis by reciprocal## ',
     .    escl,isw(dosscl /= 1))
      endif
      forall (ie = 1:ne) e(ie) = emin + (ie-1d0)*(emax-emin)/(ne-1d0)

C --- Select the next linear combination of dos ---
C     ifr should be zero, or 1 if reference dos is 1st dos
      md = 0   ! panel index
      npseq = 0 ! number of panels obtained through seq=, needed for wordg below
      dssum = 0

C ... Reference case: dref becomes the first DOS
      if (ifr == 1) then
        call dcopy(ne,dref,1,d,1)
        call idmxmn(ne,dref,1,i,j)
        dmx(1) = dref(j)
        dssum = dssum+dmx(md+ifr)
      endif

   90 continue
        print *, ' '
        if (iseq /= 0) goto 91
        call awrit1(' dos %i:  Enter list of dos to be added:',strn,80,0,md+1)

        if (cmdopt('-lst=',5,0,strn)) then
          call wordg(strn,1,'; ',md-npseq+1,is1,is2)
          doslst = ' '
          is1 = max(is1,6)
          if (is2 >= is1) doslst = strn(is1:is2)
        else
          call skpblb(strn,len(strn),i)
          call cwrite(strn,0,i+1,0)
          read(*,'(a512)') doslst
        endif
        if (cmdopt('-lst2=',6,0,strn)) then
          call wordg(strn,1,'; ',md-npseq+1,is1,is2)
          dosls2 = ' '
          is1 = max(is1,7)
          if (is2 >= is1) dosls2 = strn(is1:is2)
        elseif (cmdopt('-lst2',5,0,strn)) then
          dosls2 = '*'
        else
          dosls2 = ' '
        endif

C       Get list of dos channels and number
   91   continue
        is = 0
        nnum = 0
        num2 = 0
        call skipbl(doslst,512,is)
        if (is >= 512) goto 32
        if (iseq == 0 .and. doslst(1:4) == 'seq=') then
          call icopy(nnumsav(1)+nnumsav(2),num,1,numsav,1)
          nnseq = mkilsd(doslst(5:),-1,nseq)
          if (nnseq > size(num)) call rx('pldos: increase maxnum')
          if (nnseq <= 0) call rx('pldos: bad sequence')
          nnseq = mkilsd(doslst(5:),nnseq,nseq)
          iseq = 1
        endif
        if (iseq /= 0) then
          nnum = nnumsav(1)+nnumsav(2)
          num2 = nnumsav(2)
          if (num2 > 0) l2chan = .true.
          num(1:nnum) = numsav(1:nnum) + nseq(iseq) - numsav(1)
          iseq = mod(iseq+1,nnseq+1)
          if (iseq > 0) npseq = npseq + 1 ! Used in finding the next string, see wordg above
          goto 33  ! lists already made ... skip parsing for them
        endif
        nnum = mkilsd(doslst,-1,num)
        if (nnum > size(num)) call rx('pldos: increase maxnum')
        if (nnum <= 0) call rx('pldos: bad list')
        nnum = mkilsd(doslst,nnum,num)
   32   continue
        if (dosls2 == '*') then
          if (nnum == 0) goto 33
          num2 = nnum
          do  i = 1, nnum
            num(nnum+i) = num(i)+1
          enddo
          nnum = 2*nnum
          l2chan = .true.
          goto 33
        endif
        call skipbl(dosls2,512,is)
        if (is >= 512) goto 33
        num2 = mkilsd(dosls2,-1,num(1+nnum))
        if (num2+nnum > size(num)) call rx('pldos: increase maxnum')
        if (num2 <= 0) call rx('pldos: bad list')
        num2 = mkilsd(dosls2,num2,num(1+nnum))
        nnum = nnum + num2
        if (num2 > 0) l2chan = .true.
   33   continue
        if (num2 /= 0 .and. ifr /= 0) call rx('pldos: -ref not set up for 2nd channel')
        num1 = nnum - num2
        nnumsav(1) = num1; nnumsav(2) = num2
        if (cmdopt('-fplot',6,0,strn) .and. num2 > 0) then
          ifi = fopng('dosp2.'//extnam,11,0)
        endif
        if (nnum == 0) goto 98
        do  i = 1, nnum
          if (num(i) > nd)
     .      call fexit2(-1,1,' pldos (abort) dos %i exceeds number read (%i)',num(i),nd)
          wgt(i) = 1d0
        enddo
        if (cmdopt('-lst=',5,0,strn)) then
          call awrit3(' combining %i dos: %n:1i',' ',80,i1mach(2),nnum,nnum,num)
          if (minval(num(1:nnum)) < 0) call rx('pldos: bad list')
          strn = ' '
        else
          i = awrite(' Weights for dos:%n:1i?',strn,80,0,nnum,num,0,0,0,0,0,0)
          call cwrite(strn,0,i+1,0)
          read(*,'(a512)') strn
        endif
        if (strn == ' ' .or. strn == '/') goto 42
        is = 0
        ist = a2vec(strn,len(strn),is,4,', ',2,3,nnum,ix,wgt)
        if (ist <= 0) then
          print *, 'bad weights, try again'
          goto 90
        endif
        if (ist < nnum) call awrit2('w=%n:1d',' ',80,i1mach(2),nnum,wgt)
   42   continue
C       call awrit1('%Nsum w = %d',' ',80,i1mach(2),dsum(nnum,wgt,1))
        md = md+1
        if (md > ndmx) call fexit(-1,1,' Exit -1:  no. dos exceeds ndmx=%i',ndmx)
C   ... Accumulate the linear combination of dos
        call dpzero(d(1,1,md+ifr),ne)
        call dpzero(d(1,2,md+ifr),ne)
        do  44  i = 1, num1
   44   call daxpy(ne,dosscl*wgt(i),din(1,num(i)),1,d(1,1,md+ifr),1)
        do  45  i = num1+1, nnum
   45   call daxpy(ne,-dosscl*wgt(i),din(1,num(i)),1,d(1,2,md+ifr),1)

        if (cmdopt('-dif',4,0,strn)) then
          if (num2 /= 0) call rx('pldos: -dif not set up for num2')
          npdif = 2
          lrat = .false.
          if (strn(5:5) == ':') then
            is = 5
            ist = a2vec(strn,len(strn),is,4,', ',2,2,2,ix,xxv)
            if (ist < 0) goto 199
            npdif = nint(xxv(1))
            if (ist > 1 .and. xxv(2) > 0) lrat = .true.
          endif
          call dcopy(ne,d(1,1,md+ifr),1,dtot,1)
          call poldvm(e,dtot,ne,npdif,lrat,1d-6,lerr,d(1,1,md+ifr))
          if (lerr /= 0) call rx('pldos: error in rational poly interp')
C   ... Integration by the trapezoidal rule
        elseif (cmdopt('-int',4,0,strn)) then
          if (num2 /= 0) call rx('pldos: -int not set up for num2')
          if (ifr /= 0) call rx('pldos: -int not set up for ifr')
          call dcopy(ne,d(1,1,md+ifr),1,dtot,1)
          d(1,1,md+ifr) = 0d0
          de = e(2) - e(1)
          do  47  ie = 2, ne
            d(ie,1,md+ifr) = d(ie-1,1,md+ifr) +
     .                       (dtot(ie)+dtot(ie-1))*de/2
   47     continue
        endif
        dmx(md+ifr) = 0d0
        do  46  ie = 1, ne
   46   dmx(md+ifr) =
     .      dmax1(dmx(md+ifr),abs(d(ie,1,md+ifr)),abs(d(ie,2,md+ifr)))
        dssum = dssum+dmx(md+ifr)
      goto 90

C --- Find the first nonzero dos ---
   98 continue
      mdi = md+ifr
      if (mdi == 0) call rx('pldos: no dos to plot')
      do  ie = 1, ne
        istart = ie
        top = 0d0
        do  id = 1, mdi
          top = dmax1(top,d(ie,1,id))
        enddo
        if (top > 1d-5) exit
      enddo
   82 estart = emin + (istart-1)*(emax-emin)/(ne-1d0)

C --- Do the plotting ---
      print *, ' '
      call awrit5(' %i plots:  emin=%d  estart=%d  emax=%d  ef=%d',
     .  strn,72,i1mach(2),mdi,emin,estart,emax,eferm)
      call awrit3(' dmx=%n:1,2;2d  sum=%,2;2d',' ',80,i1mach(2),mdi,dmx,dssum)
      dtopmx = 0d0
      do  id = 1, mdi
        dtopmx = dmax1(dtopmx,dmx(id))
      enddo
      cmhgt = min(ymax/dssum*dtopmx,dtoplo)
      e1 = emin
      e2 = emax
C     print 335, 'Enter dmax, height of dmax (cm), e1, e2: '
      write(*,335,advance='no') 'Enter dmax, height of dmax (cm), e1, e2: '
  335 format(1x,a)
C     call cwrite(' Enter dmax, height of dmax (cm), e1, e2: ',0,41,0)
      read (5,*)  dtopmx,cmhgt,e1,e2
      ne1 = max0(nint((e1-emin)*ne/(emax-emin)+1),1)
      ne2 = min0(nint((e2-emin)*ne/(emax-emin)),ne)
      do  id = 1, mdi
        dmx(id) = 0d0
        do  ie = ne1, ne2
          dmx(id) = dmax1(dmx(id),abs(d(ie,1,id)),abs(d(ie,2,id)))
        enddo
      enddo
      if (e1 /= emin .or. e2 /= emax) call awrit4(
     .  'dmx[%;3d,%;3d] =%n:1,2;2d',' ',80,i1mach(2),e1,e2,mdi,dmx)
      print *, ' '

      ypos1 = 0.97d0

C ... Write files for fplot
      if (cmdopt('-fplot',6,0,strn)) then

C       Scale columns 2...md+ifr
        if (ifr == 1 .and. iscale /= 0) then
          do  ie = ne1, ne2
            if (d(ie,1,1) /= 0) then
              do  id = 2, mdi
                d(ie,1,id) = d(ie,1,id)/d(ie,1,1)
              enddo
            else
              do  id = 2, mdi
                d(ie,1,id) = 1
              enddo
            endif
          enddo
        endif

        ifi = fopng('dosp.'//extnam,-1,0)
        ypos1 = 1.4d0
        call awrit2('%% rows %i cols %i',strn,80,ifi,ne2-ne1+1,mdi+1)
        do  ie = ne1, ne2
          if (prlong) then
            write(ifi,334) e(ie), (d(ie,1,id), id=1, mdi)
  334       format(f9.5,6f15.6:/(9x,6f15.6))
          else
            write(ifi,333) e(ie), (d(ie,1,id), id=1, mdi)
  333       format(f9.5,6f11.4:/(9x,6f11.4))
          endif
        enddo
        call fclose(ifi)
        if (l2chan) then
          ifi = fopng('dosp2.'//extnam,-1,0)
          call awrit2('%% rows %i cols %i',strn,80,ifi,ne2-ne1+1,mdi+1)
          do  ie = ne1, ne2
            if (prlong) then
              write(ifi,334) e(ie), (d(ie,2,id), id=1, mdi)
            else
              write(ifi,333) e(ie), (d(ie,2,id), id=1, mdi)
            endif
          enddo
          call fclose(ifi)
        endif
      else
        call pltini('ps',0,.false.)
        call setfnt(font,fontsz)
      endif

      ltdos = '1,bold=3,col=0,0,0'
      if (iscale /= 0) then
        ltdos = '1,bold=3,col=0,0,0,colw=1,0,0'
        if (md == 2) ltdos='1,bold=3,col=0,0,0,colw=1,0,0,colw2=0,1,0'
      endif

C ... Header for plot.dos
      outs(1) = ' '
      if (cmdopt('-fplot',6,0,strn)) then
        print *
        lgcmd = fopng('plot.dos',-1,0); rewind lgcmd
        j1 = cmdoptswx('-fplot','inc=','')+4 ! Look for include file
        if (j1-4 > 0) then
          call nwordg(fpstrn,0,dcfp,1,j1,j2)
          k = fopnx(trim(fpstrn(j1:j2)),172,-1,-1)
          if (k == 0) call rxs2('pldos: missing file "',trim(fpstrn(j1:j2)),'"')
          print *, 'prepending contents of file ', trim(fpstrn(j1:j2))
          ifinc = fopng(trim(fpstrn(j1:j2)),-1,1)
          do  while (rdstrn(ifinc,strn,len(strn),.false.))
            write(lgcmd,'(a)') trim(strn)
          enddo
        endif
        if (ltdos /= ' ')
     .    write(lgcmd,'(''%char0 ltdos="'',a,''"''/)') trim(ltdos)
        call awrit0('fplot',' ',-80,lgcmd)
      endif
      if (iscale /= 0) mdi = 1

C ... For each panel, do
      do  id = 1, mdi
C     Round off dtop to a multiple of nrnde
      idtop = int(dmx(id)*1.05d0/nrnde + 1d0)
      dtop = idtop*nrnde
      if (esclxy) dtop = dmx(id)
C     Handle -fplot~dmin=#
      ii = cmdoptswx('-fplot','dmin=','') + 4
      if (ii > 4) then
        if (a2vec(fpstrn,len_trim(fpstrn),ii,4,dcfp//' ',2,3,1,iv,dtoplo) < 0)
     .    call rxs('failed to parse string ',fpstrn)
        dtop = dmin1(dmax1(dmx(id),dtoplo),dtopmx)
      else
        dtop = dmin1(dmax1(dtop,dtoplo),dtopmx)
      endif

      dbot = 0
      if (l2chan) dbot = -dtop
      dcm = dtop*cmhgt/dtopmx
      dy = dcm/ymax
      ypos2 = ypos1
      ypos1 = ypos2-dy
      if (ypos1 < 0)
     .  print *, 'pldos (warning): frame past bottom of page ...'
      call awrit6(' plot %i  dmx=%1,3;3d  dtop=%1,3;3d  dy=%1,3;3d'//
     .  '  y=%1,3;3d %1,3;3d',' ',80,i1mach(2),id,dmx(id),dtop,dy,ypos1,ypos2)
      if (cmdopt('-fplot',6,0,strn)) then
C       First frame: set up variables
        panel1 = .false.
        if (iouts == 1) then
          panel1 = .true.
          iouts = iouts+1
          call awrit5('%% var ymax=%1;3d dy=%1;3d dw=.00 ymax+=dy '//
     .      'emin=%1;3d emax=%1;3d ef=%d',
     .      outs(iouts),80,0,ypos2,ypos2-ypos1,e1,e2,eferm)
          iouts = iouts+1
        endif
        iouts = iouts+1 ! make a blank line
        call awrit3('%% var ymax-=dy+dw dy=%1;3d dmin=%1;3d dmax=%1;3d',
     .    outs(iouts),80,0,ypos2-ypos1,dbot,dtop)
        iouts = iouts+1
        call awrit0('%a -frme 0,1,{ymax-dy},{ymax}',outs(iouts),80,0)
C       Handle -fplot~open
        if (cmdoptswx('-fplot','open','') > 0 .and. panel1) call awrit0('%a -frmt th=3,0,0',outs(iouts),80,0)
        if (cmdoptswx('-fplot','open','') > 0 .and. .not. panel1) call awrit0('%a -frmt th=3,1,0',outs(iouts),80,0)

C       Handle -fplot~xl=strn
        if (id == mdi .and. cmdoptswx('-fplot','xl=','') > 0) then
          j1 = cmdoptswx('-fplot','xl=','') + 3
          call nwordg(fpstrn,0,dcfp,1,j1,j2)
          call awrit0('%a -xl '//trim(fpstrn(j1:j2)),outs(iouts),80,0)
        endif

C       Handle -fplot~tmy=#
        ii = cmdoptswx('-fplot','tmy=','') + 3
        if (dtop/dtopmx*cmhgt < 2) then  ! No tic marks if window too small
          call awrit0('%a -tmy 100',outs(iouts),80,0)
        elseif (ii-3 > 0) then  ! User specified tic marks
          if (a2vec(strn,len_trim(strn),ii,4,dcfp//' ',2,3,1,iv,tmy) < 0)
     .      call rxs('failed to parse string ',strn)
          call awrit1('%a -tmy %d',outs(iouts),80,0,tmy)
        endif
        iouts = iouts+1
        call awrit1('%a -p0 -x {emin},{emax}'
     .    //' -y {dmin},{dmax} -1p%?#n# -noxn##',outs(iouts),80,0,id-md)
        iouts = iiouts(outs,iouts,72)
        if (iscale /= 0) then
          call awrit1('%a -colsy %i -colsw 3'//
     .      ' -lt {ltdos} -ord y -qr dosp.'
     .      //trim(extnam),outs(iouts),250,0,id+1)
        else
          call awrit1('%a -colsy %i -lt {ltdos} -ord y -qr dosp.'//
     .      trim(extnam),outs(iouts),250,0,id+1)
        endif
        iouts = iiouts(outs,iouts,72)
        if (iouts > 99-1) call rx('pldos output string too long')
        if (l2chan) then
        call awrit1('%a -colsy %i -lt {ltdos} -ord y -qr dosp2.'//
     .      trim(extnam),outs(iouts),250,0,id+1)
        iouts = iiouts(outs,iouts,72)
        if (iouts > 99-1) call rx('pldos output string too long')
        endif
        if (eferm < e2 .and. eferm > e1) then
          call awrit0('%a -lt 2,bold=3,col=0,0,1,2,.5,.05,.5 '//
     .      '-tp 2~{ef},{dmin},{ef},{dmax}',outs(iouts),250,0)
          iouts = iiouts(outs,iouts,72)
        endif
      else
        call pltstp(5d-2,0.98d0,ypos1,ypos2)
        call pltstu(e1,e2,0d0,dtop,0d0)
        dytk = dytik
        if (dtop < 10.01d0) dytk = 1d0
        if (dtop > 200d0) dytk = 50d0
        xnum = 0
        ynum = 0
        if (id == mdi) then
          xnum = 1
          ynum = 1
        endif
        xlab = ' '
        ylab = ' '
        title = ' '
        xx = NULLI
        call frme(-1d30,dxtik,1,.025d0,.6d0,0,xnum,xnfmt,xlab,xx,
     .            -1d30,dytk ,1,.025d0,.6d0,0,ynum,ynfmt,ylab,xx,
     .             ' ',xxv,xxv,xxv,title,ltfrme,'  ',frmcs,fillcs)
        call plntyp(1,3,2d0,0.5d0,0d0,0d0)
        call plcrv(e(ne1),d(ne1,1,id),e,e,e,e,ne2-ne1+1,e1,e2,0d0,dtop,
     .    0,-2d0,csw,11)
        if (eferm < e2 .and. eferm > e1) then
          xd(1) = eferm
          xd(2) = eferm
          yd(1) = 0
          yd(2) = dtop
          call plntyp(2,3,2d0,1d0,0d0,0d0)
          call plcrv(xd,yd,xd,xd,xd,xd,2,xd(1)-0.1d0,xd(2)+0.1d0,yd(1),
     .      yd(2),0,-2d0,csw,11)
        endif
        print *, ' '
      endif
      if (iscale /= 0) goto 21
      enddo
   21 continue
      if (cmdopt('-fplot',6,0,strn)) then
C       call awrit1('%a',outs,80,-i1mach(2),0)
        lgcmd = fopng('plot.dos',-1,0)
        do  25  is = 1, iouts
          call skpblb(outs(is),80,i)
C         call cwrite(outs(is),0,i,1)
          write(lgcmd,'(a)') outs(is)(1:i+1)
   25   continue
        call fclose(lgcmd)
        if (cmdoptswx('-fplot','sh','') > 0) then
          call info0(2,1,0,' pldos: invoking system call: fplot -pr10 -f plot.dos')
          call fsystm('fplot -pr10 -f plot.dos',i)
        else
          print *, ' ... to make a postscript file, invoke:'
          print *, ' fplot -f plot.dos'
        endif
      else
        call pltdmp(0)
      endif
C#ifdef F90
      deallocate(e,dtot,din,d)
C#endif
      call cexit(0,1)

C   99 continue
C      call rx('pldos encountered error in attempting to open file')
  199 print *, 'usage:  pldos [switches] filename'
      print *, '        pldos reads from standard input 4 numbers:    dosmax, height of dosmax (cm), emin, emax'
      print *, ''
      print *, 'Switches:'
  399 format(/4x,a:t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a:/t28,a)
      print 399,'-fplot[~options]',
     .  'makes data files and script for fplot command.  Options are separated',
     .  'by the first character following -fplot (taken to be ~ below)',
     .  '~long    writes fplot data files with 6 decimals, rather than 4',
     .  '~open    fplot script adapted to give open feel to figure',
     .  '~dmin=#  Restricts the minimum dosmax in any panel to be at least #',
     .  '~tmy=#   Sets the tic mark spacing on the DOS axis',
     .  '~xl=strn assign an x axis label',
     .  '~inc=fn  prepends contents of file fn to fplot script file',
     .  '~sh      invoke fplot in a subshell to make fplot.ps',
     .  '~ext=nam writes file dosp.nam instead of dosp.dat'
      print 399,'-lst=list[;list;...]',
     .  'specify groups of partial DOS to sum',
     .  'list consists of numbers or ranges separated by comma,',
     .  'e.g. 11,2,5+1:8+2:2 combines channels 11 2 6 8 10',
     .  'while semicolon starts a new list in a new panel',
     .  '',
     .  'Special seq tag: following an existing panel, you can use seq=seq-list.',
     .  'For each element in this list the preceding list will be replicated,',
     .  'but uniformly shifted so that the first member in the DOS list is',
     .  'the current element in seq-list.',
     .  'Example : -lst="1:4;seq=5,9" creates three panels with',
     .  '  panel 1 (1:4)  panel 2 (5:8)  panel 3 (9:12)'
      print 399, '-lst2|-lst2=list...',
     .  'same as -lst,but for second spin',
     .  '-lst2 without args: use args of -lst, incrementing by 1'
      print 399, '-ref:fn=filename[:chan=ref-channel][:scale]','',
     .  'single DOS with color weights',
     .  'reference DOS is total DOS file, standard format',
     .  'filename is partial DOS file, standard format, e.g.',
     .  "pldos -ef=0 -escl=13.6 -fplot '-lst=3;1'"//
     .  ' -ref:fn=dos.cr:scale dosl.cr'
      print 399, '-dif:#[,rat]] | -int','differentiate or integrate DOS'
      print 399, '-ef=#|-shiftef=#','shift Fermi level to or by #'
      print 399, '-escl=# | -esclxy=#','scale energy axis by #',
     .  '-esclxy scales E by # and DOS by 1/#'
      if (cmdopt('-h',2,0,strn)) call cexit(0,1)
      call cexit(-1,1)
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
      subroutine dosio(dos,nemx,ndmx,ndos,nld,emin,emax,nspin,
     .   eferm,del,ifmt,ifile)
C- I/O for dos, MSM's format
      implicit none
C Passed parameters
      integer ifile,ifmt,ndmx,ndos,nemx,nld,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,1)
C Local parameters
      integer ie,ild,iedum

C --- Write branch ---
      if (ifile < 0) then
        write(-ifile,760) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ifmt == 0) then
          do  10  ie = 1, ndos
   10     write(-ifile,761) ie,(dos(ie,ild),ild=1,nld*nspin)
        else
          do  11  ild = 1, nld*nspin
   11     write(-ifile,762) (dos(ie,ild),ie=1,ndos)
        endif
  761   format(i5,6f12.5/(5x,6f12.5))
  760   format(2f10.5,3i5,2f10.5,i5)
  762   format(5f14.6)
      endif
C --- Read branch ---
      if (ifile > 0) then
        read(ifile,*) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ndos > nemx) call rxi('dosio: ndos > nemx =',nemx)
        if (nld*nspin > ndmx) call rxi('dosio: nld > ndmx=',ndmx)
        if (ifmt == 0) then
          do  20  ie = 1, ndos
   20     read(ifile,761) iedum,(dos(ie,ild),ild=1,nld*nspin)
        elseif (ifmt == 1) then
          do  21  ild = 1, nld*nspin
   21     read(ifile,762) (dos(ie,ild),ie=1,ndos)
        else
          call rx('dosio: bad fmt')
        endif
      endif
      end
