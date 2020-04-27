      subroutine fmain
      implicit none
      integer iarg,nbas,nvar,i,j,ifi,garg,i1mach,iprint,chrln,fopna,
     .        fopnx,ifip,ifif,ifih,ifii,ifiw,lenp,lenf,lenh,leni,lenw,
     .        lenprg,pdim,nkill,it,ir,ierr,rdm,isw,idamax,nd,ns,icom
      logical T,F,cmdstr,lsequ,a2bin,lbroy,lfp,lsd,lhess

      parameter (nd=4,ns=6)
      parameter (chrln=2048)
      character*(chrln) prgnam,outs,first,fnamp,fnamf,fnamh,fnami,fnamw
      double precision xtol,gtol,step,gfac,wkg(28),gnorm,gmax
      double precision dsqrt,ddot
      double precision , allocatable::pos(:,:),force(:,:),p(:,:),
     .                               hess(:,:)
      data T, F / .True. , .False. /
      integer wksize
      parameter (wksize=10000)
      integer w(wksize)
      common /w/ w
      call pshpr(0)
      call wkinit(wksize)
      call poppr

      call finits(1,0,0,i)
      call pshpr(61)
      prgnam = 'm_gradzr'
      call skpblb(prgnam,chrln,lenprg)
C     Defaults ..
      lfp = T
      lbroy = F
      lsd = F
      lhess = T
      xtol = 1d-3
      gtol = 1d-3
      step = 0.01
      gfac = 1.2d0
      fnamp = ''
      fnamf = ''
      fnamh = ''
      fnami = ''
      fnamw = ''
      lenp = 0
      lenf = 0
      lenh = 0
      leni = 0
      lenw = 0
      nkill = 0
      nbas = 0
      it = 0
      icom = 0
      call dcopy(28,0d0,0,wkg,1)

C     no arguments: print usage
      if (.not. cmdstr(1,first)) goto 3
      iarg = 1
    1 continue
C     finished all switches
      if (.not. cmdstr(iarg,first)) goto 2
      if (lsequ(first,'-pr',3,' ',i)) then
        i = 0
        if (.not. a2bin(first(4:40),j,2,0,' ',i,-1)) goto 3
        call pshpr(j)
      else if (first == '-broy') then
        lfp = F
        lbroy = T
      else if (first == '-sd') then
        lfp = F
        lsd = T
      else if (first == '-fp') then
      else if (garg('-xtol=',iarg,4,' ',1,1,i,j,xtol) /= 0) then
      else if (garg('-gtol=',iarg,4,' ',1,1,i,j,gtol) /= 0) then
      else if (garg('-step=',iarg,4,' ',1,1,i,j,step) /= 0) then
      else if (garg('-gfac=',iarg,4,' ',1,1,i,j,gfac) /= 0) then
      else if (garg('-nkill=',iarg,2,' ',1,1,i,j,nkill) /= 0) then
      else if (garg('-nbas=',iarg,2,' ',1,1,i,j,nbas) /= 0) then
      else if (first(1:7) == '-pfile=') then
        fnamp = first(8:)
        call skpblb(fnamp,chrln,lenp)
      else if (first(1:7) == '-ffile=') then
        fnamf = first(8:)
        call skpblb(fnamf,chrln,lenf)
      else if (first(1:7) == '-hfile=') then
        fnamh = first(8:)
        call skpblb(fnamh,chrln,lenh)
      else if (first(1:7) == '-ifile=') then
        fnami = first(8:)
        call skpblb(fnami,chrln,leni)
      else if (first(1:7) == '-wfile=') then
        fnamw = first(8:)
        call skpblb(fnamw,chrln,lenw)
      endif
      iarg = iarg+1
      goto 1
    2 continue
      if (nbas == 0) then
        print *, prgnam(1:lenprg+1), ' requires nbas'
        goto 3
      endif
      if (fnamp == '' .or. fnamf == '' .or.
     .    fnami == '' .or. fnamw == '') then
        print *, prgnam(1:lenprg+1), ' requires filenames'
        goto 3
      endif
      if (fnamh == '') lhess = F
      call awrit6(prgnam(1:lenprg+1)//
     .           ' starting .. nbas=%i, xtol=%,3g, '//
     .           'gtol=%,3g, step=%,3g, gfac=%,2d files: '
     .            //fnamp(1:lenp+1)//' '
     .            //fnamf(1:lenf+1)//' '//fnami(1:leni+1)//
     .            ' '//fnamh(1:lenh+1)//' '//fnamw(1:lenw+1)//
     .            ', nkill=%i',' ',256,
     .            i1mach(2),nbas,xtol,gtol,step,gfac,nkill)
      if (fopnx(fnamp,170,-1,-1) /= 0) then
        ifip = fopna(fnamp,-1,1)
      else
        call awrit0(' file '//fnamp(1:lenp+1)//' not found',' ',
     .              128,i1mach(2))
        goto 3
      endif
      if (fopnx(fnamf,170,-1,-1) /= 0) then
      else
        call awrit0(' file '//fnamf(1:lenf+1)//' not found',' ',
     .              128,i1mach(2))
        goto 3
      endif
      if (lhess) then
        ifih = fopna(fnamh,-1,4)
      endif
      ifiw = fopna(fnamw,-1,4)
C allocate space for pos, force, hessian and workspace
      nvar = 3 * nbas
      allocate (pos(3,nbas),stat=ierr)
      if (ierr /= 0)
     .      call rx0(prgnam(1:lenprg+1)//': failed to allocate pos')
      allocate (force(3,nbas),stat=ierr)
      if (ierr /= 0)
     .      call rx0(prgnam(1:lenprg+1)//': failed to allocate force')
      allocate (hess(nvar,nvar),stat=ierr)
      if (ierr /= 0)
     .      call rx0(prgnam(1:lenprg+1)//': failed to allocate hess')
      if (lsd) pdim = 7
      if (lfp) pdim = 7+nvar
      if (lbroy) pdim = 12+2*nvar
      allocate (p(nvar,pdim),stat=ierr)
      if (ierr /= 0)
     .      call rx0(prgnam(1:lenprg+1)//': failed to allocate p')

C open files and read data
      if (fopnx(fnami,170,-1,-1) /= 0) then
        ifii = fopna(fnami,-1,1)
        if (iprint() >= 20) then
          print *,' Found iterations file'
        endif
      else
        it = 1
        ir = 0
        ifii = fopna(fnami,-1,0)
        if (iprint() >= 20) then
          print *,' Iterations file not found, start new minimisation'
        endif
      endif
      if (it == 0) then
        call ioiter(0,ifii,it,ir,icom,wkg,gnorm,gmax,ierr)
        if (iprint() >= 20) then
          call awrit4(' Read iterations file, ierr = %i, it = %i,'//
     .                ' ir = %i, icom = %i',' ',128,i1mach(2),
     .                ierr,it,ir,icom)
        endif
      endif
      if (ierr == 1) then
        it = 1
        if (iprint() >= 20) then
          print *, 'Iteration file malformed or empty. Setting it = 1'
        endif
      endif
      if (it > 1) then
        call iopwk(0,ifiw,p,nvar*pdim,ierr)
        if (ierr == 1) then
          print *, 'Workfile empty or corrupted. Exit'
          call cexit(0,1)
        endif
      endif

C read coordinate and force files
      call pshpr(0)
      call iopos(F,0,fnamp,nbas,pos)
      call iopos(F,0,fnamf,nbas,force)
      call poppr

C make workarray for gradzr
      call dcopy(nvar,pos,1,p,1)
      call dcopy(nvar,force,1,p(1,2),1)
      call dscal(nvar,-1d0,p(1,2),1)

C set or read Hessian matrix
      if (it == 1 .or. .not. lhess) then
        call dcopy(nvar*nvar,0d0,0,hess,1)
        call dcopy(nvar,1d0,0,hess,nvar+1)
      else
        if (rdm(ifih,2,nvar*nvar,' ',hess,nvar,nvar) < 0) then
C ... attempt to read Hessian failed.  Make it unity
          call dcopy(nvar*nvar,0d0,0,hess,1)
          call dcopy(nvar,1d0,0,hess,nvar+1)
          if (iprint() >= 20) print *, 'Failed to read hessian'
        endif
        call fclose(ifih)
      endif

C set relaxation method
      if (lsd)   isw = 00021 + 40
      if (lfp)   isw = 00121 + 40
      if (lbroy) isw = 00221 + 00
      call gradzr(nvar,p,hess,xtol,step,xtol,gtol,gfac,wkg,isw,ir)
      gnorm = dsqrt(ddot(nvar,p(1,2),1,p(1,2),1))
      gmax = p(idamax(nvar,p(1,2),1),2)
      if (ir > 0) then
        call rx1(prgnam(1:lenprg+1)//'
     .           : gradzr returned ir=%i ... aborting',ir)
      endif
      if ((lsd .or. lfp) .and. iprint() >= 20) then
        if (iprint() >= 20) call awrit4(' RELAX line %d:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n==-1#new line minimization##'//
     .    '%-1j%?#n==-2#bracketed root this line##'//
     .    '%-1j%?#n==-3#bracketed root this line##'//
     .    '%-1j%?#n==-4#extrapolated along this line##'//
     .    '%-1j%?#n==-5#extrapolated along this line##'//
     .    '%-1j%?#n==-6#is in trouble##'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,i1mach(2),wkg(19),ir,
     .    wkg(1)*p(idamax(nvar,p(1,nd),1),nd),gnorm)
      elseif (lbroy) then
        if (iprint() >= 20) call awrit4(' RELAX:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n<0#completed Broyden iter %i#%j#'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,i1mach(2),
     .    ir,-ir,p(idamax(nvar,p(1,ns),1),ns),gnorm)
      endif
      if (ir == 0) icom = 1
      if (icom == 0 .and. wkg(28) < 0) icom = -1

C printout
      if (iprint() >= 40) then
        call info0(-40,0,0,'        Gradients:')
        print 100, (p(i,2), i = 1, nvar)
        if (.not. lsd .and. (it /= 1 .or. xtol /= 0)) then
          call info0(-40,0,0,'      Diagonal inverse Hessian:')
          print 100, (hess(i,i), i = 1,nvar)
        endif
      endif
  100 format(10f8.3)

C write Hessian to disc
      if (lhess) then
        rewind ifih
        write(ifih) nvar,nvar,11
        write(ifih) hess
        call fclose(ifih)
      endif

C write workarray to disc
      call iopwk(1,ifiw,p,nvar*pdim,ierr)

C update pos file (new positions in workarray p)
      call pshpr(0)
      call iopos(T,0,fnamp,nbas,p)
      call poppr

C update iterations file
      it = it + 1
      call ioiter(1,ifii,it,ir,icom,wkg,gnorm,gmax,ierr)

      call cexit(0,1)

    3 continue
      print 300, prgnam(1:lenprg+1)
  300 format (
     .  ' usage: ',a,' [switches] ' /
     .  '     -pr#            verbosity (default 30) ' /
     .  '     -nbas=#         number of atoms ',
     .                        '(all atoms relax at present)' /
     .  '     -pfile=filename positions file '/
     .  '     -ffile=filename forces file '/
     .  '     -hfile=filename Hessian file '/
     .  '     -ifile=filename iterations file '/
     .  '     -wfile=filename workspace file (for p array) '/
     .  '     -fp             use Fletcher-Powell minimisation',
     .                        ' (default) ' /
     .  '     -broy           use Broyden minimisation ' /
     .  '     -sd             use Steepest descent minimisation ' /
     .  '     -xtol=#         tolerance in displacements ',
     .                        '(default 1d-3)' /
     .  '     -gtol=#         tolerance in gradients (default 1d-3) ' /
     .  '     -step=#         maximum step in displacements ',
     .                        '(default 1d-2)' /
     .  '     -nkill=#        reset Hessian after # iterations ' /
     .  '     -gfac=#         Extrapolation growth factor for line ' ,
     .                       'minimisations (default 1.2) ' )
      call cexit(0,1)

      end

      subroutine ioiter(mode,ifii,it,ir,icom,wkg,gnorm,gmax,ierr)
C- I/O info file for m_gradzr
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode : 0 read, 1 write
Ci   ifii : file handle
Ci   it, ir, wkg
Co Outputs:
Co   it, ir, wkg
Co   ierr : 0 if read successful, 1 otherwise
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer mode,ifii,it,ir,icom,ierr
      double precision wkg(28),gnorm,gmax
      ierr = 0
      rewind ifii
      if (mode == 0) then
        read (ifii, 100, end=1, err=1) it, ir, icom, wkg, gnorm, gmax
      else
        write (ifii, 100) it, ir, icom, wkg, gnorm, gmax
      endif
      return
    1 continue
      ierr = 1
  100 format (3i7/ 30(g28.16/))
      end

      subroutine iopwk(mode,ifiw,p,pdim,ierr)
C- I/O work file for m_gradzr
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode : 0 read, 1 write
Ci   ifiw : file handle
Co Outputs:
Co   ierr : 0 if read successful, 1 otherwise
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer mode,ifiw,ierr,pdim
      double precision p(pdim)
      ierr = 0
      rewind ifiw
      if (mode == 0) then
        read (ifiw, end=1, err=1) p
      else
        write (ifiw) p
      endif
      return
    1 continue
      ierr = 1
      end
