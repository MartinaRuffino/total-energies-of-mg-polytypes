      subroutine mcrelx(lrlx,it,nit,nbas,atid,step,xtol,gtol,f,bas,p,
     .                  w,icom)
C- Relax atomic positions
C ----------------------------------------------------------------------
Ci Inputs:
Ci   it:        iteration number
Ci   nit:       total number of iterations called for
Ci   nbas
Ci   w:         work array for the hessian
Ci   p:         work array for gradzr (dimensioned in rlxstp)
Ci   delta,nelts used for printout only, if U[L] turned on (TBE)
Co Outputs:
Co   bas        updated
Co   icom       1 if relaxation is complete
Co              0 otherwise
Cr Remarks:
Cr   Adapted from lm-6.10c relax.f
Cr   This is a rather stripped down version. For now, all atoms are
Cr   relaxed and the Hessian is read and used if found, by default.
Cr   NB iomc reads ifrlx, but I can't imagine a context where you'd
Cr   want to relax only a part of a molecule.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer lrlx,it,nit,nbas,icom
      character*8 atid(nbas)
      double precision step,xtol,gtol,f(3,nbas),w(3*nbas,3*nbas),
     .                 p(3*nbas,6),bas(3,nbas)
C Local Variables
      character strn*128,dumstr*6
      integer nvar,i,j,ipr,ifi,i1mach,ix,ir,iprint,isw,rdm,idamax
      logical lpp,cmdopt,a2bin,rdhess
      double precision xtoll,grfac,wkg(28)
      double precision ddot
      save ir,wkg

      rdhess = .true.
      call pshpr(iprint()+30)
      ipr = iprint()
      nvar = 3*nbas
      ifi = 74

C --- Make vector of positions and gradients ---
      call dcopy(nvar,bas,1,p,1)
      call dcopy(nvar,f,1,p(1,2),1)
      call dscal(nvar,-1d0,p(1,2),1)

C --- Initialization ---
      if (it == 1) then
        ir = 0
        call dcopy(nvar*nvar,0d0,0,w,1)
        call dcopy(nvar,1d0,0,w,nvar+1)
C   ... Read Hessian from disc
        if (rdhess) then
          call mopen(ifi,'hssn','u')
          rewind ifi
C     ... Hessian written in rdm-compatible format
          if (rdm(ifi,2,nvar*nvar,' ',w,nvar,nvar) < 0) then
C       ... Attempt to read Hessian failed.  Make it unity
            if (ipr > 30) print *,
     .        'RELAX: (warning) failed to read Hessian; set to unity'
            call dcopy(nvar*nvar,0d0,0,w,1)
            call dcopy(nvar,1d0,0,w,nvar+1)
          else
            if (ipr >= 30) print *,' RELAX: Hessian read from disc'
          endif
        endif
        call fclose(ifi)
      endif

C --- Relax ---
      grfac = 1.2d0
      if (cmdopt('-grfac=',6,0,strn)) then
        j = 6
        if (.not. a2bin(strn,grfac,4,0,' ',j,len(strn))) then
          print *, 'RELAX: Ignored command line value of grfac'
          grfac = 1.2d0
        endif
      endif
      xtoll = xtol
      if (cmdopt('-xtoll=',6,0,strn)) then
        j = 6
        if (.not. a2bin(strn,xtoll,4,0,' ',j,len(strn))) then
          print *, 'RELAX: Ignored command line value of xtoll'
          xtoll = xtol
        endif
      endif
C ... the user sets xtol _and_ gtol:
      if (lrlx == 1) isw = 00021 + 40
      if (lrlx == 2) isw = 00121 + 40
      if (lrlx == 3) isw = 00221 + 00
      if (xtol == 0 .and. gtol == 0)
     .  call rx('RELAX: both xtol and gtol are zero')
      if (gtol == 0) isw = isw-10
      if (xtol == 0) isw = isw-20
      call gradzr(nvar,p,w,xtoll,step,xtol,gtol,grfac,wkg,isw,ir)
      if (ir > 0) then
        call rx1('RELAX: gradzr returned ir=%i ... aborting',ir)
      endif
      if (ipr >= 20 .and. (lrlx == 1 .or. lrlx == 2)) then
        call awrit4(' RELAX line %d:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n==-1#new line minimization##'//
     .    '%-1j%?#n==-2#bracketed root this line##'//
     .    '%-1j%?#n==-3#bracketed root this line##'//
     .    '%-1j%?#n==-4#extrapolated along this line##'//
     .    '%-1j%?#n==-5#extrapolated along this line##'//
     .    '%-1j%?#n==-6#is in trouble##'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,i1mach(2),wkg(19),ir,
     .    wkg(1)*p(idamax(nvar,p(1,3),1),3),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
      elseif (ipr >= 20 .and. lrlx == 6) then
        call awrit4(' RELAX:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n<0#completed Broyden iter %i#%j#'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,i1mach(2),
     .    ir,-ir,p(idamax(nvar,p(1,5),1),5),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
      endif
      if (ir == 0) icom = 1

C --- Printout  ---
      if (ipr >= 40) then
        print *,'        Gradients:'
        print 100, (p(i,2), i = 1, nvar)
        if (lrlx /= 4) then
          print *,'      Diagonal inverse Hessian:'
          print 100, (w(i,i), i = 1,nvar)
        endif
      endif
  100 format(10f8.3)

C --- Update atom positions ---
      call dcopy(nvar,p,1,bas,1)

C --- Write Hessian to disc ---
      if (rdhess .and. (icom == 1 .or. it == nit) .or. .true.) then
        call mopen(ifi,'hssn','u')
C   ... Hessian written in rdm-compatible format
        rewind ifi
        write(ifi) nvar,nvar,11
        write(ifi) w
        call fclose(ifi)
      endif

      call poppr
      ipr = iprint()

C --- Printout ---
      ifi = 71
      lpp = (icom == 1 .or. it == nit) .and. ipr >= 30 .or.
     .      ipr >= 31
      if (lpp) then
        print 120
        print *,' Site   Species                      Position'
        do  70  j = 1, nbas
          write (*,130) j,atid(j),(bas(ix,j),ix=1,3)
   70   continue
      endif

C ... Write new atom positions to LOG file, for use in CTRL file
      dumstr = 'SITE'
      do  80  j = 1, nbas
        if (j == 2) dumstr = ' '
        call awrit1(dumstr//'ATOM='//atid(j)//
     .              ' POS=%3;12,7D  ',' ',80,ifi,bas(1,j))
   80 continue

      call ftflsh(i1mach(2))
      call ftflsh(ifi)

  120 format(/' Updated atom positions:')
  130 format(i4,6x,a4,3x,3(f14.8))
      end
