      subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .  metal,ltet,norder,npts,width,rnge,wtkp,eb,efmax,efermi,sumev,
     .  wtkb,dosef,qval,ent,egap)
C- BZ integration for fermi level, band sum and qp weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of ltethedron;
Ci         :idtet(0,i) number of ltethedra of the i'th kind
Ci   zval  :valence charge
Ci   metal :0 => nonmetal
Ci         :1 => metal
Ci         :11 => Compute weights using input efermi as Fermi level
Ci   ltet  :0 => sampling integration
Ci         :1 => tetrahedron integration
Ci         :11 => tetrahedron integration for bands, sampling for weights
Ci   norder:(sampling) polynomial order in Methfessel-Paxton integration
Ci         :100s digit norder flags that metals treatment should apply
Ci         :     regardless of whether a gap is present or not
Ci   width :(sampling) gaussian width in Methfessel-Paxton integration
Ci   npts  :(sampling) number of points in DOS mesh
Ci   rnge  :(sampling) range over which sampling delta function is assumed
Ci         :to vanish, in units of width
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which eigenvectors are known
Co Outputs
Co   efermi:Fermi energy (if metal is 11, efermi is input and not altered)
Co   sumev :sum of eigenvalues
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   dosef :DOS at Fermi level
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Co   ent   :entropy term (actually TS)
Co   egap  :energy gap (returned if metal=T and insulator found)
Co         :If none sought or metal found, egap returned as NULLI
Cs Command-line switches
Cs   --cv:    : Calculate electronic specific heat, eV
Cs   --cvK:   : Calculate electronic specific heat, Ry
Cs   --efrnge : Print out indices to bands that bracket the Fermi level
Cs   --oldbz  : Not documented
Cl Local variables
Cl   lfill :true => insulator
Cu Updates
Cu   15 Apr 16 Added --efrnge
Cu   25 Apr 14 Bug fix case metal=11 and insulator
Cu   11 Oct 13 metal now an integer, with allowed values 0,1,11
Cu   09 Aug 13 ltet now an integer, with allowed values 0,1,11
Cu   02 Jun 13 Replace f77 pointers with f90 ones
Cu   31 Aug 12 Returns egap
Cu   12 Jul 08 (ATP) bzwts now returns entropy term (actually kTS)
Cu   04 Jun 08 (ATP) Handles metal case when nkp=1
Cu    4 Aug 07 bzwts can make and tabulate specific heat (F-D statistics)
Cu   29 Jul 07 (ATP) Find E_F using weights by bisection, not INTNOS
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   17 Jan 05 Use 100s digit norder as flag to treat all cases as metal,
Cu             whether or not a gap is present
Cu    1 May 04 When insulator, write gap to log file
Cu   09 May 04 When insulator, write gap to log file
Cu   01 Jul 03 When insulator, prints highest occ and lowest unocc state
Cu   24 Oct 02 Patch for weird cases when idos doesn't encompass
Cu             zval, where emin, emax found by efrang.
Cu   22 Sep 01 Returns dosef now.  Altered argument list.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .  idtet(5,ntet),ltet,metal
      double precision zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),
     .  wtkb(nevx,nsp,nkp),efmax,efermi,sumev,dosef(2),qval(2),ent
C ... Dynamically allocated arrays
      integer, allocatable :: bmap(:)
      real(8), allocatable :: bot(:), top(:), dos(:), wk(:), tlst(:)
C ... Local parameters
      integer fopna,i,i1mach,ib,ifi,ikp,ipr,isw,it,itmax,j,job,lRy,
     .  master,mkdlst,mpipid,n,nbmxx,nbpw,nev,nevxx,nptdos,nspx,
     .  nulli,procid,stdl,stdo,nsk1(4),nsk2(4)
      double precision amom,cv,dmax,dmin,dsum,dum,e,e1,e2,egap,ehi,elo,
     .  emax,emin,sumwt,tRy,tol,xx,sumevs,efermis,dosefs(2),efloc
      character outs*100
      logical lfill
      procedure(logical) cmdopt,efrng2
      procedure(integer) idnear,nglob
C     procedure(real(8))
      parameter (nulli=-99999)
C ... External calls
      external awrit1,awrit2,awrit3,awrit4,awrit5,bzints,dcopy,dpsadd,
     .         dscal,ebcpl,fclose,fermi,fexit2,getpr,iinit,info0,info2,
     .         info5,intnos,maknos,mark1,poppr,pshprt,rxx,splwts,tcn,
     .         tcx,word,zcopy

C --- Locate band limits: lfill => insulator ---
      procid = mpipid(1)
      master = 0
      call tcn('bzwts')
      call getpr(ipr)
      qval(1) = 0
      qval(2) = 0
      ent = 0
      n = isign(1,norder) * mod(iabs(norder),100)
      stdo = nglob('stdo')
      stdl = nglob('stdl')
      allocate(bot(nevx*nsp),top(nevx*nsp))
      nspx  = 1
      nevxx = nevx*nsp
      nbmxx = nbmx*nsp
C     job = 1 for non spin pol, -1 for spin pol
      job = 3-2*nsp
      dosef(1) = 0
      dosef(nspx) = 0
      egap = nulli
      if (nsp == 2 .and. nspc == 1 .and. cmdopt('--oldbz',7,0,outs))
     .  then
        nspx  = nsp
        nevxx = nevx
        nbmxx = nbmx
        job = 1
      endif

C ... Force coupled spins: find range
      if (nspx /= nsp .and. nspc == 1) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        allocate(bmap(nevx*nsp*nkp/nbpw+1))
        call iinit(bmap,nevx*nsp*nkp/nbpw+1)
        allocate(wk(nevx*nsp))
        call ebcpl(0,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap,wk,eb)
        lfill = efrng2(nspx,nkp,nbmxx,nevxx,zval*2,eb,bot,top,
     .    elo,ehi,emin,emax)
        call ebcpl(1,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap,wk,eb)
        deallocate(bmap,wk)
C     Spins not coupled: find range
      else
        lfill = efrng2(nspx,nkp,nbmxx,nevxx,nspc*zval,eb,bot,
     .    top,elo,ehi,emin,emax)
      endif

C ... Find k-points corresponding to emin, emax calculated by efrng2
      if (cmdopt('--efrnge',8,0,outs)) then
        i = idnear(nbmx*nkp*nsp,emin,eb,1)
        call m4indx(nbmx,nsp,0,i,nsk1)
        call info5(10,1,-1,
     .    ' BZWTS:  e1=%d at band=%i ik=%i %?#(n==1)##isp=%i#',
     .    eb(nsk1(1),nsk1(2),nsk1(3)),nsk1(1),nsk1(3),nsp,nsk1(2))
        i = idnear(nbmx*nkp*nsp,emax,eb,1)
        call m4indx(nbmx,nsp,0,i,nsk2)
        call info5(10,0,0,'  e2=%d at band=%i ik=%i %?#(n==1)##isp=%i#',
     .    eb(nsk2(1),nsk2(2),nsk2(3)),nsk2(1),nsk2(3),nsp,nsk2(2))
      endif

C ... Bands never filled if 100s digit norder set
      if (ltet /= 1 .and. iabs(norder) >= 100) lfill = .false.
      deallocate(bot,top)
C ... Case an insulator: put efermi at emin + tiny number
      if (metal == 11) then
      elseif (lfill) then
        efloc = emin + 1d-10
        efermi = efloc
C ... Do the best we can should be a metal, but assumption that it isn't
      elseif (metal == 0) then
        efloc = (emin + emax) / 2
        efermi = efloc
      endif

C ... Pretend as though spin-pol bands are coupled to find E_f
C     Hang onto bmap, wk; they are used later
      if (nsp == 2 .and. nspc == 1 .and. job == -1) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        allocate(bmap(nevx*nsp*nkp/nbpw+1))
        call iinit(bmap,nevx*nsp*nkp/nbpw+1)
        allocate(wk(nevx*nsp))
        call ebcpl(0,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap,wk,eb)
      endif

C --- BZ weights, sumev and E_f for an insulator  ---
      if (metal == 0) then
        if (.not. lfill .and. ipr > 10) then
          write(stdo,*)' BZWTS : partly filled bands encountered;'
          write(stdo,*)' expect serious errors from Fermi cut-off.'
          write(stdo,*)' **** re-start with METAL=T in control file'
        endif
        sumwt = 0d0
        sumev = 0d0
        nev = nint(zval/2)
        call rxx(nev > nbmx,'BZWTS: zval too big')
        do  ikp = 1, nkp
          do  ib = 1, nev*nsp
          e = eb(ib,1,ikp)
          if (e <= efermi) then
            sumwt = sumwt +   abs(wtkp(ikp))/nsp
            sumev = sumev + e*abs(wtkp(ikp))/nsp
          endif
          enddo
        enddo
        egap = emax-emin
        if (ipr >= 20) then
          write(stdo,'(/'' BZWTS : --- Non-metal sampling ---'')')
          call awrit3(' Fermi energy: %,6;6d;  %d electrons;  Sum '//
     .      'occ. bands: %,6;6d',' ',80,stdo,efermi,sumwt,sumev)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .      '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .      (emax-emin)*13.6d0,0)
        endif

C --- BZ weights, sumev and E_f by tetrahedron method (Blochl wts) ---
      else if (ltet /= 0) then
        if (ipr >= 30) write (stdo,1)
    1   format(/' BZWTS : --- Tetrahedron Integration ---')
        if (lfill) then
          egap = emax-emin
          call info2(10,0,0,' ... only filled or empty bands'//
     .      ' encountered:  ev=%,6;6d  ec=%,6;6d',emin,emax)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .      '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .      (emax-emin)*13.6d0,0)
          if (metal == 11) goto 11
          call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,nspx,emin,emin,
     .                dosef,1,efloc,2*job,ntet,idtet,sumev,qval)
          goto 10
        endif
        nptdos = 101
        allocate(dos(nspx*nptdos))
        tol = 1d-6
C       Preliminary check that dos lies within emin,emax
C       Widen emin,emax if not
        if (.not. lfill) then
          call bzints(n1,n2,n3,eb,dum,nkp,nevxx,nbmxx,nspx,emin,emax,
     .                dos,nptdos,efloc,job,ntet,idtet,sumev,qval)
          dmin = dos(1)
          if (nspx == 2) dmin = dmin + dos(nptdos+1)
          dmax = dos(nptdos)
          if (nspx == 2) dmax = dmax + dos(nptdos+nptdos)
          if (dmin > zval) then
            emin = 3*emin-2*emax
            call info5(10,0,0,
     .        ' (warning): initial NOS ( %1,6;6g %1,6;6g ) does'//
     .        ' not encompass Q=%1;6g' ,dmin,dmax,zval,0,0)
          elseif (dmax < zval) then
            emax = 3*emax-2*emin
            call info5(10,0,0,
     .        ' (warning): initial NOS ( %1,6;6g %1,6;6g ) does'//
     .        ' not encompass Q=%1;6g' ,dmin,dmax,zval,0,0)
          endif
        endif
        if (ipr >= 35) write(stdo,2)
    2   format(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'DOS(E_f)')

   11   continue                ! Re-entry point for lfill = .true. and metal = 11
        if (metal == 11) then
          call info2(10,1,0, ' BZWTS : Fermi level externally '//
     .      'given as %d ...',efermi,0)
          efloc = efermi
          call pshpr(min(ipr,15))
          call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,nspx,emin,emin,
     .                dosef,0,efloc,2*job,ntet,idtet,sumev,qval)
          call poppr
        elseif (.not. lfill) then
          itmax = 5
          do   it = 1, itmax
          call bzints(n1,n2,n3,eb,dum,nkp,nevxx,nbmxx,nspx,emin,emax,
     .                dos,nptdos,efloc,job,ntet,idtet,sumev,qval)
          call fermi(zval,dos,nptdos,emin,emax,nspx,efloc,emin,
     .               emax,dosef)
          if (ipr >= 35) write (stdo,3) efloc,emin,emax,emax-emin,
     .                            dosef(1)
    3     format(7x,6(f10.6,1x))
          if (emax-emin < tol) goto 5
          enddo
          if (ipr > 10)
     .    call awrit2(' BZWTS (warning): Fermi energy not '//
     .    'converged: %g > tol=%g',' ',80,i1mach(2),emax-emin,tol)
    5     continue
          deallocate(dos)
          call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,nspx,emin,emin,
     .                dosef,1,efloc,2*job,ntet,idtet,sumev,qval)
          efermi = efloc
        endif
   10   continue                ! Re-entry point for lfill = .true.
      endif  ! end of tetrahedron branch

C --- BZ weights, sumev and E_f by Methfessel-Paxton sampling ---
      if (metal /= 0 .and. (ltet == 0 .or. ltet == 11)) then
        if (ltet == 11) then
          sumevs = sumev
          efermis = efermi
          dosefs = dosef
        endif

        call info2(31,1,0,' BZWTS : --- Brillouin Zone sampling; '//
     .    'N=%i, W=%10z%;4g%0z ---',n,width)
C   ... Temporarily remove spin degeneracy if spins are coupled
        if (nsp == 2 .and. nspx == 1) call dscal(nkp,.5d0,wtkp,1)
C   ... Find Fermi level, sampling
        if ((.not. lfill) .or. (metal /= 0 .and. (nkp == 1))) then
          e1 = elo - rnge*width/2
          e2 = ehi + rnge*width/2
          efloc = 0.5d0*(e1 + e2)
          itmax = 1000
          do  it = 1, itmax
            call pshprt(0)
            call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efloc,
     .                  .true.,sumev,wtkb,qval,ent,dosef,cv)
            call poppr
            if (dabs(zval - qval(1)) < 1d-12) then
              call info5(41,0,0,
     .         ' Fermi energy, %d found after %i bisections,'//
     .         ' %d electrons, DOS(E_f)=%;4g',efloc,it,qval(1),dosef,0)
              if (metal == 11) then
                call info2(10,1,0, ' BZWTS : Fermi level externally '//
     .            'given as %d ... remake band weights',efermi,0)
                efloc = efermi
                call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,
     .            .true.,sumev,wtkb,qval,ent,dosef,cv)
              else
                efermi = efloc
              endif
              goto 15
            endif
            if (qval(1) > zval) then
              e2 = efloc
            else
              e1 = efloc
            endif
            efloc = 0.5d0*(e1 + e2)
          enddo
          call info0(10,0,0,' BZWTS (warning): '//
     .      'cannot find E_F by bisection, using INTNOS')
          allocate(dos(nspx*npts))
          emin = elo - rnge*width/2
          emax = emax + rnge*width/2
          call maknos(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,-rnge,
     .                emin,emax,npts,dos)
          if (nspx == 2) call dpsadd(dos,dos,npts,1,npts+1,1d0)
          call intnos(npts,dos,emin,emax,zval,efloc,dosef,sumev)
          if (metal == 11) call rx(
     .      'bzints: sampling integration not ready for metal=11')
          efermi = efloc
          deallocate(dos)
        else
          dosef(1) = 0
          dosef(2) = 0
          egap = emax-emin
          call info2(10,0,0,' ... only filled or empty bands'//
     .      ' encountered:  ev=%,6;6d  ec=%,6;6d',emin,emax)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .      '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .      (emax-emin)*13.6d0,0)
        endif

C   ... (optional) Tabulate specific heat in file for list of T's
   15   continue
        if ((cmdopt('--cv:',5,0,outs) .or. cmdopt('--cvK:',6,0,outs))
     .    .and. n < 0 .and. metal /= 0) then
        if (procid == master) then
          if (cmdopt('--cvK:',6,0,outs)) then
            lRy = 0
            i = 7
          else
            lRy = 1
            i = 6
          endif
          itmax = mkdlst(outs(i:),1d-8,-1,xx)
          if (itmax > 0) then
          allocate(tlst(itmax))
          call word(outs,1,it,j)
          call info2(20,0,0,' Writing CV(T) to file for %i vals of T: '
     .      //outs(i:j)//' %?#(n==1)#(Ry)#(K)#',itmax,lRy)
          it = mkdlst(outs(i:),1d-8,itmax+1,tlst)
          if (it /= itmax) call rx('bzwts: bug in mkdlst')
          ifi = fopna('cv',-1,0)
          rewind ifi
          call awrit1('%% rows %i cols 4%N'//
     .      '#    T(K)    T(Ry)     S(k_B)   TdS/dT(k_B)',' ',80,ifi,it)
          do  it = 1, itmax
            tRy = tlst(it)
            if (lRy == 0) then
              tRy = tRy/0.1579d6
            endif
            call pshpr(1)
            call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,tRy,efermi,
     .        metal /= 0,sumev,wtkb,qval,ent,dosef,cv)
            call poppr
            call awrit4('%;9,1D  %;9F  %;9F  %;9F',' ',80,ifi,0.1579d6*tRy,tRy,ent,cv)
          enddo
          deallocate(tlst)
          call fclose(ifi)
          endif
        endif
        endif

C   ... Make weights, sampling
        call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,
     .              (.not. lfill) .or. (metal /= 0 .and. nkp == 1),
     .              sumev,wtkb,qval,ent,dosef,cv)

C   ... Put back spin degeneracy if removed
        if (nsp == 2 .and. nspx == 1) call dscal(nkp,2d0,wtkp,1)

        if (ltet == 11) then
          sumev = sumevs
          efermi = efermis
          dosef = dosefs
        endif

      endif  ! End of sampling branch

C ... Restore to uncoupled bands; ditto with weights
      if (nsp == 2 .and. nspc == 1 .and. job == -1) then
        call ebcpl(1,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap,wk,eb)
        if (metal /= 0)
     .  call ebcpl(1,nevx,nevx,nsp,nspc,nkp,nbpw,bmap,wk,wtkb)
        deallocate(bmap,wk)
      endif

C ... Magnetic moment
      amom = 0
      if (nsp == 2 .and. nspc /= 2 .and. metal /= 0) then
        do  ikp = 1, nkp
          amom = amom + dsum(nevx,wtkb(1,1,ikp),1) -
     .                  dsum(nevx,wtkb(1,2,ikp),1)
        enddo
C        if (ipr > 1) write(stdo,922) amom
C  922   format(9x,'Mag. moment:',f15.6)
C        if (ipr > 1) then
C          call awrit1(' BZWTS: magnetic moment = %d',' ',128,i1mach(2),
C     .                 amom)
C        endif
      endif
      qval(2) = amom

C ... log file
      if (ipr > 1) then
      outs = ' '
      call awrit5(
     .  'bz met %l  tet %l  ef %,6;6d  sev %,9;9d  zval %,1;6d',
     .  outs,len(outs),0,metal /= 0,ltet /= 0,efermi,sumev,zval)
      call awrit5('%a qval %,1;6d%?#n# amom %,1;6d#%j#'//
     .  '%?#n# gap %,4;4d eV##',outs,len(outs),0,
     .  qval,int(amom*10000),amom,isw(egap > 0),egap*13.6d0)
      if (metal /= 0)call awrit1('%a dosef %;4g',outs,len(outs),0,dosef)
      call awrit1('%a',outs,len(outs),-stdl,dosef)
      endif

C ... Check that efmax is above efermi ...
      e = efermi
      if (.not. lfill .and. ltet /= 1) e = efermi + rnge*width/2
      if (e > efmax) call fexit2(-1,111,' Exit -1 BZWTS: '//
     .  'efmax=%1;6d < %1;6d',efmax,e)

      call tcx('bzwts')

      end

      subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb)
C- Gather spin-polarized bands into a single group, or redistribute
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0, make bmap and gather; 1, scatter
Ci         :2  like 0 (make bmap) but do not gather
Ci   nbmx  :leading dimension of b and dimensions wk
Ci   nevx  :number of eigenvalues
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nq    :number of k-points for which eb are calculated
Ci   nbpw  :a number no larger than the number of bits per integer word
Ci   bmap  :an integer work array with dimension at nbmx*nsp*nq/nbpw
Ci         :see Remarks
Ci   wk    :a work array with dimension nbmx*nsp
Ci         :Not touched if mode=2
Cio Inputs/Outputs
Cio   eb    :energy bands:
Cio         :mode=0 input spin-split, output merged to a single vector
Cio         :mode=1 input merged to a single vector, output spin-split
Cio         :Not touched if mode=2
Cr Remarks
Cio   Call ebcpl with mode=1 to undo call of ebcpl with mode=0.
Cio   bmap must be preserved for mode=1 call.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(*)
      double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
C ... Local parameters
      integer ib,iq,ib1,ib2,iqb,iget

      if (nsp == 1 .or. nspc == 2) return

C --- Gather bands at each qp ---
      if (mode == 0 .or. mode == 2) then
        iqb = 0
        do  iq = 1, nq

C   ... Gather and order +,- bands at this qp into one column
        ib1 = 1
        ib2 = 1
          do  ib = 1, nevx*nsp
          iqb = iqb+1
          if (eb(min(ib1,nevx),1,iq) < eb(min(ib2,nevx),2,iq)
     .       .and. ib1 <= nevx .or. ib2 > nevx) then
            if (mode == 0) wk(ib) = eb(ib1,1,iq)
            ib1 = ib1+1
          else
            if (mode == 0) wk(ib) = eb(ib2,2,iq)
            ib2 = ib2+1
            call mark1(bmap, nbpw, iqb)
          endif
C         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
C    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),
C    .      iget(bmap,nbpw,iqb),ib1,ib2)
          enddo
          if (mode == 0) call dcopy(nevx*nsp,wk,1,eb(1,1,iq),1)
          if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
        enddo
      endif

C --- Disperse bands at each qp ---
      if (mode == 1) then
        iqb = 0
        do  iq = 1, nq

C   ... Disperse bands into +,- for this qp according to bmap
        ib1 = 1
        ib2 = 1
        call dcopy(nevx*nsp,eb(1,1,iq),1,wk,1)
          do  ib = 1, nevx*nsp
          iqb = iqb+1
          if (iget(bmap,nbpw,iqb) == 0) then
            eb(ib1,1,iq) = wk(ib)
C            call awrit4(' iq=%,2i  ib=%,2i  evl=%d  down ib1=%i',
C     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib1)
            ib1 = ib1+1
          else
            eb(ib2,2,iq) = wk(ib)
C            call awrit4(' iq=%,2i  ib=%,2i  evl=%d    up ib2=%i',
C     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib2)
            ib2 = ib2+1
          endif
C         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
C    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),iget(bmap,nbpw,iqb),
C    .      ib1,ib2)
          enddo
        if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
        enddo
      endif

      end
      subroutine ebcpl2(mode,ndz,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,d,z)
C- Gather spin-polarized evecs according to bmap
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0, gather; 1, scatter
Ci         :Add 10 to work with z instead of d
Ci   ndz   :leading dimension of z; size of eigenvector
Ci   nbmx  :dimensions of z and wk
Ci   nevx  :number of eigenvalues
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nq    :number of k-points for which eb are calculated
Ci   nbpw  :a number no larger than the number of bits per integer word
Ci   bmap  :a map array generated by ebcpl, which see
Ci   wk    :a work array with dimension nbmx*nsp
Cio Inputs/Outputs
Cio   d    :eigenvector (real)
Cio   z    :eigenvector (complex)
Cio        :mode=0 input spin-split, output merged to a single vector
Cio        :mode=1 input merged to a single vector, output spin-split
Cr Remarks
Cio   Call ebcpl with mode=1 to undo call of ebcpl with mode=0.
Cio   bmap must be preserved for mode=1 call.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ndz,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(*)
C     double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
      double precision d(ndz,nbmx,nsp,nq)
      double complex z(ndz,nbmx,nsp,nq)
C ... Local parameters
      integer ib,iq,ib1,ib2,iqb,iget
      complex(8), allocatable :: zwk(:,:)

      if (nsp == 1 .or. nspc == 2) return
      allocate(zwk(ndz,nbmx*nsp))
      if (mode < 10) call rx('ebcpl2 not ready')

C --- Scatter/Gather evecs at each qp ---
      if (mod(mode,10) == 0) then
        iqb = 0
        do  iq = 1, nq

C   ... Gather and order +,- bands at this qp into one column
        ib1 = 1
        ib2 = 1
        do  ib = 1, nevx*nsp
          iqb = iqb+1
          if (iget(bmap,nbpw,iqb) == 0) then
            zwk(:,ib) = z(:,ib1,1,iq)
            ib1 = ib1+1
          else
            zwk(:,ib) = z(:,ib2,2,iq)
            ib2 = ib2+1
          endif
        enddo
        call zcopy(ndz*nevx*nsp,zwk,1,z(1,1,1,iq),1)
        if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
      enddo
      endif

C --- Disperse bands at each qp ---
      if (mod(mode,10) == 1) then
        iqb = 0
        do  iq = 1, nq

C   ... Disperse bands into +,- for this qp according to bmap
        ib1 = 1
        ib2 = 1
        call dcopy(ndz*nevx*nsp,z(1,1,1,iq),1,zwk,1)
        do  ib = 1, nevx*nsp
          iqb = iqb+1
          if (iget(bmap,nbpw,iqb) == 0) then
            zwk(:,ib) = z(:,ib1,1,iq)
            ib1 = ib1+1
          else
            zwk(:,ib) = z(:,ib2,2,iq)
            ib2 = ib2+1
          endif
          enddo
        if (ib1-1 /= nevx .and. ib2-1 /= nevx) call rx('bug')
        enddo
      endif
      deallocate(zwk)

      end
      subroutine mark1(bitmap, nbpw, n)
C- put a one in the nth bit of bitmap.
C ----------------------------------------------------------------
Ci Inputs
Ci   bitmap, n
Cr Remarks
Cr    nbpw: a number no larger than the number of bits per integer word
C ----------------------------------------------------------------
      implicit none
      integer bitmap(*), nbpw, n
C Local parameters
      integer nword,nbit,i

      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      i = 2**(nbpw-nbit-1)
      bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
      end

      integer function iget(bitmap, nbpw, n)
C- Return 0 or 1, depending on the value of the nth bit of bitmap
C ----------------------------------------------------------------
Cr Remarks
Cr   See mark1
C ----------------------------------------------------------------
      implicit none
      integer bitmap(*), nbpw, n
C Local parameters
      integer nword,nbit
      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
      end

      subroutine m4indx(n1,n2,n3,i,n123)
C- Unpack the indices of a 3D or 4D array given the absolute position
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1    : leading dimension of array
Ci   n2    : second dimension of array
Ci   n3    : third dimension of array.  Set to 0 if array has only 3D.
Ci   i     : absolute position of element relative to start of array
Co Outputs
Co   n123  : 3D case (n3=0). If A = A(n1,n2,n3) then
Co         : A(i,1,1) points to A(n123(1),n123(2),n123(3))
Co         : 4D case not implemented
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   15 Apr 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,i,n123(3)
C ... Local parameters
      integer k

      if (n3 == 0) then
        n123(3) = (i-1)/(n1*n2) + 1
        k = i - n1*n2*(n123(3)-1)
        n123(2) = (k-1)/(n1) + 1
        k = k - n1*(n123(2)-1)
        n123(1) = k
      else
        call rx('m4indx: not implemented')
      endif
      end
