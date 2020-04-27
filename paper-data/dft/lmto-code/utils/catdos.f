C  Combines sets of densities of states into a single file, standard DOS format
Cu Updates
Cu   28 May 17 scale can be adjusted for each dos
      subroutine fmain
      implicit none
      integer nemx,ndmx,nspin,ndos,nld,ndos2,nspin2,nld2
      parameter( nemx=5001, ndmx=300)
      double precision d(nemx,ndmx),d2(nemx)
      double precision emin,emax,eferm,emin2,emax2,eferm2,del2,del,de,ef,scale,dmax
      integer iarg,ifmt,ifi,i,j,n,i1mach
      character*20 fname
      character*1 f2(20)
      equivalence (f2,fname)
      logical cmdstr,lsequ,a2bin,lef
      procedure(real(8)) :: dval
      common /zz/ d
      data emin /-99/, emax /-99/, nspin /-1/, ndos /-99/

      call dpzero(d,size(d))
      lef = .false.
      scale = 1
      iarg = 0
      goto 12
   10 print 332
  332 format(
     .  ' usage: catdos [-l#] [-h#] [-n#]',
     .  ' [-s#] [-ef=#] file [[-s#] [-ef=#] file ...]'/
     .  '        combines DOS files into one file (standard DOS format)'/
     .  '        The following three switches cause DOS to be interpolated'/
     .  '        and should be appear before any DOS file'/
     .  '        -l#   set lower energy bound for DOS'/
     .  '        -h#   set upper energy bound for DOS'/
     .  '        -n#   set number of DOS points'/
     .  '        The following apply to any DOS read following switch'/
     .  '        -ef=# aligns DOS to specified Fermi level=#'/
     .  '        -s#   scales DOS subsequently read by #')
      call cexit(-1,1)

   12 iarg = iarg+1
      if (.not. cmdstr(iarg,fname)) goto 10
      if (lsequ(fname,'-l',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),emin,4,0,' ',i,-1)) goto 10
        goto 12
      else if (lsequ(fname,'-h',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),emax,4,0,' ',i,-1)) goto 10
        goto 12
      else if (lsequ(fname,'-n',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),ndos,2,0,' ',i,-1)) goto 10
        goto 12
      else if (fname(1:4) == '-ef=') then
        lef = .true.
        i = 0
        if (.not. a2bin(f2(5),ef,4,0,' ',i,-1)) goto 10
        goto 12
      else if (lsequ(fname,'-s',2,' ',n)) then
        i = 0
        if (.not. a2bin(f2(3),scale,4,0,' ',i,-1)) goto 10
        goto 12
      endif

C --- Read until files exhausted ---
      print 345
  345 format('   File',18x,
     .  'ndos  nstate     scale    max D       emin      emax     old-ef     shift')
      nld = 0
      iarg = iarg-1
   15 iarg = iarg+1
        if (.not. cmdstr(iarg,fname)) goto 20
        if (fname(1:4) == '-ef=') then
          lef = .true.
          i = 0
          if (.not. a2bin(f2(5),ef,4,0,' ',i,-1)) goto 10
          goto 15
        endif
        if (lsequ(fname,'-s',2,' ',n)) then
          i = 0
          if (.not. a2bin(f2(3),scale,4,0,' ',i,-1)) goto 10
          goto 15
        endif

        ifi = 10
        open(ifi, file=fname, status='OLD', err=10)
        call dosio(d(1,nld+1),nemx,ndmx-nld,ndos2,nld2,emin2,emax2,
     .    nspin2,eferm2,del2,ifmt,ifi)
        call dscal(nemx*(ndmx-nld),scale,d(1,nld+1),1)
        call idmxmn(nemx*(ndmx-nld),d(1,nld+1),1,i,j)
        dmax = dval(d(1,nld+1),j)
C   --- Get default values for emin, emax, ndos from first dos ---
C       nspin=-1 flags that no dos read yet
        if (nspin == -1) then
          if (emin == -99) emin = emin2
          if (emax == -99) emax = emax2
          if (ndos == -99) ndos = ndos2
          del  = del2
          nspin = nspin2
          eferm = eferm2
C       If this happens, treat list as though all nspin=1
        elseif (nspin /= nspin2) then
          print *, ' warning ... mixing files w/ and w/out spin pol'
          nspin = 1
        endif

C   --- Align DOS to the specified fermi level if requested ---
        if (lef) then
          de = ef - eferm2
          emin2  = emin2+de
          emax2  = emax2+de
          eferm2 = eferm2+de
C          do  14  i = nld+1, nld+nld2*nspin2
C            call dcopy(ndos2,d(1,i),1,d2,1)
C            call dositp(d(1,i),ndos,emin,emax,ndos2,emin2,emax2,d2)
C   14     continue
        endif

C   --- Interpolate DOS to mesh specified by ndos,emin,emax ---
        nld = nld + nld2*nspin2
        if (lef) then
          print 334, fname, nld2, nspin2, ndos2, scale, dmax, emin2, emax2,
     .      eferm2-de, de
  334     format(2x,a20,i5,'(',i1,')',i8,7f10.5)
        else
          print 337, fname, nld2, nspin2, ndos2, scale, dmax, emin2, emax2
  337     format(2x,a20,i5,'(',i1,')',i8,4f10.5,5x,'---')
        endif
        if (ndos /= ndos2 .or. emin /= emin2 .or. emax /= emax2)
     .    then
          call awrit5(' ... interpolate dos %i..%i to ndos,emin,emax='//
     .      ' %i %1;6d %1;6d',' ',80,i1mach(2),
     .      nld+1-nld2*nspin2,nld,ndos,emin,emax)
          do  i = nld+1-nld2*nspin2, nld
            call dcopy(ndos2,d(1,i),1,d2,1)
            call dositp(d(1,i),ndos,emin,emax,ndos2,emin2,emax2,d2)
          enddo
        endif
        close(ifi)
        de = 0
      goto 15

   20 continue

C --- Cleanup ---
      call info2(0,0,0,' Writing %i(%i) dos to file dos.dat',
     .  nld/nspin,nspin)
      open(ifi, file='dos.dat', err=10)
      if (.not. lef) ef = eferm
      call dosio(d,nemx,ndmx,ndos,nld/nspin,emin,emax,nspin,
     .  ef,del,ifmt,-ifi)

      end
      subroutine dositp(dos,ndos,emin,emax,ndos2,emin2,emax2,d2)
C- Interpolate density (or number) of states to another mesh
C ----------------------------------------------------------------
Ci Inputs
Ci   d2,ndos2,emin2,emax2: starting number of states,
Ci     number of mesh points and energy range of dos
Ci   ndos,emin,emax: output number of states,
Ci     number of mesh points and energy range of dos
Co Outputs
Co   dos: interpolated dos
Cr Remarks
C ----------------------------------------------------------------
C Passed parameters
      integer ndos,ndos2
      double precision dos(1),d2(1),emin,emax,emin2,emax2
C Local parameters
      integer i,iold
      double precision xiold,enew,sloold,slonew

c      if (emin < emin2) stop 'dositp: cannot extrapolate emin'
c      if (emax > emax2) stop 'dositp: cannot extrapolate emax'

      sloold = (emax2-emin2)/(ndos2-1)
      slonew = (emax -emin )/(ndos -1)
      do  10  i = 1, ndos
        enew = slonew * (i-1) + emin
        xiold = (enew-emin2+.0000001)/sloold + 1
        iold  = xiold
        if (iold <= 0) iold=1
        if (iold >= ndos2) iold=ndos2-1
        xiold = xiold-iold
        dos(i) = d2(iold) + xiold*(d2(iold+1)-d2(iold))
   10 continue
      end
      subroutine dosio(dos,nemx,ndmx,ndos,nld,emin,emax,nspin,
     .   eferm,del,ifmt,ifile)
C- I/O for DOS, MSM's format
C ----------------------------------------------------------------
Ci Inputs
Ci   ifile
Cr Remarks
Cr   dos is assumed to be declared as dos(nemx,nemx), but
Cr   dos(1..ndos,1..nld) are input/output
Cr
Cr   ifmt=0 for old format, ifmt=1 for new format
C ----------------------------------------------------------------
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
        if (ndos > nemx) call rx('dosio: ndos gt nemx')
        if (nld*nspin > ndmx) call fexit2(-1,001,' Exit -1 DOSIO: '//
     .    'nld (=%i) > ndmx (=%i)',nld*nspin,ndmx)
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
