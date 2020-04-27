C#define PRTNOCR
      subroutine rsedita(sopts,mode,s_ctrl,s_site,s_spec,s_lat,s_pot,
     .  s_bz,nbas,nat,nspec)
C- Restart file editor for ASA restart files
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: zbak lrs
Co     Stored:    lrs
Co     Allocated: *
Cio    Elts Passed:lrs
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat nabc ag bgv cg cy dlv gv gvq indxcg ips0
Ci                istab jcg kv igv igv2 kv2 pos qlv symgr
Co     Stored:    alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                jcg kv igv igv2 kv2 pos qlv symgr
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: nlml nlma ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf
Ci                dpfr gibbs gma gmar grrme mad mxy orhat palp papg pf
Ci                pfnc pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp
Ci                qt rhos rhrmx sop thetcl vdif vintr vrmax vshft smpot
Ci                smrho smrout
Co     Stored:    ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf dpfr
Co                gibbs gma gmar grrme mad mxy orhat palp papg pf pfnc
Co                pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Co                rhos rhrmx sop thetcl vdif vintr vrmax vshft smpot
Co                smrho smrout
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: iorsa bcast_strx iinit mpibc1
Ci Inputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Ci   mode  :0 ASA I/O
Ci         :1 FP  I/O, mesh density
Ci         :-1 (input only) read nbas,nat,nspec from file
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Co Outputs
Co   rsedit never returns.
Co   rst file can be written.
Cr Remarks
Cr   The density consists of a smooth part (smrho) plus
Cr   nbas atom-centered densities inside the MT spheres.
Cr   Their sum is the full charge density.
Cr   The local density is represented as the difference of the
Cr   two valence components in orhoat, plus the core density.
Cr   Density in the MT spheres:
Cr      mesh parameters rmt,nr,a;
Cr      total density rho (times r**2) to lmxl;
Cr      a spherical potential v0 defining the wave functions within rmt
Cr      pnu and idmod to lmxa
Cr   Smooth density
Cr      real part of complex*16 array smrho contains the density
Cr      k1,k2,k3 are the physical dimensions of the array
Cr      n1,n2,n3 are the dimensions of the mesh.
Cl Local variables
Cu Updates
Cu   25 Aug 16 Bug fix: pnu is updated properly in 'set' command
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   29 Jun 11 New allnc option
Cu   10 Mar 11 new rmsite and ssave options
Cu   19 Dec 08 New exch option
Cu   01 Jul 08 First created
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer mode,nbas,nat,nspec,ifi,n0
      parameter (n0=10)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      integer,allocatable:: ilst(:,:)
C     real(8),allocatable:: pnu1(:,:,:),pnu2(:,:,:)
      real(8),allocatable:: qnu1(:,:,:,:),qnu2(:,:,:,:)
C ... Local parameters
      type(str_lat)::   s_lat2
      type(str_pot)::   s_pot2
      type(str_site),allocatable::  s_site1(:),s_site2(:)
      type(str_spec),allocatable::  s_spec1(:),s_spec2(:)
      integer nglob,fopna,fopng,iorsa,a2vec
      integer i,j,k,l,im,m,n,j1,j2,js1,js2,nit(2),nspeck(2),isync,nw,
     .  lexch,nspec10,nbas10,specsize,ssizei
      logical lnsave,lbin,lsopts,have2s,havers,lfail
      integer nsp,ngabc(3),n1,n2,n3,k1,k2,k3,lrssav,nvl,nchan
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer nlst(2),ix(10),nbask(2),natk(2),isw
C      integer idup(3)
      double precision qbg,fac(2),alat,plat(3,3),xx
C     double precision plat2(3,3)
      character dc*1, fn*120, fileid*68, outs*150
C     double precision dval
C ... data statements
C     data vec0 /0d0,0d0,0d0/

C     Defaults
      fn = 'rst1'
      havers = .false.
      have2s = .false.
C     # spec and # sites in spec 2nd atom density
      nbask(1) = nbas
      natk(1) = nat
      nspeck(1) = nspec
      nbask(2) = nbas
      nspeck(2) = nspec
      natk(2) = nat
C     true if rs data available, but not saved
      lnsave = .false.
      nsp  = nglob('nsp')
      qbg  = s_ctrl%zbak(1)
      nvl = s_pot%nlml
      nchan = s_pot%nlma
      alat = s_lat%alat
      plat = s_lat%plat

C     Allocate local species structures; copy original to them.
C     Make space for extra species
      nspec10 = nspec+10
      allocate(s_spec1(nspec10))
      allocate(s_spec2(nspec10))
      specsize = ssizei(s_spec2(1),s_spec2(2))
      s_spec1(1:nspec) = s_spec(1:nspec)
      s_spec2(1:nspec) = s_spec(1:nspec)
C     Allocate local site structures; copy original to them
C     Make space for extra sites
      nbas10 = nbas+10
      allocate(s_site1(nbas10))
      allocate(s_site2(nbas10))
      s_site1(1:nbas) = s_site(1:nbas)
      s_site2(1:nbas) = s_site(1:nbas)
C     Allocate second lat structure; copy original to it
      s_lat2 = s_lat
C     Allocate second pot structure; copy original to it
      s_pot2 = s_pot

C     allocate(pnu1(n0,nsp,nbas10),pnu2(n0,nsp,nbas10))
      allocate(qnu1(3,n0,nsp,nbas10),qnu2(3,n0,nsp,nbas10))

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the ASA restart file editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the ASA restart file editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif
      if (mode /= 1) call rx('rsedit not ready for mode ne 1')

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 > len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) == dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 < js1) lsopts = .false.
      endif

C 306 format(' Failed to parse string ',a,' ... try again.')
C 100 continue
C#ifdef PRTNOCR
      print '(/'' Option : '',$)'
C#elseC
C      print '(/'' Option : '')'
C#endif
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .    ' ''?'' to see menu')
        goto 10

C ... Read 2nd density from restart file
      elseif (outs(1:5) == 'read ' .or. outs(1:6) == 'reada ') then
        call word(outs,2,j1,j2)
        if (outs(1:6) == 'reada ' .or. .true.) then
          if (j2 < j1) then
            fn = 'rsta'
            ifi = fopna('rsta',-1,1)
          else
            fn = outs(j1:j2)
            ifi = fopng(outs(j1:j2),-1,1)
          endif
          lbin = .false.
C        else
C          if (j2 < j1) then
C            fn = 'rst'
C            ifi = fopna('rst',-1,4+1)
C          else
C            fn = outs(j1:j2)
C            ifi = fopng(outs(j1:j2),-1,4+1)
C          endif
C          lbin = .true.
        endif
        lrssav = s_ctrl%lrs
        s_ctrl%lrs = 32*0+64+128*0
C       Get the number of atoms and species
        k = iorsa(-1,s_ctrl,xx,xx,s_lat2,s_pot2,s_bz,
     .    fileid,nbask(2),natk(2),nspeck(2),xx,nit(2),lbin,ifi)
        call info5(10,0,0,' rsedit: '//
     .    'read 2nd density from file "'//trim(fn)//'" ...  '//
     .    'nbas = %i  nat = %i  nspec = %i',nbask(2),natk(1),nspeck(2),
     .    0,0)
C       Allocate species and site structures, fill with available info
        deallocate(s_spec2,s_site2)
        allocate(s_spec2(nspeck(2)))
        allocate(s_site2(nbask(2)))
        call bcast_strx(-1-2,xx,xx,xx,xx,xx,xx,s_spec2,s_site2,xx,
     .    nspeck(2),nbask(2))
        do  i = 1, nspeck(2)
          s_spec2(i)%p = 0
          s_spec2(i)%pz = 0
        enddo
        call pshpr(1)
        k = iorsa(-4,s_ctrl,s_site2,s_spec2,s_lat2,s_pot2,s_bz,
     .    fileid,nbask(2),natk(2),nspeck(2),qnu2,nit(2),lbin,ifi)
C       Check for unassigned species => reduce basis
        k = 0
        do  i = 1, nbask(2)
          j = s_site2(i)%spec
          if (j == 0) exit
          k = i
        enddo
        if (k /= nbask(2)) then
          call info2(10,0,0,'%9f(warning) only %i sites had '//
     .      'associated species ... reducing nbas(2)',k,0)
        endif
        nbask(2) = k
        call poppr
        k = iorsa(0,s_ctrl,s_site2,s_spec2,s_lat2,s_pot2,s_bz,
     .    fileid,nbask(2),natk(2),nspeck(2),qnu2,nit(2),lbin,ifi)
        call fclose(ifi)
        have2s = .true.

        call info0(2,1,0,' Second Density:')
        call prsed2a(0,nbask(2),nsp,s_site2,s_spec2,qnu2)

        goto 10

C ... Read restart file
      elseif (outs(1:3) == 'rs ' .or. outs(1:4) == 'rsa ') then
        call word(outs,2,j1,j2)
        if (outs(1:4) == 'rsa ' .or. .true.) then
          if (j2 < j1) then
            fn = 'rsta'
          else
            fn = outs(j1:j2)
          endif
          ifi = fopna(fn,-1,1)
          lbin = .false.
C        else
C          if (j2 < j1) then
C            fn = 'rst'
C          else
C            fn = outs(j1:j2)
C          endif
C          ifi = fopna(fn,-1,4+1)
C          lbin = .true.
        endif
        call info0(10,0,-1,
     .    ' rsedit: read restart file "'//trim(fn)//'"')
        k = iorsa(0,s_ctrl,s_site1,s_spec1,s_lat,s_pot,s_bz,fileid,
     .    nbask(1),natk(1),nspeck(1),qnu1,nit,lbin,ifi)
C        do  k = 1, nbask(1)
C          pnu1(1:n0,1:2,k) = s_site(k)%pnu(1:n0,1:2)
C        enddo
        call info5(10,0,0,'%9pfile has nbas = %i  nat = %i  '//
     .    'nspec = %i',nbask(1),natk(1),nspeck(1),0,0)

        call fclose(ifi)

        havers = .true.
        lnsave = .true.
        call info2(2,1,0,' Density from '//
     .    '%?#n==0#rs file:#overlapping free atoms:#',
     .    isw(outs(1:5) == 'rsfa '),0)
        call prsed2a(0,nbask(1),nsp,s_site1,s_spec1,qnu1)

C ... show
      elseif (outs(1:5) == 'show ') then
        if (havers) then
          call info2(2,1,0,' Density from rs file:  '//
     .      '%i sites, %i species',nbask(1),nspeck(1))
          call prsed2a(1,nbask(1),nsp,s_site1,s_spec1,qnu1)
        endif

        if (have2s) then
          call info2(2,1,0,' Second Density:  '//
     .      '%i sites, %i species',nbask(2),nspeck(2))
          call prsed2a(1,nbask(2),nsp,s_site2,s_spec2,qnu2)
        endif

C ... Spin averaging of densities, or  Exchange l=0 parts of
C     set object 1|2 list1 list2
C     exch 1|2 site1 site2  or exch 1|2 flip site1 site2
      elseif (outs(1:4) == 'set ' .or. outs(1:5) == 'exch ') then

      lexch = 0
      if (outs(1:5) == 'exch ') lexch = 1
      if (lexch == 1) stop 'not ready for exch'
      call words(outs,nw)
      if (nw < 4) goto 98

C     Switches for exchange
      if (lexch /= 0) then

C     Switches for set
      else
        call word(outs,2,j1,j2)
        im = -1
        if (outs(j1:j2) == 'all') im = 2+1
        if (outs(j1:j2) == 'qnu') im = 2
        if (outs(j1:j2) == 'pnu') im = 1
        if (im == -1) then
          call info0(0,0,0,
     .      '     ... mode "'//outs(j1:j2)//'" not recognized')
          goto 98
        endif

C       Fourth argument same for each density
        call word(outs,4,j1,j2)
        if (outs(j1:j2) == 'zers') then
          k = im*10+0
        elseif (outs(j1:j2) == 'zerq') then
          k = im*10+1
        elseif (outs(j1:j2) == 'flip') then
          k = im*10+2
        else
          call info0(0,0,0,' failed to recognize 4th argument in: '//
     .      trim(outs))
          goto 98
        endif

C       Fifth (site list) argument same for each density
        call word(outs,5,j1,j2)
        nlst(1) = 0
        if (j2 >= j1) then
          call mkils0(outs(j1:j2),nlst(1),j)
          if (nlst(1) <= 0) then
            call info0(0,0,0,'%6p... Bad or null list : '//outs(j1:j2))
            goto 98
          endif
          allocate(ilst(nlst(1),1))
          call mkilst(outs(j1:j2),nlst(1),ilst)
        endif

C       Carry out operation on 1st or second density
        call word(outs,3,j1,j2)
        if (outs(j1:j2) == '1') then
          if (havers) then
            call prsed4a(k,nbask(1),nsp,s_site1,qnu1,nlst,ilst)
            if (allocated(ilst)) deallocate(ilst)
          else
            call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
            goto 10
          endif
        elseif (outs(j1:j2) == '2') then
          if (have2s) then
            if (.not. allocated(ilst)) allocate(ilst(1,1)) ! Avoid compiler complaints
            call prsed4a(k,nbask(2),nsp,s_site2,qnu2,nlst,ilst)
            if (allocated(ilst)) deallocate(ilst)
          else
            call info0(0,0,0,' no 2nd density... cannot do: '//
     .        trim(outs))
            goto 10
          endif
        else
          call info0(0,0,0,' failed to recognize 3rd argument in '//
     .      trim(outs))
          goto 98
        endif
      endif

      lnsave = .true.

C ... Linear combination of densities
      elseif (outs(1:4) == 'add ') then
        if (.not. (havers .and. have2s)) goto 97

        call words(outs,nw)
        if (nw < 3) goto 98
        call word(outs,2,j1,j2)
        im = -1
        if (outs(j1:j2) == 'all') im = 2+1
        if (outs(j1:j2) == 'pnu') im = 2
        if (outs(j1:j2) == 'pnu') im = 1
        if (im == -1) then
            call info0(0,0,0,
     .      '     ... mode "'//outs(j1:j2)//'" not recognized')
            goto 98
        endif
        call word(outs,3,j1,j2)
        if (outs(j1:j2) == 'sync') then
          isync = 0
C         Require 3 words or 5 words
          if (im >= 4 .and. nw /= 3 .or.
     .        im < 4 .and. nw /= 5) goto 98
          nw = 3
        else
          isync = 1
C         Require 4 words or 6 words
          if (im >= 4 .and. nw /= 4 .or.
     .        im < 4 .and. nw /= 6) goto 98
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,2,ix,fac)
          if (j /= 2) then
            call word(outs,4,j,j2)
            call info0(0,0,0,
     .        '     ... Could not parse arguments : '//outs(j1:j2))
            goto 98
          endif
          nw = 4
        endif
        call numsyv(n)
        if (im < 4) then
          do  k = 1, 2
            call word(outs,nw+k,j1,j2)
            call lodsyv('nbas',0,dble(nbask(1)),i)
            if (k == 2) call lodsyv('nbas',0,dble(nbask(2)),i)
C           call shosyv(0,0,0,6)
            call mkils0(outs(j1:j2),nlst(k),j)
            if (nlst(k) <= 0) then
              call info0(0,0,0,
     .          '%6p... Bad or null list : '//outs(j1:j2))
              goto 98
            endif
            if (k == 1) then
              allocate(ilst(nlst(k),2))
              call mkilst(outs(j1:j2),nlst(k),ilst)
              call imxmn(nlst(k),ilst,1,l,m)
            elseif (k == 2) then
              if (nlst(2) == 1) then
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
                call ivset(ilst(1,2),2,nlst(1),ilst(1,2))
              elseif (nlst(2) /= nlst(1)) then
                call info2(0,0,0,'%5p... 2nd list: '//outs(j1:j2)//
     .            ' has %i entries when %i are needed',nlst(2),
     .            nlst(1))
                goto 98
              else
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
              endif
C             Check that elements in list are within range
              call imxmn(nlst(k),ilst(1,k),1,l,m)
            endif
            if (l < 1 .or.  k == 1 .and. m > nbask(1)
     .        .or.  k == 2 .and. m > nbask(2)) then
              call info0(0,0,0,
     .          '     ... Entry in list out of range: '//outs(j1:j2))
              goto 98
            endif
          enddo
        endif
        call clrsyv(n)

        call prsed1a(im*10+isync,nlst,ilst,nbask,nsp,
     .    s_site1,qnu1,s_site2,qnu2,fac,lfail)
        if (lfail) then
          call info0(0,0,0,' Densities failed to synchronize')
        elseif (isync == 0) then
          call info0(0,0,0,' Densities synchronized')
        elseif (isync == 1) then
          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
          call info0(2,1,0,' New density from rs file:')
          call prsed2a(0,nbask(1),nsp,s_site1,s_spec1,qnu1)
        endif

        deallocate(ilst)
        lnsave = .true.
        goto 10

C ... Save rst file
      elseif (outs(1:5) == 'save ' .or. outs(1:6) == 'savea ') then
        lbin = outs(1:5) == 'save '
C       For ASA, never a binary file
        lbin = .false.
        call word(outs,2,j1,j2)
        if (j2 >= j1) fn = outs(j1:j2)
        if (lbin) then
          if (j2 < j1) fn = 'rst'
          ifi = fopna(fn,-1,4)
        else
          if (j2 < j1) fn = 'rsta'
          ifi = fopna(fn,-1,0)
        endif
        call info0(10,0,-1,
     .    ' rsedit: write restart file "'//trim(fn)//'"')
        k = iorsa(0,s_ctrl,s_site1,s_spec1,s_lat,s_pot,s_bz,fileid,
     .    nbask(1),natk(1),nspeck(1),qnu1,nit,lbin,-ifi)

        lnsave = .false.

C ... abort
      elseif (outs(1:2) == 'a ') then
        call rx0('aborting rs editor ... no file written')

C ... quit
      elseif (outs(1:2) == 'q '. or. outs(1:5) == 'quit ') then
        if (lnsave) then
          print '('' rs file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) == 'y' .or. outs(1:1) == 'q'))
     .      goto 10
        endif
        call rx0('exit rs editor')

C ... help
      elseif (outs == '?') then
        print 310
        print 311
        print 312
        print 313
        print 314
  310   format(
     .    ' Select one of these options:'/
     .  t4,'rs [fn]',t15,
     .    'read density and associated parameters from ASA restart ',
     .    'file.'/t15,'Use "fn" for file name, ',
     .    'or "rst.ext" if "fn" is not supplied')
  311   format(/
     .  t4,'read [fn]',t15,
     .    'read 2nd density from binary restart file.  Use "fn" for ',
     .    'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
C     .  t4,'addspec #',t15,
C     .    'append species # from 2nd density to 1st density.'/t15,
C     .    'If # is zero, append a generic floating orbital species'/
C     .  t4,'addflt #1 #2  posx posy posz'/
C     .    t15,'insert a floating orbital at site #1.'/t15,
C     .    '#2 = species index; it must point to a floating orbital.'//
     .  t4,'show',t15, 'Show summary information about densities'//
     .  t4,'...The following affect the charge or magnetic moment',
     .    ' of a density.'/
     .  t4,'set all|pnu|qnu 1|2  zers|zerq|flip'/
     .  t14,'* The first argument specifies which objects are affected.'
     .    /
     .  t16,'"all" applies to all quantities,'/
     .  t16,'"pnu" to the linearization parameters'/
     .  t16,'"qnu" to the moments'/
     .  t14,'* The second argument refers to the first or second ',
     .    'density'/
     .  t14,'* Third argument = zers: zeros spin part of density:   ',
     .    'n+ - n- = 0'/
     .  t14,'  Third argument = zerq: zeros charge part of density: ',
     .    'n+ + n- = 0'/
     .  t14,'  Third argument = flip: exchange n+, n- ')

  312   format(/
     .  t4,'...The following put a linear combination of the ',
     .    'two densities into the first, or'/
C    .  t4,'add all|pnu|qnu',t22,'fac1 fac2'/
     .  t4,'add all|pnu|qnu',t22,'fac1 fac2 lst1 lst2'/
     .  t14,'* The first argument''s meaning is the same as the ',
     .    'first argument to "set"'/
     .  t14,'* The objects specified will be replaced by'/
     .  t16,'  fac1*[1st type] + fac2*[second type].'/
     .  t16,'"lst2" is a matching site list for the 2nd density, OR'/
     .  t16,'"lst2" is single site (same site added to all sites)'/)

  313   format(
     .  t4,'...The following exchanges the l=0 parts of two site',
     .    ' densities '/
     .  t4,'exch 1|2',t22,'[flip] site1 site2'/
     .    )

  314   format(
     .    t4,'save  [fn]',t15,'saves restart data in ',
     .    'binary restart file (name="rst" unless fn supplied)'/
     .    t4,'savea [fn]',t15,'saves restart data in ',
     .    'ascii restart file (name="rsta" unless fn supplied)'/
     .    t4,'ssave  fn',t15,'saves site data in file fn'/
     .    t4,'q',t15,'to quit the editor'/
     .    t4,'a',t15,'to abort')

      else
        print '(1x,''unrecognized option: '',a)', trim(outs)

      endif
      goto 10

   97 call info0(0,0,0,'%10pBoth densities must be supplied before'
     .  //'invoking this command ')
   98 call info0(0,0,0,' rsedit:  improper usage of '//trim(outs)//
     .  ' ... nothing done')
      goto 10

      end

      subroutine prsed1a(mode,nlst,ilst,nbask,nsp,
     .  s_site1,qnu1,
     .  s_site2,qnu2,
     .  fac,lfail)
C- Add fac*2nd density to density
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 check synchronization
Ci         :1 add the density
Ci         :2 swap moments l=0 part of 2 local densities.  nlst should be 1
Ci         :4 like 2, but exchange spin in the swap
Ci         :10s digit:
Ci         :0 work with true density only
Ci         :1 work with pnu
Ci         :2 work with v0,v1
Ci         :4 work with n1,n2,n0
Ci         :Any combination of the above is allowed
Ci         :8 same as 7, but exclude core density
Ci   nlst  :number of site densities to which 2nd density is added
Ci   ilst  :list of indices to site densities
Ci         :ilst(:,1) is list of sites to be modified
Ci         :ilst(:,2) is list from which 2nd density is to be found
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci
Ci   site1 :struct for site-specific information; see routine usite
Ci   qnu1  :moments of 1st density
Ci   site2 :struct for site-specific information; see routine usite
Ci   qnu2  :moments of 2nd density
Ci   fac   :Add fac(1) * 1st density + fac(2) * 2nd density into density
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Jul 11
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlst,ilst(nlst,2),nbask(2),nsp,n0
      parameter (n0=10)
      double precision fac(2)
      double precision qnu1(3,n0,nsp,*),qnu2(3,n0,nsp,*)
      logical lfail
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site1(*),s_site2(*)
C ... Local parameters
      integer il,ib,jb,is,js
      integer nglob,stdo,ipr,isync,im
      integer ,allocatable:: flgspc(:)
      double precision pnu(n0,2),pnu2(n0,2)
      logical lfails

      stdo = nglob('stdo')
      call getpr(ipr)
      isync = mod(mode,10)
      im  = mod(mode/10,10)

      lfails = .false.
      if (isync >= 2) then
        call info5(0,1,0,' ... Swap '//
     .    '%?#n==4#and spin flip ##'//
     .    '%?#(n==2|n==3)# qnu##'//
     .    '%?#n==1# pnu##',
     .    isync,im,mod(im,2),0,0)
        if (isync == 3 .or. isync > 4 .or. nlst > 1) then
          call info0(0,1,0,' ... illegal parameters (aborting)')
          return
        endif
      elseif (isync == 1) then
        call info5(0,1,0,' ... Replace '//
     .    '%?#(n==2|n==3)# qnu##'//
     .    '%?#(n==1)# pnu##'//
     .    '  by linear combination'//
     .    ' %d*[1st] + %d*[2nd]',
     .    im,mod(im,2),fac(1),fac(2),0)
      endif
      call info2(0,0,0,' site 1 : %n:1,3i',nlst,ilst(1,1))
      call info2(0,0,0,' site 2 : %n:1,3i',nlst,ilst(1,2))

      lfail = .false.
      allocate(flgspc(nbask(1)))
      call iinit(flgspc,nbask(1))
      do  il = 1, nlst

        ib = ilst(il,1)
        if (ib > nbask(1))
     .    call rxi('illegal site index for density:',ib)
        jb = ilst(il,2)
        if (jb > nbask(2))
     .    call rxi('illegal index for 2nd density:',jb)

        is = s_site1(ib)%spec
        pnu = s_site1(ib)%pnu
        js = s_site2(jb)%spec
        pnu2 = s_site2(jb)%pnu

        if (isync == 0) cycle

C       True local density: modified for im=2,3
        if (im >= 2) then
          if (isync < 2) then
            call lcrho(3*n0,nsp,1,1,fac(1),fac(2),
     .        qnu1(1,1,1,ib),qnu2(1,1,1,jb))
          else
            call swrho(isync/4,3*n0,nsp,1,1,1,
     .        qnu1(1,1,1,ib),qnu2(1,1,1,jb))
          endif
        endif
C       Pnu modified if 1's bit im set
        if (mod(im,2) >= 1) then
          if (isync < 2) then
            call lcrho(n0,nsp,1,1,fac(1),fac(2),pnu,pnu2)
            s_site1(ib)%pnu = pnu
          else
            call rx('not ready for this swap')
          endif
        endif

      enddo

C      if (lfail) then
C        call info0(0,0,0,' ... Radial meshes failed to synchronize')
C      elseif (isync == 0) then
C        call info0(0,0,0,' ... Radial meshes synchronized')
C      elseif (isync == 1) then
C        call info0(0,0,0,' ... Linear combination successful')
C      endif

      lfail = lfail .or. lfails

      end

      subroutine prsed2a(mode,nbas,nsp,s_site,s_spec,qnu)
C- Print out density information
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 print out sphere charges and and magnetic moment
Ci         :1 print out sphere P,Q
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pnu
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   28 Jun 08
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nsp,n0
      parameter (n0=10)
      double precision qnu(3,n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer nglob,nl,ib,il,isp,is
      double precision pnu(n0,2),qm(2),dsum

      nl = nglob('nl')
      qm(2) = 0
      if (mode == 1) print 332
  332 format('  ib   l',7x,'P',25x,'Q')
      do  ib = 1, nbas
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu

        do isp = 1, nsp
          if (mode == 1) then
          do il = 1, nl
            print 333, ib, il-1, pnu(il,isp),
     .        qnu(1,il,isp,ib),qnu(2,il,isp,ib),qnu(3,il,isp,ib)
  333       format(2i4,f12.6,2x,3f12.6)
          enddo
          endif
          qm(isp) = dsum(nl,qnu(1,1,isp,ib),3)
        enddo

        if (nsp == 1) then
          print 334, ib, qm(1)
        elseif (nsp == 2) then
          print 334, ib, qm(1)+qm(2), qm(1)-qm(2)
  334     format(' site',i4,':  q=',f12.6:'  mom=',f12.6)
        endif
      enddo

      end

      subroutine prsed4a(mode,nbas,nsp,s_site,qnu,nlst,ilst)
C- Remove spin or charge part of ASA moments
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 remove spin part, 1 remove charge part, 2, flip spins
Ci         :10s digit:
Ci         :1 work with pnu
Ci         :2 work with qnu
Ci         :Any combination is allowed
Ci   nbas  :number of atoms
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Jul 11
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nsp,n0,nlst,ilst(nlst)
      parameter (n0=10)
      double precision qnu(3,n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,kb,is,ipr,stdo,im,mode0,nlist
      double precision pnul(n0,2),xx,xxv(1)
      procedure(integer) :: nglob

      stdo = nglob('stdo')
      call getpr(ipr)
      mode0 = mod(mode,10)
      im  = mod(mode/10,10)
      if (im > 7) call rx('prsed4 not ready for this mode')
      if (mode0 == 0 .and. nsp == 1) return
      if (mode0 == 2 .and. nsp == 1) return

C      if (nsp == 1) write(stdo,200)
C      if (nsp == 2) write(stdo,201)
C  200 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  6x,'q2')
C  201 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  5x,'mtrue',5x,'q2        m2')

      nlist = nlst
      if (nlist == 0) nlist = nbas
      do  kb = 1, nlist

        ib = kb
        if (nlst /= 0) ib = ilst(kb)

        is = s_site(ib)%spec
        pnul = s_site(ib)%pnu

C       qnu :  modified if im>=2
        if (im >= 2) then
          call splrho(30,nsp,3*n0,1,qnu(1,1,1,ib),xxv,xxv)
          call prsed5(mode0,qnu(1,1,1,ib),3*n0,1,nsp)
          call splrho(31,nsp,3*n0,1,qnu(1,1,1,ib),xx,xxv)
        endif

C       pnu v0,v1: modified if 1's bit im set
        if (mod(im,2) >= 1) then
C         call prmx('pnu',pnul,n0,n0,2)
          call splrho(30,nsp,n0,1,pnul,xxv,xxv)
          call prsed5(mode0,pnul,n0,1,nsp)
          call splrho(31,nsp,n0,1,pnul,xxv,xxv)
          s_site(ib)%pnu = pnul
C         call prmx('pnu',pnul,n0,n0,2)
        endif

      enddo

      call info0(2,1,0,' New density:')
      call prsed2a(0,nbas,nsp,s_site,xx,qnu)
      end

      subroutine prsed5(mode,rho,nr,nlml,nsp)
C- Zero out spin part (mode=0) or charge part (mode=1) of rho
      implicit none
      integer mode,nr,nlml,nsp
      double precision rho(nr,nlml,nsp)

      if (mode == 2) then
        call dscal(nr*nlml,-1d0,rho(1,1,2),1)
        return
      endif
      if (mode == 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
      call dpzero(rho,nr*nlml)
      if (mode == 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
      end

      subroutine splrho(mode,nsp,nr,nlml,rho1,rho2,rhoc)
C- Overwrite spin pol local rho+,rho- with rho,rho+ - rho-, or reverse
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 input (rho+,rho-) -> (rho+ + rho-, rho+ - rho-)
Ci         :1 input (rho+ + rho-, rho+ - rho-) -> (rho+,rho-)
Ci         :10s digit
Ci         :1 suppress splitting of rho2
Ci         :2 suppress splitting of rhoc
Ci         :3 suppress both
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Ci   rhoc  :core density
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer mode,nsp,nr,nlml
      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rhoc(nr,nsp)
      double precision fac

      if (nsp == 1) return
      fac = 1
      if (mod(mode,10) /= 0) fac = .5d0

      call dsumdf(nr*nlml,fac,rho1,0,1,rho1(1,1,2),0,1)
      if (mod(mod(mode/10,10),2) == 0)
     .  call dsumdf(nr*nlml,fac,rho2,0,1,rho2(1,1,2),0,1)
      if (mod(mod(mode/10,10)/2,2) == 0)
     .  call dsumdf(nr,fac,rhoc,0,1,rhoc(1,2),0,1)
      end

      subroutine lcrho(nr,nsp,nlml1,nlml2,fac1,fac2,rho1,rho2)
C- First density is overwritten by linear combination of two densities
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml1 :number of L channels for first density
Ci   nlml2 :number of L channels for second density
Ci   fac1  :scales rho1
Ci   fac2  :scales rho2
Ci   rho2  :second density, dimensioned (nr,nlm2,nsp)
Cio Inputs/Outputs
Ci   rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
Cio        :    rho1 <- fac1 * rho1 + fac2 * rho2
Cio        : fac1 scales rho1 in all nlm1 channels
Cio        : fac2 * rho2 is added into rho1 for min(nlm1,nlm2) channels
Cl Local variables
Cl         :
Cr Remarks
Cr   rho1 <- fac1 * rho1 + fac2 * rho2
Cu Updates
Cu   01 Jul 08 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,nlml1,nlml2
      double precision fac1,fac2
      double precision rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
C ... Local parameters
      integer isp,nlml

      if (nlml1 <= 0) return
      nlml = min(nlml1,nlml2)
      do  isp = 1, nsp
        if (fac1 /= 1) then
          call dscal(nr*nlml1,fac1,rho1(1,1,isp),1)
        endif
        if (fac2 /= 0) then
          call daxpy(nr*nlml,fac2,rho2(1,1,isp),1,rho1(1,1,isp),1)
        endif
      enddo
      end

      subroutine swrho(mode,nr,nsp,nlml,nlml1,nlml2,rho1,rho2)
C- Swap two local densities, possibly with spin flip
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :1 swap spin up 1st density w/ spin down second
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :number of L channels to swap
Ci   nlml1 :number of L channels for first density
Ci   nlml2 :number of L channels for second density
Ci   rho2  :second density, dimensioned (nr,nlm2,nsp)
Cio Inputs/Outputs
Cio  rho1  :first density, dimensioned (nr,nlm1,nsp).  On output
Cio  rho2  :second density, dimensioned (nr,nlm2,nsp).  On output
Cio        :densities are exchanged, and spins possibly swapped.
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   19 Jul 08 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nr,nsp,nlml,nlml1,nlml2
      double precision rho1(nr,nlml1,nsp),rho2(nr,nlml2,nsp)
C ... Local parameters
      integer isp,jsp,nlmll

      nlmll = min(nlml,min(nlml1,nlml2))
      if (nlmll <= 0) return

      do  isp = 1, nsp
        jsp = isp
        if (mode == 1 .and. nsp == 2) jsp = 3-isp
        call dswap(nr*nlmll,rho2(1,1,isp),1,rho1(1,1,jsp),1)
      enddo
      end
