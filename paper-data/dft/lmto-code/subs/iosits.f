      integer function iosits(lio,vn,errh,filnam,ifi,slabl,alat,plat,
     .  nbas,nspec,s_spec,s_site)
C- File I/O of site data into ssepc and ssite strux
C ----------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  (not used)
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos vel eula pl relax vshft ndelta delta bxc
Ci                 cpawt omg omgn domg dmat gc gcu gcorr gii sfvrtx j0
Ci                 pdos rho1 rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl
Ci                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Ci                 taukk pihh pihk pikk sohh sohk sokk sighhx sighkx
Ci                 sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet
Ci                 v0 v1
Co     Stored:     spec pos vel eula pl relax vshft bxc cpawt omg omgn
Co                 domg dmat gc gcu gcorr gii sfvrtx j0 pdos rho1 rho2
Co                 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl
Co                 eqkkl sighh sighk sigkk tauhh tauhk taukk pihh pihk
Co                 pikk sohh sohk sokk sighhx sighkx sigkkx tauhhx
Co                 tauhkx taukkx pihhx pihkx pikkx thet v0 v1
Ci Inputs
Ci   lio  :controls file handling and header I/O
Ci         1s digit: file handling
Ci           0 for read
Ci           1 for write
Ci           2 add to use given ifi instead of obtaining it from fnam
Ci           4 add to use s_spec%name instead of slabl
Ci         10s-100s digit: (file read) bits set conditions for read
Ci           1    file alat must match passed array
Ci           2    file plat must match passed array
Ci           4    file nbas must match passed value
Ci           8    Must read nbas data elements from file
Ci          16    version-id must match file
Ci         10s-100s digit: (file write)
Ci           1    write site positions as multiples of plat
Ci        1000s to 10000s digit: bits indicate what info is read in
Ci           1    I/O alat
Ci                ... alat is optional, if sought; no value of alat
Ci                    is assigned if missing (standard format)
Ci           2    I/O plat required, if sought
Ci           4    I/O nbas required, if sought
Ci           8    I/O species, pos information
Ci                    If bit is not set, pos array is untouched
Ci                ... these data are optional, if sought, that is
Ci                    on file read, arrays are assigned to zero if
Ci                    missing from file data.
Ci                (version 3)
Ci          16    I   return local io as functional value
Ci          32    I/O vel,eula,PL,rlx,vshft (standard format only)
Ci                If data sought, but not available, parameters
Ci                are assigned to zero (rlx set to 111)
Ci                If data not sought, parameters are not assigned
Ci          64    I/O (empirical TB only)
Ci                The following conventions apply to version 2 only:
Ci          16    I/O vel,eula,PL,rlx,vshft (standard format only)
Ci          32    I/O (empirical TB only)
Ci
Ci   errh :controls error handling when data is missing or a required
Ci         match (described in lio, above) is not satisfied.
Ci         0 flags iosits to abort
Ci        >0 iosits returns <0 and prints a warning if verbosity > errh
Ci
Ci  filnam :file name.  See description of argument ifi.
Ci
Ci   vn    : version id.
Ci   slabl : Species labels
Ci
Cio  alat  :Scaling of lattice and position vectors
Cio  plat  :Lattice vectors (dimensionless)
Ci
Cio  ifi   :File logical unit.
Cio       *If 1s digit of lio contains 2's bit,
Cio        iosits uses ifi as the file handle and fnam is not used
Cio       *Otherwise, the file 'filnam' is opened and assigned to ifi.
Cio nbas:   size of basis
Cio sspec
Cl Local variables
Cl  lxpos  :T, site positions expressed as multiples of plat
Cr Remarks
Cr  *File consists of 1st record specifying version number
Cr     % version-id ...
Cr   The file version-id consists of a floating-point number,
Cr   followed by an optional string which marks file types different
Cr   from the standard. See below for possible formats.
Cr   A check is made that caller version no matches file.
Cr
Cr  *Standard format
Cr   First line has nbas and plat data, in this format:
Cr     % version-id [options] io=# nbas=# plat=# # ..
Cr   where the tokens specify:
Cr     version-id is the current version number (2.0)
Cr     io=  indicates file's contents (1000-10000s digits of lio)
Cr     nbas= is the size of the basis
Cr     plat= are the lattice vectors
Cr     ... and there are the following optional tokens:
Cr     fast => data is not parsed through preprocessor and
Cr             may be read using fortran read (no algebra)
Cr     xpos => site positions are expressed as multiples of plat
Cr             as opposed to Cartesian coordinates
Cr
Cr   Following the header line is one single line for each site.
Cr   Its structure is (version 2.0)
Cr     spid   x y z   vx vy vz   euler-angles  PL  rlx(xyz)
Cr   Example site file:
Cr     % vn=3.0 fast io=14 nbas=1 plat=-.5 .5 .5 .5 -.5 .5 .5 .5 -.5
Cr      FE    0 0 1
Cr
Cr  *Kotani format
Cr   Any line beginning with '#' is ignored (comment line)
Cr   First line has
Cr     % [#:]kotani   where `#' is version # (2.0)
Cr   If no `#' is specified, no check is made on version number
Cr   Then follows these lines:
Cr       alat                     (lattice constant, in a.u.)
Cr       plat(1:3,1)              (first lattice vector)
Cr       plat(1:3,2)              (second lattice vector)
Cr       plat(1:3,3)              (third lattice vector)
Cr   Next follow site data, in this format, one for each site
Cr       ibas, iclass, species, pos(ibas)
Cr   Example: CaO
Cr      % vn=kotani
Cr      10.26                ! alat
Cr       0   .5   .5         ! plat(1:3,1)
Cr      .5   .0   .5         ! plat(1:3,2)
Cr      .5   .5   .0         ! plat(1:3,3)
Cr      1 1 Ca  0.0 0.0 0.0  ! ibas, iclass, species, pos
Cr      2 2 Ca  0.3 0.3 0.3  ! etc
Cr      3 3 O   0.5 0.5 0.5  !
Cr      4 3 O   1.0 1.0 1.0  !
Cr   NB: iosits ignores 'class' column.
Cb Bugs
Cb   No check is made on version number, Kotani style input.
Cb   s_spec contains slabl information; should not be passed twice
Cm MPI
Cm   Master does I/O. If read then nbas,alat,plat,ssite are broadcast.
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   30 Oct 07 New version 3, compatible with version 3 iosite.f
Cu   20 Apr 07 Bug fix, blank lines at EOF
Cu   14 Apr 03 MPI enabled
Cu   24 May 02 standard format will read alat (optional)
Cu   11 Jan 02 xpos option implemented
Cu   11 Jan 02 iosits can read input using a2vec
C ----------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
      integer lio,errh,nbas,nspec,ifi
      character*8 slabl(*)
      character*(*) filnam
      double precision alat,plat(3,3),vn
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      logical mlog,cmdopt
      integer ierr
C ... These are for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),aa*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
C Local variables
      character locnam*8, spid*8, readin*6
      double precision alatx,posl(3),vell(3),eulal(3),plx(9),ddot,
     .  vshftl,delta(6),vnx,plati(3,3),xx(1)
      logical bittst,ltmp,parstr,a2bin,lfast,lxpos,rdstrn
      integer a2vec,ipr,lio0,lio12,lio345,llio345,fopnx,ib,is,irlx(3),
     .  ipll,lgunit,i,j,nbasl,ix(9),ndel,iirlx,j1,j2,getdig,bitand,isw
      integer procid,master
      data procid /0/ master /0/
      data rdarg /'#{}% c'/

C ... Setup
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      mlog = cmdopt('--mlog',6,0,spid)
      iosits = -1
      call getpr(ipr)
      locnam = filnam
      lio0   = mod(lio,10)
      lio12 = mod(lio/10,100)
      lio345 = mod(lio/1000,1000)
C     Remake convention for lio345 if input is earlier than v2
      if (vn < 3) then
        lio345 = mod(lio345,16) + 32*(lio345/16)
      endif
      if (procid == master) then
        if (.not. bittst(lio0,2)) ifi = fopnx(filnam,100,16+8+0,-1)
        rewind ifi
      endif

C --- File write ---
      if (mod(lio0,2) == 1 .and. procid == master) then

        lxpos = bittst(lio12,1)
        if (lxpos) call dinv33(plat,0,plati,alatx)

C   ... Write header
        if (alat > 0) then
          call awrit6('%% site-data vn=%,1d %?!n!xpos !!fast io=%i'//
     .      ' nbas=%i alat=%;8,1d plat=%9:1;8,1d',a,len(a),ifi,vn,
     .      isw(lxpos),lio345,nbas,alat,plat)
        else
          call awrit4(
     .    '%% site-data vn=2.0 %?!n!xpos !!fast io=%i nbas=%i '//
     .      'plat=%9:1;8,1d',a,len(a),ifi,isw(lxpos),lio345,nbas,plat)
        endif
        if (lio345 >= 32) write (ifi,1)
        if (lio345 < 32) write (ifi,2)
    1   format('#',t26,'pos',t64,'vel',t103,'eula',t126,'vshft',t133,'PL rlx')
    2   format('#',t26,'pos')
        do  ib = 1, nbas
          is = s_site(ib)%spec
          posl = s_site(ib)%pos
          vell = s_site(ib)%vel
          eulal = s_site(ib)%eula
          ipll = s_site(ib)%pl
          irlx = s_site(ib)%relax
          vshftl = s_site(ib)%vshft
          if (lxpos) then
            call dcopy(3,posl,1,plx,1)
C           Forward: pos+ = plat posp+
C           call dgemm('N','N',3,1,3,1d0,plat,3,plx,3,0d0,posl,3)
C           Reverse:  posp+ = (plat)^-1 pos+
            call dgemm('N','N',3,1,3,1d0,plati,3,plx,3,0d0,posl,3)
          endif

          if (lio0/4 > 0) then
            spid = s_spec(is)%name
          else
            spid = slabl(is)
          endif
          a = ' '//spid
C          call awrit6('%9p%3:1;10F %3:1;8F %3:1;8F  %;8F%,3i %3,1i',
C     .      a,len(a),-ifi,posl,vell,eulal,vshftl,ipll,irlx)
          if (lio345 >= 32) then
           call awrit6('%8p%3;12,7D%3:1;12,7D%3:1;12,7D %;8F%,3i %3,1i',
     .      a,len(a),-ifi,posl,vell,eulal,vshftl,ipll,irlx)
          else
            call awrit1('%8p%3;12,7D',a,len(a),-ifi,posl)
          endif
C         Write tbe data
          if (bittst(lio345,64)) then
            call dpzero(delta,6)
            ndel = s_site(ib)%ndelta
            delta = s_site(ib)%delta
            call awrit1(' tbe:%9p%6:1;9F',' ',80,ifi,delta)
          endif
        enddo

        call info2(30,1,0,' iosits: wrote to file '''//
     .    trim(locnam)//''', %i sites',nbas,0)

        iosits = 0
C --- File read ---
      elseif (procid == master) then
        ib = 0
        readin = 'header'
        call getpr(ipr)
        nr = 0
        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
        if (nr == 0) goto 99

C   ... Check for alternative formats
        j1 = 0
        if (.not. parstr(a,'vn=',len(a)-3,3,'=',j1,j2)) goto 98
        if (a(j1:j1) /= ' ') goto 40
        call nword(a,1,j1,j2)
        i = 0
        if (.not. parstr(a(j1+3:),'kotani',j2-j1-6,6,'=',i,j)) goto 40

C   --- File read, Kotani's format : no algebraic expressions ---
        llio345 = 15
        read(ifi,*) alatx
        if (bittst(lio345,1)) alat = alatx
        if (bittst(lio12,1)) then
          call fsanrg(alatx,alat,alat,1d-6,'iosits:','alat',.true.)
        endif

C   ... Read plat
        read(ifi,*,err=98,end=98) plx(1),plx(2),plx(3)
        read(ifi,*,err=98,end=98) plx(4),plx(5),plx(6)
        read(ifi,*,err=98,end=98) plx(7),plx(8),plx(9)
        if (bittst(lio345,2)) call dcopy(9,plx,1,plat,1)
        if (bittst(lio12,2)) then
          call daxpy(9,-1d0,plat,1,plx,1)
          if (ddot(9,plx,1,plx,1) > 1d-14) goto 99
        endif

C   ... Read nbas.  Find indirectly by counting lines.
C       First char = `#' => comment line
        if (bittst(lio345,4)) then
C         i = 0
          nbasl = 0
   31     continue
          if (rdstrn(ifi,a,recl,.false.)) then
C           Skip comment lines
            if (a(1:1) == '#') goto 31
C           Formatted read ibas as a float, for portability
            call word(a,1,j1,j2)
            read(a(j1:j2),'(g72.0)',err=98) alatx
C           Error unless ibas is an integer
            ib = alatx
            if (dble(ib)-alatx /= 0) goto 98
C           Keep track of largest value of ib
C           i = max(i,ib)
            nbasl = nbasl+1
C           For now, sites must be ordered
            if (ib /= nbasl) goto 98
            goto 31
          endif
          if (bittst(lio345,4)) nbas = nbasl
          if (bittst(lio12,4))
     .      call sanrg(.true.,nbasl,nbas,nbas,'iosits:','file''s nbas')
C         file pointer now at EOF.  Restore position
          rewind ifi
          nr = 0
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .      ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          do  i = 1, 4
            ltmp = rdstrn(ifi,a,len(a),.false.)
          enddo
        endif

C  ...  Read site data, Kotani format
        if (lio345 >= 8) then
        do  ib = 1, nbas
          if (.not. rdstrn(ifi,a,len(a),.false.)) goto 98
C         First and second columns (site and class index) ignored
          call word(a,3,j1,j2)
C         Column 3 is species index
          spid = ' '
          read(a(j1:j2),'(a)',err=98,end=98) spid
          do  i = 1, nspec
            is = i
            if (lio0/4 > 0) then
              if (spid == s_spec(i)%name) exit
            else
              if (spid == slabl(i)) exit
            endif
            if (i == nspec) goto 98
          enddo
C          call tokmat(spid,slabl,nspec,8,' ',is,j,.false.)
C          is = is+1
CC         Error if species not found in list
C          if (is <= 0) goto 98
C         Columns 4-6 are site position
          read(a(j2+1:),*,err=98,end=98) posl
          call dpzero(vell,3)
          call dpzero(eulal,3)
          call ivset(irlx,1,3,1)
          ipll = 0
          vshftl = 0
          s_site(ib)%spec = is
          s_site(ib)%pos = posl
          s_site(ib)%vel = vell
          s_site(ib)%eula = eulal
          s_site(ib)%pl = ipll
          s_site(ib)%relax = irlx
          s_site(ib)%vshft = vshftl
        enddo
        endif
        goto 999

C   ... Read, std format.  Return here to resume parsing for arguments
C   ... Read version and check compatibility with calling program
   40   continue
        i = j1+2
        if (.not. a2bin(a,vnx,4,0,' ',i,len(a)-5)) goto 99
        call sanrg(.true.,int(vnx),2,3,'IOSITS:','file''s version number')
C   ... Read lio345
        i = 0
        readin = ' "io" '
        if (.not. parstr(a,'io=',len(a)-3,3,' ',i,j)) goto 99
        i = j-1
        if (.not. a2bin(a,llio345,2,0,' ',i,len(a)-3)) goto 99
C       Remake convention for llio345 if file is v2
        if (vnx < 3) then
          llio345 = mod(llio345,16) + 32*(llio345/16)
        endif
C   ... Read nbas
        if (bittst(lio12,4) .or. bittst(lio345,4)) then
          i = 0
          readin = '"nbas" '
          if (.not. parstr(a,'nbas=',len(a)-5,5,' ',i,j)) goto 99
          i = j-1
          if (.not. a2bin(a,nbasl,2,0,' ',i,len(a)-5)) goto 99
          if (bittst(lio345,4)) nbas = nbasl
          if (bittst(lio12,4))
     .      call sanrg(.true.,nbasl,nbas,nbas,'iosits:','file''s nbas')
        endif
C   ... Read alat (optional if sought), plat (reqd if sought)
        if (bittst(lio12,1) .or. bittst(lio345,1)) then
          i = 0
          readin = '"alat" '
          if (parstr(a,'alat=',len(a)-5,5,' ',i,j)) then
          i = j-1
          if (a2vec(a,len(a)-5,i,4,', ',2,-3,1,ix,plx) == 1) then
            if (bittst(lio345,1)) alat = plx(1)
          endif
          endif
        endif
        if (bittst(lio12,2) .or. bittst(lio345,2)) then
          i = 0
          readin = '"plat" '
          if (.not. parstr(a,'plat=',len(a)-5,5,' ',i,j)) goto 99
          i = j-1
          if (a2vec(a,len(a)-5,i,4,', ',2,-3,9,ix,plx) /= 9) goto 99
          if (bittst(lio345,2)) call dcopy(9,plx,1,plat,1)
          if (bittst(lio12,2)) then
            call daxpy(9,-1d0,plat,1,plx,1)
            if (ddot(9,plx,1,plx,1) > 1d-14) then
              call info(0,0,0,
     .          ' IOSITS: file plat does not match passed data',0,0)
              call rxs('IOSITS: file mismatch, file ',filnam)
            endif
          endif
        endif
C       See whether input file can be read with fortran read
        i = 0
        lfast = parstr(a,'fast ',len(a)-5,5,' ',i,j)
        i = 0
        lxpos = parstr(a,'xpos ',len(a)-5,5,' ',i,j)

        readin = 'data'
        if (lio345-bitand(lio345,16) >= 8) then
C         Exit if file doesn't contain enough data
          if (mod(lio345,16) > mod(llio345,16)) goto 99
    4     continue
          a(1:1) = ' '
          aa(1:1) = ' '
          if (lfast) then
    5       continue
            ltmp = rdstrn(ifi,a,len(a),.false.)
            nr = nr+1
            if (ltmp .and. a(1:1) == '#' .or. a == ' ' .and.
     .          ib < nbas) goto 5
            if (lio345 >= 64) then
    6         continue
              ltmp = rdstrn(ifi,aa,len(aa),.false.)
              nr = nr+1
              if (ltmp .and. aa(1:1) == '#' .or. aa == ' ') goto 6
            endif
          else
            call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .        ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
            ltmp =  nr > 0
            if (ltmp .and. lio345 >= 64)
     .        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .        ilist,nlist,vnam,ctbl,mxchr,aa,recrd,recl,nr)
            ltmp =  nr > 0
          endif

C         If all sites have been read, exit
          if (ib >= nbas) then
            call info2(40,0,0,' IOSITS: read pos'//
     .        '%?!n>15!,vel,eula,vshft,ipl,irlx!! from file unit %i',
     .        llio345,ifi)
            goto 999
          endif

          if (.not. ltmp) then
            if (ib < nbas .and. bittst(lio12,8)) then
              if (errh <= ipr .or. errh == 0) print *,
     .          'IOSITS: file does not contain',nbas,' sites'
              goto 99
            endif
            goto 999
          endif
C         Read data for this site
          ltmp = .false.
          call word(a,1,j1,j2)
          spid = ' '
          read(a(j1:j2),'(a)',err=98,end=98) spid
          do  i = 1, nspec
            is = i
            if (lio0/4 > 0) then
              if (spid == s_spec(i)%name) exit
            else
              if (spid == slabl(i)) exit
            endif
            if (i == nspec) goto 98
          enddo
C          call tokmat(spid,slabl,nspec,8,' ',is,j,.false.)
C          is = is+1
C          if (is <= 0) goto 98
C         No vel,eula,vshft,ipl,irlx available on disk
          if (llio345 < 32) then
            call dpzero(posl,3)
            if (lio345 >= 32) then  ! Default values if sought but not available
              call dpzero(vell,3)
              call dpzero(eulal,3)
              vshftl = 0
              ipll = 0
              iirlx = 111
            endif
              if (lfast) then
                read(a(j2+1:),*,err=98,end=98) posl
              else
                i = 0
                if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,posl)
     . /= 3) goto 98
              endif
C          vel,eula,vshft,ipl,irlx are available on disk
          elseif (lfast) then
            read (a(j2+1:),*,err=98,end=98) posl,vell,eulal,vshftl,ipll,
     .            iirlx
            else
              i = 0
              if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,posl)
     . /= 3) goto 98
              if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,vell)
     . /= 3) goto 98
              if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,3,ix,eulal)
     . /= 3) goto 98
              if (a2vec(a(j2+1:),len(a)-j2,i,4,', ',2,3,1,ix,vshftl)
     . /= 1) goto 98
              if (a2vec(a(j2+1:),len(a)-j2,i,2,', ',2,3,1,ix,ipll)
     . /= 1) goto 98
              if (a2vec(a(j2+1:),len(a)-j2,i,2,', ',2,3,1,ix,iirlx)
     . /= 1) goto 98
            endif
           if (lxpos) then
             call dcopy(3,posl,1,plx,1)
             call dgemm('N','N',3,1,3,1d0,plat,3,plx,3,0d0,posl,3)
           endif

           ib = ib+1
           s_site(ib)%spec = is
           s_site(ib)%pos = posl

C          Copy vel,eula,vshft,ipl,irlx if sought
           if (lio345 >= 32) then
             do  i = 1, 3
               irlx(i) = getdig(iirlx,3-i,10)
             enddo
             s_site(ib)%vel = vell
             s_site(ib)%eula = eulal
             s_site(ib)%pl = ipll
             s_site(ib)%relax = irlx
             s_site(ib)%vshft = vshftl
           endif

           if (lio345 >= 64) then
             call word(a,1,j1,j2)
             locnam = ' '
             read(a(j1:j2),'(a)',err=98,end=98) locnam
             read(aa(j2+1:),*,err=98,end=98) delta
           endif
C          We can quit now if ib eq nbas but read one more line in case
C          site file has comments at the end of file
C          if (ib < nbas) goto 41
           goto 4
         endif

      endif

      if (mod(lio0,2) == 1) then
        call mpi_bcast(iosits, 1, mpi_integer, master, mpi_comm_world, ierr)
        return
      end if
C ... Exit point for read
  999 continue

      call mpibc1(nbas,1,2,.false.,'iosits','nbas')
      call mpibc1(alat,1,4,.false.,'iosits','alat')
      call mpibc1(plat,9,4,.false.,'iosits','plat')
      if (lio345 >= 8)  then
        call bcast_strx(2**1,xx,xx,xx,xx,xx,xx,xx,s_site,xx,0,nbas)
      endif

      iosits = 0
C     Copy llio345 if sought
      if (bittst(lio345,16)) iosits = llio345
      return

C ... Error handling
   98 if (errh <= ipr .or. errh == 0) then
        call awrit1(' IOSITS line %i: missing or incompatible data:',
     .    ' ',80,lgunit(1),nr)
        a(70:) = ' ...'
        call awrit0('%a',a,len(a),-lgunit(1))
        call rx('failed to read site data')
      endif

   99 continue
C     Copy llio345
      if (bittst(lio345,16)) lio = llio345

      if (errh == 0) then
        print 3,':',readin,filnam
       elseif (errh <= ipr) then
        print 3,'(warning)',readin,filnam
       endif
      if (errh == 0) call rxs('IOSITS: file mismatch, file ',filnam)
    3 format(1x,'IOSITS ',a,' failed to read ',a,' from file ',a)

      end
