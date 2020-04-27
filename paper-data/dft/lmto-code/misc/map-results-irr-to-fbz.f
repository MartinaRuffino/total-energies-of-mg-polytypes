C#define F90
C Expands data from `results' file into full BZ
C Requires files points, allpoints, results
C Syntax:
C   map-results-irr-to-fbz [ -ehf|-ehk|-eavg ] [-s#] [-ft] [-nozero]
C or for Takao's mode3:
C   map-results-irr-to-fbz -takao3[:n] -s-1
C The output is file jr.dat
C
C To generate file 'allpoints', invoke
C   lm|lmf --pr80  ...
C For file 'allpoints' cut table beginning after these lines
C    BZMESH: qp mapping
C  i1..i3                        qp                  iq   ig g
C The table as lines such as
C (1,1,1)          0.000000   0.000000   0.000000     1    1 i*i
C (2,1,1)         -0.125000   0.125000   0.125000     2    1 i*i
C (1,2,1)          0.125000  -0.125000   0.125000     2    2 r3(1,1,-1)
C   .
C   .
C   .
C (5,7,3)         -0.500000   1.000000   0.000000    43   10 r3(-1,1,1)
C and ends with a line such as:
C BZMESH:  43 irreducible QP from 512 ( 8 8 8 )  shift= F F F
C The points table is the next table, following the line
C              Qx          Qy          Qz      Multiplicity    Weight
C and looks like, e.g.
C    1      0.000000    0.000000    0.000000         1        0.003906
C   .
C   .
C   .
C   43      0.000000    0.500000    1.000000         6        0.023438
C
C File 'results' may be generated in different ways.
C Takao's 2005 exchange calculations created jlist files.
C Those can be used as follows:
C   cp jlist results
C   map-results-irr-to-fbz -takao2 -s-1
C
C Takao's 2006 exchange calculations created Jmat files.
C They contain interactions for 1 or more atoms
C   cp Jmat results
C The qp's in 'Jmat' must match the list in 'points'
C Invoke with:
C   map-results-irr-to-fbz -takao3:nat:mom1:mom2... -pos:px1,py1,pz1,px2,py2,pz2,... -ft
C   nat is number of atoms; mom1,mom2,... are the (absolute values) of magnetic moments.
C Example for AFM MnO:
C   map-results-irr-to-fbz -takao3,2,5.36,5.36 -pos:0,0,0,0,-.5,-.5  -ft
C
C   (The -pos switch is needed to scale the J matrix by phase factors)
Cu Updates
Cu   31 Jul 15  Added -magat= switch
C ----------------------------------------------------------------------
      subroutine fmain
      implicit none
      character first*120,dc*1
      logical cmdstr,rdstrn,a2bin
      double precision tol,qp(3),ehfj,ehkj,dq,scl,eout,e0,amom(100),pi,
     .  tpi,sp,spirr,topime
      integer, allocatable:: iqmult(:),ipq(:,:,:)
      real(8), allocatable:: qirr(:,:),qfbz(:,:,:,:),pos(:,:)
      real(8), allocatable:: etot(:,:,:),ehf(:,:,:),ehk(:,:,:),xv(:)
      complex(8), allocatable:: eqij(:,:)
      complex(8), allocatable:: zeq(:,:,:),zetot(:,:,:)
      double precision e(10),bbb,ccc
      integer nirr,fopng,fopna,j1,j2,ix(100),itmp(3),nabc(3),a2vec,iarg,
     .  i,j,k,ifi,n1,n2,n3,nqtot,i1,i2,i3,iter,istart,istop,mode,modew,
     .  pvgfe7,ltakao,nat,lcplx,nblst,iblst(200)
      logical ltmp,lphase,cmdopt
      character outs*120
      equivalence (n1,nabc(1)),(n2,nabc(2)),(n3,nabc(3))
      parameter (tol=1d-5)

      integer wksize
      parameter(wksize= 80 000 000)
      integer w(wksize)
      common /w/ w

      call pshpr(0)
      call wkinit(wksize)
      call wkfast(.true.)
      call poppr

      iarg = 0
      mode = 0
      scl = 1
      modew = 21
      e0 = 0
      ltakao = 0
      nat = 1
      allocate(pos(3,1))
      pos = 0
      lphase = .false.
      pi = 4*datan(1d0)
      tpi = 2*pi
      topime = 0
C     lcplx fixes cast of Jin, Jout
C     1s digit:  0 if J(in) is real, 1 if complex
C     10s digit: 0 if cast of J(out) is determined from Jin
C                1 if Jout is real
C                2 if Jout is complex
      lcplx = 0
      do  j = 1, size(iblst)
        iblst(j) = j
      enddo
   15 continue
      iarg = iarg+1
      if (.not. cmdstr(iarg,first)) goto 30
      if (.false.) then
      else if (first(1:5) == '-help' .or.
     .         first(1:6) == '--help' .or.
     .         first(1:3) == '--h') then
        goto 99
      else if (first(1:7) == '--check') then
      else if (first(1:4) == '-ehf') then
        mode = 0
      else if (first(1:4) == '-ehk') then
        mode = 1
      else if (first(1:4) == '-eav') then
        mode = 2
      else if (first(1:4) == '-ftz') then
        modew = 121
        lcplx = 20
      elseif (first(1:7) == '-magat=')  then
        call mkils0(first(8:),nblst,iblst)
        if (nblst /= nat) call rx('size of list must match nat')
        call mkilst(first(8:),nblst,iblst)
      else if (first(1:3) == '-ft') then
        modew = 121
        lcplx = 10
      else if (first(1:7) == '-takao3') then
        dc = first(8:8)
        if (dc == ' ') call rxs('missing parameters, switch ', first)
        ltakao = 3
        e0 = -1
        k = 8
        j = a2vec(first,len(first),k,2,dc//' ',2,2,1,ix,nat)
        if (j < 0) call rxs('failed to parse switch ',first)
        j = a2vec(first,len(first),k,4,dc//' ',2,2,nat,ix,amom)
        if (j < nat) call rxs('failed to parse switch ',first)
        deallocate(pos)
        allocate(pos(3,nat))
        pos = 0
      else if (first(1:4) == '-pos') then
        dc = first(5:5)
        if (nat <= 1) call
     .    rx('-pos must be specified after multiple-atom spec')
        k = 5
        j = a2vec(first,len(first),k,4,dc//', ',3,3,3*nat,ix,pos)
        if (j < 3*nat) call rxs('failed to parse switch ',first)
        lphase = .true.
      else if (first(1:7) == '-takao2') then
        ltakao = 2
      else if (first(1:6) == '-takao') then
        ltakao = 1
      else if (first(1:7) == '-nozero') then
        e0 = -1
      else if (first(1:2) == '-s') then
        i = 2
        if (.not. a2bin(first,scl,4,0,' ',i,-1)) goto 99
      else
        call rxs('failed to parse switch ',first)
      endif
      goto 15
C ... End of switches
   30 continue
C     At this point, 10s digit of lcplx is 0,1,2
C     Output J is assumed real
      if (lcplx/10 == 2) then
        if (ltakao == 3) lcplx = 21
C     Output J is complex
      elseif (lcplx/10 == 1) then
        if (ltakao == 3) lcplx = 11
C     Else, input and output J takes the same cast
      else
        if (ltakao == 3 .or. lcplx == 1) lcplx = 11
      endif

C --- Read list of irreducible q-points ---
      ifi = fopng('points',-1,1)
      rewind ifi
      do  i = 1, 9999999
        read(ifi,*,END=31,ERR=31) j
        if (j /= i) call rx('bad points file')
      enddo
      call rx('did not find eof!')

   31 continue
      nirr = j
      allocate(qirr(3,nirr))
      allocate(iqmult(nirr))
      rewind ifi
      do  i = 1, nirr
        read(ifi,*) j,qirr(:,i),iqmult(i)
        if (j /= i) call rx('bad points file')
      enddo
      call fclr(' ',ifi)
      call info2(0,0,0,' map-results-irr-to-fbz: read %i points from '//
     .  'file points',nirr,0)

C --- Read data from results file ---
      ifi = fopng('results',-1,1)
      rewind ifi
      if (mod(lcplx,10) > 0) then
        allocate(zeq(nirr,nat,nat))
        allocate(eqij(nat,nat))
        allocate(xv(2*nat*nat))
        zeq = 0
      endif
      allocate(ehf(nirr,nat,nat))
      allocate(ehk(nirr,nat,nat))
      ehf = 0
      ehk = 0
      do  i = 1, 9999999
        if (ltakao == 3) then
C         Deal with bug in read: read eqij into real xv; then copy
          read(ifi,*,END=41,ERR=41) qp, xv(1:2*nat*nat)
C         Exchange integrals, Takao's convention
          call dcopy(2*nat*nat,xv,1,eqij,1)
C         bug: last qp in Takao's 'results' must match last in 'points'
C         Used as a sanity check
          call getirr(nirr,qirr,qp,j)

           if (cmdopt('--check',7,0,outs)) then
           if (i == 7 .or. i == 19 .or. i == 115) then
            bbb  = eqij(1,1) + eqij(2,2)
            ccc  = eqij(1,1)*eqij(2,2) - eqij(1,2)*eqij(2,1)
            e(1) = (-bbb - sqrt(bbb**2 - 4*ccc))/2d0
            print *, i
            print *, qp,e(1),amom(1)*e(1)*13.6d0*1000
            call zprm('h',2,eqij,nat,nat,nat)
C            call zhevx(nat,nat,eqij,eqij,0,.true.,nat,9d9,j1,wk,.true.,
C     .        e,nat,zeq)
          endif
          endif

C         Scale J to Mark's convention: Jij(takao) S_i S_j = Jij(mark) e_i e_j
C         Jij(mark) = Jij(takao) S_i S_j = Jij(takao) (amom(i)/2) (amom(j)/2)
          do  j1 = 1, nat
          do  j2 = 1, nat
            eqij(j1,j2) = eqij(j1,j2)*(amom(j1)/2)*(amom(j2)/2)
          enddo
          enddo
        else if (ltakao == 2) then
          read(ifi,*,END=41,ERR=41) qp, ehfj, amom(1), ehkj
C         bug: last qp Takao's 'results' must match last in 'points'
          call getirr(nirr,qirr,qp,j)
          if (j == -1) then
            call rx1('cannot match qp = %3;8,5D to any qp in points',qp)
          endif
          ehfj = ehkj / 13.6d3 / 4 * amom(1)
          ehkj = ehfj
        else if (ltakao == 1) then
          read(ifi,*,END=41,ERR=41) qp, ehfj
          ehkj = ehfj
          j = i
        else
          read(ifi,*,END=41,ERR=41) qp, j, ehfj, ehkj, iter, dq
          if (j > nirr)
     .      call rxi('results file contains funny point, line',i)
        endif
        do  k = 1, 3
          if (abs(qp(k)-qirr(k,j)) > tol) then
            call rxi('results file contains wrong qp, line',i)
          endif
        enddo

C   ... Copy J to irr part
        if (mod(lcplx,10) > 0) then

C         Scale J by phase factors
C          if (lphase) then
C            do  j1 = 1, nat
C            do  j2 = 1, nat
C              sp = tpi* (qp(1)*(pos(1,j1)-pos(1,j2)) +
C     .                   qp(2)*(pos(2,j1)-pos(2,j2)) +
C     .                   qp(3)*(pos(3,j1)-pos(3,j2)))
CC              print *, j1,j2, eqij(j1,j2)*cdexp(dcmplx(0,-sp))
C              eqij(j1,j2) = eqij(j1,j2)*cdexp(dcmplx(0,-sp))
CC              print *, j1,j2, eqij(j1,j2)
C            enddo
C            enddo
C          endif

          do  j2  = 1, nat
          do  j1  = 1, nat
C           topime = max(topime,abs(dimag(eqij(j1,j2))))
            zeq(j,j1,j2) = eqij(j1,j2)
C           ehf(j,j1,j2) = dble(eqij(j1,j2))
          enddo
          enddo

        else
          ehf(j,1,1) = ehfj
          ehk(j,1,1) = ehkj
        endif

      enddo
   41 continue
      call fclr(' ',ifi)
      if (j /= nirr) then
        call rxi('funny results file ... abort after reading point',j)
      endif

      if (lphase) then
        call info0(10,0,0,' scaled J by phase factor')
      endif
C      if (topime >= 0) then
C        call info2(10,0,0,' max Im part of J = %,6;6d',topime,0)
C      endif
C      if (topime >= 1d-6) then
C        call rx('Im(J) too large ... aborting')
C      endif

C --- Sanity check: count how many nonzero energies were read ---
      istart = 0
      do  i = 1, nirr
        if (mod(lcplx,10) > 0) then
          ltmp = zeq(i,1,1) == 0
        else
          ltmp = ehf(i,1,1) == 0
        endif
C       ltmp = ehf(i,1,1) == 0
        if (ltmp) then
          if (istart == 0) then
            istart = i
          endif
          istop = i
        else
          if (istart /= 0) then
            call info2(10,0,0,
     .        ' missing results for points %i .. %i',istart,istop)
            istart = 0
          endif
        endif
      enddo
      if (istart /= 0) call info2(10,0,0,
     .  ' missing results for points %i .. %i',istart,istop)

C --- Read list q-points in full BZ ---
      ifi = fopng('allpoints',-1,1)
      rewind ifi
      nabc = 0
      do  i = 1, 9999999
        first = ' '
        if (.not. rdstrn(ifi,first,len(first),.false.)) goto 51
        nqtot = i
        call word(first,1,j1,j2)
        k = 1
        j = a2vec(first(j1:),j2-j1,k,2,',)',2,2,3,ix,itmp)
        if (j /= 3) call rx('failed to read file allpoints')
        do  k = 1, 3
          nabc(k) = max(nabc(k),itmp(k))
        enddo
      enddo
   51 continue
      call info2(0,0,0,' read %i points from file'//
     .  ' allpoints, nabc =%3:1i',nqtot,nabc)

C --- Create the mapping ipq(i1,i2,i3) to qfbz ---
      allocate(ipq(n1,n2,n3))
      allocate(qfbz(3,n1,n2,n3))
      rewind ifi
      do  i = 1, nqtot
        first = ' '
        if (.not. rdstrn(ifi,first,len(first),.false.)) goto 52
        j = i
        call word(first,1,j1,j2)
        k = 1
        j = a2vec(first(j1:),j2-j1,k,2,',)',2,2,3,ix,itmp)
        if (j /= 3) call rx('failed to read file allpoints')
        i1 = itmp(1)
        i2 = itmp(2)
        i3 = itmp(3)
        call word(first,2,j1,j2)
C       print *, i, first(j1:j1+40)
        read(first(j1:),*) qfbz(:,i1,i2,i3),ipq(i1,i2,i3)
      enddo
   52 continue
      call fclr(' ',ifi)

C --- Copy data from irr BZ to full BZ ---
      if (mod(lcplx,10) > 0) then
        allocate(zetot(n1,n2,n3))
        allocate(etot(n1,n2,n3))
      else
        allocate(etot(n1,n2,n3))
      endif
C     allocate(etot(n1,n2,n3))

C     debugging
C      print *, 'debuging ... make JmatTest ...'
C      call snot(n1,n2,n3,qfbz,etot)

C     Reset modew if input J is complex
      if (mod(lcplx,10) > 0) then
C       Output J is complex
        if (lcplx/10 == 2) then
          modew = modew + 2000
C       Output J is real
        else
          modew = modew + 1000
        endif
      else
C       Output J is complex
        if (lcplx/10 == 2) then
          modew = modew + 2000
        endif
      endif

C --- Loop over atom pairs in lmgf compatible order ---
      topime = 0
      do  j2 = 1, nat
      do  j1 = 1, nat

C --- Loop over FBZ ---
      i = 0
      do  i1 = 1, n1
      do  i2 = 1, n2
      do  i3 = 1, n3
        j = ipq(i1,i2,i3)
        if (mod(lcplx,10) > 0) then

C         Scale J by phase factors
          if (lphase) then
            qp = qirr(:,j)
            spirr = tpi* (qp(1)*(pos(1,j2)-pos(1,j1)) +
     .                    qp(2)*(pos(2,j2)-pos(2,j1)) +
     .                    qp(3)*(pos(3,j2)-pos(3,j1)))
            qp = qfbz(:,i1,i2,i3)
            sp = tpi* (qp(1)*(pos(1,j2)-pos(1,j1)) +
     .                 qp(2)*(pos(2,j2)-pos(2,j1)) +
     .                 qp(3)*(pos(3,j2)-pos(3,j1)))
C            print *, j1,j2, zeq(j,j1,j2)
C            print *, j1,j2, zeq(j,j1,j2)*cdexp(dcmplx(0,-spirr))
            sp = -sp
            spirr = -spirr
            topime = max(topime,
     .        abs(dimag(zeq(j,j1,j2)*cdexp(dcmplx(0,-spirr)))))
C            print *, j1,j2, zeq(j,j1,j2)*cdexp(dcmplx(0,sp-spirr))
          endif
          zetot(i1,i2,i3) = zeq(j,j1,j2)*cdexp(dcmplx(0,sp-spirr))*scl

        else
          ehfj = ehf(j,j1,j2)
          ehkj = ehk(j,j1,j2)
          if (mode == 0) then
            eout = ehfj
          else if (mode == 1) then
            eout = ehkj
          else if (mode == 2) then
            eout = (ehfj+ehkj)/2
          endif
          if (eout /= 0) i = i+1
          etot(i1,i2,i3) = eout*scl
        endif
      enddo
      enddo
      enddo

C      if (topime /= 0) then
C        call info5(10,0,0,
C     .    ' ib=%i  jb=%i   max Im part of J(rot) = %,6;6d',j1,j2,topime,
C     .    0,0)
C      endif

      if (e0 /= -1) then
        e0 = etot(1,1,1)
        do  i1 = 1, n1
        do  i2 = 1, n2
        do  i3 = 1, n3
          if (etot(i1,i2,i3) /= 0) etot(i1,i2,i3) = etot(i1,i2,i3)-e0
        enddo
        enddo
        enddo
        call info2(0,0,0,' ... copied %i points to full BZ;'//
     .    ' subtracting e0=%d',i,e0)
      endif

      if (scl /= 1 .or. mod(modew,1000) == 121)
     .  call info2(0,0,0,' ... scale by %d '//
     .  '%?;(n==121); ... ft data;',scl,mod(modew,1000))

C --- Map plat to -plat ---
C      if (ltakao == 3) then
C        call fftz3(zetot,n1,n2,n3,n1,n2,n3,1,0,-1)
C        call platsw(n1,n2,n3,zetot)
C        call fftz3(zetot,n1,n2,n3,n1,n2,n3,1,0,1)
C      endif

C --- Write available points to full BZ ---
      call info0(0,0,0,' ... writing file jr.dat')
      ifi = fopna('jr',-1,0)
      if (j1 == 1 .and. j2 == 1) rewind ifi
      call pshpr(41)
      if (mod(lcplx,10) > 0 .and. modew /= 1021) then
        i = pvgfe7(zetot,modew,ifi,iblst(j2),iblst(j1),n1,n2,n3)
      else
        if (modew == 1021) then
          call dcopy(n1*n2*n3,zetot,2,etot,1)
          i = pvgfe7(etot,21,ifi,iblst(j2),iblst(j1),n1,n2,n3)
        else
          i = pvgfe7(etot,modew,ifi,iblst(j2),iblst(j1),n1,n2,n3)
        endif
      endif

      enddo
      enddo

      if (topime /= 0) then
        call info2(10,0,0,' max Im part of J(rot) = %,3;3g',topime,0)
      endif

      call fclr(' ',ifi)


      return

   99 continue
      print 333
  333 format(
     .  ' usage: '/
     .  ' map-results-irr-to-fbz  [-options]',/2x,
     .  ' Requires files points, allpoints, results to execute.'//
     .  ' Options: '/
     .  '   -ehf | -ehk | -eav'/
     .  '   -takao | -takao2 | -takao3:nat,|mom1|,|mom2|,...'/
     .  '   -pos:px1,py1,pz1,px2,py2,pz2,...'/
     .  '   -s#'/
     .  '   -ft'/
     .  '   -magat=ib1,ib2,...'/
     .  )

      end

      integer function pvgfe7(s,lio,ifi,ib,jb,n1,n2,n3)
C- I/O for exchange integrals J
C ----------------------------------------------------------------------
Ci  lio  0  read, return s in real space (see Remarks)
Ci       1  write s in real space
Ci       2  read, find start of array containing (ib,jb) pair
Ci      10  read s, return s in recip space
Ci      11  write real s(qs) from input real s(rs)
Ci      21  write real s(qs) from input real s(qs)
Ci     111  write real s(rs) from input real s(rs)
Ci     121  write real s(rs) from input real s(qs)
Ci     ...  1000s digit =1,2 implies input s is complex
Ci     ...  1000s digit =2   implies output  s is complex
Ci    1121  write real s(rs)    from complex input s(qs)
Ci    2121  write complex s(rs) from complex input s(qs)
Ci  ifi     file logical unit number
Ci  ib,jb   pair for which to read/write J
Ci  n1..3   dimensions s
Cio Inputs/Outputs
Cio   s    array containing exchange integrals
Cr  Remarks
Cr
Cr    File read (1s digit mode 0):
Cr    pvgfe7 assumes that (ib,jb) pairs are ordered by increasing ib,jb,
Cr    with jb=fast index.  If file ib>passed ib, or if
Cr    file jb>passed jb and file ib>passed ib, pvgfe7 returns with -1
Cr    If file end-of-file encountered, pvgfe7 returns with -2
Cr
Cr    File read (1s digit mode 2):
Cr    Looks in file for start of array that contains (ib,jb) pair.
Cr    No assumption is made about file order; no attempt is made
Cr    to read the array if it is found.  pvgfe7 returns 0 if start
Cr    is found, or -2 if end-of-file encountered.
Cu Updates
Cu   05 Mar 01 Added read mode 2
Cu   22 Dec 00 turned into a function call
C ----------------------------------------------------------------------
      implicit none
      integer n1,n2,n3,ib,jb,lio,ifi
      double precision s(n1,n2,n3)
C ... Heap
      integer w(1)
      common /w/ w
      logical a2bin
      integer i,j,ib0,jb0,k,nfbz,oh,lgunit,ip,recl,iprint,lio0
      parameter (recl=80)
      character*(20) fmt, space*2, recrd*(recl), dir(4)*7
C#ifdef F90
      complex(8), allocatable:: zs(:,:,:)
C#endif
      data dir /'rows','cols','rs','qs'/

      lio0 = mod(lio,10)
      pvgfe7 = 0
C ... Write
      if (lio0 == 1) then
        fmt = '(9f15.10)'
        space = 'rs'
        if (mod(lio/100,10) == 0) space = 'qs'
        if (lio == 11 .or. mod(lio,1000) == 121) then
          nfbz = n1*n2*n3
          call defcc(oh,-nfbz)
          if (lio/1000 == 1 .or. lio/1000 == 2) then
            call dcopy(nfbz*2,s,1,w(oh),1)
          else
            call dcopy(nfbz,s,1,w(oh),2)
          endif
          i = -1
          if (mod(lio,1000) == 11) i = 1
          call fftz3(w(oh),n1,n2,n3,n1,n2,n3,1,0,i)
          if (lio/1000 == 2) then
            call dcopy(nfbz*2,w(oh),1,s,1)
          else
            call dcopy(nfbz,w(oh),2,s,1)
          endif
          call rlse(oh)
        endif
        if (iprint() > 40)
     .    call awrit2(' pvgfe7:  writing '//space//
     .    ' J for ib=%i, jb=%i',' ',80,lgunit(1),ib,jb)
        if (lio/1000 == 2) then
C#ifdef F90
          call awrit4('%% rows %i cols %i complex  '//space//
     .      '  ib=%i  jb=%i',' ',80,ifi,n1*n2,n3,ib,jb)
          allocate (zs(n1,n2,n3))
          nfbz = n1*n2*n3
          call dcopy(nfbz*2,s,1,zs,1)
          do  i = 1, n1
          do  j = 1, n2
            write(ifi,fmt) (dble(zs(i,j,k)), k=1,n3)
          enddo
          enddo
          write(ifi,"(1x)")
          do  i = 1, n1
          do  j = 1, n2
            write(ifi,fmt) (dimag(zs(i,j,k)), k=1,n3)
          enddo
          enddo
          deallocate (zs)
C#elseC
C         call rx('this branch pvgfe7 requires F90 compiler')
C#endif
        else
          call awrit4('%% rows %i cols %i real  '//space//
     .      '  ib=%i  jb=%i',' ',80,ifi,n1*n2,n3,ib,jb)
          do  10  i = 1, n1
          do  10  j = 1, n2
   10     write(ifi,fmt) (s(i,j,k), k=1,n3)
        endif
C ... Read
      else
C   ... Entry point for new array to read, modes 0 and 2
   31   continue
        read (ifi, '(a)',END=99) recrd
        space = ' '
        if (recrd(1:1) == '%') then
          ip = 1
          ib0 = 0
          jb0 = 0
   30     call skipbl(recrd,recl,ip)
          if (ip >= recl) goto 32
          k = ip-1
          call tokmat(recrd(ip+1:recl),dir,4,7,' ',i,ip,.false.)
          ip = ip+k
          if (i < 0) then
            if (recrd(ip+1:ip+3) == 'ib=') then
              i = 0
              if (.not. (a2bin(recrd(ip+4:),ib0,2,0,' ',i,recl)))
     .          call rxs('pvgfe7: failed to parse ib in',recrd)
              if (lio0 == 2 .and. ib == ib0 .and. jb == jb0) then
                backspace ifi
                return
              endif
C             ib < ib0 => already past this ib ...
              if (lio0 == 0 .and. ib < ib0) then
                backspace ifi
                pvgfe7 = -1
                return
               call rx('pvgfe7 not ready')
              endif
            elseif (recrd(ip+1:ip+3) == 'jb=') then
              i = 0
              if (.not. (a2bin(recrd(ip+4:),jb0,2,0,' ',i,recl)))
     .          call rxs('pvgfe7: failed to parse jb in',recrd)
C             jb < jb0 and ib < ib0 => already past this pair ...
              if (lio0 == 2 .and. ib == ib0 .and. jb == jb0) then
                backspace ifi
                return
              endif
              if (lio0 == 0 .and. jb < jb0 .and. ib <= ib0) then
                backspace ifi
                pvgfe7 = -1
                return
              endif
            endif
            call skp2bl(recrd,recl,ip)
          elseif (i == 0 .or. i == 1) then
            if (.not. (a2bin(recrd,j,2,0,' ',ip,recl))) return
            if (i == 0 .and. j /= n1*n2 .or. i == 1 .and. j /= n3)
     .        call rx('pvgfe7:  file mismatch')
          elseif (i == 2) then
            space = 'rs'
          elseif (i == 3) then
            space = 'qs'
          endif
        else
          call rx('pvgfe7: file jr not in standard format')
        endif
        goto 30
   99   continue
        pvgfe7 = -2
        return
   32   continue
        if (space == ' ') then
          call rx('pvgfe7: file jr not in standard format')
        endif
        do  20  i = 1, n1
        do  20  j = 1, n2
   20   read(ifi,*) (s(i,j,k), k=1,n3)
        if ((ib0 /= 0 .and. ib /= ib0) .or.
     .      (jb0 /= 0 .and. jb /= jb0)) goto 31
        recrd = '%N pvgfe7: read J('//space//'), returning J(  )'
        if (ib0 /= 0) call awrit1('%a, ib=%i',recrd,80,0,ib)
        if (jb0 /= 0) call awrit1('%a, jb=%i',recrd,80,0,jb)
        if (lio == 0 .and. space == 'qs' .or.
     .    lio /= 0 .and. space == 'rs') then
          nfbz = n1*n2*n3
          call defcc(oh,-nfbz)
          call dcopy(nfbz,s,1,w(oh),2)
          i = -1
          space = 'rs'
          if (lio /= 0) i = 1
          if (lio /= 0) space = 'qs'
          call fftz3(w(oh),n1,n2,n3,n1,n2,n3,1,0,i)

C          call zprm3('jq',2,w(oh),n1,n2,n3)
C          call fftz3(w(oh),n1,n2,n3,n1,n2,n3,1,0,-1)
C          call zprm3('jr',2,w(oh),n1,n2,n3)
C          call fftz3(w(oh),n1,n2,n3,n1,n2,n3,1,0,i)
C          call zprm3('jq',2,w(oh),n1,n2,n3)
C
          call dcopy(nfbz,w(oh),2,s,1)
          call rlse(oh)
        endif
        recrd(36:37) = space
        if (iprint() >= 30) call awrit0(recrd,' ',-80,lgunit(1))
C       call prm3(recrd,0,s,n1,n2,n3)
      endif
      end

      subroutine getirr(nirr,qirr,qp,iq)
C- Find iq in qirr that matches qp
C ----------------------------------------------------------------------
Ci Inputs
Ci   nirr
Ci   qirr
Ci   qp    :k-point
Co Outputs
Co   iq    :point iq that matches qirr.  No match -> iq = -1
Cl Local variables
Cl         :
Cr Remarks
Cr
Cb Bugs
Cb   should call latvec to compare
Cu Updates
Cu   08 Sep 05
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nirr,iq
      double precision qirr(3,nirr),qp(3)
C ... Local parameters
      integer i
      double precision tol
      parameter (tol=1d-5)

      do  i = 1, nirr
        if (abs(qp(1)-qirr(1,i)) < tol) then
        if (abs(qp(2)-qirr(2,i)) < tol) then
        if (abs(qp(3)-qirr(3,i)) < tol) then
          iq = i
          return
        endif
        endif
        endif
      enddo
      iq = -1

      end

      subroutine snot(n1,n2,n3,qfbz,etot)
      implicit none
      integer n1,n2,n3
      double precision qfbz(3,n1,n2,n3),etot(n1,n2,n3)

      integer i1,i2,i3,i,j,k,ifi,fopna,modew,iffi,ifmat
C     NN, MnO lattice
      integer nn
      double precision latvec(3,6)
      double precision exnn,pi,tpi,TdotK,qp(3)
      double complex ej,eikt,img,ez(n1,n2,n3)

C     Lattice translation vectors for A-A nearest neighbors
C      latvec(:,1) = (/ 0.5000d0,  0.0000d0, -0.5000d0 /)
C      latvec(:,2) = (/ 0.0000d0,  0.5000d0, -0.5000d0 /)
C      latvec(:,3) = (/ 0.5000d0, -0.5000d0,  0.0000d0 /)
C      latvec(:,4) = (/-0.5000d0,  0.5000d0,  0.0000d0 /)
C      latvec(:,5) = (/ 0.0000d0, -0.5000d0,  0.5000d0 /)
C      latvec(:,6) = (/-0.5000d0,  0.0000d0,  0.5000d0 /)

C     Actual connecting vectors for A-B nearest neighbors
      latvec(:,1) = (/ 0.0000d0, -0.5000d0, -0.5000d0 /)
      latvec(:,2) = (/-0.5000d0,  0.0000d0, -0.5000d0 /)
      latvec(:,3) = (/ 0.5000d0,  0.5000d0,  0.0000d0 /)
      latvec(:,4) = (/-0.5000d0, -0.5000d0,  0.0000d0 /)
      latvec(:,5) = (/ 0.5000d0,  0.0000d0,  0.5000d0 /)
      latvec(:,6) = (/ 0.0000d0,  0.5000d0,  0.5000d0 /)
C     Lattice translation vectors corresponding to connecting vectors
      latvec(:,1) = (/ 0.00000d0, 0.00000d0, 0.00000d0/)
      latvec(:,2) = (/-0.50000d0, 0.50000d0, 0.00000d0/)
      latvec(:,3) = (/ 0.50000d0, 1.00000d0, 0.50000d0/)
      latvec(:,4) = (/-0.50000d0, 0.00000d0, 0.50000d0/)
      latvec(:,5) = (/ 0.50000d0, 0.50000d0, 1.00000d0/)
      latvec(:,6) = (/ 0.00000d0, 1.00000d0, 1.00000d0/)

      nn = 6
      exnn = 1
      pi = 4*datan(1d0)
      tpi = 2*pi
      img = (0d0,1d0)

C ... Use this block to generate etot in the full BZ ...
      modew = 121
      do  i1 = 1, n1
      do  i2 = 1, n2
      do  i3 = 1, n3
        ej = 0
        do  j = 1, nn
          TdotK = tpi* (qfbz(1,i1,i2,i3)*(latvec(1,j)) +
     .                  qfbz(2,i1,i2,i3)*(latvec(2,j)) +
     .                  qfbz(3,i1,i2,i3)*(latvec(3,j)))
          eikt = dcmplx(dcos(TdotK),dsin(TdotK))
          ej = ej + exnn*eikt
        enddo

C        if (abs(dimag(ej)) > 1d-10) stop 'oops'
C        etot(i1,i2,i3) = dble(ej)
        ez(i1,i2,i3) = ej

      enddo
      enddo
      enddo

      call info0(0,0,0,' ... writing file jr.dat')
      ifi = fopna('jr',-1,0)
      call awrit4('%% rows %i cols %i complex qs'//
     .  '  ib=%i  jb=%i',' ',80,ifi,n1*n2,n3,1,1)
      do  i = 1, n1
      do  j = 1, n2
        write(ifi,"(9f15.10)") (dble(ez(i,j,k)), k=1,n3)
      enddo
      enddo
      write(ifi,"(1x)")
      do  i = 1, n1
      do  j = 1, n2
        write(ifi,"(9f15.10)") (dimag(ez(i,j,k)), k=1,n3)
      enddo
      enddo
C     j = pvgfe7(etot,modew,ifi,1,1,n1,n2,n3)
      call rx0('done')

C ... Use this block to mimic Takao's output, irreducible BZ ...
C     Takes list of q-points from his file qlist
      iffi=30
      ifmat= 1017
      open(iffi,file='qlist')
      open(ifmat, file="JmatTest")
      i = 0
      do
        read(iffi,*,end=1010) qp
        i = i+1
        ej = 0
        do  j = 1, nn
          TdotK = tpi* (qp(1)*(latvec(1,j)) +
     .                  qp(2)*(latvec(2,j)) +
     .                  qp(3)*(latvec(3,j)))
          eikt = dcmplx(dcos(TdotK),dsin(TdotK))
          ej = ej + exnn*eikt




C          print "(i4,3f12.5)", j,
C     .      sngl (latvec(:,j) + (/0d0,-0.5d0,-0.5d0/))

        enddo
C        phase = exp( tpi*img*sum(qp*(/0d0,-0.5d0,-0.5d0/)) )
C        print *, ej
C        print *, phase
C        ej = ej * phase
C        print *, ej*phase
C        if (abs(dimag(ej)) > 1d-10) stop 'oops'
C       One atom
        write(ifmat,"(3d18.10, 3x, 255d18.10)")
     .        qp, ej
C       Two atoms
C        write(ifmat,"(3d18.10, 3x, 255d18.10)")
C     .        qp, 0d0,0d0, ej, dconjg(ej), 0d0,0d0

      enddo
 1010 continue
      print *, 'generated file for',i,' sites'
      call rx0('done')

      end

      subroutine platsw(n1,n2,n3,zetot)
C- Overwrites zetot(i1,i2,i3) with zetot(n1-i1,n2-i2,n3-i3) (adjusted for offsets)
      implicit none
      integer n1,n2,n3
      double complex zetot(n1,n2,n3)
      integer i1,i2,i3,i1m,i2m,i3m
      complex(8), allocatable:: wk(:,:,:)

      allocate(wk(n1,n2,n3))
      call dcopy(n1*n2*n3*2,zetot,1,wk,1)

      call zprm3('starting zetot',2,zetot,n1,n2,n3)

      do  i3 = 1, n3
      i3m = n3+2-i3
      if (i3 == 1) i3m = 1

        do  i2 = 1, n2
        i2m = n2+2-i2
        if (i2 == 1) i2m = 1

          do  i1 = 1, n1
          i1m = n1+2-i1
          if (i1 == 1) i1m = 1

          zetot(i1m,i2m,i3m) = wk(i1,i2,i3)
      enddo
      enddo
      enddo

      call zprm3('ending zetot',2,zetot,n1,n2,n3)

      deallocate(wk)
      end

