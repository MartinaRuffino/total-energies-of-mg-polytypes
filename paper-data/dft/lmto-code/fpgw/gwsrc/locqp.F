      subroutine q2octant(mode,tolq,qin,qout)
C- Adds smallest integer to qin to make qout, so that 0<=qout<1 for each component
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 1 => if qout within tolq of 0, set to 0
Ci   tolq  :precision to which qp is known
Ci   qin   : qp to be shifted
Co Outputs
Co   qout  : result of shift
Cr Remarks
Cr   qin is translated to place it in the first octant of Wigner-Seitz cell.
Cr   qin and qout are in multiples of the lattice vectors.
Cr   qin and qout can occupy the same address space
Cu Updates
Cu   30 Mar 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode
      real(8):: qin(3),qout(3),tolq
C ... Local parameters
      double precision qx(3)

C     Add an integer large enough to render all elements of qx positive
C     components within tolq of boundary will get shifted to origin
      qx   = qin + dble(int(abs(qin))) + 1 + tolq
C     Remove integer multiples; restore tolq shift
      qout = qx - int(qx) - tolq

C     Set tiny elements to zero, to support routines requiring q=0 exactly
      if (mode==1) then
        if (abs(qout(1))<tolq .and. abs(qout(2))<tolq .and. abs(qout(3))<tolq) qout = 0
      endif
      end

C     Test q2octant
C      subroutine fmain
C      real(8) qp(3),qout(3)
C      real(8),parameter:: tolq=1d-7
C
C      qp = [0d0,0d0,0d0]
C      call q2octant(1,tolq,qp,qout)
C      call info5(10,0,0,'qin = %3;11,7D  qout = %3;11,7D diff = %3;12,8D',qp,qout,qout-qp,4,5)
C      qp = [-1.01d0,2.02d0,-4.0000001d0]
C      call q2octant(1,tolq,qp,qout)
C      call info5(10,0,0,'qin = %3;11,7D  qout = %3;11,7D diff = %3;12,8D',qp,qout,qout-qp,4,5)
C      qp = [-.99999d0,1.9999d0,-3.999999d0]
C      call q2octant(1,tolq,qp,qout)
C      call info5(10,0,0,'qin = %3;11,7D  qout = %3;11,7D diff = %3;12,8D',qp,qout,qout-qp,4,5)
C      qp = [-.99999d0,1.999999d0,-3.99999991d0]
C      print *, '... last two check mode'
C      call q2octant(1,tolq,qp,qout)
C      call info5(10,0,0,'qin = %3;11,7D  qout = %3;11,7D diff = %3;12,8D',qp,qout,qout-qp,4,5)
C      qp = [-.99999991d0,1.99999991d0,-3.99999991d0]
C      call q2octant(1,tolq,qp,qout)
C      call info5(10,0,0,'qin = %3;11,7D  qout = %3;11,7D diff = %3;12,8D',qp,qout,qout-qp,4,5)
C      end

      subroutine locqp(mode,nqlst,qpp,plat,qin,tolq,iq,qout)
C- Fast search for a qp connected by a lattice vector to a some member of a list
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 qplst are given in Cartesian coordinates
Ci         :1 qplst are given as multiples of reciprocal vectors to plat
Ci         :10s digit
Ci         :0 Return original q that matched qin in form originally given
Ci         :1 Return shortened form of q in Cartesian coordinates
Ci         :2 Return shortened form of q in multiples of RLV
Ci         :100s digit
Ci         :0 Return iq=0 if no qp found
Ci         :1 Abort with error if qp found
Ci         :1000s digit
Ci         :1 If each element is less than tolq set q to 0
Ci   nqlst :size of q list
Ci   qpp   :Tailored information about the qp list, created by locqpi, which see
Ci   plat  :primitive lattice vectors.  Used if 1s digit mode is zero
Ci   qin   :qp to match up to a lattice translation vector
Ci   tolq  :precision to which qp is known
Co Outputs
Co   iq    :points to entry in list matching q, if one exists
Co         :iq=0 if no entry found
Ci   qout  :qp matching qin with, shifted to first octant in BZ.
Cr Remarks
Cr   locqp provides a fast lookup to locate a qp connected by a lattice vector
Cr   to a member of a list of points.
Cr   Call initialization routine locqpi before calling locqp, which generates qpp.
Cr   qpp should not be altered after it is made
Cu Updates
Cu   30 Mar 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nqlst,iq
      real(8) qpp(nqlst,7),plat(3,3),qin(3),tolq,qout(3)
C ... Local parameters
      integer i
      real(8) qloc(3),qlat(3,3),xx

C     qloc must be in units of the RLV
      if (mod(mode,10) == 0) then
        call dgemm('T','N',3,1,3,1d0,plat,3,qin,3,0d0,qloc,3)
      else
        qloc = qin
      endif
      call q2octant(mod(mode/1000,10),tolq,qloc,qloc)

C     Quick search for first element
      i = 0
      call huntx(qpp(1,3),nqlst,qloc(3),[0],i)
      i = max(i,1)
      do  while (qpp(i,3)+tolq > qloc(3) .and. i > 1)
        i = i-1
      enddo

C     First entry in third column matching q(3)
      do  while (qpp(i,3)+tolq < qloc(3))
        i = i+1
        if (i > nqlst) goto 999
      enddo

C     First entry in second column matching q(2)
      do  while (qpp(i,2)+tolq < qloc(2))
        i = i+1
        if (i > nqlst) goto 999
      enddo

C     First entry in first column matching q(1)
      do  while (qpp(i,1)+tolq < qloc(1))
        i = i+1
        if (i > nqlst) goto 999
      enddo
      if (abs(qloc(1)-qpp(i,1))>tolq) goto 999 ! All three components must match
      if (abs(qloc(2)-qpp(i,2))>tolq) goto 999
      if (abs(qloc(3)-qpp(i,3))>tolq) goto 999

      iq = nint(qpp(i,4))
      if (mod(mode/10,10) == 0) then
        qout(1:3) = qpp(i,5:7)
      elseif (mod(mode/10,10) == 1) then
!       call dinv33(plat,1,qlat,xx)
        call minv33tp(plat,qlat)
        call dgemm('N','T',3,1,3,1d0,qlat,3,qpp(i,1),nqlst,0d0,qout,3)
      else
        qout(1:3) = qpp(i,1:3)
      endif

      return

C     Handle case no qp found
  999 continue
      iq = 0
      if (mod(mode/100,10) == 0) return
C     call rx2('locqp: no qp matching q=%s,(%3;6,6d) within tol=%;3g',qin,tolq)
      call rxf('locqp: no qp matching qlst within tol',tolq)
      end

      subroutine locqpi(mode,nqlst,qplst,tolq,plat,qpp)
C- Setup locqp (for fast search for qp contained in a host list)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 qplst are in Cartesian coordinates
Ci         :1 qplst are in multiples of reciprocal vectors to plat
Ci   nqlst :size of qp list
Ci   qplst :vector of qp
Ci   plat  :primitive lattice vectors.  Used if 1s digit mode is zero
Co Outputs
Co   qpp   :specially prepared list consisting of qplst in RLV units, and permutation table
Co         :qpp(i,1:3) contains a qp in multiples of qlat.  The qp are ordered for fast search
Co         :qpp(i,4) is an index pointing to which entry this point originated from.
Co         :qpp(i,5:7) contains the original q list
Cr Remarks
Cr   locqpi performs initialization steps for the fast lookup routine locqp.
Cr   It generates an internal array qpp dimensioned (4,nqlst), which locqp requires.
Cr   Caller must allocate qpp before calling locqpi.
Cu Updates
Cu   30 Mar 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nqlst
      real(8) qplst(3,nqlst),plat(3,3),qpp(nqlst,7),tolq
C ... Dynamically allocated arrays
      integer, allocatable ::  iprm(:)
      real(8), allocatable ::  qloc(:,:)
C ... Local parameters
      integer iq,i
      real(8) ::  qx(3)

C     call prmx('plat',plat,3,3,3)
C     call prmx('qplst',qplst,3,3,nqlst)

C ... qplst as multiples of reciprocal lattice vectors
      allocate(qloc(3,nqlst),iprm(nqlst))
      if (mod(mode/10,10) == 0) then
        call dgemm('T','N',3,nqlst,3,1d0,plat,3,qplst,3,0d0,qloc,3)
      else
        call dcopy(3*nqlst,qplst,1,qloc,1)
      endif
C     call prmx('qplst * plat',qloc,3,3,nqlst)

C ... Render all qp to first octant; swap 1st and third colums for sorting
      do  iq = 1, nqlst
        call q2octant(mod(mode/1000,10),tolq,qloc(1,iq),qx)
        qloc(:,iq) = [qx(3), qx(2), qx(1)]
      enddo

C ... Sort the qp, preserve points and permutation array
      call dvheap(3,nqlst,qloc,iprm,tolq,0)
      do  iq = 1, nqlst
        i = iprm(iq)
        qpp(iq,1:7) = [qloc(3,iq), qloc(2,iq), qloc(1,iq), dble(i), qplst(1,i), qplst(2,i), qplst(3,i)]
      enddo
C     call prmx('qplst * plat sorted',qpp,nqlst,nqlst,4)

      deallocate(qloc,iprm)

      end

C     Tests locqp.  Requires files plat and qin
C      subroutine fmain
C      implicit none
C
C      integer ifi,i,j,iq,rdops,nr,nc
C      real(8) xx,plat(3,3),qlat(3,3),qout(3),qx(3)
CC     double precision platt(3,3)
C      real(8), allocatable :: qin(:,:),pqin(:,:),qpp(:,:)
C      real(8),parameter:: tolq=1d-7
C      procedure(logical) :: latvec
C      procedure(integer) :: fopng,rdm
C      procedure(real) :: ran1
C
C      print *, 'reading files plat and qin'
C
C      rdops = 0; nr = 0; nc = 0
C      ifi = fopng('plat',-1,1)
C      i = rdm(ifi,rdops,9,' ',plat,3,3)
C
C      rdops = 0; nr = 0; nc = 0
C      ifi = fopng('qin',-1,1)
C      i = rdm(ifi,rdops,0,' ',[xx],nr,nc)
C      if (nc /= 3) call rx('file qin must have three columns')
C      allocate(qin(nc,nr),pqin(nr,nc))
C      rewind ifi
C      i = rdm(ifi,rdops,nr*nc,' ',pqin,nr,nc)
C      call fclose(ifi)
C      print *, 'qin has',nr,'vectors'
C      forall (i=1:nr,j=1:nc) qin(j,i) = pqin(i,j)
CC     call prmx('qplst',qin,3,3,nr)
CC     deallocate(pqin); allocate(pqin(nc,nr))
C
C      allocate(qpp(nr,7))
C      call locqpi(1000,nr,qin,tolq,plat,qpp)
C
CC     call mkqlat(plat,qlat,xx)
C      call dinv33(plat,1,qlat,xx)
C
CC     platt = transpose(plat)
C      call ran1in(1)
C
C      print *, 'qp?'
C      qx = [0d0,0d0,0.077310020984095812d0]
C      qx = [0d0,0d0,0.038655d0]
C      qx = [0.01d0,0d0,0.d0]
C      read(*,*) qx
C      call info2(10,0,0,'find connection to q = %3;11,7D',qx,2)
C      call locqp(000,nr,qpp,plat,qx,tolq/10,iq,qout)
C      call info2(10,0,0,'iq with tolq = 1d-8 %i',iq,2)
C      call locqp(000,nr,qpp,plat,qx,tolq*10,iq,qout)
C      call info2(10,0,0,'iq with tolq = 1d-7 %i  mismatch = %3;11,7D',iq,qout-qx)
C
C      do  i = 1, nr
C        qx(:) = qin(:,i)
C        xx = nint(10*ran1()-5)
C        qx(:) = qx(:) + xx*qlat(:,1)
C        xx = nint(10*ran1()-5)
C        qx(:) = qx(:) + xx*qlat(:,2) + xx/10*1d-8
C        xx = nint(10*ran1()-5)
C        qx(:) = qx(:) + xx*qlat(:,3)
C
C        print *, i
C        call locqp(100,nr,qpp,plat,qx,tolq,iq,qout)
C        if (iq /= i) stop 'oops'
C        if (.not. latvec(1,tolq,plat,qx(:)-qout(:))) stop 'oops'
C
C      enddo
C      print *, 'OK for given points'
C
C      end
