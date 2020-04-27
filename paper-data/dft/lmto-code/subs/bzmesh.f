C#define DEBUG
      subroutine bzmshp(prgnam,mode,nkabc,lshft,plat,g,ng,nsgrp,
     .  lfbz,lgstar,lpbc,qb,nqp,qp,wtkp,ipq,igstar)
C- Generate irreducible k-points and attendant arrays for a uniform mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit
Ci         :-1 return nqp only; no other arguments are modified
Ci         :0 return nqp,qb; leave qp,wtkp,ipq,igstar untouched
Ci         :1 return nkp,qb,qp; leave wtkp,ipq,igstar untouched
Ci         :2 return nkp,qb,qp,wtkp; leave ipq,igstar untouched
Ci         :3 return nkp,qb,qp,wtkp,ipq; leave igstar untouched
Ci         :4 return nkp,qb,qp,wtkp,ipq,igstar (igstar according to lgstar)
Ci         :10s digit
Ci         :1 cause bzmesh to identify k-point closest to qp(:,1)
Ci         :  mode must be nonegative
Ci  nkabc  :Mesh of k-points
Ci  lshft  :lshft(i)=1 if qp is offset along axis i
Ci  plat   :primitive lattice vectors, in units of alat
Ci  g      :point group operations
Ci  ng     :number of point group operations
Ci  nsgrp  :number of space group operations.
Ci         :Choose 0 or ng not to distingish from ng.
Ci         :The effect of 0<nsgrp<ng is to affect wtkp
Ci         :See bzmesh.f and note remarks about wgt(1).
Ci  lfbz   :T => make qp for full BZ
Ci  lgstar : nozero, generate igstar according to bzmesh, which see
Ci         : 0 igstar is not made
Ci         : 2 igstar contains inverse mapping of ipq
Ci         :-2 igstar contains group ops rotating irreducible
Ci         :   to full BZ.
Ci  lpbc   :Used to restrict to 2D mesh or k-point order.
Ci         :See bzmesh.f
Co Outputs
Co   nqp   :number of irreducible k-points
Co   qb    :lattice vectors for reduced mesh
Co         :Returned for mode>=0
Co   qp    :irreducible k-points.
Co         :qp needs be dimensioned at least (3,nqp)
Co         :Returned for mode>=1
Co   wtkp  :k-point weights
Co         :wtkp needs be dimensioned at least (nqp)
Co         :Returned for mode>=2
Co   ipq   :ipq(i1,i2,i3) points to irreducible qp corresponding to
Co         :qp labelled (i1,i2,i3), where i1=1..n1, i2=1..n2, i3=1..n3
Co         :ipq needs be dimensioned at nkabc(1)*nkabc(2)*nkabc(3)
Co         :Returned for mode>=3
Co   igstar:information needed for the mapping of qp to another set.
Co         :See bzmesh for details.
Co         :igstar needs be dimensioned at nkabc(1)*nkabc(2)*nkabc(3)+1
Co         :Returned for mode>=4
Cr Remarks
Cu Updates
Cu   25 Aug 18 New 10s digit mode
Cu   18 Jun 13 Repackaged from mkqp
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical, intent(in) :: lfbz
      character prgnam*(*)
      integer mode,nqp,ng,nsgrp,lpbc,lgstar
      integer nkabc(3),lshft(3),ipq(*),igstar(*)
      double precision qp(3,*),wtkp(*),plat(3,3),g(9,ng),qb(3,3)
C ... Local parameters
      character outs*(180)
      integer mxkp,i,ngl,isw
      logical llshft(3)
      double precision qbl(3,3),qpsave(3)
      integer, allocatable::  igstarl(:),ipql(:)
      real(8), allocatable::  qpl(:),wtkl(:)
C     Debugging iqstar
C      integer i1,i2,i3,i123(3),ifac(3),iq,iqs
C      equivalence (i1,i123(1)),(i2,i123(2)),(i3,i123(3))
C      double precision q1(3)

      outs = ' '//prgnam//':  '//
     .  'Generate qp for %?;n;full;irr; BZ with:*'//
     .  '%?;n>0;%b qp weights.;;%-1j'
     .  //'%?;(n>1);%b, ipq.;;%-1j'
     .  //'%?;(n>2);%b, igstar.;;%-1j'
      call info2(45,1,0,trim(outs),isw(lfbz),mod(iabs(max(mode,0)),10))
      if (mode >= 10) qpsave(1:3) = qp(1:3,1)

      if (mod(lpbc,10) == 1) then
C        outs = ' ' // prgnam
C        if (nkabc(3) > 1 .and. iprint() >= 10) then
C          write(stdo,*) ' '
C          call awrit2('%a (warning): nk3=%i, shft3=%i; reset to 1,0',
C     .      outs,80,-stdo,nkabc(3),lshft)
C        endif
        lshft(3) = 0; nkabc(3) = 1
      endif
      mxkp = nkabc(1)*nkabc(2)*nkabc(3)
      if (mxkp <= 0) call rx('bzmesh: illegal or missing k-mesh')

      forall (i=1:3) llshft(i) = lshft(i) /= 0 ! Make logical llshft
      allocate(igstarl(mxkp+1),ipql(mxkp))
      allocate(qpl(3*mxkp),wtkl(mxkp))
      call iinit(igstarl,mxkp+1)
      igstarl(1) = lgstar; wtkl(1) = nsgrp
      ngl = ng
      if (lfbz) ngl = 0
      i = mxkp
      if (mode >= 10) then
        i = -mxkp
        qpl(1:3) = qpsave(1:3)
      endif
      call bzmesh(plat,qbl,nkabc(1),nkabc(2),nkabc(3),llshft,g,ngl,ipql,qpl,wtkl,nqp,i,igstarl,lpbc)
      if (mode >= 0) call dcopy(9,qbl,1,qb,1)
      if (mode >= 1) call dcopy(3*nqp,qpl,1,qp,1)
      if (mode >= 2) call dcopy(nqp,wtkl,1,wtkp,1)
      if (mode >= 3) call icopy(mxkp,ipql,1,ipq,1)
      if (mode >= 4 .and. lgstar /= 0) call icopy(mxkp+1,igstarl,1,igstar,1)
      deallocate(ipql,igstarl,qpl,wtkl)

C     Debugging test of iqstar
C      if (mode < 3) return
C      do  i = 1,3
C        ifac(i) = 1
C        if (llshft(i)) ifac(i) = 2
C      enddo
C      do  iq = 1, nqp
C        iqs = 0
C        do  while (.true.)      ! loop over star of q
CC         call iqstar(iq,nkabc(1),nkabc(2),nkabc(3),ipq,ifac,qb,iqs,q1)
C          call iqstarx(iq,nkabc(1),nkabc(2),nkabc(3),ipq,ifac,qbl,iqs,q1,i123)
C          if (iqs == 0) exit
C          call info8(2,0,0,' (%i,%i,%i)%15p%3:1,6;11D %,5i  %,3i',
C     .      i1,i2,i3,q1,iq,igstar(1+iqs),0,0)
C        enddo
C      enddo
C      stop

      end

      subroutine bzmesh(plat,qb,n1,n2,n3,lshft,g,ng,ipq,qp,wgt,nq,nqmx,igstar,lpbc)
C- Divides the reciprocal lattice into microcells
C-----------------------------------------------------------------------
Ci Inputs:
Ci  plat     :primitive lattice vectors
Ci  n1,n2,n3 :no. divisions for the 3 recip. latt. vecs; (see Remarks)
Ci  lshft    :logical switch, for each recip. latt. vec. :
Ci           :center the mesh through the origin (k=0)
Ci           :center the mesh straddling the origin --- no point at 0.
Ci  g,ng     :symmetry point group operations, and number
Ci  wgt(1)   :if nonzero, wgt(1) holds number of space group ops nsgrp.
Ci           :Setting wgt(1)>0 will cause bzmesh to flag which irreducible
Ci           :k-points have a star originating from symops g(nsgrp+1:ng)
Ci           :Note nsgrp should be less than or equal to ng.
Ci           :Normally these are symops which obtain not from the space
Ci           :group but from time-reversal symmetry.  Thus the
Ci           :calling program can distinguish which irreducible points
Ci           :use this special symmetry.  When switch is set (ie wgt(1)>0)
Ci           :bzmesh returns wgt(i) as a negative number those points
Ci           :which are further reduced by extra symmetries.
Ci  nqmx     :abort if number of k-points exceeds abs(nqmx)
Ci           :nqmx<0 => identify qp that most closely approaches input qp(:,1)
Ci  igstar(0):if nonzero, make igstar.   See remarks for what is made
Ci  lpbc     :1s + 10s digit:
Ci           :  0 make usual 3D mesh
Ci           :  1 or 11 make 2D mesh
Ci           :100s digit
Ci           :  0 standard order of generation of q-points:
Ci                outer loop : i3=1..n3 inner loop :i1=1..n1
Ci           :  1 reverse order of generation of q-points:
Ci                outer loop : i1=1..n1 inner loop :i3=1..n3
Co Outputs:
Co   ipq   :ipq(i1,i2,i3) points to irreducible qp corresponding to
Co         :qp labelled (i1,i2,i3), where i1=1..n1, i2=1..n2, i3=1..n3
Co         :To get the indices i1',i2',i3' corresponding to point
Co         :ipq(i1,i2,i3), call cmpi123.
Co   qp    :qp(1:3,nq) = vector of irreducible qp.
Co         :Note: if nqmx<0 qp(1:3,1) is input
Co   wgt   :wgt(nq) weight assigned to this qp (see Remarks)
Co         :Note the sign of wgt may be used to flag whether the star
Co         :of this irr qp comes at least in part from time-reversal
Co         :symmetry (see wgt(1) above).  Thus the true weight assigned
Co         :to this qp is abs(wgt(i)).  The sign of wgt(i) flags whether
Co         :this point has extra weighting from time-reversal symmetry.
Co   nq    :number of irreducible qp
Co   igstar:information needed for mapping of qp to another set.
Co         :*For a qp specified by the triplet (i1,i2,i3), let
Co         : let i = some irreducible qp with triplet (i1,i2,i3).
Co         : Define 'compressed index' j = i1 + ndmx(1)*i2 + ndmx(2)*i3
Co         : where ndmx is defined in cmpi123
Co         : Thus j contains triplet (i1,i2,i3) in compressed form
Co         :*If input igstar(0)=0, igstar is not touched.
Co         :*If input igstar(0)=2, igstar(i) = compressed index j
Co         : Thus igstar contains data for the 'inverse' of the ipq:
Co         : j = igstar(ipq(i1,i2,i3)) = compressed form of (i1,i2,i3)
Co         : In this mode igstar(1..nq) is generated.
Co         : Also igstar(0) is overwritten with n1+ndmx(1)*n2+ndmx(2)*n3
Co         :*If input igstar(0)=-2, igstar(i) contains the group
Co         : operation ig that maps points in the full BZ to the one of
Co         : the nq irreducible points, that is g(igstar(i)) rotates
Co         : qp at ipq(i1,i2,i3) into qp(i1,i2,i3).
Co         : In this mode igstar(1..n1*n2*n3) is generated
Co         : with igstar(i1+(i2-1)*n1+(i3-1)*n1*n2) = ig.
Co         : Also igstar(0) is overwritten with -(n1+ndmx(1)*n2+ndmx(2)*n3)
Co   qb    :vectors of first microcell for input to BZINTS (see bzmsh0)
Cl Local variables
Cl  lsgrp  :if true, a symmetry op used to map this qp to another qp
Cl         :exceeded input wgt(1)=>start of q contains symmetry from
Cl         :not from space group but from time-reversal symmetry
Cr Remarks:
Cr  The reciprocal lattice is divided into n1*n2*n3 microcells which
Cr  are parallelipipeds with 8 corners.  The corners are nodes of the
Cr  k-space mesh in the whole reciprocal lattice unit cell.
Cr  Thus, for i1=1..n1, i2=1..n2, i3=1..n3 the qp(i1,i2,i3) are
Cr    q_k = (i1*ifac(1)-1)*qb(k,1) +
Cr          (i2*ifac(2)-1)*qb(k,2) +
Cr          (i3*ifac(3)-1)*qb(k,3)
Cr  where ifac is 1 or 2; see bzmsh0.
Cr
Cr  Some of the qp will be symmetry-related, leaving nq irreducible
Cr  k-points, which are returned in qp(1..3,j), j=1,nq.
Cr
Cr  wgt(j) contains the sampling weight associated with qp(j), i.e.
Cr    2/(n1*n2*n3) * no. points equivalent to qp(j); factor 2 for spin.
Cr
Cr  ipq(i1,i2,i3) marks which irreducible qp to which each point in the
Cr  full BZ belongs: point (i1,i2,i3) is equivalent to irreducible
Cr  point ipq(i1,i2,i3).
Cr
Cr  qp gets mapped under symop from qp to reference.
Cr  Thus, if g = rotation, reference is q, q' = mapped qp is obtained from q as
Cr    q' = g^-1 q
Cu Updates
Cu   25 Aug 18 New option to locate closest approach to given qp (nqmx<0)
Cu   09 Oct 17 mxxyz replaced with cmpi123
Cu   09 Jan 09 Package calculation of ndmx into mxxyz, to circumvent
Cu             compiler bugs
Cu   09 Jan 03 Can pass ng=0
Cu   15 Sep 02 Use sign of wgt to flag which irr points contain
Cu             equivalent points from time-reversal symmetry
Cu   21 Jul 02 Further changes to 18 Jun revisions.
Cu   18 Jun 02 New 100s digit in lpbc option to change order in
Cu             generation of qpoints.  Stripped some options in igstar.
Cu   23 Nov 01 Bug fix for long unit cells
Cr   19 Nov 97 (WRL) added lpbc option, projecting qp to 2D
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lshft(3)
      integer n1,n2,n3,nqmx,ng,nq,igstar(0:*),ipq(n1,n2,n3),lpbc
      real(8) qb(3,3),wgt(abs(nqmx)),plat(3,3),qp(3,abs(nqmx)),g(9,1)
C ... Local parameters
      logical lsgrp
      integer i1,i2,i3,ifac(3),ig,igcnt,ii,ii1,ii2,ii3,ipr,iq,is(3),iwgt,
     .  j1,j2,j3,lpbc01,lpbc2,lstar,m1,m2,m3,ndmx(2),nn1,nn2,nn3,stdo,nsgrp
      integer iifind(3)
      real(8) w0,x1,x2,x3,swgt,v(3),v1(3),rb(3,3),qpfind(3),dqmin,
     .  vfind(3),qlat(3,3),qwk(3)
      procedure(integer) :: lgunit
      procedure(real(8)) :: dlength
C ... External calls
      external asymop,awrit5,awrit6,bzmsh0,dinv33,dpzero,fexit,getpr,
     .         grpop,iinit,projql,rx,rxx
      character*1 chr(0:2)
C#ifdef DEBUG
      character*30 sg
      double precision ag(3)
C#endif
C#ifdefC QP0     --- use only points centered within dq from q0
C      integer ifi,fopn
C      double precision q0(3),dq,ddot
CC      data q0 /0.370925d0,-0.401855d0,0.370925d0/,dq/.15d0/
C      data q0 /0.339995d0,0.339995d0,0.339995d0/,dq/.15d0/
C      integer k,m
C#endif

      call getpr(ipr)
      stdo = lgunit(1)
      call cmpi123(-1,n1,n2,n3,ii)
      call cmpi123(2,ndmx(1),ndmx(2),ii,ii)
      qpfind = 0
      if (nqmx < 0) then
        qpfind(1:3) = qp(1:3,1)
        call mkqlat(plat,qlat,w0)
      endif

      lstar = igstar(0)
      lpbc01 = mod(lpbc,100)
      lpbc2  = lpbc/100
      nsgrp = nint(wgt(1))
      chr(2) = ' '
      chr(0) = '*'
      if (nsgrp /= wgt(1)) call rx('bzmesh: invalid input wgt(1)')
C     lpbc2 = 1
      if (min(n1,n2,n3) < 1) call rx('bzmesh: improper specification of k-mesh')
      dqmin = 9d9

      call bzmsh0(plat,lshft,lpbc01,n1,n2,n3,is,ifac,rb,qb)
C#ifdef DEBUG
      if (ipr >= 80) write(stdo,1)
    1 format(/' BZMESH: qp mapping'/'  i1..i3',25x,'qp',20x,'iq   ig  g')
      call dpzero(ag,3)
C#endif
      m1 = n1*ifac(1)
      m2 = n2*ifac(2)
      m3 = n3*ifac(3)
      call iinit(ipq,n1*n2*n3)
      w0 = 2d0/(n1*n2*n3)
      nq = 0
      swgt = 0d0
      igcnt = 0
      nn1 = 6*m1
      nn2 = 6*m2
      nn3 = 6*m3

C --- For each of (n1*n2*n3) qp, find irreducible set ---
C     Loops written without 'do' construct for freedom in choosing
C     loop order
C     do  20  i3 = 0, n3-1
C     do  20  i2 = 0, n2-1
C     do  20  i1 = 0, n1-1
      i1 = 0
      i3 = 0
 10   continue
      i2 = 0
 20   continue
      if (lpbc2 == 0) then
        i1 = 0
      else
        i3 = 0
      endif
   21 continue

C     print 973, i1,i2,i3,igcnt
C 973 format(5i4)

C   ... Add qp to list if not flagged as symmetry-related to a prior
 30   continue
        if (ipq(i1+1,i2+1,i3+1) == 0 .or. nqmx < 0) then
          ii1 = i1*ifac(1)+is(1)
          ii2 = i2*ifac(2)+is(2)
          ii3 = i3*ifac(3)+is(3)
          v1(1) = ii1*qb(1,1) + ii2*qb(1,2) + ii3*qb(1,3)
          v1(2) = ii1*qb(2,1) + ii2*qb(2,2) + ii3*qb(2,3)
          v1(3) = ii1*qb(3,1) + ii2*qb(3,2) + ii3*qb(3,3)
          if (ipq(i1+1,i2+1,i3+1) == 0) v = v1
          if (nqmx < 0) then
            call shorbz(qpfind-v1,qwk,qlat,plat)
            x1 = dlength(3,qwk,1)
            if (x1 < dqmin) then
              iifind(1:3) = [ii1, ii2, ii3]
              vfind(1:3) = v1
              dqmin = x1
            endif
          endif
        endif

        if (ipq(i1+1,i2+1,i3+1) == 0) then
C         call prmx('kpt',v,3,3,1)
C#ifdefC QP0
C          call dcopy(3,v,1,v1,1)
C          call daxpy(3,-1d0,q0,1,v1,1)
CC          print 156, ddot(3,v1,1,v1,1), v1
CC  156     format(4f15.6)
C          if (ddot(3,v1,1,v1,1) > dq*dq) goto 20
C#endif

C     --- Mark each qp in the star of q as equivalent to this qp ---
          iwgt = 0
          lsgrp = .false.
          call dcopy(3,v,1,v1,1)
        do  ig = 1, max(ng,1)
            if (ng > 0) call grpop(v,v1,g,ig)
            x1 = v1(1)*rb(1,1) + v1(2)*rb(2,1) + v1(3)*rb(3,1) - is(1)
            x2 = v1(1)*rb(1,2) + v1(2)*rb(2,2) + v1(3)*rb(3,2) - is(2)
            if (lpbc01 == 0) then
            x3 = v1(1)*rb(1,3) + v1(2)*rb(2,3) + v1(3)*rb(3,3) - is(3)
            else
            x3 = 0
            endif
            j1 = idnint(x1)
            j2 = idnint(x2)
            j3 = idnint(x3)
            if (max(dabs(x1-j1),dabs(x2-j2),dabs(x3-j3)) > 1d-4) then
C              print *, 'ii1,ii2,ii3=',ii1,ii2,ii3
C              print *, 'x1,x2,x3=',sngl(x1),sngl(x2),sngl(x3)
C              print *, 'ig,qb,rb,v,g,v1=',ig
C              call prm('(3f12.6)',qb,3,3)
C              call prm('(3f12.6)',rb,3,3)
C              call prm('(3f12.6)',v,3,1)
C              call prm('(3f12.6)',g(9*ig-8),3,3)
C              call prm('(3f12.6)',v1,3,1)
              call awrit2(' qp%3:1,3;3d -> %3:1,3;3d is not'//
     .          ' on k-mesh',' ',80,stdo,v,v1)
              call rx('BZMESH: symops incompatible with this mesh')
            endif
C           call prmx('rotated kpt',v1,3,3,1)
C            if (ig /= 1 .and. ipr >=
C            call awrit2(' qp%3:1,3;3d rotated from %3:1,3;3d',
C     .        ' ',80,stdo,v1,v)
C ..        scale shifted point or discard if shifted off mesh
            if (lshft(1) .and. mod(abs(j1),2) == 1) cycle
            if (lshft(2) .and. mod(abs(j2),2) == 1) cycle
            if (lshft(3) .and. mod(abs(j3),2) == 1) cycle
            if (lshft(1)) j1 = j1/2
            if (lshft(2)) j2 = j2/2
            if (lshft(3)) j3 = j3/2
C ...       Ensure (j1,j2,j3) in first quadrant of Q
            j1 = mod(j1+2*nn1,n1) + 1
            j2 = mod(j2+2*nn2,n2) + 1
            j3 = mod(j3+2*nn3,n3) + 1
            call rxx(j1 <= 0.or.j2 <= 0.or.j3 <= 0,'neg j in bzmesh')

            if (ipq(j1,j2,j3) == 0) then
              ipq(j1,j2,j3) = nq+1
              iwgt = iwgt+1
              if (ig > nsgrp .and. nsgrp > 0) lsgrp = .true.
              igcnt = igcnt+1
              if (lstar == 2) then
                igstar(nq+1) = j1 + ndmx(1)*j2 + ndmx(2)*j3
              elseif (lstar == -2) then
                ii = j1 + (j2-1)*n1 + (j3-1)*n1*n2
                igstar(ii) = ig
              elseif (lstar /= 0) then
                call rx('bzmesh: bad igstar(0)')
              endif
C              OLD conventions
C              if (lstar == 1) then
C                igstar(igcnt) = j1 + ndmx(1)*j2 + ndmx(2)*j3
C              elseif (lstar == 2) then
C                igstar(nq+1) = j1 + ndmx(1)*j2 + ndmx(2)*j3
C              elseif (lstar == -1) then
C                igstar(igcnt) = ig
C              elseif (lstar == -2) then
C                ii = j1 + (j2-1)*n1 + (j3-1)*n1*n2
C                igstar(ii) = ig
C              endif
C#ifdef DEBUG
              if (ipr >= 80 .and. ng > 0) then
                call dpzero(ag,3)
                call asymop(g(1,ig),ag,' ',sg)
                call awrit6(' (%i,%i,%i)%15p%3:1,6;11D'//
     .          ' %,5i  %,3i '//sg//'%a',' ',80,stdo,j1,j2,j3,
     .          v1,ipq(j1,j2,j3),ig)
              endif
C#endif
            endif
            call rxx(j1 < 0.or.j2 < 0.or.j3 < 0,'neg j in bzmesh')
          enddo
          nq = nq+1
          qp(:,nq) = v(:)
          wgt(nq) = iwgt*w0
          if (lsgrp) wgt(nq) = -iwgt*w0
          swgt = swgt + abs(wgt(nq))
        endif

C     End-of-loop for i1,i2,i3
C  20 continue
      if (lpbc2 == 0) then
        i1 = i1+1
        if (i1 < n1) goto 30
      else
        i3 = i3+1
        if (i3 < n3) goto 30
      endif
      i2 = i2+1
      if (i2 < n2) goto 20
      if (lpbc2 == 0) then
        i3 = i3+1
        if (i3 < n3) goto 10
      else
        i1 = i1+1
        if (i1 < n1) goto 10
      endif
C ... Done accumulating inequivalent qp
      if (ipr >= 20) then
        call awrit6(' BZMESH: %i irreducible QP from %i'//
     .    ' ( %i %i %i )  shift=%3:1l',' ',80,
     .    stdo,nq,n1*n2*n3,n1,n2,n3,lshft)
        call awrit5(' bz nq %i ( %i %i %i )  shift %3:1l',' ',80,
     .    lgunit(2),nq,n1,n2,n3,lshft)

      if (nqmx < 0) then
        call info2(20,0,0,
     .    '%9fpoint closest to q=(%s,%3;6,6d)%N i1..i3%26fqp%23fd',qpfind,2)
        call info5(20,0,0,' (%i,%i,%i)%15p%3:1,6;11D    %;4g',
     .    iifind(1)+1,iifind(2)+1,iifind(3)+1,vfind,dqmin)
      endif

      endif
C#ifdefC  QP0
C      call awrit3('          %i points within dq=%1;6d  swgt=%1;6d',
C     .  ' ',80,stdo,igcnt,dq,swgt)
C      ifi = fopn('QPTS')
C      write(ifi,333) nq, n1,n2,n3,0
C  333 format(4i5,i8)
C      do  70  iq = 1, nq
C   70 write(ifi,334) iq,(qp(m,iq),m=1,3), wgt(iq)
C  334 format(i5,1p,4d18.10)
C#else
      if (lstar /= 0) igstar(0) = n1 + ndmx(1)*n2 + ndmx(2)*n3
      if (lstar < 0) igstar(0) = -igstar(0)
      if (igcnt /= n1*n2*n3) call rx('bug in bzmesh')
      if (dabs(swgt-2) > 1.d-9) call fexit(-1,111,
     .' Exit -1 BZMESH: QP weights sum to %1;6d but must sum to 2',swgt)
C#endif


      if (ipr >= 50) then
        write (stdo,2)
    2   format(14x,'Qx',10x,'Qy',10x,'Qz',6x,'Multiplicity    Weight')
        do  iq = 1, nq
        ii = 1+dsign(1d0,wgt(iq))
        iwgt = abs(wgt(iq)/w0) + .1d0
          write (stdo,3) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,chr(ii),
     .                   abs(wgt(iq))
    3     format(i5,2x,3F12.6,i10,1x,a,f14.6)
        enddo
      endif
      end

      subroutine bzmsh00(plat,lshft,lpbc,n123,ifac,qb)
C- Setup for a uniform mesh in the Brillouin zone
C-----------------------------------------------------------------------
Ci   plat:      primitive lattice vectors
Ci   lshft:     0 center mesh points on the origin for i'th R.L.V
Ci              1 set mesh points off center from the origin
Ci   lpbc:      0 make 3D mesh
Ci              1 or 11 make 2D mesh
Ci   n123:      number of divisions along the three R.L.V
Co Outputs:
Co   ifac(i)    1 if not shifted, 2 if shifted for i'th axis; see qb.
Co   qb:        a microcell in the Brillouin zone
Co              This, together with ifac provides information how to
Co              generate the actual q-point from a triplet of integers
Co              specifying a point on the mesh.  Given a
Co              triplet (j_1,j_2,j_3) of ipq, the i'th component q_i is
Co                 q_i(j_1,j_2,j_3) = sum_n (j'_n)*qb(i,n),  with
Co                 j'_n = j_n*ifac(n)-1
Cu Updates
Cu   10 Aug 15  simplified bzmsh0, returning sufficient info to make qp mesh
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lshft(3),n123(3),ifac(3),lpbc
      double precision plat(3,3),qb(3,3)
C ... Local parameters
      logical llshft(3)
      double precision rb(3,3)
      integer i,is(3)

      do  i = 1, 3
        llshft(i) = lshft(i) /= 0
      enddo
      call bzmsh0(plat,llshft,lpbc,n123(1),n123(2),n123(3),is,ifac,rb,qb)

      end

      subroutine bzmsh0(plat,lshft,lpbc,n1,n2,n3,is,ifac,rb,qb)
C- Setup for a uniform mesh in the Brillouin zone
C-----------------------------------------------------------------------
Ci   plat:      primitive lattice vectors
Ci   n1,n2,n3:  number of divisions along the three R.L.V
Ci   lshft:     F center mesh points on the origin for i'th R.L.V
Ci              T set mesh points off center from the origin
Ci   lpbc:      0 make 3D mesh
Ci              1 or 11 make 2D mesh
Co Outputs:
Co   is(i)      0 mesh points centered at the origin for i'th axis
Co              1 mesh points off-centered from the origin for i'th axis
Co   ifac(i)    1 if not shifted, 2 if shifted for i'th axis; see qb.
Co   rb:        a microcell in the Wigner Seitz cell.
Co   qb:        a microcell in the Brillouin zone
Co              This, together with ifac provides information how to
Co              generate the actual q-point from a triplet of integers
Co              specifying a point on the mesh.  Given a
Co              triplet (j_1,j_2,j_3) of ipq, the i'th component q_i is
Co                 q_i(j_1,j_2,j_3) = sum_n (j'_n)*qb(i,n),  with
Co                 j'_n = j_n*ifac(n)-1
Cu Updates
Cu   19 Nov 97 added lpbc option
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lshft(3)
      integer n1,n2,n3,is(3),ifac(3),lpbc
      double precision plat(3,3),rb(3,3),qb(3,3)
C ... Local parameters
      double precision qlat(3,3),vol
      integer k,m,m1,m2,m3,iprint,lgunit,stdo

      stdo = lgunit(1)
      is(1)   = 0
      is(2)   = 0
      is(3)   = 0
      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      if (lshft(1)) then
        is(1) = 1
        ifac(1) = 2
      endif
      if (lshft(2)) then
        is(2) = 1
        ifac(2) = 2
      endif
      if (lshft(3)) then
        is(3) = 1
        ifac(3) = 2
      endif
      m1 = n1*ifac(1)
      m2 = n2*ifac(2)
      m3 = n3*ifac(3)
      call dinv33(plat,1,qlat,vol)
      if (lpbc == 1 .or. lpbc == 11) call projql(qlat)
      if (iprint() > 60) then
        write (stdo,1)
    1   format(' BZMESH : ',5x,'Plat',31x,'Qlat')
        do  k = 1, 3
          write (stdo,2) (plat(m,k),m=1,3),(qlat(m,k),m=1,3)
    2     format(3F10.5,5x,3F10.5)
        enddo
      endif

      do  m = 1, 3
        if (m1*m2*m3 /= 0) then
          qb(m,1) = qlat(m,1)/m1
          qb(m,2) = qlat(m,2)/m2
          qb(m,3) = qlat(m,3)/m3
        end if
        rb(m,1) = plat(m,1)*m1
        rb(m,2) = plat(m,2)*m2
        rb(m,3) = plat(m,3)*m3
      enddo

      end
C      integer function mxxyz()
CC- Return maximum integer whose cube fits into integer word
CC  Package as a subroutine to avoid compiler bugs, e.g. result changing
CC  as compiler switches change.
CC  e.g. intel ifort v 11, for x86_64, these should produce same results:
CC  ifort -g -cm -axW -WB -c bzmesh.f
CC  ifort -g -cm -c bzmesh.f
C      implicit none
C      integer i1mach
C
C      call rx('bzmesh replace call to mxxyz')
C
C      mxxyz = (dble((i1mach(9)))/2.01d0) ** (1d0/3d0)
C      end

      subroutine iqstar(mode,iq,nk1,nk2,nk3,ipq,ifac,qb,iqs,q,i123)
C- Return next q-point in the star of q, associated with irreducible point iq
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :-1 return number of points in the star of ipq
Ci         : 0 return iqs only
Ci         : 1 return iqs and q
Ci         : 2 return iqs and q and i123
Ci   iq    : index to irreducible q-point
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   ipq   : ipq(i1,i2,i3) points to the irreducible qp into which
Ci         : mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   ifac  : used with qb to recover a wave vector from three indices
Ci         : specifying position in a uniform mesh (bzmsh0.f)
Ci         : used only if mode>0.
Ci   qb    : vectors of a microcell in the Brillouin zone
Ci         : used only if mode>0.
Cio Inputs/Outputs
Cio  iqs   : On input, index to prior star-of-q found in full BZ.
Cio        : iqs should be 0 on whenever iqstar is called for the
Cio        : first time for a new iqs.
Cio        : On output, index to next star-of-q found in full BZ.
Cio        : If there are no more elements in the star, iqs is returned as 0
Cio        : Special case mode -1 :
Cio        : On output, iqs = number of points in star of iq
Co Outputs
Co   q     : (mode>0) qp for the current star-of-q
Co   i123  : (mode>1) indices i1,i2,i3, to microcell in full BZ (bzmesh)
Cl Local variables
Cl         :
Cr Remarks
Cr   To loop over all points in the star of an irreducible point iq, do:
Cr     do  while (.true.)        ! loop over star of q
Cr       call iqstar(iq,nk1,nk2,nk3,ipq,ifac,qb,iq1,q1)
Cr       if (iq1 == 0) exit
Cr        ..
Cr     enddo
Cr
Cr   Given a loop in the full BZ, resolve k as follows
Cr   This loop requires a routine that returns irreducible iq and ig for each qp in the fbz
Cr      do  iqfbz = 1, nkfbz
Cr        call routine-to-get-iq(iq,ig)
Cr        if (ig == 1) iqs = 0  ! First point of new star
Cr        call iqstar(2,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,iqs,qpr,i123) ! recoversb
Cr        call iqstar(2,iq,nk1,nk2,nk3,ipq,ifac,qb,iqs,q,i123)
Cr        print '(3i4,2x,3i3,3f12.6)', iq,ig,iqs,iv(1:3),qpr
Cr      enddo
Cr
Cr   Vector q1 is returned UNSHORTENED.
Cr
Cr   iqs is an index to a point in the full Brillouin zone, that is,
Cr   the position of this q had no symmetry operations been used.
Cr   See discussion for igstar in subroutine bzmesh.
Cr
Cr   The corresponding point group operation is given by ig=igstar(iqs).
Cr   g(ig) rotates point q(iqs) in start of q to q.
Cr   igstar is not needed by this routine, so ig is not returned.
Cr
Cr   Should they be needed in future, the points i1,i2,i3 in the full BZ
Cr   could be returned
Cu Updates
Cu   05 Aug 15 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,iq,nk1,nk2,nk3,iqs,ifac(3),ipq(nk1*nk2*nk3),i123(3)
      double precision q(3),qb(3,3)
C ... Local parameters
      integer i1,i2,i3,k,j1,j2,j3,nq,nstar
C     integer k1,k2,k3
      double precision qk
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,j1,j2,j3) = (j1*ifac(1)-1)*qb(k,1) +
     .                 (j2*ifac(2)-1)*qb(k,2) +
     .                 (j3*ifac(3)-1)*qb(k,3)



      if (mode == -1) then
        iqs = 0
        nstar = 0
      endif

      nq = nk1*nk2*nk3
      do  while (iqs < nq)

        iqs = iqs+1

C   ... skip this qp unless it is related to iq
        if (ipq(iqs) /= iq) cycle

        if (mode == 0) return
        if (mode == -1) then
          nstar = nstar+1
          cycle
        endif

C   ... indices to microcell in full BZ
        i1 = mod(iqs-1,nk1)+1
        i2 = mod((iqs-1)/nk1,nk2)+1
        i3 = mod((iqs-1)/(nk1*nk2),nk3)+1

C   ... q into which h or z is rotated
        q(1) = qk(1,i1,i2,i3)
        q(2) = qk(2,i1,i2,i3)
        q(3) = qk(3,i1,i2,i3)
        if (mode == 1) return

C   ... copy indices to i123
        i123(1) = i1; i123(2) = i2; i123(3) = i3
        return

      enddo

C ... No additional q-points in star ... reset iqs
      iqs = 0
      if (mode == -1) iqs = nstar

      end

      subroutine alignqp(nqlst,qplst,nqibz,qibz,plat,tolq,iphost,nomatch)
C- Determine if a list of qp is contained in a host list
C ----------------------------------------------------------------------
Ci Inputs
Ci   nqlst :size of list
Ci   qplst :vector of qp in list
Ci   nqibz :size of host list
Ci   qibz  :vector of qp in host list
Ci   plat  :primitive lattice vectors, in units of alat
Ci   tolq  :qp must be equal to within tolerance tolq
Co Outputs
Co   iphost:iphost(iq) points to entry in qibz matching qplst(iq), if one exists
Co         :iphost(iq)=0 if no entry found
Co  nomatch:returns index to first iq for which no match was found
Co         :If all points matched, nomatch returns 0
Cr Remarks
Cu Updates
Cu   12 Mar 18 Adapated from sugws.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nqlst,nqibz,iphost(nqlst),nomatch
      real(8) qplst(3,nqlst),qibz(3,nqibz),plat(3,3),tolq
C ... Local parameters
      integer iq,i
      procedure(logical) :: latvec

      call iinit(iphost,nqlst)
      do  iq = 1, nqlst
        nomatch = iq
        do  i = 1, nqibz
          if (latvec(1,tolq,plat,qplst(1:3,iq)-qibz(1:3,i))) then
            iphost(iq) = i
            exit
          endif
        enddo
        if (iphost(iq) == 0) exit
        if (iq == nqlst) nomatch = 0 ! points are all contained in host
      enddo
      end

C      subroutine iqstarx(iq,nk1,nk2,nk3,ipq,ifac,qb,iqs,q,i123)
CC- Return next q-point in the star of q, associated with irreducible point iq
CC ----------------------------------------------------------------------
CCi Inputs
CCi   iq    : index to irreducible q-point
CCi   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
CCi   ipq   : ipq(i1,i2,i3) points to the irreducible qp into which
CCi           mesh point (i1,i2,i3) is mapped (bzmesh.f)
CCi   ifac  : used with qb to recover a wave vector from three indices
CCi           specifying position in a uniform mesh (bzmsh0.f)
CCi   qb    : vectors of a microcell in the Brillouin zone
CCio Inputs/Outputs
CCio  iqs   : On input, index to last star-of-q found in full BZ.
CCio        : iqs should be 0 on whenever iqstar is called for the
CCio        : first time for a new iqs.
CCio        : On output, index to next star-of-q found in full BZ.
CCio        : If there are no more elements in the star, iqs is returned as 0
CCo Outputs
CCo   q     : qp for the current star-of-q
CCl Local variables
CCl         :
CCr Remarks
CCr   This is a debugging version of iqstar
CCu Updates
CCu   05 Aug 15 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer iq,nk1,nk2,nk3,iqs,ifac(3),ipq(nk1*nk2*nk3),i123(3)
C      double precision q(3),qb(3,3)
CC ... Local parameters
C      integer i1,i2,i3,k,j1,j2,j3,nq
C      integer k1,k2,k3
C      double precision qk
CC Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      qk(k,j1,j2,j3) = (j1*ifac(1)-1)*qb(k,1) +
C     .                 (j2*ifac(2)-1)*qb(k,2) +
C     .                 (j3*ifac(3)-1)*qb(k,3)
C
CC     debugging : slow way that checks algorithm below
C      nq = iqs
C      iqs = 0
C      do  k3 = 1, nk3
C        do  k2 = 1, nk2
C          do  k1 = 1, nk1
C
C            iqs = iqs+1
CC       ... Look for the first iqs > nq
C            if (iqs <= nq) cycle
C
CC       ... skip this qp unless it is related to iq
C            if (ipq(iqs) /= iq) cycle
C
C            i1 = mod(iqs-1,nk1)+1
C            i2 = mod((iqs-1)/nk1,nk2)+1
C            i3 = mod((iqs-1)/(nk1*nk2),nk3)+1
C
C            if (i1 /= k1 .or. i2 /= k2  .or. i3 /= k3) stop 'oops'
C            i123(1) = i1; i123(2) = i2; i123(3) = i3
C
CC       ... q into which h or z is rotated
C            q(1) = qk(1,i1,i2,i3)
C            q(2) = qk(2,i1,i2,i3)
C            q(3) = qk(3,i1,i2,i3)
C
C            return
C
C          enddo
C        enddo
C      enddo
C
CC ... No additional q-points in star ... reset iqs
C      iqs = 0
C
C      end
C      subroutine nirrbz(plat,n1,n2,n3,lshft,g,ng,lpbc,nq)
CC- Find the number of irreducible k-points for a given uniform mesh
CC ----------------------------------------------------------------------
CCi Inputs
CCi  plat     :primitive lattice vectors
CCi  n1,n2,n3 :no. divisions for the 3 recip. latt. vecs; (see Remarks)
CCi  lshft    :logical switch, for each recip. latt. vec. :
CCi           :center the mesh through the origin (k=0)
CCi           :center the mesh straddling the origin --- no point at 0.
CCi  g,ng     :symmetry point group operations, and number
CCi  lpbc     :1s + 10s digit:
CCi           :  0 make usual 3D mesh
CCi           :  1 or 11 make 2D mesh
CCi           :100s digit
CCi           :  0 standard order of generation of q-points:
CCi                outer loop : i3=1..n3 inner loop :i1=1..n1
CCi           :  1 reverse order of generation of q-points:
CCi                outer loop : i1=1..n1 inner loop :i3=1..n3
CCo Outputs
CCo   nq      :number of irreducible qp
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   17 Jul 12 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      logical lshft(3)
C      integer n1,n2,n3,ng,nq,lpbc
C      double precision plat(3,3),g(9,ng)
CC ... Local parameters
C      integer:: nqmx,igstar(1)=0
C      double precision qb(3,3)
C      integer,allocatable:: ipq(:,:,:)
C      real(8),allocatable:: qibz(:,:),wibz(:)
C
C      nqmx = n1*n2*n3
C      allocate(qibz(3,nqmx),wibz(nqmx),ipq(n1,n2,n3))
C      wibz(1) = 0
C      call pshpr(1)
C      call bzmesh(plat,qb,n1,n2,n3,lshft,g,ng,ipq,qibz,wibz,nq,nqmx,
C     .  igstar,lpbc)
C      call poppr
C      deallocate(qibz,wibz,ipq)
C
C      end
C      subroutine fmain
C      implicit none
C
C      integer ifi,i,rdops,nr,nc
C      real(8) xx
C      procedure(integer) :: fopng,rdm
C
C      rdops = 0
C      ifi = fopng('qin',-1,1)
C      i = rdm(ifi,rdops,0,' ',[xx],nr,nc)
C
C      nc = 4
C      if (nc .ne. 3) call rx('file qin must have three columns')
C
C      end
