      subroutine spcgrp(mode,gen,ngen,plat,nbas,bas,ips,ngmx,g,ag,ng)
C- Sets up space group given point group generators gen
C ----------------------------------------------------------------------
Ci  mode   :a compound set of switches
Ci         :1s digit add inversion to list of generators
Ci         :(from time-reversal symmetry); assume ag is zero.
Ci         :THIS OPTION COMMENTED OUT
Ci         :10s digit 0 stop if any gen is missing a translation.
Ci         :          1 continue anyway, making the point group.
Ci  gen,ngen: generators of the point group, and number
Ci  plat,nbas,bas: primitive lattice vectors and basis vectors
Ci  ips    :species pointer table
Ci  ngmx   :maximum allowed number of group operations
Co  g,ag,ng: space group that defined by operation on vector v
Co         :(g,ag)(v) as g*v + ag, with g a 3x3 matrix and vector ag.
Co         :ng < 0 spcgrp could not find translations for generator -ng
Cu Updates
Cu   04 Apr 03 To determine whether space groups are equal,
Cu             compare only point group part. This eliminates
Cu             space groups that that have the same point group
Cu             but differing translational symmetry which can
Cu             occur if an artifically large supercell is made
C ----------------------------------------------------------------------
      implicit none
      integer mode,ips(*),ngmx,ngen,ng,nbas
      double precision plat(3,3),gen(9,ngen),g(9,ngmx),ag(3,ngmx),
     .  bas(3,nbas)
      integer lok,ig,ib,jb,i1mach,iprint
      double precision xx,qb(3,3),trans(3,50)
C      integer awrite
C      character*80 outs

      if (ngen > 50) call rx('spcgrp: ngen gt 50')

C --- Find non-primitive translations for each generator ---
      call pshpr(iprint()-30)
      call dinv33(plat,1,qb,xx)
      do  ig = 1, ngen
        do  ib = 1, nbas
          do  jb = 1, ib
            call trysop(ib,jb,gen(1,ig),nbas,ips,bas,plat,qb,lok,
     .        trans(1,ig))
            if (lok == 1 .and. iprint() > 20) call awrit2(
     .        ' spcgrp gen %i: found trans =%3:1;6d',' ',80,i1mach(2),
     .        ig,trans(1,ig))
            if (lok == 1) goto 10
          enddo
        enddo
        if (mod(mode/10,10) == 0) call fexit(-1,111,
     .    '%N  SPCGRP: no translation found for generator %i',ig)
        call dpzero(trans(1,ig),3)
   10 enddo
      call poppr

C --- Add inversion ---
C      if (mod(mode,10) == 1) then
C        ngen = ngen+1
C        call dpzero(trans(1,ngen),3)
C        call dpzero(gen(1,ngen),9)
C        gen(1,ngen) = -1
C        gen(5,ngen) = -1
C        gen(9,ngen) = -1
C        ib = awrite('%x%N SPCGRP: added inversion to group operations'
C     .    ,outs,80,0,0,0,0,0,0,0,0,0)
C        if (iprint() >= 30) call cwrite(outs,0,ib-1,0)
C      endif

C --- Make the space group ---
      call sgroup(0,gen,trans,ngen,g,ag,ng,ngmx,qb)
      end
      subroutine sgroup(mode,gen,agen,ngen,g,ag,ng,ngmx,qb)
C- Sets up space group given generators (gen,agen).
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :a compound set of switches
Ci         :1s digit
Ci         :0 two groups compare to equal when both their point
Ci         :  and space parts compare equal
Ci         :1 two groups compare to equal when their point
Ci         :  group compares equal.  This eliminates
Ci         :  space groups that that have the same point group
Ci         :  but differing translational symmetry, which can
Ci         :  occur for artifically large supercells
Ci         :10s digit
Ci         :0 if ng>ngmx, abort with error message
Ci         :1 if ng>ngmx, return with ng=ngmx+1
Ci   gen   :rotation part of generators of the group
Ci   agen  :translation part of space group generator
Ci   ngen  :number of generators
Ci   ngmx  :maximum allowed number of group operations
Ci   qb    :vectors of a microcell in the Brillouin zone
Co Outputs
Co   g     :point group operations
Co   ag    :translation part of space group
Co   ng    :number of group operations
Cr Remarks
Cr   Operations are defined as (g,a)(v):=g*v+a
Cr   where g is a (3x3) matrix, a is a vector.
Cr   Always returns the identity operation as one group operation
Cu Updates
Cu   04 Jan 06 Added 10s digit mode
Cu   14 Mar 03 Added mode
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ngen,ng,ngmx
      double precision gen(9,ngen),g(9,ngmx),qb(3,3),
     .  agen(3,ngen),ag(3,ngmx)
C ... Local parameters
      double precision h(9),hh(9),e(9),sig(9),asig(3),ah(3),ahh(3),ae(3)
      integer ipr,igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,is,nnew,n,m,
     .  stdo,mode0,mode1
      procedure(integer) :: nglob
      logical spgeql
      character sout*80,sg*35
      data e/1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/, ae/0d0,0d0,0d0/

      call getpr(ipr)
      stdo = nglob('stdo')
      sout = ' '
      call spgcop(e,ae,g,ag)
      ng = 1
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

C --- For each generator, do ---
      do  igen = 1, ngen
      call spgcop(gen(1,igen),agen(1,igen),sig,asig)

C --- Extend the group by all products with sig ----
      do  ig = 1, ng
      if (spgeql(mode0,g(1,ig),ag(1,ig),sig,asig,qb)) then
        if (ipr > 30) call awrit2(' Generator %i already in group '//
     .    'as element %i',' ',80,stdo,igen,ig)
C        write(stdo,650) igen,ig
C  650   format(' generator',i3,'  is already in group as element',i3)
            goto 10
      endif
      enddo

C ... Determine order (= power of sig that restores unit operation)
      call spgcop(sig,asig,h,ah)
      do  itry = 1, 100
        iord = itry
        if (spgeql(mode0,h,ah,e,ae,qb)) exit
        call spgprd(sig,asig,h,ah,h,ah)
      enddo

C ... Products of type  g1 sig**p g2
      nnow = ng
      if (ipr >= 40) call awrit2('%a  %i is %i,',sout,80,0,igen,iord)
      do  j = 1, ng
      call spgcop(g(1,j),ag(1,j),h,ah)
      do  ip = 1, iord-1
        call spgprd(sig,asig,h,ah,h,ah)
        do  i = 1, ng
          call spgprd(g(1,i),ag(1,i),h,ah,hh,ahh)
          do  k = 1, nnow
            if (spgeql(mode0,g(1,k),ag(1,k),hh,ahh,qb)) goto 1
          enddo
C         call asymop(hh,ahh,' ',sg)
C         write(stdo,'('' sgroup adding'',i3,2x,a)') nnow+1,sg
          nnow = nnow+1
          if (nnow > ngmx) goto 99
          call spgcop(hh,ahh,g(1,nnow),ag(1,nnow))
    1   enddo
      enddo
      if (j == 1) n2 = nnow
      enddo

C ... Products with more than one sandwiched sigma-factor
      m1 = ng+1
      m2 = nnow
      do  is = 2, 50
      nnew = 0
      do  n = ng+1, n2
      do  m = m1, m2
       call spgprd(g(1,n),ag(1,n),g(1,m),ag(1,m),h,ah)
       do  k = 1, nnow
         if (spgeql(mode0,g(1,k),ag(1,k),h,ah,qb)) goto 2
       enddo
       nnew = nnew+1
       nnow = nnow+1
       if (nnow > ngmx) goto 99
       call spgcop(h,ah,g(1,nnow),ag(1,nnow))
    2 enddo
      enddo
      m1 = m2+1
      m2 = nnow
      if (nnew == 0) exit
      enddo
      ng = nnow
   10 enddo

C --- Printout ---
      if (ipr >= 30) then
        if (sout /= ' ' .and. ipr >= 60) call awrit0
     .    (' Order of generator'//sout//'%a%b',' ',80,stdo)
        call awrit2(' SGROUP: %i symmetry operations from %i '//
     .    'generators',' ',80,stdo,ng,ngen)
        if (ipr >= 60 .and. ng > 1) then
          write(stdo,'('' ig  group op'')')
          do  ig = 1, ng
            call asymop(g(1,ig),ag(1,ig),' ',sg)
            write(stdo,'(i4,2x,a)') ig,sg
          enddo
        endif
      endif
      return

   99 continue
      if (mode1 == 0) call rx1(
     . 'SGROUP: ng greater than ngmx=%i: probably bad translation',ngmx)
C      call info2(2,0,0,
C     .  ' SGROUP (warning) ng greater than ngmx=%i ... exiting',ngmx,0)
      ng = ngmx+1
      end
      logical function spgeql(mode,g1,a1,g2,a2,qb)
C- Determines whether space group op g1 is equal to g2
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 -> space group is compared
Ci         :1 -> only point group is compared
Ci   g1,a1 :first space group
Ci   g2,a2 :second space group
Ci   qb    :reciprocal lattice vectors
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode
      double precision g1(9),g2(9),a1(3),a2(3),qb(3,3)
C ... Local parameters
      integer m,iq,iac
      double precision c,ca,dc

      spgeql=.true.
      do  m = 1, 9
      if (dabs(g1(m)-g2(m)) > 1d-5) then
        spgeql=.false.
        return
      endif
      enddo
      if (mode == 1) return
      do  iq = 1, 3
        c = (a1(1)-a2(1))*qb(1,iq)+(a1(2)-a2(2))*qb(2,iq)
     .     +(a1(3)-a2(3))*qb(3,iq)
        ca = dabs(c)
        iac = ca+0.5d0
        dc = ca-iac
        if (dabs(dc) > 1.d-5) then
          spgeql = .false.
          return
        endif
      enddo
      return
      end
      subroutine spgprd(g1,a1,g2,a2,g,a)
      implicit none
      double precision
     .  g1(3,3),g2(3,3),g(3,3),sum,a1(3),a2(3),a(3),h(3,3),ah(3)
      integer i,j,k

      do  i = 1, 3
        do  j = 1, 3
          sum=0d0
          do  k = 1, 3
            sum = sum + g1(i,k)*g2(k,j)
          enddo
          h(i,j) = sum
        enddo
      enddo
      do  j = 1, 3
        do  i = 1, 3
          g(i,j) = h(i,j)
        enddo
      enddo
      do  i = 1, 3
        ah(i) = a1(i)
        do  j = 1, 3
          ah(i) = ah(i)+g1(i,j)*a2(j)
        enddo
      enddo
      do  i = 1, 3
        a(i) = ah(i)
      enddo
      end

      subroutine spgcop(g,ag,h,ah)
      integer i
      double precision h(9),g(9),ag(3),ah(3)
      do  i = 1, 9
        h(i) = g(i)
        if (dabs(h(i)) < 1.d-10) h(i) = 0d0
      enddo
      do  i = 1, 3
        ah(i) = ag(i)
        if (dabs(ah(i)) < 1.d-10) ah(i) = 0d0
      enddo
      end
      subroutine gpfndx(g,ag,ia,ja,pos,nrc,rb,qb)
C- Finds atom ja which is transformed into ia by group operation g,ag.
      implicit none
      integer ia,ja
      double precision g(3,3),ag(3),pos(3,1),d(3),rb(3,3),qb(3,3)
      integer ka,nrc,m,k
C     integer mode(3)
C      mode(1) = 2
C      mode(2) = 2
C      mode(3) = 2
      ja = 0
      do  ka = 1, nrc
        do  m = 1, 3
          d(m) = ag(m) - pos(m,ia)
          do  k = 1, 3
            d(m) = d(m) + g(m,k)*pos(k,ka)
          enddo
        enddo
        call shorbz(d,d,rb,qb)
C       call shorps(1,rb,mode,d,d)
        if (abs(d(1))+abs(d(2))+abs(d(3)) < 1d-4) then
          ja = ka
          return
        endif
      enddo
      end
