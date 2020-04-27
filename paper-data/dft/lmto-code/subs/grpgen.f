      subroutine grpgen(gen,ngen,symops,ng,ngmx)
C- Create or extend the point group from the supplied generators
C ----------------------------------------------------------------
Ci Inputs
Ci   gen,ngen,ngmx
Ci   ng  (12 Sep 96): if>0 , add symops to the ng already in list.
Co Outputs
Co   symops,ng
Cr Remarks
Cr   This works for point groups only and is set up for integer
Cr   generators.
C ----------------------------------------------------------------
      implicit none
      integer ngen,ng,ngmx
      double precision gen(9,ngen),symops(9,ngmx)
      double precision h(9),hh(9),e(9),sig(9),ae(3)
      integer igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,n,m
      integer ipr,stdo,nglob
      logical grpeql
      character*80 sout
C     double precision plat(3,3)
C     character*80 sout2
C     logical ltmp
      data e /1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/, ae/0d0,0d0,0d0/

      stdo = nglob('stdo')
      call getpr(ipr)
      sout = ' '

C      do j = 1, ng
C      print 333, (symops(k,j), k=1,9), 0,j
C      enddo

      call grpcop(e,symops)  ! Always the identity: park into 1st op
      if (ng < 1) ng = 1
C --- For each generator, do ---
      do  igen = 1, ngen
        call grpcop(gen(1,igen),sig)
C   ... Skip this generator if it is already in the group
        do  ig = 1, ng
          if (grpeql(symops(1,ig),sig) .and. ipr > 30)
     .      call awrit2(' Generator %i already in group as element %i',
     .      ' ',80,stdo,igen,ig)
          if (grpeql(symops(1,ig),sig)) goto 10
        enddo

C   --- Determine order ---
        call grpcop(sig,h)
        do  itry = 1, 100
          iord = itry
          if (grpeql(h,e)) exit
          call grpprd(sig,h,h)
        enddo

C   --- Products of type  g1 sig**p g2 ---
        nnow = ng
        if (ipr >= 40) call awrit2('%a  %i is %i,',sout,80,0,igen,iord)
        do  j = 1, ng
          call grpcop(symops(1,j),h)
          do  ip = 1, iord-1
C       ... h = sig**ip
            call grpprd(sig,h,h)
            do  i = 1, ng
C         ... hh = symops_i sig**ip
              call grpprd(symops(1,i),h,hh)
              do  k = 1, nnow
                if (grpeql(symops(1,k),hh)) goto 1
              enddo
              nnow = nnow+1
C             call asymopn(1,hh,ae,plat,sout,sout2)
C             print *, j,nnow,trim(sout)
C             print 333, hh, j,nnow
C 333         format(9f12.6,2i3)
              if (nnow > ngmx) goto 20
              call grpcop(hh,symops(1,nnow))
    1       continue
            enddo
          enddo
          if (j == 1) n2 = nnow
        enddo

C   --- Products with more than one sandwiched sigma-factor ---
        m1 = ng+1
        m2 = nnow
        do  i = 2, 50
          do  n = ng+1, n2
            do  m = m1, m2
              call grpprd(symops(1,n),symops(1,m),h)
              do  k = 1, nnow
                if (grpeql(symops(1,k),h)) goto 2
              enddo
              nnow = nnow+1
              if (nnow > ngmx) goto 20
              call grpcop(h,symops(1,nnow))
    2       enddo
          enddo
          if (m2 == nnow) exit
          m1 = m2 + 1
          m2 = nnow
        enddo
        ng = nnow
   10   continue
      enddo

C --- Printout ---
      if (ipr >= 30) then
        if (sout /= ' ' .and. ipr >= 60) call awrit0
     .    (' Order of generator'//sout//'%a%b',' ',80,stdo)
        call awrit2(' GRPGEN: %i symmetry operations from %i '//
     .    'generator(s)',' ',80,stdo,ng,ngen)
      endif
      if (ipr >= 80 .and. ng > 1) then
        write(stdo,'('' ig  group op'')')
        do  ig = 1, ng
          call asymop(symops(1,ig),ae,' ',sout)
          write(stdo,'(i4,2x,a)') ig,sout(1:35)
        enddo
      endif

C      if (ipr >= 110) then
C        print *, 'group operations:'
C        call ywrm(0,' ',1,i1mach(2),'(5f12.6)',symops,1,9,9,ng)
C      endif
      return
   20 continue
      call rx('GRPGEN: too many elements')
      end
      subroutine grpcop(g,h)
C- Copy matrix
      implicit none
      double precision h(9),g(9)
      integer i
      do  i = 1, 9
        h(i) = g(i)
      enddo
      end
      logical function grpeql(g1,g2)
C- Checks if G1 is equal to G2
C  8 Oct 2015 increase tolerance to manage site data with lower resolution
      implicit none
      double precision g1(9),g2(9),dabs,tol,x1,x2
      parameter (tol = 1d-6)
      logical ddif
      integer i
      ddif(x1,x2) = dabs(x1-x2) > tol
      grpeql = .false.
      do  i = 1, 9
        if (ddif(g1(i),g2(i))) return
      enddo
      grpeql = .true.
      end
      subroutine grpprd(g1,g2,g1xg2)
C- Returns the product of two point group operations
      implicit none
      double precision g1(3,3),g2(3,3),g1xg2(3,3),h(3,3),sum
      integer i,j,k

      do  i = 1, 3
        do  j = 1, 3
          sum = 0d0
          do  k = 1, 3
            sum = sum + g1(i,k)*g2(k,j)
          enddo
          h(i,j) = sum
        enddo
      enddo
      do  j = 1, 3
        do  i = 1, 3
          g1xg2(i,j) = h(i,j)
        enddo
      enddo

      end
