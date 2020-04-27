      subroutine cmpi123(mode,n1,n2,n3,n123)
C- Packs/unpacks three integers into a single integer
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :-1 setup mode to set up internal ndmx.
Ci         :   Supply max values for n1,n2,n3.
Ci         :   See remarks for limits to their collective size;
Ci         :0  unpack mode: unpacks n1,n2,n3 from given n123
Ci         :1  pack mode: packs given n1,n2,n3 into n123
Ci         :2  Returns nd(1) in n1  and nd(2) in n2
Co Inputs/Outputs
Ci n1,n1,n3: triplet of integers
Ci         : Input for modes -1 and 1; output for mode 0
Ci         : mode 2: nd(1:2) returned in n1,n2
Ci   n123  : packed triplet.
Ci         : Output for modes -1 and 1; input for mode 0
Cr Remarks
Cr   Compacted number = n1 + nd(1)*n2 + nd(2)*n3
Cr   nd(1) and nd(2) are the smallest powers of 2 that encompass n1 or (n1,n2)
Cr   cmpi123 returns an error if the largest compacted number exceeds the largest integer.
Cu Updates
Cu   09 Oct 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,n1,n2,n3,n123
C ... Local parameters
      integer i,ipow,n,nd(2),nmax(3)
      procedure(integer) :: i1mach
      data nd /0,0/, nmax /0,0,0/
      save nd,nmax

      if (mode < 0) then
        nmax = [n1,n2,n3]
        do  i = 1, 2
          ipow = 0
          n = n1; if (i == 2) n = n2
          do  while (n>0)
            ipow = ipow+1
            n = n/2
          enddo
          nd(i) = 2**ipow
        enddo
        i = i1mach(9)
        if (dble(nd(1))*nd(2) > i)
     .    call rx1('cmpi123: triplet too large to compact: %3:1i',nmax)
        nd(2) = nd(1)*nd(2)
        if (n1+dble(nd(1))*n2+dble(nd(2))*n3 > i)
     .    call rx1('cmpi123: triplet too large to compact: %3:1i',nmax)
      elseif (mode == 2) then
        n1 = nd(1)
        n2 = nd(2)
      elseif (mode > 0) then
        if (n1>nmax(1) .or. n2>nmax(2) .or. n3>nmax(3))
     .    call rx2('cmpi123: triplet %3:1i exceeds nmax: %3:1i', [n1,n2,n3],nmax)
        n123 = n1 + nd(1)*n2 + nd(2)*n3
      else
        n3 = n123/nd(2); i = n123-n3*nd(2)
        n2 = i/nd(1); n1 = i-n2*nd(1)
      endif

      end
C      subroutine fmain
C      integer n(3),m(3),d(3),n123
C      logical lok
C
CC      n(1) = 5000
CC      n(2) = 2048
CC      n(3) = 11
C
C      n(1) = 129
C      n(2) = 2048
C      n(3) = 1023
C
C      call cmpi123(-1,n(1),n(2),n(3),n123)
C
C      m = 0
C      call cmpi123(2,m(1),m(2),m(3),n123)
C      call info2(1,0,0,' nd(1) = %i  nd(2) = %i',m(1),m(2))
C
C      call cmpi123(1,n(1),n(2),n(3),n123)
C      call cmpi123(0,m(1),m(2),m(3),n123)
C      d = m - n;
C      lok = d(1)==0 .and. d(2)==0 .and. d(3)==0
C      call info5(1,0,0,' compacting %3:1i -> %i decompressing to %3:1i ok=%l',n,n123,m,lok,5)
C
C      n = [40,20,10]
C      call cmpi123(1,n(1),n(2),n(3),n123)
C      call cmpi123(0,m(1),m(2),m(3),n123)
C      d = m - n;
C      lok = d(1)==0 .and. d(2)==0 .and. d(3)==0
C      call info5(1,0,0,' compacting %3:1i -> %i decompressing to %3:1i ok=%l',n,n123,m,lok,5)
C
C      n = [130,20,10]
C      call cmpi123(1,n(1),n(2),n(3),n123)
C      call cmpi123(0,m(1),m(2),m(3),n123)
C      d = m - n;
C      lok = d(1)==0 .and. d(2)==0 .and. d(3)==0
C      call info5(1,0,0,' compacting %3:1i -> %i decompressing to %3:1i ok=%l',n,n123,m,lok,5)
C
C      end
