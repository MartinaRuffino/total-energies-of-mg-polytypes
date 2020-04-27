      subroutine clcoin(mode,off1,off2,rtab1,rtab2,fuzz,iax1,iax2,
     .  leng1,leng2,match,nmatch)
C- Return a list of coincident lattice points in two clusters
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit governs permutations of lists
Ci          1's bit:  rtab1(:,i1) corresponds to iax1(:,ip1)
Ci            If 0, ip1 = i1
Ci            If 1, ip1 = iax1(7,i1)), i.e. rtab1 is permuted by iax1(7)
Ci          2's bit: rtab2(:,i2) corresponds to iax1(:,ip2)
Ci            If 0, ip2 = i2
Ci            If 1, ip2 = iax2(7,i2)), i.e. rtab2 is permuted by iax2(7)
Ci          3 options 1+2
Ci         10s+100s digit governs matching criterion to be used (may be combined)
Ci          1 check for coincident lattice
Ci          2 require iax1(6,ip1) = iax2(6,ip2)
Ci          4 require iax1(8,ip1) = iax2(8,ip2)
Ci         100s digit (see match, outputs, below)
Ci          0 Look for cluster 1 sites in cluster 2
Ci          1 Look for cluster 2 sites in cluster 1
Ci   off1  :offset in iax1 table to start of 1st cluster
Ci   off2  :offset in iax2 table to start of 2nd cluster
Ci   rtab1 :table of positions in first cluster; see Remarks
Ci   rtab2 :table of positions in second cluster; see Remarks
Ci   fuzz  :tolerance to which positions might deviate but
Ci          still count as coincident
Ci   leng1 :number of pairs in first cluster
Ci   leng2 :number of pairs in second cluster
Ci   iax1  :iax table for first cluster; see pairc.f for description
Ci   iax2  :iax table for second cluster; see pairc.f for description
Co Outputs
Co   match :If site ip1 in iax1 is coincident with site ip2 in iax2,
Co         :match(ip1) = ip2
Co         :If no coincident site is found, match(ip1) = 0
Co  nmatch :number of pairs matched
Cr Remarks
Ci   Effective site index, iax(7,ip), must be set (see pairc.f) and
Ci   rtab1,rtab2 must be sorted by it.
Cb Bugs
C    make fuzz same convention as in ppair4
c    need fancier comparison as in siteid
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C     Passed parameters
      integer mode,off1,off2,leng1,leng2,match(leng1),niax,nmatch
      parameter (niax=10)
      integer iax1(niax,leng1),iax2(niax,leng2)
      double precision rtab1(3,leng1),rtab2(3,leng2)
      double precision fuzz
C     Local parameters
      logical lprm1,lprm2,lmata,lmatb
      integer mode0,mode1,mode2,ip,ip1,ip2,i1,i2

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      lprm1 = mod(mode0,2) /= 0
      lprm2 = mode0 >= 2
      lmata = mod(mode1,2) /= 0
      lmatb = mode1 >= 2

C     i1,i2 = index to next site in clusters 1,2
C     ip1,ip2 = index to next site in permuted order
      i1 = 1 + off1
      i2 = 1 + off2

C --- Match coincidences in permuted order --
      nmatch = 0
      if (mode2 == 0) then
        call iinit(match(i1),leng1)
      else
        call iinit(match(i1),leng2)
      endif
      do  ip = 1, leng1+leng2

C       We are done when every member of either cluster has been looked at
        if (i1 > leng1+off1 .or. i2 > leng2+off2) return

C       ip1,ip2 = index to next site in permuted order
        ip1 = i1
        ip2 = i2
        if (lprm1) ip1 = iax1(7,i1) + off1
        if (lprm2) ip2 = iax2(7,i2) + off2
        if (ip1 > i1+off1) call rx('clcoin: faulty iax(7)')
        if (ip2 > i2+off1) call rx('clcoin: faulty iax(7)')

C       Increment 1 if x(2)>x(1)+fuzz
        if (rtab1(1,ip1)+fuzz < rtab2(1,ip2)) goto 24
C       Increment 2 if x(1)>x(2)+fuzz
        if (abs(rtab1(1,ip1)-rtab2(1,ip2)) > fuzz) goto 25
C       Increment 1 if y(2)>y(1)+fuzz
        if (rtab1(2,ip1)+fuzz < rtab2(2,ip2)) goto 24
C       Increment 2 if y(1)>y(2)+fuzz
        if (abs(rtab1(2,ip1)-rtab2(2,ip2)) > fuzz) goto 25
C       Increment 1 if z(2)>z(1)+fuzz
        if (rtab1(3,ip1)+fuzz < rtab2(3,ip2)) goto 24
C       Increment 2 if z(1)>z(2)+fuzz
        if (abs(rtab1(3,ip1)-rtab2(3,ip2)) > fuzz) goto 25

C   ... Sites are coincident to within tolerance fuzz.
C       Check for points coincident but miss other matching criteria
        if (lmata) then
          if (iax1(6,ip1) /= iax2(6,ip2)) then
            match(ip1) = 0
            goto 22
          endif
        endif
        if (lmatb) then
          if (iax1(8,ip1) /= iax2(8,ip2)) then
            match(ip1) = 0
            goto 22
          endif
        endif
        if (mode2 == 0) then
          match(ip1) = ip2
        else
          match(ip2) = ip1
        endif
c       if (match(ip1) == 0) match(ip1) = ip2
        nmatch = nmatch+1

   22   continue
        i1 = i1 + 1
        i2 = i2 + 1
        ip1 = i1
        ip2 = i2
        if (lprm1) ip1 = iax1(7,i1) + off1
        if (lprm2) ip2 = iax2(7,i2) + off2
        cycle

C  ...  Case rtab(1) is the smaller one: increment its index
   24   continue
        i1 = i1 + 1
        ip1 = i1
        if (lprm1) ip1 = iax1(7,i1) + off1
        cycle

C  ...  Case rtab(2) is the smaller one: increment its index
   25   continue
        i2 = i2 + 1
        ip2 = i2
        if (lprm2) ip2 = iax2(7,i2) + off2
C       cycle

      enddo
      end
