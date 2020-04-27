      subroutine fadflb(mode,lmxa,z,lmxb,iqocc)
C- Autosets LMTO l-cutoff for 1st and 2nd kappa,
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : Choice of lmxb, depending on mode
Ci         : Larger numbers make default size larger.
Ci         : Let lq = l for highest highest state in free atom
Ci         : 0 standard minimal basis, typically spd; same as 3
Ci         : 1 hyperminimal basis: lmxb=lq
Ci         :   lmxb2 = lq
Ci         : 2 lmxb=lq, but increase lxmb by 1 if lq<2
Ci         :   lmxb2 = lq
Ci         : 3 standard minimal basis, typically spd, same as 0
Ci         :   Rules:
Ci             a) start with lmxb=min(lq+1,3)
Ci             b) restrict lmxb to lmxb<=2 if Z<=18 (Ar)
Ci             c) set lmxb=3 if Z>36
Ci         :   lmxb2 = max(lq,lmxb-1)
Ci         : 4 standard basis, similar to mode=3, except:
Ci             b) restrict lmxb to lmxb<=2 if Z<=36 (Kr)
Ci         :   lmxb2 = max(lq,lmxb-1)
Ci         : 5 large basis: increment standard basis by adding 1 l
Ci             Rules:
Ci             a) start with lmxb=max(lq+1,2)
Ci             b) If Z>=6 (C) lmxb=max(lmxb+1,3)
Ci         :   lmxb2 = max(lq,lmxb-1)
Ci   lmxa  :augmentation l-cutoff (lmxb cannot exceed lmxa)
Ci   z     :nuclear charge
Co Outputs
Co   lmxb  :basis l-cutoff, for RSMH,RSMH2
Co   iqocc :
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmxa,lmxb(2),iqocc(5)
      double precision z
C ... Local parameters
      integer n0,lmxq,i,iprint,iz
      parameter (n0=10)
C     parameter (n0=10, NULLI=-99999)
      double precision p(n0,2),q(n0,2)
      character slabl*2,sout*100

C#ifdefC DEBUG
C      integer iiz
C      call setpr(90)
C      call info2(0,0,0,'-------',i,lmxq)
C      do  iiz  = 1, 100
C      z = iiz
C#endif

C ... Find highest l which contains valence charge
      call dpzero(p,2*n0)
      call dpzero(q,2*n0)
      call pshpr(1)
      call defpq(1,z,4,1,p,q)
      call poppr
      lmxq = -1
      do  i = 1, 5
        iqocc(i) = q(i,1)
        if (q(i,1) /= 0) lmxq = i-1
      enddo
      iz = nint(z)

C     Starting defaults: hyperminimal cases (mode 1)
      lmxb(1) = lmxq
      lmxb(2) = lmxq
C     Other cases
      if (mode == 1) then
        lmxb(2) = min(lmxq,1)
      elseif (mode == 2) then
        if (lmxq < 2) lmxb(1) = lmxq+1
      elseif (mode == 0 .or. mode == 3 .or. mode == 4) then
        lmxb(1) = min(lmxq+1,3)
        if (iz <= 36 .and. (mode == 0 .or. mode == 3) .or.
     .      iz <= 18 .and. (mode == 4)) lmxb(1) = min(lmxq+1,2)
C       Large atoms
        if (iz > 36) lmxb(1) = 3
        if (iz > 3) lmxb(1) = max(lmxb(1),2)
        lmxb(2) = max(lmxq,lmxb(1)-1)
      elseif (mode == 5) then
        lmxb(1) = max(lmxq+1,2)
        if (iz >= 6) lmxb(1) = max(lmxb(1),3)
        lmxb(2) = max(lmxq,lmxb(1)-1)
      endif

C ... Printout
      if (iprint() >= 80) then
        call zslabl(-1,slabl,iz)
        if (iprint() <= 90) then
          call info5(0,0,0,' fadflb mode %i, spec '//slabl//
     .      ':  Z=%i lmxb = %i, %i',mode,iz,lmxb(1),lmxb(2),0)
        else
          call info5(0,0,0,' '//slabl//': Z=%i lmxq=%i   %n:1d',iz,lmxq,
     .      lmxq+1,q,0)
        endif
      endif

      if (lmxb(1) > lmxa) then
        call zslabl(-1,slabl,iz)
        call awrit3('%xfadflb spec '//slabl//'%a, Z=%i: lmxb(=%i) '//
     .    'exceeds lmxa(=%i).  Increase lmxa for this species',sout,
     .    len(sout),0,iz,lmxb(1),lmxa)
        call rx(sout)
      endif

C#ifdefC DEBUG
C      enddo
C      stop
C#endif

      end
