      subroutine pgglst(mode,npl,pgplp,lstii)
C- Make a list of layers for which to calculate GF
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode 1s digit
Ci        0  exclude -1 from the list
Ci        1  include -1 in the list
Ci       10s digit
Ci        0  exclude npl from the list
Ci        1  include npl in the list
Ci      100s digit concerns how repeating layers are incorporated
Ci        0: lstii contains only rightmost GF for for repeating layers
Ci        1: lstii contains left- and rightmost- GF
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Co Outputs
Co   lstii : lstii(-1) number of PL
Co         : lstii(0)  = -1
Co         : lstii(1..)  = PL to be included in list
Cr Remarks
Cb Bugs
Cb   This routine needs cleaning up!
Cu Updates
Cu   15 Sep 04 Printout more compact
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,npl,lstii(-1:npl+2),pgplp(6,-1:npl)
C ... Local parameters
      integer kpl,ipl,ivl,nrpt,i,iprint,mode0,mode1,mode2,iil
      integer nrpt0
      character strn*128
      parameter (nrpt0=3)

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)

C --- Pick off which PL to calculate ---
      kpl = 0
      ipl = 0
      iil = 1
      lstii(-1) = 0
      if (mode0 == 1) then
        lstii(-1) = 1
        lstii(0) = -1
        iil = 0
      endif

   10 continue
      if (ipl >= npl) goto 11

C ... Number of repetitions of potential for next PL's
      nrpt = 1
      ivl = pgplp(2,ipl)
      do  12  i = ipl+1, npl-1
        if (pgplp(2,i) == 0 .or. pgplp(2,i) /= ivl) goto 13
        nrpt = nrpt+1
   12 continue
   13 continue
C ... If the repeating layers propagate all the way to rhs
      if (ipl+nrpt == npl) nrpt = nrpt-1
C ... Minimum number of repeating layers to treat specially
      if (nrpt <= nrpt0) nrpt = 1

      kpl = kpl+1
      if (mode2 == 1 .and. nrpt > 1) then
        lstii(kpl-iil) = ipl
        kpl = kpl+1
      endif
      lstii(kpl-iil) = ipl+nrpt-1
      ipl = ipl+nrpt
      lstii(-1) = kpl
      goto 10
   11 continue

      if (mode1 == 1) then
        lstii(-1) = lstii(-1) + 1
        lstii(lstii(-1)) = npl
      endif

      lstii(-1) = lstii(-1) + 1 - iil

      if (iprint() >= 10) then
C        stdo = nglob('stdo')
        call ilst2a(lstii(0),lstii(-1),strn)
        call word(strn,1,i,iil)
        call info2(10,0,0,' pgglst mode %i %19p%i PL %26p('
     .    //strn(1:iil)//')',mode,lstii(-1))

      endif
      end
