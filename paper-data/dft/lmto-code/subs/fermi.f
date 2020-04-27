      subroutine fermi(qval,dos,ndos,emin,emax,nsp,eferm,e1,e2,dosef)
C- Makes fermi energy from integrated density
C ----------------------------------------------------------------------
Ci Inputs
Ci   qval, number of electrons to fermi level; dos(i) integrated
Ci   density at bin i; ndos, number of bins + 1; emin, emax, energy
Ci   window.
Co Outputs
Co   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
Co   i.e., Fermi energy lies between e1 and e2.
Co   dosef:  density of states at fermi level
Cr Remarks
Cr   emin and e1 (and emax and e2) may point to the same address.
Cr   This version uses idos decomposed into spin up, down for nsp=2
Cu Updates
Cu   29 Mar 11 Handle pathological case DOS doesn't encompass qval,
Cu             but is within fuzz of it
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer ndos,nsp
      double precision qval,dos(ndos,1),emin,emax,eferm,e1,e2,dosef
C Local parameters
      integer i1,ie,i2
      double precision de,q1,q2,d1mach,wt
      double precision fuzz
      parameter (fuzz=1d-8)

C --- Check bounds of DOS ---
      de = (emax-emin)/(ndos-1)
      wt = 1d0/2
      i2 = 1
      if (nsp == 2) then
        wt = 1
        i2 = 2
      endif
      q1 = (dos(1,1)+dos(1,i2))*wt
      q2 = (dos(ndos,1)+dos(ndos,i2))*wt
      if (q2 < q1) call rx('FERMI: q2<q1')

C      if (q1 > qval .or. q2 < qval) then
C        print *, 'hi'
C        q1 = 151.99999999800159161572810249985d0
C        q2 = 151.99999999900159161572810256213d0
C      endif

C     q1 > qval, but within fuzz of it
      if (q1 > qval .and. q1 < qval+fuzz) q1 = qval
C     q2 < qval, but within fuzz of it
      if (q2 < qval .and. q2 > qval-fuzz) q2 = qval

      if (q1 > qval .or. q2 < qval) call fexit3(-1,111,
     .  ' Exit -1 FERMI: NOS ( %1,6;6g %1,6;6g ) does not '//
     .  'encompass Q = %1;6g',q1,q2,qval)

C     If qval at an endpoint, skip bin search
      if (q2 == qval) then
        e1 = emax
        e2 = e1 + de
        eferm = emax
        dosef = (q2-q1)/de
        return
      endif
      if (q1 == qval) then
        e2 = emin
        e1 = e1 - de
        eferm = emin
        dosef = (q2-q1)/de
        return
      endif

C --- Find bin that boxes E_f ---
      i1 = 1
      q1 = (qval + d1mach(3))/wt
      do  ie = 2, ndos
        if (dos(ie,1)+dos(ie,i2) > q1) goto 10
        i1 = ie
      enddo
      call rx('bug in FERMI')

C --- Linear interpolation for the Fermi level ---
   10 continue
      q1 = (dos(i1,1)+dos(i1,i2))*wt
      q2 = (dos(i1+1,1)+dos(i1+1,i2))*wt
      e1 = emin + de*(i1-1)
      e2 = e1 + de
      eferm = e1 + (qval-q1)/(q2-q1)*de
      dosef = (q2-q1)/de

      end
