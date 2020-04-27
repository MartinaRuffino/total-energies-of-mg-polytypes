      subroutine dlmq(nclass,s_spec,ics,qt,vrmax,nth,wts)
C- Compute average charge for class ic
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   nth   :number of DLM angles
Ci   wts   :DLM weights
Cio Inputs/Outputs
Cio  qt    :qt(ic) = average charge in class ic
Cio  vrmax :average potential at wsr
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass,nth(nclass)
      integer ics(nclass)
      double precision qt(*),vrmax(2,*),wts(*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ic,ith,ioff
      double precision q,vrm(2)
      character*8 sname
      integer mpipid,procid,master,iprint
      integer stdo,nglob

      procid = mpipid(1)
      master = 0
      stdo = nglob('stdo')

      ioff = nclass
      do ic = 1, nclass
        sname = s_spec(ics(ic))%name
        if (nth(ic) < 2) cycle
        q = 0d0
        vrm(1) = 0d0
        vrm(2) = 0d0
        do ith = 1,nth(ic)
          q = q + qt(ioff+ith)*wts(ioff+ith)
          vrm(1) = vrm(1) + vrmax(1,ioff+ith)*wts(ioff+ith)
          vrm(2) = vrm(2) + vrmax(2,ioff+ith)*wts(ioff+ith)
        enddo
        ioff = ioff + nth(ic)
        qt(ic) = q
        vrmax(1,ic) = vrm(1)
        vrmax(2,ic) = vrm(2)
        if (iprint() >= 40 .and. procid == master) then
        write(stdo,901) sname,q
 901    format(1X,'DLMQ:   ',A8,' has average charge ',F9.6)
        endif
      enddo

      end

