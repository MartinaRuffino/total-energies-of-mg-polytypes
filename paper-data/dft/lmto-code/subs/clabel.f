      subroutine uniqclabl(mode,nclass,nbas,ips,first,slabl,clabl)
C- Uniquifies a class or species label
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 loop over nclass species identified by slabl(1:nclass)
Ci         :  In this mode nclass is required, neither nbas nor ips is used
Ci         :1 loop over nbas species identified by slabl(ips(1:nbas))
Ci         :  In this mode nbas and ips are required
Ci   nclass:number of inequivalent classes.  Not used if mode=1.
Ci   nbas  :size of basis.  Not used if mode=0
Ci   ips   :species table: site ib belongs to species ips(ib). Not used if mode=0
Ci   slabl :vector of species labels or classes
Ci   first :first integer to append
Cio Inpus/Outputs
Cio  clabl :class name
Cio        :  On input, a name which may or may not be unique
Cio        :  On output, name is make unique by appending a number
Cr Remarks
Cr
Cu Updates
Cu   12 Nov 17 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nclass,ips(nbas),first
      character*8 slabl(*),clabl
C ... Local parameters
      integer n,ib,j
      character*8 clabl0

      do  n = first-1, 99999    ! Loop until unique label found
        clabl0 = clabl          ! Original label
        if (n >= first) call awrit1('%a%i',clabl0,8,0,n) ! append a digit
C       Check for uniqueness
        if (mode == 0) then
          do  ib = 1, nclass
            j = ib
            if (slabl(j) == clabl0) exit ! Match found
          enddo
        else
          do  ib = 1, nbas      ! Loop over sites prior to insert instruction
            j = ips(ib)
            if (slabl(j) == clabl0) exit ! Match found
          enddo
        endif
C       At this point, j is either a matching label or the last member of the list
        if (slabl(j) /= clabl0) exit ! If no match, we are done
      enddo
      clabl = clabl0

      end

      subroutine clabel0(slabl,is,idx,clabl)
C- Make class label from species label
C ----------------------------------------------------------------------
Ci Inputs
Ci   slabl :vector of species labels
Ci   is    :current species for which to make class label
Ci   idx   :idx=1 (first occurence of species): assign clabl=slabl(is)
Ci         :idx>1 class label is species label, with idx appended.
Co Outputs
Co   clabl:    class label
Cb Bugs
Cb   No check is made to ensure uniqueness of clabl.  See uniqclabl
Cu Updates
Cu   05 Apr 18 subroutine renamed from 'subroutine clabel'
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer is,idx
      character*8 slabl(is),clabl

      clabl = slabl(is)
      if (idx > 1) call awrit1('%a%i',clabl,8,0,idx)

      end

      subroutine dlmclbl(nclass,nccomp,s_spec,ics,dclabl)
C- Class label for extra DLM atom files
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: ncomp
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   dclabl:class name, packed as a real number
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   14 Mar 13 Improved printout
Cu   28 May 12 (Kirill) extensions for disordered local moments
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass,nccomp
      integer ics(nclass+nccomp)
      double precision dclabl(*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      character*8 clabl,alabl
      integer ic,ncomp,j,icd

      if (nccomp == 0) return

      call info2(20,0,0,' Generate class labels for %i CPA components',
     .  nccomp,0)
      call info0(31,0,0,' Parent%10pClass%20pSpecies')

      icd = nclass
      do  ic = 1, nclass
        call r8tos8(dclabl(ic),clabl)
        ncomp = s_spec(ics(ic))%ncomp
        if (ncomp > 1) then
C          call info2(20,0,0,'%10pGenerating'//
C     .      ' %i class labels for class '//clabl,ncomp,0)
          do  j = 1, ncomp
            icd = icd + 1
            alabl = clabl
            call awrit1('%a#%i',alabl,8,0,j)
            call s8tor8(alabl,dclabl(icd))
            call info2(31,0,0,' '//clabl//
     .        '%(n>9?9:10)p%i:'//alabl//'   %i',icd,ics(icd))
          enddo
        endif
      enddo

      end
