      subroutine dlmidx(nclass,s_spec,ics,nccomp,ndlmcl,idcc,ncomp,nrcp)
C- Set up indices for CPA
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  ncpa nthet iscpa
Co     Stored:     ncomp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nclass:number of inequivalent classes
Cio Inputs/Outputs
Cio  nccomp:total number of extra CPA classes for all sites (see Remarks)
Cio        :Input nccomp=0 => init mode, routine makes:
Cio        :  nccomp ndlmcl ncomp(:) s_spec(:)%ncomp idcc(1:nclass)
Cio        :Input nccomp>0 => normal mode, routine also makes:
Cio        :  ics(extended) idcc(extended) nrcp
Cio  ics   :species table
Cio        :On input, class ic belongs to species ics(ic), ic=1:nclass
Cio        :nccomp extra classes are tacked to the end of ics
Cio        :On output, ics updated for CPA-extended class list
Cio        :ics(i>ic) index to parent species.  nclass<i<=nclass+nccomp
Co Outputs
Co   ndlmcl:number of classes subject to CPA/DLM treatment
Co   idcc  :points to links between CPA and parent classes
Co         :idcc(ic=1:nclass) points to child classes:
Co         :If class ic is NOT a CPA class, idcc(ic) points to itself
Co         :If ic IS a CPA class, idcc(ic) points to the first of
Co         :of the CPA or DLM classes associated with ic.
Co         :idcc(i=nclass+1:nclass+nccomp) (CPA classes):
Co         :index to the parent class for this CPA class
Co   ncomp :ncomp(ic) = total #  CPA components for class ic (ic=1:nclass)
Co         :This number includes number of chemical species, and
Co         :number of angles within each chemical species.
Co         :ncomp(ic)=0 => no disorder for that site
Co   nrcp  :nrcp(i) = number of sites belonging to (extended) class i,
Co         :i=1:nclass+nccomp
Cl Local variables
Ci   first :T initialization mode, where dimensions are determined
Cl   ncomp :Local version of ncomp(ic)
Cr Remarks
Cr   In the CPA case, classes are extended by nccomp classes.
Cr   Each element treated with chemical or spin disorder is assigned
Cr   new classes; thus a binary alloy with one atom/cell has 3 classes.
Cr   The total number of classes is nclass+nccomp.
Cr   The extended class list looks like:
Cr     1 Original class (some CPA average of the next two)
Cr     2 class for first component (true atom)
Cr     3 class for second component (true atom)
Cu Updates
Cu   14 Mar 13 correct treatment of classes; improved printout
Cu   25 Apr 12 (K Belashchenko) additions for CPA
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass,ndlmcl,nccomp
      integer idcc(nclass+nccomp),nrcp(nclass+nccomp),ics(*),
     .  ncomp(nclass)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical first
      integer lgunit,i1,i2,maxcpa
      parameter (maxcpa=30)
      character*(256) strn,strn2,strn3
      integer is,ic,nthet,icd,i,j,ncpa,ncompi,iscpa(maxcpa)

      integer mpipid,procid,master
      procid = mpipid(1)
      master = 0

      first = .true.
      if (nccomp /= 0) first = .false.

      nccomp = 0
      ndlmcl = 0
      call ivset(ncomp,1,nclass,0)
      icd = nclass
      do  ic = 1, nclass
        is = ics(ic)
        ncompi = 0
        ncpa = s_spec(is)%ncpa
C ...   If not chemical CPA, check if DLM
        if (ncpa == 0) then
          nthet = s_spec(is)%nthet
          if (nthet < 2) nthet = 0
          if (.not. first .and. nthet /= 0)
     .      call ivset(ics,icd+1,icd+nthet,is)
          ncompi = ncompi + nthet
          icd = icd + nthet
C ...   Chemical CPA class
        else
          iscpa(1:ncpa) = s_spec(is)%iscpa(1:ncpa)
          do  j = 1, ncpa
            if (iscpa(j) == 0) exit
            nthet = s_spec(iscpa(j))%nthet
            if (nthet < 2) nthet = 1
            if (.not. first) call ivset(ics,icd+1,icd+nthet,iscpa(j))
            ncompi = ncompi + nthet
            icd = icd + nthet
          enddo
        endif
        ncomp(ic) = ncompi
        s_spec(is)%ncomp = ncompi
        if (ncompi >= 1) then
          idcc(ic) = nclass + 1 + nccomp
          nccomp = nccomp + ncompi
          ndlmcl = ndlmcl + 1
        else
          idcc(ic) = ic
        endif
      enddo

      if (first) return

      do  ic = 1, nclass
        icd = idcc(ic)
        if (icd /= ic) then
          do  i = 0, ncomp(ic)-1
            idcc(icd + i) = ic
            nrcp(icd + i) = nrcp(ic)
          enddo
        else
          idcc(ic) = ic
        endif
      enddo

      if (procid == master) then
        strn = ' CPA class mapping:'
        strn2 = ' Sites per class:'
        strn3 = ' Species:'
        do  ic = 1, nclass
          if (idcc(ic) == ic) then
            call awrit1('%a  %i',strn,len(strn),0,idcc(ic))
          else
            call awrit3('%a  %i->%i:%i',strn,len(strn),0,ic,idcc(ic),
     .        idcc(ic)+ncomp(ic)-1)
          endif
          call wordg(strn,1,' ',ic+3,i1,i2)
          call awrit2('%np%,2i',strn2,len(strn2),0,i1-2,nrcp(ic))
          call awrit2('%np%,2i',strn3,len(strn3),0,i1-2,ics(ic))
        enddo
        call info0(30,0,0,trim(strn))
        call info0(30,0,0,trim(strn2))
        call info0(30,0,0,trim(strn3))
        write(lgunit(2),506) trim(strn)
        write(lgunit(2),506) trim(strn2)
        write(lgunit(2),506) trim(strn3)
  506   format(a)

C        do  i = 1, 2
C          if (iprint() >= 20 .or. i == 2) then
C          write(lgunit(i),500)
C          write(lgunit(i),505) trim(strn)
C          write(lgunit(i),501)(idcc(ic),ic=1,nclass)
C          write(lgunit(i),502)(idcc(ic),ic=nclass+1,nclass+nccomp)
C          write(lgunit(i),503)(nrcp(ic),ic=1,nclass+nccomp)
C          write(lgunit(i),504)(ics(ic),ic=1,nclass+nccomp)
C          endif
C        enddo
      endif

C  500 format(' DLMIDX: index arrays for DLM')
C  501 format(' Parent class -> first DLM class:',50i4)
C  502 format(' Child class  -> parent class:   ',50i4)
C  503 format(' Debug: ## of sites per class:   ',50i4)
C  504 format(' Debug: Class -> species:        ',50i4)
C  505 format(' CPA class mapping:              ',a)

      end
