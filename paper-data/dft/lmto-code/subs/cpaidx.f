      subroutine cpaidx(mode,s_spec,s_site,s_cpasite,
     .  nbas,nlst,iblst,idcc,nsite,nbcomp,mxcomp,ibcomp)
C- Table of indices for site list, including CPA sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: ncomp ncpa iscpa nthet
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec class
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_cpasite
Ci     Elts read:  *
Co     Stored:     ib idcc ncomp angles
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do nothing
Ci         :1 Include any normal class or CPA class with angle 0 or pi
Ci         :10s digit
Ci         :0 returns nsite,nbcomp only
Ci         :1 returns s_cpasite also.
Ci         :  It is the caller's responsibility to ensure that
Ci         :  s_cpasite must be dimensioned at least as large as nsite
Ci         :2 returns ibcomp also
Ci         :  ibcomp must be dimensioned as large as nbcomp
Ci         :1+2 digits can be combined
Ci         :100s digit
Ci         :0 no not enforce condition that angles be 0 or pi
Ci         :1 enforce condition that angles be 0 or pi
Ci   nbas  :size of basis
Ci   nlst  :0  => site list is 1,2,...nbas
Ci         :>0 => use site list iblst(1),iblst(2)...
Ci   iblst :Used only if nlst>0
Ci         :Site indices to include in list this routine generates
Ci         :Any site ib not in iblst is excluded.
Ci         :NB: iblst is assumed to be ordered.
Co   idcc  :points to links between CPA and parent classes (see dlmidx.f)
Co Outputs
Co  nsite  :number of sites in list:
Co         :nsite = nbas if nlst=0
Co         :nsite = nlst if nlst>0
Co  nbcomp :number of components in site list
Co         :Unless some members of the site list are CPA sites,
Co         :nbcomp is returned as nsite.
Co         :Otherwise, nbcomp is expanded: CPA species is has
Co         :s_spec(is)%ncomp components
Co  ibcomp :ibcomp(:,1:nbcomp) contains information about components
Ci         :ibcomp(1,k) = site index associated with component k
Ci         :ibcomp(2,k) = class index associated with this component
Ci         :ibcomp(3,k) = 0 if component is a normal site
Ci                      = 1 if component is a CPA site with angle 0
Ci                      = 2 if component is a CPA site with angle pi
Ci                       -1 otherwise
Ci         :ibcomp(4,k) = component index within site ib:
Ci                        ibcomp(4,k) = 1 for first component at ib,
Ci                                      2 for second, etc
Co  mxcomp :max number of components associated with any site
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   First CPA class associated with site ib: s_site(ib)%dlmcl
Cu Updates
Cu   16 Jan 13 Redesigned
Cu   14 Nov 12 First created
C ----------------------------------------------------------------------

      use pointers

      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nlst,nsite,nbcomp,mxcomp,ibcomp(4,nbcomp)
      integer iblst(nlst),idcc(*)
C ... Structures
!      include 'structures.h'

      type(str_spec):: s_spec(*)
      type(str_site):: s_site(nbas)
      type(str_cpasite):: s_cpasite(*)
C ... Local parameters
      integer maxcpa
      parameter (maxcpa=30)
      integer iscpa(maxcpa),icomp
      integer ilst,ib,is,ic,ncomp,mode0,mode1,mode2,ncpa,icpa,nthet,ith
      double precision pi,xx

C --- Setup
      pi = 4d0*datan(1d0)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      if (mode0 == 0) return
      call sanrg(.true.,mode0,1,1,'abc','def')

C --- Loop over all sites, excluing sites not in iblst
C     ilst = index to current site in iblst
      ilst = 0
      mxcomp = 0
      nbcomp = 0
      nsite = 0
      do  ib = 1, nbas
        if (nlst > 0) then
   11     if (ilst >= nlst) cycle
          if (iblst(ilst+1) < ib) then
            ilst = ilst+1
            goto 11
          endif
          if (iblst(ilst+1) /= ib) cycle
        else
          ilst = ib
        endif

        is = s_site(ib)%spec
        ic = s_site(ib)%class
C       ncomp = s_site(ib)%ncomp
        ncomp = s_spec(is)%ncomp
        nsite = nsite+1
        mxcomp = max(mxcomp,ncomp)
        if (mod(mode1,2) /= 0) then
          s_cpasite(nsite)%ib = ib
          s_cpasite(nsite)%idcc = ic
          s_cpasite(nsite)%ncomp = ncomp
        endif
C   ... Case normal class
        if (idcc(ic) == ic) then
          nbcomp = nbcomp+1
          if (mod(mode1/2,2) /= 0) then
            ibcomp(1,nbcomp) = ib
            ibcomp(2,nbcomp) = ic
            ibcomp(3,nbcomp) = 0
            ibcomp(4,nbcomp) = 0
          endif
C   ... Case CPA class
        else
          ncpa = s_spec(is)%ncpa
          if (ncpa == 0) then
            iscpa(1) = is       ! Just one CPA element -> true element
            ncpa = 1
          else
            iscpa(1:ncpa) = s_spec(is)%iscpa(1:ncpa)
          endif
          if (mod(mode1,2) /= 0) then
            s_cpasite(nsite)%ib = ib
            s_cpasite(nsite)%idcc = ic
            s_cpasite(nsite)%ncomp = ncomp
            allocate(p_d1(ncomp)); s_cpasite(nsite)%angles => p_d1
          endif
          icomp = 0
          do  icpa = 1, ncpa
            ic = icpa + idcc(ic)-1
            is = iscpa(icpa)
            if (is == 0) call rx('wrong species for CPA class')
            nthet = s_spec(is)%nthet
            if (nthet < 2) then
              nthet = 1
            endif
            do  ith = 1, nthet
              icomp = icomp+1
C             Is this correct?
              xx = s_site(ib)%thet(ith,1)
              if (mode2 == 1) then
                if (xx /= 0 .and. xx /= pi)
     .            call rx('cpaidx: only collinear angles allowed')
              endif
              if (mod(mode1,2) /= 0) then
                s_cpasite(nsite)%angles(icomp) = xx
              endif
C             Fill out ibcomp for components in this site
              if (mod(mode1/2,2) /= 0) then
                ibcomp(1,nbcomp+ith) = ib
                ibcomp(2,nbcomp+ith) = ic+ith-1
                ibcomp(4,nbcomp+ith) = icomp
                if (s_site(ib)%thet(ith,1) == 0) then
                  ibcomp(3,nbcomp+ith) = 1
                elseif (s_site(ib)%thet(ith,1) == pi) then
                  ibcomp(3,nbcomp+ith) = 2
                else
                  ibcomp(3,nbcomp+ith) = -1
                endif
              endif
            enddo
            nbcomp = nbcomp+nthet
          enddo
        endif
      enddo
      end
