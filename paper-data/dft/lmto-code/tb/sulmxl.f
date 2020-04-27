      subroutine sulmxl(force,pv,point,nl,nclass,qpol,lmxl,nlmq,nlmq1)
C- Adjust lmxl and set dimensions for multipole-related arrays
C ----------------------------------------------------------------------
Ci Inputs
Ci   force  :T for calculating forces
Ci   pv     :T for calculating pressure
Ci   point  :T if only point charges shoud be used
Ci   nl     :lmax+1 for basis set
Ci   nclass :number of different species (classes?)
Ci   qpol   :polarisation parameters
Cio Inputs/Outputs
Cio  lmxl   :lmax for multipole charges Q_L
Co Outputs
Co   nlmq   :max L for multipole charges Q_L, all species (for dimensioning)
Co   nlmq1  :if force or pv, nlmq1 = (ll(nlmq)+1)^2, otherwise nlmq1 = nlmq
Cr Remarks
Cr   The program checks if the input lmxl can be decreased by comparing them
Cr   to the max l for which polarisation parameters are not zeroes.
Cr   Dimensioning parameter nlmq corresponds to maximal lmxl for all species
Cr   after the above adjustement. For forces (and provisionally pressure)
Cr   one needs to go one up in angular momentum. Parameter nlmq1 takes care
Cr   of this.
Cu Updates
Cu    25 Feb 2010 (SL) first created
C ----------------------------------------------------------------------
      implicit none

C ... Passed parameters
      logical, intent(in) :: force, pv, point
      integer, intent(in) :: nl,nclass
      double precision, intent(in) :: qpol(10,nclass)
      integer, intent(inout) :: lmxl(nclass)
      integer, intent(out) :: nlmq,nlmq1
C ... Local parameters
      integer il,il1,il2,ic,lmax,lmxl0,nlq,ll
      double precision M

C ...   Check lmxl against qpol and reduce if possible
      if (.not. point) then
        lmax = nl-1
        do ic = 1, nclass
          lmxl0 = 0
          do il1 = 0, lmax
            do il2 = 0, lmax
              do il = 0, 2*lmax
                call getM((il+1)**2,(il1+1)**2,(il2+1)**2,
     .                     qpol(1,ic),M)
                if (dabs(M) > 1d-8) lmxl0 = max(lmxl0,il)
              enddo
            enddo
          enddo
          if (lmxl0 < lmxl(ic)) then
            call info5(10,0,0,' SULMXL: lmxl for species %i'//
     .      ' is reduced from %i to %i',ic, lmxl(ic),lmxl0,0,0)
            lmxl(ic) = lmxl0
          endif
        enddo

C ... find largest lmxl --> nlq (for dimensioning)
        nlq = 0
        do ic = 1, nclass
          nlq = max(nlq,lmxl(ic))
        enddo
        nlq = nlq + 1
      else
C ... if 'point charges' option is on
        lmxl(1:nclass) = 0
        nlq = 1
      endif

C ... need lmax+1 for forces and pressure (provisionally)
      nlmq = nlq*nlq
      if (force .or. pv) then
        nlmq1 = (nlq+1)**2
      else
        nlmq1 = nlmq
      endif

C ... Just to make sure
      if(ll(nlmq) > 2*(nl-1)) then
        call info2(10,1,0,' sulmxl: ll(nlmq) exceeds 2*(nl-1),'//
     .    ' ll(nlmq) = %i, nl = %i',ll(nlmq),nl)
        call rx0(' sulmxl: something is wrong!')
      endif

      end


