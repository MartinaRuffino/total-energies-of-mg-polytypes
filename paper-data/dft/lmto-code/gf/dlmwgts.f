      subroutine dlmwgts(nclass,ncomp,s_spec,ics,lso,angle,wts,entropy)
C- Make angles and weights for DLM classes
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name ncpa nthet iscpa beff
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:xcpa beff
Cio    Passed to:  *
Ci Inputs
Ci   nclass:number of inequivalent classes
Ci   ncomp :total number of CPA components
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   lso   :.true. for fully-relativistic or spin-orbit coupling
Co Outputs
Co   angle :DLM angles
Co         :For each species,
Co         : nthet  # of angles  angles
Co         :   <2       1          0
Co         :    2       2         0,pi
Co         :   >2       nthet     points on gaussian quadrature
Co   wts   :corresponding weights for angles
Co  entropy:alloy+spin contribution to entropy
Cl Local variables
Cl   icp   :index to current component
Cr Remarks
Cr   angle includes all components in each class ic:
Ci   For each class there is a loop over components;
Ci   the number of components is given by spec(ics(ic))%ncpa
Ci   spec(ics(ic))%ncpa = 0 => no chemical CPA => 1 component
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   13 Oct 12 Some modifications for new structures
Cu   01 Nov 11 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclass,ncomp,ics(*)
      double precision wts(*),angle(ncomp,2),entropy
      logical lso
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: w1(:),w2(:),w3(:),angfil(:,:)
C ... Local parameters
      integer maxth,maxcpa,ntbig
      parameter (maxth=1000,maxcpa=30,ntbig=1000)
      integer ncpa,ic,is,nthet,nang(2),icp,icpa,ith,iscpa(maxcpa),i
      integer iprint,lgunit,stdo,ifi,fxst,fopn
      double precision xcpa(maxcpa),wk1(maxth),wk2(maxth),xx(3,maxth),
     .  pi,bf,s,wt,wtq,dval
      logical lpr,lsphere
      character*8 sname
      integer mpipid,procid,master

C --- Setup ---
      procid = mpipid(1)
      master = 0
      stdo = lgunit(1)
      entropy = 0
      pi = 4*datan(1d0)
      lpr = .false.
      if (procid == master .and. iprint() >= 30) then
        write(stdo,900)
  900   format(' DLMWGTS: Angles and weights for CPA classes:')
        lpr = .true.
      endif

C --- Read angles file for polar angles ---
      allocate(angfil(3,nclass+ncomp)) ; angfil = 0
      ifi = 0
      if (procid == master) then
        if (fxst('angles') == 1) then
          ifi = fopn('angles')
          rewind ifi
C         call info0(20,0,0,'          Found angles file, reading')
C         call iomagf(2,nclass+ncomp,1,angfil,angfil,1,ifi)
          call ioextf(2,'angles',nclass+ncomp,1,angfil,angfil,3,1,ifi)
          call fclr('angles',ifi)
        endif
      endif
      call mpibc1(ifi,1,2,.false.,'dlmwgts','ifi')
      if (ifi /= 0) call mpibc1(angfil,
     .  (nclass+ncomp)*3,4,.false.,'dlmwgts','shfac')

      angle = 0d0
      icp = 0  ! Running index counting CPA classes

C --- For each class do ---
      do  ic = 1, nclass
        is = ics(ic) ; wts(ic) = 1d0
        sname = s_spec(is)%name
        ncpa = s_spec(is)%ncpa
C ...   If not chemical CPA, check for DLM
        if (ncpa == 0) then
          nthet = s_spec(is)%nthet
          if (nthet < 2) cycle  ! No DLM; nothing to do
          iscpa(1) = is ! Just one CPA element -> true element
          xcpa(1) = 1d0 ; ncpa = 1
        else
          iscpa(1:ncpa) = s_spec(is)%iscpa(1:ncpa)
          call dcopy(ncpa,s_spec(is)%xcpa,1,xcpa,1)
          if (abs(sum(xcpa(1:ncpa))-1) > 1d-8) then
            call fexit(-1,111,' Exit -1 dlmwgts: concentrations do not'
     .        //' add up to 1 for species',is)
          endif
        endif

C  ...  Non-CPA sites have ncpa = 1 here
        do  icpa = 1, ncpa
          if(iscpa(icpa) == 0) exit
          nthet = s_spec(iscpa(icpa))%nthet
          if (nthet > maxth)
     .      call fexit(-1,111,' Exit -1 dlmwgts: increase maxth',maxth)
          nang(1:2) = s_spec(iscpa(icpa))%nang(1:2)
          lsphere = nang(1) /= 0
          if (.not. lsphere .and. nthet > 2 .and. lso)
     .      call rx('NANG is required for vector DLM sites for SOC')
C     ... Explicit integration over a sphere for lso=.true.
          if (lsphere) then
            if (.not. lso) call rx('NANG requires lrel=2 or SO')
            call fpiint(nang(1),nang(2),nthet,xx,wk2)
            if (nthet /= s_spec(iscpa(icpa))%nthet)
     .        call rx2('dlmwgts: wrong nthet=%i for class %i',nthet,ic)
            is = ics(nclass+icp+1)
            if (lpr) write(stdo,902) sname,(s_spec(is)%beff(i),i=1,10)
            wk1(1:nthet) = xx(3,1:nthet)
            call gibbswts(nthet,wk1,s_spec(is)%beff,wk2)
            angle(icp+1:icp+nthet,1) = dacos(xx(3,1:nthet))
            angle(icp+1:icp+nthet,2) =atan2(xx(2,1:nthet),xx(1,1:nthet))
            wts(nclass+icp+1:nclass+icp+nthet) = wk2(1:nthet)*xcpa(icpa)
            if (lpr) then
              do  ith = 1, nthet
                write(stdo,901) sname,icp+ith,angle(icp+ith,1),
     .            angle(icp+ith,2),wts(nclass+icp+ith)
              enddo
            endif
            icp = icp + nthet ; wts(ic) = 0
C     ... No DLM for this CPA component - assign angles from file if available
          elseif (nthet < 2) then
            angle(icp+1,1) = angfil(1,nclass+icp+1)
            wts(nclass+icp+1) = xcpa(icpa)
            if (lpr) write(stdo,901)sname,icp+1,angle(icp+1,1),
     .        wts(nclass+icp+1)
            icp = icp + 1 ; wts(ic) = 0
C     ... Paramagnetic DLM for this component (spin up and down)
          elseif (nthet == 2) then
            angle(icp+1,1) = 0 ; angle(icp+2,1) = pi
            bf = s_spec(ics(nclass+icp+1))%beff(1)
            wk2(1) = dexp(-bf) ; wk2(2) = dexp(bf)
            wk2(1) = wk2(1)/(wk2(1)+wk2(2)) ; wk2(2) = 1d0 - wk2(1)
            wts(nclass+icp+1) = wk2(1)*xcpa(icpa)
            wts(nclass+icp+2) = wk2(2)*xcpa(icpa)
            if (lpr) then
             write(stdo,901)sname,icp+1,angle(icp+1,1),wts(nclass+icp+1)
             write(stdo,901)sname,icp+2,angle(icp+2,1),wts(nclass+icp+2)
            endif
            icp = icp + 2 ; wts(ic) = 0d0
C     ... Axially symmetric DLM for this component (Gaussian quadrature)
          elseif (nthet > 2) then
            call gausq(nthet,-1d0,1d0,wk1,wk2,0,0)
            is = ics(nclass+icp+1)
            if (lpr) write(stdo,902)sname,(s_spec(is)%beff(i),i=1,10)
            wk1(1:nthet) = - wk1(1:nthet)
            wk2(1:nthet) = wk2(1:nthet)/2
            call gibbswts(nthet,wk1,s_spec(is)%beff,wk2)
            do  ith = 1, nthet
              angle(icp+ith,1) = dacos(wk1(ith))
              wts(nclass+icp+ith) = wk2(ith)*xcpa(icpa)
              if (lpr) then
                write(stdo,901)sname,icp+ith,angle(icp+ith,1),
     .            wts(nclass+icp+ith)
              endif
            enddo
            icp = icp + nthet ; wts(ic) = 0
          endif
C     ... Contribution to entropy from this component
          if (lsphere .or. nthet > 2) then
            allocate(w1(ntbig),w2(ntbig),w3(ntbig))
            call gausq(ntbig,-1d0,1d0,w1,w2,0,0)
            call dpscop(w2,w3,ntbig,1,1,0.5d0)
            call gibbswts(ntbig,w1,s_spec(is)%beff,w2)
            s = 0
            do  i = 1, ntbig
              wt  = dval(w2,i) ; wtq = dval(w3,i)
              s = s - wt * dlog(wt/wtq)
            enddo
            entropy = entropy + xcpa(icpa) * s
            deallocate(w1,w2,w3)
C ...       end entropy part
          endif
        enddo
      enddo

      deallocate(angfil)
  901 format(10X,A8,'#',i3,1X,3F11.6)
  902 format(' BEFF for ',A8,10F8.3)

      end
