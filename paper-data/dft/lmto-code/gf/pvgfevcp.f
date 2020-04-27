      subroutine pvgfevcp(s_site,mode,ib,icomp,pfun,ldpf,norb,norbx,sfvrtx,gamma1,gamma2,pfloc)
C- Vertex associated with a single site for exchange interactions
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 Potential functions read from s_site(ib)%pfr, stored in vector form, pfr(norb,1:2)
Ci         :1 Potential functions read from s_site(ib)%pfr, stored in matrix form, pfr(norb,2,norb,2)
Ci         :2 Potential functions read from pfun, stored in vector form pfun(:,1:2)
Ci         :10s digit
Ci         :0 Return interaction gamma1(norb) in vector form; gamma2 not used
Ci         :1 Return interaction gammmacp(norb,norb) in matrix form; gamma1 not used
Ci         :100s digit
Ci         :0 Generate gamma2 from pfr
Ci         :1 Generate gamma2 CPA 12 vertex
Ci         :  Note: in this mode the remaining options in mode are ignored!
Ci         :2 Same as 1
Ci         :  Note: in this mode the remaining options in mode are ignored!
Ci         :1000s digit
Ci         :0 Do not rotate pfr
Ci         :1 rotate pfr from (mu,lambda,lambda') representation to lms representation
Ci         :2 rotate pfr from lms representation to (mu,lambda,lambda') representation
Ci   icomp :component to CPA vertex. icomp=0 if not CPA
Ci   ib    :site containing CPA vertex
Ci   norb  :dimensions s_site(ib)%pfr and sfvrtx
Ci   norbx :dimensions
Ci   sfvrtx:spin flip vertex for all components in site ib
Co Outputs
Ci  gamma1 :(Non CPA sites only) Pfr(up)-Pfr(dn), as a vector of dimension norb
Ci  gamma2 :Non CPA sites :      Pfr(up)-Pfr(dn), as (norb,norb) matrix
Ci         :CPA sites: vertex for given component
Cl Local variables
Cl         :
Cr Remarks
Cr  *This routine uses EITHER
Cr     potential functions stored in s_site%pfr (icomp<1)
Cr     a CPA vertex stored in sfvrtx (icomp>0)
Cu Updates
Cu   04 Sep 15 Re-designed for a general case
Cu   13 Jul 15 (Vishina) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,icomp,ib,norb,norbx,ldpf
      double complex sfvrtx(norb,norb,2,*)
      double complex gamma2(norbx,norbx),gamma1(norb),pfun(ldpf,2),pfloc(norb,2,norb,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      integer, allocatable :: idxkm(:,:,:)
C     complex(8), pointer :: pf
C ... Local parameters
      integer i,j,mode2,mode3,nl,ll
      double complex wk(norb,2,norb,2),wk2(norb,2)

C     mode0 = mod(mode,10)
C     mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)

C --- Handle CPA vertex as a special case ---
      if (mode2 == 1 .or. mode2 == 2) then
        do  i = 1, norb
          do  j = 1, norb
            gamma2(i,j) = sfvrtx(i,j,mode2,icomp)
          enddo
        enddo
C       call zprm('gamma2',2,gamma2,norbx,norb,norb)
        return
      endif

C ... Copy pfr into matrix pfloc applicable to all modes below
C     Read potential functions from pfun
      if (mod(mode,10) == 2) then
        call rxx(mode3/=0,'mode 3 not allowed with mode0=2')
        call dpzero(pfloc,2*size(pfloc))
        do  i = 1, norb
          pfloc(i,1,i,1) = pfun(i,1)
          pfloc(i,2,i,2) = pfun(i,2)
        enddo
C     Read potential functions from s_site(ib)%pfr, vector form
      elseif (mod(mode,10) == 0) then
        call rxx(mode3/=0,'mode 3 not allowed with mode0=2')
        call dpzero(pfloc,2*size(pfloc))
        call zcopy(norb*2,s_site(ib)%pfr,1,wk2,1)
        do  i = 1, norb
          pfloc(i,1,i,1) = wk2(i,1)
          pfloc(i,2,i,2) = wk2(i,2)
        enddo
C     Read potential functions from s_site(ib)%pfr, matrix form
      else ! 1s digit mode =  1
C        call yprm0('(1p,9e18.10)')
C        call yprmi('Pfr for ib=%i',ib,0,3,s_site(ib)%pfr,0,norb*2,norb*2,norb*2)

        call zcopy(norb*norb*4,s_site(ib)%pfr,1,pfloc,1)
C       Rotate (mu,lambda,lambda') repsn to lms repsn or vice-versa
        if (mode3 /= 0) then
          nl = ll(norb)+1
          allocate(idxkm(2,0:nl-1,2*(nl+1)))
          call mstokm(2-mode3,nl,1,norb,wk,pfloc,idxkm)
          call zcopy(size(pfloc),wk,1,pfloc,1)
          deallocate(idxkm)
        endif
      endif

C      call yprm0('(1p,9e18.10)')
C      call yprmi('Pfloc for ib=%i',ib,0,3,pfloc,0,norb*2,norb*2,norb*2)

C ... Return site interaction as vector gamma1
      if (mod(mode/10,10) == 0) then
        do  i = 1, norb
          gamma1(i) = pfloc(i,1,i,1) - pfloc(i,2,i,2)
        enddo
C       call zprm('gamma1',2,gamma1,norbx,norb,1)

C ... Return site interaction as matrix gamma2
      else
        do  i = 1, norb
          do  j = 1, norb
            gamma2(i,j) = pfloc(i,1,j,1) - pfloc(i,2,j,2)
          enddo
        enddo
C       call zprm('gamma2',2,gamma2,norbx,norb,norb)
      endif

      end

C      subroutine pvgfevcp(s_site,mode,icomp,ib,ldpf,pfun,norb,norbx,nl,nbas,lrel,nspc,indxsh,
C     .            sfvrtx,gamiicp,gamiirel,mskm)
CC- Vertex associated with a single site for exchange interactions
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1s digit
CCi         :0 Not used here
CCi         :1 fill gamii with pfun(1)-pfun(2)
CCi         :10s digit
CCi         :0 make no contribution to gamii from a CPA vertex
CCi         :1 Copy diagonal part of CPA 12 vertex within component icomp, site ib
CCi         :  to appropriate place in gamii (nothing done if icomp<=0)
CCi         :2 Copy diagonal part of CPA 21 vertex within component icomp, site ib
CCi         :  to appropriate place in gamii (nothing done if icomp<=0)
CCi   icomp :component to CPA vertex. icomp=0 if not CPA
CCi   ib    :site containing CPA vertex
CCi   ldpf  :dimensions potential functions P
CCi   pfun  :vector of potential functions
CCi   norb  :dimensions sfvrtx
CCi   sfvrtx:spin flip vertex for all components in site ib
CCo Outputs
CCi   gamii :Pfun(up)-Pfun(dn) in the absence of vertices,
CCi         :CPA case: diagonal part of (spin) vertex is substituted
CCi         :at site ib for given component
CCl Local variables
CCl         :
CCr Remarks
CCr    For sites other than CPA sites, it is assumed that s_site%pfr has been filled
CCu Updates
CCu   02 Sep 15 (Vishina) extended to relativistic case
CCu   13 Jul 15 (Vishina) first created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,icomp,ib,ldpf,norb,norbx,nspc,lrel,nl,indxsh(*),nbas
C      double complex pfun(ldpf,2),sfvrtx(norb,norb,2,*)
C      double complex gamiicp(norbx,norbx)
C!      double complex gamiirel(norbx*2,norbx*2)
C      double complex gamiirel(ldpf)
CC ... For structures
C      include 'structures.h'
C      type(str_site)::  s_site(*)
CC ... Local parameters
C      logical ldsite
C      double complex wk(norb,nspc,norb,nspc),wk2(norb,2)
C      double complex wkkmu(norb*2,norb*2),wkms(ldpf,2,2)
C      double complex fun(nspc,nspc),xx
C      double complex wk3(norb,norb,nspc,nspc),wk4(norb,norb,nspc,nspc)
C      double complex mskm(0:nl-1,2*nl)
C      integer idx(2,0:nl-1,2*(nl+1))
C      integer i,j,mode1,ncomp,l,imu,i1,i2,pfdim,ibas,ib1,ib2
C      integer m11,m12,m21,m22,lmr11,lmr22,lmrd11,lmrd22,idxn,io
C
CC     mode0 = mod(mode,10)
C      ncomp = s_site(ib)%ncomp
C      mode1 = mod(mode/10,10)
C      pfdim = ldpf
C
C      if (icomp < 1 .and. nspc == 1) then    ! Non CPA branch, nonrelativistic.  P is a diagonal matrix
C        gamiicp = 0
C        call zcopy(norb*2,s_site(ib)%pfr,1,wk2,1)
C        do  i = 1, norb
C          gamiicp(i,i) = wk2(i,1) - wk2(i,2)
C        enddo
C
CC ... non-CPA branch
C      elseif (icomp < 1) then ! Non CPA branch
C
C        if (lrel==2) then
C
C
C
C
C          gamiirel(:) = 0
C          call zcopy(norb*norb*4,s_site(ib)%pfr,1,wk,1)
C
C          print *, wk(1,1,1,1),wk(2,1,2,1),wk(3,1,3,1),wk(4,1,4,1),wk(5,1,5,1)
C          stop
C
C          call zcopy(ldpf*4,s_site(ib)%pfr,1,wkms,1)
C          do i = 1, ldpf
C            gamiirel(i) = wkms(i,1,1) - wkms(i,2,2)
C          enddo
C        else                    !lrel ne 2
C          if (mode1 == 1 .or. mode1 == 2) then
CC          print*,'pvgfevcp mode1 1 or 2'
C            call zcopy(norb*norb*4,s_site(ib)%pfr,1,wk,1)
C            do  i = 1, norb
C              do  j = 1, norb
C                gamiicp(i,j) = wk(i,1,j,1) - wk(i,2,j,2)
C              enddo
C            enddo
C          endif
C        endif                   ! lrel eq 2
CC ... CPA branch
C      else
C        if (mode1 == 1 .or. mode1 == 2) then
C          do  i = 1, norb
C            do  j = 1, norb
C              gamiicp(i,j) = sfvrtx(i,j,mode1,icomp)
C            enddo
C          enddo
C        endif
C      endif
C
C      end
