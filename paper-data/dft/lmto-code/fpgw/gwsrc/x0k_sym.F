      subroutine x0k_sym(npm,q,ginv,nbloch,nwt,ngbb,ngc,nolfco,zzr,nmbas,zcousq,chipmzzr,ngrp,eibzsym,
     o  rcxq)
C- Symmetrize polarization matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci  npm    : 1 if time-reversal symmetry, 2 if not
Ci  q      : Symmetrize chi0(q)
Ci  ginv   : Inverse of reciprocal lattice vectors (?)
Ci  nbloch : number of product basis functions B
Ci           RENAME to ngprd
Ci  nwt    : number of energy points for which rcxq is given
Ci  ngbb   : Leading dimension of zcousq
Ci         : RENAME to ldz
Ci  ngc    : number of G vectors in mixed basis
Ci         : RENAME to ngprd
Ci  nolfco : If T => no local fields
Ci  zzr    : <eiqr|B_I> ... used for magnetic susceptibility calculations (hvccfp0.m.F)
Ci         : Used for magnetic susceptibility case, nolfco = T
Ci  nmbas  : Number of mixed basis functions : rcxq = rcxq(nmbas,nmbas,:,:)
Ci  zcousq : Eigenvectors of the diagonalized bare coulomb matrix
Ci         : Used to rotate rcxq from M~ basis to E- basis
Ci  chipmzzr : If T => rcxq is magnetic susceptibility
Ci  ngrp   : Number of point group operations
Ci  eibzsym: ? see subroutine eibzgen
Co Outputs
Co  rcxq   : polarization matrix is symmetrized, and rotated from M~ to E- basis
Cl Local variables
Cr Remarks
Cr  Symmetrizer of EIBZ PRB.81,125102(2010) Eq.(51) july2012takao
Cr  This may be not so effective ---> only for limited cases?
Cr  zrotm(J,J') = <M~^k_J| \hat{A}^k_i M~^k_J'>
Cr  We do \sum_i T_alpha_i [ zrotm_i^dagger (I,I') P_I'J' zrotm_i(J'J) ]
Cr
Cr  Other comments from Takao (?)
Cr  (exactly speaking, we insert conversion matrix between Enu basis and M_I basis).
Cr  input qin = q
Cr  \hat{A}^k_i  is specified by symops(:,:,igx),and igxt (-1 for time-reversal).
Cr  Note that k= \hat{A}^k_i(k) (S_A^k)
Cr  See Eq.(51) around in PRB81 125102(2010)
C ----------------------------------------------------------------------
      use keyvalue,only:  getkeyvalue
      use m_rotMPB,only:  rotMPB2
      use m_readqgcou,only: qtt_,nqnum
      implicit none
C ... Passed parameters
      logical :: nolfco,chipmzzr
      integer :: npm,ngbb,nwt,nbloch,jpm,iw,ngb,ngc,nmbas
      real(8):: q(3),ginv(3,3)
      complex (8):: rcxq(nmbas,nmbas,nwt,npm),zcousq(ngbb,ngbb)
C ... Local parameters
      logical,parameter:: debug=.false.
      integer :: verbose
      integer :: ix,iy
      real(8):: quu(3)
      complex(8),target :: zzr(ngbb,nmbas)
      integer:: ngrp
      complex(8),allocatable:: rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
      integer:: eibzsym(ngrp,-1:1),neibz,icc,ig,ikp,itimer,icount
      integer:: irotm,nrotmx,itt,ntimer, nccc,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:),zrrx(:,:),rcxq_core(:,:)

      integer :: i,j
      integer(8) :: cr,cm,c1,c2
      call system_clock(count_rate = cr, count_max = cm)

      if (ngrp == 0) return   ! Neither symmetrize nor convert to E-mu basis

      call system_clock(c1)

      ngb = nbloch + ngc        ! This is not ngbb for smbasis()=T. oct2005
      if (ngb/=ngbb) call rx('x0k_sym: ngb/=ngbb')
      if (ngrp == 1) goto 100

!     print *, 'x0k_sym !!'; goto 100

      call iqindx2(q,ginv,qtt_,nqnum,ikp,quu) !to get ikp for timereversal mode
      if (sum(abs(q-quu))>1d-8) call rx( 'x0kf_v2h: eibz 111 q/quu')
      neibz = sum(eibzsym(:,1))+sum(eibzsym(:,-1))
                                !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
      write(6,"(' x0k_sym: start symmetrization  ikp neibz q=',2i3,3f12.8)") ikp,neibz,q
      ntimer=1
      if (sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case
      allocate(zrotm(ngb,ngb),nrotm(ngrp*2))

C  --- Assemble rotation matrx zrr,zrrc ---
C      Rotation matrix zrrx is typically a very sparse matrix.
C      Thus it is stored to  "i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,icc),nrotm(icc)".
C      In the following: matmul(rcxqwww,zrrx) is given by
C         do irotm1 = 1,nrotm(icc)
C           rcxq0(:,i2(irotm1,icc)) = rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,i2(irotm1,icc))
      allocate(zrrx(nmbas,nmbas))
      nrotmx = 10000            !trial value
      do 1011                   !this loop only finds large enough nrotmx.
        if (allocated(i1)) deallocate(i1,i2,zrr,zrrc) !,zrr_,zrrc_)
        nccc=ngrp*2
        allocate(i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,nccc),zrrc(nrotmx,nccc))
        i1=-99999
        i2=-99999
        zrr=-99999d0
        zrrc=-99999d0
        icc=0
        do itimer=1,-1,-2
          if (ntimer==1.and.itimer==-1) exit
          if (itimer==1 ) itt=1
          if (itimer==-1) itt=2
          do ig=1,ngrp
            if (eibzsym(ig,itimer)==1) then
              icc=icc+1

C         ... Get rotation matrix zrrx, which can be a sparse matrix
C             zrotm(J,J') is the rotation matrix = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>
C             See rotMPB2 defined in readeigen.F.
              call rotMPB2(nbloch,ngb,q,ig,itimer,ginv,zrotm)

              if (nolfco.and.chipmzzr) then
C               Assume <svec_I | svec_J >= \delta_IJ, In addition, we use fact that we have no IPW parts in svec.
C               If IPW part exist, we may have to take into account <IPW|IPW> matrix, e.g. as in ppovlz.
C               svec --> zzr
                if (itimer==1) then
                  zrrx= matmul(transpose(dconjg(zzr)), matmul(zrotm, zzr))
!                 call zprm('zzr',zzr,ngbb,ngbb,nmbas)
                else
                  zrrx= matmul(transpose(zzr), matmul(dconjg(zrotm), zzr))
                endif
C              elseif (nolfco) then
C                call rx( 'x0kf_v4h: this case is not implemented')

              else
C               zrrx(mu nu)= dconjg(Zcousq(I, mu)) *zrotm(I,J)* Zcousq(J, nu)
C               zrrx is very sparse matrix. Size is \sim ngb or something.

c$$$                if (itimer==1) then
c$$$                  call matmmsparse(zcousqinv,zrotm,zcousq,zrrx,ngb,1d-8,iele)
c$$$                  ! this means zrrx= matmul(zcousqinv,matmul(zrotm, zcousq))
c$$$                else
c$$$                  call matmmsparse(dconjg(zcousqinv),dconjg(zrotm),zcousq,zrrx,ngb,1d-8,iele)
c$$$                  ! this means zrrx= matmul(dconjg(zcousqinv),matmul(dconjg(zrotm), zcousq))
c$$$                endif

                if (itimer==1) then
                  zrrx=zrotm
c                  call matmmsparse(zcousqinv,zrotm,zcousq,zrrx,ngb,1d-8,iele)
                  ! this means zrrx= matmul(zcousqinv,matmul(zrotm, zcousq))
                else
                  zrrx=dconjg(zrotm)
c                  call matmmsparse(dconjg(zcousqinv),dconjg(zrotm),zcousq,zrrx,ngb,1d-8,iele)
                  ! this means zrrx= matmul(dconjg(zcousqinv),matmul(dconjg(zrotm), zcousq))
                endif

              endif
              i1(:,icc) = 0
              i2(:,icc) = 0
              irotm = 0
              iagain = 0
              do ix = 1,ngb
                do iy = 1,ngb
                  if (abs(zrrx(ix,iy))>1d-8) then
                    irotm = irotm+1
                    if (irotm>nrotmx) then
                      iagain = 1
                    endif
                    if (iagain/=1) then
                      i1(irotm,icc) = ix
                      i2(irotm,icc) = iy
                      zrr(irotm,icc) = zrrx(ix,iy)
                      zrrc(irotm,icc) = dconjg(zrr(irotm,icc))
                    endif
                  endif
                enddo
              enddo
              if (iagain==1) then
                nrotmx=irotm    !enlarge allocation and do things again.
                print *, 'x0k_sym warning (slow speed): repeat loop 1011'
                goto 1011       !enlarge nrotmx and try again.
              endif
              nrotm(icc)=irotm
C             print *,'ig itimer icc nrotm=',ig,itimer,icc,nrotm(icc) ,iele
            endif
          enddo
        enddo
        exit
 1011 continue                  !only when nrotmx overflow.

         call system_clock(c2)
         write(*,'("x0k_sym 1:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

         call system_clock(c1)
C  --- Main part : obtain symmetrized rcxq ---
C     neibz is total number of symmetrization operation.
C     rcxq is rotated and accumulated; finally divided by neibz
!     zcousqc = dconjg(transpose(zcousq))
!$OMP parallel private(rcxq000,icc,itt,icount,rcxqwww,rcxq00,rcxq0,rcxq_core)
      allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb),rcxq_core(ngb,ngb))
!$OMP do
      do iw=1,nwt
        do jpm=1,npm
!           rcxq000 = 0d0
          icc=0
          do itimer=1,-1,-2
            if (itimer==1 ) itt=1
            if (itimer==-1) itt=2
            icount=0
            if (itimer==1) then
              call zcopy(ngb*ngb,rcxq(1,1,iw,jpm),1,rcxqwww,1)
            else
              forall(i=1:ngb,j=1:ngb) rcxqwww(i,j) = rcxq(i,j,iw,jpm)
            endif
            rcxq00 = 0d0
            do ig=1,ngrp
              if (eibzsym(ig,itimer)==1) then
                icount=icount+1
                icc=icc+1
                rcxq0 =0d0

C               The following is equivalent to
C               rcxq00= rcxq00 + matmul(zrrc_(:,:,icc),matmul(rcxqwww,zrr_(:,:,icc)))
                do irotm1 = 1,nrotm(icc)
c                 if (abs(zrr(irotm1,icc))<1d-8) cycle
!                   rcxq0(:,i2(irotm1,icc)) = rcxq0(:,i2(irotm1,icc)) + rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,icc)
                  call zaxpy(ngb,zrr(irotm1,icc),rcxqwww(1,i1(irotm1,icc)),1,rcxq0(1,i2(irotm1,icc)),1)
                enddo
                do irotm2 = 1,nrotm(icc)
c                 if (abs(zrrc(irotm2,icc))<1d-8) cycle
!                   rcxq00(i2(irotm2,icc),:) = rcxq00(i2(irotm2,icc),:) + zrrc(irotm2,icc) * rcxq0(i1(irotm2,icc),:)
                  call zaxpy(ngb,zrrc(irotm2,icc),rcxq0(i1(irotm2,icc),1),ngb,rcxq00(i2(irotm2,icc),1),ngb)
                enddo
              endif
            enddo

            if (itimer==1) then
!               rcxq000=rcxq00
              call zcopy(ngb*ngb,rcxq00,1,rcxq000,1)
            else
!               rcxq000=rcxq000+rcxq00
              call zaxpy(ngb*ngb,(1d0,0d0),rcxq00,1,rcxq000,1)
            endif
          enddo
          rcxq(:,:,iw,jpm) = rcxq000/neibz
        enddo
      enddo
!$OMP end  do
      deallocate(rcxq00,rcxq000,rcxq0,rcxqwww)
!$OMP end parallel
      call system_clock(c2)
      write(*,'("x0k_sym 2:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      deallocate(zrotm,i1,i2)
      if (verbose()>=40) call cputid2(' complete symmetrization',0)

      call system_clock(c1)

C     Entry point for conversion to E-mu basis
  100 continue

!     call zprm('rcxq(iw=20) before rotation',rcxq(1,1,20,1),nmbas,nmbas,nmbas)

      call rcxq2Ebas(0,nmbas,nmbas,nwt,npm,ngb,zcousq,rcxq)

!     call zprm('rcxq(iw=20) after rotation',rcxq(1,1,20,1),nmbas,nmbas,nmbas)

      if (verbose()>=40) call cputid2(' complete conversion to E-mu basis',0)

      write(*,"(' --- x0k_sym: end')")
      end subroutine x0k_sym

      subroutine rcxq2Ebas(mode,nmbas1,nmbas,nwt,npm,ngb,zcousq,rcxq)
C- Rotates rcxq from M~ to E basis or the inverse
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 forward rotation
Ci         :1 reverse rotation
Ci   nmbas1:leading dimension of zcousq
Ci   nmbas :second dimension of zcousq
Ci   nwt   :number of energy points for which rcxq is given
Ci   npm   :1 if time-reversal symmetry, 2 if not
Ci   ngb   :dimensions zcousq; number of product basis functions
Ci   zcousq:(mode=0) eigenvectors of the coulomb matrix ; see x0kr.F
Ci         :(mode=1) Use (inverse)+ of zcousq, aka ppovlz
Cio Inputs/Outputs
Cio  rcxq  : rxcq is rotated : rxcq -> zcousq+ rcxq zcousq
Cr Remarks
Cr   zcousq+ ppovlz = ppovlz+ zcousq = 1
Cr   See Remarks in x0kr.F
Cu Updates
Cu   21 Nov 17 Extracted from x0k_sym
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: nmbas1,nmbas,nwt,npm,ngb,mode
      real(8):: q(3),ginv(3,3)
      complex (8):: rcxq(nmbas1,nmbas,nwt,npm),zcousq(ngb,ngb)
C ... Dynamically allocated arrays
      complex(8),allocatable:: rcxq0(:,:)
C ... Local parameters
      integer iw,jpm

      if (nmbas1 /= nmbas) call rx('rcxq2Ebas not ready for nmbas1 /= nmbas')

C!$OMP parallel private(rcxq0)
      allocate(rcxq0(nmbas,nmbas))
C!$OMP do
      do iw=1,nwt
        do jpm=1,npm
          if (mode == 0) then
            call zgemm('N','N',nmbas,nmbas,nmbas,(1d0,0d0),rcxq(1,1,iw,jpm),nmbas,zcousq,ngb,(0d0,0d0),rcxq0,nmbas)
            call zgemm('C','N',nmbas,nmbas,nmbas,(1d0,0d0),zcousq,ngb,rcxq0,nmbas,(0d0,0d0),rcxq(1,1,iw,jpm),nmbas)
          else
            call zgemm('N','C',nmbas,nmbas,nmbas,(1d0,0d0),rcxq(1,1,iw,jpm),nmbas,zcousq,ngb,(0d0,0d0),rcxq0,nmbas)
            call zgemm('N','N',nmbas,nmbas,nmbas,(1d0,0d0),zcousq,ngb,rcxq0,nmbas,(0d0,0d0),rcxq(1,1,iw,jpm),nmbas)
          endif
        enddo
      enddo

C!$OMP end  do
      deallocate(rcxq0)
C!$OMP end parallel

      end
