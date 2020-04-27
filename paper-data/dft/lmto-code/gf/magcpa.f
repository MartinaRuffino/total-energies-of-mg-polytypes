      subroutine magcpa(s_ctrl,s_site,nl,sop,amag)
C- Magnetic and orbital moments, SOC energy terms for all sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lrel nbas
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lncol
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp norb thet cpawt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dmat dlmcl
Cio    Passed to:  *
Ci Inputs
Ci   nl    :dimensons sop
Ci   sop   :spin-orbit parameters (atomsr.f)
Co Outputs
Co   amag  :system magnetic moment
Cs Command-line switches
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   26 Nov 17 (Belashchenko) Added magnetic torques
Cu   08 Jun 14 (Belashchenko) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl
      double precision sop(0:nl-1,2,2,9,*),amag(3)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
C ... Local parameters
      integer nbas,ib,ncomp,icomp,norb
      real(8) eula(3),rotg(3,3),amloc(3),amg(3),orbl(3),orbg(3),trql(3),trqg(3),pi,f
      real(8), allocatable :: mag(:,:),orbm(:,:),socterms(:,:),torq(:,:)
      logical lso
      integer mpipid,procid,master
      procedure(logical) :: bittst

      pi = 4*datan(1d0)
      procid = mpipid(1)
      master = 0
      if (procid /= master) return
      lso = .false.
      if (mod(s_ctrl%lrel,10) == 2 .or. bittst(s_ctrl%lncol,4)) lso = .true.

      nbas = s_ctrl%nbas
      eula = 0d0
      amag = 0

C     write(lgunit(1),901)
      call info0(30,1,0,
     .  ' MAGCPA: spin and orbital moments, SOC (mRy) and torques (Ry) from density matrices'//
     .  '%N%24fglobal coordinates%18flocal coordinates'//
     .  '%N   ib  ic%9f x%10f y%10f z%12f x%10f y%10f z')

      do ib = 1, nbas
        ncomp = s_site(ib)%ncomp
        allocate(mag(3,ncomp),orbm(3,ncomp),socterms(5,ncomp),torq(2,ncomp))
        norb = s_site(ib)%norb
        call dmatav(lso,nl,norb,ncomp,s_site(ib)%dmat,
     .    sop(0,1,1,1,s_site(ib)%dlmcl),mag,orbm,socterms,torq)
        eula(1:3) = s_site(ib)%eula(1:3)
        do icomp = 1, ncomp
          if (s_site(ib)%thet(icomp,1) /= 0d0) then
            eula(1) = s_site(ib)%thet(icomp,2) 
            eula(2) = s_site(ib)%thet(icomp,1) ; eula(3) = - eula(1)
          endif
          call eua2rm(eula(1),eula(2),eula(3),rotg) 
C          call prmx('rotg',rotg,3,3,3)
C          call eua2rm(-eula(3),-eula(2),-eula(1),rotg) 
C          call prmx('rotg',rotg,3,3,3)
          amloc(:)  = mag(:,icomp) ; amg = matmul(rotg,amloc)
          amag = amag + s_site(ib)%cpawt(icomp) * amg/nbas
C         write (lgunit(1),902) ib,icomp,(amg(i),i=1,3),(amloc(i),i=1,3)
          call info5(30,0,0,'%,5i%,4i spin%3:2;10F  %3:2;10F',ib,icomp,amg,amloc,5)
C         call info5(30,0,0,'%,5i%,4i spin%3;12,8D  %3;12,8D',ib,icomp,amg,amloc,5)
          if (.not. lso) cycle
          orbl(:) = orbm(:,icomp) ; orbg = matmul(rotg,orbl)
          trql(3) = 0 ; trql(1:2) = torq(1:2,icomp) ; trqg = matmul(rotg,trql)
C         write (lgunit(1),903) (orbg(i),i=1,3),(orbl(i),i=1,3)
          call info2(30,0,0,'%9f orb %3;12,8D  %3;12,8D',orbg,orbl)
          call info2(30,0,0,'%9f Torq%3:2;10F  %3:2;10F',trqg,trql)
C         call info2(30,0,0,'%9f Torq%3;12,8D  %3;12,8D',trqg,trql)
C         write (lgunit(1),904) (1000*socterms(i,icomp),i=1,5)
          f = 1000
          call info5(30,0,0,'%9f SOC  11:%;10,6D 12:%;10,6D 21:%;10,6D 22:%;10,6D total:%;10,6D',
     .    f*socterms(1,icomp),f*socterms(2,icomp),f*socterms(3,icomp),f*socterms(4,icomp),f*socterms(5,icomp))
        enddo
        deallocate(mag,orbm,socterms,torq)
      enddo

C  901 format(/' MAGCPA: magnetic and orbital moments from density '
C     .  'matrices, includes CPA sites'/9x,'SOC energy terms are in mRy'/
C     .  21x,'global coordinates',15x,
C     .  'local coordinates'/'   ib  ic',9x,' x',8x,' y',8x,
C     .  ' z',10x,' x',8x,' y',8x,' z')
C  902 format(i5,i4,' spin',3F10.6,2x,3F10.6)
C  903 format(9x' orb ',3F10.6,2x,3F10.6)
C  904 format(9x' SOC, 11:',F10.6,' 12:',F10.6,' 21:',F10.6,' 22:',F10.6,
C     .  ' total:',F10.6)
      
      end

      subroutine dmatav(lso,nl,norb,ncomp,dmat,sop,mag,orbm,soc,torq)
C- Calculate spin and orbital moments, SOC energy terms, and magnetic torque for one site
Ci    lso = T: calculate orbital moments, SOC terms, and torques
      implicit none
      logical lso
      integer nl,norb,ncomp
      double complex dmat(norb,norb,2,2,4,ncomp)
      double precision sop(0:nl-1,2,2,9,ncomp),
     .  mag(3,ncomp),orbm(3,ncomp),soc(5,ncomp),torq(2,ncomp)

      double complex ci,dmxx(norb,norb),lpl,lmi,trprod
      double precision lp(norb,norb),lm(norb,norb),lz(norb,norb)
      integer i,j,jx,l(norb),icomp
      real(8) sopl(0:nl-1)

      ci = dcmplx(0,1)

C --- Magnetic moment (explicit calculation)
      mag = 0
      do i = 1, norb
        mag(1,:) = mag(1,:) + dmat(i,i,1,2,1,:) + dmat(i,i,2,1,1,:)
        mag(2,:) = mag(2,:) + ci*(dmat(i,i,1,2,1,:) - dmat(i,i,2,1,1,:))
        mag(3,:) = mag(3,:) + dmat(i,i,1,1,1,:) - dmat(i,i,2,2,1,:)
      enddo        
      if (.not. lso) return
C --- Orbital moments 
      orbm = 0
      call orbmop(norb,l,lz,lp,lm)
      do icomp = 1, ncomp
        dmxx(:,:) = dmat(:,:,1,1,1,icomp) + dmat(:,:,2,2,1,icomp)
        lpl = trprod(norb,lp,dmxx)
        lmi = trprod(norb,lm,dmxx)
        orbm(1,icomp) = (lpl + lmi)/2
        orbm(2,icomp) = (lpl - lmi)/2/ci
        orbm(3,icomp) = trprod(norb,lz,dmxx)
      enddo
C ... This follows the convention (spin 1 means -1/2, etc.)
c     orbm = - orbm   (removed)

C --- Spin-orbit coupling energy (4 terms)
      soc = 0
      do icomp = 1, ncomp
        do j = 1, 4
          jx = 2 ; if (j == 1) jx = 1 ; if (j == 4) jx = 3
          dmxx(:,:) = dmat(:,:,1,1,j,icomp)
          call dmscal(norb,l,sop(:,1,1,jx,icomp),dmxx)
c         soc(1,icomp) = soc(1,icomp) - trprod(norb,lz,dmxx)/2 ! spin #1 is s=-1/2, hence factor -1/2
          soc(1,icomp) = soc(1,icomp) + trprod(norb,lz,dmxx)/2 ! spin #1 is s=-1/2, hence factor -1/2
          dmxx(:,:) = dmat(:,:,2,2,j,icomp)
          call dmscal(norb,l,sop(:,2,2,jx,icomp),dmxx)
c         soc(4,icomp) = soc(4,icomp) + trprod(norb,lz,dmxx)/2 ! spin #2 is s=+1/2, hence factor 1/2
          soc(4,icomp) = soc(4,icomp) - trprod(norb,lz,dmxx)/2 ! spin #2 is s=+1/2, hence factor 1/2
          dmxx(:,:) = dmat(:,:,1,2,j,icomp)
          if (j == 1) sopl = sop(:,1,2,1,icomp)
          if (j == 4) sopl = sop(:,1,2,3,icomp)
          if (j == 2) sopl = sop(:,2,1,2,icomp)
          if (j == 3) sopl = sop(:,1,2,2,icomp)
          call dmscal(norb,l,sopl,dmxx)
c         soc(2,icomp) = soc(2,icomp) + trprod(norb,lm,dmxx)/2 ! dmat(1,2) goes with l-
          soc(2,icomp) = soc(2,icomp) + trprod(norb,lp,dmxx)/2 ! dmat(1,2) goes with l+
          dmxx(:,:) = dmat(:,:,2,1,j,icomp)
          if (j == 2) sopl = sop(:,1,2,2,icomp)
          if (j == 3) sopl = sop(:,2,1,2,icomp)
          call dmscal(norb,l,sopl,dmxx)
c         soc(3,icomp) = soc(3,icomp) + trprod(norb,lp,dmxx)/2 ! dmat(2,1) goes with l+
          soc(3,icomp) = soc(3,icomp) + trprod(norb,lm,dmxx)/2 ! dmat(2,1) goes with l-
        enddo
      enddo
      soc(5,:) = soc(1,:) + soc(2,:) + soc(3,:) + soc(4,:)

C --- Magnetic torque (2 components)
      torq = 0
      do icomp = 1, ncomp
        do j = 1, 4
          if (j == 1) sopl = sop(:,1,2,7,icomp)
          if (j == 4) sopl = sop(:,1,2,9,icomp)
          if (j == 2) sopl = sop(:,2,1,8,icomp)
          if (j == 3) sopl = sop(:,1,2,8,icomp)
          do i = 1, norb
c           torq(1,icomp) = torq(1,icomp) - ci*sopl(l(norb))*(dmat(i,i,1,2,j,icomp))
c           torq(2,icomp) = torq(2,icomp) +    sopl(l(norb))*(dmat(i,i,1,2,j,icomp))
            torq(1,icomp) = torq(1,icomp) + ci*sopl(l(i))*(dmat(i,i,1,2,j,icomp))
            torq(2,icomp) = torq(2,icomp) -    sopl(l(i))*(dmat(i,i,1,2,j,icomp))
          enddo        
          if (j == 2) sopl = sop(:,1,2,8,icomp)
          if (j == 3) sopl = sop(:,2,1,8,icomp)
          do i = 1, norb
c           torq(1,icomp) = torq(1,icomp) + ci*sopl(l(norb))*(dmat(i,i,2,1,j,icomp))
c           torq(2,icomp) = torq(2,icomp) +    sopl(l(norb))*(dmat(i,i,2,1,j,icomp))
            torq(1,icomp) = torq(1,icomp) - ci*sopl(l(i))*(dmat(i,i,2,1,j,icomp))
            torq(2,icomp) = torq(2,icomp) -    sopl(l(i))*(dmat(i,i,2,1,j,icomp))
          enddo        
        enddo
      enddo
      end

      double complex function trprod(n,a,b)
      implicit none
      integer n
      double precision a(n,n)
      double complex b(n,n)

      double complex prod(n,n)
      integer i

      trprod = 0
      prod = matmul(a,b)
      do i = 1, n
        trprod = trprod + prod(i,i) 
      enddo 
      end

      subroutine zmscal(norb,l,sop,dmat)
      implicit none
      integer norb,l(norb)
      complex(8) sop(0:*)
      complex(8) dmat(norb,norb)

      integer i,j
      do i = 1, norb
        do j = 1, norb
          if (l(i) .ne. l(j)) cycle
          dmat(i,j) = dmat(i,j) * sop(l(i))
        enddo
      enddo
      end

      subroutine dmscal(norb,l,sop,dmat)
      implicit none
      integer norb,l(norb)
      double precision sop(0:*)
      double complex dmat(norb,norb)

      integer i,j
      do i = 1, norb
        do j = 1, norb
          if (l(i) .ne. l(j)) cycle
          dmat(i,j) = dmat(i,j) * sop(l(i))
        enddo
      enddo
      end

      subroutine orbmop(norb,l,lz,lp,lm)
C-Set up orbital moment operator matrices
C-----------------------------------------------------------------------
Cu Updates
Cu   18 Jun 18 Synchronize with updated spherical harmonics
C-----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer norb,nl
      integer l(norb)
      double precision, dimension(norb,norb) :: lz,lp,lm
C ... Local parameters
      integer m(norb),i,j
      procedure(integer) :: ll

      nl = ll(norb)+1
      if (nl*nl /= norb) call rx1('orbmop: invalid value %i for norb',norb)
      call sanrg(.true.,nl,1,4,'orbmop:','nl')
      
C ... l and m for each orbital in order of appearance. opt=1 orders m as -l, -l+1, ..., l
!     call setmorderdefault(2)
      call lmorder(-1,nl-1,l,m)

C ... Orbital moment operators 
      lp = 0 ; lm = 0 ; lz = 0
      do i = 1, norb
        lz(i,i) = m(i)
        do j = 1, norb
          if (l(i) /= l(j)) cycle
          if (m(j) == m(i) - 1) then
C ... lp and lm are actually transposes of l+ and l- 
            lp(i,j) = sqrt(dble((l(i)+m(i))*(l(i)-m(j))))
            lm(j,i) = lp(i,j)
          endif
        enddo
      enddo

C     call prmx('lp',lp,norb,norb,norb)
C     call prmx('lm',lm,norb,norb,norb)
C     call prmx('lz',lz,norb,norb,norb)

      end
