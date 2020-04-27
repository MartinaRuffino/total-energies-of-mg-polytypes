C      module bzdata
C      implicit none
C      type str_bzdata
C
C      integer    ::   nqbz     ! Number of k-points in the full BZ
C      integer    ::   nqibz    ! Number of k-points in the irreducible part
C      integer    ::   nqbzw    ! Dimensions qbzw,ib1bz.  Should be (n1+1)*(n2+1)*(n3+1)
C      integer    ::   ntetf    ! Number of tetrahedra in the full BZ
C      integer    ::   nteti    ! Number of tetrahedra in the irreducible part
C      integer    ::   ngrp     ! Number of symmetry operations
C      integer    ::   n123(3)  ! Number of k-divisions along the reciprocal lattice vectors
C      integer,pointer::
C     .  idtetf(:,:),           ! corners of tetrahedra for full BZ
C     .  ib1bz(:),              ! maps k-points in the padded qp list (qbzw) to the original list
C     .  idteti(:,:),           ! corners of tetrahedra for irreducible BZ
C     .  nstar(:),              ! star of k
C     .  irk(:,:),              ! irreducible
C     .  nstbz(:)               !
C      real(8),pointer::
C     .  qbz(:,:),wbz(:),       ! k-points in the full Brillouin zone, and weights
C     .  qibz(:,:),wibz(:),     ! k-points in the irreducible Brillouin zone, and weights
C     .  qbzw(:,:)             ! k-points on the regular mesh, padded to repeat boundary points
C      end type str_bzdata
C      end module bzdata
      subroutine read_bzdata(job,s_bzdata)
C- Read header information from BZDATA, written by the GW code
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bzdata
Ci     Elts read:  all
Co     Stored:     all
Co     Allocated:  *
Cio    Elts passed:nqbz nqibz nqbzw ntetf nteti ngrp qbz wbz qibz wibz
Cio                nstbz nstar irk idtetf ib1bz qbzw idteti_r_r
Ci Inputs
Ci   job   : 1s digit not used now
Ci         : 10s digit
Ci         : 1 render s_bzdata%idtetf in standard Questaal s_bzdata%idteti, i.e.:
Ci         :   use regular (not padded) q mesh, and dimension it as (0:4,ntetf)
Cu Updates
Cu   28 Mar 17 First created
C ----------------------------------------------------------------------
      implicit none
C ... For structures
      type str_bzdata

      integer    ::   nqbz     ! Number of k-points in the full BZ
      integer    ::   nqibz    ! Number of k-points in the irreducible part
      integer    ::   nqbzw    ! Dimensions qbzw,ib1bz.  Should be (n1+1)*(n2+1)*(n3+1)
      integer    ::   ntetf    ! Number of tetrahedra in the full BZ
      integer    ::   nteti    ! Number of tetrahedra in the irreducible part
      integer    ::   ngrp     ! Number of symmetry operations
      integer    ::   n123(3)  ! Number of k-divisions along the reciprocal lattice vectors
      integer,pointer::
     .  idtetf(:,:),           ! corners of tetrahedra for full BZ
     .  ib1bz(:),              ! maps k-points in the padded qp list (qbzw) to the original list
     .  idteti(:,:),           ! corners of tetrahedra for irreducible BZ
     .  nstar(:),              ! star of k
     .  irk(:,:),              ! irreducible
     .  nstbz(:)               ! ? see genqbz
      real(8),pointer::
     .  qbz(:,:),wbz(:),       ! k-points in the full Brillouin zone, and weights
     .  qibz(:,:),wibz(:),     ! k-points in the irreducible Brillouin zone, and weights
     .  qbzw(:,:)              ! k-points on the regular mesh, padded to repeat boundary points
      end type str_bzdata
      type(str_bzdata)::  s_bzdata
C ... Passed parameters
      integer job
C ... Local parameters
      integer ifi,fopng,nqbz,ntetf,i,j
      real(8):: qbas(3,3),ginv(3,3),dq_(3)
      integer, pointer :: idtetf(:,:)

      ifi = fopng('BZDATA',-1,1)
      rewind ifi
      read(ifi,*)  s_bzdata%nqbz,s_bzdata%nqibz,s_bzdata%nqbzw,
     .             s_bzdata%ntetf,s_bzdata%nteti,s_bzdata%ngrp
      read(ifi,*)  s_bzdata%n123(1:3)
!     if (job == 0) return

      nqbz = s_bzdata%nqbz
      ntetf = s_bzdata%ntetf

      allocate(s_bzdata%qbz(3,nqbz),s_bzdata%wbz(nqbz))
      allocate(s_bzdata%qibz(3,nqbz),s_bzdata%wibz(nqbz),s_bzdata%nstbz(nqbz))
      allocate(s_bzdata%nstar(s_bzdata%nqibz),s_bzdata%irk(s_bzdata%nqibz,s_bzdata%ngrp))
      if (ntetf>0 .and. mod(job/10,10) == 1) then
        allocate(idtetf(4,ntetf))
        allocate(s_bzdata%idtetf(0:4,ntetf),s_bzdata%ib1bz(s_bzdata%nqbzw),s_bzdata%qbzw(3,s_bzdata%nqbzw))
        s_bzdata%idtetf(0,1:ntetf) = 1
      elseif (ntetf>0) then
        allocate(s_bzdata%idtetf(1:4,ntetf),s_bzdata%ib1bz(s_bzdata%nqbzw),s_bzdata%qbzw(3,s_bzdata%nqbzw))
        idtetf => s_bzdata%idtetf
      else
        allocate(s_bzdata%idtetf(1,1),s_bzdata%ib1bz(1),s_bzdata%qbzw(1,1))
      endif
      if (s_bzdata%nteti>0) then
        allocate(s_bzdata%idteti(0:4,6*nqbz))
      endif

      call rwbzdata(ifi,1,
     .  s_bzdata%ngrp,qbas,ginv,
     .  s_bzdata%qbz,s_bzdata%wbz,s_bzdata%nstbz,nqbz,
     .  s_bzdata%qibz,s_bzdata%wibz,s_bzdata%nstar,s_bzdata%irk,s_bzdata%nqibz,
     .  idtetf,ntetf,s_bzdata%qbzw,s_bzdata%ib1bz,s_bzdata%nqbzw,
     .  s_bzdata%idteti,s_bzdata%nteti,dq_)

      if (ntetf>0 .and. mod(job/10,10) == 1) then
        do  i = 1, ntetf
          forall (j = 1:4) s_bzdata%idtetf(j,i) = s_bzdata%ib1bz(idtetf(j,i))
        enddo
        deallocate(idtetf)
      endif

      end

      subroutine rwbzdata(ifbz,job,
C- Taken from the GW code
     .  ngrp,qbas,ginv,
     .  qbz,wbz,nstbz,nqbz,
     .  qibz,wibz,nstar,irk,nqibz,
     .  idtetf,ntetf,qbzw,ib1bz,nqbzw,
     .  idteti,nteti,dq_)
C- ReadWrite  BZ mesh data reuired for GW ---
C----------------------------------------------------------------
      implicit none
      integer:: nqbz,ntetf,nteti,nqbzw,iqbz,ifbz
     & ,nqibz,iqibz,itet,ngrp,job
      real(8):: qbz(3,nqbz),wbz(nqbz),qibz(3,nqibz),wibz(nqibz)
     &         ,qbzw(3,nqbzw),qbas(3,3),ginv(3,3),dq_(3)
      integer:: idtetf(0:3,ntetf),ib1bz(nqbzw),idteti(0:4,nteti)
     &       ,irk(nqibz,ngrp),nstar(nqibz),nstbz(nqbz)
      character(len=16) :: irk_fmt

      write(irk_fmt, "('(',i0,'(x,i8))')") ngrp

      if(job<0) write(ifbz,"(3d24.16)") qbas,ginv
      if(job>0) read (ifbz,"(3d24.16)") qbas,ginv

      do iqibz = 1,nqibz
        if(job<0) then
          write(ifbz,"(4d24.16,i9)")
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          write(ifbz,irk_fmt) irk(iqibz,1:ngrp)
        else
          read(ifbz,"(4d24.16,i9)")
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          read(ifbz,irk_fmt) irk(iqibz,1:ngrp)
        endif
      enddo

      do iqbz = 1,nqbz
        if(job<0) then
          write(ifbz,"(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
        else
          read(ifbz, "(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
        endif
      enddo

      if (ntetf>0) then
        if(job<0) then
          write(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          write(ifbz,"(i9,3d24.16)") (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        else
          read(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          read(ifbz,"(i9,3d24.16)") (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        endif
      endif

      if(nteti>0) then
        if(job<0) then
          write(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        else
          read(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        endif
      endif

      if(job<0) write(ifbz,"(3d24.16,' !dq_')") dq_
      if(job>0) read (ifbz,"(3d24.16)") dq_

      end subroutine rwbzdata
