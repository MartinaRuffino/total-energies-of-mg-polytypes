    subroutine rho2frce(tbc,lmol,lgamma,plat,nbas,nl,nsp,isp,      &
       lmx,ipc,indxsh,iwk,ldim,bk,nsite,iax,npr,typesz,h0,dh,charge, &
       rl,trh,pv,force,f,thrpv,e,rho,rholm,rhoc,zr)

      use mpi
      use tbprl, only : tbc_t
      use mod_ctx, only : cart_ctx_t

      implicit none
!C Passed parameters
      type(tbc_t), intent(in), target :: tbc
      integer, intent(in) :: nbas,nl,nsp,isp,ldim,nsite,typesz,iwk(nbas), &
                             iax(0:9,nsite),npr(0:1,nbas),lmx(*),ipc(*),indxsh(*)
      logical, intent(in) :: lmol,lgamma,charge,rl,trh,pv,force


      real(8), intent(in) :: plat(3,3),bk(3),             &
                             h0(nl**4*nsite*nsp),         &
                             dh(nl**4*nsite,4),           &
                             zr(typesz,ldim,ldim)

      real(8), intent(inout) :: thrpv, f(3,nbas), e(nbas),   &
                                rho(0:nl-1,nsp,nbas),       &
                                rholm(nl**2,nsp,nbas),      &
                                rhoc(0:nl**2-1,0:nl**2-1,nbas)

      integer :: l,m,i,ib,ixa,ixb,nlma,nlmb,la,lb,lma,lmb,indxH(nbas),ia,ipa,nrange,nend,ldz
      real(8) :: eloc, floc(3)
      real(8), allocatable :: hk0(:,:,:), dhk(:,:,:,:), sk(:,:,:), dsk(:,:,:,:)
      logical ::  sqt !, walked(nbas)
!       logical, parameter :: naive =.true.
      logical, parameter :: naive =.false.
      real(8) :: rescale
      integer :: tmp(nbas), pid
      real(8) :: rtmp(1)

!       real(8) :: wk(ldim,ldim,2),wk2(ldim,ldim,2),wkcf(nbas)

      integer, parameter :: lix(0:8,0:8) = reshape([ &
     &     0,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 1
     &    -1,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 2
     &     1, 2, 3,-1,-1,-1,-1,-1,-1,            & ! 3
     &     0, 1, 2, 3,-1,-1,-1,-1,-1,            & ! 4
     &     4, 5, 6, 7, 8,-1,-1,-1,-1,            & ! 5
     &     0, 4, 5, 6, 7, 8,-1,-1,-1,            & ! 6
     &    -1,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 7
     &     1, 2, 3, 4, 5, 6, 7, 8,-1,            & ! 8
     &     0, 1, 2, 3, 4, 5, 6, 7, 8 ], [9,9])     ! 9

      integer, parameter :: llx(0:8) = [((l,m = -l,l),l=0,2)] ! [0, 1, 1, 1, 2, 2, 2, 2, 2]

      integer :: nl2
      real(8) :: rm(typesz,0:nl*nl-1,0:nl*nl-1)
      real(8), dimension(0:nl*nl-1,0:nl*nl-1) :: rmc
      real(8) :: rv(0:nl*nl-1), rc(0:nl-1)

      real(8) ::  nul(2) = [0.0_8, 0.0_8], one(2) = [1.0_8, 0.0_8]

      logical :: lroot, offsite
      type(cart_ctx_t), pointer :: c2d


      procedure() :: s_ddot, dgemm, dgemv, s_zdotc, zgemm, zgemv
      procedure(), pointer :: gemm => null(), gemv => null(), s_dot => null()

      character :: mt

      integer :: lsz, err, mtype, u

      call tcn('tbfrce')

      if (lgamma) then
         mt = 't'
         gemm => dgemm
         gemv => dgemv
      else
         mt = 'c'
         gemm => zgemm
         gemv => zgemv
      end if


      lroot  = tbc % c3d % lrt
      pid    = tbc % c2d % id

      c2d => tbc % c2d

!       lz = tbc % globsz(1)

      nl2 = nl*nl

      mtype = mpi_real8
      if (.not. lgamma) mtype = mpi_complex16

      call tcn('Traces')


!       if (naive .and. .not. aprox) call dgemm('t','n',ldim, ldim, nstate, 1.0_8, zt, ldim, zt, ldim, 0.0_8, zr, ldim)
!       if (naive) call gemm(mt,'n',ldim, ldim, nstate, one, z, lz, zt, lz, nul, zr, lz)



      indxH(1) = 0
      do  ib = 2, nbas
         indxH(ib) = indxH(ib-1) + iwk(ib-1)
      enddo


      if (trh) allocate(hk0(typesz,ldim,ldim))

      if (force .and. pv) then
         allocate(dhk(typesz,ldim,ldim,4))
      else if (force .and. (.not. pv)) then
         allocate(dhk(typesz,ldim,ldim,3))
      else if ((.not. force) .and. pv) then
         allocate(dhk(typesz,ldim,ldim,4:4))
      end if




      if (trh .or. force .or. pv) then
         call tcn('blt')
         call pshprt(0)

         if (trh) call ztbloch(lgamma,bk,nl,1,nsp,isp,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,h0,rtmp,rtmp,.false., &
                                                                                          & ldim,hk0,ldim,typesz)
! --- spin-independent: ispu = 1
         if (force) then
            do i = 1, 3
               call ztbloch(lgamma,bk,nl,1,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,dh(1,i),rtmp,rtmp,.false., &
                                                                                    & ldim,dhk(1,1,1,i),ldim,typesz)
            end do
         end if
         if (pv) call ztbloch(lgamma,bk,nl,1,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,dh(1,4),rtmp,rtmp,.false., &
                                                                                    & ldim,dhk(1,1,1,4),ldim,typesz)
         call popprt()
         call tcx('blt')
      end if


      offsite = trh .or. force .or. pv

      call tcn('walk')

      do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)

!          nrange = 1
!          if (offsite) nrange = npr(0,ia)
         nend = tbc % neighmap(ia)+1
         if (offsite) nend = tbc % neighmap(ia+1)

!          walked = .false.

         if (trh)   eloc = 0.0_8
         if (force) floc = 0.0_8

         do ipa = tbc % neighmap(ia)+1, nend
            ib = tbc % neighidx(ipa)

!          do ipa = 1, nrange
!             ib = iax(1, npr(1, ia)+ipa)
!             if (walked(ib)) cycle

            ixa = indxH(ia)+1
            ixb = indxH(ib)+1

            nlma = iwk(ia)
            nlmb = iwk(ib)

            rm(1:typesz, 0:nlma-1,0:nlmb-1) = zr(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1)

            if (ia == ib .and. (rl .or. charge)) then
               if (ia /= ib) call rxx(.true., 'tbfrce: npr <--> iax mismatch, quite likely a bug ..')
               if (rl) then
                  rmc = 0.0_8
                  do lb = 0, nlma - 1
                     lmb = lix(lb, nlma-1)
                     do la = 0, nlma - 1
                        rmc(lix(la, nlma-1), lmb) = rm(1,la,lb)
                     end do
                  end do
                  rhoc(0:nl2-1, 0:nl2-1, ia) = rhoc(0:nl2-1, 0:nl2-1, ia) + rmc
               end if

               if (charge) then
                  rv = 0.0_8
                  rc = 0.0_8
                  do la = 0, nlma - 1
                     lma = lix(la, nlma-1)
                     rv(lma) = rm(1,la,la)
                     rc(llx(lma)) = rc(llx(lma)) + rm(1,la,la)
                  end do
                  rho(0:nl-1, isp, ia) = rho(0:nl-1, isp, ia) + rc
                  rholm(1:nl2, isp, ia) = rholm(1:nl2, isp, ia) + rv
               end if
            end if

            if (trh) eloc = eloc + sum(hk0(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1) * rm(1:typesz,0:nlma-1,0:nlmb-1))

            if (force) then
               do i = 1, 3
                  floc(i) = floc(i) + sum(dhk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,i) * rm(1:typesz,0:nlma-1,0:nlmb-1))
               end do
            end if

            if (pv) thrpv = thrpv - sum(dhk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,4) * rm(1:typesz,0:nlma-1,0:nlmb-1))
!             if (pv) then
!                write(980,'(/,2(x,i5))') ia, ib
!                do lb = 0, nlmb-1
!                   do la = 0, nlma-1
!                      write(980,'(f18.12)',advance='no') dhk(1:typesz,ixa+la,ixb+lb,4)
!                   end do
!                   write(980,'()')
!                end do
!             end if
!             walked(ib) = .true.
         end do

         if (trh) e(ia) = e(ia) + eloc
         if (force) f(1:3,ia) = f(1:3,ia) - 2*floc
      end do

      call tcx('walk')

      if (force .or. pv) deallocate(dhk)
      if (trh) deallocate(hk0)

      call tcx('Traces')
      call tcx('tbfrce')


!       write(426,'(2(9(9(x,f16.8),/),/))') rhoc

      end subroutine rho2frce



