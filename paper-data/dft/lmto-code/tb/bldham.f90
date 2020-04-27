   subroutine bldham(nsp,nspc,s_ctrl,s_lat,tbc,lscale,fitpar,nvar,&
      & ip1,ivar,nlmesp,memodk,decay,decov,poly,cutmod,cut,cutov,iam,  &
      & npm,nterm,nset,tabme,tabov,nsites,npr,iax,lmx,ham,dhm,ovl,dov)

      use tbprl
      use structures, only : str_ctrl, str_lat
      use mpi, only : mpi_in_place, mpi_real8, mpi_comm_world

      implicit none

      integer, parameter :: dp = 8

      type(str_ctrl), intent(in) ::  s_ctrl
      type(str_lat), intent(in) ::   s_lat
      type(tbc_t), intent(inout) :: tbc
      logical, intent(in) :: fitpar,lscale
      integer, intent(in) :: nsp, nspc, nvar, nlmesp, nset, memodk(nset), poly(*), cutmod(*), nterm, nsites, lmx(*)
      integer, intent(in) :: ip1(*), ivar(*), iam(3, *), npm(2, *), npr(0:1, *), iax(0:9, nsites)
      real(dp), intent(in) :: decay(nlmesp,*), cut(2,nlmesp,*), cutov(2, nlmesp, *), &
                            & decov(nlmesp,*), tabme(nterm, nlmesp, *), tabov(nterm, nlmesp, *)
      real(dp), intent(out) ::  ham(s_ctrl%nl**2, s_ctrl%nl**2, nsites, nspc**2   ), &
                              & ovl(s_ctrl%nl**2, s_ctrl%nl**2, nsites, nspc**2   ), &
                              & dhm(s_ctrl%nl**2, s_ctrl%nl**2, nsites, nspc**2, 4), &
                              & dov(s_ctrl%nl**2, s_ctrl%nl**2, nsites, nspc**2, 4)
!       real(dp), allocatable ::  ham(:,:,:,:), dhm(:,:,:,:,:), ovl(:,:,:,:), dov(:,:,:,:,:)

      real(dp) :: dcs(0:3)
      real(dp), allocatable :: v(:), dv(:)

      integer :: nl, nl2, ca, cb, ia, ib, isite, rsite, isp1, isp2, ilme, tblid, spidx, k, ltb, crd, d, memode
      real(dp) :: alat
      logical :: hdiffs, cryf, onsite, lov, ocryf
      integer :: comm, err
      integer, allocatable :: idx(:), sitemap(:), displs(:), recvcnt(:)
      procedure(logical) :: bittst

!       real(dp), pointer :: pham(:,:), povl(:,:)

      integer, parameter :: itab(2,2) = reshape([0,2,2,1],[2,2])

      call tcn('newham')

      comm = tbc % c3d % comm

      ltb = s_ctrl%ltb
      lov = bittst(ltb,1)
      cryf = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      hdiffs = bittst(ltb,16) .or. bittst(ltb,128)

      if (cryf .or. ocryf) call rx('CRYF removed')

      nl = s_ctrl % nl
      nl2 = nl*nl
      alat = s_lat % alat

      ilme = (nl*(nl+1)*(nl+2))/6
      allocate(v(ilme))

      d = 0

      if (hdiffs) then
         d = 1
         allocate(dv(ilme))
      else
         allocate(dv(1)) ! So compiler doesn't complain
      end if

! comment until the beginning of the nspc loops to stop any parallelism.
! Also change do isite = sitemap(aproc)+1, sitemap(aproc+1) to do isite = 1, nsites
! !
!       call tcn('pham prep')
!       allocate(idx(0:nsites))
!       idx(0) = 0
!       do isite = 1, nsites
!          ia = iax(0, isite)
!          ib = iax(1, isite)
!          ca = s_ctrl % ipc(ia)
!          cb = s_ctrl % ipc(ib)
!          idx(isite) = idx(isite-1) + ((lmx(ca)+1)*(lmx(cb)+1))**2
!       end do
!
!
!       allocate(sitemap(0:naproc), displs(0:naproc-1), recvcnt(0:naproc-1))
!       call vbdist(nsites, idx(1:nsites), naproc, sitemap)
!       deallocate(idx)
!       displs = nl2*nl2*sitemap(0:naproc-1)
!       recvcnt = nl2*nl2*(sitemap(1:naproc) - sitemap(0:naproc-1))
!
!       call tcx('pham prep')

      do isp2 = 1, nspc
         do isp1 = 1, nspc

            tblid = 1 + itab(isp1,isp2)*ilme
            spidx = (isp2-1)*nspc + isp1

            call tcn('ham etc..')
            do isite = 1, nsites
!             do isite = sitemap(aproc)+1, sitemap(aproc+1)

               ia = iax(0, isite)
               ib = iax(1, isite)
               onsite = isite == npr(1, ia) + 1
               if (onsite) cycle ! may be do something more here later and get rid of tbdiag?

               ca = s_ctrl % ipc(ia)
               cb = s_ctrl % ipc(ib)

               call meptr(ca,cb,iam,npm,k)
               if (k == 0) cycle

               call dlmn(s_ctrl % nbas, s_lat % plat, s_lat % pos, iax(0, isite), dcs)
               memode = memodk(k)
               call makvme(memode,d,tblid,ilme,nterm,tabme(1,1,k),decay(1,k), &
                         & alat,lscale, alat*dcs(0),poly(k),cutmod(k),cut(1,1,k),v,dv)
               call skham(dcs(1), dcs(2), dcs(3), v, nl, nl2, ham(1:nl2, 1:nl2, isite, spidx))
               if (memode == 6) call gskham(dcs(1), dcs(2), dcs(3), v, nl, nl2, ham(1:nl2, 1:nl2, isite, spidx))

               if (hdiffs) call dskham(dcs(1),dcs(2),dcs(3),dcs(0)*alat,v,dv,nl,nl2,cryf, &
                                          & dhm(1:nl2, 1:nl2, isite, spidx, 1), &
                                          & dhm(1:nl2, 1:nl2, isite, spidx, 2), &
                                          & dhm(1:nl2, 1:nl2, isite, spidx, 3), &
                                          & dhm(1:nl2, 1:nl2, isite, spidx, 4))

               if (lov) then
                  call makvme(memode,d,tblid,ilme,nterm,tabov(1,1,k),decov(1,k), &
                           & alat,lscale, alat*dcs(0),poly(k),cutmod(k),cutov(1,1,k),v,dv)
                  call skham(dcs(1), dcs(2), dcs(3), v, nl, nl2, ovl(1:nl2, 1:nl2, isite, spidx))
                  if (memode == 6) call gskham(dcs(1), dcs(2), dcs(3), v, nl, nl2, ovl(1:nl2, 1:nl2, isite, spidx))

                  if (hdiffs) call dskham(dcs(1),dcs(2),dcs(3),dcs(0)*alat,v,dv,nl,nl2,ocryf, &
                                          & dov(1:nl2, 1:nl2, isite, spidx, 1), &
                                          & dov(1:nl2, 1:nl2, isite, spidx, 2), &
                                          & dov(1:nl2, 1:nl2, isite, spidx, 3), &
                                          & dov(1:nl2, 1:nl2, isite, spidx, 4))
               end if

            end do

            call tcx('ham etc..')
!
!             call tcn('ham dist')
!
! !             Distribute
!             if (naproc > 1) then
!                call mpi_allgatherv(mpi_in_place, 0, mpi_real8, &
!                                   & ham(1,1,1,spidx), recvcnt, displs, mpi_real8, comm, err)
!                if (hdiffs) then
!                   do crd = 1, 4
!                      call mpi_allgatherv(mpi_in_place, 0, mpi_real8, &
!                                         & dhm(1,1,1,spidx,crd), recvcnt, displs, mpi_real8, comm, err)
!                   end do
!                end if
!
!                if (lov) then
!                   call mpi_allgatherv(mpi_in_place, 0, mpi_real8, &
!                                      & ovl(1,1,1,spidx), recvcnt, displs, mpi_real8, comm, err)
!                   if (hdiffs) then
!                      do crd = 1, 4
!                         call mpi_allgatherv(mpi_in_place, 0, mpi_real8, &
!                                            & dov(1,1,1,spidx,crd), recvcnt, displs, mpi_real8, comm, err)
!                      end do
!                   end if
!                end if
!             end if
!
!             call tcx('ham dist')

            call tcn('ham shuffle')

            if (nl > 1) then
               do isite = 1, nsites
                  ia = iax(0, isite)
                  onsite = isite == npr(1, ia) + 1
                  if (onsite) cycle

                  rsite = iax(5, isite)

!          take sp for ps etc.. later try to use breadcrumbs and reset both at the same time halving the accesses
!                 ps
                  ham(2:4,   1, isite, spidx) =  ham(  1, 2:4, rsite, spidx)

                  if (nl > 2) then
!                    ds
                     ham(5:9,   1, isite, spidx) =  ham(  1, 5:9, rsite, spidx)
!                    dp
                     ham(5:9, 2:4, isite, spidx) =  transpose(ham(2:4, 5:9, rsite, spidx)) ! may also get the nontransposed after the breadcrumbs are implemented
                  end if

                  if (hdiffs) then
                     do crd = 1, 3 ! do not include dr!
                        dhm(2:4,   1, isite, spidx, crd) =  -dhm(  1, 2:4, rsite, spidx, crd)
                        if (nl > 2) then
                           dhm(5:9,   1, isite, spidx, crd) =  -dhm(  1, 5:9, rsite, spidx, crd)
                           dhm(5:9, 2:4, isite, spidx, crd) =  -transpose(dhm( 2:4, 5:9, rsite, spidx, crd))
                        end if
                     end do
                     crd = 4
                     dhm(2:4,   1, isite, spidx, crd) = dhm(  1, 2:4, rsite, spidx, crd)
                     if (nl > 2) then
                        dhm(5:9,   1, isite, spidx, crd) = dhm(  1, 5:9, rsite, spidx, crd)
                        dhm(5:9, 2:4, isite, spidx, crd) = transpose(dhm( 2:4, 5:9, rsite, spidx, crd))
                     end if
                  end if

                  if (lov) then
                     ovl(2:4,   1, isite, spidx) =  ovl(  1, 2:4, rsite, spidx)
                     if (nl > 2) then
                        ovl(5:9,   1, isite, spidx) =  ovl(  1, 5:9, rsite, spidx)
                        ovl(5:9, 2:4, isite, spidx) =  transpose(ovl(2:4, 5:9, rsite, spidx)) ! may also get the nontransposed after the breadcrumbs are implemented
                     end if
                     if (hdiffs) then
                        do crd = 1, 3 ! do not include dr!
                           dov(2:4,   1, isite, spidx, crd) =  -dov(  1, 2:4, rsite, spidx, crd)
                           if (nl > 2) then
                              dov(5:9,   1, isite, spidx, crd) =  -dov(  1, 5:9, rsite, spidx, crd)
                              dov(5:9, 2:4, isite, spidx, crd) =  -transpose(dov( 2:4, 5:9, rsite, spidx, crd))
                           end if
                        end do
                        crd = 4
                        dov(2:4,   1, isite, spidx, crd) = dov(  1, 2:4, rsite, spidx, crd)
                        if (nl > 2) then
                           dov(5:9,   1, isite, spidx, crd) = dov(  1, 5:9, rsite, spidx, crd)
                           dov(5:9, 2:4, isite, spidx, crd) = transpose(dov( 2:4, 5:9, rsite, spidx, crd))
                        end if
                     end if
                  end if
               end do
            end if
            call tcx('ham shuffle')
         end do
      end do

      call tcx('newham')
!
!       call print_ham(nspc,nl,s_ctrl%ipc,nsites,iax,ham,'ham_new')
!       call print_ham(nspc,nl,s_ctrl%ipc,nsites,iax,dhm(:,:,:,:,1),'dhamx_new')
!       call print_ham(nspc,nl,s_ctrl%ipc,nsites,iax,dhm(:,:,:,:,2),'dhamy_new')
!       call print_ham(nspc,nl,s_ctrl%ipc,nsites,iax,dhm(:,:,:,:,3),'dhamz_new')
!       call print_ham(nspc,nl,s_ctrl%ipc,nsites,iax,dhm(:,:,:,:,4),'dhamr_new')
!
!       deallocate(ham)
!       call rx0('newhams PRINTED')

   end subroutine bldham




!    subroutine print_ham(nspc,nl,ipc,nsites,iax,ham,fname)
!       implicit none
!       integer, parameter :: dp = 8
!
!       integer, intent(in) :: nl, nspc, nsites, iax(0:9,nsites), ipc(*)
!       real(dp), intent(inout) :: ham(nl**2, nl**2, nsites, nspc**2)
!       character(len=*), intent(in) :: fname
!
!
!       integer :: isp1, isp2, spidx, ia, ib, ca, cb, isite, rsite, nl2, u
!
!       character(len=100) :: frm
!
!       nl2 = nl*nl
!
!       where( -1.0e-16_dp < ham .and. ham < 1.0e-16_dp) ham = 0.0_dp
!
!       write(frm, '(a,i0,a,i0,a)' ) '(',nl2,'(',nl2,'(x,f16.8),/))'
!
!       open(newunit=u, file=fname, action='write')
!
!       do isp2 = 0, nspc-1
!          do isp1 = 1, nspc
!
!             spidx = isp2*nspc + isp1
!
!             do isite = 1, nsites
!                ia = iax(0, isite)
!                ib = iax(1, isite)
!                ca = ipc(ia)
!                cb = ipc(ib)
!
!                rsite = iax(5,isite)
!
!                write(u, '(2(x,i8),4(x,i5))') isite, rsite, ia, ca, ib, cb
!                write(u, frm) ham(1:nl2, 1:nl2, isite, spidx)
!
!             end do
!          end do
!       end do
!
!       close(u)
!
!
!    end subroutine print_ham

