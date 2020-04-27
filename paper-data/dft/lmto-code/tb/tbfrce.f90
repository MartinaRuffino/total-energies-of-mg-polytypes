      subroutine tbfrce(tbc, mode,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp, &
     &  lmx,ipc,idxdn,indxsh,iwk,ldim,nev,zval,bk,wt,wts,                     &
     &  norder,width,metal,efermi,nsite,iax,npr,xyzfrz,eband,zdist,ldz,typesz,&
     &  h,h0,dh, dhcf,vso,hso,s,ds,dscf,charge,rl,rn,trh,pv,force,lso,lov,    &
     &  cryf,ocryf,sumev,entrpy,f,thrpv,e,rho,rholm,rhoc,rhon,drhos)
!C- TB band forces, energies, charges, and pressure for one k-point
!C-----------------------------------------------------------------------
!Ci Inputs
!Ci   In mode 0 (sampling or insulators only) the Fermi energy or no. of
!Ci   states is input; in mode 1 the k-point weights are input in wts;
!Ci   lmol: molecule (cluster) branch
!Ci   lgamma: gamma point only (hamiltonian real); must be TRUE if lmol
!Ci   plat,nbas,nl,lmx,ipc,indxsh
!Ci   nspc = 2 for coupled spins (empirical S-O)
!Ci   nsp1 = 2 if S-O or TB+U (dimensions eband and BZ integration)
!Ci   nsp  = 2 if TB+U or TB-L spin polarised, isp is then current spin
!Ci   iwk, number of orbitals for each atom
!Ci   ldim, dimension of "lower" block of hamiltonian matrix
!Ci   nev, number of eigenvectors found from diagno
!Ci   zval, total no. of valence electrons
!Ci   bk,wt: k-point, and k-point degeneracy weight
!Ci   wts, band weights for all k-points
!Ci   norder, width, BZ sampling parameters (metals only); metal,efermi;
!Ci   nsite, total number of neighbors in all clusters;
!Ci   iax, neighbor lists; npr, see tbham;
!Ci   eband, bands this k; z, eigenvectors this k;
!Ci   h, h0, dh(x,y,z,r), real space hamiltonian and derivatives;
!Ci   dhcf(x,y,z,r), real space crystal field derivatives
!Ci   vso,hso: table of spin-orbit parameters and the hamiltonian
!Ci   s, ds(x,y,z,r), real space overlap and derivatives;
!Ci   dscf(x,y,z,r), real space overlap crystal field derivatives
!Ci   Switches: charge, return s,p,d charges on each atom;
!Ci             rl, return on-site c*_RL c_RL' in rhoc (see Outputs)
!Ci             rn, return on-site c*_RL S c_RL' in rhon (see Outputs)
!Ci             trh, return band energy on each atom;
!Ci             pv, calculate 3PV;
!Ci             force, calculate forces;
!Ci             lso, include spin-orbit interactions
!Ci             lov, include overlap matrix (non-orthogonal TB)
!Ci   cryf, true if crystal field terms in Hamiltonian
!Ci   ocryf, true if overlap crystal field terms
!Ci   xm0,wk,wk2,wkcf: work arrays
!Co Outputs
!Co   sumev, sum occupied levels if lmol (i.e., a molecule)
!Co   entrpy, entropy term (actually TS)
!Co   f(3,nbas) forces on each atom in the basis accumulated for this k;
!Co   thrpv is 3PV for this k;
!Co   e(nbas,nsp1) band energy for each atom for this k
!Co   rho(nl,2,nbas) s,p,d charges for each atom
!Co                                       accumulated for this k, spin
!Co   rholm(nl**2,2,nbas) {lm} Mulliken charges for each atom
!Co                                       accumulated for this k, spin
!Co                         non spin pol for now
!Co   rhoc(nl**2,nl**2,nbas) s,p,d eigenvector products  c*_RL c_RL'
!Co                                       accumulated for this k and spin
!Co                                       and summed over spin
!Co   rhon(nl**2,nl**2,nbas,nsp=2) occupation numbers for each atom
!Co                                       accumulated for this k and spin
!Co   drhos(ldim,ldim,3) c*_RL dS/dR c_RL' + cc
!Co                                       accumulated for this k
!Cr Remarks
!Cr   Forces, rho, rho[ln] are symmetrized after k-point sum, see symfor
!Cr   symr and symrtb
!Cu Updates
!Cu         2013 (DP) rewrite to use sparsity. Bands loop is innermost.
!Cu   04 Jun 08 (ATP) new molecule (cluster) mode
!C-----------------------------------------------------------------------


      use mpi
      use tbprl
      use mod_ctx

      implicit none
!C Passed parameters
      type(tbc_t), intent(in), target :: tbc
      integer mode,nbas,nl,nspc,nsp,isp,nsp1,ldim,nev,norder,nsite,ldz, typesz
      integer lmx(*),ipc(*),indxsh(*),iwk(nbas),iax(0:9,nsite), npr(0:1,nbas),idxdn(0:nl-1,*)
      real(8) :: zval,wt,width,efermi,thrpv,sumev,entrpy
      real(8), intent(inout), target :: zdist(typesz,ldz,tbc%loclsz(2))
      real(8) :: plat(3,3),bk(3),wts(nev),eband(nev),                  &
     &  h0(nl**4*nsite*nspc**2*nsp),                                     &
     &  h(nl**4*nsite*nspc**2*nsp),                                      &
     &  dh(nl**4*nsite*nspc**2,4),                                       &
     &  dhcf(nl**4*nsite*nspc**2,4),vso(*),hso(*),                       &
     &  s(nl**4*nsite*nspc**2*nsp),                                      &
     &  ds(nl**4*nsite*nspc**2,4),dscf(nl**4*nsite*nspc**2,4),           &
     &  f(3,nbas), e(nbas),rho(0:nl-1,nsp1,nbas),rholm(nl**2,nsp1,nbas), &
     &  drhos(3,nbas,nbas),                                              &
     &  rhon(0:nl**2-1,0:nl**2-1,nbas),                                  &
     &  rhoc(0:nl**2-1,0:nl**2-1,nbas)
      logical lmol,lgamma,metal,charge,rl,rn,trh,pv,force,lso,lov,cryf,ocryf,xyzfrz(3)
      integer nstate,nst0,ntry,itry,n,n2,l,m,i,j,ib,id,ibas,indx,iprint,lsp,i1mach,is,k,ilm
      double precision dr,zv0,e1,e2,dosef,cv
      character*80 outs
      logical cmdopt
      logical, parameter :: Fls = .false.,T = .true.

      integer :: isite,istate,ixa,ixb,nlma,nlmb,la,lb,lma,lmb,indxH(nbas),ia,ipa,nrange,nend
      real(8) :: eloc, floc(3)
      real(8), allocatable :: hk0(:,:,:), dhk(:,:,:,:), zr(:,:,:), sk(:,:,:), dsk(:,:,:,:)
      real(8), pointer :: z(:,:,:), zt(:,:,:)
      logical ::  sqt !, walked(nbas)
!       logical, parameter :: naive =.true.
      logical, parameter :: naive =.false.
      real(8) :: rescale
      integer :: tmp(nbas), pid

!       real(8) :: wk(ldim,ldim,2),wk2(ldim,ldim,2),wkcf(nbas), xm0(ldim*nspc)

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
      real(8) :: rm(typesz,0:nl*nl-1,0:nl*nl-1),  rhos(typesz,0:nl*nl-1,0:nl*nl-1)
      real(8), dimension(0:nl*nl-1,0:nl*nl-1) :: rmc
      real(8) :: rv(0:nl*nl-1), rc(0:nl-1)

      real(8) ::  nul(2) = [0.0_8, 0.0_8], one(2) = [1.0_8, 0.0_8]

      logical :: lroot, offsite
      type(cart_ctx_t), pointer :: c2d


      procedure() :: s_ddot, dgemm, dgemv, s_zdotc, zgemm, zgemv
      procedure(), pointer :: gemm => null(), gemv => null(), s_dot => null()

      character :: mt

      integer :: destrank, lsz, err, mtype, bsz, tag, li, lj, u, lz, nlcblks, reqsz
!       real(8), allocatable ::  ztdist(:,:,:)
      integer, allocatable :: rreqs(:), sreqs(:) !, statuses(:,:)
      integer, dimension(2) :: lbcrds, gbcrds, lcrds, remcrd, lendc, tgbcrd, destcrd
      real(8), allocatable :: sblocks(:,:,:,:,:), rblocks(:,:,:,:,:), focc(:)
      integer, pointer :: nlclus, lclus(:)

      logical :: ptrans = .true.

      call tcn('tbfrce')

      if (lgamma) then
         mt = 't'
         s_dot  => s_ddot
         gemm => dgemm
         gemv => dgemv
      else
         mt = 'c'
         s_dot  => s_zdotc
         gemm => zgemm
         gemv => zgemv
      end if


      lroot  = tbc % c3d % lrt
      pid    = tbc % c2d % id

      c2d => tbc % c2d

      lz = tbc % globsz(1)
      if (naive) allocate(zr(typesz,lz,lz))


!       if (.not. aprox) then
      n = ldim
      n2 = n*n
      lsp = ldim / nspc
      call rxx(lgamma.and.nspc==2, 'TBFRCE not set up for S-O and gamma point only')
      call rxx((.not.trh).and.lov.and.(rn.or.rl),'TBFRCE: for s-c TB with overlap restart with TRH=T in ctrl')
!       call defi(oiwk, nbas)
!       if (lgamma .and. iprint() >= 30) then
!         if ((iprint() > 30 .and. ldim < 10).or. iprint() > 60) then
!           print *, 'TBFRCE: real eigenvectors ..'
!           do  i = 1, ldim
!             write(*,'(1028f8.4)') (z(1,i,j),j=1,ldim)
!           enddo
!         endif
!       endif

!C --- Set k-point weights if insulator or mode 0 ---
      nstate = nev
      if (.not. metal .and. .not. lmol) then
        nstate = (zval + 1d-4) / (3 - nspc)
        wts(1:nstate) = wt
        entrpy = 0d0
        efermi = 0.5_8*(eband(nstate) + eband(nstate+1))
      elseif (mode == 0) then
        if (lmol) then
          ntry = 100
          e1 = eband(1)
          e2 = eband(ldim)
          nst0 = (zval + 1d-4) / (3 - nspc)
          efermi = eband(nst0)
          do  itry = 1, ntry
            call pshprt(0)
            call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,T,sumev,wts(1),zv0,entrpy,dosef,cv)
            call popprt
            if (zv0 > zval) e2 = efermi
            if (zv0 <= zval) e1 = efermi
!C           ... could use qtol here:
            if (dabs(zval - zv0) < 1d-12) goto 1
            efermi = 0.5d0*(e1 + e2)
          enddo
          if (iprint() >10) print *, 'TBFRCE: ***warning*** cannot find HOMO/LUMO'
    1     continue

          if (iprint() >= 30) then
            call awrit1(' TBFRCE: locate molecule ''Fermi'' energy ... %i tries,',' ',128,i1mach(2),itry)
            if (norder >= 0) then
              call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term: %d, %d electrons',' ', &
                                             &   256,i1mach(2),norder,width,efermi,sumev,entrpy,zv0)
            else
              call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%d, %d electrons',' ', &
                                             &   256,i1mach(2),0.1579d6*width,efermi,sumev,entrpy,zv0)
            endif
            if (iprint() > 30) &
               & call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,T,sumev,wts,zv0,entrpy,dosef,cv)
          endif
        else
          call pshprt(0)
          call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,metal,sumev,wts,zv0,entrpy,dosef,cv)
          call popprt
        endif
      endif


      nl2 = nl*nl

      mtype = mpi_real8
      if (.not. lgamma) mtype = mpi_complex16

      sqt = .not. (lov .or. any(wts(1:nstate) < 0.0_8))

!       if ( .not. sqt) then !iprint() > 20 .and.
!          write(*,"(' TBFRCE: Not using the fast sqrt(focc) mode due to ')",advance='no')
!          if (lov) then
!             write(*,"('nonorthogonality.')")
!          elseif (any(wts(1:nstate) < 0.0_8)) then
!             write(*,"('presence of negative occupancies.')")
!             print *, wts
!          end if
!       end if

      if (.not. tbc % sl) ptrans = .false.

      if (tbc%sl .and. tbc%d2count(pid) > 0) then
         allocate(z(typesz,lz,tbc % globsz(2)))
      else
         z => zdist
      end if

      if (sqt) then
         if (tbc%sl .and. ptrans) then
            allocate(focc(tbc % loclsz(2)))
            focc = 0.0_8
            do j = 0, (nstate/tbc%blcksz(2))/c2d%szs(1)-1
               i = j*c2d%szs(1) + c2d%crd(1)
               focc(j*tbc%blcksz(2)+1 : (j+1)*tbc%blcksz(2)) = sqrt(wts(i*tbc%blcksz(2)+1 : (i+1)*tbc%blcksz(2)))
            end do
            k = mod(nstate,tbc%blcksz(2))
            if (k>0) then
               j = (nstate/tbc%blcksz(2))/c2d%szs(1)
               i = j*c2d%szs(1) + c2d%crd(1)
               focc(j*tbc%blcksz(2)+1 : j*tbc%blcksz(2)+k) = sqrt(wts(i*tbc%blcksz(2)+1 : i*tbc%blcksz(2)+k))
            end if
         else
            allocate(focc(lz))
            focc = 0.0_8
            focc(1:nstate) = sqrt(wts)
         end if
      end if
!
!       write(100+c2d%id, '(g0)') focc
!       flush(100+c2d%id)
!       call mpi_barrier(c2d%comm, err)
!       stop 'focc'

      call tcn('Traces')

      if (tbc%sl) then
         if (ptrans) then

            call tcn('pTrans')
            lsz = product(tbc%loclsz)

            if (c2d%szs(0) == c2d%szs(1)) then
               call tcn('trans')
               if (sqt) then
                  call mul_ctran_inplace(typesz, tbc % loclsz(1), zdist, ldz, focc)
               else
                  call     ctran_inplace(typesz, tbc % loclsz(1), zdist, ldz)
               end if
               call tcx('trans')

               call tcn('sqwap')
               destcrd = [c2d%crd(1), c2d%crd(0)]
               call mpi_cart_rank(c2d%comm, destcrd, destrank, err)
               call mpi_sendrecv_replace(zdist, lsz, mtype, destrank, destrank, destrank, c2d%id, c2d%comm,mpi_status_ignore,err)
               call tcx('sqwap')
            else

               call tcn('srambled swap')
               nlcblks = product(tbc%lcblks)
               allocate(rblocks(typesz,tbc%blcksz(1),tbc%blcksz(2),0:tbc%lcblks(1)-1,0:tbc%lcblks(2)-1),rreqs(0:nlcblks-1))
               rreqs = mpi_request_null
               bsz = product(tbc%blcksz)
               do j = 0, tbc%lcblks(2)-1
                  do i = 0, tbc%lcblks(1)-1
                     lbcrds = [i,j]
                     gbcrds = lbcrds*c2d%szs(0:1) + c2d%crd(0:1)
                     tgbcrd = [gbcrds(2), gbcrds(1)]
                     if (any(gbcrds >= tbc%blocks) .or. any(tgbcrd >= tbc%blocks)) cycle
                     remcrd = mod(tgbcrd,c2d%szs(0:1))
                     call mpi_cart_rank(c2d%comm, remcrd, destrank, err)
                     tag = gbcrds(2)*tbc%blocks(1)+gbcrds(1)
                     call mpi_irecv(rblocks(1,1,1,i,j), bsz, mtype, destrank, tag, c2d%comm, rreqs(i+j*tbc%lcblks(1)), err)
!                      write(350+c2d%id,'(a,2(x,i2),3(a,i3))') 'recv bl: ',gbcrds,' from: ',destrank,' to: ',c2d%id,' tag: ',tag
                  end do
               end do

               allocate(sblocks(typesz,tbc%blcksz(1),tbc%blcksz(2),0:tbc%lcblks(1)-1,0:tbc%lcblks(2)-1),sreqs(0:nlcblks-1))
               sreqs = mpi_request_null
               call mpi_barrier(c2d%comm, err)

               do j = 0, tbc%lcblks(2)-1
                  do i = 0, tbc%lcblks(1)-1
                     lbcrds = [i,j]
                     lcrds = lbcrds*tbc%blcksz
                     lendc = lcrds+tbc%blcksz

                     gbcrds = lbcrds*c2d%szs(0:1) + c2d%crd(0:1)
                     tgbcrd = [gbcrds(2), gbcrds(1)]
                     if (any(gbcrds >= tbc%blocks) .or. any(tgbcrd >= tbc%blocks)) cycle

                     if (sqt) then
                        call mul_ctran_ofplace(typesz, tbc % blcksz(1), tbc % blcksz(2), &
                           & zdist(1,lcrds(1)+1,lcrds(2)+1), ldz, sblocks(1,1,1,i,j), tbc%blcksz(1), focc(lcrds(2)+1))
                     else
                        call ctran_ofplace(typesz, tbc % blcksz(1), tbc % blcksz(2), &
                           & zdist(1,lcrds(1)+1,lcrds(2)+1), ldz, sblocks(1,1,1,i,j), tbc%blcksz(1))
                     end if

                     remcrd = mod(tgbcrd,c2d%szs(0:1))
                     call mpi_cart_rank(c2d%comm, remcrd, destrank, err)
                     tag = tgbcrd(2)*tbc%blocks(1)+tgbcrd(1)
                     call mpi_irsend(sblocks(1,1,1,i,j), bsz, mtype, destrank, tag, c2d%comm, sreqs(i+j*tbc%lcblks(1)), err)
!                      write(350+c2d%id,'(a,2(x,i2),3(a,i3))') 'send bl: ',tgbcrd,' from: ',c2d%id,' to: ',destrank,' tag: ',tag
                  end do
               end do


!                flush(350+c2d%id)
!                call mpi_barrier(c2d%comm, err)
!                stop
!                allocate(statuses(mpi_status_size,nlcblks))

               call mpi_waitall(nlcblks, rreqs, mpi_statuses_ignore, err)
               deallocate(rreqs)
               call mpi_waitall(nlcblks, sreqs, mpi_statuses_ignore, err)
               deallocate(sblocks,sreqs)

!                deallocate(statuses)

               do j = 0, tbc%lcblks(2)-1
                  do i = 0, tbc%lcblks(1)-1
                     lbcrds = [i,j]
                     lcrds = lbcrds*tbc%blcksz
                     lendc = lcrds+tbc%blcksz
                     zdist(1:typesz,lcrds(1)+1:lendc(1),lcrds(2)+1:lendc(2)) = rblocks(1:typesz,1:tbc%blcksz(1),1:tbc%blcksz(2),i,j)
                  end do
               end do

               deallocate(rblocks)
               call tcx('srambled swap')
            end if

!       do j = 0, tbc%lcblks(2)-1
!          do i = 0, tbc%lcblks(1)-1
!             lbcrds = [i,j]
!             lcrds = lbcrds*tbc%blcksz
!             lendc = lcrds+tbc%blcksz
!             gbcrds = lbcrds*c2d%szs(0:1) + c2d%crd(0:1)
!             do lj = 1, tbc % blcksz(2)
!                do li = 1, tbc % blcksz(1)
!                   zdist(1:typesz,lcrds(1)+li,lcrds(2)+lj) = ztdist(1:typesz,lcrds(1)+li,lcrds(2)+lj)*wts(gbcrds(1)+li)
!                end do
!             end do
!          end do
!       end do
         call tcx('pTrans')
         end if

         call tcn('gthr')
         call darray_gather(mtype, z, zdist, tbc%desc, c2d)
         call tcx('gthr')
         if (tbc%d2count(c2d%rt) /= nbas) then
             call tcn('bcast')
             call mpi_bcast(z, lz*lz, mtype, c2d%rt, c2d%comm, err)

!          call r8_printu(tbc%realsz(1), tbc%realsz(2), z, lz, 360+c2d%id)
!          flush(360+c2d%id)
!          call mpi_barrier(c2d%comm, err)
!          stop

!             call tcn('a2aw')
!             call darray_gathera(mtype, z, zdist, tbc%desc, tbc % c2d)
!             call tcn('a2aw')


! ++++++++++++++++ in this branch the root sends out only the necessary columns
! ++++++++++++++++    to the processes with nonblocking ready mode, but it can
! ++++++++++++++++    only be faster if the atoms on each process are fairly
! ++++++++++++++++    isolated and even then seems dubious since the communication is per atom.
! ++++++++++++++++    With some buffer it may be better but then the size requirements may go into lalaland
! ++++++++++++++++    if the isolations condition is not really met.
!          if (.not. c2d % lrt) then
!             nlclus => tbc % nlclus(pid)
!             lclus => tbc % lclus(1:nlclus,pid)
!             reqsz = nlclus
!             allocate(rreqs(nlclus,1))
!             do ipa = 1, nlclus
!                ia = lclus(ipa)
!                ixa = indxH(ia)+1
!                tag = pid*nbas+ia
!                call mpi_irecv(zt(1,1,ixa), iwk(ia)*lz, mtype, c2d%rt, tag, c2d%comm, rreqs(ipa,1), err)
!             end do
!          end if
!
!          call mpi_barrier(c2d % comm, err)
!
!          if (c2d % lrt) then
!             allocate(rreqs(sum(tbc%nlclus),1))
!             reqsz = 1
!             do j = 0, c2d % sz-1
!                nlclus => tbc % nlclus(j)
!                lclus => tbc % lclus(1:nlclus,j)
!                do ipa = 1, nlclus
!                   ia = lclus(ipa)
!                   ixa = indxH(ia)+1
!                   tag = j*nbas+ia
!                   call mpi_irsend(zt(1,1,ixa), iwk(ia)*lz, mtype, j, tag, c2d%comm, rreqs(reqsz,1), err)
!                   reqsz = reqsz + 1
!                end do
!             end do
!             reqsz = reqsz - 1
!          end if
!
!          allocate(statuses(mpi_status_size,reqsz))
!          call mpi_waitall(reqsz, rreqs, statuses, err)
!          deallocate(statuses, rreqs)
!+++++++++++++++++++++++++++++++++++++++++++++++++ end

             call tcx('bcast')
   !       deallocate(zdist)

!
!          call tcn('a2aw')
!          call darray_gathera(mtype, zt,  ztdist, tbc%desc, tbc % c2d)
!          call tcx('a2aw')
!          deallocate(ztdist)
         end if

         if (.not. ptrans) then
            if (sqt) then
               call mul_ctran_inplace(typesz, lz, z, lz, focc)
            else
               call     ctran_inplace(typesz, lz, z, lz)
            end if
         end if

      else
         if (sqt) then
            call mul_ctran_inplace(typesz, lz, z, lz, focc)
         else
            call     ctran_inplace(typesz, lz, z, lz)
         endif

!          call r8_printu(tbc%realsz(1), tbc%realsz(2), z, lz, 360+c2d%id)
!          flush(360+c2d%id)
!          call mpi_barrier(c2d%comm, err)
!          stop
      end if


      if (tbc % d2count(pid) > 0) then
!       call dcopy(typesz*lz*lz, zt, 1, z, 1)
! This will cause problems if the passed Z is smaller than (lz,lz).
! A potential solution is to receive it with dimension (*) and dynamically reshape by calculating the appropriate indices
          if (sqt) then
             zt => z
          else
             allocate(zt(typesz,lz,lz))
             call tcn('occmult')
             do j = 1, ldim
                do i = 1, nstate
                   zt(1:typesz, i, j) = z(1:typesz, i, j) * wts(i)
                end do
             end do
             call tcx('occmult')
          end if

!       if (naive .and. .not. aprox) call dgemm('t','n',ldim, ldim, nstate, 1.0_8, zt, ldim, zt, ldim, 0.0_8, zr, ldim)
!       if (naive) call gemm(mt,'n',ldim, ldim, nstate, one, z, lz, zt, lz, nul, zr, lz)



          indxH(1) = 0
          do  ibas = 2, nbas
             indxH(ibas) = indxH(ibas-1) + iwk(ibas-1)
          enddo


          if (trh) allocate(hk0(typesz,lz,lz))

          if (force .and. pv) then
             allocate(dhk(typesz,lz,lz,4))
          else if (force .and. (.not. pv)) then
             allocate(dhk(typesz,lz,lz,3))
          else if ((.not. force) .and. pv) then
             allocate(dhk(typesz,lz,lz,4:4))
          end if


          if (lov) then
             allocate(sk(typesz,lz,lz))
             if (force .and. pv) then
                allocate(dsk(typesz,lz,lz,4))
             else if (force .and. (.not. pv)) then
                allocate(dsk(typesz,lz,lz,3))
             else if ((.not. force) .and. pv) then
                allocate(dsk(typesz,lz,lz,4:4))
             end if
          end if


          if (trh .or. force .or. pv .or. lov) then
             call tcn('blt')
             call pshprt(0)

             if (trh) call ztbloch(lgamma,bk,nl,nspc,nsp,isp,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,h0,vso,hso,Fls, &
                                                                                              & ldim,hk0,lz,typesz)
    ! --- spin-independent: ispu = 1
             if (force) then
                do i = 1, 3
                   call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,dh(1,i),vso,hso,Fls, &
                                                                                        & ldim,dhk(1,1,1,i),lz,typesz)
                end do
             end if
             if (pv) call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,dh(1,4),vso,hso,Fls, &
                                                                                        & ldim,dhk(1,1,1,4),lz,typesz)

             if (lov) then
                call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,s,vso,hso,Fls, &
                                                                                                    & ldim,sk,lz,typesz)
          ! --- spin-independent: ispu = 1
                if (force) then
                   do i = 1, 3
                      call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,ds(1,i),vso,hso,Fls, &
                                                                                           & ldim,dsk(1,1,1,i),lz,typesz)
                   end do
                end if
                if (pv) call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,ds(1,4),vso,hso,Fls, &
                                                                                        & ldim,dsk(1,1,1,4),lz,typesz)
             end if

             call popprt()
             call tcx('blt')
          end if


          offsite = trh .or. force .or. pv .or. (lov .and. (rn .or. charge))

          call tcn('walk')

          do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)

!          nrange = 1
!          if (offsite) nrange = npr(0,ia)
             nend = tbc % neighmap(ia)+1
             if (offsite) nend = tbc % neighmap(ia+1)

!          walked = .false.

             if (trh)   eloc = 0.0_8
             if (force) floc = 0.0_8
             if (lov .and. (rn .or. charge)) rhos = 0.0_8

             do ipa = tbc % neighmap(ia)+1, nend
                ib = tbc % neighidx(ipa)

!          do ipa = 1, nrange
!             ib = iax(1, npr(1, ia)+ipa)
!             if (walked(ib)) cycle

                ixa = indxH(ia)+1
                ixb = indxH(ib)+1

                nlma = iwk(ia)
                nlmb = iwk(ib)

                if (naive) then
                   rm(1:typesz, 0:nlma-1,0:nlmb-1) = zr(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1)
                else
                   if ( nlma == 1 .and. nlmb == 1) then
!                     call pddot(nstate, rm(1,0,0), z(1,1, ixa), IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
                      call s_dot(nstate, z(1,1, ixa), 1, zt(1,1, ixb), 1, rm(1, 0, 0))
!                   print *, 'rm(0,0)', rm(0,0), nstate, zdotu(nstate, zt(1, ixa), 1, zt(1, ixb), 1)
!                   call r8_print(2,nstate,zt(1,ixa),'zt-ixa')
!                   call r8_print(2,nstate,zt(1,ixa),'zt-ixb')
!                   stop
                   else if (nlma == 1) then
                      call gemv(mt, nstate, nlmb, one, z(1,1,ixb), lz, zt(1,1,ixa), 1, nul, rm, 1)
                      rm(1,0, 0:nlmb-1) = rm(1,0:nlmb-1, 0)
                      if (.not. lgamma) rm(2, 0, 0:nlmb-1) = -rm(2,0:nlmb-1, 0)
                   else if (nlmb == 1) then
                      call gemv(mt, nstate, nlma, one, z(1,1,ixa), lz, zt(1,1,ixb), 1, nul, rm, 1)
                   else
                      call gemm(mt,'n', nlma, nlmb, nstate, one, z(1,1,ixa), lz, zt(1,1,ixb), lz, nul, rm, nl2)
                   end if
                end if

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

                   if (charge .and. (.not. lov)) then
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


!             if (lov .and. charge .and. (.not. rn)) rhos(0:nlma-1,0) = rhos(0:nlma-1,0) &
!                   & + sum(sum(sk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1) * rm(1:typesz,0:nlma-1,0:nlmb-1), dim=1), dim=2)
                if (lov .and. (rn .or. charge)) &
                   & call gemm('n','n',nlma, nlma, nlmb, one, rm, nl2, sk(1,ixb,ixa), lz, one, rhos, nl2)


                if (trh) eloc = eloc + sum(hk0(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1) * rm(1:typesz,0:nlma-1,0:nlmb-1))

                if (force) then
                   do i = 1, 3
                      floc(i) = floc(i) + sum(dhk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,i) * rm(1:typesz,0:nlma-1,0:nlmb-1))
                   end do

                   if (lov .and. (rl .or. rn)) then
                      do i = 1, 3
                         drhos(i,ib,ia) = drhos(i,ib,ia) &
                                        & + sum(dsk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,i) * rm(1:typesz,0:nlma-1,0:nlmb-1))
                      end do
                   end if
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

             if (lov .and. (rn .or. charge)) then
                if (rn) then
                   rmc = 0.0_8
                   do lb = 0, nlma - 1
                      lmb = lix(lb, nlma-1)
                      do la = 0, nlma - 1
                         rmc(lix(la, nlma-1), lmb) = rhos(1,la,lb)
                      end do
                   end do
                   rhon(0:nl2-1, 0:nl2-1, ia) = rhon(0:nl2-1, 0:nl2-1, ia) + rmc
                end if
                if (charge) then
                   rv = 0.0_8
                   rc = 0.0_8
                   do la = 0, nlma - 1
                      lma = lix(la, nlma-1)
                      rv(lma) = rhos(1,la,la)
                      rc(llx(lma)) = rc(llx(lma)) + rhos(1,la,la)
                   end do
                   rho(0:nl-1, isp, ia) = rho(0:nl-1, isp, ia) + rc
                   rholm(1:nl2, isp, ia) = rholm(1:nl2, isp, ia) + rv
                end if
             end if
          end do

          if (lov .and. (force .or. pv)) then
             call tcn('ZL')
             do j = 1, ldim
                do istate = 1, nstate
                   z(1:typesz,istate,j) = z(1:typesz,istate,j)*eband(istate)
                end do
             end do
             call tcx('ZL')

             call tcn('Swalk')

             if (naive) call gemm(mt,'n',ldim, ldim, nstate, one, z, lz, zt, lz, nul, zr, lz)

             do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)

!             nrange = 1
!             if (offsite) nrange = npr(0,ia)
                nend = tbc % neighmap(ia)+1
                if (offsite) nend = tbc % neighmap(ia+1)

!             walked = .false.

                floc = 0.0_8

!             do ipa = 1, nrange
                do ipa = tbc % neighmap(ia)+1, nend
!                ib = iax(1, npr(1, ia)+ipa)
                   ib = tbc % neighidx(ipa)

!                if (walked(ib)) cycle

                   ixa = indxH(ia)+1
                   ixb = indxH(ib)+1

                   nlma = iwk(ia)
                   nlmb = iwk(ib)

                   if (naive) then
                      rm(1:typesz, 0:nlma-1,0:nlmb-1) = zr(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1)
                   else
                      if ( nlma == 1 .and. nlmb == 1) then
                         call s_dot(nstate, z(1,1, ixa), 1, zt(1,1, ixb), 1, rm(1, 0, 0))
   !                   print *, 'rm(0,0)', rm(0,0), nstate, zdotu(nstate, zt(1, ixa), 1, zt(1, ixb), 1)
   !                   call r8_print(2,nstate,zt(1,ixa),'zt-ixa')
   !                   call r8_print(2,nstate,zt(1,ixa),'zt-ixb')
   !                   stop
                      else if (nlma == 1) then
                         call gemv(mt, nstate, nlmb, one, z(1,1,ixb), lz, zt(1,1,ixa), 1, nul, rm, 1)
                         rm(1,0, 0:nlmb-1) = rm(1,0:nlmb-1, 0)
                         if (.not. lgamma) rm(2, 0, 0:nlmb-1) = -rm(2,0:nlmb-1, 0)
                      else if (nlmb == 1) then
                         call gemv(mt, nstate, nlma, one, z(1,1,ixa), lz, zt(1,1,ixb), 1, nul, rm, 1)
                      else
                         call gemm(mt,'n', nlma, nlmb, nstate, one, z(1,1,ixa), lz, zt(1,1,ixb), lz, nul, rm, nl2)
                      end if
                   end if

                   if (force) then
                      do i = 1, 3
                         floc(i) = floc(i) + sum(dsk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,i) * rm(1:typesz,0:nlma-1,0:nlmb-1))
                      end do
                   end if

                   if (pv) thrpv = thrpv + sum(dsk(1:typesz,ixa:ixa+nlma-1,ixb:ixb+nlmb-1,4) * rm(1:typesz,0:nlma-1,0:nlmb-1))

!                walked(ib) = .true.
                end do

                if (force) f(1:3,ia) = f(1:3,ia) + 2*floc
             end do

             call tcx('Swalk')
          end if



          call tcx('walk')

!       stop

          if (force .or. pv) deallocate(dhk)
          if ((force .or. pv) .and. lov) deallocate(dsk)
          if (lov) deallocate(sk)
          if (trh) deallocate(hk0)
          if (tbc%sl) deallocate(z)
          if (.not. sqt) deallocate(zt)
      end if

      if (naive) deallocate(zr)

      call tcx('Traces')
      call tcx('tbfrce')


!       write(426,'(2(9(9(x,f16.8),/),/))') rhoc

      end subroutine tbfrce




      subroutine s_ddot(n, dx, incx, dy, incy, r)
         implicit none
         integer, intent(in) :: n, incx, incy
         real(8), intent(in) :: dx(*), dy(*)
         real(8), intent(out) :: r
         procedure(real(8)) :: ddot
         r = ddot(n, dx, incx, dy, incy)
      end subroutine s_ddot

      subroutine s_zdotc(n, zx, incx, zy, incy, r)
         implicit none
         integer, intent(in) :: n, incx, incy
         complex(8), intent(in) :: zx(*), zy(*)
         complex(8), intent(out) :: r
         procedure(complex(8)) :: zdotc
         r = zdotc(n, zx, incx, zy, incy)
      end subroutine s_zdotc





!       subroutine s_pddot(N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )


!
!



      subroutine ctran_ofplace(typesz, n, m, a, lda, b, ldb)
      implicit none
      integer, parameter :: dp = 8

      integer, intent(in) :: n, m, typesz, lda, ldb
      real(8), intent(in)  :: a(typesz,lda,*)
      real(8), intent(out) :: b(typesz,ldb,*)

      integer :: i,j
      real(8) :: t,ct(2)

      if (typesz == 1) then
         do j = 1, m
            do i = 1, n
               b(1,j,i) = a(1,i,j)
            end do
         end do
      else
         do j = 1, m
            do i = 1, n
               b(1:2,j,i) = [a(1,i,j), -a(2,i,j)]
            end do
         end do
      end if

      end subroutine ctran_ofplace


      subroutine mul_ctran_ofplace(typesz, n, m, a, lda, b, ldb, d)
      implicit none
      integer, parameter :: dp = 8

      integer, intent(in) :: n, m, typesz, lda, ldb
      real(8), intent(in)  :: a(typesz,lda,*), d(*)
      real(8), intent(out) :: b(typesz,ldb,*)

      integer :: i,j
      real(8) :: t,ct(2)

      if (typesz == 1) then
         do j = 1, m
            do i = 1, n
               b(1,j,i) = a(1,i,j)*d(j)
            end do
         end do
      else
         do j = 1, m
            do i = 1, n
               b(1:2,j,i) = [a(1,i,j), -a(2,i,j)]*d(j)
            end do
         end do
      end if

      end subroutine mul_ctran_ofplace







      subroutine mul_ctran_inplace(typesz, n, a, lda, d)
!       Conjugate transpose in place.
      implicit none
      integer, parameter :: dp = 8

      integer, intent(in) :: n, typesz, lda
      real(8), intent(inout) :: a(typesz,lda,*), d(*)

      integer :: i,j
      real(8) :: t,ct(2)

      if (typesz == 1) then
         do j = 1, n
            do i = j, n
               t        = a(1,i,j)*d(j)
               a(1,i,j) = a(1,j,i)*d(i)
               a(1,j,i) = t
            end do
         end do
      else
         do j = 1, n
            do i = j, n
               ct         = [a(1,i,j), -a(2,i,j)]*d(j)
               a(1:2,i,j) = [a(1,j,i), -a(2,j,i)]*d(i)
               a(1:2,j,i) = ct
            end do
         end do
      end if

      end subroutine mul_ctran_inplace



      subroutine ctran_inplace(typesz, n, a, lda)
!       Conjugate transpose in place.
      implicit none
      integer, parameter :: dp = 8

      integer, intent(in) :: n, typesz, lda
      real(8), intent(inout) :: a(typesz,lda,*)

      integer :: i,j
      real(8) :: t,ct(2)

      if (typesz == 1) then
         do j = 1, n-1
            do i = j+1, n
               t        = a(1,i,j)
               a(1,i,j) = a(1,j,i)
               a(1,j,i) = t
            end do
         end do
      else
         do j = 1, n-1
            a(2,j,j) = -a(2,j,j) ! conjugate the (pseudo)diagonal
            do i = j+1, n
               ct         = [a(1,i,j), -a(2,i,j)]
               a(1:2,i,j) = [a(1,j,i), -a(2,j,i)]
               a(1:2,j,i) = ct
            end do
         end do
         a(2,n,n) = -a(2,n,n) ! and the (pseudo)last piece of it
      end if

      end subroutine ctran_inplace

      subroutine r8_printu(n,m,a,lda,u)
      implicit none
      integer, intent(in) :: n,m,u,lda
      real(8), intent(in) :: a(lda,*)
      real(8) :: r

      integer :: v, i,j

      write(u, '("#",2(x,i0))') m, n
      do i = 1, m
         do j = 1, n
            r = a(j,i)
            if (abs(r) < 1.0e-8_8) r = abs(r)
            write(u, '(x,es15.8)' ,advance='no') r
         end do
         write(u,'("")')
      end do

      end subroutine r8_printu


!       subroutine r8_print(n,m,a,fln)
!       implicit none
!       integer, intent(in) :: n,m
!       real(8), intent(in) :: a(n,m)
!       character(len=*) :: fln
!
!       integer :: u, v, i,j
!
!       open(newunit=u, file=fln, action='write')
!       write(u, '(2(x,i0))') m, n
!       do i = 1, m
!          do j = 1, n
!             write(u, '(x,es20.12)' ,advance='no') a(j,i)
!          end do
!          write(u,'("")')
!       end do
!       close(u)
!
!       end subroutine r8_print
!
!
!
!       subroutine r8_ppm(n,m,a,fln)
!       implicit none
! !       a(i,j) in [-1:1]
!       integer, intent(in) :: n,m
!       real(8), intent(in) :: a(n,m)
!       character(len=*) :: fln
!
!       integer :: u, v, i,j
!
!
!       open(newunit=u, file=fln, action='write')
!       write(u, '(a,/,a,/,i5,x,i5,/,i3)') 'P3', "# transposed grayscale matrix representation ", n, m, 256
!       do i = 1, m
!          do j = 1, n
!             v = int(256*a(j,i))
!             if (v < 0) then
!                write(u, '(3(x,i3))' ,advance='no') 255+v,255+v,255
!             else if (v > 0) then
!                write(u, '(3(x,i3))' ,advance='no') 255,255-v,255-v
!             else
!                write(u, '(3(x,i3))' ,advance='no') 255,  255, 255
!             endif
!          end do
!          write(u,'("")')
!       end do
!       close(u)
!
!
!
!       end subroutine r8_ppm

!
!       subroutine subtest(wts,eband,z,ldim,lz,typesz,nstate)
!       implicit none
!       integer :: ldim, lz, typesz, nstate, istate
!       real(8) :: z(typesz,lz,lz), eband(ldim), wts(ldim)
!       real(8) ::  nul(2) = [0.0_8, 0.0_8], one(2) = [1.0_8, 0.0_8]
!
!       real(8) :: zt(typesz,lz,lz)
!
!
!       call tcn('msq')
!       if (all(wts(1:nstate) == wts(1))) then
!          z(1:typesz,1:ldim,1:nstate) = z(1:typesz,1:ldim,1:nstate)*sqrt(wts(1))
!       else
!          do istate = 1, nstate
!             z(1:typesz,1:ldim,istate) = z(1:typesz,1:ldim,istate)*sqrt(wts(istate))
!          end do
!       end if
!       call tcx('msq')
!
!       call zgemm('c','n',ldim, ldim, nstate, one, z, lz, z, lz, nul, zt, lz)
!
!       end subroutine subtest
!






!
!       subroutine tbfrce_naive(tbc, mode,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp, &
!      &  lmx,ipc,idxdn,indxsh,iwk,ldim,nev,zval,bk,wt,wts,                     &
!      &  norder,width,metal,efermi,nsite,iax,npr,xyzfrz,eband,z,lz,typesz,    &
!      & h,h0,dh, dhcf,vso,hso,s,ds,dscf,charge,rl,rn,trh,pv,force,lso,lov,     &
!      &  cryf,ocryf,sumev,entrpy,f,thrpv,e,rho,rholm,rhoc,rhon,drhos)
! !C- TB band forces, energies, charges, and pressure for one k-point
! !C-----------------------------------------------------------------------
! !Ci Inputs
! !Ci   In mode 0 (sampling or insulators only) the Fermi energy or no. of
! !Ci   states is input; in mode 1 the k-point weights are input in wts;
! !Ci   lmol: molecule (cluster) branch
! !Ci   lgamma: gamma point only (hamiltonian real); must be TRUE if lmol
! !Ci   plat,nbas,nl,lmx,ipc,indxsh
! !Ci   nspc = 2 for coupled spins (empirical S-O)
! !Ci   nsp1 = 2 if S-O or TB+U (dimensions eband and BZ integration)
! !Ci   nsp  = 2 if TB+U or TB-L spin polarised, isp is then current spin
! !Ci   iwk, number of orbitals for each atom
! !Ci   ldim, dimension of "lower" block of hamiltonian matrix
! !Ci   nev, number of eigenvectors found from diagno
! !Ci   zval, total no. of valence electrons
! !Ci   bk,wt: k-point, and k-point degeneracy weight
! !Ci   wts, band weights for all k-points
! !Ci   norder, width, BZ sampling parameters (metals only); metal,efermi;
! !Ci   nsite, total number of neighbors in all clusters;
! !Ci   iax, neighbor lists; npr, see tbham;
! !Ci   eband, bands this k; z, eigenvectors this k;
! !Ci   h, h0, dh(x,y,z,r), real space hamiltonian and derivatives;
! !Ci   dhcf(x,y,z,r), real space crystal field derivatives
! !Ci   vso,hso: table of spin-orbit parameters and the hamiltonian
! !Ci   s, ds(x,y,z,r), real space overlap and derivatives;
! !Ci   dscf(x,y,z,r), real space overlap crystal field derivatives
! !Ci   Switches: charge, return s,p,d charges on each atom;
! !Ci             rl, return on-site c*_RL c_RL' in rhoc (see Outputs)
! !Ci             rn, return on-site c*_RL S c_RL' in rhon (see Outputs)
! !Ci             trh, return band energy on each atom;
! !Ci             pv, calculate 3PV;
! !Ci             force, calculate forces;
! !Ci             lso, include spin-orbit interactions
! !Ci             lov, include overlap matrix (non-orthogonal TB)
! !Ci   cryf, true if crystal field terms in Hamiltonian
! !Ci   ocryf, true if overlap crystal field terms
! !Ci   xm0,wk,wk2,wkcf: work arrays
! !Co Outputs
! !Co   sumev, sum occupied levels if lmol (i.e., a molecule)
! !Co   entrpy, entropy term (actually TS)
! !Co   f(3,nbas) forces on each atom in the basis accumulated for this k;
! !Co   thrpv is 3PV for this k;
! !Co   e(nbas,nsp1) band energy for each atom for this k
! !Co   rho(nl,2,nbas) s,p,d charges for each atom
! !Co                                       accumulated for this k, spin
! !Co   rholm(nl**2,2,nbas) {lm} Mulliken charges for each atom
! !Co                                       accumulated for this k, spin
! !Co                         non spin pol for now
! !Co   rhoc(nl**2,nl**2,nbas) s,p,d eigenvector products  c*_RL c_RL'
! !Co                                       accumulated for this k and spin
! !Co                                       and summed over spin
! !Co   rhon(nl**2,nl**2,nbas,nsp=2) occupation numbers for each atom
! !Co                                       accumulated for this k and spin
! !Co   drhos(ldim,ldim,3) c*_RL dS/dR c_RL' + cc
! !Co                                       accumulated for this k
! !Cr Remarks
! !Cr   Forces, rho, rho[ln] are symmetrized after k-point sum, see symfor
! !Cr   symr and symrtb
! !Cu Updates
! !Cu   04 Jun 08 (ATP) new molecule (cluster) mode
! !C-----------------------------------------------------------------------
! !       use mpi
!       use tbprl
!
!       implicit none
! !C Passed parameters
!       type(tbc_t), intent(in) :: tbc
!       integer mode,nbas,nl,nspc,nsp,isp,nsp1,ldim,nev,norder,nsite,lz, typesz
!       integer lmx(*),ipc(*),indxsh(*),iwk(nbas),iax(0:9,nsite), npr(0:1,nbas),idxdn(0:nl-1,*)
!       real(8) :: zval,wt,width,efermi,thrpv(*),sumev,entrpy
!       real(8), intent(inout), target :: z(typesz,lz,lz)
!       real(8) :: plat(3,3),bk(3),wts(nev),eband(nev),                  &
!      &  h0(nl**4*nsite*nspc**2*nsp),                                     &
!      &  h(nl**4*nsite*nspc**2*nsp),                                      &
!      &  dh(nl**4*nsite*nspc**2,4),                                       &
!      &  dhcf(nl**4*nsite*nspc**2,4),vso(*),hso(*),                       &
!      &  s(nl**4*nsite*nspc**2*nsp),                                      &
!      &  ds(nl**4*nsite*nspc**2,4),dscf(nl**4*nsite*nspc**2,4),           &
!      &  f(3,nbas), e(nbas),rho(0:nl-1,2,nbas),rholm(nl**2,2,nbas),       &
!      &  drhos(3,nbas,nbas),                                              &
!      &  rhon(0:nl**2-1,0:nl**2-1,nbas,2),                                &
!      &  rhoc(0:nl**2-1,0:nl**2-1,nbas)
!       logical lmol,lgamma,metal,charge,rl,rn,trh,pv,force,lso,lov,cryf,ocryf,xyzfrz(3)
! !C Local variables
! !       integer oiwk
!       integer nstate,nst0,ntry,itry,n,n2,l,m,i,j,ib,id,ibas,indx,iprint,lsp,i1mach,is,k,ilm
!       double precision dr,zv0,e1,e2,dosef,cv
!       character*80 outs
!       logical cmdopt
!       logical, parameter :: Fls = .false.,T = .true.
!
!       integer :: isite,istate,ixa,ixb,nlma,nlmb,la,lb,lma,lmb,indxH(nbas),ia,ipa,nrange
!       real(8) :: eloc, floc(3)
!       real(8), allocatable :: zt(:,:,:), hk0(:,:,:), dhk(:,:,:,:), zr(:,:,:), sk(:,:,:), dsk(:,:,:,:)
!       logical :: walked(nbas)
! !       logical, parameter :: naive =.true.
!       logical, parameter :: naive =.false.
!       real(8) :: rescale
!       integer :: tmp(nbas)
!       integer :: pid
!
!
! !       real(8) :: wk(ldim,ldim,2),wk2(ldim,ldim,2),wkcf(nbas), xm0(ldim*nspc)
!
!
!
!       integer, parameter :: lix(0:8,0:8) = reshape([ &
!      &     0,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 1
!      &    -1,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 2
!      &     1, 2, 3,-1,-1,-1,-1,-1,-1,            & ! 3
!      &     0, 1, 2, 3,-1,-1,-1,-1,-1,            & ! 4
!      &     4, 5, 6, 7, 8,-1,-1,-1,-1,            & ! 5
!      &     0, 4, 5, 6, 7, 8,-1,-1,-1,            & ! 6
!      &    -1,-1,-1,-1,-1,-1,-1,-1,-1,            & ! 7
!      &     1, 2, 3, 4, 5, 6, 7, 8,-1,            & ! 8
!      &     0, 1, 2, 3, 4, 5, 6, 7, 8 ], [9,9])     ! 9
!
!       integer, parameter :: llx(0:8) = [((l,m = -l,l),l=0,2)] ! [0, 1, 1, 1, 2, 2, 2, 2, 2]
!
!       integer :: nl2
!       real(8) :: rm(typesz,0:nl*nl-1,0:nl*nl-1),  rhos(typesz,0:nl*nl-1,0:nl*nl-1)
!       real(8), dimension(0:nl*nl-1,0:nl*nl-1) :: rmc
!       real(8) :: rv(0:nl*nl-1), rc(0:nl-1)
!
!       real(8) ::  nul(2) = [0.0_8, 0.0_8], one(2) = [1.0_8, 0.0_8]
!
!       logical :: lroot, offsite
!
!
!       procedure() :: s_ddot, dgemm, dgemv, s_zdotc, zgemm, zgemv
!       procedure(), pointer :: gemm => null(), gemv => null(), s_dot => null()
!
!       character :: mt
!
! !       print *, 'im tbfrce'
!
!       if (lgamma) then
!          mt = 't'
!          s_dot  => s_ddot
!          gemm => dgemm
!          gemv => dgemv
!       else
!          mt = 'c'
!          s_dot  => s_zdotc
!          gemm => zgemm
!          gemv => zgemv
!       end if
!
!
!       lroot  = tbc % c3d % lrt
!       pid    = tbc % c2d % id
!
!       allocate(zt(typesz,lz,lz))
!       allocate(zr(typesz,lz,lz))
!
!       call tcn('tbfrce')
!
!
! !       if (.not. aprox) then
!       n = ldim
!       n2 = n*n
!       lsp = ldim / nspc
!       call rxx(lgamma.and.nspc==2, 'TBFRCE not set up for S-O and gamma point only')
!       call rxx((.not.trh).and.lov.and.(rn.or.rl),'TBFRCE: for s-c TB with overlap restart with TRH=T in ctrl')
! !       call defi(oiwk, nbas)
! !       if (lgamma .and. iprint() >= 30) then
! !         if ((iprint() > 30 .and. ldim < 10).or. iprint() > 60) then
! !           print *, 'TBFRCE: real eigenvectors ..'
! !           do  i = 1, ldim
! !             write(*,'(1028f8.4)') (z(1,i,j),j=1,ldim)
! !           enddo
! !         endif
! !       endif
!
! !C --- Set k-point weights if insulator or mode 0 ---
!       nstate = nev
!       if (.not. metal .and. .not. lmol) then
!         nstate = (zval + 1d-4) / (3 - nspc)
!         call dcopy(nstate,wt,0,wts,1)
!         entrpy = 0d0
!         efermi = 0.5_8*(eband(nstate) + eband(nstate+1))
!       elseif (mode == 0) then
!         if (lmol) then
!           ntry = 100
!           e1 = eband(1)
!           e2 = eband(ldim)
!           nst0 = (zval + 1d-4) / (3 - nspc)
!           efermi = eband(nst0)
!           do  itry = 1, ntry
!             call pshprt(0)
!             call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,T,sumev,wts(1),zv0,entrpy,dosef,cv)
!             call popprt
!             if (zv0 > zval) e2 = efermi
!             if (zv0 <= zval) e1 = efermi
! !C           ... could use qtol here:
!             if (dabs(zval - zv0) < 1d-12) goto 1
!             efermi = 0.5d0*(e1 + e2)
!           enddo
!           if (iprint() >10) print *, 'TBFRCE: ***warning*** cannot find HOMO/LUMO'
!     1     continue
!
!           if (iprint() >= 30) then
!             call awrit1(' TBFRCE: locate molecule ''Fermi'' energy ... %i tries,',' ',128,i1mach(2),itry)
!             if (norder >= 0) then
!               call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term: %d, %d electrons',' ', &
!                                              &   256,i1mach(2),norder,width,efermi,sumev,entrpy,zv0)
!             else
!               call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%d, %d electrons',' ', &
!                                              &   256,i1mach(2),0.1579d6*width,efermi,sumev,entrpy,zv0)
!             endif
!             if (iprint() > 30) &
!                & call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,T,sumev,wts,zv0,entrpy,dosef,cv)
!           endif
!         else
!           call pshprt(0)
!           call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,metal,sumev,wts,zv0,entrpy,dosef,cv)
!           call popprt
!         endif
!       endif
! !
! !       print *, 'nev:', nev, nstate
! !       print *, 'evals, wts:'
! !       do i = 1, nstate
! !          print *,  eband(i), wts(i)
! !       end do
!
! !
!       nl2 = nl*nl
!
! !       allocate(zt(ldim,ldim))
!
!
!       indxH(1) = 0
!       do  ibas = 2, nbas
!          indxH(ibas) = indxH(ibas-1) + iwk(ibas-1)
!       enddo
!
!
!       if (trh) allocate(hk0(typesz,lz,lz))
!
!       if (force .and. pv) then
!          allocate(dhk(typesz,lz,lz,4))
!       else if (force .and. (.not. pv)) then
!          allocate(dhk(typesz,lz,lz,3))
!       else if ((.not. force) .and. pv) then
!          allocate(dhk(typesz,lz,lz,4:4))
!       end if
!
!
!       if (lov) then
!          allocate(sk(typesz,lz,lz))
!          if (force .and. pv) then
!             allocate(dsk(typesz,lz,lz,4))
!          else if (force .and. (.not. pv)) then
!             allocate(dsk(typesz,lz,lz,3))
!          else if ((.not. force) .and. pv) then
!             allocate(dsk(typesz,lz,lz,4:4))
!          end if
!       end if
!
!
!       call tcn('blt')
!       call pshprt(0)
!
!       if (trh) call ztbloch(lgamma,bk,nl,nspc,nsp,isp,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,h0,vso,hso,Fls, &
!                                                                                           & ldim,hk0,lz,typesz)
!
!       if (lov) then
!          call ztbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,nsite,iax,npr,s,vso,hso,Fls, &
!                                                                                              & ldim,sk,lz,typesz)
!       end if
!
!
!       call popprt
!       call tcx('blt')
!
!       do i = 1, nstate
!          do j = 1, ldim
!             zt(1:typesz, j, i) = z(1:typesz, j, i) * wts(i)
!          end do
!       end do
!
!
!       call gemm('n',mt, ldim, ldim, nstate, one, z, lz, zt, lz, nul, zr, lz)
!
!       offsite = trh .or. force .or. pv .or. lov .and. (rn .or. charge)
!
!
!       if (rl) then
!          do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
!             ixa = indxH(ia)+1
!             nlma = iwk(ia)
!             rhoc(lix(0:nlma-1, nlma-1), lix(0:nlma-1, nlma-1), ia) = zr(1,ixa:ixa+nlma-1,ixa:ixa+nlma-1)
!          end do
!       end if
!
!       if (lov) then
!          call gemm('n','n', ldim, ldim, ldim, one, zr, lz, sk, lz, nul, zt, lz)
!
!          if (rn) then
!             do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
!                ixa = indxH(ia)+1
!                nlma = iwk(ia)
!                rhon(lix(0:nlma-1, nlma-1), lix(0:nlma-1, nlma-1), ia, isp) = zt(1,ixa:ixa+nlma-1,ixa:ixa+nlma-1)
!             end do
!          end if
!
!          if (charge) then
!             do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
!                ixa = indxH(ia)+1
!                nlma = iwk(ia)
!                do la = 0, nlma-1
!                   lma = lix(la,nlma-1)
!                   rholm(lma,isp,ia) = rholm(lma,isp,ia) + zt(1,ixa+la,ixa+la)
!                   rho(llx(lma),isp,ia) = rho(llx(lma),isp,ia) + zt(1,ixa+la,ixa+la)
!                end do
!             end do
!          end if
!       end if
!
!       if (trh) then
!          call gemm('n','n', ldim, ldim, ldim, one, zr, lz, hk0, lz, nul, zt, lz)
!          do ia = tbc%d2amap(pid)+1, tbc%d2amap(pid+1)
!             ixa = indxH(ia)+1
!             nlma = iwk(ia)
!             e(ia) = 0.0_8
!             do la = 0, nlma-1
!                e(ia) = e(ia) + zt(1,ixa+la,ixa+la)
!             end do
!          end do
!       end if
!
!       if (lov) deallocate(sk)
!       if (trh) deallocate(hk0)
!       deallocate(zt)
!       deallocate(zr)
!
!       call tcx('tbfrce')
!
!       end subroutine tbfrce_naive






!                write(u,'("(",2(x,i4),")","[",2(x,i4),"]", &
!                      & x,i4,"[",2(x,i4),"]",x,i5,"[",2(x,i4),"]",x,i5)') &
!                      & lbcrds,gbcrds, pid, c2d%crd(0:1), destrank, remcrd, tag
!       call tcn('msq')
!       if (all(wts(1:nstate) == wts(1))) then
!          zt(1:ldim,1:nstate) = z(1:ldim,1:nstate)*sqrt(wts(1))
!       else
!          do istate = 1, nstate
!             zt(1:ldim,istate) = z(1:ldim,istate)*sqrt(wts(istate))
!          end do
!       end if
!       call tcx('msq')
!
!       call tcn('transpose')
!       zt(1:nstate,1:ldim) = transpose(zt(1:ldim,1:nstate))
!       call tcx('transpose')


!       call tcn('msq')
!       if (all(wts(1:nstate) == wts(1))) then
!          z(1:typesz,1:ldim,1:nstate) = z(1:typesz,1:ldim,1:nstate)*sqrt(wts(1))
!       else
!          do istate = 1, nstate
!             z(1:typesz,1:ldim,istate) = z(1:typesz,1:ldim,istate)*sqrt(wts(istate))
!          end do
!       end if
!       call tcx('msq')

! !          write(u,'(" from lbcrds [gcrds] at pid, [pcrds] to ")')
!                tgbcrd = gbcrds
!                write(u,'("recv ",i4,4(x,"(",2(x,i4),")"),x,i4,x,i5)') pid,c2d%crd(0:1), lbcrds,gbcrds,remcrd,destrank,tag
!           write(u,'(" from lbcrds [gcrds] at pid, [pcrds] to dest [dpcrds], tag")')
! !                sblocks(1:tbc%blcksz(1),1:tbc%blcksz(2),i,j) = transpose(zdist(lcrds(1):lendc(1),lcrds(2):lendc(2)))

!                sblocks(1:typesz, 1:tbc%blcksz(1),1:tbc%blcksz(2),i,j) &
!                & = zdist(1:typesz,lcrds(1)+1:lcrds(1)+tbc%blcksz(1),lcrds(2)+1:lcrds(2)+tbc%blcksz(2))
!
!                ztdist(1:typesz,lcrds(1)+1:lendc(1),lcrds(2)+1:lendc(2)) = sblocks(1:typesz,1:tbc%blcksz(1),1:tbc%blcksz(2),i,j)
!
!                  ztdist(1:typesz,lcrds(1)+1:lendc(1),lcrds(2)+1:lendc(2)) &
!               &=  zdist(1:typesz,lcrds(1)+1:lendc(1),lcrds(2)+1:lendc(2))
!
!                tgbcrd = gbcrds
!
!       write(u+500,'(768(768(x,f6.2),/))') z
!       write(u+300,'(768(768(x,f6.2),/))') zt
!
!       flush(u+300)
!       flush(u+500)

!       call mpi_barrier(c2d%comm, err)

!       flush(u)
!

!       write(u+400,'(768(768(x,f6.2),/))') zt
!       flush(u+400)

!       call mpi_barrier(c2d%comm, err)

!       call rx('stop')
