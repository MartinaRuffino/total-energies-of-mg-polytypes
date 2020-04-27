      subroutine mkstrx(tbc,ldip,nbas,nlmq1,nlmq,bas,ipc, nclas, lmxl,awld,alat,   &
                      & vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy,pv,MOL,strx,dstrx,plat)
!C- Make lookup table of strux and radial derivatives
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
!Ci        : 2 include 'slab' dipole correction to Ewald
!Ci        : any other number - skip the dipole correction
!Ci   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
!Ci   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
!Ci   nlmq : max L-cutoff for multipoles, leading dimensions of dstrx,
!Ci           second dimension of strx
!Ci   lmxl : l-cutoff for multipoles for each class of species
!Ci   pv   : true if radial derivatives required
!Ci   MOL  : true for molecule (cluster) branch
!Co Outputs:
!Co  strx  : coefficients B_LL'(R'-R) (structure constants)
!Co  dstrx : derivatives of structure constants (x dB_LL'(x)/dx) at x = |R'-R|
!Co          if pv = F dstrx is not touched (see hstra.f)
!Cr Remarks
!Cr   Instead of making B on each sc iteration, the program prepares the whole B
!Cr   in the beginning of the self-consistency cycle to be used as a lookup table.
!Cr   This results in faster calculation at an expence of memory.
!Cr
!Cr   Calling mkstrx is set as the default option. To call instead hstr on each
!Cr   iteration (and therefore to reduce memory but increase computational time),
!Cr   use switch --sfly
!Cr
!Cr   Efficiency issues: there are two symmetry properties of zero energy
!Cr   structure constants that are used in mkstrx:
!Cr     B_LL'(R'-R) = B_L'L(R-R')                (1)
!Cr     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
!Cr   same properties hold for the radial derivative of B.
!Cr
!Cr   According to (1), mkstrx calculates strx(:,:,ib,jb) only for the
!Cr   lower triangle {ib = 1:nbas, jb = 1:ib}
!Cr   property (2) is used in hstra in the same way, ie strx(ilm,jlm,:,:)
!Cr   is sought only for the triangle {ilm = 1:Lmax, jlm=1:ilm}.
!Cr   hstra should therefore be called with Lmax >= Lmax', otherwise it stops.
!Cr
!Cr   mkstrx fills only lower triangles of strx and dstrx in the above sense;
!Cr   properties (1) and (2) are used later, when strx and dstrx are copied into
!Cr   'local' arrays fstrx and fdstrx in tbesel.f We decided to keep it this way
!Cr   because in future we plan to reduce the size of strx and dstrx by cutting of
!Cr   the unused parts of arrays strx and dstrx.
!Cr
!Cr   hstra vs hstr: a call to hstr (a routine that actually makes strux) is
!Cr   replaced with a call to hstra, an 'optimised' version of hstr. The difference
!Cr   between hstr and hstra is that hstra computes only the lower triangle (see above)
!Cr   of strux but is restricted to Lmax >= Lmax'. hstr does not have any restrictions
!Cr   and is called if the --sfly switch is on.
!Cr
!Cb Bugs
!Cb   mkstrx is not written in the most efficient way and could be further refined wrt to
!Cb   amount of calculated elements of strux. At the moment mkstrx is a result of certain
!Cb   trade off between performance and clarity of the program.
!Cb
!Cu Updates
!Cu    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
!Cu                      jb<=ib only)
!Cu    19 Jan 2010 (ATP) first written
!C ----------------------------------------------------------------------
      use tbprl
      use mpi
      implicit none

!C Passed Parameters
      type(tbc_t), intent(in) :: tbc
      logical, intent(in) :: pv,MOL
      integer, intent(in) :: nbas,nlmq,nlmq1,ipc(nbas),ldip,indxcg(*),jcg(*),nkd,nkg,nclas,lmxl(nclas)
      double precision, intent(in)  :: bas(3,nbas),alat,awld,vol,dlat(3,nkd),glat(3,nkg),cy(*),cg(*),plat(3,3)
      double precision, intent(inout) :: strx(nlmq,nlmq1,nbas,*),dstrx(nlmq,nlmq,nbas,*)
!C Local Variables
      integer ib,jb,lmax,lmxf,lmxst,nlx,nlf,i1mach,iprint,ibr
      integer li,li1,lj, ilm,jlm, i
      double precision tau(3),taua(3)
      integer, parameter :: dp = 8
      real(dp) :: blk(nlmq1,nlmq), dblk(nlmq,nlmq), hl(100),bl(100)
      real(dp) :: qlat(3,3), det
      integer :: cmol, cpv
      integer :: nbas1, u
      integer :: pid
      logical :: lroot
! ! #define CUDAM
!#ifdef! CUDAM
!      interface
!!          subroutine rcnsl0(tau,a,lmax,alat,rlat,nkr,dlat,nkd,vol,cy,dl) bind(c, name='rcnsl0_hcu')
!!             use iso_c_binding, only : c_int, c_double
!!             implicit none
!!             integer(c_int), intent(in), value :: lmax, nkd, nkr
!!             real(c_double), intent(in), value :: a,alat,vol
!!             real(c_double), intent(in) :: tau(3), dlat(3,nkd), rlat(3,nkr), cy(*)
!!             real(c_double), intent(out) :: dl(*)
!!          end subroutine
!
!!          void hstr_hcu(bool mol, bool pv, int ldip,
!!                         double *strx, double *drstrx,
!!                         int nlmf, int nlm, int nlmq1, int nlmq,
!!                         double *hl, double *cg, int *indxcg, int *jcg,
!!                         double vol)
!!          subroutine hstra(mol,pv,ldip,strx,drstrx,nlmf,nlm,nlmq1,nlmq,hl,cg,indxcg,jcg,vol) bind(c, name='hstr_c')
!!             use iso_c_binding, only: c_bool, c_int, c_double
!!             implicit none
!!             logical(c_bool), intent(in), value :: mol, pv
!!             integer(c_int), intent(in), value :: ldip, nlmf, nlm, nlmq1, nlmq
!!             integer(c_int), intent(in) :: indxcg(*), jcg(*)
!!             real(c_double), intent(in) :: cg(*), hl(*)
!!             real(c_double), intent(in), value :: vol
!!             real(c_double), intent(out) :: strx(*), drstrx(*) !strx(nlmq1,nlm),drstrx(nlmq,nlm)
!!          end subroutine hstra
!         subroutine mkstrx_c(ldip, nbas, nbas1, ib0, nlmq1, nlmq,    &
!                         & bas, ipc, nclas, lmxl,  awld,  alat,  vol,  &
!                         & dlat, nkd,  glat, nkg,                      &
!                         & indxcg, jcg,  cg,  cy,                      &
!                         & pv, mol,  strx,  dstrx,           &
!                         & plat, qlat) bind(c, name='mkstrx_hcu')
!            use iso_c_binding, only : c_int, c_double
!            implicit none
!            integer(c_int), intent(in), value :: pv, mol
!            integer(c_int), intent(in), value :: ldip, nbas, nbas1, ib0, nlmq1, nlmq, nkd, nkg, nclas
!            integer(c_int), intent(in), dimension(*) :: ipc, lmxl, indxcg, jcg
!            real(c_double), intent(in), value :: awld, alat, vol
!            real(c_double), intent(in), dimension(*) :: bas, dlat, glat, cg, cy, plat, qlat
!            real(c_double), intent(out), dimension(*) :: strx, dstrx
!         end subroutine mkstrx_c
!
!      end interface
!#endif

!       integer :: test_i, err
!       if (pid == 0) then
!          write (*,"('waiting')")
!          read(*,*) test_i
!       end if
!
!       call mpi_barrier(mpi_comm_world, err)


      call tcn('mkstrx')

      lroot  = tbc % c3d % lrt
      pid    = tbc % c3d % id

      if (pid == 0) print *, 'nkg,nkd: ', nkg,nkd
!       print *, 'dlat:'
!       write (857, '(i0,/,"mkstrx")') 2*768*768 + nkd
!       do i = 0, nkd-1
!          print '(3(x,f18.12))', dlat(i*3+1:(i+1)*3)
!          write(857, '("Cl ", 3(x,f18.12))'), dlat(i*3+1:(i+1)*3)
!       end do

      call dinv33(plat,0,qlat,det)
!
!       write(*, '(a,3(/,3(x,f12.6)))') 'plat:', plat
!       write(*, '(a,3(/,3(x,f12.6)))') 'qlat:', qlat
!       write(*, '(a,3(/,3(x,f12.6)))') 'glat:', glat(4:12)
!       stop 'lat'

      print *, 'enter mkstrx'
!#ifdef! CUDAM
!      cmol = 0; cpv = 0
!      if (mol) cmol = 1
!      if (pv)  cpv = 1
!
!      nbas1 = tbc%esamap(pid+1)-tbc%esamap(pid)
!      call mkstrx_c(ldip, nbas, nbas1, tbc%esamap(pid), nlmq1, nlmq,   &
!                         & bas, ipc, nclas, lmxl, awld, alat, vol,         &
!                         & dlat, nkd, glat, nkg,                      &
!                         & indxcg, jcg, cg, cy,                      &
!                         & cpv, cmol, strx, dstrx,           &
!                         & plat, qlat)
!
!
!#else
      do  ib = 1, tbc%esamap(pid+1)-tbc%esamap(pid)
         ibr = ib+tbc%esamap(pid)
         li = lmxl(ipc(ibr))
!c       nlmi = (li+1)**2
!C...  nlmq1 == nlmq if no force or pressure are required.
!C     Otherwise need to go up to lmax+1 in the first index of B_LL'
!          if (nlmq1 > nlmq) then
!             li1 = li + 1
!          else
!             li1 = li
!          endif
         if (nlmq1 > nlmq) li = li + 1
         nlf = (li+1)*(li+1)

         do  jb = 1, nbas
            lj = lmxl(ipc(jb))

!C...  The lines below are because we shall be filling the jb > ib triangle
!C     using symmetry properties of B (see tbesel.f) AND because we need l+1
!C     for forces. No complications if forces are not required.
!             if (nlmq1 > nlmq) then
!                lmax = max(li,lj)
!                lmxf = lmax+1
!                nlx = (lmax+1)*(lmax+1)
!                nlf = (lmax+2)*(lmax+2)
!                lmxst = 2*lmax+1
!             else
!                nlx = (lj+1)*(lj+1)
!                nlf = (li+1)*(li+1)
!                nlf = max(nlf,nlx)
!                lmxst = li + lj
!             endif
!
!             nlx = (lj+1)*(lj+1)
!             nlf = (li+2)**2
            nlx = (lj+1)**2
            lmxst = li + lj
!             lmxst = li+lj+1

!          print "(9(x,a,x,i0))", 'ib:',ib,'li',li,'nlmq1',nlmq1,'nlf',nlf, &
!                      & 'jb',jb,'lj',lj,'nlmq',nlmq,'nlx',nlx,'lmxst',lmxst


            tau = bas(1:3,jb)-bas(1:3,ibr)
            if (.not. mol) then
!                write(*, '(2(x,i3),3(x,f18.12))', advance='no') ib, jb, tau
!                write(857, '("Br ",3(x,f18.12))') tau
!                call shortn(tau,taua,dlat,nkd)
!                taua = tau
!                call directshortn(taua,plat,qlat)
!                call simpleshortn(tau,taua,dlat,nkd)
!                write(*, '(2(x,3(x,f18.12)))', advance='yes') taua, taua-tau
!                write(857, '("I ",3(x,f18.12))') taua
!                tau = taua
!                call simpleshortn(tau,tau,dlat,nkd)

               call directshortn(tau,plat,qlat)
!                call tcn('rcnsl0')
!                hl(1:16) = 0.0_8
               call rcnsl0(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl)
!                write(*, '(3(x,i3),3(x,f18.12))') ib, jb, lmxst, tau
!                write(*, '(2(x,i3),16(x,f18.12))') ib, jb, hl
!                call tcx('rcnsl0')
!                if (jb > 1) stop 'MKSTRX'
            else
               taua = tau * alat
!                write(*,"('f ib,jb: ',i4,x,i4)") ib, jb
               call soldhj(taua,0d0,0,lmxst,hl,bl)
!                write(*,"('hl:',/,20(5(x,f16.8),/))") hl
            endif

            blk = 0.0_dp; dblk = 0.0_dp

            if (iprint() > 60) call awrit2('tau=%3:1d, rstrx=%100:1d',' ',1028,i1mach(2),tau,hl)

!             call hstra(mol,pv,ldip,blk,dblk,nlf,nlx,nlmq1,nlmq,hl,cg,indxcg,jcg,vol)
            call hstr (mol,pv,ldip,blk,dblk,nlf,nlx,nlmq1,nlmq,hl,cg,indxcg,jcg,vol)


            strx(1:nlmq,1:nlmq1,jb,ib) = transpose(blk)

            if (pv) dstrx(1:nlmq,1:nlmq,jb,ib) = transpose(dblk)

!             write(2000+pid,*) 'ib,jb:', ibr, jb
!             do ilm = 1, nlmq1
!                do jlm = 1, nlmq
!                   write(2000+pid,'(f16.8)',advance='no') blk(ilm,jlm)
!                end do
!                write(2000+pid,'(" ")')
!             end do
!          if (jb > 2) stop 'jb>5'
         enddo
      end do


!       open(newunit=u, file='strx-f', action='write')
!
!       do ib = 1, tbc%esamap(pid+1)-tbc%esamap(pid)
!          do jb = 1, nbas
!             write(u,'("jb, ib: ",i4,x,i4)') jb -1, ib - 1
!             do ilm = 1, nlmq1
!                do jlm = 1, nlmq
!                   write(u,'(x,f16.8)',advance='no') strx(jlm,ilm,jb,ib)
!                end do
!                write(u,'()')
!             end do
!          end do
!       end do
!       close(u)
!

!       stop 'mkstrx'

!#endif
      print *, 'end mkstrx'
      call tcx('mkstrx')

!       call rx0("MKSTRX")
      end subroutine mkstrx

