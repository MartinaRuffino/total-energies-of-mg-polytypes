subroutine mkstrxd(s_ctrl, ipc, s_lat, tbc, nlmq, nlmq1, lmxl, ldip, dlat, nkd, glat, nkg, &
                     & indxcg, jcg, cg, cy, struxd, dstrxd, struxidx, dstrxidx, struxsize, dstrxsize)
!- Make lookup table of strux and radial derivatives
! ----------------------------------------------------------------------
!i Inputs:
!i  ldip  : 3 include 'spherical' dipole correction to Ewald
!i        : 2 include 'slab' dipole correction to Ewald
!i        : any other number - skip the dipole correction
!i   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
!i   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
!i   nlmq : max L-cutoff for multipoles, leading dimensions of dstrx,
!i           second dimension of strx
!i   lmxl : l-cutoff for multipoles for each class of species
!i   pv   : true if radial derivatives required
!i   MOL  : true for molecule (cluster) branch
!i   struxidx    : index for struxd, shall be allocated to (nbas,tbc%escount(pid)),
!i                 tbc%escount(pid) the number of atoms allocated to process 'pid'
!i   dstrxidx    : index for dstrxd
!
!o Outputs:
!o  struxd : coefficients B_LL'(R'-R) (structure constants)
!o  dstrxd : derivatives of structure constants (x dB_LL'(x)/dx) at x = |R'-R|
!o          if pv = F dstrx is not touched
!r Remarks
!r   Instead of making B on each sc iteration, the program prepares the whole B
!r   in the beginning of the self-consistency cycle to be used as a lookup table.
!r   This results in faster calculation at an expence of memory.
!r
!r   Calling mkstrx is set as the default option. To call instead hstr on each
!r   iteration (and therefore to reduce memory but increase computational time),
!r   use switch --sfly
!r
!r   Efficiency issues: there are two symmetry properties of zero energy
!r   structure constants that are used in mkstrx:
!r     B_LL'(R'-R) = B_L'L(R-R')                (1)
!r     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
!r   same properties hold for the radial derivative of B.
!r
!r This symmetry however is not employed because the complexity that a
!r  parallel implementation will require and the diminishing gains + increased overhead with.
!r Curently the matrices and the indices are distributed in contiguous blocks across the long
!r side and no communication whatsoever is needed. For a process with id pid, the block starts
!r at atom tbc%esamap(pid)+1 and ends at tbc%esamap(pid+1) inclusive. For atom ib from this range
!r there is a column of neighbours in struxidx(:,ib-tbc%esamap(pid)). For each of pair (jb,ib)
!r a block of size (nlmj,nlmi1) is preallocated starting at struxd(struxidx(jb,ib-tbc%esamap(pid)))
!r The block size for dstrx is (nlmj,nlmi)
!r
!b Bugs
!b   mkstrx is not written in the most efficient way and could be further refined wrt to
!b   amount of calculated elements of strux. At the moment mkstrx is a result of certain
!b   trade off between performance and clarity of the program.
!b
!u Updates
!u           2013 (DP)  strx moved to indexed linear array with auxiliary index for better compatibility with C
!u           2013 (DP)  Fortran pointer based strx for memory efficiency
!u           2012 (DP)  Parallell distributed strx & dstrx
!u    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
!u                      jb<=ib only)
!u    19 Jan 2010 (ATP) first written
! ----------------------------------------------------------------------

   use tbprl
   use structures
   implicit none
   integer, parameter :: dp = 8

   type(str_ctrl ), intent(in) :: s_ctrl
   type(str_lat  ), intent(in) :: s_lat
   type(tbc_t), intent(in) :: tbc

   integer, intent(in) :: lmxl(*), ldip, nkd, nkg, indxcg(*),jcg(*),ipc(*), struxsize, dstrxsize
   real(dp), intent(in) :: dlat(3,nkd), glat(3,nkg), cy(*),cg(*)
   integer, intent(in) :: struxidx(s_ctrl%nbas,*), dstrxidx(s_ctrl%nbas,*)
   real(dp), intent(out) :: struxd(struxsize), dstrxd(dstrxsize)


   integer :: pib,jb,lmax,lmxf,lmxst,nlmj,nlmi,nlmi1,i1mach,iprint,ib,nlmq,nlmq1
   integer :: li,li1,lj, ilm,jlm, i, nbas1, u, pi, pj, sz, dsz, lmaxl
   real(dp) :: tau(3),taua(3), hl(100),bl(100), plat(3,3), qlat(3,3), alat, &
                              & vol, awald, det, struxl(nlmq*nlmq1), dstrxl(nlmq*nlmq)
   logical :: pv, mol, force
   procedure(logical) :: bittst
   integer :: pid
   integer :: cmol, cpv, cforce
   logical :: lroot

   interface
      subroutine mkstrx_c(ldip, nbas, nbas1, ib0, lmaxl,    &
                        & bas, ipc, nclas, lmxl,  awld,  alat,  vol,  &
                        & dlat, nkd,  glat, nkg,                      &
                        & indxcg, jcg,  cg,  cy,                      &
                        & pv, mol, force, strx,  dstrx,  struxidx, dstrxidx, struxsize, dstrxsize,  &
                        & plat, qlat) bind(c, name='mkstrxd_hcu')
         use iso_c_binding, only : c_int, c_double
         implicit none
         integer(c_int), intent(in), value :: pv, mol, force
         integer(c_int), intent(in), value :: ldip, nbas, nbas1, ib0, lmaxl, nkd, nkg, nclas, struxsize, dstrxsize
         integer(c_int), intent(in), dimension(*) :: ipc, lmxl, indxcg, jcg, struxidx, dstrxidx
         real(c_double), intent(in), value :: awld, alat, vol
         real(c_double), intent(in), dimension(*) :: bas, dlat, glat, cg, cy, plat, qlat
         real(c_double), intent(out), dimension(*) :: strx, dstrx
      end subroutine mkstrx_c
   end interface




   call tcn('mkstrxp')

   lroot  = tbc % c3d % lrt
   pid    = tbc % c3d % id

   pv  = bittst(s_ctrl%ltb,128)
   mol = bittst(s_ctrl%ltb,2**18)
   force = bittst(s_ctrl%ltb,16)

   plat  = s_lat%plat
   alat  = s_lat%alat
   vol   = s_lat%vol
   awald = s_lat%awald


   cmol = 0; cpv = 0; cforce = 0
   if (mol) cmol = 1
   if (pv)  cpv = 1
   if (force) cforce = 1

   call dinv33(plat,0,qlat,det)

   if (.not. pv) pj = 1

   nbas1 = tbc%escount(pid)
   lmaxl = nint(sqrt(real(nlmq1,8))) - 1
   if (nlmq1>nlmq) lmaxl = lmaxl - 1

   call mkstrx_c(ldip, s_ctrl%nbas, nbas1, tbc%esamap(pid), lmaxl,    &
               & s_lat%pos, ipc, s_ctrl%nclass, lmxl, awald, alat, vol,    &
               & dlat, nkd, glat, nkg,                       &
               & indxcg, jcg, cg, cy,                        &
               & cpv, cmol, cforce, struxd, dstrxd, struxidx, dstrxidx, struxsize, dstrxsize,&
               & plat, qlat)

   call tcx('mkstrxp')
end subroutine mkstrxd