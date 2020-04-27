subroutine mkstrxidx(tbc, pv, force, nbas, ipc, lmxl, ssize, dsize, struxidx, dstrxidx)
!- Make structure constants index arrays
!--------------------------------------------------------------------------------
!i Inputs:
!i   tbc%esamap  : Offsets for the structure constants distribution. A contiguous block is given to each process.
!i   pv          : whether to make the radial derivatives and go to a higher harmonic in strux
!i   force       : whether to go to a higher harmonic in strux
!i   nbas        : number of atoms
!i   ipc         : atom classes
!i   lmxl        : Lmax for each atom class
!
!o Outputs:
!o   ssize       : size of the compact structure constants array in units of real(8)
!o   dsize       : size of the compact derivatives of the structure constants in units of real(8)
!o   struxidx    : index for struxd, shall be allocated to (nbas,tbc%escount(pid)),
!o                 tbc%escount(pid) the number of atoms allocated to process 'pid'
!o   dstrxidx    : index for dstrxd
!
!r Remark:
!r The sizes of each block are calculated and local offsets generated.
!r The indexing is transposed for cache friendliness and also ease of distribution in the parallel case
!r Neither the indices nor the structure constants are ever shared.
!r There is a 'real and present' danger the default integer(4) overflows for largeish systems,
!r  If this happens switching to integer(8) shall help but I guess there
!r  is a lot more to wory about in other places so I am leaving to int4 for now.


   use tbprl
   implicit none

   type(tbc_t), intent(in) :: tbc

   integer, intent(in) :: lmxl(*), ipc(*), nbas
   integer, intent(out) :: struxidx(nbas,*), dstrxidx(nbas,*), ssize, dsize
   logical :: pv, force
   procedure(logical) :: bittst
   integer :: pid
   logical :: lroot
   integer :: sidx, didx, pib, ib, jb, li, lj, nlmi, nlmi1, nlmj

   call tcn('mkstrxidx')

   lroot  = tbc % c3d % lrt
   pid    = tbc % c3d % id

   sidx = 1
   didx = 1

   do  pib = 1, tbc%esamap(pid+1)-tbc%esamap(pid)
      ib = pib + tbc%esamap(pid)
      li = lmxl(ipc(ib))
      nlmi  = (li+1)*(li+1)

      if (force .or. pv) then
         nlmi1  = (li+2)*(li+2)
      else
         nlmi1  = nlmi
      endif

      do  jb = 1, nbas
         lj = lmxl(ipc(jb))
         nlmj = (lj+1)*(lj+1)

         struxidx(jb,pib) = sidx
         sidx = sidx + nlmj*nlmi1
         if (sidx < struxidx(jb,pib)) call rx('MKSTRIDX: Integer overflow. Converto to int8 or address type')

         if (pv) then
            dstrxidx(jb,pib) = didx
            didx = didx + nlmj*nlmi
         end if
      end do
   end do

   ssize = sidx
   dsize = didx

   call tcx('mkstrxidx')

end subroutine mkstrxidx
