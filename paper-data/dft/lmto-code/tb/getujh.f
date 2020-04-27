      subroutine getujh(s_spec,nl,nbas,ipc,dclabl,idu,uh,jh)
C- Unpack Hubbard U's and J's from spec structure
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sspec,nl,nbas,ipc
Co Outputs:
Co   idu,uh,jh
Cr Remarks
Cr   self explanatory
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   02 Mar 01 Added scissors operator
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,ipc(*),idu(4,nbas)
      double precision dclabl(*),uh(4,nbas),jh(4,nbas)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,ipr,iprint,stdo,nglob
      character*8 clabl

      ipr = iprint()
      stdo = nglob('stdo')
      if (ipr > 30) then
        call awrit0(' GETUJH: unpacking U and J ..',' ',120,stdo)
      endif
      do  ib = 1, nbas
        is = ipc(ib)
        idu(1:4,ib) = s_spec(is)%idu(1:4)
        uh(1:4,ib) = s_spec(is)%uh(1:4)
        jh(1:4,ib) = s_spec(is)%jh(1:4)
        if (ipr > 30) then
          call r8tos8(dclabl(is),clabl)
          call awrit7('Atom %#5i '//clabl//' IDU=%n:2i U=%n:2d J=%n:2d',
     .      ' ',256,stdo,ib,nl,idu(1,ib),nl,uh(1,ib),nl,jh(1,ib))
        endif
      enddo
      end
