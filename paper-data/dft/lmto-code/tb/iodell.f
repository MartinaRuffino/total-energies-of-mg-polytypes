      subroutine iodelL(ltb,io,nl,nbas,nsp,del,delL)
C- I/O on-site hamiltonian increments for l >= 0
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks:
Cr  If io=0 write the delta's to disc. If io=1 look for delta's on disc
Cr  if not found, set to zero
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer ltb,io,nl,nbas,nsp
      double precision del(0:nl-1,nbas),delL(nl**2,nl**2,nbas,nsp)
C Local Variables
      integer ifi,isp,ib,ilm,fopn,iprint
      logical bittst

      ifi = fopn('DELTA')
      rewind ifi
      if (io == 0) then
        if (nl == 3) then
          write (ifi,300) delL
        else
          write (ifi,310) delL
        endif
      elseif (io == 1) then
        if (.not. bittst(ltb,2**16)) goto 100
        if (nl == 3) then
          read (ifi,300,end=200,err=200) delL
        else
          read (ifi,310,end=200,err=200) delL
        endif
      else
        call rx('iodelL: bad io')
      endif
      return
  200 continue
      if (iprint() > 20) then
        print*, 'IODELL: error reading delta''s from delta file'
        print*, '        setting delta''s to zero ..'
      endif
  100 continue
!   This is incomplete!! The offdiagonal elements are left to chance.
!       do   ib = 1, nbas
!         do  isp = 1, nsp
!           do   ilm = 1, nl**2
!             delL(ilm,ilm,ib,isp) = 0d0
!           enddo
!         enddo
!       enddo
      delL = 0.0_8
  300 format (9f15.9)
  310 format (4f15.9)
      end
