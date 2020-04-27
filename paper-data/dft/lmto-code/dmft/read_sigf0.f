       subroutine read_sigf0(ndsig,nicix,sigf0)
       ! the sig.inp.f0 file is read in Ry
       implicit none
       real, parameter :: ry2ev=13.60569193d0
       integer ,   intent(in)  :: ndsig,nicix
       real(8), intent(out) :: sigf0(2,ndsig,nicix)
       real(8)  :: rsig(ndsig*2)
       integer  :: l,uin,icix

       open(unit=uin,file='sig.inp.f0')
       read(uin,*) rsig
       rsig = rsig/ry2ev    ! convert in Ry
       do icix=1,nicix
        do l=1,ndsig
         sigf0(1,l,icix) = rsig(l*2-1)
         sigf0(2,l,icix) = rsig(l*2  )
        enddo
       enddo
       close(uin)
       end subroutine read_sigf0
