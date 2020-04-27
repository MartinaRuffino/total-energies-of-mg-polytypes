       subroutine print_dmftu(file,fmt,dmft,nindo,ndham,nlmax,nicix,iq,qp)
       implicit none
       character :: fmt
       character(len=*), intent(in) :: file
       integer, intent(in) :: ndham, nlmax, nicix, iq
       integer, intent(in) :: nindo(nicix)
       complex(8), intent(in) :: dmft(ndham,nlmax,nicix)
       real(8), intent(in) :: qp(3)

       integer :: u, iband, ind, ib, fopna

       ib = 0; if (fmt == 'b') ib = 4
       u = fopna(trim(file),-1,ib)  ! Use extension

       if (fmt == 'a') then
           write(u,'(a,1x,i0)') '# DMFT nicix: ', nicix
           write(u,'("# kpt n.",i0,1x,3(1x,f12.8))') iq, qp
       end if

       do iband = 1, ndham
         if (fmt == 'a') write(u,'(1x,i5)',advance='no') iband
         do ib = 1, nicix
           do ind = 1, nindo(ib)
             if (fmt == 'a') then
               write(u,'(2x,2(1x,f12.8))', advance='no')
     .           dble (dmft(iband,ind,ib)),aimag(dmft(iband,ind,ib))
             else if (fmt == 'b') then
               write(u) dble (dmft(iband,ind,ib)),
     .           aimag(dmft(iband,ind,ib))
             end if
           end do
         end do
         if (fmt == 'a') write(u,'("")')
       end do

C      call fclose(u)
       end subroutine print_dmftu
