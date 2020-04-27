    subroutine findmolidc(nbas, iax, niax, nsites, npr, molidc)
    implicit none
! Find the molecule index of each atom using the neighbourtable.
! As this is first stab, it scraps the old indices (if available) and starts anew. It will not be too complicated to only update the present list in later stage

    integer, intent(in) :: nbas, niax, nsites, iax(0:niax-1,nsites), npr(0:1,nbas)
    integer, intent(inout) :: molidc(nbas)

    integer, parameter :: maxchain = 128
    integer :: i, j, ij, nr, midx, ch(nbas)
    logical :: deadend

    call tcn('Find mol idcs')

    midx = 1
    do i = 1, nbas
        if (molidc(i) /= 0) cycle
        molidc(i) = midx
        nr = 1
        ch(nr) = i
        do ! recursion to find all atoms of the molecule containing atom i
            deadend = .true.
!           Probably worth it to save the starting ending position as well instead of repeating the already walked neighbours just to discover that they are set.
!           Also will be useful to filter the iax/npr table through the specific bondtype
!           lengths and remove any repeating neighbours for diff translations, which are kept in iax. May be tbc%a2map could be used instead of iax since.
            do ij = 1, npr(0, ch(nr))
                j = iax(1,npr(1, ch(nr))+ij)
                if (molidc(j) == 0) then
!                     print *, 'ch:', ch(:nr), '[',j,']'
                    deadend = .false.
                    molidc(j) = midx
                    nr = nr + 1
                    if (nr > maxchain) &
                        call rx('FINDMOLIDC: There is a path, longer than envisioned (128)! Increase "maxchain" accordingly!')
                    ch(nr) = j
                    exit
                end if
            end do
            if (deadend) then
                if (nr == 1) exit
                nr = nr - 1
            end if
        end do
        midx = midx + 1
    end do

    do i = 1, nbas
        print *, i, molidc(i)
    end do

    call tcx('Find mol idcs')

    call rx('molecules found and printed! stopping nowi to compare')

    end subroutine findmolidc
