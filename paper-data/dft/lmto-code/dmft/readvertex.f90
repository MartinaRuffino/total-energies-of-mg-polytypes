    program main
    use vertex, only : read_vertex
    implicit none
    integer :: mode, nvfl, nomv, nOm, norb, nargs
    real(8) :: beta
    character(len=20):: infile, outfile
    nargs=command_argument_count()
    if(nargs /= 2) then
       write(*,*) 'readvertex inputfile outputfile'
       write(*,*)('wrong number of argument')
       stop
    else
       call get_command_argument(1,infile)
       call get_command_argument(2,outfile)
    endif

    call read_vertex(mode,beta,nvfl,nomv,nOm,norb,infile,outfile)


    end program main
