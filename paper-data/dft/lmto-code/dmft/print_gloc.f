      subroutine print_gloc(s_dmft,gpmode,nomg,ldcix,nsp,ncix,omega,gloc,nomgp)
C- Print Gloc and Gimp on a file. Optional check for convergence.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   nomg   :dimensioning parameter, number of frequency points
Ci   gloc   :local interacting Green's function
Ci   omega  :frequency in Ry
Ci   nomgp  :number of frequencies to print
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C     ----------------------------------------------------------------------
      use structures
      implicit none
      type(str_dmft)      ::  s_dmft
C ... Passed parameters
      integer, intent(in) :: ldcix,nsp,nomg,ncix,nomgp,gpmode
      real(8), intent(in) :: omega(nomg)
      complex(8), intent(in)  :: gloc(ldcix,ldcix,nsp,ncix,nomg)
C ... Local parameters
      real(8), parameter :: ry2eV = 13.60569193d0
      integer :: isp,iomg,icix,ichan,jchan
      integer :: ifo,isig,i,j,cix
      character(len=256) :: strn,fname
      integer, allocatable :: nicixi(:)
      procedure(integer) :: fopna
      procedure(logical) :: cmdopt
      allocate(nicixi(s_dmft%ncix))    ! Loop over inequivalent cix only
      call ineqcix(s_dmft,ncix,nicixi)

      if(gpmode == 5) then
         do icix = 1, ncix
            call awrit1('%x_%i',strn,len(strn),0,icix)
            ifo=fopna('gloc'//strn,-1,2) ; rewind ifo
            do  iomg = 1, nomgp
               write(ifo,'(f14.8)',advance='no') omega(iomg)*ry2ev
               do isp=1,nsp
                  do ichan=1,ldcix
                     if (isp == 2) then
                        do jchan = 1, ldcix
                           write(ifo,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                        end do
                     end if

                     do jchan=1,ldcix
                        write(ifo,'(2(x,f14.8))',advance='no') gloc(jchan,ichan,isp,icix,iomg)/ry2ev
                     end do

                     if (isp < nsp) then
                        do jchan = 1, ldcix
                           write(ifo,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                        end do
                     end if
                  enddo
               enddo

               write(ifo,*)
            enddo               ! nomgp
            call fclose(ifo)

      enddo
      else
         ifo=fopna('gloc',-1,2) ; rewind ifo
         do icix = 1, ncix
            if (.not. cmdopt('--fullg',7,0,strn)) then
               do  iomg = 1, nomgp
                  write(ifo,'(f14.8)',advance='no') omega(iomg)*ry2ev
                  do isp=1,nsp
                     do ichan=1,ldcix
                        write(ifo,'(2(x,f14.8))',advance='no') gloc(ichan,ichan,isp,icix,iomg)/ry2ev
                     enddo
                  enddo
                  write(ifo,*)
               enddo
            else
               do  iomg = 1, nomgp
                  write(ifo,'(f14.8)',advance='no') omega(iomg)*ry2ev
                  do isp=1,nsp
                     do ichan=1,ldcix
                        if (isp == 2) then
                           do jchan = 1, ldcix
                              write(ifo,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                           end do
                        end if

                        do jchan=1,ldcix
                           write(ifo,'(2(x,f14.8))',advance='no') gloc(jchan,ichan,isp,icix,iomg)/ry2ev
                        end do

                        if (isp < nsp) then
                           do jchan = 1, ldcix
                              write(ifo,'(2(x,f14.8))',advance='no') 0.0_8, 0.0_8
                           end do
                        end if
                     enddo
                  enddo
                  write(ifo,*)
               enddo
            end if
         enddo
      endif
      call fclose(ifo)
      call info0(10,0,0,' gloc(N)   recorded in gloc file.')

      end subroutine print_gloc
