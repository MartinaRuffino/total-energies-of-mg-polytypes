      subroutine chk_selfcons(ldcix,nsp,ncix,nomg,gloc)
C- Check for SC condition Gloc=Gimp. Optional print Gloc and Gimp on a file.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci Inputs
Ci   ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   nsp    :dimensioning parameter, number of spins
Ci   ncix   :dimensioning parameter, number of correlated species
Ci   nomg   :dimensioning parameter, number of frequency points
Ci   gloc   :local interacting Green's function
Ci   omega  :frequency in Ry
Ci   lprtg  :logical to print out Gimp and Gloc
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... For structures
      integer,    intent(in) :: ldcix,nsp,nomg,ncix
      complex(8), intent(in) :: gloc(ldcix,ldcix,nsp,ncix,nomg)
C ... Local parameters
      real(8), parameter :: Ry2eV = 13.60569193d0
      procedure(integer) :: fopna,fxst
      integer    :: isp,iomg,icix,idev
      real(8)    :: devsc(ldcix*nsp*ncix*nomg)
      integer    :: ichan,i,ifi
      logical    :: lgimp
      complex(8) :: gimp(ldcix,nsp,ncix,nomg)
      real(8)    :: romg, rgimp(ldcix*nsp*2)  ! read from file


      lgimp = fxst('gimp.prev') /= 0
      icix=1

C ... Makes the check only if there is the gimp file
      if(lgimp)then
       call info0(10,1,0,' Checking SC condition (TO DO ncix/=1)...')

!  --- Read gimp of previous iteration
!  --- and compute the deviation from actual gloc
       devsc=0d0 ; idev=0
       ifi=fopna('gimp.prev',-1,1)  ; rewind ifi ! input
       do iomg=1,nomg
        read(ifi,*) romg, rgimp
        do isp=1,nsp
         do ichan=1,ldcix
          i=(isp-1)*ldcix*2+(ichan*2-1)
          gimp(ichan,isp,icix,iomg)=cmplx(rgimp(i),rgimp(i+1))

!     --- Cumulate the deviation at each element
          idev=idev+1
          devsc(idev) = abs(gloc(ichan,ichan,isp,icix,iomg)/Ry2eV-gimp(ichan,isp,icix,iomg))

         enddo
        enddo
       enddo
       call fclose(ifi)
       write(*,'(" Deviation from SC condition Gloc-Gimp=0 (max,avg) [1/eV]"
     .          ,2e14.5)') maxval(devsc), sum(devsc)/size(devsc)

C ... In case the files from the previous iteration are not present.
      else
        call info0(10,1,0,' Check SC condition skipped (missing information about previous iteration)')
      endif

      end subroutine chk_selfcons
