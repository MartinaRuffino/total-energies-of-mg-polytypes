      subroutine cmp_overlaptest(dmftu,nindo,nicix,nlmax,ndham,nev,Olapmk)
C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   dmftu
Ci   nindo
Ci   nicix
Ci   nlmax
Ci   ndham :dimensioning parameter, largest hamiltonian dimension
Ci   nev   :actual number of eigenvectors generated
Ci   iq
Ci   Olapmk
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   19 Oct 14
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: nicix,nlmax,ndham,nev
      complex(8), intent(in) :: dmftu(ndham,nlmax,nicix)
      integer, intent(in) :: nindo(nicix)
      complex(8), intent(out) :: Olapmk(nlmax,nlmax,nicix,nicix)
C ... Local parameters
      integer :: iorb1,iorb2,nind1,nind2
      complex(8) :: olp(nlmax,nlmax)
      character(len=128) :: qlbl
      real(8) :: rqlb
      logical, parameter :: debug=.true.

      rqlb = sqrt(2d0)
      call dpzero(olp,2*size(olp))
      do iorb1 = 1, nicix
        nind1 = nindo(iorb1)
        do iorb2 = 1, nicix
          nind2 = nindo(iorb2)
 !        call zgemm('C','N', nind1, nind2, nbands, (1d0,0d0), DMFTU(1,1,iorb1),nbands, DMFTU(1,1,iorb2),nbands, (0d0,0d0), olp,maxdim2)
          call zgemm('C','N',nind1,nind2,nev,(1d0,0d0),dmftu(1,1,iorb1),ndham,dmftu(1,1,iorb2),ndham,(0d0,0d0),
     .      Olapmk(1,1,iorb1,iorb2),nlmax)
        enddo
      enddo

      if (debug) then
        write(qlbl,'(f12.8)') rqlb
        call zwrmdmft(0,qlbl,'olapmkprint',' ',Olapmk,nlmax,5,5)
      endif

      end subroutine cmp_overlaptest
