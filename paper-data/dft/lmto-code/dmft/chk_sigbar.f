      subroutine chk_sigbar(s_dmft,nicix,orbmax,ndham,nev,nsp,nindo,
     .                      dmftu,sigbarinftyk,projsbark)
C- Compute the projection of sigbar and compares with sig.inp and sig.dc (eq. 5.15 of Sponza's notes)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci Inputs
Ci   dmftu  :projector to local correlated Hilbert space
Ci   ndham  :dimensioning parameter, largest hamiltonian dimension
Ci   nlohi  :parameter controlling first and last band to include
Ci   orbmax :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   nsp    :dimensioning parameter,2 for spin-polarized case, otherwise 1
Ci   nindo  :number of orbitals per channel
Ci   nicix  :dimensioning parameter, number of correlated l/site channels
Ci   sigbarinftyk :high-frequency limit of sigbar in lattice representation {ijk}
Cr Remarks
Cu Updates
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu    9 Nov 15  L.Sponza created
C ----------------------------------------------------------------------
      use structures
      implicit none
!      include 'structures.h'
      type(str_dmft)            :: s_dmft
C ... Passed parameters
      integer,    intent(in)    :: nicix,orbmax,ndham,nev,nsp
      integer,    intent(in)    :: nindo(nicix)
      complex(8), intent(in)    :: dmftu(ndham,orbmax,nicix)
      complex(8), intent(in)    :: sigbarinftyk(nev,nev,nsp)
      complex(8), intent(out)   :: projsbark(orbmax,orbmax,nicix,nicix)
C ... Local parameters
      complex(8),allocatable :: tmp1(:,:)
      integer    :: isp,iorb1,nind1,iorb2,nind2,nl1,nl2
C ... External calls
      external dpzero,zgemm

      call rx('chk_sigbar: check new argument list')

C --- Local projection of sigbar(infty) ---
      call dpzero(projsbark,2*size(projsbark))
      allocate(tmp1(orbmax,nev))
      do  isp = 1, nsp
        do  iorb1 = 1, nicix
          nind1 = nindo(iorb1)
          do  iorb2 = 1, nicix
            nind2 = nindo(iorb2)
            call dpzero(tmp1,2*size(tmp1))
            nl2 = nind2/nsp
            nl1 = nind1/nsp
C       ... Local Eimp ... zgemm call generating tmp1 : a very expensive way to make eval * dmftu ... rewrite as level 2 BLAS
            call zgemm('C','N',nl1,nev,nev,(1d0,0d0),
     .        dmftu(1,1+(isp-1)*nl1,iorb1),ndham,
     .        sigbarinftyk(1,1,isp),nev,(0d0,0d0),
     .        tmp1,orbmax)
            call zgemm('N','N',nl1,nl2,nev,(1d0,0d0),
     .        tmp1,orbmax,
     .        dmftu(1,1+(isp-1)*nl2,iorb2),ndham,(0d0,0d0),
     .        projsbark(1+(isp-1)*nl1,1+(isp-1)*nl2,iorb1,iorb2),orbmax)
          enddo
        enddo
      enddo
      deallocate(tmp1)
      end subroutine chk_sigbar

      subroutine chk_sigbar_prt(s_dmft,orbmax,nicix,nindo,ndsig,
     .                          sbarfty,sinpfty,sdcfty,mode)
C- Print files with checks on P+E(sigbar)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci Inputs
Ci   orbmax :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   nindo  :number of orbitals per channel
Ci   nicix  :dimensioning parameter, number of correlated l/site channels
Ci   ndsig  :maximum number of matrix elements in any one DMFT block
Ci   sbarfty :high-frequency limit of sigbar in local representation {LL'}
Ci   sinpfty :high-frequency limit of sig.inp in local representation {LL'}
Ci   sdcfty  :high-frequency limit of sig.dc in local representation {LL'}
Ci   mode   :file produced: (1) File 1 , (2) File 2, (3) both
Cr Remarks
Cu Updates
Cu   10 Nov 15  L.Sponza created
C ----------------------------------------------------------------------
      use structures
      implicit none
!      include 'structures.h'
      type(str_dmft)         :: s_dmft
C ... Passed parameters
      integer,    intent(in) :: orbmax,nicix,mode,ndsig
      integer,    intent(in) :: nindo(nicix)
      complex(8), intent(in) :: sbarfty(orbmax,orbmax,nicix,nicix)
      complex(8), intent(in) :: sdcfty(orbmax,orbmax,nicix,nicix)
      real(8),    intent(in) :: sinpfty(2*ndsig,nicix)
      !internal
      integer, parameter :: u1=8989 , u2=8990
      integer            :: i,il,jl,it
      complex(8)         :: sinp,sdc(ndsig,nicix)

      call rx('chk_sigbar_prt: redimension sigdc')

      ! File 1 : full matrix P+E(Sigbar(infty))_{LL'}
      if(mode==1.or.mode==3)then
        open(unit=u1,file='chk_projsbar_full.dat')
        do i=1,nicix
          write(u1,'(" ++++ icix = ", i4)') i
          do il=1,nindo(i)
            do jl=1,nindo(i)
              write(u1,'(2f15.8)',advance='no') sbarfty(il,jl,i,i)
            enddo
            write(u1,'("")')
          enddo
        enddo
        close(u1)
      endif



      ! File 2 : comparison of P+E(Sigbar)), Sig.inp and SigDC in the high-frequency limit
      if(mode==2.or.mode==3)then
        ! first extracts only diagonal values of sigDC
        call dpzero(sdc,2*size(sdc))
        do i = 1, nicix
          do il=1,nindo(i)
            do jl=1,nindo(i)
              it = abs(s_dmft%sigind(il,jl,i,1))
              if (it > 0) then
                sdc(it,i)=sdcfty(il,jl,i,i)
              endif
            enddo
          enddo
        enddo

        ! then writes the file
        open(unit=u2,file='chk_projsbar_comp.dat')
        write(u2,'(5x,"__ P+E(Sbar(infty)) __",5x,
     .             7x,"__ Sinp(infty) __",8x,
     .             8x,"__ Sdc(infty) __",8x,
     .             2x," DEVIATION FROM P+E=1")')
        do i=1,nicix
          write(u2,'(" icix ", i2)') i
          do il=1,nindo(i) !dsig
            sinp=dcmplx( sinpfty(2*il-1,i) , sinpfty(2*il,i) )
!            write(999,'(3i4,3x,2f10.5)')il,ndsig,nindo(i),sinp
           !write(u2,'(3(2f15.8,2x),e20.12)') sbarfty(il,il,i,i),sinp,sdc(il,i),
           write(u2,'(3(2f15.8,2x),e20.12)') sbarfty(il,il,i,i),sinp,sdc(il,i),
     .                          abs(sbarfty(il,il,i,i)-sinp+sdc(il,i))
          enddo
        enddo
        close(u2)
      endif

      end subroutine chk_sigbar_prt
