!       subroutine diagc_bc(H2p,Hmax,eigvalH)
!
!       implicit none
!
!       integer :: Hmax
!       complex*16 :: H2p(Hmax,Hmax),eigvalH(Hmax)
!       integer :: info,LDVL,LDVR,LWORK ! for diagonalising
!       real(8),allocatable :: RWORK(:),evals(:)
!
!       complex*16,allocatable:: work(:), VL(:,:)! maybe complex?
! C
! C      use zgeev to calculate eigvals and right eigvecs? or dgeev????
! C
!       LDVL=Hmax
!       LDVR= Hmax
!       LWORK= 2*Hmax ! if using ZGEEV then LWORK=2*Max
!       allocate(VL(LDVL,Hmax))
!       allocate(WORK(LWORK))
!       allocate(RWORK(2*Hmax)) ! if using ZGEEV
! C
! CC      Do we need to copy H2p before passing it????
!       call ZGEEV('N','V',Hmax,H2p,Hmax,eigvalH,VL,
!      &                          LDVL,eigvecH,LDVR,WORK,LWORK,RWORK,info)
!
! CC       eigvalH = eigvalH_real! + zimag * eigvalH_imag ! need, as complex?
! CC
! CC      For dgeev: if eigenvalue is complex, its neighbour will be
! CC      conjugate.  The eigenvectors(columns) are expressed differently
! CC      (see description in the lapack routine).  If the jth eigenvalue
! CC      is complex,j+1th its conjugate; the jth column contains the real
! CC      part of the eigenVECTOR for both j and j+1. The j+1th column
! CC      contains the imaginary part, which is added to jth eigvec and
! CC      subtracted from the j+1th eigvec.
!
!       deallocate(VL)
!       deallocate(WORK)
!       deallocate(RWORK)
!
!       end

      subroutine diagc_bcHD(H2p,Hmax,eigvalH,ind1)
      use keyvalue,only: getkeyvalue

      implicit none

      integer :: Hmax,ind1
      complex*16 :: H2p(Hmax,Hmax),eigvalH(Hmax)
      integer :: info,LWORK,lrwork,liwork ! for diagonalising
      real(8),allocatable :: RWORK(:),evals(:),eigvalHreal(:)
      integer, allocatable :: iWORK(:)
      complex*16,allocatable :: work(:)
      real(8) :: max_omega,min_omega

       call getkeyvalue("GWinput","MaxOmega_diag",max_omega,
     & default=1d90)     ! In Rydbergs?
      call getkeyvalue("GWinput","MinOmega_diag",min_omega,
     & default=0.0d0)

      if(max_omega<1d10.or.min_omega /= 0.0d0.or.Hmax>32700)
     . then; call diagc_bcHR(H2p,Hmax,eigvalH,ind1)
      return; endif
C
C      use zgeev to calculate eigvals and right eigvecs? or dgeev????
      LWORK= 2*Hmax!-1 ! if using ZGEEV then LWORK=2*Max
      LWORK= -1 ! if using ZGEEV then LWORK=2*Max
      allocate(WORK(2*Hmax+Hmax**2))
      LRWORK=-1
      allocate(RWORK(1+5*Hmax+2*Hmax**2)) ! if using ZGEEV
      LIWORK =-1
      allocate(IWORK(3+5*Hmax))
      allocate(eigvalHreal(Hmax))
C
CC      Do we need to copy H2p before passing it????
      call ZHEEVD('V','U',Hmax,H2p,Hmax,eigvalHreal,WORK,LWORK,RWORK,
     .           LRWORK,IWORK,LIWORK,info)
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        WRITE(*,*)'minimum values for LWORK,LRWORK,LIWORK'
        WRITE(*,*)2*Hmax+Hmax**2,1+5*Hmax+2*Hmax**2,3+5*Hmax
        WRITE(*,*)'suggested values + info'
        write(*,*)LWORK,LRWORK,LIWORK,info
        deallocate(WORK,RWORK,IWORK)
        allocate(WORK(LWORK),RWORK(LRWORK),IWORK(LIWORK))
      call ZHEEVD('V','U',Hmax,H2p,Hmax,eigvalHreal,WORK,LWORK,RWORK,
     .           LRWORK,IWORK,LIWORK,info)

      eigvalH = eigvalHreal
CC       eigvalH = eigvalH_real! + zimag * eigvalH_imag ! need, as complex?
CC
CC      For dgeev: if eigenvalue is complex, its neighbour will be
CC      conjugate.  The eigenvectors(columns) are expressed differently
CC      (see description in the lapack routine).  If the jth eigenvalue
CC      is complex,j+1th its conjugate; the jth column contains the real
CC      part of the eigenVECTOR for both j and j+1. The j+1th column
CC      contains the imaginary part, which is added to jth eigvec and
CC      subtracted from the j+1th eigvec.

      deallocate(WORK)
      deallocate(RWORK,eigvalHreal,IWORK)

      end

      subroutine diagc_bcHR(H2p,Hmax,eigvalH,ind1)
      use keyvalue,only: getkeyvalue

      implicit none

      integer :: Hmax,Mint,ind1
      complex*16 :: H2p(Hmax,Hmax),eigvalH(Hmax)
      integer :: info,LWORK,lrwork,liwork,IL,IU ! for diagonalising
      real(8),allocatable :: RWORK(:),evals(:),eigvalHreal(:)
      integer, allocatable :: iWORK(:),ISUPPZ(:)
      complex*16,allocatable :: work(:),Zmat(:,:)
      real(8) :: VL,VU,ABSTOL,DLAMCH,min_omega,max_omega
      character(1) :: charinput
C
C      use zgeev to calculate eigvals and right eigvecs? or dgeev????
      LWORK= 2*Hmax!-1 ! if using ZGEEV then LWORK=2*Max
      allocate(WORK(LWORK))
      LWORK=-1
      LRWORK=24*Hmax
      allocate(RWORK(LRWORK)) ! if using ZGEEV
      LRWORK=-1
      LIWORK =10*Hmax
      allocate(IWORK(LIWORK))
      LIWORK=-1
      allocate(eigvalHreal(Hmax))
      allocate(Zmat(Hmax,Hmax))
      ABSTOL=DLAMCH( 'Safe minimum' )
        write(*,*)'ABSTOL=',ABSTOL
      allocate(ISUPPZ(2*Hmax))

       call getkeyvalue("GWinput","MaxOmega_diag",max_omega,
     & default=1d90)     ! In Rydbergs?
      call getkeyvalue("GWinput","MinOmega_diag",min_omega,
     & default=0.0d0)
      charinput='V'
      if(max_omega>1d10)charinput='A'
      if(ind1==1)then
      VL=min_omega; VU=max_omega
      else
      VL=-max_omega; VU=-min_omega
      endif
      IL=1; IU=Hmax

      call getkeyvalue("GWinput","IU_diag",IU,
     & default=Hmax)     ! In Rydbergs?
      call getkeyvalue("GWinput","IL_diag",IL,
     & default=1)
        if(IU /= Hmax.or.IL /= 1)charinput='I'
C
CC      Do we need to copy H2p before passing it????
      write(*,*)charinput,'VL,VU,IL,IU=',VL,VU,IL,IU
      call ZHEEVR('V',charinput,'U',Hmax,H2p,Hmax,VL,VU,IL,IU,ABSTOL,
     .        Mint,eigvalHreal,Zmat,Hmax,ISUPPZ,WORK,LWORK,RWORK,
     .           LRWORK,IWORK,LIWORK,info)
C        write(*,*)'Number of eigenvalues found=',Mint
        LWORK=WORK(1)
        LRWORK=RWORK(1)
        LIWORK=IWORK(1)
        write(*,*)'LWORK,LRWORK,LIWORK=',2*Hmax,24*Hmax,10*Hmax
        write(*,*)'LWORK,LRWORK,LIWORK=',LWORK,LRWORK,LIWORK
        deallocate(WORK,RWORK,IWORK)
        allocate(WORK(LWORK),RWORK(LRWORK),IWORK(LIWORK))

      call ZHEEVR('V',charinput,'U',Hmax,H2p,Hmax,VL,VU,IL,IU,ABSTOL,
     .        Mint,eigvalHreal,Zmat,Hmax,ISUPPZ,WORK,LWORK,RWORK,
     .           LRWORK,IWORK,LIWORK,info)
      write(*,*)'number of eigenvalues found=',Mint
      eigvalH = eigvalHreal
CC       eigvalH = eigvalH_real! + zimag * eigvalH_imag ! need, as complex?
CC
CC      For dgeev: if eigenvalue is complex, its neighbour will be
CC      conjugate.  The eigenvectors(columns) are expressed differently
CC      (see description in the lapack routine).  If the jth eigenvalue
CC      is complex,j+1th its conjugate; the jth column contains the real
CC      part of the eigenVECTOR for both j and j+1. The j+1th column
CC      contains the imaginary part, which is added to jth eigvec and
CC      subtracted from the j+1th eigvec.
      H2p=Zmat
      deallocate(WORK,Zmat)
      deallocate(RWORK,eigvalHreal,IWORK,ISUPPZ)

      end


      subroutine diagc_bcH(H2p,Hmax,eigvalH)

      implicit none

      integer :: Hmax
      complex*16 :: H2p(Hmax,Hmax),eigvalH(Hmax)
      integer :: info,LWORK,iwork ! for diagonalising
      real(8),allocatable :: RWORK(:),evals(:),eigvalHreal(:)

      complex(8), allocatable:: work(:)! maybe complex?
C
C      use zgeev to calculate eigvals and right eigvecs? or dgeev????
      LWORK= 2*Hmax!-1 ! if using ZGEEV then LWORK=2*Max
      allocate(WORK(LWORK))
      allocate(RWORK(3*Hmax-2)) ! if using ZGEEV
      allocate(eigvalHreal(Hmax))
C
CC      Do we need to copy H2p before passing it????

      lwork = -1

      call ZHEEV('V','U',Hmax,H2p,Hmax,eigvalHreal,WORK,LWORK,RWORK,info)

      lwork = int(WORK(1))

      deallocate(work)
      allocate(work(lwork))


      call ZHEEV('V','U',Hmax,H2p,Hmax,eigvalHreal,WORK,LWORK,RWORK,info)
!      call zheevd('v', 'u', hmax, h2p, h2p, eigvalHreal, WORK, LWORK, RWORK, LRWORK,IWORK, LIWORK, INFO)
!      eigvecH = H2p
!      call zcopy(Hmax*Hmax,H2p,1,eigvecH,1)
      eigvalH = eigvalHreal
CC       eigvalH = eigvalH_real! + zimag * eigvalH_imag ! need, as complex?
CC
CC      For dgeev: if eigenvalue is complex, its neighbour will be
CC      conjugate.  The eigenvectors(columns) are expressed differently
CC      (see description in the lapack routine).  If the jth eigenvalue
CC      is complex,j+1th its conjugate; the jth column contains the real
CC      part of the eigenVECTOR for both j and j+1. The j+1th column
CC      contains the imaginary part, which is added to jth eigvec and
CC      subtracted from the j+1th eigvec.

      deallocate(WORK)
      deallocate(RWORK,eigvalHreal)

      end
