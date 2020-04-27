program broad_sig
!==============================================================
! This program manipulates the self-energy file
!  produced by CTQMC (Sig.out).
!
! First it does an average between channels of Sig.out
! The input average array has a specific form,
!  for instance " 1 1 2 3 2 " = electron d system
!  means that first and second channels are averaged (both for
!  spin up and spin down), the third and the fifth also are
!  averaged, while the fourth is independent.
! Any number between 1 and 9 can be chosen to specify the grouped channels.
!
! Then it performs a broadening of the averaged self-energy.
! Every channels is convoluted with a gaussin distribution
!  whose width is frequency-dependent
!     W(x) = Wmax*[ erf1(x) - erf2(x) ]                     ! E or e mode
!     W(x) = Wmax*[ erf2(x) - f(x) ] wiht f(x)=mlin*x+qlin  ! L or l mode
! In the E mode, the parameters of the two erf functions
!  are computed on the basis of 'nom' and 'wmax'.
! The L mode requires two additional parameters startfrq, stopfrq
!
! Optionally it can also set Im[Sigma]=0 whenever Im[Sigma]>0.
!==============================================================
!
! An example of input from the command line is
! echo 'Sig.out 250 1 l " 35  20 200" k "1 1 2 1 2" ' | broad_sig.out
! ! echo 'Sig.out nom nsp l " Wmax  startfrq stopfrq" k "1 1 2 1 2" ' | broad_sig.out
!
! An example that doesn't make any change to Sig.out is
! echo 'Sig.out 250 1 n n n n ' | broad_sig.out
!==============================================================
  implicit none
  character(len=128)      :: fname ,aver_char, brdpar_char
  character(len=1028)     :: line
  integer                 :: nomg,nchan,nhighfrq,header,npar,nom,nsp
  character(len=1)        :: killim,brdmode
  logical                 :: lkillim
  real(8), allocatable    :: omg(:), siggauss(:)
  complex(8), allocatable :: sig(:,:), sigbrd(:,:), av_sig(:,:)
  integer, allocatable    :: brdpar(:)
  integer                 :: io,il,ich,ounit=90


  write(*,'("")')
  write(*,'("+++ Program broad_sig.f90 (6 Nov 2015) +++" )')
  write(*,'("+")')
  ! read input and initialize parameters
  call read_input(ounit,brdmode, nom, nsp, brdpar_char ,fname, killim,aver_char,lkillim,nhighfrq, npar)
  allocate(brdpar(npar))
  call set_parameters(npar,brdmode, brdpar_char , brdpar)


  ! read the file Sig.out
  call count_freq(fname,nomg,nchan,header)
  allocate(omg(nomg))          ; omg(:)      = 0.0
  allocate(siggauss(nomg))     ; siggauss(:) = 0.0
  allocate(sig(nomg,nchan))    ; sig(:,:)    = cmplx(0d0,0d0)
  allocate(av_sig(nomg,nchan)) ; av_sig(:,:) = cmplx(0d0,0d0)
  allocate(sigbrd(nomg,nchan)) ; sigbrd(:,:) = cmplx(0d0,0d0)
  call read_sigfile(fname,nomg,nchan,omg,sig)


  !Compute the average
  call average_sigma(nchan,nomg,nsp,aver_char,sig,av_sig)
  if(trim(aver_char)/='n' .or. trim(aver_char)/='N') fname=trim(fname)//'.avg'
  ! print averaged Sigma
! if(.true.) then
!  write(*,'("  >  ",a)',advance='no') trim(fname)
!  open(60,file=trim(fname))
!  do io=1,nomg
!   write(60,'(f10.5,10(3x,2f12.5))') omg(io),(av_sig(io,ich),ich=1,nchan)
!  enddo
!  close(60)
! endif
  write(*,'("")') ; write(*,'("+")')


  ! constructs the energy-dependent width of the gaussian broadening
  write(*,'("+ Constructing gaussian width, method ",a)',advance='no') trim(brdmode)
  if(brdmode=='E'.or.brdmode=='e') siggauss= construct_width_erf(nomg,omg,nom,brdpar(1))  ! erf1 and erf2
  if(brdmode=='L'.or.brdmode=='l') siggauss= construct_width_lin(nomg,omg,nom,brdpar)     ! linear and erf2
  ! print frequency-dependent width
  if(.true.) then
   write(*,'("  >  width.dat")',advance='no')
   open(unit=60,file='width.dat')
   write(60,'(a)') '# Gaussian width produced in mode '//trim(brdmode)//' with parameters = '//trim(brdpar_char)
   do io=1,nomg
    write(60,'(2f10.5)') omg(io), siggauss(io)
   enddo
  endif
  write(*,'("")')


  ! apply gaussian broadening
  if(brdmode/='N' .or. brdmode/='n') then
    call apply_gauss_brd(nomg,nchan,omg,siggauss,av_sig,sigbrd)
    fname=trim(fname)//'.brd'
  endif


  !kill positive Im[Sigma]
  if(lkillim) call kill_positive_imaginary(nomg,nchan,sigbrd)


  ! print on output file
  write(*,'("+ WRITING OUTPUT SIGMA  >  ",a)') trim(fname)
  do io=1,nomg
   write(ounit,'(f14.8)',advance='no') omg(io)
   do ich=1,nchan
    write(ounit,'(2x,2(x,f14.8))',advance='no') sigbrd(io,ich)
   enddo
   write(ounit,'("")')
  enddo
  close(ounit)

  write(*,'("+")')
  write(*,'("++++++++++++++++++++++++++++++++++++++++++" )')
  write(*,'("")')
  deallocate(omg,sig,sigbrd,av_sig,siggauss)
contains



!===========================================================================
!==============   SUBROUTINES AND FUNCTIONS   ==============================
!===========================================================================

  subroutine read_input(ounit,broadmode, nom, nsp, broadparams, fname, killim,aver_char,lkillim, nhighfrq, nparams)
   implicit none
   integer,             intent(in)  :: ounit
   character(len=128) , intent(out) :: fname ,aver_char,broadparams
   integer            , intent(out) :: nhighfrq,nparams,nom,nsp
   character(len=1)   , intent(out) :: killim,broadmode
   logical            , intent(out) :: lkillim


   read(*,*) fname, nom, nsp, broadmode, broadparams, killim, aver_char

   ! logical to kill positive imaginary part
   if (killim=='k') then
    lkillim=.true.
   else
    lkillim=.false.
   endif

   ! number of frequencies in the tail
   nhighfrq=nomg-nom

   ! number of broading parameters
   if(broadmode=='n' .or. broadmode=='N') nparams=1
   if(broadmode=='e' .or. broadmode=='E') nparams=1
   if(broadmode=='l' .or. broadmode=='L') nparams=3

   open(ounit,file=trim(fname)//'.brd')
   write(ounit,'("# broad_sig.f90 : input parameters = ",a,2x,i3,2x,a1,2x,"(",a,")",2x,a1,2x,"(",a,")")')&
    & trim(fname), nom, broadmode, trim(broadparams), killim, trim(aver_char)
   write(*,'("+ kill poitive imaginary part = ",l2)') lkillim

  end subroutine


!================================================================================


  subroutine set_parameters(npar,broadmode, broadparams, params)
   implicit none
   integer            , intent(in)  :: npar
   character(len=1)   , intent(in)  :: broadmode
   character(len=128) , intent(in)  :: broadparams
   integer            , intent(out) :: params(npar)
   integer :: ip


   ! erf1 and erf2 mode
   if(broadmode=='e' .or. broadmode=='E') then
    read(broadparams,'(1i4)') params  ! wmax
    if(params(1)==0) params(1)=int(nom/20.0)+1

   ! linear and erf2 mode
   else if(broadmode=='l' .or. broadmode=='L') then
    read(broadparams,'(3i4)') params  ! wmax , startfrq, endfrq
    if(params(1)==0) params(1)=int(nom/20.0)+1

   ! no broeading
   else if(broadmode=='n' .or. broadmode=='N') then
    params(1)=0
   else
    write(*,'("++ Broading mode not supported ++")')
    write(*,'("++     PROGRAM TERMINATED      ++")')
    stop
   endif

   write(*,'("+ broading mode and parameters = ",a2)',advance='no') broadmode
   do ip=1,npar
    write(*,'(1x,i3)',advance='no') params(ip)
   enddo
   write(*,'("")') ; write(*,'("+")')
  end subroutine set_parameters


!================================================================================


  subroutine count_freq(fn,numfreq,nchan,headl)
    implicit none
    character(len=128),   intent(in)    :: fn
    integer,              intent(out)   :: numfreq,nchan,headl
    integer :: uin=20
    integer :: nl
    character(len=1028) :: line
    nl=0
    headl=0
    open(uin,file=trim(fn))
    do
     read(uin,'(A)',end=10) line
     if (line(1:1)=='#') then
      headl=headl+1
     else
      nchan=(count_words(line)-1)/2
     endif
     nl=nl+1
    enddo
10  close(uin)
    numfreq=nl-headl
  end subroutine count_freq



!======================================================================



  function count_words(string)
   implicit none
   character(len=1028) , intent(in)  :: string
   integer                        :: count_words
   integer :: pos, i
   pos = 1
   count_words = 0
   do
     i = verify(string(pos:), ' ')
     if (i == 0) exit
     count_words = count_words + 1
     pos = pos + i - 1
     i = scan(string(pos:), ' ')
     if (i == 0) exit
     pos = pos + i - 1
   end do
   return
   end function  count_words



! ==========================================================================



  subroutine read_sigfile(fn,numfreq,nchan,omega,sigma)
    implicit none
    character(len=128), intent(in)    :: fn
    integer,            intent(in)    :: numfreq,nchan
    real(8),            intent(inout) :: omega(numfreq)
    complex(8),         intent(inout) :: sigma(numfreq,nchan)
    integer             :: uin=20
    integer             :: lines2skip=1 , il, ich, io
    real                :: rdsig(numfreq,nchan*2)

    write(*,'("+ Reading file ",a)')trim(fn)
    open(unit=uin,file=trim(fn))
    do il=1,lines2skip
     read(uin,*)
    enddo
    do io=1,numfreq
     read(uin,*) omega(io) , rdsig(io,:)
     do ich=1,nchan
      sigma(io,ich)=cmplx( rdsig(io,ich*2-1) , rdsig(io,ich*2) )
     enddo
    enddo
    close(uin)
    write(*,'("+")')
  end subroutine read_sigfile



! ==========================================================================



  subroutine average_sigma(nchan,nomg,nsp,aver_char,sig,av_sig)
   implicit none
   integer,           intent(in) :: nomg,nchan
   character(len=128),intent(in) :: aver_char
   complex(8),        intent(in) :: sig(nomg,nchan)
   complex(8),        intent(out):: av_sig(nomg,nchan)
   !internal
   integer               :: navch, nsp, ich,jch,group_index,group_number,ig
   character(len=100)    :: formt1
   integer, allocatable  :: aver_rd(:) ! input vector
   integer, allocatable  :: aver(:)    ! averaging groups
   integer, allocatable  :: aver_group(:,:)
   complex(8)            :: av_sig_tmp



   if(trim(aver_char)=='n' .or. trim(aver_char)=='N') then ! it corresponds to av_char= " 1 2 3 4 5 "
     write(*,'("+ No average performed ")',advance='no')
     av_sig=sig
   else
       ! stores av_char into a vector
       navch=nchan/nsp
       allocate(aver_rd(navch))          ; aver_rd=0
       allocate(aver(nchan))             ; aver=0
       write(formt1,'("(",i1,"i2)")') navch
       read(aver_char,trim(formt1))   aver_rd

       ! create averaging groups
       if(nsp==1) then
        aver(1:navch)         = aver_rd
       else
        aver(1:navch)         = aver_rd
        aver(navch+1:2*navch) = aver_rd*10
       endif
       write(*,'("+ averaging channels = ")',advance='no')
       do ig=1,nchan
        write(*,'(i3)',advance='no') aver(ig)
       enddo


       ! do the average
       do io=1,nomg
        do ich=1,nchan
         group_number=0
         group_index=aver(ich)
         av_sig_tmp=cmplx(0.0,0.0)
         !write(*,'("for channel ",i3," I take channels ")',advance='no') ich
         do jch=1,nchan
           if (aver(jch)==group_index) then
            group_number=group_number+1
            av_sig_tmp=av_sig_tmp+sig(io,jch)
            !write(*,'(i4," , ")',advance='no') jch
           endif
         enddo
         !write(*,'("")')
         av_sig(io,ich) = av_sig_tmp/group_number
        enddo
       enddo
   endif
  end subroutine average_sigma



! ==========================================================================


  subroutine apply_gauss_brd(nomg,nch,omg,wg,sig,sigbrd)
    ! simple convolution integral
    implicit none
    integer,    intent(in) :: nomg,nch
    real(8),    intent(in) :: omg(nomg), wg(nomg)
    complex(8), intent(in) :: sig(nomg,nch)
    complex(8), intent(out):: sigbrd(nomg,nch)
    integer :: io,jo,ich
    real(8) :: norm,gauss
    complex(8) :: convol(nch)

    do io=1,nomg
     norm=0.0 ; convol(:)=0.0
     do jo=1,nomg
      if (abs(wg(io)) < 0.2 ) then
       ! if wdith==0 , then gauss=delta function
       if( abs(io-jo)<1.0e-3 ) then
        gauss=1.0
       else
        gauss=0.0
       endif
      else
       ! otherwise it's a proper gaussian distriution
       gauss=exp(-0.5*( (io-jo)/wg(io))**2 )
      endif
      norm=norm+gauss
      do ich=1,nch
       convol(ich)=convol(ich)+gauss*sig(jo,ich)
      enddo
     enddo
     sigbrd(io,:)=convol(:)/norm
    enddo
  end subroutine apply_gauss_brd

! ==========================================================================



  function sigma_gauss(nomg,omg,cntr_omg,w0,wk)
   implicit none
   integer, intent(in)     :: nomg
   real(8), intent(in)     :: w0,wk,omg(nomg),cntr_omg
   real(8)                 :: gaussian_weight,norm,maxw
   integer                 :: io
   real(8),dimension(nomg) :: sigma_gauss

   maxw=80
   do io=1,nomg
    gaussian_weight= maxw*exp(-0.5*( (omg(io)-cntr_omg)/wk )**2 )
    sigma_gauss(io) = w0 + gaussian_weight
   enddo
   return
  end function



! ==========================================================================



  function construct_width_erf(nomg,omg,nom,wmax)
   implicit none
   integer, intent(in)      :: nomg,nom,wmax
   real(8), intent(in)      :: omg(nomg)
   integer                  :: io
   ! parameters and arguments of erf1, erf2 and final width
   real(8)                  :: dilat1,dilat2,x
   integer                  :: wshift1, wshift2
   ! final width
   real(8)                  :: w
   real(8), dimension(nomg) :: construct_width_erf

   ! First compute the parameteres of the two erf functions
   !wmax=!int(nom/20d0)+1
   wshift1 = int(2.5*wmax)+1
   wshift2 = int(nom+(nomg-nom)/3d0)+1
   dilat1  = 0.3333*nom
   dilat2  = 1d0

   ! Then compute the value of the width
   do io=1,nomg
    ! record first erf1
    x=( real(io-wshift1)/real(nom) )*dilat1   ! scaled and shifted argument
    w=0.5*erf(x)
    ! subtract second erf
    x=( real(io-wshift2)/real(nom) )*dilat2   ! scaled and shifted argument
    w=w-0.5*erf(x)
    ! shift , scale and adjust the result
    w=w*wmax
!    if (w<=1.0) then
!     w=0.0
!    endif
    construct_width_erf(io) = w
   enddo
   return
  end function



! =============================================================================



  function construct_width_lin(nomg,omg,nom,params)
   implicit none
   integer, intent(in)      :: nomg,nom,params(3)
   real(8), intent(in)      :: omg(nomg)
   integer                  :: io
   ! parameters and arguments of erf2 part
   real(8)                  :: dilat2
   integer                  :: wshift2
   ! parameters and arguments of the liner part
   integer                  :: startlin , stoplin, wmax
   real(8)                  :: mlin, qlin
   ! final width
   real(8)                  :: w,x
   real(8), dimension(nomg) :: construct_width_lin

   ! broading parameteres
   wmax = params(1) ; startlin = params(2) ; stoplin = params(3)
   mlin = -dble(wmax)/dble(stoplin-startlin) ; qlin = -dble(mlin*stoplin)
   wshift2 = int(nom+(nomg-nom)/3d0)+1 ; dilat2  = 1d0

   ! construct the frequency-dependent gaussian width
   do io=1,nomg
     ! construct erf on tails
     x=( real(io-wshift2)/real(nom) )*dilat2   ! scaled and shifted argument
     w=-0.5*wmax*(erf(x)-1d0)

     ! construct low-energy (linear) part
     if(io<startlin) then
       w=w-wmax
     else if(io>=startlin .and. io<stoplin) then
       w=w-(mlin*io+qlin)
     endif

     ! record width
     construct_width_lin(io) = w
   enddo
   return
  end function



! =============================================================================



  subroutine kill_positive_imaginary(no,nc,s)
   implicit none
   integer,                         intent(in) :: no,nc
   complex(8), dimension(no,nc), intent(inout) :: s
   integer :: io,ic,nkill=0

   do io=1,no
    do ic=1,nc
     if(aimag(s(io,ic)) > 0.0 ) then
      s(io,ic)=cmplx(dble(s(io,ic)),0.0)
      nkill=nkill+1
     endif
    enddo
   enddo

   write(*,'("+")')
   write(*,'("+ N. of positive Im[Sig.brd] = ",i3)')nkill
   write(*,'("+")')
  end subroutine


!
! subroutine fast_gauss_brd(sig,nomg,nch,dos,omg,dosbrd)
!   ! like apply_gauss_brd, but the internal loop (wj) is
!   !  performed in a smaller interval, so to speed up the convolution integral.
!   implicit none
!   integer, intent(in) :: nomg,nch
!   real(8), intent(in) :: sig,dos(nomg,nch), omg(nomg)
!   real(8), intent(out):: dosbrd(nomg,nch)
!   integer :: wi,wj,ich
!   integer :: wjmin,wjmax,fivesigma
!   real(8) :: norm,convol(nch),gauss,dom

!   ! compute fivesigma in indexes of omega
!   !  ASSUMING CONSTANT STEP ON FREQUENCY GRID
!   dom=omg(2)-omg(1)
!   fivesigma=int(5.0*sigma/dom)+1

!   ! Make fast convolution
!   do wi=1,nomg
!    norm=0.0 ; convol(:)=0.0
!    wjmin=max(1,wi-fivesigma)
!    wjmax=min(nomg,wi+fivesigma)
!    do wj=wjmin,wjmax
!     gauss=exp(-0.5*( (omg(wj)-omg(wi))/sig )**2 )
!     norm=norm+gauss
!     do ich=1,nch
!      convol(ich)=convol(ich)+gauss*dos(wj,ich)
!     enddo
!    enddo
!    dosbrd(wi,:)=convol(:)/norm
!   enddo

! end subroutine fast_gauss_brd
!





  subroutine write_sigbrd(fn,nch,nomg,omg,sigbrd)
    implicit none
    character(len=7),   intent(in)  :: fn
    integer,            intent(in)  :: nch,nomg
    real(8),            intent(in)  :: omg(nomg)
    complex(8),         intent(in)  :: sigbrd(nomg,nch)
    character(len=11)  :: fnout
    character(len=256) :: head1
    integer  :: uout=21, uin=20
    integer  :: io,ich

    ! copy header
    open(unit=uin,file=fn)
    read(uin,'(a)') head1
    close(uin)

   ! write output file
    fnout=fn//'.brd'
    open(unit=uout,file=fnout)
    write(uout,'(a)') trim(head1)
    do io=1,nomg
     write(uout,'(f10.5,8(3x,2f12.5))') omg(io),(sigbrd(io,ich),ich=1,nch)
    enddo
    close(uout)

  end subroutine write_sigbrd




end program broad_sig
