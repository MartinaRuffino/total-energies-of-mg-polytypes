C#define OPENMP
      subroutine mkzxq_ml(nqbz,nbnbx,npm,k,nband,nt0,ntp0,nctot,ncc,nbmx,
     .  nbcut,nbcut2,ngb,ngbb,itps,nkmax,nkqmax,nbnb,n1b,n2b,ihw,nhw,
     .  jhw,whw,nwt,fack,zmelt,rcxq,nmbas1,nmbas2)

C- Imaginary part of bare response function
C ----------------------------------------------------------------------
Ci   nqbz  :
Ci   nbnbx :
Ci   npm   :1 if time reversal symmetry, 2 if not
Ci   k     :index to current kpoint
Ci   nband :number of energy bands
Ci   nt0   :number of eigenvectors associated w/ n1b;
Ci         :used here for dimensioning only
Ci   ntp0  :number of eigenvectors associated w/ n12b;
Ci         :used here for dimensioning only
Ci   nctot :number of core states for n
Ci   ncc   :number of core states for n'
Ci   nbmx  :max number of bands to include in calculation
Ci         :(value of nband_chi0 in GWinput)
Ci   nbcut :(value of nbcutlow in GWinput)
Ci   nbcut2:(value of nbcutlowto in GWinput)
Ci   ngb   : total number of product functions
Ci   ngbb  :dimensions rcxq for now.  Normally ngb, but
C          :in future, contracted size of product basis.
Ci   itps  :nkqmin, where nkqmin = min(n2b(ibib,k,jpm), or 1
Ci   nkmax :max value of n1b
Ci   nkqmax:max value of n2b
Co   nbnb  :nbnb(k) = number of nonzero (ib,jb) pairs for kth qp
Co   n1b   :n1b(ix,k)= ib for ixth nonzero (ib,jb) pair
Co   n2b   :n2b(ix,k)= jb for ixth nonzero (ib,jb) pair
Ci   ihw   :ihw(nbnb(k),k) = index to first histogram bin within
Ci                           (demin,demax) for a given nbnb pair and k
Ci   nhw   :nhw(nbnb(k),k) = number of histograms bins witin
Ci                           (demin,demax) for a given nbnb pair and k
Ci   jhw   :jhw(nbnb(k),k) = pointer to whw; see whw
Ci   whw   :whw(i:j) histogram weights in bins i:j for given ib,jb,kx
Ci         :i = jhw(nbnb,kx)
Ci         :j = jhw(nbnb,kx)+nhw(nbnb),kx)-1
Ci         :nbnb = nbnb(ib,jb,kx)
Ci   nwt   :number of histogram points (aka as nwhis)
Ci   fack  :BZ weighting from this k-point
Ci   zmelt :<psi(k+q,t') | psi(k,t) M(R,i)>
Co Outputs
Co   rcxq  :noninteracting response function on the positive real axis zxq
Co         :A portion of rxcq is accumulated.
Co         :Not summed over spin here (done by hx0fp0)
Cr Remarks
Co   zxq(i,j,iw) =  (-1/pi) \int_{freqr(iw)}^{freqr(iw+1)} d\omega
Co                           Im <M(i)|x0(q,\omega)|M(j)>
Co   Evaluate and store in rcxq as follows:
Co    z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
Co    rcxq(iw,ibg1,ig2) = sum_ibib wwk(iw,ibib)* z1p(ibib, ig1,ig2)
Cu Updates
Cu   15 May 13 Adapted from x0kf_v4h.F
C ----------------------------------------------------------------------

      use mpi_mod

      implicit none

C ... Passed parameters

      integer nqbz,npm,nbnbx,nband,nt0,ntp0,nctot,ncc,nbmx,nbcut,nbcut2,
     .  ngb,ngbb,itps,nkmax,nkqmax,nwt,k
      integer :: nmbas1, nmbas2
      integer :: nbnb(nqbz,npm),n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm),
     .  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8) :: whw(*),fack
      complex(8):: zmelt(ngb,nctot+nt0,ncc+ntp0)
!ML old      complex(8):: rcxq(ngbb,ngbb,nwt,npm)
      complex(8):: rcxq(nmbas1,nmbas2,nwt,npm)

C ... Local parameters

      integer ibib,jpm,it,itp,ig1,ig2,iw,ilo,ihi
      logical :: iww1=.true.,iww2=.true.
      complex(8) :: zmelt2
      complex(8),allocatable::  z1p(:,:)
      real(8),allocatable::  wgti(:)
      real(8) :: imagweight
      real(8):: delwc

C ... for OPENMP
C#ifdef OPENMP
      integer nthreads,omppid,tid,iwg
      integer,allocatable:: firstg(:),lastg(:)
C#ifdefC OPENMP2
C      integer,allocatable,target:: first0(:),last0(:),idx0(:,:)
C      integer,allocatable,target:: first1(:),last1(:),idx1(:,:)
C      integer,allocatable,target:: first2(:),last2(:),idx2(:,:)
C      integer,allocatable,target:: first3(:),last3(:),idx3(:,:)
C      integer,pointer :: first(:),last(:),idx(:,:)
C#endif

C ... OMP setup
C$OMP parallel shared(nthreads)
      nthreads = omppid(0)
C$OMP END PARALLEL
      nthreads = max(nthreads,1)

      allocate(firstg(0:nthreads-1),lastg(0:nthreads-1))



!ML old      call thdist1(1,1d0,2d0,nthreads,1,ngb,firstg,lastg)
      call thdist1(1,1d0,2d0,nthreads,1,nmbas2,firstg,lastg)

!ML:  this distributes the loop over the second index igb2, which runs to nmbas2

C This version merges iw and ig2 loops
C But doesn't help: seems to be limited by memory access to shared rxcq
C#ifdefC OPENMP2
C      allocate(first0(0:nthreads-1),last0(0:nthreads-1),idx0(2,1*ngb))
C      allocate(first1(0:nthreads-1),last1(0:nthreads-1),idx1(2,2*ngb))
C      allocate(first2(0:nthreads-1),last2(0:nthreads-1),idx2(2,3*ngb))
C      allocate(first3(0:nthreads-1),last3(0:nthreads-1),idx3(2,4*ngb))
C      call thdist2(1,1d0,2d0,nthreads,0,0,1,ngb,first0,last0,idx0)
C      call thdist2(1,1d0,2d0,nthreads,0,1,1,ngb,first1,last1,idx1)
C      call thdist2(1,1d0,2d0,nthreads,0,2,1,ngb,first2,last2,idx2)
C      call thdist2(1,1d0,2d0,nthreads,0,3,1,ngb,first3,last3,idx3)
C#endif

C      print 101, '0l:',last0
C      print 101, '0s:',last0-first0
C      print 101, '1l:',last1
C      print 101, '1s:',last1-first1
C      print 101, '2l:',last2
C      print 101, '2s:',last2-first2
C      print 101, '3l:',last3
C      print 101, '3s:',last3-first3
C  101 format(a,20i4)
C      call rx('done')
C#endif

C     print *, 'mkzxq starting time',nbnb(k,1),ngb,delwc()
      allocate(z1p(ngb,ngb),wgti(nwt))

      z1p = 0.d0
      wgti = 0.d0

C --- Loop over jpm=1,2 if no TR symmetry ---
      do  jpm = 1, npm

C   --- Sum over (n, n') pairs for this k ---
        do  ibib = 1, nbnb(k,jpm) !---  ibib loop


!           write(0,'("[",I4,"] mkzxq: jhw(",I3,",",I3,",",I3,") = ",I5)') mpi_rank,ibib,k,jpm,jhw(ibib,k,jpm)

          if (n1b(ibib,k,jpm) <= nband) then
            it = nctot + n1b(ibib,k,jpm) !valence
            if (it > nctot + nkmax ) cycle
          else
            it = n1b(ibib,k,jpm) - nband !core
          endif
          if (n2b(ibib,k,jpm) <= nband) then
            itp = ncc + n2b(ibib,k,jpm) - itps + 1 !val
            if (itp > ncc + nkqmax-itps+1 ) cycle
          else
            itp = n2b(ibib,k,jpm) - itps + 1 - nband !core
          endif

C     ... Check for truncation of states
C         jpm = 1 : whether or not TR symmetry
          if (jpm==1) then
            if (n2b(ibib,k,jpm)>nbmx)  then !nbmx
              if (iww1) then
                print *,' nband_chi0 nbmx=',nbmx
                iww1 = .false.
              endif
              cycle
            endif
            if (n1b(ibib,k,jpm)<=nbcut.and.nbcut2<n2b(ibib,k,jpm)) then
              if (iww2) then
!                write(*,"('[',I4,'] nband_chi0 nbcut nbcut2 n2b n1b=',4i6)")
!     .            mpi_rank, nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                iww2 = .false.
              endif
              cycle
            endif

C         jpm = 2 : Executed only if no TR symmetry
          else
            if (n1b(ibib,k,jpm) > nbmx) then !nbmx
              if (iww1) then
                print *,' nband_chi0 nbmx=',nbmx
                iww1 = .false.
              endif
              cycle
            endif
            if (n2b(ibib,k,jpm)<=nbcut.and.nbcut2<n1b(ibib,k,jpm)) then
              if (iww2) then
!                write(*,"('[',I4,'] nband_chi0 nbcut nbcut2 n2b n1b=',4i6)")
!     .            mpi_rank, nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                iww2 = .false.
              endif
              cycle
            endif
          endif

!ML  debug

!                write(*,"('[',I4,'] k ibib jpm nbcut n2b n1b=',6i6)")
!     .            mpi_rank, k, ibib, jpm, nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)


!ml

! remember that it and itp are determined from ibib
!
!

C#ifdef OPENMP
C$OMP parallel private(tid,ig2,ig1,zmelt2)
          tid = omppid(1)
          do  ig2 = firstg(tid), lastg(tid)
            zmelt2 = zmelt(ig2,it,itp)
            do  ig1 = 1, min(ig2,nmbas1)
              z1p(ig1,ig2) = dconjg(zmelt(ig1,it,itp)) * zmelt2
            enddo
          enddo
C$OMP end parallel

!          do ig2=1,min(nmbas2,10)
!             write(*,'("[",I4,"] k = ",I3,", ibib = ",I3,", z1p(:,",I3,") = ",5("("E11.2,",",E11.2,") "))')
!     &            mpi_rank, k, ibib, ig2,(z1p(ig1,ig2),ig1=1,min(5,nmbas1))
!          enddo




C#elseC
C          do  ig2 = 1, ngb
C            zmelt2 = zmelt(ig2,it,itp)
C            do  ig1 = 1, ig2
C              z1p(ig1,ig2) = dconjg(zmelt(ig1,it,itp)) * zmelt2
C            enddo
C          enddo
C#endif

          ilo = ihw(ibib,k,jpm)
          ihi = ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
          if (ihi > nwt) call rx('x0kf_v4hz: iw>nwt')
          do  iw = ilo, ihi
            wgti(iw) = fack * whw(jhw(ibib,k,jpm)+iw-ilo)
          enddo

C ... version that merges iw,ig2 loops
C#ifdef OPENMP
C#ifdefC OPENMP2
C          select case (ihi-ilo)
C            case (0); first => first0; last => last0; idx => idx0
C            case (1); first => first1; last => last1; idx => idx1
C            case (2); first => first2; last => last2; idx => idx2
C            case (3); first => first3; last => last3; idx => idx3
C            case default ;
C              allocate(first(0:nthreads-1),last(0:nthreads-1))
C              allocate(idx(2,(ihi-ilo+1)*ngb))
C              call thdist2(1,1d0,0d0,nthreads,ilo-ilo,ihi-ilo,1,ngb,
C     .          first,last,idx)
C          end select
C
CC$OMP parallel private(tid,iwg,iw,ig2,ig1,imagweight) shared(rcxq)
C          tid = omppid(1)
CC         print *, 'tid',tid,first0(tid),last0(tid)
C          do  iwg = first(tid), last(tid) ! compound index
C            iw  = idx(1,iwg)+ilo   ! outer loop index
C            ig2 = idx(2,iwg)  ! inner loop index
C            imagweight = wgti(iw)
CC           print *, 'tid',tid,iw,ig2,imagweight
C            do  ig1 = 1, ig2
C              rcxq(ig1,ig2,iw,jpm) = rcxq(ig1,ig2,iw,jpm) +
C     .          z1p(ig1,ig2)*imagweight
C            enddo               ! ig1
C          enddo                 ! ig2+iw
CC$OMP end parallel
C          if (ihi-ilo > 3) then
C            deallocate(first,last,idx)
C          endif
C#else



C$OMP parallel private(tid,iw,ig2,ig1,imagweight)
          tid = omppid(1)
          do  iw = ilo, ihi
            imagweight = fack * whw(jhw(ibib,k,jpm)+iw-ilo)
            do  ig2 = firstg(tid), lastg(tid)
              do  ig1 = 1, min(ig2,nmbas1)
                rcxq(ig1,ig2,iw,jpm) =
     .            rcxq(ig1,ig2,iw,jpm) + z1p(ig1,ig2)*imagweight
              enddo             ! ig1
            enddo               ! ig2
          enddo                 ! iw
C$OMP end parallel
C#endif
C#elseC
CC generic routine
CC$OMP     parallel do private(imagweight,ig1,ig2,iw) shared(rcxq)
C          do  iw = ilo, ihi
C            imagweight = fack * whw(jhw(ibib,k,jpm)+iw-ilo)
C            do  ig2 = 1, ngb
C              do  ig1 = 1, ig2
C                rcxq(ig1,ig2,iw,jpm) =
C     .            rcxq(ig1,ig2,iw,jpm) + z1p(ig1,ig2)*imagweight
C              enddo             ! ig1
C            enddo               ! ig2
C          enddo                 ! iw
CC$OMP     end parallel do
C
C#endif

        enddo                   ! ibib
      enddo                     ! jpm

!      do ig1=1,ngb
!         write(*,'("[",I4,"] sum abs zmelt( k= ",I4," ) = "d13.6)')
!     &        mpi_rank, k, abs(zmelt(ig1,1,1))
!      enddo


!      do ig1=1,ngb
!         write(*,'("[",I4,"] sum abs z1p( ",I4,I3," ) = "d13.6)')
!     &        mpi_rank, k, ig1,  sum(abs(z1p(ig1,:)))
!      enddo

!      do ig2=1,ngb
!         write(*,'("[",I4,"] sum abs rcxq(",I4,I3," ) = "d13.6)')
!     &        mpi_rank, k, ig2, sum(abs(rcxq(1:ig2,ig2,ilo:ihi,1)))
!      enddo

      deallocate(z1p,wgti)
C#ifdef OPENMP
      deallocate(firstg,lastg)
C#ifdefC OPENMP2
C      deallocate(first0,last0,idx0)
C      deallocate(first1,last1,idx1)
C      deallocate(first2,last2,idx2)
C      deallocate(first3,last3,idx3)
C#endif
C#endif

C      print *, 'mkzxq exit',k,sum(rcxq)
C      print *, 'mkzxq ending time',delwc() ; stop

      end
