      subroutine makedelta3(s_dmft,ldcix,nsp,ncix,nomg,gloc,eimpm,omega,omfac,delta,deltah5,g0h5,sigdc)
C- Make DMFT hybridization function delta
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci     Elts read:  l icix
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:icix sig
Cio    Passed to:  ineqcix sigp2sigij
Ci Inputs
Cl   ldcix :dimension of largest cix block, for dimensioning arrays
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ncix  :Number of correlated blocks
Ci   nomg  :dimensioning parameter, number of frequency points
Ci   gloc  :local interacting Green's function
Ci   eimpm :local impurity level
Ci   omega :frequency in Ry
Ci   omfac :factor accounting for real or Matsubara frequency (either i or 1)
Co Outputs
Co   delta  :hybridization function, stored in compressed format
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   18 May 18 (MvS) Enable equivalent cix blocks with spin flip
Cu   30 Apr 18 (MvS) First cut at equivalent cix blocks
Cu   27 Feb 16 (MvS) Adapted form makedelta2: redesign for more general cix blocks
Cu   13 Nov 15 Sponza changed to handle new definition of Eimp
Cu   22 Aug 15 (Pisanti) spin polarized
Cu   23 Apr 15  P Pisanti added imp and dc sigma
Cu   13 Oct 14  P Pisanti first created
C ----------------------------------------------------------------------
      use structures
      use h5
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft) ::  s_dmft
      type(h5space)  ::  fs_h5,ms_h5

C ... Passed parameters
      integer, intent(in) :: ldcix,nsp,nomg,ncix
      real(8), intent(in) :: omega(nomg)
      complex(8), intent(in)  :: omfac
      complex(8), intent(in)  :: gloc(ldcix,ldcix,nsp,ncix,nomg)
      complex(8), intent(in)  :: eimpm(ldcix,ldcix,nsp,*)
      real(8),    intent(out) :: delta(nomg,*)
      complex(8), intent(out) :: deltah5(nomg,*)
      complex(8), intent(out) :: g0h5(nomg,*)
      logical :: debug=.false.
C ... Dynamically allocated local arrays
C ... Local parameters
      real(8), parameter :: ry2eV = 13.60569193d0
      integer :: isp,ksp,iomg,cix,icix,cixi,m1,m2,nm,j,nicix(ncix)
      complex(8) :: glocinv(ldcix,ldcix,nsp),deltai(ldcix,ldcix,nsp),deltas(ldcix,ldcix,nsp),g0(ldcix,ldcix,nsp)
      complex(8) :: sigi(ldcix,ldcix,nsp),DC(ldcix,ldcix,nsp)
      real(8) :: sigdc(s_dmft%ndsig)
      integer ::il1,il2
      call fs_h5 % init( [ldcix,ldcix,nsp,nomg,ncix])
      call ms_h5 % init( [ldcix,ldcix,nsp])

C ... Count the number of cix per inequivalent cix
      call ineqcix(s_dmft,ncix,nicix)

C --- Make delta for each inequivalent cix ---
      do  iomg = 1, nomg
         do  cixi = 1, ncix     ! Loop all inequivalent cixi.  Inner loop will group equivalent ones
            if (nicix(cixi) >= 0) cycle ! nicix(:)<0 => first of equivalent cix
            icix = iabs(s_dmft%icix(cixi))
            call dpzero(deltas,2*size(deltas))
            do  cix = 1, ncix   ! Loop over cix equivalent to cixi: gather average
               if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
               nm = 2*s_dmft%l(icix)+1
C     if (nm > ldcix) call rx('bug in sudmft')

C     ... Gloc^-1
               do  isp = 1, nsp
                  do  m1 = 1, nm
                     do  m2 = 1, nm
                        glocinv(m1,m2,isp) = gloc(m1,m2,isp,cix,iomg)
                     enddo
                  enddo
                  call zinv(nm,glocinv(1,1,isp),ldcix)
C     call zprm('ginv',2,glocinv(1,1,isp),ldcix,nm,nm)
               enddo            ! spin

C     ... Decompress sigma for this cix block and omega, orbital repsn
               call sigp2sigij(s_dmft,1000,nomg,nsp,iomg,cix,ldcix,s_dmft%sig,s_dmft%sig,sigi)

C     ... Hybridization function: Delta = omega - Gloc^-1 - E_{imp} - Sigma
               deltai(:,:,:) =  - glocinv(:,:,:) - eimpm(:,:,:,cix) - sigi(:,:,:)
               g0(:,:,:) = glocinv(:,:,:) + sigi(:,:,:)
               do isp = 1 , nsp
                  do m1 = 1, nm
                     do m2 = 1, nm
                        if (s_dmft%sigind(m1,m2,icix,isp) >0 ) then
                           g0(m1,m2,isp)  = g0(m1,m2,isp)    + sigdc(s_dmft%sigind(m1,m2,icix,isp))
                        endif
                     enddo
                  enddo
                  call zinv(nm,g0(1,1,isp),ldcix)
               enddo
               call zdscal(size(g0), 1./ry2eV, g0, 1)
               call fs_h5 % select_hyperslab(start=[0,0,0,iomg-1,cixi-1], count=[ldcix,ldcix,nsp,1,1])
               call h5_write('g0.h5:/g0', g0(1,1,1), find_h5dtype(g0), file_space = fs_h5, mem_space = ms_h5)
               call zdscal(size(g0), ry2eV, g0, 1)

               do  m1 = 1, nm
                  deltai(m1,m1,1:nsp) = deltai(m1,m1,1:nsp) + omega(iomg)*omfac
               enddo
C     do  isp = 1, nsp
C     call yprmi('deltai, cix=%i icix=%i',cix,icix,3,deltai(1,1,isp),0,ldcix,ldcix,ldcix)
C     enddo
               do  isp = 1, nsp

                  if(s_dmft% ig(cix) /= 1) then

                  endif
               enddo


               j = ldcix**2
               do  isp = 1, nsp
                  ksp = isp
                  if (s_dmft%icix(cix) < 0 .and. nsp == 2) ksp = nsp+1-isp ! spin flip for AFM symmetry

                  call daxpy(2*j,1/dble(iabs(nicix(cixi))),deltai(1,1,ksp),1,deltas(1,1,isp),1)

               enddo

            enddo               ! loop over all blocks with equivalent cix

C     ... Store delta in compressed form
            call sigp2sigij(s_dmft,2,nomg,nsp,iomg,cixi,ldcix,delta,delta,deltas)
            call sigp2sigij(s_dmft,12,nomg,nsp,iomg,cixi,ldcix,delta,deltah5,deltas)
            call sigp2sigij(s_dmft,12,nomg,nsp,iomg,cixi,ldcix,g0,g0h5,g0)
C     call prmx('delta',delta(iomg,1),nomg,1,ldcix*nsp*2)

         enddo                  ! cix loop
      enddo                     ! omega loop

      end subroutine makedelta3

      subroutine ineqcix(s_dmft,ncix,nicix)
C     - Returns nicix with information about which cix are equivalent
C     ----------------------------------------------------------------------
C     i Inputs
C     io Structures
C     io  s_dmft
C     i     Elts read:  icix
C     o     Stored:     *
C     o     Allocated:  *
C     io    Elts passed:*
C     io    Passed to:  *
C     i   ncix
C     o Outputs
C     o   nicix :nicix(cix) <0 cix is the first occurence of this inequivalent cix
C     o                        |nicix(cix)| is the number of occurences
C     o                      0 => no occurence of this cix
C     o                     >0 => points to first inequivalent occurence of this cix
C     o         :nicix(icix,2) = number of inequivalent cix
C     s Command-line switches
C     l Local variables
C     l         :
C     r Remarks
C     r  To loop over all equivalent cix group together use a looping contstruct like this:
C     r    do  cixi = 1, ncix        ! Loop over all inequivalent cixi
C     r      if (nicix(cixi) >= 0) cycle ! nicix(:)<0 => first of equivalent cix
C     r      icix = iabs(s_dmft%icix(cixi)) ! inequivalent cix index
C     r      do  cix = 1, ncix       ! Loop over cix equivalent to cixi
C     r        if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
C     r        print *, 'body of calculation here',cix,cixi
C     r      enddo
C     r    enddo
C     u Updates
C     u   27 Apr 18 First created
C     ----------------------------------------------------------------------
      use structures
      implicit none
C     ... For structures
!     include 'structures.h'
      type(str_dmft)::  s_dmft
C     ... Passed parameters
      integer, intent(in)  :: ncix
      integer, intent(out)  :: nicix(ncix)
C     ... Local parameters
      integer cix,icix,cix2

      call iinit(nicix,size(nicix))
      do  cix = 1, ncix
         if (nicix(cix) /= 0) cycle ! Already marked; skip
         icix = iabs(s_dmft%icix(cix))
         nicix(cix) = -1        ! Initial size of this group
         do  cix2 = cix+1, ncix
            if (iabs(s_dmft%icix(cix2)) /= icix) cycle ! block not equivalent
            nicix(cix2) = cix   ! Point to first member
            nicix(cix) = nicix(cix) - 1 ! Increment the number of members
         enddo
      enddo

      end

C#ifdefC TEST
C      subroutine fmain
C
C      use structures
C      implicit none
CC ... For structures
C!      include 'structures.h'
C      type(str_dmft)::  s_dmft
C      integer ncix,nicix(5),cixi,cix,icix
C
C      ncix = 5
C      allocate(s_dmft%icix(ncix))
C
C      s_dmft%icix(1) = 1
C      s_dmft%icix(2) = 3
C      s_dmft%icix(3) = 1
C      s_dmft%icix(4) = 1
C      s_dmft%icix(5) = 3
C      print *, s_dmft%icix
C
C      call ineqcix(s_dmft,ncix,nicix)
C
CCr  To loop over all equivalent cix group together use a looping contstruct like this:
C      do  cixi = 1, ncix        ! Loop over all inequivalent cixi
C        if (nicix(cixi) >= 0) cycle ! nicix(:)<0 => first of equivalent cix
C        icix = iabs(s_dmft%icix(cixi)) ! inequivalent cix index
C        do  cix = 1, ncix       ! Loop over cix equivalent to cixi
C          if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
C          print *, 'body of calculation here',cix,cixi
C        enddo
C      enddo
C
C      end
C#endif
