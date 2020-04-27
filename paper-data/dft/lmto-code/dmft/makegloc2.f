      subroutine makegloc2(s_dmft,ncix,ldcix,ndhamx,nlo,nhi,ikp,nkp,nsp,nspx,iomg,
     .  ef,broad,evlk,omega,dmftu,sigbark,gpmode,nompr,gloc,gdiag,ifgk,gkloc)
C- Compute the local k-resolved Green's function according to Pisanti's notes
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Ci     Elts read:  lsigim icix l
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lsigim
Cio    Passed to:  *
Ci Inputs
Ci   ncix   :dimensioning parameter, number of correlated l/site channels
Ci   ldcix  :dimensioning parameter, max number of correlated orbitals (lms character)
Ci   ndhamx :dimensioning parameter, largest hamiltonian dimension
Ci   nev0   :offset to first band in evlk to be included in G
Ci   nev    :number of bands to include in g
Ci   nkp    :dimensions gdiag
Ci   nsp    :dimensioning parameter,2 for spin-polarized case, otherwise 1
Ci   iomg   :index to current frequency
Ci   ef     :chemical potential (updated for current iteration)
Ci   broad  :Gaussian broadening, in Ry
Ci   evlk   :eigenvalues from the QSGW/LDA Hamiltonian
Ci   omega  :current frequency in Ry
Ci   dmftu  :projector to local correlated Hilbert space
Ci   sigbark:(Sig.inp - SigDC) in lattice representation {ij} for one k and omega (embed_sigma)
Ci   gpmode :>1 and <8 write to ifgk local Green's function for this k and omega
Ci          :<-1       write to ifgk Hamiltonian for this k, H = omega-1/g
Ci          :8         write to ifgk Green's function for this k and omega (binary)
Ci          :9         write to ifgk Sigma for this k and omega (binary)
Ci          :18        return diagonal part of G, 1-particle basis
Ci          :20        return effective diagonal sigma which generates diagonal g
Ci   nompr  :number of frequencies for which to write gloc(k) (|gpmode|>1)
Ci   ifgk   :file logical unit for writing gloc(k), hloc, or G
Co Outputs
Co   gdiag  :Diagonal part of G, 1-particle basis
Co   gloc   :local interacting Green's function, add contribution from one k
Cr Remarks
Cu Updates
Cu   09 Feb 18 (MvS) new option gpmode=20
Cu   17 Jan 18 (MvS) First cut at SO coupled case
Cu   31 Mar 17 (MvS) new option gpmode=8
Cu   25 Oct 16 (MvS) gpmode<-1 print hamiltonian instead of g
Cu   05 Oct 16 (Sponza)  New gpmode
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu    9 Nov 15 (Sponza ) modified version taking sigbar in input
Cu   22 Aug 15 (Pisanti) spin polarized
Cu   23 Apr 15 (Pisanti) g can include DMFT sigma and DC
Cu   19 Oct 14 (Pisanti) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft)      ::  s_dmft
C ... Passed parameters
      integer,    intent(in)  :: ncix,ldcix,ndhamx,nlo,nhi,ikp,nkp,nsp,nspx
      integer,    intent(in)  :: ifgk,gpmode,nompr,iomg
      real(8),    intent(in)  :: ef,broad
      real(8),    intent(in)  :: evlk(ndhamx,nspx) ! QSGW/LDA eigenvalues
      complex(8), intent(in)  :: omega
!!
      complex(8), intent(in)  :: dmftu(nlo:nhi,ldcix,nsp,ncix)
      complex(8), intent(in)  :: sigbark(nlo:nhi,nlo:nhi,nspx)
      complex(8), intent(out) :: gloc(ldcix,ldcix,nsp,ncix)
      complex(8), intent(out) :: gdiag(nompr,nlo:nhi,nkp,nspx)
      complex(8), intent(out) :: gkloc(ldcix,ldcix,nsp,ncix)
C ... Dynamically allocated local arrays
      complex(8),allocatable :: tmp1(:,:),gij(:,:,:)

C ... Local variables
      real(8)    :: xx
      complex(8) :: hkloc
      real(8), parameter :: ry2eV = 13.60569193d0
      integer    :: cix,icix,i,j,isp,ksp,nl1,nev,nspc,ispc
      character(len=256) :: strn
C ... External calls
      external dpzero,zgemm,zinv
      procedure(logical) :: cmdopt

      nev = nhi-nlo+1
      allocate(tmp1(ldcix,nev),gij(nlo:nhi,nlo:nhi,nspx))
      call dpzero(gij,2*size(gij))
      nspc = 1; if (nsp > nspx) nspc = 2


      call dpzero(gkloc,2*size(gkloc))
      do  isp = 1, nspx

C   --- Full Green's function (eigenfunction basis) and eigenvalues in energy subspace ---
C       One-body part of g is diagonal
C       [g_ij]^-1 = delta_ij [omega + ef - evl] - sigbark_ij
        do  i = nlo, nhi
          do  j = nlo, nhi
            gij(i,j,isp) = -sigbark(i,j,isp)
            ! add omega+mu-eval (and broadening) on the diagonal
            if (i==j) then
              gij(i,i,isp) = gij(i,i,isp) + omega + ef - evlk(i,isp)
              ! if real frequencies add broadening.
              !!! CHECK IT: SIGN SHOULD DEPEND ON THE OCCUPANCY !!!
 !             if (.not. s_dmft%lsigim .and. broad /= 0) call rx('makegloc2: check broadening in g')
              if (.not. s_dmft%lsigim) gij(i,i,isp) = gij(i,i,isp) + (0d0,1d0)*broad
            endif
          enddo
        enddo
        ! perform inversion on (nlo:nhi) subblock (It's not efficient at the 1st iteration)
        call zinv(nev,gij(nlo,nlo,isp),nev)
C       call zprm('G',2,gij(nlo,nlo,isp),ndhamx,nlo:nhi,nlo:nhi)
        if (gpmode == 8) then
          if (iomg > nompr) cycle
          write(ifgk) gij(nlo:nhi,nlo:nhi,isp)
          cycle
        elseif (gpmode == 18) then
          if (iomg > nompr) cycle
          forall (i = nlo:nhi) gdiag(iomg,i,ikp,isp) = gij(i,i,isp)
          cycle
        elseif (gpmode == 20) then
          if (iomg > nompr) cycle
          do  i = nlo, nhi
            hkloc = omega + ef - evlk(i,isp) - 1/gij(i,i,isp)
            gdiag(iomg,i,ikp,isp) = hkloc
          enddo
          cycle
        endif

C   --- Local projections U+ G U ---

        do  cix = 1, ncix
C     In noncollinear case, isp=1 always => need internal ispc=1..2
C     ksp is the current spin index in both cases:
C     ksp = isp  in the collinear case
C     = ispc in the noncollinear case
C     ispc=1 for independent spins, and spin index when nspc=2
           do  ispc = 1, nspc
              ksp = min(max(ispc,isp),nsp)

              icix = iabs(s_dmft%icix(cix))
              nl1 = 2*s_dmft%l(icix)+1
              if (nl1 > ldcix) call rx('bug in sudmft')
              call zgemm('C','N',nl1,nev,nev,(1d0,0d0),
     .             dmftu(nlo,1,ksp,cix),nev,
     .             gij(nlo,nlo,isp),nev,(0d0,0d0),
     .             tmp1,ldcix)
              call zgemm('N','N',nl1,nl1,nev,(1d0,0d0),
     .             tmp1,ldcix,
     .             dmftu(nlo,1,ksp,cix),nev,(0d0,0d0),
     .             gkloc(1,1,ksp,cix),ldcix)
              call daxpy(2*ldcix*ldcix,1d0,gkloc(1,1,ksp,cix),1,gloc(1,1,ksp,cix),1)
           enddo                ! ksp
        enddo                   ! cix
      enddo                     ! isp

      deallocate(tmp1,gij)
      end subroutine makegloc2
