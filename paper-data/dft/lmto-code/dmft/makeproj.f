      subroutine makeproj(s_dmft,nbas,nsp,isp,ndhamx,ldcix,nlo,nhi,
     .  nlmax,nphimx,aus,ppn,nppn,n0,iq,qp,dmftu)
C- Add contribution from this qp to Haule's DMFT projectors
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft
Ci     Elts read:  nicix l ib
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iproj
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   isp   :current spin channel (1 or 2)
Ci   ndhamx:dimensioning parameter, largest hamiltonian dimension * nspc
Ci   nphimx:dmensions aus: global max # of partial waves of a given l
Ci   ldcix :dimension of largest cix block, dimensions dmftu
Ci   nlo,nhi :first, last eigenstate to include in projector
Ci   nlmax :largest (lmax+1)**2 over all sites : dimensions ausp
Ci   aus   :coefficients to (phi,phidot) projection of w.f
Ci   ppn   :NMTO pot pars; see potpus.f
Ci   nppn  :leading dimension of ppn
Ci   n0    :dimensioning parameter
Ci   iq    :index to current k-point
Ci   qp    :current k-point
Co Outputs
Co   dmftu :dmftu(ib,ichan,isp,icix) = factorized dmft projector for one spin, all correlated sites
Co         :ib = nlo:nhi
C          :ichan = 1:nlm
Co         :See Eq. 4.6 in Pisanti's PhD thesis
Cs Command-line switches
Cl Local variables
Cl   nicix : Number of (sites,l) to be singled out as correlated orbitals
Cr Remarks
Cr   Projectors and quantities use to make it (aux,ppn) are still locally collinear
Cu Updates
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu   22 Aug 15 (Pisanti) spin polarized
Cu   05 Oct 14 (Pisanti) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, intent(in) :: nbas,nsp,isp,ndhamx,ldcix,nlmax,nphimx,iq,n0,nppn,nlo,nhi
      real(8), intent(in) :: qp(3),ppn(nppn,n0,nbas)
      complex(8), intent(in)  :: aus(nlmax,ndhamx,nphimx,nsp,nbas)
      complex(8), intent(inout) :: dmftu(nlo:nhi,ldcix,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_dmft):: s_dmft
C ... Local parameters
      integer l,ib,cix,icix,m,nglob,stdo,nev,ncix
      double precision xv(nlmax*2)
      integer, parameter :: PRTU=50

C --- Setup ---
      stdo = nglob('stdo')
      ncix = s_dmft%ncix
      nev = nhi-nlo+1

C --- Make projector for each correlated block ---
      do  cix = 1, ncix

        ib   = s_dmft%ib(cix)
        icix = iabs(s_dmft%icix(cix))
        l    = s_dmft%l(icix)

        call dmftprojection(aus,ppn,nppn,n0,s_dmft%iproj,l,ib,nlmax,nphimx,
     .    isp,nsp,nbas,ndhamx,nlo,nhi,ldcix,cix,dmftu)
        do  m = 1, (2*l+1)
          xv(m) = sum(abs(dmftu(nlo:nhi,m,isp,cix))**2)
C         write(*,"(' Sum U^2(1:nev) after dmftprojection m : ', i4, f12.6)") m, xv(m)
        enddo
        call info5(PRTU,0,0,' ib=%,2i spin=%i cix=%i  Sum U^2%n;12,6D',ib,isp,cix,(2*l+1),xv)

C .... Write dosdmftu for density of states
C        ifi = fopna('dosdmftu',-1,0)
C        write(ifi,'("# kpt n.",i0,1x,3(1x,f12.8))') iq, qp
C        do iband = nlo, nhi
C          write(ifi,'(1x,i5)', advance='no') iband
C          do m = 1, nchan
C            write(ifi,'(x,(1x,f12.8))', advance='no') abs(dmftu(iband,m,isp,cix))**2
C          enddo
C          write(ifi,'("")')
C        enddo

      enddo                     !  loop over correlated orbitals
C     call zprm('dmftu',2,dmftu,nev,nev,ldcix*nsp)

C     call print_dmft(s_dmft,'proj_plot_dmftu','a',dmftu,nhi-nlo+1,ldcix,isp,nsp,nev,ncix,iq,qp)
!       call print_dmft(s_dmft,'proj_plot_dmftub','b',dmftu,nhi-nlo+1,ldcix,isp,nsp,nev,ncix,iq,qp)

      end subroutine makeproj

! !----Internal subroutines
!       subroutine print_dmft(s_dmft,file,fmt,dmft,ndhamx,ldcix,isp,nsp,nev,ncix,iq,qp)
!       implicit none
! C ... For structures
!       include 'structures.h'
!       type(str_dmft):: s_dmft
!
!       character :: fmt
!       character(len=*), intent(in) :: file
!       integer, intent(in) :: ndham,nev,ncix,iq,isp,nsp,ldcix
!       complex(8), intent(in) :: dmft(ndham,ldcix,nsp,ncix)
!       real(8), intent(in) :: qp(3)
!
!       integer :: u,iband,ib,ichan,fopna,cix,icix,l,ilm1
!
!       ib = 0; if (fmt == 'b') ib = 4
! C     u = fopnx(trim(file),2,ib,-1)
!       u = fopna(trim(file),-1,ib)  ! Use extension
!
!       if (fmt == 'a') then
!         write(u,'(a,1x,i0)') '# DMFT isp: ', isp
!         write(u,'("# kpt n.",i0,1x,3(1x,f12.8))') iq, qp
!       end if
!
!       do iband = 1, nev
!         if (fmt == 'a') write(u,'(1x,i5)', advance='no') iband
!         do cix = 1, ncix
!           ib   = s_dmft%ib(cix)
!           icix = iabs(s_dmft%icix(cix))
!           l    = s_dmft%l(icix)
!           ichan = 0
!           do ilm1 = l**2+1, (l+1)**2  ! Loop over m
!             ichan = ichan+1
!             if (fmt == 'a') then
!               write(u,'(2x,2(1x,f12.8))', advance='no')
!      .          dble(dmft(iband,ichan,isp,cix)),dimag(dmft(iband,ichan,isp,cix))
!             else if (fmt == 'b') then
!               write(u) real (dmft(iband,ichan,isp,cix)),aimag(dmft(iband,ichan,isp,cix))
!             end if
!           end do
!         end do
!         if (fmt == 'a') write(u,'("")')
!       end do
!
!       end subroutine print_dmft

C       subroutine print_aus(file, aus, ncorrch, ndham, nev, nlmax, ncix)
C       implicit none
C       character(len=*), intent(in) :: file
C       integer, intent(in) :: ndham, nev, nlmax, ncix
C       integer, intent(in) :: ncorrch(ncix)
C       complex(8), intent(in) ::  aus(nlmax,ndham)
C
C       integer :: u, iband, ichan, icix, l2, fopna
C
CC      open(unit=u, file=trim(file), action='write')
CC      u = fopnx(trim(file),2,0,-1)
C       u = fopna(trim(file),-1,0) ! append extension
C
C       write(u,'(2(a,1x,i0))')'# AUS ncix: ', ncix
C
C       do iband = 1, nev
C           write(u,'(1x,i5)', advance='no') iband
C           do icix = 1, ncix
C               l2 = ((ncorrch(icix)-1)/2)**2
C               do ichan = 1, ncorrch(icix)
C                   write(u,'(2x,2(1x,f12.8))', advance='no')real (aus(l2+ichan,iband)), aimag(aus(l2+ichan,iband))
C               end do
C           end do
C           write(u,'("")')
C       end do
C
C       call fclose(u)
C       end subroutine print_aus
