      subroutine dmftprojection(aus,ppn,nppn,n0,iproj,l,site,
     .  nlmax,nphimx,isp,nsp,nbas,ndhamx,nlo,nhi,ldcix,cix,dmftu)
C- Factorized DMFT projector for one l
C ----------------------------------------------------------------------
Ci Inputs
Ci   aus   :val,slo of w.f. at MT sphere surface (mkusq.f)
Ci   ppn   :potential parameters (potpus.f)
Ci   nppn  :dimensions ppn
Ci   n0    :dimensioning parameter
Ci   iproj :Projector type
Ci   l     :l quantm number for this correlated orbital
Ci   site  :site index for this correlated orbital
Ci   nlmax :largest (lmax+1)**2 over all sites : dimensions ausp
Ci   nphimx:maximum number of radial functions that make up projector
Ci   nlohi :first, last band to include into projector
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   isp   :current spin index
Ci   nbas  :size of basis
Ci   ndhamx:dimensioning parameter, largest hamiltonian dimension * nspc
Ci   ldcix :dimensioning parameter, max number of correlated orbitals in a block
Ci   cix   :index to current correlated block
Co Outputs
Co   dmftu :Factorized projector for one spin.  dmftu assumed to be zero on entry.
Co         :Data is into dmft(1:nev,1:nlm,isp,site)
Co         :where nev = nhi-nlo+1 and nlm = ilmax-ilmin+1
Co         :Thus first band index (nlo) and first orbital index (l**2+1)
Co         :are written to dmft(1,1,isp,cix)
Cl Local variables
Cl         :
Cr Remarks
Cr   Projectors and quantities use to make it (aux,ppn) are still locally collinear
Cu Updates
Cu   03 Oct 17 Rationalize dimensions for dmftu
Cu   27 Feb 16 (MvS) 1st attempt at redesign for more general cix blocks
Cu      Oct 14 (Pisanti) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: l,site,nphimx,iproj,isp,nsp,ndhamx,nlo,nhi,
     .                       nppn,n0,nlmax,nbas,ldcix,cix
      real(8), intent(in) :: ppn(nppn,n0,nsp,nbas)
      complex(8), intent(in)    :: aus(nlmax,ndhamx,nphimx,nsp,nbas)
      complex(8), intent(inout) :: dmftu(nlo:nhi,ldcix,nsp,cix)
C ... Local parameters
      integer k(3)  ! Points to normalization of appropriate partial wave in pp
      integer :: ivec,ilm1,ilmax,ilmin,ichan,ipqn,nev
      complex(8) :: sm1,sm2,ff

      k(1) = 2 ! Normalization of phi
      k(2) = 7 ! Normalization of phidot
      k(3) = 8 ! Normalization of local orbital.  Should be 0 if no LO

      nev = nhi-nlo+1
      ilmin = l**2+1; ilmax = ilmin+2*l  ! Index to 1st, last m for this l
      do  ivec = nlo, nhi ! over bands to be included in sum
C   ... Projection 1, Eq. 17 in arXiv:0907.0195v2
        if (abs(iproj) == 1) then
          ichan = 0              ! ilm index, shifted to start at 1
          do  ilm1 = ilmin, ilmax ! Loop over m
            ichan = ichan+1
            dmftu(ivec,ichan,isp,cix) = aus(ilm1,ivec,1,isp,site)*ppn(2,l+1,isp,site)
          enddo
C   ... Projection 2, Eq. 21 in arXiv:0907.0195v2
        elseif (abs(iproj) == 2) then
          ichan = 0
          do  ilm1 = ilmin, ilmax ! Loop over m
            ichan = ichan+1
            sm1 = 0
            sm2 = 0
            ff = aus(ilm1,ivec,1,isp,site)
            if (abs(ff) > 0) ff = ff/abs(ff)
            do  ipqn = 1, nphimx
              sm1 = abs(aus(ilm1,ivec,ipqn,isp,site))**2
              sm2 = sm2 + sm1*ppn(k(ipqn),l+1,isp,site)
            enddo
            dmftu(ivec,ichan,isp,cix) = ff*sqrt(sm2)
          enddo ! Loop over m
        endif ! projector type 2
      enddo ! Loop over eigenstates

      end subroutine dmftprojection
