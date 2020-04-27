C#define DEBUG
      subroutine gworthphi(mode,s_site,s_spec,nl,nsp,nbas,nat,nphimx,nrmx,nrphi,phiv,phio)
C- Makes an orthonormal set of partial waves for matrix elements of the product basis
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa nr a rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt
Cio    Passed to:  *
Ci   mode  :1s digit
Ci         : 0 make/print overlap of partial waves, but do not calculate phio
Ci         : 1 generate phio
Ci         :10s digit
Ci         : 0 Orthonormalize using Cholesky decomposition
Ci         : 1 Orthonormalize by diagonalizing overlap matrix
Ci         : 2 Same as 1, but scale Kotani rdata4gw style for compatibility w/ old GW
Ci   nl    :(maximum lmxa) + 1 --- here a dimensioning parameter
Ci         :Old GW code used lmxamx
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   nat   :number of sites in basis with augmentation spheres
Ci   nphimx:dimensioning parameter : max number of valence partial waves of a particular l
Ci   nrmx  :leading dimension of phiv,phio.  Cannot be smaller than max(s_spec(:)%nr)
Ci   nrphi :number of valence partial waves for a particular l
Ci         :formerly named nvmax in old GW code
Ci   phiv  :valence partial waves.
Co Outputs
Co   phio  :orthonormalized form phiv
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   02 Jul 18
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer nl,nphimx,nrmx,nsp,nbas,nat,mode
      integer :: nrphi(0:nl-1,nat)
!     double precision phime(nrmx,0:nl-1,nphimx,nsp)
      double precision phiv(nrmx,0:nl-1,nphimx,nsp,nat),phio(nrmx,0:nl-1,nphimx,nsp,nat)
C ... Dynamically allocated arrays
C      real(8),allocatable:: rprodl(:,:,:),rprod(:,:)
C      real(8),allocatable:: eb(:),ebx(:),z(:,:),ovb(:,:)
C ... Local parameters
      integer,parameter:: nblmx=300
      integer ib,is,iat,isp,n1,l,n2,nr,i,lmxa,ir,nphil
      double precision ovv,xx(1)
      double precision ovphi(nphimx,nphimx,0:nl-1,nsp),ovphi2(nphimx,nphimx,0:nl-1)
      double precision rofi(nrmx),rwgt(nrmx)
      procedure(integer) :: nglob,iprint
      procedure(real(8)) :: dot3

      if (mod(mode,10) > 0) call dpzero(phio,size(phio))

      call info2(30,1,0,' gworthphi: orthogonalize partial waves by '//
     .  '%?#(n==0)#Cholesky Decomposition##%-1j%?#(n>0)#diagonalizing overlap##%-1j'//
     .  '%?#(n==2)# (Kotani scaling)##',mod(mode/10,10),2)

      iat = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        if (lmxa < 0) cycle
        iat = iat+1
        call rmeshprm(2+4,s_spec,nrmx,is,is,20,xx,xx,xx,xx,rofi,rwgt)
        nr = s_spec(is)%nr
        call dpzero(ovphi,size(ovphi))
        nphil = sum(nrphi(0:lmxa,iat))

        call info5(50,1,0,' start ib=%i  lmax=%i  nrphi=%i',ib,lmxa,nphil,4,5)

C       Return ovphi to orthogonalize. ovphi depends on mode
        i = 20001; if (mod(mode/10,10) > 0) i = 30001; if (mod(mode/10,10) == 2) i = 40001
        call prodphi(i,nl,nsp,nphimx,nrmx,nr,rofi,rwgt,nrphi(0,iat),
     .    phiv(1,0,1,1,iat),phiv(1,0,1,1,iat),ovphi,ovphi) ! last arg is summy

        if (mod(mode,10) == 0) goto 99 ! Skip orthogonalization

C   --- Orthogonalize partial waves ---
        do  isp = 1, nsp
          do  l = 0, lmxa

C       ... Overlap matrix diagonalized
            if (mod(mode/10,10) > 0) then

              if (mod(mode/10,10) == 2) then
                do  n2 = 1, nrphi(l,iat)
                do  n1 = 1, nrphi(l,iat)
                  ovphi(n1,n2,l,isp) = sqrt(1d0+0.1d0*n1) * ovphi(n1,n2,l,isp)
                enddo
                enddo
              endif

C             if (l == 0) call prmx('ovphi',ovphi(1,1,0,isp),nphimx,nrphi(l,iat),nrphi(l,iat))

              n1 = nrphi(l,iat)
              do  n2 = 1, n1
                forall (ir=1:nr) phio(ir,l,n2,isp,iat) = sum(phiv(ir,l,1:n1,isp,iat) * ovphi(1:n1,n2,l,isp))
              enddo

C       ... Cholesky decomposition of overlap matrix
            else

C             phio = phiv * transpose(L^-1)
              n1 = nrphi(l,iat)
              do  n2 = 1, n1
                forall (ir=1:nr) phio(ir,l,n2,isp,iat) = sum(phiv(ir,l,1:n1,isp,iat) * ovphi(n2,1:n1,l,isp))
              enddo

            endif

C#ifdef DEBUG
C       ... Debugging: overlap matrix of phio (should be orthonormal)
            call dpzero(ovphi2(1,1,l),nphimx**2)
            do  n1 = 1, nrphi(l,ib)
            do  n2 = n1, nrphi(l,ib)
              ovv = dot3(nr,phio(1,l,n1,isp,iat),phio(1,l,n2,isp,iat),rwgt)
              if (abs(ovv) < 1d-10) ovv = 0
              ovphi2(n1,n2,l) = ovv
              ovphi2(n2,n1,l) = ovv
            enddo
            enddo
C#endif

          enddo ! loop over l

C     ... Printout overlap of phio; confirm orthonormality
C#ifdef DEBUG
          if (iprint() >= 60) then
            call info0(55,1,0,' Overlap matrix of orthonormalized phi%N'//
     .        '  ib isp n1  n2   l=0        l=1 ...')
            do  l = 0, lmxa
              do  n1 = 1, nrphi(l,iat)
              do  n2 = n1, nrphi(l,iat)
                call info5(60,0,0,' %2,3i %,3i %,3i %n:2;9,9F',[ib,isp],n1,n2,nl,ovphi2(n1,n2,:))
              enddo
              enddo
            enddo
          endif
C#endif

        enddo ! Loop over spins

   99   continue
      enddo ! loop over basis

      end
