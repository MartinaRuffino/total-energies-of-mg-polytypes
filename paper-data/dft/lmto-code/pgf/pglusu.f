      subroutine pglusu(s_ham,nbas,npl,pgplp,ntab,iax,offH)
C- Setup for LU decomposition in layer geometry
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     bandw
Co     Allocated:  *
Cio    Elts passed:lncol
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see subas.f
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   28 Apr 00 made noncollinear
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,niax,npl,pgplp(6,-1:*)
      parameter (niax=10)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer ntab(nbas),iax(niax,1),offh(n0H,nkap0,nbas+1)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Local parameters
      integer kl,ku,kli,ib1,ib2,lgunit,iprint,nspc,lrel
      logical lplh
      procedure(integer) nglob

      nspc = min(1+s_ham%lncol,2)
      lrel = mod(nglob('lrel'),10)

      lplh = .false.
      lplh = .true.
      if (lplh) then
C   ... Get kl assembling hamiltonian by layers
        call pglus1(0,npl-1,npl,pgplp,nspc,lrel,kl,ku)
      else
        if (nspc == 2) call rx('branch not implemented for nspc=2')
C   ... Get kl from straight 2D Bloch sum of layers 0..npl-1
        call spbndw(nbas,ntab,iax,offH,kl,ku)
C       kl,ku at least as large as size of 0 PL (left BC)
        ib1 = 1
        ib2 = pgplp(1,0)
        kli = offh(1,1,ib2+1) - offh(1,1,ib1) - 1
        kl = max(kl,kli)
        ku = max(ku,kli)
C       kl,ku at least as large as size of n PL (right BC)
        ib1 = pgplp(1,npl-2)+1
        ib2 = pgplp(1,npl-1)
        kli = offh(1,1,ib2+1) - offh(1,1,ib1) - 1
        kl = max(kl,kli)
        ku = max(ku,kli)
      endif
      kl = max(kl,ku)

      if (iprint() >= 30) call
     .  awrit1('%N pglusu: band subdiagonal = %i',' ',80,lgunit(1),kl)
      s_ham%bandw = kl
      end
      subroutine pglus1(ip1,ip2,npl,pgplp,nspc,lrel,kl,ku)
C- Return width of band, PL representation
C ----------------------------------------------------------------------
Ci Inputs
Ci   ip1,ip2: range of PL to consider
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f,subas.f)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   lrel  :0 for non-relativistic
Ci          1 for scalar relativistic
Ci          2 for Dirac equation
Co Outputs
Co   kl    :size of subdiagonal
Co   ku    :size of superdiagonal
Cr Remarks
Cu Updates
Cu   09 Jul 03 Add one to kl,ku if lrel=2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer ip1,ip2,kl,ku,npl,nspc,lrel,pgplp(6,-1:*)
C ... Local parameters
      integer kli,kui,ipl

      kl = 0
      ku = 0
      do  10  ipl = ip1, ip2
        kui = nspc*pgplp(4,ipl) - 1
        kli = kui
        if (ipl > 0) kui = kui + pgplp(4,ipl-1)
        if (ipl < npl-1) kli = kli + pgplp(4,ipl+1)

        kl   = max(kl,kli)
        ku   = max(ku,kui)
        if (lrel == 2) then
          kl = kl+1
          ku = ku+1
        endif

C        print 333, ipl,kli,kui
C  333   format(i4,2x,2i4,2x,2i4,2x,2i4)
   10 continue
      end
