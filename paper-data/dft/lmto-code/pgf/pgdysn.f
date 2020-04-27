      subroutine pgdysn(s_ham,s_pot,s_str,s_site,plat,isp,ip0,npl,noL,noR,nrpt,
     .  vpl,mode,pgplp,qp,lgii,gii)
C- Make crystal GF from surface GF and neighboring crystal GF
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham lncol neula lgen3 lham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plhamnc plham plhamso
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfnc palp
Cio    Passed to:  plhamnc plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plhamnc plham
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plhamnc plham plhamso
Ci Inputs:
Ci   plat  : primitive lattice vectors, in units of alat
Ci   isp   : current spin channel (1 or 2)
Ci   ip0   : current PL to embed :
Ci   npl   : number of principal layers
Ci   noL   : max number of left off-diagonal GF g(ip0..ip0+nol,ip0)
Ci   noR   : max number of right off-diagonal GF g(ip0,ip0..ip0+noR)
Ci   nrpt  : number of repeat layers
Ci   vpl   : used in adiabatic approximation; see pvemb and pginns
Ci   mode  : specifies type of GF to be created
Ci        1s digit signifies:
Ci         0 create GF
Ci         1 create GF^-1
Ci       10s digit signifies:
Ci         0 do not calculate any GF
Ci         4 convert gii~(ip0) to gii(ip0) using gii(ipR)
Ci           if ip0 = -2 make off-diagonal (*)(-2) only
Ci         5 convert gii~(ip0) to gii(ip0) using gii(ipL)
Ci   pgplp : table of PL-dependent quantities (see pgfset)
Ci   qp    : k-point (NB: projection of qp onto PL direction not used).
Cio Inputs/Outputs:
Cio  gii   : Diagonal for new principal layer, at ip0
Coo   lgii  : attributes of gii generated; see lgupac for conventions
Cl Lcal variables
Cl   ipL   : layer ip0-nrpt+1, mode=40
Cl           layer ip0-1, mode=50
Cl   ipR   : layer ip0+1, mode=40
Cl           not used, mode=50
Cl   iprt  : print verbosity
Cr Remarks
Cr   Partitioning into lower and higher blocks not implemented.
Cu Updates
Cu   02 Dec 15 mode 5 extended to noncollinear case
Cu   27 Sep 15 First cut at extending to the noncollinear case
Cu   02 Sep 15 complete migration to f90 pointers
Cu   08 Jul 13 replaced f77 pointers with f90 ones except passed pointers
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,npl,ip0,noL,noR,nrpt,mode,pgplp(6,-1:*),lgii(-2:*)
      double precision plat(3,3),qp(3),vpl(-1:npl)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wk2(:)
      complex(8), allocatable :: wk3(:)
      complex(8), allocatable :: sii(:)
      complex(8), allocatable :: r(:)
      complex(8), allocatable :: s00(:,:,:,:),wk(:,:,:,:)
      complex(8), pointer :: gLL(:,:),gRR(:,:)
C ... Local parameters
      logical, parameter :: sanity=.true.
      integer, parameter :: PRTG=80
      integer lds,ld0,ldl,ldr,ldw,mode2,iprt,pgdim,nspc,lncol
      integer ndg0,ndgl,ndgr,ofl,ofr,of0,ndim,kpl,ipL,ipR
      integer ld0x,ldlx,ldrx,ldsx,ndg0x,ndgrx,ndglx,ldwx
      double precision xv(1)
      character outs*80

C --- Setup ---
      call tcn('pgdysn')
      mode2 = mod(mode/10,10)
      if (mode2 == 0) return
      if (mode2 == 4) then
        ipL = ip0-nrpt+1
      elseif (mode2 == 5) then
        ipL = ip0-1
      endif
      ipR = min(ip0+1,npl)
      gLL => gii(ipL)%gll
      gRR => gii(ipR)%gll
      ndgL = pgplp(3,max(ipL,-1))
      ndgR = pgplp(3,min(ipR,npl))
      ndg0 = pgplp(3,max(ip0,-1))
      ofL = pgplp(6,max(ipL,-1))
      ofR = pgplp(6,min(ipR,npl))
      of0 = pgplp(6,max(ip0,-1))
      nspc = s_ham%ldham(4)
      lncol = s_ham%lncol

      call getpr(iprt)
      if (ip0 < -2 .or. ip0 > npl+1) call fexit2(-1,111,
     .  ' Exit -1 PGDYSN: expected ipl to be >=0 and <%i'//
     .  ' but found ipl=%i',npl,ip0)
      ldl = pgplp(4,max(ip0-1,-1))
      ld0 = pgplp(4,max(ip0,-1))
      ldr = pgplp(4,min(ip0+1,npl))
      lds = ldl+ld0+ldr
      ld0x = ld0*nspc
      ldlx = ldl*nspc
      ldrx = ldr*nspc
      ldsx = lds*nspc
      ndg0x = ndg0*nspc
      ndgRx = ndgR*nspc
      ndgLx = ndgL*nspc

C --- Strux connecting ip0 to adjacent layers ---
      allocate(s00(ld0,nspc,lds,nspc))
      if (nspc == 2) then
        allocate(wk(ld0,nspc,lds,nspc))
      else
        allocate(wk(1,1,1,1))
      endif
      call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ip0,npl,ld0,lds,pgplp,qp,wk,s00)
      deallocate(wk)

C      if (lncol /= 0) then
C        kcplxl = 1              ! plham requires kcplx=1
C        allocate(wk(ld0,nspc,lds,nspc))
C        call plham(1,s_ham,s_pot,s_str,s_site,plat,isp,kcplxl,max(ip0,-1),pgplp,qp,lds,ld0,wk)
CC       call awrit1(' P-S(%i) from plham:',outs,len(outs),0,max(ip0,-1))
CC       call zprm(outs,2+0,wk,ld0x,ld0x,ldsx)
C
CC   ... Rearrange (spin,orbital) indexing
C        call dpzero(s00,2*size(s00))
C        call pgflu2(0,kcplxl,nspc,ld0,ldl,ldr,lds,0,0,2,2,ld0*nspc,wk,s00)
CC       call awrit1(' spin reordered P-S(%i):',outs,len(outs),0,max(ip0,-1))
CC       call zprm(outs,2+0,s00,ld0x,ld0x,ldsx)
C
CC       Convert to kcplx=0 format
C        call ztoyy(s00,ld0x,ldsx,ld0x,ldsx,kcplxl,0)
C        deallocate(wk)
C      else
C        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,max(ip0,-1),pgplp,qp,lds,ld0,s00)
C      endif

C --- Make some GF at ip0, depending on mode and what is available ---
C ... Case left sif, repeating potential layers ...
      if (mode2 == 4) then
        ldw = max(ld0,ldr); ldwx = ldw*nspc
        allocate(wk3(ld0x*ldwx))
C   ... Make g(ip0+1..ipl+noL)(ip0) and store into appropriate gii
        if (noL > 0) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          call pgoff3(sanity,ip0,noL,pgplp,wk3,s00,lgii,gii)
        endif
C   ... Make g(ip0)(ip0+1..noR)
        if (noR > 0) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          ndim = pgdim(0,ip0,npl,ip0+1,ip0+noR,pgplp)
          call pgdys2(sanity,2,'N',lds,ld0,ldr,ld0,ndg0,ndgR,ndg0,of0,
     .      ofR,of0+ld0,ndim,ld0+ldl,wk3,s00,lgii(ip0),
     .      gii(ip0)%gll,lgii(ipR),gii(ipR)%gll,gii(ip0)%gll)
        endif
C   ... Convert g~(ip0)(ip0) -> g(ip0)(ip0)
        if (nrpt > 1) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          allocate(wk2(ld0*ld0))
          if (ndg0 /= ld0) call rx('pgdysn fix copy gii to wk2')
          call dcopy(ld0*ld0*2,gii(ip0)%gll,1,wk2,1)
        endif
C       Use g(ip0)(ip0+1) if we already have it
        if (noR > 0) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          call pgdys1(2,lds,ld0,ldr,ndg0,ndgr,ndg0,of0,ofR,of0+ld0,
     .      ld0+ldl,ldw,wk3,s00,lgii(ipR),gRR,gii(ip0)%gll,
     .      lgii(ip0),gii(ip0)%gll)
C       Otherwise use wk to create g(ip0)(ip0+1)
        else
          allocate(wk(ld0,nspc,lds,nspc))
          call pgdys1(0,ldsx,ld0x,ldrx,ndg0x,ndgrx,ldsx,of0,ofR,0,ld0x+ldlx,ldwx,
     .      wk3,s00,lgii(ipR),gRR,wk,lgii(ip0),gii(ip0)%gll)
          deallocate(wk)
        endif
C   ... Case repeating potential layers ...
        if (nrpt > 1) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
!         orvec = 2*ld0*ld0*i1mach(18)
          allocate(sii(ld0*ld0*3))
          allocate(r(2*ld0))
          call pginns(ld0,nrpt-1,r,gii(ip0)%gll(1,ld0+1),xv,vpl,
     .      lgii(ipL),gLL,xv,sii,2,wk2,lgii(ip0),gii(ip0)%gll)
          deallocate(sii,r,wk2)
        endif
        deallocate(wk3)

C --- Convert gii(ip0,ip0)~ to gii(ip0,ip0) using gii(ipL,ipL) ---
      elseif (mode2 == 5) then
        ldw = max(ld0,ldl)
        allocate(wk3(ld0x*ldw*nspc))
C   ... Make g(ipL)(ipL+1..ipL+noR)
        if (noR > 0) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          ndim = pgdim(0,ipL,npl,ipL+1,ipL+noR,pgplp)
          call pgdys2(sanity,2,'N',lds,ldl,ld0,ld0,ndgL,ndg0,ndgL,ofL,
     .      of0,ofL+ld0,ndim,ld0+ldl,wk3,s00,lgii(ipl),
     .      gii(ipL)%gll,lgii(ip0),gii(ip0)%gll,gii(ipL)%gll)
        endif
C   ... Make g(ip0-1)(ip0) and g(ip0)(ip0)
        if (noL > 0) then
          if (nspc == 2) call rx('pgdysn not ready for nspc == 2')
          if (noL > 1) call rx('pgdysn nol>1 not implemented')
          call pgdys1(0,lds,ld0,ldl,ndg0,ndgL,ndg0,of0,ofL,of0-ldl,ld0,
     .      ldw,wk3,s00,lgii(ipL),gLL,gii(ip0)%gll,
     .      lgii(ip0),gii(ip0)%gll)
C       Otherwise, use wk to hold gii(ip0-1,gii(ip0))
        else
          allocate(wk(ld0,nspc,lds,nspc))
          call pgdys1(0,ldsx,ld0x,ldlx,ndg0x,ndgLx,ldsx,of0,ofL,0,ld0x,ldw*nspc,
     .      wk3,s00,lgii(ipL),gLL,wk,lgii(ip0),gii(ip0)%gll)
          deallocate(wk)
        endif
        if (nrpt > 1) call rxi('pgdysn: mode=%i, nrpt>1 not implemented',mode)
      else
        call rx('pgdysn: bad mode')
      endif

C --- (debugging) print gii ---
      if (iprt >= PRTG/1) then
        call awrit2('g for layer %i by Dyson, mode=%i',outs,len(outs),0,ip0,mode)
        call yprm(outs,2,gii(ip0)%gll,ld0x*ndg0x,ld0x,ld0x,ndg0x)
        if (mode2 == 4 .and. nrpt > 1) then
          call awrit1('gii for layer %i',outs,len(outs),0,ip0-nrpt+1)
          call rx('pgdysn fix pointer gLL in call to pginns')
          call yprm(outs,2,gLL(1,of0),ld0x*ld0x,ld0x,ld0x,ld0x)
        endif

        if (mode2 == 4) then
          do  kpl = ip0+1, min(ip0+noL,npl)
            call awrit2('updated offg for layer %i, mode=%i',outs,len(outs),0,kpl,mode)
            ndg0 = pgplp(3,kpl)
            ld0 = pgplp(4,kpl)
            call yprm(outs,2,gii(kpl)%gll,ld0*ndg0,ld0,ld0,ndg0)
          enddo
        elseif (noL > 0) then
          kpl = 6
          call awrit2('updated offg for layer %i, mode=%i',outs,
     .      len(outs),0,kpl,mode)
          ndg0 = pgplp(3,kpl)
          ld0 = pgplp(4,kpl)
          call yprm(outs,2,gii(kpl)%gll,ld0*ndg0,ld0,ld0,ndg0)
        endif
      endif

      deallocate(s00)
      call tcx('pgdysn')
      end
      subroutine pgdys1(mode,lds,lda,ldb,nda,ndb,ndab,offa,offb,offab,
     .  offs,ldw,wk,sab,lgbb,gbb,gab,lgaa,gaa)
C- Make interfacial g_ii from one semiinfinite and one proper GF
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors, in units of alat
Ci   mode  :0 calculate gab and convert gaa~ to gaa
Ci          1 calculate gab only
Ci          2 gab is input ; convert gaa~ to gaa
Ci   sab   :strx connecting layers a and b (in sab(1,offs))
Ci   gaa   :semi-infinite GF gaa~
Ci   gbb   :proper diagonal GF, b layer
Ci   nda   :second dimension for gaa
Ci   ndab  :second dimension for gab
Ci   ndb   :second dimension for gbb
Ci   wk,ldw:workspace, size (lda,ldw,2); ldw must be >= max(lda,ldb)
Co Outputs:
Co   gab   :off-diagaonal GF gab
Co   gaa   :(mode ne 1) proper crystal diagonal GF gaa
Co   sab   :is OVERWRITTEN
Cr Remarks
Cr   Uses special case of general formulas from Dyson equation
Cr     g_AB = g_Aa~ Sab gbB = g_Aa Sab gbB~
Cr     g_BA = g_Bb~ Sba gaA = g_Bb Sba gaA~
Cr     g_AA = g_AA~  +  g_Ab  Sba g_aA~
Cr            g_AA~  +  g_Aa~ Sab g_bA
Cr   Whether gaa is right- or left- semi-infinite and gbb is the proper
Cr   GF depends only on which s is used: s0r when gaa~ is on the left;
Cr   use s0l when gaa~ is the on the right.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lds,lda,ldb,offs,ldw,lgaa,lgbb,nda,ndab,ndb,mode,
     .  offa,offb,offab
      double precision sab(lda,lds,2),wk(lda,ldw,2),
     .  gaa(lda,nda,2),gbb(ldb,ndb,2),gab(lda,ndab,2)
C ... Local parameters
      logical, parameter :: sanity=.true.

      call tcn('pgdys1')

C --- Make gab ---
      if (mode < 2)
     .  call pgdys2(sanity,2,'N',lds,lda,ldb,lda,nda,ndb,ndab,offa,
     .  offb,offab,ldb,offs,wk,sab,lgaa,gaa,lgbb,gbb,gab)

      if (mode == 1) return

C --- S_ab(1) <- delta g_aa = g_ab S_ba g_aa~ = g_ab S_ab+ g_aa~ ---
      call yygemm('N','C',lda,lda,ldb,1d0,gab(1,1+offab,1),gab(1,1+offab,2),lda,
     .  sab(1,1+offs,1),sab(1,1+offs,2),lda,0d0,wk,wk(1,1,2),lda)
      call yygemm('N','N',lda,lda,lda,1d0,wk,wk(1,1,2),lda,
     .  gaa(1,1+offa,1),gaa(1,1+offa,2),lda,0d0,sab(1,1,1),sab(1,1,2),lda)
C     call yprm('g_ab S_ba gaa~',2,sab(1,1,1),lda*nda,lda,lda,lda)

C --- Make g_aa = gaa~ + g_ab s_ba g_aa~ ---
      call daxpy(lda*lda,1d0,sab(1,1,1),1,gaa(1,1+offa,1),1)
      call daxpy(lda*lda,1d0,sab(1,1,2),1,gaa(1,1+offa,2),1)
      call lgupac(lgaa,'pt',0,3,0,0,0,0)
      call lgupac(lgaa,' pgdys1 created',6,0,0,0,0,0)

      call tcx('pgdys1')
      end
      subroutine pgdys2(sanity,mode,cc,lds,lda,ldb,ldab,nda,ndb,ndab,
     .  offa,offb,offab,n,offs,wk,sab,lgaa,gaa,lgbb,gbb,gab)
C- Make an off-diagonal GF a surface diagonal and a crystal diagonal GF
C ----------------------------------------------------------------
Ci Inputs:
Ci   sanity   if T, make sanity checks
Ci   sab:     strx connecting layers a and b (in sab(1,offs))
Ci   cc:      'N' or 'C':  use Sab or Sab+
Ci   lds,offs:leading dimension and offset to s connecting a to b
Ci   gaA:     aA GF on a side of ab interface
Ci            Requires gaa(1..lda,1+offa..1+lda)
Ci   gbB:     bB GF on b side of ab interface
Ci            Requires gbb(1..ldb,1+offb..1+n)
Ci   lda,nda: first and second dimension of gaA
Ci            Also lda is leading dimension of Sab and
Ci            is dimension of product gaa Sab
Ci   ldb,ndb: first and second dimension of gbB
Ci            Also ldb is dimension of product S_ab gbb
Ci   offa,offb: offset to diagonal parts of gaA and gbB
Ci   ldab,ndab:first and second dimension of gAB
Ci   mode:    0 make nothing
Ci            1 wk = gaa S_ab and exit
Ci            2 Make wk as in (1) and gaB = (g_aa S_ab) gbB
Ci            3 Take wk as input; make gaB = wk gbB
Ci            4 wk = S_ab gbb and exit
Ci            5 Make wk as in (4) and gAb = g_Aa (S_ab gbb)
Ci            6 Take wk as input; gAb = g_Aa wk
Ci   n:       calculate gaB(*,1+offab..n+offab)
Ci   wk       workspace, size (lda*ldb*2)
Co Outputs:
Co   gaB:    gaB(1..lda,1+offab..n+offab) is calculated
Co   wk:     g_aa S_ab
Cr Remarks
Cr   From the Dyson equation, the general formulas apply
Cr     g_i,j = g~i,k S_k,k+1 g_k+1,j  = g_i,k S_k,k+1 g~k+1,j
Cr   pgdys2 calculates
Cr     g_aB = g~_aa Sab gbB or g_aa Sab g~bB
Cr   All offsets offa,offb,offab,offs are numbered from zero.
C ----------------------------------------------------------------
      implicit none
      logical sanity
      character*1 cc
      integer lds,lda,ldb,ldab,offs,lgaa,lgbb,nda,ndab,ndb,
     .  offa,offb,offab,n,mode
      double precision sab(lda,lds,2),wk(lda,ldb,2),
     .  gaa(lda,nda,2),gbb(ldb,ndb,2),gab(ldab,ndab,2)
      logical lok

C     print *, mode,lds,lda,ldb,nda,ndb,ndab,offa,offb,offab,n,offs

C ... Error checks
      if (sanity) then
        call lgupac(lgaa,' pgdys2 using gaa:',6,0,0,0,0,0)
        call lgupac(lgbb,' pgdys2 using gbb:',6,0,0,0,0,0)
        lok = mod(lgaa/2,4) == 1 .or. mod(lgaa/2,4) == 2 .or.
     .        mod(lgbb/2,4) == 1 .or. mod(lgbb/2,4) == 2
        call rxx(.not.lok,'pgdys2: gaa or gbb must be semi-infinite GF')
C       call yprm('g_aa',2,gaa(1,1+offa,1),lda*nda,lda,lda,lda)
C       call yprm('g_bB',2,gbb(1,1+offb,1),ldb*ndb,ldb,ldb,n)
C       call yprm('S_ab',2,sab(1,1+offs,1),lda*lds,lda,lda,ldb)
      endif

C --- Make gaB and/or similar quantities ---
C ... wk = g_aa S_ab
      if (mode == 1 .or. mode == 2)
     .  call yygemm('N',cc,lda,ldb,lda,1d0,gaa(1,1+offa,1),gaa(1,1+offa,2),lda,
     .  sab(1,1+offs,1),sab(1,1+offs,2),lda,0d0,wk,wk(1,1,2),lda)
C     call yprm('g_aa s_ab',2,wk,lda*ldb,lda,lda,ldb)

C ... gaB = (g_aa S_ab) gbB
      if (mode == 2 .or. mode == 3) then
        call yygemm('N','N',lda,n,ldb,1d0,wk,wk(1,1,2),lda,gbb(1,1+offb,1),
     .    gbb(1,1+offb,2),ldb,0d0,gab(1,1+offab,1),gab(1,1+offab,2),ldab)
C       call yprm('g_aB',2,gab(1,1+offab,1),ldab*ndab,ldab,ldab,n)
      endif

C --- Make gAb and/or similar quantities ---
C ... wk = S_ab gbb
      if (mode == 4 .or. mode == 5)
     .  call yygemm(cc,'N',lda,ldb,ldb,1d0,sab(1,1+offs,1),sab(1,1+offs,2),lda,
     .  gbb(1,1+offb,1),gbb(1,1+offb,2),ldb,0d0,wk,wk(1,1,2),lda)

C ... gAb = g_Aa S_ab gbb
      if (mode == 5 .or. mode == 6) then
        call yygemm('N','N',n,ldb,lda,1d0,gaa(1,1+offa,1),gaa(1,1+offa,2),lda,
     .    wk,wk(1,1,2),lda,0d0,gab(1,1+offab,1),gab(1,1+offab,2),ldab)
C       call yprm('g_Ab',2,gab(1,1+offab,1),ldab*ndab,ldab,n,ldb)
      endif

      end
      subroutine pgoff3(sanity,ipl,nol,pgplp,wk,sab,lgii,gii)
C- Off-diagonal g(ipl+1..ipl+noL)(ipl), layer by layer from Dyson
C ----------------------------------------------------------------
Ci Inputs:
Ci   sanity   if T, make sanity checks
Ci   ipl,nol  make g(ipl+1..ipl+nol,ipl)
Ci   sab:     strux s(ipl+1,ipl).
Ci            If ipl>-2, have only s(ipl,ipl+1); use that instead.
Ci   wk:      dimensioned at least (lda,ldb)*2, with
Ci            lda=max(lda) for ipl+1..ipl+nol
Co Outputs
Co   gii     blocks to the off-diagonal gii(ipl+1,ipl), gii(ipl+nol,ipl)
Co            updated
Cr Remarks
C    Uses g_i+n,i+1 S_i+1,i g~i,i
C ----------------------------------------------------------------
      use structures, only : s_lgf
      implicit none
C ... Passed parameters
      logical sanity
      integer ipl,nol,pgplp(6,-1:*),lgii(-2:*)
C     double precision sab(lda,lds,2),wk(lda,ldb,2)
      double precision sab(*),wk(*)
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*)
C ... Local parameters
      integer ldl,lda,ldb,nda,ndb,lmode,offa,offb,offab,ipn,nan,lds,offs
      character*80 outs
      character *1 cc

      lmode = 5
      ldb  = pgplp(4,max(ipl,-1))
      ndb  = pgplp(3,max(ipl,-1))
      offb = pgplp(6,max(ipl,-1))
      ldl  = pgplp(4,max(ipl-1,-1))
      lda  = pgplp(4,ipl+1)
      lds  = ldl+ldb+lda
C ... Want sab = s(ipl+1,ipl).  if ipl>-2, have only s(ipl,ipl+1)
      offs = ldl+ldb
      cc = 'C'
C ... if ipl is -2, have directly s(ipl+1,ipl) (not needed)
C      if (ipl == -2) then
C        offs = ldl
C        cc = 'N'
C      endif

C ... nan is the number of rows between PL i+n and i
      nan = 0
      do  10  ipn = ipl+1, ipl+noL

C   ... Dimensioning parameters for row i+n
        lda  = pgplp(4,ipn)
        nan  = nan + lda
        nda  = pgplp(3,ipn)
        offab= pgplp(6,ipn) - nan
        offa = offab + lda

C        print *, sanity,lmode,lds,lda,ldb,lda,nda,ndb,nda,
C     .    offa,offb,offab,offs,wk,sab

        call pgdys2(sanity,lmode,cc,lds,lda,ldb,lda,nda,ndb,nda,
     .    offa,offb,offab,ldb,offs,wk,sab,lgii(ipn),gii(ipn)%gll,
     .    lgii(ipl),gii(ipl)%gll,gii(ipn)%gll)

        call awrit2(' pgoff3 created g(%i)(%i)',outs,len(outs),0,ipn,ipl)
        call lgupac(lgii(ipn),outs,6,0,0,0,0,0)
C       call yprm(outs,2,gii(ipn)%gll,lda*nda,lda,lda,nda)

        lmode = 6

   10 continue

      end

      subroutine pgofLR(sanity,ld0,ldl,ldr,ldm1,sab,wk,lgii,
     .  ldji,gji,ldij,gij,ldRL,ndRL,giL,gLi)
C- Make g-1,i and gi,-1 from g~_-1,j-1, g_j,i and g~_j-1,-1, g_i,j
C ----------------------------------------------------------------
Ci Inputs:
Ci   sanity: T, make sanity checks
Ci   ld0   : Dimension of g at current layer
Ci   ldl   : Dimension of layer left adjacent to current layer
Ci   ldr   : Dimension of layer right adjacent to current layer
Ci   ldm1  : Dimension of g at -1 layer
Ci   sab   : strux s(ipl,ipl-1..ipl+1), ipl = current layer
Ci   wk    : workspace dimensioned at least (lda,ldL)*2
Ci   lgii  : attributes of gii?
Ci   ldji  : first and second dimensions of gji
Ci   gji   : Green's function connecting layers j and i
Ci         : (typically layer j = layer i; see Remarks)
Ci   ldij  : first and second dimensions of gij
Ci   gij   : Green's function connecting layers i and j
Ci         : (typically layer j = layer i; see Remarks)
Ci   ldRL,ndRL leading two dimensions of giL and gLi
Cio Inputs/Outputs
Cio  giL       On input, surface GF g~_j-1,-1, see Remarks
Cio            On output, interface g_i,-1, see Remarks
Cio  gLi       On input, surface GF g~_-1,j-1; see Remarks
Cio            On output, g_-1,i; see Remarks
Cr Remarks
Cr   pgofLR uses Dyson equation:
Cr     g_i,-1 = g~i,j  S_j,j-1  g_j-1,-1
Cr            = gi,j   S_j,j-1  g~_j-1,-1
Cr     g_-1,i = g~-1,i S_j-1,j  g_j,i
Cr            = g-1,i  S_j-1,j  g~_j,i
Cr   It is assumed for convenience layers j and i have the same
Cr   matrix elements; then j's value need not be explicit.
Cr   If g_ij is xtal g, resulting g_iL is xtal g_iL
Cr   If g_ij is g~_ii,  resulting g_iL is surface g~_iL
Cr   Example:
Cu Updates
Cu   20 Feb 02 bug fix when ipl=npl.  Altered argument list.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical sanity
      integer ldRL,ndRL,lgii,ldji,ldij,ld0,ldl,ldr,ldm1
      double precision gji(ldji,ldji,2),giL(ldRL,ndRL,2)
      double precision gij(ldij,ldij,2),gLi(ndRL,ldRL,2)
C     double precision sab(*),wk(*)
      double precision sab(ld0,ld0+ldl+ldr,2),wk(*)
C ... Local parameters
      integer offs,lds

C ... Parameters for i layer -> ld0
C     ldi = ld0 ! pgplp(4,ipl)
C ... Parameters for j-1 layer
C     ldl  = pgplp(4,max(ipl-1,-1))
C ... Parameters for -1 layer
C     ldm1  = pgplp(4,-1)
C ... Parameters for Sab -> lds
C     lds  = ldl + pgplp(4,min(ipl+1,npl)) + ld0
      lds = ld0+ldl+ldr
      offs = ld0

C ... Error checks
      if (sanity) then
        call lgupac(lgii,' pgofLR using gaa:',6,0,0,0,0,0)
C       print *, 'pgofLR: ipl is', ipl
C       call yprm('g_ij',2,gij,ldij*ldij,ldij,ld0,ld0)
C       call yprm('g_ji',2,gji,ldji*ldji,ldji,ld0,ld0)
C       call yprm('S_ab',2,sab(1,1+offs,1),ld0*lds,ld0,ld0,ldl)
C       call yprm('g_j-1L',2,giL,ldRL*ndRL,ldRL,ldl,ldm1)
C       call yprm('g_Lj-1',2,gLi,ndRL*ldRL,ndRL,ldm1,ldl)
      endif

C ... wk(ld0,ldm1) = Sab initial-giL; sab(ld0,lds) = s(ipl,ipl-1)
      call yygemm('N','N',ld0,ldm1,ldl,1d0,
     .  sab(1,1+offs,1),sab(1,1+offs,2),ld0,
     .  giL,giL(1,1,2),ldRL,0d0,wk,wk(1+ld0*ldm1),ld0)
C     call yprm('Sj,j-1 gj-1,-1',2,wk,ld0*ldm1,ld0,ld0,ldm1)

C ... New gi,-1 = gij Sab initial-gi,-1.  gij has dim. (ldij,ldij)
      call yygemm('N','N',ld0,ldm1,ld0,1d0,gij,gij(1,1,2),ldij,
     .  wk,wk(1+ld0*ldm1),ld0,0d0,giL,giL(1,1,2),ldRL)
C     call yprm('gi,-1',2,giL,ldRL*ndRL,ldRL,ld0,ldm1)

C ... wk(ldR,ld0) = initial-g-1,i Sba, Sba = S(jpl-1,jpl)
C     Use for s(ipl-1,ipl) : sba(ld0,lds)+ = S(ipl,ipl-1)+
      call yygemm('N','C',ldm1,ld0,ldl,1d0,gLi,gLi(1,1,2),ndRL,
     .  sab(1,1+offs,1),sab(1,1+offs,2),ld0,
     .  0d0,wk,wk(1+ldm1*ld0),ldm1)
C     call yprm('gl,j-1 sj-1,j',2,wk,ldm1*ld0,ldm1,ldm1,ld0)

C ... New g-1,i = initial-g-1,i Sba gji.  gji has dim. (ldji,ldji)
      call yygemm('N','N',ldm1,ld0,ld0,1d0,wk,wk(1+ldm1*ld0),ldm1,
     .  gji,gji(1,1,2),ldji,0d0,gLi,gLi(1,1,2),ndRL)
C     call yprm('g-1,i',2,gLi,ndRL*ldRL,ndRL,ldm1,ld0)

      end
