      subroutine pgflu(s_ctrl,s_lat,s_ham,s_pot,s_str,s_site,isp,izp,zp,qp,
     .  semsh,lsprse,lbiglu,pgplp,lpgf,lgii,gii,strRL,gend,gij)
C- Layer Green's function by LU decomposition of banded matrix
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl nbas npl lpgf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  bandw ldham lncol neula offH lham lgen3 eula
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula lham neula
Cio    Passed to:  sblham plham plhamso plhamnc rotheu
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  palp
Co     Stored:     *
Co     Allocated:  pfnc
Cio    Elts passed:pfnc palp
Cio    Passed to:  sblham plham plhamnc
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax s nds npr
Cio    Passed to:  sblham plham plhamnc
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfra
Cio    Passed to:  plham plhamso plhamnc
Ci Inputs
Ci   isp   :current spin channel (1 or 2)
Ci   izp   :current energy index
Ci   zp    :complex energy : zp(izp) is current energy
Ci   qp    :k-point
Ci   semsh :parameters defining energy contour (emesh.f)
Ci   lsprse:not used.  For now, lsprse must be unity
Ci   lbiglu:include  -1 and npl layers in LU decomposition
Ci         :If lbiglu is set, pgflu generates xtl GF g(-1) and g(npl)
Ci         :What GF are computed depend on lpgf; see gii in Outputs
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   lpgf  :1 compute diagonal GF for all layers
Ci         :5 surface GF for layers -1, npl and off-diagonal
Ci         :  g(-1,npl) and g(npl,-1)
Ci   strRL :(lpgf=5) holds information about g(-1,npl) and g(npl,-1)
Ci         :strRL(1) leading dimension of gRL, second dim. of gRL
Ci         :strRL(2) leading dimension of gLR, second dim. of gLR
Co Outputs
Co   lgii  :attributes of the gii generated; see lgupac for conventions
Co    gii  :digaonal GF by layer: g_ii = gii(i)
Co         :Memory is allocated by pgflu.
Co         :What is generated depends on lpgf:
Co         :lpgf=1 : gii(0..npl-1) contain g for those PL.
Co         :         If lbiglu=0, gii(-1),gii(npl) are returned as gsLL and gsRR
Co         :         If lbiglu=1, gii(-1),gii(npl) are computed as gLL and gRR
Co         :lpgf=5 : the left- and right- surface Green's functions
Co                   are returned in gii(-1),gii(npl)
Co   gend  :some GF at endpoints or connecting endpoints
Co         :Which of the following are made depends on the mode
Co         :gend(1) = g(npl,-1)  Returned if lpgf=5
Co         :gend(2) = g(-1,npl)  Returned if lpgf=5
Co         :gend(3) = gsLL (left surface GF at -1) returned if gend(3)%gll is associated
Co         :gend(4) = gsRR (right surface GF at npl) returned if gend(4)%gll is associated
Co         :Note: gend(3)%gll and gend(4)%gll should be passed as NULL if you
Co         :      don't want them populated.
Co   gij   :(lpgf=8) off-diagonal GF connecting layer i to layers -1 and npl
Co         :gij(1,i) = gij(-1,i)
Co         :gij(2,i) = gij(i,-1)
Co         :gij(3,i) = gij(npl,i)
Co         :gij(4,i) = gij(i,npl)
Cl Local variables
Cl   ldisk :0 keep everything in memory
Cl          1 save gii to disk, reload after releasing h
Cl   ld0mx :maximum dimension of any PL
Cl   ldimx :sum of dimensions of all PL (hamiltonian dimensions)
Cl    ld0x :dimension of current PL
Cl   ldm1x :dimension of leftmost PL
Cl    ofb  :hamiltonian offset for current PL (layer for which g is calculated)
Cl    of1b :of1b+1 = first row for which g = (LU)^-1 rhs is calculated for for one PL
Cl    ldpl :dimension of hamiltonian between leftmost and current PL
Cl         :last row for which g = (LU)^-1 rhs is calculated for for one PL
Cl    b    :rhs to LU backsubstitution, portions of which get copied to g
Cl    pivc :pivc(kpl) = offset to the start of layer kpl in hamiltonian and thus:
Cl         :pivc(kpl) = offset to first element in layer in vector b for layer kpl
Cl         :            In noncollinear case, dimension is doubled.
Cl         :pivc(kpl+1)-pivc(kpl) is hamiltonian dimension of layer kpl
Cl    ndg  :source (column) dimension of GF for current PL, same as pgflp(3)
Cl    ld0  :field (row) matrix dimension for this PL, same as pgflp(4)
Cr Remarks
Cr    The hamiltonian has the block structure
Cr
Cr          (P0-S00    -S01                              )
Cr          (  -S10  P1-S11    -S12                      )
Cr          (          -S21  P2-S22    -S23              )
Cr          (                  -S32  P3-S33    -S34  ... )
Cr
Cr
Cr    In the noncollinear case the ii block is structured
Cr
Cr                      ( P++ - S     P+-     )
Cr           (P-S)_ii = (                     )
Cr                      ( P-+         P-- - S )_ii
Cr
Cr    Before exiting, g is rotated to the spinor representation which
Cr    makes P diagonal.
Cr
Cb Bugs
Cb   gii should be returned in kcplx=2 mode?
Cb   mode 8 has not been throughly checked!
Cu Updates
Cu   15 Jun 17 Noncollinear branch now works in global axis.
Cu             Rotation to local axis moved to pgfasa, after g's have been made
Cu             This renders LU and embedding modes consistent.
Cu   29 Apr 17 pgflu rotates both spin and orbital components to local quantization axis
Cu   06 Jun 16 Leads can be noncollinear
Cu      May 16 Can return off-diagonal GF connecting leads
Cu   24 Nov 15 Returns left surface gf in glr(3)%gll, if that pointer is associated
Cu                and right surface gf in glr(4)%gll, if that pointer is associated
Cu   02 Sep 15 complete migration to f90 pointers
Cu   08 Aug 13 replaced f77 pointers with f90 ones except passed pointers
Cu   10 Nov 11 Begin migration to f90 structures
Cu   25 Apr 05 Self-energy can have different Im Z than active region
Cu             for contours close to real axis.
Cu   03 Nov 03 lpgf=5 works for noncollinear branch
Cu   26 Feb 02 pgflu now works for lbiglu=1 for both lpgf=1 and lpgf=5
Cu   12 Dec 01 pgflu can compute off-diagonal GF for
Cu             transmission mode (lpgf=5)
Cu   28 Apr 00 made noncollinear
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lsprse,lbiglu,isp,izp,lpgf,pgplp(6,-1:*),lgii(-2:*),strRL(2)
      double precision zp(2,izp),qp(3),semsh(10)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: pivc(:),piv(:)
      complex(8), allocatable, target :: s(:)
      complex(8), allocatable,target :: sgs(:)
      complex(8), allocatable :: sll(:)
      complex(8), allocatable :: h(:)
      complex(8), allocatable :: b(:,:)
      complex(8), pointer :: wk(:),pfa(:,:)
      complex(8), pointer :: g(:,:)
C     double complex, allocatable:: gii(:,:,:)
C     double precision, allocatable:: gii(:,:,:)
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*),gij(4,-1:*),gend(4)
C ... Local parameters
      logical bittst,lso
      integer, parameter :: nkap0=4,n0H=5, PRTG=80,ISH=1,MZP=64,LGAM=128
      integer i,i1,i2,idim,ipl,ipr,itr,iwk(10),kcplx,kl,kpl,ku,lLR,lRL,
     .  lbloch,ld0,ld0mx,ld0x,ldh,ldham(16),ldim,ldimx,ldisk,ldkx,ldl,ldm1,ldm1x,
     .  ldn,ldnx,ldpl,ldr,lds,lhdimp,lidimp,lncol,lpack,lspc,modegg,mxorb,nbas,
     .  nbaslu,nbasp,ndg,ndgx,neul,nfilet,nl,npadl,npadr,npl,nspc,of1b,ofb,offL,
     .  offR,offa,offi,offp,ofg
      equivalence (lidimp,ldham(2)),(lhdimp,ldham(3))
      double precision xx
      double precision plat(3,3)
C     character*40 outs
      procedure(integer) :: isw,isum,fopn,fopnT,nglob,pgdim,rotspd

C --- Setup ---
      call tcn('pgflu')
C     call pshpr(70)
      call getpr(ipr)
      plat = s_lat%plat
      nl = s_ctrl%nl
      nbas = s_ctrl%nbas
      npl = s_ctrl%npl
      kl = s_ham%bandw
C     offH = s_ham%ooffH
      ldham = s_ham%ldham
      lncol = s_ham%lncol
      neul = s_ham%neula
      nspc = min(1+lncol,2)
      call sanrg(.true.,lsprse,1,1,'pgflu:','lsprse')
      npadl = pgplp(1,0)
      npadr = pgplp(1,npl-1) - pgplp(1,npl-2)
      nbasp = nbas + npadl + npadr
      mxorb = nglob('mxorb')
      lspc = isw(s_ctrl%lpgf(2) >= 10 .and. nspc == 2)

      ldisk = 0
      if (lpgf >= 5) ldisk = 0
      modegg = 0
      if (nspc == 2) modegg = 1
      if (ldisk /= 0) call rx('update pgflu for ldisk>0')

C --- Allocate gLL, gRR ---
      call pggmal(modegg,-1,-1,pgplp,-1,gii)
      call pggmal(modegg,npl,npl,pgplp,npl,gii)
      if (ldisk == 0 .and. lpgf /= 5) then
        call pggmal(modegg,0,npl-1,pgplp,0,gii)
      else
        nfilet = fopnT('TMP',-1,4,0)
        rewind nfilet
      endif

C     lbiglu: larger dimensions
      if (lbiglu /= 0) then
        ldim = s_ham%offH(nkap0*n0H,nbasp)
        nbaslu = nbasp
      else
        ldim = s_ham%offH(nkap0*n0H,nbas)
        nbaslu = nbas
      endif

      ldimx = nspc*ldim
      idim  = 0
      ku = kl
      ldh = 2*kl+ku+1
      lso = bittst(lncol,4)

C --- Generate noncollinear potential functions ---
      if (lso) then  ! plhamso does rotation internally; not needed here
      elseif (lncol /= 0) then
        call ptr_pot(s_pot,8+1,'pfnc',lhdimp,nspc**2,[0d0])
        pfa => s_pot%palp
C       nsp = 2
C       allocate(pfa(lhdimp*nsp))
C       call gf1kps(nsp,ldim,lhdimp,isw(.not.ldiwav),s_pot%palp,pfa)
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1,1,1d0)
        offp = lhdimp
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1+2*offp,1+2*3*lhdimp,1d0)
        i = rotspd(0)           ! local -> global z
        call rotspn(230010+100*i,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,xx,
     .    xx,lhdimp,lidimp,lidimp,lhdimp,1,xx,s_pot%pfnc)
C       call zprm('pfnc',2,s_pot%pfnc,lhdimp,lhdimp,4)
      endif

C --- Assemble hamiltonian in packed storage ---
      offa = 0
      lpack = 1  ! Packed storage
      kcplx = 1  ! complex*16 arithmetic
      if (lpack == 0) ldh = ldimx
      allocate(h(ldh*ldimx)); call dpzero(h,2*ldh*ldimx)

C ... Assemble hamiltonian in one fell swoop
      if (lncol == 0) then
        lbloch = lpack + 10*kcplx
C        if (bittst(s_ctrl%lham,LGAM)) then
C          pfa => s_pot%palp
C        else
C          pfa => s_pot%pf
C        endif
        call sblham(s_ham,s_pot,s_str,1,nbaslu,2*nbasp-nbas,
     .    lbloch,plat,ldim,idim,ldh,idim,ldim,kl,qp,h,xx,xx)

C ... Assemble hamiltonian PL by PL
C     Global spin quantization axis in the noncollinear case
      else
        call dpzero(h,ldh*ldimx*2)
        i1 = 0
        i2 = npl-1
        if (lbiglu /= 0) then
          i1 = -1
          i2 = npl
        endif
        do  ipl = i1, i2
          ldl = pgplp(4,max(ipl-1,-1))
          ld0 = pgplp(4,ipl)
          ldr = pgplp(4,min(ipl+1,npl))
          lds = ldl+ld0+ldr
          allocate(s(lds*nspc*ld0*nspc))
          call plham(1,s_ham,s_pot,s_str,s_site,plat,isp,kcplx,ipl,pgplp,qp,lds,ld0,s)
C         call pgflu2(0,1,nspc,ld0,ldl,ldr,lds,kl,offa,.false.,.true.,ld0*nspc,s,h)
          call pgflu2(lpack,kcplx,nspc,ld0,ldl,ldr,lds,kl,offa,isw(ipl>i1),isw(ipl<i2),ldh,s,h)
          deallocate(s)
          offa = offa + ld0*nspc
        enddo
      endif

C      allocate(s(ldimx*ldimx)); call dpzero(s,2*ldimx*ldimx)
C      call zpack(0,lpack,0,1,ldimx,1,ldimx,kl,ku,ldh,ldimx,h,s)
C      print *,'z=',zp(1,izp),zp(2,izp)
C      call yprm('after h is assembled h no bc(unp)',kcplx+2,s,ldimx**2,ldimx,ldimx,ldimx)
C     Alternatively write a binary file
C      i = fopna('pms',-1,4)
C      call ywrm(1,' ',kcplx+2,i,'(5f15.9)',s,ldimx**2,ldimx,ldimx,ldimx)
C      deallocate(s)
C      call fclose(i)
C      print *, 'wrote binary pms'
C      call rx('done')

C --- Add boundary condition connecting to left SIF GF ---
      ld0 = pgplp(4,0)
      ld0x = nspc*ld0
      allocate(sgs(ld0x**2))

C ... Read left surface g into sgs, possibly store into gii(-1) and gend(3), kcplx=0
      ipl = -1
      call rdsurfg(ipl,1,isp,lspc,nspc,izp,npl,0,semsh,pgplp,zp,qp,plat,lgii,sgs)
      call dcopy(2*size(sgs),sgs,1,gii(ipl)%gll,1)
C     if (ipr >= PRTG/1) call yprm('gsLL',2,gii(ipl)%gll,ld0x**2,ld0x,ld0x,ld0x)
C     Copy left surface gLL to gend(3)%gll
      if (associated(gend(3)%gll)) then
        if (size(gend(3)%gll) /= size(sgs)) call rx('pgflu: size mismatch')
        call dcopy(2*size(sgs),sgs,1,gend(3)%gll,1)
      endif

C ... Make Sigma(L)
      allocate(sll(3*ld0x**2),wk(3*ld0x**2))
      if (nspc == 2) then
        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ipl,npl,ld0,3*ld0,pgplp,qp,wk,sll)
      else
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,ipl,pgplp,qp,3*ld0,ld0,sll)
      endif
      call pgsif2(100,sll,ld0x,sgs)
C     if (ipr >= PRTG/1) call zprm('sigma(L)',2,sgs,ld0x,ld0x,ld0x)
      deallocate(sll,wk)

C ... Add sigma(L) to h
!      print *, 'omit sigma(L)'
      call pgflu1(lpack,kcplx,ld0*nspc,kl,ldh,ld0*nspc,0,0,0,sgs,sgs,h,h)
C      allocate(s(ldimx*ldimx)); call dpzero(s,2*ldimx*ldimx)
C      call zpack(0,lpack,0,1,ldimx,1,ldimx,kl,ku,ldh,ldimx,h,s)
C      print *,'z=',zp(1,izp),zp(2,izp)
C      call yprm('after h+lbc is assembled',kcplx+2,s,ldimx**2,ldimx,ldimx,ldimx)
      deallocate(sgs)

C --- Add boundary condition connecting to right SIF GF ---
      ld0 = pgplp(4,npl-1)
      ld0x = nspc*ld0
      allocate(sgs(ld0x**2))

C ... Read right surface g into sgs, possibly store into gii(npl) and gend(4)
      ipl = npl
      call rdsurfg(ipl,1,isp,lspc,nspc,izp,npl,0,semsh,pgplp,zp,qp,plat,lgii,sgs)
      call dcopy(2*size(sgs),sgs,1,gii(ipl)%gll,1)
C     if (ipr >= PRTG/1) call yprm('gsRR',2,gii(ipl)%gll,ld0x**2,ld0x,ld0x,ld0x)
C     Copy right surface gRR to gend(4)%gll
      if (associated(gend(4)%gll)) then
        if (size(gend(4)%gll) /= size(sgs)) call rx('pgflu: size mismatch')
        call dcopy(2*size(sgs),sgs,1,gend(4)%gll,1)
      endif

C ... Make Sigma(R)
      allocate(sll(3*ld0x**2),wk(3*ld0x**2))
      if (nspc == 2) then
        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ipl,npl,ld0,3*ld0,pgplp,qp,wk,sll)
      else
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,ipl,pgplp,qp,3*ld0,ld0,sll)
      endif
      call pgsif2(101,sll,ld0x,sgs)
C     if (ipr >= PRTG/1) call zprm('sigma(R)',2,sgs,ld0x,ld0x,ld0x)
      deallocate(sll,wk)

C     Add sigma(R) to h
      offa = nspc*isum(npl-1,pgplp(4,0),6)
      if (lbiglu /= 0) offa = nspc*isum(npl+1,pgplp(4,-1),6)
!     print *, 'omit sigma(R)'
      call pgflu1(lpack,kcplx,ld0*nspc,kl,ldh,ld0*nspc,0,offa,offa,sgs,sgs,h,h)
      deallocate(sgs)

C  ... Printout
C      if (ipr >= PRTG/100 .and. izp == 17) then
CC      This generates an unpacked version from a packed one
C       allocate(s(ldimx*ldimx)); call dpzero(s,2*ldimx*ldimx)
C       call zpack(0,lpack,0,1,ldimx,1,ldimx,kl,ku,ldh,ldimx,h,s)
C       call yprm('h+bc(unp)',kcplx+2,s,ldimx**2,ldimx,ldimx,ldimx)
C       stop
C       deallocate(s)
CC      call yprm('pgflu: h+bc',kcplx+2,h,ldh*ldimx,ldh,ldh,ldimx)
C      endif

C --- Table of hamiltonian offsets ---
C     pivc is used for pivot constraints, and to find off-diagonal g in b
C     pivc(kpl) = offset to the start of layer kpl in hamiltonian
      i1 = 0
      i2 = npl-1
      if (lbiglu /= 0) then
        i1 = -1
        i2 = npl
      endif
      allocate(pivc(i1:i2+1))
      pivc(i1) = 0
      do  ipl = i1, i2
        call pghoff(ipl,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
        pivc(ipl+1) = nspc*(offi+pgplp(4,ipl))
      enddo

C --- Initial LU decomposition ---
      call tcn('LU-decomp')
      allocate(piv(ldimx))
      if (lpack == 0) then
        call zgetrf(ldimx,ldimx,h,ldh,piv,i)
      else
C       call zgbtrf(ldimx,ldimx,kl,ku,h,ldh,piv,i)
        if (pivc(i2+1) /= ldimx) call rx('pgflu: bug in pivc')
        call z1btrf(ldimx,ldimx,kl,ku,h,ldh,piv,pivc(i1+1),i)
      endif
      call tcx('LU-decomp')
      call rxx(i/=0,'pgflu : failed to decompose h')

C --- For each layer, diagonal g's by back substitution ---
C     For lpgf=5, only off-diagonal GF needed.
C     For lpgf=8, g(npl,all), g(all,npl), g(-1,all), g(all,-1) are calculated
C     Also if lpgf=5 or lpgf=8, gii(-1) and gii(npl) are left untouched
      ld0mx = nspc*pgdim(3,-1,npl,-1,npl,pgplp)  ! Largest PL dimension
C     b=rhs to backsubstitution, performed layer by layer to preserve
C     O(N) scaling.  dimensioned b(ldimx,nspc*(max PL dimension))
      allocate(b(ldimx,ld0mx)); call dpzero(b,2*ldimx*ld0mx)
      ldpl = ldim*nspc
      itr = 0
      do  ipl = i2, i1, -1
        ld0  = pgplp(4,ipl)
        ld0x = nspc*ld0
        ldl  = pgplp(4,max(ipl-1,-1))
        ldr  = pgplp(4,min(ipl+1,npl))
        lds  = ldl+ld0+ldr
        ndg  = pgplp(3,ipl)
        ndgx = nspc*ndg
        ofg  = pgplp(6,ipl)
C       ibR  = pgplp(1,ipl-1)
        call pghoff(ipl,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
C       Offset to diagonal part of hamiltonian for this PL
        ofb = nspc*offi
C       Offset to diagonal part of hamiltonian for prior PL
C       This guarantees g(ipl-1,ipl) and g(ipl,ipl) are both made.
        of1b = nspc*min(max(offL,0),offi)
C       Transmission mode requires g_LR; spin torque modes require g_Li
        if (lpgf == 5 .or. lpgf == 8) of1b = 0

C   ... Set up boundary conditions for this PL
        if (ipl < i2 .and. lpack /= 0) then
          allocate(s(lds*nspc*ld0x))
          call plham(1,s_ham,s_pot,s_str,s_site,plat,isp,kcplx,ipl,pgplp,qp,lds,ld0,s)
          if (nspc == 2) then
            allocate(wk(ld0x**2))
          else
            wk => s
          endif
!          print *, '!! remove bc in call to pgflu5'
!          call pgflu5(.false.,kcplx,itr,ofb,ld0,ldl,ldr,nspc,s,wk,ldimx,b)
          call pgflu5(lpgf>0,kcplx,itr,ofb,ld0,ldl,ldr,nspc,s,wk,ldimx,b)
          if (nspc == 2) deallocate(wk)
          deallocate(s)
        endif

C   ... Make b = gji connecting layers (jpl..ipl,ipl) to ipl by back-substitution
C       jpl (specified through of1b) is typically -1 or ipl-1
        call pgflu3(lpack,kcplx,itr,ldpl,ofb,of1b,ld0x,kl,ku,ldh,h,piv,ldimx,b)

C   ... Transmission mode: make g(-1,npl)  and  g(npl,-1)
C   ... Spin torque mode:  make g(ipl,npl),g(npl,ipl),g(-1,ipl),g(ipl,-1)
        if (lpgf == 5 .or. lpgf == 8) then

          ldm1 = pgplp(4,-1)
          ldn  = pgplp(4,npl)
          ldm1x = nspc*ldm1
          ldnx = nspc*ldn

C         call zprm('b = g(all,npl)',2,b,ldimx,ldimx,ldnx)

C         Save g(-1,npl)
          lRL = strRL(1)
          lLR = strRL(2)
          call zmcpy('N',b,ldimx,1,gend(2)%gll,lLR,1,ldm1x,ldnx)
C          if (ipr >= PRTG/10) then
C            call zprm('g_-1,npl',2,gend(2)%gll,lLR,ldm1x,ldnx)
C          endif

C         Save g(kpl,npl), kpl = -1..npl
          if (lpgf == 8) then
          do  kpl = -1, npl
            ldkx = pivc(kpl+1)-pivc(kpl)
            iwk(1:2) = shape(gij(4,kpl)%gll)
            if (iwk(1)/=ldkx .or. iwk(2)/=ldnx) call rx('pgflu: dimensioning error')
            call zmcpy('N',b(1+pivc(kpl),1),ldimx,1,gij(4,kpl)%gll,ldkx,1,ldkx,ldnx)
            call info5(PRTG,0,0,' g_%i,%i(1,1): %;12,7D',kpl,i2,b(1+pivc(kpl),1),0,0)
C            if (ipr >= PRTG/1) then
C              call awrit2(' g_%i,%i:',outs,len(outs),0,kpl,i2)
C              call zprm(outs,2,gij(4,kpl)%gll,ldkx,ldkx,ldnx)
C            endif
          enddo
          endif

C         Create b = g+(npl,all).  of1b=0 => Elements b(1:ldpl,1:ldnx) are generated
          call dpzero(b, ldimx*ld0mx*2)
          call pgflu3(lpack,kcplx,1,ldpl,ofb,of1b,ldnx,kl,ku,ldh,h,piv,ldimx,b)
C         call zprm('b(all,npl) = g+(npl,all)',2,b,ldimx,ldimx,ldnx)

C         Save g(npl,-1)
          call zmcpy('C',b,ldimx,1,gend(1)%gll,1,lRL,ldnx,ldm1x)
C          if (ipr >= PRTG/1) then
C           call zprm('g_npl,-1',2,gend(1)%gll,lrl,ldnx,ldm1x)
C          endif

C     ... All off-diagonal g's connecting to layer -1 or npl
          if (lpgf == 8) then
C         Save g(npl,kpl), kpl = -1..npl
          do  kpl = -1, npl
            ldkx = pivc(kpl+1)-pivc(kpl)
            iwk(1:2) = shape(gij(3,kpl)%gll)
            if (iwk(1)/=ldnx .or. iwk(2)/=ldkx) call rx('pgflu: dimensioning error')
            call zmcpy('C',b(1+pivc(kpl),1),ldimx,1,gij(3,kpl)%gll,1,ldnx,ldnx,ldkx)
            call info5(PRTG,0,0,'g_%i,%i(1,1): %;12,7D',i2,kpl,b(1+pivc(kpl),1),0,0)
C            if (ipr >= PRTG/10) then
C              call awrit2('g_%i,%i:',outs,len(outs),0,i2,kpl)
C              call zprm(outs,2,gij(3,kpl)%gll,ldnx,ldnx,ldkx)
C            endif
          enddo

C         Create b = g(all,-1).   Elements b(1:ldpl,1:ldnx) are generated
          call dpzero(b, ldimx*ldm1x*2)
          call pgflu3(lpack,kcplx,0,ldpl,0,0,ldm1x,kl,ku,ldh,h,piv,ldimx,b)
C         call zprm('b = g(all,-1)',2,b,ldimx,ldimx,ldm1x)
          do  kpl = -1, npl
            ldkx = pivc(kpl+1)-pivc(kpl)
            iwk(1:2) = shape(gij(4,kpl)%gll)
            if (iwk(1)/=ldkx .or. iwk(2)/=ldm1x) call rx('pgflu: dimensioning error')
            call zmcpy('N',b(1+pivc(kpl),1),ldimx,1,gij(4,kpl)%gll,ldkx,1,ldkx,ldm1x)
            call info5(PRTG,0,0,' g_%i,%i(1,1): %;12,7D',kpl,i1,b(1+pivc(kpl),1),0,0)
C            if (ipr >= PRTG/10) then
C              call awrit2('g_%i,%i:',outs,len(outs),0,kpl,i1)
C              call zprm(outs,2,gij(4,kpl)%gll,ldkx,ldkx,ldm1x)
C            endif
          enddo

C         Create b = g+(-1,all).  of1b=0 => Elements b(1:ldpl,1:ldm1x) are generated
          call dpzero(b, ldimx*ldm1x*2)
          call pgflu3(lpack,kcplx,1,ldpl,0,0,ldm1x,kl,ku,ldh,h,piv,ldimx,b)
C         call zprm('b = g+(-1,all)',2,b,ldimx,ldimx,ldm1x)
          do  kpl = -1, npl
            ldkx = pivc(kpl+1)-pivc(kpl)
            iwk(1:2) = shape(gij(3,kpl)%gll)
            if (iwk(1)/=ldm1x .or. iwk(2)/=ldkx) call rx('pgflu: dimensioning error')
            call zmcpy('C',b(1+pivc(kpl),1),ldimx,1,gij(3,kpl)%gll,1,ldm1x,ldm1x,ldkx)
            call info5(PRTG,0,0,'g_%i,%i(1,1): %;12,7D',i1,kpl,b(1+pivc(kpl),1),0,0)

C            if (ipr >= PRTG/10) then
C              call awrit2('g_%i,%i:',outs,len(outs),0,i1,kpl)
C              call zprm(outs,2,gij(3,kpl)%gll,ldm1x,ldm1x,ldkx)
C              call yprmi('g_%i,%i:',i1,kpl,3,gij(3,kpl)%gll,0,ldm1x,ldm1x,ldkx)
C            endif
          enddo

          endif ! lpgf == 8

          call ztoyy(gend(1)%gll,lRL,lLR,lRL,lLR,1,0)
          call ztoyy(gend(2)%gll,lLR,lRL,lLR,lRL,1,0)

          exit   ! We have g's connecting L and R, which is all mode 5 or 8 needs
        endif

C   ... Poke into in gii, else create temporary array
        g => gii(ipl)%gll
        if (associated(gii(ipl)%gll)) then
          g => gii(ipl)%gll
        else
          allocate(g(ld0x,ndgx))
        endif
        if (lpgf<0) call lgupac(lgii(ipl),'p',1,1,0,0,0,0)
        if (lpgf>0 .or. ipl >= i2) call lgupac(lgii(ipl),'p',1,3,0,0,0,0)

        call pgflu4(kcplx,itr,ofb,of1b,ld0x,ndgx,ofg,ldimx,b,g) ! g in kcplx=0

C        if (izp == 17) then
C          call yprmi('g_%i,%-1j%i global spin axis',ipl,ipl,2+0,g,ld0x*ndgx,ld0x,ld0x,ndgx)
C        endif

C   ... Rotate noncollinear (P-S)^-1 to local spin quantization axis
C       This has been moved to pgfasa
C        if (nspc == 2) then
C          if (ndg /= ld0) call rx('implement spinor rotation for off-diagonal g')
CC         Shift indxsh to starting at 1 for ipl
C          offa = offi
C          call pblch3(0,1,offa,1,nbasp,mxorb,ldham,s_ham%iprmb)
CC         call pblch2(-offa,1,nbasp,mxorb,s_ham%iprmb)
C          optrot = 100*kcplx*0 + rotspd(1-0) ! Backwards rotation
C          if (bittst(s_ham%lham,256)) optrot = optrot + 10 ! Spherical harmonics
C          if (mod(s_lat%lsym/8,2) /= 0) optrot = optrot + 1000 ! Suppress orbital rotation
C          call gtibpl(ipl,npl,pgplp,ib1,ib2)
C          call rotheu(optrot,s_ham,nl,ib1,ib2,nbasp,0,ld0,ld0,ndg,g)
CC          ncsw = 30000 + 10*kcplx*0  g is already noncollinear
CC          call rotspn(ncsw,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,
CC     .      xx,xx,0,ld0,ndg,ld0,ndg,xx,g)
CC         Restore indxsh
C          call pblch3(11,1,offa,1,nbasp,mxorb,ldham,s_ham%iprmb)
CC         call pblch2(offa,1,nbasp,mxorb,s_ham%iprmb)
C
CC         call priprm(' ham offsets for current PL',
CC     .      1,nbasp,mxorb,ldham,ldham,ldham,s_ham%iprmb)
C
CC          call awrit1(' g_%i,%-1j%i local spin axis',outs,len(outs),0,ipl)
CC          call yprm(outs,2+0,g,ld0x*ndgx,ld0x,ld0x,ndgx)
C
C!       if (ipl == 6) then
C!        if (ipr >= PRTG/10) then
C!          call yprmi('g_%i,%-1j%i local spin axis',ipl,ipl,2+0,g,ld0x*ndgx,ld0x,ld0x,ndgx)
C!        endif
C
C        endif

        if (ldisk /= 0) then
          call ywrm(1,' ',2,nfilet,' ',g,ld0x*ndgx,ld0x,ld0x,ndgx)
        endif

        ldpl = ldpl - ld0x

      enddo
      deallocate(b,piv,pivc,h)

C     OLD: When gii(-1),gii(npl) = surface GF, put in kcplx=0 mode.
C     Notw read from rdsurfg, already in kcplx=0 mode
C      if (lbiglu==0 .or. lpgf==5 .or. lpgf==8) then
C        if (lgii(-1) > 0) call ztoyy(gii(-1)%gll,ld0x,ld0x,ld0x,ld0x,1,0)
C        if (lgii(npl) > 0) call ztoyy(gii(npl)%gll,ld0x,ld0x,ld0x,ld0x,1,0)
C      endif

C      if (ldisk /= 0) then
C        call rx('ldisk branch needs checking')
C        rewind nfilet
C        call pggmal(modegg,0,npl-1,pgplp,0,gii)
C        if (nspc == 2) call rx('pgflu: check disk branch')
C        do  ipl = npl-1, 0, -1
C          ld0 = pgplp(4,ipl)
C          ndg = pgplp(3,ipl)
C          call rdm(nfilet,32,ld0*ndg*2,' ',gii(ipl)%gll,ld0,ndg)
C          if (ipr >= PRTG) then
C            call awrit1(' g_%i:',outs,len(outs),0,ipl)
C            call yprm(outs,2+0,gii(ipl)%gll,ld0*ndg,ld0,ld0,ndg)
C          endif
C        enddo
C      endif

      call tcx('pgflu')
      return

   99 call rx('pgflu: failed to read bc file')

      end

      subroutine pgflu1(lpack,kcplx,ld0,kl,ldh,ncol,offcs,offa,offb,splc,splr,hplc,hplr)
C- Kernel called by pgflu, to add subblocks of PL to larger array
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack :0,h is stored in unpacked form
Ci         :1,h is stored in banded form (see Remarks)
Ci   kcplx :1: spl,h in complex*16 format
Ci          2: spl,h have real, imaginary separated by columns
Ci   ld0   :dimension of current layer
Ci   kl    :number of sub- and super-diagonals within the band of h
Ci         :Used only if lpack=1
Ci   ldh   :leading dimension of h
Ci   ncol  :number columns to copy.  This is:
Ci          ld0 when copying (PL,PL) block
Ci          ldl when copying (PL,PL-1) block
Ci          ldr when copying (PL,PL+1) block
Ci          where ldl and ldr = dimension of PL-1 and PL+1 hamiltonians
Ci   offcs :column offset to (PL,PL') block of spl.  This is:
Ci          0       when copying (PL,PL) block
Ci          ld0     when copying (PL,PL-1) block
Ci          ld0+ldl when copying (PL,PL+1) block
Ci          For noncollinear second spin channel, add ld0+ldl+ldr
Ci          (Row offset is always zero)
Ci   offa  :row offset in h
Ci         :For noncollinear second spin channel, add ld0
Ci   offb  :column offset in h.  This is:
Ci          offa          when copying (PL,PL) block
Ci          offa-ldl*nspc when copying (PL,PL-1) block
Ci          offa+ld0*nspc when copying (PL,PL+1) block
Ci   splc  :strux for this layer, and coupling to (left,right)
Ci          in complex*16 format (case kcplx=1)
Ci   splr  :strux for this layer, and coupling to (left,right)
Ci          with real, imaginary separated by columns (case kcplx=2)
Co Outputs
Co   hplc  :hamiltonian accumulated for this layer
Co          in complex*16 representation  (case kcplx=1)
Co   hplr  :hamiltonian accumulated for this layer
Co          with real, imaginary sep. by columns.  (case kcplx=2)
Cr Remarks
Cr   pgflu1 pokes one block of principal layer structure constants (or
Cr   hamiltonian) into the appropriate place in full hamiltonian array.
Cr   Three calls are necessary to copy all blocks connecting a PL to
Cr   its neighbors; see ncol,offcs,offb.
Cr
Cr   hpl is stored in normal (lpack=0) or band (lpack=1) format.
Cr   The band format follows LAPACK band storage conventions:
Cr      h(i,j) is stored in location (kl+ku+1+i-j,j)
Cr   We assume here ku=kl
Cr
Cr   For complex*16 representation, splr and hplr are not used.
Cr   For separated real,imag representation, splc and hplc are not used.
Cr   Routine ztoy converts one representation to another.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer ld0,offa,ldh,lpack,kcplx,kl,ncol,offcs,offb
      double precision splc(2,ld0,*),hplc(2,ldh,ncol+offb)
      double precision splr(ld0,2,*),hplr(ldh,2,ldh)
C ... Local parameters
      integer i,j,ofc

      call tcn('pgflu1')

      if (kcplx == 1) then
        call zpack(1,0,lpack,1+offa,ld0+offa,1+offb,ncol+offb,kl,kl,
     .    ld0,ldh,splc(1,1-offa,1+offcs-offb),hplc)
      else
        if (lpack == 0) then
C         Unpacked storage case
          call ymscop(1,ld0,ncol,ld0*2,ldh*2,0,offcs,offa,offb,
     .      splr,ld0,hplr,ldh)
        elseif (lpack == 1) then
C         Packed storage case
          ofc = kl+kl+1
          do  12  j = 1, ncol
          do  12  i = 1, ld0
          hplr(i-j+ofc,1,j+offb) = hplr(i-j+ofc,1,j+offb) +
     .                             splr(i,1,j+offcs)
   12     hplr(i-j+ofc,2,j+offb) = hplr(i-j+ofc,2,j+offb) +
     .                             splr(i,2,j+offcs)
        else
          call rxi('pgflu1: bad argument lpack,',lpack)
        endif
      endif

      call tcx('pgflu1')
      end

      subroutine pgflu2(lpack,kcplx,nspc,ld0,ldl,ldr,lds,kl,offa,lcpl,
     .  lcpr,ldh,spl,h)
C- Assembles 1 PL of hamiltonian for LU decomposition
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack :0,h is unpacked
Ci         :1,h is stored in packed storage (see Remarks
Ci   kcplx :1: spl,h in complex*16 format
Ci          2: spl,h have real, imaginary separated by columns
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ld0   :dimension of current layer
Ci   ldl   :dimension of left adjacent layer
Ci   ldr   :dimension of right adjacent layer
Ci   lds   :second dimension of s
Ci   kl    :number of sub- and super-diagonals within the band of h
Ci   offa  :offset to hamiltonian for this PL
Ci   lcpl  :0 omit coupling to left layer h(i,i-1)
Ci         :1 copy h(i,i-1) in the normal hamiltonian order; see Remarks
Ci         :2 copy h(i,i-1) to the right of h(i,i); see Remarks
Ci   lcpr  :0 omit coupling to right layer h(i,i+1)
Ci         :1 copy h(i,i+1) in the normal hamiltonian order; see Remarks
Ci         :2 copy h(i,i+1) to the left of h(i,i); see Remarks
Ci   ldh   :leading dimension of h
Ci   spl   :strux for this layer, and coupling to (left,right)
Ci          in complex*16 format, when kcplx=1 or
Ci          with real, imaginary separated by columns when kcplx=2
Co Outputs
Co    h    :hamiltonian accumulated for this layer
Co          in complex*16 format, when kcplx=1 or
Co          with real, imaginary separated by columns when kcplx=2
Cr Remarks
Cr   Let  j=layer i-1, k=layer i+1.  Let ldim = ld0+ldl+ldr
Cr   * Collinear case.  spl is ordered as:
Cr
Cr                0      ld0   ld0+ldl
Cr                |      |       |
Cr       0  ->    \/     \/      \/
Cr                ii     ij     ik
Cr
Cr   The three blocks are copied to h in the normal hamiltonian order, which looks like:
Cr           offa-2*ldl  offa  offa+ld0
Cr                |       |       |
Cr     offa ->    \/      \/      \/
Cr                ij      ii      ik
Cr
Cr   Special mode: lcpl=2 or lcpr=2.  ij and ik are exchanged:
Cr           offa-2*ldr  offa  offa+ld0
Cr                |       |       |
Cr     offa ->    \/      \/      \/
Cr                ik      ii      ij
Cr
Cr   * Noncollinear case. spl for layer i looks as follows, generated
Cr     from plhamso or plham2
Cr
Cr                                      ldim+  ldim+  ldim+
Cr                0      ld0   ld0+ldl    0    ld0   ld0+ldl
Cr                |      |       |        |       |      |
Cr       0  ->   \/     \/      \/       \/      \/     \/
Cr                ii++   ij++   ik++      ii+-    ij++   ik+-
Cr                ii-+   ij-+   ik-+      ii--    ij--   ik--
Cr
Cr   spl is copied to h in the normal hamiltonian order, which looks like:
Cr
Cr               offa-   offa-                   offa+  offa+
Cr               2*ldl   ldl    offa  offa+ld0   2*ld0  2*ld0+ldr
Cr                |      |       |     |         |      |
Cr      offa ->   \/     \/      \/    \/        \/     \/
Cr                ij++   ij+-    ii++  ii+-      ik++   ik+-
Cr                ij-+   ij--    ii-+  ii--      ik-+   ik--
Cr
Cr   Special mode: lcpl=2 or lcpr=2.  ij and ik are ordered like spl, but
Cr   spin components are grouped together for a particular (ij) PL pair.
Cr   pgembs and pgsif require this ordering.
Cr
Cr                                offa+  offa+       offa+   offa+
Cr               offa  offa+ld0   2*ld0  2*ld0+ldl   2*ld2   2*ld2+ldr   where ld2=ld0+ldl
Cr                |     |         |      |            |      |
Cr     offa ->    \/    \/        \/     \/           \/     \/
Cr                ii++  ii+-      ij++   ij+-         ik++   ik+-
Cr                ii-+  ii--      ij-+   ij--         ik-+   ik--
Cr
Cr   All rows are copied in one block.
Cr   Columns are copied in six blocks.
Cu Updates
Cu   27 Sep 15 lcpl and lcpr made into integers with extra functionality (lcpl=2)
Cu   27 Apr 00 Adapted to handle the noncollinear case.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lcpl,lcpr
      integer ld0,ldl,ldr,lds,offa,ldh,lpack,kcplx,kl,nspc
C     this is for complex*16 format
      double precision spl(2,ld0,nspc,lds,nspc),h(2,ldh,*)
C ... Local parameters
      integer isp,offcs,offci,offb,ld2

      call sanrg(.true.,kcplx,1,2,'pgflu:','kcplx')

      do  isp = 1, nspc

C   ... Coupling (ipl,ipl)
        offci = (isp-1)*(ld0+ldl+ldr)
        offcs = offci
        offb  = offa + (isp-1)*ld0
        call pgflu1(lpack,kcplx,ld0*nspc,kl,ldh,ld0,offcs,offa,offb,spl,spl,h,h)

C   ... Coupling (ipl,ipl-1)
        if (lcpl /= 0) then
          offcs = offci + ld0
          offb  = offa - ldl*(nspc-isp+1)
          if (lcpl == 2 .or. lcpr == 2) offb  = offa + nspc*ld0 + (isp-1)*ldl
          call pgflu1(lpack,kcplx,ld0*nspc,kl,ldh,ldl,offcs,offa,offb,spl,spl,h,h)
        endif

C   ... Coupling (ipl,ipl+1)
        if (lcpr /= 0) then
          offcs = offci + ld0+ldl
          offb  = offa + nspc*ld0 + (isp-1)*ldr
C         if (lcpl == 2 .or. lcpr == 2) offb  = offa - ldr*(nspc-isp+1)
          if (lcpl == 2 .or. lcpr == 2) then
            ld2 = ld0+ldl
            offb  = offa + nspc*ld2 + (isp-1)*ldr
          endif
          call pgflu1(lpack,kcplx,ld0*nspc,kl,ldh,ldr,offcs,offa,offb,spl,spl,h,h)
        endif

      enddo

C     Only for the unpacked case, kcplx=1
C     call zprm('h',2,h,ldh,ld0*nspc,lds*nspc)

      end

      subroutine pgflu3(lpack,kcplx,itr,ldpl,ofb,of1b,ld0x,kl,ku,ldh,h,ipiv,ldb,b)
C- Layer Green's function for 1 PL by back-substitution of LU decomp'sn
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack :0,h is stored in unpacked form
Ci         :1,h is stored in banded form (see Remarks)
Ci   kcplx :1: spl,h in complex*16 format
Ci          2: spl,h have real, imaginary separated by columns
Ci   itr   :0 calculate g solving LU g = b
Ci          1 calculate g solving (LU)+ g+ = b
Ci          The former produces g_ki with k = j..i (j<=i)
Ci          The latter produces (g_ik)+ with k = j..i (j<=i)
Ci          i is current layer index spec'd by ofb (see below)
Ci          j is is index specified through of1b (see below)
Ci          See description for b, below
Ci   ldpl  :dimension of hamiltonian for PL 0...current PL
Ci         :last row for which b is made.  ldpl is normally ofb+ld0x.
Ci   ofb   :hamiltonian offset for this PL, called i in these notes
Ci         :-1 is added to b(i+ofb,i) before back-substitution
Ci   of1b  :hamiltonian offset for PL j: g_j:i,i are returned in b
Ci         :of1b is starting row index for layer j
Ci   ld0x  :dimension of current PL, and number of columns of matrix b
Ci   kl,ku :number of band sub- and super-diagonals
Ci   ldh   :leading dimension of band hamiltonian h
Ci   h     :H-z or S-P, decomposed by LU decomposition
Ci   ipiv  :work array holding pivoting information
Ci   ldb   :leading dimension of b
Cio Inputs/Outputs
Cio  b     :On input, b = 0 or if adding correction from layer k,
Cio         b(ofb+1:ofb+ld0x,1:ld0x) = -(gik Ski)+ if itr=0
Cio         b(ofb+1:ofb+ld0x,1:ld0x) = - Sik gki   if itr=1
Cio         where i is PL of interest, and k is right adjacent.
Cio         (See Remarks to pgflu5).
Cio        :On output, b holds subblock of g for PL's j..i connected to i.
Cio         Typically of1b corresponds to offset of PL j with j<=i.
Cio         b(1+of1b:ldpl,1:ld0x) =  g_j..i,i   (itr = 0)
Cio         b(1+of1b:ldpl,1:ld0x) = (g_i,j..i)+ (itr = 1)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lpack,kcplx,itr,ldpl,ldb,ofb,of1b,kl,ku,ldh,ipiv(ldb),ld0x
      double precision h(2,ldh,ldb),b(2,ldb,ld0x)
C ... Local parameters
      integer i,PRTG,iprint
      parameter (PRTG=100)
      character*1 trans,ch*2

      call tcn('LU back-subst')

C     Only complex*16 implemented for LU decomposition
      call sanrg(.true.,kcplx,1,1,'pgflu3:','kcplx')
      call sanrg(.true.,lpack,0,1,'pgflu3:','lpack')

      if (lpack == 0) call dpzero(b,2*ldb*ld0x)
C     Not needed if use z1btrs -> z1btrs
C     do  5  j = 1, ld0x
C     do  5  i = 1, ofb
C       b(1,i,j) = 0
C       b(2,i,j) = 0
C   5 continue
C     endif

      do  i = 1, ld0x
        b(1,i+ofb,i) = b(1,i+ofb,i) - 1
      enddo

      if (iprint() >= PRTG) call yprm('su rhs',2+kcplx,b,ldb*ld0x,ldb,ldpl,ld0x)
      trans = 'N'
      if (mod(itr,10) == 1) trans = 'C'

      if (lpack == 0) then
C       zgetrs makes the entire b(:,1:ld0x), not just a subblock
        call zgetrs(trans,ldb,ld0x,h,ldh,ipiv,b,ldb,i)
        if (iprint() >= PRTG/1)
     .  call yprm('g from zgetrs',2+kcplx,b,ldb*ld0x,ldb,ldb,ld0x)
      else
C       z1btrs makes the subblock b(1+of1b:ldpl,1:ld0x)
C       call zgbtrs(trans,ldpl,kl,ku,ld0x,h,ldh,ipiv,b,ldb,i)
        call z1btrs(trans,ldpl,kl,ku,ld0x,h,ldh,ipiv,b,ldb,1+of1b,i)
        if (iprint() >= PRTG/1) then
         ch = 'g'
         if (trans == 'C') ch = 'g+'
         call yprm(ch//' from z1btrs',2+kcplx,b,ldb*ld0x,ldb,ldpl,ld0x)
        endif
      endif
      call rxx(i /= 0,'pgflu : failed to decompose h')

      call tcx('LU back-subst')
      end

      subroutine pgflu4(kcplx,itr,ofb,of1b,ld0x,ndgx,ofg,ldb,b,gii)
C- Copy g as calculated from LU decomposition to gii
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx :1: b in complex*16 format
Ci          2: b have real, imaginary separated by columns
Ci   itr   :0 b holds g_ki with k<=i
Ci          1 b holds (g_ik)+ with k<=i.
Ci   ofb   :hamiltonian offset for this PL.
Ci         :b(1+ofb,1) corresponds to the first element of the
Ci         :diagonal part of g
Ci   of1b  :offset to first row index for which b is calculated
Ci         :of1b<ofb means that b contains some (off-digaonal parts) of
Ci         :the GF connecting current layer to preceding layers.
Ci   ld0x  :leading dimension of gii
Ci   ndgx  :second dimension of gii
Ci   ofg   :offset in gii to diagonal part of gii
Ci         : ofg and ofb may differ because b is dimensioned for
Ci         : the entire GF connected to the current layer, but
Ci         : only a portion will be stored in gii.  Specifically:
Ci         : b = b(ldb,ndgx), with b(1,1) = first element of
Ci         : off-diagonal GF connecting layer 0 with current layer
Ci         : and b(1+ofb,1) = 1st element of diagonal g
Ci         : NB: Not all of b may be calculated; see of1b.
Ci         : g = g(ld0x,ndgx) with g(1+ofg,1) = 1st elemeent of
Ci         : diagonal g for current layer.
Ci         : Thus, g(1+ofg,1) = b(1+ofb,1)
Ci   ldb   :leading dimension of b and
Ci   b     :b holds subblock g connected to PL i:
Ci          (itr = 0) b(1+of1b:ldb,1:ld0x) = g_ki with k<=i
Ci          (itr = 1) b(1+of1b:ldb,1:ld0x) = (g_ik)+ with k<=i
Co Outputs
Co   gii   :gii(ld0x,ndgx) is calculated (in kcmplx=0 format for now)
Cb Bugs
Cb   For now, gii in kcmplx=0 format
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kcplx,itr,ldb,ofb,of1b,ld0x,ndgx,ofg
      double precision b(2,ldb,*),gii(ld0x,ndgx,2)
C ... Local parameters
      character*1 trans

      trans = 'N'
      if (mod(itr,10) == 1) trans = 'C'

      call sanrg(.true.,kcplx,1,1,'pgflu4:','kcplx')

      if (trans == 'N') then
        call sanrg(.true.,ofg,0,0,'pgflu4:','ofg')
        call ymcpy(b(1,1+ofb-ofg,1),2*ldb,2,1,gii,ld0x,1,ld0x*ndgx,ld0x,ndgx)
      else
        call sanrg(.true.,of1b,0,ofb-ofg,'pgflu4:','of1b')
        call ymcpy(b(1,1+ofb-ofg,1),2,2*ldb,1,gii,ld0x,1,ld0x*ndgx,ld0x,ndgx)
        call dscal(ld0x*ndgx,-1d0,gii(1,1,2),1)
      endif
      end

      subroutine pgflu5(lbc,kcplx,itr,ofb,ld0,ldl,ldr,nspc,spl,wk,ldb,b)
C- Setup rhs to generate GF for 1 PL by back-substitution of LU decomp'n
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbc   :F : Set rhs to 0
Ci         :T : Set rhs to (g_ab sba)+ or (Sab g_ba)
Ci   kcplx :1: spl,b in complex*16 format
Ci          2: spl,b have real, imaginary separated by columns
Ci   itr   :0 b holds g_ab
Ci          1 b holds (g_ba)+
Ci   ofb   :offset to diagonal part of hamiltonian for prior PL
Ci   ld0   :dimension of current layer
Ci   ldl   :dimension of left adjacent layer
Ci   ldr   :dimension of right adjacent layer
Ci   spl   :strux for this layer, and coupling to (left,right)
Ci   ldb   :leading dimension of b
Cio Inputs/Outputs
Cio  itr   :specifies whether b holds gab+ or gba
Cio        :Output itr specifies which form of back-substitution is
Cio        :needed to generate gaa for this PL; see Remarks
Cio        :itr is set to 1 if input itr = 0
Cio        :itr is set to 0 if input itr = 1
Cio  b     :On input, b(ofb+1:ofb+ld0,1:ldr) contains gab+ or gba,
Cio        :(depending on itr) connecting layer a to right adjacent.
Cio        :b is overwritten to create boundary conditions for
Cio        :calculation of gaa; see Remarks.
Cio        :b(ofb+1:ofb+ld0,1:ld0) = -(gab Sba)+ if input itr=0
Cio        :b(ofb+1:ofb+ld0,1:ld0) = - Sab gba   if input itr=1
Cr Remarks
Cr   Given LU decomposition LU = -(g_aa~)^-1, g_aa is obtained via
Cr   Dyson's equation.  Which form of Dyson's equation is used
Cr   depends on whether gab is or gba is input:
Cr     g_aa = (1 + g_ab Sba) g_aa~  = (LU)+^-1 (-1 - g_ab Sba)+
Cr          = g_aa~ (1 + Sab g_ba)  = (LU)^-1  (-1 - Sab g_ba)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lbc
      integer kcplx,itr,ldb,ofb,ld0,ldl,ldr,nspc
      double precision b(2,ldb,*),spl(2,ld0,nspc,ld0+ldl+ldr,nspc)
      double precision wk(2,ld0,nspc,ld0,nspc)
C ... Local parameters
      integer i,j,PRTG,iprint,ld0x
      parameter (PRTG=100)
      double precision alpha(2),beta(2)
      character*20 strn

      call tcn('pgflu5')
      ld0x = ld0*nspc

C     Only complex*16 implemented now
      call sanrg(.true.,kcplx,1,1,'pgflu5:','kcplx')

      alpha(1) = 1
      alpha(2) = 0
      beta(1) = 0
      beta(2) = 0

      if (lbc) then

C     spl <- (gab Sba)^+ = (gab Sab+)^+ = Sab gab+ (itr=0)
C         <- (Sab gba) (itr = 1)
      call zgemm('N','C',ld0,ld0x,ldr,alpha,spl(1,1,1,1+ld0+ldl,1),
     .  ld0x,b(1,1+ofb,1),ldb,beta,wk,ld0x)

C     NC case: 11+12 blocks with s11 b1 + s12 b2 = s11 b1 since s12=0
C     (zgemm above). Still need  s21 b1 + s22 b2 = s22 b2
      if (nspc == 2) then
        call zgemm('N','C',ld0,ld0x,ldr,alpha,spl(1,1,2,1+ld0+ldl,2),
     .    ld0x,b(1,1+ofb,1+ldr),ldb,beta,wk(1,1,2,1,1),ld0x)
      endif

      else
        call dpzero(wk,size(wk))
      endif

      do  j = 1, ld0x
      do  i = 1, ld0x
        b(1,i+ofb,j) = -wk(1,i,1,j,1)
        b(2,i+ofb,j) = -wk(2,i,1,j,1)
      enddo
      enddo
      if (iprint() >= PRTG) then
        strn = '-(gab Sba)^+'
        if (itr == 1) strn = '-(Sab gba)'
        call yprm(strn,kcplx+2,b(1,1+ofb,1),ldb*ldb,ldb,ld0x,ld0x)
      endif
      itr = 1 - itr

      call tcx('pgflu5')
      end

C      subroutine pgflu6(nspc,lRL,lLR,ldL,ldR,spl,spr,wk,gsLL,gsRR,
C     .  gRL,gLR)
CC- Make off-diagonal g_n-1, g_-1,n from g_n-1,0 and g_0,n-1
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
CCi         :(not implemented)
CCi   lRL   :leading dimension of gRL, second dim. of gRL
CCi   lLR   :leading dimension of gLR, second dim. of gLR
CCi   ldL   :hamiltonian dimension of layers 0 and -1
CCi   ldR   :hamiltonian dimension of layers npl-1 and npl
CCi   spl   :strux connecting layers npl-1 and npl (plham)
CCi   spr   :strux connecting layers 0 and -1 (plham)
Cci   wk    :work array of dimension 2*lLR*lRL
CCi   gsLL  :left surface GF for layer -1
CCi   gsRR  :right surface GF for layer npl
CCio Inputs/Outputs
CCi   gRL   :On input, crystal GF connecting layers 0 and npl-1
CCi         :On output crystal GF connecting layers -1 and npl
CCi   gLR   :On input, crystal GF connecting layers npl-1 and 0
CCi         :On output crystal GF connecting layers npl and -1
CCr Remarks
CCr   From Dyson's equation,
CCr     g(0,npl)  = g(0,npl-1) S(npl-1,npl) gsRR
CCr     g(-1,npl) = gsLL S(-1,0) g(0,npl)
CCr     g(npl,0)  = gsRR S(npl,npl-1) g(npl-1,0)
CCr     g(npl,-1) = g(npl,0) S(0,-1) gsLL
CCu Updates
CCu   12 Dec 01 First created
CC ----------------------------------------------------------------------
C      implicit none
C      integer nspc,ldL,ldR,lLR,lRL
C      double precision wk(2,lLR,lRL),
C     .  spl(2,ldL,nspc,ldL,3,nspc),gLR(2,lLR,lRL),gsLL(2,ldL,ldL),
C     .  spr(2,ldR,nspc,ldR,3,nspc),gRL(2,lRL,lLR),gsRR(2,ldR,ldR)
CC ... Local parameters
CC     integer ldLx,ldRx
C      double precision alp(2),beta(2)
C
C      call tcn('pgflu6')
CC     ldLx = ldL*nspc
CC     ldRx = ldR*nspc
C      if (nspc == 2) call rx('pgflu6 not ready for nspc=2')
C
CC      call yprm('gLR',2+1,gLR,0,nspc*lLR,ldL,ldR)
CC      call yprm('g+RL',2+1,gRL,0,nspc*lRL,ldR,ldL)
C
C      alp(1) = 1
C      alp(2) = 0
C      beta(1) = 0
C      beta(2) = 0
C
CC     spl(1) <- gsLL S(-1,0)
C      call zgemm('N','N',ldL,ldL,ldL,alp,gsLL,ldL,
C     .  spl(1,1,1,1,3,1),ldL,beta,spl,ldL)
CC     call yprm('gsLL S',2+1,spl,0,ldL,ldL,ldL)
C
CC     spr(1) <- S(npl-1,npl) gsRR
C      call zgemm('N','N',ldR,ldR,ldR,alp,spr(1,1,1,1,3,1),ldR,
C     .  gsRR,ldR,beta,spr,ldR)
CC     call yprm('S gsRR',2+1,spr,0,ldR,ldR,ldR)
C
CC     wk <- g(0,npl) = g(0,npl-1) S(npl-1,npl) gsRR
C      call zgemm('N','N',ldL,ldR,ldR,alp,gLR,lLR,spr,ldR,beta,wk,ldL)
CC     call yprm('g_0n',2+1,wk,0,ldL,ldL,ldR)
C
CC     gLR <- g(-1,npl) = gsLL S(-1,0) g(0,npl)
C      call zgemm('N','N',ldL,ldR,ldR,alp,spl,ldL,wk,ldL,beta,gLR,lLR)
CC     call yprm('g_-1n',2+1,gLR,0,lLR,ldL,ldR)
C
CC     spl(1) <- S(0,-1) gsLL
C      call zgemm('N','N',ldL,ldL,ldL,alp,spl(1,1,1,1,2,1),ldL,
C     .  gsLL,ldL,beta,spl,ldL)
CC      call yprm('S gsLL',2+1,spl,0,ldL,ldL,ldL)
C
CC     spr(1) <- gsRR S(npl,npl-1)
C      call zgemm('N','N',ldR,ldR,ldR,alp,gsRR,ldR,
C     .  spr(1,1,1,1,2,1),ldR,beta,spr,ldR)
CC      call yprm('gsRR S',2+1,spr,0,ldR,ldR,ldR)
C
CC     wk <- g(npl,0) = gsRR S(npl,npl-1) g(npl-1,0)
C      call zgemm('N','C',ldR,ldL,ldL,alp,spr,ldL,gRL,lRL,beta,wk,ldR)
CC     call yprm('g_n0',2+1,wk,0,ldR,ldR,ldL)
C
CC     gRL <- g(npl,-1) = g(npl,0) S(0,-1) gsLL
C      call zgemm('N','N',ldR,ldL,ldL,alp,wk,ldR,spl,ldL,beta,gRL,lRL)
CC     call yprm('g_n-1',2+1,gRL,0,lRL,ldR,ldL)
C
C      call tcx('pgflu6')
C      end
      subroutine pgflu7(nbas,npadl,npadr,nl,s_spec,ips,ipa,offH,iprmb,ldham)
C- Makes hamiltonian offsets, including padded layers in hamiltionaian
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makidx
Ci Inputs
Ci   nbas  :size of basis
Ci   npadl :number of padding layers on left
Ci   npadr :number of padding layers on right
Ci   nl    :(global maximum l) + 1
Ci   ips   :species table: site ib belongs to species ips(ib)
Co Outputs
Ci   ipa   :indices to basis, included paadded parts
Co   offH  :Offsets to hamiltonian matrix (makidx.f)
Co   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co   ldham :ldham (i) holds size of block i, plus all the preceding blocks.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,npadl,npadr,nl,ips(*),ipa(nbas+2*npadl+2*npadr+1),iprmb(*),ldham(*)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,nkaph,nbasp

      nbasp = nbas + npadl + npadr

C      do  ib = 1, nbas
C        ipa(ib) = ib+npadl
C      enddo
C      do  ib = 1, npadl
C        ipa(nbas+ib) = ib
C      enddo
      do  ib = 1, npadl
        ipa(ib) = nbas+ib
      enddo
      do  ib = 1, nbas
        ipa(ib+npadl) = ib
      enddo
      do  ib = nbas+npadl+1,  nbas+npadl+npadr+1
        ipa(ib) = ib
      enddo

      nkaph = 1
      call makidx(nl,nkaph,1,nbasp,ipa,s_spec,ips,offH,iprmb,ldham)

      end
