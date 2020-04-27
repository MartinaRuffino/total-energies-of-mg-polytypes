      subroutine pgemb(s_ham,s_pot,s_str,s_site,plat,isp,ipL,ip0,ipR,
     .  nrpt,vpl,mode,pgplp,qp,iprt,lgii,gii,strRL,gij)
C- Embeds a new layer between one or two semi-infinite Green's functions
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham lncol neula lgen3 lham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plhamnc plham plhamso pvemb pgbevl
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  palp
Co     Stored:     *
Co     Allocated:  pfnc
Cio    Elts passed:pfnc palp
Cio    Passed to:  plhamnc plham pvemb pgbevl
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plhamnc plham pvemb pgbevl
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plhamnc plham plhamso pvemb pgbevl
Ci Inputs:
Ci   plat:   primitive 2D lattice translation vectors
Ci   isp:    current spin channel (1 or 2)
Ci   ipL:    layer to left of current PL, ip0-1 mode=1,3,5)
Ci           or current PL (mode=4)
Ci   ip0:    current PL to embed
Ci   ipR:    layer to right of current PL (mode=2,3,4,6)
Ci           or ipl+nrpt-1 (mode=10, nrpt>1)
Ci           or ipL (mode=5)
Ci   nrpt:   number of layers for which repeats
Ci   vpl:    used in adiabatic approximation; see pvemb.
Ci   mode:   specifies type of GF to be created.  See Remarks
Ci          1s digit signifies:
Ci           0 do not calculate any GF
Ci           1 (nrpt=0)
Ci             gii(ip0) <- left surface GF created from gii~(ipL)
Ci             i.e.  gii(ip0) = g~(ip0,ip0)
Ci             (nrpt>1)
Ci             off-diagonal GF also made ??
Ci             gii(ip0+nrpt-1) <- g~(ip0+nrpt-1,ip0+nrpt-1)
Ci             nrpt>0 requires rvec to be tacked to the end of gii(ipR)%gll,
Ci             rvec contained in gii(ip0)%gll(1,ld0+1)
Ci           2 Analogous to mode 1, except surface gii~(ip0) is created
Ci             from right gii~(ipR), and is right surface GF
Ci           3 gii(ip0) <- crystal GF created by embedding both
Ci             left and right terms from gii(ipL)~,gii(ipR)~.
Ci           6 gii(ip0) <- surface gii~(ip0,ip0) created from
Ci             gii~(ipL,ipL) as in mode 1.  Also make from the off-diagonal
Ci             gii~(-1,ip0),gii~(ip0,-1) from gii~(-1,ipL),gii~(ipL,-1)
Ci           7 Last layer, transport mode.
Ci             Make gii(-1,ip0),gii(ip0,-1) corresponding to the true xtal GF.
Ci             Requires:
Ci             gii~(-1,ipL),gii~(ipL,-1),gii~(ipR,ipR).
Ci             Also gRR~ is saved in g(npl,npl), instead of gRR; see Remarks
Ci          16 like 6, but weakly varying potential
Ci   pgplp:  table of PL-dependent quantities (see pgfset)
Ci   qp:     k-point (NB: projection of qp onto PL direction not used).
Ci   iprt    printout verbosity
Ci   strRL (mode=6,7) strux containing info about RL GF:
Ci         (1) = ldRL, (2) = ndRL, (3) = offRL, (4) = offLR
Cl Local and global variables
Cl   npl   : number of principal layers for layer programs (pgfset.f)
Cio Inputs/Outputs:
Cio  gii   : Vector of layer diagonal GF : gij(layer k) = gii(k)%gll
Cio        : On output, gii(ip0), and possibly some off-diagonal gf
Cio        : connecting adjacent layers to ip0.
Cio        : In modes 1 and 6, gii~  is stored in gii(ip0)
Cio        : If nrpt>0, gii(ipR) is stored in gii(ipR)
Cio  gij   : Holds certain off diagonal G's
Cio        :     mode      gij(1)       gij(2)   gij(3)   gij(4)
Cio        :       1         -            -        -        -
Cio        :       3         -            -        -        -
Cio        :       6  In  g~(ipL,-1)  g~(-1,ipL)   -        -
Cio        :       6  Out g~(ip0,-1)  g~(-1,ip0)   -        -
Cio  lgii  : attributes of gii generated; see lgupac for conventions
Cr Remarks
Cr  *Embedding proceeds layer-by-layer, in a left-to-right pass
Cr   in which a succession of surface GF are created, followed
Cr   a right-to-left pass in which the surface GF are converted into
Cr   crystal GF.  pgemb handles on step in the left-to-right pass.
Cr   The diagonal surface GF at layer ip0, g~(ip0,ip0) is built up from
Cr   the diagonal surface GF at layer ipL.  Thus,
Cr     g~(-1,-1) -> g~(0,0), -> g~(1,1) -> ... -> g~(ip0,ip0)
Cr   In modes 1,6, or 7 pgemb creates g~(ip0,ip0) from g~(ipL,ipL)
Cr   (the differences between modes 1,6,7 are described below).
Cr
Cr   In the special case that the right surface GF to the right of ip0
Cr   is also known, the proper crystal GF at layer ip0 can be created
Cr   by adding embedding terms from both left and right surface GF:
Cr     g~(ipL,ipL) + g~(ipR,ipR) -> g(ip0,ip0)
Cr   This is done in mode 3.
Cr
Cr   The right-to-left passes in which the surface GF are converted into
Cr   crystal GF are also accomplished layer-by-layer.  See pgdysn for
Cr   a description.
Cr
Cr   Some applications require other GF besides the diagonal ones.
Cr   In particular, the off-diagonal g(-1,npl), g(npl,-1) are useful
Cr   for computing transmission probabilities.  Mode 6 is similar
Cr   to mode 1, except that in addition to creating g~(ip0,ip0),
Cr   the off-diagonal g~(-1,ip0), g~(ip0,-1) are created from
Cr   the g~(-1,ipL), g~(ipL,-1) (created from a prior call to pgemb).
Cr
Cr   Mode 7 is like mode 6, except that now the right GF g~(ipR,ipR)
Cr   must be available, and the proper g(-1,npl), g(npl,-1) are created.
Cr   Since the diagonal GF are not used, the proper crystal GF for
Cr   diagonal blocks, e.g. g(npl,npl) are not saved in gii(npl); rather
Cr   we save gsR = gRR~ in gii(npl), as needed by pgcurr.
Cr
Cr   To summarize the modes:
Cr   (mode 1,6) left semi-infinite GF, embedding g~(ipL,ipL)
Cr              mode 6 also creates off-diagonal g~(-1,ip0), g~(ip0,-1)
Cr   (mode 2)   right semi-infinite GF,  embedding g(ipR,ipR)
Cr   (mode 3)   crystal GF, embedding g(ipL,ipL) and g(ipR,ipR)
Cr   (mode 7)   like mode 6, except the off-diagonal g(-1,ip0) and
Cr              g(ip0,-1) correspond to the proper crystal GF.
Cb Bugs
Cb   Partitioning into lower and higher blocks not implemented.
Cu Updates
Cu   02 Dec 15 modes 5,6,7 extended to noncollinear case
Cu    3 Oct 15 mode 1 extended to the noncollinear case
Cu   02 Sep 15 complete migration to f90 pointers
Cu   08 Jul 13 replaced f77 pointers with f90 ones except passed pointers
Cu   10 Nov 11 Begin migration to f90 structures
Cu    3 Jan 03 bug fix (K. Belashchenko)
Cu   20 Feb 02 bug fix for lpgf=5, last layer.
Cu   13 Dec 01 small changes to accomodate calc. of transmission
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,ipL,ip0,ipR,nrpt,mode,iprt,pgplp(6,-1:10),lgii(-2:10),strRL(2),nl2
      double precision plat(3,3),qp(3),vpl(-1:*)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C     Pointer arrays to different layer GF's
      type (s_lgf) :: gii(-2:*),gij(4)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: s00(:,:,:,:)
      complex(8), allocatable :: gnn(:,:)
      complex(8), allocatable :: wk(:,:,:,:)
      complex(8), pointer :: gLL(:,:),gRR(:,:),pfa(:,:)
C ... Local parameters
      logical, parameter :: sanity=.true.
      integer, parameter :: PRTG=80/1
      logical bittst,lso
      integer lds,ld0,ldl,ldr,mode0,mode1,lgnn,lgbii,lncol,lhdimp,lidimp,i,
     .  npl,nglob,ndg0,ndgl,ndgr,ofl,ofr,of0,offp,nspc,neul,nbasp,nl,mxorb
      integer ld0x,ldlx,ldrx,ldsx,ldm1x,ndg0x,ndgrx,ndglx
      double precision xv(1)
      character outs*80
      procedure(integer) :: rotspd
      procedure(real(8)) :: dlength

C --- Setup ---
      call tcn('pgemb')
      npl = nglob('npl')
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
      lhdimp = s_ham%ldham(3)
      lidimp = s_ham%ldham(2)
C      npadl = pgplp(1,0)
C      npadr = pgplp(1,npl-1) - pgplp(1,npl-2)
C      nbasp = nbas + npadl + npadr
      nl = nglob('nl')
      nl2 = nl*nl
      nbasp = nglob('nbasp')
      neul = s_ham%neula
      mxorb = nglob('mxorb')

C      if (bittst(s_ham%lncol,4)) then ! SO is turned on
C        if (dlength(size(s_ham%eula),s_ham%eula,1) /= 0)
C     .    call rx('SO not yet implemented with noncollinear in Dyson mode')
C      endif

      call getpr(iprt)
      if (ip0 < -2 .or. ip0 > npl+1) call fexit2(-1,111,
     .  ' Exit -1 PGEMB: expected ipl to be >=0 and <%i but found ipl=%i',npl,ip0)
      ldl = pgplp(4,max(ip0-1,-1))
      ld0 = pgplp(4,max(ip0,-1))
      ldr = pgplp(4,min(ip0+1,npl))
      lds = ldl+ld0+ldr
      ld0x = ld0*nspc
      ldlx = ldl*nspc
      ldrx = ldr*nspc
      ldsx = lds*nspc
      ldm1x = pgplp(4,-1)*nspc
      ndg0x = ndg0*nspc
      ndgRx = ndgR*nspc
      ndgLx = ndgL*nspc
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode0 == 0) return

C ... Generate noncollinear potential functions.  Rotate P; S generated in global axis
      lso = bittst(lncol,4)
      if (lso) then ! plhamso does rotation of P internally
      elseif (lncol /= 0) then
        call ptr_pot(s_pot,8+1,'pfnc',lhdimp,nspc**2,xv)
        pfa => s_pot%palp
C       nsp = 2
C       allocate(pfa(lhdimp*nsp))
C       call gf1kps(nsp,ldim,lhdimp,isw(.not.ldiwav),s_pot%palp,pfa)
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1,1,1d0)
        offp = lhdimp
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1+2*offp,1+2*3*lhdimp,1d0)
        i = rotspd(0)           ! local -> global z
        call rotspn(230010+100*i,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xv,xv,
     .    xv,lhdimp,lidimp,lidimp,lhdimp,1,xv,s_pot%pfnc)
C       call zprm('pfnc',2,s_pot%pfnc,lhdimp,lhdimp,4)
      endif

C --- Strux connecting ip0 to adjacent layers, global spin quantization axis ---
      allocate(s00(ld0,nspc,lds,nspc),wk(ld0,nspc,lds,nspc))
      call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ip0,npl,ld0,lds,pgplp,qp,wk,s00)
C      if (lncol /= 0) then
C        kcplxl = 1              ! plham requires kcplx=1
C        call plham(1,s_ham,s_pot,s_str,s_site,plat,isp,kcplxl,max(ip0,-1),pgplp,qp,lds,ld0,wk)
C        call awrit1(' P-S(%i) from plham:',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,wk,ld0x,ld0x,ldsx)
C
CC       Debugging: 1,2,3 blocks spin 1
CC       mch out.copt -split g 1,36+1,72+1 1,36+1,72+1,108+1,144+1,180+1,216+1 g11
CC       mch out.copt -split g 1,36+1,72+1 1,36+1,72+1,108+1,144+1,180+1,216+1 g12
CC       mch out.copt -split g 1,36+1,72+1 1,36+1,72+1,108+1,144+1,180+1,216+1 g13
C
CC   ... Rotate P-S by Euler angles
CC       Convert raw hamiltonian offsets offi,offL,offR into layer form
CC       offi -> 0, offL -> ld0, offR -> ld0+ldl
C        call pghoff(ip0,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
C        call gtibpl(ip0,npl,pgplp,ib01,ib02)
C        call gtibpl(max(ip0-1,-1),npl,pgplp,ibL1,ibL2)
C        call gtibpl(min(ip0+1,npl),npl,pgplp,ibR1,ibR2)
CC       Modify iprmb by adding largest index to exclude all pairs but restore PL, PL-1, PL+1
C        offi  = offi + lhdimp
C        offL  = offL - ld0 + lhdimp
C        offR  = offR - ld0 - ldl + lhdimp
C        allocate(iprmb(size(s_ham%iprmb))); call icopy(size(iprmb),s_ham%iprmb,1,iprmb,1)
C        call pblch2(lhdimp,1,nbasp,nl2,iprmb)
C        call pblch2(-offi,ib01,ib02,nl2,iprmb)
C        call pblch2(-offL,ibL1,ibL2,nl2,iprmb)
C        call pblch2(-offR,ibR1,ibR2,nl2,iprmb)
CC       Rotate the three blocks
C        ncsw = 30000 + 10*kcplxl + 100*rotspd(1-1)
C        call rotspn(ncsw,1,nbasp,nbasp,nl,iprmb,s_ham%eula,neul,xv,xv,xv,ld0,ld0,lds,ld0,lds,wk,wk)
C        deallocate(iprmb)
C        call awrit1(' rot P-S(%i):',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,wk,ld0x,ld0x,ldsx)
C
CC   ... Rearrange (spin,orbital) indexing
C        call dpzero(s00,2*size(s00))
C        call pgflu2(0,kcplxl,nspc,ld0,ldl,ldr,lds,0,0,2,2,ld0*nspc,wk,s00)
C        call awrit1(' spin reordered P-S(%i):',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,s00,ld0x,ld0x,ldsx)
C
CC       Convert to kcplx=0 format
C        call ztoyy(s00,ld0x,ldsx,ld0x,ldsx,kcplxl,0)
C      else
C        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,max(ip0,-1),pgplp,qp,lds,ld0,s00)
C      endif

      if (iprt >= PRTG/1) then
        call awrit1(' P-S(%i):',outs,len(outs),0,max(ip0,-1))
        call yprm(outs,2+0,s00,ld0x*ldsx,ld0x,ld0x,ldsx)
      endif

C --- Make some GF at ip0, depending on mode and what is available ---
C ... Embed nrpt layers at ip0..ip0+nrpt-1 from left SIF
      if (mode0 == 1 .and. nrpt > 1) then
        if (nspc == 2) call rx('pgemb not ready for nspc == 2')
C       Eigenvectors of bulk GF, and:
C       gii(ip0) = gf(ip0) for repeating pot at ip0,ip0+1,ip0+2,...
        deallocate(wk)
        if (pgplp(3,ip0) /= ld0) call rx('repeating branch not ready')
        call pvemb(s_ham,s_pot,s_str,s_site,sanity,plat,isp,nspc,1,vpl,ip0,npl,
     .    pgplp,qp,ndg0,ndgl,ndgr,nrpt,ofL,of0,ofR,lgbii,lds,ldl,ldr,
     .    lgii,ipL,ipR,s00,gLL,gii(ip0)%gll,gRR,ld0,strRL,gij)

C ... Embed nrpt layers at ip0..ip0-nrpt+1 from right SIF
      elseif (mode0 <= 3 .and. nrpt > 1) then
        if (nspc == 2) call rx('pgemb not ready for nspc == 2')
        call rxi('pgemb not ready for nrpt gt 1 and mode',mode)

C ... Embed layer at ip0 from left and/or right SIF
C     mode0 = 1 -> gii~(ip0) from gii~(ipL)
C     mode0 = 2 -> gii~(ip0) from gii~(ipR)
C     mode0 = 3 -> gii(ip0) from gii~(ipL) and gii~(ipR)
      elseif (mode0 <= 3) then
        call pgembs(ldsx,ld0x,ldlx,ldrx,ndg0x,ndgLx,ndgrx,of0,ofL,ofR,s00,
     .    10*mode0,lgii(ipL),gLL,lgii(ipR),gRR,wk,lgii(ip0),gii(ip0)%gll)

C ... Embed nrpt layers at ip0..ip0+nrpt-1 from left SIF
C     and also make off-diagonal g(-1,ip0+nrpt-1)
      elseif (mode0 == 6 .and. nrpt > 1) then
        if (nspc == 2) call rx('pgemb not ready for nspc == 2')
C   ... Make gs(jpl,jpl),gs(ipl,jpl),gs(jpl,ipl),gs(-1,jpl),gs(jpl,-1)
        deallocate(wk)
        i = 3 + 10*mode1
        call pvemb(s_ham,s_pot,s_str,s_site,sanity,plat,isp,nspc,i,vpl,ip0,npl,
     .    pgplp,qp,ndg0,ndgl,ndgr,nrpt,ofL,of0,ofR,lgbii,lds,ldl,ldr,
     .    lgii,ipL,ipR,s00,gLL,gii(ip0)%gll,gRR,ld0,strRL,gij)

C ... Embed layer at ip0 from left SIF
C     and also make off-diagonal g(-1,ip0), g(ip0,-1)
      elseif (mode0 == 6) then
        if (lds < strRL(1)) then
          deallocate(wk)
          allocate(wk(ld0,nspc,strRL(1),nspc))
        endif
C       Make  gii~(ip0)
        call pgembs(ldsx,ld0x,ldlx,ldrx,ndg0x,ndglx,ndgrx,of0,ofL,ofR,s00,
     .    10,lgii(ipL),gLL,lgii(ipR),gRR,wk,lgii(ip0),gii(ip0)%gll)
C       Make  gij(1)<-gi-1~, gij(2)<-g-1i~ where i = ip0
        call pgofLR(sanity,ld0x,ldlx,ldrx,ldm1x,s00,wk,lgii(ip0),ld0x,
     .    gii(ip0)%gll,ld0x,gii(ip0)%gll,strRL(1),strRL(2),gij(1)%gll,gij(2)%gll)

      elseif (mode0 == 7 .and. nrpt > 1) then
        if (nspc == 2) call rx('pgemb not ready for nspc == 2')
        call rx('pgemb : cannot combine mode 7 with nrpt>1')

C ... Last step, transport mode: make g(-1,npl), g(npl,-1)
C     In this (transport) mode, ipL=npl-1, ip0=npl, ipR=npl
C     We make g(npl,npl) first, from which we make g(-1,npl), g(npl,-1)
C     Instead of saving g(npl,npl) in gii(npl) as would be customary
C     we save gRR~ (see Remarks).  Since we need to make g(npl,npl)
C     temporary space is allocated for it.
      elseif (mode0 == 7) then
C       Make xtal g(npl,npl) in a temporary array and discard it.
        allocate(gnn(ld0x,ld0x))
C       2 Jan 03 bug fix: next line added (Kirill Belashchenko)
        lgnn = 1
C       Make g(npl,npl)
        call pgembs(ldsx,ld0x,ldlx,ldrx,ndg0x,ndglx,ndgrx,of0,ofL,ofR,s00,
     .    30,lgii(ipL),gLL,lgii(ipR),gRR,wk,lgnn,gnn)
C       Make g(npl,-1) and g(-1,npl)
        call pgofLR(sanity,ld0x,ldlx,ldrx,ldm1x,s00,wk,lgnn,ld0x,
     .    gnn,ld0x,gnn,strRL(1),strRL(2),gij(1)%gll,gij(2)%gll)
        deallocate(gnn)
      else
        call rx('pgemb: bad mode')
      endif
      if (allocated(wk)) deallocate(wk)

C --- (debugging) print gii ---
      if (iprt >= PRTG/1) then
      if (mode0 /= 6 .and. mode0 /= 7) then
        call awrit2('embedded g for layer %i, mode=%i',outs,len(outs),0,ip0,mode)
        call yprm(outs,2,gii(ip0)%gll,ld0x*ndg0x,ld0x,ld0x,ndg0x)
C        call ymscop(0,ld0,ld0,ld0,ld0,0,of0,0,0,
C     .    gii(ip0)%gll,ld0*ndg0,s00,ld0*ld0)
C        call awrit2('embedded g_ii for layer %i, mode=%i',outs,
C     .    len(outs),0,ip0,mode)
C        call yprm(outs,2,s00,ld0*ld0,ld0,ld0,ld0)
      else
        call awrit1('g_%i,-1',outs,len(outs),0,ipR)
        call yprm(outs,2,gij(1)%gll,strRL(1)*strRL(2),strRL(1),ld0x,pgplp(4,-1)*nspc)
        call awrit1('g_-1,%i',outs,len(outs),0,ipR)
        call yprm(outs,2,gij(2)%gll,strRL(2)*strRL(1),strRL(2),pgplp(4,-1)*nspc,ld0x)
      endif
      if (mode0 == 1 .and. nrpt > 1) then
        call awrit1('gii for layer %i',outs,len(outs),0,ip0+nrpt-1)
        call yprm(outs,2,gRR(1,of0),ld0x**2,ld0x,ld0x,ld0x)
      endif
      if (mode0 == 4 .and. nrpt > 1) then
        call awrit1('gii for layer %i',outs,len(outs),0,ip0-nrpt+1)
        call yprm(outs,2,gLL(1,of0),ld0*ld0x**2,ld0x,ld0x,ld0x)
      endif
      endif

      deallocate(s00)
      call tcx('pgemb')
      end
      subroutine pgembs(lds,ld0,ldl,ldr,ndg0,ndgl,ndgr,of0,ofL,ofR,
     .  sii,mode,lgLL,gLL,lgRR,gRR,wk,lgii,gii)
C- Make surface g_ii from semiinfinite g_LL and/or semiinfinite g_RR
C ----------------------------------------------------------------
Ci Inputs:
Ci   sii:   strx connecting layer i with L and R
Ci   ld0,lds       dimensions sii and wk
Ci   ld0,ndg0,of0  diagonal dim, tot. dim. and diag. offset to gii
Ci   ldl,ndgl,ofl  diagonal dim, tot. dim. and diag. offset to gLL
Ci   ldr,ndgr,ofr  diagonal dim, tot. dim. and diag. offset to gRR
Ci   mode 10s digit
Ci          1 Create gii~ from gll~  and (P-S)ii
Ci          2 Create gii~ from grr~  and (P-S)ii
Ci          3 Create gii  from gll~ and grr~(ipl+1) and (P-S)ii
Ci         1s digit
Ci          1 create gii^-1 instead of gii
Ci   lgLL  : attributes of gLL generated; see lgupac for conventions
Ci   gLL   :left surface gf, needed if 10s digit mode = 1 or 3
Ci   lgRR  :attributes of gRR generated; see lgupac for conventions
Ci   gRR   :right surface gf, needed if 10s digit mode = 2 or 3
Ci   wk    :a work array of the same dimension as s...
Ci          maybe doesn't need to be that big?
Co Outputs:
Co   gii,lgii: some kind diagonal GF at layer i, and information
Co          describing what g contains (see pgemb)
Cr Remarks
Cr   sii   is OVERWRITTEN (diagonal part only)
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lds,ld0,ldl,ldr,mode,lgLL,lgRR,lgii,ndg0,ndgl,ndgr,
     .  of0,ofL,ofR
      double precision gii(ld0,ndg0,2),gLL(ldl,ndgl,2),gRR(ldr,ndgr,2),
     .  sii(ld0,lds,2),wk(ld0,lds,2)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wkr(:)
      real(8), allocatable :: gRRr(:)
      real(8), allocatable :: gRRi(:)
      real(8), allocatable :: gLLr(:)
      real(8), allocatable :: gLLi(:)
C ... Local parameters
      integer i,j,ids,ltype,lmode,lgen,ltyp,lmod,lbulk
      character outs*80

C     call yprm('pgembs: gll',2,gll(1,1+ofL,1),ldl*ndgl,ldl,ldl,ldl)

      call awrit1('%x pgembs: mode=%i',outs,80,0,mode)
      call rxx(lgii == 0,'pgembs: no memory allocated for gii')
      ltype = 0
      lmode = 1
      lgen = 0
      if (lds /= ld0+ldl+ldr) call rx('pgemb: lds mismatch')
C     call yprm('Layer strux in pgembs',2,sii,ld0*lds,ld0,ld0,lds)

C --- g_ii <- -(S-P) ---
      do  j = 1, ld0
      do  i = 1, ld0
        gii(i,j+of0,1:2) = -sii(i,j,1:2)
      enddo
      enddo

C --- Add embedding term from left semi-infinite GF ---
      if (mod(mode/10,4) /= 2) then
        call lgupac(lgll,' pgembs using gll:',6,0,0,0,0,0)
        call awrit1('%a  lgll=%i',outs,80,0,lgll)
        ids = 1+ld0
        if (mod(lgLL/2,4) == 0) call rx('pgembs: missing gLL')
C   ... gll must be for surface or xtal (ltyp=1 or 3)
        call lgupac(lgll,'utmb',0,ltyp,lmod,lbulk,0,0)
        if (mod(lbulk,2) == 1) lgen = lgen + 1
C   ... Create and stored embedding term S_iL g_LL~ S_Li in sii
C       wk used as a work array
        if (ltyp == 1) then   ! g_LL contains left semi-infinite or related
          if (lmod == 0) then ! g_LL is GF
            call yygemm('N','N',ld0,ldl,ldl,1d0,sii(1,ids,1),sii(1,ids,2),ld0,
     .        gLL(1,1+ofL,1),gLL(1,1+ofL,2),ldl,0d0,wk,wk(1,1,2),ld0)
C     ... Subtract S_iL g_LL~ S_Li = S_iL (g_LL~^-1)^-1 S_iL+ from gii
          elseif (lmod == 1) then ! g_LL is GF^-1
C           Make local copy of g, since yyqnvb destroys it
            allocate(gLLr(ldl*ldl),gLLi(ldl*ldl))
            call dcopy(ldl*ldl,gLL(1,1+ofL,1),1,gLLr,1)
            call dcopy(ldl*ldl,gLL(1,1+ofL,2),1,gLLi,1)
            call dcopy(ldl*ld0,sii(1,ids,1),1,wk,1)
            call dcopy(ldl*ld0,sii(1,ids,2),1,wk(1,1,2),1)
            allocate(wkr((ldl+1)*max(ldl,ld0)))
            call yyqnvb('t',gLLr,gLLi,ldl,ldl,ld0,wkr,ldl,wkr,wk,wk(1,1,2),ld0,i)
            call rxx(i /= 0,'pgemb: GF singular')
            deallocate(wkr,gLLr,gLLi)
          elseif (lmod == 3) then
            print *, 'pgemb not ready for lmod=',lmod
            goto 999
          else
            goto 999
          endif
C     ... This could go straight into gii, but it's faster this way
          call yygemm('N','C',ld0,ld0,ldl,1d0,wk,wk(1,1,2),ld0,
     .      sii(1,ids,1),sii(1,ids,2),ld0,0d0,sii,sii(1,1,2),ld0)
          ltype = ltype + 1
C         call yprm('S_iL g_LL~ S_Li',2,sii,ld0*ld0,ld0,ld0,ld0)
        elseif (ltyp == 3) then
          print *, 'pgemb not ready for ltyp=3'
          goto 999
        else
          goto 999
        endif

C       call yprm('S_iL g~ S_Li',2,sii,ld0*lds,ld0,ld0,ld0)
C   ... Add embedding term S_iL g_LL~ S_Li to gii(of0)
        do  j = 1, ld0
        do  i = 1, ld0
          gii(i,j+of0,1:2) = gii(i,j+of0,1:2) - sii(i,j,1:2)
        enddo
        enddo
      endif

C --- Add embedding term from right semi-infinite GF ---
      if (mod(mode/10,4) /= 1) then
        call lgupac(lgrr,' pgembs using grr:',6,0,0,0,0,0)
        call awrit1('%a  lgrr=%i',outs,80,0,lgrr)
        ids = 1+ld0+ldl
        if (mod(lgRR/2,4) == 0) call rx('pgembs: missing gRR')
C   ... grr must be for surface or xtal (ltyp=2 or 3)
        call lgupac(lgrr,'utmb',0,ltyp,lmod,lbulk,0,0)
        if (lbulk/2 == 1) lgen = lgen + 4
C   ... Create and stored embedding term S_iR g_RR~ S_Ri in sii
C       wk used as a work array
        if (ltyp == 2) then
C     ... Subtract S_iR g_RR~ S_Ri = S_iR g_RR~ S_iR+ from gii
          if (lmod == 0) then
            call yygemm('N','N',ld0,ldr,ldr,1d0,sii(1,ids,1),sii(1,ids,2),ld0,
     .        gRR(1,1+ofr,1),gRR(1,1+ofr,2),ldr,0d0,wk,wk(1,1,2),ld0)
C     ... Subtract S_iR g_RR~ S_Ri = S_iR (g_RR~^-1)^-1 S_iR+ from gii
          elseif (lmod == 1) then
C           Make local copy of g, since yyqnvb destroys it
            allocate(gRRr(ldr*ldr))
            allocate(gRRi(ldr*ldr))
            call dcopy(ldr*ldr,gRR(1,1+ofr,1),1,gRRr,1)
            call dcopy(ldr*ldr,gRR(1,1+ofr,2),1,gRRi,1)
            call dcopy(ldr*ld0,sii(1,ids,1),1,wk,1)
            call dcopy(ldr*ld0,sii(1,ids,2),1,wk(1,1,2),1)
            allocate(wkr((ldr+1)*max(ldr,ld0)))
            call yyqnvb('t',gRRr,gRRi,ldr,ldr,ld0,wkr,ldr,wkr,wk,wk(1,1,2),ld0,i)
            call rxx(i /= 0,'pgemb: GF singular')
            deallocate(wkr,gRRr,gRRi)
          elseif (lmod == 2) then
            print *, 'pgemb not ready for lmod=',lmod
            goto 999
          else
            goto 999
          endif
C     ... This could go straight into gii, but it's faster this way
          call yygemm('N','C',ld0,ld0,ldr,1d0,wk,wk(1,1,2),ld0,
     .      sii(1,ids,1),sii(1,ids,2),ld0,0d0,sii,sii(1,1,2),ld0)
          ltype = ltype + 2
        elseif (ltyp == 3) then
          print *, 'pgemb not ready for ltyp=3'
          goto 999
        else
          goto 999
        endif
C       call yprm('S_iR g~ S_Ri',2,sii,ld0*lds,ld0,ld0,ld0)
        do  j = 1, ld0
        do  i = 1, ld0
          gii(i,j+of0,1:2) = gii(i,j+of0,1:2) - sii(i,j,1:2)
        enddo
        enddo
      endif

C --- gii = (P-S  -  S_iL g_LL~ S_iL+  -  S_iR g_RR~ S_iR+ )^-1 ---
      if (mod(mode,10) == 0) then
        call yyqinv('n',gii(1,1+of0,1),gii(1,1+of0,2),ld0,2,ld0,wk,ld0,i)
        call rxx(i /= 0,'pgemb: GF singular')
C       call yygefa(gii(1,1+of0,1),gii(1,1+of0,2),ld0,ld0,wk,i)
C       if (i /= 0) call fexit(-1,111,' Exit -1 PGEMB: '//
C     .    'Green''s function matrix singular',0)
C        call yygedi(gii(1,1+of0,1),gii(1,1+of0,2),ld0,ld0,wk,wk,
C     .    wk(1,2,1),wk(1,3,1),1)
        lmode = 0
      endif

      call lgupac(lgii,'p',1,ltype,lmode,0,0,lgen)
      call lgupac(lgii,' pgembs created',6,0,0,0,0,0)
C     call yprm('Layer strux in pgembs*',2,sii,ld0*lds,ld0,ld0,lds)

      return

C --- Error exit ---
  999 continue
      call rx(outs//'%a,  not implemented')

      end
      subroutine pvemb(s_ham,s_pot,s_str,s_site,sanity,plat,isp,nspc,mode,vpl,
     .  ip0,npl,pgplp,qp,ndg0,ndgl,ndgr,nrpt,ofL,of0,ofR,lgbii,lds,ldl,
     .  ldr,lgii,ipL,ipR,s00,gLL,gii,gRR,ld0,strRL,gij)
C- Kernel to make gii(ip0+nrpt-1)~
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham lncol lham neula offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  pgbevl plham plhamso
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfnc palp
Cio    Passed to:  pgbevl plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  pgbevl plham
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  pgbevl plham plhamso
Ci Inputs
Ci   sanity:
Ci   plat  :primitive lattice vectors, in units of alat
Ci   isp   :current spin channel (1 or 2)
Ci   mode  :passed to pginns
Ci         :1s digit does: (pginns)
Ci         : 0
Ci         : 1 put info into rvec needed for backward updates
Ci         : 2 Use info from rvec to make backward updates
Ci         :   (see mode, Outputs)
Ci         : 3 generate information needed for g0n,gn0
Ci         :   Also, if mode=3, call to pgofLR (pvemb)
Ci         :10s digit
Ci         :  1 adiabatic V
Ci   vpl   :Adiabatic potential shift (10s digit mode)
Ci   ip0   :make gii(ip0+nrpt-1) from gii(ipL)
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see subas.f
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   qp    :k-point
Ci   ld0,ndg0,of0  diagonal dim, tot. dim. and diag. offset to gii
Ci   ldl,ndgl,ofL  diagonal dim, tot. dim. and diag. offset to gLL
Ci   ldr,ndgr,ofR  diagonal dim, tot. dim. and diag. offset to gRR
Ci   nrpt  :Number of repeat layers
Ci   lds   :ldl+ld0+ldr
Ci   lgii  :attributes of the gii generated; see lgupac for conventions
Ci   ip0   :current PL to embed
Ci   ipL   :layer to left of current PL
Ci   ipR   :layer index to gRR
Ci         :if nrpt>0, should be ipl = ip0+nrpt
Ci   s00   :Structure contants for layer ip0
Ci   gLL   :gii(ipl)%gll
Co Outputs
Co   lgbii :contents of gii
Co   gii   :
Co   gRR   :Make gRR
Co         :Also, rvec-derived parameters are tacked to the end of gRR,
Co         :at gRR(1,ld0+1) when gRR is defined as a complex array (pginns)
Co         :(Note for now, gRR is defined locally as REAL)
Co  strrl  :(only if mode=3)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   08 Jul 13 replaced f77 pointers with f90 ones except passed pointers
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical sanity
      integer lds,ldl,ldr,ld0,ipL,ipR,ndg0,ndgl,ndgr,nrpt,ofL,of0,ofR,
     .  lgbii,lgii(-2:ipL),npl,ip0,mode,isp,nspc,pgplp(6,-1:10),strRL(4)
      double precision plat(3,3),qp(3),
     .  gLL(1),gii(ld0,ndg0,2),gRR(*),s00(*),vpl(-1:npl)
C     Pointer array to different layer GF's
      type (s_lgf) :: gij(4)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wkc(:),gbii(:)
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: wt(:)
      complex(8), allocatable :: rn(:)
      complex(8), allocatable :: sii(:)
      complex(8), allocatable :: r(:)
C ... Local parameters
      integer jpl,offRR,off21,off22,ip1,ip2
      integer ld0x,ldlx,ldrx,ldsx,ldm1x,ndg0x,ndgrx,ndglx
      double precision x

      if (nspc == 2) call rx('pvemb not ready for nspc=2')
      ld0x = ld0*nspc
      ldlx = ldl*nspc
      ldrx = ldr*nspc
      ldsx = lds*nspc
      ldm1x = pgplp(4,-1)*nspc
      ndg0x = ndg0*nspc
      ndgRx = ndgR*nspc
      ndgLx = ndgL*nspc
      allocate(r(2*ld0),sii(ld0*ld0*3))
      offRR = 1+2*ld0*ld0

C ... for now
      if (mode/10 == 1) then
        allocate(rn(2*ld0))
        ip2 = ip0+nrpt-2
        call pgbevl(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ip2,npl,pgplp,qp,100,-1,gRR(offRR),rn,sii,x,ld0,0,x,1)
      else
        allocate(rn(1))
      endif
      ip1 = ip0+1
C     gRR(offRR)<-rvec; sii<- P-S(ip1)
      call pgbevl(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ip1,npl,pgplp,qp,100,-1,gRR(offRR),r,sii,x,ld0,0,x,1)
      allocate(wt(nrpt))
      if (mode/10 == 1) then
        call pvemb2(ip0,ip0+nrpt-1,ip1,ip2,ld0,r,rn,vpl,wt)
      endif

C --- Make gii(ip0+nrpt-1)~ ---
C ... gbii<- Inverse right GF and sii(1..3) <- S0R, PrP^-1 and QrQ^-1
      allocate(wkc(lds*ld0),gbii(ld0**2))
      lgbii = 1
      call pgbulk(ld0,gRR(offRR),wkc,r,sii,9110,gRR,ld0,gbii,lgbii)
C ... gii(1,1+of0)<-g(ip0,ip0) for repeating pot. at ip0,ip0+1,ip0+2,...
      call pgembs(ldsx,ld0x,ldlx,ldrx,ndg0x,ndglx,ndgrx,of0,ofL,ofR,s00,
     .  30,lgii(ipL),gLL,lgbii,gbii,wkc,lgii(ip0),gii(1,1+of0,1))
      deallocate(gbii,wkc)
C     call yprm('gii',2,gii(1,1+of0,1),ld0*ndg0,ld0,ld0,ld0)
C ... gRR<-gii(ip0+nrpt-1)~ and gRR(offRR)<-rvec-related parameters for embedding (see pginns)
      allocate(wkc(ld0*ld0))
      call pginns(ld0,nrpt-1,r,gRR(offRR),rn,wt,lgii(ip0),
     .  gii(1,1+of0,1),gRR,sii,mode,wkc,lgii(ipR),gRR)
      deallocate(wkc,wt,rn,r,sii)

      if (mod(mode,10) /= 3) return

C --- Make  gj-1, g-1j ---
      allocate(wk(strRL(2)*ld0*2))
      off21 = 1 + 2*ld0*ld0 + ld0
      off22 = off21 + 2*ld0*ld0
      jpl = ip0+nrpt-1
      if (nspc == 2) call rx('pvemb not ready for nspc=2')
      call pgofLR(sanity,ld0x,ldlx,ldrx,ldm1x,s00,wk,lgii(jpl),ld0*2,
     .  gRR(off21),ld0*2,gRR(off22),strRL(1),strRL(2),gij(1)%gll,gij(2)%gll)
      deallocate(wk)

      end
      subroutine pvemb2(ipL,ipR,ip1,ip2,ld,r,rn,vpl,wt)
      implicit none
      integer ipl,ipr,ip1,ip2,ld
      double precision vpl(-1:ipR),wt(ipL:ipR)
      double precision r(ld,2,2),rn(ld,2,2)
      double precision dv
      double complex dlnr
      integer i,ip

C ... No effect if no delta v
      dv = vpl(ip2) - vpl(ip1)
      if (dv == 0) then
        call dpzero(rn,ld*4)
        call dpzero(wt,ipR-ipL+1)
        return
      endif

C      call yprm0('(1p,9e18.10)')
C      call yprm('r',2,r,ld*2,ld*2,ld*2,1)
C      call yprm('rn',2,rn,ld*2,ld*2,ld*2,1)

C ... Overwrite rn with delta ln r
      do  10  i = 1, ld*2
        dlnr = cdlog(dcmplx(rn(i,1,1),rn(i,1,2))) -
     .          cdlog(dcmplx(r(i,1,1),r(i,1,2)))
C       print *, i,dlnr
        rn(i,1,1) = dble(dlnr)
        rn(i,1,2) = dimag(dlnr)
   10 continue

C     call yprm('d lnr',2,rn,ld*2,ld*2,ld*2,1)

C ... Make weights
      do  20  ip = ipL, ipR
        wt(ip) = (vpl(ip) - vpl(ip1))/dv
   20 continue

      end
C      subroutine pgembi(mode,lds,lda,ldb,nda,ndb,ndab,offa,offb,offab,
C     .  offs,ldw,wk,sab,lgbb,gbb,gab,lgaa,gaa)
CC- Make interfacial g_ii from one semiinfinite and one proper GF
CC ----------------------------------------------------------------
CCi Inputs:
CCi   mode:    0 calculate gab and convert gaa~ to gaa
CCi            1 calculate gab only
CCi            2 gab is input ; convert gaa~ to gaa
CCi   sab(offs):strx connecting layers a and b (in sab(1,offs))
CCi   gaa:     semi-infinite GF gaa~
CCi   gbb:     proper diagonal GF, b layer
CCi   nda:     second dimension for gaa
CCi   ndab:    second dimension for gab
CCi   ndb:     second dimension for gbb
CCi   wk,ldw   workspace, size (lda,ldw,2); ldw must be >= max(lda,ldb)
CCo Outputs:
CCo   gab:     off-diagaonal GF gab
CCo   gaa:     (mode ne 1) proper crystal diagonal GF gaa
CCo   sab(1)   is OVERWRITTEN
CCr Remarks
CCr   Uses special case of general formulas from Dyson equation
CCr     g_AB = g_Aa~ Sab gbB = g_Aa Sab gbB~
CCr     g_BA = g_Bb~ Sba gaA = g_Bb Sba gaA~
CCr     g_AA = g_AA~  +  g_Ab  Sba g_aA~
CCr            g_AA~  +  g_Aa~ Sab g_bA
CCr   Whether gaa is right- or left- semi-infinite and gbb is the proper
CCr   GF depends only on which s is used: s0r when gaa~ is on the left;
CCr   use s0l when gaa~ is the on the right.
CC ----------------------------------------------------------------
C      implicit none
C      integer lds,lda,ldb,offs,ldw,lgaa,lgbb,nda,ndab,ndb,mode,
C     .  offa,offb,offab
C      double precision sab(lda,lds,2),wk(lda,ldw,2),
C     .  gaa(lda,nda,2),gbb(ldb,ndb,2),gab(lda,ndab,2)
C      logical sanity
C      parameter (sanity=.true.)
C
CC --- Make gab ---
C      if (mode < 2)
C     .  call pgdys2(sanity,2,'N',lds,lda,ldb,lda,nda,ndb,ndab,offa,
C     .  offb,offab,ldb,offs,wk,sab,lgaa,gaa,lgbb,gbb,gab)
C
C      if (mode == 1) return
C
CC --- S_ab(1) <- delta g_aa = g_ab S_ba g_aa~ = g_ab S_ab+ g_aa~ ---
C      call yygemm('N','C',lda,lda,ldb,1d0,gab(1,1+offab,1),
C     .  gab(1,1+offab,2),lda,sab(1,1+offs,1),sab(1,1+offs,2),lda,0d0,
C     .  wk,wk(1,1,2),lda)
C      call yygemm('N','N',lda,lda,lda,1d0,wk,wk(1,1,2),lda,gaa(1,1+offa,
C     .  1),gaa(1,1+offa,2),lda,0d0,sab(1,1,1),sab(1,1,2),lda)
CC     call yprm('g_ab S_ba gaa~',2,sab(1,1,1),lda*nda,lda,lda,lda)
C
CC --- Make g_aa = gaa~ + g_ab s_ba g_aa~ ---
C      call daxpy(lda*lda,1d0,sab(1,1,1),1,gaa(1,1+offa,1),1)
C      call daxpy(lda*lda,1d0,sab(1,1,2),1,gaa(1,1+offa,2),1)
C      call lgupac(lgaa,'pt',0,3,0,0,0,0)
C      call lgupac(lgaa,' pgembi created',6,0,0,0,0,0)
C
C      end
