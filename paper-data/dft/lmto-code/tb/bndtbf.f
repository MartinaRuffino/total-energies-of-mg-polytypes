      subroutine bndtbf(s_ctrl,s_lat,nterm,nlmesp,nset,nclass,
     .  nl,nsp,nspc,nkp,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,cut,cutov,
     .  iam,npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .  dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,
     .  eband,epar,dpar)
C- Get eigenvalues and their parameter derivatives for each k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nterm,nlmesp,nset,nclass,nl,nsp,nkp,ldim,mxnbr,memode
Ci   nelts: first dimension of delta
Ci   nfit: number of eigenvalues to fit out of ldim*nkp total bands
Ci   nvar: number of parameters to vary out of all parameters
Ci   iam,npm,lmx,indxsh
Ci   ip1(nterm+1,nlmesp,nset,4): pointer to locations in full list of
Ci     variables for parameter types 1-4 (Hamiltonian and overlap
Ci     parameters).  A first index of nterm+1 corresponds to the decay
Ci     parameter (memode=2 or 3 only).
Ci   ip2(nl,nsp,nclass,5:6): pointer to locations in full list of
Ci     variables for parameter types 5 and 6 (diagonal Hamiltonian and
Ci     spin-orbit parameters).  For the spin-orbit parameters a first
Ci     index of 1 corresponds to p, 2 to d, etc.  A first index of nl
Ci     is not used (there is no spin-orbit coupling for s-states).
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   ivar(1,i): points to the position in the full list of ith variable
Ci   ivar(2,i): paramter type of ith variable, with types 1 to 6:
Ci     1: Hamiltonian parameter
Ci     2: Hamiltonian crystal field parameter
Ci     3: overlap parameter
Ci     4: overlap crystal field parameter
Ci     5: diagonal Hamiltonian parameter
Ci     6: spin-orbit parameter
Ci   rmaxh,efmax,qp,hso
Ci   decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,vso: parameters
Ci   delta: holds electrostatic increments to the diagonal matrix
Co Outputs
Co   eband: energy bands for all k-points
Co   epar: subset of eband, evals to be included in fit (all k-points)
Co   dpar: derivatives of epar wrt each of 1..nvar TB parameters
Cr Remarks
Cr   The fitting works with ME type independent on a particular pair of
Cr   species, whereas tbham can cope with different ME for different
Cr   species pairs. Consequently, tbham expects memode as an array, not
Cr   as a scalar. A quick fix for this mismatch employed here is to make
Cr   an array memodk with all elements equal and pass it into tbham.
Cr   Cutoff-related arrays cutmd and poly are dealt with similarly
Cr   (see a few lines before the tbham call).
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nterm,nlmesp,nset,nclass,nl,nsp,nkp,nelts,ldim,nfit,mxnbr,
     .  memode,nvar
      integer ip1(nterm+1,nlmesp,nset,4),iam(3,1),npm(2,nclass),
     .  ip2(nl,nsp,nclass,5:6),lmx(nclass),ifit(2,nkp),indxsh(1),
     .  ivar(2,nvar)
      logical rl
      double precision rmaxh,efmax
      double precision decay(nlmesp,nset),deccf(nlmesp,nset),
     .  decov(nlmesp,nset),dcocf(nlmesp,nset),tabme(nterm,nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),tabov(nterm,nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),qnu(3,nl,nsp,nclass),delta(nelts,1),
     .  qp(3,nkp),vso(nl,nclass),hso(nl**2,nl**2,4,2),
     .  eband(ldim*nsp,nkp),epar(nfit),dpar(nvar,nfit),pot0(*),
     .  cut(2,nlmesp,1),cutov(2,nlmesp,1)
C ... Dynamically allocated arrays
      integer, allocatable :: ntab(:),iax(:),npr(:),pgfsl(:)
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: zll(:),zwk(:)
      real(8), allocatable, target :: h(:),dh(:),h0(:),ov(:)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
C ... Local parameters
      integer i,j,iprint,mxcsiz,nttab,nsite,i1mach,ipr,nspc,ldimx,iev,
     .  ikp,iprx,nev,ltb,lncol,lpgf(2),nbas,npadl,npadr,nclasp,nbasp,
     .  fopn,isp,niax,os
      parameter (niax=10)
C     integer ontab,oiax,onpr,oh,oh0,odh,oov,owk,os,ozll
      double precision ckbas,cksumf,plat(3,3),kap2(20)
      logical lov,cryf,ocryf,addsll,lso,fitpar,bittst,ldum1,ldum2,
     .  iostrx,cmdopt
      character*120 outs
C ... need only for the tbham call
      integer poly(nset),cutmd(nset),memodk(nset)

      type(str_str0) :: strh
      integer :: inull
      real(8) :: rnull

C --- Setup ---
      ltb = s_ctrl%ltb
      lncol = s_ctrl%lncol
      lpgf = s_ctrl%lpgf
      nbas = s_ctrl%nbas
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = s_ctrl%nclasp
      nbasp = nbas + npadl + npadr
      plat = s_lat%plat
      ckbas = cksumf(s_lat%pos,3*nbas)

C --- Initialization ---
      lov    = bittst(ltb,1)
      cryf   = bittst(ltb,2)
      ocryf  = bittst(ltb,4)
      addsll = bittst(ltb,8) .and. lov
      lso    = bittst(lncol,4)

      call rxx(lso .and. nsp /= 2,'BNDTBF: must set NSPIN=2 for SO=T')

C --- Get neighbor table iax for each atom in the cluster ---
      if (lpgf(1) /= 0) then
        i = 2
        j = 1
      else
        i = 3
        j = -1
      endif
C ... Make nttab,ontab,oiax
      mxcsiz = mxnbr
C     call defi(ontab, nbasp+1)
C     call defi(oiax, nbasp+1)
      allocate(ntab(nbasp+1))
      allocate(pgfsl(nbasp+1))
      call pairs(nbas,nbasp,1d0,plat,[rmaxh/2],s_lat%pos,
     .  [-1],i,j,pgfsl,nttab,ntab,iax,mxcsiz)
      deallocate(pgfsl)

C ... Patch iax(6) for the padded basis
      if (nbasp > nbas) call pairp6(nbas,npadl,npadr,iax,ntab)
      nsite = nttab
      call mkiaxd(nsite,lmx,s_ctrl%ips,iax)

C ... tb programs not converted to ntab from npr yet ...
      allocate(npr(2*nbas))
      call npr2tb(1,nbasp,npr,ntab)

C --- Slater-Koster Hamiltonian ---
      fitpar = .true.
c ... no cutoff for fitting
      poly(1:nset) = 5
      cutmd(1:nset) = 0
      memodk(1:nset) = memode
      allocate(h(nl**4*nsite*min(nspc**2*nsp,4))); h = 0d0
      allocate(h0(nl**4*nsite*min(nspc**2*nsp,4))); h0 = 0d0
      allocate(ov(nl**4*nsite*nspc**2)); ov = 0d0
      allocate(dh(nvar*nl**4*nsite*nspc**2)); dh = 0d0

      call tbham(nsp,nspc,s_ctrl,s_lat,.true.,fitpar,nvar,
     .  ip1,ivar,nlmesp,memodk,decay,deccf,decov,dcocf,
     .  poly,cutmd,cut,cutov,iam,npm,nterm,nset,
     .  tabme,tabcf,tabov,tbocf,nsite,npr,iax,h,h0,dh,ov,1,1,1)

C ... diagonal elements
      ldum1 = bittst(ltb,2**15)
      if (ldum1) ltb = ltb - 2**15
      if (nspc == 2) call rx('update tbdiag; eliminate nsp')
      call tbdiag(ltb,nbasp,nl,nsp,nspc,nsp,s_ctrl%ipc,nsite,npr,
     .  s_ctrl%initc,qnu,nelts,delta,rnull,h,h0,ov)
      if (ldum1) ltb = ltb + 2**15

C --- Accumulate parameter derivatives for diagonal hamiltonian MEs ---
      call tbfdia(nbas,nl,nsp,nclass,nsite,npr,ip2,s_ctrl%ipc,dh)

C ... What about fit for spin-orbit parameters?

C --- Add to Hamiltonian MEs: hLL' -> hLL' + ebarLL' * sLL' ---
      if (addsll) then
C ... Fit does not work here yet, should it?
        ldum1 = bittst(ltb,16)
        ldum2 = bittst(ltb,128)
        if (ldum1) ltb = ltb - 16
        if (ldum2) ltb = ltb - 128
        allocate(wk(nl*nbas))
        call addes(nsite,nl,nsp,nbasp,ltb,npr,iax,ov,0,wk,h,0)
        deallocate(wk)
        if (ldum1) ltb = ltb + 16
        if (ldum2) ltb = ltb + 128
      endif

C --- Write Hamiltonian to STR file ---
      if (lpgf(1) /= 0 .or. cmdopt('-dumph',6,0,outs)) then
        call npr2tb(0,nbasp,npr,npr)
        allocate(strh % a(nl*nl*nbasp)); strh % a = 0.0_8
        allocate(wk(nbasp*nl**2)); call dpzero(wk,nbasp*nl**2)
        strh % s => h

        call iostrt(1,'STR',nl,nbasp,1,kap2,0,ckbas,-1,nsite,strh)
        nsite = strh % n(1+nbasp)
        call fclose(fopn('STR'))
        call npr2tb(1,nbasp,npr,npr)
        call free_str0(strh)
      endif

C --- Printout ---
      if (iprint() >= 50) then
        allocate(wk(nbasp*nl**2)); call dpzero(wk,nbasp*nl**2)
C   ... Printout Hamiltonian
        print *
        print *,'Real-space TB hamiltonian :'
        do  10  j = 1, nsp**2
          if (nsp == 2) print *
          if (j == 1 .and. nsp == 2)
     .      print *,'Spin: Up-Up, w/o spin-orbit'
          if (j == 2) print *,'Spin: Up-Down, w/o spin-orbit'
          if (j == 3) print *,'Spin: Down-Up, w/o spin-orbit'
          if (j == 4) print *,'Spin: Down-Down, w/o spin-orbit'
          os = 1 + (j-1)*nl**4*nsite
          call shostr(nl**2,nsite,nbasp,plat,s_lat%pos,0,wk,iax,
     .      ntab,h(os),1,1,0d0,1,1d0)
C     ... Printout Overlap
          if (lov) then
            print *,'overlap matrix'
            os = 1 + (j-1)*nl**4*nsite
            call shostr(nl**2,nsite,nbasp,plat,s_lat%pos,0,wk,
     .        iax,ntab,ov(os),1,1,0d0,1,1d0)
          endif
   10   continue
C   ... Printout parameter derivatives
        do  30  i = 1, nvar
          print *
          if (ivar(2,i) == 1 .or. ivar(2,i) == 5
     .                        .or. (ivar(2,i) == 6 .and. lso)) then
            write(*,500) ivar(1,i)
  500       format(' Derivatives for Hamiltonian parameter: ',i4)
          elseif (ivar(2,i) == 3 .and. lov) then
            write(*,510) ivar(1,i)
  510       format(' Derivatives for overlap parameter: ',i4)
          elseif (ivar(2,i) == 2 .and. cryf) then
            write(*,520) ivar(1,i)
  520       format(' Derivatives for crystal field parameter: ',i4)
          elseif (ivar(2,i) == 4 .and. ocryf) then
            write(*,530) ivar(1,i)
  530       format(' Derivatives for overlap crystal field parameter: ',
     .        i4)
          endif
          do  20  j = 1, nsp**2
            if (nsp == 2) print *
            if (j == 1 .and. nsp == 2) print *,'Spin: Up-Up'
            if (j == 2) print *,'Spin: Up-Down'
            if (j == 3) print *,'Spin: Down-Up'
            if (j == 4) print *,'Spin: Down-Down'
            os = 1 + (nsp**2*(i-1)+(j-1))*nl**4*nsite
            call shostr(nl**2,nsite,nbasp,plat,s_lat%pos,0,wk,
     .        iax,ntab,dh(os),1,1,0d0,1,1d0)
   20     continue
   30   continue
        deallocate(wk)
      endif

C --- Determine verbosity for secmtb ---
      ipr = 0
      if (iprint() >= 40) ipr = 1
      if (iprint() >= 45) ipr = 2
      if (iprint() >= 50) ipr = 3

C --- Loop over k-points ---
      ldimx = ldim*nspc
      allocate(zll(ldimx*ldimx))
      call dpzero(dpar,nvar*nfit)
      iev = 0
      do  40  ikp = 1, nkp
        if (ifit(1,ikp) > ifit(2,ikp)) goto 40
        iprx = ipr
        if (ikp >= 1) iprx = max(ipr-1,0)
        if (iprint() >= 40 .and. mod(ikp,100) == 0) iprx=max(ipr,1)

C --- Make bands and eigenvectors ---
        if (nspc == 2 .or. nsp == 2) call rx('need isp here')
        isp = 1
        call secmtbf(s_ctrl,
     .    plat,nbas,nl,nspc,nsp,isp,lmx,s_ctrl%ipc,indxsh,
     .    ldimx,ifit(2,ikp),efmax,ikp,nkp,qp(1,ikp),nsite,iax,
     .    npr,h,vso,hso,ov,pot0,rl,iprx,.true.,nev,
     .    zll,eband(1,ikp),rnull,rnull)

C --- Make derivatives of eigenvalues wrt TB parameters ---
        allocate(zwk(ldimx*ldimx))
        call derivp(nkp,nbas,nclass,nsite,nl,nspc,ldimx,nvar,nfit,iprx,
     .    ikp,nev,ifit,s_ctrl%ipc,ivar,lmx,indxsh,iax,npr,
     .    qp(1,ikp),plat,dh,zll,eband(1,ikp),lso,lov,cryf,ocryf,
     .    zwk,iev,dpar)
        deallocate(zwk)
   40 continue
      deallocate(npr,ntab,iax,zll)

C --- Copy bands to be fit into a single array ---
      call xcopb(nkp,ldimx,nfit,ifit,eband,epar)

      deallocate(h,h0,dh,ov)

      end
