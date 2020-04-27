      subroutine gfg2g(s_ham,s_pot,mode,lrel,isp,ib1,ib2,jb1,jb2,
     .  ldg1,ldg2,gf,ghh,npos)
C- Convert g_ij = to G_ij by energy scaling
C-----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 lncol ldham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:papg gmar palp pf gma dpfr ddpfr dpf ddpf
Cio    Passed to:  *
Ci Inputs
Ci   mode :1s digit pertains to storage (not used now)
Ci        10s digit distinguishes how complex arithmetic is handled
Ci          0: gf has real, imaginary separated
Ci             gf = gf(ldg1,ldg2,2), with gf(*,*,1..2) = real..imag
Ci          1: gf is returned complex*16 format:
Ci             gf = gf(2,ldg1,ldg2), with gf(1..2,*,*) = real..imag
Ci          2: gf has real, imaginary separated by columns
Ci             gf = gf(ldg1,2,ldg2), with gf(*,1..2,*) = real..imag
Ci       100s digit
Ci         0  no scaling
Ci         1  scale g -> G
Ci         2  alp-to-gam transformation g^alp -> g^gam
Ci         3  both 2 + 1
Ci         4  non-equilibrium mode, scale g^< to G^<, alpha repns only
Ci       1000s digit
Ci         1  reverse sense of scaling: G -> g
Ci         2  suppress diagonal contribution P-dot-dot
Ci         3  1 + 2
Ci   lrel  :0 for non-relativistic
Ci          1 for scalar relativistic
Ci          2 for Dirac equation
Ci ib1,ib2:scale subblock in range ib1..ib2, row dimension
Ci jb1,jb2:scale subblock in range jb1..jb2, column dimension
Ci   ldg1  :leading dimension of gf.  In the noncollinear case
Ci         :the leading dimension of gf for one spin channel.
Ci   ldg2  :second dimension of gf.  In the noncollinear case
Ci         :the second dimension of gf for one spin channel.
Cio Inputs/Outputs
Cio  gf    :On input gf is unscaled g, eg (P-S)^-1 in 2nd gen LMTO
Cio        :On output gf is scaled g
Cio  ghh   :On input gf is unscaled hh block of g
Cio        :On output gf is scaled hh block of g
Co Outputs
Co   npos  :number of elements for which g_ij > 0
Cl Local variables
Cl   ldpf  :leading dimension of pfun (suhamz.f)
Cl   ofhi  :offset to imaginary part of ghh
Cl   offgi :offset to imaginary part of gf
Cr Remarks
Cr   This hamiltonian-dependent routine scales the dimensionless GF
Cr   the proper GF.
Cr   For 2nd-generation LMTO:
Cr       g_in  = (P-S)^-1
Cr       g_out  scaled g; see gfg2g2
Cu Updates
Cu  18 Jun 18 Synchronize with updated spherical harmonics
Cu  05 Jun 16 no special treatment for nonequilibrium branch
Cu     May 16 redesign of scaling for SO and FR cases
Cu  10 Nov 11 Begin migration to f90 structures
Cu   1 Jun 05 New options 1000s digit mode
Cu  16 Jan 04 (S.Faleev) Changes to handle non-equilibrium mode
Cu  18 Mar 03 (A Chantis) first cut at relativistic case
Cu             Altered argument list.
Cu  26 Feb 02 Adapted to handle case hamiltonian permuted by site
Cu  28 Apr 00 changes to adapt to pgf
Cu            Argument list changed
Cu  15 Mar 00 redesigned the iwaves downfolding, added h-waves
Cu  21 Dec 99  added option to transform GF to gamma-representation
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isp,ib1,ib2,jb1,jb2,ldg1,ldg2,npos,lrel
      double precision gf(ldg1,ldg2),ghh(*)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer PRTG,hdim,hdimx,i,ib,is1,is2,ixx,kcplx,ld11,ld21,ldgx,ldham(16),
     .  ldi,ldj,ldpf,ldr,ldtm,lgen3,lhdim,lidim,lidimx,lncol,mode2,mode3,
     .  mxorb,nspc,offgi,offi,offpi,offpj,ofgh,ofgl,ofhi,pfdim,nl,morder
      logical lso
      equivalence (lidim,ldham(2)),(lhdim,ldham(3))
      equivalence (nspc,ldham(4)),(lidimx,ldham(6))
      equivalence (pfdim,ldham(9))
      parameter (PRTG=100)
      procedure(logical) bittst
      procedure(integer) :: iprint,nglob
C ... External calls
      external cplxdm,gf2gam,gfg2g2,gfr2gam,info0,isanrg,offsHp,rxi,rxx,soscale,yprm,ztoyy

      call sanrg(.true.,mod(mode,100),0,0,'gfg2g','mode')
      kcplx = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      if (mode2 == 0) return
      lgen3 = s_ham%lgen3
      lncol = s_ham%lncol
      ldham = s_ham%ldham
C     Never an independent spin channel in noncollinear case
      if (isp == 2) call sanrg(.true.,nspc,1,1,'gfg2g','nspc')
      ldgx  = ldg1*nspc
      nl = nglob('nl')
      mxorb = nglob('mxorb')
      lso = bittst(lncol,4)
      call lmorder(0,morder,[0],[0])

C      if (iprint() >= PRTG) then
C        print *, 'ib1..ib2',ib1,ib2
C        ldtm = ldg2*nspc
C        call yprmi('gf before scaling (ib1,ib2)=%s,(%2i)',[ib1,ib2],0,2,gf,ldgx*ldtm,ldgx,ldgx,ldtm)
C        if (lhdim-lidim /= 0)
C     .    call yprm('gh before scaling',2,ghh,ofhi,hdimx,hdimx,mxorb)
C      endif


C --- 3nd generation NMTO --
      if (lgen3 /= 0) then
        call rx('gfg2g not ready for 3rd generation')

C --- 2nd generation ASA : Scale (P-S)^-1 in gf ---
      else
C       ldpf  = lhdim
        ldpf  = pfdim
        hdim = lhdim - lidim
        hdimx = hdim*nspc

C       call yprm('offH',0,s_ham%offH,0,n0H*nkap0,n0H,32)
C       call snot(s_pot%dpf,ldpf*nspc)
C       call zprm('dotP',2,s_pot%dpf,ldpf*nspc,ldpf*nspc,1)
C       call zprm('dotdotP',2,s_pot%ddpf,ldpf*nspc,ldpf*nspc,1)

C       Loop over four spin channels in noncollinear case
C       In collinear case, is1=is2=isp
C       Both pf and gf are in downfolding order; thus no permutations
        do  is1 = isp, max(isp,nspc)
        do  is2 = isp, max(isp,nspc)

C         Fast version, but no need to bother
C         call offsHp(s_ham%iprmb,ib1,ib2,1,ldham(1),lidim,offpi,ixx)
C         call offsHp(s_ham%iprmb,ib1,ib2,1,0,ldham(1),offpi,ldi)
C         ldi = ldi+ixx
          call offsHp(s_ham%iprmb,ib1,ib2,0,0,lidim,offpi,ldi)
          offpi = offpi + (is1-1)*ldpf
          call offsHp(s_ham%iprmb,jb1,jb2,0,0,lidim,offpj,ldj)
          offpj = offpj + (is2-1)*ldpf

C         print 346, 'ib1,ib2',ib1,ib2,offpi,ldi
C 346     format(a,4i5)

          call cplxdm(kcplx,ldg1*nspc,ldg1*nspc,ld11,ld21,ldr,offgi)
C         Case g contains only spin channel isp
          if (nspc == 1) then
            ofgl = 0
            ofgh = 0
          else
            ofgl = ldg1*(is1-1) + ldg1*nspc*ldg2*(is2-1)
            ofgh = hdim*(is1-1)
          endif

C     ... Transform to gamma representation
C         call yprm('gf',kcplx+2,gf,(ldg1*nspc)**2,ldg1*nspc,ldg1*nspc,ldg1*nspc)

          if (mode2 == 2 .or. mode2 == 3) then
            call info0(PRTG,0,0,' gfg2g: g(alp)->gamma repsn')
            if (mode3 /= 0) call rxi('gfg2g not ready mode=',mode)
            if (lrel == 2 .or. lso) then
              if (ldg1/=ldi .or. ldg1/=ldj)
     .          call rx('gfr2gam not ready for subblock dimensioning')
              i = 100*morder + kcplx
              call gfr2gam(i,ib1,ib2,s_ham%iprmb,ldg1,ldi,ldj,offpi,offpj,pfdim,
     .          s_pot%papg,s_pot%gmar,gf,npos)
            else
              call gf2gam(0,ldgx,offgi,ldi,ldj,offpi,offpj,s_pot%palp,
     .          s_pot%pf,s_pot%gma,gf(1+ofgl,1))
            endif
C           Transform the higher block
            if (lhdim > lidim .and. offpi == offpj) then
              if (lso) call rx('gfg2g not ready for higher block')
              mxorb = nglob('mxorb')
              call cplxdm(kcplx,hdim*nspc,mxorb,ixx,ixx,hdimx,ofhi)
              do  ib = ib1, ib2
                call offsHp(s_ham%iprmb,ib,ib,0,lidim,lhdim,offi,ldi)
                offi  = offi - lidim
                offpi = offi + lidim + (is1-1)*ldpf
                call gf2gam(0,hdimx,ofhi,ldi,ldi,offpi,offpi,
     .            s_pot%palp,s_pot%pf,s_pot%gma,ghh(1+offi+ofgh))
              enddo
            endif

C            if (iprint() >= PRTG) then
C              ldtm = ldg2*nspc
C              call yprm('gf after a2g',2,gf,ldgx*ldtm,ldgx,ldgx,ldtm)
C              if (hdim /= 0)
C     .        call yprm('gh after a2g',2,ghh,ofhi,hdimx,hdimx,mxorb)
C            endif

          endif

C     ... Transform g to G by energy scaling (mode2=1,3,4,5,6)
          if (mod(mode2,2) /= 0 .or. mode2 == 4) then
            call info0(PRTG,0,0,' gfg2g: g->G')
            i = 0
            if (mode3 == 1) i = 1
            if (mode3 == 2) i = 10
            if (mode3 == 3) i = 11

            if (lrel == 2 .or. lso) then
              i = 100*morder + kcplx
              if (mode3 /= 0) call rxi('gfg2g not ready mode=',mode)
              call gfr2gam(i,ib1,ib2,s_ham%iprmb,ldg1,ldi,ldj,offpi,offpj,pfdim,
     .          s_pot%dpfr,s_pot%ddpfr,gf,npos)
            else
              call gfg2g2(i,ldgx,offgi,ldi,ldj,offpi,offpj,s_pot%dpf,
     .          s_pot%ddpf,gf(1+ofgl,1),npos)
            endif
C           Transform the higher block
            if (lhdim > lidim .and. offpi == offpj) then
              if (lso) call rx('gfg2g not ready for higher block')
              mxorb = nglob('mxorb')
              call cplxdm(kcplx,hdim*nspc,mxorb,ixx,ixx,hdimx,ofhi)
              do  ib = ib1, ib2
              offi = s_ham%offH(1+2+nkap0*n0H*(ib-1),1)
              ldi  = s_ham%offH(1+2+nkap0*n0H*ib,1) - offi
C             offi  = w(ooffH+2+nkap0*n0H*(ib-1))
C             ldi   = w(ooffH+2+nkap0*n0H*ib) - offi
              call offsHp(s_ham%iprmb,ib,ib,0,lidim,lhdim,offi,ldi)
              offi  = offi - lidim
              offpi = offi + (is1-1)*ldpf + lidim
              call gfg2g2(i,hdimx,ofhi,ldi,ldi,offpi,offpi,s_pot%dpf,
     .          s_pot%ddpf,ghh(1+offi+ofgh),npos)
              enddo
            endif

C         Non-equilibrium mode
C          elseif (mode2 == 4)  then
C            if (mode3 /= 0) call rxi('gfg2g not ready mode=',mode)
C            if (lrel >= 2)
C     .      call rx('gfg2g not ready for lrel>1 for non-equilibr. mode')
C            if (lhdim > lidim  .or. ldgx /= ldi)
C     .      call rx('gfg2g not ready for high block in non-equil. mode')
C            call gfg2gn(ldi,offgi,offpi,s_pot%dpf,gf(1+ofgl,1),npos)
          endif

          if (lrel == 2 .or. lso) goto 17
        enddo  ! loop over is2
        enddo  ! loop over is1
   17   continue
      endif

      if (iprint() >= PRTG/1) then
        ldtm = ldg2*nspc
        call yprmi('gf after scaling (ib1,ib2)=%s,(%2i)',[ib1,ib2],0,2,gf,ldgx*ldtm,ldgx,ldgx,ldtm)
        if (hdim /= 0)
     .    call yprm('gh after scaling',2,ghh,ofhi,hdimx,hdimx,mxorb)
      endif
      end

      subroutine gfg2g2(mode,ldg1,offi,ni,nj,offpi,offpj,dpf,ddpf,gf,npos)
C- Convert 2nd generation g_ij = (P-S)^-1 to G_ij by energy scaling
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :  0 convert g_ij to G_ij
Ci         :  1 convert G_ij to g_ij
Ci         :10s digit
Ci         :  0 include diagonal contribution ddpf
Ci         :  1 do not include diagonal contribution ddpf
Ci   ldg1  :leading dimension of g_ij
Ci   offi  :separation between real and imaginary parts
Ci   ni    :row dimension of gf to scale
Ci   nj    :column dimension of gf to scale
Ci   offpi :offset to first entry in row index of dpf,ddpf
Ci   offpj :offset to first entry in column index of dpf,ddpf
Ci   dpf   :1st energy derivative of potential function P-dot (mkptfp.f)
Ci          dpf is in complex*16 format
Ci   ddpf  :2nd energy derivative  -1/2 P-dotdot/P-dot (mkptfp.f)
Ci          ddpf is in complex*16 format
Cio Inputs/Outputs
Cio  gf    :On input gf is unscaled g eg (P-S)^-1 in 2nd gen LMTO
Cio        :On output gf(i,j) is scaled by
Cio        :    sqrt(dpf(i+offpi)) sqrt(dpf(j+offpj)) +
Cio        :          ddpf(i+offpi) delta_ij
Co Outputs
Co   npos  :(mode=0) number of diagonal elements for which G_ij > 0
Cr Remarks
Cr   In the two-center, gamma repsn, sqrt(P-dot) = 1/sqrt(delta)
Cu Updates
Cu   1 Jun 05 Re-definition of mode
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ldg1,offi,ni,nj,offpi,offpj,npos
      double precision gf(ldg1,*),ddpf(2,nj*2)
      double complex dpf(nj*2)
C ... Local parameters
      double complex xxc
      integer i,j

      npos = 0

C --- G -> g: subtract diagonal part -1/2 P^dotdot_i / P^dot_i
      if (mode == 1 .and. offpi == offpj) then
        do  i = 1, ni
          gf(i,i)      = gf(i,i)      - ddpf(1,i+offpi)
          gf(i+offi,i) = gf(i+offi,i) - ddpf(2,i+offpi)
        enddo
      endif

C --- g_ij <- sqrt(P^dot_i) g_ij sqrt(P^dot_j) ---
      if (mod(mode,10) == 0) then
        do  j = 1, nj
        do  i = 1, ni
          xxc = sqrt(dpf(i+offpi)*dpf(j+offpj))
     .         *dcmplx(gf(i,j),gf(i+offi,j))
          gf(i,j)      = dble(xxc)
          gf(i+offi,j) = dimag(xxc)
        enddo
        enddo
      else
        do  j = 1, nj
        do  i = 1, ni
          xxc = dcmplx(gf(i,j),gf(i+offi,j))
     .         /sqrt(dpf(i+offpi)*dpf(j+offpj))
          gf(i,j)      = dble(xxc)
          gf(i+offi,j) = dimag(xxc)
        enddo
        enddo
      endif

C --- g_ij <- g_ij - 1/2 P^dotdot_i / P^dot_i ---
      if (mode == 0 .and. offpi == offpj) then
        do  i = 1, ni
          gf(i,i)      = gf(i,i)      + ddpf(1,i+offpi)
          gf(i+offi,i) = gf(i+offi,i) + ddpf(2,i+offpi)
          if (gf(i+offi,i) > 0) npos = npos+1
        enddo
      endif

      end

      subroutine gfr2gam(mode,ib1,ib2,iprmb,ldg,ldi,ldj,offpi,offpj,pfdim,paopg,gma,gf,npos)
C- 2nd generation scaling transformation alpha -> gamma, fully relativistic version
C- Or alternatively g = (P-S)^-1 -> G
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :kcplx distinguishes how gf stores complex matrices
Ci          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
Ci          1: complex*16: gf = gf(ldg,ldg)
Ci          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
Ci         :10s digit
Ci          0: SR convention for m ordering
Ci          1: FR convention for m ordering
Ci         :100s digit
Cu         :1: internally permute m ordering before transformation
Ci   ldg   :leading dimension of g_ij
Ci   ldi   :row dimension of gf to scale.  Not used now
Ci   ldj   :column dimension of gf to scale.  Not used now
Ci   offpi :offset to first entry in row index of dpf,ddpf
Ci   offpj :offset to first entry in column index of dpf,ddpf
Ci   pfdim :leading dimension of paopg and gma
Ci   paopg :g^alp->g^gam : palp pgam^-1 (ratios of P in 2 representations)
Ci         :g -> G       : A^T  where A^T A = energy derivative Pdot of pgam
Ci   gma   :g^alp->g^gam : (gamma-alpha) palp/pgam
Ci         :g -> G       : -1/2 P-dotdot/P-dot in gamma repsn
Cio Inputs/Outputs
Cio  gf    :On output gf(i,j) is scaled by
Cio        :  gma(i)  +  paopg^T (i+offpi)  gf(i,j)  paopg (j+offpj)
Co   npos  :number of diagonal elements for which G_ij > 0
Cr Remarks
Cu Updates
Cu   11 May 18 Added 100s digit option for m ordering
Cu   19 Mar 16 (MvS) redesigned to work with SO coupling.  gfg2g3 subsumed into this routine
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ib1,ib2,ldg,ldi,ldj,offpi,offpj,pfdim,npos,iprmb(pfdim)
      double complex gf(ldg,2,ldg,2),paopg(pfdim,2,2),gma(pfdim,2,2)
C ... Local parameters
      integer i,is,ldgx,kcplx,mod1,mod2,nl
      procedure(integer) :: nglob
C     double complex temp(ldg,2,ldg,2)
C     logical :: debug=.false.

      ldgx = ldg*2
      kcplx = mod(mode,10)
      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)
      call rxx(ldg/=ldi .or. ldg/=ldj .or. offpi/=offpj,
     .  'gfr2gam not ready for off-diagonal subblocks')

      call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)

C     if (debug) call yprmi('(P-S)^-1-rel, alpha ib1,2 = %2:1i',[ib1,ib2],0,3,gf,0,ldgx,ldgx,ldgx)

C     For now, just permute m order for socscl.
C     Better to modify soscale to handle reversed m order.
      if (mod2 > 0) then
        nl = nglob('nl')
C       call zprm('dpfr',2,paopg,pfdim,pfdim,4)
        call pmorderx(1,ib1,ib2,1,nl,iprmb,0,pfdim,pfdim,4,pfdim,4,0,paopg)
C       call zprm('dpfr permuted',2,paopg,pfdim,pfdim,4)
C       call zprm('ddpfr',2,gma,pfdim,pfdim,4)
        call pmorderx(1,ib1,ib2,1,nl,iprmb,0,pfdim,pfdim,4,pfdim,4,0,gma)
C       call zprm('ddpfr permuted',2,gma,pfdim,pfdim,4)
C       call zprm('(P-S)^-1-rel, alpha',2,gf,ldgx,ldgx,ldgx)
        call pmorderx(11,ib1,ib2,2,nl,iprmb,0,ldg,ldg,ldg,ldi,ldj,offpi,gf)
C       call zprm('(P-S)^-1-rel, permuted',2,gf,ldgx,ldgx,ldgx)
C       if (debug) call yprmi('(P-S)^-1-rel, permuted ib1,2 = %2:1i',[ib1,ib2],0,3,gf,0,ldgx,ldgx,ldgx)
      endif

C ... gf <- (paopg)^T gf (paopg)
C     call zprm('paopg calling soscale',2,paopg,pfdim,pfdim,4)
      call soscale(mode,ldi,ldj,pfdim,offpi,offpj,paopg,gf)
C     call frscale(ldg,pfdim,paopg,gf)  ! superseded
C     if (debug) call yprmi('after soscale Pa/Pg (P-S)^-1 Pa/Pg ib1,2 = %2:1i',[ib1,ib2],0,3,gf,0,ldgx,ldgx,ldgx)
C     stop

C     debugging
C     print *, '!! gfr2gam return early'; return
C     print *, '!! gfr2gam'; gf=0

C ... gf += gma  (diagonal part of g)
      npos = 0
      if (offpi == offpj) then
      do  i = 1, ldg  ! spin diagonal parts for the lower block
        do  is = 1, 2
          gf(i,is,i,is) = gf(i,is,i,is) + gma(i+offpi,is,is)
          if (dimag(gf(i,is,i,is)) > 0d0) npos = npos + 1
        enddo
      enddo
      do  i = 2, ldg  ! spin off-diagonal parts for the lower block
        if (mod1 == 1) then
          gf(i-1,1,i,2) = gf(i-1,1,i,2) + gma(i+offpi-1,1,2)
          gf(i,2,i-1,1) = gf(i,2,i-1,1) + gma(i+offpi,2,1)
        else
          gf(i,1,i-1,2) = gf(i,1,i-1,2) + gma(i+offpi-1,1,2)
          gf(i-1,2,i,1) = gf(i-1,2,i,1) + gma(i+offpi,2,1)
        endif
      enddo
      endif

      if (mod2 > 0) then
        call pmorderx(1,ib1,ib2,1,nl,iprmb,0,pfdim,pfdim,4,pfdim,4,0,paopg)
        call pmorderx(1,ib1,ib2,1,nl,iprmb,0,pfdim,pfdim,4,pfdim,4,0,gma)
        call pmorderx(11,ib1,ib2,2,nl,iprmb,0,ldg,ldg,ldg,ldi,ldj,offpi,gf)
C       call zprm('dpfr orig',2,paopg,pfdim,pfdim,4)
C       call zprm('ddpfr orig',2,gma,pfdim,pfdim,4)
      endif

C     if (debug) call yprmi('(P-S)^-1-rel, gamma ib1,2 = %2:1i',[ib1,ib2],0,3,gf,0,ldgx,ldgx,ldgx)
      call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)

      end

C      subroutine gfg2g3(mode,ldg,ni,nj,offpi,offpj,pfdim,dpffr,ddpffr,gf,npos)
CC- 2nd generation g_ij = (P-S)^-1 to G_ij by energy scaling, fully relativistic version
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1s digit
CCi         :kcplx distinguishes how gf stores complex matrices
CCi          0: real, imaginary separated: gf = gf(ldg,ldg,1..2)
CCi          1: complex*16: gf = gf(ldg,ldg)
CCi          2: real, imaginary in columns : gf = gf(ldg,1..2,ldg)
CCi         :10s digit
CCi          0: SR convention for m ordering
CCi          1: FR convention for m ordering
CCi   ldg   :leading dimension of g_ij
CCi   ni    :row dimension of gf to scale ... not used now
CCi   nj    :column dimension of gf to scale ... not used now
CCi   offpi :offset to first entry in row index of gf
CCi   offpj :offset to first entry in column index of gf
CCi   pfdim :leading dimension of the potential parameter arrays
CCi   dpffr :1st energy derivative of potential function P-dot (mkfrpf.f)
CCi          dpf is in complex*16 format
CCi   ddpffr:2nd energy derivative  -1/2 P-dotdot/P-dot (mkfrpf.f)
CCi          ddpf is in complex*16 format
CCio Inputs/Outputs
CCio  gf    :On input gf is unscaled g eg (P-S)^-1 in 2nd gen LMTO
CCio        :On output gf(i,j) is scaled by
CCio        :    sqrt(dpf(i+offpi)) sqrt(dpf(j+offpj)) +
CCio        :          ddpf(i+offpi) delta_ij
CCo Outputs
CCo   npos  :number of diagonal elements for which g_ij > 0
CCr Remarks
CCr   dpffr and ddpffr are stored in kappa rep
CCr   (for each orbital: {{{mu+1/2,mu+1/2},{mu+1/2,mu-1/2}},
CCr                       {{mu-1/2,mu+1/2},{mu-1/2,mu+1/2}},
CCr   spins numbered in order -1/2,1/2; thus m+s=mu for all elements)
CCr   offpi,offpj seem to be unnecessary
CCr   ni,nj seem to be unnecessary
CCr   mode is not used
CCr   Further cleanup required, perhaps when downfolding is enabled
CCb Bugs
CCb   Scaling assumes lambda and kappa indices are equivalent, which is wrong
CCb   Complete revision is necesssary
CCu Updates
CCu  19 Mar 16 merged into gfrgam
CCu  13 Aug 13 (K Belashchenko) Simplified
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,ldg,ni,nj,offpi,offpj,npos,pfdim
C      double complex dpffr(pfdim,2,2),ddpffr(pfdim,2,2),gf(ldg,2,ldg,2)
CC ... Local parameters
C      integer i,is,ldgx,kcplx,mod1
C
CC     call sanrg(.true.,mode,0,0,'gfg2g3','mode')
C
C      ldgx = ldg*2
C      kcplx = mod(mode,10)
C      mod1 = mod(mode/10,10)
C      call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,kcplx,1)
C
CC     call zprm('(P-S)^-1',2,gf,ldgx,ldgx,ldgx)
C
C      if (offpi /= 0 .or. offpj /= 0) call rx('Unexpected offsets')
C
CC --- g_ij <- sqrt(P^dot_i) g_ij sqrt(P^dot_j) , 2x2 matrices ---
C      call frscale(ldg,pfdim,dpffr,gf)
C
CC     call zprm('sqrt(pd)(P-S)sqrt(pd)',2,gf,ldgx,ldgx,ldgx)
C
CC ... Add (-1/2)(P_dot_dot)/P_dot
C      do  i = 1, ldg  ! spin diagonal parts for the lower block
C        do  is = 1, 2
C          gf(i,is,i,is) = gf(i,is,i,is) + ddpffr(i+offpi,is,is)
C        enddo
C        if (i == 1) cycle
C        gf(i-1,1,i,2) = gf(i-1,1,i,2) + ddpffr(i+offpi-1,1,2)
C        gf(i,2,i-1,1) = gf(i,2,i-1,1) + ddpffr(i+offpi,2,1)
C      enddo
C
CC     call zprm('Green-function',2,gf,ldgx,ldgx,ldgx)
C
C      call ztoyy(gf,ldgx,ldgx,ldgx,ldgx,1,kcplx)
C
C      end

C      subroutine gfg2gn(ldi,offi,offpi,dpf,gf,npos)
CC- Convert 2nd generation diagonal g^< to G^< by energy scaling
CC ----------------------------------------------------------------------
CCi   ldi   :leading dimension of g_ii
CCi   offi  :separation between real and imaginary parts
CCi   offpi :offset to first entry in row index of gf
CCi   dpf   :1st energy derivative of potential function P-dot (mkptfp.f)
CCi          dpf is in complex*16 format; alpha repsn.
CCio Inputs/Outputs
CCio  gf    :On input diagonal gf is unscaled crystal g^< in 2nd gen LMTO
CCio        :On output diagonal gf is scaled by
CCio        :    sqrt( dpf(i+offpi)) * conjug[sqrt(dpf(i+offpi)] )
CCo Outputs
CCo   npos  :number of diagonal elements for which g_ii > 0
CCr Remarks
CCr   No longer needed ... dpf is now generated with Im z = 0
CCu Updates
CCu   16 Jan 04 (S.Faleev) First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer ldi,offi,offpi,npos,i
C      double precision gf(ldi,*)
C      double complex   dpf(ldi)
CC ... Local parameters
C      double complex xxc
C
C      npos = 0
C      do i = 1, ldi
C         xxc = sqrt( dpf(i+offpi)* dconjg(dpf(i+offpi)) )
C     .         *dcmplx(0d0,gf(i+offi,i))
CC        print *, gf(i,i), dpf(i+offpi)
C         gf(i,i) = 0d0
C         gf(i+offi,i) = dimag(xxc)
CC         if (abs(sqrt(dpf(i+offpi)* dconjg(dpf(i+offpi)))-  !TEMP!!!!!!!
CC     .        abs(dpf(i+offpi))) > 1d-10 .or. offi /= ldi**2)
CC     .         call rx('gfg2gn: wrong something')         !TEMP!!!!!!!
C         if (gf(i+offi,i) > 0d0) npos = npos + 1
C      enddo
C      end

C      subroutine frscale(ldg,pfdim,del,g)
CC- Scale 2x2 matrix g <- transpose(del) * g * del in fully relativistic case
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ldg   :leading dimension of g
CCi   pfdim :leading dimension of del
CCi   del   :1st energy derivative of potential function P-dot (mkptfp.g)
CCio Inputs/Outputs
CCio  g     :function to scale (overwritten)
CCr Remarks
CCu Updates
CCu  19 Mar 16 (MvS) superseded by soscale
CCu  13 Aug 13 (K Belashchenko) First created
CC ----------------------------------------------------------------------
C      implicit none
C      integer ldg,pfdim
C      double complex del(pfdim,2,2),g(ldg,2,ldg,2)
C      double complex wk(ldg,2,ldg,2),delT(ldg,2,2)
C      integer i,j,k,ip,jp,kp,s,sp
C
C      if (ldg > pfdim) call rx('frscale: wrong dimensions')
C
C      wk = dcmplx(0d0,0d0)
C
C      delT = dcmplx(0d0,0d0)
C      do j = 1, 2
C        delT(1:ldg,j,j) = del(1:ldg,j,j)
C      enddo
C
C      delT(3:ldg,2,1)   = del(2:ldg-1,1,2)
C      delT(2:ldg-1,1,2) = del(3:ldg,2,1)
C
C      do s = 1, 2
C      do sp = 1, 2
C        do k = 1, 2
C        do kp = 1, 2
C          do i = 1, ldg
CC      ...  Sum rule (conserved mu: i-s = j1-k)
CC      ...  (spins are numbered -1/2, 1/2; m's are l,l-1,...,-l)
C            ip = k - s + i
CC      ...  check we are within bounds
C            if (ip < 1 .or. ip > ldg) cycle
C            do j = 1, ldg
C              jp = kp - sp + j
C              if (jp < 1 .or. jp > ldg) cycle
CC        ...  accumulate in wk array
C              wk(i,s,j,sp) = wk(i,s,j,sp) + delT(i,s,k)
C     .          * g(ip,k,jp,kp) * del(jp,kp,sp)
C
C              if (i == 6 .and. s == 1 .and. j == 17 .and. sp == 2) then
C                print 321, i,s,j,sp,delT(i,s,k),ip,k,jp,kp,g(ip,k,jp,kp),del(jp,kp,sp),g(ip,k,jp,kp)*del(jp,kp,sp)
C  321           format(4i3,2f15.10,2x,4i3,6f15.10)
C              endif
C
C
C            enddo !j
C          enddo !i
C        enddo !kp
C        enddo !k
C      enddo !sp
C      enddo !s
C
C      g = wk
C
C      end

C      subroutine frscale(ldg,pfdim,del,f)
CC- Scale 2x2 matrix f <- transpose(del) * f * del in fully relativistic case
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ldg   :leading dimension of f
CCi   pfdim :leading dimension of del
CCi   del   :1st energy derivative of potential function P-dot (mkptfp.f)
CCio Inputs/Outputs
CCio  f     :function to scale (overwritten)
CCr Remarks
CCu Updates
CCu  13 Aug 13 (K Belashchenko) First created
CC ----------------------------------------------------------------------
C      integer ldg,pfdim
C      double complex del(pfdim,2,2),f(ldg,2,ldg,2)
C
C      double complex wk(ldg,2,ldg,2),delT(ldg,2,2)
C      integer s1,s2,sp1,sp2,i,j,i1,i2
C
C      if (ldg > pfdim) call rx('frscale: wrong dimensions')
C
C      wk = dcmplx(0d0,0d0)
C
C      delT = dcmplx(0d0,0d0)
C      do j = 1, 2
C        delT(1:ldg,j,j) = del(1:ldg,j,j)
C      enddo
C
C      delT(3:ldg,2,1)   = del(2:ldg-1,1,2)
C      delT(2:ldg-1,1,2) = del(3:ldg,2,1)
C
C      do s1 = 1, 2
C      do s2 = 1, 2
C        do sp1 = 1, 2
C        do sp2 = 1, 2
C          do i = 1, ldg
CC      ...  Sum rule (conserved mu: i-s1 = j1-sp1)
CC      ...  (spins are numbered -1/2, 1/2; m's are l,l-1,...,-l)
C            i1 = sp1 - s1 + i
CC      ...  check we are within bounds
C            if (i1 < 1 .or. i1 > ldg) cycle
C            do j = 1, ldg
C              j1 = sp2 - s2 + j
C              if (j1 < 1 .or. j1 > ldg) cycle
CC        ...  accumulate in wk array
C              wk(i,s1,j,s2) = wk(i,s1,j,s2) + delT(i,s1,sp1)
C     .          * f(i1,sp1,j1,sp2) * del(j1,sp2,s2)
C            enddo !j
C          enddo !i
C        enddo !sp2
C        enddo !sp1
C      enddo !s2
C      enddo !s1
C
C      f = wk
C
C      end

C      subroutine frscale(ldg,pfdim,del,f)
CC- Scale 2x2 matrix f <- transpose(del) * f * del in fully relativistic case
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ldg   :leading dimension of f
CCi   pfdim :leading dimension of del
CCi   del   :1st energy derivative of potential function P-dot (mkptfp.f)
CCio Inputs/Outputs
CCio  f     :function to scale (overwritten)
CCr Remarks
CCu Updates
CCu  13 Aug 13 (K Belashchenko) First created
CC ----------------------------------------------------------------------
C      integer ldg,pfdim
C      double complex del(pfdim,2,2),f(ldg,2,ldg,2)
C
C      double complex wk(ldg,2,ldg,2),delR(ldg,2,2),delT(ldg,2,2)
C      integer s1,s2,sp1,sp2,i,j,i1
C      if (ldg > pfdim) call rx('frscale: wrong dimensions')
C
C      mod1 = 0
C
CC     call dpzero(delR,2*size(delR))
C
CC      do  i = 2, ldg
CC        if (mod1 == 1) then
CC          gf(i-1,1,i,2) = gf(i-1,1,i,2) + del(i-1,1,2)
CC          gf(i,2,i-1,1) = gf(i,2,i-1,1) + del(i,2,1)
CC        else
CC          gf(i,1,i-1,2) = gf(i,1,i-1,2) + del(i-1,1,2)
CC          gf(i-1,2,i,1) = gf(i-1,2,i,1) + del(i,2,1)
CC        endif
CC      enddo
CC
C
C      wk = dcmplx(0d0,0d0)
C      delT = dcmplx(0d0,0d0)
C      do  j = 1, 2
C        delR(1:ldg,j,j) = del(1:ldg,j,j)
C        delT(1:ldg,j,j) = del(1:ldg,j,j)
C      enddo
C
C      if (mod1 == 0) then
C        delR(1:ldg-1,2,1) = del(2:ldg,2,1); delR(ldg,2,1) = 0
C      else
C        delR(1:ldg,2,1) = del(1:ldg,2,1)
C
C        delT(3:ldg,2,1)   = delR(2:ldg-1,1,2)
C        delT(2:ldg-1,1,2) = delR(3:ldg,2,1)
C
C      endif
C
C      delT(3:ldg,2,1)   = del(2:ldg-1,1,2)
C      delT(2:ldg-1,1,2) = del(3:ldg,2,1)
C
CC      print *, '!!'; f = 0
CC      do  i1 = 1, 2*ldg
CC        f(i1,1,i1,1) = 1
CC      enddo
C
C      do  s = 1, 2
C      do  sp = 1, 2
C        do  k = 1, 2
C        do  kp = 1, 2
C          do  i = 1, ldg
CC      ...  Sum rule (conserved mu: i-s = jp-k)
CC      ...  (spins are numbered -1/2, 1/2; m's are l,l-1,...,-l)
C            ip = k - s + i
CC      ...  check we are within bounds
C            if (ip < 1 .or. ip > ldg) cycle
C            do  j = 1, ldg
C              jp = kp - sp + j
C              if (jp < 1 .or. jp > ldg) cycle
CC        ...  accumulate in wk array
CC              wk(i,s,j,sp) = wk(i,s,j,sp) +
CC     .          delT(i,s,k) * f(ip,k,jp,kp) * delR(jp,kp,sp)
C
C              if (ip == 6 .and. jp == 5) then
CC                print *, 'hi',k,kp,s,sp,ip,jp,delR(jp,kp,sp)
C                print *, 'hi',k,kp,kp,s,ip,jp,delR(jp,kp,sp)
C              endif
C              wk(i,s,j,sp) = wk(i,s,j,sp) +
C     .          f(ip,k,jp,kp) * delR(jp,kp,sp)
C
C
C
C            enddo !j
C          enddo !i
C        enddo !kp
C        enddo !k
C      enddo !sp
C      enddo !s
C
C      f = wk
C
C      end


      subroutine gf2gam(mode,ldg1,offi,ni,nj,offpi,offpj,palp,pgam,gma,gf)
C- 2nd generation scaling transformation alpha -> gamma
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit = 1 scale the diagonal elements of g only
Ci         :10s digit = 1 gf is a one-dimensional matrix
Ci         :Must be used in conjuction with 1s digit = 1
Ci   ldg1  :leading dimension of g_ij
Ci   offi  :separation between real and imaginary parts
Ci   ni    :row dimension of gf to scale
Ci   nj    :column dimension of gf to scale
Ci   offpi :offset to first entry in row index of gf
Ci   offpj :offset to first entry in column index of gf
Ci   palp  :potential function P in alpha rep'sn
Ci   pgam  :potential function P in gamma rep'sn
Ci   gma   :gamma - alpha
Cio Inputs/Outputs
Cio  gf    :On input gf is unscaled g eg (P-S)^-1 in 2nd gen LMTO
Cio        :On output gf(i,j) is scaled by
Cio        :  gma(i) + palp/pgam (i+offpi) gf(i,j)  palp/pgam (j+offpj)
Cr Remarks
Cr   Testing: need files palp-nc pgam-nc gma-nc galp-nc ggam-nc
Cr   mch palp-nc pgam-nc -de -a:nr=64 ratio ratio -v2dia -a pabyg gma-nc -a:nr=64 gma
Cr   pabyg galp-nc pabyg -x -x gma -v2dia pabyg -x -+ ggam-nc  -- -px
Cr
Cr   Testing salp->sgam: gma-nc need files gma-nc msalp-nc msgam-unrot
Cr   mch gma-nc -a:nr=64 gma  msalp-nc -i gma -v2dia -+ -i msgam-unrot -- -px
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer mode,ldg1,offi,ni,nj,offpi,offpj
      double precision gf(ldg1,*),gma(nj)
      double complex palp(nj),pgam(nj)
C ... Passed parameters
      double complex xx
      integer i,j

C      call zprm('(P-S)^(-1)scal',1,gf,ldg1,ldg1,ldg1)

C --- g_ij <- P^alp_i/P^gam_i  g_ij   P^alp_j/P^gam_j
      if (mode == 0) then
      do  10  j = 1, nj
      do  10  i = 1, ni
        xx = palp(i+offpi)/pgam(i+offpi)
     .      *dcmplx(gf(i,j),gf(i+offi,j))
     .      *palp(j+offpj)/pgam(j+offpj)
        gf(i,j) = dble(xx)
        gf(i+offi,j) = dimag(xx)
   10 continue
      else
        do  11  i = 1, ni
          if (mode < 10) then
            xx = palp(i+offpi)/pgam(i+offpi)
     .        *dcmplx(gf(i,i),gf(i+offi,i))
     .        *palp(i+offpi)/pgam(i+offpi)
            gf(i,i) = dble(xx)
            gf(i+offi,i) = dimag(xx)
          else
            xx = palp(i+offpi)/pgam(i+offpi)
     .        *dcmplx(gf(i,1),gf(i+offi,1))
     .        *palp(i+offpi)/pgam(i+offpi)
            gf(i,1) = dble(xx)
            gf(i+offi,1) = dimag(xx)
          endif
   11   continue
      endif

C      call zprm('Pa/Pg(P-S)Pa/Pg-scal',1,gf,ldg1,ldg1,ldg1)

C --- g_ij += (gamma-alpha) Palp/Pgam ---
      if (offpi == offpj) then
        do  20  i = 1, ni
          xx = gma(i+offpi) * palp(i+offpi)/pgam(i+offpi)
          if (mode < 10) then
            gf(i,i) = gf(i,i) + dble(xx)
            gf(i+offi,i) = gf(i+offi,i) + dimag(xx)
          else
            gf(i,1) = gf(i,1) + dble(xx)
            gf(i+offi,1) = gf(i+offi,1) + dimag(xx)
          endif
   20   continue
      endif

C      call zprm('gf-scalar',1,gf,ldg1,ldg1,ldg1)

      end

C     subroutine offsH(ipa,irow,ikap,ib,offH,offset)
C      integer nkap0,n0H,offset
C      parameter (nkap0=4,n0H=5)
C      integer irow,ikap,ib,ipa(*)
C      integer offH(n0H,nkap0,ib)
C      integer kb
C
C      kb = ib
C      if (ipa(1) /= 0) then
C        kb = ipa(ib)
C      endif
C      offset = offH(irow,ikap,kb)
