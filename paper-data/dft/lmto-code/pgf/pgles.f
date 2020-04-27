      subroutine pgles(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ip0,npl,pgplp,ldmx,
     .  ld0n,ndg0,ld0,qp,bnn,gin,lgii,gii)
C- Make crystal diagonal g^<(ip0,ip0) from crystal non-diagonal retarded
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham lncol lham neula offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plham
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:palp pfnc
Cio    Passed to:  plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plham
Ci Inputs:
Ci   plat:  primitive lattice vectors, in units of alat
Ci   isp:   current spin index
Ci   ip0:   current PL
Ci   npl:   number of principal layers
Ci   pgplp: table of PL-dependent quantities (see pgfset)
Ci   ldmx:  maximum leading dimension of GF for -1,0,..,npl layers
Ci   ld0n:  leading dimension of pnl layer (the same as npl-1 layer)
Ci   ndg0:  column  dimension of ip0 layer; normally same as ld0
Ci   ld0:   dimension of lower block, layer ip0
Ci   qp:    k-point (NB: projection of qp onto PL direction not used)
Ci   bnn:   self energy for the right lead, bnn[npl-1,npl-1] =
Ci          S[npl-1,npl]*(g^r-g^a)[right semi-infin.]*S[npl,npl-1]
Ci   gii:   Diagonal retarded left semi-infinite GF for ip0 PL
Ci   lgii:  attributes of gii
Co Input/Output:
Co   gin:   non-diagonal retarded crystal g(ip0+1,npl-1) on input;
Co          g(ip0,npl-1) on output
Co Output:
Co   gii:   Diagonal crystal lesser (-1/2)*g^<(ip0,ip0) for new ip0 PL
Co   lgii:  attributes of gii generated; see lgupac for conventions
Cr Remarks
Cr  This routine is based on this paer:
Cr  S. Faleev et al Phys. Rev. B71}, 195422 (2005).
Cr
Cr  gin(ip0,npl-1) and right self-energy bnn, using formula
Cr  g^<(ip0,npl-1) = gin(ip0,npl-1)*bnn*hermit.conjug[gin(ip0,npl-1)].
Cr  gin(ip0,npl-1) is calculated from Dyson equation (for i < npl-1)
Cr  using gin(i,npl-1) = g~(i,i)*S(i,i+1)*gin(i+1,npl-1)       (1)
Cr  where g~(i,i) = g(left semi-inf.)(i,i)
Cr
Cr  Noncollinear version:
Cr  Suppressing spin indices, write Eq (1) as:
Cr    gin = g~L S g
Cr  Include spin indices and assume S is spin diagonal.
Cr    gin_1,2 = sum_3 g~L_1,3 S_3,3 g_3,2                      (2)
Cr            = sum_3 g~L_1,3 wk_3,2 where wk_3,2 = S_3,3 g_3,2
Cb Bugs
Cb   S(i,i+1) is assumed to be collinear and independent of spin
Cb   Thus spin quantization axis is assumed to be fixed along z.
Cb mch -f9f15.10 -qr out.copt $s22 s11 -sub 1,36,1,36 tmpx/g4nu -- -abs -max:g
Cb mch -f9f15.10 -qr out.copt $s22 s22 -sub 1,36,1,36 tmpx/g4nd -- -abs -max:g
Cu Updates
Cu   28 Aug 16 Rewritten to handle noncollinear case
Cu   18 May 16 First cut at extension to noncollinear case
Cu   10 Nov 11 Begin migration to f90 structures
Cu   16 Jan 04 (S.Faleev) First created
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nspc,npl,ip0,pgplp(6,-1:*),ldmx,ld0n,ld0,lgii,ndg0
      double precision plat(3,3),qp(3)
      double precision bnn(ld0n,nspc,ld0n,nspc,2),gii(ld0,nspc,ndg0,nspc,2)
      real(8), target :: gin(ldmx,nspc,ld0n,nspc,2)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      real(8), pointer ::    gwk(:,:,:,:,:)
      real(8),allocatable :: sij(:,:,:),wk(:,:,:,:,:)
C ... Local parameters
      integer off0,ldr,ldl,i,i1,ld0x,ldlx,ldrx,ld0nx,lds,ldgwk
      double precision sum

      call tcn('pgles')

C ... Dimensioning parameters
      off0 = pgplp(6,max(ip0,-1))
      ldr  = pgplp(4,min(ip0+1,npl))
      ldl  = pgplp(4,max(ip0-1,-1))
      lds  = ldl+ld0+ldr
      ld0x = ld0*nspc           ! dimension of noncollinear lower block, layer ip0
      ldlx = ldl*nspc           ! dimension of noncollinear lower block, layer ip0-1
      ldrx = ldr*nspc           ! dimension of noncollinear lower block, layer ip0+1
      ld0nx = ld0n*nspc         ! dimension of noncollinear lower block, layer npl

      if (ld0 /= ndg0 .or. off0 /= 0) call rx('pgles: wrong PL dimension')

C      if (isp == 1) then
C      print *, 'pgles : g connecting layers', min(ip0+1,npl-1), npl-1
C      call yprm('gin',2,gin,ldmx*nspc*ld0nx,ldmx*nspc,ldmx*nspc,ld0nx)
C      call yprm('bnn',2,bnn,(ld0nx)**2,ld0nx,ld0nx,ld0nx)
C      endif

C --- Make gin(i,npl-1) = g(left semi-inf.)(i,i)*S(i,i+1)*gin(i+1,npl-1) ---
      if (ip0 <= npl-2 .and. ip0 >= -1) then
        if (mod(lgii/2,4) /= 1) call rx('pgles: gii is not left semi-infinite')

C   ... Make S(i,i+1) and pack gin(i+1,npl) into contiguous gwk
        allocate(sij(ld0x,lds*nspc,2),wk(ld0,nspc,lds,nspc,2))

        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ip0,npl,ld0,lds,pgplp,qp,wk,sij)
        if (nspc == 2) then
          allocate(gwk(ldr,nspc,ld0n,nspc,2))
          call pgsetg(1,ldmx,ld0n,ldr,ld0n,0,nspc,gwk,gin)
          ldgwk = ldrx
        else
          gwk => gin
          ldgwk = ldmx
        endif
C       call yprm('Sij',2,sij,ld0x*lds*nspc,ld0x,ld0x,lds*nspc)
        deallocate(wk)

C   ... wk <-- S(i,i+1) * gin(i+1,npl-1)
C       Multiply using dimensions wk(ld0x,ld0nx) <- s(ld0x,ldrx) * g(ldrx,ld0nx)
C set s22 = '-split s 1,36+1,72+1 1,36+1,72+1'
C set s22x = '-split s 1,54+1,108+1 1,36+1,72+1'
C set s22y = '-split s 1,54+1,108+1 1,54+1,108+1'
C set s34 = '-split s 1,108+1 1,108+1,216+1,288+1'
C mch out.copt $s22x s11 tmpx/g4nu -- -px
C mch out.copt $s22 s11 tmpx/bnnu -- -px:6
C matches coll because 15 doesn't couple to any of sites in ipl=4
        allocate(wk(ld0,nspc,ld0n,nspc,2))
C       call yprm('Sii+1',2,sij(1,ldlx+ld0x+1,1),ld0x*lds*nspc,ld0x,ld0x,ldrx)
        call yygemm('N','N',ld0x,ld0nx,ldrx,1d0,
     .    sij(1,ldlx+ld0x+1,1),sij(1,ldlx+ld0x+1,2),ld0x,
     .    gwk,gwk(1,1,1,1,2),ldgwk,0d0,wk,wk(1,1,1,1,2),ld0x)
C       call yprm('S*gin',2,wk,ld0x*ld0nx,ld0x,ld0x,ld0nx)
        deallocate(sij)

C   ... gin(i,npl-1) <-- g~(i,i) * S(i,i+1) * gin(i+1,npl-1) ... dimensioned (ld0x,ld0nx)
        if (nspc == 2) then
          deallocate(gwk)
          allocate(gwk(ld0,nspc,ld0n,nspc,2))
          ldgwk = ld0x
        endif

C       Multiply using dimensions gin(ld0x,ld0nx) <- gii~(ld0x,ld0x) * wk(ld0x,ld0nx)
C       call yprm('gii~',2,gii,ld0x*ld0x,ld0x,ld0x,ld0x)
        call yygemm('N','N',ld0x,ld0nx,ld0x,1d0,
     .    gii(1,1,1+off0,1,1),gii(1,1,1+off0,1,2),ld0x,
     .    wk,wk(1,1,1,1,2),ld0x,0d0,gwk,gwk(1,1,1,1,2),ldgwk)
C       call yprm('gin(i,npl-1)',2,gwk,ld0x*ld0nx,ld0x,ld0x,ld0nx)
        deallocate(wk)

C   ... Pack contiguous gwk into gin(i+1,npl)
        if (nspc == 2) then
          call pgsetg(0,ldmx,ld0n,ld0,ld0n,0,nspc,gwk,gin)
        endif

      elseif (ip0 == npl-1) then  ! End layer already has gin
C       For ip0=npl-1 gin(ip0,npl-1) is already given
C       Noncoll case: copy gin into contiguous gwk
        if (nspc == 2) then
          allocate(gwk(ld0,nspc,ld0n,nspc,2))
          call pgsetg(1,ldmx,ld0n,ld0,ld0n,0,nspc,gwk,gin)
          ldgwk = ld0x
        else
          gwk => gin
          ldgwk = ldmx
        endif

C       call yprm('gwk',2,gwk,ld0x*ld0nx,ld0x,ld0x,ld0nx)

      else
        call rx('pgles: wrong ipl')
      endif

C --- Make (-1/2)*g^<(ip0,ip0) ---
C     wk <-- gin(i,npl-1) * bnn (npl-1,npl-1)
C     Multiply using dimensions wk(ld0x,ld0nx) <- gin(ld0x,ld0nx) * bnn(ld0nx,ld0nx)
      allocate(wk(ld0,nspc,ld0n,nspc,2))
      call yygemm('N','N',ld0x,ld0nx,ld0nx,1d0,
     .  gwk,gwk(1,1,1,1,2),ldgwk,
     .  bnn,bnn(1,1,1,1,2),ld0nx,0d0,wk,wk(1,1,1,1,2),ld0x)
C     call yprm('gin*bnn',2,wk,ld0x*ld0nx,ld0x,ld0x,ld0nx)

C     gii <-- -1/2 gin(i,npl-1) * bnn (npl-1,npl-1) * [gin(i,npl-1)]+
C ... gii contains (-1/2)*g^<(ip0,ip0), diagonal over orbit and PL index
c     There is no "-" because two i in definition of self-energy and
c     lesser GF g^<, see paper by S.Faleev cited in Remarks
C     Multiply using dimensions gii(ld0x,ld0x) <- wk(ld0x,ld0nx) * [gin(ld0x,ld0nx)]+
      call yygemm('N','C',ld0x,ld0x,ld0nx,1d0/2,
     .  wk,wk(1,1,1,1,2),ld0x,
     .  gwk,gwk(1,1,1,1,2),ldgwk,
     .  0d0,
     .  gii,gii(1,1,1,1,2),ld0x)
C     call yprm('1/2 gin*bnn*gin+',2,wk,ld0x*ld0x,ld0x,ld0x,ld0x)

      if (nspc == 2) deallocate(gwk)

      sum = 0d0
      do  i1 = 1, nspc
        do  i = 1, ld0
          sum = sum + abs(gii(i,i1,i,i1,1))
          gii(i,i1,i,i1,1) = 0
        enddo
      enddo
      if (sum > 1d-8) call rx('pgles: Re[g^<] not zero')

C     Label crystal (-1/2)*g^<  as 'interface GF' (analogous to crystal g^r)
      call lgupac(lgii,'pt',0,3,0,0,0,0)

C     print *, 'pgles : g connecting layers', ip0, npl-1
C     call yprm('gin',2,gin,ldmx*nspc*ld0nx,ldmx*nspc,ldmx*nspc,ld0nx)

      deallocate(wk)
      call tcx('pgles')

      end
