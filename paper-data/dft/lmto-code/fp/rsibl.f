      subroutine rsibl(s_site,s_spec,s_lat,s_optic,lfrce,nbas,isp,q,iq,
     .  ndimh,nspc,napw,igapw,iprmb,numq,nevec,evec,ewgt,k1,k2,k3,smpot,
     .  smrho,f)
C- Add smooth part of output density into smrho and forces.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rsibl1
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  tbhsi uspecb rsibl1
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat gmax nabc ng vol kv kv2 igv igv2
Co     Stored:     igv igv2
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:igv gv kv qlat igv2
Cio    Passed to:  sugvec
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic ocrng unrng optme
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:optme ocrng unrng
Cio    Passed to:  opt_nstate
Ci Inputs
Ci   lfrce :1s digit
Ci         :if nonzero, accumulate contribution to force
Ci         :10s digit [NOT IMPLEMENTED]
Ci         :if nonzero, do not make psi, or accumulate smrho or force contribution
Ci         :   (useful if, e.g. only optical matrix elements are sought)
Ci   nbas  :size of basis
Ci   q     :Bloch vector
Ci   iq    :index to current k-point (used by MPI)
Ci   ndimh :dimension of hamiltonian
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial Fermi levels
Ci   nevec :number of eigenvectors with nonzero weights
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   k1..3 :dimensions smpot,smrho
Ci   smpot :smooth potential on uniform mesh, needed for forces
Co Outputs
Co   smrho :smooth density accumulated for this qp
Co         :Noncollinear case:
Co         :smrho(:,:,:,:,1:2,1) are the spin-diagonal parts
Co         :smrho(:,:,:,:,1:2,2) are the spin (1,2) and spin (2,1) blocks
Co   f     :force contribution accumulated for this qp
Co   optme(m,i,j) (an element in s_optic). If loptic>0, adds the interstitial part of
Co         :<i|grad|j> connecting occ i with unocc j, polarization m
Co         :The total <i|grad|j> has three terms:
Co         : + matrix element of the true w.f. (calculated here)
Co         : - matrix element of the local rep of the envelope w.f. (calculated here)
Co         : + matrix element of the envelope w.f. (calculated in rsibl)
Cl Local variables
Cl  loptic :if nonzero, calculate optical matrix elements
Cl  optl   :Interstitial part of matrix element <i | grad | j> connecting
Cl         :all occ states i, ocrng(1)<=i<=ocrng(2) to
Cl         :to unocc states j, unrng(1)<=j<=unrng(2)
Cr Remarks
Cm MPI
Cm   Parallelise over the eigenvector loop. The vector block size is
Cm   chosen (in the range 6-16, by dstrbp.f) so as to distribute the
Cm   work optimally across processes. Two work arrays of the size of
Cm   smrho are allocated from the heap as buffers. Only one will be
Cm   needed under MPI-2. See comments in hsibl.
Cb Bugs
Cb    openmp doesn't work ... crashes when calling fftz3 in rsibl2.
Cb    Possibly a problem with calling fftw within a rsibl2 which
Cb    is itself a thread?
Cb    Replace call to gvgvcomp and pass ivp as input
Cu Updates
Cu   21 Aug 14 Optics reworked to reflect the proper three-fold representation
Cu             of optical matrix elements
Cu   10 Feb 14 Optical matrix elements for envelopes generated if loptic>0
Cu   14 May 13 Complete migration to f90 structures
cu             Unsuccessful attempt (again) to get openmp to work
Cu   05 May 12 Repackaged G vector generation into sugvec
Cu   10 Nov 11 Begin migration to f90 structures
Cu   09 Apr 10 First attempt at noncollinear smrho.
Cu             Not checked; force terms not implemented
Cu   29 Dec 08 Unsuccessful attempt to make work with openmp
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu   10 Sep 06 Added MPI parallelization in the spin-coupled case
Cu   23 Dec 04 Extended to spin-coupled case
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   12 Oct 00 Use q-dependent list of G vectors
Cu    6 Jul 00 attempt to vectorize by grouping eigenvectors in blocks
Cu   17 Jun 00 Spin polarized
Cu   23 May 00 Adapted from nfp rsif_q.f
C ----------------------------------------------------------------------
      use structures
      implicit none
      integer lfrce,nbas,isp,k1,k2,k3,ndimh,nevec,numq,iprmb(*),iq,nspc
      integer napw,igapw(3,napw)
      double precision q(3),ewgt(numq,nevec),f(3,nbas,numq)
      double complex evec(ndimh,nspc,nevec),smrho(k1,k2,k3,numq,2,nspc),
     .  smpot(k1,k2,k3,isp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_optic):: s_optic
C ... Dynamically allocated local arrays
      real(8), allocatable :: gq(:),g2(:),yl(:),ylw(:),he(:),hr(:)
      complex(8), allocatable :: optl(:,:,:)
C     Local for parallel threads
      integer,allocatable:: ivp(:)
      complex(8),allocatable::psi(:,:,:),psir(:,:,:,:),vpsi(:,:,:),wk(:,:,:)
      real(8),allocatable:: cosi(:),sini(:),wk2(:)
C ... Local parameters
      integer procid, nproc
      integer n0,nkap0,nermx,npmx,nblk,nlmto,loptic
      parameter (n0=10,nkap0=4,nermx=100,npmx=128)
      integer nspec,ngabc(3),n1,n2,n3,nrt,net,ng,ltop,nlmtop
      integer nfilo,nfiup,nemlo,nemup,nfilm,nempm
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
      double precision alat,qlat(3,3),plat(3,3),q0(3),gmax,xx
      procedure(integer) :: iprint,nglob,mpipid
C     Shared variables (openmp)
      integer ivec,nvec
      double precision vol,etab(nermx),rtab(nermx)

      if (nevec <= 0) return

C ... First setup
      call tcn('rsibl')
      nproc  = mpipid(0)
      procid = mpipid(1)
      nbas  = nglob('nbas')
      nspec = nglob('nspec')
      nlmto = ndimh-napw
      loptic = s_optic%loptic  ! Set loptic=0 to suppress interstital part of optical matrix element
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      gmax = s_lat%gmax
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol

C ... Setup for q-dependent gv ... also makes kv, gv+q and iv
      call pshpr(iprint()-30)
      call sugvec(s_lat,5009,q,0d0,gmax,ngabc,0,ng)
      call poppr

C     For PW basis ... for now.
      if (napw > 0) then
        allocate(ivp(napw))
        call gvgvcomp(ng,s_lat%igv,napw,igapw,ivp)
      else
        allocate(ivp(1))
      endif

C --- Tables of energies, rsm, indices to them ---
      call tbhsi(s_spec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)

C --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2
      allocate(gq(ng*3))
      allocate(g2(ng))
      allocate(yl(ng*nlmtop))
      allocate(he(ng*net))
      allocate(hr(ng*nrt))
      call dpzero(q0,3)
      if (nlmto > 0) then
      call hsibq1(net,etab,nrt,rtab,ltop,alat,q0,ng,s_lat%gv,gq,
     .  g2,yl,he,hr)
      endif

      if (loptic > 0) then  ! for now
        nblk = nevec
      else
        nblk = 16
      endif

C  --- Loop over eigenstates ---
C$OMP PARALLEL PRIVATE (psi,vpsi,psir,cosi,sini,wk,wk2,ylw)
C$OMP& REDUCTION(+:smrho,f)
      allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
      allocate(psir(k1,k2,k3,nspc),cosi(ng),sini(ng),wk2(ng))
      allocate(ylw(ng*nlmtop))
C$OMP DO schedule(static,1) PRIVATE (ivec,nvec)
      do  ivec = 1, nevec, nblk
        nvec = min(nblk, nevec-ivec+1)

C   ... Add together total smooth wavefunction
C       i = 0; if (loptic > 0) i = 2
        call rsibl1(0,s_site,s_spec,q,nbas,iprmb,ng,gq,s_lat%igv,n1,n2,
     .    n3,qlat,cosi,sini,yl,ylw,he,hr,wk,wk2,vol,iprt,ipet,etab,rtab,
     .    ndimh,nlmto,nspc,numq,ewgt,ivec,nvec,evec,xx,psi,xx)
C   ... Add smooth contribution to optical matrix element
        if (loptic > 0) then
          call opt_nstate(s_optic,3,nevec,nfilo,nfiup,nemlo,nemup,nfilm,nempm)
          if (size(s_optic%optme,1) /= nfilm .or.
     .        size(s_optic%optme,2) /= nempm .or.
     .        size(s_optic%optme,3) /= 3)
     .      call rx('rsibl: array mismatch')
          if (nfilm*nempm > 0) then  !Do only if some transitions in window
          allocate(optl(nfilo:nfiup,nemlo:nemup,3))
          call rsibl8(vol,ng,nspc,nfilo,nfiup,nemlo,nemup,gq,psi,optl)
C         call zprm('optme(istl)',2,optl,nfilm,nfilm,nempm)
          call daxpy(2*nfilm*nempm*3,1d0,optl,1,s_optic%optme,1)
C         call zprm('optme(1-2+istl)',2,s_optic%optme,nfilm,nfilm,nempm)
          deallocate(optl)
          endif
        endif
        call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol),
     .    evec(1,1,ivec),psi)

C   ... Add to real-space mesh, optionally make smpot*psi for forces
        call rsibl2(ng,nspc,nvec,psi,n1,n2,n3,k1,k2,k3,s_lat%kv,numq,
     .    ewgt(1,ivec),mod(lfrce,10),smpot(1,1,1,isp),
     .    psir,smrho(1,1,1,1,isp,1),vpsi)

C    --- Add to forces, or calculate optical matrix elements ---
        if (mod(lfrce,10) /= 0) then
          call rsibl1(1,s_site,s_spec,q,nbas,iprmb,ng,gq,s_lat%igv,
     .      n1,n2,n3,qlat,cosi,sini,yl,ylw,he,hr,wk,wk2,vol,iprt,ipet,
     .      etab,rtab,ndimh,nlmto,nspc,numq,ewgt,ivec,nvec,evec,vpsi,
     .      psi,f)
        endif

      enddo
C$OMP END DO
      deallocate(psi,vpsi,wk,psir,cosi,sini,wk2,ylw)
C$OMP END PARALLEL

      deallocate(gq,g2,yl,he,hr)
      deallocate(ivp)
C     call zprm3('smrho after rsibl',0,smrho,k1,k2,k3)
      call tcx('rsibl')
      end

      subroutine rsibl1(mode,s_site,s_spec,q,nbas,iprmb,ng,gq,iv,n1,n2,
     .  n3,qlat,cosgp,singp,yl,ylw,he,hr,psi0,wk2,vol,iprt,ipet,etab,
     .  rtab,ndimh,nlmto,nspc,numq,ewgt,ivec,nvec,evec,vpsi,psi,f)
C- Make wave function for a block of evecs, or add contr. to forces
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  ngcut lmxa lmxb pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Ci Inputs
Ci   mode  :0 make wave function
Ci         :1 Add 2*Re( (v psi+) grad(psi) ) to f
Ci   q     :Bloch wave number
Ci   nbas  :size of basis
Ci   ng    :number of G-vectors
Ci   gq    :2*pi/alat * (q+G) for all G-vectors
Ci   iv    :g-vectors as integer multiples of qlat (suphs0)
Ci   n1..3 :size uniform mesh for smooth density and potential
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   cosgp :cos(phase) for each g-vector
Ci   singp :sin(phase) for each g-vector
Ci   yl    :spherical harmonics for ng vectors
Ci   ylw   :work array of same dimension as yl
Ci   he    :table of energy factors
Ci   hr    :table of smoothing radius factors
Ci   psi0  :work array (dim ng*2*nspc*nev): psi sans phase factors
Ci   wk2   :work array of dimension ng
Ci   vol   :cell volume
Co   iprt  :index to which entry in rt a given orbital belongs
Ci   ipet  :index to which entry in etab a given orbital belongs
Ci   etab  :table of all inequivalent energies
Ci   rtab  :table of all inequivalent smoothing radii
Ci   ndimh :dimensions evec
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   numq  :number of trial fermi levels
Ci   ewgt  :weights for each of the trial fermi levels
Ci   ivec  :first of current block of eigenvectors
Ci   nvec  :number of eigenstates to generate
Ci   evec  :eigenvectors
Ci   vspi  :potential * wave function, needed only for mode=1
Co Outputs
Co   psi   :wave function (mode=0); work area (mode=1)
Co   f     :term added to forces (mode=1)
Cr Remarks
Cu Updates
Cu   14 Jan 14 (B. Kaube) option for calculating optmt
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Apr 09 (T. Kotani) Added call to ncutcorrect
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision q(3)
      integer mode,nbas,ng,ndimh,nlmto,nspc,numq,ivec,nvec,iv(ng,3),
     .  n1,n2,n3,n0,nkap0,iprmb(*)
      parameter (n0=10,nkap0=4)
      integer iprt(n0,nkap0,*),ipet(n0,nkap0,*)
      double precision vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*),
     .  psi0(ng,2,nspc,*),wk2(ng),cosgp(ng),singp(ng),etab(*),
     .  rtab(*),gq(ng,3),f(3,nbas,numq),ewgt(numq,*),qlat(3,3)
      double complex psi(ng,nspc,*),evec(ndimh,nspc,*),vpsi(ng,nspc,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
C ... Local parameters
      logical lmv6,cmdopt
      character*40 strn
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .  blks(n0*nkap0),ntab(n0*nkap0),ncut(n0,nkap0),lh(nkap0),nkapi
      double precision e,rsm,eh(n0,nkap0),rsmh(n0,nkap0),f0(3)
      double precision xx(n0),wt,p(3)
      integer ib,is,io,jo,l2,kp,ie,ir,ioff,nlm1,nlm2,iq,kb,lt,i

      call dpzero(psi, 2*ng*nspc*nvec)

      if (nlmto == 0) return

      lmv6 = cmdopt('--lmv6',6,0,strn)

      do  ib = 1, nbas
        is = s_site(ib)%spec
        p = s_site(ib)%pos
        ncut = s_spec(is)%ngcut
        call suphas(q,p,ng,iv,n1,n2,n3,qlat,cosgp,singp)
C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C       Block into groups with consecutive l and common (e,rsm)
        call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkapi)
        call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)

        call dpzero(psi0,ng*2*nspc*nvec)
        if (mode == 1) call dpzero(psi, 2*ng*nspc*nvec)
        do  io = 1, norb
        if (blks(io) /= 0) then
          jo = ntab(io)
          l2 = ltab(io)
          lt = ltab(jo)
          kp = ktab(io)
          ie = ipet(l2+1,kp,is)
          ir = iprt(l2+1,kp,is)
          ioff = offl(io)
          nlm1 = l2**2+1
          nlm2 = nlm1 + blks(io)-1
          rsm = rtab(ir)
          e   = etab(ie)
          if (.not.lmv6) call ncutcorrect(ncut(lt+1,kp),1,gq,ng) ! TK April 2009
          call rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,min(ng,ncut(lt+1,kp))
     .      ,yl,ylw,he,hr,wk2,ioff,evec(1,1,ivec),ndimh,nspc,nvec,psi0)
        endif
        enddo
        call rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)

C   --- Forces --
        if (mode == 1) then
          do  i = 1, nvec
            call rsibl4(vol,ng,nspc,gq,vpsi(1,1,i),psi(1,1,i),f0)
            do  iq = 1, numq
              wt = ewgt(iq,i+ivec-1)
              f(1,ib,iq) = f(1,ib,iq) + wt*f0(1)
              f(2,ib,iq) = f(2,ib,iq) + wt*f0(2)
              f(3,ib,iq) = f(3,ib,iq) + wt*f0(3)
C             This shouldn't be necessary
              do  kb = 1, nbas
                f(1,kb,iq) = f(1,kb,iq) - wt*f0(1)/nbas
                f(2,kb,iq) = f(2,kb,iq) - wt*f0(2)/nbas
                f(3,kb,iq) = f(3,kb,iq) - wt*f0(3)/nbas
              enddo
            enddo
          enddo
        endif
      enddo

      end

      subroutine rsibl2(ng,nspc,nev,psi,n1,n2,n3,k1,k2,k3,kv,numq,ewgt,
     .  lfrce,smpot,psir,smrho,vpsi)
C- FT wave function to real space and add square into mesh density
C  and optionally make smpot * psi
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nev   :number of wave functions
Ci   psi   :wave function in reciprocal space
Ci   n1..3 :size of FT mesh
Ci   k1..3 :dimensions smpot,smrho
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   numq  :number of trial fermi levels
Ci   ewgt  :weights for each of the trial fermi levels
Ci   lfrce :if nonzero, make vpsi
Ci   smpot :potential, needed if lfrce is nonzero
Co Outputs
Co   psir  :psi in real space
Co   smrho :ewgt (psir+)psir added to smooth density
Co   vpsi  :FT (smpot * r.s. wave function) if lfrce is nonzero
Cr Remarks
Cu Updates
Cu   09 Apr 10 First attempt at noncollinear smrho.
Cu             Not checked; force terms not implemented
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,ng,nspc,nev,numq,kv(ng,3),lfrce
      double precision ewgt(numq,nev)
      double complex psi(ng,nspc,nev),vpsi(ng,nspc,nev),
     .  psir(k1,k2,k3,nspc)
      double complex smrho(k1,k2,k3,numq,nspc*nspc),smpot(k1,k2,k3,nspc*nspc)
C ... Local parameters
      integer i1,i2,i3,iq,i,ispc
      double precision wgt1

      call tcn('rsibl2')
      do  ispc = 1, nspc
      do  i = 1, nev
      call gvputf(ng,1,kv,k1,k2,k3,psi(1,ispc,i),psir)
      call fftz3(psir,n1,n2,n3,k1,k2,k3,1,0,1)
      do  iq = 1, numq
        wgt1 = ewgt(iq,i)
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              smrho(i1,i2,i3,iq,ispc) = smrho(i1,i2,i3,iq,ispc) +
     .           wgt1*dconjg(psir(i1,i2,i3,1))*psir(i1,i2,i3,1)
            enddo
          enddo
        enddo
      enddo
C ... Make the spin off-diagonal parts of smrho
      if (nspc == 2) then
        call gvputf(ng,1,kv,k1,k2,k3,psi(1,3-ispc,i),psir(1,1,1,2))
        call fftz3(psir(1,1,1,2),n1,n2,n3,k1,k2,k3,1,0,1)
        do  iq = 1, numq
          wgt1 = ewgt(iq,i)
          do  i3 = 1, n3
            do  i2 = 1, n2
              do  i1 = 1, n1
                smrho(i1,i2,i3,iq,ispc+2) = smrho(i1,i2,i3,iq,ispc+2) +
     .             wgt1*dconjg(psir(i1,i2,i3,1))*psir(i1,i2,i3,2)
              enddo
            enddo
          enddo
        enddo
      endif

      if (lfrce /= 0) then
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              psir(i1,i2,i3,1) = psir(i1,i2,i3,1)*smpot(i1,i2,i3,ispc)
            enddo
          enddo
        enddo
        call fftz3(psir,n1,n2,n3,k1,k2,k3,1,0,-1)
        call gvgetf(ng,1,kv,k1,k2,k3,psir,vpsi(1,ispc,i))
      endif
      enddo
      enddo
      call tcx('rsibl2')
      end

C      subroutine rsibl3(ng,n1,n2,n3,k1,k2,k3,kv,smpot,f,vpsi)
CC- Make f*smpot and transform to reciprocal space
C      implicit none
CC ... Passed parameters
C      integer k1,k2,k3,n1,n2,n3,ng,kv(ng,3)
C      double complex vpsi(ng),f(k1,k2,k3),smpot(k1,k2,k3)
CC ... Local parameters
C      integer i1,i2,i3
C
C      call tcn('rsibl3')
C      do  i3 = 1, n3
C        do  i2 = 1, n2
C          do  i1 = 1, n1
C            f(i1,i2,i3) = f(i1,i2,i3)*smpot(i1,i2,i3)
C          enddo
C        enddo
C      enddo
C      call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,kv,k1,k2,k3,f,vpsi)
C      call tcx('rsibl3')
C      end

      subroutine rsibl4(vol,ng,nspc,gq,vpsi,psi,f0)
C- Force term 2*Re( (psi_nu+) vsm grad(psi_nu) )
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   gq    :2*pi/alat * (q+G) for all G-vectors
Ci   vpsi  :(psi vsm) in reciprocal space
Ci   psi   :portion of wave function associated with one site ib
Co Outputs
Co   f0    :2*Re( (psi_nu+) vsm grad(psi_ib,nu) ) put in f0
Cr Remarks
Cr   gradient operator is i*G
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nspc
      double precision vol,gq(ng,3),f0(3)
      double precision vpsi(2,nspc,ng),psi(2,nspc,ng)
C ... Local parameters
      integer ispc,ig
      double precision sum1,sum2,sum3,xx
C     double complex cc,ovl

      call tcn('rsibl4')
      sum1 = 0
      sum2 = 0
      sum3 = 0
C     ovl = 0
      do  ispc = 1, nspc
      do  ig = 1, ng
C        cc  = dcmplx(psi(1,ispc,ig),-psi(2,ispc,ig))
C     .       *dcmplx(vpsi(1,ispc,ig),vpsi(2,ispc,ig)) * vol
C        ovl = ovl + cc
C     .      * dcmplx(0d0,1d0)*gq(ig,1)
        xx = vpsi(2,ispc,ig)*psi(1,ispc,ig)
     .     - vpsi(1,ispc,ig)*psi(2,ispc,ig)
        sum1 = sum1 + gq(ig,1)*xx
        sum2 = sum2 + gq(ig,2)*xx
        sum3 = sum3 + gq(ig,3)*xx
      enddo
      enddo

      f0(1) = 2*vol*sum1
      f0(2) = 2*vol*sum2
      f0(3) = 2*vol*sum3

C     print *, ovl

      call tcx('rsibl4')

      end

      subroutine rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,ncut,yl,ylw,he,hr,
     .  wk,ioff,evec,ndimh,nspc,nvec,psi)
C- Add contribution to wave function from one block of orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci   ie    :index to appropriate entry in energy factor table he
Ci   ir    :index to appropriate entry in sm radius factor table hr
Ci   e     :hankel energy
Ci   rsm   :smoothing radius
Ci   vol   :cell volume
Ci   nlm1  :starting orbital L for which to accumulate wave function
Ci   nlm2  :final orbital L for which to accumulate wave function
Ci   ng    :number of G-vectors
Ci   ncut  :G-cutoff for wave function
Ci   yl    :spherical harmonics for ng vectors
Ci   ylw   :work array dimensioned same as yl
Ci   he    :table of energy factors
Ci   hr    :table of smoothing radius factors
Ci   wk    :work array of dimension at least ncut
Ci   ioff  :offset to hamiltonian (eigenvector) for this orbital block
Ci   ndimh :dimension of hamiltonian
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   evec  :eigenvectors
Ci   nvec  :number of eigenvectors
Co Outputs
Co   psi   :contribution to psi from this block accumulated
Cr Remarks
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ie,ioff,ir,ncut,ng,nlm1,nlm2,ndimh,nspc,nvec
      double precision e,rsm,vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*),
     .  wk(ncut)
      double complex evec(ndimh,nspc,nvec),psi(ng,nspc,nvec)
C ... Local parameters
      integer i,ii,ilm,l,ll,lmax,m,iv,nlmx,ispc
      parameter (nlmx=100)
      double complex cfac(nlmx),cc,evp(nlmx,nspc,nvec)
      double precision pi,fac
      parameter (pi=3.1415926535897931d0)

      if (nlm2 == 0) return
      call tcn('rsibl5')

C ... Phase and other factors
      lmax = ll(nlm2)
      fac = -4d0*pi*dexp(e*rsm*rsm*0.25d0)/vol
      cc = (0d0,1d0)*fac
      ilm = 0
      do  l = 0, lmax
        cc = cc*(0d0,-1d0)
        do m = -l,l
          ilm = ilm+1
          cfac(ilm) = cc
        enddo
      enddo

C ... Combine G-dependent energy, rsm and YL factors
      do  i = 1, ncut
        wk(i) = he(i,ie)*hr(i,ir)
      enddo
      do  ilm = nlm1, nlm2
      do  i = 1, ncut
        ylw(i,ilm) = wk(i)*yl(i,ilm)
      enddo
      enddo

C ... Make vector evec*phase
      do  ispc = 1, nspc
      do  ilm = nlm1, nlm2
        ii = ilm-nlm1+ioff+1
        do  iv = 1, nvec
C          cc = evec(ii,ispc,iv)*cfac(ilm)
C          evpr(ilm,ispc,iv) = dble(cc)
C          evpi(ilm,ispc,iv) = dimag(cc)
          evp(ilm,ispc,iv) = evec(ii,ispc,iv)*cfac(ilm)
        enddo
      enddo

C ... For each orbital and evecs 1..nvec, accumulate psi
C     ii = ilm-nlm1+ioff+1
      do  ilm = nlm1, nlm2
        do  iv = 1, nvec
          do  i = 1, ncut
            psi(i,ispc,iv) = psi(i,ispc,iv)+ylw(i,ilm)*evp(ilm,ispc,iv)
          enddo
        enddo
      enddo
      enddo

      call tcx('rsibl5')
      end

      subroutine rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)
C- Multiply by phase to make final FT of partial wave function
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nvec  :number of eigenvectors
Ci   cosgp :cos(phase)
Ci   singp :sin(phase)
Ci   wr    :real part of psi, unscaled by phase
Ci   wi    :imaginary part of psi, unscaled by phase
Co Outputs
Co   psi   :(wr,si)*(cosgp,singp) is added into psi
Cr Remarks
Cu Updates
Cu   23 Dec 04 Extended to noncollinear case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nspc,nvec
      double precision cosgp(ng),singp(ng)
      double complex psi0(ng,nspc,nvec),psi(ng,nspc,nvec)
C ... Local parameters
      integer i,iv,ispc

      call tcn('rsibl6')

      do  iv = 1, nvec
      do  ispc = 1, nspc
      do  i = 1, ng
        psi(i,ispc,iv) = psi(i,ispc,iv)
     .                 + psi0(i,ispc,iv)*dcmplx(cosgp(i),singp(i))
      enddo
      enddo
      enddo

      call tcx('rsibl6')
      end

      subroutine rsibl7(xsmrho,k1,k2,k3,numq,smrho)
C- Combine rho from separate parallel threads
      implicit none
C ... Passed parameters
      integer k1,k2,k3,numq
      double complex smrho(k1,k2,k3,numq),xsmrho(k1,k2,k3,numq)
C ... Local parameters
      integer ik1,ik2,ik3,iq
      do  iq = 1, numq
        do  ik3 = 1, k3
          do  ik2 = 1, k2
            do  ik1 = 1, k1
              smrho(ik1,ik2,ik3,iq) = smrho(ik1,ik2,ik3,iq) +
     .                               xsmrho(ik1,ik2,ik3,iq)
            enddo
          enddo
        enddo
      enddo

      end

      subroutine rsibl8(vol,ng,nspc,nfilo,nfiup,nemlo,nemup,gq,psi,opt)
C- Matrix element of the gradient operator (psi_i)+ grad grad(psi_j)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nfilo :Matrix elements for occupied bands nfilo, nfiup
Ci   nfiup :Matrix elements for occupied bands nfilo, nfiup
Ci   nemlo :Matrix elements for unoccupied bands nemlo, nemup
Ci   nemup :Matrix elements for unoccupied bands nemlo, nemup
Ci   gq    :2*pi/alat * (q+G) for all G-vectors
Ci   psi   :wave function on a mesh of points (G+k)
Co Outputs
Co   opt   :(psi_i)+ grad_(1,2,3) grad(psi_j) accumulated into opt(i,j,1:3)
Cr Remarks
Cr   gradient operator is i*(G+q) in k-space
Cu Updates
Cu   04 Jun 14 first cut at the noncollinear case
Cu   06 Feb 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nspc
      integer nfilo,nfiup,nemlo,nemup
      complex(8) :: psi(ng,nspc,nemup),opt(nfilo:nfiup,nemlo:nemup,3)
      real(8) :: vol,gq(ng,3)
C ... Local parameters
      complex(8), allocatable :: psioc(:,:,:),psiun(:,:,:)
      integer ig,i,k,j,isp,nfilm,nempm
C     double complex zdotc

      call tcn('rsibl8')

      allocate(psioc(ng,nspc,nfilo:nfiup),psiun(ng,nspc,nemlo:nemup))
      nfilm = nfiup-nfilo+1
      nempm = nemup-nemlo+1

C ... Set up <psi_i|
      do  isp = 1, nspc
      do  i = nfilo, nfiup
        call zcopy(ng,psi(1,isp,i),1,psioc(1,isp,i),1)
      enddo
      enddo
C      call zprm('psi',2,psi,ng*nspc,ng*nspc,nfilm)
C      call prmx('G+q',gq,ng,ng,3)

C --- For each x,y,z, do ---
      do  k = 1, 3

C   ... Set up i(G+q)|psi_j>
        do  isp = 1, nspc
        do  j = nemlo, nemup
          do  ig = 1, ng
            psiun(ig,isp,j) = dcmplx(0d0,gq(ig,k))*psi(ig,isp,j)
C           psiun(ig,isp,j) = dcmplx(1d0,0d0)*psi(ig,isp,j)
          enddo
        enddo
        enddo

C   ... <i | i(G+q) | j >
        call zgemm('C','N',nfilm,nempm,ng*nspc,(1d0,0d0),psioc,ng*nspc,
     .    psiun,ng*nspc,(0d0,0d0),opt(nfilo,nemlo,k),nfilm)

      enddo ! Loop over x,y,z

C     Debugging
c     mch -qr psi -a p p gq -coll 1 -v2dia -s0,1 -a gq psi -cc -t gq -x psi -x out.ag  -- -px
C     call zprm('opt(x)',2,opt(nfilo,nemlo,1),nfilm,nfilm,nempm)
C     call zprm('opt(y)',2,opt(nfilo,nemlo,2),nfilm,nfilm,nempm)
C     call zprm('opt(z)',2,opt(nfilo,nemlo,3),nfilm,nfilm,nempm)

C     Psi defined as 1/sqrt(vol) * exp(..)
      call dscal(2*nfilm*nempm*3,vol,opt,1)

      deallocate(psioc,psiun)
      call tcx('rsibl8')

      end

      subroutine gvgvcomp(ng,igv,napw,igapw,ivp)
C- Find pointer ivp that maps igv to igapw.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Co   igv   :list of reciprocal lattice vectors G (gvlist.f)
Ci   napw   :number of R.L.V for PW basis
Ci   igapw :reciprocal lattice vectors for PW basis.
Co Outputs
Co   ivp   :if ig = ivp(jg), igv(jg) and nvec(ig) are same vector
Cr Remarks
Cr  This routine should be be cleaned up and ivp
Cr  used by rest of program in place of igapw
Cu Updates
Cu   05 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,igv(ng,3),napw,igapw(3,napw),ivp(napw)
C ... Local parameters
      integer jg,ig

C     Redesign ... inefficient.
      do  ig = 1, napw
        do  jg = 1, ng
          if(sum(abs(igv(jg,:)-igapw(:,ig))) == 0) then
            ivp(ig) = jg
            goto 333
          endif
        enddo
        ivp(ig) = -9999
        call rx('gvgvcomp wrong 111! maybe enlarge GMAX or so')
  333   continue
        if( sum(abs( igapw(:,ig)-igv(ivp(ig),:) )) /=0)
     .    call rx('bug in gvgvcomp.  Cannot find ivp')
      enddo
      end

      subroutine rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,sqv,evec,psi)
C- Plane wave part of evec
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nvec  :number of eigenstates to generate
Ci   evec  :eigenvectors
Ci   vspi  :potential * wave function, needed only for mode=1
Ci   sqv   :square root of volume
Co Outputs
Co   psi   :wave function
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,ndimh,nlmto,nspc,nvec
      integer napw,ivp(napw)
      double precision sqv
      double complex psi(ng,nspc,nvec),evec(ndimh,nspc,nvec)
C ... Local parameters
      integer i,ispc,igv

      if (napw <= 0) return

      do  ispc = 1, nspc
        do  i = 1, nvec
          do  igv = 1, napw
            psi(ivp(igv),ispc,i) = psi(ivp(igv),ispc,i)
     .        + evec(nlmto+igv,ispc,i)/sqv
          enddo
        enddo
      enddo

      end
