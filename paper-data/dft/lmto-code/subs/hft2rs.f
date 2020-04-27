      subroutine hft2rs(mode,s_ham,n1,n2,n3,k1,k2,k3,qoff,isp,nsp,nbas,
     .  g,ag,nsgrp,range,ia1,ia2,idim,jdim,hq,plat,pos,ndhrs,hrs)
C- Fourier transform of a matrix to real-space representation
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :See also Remarks
Ci           :1s digit specifies what is created.
Ci           :  0 (setup) create tables iax and ntab
Ci           :  1 (FT)    convert hq to h(T) by FFT
Ci           :  These modes are distinct so the caller can
Ci           :  allocate the appropriate memory for hrs.
Ci           :  Note: for converting hq to h(T) by standard FT,
Ci              see Remarks
Ci           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Ci           :  1 Assume hrs is real (correct if hq is hermitian)
Ci           :100s digit specifies method for FT:
Ci           :  0 convert hq to h(T) by FFT.
Ci           :    hq is supplied on a uniform mesh of q-points
Ci           :    It is converted to R.S. by fft and copied to hrs.
Ci           :  1 convert hq to h(T) standard inverse Bloch transform
Ci           :    Note: hft2rs only does the setup for the
Ci           :    standard inverse Bloch transform; see Remarks
Ci           :1000s digit concerns generation of neighbor table:
Ci           :      (applicable only when 1s digit mode is zero)
Ci           :      See Note 6 in Remarks.
Ci           :  0 table is not altered by symmetry considerations
Ci           :  1 table is enlarged so that any pair connecting ri-rj
Ci           :    has a corresponding entry rj-ri.
Ci           :  2 table is reduced so that any pair connecting ri-rj
Ci           :    has a corresponding entry rj-ri.
Ci           :10000s digit specifies indexing of hrs_RL,R'L'
Ci           :  0 poke h_RL,R'L' (T) into hrs_RL,R'L'
Ci           :  1 poke h_RL,R'L' (T) into hrs_R'L',RL
Ci           :    (used when h(T) constructed for rj-ri rather than ri-rj)
Ci   n1,..n3 :number divisions in QP mesh
Ci   k1,..k3 :dimensions hq
Ci   qoff    :offset of gamma-point = q(1,1,1)
Ci   isp     :current spin channel (1 or 2): poke hq into hrs(isp)
Ci           :NB: hq has no spin index
Ci   nsp     :2 for spin-polarized case, otherwise 1
Ci   nbas    :size of basis
Ci   g       :point group operations
Ci   ag      :translation part of space group
Ci   nsgrp   :number of group operations
Ci   range   :maximum length for connecting vectors
Ci   idim    :hamiltonian dimension of ia1..ia2
Ci   jdim    :hamiltonian dimension of ib1..ib2
Ci   hq      :hamiltonian on uniform mesh of q-points
Ci           :or translation vectors (related to q-points by FT)
Ci   plat    :primitive lattice vectors, in units of alat
Ci   pos     :basis vectors
Ci   ndhrs   :leading dimensions of hrs: must be at least as large as
Ci           :the total number of orbitals of any atom.
Cio Inputs/Outputs
Cio  s_ham%nprs: offset ntab(ib) to neighbor table for cluster ib (pairc.f)
Cio           If 1s digit mode=0, memory is allocated for ntab and array
Cio           is created.  If 1s digit mode>1, ntab is an input.
Cio  s_ham%iaxs: neighbor table iax containing pair information (pairc.f)
Cio           If 1s digit mode=0, memory is allocated for iax and array
Cio           is created.  If 1s digit mode>1, iax is an input.
Co   hrs     :real-space h in neighbor-list format (see bloch.f)
Co           :hrs is created when 1s digit mode>0 only.
Cr Remarks
Cr   hft2rs is designed to generate a real-space h_ij(T) from h_ij(q).
Cr   h_ij(T) is kept in standard tight-binding format, which entails
Cr   the following:
Cr     1. a neighbor table iax that contains a list of all connecting
Cr        vectors connecting sites i and j for the r.s. hq.
Cr     2. h_ij(T) is stored as a list of matrices; one matrix
Cr        connects all orbitals at site i to all orbitals at site j+T.
Cr     3. the forward Bloch transform of h_ij(T) should generate
Cr        h_ij (q); this routine generates h_ij(T) in the form that
Cr        routine Bloch will make h_ij(q) from h_ij(T).
Cr     4. The inverse Bloch transform may be done by FFT or by
Cr        a standard transform.  In the former case, the forward
Cr        transform from the hq(T) computed here should be exact on
Cr        the mesh of q-points from which hq(T) was made.
Cr        For the standard inverse transform it need not.
Cr     5. This routine is called once in a setup mode to create the
Cr        pair table (1s digit mode=0).  For the inverse transform by
Cr        FFT it is called a second time to do the mapping.  hft2rs will
Cr        not do the inverse transform using the standard approach,
Cr        because of complications connected with symmetrization.
Cr
Cr        For FFT, you will want to use the following procedure:
Cr          a) call hft2rs in setup mode (1s digit mode = 0)
Cr             Note: after this call, you may want to adjust the
Cr             neigbor table to accomodate symmetry operations;
Cr             see Note 7 below.
Cr          b) Create hq(q) in the full BZ.  If you have only hq(q)
Cr             for irreducible points, see hamfb3.f
Cr          c) call hft2rs again, with 1s digit mode = 1
Cr             At this point the FT is complete.  However it may
Cr             have someone undesirable properties as regards
Cr             symmetry. You may want to symmetrize the result;
Cr             see Note 6 below.
Cr        The standard inverse transform differs in the following ways:
Cr          * the inverse transform need to be exact
Cr          * the inverse transform uses hq(q) at the irreducible
Cr            qpoints
Cr        For this method use the following procedure:
Cr          a) call hft2rs in setup mode (1s digit mode = 0)
Cr             Note: after this call, you will want to adjust the
Cr             neigbor table to accomodate symmetry operations;
Cr             see Note 6 below.
Cr          b) make an unsymmetrized inverse Bloch transforming
Cr             by calling ibloch for each irreducible qp.
Cr          c) Symmetrize the inverse transform by calling rsmsym.
Cr
Cr   6.  Symmetry.  The pair table may be at variance with the symmetry
Cr       of the system.  This can happen for two reasons: points that
Cr       are symmetry-related may not all find their way into the pair
Cr       table because their connecting vector is within a tight
Cr       tolerance of the r.s. cutoff radius, or in the case of the FFT,
Cr       the list of points is reduce to make a list whose size
Cr       corresponds to the mesh of q-points.  You can cause hft2rs to
Cr       make the table compatible with symmetry operations by setting
Cr       the 1000s digit of mode in the setup call.  The symmetry
Cr       operations are passed as the triplet (g,ag,nsgrp)
Cr       hft2rs actually changes the table in two ways:
Cr
Cr       * by adding to the table all missing vectors that are
Cr         symmetry-related to vectors already in the table
Cr
Cr       * by adding to the table all vectors not in the table that
Cr         correspond to the inverse of vectors in the table.
Cr
Cr       Thus, even if you tell hft2rs to change the table with only
Cr       on group operation, it may affect the table.
Cr
Cr       Note, however, the FFT inverse will only poke hq(q) into
Cr       hq(T) for mesh of T it chooses.  If you want the transform
Cr       to be symmetrize, call rsmsym to do it.
Cr
Cr   Notes on the meaning of row and column dimensions:
Cr     ib=field(hq=gf)   or augmentation(hq=ham)
Cr     jb=source(hif=gf) or basis(hq=ham)
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   25 Sep 04 Handles offset BZ mesh
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu    3 Oct 02 redesigned
Cu   23 Jun 02 First written
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,idim,jdim,nbas,n1,n2,n3,k1,k2,k3,isp,nsp,ndhrs,nsgrp
      integer ia1,ia2
      integer nkap0,n0H,niax
      parameter (nkap0=4,n0H=5,niax=10)
      double precision plat(3,3),pos(3,*),range,g(9,*),ag(3,*),qoff(3)
      double complex hq(k1,k2,k3,idim,jdim,nsp),hrs(ndhrs,ndhrs,nsp,1)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
C#ifndef LINUXF
      integer, allocatable, target :: iax0(:)
C#elseC
C      integer, pointer :: iax0(:) => null()
C#endif
      integer, allocatable :: iwk(:)
      integer, pointer :: iax(:)
C ... Local parameters
      integer i,ipr,iset,mode0,mode1,mode2,mode3,mxcsiz,nttab,nviax,
     .  nvtot
      double precision dmx(2)
C#ifdefC DEBUG
C      integer iaxi(5)
C      double precision tau(3)
C#endif
C ... External calls
      external dinv33,dpsadd,dpscop,fftz3s,getpr,gvlist,hft2r1,hft2r2,
     .         hft2r3,icopy,info5,isanrg,ivheap,ivprm,poppr,
     .         ppair1,pshpr,ptr_ham,rxi,symiax,tcn,tcx

      call tcn('hft2rs')
      call getpr(ipr)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      call sanrg(.true.,mode1,0,1,'hft2rs:','10s digit mode')
      call sanrg(.true.,mode2,0,1,'hft2rs:','100s digit mode')

C --- Mode 0 (setup) make and reduce iax table ---
      if (mode0 == 0) then

        call info(30,0,0,' hft2rs: make neighbor table for '//
     .    'r.s. hamiltonian using range = %;4d * alat',range,0)

C   ... Pair table for rs hamiltonian
        call pshpr(min(ipr-10,30))
        mxcsiz = 0
        call pairs(nbas,nbas,1d0,plat,[range/2],pos,
     .    [-1],3,-1,[0],nttab,s_ham%nprs,iax0,mxcsiz)
        call poppr

C   ... Standard Bloch transform doesn't require mapping to mesh
        if (mode2 == 0) then

C     ... Mark all members in iax table with vector in sfz (for FFT)
          nviax = 0
          nvtot = 0
          call hft2r1(mode,n1,n2,n3,k1,k2,k3,isp,nsp,ia1,ia2,1,nbas,
     .      s_ham%offH,idim,jdim,hq,plat,pos,s_ham%nprs,iax0,qoff,
     .      ndhrs,hrs,nviax,nvtot,dmx)
          call info(30,0,0,' hft2rs: found %i connecting vectors out'//
     .      ' of %i possible for FFT',nviax,nvtot)

C     ... Reduce iax table to elements on mesh in common with pairc
          allocate(iwk((niax+1)*mxcsiz))
          call hft2r3(1,nbas,iwk,s_ham%nprs,iax0,nttab,mxcsiz)
          deallocate(iwk)
        endif

C   ... Add pairs so that any ri-rj has corresponding rj-ri
        iax => iax0
        if (mode3 /= 0) then
          i = 110
          if (mode3 == 2) i = 111
          allocate(iax(niax*nttab*3))  ! maximum size of redimensioned iax
          call icopy(niax*nttab,iax0,1,iax,1) ! copy original
          deallocate(iax0) ! original no longer needed
          nttab = nttab*3
          call pshpr(ipr)
          call symiax(i,plat,nbas,pos,g,ag,nsgrp,s_ham%nprs,iax,
     .      nttab,mxcsiz)
          call poppr

C     ... Fill out iax(6,7,10)
          call pshpr(1)
          call ppair1(0,1,nbas,nbas,-1,1d0,plat,pos,range/2,
     .      nttab,s_ham%nprs,iax,mxcsiz)
          call poppr
        endif

C   ... iax is now complete ... copy array to strux
!       call ptr_ham(s_ham,4+1,'iaxs',niax*nttab,0,iax)
        if (associated(s_ham%iaxs)) deallocate(s_ham%iaxs)
        allocate(s_ham%iaxs(niax*nttab))
        s_ham%iaxs(:) = iax(1:niax*nttab)
        deallocate(iax)

C       Printout of pairs table
C       ipr = 70
C       call pshpr(ipr)
        if (ipr > 40) then
          call pshpr(ipr-10)
          call ppair1(20,1,nbas,nbas,-1,1d0,plat,pos,range/2,
     .    nttab,s_ham%nprs,s_ham%iaxs,mxcsiz)
          call poppr
        endif
        goto 999
      endif

C ... hq to real-space representation by FFT; copy hq to hrs
      if (mode0 == 1 .and. mode2 == 0) then
        iset = 0
        call fftz3s(n1,n2,n3,k1,k2,k3,iset)
C       call yprm('hfbz',3,hq,0,k1*k2*k3,k1*k2*k3,idim*jdim)
        call fftz3(hq,n1,n2,n3,k1,k2,k3,idim*jdim,iset,-1)
C       call yprm('hfbz(FT)',3,hq,0,k1*k2*k3,k1*k2*k3,idim*jdim)
C       Scale by phase for offset k-mesh
C        call hft2rp(0,lshft,n1,n2,n3,k1,k2,k3,idim,jdim,hq)
C        call yprm('hfbz(FT)+phase',3,hq,0,k1*k2*k3,k1*k2*k3,idim*jdim)

        call hft2r1(mode,n1,n2,n3,k1,k2,k3,isp,nsp,ia1,ia2,1,nbas,
     .    s_ham%offH,idim,jdim,hq,plat,pos,s_ham%nprs,s_ham%iaxs,qoff,
     .    ndhrs,hrs,nviax,nvtot,dmx)
        call info5(30,0,0,' hft2rs created hrs:  ndhrs=%i  '//
     .    'max Re(hrs) = %;3g  max Im(hrs) = %;3g',ndhrs,dmx,dmx(2),0,0)
      elseif (mode0 == 1 .and. mode2 == 1) then
        call rxi('hft2rs: not implemented',mode)
C        call ibloch(11-mode1,isp,nsp,w(oqp),w(owtkp),iq,plat,w(s_ham%offH),
C     .    w(oiprmb),1,nitab,s_ham%iaxs,ldim,ldim,w(oh2),w(osigr),nl**2)
      else
        call rxi('hft2rs: bad mode',mode)
      endif

  999 continue
      call tcx('hft2rs')

C#ifdefC DEBUG
C      if (mode0 == 0) return
C  998 print *, 'is?'
C      read(*,*) iset
C      if (iset == 0) return
C      call icopy(5,w(oiax-niax+niax*iset),1,iaxi,1)
C      do  i = 1, 3
C        tau(i) = pos(i,iaxi(2)) - pos(i,iaxi(1)) +
C     .    plat(i,1)*iaxi(3) + plat(i,2)*iaxi(4) + plat(i,3)*iaxi(5)
C      enddo
C      call info2(0,0,0,' iax(1:5,iset)=%5:2i %3:1,6;12D',iaxi,tau)
C      do  i = 1, min(ndhrs,9)
C        print 393, (dble(hrs(i,ipr,1,iset)), ipr=1,min(ndhrs,9))
C      enddo
C      print 393
C      do  i = 1, min(ndhrs,9)
C        print 393, (dimag(hrs(i,ipr,1,iset)), ipr=1,min(ndhrs,9))
C      enddo
C  393 format(9f12.6)
C      goto 998
C#endif

      end
      subroutine hft2r1(mode,n1,n2,n3,k1,k2,k3,isp,nsp,ia1,ia2,ib1,ib2,
     .  offH,idim,jdim,hij,plat,pos,ntab,iax,qoff,ndhrs,hrs,nviax,nvtot,
     .  dmx)
C- Copy a R.S. hamiltonian subblock on a uniform mesh of translation
C- vectors to a form specified by a neighbor list,
C  or just mark entries in the table that correspond to mesh points
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :1s digit specifies what is created.
Ci           :  0 marks entries in row 8 of iax table; see iax below
Ci           :    In this mode, hij and hrs are not used
Ci           : >0 copy hij to hrs
Ci           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Ci           :  1 Assume hrs is real (correct if hq is hermitian)
Ci           :10000s digit specifies indexing of hrs_RL,R'L'
Ci           :  0 poke h(T) into hrs_RL,R'L'
Ci           :  1 poke h(T) into hrs_R'L',RL
Ci           :    (used when h(T) constructed for rj-ri rather than ri-rj)
Ci   n1,..n3 :number divisions in QP mesh
Ci   isp     :current spin channel (1 or 2): poke hij into hrs(isp)
Ci           :NB: hij has no spin index
Ci   nsp     :2 for spin-polarized case, otherwise 1
Ci   ia1,ia2 :range of field point sites for which hij is tabulated
Ci   ib1,ib2 :range of source point sites for which hij is tabulated
Ci   offH    :Offsets to hamiltonian matrix (makidx.f)
Ci   idim    :hamiltonian dimension of ia1..ia2
Ci   jdim    :hamiltonian dimension of ib1..ib2
Ci   hij     :hamiltonian on uniform mesh of q-points
Ci   plat    :primitive lattice vectors, in units of alat
Ci   pos     :basis vectors
Ci   ntab    :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax     :neighbor table containing pair information (pairc.f)
Ci   ndhrs   :leading dimensions of hrs: must be at least as large as
Ci           :the total number of orbitals of any atom.
Co Outputs
Co   hrs     :diagonal GF for this principal layer
Co   nviax   :number of connecting vectors found in iax table added to nviax
Co   nvtot   :total number of connecting vectors added to nvtot
Co   dmx     :dmx(1) = largest real element found in h
Co           :dmx(2) = largest imaginary element found in h
Cr Remarks
Cr   ib=field(hij=gf)  or augmentation(hij=ham)
Cr   jb=source(hif=gf) or basis(hij=ham)
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   01 May 06 dmx is properly global maximum
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu   23 Jun 02 Adapted from gfq2rs.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,idim,jdim,ia1,ia2,ib1,ib2,n1,n2,n3,k1,k2,k3,isp,nsp,
     .  niax,ndhrs,nviax,nvtot
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5,niax=10)
      integer offH(n0H,nkap0,1),iax(niax,*),ntab(ib2+1)
      double precision plat(3,3),pos(3,*),dmx(2),qoff(3)
      double complex hij(k1,k2,k3,idim,jdim),hrs(ndhrs,ndhrs,nsp,1)
C ... Dynamically allocated local arrays
      integer, allocatable :: kv(:)
      real(8), allocatable :: gv(:)
      real(8), allocatable :: rtab(:)
C ... Local parameters
      logical ltrans
      integer id,jd,i,ia,ib,ngmx,ng,is1,is2,nlma,nlmb,offi,offj,
     .  nhit,nmiss
      double precision gmax,qlat(3,3),tau(3),dmx1(2),xx
C ... External calls
      external dinv33,dpsadd,dpscop,gvlist,hft2r2,isanrg,poppr,pshpr

      id = offH(1,1,ia2+1)-offH(1,1,ia1)
      call sanrg(.true.,idim,id,id,'hft2rs:','idim')
      jd = offH(1,1,ib2+1) - offH(1,1,ib1)
      call sanrg(.true.,jdim,jd,jd,'hft2rs:','jdim')
      ltrans = mod(mode/10000,10) /= 0
      dmx(1) = 0
      dmx(2) = 0

      offi = 0
      do  ia = ia1, ia2
      offj = 0
      do  ib = ib1, ib2

      is1 = ntab(ib)+1
      is2 = ntab(ib+1)

      nlma = offH(1,1,ia+1) - offH(1,1,ia)
      nlmb = offH(1,1,ib+1) - offH(1,1,ib)

C --- Make list for a dimensionless lattice (alat=1) ---
      call dinv33(plat,1,qlat,gmax)
      call dpscop(pos,tau,3,3*ib-2,1,1d0)
      call dpsadd(tau,pos,3,1,3*ia-2,-1d0)
c     call gvctof(0,1d0,qlat,tau,n1,n2,n3,gmax,ngmx)
      ngmx = n1*n2*n3
      allocate(gv(ngmx*3),kv(ngmx*3))
C     print *, 'ia,ib=',ia,ib
C     call ogvlst(1d0,qlat,tau,n1,n2,n3,9d9,ngmx,ng,gv,kv)
      call pshpr(80*0)
      i = 208 + 400*0
      call gvlist(1d0,qlat,tau,n1,n2,n3,1d5,i,ngmx,ng,kv,gv,xx,xx)
      call poppr

C --- For each connecting vector, mark iax table or poke into hrs ---
      allocate(rtab(3*(is2-is1+1)))
      if (ltrans) then
        call hft2r2(mode,k1,k2,k3,ia,ib,jdim,idim,iax,qoff,is1,is2,
     .    ndhrs,isp,nsp,nlmb,nlma,hij(1,1,1,1+offj,1+offi),gv,
     .    kv,ng,tau,plat,pos,rtab,hrs,hrs,nhit,nmiss,dmx1)
      else
        call hft2r2(mode,k1,k2,k3,ia,ib,idim,jdim,iax,qoff,is1,is2,
     .    ndhrs,isp,nsp,nlma,nlmb,hij(1,1,1,1+offi,1+offj),gv,
     .    kv,ng,tau,plat,pos,rtab,hrs,hrs,nhit,nmiss,dmx1)
      endif
      dmx(1) = max(dmx(1),dmx1(1))
      dmx(2) = max(dmx(2),dmx1(2))
      nviax = nviax + nhit
      nvtot = nvtot + nhit + nmiss
      deallocate(rtab,gv,kv)
      offj = offj + offH(1,1,ib+1) - offH(1,1,ib)
      enddo
      offi = offi + offH(1,1,ia+1) - offH(1,1,ia)
      enddo
      end

      subroutine hft2r2(mode,n1,n2,n3,ia,ib,idim,jdim,iax,qoff,is1,is2,
     .  ndhrs,isp,nsp,nlma,nlmb,hij,gv,kv,ng,tau,plat,pos,rtab,hrs,hrsr,
     .  nhit,nmiss,dmx)
C- For each connecting vector in a list taken from a uniform mesh of
C- translation vectors find corresponding place r.s neighbor table, and:
C  either mark iax table, or copy into hrs h generated from uniform mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :1s digit specifies what is created.
Ci           :  0 marks entries in row 8 of iax table; see iax below
Ci           :    In this mode, hij and hrs are not used
Ci           : >0 copy hij to hrs
Ci           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Ci           :  1 Assume hrs is real
Ci           :Other digits are not used.
Ci   n1,..n3:number divisions in QP mesh
Cr   ia    :site index for field point (or row or augmentation dim)
Ci   ia    :site index for source point (or col or basis dimension)
Ci   idim  :augmentation dimension of hij
Ci   jdim  :basis dimension of hij
Ci   is1,is2:Look in iax table for pairs is1..is2
Ci   ndhrs :leading dimensions of hrs: must be at least as large as
Ci         :the total number of orbitals of any atom.
Ci   isp   :current spin channel (1 or 2): poke hij into hrs(isp)
Ci         :NB: hij has no spin index
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlma  :number of orbitals in field  (row) dimension
Ci   nlmb  :number of orbitals in source (col) dimension
Ci   hij   :R.S hamiltonian on uniform mesh of translation vectors
Ci   gv    :list of r.s. lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f) mapping gv
Ci         :to the uniform mesh.  Not used if 1s digit mode is zero
Ci   ng    :number of G-vectors
Ci   tau   :offset pos(ib)-pos(ia) added to lattice vectors gv to make
Ci         :true connecting vectors.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Cio Inputs/Outputs
Cio   iax  :neighbor table containing pair information (pairc.f)
Cio        :Output if 1s digit mode is zero:
Cio        :iax(8,k) set to -1 for pairs k that correspond
Cio        :         to a connecting vector
Cio        :if 1s digit mode is >0, iax is not changed.
Co   hrs   :see hrsr
Co   hrsr  :copy hij to hrs or hrsr for pairs found in aix table.
Co         :hij copied EITHER into hrs (complex) OR hrsr (real)
Co         :depending on 10s digit mode
Co         :Also copy only takes place if 1s digit mode>0
Co   nhit  :number of connecting vectors found in aix table
Co   nmiss :number of connecting vectors not found in aix table
Co   dmx   :dmx(1) = largest real element found in h
Co         :dmx(2) = largest imaginary element found in h
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu   23 Jun 02  Adapted from pgfq2r
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,idim,jdim,ia,ib,n1,n2,n3,ndhrs,isp,nsp,nlma,nlmb,is1,
     .  is2,mode,nhit,nmiss
      parameter (niax=10)
      integer ng,kv(ng,3),iax(niax,*)
      double precision tau(3),gv(ng,3),rtab(3,*),pos(3,*),plat(3,3)
      double precision hrsr(ndhrs,ndhrs,nsp,is2),dmx(2),qoff(3)
      double complex hrs(ndhrs,ndhrs,nsp,is2),hij(n1,n2,n3,idim,jdim)
C ... Local parameters
      logical ltrans
      integer ig,i1,i2,i3,is,ix,offx,irow,icol,mode0,mode1,j,k
      double precision v(3),tol,dimx,dimy,TdotK,twopi
      double complex phase
C     character*20 strn
      parameter (tol=1d-5)
C ... External calls
      external isanrg

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      ltrans = mod(mode/10000,10) /= 0
      dimx = 0
      dimy = 0
      offx = is1-1
      twopi = 8*datan(1d0)

C --- Make -(connecting vector) = rtab = (pos(ia=field)-pos(ib=src)) ---
C     NB : Connecting vector is pos(src)-pos(field) = -rtab
      do  is = is1, is2
      do  ix = 1, 3
        rtab(ix,is-offx) = pos(ix,ia) - pos(ix,ib)
     .    + plat(ix,1)*iax(3,is)
     .    + plat(ix,2)*iax(4,is)
     .    + plat(ix,3)*iax(5,is)
      enddo
      enddo
C     call prmx('rtab',rtab,3,3,is2-is1+1)

C --- Loop through list of R vectors ---
      nhit = 0
      nmiss = 0
      do  ig = 1, ng
        if (mode0 /= 0) then
          i1 = kv(ig,1)
          i2 = kv(ig,2)
          i3 = kv(ig,3)
        endif
        v(1) = -(gv(ig,1) - tau(1))
        v(2) = -(gv(ig,2) - tau(2))
        v(3) = -(gv(ig,3) - tau(3))

C   ... If a match in the iax table is found, copy hij there
        do  is = is1, is2
          if (ia /= iax(2,is)) cycle
          if (abs(v(1)+rtab(1,is-offx)) < tol) then
          if (abs(v(2)+rtab(2,is-offx)) < tol) then
          if (abs(v(3)+rtab(3,is-offx)) < tol) then

            nhit = nhit+1

C            if (iax(1,is) == 1 .and. iax(2,is) == 1 .and. mode0 /= 0)
C     .        then
CC              print 335, 'ok',is,ig,i1,i2,i3,v,dble(hij(i1,i2,i3,1,1))
CC  335         format(1x,a,' connecting vector ib-ia ',2i4,2x,3i3,6f12.6)
C              print 336,
C     .          'ok',is,ig,i1,i2,i3,iax(3,is),iax(4,is),iax(5,is),
CC     .          mod(i1,6)-1-iax(3,is),
CC     .          mod(i2,6)-1-iax(4,is),
CC     .          mod(i3,6)-1-iax(5,is)
C     .          mod(iax(3,is)+6,6)+1-i1,
C     .          mod(iax(4,is)+6,6)+1-i2,
C     .          mod(iax(5,is)+6,6)+1-i3,
C     .          hij(i1,i2,i3,6,6)
C  336         format(1x,a,2i4,2x,3i3,2x,3i3,2x,3i3,6f12.6)
C            endif
C            if (is >= 65 .and. is <= 65+9 .and. iax(2,is) == 1 .and.
C     .        mode0 /= 0) then
C              print 335, 'ok',is,ig,i1,i2,i3,v,dble(hij(i1,i2,i3,1,1))
C            endif

            if (mode0 == 0) then
              iax(8,is) = -1

            else

              TdotK = 0
              do  j = 1, 3
                do  k = 1, 3
                  TdotK = TdotK - twopi*qoff(j)*plat(j,k)*iax(2+k,is)
                enddo
              enddo
              phase = cdexp(dcmplx(0d0,TdotK))

              if (mode1 == 0) then
                if (ltrans) then
                  do  icol = 1, nlmb
                    do  irow = 1, nlma
                      hrs(icol,irow,isp,is) = hij(i1,i2,i3,irow,icol)*phase
                      dimx = max(dimx,abs(dble(hij(i1,i2,i3,irow,icol))))
                      dimy = max(dimy,abs(dimag(hij(i1,i2,i3,irow,icol))))
                    enddo
                  enddo
                else
                  do  icol = 1, nlmb
                    do  irow = 1, nlma
                      hrs(irow,icol,isp,is) = hij(i1,i2,i3,irow,icol)*phase
                      dimx = max(dimx,abs(dble(hij(i1,i2,i3,irow,icol))))
                      dimy = max(dimy,abs(dimag(hij(i1,i2,i3,irow,icol))))
                    enddo
                  enddo
                endif
              elseif (ltrans) then
                do  icol = 1, nlmb
                  do  irow = 1, nlma
                    hrsr(icol,irow,isp,is) = hij(i1,i2,i3,irow,icol)*phase
                    dimx = max(dimx,abs(dble(hij(i1,i2,i3,irow,icol))))
                    dimy = max(dimy,abs(dimag(hij(i1,i2,i3,irow,icol))))
                  enddo
                enddo
              else
                do  icol = 1, nlmb
                  do  irow = 1, nlma
                    hrsr(irow,icol,isp,is) = hij(i1,i2,i3,irow,icol)*phase
                    dimx = max(dimx,abs(dble(hij(i1,i2,i3,irow,icol))))
                    dimy = max(dimy,abs(dimag(hij(i1,i2,i3,irow,icol))))
                  enddo
                enddo
              endif
            endif
            goto 10
          endif
          endif
          endif
        enddo
C       print 335, 'NO',is,ig,i1,i2,i3,v
        nmiss = nmiss+1

C        call cwrite(' quit (q) ',0,9,0)
C        read(*,'(a20)',err=99,end=99) strn
C        if (strn == 'x') stop
C        if (strn == 'q') goto 99
   10 enddo
C  99 continue

      dmx(1) = dimx
      dmx(2) = dimy

      end

      subroutine hft2r3(ib1,ib2,iwk,ntab,iax,nttab,mxcsiz)
C- Reduce iax table to those marked
C ----------------------------------------------------------------------
Ci Inputs
Ci   ib1   :compact iax starting at site ib1
Ci   ib2   :compact iax ending at site ib2
Ci   iwk   :integer work array of dimension max cluster size for
Ci         :any site.
Cio Inputs/Outputs
Cio  ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio        :Pairs for which iax(8,is) is set to -1 are retained;
Cio        :the remaining pairs are purged from the list.
Cio        :On output iax(8,*) is reset to zero.
Co   nttab :number of pairs in reduced table
Co   mxcsiz:maximum cluster size (number of pairs connected to any site)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   23 Jun 02  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,ib1,ib2,ntab(*),iwk(*),nttab,mxcsiz
      parameter (niax=10)
      integer iax(niax,*)
C ... Local parameters
      integer is1,is2,ib,nstart,i,nib
C     integer iprint,lgunit,stdo
C ... External calls
      external icopy,ivheap,ivprm

C --- Reorder iax table site-by-site ---
      do  ib = ib1, ib2

        is1 = ntab(ib)+1
        is2 = ntab(ib+1)
        nstart = is2-is1+1
C       Copy iax(8) to iwk(1+nstart)
        call icopy(nstart,iax(8,is1),niax,iwk(1+nstart),1)
C       Sort, using iwk(1) as permutation array
        call ivheap(1,nstart,iwk(1+nstart),iwk,101)
C       Rearrange iax, putting marked ones at top
        call ivprm(niax,nstart,iax(1,is1),iwk(1+nstart),iwk,1)
      enddo

C --- Condense iax table and ntab ---
      mxcsiz = 0
      nttab = ntab(ib1)
      do  ib = ib1, ib2
C       Count nib = number of marked pairs for this site
        is1 = ntab(ib)+1
        is2 = ntab(ib+1)
        nib = 0
        do  i = is1, is2
          nib = i-is1
          if (iax(8,i) == 0) goto 5
        enddo
        nib = is2-is1+1
    5   continue
        mxcsiz = max(mxcsiz,nib)
C       Append iax for this ib to the end of now shortened table
        if (is1 /= nttab+1) then
C         print *, 'copy iax for site',ib
          call icopy(niax*nib,iax(1,is1),1,iax(1,nttab+1),1)
        endif
        ntab(ib) = nttab
        nttab = nttab+nib
C       print 369, ib, nib, ntab(ib)
C 369   format(3i5)
      enddo
      ntab(ib2+1) = nttab
C --- Zero out iax(6..10) ---
      do  i = ntab(ib1)+1, nttab
        iax(6,i) = 0
        iax(7,i) = 0
        iax(8,i) = 0
        iax(9,i) = 0
        iax(10,i) = 0
      enddo
C     print 369, ib2+1, ntab(ib2+1)

C ... Printout
C      if (iprint() >= 30) then
C        stdo = lgunit(1)
C        call awrit4(' hft2rs (%i sites): neighbor table reduced from '//
C     .    '%i to %i pairs (%i/site)',' ',80,stdo,ib2-ib1+1,is2,
C     .    ntab(ib2+1)-ntab(ib1),(ntab(ib2+1)-ntab(ib1))/(ib2-ib1+1))
C        call info(70,0,0,' ntab: %n:1i',ib2-ib1+2,ntab)
C      endif

      end
