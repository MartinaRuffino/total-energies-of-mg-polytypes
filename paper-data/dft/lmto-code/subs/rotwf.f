      subroutine rotwf(opt,nl,nbas,nspc,bas,offH,indxsh,istab,g,ag,q,rmat,
     .  ldmpa,ldima,ldz,nev,z0,z)
C- Rotate a wave function z0 into z
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :1     Rotate z using g^-1
Ci         :10s digit
Ci         :0     z is in real spherical harmonics
Ci         :Add 2 if z is in true (complex) spherical harmonics
Ci         :      See rcub2sph.f for definition of spherical harmonics
Ci         :Add 4 to choose other sign of phase convention, i.e.
Ci         :      phase = q * [(g R_j + a) - R_i]
Ci         :100s digit distinguishes how complex arithmetic is handled
Ci         : 0: z,z0 have real, imaginary separated
Ci         :    z = z(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci         : 1: z,z0 are in complex*16 format
Ci         :    z = z(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci         :1000s digit
Ci         : 1: return in z rotation matrix R instead of R z
Ci   nl    :(global maximum l) + 1
Ci   nbas  :number of atoms in the basis (input)
Ci   bas   :basis vectors (input)
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci         :Site istab(i,ig) is transformed into site i by grp op ig.
Ci         :Negative sign of istab(1) => flags that group op comes from time reversal
Ci   g,ag  :space group operation; see Remarks
Ci   q     :qp for which z is sought.
Ci   rmat  :a work array of dimension at least nl**4
Ci         :If wave functions use true spherical harmonics, rmat
Ci         :must be dimensioned at least nl**4*2
Ci   ldmpa :offset to current downfolding block
Ci   ldima :last orbital in current downfolding block
Ci   ldz   :dimensions z
Ci   nev   :number eigenvectors to rotate
Ci   z0    :unrotated wave function; see Remarks
Co Outputs
Co   rmat  :rotation matrix for this g
Co   z     :is returned at q
Cb Bugs
Cb   The 2nd dimension of z0,z must be the same as ldz.  Matters if kcplx=0 or 2.
Cr Remarks
Cr   Input vectors z0 are rotated to z = R z0.
Cr   Rotating the crystal axes r -> g*r, so that Y_L(rhat) -> Y_L(g*rhat)
Cr   An eigenfunction is chi(\r) = sum_L\a z_l(|\r-\a|) YL((|\r-\a|) \phi_l(|\r-\a|)
Cr   Let R = rotation matrix mapping harmonics Y_L(g*rhat) = R Y_L(rhat)
Cr
Cr   Program computes z = R z0
Cr
Cr   Analogous to operation that maps starting qp to q, i.e.
Cr     z(q) = R z0 (g^-1 q).
Cr   In one case the crystal is rotated; in the other q is rotated.
Cr
Cr   There is also a phase shift because position a'_j = g a_j + ag
Cr   can differ from the corresponding
Cr   The Bloch sum picks up a phase shift
Cr     exp(i q (g a_j + a) - a_i)) where  (g a_j + ag) - a_i
Cu Updates
Cu   14 Jan 18 Extended to the noncollinear case
Cu   10 Aug 15 Allows additional rotations from time-reversal symmetry.
Cu   04 Jan 04 Allow rotation of w.f. represented in true
Cu             spherical harmonics
Cu   25 Apr 02 first cut at applicability to multiple-kappa basis.
Cu             bug: fp code has different phase convention.
Cu   02 Oct 01 Replace by call to prothl in roth, so downfolding works
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: nkap0=4,n0H=5
      integer opt,ldmpa,ldima,ldz,nev,nl,nbas,nspc,istab(nbas),offH(n0H,nkap0,*),indxsh(*)
      double precision bas(3,nbas),g(3,3),ag(3),q(3),rmat(nl*nl,nl*nl,2)
C     z,z0 dimensioned as if kcplx=2.  rotwf internally converts to kcplx=2 to perform rotation.
      double precision z(ldz,nspc,2,ldz),z0(ldz,nspc,2,ldz)
C ... Local parameters
      logical lroti
      integer ofzi,kcplx,ixx,ldzr,opt0,opt1,opt3,nl2,ldzx,isp,morder
      double precision wk(nl*nl,nl*nl,2),gt(3,3),xx
C     procedure(real(8)) :: det33

      call tcn('rotwf')

C     lunitg is true if unit rotation matrix
C      lunitg = ddot(9,g,1,g,1)-3 < 1d-10 .and.
C     .             abs(g(1,1)-1) < 1d-10 .and.
C     .             abs(g(2,2)-1) < 1d-10 .and.
C     .             abs(g(3,3)-1) < 1d-10

C     opt1 = mod(mod(opt/10,10),4)
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt3 = mod(opt/1000,10)
      nl2 = nl*nl
      ldzx = ldz*nspc
      call lmorder(0,morder,[0],[0])

C --- Setup for pseudorotation (rotation + inversion) ---
      lroti = .false.
      if (istab(1) < 0) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
C       TR symmetry : spin up and spin down must be flipped
        if (nspc == 2) call rx('rotwf needs still missing TR symmetry, nc case')
        lroti = .true.
      endif

C     call zprm('z0',2,z0,ldzx,ldzx,nev)
      kcplx = mod(opt/100,10)
      if (kcplx == 1) then
        kcplx = 2
        call ztoy(z0,ldzx,ldzx,nev,0)
      endif
      call cplxdm(kcplx,ldzx,ldzx,ixx,ixx,ldzr,ofzi)

C     call yprm('z0',kcplx+2,z0,ofzi,ldzx,ldzx,nev)

      gt = g; if (opt0 /= 0) call dinv33(g,0,gt,xx)
      if (mod(opt1,4) == 0) then
        call ylmrtg(nl2,gt,rmat)
      else
        call ylmrtg(nl2,gt,wk)
        call s2sph(kcplx+100*morder,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C       call yprm('rmat',kcplx+2,rmat,nl2*nl2,nl2,nl2,nl2)
      endif

      do  isp = 1, nspc
        call prothl(10*opt1+100*opt3,nl*nl,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,
     .    q,rmat,1,nev,1,nbas,0,ldzr,-1,0,z(1,isp,1,1),ofzi,ldzr,ixx,z0(1,isp,1,1),ofzi)
      enddo

C     call yprm('z(r)',kcplx+2,z,ofzi,ldzx,ldzx,nev)

      if (mod(opt/100,10) == 1) then
        call ztoy(z0,ldzx,ldzx,nev,1)
        call ztoy(z,ldzx,ldzx,nev,1)
      endif

C --- Cleanup for pseudorotation (rotation + inversion) ---
      if (lroti) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
        call protwft(mod(opt/100,10),z,z,z,ldzx,ldzx,ldima,nev)
      endif

      call tcx('rotwf')
C#endif


C --- OLD ... does not work with downfolding ---
C#ifdefC OLD
C      call tcn('rotwf')
C      twopi = 8d0*datan(1d0)
C      nl2 = nl*nl
C      call ylmrtg(nl2,g,rmat)
C      call dpzero(z,ldz*ldz*2)
C
C      do  10  ibas = 1, nbas
C        jbas = istab(ibas)
CC   ... tbas = (g R_j + a) - R_i; tbas should be a lattice vector
C        do  14  i = 1, 3
C        tbas(i) = ag(i) - bas(i,ibas)
C        do  16  j = 1, 3
C   16   tbas(i) = tbas(i) + g(i,j)*bas(j,jbas)
C   14   continue
CC   ... exp(i q tbas)
C        sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
C        cosP = dcos(sp)
C        sinP = dsin(sp)
C
CC   ... Offset to first orbital in ibas,jbas
C        offi = offH(1,1,ibas)
C        offj = offH(1,1,jbas)
C        nlm  = offH(1,1,ibas+1) - offH(1,1,ibas)
C
CC   ... For each l, do rmat z0(jbas,lm)
C        offr = 1
C        do  20  li = 0, nl-1
C          offr = li**2 + 1
C          nlmi = 2*li+1
C          if (indxsh(offj+offr) > ldz) goto 20
C
C          if (nlmi > ncut) then
C            call dgemm('N','N',nlmi,ldz,nlmi,1d0,rmat(offr,offr),nl2,
C     .        z0(offj+offr,1,1),ldz,0d0,z(offi+offr,1,1),ldz)
C            call dgemm('N','N',nlmi,ldz,nlmi,1d0,rmat(offr,offr),nl2,
C     .        z0(offj+offr,1,2),ldz,0d0,z(offi+offr,1,2),ldz)
C          else
C            offr = offr - 1
C            ir = offi+offr
C            jr = offj+offr
C            do  22  k = 1, ldz
C            do  22  j = 1, nlmi
C            do  22  i = 1, nlmi
C            z(i+ir,k,1) = z(i+ir,k,1) + rmat(i+offr,j+offr)*z0(j+jr,k,1)
C            z(i+ir,k,2) = z(i+ir,k,2) + rmat(i+offr,j+offr)*z0(j+jr,k,2)
C   22       continue
C          endif
C
C   20   continue
C
CC        if (abs(sp) > 1d-8) print 333, 'tbas=',ibas,tbas,istab,sp
CC  333   format(a,i3,3f14.8,3i3,f14.8)
C
C
CC   ... Multiply by phase
C        if (dabs(sp) > 1d-8) then
C          do  30  i = 1, nlm
C   30     call yscal(ldz,cosP,sinP,z(offi+i,1,1),z(offi+i,1,2),ldz)
C        endif
C   10 continue
C#endif

C#ifdefC DEBUG
C      ofzi = ldz*ldz
CC     zz = 0
C      call prothl(opt1,nl2,nbas,bas,0,ldz,indxsh,istab,g,ag,q,rmat,1,
C     .  ldz,1,nbas,0,ldz,-1,0,zz,ofzi,ldz,ldz,z0,ofzi)
C
C      do  i = 1, ldz
C      do  j = 1, ldz
C        if (cdabs2(z(i,j,1)-zz(i,j,1),z(i,j,2)-zz(i,j,2)) >= 1d-12) then
C          print *, i,j
C          call rx('oops')
C        endif
C      enddo
C      enddo
C      call tcx('rotwf')
C#endif
      end

      subroutine protwft(kcplx,z0,z1,z2,ldz,ldz2,ldim,nev)
C- Take the complex conjugate of a matrix z
      implicit none
C ... Passed parameters
      integer kcplx,ldz,ldz2,ldim,nev
      double precision z0(ldz,ldz2,2),z1(2,ldz,ldz2),z2(ldz,2,ldz2)
C ... Local parameters
      integer ic

C     call zprm('z before cc',2,z1,ldz,ldim,nev)

      do  ic = 1, nev

      if (kcplx == 0) then
        call dscal(ldim,-1d0,z0(1,ic,2),1)
      elseif (kcplx == 2) then
        call dscal(ldim,-1d0,z2(1,2,ic),1)
      else
        call dscal(ldim,-1d0,z1(2,1,ic),2)
      endif

      enddo

C     call zprm('z after cc',2,z1,ldz,ldim,nev)

      end
