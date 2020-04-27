      subroutine optdme(opt,nfc,nrfn,nl,nbas,isp,nsp,nspc,nev,ipc,zgradm,
     .  nfilo,nfiup,nemlo,nemup,fac,ccp,optme)
C- Sphere contributions to optical matrix element <phi_j (fac*grad) phi_i>
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 Add to optme :      <j | (fac*grad) | i>^transpose
Ci         :1 Add to optme :  cc [<j | (fac*grad) | i>^transpose]
Ci   nfc   :Dimensions ccp,gradm : nfc must be at least max(nf1,nf2)
Ci   nrfn  :site-dependent number of radial functions
Ci   nl    :1 + l-cutoff for partial waves
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2) used only as index to gradm
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nev   :number of eigenvectors for this k
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci         :if zgradm is organized by site instead of class, :pass ipc(1) <= 0
Ci   zgradm:zgradm(ilm1,ilm2,1:3,if2,if1,isp,ic):
Ci         :<g2 grad g1> for: 1:3 corresponding to x,y,z  or (r+, r-, z) components.
Ci         :g1 = g1(r,if1,isp) Y_L(ilm1)  and  g2 = g2(r,if2,isp) Y_L(ilm2)
Co   nfilo,nfiup :Matrix elements for occupied bands restricted to nfilo:nfiup
Co   nemlo,nemup :Matrix elements for unoccupied bands restricted to nemlo:nemup
Ci   fac   :scale contribution to optme by fac
Ci   ccp   :Coefficients for projection of eigenfunction onto partial waves
Ci         :ccp(ilm,ib,jsp,ifn,j) =  Portion of eigenfunction j projected onto
Ci         :                         site ib, radial function ifn, L=ilm
Ci         :                         jsp ranges from 1..nspc.
Co Outputs
Co   optme :Augmented wave contribution to matrix element  <j | (fac*grad) | i>
Co         :connecting unocc states j to occ states i
Co         :Note optdme returns transpose or complex conjugate of <j | (fac*grad) | i>
Cd Debugging
Cd   Uncomment lines dumping ccp, grad1, grad2, grad3, copy to:
Cd     ccp(if1=1) -> zr1  ccp(if1=2) -> zr2   grad1 -> rcgx21  grad2 -> rcgy21  grad3 -> rcgz21
Cd   Also uncomment lines dumping 'z^+ zgradm z
Cd   Checks: use for test 6:9 for occ and 10:16 for unocc, write ME to file out.ogan
Cd   Compare z+ zgradm z to out.ogan:
Cd   mc -f9f15.10 zr2 -coll 10:16 -cc -t rcgx21 -x zr1 -coll 6:9 -x out.ogan -- -px
Cd   mch -f9f12.6 z -coll 6:9 -a zo z -coll 10:16 -a zu zu -cc -t grad1 -x zo -x -t  out.ogan --
Cr Remarks
Cr   Optics code adapted by Sergey Rashkeev with Walter Lambrecht from V. Antropov
Cr   Rewritten by MvS to use either real or spherical harmonics
Cr   ccp and zgradm must be in consistent representations (real or spherical harmonics)
Cr   zgradm may correspond to (x, y, z) or (r+, r-, z) components.
Cu Updates
Cu   02 May 18 Revised treatment of spherical harmonics
Cu   21 Aug 14 Revisited with simplified treatment of range in occ, unocc states
Cu   22 Jun 14 Rewritten for the noncollinear case: gradm is dimensioned differently
Cu   12 Jun 14 Symmetrization part split out of optdme to handle FP case.
Cu             See symdme for symmetrization part.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nev,isp,nfc,nl,nbas,nfilo,nfiup,nemlo,nemup,nsp,nspc
      integer ipc(nbas),nrfn(nbas)
      double precision fac
      double complex zgradm(nl**2,nl**2,3,nfc,nfc,nsp*nspc,nbas)
      double complex ccp(0:nl*nl-1,nbas,nspc,nfc,nev)
      complex(8) :: optme(nfilo:nfiup,nemlo:nemup,3) ! Output
C ... Local parameters
C     logical :: debug=.false.
      integer ibas,i,j,k,if1,if2,ndc,nl2,nocc,nunocc,ic
      integer ispc,jspc,ksp
      complex(8), allocatable :: zg(:,:),zgz(:,:,:)

C      if (debug) then
C      ibas = 1; if1 = 1; if2 = 1
C      call yprmi('ccp ib=%i if1=1',ibas,if1,3,ccp(0:nl*nl-1,ibas,1,1,1:nev),0,nl*nl,nl*nl,nev)
C      call yprmi('ccp ib=%i if1=2',ibas,if1,3,ccp(0:nl*nl-1,ibas,1,2,1:nev),0,nl*nl,nl*nl,nev)
C      call yprmi('grad1 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],3,zgradm(1,1,1,if2,if1,isp,ibas),0,nl*nl,nl*nl,nl*nl)
C      call yprmi('grad2 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],3,zgradm(1,1,2,if2,if1,isp,ibas),0,nl*nl,nl*nl,nl*nl)
C      call yprmi('grad3 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],3,zgradm(1,1,3,if2,if1,isp,ibas),0,nl*nl,nl*nl,nl*nl)
C      endif

C ... Accumulate sum_(sites,partial waves) < i | grad | j > for each (i,j,polarization)
C     Debugging: put evec into zc and grad into grad1
C     mch -f9f15.10 -vfac=-1 zc -coll 6:9 -a zo zc -coll 10:16 -a zu zu -cc -t grad1 -sfac -x zo -x -t out.zbgan -cc --
      nl2 = nl*nl; ndc = nl2*nbas*nspc*nfc
      nocc = nfiup-nfilo+1; nunocc = nemup-nemlo+1
      allocate(zg(nunocc,nl2),zgz(nemlo:nemup,nfilo:nfiup,3))
      call dpzero(zgz,2*size(zgz))

      do  k = 1, 3
        do  ibas = 1, nbas
          ic = ibas; if (ipc(1) > 0) ic = ipc(ibas)

          do  if1 = 1, nrfn(ibas)
          do  if2 = 1, nrfn(ibas)

C           if (ibas/=1) cycle
C           if (if1/=1 .or. if2/=1 .or. ibas/=1 .or. k/=1) cycle
C           if (if1/=1 .or. if2/=1) cycle
C           if (if1/=1) cycle

            do  ispc  = 1, nspc
              jspc = ispc
              ksp = max(isp,ispc)

C             z^+ * (fac * zgradm) for this (if2,if1,ibas,ksp) block
              call zgemm('C','N',nunocc,nl2,nl2,fac*(1d0,0d0),
     .          ccp(0,ibas,ispc,if2,nemlo),ndc,zgradm(1,1,k,if2,if1,ksp,ic),nl2,
     .          (0d0,0d0),zg,nunocc)
C              call zprm('fac * z^+ zgradm',2,zg,nunocc,nunocc,nl2)

C             Add contribution z^+ (fac * zgradm) z from this (if2,if1,ibas,ksp) block
              call zgemm('N','N',nunocc,nocc,nl2,(1d0,0d0),
     .          zg,nunocc,ccp(0,ibas,jspc,if1,nfilo),ndc,
     .          (1d0,0d0),zgz(nemlo,nfilo,k),nunocc)
C              call yprmi('z^+ zgradm z  ib k = %2:1i  if2 if1 = %2:1i',
C     .          [ibas,k],[if2,if1],3,zgz(nemlo,nfilo,k),0,nunocc,nunocc,nocc)
C              print *, '!!'; zgz = 0

          enddo                   ! loop over ispc

          enddo                   ! loop over over if2
          enddo                   ! loop over over if1
  991     continue
        enddo                   ! loop over ibas

C       Matrix element of grad operator
        if (mod(opt,10) == 0) then
          forall (i=nfilo:nfiup, j=nemlo:nemup) optme(i,j,k) = optme(i,j,k) + zgz(j,i,k)
C       Hermitian conjugate of matrix element
        else
          forall (i=nfilo:nfiup, j=nemlo:nemup) optme(i,j,k) = optme(i,j,k) + dconjg(zgz(j,i,k))
        endif

C       debugging
C        if (debug) then
C        call yprmi('optdme: optme k=%i isp=%i',k,isp,3,optme(nfilo,nemlo,k),0,nocc,nocc,nunocc)
C        endif

      enddo                       ! loop over over (xyz)
      deallocate(zg)

C     debugging
C      if (debug) call info5(1,0,0,' sumcheck optme %2;11,6D%2;11,6D %2;11,6D',
C     .  sum(optme(:,:,1)),sum(optme(:,:,2)),sum(optme(:,:,3)),4,5)

      end

      subroutine symdme(loptme,g,ngrp,nfiloe,nfiupe,nemloe,nemupe,
     .        nfilo,nfiup,nemlo,nemup,fac,optme,optmt,optmc)
C- Symmetrize optical matrix elements
C ----------------------------------------------------------------------
Ci Inputs
Ci   loptme: 0 do nothing
Ci         : 1 return optmt = symmetrized |optmt|^2
Ci         : 2 return optmc = symmetrized optmt
Ci         : 3 combination of 1+2
Ci         : 4 return optmc = optmt(k), not symmetrized
Ci         : 5 combination of 1+4
Ci   g,ngrp:group operations and number
Ci   nfiloe,nfiupe :Matrix elements for occupied bands in range nfilo:nfiup
Ci   nemloe,nemupe :Matrix elements for unoccupied bands in range nemlo:nemup
Ci   nfilo,nfiup :Dimensions optmt
Ci   nemlo,nemup :Dimensions optmt
Ci   fac   :scale optical matrix elements by fac
Ci   optme : <i | grad | j> connecting
Ci         :all occ states i, ocrng(1)<=i<=ocrng(2) to
Ci         :to unocc states j, unrng(1)<=j<=unrng(2)
Co Outputs
Co   optmt : optmt is returned if loptme is 1 or 3 or 5:
Co         : optmt(1:3,:,:) += fac * symmetrized |optme(:,:,1:3)|^2
Co   optmc : if loptme is 2 or 3 (probably not useful!)
Co         : optmc(1:3,:,:) += fac * symmetrized optme(:,:,1:3)
Co         : if loptme is 4 or 5:
Co         : optmc(1:3,:,:) += fac optme(:,:,1:3)
Cu Updates
Cu   21 Aug 14 Revisited with corrected form of matrix elements
Cu   12 Jun 14 Symmetrization part split out of optdme to handle FP case.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer loptme
      integer nfiloe,nfiupe,nemloe,nemupe,nfilo,nfiup,nemlo,nemup,ngrp
      double precision g(3,3,48),fac
      double complex optme(nfiloe:nfiupe,nemloe:nemupe,3)
C     Outputs
      double precision optmt(3,nfilo:nfiup,nemlo:nemup)
      double complex optmc(3,nfilo:nfiup,nemlo:nemup)
C ... Local parameters
      logical :: debug=.false.
      integer igrp,i,j,k
      double precision opt2(3)
      double complex opts(3),optc(3)

      if (nfiloe<nfiloe .or. nfiupe>nfiup .or.
     .    nemloe<nemlo  .or. nemupe>nemup)
     .  call rx('SYMDME: bad dimensions for optmt')

C     debugging
      if (debug) call info5(1,0,0,' sumcheck symdme %2;11,6D%2;11,6D %2;11,6D',
     .  sum(optme(:,:,1)),sum(optme(:,:,2)),sum(optme(:,:,3)),4,5)
C      do  k = 1, 3
C        call yprmi('symdme optme,  k = %i',
C     .    k,0,3,optme(nfilo,nemlo,k),0,nfiup-nfilo+1,nfiup-nfilo+1,nemup-nemlo+1)
C      enddo


C --- Loop over all (occ, unocc) pairs ---
      do  i = nfiloe, nfiupe
        do  j = nemloe, nemupe

        opt2 = 0; optc = 0
        do  igrp = 1, ngrp
          call zgrpop(optme(i,j,1:3),opts,g,igrp)
          opts = fac*opts
          do  k = 1, 3
            opt2(k) = opt2(k) + (dble(opts(k))**2+dimag(opts(k))**2)
            optc(k) = optc(k) + opts(k)
          enddo
        enddo
        if (mod(loptme,2) == 1) optmt(:,i,j) = optmt(:,i,j) + opt2(:)/ngrp
        if (loptme == 2 .or. loptme == 3) optmc(:,i,j) = optmc(:,i,j) + optc(:)/ngrp
        if (loptme == 4 .or. loptme == 5) optmc(:,i,j) = optmc(:,i,j) + fac*optme(i,j,1:3)

      enddo ! unoccupied states
      enddo ! occupied states

C     debugging
      if (debug) call info5(1,0,0,' sumcheck optme %;12,6D %;12,6D %;12,6D',
     .  sum(optmt(1,:,:)),sum(optmt(1,:,:)),sum(optmt(1,:,:)),4,5)
C      do  k = 1, 3
C        call yprmi('symdme optmt,  k = %i',
C     .    k,0,1,optmt(k:k,:,:),0,nfiup-nfilo+1,nfiup-nfilo+1,nemup-nemlo+1)
C      enddo

      end
