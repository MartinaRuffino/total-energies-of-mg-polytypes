      subroutine ooptdme(nfc,nrfn,nl,nlm,nbas,isp,nsp,nspc,nev,ipc,gradm,
     .  nfilo,nfiup,nemlo,nemup,fac,ccp,optme)
C- Sphere contributions to optical matrix element <phi_j (-grad) phi_i>^*
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfc   :Dimensions ccp,gradm : nfc must be at least max(nf1,nf2)
Ci   nrfn  :site-dependent number of radial functions
Ci   nl    :1 + l-cutoff for partial waves
Ci   nlm   :site-dependent cutoff (1+lmax)**2 for partial waves
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2) used only as index to gradm
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nev   :number of eigenvectors for this k
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   gradm :gradm(1:3,if1,if2,ilm1,ilm2,isp,ic):
Ci         :<g2 grad g1> for: 1:3 corresponding to x,y,z components and
Ci         :g1 = g1(r,if1,isp) Y_L(ilm1)  and  g2 = g2(r,if2,isp) Y_L(ilm2)
Co   nfilo,nfiup :Matrix elements for occupied bands in range nfilo:nfiup
Co   nemlo,nemup :Matrix elements for unoccupied bands in range nemlo:nemup
Ci   fac   :scale contribution to optme by fac
Ci   ccp   :Coefficients fixing projection of eigenfunction onto partial waves
Ci         :ccp(ilm,ib,jsp,ifn,j) =  Portion of eigenfunction j projected onto
Ci         :                         site ib, radial function ifn, L=ilm
Ci         :                         jsp ranges from 1..nspc.
Ci         :ccp should correspond to partial waves in spherical harmonics
Co Outputs
Co   optme :Augmented wave contribution to matrix element [<j | (fac*grad) | i>]^*
Co         :connecting unocc states j to occ states i
Co         :is added to optme with scale factor fac
Cr Remarks
Cr   Optics code adapted by Sergey Rashkeev with Walter Lambrecht from V. Antropov
Cu Updates
Cu   21 Aug 14 Revisited with simplified treatment of range in occ, unocc states
Cu   22 Jun 14 Rewritten for the noncollinear case: gradm is dimensioned differently
Cu   12 Jun 14 Symmetrization part split out of optdme to handle FP case.
Cu             See symdme for symmetrization part.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nev,isp,nfc,nl,nbas,nfilo,nfiup,nemlo,nemup,nsp,nspc
      integer ipc(nbas),nrfn(nbas),nlm(nbas)
      double precision fac
      double precision gradm(3,nfc,nfc,nl**2,nl**2,nsp*nspc,nbas)
      double complex ccp(0:nl*nl-1,nbas,nspc,nfc,nev)
C     Outputs
      complex(8) :: optme(nfilo:nfiup,nemlo:nemup,3)
C ... Local parameters
      integer i,j,ibas,l1,m1,k,l2,m2,li,lf,if1,if2
      integer ispc,jspc,kspc
      double precision opt1(3,2),optrr(3),optii(3),xx
      double complex opttt

C      ibas = 1; if1 = 1; if2 = 1
C      call yprmi('ccp before rotation ib=%i if1=1',ibas,if1,3,ccp(0:nl*nl-1,ibas,1,1,1:nev),0,nl*nl,nl*nl,nev)
C      call yprmi('grad1 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],1,gradm(1,if1,if2,1:nl*nl,1:nl*nl,isp,ibas),0,nl*nl,nl*nl,nl*nl)
C      call yprmi('grad2 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],1,gradm(2,if1,if2,1:nl*nl,1:nl*nl,isp,ibas),0,nl*nl,nl*nl,nl*nl)
C      call yprmi('grad3 ib = %i  if2 if1 = %2:1i',ibas,[if2,if1],1,gradm(3,if1,if2,1:nl*nl,1:nl*nl,isp,ibas),0,nl*nl,nl*nl,nl*nl)

C --- Loop over all (occ, unocc) pairs ---
      do  i = nfilo, nfiup
        do  j = nemlo, nemup

C   ... Make < i | grad | j > for each (i,j,polarization)
        do  k = 1, 3
          opttt = (0d0,0d0)
          do  l1 = 0, nl-1
          do  m1 = -l1, l1
          do  l2 = 0, nl-1

            if (l2 /= (l1+1) .and. l2 /= (l1-1) ) cycle
            if (k == 1) then
              m2 = m1+1
              if (m2 > l2) cycle
            endif
            if (k == 2) then
              m2 = m1-1
              if (m2 < -l2) cycle
            endif
            if (k == 3) then
              m2 = m1
              if (m2 < -l2 .or. m2 > l2) cycle
            endif

            li = l1*l1 + l1 + m1
            lf = l2*l2 + l2 + m2
C           if (li < 0 .or. lf < 0) cycle
C           if (lf < 0) call rx('optdme: should not happen')
            do  ibas = 1, nbas
            if (li > nlm(ibas) .or. lf > nlm(ibas)) cycle

            if (nspc == 1) then
              do  if2 = 1, nrfn(ibas)
              do  if1 = 1, nrfn(ibas)

C               if (if1 /= 1 .or. if2 /= 1) cycle
C               if (if1 /= 2 .or. if2 /= 1) cycle
C               if (if1 /= 1) cycle

                xx = gradm(k,if2,if1,lf+1,li+1,isp,ipc(ibas))
                if (k == 2) xx = -xx ! Old gradm uses opposite sign for component 2
                opttt = opttt +
     .          dconjg(ccp(li,ibas,1,if1,i))*ccp(lf,ibas,1,if2,j)*xx


                if (i == 1 .and. j == 1 .and. k == 1) then
                  call info8(1,0,0,' l2,m2=%2:1,2i  l1,m1=%2:1,2i ilm1,2=%2:1,2i %4;12,6D  %;11,6D  %2;12,6D',
     .              [l2,m2],[l1,m1],[li+1,lf+1],
     .              [ccp(lf,ibas,1,if2,j),dconjg(ccp(li,ibas,1,if1,i))],
     .              xx,
     .              opttt,7,8)
                endif

              enddo             ! if1
              enddo             ! if2

            else
              do  ispc = 1, nspc
C             do  jspc = 1, nspc  ! spin diagonal part for optics
              jspc = ispc
              kspc = ispc + 2*(jspc-1)
              do  if2 = 1, nrfn(ibas)
              do  if1 = 1, nrfn(ibas)

                xx = gradm(k,if2,if1,lf+1,li+1,kspc,ipc(ibas))
                if (k == 2) xx = -xx ! Old gradm uses opposite sign for component 2
                opttt = opttt +
     .          dconjg(ccp(li,ibas,ispc,if1,i))*ccp(lf,ibas,jspc,if2,j)*xx

              enddo             ! if1
              enddo             ! if2
C             enddo             ! Loop over noncollinear spins
              enddo             ! Loop over noncollinear spins
            endif
          enddo                 ! Loop over sites

          enddo ! l2
          enddo ! m1
          enddo ! l1

          opt1(k,1) = dble(opttt)
          opt1(k,2) = dimag(opttt)
        enddo  ! Polarization loop

C   ... Real, imaginary parts of grad in Cartesian coordinates
C       optx = (-opt1-opt2)/sqrt(2)  opty = (-i*opt1 + i*opt2)/sqrt(2)  opt3 = -optz
C                  (-1  -i)                  (-1 -1)
C       ue = 1/sr2*(      )   inverse= 1/sr2*(     )
C                  (-1   i)                  (i  -i)
C Debugging:
C mc -vfac=1 -vsr2='sqrt(2)' zgz1 -s-1/sr2 zgz2 -s-1/sr2 -+ -sfac zgzx --
C mc -vfac=1 -vsr2='sqrt(2)' zgz1 -s0,-1/sr2 zgz2 -s0,1/sr2 -+ -sfac zgzy --
C       This is (-grad)
        optrr(1) = -(opt1(1,1)+opt1(2,1))/dsqrt(2d0)
        optrr(2) = +(opt1(1,2)-opt1(2,2))/dsqrt(2d0)
        optrr(3) = -opt1(3,1)

        optii(1) = -(opt1(1,2)+opt1(2,2))/dsqrt(2d0)
        optii(2) = -(opt1(1,1)-opt1(2,1))/dsqrt(2d0)
        optii(3) = -opt1(3,2)

C   ... Add -fac*<-grad> into optme
        do  k = 1, 3
          optme(i,j,k) = optme(i,j,k) - fac*dcmplx(optrr(k),optii(k))
        enddo

C   ... Write as spherical vectors for debugging
C        print *, '!!'
C        do  k = 1, 3
C          optme(i,j,k) = optme(i,j,k) + fac*dcmplx(opt1(k,1),opt1(k,2))
C        enddo

      enddo ! unoccupied states
      enddo ! occupied states

C     Debugging printout
C      do  k = 1, 3
C      call yprmi('ooptdme : optme k = %i',
C     .  k,0,3,optme(nfilo,nemlo,k),0,nfiup-nfilo+1,nfiup-nfilo+1,nemup-nemlo+1)
C      enddo
C
C      call info5(1,0,0,' sumcheck optdme %2;11,6D%2;11,6D %2;11,6D',
C     .  sum(optme(:,:,1)),sum(optme(:,:,2)),sum(optme(:,:,3)),4,5)
C

      end
