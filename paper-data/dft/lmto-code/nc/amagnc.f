      subroutine amagnc(nbas,nl,ipc,rhos,nrhos,qnu,eula,neul,mode,
     .  amag,aamom,bxc)
C- Printout magnetic moments in unit cell
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0)
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :amagnc uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   mode  :0 use spin density matrix to make moments along z
Ci         :1 use qnus to make moments along qnu
Ci         :Add 2 to generate average magnetization direction bxc
Co Outputs
Co   amag(1..3): net system magnetic moment
Co   aamom :local magnetic moments
Co   bxc   :average magnetization direction for each class
C          :(made when mode>=2)
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma:
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr     M_z =  (rho11)-(rho22)
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be properly hermitian, e.g. when rhos is generated
Cr   by a Green's function technique.
Cb Bugs
Cb   This routine does not symmetrize over class.
Cb   The first member of the class is assumed to be representative.
Cu Updates
Cu   09 Apr 11 can generate bxc; also small bug fixes
Cu   21 Apr 04 Revised to properly accomodate m-dependent Euler angles
Cu   17 Feb 03 Revised amagnc; cleaner and bug fixes.
C ----------------------------------------------------------------------
      implicit none
C Passed variables
      integer nbas,neul,nl,nrhos,ipc(nbas),mode
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,*),
     .  qnu(3,nl,2,*),amag(3),aamom(nbas),bxc(3,*)
C Local variables
      logical lrhol
      integer ib,ic,lp1,lgunit,ipr,k,ilm,l,m,mode0,iclbsj
      integer stdo,nglob
      double precision alphan,betan,gamman,arhol,arhom,sarho,bxcl(3),
     .  rotg(3,3),amom,amlm(3),amgm(3),saml(3),samg(3),amgl(3),aml(3),
     .  fac
      integer PRT1,PRT2
      parameter (PRT1=30,PRT2=35)

C     character*1 cmode(2)
C     data cmode /'o','i'/


C      print *, '!!'; call pshpr(50)

C --- Setup ---
      call getpr(ipr)
      stdo = nglob('stdo')
      call sanrg(.true.,mode,0,3,' amagnc','mode')
      mode0 = mod(mode,2)
      amag(1) = 0
      amag(2) = 0
      amag(3) = 0
      sarho = 0
      lrhol = nrhos == nl .or. mode0 == 1
      if (nrhos /= nl .and. nrhos /= nl**2 .and. mode0 == 0)
     .  call rxi('amagnc: bad nrhos',nrhos)

      if (ipr >= 20) then
        if (mode0 == 0)write(stdo,332) 'density matrix'
        if (mode0 == 1)write(stdo,332) 'Euler angles and sphere charges'
        if (lrhol .and. ipr >= PRT1) write(stdo,335) '   l'
        if (.not. lrhol .and. ipr >= PRT1) write(stdo,335) ' ilm'
      endif

  332 format(/' AMAGNC: magnetic moments from ',a,':')
  335 format(16x,'global coordinates',15x,'local coordinates'/
     .        3x,'ib',a,'    Mx        My        Mz',
     .       10x,           'Mx        My        Mz        |M|')

C --- Magnetization vector for each site, orbital ---
      do  ib = 1, nbas
        call dpzero(bxcl,3)
        call dpzero(saml,3)
        call dpzero(samg,3)
        ic = ipc(ib)
        aamom(ib) = 0
        ilm = 0
        do  l = 0, nl-1
        lp1 = l+1
        call dpzero(amgl,3)
        call dpzero(aml,3)
        arhol = 0
        do  m = -l, l
          ilm = ilm+1
C         Rotation matrix for this ilm
          k = 1
          if (neul == nl) k = lp1
          if (neul == nl**2) k = ilm
          alphan = eula(ib,k,1)
          betan  = eula(ib,k,2)
          gamman = eula(ib,k,3)
          call eua2rm(alphan,betan,gamman,rotg)
C         print 336, 'rotg',((rotg(i,k),k=1,3),i=1,3)
C 336     format(' rot matrix:',a,/(3f15.10))

C         Magnetization for this ilm in local coordinates
          if (mode0 == 0) then
            if (.not. lrhol) then
              k = ilm
              fac = 1
            else
              fac = 1/dble(2*l+1)
              k = l+1
            endif
C           This is M in local coordinates, symmetrized form
            amlm(1) = (rhos(1,0,k,1,2,ic) + rhos(1,0,k,2,1,ic))*fac
            amlm(2) = (rhos(2,0,k,2,1,ic) - rhos(2,0,k,1,2,ic))*fac
            amlm(3) = (rhos(1,0,k,1,1,ic) - rhos(1,0,k,2,2,ic))*fac
          else
            amlm(1) = 0
            amlm(2) = 0
            amlm(3) = qnu(1,lp1,1,ic) - qnu(1,lp1,2,ic)
            call dscal(3,1/dble(2*l+1),amlm,1)
          endif
          arhom = dsqrt(amlm(1)**2 + amlm(2)**2 + amlm(3)**2)

C         Rotate amlm to global coordinates (amgm)
          call dgemm('T','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)
C         call dgemm('N','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)

C         l-resolved moment, global coordinates
          call daxpy(3,1d0,amgm,1,amgl,1)
C         ditto, local coordinates; meaningful only for lrhol=T
          call daxpy(3,1d0,amlm,1,aml,1)
          arhol = arhol + arhom
C         aamom(ib) = aamom(ib) + arhom

C          if (.not. lrhol .and. ipr >= PRT2) write(stdo,333) ib,ilm,amgm,amlm,arhom
C  333     format(i5,i4,3f10.6,2x,3f10.6,f10.6)
          if (.not. lrhol .and. ipr >= PRT2) then
            call info5(2,0,0,'%,5i%,4i%3;10,6D  %3;10,6D%;10,6D',ib,ilm,amgm,amlm,arhom)
          endif
C         End of loop over m
          enddo

          call daxpy(3,1d0,aml,1,saml,1)
          call daxpy(3,1d0,amgl,1,samg,1)
C          if (lrhol .and. ipr >= PRT2) write(stdo,333) ib,l,amgl,aml,arhol
          if (lrhol .and. ipr >= PRT2)
     .      call info5(2,0,0,'%,5i%,4i%3;10,6D  %3;10,6D%;10,6D',ib,l,amgl,aml,arhol)
        enddo
C       aamom(ib) = dabs(aamom(ib))
        amom = dsqrt(samg(1)**2 + samg(2)**2 + samg(3)**2)
        aamom(ib) = amom
        sarho = sarho + dabs(amom)
C        if (ipr >= PRT1) write(stdo,334) ib,samg,saml,amom
C  334   format(i5,'   *',3f10.6,2x,3f10.6,f10.6)
        if (ipr >= PRT1)
     .    call info5(2,0,0,'%,5i   *%3;10,6D  %3;10,6D%;10,6D',ib,samg,saml,amom,0d0)
        if (ipr >= PRT2) write(stdo,'(1x)')
        amag(1) = amag(1) + samg(1)/nbas
        amag(2) = amag(2) + samg(2)/nbas
        amag(3) = amag(3) + samg(3)/nbas
        if (mode >= 2) call daxpy(3,1d0,samg,1,bxcl,1)

C       If ib is first member of class, generate bxc(:,ic)
        if (mode >= 2 .and. ib == iclbsj(ic,ipc,nbas,1)) then
C         Normalize bxc, if it is nonzero
          amom = dsqrt(bxcl(1)**2 + bxcl(2)**2 + bxcl(3)**2)
          if (amom /= 0) call dscal(3,1/amom,bxcl(1),1)
C         If no orbital-dependence => use Euler angles
          if (neul == 1 .and. mode0 == 1 .or. amom == 0) then
            bxcl(1) = dcos(alphan)*dsin(betan)
            bxcl(2) = dsin(alphan)*dsin(betan)
            bxcl(3) = dcos(betan)
          endif
          call dcopy(3,bxcl,1,bxc(1,ic),1)
        endif

      enddo
      amom = dsqrt(amag(1)**2 + amag(2)**2 + amag(3)**2)
      call info5(20,0,0,' <Mx>=%,6;6d  <My>=%,6;6d  <Mz>=%,6;6d'//
     .  '  <M>=%,6;6d  <|M|>=%,6;6d',amag(1),amag(2),amag(3),amom,
     .  sarho/nbas)
      if (lgunit(2) > 0 .and. ipr > 1) then
        call awrit5(' <Mx>=%,6;6d  <My>=%,6;6d  <Mz>=%,6;6d'//
     .    '  <M>=%,6;6d  <|M|>=%,6;6d',
     .    ' ',80,lgunit(2),amag(1),amag(2),amag(3),amom,sarho/nbas)
      endif

      end
      subroutine amagn2(nl,nclass,nbas,ipc,eula,neul,bxc,qnu,
     .  nrhos,rhos)
C- Renormalize magnetic part of qnu, corresponding to average direction
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   nrhos :number of channels which spin density-matrix is stored
Ci   bxc   :average direction for XC field (amagnc.f)
Ci   rhos  :spin density-matrix.  rhos++ and rhos-- are modified.
Co Outputs
Co   qnu   :energy-weighted moments of the sphere charges
Cl Local variables
Cl         :
Cr Remarks
Cr   If bxc is scaled by -1, the magnetic part of qnu is scaled by -1.
Cr   Only the product is relevant.  This routine scales bxc by -1 if
Cr   the moment from the given qnu is negative.
Cb Bugs
Cb   This routine does not symmetrize over class.
Cb   The first member of the class is assumed to be representative.
Cu Updates
Cu   09 Apr 11 Completely rewritten
Cu   09 Apr 04 First created.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nclass,nbas,neul,nrhos,ipc(nbas)
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,nclass)
      double precision qnu(3,nl,2,*),bxc(3,*)
C ... Local parameters
      logical lrhol
      integer ilm,l,m,k,i,ic,ib,iclbsj,lp1
C      integer nglob,stdo,ipr,PRT1,PRT2
C      parameter (PRT1=40,PRT2=50)
      double precision ddot,rotg(3,3),eulat(neul,3),
     .  alphan,betan,gamman,qp,qm,dsum,bhat(3),qli(0:2),fac
C      double precision pi
      double precision amgm(3),amlm(3),amdotb(0:2)

C     Even if input Euler angles all identical, output need not be
C     We will ignore this in the spirit of the ASA.
      if (neul == 1) return

      if (nrhos < neul .or. nrhos /= nl .and. nrhos /= nl*nl) call
     .  rx('amagn2: illegal dimension for dens matrix')
C      call getpr(ipr)
C      stdo = nglob('stdo')
C      pi = 4*datan(1d0)
      lrhol = nrhos == nl

C ... Determine whether to scale bxc by -1
      do  ic = 1, nclass
        fac = dsum(nl,qnu(1,1,1,ic),3) - dsum(nl,qnu(1,1,2,ic),3)
        if (fac < 0) call dscal(3,-1d0,bxc(1,ic),1)
      enddo

      do  ic = 1, nclass
      ib = iclbsj(ic,ipc,-nbas,1)
      if (ib > 0) then
        call dmscop(eulat,1,eula,nbas,ib,ib,1,neul*3,1,1,1d0)
        call dcopy(3,bxc(1,ic),1,bhat,1)

C --- Projection of rhos onto bxc; store to magnetic part of qnu ---
        ilm = 0
        do  l = 0, nl-1
        lp1 = l+1
        call dpzero(amdotb,3)
        call dpzero(qli,3)
        do  m = -l, l
          ilm = ilm+1
          k = 1
          if (neul == nl) k = lp1
          if (neul == nl**2) k = ilm
C         Rotation matrix for this ilm
          alphan = eula(ib,k,1)
          betan  = eula(ib,k,2)
          gamman = eula(ib,k,3)
          call eua2rm(alphan,betan,gamman,rotg)
C         print 336, 'rotg',((rotg(i,k),k=1,3),i=1,3)
C 336     format(' rot matrix:',a,/(3f15.10))

C         Magnetization for this ilm in local coordinates
          if (.not. lrhol) then
            k = ilm
            fac = 1
          else
            fac = 1/dble(2*l+1)
            k = l+1
          endif
C         This is M in local coordinates, symmetrized form
          do  i = 0, 2
            amlm(1) = (rhos(1,i,k,1,2,ic) + rhos(1,i,k,2,1,ic))*fac
            amlm(2) = (rhos(2,i,k,2,1,ic) - rhos(2,i,k,1,2,ic))*fac
            amlm(3) = (rhos(1,i,k,1,1,ic) - rhos(1,i,k,2,2,ic))*fac
C           Rotate amlm to global coordinates (amgm)
            call dgemm('T','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)
C           call dgemm('N','N',3,1,3,1d0,rotg,3,amlm,3,0d0,amgm,3)
C           Projection onto bhat
            qm = ddot(3,amgm,1,bhat,1)
C           Accumulate into l-resolved moment
            amdotb(i) = amdotb(i) + qm
C           m-averaged rho+ + rho-, rho+ - rho- in Bxc frame
            qli(i) = qli(i) +(rhos(1,i,k,1,1,ic)+rhos(1,i,k,2,2,ic))*fac
C         End of loop over i
          enddo
C       End of loop over m
        enddo

C   ... Accumulate qnu
        do  i = 0, 2
          qp = (qli(i)+amdotb(i))/2
          qm = (qli(i)-amdotb(i))/2
C         print *, sngl(qp),sngl(qm)
C         These should be equivalent to the above
          qp = (qnu(i+1,l+1,1,ic)+qnu(i+1,l+1,2,ic)+amdotb(i))/2
          qm = (qnu(i+1,l+1,1,ic)+qnu(i+1,l+1,2,ic)-amdotb(i))/2
C         print *, sngl(qp),sngl(qm)

          qnu(i+1,l+1,1,ic) = qp
          qnu(i+1,l+1,2,ic) = qm
        enddo

C       End of loop over l
        enddo

C     End of loop over class
      endif
      enddo

      end
C     this is old amagn2
C      subroutine amagn3(mode,nl,nlo,nclass,nbas,ipc,eula,neul,bxc,qnu,
C     .  nrhos,rhos)
CC- When rhos++ and rhos-- is only l-dependent, approximately adjust
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 use spin density matrix to make moments m3 along zhat
CCi         :1 use qnus to make moments along qnu
CCi   nl    :(global maximum l) + 1
CCi   nlo   :number of channels for which to accumulate ASA moments.
CCi         :nl0=nl when moments are contracted by m,
CCi         :nl0=nl**2 when moments are not contracted
CCi   nbas  :size of basis
CCi   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
CCi   eula  :Euler angles for noncollinear spins
CCi   neul  :1 if Euler angles are l-independent, nl otherwise
CCi   qnu   :energy-weighted moments of the sphere charges
CCi   nrhos :number of channels which spin density-matrix is stored
CCi   bxc   :average direction for XC field
CCio Inputs/Outputs
CCi   rhos  :spin density-matrix.  rhos++ and rhos-- are modified.
CCl Local variables
CCl         :
CCr Remarks
CCr   If because rhos++ and rhos-- were averaged over m (which can happen
CCr   if makwts accumulates doswts by l only; also for the ASA moments
CCr   they are always averaged over m), no account is
CCr   taken of their m-dependent orientation.  A better average to
CCr   rhos++ and rhos-- would be to make
CCr     Mbar = 1/(2l+1) sum_m Rm^-1 Mm(loc)
CCr   where Mm(loc) = m-dependent local magnetization in its local
CCr   coordinate system and Rm is the rotation matrix to loc. coords.
CCr   This is what would obtain if the entire rhos were averaged over m.
CCr   Finally, from Mbar we obtain the z-part of Mm in its local
CCr   coordinate system as Mmz(loc) = Rm Mbar.  We leave the x and y
CCr   parts untouched as they were explicitly resolved.
CCr   This is of course an approximation.  Better to allow the entire
CCr   rhos to be accumulated by m, which is accomplished by resolving
CCr   doswts by m (nlo = nl**2).
CCr
CCr   Correcting the ASA moments qnu+-qnu-.  A similar principle applies.
CCr   qnu+-qnu- has less information than the full rhos, because only
CCr   the z component is retained.  We construct qnu+-qnu- by requiring
CCr   that the rotation of qnu+-qnu- to the B-field axis matches
CCr   the rhos++ - rhos-- in this same reference.
CCr   If qnu+ - qnu- = (rho++ - rho--) in the B-field reference frame,
CCr   we assign qnu+ - qnu-  = (rho++ - rho--)/R_33 in the coordinates
CCr   of the orbital, where R is the rotation from the orbital frame
CCr   to the Bxc frame: R=Rb*Rm+ (see below).
CCr   When the orbital axis (spec by rm) and and Bxc align then
CCr   (Rb*Rm+)_33=1 and the correction vanishes.  This should always
CCr   occur in the rigid-spin approximation, when the each orbital
CCr   is collinear with Bxc.
CCu Updates
CCu   09 Apr 04  First created.
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,nl,nlo,nclass,nbas,neul,nrhos,ipc(nbas)
C      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,nclass)
C      double precision qnu(3,nl,2,*),bxc(3,*)
CC ... Local parameters
C      integer ilm,l,j,m,k,i,ic,ib,iclbsj
C      double precision mbar(3,0:2),amg(3,0:2),ammi(3),ddot,rotm(3,3),
C     .  alpha,beta,gamma,ql,amli(0:2),qp,qm,eulat(neul,3),ph,th,
C     .  hatbxc(3),pi,rotb(3,3),rotmb(3,3),rotmb3,rhom(3),qli(0:2),fac
C
C      pi = 4*datan(1d0)
C      do  ic = 1, nclass
C      ib = iclbsj(ic,ipc,-nbas,1)
C      if (ib > 0) then
C      call dmscop(eulat,1,eula,nbas,ib,ib,1,neul*3,1,1,1d0)
C      call dcopy(3,bxc(1,ic),1,hatbxc,1)
C      ilm = 0
C      do  l = 0, nl-1
C
CC       Retain starting ilm for this l
C        j = ilm
C
CC   --- Adjust rhos when ++,-- part averaged over m ---
C        if (nrhos > nlo) then
C
CC   ... Determine average M for this l
C        call dpzero(mbar,3*3)
C        do  m = -l, l
C          ilm = ilm+1
CC         Euler angles for this channel
C          k = 1
C          if (neul == nl) k = l+1
C          if (neul == nl*nl) k = ilm
C          alpha = eulat(k,1)
C          beta  = eulat(k,2)
C          gamma = eulat(k,3)
C          call eua2rm(alpha,beta,gamma,rotm)
C          do  i = 0, 2
C            if (mode == 0) then
C              ammi(1) = rhos(1,i,ilm,1,2,ic) + rhos(1,i,ilm,2,1,ic)
C              ammi(2) = rhos(2,i,ilm,2,1,ic) - rhos(2,i,ilm,1,2,ic)
C              ammi(3) = rhos(1,i,ilm,1,1,ic) - rhos(1,i,ilm,2,2,ic)
C            elseif (mode == 1) then
C              ammi(1) = 0
C              ammi(2) = 0
C              ammi(3) = (qnu(i+1,l+1,1,ic) - qnu(i+1,l+1,2,ic))/(2*l+1)
C            endif
C            call dgemm('T','N',3,1,3,1d0,rotm,3,ammi,3,0d0,amg(1,i),3)
C            call daxpy(3,1d0/(2*l+1),amg(1,i),1,mbar(1,i),1)
C          enddo
C        enddo
C
CC   ... Distribute z part over rhom
C        ilm = j
C        do  m = -l, l
C          ilm = ilm+1
CC         Euler angles for this channel
C          k = 1
C          if (neul == nl) k = l+1
C          if (neul == nl*nl) k = ilm
C          alpha = eulat(k,1)
C          beta  = eulat(k,2)
C          gamma = eulat(k,3)
C          call eua2rm(alpha,beta,gamma,rotm)
C          do  i = 0, 2
C            if (mode == 0) then
C              ql  = rhos(1,i,ilm,1,1,ic) + rhos(1,i,ilm,2,2,ic)
CC             rhom(3) = rhos(1,i,ilm,1,1,ic) - rhos(1,i,ilm,2,2,ic)
C            elseif (mode == 1) then
C              ql  = (qnu(i+1,l+1,1,ic) + qnu(i+1,l+1,2,ic))/(2*l+1)
CC             rhom(3) = (qnu(i+1,l+1,1,ic) - qnu(i+1,l+1,2,ic))/(2*l+1)
C            endif
CC           rhoz_m(loc) = rotm * M
C            rhom(3) = ddot(3,rotm(3,1),3,mbar(1,i),1)
C            qp = (ql+rhom(3))/2
C            qm = (ql-rhom(3))/2
CC            print 533, i,ilm,rhos(1,i,ilm,1,1,ic)-qp,
CC     .                       rhos(1,i,ilm,2,2,ic)-qm
CC  533       format(2i3,2f15.10)
C            rhos(1,i,ilm,1,1,ic) = qp
C            rhos(1,i,ilm,2,2,ic) = qm
C          enddo
C        enddo
C        endif
C
CC   --- Generate qnu+-qnu- from projection of rhos onto Bxc-field ---
CC   ... Projection of m-averaged rhos onto Bxc-field
CC       The coordinate system for rhos for this orbital is specified by
CC       Euler angles (which generate Rm).  Then rhos(loc)=Rm*rhos(glob).
CC       Let Rb = rotation of of Bxc to zhat.  (NB: |Bxc|=1).
CC       Bxc.rhos = Rb*Bxc . Rb*rhos = zhat . (Rb*Rm+) rhos(loc)
CC       The z component of (Rb*(Rm+)) (rho++ - rho--) is qnu+ - qnu-
CC       in the Bxc coordinate system.  Because qnu+ - qnu- is stored
CC       in the local coordinate system of Euler angles, we compute
CC       qnu+ - qnu-  = m-averaged rotation of qnu+ - qnu- in the
CC       Bxc frame: qnu+ - qnu- (loc-frame) = R qnu+ - qnu- (Bxc-frame)
CC       where  R = 1/(2l+1) sum m (Rb*(Rm+))^-1 = 1/(2l+1) sum m Rm*Rb+
CC       Branch has no effect if Bxc || all orbitals in this l channel
C        if (neul /= 1 .or. nrhos > nlo) then
C        ilm = j
C        call dpzero(amli,3)
C        call dpzero(qli,3)
C        rotmb3 = 0
CC       Determine average M for this l, neglecting rhox and rhoy
C        do  m = -l, l
C          ilm = ilm+1
CC         Euler angles for magnetic moment this channel
C          k = 1
C          if (neul == nl) k = l+1
C          if (neul == nl*nl) k = ilm
C          alpha = eula(ib,k,1)
C          beta  = eula(ib,k,2)
C          gamma = eula(ib,k,3)
C          call eua2rm(alpha,beta,gamma,rotm)
CC         Rotation matrix Rb rotates Bxc to zhat
C          if (hatbxc(1)**2+hatbxc(2)**2 /= 0) then
C            ph = datan2(hatbxc(2),hatbxc(1))
C            th = datan2(dsqrt(hatbxc(1)**2+hatbxc(2)**2),hatbxc(3))
C          elseif (hatbxc(3) < 0) then
C            ph = 0
C            th = pi
C          else
C            ph = 0
C            th = 0
C          endif
C          call rotma(ph+pi/2,pi/2,th,rotb)
CC         Rb * Rm+
C          call dgemm('N','T',3,3,3,1d0,rotb,3,rotm,3,0d0,rotmb,3)
CC         ammi: j=3 component of rhom = qnu+-qnu- in Bxc coord sys.
C          if (nrhos == nl*nl) then
C            k = ilm
C            fac = 1
C          else
C            fac = 1/dble(2*l+1)
C            k = l+1
C          endif
C          do  i = 0, 2
C            rhom(1) = (rhos(1,i,k,1,2,ic) + rhos(1,i,k,2,1,ic))*fac
C            rhom(2) = (rhos(2,i,k,2,1,ic) - rhos(2,i,k,1,2,ic))*fac
C            rhom(3) = (rhos(1,i,k,1,1,ic) - rhos(1,i,k,2,2,ic))*fac
C            call dgemm('N','N',3,1,3,1d0,rotmb,3,rhom,3,0d0,ammi,3)
C            if (i == 0) print 521, k, rhom, ammi
C  521       format(i4,3f12.6,2x,3f12.6)
CC           m-averaged rho+ + rho-, rho+ - rho- in Bxc frame
C            qli(i) = qli(i) +(rhos(1,i,k,1,1,ic)+rhos(1,i,k,2,2,ic))*fac
C            amli(i) = amli(i) + ammi(3)
C          enddo
CC         Accumulate sum_m (Rb * Rm+) zhat to rotate amli->loc frame
C          rotmb3 = rotmb3 + rotmb(3,3)/(2*l+1)
C        enddo
CC       Make (qp-qm)_i = amli_i/R_33.  Then
CC       sum Rb Rm+ (qnu+ - qnu-) = m-averaged rhos++-rhos-- in Bxc frame
C        do  i = 0, 2
C          amli(i) = amli(i)/rotmb3
C          qp = (qli(i)+amli(i))/2
C          qm = (qli(i)-amli(i))/2
CC          if (i == 0) print533,i,l,ic,
CC     .      qnu(i+1,l+1,1,ic)+qnu(i+1,l+1,2,ic)-qp-qm,
CC     .      qnu(i+1,l+1,1,ic)-qnu(i+1,l+1,2,ic)-(qp-qm),
CC     .      amli(i)*rotmb3,amli(0),qp,qm
CC  533     format(3i3,6f13.8)
C          qnu(i+1,l+1,1,ic) = qp
C          qnu(i+1,l+1,2,ic) = qm
C        enddo
C      endif
C      enddo
C      endif
C      enddo
C
CC      stop
C
C      end
