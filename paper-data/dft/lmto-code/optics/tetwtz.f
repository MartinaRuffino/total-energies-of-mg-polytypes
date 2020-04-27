      subroutine tetwtz(omg,chipm,npm,linme,esciss,qlat,
     .  efermi,eband1,eband2,nctot,ncc,ecore,
     .  ntetf,nqbzw,nband,nqbz,idtetf,qbzw,ib1bz,
     x  q,iq,isp1,isp2,nqibz,
     o  nwgt,wgt)
C- Weights for dielectric function by tetrahedron method.
C ----------------------------------------------------------------------
Ci Inputs
Ci   chipm :T : to calculate magnetic susceptibility
Ci   npm   :1  use time reversal symmetry
Ci         :2  no time reversal symmetry
Ci         :   swap (ib,jb) => (jb,ib)
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Ci             Not used now
Ci   esciss:energy for scissors operator, unoccupied states
Ci   qlat  :reciprocal lattice vectors
Ci   efermi:Fermi energy  (see Bugs below)
Ci   eband1:eigenvalues at k, in Ry units (see Bugs below)
Ci   eband2:eigenvalues at k+q (see Bugs below)
Ci   nctot :number of cores
Ci   ncc   :number of cores for jb in (ib,jb) pair.  Needed for dimensioning.
Ci         :Should be 0 if time-reversal, nctot otherwise
Ci   ecore :core eigenvalues. (see Bugs below)
Ci         :Not used if nctot=0
Ci   ntetf :number of tetrahedra
Ci   nqbzw :dimensions qbzw,ib1bz. Should be (n1+1)*(n2+1)*(n3+1)
Ci   nband :number of bands
Ci   nqbz  :number of k-points in the 1st BZ
Ci   idtetf:idtetf(1:4,i) points to the 4 k-points defining the ith
Ci         :tetrahedron.  Note: idtetf points to the padded qbzw
Ci   qbzw  :k-points on the regular mesh, padded to repeat points
Ci         :on the boundary lines (tetfbz.f)
Ci   ib1bz :maps the points in the padded qp list to the original list
Ci   q     :(for mtet mode only) q-vector separating occ, unocc states
Ci         :True q = 2*pi*q(1:3)/alat
Ci         :Not used now
Ci   iq    :(mtet mode) used to locate q vector in 1st BZ
Ci         :Not used now
Ci   isp1  :(mtet mode) spin index for occ   state (evals from disk)
Ci         :Also used in chipm
Ci   isp2  :(mtet mode) spin index for unocc state (evals from disk)
Ci         :Not used now
Ci   nqibz :(mtet mode) # kpoints read from disk?
Ci         :Not used now
Co Outputs
Co     wgt :(complex) integration weight for the Lindhard function. complex
Co         :wgt=wgt(band indexfor k, band indexfor q+k,k)
Co    nwgt :the number of non-zero weights for each k.
Cl Local variables
Cl         :
Cr Remarks
Cr   This is an implementation of J.Rath&A.J.Freeman PRB11, 2109 (1975)
Cr   lindtet3 calculates the weights for one tetrahedron.
Cr
Cr *Kotani's notes from tetwt4.F
Cr   The numbering of band index is not unique due to degeneracy.
Cr   It affects how to choose the tetrahedron.  As a result, it affects
Cr   the integration weight, and it might break the crystal symmetry.
Cr   Symmetrization is added at the end to recover the symmetry.
Cr
Cr   subroutine lindtet3 is the kernel, and calculates the microcell integral
Cr    \int d^3k f(e(k)) (1-f(e(q+k))) / (omg - e(q+k) + e(k)),
Cr   where e(q+k) is unoccupied, e(k) is occupied
Cr   and f(E) denotes the Fermi distribution function at T=0.
Cr   omg is complex. Use omg=+1d-8 or so to choose the correct log branch.
Cr
Cr   The Rath&Freeman algorithm requires rather tedious classifications
Cr   because of different possible orderings of the eigenvalues.
Cr   The smoothness might be not fine.
Cr   There is some room to improve this routine, especially in intttvc,
Cr   which is the core part of the tetrahedron method.
Cr   Parameters eps in intttvc might not be optimum.
Cr   These parameters (and classification) are necessary to avoid
Cr   out-of-range arguments, e.g., to avoid divergent denominators in intv4c.
Cr
Cr *mtet mode: doubles the k-mesh in the integration of the delta-function
Cr              (mesh has 2x2x2 more tetrahedra than normal mode).
Cr   Normal  mode requires mtet=(/1,1,1/)
Cr   Doubled mode requires mtet=(/2,2,2/)
Cr   mtet mode is disabled at present
Cb Bugs
Cb   omg  is in Hartree units!
Cu Updates
Cu   25 Apr 14  Adapted from Kotani's tetwt4.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme,chipm
      integer    :: job,nband,nqbz,nctot,ntetf,nqbzw,npm,ncc,
     .  idtetf(0:3,ntetf),ib1bz(nqbzw),nwgt(nqbz,npm)
      real(8) :: efermi,esciss,qlat(3,3),qbzw(3,nqbzw),
     .  eband1(nband,nqbz),eband2(nband,nqbz),ecore(nctot)
      complex(8):: omg, wgt(nband+nctot,nband+ncc,nqbz,npm)
C     mtet mode
      integer   ::isp1,isp2,iq,nqibz
      real(8) :: q(3)
C ... Dynamically allocated arrays
      complex(8),allocatable:: wgt1(:,:)
      real(8),pointer:: eocc(:),eunocc(:)
C ... Local parameters
      integer   :: jpm,ibxmx,jbxmx,jbx,nrankc1,nrankc2,nnn1,nnn2,
     .  itet,ib,jb,iprint,i,j,ik,ibx,nbnc,nb1,nnni,nnnj,
     .  kk(0:3),kkm(0:3),
     .  nbandmx_kq,nbandmx_k,mxbnds,jobl
      integer   :: noccx_k,noccx_kq,irnk1,nrank1,irnk2,nrank2,kx,
     .  im,ini1(nband+nctot),ied1(nband+nctot),ixi1,ixi2,ixe1,ixe2,
     .  ini2(nband+nctot),ied2(nband+nctot)
      real(8) :: det33,kvec(3,0:3),x(0:3),voltot,volt
      real(8),target :: ek(nband+nctot,0:3),ekq(nband+nctot,0:3)
C     integer   ::  noccx1
C     complex(8)::  wmean,a1,a2,wgt2
      logical :: ipr=.false.
      real(8) :: piofvoltot,kkv(3,0:3)
      logical :: chkwrt
      real(8),parameter:: pi=3.1415926535897932d0
      complex(8)::  wgtA,wgtB,wmean,wttx
      real(8),parameter:: eps=0d0 !1d-12 ! cutoff check to determine cancellation.
C ... Needed to synchronize with mtet mode
      logical ::mtett
      integer   :: mtet(3),nmtet,nqbzwm,nqbzm
      integer   ,allocatable:: idtetfm(:,:,:),ib1bzm(:)
      real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
      real(8),allocatable:: ekzz1(:,:),ekzz2(:,:),wtet(:,:,:)

C ... External calls
      external addsciss,chkdgn,lindtet6

C --- Setup ---
      chkwrt = iprint()>=150
      voltot = abs(det33(qlat))
      piofvoltot = pi/4d0/voltot
      jobl = mod(job,10)
      call dpzero(wgt,2*(nband+nctot)*(nband+ncc)*nqbz*npm)

C ... For Divided tetrahedron method
      mtet = (/1,1,1/)
      if (sum(abs(mtet))/=3) then
        mtett = .true.
        call info2(30,0,0,' TETWTZ:  divided-tetrahedron scheme: '//
     .    'mtet = %3;2i',mtet,0)
      else
        mtett = .false.
      endif

      if (mtett) then
        call rx('divided tetrahedron not implemented')
C        ifmtet=501
C        open (ifmtet, file='mtet',form='unformatted')
C        read(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
C        allocate(
C     .       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
C        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
C     .       ib1bzm(nqbzwm), qbzm(3,nqbzm) !,index_qbzm(n_index_qbzm,n_index_qbzm,n_index_qbzm)
C        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
C     .       , wtet(0:3,nmtet,ntetf) )
C        read(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
C        close(ifmtet)
C        ifeig=501
C        open (ifeig,file='eigmtet',form='unformatted')
C        read(ifeig) nband_x,nqnumm,nsp
C        print *,'readin eigmtet ', nband_x,nqnumm,nsp
C        if (nband_x /=nband) then
C          stop 'tetwt5: nband_x /=nband'
C        endif
C        allocate( eigtet(nband,nqnumm,nsp) )
C        do  iqx = 1, nqnumm
C          do  ispx = 1, nsp
C            read (ifeig) eigtet(1:nband,iqx,ispx)
C            if (iprint()>150) write(6,"('iq=',i3,'  eig(1:5)=',5f10.4)")iqx,eigtet(1:5,iqx,ispx)
C          enddo
C        enddo
C        close(ifeig)
      else
        nmtet = 1
        nqbzwm = nqbzw
        nqbzm = nqbz
        allocate(wtet(0:3,nmtet,ntetf),qbzm(3,nqbzm),qbzwm(3,nqbzwm),
     .    idtetfm(0:3,nmtet,ntetf),ib1bzm(nqbzwm))
        idtetfm(:,1,:) = idtetf
        qbzwm = qbzw
        ib1bzm = ib1bz
        wtet = 0.25d0
      endif

c-----------------------------------------------------------------------
      allocate(wgt1(nband+nctot,nband+ncc))
      nbnc = nband+nctot
      nb1 = nband+1
      allocate(ekxx1(nband+nctot,nqbz),ekxx2(nband+nctot,nqbz))
      do  kx = 1, nqbz
        ekxx1( 1:nband, kx) = eband1(1:nband,kx)
        ekxx2( 1:nband, kx) = eband2(1:nband,kx)
        ekxx1( nband+1: nband+nctot, kx) = ecore(1:nctot)
        ekxx2( nband+1: nband+nctot, kx) = ecore(1:nctot)
      enddo

C     print *, 'sumcheck',sum(abs(ekxx1)),sum(abs(eband1))
C     print *, 'sumcheck',sum(abs(ekxx2)),sum(abs(eband2))

C --- Eigenvalues at q and q+k, mtet mode ----
C     At the end of this block:
C     ekzz1 holds bands at k,   for mtett either true or false
C     ekzz2 holds bands at q+k, for mtett either true or false
      allocate(ekzz1(nband+nctot,nqbzm),ekzz2(nband+nctot,nqbzm))
      if (mtett) then
C        call mkqlat(qlat,ginv,volt)
CC        call dinv33x(qlat,ginv)
C        print *, ' mtett mode nqbzm nqibz=',nqbzm,nqibz
C        do  kx = 1, nqbzm
C          ekzz1( 1:nband, kx)              = eigtet(1:nband,kx,isp1)
C          ekzz1( nband+1: nband+nctot, kx) = ecore(1:nctot)
C          ekzz2( nband+1: nband+nctot, kx) = ecore(1:nctot)
C
Cccccccccccccccccccccc
Cc          do  ib = 1, nband
Cc          if (eigtet(ib,kx,isp1)>20) then
Cc            write(6,"('eigtet=',3i5,f13.5)")
Cc     .      ib,kx,isp1,eigtet(ib,kx,isp1)
Cc          endif
Cc          enddo
Cccccccccccccccccccccc
C
C        enddo
C
C         if (iq<=nqibz) then
C          do  kx = 1, nqbzm
C!          q(1:3)+qbzm(1:3,kx) ----> kqxxm
Cc      print *,' xxx3 ginv=',ginv
Cc      print *,' q+qbzm= ', q(1:3) + qbzm(1:3,kx)
Cc          call fbz2 ( q(1:3) + qbzm(1:3,kx),
Cc     i      ginv,index_qbzm,n_index_qbzm,qbzm,nqbzm,
Cc     o      qdummy,kqxxm)
C            kqxxm = iqindx(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm)
Cc          write(6,*) 'kqxxm2 kqxxm=',kqxxm2,kqxxm
Cc          if (kqxxm/=kqxxm2) then
Cc             print *,' q=',q
Cc             stop 'kqxxm/=kqxxm2'
Cc          endif
C            ekzz2(1:nband, kx) = eigtet(1:nband,kqxxm,isp2)
C          enddo
C        else
C          do  kx = 1, nqbzm
C            kp  = nqbzm *(iq - nqibz) + kx
C            ekzz2(1:nband, kx)= eigtet(1:nband,kp,isp2)
C          enddo
Ccccccccccccccccccccccccccccccccccccccccc
CC        print *, 'sumcheck ssssss1=',sum(abs(ekzz2(1:nband,1:nqbzm)))
CC        kp  = nqbzm *(iq - nqibz) + 1
CC        print *,'sumcheck ssssss2=',sum(abs(eigtet(1:nband,kp:kp+nqbzm-1,isp)))
Ccccccccccccccccccccccccccccccccccccccccc
C        endif
C        deallocate( eigtet)
      else
        ekzz1 = ekxx1
        ekzz2 = ekxx2
      endif

C ... Scissors Operator
      if (esciss /= 0d0) then
        call addsciss(esciss,efermi,(nband+nctot)*nqbzm,ekzz1)
        call addsciss(esciss,efermi,(nband+nctot)*nqbzm,ekzz2)
      endif

C      call info5(10,0,0,' tetwtz: nqbz=%i vol=%;6d qlat=%9:1;6d '//
C     .  ' Ef=%;6d nband=%i',nqbz,voltot,qlat,efermi,nband)

C ... Sanity check: tetrahedra volumes must add to cell volume
      volt = 0d0
      do  itet = 1, ntetf
        do  im = 1, nmtet
          kvec(1:3,0:3) = qbzwm(1:3,idtetfm(0:3,im,itet) )
          do  i = 1, 3
            kvec(1:3,i) = kvec(1:3,i) - kvec(1:3,0)
          enddo
          volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
        enddo
      enddo
      if (abs(volt-voltot)>1d-10) then
        call rx('tetwtq: tetrahedra volumes do not '//
     .    'sum to unit cell volume')
      endif

C --- Large loop over tetrahedra ---
      do  itet = 1, ntetf
        kk(0:3) =      ib1bz(idtetf(0:3,itet))     ! Indices to 4 corners
        kkv(1:3,0:3) = qbzw(1:3,idtetf(0:3,itet))  ! k-points at corners
C   ... Loop over multitet
        do  im = 1, nmtet
C         kq (0:3)      = kqxx (kk(0:3)) !  k+q in micro-tet
          kkm(0:3)      = ib1bzm(idtetfm(0:3,im,itet)) !  k   in micro-tet
          kvec(1:3,0:3) = qbzwm ( 1:3, idtetfm(0:3,im,itet) )

          ek ( 1:nband+nctot, 0:3) = ekzz1( 1:nband+nctot, kkm(0:3)) ! k
          ekq( 1:nband+nctot, 0:3) = ekzz2( 1:nband+nctot, kkm(0:3)) ! k+q
C         noccx_k  = noccx1( ek(1:nband, 0:3),4,nband, efermi)  !the highest number of occupied states (HNOS)
C         noccx_kq = noccx1(ekq(1:nband, 0:3),4,nband, efermi)
          noccx_k  = mxbnds(2,ek(1:nband,0:3),nband,nband,4,efermi,i) ! the most number of occupied states
          noccx_kq = mxbnds(2,ekq(1:nband,0:3),nband,nband,4,efermi,i)

C     ... Exclude bands with energies 9999 and higher
          nbandmx_k  = nband
          nbandmx_kq = nband
          do  i = 1, nband
            if (maxval(ek(i, 0:3)) >= 9999) then
              nbandmx_k = i-1
              exit
            endif
          enddo
          do  i = 1, nband
            if (maxval(ekq(i, 0:3)) >= 9999) then
              nbandmx_kq = i-1
              exit
            endif
          enddo

C ...     If no time reversal, loop twice, exchanging i,j second pass
          do  jpm = 1, npm      !feb2006
            if (jpm == 1) then
              ibxmx = noccx_k + nctot
              jbxmx = nbandmx_kq !nband  Apr2009takao
            else
              call rx('tetwtz npm=2 never checked')
              ibxmx = nbandmx_k  !nband  Apr2009takao
              jbxmx = noccx_kq + nctot
            endif
C       --- Loop over all potential (occ, unocc) pairs ---
            do  ibx  = 1, ibxmx  !noccx_k + nctot  !   occupied
            do  jbx  = 1, jbxmx  !nband            ! unoccupied
C             Get true ib,jb for this jpm
              if (ibx<=noccx_k .or. jpm==2) then
                ib = ibx
              else
                ib = ibx - noccx_k + nband
              endif
              if (jbx<=noccx_kq .or. jpm==1) then
                jb = jbx
              else
                jb = jbx - noccx_kq + nband
              endif
C             job=10 case: jb must be 1
              if (job >= 10 .and. jb /= 1) cycle

              if (jpm == 1) then
C               ibxmx = noccx_k + nctot
C               jbxmx = nband
                eocc   => ek (ib,0:3)
                eunocc => ekq(jb,0:3)
              else
c               ibxmx = nband
c               jbxmx = noccx_kq + nctot
                eunocc => ek (ib,0:3)
                eocc   => ekq(jb,0:3)
              endif

              if (minval(eocc)   <= efermi .and.
     .            maxval(eunocc) >= efermi) then
                if (maxval(eunocc(:)-eocc(:)) <0) cycle ! eliminate cases with no contribution

                x(0:3) = .5d0*(eocc-eunocc) ! + omg !Denominator. unit in Hartree.
                if (chkwrt) then
                  write(6,"(' ... Before lindtet3: itet ib jb Ef='
     .              ,i8,2i5,d12.4,' ###')") itet,ib,jb,efermi
                  write(6,"('  eocc  - Ef= ',4f10.3)") eocc  -efermi
                  write(6,"('  eunocc- Ef= ',4f10.3)") eunocc-efermi
                  write(6,"('  -x(a.u.)= ',4f10.3)") -x(0:3)
                endif

                if (chipm .and. isp1==2) then
                  wgtA = 0d0
                else
                  call lindtet3(omg,kvec,ek(ib,0:3),ekq(jb,0:3),x,
     .              efermi,wttx) ! kvec(1:3, 0:3), ea(0:3), x(0:3)
                  wgtA = wttx/4d0/voltot
                endif

                if (chipm .and. isp1==1) then
                  wgtB = 0d0
                elseif (abs(omg)<1d-10.and.(.not.chipm)) then ! static case
                  wgtB = wgtA
                else
                  call lindtet3(-omg,kvec,ek(ib,0:3),ekq(jb,0:3),x,
     .              efermi,wttx)
                  wgtB = wttx/4d0/voltot
                endif
                wgt1(ib,jb) = wgtA + wgtB

C              if (chkwrt) then
C                write(6,"('ttt',3i4,f10.5,d12.4)") itet,ib,jb,dble(omg),-dimag(wttx/pii)
C              endif


              endif

            enddo
            enddo

          wgt(1:noccx_k,:,kk(0),jpm) = wgt(1:noccx_k,:,kk(0),jpm)
     .        + wgt1(1:noccx_k,:)* 4*wtet(0,im,itet)
          wgt(nb1:nbnc,: ,kk(0),jpm) = wgt(nb1:nbnc,: ,kk(0),jpm)
     .      + wgt1(nb1:nbnc, :)* 4*wtet(0,im,itet)
          wgt(1:noccx_k,:,kk(1),jpm) = wgt(1:noccx_k,:,kk(1),jpm)
     .      + wgt1(1:noccx_k,:)* 4*wtet(1,im,itet)
          wgt(nb1:nbnc,: ,kk(1),jpm) = wgt(nb1:nbnc,: ,kk(1),jpm)
     .      + wgt1(nb1:nbnc, :)* 4*wtet(1,im,itet)
          wgt(1:noccx_k,:,kk(2),jpm) = wgt(1:noccx_k,:,kk(2),jpm)
     .      + wgt1(1:noccx_k,:)* 4*wtet(2,im,itet)
          wgt(nb1:nbnc,: ,kk(2),jpm) = wgt(nb1:nbnc,: ,kk(2),jpm)
     .      + wgt1(nb1:nbnc, :)* 4*wtet(2,im,itet)
          wgt(1:noccx_k,:,kk(3),jpm) = wgt(1:noccx_k,:,kk(3),jpm)
     .      + wgt1(1:noccx_k,:)* 4*wtet(3,im,itet)
          wgt(nb1:nbnc,: ,kk(3),jpm) = wgt(nb1:nbnc,: ,kk(3),jpm)
     .      + wgt1(nb1:nbnc, :)* 4*wtet(3,im,itet)
          enddo  ! loop over jpm
        enddo ! loop over multitet
      enddo ! loop over tetrahedra
      deallocate(idtetfm, qbzwm,ib1bzm, qbzm)

c      print *, 'sumcheck wgt=', sum(wgt)

C --- Symmetrization of wgt ---
C     Average contributions from degenerate levels.
      if (npm == 2)
     .  call rx('tetwtz symmetrization not ready for npm=2')
      do  kx = 1, nqbz
C       Group occ,unocc bands into families of degenerate levels
        call chkdgn(ekxx1(:,kx),nband,nrank1,ini1,ied1,0,ipr)
        call chkdgn(ekxx2(:,kx),nband,nrank2,ini2,ied2,0,ipr)
C       Ditto for the core levels
        nrankc1 = 0
        if (nctot /= 0) then
          call chkdgn(ecore,nctot,nrankc1,ini1(nrank1+1),
     .      ied1(nrank1+1),nband,ipr)
        endif
        nrankc2 = 0
        if (nctot /= 0) then
          call chkdgn(ecore,nctot,
     .    nrankc2,ini2(nrank2+1),ied2(nrank2+1),nband,ipr)
        endif
        do  jpm = 1, npm
          if (jpm == 1) then
            nnn1 = nrankc1
            nnn2 = 0
          else
            nnn1 = 0
            nnn2 = nrankc2
          endif

C     ... Loop over nrank + nnn distinct valence + core levels
          do  irnk1 = 1, nrank1 + nnn1
            do  irnk2 = 1, nrank2 + nnn2
              ixi1  = ini1(irnk1);  ixe1 = ied1(irnk1)
              ixi2  = ini2(irnk2);  ixe2 = ied2(irnk2)

              wmean = sum(wgt(ixi1:ixe1,ixi2:ixe2,kx,jpm))
     .              / ((ixe1-ixi1+1)*(ixe2-ixi2+1))
              wgt(ixi1:ixe1,ixi2:ixe2,kx,jpm) = wmean
            enddo
          enddo
        enddo ! end of jpm loop
      enddo ! end of kx loop

C --- Make nwgt = number of (ib,jb) pairs for a given k ---
C      if (jobl==0) then
        do  jpm = 1, npm
          do  ik = 1, nqbz
            nwgt(ik,jpm) = 0
            if (jpm == 1) then
              nnni = nctot
              nnnj = 0
            else
              nnni = 0
              nnnj = nctot
            endif
            do  i = 1, nband+nnni
              do  j = 1, nband+nnnj
C               if (iwgt(i,j,ik,jpm)) then
                if (abs(wgt(i,j,ik,jpm)) > eps) then
                  nwgt(ik,jpm) = nwgt(ik,jpm)+1
                endif
              enddo
            enddo
          enddo
        enddo
c       write(6,*)' max num of nonzero wgt(k)=',maxval(nwgt)
c       write(6,*)' tot num of nonzero wgt   =',sum(nwgt)
c       write(6,*)' sum of wgt   =',sum(wgt)
c       write(6,*)' tetwt4: end '; call cputid  (0)
        call info2(20,0,0,' TETWTZ:  total number of (occ,unocc) '//
     .  'pairs = %i, %i or fewer per qp',
     .    sum(nwgt),maxval(nwgt(1:nqbz,1:npm)))
        call rx('another check?')
C        if (sum(nwgt) /= count(iwgt))
C     .    call rx(' bug in tetwtz: pairs improperly counted')
C      endif

C --- Cleanup ---
      deallocate(wgt1)
      if (allocated(idtetfm)) deallocate(idtetfm)
      if (allocated(ib1bzm)) deallocate(ib1bzm)
      if (allocated(qbzm)) deallocate(qbzm)
      if (allocated(qbzwm)) deallocate(qbzwm)
      if (allocated(ekxx1)) deallocate(ekxx1)
      if (allocated(ekxx2)) deallocate(ekxx2)
      if (allocated(ekzz1)) deallocate(ekzz1)
      if (allocated(ekzz2)) deallocate(ekzz2)
      if (allocated(wtet)) deallocate(wtet)
      if (allocated(idtetfm)) deallocate(idtetfm)
      if (allocated(ib1bzm)) deallocate(ib1bzm)
      if (allocated(qbzm)) deallocate(qbzm)
      if (allocated(qbzwm)) deallocate(qbzwm)
      if (allocated(ekxx1)) deallocate(ekxx1)
      if (allocated(ekxx2)) deallocate(ekxx2)
      if (allocated(ekzz1)) deallocate(ekzz1)
      if (allocated(ekzz2)) deallocate(ekzz2)
C     if (allocated(eigtet)) deallocate(eigtet)
      if (allocated(wtet)) deallocate(wtet)
      end

      subroutine lindtet3(omg,kvec,ea,eb,x,efermi,wttx)
C- Calculate the integral of \int dk1 dk2 dk3 f(ea)(1-f(eb))/x
C ----------------------------------------------------------------------
Ci Inputs
Ci   omg   :frequency (complex)
Ci   kvec  :wave numbers at the four corners of the tetrahedron
Ci   ea    :energy at four corners of tetrahedron at k
Ci   eb    :energy at four corners of tetrahedron at k+q
Ci   x     :energy demoninator entering into delta-function integral
Ci   efermi:Fermi energy
Co Outputs
Co   wttx  :weight from this tetrahedron
Cs Command-line switches
Cl Local variables
Cr Remarks
Cr   f(E) denotes the Fermi distribution function.
Cr   For T=0 only
Cr   Tetrahedron is specified by values at 4 corners
Cr     kvec(1:3,1:4), ea(1:4), eb(1:4), x(1:4)
Cr   This code is based on J.Rath & A.J.Freeman PRB11(6) p.2109(1975).
Cu Updates
Cu   25 Apr 14  Adapted from gwsrc/tetwt4.F
C---------------------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8)    :: kvec(1:3,1:4),ea(1:4),eb(1:4),x(1:4),efermi
      complex(8) :: omg,wttx
C ... Local parameters
      integer    :: ieaord(1:4),i,isig,n,itmp,ix
      real(8) :: kk(3,4),xx(4),ee(4),ebf(4),ebfKx(4),Kx(3,4),xKx(4)
      complex(8) :: inttetra3
      integer :: iprint
      logical :: debug
C     integer stdo,nglob

      debug = .false.
C     if (iprint()>=110) debug = .true.
C     stdo = nglob('stdo')
      if (debug) then
        write(996,*)' lindtet3: *****************************'
        write(996,"(' i=',i3,' x ea eb =',3f10.5)")
     .  (i, x(i), ea(i)-efermi, eb(i)-efermi,i=4,1,-1)
      endif

C ... Order Eigenvalues
      n = 4
      isig = 1
      do  i = 1, n
        ieaord(i) = i
      enddo
      do  ix = 2, n
        do  i = ix, 2,-1
          if (ea(ieaord(i-1)) > ea(ieaord(i))) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
            isig = -isig
            cycle
          endif
          exit
        enddo
      enddo
      ieaord(1:4) = ieaord(4:1:-1)

      ! ordering of E4,E3,E2,E1 This suits for Rath&Freeman.
      kk(1:3,1:4) = kvec(1:3,ieaord(1:4))
      ! 4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
      ee (1:4) = ea (ieaord(1:4)) - efermi
      xx (1:4) = x  (ieaord(1:4))
      ebf(1:4) = eb (ieaord(1:4)) - efermi
ccccccccccccc
      if (debug) then
        write(996,"(' i iea=',2i4)") (i,ieaord(i),i=1,4)
        write(996,"(' i=',i3,' xx ee ebf =',3f10.5,'  kk=',3f10.5)")
     .     (i,xx(i),ee(i),ebf(i),kk(1:3,i),i=4,1,-1)
      endif
ccccccccccccc
      if (0d0 <= ee(4)) then
        wttx = (0d0,0d0)
      elseif (ee(4)<0d0 .and. 0d0<=ee(3)) then   !!! Fig 1.
        if (debug) write(996,*) 'lindtet3: fig 1 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  Kx(1,1),xKx(1),ebfKx(1))
        !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk3(kk,ee,xx,ebf, 4,1,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 4,3,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4,10,20,30/) )  ! k4,K1,K2,K3
      elseif (ee(3)<0d0 .and. 0d0<=ee(2)) then   !!! Fig 2.
        if (debug) write(996,*) 'lindtet3: fig 2 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  Kx(1,1),xKx(1),ebfKx(1)) !K1
        call midk3(kk,ee,xx,ebf, 4,1,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 3,1,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        call midk3(kk,ee,xx,ebf, 3,2,  Kx(1,4),xKx(4),ebfKx(4)) !K4
ccccccccccc
c         write(stdo,*)'xxxx inttetra='
c     .       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 3,10,20/))  ! k4,k3,K1,K2
c         write(stdo,*) omg,kk,xx,ebf
c         write(stdo,*) ' Kx=',Kx,xKx,ebfKx
c         write(stdo,*)'xxxx inttetra='
c     .       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,20,30,10/))  ! k3,K2,K3,K1
c         write(stdo,*) omg,kk,xx,ebf
c         write(stdo,*) ' Kx=',Kx,xKx,ebfKx
c         write(stdo,*)'xxxx inttetra='
c     .       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,10,30,40/))  ! k3,K1,K3,K4
c         write(stdo,*) omg,kk,xx,ebf
c         write(stdo,*) ' Kx=',Kx,xKx,ebfKx
c     stop ' ***************** end test ttt *****************'
ccccccccccc
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 3,10,20/))  ! k4,k3,K1,K2
     .       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,20,30,10/))  ! k3,K2,K3,K1
     .       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,10,30,40/))  ! k3,K1,K3,K4
      elseif (ee(2)<0d0 .and. 0d0<=ee(1)) then   !!! Fig 3.
        if (debug) write(*,*) 'lindtet3: fig 3 xxx'
        call midk3(kk,ee,xx,ebf, 1,4,  Kx(1,1),xKx(1),ebfKx(1)) !K1
        call midk3(kk,ee,xx,ebf, 1,2,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 1,3,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3, 4,30, 2/))  ! k3,k4,K3,k2
     .       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4,10,20,30/))  ! k4,K1,K2,K3
     .       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 2,20,30/))  ! k4,k2,K2,K3
cccccccccc test fulled cccccccccccccccccccccccccccccccccccc
cc     .        + inttetra(kk,xx,Kx,xKx, (/1, 10,20,30/) )  ! test fulling
cc           print *,' fulled test'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
        if (debug) write(996,*) 'lindtet3: fig 4 xxx'
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfkx,(/1,2,3,4/))  ! k1,k2,k3,k4
      endif
      end

C --- Function inttetra3 ---
      complex(8) function
     .  inttetra3(omg,kk_,xx_,ebf,Kx_,xKx_,ebfKx,itetx)
C- Calculate tetrahedron integral Rath & Freeman Eq.(16)
C---------------------------------------------------------------------------------
Ci   omg   :frequency (complex)
Ci    kk_  :(k1-k4) in Rath & Freeman paper
Ci    Kx_  :(K1-K4) in Rath & Freeman paper
Ci    xx_  :xx_ + omg = energy denominator
Ci   xKx_  :xx_ corresponding to Kx_
Ci  itetx  :denote the four corners of the tetrahedron
Cr Remarks
Cr The four corners are selected from 8 points. It is specified by itetx.
C---------------------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: itetx(4)
      complex(8) :: omg
      real(8) :: kk_(3,4),xx_(4),ebf(4),Kx_(3,4),xKx_(4),ebfKx(4)
C ... Local parameters
      integer ix,i,ieaord(4),isig,n,itmp
      real(8) :: kk(3,4),xx(4),ebfin(4),Kx(3,4),xKx(4),ee(4),
     .  kin(3,4),xin(4)
      complex(8) :: wttx,inttetrac
      logical :: debug=.false.

C ... kk xin ebfin
      do  i = 1, 4
        if (itetx(i) < 10) then   ! kk (k1-k4 in Paper)
          ix = itetx(i)
          xin(    i) = xx_(    ix)
          ebfin(i)   = ebf(ix)
          kin(1:3,i) = kk_(1:3,ix)
        else                       ! Kx (K1-K4 in Paper)
          ix = itetx(i)/10
          xin(    i) = xKx_(    ix)
          ebfin(i)   = ebfKx(ix)
          kin(1:3,i) = Kx_ (1:3,ix)
        endif
      enddo
c      print *, ' ebf  =', ebf
c      print *, ' ebfkx=', ebfkx
c      write(*,"(' i=',i3,' xin ebfin =',2f10.5,'  kin=',3f10.5)")
c     .  (i,xin(i),ebfin(i),kin(1:3,i),i=4,1,-1)

C ... Order Eigenvalues in descending order ee(4)>ee(3)>ee(2)>ee(1)
      n = 4
      isig = 1
      do  i = 1, n
        ieaord(i) = i
      enddo
      do  ix = 2, n
        do  i = ix, 2,-1
          if (ebfin(ieaord(i-1)) > ebfin(ieaord(i))) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      kk(1:3,1:4) = kin  (1:3,ieaord(1:4))   ! 4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
      ee(    1:4) = ebfin(    ieaord(1:4))
      xx(    1:4) = xin  (    ieaord(1:4))
cccc
      if (debug) then
        write(996,*)' inttetra3: ***** '
        write(996,"(' i=',i3,' xx ee =',2f10.5,'  kk=',3f10.5)")
     .  (i,xx(i),ee(i),kk(1:3,i),i=4,1,-1)
      endif
cccc
      if (0d0>=ee(4)) then
        wttx = (0d0,0d0)
      elseif (ee(4)>0d0 .and. 0d0>=ee(3)) then   !!! Fig 1.
        if (debug)  write(996,*)' intterra3: fig1'
        call midk(kk,ee,xx, 4,2,  Kx(1,1),xKx(1)) !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk(kk,ee,xx, 4,1,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 4,3,  Kx(1,3),xKx(3)) !K3
        wttx = inttetrac(omg,kk,xx,Kx,xKx,  (/4,10,20,30/))  ! k4,K1,K2,K3
      elseif (ee(3) > 0d0 .and. 0d0>= ee(2)) then   !!! Fig 2.
        if (debug)  write(996,*)' intterra3: fig2'
        call midk(kk,ee,xx, 4,2,  Kx(1,1),xKx(1)) !K1
        call midk(kk,ee,xx, 4,1,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 3,1,  Kx(1,3),xKx(3)) !K3
        call midk(kk,ee,xx, 3,2,  Kx(1,4),xKx(4)) !K4
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/4, 3,10,20/))  ! k4,k3,K1,K2
     .        + inttetrac(omg,kk,xx,Kx,xKx, (/3,20,30,10/))  ! k3,K2,K3,K1
     .        + inttetrac(omg,kk,xx,Kx,xKx, (/3,10,30,40/))  ! k3,K1,K3,K4
      elseif (ee(2) > 0d0 .and. 0d0>= ee(1)) then   !!! Fig 3.
        if (debug)  write(996,*)' intterra3: fig3'
        call midk(kk,ee,xx, 1,4,  Kx(1,1),xKx(1)) !K1
        call midk(kk,ee,xx, 1,2,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 1,3,  Kx(1,3),xKx(3)) !K3
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/3, 4,30, 2/))  ! k3,k4,K3,k2
     .        + inttetrac(omg,kk,xx,Kx,xKx, (/4,10,20,30/))  ! k4,K1,K2,K3
     .        + inttetrac(omg,kk,xx,Kx,xKx, (/4, 2,20,30/))  ! k4,k2,K2,K3
      else
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/1,2,3,4/) )  ! k1,k2,k3,k4
        if (debug)  write(996,*)' intterra3: fig4'
      endif
      inttetra3 = wttx
      end

C --- Function inttetrac ---
      complex(8) function inttetrac(omg,kk,xx,Kx,xKx,itetx)
C- Calculate tetrahedron integral Rath & Freeman Eq.(16).
C---------------------------------------------------------------------------------
Ci Inputs
Ci   omg   :frequency (complex)
Ci    kk   :(k1-k4) and xx (value of denominator)
Ci    Kx   :(K1-K4) and xkx
Ci  itetx  :denote the four corners of the tetrahedron
Cr Remarks
Cr The four corners are selected from 8 points. It is specified by itetx.
C---------------------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: itetx(4)
      real(8) :: kk(3,4),xx(4),Kx(3,4),xKx(4)
      complex(8) :: omg
C ... Local parameters
      integer   :: ix,i
      real(8) ::  am(3,3),kin(3,4),xin(4),det33
      complex(8) :: intttvc2
      logical:: debug

      debug = .false.

      if (omg==1d30 .or. omg==-1d30) then
        inttetrac = 1d0 !this mode is only for cont the non-zero pairs, occu and unocc
        return
      endif
      do  i = 1, 4
        if (itetx(i) < 10) then   ! kk (k1-k4 in Paper)
          ix = itetx(i)
          xin(    i) = xx(    ix)
          kin(1:3,i) = kk(1:3,ix)
        else                      ! Kx (K1-K4 in Paper)
          ix = itetx(i)/10
          xin(    i) = xKx(    ix)
          kin(1:3,i) = Kx (1:3,ix)
        endif
      enddo

      do  i = 1, 3
        am(1:3,i) = kin(1:3,i) - kin(1:3,4)
      enddo
      inttetrac = intttvc2(omg,xin) * abs(det33(am)/6d0) ! \omega (volume of tetrahedra) = abs(det33(am)/6d0) See Eq. (17).

      if (debug) then
        write(996,"(' omg  abs(det33(am)/6d0) =',2d13.5,3x,4d13.5)")
     .     omg,abs(det33(am)/6d0)
        write(996,"(' xin       =',4d13.5)") xin
        write(996,"(' inttetrac =',4d13.5)") inttetrac
        write(996,*)
      endif

      end

C --- Function intttvc2 ---
      complex(8) function intttvc2(omg,v)
C- Tetrahedron integral Rath & Freeman Eq.(16) (except volume factor)
C---------------------------------------------------------------------------------
Ci   omg   :frequency (complex)
Ci   v(1:4):v+omg = energy demoninator
Cr Remarks
Cr   Include imaginary part
C---------------------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8) :: v(4)
      complex(8):: omg
C ... Local parameters
      complex(8) :: vvvc,intvvc,intv4c,intv3c
C     complex(8) :: intv2c
      integer    :: ieaord(1:4),i,isig,idf,ix,n,itmp
      real(8)    :: vv(4),v3,v4,vvd,vvas,vvv,diffc
      logical ::debug=.false.
      real(8),parameter :: eps = 1d-3  ! Judge coincidence. Safer than 1d-4

C ... Order eigenvalues
      n=4
      isig = 1
      do  i = 1, n
        ieaord(i) = i
      enddo
      do  ix = 2, n
        do  i = ix, 2,-1
          if ( v(ieaord(i-1)) > v(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      vv(1:4) = v(ieaord(1:4)) !vv is ordered

      vvas = sum(abs(vv(1:4)))   ! ;print *,'vvas=',vvas
      diffc = max(eps,eps*vvas)
      idf = 0
      do  i = 1, 3
        if ( abs( vv(i)-vv(i+1) ) < diffc ) idf = idf+10**(3-i)
      enddo

      if (debug) then
        write(996,"(' xxx vv=',4d23.15,' idf= ',i3)")
     . vv(1),vv(2),vv(3),vv(4),idf
      endif

c      if (   abs( vv(1)-vv(4)   ) < eps2*vvas+1d-4)  idf = 111
      if (    idf==111 ) then  !2002 apr
        vvvc = (vv(1)+vv(2)+vv(3)+vv(4))/4
        vv(1) = vvvc -.5d0*diffc
        vv(2) = vvvc -.5d0*diffc
        vv(3) = vvvc +.5d0*diffc
        vv(4) = vvvc +.5d0*diffc
        idf = 101
      elseif (idf==011) then
        vv(4) = vv(3) + diffc
        idf = 010
      elseif (idf==110) then
        vv(1) = vv(2) -diffc
        idf = 010
      endif

      if (debug) then
        write(996,"(' yyy vv=',4d23.15,' idf= ',i3)")
     . vv(1),vv(2),vv(3),vv(4),idf
      endif

! Possible values of idf are denoted by ijk, where i,j and k can be 0 or 1.
! 1 means corresponding vv are the same.
      intttvc2 = (0d0,0d0)
      if (    idf==  000 ) then
        intttvc2 = intvvc(omg,vv)
      elseif (idf==  100 ) then
        vvv = (vv(1)+vv(2))/2
        v3  = vv(3)
        v4  = vv(4)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif (idf==  010 ) then
        vvv = (vv(2)+vv(3))/2
        v3  = vv(1)
        v4  = vv(4)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif (idf==  001 ) then
        vvv = (vv(3)+vv(4))/2
        v3  = vv(1)
        v4  = vv(2)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif (idf==  101 ) then
        vvv = (vv(1)+vv(2))/2
        vvd = (vv(3)+vv(4))/2
        intttvc2 = intv3c(omg,vvv,vvd)
      else
        stop 'intttvc2: idf wrong...'
c      elseif (idf==  110 ) then
c        vvv = (vv(1)+vv(2)+vv(3))/3
c        v4  = vv(4)
c        if ( abs(vvv+omg)>eps3 ) intttvc2 = intv2c(omg,vvv,v4)
c      elseif (idf==  011 ) then
c        vvv = (vv(2)+vv(3)+vv(4))/3
c        v4  = vv(1)
c        if ( abs(vvv+omg)>eps3 ) intttvc2 = intv2c(omg,vvv,v4)
c      elseif (idf==  111 ) then
c        vvvc = (vv(1)+vv(2)+vv(3)+vv(4))/4 + omg
c        if ( abs(vvvc)>eps3 ) intttvc2 = 1d0/vvvc
      endif
      end

C --- Function intvvc ---
      complex(8) function intvvc(omg,vv)
C- ?
      implicit none
C ... Passed parameters
      real(8) :: vv(4)
      complex(8):: omg
C ... Local parameters
      complex(8):: v1c,v2c,v3c,v4c,intvvcc
      real(8) :: v1,v2,v3,v4,d1,d2,d3

      v1=vv(1); v2=vv(2); v3=vv(3); v4=vv(4)
      if (v1 == 0) v1 = 1d-15
      if (v2 == 0) v2 = 1d-15
      if (v3 == 0) v3 = 1d-15
      if (v4 == 0) v4 = 1d-15
      d1 = (v1-v4)*(v1-v2)*(v1-v3)
      d2 = (v2-v4)*(v2-v3)*(v2-v1)
      d3 = (v3-v4)*(v3-v1)*(v3-v2)
      v1c = v1+omg
      v2c = v2+omg
      v3c = v3+omg
      v4c = v4+omg
      intvvcc = v1c**2/d1*log(v1c/v4c)  ! rath&freeman eq.(17)
     .        + v2c**2/d2*log(v2c/v4c)
     .        + v3c**2/d3*log(v3c/v4c)
      intvvc = 3*intvvcc
      end

C --- Function intv4c ---
      complex(8) function intv4c(omg,v,v3,v4) !case (iv) at page 2113.
C-
      implicit none
C ... Passed parameters
      complex(8):: omg
      real(8) :: v,v3,v4
C ... Local parameters
      complex(8):: v3c,v4c,vc,xfun,x3,x4
C     complex(8):: intv4cc

      if (v  == 0) v  = 1d-15
      if (v3 == 0) v3 = 1d-15
      if (v4 == 0) v4 = 1d-15
      v3c = V3+omg
      v4c = V4+omg
      vc  = V +omg

C     case4
      x3 = ( V -V3 )/V3c
      x4 = ( V -V4 )/V4c
      intv4c = 3d0*(xfun(x3) - xfun(x4))/(v3-v4)

c      write(996,"('intv4c: case4 x3',2d13.5)") x3
c      write(996,"('              x4',2d13.5)") x4

c--- case1
cc      write(996,"('intv4c: case1'))")
c      intv4cc= V3c**2/(v3-v)**2/(v3-v4)*log(V3c/Vc)
c     .     +   V4c**2/(v4-v)**2/(v4-v3)*log(V4c/Vc)
c     .     +   Vc/(v3-v)/(v4-v)
c      intv4c = 3*intv4cc

c--- case2
c      write(996,"('intv4c: case2'))")
c      intv4cc= 1d0/(v3-v4) * (
c     .     1d0/(v3-v) * ( V3c**2* log(V3c/Vc)/(v3-v) - Vc )
c     .   - 1d0/(v4-v) * ( V4c**2* log(V4c/Vc)/(v4-v) - Vc )
c     .         )
c      intv4c = 3*intv4cc

c---case3
c      write(996,"('intv4c: case3'))")
c      intv4c= 3d0/(v3-v4) * (
c     .   V3c**2/(v3-v)**2 * (log(V3c/Vc) -1d0 + vc/v3c)
c     . - V4c**2/(v4-v)**2 * (log(V4c/Vc) -1d0 + vc/v4c)
c     .         )
c      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(996,"('intv4c-1:', 2d13.5))")
c     .   1d0/(v3-v4) *( 1d0/(v3-v) * ( V3c**2* log(V3c/Vc)/(v3-v) - Vc))
c      write(996,"('intv4c-2:', 2d13.5))")
c     .   1d0/(v3-v4) *( 1d0/(v4-v) * ( V4c**2* log(V4c/Vc)/(v4-v) - Vc))
c      write(996,"('intv4c:', 2d13.5))")
c     .        V4c**2/(v4-v)**2/(v4-v3)*log(V4c/Vc)
c      write(996,"('intv4c:', 2d13.5))")
c     .        Vc/(v3-v)/(v4-v)
cccccccccccccccccccccccccccccccccccccccccccccdddddddd
      end

C --- Function intv3c ---
      complex(8) function intv3c(omg,v,vd)
C-
      implicit none
C ... Passed parameters
      complex(8):: omg
      real(8) :: v,vd
C ... Local parameters
      complex(8):: vc,vdc,x,xfun2
C     complex(8):: intv3cc,intv3cy,xfun3
c     integer   ,save :: ic=0

      if (v  == 0) v =1d-15
      if (vd == 0) vd=1d-15

C     case2 numerically better
      vc  = V +omg
      vdc = Vd+omg
      x = vdc/vc - 1d0
      intv3c = 3*vc*xfun2(x)/(v-vd)**2

C      case1
c      vc  = V +omg
c      vdc = Vd+omg
c      intv3cc = (2*vc*vdc/(v-vd)*log(vdc/vc) + (vc+vdc))
c     .          /(v-vd)**2
c      intv3c = 3*intv3cc

cccccccccccccccccccccccccccccccccccccccccccccccccc
c      ic=ic+1
c      intv3cy = 3*vc*xfun3(x)/(v-vd)**2
c      write(996,"('3',i7,2d12.4,2f10.5))") ic,omg,v,vd
c      write(996,"('3',i7,2d16.8))") ic,intv3c
c      write(996,"(8x,2d16.8))") intv3cx
c      write(996,"(8x,2d16.8))") intv3cy
c      write(996,*)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      end

C --- Function intv2c ---
      complex(8) function intv2c(omg,v,v4)
C-
      implicit none
C ... Passed parameters
      complex(8):: omg
      real(8) :: v,v4
C ... Local parameters
      complex(8):: vc,v4c,intv2cc

      if (v  == 0) v =1d-15
      if (v4 == 0) v4=1d-15
      vc  = V +omg
      v4c = V4+omg
      intv2cc = v4c**2/(v-v4)**3*log(vc/v4c)
     .     + (1.5d0*v4c**2 + .5d0*vc**2 - 2*vc*v4c) / (v-v4)**3
      intv2c = 3*intv2cc
      end

C --- Function xfun ---
      complex(8) function xfun(x)
C-
      implicit none
C ... Passed parameters
      complex(8) :: x
C ... Local parameters
      integer n
      complex(8) :: a

      if (abs(x)<1d-1) then
        a = -x
        xfun = 0d0
        do  n = 20, 3,-1
          xfun = xfun*a + 1d0/n !honar procedure
        enddo
        xfun = xfun*a
      else
        xfun = ( (1d0-x/2d0)*x - log(1d0+x) )/x**2
      endif
      end

C --- Function xfun2 ---
      complex(8) function xfun2(x)
C-
      implicit none
      complex(8) x
      integer n
      complex(8) a

      if (abs(x)<1d-1) then
        a = -x
        xfun2 = 0d0
        do  n = 20, 3, -1
          xfun2 = xfun2*a + 1d0/n !honar procedure
        enddo
        xfun2 = (-xfun2* 2*(1+x) +1d0 )* x**2
      else
        xfun2 =  x + 2 - 2*(1+x)*log(1+x)/ x
      endif
      end

