      subroutine tetwtq(job,npm,ncc,linme,esciss,qlat,fac,
     .  efermi,eband1,eband2,nctot,ecore,
     .  ntetf,nqbzw,nband,nqbz,idtetf,qbzw,ib1bz,
     .  frhis,nwhis,nwgtx,ibjb,nhwtot,ihw,nhw,jhw,
     x  q,iq,isp1,isp2,nqibz,
     o  iwgt,nwgt,demin,demax,
     o  whw)
C- Weights for dielectric function by the tetrahedron method.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0: return iwgt,nwgt,demin,demax only
Ci         :1: make (ib,jb) and k-resolved weights whw
Ci         :   Add 10 to restrict jb to 1 (for single DOS)
Ci   npm   :1 : use time reversal symmetry
Ci         :2 : no time reversal symmetry
Ci   ncc   :number of cores for jb in (ib,jb) pair:
Ci         :Should be 0 if time-reversal, nctot otherwise
Ci   eband1:eigenvalues at k, in Ry units (see Bugs below)
Ci   eband2:eigenvalues at k+q (see Bugs below)
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Ci   esciss:energy for scissors operator, unoccupied states
Ci   qlat  :reciprocal lattice vectors
Ci   efermi:Fermi energy  (see Bugs below)
Ci   ntetf :number of tetrahedra
Ci   nqbzw :dimensions qbzw,ib1bz. Should be (n1+1)*(n2+1)*(n3+1)
Ci   nband :number of bands
Ci   nqbz  :number of k-points in the 1st BZ
Ci   nctot :number of cores
Ci   ecore :core eigenvalues. (see Bugs below)
Ci         :Not used if nctot=0
Ci   idtetf:idtetf(1:4,i) points to the 4 k-points defining the ith
Ci         :tetrahedron.  Note: idtetf points to the padded qbzw
Ci   qbzw  :k-points on the regular mesh, padded to repeat points
Ci         :on the boundary lines (tetfbz.f)
Ci   ib1bz :maps the points in the padded qp list to the original list
Ci   iwgt  :T: ib(occ) jb(unocc) pair makes nonzero contribution
Ci         :Output when job=0
Ci         :Input  when job=1
Ci   nwgt  :number of ib(occupied) jb(unoccupied) pairs for this kpoint
Ci         :Output when job=0
Ci         :Input  when job=1
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram (See Bugs below)
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci         :Not used when 1s digit job=0
Ci         :Note!! frhis is in Hartree, while eband are in Ry
Ci   nwhis :number of energies in histogram
Ci         :Not used when job=0
Ci   nwgtx :leading dimension of ihw,jhw,nhw
Ci         :should be maxval(nwgt)
Ci         :Not used when job=0
Ci   ibjb  :ibjb(nctot+nband,nband,nqbz) = ibjb index for given ib jb k.
Ci   nhwtot:dimension of whw
Ci   ihw   :ihw(ibjb,k) = index to first histogram bin encompassing
Ci                        (demin,demax) for a given ibjb pair and k
Ci         :Not used when job=0
Ci   nhw   :nhw(ibjb,k) = number of histogram bins encompassing
Ci                        (demin,demax) for a given ibjb pair and k
Ci         :Not used when job=0
Ci   jhw   :jhw(ibjb,kx) = pointer to whw
Ci         :Not used when job=0
Ci   q     :(mtet mode) q-vector separating occ, unocc states
Ci         :True q = 2*pi*q(1:3)/alat
Ci         :Not used now
Ci   iq    :(mtet mode) used to locate q vector in 1st BZ
Ci         :Not used now
Ci   isp1  :(mtet mode) spin index for occ   state (evals from disk)
Ci         :Not used now
Ci   isp2  :(mtet mode) spin index for unocc state (evals from disk)
Ci         :Not used now
Ci   nqibz :(mtet mode) # kpoints read from disk?
Ci         :Not used now
Co Outputs
Co   demin :minimum excitation energy for ib(occupied) jb(unoccupied) pair
Co         :Output when job=0 (see Bugs below)
Co         :Otherwise, not used
Co   demax :maximum excitation energy for ib(occupied) jb(unoccupied) pair
Co         :Output when job=0 (see Bugs below)
Co         :Otherwise, not used
Co   whw   :whw(i:j) histogram weights in bins i:j for given ib,jb,kx
Co         :i = jhw(ibjb,kx)
Co         :j = jhw(ibjb,kx)+nhw(ibjb),kx)-1
Co         :ibjb = ibjb(ib,jb,kx)
Co         :Generated when job=1
Cl Local variables
Cl         :
Cr Remarks
Cr  *tetwtq calculates by the tetrahedron method this integral:
Cr     whw(i,jb,ib,ik) =
Cr        \int_om_i^om_i+1 d\omg \times
Cr        \int d^3k f(e(k))(1-f(e(q+k)))\delta(omg-e2_jb(q+k)+e1_ib(k))
Cr   on a mesh of k-points ik for a sequence of histograms i=1,2,3,...
Cr   Here:
Cr     om_i   = lower energy bound of histogram, frhis(i)
Cr     om_i+1 = upper energy bound of histogram, frhis(i+1)
Cr     f(E)   = Fermi distribution function for T=0 (efermi is input)
Cr     e2_jb  = jbth state eband2 (unoccupied)
Cr     e2_ib  = ibth state eband1 (occupied)
Cr
Cr   Note: whw is stored in a compressed format; see Outputs
Cr   The total joint DOS in the interval (om_i,om_i+1) is
Cr     JNOS(i) = sum_(ib,jb,ik) whw(i,ib,jb,ik)
Cr
Cr   JNOS is not calculated here, but in tetwtt after call to tetwtq.
Cr
Cr   This code is based on: J.Rath&A.J.Freeman PRB11, 2109 (1975)
Cr   lindtet6 calculates the weights for one tetrahedron.
Cr
Cr  *Note that \delta(omg-e2_jb(q+k)+e1_ib(k)) is only the imaginary
Cr   part of the energy denominator 1(omg-e2_jb(q+k)+e1_ib(k) + idelta)
Cr   entering into Im x0.
Cr   Thus whw are real and corresponds to weights the imaginary part only.
Cr   The real part must be obtained from the imaginary part by a Hilbert
Cr   transformation. To calculate the full (complex) weight, use tetwtz.f.
Cr
Cr *Kotani's notes from tetwt5.F
Cr   The numbering of band index is not unique due to degeneracy.
Cr   It affects how to choose tetrahedron.  As a result, it affect on
Cr   the integration weight, and it might break the crystal symmetry.
Cr   So I add symmetrization at the last in this routine so as to recover the symmetry.
Cr
Cr *mtet mode: doubles the k-mesh in the integration of the delta-function
Cr              (mesh has 2x2x2 more tetrahedra than normal mode).
Cr   Normal  mode requires mtet=(/1,1,1/)
Cr   Doubled mode requires mtet=(/2,2,2/)
Cr   mtet mode is disabled at present
Cb Bugs
Cb   Input frhis is in Hartree units, while other input energies
Cb   (efermi, esciss, eband1, eband2, ecore) are in atomic Rydberg units.
Cb   demin, demax (internally generated) are in Hartree units.
Cu Updates
Cu   07 Sep 09  Adapted from Kotani's tetwt5.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer(4) :: job,nband,nqbz,nctot,ntetf,nqbzw,npm,ncc,
     .  idtetf(0:3,ntetf),ib1bz(nqbzw),nwgt(nqbz,npm)
      integer(4):: nwgtx,nhwtot,nwhis
      real(8) :: fac,efermi,esciss,qlat(3,3),qbzw(3,nqbzw),
     .  eband1(nband,nqbz),eband2(nband,nqbz),ecore(nctot)
      real(4) :: demax(nband+nctot,nband+ncc,nqbz,npm),
     .           demin(nband+nctot,nband+ncc,nqbz,npm)
      logical :: iwgt(nband+nctot,nband+ncc,nqbz,npm)
      real(8) :: whw(nhwtot),frhis(nwhis+1)
      integer(4):: ihw(nwgtx,nqbz,npm), ! omega pointer
     .             nhw(nwgtx,nqbz,npm), ! number of bins
     .             jhw(nwgtx,nqbz,npm), ! index to whw
     .             ibjb(nctot+nband,nband+ncc,nqbz,npm)
C     mtet mode
      integer(4)::isp1,isp2,iq,nqibz
      real(8) :: q(3)
C ... Local parameters
      integer(4):: jpm,ibxmx,jbxmx,jbx,nrankc1,nrankc2,nnn1,nnn2,
     .  nnni,nnnj,itet,ib,jb,iprint,i,j,ik,ibx,nbnc,nb1,ixx,ihis,
     .  ikx,ibib,ioff,isum,ini,ied,jini,iini,nnn,kk(0:3),kkm(0:3),
     .  nbandmx_kq,nbandmx_k,mxbnds,jobl
      integer(4):: noccx_k,noccx_kq,irnk1,nrank1,irnk2,nrank2,kx,
     .  im,ini1(nband+nctot),ied1(nband+nctot),ixi1,ixi2,ixe1,ixe2,
     .  ini2(nband+nctot),ied2(nband+nctot)
      real(8) :: det33,kvec(3,0:3),x(0:3),voltot,volt,ddw,wtthis(nwhis)
      real(8),target :: ek(nband+nctot,0:3),ekq(nband+nctot,0:3)
C     integer(4)::  noccx1
C     complex(8)::  wmean,a1,a2,wgt2
      logical :: ipr=.false.
      real(4) :: demax_,demaxx,demin_,deminn
      real(8),pointer:: eocc(:),eunocc(:)
      real(8) :: wtthis2(nwhis,0:3),piofvoltot,kkv(3,0:3)
      logical :: chkwrt,wxx
      real(8),parameter:: pi=3.1415926535897932d0

C ... Needed to synchronize with mtet mode
      logical ::mtett
      integer(4):: mtet(3),nmtet,nqbzwm,nqbzm
      integer(4),allocatable:: idtetfm(:,:,:),ib1bzm(:)
      real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
      real(8),allocatable:: ekzz1(:,:),ekzz2(:,:),wtet(:,:,:)
C ... for mtet mode
C      integer(4):: ntetfm,ifeig,ifmtet,nqnumm,kqxxm,nsp,nband_x,kp,im,
C     .  isp1,isp2,iq,nqibz,iqx,ispx,iqindx
C      integer(4),allocatable:: idtetfm(:,:,:),ib1bzm(:)
C      real(8),allocatable:: eigtet(:,:,:),ginv(3,3)
CC     real(8)::qdummy(3)
C      real(8) :: rk(3,nqbz)
C      complex(8),allocatable:: wgt1(:,:)

C ... External calls
      external addsciss,chkdgn,lindtet6

C --- Setup ---
      chkwrt = iprint()>=150
      voltot = abs(det33(qlat))
      piofvoltot = pi/4d0/voltot
      jobl = mod(job,10)

C ... For divided tetrahedron method
      mtet = (/1,1,1/)
      if (sum(abs(mtet))/=3) then
        mtett = .true.
        call info2(30,0,0,' TETWTQ:  divided-tetrahedron scheme: '//
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
C        if(nband_x /=nband ) then
C          stop 'tetwt5: nband_x /=nband'
C        endif
C        allocate( eigtet(nband,nqnumm,nsp) )
C        do iqx=1,nqnumm
C          do ispx=1,nsp
C            read (ifeig) eigtet(1:nband,iqx,ispx)
C            if (iprint()>150) write(6,"('iq=',i3,'  eig(1:5)=',5f10.4)")iqx,eigtet(1:5,iqx,ispx)
C          enddo
C        enddo
C        close(ifeig)
      else
        nmtet = 1
        nqbzwm= nqbzw
        nqbzm = nqbz
        allocate(wtet(0:3,nmtet,ntetf),qbzm(3,nqbzm),qbzwm(3,nqbzwm),
     .    idtetfm(0:3,nmtet,ntetf),ib1bzm(nqbzwm))
        idtetfm(:,1,:)=idtetf
        qbzwm = qbzw
        ib1bzm = ib1bz
        wtet = 0.25d0
      endif

c-----------------------------------------------------------------------
C     allocate( wgt1(nband+nctot,nband) )
      nbnc = nband+nctot
      nb1  = nband+1
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
C        do kx = 1, nqbzm
C          ekzz1( 1:nband, kx)              = eigtet(1:nband,kx,isp1)
C          ekzz1( nband+1: nband+nctot, kx) = ecore(1:nctot)
C          ekzz2( nband+1: nband+nctot, kx) = ecore(1:nctot)
C
Cccccccccccccccccccccc
Cc          do ib=1,nband
Cc          if(eigtet(ib,kx,isp1)>20) then
Cc            write(6,"('eigtet=',3i5,f13.5)")
Cc     .      ib,kx,isp1,eigtet(ib,kx,isp1)
Cc          endif
Cc          enddo
Cccccccccccccccccccccc
C
C        enddo
C
C        if(iq<=nqibz) then
C          do kx = 1, nqbzm
C!          q(1:3)+qbzm(1:3,kx) ----> kqxxm
Cc      print *,' xxx3 ginv=',ginv
Cc      print *,' q+qbzm= ', q(1:3) + qbzm(1:3,kx)
Cc          call fbz2 ( q(1:3) + qbzm(1:3,kx),
Cc     i      ginv,index_qbzm,n_index_qbzm,qbzm,nqbzm,
Cc     o      qdummy,kqxxm)
C            kqxxm = iqindx(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm)
Cc          write(6,*) 'kqxxm2 kqxxm=',kqxxm2,kqxxm
Cc          if(kqxxm/=kqxxm2) then
Cc             print *,' q=',q
Cc             stop 'kqxxm/=kqxxm2'
Cc          endif
C            ekzz2(1:nband, kx) = eigtet(1:nband,kqxxm,isp2)
C          enddo
C        else
C          do kx = 1, nqbzm
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

C ... Scissors operator
      if (esciss /= 0d0) then
        call addsciss(esciss,efermi,(nband+nctot)*nqbzm,ekzz1)
        call addsciss(esciss,efermi,(nband+nctot)*nqbzm,ekzz2)
      endif

C      call info5(10,0,0,' tetwtq: nqbz=%i vol=%;6d qlat=%9:1;6d '//
C     .  ' Ef=%;6d nband=%i',nqbz,voltot,qlat,efermi,nband)

C ... Sanity check: tetrahedra volumes must add to cell volume
      volt = 0d0
      do  itet = 1, ntetf
        do  im = 1, nmtet
          kvec(1:3,0:3) = qbzwm(1:3,idtetfm(0:3,im,itet) )
          do  i = 1,3
            kvec(1:3,i) = kvec(1:3,i) - kvec(1:3,0)
          enddo
          volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
        enddo
      enddo
      if (abs(volt-voltot)>1d-10) then
        call rx('tetwtq: sum of tetrahedra volumes does not '//
     .    'match unit cell volume')
      endif

      if (jobl==0) then
        demin =  1d10
        demax = -1d10
        iwgt = .false.
      endif

C --- Large loop over tetrahedra ---
      do  itet = 1, ntetf
        kk (0:3)      = ib1bz( idtetf(0:3,itet) )     ! Indices to 4 corners
        kkv(1:3, 0:3) = qbzw (1:3, idtetf(0:3,itet) ) ! k-points at corners
C   ... Loop over multitet
        do  im = 1, nmtet
C         kq (0:3)      = kqxx (kk(0:3)) !  k+q in micro-tet
          kkm (0:3)     = ib1bzm( idtetfm(0:3,im,itet) ) !  k   in micro-tet
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
          do i = 1, nband
            if (maxval(ekq(i, 0:3)) >= 9999) then
              nbandmx_kq = i-1
              exit
            endif
          enddo

C     ... If no time reversal, loop twice, exchanging i,j second pass
          do  jpm = 1, npm      !feb2006
            if (jpm==1) then
              ibxmx = noccx_k + nctot
              jbxmx = nbandmx_kq !nband  Apr2009takao
            else
              ibxmx = nbandmx_k  !nband  Apr2009takao
              jbxmx = noccx_kq + nctot
            endif
C       --- Loop over all potential (occ, unocc) pairs ---
            do  ibx  = 1, ibxmx  !noccx_k + nctot  !   occupied
            do  jbx  = 1, jbxmx  !nband            ! unoccupied
C         ... Get true ib,jb for this jmp
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

              if (jpm==1) then
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

C           ... Setup mode (jobl=0)
                if (jobl==0) then
                  iwgt(ib,jb,kk(0:3),jpm)= .true.
                  x(0:3) = .5d0*(eocc-eunocc) ! + omg !Denominator. unit in Hartree.
                  demax_ =  maxval(-x(0:3)) ! in Hartree
                  demin_ =  minval(-x(0:3)) ! in Hartree
                  do  ixx = 0, 3
                    demax(ib,jb,kk(ixx),jpm) =
     .                max(demax_,demax(ib,jb,kk(ixx),jpm))
                    demin(ib,jb,kk(ixx),jpm) =
     .                min(demin_,demin(ib,jb,kk(ixx),jpm))
                  enddo
C             ... Integration mode (jobl=1)
                else
                  x(0:3) = .5d0*(eocc-eunocc) ! + omg !Denominator. unit in Hartree.
                  wtthis2 = 0d0
                  if (chkwrt) then
                    write(6,"('### Goto lindtet6: itet ib jb Ef='
     .                ,i8,2i5,d12.4,' ###')") itet,ib,jb,efermi
                    write(6,"('  eocc  - Ef= ',4f10.3)") eocc  -efermi
                    write(6,"('  eunocc- Ef= ',4f10.3)") eunocc-efermi
                    write(6,"('  -x(a.u.)= ',4f10.3)") -x(0:3)
                  endif
! kvec(1:3, 0:3), ea(0:3), x(0:3) wgt(ib,jb,kk(0:3)) = wgt(ib,jb,kk(0:3)) + wttx/4d0/voltot
                  call lindtet6(kkv,kvec,eocc,eunocc,x,efermi,
     i              frhis,nwhis,linme,
     o              wtthis2)
                  if (chkwrt) then
                    write(6,"('  sum nwhis wtthis2=',i5,d13.5)")
     .                nwhis,sum(wtthis2)
                  endif

C               ... Debugging printouts
                  if (chkwrt) then
                    write(6,"(
     .                '  === ihis [range(a.u.)] wtthis2(0:3) ===')")
                    do  ihis = 1, nwhis
                      ddw = frhis(ihis+1)-frhis(ihis)
                      if (chkwrt.and.sum(abs(wtthis2(ihis,:)))/=0d0)
     .                  then
                        write(6,"('ttt',3i4,f10.5,4d12.4)")
     .                    itet,ib,jb,(frhis(ihis)+frhis(ihis+1))/2d0,
     .                    wtthis2(ihis,:)/ddw
                      endif
                      do ikx = 0, 3
                        ibib = ibjb(ib,jb,kk(ikx),jpm)
                        if (ihis< ihw(ibib,kk(ikx),jpm).or.
     .                    ihis> ihw(ibib,kk(ikx),jpm)- 1 +
     .                    nhw(ibib,kk(ikx),jpm)) then
                          if (sum(abs(wtthis2(ihis,0:3)))/=0d0)
     .                      stop 'tetwt5:wtthis2/=0 for out of range'
                        endif
                      enddo
                    enddo
                  endif

                  do  ikx = 0, 3
                    ibib = ibjb(ib,jb,kk(ikx),jpm)
                    jini = jhw(ibib,kk(ikx),jpm)
                    iini = ihw(ibib,kk(ikx),jpm)
                    nnn =  nhw(ibib,kk(ikx),jpm)

                    if (linme) then
                      whw(jini:jini+nnn-1) = whw(jini:jini+nnn-1) +
     .                  wtthis2(iini:iini+nnn-1,ikx) * piofvoltot
                    else
                      whw(jini:jini+nnn-1) = whw(jini:jini+nnn-1) +
     .                  wtthis2(iini:iini+nnn-1,0) * piofvoltot*4*
     .                  wtet(ikx,im,itet) ! piofvoltot= pi/voltot/4
                    endif
                  enddo

                endif
              endif

            enddo
          enddo
          enddo

C1100   continue
        enddo
      enddo
      deallocate(idtetfm, qbzwm,ib1bzm, qbzm)

c      print *, 'sumcheck wgt=', sum(wgt)

C --- Symmetrization of wgt and whw ---
C     Average contributions from degenerate levels.
      do  kx = 1, nqbz
C       Group occ,unocc bands into families of degenerate levels
        call chkdgn(ekxx1(:,kx),nband,nrank1,ini1,ied1,0,ipr)
        call chkdgn(ekxx2(:,kx),nband,nrank2,ini2,ied2,0,ipr)
C       Ditto for the core levels
        nrankc1 = 0
        if (nctot/=0) then
          call chkdgn(ecore,nctot,
     .      nrankc1,ini1(nrank1+1),ied1(nrank1+1),nband,ipr)
        endif
        nrankc2 = 0
        if (nctot/=0) then
          call chkdgn(ecore,nctot,
     .    nrankc2,ini2(nrank2+1),ied2(nrank2+1),nband,ipr)
        endif
        do  jpm = 1, npm
          if (jpm==1) then
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

C         ... Setup mode (jobl=0) : make demin, demax
              if (jobl==0) then
C               wmean = sum( wgt(ixi1:ixe1,ixi2:ixe2, kx))
C    .            / ((ixe1-ixi1+1)*(ixe2-ixi2+1))
C                wxx = .false.
C                if (count(iwgt(ixi1:ixe1,ixi2:ixe2,kx,jpm))>0)
C     .            wxx = .true.
                wxx = count(iwgt(ixi1:ixe1,ixi2:ixe2,kx,jpm)) > 0
                iwgt(ixi1:ixe1,ixi2:ixe2,kx,jpm) = wxx
                demaxx = maxval(demax(ixi1:ixe1,ixi2:ixe2,kx,jpm))
                demax(ixi1:ixe1,ixi2:ixe2,kx,jpm) = demaxx
                deminn = minval(demin(ixi1:ixe1,ixi2:ixe2,kx,jpm))
                demin(ixi1:ixe1,ixi2:ixe2,kx,jpm) = deminn

C         ... Integration mode (jobl=1): symmetrize
              else
C               Sanity check: every state in this group of (ib,jb) pairs
C               should be included, or every state left out
                isum = 0
                do  ib = ixi1, ixe1
                  do  jb = ixi2, ixe2
                    if (ibjb(ib,jb,kx,jpm)/=0) isum = isum+1
                  enddo
                enddo
                if (isum/=0 .and. isum /= (ixe1-ixi1+1)*(ixe2-ixi2+1))
     .            then
                  call rx('bug in tetwtq: inconsistent isum')
                endif
                if (isum==0) cycle

C               wtthis = whw averaged over (ib,jb) pairs in the group
                wtthis = 0d0
                do  ib = ixi1, ixe1
                  do  jb = ixi2, ixe2
                    ibib = ibjb(ib,jb,kx,jpm)
                    ini = ihw(ibib,kx,jpm)
                    ied = ihw(ibib,kx,jpm)+ nhw(ibib,kx,jpm)-1
                    ioff= jhw(ibib,kx,jpm)
                    wtthis(ini:ied) = wtthis(ini:ied) +
     .                whw(ioff:ioff+nhw(ibib,kx,jpm)-1)
                  enddo
                enddo
                wtthis = wtthis/((ixe1-ixi1+1)*(ixe2-ixi2+1))
C               Distribute average into whw
                do  ib = ixi1, ixe1
                  do  jb = ixi2, ixe2
                    ibib = ibjb(ib,jb,kx,jpm)
                    ini = ihw(ibib,kx,jpm)
                    ied = ihw(ibib,kx,jpm)+ nhw(ibib,kx,jpm)-1
                    ioff= jhw(ibib,kx,jpm)
                    whw(ioff:ioff+nhw(ibib,kx,jpm)-1) = wtthis(ini:ied)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo ! end of jpm loop
      enddo ! end of kx loop

C --- Make nwgt = number of (ib,jb) pairs for a given k ---
      if (jobl==0) then
        do  jpm = 1, npm
          do  ik = 1, nqbz
            nwgt(ik,jpm) = 0
            if (jpm==1) then
              nnni = nctot
              nnnj = 0
            else
              nnni = 0
              nnnj = nctot
            endif
            do  i = 1, nband+nnni
              do  j = 1, nband+nnnj
                if (iwgt(i,j,ik,jpm)) then
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
        call info2(20,0,0,' TETWTQ:  total number of (occ,unocc) '//
     .  'pairs = %i, %i or fewer per qp',
     .    sum(nwgt),maxval(nwgt(1:nqbz,1:npm)))
        if (sum(nwgt) /= count(iwgt))
     .    call rx(' bug in tetwt5: pairs improperly counted')
      endif

C --- Cleanup ---
C     deallocate(wgt1)
C     deallocate(idtetfm, qbzwm,ib1bzm, qbzm)
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

      subroutine hisrange(frhis,nwhis,demin,demax,ihw,nhw)
C- Find the histograms encompassing an energy range
C ----------------------------------------------------------------------
Ci Inputs
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci   nwhis :number of energies in histogram
Ci   demin :minimum excitation energy for ib(occupied) jb(unoccupied) pair
Ci   demax :maximum excitation energy for ib(occupied) jb(unoccupied) pair
Co Outputs
Co   ihw   :index to first bin encompassing demin
Co   nhw   :number of histograms needed encompass (demin,demax)
Cl Local variables
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer(4):: nwhis,ihw,nhw
      real(8):: frhis(nwhis+1)
      real(4):: demin,demax
C ... Local parameters
      integer(4):: ihis

      ihw = -9999
      nhw = -9999
C ... Find last bin whose lower bound is below demin
      do  ihis = 1, nwhis
        if (demin < frhis(ihis+1)) then
          ihw = ihis
          exit
        endif
      enddo
      if (ihw==-9999) then
        call rx2(' HISRANGE: energy %d exceeds upper bound '
     .    //'in histogram (%;6d)',dble(demin),frhis(nwhis+1))
      endif

C ... Find last bin whose lower bound is below demax
      do ihis = ihw, nwhis
        if (demax < frhis(ihis+1)) then
          nhw = ihis - ihw +1
          exit
        endif
      enddo
      if (nhw==-9999) then
        call rx2(' HISRANGE: energy %d exceeds upper bound '
     .    //'in histogram (%;6d)',dble(demax),frhis(nwhis+1))

      endif
      end

C      integer function mxbnds(mode,eb,nde,nband,nk,ef)
CC- Find the maximum number of states below Ef among a set of k-points
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0, find the largest  number of bands under Ef for any k
CCi         :1, find the smallest number of bands under Ef for any k
CCi   eb    :eigenvalues k-points 1..nk
CCi   nde   :leading dimension of eb
CCi   nband :maximum number of states to search
CCi   nk    :number of k-points
CCi   ef    :Fermi level
CCo Outputs
CCo   mxbnds: highest number of occupied states
CCl Local variables
CCr Remarks
CCr   Adapted from noccx1 in tetwt4.F
CCu Updates
CC ----------------------------------------------------------------------
CC      implicit none
CC ... Passed parameters
C      integer mode,nde,nband,nk
C      double precision eb(nde,nk),ef
CC ... Local parameters
C      integer nmax,k,it
C
C      nmax = 0
C      if (mode == 1) nmax = nband
C      do  k = 1, nk
C        do  it = 1, nband
C          if (eb(it,k) > ef) goto 1111
C        end do
C        it = nband
C 1111   continue
C        if (mode == 0) then
C          if (it > nmax) nmax = it
C        else
C          if (it < nmax) nmax = it
C        endif
CC       print *, it,nmax,k
C      end do
C      mxbnds = nmax
C
C      end

      subroutine lindtet6(kkv,kvec,ea,eb,x,efermi,frhis,nwhis,linme,
     .  wtthis)
C- Imaginary part of \int dk1 dk2 dk3 f(ea)(1-f(eb))/(\omega+x(k))
C ----------------------------------------------------------------------
Ci Inputs
Ci   kkv
Ci   kvec
Ci   ea    :energy at four corners of tetrahedron at k
Ci   eb    :energy at four corners of tetrahedron at k+q
Ci   x     :energy entering into delta-function integral
Ci   efermi:Fermi energy
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci   nwhis :number of energies in histogram
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Co Outputs
Co   wtthis: Weights are added into wtthis at each of 4 corners of tetrahedron
Cl Local variables
Cl         :
Cr Remarks
Cr  Accumulates Im \int dk1 dk2 dk3 f(ea)(1-f(eb))/(\omg+x)
Cr  by a histogram for one tetrahedron.
Cr  wtthis(ihis) += \int_lower(ihis)^upper(ihis) d\omg \int d^3k
Cr                   f(ea(k)) (1-f(eb(q+k))) \delta (omg + x(k))
Cr  where:
Cr    f(E) = Fermi distribution function for T=0
Cr  Tetrahedon is specified by values at 4 corners:
Cr    kvec(1:3, 1:4), ea(1:4), eb(1:4), x(1:4)
Cr  This code is based on J.Rath & A.J.Freeman PRB11, 2109 (1975).
Cr  Note:
Cr  \sum_ihis wthis(ihis) = volume of the microcell = \int d^3 k
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer nwhis
      real(8) :: kkv(3,4),kvec(1:3,1:4),x(1:4),ea(1:4),eb(1:4)
      real(8) :: frhis(nwhis+1),wtthis(nwhis,4),efermi
C ... Local parameters
C     integer(4) :: isig
      integer(4) :: ieaord(1:4),i,n,itmp,ix
      real(8)    :: kk(3,1:8),xx(1:8),ee(1:4),ebf(1:8)
!     real(8)    ::  ebfKx(1:4),Kx(1:3,4),xKx(4),
      logical:: chkwrt=.false.

      if (chkwrt) then
        print *, ' lindtet6: *****************************'
        write(6,"(' i=',i3,' x ea eb =',3f10.5)")
     .  (i, x(i), ea(i)-efermi, eb(i)-efermi,i=4,1,-1)
      endif

C ... Order points E_4>E_3>E_2>E_1, suitable  Rath&Freeman.
      n = 4
c      isig = 1
      do i = 1,n
        ieaord(i) = i
      enddo
      do  ix = 2, n
        do  i = ix, 2, -1
          if ( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
c            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      ieaord(1:4) = ieaord(4:1:-1)

C     4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
      kk(1:3,1:4) = kvec(1:3,ieaord(1:4))
      ee (1:4)    = ea (ieaord(1:4)) - efermi
      xx (1:4)    = x  (ieaord(1:4))
      ebf(1:4)    = eb (ieaord(1:4)) - efermi
ccccccccccccc
c      write(6,"(' i iea=',2i4)") (i,ieaord(i),i=1,4)
c      write(6,"(' i=',i3,' xx ee ebf =',3f10.5,'  kk=',3f10.5)")
c     .     (i,xx(i),ee(i),ebf(i),kk(1:3,i),i=4,1,-1)
ccccccccccccc

C ... Figure 0
      if ( 0d0<=ee(4) ) then
        if (chkwrt) write(6,*) 'lindtet6: fig 000'
C       wttx = (0d0,0d0)
C ... Figure 1
      elseif ( ee(4) < 0d0 .and. 0d0<= ee(3) ) then
        if (chkwrt) write(6,*) 'lindtet6: fig 1 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1 !K1 is on the line k4---k2.
        call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 4,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3
C       k4,K1,K2,K3
        call inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/),frhis,nwhis,
     .    linme,wtthis)
C ... Figure 2
      elseif ( ee(3) < 0d0 .and. 0d0<= ee(2) ) then
        if (chkwrt) write(6,*) 'lindtet6: fig 2 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
        call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 3,1,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3
        call midk3(kk,ee,xx,ebf, 3,2,  kk(1,4+4),xx(4+4),ebf(4+4)) !K4
        if (chkwrt) write(6,*) 'lindtet6: fig 2 xxx2'
C       k4,k3,K1,K2
        call inttetra6(kkv,kk,xx,ebf,(/4,3,1+4,2+4/),frhis,nwhis,
     .    linme,wtthis)
        if (chkwrt) write(6,*) 'lindtet6: fig 2 xxx3'
C       k3,K2,K3,K1
        call inttetra6(kkv,kk,xx,ebf,(/3,2+4,3+4,1+4/),frhis,nwhis,
     .    linme,wtthis)
        if (chkwrt) write(6,*) 'lindtet6: fig 2 xxx4'
C       k3,K1,K3,K4
        call inttetra6(kkv,kk,xx,ebf,(/3,1+4,3+4,4+4/),frhis,nwhis,
     .    linme,wtthis)
        if (chkwrt) write(6,*) 'lindtet6: fig 2 xxx5'

C ... Figure 3
      elseif ( ee(2)<0d0 .and. 0d0<=ee(1) ) then
        if (chkwrt) write(6,*) 'lindtet6: fig 3 xxx'
        call midk3(kk,ee,xx,ebf, 1,4,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
        call midk3(kk,ee,xx,ebf, 1,2,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 1,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3
C       k3,k4,K3,k2
        call inttetra6(kkv,kk,xx,ebf,(/3,4,3+4,2/),frhis,nwhis,
     .    linme,wtthis)
C       k4,K1,K2,K3
        call inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/),frhis,nwhis,
     .    linme,wtthis)
C       k4,k2,K2,K3
        call inttetra6(kkv,kk,xx,ebf,(/4,2,2+4,3+4/),frhis,nwhis,
     .    linme,wtthis)
      else
        if (chkwrt) write(6,*) 'lindtet6: fig 4 xxx'
C       k1,k2,k3,k4
        call inttetra6(kkv,kk,xx,ebf,(/1,2,3,4/),frhis,nwhis,
     .    linme,wtthis)
      endif
      if (chkwrt) write(6,*) 'end of lindtet6',wtthis
c      stop 'yyyyyyyyyyyyyyyyyyyyyy'
      end

      subroutine inttetra6(kkv,kk_,xx_,ebf,itetx,frhis,nwhis,linme,
     o  wtthis)
C- Calculate tetrahedron integral Rath&Freeman Eq.(16).
C ----------------------------------------------------------------------
Ci Inputs
Ci   kkv   :
Ci   kk_   :
Ci   xx_   :
Ci   ebf   :
Ci   itetx : the four corners of the tetrahedron
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci   nwhis :number of energies in histogram
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Co Outputs
Co   wtthis: Weights are added into wtthis at each of 4 corners of tetrahedron
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   07 Sep 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer(4):: nwhis,itetx(4)
      real(8) ::  kk_(3,1:8),xx_(1:8),kkv(3,4)
!     real(8) ::  Kx(3,1:4),xKx(1:4)
      real(8):: frhis(nwhis+1),wtthis(nwhis,4),ebf(1:8)
C ... Local parameters

C calculate tetrahedron integral Eq.(16).
C the four corners and denoted by itetx.
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx.
Cr wtthis is accumulated
C     integer(4):: isig
      integer(4):: ix,i,ieaord(4),n,itmp
      real(8) ::
     .  kk(3,1:8),xx(1:8),ebfin(1:8) ,ee(4)  ,kin(3,4), xin(4)
!     real(8) ::  Kx_(3,1:4),xKx_(1:4),ebfKx(1:4)
      logical :: chkwrt=.false.

      if (chkwrt) then
        print *,' inttetra6: ---------------------------'
c        print *,' ebf  =', ebf(1:4)
c        print *,' ebfKx=', ebfKx(1:4)
      endif

c ... kk xin ebfin
c      do i = 1,4
c        if ( itetx(i) < 10 ) then   ! kk (k1-k4 in Paper)
c          ix = itetx(i)
c          xin(    i) = xx_(    ix)
c          ebfin(i)   = ebf(ix)
c          kin(1:3,i) = kk_(1:3,ix)
c        else                       ! Kx (K1-K4 in Paper)
c          ix = itetx(i)/10
c          xin(    i) = xKx_(    ix)
c          ebfin(i)   = ebfKx(ix)
c          kin(1:3,i) = Kx_ (1:3,ix)
c        endif
c      enddo

c ... kk xin ebfin
      xin(    1:4) = xx_(    itetx(1:4))
      kin(1:3,1:4) = kk_(1:3,itetx(1:4))
      ebfin(  1:4) = ebf(    itetx(1:4))
c
      if (chkwrt) then
        write(6,"(' i=',i3,' xin ebfin =',2f10.5,'  kin=',3f10.5)")
     .  (i,xin(i),ebfin(i),kin(1:3,i),i=4,1,-1)
      endif

c  ... ee decending order ee(4)>ee(3)>ee(2)>ee(1)
      n = 4
C     isig = 1
      do i = 1,n
        ieaord(i) = i
      enddo
      do  ix =  2,n
        do  i = ix,2,-1
          if (ebfin(ieaord(i-1)) > ebfin(ieaord(i))) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
C           isig= -isig
            cycle
          endif
          exit
        enddo
      enddo

      kk(1:3,1:4) = kin  (1:3,ieaord(1:4))   ! 4 corners denoted by kkv(:,1:4), ee(1:4), and xx(1:4)
      ee(    1:4) = ebfin(    ieaord(1:4))
      xx(    1:4) = xin  (    ieaord(1:4))
c
      if (chkwrt) write(6,"(' i=',i3,' xx ee =',2f10.5,'  kk=',3f10.5)")
     .  (i,xx(i),ee(i),kk(1:3,i),i=4,1,-1)
c
      if ( 0d0>=ee(4) ) then
!        wttx = (0d0,0d0)
      elseif ( ee(4) > 0d0 .and. 0d0 >= ee(3) ) then   !!! Fig 1.
        if (chkwrt) write(6,*) 'inttetra5: fig 1'
        call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 4,3,  kk(1,3+4),xx(3+4)) !K3
C       k4,K1,K2,K3
        call inttetrac6(kkv,kk,xx,(/4,1+4,2+4,3+4/),frhis,nwhis,linme,
     o         wtthis)
      elseif ( ee(3) > 0d0 .and. 0d0>= ee(2) ) then   !!! Fig 2.
        if (chkwrt) write(6,*) 'inttetra5: fig 2'
        call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1
        call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 3,1,  kk(1,3+4),xx(3+4)) !K3
        call midk(kk,ee,xx, 3,2,  kk(1,4+4),xx(4+4)) !K4
C       k4,k3,K1,K2
        call inttetrac6(kkv,kk,xx,(/4,3,1+4,2+4/),frhis,nwhis,linme,
     o         wtthis)
C       k3,K2,K3,K1
        call inttetrac6(kkv,kk,xx,(/3,2+4,3+4,1+4/),frhis,nwhis,linme,
     o         wtthis)
C       k3,K1,K3,K4
        call inttetrac6(kkv,kk,xx,(/3,1+4,3+4,4+4/),frhis,nwhis,linme,
     o         wtthis)
      elseif ( ee(2) > 0d0 .and. 0d0>= ee(1) ) then   !!! Fig 3.
        if (chkwrt) write(6,*) 'inttetra5: fig 3'
        call midk(kk,ee,xx, 1,4,  kk(1,1+4),xx(1+4)) !K1
        call midk(kk,ee,xx, 1,2,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 1,3,  kk(1,3+4),xx(3+4)) !K3
C       k3,k4,K3,k2
        call inttetrac6(kkv,kk,xx,(/3,4,3+4,2/),frhis,nwhis,linme,
     o         wtthis)
C       k4,K1,K2,K3
        call inttetrac6(kkv,kk,xx,(/4,1+4,2+4,3+4/),frhis,nwhis,linme,
     o         wtthis)
C       k4,k2,K2,K3
        call inttetrac6(kkv,kk,xx,(/4,2,2+4,3+4/),frhis,nwhis,linme,
     o         wtthis)
      else
C       k1,k2,k3,k4
        if (chkwrt) write(6,*) 'inttetra5: fig 4'
        call inttetrac6(kkv,kk,xx,(/1,2,3,4/),frhis,nwhis,linme,
     o    wtthis)
      endif
      end

      subroutine inttetrac6(kkv,kk,xx,itetx,frhis,nwhis,linme,
     o         wtthis)
C- Kernel called by inttetra6
C ----------------------------------------------------------------------
Ci Inputs
Ci   kkv   :
Ci   kk    :
Ci   xx    :
Ci   itetx : the four corners of the tetrahedron
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci   nwhis :number of energies in histogram
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Co Outputs
Co   wtthis: Weights are added into wtthis at each of 4 corners of tetrahedron
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   07 Sep 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer(4):: nwhis,itetx(4)
      real(8) ::  kk(3,1:8),xx(1:8),kkv(3,4)
!     real(8) ::  Kx(3,1:4),xKx(1:4)
      real(8):: frhis(nwhis+1),wtthis(nwhis,4)
C ... Local parameters
C
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx(1:4).
Cr wtthis is accumulating.
      integer(4):: i
      real(8) ::  am(3,3),bm(3,3),kin(3,4),kkvkin(4,4),xin(4),det33
      real(8):: voltet,kkv4bm(3),xtmp
      logical:: chkwrt=.false.

c
      if (chkwrt) print *, ' inttetrac6: === '

      xin(    1:4) = xx(    itetx(1:4))
      kin(1:3,1:4) = kk(1:3,itetx(1:4))
      do  i = 1,3
        am(1:3,i) = kin(1:3,i) - kin(1:3,4)
      enddo
      voltet = abs(det33(am)/6d0) ! \omega (volume of tetrahedra) = abs(det33(am)/6d0) See Eq. (17).

c ... kkvkin: kin is decomplosed into kkv
C     e.g. kin (:,1) =  \sum_i kkv(:,i) * kkvkin(i,1)
      if (linme) then
        do  i = 1, 3
          am(1:3,i) = kkv(1:3,i) - kkv(1:3,4)
        enddo
        call mkqlat(am,bm,xtmp)
        kkv4bm(1) = sum( kkv(:,4)*bm(:,1))
        kkv4bm(2) = sum( kkv(:,4)*bm(:,2))
        kkv4bm(3) = sum( kkv(:,4)*bm(:,3))
        do i = 1, 4
          kkvkin(1,i) = sum( kin(:,i)*bm(:,1) ) - kkv4bm(1)
          kkvkin(2,i) = sum( kin(:,i)*bm(:,2) ) - kkv4bm(2)
          kkvkin(3,i) = sum( kin(:,i)*bm(:,3) ) - kkv4bm(3)
          kkvkin(4,i) = 1d0 - sum(kkvkin(1:3,i))
        enddo
      endif

      call intttvc6(kkvkin,xin,voltet,linme,frhis,nwhis,wtthis)
      end

      subroutine intttvc6(kkvkin,v,voltet,linme,frhis,nwhis,
     o         wtthis)
C- Weights for each bin in the histogram.
C ----------------------------------------------------------------------
Ci Inputs
Ci   kkvkin:
Ci   v     :energy entering into delta-function: delta(omega+v)
Ci   voltet:
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci   nwhis :number of energies in histogram
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Co Outputs
Co   wtthis: Weights are added into wtthis at each of 4 corners of tetrahedron
Cl Local variables
Cl   ieaord:orders 4 corners of tetrahedron by increasing energy
Cl   ww    :v, ordered by increasing energy
Cr Remarks
Cr   The i-th bin has fermi function weights [frhis(i) frhis(i+1)].
Cr    wtthis(ihis) += \int_{frhis(ihis)^{frhis(ihis+1)} d \omega *
Cr                    \int d^3k \delta(\omega + v(k) )
Cr  Total sum of the histogram weight is equal to voltet.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer nwhis
      real(8):: kkvkin(4,4),wtthis(nwhis,4)
      double precision v(4),voltet,frhis(nwhis+1)
C ... Local parameters
      logical   :: chkwrt=.false.
      integer(4):: i, inihis,iedhis,ichk=2313
      integer(4):: ieaord(1:4),ix,n,itmp,ihis
      real(8)   :: ww(4),norm,pvn,intega,integb
      real(8)   :: www(1:4)=.25d0,ec,wcg(4),wttt

      n=4
c     isig = 1
      do   i = 1, n
        ieaord(i) = i
      enddo
      do  ix = 2, n
        do  i = ix, 2, -1
          if (-v(ieaord(i-1)) > -v(ieaord(i))) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
c            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo

      ww(1:4) = -v( ieaord(1:4) )
C     Check that ww(1)<ww(2)<ww(3)<ww(4)
      if ((.not.(WW(1)<=WW(2))) .or. (.not.(WW(2)<=WW(3))) .or.
     .    (.not.(WW(3)<=WW(4))) ) then
        write(6,"(/,' --- intttvc6: wrong order WW=',4d14.6)") WW
        stop 'intttvc6: wrong order of WW'
      endif

      if (chkwrt) then
        write(ichk,"(/,' --- intttvc6: e=',4d23.16)") WW
      endif

C ?
      inihis =  -999
      iedhis =  -999
      ix = 1
      do  ihis  =  1, nwhis
        if (ix == 1 .and. WW(ix)<frhis(ihis+1)) then
          inihis = ihis
          ix = 4
        endif
        if (ix == 4 .and. WW(ix)<frhis(ihis+1)) then
          iedhis = ihis
          exit
        endif
      enddo

      if (iedhis == -999 .or. inihis == -999) then
        print *,' intttvc6: can not find inihis iedhis'
        stop    ' intttvc6: can not find inihis iedhis'
      endif

c      if (chkwrt) print *,' inihis iedhis=', inihis,iedhis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' ###########integtetn test #########'
c      ww(1:4)=(/1d0,1.5d0,9.5d0,10d0/)
c      do ix=1,101
c        wx = (ix-1)/10d0
c        call integtetn(WW, Wx, xxx)
c        write(61,"(i3,2d23.16)")ix,wx,xxx
c      enddo
c      stop 'test end integtetn---'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call integtetn(WW,WW(4),norm)
C     if (chkwrt) write(ichk,"(' norm voltet=',2d24.16)") norm,voltet
      pvn = voltet/norm
      if (chkwrt) write(ichk,"(' norm voltet pvn=',3d14.6)")
     .  norm,voltet,pvn

      intega = 0d0
      do  ihis = inihis, iedhis
        if ( frhis(ihis+1)>ww(4) ) then
          integb = norm
        else
          call integtetn(WW,frhis(ihis+1),integb)
        endif

        if (linme) then !-- Weight for each corner sum(wcg)=1d0
          ec = (frhis(ihis)+frhis(ihis+1))/2d0
          call mkwcg(WW, ec, wcg)
          do  ix = 1, 4
            www(ix)= sum(kkvkin(ix,ieaord(1:4))*wcg(1:4))
          enddo
c          write(6,"('sum(www)= ',d13.4,' www(1:4)=',4f10.4)") sum(www),www(1:4)
          wttt = pvn*(integb - intega)*4d0
          wtthis(ihis,:) = wtthis(ihis,:) + wttt * www(:)
        else
          wtthis(ihis,1) = wtthis(ihis,1) + pvn*(integb - intega)
        endif

        if (chkwrt) then
          write(ichk,"(' ihis [Init End] wtt=', i5,3f11.6)")
     .    ihis, frhis(ihis), frhis(ihis+1), pvn*(integb-intega)
        endif

        intega = integb
      enddo

c      if (chkwrt) then
c        dini = 0d0
c        do ihis=inihis,iedhis
c          if (ihis>1) dini=wtthis(ihis-1)
c          write(6,"(' ihis [Init  End]', i4,2f13.6,' wtt=',f13.6)")
c     .    ihis,dini,frhis(ihis),wtthis(ihis)
c          dini = frhis(ihis)
c        enddo
c      endif
      if (chkwrt) write(ichk,*) ' end of intttvc6'
      end

      subroutine integtetn(e,ee,integb)
C- Calculate primitive integral integb = 1/pi Imag[\int^ee dE' 1/(E' -e(k))]
C-                                     = \int^ee dE' S[E']
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :energies at 4 corners of tetrahedron
Ci   ee    :plane of constant energy
Co Outputs
Co   integb:integral S(ee)
Cl Local variables
Cl         :
Cr Remarks
Cr   S[E] : is area of the cross-section between the omega-constant plane and the tetrahedron.
Cr   This routine assumes e1<e2<e3<e4.
Cr   Normalization is not considered here.
Cr   See Rath&Freeman integration of imaginary part, Eq.(17)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8) ::  e(1:4),ee,integb
C ... Local parameters
      real(8) ::  a,b,V1,V2,V3,V4,D1,D2,D3,D4,e1,e2,e3,e4
C     integer(4)::i1,i2,i3

c-------------------------------------------------------------
c      if ((.not.(e(1)<=e(2))).or.(.not.(e(2)<=e(3)))
c     . .or.(.not.(e(3)<=e(4))) ) then
c        write(6,"(/,' --- integtetn wrong order of e=',4d14.6)") e
c        stop 'integtetn: wrong order of e'
c      endif

c      write(6,"(' --- integtetn e =',4d23.16)") e
c      write(6,"(' --- integtetn ee=',d23.16)") ee

      if ( ee<e(1) ) then    !jan2008 ee<=e(1) ---> ee<e(1)
        integb = 0d0
        return
      elseif (ee>e(4) ) then
        call rx(' integtetn: ee>e(4)')
      endif

C --- Case 1. poor numerical accuracy for cases.
      if (.false.) then
        e1 = e(1)-3d-6
        e2 = e(2)-2d-6
        e3 = e(3)-1d-6
        e4 = e(4)
        V1= ee - e1
        V2= ee - e2
        V3= ee - e3
        V4= ee - e4
        D2 = (V2-V4)*(V2-V3)*(V2-V1) !<0
        D3 = (V3-V4)*(V3-V2)*(V3-V1) !>0
        D1 = (V1-V4)*(V1-V3)*(V1-V2) !>0
        if ( e1<=ee ) integb =          V1**3/D1
        if ( e2<=ee ) integb = integb + V2**3/D2
        if ( e3<=ee ) integb = integb + V3**3/D3
        return
      endif

C --- Case 2 ---
      e1 = e(1)-3d-8
      e2 = e(2)-2d-8
      e3 = e(3)-1d-8
      e4 = e(4)
      V1= ee - e1
      V2= ee - e2
      V3= ee - e3
      V4= ee - e4
c      D1 = (V1-V4)*(V1-V3)*(V1-V2) !>0
c      D2 = (V2-V4)*(V2-V3)*(V2-V1) !<0
c      D3 = (V3-V4)*(V3-V2)*(V3-V1) !>0
c      D4 = (V4-V1)*(V4-V2)*(V4-V3) !<0
      if    ( e1<=ee .and. ee<e2 ) then
        D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
        integb =  V1**3/D1
      elseif ( e2<=ee .and. ee<e3 ) then
        D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
        D2 = (e4-e2)*(e3-e2)*(e1-e2) !<0
        a  =  V1/   D1**(1d0/3d0)
        b  =  V2/(-D2)**(1d0/3d0)
        integb = (a-b) * (a**2+a*b+b**2)
      elseif ( e3<=ee .and. ee<e4 ) then
        D4 = (e1-e4)*(e2-e4)*(e3-e4) !<0
        integb = 1d0 - V4**3/D4
      elseif ( ee==e4 ) then
        integb = 1d0
      endif

c-----------------------------------------------------------
c      D2 = (V2-V4)*(V2-V3)*(V2-V1) !>0
c      D3 = (V3-V4)*(V3-V2)*(V3-V1) !<0
c      integb = 0d0
c      if ( e(1)<=ee ) integb =          V1**3* D2*D3
cc      write(6,*) ' integb1=',integb
c      D1 = (V1-V4)*(V1-V3)*(V1-V2) !<0
c      if ( e(2)<=ee ) integb = integb + V2**3* D3*D1
cc      write(6,*) ' integb2=',integb
c      if ( e(3)<=ee ) integb = integb + V3**3* D1*D2
cc      write(6,*) ' integb3=',integb
c--
c      if (V4==0d0) then
c       write(6,*)
c       write(6,"(' V  =',4d18.10)") V1, V2, V3, V4
c       write(6,"(' int=',3d13.5,16x,d13.5)") V1**3/D1,V2**3/D2,V3**3/D3,
c     .             max(abs(V1**3/D1),abs(V2**3/D2),abs(V3**3/D3))
c         write(6,*) ' integb=',integb
c      endif
      end

      subroutine mkwcg(e,ee,wcg)
C- Calculate weight for each corner
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :e(1:4) at corners
Ci   ee    :energy
Co Outputs
Co   wcg   :weights
Cl Local variables
Cl         :
Cr Remarks
Cr   Normalization is sum(scg(1:4))=1d0
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8) ::  e(1:4), ee, wcg(1:4)
C ... Local parameters
      real(8) :: e1,e2,e3,e4,w14_1,w14_4,w12_1,w12_2,w13_1,w13_3,w23_2,
     .  w23_3,w24_2,w24_4,w34_3,w34_4

C This is for original case; changed Dec 2003.
c      wcg=.25d0
c      return

      e1 = e(1)-3d-8
      e2 = e(2)-2d-8
      e3 = e(3)-1d-8
      e4 = e(4)
      if    ( ee<=e(1) ) then
        wcg(1)  =1d0
        wcg(2:4)=0d0
      elseif ( e1<=ee .and. ee<e2 ) then
        call wab(e1,e2,ee,w12_1,w12_2)
        call wab(e1,e3,ee,w13_1,w13_3)
        call wab(e1,e4,ee,w14_1,w14_4)
        wcg(1) = (w12_1 +w13_1 + w14_1)/3d0
        wcg(2) =  w12_2/3d0
        wcg(3) =  w13_3/3d0
        wcg(4) =  w14_4/3d0
      elseif ( e2<=ee .and. ee<e3 ) then !Is this correct?
        call wab(e1,e3,ee,w13_1,w13_3)
        call wab(e1,e4,ee,w14_1,w14_4)
        call wab(e2,e3,ee,w23_2,w23_3)
        call wab(e2,e4,ee,w24_2,w24_4)
        wcg(1) = (w13_1 + w14_1)/4d0
        wcg(2) = (w23_2 + w24_2)/4d0
        wcg(3) = (w13_3 + w23_3)/4d0
        wcg(4) = (w14_4 + w24_4)/4d0
      elseif ( e3<=ee .and. ee<e4 ) then
        call wab(e1,e4,ee,w14_1,w14_4)
        call wab(e2,e4,ee,w24_2,w24_4)
        call wab(e3,e4,ee,w34_3,w34_4)
        wcg(1) =  w14_1/3d0
        wcg(2) =  w24_2/3d0
        wcg(3) =  w34_3/3d0
        wcg(4) = (w24_4 +w34_4 + w14_4)/3d0
      elseif ( ee> e(4) ) then
        wcg(4)  =1d0
        wcg(1:3)=0d0
      endif
      end

      subroutine wab(ea,eb,ee,wa,wb)
      implicit none
      real(8)::ea,eb,wa,wb,eet,ee
      eet= eb - ea
      wa= (eb-ee)/eet
      wb= (ee-ea)/eet
      end


      subroutine pair_index(jpm,iwgt,nqbz,nband,nctot,ncc,nwgtx,
     .  nib,njb,noccxv,nbnb)
C- Get indices ib,jb for (ib,jb) pairs with nonzero weight
C ----------------------------------------------------------------------
Ci Inputs
Ci   jpm   :1 for normal case, 2 when swap (ib,jb) (no time reversal)
Ci   iwgt  :T ib(occ) jb(unocc) pair makes nonzero contribution
Ci   nqbz  :Number of k-points in 1st BZ
Ci   nband :number of bands
Ci   nctot :number of cores
Ci   ncc   :number of cores for jb in (ib,jb) pair:
Ci         :Should be 0 if time-reversal, nctot otherwise
Ci   nwgtx :dimensions nib,njb.  Must be at least as large as the
Ci         :maximum number of ib(occupied) jb(unoccupied) pairs for any k
Co Outputs
Co   nib   :nib(ix,ik)= ib for ixth nonzero (ib,jb) pair
Co   njb   :njb(ix,ik)= jb for ixth nonzero (ib,jb) pair
Co   noccxv:?
Co   nbnb  :nbnb(ik) = number of nonzero (ib,jb) pairs for this k
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Sep 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer(4) :: jpm,nqbz,nband,nctot,ncc,nwgtx
      integer(4) :: nib(nwgtx,nqbz),njb(nwgtx,nqbz),noccxv,nbnb(nqbz)
      logical    :: iwgt(nband+nctot,nband+ncc,nqbz)
C ... Local parameters
      integer(4) :: ib,jb,ik,ix

      noccxv = 0
      do  ik  = 1, nqbz
        ix  = 0
        do  ib  = 1, nband + nctot
          do  jb  = 1, nband + ncc
            if (iwgt(ib,jb,ik)) then
              ix          = ix+1
              nib(ix, ik) = ib
              njb(ix, ik) = jb
              if (jpm==1 .and. ib<=nband .and. ib>noccxv) noccxv = ib
              if (jpm==2 .and. jb<=nband .and. jb>noccxv) noccxv = jb
            endif
          enddo
        enddo
        nbnb(ik) = ix
      enddo
      end

C      subroutine dinv33x(plat,qlat)
CC- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
CCr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
C      implicit none
C      double precision plat(3,3),qlat(3,3),det
C      call cross(plat(1,2),plat(1,3),qlat)
C      call cross(plat(1,3),plat,qlat(1,2))
C      call cross(plat,plat(1,2),qlat(1,3))
C      det  = sum( plat(1:3,1)*qlat(1:3,1) )
C      qlat = qlat/det
C      end

C      integer function noccx1(ekt,nk,nt,ef)
CC- Find the maximum number of states below Ef amnong a set of k-points
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ekt   :eigenvalues for all k-points and states
CCi   nk    :number of k-points
CCi   nt    :number of states, and leading dimension of ekt
CCi   ef    :Fermi level
CCo Outputs
CCo   noccx1: highest number of occupied states
CCl Local variables
CCr Remarks
CCr   Adapted from noccx1 in tetwt4.F
CCu Updates
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nk,nt
C      double precision ekt(nt,nk),ef
CC ... Local parameters
C      integer noccx,k,it
C
C      noccx = 0
C      do  k = 1, nk
C        do  it = 1, nt
C          if (ekt(it,k) > ef) goto 1111
C        end do
C 1111   if (it > noccx) noccx = it
C      end do
C      noccx1 = noccx
C
C      end
