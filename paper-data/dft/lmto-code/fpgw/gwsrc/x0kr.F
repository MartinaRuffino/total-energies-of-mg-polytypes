      subroutine x0kr(npm,ihw,nhw,jhw,whw,nhwtot,n1b,n2b,nbnbx,nbnb,q,isp1,isp2,qlat,rk,
     i  mdim,nlnmv,nlnmc,iclass,ppb,icore,ncore,nlmto,nqbz,noccx,noccxv,
     i  nl,nclass,natom,nnc,nlnmx,mdimx,nbloch,nw,iq,nmb,ngc,nband,
     i  nolfco,zzr,nmbas1,nmbas,ppovlz,ppovl,chipmzzr,eibzmode,nwgt,ebasis,
     o  rcxq)
C- Calculate the imaginary part of chi0 or transverse magnetic chi0 along the real axis.
C ----------------------------------------------------------------------
Ci Inputs
Ci  npm    : 1 if time-reversal symmetry, 2 if not
Ci  ...      See tetwt5 for parameters isw,nhw,jhw,whw,nhwtot,n1b,n2b,nbnbx,nbnb
Ci  ihw    : ihw(ibjb,k): index to first histogram bin which discretizes interval (demin,demax)
Ci         : for a given ibjb pair and k.
Ci         : ibjb is an index to an (occ,unocc) pair of states and
Ci         : (demin,demax) = energy window of histogram
Ci         :               = (min,max) excitation energy to include in polarizability
Ci  nhw    : nhw(ibjb,k) = number of histogram bins for a given ibjb pair and k
Ci  jhw    : jhw(ibjb,k) = pointer to whw
Ci  whw    : whw(i:j) histogram weights in bins i:j for given ibjb = (occ,unocc) pair and k
Co         : i = jhw(ibjb,kx)
Co         : j = jhw(ibjb,kx)+nhw(ibjb),kx)-1
Co         : ibjb = ibjb(ib,jb,kx)
Ci  nhwtot : dimensions whw
Ci  n1b    : n1b(ibjb,1) = band index to occupied state ib for (occ,unocc) pair ibjb
Ci         : n1b(ibjb,2) = band index to unoccupied state ib for (unocc,occ) pair ibjb
Ci         :               used only if npm=2
Ci  n2b    : n2b(ibjb,2) = band index to unoccupied state jb for (occ,unocc) pair ibjb
Ci         : n2b(ibjb,2) = band index to occupied state jb for (unocc,occ) pair ibjb
Ci         :               used only if npm=2
Ci  nbnbx  : Dimensions n1b and n2b
Ci  nbnb   : nbnb(k,jpm) = number of (occ,unocc) pairs which contribute to polarizability
Ci  q      : Make chi0(q)
Ci  isp1   : Spin index for wave functions at k
Ci  isp2   : Spin index for wave functions at k+q.  Different from isp1 for transverse chi
Ci  qlat   : Reciprocal lattice vectors
Ci         : NEW: qlat
Ci  rk     : k-points in the 1st BZ CHECK
Ci  mdim   : nrprdc(i) = number of product basis function for class i
Ci         : NEW: nrprdc
Ci  nlnmv  : For each class, number of partial waves for valence
Ci  nlnmc  : For each class, number of partial waves for core
Ci  iclass : the jth atom belongs to class iclass(j)
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'n') B(R,I)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Ci  icore  : index for core states
Ci  ncore  : no. core states in each class
Ci  naugm  : number of augmentation channels ... called ndima in sugw.f
Ci         : OLD: nlmto
Ci  nqbz   : number of k-points in the 1st BZ
Ci  noccxv : maximum number of occupied valence states at any k
Ci  noccx  : maximum number of occupied states = noccxv + number of core states
Ci  nl     : (global maximum l) + 1
Ci  nclass : number of inequivalent classes
Ci  nbas   : number of sites
Ci         : OLD natom
Ci  nnc    : dimensions icore
Ci  nlnmx  : maximum number of product basis functions; dimensions ppb
Ci  nrprdcx: dimensions ppb
Ci         : OLD: mdimx
Ci  nrprd  : number of product basis functions B.
Ci           OLD nbloch
Ci  ngprd  : number of G vectors in mixed basis
Ci         : OLD: ngc
Ci  nmbas1 : Number of mixed basis functions : rcxq = rcxq(nmbas1,nmbas,:,:)
Ci  nmbas  : Number of mixed basis functions : rcxq = rcxq(nmbas1,nmbas,:,:)
Ci  nmb    : total number of mixed product basis functions I,J in which <I|v|J> is calculated
Ci         : (formerly ngb or ngbb)
Ci  nw     : number of frequency points at which chi is calculated
Ci  iq     : index to current q
Ci  ndham  : dimensioning parameter: maximum number of bands at any k.
Ci         : OLD: nband
Ci  nolfco : T => calculate x0 without local fields
Ci  zzr    : <eiqr|B_I> ... used for magnetic susceptibility calculations (hvccfp0.m.F)
Ci  ppovlz : ppovlz(I,mu) = \sum_J O_IJ Zcousq(J, mu)
Ci         : I = mixed basis function index
Ci         : mu = index to eigenfunction of the bare coulomb interaction at this q
Ci  ppovl  :<I|J> = Overlap in product basis
Ci         : should be replaced with what is called here ppovl0 -- used only when eibzmode=T
Ci  chipmzzr:T => calculate transverse chipm
Ci  eibzmode:T => use symmetry to reduce the number of k in convolution G(k) * G(k+q)
Ci  nwgt   : nwgt(k) = number of k-points symmetry equivalent to k
Ci         : nwgt(k) = 0 => this k is subsumed under another, equivalent k point
Ci         : Not used unless eibzmode = T.
Ci  ebasis : If nonzero, make rcxq in E-Basis.
Ci         : This occurs anyway if nmbas1 or nmbas differs from nmb
Ci         : and it is not applicable if using spin-basis
Co Outputs
Co  rcxq   : rcxq(I,J,iw,:) = (I,J) matrix element of x0(I,J;iw)(q)
Co           x0(I,J;w)(q) = S[k=FBZ] S[t=occ] S[t'=unocc]
Co                          <M_I(q) psi(k,t) |psi(k+q,t')> <psi(k+q,t')| psi(k,t) M_J(q)>
Co                         { 1/[w-e(k+q,t')+e(k,t)+i*delta] - 1/[w+e(k+q,t')-e(k,t)-i*delta] }
Co           When npm==2 x0 is also calculated for omega<0 and stored in rcxq(:,:,:,2)
Co           rcxq is usually calculated in the M~ basis (see Remarks)
Cs Command-line switches
Cl Local variables
Cl  cphik  : amount of partial wave contributing to occ valence eigenfunction at k
Ci  cphikq : amount of partial wave contributing to unocc valence eigenfunction at k+q
Cl  ... Variables once declared but no longer used
Cl  igx    : \hat{A}^k_i is specified by symops(:,:,igx), and igxt
Cl  ppovl  : overlap of mixed basis <I|J>
Cl  nfilo  : lower bound for occupied states at k, n in in Eq. 32, PRB76, 165106
Cl  nfiup  : upper bound for occupied states at k
Cl         : OLD: nkmin,nkmax
Cl  nfilk  : no. occupied valence states at k, n in in Eq. 32, PRB76, 165106
Cl         : OLD: nt0
Cl  nemqk  : no. unoccupied states at q+k
Cl         : OLD: ntp0
Cl  nemlo  : lower bound for unoccupied states at q+k
Cl         : OLD: nkqmin (also itps)  If npm=2, nemlo should be 1.
Cl  nemup: : upper bound for unoccupied states at k+q
Cl         : OLD: nkqmax
Cl  zmelt1 : Matrix elements of mixed basis, <psi(q+k,t') | psi(k,t) M > .. Eq. 32, PRB76, 165106
Cl         : Stored as zmelt1(I(q), it1(k), it2(q+qrk))
Cl         : I  : mixed basis index
Cl         : it1: band index for occupied states t at k (more precisely rk)
Cl         : it2: band index for unoccupied t' states at q+k
Cl         : Set mkmelt.F
Cl         : RENAME as zmel0
Cl  zmelt  : Matrix elements scaled by overlap or transformed into some other
Cl         : representation, e.g E_\mu
Cl  nbmx   : max number of bands (excluding core) to include in susceptibility calculation
Cl  ncore1 : number of core states for n (1st band index of zmelt)
Cl         : and offset to first valence state in zmelt(:,ncore1+1:,:)
Cl         : OLD: nctot
Cl  ncore2 : 0 if npm is 1, ncore1 if npm is 2.
Cl         : OLD: ncc
Cr Remarks
Cr   Matrix elements and overlerlap of bare coulomb interaction are
Cr     V_IJ = <M_I|v|M_J>, O_IJ = <M_I|M_J> = ppovl
Cr     V_IJ is diagonalized by calling routine:
Cr     (V_IJ - vcoud_mu O_IJ) Z_J,mu = 0   (1)
Cr   Which implies
Cr     Z+ V Z = vcoud     (2)
Cr     Z+ O Z = 1         (3)
Cr     V Z = O Z vcoud    (4)
Cr     Z^-1 = Z+ O        (5)
Cr     V = O Z vcoud Z+ O (6)
Cr   Names of arrays:
Cr     Z_J,mu = zcousq(J, mu)
Cr     O_IJ   = <M_I|M_J> = ppovl
Cr     ppovlz_I,mu = \sum_J O_IJ Z_J,mu = O Z
Cr   so that (3) can be written
Cr     zcousq+ ppovlz = ppovlz+ zcousq = 1  (7)
Cr   and (6) be written
Cr     V_IJ  = \sum_mu ppovl(I,:)*zcousq(:;mu) vcoud_mu (Zcousq^*(:,mu) ppovl(:,J)
Cr   Kinds of mixed basis sets:
Cr     M_I  Matrix element of the as-given product basis <psi1 | psi2 M_I>.
Cr          with overlap O_IJ.  O_IJ deviates from unity only for the (PW,PW) bloch
Cr     M~_I = \sum_J O_IJ M_J  (Takao calls M~ Mbar in x0k_sym)
Cr     E_mu = is the basis which makes V diagonal:
Cr            <E_mu|v|E_nu> = delta_mu,nu vcoud_mu
Cr     E and M are related (not the funny transpose, not cc (!?)
Cr            \sum_I Z^T_I,mu M~_I = (zcousq^T M~) = (ppovlz^T M)
Cr
Cr   Susceptibility rxcq is ususally generated in the M~ basis: rcxq = rcxq_IJ
Cr   But if ebasis is nonzero, it is generated in the E basis: rcxq = rcxq_munu
Cr   It be rotated from the M~ to the E- basis
Cr     rcxq_munu = zcousq+ rcxq_IJ zcousq
Cr   The inverse rotation is, using (7)
Cr     rcxq_IJ = ppovlz rcxq_munu ppovlz+
Cr   For spin susceptibilities, a special basis is used.  In that case nmbas1 is not nmbas.
Cr
Cr   Calculation of rcxq in practice:
Cr     rcxq(iw,I,J) = sum_ibib wwk(iw,ibib)* z1p(ibib,I,J)
Cr   where
Cr     z1p = <M_I psi_it | psi_itp> < psi_itp | psi_it M_J >
Cu Updates
Cu   22 Nov 17  rcxq can be calculated in the E basis
Cu   10 Nov 14  Adapted from x0kf_v4h, which was in turn taken from Ferdi's x0k.
C ----------------------------------------------------------------------
      use m_readeigen,only: readeval
      use keyvalue,only:  getkeyvalue
!     use m_rotMPB,only:  rotMPB2
      use m_readqgcou,only: qtt_,nqnum

!     Martin Lueder's implementation of MPI
      use mpi_mod, q_comm => comm_1, q_inter => inter_1,
     &             k_comm => comm_2, k_inter => inter_2,
     &             rank => mpi_rank

      implicit none
C ... Passed parameters
      logical :: eibzmode,chipmzzr,nolfco,BSE_x0
      integer :: npm,nmb,natom,nw,isp1,isp2,nlmto,noccx,noccxv,nl,nclass,nnc,nlnmx,mdimx,
     .  nbloch,iq,nbnbx,nqbz,nmbas1,nmbas,nband,ngc,numBSE
      integer :: mdim(natom),nwgt(nqbz),nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*)
      real(8) :: q(3),qlat(3,3),rk(3,nqbz),ebasis
      real(8), intent(IN) :: ppb(nlnmx*nlnmx*mdimx*nclass)  ! ppb should not be changed in this routine.
C     indexing for occupied, unoccupied pairs
      integer :: nbnb(nqbz,npm),n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm)
C     Tetrahedron integration
      integer :: nhwtot,ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8) :: whw(nhwtot)
      complex(8) :: rcxq(nmbas1,nmbas,nw,npm)
      complex(8) :: ppovl(ngc,ngc)
      complex(8),target :: ppovlz(nmb,nmb),zzr(nmb,nmbas)

C ... Local parameters
      logical :: lprt,smbasis,onceww
      logical, parameter :: debug=.false.
C     integer ifpomat,iopen,iclose,iqxdummy,isx
      integer iatom,ibib,ipb1,ipb2,imb1,imb2,it,itp,iw,jpm,k,nbcut,nbcut2,
     .  nbmx,nctot,nfiup,nfilo,nemup,nemlo,nocc,nt0,ntp0,n1top,n2top,nbmxc,ncc
      integer :: verbose,idalloc,iverbose  ! Parameters for printout
      integer imdim(natom),iatomp(natom)
      real(8) :: imagweight,det,ebmx,fack,eband(nband)
      real(8) :: qq(3),qlatinv(3,3),symope(3,3),quu(3),shftv(3)=(/0d0,0d0,0d0/)

      complex(8) :: zq01,zq02,zq0zq0 !?This is a bug for the case of two atoms per cell!!! oct2006
      complex(8) :: phase(natom)
      complex(8),allocatable :: cphik(:,:),cphikq(:,:)
!     complex(8),allocatable :: pomat(:,:),zmeltn(:,:,:)
      complex(8),allocatable,target :: ppovl0(:,:),zmelt1(:,:,:)
      complex(8),pointer :: zebas(:,:) ! matrix elts of M basis with E- or other basis
      complex(8),pointer :: zmelt(:,:,:)

      integer :: i, j
C ... For MPI
      integer k_start, k_step, k_end
#ifdef USE_MPI
      integer num_k_groups, rcxq_size
#endif
      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/

C --- Setup ---
      lprt = verbose()>=40
      if (q_comm%ID == 0) write(*,'(" [",I4,"] x0kr begin: q=",3f8.4,$)') rank, q
      call cputid(0)
      imdim(1) = 1
      do iatom = 1,natom
        iatomp(iatom) = iatom
        if (iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
      enddo
      if (mdimx /= maxval(mdim) ) call rx('x0kr: wrong mdimx')
      if (sum(mdim(iclass(1:natom))) /= nbloch) call rx('mkmelt: product basis mismatch')
      nctot = noccx - noccxv
      ncc = 0; if (npm == 2) ncc = nctot
!     call dinv33(qlat,0,qlatinv,det)
      call minv33(qlat,qlatinv)
      phase= (1d0,0d0)
      allocate(cphik(nlmto,nband),cphikq(nlmto,nband))
      if (nbloch+ngc /= nmb) call rx('x0kr: inconsistent number of product functions')
C ... read in nbcut, nbcut2
      if (rank == 0) then
        call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
        call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      endif
#ifdef USE_MPI
      call MPI_Bcast(nbcut, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nbcut2,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif
      call getnemx(nbmx,ebmx,7,.true.) ! nbmx = band cutoff excluding cores
      nbmxc = nctot+nbmx ! nbmxc = band cutoff including cores

C     Reorder n1b, n2b so that core states precede valence states.
C     Comment out : reordering is now done external to x0kr
C     call alignn12(0,nbnbx,nqbz,0,npm,nband,nctot,ncc,nbnb,n1b,n2b)

C ... smoothed mixed basis !oct2005
C     Replace original zmelt with new zmelt based on smoothed mixed basis.
      if (smbasis()) then !need to check again, when we will make smbasis on.
        call rx('smbasis no longer supported, sorry')
C        if (ncc/=0) then
C          print *,"Timereversal=F(ncc/=0) not yet implemented for smbasis."
C          print *," pomat should be generated correctly                  ."
C          call rx( "Timereversal=F(ncc/=0) not yet implemented for smbasis.")
C        endif
C        ifpomat = iopen('POmat',0,-1,0) !oct2005
C        do
C          read(ifpomat) q_r,nn,no,iqxdummy !readin reduction matrix pomat
Cc          write(6,"('ttt: q  =',3f12.5)") q
Cc          write(6,"('ttt: q_r=',3f12.5)") q_r
C          allocate( pomat(nn,no) )
C          read(ifpomat) pomat
C          if (sum(abs(q-q_r))<1d-10) then
C            write(6,*) 'ok find the section for give qibz_k'
C            exit
C          endif
C          deallocate(pomat)
C        enddo
C        if (sum(abs(q-q_r))>1d-10 ) then
C          write(6,"('q  =',3f12.5)") q
C          write(6,"('q_r=',3f12.5)") q_r
C          call rx('POmat reading err q/=q_r')
C        endif
C        isx = iclose('POmat')
      endif

      if (q_comm%ID == 0) write(*,"(' x0kf using eibzmode =',L2)") eibzmode

C ... Big loop over k
#ifdef USE_MPI

      num_k_groups = mpi_N2

C     Split each processor group for the q-loop into num_k_groups smaller
C     groups for the k-loop.

      write(*,*) "MPI in x0kf: ",q_comm%ID, q_comm%communicator, nqbz, num_k_groups

      call mpi_make_groups(q_comm, nqbz, num_k_groups, k_comm, k_inter)

      k_start = 1 + k_comm%group_ID
      k_end   = nqbz
      k_step  = num_k_groups

!      write(0,'(" [",I4,"] q_comm%ID = ",I3,", q_comm%group_ID = ",
!     &          I3,", k_comm%ID = ",I3,", k_comm%group_ID = ",I3)')
!     &         rank, q_comm%ID, q_comm%group_ID, k_comm%ID, k_comm%group_ID

!
!     Here we can try to be more efficient by using the thdist1 routine !!
!

#else

      k_start  = 1
      k_end    = nqbz
      k_step   = 1

#endif

!ML   as this routine is summing rcxq over all k with an ALL_REDUCE,
!ML   we have to make sure that only the first processor holds the
!ML   previous (input) value:
      if (k_comm%group_ID /= 0) then
        rcxq = 0d0
      endif

      qq = 0d0
      iverbose = -1
      do  1000 k = k_start, k_end, k_step ! Perform k-integral to obtain chi0(q)
        if (eibzmode) then
          if (nwgt(k) == 0) cycle ! Equivalent points will be handled by the symmetrizer (x0k_sym)
        endif
        if (iverbose > 0) iverbose = 0
        if (iverbose < 0) iverbose = (verbose()/30)*10
        if (mod(k,10) == 1 .or. verbose() >= 40) then
#ifdef USE_MPI
        write(0,'(" [",I4,"](",I3,",",I3,") do 1000  k=",i0, " of",i0)') rank, q_comm%ID, q_comm%group_ID, k,nqbz
#else
        write(*,'(" do 1000  k=",i0, " of",i0)') k,nqbz
#endif
        endif
C       ipr = (k<5.or.k==nqbz.or.debug)
        if (sum(nbnb(k,1:npm))==0) cycle
        if (k<=5 .or. debug) then
          write(6,"(' [',I4,'] x0kr: iq, ik rk=',2i7,3f10.4,$)") rank, iq, k, rk(:,k)
          call cputid(0)
C        elseif (mod(k,max(10,nqbz/20))==1.or.k>nqbz-10) then
C          write(6,"(' x0kr: k rk=',i4,3f10.4)")k, rk(:,k)
        endif

C       Augmentation part of eigenfunctions for k and k+q
C       call readcphi(  rk(:,k)-qq,nlmto,isp1,  quu,cphik) !quu is used q for eigenfunctions.
C       call readcphi(q+rk(:,k)-qq,nlmto,isp2,  quu,cphikq)
C       call zprm('cphik',cphik,nlmto,nlmto,nband)
C       call zprm('cphikq',cphikq,nlmto,nlmto,nband)
        call getcphi(  rk(:,k)-qq,nlmto,1,nband,isp1, quu,cphik) !quu is used q for eigenfunctions.
        call getcphi(q+rk(:,k)-qq,nlmto,1,nband,isp2, quu,cphikq)
C       call zprm('cphik',cphik,nlmto,nlmto,nband)
C       call zprm('cphikq',cphikq,nlmto,nlmto,nband)

!ML debug output for  cphik (nlmto,nband) --------------------------------------------------
!         do i=1,nband
!
!            write(*,'(" [",I4,"](",I3,",",I3,")  cphik-l(j,",I4,") = ",5("(",2E11.3,") "))')
!     &           rank, iq, k, i, (cphik(j,i), j=1,min(5,nlmto))
!            write(*,'(" [",I4,"](",I3,",",I3,")  cphik-h(j,",I4,") = ",5("(",2E11.3,") "))')
!     &           rank, iq, k, i, (cphik(j,i), j=max(1,nlmto-5),nlmto)
!
!         enddo
!
!         do i=1,nband
!
!            write(*,'(" [",I4,"](",I3,",",I3,")  cphikq-l(j,",I4,") = ",5("(",2E11.3,") "))')
!     &           rank, iq, k, i, (cphik(j,i), j=1,min(5,nlmto))
!            write(*,'(" [",I4,"](",I3,",",I3,")  cphikq-h(j,",I4,") = ",5("(",2E11.3,") "))')
!     &           rank, iq, k, i, (cphik(j,i), j=max(1,nlmto-5),nlmto)
!
!         enddo
!ML
!ML  both cphik and cphikq are ok
!ML
!ML -----------------------------------------------------------------------------------------

C   --- Range, number of occupied states at k and unoccupied states at q+k ---
C       Make
C       nfilo,nfiup: lower and upper bounds for occupied states
C       nemlo,nemup: lower and upper bounds for unoccupied states
C       nt0:         "Occupied"   valence states range from 1:nt0
C       ntp0:        "Unoccupied" valence states  range from nemlo:nemlo-1+ntp0
        nfilo =  999999
        nfiup = -999999
        nemlo =  999999
        nemup = -999999
        do  jpm = 1, npm
          do  ibib = 1, nbnb(k,jpm)
C           Use these lines if cores follow valence in n1b; see alignn12
C           nfilo  = min(n1b(ibib,k,jpm),nfilo)
C           nemlo = min(n2b(ibib,k,jpm),nemlo)
C           if (n1b(ibib,k,jpm)<=nband)  nfiup = max(n1b(ibib,k,jpm),nfiup)
C           if (n2b(ibib,k,jpm)<=nband) nemup = max(n2b(ibib,k,jpm),nemup)
C           Use these lines if valence follow cores; see Remarks in alignn12
            if (n1b(ibib,k,jpm)>nctot) nfilo = min(n1b(ibib,k,jpm)-nctot,nfilo)
            if (n2b(ibib,k,jpm)>ncc) nemlo = min(n2b(ibib,k,jpm)-ncc,nemlo)
            nfiup  = max(n1b(ibib,k,jpm)-nctot,nfiup)
            nemup = max(n2b(ibib,k,jpm)-ncc,nemup)
          enddo
        enddo
!        print *, 'x0kr lims',k,nfilo,nfiup,nemlo,nemup
        if (nfilo/=1) call rx('nfilo should be 1')
        if (npm==2 .and. nemlo/=1) call rx('nemlo should be 1 (calculation has no TR symmetry)')
        call readeval(q+rk(:,k)-qq,isp2,eband)
        nemup = nocc(eband,ebmx,.true.,nband)
        if (npm==2) then
          call readeval(rk(:,k)-qq,isp1,eband)
          nfiup = nocc(eband,ebmx,.true.,nband)
        endif
        nt0  = nfiup  ! First occupied is 1
        ntp0 = nemup - nemlo + 1  ! Number of unoccupied states

C --- Calculate the matrix elements in the mixed basis <psi(k+q,t') | psi(k,t) M_I> ---
C     Product basis part:
C     zmelt1(1:nbloch,ibk,jbqk)       where:
C       ibk  = [1:nctot]               core
C       ibk  = [nctot+1:nctot+nfiup]   occ valence [nfilo,nfiup]; nfilo=1 and nfiup = highest occ state
C       jbqk = [1:nctot]               core
C       jbqk = [ncc+nemlo:ncc+nemup]   unocc valence [nemlo,nemup] ncc=0 and nemlo = 1st unocc state
C               EXCEPT if jpm=2. Then: ncc=nctot and nemlo should be 1.
C     Note for zmelt1 :
C        core states precede valence states; indexing for occ and unocc valence states start at nctot+1 and ncc+1
C     For n1b and n2b
C        reorganized so indexing consistent with of zemlt1 apart from offset nemlo in unocc channel
C   ... core
        if (debug) write(6,*) 'nbloch,nctot,nt0,ncc,ntp0=',nbloch,nctot,nt0,ncc,ntp0
        iw = idalloc('zmelt1',iverbose+2,nmb,(nctot+nt0)*(ncc+ntp0)*2)
C                       q      k        q+k
        allocate(zmelt1(nmb,ncc+ntp0,nctot+nt0)); zmelt1 = 0d0
        if (onceww(5)) print *,'nctot ncc=',nctot,ncc

        call mkmelt(0,isp2,isp1,qlat,q+rk(:,k),q,q,cphikq(1,nemlo),cphik(1,nfilo),ppb,icore,ncore,
     .    symope,shftv,[0d0],nl,natom,nclass,iclass,iatomp,nlnmx,nlnmv,nlnmc,imdim,mdim,mdimx,nlmto,
     .    ngc,nbloch,nemlo,ntp0,[0],npm,1,nt0,1,nmb,1,nctot+nt0,1,ncc+ntp0,nnc,lprt,
     o    zmelt1)

         call getkeyvalue("GWinput","x0_small_BSE",BSE_x0,default=.false.)
         if(BSE_x0)then
         write(*,*)'nctot=',nctot
         call getkeyvalue("GWinput","NumValBSE_no",numBSE,default=0)
         zmelt1(:,:,1:numBSE)=0.0d0
         endif

        if (debug) write(*,"('sum zmelt, abszmelt=',1p,4d23.16)") sum(zmelt1),sum(abs(zmelt1))
C       call zprm('zmelt1',zmelt1,nmb,nmb,(nctot+nt0)*(ncc+ntp0))

C   ... Smooth mixed basis Oct 2005 (not maintained now)
C       Replace original zmelt with new zmelt based on smoothed mixed basis.
        if (smbasis()) then
          call rx('smbasis not maintained recently need to check again')
C          if (nn /= nmb) call rx('x0kr: nmb is not nmb')
C          if (no /= nmb)  call rx('x0kr: POmat err no/=nmb')
C          nmb = nn ! Renew nmb
C          iw = idalloc('zmeltn',iverbose+2,nmb,(nctot+nt0)*ntp0*2)
C          allocate(zmeltn(nmb, 1:nctot+nt0, 1:ntp0))
C          call matm(pomat, zmelt1, zmeltn, nn,no,(nctot+nt0)*ntp0)
C          deallocate(zmelt1)
C          allocate(zmelt1(nmb, 1:nctot+nt0, 1:ntp0))
C          zmelt1 = zmeltn
C          deallocate(zmeltn)
C          iw = idalloc('zmeltn',iverbose+4,nmb,(nctot+nt0)*ntp0*2)
c        write(6,"('smbasis: nmb nmb nn no=',4i5)") nmb,nmb,nn,no
        endif

C   --- Convert zmelt to E-\mu-basis representation, or some other basis ---
        iw = idalloc('zmelt',iverbose+2,nmbas,(nctot+nt0)*(ncc+ntp0)*2)
        if (chipmzzr .and. nolfco) then
          zebas => zzr          ! to rotate to spin moment basis
        elseif (chipmzzr) then
          call rx('not ready for chiPM without local fields')
C       Rotates M to E-mu basis; see Remarks
        elseif (nmbas1 /= nmb .or. nmbas /= nmb .or. ebasis /= 0d0) then
          zebas => ppovlz
C         call zprm('ppovlz',ppovlz,nmb,nmb,nmb)
C       Rotates M to M~ basis; see Remarks
        else
          nullify(zebas)
          allocate(zmelt(nmbas,ncc+ntp0,nctot+nt0))

          do j = 1, nctot+nt0
            do i = 1, ncc+ntp0
              call zcopy(nbloch, zmelt1(1,i,j), 1, zmelt(1,i,j), 1)
            end do
          end do
!           zmelt(1:nbloch,:,:) = zmelt1(1:nbloch,:,:)
          call zgemm('T','N',ngc,(nctot+nt0)*(ncc+ntp0),ngc,(1d0,0d0),ppovl,ngc,
     .      zmelt1(nbloch+1,1,1),nmb,(0d0,0d0),zmelt(nbloch+1,1,1),nmbas)

C         Also accomplishes rotation to M~ but less efficient since zgemm spans more functions
C         call zprm('ppovl',ppovl,ngc,ngc,ngc)
C         call zprm('ppovlz',ppovlz,nmb,nmb,nmb)
C         allocate(ppovl0(nmb,nmb))
C         call dinit(ppovl0,2*nmb*nmb)
C         do  ix = 1, nbloch; ppovl0(ix,ix) = 1; enddo
C         ppovl0(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc) = ppovl
C         zebas => ppovl0       ! Standard mixed basis
        endif
!       call zprm('zebas',zebas,nmb,nmb,nmb)

C       Make zmelt in E basis = ppovlz^T * M or in spin moment basis = zzr^T * M
        if (associated(zebas)) then ! Do basis conversion zmelt = zebas * M
          allocate(zmelt(nmbas,ncc+ntp0,nctot+nt0))
C         call zprm('zebas',zebas,nmb,nmb,nmb)
          call zgemm('T','N',nmbas,(nctot+nt0)*(ncc+ntp0),nmb,(1d0,0d0),zebas,nmb,zmelt1,nmb,(0d0,0d0),zmelt,nmbas)
C         call zprm('E = ppovlz^T * M',zmelt,nmb,nmb,(nctot+nt0)*(ncc+ntp0))
        endif
        deallocate(zmelt1)
        iw = idalloc('zmelt1',iverbose+4,nmb,(nctot+nt0)*(ncc+ntp0)*2)
        if (verbose()>=30 .and. k == 1 .or. lprt) call cputid2(' complete zmelt * overlap',0)
!       debug = .true.
        if (debug) write(*,'(" sum abs zmelt*ovlp",3d13.6)') sum(abs(zmelt))
        if (associated(zebas,ppovl0)) then
            nullify(zebas)
            deallocate(ppovl0)
        end if

C   --- Calculate chi0 without local fields ---
        if (nolfco) then !iepsmode==202) then ! just for <exp(iq r)|x0(q,\omega)|exp(iq r)>
          n1top = nctot+nfiup ! Last occupied state to include in chi0
          n2top = ncc+nemup  ! Last unoccupied state to include in chi0
          do  jpm = 1,npm !
            do  ibib = 1, nbnb(k,jpm) !---  ibib loop

C              if (jpm == 1) then
C                if (n1b(ibib,k,jpm)<=nbcut .and. n2b(ibib,k,jpm)>nbcut2) cycle
C              else
C                if (n2b(ibib,k,jpm)<=nbcut .and. n1b(ibib,k,jpm)>nbcut2) cycle
C              endif
C
C              if (jpm==1 .and. n2b(ibib,k,jpm)>nbmx) cycle
C              if (jpm==2 .and. n1b(ibib,k,jpm)>nbmx) cycle
C
C              if (n1b(ibib,k,jpm) <= nband) then
C                it = nctot + n1b(ibib,k,jpm) !valence
C              else
C                it = n1b(ibib,k,jpm) - nband !core
C              endif
C
C              if (n2b(ibib,k,jpm) <= nband) then
C                itp = ncc + n2b(ibib,k,jpm) - nemlo + 1 ! valence
C                if (itp > ncc + nemup-nemlo+1) cycle
C              else
C                itp = n2b(ibib,k,jpm) - nemlo + 1 - nband !core
C              endif

C             ... This variant assumes valence states precede cores in n1b and n2b
              if (n1b(ibib,k,jpm) > n1top .or. n2b(ibib,k,jpm) > n2top) cycle
              if (jpm==1 .and. n2b(ibib,k,jpm)>nbmxc) cycle ! Past highest unocc state
              if (jpm==2 .and. n1b(ibib,k,jpm)>nbmxc) cycle
              if (nbcut > 0) then ! Special handling if nbcut has been set
                if (jpm == 1 .and. n1b(ibib,k,jpm)<=nbcut .and. n2b(ibib,k,jpm)>nbcut2) cycle
                if (jpm == 2 .and. n2b(ibib,k,jpm)<=nbcut .and. n1b(ibib,k,jpm)>nbcut2) cycle
              endif
              it = n1b(ibib,k,jpm) !Index to matrix elements of occ states, with first occ state = nctot+1
              itp = n2b(ibib,k,jpm) - nemlo + 1 !Index to matrix elements start with first unocc state = ncc+1

              do imb1 = 1, nmbas1
              do imb2 = 1, nmbas
                zq01 = zmelt(imb1,itp,it)
                zq02 = zmelt(imb2,itp,it)
                zq0zq0 = dconjg(zq01)*zq02
                do  iw = ihw(ibib,k,jpm), ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
                  if (iw > nw) call rx('x0kr: iw > nw')
                  imagweight   = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
                  if (eibzmode) imagweight = nwgt(k)*imagweight
                  rcxq(imb1,imb2,iw,jpm) = rcxq(imb1,imb2,iw,jpm) + zq0zq0*imagweight ! sum over ibib (i.e n, n') and k.
                enddo           ! iw
              enddo ! imb2
              enddo ! imb1
            enddo ! ----- ibib loop
          enddo ! ----- jpm loop
          iw = idalloc('zmelt',iverbose+4,nmb,(nctot+nt0)*(ncc+ntp0)*2)
          deallocate(zmelt) !,z1p)
          if (lprt) call cputid2(' complete generation of rcxq',0)
          cycle  !to do 1000
        endif

C   --- Make x0(I,J;iw)(q) ---
        if (nmbas /= nmb .or. nmbas1 /= nmb) call rx('x0kr: product basis dimension mismatch')
        fack = 1; if (eibzmode) fack = nwgt(k)
        call mkzxq(nqbz,nbnbx,npm,k,nband,nctot,ncc,nbmx,nbcut,nbcut2,nmb,nemlo,nfiup,nemup,nbnb,
     .    n1b,n2b,ihw,nhw,jhw,whw,nw,fack,1,nmb,1,nctot+nt0,1,ncc+ntp0,zmelt,rcxq)

        deallocate(zmelt)
        if (lprt) call cputid2(' complete generation of rcxq',0)
        iw = idalloc('zmelt',iverbose+4,nmb,(nctot+nt0)*(ncc+ntp0)*2)

        if (debug) write(*,'(" sum rcxq",3d13.6,$)') sum(abs(rcxq))
        if (debug) call cputid(0)
        if (debug) write(6,*)' end of kloop k jpm=',k,jpm
C       call zprm('rcxq',rcxq,nmbas1*nmbas1,nmbas*nmbas,nw)

 1000 continue
      if (verbose()>=30) call cputid2(' completed k-loop',0)

!------------------------------------------------------------------------------------------
! MPI remark:
!
! At this point, each processor holds the sum over its set of k-contributions to rcxq.
! After the k-loop, we need to sum each element (ipb1,ipb2,iw) over the different k-sets!!
!------------------------------------------------------------------------------------------
#ifdef USE_MPI

!      write(0,'(" [",I4,"](",I3,",",I3,";",I3,",",I3,") Before MPI_Barrier 1000")')
!     &     rank, q_comm%group_ID, q_comm%ID, k_comm%group_ID, k_comm%ID
!
!
!      CALL MPI_BARRIER(k_inter%communicator,ierror)
!
!      write(0,'(" [",I4,"](",I3,",",I3,";",I3,",",I3,") After MPI_Barrier 1000")')
!     &     rank, q_comm%group_ID, q_comm%ID, k_comm%group_ID, k_comm%ID


!      write(0,'(" [",I4,"](",I3,",",I3,";",I3,",",I3,") , k-inter = ", 2I5)')
!     &     rank, q_comm%group_ID, q_comm%ID, k_comm%group_ID, k_comm%ID,
!     &     k_inter%group_ID, k_inter%ID

      rcxq_size = SIZE(rcxq,1)*SIZE(rcxq,2)*SIZE(rcxq,3)*SIZE(rcxq,4)

!      write(0,'(" [",I4,"](",I3,",",I3,";",I3,",",I3,") Before MPI_ALLREDUCE: rcxq ",3I10)')
!     &     rank, q_comm%group_ID, q_comm%ID, k_comm%group_ID, k_comm%ID,
!     &     nmb*nmb*nw*npm, nw*npm*nmbas1*nmbas2, rcxq_size


      CALL MPI_ALLREDUCE(MPI_IN_PLACE, rcxq, rcxq_size,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM,
     &     k_inter%communicator, ierror)


      if (ierror /= 0) write(*,*) "Something wrong with AllReduce rcxq:",rank
!     write(*,*) "rcxq(1,1,1,1) = ",rcxq(1,1,1,1)

#endif

C --- Hermitianize. jun2012takao moved from dpsion5 ---
      if (.not.nolfco) then
        if (nmbas1/=nmbas) then
          call rxi('x0kf: symmetrization requires nmbas1=nmbas but nmbas1=',nmbas1)
        endif
        do jpm = 1, npm
          do iw = 1, nw
            do ipb2 = 1, nmbas
              do ipb1 = 1, ipb2-1
                rcxq(ipb2,ipb1,iw,jpm) = dconjg(rcxq(ipb1,ipb2,iw,jpm))
              enddo
            enddo
          enddo
        enddo
      endif

      deallocate(cphik,cphikq)
!     if (smbasis()) deallocate(pomat)
      if (lprt) then
        iw = idalloc(' ',11,1,1)
      else
        iw = idalloc(' ',1,1,1)
        write(*,"(/' Max dynamic allocation x0kr:',i7,' MB'/)") iw
      endif

C      iw = idalloc(' ',1,1,1)
C      write(*,"(/' Max dynamic allocation x0kr:',i7,' MB'/)") iw

      print *,'x0kr: end'
      end subroutine x0kr

      subroutine dpsion5(frhis,nwhis,freqr,nw_w,freqi,niwt,realomega,imagomega,
     i  rcxq,npm,nw_i,nmbas1,nmbas2,
     o  zxq,zxqi,
     i  chipm,schi,isp,ecut,ecuts)
C- Make full noninteracting susceptibility on the real or imaginary axis from the spectrum weight rcxq
C ----------------------------------------------------------------------
Ci Inputs
Ci   frhis : histogram bins for frequencies on the real axis.
Ci         : The i-th bin is the interval [frhis(i), frhis(i+1)].
Ci         : (Routine hilbertmat is not itself limited by this condition).
Ci   nwhis : number of bins.  frhis has nwhis+1 elements.
Ci   freqr : freqr(0:nw_w) = table of real energies for which to calculate zxq
Ci         : We take freqr(i) to be the middle of i-th bin.  Note freqr(0)=0.
Ci   nw_w  : dimension of freqr
Ci   freqi : freqi(1:niwt) = table of imagiary energies for which to calculate zxqi
Ci   niwt  : dimension of freqi
Ci   realomega : T => calculate zxq
Ci   imagomega : T => calculate zxqi
Ci   rcxq  : spectrum weight for given in frhis along the real-axis (see x0kr)
Ci         : i.e. Im (chi0(omega)) for each of the frequencies in frhis
Ci         : rcxq may be OVERWRITTEN by this routine
Ci   npm   : 1 => timereversal=T; 2 => timereversal=F
Ci   nw_i  : Make zxq for frequencies nw_i:nw_w
Ci   nmbas1: Number of mixed basis functions : rcxq = rcxq(nmbas1,nmbas2,nwhis,npm)
Ci   nmbas2: Number of mixed basis functions : rcxq = rcxq(nmbas1,nmbas2,nwhis,npm)
Ci   chipm : If T => caculating transverse susceptibility
Ci   schi  : ?
Ci   isp   :current spin channel (1 or 2)
Ci   ecut  : ?
Ci   ecuts : ?
Co Outputs
Co   zxq   : (if realomega=T) W-v along the real axis on freqr(0:nw_w)
Co         : If chipm is also true, result is added (not stored) in zxq
Ci   zxqi  : (if imagomega=T) W-v along the imaginary axis on freqi(niwt)
Cr Remarks
Cr   v4 works for timereversal=F (npm=2 case).
Cu Updates
Cu   Feb2006:  v4 for timereversal=F
Cu   July2005: v3Add spin chipm mode
Cu   sergey faleev Apr 2002 ; Rebuilt by takao
C------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer:: nw_i,nw_w,niwt,nwhis,npm,nmbas1,nmbas2
      real(8) :: freqi(niwt),frhis(nwhis+1),freqr(0:nw_w)
      logical :: realomega,imagomega
      complex(8):: rcxq(nmbas1,nmbas2,nwhis,npm)
      logical :: chipm
      integer::isp
      real(8):: ecut,ecuts,schi
      complex(8):: zxq(nmbas1,nmbas2,nw_i:nw_w) !iw=0 corresponds to omg=0
      complex(8):: zxqi(nmbas1,nmbas2,niwt)
C ... Local parameters
      logical lerr
      integer iw,ispx,it,jpm,i,j
      real(8) :: pi,aaa,wcutef
      logical:: debug=.false.
      real(8),allocatable :: his_L(:),his_R(:),his_C(:)  !,ebb(:)
      real(8),allocatable ::rmat(:,:,:),rmati(:,:,:),rmatt(:,:,:),imatt(:,:,:)
      complex(8),allocatable :: rmatiC(:,:,:),imattC(:,:,:),rrr(:)
      complex(8):: wfac,zz,img,alpha,beta

      if (debug) then
        write(6,*)' dpsion5 nmbas1 nmbas2 nwhis npm =',nmbas1,nmbas2,nwhis,npm
        write(6,*)' sumchk rcxq=',sum(abs(rcxq))
      endif
      pi  = 4d0*datan(1d0)
      img = (0d0,1d0); alpha=1d0; beta=0d0
      call cputid(0)
      ispx = isp
      if (schi<0) then
        ispx = 3-isp !flip
      endif

C ... Check freqr = frhis_m.
C     But I think now this is not necessary. You can supply any freqr(iw). But be careful.
      if (realomega) then
        if (nwhis < nw_w) then
          print *,nwhis,nw_w
          call rx(' dpsion5: nwhis<nw_w')
        endif
        if (freqr(0)/=0d0) call rx(' dpsion5: freqr(0)/=0d0')
        aaa = 0d0
        if (nw_w>0) then
          do iw = 1,nw_w
            aaa = aaa + abs(freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0)
            if (debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)")
     &        iw,freqr(iw),(frhis(iw)+frhis(iw+1))/2d0
          enddo
          if (aaa>1d-10) call rx('dpsion5: freqr/=frhis_m is not implemented yet')
        endif
      endif                     !realomega

C--------------------------------------------------------------
! Each bin  [his_Left, his_Right]  his_Center is middle.
! his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
      allocate(his_L(-nwhis:nwhis),his_R(-nwhis:nwhis),his_C(-nwhis:nwhis))
      his_L(1:nwhis) = frhis(  1:  nwhis)
      his_R(1:nwhis) = frhis(1+1:1+nwhis)
      his_C(1:nwhis) = (his_L(1:nwhis) + his_R(1:nwhis))/2d0
      do iw= 1,nwhis
        his_L(-iw) = -his_R(iw)
        his_R(-iw) = -his_L(iw)
        his_C(-iw) = -his_C(iw)
      enddo
      his_C(0) = 0d0; his_R(0)=-999; his_L(0)=-999

      do iw= 1, nwhis
        if (ecut<1d9) then
          wfac = wcutef(his_C(iw),ecut,ecuts)
        else
          wfac = 1d0
        endif
!       rcxq is used as work---> rcxq= Average value of Im chi.
!       Note scale rcxq by -1 : "negative" of physical quantity
        do jpm = 1, npm
          call dscal(2*nmbas1*nmbas2,-wfac/(his_r(iw)-his_l(iw)),rcxq(1,1,iw,jpm),1)
        enddo
      enddo

C      if (nmbas1==nmbas2) then
C        write(6,"('--- dpsion5 : EigenValues for rcxq --------')")
C        allocate(ebb(nmbas1))
C        do jpm= 1,npm
C          do iw = 1, nwhis
C            call diagcvh2(rcxq(:,:,iw,jpm),nmbas1,ebb)
C            do ii=1,nmbas1
C              write(6,"('hhh1: xxxxxxxxxxxxxxxxx',2i4)") jpm,iw
C              if (abs(ebb(ii))>1d-8.and.ebb(ii)>0)
C     &      write(6,"('hhh1: jpm iw eb=',2i4,d13.5)") jpm,iw,ebb(ii)
C            enddo
C          enddo
C        enddo
C        deallocate(ebb)
C      endif

C --- realomega case ---
      if (realomega)then
        write(6,*) " --- realomega --- "
        if (npm==1) then
          allocate(rmat(0:nw_w,-nwhis:nwhis,npm), rrr(-nwhis:nwhis))
          rmat  = 0d0
          do it = 0, nw_w
            zz = freqr(it)      !his_C(it)
            call hilbertmat(zz,nwhis,his_L,his_C,his_R,rrr)
            rmat(it,:,1) = dreal(rrr)
          enddo
          allocate(rmatt(0:nw_w,nwhis,npm))
          if (chipm .and. ispx==1) then
            rmatt(:,:,1) = rmat(:,1:nwhis,1)
          elseif (chipm .and. ispx==2) then
            do iw= 1,nwhis
              rmatt(:,iw,1) = -rmat(:,-iw,1)
            enddo
          else
            do iw= 1,nwhis
              rmatt(:,iw,1) = rmat(:,iw,1) - rmat(:,-iw,1)
            enddo
          endif
          deallocate(rmat,rrr)
        else                    ! npm==2 case
          allocate(rmatt(-nw_w:nw_w,nwhis,npm),rrr(-nwhis:nwhis))
          rmatt = 0d0
          do it = -nw_w, nw_w
            if (it<0) then
              zz = -freqr(-it)  !his_C(it)
            else
              zz = freqr(it)    !his_C(it)
            endif
            call hilbertmat(zz,nwhis,his_L,his_C,his_R,rrr)
            rmatt(it,:,1) =  dreal(rrr  (1:nwhis))
            rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))
          enddo
          deallocate(rrr)
        endif
        rmatt = rmatt/pi

        if (npm==2 .and. chipm) call rx('x0kr: npm==2.and.chipm is not meaningful probably')

!!      Note rcxq is now negative (converted at the top of this routine)
        if (chipm .and. ispx==2) then
                                !nothing here
                                !Since the range of zxq is nw_i=0, we have no area to store negative energy part of chipm.
        elseif (chipm) then ! Add i*rcxq to zxq(:,:,1:nw_w)
C         Product nmbas1*nmbas2*nw_w may overflow: replace this call with a loop over iw
C         call zaxpy(nmbas1*nmbas2*nw_w,img,rcxq,1,zxq(1,1,1),1)
          do  iw = 1, nw_w
            call zaxpy(nmbas1*nmbas2,img,rcxq(1,1,iw,1),1,zxq(1,1,iw),1)
          enddo
        else ! Copy i*rcxq to zxq
          zxq = 0d0
C         Product nmbas1*nmbas2*nw_w may overflow: replace this call with a loop over iw
C         call zaxpy(nmbas1*nmbas2*nw_w,img,rcxq,1,zxq(1,1,1),1)
          do  iw = 1, nw_w
            call zaxpy(nmbas1*nmbas2,img,rcxq(1,1,iw,1),1,zxq(1,1,iw),1)
          enddo
        endif
        if (npm==2) then  ! No TR symmetry: copy rcxq(:,:,:,2) to zxq(-iw)
          do  iw = 1, nw_w
            call zaxpy(nmbas1*nmbas2,img,rcxq(1,1,iw,2),1,zxq(1,1,-iw),1)
          enddo
        endif

        if (npm==1) then
          call dgemm('N','T',2*nmbas1*nmbas2,nw_w+1,nwhis,1d0,rcxq,
     .      2*nmbas1*nmbas2,rmatt,nw_w+1,1d0,zxq,2*nmbas1*nmbas2)
        elseif (npm==2) then
          call dgemm('N','T',2*nmbas1*nmbas2,npm*nw_w+1,nwhis,1d0,rcxq(1,1,1,1),
     .      2*nmbas1*nmbas2,rmatt(:,:,1),npm*nw_w+1,1d0,zxq,2*nmbas1*nmbas2)
          call dgemm('N','T',2*nmbas1*nmbas2,npm*nw_w+1,nwhis,1d0,rcxq(1,1,1,2),
     .      2*nmbas1*nmbas2,rmatt(:,:,2),npm*nw_w+1,1d0,zxq,2*nmbas1*nmbas2)
        else
          call rx('dpsion5: npm must be 1 or 2')
        endif
        deallocate(rmatt)
      endif

C --- imagomega case ---
      if (imagomega) then
        allocate(rrr(-nwhis:nwhis))
        if (npm==1) then
          allocate(rmati (niwt,-nwhis:nwhis,npm))
          rmati= 0d0
        else
          allocate(rmatiC(niwt,-nwhis:nwhis,npm))
          rmatiC = 0d0
        endif
        do  it = 1, niwt
          zz = img*freqi(it)    !his_C(it)
          call hilbertmat(zz,nwhis,his_L,his_C,his_R,rrr) !Im(zz)>0
          if (npm==1) then
            rmati (it,:,1) = dreal(rrr)
          else
            rmatiC(it,:,1) = rrr
          endif
        enddo

        if (npm==1) then
          allocate(imatt(niwt,nwhis,npm))
          do  iw = 1, nwhis
            imatt(:,iw,1) = rmati(:,iw,1) - rmati(:,-iw,1)
          enddo
          deallocate(rmati,rrr)
          imatt = imatt/pi
          call dgemm('N','T',2*nmbas1*nmbas2,niwt,nwhis,1d0,rcxq,
     .      2*nmbas1*nmbas2,imatt,niwt,0d0,zxqi,2*nmbas1*nmbas2)
          lerr = .false.
          do  i = 1, nmbas1
            do  j = 1, nmbas2
              if (lerr) exit
              do it = 1, niwt
                lerr = lerr .or. (zxqi(i,j,it)+0 /= zxqi(i,j,it))
C                if (lerr)
C                  print *, i,j,it
C                  call rx('found error')
C                endif
              enddo
            enddo
          enddo
          if (lerr) then
            write(*,'(" !!! dpsion5 error in dgemm call ... redo with zgemm")')
            allocate(imattC(niwt,nwhis,npm))
            imattC = imatt
            call zgemm('N','T',nmbas1*nmbas2,niwt,nwhis,alpha,rcxq,
     .        nmbas1*nmbas2,imattC,niwt,beta,zxqi,nmbas1*nmbas2)
            deallocate(imattC)
          endif
          deallocate(imatt)
        else
          allocate(imattC(niwt,nwhis,npm))
          do iw= 1,nwhis
            imattC(:,iw,1) =   rmatiC(:,iw,1)
            imattC(:,iw,2) = - rmatiC(:,-iw,1)
          enddo
          deallocate(rmatiC,rrr)
          imattC = imattC/pi
          call zgemm('N','T',nmbas1*nmbas2,niwt,nwhis,alpha,rcxq(1,1,1,1),
     .      nmbas1*nmbas2,imattC(1,1,1),niwt,beta,zxqi,nmbas1*nmbas2)
          call zgemm('N','T',nmbas1*nmbas2,niwt,nwhis,alpha,rcxq(1,1,1,2),
     .      nmbas1*nmbas2,imattC(1,1,2),niwt,alpha,zxqi,nmbas1*nmbas2)
          deallocate(imattC)
        endif
      endif

      deallocate(his_L,his_C,his_R)
      write(*,'(" ... end dpsion5 ",$)')
      call cputid(0)
      end

C      logical function checkbelong(qin,qall,nq,ieibz) !ieibz is also returned
C      integer:: nq,ieibz
C      real(8):: qin(3),qall(3,nq)
C
C      checkbelong=.false.
C      do i=1,nq
C        if (sum(abs(qin-qall(:,i)))<1d-8) then
C          ieibz=i
C          checkbelong=.true.
C          return
C        endif
C      enddo
C      end
