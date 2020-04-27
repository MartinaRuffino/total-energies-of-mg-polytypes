C#define OPENMP
      subroutine mkzxq(nqbz,nbnbx,npm,k,ndham,ncore1,ncore2,nbmx,
     .  nbcut,nbcut2,nmb,nkqmin,nkmax,nkqmax,nbnb,n1b,n2b,ihw,nhw,
     .  jhw,whw,nw,fack,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,zmelt,rcxq)
C- Imaginary part of bare response function
C ----------------------------------------------------------------------
Ci  nqbz   :number of k-points in the 1st BZ
Ci  nbnbx  :dimensions n1b,n2b,ihw,nhw,jhw
Ci         :Must be at least as large as the max value of nbnb
Ci  npm    :1 if time reversal symmetry, 2 if not
Ci  k      :index to current kpoint
Ci  ndham  :number of energy bands (no core states)
Ci  ncore1 :number of core states for n (1st band index of zmelt)
Ci         :and offset to first valence state in zmelt(:,ncore1+1:,:)
Ci         :OLD: nctot
Ci  ncore2 :Number of core states for n' (1nd band index of zmelt) and should be eiher:
Ci         :0      => do not calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :ncore1 => do calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :OLD: ncc
Ci  nbmx   :max number of bands (excluding core) to include in calculation
Ci         :(value of nband_chi0 in GWinput)
Ci  nbcut  :exclude contribution to chi0 from pairs with
Ci         :occ valence state <= nbcut AND unoccupied valence state > nbcut2
Ci         :(nbcut read as nbcutlow in GWinput)
Ci  nbcut2 :exclude contribution to chi0 from pairs with
Ci         :occ valence state <= nbcut AND unoccupied valence state > nbcut2
Ci         :(nbcut2 read as nbcutlowto in GWinput)
Ci  nmb    :number of mixed product basis functions entering into calculation of chi0
Ci         :nmb also dimensions rcxq and zmelt
Ci         :(formerly ngb)
Ci  nfilo  :include contribution from occupied states nfilo:nfiup at k
Ci  nfiup  :include contribution from occupied states nfilo:nfiup at k
Ci         :OLD: nkmin,nkmax
Ci  nemlo  :lower bound for unoccupied states at k+q
Cl         :OLD: nkqmin (also itps)  If npm=2, nkqmin should be 1.
Cl  nemup: :upper bound for unoccupied states at k+q
Cl         :OLD: nkqmax
Ci  nbnb   :nbnb(k) = number of nonzero (it,jb) pairs for kth qp
Ci  n1b    :n1b(ix,k) = ib for ix^th nonzero (ib,jb) pair
Ci         :ib=1 = first valence state.
Ci         :When core states are included they follow the valence states
Ci  n2b    :n2b(ix,k)= jb for ix^th nonzero (ib,jb) pair
Ci         :jb=1 = first valence state.
Ci         :When core states are included they follow the valence states
Ci  ihw    :ihw(nbnb(k),k) = index to first histogram bin within
Ci                           (demin,demax) for a given nbnb pair and k
Ci  nhw    :nhw(nbnb(k),k) = number of histograms bins witin
Ci                           (demin,demax) for a given nbnb pair and k
Ci  jhw    :jhw(nbnb(k),k) = pointer to whw; see whw
Ci  whw    :whw(i:j) histogram weights in bins i:j for given ib,jb,kx
Ci         :i = jhw(nbnb,kx)
Ci         :j = jhw(nbnb,kx)+nhw(nbnb),kx)-1
Ci         :nbnb = nbnb(ib,jb,kx)
Ci  nw     :number of frequencies where x0 is calculated.
Ci         :(formerly nwhis)
Ci  fack   :BZ weighting for this k-point
Ci  zmelt  :<psi(q+k,t') | psi(k,t) M(I)>, I=1..nmb
Co Outputs
Co   rcxq  :Noninteracting response function x0 on the positive real axis
Co         :A portion of rxcq is accumulated for one spin.
Cr Remarks
Cr   zxq(i,j,iw) =  (-1/pi) \int_{freqr(iw)}^{freqr(iw+1)} d\omega
Cr                           Im <M(i)|x0(q,\omega)|M(j)>
Cr   Evaluate and store in rcxq as follows:
Cr    z1p = <M_I psi_it | psi_itp> < psi_itp | psi_it M_J >
Cr    For a given (occ,unocc) pair.
Cr    z1p(I,J) =  dconjg(zmelt(I,occ,unocc)) * zmelt(J,occ,unocc)
Cr    Index a particular (occ,unocc) pair by ibib.  Then
Cr    rcxq(iw,I,J) = sum_ibib weight(iw,ibib) * dconjg(zmelt(I,occ,unocc)) * zmelt(J,occ,unocc)
Cr
Cr   This routine first identifies all pairs ibib belonging to a particular bin iw (pmkzxq)
Cr   For each iw, it evaluates
Cr      sum_ibib weight(iw,ibib) * dconjg(zmelt(I,[ibib])) * zmelt(J,[ibib])
Cr   using the BLAS3 library routine zherk. Note that zherk can be multithreaded.
Cb Bugs
Cb   Can we arrange for n1b and n2b to have the same core-valence ordering
Cb   as in the matrix elements?
Cu Updates
Cu   04 Jun 16 (DMT) rewrote loop order to enable BLAS3 zherk call
Cu   15 May 13 Adapted from x0kf_v4h.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,intent(in) :: ipbhi,ipblo,it1hi,it1lo,it2hi,it2lo,nkqmin,k,ndham,
     .  nbcut,nbcut2,nbmx,nbnbx,ncore2,ncore1,nkmax,nkqmax,nmb,npm,nqbz,nw
      integer,intent(in) :: nbnb(nqbz,npm),n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm),
     .  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8),intent(in) :: whw(*),fack
      complex(8) :: rcxq(nmb,nmb,nw,npm)
C     Original index order
C     complex(8) :: zmelt(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
C     Reversed index order
      complex(8) :: zmelt(ipblo:ipbhi,it2lo:it2hi,it1lo:it1hi)

C ... Allocatable arrays
      integer,allocatable :: pibiw(:,:,:)
      complex(8),allocatable :: ztmp(:,:)
C ... Local parameters
      integer ibib,jpm,it,itp,iw,ndibiw,ipb
      real(8) :: imagweight,minwt,maxwt
      real(8), parameter :: wcut = 1d-10
C     integer(8) :: c1, c2, cr, cm
      complex(8) :: zfac
      integer :: nibiw(nw,npm), itc, iwa(1)

C#ifdefC TIMING
C      call system_clock(count_rate = cr, count_max = cm)
C      call system_clock(c1)
C#endif

C --- Identify ibib pairs belonging to each frequency bin ---
C     iwa(1) = iw
      ndibiw = 1
      call pmkzxq(0,nqbz,nbnbx,npm,k,ncore1,ncore2,nbmx,nbcut,nbcut2,
     .  nkqmin,nkmax,nkqmax,nbnb,n1b,n2b,ihw,nhw,nw,ndibiw,nibiw,iwa)
      ndibiw = maxval(nibiw)
C     Original index order
C     allocate(pibiw(ndibiw,nw,npm),ztmp(ndibiw,nmb))
C     Reversed index order
      allocate(pibiw(ndibiw,nw,npm),ztmp(nmb,ndibiw))
      call pmkzxq(1,nqbz,nbnbx,npm,k,ncore1,ncore2,nbmx,nbcut,nbcut2,
     .  nkqmin,nkmax,nkqmax,nbnb,n1b,n2b,ihw,nhw,nw,ndibiw,nibiw,pibiw)

C --- jpm=1 if TR symmetry; Loop over jpm=1,2 if not ---
      maxwt = -99999; minwt = 0
      do  jpm = 1, npm
C   --- Sum over (n, n') pairs for this k ---
        do  iw = 1, nw
          if (nibiw(iw,jpm) == 0) cycle

C         Original index order
C         ztmp = 0
          do  itc = 1, nibiw(iw,jpm)
            ibib = pibiw(itc,iw,jpm)
            it = n1b(ibib,k,jpm) !valence, with first occ state = ncore1+1
            itp = n2b(ibib,k,jpm) - nkqmin + 1 !valence, with first unocc state = ncore2+1

            imagweight = fack * whw(jhw(ibib,k,jpm) - ihw(ibib,k,jpm) + iw)
            maxwt = max(imagweight,maxwt); minwt = min(imagweight,minwt)
            if (imagweight < 0d0) then
              imagweight = 0d0  ! For debugging
              if (-imagweight < wcut) imagweight = 0d0
              if (imagweight < 0d0) then
                print *, 'mkzxq: negative weight', imagweight, ibib, it, itp, iw
                call rx('mkzxq: negative weight')
              endif
            endif
            zfac = sqrt(imagweight)

C           Original index order. We could use unit stride here if we generate transpose of rcxq.
C           call zaxpy(nmb, zfac, zmelt(ipblo,it,itp), 1, ztmp(itc,1), ndibiw)
C           Reversed index order
            forall (ipb=1:nmb) ztmp(ipb,itc) = zfac*conjg(zmelt(ipblo+ipb-1,itp,it))
        end do

C         Original index order
C         call zherk('u','c',nmb,nibiw(iw,jpm),1d0,ztmp,ndibiw,1d0,rcxq(1,1,iw,jpm),nmb)
C         Reversed index order
          call zherk('u','n',nmb,nibiw(iw,jpm),1d0,ztmp,nmb,1d0,rcxq(1,1,iw,jpm),nmb)

        end do                  ! iw
      enddo                     ! jpm

      if (minwt < 0)
     .  write(*,'(" mkzxq (warning!) negative weight",1pe10.2," as fraction of max",1pe10.2)') minwt,minwt/maxwt

C#ifdefC TIMING
C      call system_clock(c2)
C      write(*,'("mkzxq:",x,g20.12)') real(c2 - c1,8)/real(cr,8)
C      write(*,'("sum abs rcxq",1pd15.6,0p,2f15.10)') sum(abs(rcxq)), sum(rcxq(:,2,2,1))
C#endif

      end

      subroutine alignn12(job,nbnbx,nqbz,kshft,npm,ndham,ncore1,ncore2,nbnb,n1b,n2b)
C- Reorder precedence of (core,valence) states in n1b and n2b for 1 k, or all k's
C ----------------------------------------------------------------------
Ci Inputs
Ci  job    :reverse (core,valence) ordering in n1b and n2b, with :
Ci         :0: valence states precede core states on input
Ci         :1: core states precede valence states on input
Ci  nbnbx  :dimensions n1b,n2b,ihw,nhw,jhw
Ci         :Must be at least as large as the max value of nbnb
Ci  nqbz   :number of k-points in the 1st BZ
Ci  kshft  :which k-point in n1b, n2b to effect the reordering
Ci         :Use kshft=0 to reorder pairs for all 1:nqbz k-points
Ci  npm    :1 if time reversal symmetry, 2 if not
Ci  ndham  :number of energy bands (no core states)
Ci  ncore1 :number of core states for n (1st band index of zmelt)
Ci         :and offset to first valence state in zmelt(:,ncore1+1:,:)
Ci         :OLD: nctot
Ci  ncore2 :Number of core states for n' (1nd band index of zmelt) and should be eiher:
Ci         :     0 => do not calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :ncore1 => do calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :OLD: ncc
Ci  nbnb   :nbnb(k) = number of nonzero (it,jb) pairs for each k
Cio Inputs/Outputs
Cio n1b    :n1b(ix,k,jpm) = ib for ix^th nonzero (ib,jb) pair
Cio        :job=0 Input or job=1 Output :
Cio        : ib=1 = first valence state
Cio        : ib=ndham+1 = first core state
Cio        :job=1 Input and job=0 Output
Cio        : ib=1..ncore1  core states
Cio        : ib=ncore1+1 = first valence state
Cio n2b    :n2b(ix,k,jpm)= jb for ix^th nonzero (ib,jb) pair
Cio        : jb=1 = first valence state
Cio        : jb=ndham+1 = first core state
Cio        :job=1 Input and job=0 Output
Cio        : jb=1..ncore2  core states
Cio        : jb=ncore2+1 = first valence state
Cl Local variables
Cl  nfilo  :1st occ band state (1 corresponds to first band)
Cl  nfiup  :last occ band state
Cl  nemlo  :1st unocc band state (1 corresponds to first band)
Cl  nemup  :last unocc band state
Cr Remarks
Cr   Remarks on structure of n1b and n2b, normal case (jpm=1) and special case (jpm=2)
Cr   The table below shows the ordering of states
Cr      array  jpm  type       Ordering,core follows val          Ordering,core precedes val
Cr       n1b    1   occ     nfilo:nfiup,ndham+1:ndham+ncore1     1:ncore1,nfilo+ncore1:nfiup+ncore1
Cr       n2b    1   unocc   nemlo:nemup,ndham+1:ndham+ncore2     1:ncore2,nemlo+ncore2:nemup+ncore2
Cr       n1b    2   unocc   nemlo:nemup,ndham+1:ndham+ncore1     1:ncore1,nemlo+ncore1:nemup+ncore1
Cr       n2b    2   occ     nfilo:nfiup,ndham+1:ndham+ncore2     1:ncore2,nfilo+ncore2:nfiup+ncore2
Cr   Note: in the last two lines ncore1 and ncore2 should be reversed.
Cr   However, if ncore2 = 0 when npm=1  and  ncore2=ncore1 when npm=2 (as this routine assumes)
Cr   the inconsistency may be safely ignored.
Cu Updates
Cu   14 Dec 14  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: job,nbnbx,nqbz,npm,ncore1,ncore2,kshft,ndham
      integer :: nbnb(nqbz,npm),n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm)
C ... Local parameters
      integer ibib,jpm,k,it,itp,k1,k2

      if (ncore1 == 0 .and. ncore2 == 0) return

      if (kshft == 0) then
        k1=1; k2=nqbz
      else
        k1=kshft; k2=kshft
      endif


      do  jpm  = 1,npm
        do  k = k1, k2
          do  ibib = 1, nbnb(k,jpm)

            if (job == 0) then

              if (n1b(ibib,k,jpm) > ndham) then
                it = n1b(ibib,k,jpm) - ndham !core, starting at 1
                if (jpm == 2) call rx('alignn12: core states in unoccupied channel (jpm=2)')
              else
                it = ncore1 + n1b(ibib,k,jpm) !valence, first state at it=ncore1+1
              endif
              n1b(ibib,k,jpm) = it

              if (n2b(ibib,k,jpm) > ndham) then
                if (jpm == 1) call rx('alignn12: core states in unoccupied channel (jpm=1)')
                itp = n2b(ibib,k,jpm) - ndham !core
              else
                itp = ncore2 + n2b(ibib,k,jpm) !valence, first state = ncore2+1, 1st unocc state = nkqmin
              endif
              n2b(ibib,k,jpm) = itp

            else

              if (n1b(ibib,k,jpm) <= ncore1) then
                it = n1b(ibib,k,jpm) + ndham !Now core starts at (last valence)+1
              else
                it = n1b(ibib,k,jpm) - ncore1 !valence, first state at it=1
              endif
              n1b(ibib,k,jpm) = it

              if (n2b(ibib,k,jpm) <= ncore2) then
                itp = n2b(ibib,k,jpm) + ndham !Now core starts at (last valence)+1
              else
                itp = n2b(ibib,k,jpm) - ncore2 !valence, first state at itp=1
              endif
              n2b(ibib,k,jpm) = itp

            endif
          enddo
        enddo
      enddo
      end


      subroutine pmkzxq(job,nqbz,nbnbx,npm,k,ncore1,ncore2,nbmx,nbcut,nbcut2,
     .  nkqmin,nkmax,nkqmax,nbnb,n1b,n2b,ihw,nhw,nw,ndibiw,nibiw,pibiw)
C- Count number of pairs in each frequency bin for imaginary part of response function
C ----------------------------------------------------------------------
Ci  job    :0 make nibiw only
Ci         :1 make nibiw and pibiw
Ci  nqbz   :number of k-points in the 1st BZ
Ci  nbnbx  :dimensions n1b,n2b,ihw,nhw
Ci         :Must be at least as large as the max value of nbnb
Ci  npm    :1 if time reversal symmetry, 2 if not
Ci  k      :index to current kpoint
Ci  ncore1 :number of core states for n (1st band index of zmelt)
Ci         :and offset to first valence state in zmelt(:,ncore1+1:,:)
Ci         :OLD: nctot
Ci  ncore2 :Number of core states for n' (1nd band index of zmelt) and should be eiher:
Ci         :0      => do not calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :ncore1 => do calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         :OLD: ncc
Ci  nbmx   :max number of bands (excluding core) to include in calculation
Ci         :(value of nband_chi0 in GWinput)
Ci  nbcut  :exclude contribution to chi0 from pairs with
Ci         :occ valence state <= nbcut AND unoccupied valence state > nbcut2
Ci         :(nbcut read as nbcutlow in GWinput)
Ci  nbcut2 :exclude contribution to chi0 from pairs with
Ci         :occ valence state <= nbcut AND unoccupied valence state > nbcut2
Ci         :(nbcut2 read as nbcutlowto in GWinput)
Ci  nbnb   :nbnb(k) = number of nonzero (it,jb) pairs for kth qp
Ci  n1b    :n1b(ix,k) = ib for ix^th nonzero (ib,jb) pair
Ci         :ib=1 = first valence state.
Ci         :When core states are included they follow the valence states
Ci  n2b    :n2b(ix,k)= jb for ix^th nonzero (ib,jb) pair
Ci         :jb=1 = first valence state.
Ci         :When core states are included they follow the valence states
Ci  ihw    :ihw(nbnb(k),k) = index to first histogram bin within
Ci                           (demin,demax) for a given nbnb pair and k
Ci  nhw    :nhw(nbnb(k),k) = number of histograms bins witin
Ci                           (demin,demax) for a given nbnb pair and k
Ci  nw     :number of frequencies where x0 is calculated.
Ci         :(formerly nwhis)
Ci  nibiw  :leading dimension of pibiw
Co Outputs
Co   nibiw :number of (occ,unocc) pairs in each bin
Co   pibiw :(job=1) pibiw(i,iw) points to ibib for ith (occ,unocc) pair
Cr Remarks
Cr   This is a help routine to count (occ,unocc) pairs for a given frequency.
Cr   This routine identifies all ibib pairs modifying rcxq in frequency bin iw.
Cr   Evaluating rcxq using I as an inner product greatly improves efficiency
Cu Updates
Cu   04 Jun 16 (DMT) redesign enabling inner product over product functions I
Cu   03 Jun 16 Adapted from mkzxq.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,intent(in) :: job,nkqmin,k,nbcut,nbcut2,nbmx,nbnbx,
     .  ncore2,ncore1,nkmax,nkqmax,npm,nqbz,nw,ndibiw
      integer,intent(in) :: nbnb(nqbz,npm),n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm),
     .  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm)
      integer,intent(out):: nibiw(nw,npm),pibiw(ndibiw,nw,npm)
C ... Local parameters
      integer ibib,jpm,it,itp,iw,ilo,ihi,n1top,n2top,nbmxc

      nibiw = 0

C --- jpm=1 if TR symmetry; Loop over jpm=1,2 if not ---
      do  jpm = 1, npm

C   --- Sum over (n, n') pairs for this k ---
        do  ibib = 1, nbnb(k,jpm) !---  ibib loop

          n1top = ncore1+nkmax
          n2top = ncore2+nkqmax
          nbmxc = ncore1+nbmx
          if (n1b(ibib,k,jpm) > n1top .or. n2b(ibib,k,jpm) > n2top) cycle
          if (jpm==1 .and. n2b(ibib,k,jpm)>nbmxc) cycle  ! Past highest unocc state
          if (jpm==2 .and. n1b(ibib,k,jpm)>nbmxc) cycle
          if (nbcut > 0) then ! Special handling if nbcut has been set
            if (jpm == 1 .and. n1b(ibib,k,jpm)<=nbcut .and. n2b(ibib,k,jpm)>nbcut2) cycle
            if (jpm == 2 .and. n2b(ibib,k,jpm)<=nbcut .and. n1b(ibib,k,jpm)>nbcut2) cycle
          endif
          it = n1b(ibib,k,jpm)   !valence, with first occ state = ncore1+1
          itp = n2b(ibib,k,jpm) - nkqmin + 1 !valence, with first unocc state = ncore2+1

C         print *, 'ibib',ibib,k,jpm,it,itp

          ilo = ihw(ibib,k,jpm)
          ihi = ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
          if (ihi > nw) call rx('mkzxq: histogram point out of range')

          do  iw = ilo, ihi
            nibiw(iw,jpm) = nibiw(iw,jpm) + 1
            if (job == 1) pibiw(nibiw(iw,jpm),iw,jpm) = ibib
          enddo                 ! iw
        enddo                   ! ibib
      enddo                     ! jpm

      end
