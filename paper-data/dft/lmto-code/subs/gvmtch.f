      subroutine gvmtch(ngs,gvs,kvs,ngb,gvb,kvb)
C- Align block of big list of gvecs to small list, purging mismatches
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngs   :number of G vectors in small list
Ci   gvs   :small list of G vectors. gvs must be sorted by length
Ci   kvs   :indices for gather/scatter operations in gvs; see gvlist.f
Ci   ngb   :number of G vectors in big list.  ngb must be >= ngs
Cio Inputs/Outputs
Cio  gvb   :big list of G vectors. gvb must be sorted by length
Cio        :All vectors in gvs should be contained somewhere in gvb.
Cio        :On output, gvb is permuted so that first ngs vectors
Cio        :match gvs.  Only gvb(1:ngs,:) has meaning.
Cio  kvb   :indices for gather/scatter operations in gvb; see gvlist.f
Cio        :On output, kvb is permuted in the same manner as gvs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngs,ngb
      integer kvs(ngs,3),kvb(ngb,3)
      double precision gvs(ngs,3),gvb(ngb,3)
C ... Dynamically allocated local arrays
      real(8), allocatable :: gg(:)
      real(8), allocatable :: kk(:)
      real(8), allocatable :: iwk(:)
C ... Local parameters
      integer ig,low,high,jg,mm
      double precision xx,tol,tol2
      parameter (tol=1d-6,tol2=1d-9)
C ... External calls
      external huntx,rxi

      allocate(gg(ngs),kk(ngs),iwk(ngs))

C      call prm('gvs',gvs,ngs,ngs,3)
C      call prm('gvb',gvb,ngb,ngb,3)

C ... Generate length of g and ensure length matches for both
      do  ig = 1, ngs
        gg(ig) = gvs(ig,1)**2+gvs(ig,2)**2+gvs(ig,3)**2
      enddo

C     call prm('gg',gg,ngs,ngs,1)

C ... Sanity check.  Make sure |G| is same for each vector in both lists
C      do  20  ig = 1, ngs
C        xx = gvb(ig,1)**2 + gvb(ig,2)**2 + gvb(ig,3)**2
C        if (abs(xx-gg(ig)) > tol) then
C          print *, xx,gg(ig)
C          call rxi('bug in gvmatch, ig=',ig)
C        endif
C   20 continue

C --- For each G in small, find matching G in big list ---
      xx = -1d0
      low = 0
      high = 0
      do  ig = 1, ngs
        iwk(ig) = -1
C   ... Find first and last g-vector list with same length
        if (abs(xx-gg(ig)) > tol) then
          call huntx(gg,ngs,gg(ig)+tol,0,high)
          low = ig
          high = min(high,ngs)
          xx = gg(ig)
        endif

        do  jg = low, high
          do  mm = 1, 3
            if (abs(gvb(jg,mm)-gvs(ig,mm)) > tol2) goto 5
          enddo
C     ... Found a match
          iwk(ig) = jg
    5   enddo
C   ... Sanity check
        if (iwk(ig) == -1) call rxi('bug in gvmatch, ig=',ig)
      enddo

C ... Rearrange gvb, kvb
      do  mm = 1, 3
        do  ig = 1, ngs
          jg = iwk(ig)
          gg(ig) = gvb(jg,mm)
          kk(ig) = kvb(jg,mm)
        enddo
        do  ig = 1, ngs
          gvb(ig,mm) = gg(ig)
          kvb(ig,mm) = kk(ig)
        enddo
      enddo

c      call prm('gvb',gvb,ngb,ngb,3)

      deallocate(gg,kk,iwk)
      end

      subroutine gvmtch2(ngd,ngs,gvs,ngb,gvb,kvb)
C- Align block of big list of gvecs to small list, purging mismatches
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngd   :leading dimension of gvs,gvb
Ci   ngs   :number of G vectors in small list
Ci   gvs   :small list of G vectors. gvs must be sorted by length
Ci   ngb   :number of G vectors in big list.  ngb must be >= ngs
Cio Inputs/Outputs
Cio  gvb   :big list of G vectors. gvb must be sorted by length
Cio        :All vectors in gvs should be contained somewhere in gvb.
Cio        :On output, gvb is permuted so that first ngs vectors
Cio        :match gvs.  Only gvb(1:ngs,:) has meaning.
Cio  kvb   :indices for gather/scatter operations in gvb; see gvlist.f
Cio        :On output, kvb is permuted in the same manner as gvs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   Routine finds gvb(jg,:) that match gvs(ig,:) and
Cr   reorders gvb to so the first ngs vectors match gvs.
Cr   The reordering is also applied to kvb
Cu Updates
Cu   25 Jun 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngd,ngs,ngb
      integer kvb(ngd,3)
      double precision gvs(ngd,3),gvb(ngd,3)
C ... Dynamically allocated local arrays
      real(8), allocatable :: ggs(:),ggb(:)
      real(8), allocatable :: kk(:)
      real(8), allocatable :: iwk(:)
C ... Local parameters
      integer ig,low,high,jg,mm
      double precision xx,tol,tol2,dlength
      parameter (tol=1d-6,tol2=1d-9)
C ... External calls
      external huntx,rxi

      if (ngb < ngs) call rx('gvmtch2: ngb < ngs')

      allocate(ggs(ngs),ggb(ngb),kk(ngs),iwk(ngs))

C      call prm('gvs',gvs,ngs,ngs,3)
C      call prm('gvb',gvb,ngb,ngb,3)

C ... Generate length of gvs, gvb
      do  ig = 1, ngs
        ggs(ig) = gvs(ig,1)**2+gvs(ig,2)**2+gvs(ig,3)**2
      enddo
      do  ig = 1, ngb
        ggb(ig) = gvb(ig,1)**2+gvb(ig,2)**2+gvb(ig,3)**2
      enddo

C     call prm('ggs',ggs,ngs,ngs,1)

C ... Sanity check.  Make sure |G| is same for each vector in both lists
C      do  20  ig = 1, ngs
C        xx = gvb(ig,1)**2 + gvb(ig,2)**2 + gvb(ig,3)**2
C        if (abs(xx-ggs(ig)) > tol) then
C          print *, xx,ggs(ig)
C          call rxi('bug in gvmatch, ig=',ig)
C        endif
C   20 continue

C --- For each G in small, find matching G in big list ---
      xx = -1d0
      low = 0
      high = 0
      do  ig = 1, ngs
        iwk(ig) = -1
C   ... Find first and last g-vector in big list with matching length
        if (abs(xx-ggs(ig)) > tol) then
          call huntx(ggb,ngb,ggs(ig)-tol,0,low)
          low = max(min(low,ngb),1)
          high = low
          call huntx(ggb,ngb,ggs(ig)+tol,0,high)
          high = min(high,ngb)
          xx = ggs(ig)
        endif

        do  jg = low, high
          if (dlength(3,gvb(jg,1:3)-gvs(ig,1:3),1) < tol) then
            iwk(ig) = jg
            exit
          endif
        enddo
        if (iwk(ig) == -1) call rxi('bug in gvmtch2, ig=',ig)
    5   continue
      enddo

C ... Rearrange gvb, kvb
      do  mm = 1, 3
        do  ig = 1, ngs
          jg = iwk(ig)
          ggs(ig) = gvb(jg,mm)
          kk(ig) = kvb(jg,mm)
        enddo
        do  ig = 1, ngs
          gvb(ig,mm) = ggs(ig)
          kvb(ig,mm) = kk(ig)
        enddo
      enddo

C      call prmx('gvs',gvs,ngd,ngs,3)
C      call prmx('gvb',gvb,ngd,ngs,3)

      deallocate(ggs,kk,iwk)
      end

