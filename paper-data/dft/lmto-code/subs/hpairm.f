      subroutine hpairm(nbas,ips,alat,plat,pos,rham,nttab,ntab,iax,rtab,
     .                 mxcsiz,ordrd,niax)
C- Neighbor table for Hamiltonian matrix (molecules version)
C ----------------------------------------------------------------------
C  Inputs:  nbas,ips,alat,plat,pos,rham
C  Outputs:
C    nttab: total number of pairs within range dr <= rham(i)+rham(j)
C    iax(1:2,*) sites ib,jb of each pair
C    iax(3:5,*) multiples of plat connecting sites (zero for molecules)
C    iax(6,*)   points to conjugate (jb,ib) pair
C    iax(7,*)   permutation table ordering cluster by increasing ib
C    ntab(1:nbas+1):  offset to first element in iax and rtab for ib
C    rtab(1:3,*): connecting vector defined as pos(jb)-pos(ib)
C    mxcsiz: number of atoms in largest cluster (for dimensioning)
C    ordrd: 0 order each cluster by increasing jb; 1 by increasing dist.
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nttab,ips(nbas),ntab(nbas+1),mxcsiz,ordrd,niax
      double precision alat,plat(3,3),pos(3,nbas),rham(nbas),rtab(3,*)
C#ifndef LINUXF
      integer, allocatable, intent(inout) :: iax(:)
C#elseC
C      integer, pointer :: iax(:)
C      include 'interfaces77.h'
C#endif
C ... Dynamically allocated local arrays
      integer, allocatable :: wk1(:),wk2(:),iaxl(:,:)
      real(8), allocatable :: wk3(:)
C ... Local parameters
      double precision r1,r2,rr
      integer ib,is,jb,js,ipr,nt,mtab,i,it,jt,stdo,nglob

      call getpr(ipr)
      stdo = nglob('stdo')
      if (ipr > 40) print 1
    1 format(/' hpair:     ib   jb',9x,'------- dr --------',11x,
     .        'd      r1+r2')

      allocate(iaxl(niax,300*nbas))
C --- Accumulate all pairs for which dr < ri + rj ---
      nttab = 1
      ntab(1) = 0
      mtab = 1
      mxcsiz = 0
      do  ib = 1, nbas
      is = ips(ib)
      r1 = rham(is)
        do  jb = 1, nbas
         js = ips(jb)
         r2 = rham(js)
C ...   Add to list if connecting vector < r1+r2
         rtab(1:3,nttab) = alat*(pos(1:3,jb)-pos(1:3,ib))
         rr = sum(rtab(1:3,nttab)*rtab(1:3,nttab))
         if (rr < (r1+r2)*(r1+r2)) then
            iaxl(1,nttab) = ib
            iaxl(2,nttab) = jb
            iaxl(3,nttab) = 0
            iaxl(4,nttab) = 0
            iaxl(5,nttab) = 0
            nttab = nttab+1
         endif
        enddo

C --- Sort table by increasing length ---
        if (ordrd == 1) then
          allocate(wk1((nttab-mtab)*6))
          allocate(wk2(nttab-mtab))
          allocate(wk3((nttab-mtab)*3))
          call xhpair(nttab-mtab,iaxl(1,mtab),rtab(1,mtab),wk1,wk2,wk3)
          deallocate(wk1,wk2,wk3)
        endif

C --- Printout for this ib ---
        if (ipr >= 30) print 2,ib,nttab-mtab
    2   format(' hpair, ib=',i3,':',i4,' neighbors')
        if (ipr > 40) then
          do  nt = mtab, nttab-1
            rr = sqrt(sum(rtab(1:3,nt)*rtab(1:3,nt)))
            r1 = rham(ips(iaxl(1,nt)))
            r2 = rham(ips(iaxl(2,nt)))
C    .    ,iax(3,nt),iax(4,nt),iax(5,nt)
            write(stdo,3) iaxl(1,nt),iaxl(2,nt),
     .        (rtab(i,nt),i=1,3),rr,r1+r2,iaxl(7,nt)
    3       format(i14,i5,3F11.6,2F9.4,2x,3I3)
          enddo
        endif

C --- Cleanup for this ib ---
        mtab = nttab
        ntab(ib+1) = nttab-1
        mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
      enddo
      nttab = nttab-1
C     call awrit2('xx ntab %n:1i',' ',80,6,nbas+1,ntab)
      if (ipr >= 20 .and. ipr <= 40) print
     .  '('' hpair:'',i5,'' pairs total'',i6,'' is max cluster size'')',
     .  nttab, mxcsiz
C --- For each pair, find matching pair ---
      do  it = 1, nttab
        iaxl(6,it) = 0
      enddo
      do  ib = 1, nbas
      do  it = ntab(ib)+1, ntab(ib+1)
        if (iaxl(6,it) /= 0) cycle
        jb = iaxl(2,it)
        do  jt = ntab(jb)+1, ntab(jb+1)
          if (iaxl(2,jt) == ib .and.
     .        iaxl(3,it) == -iaxl(3,jt) .and.
     .        iaxl(4,it) == -iaxl(4,jt) .and.
     .        iaxl(5,it) == -iaxl(5,jt))  then
            iaxl(6,it) = jt
            iaxl(6,jt) = it
            goto 5
          endif
        enddo
        call rx('bug in hpair')
    5   continue
      enddo
      enddo

C#ifndef LINUXF
      if (allocated(iax)) deallocate(iax)
      allocate(iax(niax*nttab))
C#elseC
C      if (associated(iax)) deallocate(iax)
C      call i1alloc(iax,niax*nttab,0)
C#endif
      call icopy(niax*nttab, iaxl, 1, iax, 1)
      deallocate(iaxl)

C     do  80  ib = 1, nbas
C  80 print '(99i3)', (iax(6,i), i=ntab(ib)+1, ntab(ib+1))
c     call rx('hpairm')
      end
      subroutine xhpair(nttab,iax,rtab,iwk,iwk2,rwk)
C- Sort neighbor table by distance
      implicit none
      integer nttab,iax(10,nttab),iwk(6,nttab),iwk2(nttab),i,j,k
      double precision rtab(3,nttab),rwk(3,nttab)

      do  i = 1, nttab
        rwk(1,i) = rtab(1,i)
        rwk(2,i) = rtab(2,i)
        rwk(3,i) = rtab(3,i)
        do  k = 1, 6
          iwk(k,i) = iax(k,i)
        enddo
      enddo
      call dvshel(3,nttab,rtab,iwk2,11)
      do  i = 1, nttab
        j = iwk2(i)+1
        rtab(1,i) = rwk(1,j)
        rtab(2,i) = rwk(2,j)
        rtab(3,i) = rwk(3,j)
        do  k = 1, 6
          iax(k,i) = iwk(k,j)
        enddo
        iax(7,j) = i
      enddo
      end
      subroutine pairmc(kb,nclus3,iax,iax3,ntab,rtab,ipsite)
C- Table of pointers to pair table for cluster around a third site
C  Output: ipsite(ic,jc)
C  Three-center terms arising from augmentation site kb contribute to
C  the hamiltonian for each pair (ic,jc) that connects to kc.
C  is=ipsite(ic,jc) locates the iax(*,is) corresponding to (ic,jc).
C
C  For efficient implementation, this routine relies on a sorted iax
C  table.  iax(7,*) is a table of permutations ip(*) that sorts each
C  local cluster by basis atom ib (see hpairm).  A sorted
C  iax table guarantees that ib(ip(i+1)) > ib(ip(i)); therefore
C  ipsite(ib(ip(ic)),jb(ip(jc+1))) > ipsite(ib,jb(ip(jc))).
      implicit none
      integer kb,nclus3,iax(10,*),iax3(nclus3),ntab(*),
     .  ipsite(nclus3,nclus3)
      double precision rtab(3,1)
      integer kc,i,j,ii,jj,ip,jp,ib,jb,nc,ic,icp,ic0,is,nclus,ipr
      integer ipx(200),ipxi(200)

      call getpr(ipr)
      if (nclus3 > 200) call rx('pairmc: increase dim of ipx')
      do  ii = 1, nclus3**2
        ipsite(ii,1) = 0
      enddo
      kc = 1+ntab(kb)

C --- Make permutation table for all elts in 3C range ---
      nclus = ntab(kb+1)-ntab(kb)
      i = 0
      do  ii = 1, nclus
        ip = iax(7,kc-1+ii)
C   ... if ip is a member of iax(3), keep it; else go to next
        do  j = 1, nclus3
          if (ip == iax3(j)) then
            i = i+1
            ipx(i) = ip
            ipxi(i) = j
            exit
          endif
        enddo
      enddo

C      print 345, (iax3(i), i=1,nclus3)
C      print 345, (iax(7,kc-1+i), i=1,nclus)
C      print 345, (ipx(i), i=1,nclus3)
C      print 345, (ipxi(i), i=1,nclus3)
C      print 345, (iax(2,kc-1+ipx(i)), i=1,nclus3)
C  345 format(34i4)

C --- Loop over 3C pairs i,j in permuted order ---
      do  ii = 1, nclus3
        ip = ipx(ii)
        ib = iax(2,kc-1+ip)
C   ... ic is offset corresponding to iax(1,ic) = 1st pair for ib
        ic0 = ntab(ib)
        ic  = ic0
        nc  = ntab(ib+1)
        do  jj = 1, nclus3
          jp = ipx(jj)
C     ... Do upper triangle only
          if (jp < ip) cycle
          jb = iax(2,kc-1+jp)
C     ... Increment ic until this (ib,jb) pair is found
   42     ic = ic+1
          icp = ic0+iax(7,ic)
C     ... Skip if missing pair corresponding to this 3C (ib,jb)
          if (ic > nc) cycle
C          call awrit5('want ib,jb=%i %i, site %i has %i %i',
C     .      ' ',80,6,ib,jb,icp,ib,iax(2,icp))
C          if (ic > nc .or. icp < 0) then
C            call awrit2('no pair ib,jb=%i %i ...',' ',80,6,ib,jb)
C            pause
C            cycle
C          endif
          if (jb /= iax(2,icp)) goto 42
C          call awrit7('found ip,jp=%i %i  ib,jb=%i %i ic,jc= %i %i'//
C     .      ' site %i',' ',80,6,ip,jp,ib,jb,ipxi(ii),ipxi(jj),icp)
          ipsite(ipxi(ii),ipxi(jj)) = icp
        enddo
      enddo

C --- Construct lower triangle of table ---
      do  i = 1, nclus3
        do  j = i+1, nclus3
          is = ipsite(i,j)
          if (is > 0) ipsite(j,i) = iax(6,is)
        enddo
      enddo

C --- Check ---
      if (ipr < 70) return
C      print '(''% rows'',i4,'' cols'',i4)', nclus3, nclus3
C      do  60  i = 1, nclus3
C   60 print '(35i4)', (ipsite(i,j), j=1, nclus3)
      do  ii = 1, nclus3
      do  jj = 1, nclus3
        i = iax3(ii)
        j = iax3(jj)
        ib = iax(2,kc-1+i)
        jb = iax(2,kc-1+j)
        is = ipsite(ii,jj)
        if (is == 0) then
          print 356, kb,ib,jb
  356     format(' pairmc, kb=',i3,':  missing 2C pair for ib,jb=',2i4)
        elseif (iax(1,is) /= ib .or. iax(2,is) /= jb) then
          print *, 'problem with i,j,is=',i,j,ib,jb,is
          call rx('bug in pairmc')
        endif
      enddo
      enddo
      end
C Testing ...
C      subroutine fmain
C      use structures, only : hpairm
C      implicit none
C      integer nbas
C      parameter (nbas=3)
C      double precision alat,plat(3,3),pos(3,nbas),rham(nbas)
C      integer ips(nbas)
C      integer mxcsiz,nttab
CC     double precision w(10000)
C      integer, allocatable, target :: iax(:)
C      integer ntab(nbas+1)
C      double precision rtab(3,300*nbas),xx,dglob
C
C      data alat/1d0/ plat/10d0,0d0,0d0,0d0,10d0,0d0,0d0,0d0,10d0/,
C     .  pos/0d0,0d0,0d0,1.3804348d0,0d0,1.05924488d0,
C     .  -1.3804348d0,0d0,1.05924488d0/,
C     .  rham/2d0,3d0,0d0/, ips/2,1,1/
C
C      call pshpr(100)
C      xx = dglob('stdo',6d0,1)
C
CC --- Obtain neighbor table for hamiltonian ---
C      call hpairm(nbas,ips,alat,plat,pos,rham,nttab,ntab,
C     .  iax,rtab,mxcsiz,1,10)
C
C      end
