      subroutine clsset(mode,s_lat,nbas,pos,ncl,clssl,cllst,clord,clp)
C- Groups all sites in a crystal into clusters
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit determines how clssl is made.
Ci           0  clssl is input
Ci           1  divide lattice into microcells of size plate
Ci          10s digit set rationalizes cluster index numbering
Ci   slat  :struct containing information about the lattice
Ci   nbas  :size of basis (input)
Ci   pos   :basis vectors (input)
Co Outputs
Co   ncl   :number of clusters
Co   clssl :cluster index: site ib belongs to cluster clssl(ib)
Co   cllst :list of sites belonging to each cluster.
Co          The range of sites is clp(1,ic)+1 ... clp(1,ic+1)
Co   clord :permutation table of indices that sorts basis by cluster
Co   clp   :CL-dependent indices.  clsset returns for each cluster:
Co          1: cumulative number of basis atoms in CL=1 ... this CL
Co          (the following are not made here)
Co          2: not used now
Co          3: lower matrix dimension for diagonal part of this CL
Co          4: downfolded matrix dimension for diagonal part of this CL
Co          5: matrix dimension for CL, incl. off-diagonal
Co          6: not used now
Cr Remarks
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nclp
      parameter (nclp=9)
      integer mode,nbas,clssl(nbas),cllst(nbas),clord(nbas),clp(nclp,*)
      double precision pos(3,nbas)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      real(8), allocatable :: posl(:),pose(:),clctr(:)
      integer, allocatable :: ncltb(:),iwk(:)
C ... Local parameters
      logical lrat
      integer lgunit,m,k,i,ib,ipr,ipose(3),i1mach,ndig,ndigi,
     .  ifac,ncl,idmx,ngmx,ngmn,ncol,isw,ix(3),j,i1,i2
      double precision plat(3,3),plate(3,3),qlat(3,3),qlate(3,3)
      double precision vol,volm,xx
      character*80 outs,fmt

C ... Setup
      plat = s_lat%plat
      lrat = mod(mode/10,10) /= 0
      call getpr(ipr)

C --- Divide total cell into microcells ---
      if (mod(mode,10) == 0) then
      elseif (mod(mode,10) == 1) then
        plate = s_lat%plate

C   ... Determine supercell size
        call dinv33(plat,1,qlat,vol)
        call dinv33(plate,1,qlate,volm)

        print 1
    1   format(/15x,'Plat',31x,'Mlat')
        print 2,((plat(m,k),m=1,3),(plate(m,k),m=1,3),k=1,3)
    2   format(3F10.5,5x,3F10.5)
        call awrit3('   Cell vol: %,3d   Cluster vol: %;3d'//
     .    '  ratio: %;3g',
     .    ' ',80,lgunit(1),abs(vol),abs(volm),abs(vol/volm))

        ndig = int(dlog(dble(i1mach(9)))/dlog(2d0)/3)
        ifac = 2**ndig
        if (ifac == 1024) ifac = 1000
        allocate(pose(3*nbas))
        ix(1) = 1
        ix(2) = 1
        ix(3) = 1
        call shorps(nbas,plat,ix,pos,pos)
c       call prmx('shortened pos',pos,3,3,nbas)
        allocate(posl(3*nbas))
        call dcopy(3*nbas,pos,1,posl,1)
C   ... Shift pos off axis by a tiny amount, avoiding boundary problems
*       call daxpy(3*nbas,1d0,1d-5,0,posl,1)
c       call prmx('plat',plat,3,3,3)
c       call prmx('plate',plate,3,3,3)
C   ... posl = positions as multiples of plat
        call prjpos(2,ix,xx,plat,nbas,pos,posl)
C       call prmx('pos, units of plat',posl,3,3,nbas)
C   ... Positions as multiples of plate
        call prjpos(63,ix,plat,plate,nbas,posl,pose)
C       call prmx('pos, units of plate',pose,3,3,nbas)
C   ... index as ifac**2*pose(1) + ifac*pose(2) + pose(3), integer parts
        do  ib = 1, nbas
          call discop(pose,ipose,3,3*ib-2,1,0)
          clssl(ib) = ipose(1)*ifac*ifac + ipose(2)*ifac + ipose(3)
        enddo
        deallocate(posl)

C   ... Attempt to combine clusters
C  30   continue
        call ngroup(nbas,clssl,clord,ncl,ngmn,ngmx)
C        if (ncl < nbas) then
C          if (pvclst(nbas,pose,clssl,clord)) goto 30
C        endif
        deallocate(pose)
      else
        call rx1('clsset: unrecognized mode %i',mode)
      endif

C --- Count and sort cluster list ---
      call ngroup(nbas,clssl,clord,ncl,ngmn,ngmx)

C --- Rationalize clusters and printout
C     call awrit2(' clsset: %n:1i',' ',255,lgunit(1),nbas,clssl)
      idmx = clssl(clord(nbas))
      ndigi = max(int(dlog10(dble(max(idmx,1))))+2,5)
      ndig = int(dlog10(dble(nbas))) + 2
      if (ipr >= 20)
     . call awrit7('%N CLSSET: %i clusters.'//
     .    '%?#n#%b (id rationalized).##'//
     .    '  Min, max, average sizes= %i, %i, %;1d'//
     .    '%?#n>=30#%N Cluster%(n-4)fid   sites ...',
     .    ' ',120,lgunit(1),ncl,isw(lrat),ngmn,ngmx,dble(nbas)/ncl,
     .    ipr,ndigi)
      i2 = 0
      outs = ' '
      clp(1,1) = 0
      do  i = 1, ncl
        i1 = i2+1
        call nxtgrp(nbas,clssl(clord(i1)),clssl,clord,i1,i2)
        clp(1,i+1) = clp(1,i) + i2-i1+1
        if (lrat) clssl(clord(i1)) = i
        if (ipr >= 30) call awrit4('  %,ni %,ni:',outs,len(outs),
     .    0,ndig,i,ndigi,clssl(clord(i1)))
        k = 1
        do  j = i1, i2
          if (lrat) clssl(clord(j)) = i
          if (ipr >= 30)
     .      call awrit2('%a%,ni',outs,len(outs),0,ndig+k,clord(j))
          k=0
        enddo
        call awrit0(outs,' ',-80,lgunit(1))
      enddo

C ... Make cllst
      allocate(ncltb(ncl+1),clctr(3*nbas))
      call clsctr(nbas,pos,clssl,clord,ncltb,cllst,clctr)
      deallocate(ncltb,clctr)

C ... Print re-ordering of sites
C     call awrit2(' clord: %n:1i',' ',255,lgunit(1),nbas,clord)
C     call awrit2(' clssl: %n:1i',' ',255,lgunit(1),nbas,clssl)
C     call defi(oiwk, nbas)
      allocate(iwk(nbas))
      call ivprm(1,nbas,clssl,iwk,clord,0)
      call awrit1(' ... sites reordered by cluster%?#n>30#:',' ',80,
     .  lgunit(1),ipr)
      if (ipr > 30) then
        idmx = clssl(clord(nbas))
        ndigi = max(int(dlog10(dble(idmx)))+2,5)
        call awrit3(' %%,%ii%%,%ii%%,%ii',fmt,len(fmt),0,ndig+1,ndig+1,
     .    ndigi+1)
        call word(fmt,1,i,k)
        ncol = 84/(2*ndig+ndigi+7)
        outs = ' Old New   id'
        call arrprt(outs,fmt(i:k),'iIi',nbas,0,ncol,
     .    0,'  | ',clord,xx,iwk,xx,xx,xx,xx,xx)
      endif

C ... Re-order clssl
      call icopy(nbas,iwk,1,clssl,1)
      deallocate(iwk)

C     call awrit2(' clssl: %n:1i',' ',255,lgunit(1),nbas,clssl)

      end
      logical function pvclst(nbas,pose,clssl,clord)
C- Attempt to combine clusters
      implicit none
      integer nbas,clssl(nbas),clord(nbas)
      double precision pose(3,nbas)
      integer i1,i2,i,il,ib,ic,lgunit
      integer j1,j2,j,jl,jb,jc
      double precision tol
      parameter (tol=1d-5)

      pvclst = .false.
   10 continue
      i2 = 0
      do  il = 1, nbas
        i1 = i2+1
        call nxtgrp(nbas,clssl(clord(i1)),clssl,clord,i1,i2)
        if (i1 > i2) return

        j2 = i2
        do  jl = il, nbas
          j1 = j2+1
          call nxtgrp(nbas,clssl(clord(j1)),clssl,clord,j1,j2)
          if (j1 <= j2) then
C     ... Check whether every member in i close to every member in j
            do  i = i1, i2
            ib = clord(i)
              do  j = j1, j2
              jb = clord(j)
              if (abs(pose(1,ib)-pose(1,jb)) > 1+tol .or.
     .            abs(pose(2,ib)-pose(2,jb)) > 1+tol .or.
     .              abs(pose(3,ib)-pose(3,jb)) > 1+tol) goto 15
              enddo
            enddo
          ic = clssl(clord(i1))
          jc = clssl(clord(j1))
          pvclst = .true.
C     ... Merge groups i and j
          call awrit2(' merge id=%i with id=%i',' ',80,lgunit(1),ic,jc)
            do  i = i1, i2
              clssl(clord(i)) = jc
            enddo
C     ... resort the group
          call ivshel(1,nbas,clssl,clord,.true.)
            do  i = 1, nbas
              clord(i) = clord(i)+1
            enddo
            goto 10
          endif
 15     enddo

      enddo

      call rx('bug in pvclst')

      end
      subroutine clsctr(nbas,bas,clssl,clord,ncltab,cllst,ctr)
C- Makes cluster list and center of gravity
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   bas   :basis vectors, in units of alat (input)
Ci   clssl :cluster index: site ib belongs to cluster clssl(ib) (clsset.f)
Co Outputs
Co   ncltab:ncltab(ic)=# sites in clsst preceding cluster ic
Co   cllst :cllst(ncltab(ic+1)..ncltab(ic+1)) are the list of
Co          sites in cluster ic.
Co   ctr   :positions of cluster centers, for each site
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,clssl(nbas),clord(nbas),ncltab(*),cllst(nbas)
      double precision bas(3,nbas),ctr(3,nbas)
C ... Local parameters
      integer ncl,ib,kb,icl,kcl,ncli,ix,lgunit,ipr
      double precision cg

C --- Make cluster list cllst ---
      icl = clssl(clord(1))
      ncl = 0
      ncltab(1) = 0
      do  ib = 1, nbas
        kb = clord(ib)
        kcl = clssl(kb)
        cllst(ib) = kb
        if (kcl /= icl) then
          ncl = ncl+1
          ncltab(ncl+1) = ib-1
          icl = kcl
        endif
      enddo
      ncl = ncl+1
      ncltab(ncl+1) = nbas

C --- Make cluster center of gravity ---
      do  ix = 1, 3
        do  icl = 1, ncl
        cg = 0
        ncli = ncltab(icl+1) - ncltab(icl)
        cg = 0
          do  ib = ncltab(icl)+1, ncltab(icl+1)
          kb = cllst(ib)
          cg = cg + bas(ix,kb)
          enddo
        cg = cg/ncli
          do  ib = ncltab(icl)+1, ncltab(icl+1)
          kb = cllst(ib)
          ctr(ix,kb) = cg
          enddo
        enddo
      enddo

C --- Printout ---
      call getpr(ipr)
      if (ipr < 40) return
      if (ipr > 40) then
        print '(1x)'
        print 1
    1   format(/'  site',9x,'------- pos --------',14x,
     .    '------- pos-cg --------')
      endif
      do  icl = 1, ncl
        ncli = ncltab(icl+1) - ncltab(icl)
c       print 332, icl, ncli, (ctr(ix,kb), ix=1,3)
        kb = cllst(ncltab(icl)+1)
        call awrit3(' cluster%,3i  has %i elements.  %33pcenter ='//
     .    ' %3;11,6D',' ',80,lgunit(1),icl,ncli,ctr(1,kb))
        if (ipr > 40) then
          do  ib = ncltab(icl)+1, ncltab(icl+1)
          kb = cllst(ib)
            print 2,kb,(bas(ix,kb),ix=1,3),
     .                   (bas(ix,kb)-ctr(ix,kb), ix=1,3)
    2       format(i4,2x,3F11.6,3x,3F11.6)
          enddo
        endif
      enddo

      end
