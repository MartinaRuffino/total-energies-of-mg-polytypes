      subroutine baslst(lstyld,iopt,slst,ilast,ips,nbas,slabl,z,nopt,
     .  optlst,optarg,nlist,list)
C- Generates a list of sites from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   lstyld:default list style
Ci   iopt  :1s digit
Ci         : 1  print out list
Ci         :10s digit
Ci         : 0, do not sort list
Ci         : 1, sort list, paring duplicates
Ci   slst  :string bearing information for list
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   nbas  :size of basis
Ci   z     :Table of nuclear charges by species
Ci     ... baslst is designed to convert switches (substring of slst)
Ci         into numbers; not implemented
Ci   nopt  :number of special modifiers
Ci   optlst:string of nopt tokens for modifiers
Co Outputs
Co   nlist :number of elements in list
Co   list  :list of sites
Co   ilast :index to last char in string parsed by baslst
Co   optarg:numbers possibly generated from optional switches
Cr Remarks
Cr   See slist for styles and their syntax
Cb Bugs
Cb  nlist may exceed size of list. Solution:
Cb  Require nlist as input, giving size of allocated list
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   07 Apr 03 New arguments lstyld and ilast
Cu   13 Sep 01 Supersedes old optget.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lstyld,iopt,nbas,ips(*),nlist,list(*),nopt,ilast
      character slst*(*),optlst(nopt)*(*)
      character*8 slabl(*)
      double precision z(*),optarg(*)
C ... Dynamically allocated arrays
      integer, allocatable :: llist(:)
C ... Local parameters
      integer i,j,j1,j2,ls,m,lstyle,k1,k2,
     .  iv,nlists,mxint,nspec,ib,is,iclbsj
      procedure(integer) :: mkilsd,lgunit,parg
      character dc*1
C ... External calls
      external cplxdm,dpzero,pblch1,rx,s2sph,tcn,tcx,ymscop,ymtrns,yprm,
     .         zmscop,ztoyy

      nlist = 0
      if (slst == ' ') return
      ls = len(slst)
      j1 = 1
      dc = slst(j1:j1)
      j1 = j1+1

C ... Return here to resume parsing for arguments
C  40 continue
      call nwordg(slst,0,dc//' ',1,j1,j2)

C ... Parse special arguments (not implemented)
      if (slst(j2+1:j2+1) /= ' ')  then
        do  i = 1, nopt
          call word(optlst(i),1,k1,k2)
          m = j1-1
          j = parg(optlst(i)(k1:k2),4,slst,m,ls,dc,1,1,iv,optarg(i))
        enddo
      endif

C --- Parse style for specifiying sites list ---
      lstyle = lstyld
      m = j1-1
      i = parg('style=',2,slst,m,ls,dc,1,1,iv,lstyle)
      if (i > 0) then
        j1 = m+2
        call nwordg(slst,0,dc//' ',1,j1,j2)
      endif

C --- List of sites ---
      if (lstyle >= 1) then
        nspec = mxint(nbas,ips)
        allocate(llist(nspec))
        if (lstyle == 2) j2=j2+1
        call slist(lstyle,slst(j1:j2),slabl,z,nspec,nlists,llist)
        if (lstyle /= 2) j2=j2+1
        nlist = 0
        do  i = 1, nlists
          is = llist(i)
          do  j = 1, nbas
            ib = iclbsj(is,ips,-nbas,j)
            if (ib >= 0) then
              nlist = nlist+1
              list(nlist) = ib
            endif
          enddo
        enddo
        if (nlist == 0 .or. mod(iopt/10,10) == 0) goto 99
        call ishell(nlist,list)
C        j = 1
C        do  i = 2, nlist
C          if (list(i) > list(j) .and. list(i) <= nbas) then
C            j = j+1
C            list(j) = list(i)
C          endif
C        enddo
C        nlist = j
C        j2=j2+1
        deallocate(llist)
C     Special style: sites with z=0
      elseif (slst(j1:j1+1) == 'z ' .or.
     .        slst(j1:j1+1) == 'Z ') then
        nlist = 0
        do  ib = 1, nbas
          is = ips(ib)
          if (z(is) == 0) then
            nlist = nlist+1
            list(nlist) = ib
          endif
        enddo

        j2=j2+1

      else
        j = j1
        call nwordg(slst,0,dc//' ',1,j,j2)
C       call mkilst(slst(j1:j2),nlist,list)
        nlist = -1
        nlist = mkilsd(slst(j1:j2),-1,list)
        if (nlist < 0) call rxs('baslst: failed to parse site list ',slst(j1:))
        nlist = mkilsd(slst(j1:j2),nlist,list)
        if (nlist == 0 .or. mod(iopt/10,10) == 0) goto 99
        call ishell(nlist,list)
        j = 1
        do  i = 2, nlist
          if (list(i) > list(j) .and. list(i) <= nbas) then
            j = j+1
            list(j) = list(i)
          endif
        enddo
        nlist = j
        j2=j2+1
      endif

      ilast = j2 ! parsing ended at ilast
   99 if (mod(iopt,10) == 0) return
      call awrit2(' baslst: %n:1i',slst,255,lgunit(1),nlist,list)
      end
