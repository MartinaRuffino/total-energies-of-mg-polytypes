      subroutine pgfset(s_spec,nbas,bas,plat,lorder,lrat,ips,
     .  vshft,pgfsl,pgfvl,pgord,pgfn,npl,npadl,npadr,pgplp)
C- Make arrays related to principal layers
C ----------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: name
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   plat  :principal lattice vectors
Ci   lorder:if .false., no final check is made on the ordering of sites
Ci         :within each PL.  Else, a check is made.  pgfset aborts if
Ci         :the sites are badly ordered.
Ci   lrat  :if T, rationalize pgfsl, numbering PL 0, 1, 2, ...
Ci   ips   :ips(ib) is species corresponding to site ib
Ci   vshft :(for printout only)
Ci   pfgsl :pgfsl(ib) is PL to which site ib corresponds
Ci   pgfvl :pgfsl(ib) is the PL potential group
Co Outputs
Co   npl   :number of principal layers
Co   npadl,npadr: number of padding sites L and R
Co   pgord :permutataion table that orders basis; see Remarks (dimensioned 2*nbas)
Co   pgfn:  projection normal to transverse vectors, increasing with PL.
Co   pgplp: PL-dependent indices.  pgfset returns:
Co          1: cumulative number of basis atoms in this PL=0 ... this PL
Co          2: index to PL with inequivalent potential
Co          (the following are not made here; see supgh)
Co          3: source (column) dimension of GF for this PL
Co          4: field (row) matrix dimension for this PL
Co          5: matrix dimension for this PL including i-waves
Co          6: offset to diagonal part of g
Co          PL are numbered 0...npl-1 here
Cio Inputs/Outputs
Cr Local variables
Cl   pgfn   :vector normal to plane Plat(1) x Plat(2)
Cr Remarks
Cr   The PL must be ordered so that their projection along plat(3)
Cr   increases with increasing PL, and the bulk left and right layers
Cr   have the smallest and largest projections.
Cr
Cr   Also sites should be ordered by increasing PL.
Cr   Then sites within one PL are contiguous.
Cr
Cr   pgfset does the following:
Cr     1. Finds the permutation table pgord that orders basis by increasing PL
Cr     2. Determines pgfn = vector normal to the plane P1xP2,
Cr        and the rightmost projection onto pgfn in the L layer
Cr     3. Shift all PL that are ordered improperly relative to layer L
Cr     4. Optionally check that the set of PL are properly ordered
Cu Updates
Cu   01 Sep 11 Begin migration to f90 structures
Cu   09 Aug 05 Better check of sites ordering within PL
Cu   20 Sep 04 Better printout when PL are badly ordered
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lorder,lrat
      integer nbas,npl,npadl,npadr,pgfsl(nbas),pgfvl(nbas),
     .  pgord(nbas,2),pgplp(6,0:*),ips(nbas)
      double precision vshft(-7:nbas)
C ... For structures
!      include 'structures.h'
      type(str_spec):: s_spec(*)
C local variables
      logical flag
      integer i,j,ib,il,jl,jb,ipr,iv,nv,ib1,isw,isum,ibs,
     .  ils,ilL,stdo,ibold,ilold
      double precision bas(3,nbas),pgfn(3),plat(3,3),
     .  xx(4),xl,xpl,xr,xpr,x0,ddot,xold,xls,xlmax
      equivalence (xx(1),xl), (xx(2),xpl), (xx(3),xpr), (xx(4),xr)
      character*80 outs, fmt*20, spid*8
      procedure(integer) :: nglob

      call getpr(ipr)
      stdo = nglob('stdo')
      if (ipr >= 20) print '(1x)'

C --- Sort basis by increasing PL, preserving order for constant PL ---
      call ivshel(1,nbas,pgfsl,pgord,.true.)
      do  i = 1, nbas
        pgord(i,1) = pgord(i,1)+1
      enddo
      j = 1
      jl = pgfsl(pgord(1,1))
      do  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        if (il /= jl) then
          call ivshel(1,i-j,pgord(j,1),pgord(1,2),.false.)
C         call awrit1('%32:1i',' ',80,stdo,pgord)
          j  = i
          jl = il
        endif
      enddo
      call ivshel(1,nbas+1-j,pgord(j,1),pgord(1,2),.false.)

C --- Vector pgfn, and xl=Rightmost proj. in L layer ---
      call cross(plat,plat(1,2),pgfn)
      x0 = ddot(3,pgfn,1,plat(1,3),1)
      if (x0 < 0) call dscal(3,-1d0,pgfn,1)
      x0 = ddot(3,pgfn,1,pgfn,1)
      call dscal(3,1/sqrt(x0),pgfn,1)
      xl = -9999d0
      ilL = pgfsl(pgord(1,1))
      do  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        if (il == ilL) then
          xl = max(xl,ddot(3,pgfn,1,bas(1,ib),1))
        endif
      enddo

C --- Add plat(3) to any layer not to the right of the L layer ---
C     call prmx('starting pos',bas,3,3,nbas)
      flag = .false.
      do  jl = ilL+1, pgfsl(pgord(nbas,1))
        xr =  9999d0
        do  i = 1, nbas
          ib = pgord(i,1)
          il = pgfsl(ib)
          if (il == jl) then
            if (ddot(3,pgfn,1,bas(1,ib),1) < xr) then
              xr = ddot(3,pgfn,1,bas(1,ib),1)
              jb = ib
            endif
          endif
        enddo

        if (xr < xl) then
C         First-time printout if needed
          if (.not. flag) then
            call info2(2,0,0,' PGFSET (warning): PL %i%-1j '//
     .        'extends beyond sites in other PL:'//
     .        '%N max P1xP2 for PL %i is %,6d ... add plat(3) for '//
     .        'these layers:%N'//'   PL  min P1xP2  ib   New P1xP2',
     .        ilL,xl)
            flag = .true.
          endif
          do  i = 1, nbas
            ib = pgord(i,1)
            il = pgfsl(ib)
            if (il == jl) then
              call daxpy(3,1d0,plat(1,3),1,bas(1,ib),1)
            endif
          enddo
          xold = xr
          xr = ddot(3,pgfn,1,bas(1,jb),1)
          if (ipr > 1) write(stdo,1) jl,xold,jb,xr
    1     format(i4,f12.6,i4,f12.6)
        endif
      enddo
C     call prmx('pos after shift bulk',bas,3,3,nbas)

C --- Check for increasing projection with increasing PL  ---
      flag = .false.
      ibs = pgord(1,1)
      ils = pgfsl(ibs)
      xold  = -1d10
      ibold = 0
      ilold = 0
      xlmax = -1d10
      jl = pgfsl(pgord(1,1))
      npl = 1
      pgplp(1,npl) = 0
      pgplp(2,npl) = pgfvl(pgord(1,1))
      do  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        xl = ddot(3,pgfn,1,bas(1,ib),1)
C       Stay within a layer ... update xlmax
        if (il == jl) then
          pgplp(1,npl) = pgplp(1,npl)+1
          xls = xlmax
          xlmax = max(xlmax,xl)
C         for error tracking only ... save site index that changes xlmax
          if (xlmax /= xls) ibs = ib
          if (xlmax /= xls) ils = il
C       New layer: update xold and jl
        elseif (il > jl) then
          xold = xlmax
          ibold = ibs
          ilold = ils
          xlmax = xl
          npl = npl+1
          pgplp(1,npl) = 1
          pgplp(2,npl) = pgfvl(ib)
          jl = il
        else
          call rx('bug in pgfset')
        endif
        if (lorder .and. xl < xold) then
          call info8(10,1,0,' PGFSET:  program will abort ...%N'
     .    //' site %i (in PL %i) has projection along P1 x P2 = %;4d%N'
     .    //' site %i (in PL %i) has projection along P1 x P2 = %;4d',
     .      ibold,ilold,xold,ib,il,xl,0,0)
          flag = .true.
        endif
        if (lrat) pgfsl(ib) = npl
      enddo
      npadl = pgplp(1,1)
      npadr = pgplp(1,npl)

C --- Printout ---
      if (ipr >= 30) then
        call awrit5(' PGFSET: %i principal layers, %i+%i+%i sites.'//
     .    '  Normal:%3:1;5d',' ',80,stdo,npl,nbas,npadl,npadr,pgfn)
        j = 18
        fmt = '(1x,a,20I4)'
        do  i = 0, npl-1,j
          if (i > 0) print *, ' '
          print fmt,'  PL |',(ib, ib = i, min(j+i-1,npl-1))
          print fmt,'size |',(pgplp(1,ib+1), ib = i, min(j+i-1,npl-1))
          print fmt,' PLV |',(pgplp(2,ib+1), ib = i, min(j+i-1,npl-1))
        enddo
        do  34  i = 2, nbas
   34   if (pgord(i,1) <= pgord(i-1,1)) goto 36
        goto 37
   36   if (lrat) print '('' Basis reordered by increasing PL:'')'
   37   continue
C       j = 18
C       fmt = '(1x,a,20I4)'
        if (nbas < 100) fmt = '(1x,a,26I3)'
        if (ipr >= 40) then
          print 2
    2     format(' New ib Old   Spec      PL  PLV',15x,'Pos',16x,
     .            'Pos.h',4x,'Vshift')
          il = -1
          do  jb = 1, nbas
            ib = pgord(jb,1)
            outs = ' '
C           call spacks(0,'spec name',sspec,spid,ips(ib),ips(ib))
            spid = s_spec(ips(ib))%name
            write(outs,3) jb, ib, spid, pgfsl(ib)-1,
     .        pgfvl(ib), (bas(ib1,ib), ib1=1,3),
     .        ddot(3,pgfn,1,bas(1,ib),1)
    3       format(i4,i6,4x,a8,2i4,1x,3f10.5,f10.5,f10.6)
            call awrit3('%?#n#%23p        ##%a%?#n#%;10,6D##',outs,
     .        len(outs),-stdo,isw(pgfsl(ib) == il),
     .        vshft(ib),vshft(ib))
            il = pgfsl(ib)
          enddo
        else
          do  jb = 1, nbas,j
            print fmt, '  PL',
     .      (pgfsl(pgord(ib,1)), ib = 1+jb-1, min(j+jb-1,nbas))
            print fmt,'Site',(pgord(ib,1),ib=1+jb-1,min(j+jb-1,nbas))
          enddo
        endif
      endif

      if (flag) call fexit(-1,111,' PGFSET: badly ordered sites',0)

C --- Convert pgplp(1) into accumulated number.  isum avoids compiler bug
      do  jl = npl, 1,-1
        pgplp(1,jl) = isum(jl,pgplp(1,1),6)
      enddo
      pgplp(1,npl+1) = pgplp(1,npl) + npadr
C     pgplp(1,-1) = 0
C      pgplp(1,npl+2) = pgplp(1,npl+1) + npadr
C      pgplp(1,npl+2) = pgplp(1,npl+1) + npadr
C      print '(a,20i3)','size |',(pgplp(1,jl), jl = 0, npl+1)
C      stop

C --- Map pgfvl(isite) into pgplp(2,PL) ---
      if (lrat) then
        il = 0
        iv = pgfvl(pgord(1,1))
C ...   Check that all sites with a given pgfsl have constant pgfvl
        do  jl = 1, nbas
          ib = pgord(jl,1)
C         print *, pgfsl(ib),il, pgfvl(ib),iv
          if (pgfsl(ib) == il .and. pgfvl(ib) /= iv)
     .      call fexit3(-1,111,' Exit -1 PGFSET: improper pgfvl(%i):'//
     .      ' expected %i but found %i',ib,iv,pgfvl(ib))
          il = pgfsl(ib)
          iv = pgfvl(ib)
        enddo
        il = -1
        nv = 0
        do  jl = 1, nbas
          ib = pgord(jl,1)
          if (pgfsl(ib) /= il) then
            nv = nv+1
C           pgplp(2,nv) = pgfvl(ib)
          endif
          il = pgfsl(ib)
        enddo
      endif

C     call prmx('end of pgfset, pos',bas,3,3,nbas)

      end
