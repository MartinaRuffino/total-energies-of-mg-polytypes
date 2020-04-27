      subroutine groupg(mode,ng,g,ag,plat,ngen,gen,agen,gens,ngout)
C- Finds a set of generators for the symmetry group
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode  :0 two groups compare to equal when both their point
Ci         :  and space parts compare equal
Ci         :1 two groups compare to equal when their point
Ci         :  group compares equal.  This eliminates
Ci         :  space groups that that have the same point group
Ci         :  but differing translational symmetry, which can
Ci         :  occur for artifically large supercells.
Ci   plat  :primitive translation vectors in real space
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   g:symmetry operation symbol
Ci   ng:number of symmetry operations as supplied by the generators
Co Outputs:
Co   gen,ngen:generators, and number needed to produce g
Co   ngout :number of group ops generated by (gen,ngen)
Co         :usually ngout=ng unless artificial translations
Co   gens  :ascii representation of generators
Cr Remarks:
Cr   The smallest set of generators is sought.
Cr   This subroutine performs the inverse function of sgroup.
Cu Updates
Cu   09 Jul 08 Extra check to find new generators beyond
Cu             the given ones.
Cu   12 May 07 Always returns gens, independent of verbosity
Cu   04 Jan 06 Returns ngout
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer mode,ngen,ng,ngout
      double precision plat(9)
      double precision gen(9,*),agen(3,*),g(9,*),ag(3,*)
      character*(*) gens
C Local parameters:
      integer imax,isop,ngloc,ngmax,iprint,ngen0,ngmx
      integer stdo,nglob
      parameter (ngmx=48*64)
      character*100 sout*120,sout2*120
      double precision gloc(3,3,ngmx),agloc(3,ngmx),qlat(9),xx

      stdo = nglob('stdo')

C --- Starting number of group ops ---
      call mkqlat(plat,qlat,xx)
      call pshpr(1)
      ngen0 = ngen
      call sgroup(0,gen,agen,ngen,gloc,agloc,ngout,ngmx,qlat)

   10 continue
C --- Do until enough generators added to make whole group ---
      if (ngout < ng) then
C   ... Run through all symops, choosing whichever adds the most ops
        imax = 0
        ngmax = 0
        do  isop = 1, ng
          call dcopy(9,g(1,isop),1,gen(1,ngen+1),1)
          call dcopy(3,ag(1,isop),1,agen(1,ngen+1),1)
C         call pshpr(61)
          call sgroup(mode,gen,agen,ngen+1,gloc,agloc,ngloc,ngmx,qlat)
C         call poppr
          if (ngloc > ngmax) then
            imax = isop
            ngmax = ngloc
            ngout = ngloc
          endif
        enddo
        ngen = ngen+1
        call dcopy(9,g(1,imax),1,gen(1,ngen),1)
        call dcopy(3,ag(1,imax),1,agen(1,ngen),1)
        goto 10
      endif

C     One last pass in case extra generators
      if (.true.) then
C   ... Run through all symops, choosing whichever adds the most ops
        imax = 0
        ngmax = ngout
        do  isop = 1, ng
          call dcopy(9,g(1,isop),1,gen(1,ngen+1),1)
          call dcopy(3,ag(1,isop),1,agen(1,ngen+1),1)
C         call pshpr(61)
          call sgroup(mode,gen,agen,ngen+1,gloc,agloc,ngloc,ngmx,qlat)
C         call poppr
          if (ngloc > ngmax) then
            imax = isop
            ngmax = ngloc
            ngout = ngloc
          endif
        enddo
        if (ngout > ngmax) then
          ngen = ngen+1
          call dcopy(9,g(1,imax),1,gen(1,ngen),1)
          call dcopy(3,ag(1,imax),1,agen(1,ngen),1)
        endif
      endif

      call poppr

C --- Create gens, optionally printout ---
      if (ngen0 == 0) then
        call info0(20,0,0,' GROUPG: the following '//
     .  'are sufficient to generate the space group:')
      else
        call info2(20,0,0,' GROUPG: %i generator(s) were added to '//
     .    'complete the group%?#n#:',ngen-ngen0,ngen-ngen0)
      endif
      sout = ' '
      sout2 = ' '
      call asymopn(ngen,gen,agen,plat,sout(9:),sout2(9:))
      if (ngen > ngen0 .and. iprint() >= 20) then
        call awrit0('%a',sout,len(sout),-stdo)
        call awrit0('%a',sout2,len(sout2),-stdo)
      endif
      gens = sout2
      if (ngout > ng) then
        call info2(20,0,0,
     .    '%9f(warning) %i group ops supplied but generators create'//
     .    ' %i ops',ng,ngout)
      endif

      end