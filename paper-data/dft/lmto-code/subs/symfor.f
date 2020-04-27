      subroutine symfor(nbas,mode,g,ng,istab,f)
C- Symmetrize forces
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas:  number of atoms in basis
Ci   mode:  1 symmetrize by f(  ib ) = sum_g(ig) f(R(ib))
Ci          2 symmetrize by f(R(ib)) = sum_g(ig) f(  ib)
Ci   g,ng:  symmetry operations, and number
Ci   istab: site into which g,ag transforms site i
Co Outputs
Co   f:  forces are symmetrized
Cu Updates
Cu        2011 (DP) eliminated work array as passed argument
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ng
      integer istab(nbas,*)
      double precision g(3,3,1),f(3,nbas)
C ... Local parameters
      integer i,j,ib,ig,jb,mode
      real(8) :: fwk(3,nbas)
C ... External calls
      external dpcopy,dpzero

C     call prmx('f-in',f,3,3,nbas)
      call dpcopy(f,fwk,1,3*nbas,1d0/ng)
      call dpzero(f,3*nbas)

      if (mode == 1) then
        do  ig = 1, ng
          do  ib = 1, nbas
          jb = istab(ib,ig)
            do  i = 1, 3
              do  j = 1, 3
                f(i,ib) = f(i,ib) + g(i,j,ig)*fwk(j,jb)
              enddo
            enddo
          enddo
        enddo
      else
        do  ig = 1, ng
          do  ib = 1, nbas
          jb = istab(ib,ig)
            do  i = 1, 3
              do  j = 1, 3
                f(i,jb) = f(i,jb) + g(i,j,ig)*fwk(j,ib)
              enddo
            enddo
          enddo
        enddo
      endif

C     call prmx('f-sym',f,3,3,nbas)

      end
