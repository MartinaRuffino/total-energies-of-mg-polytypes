      subroutine ncutcorrect(ncut,n,gv,ng)
C- Adjust number of G vectors to include vecs of equivalent length near cutoff
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of orbitals for which to adjust cutoff
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   ng    :number of G-vectors
Cio Inputs/Outputs
Cio  ncut  :orbital-dependent maximum number of G vectors to use
Cr Remarks
Cr   Input ncut was generated for the Gamma point; but the G vectors
Cr   are q-dependent.  This routine increases the G vector cutoff
Cr   to include all G vectors of the same length.
Cu Updates
Cu   02 Apr 09 (TK) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,ng,ncut(n)
      double precision gv(ng,3)
C ... Local parameters
      integer:: i,ni,ig
      real(8):: gg, eps=1d-6

C     print *, ncut
      do  i = 1, n
        ni = ncut(i)
        if (ni < ng) then
          gg = sum(gv(ni,:)**2)
          do  ig = ni+1, ng
C           print *, i, ig, gg,sum(gv(ig,:)**2)
            if (abs(sum(gv(ig,:)**2)/gg-1d0) > eps) then
              ncut(i) = ig-1
              goto 99
            endif
          enddo
          ncut(i) = ng
   99     continue
        else
          ncut(i) = ng
        endif
      enddo
C     print *, ncut
      end
