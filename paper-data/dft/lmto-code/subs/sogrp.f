      subroutine sogrp(symops,ag,ng,imode)
C- Generate all point symmetry operations from the generation group
C ----------------------------------------------------------------
Ci Inputs
Ci   symops:rotation part of symmetry operation matrices
Ci   ag    :translation part of space group
Ci   ng    :number of group operations
Ci   imode: 1. Only take symmetry operations make spin quantization axis invariant.
Ci         -1. Also take symmetry operations flip the axis direction.
Cio Inputs/Outputs
Ci   symops:rotation part of symmetry operation matrices
Cr Remarks
Cr   This works for point groups only and is set up for integer generators.
Cr   liqin: About imode
cr      From what I test so far, it seems ASA is ok with -1 option, but full potential
Cr      probably need imode=1 for safe
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,imode
      double precision symops(9,ng),ag(3,ng)
C ... Local parameters
      double precision gi(3,3),symops2(9,ng),ag2(3,ng)
      double precision x0(3),y0(3),x2(3),y2(3),z2(3),z0(3)
      integer isymsoc(ng),ic1,ic2,ic
      integer ig,j,i,stdo,ng0
      procedure(integer) :: nglob,iprint
      character sg*70

      stdo = nglob('stdo')

      x0 = (/ 1d0,0d0,0d0 /)
      y0 = (/ 0d0,1d0,0d0 /)
      z0 = (/ 0d0,0d0,1d0 /)
      isymsoc = 0
      ic1 = 0
      ic2 = 0
C     For each group op, set isymsoc to
C       1 if  xrot x yrot points along z
C      -1 if  xrot x yrot points along -z
      do  i = 1, ng
        gi = reshape(symops(:,i), (/ 3, 3 /))
        x2 = matmul(gi,x0)  ! Rotated x axis
        y2 = matmul(gi,y0)  ! Rotated y axis
        call cross(x2,y2,z2)

        if (abs(z2(3)-1d0) < 1d-5) then
          isymsoc(i) = 1
          ic1 = ic1+1
        elseif (abs(z2(3)+1d0) < 1d-5) then
          isymsoc(i) = -1
          ic2 = ic2+1
        endif

        if (iprint() > 70) then
          write(stdo,"('!',i3 )")  i
          do  j = 1, 3
            write(stdo,"( 3f13.7, 5x,2f13.7 )") gi(j,1:3),ag(j,i),z2(j)
          enddo
        endif
      enddo

      ng0 = ng

C ... Return a reduced list of symops
      ic = 0
      if (imode == 1) then
        do  ig = 1, ng
          if (isymsoc(ig) == 1) then
            ic = ic+1
            symops2(:,ic) = symops(:,ig)
            ag2(:,ic) = ag(:,ig)
          endif
        enddo
        ng = ic
        symops(:,1:ng) = symops2(:,1:ng)
        ag(:,1:ng) = ag2(:,1:ng)
      elseif (imode == -1) then
        do  ig = 1, ng
          if ((isymsoc(ig)==1) .or. (isymsoc(ig) == -1)) then
            ic = ic+1
            symops2(:,ic) = symops(:,ig)
            ag2(:,ic) = ag(:,ig)
          endif
        enddo
        ng = ic
        symops(:,1:ng) = symops2(:,1:ng)
        ag(:,1:ng) = ag2(:,1:ng)
      endif

C --- Printout ---
      call info5(20,0,0,' SOGRP:  of %i symops, %i retained with z '//
     .  'preserved%?;(n==-1); + %i with z -> -z;;',
     .  ng0,ic1,imode,ic2,0)
      if (iprint() < 60 .or. ng < 2) return

      write(stdo,'('' ig  group op'')')
      do  ig = 1, ng
        call asymop(symops(1,ig),ag(1,ig),' ',sg)
        write(stdo,'(i4,2x,a)') ig,trim(sg)
      enddo

      end      subroutine sogrp
