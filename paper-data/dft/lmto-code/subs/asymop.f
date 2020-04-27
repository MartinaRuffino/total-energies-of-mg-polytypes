      subroutine asymop(grp,ag,asep,sg)
C- Generate the symbolic representation of a group operation
C ----------------------------------------------------------------------
Ci Inputs:
Ci  grp,ag :  space group rotation + translation matrix
Ci  asep:     If first character is 'f', rot. and trans. vectors
Ci            are not converted into small algebraic expressions.
Ci            If first char. is 'f', The second-* character are
Ci            used for the separators.
Co Outputs:
Co  sg  :  symbolic representation of group op
Cb Bugs
Cb  No check is made on the length of sg
C ----------------------------------------------------------------------
      implicit none
      double precision grp(3,3),ag(3)
      character*(*) sg,asep
C Local variables
      double precision vecg(3),dasum,tiny
      integer nrot,ip,isw,awrite,i1,i2,fmtv
      logical li,parsvc
      parameter(tiny=1d-4)

C --- Get consitutents of grp ---
      call csymop(1,grp,li,nrot,vecg)

C --- Rotational part ---
      i1 = 1
      fmtv = 0
      if (asep(1:1) == 'f') then
        fmtv = 4
        i1 = 2
      endif
      sg = ' '
      if (nrot == 1) then
        sg = 'i*i'
        ip = 3
        if (li) sg = 'i'
        if (li) ip = 1
      else
        if (li .and. nrot == 2) then
          sg = 'm'
          ip = 1
        else
          ip = awrite('%?#n#i*##r%i',sg,len(sg),0,isw(li),nrot,
     .      0,0,0,0,0,0)
        endif
        call rxx(.not. parsvc(2+fmtv,sg,ip,vecg),'bug in asymop')
      endif

C --- Translational part ---
      if (dasum(3,ag,1) > tiny) then
        if (asep(i1:i1) /= ' ') then
          call nword(asep,1,i1,i2)
          sg(ip+1:) = asep(i1:i2)
          ip = ip+i2-i1+1
        endif
        call rxx(.not. parsvc(1+fmtv,sg,ip,ag),'bug in asymop')
      endif

      end

      subroutine asymopn(ngen,gen,agen,plat,sgenc,sgenp)
C- Generate the symbolic representation of a family of group operations
C ----------------------------------------------------------------------
Ci Inputs
Ci   ngen  :number of operations
Ci   gen   :point group for each group op
Ci   agen  :translation for each group op, Cartesian coordinates
Ci          (without alat)
Ci   plat  :primitive lattice vectors
Co Outputs
Co   sgenc :ascii representation, agen in Cartesian coordinates
Ci   sgenp :ascii representation, agen in units of plat
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   26 Feb 14
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngen
      double precision gen(9,ngen),agen(3,ngen),plat(3,3)
      character sgenc*(*), sgenp*(*)
C ... Local parameters
      integer isop,ngmx
      integer i1,i2,j1,j2
      parameter (ngmx=48*64)
      character*100 sg,sg1
      double precision qlat(3,3),vec(3)

      call mkqlat(plat,qlat,vec)
      sgenc = ' '; sgenp = ' '
      do  isop = 1, ngen
        call asymop(gen(1,isop),agen(1,isop),':',sg)
        call awrit0('%a '//sg,sgenc,len(sgenc),0)
        call dcopy(3,agen(1,isop),1,vec,1)
        call dgemm('N','N',1,3,3,1d0,agen(1,isop),1,qlat,3,0d0,vec,1)
        call asymop(gen(1,isop),vec,'::',sg1)
        call word(sg1,1,i1,i2)
        call shorbz(vec,vec,plat,qlat)
        call asymop(gen(1,isop),vec,'::',sg)
        call word(sg,1,j1,j2)
        if (i2-i1 < j2-j1) sg = sg1
        call awrit0('%a '//sg,sgenp,len(sgenp),0)
      enddo
      end

      subroutine csymop(iopt,grp,li,nrot,vecg)
C- Decomposes a group operation into its consitutents, or vice-versa
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt  := -1 to convert (nrot,vecg,li) to grp
Ci          =  1 to convert grp to to (nrot,vecg,li)
Co Inputs/Outputs:
Cio grp   :group operation matrix
Cio li    :if T: inversion or rotoinversion
Cio nrot  :rotation angle = 2*pi/nrot
Cio vecg  :rotation axis
Cr Remarks
Cr   for nrot > 2 the matrix is non-symmetric and the rotation
Cr   axis can be calculated from the antisymmetric part.
Cr   For nrot = 2 this not possible.  However, the squared vector
Cr   components are given by:  mat(i,i) = 2 v_i * v_i - 1.
Cr   This is used for the largest component. The others are taken
Cr   from: mat(i,j) = 2 v_i * v_j for i ne j.  This way we also
Cr   get the right phases between the components.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer nrot,iopt
      double precision vecg(3),grp(3,3)
      logical li
C Local parameters:
      integer i,idamax,j,in
      double precision costbn,detop,ddet33,dnrm2,sinpb3,tiny,twopi,vfac,
     .  wk(9),sintbn,omcos,ddot
      parameter(tiny=1d-5)
C External calls:
      external daxpy,dcopy,ddet33,dpzero,dnrm2,dscal,idamax,ddot

      twopi = 8*datan(1d0)
C --- Make grp from (nrot,vecg,li) ---
      if (iopt == -1) then
        call dpzero(grp,9)
        in = iabs(nrot)
        if (in <= 0.or.in == 5.or.in > 6)
     .    call fexit(-1,111,'%N Exit -1 CSYMOP: '//
     .    'abs(nrot) must 1,2,3,4 or 6, but is %i',in)
        if (in == 1) then
          call dcopy(3,1d0,0,grp,4)
        else
          sintbn = dsin(twopi/nrot)
          costbn = dcos(twopi/nrot)
          omcos  = 1d0-costbn
          call rxx(dnrm2(3,vecg,1) < tiny,
     .             'CSYMOP: zero rotation vector')
          call dscal(3,1/sqrt(ddot(3,vecg,1,vecg,1)),vecg,1)
          grp(1,1) = omcos*vecg(1)*vecg(1) + costbn
          grp(1,2) = omcos*vecg(1)*vecg(2) - sintbn*vecg(3)
          grp(1,3) = omcos*vecg(1)*vecg(3) + sintbn*vecg(2)
          grp(2,1) = omcos*vecg(2)*vecg(1) + sintbn*vecg(3)
          grp(2,2) = omcos*vecg(2)*vecg(2) + costbn
          grp(2,3) = omcos*vecg(2)*vecg(3) - sintbn*vecg(1)
          grp(3,1) = omcos*vecg(3)*vecg(1) - sintbn*vecg(2)
          grp(3,2) = omcos*vecg(3)*vecg(2) + sintbn*vecg(1)
          grp(3,3) = omcos*vecg(3)*vecg(3) + costbn
        endif
        if (li) call dscal(9,-1d0,grp(1,1),1)

C --- Make (nrot,vecg,li) from grp ---
      else if (iopt == 1) then

C ... Require |determinant=1|
        call dinv33(grp,0,wk,detop)
        if (dabs(dabs(detop)-1.0d0) > tiny)
     .    call fexit(-1,111,'%N Exit -1 ASYMOP: '//
     .    'determinant of group op must be +/- 1, but is %d',detop)
        detop = dsign(1.d0,detop)
C   ... li is T if to multiply by inversion
        li = detop < 0d0
C   ... Multiply operation grp with detop to guarantee pure rotation
        call dscal(9,detop,grp(1,1),1)
C   --- Calculate rotation angle from the normalization of v ---
C       sum_i grp(i,i) = sum_i (1-cos) v_i*v_i + 3*cos = 1 + 2 * cos
C       costbn = -0.5d0
C       call daxpy(3,0.5d0,grp(1,1),4,costbn,0)
        costbn = 0.5d0*(-1 + grp(1,1) + grp(2,2) + grp(3,3))

        if (dabs(costbn-1d0) < tiny) then
          nrot = 1
          call dpzero(vecg,3)
        else
C     ... See Remarks
          nrot = idnint(twopi/dacos(dmax1(-1d0,costbn)))
          if (nrot == 2) then
            do  i = 1, 3
              vecg(i) = 0.5d0*(grp(i,i)+1.0d0)
            enddo
            j = idamax(3,vecg,1)
            if (vecg(j) < 0d0)
     .        call fexit2(-1,111,' Exit -1 ASYMOP:  bad component %i'//
     .        ' of operation.  Diagonal element = %d',j,grp(j,j))
            vecg(j) = dsqrt(vecg(j))
            vfac = 0.5d0/vecg(j)
            do  i = 1, 3
              if (i /= j) vecg(i) = vfac*grp(i,j)
            enddo
          else
            vecg(1) = grp(3,2)-grp(2,3)
            vecg(2) = grp(1,3)-grp(3,1)
            vecg(3) = grp(2,1)-grp(1,2)
          endif

C     --- Renormalize at least one component to 1 ---
C         to allow for abbreviations as 'D', 'X', 'Y' or 'Z'
          sinpb3 = dsqrt(.75d0)
          if (dabs((sinpb3-dabs(vecg(1)))*(sinpb3-dabs(vecg(2)))*
     .            (sinpb3-dabs(vecg(3)))) > tiny) then
            do  j = 3, 1, -1
             vfac = dabs(vecg(j))
             if(vfac > tiny) call dscal(3,1.d0/vfac,vecg,1)
            enddo
          endif
        endif
        call dscal(9,detop,grp(1,1),1)

      endif
      end
