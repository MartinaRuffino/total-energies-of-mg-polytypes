      subroutine mixpqc(nclass,nrclas,nl,nsp,igroup,pnu,qnu)
C- Averages together classes of equivalent groups
C ----------------------------------------------------------------
Ci Inputs
Ci   nclass,nrclas,nl,nsp,igroup
Ci   jgroup,pwk,qwk: work arrays
Co Outputs
Co   pnu,qnu are averaged by group
Cr Remarks
Cl Variables
Cl   cgrp is current group; grps is size of current group;
Cl   ngrp is number of groups; ic is class
C ----------------------------------------------------------------
      implicit none
      integer nclass,nrclas(0:*),nl,nsp,igroup(0:nclass)
      double precision pnu(nl,nsp,0:nclass-1),
     .                 qnu(3,nl,nsp,0:nclass-1)
C Local variables:
      integer ncmx
      parameter (ncmx=4096)
      integer ngrp,ic,i,grps,grpw,cgrp,iprint,jgroup(0:2*ncmx),isp,ispc,
     .  i0
      double precision pwk(10,2), qwk(3,10,2)

C --- Copy iabs(igroup) to jgroup ---
      if (nclass > ncmx) call rx('mixpqc: increase ncmx')
      call icopy(nclass,igroup,1,jgroup,1)
      do  ic = 0, nclass-1
        jgroup(ic) = iabs(jgroup(ic))
      enddo

C --- Sort groups and find first with group>0 ---
      call ivshel(1,nclass,jgroup,jgroup(nclass),.true.)
      do  i = 0, nclass-1
        i0 = i
        if (jgroup(jgroup(nclass+i)) /= 0) exit
      enddo
      if (i0 == nclass-1) return

C --- Accumulate P, Q by group ---
      if (iprint() >= 30) print '(/'' Mixpqc: ic  ngrp  group'')'
      ngrp = 0
      grps = 0
      grpw = 0
      cgrp = 0
      call dpzero(pwk,nl*nsp)
      call dpzero(qwk,nl*nsp*3)
C --- Loop through classes in order of permutation table ---
      do  i = i0, nclass-1
C   ... New group if group ID changes
        if (cgrp /= jgroup(jgroup(nclass+i))) then
          cgrp = jgroup(jgroup(nclass+i))
          if (grps /= 0)
     .    call mixpqx(i,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .      jgroup,pwk,qwk,pnu,qnu)
          grps = 0
          grpw = 0
          ngrp = ngrp+1
          call dpzero(pwk,20)
          call dpzero(qwk,60)
        endif
        grps = grps + 1
        ic = jgroup(nclass+i)
        grpw = grpw + nrclas(ic)
C   ... Sum P,Q of this class into its group
        do  isp = 1, nsp
          ispc = isp
          if (igroup(ic) < 0) ispc = nsp+1-isp
          call daxpy(nl,dble(nrclas(ic)),pnu(1,isp,ic),1,pwk(1,ispc),1)
          call daxpy(3*nl,dble(nrclas(ic)),qnu(1,1,isp,ic),1,
     .      qwk(1,1,ispc),1)
        enddo
      enddo
      call mixpqx(nclass,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .  jgroup,pwk,qwk,pnu,qnu)

      end
      subroutine mixpqx(i,cgrp,ngrp,grps,grpw,nclass,nl,nsp,igroup,
     .  jgroup,pwk,qwk,pnu,qnu)
C- Kernel called by mixpqc to redistribute group avg into each class
      implicit none
      integer i,ngrp,grpw,grps,nl,nsp,igroup(0:*),jgroup(0:*),
     .  cgrp,nclass
      integer, parameter :: n0=10
      double precision pnu(nl,nsp,0:*),qnu(3,nl,nsp,0:*),
     .  pwk(n0,2),qwk(3,n0,2),q
      integer i1,iprint,ic,isp,ispc,stdo,nglob
      double precision dsum


      stdo = nglob('stdo')
      call dscal(n0*nsp,1/dble(grpw),pwk,1)
      call dscal(3*n0*nsp,1/dble(grpw),qwk,1)
      q = dsum(nl,qwk,3)
      if (nsp == 2) q = q + dsum(nl,qwk(1,1,2),3)

C --- Redistribute P,Q of this group into each class ---
      do  i1 = i-grps, i-1
        ic = jgroup(nclass+i1)
        if (iprint() >= 40) write(stdo,1) ic+1,cgrp,ngrp
    1   format(i11,i5,i7)
        do  isp = 1, nsp
          ispc = isp
          if (igroup(ic) < 0) ispc = nsp+1-isp
          call dcopy(nl,pwk(1,ispc),1,pnu(1,isp,ic),1)
          call dcopy(3*nl,qwk(1,1,ispc),1,qnu(1,1,isp,ic),1)
        enddo
      enddo
      if (iprint() >= 30) write(stdo,2) grps,ngrp,q
    2 format(i6,' elts in group',i3,'  Q=',f9.6)

      end
