      subroutine cpadlm(mode,nspc,lrel,lso,s_site,norb,pf,omega,cp,
     .  pmomg,lidim,iprm)
C- Coherent potential generator for ASA GF
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 insert coherent potential cohP in iprm order in cp
Ci         :2 return cohP-Omega for one site in orbital order in pmomg
Ci         :  (ignore iprm)
Ci         :3 do both 1 and 2
Ci         :10s digit
Ci         :0 do not include constraining fields
Ci         :1 include constraining fields
Ci   nspc  :number of coupled spins
Ci   lrel  :0 for nonrelativistic
Ci         :1 for scalar relativistic
Ci         :2 for fully  relativistic
Ci   lso   :.true. for spin-orbit coupling (should not occur with lrel=2)
Ci   norb  :dimension of the basis
Ci   pf    :array containing potential parameters P(z)
Ci   omega :coherent interactor (gfomg)
Ci   lidim :dimensions cp.  Also if 2s digit mode is zero, it is the
Ci         :dimension of the crystal basis (l+i waves)
Ci   iprm  :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci         :passed in for the given basis site
Co Outputs
Co   cp    :coherent potential for given site inserted in iprm order
Co   pmomg :cohP-Omega in orbital order for one site (ignore iprm)
Cu Updates
Cu   08 Jun 14 (Belashchenko) extended to spin-coupled case
Cu   06 Jan 13 New mode, and associated options
Cu   18 Dec 12 Completed migration to F90 structures
Cu   08 Oct 08 (Kirill) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nspc,lrel,norb,lidim,iprm(norb)
      logical lso
      complex(8) pf(norb,2,1+lrel/2,*),omega(norb,nspc,norb,2)
      complex(8) cp(lidim,nspc,lidim,2),pmomg(*)

c     complex(8), allocatable :: pf(:,:,:,:),omg(:,:,:,:)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site
C ... Dynamically allocated arrays
      real(8),allocatable    :: wrk(:)
      complex(8),allocatable :: wk(:,:,:,:),wk1(:,:,:,:)
C ... Local parameters
      integer job,icp,ierr,job1,job2,mod1,nl,ncomp
C     integer n1,n2,ip1,ip2
      real(8) eula(3),wt,pi,th
      real(8), allocatable :: bxc(:,:)

      pi = 4*datan(1d0)
      job = mod(mode,10)
      mod1 = mod(mode/10,10)
      call sanrg(.true.,lrel,0,2,'cpadlm:','lrel')
      call sanrg(.true.,job,0,3,'cpadlm:','job')
      call sanrg(.true.,mod1,0,1,'cpadlm:','mod1')
      if (job == 0) return
      job1 = mod(job,2)
      job2 = mod(job/2,2)
      ncomp = s_site%ncomp

c     if (job2 /= 0 .and. lidim /= norb)
c    .  call rx('cpadlm: job2 and lidim != norb')

      nl = sqrt(norb + 0.1)
      if (mod1 == 1) then
        allocate(bxc(3,ncomp))
        call dcopy(ncomp*3,s_site%bxc,1,bxc,1)
      endif

C     call yprm('pf in mkcpa',3,pf,1,norb,norb,ncomp*2)

C --- Make cohP-Omega in wk1
      allocate(wk(norb,2,norb,2),wk1(norb,nspc,norb,2),wrk(66*2*norb))
      wk1 = 0
      do  icp = 1, ncomp
        wt = s_site%cpawt(icp)
        if (wt < 1d-9) cycle
C   ... Record rotated Pa into wk
        if (lrel == 2 .or. lso) then
          eula(1) = s_site%thet(icp,2)
          eula(2) = s_site%thet(icp,1) ; eula(3) = - eula(1)
          if (mod1 == 1) eula(2) = eula(2) + datan(bxc(1,icp))
          if (lso) then
            call zcopy(norb*norb*4,s_site%pfr(1,icp),1,wk,1)
            call rot_LandS(10,eula,nl,norb,1,wk)
          else
            call rx('cpadlm: check rotation convention in call to rotpfr')
            call rotpfr(10,eula,nl,0,0,0,0,norb,pf(1,1,1,icp),norb,1,wk)
          endif
        else
          th = s_site%thet(icp,1)
          if (mod1 == 1) th = th + datan(bxc(1,icp))
          call adrotp(1,norb,pf(1,1,1,icp),pf(1,2,1,icp),th,wk)
        endif
C   ... Pa-Omega in global coordinate system
        if (nspc == 2) then
          wk = wk - omega(:,:,:,:)
        else
          wk(:,1,:,1) = wk(:,1,:,1) - omega(:,1,:,1)
          wk(:,2,:,2) = wk(:,2,:,2) - omega(:,1,:,2)
        endif
C   ... Add wa * (Pa-Omega)^-1 to wk1 array
        call zqinv('N',wk,2*norb,-66*2*norb,2*norb,wrk,2*norb,ierr)
        if (ierr /= 0) call rx('cpadlm:  wk is singular')
        if (nspc == 2) then
          wk1 = wk1 + wt * wk
        else
          wk1(:,1,:,1) = wk1(:,1,:,1) + wt * wk(:,1,:,1)
          wk1(:,1,:,2) = wk1(:,1,:,2) + wt * wk(:,2,:,2)
        endif
      enddo
C      call zprm('gii (up)',2,wk1,norb,norb,norb)
C      call zprm('gii (dn)',2,wk1(1,1,2),norb,norb,norb)

      if (allocated(bxc)) deallocate(bxc)

C ... Invert to make (Pcoh-Omega)
      if (nspc == 2) then
        call zqinv('N',wk1,2*norb,-66*2*norb,2*norb,wrk,2*norb,ierr)
      else
C ...   Each spin component separately if nspc = 1
        call zqinv('N',wk1,norb,-66*2*norb,norb,wrk,norb,ierr) ! spin 1
        if (ierr == 0) call zqinv('N',wk1(1,1,1,2),norb,-66*2*norb,
     .    norb,wrk,norb,ierr) ! spin 2
      endif
      if (ierr /= 0) call rx('cpadlm:  wk1 is singular')

C --- Record wk1 in pmomg (orbital order, ignore iprm)
      if (job2 == 1) call zcopy(norb*nspc*norb*2,wk1,1,pmomg,1)

C --- Record coherent potential in cp (iprm order)
      if (job1 == 1) then
        wk1 = wk1 + omega
        call pokeg0(0,norb,iprm,lidim,lidim,lidim,nspc,2,wk1,cp)
      endif

      deallocate(wk,wk1,wrk)

      end


      subroutine prtarr(wk1,norb)
      integer norb
      complex(8) wk1(norb,2,norb,2)
      integer n1
      print *,'block 11:'
      do  n1 = 1, norb
        write(6,905) wk1(n1,1,:,1)
      enddo
      print *,'block 22:'
      do  n1 = 1, norb
        write(6,905) wk1(n1,2,:,2)
      enddo
      print *,'block 12:'
      do  n1 = 1, norb
        write(6,905) wk1(n1,1,:,2)
      enddo
      print *,'block 21:'
      do  n1 = 1, norb
        write(6,905) wk1(n1,2,:,1)
      enddo
 905  format(16(2f7.3,1x))
      end
