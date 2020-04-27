!      module specfun
!      contains
      subroutine prepsf0(norb,nspc,s_site,lidim,iprm,pfdim,dpf,ddpf,am,
     .  amt,al)
C- Make L,M,M~ matrices for non-CPA sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   norb  :number of orbitals for this site
Ci   nspc  :number of coupled spins (if =2 then input matrices in s_site)
Ci   lidim :number of local+intermediate orbitals in the Hamiltonian
Ci   iprm  :permutations ordering orbitals for given site (makidx.f)
Ci   pfdim :leading dimension of the potential parameter arrays
Ci   dpf   :1st energy derivative of potential function P-dot (mkptfp.f)
Ci   ddpf  :2nd energy derivative  -1/2 P-dotdot/P-dot (mkptfp.f)
Co Outputs
Co   al    :Matrices needed for spectral function
Co   am,amt:al = Lambda, am = M, amt = M-tilde: see Turek's book (4.26)
Cu Updates
Cu   23 Jan 14 (Belashchenko) Cleaned up
Cu   00 Feb 13 (Pujari) First created
C ----------------------------------------------------------------------
C!DEC$ ATTRIBUTES NO_ARG_CHECK :: s_site
      use structures
      implicit none
C Passed parameters:
      integer norb,nspc,lidim,pfdim,iprm(norb)
      complex(8),dimension(pfdim,2)          :: dpf,ddpf
      complex(8),dimension(norb,nspc,norb,2) :: al,am,amt
C ..  For structures
!      include 'structures.h'
      type(str_site)::  s_site
C Local parameters:
      integer n,indx

      if (nspc == 1) then
        al = 0 ; am = 0 ; amt = 0
        do  n = 1, norb
          indx = iprm(n) ; if (indx > lidim) cycle
          al(n,1,n,:) = ddpf(indx,:)
          am(n,1,n,:) = sqrt(dpf(indx,:))
          amt(n,1,n,:) = am(n,1,n,:)
        enddo
      elseif (nspc == 2) then
        call zcopy(norb*norb*4,s_site%ddpfr,1,al,1)
        call zcopy(norb*norb*4,s_site%dpfr,1,am,1)
        amt(:,1,:,1) = transpose(am(:,1,:,1))
        amt(:,1,:,2) = transpose(am(:,2,:,1))
        amt(:,2,:,1) = transpose(am(:,1,:,2))
        amt(:,2,:,2) = transpose(am(:,2,:,2))
      else
        call rx1('Unexpected value nspc=',nspc)
      endif

      end subroutine prepsf0

      subroutine prepsf(mode,ncomp,nspc,s_site,th,wt,norb,pf,dpf,ddpf,
     .  omg,pmomg,am,amt,al)
C- Make L,M,M~ matrices for one CPA site, see Turek (4.26)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do not include constraining fields
Ci         :1 include constraining fields
Ci   ncomp :number of CPA components for this site
Ci   nspc  :number of coupled spins
Ci   th    :array of angles in radians
Ci   wt    :CPA weights, assumed normalized
Ci   norb  :number of orbitals for this site
Ci   pf    :potential functions
Ci   dpf   :1st energy derivative of potential function P-dot (mkptfp.f)
Ci   ddpf  :2nd energy derivative  -1/2 P-dotdot/P-dot (mkptfp.f)
Ci   omg   :Omega(z) for this site
Ci   pmomg :cohP-Omega for this site (cohP = coherent potential)
Ci         :(note it is overwritten in this routine)
Co Outputs
Co   al    :Matrices needed for spectral function
Co   am,amt:al = Lambda, am = M, amt = M-tilde: see Turek's book (4.26)
Cu Updates
Cu   23 Jan 14 (Belashchenko) Cleaned up
Cu   00 Feb 13 (Pujari) First created
C ----------------------------------------------------------------------
C!DEC$ ATTRIBUTES NO_ARG_CHECK :: s_site
      use gcpa
      use structures
      implicit none
C Passed parameters:
      integer,intent(in) :: mode,norb,nspc,ncomp
      real(8),intent(in) :: wt(ncomp),th(ncomp,nspc)
      complex(8),intent(in),dimension(norb,2,nspc,ncomp) :: pf,dpf,ddpf
      complex(8),intent(in) :: omg(norb,nspc,norb,2)
      complex(8),intent(inout),dimension(norb,nspc,norb,2) :: pmomg,
     .  al,am,amt
C ..  For structures
!      include 'structures.h'
      type(str_site)::  s_site
C Local parameters:
      integer icp,n1,n2,is1,is2,ierr,lrel,nl
      real(8) pi,wei,thi(nspc),eula(3),bxc(3,ncomp)
      complex(8) mu1(2),mu2(2)
      complex(8), dimension(norb,2,norb,2) :: wk,wk2,wk3
      real(8), allocatable :: wrk(:)
      complex(8), allocatable :: wkpmo(:,:)
      logical lso

      call sanrg(.true.,mode,0,1,'prepsf:','mode')
      pi = 4*datan(1d0)
      lrel = 1 ; lso = .false.
      if (nspc == 2) lso = .true.

      al = 0 ; am = 0 ; amt = 0

      if (mode == 1) call dcopy(3*ncomp,s_site%bxc,1,bxc,1)

      do  icp = 1, ncomp
        thi = th(icp,:) ; wei = wt(icp)
        if (mode == 1) thi(1) = thi(1) + datan(bxc(1,icp))
        if (lso) then
          call gc00(0,norb,nspc,lrel,lso,thi,pfr=s_site%pfr(:,icp),
     .      omg=omg,g00=wk)
        else
          call gc00(0,norb,nspc,lrel,lso,thi,pfi=pf(:,:,:,icp),
     .      omg=omg,g00=wk)
        endif

C   --- mu*(P - Omega)^-1*mu
        if (nspc == 1) then
          wk2 = 0 ; wk3 = 0
          do  n2 = 1, norb
            mu2(1) = sqrt(dpf(n2,1,1,icp))
            mu2(2) = sqrt(dpf(n2,2,1,icp))
            do  n1 = 1,norb
              mu1(1) = sqrt(dpf(n1,1,1,icp))
              mu1(2) = sqrt(dpf(n1,2,1,icp))
              do  is1 = 1, 2
              do  is2 = 1, 2
                wk2(n1,is1,n2,is2) =          wk(n1,is1,n2,is2)*mu2(is2)
                wk3(n1,is1,n2,is2) = mu1(is1)*wk(n1,is1,n2,is2)
                wk (n1,is1,n2,is2) = mu1(is1)*wk(n1,is1,n2,is2)*mu2(is2)
              enddo
              enddo
            enddo
            do  is1 = 1, 2
              wk(n2,is1,n2,is1) = wk(n2,is1,n2,is1) + ddpf(n2,is1,1,icp)
            enddo
          enddo

C     ... Rotate full \Lambda and M to the global frame
          call rotm(norb,thi(1),wk)
          call rotm(norb,thi(1),wk2)
          call rotm(norb,thi(1),wk3)

C     ... Sum and extract the diagonal blocks after multiplying them by wei
          do  is1 = 1, 2
            al(:,1,:,is1) = al(:,1,:,is1) + wei * wk(:,is1,:,is1)
            am(:,1,:,is1) = am(:,1,:,is1) + wei * wk2(:,is1,:,is1)
            amt(:,1,:,is1) = amt(:,1,:,is1) + wei * wk3(:,is1,:,is1)
          enddo
        elseif (lso) then
          n2 = 2*norb
          wk2 = reshape(matmul(reshape(wk,(/n2,n2/)),
     .                         reshape(s_site%dpfr(:,icp),(/n2,n2/))),
     .          (/norb,2,norb,2/))
          wk3 = reshape(matmul(
     .               transpose(reshape(s_site%dpfr(:,icp),(/n2,n2/))),
     .                         reshape(wk,(/n2,n2/))),
     .          (/norb,2,norb,2/))
          wk  = reshape(matmul(
     .               transpose(reshape(s_site%dpfr(:,icp),(/n2,n2/))),
     .                         reshape(wk2,(/n2,n2/))),
     .          (/norb,2,norb,2/))
     .        + reshape(s_site%ddpfr(:,icp),(/norb,2,norb,2/))
          eula(1) = thi(2)
          eula(2) = thi(1) ; eula(3) = - eula(1)
          nl = sqrt(norb + 0.1)
          call rot_LandS(10,eula,nl,norb,1,wk)
          call rot_LandS(10,eula,nl,norb,1,wk2)
          call rot_LandS(10,eula,nl,norb,1,wk3)
          al = al + wei*wk ; am = am + wei*wk2 ; amt = amt + wei*wk3
        else
          call rx('prepsf: nspc=2 but lso=F')
        endif
      enddo  ! icp loop

C ... Multiply by (cohP-Omg)
      if (nspc == 1) then
        do  is1 = 1, 2
          am(:,1,:,is1) = matmul(pmomg(:,1,:,is1),am(:,1,:,is1))
          amt(:,1,:,is1) = matmul(amt(:,1,:,is1),pmomg(:,1,:,is1))
        enddo
      else
        n2 = 2*norb
        allocate(wkpmo(n2,n2))
        wkpmo = reshape(pmomg,(/n2,n2/))
        am = reshape(matmul(wkpmo,
     .                      reshape(am,(/n2,n2/))),(/norb,2,norb,2/))
        amt= reshape(matmul(reshape(amt,(/n2,n2/)),wkpmo),
     .                                             (/norb,2,norb,2/))
      endif

C ... Subtract ~M * (CohP - Omega)^-1* M from Lambda, see Turek (4.25)
C ... Note p_omg is overwritten
      allocate(wrk(66*2*norb))
      if (nspc == 1) then
        do  is1 = 1, 2
          call zqinv('N',pmomg(:,1,:,is1),norb,-66*norb,norb,wrk,norb,
     .      ierr)
          if (ierr /= 0) call rx('prepsf: singular matrix')
          al(:,1,:,is1) = al(:,1,:,is1) - matmul(amt(:,1,:,is1),
     .      matmul(pmomg(:,1,:,is1),am(:,1,:,is1)))
        enddo
      else
        call zqinv('N',wkpmo,n2,-66*n2,n2,wrk,n2,ierr)
        if (ierr /= 0) call rx('prepsf: singular matrix')
        al = al - reshape(matmul(reshape(amt,(/n2,n2/)),
     .      matmul(wkpmo,reshape(am,(/n2,n2/)))),(/norb,2,norb,2/))
        deallocate(wkpmo)
      endif
      deallocate(wrk)
      end subroutine prepsf

      subroutine cpasf(mode,nspc,s_site,zp,wz,nl,nbas,lidim,iprm,nkp,qp,
     .  theta,sf,sfint,sfr,sfrl)
C- Calculate the CPA spectral function
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb alr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:amrt amr
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 no site-resolved SF
Ci         :1 make site-resolved SF and return in sfr
Ci         :2 also make orbital resolved SF and return in sfrl
Ci         :10s digit is ignored
Ci         :100s digit
Ci         :1 assume integration over contour, use wz
Ci   nspc  :number of coupled spins
Ci   zp    :complex energy
Ci   wz    :weight of this z-point (if integration requested)
Ci   nbas  :number of basis sites
Ci   lidim :number of local+intermediate orbitals in the Hamiltonian
Ci   iprm  :permutations ordering orbitals for given site (makidx.f)
Ci   nkp   :number of q-points
Ci   qp    :array of q-points
Ci   theta :rotate the GF in spin space by this angle around the y axis
Co Outputs
Co   sf    :spectral function for this energy point
Co   sfr   :site-resolved spectral function for this energy point
Co   sfint :contribution to integrated SF (mode 100s digit)
Cu Updates
Cu   23 Jan 14 (Belashchenko) Cleaned up
Cu   00 Feb 13 (Pujari) First created
C ----------------------------------------------------------------------
C!DEC$ ATTRIBUTES NO_ARG_CHECK :: s_site
      use structures
      implicit none
C ... Passed parameters
      integer,intent(in) :: mode,nspc,nl,nbas,lidim,nkp,iprm(*)
      real(8),intent(in) :: qp(3,nkp),theta
      real(8) sf(2,nkp),sfr(*),sfrl(*),sfint(2,nkp)
      complex(8),intent(in) :: wz,zp
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(nbas)
C Local parameters:
      logical lsite,lorb
      integer mod0,mod2,fopnx,ifi,iogfrs,lio,clp(18)
      integer nglob,mxorb
      integer isp,jsp,ksp,iq,ib,norb,n,n2,idx
      real(8) plat(3,3),xx
      complex(8) trg(nspc)
      complex(8), dimension(:,:,:,:), allocatable :: gk,gkii,Gbb
      integer, parameter:: ISH=1,MCD=32,MZP=64,MNC=1,MNB=2,MNT=4,MLB=8
      real(8), parameter:: pi=3.14159265359
      character*20 fnam
      real(8),allocatable :: wksf(:,:,:),wksfl(:,:,:,:),wksfi(:,:)
C ... Parallel mode
      integer mpipid,procid,master,numprocs
      integer, allocatable :: kpproc(:)

      mxorb = nglob('mxorb')
      mod0 = mod(mode,10)
      mod2 = mod(mode/100,10)
      call sanrg(.true.,mod0,0,2,'cpasf:','mod0')
      call sanrg(.true.,mod2,0,1,'cpasf:','mod2')
      procid = mpipid(1) ; numprocs = mpipid(0) ; master = 0

      lsite = .false. ; lorb = .false.
      if (mod0 >= 1) lsite = .true.
      if (mod0 == 2) lorb = .true.
      plat = 0 ; clp = 0 ; clp(3) = lidim*nspc ; clp(4) = lidim*nspc

      allocate(gk(lidim,nspc,lidim,nspc))
      if (lsite) allocate(wksf(2,nbas,nkp))
      if (lorb) then
        allocate(wksfl(2,nl*nl,nbas,nkp))
        wksfl = 0
      endif
      if (mod2 == 1) then
        allocate(wksfi(2,nkp))
        wksfi = 0
      endif

      sf = 0
C --- Select k-point range for this process
      allocate (kpproc(0:numprocs))
      if (numprocs > 1) then
        call pshpr(1)
        call dstrbp(nkp,numprocs,1,kpproc(0))
        call poppr
      else
        kpproc(0) = 1 ; kpproc(1) = nkp+1
      endif

      do  isp = 1, 3-nspc
C  ...  Read the header of gfqp file
        fnam = 'gfqp1' ; if (isp == 2) fnam = 'gfqp2'
        if (nspc == 2) fnam = 'gfqp'
        ifi = fopnx(fnam,100,16+8+4+0,-1) ; rewind ifi
        call pshpr(0)
        lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MCD*1+MZP*0)
        if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zp,qp,plat,xx,clp,xx,xx,
     .    0,0,0,xx) /= 0) call rx('cpasf failed to read file header')
        call poppr
C  ...  For each q-point, read g(k) from disk, make G_BB(k), and make the spectral function
        do  iq = kpproc(procid), kpproc(procid+1)-1
          call pshpr(0)
          lio = 10000*ISH+0*MZP+10*MCD+6
          if (iogfrs(lio,1,0,' ',ifi,1,1,0,0,zp,qp(1,iq),xx,xx,clp,xx,
     .      xx,0,0,1,gk) /= 0)
     .      call rxi('cpasf failed to read gf for qp no.',iq)
          call poppr
          do  ib = 1, nbas
            norb = s_site(ib)%norb ; n2 = 2*norb
            allocate(gkii(norb,nspc,norb,nspc),Gbb(norb,nspc,norb,nspc))
C      ...  Extract GF for current site
            idx = mxorb*(ib-1)
C           call pokeg0(1,norb,iprm(1+idx),lidim,lidim,lidim,nspc,2,gkii,gk) ! Here gk is for one spin, pokeg0 can't process
            gkii = gk(iprm(idx+1:idx+norb),:,iprm(idx+1:idx+norb),:)
C      ...  Make G_BB = L - ~M*(CohP-Omg)^-1*M + ~M*g(k)*M
            if (nspc == 1) then
              Gbb(:,1,:,1) = matmul(s_site(ib)%amrt(:,1,:,isp),
     .          matmul(gkii(:,1,:,1),s_site(ib)%amr(:,1,:,isp)))
     .          + s_site(ib)%alr(:,1,:,isp)
            else
              Gbb = reshape(
     .          matmul(reshape(s_site(ib)%amrt,(/n2,n2/)),
     .                 matmul(reshape(gkii,(/n2,n2/)),
     .                        reshape(s_site(ib)%amr,(/n2,n2/))
     .                      )
     .                ),   (/norb,2,norb,2/)) + s_site(ib)%alr
            endif
            if (theta /= 0) call rotm(norb,theta,Gbb)
C       ... Store the spectral functions according to mode
            trg = 0
            do n = 1, norb
              do  jsp = 1, nspc
                trg(jsp) = trg(jsp) + Gbb(n,jsp,n,jsp)
                if (lorb) then
                  wksfl(max(isp,jsp),n,ib,iq) =
     .              - dimag(Gbb(n,jsp,n,jsp))/pi
                endif
              enddo
            enddo
            do  jsp = 1, nspc
              ksp = max(isp,jsp)
              if (lsite) wksf(ksp,ib,iq) = - dimag(trg(jsp))/pi
              sf(ksp,iq) = sf(ksp,iq) - dimag(trg(jsp))/pi
              if (mod2 == 1) then
                wksfi(ksp,iq) = wksfi(ksp,iq) - dimag(wz*zp*trg(jsp))/pi
              endif
            enddo
            deallocate(gkii,Gbb)
          enddo
        enddo  ! k loop
        call fclr(fnam,ifi)
      enddo  ! spin loop
      deallocate(gk)
      call mpibc4(sf,sf,kpproc,2,4)
      if (lsite) then
        call mpibc4(wksf,wksf,kpproc,2*nbas,4)
        call dcopy(2*nbas*nkp,wksf,1,sfr,1)
        deallocate(wksf)
      endif
      if (lorb) then
        call mpibc4(wksfl,wksfl,kpproc,2*nl*nl*nbas,4)
        call dcopy(2*nl*nl*nbas*nkp,wksfl,1,sfrl,1)
        deallocate(wksfl)
      endif
      if (mod2 == 1) then
        call mpibc4(wksfi,wksfi,kpproc,2,4)
        call dpsadd(sfint,wksfi,2*nkp,1,1,1d0)
        deallocate(wksfi)
      endif
      deallocate(kpproc)
 901  format(f10.5,i5,2f14.6)
      end subroutine cpasf

      subroutine iospf(mode,nl,ef,nzp,zp,nkp,qp,sf,nbas,sfr,sfrl)
C- Write the spectral function to file spf.ext and other files
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 1s digit
Ci         : 0 no site-resolved SF
Ci         : 1 make site-resolved SF and return in sfr
Ci         : 2 also make orbital resolved SF and return in sfrl
Ci         : 10s digit
Ci         : 0 assume k-points are along lines, arrange on one line
Ci         : 1 assume k-points are in BZ, list them in 3d format
Ci   nl    :(global maximum l) + 1
Ci   ef    : Shift spectral function by -ef (choose ef to be the zero)
Ci   nzp   : number of complex energy points
Ci   zp    : points for complex energy contour
Ci   nkp   : number of k-points
Ci   qp    : array of k-points
Ci   sf    : spectral function to write to file
Ci   sfr   : site-resolved spectral function to write to file
Cu Updates
Cu   23 Jan 14 (Belashchenko) Cleaned up
Cu   00 Feb 13 (Pujari) First created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer mode,nl,nzp,nkp,nbas
      real(8) qp(3,nkp),zp(2,nzp),sf(2,nkp,nzp),sfr(*),sfrl(*),ef
C Local parameters
      integer mod0,mod1,izp,iq,ib,ibmax,ncut,icut,ifi(0:nbas),stdo
      real(8) dq(3),qdist,distprev,sfl(2)
      real(8), allocatable :: qcut(:),qxlen(:)
      real(8) sfll(2,nl*nl),dos(2)
      real(8),allocatable :: wksf(:,:,:,:),wksfl(:,:,:,:,:)
      logical lsite,lorb,ldebug
      character*8 fnam
      procedure(integer) :: fopnx,nglob

      ldebug = .true.

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      call sanrg(.true.,mod0,0,2,'iospf:','mod0')
      call sanrg(.true.,mod1,0,1,'iospf:','mod1')
      stdo = nglob('stdo')

      lsite = .false. ; lorb = .false.
      if (mod0 >= 1) lsite = .true.
      if (mod0 == 2) lorb = .true.

C --- Find the path cuts in the set of k-points (heuristic)
      if (mod1 == 0) then
        allocate(qcut(nkp),qxlen(nkp))
        ncut = 0 ; icut = 1 ; qcut = 0 ; qxlen = 0 ; distprev = 0
        if (nkp > 2) then
          do iq = 2, nkp
            dq = qp(:,iq)-qp(:,iq-1) ; qdist = sqrt(dot_product(dq,dq))
            if (qdist < 1d-6) then
              ncut = ncut + 1 ; qcut(ncut) = qxlen(iq-1)
              icut = 1 ; qxlen(iq) = qxlen(iq-1) + qdist
            elseif (qdist > distprev*1.01 .and. icut == 0) then
              qxlen(iq) = qxlen(iq-1) ! + 0.02
              ncut = ncut + 1 ; qcut(ncut) = qxlen(iq-1)
              ncut = ncut + 1 ; qcut(ncut) = qxlen(iq)
              icut = 1
            else
              qxlen(iq) = qxlen(iq-1) + qdist ; icut = 0
            endif
            distprev = qdist
          enddo
          if (icut == 0) then
            ncut = ncut + 1 ; qcut(ncut) = qxlen(nkp)
          endif
        endif
      endif

C --- Write spectral function to file. Header includes locations of path cuts

      ibmax = 0
      if (lsite) then
        ibmax = nbas
        allocate(wksf(2,nbas,nkp,nzp))
        call dcopy(2*nbas*nkp*nzp,sfr,1,wksf,1)
        if (lorb) then
          allocate(wksfl(2,nl*nl,nbas,nkp,nzp))
          call dcopy(2*nl*nl*nbas*nkp*nzp,sfrl,1,wksfl,1)
        endif
      endif

C     ib = 0 writes the entire spectral function
C     ib > 0 writes site- or site+orbital resolved spectral function
      do  ib = 0, ibmax
        fnam = 'spf'
        if (ib > 0)  call awrit1('%a%i',fnam,8,0,ib)
        ifi(ib) = fopnx(fnam,100,16+8+2+0,-1)
        if (mod1 == 0) write(ifi(ib),901) qcut(1:ncut)
        do izp = 1, nzp
          if (mod1 == 1 .and. ldebug) then
            dos(1) = sum(sf(1,:,izp))/nkp
            dos(2) = sum(sf(2,:,izp))/nkp
            write(stdo,951) zp(1,izp)-ef,dos
 951        format(' DOS from SF:',3F12.6)
          endif
          do iq = 1, nkp
            if (ib == 0) then
              sfl = sf(1:2,iq,izp)
            else
              sfl = wksf(1:2,ib,iq,izp)
              if (lorb) sfll = wksfl(1:2,:,ib,iq,izp)
            endif
            if (ib == 0) then
              if (mod1 == 0) then
                write(ifi(ib),902) zp(1,izp)-ef,qxlen(iq),sfl
              else
                write(ifi(ib),903) zp(1,izp)-ef,qp(:,iq),sfl
              endif
            else
              if (lorb) then
                if(mod1 == 0) then
                  write(ifi(ib),902) zp(1,izp)-ef,qxlen(iq),sfl,sfll
                else
                  write(ifi(ib),903) zp(1,izp)-ef,qp(:,iq),sfl,sfll
                endif
              else
                if(mod1 == 0) then
                  write(ifi(ib),902) zp(1,izp)-ef,qxlen(iq),sfl
                else
                  write(ifi(ib),903) zp(1,izp)-ef,qp(:,iq),sfl
                endif
              endif
            endif
          enddo
        enddo
C       call fclose(ifi(ib))
        call fclr(' ',ifi(ib))
      enddo

      if (allocated(qcut)) deallocate(qcut,qxlen)

 901  format('#',100f9.5)
 902  format(2f10.5,50f14.6)
 903  format(4f10.5,50f14.6)
      if (lsite) deallocate(wksf)
      if (lorb)  deallocate(wksfl)
c     call fclose(ifi)
c     do  ib = 1, nbas
c       call fclose(ifi(ib))
c     enddo
      end subroutine iospf

      subroutine iospf2(nkp,qp,sf)
C- Write the integrated spectral function to file sfint.ext
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkp   : number of k-points
Ci   qp    : array of k-points
Ci   sf    : k-resolved quantity to write to file
Cu Updates
Cu   29 Jul 14 (Belashchenko) First created
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nkp
      real(8) qp(3,nkp),sf(2,nkp)
C Local parameters
      integer ifi,fopnx,iq
      real(8) sumev(2)

      ifi = fopnx('sfint',100,16+8+2+0,-1)
      write(ifi,901)
      sumev = 0
      do iq = 1, nkp
        write(ifi,902) qp(:,iq),sf(:,iq)
        sumev = sumev + sf(:,iq)
      enddo
      call fclose(ifi)

      sumev = sumev/nkp
      write(6,903) sumev(1)+sumev(2)

 901  format('#   Qx       Qy        Qz         sfup           sfdn')
 902  format(3f10.5,2f14.6)
 903  format(' sumev from integrated SF:',F12.6)
      end

!      end module specfun
