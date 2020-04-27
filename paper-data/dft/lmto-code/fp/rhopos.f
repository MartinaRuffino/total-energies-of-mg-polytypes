      subroutine rhopos(mode,s_site,s_spec,s_lat,s_rhat,smrho)
C- Render density everywhere positive
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Cio    Passed to: siterhopos
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: a nr rmt lmxl
Co     Stored:    *
Cio    Passed to: siterhopos
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: nabc vol
Co     Stored:    *
Cio    Passed to: *
Cio  s_rhat
Ci     Elts read: rho1 rho2
Co     Stored:    Site densities set to zero where negative.
Co                Difference is put back in one-center expansion
Cio    Passed to: siterhopos
Ci Inputs
Ci   mode  :0 attempt to render rho1 and rho2 > 0 independently
Ci         :1 attempt to render rho1 and rho2 > 0 by shifting rho1 and rho2 in tandem
Ci Inputs/Outputs
Cio  smrho :smooth density on uniform mesh
Cio        :negative points of smrho are set to zero
Cio        :smrho is scaled so that net change in charge is zero
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   19 Nov 14 New mode
Cu   02 Aug 13 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode
      double complex smrho(*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_rhat)::  s_rhat(*)
C ... Local parameters
      integer ngabc(3),n1,n2,n3,k1,k2,k3,nsp,nglob
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision vol,dqm,dq1,dq2,qm,fac
      double precision tol ! negligible density correction
      parameter(tol=1d-6)
C ... External calls
      external daxpy,dcopy,dgemm,dmcpy,dpzero,dscal,fftz30,fpiint,getpr,
     .         info0,info2,info5,info8,mshint,prhoposx,radmsh,radwgt,
     .         ropyln,rxi,siterhopos,smrpos

      call info0(30,1,-1,' rhopos: try to render density everywhere positive ...')

C ... Site densities
      call siterhopos(mode,s_site,s_spec,s_rhat,dq1,dq2)

C ... Mesh density
      ngabc = s_lat%nabc
      vol = s_lat%vol
      call fftz30(n1,n2,n3,k1,k2,k3)
      nsp = nglob('nsp')
      call smrpos(smrho,vol,k1,k2,k3,n1,n2,n3,nsp,dqm)

C ... No adjustments; nothing to do
      if (abs(dqm)+dabs(dq1)+dabs(dq2) < tol) then
        call info0(30,0,0,' nothing to change')
        return
      endif

C ... Scale mesh density to ensure system charge neutrality
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,qm,fac)
      fac = 1 - (dqm+dq1-dq2)/qm
      call info5(30,0,0,' rhopos: scale mesh density by '//
     .    '1-(%,3;3e+%,3;3e-%,3;3e)/%;5F = %,6;6d',
     .  dqm,dq1,dq2,qm,fac)

C     call zprm3('smrho before scaling',0,smrho,k1,k2,k3)
      call dscal(2*k1*k2*k3*nsp,fac,smrho,1)
C     call zprm3('smrho after scaling',0,smrho,k1,k2,k3)
C     call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,dqm,fac)

      end
      subroutine smrpos(smrho,vol,k1,k2,k3,n1,n2,n3,nsp,dqm)
C- Make smrho real and positive
      implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,nsp
      double precision vol,dqm
      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer stdo,nglob,i1,i2,i3,nneg,isp
      double precision rmin,xx,qs(2),qsp(2)
      double precision tol ! 'small' density
      parameter(tol=1d-15)

      stdo = nglob('stdo')
      nneg = 0
      rmin = 999
      qs = 0 ; qsp = 0
      do  isp = 1, nsp
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              xx = dble(smrho(i1,i2,i3,isp))
              rmin = min(rmin,xx)
              qs(isp) = qs(isp) + xx
              if (xx < 0) then
                qsp(isp) = qsp(isp) - xx
                nneg = nneg+1
                smrho(i1,i2,i3,isp) = tol
              endif
            enddo
          enddo
        enddo
      enddo

      dqm = 0
      if (nneg > 0) then
        qs(1:nsp) = qs(1:nsp)*vol/(n1*n2*n3)
        qsp(1:nsp) = qsp(1:nsp)*vol/(n1*n2*n3)
        call info5(30,1,0,
     .    ' rhopos: mesh density negative at %i points.  '//
     .    'rhomin=%,3;3e  dQ=%,3;3e%?#n==2#  dQ(2)=%,3;3e##',
     .    nneg,rmin,qsp(1),nsp,qsp(2))
        dqm = qsp(1)+qsp(2)
      endif

      end

      subroutine siterhopos(mode,s_site,s_spec,s_rhat,dq1,dq2)
C- Render local site densities positive
      use structures
      implicit none
C ... Passed parameters
      integer mode
      double precision dq1,dq2
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated arrays
      real(8), allocatable :: rofi(:),rwgt(:)
C ... Local parameters
      integer nbas,nglob,ib,is,nr,lmxl,intopt,nlml,nsp,p1
      double precision a,rmt,dq1s,dq2s

      nbas = nglob('nbas')
      nsp = nglob('nsp')

      dq1 = 0 ; dq2 = 0; p1 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        lmxl = s_spec(is)%lmxl
        if (lmxl == -1) cycle
        nlml = (lmxl+1)**2
        intopt = 10*nglob('lrquad')
        allocate(rofi(nr),rwgt(nr))
        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)

        if (mode == 0) then
          call prhopos1(ib,1,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho1,dq1s)
          call prhopos1(ib,2,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho2,dq2s)
        else
          call prhopos2(ib,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho1,s_rhat(ib)%rho2,dq1s,dq2s)
        endif

        dq1 = dq1 + dq1s
        dq2 = dq2 + dq2s

        if (mode == 0) then
          call prhopos1(ib,1,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho1,dq1s)
          call prhopos1(ib,2,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho2,dq2s)
        else
          call prhopos2(ib,p1,rofi,nr,rwgt,nlml,nsp,s_rhat(ib)%rho1,s_rhat(ib)%rho2,dq1s,dq2s)
        endif

        deallocate(rofi,rwgt)

      enddo

      end
      subroutine prhopos1(ib,it,p1,ri,nr,rwgt,nlm,nsp,rhol,dq)
C- Render a local density positive everywhere in space
C ----------------------------------------------------------------------
Ci Inputs
Ci   ib    :site index (for printout only)
Ci   it    :1 for true rho, 2 for sm rho (for printout only)
Ci   p     :0 => head not yet printed
Ci         :1 => header printed out
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rhol  :full charge density * r**2
Co Outputs
Cl Local variables
Cl   yl    :YL(rp)
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib,it,nr,nlm,nsp,p1
      double precision ri(nr),rwgt(nr),rhol(nr,nlm,nsp),dq
C ... Dynamically allocated arrays
      real(8), allocatable :: yl(:,:)
      real(8), allocatable :: rp(:,:,:),rps(:,:,:),rl(:,:,:)
C ... Local parameters
      integer nnn  ! Maximum number of angular mesh points
      integer maxit  ! Maximum number of iterations
      double precision tol ! 'small' density
      parameter(nnn=122,tol=1d-12,maxit=5)
      integer stdo,ipr,ll,lmax,np,nph,nth,nxl(0:7),nglob
      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn),fpi,weight
      character cct
      integer ip,ir,i,n1,ineg(2),iter,inegi(2)
      double precision rpneg,rpnegi,qsi(2),qspi(2),dqi(2),q0
      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

      stdo = nglob('stdo')
      call getpr(ipr)
      fpi = 16d0*datan(1d0)
C ... Create angular integration mesh
      lmax = ll(nlm)
      if (lmax > 6) then
        nth = 2*lmax+2
        nph = 0
      else
        nth = nxl(lmax)
C       nth = nxl(7)
        nph = 0
      endif
      call fpiint(nth,nph,np,p,wp)
      dq = 0

      n1 = nr*np
C      if (ipr >= 30) write (stdo,1) nth,nph,np,nr
C    1 format(' mesh:   nth,nph=',2I4,'   gives',i4,'  angular points,',
C     .    '   nrad=',i4)
      if (np > nnn) call rxi('rhopos: increase nnn, need',np)

C ... Scale rhol to true density
      allocate(rl(nr,nlm,nsp))
      call dcopy(nr*nlm*nsp,rhol,1,rl,1)
      call prhoposx(nr,nlm,nsp,ri,rl,1)

C ... Spherical harmonics (unit radius) ... to be overwritten by Yl * wp
      call dmcpy(p,1,3,p2,np,1,np,3)
      allocate(yl(np,nlm))
      call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,yl,r2)

C ... Setup
      allocate(rp(nr,np,nsp),rps(nr,np,nsp))
      call dpzero(rp,n1*nsp)
      ineg = 0; rpneg = 0; dqi = 0

C --- Generate density point-wise through sphere, render > 0 ---
      do  iter = 1, maxit
      inegi = 0; qsi = 0; qspi = 0; rpnegi = 0
      call fyl2p(nr,nlm,nsp,np,yl,rl,rps)
      call dcopy(n1*nsp,rps,1,rp,1) ! Save original
      do  i = 1, nsp

C   ... Counts number of points with negative density; reset to pos.
        do  ip = 1, np
          do  ir = 1, nr
            weight = (ri(ir)**2*rwgt(ir))*wp(ip)
            qsi(i)  = qsi(i) + rps(ir,ip,i)*weight
            if (rps(ir,ip,i) < 0d0) then
              rpneg = min(rpneg,rps(ir,ip,i))
              rpnegi = min(rpnegi,rps(ir,ip,i))
              inegi(i) = inegi(i) + 1
              rps(ir,ip,i) = tol
C             if (rps(ir,ip,i) < -.002) then
C               print *, 'hi',ir,ip,rps(ir,ip,i)
C             endif
            endif
            qspi(i)  = qspi(i) + rps(ir,ip,i)*weight
          enddo
        enddo
      enddo

      ineg(1) = max(ineg(1),inegi(1))
      ineg(2) = max(ineg(2),inegi(2))
      if (ineg(1)+ineg(2) == 0) goto 99 ! no negative points

      if (iter == 1) q0 = qsi(1) + qsi(2)
      dqi(1) = dqi(1) + qspi(1)-qsi(1)
      dqi(2) = dqi(2) + qspi(2)-qsi(2)

C     rps <- change in rp
      call daxpy(n1*nsp,-1d0,rp,1,rps,1)
C     Add change to rl
      call fp2yl(nr,nlm,nsp,np,wp,rps,yl,1d0,rl)

      if (iter == maxit .or. inegi(1)+inegi(2) == 0) then

      if (p1 == 0) then
        call info2(30,1,0,'  ib%6fnneg   min val%4flast'//
     .  '%7fdQ%6fdQ/Q%?#n==2#     spin1    spin2##',nsp,0)
        p1 = 1
      endif

C      call info8(30,0,0,'%,4i%?#n==1#t#s#%,8i  '//
C     .  '%:-1,3;3e %;12,7D %:-1,3;3e%?#n==2# %:-1,3;3e##',ib,it,
C     .  inegi(1)+inegi(2),rpneg,(dqi(1)+dqi(2))/q0,dqi(1),nsp,dqi(2))
C
C      cct = 's'; if (it == 1) cct = 't'
C      call info8(30,0,0,'%,4i%?#n==1#t#s#%,8i  '//
C     .  '%:-1,3;3e%:-2;8F %:-1;6F%?#n==2# %:-1;6F##',ib,it,
C     .  inegi(1)+inegi(2),rpneg,(dqi(1)+dqi(2))/q0,dqi(1),nsp,dqi(2))

      cct = 's'; if (it == 1) cct = 't'
C     dq = qspi(1)-qsi(1) + qspi(2)-qsi(2)
      if (nsp == 1) then
        call info8(30,0,0,'%,4i '//cct//'%,8i%:-3,3;3e%:-3,3;3e'//
     .    '%:-3,3;3e%:-3;6F',ib,ineg(1)+ineg(2),rpneg,
     .    rpnegi,(dqi(1)+dqi(2)),(dqi(1)+dqi(2))/q0,0,0)
      else
        call info8(30,0,0,'%,4i '//cct//'%,8i%:-3,3;3e%:-3,3;3e'//
     .    '%:-3,3;3e%:-3;6F%:-3;6F%:-3;6F',ib,ineg(1)+ineg(2),rpneg,
     .    rpnegi,(dqi(1)+dqi(2)),(dqi(1)+dqi(2))/q0,dqi(1),dqi(2))
      endif

      if (inegi(1)+inegi(2) == 0) exit

      endif ! Last iteration

      enddo ! Loop over iterations
C     Net change in charge
      dq = (dqi(1)+dqi(2))

C --- debugging : confirm that rl -> positive rps.  Comment out daxpy above ---
C      call xxx(ib,it,ri,nr,rwgt,nlm,nsp,rl)
C      stop

C ... Undo the r^2 scaling; add to rhol
      call prhoposx(nr,nlm,nsp,ri,rl,0)
      call dcopy(nr*nlm*nsp,rl,1,rhol,1)

   99 continue
      deallocate(rl,yl,rp,rps)
      end

      subroutine prhopos2(ib,p1,ri,nr,rwgt,nlm,nsp,rhol1,rhol2,dq1,dq2)
C- Render a a pair of local densities positive everywhere in space
C ----------------------------------------------------------------------
Ci Inputs
Ci   ib    :site index (for printout only)
Ci   p     :0 => head not yet printed
Ci         :1 => header printed out
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rhol1 :true full charge density * r**2
Ci   rhol2 :smooth charge density * r**2
Co Outputs
Cl Local variables
Cl   yl    :YL(rp)
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib,nr,nlm,nsp,p1
      double precision ri(nr),rwgt(nr),rhol1(nr,nlm,nsp),dq1,rhol2(nr,nlm,nsp),dq2
C ... Dynamically allocated arrays
      real(8), allocatable :: yl(:,:)
      real(8), allocatable :: rp1(:,:,:),rps1(:,:,:),rl1(:,:,:)
      real(8), allocatable :: rp2(:,:,:),rps2(:,:,:),rl2(:,:,:)
C ... Local parameters
      logical lneg1,lneg2
      integer nnn  ! Maximum number of angular mesh points
      integer maxit  ! Maximum number of iterations
      double precision tol ! 'small' density
      parameter(nnn=122,tol=1d-12,maxit=5)
      integer stdo,ipr,ll,lmax,np,nph,nth,nxl(0:7),nglob
      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn),fpi,weight
      integer ip,ir,i,n1,iter,ineg(2,2),inegi(2,2),iv(2)
      double precision rpneg1,rpnegi1,qsi1(2),qspi1(2),dqi1(2),q01,add1
      double precision rpneg2,rpnegi2,qsi2(2),qspi2(2),dqi2(2),q02,add2
      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

      stdo = nglob('stdo')
      call getpr(ipr)
      fpi = 16d0*datan(1d0)
      dq1 = 0; dq2 = 0
C ... Create angular integration mesh
      lmax = ll(nlm)
      if (lmax > 6) then
        nth = 2*lmax+2
        nph = 0
      else
        nth = nxl(lmax)
C       nth = nxl(7)
        nph = 0
      endif
      call fpiint(nth,nph,np,p,wp)
      dq1 = 0

      n1 = nr*np
      if (np > nnn) call rxi('rhopos: increase nnn, need',np)

C ... Scale rhol1,rhol2 to true density
      allocate(rl1(nr,nlm,nsp),rl2(nr,nlm,nsp))
      call dcopy(nr*nlm*nsp,rhol1,1,rl1,1)
      call prhoposx(nr,nlm,nsp,ri,rl1,1)
      call dcopy(nr*nlm*nsp,rhol2,1,rl2,1)
      call prhoposx(nr,nlm,nsp,ri,rl2,1)

C ... Spherical harmonics (unit radius) ... to be overwritten by Yl * wp
      call dmcpy(p,1,3,p2,np,1,np,3)
      allocate(yl(np,nlm))
      call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,yl,r2)

C --- Generate density point-wise through sphere, render > 0 ---
      allocate(rp1(nr,np,nsp),rps1(nr,np,nsp),rp2(nr,np,nsp),rps2(nr,np,nsp))
      call dpzero(rp1,n1*nsp); call dpzero(rp2,n1*nsp)
      ineg = 0; rpneg1 = 0; dqi1 = 0
      ineg = 0; rpneg2 = 0; dqi2 = 0
      do  iter = 1, maxit
      inegi = 0; qsi1 = 0; qspi1 = 0; rpnegi1 = 0
      inegi = 0; qsi2 = 0; qspi2 = 0; rpnegi2 = 0
      call fyl2p(nr,nlm,nsp,np,yl,rl1,rps1)
      call dcopy(n1*nsp,rps1,1,rp1,1) ! Save original
      call fyl2p(nr,nlm,nsp,np,yl,rl2,rps2)
      call dcopy(n1*nsp,rps2,1,rp2,1) ! Save original

      do  i = 1, nsp

C   ... Counts number of points with negative density; reset to pos.
        do  ip = 1, np
          do  ir = 1, nr
            weight = (ri(ir)**2*rwgt(ir))*wp(ip)
            qsi1(i)  = qsi1(i) + rps1(ir,ip,i)*weight
            qsi2(i)  = qsi2(i) + rps2(ir,ip,i)*weight
            lneg1 = rps1(ir,ip,i) < 0d0
            lneg2 = rps2(ir,ip,i) < 0d0
            add1 = 0; add2 = 0
            if (lneg1) then
              rpneg1 = min(rpneg1,rps1(ir,ip,i))
              rpnegi1 = min(rpnegi1,rps1(ir,ip,i))
              inegi(i,1) = inegi(i,1) + 1
              add1 = -rps1(ir,ip,i)*(1+tol)
            endif
            if (lneg2) then
              rpneg2 = min(rpneg2,rps2(ir,ip,i))
              rpnegi2 = min(rpnegi2,rps2(ir,ip,i))
C              if (rps2(ir,ip,i) < -.002) then
C                print *, 'hi',ir,ip,rps2(ir,ip,i)
C              endif
              inegi(i,2) = inegi(i,2) + 1
              add2 = -rps2(ir,ip,i)*(1+tol)
            endif
            rps1(ir,ip,i) = rps1(ir,ip,i) + max(add1,add2)
            rps2(ir,ip,i) = rps2(ir,ip,i) + max(add1,add2)
            qspi1(i)  = qspi1(i) + rps1(ir,ip,i)*weight
            qspi2(i)  = qspi2(i) + rps2(ir,ip,i)*weight
          enddo
        enddo
      enddo

      ineg(1,1) = max(ineg(1,1),inegi(1,1))
      ineg(2,1) = max(ineg(2,1),inegi(2,1))
      ineg(1,2) = max(ineg(1,2),inegi(1,2))
      ineg(2,2) = max(ineg(2,2),inegi(2,2))
      if (sum(ineg(:,:)) == 0) goto 99 ! no negative points

      if (iter == 1) then
        q01 = qsi1(1) + qsi1(2)
        q02 = qsi2(1) + qsi2(2)
      endif
      dqi1(1) = dqi1(1) + qspi1(1)-qsi1(1)
      dqi1(2) = dqi1(2) + qspi1(2)-qsi1(2)
      dqi2(1) = dqi2(1) + qspi2(1)-qsi2(1)
      dqi2(2) = dqi2(2) + qspi2(2)-qsi2(2)

C     rps(1,2) <- change in rp(1,2)
      call daxpy(n1*nsp,-1d0,rp1,1,rps1,1)
      call daxpy(n1*nsp,-1d0,rp2,1,rps2,1)
C     Add change to rl(1,2)
      call fp2yl(nr,nlm,nsp,np,wp,rps1,yl,1d0,rl1)
      call fp2yl(nr,nlm,nsp,np,wp,rps2,yl,1d0,rl2)

C      rl2 = 0
C      call fp2yl(nr,nlm,nsp,np,wp,rps2,yl,1d0,rl2)
C      call fyl2p(nr,nlm,nsp,np,yl,rl2,rps2)

      if (iter == maxit .or. sum(inegi(:,:)) == 0) then

        if (p1 == 0) then
          call info2(30,1,0,'  ib%3fnneg%9fmin val%4flast'//
     .      '%7fdQ%6fdQ/Q%?#n==2#     spin1    spin2##',nsp,0)
          p1 = 1
        endif

C       iv(1) = sum(inegi(:,1)); iv(2) = sum(inegi(:,2))
        iv(1) = sum(ineg(:,1)); iv(2) = sum(ineg(:,2))
        if (nsp == 1) then
          call info8(30,0,0,'%,4i  %2:1i %17p%:-3,3;3e%:-3,3;3e%:-3,3;3e%:-3;6F',
     .      ib,iv,min(rpneg1,rpneg2),min(rpnegi1,rpnegi2),
     .      (dqi1(1)+dqi1(2)),(dqi1(1)+dqi1(2))/q01,dqi1(1),dqi1(2))
        else
          call info8(30,0,0,'%,4i  %2:1i %17p%:-3,3;3e%:-3,3;3e%:-3,3;3e%:-3;6F%:-3;6F%:-3;6F',
     .      ib,iv,min(rpneg1,rpneg2),min(rpnegi1,rpnegi2),
     .      (dqi1(1)+dqi1(2)),(dqi1(1)+dqi1(2))/q01,dqi1(1),dqi1(2))

        endif

        if (sum(inegi(:,:)) == 0) exit ! Last iteration
      endif                     ! printout

      enddo ! Loop over iterations
C     Net change in charge
      dq1 = (dqi1(1)+dqi1(2))
      dq2 = (dqi2(1)+dqi2(2))

C --- debugging : confirm that rl1 -> positive rps1.  Comment out daxpy above ---
C      call xxx(ib,it,ri,nr,rwgt,nlm,nsp,rl1)
C      stop

C ... Undo the r^2 scaling; add to rhol(1,2)
      call prhoposx(nr,nlm,nsp,ri,rl1,0)
      call dcopy(nr*nlm*nsp,rl1,1,rhol1,1)
      call prhoposx(nr,nlm,nsp,ri,rl2,0)
      call dcopy(nr*nlm*nsp,rl2,1,rhol2,1)

   99 continue
      deallocate(yl,rl1,rp1,rps1,rl2,rp2,rps2)
      end

      subroutine prhoposx(nr,nlm,nsp,ri,rl,isgn)
C- Scales rho by r**2, or undoes scaling
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ri    :radial mesh points
Ci   isgn  :1, scale rl by 1/r**2
Ci         :else scale rl by r**2
Co Outputs
Co   rl   :rl is scaled by r**2 or 1/r**2, depending on isgn
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer nr,nlm,nsp,isgn,i,ilm,ir
      double precision rl(nr,nlm,nsp),ri(nr),rho2,rho3

C --- Scale rho by 1/r**2 ---
      if (isgn == 1) then
        do  i = 1, nsp
          do  ilm = 1, nlm
            rl(1,ilm,i) = 0d0
            do  ir = 2, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)/ri(ir)**2
            enddo
          enddo
        enddo
C  ...  Extrapolate rho to origin
        do  i = 1, nsp
          rho2 = rl(2,1,i)
          rho3 = rl(3,1,i)
          rl(1,1,i) = (rho2*ri(3)-rho3*ri(2))/(ri(3)-ri(2))
        enddo
      else
        do  i = 1, nsp
          do  ilm = 1, nlm
            do  ir = 1, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)*ri(ir)**2
            enddo
          enddo
        enddo
      endif
      end

C      subroutine xxx(ib,it,ri,nr,rwgt,nlm,nsp,rhol)
CC- Render local density positive everywhere in space
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ib    :site index (for printout only)
CCi   it    :1 for true rho, 2 for sm rho (for printout only)
CCi   ri    :radial mesh points
CCi   nr    :number of radial mesh points
CCi   rwgt  :radial mesh weights
CCi   nlm   :L-cutoff for density expansion
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi   rhol    :full charge density * r**2
CCo Outputs
CCl Local variables
CCl   yl    :YL(rp)
CCr Remarks
CCu Updates
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer ib,it,nr,nlm,nsp
C      double precision ri(nr),rwgt(nr),rhol(nr,nlm,nsp)
CC ... Dynamically allocated arrays
C      real(8), allocatable :: yl(:,:)
C      real(8), allocatable :: rp(:,:,:),rps(:,:,:),rl(:,:,:)
CC ... Local parameters
C      integer nnn  ! Maximum number of angular mesh points
C      double precision tol ! 'small' density
C      parameter(nnn=122,tol=1d-15)
C      integer stdo,ipr,ll,lmax,np,nph,nth,nxl(0:7),nglob
C      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn),fpi
CC ... Local parameters
C      integer ip,ir,i,n1,ineg(2)
C      double precision weight
C      double precision qs(2),qsp(2),rpneg
C      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/
C
C      stdo = nglob('stdo')
C      call getpr(ipr)
C      fpi = 16d0*datan(1d0)
CC ... Create angular integration mesh
C      lmax = ll(nlm)
C      if (lmax > 6) then
C        nth = 2*lmax+2
C        nph = 0
C      else
C        nth = nxl(lmax)
C        nph = 0
C      endif
C      call fpiint(nth,nph,np,p,wp)
C
C      n1 = nr*np
CC      if (ipr >= 30) write (stdo,1) nth,nph,np,nr
CC    1 format(' mesh:   nth,nph=',2I4,'   gives',i4,'  angular points,',
CC     .    '   nrad=',i4)
C      if (np > nnn) call rxi('rhopos: increase nnn, need',np)
C
CC ... Scale rhol to true density
C      allocate(rl(nr,nlm,nsp))
C      call dcopy(nr*nlm*nsp,rhol,1,rl,1)
C      call prhoposx(nr,nlm,nsp,ri,rl,1)
C
CC ... Spherical harmonics (unit radius) ... to be overwritten by Yl * wp
C      call dmcpy(p,1,3,p2,np,1,np,3)
C      allocate (yl(np,nlm))
C      call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,yl,r2)
C
CC --- Generate density point-wise through sphere, render > 0 ---
C      allocate (rp(nr,np,nsp),rps(nr,np,nsp))
C      call dpzero(rp,n1*nsp)
C      ineg = 0; rpneg = 0; qs = 0; qsp = 0
C      do  i = 1, nsp
C        call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0,
C     .    rps(1,1,i),nr)
C        call dcopy(n1,rps(1,1,i),1,rp(1,1,i),1) ! Save original
C
CC   ... Counts number of points with negative density; reset to pos.
C        do  ip = 1, np
C          do  ir = 1, nr
C            weight = (ri(ir)**2*rwgt(ir))*wp(ip)
C            qs(i)  = qs(i) + rps(ir,ip,i)*weight
C            if (rps(ir,ip,i) < 0d0) then
C              rpneg = min(rpneg,rps(ir,ip,i))
C              ineg(i) = ineg(i) + 1
C              rps(ir,ip,i) = tol
C            endif
C            qsp(i)  = qsp(i) + rps(ir,ip,i)*weight
C          enddo
C        enddo
C      enddo
C
C      call info2(30,1,0,'  ib%5fnneg  min val     Qnew/Qold'//
C     .  '    dQ%?#n==2#%7fdQ(2)##',nsp,0)
C
C      call info8(30,0,0,'%,4i%?#n==1#t#s#'//
C     .  '%,8i  %:-1,3;4e %;12,7D %:-1,3;4e'//
C     .  '%?#n==2# %:-1,3;4e##',ib,it,ineg(1)+ineg(2),rpneg,
C     .  (qs(1)+qs(2))/(qsp(1)+qsp(2)),qsp(1)-qs(1),nsp,qsp(2)-qs(2))
C
C      print *, qs
C      print *, qsp
C      end
