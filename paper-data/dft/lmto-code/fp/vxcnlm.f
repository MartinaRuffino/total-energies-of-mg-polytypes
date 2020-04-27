      subroutine vxcnlm(lxcfnl,lxcg,nsp,k1,k2,k3,s_lat,smrho,
     .  rex,rec,rexc,rvx,rvc,rvxc,vavg,smvx,smvc,smvxc,smexc)
C- Gradient correction to smoothed rho(q) tabulated on a mesh
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc ng alat vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv kv
Cio    Passed to:  *
Ci Inputs
Ci   lxcfnl: if nonzero, full GGA to be evaluated by libxc
Ci         : if zero, GGA is added to an LDA calculated elsewhere
Ci         :          GGA is specified by lxcg
Ci   lxcg  : Specifies GGA, if lxcfnl is zero
Ci         :  0    LSDA
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Ci   nsp   : 2 for spin-polarized case, otherwise 1
Ci   k1..k3: dimensions of smrho,vnl for smooth mesh density
Ci   smrho :smooth density on uniform mesh
Co Outputs
Co   rex   :integral [smrho * ex[smrho]]
Co         :calculated only if lxcfnl>0
Co   rec   :integral [smrho * ec[smrho]]
Co         :calculated only if lxcfnl>0
Co   rexc  :integral [smrho * exc[smrho]]
Co         :If lxcfnl > 0, entire contribution to sm rexc
Co         :Otherwise rexc = nonlocal contribution only
Co   rvx   :integral [smrho * vx[smrho]]
Co         :calculated only if lxcfnl>0
Co   rvc   :integral [smrho * vc[smrho]]
Co         :calculated only if lxcfnl>0
Co   rvxc  :integral [smrho * smvxc[smrho]]
Co         :If lxcfnl > 0, entire contribution to sm rvxc
Co         :Otherwise rvxc = nonlocal contribution only
Co   vavg  :average of smvxc
Co         :If lxcfnl = 0, contribution only
Co   smvx  :exchange potential on uniform mesh
Co         :calculated only if lxcfnl>0
Co   smvc  :correlation potential on uniform mesh
Co         :calculated only if lxcfnl>0
Co   smvxc :XC potential on uniform mesh added to smvxc
Co         :If lxcfnl > 0 potential copied to smvxc
Cl Local variables
Cl   ... Case lxcnl nonzero, and GGA specified by it
Cl   agr(*,1)  : |grad rhop|^2 or |grad rho|^2 if nsp=1
Cl   agr(*,2)  : |grad rhom|^2 (nsp=2)
Cl   agr(*,3)  : grad rhop . grad rhom
Cl   ggr(*,1)  : Laplacian of rhop (total rho if nsp=1)
Cl   ggr(*,2)  : Laplacian of rhom (nsp=2)
Cl   ... Case lxcnl = 0, and GGA specified by lxcg
Cl   agr(*,1)  : |grad rhop| or |grad rho| if nsp=1
Cl   agr(*,2)  : |grad rhom| (nsp=2)
Cl   agr(*,k)  : |grad total rho|. k=3 for nsp=2; else k=1
Cl   agr(*,4)  : grad rho+ . grad rho- (only for Langreth-Mehl-Hu)
Cl   ggr(*,1)  : Laplacian of rhop (total rho if nsp=1)
Cl   ggr(*,2)  : Laplacian of rhom (nsp=2)
Cl   gagr(*,k) : (grad rho).(grad |grad rho|)
Cl   gagr(*,1) : (grad rhop).(grad |grad rhop|) (total rho if nsp=1)
Cl   gagr(*,2) : (grad rhom).(grad |grad rhom|) (nsp=2)
Cl   gagr(*,k) : (grad rho).(grad |grad rho|). k=3 for nsp=2; else k=1
Cr Remarks
Cr
Cu Updates
Cu   17 Dec 16 (Jerome Jackson) return to analytic gradient scheme for
Cu             GGA potential, using grfmsh
Cu   24 Nov 16 (Jerome Jackson) correction to GGA potential (de/dgradient term included)
Cu             requires additional argument to xc_libxc call
Cu   09 Dec 13 First cut at using libxc functional for XC potential
Cu             Modified Argument list.
Cu   10 Nov 11 Begin migration to f90 structures
Cu   06 Apr 09 Adapted from vxcnlp.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lxcfnl,lxcg,k1,k2,k3,nsp
      double precision rexc(2),rvxc(2),vavg(2)
      double complex   smrho(k1,k2,k3,nsp)
      double complex   smvx(k1,k2,k3,nsp),smvc(k1,k2,k3,nsp)
      double complex   smvxc(k1,k2,k3,nsp),smexc(k1,k2,k3)
C ... For structures
!       include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8),allocatable :: ex(:,:,:),ec(:,:,:)
      real(8),allocatable :: vx(:,:,:,:),vc(:,:,:,:)
      real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
      real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:)
      complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
      double precision, allocatable:: dxdg2(:,:)
      double complex, allocatable :: dxdg2cplx(:),graddxdg2(:,:,:)
C ... Local parameters
      integer :: ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3)
      integer :: nx,px,ny,py,nz,pz,indexn,indexp
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision :: alat,vol,xx,tmp,plat(3,3)
      double precision :: rex(2),rec(2),rvx(2),rvc(2)
      double precision :: rexcnl(2),rmunl(2),vavgnl(2)
      double precision :: grad(nsp*2-1,3)

C --- Setup and memory allocation ---
      call tcn('vxcnlm')
      call getpr(ipr)
      ngabc = s_lat%nabc
      ng = s_lat%ng
      alat = s_lat%alat
      plat(:,:) = s_lat%plat(:,:)
      vol = s_lat%vol
      np = k1*k2*k3
      allocate(agr(np,3*nsp-2),gagr(np,2*nsp-1),ggr(np,nsp))
      allocate(zgrho(np,3,nsp),zggrho(np,nsp))

C --- Grad rho_i and Laplacian rho_i (complex) ---
      do  i = 1, nsp
C       call zprm3('smrho(isp=%i)',i,smrho(1,1,1,i),n1,n2,n3)
        call grfmsh(601,alat,ng,s_lat%gv,s_lat%kv,k1,k2,k3,n1,n2,n3,
     .    smrho(1,1,1,i),zgrho(1,1,i),zggrho(1,i))
C       call zprm3('gradx smrho(isp=%i)',i,zgrho(1,1,i),n1,n2,n3)
C       call zprm3('grady smrho(isp=%i)',i,zgrho(1,2,i),n1,n2,n3)
C       call zprm3('gradz smrho(isp=%i)',i,zgrho(1,3,i),n1,n2,n3)
C       call zprm3('lap smrho(isp=%i)',i,zggrho(1,i),n1,n2,n3)

C       call fftgrad(n1,n2,n3,alat,plat,smrho(1,1,1,i),zgrho(1,1,i))
C       call fftlap(n1,n2,n3,alat,plat,smrho(1,1,1,i),zggrho(1,i))
C       call zprm3('gradx smrho(isp=%i)',i,zgrho(1,1,i),n1,n2,n3)
C       call zprm3('grady smrho(isp=%i)',i,zgrho(1,2,i),n1,n2,n3)
C       call zprm3('gradz smrho(isp=%i)',i,zgrho(1,3,i),n1,n2,n3)
C       call zprm3('lap smrho(isp=%i)',i,zggrho(1,i),n1,n2,n3)
      enddo

C --- Setup for libxc call ---
      if (lxcfnl /= 0) then
C ... agr_i : |grad rho_i|^2, i=1,2 and agr(3) = grad rho+ . grad rho-
      do  i = 1, nsp
        do  ip = 1, np
          agr(ip,i) = dble(zgrho(ip,1,i))**2 +
     .                dble(zgrho(ip,2,i))**2 +
     .                dble(zgrho(ip,3,i))**2
          ggr(ip,i) = dble(zggrho(ip,i))
        enddo
      enddo
      if (nsp == 2) then
        do  ip = 1, np
          agr(ip,3) = dble(zgrho(ip,1,1)*zgrho(ip,1,2)) +
     .                dble(zgrho(ip,2,1)*zgrho(ip,2,2)) +
     .                dble(zgrho(ip,3,1)*zgrho(ip,3,2))
        enddo
      endif

      allocate(rho(np,nsp))
      do  i = 1, nsp
        call dcopy(np,smrho(1,1,1,i),2,rho(1,i),1)
      enddo
      allocate(vx(k1,k2,k3,nsp),vc(k1,k2,k3,nsp))
      allocate(ex(k1,k2,k3),ec(k1,k2,k3))
      allocate(dxdg2(2*nsp-1,np),graddxdg2(np,3,2*nsp-1))

C     Cannot poke result directly into smvxc, since smvxc is complex
      call xc_libxc(np,nsp,lxcfnl,rho,rho(1,nsp),xx,xx,xx,agr,
     .  xx,agr(1,nsp),agr(1,2*nsp-1),xx,xx,xx,ex,ec,vx,vc,
     .  vx(1,1,1,nsp),vc(1,1,1,nsp),xx,xx,xx,dxdg2)

CJJ, require d(exc)/d(gradient rho^2) from libxc call
C see, eg, equation (14) of libxc paper
C Marques, Oliveira, Burnus, CPC 183 (2012) 2272--2281
C
C dxdg2 corresponds to the libxc variable vsigma
C
      ! evaluate grad(dexc/d(grad rho)) via chain rule from grad(dexc/dsigma_i)
      ! with sigma_uu=grad(rho_u).grad(rho_u), similarly for ud,dd

      allocate(dxdg2cplx(np))
      do i = 1, 2*nsp-1
        dxdg2cplx(:) = dxdg2(i,:)!real part
        call grfmsh(201,alat,ng,s_lat%gv,s_lat%kv,k1,k2,k3,n1,n2,n3,
     .    dxdg2cplx,graddxdg2(1,1,i),xx)
C       call fftgrad(n1,n2,n3,alat,plat,dxdg2cplx,graddxdg2(1,1,i))
      end do
      deallocate(dxdg2cplx)

      ip=1
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        !spin 1
        tmp=0.d0

        do i=1,3 ! components
          tmp=tmp+2*dble(graddxdg2(ip,i,1))*dble(zgrho(ip,i,1))
          if(nsp==2)tmp=tmp+dble(graddxdg2(ip,i,2))*dble(zgrho(ip,i,2))
        end do
        tmp=tmp+2*dxdg2(1,ip)*dble(zggrho(ip,1)) ! imag(zggrho) anyway zero
        if(nsp==2)tmp=tmp+dxdg2(2,ip)*dble(zggrho(ip,2))
        vx(i1,i2,i3,1)=vx(i1,i2,i3,1)-tmp

        !spin 2
        if(nsp==2)then
          tmp=0.d0
          do i=1,3 ! components
            tmp=tmp+2*dble(graddxdg2(ip,i,3))*dble(zgrho(ip,i,2))+dble(graddxdg2(ip,i,2))*dble(zgrho(ip,i,1))
          end do
          tmp=tmp+2*dxdg2(3,ip)*dble(zggrho(ip,2))+dxdg2(2,ip)*dble(zggrho(ip,1))
          vx(i1,i2,i3,2)=vx(i1,i2,i3,2)-tmp
        end if

        ip=ip+1
      enddo
      enddo
      enddo
      deallocate(zgrho,zggrho,dxdg2,graddxdg2)
CJJ

C --- Integrated quantities ---
      do  i = 1, nsp
        rex(i) = 0
        rec(i) = 0
        rexc(i) = 0
        rvxc(i) = 0
        vavg(i) = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          smvx(i1,i2,i3,i) = vx(i1,i2,i3,i)
          smvc(i1,i2,i3,i) = vc(i1,i2,i3,i)
          smvxc(i1,i2,i3,i) = vx(i1,i2,i3,i) + vc(i1,i2,i3,i)
          smexc(i1,i2,i3) = ex(i1,i2,i3) + ec(i1,i2,i3)
          vavg(i) = vavg(i) + smvxc(i1,i2,i3,i)
          rex(i) = rex(i) + smrho(i1,i2,i3,i)*ex(i1,i2,i3)
          rec(i) = rec(i) + smrho(i1,i2,i3,i)*ec(i1,i2,i3)
          rvx(i) = rex(i) + smrho(i1,i2,i3,i)*vx(i1,i2,i3,i)
          rvc(i) = rec(i) + smrho(i1,i2,i3,i)*vc(i1,i2,i3,i)
          rexc(i) = rexc(i) + rex(i) + rec(i)
          rvxc(i)  = rvxc(i)  + rvx(i) + rvc(i)
        enddo
        enddo
        enddo
        vavg(i) = vavg(i)/(n1*n2*n3)
        rex(i)  = rex(i)*vol/(n1*n2*n3)
        rec(i)  = rec(i)*vol/(n1*n2*n3)
        rvx(i)  = rvx(i)*vol/(n1*n2*n3)
        rvc(i)  = rvc(i)*vol/(n1*n2*n3)
        rexc(i) = rex(i) + rec(i)
        rvxc(i)  = rvx(i) + rvc(i)
      enddo

      deallocate(ex,ec,vx,vc,ggr,agr,gagr,rho)

C --- Setup for other GGA calls --
      else

      allocate(gzgrho(np,3,2*nsp-1))

C ... agr_i : |grad rho_i|, i=1,2 and agr_i(3) : |grad rho|
C     and ggr_i = lap rho_i.  Also agr(4) : grad rho+ . grad rho-
      do  i = 1, nsp
        do  ip = 1, np
          agr(ip,i) = dsqrt(dble(zgrho(ip,1,i))**2 +
     .                      dble(zgrho(ip,2,i))**2 +
     .                      dble(zgrho(ip,3,i))**2)
          ggr(ip,i) = dble(zggrho(ip,i))
        enddo
      enddo
      if (nsp == 2) then
        do  ip = 1, np
          agr(ip,3) = dsqrt(dble(zgrho(ip,1,1)+zgrho(ip,1,2))**2 +
     .                      dble(zgrho(ip,2,1)+zgrho(ip,2,2))**2 +
     .                      dble(zgrho(ip,3,1)+zgrho(ip,3,2))**2)
          agr(ip,4) =       dble(zgrho(ip,1,1)*zgrho(ip,1,2)) +
     .                      dble(zgrho(ip,2,1)*zgrho(ip,2,2)) +
     .                      dble(zgrho(ip,3,1)*zgrho(ip,3,2))
        enddo
      endif
C     do  i = 1, 3*nsp-2
C       call prm3('|grad rho(isp=%i)|',i,agr(1,i),n1,n2,n3)
C     enddo

C ... gzgrho (complex) : grad |grad rho_i|, i=1,2,3 (see above for i=3)
C     Use zggrho as complex work array
      do  i = 1, 2*nsp-1
        call dpzero(zggrho,np*2)
        call dcopy(np,agr(1,i),1,zggrho,2)
C       call zprm3('|grad rho_i|',0,zggrho(1,i),n1,n2,n3)
        call grfmsh(201,alat,ng,s_lat%gv,s_lat%kv,k1,k2,k3,n1,n2,n3,
     .    zggrho,gzgrho(1,1,i),xx)
C       call zprm3('gradx |grad rho_%i|',i,gzgrho(1,1,i),n1,n2,n3)

C        call fftgrad(n1,n2,n3,alat,plat,zggrho,gzgrho(1,1,i))
C        call zprm3('gradx |grad rho_%i|',i,gzgrho(1,1,i),n1,n2,n3)
C        call zprm3('grady |grad rho_%i|',i,gzgrho(1,2,i),n1,n2,n3)
C        call zprm3('gradz |grad rho_%i|',i,gzgrho(1,3,i),n1,n2,n3)
      enddo
      deallocate(zggrho)

C ... gagr : grad rho_i . grad |grad rho_i|, i=1,2,3 (see above for i=3)
      do  i = 1, nsp
        do  ip = 1, np
          gagr(ip,i) =
     .      dble(zgrho(ip,1,i))*dble(gzgrho(ip,1,i)) +
     .      dble(zgrho(ip,2,i))*dble(gzgrho(ip,2,i)) +
     .      dble(zgrho(ip,3,i))*dble(gzgrho(ip,3,i))
        enddo
C       call prm3('grad rho . grad |grad rho_%i|',i,gagr(1,i),n1,n2,n3)
      enddo
      if (nsp == 2) then
        do  ip = 1, np
          gagr(ip,3) =
     .      dble(zgrho(ip,1,1)+zgrho(ip,1,2))*dble(gzgrho(ip,1,3)) +
     .      dble(zgrho(ip,2,1)+zgrho(ip,2,2))*dble(gzgrho(ip,2,3)) +
     .      dble(zgrho(ip,3,1)+zgrho(ip,3,2))*dble(gzgrho(ip,3,3))
        enddo
C       call prm3('grad rho . grad |grad rho_%i|',3,gagr(1,3),n1,n2,n3)
      endif

      deallocate(zgrho,gzgrho)

C --- Nonlocal potential for all points  ---
      allocate(vnl(k1,k2,k3,nsp),enl(k1,k2,k3),rho(np,nsp))
      call dpzero(vnl,np*nsp)
      call dpzero(enl,np)
      do  i = 1, nsp
        call dcopy(np,smrho(1,1,1,i),2,rho(1,i),1)
C       call zprm3('smrho_%i',i,smrho(1,1,1,i),n1,n2,n3)
C       call prm3 ('rho_%i',i,rho(1,i),n1,n2,n3)
C       call prm3 ('lap-rho_%i',i,ggr(1,i),n1,n2,n3)
      enddo

      if (lxcg > 2) then
        i = 1 ; if (nsp == 2) i = 4 ! To enable vxcnlm to pass bounds check
        call vxcgga(lxcg,np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .    ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,i),
     .    gagr(1,2*nsp-1),gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
      else
        lcut = 1
        if (lcut == 1) then
        call vxnlcc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .      ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1),
     .      gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
        else
          call vxnloc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .      ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1),
     .      gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
        endif
      endif
C      call prm3('enl',i,enl,n1,n2,n3)
C      do  i = 1, nsp
C        call prm3('vnl(isp=%i)',i,vnl(1,1,1,i),n1,n2,n3)
C      enddo

C --- Make nonlocal rexc, rvxc ---
      do  i = 1, nsp
        rexcnl(i) = 0
        rmunl(i) = 0
        vavgnl(i) = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rexcnl(i) = rexcnl(i) + dble(smrho(i1,i2,i3,i))*enl(i1,i2,i3)
          rmunl(i) = rmunl(i) + dble(smrho(i1,i2,i3,i))*vnl(i1,i2,i3,i)
          smvxc(i1,i2,i3,i) = smvxc(i1,i2,i3,i) + vnl(i1,i2,i3,i)
          vavgnl(i) = vavgnl(i) + vnl(i1,i2,i3,i)
        enddo
        enddo
        enddo
        rexc(i) = rexc(i) + rexcnl(i)*vol/(n1*n2*n3)
        rvxc(i) = rvxc(i) + rmunl(i)*vol/(n1*n2*n3)
        vavg(i) = vavg(i) + vavgnl(i)/(n1*n2*n3)
      enddo
      deallocate(rho,agr,gagr,ggr,enl,vnl)

      endif

      call tcx('vxcnlm')
      end
