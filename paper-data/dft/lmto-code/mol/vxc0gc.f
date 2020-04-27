      subroutine vxc0gc(nr,nsp,rofi,rwgt,rho,vxc,exc,rep,rmu,lxcfun)
C- Gradient-corrected part of vxc and exc in a spherical potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   rho   :(true rho)*(4*pi*r**2)
Ci   lxcfun:specifies exchange-correlation functional
Ci         :1s digit sets local xc functional
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci         :100s digit sets gradient corrections
Ci         :  0    no gradient corrections
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Co Outputs
Co   vxc   :nonlocal XC potential
Co   exc   :nonlocal XC energy
Co   rep   :int rho * exc
Co   rmu   :int rho * vxc
Cl Local variables
Cr Remarks
Cu Updates
Cu   03 Jul 13 Replace f77 pointers with f90 ones
Cu   05 Apr 09  Use rwgt for mesh weights
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lxcfun,nr,nsp
      double precision rofi(nr),rwgt(nr),vxc(nr,nsp),rho(nr,nsp),
     .  rep(nsp),rmu(nsp),exc(nr)
C ... Dynamically allocated local arrays
      real(8), allocatable :: grh(:),ggrh(:),agrh(:),grgag(:)
C ... Local parameters
      integer nrmx
      parameter (nrmx=5001)
      double precision pi,rho2,rho3,ub4pi,rp(nrmx*2),rho0(2)
      integer i,ir,isp,lxcg

      pi = 4d0*datan(1d0)
      ub4pi = 1d0/(4d0*pi)
      call dpzero(vxc, nr*nsp)
      call dpzero(exc, nr)
      call dpzero(rep, nsp)
      call dpzero(rmu, nsp)
      lxcg = iabs(lxcfun)/100
      if (lxcg == 0) return

      if (nr > nrmx) call rx('vxc0gc: nr > nrmx')

C --- Make true rho ---
      do  isp = 1, nsp
        rho2 = rho(2,isp)/rofi(2)**2
        rho3 = rho(3,isp)/rofi(3)**2
        rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
        rp(1+nr*(isp-1)) = rho0(isp)
        do  ir = 2, nr
          rp(ir+nr*(isp-1)) = rho(ir,isp)*ub4pi/rofi(ir)**2
        enddo
      enddo

C      print *, 'rl,l=0'
C      call prmr(21,rofi,rp,1)

C --- Gradient correction ---
      allocate(grh(nrmx*nsp))
      allocate(ggrh(nrmx*nsp))
      allocate(agrh(nrmx*(3*nsp-2)))
      allocate(grgag(nrmx*(2*nsp-1)))
      call vxcgr2(lxcg,nr,nsp,nr,rofi,rp,grh,ggrh,agrh,grgag,exc,vxc)
      deallocate(grh,ggrh,agrh,grgag)
      do  i = 1, nsp
        rep(i) = 0d0
        rmu(i) = 0d0
        do  ir = 1, nr
          rep(i) = rep(i) + rwgt(ir)*rho(ir,i)*exc(ir)
          rmu(i) = rmu(i) + rwgt(ir)*rho(ir,i)*vxc(ir,i)
        enddo
      enddo

      end
