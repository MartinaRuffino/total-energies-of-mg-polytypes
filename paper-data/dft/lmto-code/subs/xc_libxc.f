C#define LIBXC
      subroutine xc_libxc(n,nsp,lxcf,rho,rho2,g2rho,g2up,g2dn,grrsqr,
     .  gup2,grrsq2,gupdn,tau,tauup,taudn,ex,ec,vx,vc,vx2,vc2,exc,vxc,
     .  vxc2,dxdg2)
C- Interface to calculate exchange-correlation functionals from libxc
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     : number of points
Ci   nsp   : 2 for spin-polarized case, otherwise 1
Ci   lxcf  : a compound of 2 digits, so that two functionals are evaluated
Ci         : lxcf/2**16 = index to the first functional in libxc
Ci         : mod(lxcf,2**16) = index to the 2nd functional in libxc
Ci         : Sign of lxcf is used as a switch.
Ci         : lxcf<0 => return exc,vxc.  Otherwise leave untouched
Ci   rho   : spin 1 density, or if not spin polarized, total spin density
Ci   rho2  : spin 2 density
Ci   g2rho : nabla rho.   Used only for meta-GGAs
Ci   g2up  : nabla rhoup. Used only for meta-GGAs
Ci   g2dn  : nabla rhodn. Used only for meta-GGAs
Ci   grrsqr: |grad rhoup|^2, or if not spin polarized, |grad rho|^2
Ci         :              Used for GGA and meta-GGA
Ci   grrsq2: |grad rhodn|^2 (spin polarized only)
Ci         :              Used for GGA and meta-GGA
Ci   gupdn : (grad rhoup).(grad rhodn) (spin polarized only)
Ci         :              Used for GGA and meta-GGA
Ci   tau   : kinetic energy density        Used only for meta-GGAs
Ci   tauup : spin1 kinetic energy density. Used only for meta-GGAs
Ci   taudn : spin2 kinetic energy density. Used only for meta-GGAs
Co Outputs
Co   ex    : exchange energy
Co         : (actually, energy density of 1st functional)
Co   ec    : correlation energy
Co         : (actually, energy density of 2nd functional)
Co   vx    : exchange potential (spin1 if spin polarized case)
Co         : (actually, potential or spin1 potential of 1st functional)
Co   vc    : correlation potential (spin1 if spin polarized case)
Co         : (actually, potential or spin1 potential of 2nd functional)
Co   vx2   : exchange potential, spin2.  Not touched if nsp=1.
Co         : (actually, spin2 potential of 1st functional)
Co   vc2   : correlation potential, spin2.  Not touched if nsp=1.
Co         : (actually, spin2 potential of 2nd functional)
Co   vxc   : vx + vc  (returned only if lxcf<0).
Co   exc   : ex + ec  (returned only if lxcf<0).
Cr Remarks
Cr   Compile with :  fc -I$(LIBXCPATH)/include -c xc_libxc.f
Cr   where $(LIBXCPATH)/lib is the directory containing libxc library.
Cb Bugs
Cb   When grrsqr=0, dxdg2 need not be 0.
Cb   However, libxc may return 0 when grrsqr falls below a cutoff.
Cb   Fix: set Use a small, nonzero number in place of grrsqr when grrsqr=0.
Cu Updates
Cu   17 Dec 16 (MvS) workaround libxc bug for GGA.
Cu   26 Nov 16 (Jerome Jackson) Copy out vsigma from libxc call for GGA
Cu   30 Nov 13 First implemented.  LDA and GGA only so far.
C ----------------------------------------------------------------------
C#ifdef LIBXC
      use xc_f90_lib_m
C#endif
      implicit none
C ... Passed parameters
      integer n,nsp,lxcf
      real(8) rho(n),rho2(n),g2rho(n),g2up(n),g2dn(n),
     .  grrsqr(n),gup2(n),grrsq2(n),gupdn(n),tau(n),tauup(n),taudn(n)
      real(8) exc(n),vxc(n),vxc2(n),dxdg2(2*nsp-1,n)
      real(8), target :: ex(n),vx(n),vx2(n),ec(n),vc(n),vc2(n)

C#ifndefC LIBXC
C      call rx('this density functional requires '//
C     .  'that you link the libxc library; sorry')
C#else
C ... Local arrays and pointers
      real(8), pointer :: pe(:),pv(:),pv2(:)
      type(xc_f90_pointer_t) xc_func,info
      real(8), allocatable :: r(:,:),vrho(:,:)
      real(8), allocatable :: sigma(:,:),vsigma(:)
C ... Local parameters
      integer i,k,nspin,id,xctype
C     Bug in libxc ... patch by bounding grrsqr(:) to tol from below
      real(8), parameter :: tol = 2d0**(-99) ! Slightly larger than 1d-30

      nspin=XC_UNPOLARIZED; if (nsp > 1) nspin=XC_POLARIZED

      do  k = 1, 2

        if (k == 1) id = iabs(lxcf)/2**16
        if (k == 2) id = mod(iabs(lxcf),2**16)

        if (id == 0 .and. k == 1) then  ! no first functional
          call dpzero(ex,n); call dpzero(vx,n)
          if (nsp == 2) call dpzero(vx2,n)
          cycle
        elseif (id == 0) then   ! no second functional
          call dpzero(ec,n); call dpzero(vc,n)
          if (nsp == 2) call dpzero(vc2,n)
          cycle
        endif

        xctype = xc_f90_family_from_id(id)  ! Determines whether LDA, GGA, or Meta-LDA
        select case(xctype)

C   ... LDA
        case(XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,info,id,nspin)
          if (k. eq. 1) then    ! First functional : poke into ex,vx,vx2
            pe => ex
            pv => vx
            if (nsp == 2) pv2 => vx2
          else                  ! Second functional : poke into ec,vc,vc2
            pe => ec
            pv => vc
            if (nsp == 2) pv2 => vc2
          endif
          if (nsp == 1) then
            call xc_f90_lda_exc_vxc(xc_func,n,rho(1),pe(1),pv(1))
            call dscal(n,2d0,pe,1); call dscal(n,2d0,pv,1)
          else
            allocate(r(2,n),vrho(2,n))
            r(1,:) = rho(:); r(2,:) = rho2(:)
            call xc_f90_lda_exc_vxc(xc_func,n,r(1,1),pe(1),vrho(1,1))
            pv(:) = 2*vrho(1,:); pv2(:) = 2*vrho(2,:); pe(:) = 2*pe(:) ! => Atomic Rydberg units
C           print *, pv(100),vx(100),vc(100)
            deallocate(r,vrho)
          endif
          call xc_f90_func_end(xc_func) ! destroy functional

C    .. GGA
        case(XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,info,id,nspin)
          if (k. eq. 1) then    ! First functional : poke into ex,vx,vx2
            dxdg2=0.d0 ! zero vsigma at first entry to function
            pe => ex
            pv => vx
            if (nsp == 2) pv2 => vx2
          else                  ! Second functional : poke into ec,vc,vc2
            pe => ec
            pv => vc
            if (nsp == 2) pv2 => vc2
          endif
          if (nsp == 1) then
C            call xc_f90_lda_exc_vxc(xc_func,n,rho(1),pe(1),pv(1))
            allocate(sigma(1,n),vsigma(n))
            do  i = 1, n
              sigma(1,i) = max(grrsqr(i),tol)
            enddo
            call xc_f90_gga_exc_vxc(xc_func,n,rho(1),sigma(1,1),pe(1),pv(1),vsigma(1))
            call daxpy(n,2.d0,vsigma,1,dxdg2,1)!JJ, copy out dxdg2; alpha=2 for hartree->ryd
            deallocate(sigma,vsigma)
            call dscal(n,2d0,pe,1); call dscal(n,2d0,pv,1)
          else
            allocate(r(2,n),vrho(2,n),sigma(3,n),vsigma(3*n))
            r(1,:) = rho(:); r(2,:) = rho2(:)
            do  i = 1, n
              sigma(1,i) = max(grrsqr(i),tol)
              sigma(2,i) = gupdn(i)
              sigma(3,i) = max(grrsq2(i),tol)
            enddo
            call dpzero(vsigma,size(vsigma))
C           call prmx('r',r,size(r),size(r),1)
C           call prmx('sigma',sigma,size(sigma),size(sigma),1)
            call xc_f90_gga_exc_vxc(xc_func,n,r(1,1),sigma(1,1),pe(1),vrho(1,1),vsigma(1))
C           call prmx('vsigma',vsigma,size(vsigma),size(vsigma),1)
            pv(:) = 2*vrho(1,:); pv2(:) = 2*vrho(2,:); pe(:) = 2*pe(:) ! => Atomic Rydberg units
            call daxpy(3*n,2.d0,vsigma,1,dxdg2,1)!JJ, copy out dxdg; alpha=2 for hartree->ryd
            deallocate(r,vrho,sigma,vsigma)
          endif
          call xc_f90_func_end(xc_func) ! destroy functional
C
CC   ... Meta-GGA functionals
        case(XC_FAMILY_MGGA)
          call rx('Meta GGAs not ready in this code yet, sorry')
C          call xc_f90_func_init(xc_func,info,id,nspin)
CC         set Tran-Blaha '09 constant if required
C          if (id == XC_MGGA_X_TB09) then
C            if (present(c_tb09)) call xc_f90_mgga_x_tb09_set_par(xc_func,
C     .        c_tb09)
C          end if
C          if (k == 1) then      ! exchange
C            if (nsp == 1) then
C              allocate(vsigma(1,n),vlapl(1,n),vtau(1,n))
C              call xc_f90_mgga_vxc(xc_func,n,rho(1),grrsqr(1),g2rho(1),tau(1),
C     .          vx(1),vsigma(1,1),vlapl(1,1),vtau(1,1))
C              deallocate(vsigma,vlapl,vtau)
C            else
C              allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
C              allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
C              r(1,:)=rhoup(:); r(2,:)=rhodn(:)
C              sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=
C     .          grrsq2(:)
C              lapl(1,:)=g2up(:); lapl(2,:)=g2dn(:)
C              t(1,:)=tauup(:); t(2,:)=taudn(:)
C              call xc_f90_mgga_vxc(xc_func,n,r(1,1),sigma(1,1),lapl(1,1),t(1,
C     .          1),vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
C              vxup(:)=vrho(1,:); vxdn(:)=vrho(2,:)
C              deallocate(r,sigma,lapl,t)
C              deallocate(vrho,vsigma,vlapl,vtau)
C            end if
C          else                  ! correlation
C            if (nsp == 1) then
C              allocate(vsigma(1,n),vlapl(1,n),vtau(1,n))
C              call xc_f90_mgga_vxc(xc_func,n,rho(1),grrsqr(1),g2rho(1),tau(1),
C     .          vc(1),vsigma(1,1),vlapl(1,1),vtau(1,1))
C              deallocate(vsigma,vlapl,vtau)
C            else
C              allocate(r(2,n),sigma(3,n),lapl(2,n),t(2,n))
C              allocate(vrho(2,n),vsigma(3,n),vlapl(2,n),vtau(2,n))
C              r(1,:)=rhoup(:); r(2,:)=rhodn(:)
C              sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=
C     .          grrsq2(:)
C              lapl(1,:)=g2up(:); lapl(2,:)=g2dn(:)
C              t(1,:)=tauup(:); t(2,:)=taudn(:)
C              call xc_f90_mgga_vxc(xc_func,n,r(1,1),sigma(1,1),lapl(1,1),t(1,
C     .          1),vrho(1,1),vsigma(1,1),vlapl(1,1),vtau(1,1))
C              vcup(:)=vrho(1,:); vcdn(:)=vrho(2,:)
C              deallocate(r,sigma,lapl,t)
C              deallocate(vrho,vsigma,vlapl,vtau)
C            end if
C          end if
C          call xc_f90_func_end(xc_func) ! destroy functional

        case default
          call rxi('xc_libxc: unsupported functional',xctype)
        end select

      enddo

      if (lxcf < 0) then
        call dmadd(ex,n,1,1d0,ec,n,1,1d0,exc,n,1,n,1)
        call dmadd(vx,n,1,1d0,vc,n,1,1d0,vxc,n,1,n,1)
        if (nsp == 2) call dmadd(vx2,n,1,1d0,vc2,n,1,1d0,vxc2,n,1,n,1)
      endif
C#endif
      end

      subroutine xc_libxc_type(lxcf,nsp,xctype)
C- Returns information about the libxc functional
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcf  :compound of two digits, each referring to a libxc functional
Ci         :lxcf/2**16 = index to the first functional
Ci         :mod(lxcf,2**16) = index to the 2nd functional
Ci   nsp   :0   => return in xctype the family (LDA, GGA, Meta-GGA)
Ci         :1,2 => return in xctype the kind (exchange, correlation, or XC)
Ci         :       nsp should be the number of spins
Co Outputs
Co   xctype:if nsp = 0, returns one value from XC_FAMILY, i.e. one of:
Co         :   XC_FAMILY_UNKNOWN = -1
Co         :   XC_FAMILY_NONE = 0
Co         :   XC_FAMILY_LDA = 1
Co         :   XC_FAMILY_GGA = 2
Co         :   XC_FAMILY_MGGA = 4
Co         :   XC_FAMILY_LCA = 8
Co         :   XC_FAMILY_OEP = 16
Co         :   XC_FAMILY_HYB_GGA = 32
Co         :if nsp = 1 or 2, returns kind of functional, i.e. one of:
Co         :   XC_EXCHANGE = 0
Co         :   XC_CORRELATION = 1
Co         :   XC_EXCHANGE_CORRELATION = 2
Co         :   XC_KINETIC = 2
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Nov 13
C ----------------------------------------------------------------------
C#ifdef LIBXC
      use xc_f90_lib_m
C#endif
      implicit none
C ... Passed parameters
      integer lxcf,nsp,xctype(2)

C#ifndefC LIBXC
C      call rx('this density functional requires '//
C     .  'that you link the libxc library; sorry')
C#else
C ... Local arrays and pointers
      type(xc_f90_pointer_t) xc_func,info
C ... Local parameters
      integer k,id,nspin

      do  k = 1, 2
        if (k == 1) id = lxcf/2**16
        if (k == 2) id = mod(lxcf,2**16)

        if (id == 0) then    ! no functional
          xctype(k) = -1
        elseif (nsp == 0) then ! whether LDA, GGA, or Meta-LDA
          xctype(k) = xc_f90_family_from_id(id)
        else
          nspin=XC_UNPOLARIZED; if (nsp > 1) nspin=XC_POLARIZED
          call xc_f90_func_init(xc_func,info,id,nspin)
          xctype(k) = xc_f90_info_kind(info)
          call xc_f90_func_end(xc_func)
        endif
      enddo
C#endif

      end

      subroutine xc_libxc_name(lxcf,nsp,k,name)
C- Returns name of the libxc functional
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcf  :compound of two digits, each referring to a libxc functional
Ci         :lxcf/2**16 = index to the first functional
Ci         :mod(lxcf,2**16) = index to the 2nd functional
Ci   nsp   :0   => return in xctype the family (LDA, GGA, Meta-GGA)
Ci         :1,2 => return in xctype the kind (exchange, correlation, or XC)
Ci         :       nsp should be the number of spins
Co Outputs
Co   xctype:if nsp = 0, returns one value from XC_FAMILY, i.e. one of:
Co         :   XC_FAMILY_UNKNOWN = -1
Co         :   XC_FAMILY_NONE = 0
Co         :   XC_FAMILY_LDA = 1
Co         :   XC_FAMILY_GGA = 2
Co         :   XC_FAMILY_MGGA = 4
Co         :   XC_FAMILY_LCA = 8
Co         :   XC_FAMILY_OEP = 16
Co         :   XC_FAMILY_HYB_GGA = 32
Co         :if nsp = 1 or 2, returns kind of functional, i.e. one of:
Co         :   XC_EXCHANGE = 0
Co         :   XC_CORRELATION = 1
Co         :   XC_EXCHANGE_CORRELATION = 2
Co         :   XC_KINETIC = 2
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Nov 13
C ----------------------------------------------------------------------
C#ifdef LIBXC
      use xc_f90_lib_m
C#endif
      implicit none
C ... Passed parameters
      integer lxcf,nsp,k
      character*(*) name
C#ifndefC LIBXC
C      call rx('this density functional requires '//
C     .  'that you link the libxc library; sorry')
C#else
C ... Local arrays and pointers
      type(xc_f90_pointer_t) xc_func,info
C ... Local parameters
      character*(120) strn1, strn3, strn2*(3)
      integer id,nspin

      if (k == 1) id = lxcf/2**16
      if (k == 2) id = mod(lxcf,2**16)

      name = '*'
      if (id == 0) return

      select case(xc_f90_family_from_id(id))
      case (XC_FAMILY_LDA);      write(strn1,'(a)') "LDA"
      case (XC_FAMILY_GGA);      write(strn1,'(a)') "GGA"
      case (XC_FAMILY_HYB_GGA);  write(strn1,'(a)') "Hybrid GGA"
      case (XC_FAMILY_MGGA);     write(strn1,'(a)') "MGGA"
      case (XC_FAMILY_LCA);      write(strn1,'(a)') "LCA"
      end select

      nspin=XC_UNPOLARIZED; if (nsp > 1) nspin=XC_POLARIZED
      call xc_f90_func_init(xc_func,info,id,nspin)

      select case(xc_f90_info_kind(info))
      case(XC_EXCHANGE)
        write(strn2, '(a)') ' X'
      case(XC_CORRELATION)
        write(strn2, '(a)') ' C'
      case(XC_EXCHANGE_CORRELATION)
        write(strn2, '(a)') ' XC'
      end select
      call xc_f90_info_name(info, strn3)

      name = trim(strn1) // trim(strn2) // ' (' // trim(strn3) // ')'
      call xc_f90_func_end(xc_func)
C#endif

      end
