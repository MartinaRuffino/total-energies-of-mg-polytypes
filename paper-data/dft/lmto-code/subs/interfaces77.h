      interface
        subroutine hpairm(nbas,ips,alat,plat,pos,rham,nttab,ntab,iax,
     .                    rtab,mxcsiz,ordrd,niax)
        implicit none
        integer :: nbas,nttab,ips(nbas),ntab(nbas+1),mxcsiz,ordrd,niax
        real(8) :: alat,plat(3,3),pos(3,nbas),rham(nbas),rtab(3,*)
!#ifndef LINUXF
        integer, allocatable, intent(inout) :: iax(:)
!#else!
!        integer, pointer :: iax(:)
!#endif
        end subroutine hpairm

        subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,
     .                   iltab,pltab,nttab,ntab,iax,mxcsiz)
        implicit none
        integer :: nbas,nbasp,nttab,nd,iltab,mxcsiz
        integer :: ntab(nbasp+1),ipsp(*),pltab(*)
        real(8) :: alat,plat(9),rmax(*),baspp(3,*)
!#ifndef LINUXF
        integer, intent(inout), allocatable :: iax(:)
!#else!
!        integer, pointer :: iax(:)
!#endif
        end subroutine pairs

        subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,lrots,eula,acof,rofi,
     .  rwgt,nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,
     .  amom,rho1,rho2)
        implicit none
        logical lrots
        integer nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca
!       Noncollinear case: rho1,rho2,rhoc have 4 components
        double precision qcsm,qloc,rsmv,rsmfa,amom(0:3),eula(3),
     .    rho1(nr,nlml,nsp*2),rofi(nr),rwgt(nr),exi(1),hfc(nxi0,nsp),
     .    rho2(nr,nlml,nsp*2),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp*2)
        real(8),target :: rhofa(nr,nsp)
        end subroutine p2ovlc

        subroutine i1alloc(p_iarr,n,linit)
        implicit none
        integer n,linit
        integer, pointer :: p_iarr(:)
        end subroutine i1alloc

        subroutine d1alloc(p_darr,n,linit)
        integer n,linit
        real(8), pointer :: p_darr(:)
        end subroutine d1alloc

        subroutine unievp(ovlp, herm, c2d, itype, jobz, rng, uplo, n, a,
     &    ia, ja, desca, b, ib, jb, descb, vl, vu, il, iu, m, nz, w, z,
     &    iz, jz, descz)
        use mod_ctx, only : cart_ctx_t
        implicit none
        logical, intent(in) :: ovlp, herm
        character, intent(in) :: jobz, rng, uplo
        integer, intent(in) :: ia, ib, itype, il, iu, iz, ja, jb, jz, n,
     &                                      desca(9), descb(9), descz(9)
        integer, intent(inout) :: m, nz
        real(8), intent(in) :: vl, vu
        real(8), intent(inout) :: w(*)
        real(8), intent(in), target :: a(*), b(*)
        real(8), intent(inout), target :: z(*)
        type(cart_ctx_t), intent(in) :: c2d
        end subroutine unievp
      end interface


