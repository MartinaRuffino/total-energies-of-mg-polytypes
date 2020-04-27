      subroutine tbmpol(nbas,nsp,nl,nlmq,qnu,ipc,lmxl,gaunt,qpol,
     .                  rho,rhoc,qmpol,mmom)
C- Make multipole (and magnetic) moments from TB "density matrix"
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,nl,nsp,qnu,ipc,gaunt
Ci   nlmq : L-cutoff for multipoles, also leading dimension of qmpol
Ci   lmxl : species-resolved l-cutoffs for multipoles
Ci   qpol : polarisation parameters (delta_ll'l'')
Ci   rho, rhoc (see tbfrce); if mixrho then these are from rhomix.
Co Outputs:
Co   qmpol: electrostatic point multipoles Q_L [Finnis, Eq.(7.76)]
Co   mmom : magnetic moments
Cr Remarks
Cr   Our multipoles are defined such that the monopole moment is the
Cr   total Mulliken charge, while the higher multipoles are approximated
Cr   as an onsite sum without reference to overlap (Finnis, p. 216)
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in) :: nbas,nl,nlmq,nsp,ipc(nbas),lmxl(*)
      double precision, intent(in)  :: qpol(10,*),qnu(3,0:nl-1,nsp,*),
     .                                 gaunt(9,9,25),rho(nl,nsp,nbas),
     .                                 rhoc(nl**2,nl**2,nbas)
      double precision, intent(out) :: qmpol(nlmq,nbas),mmom(nbas)
C Local Variables
      integer ib,ic,ilm,ilmp,ilmpp,iprint,nlm
      integer ll
      double precision q0,M,dsum
c...deb
      double precision qmpt(nlmq)
c...deb
      real(8) :: qdown

      call tcn('tbmpol')

c... Checks
      if (nlmq > 25)
     .  call rxi(' tbmpol: nlmq is too big, nlmq = ',nlmq)

c     call dcopy(9*nbas,0d0,0,qmpol,1)
c     call dcopy(nbas,0d0,0,mmom,1)
      qmpol = 0d0
      mmom = 0d0

C --- get multipole moments ---
      do  ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxl(ic)+1)**2

        q0 = dsum(nl,qnu(1,0,1,ic),3)
        qmpol(1,ib) = dsum(nl,rho(1,1,ib),1)
        if (nsp /= 1) then
          q0 = q0 + dsum(nl,qnu(1,0,2,ic),3)
          qdown = dsum(nl,rho(1,2,ib),1)
          mmom(ib) = qmpol(1,ib) - qdown
          qmpol(1,ib) = qmpol(1,ib) + qdown
        end if
        qmpol(1,ib) = qmpol(1,ib) - q0

!         if (nsp == 1) then
!           q0 = dsum(nl,qnu(1,0,1,ic),3)
!           qmpol(1,ib) = dsum(nl,rho(1,1,ib),1) - q0
!         else
!           q0 = dsum(nl,qnu(1,0,1,ic),3) + dsum(nl,qnu(1,0,2,ic),3)
!           qmpol(1,ib) = dsum(nl,rho(1,1,ib),1)
!      .                + dsum(nl,rho(1,2,ib),1) - q0
!         endif

        if (iprint() > 60) write (*,100)
        if (nlm > 1) then
          do  ilm = 2, nlm
            do  ilmp  = 1, nl**2
              do  ilmpp = 1, nl**2
                call getM(ilm,ilmp,ilmpp,qpol(1,ic),M)
                qmpol(ilm,ib) = qmpol(ilm,ib) +
     .            rhoc(ilmp,ilmpp,ib) * M * gaunt(ilmp,ilmpp,ilm)

C --- verbose output ---
                if (iprint() > 60 .and. M /= 0d0
     .              .and. gaunt(ilmp,ilmpp,ilm) /= 0d0) then
                  write (*,200)
     .                ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .                M,gaunt(ilmp,ilmpp,ilm),rhoc(ilmp,ilmpp,ib)
                endif
C ---------------------

              enddo
            enddo
          enddo
        endif
C        if (iprint() > 30) then
C          call awrit3('        Q^e/e=%d, %3:1d, %5:1d',' ',180,
C     .                i1mach(2),qmpol(1,ib),qmpol(2,ib),qmpol(5,ib))
C        endif
! Functionality replaced above
! C --- make magnetic moments ---
!         if (nsp == 2) then
!           mmom(ib) = dsum(nl,rho(1,1,ib),1) - dsum(nl,rho(1,2,ib),1)
!         endif
      enddo

  100 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .      rho_L''L''''')
  200 format (6(i3,2x),2(2x,f6.2,2x),2f10.6)

c...deb
C total MP moments of the cell
      if (iprint() >= 40) then
        qmpt = 0d0
        do ib = 1, nbas
          nlm = (lmxl(ipc(ib))+1)**2
          qmpt(1:nlm) = qmpt(1:nlm) + qmpol(1:nlm,ib)
        enddo
        print *,' tbmpol: total multipole moments of the cell = '
        print '(1x,25f8.4)',qmpt(1)
        if (nlmq >= 4) print '(1x,25f8.4)',qmpt(2:4)
        if (nlmq >= 9) print '(1x,25f8.4)',qmpt(5:9)
        if (nlmq >= 16) print '(1x,25f8.4)',qmpt(10:16)
        if (nlmq >= 25) print '(1x,25f8.4)',qmpt(17:25)
      endif
c...deb

      call tcx('tbmpol')

      end subroutine tbmpol



!       write(520,"('')")
!       write(520,'("rho",/,12(x,f16.8),/)') rho
!       write(520,'("rhoc",/,2(9(9(x,f16.8),/),/))') rhoc
!       write(520,'("qmpol",/,2(25(x,f16.8),/))') qmpol
!       write(520,'("qnu",/,2(2(9(x,f16.8),/),/))') qnu(:,:,:,:2)
!       write(520,'("qpol",/,2(10(x,f16.8),/))') qpol(:,:2)
