      subroutine rstr0(nxi,lxi,exi,nlmx,np,x,y,z,lmxa,iop,hl,hd)
C- Reduced strux for a vector of energies, standard definition IV-43.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nxi   :number of energies for each connecting vector
Ci   lxi   :lmax for each energy and point
Ci   exi   :list of nxi energies
Ci   nlmx  :leading dimension of hl,hd
Ci   np    :number of connecting vectors
Ci   x,y,z :connecting vectors
Ci   lmxa  :Make hl for lxi + lmxa
Ci   iop   :1: calculate h only
Ci          any other number: calculate both h and hdot
Co Outputs
Co   hl    :reduced strux for each energy, point, to lxi + lmxa
Co   hd    :reduced energy derivative
Cr Remarks
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer np,nxi,lxi(nxi,np),lmxa,nlmx,iop
      double precision hl(nlmx,nxi,np),hd(nlmx,nxi,np),
     .  exi(nxi),x(np),y(np),z(np)
C ... Local parameters
      integer ie,ip,lmax,ilm,l,m,lmxx,lmxy
      parameter (lmxx=20)
      double precision psi(-1:lmxx),phi(-1:lmxx),yl((lmxx+1)**2),r2,
     .  rfac,xxx,fac2l(0:lmxx),psid,psil
      real(8), pointer :: ylm(:)
      real(8), pointer :: r2l(:)

C --- Make ylm for all points, up to lmxy, which find ---
      lmxy = -1
      do  ip = 1, np
        do  ie = 1, nxi
          lmxy = max(lmxy,lxi(ie,ip))
        enddo
      enddo
      lmxy = lmxy+lmxa
      fac2l(0) = 1
      do  l = 1, lmxy
        fac2l(l) = fac2l(l-1)*(2*l-1)
      enddo
      if (lmxy > lmxx) call rxi('rstr0: lmax exceeds lmxx:',lmxy)
      allocate(ylm((lmxy+1)**2*np),r2l(np))
      call ropyln(np,x,y,z,lmxy,np,ylm,r2l)

C --- For each point, do ---
      do  ip = 1, np
      call pvstr0(lmxy,ip,np,ylm,r2l,yl,r2)
      if (r2 < 1d-10) then
        call dpzero(hl(1,1,ip),nlmx*nxi)
        if (iop /= 1) call dpzero(hd(1,1,ip),nlmx*nxi)
          cycle
      endif

C --- Reduced strx hl, or hl and hd for all energies, this point ---
      if (iop == 1) then
        do  ie = 1, nxi
          lmax = lmxa+lxi(ie,ip)
C     ... Not worth vectorizing (I think)
c         call bessl2(exi(ie)*r2,0,lmax,phi(0),psi(0))
          call besslr(exi(ie)*r2,0,-1,lmax,phi,psi)
          ilm = 0
          rfac = dsqrt(r2)
          xxx = 1d0/r2
          do  l = 0, lmax
            rfac = rfac*xxx
C       ... undo fac2l scaling, to recover standard MSM IV-43.
C           psil = rfac * psi(l) * fac2l(l)
            psil = rfac * psi(l)
            do  m = -l, l
              ilm = ilm+1
              hl(ilm,ie,ip) = psil*yl(ilm)
            enddo
          enddo
        enddo
      else
        do  ie = 1, nxi
          lmax = lmxa+lxi(ie,ip)
C         call bessl2(exi(ie)*r2,-1,lmax,phi(-1),psi(-1))
          call besslr(exi(ie)*r2,0,-1,lmax,phi(-1),psi(-1))
          ilm = 0
          rfac = dsqrt(r2)
          xxx = 1d0/r2
          do  l = 0, lmax
            rfac = rfac*xxx
C           psil = rfac * psi(l) * fac2l(l)
C           psid = rfac * psi(l-1)*r2/(4*l-2) * fac2l(l)
            psil = rfac * psi(l)
            psid = rfac * psi(l-1)*r2/2
            do  m = -l, l
              ilm = ilm+1
              hl(ilm,ie,ip) = psil*yl(ilm)
              hd(ilm,ie,ip) = psid*yl(ilm)
            enddo
          enddo
        enddo
      endif
      enddo
      deallocate(ylm,r2l)
      end
      subroutine pvstr0(lmax,ip,nd,ylm,r2,yl,rsq)
      implicit none
      integer lmax,ip,nd,nlm,ilm
      double precision ylm(nd,1), yl(*), r2(ip), rsq

      nlm = (lmax+1)**2
      rsq = r2(ip)
      do  ilm = 1, nlm
        yl(ilm) = ylm(ip,ilm)
      enddo
      end
