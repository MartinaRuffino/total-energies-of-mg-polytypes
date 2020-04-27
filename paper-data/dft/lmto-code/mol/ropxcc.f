      subroutine ropxcc(npx,rhos,rhot,grh,g2rh,ggrh,
     .  agrh,grgagr,grpgrm,exc,vxc,sum1,sum2,np)
C- Smooth density's vxc and exc for a mesh of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   np     : number of points
Ci   npx    : leading dimension of the density arrays
Ci   rhot   : spin-resolved density
Ci   grh    :               grad rho
Ci   g2rh   :               matrix of second derivatives of rho
Ci   ggrh   :               laplacian rho
Co Outputs
Co   rhos   : total density: rhos = rho_up + rho_down
Co   agrh   : spin-resolved |grad rho|
Co   grgagr : grad rho . grad abs grad rho; 1 - up, 2 - down, 3 - tot
Co   grpgrm : grad rho_up . grad rho_down
Cu Updates
Cu   21 Jul 07 (S. Lozovoi) Generates gradients for GGA
C ----------------------------------------------------------------------
C
      implicit none
      integer npx,np
      double precision rhos(npx),rhot(npx,*),exc(npx),vxc(npx,*)
      double precision grh(npx,3,*),g2rh(npx,5,*),ggrh(npx,*)
      double precision agrh(npx,*),grgagr(npx,*),grpgrm(npx)
      double precision sum1,sum2
c Local variables
      integer nsp,lsp,lxcf,lxcg,lcut
      integer oexcx,oexcc,ovxcx,ovxcc,ovxcx2,ovxcc2
      integer i,ip,ineg(2)
      integer i1mach,iprint
      double precision rhomin,tol
c     double precision repnl(2),rmunl(2)
      double precision rhx,rhy,rhz,rhxx,rhxy,rhxz,rhyy,rhyz,rhzz
c for debugging:
C     integer iii,j
c
      data rhomin/1d-15/,tol/1.d-15/
      real w(1)
      common /w/ w

      call tcn('ropxcc')
      if (npx < np) call rxi('ropxcc: increase npx, needed',np)
c              print *, 'np,npx =', np,npx

      nsp=lsp()+1
C --- If negative density, set it to rhomin and produce a warning ---
      call tcn('density cleaning')
        do  i = 1, nsp
          ineg(i) = 0
          do  ip = 1, np
              if(rhot(ip,i) < 0d0) then
                ineg(i) = ineg(i) + 1
                rhot(ip,i) = rhomin
              endif
          enddo
        enddo
        if(iprint() >= 40) then
          if(nsp. eq. 1) then
            if(ineg(1) /= 0)
     .      call awrit1(' ropxcc (warning):  %i non positive density'//
     .      ' points',' ',128,i1mach(2),ineg(1))
          else
            if(ineg(1) /= 0)
     .      call awrit1(' ropxcc (warning):  %i non positive density'//
     .      ' points for the majority spin',' ',128,i1mach(2),ineg(1))
            if(ineg(2) /= 0)
     .      call awrit1(' ropxcc (warning):  %i non positive density'//
     .      ' points for the minority spin',' ',128,i1mach(2),ineg(2))
          endif
        endif
      call tcx('density cleaning')

      call tcn('make local vxc')
c ------- rhos = rhop+ + rhop-  ---------------------
      call dpzero(rhos,np)
      do  i = 1, nsp
        call daxpy(np,1d0,rhot(1,i),1,rhos,1)
      enddo

C --- Make local eps and mu ---
      call defrr(oexcx,np)
      call defrr(oexcc,np)
      call defrr(ovxcx,np)
      call defrr(ovxcc,np)
      call defrr(ovxcx2,np)
      call defrr(ovxcc2,np)
      if (lxcf() > 2) then
        call evxcp(rhot,rhot(1,nsp),np,nsp,lxcf(),w(oexcx),w(oexcc),exc,
     .    w(ovxcx),w(ovxcx2),w(ovxcc),w(ovxcc2),vxc,vxc(1,nsp))
      else
        do  i = 1, nsp
          call evxcv(rhos,rhot(1,i),np,nsp,lxcf(),
     .      exc,w(oexcx),w(oexcc),
     .      vxc(1,i),w(ovxcx),w(ovxcc))
        enddo
      endif
      call rlse(oexcx)
      call tcx('make local vxc')
c ---------- if GGA, make non-local xc ----------------
      if(lxcg(). ne. 0) then
      call tcn('make non-local vxc')
c --------- prepare |grad rho|, (grad rho) * grad |grad rho|, etc.
C --- agrh : |grad rho(r)| ---
      do  i = 1, nsp
        do  ip = 1, np
          agrh(ip,i) = dsqrt(grh(ip,1,i)**2
     .                     + grh(ip,2,i)**2
     .                     + grh(ip,3,i)**2)
        enddo
      enddo
C     call prm3('abs grad rho',agrq,n1+2,n1,n2,n3)
C --- agrh(3) : |grad tot rho|;  grpgrm : grad rho+ . grad rho- ---
      if (nsp == 2) then
        do  ip = 1, np
          agrh(ip,3) =
     .      dsqrt((grh(ip,1,1)+grh(ip,1,2))**2 +
     .            (grh(ip,2,1)+grh(ip,2,2))**2 +
     .            (grh(ip,3,1)+grh(ip,3,2))**2)
        enddo
c this term is only needed for L-M GGA
        if(lxcg() == 1) then
          do  ip = 1, np
            grpgrm(ip) =
     .             grh(ip,1,1)*grh(ip,1,2) +
     .             grh(ip,2,1)*grh(ip,2,2) +
     .             grh(ip,3,1)*grh(ip,3,2)
          enddo
        endif
      endif
C --- grgagr(1-2) : grad rho . grad |grad rho| ---
      do  i = 1, nsp
        do  ip = 1, np
c set to zero if in no-man's land or if |grad rho| = 0
          if(rhot(ip,i) < rhomin .or. agrh(ip,i) < tol) then
            grgagr(ip,i) = 0.d0
          else
            rhx=grh(ip,1,i)
            rhy=grh(ip,2,i)
            rhz=grh(ip,3,i)
            rhxx=g2rh(ip,1,i)
            rhxy=g2rh(ip,2,i)
            rhxz=g2rh(ip,3,i)
            rhyz=g2rh(ip,4,i)
            rhzz=g2rh(ip,5,i)
            rhyy=ggrh(ip,i) - rhxx - rhzz
            grgagr(ip,i) = (rhx*rhx*rhxx + rhy*rhy*rhyy + rhz*rhz*rhzz
     .      + 2.d0 * (rhx*rhy*rhxy + rhx*rhz*rhxz + rhx*rhz*rhxz))/
     .      agrh(ip,i)
          endif
        enddo
      enddo
C --- grgagr(3) : (grad tot rho) . grad |grad tot rho|  ---
      if (nsp == 2) then
        do  ip = 1, np
          if(rhos(ip) < rhomin) then
            grgagr(ip,3) = 0.d0
          else
            rhx=grh(ip,1,1) + grh(ip,1,2)
            rhy=grh(ip,2,1) + grh(ip,2,2)
            rhz=grh(ip,3,1) + grh(ip,3,2)
            rhxx=g2rh(ip,1,1) + g2rh(ip,1,2)
            rhxy=g2rh(ip,2,1) + g2rh(ip,2,2)
            rhxz=g2rh(ip,3,1) + g2rh(ip,3,2)
            rhyz=g2rh(ip,4,1) + g2rh(ip,4,2)
            rhzz=g2rh(ip,5,1) + g2rh(ip,5,2)
            rhyy=ggrh(ip,1) + ggrh(ip,2) - rhxx - rhzz
            grgagr(ip,3) = (rhx*rhx*rhxx + rhy*rhy*rhyy + rhz*rhz*rhzz
     .      + 2.d0 * (rhx*rhy*rhxy + rhx*rhz*rhxz + rhx*rhz*rhxz))/
     .      agrh(ip,3)
          endif
        enddo
      endif
c******** for debugging: selectively set terms to zero  *******
c               call dpzero(agrh,npx*nsp)
c               print *,' **** agrh is set to zero *****'
c               call dpzero(ggrh,npx*nsp)
c               print *,' **** ggrh is set to zero *****'
c               call dpzero(agrh,npx*(2*nsp-1))
c               print *,' **** agrh(i,2*nsp-1) is set to zero *****'
c               call dpzero(grpgrm,npx)
c               print *,' **** grpgrm is set to zero *****'
c               call dpzero(grgagr,npx*nsp)
c               print *,' **** grgagr is set to zero *****'
c               call dpzero(grgagr,npx*(2*nsp-1))
c               print *,' **** grgagr(i,2*nsp-1) is set to zero *****'
c*****************************************************
C --- Nonlocal exc and mu for all points  ---
c     call dpzero(vxc,np*nsp)
c     call dpzero(exc,np)
      if (lxcg() > 2) then
        call vxcgga(lxcg(),np,nsp,rhot(1,1),rhot(1,nsp),agrh(1,1),
     .    agrh(1,nsp),ggrh(1,1),ggrh(1,nsp),
     .    agrh(1,2*nsp-1),grpgrm,grgagr(1,2*nsp-1),
     .    grgagr(1,1),grgagr(1,nsp),vxc(1,1),
     .    vxc(1,nsp),exc)
      else
        lcut = 1
        if (lcut == 1) then
          call vxnlcc(np,nsp,rhot(1,1),rhot(1,nsp),agrh(1,1),
     .      agrh(1,nsp),ggrh(1,1),ggrh(1,nsp),
     .      agrh(1,2*nsp-1),grpgrm,grgagr(1,2*nsp-1),
     .      grgagr(1,1),grgagr(1,nsp),vxc(1,1),
     .      vxc(1,nsp),exc)
        else
          call vxnloc(np,nsp,rhot(1,1),rhot(1,nsp),agrh(1,1),
     .      agrh(1,nsp),ggrh(1,1),ggrh(1,nsp),
     .      agrh(1,2*nsp-1),grpgrm,grgagr(1,2*nsp-1),
     .      grgagr(1,1),grgagr(1,nsp),vxc(1,1),
     .      vxc(1,nsp),exc)
        endif
      endif
      call tcx('make non-local vxc')

C --- Make nlocal reps, rmu ---
c     do  i = 1, nsp
c       call dpdot(rhot(1,i),exc,np,repnl(i))
c       call dpdot(rhot(1,i),vxc(1,i),np,rmunl(i))
c       repnl(i) = repnl(i)*vol0
c       rmunl(i) = rmunl(i)*vol0
c     enddo

C ------- end of GGA if-loop
      endif

c ------- add to total charge and xc energy -------
        do i=1,np
          sum1 = sum1 + rhos(i)
          sum2 = sum2 + rhos(i)*exc(i)
        enddo
c ------- overwrite rhot with vxc -----------------
      call dcopy(npx*nsp,vxc,1,rhot,1)
c     do 4 isp=1,nsp
c       do 3 i=1,np
c 3     rhot(i,isp)=vxc(i,isp)
c 4   continue

      call tcx('ropxcc')
      end
