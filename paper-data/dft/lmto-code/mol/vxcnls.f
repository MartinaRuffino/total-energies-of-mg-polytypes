      subroutine vxcnls(ri,lcut,nr,np,nlm,nsp,rp,grp,ggrp,agrp,
     .  yl,gyl,ylwp,rwgt,wp,rl,agrl,lxcg,vxcnl,excnl,vl,rep,rmu)
C- Gradient correction to nspher. density on a radial and angular mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   ri    :mesh of points
Ci   lcut  :1 if cutoff exc for small rho to avoid blowup in exc
Ci   nr    :number of radial mesh points
Ci   np    :number of points for angular integration
Ci   nlm   :maximum (l+1)**2
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   yl    :Ylm's tabulated on angular mesh
Ci   gyl   :Gradient of Ylm's tabulated on angular mesh
Ci   ylwp  :Ylm*wp, where wp are weights for angular mesh
Ci   rwgt  :radial mesh weights
Ci   wp    :angular mesh weights (not needed unless debugging)
Ci   rl    :density on radial mesh
Ci   agrl  :|rl|
Ci   lxcg  :indicates which nonlocal XC functional to use
Ci         :  0    No nonlocal functional
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Co Outputs
Co   rp    :spin pol density on the combined radial and angular mesh
Co   grp   :grad rp
Co   ggrp  :Laplacian rp
Co   agrp  :|grad rp|
Co   vxcnl  :nonlocal XC potential on the combined radial * angular mesh
Co   excnl  :nonlocal XC energy on the combined radial * angular mesh
Co   vl    :vxcnl is added to vl
Co   rep   :int rho * excnl added to rep
Co   rmu   :int rho * vxcnl added to rmu
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Apr 09 reduced the calling arguments for better portability
Cu   29 Apr 05 (ATP) adaped to new vxcnsp
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lcut,nr,np,nlm,nsp,lxcg
      double precision ri(3),yl(np,nlm),gyl(np,nlm,3),ylwp(np,nlm),
     .  wp(np),rep(2),rmu(2),rwgt(nr),vxcnl(nr,nlm,nsp),excnl(nr,nlm),
     .  rl(nr,nlm,nsp),agrl(nr,nlm,nsp),agrp(nr,np,nsp),
     .  vl(nr,nlm,nsp),rp(nr,np,nsp),grp(nr,np,3,nsp),ggrp(nr,np,nsp)
C ... Local parameters
      double precision pi,vhold,tol,sumnl(0:20),repnl(4),rmunl(4),weight
      double precision wk(nr,nsp),wk2(nr,nsp)
      integer ilm,ip,ipr,ir,i,l,ll,lmax,nn,oagrp,ogagrp,lopgfl,
     .  ovxcp,oexcp,owk,owkl,lmx,jx
      integer ir0,ir1,nri,incn
      logical lz
      parameter (incn=50)
      real w(1)
      common /w/ w
      data tol/1d-15/

C     call pshpr(80)
      call tcn('vxcnls')
      if (lxcg == 0) return
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      nn = 6
      call dpzero(repnl,4)
      call dpzero(rmunl,4)
      call dpzero(excnl,nr*nlm)
      call dpzero(vxcnl,nr*nlm*nsp)
      lz = ri(1) == 0

C --- Make ylwp = yl*wp for fast multiplication ---
C      do  ilm = 1, nlm
C        do  ip = 1, np
C          ylwp(ip,ilm) = yl(ip,ilm)*wp(ip)
C        enddo
C      enddo

C --- Generate density point-wise through sphere ---
      do  i = 1, nsp
        call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0,
     .             rp(1,1,i),nr)
      enddo

C --- If negative density, set to tol ---
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              if (rp(ir,ip,i) <= 0d0) rp(ir,ip,i) = tol
            enddo
          enddo
        enddo

C --- Potential from spherical part of density (valid for small r) ---
      if (lz) then
        do  i = 1, nsp
          call dpcopy(rl(1,1,i),wk(1,i),1,nr,dsqrt(4*pi))
          do  ir = 1, nr
          wk(ir,i) = wk(ir,i)*ri(ir)**2
          enddo
        enddo
        call dpzero(wk2,nr*nsp)
        call vxc0gc(nr,nsp,ri,rwgt,wk,wk2,excnl,repnl(3),rmunl(3),
     .    100*lxcg)
        call dscal(nr,dsqrt(4*pi),excnl,1)
        do  i = 1, nsp
          if (ipr >= 30 .and. i == 1)
     .      print 4, rmunl(i+2),repnl(i+2),'  (l=0 rho)'
          if (ipr >= 30 .and. i == 2)
     .      print 5, rmunl(i+2),repnl(i+2),
     .      rmunl(3)+rmunl(4),repnl(3)+repnl(4)
          call dpcopy(wk2(1,i),vxcnl(1,1,i),1,nr,dsqrt(4*pi))
        enddo
C       call prrmsh('l=0 vxcnl in vxcnls',ri,vxcnl,nr,nr,1)
C       call prrmsh('rl in vxcnls',ri,rl,nr,nr,1)
C       call prrmsh('wl in vxcnls',ri,rwgt,nr,nr,1)
      endif

C --- Gradient of density point-wise through sphere ---
      lopgfl = 0
      if (lz) lopgfl = 10
      do  i = 1, nsp
      call gradfl(ll(nlm),nlm,nr,np,1,nr,1,lopgfl,nn,ri,yl,gyl,
     .    rl(1,1,i),grp(1,1,1,i),ggrp(1,1,i))
      enddo

C --- agrp, agrl = abs grad rho and its Yl-projection ---
      do  i = 1, nsp
      do  ip = 1, np
      do  ir = 1, nr
        agrp(ir,ip,i) =
     .    dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
      enddo
      enddo
      enddo
      do  i = 1, nsp
        call dgemm('N','N',nr,nlm,np,1d0,agrp(1,1,i),nr,ylwp,np,0d0,
     .    agrl(1,1,i),nr)
      enddo

C --- Do gradient in blocks nr-incn..nr, nr-2*incn..nr-incn ... ---
      ir1  = nr
      lmax = ll(nlm)
   10 continue
        ir0 = max(ir1-incn,1)
        if (lz) ir0 = max(ir1-incn,2)
        nri = ir1-ir0+1
        if (nri > 0) then
        vhold = (vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2

C ...   Gradient-corrected Vxc for points between ir0 and ir1
        call defrr(oagrp,nri*np*(3*nsp-2))
        call defrr(ogagrp,nri*np*nsp*3)
        call defrr(ovxcp,nri*np*nsp)
        call defrr(oexcp,nri*np*nsp)
        call defrr(owk,nr)
        call defrr(owkl,nri*nlm*nsp)
        call xxcnls(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,ylwp,
     .    wp,w(owkl),rp,ggrp,grp,agrl,w(oagrp),w(ogagrp),rl,rwgt,lcut,
     .    w(ovxcp),w(oexcp),vxcnl,excnl,sumnl)
        call rlse(oagrp)
        ir1 = ir0-1

C ... Check rmu to determine largest lmax next pass
        do  l = lmax, 0,-1
          lmx = l
          if (dabs(sumnl(l)) > 1d-7) exit
        enddo
        lmax = lmx

        if (dabs(vhold-(vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2)
     . >= 1d-6 .or. dabs(vhold) <= 1) goto 10
        endif

C --- Add nonlocal vxc into vl ----
      call daxpy(nr*nlm*nsp,1d0,vxcnl,1,vl,1)
      if (lz) then
        do  i = 1, nsp
        vl(1,1,i) = (vl(2,1,i)*ri(3)-vl(3,1,i)*ri(2))/(ri(3)-ri(2))
        jx = 1
        call polint(ri(2),vl(2,1,i),nr-1,nn,ri,0d0,0,jx,vl(1,1,i),vhold)
        if (ipr >= 50 .and. dabs(vhold) > dabs(vl(1,1,i)/100))
     .        print 1,vl(1,1,i),vhold/vl(1,1,i)*100
    1     format(' vxcnls (warning): expect error in V at origin:',
     .  'V=',1pe10.3,' est err=',0pf7.1,'%')
        enddo
      endif
C     call prrmsh('nlocal v(l=0)',ri,vxcnl,nr*nlm,nr,1)

C --- Nonlocal rho*exc, rho*vxc ---
      do  i = 1, nsp
        do  ilm = 1, nlm
        do  ir = 1, nr
          weight = ri(ir)**2*rwgt(ir)
          rmunl(i) = rmunl(i) + rl(ir,ilm,i)*vxcnl(ir,ilm,i)*weight
          repnl(i) = repnl(i) + rl(ir,ilm,i)*excnl(ir,ilm)*weight
        enddo
        enddo
      if (ipr >= 30 .and. i == 1) print 4, rmunl(i),repnl(i)
      if (ipr >= 30 .and. i == 2) print 5, rmunl(i),repnl(i),
     .  rmunl(1)+rmunl(2),repnl(1)+repnl(2)
      rep(i) = rep(i) + repnl(i)
      rmu(i) = rmu(i) + rmunl(i)
      enddo
      call tcx('vxcnls')

C --- Print out rmu by angular momentum ---
      if (ipr < 35) return
      lmax = ll(nlm)
      do  i = 1, nsp
        do  l = 0, lmax
          sumnl(l) = 0d0
        enddo
        do  ilm = 1, nlm
          l = ll(ilm)
          do  ir = 1, nr
            sumnl(l) = sumnl(l) +
     .                 ri(ir)**2*rl(ir,ilm,i)*vxcnl(ir,ilm,i)*rwgt(ir)
          enddo
        enddo
        if (i == 1) print 2, (sumnl(l),l=0,lmax)
        if (i == 2) print 3, (sumnl(l),l=0,lmax)
      enddo
    2 format(' rvnlc by L: ',f12.6,4f10.6:/(15x,4f10.6))
    3 format('     spin 2: ',f12.6,4f10.6:/(15x,4f10.6))
    4 format(' vxcnls: nlc rmu=',f11.6,'  rep=',f11.6,a)
    5 format(' spin 2:         ',f11.6,'      ',f11.6/
     .       '  total:         ',f11.6,'      ',f11.6)


C      if (ipr >= 40) print 887,
C     .  vl(1,1,1), vxcnl(1,1,1), vl(nr,1,1), vxcnl(nr,1,1)
C  887 format(' V_0(0)=',f15.6,'  nloc VXC_0(0)=',f12.6/
C     .       ' V_0(R)=',f15.6,'  nloc VXC_0(R)=',f12.6)

C     call poppr
      end

      subroutine xxcnls(lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,yl,gyl,
     .  ylwp,wp,wkl,rp,ggrp,grp,agrl,agrp,gagrp,rl,rwgt,lcut,vxcp,excp,
     .  vxcnl,excnl,sumnl)
      implicit none
      integer ir0,ir1,nr,np,nlm,nsp,nn,lmax,lcut,lxcg
      double precision rp(nr,np,nsp),ggrp(nr,np,nsp),ri(1),wp(np),
     .  grp(nr,np,3,nsp),agrp(ir0:ir1,np,nsp),vxcnl(nr,nlm,nsp),
     .  excnl(nr,nlm),agrl(nr,nlm,nsp),gagrp(ir0:ir1,np,3,nsp),
     .  yl(np,nlm),ylwp(np,nlm),gyl(np,nlm,3),wkl(ir0:ir1,nlm,nsp),
     .  vxcp(ir0:ir1,np,nsp),excp(ir0:ir1,np),rl(nr,nlm,nsp),rwgt(nr),
     .  sumnl(0:20)
      integer ll,ir,ip,i,nri,ilm,l,iprint,mlm,lopgfl

C     call pshpr(80)
      nri = ir1-ir0+1

C --- gagrp(store in vxcp) = grad rho . grad abs grad rho ---
      lopgfl = 0
      if (ri(ir0) == 0) lopgfl = 10
      do  i = 1, nsp
      call gradfl(lmax,nlm,nr,np,ir0,ir1,0,lopgfl,nn,ri,yl,gyl,
     .    agrl(1,1,i),gagrp(ir0,1,1,i),0d0)
      enddo
      do  i = 1, nsp
        do  ip = 1, np
          do  ir = ir0, ir1
            vxcp(ir,ip,i) = gagrp(ir,ip,1,i)*grp(ir,ip,1,i)
     .                     +gagrp(ir,ip,2,i)*grp(ir,ip,2,i)
     .                     +gagrp(ir,ip,3,i)*grp(ir,ip,3,i)
          enddo
        enddo
      enddo

C --- store in agrp:  grad total rho . grad abs grad total rho ---
      if (nsp == 2) then
        do  ip = 1, np
        do  ir = ir0, ir1
          agrp(ir,ip,1) =
     .      (grp(ir,ip,1,1)+grp(ir,ip,1,2))*
     .      (gagrp(ir,ip,1,1)+gagrp(ir,ip,1,2)) +
     .      (grp(ir,ip,2,1)+grp(ir,ip,2,2))*
     .      (gagrp(ir,ip,2,1)+gagrp(ir,ip,2,2)) +
     .      (grp(ir,ip,3,1)+grp(ir,ip,3,2))*
     .      (gagrp(ir,ip,3,1)+gagrp(ir,ip,3,2))
        enddo
        enddo
      endif

C --- Copy grad rho . grad abs grad rho into gagrp ---
      do  i = 1, nsp
        do  ip = 1, np
          do  ir = ir0, ir1
            gagrp(ir,ip,i,1) = vxcp(ir,ip,i)
          enddo
        enddo
      enddo
      if (nsp == 2) then
        do  ip = 1, np
          do  ir = ir0, ir1
            gagrp(ir,ip,3,1) = agrp(ir,ip,1)
          enddo
        enddo
      endif
C     call px('gr.gagr',nri,nlm,1,np,ri(ir0),wp,gagrp(ir0,1,3,1),yl,wkl)

C --- Make agrp+,agrp- for ir0 .. ir1 ---
      do  i = 1, nsp
        do  ip = 1, np
          do  ir = ir0, ir1
            agrp(ir,ip,i) =
     .      dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
          enddo
        enddo
      enddo

C --- Make agrp (total rho) agrp+.agrp-  for ir0 .. ir1 ---
      if (nsp == 2) then
        do  ip = 1, np
        do  ir = ir0, ir1
          agrp(ir,ip,3) =
     .      dsqrt((grp(ir,ip,1,1)+grp(ir,ip,1,2))**2 +
     .            (grp(ir,ip,2,1)+grp(ir,ip,2,2))**2 +
     .            (grp(ir,ip,3,1)+grp(ir,ip,3,2))**2)
          agrp(ir,ip,4) =
     .             grp(ir,ip,1,1)*grp(ir,ip,1,2) +
     .             grp(ir,ip,2,1)*grp(ir,ip,2,2) +
     .             grp(ir,ip,3,1)*grp(ir,ip,3,2)
        enddo
        enddo
C       call px('x',nri,nlm,nsp,np,ri(ir0),wp,agrp(ir0,1,3),yl,wkl)
      endif

C --- Make nonlocal potential for points ir0 .. ir1 ---
      call dpzero(vxcp,nri*np*nsp)
      call dpzero(excp,nri*np)
      do  ip = 1, np
        if (lxcg > 2) then
          call vxcgga(lxcg,nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),
     .      agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1),
     .      ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),
     .      gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1),
     .      gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp),
     .      excp(ir0,ip))
        elseif (lcut == 0) then
        call vxnloc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1),
     .  agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp),
     .  agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1),
     .  gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),
     .  vxcp(ir0,ip,nsp),excp(ir0,ip))
        else
        call vxnlcc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),agrp(ir0,ip,1),
     .  agrp(ir0,ip,nsp),ggrp(ir0,ip,1),ggrp(ir0,ip,nsp),
     .  agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,4),gagrp(ir0,ip,2*nsp-1,1),
     .  gagrp(ir0,ip,1,1),gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),
     .  vxcp(ir0,ip,nsp),excp(ir0,ip))
        endif
      enddo

C ... (test): yl projection of various quantities'
C      call px('rho',nr,nlm,nsp,np,ri,wp,rp,yl,wkl)
C      call px('ggrh',nr,nlm,nsp,np,ri,wp,ggrp,yl,wkl)
C      call px('agrh',nri,nlm,nsp,np,ri(ir0),wp,agrp,yl,wkl)
C      call px('gr.gagr',nri,nlm,nsp,np,ri(ir0),wp,gagrp,yl,wkl)
C      call px('vxc',nri,nlm,nsp,np,ri(ir0),wp,vxcp,yl,wkl)
C      call px('exc',nri,nlm,nsp,np,ri(ir0),wp,excp,yl,wkl)

C --- Yl-projection of vxc,exc into vxcnl,excnl ---
      mlm = (lmax+1)**2
      call dgemm('N','N',nri,mlm,np,1d0,excp(ir0,1),nri,ylwp,np,0d0,
     .  wkl,nri)
      do  ilm = 1, mlm
        do  ir = ir0, ir1
          excnl(ir,ilm) = wkl(ir,ilm,1)
        enddo
      enddo
      do  i = 1, nsp
        call dgemm('N','N',nri,mlm,np,1d0,vxcp(ir0,1,i),nri,ylwp,np,0d0,
     .    wkl,nri)
        do  ilm = 1, mlm
          do  ir = ir0, ir1
            vxcnl(ir,ilm,i) = wkl(ir,ilm,1)
          enddo
        enddo
      enddo
C     call prmr(nr,ri,vxcnl,nlm)

C --- Estimate rmu in this ir0 ir1 interval by angular momentum ---
      do  i = 1, nsp
        do  l = 0, lmax
          sumnl(l) = 0d0
        enddo
        do  ilm = 1, nlm
          l = ll(ilm)
          do  ir = ir0, ir1
            sumnl(l) = sumnl(l) +
     .        rl(ir,ilm,i)*vxcnl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
          enddo
        enddo
        if (iprint() >= 80) then
          if (i == 1) print 1,ri(ir0),(sumnl(l),l = 0,lmax)
          if (i == 2) print 2,(sumnl(l),l = 0,lmax)
    1     format(' R>',f8.6,': ',f12.6,4F10.6:/(15x,4F10.6))
    2     format('     spin 2: ',f12.6,4F10.6:/(15x,4F10.6))
        endif
      enddo
C     call poppr

      end
C      subroutine px(strn,nr,nlm,nsp,np,ri,wp,fp,yl,fl)
C      implicit none
C      character *(*) strn
C      integer nr,nlm,nsp,np
C      double precision ri(1),wp(np),fp(nr,np,nsp),yl(np,nlm),
C     .  fl(nr,nlm,nsp)
C
C      call dpzero(fl,nr*nlm*nsp)
C      call fp2yl(nr,nlm,nsp,np,wp,fp,yl,0d0,fl)
C      print *, fl(1,1,1), fl(1,1,2)
C      print *, fl(591,1,1), fl(591,1,2)
C      print *, strn
C      call prmr(nr,ri,fl,nlm)
C
C      end
C
C      subroutine prmr(nr,ri,f,nl)
C      implicit none
C      integer nr,nl,ir,j,fopna,ifi,ir0
C      double precision ri(nr),f(nr,nl)
C      character*(10) fmt
C      ifi = fopna('out',19,0)
C      ir0 = 1
CC ... first nonzero l=0 point ...
C      do  20  ir = 1, nr
C        ir0 = ir
C        if (dabs(f(ir,1)) > 1d-12) goto 21
C   20 continue
C   21 continue
C
C      write(ifi,*) nr-ir0+1, nl+1
C      do  10  ir = ir0, nr
CC        write(ifi,333) ri(ir), (f(ir,3*j-2), j=1, nl)
C        write(ifi,333) ri(ir), (f(ir,j), j=1, nl)
CC 333   format(f12.7,(7g18.10:/12x))
C  333   format(f12.9,(9g18.10:/12x))
C   10 continue
C      call fclose(ifi)
C      print *, 'prmr:'
C      pause
C      end
C
