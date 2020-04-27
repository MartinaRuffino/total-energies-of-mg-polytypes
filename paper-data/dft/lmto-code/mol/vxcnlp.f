      subroutine vxcnlp(n1,n2,n3,nsp,alat,plat,vol0,ecut,rhoq,
     .  repnl,rmunl,vnl)
C- Gradient correction to smoothed rho(q) tabulated on a mesh
C  rho(q) input as rho(Q); returned as rho(R)
      implicit none

      integer n1,n2,n3,nsp
      double precision alat,plat(3,3),vol0,ecut,repnl(2),rmunl(2)
      double precision rhoq(n1+2,n2,n3,nsp),vnl(n1+2,n2,n3,nsp)
*     complex*16 rhoq(n1,n2,n3,1)
      integer n,nn,id,id2,np,ipr
      integer ogrq,oagrq,ogagrq,oggrq,oexc
      real w(1)
      common /w/ w

C --- Setup and memory allocation ---
      call tcn('vxcnlp')
      call getpr(ipr)
      id = n1/2+1
      id2 = n1+2
* reset id
*     id  = n1
*     id2 = n1
      np = 2*id*n2*n3
      n = n2*n3
      nn = 2*n
      call dpzero(repnl,2)
      call dpzero(rmunl,2)
      call defrr(ogrq,   np*3*nsp)
      call defrr(oagrq,  np*(3*nsp-2))
      call defrr(ogagrq, np*3*nsp)
      call defrr(oggrq,  np*nsp)
      call defrr(oexc,   np*nsp)

C --- Gradient-corrected potential ---
      call xxcnlp(n1,n2,n3,nsp,n,nn,id,np,alat,plat,vol0,ecut,rhoq,
     .  w(ogrq),w(oggrq),w(oagrq),w(ogagrq),repnl,rmunl,vnl,w(oexc))
      call rlse(ogrq)
      call tcx('vxcnlp')
      end

      subroutine xxcnlp(n1,n2,n3,nsp,n,nn,id,np,alat,plat,vol0,ecut,
     .  rhoq,grq,ggrq,agrq,gagrq,repnl,rmunl,vxc,exc)
C- Help subroutine for vxcnlp
      implicit none
      integer n1,n2,n3,n,nn,id,np,nsp
      double precision alat,plat(3,3),vol0,ecut,repnl(2),rmunl(2)
      double precision rhoq(np,nsp),grq(np,3,nsp),ggrq(np,nsp),
     .  agrq(np,*),gagrq(np,3,*),vxc(np,nsp),exc(np)
      integer ip,ix,i,lcut,lxcg

C --- grad rho and laplacian rho in Q-space ---
      call ropgrq(n1,n2,n3,nsp,n,nn,id,3,alat,plat,ecut,rhoq,grq,ggrq)
C      call zprx('grqx',grq,id,n2,n3)
C      call zprx('grqy',grq(1,2,1),id,n2,n3)
C      call zprx('grqz',grq(1,3,1),id,n2,n3)
C      call zprx('ggrq',ggrq,id,n2,n3)
*     call zprm3('grqx',grq,n1,n2,n3)
*     call zprm3('ggrq',ggrq,n1,n2,n3)

C --- Make rho(r) by backward FT rho(q) ---
      call rfft3(rhoq,n1,n2,n3,1,1)
C      call prm3('rho(r)',rhoq,n1+2,n1,n2,n3)

C --- Make grad rho(r) by backward FT grad rho(q) ---
      do  10  i = 1, nsp
      do  10  ix = 1, 3
   10 call rfft3(grq(1,ix,i),n1,n2,n3,1,1)
C      call prm3('grhox(r)',grq(1,1,1),n1+2,n1,n2,n3)
C      call prm3('grhoy(r)',grq(1,2,1),n1+2,n1,n2,n3)
C      call prm3('grhoz(r)',grq(1,3,1),n1+2,n1,n2,n3)

C --- Make laplacian rho(r) by backward FT nabla rho(q) ---
      call rfft3(ggrq,n1,n2,n3,1,1)
C     call prm3('ggrho(r)',ggrq,n1+2,n1,n2,n3)

C --- agrq : abs grad rho(r) ---
      do  20  i = 1, nsp
      do  20  ip = 1, np
   20 agrq(ip,i) = dsqrt(grq(ip,1,i)**2+grq(ip,2,i)**2+grq(ip,3,i)**2)
C     call prm3('abs grad rho',agrq,n1+2,n1,n2,n3)

C --- vxc (temp) : Forward FFT of agrq = abs grad rho ---
      call dpcopy(agrq,vxc,1,np,1d0)
      call rfft3(vxc,n1,n2,n3,1,-1)
C     call zprx('abs grad rho(q)',vxc,id,n2,n3)

C --- gagrq : grad abs grad rho(q) ---
      call ropgrq(n1,n2,n3,nsp,n,nn,id,1,alat,plat,ecut,vxc,gagrq,1d0)
C      call zprx('gagrqx',gagrq(1,1,1),id,n2,n3)
C      call zprx('gagrqy',gagrq(1,2,1),id,n2,n3)
C      call zprx('gagrqz',gagrq(1,3,1),id,n2,n3)

C --- gagrq : backward FT to make grad abs grad rho(r) ---
      do  12  i = 1, nsp
      do  12  ix = 1, 3
   12 call rfft3(gagrq(1,ix,i),n1,n2,n3,1,1)
C      call prm3('gagrrx(r)',gagrq(1,1,1),n1+2,n1,n2,n3)
C      call prm3('gagrry(r)',gagrq(1,2,1),n1+2,n1,n2,n3)
C      call prm3('gagrrz(r)',gagrq(1,3,1),n1+2,n1,n2,n3)

C --- vxc (temp) : grad rho . grad abs grad rho ---
      do  32  i  = 1, nsp
      do  32  ip = 1, np
   32 vxc(ip,i) =
     .    gagrq(ip,1,i)*grq(ip,1,i) +
     .    gagrq(ip,2,i)*grq(ip,2,i) +
     .    gagrq(ip,3,i)*grq(ip,3,i)
C     call prm3('grad rho . grad abs grad rho',vxc,n1+2,n1,n2,n3)

C --- agrq(3) (temp) grad total rho . grad abs grad total rho ---
C NB: abs grad rho+ + abs grad rho- used in place of abs grad total rho!
      if (nsp == 2) then
        do  34  ip = 1, np
          agrq(ip,3) =
     .      (grq(ip,1,1)+grq(ip,1,2))*
     .      (gagrq(ip,1,1)+gagrq(ip,1,2)) +
     .      (grq(ip,2,1)+grq(ip,2,2))*
     .      (gagrq(ip,2,1)+gagrq(ip,2,2)) +
     .      (grq(ip,3,1)+grq(ip,3,2))*
     .      (gagrq(ip,3,1)+gagrq(ip,3,2))
   34   continue
      endif

C --- gagrq : grad rho . grad abs grad rho (move from vxc)---
      do  36  i  = 1, nsp
      do  36  ip = 1, np
   36 gagrq(ip,i,1) = vxc(ip,i)

C --- gagrq(3) : mv from agrq(3) grad tot rho . grad abs grad tot rho
      if (nsp == 2) then
        do  38  ip = 1, np
   38   gagrq(ip,3,1) = agrq(ip,3)
      endif

C --- agrq(3) : abs grad tot rho;  agrq(4) : grad rho+ . grad rho- ---
      if (nsp == 2) then
        do  42  ip = 1, np
          agrq(ip,3) =
     .      dsqrt((grq(ip,1,1)+grq(ip,1,2))**2 +
     .            (grq(ip,2,1)+grq(ip,2,2))**2 +
     .            (grq(ip,3,1)+grq(ip,3,2))**2)
          agrq(ip,4) =
     .             grq(ip,1,1)*grq(ip,1,2) +
     .             grq(ip,2,1)*grq(ip,2,2) +
     .             grq(ip,3,1)*grq(ip,3,2)
   42   continue
      endif

C --- Nonlocal potential for all points  ---
      call dpzero(vxc,np*nsp)
      call dpzero(exc,np)
      if (lxcg() > 2) then
        call vxcgga(lxcg(),np,nsp,rhoq(1,1),rhoq(1,nsp),agrq(1,1),
     .    agrq(1,nsp),ggrq(1,1),ggrq(1,nsp),
     .    agrq(1,2*nsp-1),agrq(1,4),gagrq(1,2*nsp-1,1),
     .    gagrq(1,1,1),gagrq(1,nsp,1),vxc(1,1),
     .    vxc(1,nsp),exc)
      else
        lcut = 1
        if (lcut == 1) then
          call vxnlcc(np,nsp,rhoq(1,1),rhoq(1,nsp),agrq(1,1),
     .      agrq(1,nsp),ggrq(1,1),ggrq(1,nsp),
     .      agrq(1,2*nsp-1),agrq(1,4),gagrq(1,2*nsp-1,1),
     .      gagrq(1,1,1),gagrq(1,nsp,1),vxc(1,1),
     .      vxc(1,nsp),exc)
        else
          call vxnloc(np,nsp,rhoq(1,1),rhoq(1,nsp),agrq(1,1),
     .      agrq(1,nsp),ggrq(1,1),ggrq(1,nsp),
     .      agrq(1,2*nsp-1),agrq(1,4),gagrq(1,2*nsp-1,1),
     .      gagrq(1,1,1),gagrq(1,nsp,1),vxc(1,1),
     .      vxc(1,nsp),exc)
        endif
      endif

C --- Make nlocal reps, rmu ---
      do  50  i = 1, nsp
        call dpdot(rhoq(1,i),exc,np,repnl(i))
        call dpdot(rhoq(1,i),vxc(1,i),np,rmunl(i))
        repnl(i) = repnl(i)*vol0
        rmunl(i) = rmunl(i)*vol0
   50 continue

      end

      subroutine ropgrq(n1,n2,n3,nsp,n,nn,id,lgrd,alat,plat,ecut,rhoq,
     .  grq,ggrq)
C- Makes grad rho(q) or nabla rho(q) from rho(q)
C  lgrd: 1 or 3 make grad rho (store in grq);
C        2 or 3 make nabla rho (in ggrq)
      implicit none
      real w(1)
      common /w/ w

      integer n1,n2,n3,n,nn,id,nsp,lgrd
      double precision alat,plat(3,3),ecut
      integer oqp,orho,ogrho,oidx,oibp,obp,i1,nq,nbp
      double precision rhoq(n1+2,n2,n3,nsp),grq(n1+2,n2,n3,3,nsp),
     .  ggrq(n1+2,n2,n3,nsp)
*     complex*16 rhoq(n1,n2,n3,1),grq(n1,n2,n3,3,1),ggrq(n1,n2,n3,1)

      call defrr(oqp, nn*3)
      call defcc(orho,nn+1)
      call defcc(ogrho,3*(nn+1))
      call defrr(oidx,nn)
      call defrr(oibp,nn)
      call defrr(obp, nn*3)
      do  30  i1 = 1, id
        call gmeshx(alat,plat,0,ecut,n1,n2,n3,i1,nn,
     .    nq,nbp,w(oibp),w(obp),w(oidx),w(oqp))
        if (lgrd == 1 .or. lgrd == 3) then
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,0,w(orho),1,w(ogrho),rhoq)
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,1,w(orho),1,w(ogrho),grq)
        endif
        if (lgrd == 2 .or. lgrd == 3) then
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,0,w(orho),2,w(ogrho),rhoq)
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,1,w(orho),2,w(ogrho),ggrq)
        endif
        if (nsp == 2 .and. (lgrd == 1 .or. lgrd == 3)) then
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,0,w(orho),1,w(ogrho),rhoq(1,1,1,2))
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,1,w(orho),1,w(ogrho),grq(1,1,1,1,2))
        endif
        if (nsp == 2 .and. (lgrd == 2 .or. lgrd == 3)) then
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,0,w(orho),1,w(ogrho),rhoq(1,1,1,2))
          call zrhfxc(id,n,w(oidx),nn,nq,nbp,w(oibp),
     .      w(oqp),i1,1,w(orho),1,w(ogrho),ggrq(1,1,1,2))
        endif
   30 continue
      call rlse(oqp)
      end

