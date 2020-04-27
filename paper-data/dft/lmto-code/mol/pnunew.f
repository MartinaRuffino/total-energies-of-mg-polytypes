      subroutine pnunew(nbas,nclas,ipclas,nsp,n0,lmxa,ips,
     .   rmt,qus0,hab,sab,idmod,pnu,lfrz)
C- Shift pnu's to band cg
C 6 Jan 95 spin polarized (MvS)
      implicit real*8 (a-h,p-z), integer (o)
      dimension ipclas(1),ips(1),lmxa(1),idmod(n0,1),
     .  pnu(n0,1),qus0(3,n0,nbas,2),rmt(1),sab(4,n0,1),
     .  hab(4,n0,1),ql(20),hl(20),pnew(20)
      call tcn('pnunew')
      call getpr(ipr)
      pi=4d0*datan(1d0)
      nsp=lsp()+1

c ------ if frozen: return after printout ------
      if (lfrz == 1) then
         if (ipr >= 40) then
            write(6,411)
  411       format(/' pnunew: pnu''s are frozen'/
     .           ' class     pnu''s')
            do ic=1,nclas
              do jb=1,nbas
                if (ipclas(jb) == ic) ib=jb
              enddo
              is=ips(ib)
              kx=lmxa(is)+1
              write(6,410) ic,(pnu(k,ib),k=1,kx)
  410         format(i6,8f12.6)
            enddo
            write(6,*) ' '
         endif
         return
      endif

c ------ loop over classes: add up symmetrized chden coeffs ----
      do 10 isp=1,nsp
      do 10 ic=1,nclas
      nrclas=0
      do 1 k=1,20
      ql(k)=0d0
  1   hl(k)=0d0
      do 11 ib=1,nbas
C ibs is effective index to pnu(*,isp,ib), sab(*,isp,ib)
      ibs = nsp*(ib-1)+isp
      if(ipclas(ib) == ic) then
        nrclas=nrclas+1
        is=ips(ib)
        jb=ib
        lmax=lmxa(is)
        rmax=rmt(is)
        do 20 k=1,lmax+1
          ql(k)=ql(k)+qus0(1,k,ib,isp)*sab(1,k,ibs) +
     .                qus0(3,k,ib,isp)*sab(4,k,ibs) +
     .          0.5d0*qus0(2,k,ib,isp)*(sab(2,k,ibs)+sab(3,k,ibs))
          hl(k)=hl(k)+qus0(1,k,ib,isp)*hab(1,k,ibs) +
     .                qus0(3,k,ib,isp)*hab(4,k,ibs) +
     .          0.5d0*qus0(2,k,ib,isp)*(hab(2,k,ibs)+hab(3,k,ibs))
   20   continue
      endif
   11 continue

c ------ get dos center energy and new pnu ----------
      jbs = nsp*(jb-1)+isp
      if(ipr >= 40 .and. nsp == 1) write(6,312) ic,nrclas
      if(ipr >= 40 .and. nsp == 2) write(6,312) ic,nrclas,isp
  312 format(/' pnunew:   class',i3,'   nrclas=',i3:'   (spin',i2,')')
      if(ipr >= 40) write(6,311)
      do 12 k=1,lmax+1
        ql(k)=ql(k)/nrclas
        hl(k)=hl(k)/nrclas
        if(dabs(ql(k)) > 1d-8) then
          p1=2d10
          ebar=hl(k)/ql(k)
          huu=hab(1,k,jbs)
          hus=(hab(2,k,jbs)+hab(3,k,jbs))/2d0
          hss=hab(4,k,jbs)
          suu=sab(1,k,jbs)
          sus=(sab(2,k,jbs)+sab(3,k,jbs))/2d0
          sss=sab(4,k,jbs)
          a=hss-ebar*sss
          b=2d0*(hus-ebar*sus)
          c=huu-ebar*suu
c     ... first root of quad equation ...
          x1=(-b-dsqrt(b*b-4*a*c))/(2*a)
          q1=suu+2*x1*sus+x1*x1*sss
          h1=huu+2*x1*hus+x1*x1*hss
          p1=0.5d0-datan(rmax*x1)/pi
          if(ipr >= 80) write(6,341) k-1,x1,q1,h1,p1
  341     format(' l=',i2,' 1st root  x=',f11.5,'   q,h=',2f11.5,
     .      '   p=',f11.5)
c     ... second root of quad equation ...
          x2=(-b+dsqrt(b*b-4*a*c))/(2*a)
          q2=suu+2*x2*sus+x2*x2*sss
          h2=huu+2*x2*hus+x2*x2*hss
          p2=0.5d0-datan(rmax*x2)/pi
          if(ipr >= 80) write(6,342) x2,q2,h2,p2
  342     format(5x,' 2nd root  x=',f11.5,'   q,h=',2f11.5,
     .      '   p=',f11.5)
        endif

c ------ set new pnu (use first root) ----------------
        pold=pnu(k,jbs)
        ipqn=pold
        ptry=ipqn+p1
        pnew(k)=pold
        if(mod(idmod(k,is),10) == 0.and.dabs(p1) < 1d10) pnew(k)=ptry
        if(ipr >= 40) write(6,310) k-1,mod(idmod(k,is),10),ql(k),ebar,
     .    pold,ptry,pnew(k)
  310   format(i3,i6,f12.6,4f13.6)
  311   format('  l   idmod     ql',10x,'ebar',8x,' pold',9x,
     .    'ptry',9x,'pnew')
   12 continue

c ------ poke into array pnu -----------
      do 30 ib=1,nbas
      ibs = nsp*(ib-1)+isp
c|    write(6,888) ib,(pnu(k,ibs),k=1,lmax+1)
c|888 format(i5,5f11.6)
        if(ipclas(ib) == ic) then
          do 31 k=1,lmax+1
   31     pnu(k,ibs)=pnew(k)
        endif
c|    write(6,888) ib,(pnu(k,ibs),k=1,lmax+1)
   30 continue

   10 continue

      call tcx('pnunew')
      end
