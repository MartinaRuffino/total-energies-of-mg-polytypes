      subroutine hy2gen(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,np1,lp2,ep2,np2,
     .  n0,llink,el1,tl1,el2,tl2,rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,
     .  d0,ndist,adec,wx,nalf,nb1,nb2,ri1,ri2,ncupl,ndust,ifi)
C- Global fit to all distances simultaneously, linked basis
C  The phi's and xi's are smoothed Hankel functions,
C  lp,ep,np describe the phi-basis; lx,ex,nx describe the xi-basis.
C  rmt are the muffin-tin radii for which the fit is designed
C  but do not influence the fit itself.
C  HOLINT determines integration mesh;  requires ri1,ri2,nb1,nb2
C         Fit is done for the volume outside RI1 and RI2.
C  HYFSCH generates fit functions.
C  Output: TCF info and data on file IFI.
      implicit none
      integer ndim(0:8,0:10),nrhs(0:20),ndimb(0:20),lx1(8),lx2(8),llink,
     .  ncupl,ndust
      integer lp1(1),lp2(1),lsym,mmax,nb1,nb2,nalf,ndimbx,ndimx,n0,
     .  ndist,nrhsx,nsdmbx,nsdmx,nx1,nx2,np1,np2,ifi
      double precision ep1(1),ep2(1),rmt1,rmt2,zc1,zc2,ri1,ri2,
     .  rsm1,rsm2,rsmp1,rsmp2,ep1j,ep2j,add1(25,5),add2(25,5),
     .  ex1(8),ex2(8),dist(80),wist(80),err(100),cx(300),
     .  el1,tl1(n0,5),el2,tl2(n0,5),adec,wx,d0,d
      integer oph1,oph2,ob,of,os,ocof,obot,otop,oba,obb,owk,owk2,owk3,
     .  ofit,oprd,onrhsj,oerr1,oerr2,osb,oxi,owp,oxp,ozp,oserr,osm,i,j,
     .  idist,lmxp1,lmxp2,nlmp1,nlmp2,ipr,irep,lmax,nbisi,
     .  ncof,nlm,np,nrep,nx,nxi,idxj(2,20),lsymj(20),nph,nhl,nelx,
     .  lp1j,lp2j,iph,modef,lscal,lx(20)
C     double precision dval,alfaj(20),alfph(2,20)

      integer ndCRAY
C#ifndef CRAY-DP
      parameter (ndCRAY=1)
C#elseC
C      parameter (ndCRAY=2)
C#endif
      real w(1)
      common /w/ w
      data lsymj /20*0/

C --- Setup ---
      call getpr(ipr)
      lscal = 1
      if (ipr >= 20)
     .  print 408, ri1,ri2,ndist,nb1,nb2,lscal == 1,ncupl,ndust
  408 format(/' hy2gen: ri=',2f7.3,'  ndist=',i2,
     .  '  nb1,nb2=',2i3,'  scale=',l1,'  ncupl=',i2,'  ndust=',i2)
      call sxlmnc(cx,10)
      if (ndist > 80) call rx('HY2GEN: ndist too big')
      nbisi = 2*nb1*nb2
      np = nbisi
      mmax = 0

      if (llink == 1) then
C        print 333, el1, el2
C  333   format(' linked basis,  elink=',f8.4,' (site 1)  ',4f8.4,
C     .    '  (site 2)')
C       call msadd(n0,np1,rsmp1,rmt1,llink,el1,tl1,lp1,ep1,add1)
        call msetpl(n0,np1,' 1  ',1,1,rmt1,np1,1,llink,
     .    el1,tl1,lp1,ep1,nelx,nhl,add1)
C        call msadd(n0,np2,rsmp2,rmt2,llink,el2,tl2,lp2,ep2,add2)
        call msetpl(n0,np2,' 2  ',1,1,rmt2,np2,1,llink,
     .    el2,tl2,lp2,ep2,nelx,nhl,add2)
      endif

      do  12  i = 1, np1
      do  12  j = 1, np2
   12 mmax = max(lp1(i)+lp2(j),mmax)
      call defrr(onrhsj, (mmax+1)*(np1*np2+1))
      call hy2dm(mmax,lp1,ep1,np1,lx1,ex1,nx1,lp2,ep2,np2,lx2,ex2,nx2,
     .  rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,el1,tl1,el2,tl2,
     .  ndim,nrhs,w(onrhsj),ndimx,nsdmx,nrhsx,ncof,nxi,idxj,lsymj,nph,
     .  nalf,ndimb,ndimbx,nsdmbx)
      call defrr(of,     ndist*nalf)
      call hyfsch(ndist,nalf,d0,adec,wx,dist,wist,w(of))

C --- Make b and s for each idist ---
      call defrr(ob,     ndist*ndimx*nrhsx*(mmax+1))
      call defrr(os,     ndist*nsdmx*(mmax+1))
      call hy2nrm(mmax,lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,lscal,rmt1,rmt2,
     .  ri1,ri2,rsmp1,rsmp2,llink,el1,add1,el2,add2,rsm1,rsm2,
     .  lp1,ep1,np1,lp2,ep2,np2,idxj,nph,lsymj,w(os),w(ob),
     .  ndimx,nsdmx,nrhsx,dist,wist,ndist)

C --- Set up and solve hyperfit ---
C     call hysnit(ncupl)
      print *, 'hy2gen reduce dim of of osb'
      call defrr(obb,    nalf*ndimx*nrhsx*(mmax+1))
      call dpzero(w(obb),nalf*ndimx*nrhsx*(mmax+1))
      call defrr(osb,    nsdmx*ncupl*ncupl*ndCRAY * 2)
      call defrr(osm,    ndimx*ndimx)
      call defrr(owk,    ndimbx*nrhsx*ndCRAY)
      call defrr(owk2,   ndimx*max(nrhsx,ncupl)*ndCRAY)
      call defrr(owk3,   ndimx*max(nrhsx,ncupl)*ndCRAY)
      call defrr(oerr1,  nrhsx*2)
      call hyfsbx(ncupl,ndust,ndimbx,nrhsx,mmax,nsdmbx,nsdmx,ndimx,ndim,
     .  nrhs,w(ob),w(os),ndist,nalf,w(of),wist,w(owk),w(owk2),w(owk3),
     .  w(oerr1),w(osm),w(osb),w(obb))
C ... Clean up memory by copying bb to b, reallocating b, renaming
      call dpcopy(w(obb),w(ob),1,nalf*ndimx*nrhsx*(mmax+1),1d0)
      call rlse(ob)
      call defrr(ob,     nalf*ndimx*nrhsx*(mmax+1))
      obb = ob
      call hyfusc(mmax,ndim,nrhs,nx1,lx1,ex1,nx2,lx2,ex2,
     .  nalf,rmt1,rmt2,w(obb),ndimx,nrhsx)

C --- Setup for checking fit ---
      call tcn('hy2gen')
      nrep = min(3,ndist)
      if (ipr >= 50) nrep = ndist
      call defrr(oserr,      nrhsx*(mmax+1))
      call dpzero(w(oserr),  nrhsx*(mmax+1))
      if (ipr < 30) nrep = 0
      do  70  irep = 1, nrep
      if (irep == 1) idist = 1
      if (irep == 2) idist = (ndist+1)/2
      if (irep == 3) idist = ndist
      if (ipr >= 50) idist = irep
      d = dist(idist)
      call defrr(ob,     ndimx*nrhsx*(mmax+1))
C      call hyfip2(d,d0,modef,nalf,nph,alfph,a1,mmax,
C     .  ndimx,nrhsx,ndim,w(onrhsj),nx1,w(obb),w(ob))
      call hyfipt(d,d0,adec,nalf,nph,mmax,
     .  ndimx,nrhsx,ndim,w(onrhsj),nx1,w(obb),w(ob))
      call hy2xim(lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
     .  ri1,ri2,rsmp1,rsmp2,llink,el1,add1,el2,add2,
     .  rsm1,rsm2,rmt1,rmt2,0,0,
     .  np,nbisi,lp1,ep1,np1,lp2,ep2,np2,cx,d,
     .  lmax,lmxp1,lmxp2,nx,lx,zc1,zc2,oph1,oph2,owp,oxi,oxp,ozp)

      nlm = ((lmax+1)*(lmax+2))/2
      nlmp1 = ((lmxp1+1)*(lmxp1+2))/2
      nlmp2 = ((lmxp2+1)*(lmxp2+2))/2

C --- Check fit ---
      call defrr(oerr1,      nrhsx*(mmax+1))
      call defrr(oerr2,      nrhsx*(mmax+1))
      call defrr(oprd,       nrhsx*(mmax+1))
      call defrr(obot,       nrhsx*(mmax+1))
      call defrr(otop,       nrhsx*(mmax+1))
      call defrr(owk,        np)
      call defrr(ofit,       np)
      call dpzero(err,       100)
      if (ipr >= 30)
     .  call tcferr(d,mmax,w(owp),np,w(owk),w(ofit),
     .  nx,lx,w(oxi),nlm,lmax,w(oph1),w(oph2),nlmp1,nlmp2,
     .  lp1,ep1,lp2,ep2,idxj,nph,lsymj,w(ob),
     .  ndimx,nrhsx,w(oerr1),w(oerr2),w(oprd),w(otop),err)
      if (ipr >= 50) call dpadd(w(oserr),w(oerr2),1,
     .  nrhsx*(mmax+1),wist(idist))
      call rlse(ob)
   70 continue
      call rlse(oserr)

        print *, 'ZERO err'
        CALL DPZERO(ERR,       100)

C --- For each energy pair, sort coefficients into final form; save ---
      call defrr(ocof,   ncof*nalf)
      call defrr(oba,    ndimx*nrhsx*(mmax+1)*nalf)
      do  80  iph = 1, nph
        lp1j = lp1(idxj(1,iph))
        lp2j = lp2(idxj(2,iph))
        ep1j = ep1(idxj(1,iph))
        ep2j = ep2(idxj(2,iph))
        lsym = lsymj(iph)
C ...   Sort coefficients into cof by ialf-columns
        call hyfcof(mmax,lx1,nx1,lx2,nx2,iph,lp1j,lp2j,lsym,ndim,
     .    w(onrhsj),w(oba),w(obb),ndimx,nrhsx,ndimbx,nalf,w(ocof),ncof)
        call dscal(ncof,2d0,w(ocof),1)
C ...  Output into file ifi
        call hyfout(rmt1,rmt2,rsm1,rsm2,d0,adec,adec,ndist,nalf,mmax,
     .    ncof,nb1,nb2,ri1,ri2,nxi,lp1j,ep1j,lp2j,ep2j,lsym,
     .    lx1,ex1,nx1,lx2,ex2,nx2,w(ocof),err,ifi)
   80 continue
      call rlse(onrhsj)
      call tcx('hy2gen')
      end
      subroutine hyfcof(mmax,lx1,nx1,lx2,nx2,iph,lp1j,lp2j,lsym,
     .  ndim,nrhsj,ba,bb,ndimx,nrhsx,ndimbx,nalf,cofa,ncof)
C- Sort into columns by ialf
      implicit none
      integer mmax,lx1(1),nx1,lx2(1),nx2,lp1j,lp2j,lsym,iph,
     .  ndimx,nrhsx,ndimbx,nalf,ncof,ndim(0:20),nrhsj(0:mmax,1)
      double precision cofa(1),
     .  ba(ndimx,nrhsx,0:mmax,nalf),bb(ndimx,nrhsx,0:mmax,nalf)
      integer ialf,id,ir,m,nd,ir1,nr,irhs(0:20)

      ncof = 0
C --- Copy coefficients into ba ---
      do  10  m = 0, mmax
        nd = ndim(m)
        ir1 = nrhsj(m,iph)
        nr = nrhsj(m,iph+1) - ir1
        irhs(m) = nr
        ncof = ncof + ndim(m)*nr
        do  11  ialf = 1, nalf
        do  11  ir = 1, nr
        do  11  id = 1, nd
   11   ba(id,ir,m,ialf) = bb(id,ir+ir1,m,ialf)
   10 continue

C --- Call tcfsrt for each ialf ---
      do  20  ialf = 1, nalf
        call tcfsrt(mmax,lx1,nx1,lx2,nx2,lp1j,lp2j,lsym,ndim,irhs,
     .    ba(1,1,0,ialf),ndimx,nrhsx,cofa(1+ncof*(ialf-1)))
   20 continue

      end
      subroutine hysnit(ncupl)
      integer ncupl
C#ifdefC TEST
C      logical cmdopt,a2bin
C      character strn*120
C#endif
      ncupl = 4
C#ifdefC TEST
C      if (cmdopt('-ncupl=',7,0,strn)) then
C        j=7
C        call rxx(.not. a2bin(strn,ncupl,2,0,' ',j,40))
C      endif
C#endif
      end
