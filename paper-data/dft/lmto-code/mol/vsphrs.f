      subroutine vsphrs(nbas,rmt,z,a,nr,lmxl,ips,orho,vval,ioffv0,
     .  cg,jcg,indxcg,cy,vnuc,zvn,qsum,asum,rep,rmu,rvh,sec,stc,lmxa,
     .  pnu,pnuz,lsc,n0,hab,vab,sab,ppnl,ioffp,puu,pus,pss,ov0,lfrz)
C- Sphere integrals and pot pars for both panels.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   rmt   :augmentation radius, in a.u., by species
Ci   z     :complex energy
Ci   z     :nuclear charge
Ci   a     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   lmxl  :l-cutoff for local density representation on a radial mesh
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   orho  :vector of local densities
Ci   vval  :vector of potentials at sphere boundary
Ci   ioffv0:offsets to vval for each site
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   cy    :Normalization constants for spherical harmonics
Ci   vnuc  :
Ci   lmxa  :augmentation l-cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnuz  :b.c. for core states
Ci   lsc   :1 if two-panel calculation
Ci   n0    :dimensioning parameter
Ci   hab
Ci   vab
Ci   sab
Ci   ppnl
Ci   ioffp
Ci   ov0   :potential
Ci   lfrz
Co Outputs
Ci   zvn   :
Ci   qsum
Ci   asum
Ci   rvh
Ci   rep
Ci   rmu
Ci   sec
Ci   stc
Ci   puu
Ci   pus
Ci   pss
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu  05 Apr 09 vsphrs calls vxcnsp for XC potential
Cu  17.05.04  added ppnl (NMTO pot. pars.)
Cu  10.01.95  spin polarized
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,n0,lsc,lfrz,nr(1),lmxl(1),ips(1),ioffv0(1),jcg(1),
     .  indxcg(1),lmxa(1),ioffp(1)
      integer orho(1),ov0(1)
      double precision rep,rmu,rvh,zvn,sec,stc,asum,qsum
      double precision rmt(1),a(1),z(1),cg(1),cy(1),vnuc(1),vval(1),
     .  pnu(n0,2,1),pnuz(n0,2,1),puu(1),pus(1),pss(1),
     .  hab(4,n0,1),sab(4,n0,1),vab(4,n0,1),ppnl(n0,9,2,nbas)
C ... Local parameters
C      No longer needed with vxcnsp call
C      integer nxl(0:7),nnn
C      integer oagrl,oagrp,oexc,oexcl,oexcnl,oggrp,ogrp,ogyl,
C     .  or2,orp,ovxc,ovxcnl,oxp,oyl,oylwp,oyp,ozp
C      parameter( nnn=144 )
C      double precision p(3,nnn),wp(nnn)
C      integer np,nph,nx
      double precision pp(n0,2,5)
      integer i0,i1,ib,ibsp,ipan,ipr,is,l,lmaxa,lmaxl,lsp,lxcf,lxcfun,
     .  lxcg,ncore,nlma,nlml,npan,nr1,nsp,intopt,nglob
      integer orl,ov00,orofi,orwgt,owk,ovl,ovdif,ovx,ovxcl,og,ogp,ogpp,
     .  opz0,osl,oul,ogz,orss,orus,oruu
      double precision y0,r1,a1,z1,b,rhoves,rhves1,vnucl,vsum,edif
      double precision sum(0:20),q1(2),rep1(2),rmu1(2),sme(2),smt(2),
     .  repsx(2),repsc(2),qcor(2)
C     For local orbitals ..
      double precision ehl(n0),rsml(n0),rs3,vmtz
      integer w(1)
      common /w/ w
C     data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

C --- Setup ---
      call getpr(ipr)
      call tcn('vsphrs')
      qsum = 0d0
      asum = 0d0
      rep = 0d0
      rmu = 0d0
      rvh = 0d0
      zvn = 0d0
      sec = 0d0
      stc = 0d0
      nsp = lsp()+1
      lxcfun = lxcf()
      y0 = 1/dsqrt(16*datan(1d0))

c ------ start loop over atoms ---------------------
      do 10 ib = 1,nbas
      is = ips(ib)
      lmaxl = lmxl(is)
      lmaxa = lmxa(is)
      nlml = (lmaxl+1)**2
      r1 = rmt(is)
      nr1 = nr(is)
      a1 = a(is)
      z1 = z(is)
      q1(2) = 0d0
      rep1(2) = 0d0
      rmu1(2) = 0d0
      orl = orho(ib)
      ov00 = ov0(ib)
      call defrr(orofi,       nr1)
      call defrr(orwgt,       nr1)
      call defrr(owk,         nr1*nsp)
      call defrr(ovl,         nr1*nlml*nsp)
      call defrr(ovx,         nr1*nsp)
      call defrr(ovdif,       nr1*nsp)
      call defrr(orl,         nr1*nlml*nsp)
C ... Total charge into rl for Poisson's equation
      if (nsp == 1) then
        call dpcopy(w(orho(ib)),w(orl),1,nr1*nlml,1d0)
      else
        call dpscop(w(orho(ib)),w(orl),nr1*nlml,1+nr1*nlml,1,1d0)
        call dpadd(w(orl),w(orho(ib)),1,nr1*nlml,1d0)
      endif
      call radmsh(r1,a1,nr1,w(orofi))
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,r1,a1,nr1,w(orwgt))
      b = r1/(dexp(a1*nr1-a1)-1)
      if(ipr >= 30) write(6,725) ib,nlml,nr1,r1,z1
  725 format(/' vsphrs,  atom',i3,' :   nlml=',i3,'   nr=',i4,
     .  '   rmt=',f9.5,'   z=',f7.1)

C --- Solve poisson equation in sphere ---
      i0 = ioffv0(ib)+1
      call poinsp(z1,vval(i0),nlml,a1,b,w(ovl),w(orofi),w(orl),
     .   w(owk),nr1,rhoves,rhves1,vnucl,vsum)
      vnuc(ib) = vnucl
      zvn = zvn-z1*vnucl
      call vrobyl(nr1,nlml,w(ovl),w(orl),w(orofi),w(orwgt),z1,
     .  sum,rhoves,001)

C --- X-C potential and integrals ---
C ... copy ves into spin-down; copy rho(up,down) into rho
      if (nsp == 2) then
        call dpscop(w(ovl),w(ovl),nr1*nlml,1,1+nr1*nlml,1d0)
        call dpcopy(w(orho(ib)),w(orl),1,nr1*nlml*nsp,1d0)
      endif

      call defrr(ovxcl,   -nr1*nlml*nsp)

C --- XC potential into vxcl ---
      call vxcnsp(0,w(orofi),nr1,w(orwgt),nlml,nsp,w(orl),
     .  lxcf()+100*lxcg(),w(1),w(1),w(1),w(1),w(1),rep1,repsx,repsc,
     .  rmu1,w(ovxcl),w(1),q1)

C --- Collect energy integrals ---
      if(ipr >= 30 .and. nsp == 1) write(6,955) rhoves,vnucl
      if(ipr >= 30 .and. nsp == 2) write(6,955) rhoves,vnucl,q1(2)-q1(1)
  955 format(' rhoves=',f15.6,'    vnucl=',f14.6:'    mag. mom=',f10.6)
      if(ipr >= 35) write(6,341) (sum(l),l=0,lmaxl)
  341 format(' rho*ves by l: ',f13.6,4f10.6:/(17x,4f10.6))
      qsum = qsum+q1(1)+q1(2)
      asum = asum+q1(2)-q1(1)
      rvh = rvh+rhoves
      rep = rep+rep1(1)+rep1(2)
      rmu = rmu+rmu1(1)+rmu1(2)

C --- Sum pot, maybe get difference to frozen pot ---
      call dpadd(w(ovl),w(ovxcl),1,nr1*nlml*nsp,1d0)
      call rlse(ovxcl)
      call dpcopy(w(ovl), w(ovx), 1, nr1, y0)
      if (nsp == 2) call dpscop(w(ovl),w(ovx),nr1,1+nr1*nlml,1+nr1,y0)
      if (lfrz == 0) call dpcopy (w(ovx), w(ov00), 1, nr1*nsp, 1d0)
      call dpcopy(w(ovx), w(ovdif), 1, nr1*nsp, 1d0)
      call dpadd(w(ovdif), w(ov00), 1, nr1*nsp, -1d0)

C --- Calculate core ---
C     call prmx('v00 in vsphrs',w(ov00),nr(is),nr(is),nsp)
      sme(2) = 0
      smt(2) = 0
C ibs is effective index to pnu(1,isp,ib) and hab(1,isp,ib), isp=1
C     ibs = nsp*(ib-1)+1
C --- call getcore from lmf, single panel only
      call defrr(opz0,  -10*nsp)
      qcor(1) = 0
      qcor(2) = 0
      call getcor(0,z1,a1,pnu(1,1,ib),w(opz0),nr1,lmaxa,w(orofi),
     .            w(ov00),0,-1,qcor,sme,smt,w(owk),ncore,0d0,0d0)
      call rlse(opz0)
      call ecrdif(nr1,nsp,w(owk),w(ovdif),w(orofi),w(orwgt),edif)
      if(ipr >= 30) write(6,771) sme(1)+sme(2),edif,sme(1)+sme(2)+edif
  771 format(' sum(ec0)=',f14.6, '   edif=',f12.6,'    sum(ec)=',f14.6)
      sme(1) = sme(1)+sme(2)+edif
      sec = sec+sme(1)
      stc = stc+smt(1)+smt(2)

C --- Make pot pars and pert matrices for each panel ---
C ibsp is effective index to puu(1,isp,ib,ipan) for isp=1
      npan = lsc+1
      do  20  ipan = 1, npan
        ibsp = nsp*(ib-1) + nsp*nbas*(ipan-1) + 1
        i1   = nsp*ioffp(ib) + nsp*ioffp(nbas+1)*(ipan-1) + 1
        call defrr(og,     2*nr1)
        call defrr(ogp,    2*nr1)
        call defrr(ogpp,   2*nr1)
C  ...  call potpus from lmf, single panel only, lrel in common /cc/
        call defrr(opz0,  -10*nsp)
        call potpus(z1,r1,lmaxa,w(ov00),w(ovdif),a1,nr1,nsp,0,w(orofi),
     .    pnu(1,1,ib),w(opz0),ehl,rsml,rs3,vmtz,4,n0,pp,ppnl(1,1,1,ib),
     .    hab(1,1,ibsp),vab(1,1,ibsp),sab(1,1,ibsp),w)
        call rlse(og)
C       call prmx('hab in vsphrs',hab(1,ibsp),4,4,1+lmaxa)
C       call prmx('sab in vsphrs',sab(1,ibsp),4,4,1+lmaxa)

c   ... make pertubation matrices
        nlma=(lmaxa+1)**2
        call defrr(oul,   nr1*(lmaxa+1)*nsp)
        call defrr(osl,   nr1*(lmaxa+1)*nsp)
C   ... call makusp from lmf, single panel only, lrel in common /cc/
        call defrr(opz0,  -10*nsp)
        call defrr(ogz,   nr1*(lmaxa+1)*nsp)
        call defrr(oruu,  nr1*(lmaxa+1)*2*nsp)
        call defrr(orus,  nr1*(lmaxa+1)*2*nsp)
        call defrr(orss,  nr1*(lmaxa+1)*2*nsp)
        call makusp(n0,z1,nsp,1,r1,lmaxa,w(ov00),a1,nr1,
     .    0d0,0d0,pnu(1,1,ib),w(opz0),w,w,w(oul),w(osl),w(ogz),w(oruu),
     .    w(orus),w(orss))
        call rlse(opz0)
        call pertms(w(oul),w(osl),lmaxa,nlma,nlma,w(ovl),nlml,w(orofi),
     .    a1,nr1,puu(i1),pus(i1),pss(i1),cg,jcg,indxcg,nsp)
   20 continue

      call rlse(orofi)
  10  continue
      call tcx('vsphrs')
      end
      subroutine vrobyl(nr1,nlml,vl,rl,rofi,rwgt,z,sum,tot,job)
      implicit real*8 (a-h,p-z), integer (o)
      dimension sum(0:10),rl(nr1,1),vl(nr1,1),rofi(1),rwgt(1)
      lmaxl=ll(nlml)
      do 33 l = 0,lmaxl
  33  sum(l) = 0d0
      if(job == 1) then
        xx = 2d0*z*dsqrt(16d0*datan(1d0))
        do 3 ir = 2,nr1
  3     sum(0) = sum(0)+rl(ir,1)*(vl(ir,1)-xx/rofi(ir))*rwgt(ir)
      else
        do 4 ir = 2,nr1
  4     sum(0) = sum(0)+rl(ir,1)*vl(ir,1)*rwgt(ir)
      endif
      do 30 ilm = 2,nlml
      l = ll(ilm)
      do 30 ir = 2,nr1
  30  sum(l) = sum(l)+ rl(ir,ilm)*vl(ir,ilm) * rwgt(ir)
      tot = 0d0
      do 35 l = 0,lmaxl
  35  tot = tot+sum(l)
      end
      subroutine xxph21(nr,nlm,nsp,ri,rl,isgn)
      implicit none
      integer nr,nlm,nsp,isgn,i,ilm,ir
      double precision rl(nr,nlm,nsp),ri(nr),rho2,rho3

C --- Scale rho by 1/r**2 ---
      if (isgn == 1) then
        do  10  i = 1, nsp
        do  10  ilm = 1, nlm
        rl(1,ilm,i) = 0d0
        do  10  ir = 2, nr
   10   rl(ir,ilm,i) = rl(ir,ilm,i)/ri(ir)**2
C  ...  Extrapolate rho to origin
        do  20  i = 1, nsp
        rho2 = rl(2,1,i)
        rho3 = rl(3,1,i)
   20   rl(1,1,i) = (rho2*ri(3)-rho3*ri(2))/(ri(3)-ri(2))
C  20   rl(1,1,i) = 0
      else
        do  110  i = 1, nsp
        do  110  ilm = 1, nlm
        rl(1,ilm,i) = 0d0
        do  110  ir = 2, nr
  110   rl(ir,ilm,i) = rl(ir,ilm,i)*ri(ir)**2
      endif
      end

      subroutine ecrdif(nr,nsp,rhoc,vdif,rofi,rwgt,edif)
c make core energy shift due to difference potential
      implicit real*8 (a-h,p-z), integer (o)
      dimension rhoc(nr,nsp),vdif(nr,nsp),rofi(2),rwgt(2)
      edif = 0d0
c|    q = 0d0
      do isp = 1,nsp
      do ir = 1,nr
c|      q = q+rwgt(ir)*rhoc(ir,isp)
        edif = edif+rwgt(ir)*rhoc(ir,isp)*vdif(ir,isp)
      enddo
      enddo
c|    write(6,300) edif,q
c|300 format('ecrdif:  edif=',f12.6,'   qcr=',f12.6)
      end
