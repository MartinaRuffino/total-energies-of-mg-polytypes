      subroutine vtailz(rhoi,poti,pot0,nxi,lxi,exi,n0,rmt,rint,rcut,
     .  rsm,lmxl,nbas,nbasp,lmaxp,alat,pos,ips,ioff,ioffv0,cg,jcg,
     .  indxcg,cy,nri,ixp,dxp,nvi,vval,zeta,zetxc,sxi0,qsmo,asmo,repsmo,
     .  rmusmo,qmom,upm,f1,f2,f4)
C- Contribution to zeta,zetxc from smooth potential tails
C ----------------------------------------------------------------------
Ci Inputs
Ci   rhoi  :coefficients to density in atom-centered density basis
Ci   poti  :coefficients to estat potential in atom-centered density basis
Ci   pot0  :coefficients to estat potential, homogeneous contribution
Ci   nxi   :number of functions per l for c.d. basis, by species
Ci   lxi   :l cutoffs for each function in c.d. basis, by species
Ci   exi   :Hankel energies for each function in c.d. basis, by species
Ci   n0    :dimensioning parameter
Ci   rmt   :augmentation radius, in a.u., by species
Ci   rint  :range of interstial charge density
Ci   rcut  :radius for which to correct artificial smoothing in
Ci         :construction of xc potential
Ci   rsm   :smoothing radius
Ci   lmxl  :l-cutoff for local density representation on a radial mesh
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (including classical multipoles)
Ci   lmaxp :l-cutoff for classical multipoles
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   pos   :basis vectors
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ioff  :table of offsets to rhoi and poti for each site
Ci   ioffv0:table of offsets to vval and pot0 for each site
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   cy    :Normalization constants for spherical harmonics
Ci   nri   :total number of functions in cd basis: dimensions of rhoi
Ci   ixp   :0 for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   dxp   :alat for molecules branch
Ci         :Otherwise, parameters related to Ewald summations
Ci   nvi   :dimensions of vval and pot0
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Co Outputs
Co   vval  :electrostatic potential at MT boundaries
Co   zeta  :integrals of density functions with electrostatic potential
Co   zetaxc:integrals of density functions with xc potential
Co   sxi0  :amount of istl charge in augmentation sphere?
Co   qsmo  :
Co   asmo  :
Co   repsmo:
Co   rmusmo:
Co   upm   :interaction energy between total electrostatic potential
Co         :and classical point multipoles (sites ib>nbas)
Co   f1    :force from integral( drho * eps ) + istl ring corrections
Co   f2    :force from drho*ves + rho*dves + terms from vval,qval
Co   f4    :force from derivative of istl charge
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   15 Jun 09 (ATP) parallelized
Cu   05 Apr 09 updated call to vxcnls
Cu   10.01.95  spin polarized
Cu   02.02.95  gradient corrections added
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n0,nbas,nbasp,lmaxp,nri,nvi,nxi(1),ips(1),lxi(n0,1),
     .  ioff(1),ioffv0(1),lmxl(1),jcg(1),indxcg(1),ixp(1)
      double precision alat,asmo,qsmo,repsmo,rmusmo,upm
      double precision exi(n0,1),pos(3,1),rmt(1),rint(1),cg(1),
     .  cy(1),poti(1),pot0(1),rhoi(1),vval(1),zeta(1),zetxc(1),rsm(1),
     .  sxi0(1),f1(3,nbasp),f2(3,nbasp),f4(3,nbasp),qmom(1),
     .  rcut(1),dxp(1)
C ... Local parameters
      integer i,i0,i1,ib,ib1,ib2,ie,ipr,ir,is,ix,lmax,lmaxl,lmaxx,lsp,
     .  lxcf,lxcfun,lxcg,m,nlml,nlmx,np,nph,nr,nr1,nr2,nrwk,nsp,nx,nnn,
     .  nr12,lgunit,stdo,nxl(0:7),intopt,nglob
      parameter(nnn=144, nr12=162)
      integer oagrl,oagrp,odrl2,odrp,odvxc2,oexc,oexcl,oggrp,ogrp,ogyl,
     .  or2,orhoi,orl1,orl2,orlt,orlx,orp,ovl,ovxc,ovxc1,ovxc2,ovxcl,
     .  owk1,owk2,oxi,oxi0,oxp,oyl,oylwp,oyp,ozp
      double precision a1,b1,r1,r2,vmad,vmtail,y0
      double precision p(3,nnn),rofi(nr12),wp(nnn),rwgt(nr12),
     .  rep(2),rmu(2),qsm(2),rep2(2),rmu2(2),q2(2),repsx(2),repsc(2)
C ... Heap
      integer w(1)
      common /w/ w
      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      double precision, dimension(:), allocatable :: ztxcbk,f1bk
      integer mpipid,procid,pid,master,numprocs,length,ierr
      integer i1mach,iprint
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      stdo = lgunit(1)
      call getpr(ipr)
      call tcn('vtailz')
      a1 = 0.01d0
      if (lxcg() == 0) then
        nr1 = 51
        nr2 = 51
      else
        nr1 = 81
        nr2 = 81
      endif
      intopt = 10*nglob('lrquad')
      y0 = 1d0/dsqrt(16d0*datan(1d0))
      qsmo = 0d0
      asmo = 0d0
      repsmo = 0d0
      rmusmo = 0d0
      upm  = 0
      call dpzero(sxi0,   nri)
      call dpzero(vval,   nvi)
      nr = nr1+nr2
      nsp = lsp()+1
      lxcfun = lxcf()
C ... Make rhoi+ + rhoi-
      call defrr(orhoi, nri)
      call dpcopy(rhoi,w(orhoi),1,nri,1d0)
      if (nsp == 2) call dpadd(w(orhoi),rhoi(1+nri),1,nri,1d0)

C --- start loop over atoms. Setup, make meshes ---
      if (MPI) then
        allocate (ztxcbk(1:nri*nsp), stat=ierr)
        allocate (f1bk(1:3*nbasp), stat=ierr)
C --- buffer zetxc and f1 and set locally to zero --
        call dcopy(nri*nsp,zetxc,1,ztxcbk,1)
        call dcopy(nri*nsp,0d0,0,zetxc,1)
        call dcopy(3*nbasp,f1,1,f1bk,1)
        call dcopy(3*nbasp,0d0,0,f1,1)
        allocate (bproc(0:numprocs), stat=ierr)
        call dstrbp(nbas,numprocs,1,bproc(0))
        ib1 = bproc(procid)
        ib2 = bproc(procid+1)-1
      else
        ib1 = 1
c...sl
c       ib2 = nbas
        ib2 = nbasp
c...sl
      endif
      do  10  ib = ib1, ib2
      if (MPI .and. mlog) then
        if (ib == bproc(procid)) then
        call gettime(datim)
        call awrit4(' vtailz '//datim//' Process %i of %i on '
     .    //shortname(1:length)//
     .    ' starting atoms %i to %i',' ',256,lgunit(3),
     .    procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
      endif
      qsm(1) = 0
      rmu(1) = 0
      rep(1) = 0
      q2(1) = 0
      rmu2(1) = 0
      rep2(1) = 0
      qsm(2) = 0
      rmu(2) = 0
      rep(2) = 0
      q2(2) = 0
      rmu2(2) = 0
      rep2(2) = 0
      if (ib <= nbas) then
        is = ips(ib)
        lmaxl = lmxl(is)
        lmaxx = -1
        do  4  ie = 1, nxi(is)
    4   lmaxx = max0(lmaxx,lxi(ie,is))
        r1 = rmt(is)
        r2 = rcut(is)
        call radmsh(r1,a1,nr1,rofi)
        call radwgt(intopt,r1,a1,nr1,rwgt)
        call mklegw(nr2,rofi(nr1+1),rwgt(nr1+1),000)
        do  5  ir = 1, nr2
          rofi(nr1+ir) = (r2+r1+(r2-r1)*rofi(nr1+ir))*0.5d0
          rwgt(nr1+ir) = 0.5d0*(r2-r1)*rwgt(nr1+ir)
    5   continue
      else
        lmaxl = lmaxp
        lmaxx = lmaxp
        r1 = 0
        nr = 1
        nr1 = 1
        nr2 = 1
      endif
      nlml = (lmaxl+1)**2
      nlmx = (lmaxx+1)**2
      call defrr(orl1, nr1*nlml*nsp)
      call defrr(orl2, nr2*nlml*nsp)
      call defrr(ovl,  nr1*nlml)
      call defrr(oxi,  (lmaxl+1)*nr)
      i0 = ioffv0(ib)+1
      i1 = ioff(ib)+1
      if (ipr >= 30) write(stdo,345) ib,nlml,r1,rsm(is)
  345 format(/' vtailz,  atom',i3,' :    nlml=',i3,'    rmt=',f10.6,
     .   '    rsm=',f10.6)
C --- Add together tail density and estatic potential ---
      call defrr(orlx,   nr*nlml*nsp)
      call sumtlz(rhoi,poti,pot0,nxi,lxi,exi,n0,rint,nbas,nbasp,lmaxp,
     .  alat,pos,ips,ioff,ioffv0,cg,jcg,indxcg,cy,ib,r1,nr,nr1,
     .  rofi,nlml,nsp,ixp,dxp,vval(i0),w(orlx),w(ovl),w(oxi),sxi0)
      vmtail = y0*vval(i0)
C     Skip interstitial integrals for classical multipoles
      if (ib > nbas) goto 50
      call xxcopi(nr,nr1,  0,nlml*nsp,w(orlx),w(orl1))
      call xxcopi(nr,nr2,nr1,nlml*nsp,w(orlx),w(orl2))
      call rlse(orlx)

C --- Save tail density then add smooth head density ---
      call defrr(orlt,    nr1*nlml)
      call defrr(odrl2,   nr2*nlml*nsp)
      if (nsp == 1) then
        call dpcopy(w(orl1),w(orlt),1,nr1*nlml,1d0)
      else
        call dpscop(w(orl1),w(orlt),nr1*nlml,1+nr1*nlml,1,1d0)
        call dpadd(w(orlt),w(orl1),1,nr1*nlml,1d0)
      endif
      call smrohx(nr,nr1,nr2,rofi,rwgt,nxi(is),lxi(1,is),exi(1,is),
     .   rsm(is),w(oxi),rhoi(i1),nri,nlml,nsp,w(orl1),w(orl2),w(odrl2))

C --- Yl-expansion of smooth xc potential ---
C     nx=lmaxl+4
      nx = 2*lmaxl+2
      nph = 0
      if (lmaxl <= 6) then
        nx = nxl(lmaxl)
        if (lxcfun/100 > 0) nx = nxl(lmaxl+1)
      endif
      call fpiint(nx,nph,np,p,wp)
      if (np > nnn) call rx('vtailz: increase nnn')
      if (ipr >= 30) write(stdo,200) nx,nph,np,nr1
  200 format(' mesh:    nx,nph=',2i4,'   gives',i4,'  angular points,',
     .    '   nrad=',i5)
      nrwk = max0(nr1,nr2)
      call defrr(ovxc1,   -nlml*nr1*nsp)
      call defrr(ovxc2,   -nlml*nr2*nsp)
      call defrr(odvxc2,  -nlml*nr2*nsp)
      call defrr(oyl,     (lmaxl+2)**2*np)
      call defrr(orp,     np*nrwk*nsp)
      call defrr(oexc,    np*nrwk*nsp*2)
      call defrr(ovxc,    np*nrwk*(1+nsp))
      call defrr(odrp,    np*nr2*nsp)
C      call setylm(nlml,np,p,cy,w(oyl))
C      call xcpmsh(nlml,np,nr1,rofi,rwgt,w(oyl),wp,w(orl1),
C     .  w(orp),w(oexc),w(ovxc),qsm,rep,rmu,w(ovxc1))
C      call xcpmsx(nlml,np,nr2,rofi(nr1+1),rwgt(nr1+1),w(oyl),wp,
C     .   w(orl2),w(orp),w(oexc),w(ovxc),w(ovxc2),w(odvxc2),w(odrl2),
C     .   w(odrp),q2,rep2,rmu2)
      call defrr(or2,     np)
      call defrr(oxp,     np)
      call defrr(oyp,     np)
      call defrr(ozp,     np)
      call dcopy(np,p(1,1),3,w(oxp),1)
      call dcopy(np,p(2,1),3,w(oyp),1)
      call dcopy(np,p(3,1),3,w(ozp),1)
      call ropyln(np,w(oxp),w(oyp),w(ozp),lmaxl+1,np,w(oyl),w(or2))
      call defrr(oylwp,   nlml*np)
C ... Memory allocation for gradient corrections
      if (lxcg() > 0) then
        call defrr(ogyl,     nlml*np*3)
        call ropylg(1,lmaxl,nlml,np,np,w(oxp),w(oyp),w(ozp),w(or2),
     .    w(oyl),w(ogyl))
        call defrr(ogrp,     nrwk*np*nsp*3)
        call defrr(oagrp,    nrwk*np*nsp)
        call defrr(oggrp,    nrwk*np*nsp)
        call defrr(owk1,     nrwk*nsp)
        call defrr(owk2,     nrwk*nsp)
        call defrr(oagrl,    nrwk*nlml*nsp)
        call defrr(ovxcl,    nrwk*nlml*nsp)
        call defrr(oexcl,    nrwk*nlml*nsp)
      endif
C ... vxc2 and dvxc2 beyond r1
      call pshpr(ipr-20)
      call dcopy(nlml*np,w(oyl),1,w(oylwp),1)
      call vxcns2(10,rofi(nr1+1),nr2,rwgt(nr1+1),np,wp,w(oylwp),nlml,
     .  nsp,w(orl2),w(1),lxcfun,.false.,w(1),w(1),w(1),w(1),rep,
     .  repsx,repsc,rmu,w(ovxc2),w(1),qsm)
      if (lxcg() > 0) then
        call vxcnls(rofi(nr1+1),1,nr2,np,nlml,nsp,w(orp),
     .    w(ogrp),w(oggrp),w(oagrp),w(oyl),w(ogyl),w(oylwp),rwgt(nr1+1),
     .    wp,w(orl2),w(oagrl),lxcg(),w(ovxcl),w(oexcl),w(ovxc2),rep,rmu)
      endif
      call dpadd(w(odrl2),w(orl2),1,nr2*nlml*nsp,1d0)
      call dcopy(nlml*np,w(oyl),1,w(oylwp),1)
      call vxcns2(10,rofi(nr1+1),nr2,rwgt(nr1+1),np,wp,w(oylwp),nlml,
     .  nsp,w(odrl2),w(1),lxcfun,.false.,w(1),w(1),w(1),w(1),rep2,
     .  repsx,repsc,rmu2,w(odvxc2),w(1),q2)
      if (lxcg() > 0) then
        call vxcnls(rofi(nr1+1),1,nr2,np,nlml,nsp,w(orp),
     .    w(ogrp),w(oggrp),w(oagrp),w(oyl),w(ogyl),w(oylwp),rwgt(nr1+1),
     .    wp,w(odrl2),w(oagrl),lxcg(),w(ovxcl),w(oexcl),w(odvxc2),rep2,
     .    rmu2)
      endif
      call dpadd(w(odrl2),w(orl2),1,nr2*nlml*nsp,-1d0)
      call dpadd(w(odvxc2),w(ovxc2),1,nr2*nlml*nsp,-1d0)
      call poppr
      do  28  i = 1, nsp
      q2(i) = q2(i)-qsm(i)
      rep2(i) = rep2(i)-rep(i)
   28 rmu2(i) = rmu2(i)-rmu(i)
C ... vxc inside r1
      call dcopy(nlml*np,w(oyl),1,w(oylwp),1)
      call vxcns2(10,rofi,nr1,rwgt,np,wp,w(oylwp),nlml,
     .  nsp,w(orl1),w(1),lxcfun,.false.,w(1),w(1),w(1),w(1),rep,
     .  repsx,repsc,rmu,w(ovxc1),w(1),qsm)
      if (lxcg() > 0) then
        b1 = r1/(dexp(a1*nr1-a1)-1)
        call vxcnls(rofi,0,nr1,np,nlml,nsp,w(orp),w(ogrp),
     .    w(oggrp),w(oagrp),w(oyl),w(ogyl),w(oylwp),rwgt,wp,w(orl1),
     .    w(oagrl),lxcg(),w(ovxcl),w(oexcl),w(ovxc1),rep,rmu)
      endif
      call rlse(oyl)

C      call prrmsh('vxc2 in vtailz',rofi(nr1+1),w(ovxc2),
C     .  nr2*nlml,nr2,nsp)
C      call prrmsh('dvxc2 in vtailz',rofi(nr1+1),w(odvxc2),
C     .  nr2*nlml,nr2,nsp)

C --- Subtract head contribution from zetxc ---
      call defrr(oxi0,   nr*(lmaxl+1))
      call xczetx(nr,nr1,nr2,rofi,rwgt,nxi(is),lxi(1,is),exi(1,is),
     .   rsm(is),w(oxi),w(oxi0),nlml,w(ovxc1),w(ovxc2),w(odvxc2),
     .   zetxc(i1),nri)
       call rlse(oxi0)

C --- Add on-site part of vval and sxi0 ---
      ix = i1
      lmax = -1
      do  29  ie = 1, nxi(is)
      call vtlxx1(nlml,poti(ix),lxi(ie,is),exi(ie,is),r1,vval(i0))
      call vtlxx2(exi(ie,is),r1,sxi0(ix))
      lmax = max0(lmax,lxi(ie,is))
   29 ix = ix+(lxi(ie,is)+1)**2
      call vtlxx1(nlml,pot0(i0),lmax,0d0,r1,vval(i0))

C ... Re-entry point for classical multipole sites (ib>nbas)
   50 continue
      vmad = y0*vval(i0)
      if (ipr >= 30) write(stdo,333) vmad-vmtail,vmtail,vmad
  333 format(' vmad...   head=',f11.6,'    tails=',f11.6,
     .   '    tot=',f11.6)

C --- Subtract sphere integrals from zeta and zetxc ---
C      call subtlz(nxi,lxi,exi,n0,rint,nbas,alat,pos,ips,ioff,
C     .   cg,jcg,indxcg,cy,pos(1,ib),r1,nr1,nr2,rofi,rwgt,nlml,nlmx,
C     .   ixp,dxp,
C     .   w(ovl),w(ovxc1),w(orlt),w(odvxc2),zeta,zetxc,ib,rhoi,
C     .   poti,pot0,ioffv0,vval(i0),qmom(i0),f1,f2,f4)
      call subtls(nxi,lxi,exi,n0,rint,nbas,nbasp,lmaxp,alat,pos,ips,
     .  ioff,cg,jcg,indxcg,cy,r1,nr1,nr2,rofi,rwgt,nlml,nlmx,ixp,dxp,
     .  w(ovl),w(ovxc1),w(orlt),w(odvxc2),zeta,zetxc,ib,rhoi,w(orhoi),
     .  poti,pot0,ioffv0,vval(i0),qmom(i0),upm,f1,f2,f4)

      call rlse(orl1)
      if (ipr >= 30 .and. nsp == 2) then
        write(stdo,967) 1,q2(1),rep2(1),rmu2(1)
        write(stdo,967) 2,q2(2),rep2(2),rmu2(2)
  967   format(8x,'spin ',i1,'   q2=',f10.6,'    rep2=',f10.6,
     .  '   rmu2=',f10.6)
      endif
      if (ipr >= 30) write(stdo,968) r2,q2(1)+q2(2),rep2(1)+rep2(2),
     .  rmu2(1)+rmu2(2)
  968 format(' rcut=',f8.3,'   q2=',f10.6,'    rep2=',f10.6,
     .  '   rmu2=',f10.6)
      qsmo = qsmo + qsm(1)+qsm(2)-q2(1)-q2(2)
      asmo = asmo + qsm(2)-qsm(1)-q2(2)-q2(1)
      repsmo = repsmo + rep(1)+rep(2)-rep2(1)-rep2(2)
      rmusmo = rmusmo + rmu(1)+rmu(2)-rmu2(1)-rmu2(2)
   10 continue
      call rlse(orhoi)
      if (MPI) then
        call mpibc2(vval,nvi,4,3,mlog,'vtailz','vval')
        call mpibc2(zeta,nri,4,3,mlog,'vtailz','zeta')
        call mpibc2(zetxc,nri*nsp,4,3,mlog,'vtailz','zetxc')
        call mpibc2(sxi0,nri,4,3,mlog,'vtailz','sxi0')
        call mpibc2(qsmo,1,4,3,mlog,'vtailz','qsmo')
        call mpibc2(asmo,1,4,3,mlog,'vtailz','asmo')
        call mpibc2(repsmo,1,4,3,mlog,'vtailz','repsmo')
        call mpibc2(rmusmo,1,4,3,mlog,'vtailz','rmusmo')
        call mpibc2(upm,1,4,3,mlog,'vtailz','upm')
        call mpibc2(f1,3*nbasp,4,3,mlog,'vtailz','f1')
        call mpibc2(f2,3*nbasp,4,3,mlog,'vtailz','f2')
        call mpibc2(f4,3*nbasp,4,3,mlog,'vtailz','f4')
C --- restore buffered components of zetxc and f1 --
        ierr = mpipid(2)
        call daxpy(nri*nsp,1d0,ztxcbk,1,zetxc,1)
        call daxpy(3*nbasp,1d0,f1bk,1,f1,1)
        deallocate (ztxcbk, stat=ierr)
        deallocate (f1bk, stat=ierr)
        deallocate (bproc, stat=ierr)
      endif

C --- Printout ---
      if (ipr >= 40) then
        write(stdo,289)
        do   ib = 1, nbas
          write(stdo,288) ib,(f1(m,ib),m=1,3)
        enddo
        write(stdo,279)
        do   ib = 1, nbas
          write(stdo,288) ib,(f2(m,ib),m=1,3)
        enddo
        write(stdo,275)
        do  ib = 1, nbas
          write(stdo,288) ib,(f4(m,ib),m=1,3)
        enddo
      endif
      call tcx('vtailz')
  288 format(i6,3f13.6)
  289 format(/'    ib     xc-force from interstitial density')
  279 format(/'    ib     es-force from interstitial density')
  275 format(/'    ib     f4')
      end

      subroutine vtlxx1(nlml,poti,lmax,e,r1,vval)
C- Add onsite terms to vval
      implicit none
      integer nlml,lmax
      double precision vval(nlml),poti(1),e,r1
      integer ilm,l,ll,lmaxl,lmux,nlmx
      double precision phi(0:20),psi(0:20)
      lmaxl=ll(nlml)
      lmux=min0(lmaxl,lmax)
      nlmx=(lmux+1)**2
      call bessl(e*r1*r1,lmux,phi,psi)
      do 29 ilm=1,nlmx
      l=ll(ilm)
   29 vval(ilm)=vval(ilm)+poti(ilm)*psi(l)/r1**(l+1)
      end
      subroutine vtlxx2(e,r1,sxi0)
C- Add onsite terms to sxi0
      implicit none
      double precision e,r1,sxi0(1)
      double precision phi(0:2),psi(0:2),srfpi,sumh
      srfpi=dsqrt(16d0*datan(1d0))
      call bessl(e*r1*r1,1,phi,psi)
      sumh=srfpi*psi(1)/e
      sxi0(1)=sxi0(1)-sumh
      end
      subroutine xxcopi(nr,nr1,nradd,nlm,rl,rl1)
C- copy out part of mesh
      implicit none
      integer nlm,nr,nr1,nradd
      double precision rl(nr,nlm),rl1(nr1,nlm)
      integer ilm,ir
      do 10 ilm=1,nlm
      do 10 ir=1,nr1
  10  rl1(ir,ilm)=rl(ir+nradd,ilm)
      end
