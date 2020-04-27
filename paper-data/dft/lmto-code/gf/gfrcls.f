C#define FAST
      subroutine gfrcls(mode,ic1,ic2,avw,alat,plat,bas,zkin,cy,cg,
     .  indxcg,jcg,rmax,ips,iaxb,clp,clssl,ntabb,offgH,cllst,offcH,
     .  ndofH,offH,pfun,pfj,dpfj,ddpfj,g)
C- Real space, free-electron Green's function for a group of clusters
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          1 return with real and imaginary separated
Ci         10s digit
Ci          1 make gfree, perturbed by pfun at each site
Ci   ic1   :starting cluster for which to make g
Ci   ic2   :last cluster for which to make g
Ci   avw   :a length scale, usu. average Wigner-Seitz sphere radius
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   zkin  :complex kinetic energy
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Ci   rmax  :augmentation radius, in a.u.
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   iaxb  :neighbor table for a cluster (paircl.f)
Ci   dpfj  :scaling to convert gfree to some screened representation.
Ci          usu. sqrt(P-dot) or ratio of sqrt(P-dot), viz:
Ci          to Gamma, make dpfj by calling mkpotj with iopt=130
Ci          to alpha, make dpfj by calling mkpotj with iopt=40
Ci   ddpfj :scaling to convert diagonal gfree to a screened repsn.
Ci          usu sqrt P-dot-dot, bare rep a difference in two, viz:
Ci          to Gamma, make dpfj by calling mkpotj with iopt=120
Ci          to alpha, make dpfj by calling mkpotj with iopt=50
Co Outputs
Co   g     :Free-electron GF
Cr Remarks
C ----------------------------------------------------------------
      implicit none
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer ic1,ic2,mode,niax,nclp,ndofH,ips(*),offH(n0H,nkap0,1)
      integer ntabb(ic2),offgH(1),cllst(*)
      parameter (niax=10,nclp=9)
      integer iaxb(niax,1),clp(nclp,ic2),offcH(ndofH,1),clssl(*)
      double precision avw,alat,plat(3,3),bas(3,ic2),rmax(*),
     .  pfun(2,1),pfj(2,1),dpfj(2,1),ddpfj(2,1),zkin(2),g(*)
C For Clebsch Gordan coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local variables
      integer ic,ib,ntvec,offxb,ndimg,offg,nsc,owk,ndimb,siz,offb
      integer nlmb,kb,mode0,mode1,stdo,lgunit

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
C      print *, 'setting mode1'
C      mode1 = 1
      stdo = lgunit(1)
      call awrit1(' gfrcls: mode=%i',' ',80,stdo,mode)

      do  10  ic = ic1, ic2

C       Number of sites in the cluster
        nsc  = clp(1,ic+1) - clp(1,ic)
        do  20  kb = clp(1,ic)+1, clp(1,ic+1)
          ib = cllst(kb)

          call gfornf(clp,clssl,ntabb,offgH,ib,ic,offxb,ndimg,siz,offg,
     .      offb)

          ntvec = ntabb(ib+1) - offxb
C         Column dimension for entire cluster
          ndimb = clp(3,ic)
C         Column dimension for this site
          nlmb = iaxb(9,1+offxb)

C     ... Make perturbed gfree for some pairs in this cluster
          if (mode1 == 1) then
            call pgfpcl(ib,ndimg,ndimb,nlmb,offb,ntvec,iaxb(1,1+offxb),
     .        avw,alat,plat,bas,ips,rmax,cy,cg,indxcg,jcg,zkin,
     .        offcH(1,ic),offH,pfun,pfj,dpfj,ddpfj,g(1+offg))
            mode0 = 0

C     ... Make gfree for all pairs in this cluster
          else
            call pgfrcl(ib,ndimg,ndimb,nlmb,offb,ntvec,iaxb(1,1+offxb),
     .        avw,alat,plat,bas,ips,rmax,cy,cg,indxcg,jcg,zkin,
     .        offcH(1,ic),offH,dpfj,ddpfj,g(1+offg))
          endif
   20   continue

        if (mode0 == 1) then
          call ztoy(g(1+offg),ndimg*ndimb,ndimg*ndimb,1,0)
C         call yprm('gfree',2,g(1+offg),ndimg*ndimb,ndimg,ndimg,ndimb)
        else
C         call yprm('gfree',3,g(1+offg),ndimg*ndimb,ndimg,ndimg,ndimb)
        endif
   10 continue

      end

      subroutine pgfrcl(ib,ndimg,ndimb,nlmb,offb,nvec,iaxb,avw,alat,
     .  plat,bas,ips,rmax,cy,cg,indxcg,jcg,zkin,offcH,offH,dpfj,
     .  ddpfj,gf)
C- Free-electron Green's function for all pairs connected to one cluster
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimg :leading dimension of g
Ci   nvec  :number of connecting vectors to the cluster
Ci   iaxb  :neighbor table for clster (paircl.f)
Ci   avw   :a length scale, usu. average Wigner-Seitz sphere radius
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   rmax  :augmentation radius, in a.u.
Ci   cy    :normalization constants for spherical harmonics (sylmnc.f)
Ci   cg,indxcg,jcg:Clebsch Gordan coefficients, and indices (scg.f)
Ci   dpfj  :scaling to convert gfree to some screened representation.
Ci          usu. sqrt(P-dot) or ratio of sqrt(P-dot), viz:
Ci          to Gamma, make dpfj by calling mkpotj with iopt=130
Ci          to alpha, make dpfj by calling mkpotj with iopt=40
Ci   ddpfj :scaling to convert diagonal gfree to a screened repsn.
Ci          usu sqrt P-dot-dot, bare rep a difference in two, viz:
Ci          to Gamma, make ddpfj by calling mkpotj with iopt=120
Ci          to alpha, make ddpfj by calling mkpotj with iopt=50
Ci   z     :energy
Co Outputs
Co   gf    :free-electron Green's function for all field point
Co          sites connected by the neighbor table to the source point.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nkap0,n0H,niax
      parameter (nkap0=4,n0H=5,niax=10)
      integer ndimg,nvec,ips(*),offcH(nvec),offH(n0H,nkap0,1),offb,nlmb
      integer iaxb(niax,nvec)
      double precision alat,avw,bas(3,*),plat(3,3),gf(2,ndimg,*)
      double precision rmax(*)
      double complex zkin,dpfj(*),ddpfj(*)
C For Clebsch Gordan coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      integer ipr,offa,nlma,ia,ib,ip,offPb,offPa,ndimb
      double precision dsqr,drr2,dr(3),rmaxa,rmaxb

C ... Setup and memory allocation
      call getpr(ipr)
      call tcn('pgfrcl')

      offPb = offH(1,1,ib)

      if (ips(1) > 0) then
        rmaxb = rmax(ips(ib))
      else
        rmaxb = rmax(ib)
      endif

C --- For each connecting vector, make gf^bare ---
      do  20  ip = 1, nvec
        ia = iaxb(2,ip)
        offa = offcH(ip)
        nlma = offcH(ip+1) - offa
        offPa = offH(1,1,ia)
        if (ips(1) > 0) then
          rmaxa = rmax(ips(ia))
        else
          rmaxa = rmax(ia)
        endif

C   ... Connecting vector
        dsqr = drr2(plat,bas(1,ia),bas(1,ib),iaxb(3,1)-iaxb(3,ip),
     .    iaxb(4,1)-iaxb(4,ip),iaxb(5,1)-iaxb(5,ip),dr)
        dr(1) = dr(1)*alat/avw
        dr(2) = dr(2)*alat/avw
        dr(3) = dr(3)*alat/avw

C   ... Free-electron g^alpha for this connecting vector
        call gfree(0,nlma,nlmb,ndimg,avw,dr,rmaxa,rmaxb,zkin,indxcg,
     .    jcg,cg,cy,dpfj(1+offPa),dpfj(1+offPb),ddpfj(1+offPa),
     .    gf(1,1+offa,1+offb),gf(1,1+offa,1+offb))

C   ... debugging printout
C       print 368, ip,offa,offb,ia,ib,iaxb(3,1)-iaxb(3,ip),
C     .   iaxb(4,1)-iaxb(4,ip),iaxb(5,1)-iaxb(5,ip)
C  368  format(3i3,4x,10i3)
C       print 331, ia,ib,dr(1)*avw/alat,dr(2)*avw/alat,dr(3)*avw/alat
C  331  format(' sites',2i4,' vector:',3f12.6)
C       call yprm('gfree',3,gf(1,1+offa,1+offb),
C     .    ndimg*ndimb,ndimg,nlma,nlmb)

   20 continue

C ... debugging check
C#ifdefC DEBUG
C      if (offa+nlma > ndimg) call rx('bad ndimg in pgfrcl')
C      call yprm('gfree',3,gf(1,1,1+offb),ndimg*ndimb,ndimg,ndimg,nlmb)
C#endif
      call tcn('pgfrcl')
      end

      subroutine pgfpcl(ib,ndimg,ndimb,nlmb,offb,nvec,iaxb,avw,alat,
     .  plat,bas,ips,rmax,cy,cg,indxcg,jcg,zkin,offcH,offH,pfun,pfj,
     .  dpfj,ddpfj,gf)
C- Perturbed free-electron Green's function for selected pairs in a cluster
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimg :leading dimension of g
Ci   nvec  :number of connecting vectors to the cluster
Ci   iaxb  :neighbor table for clster (paircl.f)
Ci   avw   :a length scale, usu. average Wigner-Seitz sphere radius
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   rmax  :augmentation radius, in a.u.
Ci   cy    :normalization constants for spherical harmonics (sylmnc.f)
Ci   cg,indxcg,jcg:Clebsch Gordan coefficients, and indices (scg.f)
Ci   dpfj  :scaling to convert gfree to some screened representation.
Ci          usu. sqrt(P-dot) or ratio of sqrt(P-dot), viz:
Ci          to Gamma, make dpfj by calling mkpotj with iopt=130
Ci          to alpha, make dpfj by calling mkpotj with iopt=40
Ci   ddpfj :scaling to convert diagonal gfree to a screened repsn.
Ci          usu sqrt P-dot-dot, bare rep a difference in two, viz:
Ci          to Gamma, make ddpfj by calling mkpotj with iopt=120
Ci          to alpha, make ddpfj by calling mkpotj with iopt=50
Ci   z     :energy
Co Outputs
Co   gf    :free-electron Green's function for all field point
Co          sites connected by the neighbor table to the source point.
Cr Remarks
Cu Updates
Cu   Brute force inversion can be more efficient.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer ndimg,nvec,niax,ips(*),ndimb,offb,nlmb,nkap0,n0H
      parameter (nkap0=4,n0H=5,niax=10)
      integer iaxb(niax,nvec),offcH(nvec),offH(n0H,nkap0,1)
      double precision alat,avw,bas(3,*),plat(3,3),gf(ndimg,ndimb,2)
      double precision rmax(*),pfun(2,*),pfj(2,*)
      double complex zkin,dpfj(*),ddpfj(*)
C For Clebsch Gordan coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      integer ipr,offa,nlma,ia,ib,ip,offPb,offPa,ndimx,np1,
     .  ll,lmxa,lmxb,ilma,ilmb,la,lb,ma,mb,nlmab,ierr,klma
      parameter (ndimx=81*2,np1=ndimx+1)
      double precision gc16(2,ndimx,ndimx+1),
     .  gab(ndimx,ndimx,2),gdP(ndimx,ndimx,2),gwk(2,ndimx,2)
      double precision dsqr,drr2,dr(3),rmaxa,rmaxb,xx,dPr,dPi,dpgr,dpgi
      double complex x12

C ... Setup and memory allocation
      call getpr(ipr)
      call tcn('pgfpcl')

      offPb = offH(1,1,ib)
      if (ips(1) > 0) then
        rmaxb = rmax(ips(ib))
      else
        rmaxb = rmax(ib)
      endif

C --- On-site free-electron Green's function, source point ---
      call dpzero(dr,3)
      call gfree(0,nlmb,nlmb,ndimx,avw,dr,rmaxb,rmaxb,zkin,indxcg,
     .  jcg,cg,cy,dpfj(1+offPb),dpfj(1+offPb),ddpfj(1+offPb),gc16,gc16)
      call zcopy(nlmb,gc16,np1,gwk(1,1,1),1)
C      call dcopy(nlmb,gc16,2*np1,gwk(1,1,1),2)
C      call dcopy(nlmb,gc16(2,1,1),2*np1,gwk(2,1,1),2)

      lmxb = ll(nlmb)

C --- For each connecting vector, make gf^bare ---
      do  20  ip = 1, nvec
        ia = iaxb(2,ip)
        offa = offcH(ip)
        nlma = offcH(ip+1) - offa
        lmxa = ll(nlma)
        offPa = offH(1,1,ia)
        if (ips(1) > 0) then
          rmaxa = rmax(ips(ia))
        else
          rmaxa = rmax(ia)
        endif
        nlmab = nlma + nlmb

        call dpzero(gab,ndimx*nlmab)
        call dpzero(gab(1,1,2),ndimx*nlmab)

C ---   On-site free-electron Green's function, field point ---
        call dpzero(dr,3)
        call gfree(0,nlma,nlma,ndimx,avw,dr,rmaxa,rmaxa,zkin,indxcg,
     .   jcg,cg,cy,dpfj(1+offPa),dpfj(1+offPa),ddpfj(1+offPa),gc16,gc16)
        call dcopy(nlmb,gc16,2*np1,gab(1+nlmb,1+nlmb,1),np1)
        call dcopy(nlmb,gc16(2,1,1),2*np1,gab(1+nlmb,1+nlmb,2),np1)
        call dcopy(nlmb,gwk,2,gab(1,1,1),np1)
        call dcopy(nlmb,gwk(2,1,1),2,gab(1,1,2),np1)

C   ... Free-electron g^alpha_ab
        dsqr = drr2(plat,bas(1,ia),bas(1,ib),iaxb(3,1)-iaxb(3,ip),
     .    iaxb(4,1)-iaxb(4,ip),iaxb(5,1)-iaxb(5,ip),dr)
        dr(1) = dr(1)*alat/avw
        dr(2) = dr(2)*alat/avw
        dr(3) = dr(3)*alat/avw
        if (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) < 1d-8) then
          nlmab = nlmb
          goto 28
        endif
        call gfree(0,nlma,nlmb,ndimx,avw,dr,rmaxa,rmaxb,zkin,indxcg,
     .   jcg,cg,cy,dpfj(1+offPa),dpfj(1+offPb),ddpfj(1+offPa),gc16,gc16)

C   ... Copy g_ab and g_ba
        ilmb = 0
        do  22  lb = 0, lmxa
        do  22  mb = -lb, lb
          ilmb = ilmb+1
          ilma = 0
          xx = - (-1)**lb
          do  23  la = 0, lmxa
          xx = -xx
          do  23  ma = -la, la
            ilma = ilma+1
            gab(ilma+nlmb,ilmb,1) = gc16(1,ilma,ilmb)
            gab(ilma+nlmb,ilmb,2) = gc16(2,ilma,ilmb)
            gab(ilma,ilmb+nlmb,1) = xx*gc16(1,ilma,ilmb)
            gab(ilma,ilmb+nlmb,2) = xx*gc16(2,ilma,ilmb)
   23     continue
   22   continue

C   ... Make dPb for b column
   28   continue
        offPb = offH(1,1,ib)
        offPa = offH(1,1,ia) - nlmb
        do  24  ilmb = 1, nlmab
          if (ilmb <= nlmb) then
            dPr = pfun(1,offPb+ilmb) - pfj(1,offPb+ilmb)
            dPi = pfun(2,offPb+ilmb) - pfj(2,offPb+ilmb)
          else
            dPr = pfun(1,offPa+ilmb) - pfj(1,offPa+ilmb)
            dPi = pfun(2,offPa+ilmb) - pfj(2,offPa+ilmb)
          endif

C          dpr = 0
C          dpi = 0
          do  25  ilma = 1, nlmab
            dPgr = gab(ilma,ilmb,1)*dPr - gab(ilma,ilmb,2)*dPi
            dPgi = gab(ilma,ilmb,2)*dPr + gab(ilma,ilmb,1)*dPi
            gdP(ilma,ilmb,1) = dPgr
            gdP(ilma,ilmb,2) = dPgi
   25     continue
            gdP(ilmb,ilmb,1) = gdP(ilmb,ilmb,1) + 1
   24   continue

C         print *, 'dr=',ip,dr
C         call yprm('gab',2,gab,ndimx*ndimx,ndimx,nlmab,nlmab)
C         call yprm('1+gdP',2,gdp,ndimx*ndimx,ndimx,nlmab,nlmab)

C#ifdef FAST
C   ... Fast inversion, assuming gab is diagonal in the (11,22) blocks
        if (nlmab /= nlmb) then
          do  30  klma = 1, nlma
          do  30  ilmb = 1, nlmb
            x12 = dcmplx(gdP(ilmb,nlmb+klma,1),gdP(ilmb,nlmb+klma,2))
     .          / dcmplx(gdP(ilmb,ilmb,1),gdP(ilmb,ilmb,2))
            dPr = dble(x12)
            dPi = dimag(x12)
            do  32  ilma = 1, nlma
              gdP(nlmb+ilma,nlmb+klma,1) = gdP(nlmb+ilma,nlmb+klma,1) -
     .        (gdP(nlmb+ilma,ilmb,1)*dPr - gdP(nlmb+ilma,ilmb,2)*dPi)
              gdP(nlmb+ilma,nlmb+klma,2) = gdP(nlmb+ilma,nlmb+klma,2) -
     .        (gdP(nlmb+ilma,ilmb,2)*dPr + gdP(nlmb+ilma,ilmb,1)*dPi)
   32       continue
   30     continue
C         call yprm('a22',2,gdp,ndimx*ndimx,ndimx,nlmab,nlmab)

C         Set up rhs
          do  36  ilmb = 1, nlmb
            x12 = dcmplx(gab(ilmb,ilmb,1),gab(ilmb,ilmb,2))
     .          / dcmplx(gdP(ilmb,ilmb,1),gdP(ilmb,ilmb,2))
          do  36  ilma = 1, nlma
            dPr = gdP(nlmb+ilma,ilmb,1)*dble(x12) -
     .            gdP(nlmb+ilma,ilmb,2)*dimag(x12)
            dPi = gdP(nlmb+ilma,ilmb,2)*dble(x12) +
     .            gdP(nlmb+ilma,ilmb,1)*dimag(x12)
            gab(ilma+nlmb,ilmb,1) = gab(ilma+nlmb,ilmb,1) - dPr
            gab(ilma+nlmb,ilmb,2) = gab(ilma+nlmb,ilmb,2) - dPi
   36     continue
        endif

C       Offset to gab(2,1) is nlmb, except for dr=0.  Offset in nlmab
        nlmab = nlmab - nlma
        call yygefa(gdp(nlmab+1,nlmab+1,1),gdp(nlmab+1,nlmab+1,2),
     .    ndimx,nlma,gc16,ierr)
        if (ierr /= 0) call rx(' failed to invert g')
        do  38  ilmb = 1, nlmb
          call yygesl(gdp(nlmab+1,nlmab+1,1),gdp(nlmab+1,nlmab+1,2),
     .      ndimx,nlma,gc16,gab(1+nlmab,ilmb,1),gab(1+nlmab,ilmb,2),0)
   38   continue
C#elseC
C        call yyqnvb(' ',gdP,gdP(1,1,2),ndimx,nlmab,nlmab,gc16,nlmab,
C     .    gc16,gab,gab(1,1,2),ndimx,ierr)
C        if (ierr /= 0) call rx(' failed to invert g')
CC       Offset to gab(2,1) is nlmb, except for dr=0.  Offset in nlmab:
C        nlmab = nlmab - nlma
C#endif

C        call yprm('(1+gdP)^-1 g',2,gab(1+nlmab,1,1),ndimx*ndimx,ndimx,
C     .    nlma,nlma)

        call ymscop(0,nlma,nlmb,ndimx,ndimg,0,0,0,0,
     .    gab(1+nlmab,1,1),ndimx*ndimx,gf(1+offa,1+offb,1),ndimg*ndimb)
     .

C   ... debugging printout
C       print 368, ip,offa,offb,ia,ib,iaxb(3,1)-iaxb(3,ip),
C     .   iaxb(4,1)-iaxb(4,ip),iaxb(5,1)-iaxb(5,ip)
C  368  format(3i3,4x,10i3)
C       print 331, ia,ib,dr(1)*avw/alat,dr(2)*avw/alat,dr(3)*avw/alat
C  331  format(' sites',2i4,' vector:',3f12.6)
C        call yprm('gfree',3,gf(1,1+offa,1+offb),ndimg*ndimb,ndimg,nlma,
C     .    nlmb)

   20 continue

C ... debugging check
C#ifdefC DEBUG
C      if (offa+nlma > ndimg) call rx('bad ndimg in pgfrcl')
C      call yprm('gfree',2,gf(1,1+offb,1),ndimg*ndimb,ndimg,ndimg,nlmb)
C#endif
      call tcx('pgfpcl')
      end

