      subroutine mcatm(spid,nbas,nspec,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,nxi,lxi,exi,orhoi,orhoic,orho0,nri,nvi,ioff,ioffv0,
     .  z,rsmfa,pin,qin,irchan,rc,pnu,orho,orhoc,ov0,setot)
C- Make starting density as a Mattheis construction
      implicit none
      integer n0,nbas,nspec,npan,nsp,orhoi,orhoic,orho0,
     .  nri,nvi,ips(1),ioff(1),ioffv0(1)
C ... Dimensions for species (nsx)
      character*8 spid(1)
      integer nr(1),irchan(n0,1),lmxa(1),nxi(1),lxi(n0,1),lmxl(1)
      double precision z(nspec),rmt(nspec),rsmfa(nspec),a(nspec),
     .  exi(n0,nspec),pin(n0,2,nspec),qin(n0,2,nspec),
     .  rc(2,nspec),setot
C ... Dimensions for atoms (nbx)
      integer orho(1),orhoc(1),ov0(1)
      double precision pnu(n0,1),pz(20)
C ... Local variables
      integer orofi,ov,orhoin,oidmod,is,nrmx,nrmix(2),kcor,lcor,lwf,
     .  nrrmx,orhos,orhocs,orhoa,nrlmsp,nrsp,opz0,lxcfun,lxcf,idmod(20),
     .  nglob,nrmt,orhocl,ov00
      double precision v0,qc,qv,sec,stc,etot,y0,
     .  rsm,rtab(n0,2),etab(n0,2),rfoca,ccof,ceh,eref,qcor(2),
     .  hfc(n0,2),hfct(n0,2),rs3,eh3,vmtz,rcfa(2),
     .  basopt(15),pnuopt(n0),pzopt(n0),gtop(2)

      real w(1)
      common /w/ w

      call dcopy(20,0d0,0,pz,1)
      call icopy(20,1,0,idmod,1)
      call dcopy(2*n0,0d0,0,hfct,1)
C     for now, no renormalization of free atom
      call dcopy(2,0d0,0,rcfa,1)

C --- Setup ---
C     call wkprnt(1)
      v0 = 0
      setot = 0
      y0 = 1/dsqrt(16*datan(1d0))
      call defrr(orhoi,     nri*nsp)
      call defrr(orho0,     nvi*nsp)
      call dpzero(w(orhoi), nri*nsp)
      call dpzero(w(orho0), nvi*nsp)

      do  10  is = 1, nspec
C ---   Make free atom rho for this species, attach Hankel tails ---
        nrmx = 1501
C ...   work space for orhocs, ov, and orhoa will are retained
        call defrr(orhocs,     nrmx*2)
        call defrr(ov,         nrmx*2)
        call defrr(orhoa,     -nrmx*2)
        call defrr(orofi,      nrmx*2)
        call defrr(orhoin,     nrmx*2)

        lxcfun = mod(nglob('lxcf'),100)
        rfoca = 0.4d0*rmt(is)
C --- No smoothing for now ..
        rsm = 0d0
C        rsm = rsmfa(is)
        nrmix(1) = 80
        nrmix(2) = 2
        kcor = 0
        lcor = -1
        qcor(1) = 0d0
        qcor(2) = 0d0
        lwf = 1
        nrmt = nr(is)
        nrrmx = nr(is)
        eref = 0
        gtop(1) = 0
        gtop(2) = 0
        call dpzero(basopt,10)

        call freats(spid(is),is,n0,nxi(is),exi(1,is),rfoca,rsm,
     .    kcor,lcor,qcor,0,nrmix,lwf,lxcfun,z(is),rmt(is),a(is),nrmt,
     .    pin(1,1,is),pz,qin(1,1,is),rs3,eh3,vmtz,rcfa,idmod,
     .    lmxa(is),eref,basopt,0,rtab,etab,hfc,hfct,nrrmx,w(orofi),
     .    w(orhoa),w(orhocs),qc,ccof,ceh,sec,stc,w(ov),etot,pnuopt,
     .    pzopt,gtop)
C   ... add core to valence
        call daxpy(nrrmx*nsp,1d0,w(orhocs),1,w(orhoa),1)
        nrsp   = nr(is)*nsp
        nrlmsp = nr(is)*nsp*(lmxl(is)+1)**2

C ---   Copy into arrays dimensioned (nr,..); scale by y0 ---
        call defrr(orhos,   -nrlmsp)
        call defrr(ov00,      nrsp)
        call defrr(orhocl,   nrsp)
        call dpcopy(w(orhoa),w(orhos),1,nr(is),y0)
        call dpcopy(w(ov),w(ov00),1,nr(is),1d0)
        call dpcopy(w(orhocs),w(orhocl),1,nr(is),1d0)
        if (nsp == 2) then
          call dpscop(w(orhoa),w(orhos),nr(is),1+nrrmx,1+nrlmsp/2,y0)
          call dpscop(w(ov),w(ov00),nr(is),1+nrrmx,1+nr(is),1d0)
          call dpscop(w(orhocs),w(orhocl),nr(is),1+nrrmx,1+nr(is),1d0)
        endif
        call mcatx1(is,n0,nbas,nsp,lmxa(is),ips,nrlmsp,nrsp,
     .    nri,nvi,nxi(is)+0,lxi(1,is),exi(1,is),ioff,ioffv0,
     .    etot,setot,hfc,orhos,orhocl,ov00,pin(1,1,is),orho,orhoc,ov0,
     .    pnu,w(orhoi),w(orho0))

   10 continue

C      call ppot3(nri,nxi,lxi,nvi,n0,nbas,ips,w(orhoi),w(orhoi),
C     .  w(orhoi),w(orho0),w(orho0),w(orho0),0)

      end
      subroutine mcatx1(is,n0,nbas,nsp,lmxa,ips,nrlmsp,nrsp,
     .  nri,nvi,nxi,lxi,exi,ioff,ioffv0,etot,setot,hfc,orhos,orhocs,ov,
     .  pin,orho,orhoc,ov0,pnu,rhoi,rho0)
C- Assemble contribution from one species free-atomic density
      implicit none
      integer is,nbas,ips(1),orhos,orhocs,ov,orho(1),orhoc(1),ov0(1),
     .  nrlmsp,nrsp,n0,lmxa,nri,nvi,nxi,lxi(1),ioff(1),ioffv0(1),nsp
      double precision pin(n0,2),pnu(n0,2,1),exi(1),etot,setot,
     .  rhoi(nri,1),rho0(nvi,1),hfc(n0,1)
      integer ib,off,off0,ie,lmax
      logical first
      real w(1)
      common /w/ w

      first = .true.
      do  10  ib = 1, nbas

C ...   snot for now to add a nonspherical component of the density
        if (is /= ips(ib)) goto 10
C ---   Copy rho for this species into appropriate basis ---
        setot = setot+etot
        if (first) then
          orho(ib) = orhos
          orhoc(ib) = orhocs
          ov0(ib)   = ov
          first = .false.
        else
          call defrr(orhoc(ib), nrsp)
          call defrr(ov0(ib),   nrsp)
          call defrr(orho(ib),  nrlmsp)
          call dpcopy(w(orhos),w(orho(ib)),1,nrlmsp,1d0)
          call dpcopy(w(orhocs),w(orhoc(ib)),1,nrsp,1d0)
          call dpcopy(w(ov),w(ov0(ib)),1,nrsp,1d0)
        endif
        call dpcopy(pin,pnu(1,1,ib),1,1+lmxa,1d0)
        if (nsp == 2) call dpcopy(pin(1,2),pnu(1,2,ib),1,1+lmxa,1d0)

C ---  Copy fit coffs into rhoi or rho0 ---
        off = ioff(ib)
        off0 = ioffv0(ib)
        lmax = -1
        do  20  ie = 1, nxi
          if (exi(ie) < -1d-9) then
            rhoi(1+off,1) = hfc(ie,1)
            rhoi(1+off,nsp) = hfc(ie,nsp)
            lmax = max(lmax,lxi(ie))
            off = off + (lxi(ie)+1)**2
          else
            rho0(1+off0,1) = hfc(ie,1)
            rho0(1+off0,nsp) = hfc(ie,nsp)
            off0 = off0 + (lmax+1)**2
          endif
   20   continue
   10 continue

      end
