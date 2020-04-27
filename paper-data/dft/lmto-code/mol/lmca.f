C#define MCA
C#ifdefC MCE
C      subroutine lmce(prgnam,s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
C     .  s_site,s_str,s_move,s_strn)
C#endif
C#ifdef MCA
      subroutine lmca(prgnam,s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,s_spec,
     .  s_site,s_str,s_move,s_strn)
C#endif
C#ifdefC MCFIT
C      subroutine lmcfit(prgnam,s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,
C     .  s_spec,s_site,s_str,s_move,s_strn)
C#endif
C#ifdefC MCXBS
C      subroutine lmxbs(prgnam,s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,
C     .  s_spec,s_site,s_str,s_move,s_strn)
C#endif
C#ifdefC MCRHO
C      subroutine lmrho(prgnam,s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,
C     .  s_spec,s_site,s_str,s_move,s_strn)
C#endif
C- MC Cluster / Molecule programs
C ----------------------------------------------------------------------
Cio Structures
Ci   ... lmce
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfrce quit mdprm nitmv lrel lxcf nbas nbasp nspec
Ci                 nspin zbak maxit tol defm ltb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lham lbas lcd lxcf
Cio    Passed to:  st2ary rlxstp relax
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  kmto nmto alfsi dabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z rmt a nr rg rint rcut rham name p q eref lxi exi
Ci                 colxbs radxbs lmxpb lmxa lmxl idmod mass rsmfa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary relax
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos vel mpole dpole spec relax clabel
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary rlxstp relax iopos
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  nbisi nalf ncupl ndust adec wztcf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  *
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  str_pack
Ci   ... lmca
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfrce mdprm nitmv lrel lxcf nbas nbasp nspec nspin
Ci                 zbak maxit tol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lham lbas lcd lxcf
Cio    Passed to:  st2ary
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  kmto nmto alfsi dabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z rmt a nr rg rint rcut rham name p q eref lxi exi
Ci                 colxbs radxbs lmxpb lmxa lmxl idmod mass rsmfa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos vel mpole dpole spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  nbisi nalf ncupl ndust adec wztcf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   ... lmcfit
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lfrce mdprm nitmv lrel lxcf nbas nbasp nspec nspin
Ci                 zbak maxit tol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lham lbas lcd lxcf
Cio    Passed to:  st2ary
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  kmto nmto alfsi dabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat nabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  w
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z rmt a nr rg rint rcut rham name p q eref lxi exi
Ci                 colxbs radxbs lmxpb lmxa lmxl idmod mass rsmfa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos vel mpole dpole spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  nbisi nalf ncupl ndust adec wztcf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  st2ary
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   prgnam:name of calling program
Co Outputs
Cs Command-line switches
Cs   --LUMO    :
Cs   --flush   :
Cs   --grfac=  :
Cs   --md=     :
Cs   --mlog    : (MPI) write MPI commands to log file
Cs   --mull    : Mulliken analysis; see doc/Command-line-options.html
Cs   --mv=     :
Cs   --plot    :
Cs   --psi     :
Cs   --rs      : Controls I/O with rst file; see Command-line-options.html
Cs   --st      :
Cs   --wforce= :
Cs   --xtoll=  :
Cs   --xyz=    :
Cs   -excite=  :
Cl Local variables
Cr Remarks
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   25 Oct 11 Started migration to f90 structures
Cu   14 Jan 05 Replaced call to mcsav with call to nwit
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters:
      character*(*)  prgnam
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
      type(str_mix)::   s_mix
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_str)::   s_str
      type(str_move)::  s_move
      type(str_strn) :: s_strn(*)
C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer n0,nsx,nbx,nvx
      parameter(n0=10, nsx=10, nbx=200, nvx=12)
C  dimensions for species (nsx)
      double precision z(nsx),rmt(nsx),rsm(nsx),rsmfa(nsx),rint(nsx),
     .  rcut(nsx),nphi(nsx),ephi(n0,nsx),a(nsx),
     .  exi(n0,nsx),pin(n0,nsx*2),qin(n0,nsx*2),rsm0(nsx),
     .  amass(nsx),piz(n0,nsx*2),colxbs(3,nsx),radxbs(nsx),rham(nsx),
     .  tphi(n0,5,nsx),alfa(2)
      integer lphi(n0,nsx),nr(nsx),idmod(n0,nsx),lmxa(nsx),lmxl(nsx),
     .  nxi(nsx),lxi(n0,nsx),idmoz(n0,nsx)
C  lmax of 'padded' point multipoles (lmaxp=1)
      integer lmaxp
C  dimensions for atoms (nbx)
      double precision  pos0(3,nbx),pol(3,nbx),vnuc(nbx),ipclas(nbx),
     .  pnu(n0,2,nbx),qus0(3,n0,2,nbx),qmull(nbx,2),qc(nbx),
     .  pos(3,nbx),mpole(nbx),dpole(3,nbx),psav(3,nbx),f(3,nbx),
     .  f1(3,nbx),f2(3,nbx),f3(3,nbx,2),f4(3,nbx),fsav(3,nbx),
     .  pnuz(n0,2,nbx),vel(3,nbx),dxp(10)
      integer ioff(nbx+1),ioffv0(nbx+1),ioffp(nbx+1),ips(nbx),ixp(10),
     .  ifrlx(3,nbx)
C  dimensions for Clebsch Gordans (lmax=6 .. 10))
      integer lmxcg,lnjcg,lnxcg
C      parameter(lmxcg=6, lnjcg=6409,  lnxcg=1226)
C      parameter(lmxcg=7, lnjcg=12537, lnxcg=2081)
C      parameter(lmxcg=8, lnjcg=22662, lnxcg=3322)
C      parameter(lmxcg=9, lnjcg=38492, lnxcg=5051)
       parameter(lmxcg=10,lnjcg=62153, lnxcg=7382)
C  other dimensions
      double precision cg(lnjcg),jcg(lnjcg),cy(300),el(10),hi(20),
     .  si(20),plat(3,3),amp(15),qrdiff(2),qrtol(2),apos(3),avel(3),
     .  qrdifm,qrtolm
      integer lsw,lswtch(20),indxcg(lnxcg),lrs(5),ndim(10),nstate(2)
      double precision dr(3),rchar(nsx),rcharw(nsx),
     .  tcpm(10,20),vals(15)
      integer nbisi(3),irchan(n0,nsx)

      logical cmdopt,a2bin,l,xyzfrz(3)
      character rhoid*64,outs*72,strn*128
      character*8 spid(nsx),atid(nbx)
      double precision dclabl(nsx)
      character symgrp*60
      character tabs*80
      character*3 ensbl

C Local variables
      integer nel,nfile,ntx,ndx,nft1,nft2,nft3,lsc,nsc,nbas,nbasp,nspec,
     .  ntcpm,ifi,nvario,nsp,npan,ipr,iprint,i,j,iamp,lmxst,nri,nvi,nhs,
     .  nla,nll,nelx,nhl,njj,nclas,nv,ng,nit,it,nhtab,mxcsiz,ib,ic,m,is,
     .  ipan,nbb1,ie,nlbx,isp,nev,i1mach,itdyn,icom,nvar,natrlx,itrlx,
     .  itmxrl,tskip,pdim,broy,lskip,i1,i2,mmix,str_pack
      double precision elink,beta1,beta2,del,dqval,etol,seref,
     .  alat,dx,dy,dz,as,ewtol,scale,qval,qsmc,vint0,etot0,rep0,
     .  rmu0,q0,amom0,rhovxc,qsmo,asmo,repsmo,rmusmo,qint,zvnuc,qsph,
     .  amsph,repsph,rmusph,rvhsph,sumec,sumtc,rvhi,rhoep,rhomu,rhovh,q,
     .  vint,vintf,upm,utot,ertot,sumev,wgt,sumev1,qsph1,qsum,asum,
     .  rescal,ehtot,ftop,xx,logs,zeta,zacc,zsqacc,ekin,tkin,cons,
     .  autime,etot(2),ft(3)
      integer nf,ierr
      double precision mdtime,temp,taup,step,tstep,veps,xtol,gtol,time,
     .                 time0,eps,dum
      double precision dmax1,dist0
      logical parmxp
C Switches
      integer lrx,ldyn,ltchk,lfrz,lvint,lxcqs,llink,lforce,lrdxc,lwrxc,
     .  lsp,lcor,ltail
C Pointers
      real(8), allocatable :: rtab(:,:)
C#ifndef LINUXF
      integer, allocatable, target :: iax(:)
C#elseC
C      integer, pointer :: iax(:)
C#endif
      integer orho(nbx),ov0(nbx),orhsms(nbx)
      integer orhoc(nbx)
      integer otspec,otdata,ochadd,obl,oblp,og,oag,orhoi,opoti,
     .  opot0,ozeta,ozetxc,ovval,osxi0,ohvecs,osvecs,orhoic,orho0,
     .  odummy,ontab,odvxci,orhoq,ovxc0,odvxc0,oqmom,
     .  oqmoi,ovab,osab,ohab,opuu,opus,opss,ozetf,orhsmi,oquu,oqus,oqss,
     .  orhout,ogh,ogj,ob,obp,oioffa,oioffb,oioffl,oh,os,osm,oc,osi,
     .  osub,ohvucs,osvucs,oevl,oewt,ot,owk,ov,oyk,op,ow,oppnl,oindrx

C#ifdefC MCFIT
C      integer nn,nph,js,ip,nalftc,ncupl,ndust
C      double precision ri1,ri2,adec,wftx,d0
C#endif
C#ifdef MCA
      double precision setot
C#endif
C#ifdefC RD
C      integer iov
C#endif

C ..for point multipoles
      data lmaxp /1/ ehtot /0d0/

C ..for Nose-Hoover and molecular statics
      data xyzfrz /.false.,.false.,.false./
      data itdyn /0/ zacc /0d0/ time /0d0/
      data autime /0.048377d0/

C For MPI ...
      integer mpipid,procid,master,numprocs
      logical MPI,mlog
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call tcn(prgnam)
C      call mopen(82,'pd','f')

      qval = 0d0

C  94= workfile to dump matrices, 73= eigenstates for later use
      if (procid == master) then
C        call mfdef(1)
        call mopen(94,'h1','u')
        call mopen(95,'h2','u')
C      call mopen(73,'st','u')
        call mopen(70,'mxi','u')
        call mopen(74,'mxs','u')
      endif
C      call poseof(82)
      call scg(lmxcg,cg,indxcg,jcg)
      call sylmnc(cy,6+6+1)

C --------- Input -----------------
      call st2ary(s_ctrl,s_ham,s_pot,s_lat,s_bz,s_mix,
     .  s_spec,s_site,s_str,s_move,
     .  nel,el,elink,
     .  alfa,lswtch,nit,del,dqval,
     .  etol,qrtol,seref,nfile,ntx,ndx,
     .  alat,plat,nft1,nft2,nft3,dx,dy,dz,
     .  atid,spid,z,amass,rmt,rsm,rsmfa,rint,rcut,rham,nphi,lphi,tphi,
     .  colxbs,radxbs,
     .  nxi,lxi,exi,lmxl,lmxa,pin,idmod,piz,idmoz,nsc,nbas,nbasp,nspec,
     .  ipclas,nr,a,
     .  qin,pos0,mpole,dpole,n0,nsx,nbx,
     .  as,ewtol,
     .  ntcpm,tcpm,
     .  lrx,ldyn,tstep,tskip,mdtime,temp,taup,itmxrl,
     .  step,xtol,gtol,lrs,ifrlx,vel,
     .  nvario)
      if (nbasp > nbx) call rx('nbasp gt nbx in '//prgnam)
      if (nspec > nsx) call rx('nspec gt nsx in '//prgnam)
C ... set qrtolm to the smallest positive qrtol
      if (qrtol(2) > 0d0) then
        qrtolm = qrtol(2)
        if (qrtol(1) > 0d0) qrtolm = min(qrtol(1),qrtolm)
      else
        qrtolm = qrtol(1)
      endif
C ... until we get rid of pol and amp, pad out pos from pos0
      if (nbasp > nbas) then
        do  ib = nbas+1, nbasp
          do  i = 1, 3
            pos(i,ib) = pos0(i,ib)
          enddo
        enddo
      endif
      call dpzero(amp,15)
      ltchk = lswtch(1)
      lfrz  = lswtch(3)
      lvint = 1
      lxcqs = lswtch(5)
      llink = lswtch(6)
      lforce= s_ctrl%lfrce
      lrdxc = lswtch(8)
      lwrxc = lswtch(9)
      lcor = 0
      orhoic = 1
      scale = 1d0
      nsp = lsp()+1
      npan = nsc+1
      ipr = iprint()
      zeta = 0

C ... convert spid to dclabl
      do  ic = 1, nspec
        call s8tor8(spid(ic),dclabl(ic))
      enddo

C --- version 2 is not yet setup to read SYMGRP
      symgrp = ' '

C#ifdef MCA | RD | MPRH | MSHOW | MCRHO
C ---- Setup and generate or input density ---
C ... this ensures pnuz set even though no second panel now
      nsc = 1
      call msetp2(spid,z,rmt,rsm,rint,rcut,nphi,lphi,ephi,nxi,lxi,exi,
     .   nsp,lmxl,lmxa,pin,pnu,piz,pnuz,nsc,nspec,n0,nel,el,
     .   atid,pos0,pol,pos,ips,nbas,nbasp,lmaxp,amp,ioff,ioffv0,ioffp,
     .   nri,nvi,ndim,nhs,nla,nll,qval,qsmc)
      nsc = npan-1
C#ifdef MCA
      call mcatm(spid,nbas,nspec,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,nxi,lxi,exi,orhoi,orhoic,orho0,nri,nvi,ioff,ioffv0,
     .  z,rsmfa,pin,qin,irchan,rchar,pnu,orho,orhoc,ov0,setot)
      print 718, setot,seref,setot-seref
  718 format(/' mca: sum etot=',f13.6,'  sum eref=',f13.6,
     .  '  diff=',f13.6)
      ltail = 0
      rhoid = 'MCA'
C#elseC
C      lcor = 0
C      orhoic = 1
C      if (cmdopt('--ra',4,0,outs)) then
C        call mopen(78,'fa','f')
C        rewind 78
C        if (cmdopt('--ra0',5,0,outs)) then
C          call rhdma0(rhoid,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
C     .    rmt,nr,a,lcor,ltail,lrs(4),nxi,lxi,exi,orhoi,orhoic,orho0,
C     .    scale,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .    pos,f,psav,fsav,it,lrs(3),ipr,78)
C        else
C          call rhdmpa(rhoid,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
C     .    rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
C     .    scale,vint,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .    pos,f,psav,fsav,vel,itdyn,
C     .    tstep,time,etot0,zeta,lrs,ipr,78)
C        endif
C        call mclose(78)
C      else
C        if (cmdopt('--rs',4,0,outs)) then
C          ifi = 80
C          call mopen(ifi,'p2','u')
C        else
C          ifi = 79
C          call mopen(ifi,'p1','u')
C        endif
C        rewind ifi
C        if (cmdopt('--rb0',5,0,outs)) then
C          call rhdmp0(rhoid,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
C     .    rmt,nr,a,lcor,ltail,lrs(4),nxi,lxi,exi,orhoi,orhoic,orho0,
C     .    scale,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .    pos,f,psav,fsav,it,lrs(3),ipr,ifi)
C        else
C          call rhdump(rhoid,spid,nbas,nbasp,nsp,npan,n0,ips,lmxl,lmxa,
C     .    rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
C     .    scale,vint,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .    pos,mpole,dpole,f,psav,fsav,vel,itdyn,time,etot0,zeta,lrs,
C     .    ipr,ifi)
C        endif
C        call mclose(79)
C      endif
C#endif
C ---- Manipulate and store density ---
C#ifdefC RD
C      lmxst = 5
C      call stewld(as,ewtol,alat,alat,plat,1d0,lmxst,ixp,dxp)
CC --------- Decompose or overlap charge density -----------------
C      if (cmdopt('--ov',4,0,outs)) then
C        iov = 1
C        if (ltail /= 0) print '(/'' rd: rho already overlapped ...'')'
C      else
C        iov = -1
C        if (ltail == 0) print '(/'' rd: rho already decomposed ...'')'
C      endif
C      if (ltail == 0.and.iov == 1 .or. ltail == 1.and.iov == -1) then
C        call rovlas(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,nbas,alat,pos,
C     .    ips,w(orhoi),ioff,cg,jcg,indxcg,ixp,dxp,orho,iov)
C        ltail = 1-ltail
C      endif
C#elseifC MPRH
C      lcor = 0
C      orhoic = 1
C      call maprho(spid,nbas,nsp,npan,ips,lmxl,lmxa,
C     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,
C     .  ioff,ioffv0,scale,nri,nvi,pnu,pnuz,orho,ov0,orhoc,
C     .  pos,f,psav,fsav,it)
C      rhoid = 'MPRH'
C#elseifC MSHOW
C      call distab(nbas,ips,spid,alat,pos,6)
C      call rx0('MSHOW')
C#elseifC MCRHO
C      lsw = 0
C      if (cmdopt('--psi',5,0,outs)) lsw = 1
C      call plrho(n0,lsw,vals,w(orhoi),nbas,ips,ioff,lmxl,a,nr,rsm,rmt,
C     .  lcor,nxi,lxi,exi,alat,pos,orho,cy,
C     .  nphi,ephi,lphi,nel,nhs,w(oevl),w(1))
C      call rx0('plro')
C#endif
C ---  Dump to file 'p2' ---
      vint = 0
      etot0 = 0
      if (cmdopt('--wa',4,0,outs)) then
        call mopen(78,'fa','f')
        rewind 78
        call rhdmpa(rhoid,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
     .    rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,
     .    ov0,scale,vint,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
     .    pos,f,psav,fsav,vel,itdyn,
     .    tstep,time,etot0,zeta,lrs,ipr,-78)
        call mclose(78)
      else
C#ifdef MCA
        call mopen(80,'p1','u')
C#elseC
C        call mopen(80,'p2','u')
C#endif
        rewind 80
        call rhdump(rhoid,spid,nbas,nbasp,nsp,npan,n0,ips,lmxl,lmxa,
     .    rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,
     .    ov0,scale,vint,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
     .    pos,mpole,dpole,f,psav,fsav,vel,itdyn,time,etot0,zeta,lrs,
     .    ipr,-80)
        call mclose(80)
      endif
      end
C#endif MCA | RD | MPRH

C#ifdefC MCFIT
CC --- make and store OCF ---
C      ifi = 61
C      call mopen(ifi,'t1','u')
C      rewind ifi
C      do  20  is = 1, nspec
C        nn = 0
C        rsm(is) = 0d0
C        call on2gen(n0,lxi(1,is),exi(1,is),nxi(is),lphi(1,is),el,nel,
C     .    llink,elink,tphi(1,1,is),rmt(is),rsm(is),0d0,ri1,ri2,nph,ifi)
C   20 continue
C      call mclose(ifi)
C
CC --- make and store TCF ---
C   29 call query('TCF?',-1,0)
C      ifi = 62
C      call mopen(ifi,'t2','u')
C      rewind ifi
C      do  24  is = 1, nspec
C        do  25  js = 1, is
C
CC ...     Extract TCF parms specific to this pair (of default parms)
C          do  27  i = 2, ntcpm
C          ip = i
C   27     if (is == nint(tcpm(1,i)) .and. js == nint(tcpm(2,i)) .or.
C     .        is == nint(tcpm(2,i)) .and. js == nint(tcpm(1,i))) goto 26
C          ip = 1
C   26     continue
C          nbisi(1) = nint(tcpm(3,ip))
C          nbisi(2) = nint(tcpm(4,ip))
C          nbisi(3) = nint(tcpm(5,ip))
C          nalftc   = nint(tcpm(6,ip))
C          adec     = tcpm(7,ip)
C          wftx     = tcpm(8,ip)
C          ncupl    = nint(tcpm(9,ip))
C          ndust    = nint(tcpm(10,ip))
C
C          d0 = rmt(is)+rmt(js)
CC ... No TCF unless at least 2 atoms of a species or --allp set
C          if (.not. cmdopt('--allp',6,0,outs)) then
C          i = 0
C          j = 0
C          do  28  ib = 1, nbas
C            if (atid(ib) == spid(is) .and. i == 0) i = ib
C            if (atid(ib) == spid(js) .and. (j == 0 .or. j == i)) j = ib
C   28     continue
C          if (j == i) goto 25
C          endif
C          if (nbisi(3) == 1) d0 = dsqrt((pos0(1,i)-pos0(1,j))**2+
C     .      (pos0(2,i)-pos0(2,j))**2+(pos0(3,i)-pos0(3,j))**2)
C          rsm(is) = 0d0
C          rsm(js) = 0d0
C          call hy2gen(lxi(1,is),exi(1,is),nxi(is),lxi(1,js),exi(1,js),
C     .      nxi(js),lphi(1,is),el,nel,lphi(1,js),el,nel,n0,
C     .      llink,elink,tphi(1,1,is),elink,tphi(1,1,js),rmt(is),rmt(js),
C     .      0d0,0d0,rsm(is),rsm(js),
C     .      d0,nbisi(3),adec,wftx,nalftc,
C     .      nbisi(1),nbisi(2),rmt(is),rmt(js),ncupl,ndust,ifi)
C
CC --- Uncomment following for plot ---
CC          rewind 61
CC          npwr = 0
CC          is = 1
CC          js = 1
CC          ie = 1
CC          je = 1
CC          lp1 = lphi(ie,is)
CC          lp2 = lphi(je,js)
CC          call defrr(otspec,  100*ntx)
CC          call defrr(otdata,  ndx)
CC          call hyfinp(1,w(otspec),w(otdata),ntx,ndx)
CC          print *, 'enter dir,d'
CC          read(*,*) dr,d
CC          dr1 = dsqrt(dr(1)**2+dr(2)**2+dr(3)**2)
CC          do  7  m = 1, 3
CC    7     dr(m) = dr(m)*d/dr1
CC          nlm1 = (lp1+1)**2
CC          nlm2 = (lp2+1)**2
CC          nf1 = 0
CC          do  17  k = 1, nxi(is)
CC   17     nf1 = nf1 + (lxi(k,is)+1)**2
CC          nf2 = 0
CC          do  18  k = 1, nxi(js)
CC   18     nf2 = nf2 + (lxi(k,js)+1)**2
CC          call defrr(oc1,     nf1*nlm1*nlm2)
CC          call defrr(oc2,     nf2*nlm1*nlm2)
CC          call hyfget(dr,rmt(is),el(ie),nlm1,rmt(js),el(je),nlm2,
CC     .      w(oc1),w(oc2),nf1,nf2,nxi(is),lxi(1,is),exi(1,is),
CC     .      nxi(js),lxi(1,js),exi(1,js),w(otspec),w(otdata))
CCC          call ious(is,oul1,osl1,a1,r1,nris,lmxa1,83)
CCC          call ious(js,oul2,osl2,a2,r2,nrjs,lmxa2,83)
CC          call hy2plt(lxi(1,is),exi(1,is),nxi(is),w(1),npwr,
CC     .      lxi(1,js),exi(1,js),nxi(js),w(1),npwr,dr,rmt(is),
CC     .      rmt(js),rsm(is),rsm(js),el(ie),el(je),w(oc1),w(oc2),
CC     .      nf1,nf2,nlm1,nlm2)
C
C   25 continue
C   24 continue
C      call mclose(ifi)
C      call rx0('mcfit')
C      end
C#endif MCFIT

C#ifdefC MCE | MCXBS
C* Use hyfrd to allocate exact amount of workspace
Cc --------- input tables -----------------
CC      call hyfinp(nfile,w(otspec),w(otdata),ntx,ndx)
C      tabs = ' t1 t2'
C      call hyfrd(tabs,otspec,otdata,i,j)
C
CC --------- Setup for structure constants -----------------
C      lmxst = 5
C      call stewld(as,ewtol,alat,alat,plat,1d0,lmxst,ixp,dxp)
C
CC --------- Other setup -----------------
C      if (lforce /= 1 .and. (ldyn /= 0 .or. lrx > 0))
C     .  call rx('no atomic movement without forces')
C
Cc --------- setup ----------------------
C*#    do 88 iamp=1,namp
C      iamp = 1
CCL      write(71,716) id(1:60)
CCL  716 format(' mce ',a)
C*# added nsp
C      call msetp2(spid,z,rmt,rsm,rint,rcut,nphi,lphi,ephi,nxi,lxi,exi,
C     .   nsp,lmxl,lmxa,pin,pnu,piz,pnuz,nsc,nspec,n0,nel,el,
C     .   atid,pos0,pol,pos,ips,nbas,nbasp,lmaxp,amp(iamp),ioff,ioffv0,
C     .   ioffp,nri,nvi,ndim,nhs,nla,nll,qval,qsmc)
C#ifdefC MCXBS
CC --------- make bs. file for xbs -----------------
C      if (.not. cmdopt('--rs',4,0,outs)) then
C        call mkbs(atid,spid,ips,alat,nbas,nspec,z,rmt,pos,
C     .            colxbs,radxbs)
C        if (cmdopt('--bonds',7,0,strn)) then
C          dist0 = 5
C          if (cmdopt('-d=',3,0,strn)) then
C            j = 3
C            if (.not. a2bin(strn,dist0,4,0,' ',j,len(strn)))
C     .        call rx('mcxbs: bad value -d')
C          endif
C          call bonds(dist0,nbas,pos,alat,atid)
C        endif
C        goto 1
C      endif
C#endifC
C*# increment qval for Fukui function
C      qval = qval + dqval
C*# added llink
C      if (llink == 1) then
C        call defrr(ochadd, 125*nspec)
C        call msetpl(n0,nel,spid,nspec,nbas,rmt,nphi,ips,llink,
C     .    elink,tphi,lphi,el,nelx,nhl,w(ochadd))
C      else
C        nelx = nel
C        nhl = 0
C        ochadd = 1
C        obl = 1
C        oblp = 1
C      endif
CCL      write(71,715) beta1,beta2,dx,dy,dz,alfa
CCL  715 format(' mce bet',2f8.3,'   dxyz',3f8.3,'   alf',2f8.4)
C      njj = ioffp(nbas+1)
C      call setgrp(symgrp,ng,og,oag)
C      call symcls(ng,w(og),w(oag),nbas,pos,ips,spid,nclas,ipclas)
C      call dpzero(rsm0,nspec)
C      if(ipr >= 1) call tm(' ')
C
Cc --------- check tables -----------------
C      if(ltchk == 1) then
C      write(6,444)
C  444 format(/' Checking tables ...')
C      call oncchk(nspec,rmt,rsm,nphi,lphi,ephi,nxi,lxi,exi,spid,
C     .  n0,w(otspec))
C      call hyfchk(nspec,rmt,rsm,nphi,lphi,ephi,nxi,lxi,exi,spid,
C     .  n0,w(otspec))
C      endif
C
Cc --------- def permanent arrays ------------
C      call defrr(orhoi,   nri*nsp)
C      call defrr(opoti,   nri*nsp)
C      call defrr(opot0,   nvi*nsp)
C      call defrr(ozeta,   nri*nsp)
C      call defrr(ozetxc,  nri*nsp)
C      call defrr(ovval,   nvi*nsp)
C      call defrr(osxi0,   nri*nsp)
C      nv=(nelx*(nelx+1))/2
C      if(nv > nvx) call rx('mce: nv gt nvx')
C      call defrr(ohvecs,  nla*4*nv*nsp)
C      call defrr(osvecs,  nla*4*nv*nsp)
C
CC --- Read from file 79 ---
C      ifi = 1
C      if (procid == master) then
C        if (cmdopt('--rs',4,0,outs)) then
C          ifi = 80
C          call mopen(ifi,'p2','u')
C        else
C          ifi = 79
C          call mopen(ifi,'p1','u')
C        endif
C      endif
C      call rhdump(rhoid,spid,nbas,nbasp,nsp,npan,n0,ips,lmxl,lmxa,
C     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
C     .  scale,vint0,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .  pos,mpole,dpole,f,psav,fsav,vel,itdyn,time,etot0,zeta,lrs,ipr,
C     .  ifi)
C      if (procid == master) then
C        call mclose(ifi)
C      endif
C      call shopos(nbas,nbasp,pos,mpole,dpole,atid)
C
CC --- Set up for MD ---
C      if (ldyn /= 0) then
C        if (ldyn == 1) then
C          ensbl = 'NVE'
C        elseif (ldyn == 2) then
C          ensbl = 'NVT'
C        else
C          call rx0('LMCE: DYN must be 1 (NVE) or 2 (NVT)')
C        endif
C        time0 = tstep
C        ierr = 0
C        if (.not. cmdopt('--st',4,0,outs)) then
C          ifi = 85
C          call mopen(ifi,'strt','f')
C          call iostrt(nbas,nf,itdyn,pos,vel,eps,zeta,zacc,
C     .                zsqacc,veps,time0,ifi,ierr)
C          call mclose(ifi)
C          time = time0
C          call zercmv(nbas,vel,amass,ips)
C        endif
C        if (cmdopt('--st',4,0,outs) .or. ierr == 1) then
C          if (ierr == 1 .and. iprint() >= 10) then
C            print *,' LMCE *warning* could not read strt file '//
C     .              'starting new MD'
C          endif
C          call initv(ensbl,nbas,nf,tstep,ips,amass,temp,0d0,taup,
C     .               0d0,nspec,1,0,dclabl,zeta,zacc,zsqacc,veps,vel)
C          itdyn = 1
C          time = 0d0
C        endif
C        logs = zacc * tstep
C        mdtime = time0 + mdtime - tstep
C      endif
C#ifdefC MCXBS
CC --------- make bs. file for xbs -----------------
C      call mkbs(atid,spid,ips,alat,nbas,nspec,z,rmt,pos,
C     .          colxbs,radxbs)
C      if (cmdopt('--bonds',7,0,strn)) then
C        dist0 = 5
C        if (cmdopt('-d=',3,0,strn)) then
C          j = 3
C          if (.not. a2bin(strn,dist0,4,0,' ',j,len(strn)))
C     .        call rx('mcxbs: bad value -d')
C        endif
C        call bonds(dist0,nbas,pos,alat,atid)
C      endif
C    1 continue
C#endifC
C#ifndefC MCXBS
CC --- Set up relaxation parameters ---
C      call defi(oindrx,6*nbas)
C      call rlxstp(s_ctrl,s_site,
C     .  natrlx,nvar,w(oindrx),xyzfrz,pdim)
C      call rlse(oindrx)
C      icom = 0
C      itrlx = 0
C      if (nvar /= 0) then
C        call defi(oindrx,2*natrlx)
C        call defdr(ow,nvar*nvar)
C        call defdr(op,pdim)
C      endif
C      if (ltail == 0) then
C        call rovlas(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,nbas,alat,pos,
C     .    ips,w(orhoi),ioff,cg,jcg,indxcg,ixp,dxp,orho,1)
C      endif
C
Cc ------- * start of iteration loop * ---------------
C      if (lfrz == 0 .and. iprint() > 20) write(6,177)
C  177 format(/' No extra freezing used')
C      if (lfrz == 1) write(6,178)
C  178 format(/' Freeze phi''s, phi-dot''s, core functions, vint')
Cc|       rewind 81
C      it = 0
CC --- Come back here after 1 iteration if not exceeded max or converged
C    1 continue
C      call defrr(odummy,   10)
C      call info2(10,1,0,' ----- start iteration %i of %i -----',it+1,
C     .  nit)
CCL      write(71,730) it,nit
CCL  730 format(' mce ------------  it',i6,'  of',i6,'  ------------')
C      call dpzero(f1,  3*nbasp)
C      call dpzero(f2,  3*nbasp)
C      call dpzero(f3,  3*nbx*2)
C      call dpzero(f4,  3*nbasp)
C
CC ---  Obtain neighbor table for hamiltonian ---
C      call defrr(ontab, nbas+1)
C!       call defrr(oiax,  300*nbas*10)
C!       call defrr(o,     300*nbas*3)
C      allocate(rtab(3,300*nbas))
C      call hpairm(nbas,ips,alat,plat,pos,rham,nhtab,w(ontab),iax,rtab,
C     .            mxcsiz,1,10)
C!       call rlse(oiax)
C!       call defrr(oiax, nhtab*10)
C!       call defrr(ortab, nhtab*3)
C!       call dpcopy(w(o),w(ortab),1,nhtab*3,1d0)
CC ... these needed for hsmatr
CC     call defrr(ontab3, nbas+1)
CC     call defrr(oiax3, nhtab*10)
CC     call defrr(ortab3, nhtab*3)
CC     call hpair3(nbas,ips,w(ontab),w(oiax),rham,w(ortab),
CC    .  w(ontab3),w(oiax3),w(ortab3))
C
Cc --------- interstitial xc energy and potential --------------
C      if(lrdxc == 0) then
C        if (lxcqs == 1) then
C          call dpzero(w(ozetxc),nri*nsp)
C          call defrr(odvxci,   3*nri*nsp)
C          call dpzero(w(odvxci),nri*3*nsp)
C          call defrr(orhoq, (nft1+2)*nft2*nft3*nsp)
C          ovxc0 = 1
C          odvxc0 = 1
C          call ropfxc(nft1,nft2,nft3,nbas,ips,n0,nxi,lxi,exi,rsm,
C     .      alat,plat,pos,0,ioff,w(orhoi),w(orho0),w(ozetxc),w(ovxc0),
C     .      lforce,w(odvxci),w(odvxc0),w(orhoq),rep0,rmu0,f1)
C          call rlse(odvxci)
C        else
C#ifdefC C1
C        call rx('rhoix7: not implemented for crystal')
C#endifC
C          call rhoix7(nbas,nspec,pos,rsm,rint,nxi,lxi,exi,n0,ips,
C     .      dx,dy,dz,ioff,cy,cg,jcg,indxcg,nri,w(orhoi),rep0,rmu0,
C     .      w(ozetxc),f1)
C        endif
C
Cc  read in old results from ropexc....
C      elseif(lrdxc == 1) then
C        call rx('LMCE: warning this branch line closed!')
C        write(6,*) 'read in old results from ropexc ...'
C        call mopen(83,'rsm','u')
C        rewind 83
C        call dpdump(w(ozetxc),nri*nsp,83)
C        read(83) rep0,rmu0
C        call dpdump(f1,3*nbas,83)
C        call mclose(83)
C      endif
Cc  dump out results from ropexc....
C      if(lwrxc == 1) then
C        call rx('LMCE: warning this branch line closed!')
C        write(6,*) 'dump out results from ropexc ...'
C        call mopen(83,'rsm','u')
C        rewind 83
C        call dpdump(w(ozetxc),nri*nsp,-83)
C        write(83) rep0,rmu0
C        call dpdump(f1,3*nbas,-83)
C        call mclose(83)
C        call rx0('mce')
C      endif
Cc  back to normal program
C* spin polarized
C      call qsmuz(rsm,nxi,lxi,exi,n0,nbas,ips,w(orhoi),q0,amom0)
C      call dpdot(w(orhoi),w(ozetxc),nri*nsp,rhovxc)
CC|    if(ipr >= 30) write(6,988) rep0,rmu0,rhovxc
CC|988 format(' rhoep,rhomu,rhovxc=',3f13.6)
C
Cc --------- interstitial estatic potential ----------
C      call defrr(oqmom,    nvi)
C      call defrr(oqmoi,    nvi)
C* made spin pol
C      call vesmom(rmt,z,nr,a,lmxl,nbas,ips,orho,nvi,ioffv0,w(oqmom))
C* added ixp,dxp, made spin pol
C      call vesmoz(nxi,lxi,exi,n0,rmt,rint,lmxl,nbas,alat,pos,ips,
C     .   w(orhoi),ioff,cg,jcg,indxcg,cy,ioffv0,nvi,ixp,dxp,w(oqmoi))
C* made spin pol
C      call vesint(nxi,lxi,exi,rmt,lmxl,n0,nri,nvi,w(orhoi),
C     .  nbas,ips,w(oqmom),w(oqmoi),w(opoti),w(opot0))
CC ---- Potential from point multipoles; print out multipole moments
C      call potpm(lmaxp,nbas,nbasp,n0,ips,mpole,dpole,atid,
C     .           nxi,lxi,pos,nvi,w(opot0),w(oqmom))
Cc --------- subtract sphere integrals from zeta, zetxc ------
C      call dpzero(w(ozeta),    nri*nsp)
C* added ixp,dxp, made spin pol
C      call vtailz(w(orhoi),w(opoti),w(opot0),nxi,lxi,exi,n0,rmt,rint,
C     .  rcut,rsm,lmxl,nbas,nbasp,lmaxp,alat,pos,ips,ioff,ioffv0,cg,jcg,
C     .  indxcg,cy,nri,ixp,dxp,nvi,w(ovval),w(ozeta),w(ozetxc),w(osxi0),
C     .  qsmo,asmo,repsmo,rmusmo,w(oqmom),upm,f1,f2,f4)
C      if (nbas == 0) goto 2
C      if (nsp == 2) call dpscop(w(osxi0),w(osxi0),nri,1,1+nri,1d0)
C      call dpdot(w(orhoi),w(osxi0),nri*nsp,qint)
C
Cc --------- def arrays for atomic properties ---------
C      call defrr(ovab,  4*n0*nbas*nsp*npan)
C      call defrr(ohab,  4*n0*nbas*nsp*npan)
C      call defrr(osab,  4*n0*nbas*nsp*npan)
C      call defrr(opuu,  njj*nsp*npan)
C      call defrr(opus,  njj*nsp*npan)
C      call defrr(opss,  njj*nsp*npan)
C      call defrr(oppnl, n0*18*nbas)
C
Cc --------- sphere potential, potpars ----------
C      call vsphrs(nbas,rmt,z,a,nr,lmxl,ips,orho,w(ovval),ioffv0,
C     .  cg,jcg,indxcg,cy,vnuc,zvnuc,qsph,amsph,repsph,rmusph,rvhsph,
C     .  sumec,sumtc,lmxa,pnu,pnuz,nsc,n0,w(ohab),w(ovab),w(osab),
C     .  w(oppnl),ioffp,w(opuu),w(opus),w(opss),ov0,lfrz)
C
Cc --------- combine terms of potential; choose vint ------
C      if (nsp == 2) call dpscop(w(ozeta),w(ozeta),nri,1,1+nri,1d0)
C      call dpdot(w(orhoi),w(ozeta),nri*nsp,rvhi)
C      call prpot(nri,nxi,lxi,n0,nbas,ips,w(ozeta),w(ozetxc),w(orhoi))
C      call dpadd(w(ozeta),w(ozetxc),1,nri*nsp,1d0)
C      call etsums(rep0,rmu0,q0,repsmo,rmusmo,qsmo,rvhi,
C     .  repsph,rmusph,rvhsph,qsph,rhoep,rhomu,rhovh,q)
C*# ... msm vint call vavg, save for printout
C      vint=(rvhi+rmu0-rmusmo)/(q0-qsmo)
CC     vavg=vint
C      if(lfrz == 1) vint=vint0
C      if (lvint == 0) vint=0d0
C      vintf=vint
Cc|    vintf=0d0
Cc|    scale=1d0
C      if(ipr >= 20) write(6,380) vint, scale, scale
C  380 format(/' vint=',f10.6,'  scale=',f10.6,'  scale zeta by',f10.6)
C      call defrr(ozetf,   nri*nsp)
C      call dpcopy(w(ozeta),w(ozetf),1,nri*nsp,1d0)
C      call dpadd(w(ozeta),w(osxi0),1,nri*nsp,-vint)
C      call dpadd(w(ozetf),w(osxi0),1,nri*nsp,-vintf)
C      call dpcopy(w(ozeta),w(ozeta),1,nri*nsp,scale)
C      call dpcopy(w(ozetf),w(ozetf),1,nri*nsp,scale)
C      call dpcopy(f4,f4,1,3*nbas,vintf)
C      if(ipr >= 40) then
C      write(6,275)
C      do 58 ib=1,nbas
C   58 write(6,288) ib,(f4(m,ib),m=1,3)
C  288 format(i6,3f13.6)
C  275 format(/'    ib     force from vint * interstitial charge')
C      endif
C
Cc --------- harris energy term to add to evals ------
Cc ... uncomment this line to include interaction energy with point
Cc     multipoles
CC      utot=0.5d0*(rhovh+zvnuc+upm)
C      utot=0.5d0*(rhovh+zvnuc)
C      ertot=rhoep-rhomu+utot-rhovh+sumec
C      if(ipr >= 20) write(6,890) zvnuc,rhovh,utot,rhoep,rhomu,
C     .   ertot,sumec,sumtc
C  890 format(/' zvnuc= ',f14.6,'    rhovh= ',f14.6,'    utot= ',f14.6
C     .       /' rhoeps=',f14.6,'    rhomu= ',f14.6,'    ertot=',f14.6,
C     .       /' sumec= ',f14.6,'    sumtc= ',f14.6)
C
Cc --------- loop over panels: dim chden arrays ----------
C      npan=nsc+1
C      sumev=0d0
C      qsph=0d0
C      do 16 ib=1,nbas
C      is=ips(ib)
C  16  call defrr(orhsms(ib),    nr(is)*(lmxl(is)+1)**2*nsp)
C      call defrr(orhsmi,        nri*nsp)
C      do  ipan =1, npan
Cc --------- copy habs and pert matrices calc for 2nd panel ------
C      if(ipan == 2) then
C        call dpscop(w(ohab),w(ohab),4*n0*nbas*nsp,1+4*n0*nbas*nsp,1,1d0)
C        call dpscop(w(osab),w(osab),4*n0*nbas*nsp,1+4*n0*nbas*nsp,1,1d0)
C        call dpscop(w(opuu),w(opuu),njj*nsp,1+njj*nsp,1,1d0)
C        call dpscop(w(opus),w(opus),njj*nsp,1+njj*nsp,1,1d0)
C        call dpscop(w(opss),w(opss),njj*nsp,1+njj*nsp,1,1d0)
C      endif
Cc --------- dimension output chden arrays ------
C      call defrr(oquu,     njj*nsp)
C      call defrr(oqus,     njj*nsp)
C      call defrr(oqss,     njj*nsp)
C      call defrr(orhout,   nri*nsp)
C
Cc --------- initialize output density ---------
C      call dpzero(w(oquu),    njj*nsp)
C      call dpzero(w(oqus),    njj*nsp)
C      call dpzero(w(oqss),    njj*nsp)
C      call dpzero(w(orhout),  nri*nsp)
C      call dpzero(qus0,       3*n0*nbas*nsp)
C
Cc --------- vectors for hamiltonian and overlap -----
C      call defrr(ogh,        nla*2*nelx)
C      call defrr(ogj,        nla*2*nelx)
C      call hsvecs(nelx,el,lmxa,rmt,ips,nbas,w(ohab),w(osab),
C     .  vint,nla,w(ohvecs),w(osvecs),hi,si,w(ogh),w(ogj))
Cc --------- make structure constants and e-derivatives -------
C      nbb1 = 0
C      do  ie = 1, nel
C        nbb1=nbb1+ndim(ie)**2
C      enddo
C      call defrr(ob,      nla*nhs)
C      call defrr(obp,     nbb1)
C      if (llink == 1) then
C        call defrr(obl,     nla*nhl)
C        call defrr(oblp,    nhl*nhl)
C      endif
C      call defrr(oioffa, nbas+1)
C      call defrr(oioffb, nbas*nel+1)
C      call defrr(oioffl, nbas+1)
C      call hoffs(nbas,lphi,lmxa,n0,ips,nel,w(oioffb),w(oioffa),
C     .  w(oioffl),nlbx)
C      call astrx0(nphi,lphi,nel,el,lmxa,pos,n0,nbas,ips,alat,cg,jcg,
C     .  indxcg,cy,nla,nhs,nhl,w(oioffa),w(oioffb),w(oioffl),
C     .  ixp,dxp,w(ob),w(obp),w(obl),w(oblp))
Cc --------- obtain evecs for both spins --------
C      if (nsp == 2 .and. procid == master) then
C        rewind 95
C      endif
C      do  isp = 1, nsp
Cc --------- make interstitial fix matrix and qi^-qi ----
C        call defrr(oh,      nhs*nhs)
C        call defrr(os,      nhs*nhs)
C        if (alfa(ipan) > 1d-8 .and. isp == 1) then
C          call defrr(osi,     nhs*nhs)
C          call defrr(osub,    4*n0*nbas*nsp)
C          call defrr(ohvucs,  nla*4*nv*nsp)
C          call defrr(osvucs,  nla*4*nv*nsp)
C          call dpzero(w(osub),4*n0*nbas*nsp)
C          call hsvecs(nelx,el,lmxa,rmt,ips,nbas,w(ohab),w(osub),
C     .         vint,nla,w(ohvucs),w(osvucs),hi,si,w(ogh),w(ogj))
C          call dpzero(w(osi), nhs*nhs)
C          call smatl(el,nel,w(ob),w(obp),w(obl),w(oblp),w(ochadd),
C     .     w(osvucs),si,nbas,isp,ips,w(oioffa),w(oioffb),w(oioffl),
C     .     nhs,nhl,nla,w(osi))
C          call dpzero(w(os), nhs*nhs)
C          call hmvi0(el,nel,nphi,lphi,nxi,lxi,exi,n0,rmt,ips,
C     .     nbas,isp,w(otspec),w(otdata),cg,jcg,indxcg,ioff,
C     .     nhtab,iax,rtab,w(oioffb),w(osxi0),nhs,w(os))
C          call qifix(alfa(ipan),nhs,w(osi),w(os),w(oh))
C          if (procid == master) then
C            rewind 94
C            if (iprint() > 30) then
C              print *, 'LMCE writing hamiltonian to disc'
C            endif
C            write (94) nhs,nhs,11
C            call dpdump(w(oh),nhs*nhs,-94)
C          endif
C          call rlse(osi)
C        elseif (alfa(ipan) > 1d-8 .and. isp == 2) then
C          if (procid == master) then
C            rewind 94
C            read (94)
C            call dpdump(w(oh),nhs*nhs,94)
C          endif
C        else
C          call dpzero(w(oh),  nhs*nhs)
C        endif
C        call hmatl(el,nel,w(ob),w(obp),w(obl),w(oblp),w(ochadd),
C     .   w(ohvecs),hi,nbas,isp,ips,w(ogh),w(ogj),w(opuu),w(opus),
C     .   w(opss),w(oioffa),w(oioffb),w(oioffl),ioffp,nhs,nhl,nla,
C     .   w(os),w(oh))
C        call hmvi0(el,nel,nphi,lphi,nxi,lxi,exi,n0,rmt,ips,
C     .   nbas,isp,w(otspec),w(otdata),cg,jcg,indxcg,ioff,
C     .   nhtab,iax,rtab,w(oioffb),w(ozeta),nhs,w(oh))
C        call dpzero(w(os),  nhs*nhs)
C        call smatl(el,nel,w(ob),w(obp),w(obl),w(oblp),w(ochadd),
C     .   w(osvecs),si,nbas,isp,ips,w(oioffa),w(oioffb),w(oioffl),nhs,
C     .   nhl,nla,w(os))
C
Cc --------- diagonalize ------------------------
C        call defrr(oevl,    nhs*nsp)
C        call defrr(oewt,    3*nhs*nsp)
C        call defrr(ot,      nhs*nhs)
C        call tcget(i)
C        if (i == 1) i = 6
C        if (i == 0) i = -1
C        if (cmdopt('--mull',6,0,outs)) then
C          call defrr(osm,nhs*nhs)
C          call dcopy(nhs*nhs,w(os),1,w(osm),1)
C        endif
C        call defrr(owk,     nhs*11)
C        call pshpr(iprint()+30)
C        call dsev1(nhs,w(oh),w(os),w(owk),i,.true.,.true.,.true.,
C     .             nhs,99d0,nev,w(ot),w(oevl))
C        call poppr
C        call rlse(owk)
C        if (nsp == 2) then
C          if (procid == master) then
C            call dpdump(w(oevl),nhs,-95)
C            call dpdump(w(ot),nhs*nhs,-95)
C          endif
C          if (isp == 1) then
C            call rlse(oh)
C          else
C            call dpscop(w(oevl),w(oevl),nhs,1,1+nhs,1d0)
C            if (procid == master) then
C              rewind 95
C              call dpdump(w(oevl),nhs,95)
C            endif
C            call mpibc1(w(oevl),nhs*2,4,mlog,'lmce','oevl')
C          endif
C        endif
CC end loop over isp
C      enddo
C
Cc --------- pert corr; find fermi level and weighting of states ----
C      wgt=2d0/nsp
C      call evlwss(qval,qsmc,ipan,nhs,w(oevl),del,w(oewt),nstate,nsp)
C      if (alfa(ipan) > 1d-8) then
C        do  isp = 1, nsp
C          call percos(nhs,isp,nsp,nstate(isp),w(os),w(ot),w(oevl))
C        enddo
C        call evlwss(qval,qsmc,ipan,nhs,w(oevl),del,w(oewt),nstate,nsp)
C      endif
C      if (cmdopt('--mull',6,0,outs)) then
C        if (llink == 1) then
C          print *, 'Warning: lmce not set up for mull and link'
C        else
C          call defrr(oc,nhs*nhs)
C          call mcmull(nsp,nhs,nbas,atid,z,n0,pnu,pnuz,lmxa,ips,nel,
C     .                w(oioffb),w(oewt),w(ot),w(osm),w(oc),pos,qmull,qc)
C          call rlse(oc)
C        endif
C      endif
C      call dpcopy(w(ot),w(os),1,nhs*nhs,1d0)
C      call rlse(ot)
C      ot = os
C
Cc --------- plot wave function ----
C      if (cmdopt('--plot',6,0,strn)) then
C        lsw = 0
C        if (cmdopt('--psi',5,0,outs)) lsw = 1
C        call plrho(n0,lsw,vals,w(orhoi),nbas,ips,ioff,lmxl,a,nr,rsm,rmt,
C     .    lcor,nxi,lxi,exi,alat,pos,orho,cy,
C     .    nphi,ephi,lphi,nel,nhs,w(oevl),w(ot))
C        call rx0('Done plot output.')
C      endif
C
Cc --------- accumulate forces and output rho for each spin ----
C      if (nsp == 2 .and. procid == master) then
C        rewind 95
C      endif
C      do  isp = 1, nsp
C        if (nsp == 2) then
C          if (procid == master) then
C            call dpdump(w(odummy),1,95)
C            call dpdump(w(ot),nhs*nhs,95)
C          endif
C          call mpibc1(w(ot),nhs*nhs,4,mlog,'lmce','ot')
C        endif
C
Cc --------- force from asa hamiltonian ---------
C        call hsvecs(nelx,el,lmxa,rmt,ips,nbas,w(ohab),w(osab),
C     .   vintf,nla,w(ohvecs),w(osvecs),hi,si,w(ogh),w(ogj))
C        if (lforce == 1) then
C          call hsmdvl(nbas,isp,nphi,lphi,ips,n0,nhs,nhl,nstate(isp),nla,
C     .     w(oioffa),w(oioffb),w(oioffl),el,nel,w(oevl),w(ot),w(ohvecs),
C     .     w(osvecs),hi,si,w(ochadd),w(ogh),w(ogj),w(opuu),w(opus),
C     .     w(opss),ioffp,ixp,dxp,pos,alat,cg,jcg,indxcg,cy,wgt,w(oewt),
C     .     f3(1,1,isp))
C        endif
C
Cc --------- accumulate output density --------
C        call defrr(ov,   nla)
C        call defrr(os,   nla)
C        call rhutsl(wgt,w(oewt),w(oevl),sumev1,nstate(isp),nbas,isp,
C     .   lmxa,ips,n0,nhs,nla,w(oioffa),w(oioffb),w(oioffl),ioffp,llink,
C     .   w(ob),w(obl),w(ochadd),w(ot),nel,w(ogh),w(ogj),w(ov),w(os),
C     .   w(oquu),w(oqus),w(oqss),qus0)
C        sumev=sumev+sumev1
C        call rhoti0(wgt,w(oewt),nstate(isp),el,nel,nphi,lphi,nxi,lxi,
C     .   exi,rmt,n0,nbas,isp,ips,w(otspec),w(otdata),nhtab,iax,
C     .   rtab,w(oioffb),cg,jcg,indxcg,ioff,nhs,w(ot),w(oh),
C     .   w(orhout),w(ozetf),f3(1,1,isp))
C        call rlse(ov)
C      enddo
Cc ATP added this .. (memory leak?)
C      call rlse(oh)
C
C      call rlse(ob)
C      if(ipan == 1) call dpcopy(w(orhout),w(orhsmi),1,nri*nsp,1d0)
C      if(ipan == 2) call dpadd(w(orhsmi),w(orhout),1,nri*nsp,1d0)
Cc --------- assemble output density in spheres -----------
C      call sprho2(nbas,pos,ng,w(og),w(oag),nclas,ipclas,ips,lmxa,
C     .   lmxl,w(oquu),w(oqus),w(oqss),ioffp,z,rmt,nr,a,pnu,pnuz,
C     .   ipan,nsc,n0,ov0,orhsms,qsph1,cg,indxcg,jcg,cy)
C      qsph=qsph+qsph1
C
Cc --------- shift pnu's to dos center; end loop over panels -----
C      if(ipan == 1) then
C        call pnunew(nbas,nclas,ipclas,nsp,n0,lmxa,ips,rmt,qus0,
C     .    w(ohab),w(osab),idmod,pnu,lfrz)
C      else
C        call pnunew(nbas,nclas,ipclas,nsp,n0,lmxa,ips,rmt,qus0,
C     .    w(ohab),w(osab),idmoz,pnuz,lfrz)
C      endif
C      call rlse(oquu)
CC end "loop" over panels
C      enddo
C
Cc --------- get new scale factor, mix densities -------
C      call dpcopy(w(orhsmi),w(orhsmi),1,nri*nsp,scale)
C      call dpdot(w(orhsmi),w(osxi0),nri*nsp,qint)
C      qsum=qval+qsmc
C      asum=amsph-asmo+amom0
C      rescal=(qsum-qsph)/qint
C      scale=scale*rescal
Cc|    scale=1d0
C      if(ipr > 10) write(6,858) qint,qsph,qsum,qsum-qsph,rescal,scale
C  858 format(/' qint=',f12.6,'       qsph=',f12.6,'      qval=',f12.6
C     .       /' qv-s=',f12.6,'     rescal=',f12.6,'     scale=',f12.6)
C      if(ipr > 10 .and. nsp == 2) write(6,860) amom0-asmo,amsph,asum
C  860 format( ' mint=',f12.6,'       msph=',f12.6,'  tot amom=',f12.6)
C      if(dabs(scale-1d0) > 0.05d0) call rx('mce: unacceptable scale')
CC     if(dabs(scale-1d0) > 0.01d0) call rx('mce: unacceptable scale')
C      i1 = str_pack('mix',-2,s_strn,strn)
C      call pshpr(1)
CC     Defaults in the absence of mixing parameters
C      broy = 0
C      mmix = 5
C      beta1 = 0.2d0
C      beta2 = 0.2d0
C      if (.not. parmxp(it+1,strn,len(strn),broy,mmix,dum,beta1,
C     .  dum,dum,dum,dum,dum,dum,dum,beta2,dum))
C     .  call rx('lmce: parse in parmxp failed')
C      call poppr
C      if (broy /= 0) call rx('lmce not set up for Broyden mixing')
C      call mcmixr(mmix,beta1,beta2,rescal,nbas,ips,rmt,nr,a,lmxl,orhsms,
C     .  orho,nri,w(orhsmi),w(orhoi),qrdiff)
C      call rlse(orhsms(1))
C
Cc --------- harris energy and forces --------------
C      ehtot = sumev + ertot - seref
C      if(ipr >= 05) call awrit3('%N sumev=%d ehtot=%d ehdiff=%d',' ',
C     .                          128,i1mach(2),sumev,sumev+ertot,ehtot)
CC      if(ipr >= 05) write(6,550) sumev,ertot,ehtot
CC  550 format(/' sumev=',f13.6,'    ertot=',f13.6,'    ehtot=',f13.6)
C      call roisym(nbas,pos,ng,w(og),w(oag),nclas,ipclas,ips,
C     .   nxi,lxi,n0,ioff,cy,w(orhoi),f)
C    2 continue
C      call dcopy(3,0d0,0,ft,1)
C      do  ib = 1, nbasp
C        do  m = 1, 3
C          f(m,ib)=f1(m,ib)+f2(m,ib)+f3(m,ib,1)+f3(m,ib,2)+f4(m,ib)
C          ft(m) = ft(m) + f(m,ib)
C        enddo
C      enddo
CC zero centre of mass force ..
C      if (iprint() >= 30) print *, 'Zero centre of mass force ..'
C      do  m = 1, 3
C        do  ib = 1, nbasp
C          f(m,ib) = f(m,ib) - ft(m)/nbasp
C        enddo
C      enddo
C      call dcopy(3,0d0,0,ft,1)
C      do  ib = 1, nbasp
C        do  m = 1, 3
C          ft(m) = ft(m) + f(m,ib)
C        enddo
C      enddo
C      if(ipr >= 20) write(6,389)
C      ftop = 0d0
C      do  ib = 1, nbasp
C        xx = dsqrt(f(1,ib)**2+f(2,ib)**2+f(3,ib)**2)
C        ftop = dmax1(ftop,xx)
C        if(ipr >= 30) write(6,388)
C     .    ib,(f(m,ib),m=1,3),xx,(pos(m,ib),m=1,3)
C      enddo
C      xx = dsqrt(ft(1)**2+ft(2)**2+ft(3)**2)
C      if(ipr >= 20) then
CC        write(6,387) (ft(m),m=1,3),xx
C        write(6,391) ftop, xx
C      elseif(ipr >= 10) then
C        write(6,390) ftop
C      endif
C  387 format(1x,'sum',3g10.2,g8.1)
C  388 format(i4,3f10.6,f8.4,3f10.6)
C  389 format(/'  ib',10x,'total force',10x,'abs val',14x,'pos')
C  390 format(/29x,'fmax=',f8.4)
C  391 format(23x,'fmax = ',f7.4,';',1x,'|sum f_i| = ',g9.2)
Cc --------- decompose, move atoms, overlap -------
C      lskip = 0
C      qrdifm = max(qrdiff(1),qrdiff(2))
C      if (lrx > 0 .or. ldyn > 0) then
Cc       if (qrdiff(1) > qrtol(1) .or. qrdiff(2) > qrtol(2)) then
C        if (qrdifm > qrtolm) then
C          lskip = 1
C          if (iprint() >= 10) then
C            call awrit0(' mce: Skip atom move until self consistent...',
C     .        ' ',60,i1mach(2))
C            call awrit4(' mce: it %i mv=0 e=%1;6,6d qdiff=%2:1,5;5d'//
C     .        ' qtol=%2:1;5d',' ',80,i1mach(2),it,ehtot,qrdiff,qrtol)
C          endif
C        else
Cc --------- Reset # of iterations for s-c mixing ------
C          it = 0
C          call rovlas(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,nbas,alat,pos,
C     .      ips,w(orhoi),ioff,cg,jcg,indxcg,ixp,dxp,orho,-1)
C          if (ldyn > 0) then
Cc --------- Molecular dynamics ---------
C            call verlet(ensbl,itdyn,nbas,nf,amass,tstep,f,ips,dclabl,
C     .                  temp,ehtot,0d0,0d0,0d0,.false.,0d0,taup,0d0,
C     .                  mdtime,time,1d0,pos,vel,zeta,0d0,0d0,logs,zacc,
C     .                  zsqacc,ekin,tkin,cons)
C            itdyn = itdyn + 1
C            ifi = 85
C            call mopen(ifi,'strt','u')
C            call iostrt(nbas,nf,itdyn,pos,vel,eps,zeta,zacc,
C     .                  zsqacc,veps,time,-ifi,ierr)
C            call mclose(ifi)
C          else
Cc --------- Molecular statics ---------
C            call dcopy(3*nbas,f,1,fsav,1)
C            call dcopy(3*nbas,pos,1,psav,1)
C            itrlx = itrlx + 1
C            call relax('MCE',s_ctrl,s_site,s_spec,
C     .        itrlx,w(oindrx),natrlx,nvar,f,w(op),w(ow),0,0,pos,icom)
C          endif
C          call rovlas(nxi,lxi,exi,n0,rmt,rint,a,nr,lmxl,nbas,alat,pos,
C     .      ips,w(orhoi),ioff,cg,jcg,indxcg,ixp,dxp,orho,1)
Cc ------- write to mv.ext file for xbs --------
C          if (lrx > 0) then
C            ifi = 85
C            call mopen(ifi,'mv','f')
C            if (itrlx /= 1) then
C              if (iprint() > 10) print*, 'Appending to mv file ..'
C              call poseof(ifi)
C            endif
C            call awrit2('frame ehf=%;4dRy fmax=%;4da.u.',' ',72,ifi,
C     .        ehtot,ftop)
C            write (ifi,880) ((pos(m,ib),m=1,3),ib=1,nbas)
C  880       format(8f10.5)
C            call mclose(ifi)
C          endif
C        endif
C      else
C        itdyn=0
C      endif
C      if (nbas == 0) goto 88
C
Cc ------ * end of iteration loop * ------------
C      lcor = 0
C      orhoic = 1
C      ltail = 1
C      if (mod(it,lrs(1)) == 0) then
C        if (procid == master) then
C          call mopen(80,'p2','u')
C        endif
C        call rhdump('rho-out',spid,nbas,nbasp,nsp,npan,n0,ips,lmxl,lmxa,
C     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
C     .  scale,vint,nri,nvi,pnu,pnuz,idmod,idmoz,orho,orhoc,
C     .  pos,mpole,dpole,f,psav,fsav,vel,itdyn,time,ehtot,zeta,lrs,
C     .  ipr,-80)
C        if (procid == master) then
C          call mclose(80)
C        endif
C      endif
C
CC   ... Write forces to file
C      if (cmdopt('--wforce=',9,0,strn) .and. procid == master)
C     .    call iopos(.true.,0,strn(10:),nbas,f,s_site)
C
C      call rlse(oqmom)
C      deallocate(iax,rtab)
C
C!       call rlse(oiax)
CCL      write(71,712) qint,qsph,scale,vint,sumev,ertot,ehtot
CCL  712 format(' mce  qint',f11.6,'    qsph',f11.6,'   scl',f11.6,'   vi',
CCL     .  f9.3/' mce  sumev',f13.6,'    er',f15.6,'    eh',f15.6)
CCL      do 64 ib=1,nbas
CCL  64  write(71,488) ib,(pos(m,ib),m=1,3),(f(m,ib),m=1,3)
CCL  488 format(' mce  pf',i4,6f10.5)
C      etot(1) = ehtot
C      etot(2) = 0
C      it = it+1
C      call nwit(nvario,it,nit,it == 1,0,etol,qrtolm,qrdifm,0d0,0d0,
C     .    'cxhi',asum,etot,lsc)
C
C      if (lsc == 2 .and. nit > 1) lsc = 3
C      if (lsc == 2 .and. nit == 1) lsc = 1
C      if (s_ctrl%quit == 4)
C     .  call rx0('lmmc : exit (--quit=band)')
CC
CC     it = it-1
CC     call mcsav(it,nit,ldyn+lrx,nvario,etol,qrtol,qrdiff,asum,ehtot)
CC
C
CC      if(ipr >= 10) call tm('end of iteration')
C      if (icom == 1 .and. lrx > 0) then
C        write (6,*) 'Relaxation complete'
C        goto 71
C      endif
C      if (itrlx == itmxrl .and. lrx > 0) then
C        call awrit1('Relaxation incomplete after %i iterations',' ',128,
C     .              6,it)
C        goto 71
C      endif
C      if (time > mdtime .and. ldyn > 0) then
C        call awrit1('Completed %d (a.u.) MD time',' ',128,6,time)
C        goto 71
C      endif
C      if (lsc > 2 .or. lskip /= 0) then
C        goto 1
C      endif
C
C
C   71 continue
CC     if(ipr >= 10) write(6,*) 'End of iteration loop'
C
Cc ------ end loop over amps -----------
C      call rlse(orhoi)
C      if (llink == 1) call rlse(ochadd)
CC      do 63 i=1,3
CC  63  hi(i)=0d0
CC      if(nbas >= 2) call dpdist(pos(1,1),pos(1,2),3,hi(1))
CC      if(nbas >= 3) call dpdist(pos(1,1),pos(1,3),3,hi(2))
CC      if(nbas >= 3) call dpdist(pos(1,2),pos(1,3),3,hi(3))
CCCL      write(71,788) amp(iamp),hi(1),hi(2),hi(3),ehtot
CCCL  788 format(' mce  done:  amp',f9.4,1x,3f9.4,f16.6)
CCCL      call distab(nbas,ips,spid,alat,pos,71)
CC      call distab(nbas,ips,spid,alat,pos,6)
C   88 continue
C
C      call rlse(otspec)
C#endifC !MCXBS
C
C      call tcx(prgnam)
C
C      end
C#endif MCE
C      subroutine rrrcut(e,e1)
CC- Cut real to one digit before decimal point
C      implicit real*8 (a-h,p-z)
C      ii=dabs(e)/10
C      e1=dabs(e)-ii*10d0
C      if(e < 0d0) e1=-e1
C      end

