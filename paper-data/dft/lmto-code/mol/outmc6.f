      subroutine outmc6(ipr,id,tabs,ops,symgrp,nel,el,lswtch,nit1,nit2,
     .  beta1,beta2,del,amp,namp,nfile,ntx,ndx,alfa,alat,plat,nft1,
     .  nft2,nft3,dx,dy,dz,nit0,lrx,gamrx,ldyn,tau,fric,feedb,lrs,ifi)
C  Write general data, MSM's format
      implicit real*8 (a-h,p-z), integer (o)
      dimension el(10),lswtch(20),plat(9),amp(15),lrs(4)
      character*1 id(60,2),symgrp*60,tabs*80,ops*80
      character*120 outs

      outs = ' files  '//tabs
      call awrit0('%a',outs,-80,-i1mach(2))
      call awrit6(' cntrl   {verb} %i  {nit} %i %i {beta} %d %d  '//
     .  '{smear} %d',' ',80,6,ipr,nit1,nit2,beta1,beta2,del)
      call awrit4(' tables  {nfil} %i  {dims} %i %i  {alfa} %d',
     .  ' ',80,6,nfile,ntx,ndx,alfa)
      call awrit1(' switch  {tchk,rel,freeze,tc}%9:1i',' ',80,6,lswtch)
      call awrit3(' xcmesh  {dxyz} %d %d %d',' ',80,6,dx,dy,dz)
      call awrit5(' ftmesh  {alat,plat} %d %9:1d {nabc} %i %i %i',' ',
     .  80,6,alat,plat,nft1,nft2,nft3)
      outs = symgrp
      call awrit0('%11o%1psymgrp  <%f%c%o%a >',outs,80,-6)
      call awrit7(
     .  ' move    {wait} %i  {rx} %i %d  {dyn} %i {tau,fric,feedb}'//
     .  ' %d %d %d',' ',80,6,nit0,lrx,gamrx,ldyn,tau,fric,feedb)
      call awrit4(' rsta    {write} %i  {read} %i  {these pos} %i'//
     .  ' {autosave} %i',' ',80,6,lrs(1),1,lrs(3),lrs(4))
      call awrit1(' elmto  %3:1d',' ',80,6,el)
C     print '('' amp '')'

      end
      subroutine outsp6(spid,z,amass,rmt,rsm,rint,rcut,rham,nphi,lphi,
     .   ephi,nxi,lxi,exi,lmaxl,lmxa,q,p,pz,idmod,nspec,nsx,n0,ifi)
C- Write general data, MSM's format
      implicit real*8 (a-h,p-z), integer (o)
      character*8 spid(1)
      dimension rmt(1),rsm(1),rint(1),rcut(1),rham(1),z(1),
     .  nxi(1),lxi(n0,1),exi(n0,1),nphi(1),ephi(n0,1),lphi(n0,2),
     .  lmaxl(1),p(n0,1),pz(n0,1),lmxa(1),konf(10),idmod(n0,1),
     .  amass(1),q(n0,1),pp(10)
      character s*80

      print *, ' '

      do  10  is = 1, nspec
        s = spid(is)
        do  14  i = 1, nphi(is)
   14   if (lphi(i,is) == -1) lphi(i,is)=9
        call awrit4('%8o%1pspec  <%f%c%a>  {Z} %d {basis} <%ni> '//
     .    '{mass} %d',s,80,-6,z(is),nphi(is),lphi(1,is),amass(is))
        s = spid(is)
        call awrit5('%8o%1p rad  <%f%c%a>  {Rmt,Rsm,Rint,Rcut,Rham}'//
     .    ' %d %d %d %d %d',s,80,-6,rmt(is),rsm(is),rint(is),rcut(is),
     .    rham(is))
        s = spid(is)
        do  15  i = 1, nphi(is)
   15   if (lxi(i,is) == -1) lxi(i,is)=9
        call awrit5('%8o%1p chd  <%f%c%a>  {lmxl} %i {lxi,exi}'//
     .    ' <%ni> %n:1d',s,80,-6,lmaxl(is),nxi(is),lxi(1,is),
     .    nxi(is),exi(1,is))
        s = spid(is)
        do  16  i = 1, lmxa(is)+1
        konf(i) = p(i,is)
   16   pp(i) = p(i,is)-konf(i)
        call awrit4('%8o%1p aug  <%f%c%a>  {konf}'//
     .    ' <%ni> {pnu} %n:1d',s,80,-6,lmxa(is)+1,konf,lmxa(is)+1,pp)
        s = spid(is)
        do  17  i = 1, lmxa(is)+1
        konf(i) = mod(pz(i,is),10d0)
   17   pp(i) = mod(pz(i,is),10d0)-konf(i)
        call awrit4('%8o%1p sc   <%f%c%a>  {konf}'//
     .    ' <%ni> {pnu} %n:1d',s,80,-6,lmxa(is)+1,konf,lmxa(is)+1,pp)
        s = spid(is)
        call awrit4('%8o%1p fa   <%f%c%a>  {konf}'//
     .   ' <%ni> {ql} %n:1d',s,80,-6,lmxa(is)+1,konf,lmxa(is)+1,q(1,is))
   10 continue

      end
      subroutine outatm(atid,pos,pol,nbas,nbx,ifi)
C- Write atomic positions, MSM's format
      implicit real*8 (a-h,p-z), integer (o)
      dimension pos(3,1),pol(3,1)
      character*8 atid(1)
      character s*80

      print *, ' '
      do  10  ib = 1, nbas
        s = atid(ib)
        call awrit1('%8o%1patom  <%f%c%a>%14p%3,6;10D',s,80,-6,
     .    pos(1,ib))
   10 continue
      end
