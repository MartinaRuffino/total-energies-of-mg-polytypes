      subroutine sthsv1(n0,nel,el,lmxa,rmt,ips,nbas,hab,sab,
     .  vab0,vint,nla,hvecs,svecs,hi,si,gh,gj)
C- Makes vectors for matrices: ASA part hvecs, svecs and hi,si
C  and adds ASA part to perturbation matrices
      implicit real*8 (a-h,p-z), integer(o)
      dimension r(1),lmxa(1),el(1),hi(1),si(1),ips(1),
     .  svecs(nla*4,1),hvecs(nla,4,1),sab(36,1),hab(36,1),
     .  vab0(9,nbas,1),gh(nla*2,1),gj(nla*2,1)

c --- make gh and gj ---
      do 10 ie=1,nel
   10 call svsvec(el(ie),lmxa,rmt,ips,nbas,gh(1,ie),gj(1,ie),nla)


c --------- start loop over lmto energies -----------
      iv=0
      do 20 ie=1,nel
      do 20 je=ie,nel
      iv=iv+1
      call xxsths(el(ie),el(je),lmxa,rmt,ips,nbas,sab,hab,
     .  vab0(1,1,iv),vint,nla,svecs(1,iv),si(iv),hvecs(1,1,iv),hi(iv))
   20 continue
c --------- printout -------------
      if(iprint() < 60) return
      nv=iv
      write(6,410) nv,nla,vint
  410 format(/' sethsv:   nv=',i4,'   nla=',i5,'   vint=',f10.5)
      write(6,400) (hi(iv),iv=1,nv)
  400 format(' hvecs...  hi= ',6f11.3)
      call yyyhsv(hvecs,nv,nla)
      write(6,401) (si(iv),iv=1,nv)
  401 format(/' svecs...  si= ',6f11.3)
      call yyyhsv(svecs,nv,nla)

      end
c ----------- printout subroutine ----
      subroutine yyyhsv(hvecs,nv,nla)
      implicit real*8 (a-h,p-z), integer (o)
      dimension hvecs(nla,4,nv)
      do 10 kk=1,4
      write(6,220) kk
  220 format(' -- kk=',i5,' ---')
      i0=1
      do 44 i=2,nla
      top=0d0
      do 4 iv=1,nv
  4   top=dmax1(top,dabs(hvecs(i,kk,iv)-hvecs(i-1,kk,iv)))
      if(top >= 1d-5) write(6,440) i0,i-1,(hvecs(i-1,kk,iv),iv=1,nv)
  44  if(top >= 1d-5) i0=i
      write(6,440) i0,nla,(hvecs(nla,kk,iv),iv=1,nv)
  440 format(i5,'  to',i4,2x,6f11.5)
  10  continue
      end

c ---------- sub xxsths --------------------------
      subroutine xxsths(e1,e2,lmxa,rmt,ips,nbas,sab,hab,vab0,vint,
     .   nla,svec,si,hvec,hi)
c  vectors for hamiltonian and overlap matrices, sequence kk kj jk jj
      implicit real*8 (a-h,p-z), integer(o)
      dimension rmt(1),lmxa(1),ips(1),svec(nla,4),hvec(nla,4),
     .  wkk(9),wkj(9),wjk(9),wjj(9),skk(9),skj(9),sjk(9),sjj(9),
     .  hkk(9),hkj(9),hjk(9),hjj(9),sab(36,1),hab(36,1),vab0(9,1)
      hi=e2+vint
      si=1.d0
      ila=0
      do 30 ib=1,nbas
      is=ips(ib)
      lmax=lmxa(is)
      call wronkj(e1,e2,rmt(is),lmax,wkk,wkj,wjk,wjj)
      call augskj(e1,e2,rmt(is),lmax,sab(1,ib),skk,skj,sjk,sjj)
      call augskj(e1,e2,rmt(is),lmax,hab(1,ib),hkk,hkj,hjk,hjj)

      do 25 l=0,lmax
      k=l+1
      do 25 m=1,2*l+1
      ila=ila+1
      svec(ila,1)= skk(k) + wkk(k)*si
      svec(ila,2)= skj(k) + wkj(k)*si
      svec(ila,3)= sjk(k) + wjk(k)*si
      svec(ila,4)= sjj(k) + wjj(k)*si
      hvec(ila,1)= hkk(k) + wkk(k)*hi
      hvec(ila,2)= hkj(k) + wkj(k)*hi
      hvec(ila,3)= hjk(k) + wjk(k)*hi
      hvec(ila,4)= hjj(k) + wjj(k)*hi - vab0(l+1,ib)
  25  continue
  30  continue
      end
