C# ifndefC ASA-STRX
Cc --- nstru0 -------
C      subroutine nstru0(nlma,nlmb,nlmx,npow,cg,jcg,indxcg,
C     .   ikl,jkl,ncfx,ip,cf)
CC- Set up Clebsch-Gordans in a way which permits longer vectors.
CC  ikl,jkl point to first,last element for a (ilm,ipow) pair.
CC  cf is coeff*(-1)**lk, ip is position in strux as linear array.
C      implicit real*8 (a-h,p-z)
C      dimension cg(1),jcg(1),indxcg(1),ikl(nlmx,0:npow),
C     .   jkl(nlmx,0:npow),ip(1),cf(1)
C      lmaxa=ll(nlma)
C      lmaxb=ll(nlmb)
C      do 1 ipow=0,npow
C      do 1 ilm=1,nlmx
C  1   jkl(ilm,ipow)=0
Cc ------- loop over cg coeffs: count non-zero elements -----
C      do 10 klm=1,nlmb
C      lk=ll(klm)
C      do 10 mlm=1,nlma
C      lm=ll(mlm)
C      ii=max0(mlm,klm)
C      indx=(ii*(ii-1))/2+min0(mlm,klm)
C      icg1=indxcg(indx)
C      icg2=indxcg(indx+1)-1
C      do 10 icg=icg1,icg2
C      ilm=jcg(icg)
C      ipow=(lm+lk-ll(ilm))/2
C  10  if(ilm <= nlmx.and.ipow <= npow) jkl(ilm,ipow)=jkl(ilm,ipow)+1
Cc ------- set up start pointers ikl ---------
C      ntot=0
C      do 14 ipow=0,npow
C      do 14 ilm=1,nlmx
C      ikl(ilm,ipow)=ntot+1
C      ntot=ntot+jkl(ilm,ipow)
C  14  jkl(ilm,ipow)=ikl(ilm,ipow)-1
C      if(iprint() >= 80) write(6,716) nlma,nlmb,nlmx,npow,ntot
C  716 format(' nstru0:  nlma,nlmb,nlmx,npow',4i5,'    ntot=',i8)
C      if(ntot > ncfx) write(6,*) '------ ntot=',ntot
C      if(ntot > ncfx) call rx('nstru0: increase ncfx')
Cc ------- loop over cg coeffs: store coeffs and pointers ------
C      jj=0
C      do 20 klm=1,nlmb
C      lk=ll(klm)
C      do 20 mlm=1,nlma
C      lm=ll(mlm)
C      jj=jj+1
C      ii=max0(mlm,klm)
C      indx=(ii*(ii-1))/2+min0(mlm,klm)
C      icg1=indxcg(indx)
C      icg2=indxcg(indx+1)-1
C      do 20 icg=icg1,icg2
C      ilm=jcg(icg)
C      ipow=(lm+lk-ll(ilm))/2
C      if(ilm <= nlmx.and.ipow <= npow) then
C        jkl(ilm,ipow)=jkl(ilm,ipow)+1
C        ii=jkl(ilm,ipow)
C        ip(ii)=jj
C        cf(ii)=cg(icg)*(-1d0)**lk
C        endif
C  20  continue
C      end
Cc --- nstru1 -------
C      subroutine nstru1(nlma,nlmb,nlmx,npow,klmb,jcg,indxcg,ikl,jkl)
CC- Call after nstru0, sets up pointers jkl for klmb <= nlmb
C      implicit real*8 (a-h,p-z)
C      dimension jcg(1),indxcg(1),ikl(nlmx,0:npow),jkl(nlmx,0:npow)
C      do 14 ipow=0,npow
C      do 14 ilm=1,nlmx
C  14  jkl(ilm,ipow)=ikl(ilm,ipow)-1
Cc ------- loop over cg coeffs: count elements up to klmb ----
C      do 20 mlm=1,nlma
C      lm=ll(mlm)
C      do 20 klm=1,klmb
C      lk=ll(klm)
C      ii=max0(mlm,klm)
C      indx=(ii*(ii-1))/2+min0(mlm,klm)
C      icg1=indxcg(indx)
C      icg2=indxcg(indx+1)-1
C      do 20 icg=icg1,icg2
C      ilm=jcg(icg)
C      ipow=(lm+lk-ll(ilm))/2
C      if(ilm <= nlmx.and.ipow <= npow) jkl(ilm,ipow)=jkl(ilm,ipow)+1
C  20  continue
Cc ------- inspect --------------
CC|    do 56 ilm=1,nlmx
CC|    write(6,230) ilm,(ikl(ilm,ipow),ipow=1,npow)
CC|56  write(6,231)     (jkl(ilm,ipow),ipow=1,npow)
CC|230 format(i5,'  ikl',8i5)
CC|231 format(5x,'  jkl',8i5)
C      end
Cc --- nstrux -------
C      subroutine nstrux(e,dr,nlma,nlmb,nlmx,npow,ikl,jkl,ip,cf,
C     .   cy,s)
CC- Mol strux, standard definition IV-43. Use setup from nstru0.
Cc  Exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
C      implicit real*8 (a-h,p-z), integer(o)
C      dimension s(1),cy(1),hl(200),psi(0:20),phi(0:20),dr(3),ip(1),
C     .   efac(0:20),ikl(nlmx,0:1),jkl(nlmx,0:1),cf(1)
C      fourpi=16.d0*datan(1.d0)
C      lmax=ll(nlma)+ll(nlmb)
C      nlm=(lmax+1)**2
C      if(lmax > 15) stop '*** change dimensions in nstrux'
C      call dpzero(s,nlma*nlmb)
C      call sylm(dr,hl,lmax,r2)
C      if(r2 < 1.d-10) return
C      call bessl(e*r2,lmax,phi,psi)
C      tol=1d-12
C      ilm=0
C      rfac=dsqrt(r2)
C      xxx=1.d0/r2
C      do 10 l=0,lmax
C      rfac=rfac*xxx
C      do 10 m=1,2*l+1
C      ilm=ilm+1
C  10  hl(ilm)=rfac*psi(l)*cy(ilm)*hl(ilm)
C      efac(0)=1.d0
C      do 1 l=1,lmax
C  1   efac(l)=-e*efac(l-1)
Cc ------ loop over ipow,ilm -------------
C      do 30 ipow=0,npow
C      do 30 ilm=1,nlm
C      xxx=fourpi*efac(ipow)*hl(ilm)
C      ii=ikl(ilm,ipow)
C      jj=jkl(ilm,ipow)
C      if(dabs(xxx) > tol) then
C      do 11 i=ikl(ilm,ipow),jkl(ilm,ipow)
C  11  s(ip(i))=s(ip(i))+xxx*cf(i)
C      endif
C  30  continue
C      end
Cc --- nstrug -------
C      subroutine nstrug(e,dr,nlma,nlmb,nlmx,npow,ikl,jkl,ip,cf,
C     .   cy,s,gs)
CC- Molecular structure constants and gradient
C      implicit real*8 (a-h,p-z), integer (o)
C      dimension s(1),gs(nlma*nlmb,3),cy(1),dr(3),
C     .   hl(200),efac(0:20),ghl(200,3),ikl(nlmx,0:1),jkl(nlmx,0:1),
C     .   cf(1),ip(1)
C      tol=1d-12
C      if(nlma < 0.or.nlmb < 0) return
C      fourpi=16*datan(1d0)
C      lmax=ll(nlma)+ll(nlmb)
C      if(lmax > 15) call rx('nstrug: increase dims')
C      call dpzero(s,  nlma*nlmb)
C      call dpzero(gs, nlma*nlmb*3)
C
CC --- Make solid Hankel functions HL ---
C      ndh=200
C      call solhg(e,dr,lmax,ndh,hl,ghl,cy)
C      nlm=(lmax+1)**2
C      efac(0)=1d0
C      do 1 l=1,lmax
C  1   efac(l)=-e*efac(l-1)
C
CC --  Make strx ---
C      do 30 ipow=0,npow
C      do 30 ilm=1,nlm
C      aaa=fourpi*efac(ipow)*hl(ilm)
C      if(dabs(aaa) > tol) then
C        do 11 i=ikl(ilm,ipow),jkl(ilm,ipow)
C  11    s(ip(i))=s(ip(i))+aaa*cf(i)
C      endif
C      do 34 ixyz=1,3
C      xxx=fourpi*efac(ipow)*ghl(ilm,ixyz)
C      if(dabs(xxx) > tol) then
C        do 12 i=ikl(ilm,ipow),jkl(ilm,ipow)
C  12    gs(ip(i),ixyz)=gs(ip(i),ixyz)+xxx*cf(i)
C      endif
C  34  continue
C  30  continue
C      end
Cc --- nstrpg -------
C      subroutine nstrpg(e,dr,nlma,nlmb,nlmx,npow,ikl,jkl,ip,cf,cy,
C     .   s,sd,gs,gd)
CC- Gradient of mol strux and energy derivatives
C      implicit real*8 (a-h,p-z), integer (o)
C      dimension s(nlma*nlmb),gs(nlma*nlmb,3),cy(1),dr(3),
C     .  sd(nlma*nlmb),gd(nlma*nlmb,3),ikl(nlmx,0:1),
C     .  jkl(nlmx,0:1),ip(1),cf(1),hl(200),hd(200),
C     .  efac(0:20),ghl(200,3),ghd(200,3)
C      tol=1d-12
C      if(nlma < 0.or.nlmb < 0) return
C      fourpi=16*datan(1d0)
C      lmax=ll(nlma)+ll(nlmb)
C      call dpzero(s,  nlma*nlmb)
C      call dpzero(sd, nlma*nlmb)
C      call dpzero(gs, nlma*nlmb*3)
C      call dpzero(gd, nlma*nlmb*3)
C
CC --- Make solid Hankel functions HL ---
C      ndh=200
C      call solhpg(e,dr,lmax,ndh,hl,ghl,hd,ghd,cy)
C      nlm=(lmax+1)**2
C      efac(0)=1d0
C      do 1 l=1,lmax
C  1   efac(l)=-e*efac(l-1)
C
CC --- Make strx ---
C      do 30 ipow=0,npow
C      do 30 ilm=1,nlm
C      i1=ikl(ilm,ipow)
C      i2=jkl(ilm,ipow)
C      aaa=fourpi*efac(ipow)*hl(ilm)
C      if(dabs(aaa) > tol) then
C        do 11 i=i1,i2
C  11    s(ip(i))=s(ip(i))+aaa*cf(i)
C        endif
C      do 34 ixyz=1,3
C      xxx=fourpi*efac(ipow)*ghl(ilm,ixyz)
C      if(dabs(xxx) > tol) then
C        do 12 i=i1,i2
C  12    gs(ip(i),ixyz)=gs(ip(i),ixyz)+xxx*cf(i)
C        endif
C  34  continue
C      aad=fourpi*efac(ipow)*(hd(ilm)+ipow*hl(ilm)/e)
C      if(dabs(aad) > tol) then
C        do 15 i=i1,i2
C  15    sd(ip(i))=sd(ip(i))+aad*cf(i)
C        endif
C      do 36 ixyz=1,3
C      xxd=fourpi*efac(ipow)*(ghd(ilm,ixyz)+ipow*ghl(ilm,ixyz)/e)
C      if(dabs(xxd) > tol) then
C        do 16 i=i1,i2
C  16    gd(ip(i),ixyz)=gd(ip(i),ixyz)+xxd*cf(i)
C        endif
C  36  continue
C  30  continue
C      end
Cc --- hstrux -------
C      subroutine hstrux(e,nlma,nlmb,nlmx,npow,ih,ikl,jkl,ip,cf,hl,s)
CC- Strux from reduced strx, standard def IV-43. Setup from nstru0.
Cc  Exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
CC  nlmx dimensions ikl,jkl
C      implicit real*8 (a-h,p-z), integer(o)
C      dimension s(1),hl(nlmx),ip(1),ikl(nlmx,0:1),jkl(nlmx,0:1),cf(1),
C     .  efac(0:20)
C      parameter (tol=1d-12)
C      fourpi=16d0*datan(1d0)
C      lmax=ll(nlma)+ll(nlmb)
C      nlm=(lmax+1)**2
C      if(lmax > 20) call rx('hstrux: lmax gt 20')
C      efac(0)=1d0
C      do 1 l=1,lmax
C    1 efac(l)=-e*efac(l-1)
Cc ------ loop over ipow,ilm -------------
C      call dpzero(s,nlma*nlmb)
C      do 10 ipow=0,npow
C      do 10 ilm=1,nlm
C      xxx=fourpi*efac(ipow)*hl(ilm+ih)
C      ii=ikl(ilm,ipow)
C      jj=jkl(ilm,ipow)
CC     if(dabs(xxx) > tol) then
C      do 11 i=ii,jj
C      k = ip(i)
C   11 s(k)=s(k)+xxx*cf(i)
CC     endif
C   10 continue
C      end
C# endif ASA-STRX
c --- hstrug -------
      subroutine hstrug(e,nlma,nlmb,nlmx,npow,ih,ikl,jkl,ip,cf,hl,s,gs)
C- Structure constants and gradient, given reduced strx.
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(1),gs(nlma*nlmb,3),
     .   hl(1),efac(0:20),ghl(200,3),ikl(nlmx,0:1),jkl(nlmx,0:1),
     .   cf(1),ip(1),xxx(3)
      tol=1d-12
      if(nlma < 0.or.nlmb < 0) return
      fourpi=16*datan(1d0)
      lmax=ll(nlma)+ll(nlmb)
      if(lmax > 15) call rx('hstrug: increase dims')
      call dpzero(s,  nlma*nlmb)
      call dpzero(gs, nlma*nlmb*3)

C --- Gradient of solid Hankel functions ---
      ndh=200
      call slhhg(e,lmax,ndh,hl(1+ih),ghl)
      nlm=(lmax+1)**2
      efac(0)=1d0
      do 1 l=1,lmax
  1   efac(l)=-e*efac(l-1)

C --  Make strx ---
      do 30 ipow=0,npow
      dd3 = fourpi*efac(ipow)
      do 30 ilm=1,nlm
      i1=ikl(ilm,ipow)
      i2=jkl(ilm,ipow)
      aaa=dd3*hl(ilm+ih)
      do 34 ixyz=1,3
   34 xxx(ixyz)=dd3*ghl(ilm,ixyz)
C     if(dabs(aaa) < tol) aaa=0
C     do 37 ixyz=1,3
C  37 if(dabs(xxx(ixyz)) < tol) xxx(ixyz)=0
      do 12 i=i1,i2
        k = ip(i)
        s(k)=s(k)+aaa*cf(i)
        gs(k,1)=gs(k,1)+xxx(1)*cf(i)
        gs(k,2)=gs(k,2)+xxx(2)*cf(i)
        gs(k,3)=gs(k,3)+xxx(3)*cf(i)
   12 continue
   30 continue
      end
C# ifndefC ASA-STRX
Cc --- hstrud -------
C      subroutine hstrud(e,nlma,nlmb,nlmx,npow,ih,ikl,jkl,ip,cf,hl,hd,
C     .  s,sd)
CC- Strux and energy der. from reduced strx, standard def IV-43
Cc  Exp. theorem is:  k(k,r-dr) = sum(m) s(m,k,dr)*j(m,r)
CC  nlmx dimensions ikl,jkl
C      implicit real*8 (a-h,p-z), integer(o)
C      dimension hl(nlmx),hd(nlmx),ip(1),ikl(nlmx,0:1),jkl(nlmx,0:1),
C     .  cf(1),s(nlma,1),sd(nlma,1),efac(0:20)
C      parameter (tol=1d-12)
C      fourpi=16d0*datan(1d0)
C      lmax=ll(nlma)+ll(nlmb)
C      nlm=(lmax+1)**2
C      if(lmax > 20) call rx('hstrux: lmax gt 20')
C      efac(0)=1d0
C      do 1 l=1,lmax
C    1 efac(l)=-e*efac(l-1)
Cc ------ loop over ipow,ilm -------------
C      call dpzero(s, nlma*nlmb)
C      call dpzero(sd,nlma*nlmb)
C      do 10 ipow=0,npow
C      do 10 ilm=1,nlm
C        i1=ikl(ilm,ipow)
C        i2=jkl(ilm,ipow)
C        xxx=fourpi*efac(ipow)*hl(ilm+ih)
C        xxd=fourpi*efac(ipow)*(hd(ilm+ih)+ipow*hl(ilm+ih)/e)
CC       if(dabs(xxx) > tol .or. (dabs(xxd) > tol)) then
C          do 12 i=i1,i2
C            k = ip(i)
C            s(k,1) = s(k,1)  + xxx*cf(i)
C            sd(k,1)= sd(k,1) + xxd*cf(i)
C   12     continue
CC        endif
C   10 continue
C      end
C# endif ASA-STRX
c --- hstrpg -------
      subroutine hstrpg(e,nlma,nlmb,nlmx,npow,ih,ikl,jkl,ip,cf,hl,hd,
     .  s,sd,gs,gd)
C- Gradient and energy derivatives of strux, given reduced strux.
C  nlmx dimensions ikl,jkl
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(nlma*nlmb),gs(nlma*nlmb,3),
     .  sd(nlma*nlmb),gd(nlma*nlmb,3),ikl(nlmx,0:1),
     .  jkl(nlmx,0:1),ip(1),cf(1),hl(1),hd(1),
     .  efac(0:20),ghl(200,3),ghd(200,3),xxx(3),xxd(3)
      tol=1d-12
      if(nlma < 0.or.nlmb < 0) return
      fourpi=16*datan(1d0)
      lmax=ll(nlma)+ll(nlmb)
      call dpzero(s,  nlma*nlmb)
      call dpzero(sd, nlma*nlmb)
      call dpzero(gs, nlma*nlmb*3)
      call dpzero(gd, nlma*nlmb*3)

C --- Make solid Hankel functions HL ---
      ndh=200
      call slhhpg(e,lmax,ndh,hl(1+ih),ghl,hd(1+ih),ghd)
      nlm=(lmax+1)**2
      efac(0)=1d0
      do 1 l=1,lmax
  1   efac(l)=-e*efac(l-1)

C --- Make strx ---
      do 30 ipow=0,npow
      dd3 = fourpi*efac(ipow)
      do 30 ilm=1,nlm
        i1=ikl(ilm,ipow)
        i2=jkl(ilm,ipow)
        aaa=dd3*hl(ilm+ih)
        aad=dd3*(hd(ilm+ih)+ipow*hl(ilm+ih)/e)
        do 34 ixyz=1,3
        xxx(ixyz)=dd3*ghl(ilm,ixyz)
   34   xxd(ixyz)=dd3*(ghd(ilm,ixyz)+ipow*ghl(ilm,ixyz)/e)
        do 12 i=i1,i2
          k = ip(i)
          s(k)=s(k) + aaa*cf(i)
          sd(k)=sd(k)+aad*cf(i)
          gs(k,1)=gs(k,1)+xxx(1)*cf(i)
          gd(k,1)=gd(k,1)+xxd(1)*cf(i)
          gs(k,2)=gs(k,2)+xxx(2)*cf(i)
          gd(k,2)=gd(k,2)+xxd(2)*cf(i)
          gs(k,3)=gs(k,3)+xxx(3)*cf(i)
          gd(k,3)=gd(k,3)+xxd(3)*cf(i)
   12   continue
   30 continue
      end
