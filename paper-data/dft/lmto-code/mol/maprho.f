      subroutine maprho(spid,nbas,nsp,npan,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhic,orho0,
     .  ioff,ioffv0,scale,nri,nvi,
     .  pnu,pz,orho,ov0,orhoc,pos,f,psav,fsav,nitx)
C- Assembles rho from different rho files
C  Just a skeleton for now ...
      implicit none
      integer nbx,n0
      parameter (nbx=200, n0=10)
      integer nbas,nsp,lcor,npan,nri,nvi,ltail,nitx,
     .  nr(nbas),ips(nbas),lxi(n0,nbas),nxi(nbas),lmxl(nbas),lmxa(nbas),
     .  orho(nbas),ov0(nbas),orhoc(nbas),orhoi,orhic,orho0,
     .  idmod(n0,1),idmoz(n0,1),ioff(nbas+1),ioffv0(nbas)
      double precision rmt(nbas),exi(n0,nbas),a(nbas),scale,
     .  pnu(n0,nbas),pz(n0,nbas),
     .  pos(3,nbas),psav(3,nbas),f(3,nbas),fsav(3,nbas)
      integer nbas2,nsp2,lcor2,npan2,nri2,nvi2,ltail2,lastmv2,
     .  nr2(nbx),ips2(nbx),lxi2(n0,nbx),nxi2(nbx),lmxl2(nbx),
     .  lrs2(5),lmxa2(nbx),orho2(nbx),ov02(nbx),orhoc2(nbx),
     .  orhoi2,orhic2,orho02,idmod2(n0,nbx),idmoz2(n0,nbx)
      double precision rmt2(nbx),exi2(n0,nbx),a2(nbx),
     .  pnu2(n0,nbx),pz2(n0,nbx),
     .  pos2(3,nbx),psav2(3,nbx),f2(3,nbx),fsav2(3,nbx),vsav2(3,nbx),
     .  scale2,tau2,tim2,fricfb2,vint2,etot02
      character*8 spid(nbas),spid2(nbx)
      integer ib,jb,ie,is,iri,nspec,ifi,ipr,i,j,ierr,irpos,js,
     .  nspec2,ioff2(nbx+1),ioffv2(nbx),nitx2,ifout,nlst1,nlst2,i1mach,
     .  list1(nbx),list2(nbx),nrix,ich
      logical cmdopt
      double precision gam,q,srfpi,dval,ddot,pi,y0
      character*(68) head2, sout*80,lst1*80,lst2*80
      integer w(1)
      common /w/ w

C --- Misc setup ---
      ifout = i1mach(2)
      pi = 4d0*datan(1d0)
      y0 = 1/dsqrt(4*pi)
C ... For now just read in rho file from 83 ...
      ifi = 83
      call getpr(ipr)

C --- Read in map file ---
      print 345, ifi
  345 format(/' --- maprho: read from unit',i3,' ---')
      nbas2 = 0
      lcor2 = 0
      orhic2 = 1
      call iinit(lrs2,5)
      call mopen(ifi,'p3','u')
      stop 'maprho: update call to rhdump'
C      do  12  i = 1, 2
C   12 call rhdump(head2,spid2,nbas2,nsp2,npan2,n0,ips2,lmxl2,lmxa2,
C     .  rmt2,nr2,a2,lcor2,ltail2,nxi2,lxi2,exi2,orhoi2,orhic2,orho02,
C     .  ov02,scale2,vint2,nri2,nvi2,pnu2,pz2,idmod2,idmoz2,orho2,
C     .  orhoc2,pos2,f2,psav2,fsav2,vsav2,nitx2,lastmv2,tau2,tim2,
C     .  etot02,fricfb2,lrs2,ipr,ifi)
      call mclose(ifi)

      if (nbas > nbx .or. nbas2 > nbx)
     .  call rx('maprho: nbas gt nbx')

C --- Map table ---
C ... Make defaults for map lists
      lst1 = ' '
      call awrit1('1:%i',lst1,80,0,min(nbas,nbas2))
      lst2 = lst1
      if (cmdopt('-mapfrom=',5,0,sout)) then
        i = 0
        call chrpos(sout,'=',9,i)
        if (i <= 8) lst1 = sout(i+2:len(sout))
      endif
      if (cmdopt('-mapto=',5,0,sout)) then
        i = 0
        call chrpos(sout,'=',7,i)
        if (i <= 6) lst2 = sout(i+2:len(sout))
      endif
      call mkilst(lst1,nlst1,list1)
C     if (nlst1 < 0) call rx('error parsing map list '//lst1)
      call mkilst(lst2,nlst2,list2)
C     if (nlst2 < 0) call rx('error parsing map list '//lst2)


C --- Printout ---
      if (ipr >= 30) then
        print *
        call awrit2(' maprho:  nbas=%i; map into %i sites:  '//
     .    lst1,sout,80,0,nbas,nlst1)
        call awrit0(sout,' ',-80,ifout)
        call awrit2(' map file nbas=%i; map from %i sites:  '//
     .    lst2,sout,80,0,nbas2,nlst2)
        call awrit0(sout,' ',-80,ifout)
      endif

C --- Checks ---
      if (nlst1 > nbas) then
        nlst1 = nbas
        call awrit1(' maprho: map lists exceeds nbas.'//
     .  '  Truncate to %i sites',sout,80,ifout,nlst1)
      endif
      if (nlst1 /= nlst2) then
        nlst1 = min(nlst1,nlst2)
        call awrit1(' maprho: map lists are of different lengths.'//
     .  '  Truncate to %i sites',sout,80,ifout,nlst1)
      endif
      if (ltail /= ltail2) call rx('maprho: ltail mismatch')

c --- Make ioff2 ---
      nrix = 0
      do  20  jb = 1, nbas2
        ioff2(jb) = nrix
        js = ips2(jb)
        do  22  j = 1, nxi2(js)
   22   nrix = nrix + (lxi2(j,js)+1)**2
   20 continue
      if (nrix /= nri2) call rx('inconsistent map file')


      call mprh2(ioff,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,w(orhoi),w(orhic),w(orho0),
     .  scale,nri,nvi,pnu,pz,idmod,idmoz,orho,ov0,orhoc,
     .  pos,f,psav,fsav,
     .  ioff2,spid2,nbas2,nsp2,npan2,ips2,lmxl2,lmxa2,
     .  rmt2,nr2,a2,lcor2,ltail2,nxi2,lxi2,exi2,w(orhoi2),
     .  w(orhic2),w(orho02),
     .  scale2,nri2,nvi2,pnu2,pz2,idmod2,idmoz2,orho2,ov02,orhoc2,
     .  pos2,f2,psav2,fsav2,
     .  nlst1,list1,list2)

      call ppot3(nri,nxi,lxi,nvi*0,n0,nbas,ips,w(orhoi),w(orhoi),
     .  w(orhoi),w(orho0),w(orho0),w(orho0),0)

      end
      subroutine mprh2(ioff,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,rhoi,rhic,rho0,
     .  scale,nri,nvi,pnu,pz,idmod,idmoz,orho,ov0,orhoc,
     .  pos,f,psav,fsav,
     .  ioff2,spid2,nbas2,nsp2,npan2,ips2,lmxl2,lmxa2,
     .  rmt2,nr2,a2,lcor2,ltail2,nxi2,lxi2,exi2,rhoi2,rhic2,rho02,
     .  scale2,nri2,nvi2,pnu2,pz2,idmod2,idmoz2,orho2,ov02,orhoc2,
     .  pos2,f2,psav2,fsav2,
     .  nlist,list1,list2)

      implicit none
      integer nbas,nsp,lcor,npan,nri,nvi,ltail,nitx,n0,
     .  nr(nbas),ips(nbas),lxi(n0,nbas),nxi(nbas),lmxl(nbas),lmxa(nbas),
     .  idmod(n0,1),idmoz(n0,1),ioff(nbas+1),orho(nbas),ov0(nbas),
     .  orhoc(nbas),nlist,list1(nlist),list2(nlist)
      double precision rmt(nbas),exi(n0,nbas),a(nbas),scale,
     .  pnu(n0,nbas),pz(n0,nbas),
     .  pos(3,nbas),psav(3,nbas),f(3,nbas),fsav(3,nbas),
     .  rhoi(nri,1),rhic(1),rho0(1)
      integer nbas2,nsp2,lcor2,npan2,nri2,nvi2,ltail2,
     .  nr2(nbas2),ips2(nbas2),lxi2(n0,nbas2),nxi2(nbas2),lmxl2(nbas2),
     .  idmod2(n0,nbas2),idmoz2(n0,nbas2),ioff2(nbas2+1),
     .  lmxa2(nbas2),orho2(nbas2),ov02(nbas),orhoc2(nbas2)
      double precision rmt2(nbas2),exi2(n0,nbas2),a2(nbas2),scale2,
     .  pnu2(n0,nbas2),pz2(n0,nbas2),
     .  pos2(3,nbas2),psav2(3,nbas2),f2(3,nbas2),fsav2(3,nbas2),
     .  rhoi2(nri2,1),rhic2(nri2),rho02(1),edel
      character*4 spid(nbas),spid2(nbas2)
      integer ipr,i,j,ie,je,ke,is,js,ib,jb,ifout,i1mach,iri1,iri2,n1,n2,
     .  iwk1(20),iwk2(20)
      real w(1)
      common /w/ w
      character*1 dig(-1:9), sout*80
      data dig /'-','0','1','2','3','4','5','6','7','8','9'/

C --- For each atom in list, map rhoi and rho ---
      call getpr(ipr)
      ifout = i1mach(2)
      if (ipr >= 30) print 345
  345 format(/' MAPR2: mapping rhoi:'/
     .  '     ib:spc la ll   rmt    nr   lxi     exi')

      if (nsp /= nsp2) call rx('mprh cannot mix spin pol')

      do  10  i = 1, nlist
        ib = list1(i)
        jb = list2(i)
        is = ips(ib)
        js = ips2(jb)
        write(sout,349) 'from',jb,spid2(js),lmxa2(js),lmxl2(js),
     .    rmt2(js),nr2(js),(dig(lxi2(j,js)),j=1,nxi2(js))
        print 350, sout(1:37),(exi2(j,js),j=1,nxi2(js))
        write(sout,349) 'into',ib,spid(is),lmxa(is),lmxl(is),
     .    rmt(is),nr(is),(dig(lxi(j,is)),j=1,nxi(is))
        print 350, sout(1:37),(exi(j,is),j=1,nxi(is))
  349   format(a4,i3,':',a4,2i2,f9.5,i4,2x,7a1)
  350   format(a,6f7.2)

        call dvshel(1,nxi(is),exi(1,is),iwk1,1)
        call dvshel(1,nxi2(js),exi2(1,js),iwk2,1)
        ie = nxi(is)+1
        je = nxi2(js)+1
        do  11  ke = 1, nxi(is)
          ie = ie-1
          je = je-1
C         print *,  ie,exi(1+iwk1(ie),is), je,exi2(1+iwk2(je),js)
          if (je >= 0) then
            if (dabs(exi(1+iwk1(ie),is)-exi2(1+iwk2(je),js)) < 1d-6)
     .        goto 11
          endif
          call awrit1('  ... new exi=%1;6d (absent in mapped rho)',
     .      sout,80,ifout,exi(1+iwk1(ie),is))
          je = je+1
   11   continue
C     ---  Copy rhoi from each energy map ---
        call dpzero(rhoi(ioff(ib)+1,1),ioff(ib+1)-ioff(ib))
        call dpzero(rhoi(ioff(ib)+1,nsp),ioff(ib+1)-ioff(ib))
        do  12  je = 1, nxi2(js)
C ...     Find corresponding energy in ib
          edel = 99d0
          do  14  ke = 1, nxi(is)
            if (edel > dabs(exi2(je,js)-exi(ke,is))) then
              edel = dabs(exi2(je,js)-exi(ke,is))
              ie = ke
            endif
   14     continue
C ...     Error for now if energy missing in ib
          if (edel > 1d-6) then
            call awrit1(' ... exi=%1;6d of mapped rho missing in'//
     .        ' destination',sout,80,0,exi2(je,js))
            call rx(sout)
          endif
C ...      Get proper offsets for mapping
          iri1 = ioff(ib)+1
          do  15  ke = 1, ie-1
   15     iri1 = iri1 + (lxi(ke,is)+1)**2
          n1 = (lxi(ie,is)+1)**2
          iri2 = ioff2(jb)+1
          do  16  ke = 1, je-1
   16     iri2 = iri2 + (lxi2(ke,js)+1)**2
          n2 = (lxi2(je,js)+1)**2
          call dpscop(rhoi2,rhoi,min(n1,n2),iri2,iri1,1d0)
          call dpscop(rhoi2(1,nsp),rhoi(1,nsp),min(n1,n2),iri2,iri1,1d0)
C ...     Printout
          if (lxi2(je,js) > lxi(ie,is))
     .     call awrit3('  ... lxi for energy %1;6,2d truncated from'//
     .     ' %i to %i',sout,80,ifout,exi2(je,js),lxi2(je,js),lxi(ie,is))
          if (lxi2(je,js) < lxi(ie,is))
     .     call awrit3('  ... lxi for energy %1;6,2d expanded from'//
     .     ' %i to %i',sout,80,ifout,exi2(je,js),lxi2(je,js),lxi(ie,is))
   12 continue

C ---   Copy rho0 ---
C        n2 = (lmax2 + 1)**2
C        call dpcopy(rho02(ivi2),rho0(ioffv0(ib)+1),1,n2,1d0)
C        ivi2 = ivi2 + n2
C        if (iprint() >= 30) print *, s(1:ip2)

C ---   Copy rhol ---
        if (nr(is) /= nr2(js)) then
          call awrit2(' ... rhol has %i mesh points but mapped rhol'//
     .      ' has %i points',sout,80,0,nr(is),nr2(js))
          call rx(sout)
        endif
        if (lmxl(is) < lmxl2(js)) then
          call awrit2(' ... truncating rhol to lmxl=%i from'//
     .      ' lmxl=%i',' ',80,ifout,lmxl(is),lmxl2(js))
        endif
        if (lmxl(is) > lmxl2(js)) then
          call awrit2(' ... increasing rhol to lmxl=%i from'//
     .      ' lmxl=%i',' ',80,ifout,lmxl(is),lmxl2(js))
        endif
        n1 = nr(is)*(lmxl(is)+1)**2
        n2 = nr2(js)*(lmxl2(js)+1)**2
        call dpzero(w(orho(ib)),n1*nsp)
        call dpcopy(w(orho2(jb)),w(orho(ib)),1,n1,1d0)
        if (nsp == 2)
     .  call dpscop(w(orho2(jb)),w(orho(ib)),n1,n2+1,n2+n1,1d0)
        call dpcopy(w(ov02(jb)),w(ov0(ib)),1,nr(is)*nsp,1d0)

C ---   Copy pnu, pos,psav,f,fsav ---
        call dpcopy(pnu2(1,jb),pnu(1,ib),1,n0,1d0)
        call dpcopy(pz(1,jb),pz(1,ib),1,n0,1d0)
        call dpcopy(pos2(1,jb),pos(1,ib),1,3,1d0)
        call dpcopy(psav2(1,jb),psav(1,ib),1,3,1d0)
        call dpcopy(f2(1,jb),f(1,ib),1,3,1d0)
        call dpcopy(fsav2(1,jb),fsav(1,ib),1,3,1d0)

        print *
   10 continue

      end
