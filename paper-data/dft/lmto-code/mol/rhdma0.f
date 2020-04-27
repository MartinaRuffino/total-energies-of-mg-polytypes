      subroutine rhdma0(head,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,irpnu,nxi,lxi,exi,orhoi,orhoic,orho0,scale,
     .  nri,nvi,pnu,pz,idmod,idmoz,orho,orhoc,pos,f,psav,fsav,nit,
     .  irpos,ipr,ifi)
C- I/O for the charge density, ASCII representation
C  This is an ASCII rendition of rhdump
      implicit none
      integer nbas,nsp,n0,nr(1),ips(1),lcor,ifi,ipr,npan,lxi(n0,1),
     .  nxi(1),lmxl(1),lmxa(1),orho(1),orhoc(1),orhoi,orhoic,orho0,
     .  nri,nvi,ltail,nit,irpos,irpnu,idmod(n0,1),idmoz(n0,1)
      character*(*) head
      character*8 spid(1),spid0
      double precision a(1),rmt(1),pnu(n0,nsp,1),pz(n0,nsp,1),exi(n0,1),
     .  pos(3,1),f(3,1),psav(3,1),fsav(3,1)
      character*68 head0,datim
      integer ib,ie,is,jfi,l,nlml,lcor0,lmx0,lmxa0,nlm0,nri0,nvi0,
     .  nbas0,npan0,nsp0,nxi0,lxi0(10),ip,odum,isp,awrite,i1mach,i
      logical lrdall
      double precision exi0(10),rmt0,scale
      character s*80
      real w(1)
      common /w/ w

C --- Input ---
      if (ifi > 0) then
        rewind ifi
        lrdall = nbas <= 0
        read(ifi,'(a)') head0
        head = head0
        read(ifi,'(a)') datim
        read(ifi,501) nbas0,nsp0,npan0,nri0,nvi0,lcor0,ltail,scale,nit
        if (ipr > 30) print 333, 'input ',
     .    nbas0,nsp0,npan0,lcor0,ltail,nri0,nvi0,head0
        if (ipr >= 50 .and. datim /= ' ') print 367, datim
  367   format(' written ',a)
        if (lrdall) then
          nbas = nbas0
          nsp  = nsp0
          npan = npan0
          nri  = nri0
          nvi  = nvi0
        endif
C       if (nbas /= nbas0 .or. nsp /= nsp0 .or. npan /= npan0 .or.
        if (nbas /= nbas0 .or. nsp /= nsp0 .or.
     .      nri /= nri0  .or. nvi /= nvi0)
     .    call rx('RHDMA0: mismatch in nbas,nsp,npan,nri or nvi')
        if (lcor > lcor0) call rx('RHDMA0:  missing core')
        if (npan > npan0) print 347, npan,npan0
  347   format(' rhdump : program has npan=',i1,
     .      ' but file''s contents are for npan=',i1)
        if (ipr >= 30) then
          s = ' Reading:'
          if (irpos == 0) call awrit0('%a  atom positions',s,-80,0)
          if (irpnu == 0) call awrit0('%a  pnu',s,-80,0)
          call awrit0('%a ...',s,-80,-i1mach(2))
        endif
        if (irpos == 0) then
          read(ifi,379)
          call dfdump(pos, 3*nbas,ifi)
          read(ifi,379)
          call dfdump(f,   3*nbas,ifi)
          read(ifi,379)
          call dfdump(psav,3*nbas,ifi)
          read(ifi,379)
          call dfdump(fsav,3*nbas,ifi)
          if (ipr >= 40) then
            print '(''   ib'',16x,''pos'')'
            do  9  ib = 1, nbas
    9       print '(i5, 3f12.6)', ib, (pos(l,ib),l=1,3)
          endif
        else
          call defrr(odum, 3*nbas)
          do  8  ib = 1, 4
          read(ifi,379)
    8     call dfdump(w(odum), 3*nbas,ifi)
          call rlse(odum)
        endif
C ...   Allocate memory for and read coefficient arrays
        call defrr(orhoi, nri*nsp)
        call defrr(orho0, nvi*nsp)
        read(ifi,379)
        call dfdump(w(orhoi),nri*nsp,ifi)
        read(ifi,379)
        call dfdump(w(orho0),nvi*nsp,ifi)
        read(ifi,379)
        if (lcor0 == 1) then
          call defrr(orhoic, nri*nsp)
          call dfdump(w(orhoic),nri*nsp,ifi)
          if (lcor == 0) call rlse(orhoic)
        endif
C ...   Allocate memory while reading density arrays
        if (ipr >= 50) print 364
  364   format(' ib:spc la ll   rmt    nr   lxi   exi',15x,'pnu')
        do  10  ib = 1, nbas
          read (ifi,389)
          read (ifi,502) is,spid0,lmxa0,lmx0,nr(is),rmt0,a(is),nxi0
          if (lrdall) then
            ips(ib) = is
            rmt(is) = rmt0
            nxi(is) = nxi0
            lmxa(is) = lmxa0
            lmxl(is) = lmx0
            spid(is) = spid0
          endif
          write(s,349) ib,spid0,lmxa0,lmx0,rmt0,nr(is)
  349     format(i3,':',a4,2i2,f9.5,i4)
          if (lmx0 /= lmxl(is)) print 346,spid(is),'lmxl',lmxl(is),lmx0
  346     format(' rhdma0 ',a4,': seek rho with ',a4,'=',i1,
     .      ' but file''s contents are to ',i1)
          if (dabs(rmt(is)-rmt0) > 1d-5 .or. nxi(is) /= nxi0 .or.
     .        is /= ips(ib))
     .      call rx('RHDMA0: mismatch in rmt,nxi or ips')
          if (irpnu == 0) then
            do  11  isp = 1, nsp
            call dfdump(pnu(1,isp,ib),lmxa0+1,ifi)
   11       call dfdump(pz(1,isp,ib), lmxa0+1,ifi)
            read (ifi,*) (idmod(l+1,is), l=0,lmxa0)
            read (ifi,*) (idmoz(l+1,is), l=0,lmxa0)
          else
            call rx('rhdma0 not set up to skip reading pnu')
          endif
          read (ifi,503) (lxi0(ie),exi0(ie), ie=1,nxi0)
          if (ipr >= 50) then
            i = max(49,awrite('%a %ni%32p%W%n:1;3d:',s,80,0,
     .        nxi0,lxi0,nxi0,exi0,0,0,0,0))
            ip = awrite('%np %n:1;3d',s,80,0,i,1+lmxa(is),pnu(1,1,ib),
     .        0,0,0,0,0)
            if (ip > 79) s(79:80) = '..'
            call awrit0('%a',s,-80,-i1mach(2))
            if (nsp == 2) call awrit3('%npspin 2: %n:1;3d',
     .        ' ',80,i1mach(2),i-7,1+lmxa(is),pnu(1,2,ib))
            if (npan == 2) call awrit3('%nppanel 2: %n:1;3d',
     .        ' ',80,i1mach(2),i-8,1+lmxa(is),pz(1,1,ib))
            if (npan == 2 .and. nsp == 2) call awrit3('%npspin 2 '//
     .        'panel 2: %n:1;3d',' ',80,i1mach(2),i-15,1+lmxa(is),
     .        pz(1,2,ib))
          endif
          do  12  ie = 1, nxi(is)
            if (lrdall) lxi(ie,is) = lxi0(ie)
            if (lrdall) exi(ie,is) = exi0(ie)
            if (lxi(ie,is) /= lxi0(ie) .or.
     .          dabs(exi(ie,is)-exi0(ie)) > 1d-10)
     .        call rx('RHDMA0: mismatch in lxi or exi')
   12     continue
C     ... Allocate rho, read spin up, down parts separately
          nlml = (lmxl(is)+1)**2
          nlm0 = (lmx0+1)**2
          i = nr(is)*nlml
          l = nr(is)*min(nlml,nlm0)
          call defrr(orho(ib),i*nsp)
          call dpzero(w(orho(ib)),i*nsp)
          call dfdump(w(orho(ib)),l,ifi)
          if (nsp == 2) call dfsdmp(w(orho(ib)),1+i,i+l,ifi)
          if (lcor0 == 1) then
            call defrr(orhoc,nr(is))
            if (nsp == 2)call dfsdmp(w(orhoc(ib)),1+nr(is),2*nr(is),ifi)
          endif
   10   continue
        if (ipr > 50) call ppot3(nri,nxi,lxi,nvi,n0,nbas,ips,
     .    w(orhoi),w(orhoi),w(orhoi),w(orho0),w(orho0),w(orho0),0)
C --- Output ---
      else
        jfi = -ifi
        rewind jfi
        head0 = head
        if (ipr > 30)
     .    print 333, 'output',nbas,nsp,npan,lcor,ltail,nri,nvi,head0
  333   format(/' RHDMA0: ',a6,i3,
     .  ' atoms:  nsp,npan,lcor,ltail,nri,nvi,=',4i2,2i5/' header: ',a)
        write(jfi,'(a)') head0
        datim = ' '
C ...   date and time for unix
        call ftime(datim)
        write(jfi,'(a)') datim
        if (ipr >= 40 .and. datim /= ' ') print *, ' ',datim
        write(jfi,501) nbas,nsp,npan,nri,nvi,lcor,ltail,scale,nit
  501   format(7i6,f15.10,i6)
        write(jfi,379) 'positions'
  379   format('------------- ',a,' ----------------')
        call dfdump(pos, 3*nbas,-jfi)
        write(jfi,379) 'forces'
        call dfdump(f,   3*nbas,-jfi)
        write(jfi,379) 'psav'
        call dfdump(psav,3*nbas,-jfi)
        write(jfi,379) 'fsav'
        call dfdump(fsav,3*nbas,-jfi)
        write(jfi,379) 'rhoi'
        call dfdump(w(orhoi),nri*nsp,-jfi)
        write(jfi,379) 'rho0'
        call dfdump(w(orho0),nvi*nsp,-jfi)
        write(jfi,379) 'rhoic'
        if (lcor == 1) call dfdump(w(orhoic),nri*nsp,-jfi)
        if (ipr >= 50) print 364
        do  20  ib = 1, nbas
          is = ips(ib)
          write(jfi,389) 'pnu,pz,idmod,idmoz,rho', ib
  389     format('------------- ',a,' ib=',i3,'  ----------------')
          nlml = (lmxl(is)+1)**2
          write(jfi,502) is,spid(is),lmxa(is),lmxl(is),nr(is),rmt(is),
     .      a(is),nxi(is)
  502     format(i4,a4,3i5,2f15.10,i5)
          do  21  isp = 1, nsp
          call dfdump(pnu(1,isp,ib),lmxa(is)+1,-jfi)
   21     call dfdump(pz(1,isp,ib), lmxa(is)+1,-jfi)
          write(jfi,'(10i3)') (idmod(l+1,is), l=0,lmxa(is))
          write(jfi,'(10i3)') (idmoz(l+1,is), l=0,lmxa(is))
          write(jfi,503) (lxi(ie,is),exi(ie,is), ie=1,nxi(is))
  503     format(i5,f15.10)
C     ... dump rho in spin up, down blocks separately
          i = nr(is)*nlml
          call dfdump(w(orho(ib)),i,-jfi)
          if (nsp == 2) call dfsdmp(w(orho(ib)),1+i,2*i,-jfi)
          if (lcor == 1) call dfdump(w(orhoc(ib)),nr(is),-jfi)
          if (nsp == 2 .and. lcor == 1)
     .      call dfsdmp(w(orhoc(ib)),1+nr(is),2*nr(is),-jfi)
          if (ipr >= 50) then
            write(s,349) ib,spid(is),lmxa(is),lmxl(is),rmt(is),nr(is)
            i = max(49,awrite('%a %ni%32p%W%n:1;3d:',s,80,0,
     .        nxi(is),lxi(1,is),nxi(is),exi(1,is),0,0,0,0))
            ip = awrite('%np %n:1;3d',s,80,0,i,1+lmxa(is),pnu(1,1,ib),
     .        0,0,0,0,0)
            if (ip > 79) s(79:80) = '..'
            call awrit0('%a',s,-80,-i1mach(2))
            if (nsp == 2) call awrit3('%npspin 2: %n:1;3d',
     .        ' ',80,i1mach(2),i-7,1+lmxa(is),pnu(1,2,ib))
            if (npan == 2) call awrit3('%nppanel 2: %n:1;3d',
     .        ' ',80,i1mach(2),i-8,1+lmxa(is),pz(1,1,ib))
            if (npan == 2 .and. nsp == 2) call awrit3('%npspin 2 '//
     .        'panel 2: %n:1;3d',' ',80,i1mach(2),i-15,1+lmxa(is),
     .        pz(1,2,ib))
          endif
   20   continue
      endif

      end
