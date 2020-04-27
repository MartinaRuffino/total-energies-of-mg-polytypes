      subroutine rhdmpa(head,spid,nbas,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
     .  scale,vint,nri,nvi,pnu,pz,idmod,idmoz,orho,orhoc,
     .  pos,f,psav,fsav,vsav,itdyn,tstep,time,etot0,zeta,
     .  lrs,ipr,ifi)
C- I/O for the charge density, ASCII representation
C  This is an ASCII rendition of rhdump
      implicit none
      integer nbas,nsp,n0,nr(1),ips(1),lcor,ifi,ipr,npan,lxi(n0,1),
     .  nxi(1),lmxl(1),lmxa(1),nri,nvi,ltail,idmod(n0,1),idmoz(n0,1),
     .  lrs(5),itdyn,orho(1),orhoc(1),ov0(1),orhoi,orhoic,orho0
      character*(*) head
      character*8 spid(1),spid0
      double precision tstep,time,etot0,zeta,vint,
     .  a(1),rmt(1),pnu(n0,nsp,1),pz(n0,nsp,1),exi(n0,1),
     .  pos(3,1),f(3,1),psav(3,1),fsav(3,1),vsav(3,1)
      character*68 head0,datim
      integer ib,ie,is,jfi,l,nlml,lcor0,lmx0,lmxa0,nlm0,nri0,nvi0,
     .  nbas0,npan0,nsp0,nxi0,lxi0(10),ip,odum,isp,awrite,i1mach,i,
     .  irpos,irpnu,irv0
      logical lrdall
      double precision exi0(10),rmt0,scale,fvsn
      character s*80
      real w(1)
      common /w/ w
      parameter (fvsn=1d0)

C --- Input ---
      if (ifi > 0) then
        irpos = lrs(3)
        irpnu = lrs(4)
        irv0  = lrs(5)
        lrdall = nbas <= 0
        rewind ifi
        read(ifi,'(f12.6)') rmt0
        if (dabs(rmt0-fvsn) > 1d-10)
     .    call fexit2(-1,111,' Exit -1 rhdmpa: '//
     .    'expected file v%,2d but file contains v%,2d',fvsn,rmt0)
        read(ifi,'(a)') head0
        head = head0
        read(ifi,'(a)') datim
        if (ipr > 30) then
          print *, ' '
          call awrit0(' rhdmpa:  input file, id: '//head0//'%a ('//
     .      datim//'%a)',s,-80,-i1mach(2))
        endif
        read(ifi,'(a)') head0
        read(ifi,501) nbas0,nsp0,npan0,nri0,nvi0,lcor0,ltail
        read(ifi,'(a)') head0
        read(ifi,505) scale,vint,time,etot0,zeta,itdyn
        if (ipr >= 30) then
          call awrit2(' itdyn %i, time= %;6d',
     .      ' ',80,i1mach(2),itdyn,time)
          call awrit4(' vint=%;6d%18pscale=%;6d%35pzeta=%,6;6d'//
     .      '%53pehf=%,6;6d',' ',80,i1mach(2),vint,scale,zeta,etot0)
        endif
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
     .    call rx('RHDMPA: mismatch in nbas,nsp,npan,nri or nvi')
        if (lcor > lcor0) call rx('RHDMPA:  missing core')
        if (npan > npan0) print 347, npan,npan0
  347   format(' rhdump : program has npan=',i1,
     .      ' but file''s contents are for npan=',i1)
        if (ipr >= 30) then
          s = ' Reading:'
          if (irpos == 0) call awrit0('%a  atom positions',s,-80,0)
          if (irv0 == 0) call awrit0('%a  atom velocities',s,-80,0)
          if (irpnu == 0) call awrit0('%a  pnu',s,-80,0)
          if (lrdall)       call awrit0('%a  lxi etc',s,-80,0)
          call awrit0('%a ...',s,-80,-i1mach(2))
        endif
        call defrr(odum, 3*nbas)
        if (irpos == 0) then
          read(ifi,379)
          call dfdump(pos, 3*nbas,ifi)
          read(ifi,379)
          call dfdump(f,   3*nbas,ifi)
          read(ifi,379)
          call dfdump(psav,3*nbas,ifi)
          read(ifi,379)
          call dfdump(fsav,3*nbas,ifi)
        else
          do  8  ib = 1, 4
          read(ifi,379)
    8     call dfdump(w(odum), 3*nbas,ifi)
        endif
        read(ifi,379)
        if (irv0 == 0) then
          call dfdump(vsav,3*nbas,ifi)
        else
          call dpdump(w(odum), 3*nbas,ifi)
        endif
        call rlse(odum)
        if (ipr >= 40) then
          print '(''   ib'',16x,''pos'',32x,''velocity'')'
          do  9  ib = 1, nbas
    9     print '(i5, 6f12.6)',ib,(pos(l,ib),l=1,3),(vsav(l,ib),l=1,3)
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
  346     format(' rhdmpa ',a4,': seek rho with ',a4,'=',i1,
     .      ' but file contains lmxl=',i1)
          if (dabs(rmt(is)-rmt0) > 1d-5 .or. nxi(is) /= nxi0 .or.
     .        is /= ips(ib))
     .      call rx('RHDMPA: mismatch in rmt,nxi or ips')
          if (irpnu == 0) then
            do  11  isp = 1, nsp
            call dfdump(pnu(1,isp,ib),lmxa0+1,ifi)
   11       call dfdump(pz(1,isp,ib), lmxa0+1,ifi)
            read (ifi,*) (idmod(l+1,is), l=0,lmxa0)
            read (ifi,*) (idmoz(l+1,is), l=0,lmxa0)
          else
            call rx('rhdmpa not set up to skip reading pnu')
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
     .        call rx('RHDMPA: mismatch in lxi or exi')
   12     continue
          nlml = (lmxl(is)+1)**2
          nlm0 = (lmx0+1)**2
C     ... Allocate rho, read spin up, down parts separately
          i = nr(is)*nlml
          l = nr(is)*min(nlml,nlm0)
          call defrr(orho(ib),i*nsp)
          call dpzero(w(orho(ib)),i*nsp)
          read(ifi,379)
          call dfdump(w(orho(ib)),l,ifi)
          if (nsp == 2) read(ifi,379)
          if (nsp == 2) call dfsdmp(w(orho(ib)),1+i,i+l,ifi)
          i = nr(is)
          call defrr(orhoc(ib),i*nsp)
          if (lcor0 == 1) call dfdump(w(orhoc(ib)),i*nsp,ifi)
          if (lcor /= 1) call rlse(orhoc(ib))
          call defrr(ov0(ib),i*nsp)
          read(ifi,379)
          call dfdump(w(ov0(ib)),i*nsp,ifi)
   10   continue
        if (ipr > 50) call ppot3(nri,nxi,lxi,nvi,n0,nbas,ips,
     .    w(orhoi),w(orhoi),w(orhoi),w(orho0),w(orho0),w(orho0),0)
C --- Output ---
      else
        jfi = -ifi
        rewind jfi
        head0 = head
C ...   date and time for unix
        datim = ' '
        call ftime(datim)
        if (ipr > 30) then
          print *, ' '
          call awrit0(' rhdmpa:  output file, id: '//head0//'%a ('//
     .      datim//'%a)',s,-80,-i1mach(2))
        endif
        write(jfi,'(f12.6)') fvsn
        write(jfi,'(a)') head0
        write(jfi,'(a)') datim
        write(jfi,
     .    '(''   nbas   nsp  npan  nri  nvi    lcor ltail     '')')
        write(jfi,501) nbas,nsp,npan,nri,nvi,lcor,ltail
  501   format(7i6)
        write(jfi,'('' scale,vint,time,etot0,zeta,itdyn'')')
        write(jfi,505) scale,vint,time,etot0,zeta,itdyn
  505   format(5f18.12,i4)
        write(jfi,379) 'positions'
  379   format('------------- ',a,' ----------------')
        call dfdump(pos, 3*nbas,-jfi)
        write(jfi,379) 'forces'
        call dfdump(f,   3*nbas,-jfi)
        write(jfi,379) 'psav'
        call dfdump(psav,3*nbas,-jfi)
        write(jfi,379) 'fsav'
        call dfdump(fsav,3*nbas,-jfi)
        write(jfi,379) 'vsav'
        call dfdump(vsav,3*nbas,-jfi)
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
          write(jfi,389) 'rho', ib
          call dfdump(w(orho(ib)),i,-jfi)
          if (nsp == 2) write(jfi,389) 'rho, spin 2', ib
          if (nsp == 2) call dfsdmp(w(orho(ib)),1+i,2*i,-jfi)
          i = nr(is)
          if (lcor == 1) call dfdump(w(orhoc(ib)),i*nsp,-jfi)
          write(jfi,389) 'v0', ib
          call dfdump(w(ov0(ib)),i*nsp,-jfi)
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
