      subroutine rhdump(head,spid,nbas,nbasp,nsp,npan,n0,ips,lmxl,lmxa,
     .  rmt,nr,a,lcor,ltail,nxi,lxi,exi,orhoi,orhoic,orho0,ov0,
     .  scale,vint,nri,nvi,pnu,pz,idmod,idmoz,orho,orhoc,
     .  pos,mpole,dpole,f,psav,fsav,vsav,itdyn,time,etot0,zeta,
     .  lrs,ipr,ifi)
C- I/O for the charge density (version 1)
C  On output, density consists of nbas sites, each site consisting of
C    the mesh parameters rmt,nr,a;
C    the total density (times r**2) to lmxl, and pnu to lmxa;
C    optionally (lcor=1) a core density;
C  and for the interstitial density,
C    coefficient arrays rhoi and rho0 and the (lxi,exi),
C    and coefficient array rhoic for tails of core density (lcor=1)
C  A header, a string of length 64 characters or less,
C  and the date and time of writing, is stored.
C  ltail=0: rho is not overlapped; ltail=1: rho overlapped.
C
C  Input is in one of three modes.
C  Normal mode (nbas ne 0): parameters nbas,nsp,npan,nri,nvi and
C    rmt(*),nxi(*),lxi(*),exi(*),ips(*) are already known.
C    rhdump aborts on any mismatches.  It reads:
C      scale,vint,itdyn,time,etot0,zeta
C      pos,f,psav,fsav,vsav (unless lrs(3) ne 0 --- lrs(4) for vsav)
C      rhoi,rho0            (and allocates space for them)
C      rhoic                (if lcor0 is 1)
C   .. for each ib
C      pnu,pz,idmod         (unless lrs(4) ne 0)
C      rho(ib),v0(ib)       (and allocates space for them)
C      rhoc(ib)             (if lcor0 is 1)
C  Second mode  (input ips(1) eq 0) like Normal mode but reads
C    rmt(*),nxi(*),lxi(*),exi(*),ips(*) for each ib
C  Third mode (input nbas eq 0) rhdump reads nbas,nsp,npan,nri,nvi
C    and returns.
C  Warning messages when (npan,lcor,lmxl,nxi,exi) mismatch.
      implicit none
      integer nbas,nbasp,nsp,n0,nr(1),ips(1),lcor,ifi,ipr,npan,
     .  lxi(n0,1),nxi(1),lmxl(1),lmxa(1),nri,nvi,ltail,idmod(n0,1),
     .  idmoz(n0,1),lrs(5),itdyn,orho(1),orhoc(1),ov0(1),orhoi,orhoic,
     .  orho0
      character*(*) head
      character*8 spid(1),spid0
      double precision time,etot0,zeta,vint,
     .  a(1),rmt(1),pnu(n0,2,nbasp),pz(n0,2,nbasp),exi(n0,nbasp),
     .  pos(3,nbasp),mpole(nbasp),dpole(3,nbasp),f(3,nbasp),
     .  psav(3,nbasp),fsav(3,nbasp),vsav(3,nbasp)
      character*68 head0,datim
      integer ib,ie,is,jfi,l,nlml,lcor0,lmx0,lmxa0,nlm0,nri0,nvi0,
     .  nbas0,nbasp0,npan0,nsp0,nxi0,lxi0(10),ip,odum,isp,awrite,
     .  irpos,irpnu,irv0,i1mach,i
      logical lrdall
      double precision exi0(10),rmt0,scale,fvsn
      character s*80
      real w(1)
      common /w/ w
      parameter (fvsn=1.05d0)

C For MPI ...
      integer mpipid,procid,master,numprocs
      logical MPI,mlog,cmdopt
      character*80 outs
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

C --- Input ---
      if (ifi > 0) then

        irpos = lrs(3)
        irpnu = lrs(4)
        irv0  = lrs(5)
        lrdall = ips(1) <= 0
        if (procid == master) then
          rewind ifi
          read(ifi) rmt0
        endif
        call mpibc1(rmt0,1,4,mlog,'rhdump','rmt0')
        if (dabs(rmt0-fvsn) > 1d-10)
     .    call fexit2(-1,111,' Exit -1 rhdump: '//
     .    'expected file v%,2d but file contains v%,2d',fvsn,rmt0)
        if (procid == master) then
          read(ifi) head0
        endif
        call mpibcc(head0,1,mlog,'rhdump','head0')
        head = head0
        if (procid == master) then
          read(ifi) datim
        endif
        call mpibcc(datim,1,mlog,'rhdump','datim')
        if (ipr > 30) then
          print *, ' '
          call awrit0(' rhdump:  input file, id: '//head0//'%a ('//
     .      datim//'%a)',s,-80,-i1mach(2))
        endif
        if (procid == master) then
          read(ifi) nbas0,nbasp0,nsp0,npan0,nri0,nvi0,lcor0,ltail
          read(ifi) scale,vint,itdyn,time,etot0,zeta
        endif
        call mpibc1(nbas0, 1,2,mlog,'rhdump','nbas0')
        call mpibc1(nbasp0,1,2,mlog,'rhdump','nbasp0')
        call mpibc1(nsp0,  1,2,mlog,'rhdump','nsp0')
        call mpibc1(npan0, 1,2,mlog,'rhdump','npan0')
        call mpibc1(nri0,  1,2,mlog,'rhdump','nri0')
        call mpibc1(nvi0,  1,2,mlog,'rhdump','nvi0')
        call mpibc1(lcor0, 1,2,mlog,'rhdump','lcor0')
        call mpibc1(ltail, 1,2,mlog,'rhdump','ltail')
        call mpibc1(scale, 1,4,mlog,'rhdump','scale')
        call mpibc1(vint,  1,4,mlog,'rhdump','vint')
        call mpibc1(itdyn, 1,2,mlog,'rhdump','itdyn')
        call mpibc1(time,  1,4,mlog,'rhdump','time')
        call mpibc1(etot0, 1,4,mlog,'rhdump','etot0')
        call mpibc1(zeta,  1,4,mlog,'rhdump','zeta')
        if (ipr >= 10) then
          call awrit2(' itdyn %i, time= %;6d',
     .      ' ',80,i1mach(2),itdyn,time)
          call awrit4(' vint=%;6d%18pscale=%;6d%35pzeta=%,6;6d'//
     .      '%53pehf=%,6;6d',' ',80,i1mach(2),vint,scale,zeta,etot0)
        endif
        if (nbas == 0) then
          nbas = nbas0
          nbasp= nbasp0
          nsp  = nsp0
          npan = npan0
          nri  = nri0
          nvi  = nvi0
          ips(1) = 0
          return
        endif
        if (nbas /= nbas0 .or. nbasp /= nbasp0 .or. nsp /= nsp0 .or.
     .      nri /= nri0  .or. nvi /= nvi0) then
          if (ipr > 1) then
            print 355, 'expected',nbas,nbasp,nri,nvi,nsp
            print 355, 'but read',nbas0,nbasp0,nri0,nvi0,nsp0
          endif
  355     format(1x,a8,'  nbas=',i4,'  nbasp=',i4,'  nri=',i5,
     .                 '  nvi=',i5,'  nsp=',i2)
          call rx('RHDUMP: mismatch in nbas,nbasp,nsp,nri or nvi')
        endif
        if (lcor > lcor0) call rx('RHDUMP:  missing core')
        if (npan > npan0 .and. ipr > 1) print 347, npan,npan0
  347   format(' rhdump : program has npan=',i1,
     .      ' but file''s contents are for npan=',i1)
        if (ipr >= 30) then
          s = ' Reading:'
C         call awrit1('%a (unit %i)',s,-80,0,ifi)
          if (irpos == 0) call awrit0('%a  atom positions',s,-80,0)
          if (irv0 == 0) call awrit0('%a  atom velocities',s,-80,0)
          if (irpnu == 0) call awrit0('%a  pnu',s,-80,0)
          if (lrdall)       call awrit0('%a  lxi etc',s,-80,0)
          call awrit0('%a ...',s,-80,-i1mach(2))
        endif
        call defrr(odum, 3*nbasp)
        if (irpos == 0) then
          if (procid == master) then
            call dpdump(pos,  3*nbasp,ifi)
            call dpdump(f,    3*nbasp,ifi)
            call dpdump(psav, 3*nbasp,ifi)
            call dpdump(fsav, 3*nbasp,ifi)
            call dpdump(dpole,3*nbasp,ifi)
            call dpdump(mpole,  nbasp,ifi)
          endif
          call mpibc1(pos,  3*nbasp,4,mlog,'rhdump','pos')
          call mpibc1(f,    3*nbasp,4,mlog,'rhdump','f')
          call mpibc1(psav, 3*nbasp,4,mlog,'rhdump','psav')
          call mpibc1(fsav, 3*nbasp,4,mlog,'rhdump','fsav')
          call mpibc1(dpole,3*nbasp,4,mlog,'rhdump','dpole')
          call mpibc1(mpole,  nbasp,4,mlog,'rhdump','mpole')
        else
          do  8  ib = 1, 5
          if (procid == master) then
            call dpdump(w(odum), 3*nbasp,ifi)
          endif
    8     continue
          if (procid == master) then
            call dpdump(w(odum),   nbasp,ifi)
          endif
        endif
        if (irv0 == 0) then
          if (procid == master) then
            call dpdump(vsav,3*nbasp,ifi)
          endif
          call mpibc1(vsav,3*nbasp,4,mlog,'rhdump','vsav')
        else
          if (procid == master) then
            call dpdump(w(odum), 3*nbasp,ifi)
          endif
        endif
        call rlse(odum)
        if (ipr >= 30) then
          print '(''   ib'',16x,''pos'',32x,''velocity'')'
          do   ib = 1, nbasp
            print '(i5, 6f12.6)',ib,(pos(l,ib),l=1,3),(vsav(l,ib),l=1,3)
          enddo
        endif
C ...   Allocate memory for and read coefficient arrays
        call defrr(orhoi, nri*nsp)
        call defrr(orho0, nvi*nsp)
        if (procid == master) then
          call dpdump(w(orhoi),nri*nsp,ifi)
          call dpdump(w(orho0),nvi*nsp,ifi)
        endif
        call mpibc1(w(orhoi),nri*nsp,4,mlog,'rhdump','orhoi')
        call mpibc1(w(orho0),nvi*nsp,4,mlog,'rhdump','orho0')
        if (lcor0 == 1) then
          call defrr(orhoic, nri*nsp)
          if (procid == master) then
            call dpdump(w(orhoic),nri*nsp,ifi)
          endif
          call mpibc1(w(orhoic),nri*nsp,4,mlog,'rhdump','orhoic')
          if (lcor == 0) call rlse(orhoic)
        endif
C ...   Allocate memory while reading density arrays
        if (ipr >= 50) print 364
  364   format(' ib:spc la ll   rmt    nr   lxi   exi',15x,'pnu')
        do  10  ib = 1, nbas
          if (procid == master) then
            read (ifi) is,spid0,lmxa0,lmx0,nr(is),rmt0,a(is),nxi0
          endif
          call mpibc1(is,    1,2,mlog,'rhdump','is')
          call mpibcc(spid0, 1,  mlog,'rhdump','spid0')
          call mpibc1(lmxa0, 1,2,mlog,'rhdump','lmxa0')
          call mpibc1(lmx0,  1,2,mlog,'rhdump','lmx0')
          call mpibc1(nr(is),1,2,mlog,'rhdump','nr(is)')
          call mpibc1(rmt0,  1,4,mlog,'rhdump','rmt0')
          call mpibc1(a(is), 1,4,mlog,'rhdump','a(is)')
          call mpibc1(nxi0,  1,2,mlog,'rhdump','nxi0')
          if (lrdall) then
            ips(ib) = is
            rmt(is) = rmt0
            nxi(is) = nxi0
            lmxa(is) = lmxa0
            lmxl(is) = lmx0
            spid(is) = spid0
          endif
          if (ipr > 1) then
            write(s,349) ib,spid0,lmxa0,lmx0,rmt0,nr(is)
            if (lmx0 /= lmxl(is)) then
              print 346,spid(is),'lmxl',lmxl(is),lmx0
            endif
          endif
  349     format(i3,':',a4,2i2,f9.5,i4)
  346     format(' rhdump ',a4,': seek rho with ',a4,'=',i1,
     .      ' but file contains lmxl=',i1)
          if (dabs(rmt(is)-rmt0) > 1d-5 .or. nxi(is) /= nxi0 .or.
     .        is /= ips(ib))
     .      call rx('RHDUMP: mismatch in rmt,nxi or ips')
          if (irpnu == 0) then
            do  11  isp = 1, nsp
              if (procid == master) then
                read (ifi) (pnu(l+1,isp,ib), l=0,lmxa0)
                read (ifi) (pz(l+1,isp,ib), l=0,lmxa0)
              endif
              call mpibc1(pnu(1,isp,ib),lmxa0+1,4,mlog,'rhdump','pnu')
              call mpibc1(pz (1,isp,ib),lmxa0+1,4,mlog,'rhdump','pz')
   11       continue
            if (procid == master) then
              read (ifi) (idmod(l+1,is), l=0,lmxa0)
              read (ifi) (idmoz(l+1,is), l=0,lmxa0)
            endif
            call mpibc1(idmod(1,is),lmxa0+1,2,mlog,'rhdump','idmod')
            call mpibc1(idmoz(1,is),lmxa0+1,2,mlog,'rhdump','idmoz')
          else
            if (procid == master) then
              read (ifi)
              read (ifi)
              if (nsp == 2) then
                read (ifi)
                read (ifi)
              endif
              read (ifi)
              read (ifi)
            endif
          endif
          if (procid == master) then
            read (ifi) (lxi0(ie),exi0(ie), ie=1,nxi0)
          endif
          call mpibc1(lxi0,nxi0,2,mlog,'rhdump','lxi0')
          call mpibc1(exi0,nxi0,4,mlog,'rhdump','exi0')
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
            if (lxi(ie,is) /= lxi0(ie) .or. exi(ie,is) /= exi0(ie))
     .        call rx('RHDUMP: mismatch in lxi or exi')
   12     continue
          nlml = (lmxl(is)+1)**2
          nlm0 = (lmx0+1)**2
C          if (nlm0 /= nlml) then
C            call awrit3(' RHDUMP (warning) ib=%i expected lmxl=%i'//
C     .        'but reading %i',' ',80,i1mach(2),ib,lmxl(is),lmx0)
C          endif
C     ... Allocate rho, read spin up, down parts separately
          i = nr(is)*nlml
          l = nr(is)*min(nlml,nlm0)
          call defrr(orho(ib),i*nsp)
          call dpzero(w(orho(ib)),i*nsp)
          if (procid == master) then
            call dpdump(w(orho(ib)),l,ifi)
          endif
          call mpibc1(w(orho(ib)),l,4,mlog,'rhdump','orho(ib)')
          if (nsp == 2) then
            if (procid == master) then
              call dpsdmp(w(orho(ib)),1+i,i+l,ifi)
            endif
C 2*i because w is real*4:
            call mpibc1(w(orho(ib)+2*i),l,4,mlog,'rhdump','orho(ib)')
          endif
          i = nr(is)
          call defrr(orhoc(ib),i*nsp)
          if (lcor0 == 1) then
            if (procid == master) then
              call dpdump(w(orhoc(ib)),i*nsp,ifi)
            endif
            call mpibc1(w(orhoc(ib)),i*nsp,4,mlog,'rhdump','orhoc(ib)')
          endif
          if (lcor /= 1) call rlse(orhoc(ib))
          call defrr(ov0(ib),i*nsp)
          if (procid == master) then
            call dpdump(w(ov0(ib)),i*nsp,ifi)
          endif
          call mpibc1(w(ov0(ib)),i*nsp,4,mlog,'rhdump','ov0(ib)')
C         call prmx('v0 in rhdump',w(ov0(ib)),nr(is),nr(is),nsp)
   10   continue

C --- Output ---
      elseif (procid == master) then
        jfi = -ifi
        rewind jfi
        head0 = head
C ...   date and time for unix
        datim = ' '
        call ftime(datim)
        if (ipr > 30) then
          print *, ' '
          call awrit0(' rhdump:  output file, id: '//head0//'%a ('//
     .      datim//'%a)',s,-80,-i1mach(2))
        endif
        write(jfi) fvsn
        write(jfi) head0
        write(jfi) datim
        write(jfi) nbas,nbasp,nsp,npan,nri,nvi,lcor,ltail
        write(jfi) scale,vint,itdyn,time,etot0,zeta
        call dpdump(pos,  3*nbasp,-jfi)
        call dpdump(f,    3*nbasp,-jfi)
        call dpdump(psav, 3*nbasp,-jfi)
        call dpdump(fsav, 3*nbasp,-jfi)
        call dpdump(dpole,3*nbasp,-jfi)
        call dpdump(mpole,  nbasp,-jfi)
        call dpdump(vsav, 3*nbasp,-jfi)
        call dpdump(w(orhoi),nri*nsp,-jfi)
        call dpdump(w(orho0),nvi*nsp,-jfi)
        if (lcor == 1) call dpdump(w(orhoic),nri*nsp,-jfi)
        if (ipr >= 50) print 364
        do  20  ib = 1, nbas
          is = ips(ib)
          nlml = (lmxl(is)+1)**2
          write(jfi) is,spid(is),lmxa(is),lmxl(is),nr(is),rmt(is),
     .      a(is),nxi(is)
          do  21  isp = 1, nsp
            write(jfi) (pnu(l+1,isp,ib), l=0,lmxa(is))
            write(jfi) (pz(l+1,isp,ib), l=0,lmxa(is))
   21     continue
          write(jfi) (idmod(l+1,is), l=0,lmxa(is))
          write(jfi) (idmoz(l+1,is), l=0,lmxa(is))
          write(jfi) (lxi(ie,is),exi(ie,is), ie=1,nxi(is))
C     ... dump rho in spin up, down blocks separately
          i = nr(is)*nlml
          call dpdump(w(orho(ib)),i,-jfi)
          if (nsp == 2) call dpsdmp(w(orho(ib)),1+i,2*i,-jfi)
          i = nr(is)
          if (lcor == 1) call dpdump(w(orhoc(ib)),i*nsp,-jfi)
          call dpdump(w(ov0(ib)),i*nsp,-jfi)
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
