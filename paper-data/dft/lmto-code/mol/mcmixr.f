      subroutine mcmixr(mmix,beta1,beta2,scale,nbas,ips,rmt,nr,a,lmxl,
     .   orhut,orhos,nri,rhoix,rhoi,qrdiff)
C- Mix density
C ----------------------------------------------------------------------
Ci Inputs: mmix, beta1,2 (sphere and interstitial mixing), scale, nbas
Ci         ips, rmt, nr, a, lmxl, orhut, orhos (pointers to new and
Ci         sphere density (by L), nri, rhoix, rhoi (coefficients to
Ci         new and old interstitial density)
Co Outputs:
Co         qrdiff
Cr Remarks
Cu Updates
Cu   Adapted from mixrho to include Anderson mixing
C ----------------------------------------------------------------------
      implicit none
      integer nbmx
      parameter( nbmx=2000 )
C Passed
      integer orhut(1),orhos(1)
      integer nbas,mmix,nri,ips(1),nr(1),lmxl(1)
      double precision rmt(1),a(1),rhoi(nri,2),rhoix(nri,2),
     .  qrdiff(2),beta1,beta2,scale
C Local
      integer nsp,lsp,isp,ib,is,nr1,ipr,nlml,i,iprint,nelts,adim,aptr
      integer ornew,orold,orofi,orwgt,oh,oa
      double precision diffat(nbmx),diftop,top,rdiff,diff,xx
C Heap
      real w(1)
      common /w/ w
      call query('beta1',4,beta1)
      call query('beta2',4,beta2)
      call getpr(ipr)
      nsp = lsp() + 1
      if (nbas > nbmx) call rx('mcmixr: nbas gt nbmx')

C --- Mixing of sphere density ---
      diftop = 0
      nelts = 0
      do ib = 1, nbas
C --- Get number of elements to mix in this sphere ---
        is = ips(ib)
        nr1 = nr(is)
        nlml = (lmxl(is)+1)**2
        nelts = nelts + nr1*nlml*nsp
C --- Make the sphere charge density difference ---
        ornew = orhut(ib)
        orold = orhos(ib)
        call defrr(orofi,   nr1)
        call defrr(orwgt,   nr1)
        call defrr(oh,      nr1)
        call rhodif(w(ornew),w(orold),nr1,nlml,a(is),rmt(is),diff)
        diffat(ib) = diff
        call rlse(orofi)
        diftop = dmax1(diftop,diff)
      enddo
      call defrr(oa,nelts*(mmix+2)*2)
      adim = nelts*(mmix+2)
      aptr = 1
C --- Poke elements from each sphere into work array ---
      do ib = 1, nbas
        is = ips(ib)
        nr1 = nr(is)
        nlml = (lmxl(is)+1)**2
        call dpscop(w(orhut(ib)),w(oa),nr1*nlml*nsp,1,aptr,1d0)
        call dpscop(w(orhos(ib)),w(oa),nr1*nlml*nsp,1,aptr+adim,1d0)
        aptr = aptr + nr1*nlml*nsp
      enddo
C --- Mix ---
      if (iprint() >= 30) print *, 'mcmixr: mix sphere density ..'
      call xmxrho(74,nelts,mmix,beta1,w(oa))
C --- Take new density from work array ---
      aptr = 1
      do ib = 1, nbas
        is = ips(ib)
        nr1 = nr(is)
        nlml = (lmxl(is)+1)**2
        call dpscop(w(oa),w(orhos(ib)),nr1*nlml*nsp,aptr+adim,1,1d0)
        aptr = aptr + nr1*nlml*nsp
      enddo
      if (aptr /= nelts+1) call rx(' bug in mcmixr')
      call rlse(oa)

      if (ipr >= 70) print 1
    1 format(' iri    out rhoi    old rhoi      diff       mixed')

C --- Mixing of interstitial density ---
      do i = 1,nri*nsp
        rhoix(i,1) = scale*rhoix(i,1)
      enddo
      top = 0d0
      rdiff = 0d0
      do  isp = 1,nsp
        do  i = 1,nri
          xx = rhoix(i,isp)-rhoi(i,isp)
          rdiff = rdiff + (xx/max(1d0/nsp,dabs(rhoix(i,isp))))**2
          top = dmax1(top,dabs(xx))
          if (ipr >= 70 .and. dabs(xx) > 1d-6)
     .      print 2, i,rhoi(i,isp),rhoix(i,isp),xx,
     .      beta2*rhoix(i,isp)+(1d0-beta2)*rhoi(i,isp)
    2     format(i4,4f12.6)
        enddo
      enddo
C --- Poke density in work array, mix and poke back ---
      call defrr(oa,nri*nsp*(mmix+2)*2)
      call dpcopy(rhoix,w(oa),1,nri*nsp,1d0)
      call dpscop(rhoi,w(oa),nri*nsp,1,1+nri*nsp*(mmix+2),1d0)
      if (iprint() >= 30) print*, 'mcmixr: mix interstitial density ..'
      call xmxrho(70,nri*nsp,mmix,beta2,w(oa))
      call dpscop(w(oa),rhoi,nri*nsp,1+nri*nsp*(mmix+2),1,1d0)
      call rlse(oa)

      rdiff = dsqrt(rdiff/nri/nsp)
      if(ipr >= 40) write(6,3) (diffat(ib),ib = 1,nbas)
    3 format(' at diff=',6f10.5)
      qrdiff(1) = diftop
      qrdiff(2) = rdiff
      if (ipr >= 30) print 4, diftop, rdiff
    4 format(' mcmixr: sdiff=',f12.6,'   rdiff',f12.6)

      if (iprint() >= 40) call query('abort?',-1,0)

      end
      subroutine xmxrho(ifi,nelts,mmix,beta,a)
C- Read from disk any previous iterations
      implicit none
      integer nelts,mmix,ifi
      double precision a(nelts,mmix+2,2),tj(10),norm(100),kpvt(10),
     .  beta,rms2
      integer nelts0,imix,amix,jmix,nmix,i,j,k,iprint,i1mach

      integer procid,master,numprocs,mpipid
      logical T,F,mlog,ioerr,cmdopt,ipr
      character*80 outs
      data T /.true./ F /.false./
      procid = mpipid(1)
      numprocs = mpipid(0)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)
      ipr = iprint() > 0 .and. numprocs > 1

C --- Do the mixing ---
      imix = 0
      if (mmix > 0) then
        if (procid == master) then
          ioerr = F
          rewind(ifi,err=11)
          read(ifi,err=11,end=11) imix, nelts0
          if (ipr) then
            call awrit4(' MCMIXR: read mixfile unit %i, '//
     .        'mmix=%i imix=%i nelts0=%i',
     .        ' ',128,i1mach(2),ifi,mmix,imix,nelts0)
          endif
          goto 12
   11     continue
          if (ipr) then
            call awrit1(' MCMIXR: mixfile unit %i not found.',
     .        ' ',128,i1mach(2),ifi)
          endif
          ioerr = T
   12     continue
        endif
        call mpibc1(ioerr,1,1,mlog,'mcmixr','ioerr')
        if (ioerr) goto 3
        call mpibc1(imix,1,2,mlog,'mcmixr','imix')
        call mpibc1(nelts0,1,2,mlog,'mcmixr','nelts0')
        if (nelts0 == nelts .and. imix <= mmix) then
          goto 2
        else
          goto 3
        endif
C ... Reaching this point means have file holding previous iters
    2   continue
        if (procid == master) then
          ioerr = F
          read(ifi,err=13) (((a(i,j,k), i=1,nelts), j=2, imix+1), k=1,2)
          if (ipr) then
            call awrit2(' MCMIXR: read mixfile unit %i, %i elements.',
     .        ' ',128,i1mach(2),ifi,nelts*imix*2)
          endif
          goto 14
   13     continue
          ioerr = T
          if (ipr) then
            call awrit1(' MCMIXR: mixfile on unit %i read failed.',
     .        ' ',128,i1mach(2),ifi)
          endif
   14     continue
        endif
        call mpibc1(ioerr,1,1,mlog,'mcmixr','ioerr')
        if (ioerr) goto 3
        call mpibc1(a,nelts*(mmix+2)*2,4,mlog,'mcmixr','a')
        goto 4
    3   continue
C ... Reaching this point means file empty or file mismatch
        imix = 0
        goto 4
      endif
C ... imix now set to number of previous iterations
    4 continue
      jmix = imix
      if (beta < 0) jmix = -imix

      call pshpr(1)
      nmix = amix(nelts,jmix,mmix,0,dabs(beta),iprint(),10d0,
     .            norm,kpvt,a,tj,rms2)
      call poppr
      if (iprint() >= 30) then
        call awrit5(' mcmixr: mixed %i of %i, beta=%d, %i elements, '//
     .  'rms diff: %g',' ',120,i1mach(2),nmix,imix,beta,nelts,rms2)
        if (nmix > 0) then
          print 5, (tj(i),i=1,nmix)
    5     format (5x,'tj:',10(f8.5,2x))
        endif
      endif

C --- Save this iteration into mixing file ---
      imix = min(imix+1,mmix)
      if (mmix > 0) then
        if (procid == master) then
          rewind ifi
          write(ifi) imix, nelts
          write(ifi) (((a(i,j,k), i=1,nelts), j=2, imix+1), k=1,2)
          if (ipr) then
            call awrit5(' MCMIXR: write mixfile, unit %i. %i elements.'
     .        //' mmix=%i imix=%i nelts=%i',
     .        ' ',128,i1mach(2),ifi,nelts*imix*2,mmix,imix,nelts)
          endif
        endif
      endif

      end

