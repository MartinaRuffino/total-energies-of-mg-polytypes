      subroutine rxes(s_ctrl,s_site,s_spec,s_lat,nlmax,ndham,nphimx,ndimh,
     .  nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn,wgts)
C- RXES calculation
C--------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:v0
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lmxa lmxb pz orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg jcg indxcg symgr
Cio    Passed to:  *
Ci Inputs:
Ci   nlmax :first dimension of aus; largest augmentation (l+1)^2
Ci   ndham :second dimension of aus, at least as large as ndimh
Ci   nphimx:dimensions aus: max number of radial partial waves per l
Ci   ndimh :number of eigenvalues
Ci   nq    :number of k-points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ef    :Fermi energy
Ci   evl   :energy bands at the nq k-points
Ci   aus   :values and slopes of eigenstates at MT sphere surfaces
Ci          (makusq)
Ci   nsite,isite,iclsl,iclsn see suclst
Ci    N, WW : Methfessel-Paxton order and broadening parameters;
Ci    tol  :allowed error in DOS due to truncating the gaussian,
Ci          if negative on entry, range is set to -tol*W
Ci    wgts :eigenvalue weights
Ci         :NB: only abs(wtkp) is used; see bzmesh.f
Ci    nfilo(2),nfiup(2) lowest filled band and highest filled band for
Ci                      each spin to be included in RXES calc
Ci    nemlo(2), nemup(2) lowest empty band and highest empty band
Ci    npts, number of energy points for mesh
Ci    emin,emax min and max energy for mesh
Co Outputs:
Co Makes RXES spectrum and outputs in jdos style
Cr Remarks
Cr calculates sum_n sum_n' sum_k
Cr            |m1(c,n',pol)m2(c,n,pol)*/(e(n',k-ef-exas-igam)|**2
Cr            delta(w-(en'k-enk))
Cr this is a jdos like quantity with different matrix elements
Cr the jdos is calculated by sampling as in mkjdos
Cr the core valence matrix elements are calculated as in vcdmel
Cr needs inputs from stdin: emin,emax,npts for w mesh
Cr                          exas: XAS energy relative to ef
Cr                          gam: broadening factor
Cr all of this is done for core level  and sites where core-level exists
Cr as in CLS calculation
Cr using --rxes assumes you also set OPTICS params and --cls
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   19 Jun 09 Written by WL adapted from vcdmel and mkjdos
C---------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer npts,nlmax,ndham,nphimx,ndimh,nq,nsp,nspc,nsite,n
      integer isite(nsite),iclsl(nsite),iclsn(nsite)
      double precision ef,
     .                 evl(ndham,nsp,nq),wgts(nq),tol,ww
      double precision exas(0:20),gamma
      double complex aus(nlmax,ndham,nphimx,nsp,nsite,nq)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:)
      real(8), allocatable :: rwgt(:)
      real(8), allocatable :: ul(:)
      real(8), allocatable :: sl(:)
      real(8), allocatable :: gz(:)
      real(8), allocatable :: ruu(:)
      real(8), allocatable :: rus(:)
      real(8), allocatable :: rss(:)
      real(8), allocatable :: me(:)
      real(8), allocatable :: g(:)
      real(8), allocatable :: dos(:)
C ... Local parameters
      integer n0,lmxax,nglob
      parameter (n0=10,lmxax=10)
      integer isp,ib,is,lcls,ncls,i,j,k,nr,lmxa,iq,nlma,i1mach,intopt
      integer lh(10)
      double precision pnu(n0,2),pnz(n0,2),a,rmt,z,xx,s,rsml(n0),ehl(n0)
      double precision emin,emax,step,range,e,d,x,wt,test
      double precision emesh,value
      integer ibn,jbn,jbnlo,mrange,iprint,lgunit,stdo
      integer meshpt,mesh1,mesh2,nfilo(2),nfiup(2),nemlo(2),nemup(2)
      integer npol,ipol,jpol,nexas,ixas
      double complex rxm1(3),rxm2(3),denom
      double precision rxmel,rxm1r(3),rxm1i(3),rxm2r(3),rxm2i(3),
     .  sumg1(3),sumg2(3),rxm1gr(3),rxm1gi(3),rxm2gr(3),rxm2gi(3)
      integer igr,nsgrp
      double precision ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite)
      character clabl*8
C ... External calls
      external adddos,awrit3,awrit8,delstp,dpzero,fclose,getcvmel,grpop,
     .         info0,iorxes,makusp,pvcdm1,pvcdm3,radmwt,rxi,skpblb,tcn,
     .         uspecb

      print *,'RXES nsp=',nsp, 'ef=',ef
      tol = -6d0
      n = 0
      call tcn('rxes')
      call dpzero(rsml,n0)
      call dpzero(ehl,n0)
      if (nspc == 2) call rx('rxes not implemented for nspc=2')

C ... Input from stdi
C      print *, 'give xas e1,e2,de above ef, broadening gamma,ww (in eV)'
C      read(*,*) exasmi,exasma,dexas,gamma,ww
C      nexas=(exasma-exasmi)/dexas
C      print *, 'nexas=',nexas
C      do  ixas = 0, nexas
C      exas(ixas)=(exasmi+ixas*dexas)/13.6057d0
C      enddo
      print *, ' number of xas energies, gamma, ww'
      read (*,*) nexas, gamma,ww
      read (*,*) (exas(ixas),ixas=0,nexas-1)
      print *, ' exas=', (exas(ixas),ixas=0,nexas-1)
      do  ixas = 0, nexas-1
         exas(ixas) = exas(ixas)/13.6057d0
      enddo
c      print *, ' exas=', (exas(ixas),ixas=0,nexas-1)
      gamma = gamma/13.6057d0
      ww = ww/13.6057d0
      print *, 'give nfilo,nfilup for each spin'
      print *, 'followed by nemlo,nemup for each spin'
      read (*,*) (nfilo(isp),nfiup(isp), isp=1,nsp)
      read (*,*) (nemlo(isp),nemup(isp), isp=1,nsp)
      print *, ' give npts, emin,emax (in eV)'
      read(*,*) npts,emin,emax
      emin = emin/13.6057d0
      emax = emax/13.6057d0
C --- Preliminaries for energy mesh and gaussian broadening
      npol = 3
C      print *, 'RXES: npol,npts,nsp:', npol,npts,nsp
      k = npol*npol*npts*nsp*(nexas+1)
      allocate(dos(k)); call dpzero(dos,k)
      stdo = lgunit(1)
      step = (emax - emin) / (npts - 1)
      if (emin < 0) call info0(10,1,0,
     .  ' RXES: (warning) emin<0 for RXES')
C --- Set tolerance and/or range of gaussians ---
      if ( tol > 0d0 ) then
        do  i = 0, npts-1
          x = i * step / ww
          call delstp(0,x,test,s,xx)
          if ( test < tol ) then
            mrange = i + 1
            goto 5
          endif
        enddo
        call info0(20,1,0,' *** Warning : tol too small in MKJDOS')
    5   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * ww
        mrange = range / ( 2 * step )
        call delstp(0,-tol/2,test,s,xx)
      endif

      if (iprint() > 30) call awrit8(
     .  ' RXES:  ef=%1;6d  N.W=%1;6d;  emin=%d  emax=%d  %i bins'//
     .  '%N%10frange of gaussians=%1;4d W (%i bins)  Est DOS err=%1;6d',
     .  ' ',200,stdo,ef,n+ww,emin,emax,npts-1,range/ww,2*mrange,
     .  test)

C reserving space for wave functions, matrix elements etc.
      nsgrp = s_lat%nsgrp
      print *, 'nsite=',nsite
      do  i = 1, nsite
        ib = isite(i)
        ncls = iclsn(i)
        lcls = iclsl(i)
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
C       p_v0 => s_site(ib)%v0
        pnz = s_site(ib)%pz
        clabl = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        z = s_spec(is)%z
        lmxa = s_spec(is)%lmxa
        if (lmxa > lmxax) call rxi('rxes needs lmxax ',lmxa)
        if (lmxa == -1) cycle
        allocate(rofi(nr))
        allocate(rwgt(nr))
C       call radmsh(rmt,a,nr,rofi)
        intopt = 10*nglob('lrquad')
        call radmwt(intopt,rmt,a,nr,rofi,rwgt)
C   --- Augmented wave functions u,s
        k = nr*(lmxa+1)*nsp
        allocate(ul(k),sl(k),gz(k),ruu(k*2),rus(k*2),rss(k*2))
        call uspecb(0,4,s_spec,is,is,lh,rsml,ehl,j)
        call makusp(n0,z,nsp,1,rmt,lmxa,s_site(ib)%v0,a,nr,xx,xx,
     .    pnu,pnz,rsml,ehl,ul,sl,gz,ruu,rus,rss)
C   --- Matrix elements of u,s with core
        call skpblb(clabl,8,j)
        call awrit3('%N RXES atom %i ('//clabl(1:j+1)//') n=%i l=%i',
     .              ' ',120,i1mach(2),ib,ncls,lcls)
        allocate(g(nr*2))
        call pvcdm1(ncls,lcls,g,z,lmxa,s_site(ib)%v0,a,nr,rofi,
     .    rwgt,ul,sl,nsp,lmxax,ume(0,1,i),sme(0,1,i))
        deallocate(rofi,rwgt,g)
        deallocate(ul,sl,gz,ruu,rus,rss)
      enddo

C --- For each qp, make <nk|x,y,z|core> at each site
C     array containing core-valence matrix elements
      allocate(me(3*ndimh*nsite*2))
      do  iq = 1, nq
        do  isp = 1, nsp
          call dpzero(me,3*ndimh*nsite*2)
          do  i = 1, nsite
            lcls = iclsl(i)
            ib = isite(i)
            is = s_site(ib)%spec
            lmxa = s_spec(is)%lmxa
            nlma = (lmxa+1)**2
            if (lmxa > -1) then
              call pvcdm3(i,nsite,ndham,ndimh,nlma,nlmax,
     .          aus(1,1,1,isp,i,iq),ume(0,isp,i),sme(0,isp,i),
     .          lcls,s_lat%cg,s_lat%jcg,s_lat%indxcg,me)
            endif
          enddo
c          print *,'me', (me(j),j=1,3*ndimh*nsite*2)

C now make spectrum contribution from this iq,isp at this site.
          do  i = 1, nsite
c          print *, 'filled bands for spin:', isp,nfilo(isp),nfiup(isp)
          do  ibn = nfilo(isp), nfiup(isp)
           if(evl(ibn,isp,iq) > ef) cycle
           jbnlo = max0(ibn+1,nemlo(isp))
c           print *, 'empty bands:', jbnlo,nemup(isp)
           do  jbn = jbnlo, nemup(isp)
             if (evl(jbn,isp,iq) < ef) cycle
             e = evl(jbn,isp,iq)-evl(ibn,isp,iq)
c              print *, 'eigenvalue diff:', ibn,jbn,e,
c     .        evl(jbn,isp,iq),evl(ibn,isp,iq)
              if (e > emax+step+range/2) cycle
              meshpt = (e - emin) / step
              mesh1 = max(meshpt-mrange,0)
              mesh2 = min(meshpt+mrange,npts-1)
C extract wanted matrix elements from array me
              do  jpol = 1, 3
              call getcvmel(i,ibn,jpol,nsite,ndimh,me,rxm1(jpol))
              call getcvmel(i,jbn,jpol,nsite,ndimh,me,rxm2(jpol))
c              print *, 'rxm1', jpol, rxm1(jpol)
              rxm1r(jpol) = dble(rxm1(jpol))
              rxm1i(jpol) = dimag(rxm1(jpol))
              rxm2r(jpol) = dble(rxm2(jpol))
              rxm2i(jpol)=-dimag(rxm2(jpol))
              sumg1(jpol) = 0d0
              sumg2(jpol) = 0d0
              enddo
c              print *,' xas matrix element before sym',
c     .        (rxm1r(jpol),rxm1i(jpol),jpol=1,3)
c              print *, ' xes matrix element before sym',
c     .        (rxm2r(jpol),rxm2i(jpol),jpol=1,3)
C symmetrize matrix elements
c              print *,' nsgrp=', nsgrp
              do  igr = 1, nsgrp
                 call grpop(rxm1r,rxm1gr,s_lat%symgr,igr)
                 call grpop(rxm1i,rxm1gi,s_lat%symgr,igr)
                 call grpop(rxm2r,rxm2gr,s_lat%symgr,igr)
                 call grpop(rxm2i,rxm2gi,s_lat%symgr,igr)
                 do  jpol = 1, 3
               sumg1(jpol) = sumg1(jpol)+rxm1gr(jpol)**2+rxm1gi(jpol)**2
               sumg2(jpol) = sumg2(jpol)+rxm2gr(jpol)**2+rxm2gi(jpol)**2
                 enddo
              enddo
              do  jpol = 1, 3
              sumg1(jpol) = sumg1(jpol)/nsgrp
              sumg2(jpol) = sumg2(jpol)/nsgrp
              enddo
c              print *,' xas matrix element after sym',
c     .        (sumg1(jpol),jpol=1,3)
c              print *, ' xes matrix element after sym',
c     .        (sumg2(jpol),jpol=1,3)
              do  ixas = 0, nexas-1
              denom = dcmplx(evl(jbn,isp,iq)-ef-exas(ixas),-gamma*0.5d0)
c              print *, 'resonance factor:', 1d0/(denom*dconjg(denom))
              do  ipol = 1, 3
                 do  jpol = 1, 3
              rxmel = sumg1(ipol)*sumg2(jpol)/abs(denom)**2
              wt = abs(wgts(iq))/nsp
c              print *, 'iq,wt', iq,wt
              wt = wt*rxmel
c              print *, 'i,jpol=',ipol,jpol,'matrix elements:',
c     .        sumg1(ipol),sumg2(jpol),denom,rxmel
              do  meshpt = mesh1, mesh2
                  emesh = emin + meshpt * step
c                  if(exas(ixas)-emesh > 0d0) cycle
                  x = (emesh - e) / ww
                  call delstp(n,x,d,s,xx)
                  value = wt*d/ww
                  call adddos(npts,meshpt,ipol,jpol,isp,nexas,ixas,
     .                        value,dos)
                enddo
                enddo
                        enddo
                      enddo
C   4         continue
                enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(me)

C ... Output of dos to file
      call iorxes(dos,npts,emin,emax,nexas,exas,nsp)

      deallocate(dos)
      call tcx('rxes')
      end
      subroutine getcvmel(isite,iband,jpol,nsite,ndimh,s,rxm)
C- select a core-valence matrix element out of the array
      implicit none
      integer isite,iband,jpol,nsite,ndimh
      double complex rxm
      double precision s(3,nsite,ndimh,2)
      rxm = dcmplx(s(jpol,isite,iband,1),s(jpol,isite,iband,2))
      end
      subroutine adddos(npts,meshpt,ipol,jpol,isp,nexas,ixas,v,dos)
C- Adds a value to the dos array
      implicit none
      integer npts,meshpt,ipol,jpol,isp,ixas,nexas
      double precision dos(0:npts-1,3,3,0:nexas,2),v
      dos(meshpt,ipol,jpol,ixas,isp) = dos(meshpt,ipol,jpol,ixas,isp)+v
      end
      subroutine iorxes(dos,npts,emin,emax,nexas,exas,nsp)
C- output of rxes on file
      implicit none
      integer npts,nsp,fopn,ixas,nexas
      double precision dos(0:npts-1,3,3,0:nexas,2),emin,emax
C local variables
      double precision de,e,exas(0:nexas),sum
      integer isp,ifio,ipts,ipol,jpol
      ifio = fopn('RXES')
      rewind ifio
      de=(emax-emin)/(npts-1)
      do  ipol = 1, 3
      do  jpol = 1, 3
         write(ifio,'(''ipol='',i3,2x,''jpol='',i3)') ipol,jpol
      do  ixas = 0, nexas-1
         write(ifio,'(''exas='',f12.5)') exas(ixas)*13.6057d0
        do  ipts = 1, npts-1
        e = emin+de*(ipts-1)
        e = exas(ixas)-e
        if (nsp == 2)  then
        sum = dos(ipts,ipol,jpol,ixas,1)+
     .     dos(ipts,ipol,jpol,ixas,2)
        write(ifio,550) e*13.6057d0,(dos(ipts,ipol,jpol,ixas,isp),
     .  isp=1,nsp),sum
        else
           write(ifio,551) e*13.6057d0,dos(ipts,ipol,jpol,ixas,1)
        endif
        enddo
      enddo
      enddo
      enddo
  550 format(f12.5,3(e14.6))
  551 format(f12.5,e14.6)
       call fclose(ifio)
       return
       end
      subroutine pvcdm3(isite,nsite,ndham,ndimh,nlma,nlmax,aus,ume,sme,
     .                  lcls,cg,jcg,indxcg,s)
C- Kernel called by vcmdel
C ----------------------------------------------------------------------
Ci Inputs
Ci   isite :
Ci   nsite :
Ci   ndimh :
Ci   nlma  :
Ci   nlmax :
Ci   aus   :
Ci   ume   :
Ci   sme   :
Ci   lcls  :
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Co Outputs
Co   s     :Matrix elements
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer isite,lcls,ndham,ndimh,nlma,nlmax,nsite,indxcg(*),jcg(*)
      double precision cg(*),ume(0:1),sme(0:1),s(3,nsite,ndimh,2)
      double complex aus(nlmax,ndham,2)
C ... Local parameters
      integer kk(4),mlm,lm,ll,klm,ii,indx,icg1,icg2,icg,llm,ib
      double complex cxx
C     Transposes (y,z,x) to (x,y,z)
      data kk /0,2,3,1/

C      write (*,100)
C ... Loop over lm of (u,s)
      do  mlm = 1, nlma
        lm = ll(mlm)
C       Selection rule would be handled by CG anyway:
        if (lm == lcls-1 .or. lm == lcls+1) then
C     ... Loop over m index of r
          do  klm = 2, 4
            ii = max0(mlm,klm)
            indx = (ii*(ii-1))/2 + min0(mlm,klm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do  icg = icg1, icg2
C             lm of core
              llm  = jcg(icg)
              if (ll(llm) == lcls) then
                do  ib = 1, ndimh
                  cxx =  cg(icg)*
     .              (dconjg(aus(mlm,ib,1))*ume(lm) +
     .               dconjg(aus(mlm,ib,2))*sme(lm))
C                 write (*,200) kk(klm),mlm,klm,llm,aus(mlm,ib,1),
C     .              aus(mlm,ib,2),cg(icg),cxx
                  s(kk(klm),isite,ib,1) = s(kk(klm),isite,ib,1)
     .                                  + dble(cxx)
                  s(kk(klm),isite,ib,2) = s(kk(klm),isite,ib,2)
     .                                  + dimag(cxx)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

C  100 format(2x,'n',2x,'L(val)',1x,'L(r)',1x,'L(core)',11x,'au',25x,
C     .        'as',14x,'CG',20x,'cxx')
C  200 format(1x,i3,3i5,2x,'(',g12.5,',',g12.5,')',
C     .  '(',g12.5,',',g12.5,')',g12.5,'(',g12.5,',',g12.5,')')

      end

