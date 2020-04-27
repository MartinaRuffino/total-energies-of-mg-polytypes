      subroutine mkrout(sopts,s_site,s_spec,s_lat,s_ham,nbas,nsp,nlmto,
     .  lekkl,s_rout,hab,sab,qbyl,hbyl,lrout)
C- Assembles local output densities out of the qkkl, and core states
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu v0 v1 pz qkkl qhkl qhhl eqkkl eqhkl eqhhl
Ci                 class pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  symrat
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rsma lmxa lmxl kmxt a nr rmt lmxb stc lfoca rfoca qc
Ci                 z ctail etail p pz rg name orbp coreh coreq
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:rhoc
Cio    Passed to:  corprm uspecb gtpcor symrat
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat qlat nsgrp afmt npgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg jcg indxcg istab symgr ag
Cio    Passed to:  symrat
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  eterms
Co     Stored:     eterms
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Cio  s_rout
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  symrat
Ci Inputs
Ci   sopts:string specifying partial densities to write
Ci        :If sopts(1:8)=='--wrhoat' mkrout acts in a special mode:
Ci        :atom densities are written to disk.
Ci        :Optional modifiers: --wrhoat[:l=#][:lx=#]
Ci        :  l=#  accumulates density from psi_l psi_l products only
Ci        :  lx=# accumulates density excluding psi_l psi_l products
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmto :dimension of lower block of LMTO part of hamiltonian
Cl   lekkl :1 Find band CG from eqkkl/qkkl
Ci   hab   :hamiltonian matrix elements of radial wave functions at each
Ci         :site.  See potpus for their definition.
Ci   sab   :overlap matrix elements of radial wave functions at each
Ci         :site.  See potpus for their definition.
Ci         :hab and sab are used here to find band cg
Ci   lrout :0 calculate core states part only
Co Outputs
Co   s_rout:vector of offsets containing site density
Co         :lfoc=0 => rhoc is updated (call to getcor)
Co         :else      taken from spec%rhoc
Co   sham->eterms various integrals for the total energy are stored:
Co         :(8)  sumec = sum of foca=0 core eigenvalues
Co         :(9)  sumtc = sum of all core kinetic energies
Co         :(12) sumt0 = sum of frozen core kinetic energies
Co   qbyl  :l-decomposed charge
Co   hbyl  :l-decomposed eigenvalue sum
Cr Remarks
Cr   u and s are linear combinations of and phi,phidot defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1
Cu Updates
Cu   03 Dec 18 symmetrisation of charges for printout (bugfix,JJackson)
Cu   19 Dec 13 First cut at nocollinear case (up to mkrou5)
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   22 Nov 12 Replace qkkl with structures
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   30 Jul 08 (T. Kotani) Use ekkl to set hbyl
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   28 Aug 01 Extended to local orbitals.
Cu   18 Jun 00 spin polarized
Cu   20 May 00 adapted from nfp mk_dens2
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer n0,nab,nbas,nsp,lekkl,nlmto,lrout
      parameter (n0=10,nab=9)
      double precision qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas),
     .  sab(nab,n0,nsp,nbas),hab(nab,n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_rhat)::  s_rout(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:),rwgt(:)
      real(8), allocatable :: fh(:),xh(:),vh(:),dh(:)
      real(8), allocatable :: fp(:),xp(:),vp(:),dp(:)
      real(8), allocatable :: ul(:),sl(:),gz(:)
      real(8), allocatable :: ruu(:),rus(:),rss(:)
      real(8), allocatable :: chh(:,:),chk(:,:),ckk(:,:)
      real(8), allocatable :: dmatl(:,:)
      real(8), pointer :: p_v0(:,:),p_v1(:,:),qkk(:,:),qhk(:,:),qhh(:,:)
C ... Local parameters
      character dc*1
      integer ib,intopt,ipr,is,j1,j2,k,kcor,kmax,lcor,lfoc,
     .  linrh1,linrh2,lmxa,lmxh,lmxl,lxrh1,lxrh2,ncore,nlma,nlmh,nlml,
     .  nlml1,nr,nspc,stdo
C     integer iv(2)
      double precision a,ceh,pi,res,rfoc,rmt,rsma,smec,smtc,stc0,
     .  sum1,sum2,sumec,am1(0:3),am2(0:3),sumt0,sumtc,xx,y0,z
      double precision qcor(2),rsml(n0),ehl(n0),x2(1)
      procedure(integer) :: iprint,a2vec,lgunit
C     character strn*40
C ... For basis
      integer nkap0,nkapi,nkape,nkaph,nglob
      parameter (nkap0=4)
      integer lh(nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer norb,ntab(n0*nkap0),blks(n0*nkap0),ipc(n0)
      double precision pnu(n0,2),pnz(n0,2),qcorg,qcorh,qsc,cofg,cofh
C     this ordering must match sham->eterms; see uham
      double precision eterms(22)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(12),sumt0)
C ... External calls
      external corprm,dinv22,dpcopy,dpzero,dsumdf,fradhd,fradpk,getcor,
     .         gtbsl1,gtpcor,info5,makusp,mkcfus,mkrou1,mkrou2,mkrou3,
     .         mkrou4,mkrou5,mkrou6,nwordg,orbl,poppr,pshpr,radmsh,
     .         radsum,radwgt,rx0,rxs2,symrat,tcn,tcx,uspecb,wrhomt

C ... Setup
      stdo = lgunit(1)
      ipr  = iprint()
      pi   = 4d0*datan(1d0)
      y0   = 1d0/dsqrt(4d0*pi)
      nkaph = nglob('nkaph')
      nspc = nglob('nspc')
      linrh1 = -1
      lxrh1 = -1
      linrh2 = -1
      lxrh2 = -1
      call tcn('mkrout')
      eterms = s_ham%eterms
      intopt = 10*nglob('lrquad')

C ... Dump atom files, partial density
      if (len(sopts) >= 8) then  ! So that the routine can be compiled with -check
      if (sopts(1:8) == '--wrhoat') then
        j2 = 9
        dc = sopts(j2:j2)
        if (dc /= ' ') then
C   ... Return here to resume parsing for arguments
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then
          elseif (sopts(j1:j1+1) == 'l=')  then
            k = j1+1
            k = a2vec(sopts,j2,k,2,dc//' ',2,3,1,lh,linrh1)
            if (k /= 1) call rxs2('mkrout: failed to parse:  ',
     .        '--wrhoat'//dc//'...',sopts(j1:j2))
            if (linrh2 == -1) linrh2 = linrh1
          elseif (sopts(j1:j1+2) == 'lx=')  then
            k = j1+2
            k = a2vec(sopts,j2,k,2,dc//' ',2,3,1,lh,lxrh1)
            if (k /= 1) call rxs2('mkrout: failed to parse:  ',
     .        '--wrhoat'//dc//'...',sopts(j1:j2))
            if (lxrh2 == -1) lxrh2 = lxrh1
          endif
          goto 10
        endif
      endif
      if (linrh1 >= 0) lxrh1 = -1
      if (linrh2 >= 0) lxrh2 = -1
      call info5(30,1,0,' mkrout:  write atom files with:'//
     .  '%?#n>=0#%-1j linc1=%i##'//
     .  '%?#n>=0#%-1j linc2=%i##'//
     .  '%?#n>=0#%-1j lxc1=%i##'//
     .  '%?#n>=0#%-1j lxc2=%i##',
     .  linrh1,linrh2,lxrh1,lxrh2,0)
      endif
      endif

C --- Start loop over atoms ---
      sumtc = 0d0
      sumec = 0d0
      sumt0 = 0d0
      if (ipr >= 30 .and. lrout > 0) then
        if (nsp == 1) write(stdo,201)
        if (nsp == 2) write(stdo,202)
  201   format(/' mkrout:  Qtrue      sm,loc       local')
  202   format(/' mkrout:  Qtrue      sm,loc       local',8x,
     .    'true mm   smooth mm    local mm')
      endif
      do  ib = 1, nbas
        is = s_site(ib)%spec
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl
        kmax = s_spec(is)%kmxt
        if (lmxa == -1) cycle

        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        lmxh = s_spec(is)%lmxb
        stc0 = s_spec(is)%stc
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)

        call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
        call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkape)
        call orbl(ib,0,nlmto,s_ham%iprmb,norb,ltab,ktab,xx,offl,xx)
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        p_v0 => s_site(ib)%v0
        p_v1 => s_site(ib)%v1
        pnz = s_site(ib)%pz
        call gtpcor(s_spec,is,kcor,lcor,qcor)
        nlml = (lmxl+1)**2
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        allocate(rofi(nr),rwgt(nr))
        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)
        if (lrout /= 0) then
          qkk => s_site(ib)%qkkl
          qhk => s_site(ib)%qhkl
          qhh => s_site(ib)%qhhl
C         Take average of qhk(2,1)+qhk(2,1)* since make (1,2) block only
          if (nspc == 2) then
            k = size(qhk)/6
            call dsumdf(k*2,0.5d0,qhk(1,3),0,1,qhk(1,5),0,1)

C            k = size(qhk)/6
C            call prmx('qhk(1)',qhk(1,1),k,k,1)
C            call prmx('qhk(2)',qhk(1,2),k,k,1)
C            call prmx('qhk(3)',qhk(1,3),k,k,1)
C            call prmx('qhk(4)',qhk(1,4),k,k,1)
C            call prmx('qhk(5)',qhk(1,5),k,k,1)
C            call prmx('qhk(6)',qhk(1,6),k,k,1)
C           stop
          endif

          call dpzero(s_rout(ib)%rho1,nr*nlml*nsp*nspc)
          call dpzero(s_rout(ib)%rho2,nr*nlml*nsp*nspc)

C     --- Assemble rho1 = products of augmented functions ---
C     ... Set up all radial head and tail envelope functions, and their bc's
          k =  (lmxh+1)*nkaph
          allocate(fh(nr*k),xh(nr*k),vh(k),dh(k))
          call fradhd(nkaph,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh)
          k =  (lmxa+1)*(kmax+1)
          allocate(fp(nr*k),xp(nr*k),vp(k),dp(k))
          call fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)

C     ... Augmented wave functions
          k = nr*(lmxa+1)*nsp
          allocate(ul(k),sl(k),gz(k*2),ruu(k*3),rus(k*3),rss(k*3))
          call uspecb(0,4,s_spec,is,is,lh,rsml,ehl,k)
          call makusp(n0,z,nsp,nspc,rmt,lmxa,p_v0,a,nr,xx,xx,pnu,pnz,
     .      rsml,ehl,ul,sl,gz,ruu,rus,rss)

C     ... Contracted density matrix as coffs to products of (u,s,gz)
          k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml*nsp
          allocate(chh(k,nsp*nspc)); call dpzero(chh,k*nsp*nspc)
          allocate(chk(k,nsp*nspc)); call dpzero(chk,k*nsp*nspc)
          allocate(ckk(k,nsp*nspc)); call dpzero(ckk,k*nsp*nspc)
          k = (lmxa+1)**2*nlml*9
          allocate(dmatl(k,nsp*nspc)); call dpzero(dmatl,k*nsp*nspc)

          call mkrou1(nsp,nspc,nlmh,nlma,nlml,kmax,s_lat%cg,
     .      s_lat%jcg,s_lat%indxcg,nkaph,nkapi,norb,ltab,ktab,blks,
     .      qhh,qhk,qkk,vh,dh,vp,dp,chh,chk,ckk,dmatl)

C         Restore original (1,2) and (2,1) blocks of qhk
          if (nspc == 2) then
            k = size(qhk)/6
            call dsumdf(k*2,1d0,qhk(1,3),0,1,qhk(1,5),0,1)
          endif

C          iv = shape(qhh)
C          call prmx('qhh',qhh,iv(1),iv(1),iv(2))
C          call prmx('dmatl',dmatl,(lmxa+1)**2,(lmxa+1)**2,nlml*9*nsp*nspc)
C          k = (lmxa+1)**2*nlml*9
C          call prmx('dmatl',dmatl,k,k,nsp*nspc)

C     ... True local density for this sphere
          call mkrou2(nsp,nspc,lmxa,nlml,linrh1,lxrh1,linrh2,lxrh2,pnz,
     .      dmatl,nr,ul,sl,gz,ruu,rus,rss,s_rout(ib)%rho1)

C     --- Assemble rho2 = unaugmented products Pkl*Pk'l' ---
C         H H product.  Include local orbitals w/ envelopes (nkapi->nkaph)
          call mkrou5(nsp*nspc,nr,nlml,nkaph,nkaph,fh,lmxh,
     .      nkaph,nkaph,fh,lmxh,chh,s_rout(ib)%rho2)
C         H Pkl product
          call mkrou5(nsp*nspc,nr,nlml,nkaph,nkaph,fh,lmxh,
     .      kmax+1,kmax+1,fp,lmxa,chk,s_rout(ib)%rho2)
C         Pkl Pkl product
          call mkrou5(nsp*nspc,nr,nlml,kmax+1,kmax+1,fp,lmxa,
     .      kmax+1,kmax+1,fp,lmxa,ckk,s_rout(ib)%rho2)
C         call prmx('rho2',s_rout(ib)%rho2,nr*nlml,nr*nlml,nsp*nspc)
C         call prrmsh('sm rout',rofi,s_rout(ib)%rho2,nr,nr,nlml*nsp)

C     ... Noncollinear development not finished below this point
C         Need rotate to a local quantization axis.

C     ... Site charge and eigenvalue sum decomposed by l
          call mkrou3(2,lmxa,nlml,nsp,pnz,dmatl,
     .      hab(1,1,1,ib),sab(1,1,1,ib),qbyl(1,1,ib),hbyl(1,1,ib))
C         Remake hbyl from energy-weighted density matrix
          if (lekkl /= 0) then
            qkk => s_site(ib)%eqkkl
            qhk => s_site(ib)%eqhkl
            qhh => s_site(ib)%eqhhl
            nlml1 = 1
            call dpzero(dmatl,(lmxa+1)**2*nlml1*nsp*9)
            k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml1*nsp
            call dpzero(chh,k)
            call dpzero(chk,k)
            call dpzero(ckk,k)
            call mkrou1(nsp,nspc,nlmh,nlma,nlml1,kmax,s_lat%cg,
     .        s_lat%jcg,s_lat%indxcg,nkaph,nkapi,norb,ltab,ktab,blks,
     .        qhh,qhk,qkk,vh,dh,vp,dp,chh,chk,ckk,dmatl)
            call mkrou3(1,lmxa,nlml1,nsp,pnz,dmatl,x2,
     .        sab(1,1,1,ib),hbyl(1,1,ib),x2)
          endif

C     ... Print charges for information
CJJ,30.11.2018,bug: sum1,sum2,am1,am2 differ with and without symmetry
C          call radsum(nr,nr,nlml,nsp,rwgt,s_rout(ib)%rho1,sum1)
C          call radsum(nr,nr,nlml,nsp,rwgt,s_rout(ib)%rho2,sum2)
C          sum1 = sum1/y0
C          sum2 = sum2/y0
C          if (nsp == 2) then
C            call magvec(nr,nlml,nsp,nspc,rwgt,s_rout(ib)%rho1,am1)
C            call magvec(nr,nlml,nsp,nspc,rwgt,s_rout(ib)%rho2,am2)
C            if (ipr >= 30) write(stdo,200) ib,sum1,sum2,sum1-sum2,
C     .                     am1(0),am2(0),am1(0)-am2(0)
C            if (nspc == 2) call info(30,0,0,' mag true %3;10,6D '//
C     .          ' smooth%3;10,6D',am1(1),am2(1))
C          else
C            if (ipr >= 30) write(stdo,200) ib,sum1,sum2,sum1-sum2
C          endif
C  200     format(i4,3f12.6,2x,3f12.6)
          deallocate(dmatl)

C     ... Contribution to mag outsite rmt: extrapolate tail to infinity
C          call mkrou6(rofi,s_rout(ib)%rho1,nr,nlml,nsp,xx,xx,res)
C          if (ipr >= 30 .and. res /= 0) then
C            write(stdo,211) res,res+am1(0)
C  211       format(7x,'contr. to mm extrapolated for r>rmt:',f11.6,
C     .        ' est. true mm =',f9.6)
C          endif

        endif

C   --- Make new core density and core eigenvalue sum ---
        if (lfoc == 0) then
          call pshpr(ipr-11)
          call getcor(0,z,a,pnu,pnz,nr,lmxa,rofi,p_v1,kcor,lcor,
     .      qcor,smec,smtc,s_rout(ib)%rhoc,ncore,0d0,0d0)
          call poppr
C     ... here if v0 used to make core... then effectively frozen
C         call p1kden(nr,rwgt,s_rout(ib)%rhoc,p_v1,p_v0,sum)
C         write(stdo,996) smec,sum,smec+sum
C 996     format('smec,sum,smec+sum=',3f14.6)
C         smec = smec+sum

          sumtc = sumtc + smtc
          sumec = sumec + smec
        else
          if (ipr >= 41) write(stdo,288) stc0
  288     format(' foca..  use smtc = ',f14.6)
          sumtc = sumtc + stc0
          sumt0 = sumt0 + stc0
          if (associated(s_spec(is)%rhoc)) then
            call dpcopy(s_spec(is)%rhoc,s_rout(ib)%rhoc,1,nr*nsp,1d0)
          endif
        endif

C        if (ipr >= 80) then
C          call awrit1(' output local density, ib=%i',strn,80,0,ib)
C          call shorho(3,strn,nr,nlml,nsp,w(orofi),w(orwgt),
C     .      s_rout(ib)%rho1,s_rout(ib)%rho2)
C        endif

        deallocate(rofi,rwgt)
        deallocate(fh,xh,vh,dh,fp,xp,vp,dp)
        deallocate(ul,sl,gz,ruu,rus,rss)
        deallocate(chh,chk,ckk)

      enddo

CJJ, symmetrise site charges/moments before printing
C fixes symmetry/nosym bug for printed charges/moments
C symrat (file symrho.f) acts on all bases at once (--> outside ib loop)
      if (ipr >= 30 .and. lrout > 0) then
        call sitepack(s_site,1,nbas,'class',1,ipc,x2)
        call symrat(s_site,s_spec,nbas,nsp,nspc,6,
     .  s_lat%plat,s_lat%nsgrp,s_lat%istab,s_lat%symgr,s_lat%ag,ipc,s_rout,qbyl,hbyl,x2)

        do  ib = 1, nbas
          is = s_site(ib)%spec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          lmxl = s_spec(is)%lmxl
          nlml = (lmxl+1)**2
          if (lmxl<0) cycle
          allocate(rofi(nr),rwgt(nr))
          call radmsh(rmt,a,nr,rofi)
          call radwgt(intopt,rmt,a,nr,rwgt)
          call radsum(nr,nr,nlml,nsp,rwgt,s_rout(ib)%rho1,sum1)
          call radsum(nr,nr,nlml,nsp,rwgt,s_rout(ib)%rho2,sum2)
          sum1 = sum1/y0
          sum2 = sum2/y0
          if (nsp == 2) then
            call magvec(nr,nlml,nsp,nspc,rwgt,s_rout(ib)%rho1,am1)
            call magvec(nr,nlml,nsp,nspc,rwgt,s_rout(ib)%rho2,am2)
            if (ipr >= 30) write(stdo,200) ib,sum1,sum2,sum1-sum2,
     .                     am1(0),am2(0),am1(0)-am2(0)
            if (nspc == 2) call info(30,0,0,' mag true %3;10,6D '//
     .          ' smooth%3;10,6D',am1(1),am2(1))
          else
            if (ipr >= 30) write(stdo,200) ib,sum1,sum2,sum1-sum2
          endif
  200     format(i4,3f12.6,2x,3f12.6)

C     ... Contribution to mag outsite rmt: extrapolate tail to infinity
          call mkrou6(rofi,s_rout(ib)%rho1,nr,nlml,nsp,xx,xx,res)
          if (ipr >= 30 .and. res /= 0) then
            write(stdo,211) res,res+am1(0)
  211       format(7x,'contr. to mm extrapolated for r>rmt:',f11.6,
     .          ' est. true mm =',f9.6)
          endif

          deallocate(rofi,rwgt)
        enddo
      endif
CJJ, ends symmetrisation of charges/moments

      if (len(sopts) >= 8) then  ! So that the routine can be compiled with -check
      if (sopts(1:8) == '--wrhoat') then
        call rx('mkrout : fix call to symrat ... need ipc')
C        call symrat(s_site,s_spec,nbas,nsp,nspc,.false.,
C     .    s_lat%plat,s_lat%nsgrp,s_lat%istab,s_lat%symgr,s_lat%ag,
C     .    ipc,s_rout,qbyl,hbyl,xx)
        do  ib = 1, nbas
          if (lmxa == -1) cycle
          is = s_site(ib)%spec
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          lmxl = s_spec(is)%lmxl
          nlml = (lmxl+1)**2
          allocate(rofi(nr))
          call radmsh(rmt,a,nr,rofi)
C         call prrmsh('rout',rofi,s_rout(ib)%rho1,nr,nr,nlml*nsp)
          call pshpr(40)
          call wrhomt(1,'rhoMT.','density',ib,s_rout(ib)%rho1,xx,
     .      rofi,nr,nlml,nsp)
          call poppr
          deallocate(rofi)
        enddo
        call rx0('atom files written to disk')
      endif
      endif

C --- Put sumec,sumtc into etot struct ---
      s_ham%eterms = eterms

      call tcx('mkrout')
      end

      subroutine mkrou1(nsp,nspc,nlmh,nlma,nlml,kmax,cg,jcg,indxcg,
     .  nkaph,nkapi,norb,ltab,ktab,blks,qhh,qhk,qkk,vh,dh,vp,dp,
     .  chh,chk,ckk,dmatl)
C- Contracted density matrix for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmh  :(lmxh+1)**2, where lmxh = L-cutoff in basis
Ci   nlma  :(lmxa+1)**2, where lmxa = L-cutoff in augmentation
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   kmax  :polynomial cutoff in augmentation
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   nkaph :number of types of one l-quantum number in the MTO basis
Ci   nkapi :number of the nkaph functions that are envelope functions
Ci   norb  :number of orbital types for this site; see orbl.f
Ci   ltab  :table of l-quantum numbers for each type
Ci   ktab  :table of energy index for each type
Ci   blks  :blks(iorb) = size of contiguous block of orbitals for
Ci         :orbital type iorb and possibly iorb+1,iorb+2...
Ci   qhh   :head-head component of site part of density matrix (rlocbl.f)
Ci   qhk   :head-tail component of site part of density matrix (rlocbl.f)
Ci   qkk   :tail-tail component of site part of density matrix (rlocbl.f)
Ci   vh    :values of head functions on MT sphere
Ci   dh    :slopes of tail functions on MT sphere
Ci   vp    :values of Pkl  functions on MT sphere
Ci   dp    :slopes of Pkl  functions on MT sphere
Co Outputs
Co   chh   :head-head product function coefficients
Co   chk   :head-tail product function coefficients
Co   ckk   :tail-tail product function coefficients
Co   dmatl :dmatl(l1,l2,mlm,i,j,:,:) holds coefficients to a Y_lm
Co         :expansion of the function products f_i(l1) f_j(l2)
Co         :where f_i,f_j form this table.
Co         :    (uu  us  uz)
Co         :    (su  ss  sz)
Co         :    (zu  zs  zz)
Co         :Spin indices:
Co         :  dmatl(:,:,:,:,:,1) is dmat of (1,1) spin block
Co         :  dmatl(:,:,:,:,:,2) is dmat of (2,2) spin block
Co         :  dmatl(:,:,:,:,:,3) is Re[dmat] of (1,2) spin block
Co         :  dmatl(:,:,:,:,:,4) is Im[dmat] of (1,2) spin block
Cr Remarks
Cr   Transforms density-matrix as generated by addrbl.f into
Cr   contracted density-matrix of wave function products.
Cr   Radial parts of chh,chk,ckk are
Cr     (psi1)_l (psi2)_l'
Cr   Radial parts of dmatl are:
Cr     (u, s, or gz)_l (u, s, or gz)_l'
Cu Updates
Cu   29 Jul 08 Adapted from mkrout.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nsp,nspc,nlmh,nlma,nlml,kmax,norb,nkaph,nkapi
      integer jcg(*),indxcg(*)
      integer ltab(norb),ktab(norb),blks(norb)
      double precision qhh(nkaph,nkaph,nlmh,nlmh,nsp*nspc)
      double precision qhk(nkaph,kmax+1,nlmh,nlma,nsp*nspc)
      double precision qkk(kmax+1,kmax+1,nlma,nlma,nsp*nspc)
C     double precision dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
C     double precision vh(0:lmxh,nkaph),dh(0:lmxh,nkaph),
C    .                 vp(0:lmxa,kmxa+1),dp(0:lmxa,kmxa+1)
C     chh = chh(nkaph,nkaph,0:lmxh,0:lmxh,nlml,nsp,nspc)
C     chk = chk(nkaph,kmax+1,0:lmxh,0:lmxa,nlml,nsp,2*nspc-1)
C     ckk = ckk(kmax+1,kmax+1,0:lmxa,0:lmxa,nlml,nsp,nspc)
      double precision dmatl(*),vh(*),dh(*),vp(*),dp(*),
     .  chh(*),chk(*),ckk(*)

      double precision cg(*)
      integer ll,k,lmxa,lmxh
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer ltba(n0*nkap0),ktba(n0*nkap0),blka(n0*nkap0)

C     lmxl = ll(nlml)
      lmxa = ll(nlma)
      lmxh = ll(nlmh)
      do  k = 0, kmax
        ltba(k+1) = 0
        ktba(k+1) = k+1
        blka(k+1) = nlma
      enddo

C ... Contracted density-matrix chh,chk,ckk from qhh,qhk,qkk
C     H H product
      call mkrou4(nsp*nspc,nlml,cg,jcg,indxcg,
     .  nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .  nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .  qhh,chh)
C     H Pkl product
      call mkrou4(nsp*nspc,nlml,cg,jcg,indxcg,
     .  nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .  kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .  qhk,chk)
C     Pkl Pkl product
      call mkrou4(nsp*nspc,nlml,cg,jcg,indxcg,
     .  kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .  kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .  qkk,ckk)

C ... Contracted density matrix as coffs to products of (u,s,gz)
C     H H product
      call mkcfus(nsp*nspc,lmxa,nlml,
     .  nkaph,nkapi,vh,dh,lmxh,
     .  nkaph,nkapi,vh,dh,lmxh,
     .  chh,dmatl)
C     H Pkl product
      call mkcfus(nsp*nspc,lmxa,nlml,
     .  nkaph,nkapi,vh,dh,lmxh,
     .  kmax+1,kmax+1,vp,dp,lmxa,
     .  chk,dmatl)
C     Pkl Pkl product
      call mkcfus(nsp*nspc,lmxa,nlml,
     .  kmax+1,kmax+1,vp,dp,lmxa,
     .  kmax+1,kmax+1,vp,dp,lmxa,
     .  ckk,dmatl)

      end

      subroutine mkrou2(nsp,nspc,lmxa,nlml,linrh1,lxrh1,linrh2,lxrh2,
     .  pnz,dmatl,nr,ul,sl,gz,ruu,rus,rss,rho)
C- Assemble true site density from product function coefficients
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   lmxa  :augmentation l-cutoff
Ci   nlml  :l-cutoff for charge density
Ci   linrh1:<0 : linrh1 is not used
Ci         :accumulate contributions only from l1=linrh1
Ci   lxrh1 :<0 or linrh1>=0: lxrh1 is not used
Ci         :exclude contributions only from l1=lxrh1
Ci   linrh2:<0 : linrh2 is not used
Ci         :accumulate contributions only from l2=linrh2
Ci   lxrh2 :<0 or linrh2>=0: lxrh2 is not used
Ci         :exclude contributions only from l2=lxrh2
Ci   dmatl :dmatl(l1,l2,mlm,i,j,:,:) holds coefficients to a y_lm
Ci         :expansion of the function products f_i(l1) f_j(l2)
Ci         :where f_i,f_j form this table.
Ci         :    (uu  us  uz)
Ci         :    (su  ss  sz)
Ci         :    (zu  zs  zz)
Ci         :Spin indices:
Ci         :  dmatl(:,:,:,:,:,isp)   = dmat of (isp,isp) spin block
Ci         :  dmatl(:,:,:,:,:,isp+2) = (Re[dmat],Im[dmat]  of (1,2) spin block
Ci         :  dmatl(:,:,:,:,:,isp+4) = (Re[dmat],-Im[dmat] of (2,1) spin block
Ci         :                           See prlcb2 for (-) sign in (2,1) block
Ci   nr    :number of radial mesh points
Ci   ul    :r*radial wave functions; see remarks
Ci   sl    :r*radial wave functions; see remarks
Ci   gz    :r*semicore wave functions; see remarks
Ci   ruu   :l-diagonal (uu) product including small component (makusp)
Ci   rus   :l-diagonal (us) product including small component (makusp)
Ci   rss   :l-diagonal (ss) product including small component (makusp)
Ci         :rxx(:,:,:,1:2) are spin-diagonal products
Ci         :rxx(:,:,:,3)  is a spin off-diagonal products (avg of 12, 21)
Ci         :See makusp.f
Co Outputs
Co   rho   :charge density assembled
Co         :Spin indices:
Co         :  rho(:,:,isp,1,) is rho of (isp,isp) spin block
Co         :  rho(:,:,1,2) is Re[rho] of (1,2) spin block
Co         :  rho(:,:,2,2) is Im[rho] of (1,2) spin block
Cl Local variables
Cr Remarks
Cr   This routine uses linear combinations (u,s) of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are returned as r * u and r * s, respectively.
Cr
Cr   Let phi_z be the w.f. corresponding to pnu_z.
Cr   A local orbital of the first type is defined as follows.
Cr      gz = r * ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
Cr   By construction, gz/r has both value = 0 and slope = 0 at rmax.
Cr   A local orbital of the second type is defined as gz=r*phi_z;
Cr   for r>rmax a smooth Hankel tail (spec'd by ehl,rsml) is attached.
Cr
Cr   Spherical density assembled using ruu,rus,rss which includes the
Cr   small component of the scalar relativistic wave function.
Cu Updates
Cu   28 aug 01 extended to local orbitals.  altered argument list.
C ----------------------------------------------------------------------
      implicit none
C ... passed parameters
      integer nsp,nspc,lmxa,nlml,nr,n0,linrh1,lxrh1,linrh2,lxrh2
      parameter (n0=10)
      double precision rho(nr,nlml,nsp,nspc),pnz(n0,nsp),
     .  dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp*nspc),
     .  ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,2*nsp-1),
     .  sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,2*nsp-1),
     .  gz(2*nr,0:lmxa,nsp),rss(nr,0:lmxa,2,2*nsp-1)
C ... local parameters
      logical lpz1,lpz2
      integer l,i,mlm,l1,l2,isp,lin,lx,ispc,i1,i2,is1,is2,ksp
      double precision xuu,xus,xsu,xss,xuz,xsz,xzu,xzs,xzz
C     double precision top

      call tcn('mkrou2')

C --- Full density as products (u,s) (u,s); no small component ---
      do  ispc = 1, nspc    ! Spin off-diagonal parts with ispc=2
      do  isp = 1, nsp
      do  mlm = 1, nlml
        do  l1 = 0, lmxa
          if (linrh1 >= 0 .and. l1 /= linrh1) then
            cycle
          elseif (lxrh1 >= 0 .and. l1 == lxrh1) then
            cycle
          endif
          do  l2 = 0, lmxa
            if (linrh2 >= 0 .and. l2 /= linrh2) then
              cycle
            elseif (lxrh2 >= 0 .and. l2 == lxrh2) then
              cycle
            endif

            lpz1 = pnz(l1+1,1) /= 0
            lpz2 = pnz(l2+1,1) /= 0

C            top = xuu*xuu + xus*xus + xsu*xsu + xss*xss
C            if (dsqrt(top) > 1d-6 .and. mlm == 16)
C     .        write (6,700) mlm,l1,l2,isp,ispc,xuu,xus,xsu,xss
C  700       format(5i5,4f14.8)

C           Spin indices for left and right partial waves:
C           i1,i2 = spin indices to radial functions for current spin block
C           is1,is2 = indices to spin block of rho we are writing to
C           Note: contr to (2,1) block from qhp added into rho(1,2) block
C           dmatl(:,:,:,:,:,isp)   : dmat of (isp,isp) spin block (always real)
C           dmatl(:,:,:,:,:,isp+2) : (Re[dmat], Im[dmat]) of (1,2) block
C           rho(:,:,isp,1) : spin diagonal block of rho
C           rho(:,:,1,2)   : Re[rho] of (1,2) spin block
C           rho(:,:,2,2)   : Im[rho] of (1,2) spin block
            if (ispc == 1) then
              i1 = isp; i2 = isp ! For spin diagonal block
              is1 = isp; is2 = 1 ! indices to spin diagonal part of rho(isp)
              ksp = isp
            elseif (ispc == 2) then
              i1 = 1; i2 = 2    ! For (1,2) block
              is1 = isp; is2 = 2! Indices to spin 12 block of rho:
              ksp = isp+2       ! Re rho(1,2) if isp=1, Im rho(1,2) if isp=2
            endif
            xuu = dmatl(l1,l2,mlm,1,1,ksp) ! uu
            xus = dmatl(l1,l2,mlm,1,2,ksp) ! us
            xsu = dmatl(l1,l2,mlm,2,1,ksp) ! su
            xss = dmatl(l1,l2,mlm,2,2,ksp) ! ss

            do  i = 1, nr
              rho(i,mlm,is1,is2) = rho(i,mlm,is1,is2)
     .           + xuu * ul(i,l1,i1) * ul(i,l2,i2)
     .           + xsu * sl(i,l1,i1) * ul(i,l2,i2)
     .           + xus * ul(i,l1,i1) * sl(i,l2,i2)
     .           + xss * sl(i,l1,i1) * sl(i,l2,i2)
            enddo

C       ... Local orbitals
            if (lpz1 .or. lpz2) then
              xuz = dmatl(l1,l2,mlm,1,3,ksp)
              xsz = dmatl(l1,l2,mlm,2,3,ksp)
              xzu = dmatl(l1,l2,mlm,3,1,ksp)
              xzs = dmatl(l1,l2,mlm,3,2,ksp)
              xzz = dmatl(l1,l2,mlm,3,3,ksp)
              if (xuz /= 0 .or. xsz /= 0 .or. xzu /= 0 .or.
     .            xzs /= 0 .or. xzz /= 0) then
              do  i = 1, nr
                rho(i,mlm,is1,is2) = rho(i,mlm,is1,is2)
     .            + xuz * ul(i,l1,i1) * gz(i,l2,i2)
     .            + xsz * sl(i,l1,i1) * gz(i,l2,i2)
     .            + xzu * gz(i,l1,i1) * ul(i,l2,i2)
     .            + xzs * gz(i,l1,i1) * sl(i,l2,i2)
     .            + xzz * gz(i,l1,i1) * gz(i,l2,i2)
              enddo
              endif
            endif

          enddo
        enddo
      enddo
      enddo
      enddo

C     call prmx('rho1',rho,nr*nlml,nr*nlml,nsp*nspc)

C --- Remake spherical density including small component ---
C     print *, 'skip remaking spher. density'
C     return
      lin = max(linrh1,linrh2)
      lx = max(lxrh1,lxrh2)
      do  ispc = 1, nspc        ! Spin off-diagonal parts with ispc=2
      do  isp = 1, nsp

        if (ispc == 1) then
          i1 = isp              ! ruu,rus,rss for spin diagonal block
          is1 = isp; is2 = 1    ! indices to spin diagonal part of rho(isp)
          ksp = isp
        elseif (ispc == 2) then
          i1 = 3                ! ruu,rus,rss for spin 12 block
          is1 = isp; is2 = 2    ! Indices to spin 12 block of rho:
          ksp = isp+2           ! Re rho(1,2) if isp=1, Im rho(1,2) if isp=2
        endif

        call dpzero(rho(1,1,is1,is2),nr)
        do  l = 0, lmxa

          if (lin >= 0 .and. l /= lin) then
            cycle
          elseif (lx >= 0 .and. l == lx) then
            cycle
          endif

          xuu = dmatl(l,l,1,1,1,ksp)
          xus = dmatl(l,l,1,1,2,ksp) + dmatl(l,l,1,2,1,ksp)
          xss = dmatl(l,l,1,2,2,ksp)
          xuz = dmatl(l,l,1,1,3,ksp) + dmatl(l,l,1,3,1,ksp)
          xsz = dmatl(l,l,1,2,3,ksp) + dmatl(l,l,1,3,2,ksp)
          xzz = dmatl(l,l,1,3,3,ksp)

          do  i = 1, nr
            rho(i,1,is1,is2) = rho(i,1,is1,is2) +
     .       xuu*ruu(i,l,1,i1) + xus*rus(i,l,1,i1) + xss*rss(i,l,1,i1)
          enddo

          if (pnz(l+1,1) /= 0) then
          do  i = 1, nr
            rho(i,1,is1,is2) = rho(i,1,is1,is2) +
     .       xuz*ruu(i,l,2,i1) + xsz*rus(i,l,2,i1) + xzz*rss(i,l,2,i1)
          enddo
          endif

      enddo
      enddo
      enddo

C     call prmx('rho1',rho,nr*nlml,nr*nlml,nsp*nspc)

      call tcx('mkrou2')
      end

      subroutine mkrou3(mode,lmxa,nlml,nsp,pnz,dmatl,hab,sab,qsum,hsum)
C- l-decomposed charges and eigenvalue sum
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 make qsum only
Ci         :2 make qsum and hsum both
Ci   lmxa  :augmentation l-cutoff
Ci   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
Ci         :expansion of the function products f_i(l1) f_j(l2)
Ci         :where f_i,f_j form this table.
Ci         :    (uu  us  uz)
Ci         :    (su  ss  sz)
Ci         :    (zu  zs  zz)
Ci   hab   :<u,s | h | u,s> for each pair uu, us, su, ss; see remarks
Ci   sab   :<u,s | 1 | u,s>
Co Outputs
Co   qsum  :l-decomposed sphere charge
Co   hsum  :l-decomposed one-electron energy
Cr Remarks
Cr   qsum and hsum are used to find the band centers of gravity
Cr   u and s ae linear combinations radial wave functions defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,lmxa,n0,nab,nlml,nsp
      parameter (n0=10,nab=9)
      double precision hab(nab,n0,nsp),sab(nab,n0,nsp),
     .  qsum(n0,nsp),hsum(n0,nsp),
     .  dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2)

C ... Local parameters
      integer l,isp,m
      double precision pi,srfpi,qz,hz

      if (mode == 0) return
C     call tcn('mkrou3')

      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)

      do  isp = 1, nsp
        do  l = 0, lmxa
          if (mode >= 1) qsum(l+1,isp) = 0d0
          if (mode >= 2) hsum(l+1,isp) = 0d0
        enddo
      enddo

      do  isp = 1, nsp
        do  l = 0, lmxa
          m = l+1
          qsum(m,isp) =
     .       + dmatl(l,l,1,1,1,isp)*sab(1,m,isp)*srfpi
     .       + dmatl(l,l,1,1,2,isp)*sab(2,m,isp)*srfpi
     .       + dmatl(l,l,1,2,1,isp)*sab(3,m,isp)*srfpi
     .       + dmatl(l,l,1,2,2,isp)*sab(4,m,isp)*srfpi
          if (mode >= 2) then
          hsum(m,isp) =
     .       + dmatl(l,l,1,1,1,isp)*hab(1,m,isp)*srfpi
     .       + dmatl(l,l,1,1,2,isp)*hab(2,m,isp)*srfpi
     .       + dmatl(l,l,1,2,1,isp)*hab(3,m,isp)*srfpi
     .       + dmatl(l,l,1,2,2,isp)*hab(4,m,isp)*srfpi
          endif
C         ... uz, sz, zu, zs, zz terms
          if (pnz(m,1) /= 0) then
          qz =
     .       + dmatl(l,l,1,1,3,isp)*sab(5,m,isp)*srfpi
     .       + dmatl(l,l,1,2,3,isp)*sab(6,m,isp)*srfpi
     .       + dmatl(l,l,1,3,1,isp)*sab(5,m,isp)*srfpi
     .       + dmatl(l,l,1,3,2,isp)*sab(6,m,isp)*srfpi
     .       + dmatl(l,l,1,3,3,isp)*sab(7,m,isp)*srfpi
          if (mode >= 2) then
          hz =
     .       + dmatl(l,l,1,1,3,isp)*hab(5,m,isp)*srfpi
     .       + dmatl(l,l,1,2,3,isp)*hab(6,m,isp)*srfpi
     .       + dmatl(l,l,1,3,1,isp)*hab(5,m,isp)*srfpi
     .       + dmatl(l,l,1,3,2,isp)*hab(6,m,isp)*srfpi
     .       + dmatl(l,l,1,3,3,isp)*hab(7,m,isp)*srfpi
          endif
          endif

        enddo
      enddo

C     call tcx('mkrou3')

      end

      subroutine mkrou4(nspch,nlml,cg,jcg,indxcg,
     .  nk1,norb1,ltab1,ktab1,blks1,lmx1,ndim1,
     .  nk2,norb2,ltab2,ktab2,blks2,lmx2,ndim2,
     .  qkk12,ckk)
C- Assemble contracted density-matrix (in Y_lm form)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspch :Number of spin channels
Ci         :1 in nonmagnetic case
Ci         :2 in spin-polarized case when spin channels are decoupled
Ci         :4 in spin-polarized case when spin channels are coupled
Ci   nlml  :charge density L-cutoff
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   nk1   :no. orbital types for a given L of first function
Ci         :nk1 merely dimensions qkk12,ckk
Ci   norb1 :total number of orbitals of the first type
Ci   ltab1 :table of l quantum numbers for the norb1 orbitals
Ci   ktab1 :table of k numbers (orbital type) for the norb1 orbitals
Ci   blks1 :block size for grouping orbitals into blocks (gtbls1)
Ci   lmx1  :dimensions ckk
Ci   ndim1 :dimensions qkk12; ndim1 should be (lmx1+1)**2
Ci   nk2   :no. orbital types for a given L of second function
Ci         :nk2 merely dimensions qkk12,ckk
Ci   norb2 :total number of orbitals of the first type
Ci   ltab2 :table of l quantum numbers for the norb1 orbitals
Ci   ktab2 :table of k numbers (orbital type) for the norb1 orbitals
Ci   blks2 :block size for grouping orbitals into blocks (gtbls1)
Ci   lmx2  :dimensions ckk
Ci   ndim2 :dimensions qkk12; ndim1 should be (lmx1+1)**2
Ci   qkk12 :density matrix between norb1 and norb2 orbitals
Co Outputs
Co   ckk   :product function coefficients
Cr Remarks
Cr   Wave function products have a radial part, written as products
Cr   of the radial functions (psi1)_l (psi2)_l'
Cr   The angular part Y_L Y_L' is contracted into a single index M
Cr   using Clebsch Gordan coefficients.  Thus ckk has indices
Cr      ckk(k1,k2,l1,l2,M),
Cr   where k1 ranges over the orbital types of the first function
Cr   and   k2 ranges over the orbital types of the second function.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nspch,nlml,norb1,lmx1,ndim1,norb2,lmx2,ndim2,nk1,nk2
      integer jcg(*),indxcg(*)
      integer ltab1(norb1),ktab1(norb1),blks1(norb1),
     .        ltab2(norb2),ktab2(norb2),blks2(norb2)
      double precision qkk12(nk1,nk2,ndim1,ndim2,nspch),cg(*),
     .                 ckk(nk1,nk2,0:lmx1,0:lmx2,nlml,nspch)
C ... Local parameters
      integer ilm1,io1,l1,nlm11,nlm12,k1,
     .        ilm2,io2,l2,nlm21,nlm22,k2,
     .        icg,ll,mlm,ix,isp
      double precision xx

C     call tcn('mkrou4')

      do  isp = 1, nspch
      do  io2 = 1, norb2
      if (blks2(io2) /= 0) then
C       k2,l2 = k and starting l index for this block
        l2 = ltab2(io2)
        k2 = ktab2(io2)
        nlm21 = l2**2+1
        nlm22 = nlm21 + blks2(io2)-1
        do  ilm2 = nlm21, nlm22
          l2 = ll(ilm2)
          do  io1 = 1, norb1
          if (blks1(io1) /= 0) then
C           k1,l1 = k and starting l index for this block
            l1 = ltab1(io1)
            k1 = ktab1(io1)
            nlm11 = l1**2+1
            nlm12 = nlm11 + blks1(io1)-1
            do  ilm1 = nlm11, nlm12
              l1 = ll(ilm1)

              ix = max0(ilm1,ilm2)
              ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
              do icg = indxcg(ix),indxcg(ix+1)-1
                mlm = jcg(icg)
                if (mlm <= nlml) then
                  xx = cg(icg)*qkk12(k1,k2,ilm1,ilm2,isp)
                  ckk(k1,k2,l1,l2,mlm,isp) = ckk(k1,k2,l1,l2,mlm,isp)+xx
                endif
              enddo
            enddo
          endif
          enddo
        enddo
      endif
      enddo
      enddo

C     call tcx('mkrou4')
      end

      subroutine mkcfus(nspch,lmxa,nlml,nf1,nf1s,val1,slo1,lmx1,
     .  nf2,nf2s,val2,slo2,lmx2,ckk,dmatl)
C- Assemble contracted density matrix as coffs to products of (u,s,gz)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspch :Number of spin channels
Ci         :1 in nonmagnetic case
Ci         :2 in spin-polarized case when spin channels are decoupled
Ci         :4 (or 6, for qkk -> qhk) when spin channels are coupled
Ci   lmxa  :leading dimension of dmatl
Ci   nlml  :charge density L-cutoff
Ci   nf1   :number of orbital types of first kind for each l
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there is a smooth part to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   val1  :function values at MT boundary, first function
Ci   slo1  :function slopes at MT boundary, first function
Ci   lmx1  :dimensions val1,slo1
Ci   nf2   :number of orbital types of second kind for each l
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :a smooth part is to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   val2  :function values at MT boundary, second function
Ci   slo2  :function slopes at MT boundary, second function
Ci   lmx2  :dimensions val1,slo1
Ci   ckk   :density matrix between  and  orbitals
Co Outputs
Co   dmatl :dmatl(l1,l2,mlm,i,j,:,:) holds coefficients to a Y_lm
Co         :expansion of the function products f_i(l1) f_j(l2)
Co         :where f_i,f_j form this table.
Co         :    (uu  us  uz)
Co         :    (su  ss  sz)
Co         :    (zu  zs  zz)
Co         :Spin indices:
Co         :  dmatl(:,:,:,:,:,1) : dmat of (1,1) spin block
Co         :  dmatl(:,:,:,:,:,2) : dmat of (2,2) spin block
Co         :  dmatl(:,:,:,:,:,3) : Re[dmat] of (1,2) spin block
Co         :  dmatl(:,:,:,:,:,4) : Im[dmat] of (1,2) spin block
Cr Remarks
Cr   Transforms contracted density-function matrix into coefficients
Cr   of wave function products. Radial parts are
Cr     (u, s, or gz)_l (u, s, or gz)_l'
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nspch,nlml,lmxa,lmx1,lmx2,nf1,nf2,nf1s,nf2s
      double precision ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nspch),
     .  val1(0:lmx1,nf1s),slo1(0:lmx1,nf1s),
     .  val2(0:lmx2,nf2s),slo2(0:lmx2,nf2s)
      double precision dmatl(0:lmxa,0:lmxa,nlml,3,3,nspch)
C ... Local parameters
      integer l1,k1,l2,k2,mlm,isp
      double precision xx

C     call tcn('mkcfus')

      do  isp = 1, nspch
      do  mlm = 1, nlml
        do  k2 = 1, nf2s
          do  k1 = 1, nf1s
            do  l2 = 0, lmx2
              do  l1 = 0, lmx1
                xx = ckk(k1,k2,l1,l2,mlm,isp)

                dmatl(l1,l2,mlm,1,1,isp) = dmatl(l1,l2,mlm,1,1,isp)
     .                             + xx * val1(l1,k1) * val2(l2,k2)
                dmatl(l1,l2,mlm,1,2,isp) = dmatl(l1,l2,mlm,1,2,isp)
     .                             + xx * val1(l1,k1) * slo2(l2,k2)
                dmatl(l1,l2,mlm,2,1,isp) = dmatl(l1,l2,mlm,2,1,isp)
     .                             + xx * slo1(l1,k1) * val2(l2,k2)
                dmatl(l1,l2,mlm,2,2,isp) = dmatl(l1,l2,mlm,2,2,isp)
     .                             + xx * slo1(l1,k1) * slo2(l2,k2)

              enddo
            enddo
          enddo
        enddo
      enddo
      enddo

C     call tcx('mkcfus')

C --- Products involving local orbitals ---
      if (nf1s >= nf1 .and. nf2s >= nf2) return
C     call tcn('mkcfus')

      do  isp = 1, nspch
      do  mlm = 1, nlml
        do  k2 = 1, nf2
          do  k1 = 1, nf1
            if (k1 > nf1s .or. k2 > nf2s) then
            do  l2 = 0, lmx2
              do  l1 = 0, lmx1
                xx = ckk(k1,k2,l1,l2,mlm,isp)

C               sc-sc product
                if (k1 > nf1s .and. k2 > nf2s) then
                  dmatl(l1,l2,mlm,3,3,isp) = dmatl(l1,l2,mlm,3,3,isp)+xx

C               sc-valence product
                elseif (k1 > nf1s) then
                  dmatl(l1,l2,mlm,3,1,isp) = dmatl(l1,l2,mlm,3,1,isp)
     .                                     + xx * val2(l2,k2)
                  dmatl(l1,l2,mlm,3,2,isp) = dmatl(l1,l2,mlm,3,2,isp)
     .                                     + xx * slo2(l2,k2)

C               valence-sc product
                elseif (k2 > nf2s) then
                  dmatl(l1,l2,mlm,1,3,isp) = dmatl(l1,l2,mlm,1,3,isp)
     .                                     + val1(l1,k1) * xx
                  dmatl(l1,l2,mlm,2,3,isp) = dmatl(l1,l2,mlm,2,3,isp)
     .                                     + slo1(l1,k1) * xx

                endif

C                if (l1 == 2.and.l2 == 2 .and. mlm == 1) then
C                  print 333, k1,k2,xx,
C     .              dmatl(l1,l2,mlm,1,3,isp),dmatl(l1,l2,mlm,3,1,isp)
C  333             format(2i4,3f12.6)
C                endif

              enddo
            enddo
            endif
          enddo
        enddo
      enddo
      enddo

C     call tcx('mkcfus')

C      isp = 1
C      do  mlm = 1, nlml
C        print *, ' '
C        do  l2 = 0, lmx2
C          do  l1 = 0, lmx1
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,1,2,isp),
C     .        dmatl(l2,l1,mlm,2,1,isp),
C     .        dmatl(l1,l2,mlm,1,2,isp)-dmatl(l2,l1,mlm,2,1,isp)
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,1,3,isp),
C     .        dmatl(l2,l1,mlm,3,1,isp),
C     .        dmatl(l1,l2,mlm,1,3,isp)-dmatl(l2,l1,mlm,3,1,isp)
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,2,3,isp),
C     .        dmatl(l2,l1,mlm,3,2,isp),
C     .        dmatl(l1,l2,mlm,2,3,isp)-dmatl(l2,l1,mlm,3,2,isp)
C
C  334       format(3i4,2f12.6,f14.8)
C          enddo
C        enddo
C      enddo
C      pause

      end

      subroutine mkrou5(nspch,nr,nlml,nf1,nf1s,f1,lmx1,nf2,nf2s,f2,lmx2,
     .  ckk,rho)
C- Assemble smooth site density from contracted density matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspch :Number of spin channels
Ci         :1 in nonmagnetic case
Ci         :2 in spin-polarized case when spin channels are decoupled
Ci         :4 in spin-polarized case when spin channels are coupled
Ci   nlml  :charge density L-cutoff
Ci   nf1   :number of orbital types of first kind for each l
Ci         :nf1 merely dimensions ckk
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there are functions of f1 type defined
Ci   f1    :first function on a radial mesh
Ci   lmx1  :l-cutoff for f1 functions
Ci   nf2   :number of orbital types of second kind for each l
Ci         :nf2 merely dimensions ckk
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :there are functions of f2 type defined
Ci   f2    :second function on a radial mesh
Ci   lmx2  :l-cutoff for f2 functions
Ci   ckk   :density matrix between f1 and f2 orbitals
Co Outputs
Co   rho   :density assembled on mesh
Cr Remarks
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
Cu             Envelopes f1,f2 must be zero for all channels that
Cu             have no smooth counterparts to subtract.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nspch,nr,nlml,lmx1,lmx2,nf1,nf2,nf1s,nf2s
      double precision ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nspch),
     .  f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s),rho(nr,nlml,nspch)
C ... Local parameters
      integer l1,k1,l2,k2,mlm,isp,i
      double precision xx

      call tcn('mkrou5')

      do  isp = 1, nspch
      do  mlm = 1, nlml
        do  k2 = 1, nf2s
          do  k1 = 1, nf1s
            do  l2 = 0, lmx2
              do  l1 = 0, lmx1
                xx = ckk(k1,k2,l1,l2,mlm,isp)
                if (dabs(xx) > 1d-16) then
                do  i= 1, nr
                  rho(i,mlm,isp) = rho(i,mlm,isp)
     .                           + xx*f1(i,l1,k1)*f2(i,l2,k2)
                enddo
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      enddo

      call tcx('mkrou5')
      end

      subroutine mkrou6(rofi,rho,nr,nlml,nsp,rho0,decay,res)
C- Fit tail of spin density; integrate charge beyond MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   rofi  :radial mesh points
Ci   rho   :spin-polarized charge density
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Co Outputs
Co   rho0  :fit density of form rho0*exp(-decay*r)
Co   decay :fit density of form rho0*exp(-decay*r)
Co   res   :integral of fit density from rofi(nr) to infinity
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Sep 03  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nlml,nsp
      double precision rofi(nr),rho(nr,nlml,2),rho0,decay,res
C ... Local parameters
      integer ir
      double precision norm(2,2),tnorm(2,2),rhs(2),y,dy,fac,y0,r0,pi,a,b

      res = 0
      if (nr < 10 .or. nsp == 1) return
      call dpzero(norm,4)
      call dpzero(rhs,2)
      fac = 1
      if (rho(nr,1,1) < rho(nr,1,2)) fac = -1
      do  ir = nr-5, nr
        y = fac*(rho(ir,1,1) - rho(ir,1,2))/rofi(ir)**2
C       If the spin density changes sign, nonsensical to try and fit
        if (y <= 0) return
        dy = dlog(y)

        norm(1,1) = norm(1,1) + 1
        norm(1,2) = norm(1,2) + rofi(ir)
        norm(2,1) = norm(2,1) + rofi(ir)
        norm(2,2) = norm(2,2) + rofi(ir)**2

        rhs(1) = rhs(1) + dy
        rhs(2) = rhs(2) + rofi(ir)*dy
      enddo

      call dinv22(norm,tnorm)
      a = tnorm(1,1)*rhs(1) + tnorm(1,2)*rhs(2)
      b = tnorm(2,1)*rhs(1) + tnorm(2,2)*rhs(2)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      a = fac*exp(a)/y0
      b = -b

C     Nonsensical if density not decaying fast enough
      if (b < 1) return
      r0 = rofi(nr)

C      a = 1d6
C      b = 2.5d0
C      r0 = 3

C     Integral a*exp(-b*r)*r*r = (2+2*b*r0+b**2*r0*2)/b**3*exp(-b*r0)
      res = a*(2+2*b*r0+b**2*r0**2)/b**3*exp(-b*r0)
      decay = b
      rho0 = a

      end
