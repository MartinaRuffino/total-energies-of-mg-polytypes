      subroutine vcdmel(s_ctrl,s_site,s_spec,s_lat,nlmax,ndham,nphimx,ndimh,
     .  nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn)
C- Valence-core dipole matrix elements
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu v0 pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name a nr rmt z lmxa lmxb pz orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg jcg indxcg
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
Co Outputs:
Co   weights for each channel output in iomoms style
Cr Remarks
Cu Updates
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   08 Jul 08 Dimension aus separately from ndimh
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   19 Sep 03 (ATP) Bug fixes
Cu   28 Mar 01 (MvS) rearrangement of indices to accommodate new makusq.
Cu   20 Mar 01 (ATP) extended to handle multiple core levels
Cu   20 Feb 01 Written by ATP
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nlmax,ndham,nphimx,ndimh,nq,nsp,nspc,nsite
      integer isite(nsite),iclsl(nsite),iclsn(nsite)
      double precision ef,evl(ndham,nsp,nq)
      double complex aus(nlmax,ndham,nphimx,nsp,nsite,nq)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:),wt(:)
      real(8), allocatable :: ul(:),sl(:),g(:),gz(:)
      real(8), allocatable :: ruu(:),rus(:),rss(:)
      real(8), allocatable :: s(:)
C ... Local parameters
      integer n0,lmxax
      parameter (n0=10,lmxax=10)
      integer ifi,fopna,isp,ib,is,lcls,ncls,nl,i,j,iomoms,intopt,nglob,
     .        nr,lmxa,iq,nlma,i1mach,nfstg,nchan
      integer lh(10)
C     integer orofi,owt,oul,osl,ogz,oruu,orus,orss,og,os
      double precision pnu(n0,2),pnz(n0,2),a,rmt,z,xx,rsml(n0),ehl(n0)
      double precision ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite)
      real(8), pointer    :: p_v0(:,:)
      character clabl*8

      call tcn ('vcdmel')
      call dpzero(rsml,n0)
      call dpzero(ehl,n0)

      do  i = 1, nsite
        ib = isite(i)
        ncls = iclsn(i)
        lcls = iclsl(i)
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        p_v0 => s_site(ib)%v0
        pnz = s_site(ib)%pz
        clabl = s_spec(is)%name
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        z = s_spec(is)%z

        lmxa = s_spec(is)%lmxa
        if (lmxa > lmxax) call rxi('vcdmel needs lmxax ',lmxa)
        if (lmxa == -1) cycle
        allocate(rofi(nr),wt(nr))
        intopt = 10*nglob('lrquad')
        call radmwt(intopt,rmt,a,nr,rofi,wt)
C   --- Augmented wave functions u,s
        j = nr*(lmxa+1)*nsp
        allocate(ul(j),sl(j),gz(j),ruu(j*2),rus(j*2),rss(j*2))
        call uspecb(0,4,s_spec,is,is,lh,rsml,ehl,j)
        call makusp(n0,z,nsp,1,rmt,lmxa,p_v0,a,nr,xx,xx,pnu,pnz,rsml,
     .    ehl,ul,sl,gz,ruu,rus,rss)
C   --- Matrix elements of u,s with core
        call skpblb(clabl,8,j)
        call awrit3('%N CLS atom %i ('//clabl(1:j+1)//') n=%i l=%i',
     .              ' ',120,i1mach(2),ib,ncls,lcls)
        allocate(g(nr*2))
        call pvcdm1(ncls,lcls,g,z,lmxa,p_v0,a,nr,rofi,wt,
     .              ul,sl,nsp,lmxax,ume(0,1,i),sme(0,1,i))
        deallocate(ul,sl,gz,ruu,rus,rss,g,rofi,wt)
      enddo

C --- Open CLS weights file and write first line
      allocate(s(3*ndimh*nsite*2))
      nl = s_ctrl%nl
      nfstg = 11
      nchan = 3*nsite
      ifi = fopna('cls',-1,4)
      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,0,1,0,0,0,0,0,
     .           0d0,0d0,0d0,0d0,0d0,0d0)

C --- For each qp, make <nk|x,y,z|core> at each site and save to disk in
C     iomoms-compatible format ---
      do   iq = 1, nq
        do  isp = 1, nsp
          call dpzero(s,3*ndimh*nsite*2)
          do  i = 1, nsite
            lcls = iclsl(i)
            ib = isite(i)
            is = s_site(ib)%spec
            lmxa = s_spec(is)%lmxa
            nlma = (lmxa+1)**2
            if (lmxa > -1) then
              call pvcdm2(i,nsite,ndham,ndimh,nlma,nlmax,
     .                  aus(1,1,1,isp,i,iq),ume(0,isp,i),sme(0,isp,i),
     .                  lcls,s_lat%cg,s_lat%jcg,s_lat%indxcg,s)
            endif
          enddo

C     ... Scale weights arbitrarily by 100 for plotting etc
          call dscal(3*ndimh*nsite*2,1d2,s,1)
          j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,1,1,ndimh,ndimh,
     .      nchan,nchan,ndimh,evl(1,isp,iq),0d0,s,0d0,0d0,0d0)
        enddo
      enddo
      iq = nq*nsp + 1
      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,iq,1,ndimh,ndimh,
     .  nchan,nchan,ndimh,0d0,0d0,0d0,0d0,ef,0d0)
      deallocate(s)
      call fclose(ifi)
      call tcx ('vcdmel')
      end
      subroutine pvcdm1(ncls,lcls,gcore,z,lmxa,v,a,nr,rofi,wt,
     .                  ul,sl,nsp,lmxax,ume,sme)
C- Radial matrix elements < (u,s) | r | core >
      implicit none
      integer ncls,lcls,lmxa,nr,nsp,lmxax
      double precision a,z,gcore(nr,2),rofi(nr),wt(nr),v(nr,nsp),
     .  ul(nr,0:lmxa,nsp),sl(nr,0:lmxa,nsp),ume(0:lmxax,nsp),
     .  sme(0:lmxax,nsp)
C Local
      integer nodes,l,nre,isp,ll,ir,i1mach,nglob,stdo,kc
      double precision e1,e2,slo,val,rmax,b,ecore,tol,yyy,dlml,slo1,
     .  r,wgt,uc,sc,ecor0,sum

      stdo = nglob('stdo')

      do  isp = 1, nsp
        if (nsp == 2) then
          call info2(30,0,0,' Spin %i ..',isp,0)
        endif
C   --- gcore <- core level wave function * r ---
        tol = 1.d-8
        e1 = -2.5d0*z*z - 5
        e2 = 20.d0
        val = 1.d-30
        slo = -val
        l = lcls
        rmax = rofi(nr)
        b = rmax/(dexp(a*nr-a)-1.d0)
        nodes = ncls - (l+1)
        ecore = (e1+e2)/2
        call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,v(1,isp),gcore,sum,
     .            a,b,rofi,nr,nre,kc)
        ecor0 = ecore
C   ... Correct core energy by using hankel bc's
        yyy = ecore - v(nr,isp) + 2*z/rmax
        if (nre == nr .and. yyy < 0.d0) then
          dlml = -1.d0-dsqrt(-yyy)*rmax
          do  ll = 1, l
            dlml = -yyy*rmax*rmax/dlml - (2*ll+1)
          enddo
          slo1 = val*(dlml+l+1)/rmax
          call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo1,v(1,isp),gcore,
     .              sum,a,b,rofi,nr,nre,kc)
        endif
        call awrit2(' vcdmel: ecor0=%,8d ecore=%,8d%N'
     .            //' (not including electrostatic potential shift)',
     .              ' ',120,i1mach(2),ecor0,ecore)
C      call prrmsh('gcore',rofi,gcore,nr,nr,1)
C      call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)
C      call prrmsh('sl',rofi,sl,nr,nr,1+lmxa)

C --- Matrix elements < (u,s) | r | core > ---
        write(stdo,1)
    1   format('   l',3x,'<u|core>',5x,'<s|core>',4x,'<u|r|core>',2x,
     .          '<s|r|core>')
        do  l = 0, lmxa
          ume(l,isp) = 0
          sme(l,isp) = 0
          uc = 0
          sc = 0
          do  ir = 1, nre
            r = rofi(ir)
            wgt = wt(ir)
            uc     = uc + wgt * ul(ir,l,isp) * gcore(ir,1)
            sc     = sc + wgt * sl(ir,l,isp) * gcore(ir,1)
            ume(l,isp) = ume(l,isp) + wgt*ul(ir,l,isp)*r*gcore(ir,1)
            sme(l,isp) = sme(l,isp) + wgt*sl(ir,l,isp)*r*gcore(ir,1)
          enddo
          write(stdo,2) l,uc,sc,ume(l,isp),sme(l,isp)
    2     format(i4,4F12.6)
        enddo
      enddo
      end
      subroutine pvcdm2(isite,nsite,ndham,ndimh,nlma,nlmax,aus,ume,sme,
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
      double precision cg(*),ume(0:*),sme(0:*),s(3,nsite,ndimh,2)
      double complex aus(nlmax,ndham,2)
C ... Local parameters
      integer kk(4),mlm,lm,ll,klm,ii,indx,icg1,icg2,icg,llm,ib,k
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
C                  write (*,200) kk(klm),mlm,klm,llm,aus(mlm,ib,1),
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

      do  k = 1, 3
        do  ib = 1, ndimh
        s(k,isite,ib,1) = s(k,isite,ib,1)*s(k,isite,ib,1)
     .                  + s(k,isite,ib,2)*s(k,isite,ib,2)
        enddo
      enddo
C  100 format (2x,'n',2x,'L(val)',1x,'L(r)',1x,'L(core)',11x,'au',25x,
C     .        'as',14x,'CG',20x,'cxx')
C  200 format (1x,i3,3i5,2x,'(',g12.5,',',g12.5,')',
C     .  '(',g12.5,',',g12.5,')',g12.5,'(',g12.5,',',g12.5,')')

      end
