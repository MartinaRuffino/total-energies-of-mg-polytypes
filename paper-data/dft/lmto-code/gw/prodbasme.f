      subroutine prodbasmea(opt,s_site,s_spec,nl,nsp,nbas,nat,ndrphi,
     .  nrphi,nrmx,nrpb,rprodb,gtoto,rprbme)
C- Matrix elements of partial waves and product basis for all augmented sites
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa a rmt nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr
Cio    Passed to:  *
Ci Inputs
Ci  opt    :1s digit not used now.  Should be 0
Ci         :10s digit
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nbas   :size of basis
Ci  nat    :number of sites in basis with augmentation spheres
Ci  ndrphi :(dimensioning) global maximum number of radial product functions for a particular l
Ci         :Formerly named nn in old GW code
Ci  nrphi   :number of (core states + valence partial waves) for a particular l
Ci         :nrphi formerly named nindx in old GW code
Ci  nrmx   :leading dimension of gtoto, rprodb
Ci         :nrmx formerly named nrofi in old GW code
Ci  rprodb :radial parts of orthonormalized product basis functions B for all l and sites
Ci  gtoto  :orthogonalized partial waves, or some combination of them and core states
Ci  nrpb   :indexing for rprodb.
Ci         :Functions nrpb(l,iat):nrpb(l+1,iat) are the family of radial B functions
Ci         :stored in rprodb for quantum number l and site iat.
Ci         :In former GW code, nrpb(l+1,iat)-nrpb(l,iat) was stored in nxx
Co Outputs
Co  rprbme :radial matrix elements < gtoto gtoto B> for all sites
Cl Local variables
Cl  nblr   :number of radial product basis functions for each l
Cl         :Formerly named nxx in old GW code
Cr Remarks
Cu Updates
Cu   05 Jul 18 Adapted from old GW basnfp_v2
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!     include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Passed parameters
      integer opt,nl,nsp,nbas,nat,ndrphi,nrmx
      integer :: nrphi(0:nl-1,nat)
      integer :: nrpb(0:2*(nl-1),nat+1) !dimensioned as a single vector to avoid compiler complaints
      double precision gtoto(nrmx,0:nl-1,ndrphi,nsp,nat)
      double precision rprodb(nrmx,*),rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),*)
C ... Local parameters
      double precision rofi(nrmx),rwgt(nrmx),xx(1)
      integer ib,iat,is,mxrbli,ntpbi,offpb,offpr,nblr(0:2*(nl-1))
      procedure(integer) :: nglob

      call info5(20,1,0,
     .  ' rprodbas:  radial matrix elements of partial waves and product basis, %i sites  nrphi=%i  lmax=%i',
     .  nat,ndrphi,nl-1,4,5)

      iat = 0; offpr = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        if (s_spec(is)%lmxa < 0) cycle
        iat = iat+1
        call rmeshprm(2+4,s_spec,nrmx,is,is,20,xx,xx,xx,xx,rofi,rwgt)
        call nrpb2nbl(03,nl,iat,nrpb,nblr,mxrbli,ntpbi) ! We need nblr and mxrbli
        offpb = 1+nrpb(0,iat)

C       call prmx('gtoto(nr,l=0)',gtoto(s_spec(is)%nr,0,:,1,1),ndrphi,ndrphi,1)

        call prodbasmei(0,nl,nsp,nblr,ndrphi,nrmx,s_spec(is)%nr,rofi,rwgt,nrphi,
     .    rprodb(1,offpb),gtoto(1,0,1,1,iat),gtoto(1,0,1,1,iat),rprbme(0,1,0,1,0,1+offpr))
        offpr = offpr + mxrbli
      enddo

C#ifdef DEBUG
      call info2(30,0,0,' sumcheck rprbme, all sites %;12F',sum(rprbme(:,:,:,:,:,1:offpr)),2)
C#endif

      end

      subroutine prodbasmei(mode,nl,nsp,nblr,ndrphi,nrmx,nr,ri,rwt,nrphi,rprodb,phi1,phi2,rprbme)
C- Matrix elements of partial waves and product basis functions
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :Not used now.  Should be 0
Ci  nl     :(global maximum l) + 1 --- here a dimensioning parameter
Ci  nsp    :2 for spin-polarized case, otherwise 1
Ci         :  In the spin pol case, product basis formed from spin average of partial waves
Ci  nblr   :number of radial product basis functions for each l
Ci  ndrphi :dimensioning parameter : max number of valence partial waves of a particular l
Ci         :ndrphi formerly named nn in old GW code
Ci  nrmx   :leading dimension of gtoto, rprodb
Ci         :nrmx formerly named nrofi in old GW code
Ci  nr     :Number of radial mesh points
Ci  ri     :ri(1:nr) = radial mesh
Ci  rwt    :rwt(1:nr) = radial mesh weights
Ci  nrphi   :number of (core states + valence partial waves) for a particular l
Ci         :nrphi formerly named nindx in old GW code
Ci  rprodb :radial parts of product basis functions B for this site
Ci  phi1   :(r * partial waves of 1st kind) entering into matrix elements
Ci  phi2   :(r * partial waves of 2nd kind) entering into matrix elements
Co Outputs
Co  rprbme :radial matrix elements < phi1 phi2 B> for one site
Cl Local variables
Cr Remarks
Cu Updates
Cu   05 Jul 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nl,nsp,nrmx,nr,ndrphi
      integer nblr(0:2*(nl-1)),nrphi(0:nl-1)
      double precision rprodb(nrmx,*),ri(nr),rwt(nr)
      double precision rprbme(0:nl-1,ndrphi,0:nl-1,ndrphi,0:2*(nl-1),maxval(nblr))
      double precision phi1(nrmx,0:nl-1,ndrphi,nsp),phi2(nrmx,0:nl-1,ndrphi,nsp)
C ... Dynamically allocated local arrays
      integer,allocatable:: iprlc(:),lprc(:) !,ibo(:)
C ... Local parameters
      logical lmagbas
      integer irad,nradpb,lb,nb,ntpba,l1,isp,isp1,isp2,n1,l2,n2
      double precision rphiphi(nrmx)
      procedure(real(8)) :: dot3

      lmagbas = mod(mode,10) == 1

      nradpb = sum(nblr(0:2*(nl-1)))
      if (nradpb == 0) return
      allocate(iprlc(nradpb),lprc(nradpb))
      iprlc(1:nradpb) = 0

C ... For each radial function, make offset iprlc to full (m-resolved) P.B. and l-index lprc
      irad = 0
      do  lb = 0, 2*(nl-1)
        do  nb = 1, nblr(lb)
          irad = irad+1
          lprc(irad) = lb
          if (irad==1) cycle
          iprlc(irad) = iprlc(irad-1) + 2*lprc(irad-1)+1
        enddo
      enddo
      if (irad /= nradpb) call rx('prodbasmei: counting mismatch')
      ntpba = iprlc(nradpb) + 2*lprc(nradpb)+1

C --- Calculate radial matrix elements ---
      do  isp = 1,nsp
        isp1 = isp
        isp2 = isp
        if (lmagbas .and. isp == 1) then
          isp1 = 1
          isp2 = 2
        elseif (lmagbas .and. isp == 2) then
          isp1 = 2
          isp2 = 1
        endif

        irad = 0
        do  lb = 0, 2*(nl-1)
        do  nb = 1, nblr(lb)
          irad = irad+1
          do  l1 = 0, nl-1
          do  n1 = 1, nrphi(l1)
          do  l2 = 0, nl-1
          do  n2 = 1, nrphi(l2)
            rprbme(l1,n1,l2,n2,lb,nb) = 0
            if (lb<abs(l1-l2) .or. l1+l2<lb) cycle
            rphiphi(1)    = 0d0
            rphiphi(2:nr) = phi1(2:nr,l1,n1,isp2) * phi2(2:nr,l2,n2,isp1)/ri(2:) ! phi = u = r \phi
            rprbme(l1,n1,l2,n2,lb,nb) = dot3(nr,rprodb(1,irad),rphiphi,rwt)
          enddo
          enddo
          enddo
          enddo
        enddo                   ! nb
      enddo                     ! lb

        if (irad /= nradpb) call rx('prodbasmei : counter mismatch')
C       call prmx('rprbme',rprbme,nl*ndrphi*nl*ndrphi,nl*ndrphi*nl*ndrphi,(2*nl-1)*maxval(nblr))

      enddo  ! isp

C... prodmt proddmt (value and slope of the product at MT). !oct2005
!    Stored into the tail of PPBRD* files.

C --- MixSpin= <rho_up - rho_down | B> matrix calculation. May2005
Cr1  Suppose  "ibas==iclass"--- it is already checked in hbasfp0.m.f
Cr2  ValMT.* is written with subroutine savemtval(ib,rho1,rofi,nr,nlml,nsp)
Cr    in fp/locpot.f just befor locpt2 in lmto (lmf).
C      if(lmagbas) then
C        sqrtfpi = sqrt(fpi)
C        ifv = 6301
C        ibas=ic
C
C        if(valmt) then
C          open(ifv,file='ValMT.'//charnum3(ibas),form='unformatted')
C          read(ifv) nr_r,nlml_r,nsp_r
C          write(6,"('readin nr nlml nsp=',3i5)") nr_r,nlml_r,nsp_r
C          allocate(rofi_r(nr_r),rho1(nr_r,nlml_r,nsp_r),rspin(nrmx),den(nrmx,nsp),r11(nrmx))
C          r11(1:nrmx)= 1d0
C          read(ifv) rofi_r, rho1
C          close(ifv)
C        else
Cc
C          open(ifv,file='rhoMT.'//trim(charnum3n(ibas)),form='unformatted',status='old',err=1031)
C          goto 1032
C 1031     continue !bug fix for lmfgw---remove this path in future.
C
CC         if(ibas>9) call rx( 'rhoMT indexing is over 10')
C          if(ibas>9) then
C            open(ifv,file='rhoMT.'//char(48+ibas/10)//char(48+mod(ibas,10)),form='unformatted')
C          else
C            open(ifv,file='rhoMT.'//char(48+ibas),form='unformatted')
C          endif
C
C
C 1032     continue
C          read (ifv) nr_r
C          allocate(rofi_r(nr_r))
C          read (ifv) rofi_r
C          read (ifv) nr_r,nlmlsp,kxx,kxx,nsp_r
C          write(6,*)' rho1 xxx=', nr_r,nlmlsp,kxx,kxx,nsp_r
C          nlml_r = nlmlsp/nsp_r
C          allocate( rho1(nr_r,nlml_r,nsp_r),rspin(nrmx),den(nrmx,nsp),r11(nrmx))
C          r11(1:nrmx)= 1d0
C          read (ifv) rho1
C        endif
C
Cccccccccccccccccccccc
Cc        print *, 'sum rho1=', sum (abs(rho1(1:nr_r,1:nlml_r,1:nsp_r)))
Ccccccccccccccccccccc
C
C        if(nsp_r/=nsp) call rx( " ReadinError: ValMT: nspr/= nsp")
C        if(nsp/=2    ) call rx( " This mode is only for nsp==2")
C        rho1= sqrtfpi*rho1  !rho1 is not including sqrt(fpi) Right?
C
CC
C        open(ifv,file='MixSpin.'//charnum3(ibas))
C        write(ifv,"(2i10,' ! ibas, max l of product basis' )") ibas,2*(nl-1)
C        write(ifv,"(i10,'           ! nblr(lb)'  )") nblr(0:2*(nl-1))
C        do ilmx = 1, (2*(nl-1)+1)**2
C          lb = ll(ilmx)
C          if(ilmx <=nlml_r) then
C            rspin(1) = 0d0  !rspin = rho^{true spin density} * r
C            do ir =2,nrmx
C              den(ir,1)=  polinta(r(ir), rofi_r,rho1(:,ilmx,1),nr_r)
C              den(ir,2)=  polinta(r(ir), rofi_r,rho1(:,ilmx,2),nr_r)
C              rspin(ir)  = (den(ir,1)-den(ir,2))/r(ir)
Ccccccccccccccccccccc
Cc          write(6,"(' den=',3d13.6)") r(ir),den(ir,1:2)
Ccccccccccccccccccccc
C            enddo
C            den(1,1:2)=0d0
C          else
C            rspin=0d0
C            den=0d0
C          endif
CC ... sumcheck
Cccccccccccccccccccccc
Cc        print *,' nsp=',nsp,nrmx,nr_r
Cc        print *, 'sumchk den1=', sum ( abs(den(1:nr_r,1)) ),maxval(abs(den(1:nr_r,1)))
Cc        print *, 'sumchk den2=', sum ( abs(den(1:nr_r,2)) ),maxval(abs(den(1:nr_r,2)))
Cc        stop 'xxxxxxxxxxxxxx'
Ccccccccccccccccccccc
C
Cc den = 4 pi r^2 * rho_true(r)
Cc rspin = 4 pi r * rho_true(r)
C          if( nblr(lb)/=0) then
C            do isp=1,nsp
C              call gintxx( den(1,isp), r11, aa,bb,nrmx, sumc(isp) )
C            enddo
C            write(6,"(' charge: ilm charge=',i5,2f13.6)") ilmx,sumc(1:nsp)
C          endif
Ccccccccccccccccccccccccc
Cc          if(lb==0) then
Cc            do ir=1,nrmx
Cc              write(6,"('rrr ',2f13.6)") r(ir), r(ir)*rspin(ir)
Cc            enddo
Cc          endif
Cccccccccccccccccccccc
C          bb1s=0d0
C          do nb = 1, nblr(lb)
C
CcCase1 -----------
C            call gintxx( rprodb(1,nb,lb), rspin,aa,bb,nrmx,spinvec0)
C            spinvec = spinvec0/sqrtfpi
C
C!2007
C! const = <1|B> where 1 is normalized within the sphere
C            if(lb==0) then
C              call gintxx( rprodb(1,nb,lb), r,aa,bb,nrmx,const)
CC             call gintxx( r, r,aa,bb,nrmx,const )
C            else
C              const=0d0
C            endif
C            const= const * sqrtfpi !/((fpi/3d0)*r(nrmx)**3)
C
C! Now spinvec = <B_I(\bfr) | m_true(\bfr) >
C            if(abs(spinvec)<1d-10 ) spinvec=0d0
C            write(ifv,"(     2i6,d24.16,2x,f13.10,2x,f13.10,d24.16
C     &       ' ! I=(ilm, nb), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')")
C     &       ilmx, nb, spinvec, sumc(1:nsp),const
C            write(6,"('ttt:',2i6, d24.16,2x, 2f14.10,d24.16
C     &       ' ! I=(ilm, nb), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')")
C     &       ilmx, nb, spinvec, sumc(1:nsp),const
C
CcCase2 tested--- too bad result -------
Cc            call gintxx( rprodb(1,nb,lb), r,aa,bb,nrmx,
Cc     &        bb1)
Cc            bb1= bb1*sqrt(4*pi) ! bb1= bb1/sqrt(1d0/3d0*r(nrmx)**3)
Cc            if(lb/=0) bb1=0d0
Cc            bb1s=  bb1s+bb1*bb1 !spinvec
Cc            write(6,"('ttt:',2i6, d24.16,2x, 2f14.10,
Cc     &       ' ! I=(ilm, nb), <1|B(I)>, sum of <1|B(I)> ')")
Cc     &      ilmx, nb, bb1, bb1s
Cc           write(ifv,"(2i6,d24.16,2x,f13.10,2x,f13.10,
Cc     &       ' ! I=(ilm, nb), <1|B(I)> intrho(1:nsp)')")
Cc     &       ilmx, nb, bb1, sumc(1:nsp)
Cc        print *,' test: basn: <1|B> ------ '
C
C          enddo
C        enddo
C        deallocate(rofi_r,rho1,rspin)
C        do ilmx = 1, (2*(nl-1)+1)**2
C          do nb = 1, nblr(lb)
C          enddo
C        enddo
C        close(ifv)
C      endif !ixc==8 end

      deallocate(iprlc,lprc)

C#ifdefC DEBUG
C      call info2(30,0,0,' sumcheck rprbme %;12F',sum(rprbme(:,:,:,:,:,1:maxval(nblr))),2)
C#endif

      end
