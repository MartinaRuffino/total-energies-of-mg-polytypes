      subroutine chkstr(s_ctrl,s_lat,s_spec,s_str)
C- Plots screened envelope functions for checking
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbasp nspec nl ips
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat avw pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos cy cg indxcg jcg
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxb kmxt rsma hcr rmt hsfitk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:name
Cio    Passed to:  e0parm
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  nkaps nitab npr iax kaps
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:nkaps kaps npr s iax alph nitab
Cio    Passed to:  *
Ci Inputs
Co Outputs
Cs Command-line switches
Cs   --atom=#  : Check functions centered at site #
Cs   --edge=#  : ?
Cs   --ikap=#  : Which Hankel energy to check
Cs   --ixi=#   : which neighbor in cluster relative to atom
Cs             : ixi=1 for head, 2 for 1st NN, etc
Cs   --plot:.. : See Remarks
Cs   --radius=#: Evaluate one-center expansion at radius=#
Cs             : instead of hard core radius.
Cl Local variables
Cl   mode  : plotting mode:
Cl         :  1  tabulate screened envelope in a plane
Cl         :  2  tabulate screened envelope on a line
Cl         :  3  check one-center expansion
Cl   lhsm  :  1 => transform to sm Hankels
Cr Remarks
Cr   Routine does one of:
Cr   1.  tabulate screened envelope in a plane for contour plots
Cr       Invoke with --plot:con
Cr   2.  tabulate screened envelope on a line
Cr       Invoke with --plot[:smh]:line[:v1x,v1y,v1z,v2x,v2y,v2z]
Cr   3.  check one-center expansion
Cr       Invoke with --plot[:smh]:onec
Cu Updates
Cu   10 Apr 15 Some bug fixes, new smh2 mode
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   09 Oct 11 chkstr reworked for mapping H -> smH
Cu   14 Sep 11 Reworked 1c expansion to include gradients
Cu   06 Sep 11 Started migration to f90 structures
Cu   23 Jul 08 (S Lozovoi)) Adapted to l-dependent augmentation radii
Cu             and species-dependent lmax
Cu   06 Aug 06 Revised to work with 2-kappa basis
Cu   19 May 04 Added --plot options
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_str)::   s_str
C ... Dynamically allocated arrays
      integer,pointer :: ips(:)
      integer,allocatable :: ntabg(:),lmxb(:),kmxt(:)
      real(8),pointer :: pos(:,:)
      real(8),allocatable :: hcrs(:,:),ehl(:,:),rsma(:),rtab(:,:),halp(:),s(:)
C ... Local parameters
      logical ltmp
      procedure(logical) :: isanrg,iostr1,cmdopt,a2bin
      integer, parameter :: npmax=51,nlmy=16,nvmax=1001,n0=10,niax=10
      integer a2vec,fopn,i,ii,iat,iatp,ifi,ikap,imax,ip,ixi,j,j1,j2,ll,il,ila,
     .  ilb,mode,nbasp,nl,nlf, np,nttab,nvec,lio,lfit2,nbasf,stdo,lio23,
     .  lio45,getdig,ncplx, loka,isw,is,nclus,nlma,nlmb,la,lb,nspec,nds,irot,
     .  lhsm,nttabf,offR,offRp
      double precision alat,avw,ctr(3),ddot,radius,tolke,tstmax,xx
      double precision fi(10),gi(10),gi0(10),gslj(3,25,2),gslk(3,25,2),hcr(10),
     .  plat(3,3),rsq(122),slj(25,2),slk(25,2),v1(6),v2(3),wp(122),xmrp(3),
     .  xp(3,122),xpt(122,3)
      double precision  yl(122,nlmy),
     .  v1chs(nlmy,nlmy,2),v1clhs(nlmy,nlmy,2),
     .  dv1chs(nlmy,nlmy,2),dv1clhs(nlmy,nlmy,2),
     .  dvalj(nlmy,nlmy),dvalja(nlmy,nlmy),
     .  dvaln(nlmy,nlmy),dvalna(nlmy,nlmy),
     .  vvalj(nlmy,nlmy),vvaln(nlmy,nlmy),lvalj(nlmy,nlmy)
      double precision reslin(nvmax,12),resk(-npmax:npmax,-npmax:npmax)
      equivalence (resk,reslin)
C     double precision v1chsp(nlmy,nlmy,2)
      equivalence (v1(4),v2)
      character outs*120,out2*120,dc*1,dc2*1
C     for error printing
      integer ier1,ier2,ixx
      double precision err
      logical :: lerr
      procedure(integer) :: nglob
c...deb
c debugging 1c expansion
C      integer nlmax
C      integer orsmh,orsma,okmx

C     double precision, allocatable  :: hvx(:,:,:),hlx(:,:,:)
C     double precision, allocatable  :: hsvx(:,:,:,:),hslx(:,:,:,:)
c...deb

      stdo = nglob('stdo')
      nbasp = s_ctrl%nbasp
      nspec = s_ctrl%nspec
      nl = s_ctrl%nl
      alat = s_lat%alat
      plat = s_lat%plat
      avw = s_lat%avw
      ips => s_ctrl%ips

      allocate(lmxb(nspec),kmxt(nspec))
      allocate(rsma(nspec),hcrs(nl,nspec))
      do  is = 1, nspec
        lmxb(is) = s_spec(is)%lmxb
        kmxt(is) = s_spec(is)%kmxt
        rsma(is) = s_spec(is)%rsma
        hcrs(1:nl,is) = s_spec(is)%hcr(1:nl)
      enddo
      pos => s_lat%pos

C ... Read some strux data from from disk
      call info0(10,1,0,' CHKSTR: check or plot strx from file STR ...')
      lio = 0
      ltmp = iostr1(lio,'STR',ii,nttabf,ii,nlf,nbasf,ii,ii,xx,ifi)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      ncplx = 1+getdig(lio23,0,2)
      loka = isw(lio45 < 2)

      call info8(10,0,0,
     .  ' file contains: '//
     .  '%?#(n==0)#2nd generation ##%-1j'//
     .  '%?#(n==1)#NMTO ##%-1j'//
     .  '%?#(n==2)#Screened ##'//
     .  '%?#n#sdot#strux#'//
     .  '%?#(n<2)#, OKA defs##'//
     .  '%?#(n==0)# (1-center coffs)##'//
     .  '%?#(n>1)#, %-1j%i energies##',
     .  lio45,
     .  getdig(lio23,3,2),
     .  lio45,
     .  getdig(lio23,2,2),
     .  s_str%nkaps,
     .  s_str%nkaps,
     .  s_str%kaps,
     .  0)
      lerr = isanrg(nbasf,nbasp,nbasp,'file:CHKSTR:','nbas',.true.)
      call info5(10,0,0,' %i pairs in neighbor table   nlf = %i   E =%n:1;4,4d%N',
     .  s_str%npr(nbasp+1),nlf,s_str%nkaps,s_str%kaps,0)
      nds = nlf**2

C ... Setup for rotating to value, slope.  irot=0 for no rotation
      irot = 0
      irot = 1
      call query('convert to value-slope functions (0=no)',2,irot)
      if (irot /= 0) irot = 100
      if (irot == 100) then
        i = nds*nbasp*s_str%nkaps**2
        allocate(halp(i)); call dpzero(halp,i)
        call makalp(nl,nds,nbasp,s_str%nkaps,s_str%kaps,hcrs,100*loka+2,ips,lmxb,xx,halp)
      else
        allocate(halp(1))
      endif

      if (loka == 1) then
        call dscal(nds**2*s_str%nitab*s_str%nkaps**2,-1d0,s_str%s,1)
      endif

C ... Convert 1-center expansion to strux matrix
      nttab = s_str%npr(nbasp+1)
      i = nds**2*s_str%nitab*s_str%nkaps**2
      allocate(s(i))
      call hanb2s(irot+0,nds,nbasp,s_str%nkaps,s_str%npr,s_str%iax,s_str%alph,halp,s,s_str%s)
c...deb
c     call shostr(nds,s_str%nitab,nbasp,plat,pos,20001,s_str%alph,s_str%iax,
c    .           s_str%npr,s,s_str%nkaps,s_str%nkaps,s_str%kaps,1,1d0)
c...deb
C     1-center also for (val,slo)
      if (irot == 100) then
        call hanb2s(irot+2,nds,nbasp,s_str%nkaps,s_str%npr,s_str%iax,s_str%alph,halp,s_str%s,xx)
      endif

C ... For now, set ntabg = full cluster size.
C     ntabg clusters need be only as large as the range of the Gaussians
      allocate(ntabg(nbasp))
      do  i = 1, nbasp
        ntabg(i) = s_str%npr(i+1) - s_str%npr(i)
      enddo

C ... Initialize some variables and defaults
      ltmp = .false.
      mode = 0
      lhsm = 0
      lfit2 = 2
      tstmax = alat
      if (cmdopt('--edge=',7,0,outs)) then
        j = 7
        call rxx(.not.a2bin(outs,tstmax,4,0,' ',j,len(outs)),'failed to parse '//outs)
        call info2(0,0,0,' square dimensioned %,1d au (--edge)',tstmax,0)
      endif
      v1(1) = -tstmax
      v1(2) = 0
      v1(3) = 0
      v2(1) = tstmax
      v2(2) = 0
      v2(3) = 0
      nvec = 501

C ... Get plotting(/testing) mode
      if (cmdopt('--plot',6,0,outs)) then
        out2 = outs(7:)
        dc = out2(1:1)
        if (dc /= ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 0
   50     continue
          j2 = j2+1
          if (out2(j2:j2) == dc) goto 50
          j1 = min(len(out2),j2)
          call nwordg(out2,0,dc//' ',1,j1,j2)
          if (j2 >= j1) then
            if (.false.) then
C            elseif (out2(j1:j1+2) == 'hsm')  then
C              lhsm = 1
            elseif (out2(j1:j1+3) == 'smh2')  then
              lhsm = 1; lfit2 = 2
            elseif (out2(j1:j1+3) == 'smh3')  then
              lhsm = 1; lfit2 = 3
            elseif (out2(j1:j1+3) == 'smh4')  then
              lhsm = 1; lfit2 = 4
            elseif (out2(j1:j1+3) == 'smh5')  then
              lhsm = 1; lfit2 = 5
            elseif (out2(j1:j1+2) == 'smh')  then
              lhsm = 1
            elseif (out2(j1:j1+2) == 'con')  then
              mode = 1
            elseif (out2(j1:j1+3) == 'line')  then
              mode = 0
              dc2 = out2(j1+4:j1+4)
              if (dc2 /= ' ') then
                ip = j2+1
                i = a2vec(out2,len(out2),ip,4,', '//dc,3,2,6,fi,v1)
                if (i /= 6) goto 52
                j2 = ip
              endif
            elseif (out2(j1:j1+3) == 'onec')  then
              mode = 2
            else
              goto 52
            endif
            goto 50
   52       continue
            call rxs2('chkstr: failed to parse --plot option: "',out2(1:j2),'%a ..."')
          endif
        endif
      endif

      if (lhsm == 1) then
        call info0(10,0,0,' Envelopes for sm Hankels ...  generating parameters')
        tolke = .08d0
        allocate(s_str%c01g(0:1,0:nl-1,nspec),s_str%rsmh(0:nl-1,nspec))
        j = 30 + lfit2
        call pshpr(51)
        call e0parm(s_spec,j,nl,s_str%kaps,tolke,s_str%rsmh,s_str%c01g)
        call poppr
        allocate(ehl(nl,nspec)) ! Because strg1c requires it
        call dvset(ehl,1,nl*nspec,s_str%kaps)
C       print *, '!! c01 = 0, rsmh = > small'; s_str%c01g = 0; s_str%rsmh(:,:) = .2d0
      else
        allocate(s_str%c01g(1,1,1),s_str%rsmh(1,1))
      endif

C ... Make rtab, plat in atomic units
      nttab = s_str%npr(nbasp+1)
      allocate(rtab(3,nttab))
      call dpzero(ctr,3)
      call mkrtab(100,alat,plat,pos,s_str%iax,nttab,ctr,rtab)
C     call prmx('rtab',rtab,3,3,nttab)

      ikap = 1
      if (cmdopt('--ikap=',7,0,outs)) then
        j = 7
        call rxx(.not.a2bin(outs,ikap,2,0,' ',j,len(outs)),'failed to parse '//outs)
      endif
      if (s_str%nkaps > 1) call query('ikap=',2,ikap)
      if (ikap > s_str%nkaps) call rx('ikap cannot be larger than nkap=')

      ixi = 2 ! Default augmentation sphere : first NN to head
      if (cmdopt('--ixi=',6,0,outs)) then
        j = 6
        call rxx(.not.a2bin(outs,ixi,2,0,' ',j,len(outs)),'failed to parse '//outs)
      endif

C ... Scaling by avw (applicable when OKA conventions are used)
      if (loka == 0) avw = 1
      call dscal(3*nttab,1/avw, rtab,1)
      call dscal(s_str%nkaps,avw**2,s_str%kaps,1)

      if (mode == 0) goto 70
      if (mode == 1) goto 60
      if (mode == 2) goto 20
      call rxi('chkstr: invalid mode',mode)

C --- Contour plot ---
   60 continue
      print 923
  923 format(/' 1. Contour plot of screened Hankel function')
      if (irot == 100 .and. s_str%nkaps == 1)  print 924
  924 format( '    function is normalized to 1 on head')
      if (irot == 100 .and. s_str%nkaps == 2) then
        call info2(10,0,0,'%4pGenerating %?#(n==1)#value#slope# function',ikap,0)
      endif

      call query('half length of square (a.u.)',4,tstmax)
      imax = npmax
      call query('number of mesh points=',2,imax)
      iat = 1

      call query('atom=',2,iat)
      if (iat > nbasp) call rxi('chkstr: atom cannot be larger than nbas=',nbasp)
      is = ips(iat)
      ilb = 1
      call query('orbital L to plot=',2,ilb)
      ixi = 1
      offR  = s_str%npr(iat)
      nlmb = s_str%iax(niax*offR+9)
      if (ilb > nlmb) call rxi('chkstr: L cannot be larger than nlmb=',nlmb)
      lb = ll(nlmb)
      hcr(1:nl) = s_spec(is)%hcr(1:nl)
      call info5(10,1,0,' Tabulating L=%i  '//
     .  '%?#(n==1)#smooth#2nd gen#  envelope function for atom %i:  hcr=%n:1;4,4d',
     .  ilb,lhsm,iat,lb+1,hcr)
      call info2(0,0,1,' File "plot" holds raster of points between +/- %da.u., or %;4d*alat',
     .  tstmax,tstmax/alat)

C --- Contour plot ---
      do  i = -imax, imax
        do  j = -imax, imax
          xmrp(1) =  (i*tstmax)/imax
          xmrp(2) =  (j*tstmax)/imax
          xmrp(3) =  1d-4
          call dscal(3,1/avw,xmrp,1)  ! avw=1 unless OKA conventions are used
          call strck(nlf**2,iat,ixi,rtab,s_str%kaps,ikap,s_str%nkaps,nbasp,
     .      s_str%alph,halp,s_str%iax,s_str%npr,s,s_str%s,s_str%nitab,ips,nl,
     .      s_str%rsmh,s_str%c01g,xmrp,nlf-1,lhsm*1000+irot+loka,slk,slj,gslk,gslj)
C          resk(i,j) = ddot(s_str%nkaps,slk(1,1),25,halp,nlf**2*nbasp)
           resk(i,j) = max(min(slk(ilb,1),999d0),-99d0)
         enddo
      enddo
      ifi = fopn('PLOT')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%5,6g)',resk,0,2*imax+1,2*imax+1,2*imax+1)

      deallocate(halp)

      return

C --- One-center expansion ---
   20 continue
C     Test expansion for R' corresponding to iax(ixi)
      print 926
  926 format(/
     .  ' 2. Compare numerically integrated YL expansions of H^a,'/
     .  ' its gradient and its',
     .  ' 1-center expansion on a neighboring sphere surface.'/' These',
     .  ' should be equal apart from numerical integration errors.')
      if (irot == 0) print 927
      if (irot /= 0 .and. s_str%nkaps == 1) print 929
      if (irot /= 0 .and. s_str%nkaps == 2) print 930
  927 format(
     .  ' When evaluated at radius=hcr,'/
     .  ' function should be H^0 on (head sphere)'/
     .  ' and vanish on tail spheres.'/)
  929 format(
     .  ' When evaluated at radius=hcr,'
     .  ' function value should be delta_RL,R''L''')
  930 format(
     .  ' When evaluated at radius=hcr,'/
     .  ' function value should be delta_RL,R''L'', slope 0, or:'/
     .  ' function slope should be delta_RL,R''L'', value 0 '
     .  ' (ikap=2)')

      if (loka == 1) then
        call info0(10,1,0,' CHKSTR (warning) unidentified numerical errors for loka=1')
      endif
      if (lio45 == 1) call rx('onec not implemented for NMTO')

      iat = 2  ! Choose as default 2nd atom as head sphere
      if (cmdopt('--atom=',7,0,outs)) then
        j = 7
        call rxx(.not.a2bin(outs,iat,2,0,' ',j,len(outs)),'failed to parse '//outs)
      endif
   23 call query('Show 1-c for cluster centered at atom',2,iat)
      if (iat > nbasp) call rxi('atom cannot be larger than nbas=',nbasp)

      call query('ixi= (1 for head, 2 for 1st NN, etc)',2,ixi)
      offR  = s_str%npr(iat)
      nclus = s_str%npr(iat+1) - offR
      if (ixi > nclus) call rxi('ixi cannot be larger than nclus=',nclus)
      nlma = s_str%iax(niax*(offR+ixi-1)+9)
      la = ll(nlma)
      offRp = ixi-1 + offR
      iatp = s_str%iax(niax*offRp+2)
C ... la, nlma - head of the cluster, lb, nlmb - expansion sphere
C     NB: 'a' traditionally refers to augmentation, 'b' to basis
C     This routine should swap meaning of symbols!
      nlmb = s_str%iax(niax*offRp+9)
      lb = ll(nlmb)
      hcr(1:nl) = s_spec(ips(iatp))%hcr(1:nl)
      is = ips(iat)
      call info5(10,1,0,' Expand head (ib=%i,ikap=%i,spec='//trim(s_spec(is)%name)//
     .  ' about ixi=%i (jb=%i,spec='//trim(s_spec(ips(iatp))%name)//')',
     .  iat,ikap,ixi,iatp,0)
      if (ixi == 1) then
        call info0(10,0,0,' ixi=1 => expand head')
      else
        call info5(10,0,0,' ixi=%i => expand tail  hcr=%1;3,3d ... R(dest-src) =%3:1;8,5D',
     .             ixi,hcr,rtab(:,offR+ixi)-rtab(:,offR+1),0,0)
      endif

      call info2(10,1,1,' In the table below:'//
     .  '%N Ha    1c expansion from num. int(YL*Ha) on sphere surface'//
     .  '%?#(n==1 )#'//
     .  '%N 1cs   analytic 1c expansion where H -> smH'//
     .  '%N 1cu   analytic 1c expansion for unsm H#'//
     .  '%N 1c    analytic 1-center expansion#%-1j'//
     .  "%N H'a   dHa/dr, from int[YL*grad(Ha).r] on sphere surface"//
     .  "%N H'n   dHa/dr, from numerical radial derivative of Ha"//
     .  '%?#(n==1)#'//
     .  "%N 1c's  1c dHa/dr analytic radial derivative of 1c"//
     .  "%N 1c'u  1c dHa/dr analytic radial der. of unsm. 1c"//
     .  "%N 1c'n  1c dHa/dr from num. radial der. of unsm. 1c#"//
     .  "%N 1c'a  1c dHa/dr from int[YL*grad(Ja).r]"//
     .  "%N 1c'n  1c dHa/dr from numerical radial derivative of 1c#",
     .  lhsm,lhsm)

      if (lhsm == 1) then
        call info0(10,0,1,
     .    " 1c''n Laplacian Ha (analytic)"//
     .  "%N 1c''n Lap 1cu (analytic)")
      endif

C --- One-center expansion by brute force integration ---
C     Integration takes place on sphere R' = iax(2,ixi)
c     call fpiint(-60,0,np,xp,wp)
      call fpiint(-122,0,np,xp,wp)
C     Normalized spherical harmonics for all points on sphere R'
      call dmcpy(xp,1,3,xpt,np,1,np,3)
      call ropyln(np,xpt(1,1),xpt(1,2),xpt(1,3),lb,np,yl,rsq)

C     For each head channel, do
      call dpzero(vvaln,nlmy*nlf**2)
      call dpzero(vvalj,nlmy*nlf**2)
      call dpzero(dvalna,nlmy*nlf**2)
      call dpzero(dvalja,nlmy*nlf**2)
C     Integrate H^a and its one-center expansion on sphere R'
      if (cmdopt('--radius=',9,0,outs)) then
        j = 9
        call rxx(.not.a2bin(outs,radius,4,0,' ',j,len(outs)),
     .    'failed to parse '//outs)
        hcr = radius
        call info2(0,0,0,' Evaluate 1-center at radius %d',radius,0)
      else
        call info0(0,0,0,' Evaluate 1-center at hcr')
      endif
      do  il = 0, lb
        radius = hcr(il+1)
        do  ip = 1, np
          xmrp(1) =  radius*xp(1,ip)
          xmrp(2) =  radius*xp(2,ip)
          xmrp(3) =  radius*xp(3,ip)

          call dscal(3,1/avw,xmrp,1)  ! avw=1 unless OKA conventions are used
          call strck(nlf**2,iat,ixi,rtab,s_str%kaps,ikap,s_str%nkaps,nbasp,s_str%alph,halp,
     .      s_str%iax,s_str%npr,s,s_str%s,s_str%nitab,ips,nl,s_str%rsmh,s_str%c01g,xmrp,la,
     .      lhsm*1000+irot+loka+10,slk,slj,gslk,gslj)
C         Integration over the surface of sphere a_R' accomplished by
C         H^a_RL(x on sphere R') = sum_L' vval_L,L'(a_R') Y_L'(a_R')
C         vvaln holds H^a(x-R'); vvalj holds 1-center expansion to H
C         ilb is L'; ila is L; thus vval(ilb,ila) = vval_L,L'
          do  ilb = il*il+1, (il+1)**2
            xx = wp(ip)*yl(ip,ilb)
            do  ila = 1, nlma
              vvaln(ilb,ila) = vvaln(ilb,ila) + slk(ila,1)*xx
              vvalj(ilb,ila) = vvalj(ilb,ila) + slj(ila,1)*xx
              dvalna(ilb,ila) = dvalna(ilb,ila) + ddot(3,gslk(1,ila,1),1,xmrp,1)*xx/radius
              dvalja(ilb,ila) = dvalja(ilb,ila) + ddot(3,gslj(1,ila,1),1,xmrp,1)*xx/radius
            enddo
          enddo
        enddo
      enddo

c...deb
c --- Make screened 1c hankels in lmto3 way ---
C      itab = s_str%npr(iat)+ixi
C      call info5(10,0,0,' chkstr:  iat = %i, ixi = %i'//
C     .  ' corresponds to itab = %i',iat,ixi,itab,0,0)
C      nlmax = nlf**2
C      allocate(hvx(nlmax,nlmax,nttab))
C      allocate(hlx(nlmax,nlmax,nttab))
C      allocate(hsvx(nlmax,nlmax,nttab,1))
C      allocate(hslx(nlmax,nlmax,nttab,1))
CC ... unpack Clebsch-Gordan coefficients and sm. radii
C
C      call dscal(9,1/alat,plat,1)
C      call dscal(3*nbasp,1/alat,pos,1)
C
C      call hs1cx(nlf**2,alat,plat,pos,s_str%kaps(ikap),
C     .  nbasp,nlf,ips,w(orsmh),w(orsma),w(okmx),hcrs,
C     .  s_str%iax,s_str%npr,s,s_lat%cy,s_lat%cg,s_lat%indxcg,s_lat%jcg,
C     .  hsvx,hslx,hvx,hlx)
Cc       do  ip = 1, s_str%npr(nbasp+1)
Cc         print *,' chkstr: ip, hvx(1,1,ip) = ',ip, hvx(1,1,ip)
Cc       enddo
C
C      call dscal(9,alat,plat,1)
C      call dscal(3*nbasp,alat,pos,1)

C ... Calculate corresponding radial derivatives
      call dpzero(dvaln,nlmy*nlf**2)
      call dpzero(dvalj,nlmy*nlf**2)
      do  il = 0, lb
        radius = hcr(il+1)
        do  ip = 1, np
          xmrp(1) =  (radius + .0005d0)*xp(1,ip)
          xmrp(2) =  (radius + .0005d0)*xp(2,ip)
          xmrp(3) =  (radius + .0005d0)*xp(3,ip)
          call dscal(3,1/avw,xmrp,1) ! avw=1 unless OKA conventions are used
          i = lhsm*1000+irot+loka
          call strck(nlf**2,iat,ixi,rtab,s_str%kaps,ikap,s_str%nkaps,
     .      nbasp,s_str%alph,halp,s_str%iax,s_str%npr,s,s_str%s,s_str%nitab,
     .      ips,nl,s_str%rsmh,s_str%c01g,xmrp,la,i,slk,slj,gslk,gslj)
          xmrp(1) =  (radius - .0005d0)*xp(1,ip)
          xmrp(2) =  (radius - .0005d0)*xp(2,ip)
          xmrp(3) =  (radius - .0005d0)*xp(3,ip)
          call dscal(3,1/avw,xmrp,1) ! avw=1 unless OKA conventions are used
          call strck(nlf**2,iat,ixi,rtab,s_str%kaps,ikap,s_str%nkaps,
     .      nbasp,s_str%alph,halp,s_str%iax,s_str%npr,s,s_str%s,s_str%nitab,
     .      ips,nl,s_str%rsmh,s_str%c01g,xmrp,la,i,slk(1,2),slj(1,2),gslk,gslj)
          do  ilb = il*il+1, (il+1)**2
            do  ila = 1, nlma
              dvaln(ilb,ila) = dvaln(ilb,ila) + (slk(ila,1)-slk(ila,2))/.001d0*wp(ip)*yl(ip,ilb)
              dvalj(ilb,ila) = dvalj(ilb,ila) + (slj(ila,1)-slj(ila,2))/.001d0*wp(ip)*yl(ip,ilb)
            enddo
          enddo
        enddo
      enddo

C ... make H0 = bare, unsmoothed Hankel
      if (ixi == 1 .and. irot == 0) then
        do  il = 1, lb+1
          radius = hcr(il)
          call besslr(s_str%kaps(ikap)*(radius/avw)**2,loka,0,il-1,fi,gi0)
          gi(il) = gi0(il)/(radius/avw)**il
        enddo
      endif

C --- One-center expansion by P_kL expansion ---
      if (lhsm == 1) then
C        v1chsp = 0; v1chs = 0
C        call strg1c(3,1,nlf**2,nlmy,nl,iat,ixi,ips,s_str%rsmh,s_str%c01g,
C     .    ehl,hcr+1d-3,kmxt,rsma,alat,plat/alat,s_lat%pos,s_str%iax,
C     .    s_str%npr,ntabg,s_lat%cy,s_lat%cg,s_lat%indxcg,s_lat%jcg,s,
C     .    v1chsp,dv1chs,v1clhs,dv1clhs)
C        call strg1c(3,1,nlf**2,nlmy,nl,iat,ixi,ips,s_str%rsmh,s_str%c01g,
C     .    ehl,hcr-1d-3,kmxt,rsma,alat,plat/alat,s_lat%pos,s_str%iax,
C     .    s_str%npr,ntabg,s_lat%cy,s_lat%cg,s_lat%indxcg,s_lat%jcg,s,
C     .    v1chs,dv1chs,v1clhs,dv1clhs)
C        v1chsp = (v1chsp - v1chs)/2d-3
        call info2(10,0,0,' ... calling strg1c with ivl=3 [envelope = Hsm + c01(0)*G0 + c01(1)*G1]'//
     .    ' with ixi=%i',ixi,2)
        call strg1c(3,1,nlf**2,nlmy,nl,iat,ixi,ips,s_str%rsmh,s_str%c01g,ehl,hcr,
     .    kmxt,rsma,alat,plat,s_lat%pos,s_str%iax,s_str%npr,ntabg,
     .    s_lat%cy,s_lat%cg,s_lat%indxcg,s_lat%jcg,s,v1chs,dv1chs,
     .    v1clhs,dv1clhs)
      endif

C ... Printout
      do  ila = 1, nlma
        call info(10,1,0,' Compare numerical YL expns and radial '//
     .    'deriv of H, 1C for ilm=%i, ikap=%i',ila,ikap)
        if (ixi == 1 .and. irot == 0)
     .  write(stdo,928) 'H0  ',(0d0,ilb=1,ila-1),gi(1+ll(ila)),(0d0,ilb=ila+1,min(nlmy,nlmb))
        write(stdo,928) 'Ha  ',(vvaln(ilb,ila), ilb = 1,min(nlmy,nlmb))
        if (lhsm == 1) then
          write(stdo,928) '1cs ',(v1chs(ilb,ila,1), ilb=1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(vvaln(ilb,ila)-v1chs(ilb,ila,1),ilb = 1, min(nlmy,nlmb))
          write(stdo,928) '1cu ',(vvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
        else
          write(stdo,928) '1c  ',(vvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(vvaln(ilb,ila)-vvalj(ilb,ila),ilb = 1, min(nlmy,nlmb))
        endif

        write(stdo,928) "H'n ",(dvaln(ilb,ila), ilb = 1,min(nlmy,nlmb))
        if (loka == 0) then
          write(stdo,928) "H'a ",(dvalna(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(dvalna(ilb,ila)-dvaln(ilb,ila),ilb = 1, min(nlmy,nlmb))
        endif

        if (loka == 0 .and. lhsm == 0) then
          write(stdo,928) "1c'a",(dvalja(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(dvalna(ilb,ila)-dvalja(ilb,ila),ilb = 1, min(nlmy,nlmb))
          write(stdo,928) "1c'n",(dvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(dvalja(ilb,ila)-dvalj(ilb,ila),ilb = 1, min(nlmy,nlmb))
        elseif (loka == 0) then
          write(stdo,928) "1c's",(dv1chs(ilb,ila,1), ilb=1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(dvalna(ilb,ila)-dv1chs(ilb,ila,1),ilb = 1, min(nlmy,nlmb))
          write(stdo,928) "1c'n",(dvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) "1c'u",(dvalja(ilb,ila), ilb = 1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(dvalj(ilb,ila)-dvalja(ilb,ila),ilb = 1, min(nlmy,nlmb))
        endif

        if (loka == 0 .and. lhsm /= 0) then
          lvalj(1:nlmb,1:nlma) = -s_str%kaps(ikap)*vvalj(1:nlmb,1:nlma)
          write(stdo,931) "1c''s",(v1clhs(ilb,ila,1), ilb=1,min(nlmy,nlmb))
C         write(stdo,931) "1c''u",(-s_str%kaps(ikap)*vvalj(ilb,ila), ilb=1,min(nlmy,nlmb))
          write(stdo,931) "1c''u",(lvalj(ilb,ila), ilb=1,min(nlmy,nlmb))
          write(stdo,928) '-   ',(v1clhs(ilb,ila,1)-(lvalj(ilb,ila)),ilb = 1,min(nlmy,nlmb))
        endif

      enddo
  928 format(a,1x,25F9.5)
  931 format(a,25F9.5)

      call info0(10,1,0,' Summary')
      if (lhsm == 1) then
        call prterr(nlmy,nlmb,nlma,v1chs,vvaln,err,ier1,ier2,xx,ixx,ixx)
      else
        call prterr(nlmy,nlmb,nlma,vvalj,vvaln,err,ier1,ier2,xx,ixx,ixx)
      endif
      call info5(10,0,0,' Ha vs 1c:%16p maximum difference %,3;3g'//
     .  ' at ilm = %i, ilm'' = %i',err,ier2,ier1,0,0)

      if (loka == 0) then
        call prterr(nlmy,nlmb,nlma,dvalna,dvaln,err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," H'n vs H'a:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)
      endif

      if (loka == 0 .and. lhsm == 0) then
        call prterr(nlmy,nlmb,nlma,dvalna,dvalja,err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," H'a vs 1c'a:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)
        call prterr(nlmy,nlmb,nlma,dvalj,dvalja,
     .    err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," 1c'n vs 1c'a:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)
      elseif (loka == 0) then
        call prterr(nlmy,nlmb,nlma,dvalna,dv1chs,err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," H'a vs 1c's:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)
        call prterr(nlmy,nlmb,nlma,dvalj,dvalja,err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," 1c'n vs 1c'u:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)

        call prterr(nlmy,nlmb,nlma,v1clhs,lvalj,err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0," 1c''n vs 1c''u:%16p maximum difference %,3;3g"//
     .    " at ilm = %i, ilm' = %i",err,ier2,ier1,0,0)

      endif

c...deb
C      deallocate(hsvx,hslx,hvx,hlx)
c...deb

      call getqu(ltmp)
      if (.not. ltmp) call fexit(0,0,outs,0)
      goto 23

C --- Plot screened Hankel along a specified direction vector ---
   70 continue

      call info2(0,1,0,' 3. Plot line from x =%3:1,d to x =%3:1,d (a.u.)',v1,v2)
      if (irot == 100 .and. s_str%nkaps == 1)  print 924
      if (irot == 100 .and. s_str%nkaps == 2) then
        call info2(10,0,0,'%4pGenerating %?#(n==1)#value#slope# function',ikap,0)
      endif

      call query('number of points=',2,nvec)
      if (nvec > nvmax) call rxi('number of points cannot exceed max=',nvmax)
      iat = 1
      if (cmdopt('--atom=',7,0,outs)) then
        j = 7
        call rxx(.not.a2bin(outs,iat,2,0,' ',j,len(outs)),'failed to parse '//outs)
      endif
      call query('atom=',2,iat)
      if (iat > nbasp) call rxi('atom cannot be larger than nbas=',nbasp)
      ixi = 1
      offR  = s_str%npr(iat)
      nlmb = s_str%iax(niax*offR+9)
      lb = ll(nlmb)
      hcr(1:1+lb) = s_spec(ips(iat))%hcr(1:1+lb)
      call info5(10,1,0,' Atom %i has hcr=%n:1;4,4d',iat,lb+1,hcr,0,0)

      do  i = 1, nvec
        xmrp(1) = v1(1) + dble(i-1)/max(nvec-1,1)*(v2(1)-v1(1))
        xmrp(2) = v1(2) + dble(i-1)/max(nvec-1,1)*(v2(2)-v1(2))
        xmrp(3) = v1(3) + dble(i-1)/max(nvec-1,1)*(v2(3)-v1(3))
        call dscal(3,1/avw,xmrp,1) ! avw=1 unless OKA conventions are used
        call strck(nlf**2,iat,ixi,rtab,s_str%kaps,ikap,s_str%nkaps,nbasp,
     .    s_str%alph,halp,s_str%iax,s_str%npr,s,s_str%s,s_str%nitab,ips,nl,
     .    s_str%rsmh,s_str%c01g,xmrp,nlf-1,lhsm*1000+irot+loka,slk,slj,gslk,gslj)
        reslin(i,1) = avw*xmrp(1)
        reslin(i,2) = avw*xmrp(2)
        reslin(i,3) = avw*xmrp(3)
        do  ilb = 1, min(nlmb,9)
          reslin(i,3+ilb) = min(max(slk(ilb,1),-100d0),100d0)
        enddo
      enddo

      call info2(10,0,0,' Writing line data, x,y,z + Ha(1:%i)  to file "plot" ',nlmb,0)
      ifi = fopn('PLOT')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%6,6g)',reslin,0,nvmax,nvec,min(3+nlmb,12))
      call fclose(ifi)
      deallocate(s,halp)
C     call fexit(0,0,outs,0)

c ... write line data to file 76 for gnuplot
c     do  i = 1, nvec
c       write(76,'(1x,3f12.6,1x,18g14.5)')
c    .    (reslin(i,ilb),ilb=1, 3+min(nlmb,9))
c     enddo

      end
