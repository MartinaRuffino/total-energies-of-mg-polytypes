      subroutine trothf(s_bz,s_lat,s_site,s_spec,s_ham,
     .  nl,nbas,isp,nsp,offH,istab,qp,ldim,smpot,vconst)
C- Test roth for FP hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc lshft ipq star
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat npgrp qlat pos nabc
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvtrob smhsbl
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvtrob smhsbl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pvtrob smhsbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  pvtrob
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   qp    :irreducible k-points
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst
Ci   g     :point group operations
Ci   ag    :translation part of space group
Co Outputs
Cl Local variables
Cl  lbloch :0, Bloch transform by FFT
Cl         :1, standard FT
Cl         :2, both
Cl   idim  :dimension of block of orbitals to be folded down (makidx.f)
Cb Bugs
Cb   should use s_bz%qp, s_lat%istab
Cr Remarks
Cr   Cr3Si6 test is taken from the FP standard test library (9 atoms, 2 classes)
Cr    *roth test by:
Cr     lmf cr3si6 --pr35 -vnit=1 --rs=0,0 --no-iactiv
Cr    *rotwf test by:
Cr     lmf cr3si6 --pr35 -vnit=1 --rs=0,0 --no-iactiv --evec
Cr     Note the significant deviation from (z^ H z - E) = 0:
Cr       global max error = 1.1d-4  with k-point avg = 1.1E-05
Cr     Tighten the Ewald sum tolerance and it should reduce to
Cr       global max error = 1.6E-05 with k-point avg = 3.2E-06
Cr   MnN check demonstrates additional features:
Cr   offset BZ mesh, extra group ops from TR symmetry, multiple kappa basis
Cr    *roth test by:
Cr     cp gf/test/ctrl.mnn .
Cr     lmf mnn -vlmf=t -opt=0 --no-iactiv
Cr     lmf mnn -vlmf=t -opt=1 --no-iactiv (no different from opt=0?)
Cr     lmf mnn -vlmf=t -opt=3 --no-iactiv
Cr    *rotwf test by:
Cr     lmf mnn -vlmf=t -opt=0 --no-iactiv --evec
Cr     You should get:
Cr       max error encountered = 3.5E-06 with k-point avg = 9.6E-07
Cu Updates
Cu   10 Aug 15  Re-designed for v7.11 with f90 pointers and structures
Cu   14 Jun 02  Adapted from troth.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nbas,nl,isp,nsp,ldim,offH(n0H,nkap0,nbas+1),istab(nbas,*)
      double precision qp(3,*),smpot(*),vconst
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      complex(8), allocatable :: wki2(:)
      complex(8),allocatable:: sq(:,:,:,:),srot(:,:)
      real(8), pointer :: pos(:,:),g(:,:,:)
      complex(8), pointer :: wki(:,:),wk2z(:,:)
C ... Local parameters
      logical latvec,cmdopt
      logical levec ! If true, rotated evecs rather than h
      logical lovl  ! If true, rotated overlap rather than h
      integer lshft(3)
      integer iat,i1,i2,i3,ib,ibla,iblb,idim,ig,ipq1,iq1,isw,jj1,jj2,
     .  jj3,k,kcplx,ldima,ldimb,ldimbs,mxorb,nev,nglob,niax,nl2,nlma,
     .  nlmb,nqp,opt,roth,npgrp,ifac(3),nkabc(3)
      double precision plat(3,3),qlat(3,3),qb(3,3),gerrmx,erravg,qoff(3)

      parameter (niax=10)
      character*80 strn
      integer :: napw=0 !, igapw(1)
      integer nlx
      parameter (nlx=5)
      double precision evl(ldim),rmat(nlx**5*2)
      double precision q1(3),q0(3),qk,dq(3),vavg
      double precision errmx
      double complex sq1(ldim,ldim)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)

C ... External calls
      external awrit7,bzmsh00,cwrite,dcopy,dpzero,dscal,dvset,fexit3,
     .  fftz30,hambl,info0,info2,info5,phmbls,poppr,prtev,pshpr,pvtro4,
     .  pvtro5,pvtrob,pvtroc,pvtroh,pvtroz,ran1in,rotwf,rx1,shorbz,to,
     .  yprm,zcopy,zhevx,zmcpy,zprm3,zqinv,ztoyy

      nsp = nglob('nsp')
      plat = s_lat%plat
      npgrp = s_lat%npgrp
      nqp = s_bz%nkp
      nkabc = s_bz%nkabc
      lshft = s_bz%lshft
      qlat = s_lat%qlat
      pos => s_lat%pos
      idim = 0  ! FP has no iwaves yet

C ... Make ifac and qb
      call pshpr(1)
      call bzmsh00(plat,lshft,0,nkabc,ifac,qb)
      call poppr
      allocate(g(3,3,npgrp))
      call dcopy(9*npgrp,s_lat%symgr,1,g,1)

      call ran1in(12)
      allocate(sq(ldim,ldim,nqp,2))
      if (napw /= 0) call rx('trothf not ready for apw')
      mxorb = nglob('mxorb')
      opt = 0
      if (cmdopt('-opt=1',6,0,strn) .or. cmdopt('-opt=5',6,0,strn)) then
        opt = 1
        if (cmdopt('-opt=5',6,0,strn)) opt = 5
        print "(/' trothf: testing for opt=',i4)", opt
      else if (cmdopt('-opt=3',6,0,strn)) then
        opt = 3
        iat = 3
        print "(/' trothf: testing for opt=',i4,' iat=',i4)", opt,iat
      else if (cmdopt('-opt=4',6,0,strn)) then
        opt = 4
        iat = 3
        print "(/' trothf: testing for opt=',i4,' iat=',i4)", opt,iat
      else
        print "(/' trothf: testing for opt=',i4)", opt
      endif

      levec = cmdopt('--evec',6,0,strn)
      lovl  = cmdopt('--ovl',5,0,strn)

      ldima = ldim
      ldimb = ldim

      qoff(1) = qk(1,1,1,1)
      qoff(2) = qk(2,1,1,1)
      qoff(3) = qk(3,1,1,1)

C ... complex storage mode
      kcplx = 1
C     if (cmdopt('-kcplx=2',8,0,strn)) kcplx=2
C     if (cmdopt('-kcplx=1',8,0,strn)) kcplx=1
      opt = opt + 100*kcplx + 40  ! 40 for FP phase convention

      call info5(1,1,1,' ... make '//
     .  '*Hamiltonian*'//
     .  '%?#n#%10p*Overlap*##'//
     .  '%?#n#%10p*Eigenvectors*##'//
     .  ' for %i irreducible qp'//
     .  '%?#n==2# and 2 spins## ...',
     .  isw(cmdopt('--ovl',5,0,strn)),isw(cmdopt('--evec',6,0,strn)),nqp,nsp,0)

      do  isp = 1, nsp
      do  iq1 = 1, nqp

C ... hamiltonian
      call dvset(sq(1,1,iq1,isp),1,2*ldim**2,0d0)
      vavg = vconst

      call dcopy(3,qp(1,iq1),1,q1,1)
C     call shorbz(qp(1,iq1),q1,qlat,plat)

      call pvtrob(isp,ldim,s_lat,s_site,s_spec,s_ham,smpot,
     .  vconst,levec,lovl,q1,sq(1,1,iq1,isp),evl)

      enddo
      enddo

      print *, 'done making H for irr q and spin'
C     call zprm('s(7)',2,sq(1,1,7,1),ldim,ldim,ldim)

      call info2(1,1,1,' ... check '//
     .  '*roth*%?#n#%11p*rotwf*##'//
     .  ' for all points in the irreducible BZ ...',
     .  isw(cmdopt('--evec',6,0,strn)),0)

      allocate(wk2z(ldim,ldim),srot(ldim,ldim))
      nl2 = nl*nl
      iq1 = 0
      isp = 1
      gerrmx = 0; erravg = 0
      do  i3 = 1, s_bz%nkabc(3)
      do  i2 = 1, s_bz%nkabc(2)
      do  i1 = 1, s_bz%nkabc(1)

        iq1 = iq1+1

        ipq1 = s_bz%ipq(iq1)
        ig = s_bz%star(iq1+1)

C   ... q into which h or z is rotated
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

C   ... q from which h or z is rotated
        q0(1) = g(1,1,ig)*q1(1) + g(2,1,ig)*q1(2) + g(3,1,ig)*q1(3)
        q0(2) = g(1,2,ig)*q1(1) + g(2,2,ig)*q1(2) + g(3,2,ig)*q1(3)
        q0(3) = g(1,3,ig)*q1(1) + g(2,3,ig)*q1(2) + g(3,3,ig)*q1(3)

C   ... Confirm that  g*q1 differs only a a G vector from q(ipq1)
        dq(1) = q0(1)-qp(1,ipq1)
        dq(2) = q0(2)-qp(2,ipq1)
        dq(3) = q0(3)-qp(3,ipq1)
        if (.not. latvec(1,1d-6,plat,dq)) then
          print *, q0
          print *, q1
          print *, qp(1,ipq1),qp(2,ipq1),qp(3,ipq1)
          print *, iq1,ipq1,ig
          call rx('trothf: bad ipq')
        endif

C   ... Rotate H from q0 to q1, store in srot
        call dcopy(ldima*ldimb*2,sq(1,1,ipq1,1),1,srot,1)

C       Check to make sure roth uses proper amount of memory
        if (mod(opt,10) >= 3 .and. mod(opt,10) <= 4) then
          ib = iat
          ibla = max(mod(opt/1000,10),1)
          iblb = max(mod(opt/10000,10),1)
C         For roth
          k = (ldima * (offH(iblb,1,ib+1)-offH(iblb,1,ib)) +
     .        (offH(ibla,1,ib+1)-offH(ibla,1,ib)) * ldimb)
C    .      - 2
          if (ibla == 4) k =
     .      2*ldima * (offH(ibla,1,ib+1)-offH(ibla,1,ib))
C    .      - 2
          if (mod(opt,10) == 4) k = 2*k
C         call wkprnt(1)
C         call defcc(owki,k)
          allocate(wki(k,1))
C         check that roth properly zero's subblock
          call dvset(wki,1,2*k,-99d0)
C         save srot, to make sure it is preserved
          if (mod(opt,10) == 4) then
C           call defcc(owki2,ldima*ldimb)
            allocate(wki2(ldima*ldimb))
            call dcopy(ldima*ldimb*2,srot,1,wki2,1)
          endif
        else
          wki => wk2z
        endif
        ldimbs = ldimb
        if (mod(opt,10) == 5) then
          call pvtro5(kcplx,ldima,srot,srot)
          ldimb = 1
        endif

        if (levec) then
          nev = ldim
          call rotwf(140,nl,nbas,1,s_lat%pos,s_ham%offH,s_ham%iprmb,
     .      s_lat%istab(1,ig),s_lat%symgr(1,ig),s_lat%ag(1,ig),
     .      q1,rmat,0,ldim,ldim,nev,sq(1,1,ipq1,1),srot)
        else
          if (roth(opt,nl,nbas,pos,iat,offH,s_ham%iprmb,istab(1,ig),
     .    s_lat%symgr(1,ig),s_lat%ag(1,ig),q1,ldima,ldimb,wki,srot) < 0)
     .    call rx1('roth: rotation+inv for mode %i not implemented',opt)
        endif

C   ... H or S at rotated q into sq1
        call shorbz(q1,q1,qlat,plat)
        call pvtrob(isp,ldim,s_lat,s_site,s_spec,s_ham,smpot,
     .    vconst,.false.,lovl,q1,sq1,evl)

        if (mod(opt,10) == 5) then
          call pvtro5(kcplx,ldima,sq1,sq1)
          ldimb = 1
        endif
        if (levec) then
          call pshpr(50)
          call pvtroz(ldim,iq1,ipq1,ig,sq1,evl,srot,errmx)
          call poppr
          gerrmx = max(gerrmx,errmx)
          erravg = erravg + errmx/(s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3))
        elseif (mod(opt,10) == 4) then
          nlma = offH(ibla,1,ib+1)-offH(ibla,1,ib)
          nlmb = offH(iblb,1,ib+1)-offH(iblb,1,ib)
          if (ibla == 4) nlmb = nlma
          call pvtro4(opt,nbas,iat,offH,ldima,ldimb,iq1,ipq1,ig,
     .      wki2,srot,sq1,nlma,nlmb,wki,wki)
        else
          call pshpr(50)
          call pvtroh(opt,nbas,iat,offH,ldima,ldimb,iq1,ipq1,ig,
     .      sq(1,1,ipq1,1),srot,sq1,errmx)
          gerrmx = max(gerrmx,errmx)
          erravg = erravg + errmx/(s_bz%nkabc(1)*s_bz%nkabc(2)*s_bz%nkabc(3))
          call poppr
        endif
        if (mod(opt,10) >= 3) then
        endif
        ldimb = ldimbs

      enddo
      enddo
      enddo
      deallocate(srot,wk2z)
      if (mod(opt,10) == 4) then
        deallocate(wki2)
      endif

      print 334, gerrmx,erravg
  334 format(/' trothf: max error encountered :',1pe8.1,3x,'k-point avg of maxerr :',e8.1)

      call fexit3(0,111,' Exit 0 finished trothk,'//
     .  ' opt=%i, ldim=%i, kcplx=%i',opt,ldim,kcplx)
      end

      subroutine pvtro4(opt,nbas,iat,offH,ldima,ldimb,iq,ipq,ig,
     .  su,s1,s2,nlma,nlmb,srow,scol)
C- Check for opt=4 case
      implicit none
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer opt,iat,offh(n0H,nkap0,1),ldima,ldimb,iq,ipq,ig,
     .  nbas,nlma,nlmb
      double precision su(ldima,ldimb,2),s1(ldima,ldimb,2),
     .  s2(ldima,ldimb,2),srow(nlma,ldimb,2),scol(ldima,nlmb,2)
      double precision tol
      integer i1,i2,j1,j2,opt0,ibla,off2,iprint,i3,i4,j3,j4
      integer i,j,kcplx,ioff,joff,iblb

C     to look at intermediate blocks (ibla=4)
      i3 = 0
      i4 = -1
      j3 = 0
      j4 = -1

      tol = 0
      do  50  j = 1, ldimb
      do  50  i = 1, ldima
        if (abs(su(i,j,1)-s1(i,j,1)) > tol .or.
     .      abs(su(i,j,2)-s1(i,j,2)) > tol) then
          print 332, iq,ipq,ig,i,j
  332     format('trothf: opt=4 but hamiltonian changed iq,ipq=',2i5,
     .      '  ig=',i2,'  i,j=',2i3)
          call rx('oops')
        endif
   50 continue

      kcplx = mod(opt/100,10)
      if (kcplx /= 0) then
        call ztoyy(srow,nlma,ldimb,nlma,ldimb,kcplx,0)
        off2 = nlma*ldimb*2
        call ztoyy(scol(1+off2,1,1),ldima,nlmb,ldima,nlmb,kcplx,0)
        call ztoyy(s2,ldima,ldimb,ldima,ldimb,kcplx,0)
      endif

      opt0 = mod(opt,10)
      kcplx = mod(opt/100,10)
      if (opt0 /= 4) call rx('should not call pvtro4 for this opt')
      ibla = max(mod(opt/1000,10),1)
      iblb = max(mod(opt/10000,10),1)
      tol = 1d-6

      i1 = 1+offH(ibla,1,iat)
      i2 = offH(ibla,1,iat+1)
      if (ibla == 4) then
        i1 = 1+offH(iblb,1,iat)
        i2 = offH(iblb,1,iat+1)
C       i3,i4 are i block (a second pass is made for them)
        i3 = offH(iblb,1,nbas+1) + 1+offH(2,1,iat)
        i4 = offH(iblb,1,nbas+1) + offH(2,1,iat+1)
      endif
      ioff = i1
      j1 = 1
      j2 = ldimb
      j3 = 1
      j4 = ldimb

      if (i2-i1+1+i4-i3+1 /= nlma) call rx('bad call to pvtro4')
      if (j2-j1+1 /= ldimb) call rx('bad call to pvtro4')

    6 if (iprint() >= 50)
     .  call awrit7('checking iq=%i  ig=%i  ipq=%i  range=%i %i %i %i',
     .  ' ',80,6,iq,ig,ipq,i1,i2,j1,j2)


      do  10  i  = i1, i2
      do  10  j  = j1, j2

        if (abs(srow(i-ioff+1,j-j1+1,1)-s2(i,j,1)) > tol .or.
     .      abs(srow(i-ioff+1,j-j1+1,2)-s2(i,j,2)) > tol) then

          print 333, iq,ipq,ig,i,j
  333     format(' trothf: mismatch iq,ipq=',2i5,' ig=',i2,'  i,j=',2i3)
          call yprm('unrotated H',2,su,ldima,ldima,ldimb,ldimb)
          call yprm('rotated H',2,s1,ldima,ldima,ldimb,ldimb)
          call yprm('H (rotated q)',2,s2,ldima,ldima,ldimb,ldimb)
          return
        endif
   10 continue
C     Make another pass for intermediate waves
      if (i3 /= 0 .and. j3 /= 0) then
        ioff = i3 - (i2-i1+1)
        i1 = i3
        i2 = i4
        j1 = j3
        j2 = j4
        i3 = 0
        j3 = 0
        goto 6
      endif


      i1 = 1
      i2 = offH(ibla,1,nbas+1)
      j1 = 1+offH(iblb,1,iat)
      j2 = offH(iblb,1,iat+1)

C     to look at intermediate blocks (ibla=4)
      i3 = 0
      i4 = -1
      j3 = 0
      j4 = -1

      off2 = nlma*ldimb*2
      joff = j1
      if (ibla == 4) then
        i3 = i1
        i4 = i2
        j3 = offH(1,1,nbas+1) + 1+offH(2,1,iat)
        j4 = offH(1,1,nbas+1) + offH(2,1,iat+1)
      endif
      if (i2-i1+1 /= ldima) call rx('bad call to pvtro4')
      if (j2-j1+1+j4-j3+1 /= nlmb) call rx('bad call to pvtro4')

  106 if (iprint() >= 50)
     .  call awrit7('checking iq=%i  ig=%i  ipq=%i  range=%i %i %i %i',
     .  ' ',80,6,iq,ig,ipq,i1,i2,j1,j2)

      do  110  i  = i1, i2
      do  110  j  = j1, j2

        if (abs(scol(i-i1+off2+1,j-joff+1,1)-s2(i,j,1)) > tol .or.
     .      abs(scol(i-i1+off2+1,j-joff+1,2)-s2(i,j,2)) > tol) then

          print 333, iq,ipq,ig,i,j
          call yprm('unrotated H',2,su,ldima,ldima,ldimb,ldimb)
          call yprm('rotated H',2,s1,ldima,ldima,ldimb,ldimb)
          call yprm('H (rotated q)',2,s2,ldima,ldima,ldimb,ldimb)
          return
        endif
  110 continue
C     Make another pass for intermediate waves
      if (i3 /= 0 .and. j3 /= 0) then
        joff = j3 - (j2-j1+1)
        i1 = i3
        i2 = i4
        j1 = j3
        j2 = j4
        i3 = 0
        j3 = 0
        goto 106
      endif


      if (kcplx /= 0) then
        call ztoyy(srow,nlma,ldimb,nlma,ldimb,0,kcplx)
        off2 = nlma*ldimb*2
        call ztoyy(scol(1+off2,1,1),ldima,nlmb,ldima,nlmb,0,kcplx)
        call ztoyy(s2,ldima,ldimb,ldima,ldimb,0,kcplx)
      endif

      end
      subroutine pvtroh(opt,nbas,iat,offH,ldima,ldimb,iq,ipq,ig,
     .  su,s1,s2,errmx)
C- Compares rotated and unrotated s matrix, analyzes errors
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt
Ci   nbas  :size of basis
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   ldima
Ci   ldimb
Ci   iq
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   ig
Ci   su    :matrix at unrotated qp
Ci   s1    :rotated h
Ci   s2    :h at rotated qp
Co Outputs
Co   errmx :
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   29 Jan 01
C ----------------------------------------------------------------------
      implicit none
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer opt,iat,offh(n0H,nkap0,1),ldima,ldimb,iq,ipq,ig,nbas
      double precision su(ldima,ldimb,2),s1(ldima,ldimb,2),
     .  s2(ldima,ldimb,2)
      double precision tol,err,errmx
      integer i1,i2,j1,j2,opt0,ibla,iblb,ib,iopt0,iprint
      integer i3,i4,j3,j4,awrite
      logical pass1
      integer i,j,kcplx,recl
      parameter (recl=512)
      character*(recl) lstr

      opt0 = mod(opt,10)
      kcplx = mod(opt/100,10)
      errmx = 0

      if (kcplx /= 0) then
        call ztoyy(s1,ldima,ldimb,ldima,ldimb,kcplx,0)
        call ztoyy(s2,ldima,ldimb,ldima,ldimb,kcplx,0)
      endif

C      if (ig == 12) then
C        print *, 'ig,iq=',ig,iq
C        print '(2f15.10)', s1(1,10,1),s1(1,10,2)
C        print '(2f15.10)', s2(1,10,1),s2(1,10,2)
C        stop
C      endif

C     to look at intermediate blocks (ibla=4)
      i3 = 0
      i4 = 0
      j3 = 0
      j4 = 0

      if (opt0 == 4) opt0 = 3
      ibla = max(mod(opt/1000,10),1)
      iblb = max(mod(opt/10000,10),1)
      tol = 1d-6
      if (opt0 == 5) tol = .02d0
      pass1 = .true.
    5 continue
      if (opt0 == 3) then
        if (pass1) then
          i1 = 1
          i2 = offH(ibla,1,nbas+1)
          j1 = 1+offH(iblb,1,iat)
          j2 = offH(iblb,1,iat+1)
          if (ibla == 4) then
            i3 = i1
            i4 = i2
            j3 = offH(iblb,1,nbas+1) + 1+offH(2,1,iat)
            j4 = offH(iblb,1,nbas+1) + offH(2,1,iat+1)
          endif
          pass1 = .false.
        else
          i1 = 1+offH(ibla,1,iat)
          i2 = offH(ibla,1,iat+1)
          if (ibla == 4) then
            i1 = 1+offH(iblb,1,iat)
            i2 = offH(iblb,1,iat+1)
            i3 = offH(iblb,1,nbas+1) + 1+offH(2,1,iat)
            i4 = offH(iblb,1,nbas+1) + offH(2,1,iat+1)
          endif
          j1 = 1
          j2 = ldimb
          j3 = 1
          j4 = ldimb
          pass1 = .true.
        endif
      else
        i1 = 1
        i2 = offH(ibla,1,nbas+1)
        j1 = 1
        j2 = offH(iblb,1,nbas+1)
        if (ibla == 4) then
          j1 = i1
          j2 = i2
        endif
      endif

    6 if (iprint() >= 50) then
        lstr = ' '
C       if (.not. pass1) print *, ' '
        i = awrite(' pvtroh: check |R H R^-1 - H^R| for iq=%i  ig=%i  ipq=%i  range=%i %i %i %i',
     .  lstr,512,0,iq,ig,ipq,i1,i2,j1,j2,0)
        call cwrite(lstr,0,i-1,0)
      endif
C     A dummy for debugging breaks
      if (ig /= 1) then
        iopt0 = 0
      endif

      iblb = max(mod(opt/10000,10),1)
      iopt0 = 0
c     re-entry for opt0 = 1, downfolding
   19 continue
      do  10  i  = i1, i2
        if (opt0 == 1 .or. opt0 == 5) then
C         Find basis atom corresonding to site, make sure diag
C         NB this doesn't work for sil
          do  12  ib = 1, nbas
            if (i <= offH(iblb,1,ib)+iopt0) goto 12
            if (i > offH(iblb,1,ib+1)+iopt0) goto 12
            j1 = offH(iblb,1,ib)+1 +iopt0
            j2 = offH(iblb,1,ib+1) +iopt0
            goto 14
   12     continue
          goto 10
   14     continue
        endif

        if (opt0 == 5) then
          j1 = 1
          j2 = 1
        endif
        do  j  = j1, j2

          errmx = max(abs(s1(i,j,1)-s2(i,j,1)),
     .                abs(s1(i,j,2)-s2(i,j,2)),
     .                errmx)
     .

        if (abs(s1(i,j,1)-s2(i,j,1)) > tol .or.
     .      abs(s1(i,j,2)-s2(i,j,2)) > tol) then

          err = max(abs(s1(i,j,1)-s2(i,j,1)),
     .              abs(s1(i,j,2)-s2(i,j,2)))
     .
          print 333, iq,ipq,ig,i,j,err
  333     format(/' trothk: mismatch iq,ipq=',2i5,'  ig=',i2,'  i,j=',
     .      2i3,'  err=',1pe8.1)
          if (opt0 /= 5)
     .    call yprm('unrotated H',2+kcplx,su,ldima*ldimb,ldima,ldima,ldimb)
          call yprm('rotated H',2,s1,ldima*ldimb,ldima,ldima,ldimb)
          call yprm('H (rotated q)',2,s2,ldima*ldimb,ldima,ldima,ldimb)
          return
        endif
      enddo
   10 continue
C     Make another pass for intermediate waves
      if (i3 /= 0 .and. j3 /= 0) then
        i1 = i3
        i2 = i4
        j1 = j3
        j2 = j4
        i3 = 0
        j3 = 0
C       pass1 = .not. pass1
        goto 6
      endif

C     Make second pass for opt0=1 (site-diagonal) mode
      if (iblb == 1 .and. (opt0 == 1 .or. opt0 == 5)) then
        iblb = 2
        iopt0 = offH(1,1,nbas+1)
        goto 19
      endif
C     Make second pass for opt0=3 mode
      if (.not. pass1) goto 5

      if (kcplx /= 0) then
        call ztoyy(s1,ldima,ldimb,ldima,ldimb,0,kcplx)
        call ztoyy(s2,ldima,ldimb,ldima,ldimb,0,kcplx)
      endif

      if (iprint() >= 50) then
        print 334, errmx
  334   format(' max err=',1pe8.1)
      endif

      end

      subroutine pvtroz(ldim,iq,ipq,ig,ham,evl,z,errmx)
C- Check that z+ H z - evl is small
      implicit none
      integer ldim,iq,ipq,ig
      double precision evl(ldim),errmx
      double complex ham(ldim,ldim),z(ldim,ldim)
C Local
      double complex wk(ldim,ldim),zhz(ldim,ldim)
      double precision tol,err
      integer i,j,iprint,recl,awrite
      parameter (recl=512)
      character*(recl) lstr

      tol = 1d-5

      call phmbls(1,ldim,ldim,evl,[0],wk,ham,z,z,zhz)

      do  i  = 1, ldim
        zhz(i,i) = zhz(i,i) - evl(i)
      enddo

      if (iprint() >= 50) then
        lstr = ' '
C       if (.not. pass1) print *, ' '
        i = awrite(' pvtroz: check |z+ H z - evl| for iq=%i  ig=%i  ipq=%i',
     .  lstr,512,0,iq,ig,ipq,0,0,0,0,0)
        call cwrite(lstr,0,i-1,0)
      endif

      errmx = 0
      do  j  = 1, ldim
      do  i  = 1, ldim

        err = cdabs(zhz(i,j))
        errmx = max(errmx,err)
        if (err > tol) then
          print 333, iq,ipq,ig,i,j,err
  333     format('trothf: mismatch iq=',2i5,' ig=',i2,' i,j=',2i3,'  zhz',f12.6)
C         call zprm('z+ H z - evl',2,zhz,ldim,ldim,ldim)
        endif
      enddo
      enddo

      if (iprint() >= 50) then
        print 334, errmx
  334   format(' max err=',1pe8.1)
      endif

      end
      subroutine pvtro5(kcplx,ldima,h,hd)
      implicit none
      integer kcplx,ldima
      double precision h(ldima,ldima),hd(ldima,2)

      integer ofh,i
      double precision wk(ldima,2)

      if (kcplx /= 0) call rx('not ready')

      ofh = ldima*ldima
      do  10  i = 1, ldima
        wk(i,1) = h(i,i)
        wk(i,2) = h(i+ofh,i)
   10 continue
      call dcopy(ldima*2,wk,1,hd,1)
      end

      subroutine pvtro6(i1,i2,i3,k1,k2,k3,ldima,nspc,ldimb,gfbz,gq)
      implicit none
      integer i1,i2,i3,k1,k2,k3,ldima,nspc,ldimb
      double complex gfbz(k1,k2,k3,ldima,nspc,ldimb,nspc)
      double complex gq(ldima,ldimb)
      double complex wk(ldima,ldimb)
      integer is,js,i,j

      is = 1
      js = 1

C     In line
      do  j  = 1, ldimb
      do  i  = 1, ldima
        wk(i,j) = gfbz(i1,i2,i3,i,is,j,js)
      enddo
      enddo

C     matrix copy does same job
      call zmcpy(' ',gfbz(i1,i2,i3,1,is,1,js),
     .  k1*k2*k3*ldima*nspc,k1*k2*k3,gq,ldima,1,ldima,ldimb)

C     check
      do  j  = 1, ldimb
      do  i  = 1, ldima
        if (wk(i,j) /= gq(i,j)) then
          print *, i,j
          stop 'oops'
        endif
      enddo
      enddo

      end

      subroutine pvtro7(nds,nttab,s)
      implicit none
      integer nds,nttab
      double complex s(nds,nds,nttab)
C     double precision s0(nds,nds,nttab)
      real(8), allocatable :: s0(:,:,:)
      double precision xx,xy

C     print *, 'allocate',nds*nds,nttab
      allocate(s0(nds,nds,nttab))
      s0 = dble(s)
      xx = max(abs(maxval(s0)),abs(minval(s0)))
      s0 = dimag(s)
      xy = max(abs(maxval(s0)),abs(minval(s0)))
      deallocate(s0)
      print *, 'maximum real part of s=',xx
      print *, 'maximum imaginary part of s=',xy

      end

      subroutine pvtro8(hreal,nhs,nttab,hrs,hrsc,hrss,hrssc)
C- find maximum diff(hrss,hrs)
      implicit none
      integer nhs,nttab,hreal
      double precision hrs(nhs,nhs,nttab),hrss(nhs,nhs,nttab)
      double complex hrsc(nhs,nhs,nttab),hrssc(nhs,nhs,nttab)
C     double precision h0(nhs,nhs,nttab)
      real(8), allocatable :: h0(:,:,:)
      double precision xx,xy
      integer i,j,k


C      call zprm('s(3)',2,hrsc(1,1,3),nhs,nhs,nhs)
C      call zprm('s(4)',2,hrsc(1,1,4),nhs,nhs,nhs)
C      call zprm('ss(3)',2,hrssc(1,1,3),nhs,nhs,nhs)

C     print *, 'allocate',nhs*nhs,nttab
      allocate(h0(nhs,nhs,nttab))
      if (hreal == 0) then
        h0 = dble(hrssc) - dble(hrsc)
        xx = max(abs(maxval(h0)),abs(minval(h0)))
        print *, 'maximum diff in real part of h=',xx

        do  k = 1, nttab
        do  j = 1, min(nhs,9)
        do  i = 1, min(nhs,9)
C          if (abs(h0(i,j,k)) > 1d-4) then
C            print *, i,j,k,h0(i,j,k)
C            stop
C          endif
        enddo
        enddo
        enddo

        h0 = dimag(hrssc) - dimag(hrsc)
        xy = max(abs(maxval(h0)),abs(minval(h0)))
        print *, 'maximum diff in imaginary part of h=',xy
      else
        h0 = dble(hrss) - dble(hrs)
        xx = max(abs(maxval(h0)),abs(minval(h0)))
        print *, 'maximum diff in h=',xx
      endif
      deallocate(h0)

      end

      subroutine pvtro9(lcplx,i,j,nds,s,sz)
C- print and record changes in s
      integer lcplx,i,j,nds
      double precision sold,s(nds,nds)
      double complex soldz,sz(nds,nds)
      save sold,soldz
      data sold/0d0/, soldz /(0d0,0d0)/

      return

      if (lcplx == 0) then
        if (sold /= s(i,j)) then
          print 333, ' pvtro9: s is now',
     .      s(i,j), ' diff=',s(i,j)-sold
  333     format(a,f11.6,a,f11.6)
          sold = s(i,j)
        endif
      else
        if (soldz /= sz(i,j)) then
          print 334, ' pvtro9: s is now',
     .      sz(i,j), ' diff=',sz(i,j)-soldz
  334     format(a,2f11.6,a,2f11.6)
          soldz = sz(i,j)
        endif
      endif
      end
      subroutine pvtroa(qp,plat,sq)
      double precision qp(3),plat(3,3)
      double complex sq,xxc
      integer iax(10,5)
      double precision tau(3,3),tdotk,cost,sint,twopi,taudotk,
     .  cosTau,sinTau,sumct,sumst,sumctau,sumstau
      integer j,k,isite
      data tau(1:3,1)/-0.288675d0, 0.500000d0,0.200324d0/
      data tau(1:3,2)/-0.288675d0,-0.500000d0,0.200324d0/
      data tau(1:3,3)/ 0.577350d0, 0.000000d0,0.200324d0/
      data iax(3:5,3) /-1, 0, 0/
      data iax(3:5,4) / 0, 0, 0/
      data iax(3:5,5) / 0, 1, 0/

      print *, 'qp=',qp
      twopi = 8*datan(1d0)
      sumct=0
      sumst=0
      sumctau=0
      sumstau=0

      do  isite = 3,5
        TdotK = 0
        do  30  j = 1, 3
        do  30  k = 1, 3
   30   TdotK = TdotK + twopi*qp(j)*plat(j,k)*iax(2+k,isite)
        taudotk = 0
        do  40  j = 1, 3
   40   TaudotK = TaudotK + twopi*qp(j)*tau(j,isite)

        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        sumct = sumct + cosT
        sumst = sumst + sinT

        cosTau = dcos(TaudotK)
        sinTau = dsin(TaudotK)
        sumctau = sumctau + cosTau
        sumstau = sumstau + sinTau
        print 333, isite,tdotk,cosT,sinT,taudotk,cosTau,sinTau
  333   format(i4,3f12.6,2x,3f12.6)
      enddo
      print 335, atan(sumst/sumct),sumct/3,sumst/3,
     .  atan(sumstau/sumctau),sumctau/3,sumstau/3
  335 format(' avg',3f12.6,2x,3f12.6)

      xxc = sq/6.6069767862013875d-05
      print 334, sq, xxc, abs(xxc),atan(dimag(xxc)/dble(xxc))
  334 format(2f12.6,2x,2f12.6,2x,2f12.6)

      end

      subroutine pvtrob(isp,ldim,s_lat,s_site,s_spec,s_ham,smpot,
     .  vconst,levec,lovl,q1,sq,evl)
C- Makes hamiltonian, overlap, or evec for a given q
      use structures
      implicit none
C ... Passed parameters
      logical levec,lovl
      integer isp,ldim
      double precision smpot(1),q1(3)
      double precision vconst,evl(ldim)
      double complex sq(ldim,ldim)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      real(8), allocatable :: wk2(:),ww(:)
      real(8), allocatable :: h(:)
      real(8), allocatable :: s(:)
C ... Local variables
      integer n1,n2,n3,ngabc(3),k1,k2,k3,nev,iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision plat(3,3),fac
      integer :: napw=0, igapw(1)

      if (napw /= 0) call rx('not ready for apw')

C     call upack('lat plat nabc',slat,plat,ngabc,0,0,0)
      plat = s_lat%plat
      ngabc = s_lat%nabc
      call fftz30(n1,n2,n3,k1,k2,k3)
C     call defrr(owk2,-ldim*ldim*2)
      allocate(wk2(ldim*ldim*2)); call dpzero(wk2,ldim*ldim*2)
C     call defrr(oh,-ldim*ldim*2)
      allocate(h(ldim*ldim*2)); call dpzero(h,ldim*ldim*2)
C     call defrr(os,-ldim*ldim*2)
      allocate(s(ldim*ldim*2)); call dpzero(s,ldim*ldim*2)
      fac = 1d0

      call fftz30(n1,n2,n3,k1,k2,k3)
      call hambl(0,s_site,s_spec,s_lat,s_ham,isp,q1,k1,k2,k3,
     .  smpot,vconst,0,0d0,ldim,napw,igapw,h,s,h)
      call dscal(ldim**2*2,fac,h,1)
C     call zprm('s',2,h,ldim,ldim,ldim)

      if (levec) then
        call info0(45,0,0,'pvtrob: poking evecs into s(k)')
      elseif (lovl) then
        call info0(45,0,0,'pvtrob: poking overlap into s(k)')
        call pvtroc(ldim,s)
        call zcopy(ldim**2,s,1,sq,1)
        call zcopy(ldim**2,s,1,h,1)
      else
        call info0(45,0,0,'pvtrob: poking hamiltonian into s(k)')
        call pvtroc(ldim,h)
        call zcopy(ldim**2,h,1,sq,1)
      endif

C     call zprm('s(k)',2,sq,ldim,ldim,ldim)
      call info2(45,0,0,' pvtrob: q1=%3;11,6D',q1,0)
      allocate(ww(11*ldim))
      call zhevx(ldim,ldim,h,s,1,.true.,ldim,9d9,nev,ww,.false.,evl,ldim,wk2)
      deallocate(ww)
      if (iprint() >= 45) call prtev(wk2,ldim,evl,ldim,9d9,nev)
      if (levec) then
        call zcopy(ldim**2,wk2,1,sq,1)
      endif
C     call zprm('s',2,sq,ldim,ldim,ldim)

      deallocate(wk2,h,s)

      end

      subroutine pvtroc(ndimh,s)
C- Add random number to Im s; symmetrize
      integer ndimh
      double complex s(ndimh,ndimh)
      real ran1
      double precision xx
      double complex z
      integer i1,i2
      double complex wk(ndimh,ndimh+1)

      return
      xx = maxval(abs(s))
      print *, 'pvtroc: largest value of s = ', xx

      return

C ... overwrite s by inv(s-z)
      print *, 'overwrite s by inv(s-z) ..'
C     Caution!!! call rsmsym with 10*1 ---
C     matrix is symmetric, not hermitian
      z = (0d0,.1d0)
      do  i1 = 1, ndimh
        s(i1,i1) = s(i1,i1) - z
      enddo
      call zqinv('n',s,ndimh,0,ndimh,wk,ndimh,i1)
      if (i1 /= 0) stop 'pvtroc: oops'

C ... add a random component
      return
      xx = ran1()/10
      print *,'add random component ...'

      do  i1 = 1, ndimh
        do  i2 = i1, ndimh
          xx = ran1()/10
          s(i2,i1) = s(i2,i1) + dcmplx(xx,ran1()/10)
        enddo
      enddo
      do  i1 = 1, ndimh
        do  i2 = i1, ndimh
          s(i1,i2) = dconjg(s(i2,i1))
        enddo
        s(i1,i1) = dble(s(i1,i1))
      enddo

      end

      subroutine pvtrod(sfz,k1,k2,k3,ldim)
      integer k1,k2,k3,ldim
      double complex sfz(k1,k2,k3,ldim,ldim)
      integer i1,i2

      i1 = 6
      i2 = 6
      call zprm3('gf(6,6)',0,sfz(1,1,1,i1,i2),k1,k2,k3)
      return

      call zprm3('gf(2,1)',0,sfz(1,1,1,2,1),k1,k2,k3)
      call zprm3('gf(1,2)',0,sfz(1,1,1,1,2),k1,k2,k3)
      end

      subroutine pvtrof(hrs,ndhrs,nsp,nk1,nk2,nk3,iax,is1,is2,hfz,ldim)
      implicit none
      integer ldim,nk1,nk2,nk3,ndhrs,nsp,is1,is2,niax
      parameter (niax=10)
      integer iax(niax,1)
      double complex hrs(ndhrs,ndhrs,nsp,is2)
      double complex hfz(nk1,nk2,nk3,ldim,ldim)
      integer i1,i2,ib1,ib2,is,ik1,ik2,ik3,j1,j2,isprt

C     pairs to match
      ib1 = 1
      ib2 = 2
C     orbitals within pairs
      i1 = 4
      i2 = 4
      isprt = 2

C      do  is = is1,is2
C        print 222, is, (iax(j1,is),j1=1,5)
C  222   format(i4,2x,5i3)
C      enddo
C      stop

      do  is = is1, is2
        if (iax(1,is) == ib1 .and. iax(2,is) == ib2) then
          ik1 = mod(iax(3,is)+6,6)+1
          ik2 = mod(iax(4,is)+6,6)+1
          ik3 = mod(iax(5,is)+6,6)+1
          print 336, is,ik1,ik2,ik3,hrs(i1,i2,1,is)
  336     format(' match at is=',i4,2x,'ik1,ik2,ik3=',3i3,2f12.6)
          do j1 = 1, ndhrs
          do j2 = 1, ndhrs
            hfz(ik1,ik2,ik3,j1,j2) = hrs(j1,j2,1,is)
          enddo
          enddo
        endif
      enddo

C     Prints out one orbital pair, all translation vectors
      call zprm3('hrs for 1 elt, all T',0,hfz(1,1,1,i1,i2),nk1,nk2,nk3)
      print *, 'matrix element for pair',isprt,
     .  ' ib1=',iax(1,isprt),' ib2=',iax(2,isprt)
      call zprm('hrs',2,hrs(1,1,1,isprt),ndhrs,ndhrs,ndhrs)

      end

      subroutine pvtroe(sfz,k1,k2,k3,ldim)
      integer k1,k2,k3,ldim
      double complex sfz(k1,k2,k3,ldim,ldim)
      integer i1,i2,i3
      double complex wk(ldim,ldim)

      i1 = 0
      i2 = 1
      i3 = 1
      print *, 'i1,i2,i3=?'
      read (*,*) i1,i2,i3
      if (i1 < 1 .or. i1 > k1) return
      wk(1:ldim,1:ldim) = sfz(i1,i2,i3,1:ldim,1:ldim)
      call zprm('s',2,wk,ldim,ldim,ldim)
      end


C      subroutine snot(h,ldh)
C      integer ldh
C      double complex h(ldh,ldh)
C
C      call zprm('h',2,h,ldh,ldh,ldh)
C
C      end
