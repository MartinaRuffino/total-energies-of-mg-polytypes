      subroutine mkqp(s_ctrl,s_bz,s_lat,gettet,lnoirr,lreduc,lgstar)
C- Set up k-points and related quantities for BZ integration
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: lpgf lmet
Co     Stored:    lmet
Co     Allocated: *
Cio    Elts Passed:lmet lsx lscr
Cio    Passed to: *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: nkabc lshft
Co     Stored:    star nkp nkabc ntet
Co     Allocated: qp wtkp idtet star ipq
Cio    Elts Passed:lopt lio qp wtkp idtet star
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: plat nsgrp npgrp
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:symgr
Cio    Passed to: *
Ci Inputs
Ci   gettet: T read or generate tetrahedra corners, if
Ci             tetrahedron integration set
Ci   lnoirr: T suppress generation of inequivalent tetrahedra
Ci   lreduc: 0 do not save array ipq
Ci         : 1 save array ipq
Ci         :-1 ignore symmetry operations, make qp for full BZ.
Ci   lgstar: nozero, generate igstar according to bzmesh, which see
Ci         : 0 igstar is not made
Ci         : 2 igstar contains inverse mapping of ipq
Ci         :-2 igstar contains group ops rotating irreducible
Ci         :   to full BZ.
Cio Inputs/Outputs
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Cio    Elts read: nkabc lshft lopt lio
Cio    Stored:    nkp nkabc qp wtkp star ntet idtet ipq
Cl Local variables
Cl   lipq  :T save array ipq
Cs Command-line switches
Cs  --findq=q1,q2,q3 | --findq~mq~q=q1,q2,q3
Cs         : When generating a k-mesh of points, 
Cs         : find the mesh point closest to q=(q1,q2,q3).
Cr Remarks
Cr   If s_ctrl%lwsig == LW8:
Cr   code will read qp from evecs file and rearrange the following arrays
Cr     s_bz%qp  s_bz%wtkp  s_bz%ipq  s_bz%idtet
Cr   in the qp order read from evecs.
Cu Updates
Cu  25 Aug 18 New option to locate closest approach to given qp (--findq)
Cu  18 Mar 18 Enable permutation of qp from contents of evec file
Cu  02 Jun 13 completed replacement of f77 pointers
Cu  06 Sep 11 Started migration to f90 structures
Cu  27 Jun 08 Adapt to new getqp.f
Cu  15 Sep 02 Can use sign of wgt to flag which irr points contain
Cu            equivalent points from time-reversal symmetry
Cu  21 Jul 02 Bug fix in second call to bzmesh
Cu   2 Feb 01 revised code to be consistent with comments (lreduc=0,1)
Cr   9 Oct 00 New lreduc, replacing lipq
Cr   6 Jan 98 (MvS) Split lnoirr into lnoirr+lipq options.
Cr  19 Nov 97 (WRL) added lpgf option, projecting qp to 2D
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical gettet
      integer lgstar,lreduc
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_bz)::    s_bz
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      integer, allocatable :: ipq(:),iwk(:)
      real(8), allocatable :: qp(:,:),wk(:)
C ... Local parameters
      integer, parameter :: LW8=8
      real(8),parameter:: tolq=1d-7
      character outs*72,dc*1
      logical ltet,lsx,lnoirr,lipq
      integer mxkp,nfilqp,nkp,nkxyz(3),nsgrp,npgrp,lshft(3),lpgf,
     .  lpbc,ntet,i,j,k,iq,stdo,ifi,nkabc(3),nsp,nspc
      double precision plat(3,3),qlat(3,3),xx,xv(10)
      procedure(integer) :: fopna,iprint,nglob,a2vec,wordsw
      procedure(logical) :: cmdopt

C ... Setup
      ntet = 0
      stdo = nglob('stdo')
      nkxyz = s_bz%nkabc
      lshft = s_bz%lshft
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      npgrp = s_lat%npgrp
      lpgf = s_ctrl%lpgf(1)
C     lpbc = 0 for kmesh in 3 dimensions, 1 kmesh in 2 dimensions
      lpbc = 0
      if (lpgf > 0) lpbc = 1
      ltet  = gettet .and. IAND(s_ctrl%lmet,2) == 2
C     call lsets('ctrl lmet',sctrl,ltet,2)
      if (ltet) then
        s_ctrl%lmet = IOR(s_ctrl%lmet,2)
      else
        s_ctrl%lmet = s_ctrl%lmet - IAND(s_ctrl%lmet,2)
      endif
      lsx  = IAND(s_ctrl%lsx,1) /= 0 .or. IAND(s_ctrl%lscr,1) /= 0
      lipq = lreduc == 1
      if (lreduc == -1) npgrp = 1

C ... q-points from BZMAP
      if (IAND(s_bz%lopt,2) /= 0) then
        call rx('recompile mkqp with BZMAP option')
C ... Read qp from disk
      elseif (IAND(s_bz%lio,1) /= 0) then
        call info0(30,0,0,' MKQP:   reading data from file QPTS ...')
        nfilqp = fopna('QPTS',-1,1)

        call getqp(0,nfilqp,nkp,nkxyz,lshft,ntet,xx,xx,xx)
        if (ltet) then
          if (ntet == 0) call rx('tetrahedron method specified but no tet weights given')
        else
          ntet = 0
        endif
        call ptr_bz(s_bz,8+1,'qp',3,nkp,xx)
        call ptr_bz(s_bz,8+1,'wtkp',nkp,0,xx)
        call ptr_bz(s_bz,1,'idtet',5,max(ntet,1),xx)
        call getqp(2,nfilqp,nkp,nkxyz,lshft,ntet,s_bz%qp,s_bz%wtkp,s_bz%idtet)

        call fclose(nfilqp)
        if (iprint() >= 20) call
     .    awrit1(' MKQP:   read %i qp from disc',' ',80,stdo,nkp)
        call rxx(ltet,'tet. integration with non-standard k-mesh')
        if (lgstar /= 0) then
          call rx('mkqp: lgstar not allowed with user supplied k-mesh')
        endif
C ... Make the qp list from bzmesh
      else
        ! Get number of points nkp, and printout what will be made
C       i = 3; if (lgstar /= 0) i = 4
        call pshpr(10)
        call bzmshp('  MKQP',-1,nkxyz,lshft,plat,s_lat%symgr,npgrp,
     .    nsgrp,.false.,lgstar,lpbc,qlat,nkp,xx,xx,xx,xx)
        call poppr

C       Allocate arrays
        mxkp = nkxyz(1)*nkxyz(2)*nkxyz(3)
        if (lgstar == 0) then
          call ptr_bz(s_bz,8+1,'star',1,0,xx)
        else
          call ptr_bz(s_bz,8+1,'star',mxkp+1,0,xx)
        endif
        call ptr_bz(s_bz,8+1,'wtkp',nkp,0,xx)
        call ptr_bz(s_bz,1,'qp',3,nkp,xx)
C       Generate ipq in a local array
C       Whether it is copied to s_bz%ipq depends on lreduc.
        allocate(ipq(mxkp))

C        call dinv33(plat,1,qlat,vol)
CC   ... Restrict BZ to two dimensions
C        if (lpbc == 1) then
C          outs = ' ' // prgnam
C          if (nkxyz(3) > 1 .and. iprint() >= 10) then
C            write(stdo,*) ' '
C            call awrit2('%a (warning): nk3=%i, shft3=%i; reset to 1,0',
C     .      outs,80,-stdo,nkxyz(3),lshft)
C          endif
C          lshft(3)=0
C          nkxyz(3) = 1
C          call projql(qlat)
C        endif
C
C        do  i = 1, 3
C          llshft(i) = lshft(i) /= 0
C        enddo
C        s_bz%star(1) = lgstar

C   ... Sanity check
        if (lsx .and. lshft(1)+lshft(2)+lshft(3) > 0) call
     .    rx('MKQP:  shifted BZ mesh not allowed with SX')

C   ... Process switch --findq
        k = 0
        if (cmdopt('--findq',7,0,outs)) then
          dc = outs(8:8)
          if (dc == '=') then
            k = 8
          else
            k = wordsw(outs,dc,'q=','',i) + 2
          endif
          if (a2vec(outs,len(outs),k,4,', ',2,-2,3,nkabc,xv) /= 3)
     .      call rx('failed to parse '//trim(outs))
          s_bz%qp(1:3,1) = xv(1:3)
          if (wordsw(outs,dc,'mq','',i) > 0) then
            call dgemm('N','T',1,3,3,1d0,xv,1,s_lat%qlat,3,0d0,s_bz%qp,1)
          endif
          k = 10
        endif

C   ... 2nd call stores the arrays
        i = 3+k; if (lgstar /= 0) i = 4+k
        call bzmshp('  MKQP',i,nkxyz,lshft,plat,s_lat%symgr,npgrp,nsgrp,
     .    .false.,lgstar,lpbc,qlat,nkp,s_bz%qp,s_bz%wtkp,ipq,s_bz%star)

C   ... Copy to s_bz%ipq
        if (lreduc /= 0) then
          call ptr_bz(s_bz,4+1,'ipq',mxkp,0,ipq)
        endif

C   ... Generate inequivalent tetrahedra
        if (ltet .and. .not. lnoirr) then
          call ptr_bz(s_bz,1,'idtet',5,6*mxkp,xx)
C         Require work array dimensioned 6*mxkp
          allocate(iwk(6*mxkp))
          call icopy(mxkp,ipq,1,iwk,1)
          call tetirr(qlat,nkxyz(1),nkxyz(2),nkxyz(3),iwk,ntet,s_bz%idtet)
          call ptr_bz(s_bz,2,'idtet',5,ntet,xx)
          deallocate(iwk)
        else
          call ptr_bz(s_bz,1,'idtet',1,1,xx)
        endif

        deallocate(ipq)

      endif

C --- Pack new info into structures ---
C     call prmx('qp',s_bz%qp,3,3,nkp)
      s_bz%nkp = nkp
      s_bz%nkabc = nkxyz
      s_bz%ntet = ntet

C --- Permute qp list to align with evec file ---
      if (s_ctrl%lwsig == LW8) then
        call info0(30,0,0,' MKQP:   permuting qp as read from evec file ...')
        ifi = fopna('evec',-1,4); rewind ifi
        call iosigh(2,5,nsp,nspc,i,j,nkabc(1),nkabc(2),nkabc(3),nkp,iq,lshft(1),lshft(2),lshft(3),ifi,xv)
        if (nkabc(1) /= nkxyz(1) .or. nkabc(2) /= nkxyz(2) .or. nkabc(3) /= nkxyz(3) .or. s_bz%nkp /= nkp)
     .    call rx2('evec file mesh (%i irr qp) is incompatible with given mesh (%i irr qp)',nkp,s_bz%nkp)
C     .   call fexit3(-1,111,' Exit -1 MKQP: '//
C     .    'evec file mesh (%s,%3i) is incompatible with given mesh (%s,%3i) with (%2i) irr points',
C     .    nkabc,nkxyz,[s_bz%nkp,nkp])
C    .    call rx2('evec file mesh (%s,%3i) is incompatible with given mesh (%s,%3i)',nkabc,nkxyz)
        allocate(qp(3,nkp),iwk(nkp),wk(nkp),ipq(nkp))
        do  iq = 1, nkp
        do  i = 1, nsp
          read(ifi) qp(1:3,iq)
          read(ifi)
          read(ifi)
        enddo
        enddo
        call alignqp(nkp,qp,nkp,s_bz%qp,plat,tolq,iwk,i) ! iwk = permutation list
        call rxx(i /= 0,'evec file mesh is incompatible with given mesh')
        call dcopy(nkp,s_bz%wtkp,1,wk,1)
C       Permute qp, wtkp
        do  iq = 1, nkp
          s_bz%qp(:,iq) = qp(:,iq)    ! s_bz%qp(iwk(iq)) -> s_bz%qp(iq)
          s_bz%wtkp(iq) = wk(iwk(iq)) ! s_bz%wk(iwk(iq)) -> s_bz%(iq)
          ipq(iwk(iq)) = iq           ! inverse permutation:  ipq(iwk(iq))) = iq
        enddo
C       Replace s_bz%ipq(i) -> s_bz%ipq(iwk(i))
        if (lreduc /= 0) then
          do  i = 1, mxkp
            j = s_bz%ipq(i); iq = ipq(j) ! j = original iq and iq = permuted iq
            s_bz%ipq(i) = iq
          enddo
        endif
        if (ntet > 0) then
          do  i = 1, ntet
            do  k = 2, 5
              j = s_bz%idtet(k,i); iq = ipq(j) ! k = original iq and iq = permuted iq
              s_bz%idtet(k,i) = iq
            enddo
          enddo
        endif

      endif

C --- Write q-points to disc ---
      if (IAND(s_bz%lio,2) /= 0) then
        nfilqp = fopna('QPTS',-1,0)
        call getqp(2,-nfilqp,nkp,nkxyz,lshft,ntet,s_bz%qp,s_bz%wtkp,s_bz%idtet)
        call fclose(nfilqp)
        call rx0('done writing qp to disk')
      endif

      end
