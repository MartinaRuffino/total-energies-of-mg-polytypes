      subroutine mk_hamindex(s_ham,s_lat,s_site,s_spec)
C- Writes orbital information for wave function rotation into file HAMindex
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham pwmode pwemax pwemin
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nsgrp plat qlat alat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:symgr ag pos
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:name
Cio    Passed to:  *
Cl Local variables
Cl   ipb   :index to true basis (excluding floating orbitals)
Cl         :given site index including those orbitals
cl   ipbi  :inverse mapping of ipb
Cr Remarks
Cr     iprmb is hard to understand (Mark's convension).  But, in anyway,
Cr     it is converted into clean indexing for Hamiltonian block.
Cr As you see in subroutine rotwv, the index for Hamiltonian reads as;
Cr     do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
Cr       ibas  = ibastab(iorb)
Cr        l    = ltab(iorb)
Cr        k    = ktab(iorb)        !kappa (different MTO index for each l)
Cr       init1 = offl(iorb)+1      !starting index for the block iorb
Cr       iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
Cr     enddo
Cu Updates
Cu   28 Aug 17 new ipbi that lists which atoms have augmentation spheres
Cu             ibpi is appended to HAMindex
Cu   Adapted from Takao's GW module m_hamindex
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
C     type(str_bz)::    s_bz
      type(str_ham)::   s_ham

      integer:: ibas,k,nl,stdo,iprint,ikt,nkt
      integer:: l,ndim,nglob,ldim,off,iorb,ib,is,norb
      real(8):: alat,dum
      integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:)
      integer,parameter :: n0nkap0=40
      integer:: pwmode,napwx,ig,nini
      integer,allocatable:: kv(:)
      real(8):: pwgmax, pwemax, pwgmin, pwemin, QpGcut_psi,qxx(3),
     .  qtarget(3),platt(3,3),q(3),qqx(3)
      integer:: ngp, ifiqg,iq,fopnx,nnn(3),ixx,nqbz
      integer:: ifi,ifisym,i,igg,iqq,iqi,irr,iqi_
      procedure(logical) :: latvec
      procedure(real(8)) :: dlength


c2012Sep02 kino, add for
      integer:: procid=0,ier=0
      integer,parameter::master=0

      integer:: nqtt,nqi,nqnum,ngpmx,napwmx
      integer nbas,ngrp,norbmto,kapmx,lmxbx,lmxax,lxxa,ndimham,imx,nat
      double precision plat(3,3),qlat(3,3),xx
      integer,allocatable:: ltab(:),ktab(:),offl(:),
     .  iclass(:),offlrev(:,:,:),ibastab(:),ipb(:),ipbi(:)
      integer,allocatable:: iqimap(:),iqmap(:),igmap(:),invgx(:),
     .  miat(:,:)
      integer,allocatable:: napwk(:),igv2(:,:,:),igv2rev(:,:,:,:)
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),
     .  shtvg(:,:),dlmm(:,:,:,:),qq(:,:)
      real(8),allocatable:: qtt(:,:),qtti(:,:)

C --- setup ---
C     if(.not.siginit) return
      call mpi_comm_rank( mpi_comm_world, procid, ier )
!       call mpi_comm_size( mpi_comm_world, numprocs, ier )

      nbas  = nglob('nbas')
      nl  = nglob('nl')
      stdo = nglob('stdo')
      ldim = s_ham%ldham(1)
      ngrp=s_lat%nsgrp
      plat=s_lat%plat
      qlat=s_lat%qlat
      alat=s_lat%alat
      allocate(symops(3,3,ngrp),ag(3,ngrp))
      call dcopy(9*ngrp,s_lat%symgr,1,symops,1)
      call dcopy(3*ngrp,s_lat%ag,1,ag,1)

      allocate(ltabx(n0nkap0,nbas),ktabx(n0nkap0,nbas),offlx(n0nkap0,nbas),ipb(nbas),ipbi(nbas))

C --- MTO part ---
C ... Obtain norbmto,ltabx,ktabx : orbital types for all atoms
      nat = 0     !number of sites with augmentation spheres
      norbmto = 0 !number of l's in mto part of hamiltonian
      kapmx = -1  !Largest kappa in system
      lmxbx = -1  !Largest basis lmax in system
      lmxax = -1  !Largest lmax in system
      ndimham = 0 !dimension of mto part of hamiltonian
      do  ib = 1, nbas
        is=s_site(ib)%spec
        ipb(ib) = 0
        if (s_spec(is)%lmxa >= 0) then
          nat = nat + 1
          ipb(ib) = nat
          ipbi(nat) = ib
        endif
        lmxax = max(lmxax,s_spec(is)%lmxa)
        call orbl(ib,0,ldim,s_ham%iprmb,norb,ltabx(1,ib),ktabx(1,ib),off,offlx(1,ib),ndim)
        do iorb = 1, norb
          norbmto = norbmto+1
          if(ltabx(iorb,ib)>lmxbx)  lmxbx = ltabx(iorb,ib)
          if(ktabx(iorb,ib)>kapmx)  kapmx = ktabx(iorb,ib)
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
        enddo
      enddo

C ... Make index table
C     allocate(ibasindex(ndimham))
      allocate(ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto))
      norbmto = 0 !number of l's in mto part of hamiltonian
      ndimham = 0 !dimension of mto part of hamiltonian
      if (iprint() > 40) then
        write(stdo,*) 'mk_hamindex:  Hamiltonian index '
        write(stdo,*) 'ib  l  k  start  end  spec'
      endif
      do  ib = 1, nbas
        is=s_site(ib)%spec
        call orbl(ib,0,ldim,s_ham%iprmb,norb,ltabx(1,ib),ktabx(1,ib),off,offlx(1,ib),ndim)
        do  iorb = 1, norb
          norbmto=norbmto+1
          ibastab(norbmto)= ib
          ltab(norbmto)   = ltabx(iorb,ib)
          ktab(norbmto)   = ktabx(iorb,ib)
          offl(norbmto)   = offlx(iorb,ib)
          nini = ndimham+ 1
          ndimham = ndimham+ 2*ltabx(iorb,ib)+1
C         ibasindex(nini:ndimham) = ib
          if (iprint() > 40) then
            write(stdo,"(3i3,2x,2i5,3x,a)")
     .        ib,ltab(norbmto),ktab(norbmto),
     .        offl(norbmto)+1,offl(norbmto)+2*ltab(norbmto)+1,
     .        trim(s_spec(is)%name)
          endif
        enddo
      enddo

C ... Reverse mapping of offset-index for hamiltonian
      allocate(offlrev(nbas,0:lmxbx,kapmx)); offlrev = 0
      do  iorb = 1, norbmto
        ibas = ibastab(iorb)
        l = ltab(iorb)
        k = ktab(iorb)
        offlrev(ibas,l,k) = offl(iorb)
      enddo

C --- Symmetry operations --- ---
      allocate(iclass(nbas),invgx(ngrp),
     .  miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      call sitepack(s_site,1,nbas,'class',1,iclass,xx)
C     Return space group information
      call mptauof(symops,ngrp,plat,nbas,s_lat%pos,iclass,miat,tiat,invgx,shtvg)

C --- Write SYMOPS ---- mar2012takao
      if (procid == master) then
      ifisym = fopnx('SYMOPS',2,0,-1) ! open as unknown
      write(ifisym,*) ngrp
      do  ig = 1, ngrp
        write(ifisym,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
        do  i = 1, 3
          write(ifisym,"(3e24.16)") symops(i,1:3,ig)
        enddo
      enddo
      call fclose(ifisym)
      endif

C ... Get rotation matrix dlmm
      lxxa = lmxax
      allocate( dlmm( -lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
      call dpzero(dlmm,size(dlmm))
      call rotdlmm(symops,ngrp, nl, dlmm)

!! no GW-related part
C      if(.not.llmfgw) then !.or.(llmfgw.and.jobgw==0)) then
C        print *,'gen_hamindex: no GW related part. not readin QGpsi'
C        return
C      endif
C

C --- Pw part. info for eigenfunctions are expanded as MTpart+PWpart.!feb2012takao ---
C      inquire(file='QGpsi',EXIST=qpgexist)  !feb2012takao
C      if(.not.qpgexist) then
C        return
C      endif

!! q on mesh and shortened q.
      ifiqg = fopnx('QGpsi',2,4,-1)
      read(ifiqg) nqnum, ngpmx, QpGcut_psi, nqbz, nqi,imx !,nqibz

c$$$cccccccccccccccccccccccccccccc
c$$$      nkt = 2*nqbz + 2*nkp
c$$$      if(allocated(qq)) deallocate(qq)
c$$$      allocate( qq(3,nkt) )
c$$$      print *,'gen_hamindex: nkt nqbz=',nkt,nqbz
c$$$      do  iq = 1, nqbz
c$$$        read(ifiqg)  qxx  ! q, and number of G vectors for
c$$$        read(ifiqg)
c$$$        qq(:,iq)=qxx
c$$$        call shorbz(qq(:,iq),qq(:,iq+nqbz+nkp),qlat,plat)
c$$$        write(stdo,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqbz+nkp)
c$$$      enddo
c$$$      print *
c$$$      call fclr(' ',ifiqg)
c$$$      do iq = nqbz+1, nqbz+nkp
c$$$        qq(:,iq) = qplist(:,iq-nqbz)
c$$$        call shorbz(qq(:,iq),qq(:,iq+nqbz+nkp),qlat,plat)
c$$$        write(stdo,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqbz+nkp)
c$$$      enddo
cccccccccccccccccccccccccccccc


!! === feb2012takao. ===
!! we have two set of data for original qxx in QGpsi and their shortened.
c      nqnum2 = 2*nqnum !+ 2*nkp !takao feb2012 test xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(allocated(qq)) deallocate(qq)
!! aug2012 we are removing shorbz
      nqtt = nqnum !nqnum*2 !doubled. second series nqnum+1:2*nqnum are for shortened q.
      nkt = nqtt
      allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
c      allocate( ngvecp(3,ngpmx,nqnum))
c      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
C      print *
C      print *,'gen_hamindex: Readin QGpsi. nqnum=',nqnum
      iqi = 0
      do  iq = 1, nqnum
        read(ifiqg)  qxx,ngp,irr  ! q, and number of G vectors for
        if(irr/=0) then
          iqi = iqi+1
          qtti(:,iqi) = qxx
          iqtt(iqi) = iq
        endif
        read(ifiqg)
c        read(ifiqg) ngvecp(1:3,1:ngp,iq), ngvecprev(-imx:imx,-imx:imx,-imx:imx,iq)
        qq(:,iq) = qxx
cc comment out aug2012takao
cc        call shorbz(qq(:,iq),qq(:,iq+nqnum),qlat,plat)
cc        write(stdo,"(' qq=',i5,3f10.5,' shortened qq=',3f10.5)") iq,qq(:,iq),qq(:,iq+nqnum)
C       write(stdo,"(' qq=',i5,3f10.5,'  irr ',i3)") iq,qq(:,iq),irr
      enddo
      call fclr(' ',ifiqg) !close file

!! ==== Generate info for rotwv and write ====
c     print *,' nqtt nqi=',nqtt,nqi
c     print *, ' symops=',symops
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      platt = transpose(plat) !this is inverse of qlat
      allocate(qtt(3,nqtt))
      qtt(:,1:nqtt) = qq(:,1:nqtt)
      do  i = 1, nqtt
        qtarget(:) = qtt(:,i)
c        print *,' i in nqtt=',i,qtarget
        do  iqi = 1, nqi
          q = qtti(:,iqi)
          iqq = iqtt(iqi)
          iqi_ = iqi
c          print *,' xxxxx qtti=',q
          do  ig = 1, ngrp
            if (latvec(1,1d-6,plat,qtarget-matmul(symops(:,:,ig),q))) then ! Equivalent
C             Warning if doesn't meet tighter tolerance
              if (.not. latvec(1,1d-7,plat,qtarget-matmul(symops(:,:,ig),q))) then
                xx = dlength(3,qtarget-matmul(symops(:,:,ig),q),1)
                call info2(2,1,0,' mk_hamindex (warning): small mismatch qtarget (%s,%3d), err=%;3g',
     .            qtarget,xx)
              endif
              igg = ig
C#ifdefC DEBUG
C              print *,' q qtarget     =',q,qtarget,ig
C              print *,' matmul q =',matmul(symops(:,:,ig),q)
C              print *
C#endif
              goto 2012
            endif
          enddo
        enddo
        call info2(1,1,0,' gen_hamindex: no q rotates to target %s,%3d.'//
     .    '%N Possibly inconsistent SYMGRP or QGpsi, positions require higher resolution.',
     .    qtarget,2)
        call rx('gen_hamindex: no qtarget found with given SYMOPS')
 2012   continue
        iqmap(i) = iqq
        iqimap(i) = iqi_
        igmap(i) = igg
      enddo

!! === rotation of APW. (not the IPW part for GW).===
c$$$ 1012 continue
      pwmode = s_ham%pwmode
      pwemax = s_ham%pwemax
      pwemin=s_ham%pwemin
      if(pwmode==0 .or. pwemax<1d-8) then
        if(allocated(napwk)) deallocate(napwk)
        allocate(napwk(nkt))
        napwk = 0
        napwmx = 0
        goto 900  ! put writehamindex inline
C        print *,'pwmode=0 writehamindex'
C        call writehamindex() !sep2012takao
C        return
      endif

!! for APW rotation.
C ... Get Igv2(3,iapw,ikt). pwmode>=10 only
C     print *,' goto APW part pwmode=',pwmode,pwemax !pwemin
      if(allocated(napwk)) deallocate(napwk,igv2,igv2rev)
      allocate( napwk(nkt))
!! takao is
      if(mod(pwmode,10)==0) then ! MTO basis only
        return
      endif
      pwgmax = dsqrt(pwemax)
      pwgmin = dsqrt(pwemin) !this will be removed.
      napwmx = 0
      call pshpr(1)
      do  ikt = 1, nkt
        qqx=qq(:,ikt)
        if(pwmode<10) qqx=0d0  !Is this call below to gvlst2 is consistent with other part?
        call gvlst2(alat,plat,qqx,0,0,0,pwgmin,pwgmax,0,0,0,napwx,dum,dum,dum,dum)
        napwk(ikt) = napwx
        if(napwmx<napwx) napwmx = napwx
      enddo
      call poppr
      allocate(igv2(3,napwmx,nkt),kv(3*napwmx) )

      do  ikt = 1, nkt
        qqx=qq(:,ikt)
        if(pwmode<10) qqx=0d0
        call pshpr(1)
        call gvlst2(alat,plat,qqx,0,0,0,pwgmin,pwgmax,0,2,napwmx,napwk(ikt),kv,dum,dum,igv2(1,1,ikt))
        call poppr
      enddo
      deallocate(kv)
C ... Find inverse mapping of igv2 --->igv2rev
      imx = -999
      do  ikt = 1, nkt
        ixx = maxval( abs(igv2(1:3,1:napwk(ikt),ikt)))
        if(ixx>imx) imx=ixx
      enddo
      allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
      igv2rev=999999
      do  ikt = 1, nkt
!       print *
        do  ig = 1, napwk( ikt )
          nnn  = igv2(1:3, ig, ikt)
          igv2rev( nnn(1), nnn(2),nnn(3), ikt) = ig
          !write(stdo,"(a,3f8.3,4i4,i6,i6)")'igv2rev: ',qq(:,ikt), nnn, ig, ikt
        enddo
      enddo
!     print *,' ---- nkt,napwmx norbmto= ',nkt,napwmx,norbmto
C     call writehamindex() !sep2012takao

C ... Replace call to writehamindex with inline code
C     Aug 2017 added nat, ipbi
  900 continue
      ifi = fopnx('HAMindex',2,4,-1)
      write(ifi) ngrp,nbas,kapmx,lmxbx,nqtt,nqi,nqnum,imx,ngpmx,norbmto,nat
      write(ifi) symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,
     .  iqimap
      write(ifi) lxxa
      write(ifi) dlmm
      write(ifi) ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
      write(ifi) qq !,ngvecp,ngvecprev
      write(ifi) plat,qlat,napwmx
      if(napwmx/=0) then        !for APW rotation used in rotwvigg
        write(ifi) igv2,napwk,igv2rev
      endif
      write(ifi) ipbi(1:nat)  ! Appended Aug 17
      call fclose(ifi)

      end

C     Taken from Kotani
      subroutine mptauof(symops,ng,plat,nbas,bas,iclass,miat,tiat,invg,delta)
C- Mapping of each atom by point group operations
Ci  Input
Ci     symops(1,ng),ng,plat,nbas,bas(3,nbas)
Ci     iclass(nbas); denote class for each atom
Co  Output
Co    space group operation does:
Co         bas' = matmul (symops(ig), bas) + delta(ig)     (1)
Co    symops is given; delta is calculated here
Co    delta : shifting vector for non-symmorphic group, called ag in lm package
Co    miat(ibas,ig); is atom that ibas is rotated into, but shifted by tiat.
Co    tiat(ibas,ig) shift to add to bas(:,miat(ibas,ig)) to coincide with bas', Eq (1).
Co    tiat must be some integer multiple of plat.  Thus:
Cr    bas(:,miat(ibas,ig)) + tiat(:,ibas,ig) = matmul(symop(ig),bas(:ibas)) + delta(:,ig)
Cu Updates
Cu   Coded by okuda 1994 March.
Cu   Simplified by kotani 1994 7/31.
C--------------------------------------------------------------------
      implicit none
      integer ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng),ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &                 tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &                 ep,dd1,dd2,dd3,t1,t2,t3
      integer  iprint
      external iprint
      integer ires(3, nbas, ng)
      integer:: ib1,ib2,PRT=100,rank=0
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
C     if(rank == 0) write(6,*)' MPTAUO2: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ) <= ep.and.
     &        abs( symops(2,ig)-symops(4,igd) ) <= ep.and.
     &        abs( symops(3,ig)-symops(7,igd) ) <= ep.and.
     &        abs( symops(4,ig)-symops(2,igd) ) <= ep.and.
     &        abs( symops(5,ig)-symops(5,igd) ) <= ep.and.
     &        abs( symops(6,ig)-symops(8,igd) ) <= ep.and.
     &        abs( symops(7,ig)-symops(3,igd) ) <= ep.and.
     &        abs( symops(8,ig)-symops(6,igd) ) <= ep.and.
     &        abs( symops(9,ig)-symops(9,igd) ) <= ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
        if(iprint() >= PRT) then
          if (rank == 0) then
          print *,' '
          print *,' '
          print *,' **** group ops no. ig (igd)= ', ig, invg(ig)
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
 1731     format (' ',3f9.4)
        endif
        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c       Look for translation vector ag (it is not supplied but calculated here and called tran
        do 20 ib1=1,nbas
        do 20 ib2=1,nbas
          tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))

c         (b1,b2,b3) = coordinates of ibas after rotation + shft by trial tran
          do 30 ibas=1,nbas
            b1=am(1,1)*bas(1,ibas)
     &        +am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas)
     &        +tran(1)
            b2=am(2,1)*bas(1,ibas)
     &        +am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas)
     &        +tran(2)
            b3=am(3,1)*bas(1,ibas)
     &        +am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas)
     &        +tran(3)
c
c           Find for atom mi for which (b1,b2,b3) - bas(:,mi) = lattice vector
            do 40 mi=1,nbas
              if( iclass(mi) /= iclass(ibas) ) go to 40

              do 50 i1=-3,3
              do 50 i2=-3,3
              do 50 i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1) <= ep.and.abs(t2) <= ep.and.
     &             abs(t3) <= ep) go to 60
   50         continue
   40       continue
C           Search failed. Try next (tr).
            goto 20

   60       continue
            miat(ibas,ig)  = mi
            tiat(1,ibas,ig)= dd1
            tiat(2,ibas,ig)= dd2
            tiat(3,ibas,ig)= dd3
            ires(1,ibas,ig)= i1
            ires(2,ibas,ig)= i2
            ires(3,ibas,ig)= i3
c
   30     continue
c When the do-30 loop has been completed, we get out of do-20 loop
          goto 21
   20   continue
        call rx( 'mptauo2: Can not find miat and tiat')
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

C       Translation found, debugging printout
        if (iprint() >= PRT) then
          if (rank == 0) then
          write(6,4658)tran
 4658     format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
            write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &      tiat(2,ibas,ig), tiat(3,ibas,ig),
     &      ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150     format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123       continue
          endif

          print 554, 'check formula'//
     .      ' bas(miat(ibas,ig)) = matmul(symop(ig),bas(ibas)) + delta(ig) - tiat(ibas,ig)'
  554     format(a)
          do  ibas=1,nbas
            tran = matmul(am,bas(:,ibas)) + delta(:,ig) - tiat(:,ibas,ig)
            tran = tran - bas(:,miat(ibas,ig))  ! This should be zero
            print 555, ibas, tran+bas(:,miat(ibas,ig)), tran
  555       format(i4,3f12.6,2x,3f12.6)
            if (sqrt(tran(1)**2+tran(2)**2+tran(3)**2) .gt. 1d-6) call rx('oops')
          enddo

        endif
   10 continue
      end

      subroutine rotdlmm(symops,ng,nl,dlmm)
C- Generate rotation matrix D^l_{m,m'} for L-representation, given point groups
C ----------------------------------------------------------------------
Ci Inputs
Ci   symops:point group operations
Ci   ng    :number of group operations
Ci   nl    :lmax+1
Co Outputs
Co   dlmm  :Rotation matrix D^l_{m,m'} that rotates Y_L(r) to Y_L(rotp r),
Co         :where _L are real harmonics
Co         :Dimensioned dlmm(-lmax:lmax,-lmax:lmax,0:lmax,ng), lmax = nl-1
Cu Updates
Cu   19 Aug 18 Generate D^l_{m,m'} by calling ylmrtg
Cu   04 May 16 Redesign of tolerance check
c-----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nl
      double precision symops(9,ng)
      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
C ... Local parameters
      integer ig,l,m1,m2
      real(8) rmat(nl*nl,nl*nl)

C#ifndef ORIG
      do  ig = 1, ng
        call ylmrtg(nl*nl,symops(1,ig),rmat)  ! lm code
        do  l = 0, nl-1
        do  m2 = -l, l
        do  m1 = -l, l
          dlmm(m2,m1,l,ig) = rmat(l*l+l+1+m2,l*l+l+1+m1)
          if (m2==8 .and. m1==7 .and. l==6) then
            print *, 'hi'
          endif
        enddo
        enddo
        enddo
      enddo
C#endif


C#ifdefC ORIG
C      double complex   dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
C      double precision det,igann,osq2
C      double complex   msc(0:1,2,2), mcs(0:1,2,2),dum(2)
CC     parameter(Img=(0d0,1d0))
C      integer i,j,md,m,mx,ikap,is
C      double precision am(3,3),fac1,fac2,amv(9),dlength
C      double precision cbeta,beta,sbeta,calpha,salpha,alpha,cgamma,sgamma,co2,so2,gamma,add
C      real(8),parameter:: tolg = 1d-7
C      complex(8),parameter:: Img=(0d0,1d0)
C      equivalence(amv,am)
C
C      do 10 ig =1,ng
C        do 20 i=1,3
C        do 20 j=1,3
C          am(i,j) = symops(i+3*(j-1),ig)
C   20   continue
Cc calculate determinant(signature)
C        det= am(1,1)*am(2,2)*am(3,3)
C     &      -am(1,1)*am(3,2)*am(2,3)
C     &      -am(2,1)*am(1,2)*am(3,3)
C     &      +am(2,1)*am(3,2)*am(1,3)
C     &      +am(3,1)*am(1,2)*am(2,3)
C     &      -am(3,1)*am(2,2)*am(1,3)
C        if(abs(abs(det)-1d0) >= 1d-10) then
C          print *,' rotdlmm: det/=1 ig and det=',ig,det
C          call rx(' ')
C        endif
Cc seek Euler angle   print *,' goto cbeta',ig,det
C        cbeta = am(3,3)/det
Cc added region correction so as to go beyond domain error for functions, dsqrt and acos.
C        if(abs(cbeta-1d0) <= 1d-6) cbeta= 1d0
C        if(abs(cbeta+1d0) <= 1d-6) cbeta=-1d0
C        beta = dacos(cbeta)
C        sbeta= sin(beta)
Cc beta= 0~pi
C        if(sbeta <= 1.0d-6) then
C          calpha= 1d0
C          salpha= 0d0
C          alpha = 0d0
C          cgamma= am(2,2)/det
C          sgamma= am(2,1)/det
C        else
C          salpha =  am(2,3)/sbeta/det
C          calpha =  am(1,3)/sbeta/det
C          sgamma =  am(3,2)/sbeta/det
C          cgamma = -am(3,1)/sbeta/det
C        endif
C        co2 = dcos(beta/2)
C        so2 = dsin(beta/2)
Cc         print *,' calpha=',calpha
C        if(abs(calpha-1.0d0) <= 1.0d-6) calpha= 1.0d0
C        if(abs(calpha+1.0d0) <= 1.0d-6) calpha=-1.0d0
C        if(abs(cgamma-1.0d0) <= 1.0d-6) cgamma= 1.0d0
C        if(abs(cgamma+1.0d0) <= 1.0d-6) cgamma=-1.0d0
C        alpha=dacos(calpha)
C        if(salpha < 0d0) alpha=-alpha
C        gamma=dacos(cgamma)
C        if(sgamma < 0d0) gamma=-gamma
Cc         print *,'alpha beta gamma det=',alpha,beta,gamma,det
C        do 30 l =  0, nl-1
C        do 30 md= -l, l
C        do 30 m = -l, l
Cc  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wiley and Sons.  p.52 (4.13)
C          fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
C          fac2 = 0d0
C          do 40 ikap=0,2*l
C            if(l-md-ikap >= 0 .and. l+m-ikap >= 0 .and. ikap+md-m >= 0) then
C              add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap)
C     &            *igann(ikap+md-m)*igann(ikap) )
C              if(2*l+m-md-2*ikap /= 0) add=add*co2**(2*l+m-md-2*ikap)
C              if(md-m+2*ikap /= 0)     add=add*(-so2)**(md-m+2*ikap)
C              fac2 = fac2+add
C            endif
C   40     continue
Cc l-th rep. is odd or even according to (det)**l
C          dlmmc(md,m,l,ig) = fac1*fac2*det**l*
C     &        cdexp( -Img*(alpha*md+gamma*m) )
C   30   continue
C
C        am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
C        am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
C        am(1,3)= sin(beta)*cos(alpha)
C        am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
C        am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
C        am(2,3)= sin(beta)*sin(alpha)
C        am(3,1)=-sin(beta)*cos(gamma)
C        am(3,2)= sin(beta)*sin(gamma)
C        am(3,3)= cos(beta)
C
C        if (dlength(9,amv(:)*det-symops(1:9,ig),1) > tolg)
C     .    call rx('Rotation by symgrp does not match rotation by Euler angle')
C
CC        if(abs(am(1,1)*det-symops(1,ig)) > 1.0d-8.or.
CC     &     abs(am(2,1)*det-symops(2,ig)) > 1.0d-8.or.
CC     &     abs(am(3,1)*det-symops(3,ig)) > 1.0d-8.or.
CC     &     abs(am(1,2)*det-symops(4,ig)) > 1.0d-8.or.
CC     &     abs(am(2,2)*det-symops(5,ig)) > 1.0d-8.or.
CC     &     abs(am(3,2)*det-symops(6,ig)) > 1.0d-8.or.
CC     &     abs(am(1,3)*det-symops(7,ig)) > 1.0d-8.or.
CC     &     abs(am(2,3)*det-symops(8,ig)) > 1.0d-8.or.
CC     &     abs(am(3,3)*det-symops(9,ig)) > 1.0d-8) then
CC          call rx('Rotation by symgrp does not match rotation by Euler angle')
CC        endif
C
C
Cc        if(iprint() >= 140) then
CC        if(.false.) then
CC          print *;print *;print *,' **** group ops no. ig=', ig
CC          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
CC          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
CC          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
CC          print *,' by Eular angle '
CC          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
CC          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
CC          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
CC        endif
C 1731   format (' ',3f9.4)
C
C   10 continue
Cc conversion to cubic rep. Belows are from csconvs
Cc  msc mcs conversion matrix generation 2->m 1->-m for m>0
C      osq2 = 1d0/sqrt(2d0)
C      do m = 0,1
C        Msc(m,1,1)= osq2*(-1)**m
C        Msc(m,1,2)=-osq2*Img*(-1)**m
C        Msc(m,2,1)= osq2
C        Msc(m,2,2)= osq2*Img
C
C        Mcs(m,1,1)= osq2*(-1)**m
C        Mcs(m,1,2)= osq2
C        Mcs(m,2,1)= osq2*Img*(-1)**m
C        Mcs(m,2,2)=-osq2*Img
C      enddo
Cc
C      do 23 is=1,ng
C        if(.false.) then
Cc        if(iprint() >= 150) then
C          print *; print *,' **** group ops no. ig=', is
C          write(6,1731) symops(1,is),symops(4,is),symops(7,is)
C          write(6,1731) symops(2,is),symops(5,is),symops(8,is)
C          write(6,1731) symops(3,is),symops(6,is),symops(9,is)
C        endif
Cc convert to cubic rep.
C      do 23   l =0,nl-1
C        do 33 m2=-l,l
C        do 33 m1= 1,l
C          dum(1)= dlmmc(m2, m1,l,is)
C          dum(2)= dlmmc(m2,-m1,l,is)
C          mx    = mod(m1,2)
C          dlmmc(m2,  m1,l,is)=
C     &                       dum(1)*msc(mx,1,1)
C     &                      +dum(2)*msc(mx,2,1)
C          dlmmc(m2, -m1,l,is)=
C     &                       dum(1)*msc(mx,1,2)
C     &                      +dum(2)*msc(mx,2,2)
C   33   continue
C        do 43 m2=  1,l
C        do 43 m1= -l,l
C          dum(1)=dlmmc( m2, m1,l,is)
C          dum(2)=dlmmc(-m2, m1,l,is)
C          mx=mod(m2,2)
C          dlmmc( m2, m1,l,is)=
C     &                       mcs(mx,1,1)*dum(1)
C     &                      +mcs(mx,1,2)*dum(2)
C          dlmmc(-m2, m1,l,is)=
C     &                       mcs(mx,2,1)*dum(1)
C     &                      +mcs(mx,2,2)*dum(2)
C   43   continue
C        do 53 m2=-l,l
C        do 53 m1=-l,l
C          dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
C          if( abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12 ) call rx( ''//
C     &     ' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12')
C   53   continue
Cccccccccccccccccccccc
C        if(.false.) then
Cc        if(.true.) then
Cc        if(iprint() >= 41) then
C          print *; print *,'  points ops  ig, l=', is,l,' cubic   '
C          do m2=-l,l
C            write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
Cc    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
C          enddo
C        endif
Ccccccccccccccccccccccc
C   23 continue
C
C#endif

      end
C#ifdefC ORIG
C      double precision function igann(i)
C      igann  = 1d0
C      do ix =1,i
C        igann=igann*ix
C      enddo
C      end
C#endif

C#ifdefC TEST
C      subroutine fmain
CC- Test original rotdlmm against ylmrtg
CC  Symops read from SYMOPS file
C      implicit none
C      integer ifi,ng,ig,i,l,m1,m2
C      integer, parameter :: nl=4, ngmx=48
C      procedure(integer) :: fopng
C      real(8) rmat(nl*nl,nl*nl)
C      real(8) :: symops(3,3,ngmx), dlmm(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ngmx)
C!     real(8),allocatable :: dlmm(:,:,:,:)
C
C
CC     read group ops
C      ifi = fopng('SYMOPS',-1,1); rewind ifi
C      read(ifi,*) ng
C      ng = min(ng,ngmx)
C      do  ig = 1, ng
C        read(ifi,*)
C        do  i = 1, 3
C          read(ifi,*) symops(i,1:3,ig)
C        enddo
C      enddo
C      close(ifi)
C
C      call rotdlmm(symops,ng,nl,dlmm) ! old GW code
C      do  ig = 1, ng
C        call ylmrtg(nl*nl,symops(1,1,ig),rmat)  ! lm code
C        print *
C        print *, 'group op',ig
C        do  i = 1, 3
C          call info2(2,0,0,'%3;12,7D',symops(i,:,ig),0)
C        enddo
C
CC       Compare
C        print *, 'rotation of real Ylm'
C        do  l = 0, nl-1
C        do  m2 = -l, l
C          write(*,"(25f12.7)") ( dlmm(m2,m1,l,ig), m1 = -l, l)
C        enddo
C        enddo
C        print *, 'subtract rmat'
C        do  l = 0, nl-1
C        do  m2 = -l, l
C        do  m1 = -l, l
C          dlmm(m2,m1,l,ig) = dlmm(m2,m1,l,ig) - rmat(l*l+l+1+m2,l*l+l+1+m1)
C        enddo
C        enddo
C        enddo
C        do  l = 0, nl-1
C        do  m2 = -l, l
C          write(*,"(25f12.7)") ( dlmm(m2,m1,l,ig), m1 = -l, l)
C        enddo
C        enddo
C        call info2(2,0,0,' max diff : %g',maxval(abs(dlmm(:,:,:,ig))),0)
C
C
C      enddo
C
C      end
C#endif
