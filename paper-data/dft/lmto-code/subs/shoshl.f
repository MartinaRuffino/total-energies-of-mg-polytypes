      subroutine shoshl(sopts,nbas,s_site,pos,alat,plat,mxnbr0,z,slabl,
     .  dclabl,ips,ipc,ves,eula,nclass)
C- Print nearest-neighbor shells
C ----------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  iopos
Ci Inputs
Ci   sopts: a set of modifiers, with the syntax
Ci          [:v][:e][:r=#][:sites:site-list][:pairs:pair-list] ..
Ci          [:tab[=#]][:disp=fnam][:nn][:fn=fnam]
Ci          :
Ci          :r=# sets range for shells
Ci          :v prints out electrostatic potential
Ci          :e prints out Euler angles
Ci          :r=# restricts neighbor table to range #
Ci          :mxcsiz=# Enlarges dimensioning of the neighbor table.
Ci          :sites:list print table only for sites in site-list
Ci          :pairs:pair-list print table only for pairs in pair-list
Ci          :tab prints out neighbor table
Ci          :tab=# prints out neighbor table in style #.  See ltab in psho1 for styles.
Ci          :  ...The following only apply to :tab
Ci          :  fn=fnam write table to file fnam
Ci          :  disp=fnam : read another site positions file; neighbor
Ci          :              table for both original and displaced
Ci          :              positions is written.
Ci          :  nn only print first entry for a given pair site
Ci          :     in neighbor table
Ci  OLD: doesn't work
Ci          :i[=style-#]:list  restricts neighbors in shell to list.
Ci                             This must be the last modifier.
Ci   nbas   :size of basis
Ci   pos    :basis vectors
Ci   alat   :length scale of lattice and basis vectors, a.u.
Ci   plat   :primitive lattice vectors, in units of alat
Ci   mxnbr0
Ci   z      :nuclear charge
Ci   slabl  :vector of species labels
Ci   dclabl :class name, packed as a real number
Ci   ips    :species table: site ib belongs to species ips(ib)
Ci   ipc    :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ves    :
Ci   eula   :Euler angles for noncollinear spins
Ci   nclass :number of inequivalent classes
Cb Bugs
Cb   mxnbr0 should be removed and replaced with mxcsiz
Cu Updates
Cu   14 Mar 19 Added mxcsiz
Cu   31 Aug 13 Obscure bug fix
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   19 Oct 12 Modified default maximum number of neighbors for table
Cu   19 Apr 03 Can write displaced neighbor table
Cu   12 Apr 03 Can write neigbhor table
Cr   24 Nov 97 changed modifier list
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nclass,mxnbr0,ips(nbas),ipc(nbas),nttab
      double precision alat,plat(3,3),pos(3,nbas),dclabl(nclass),ves(*),
     .  eula(*),z(*)
      character sopts*(*),slabl(nbas)*8
C ... For structures
!      include 'structures.h'
      type(str_site),intent(inout)::  s_site(*)
C ... Dynamically allocated local arrays
C#ifndef LINUXF
      integer, allocatable :: iax(:)
C#elseC
C      integer, pointer :: iax(:) => null()
C#endif
      integer, allocatable :: lst1(:),lst2(:),mark(:),ntab(:)
      real(8), allocatable :: dist(:),wk(:),rtab(:)
C ... Local parameters
      integer niax
      parameter (niax=10)
      logical lves,leula,a2bin
      double precision avwsr,avw,range,xx,stdo
      integer npr(2),mxnbr,ib,ic,j,j1,j2,lstyle,
     .  scrwid,ltab,nlstc,nlst1,nlst2,mxcsiz,ifi,ldisp,lnn,lstc(1),iv(10)
      character*8  dc*1
      parameter (scrwid=120)
      procedure(integer) :: lgunit,rdm,iprint,a2vec,cmdoptswx,wordsw,ival,fopn,iclbsj

C --- Parse modifiers ---
      stdo = lgunit(1)
      ifi = stdo
      lves  = .false.
      leula = .false.
      lstyle = 1
      range = 2
      nlstc = 0
      ltab = 0
      nlst1 = 0
      nlst2 = 0
      allocate(lst1(nbas),lst2(nbas),mark(nbas))
      allocate(dist(nbas*3))
      lnn = 0
      ldisp = 0
      call dcopy(nbas*3,pos,1,dist,1)
      mxcsiz = 0  ! Used only with :tab for now

      if (sopts /= ' ') then
C       ls = len(sopts)
        j1 = 1
        dc = sopts(j1:j1)
        j1 = j1+1

C   ... Gradually migrate options to new format. Start with mxcsiz.
        j2 = j1
        j = wordsw(sopts,dc,'mxcsiz','= ',j2) ! dimensioning for cluster table
        if (j /= 0) then
          if (a2vec(sopts,len_trim(sopts),j2,2,' '//dc,2,2,1,iv,mxcsiz) < 1) goto 999
        endif

C   ... Return here to resume parsing for arguments
   40   continue
        call nwordg(sopts,0,dc//' ',1,j1,j2)

C   ... Parse special arguments
        if (j2 >= j1) then
C         print *, sopts(j1:j2)
          if (sopts(j1:j2) == 'v')  then
            lves = .true.

          elseif (sopts(j1:j2) == 'e')  then
            leula = .true.

          elseif (sopts(j1:j1+1) == 'r=') then
            j = 0
            if (.not. a2bin(sopts(j1+2:),range,4,0,' ',j,j2-j1-2)) goto 999

          elseif (sopts(j1:j1+2) == 'fn=')  then
            if (j1+3 <= j2) ifi = fopn(sopts(j1+3:j2))

          elseif (sopts(j1:j1+3) == 'tab=')  then
            j = 0
            if (.not. a2bin(sopts(j1+4:),ltab,2,0,' ',j,j2-j1-4)) goto 999

          elseif (sopts(j1:j1+4) == 'disp=')  then
            call iopos(.false.,-1,sopts(j1+5:j2),nbas,dist,s_site)
            ldisp = 1

          elseif (sopts(j1:j1+1) == 'nn')  then
            lnn = 1

          elseif (sopts(j1:j1+2) == 'tab')  then
            ltab = 1

          elseif (sopts(j1:j1+6) == 'mxcsiz=') then ! migrated to above
C
          elseif (sopts(j1:j1+4) == 'sites')  then
            call baslst(0,11,sopts(j1+5:j2),j,ips,nbas,slabl,z,0,
     .        ' ',xx,nlst1,lst1)
            j2 = j1+5+j-2

          elseif (sopts(j1:j1+4) == 'pairs')  then
            call baslst(0,11,sopts(j1+5:j2),j,ips,nbas,slabl,z,0,' ',xx,nlst2,lst2)
            j2 = j1+5+j-2

          else
            call rxs('shoshl: failed to parse --shell switch: ',sopts(j1:j2))
            goto 999
          endif
          j1 = j2+2
          goto 40
        endif
      endif

C --- Print neighbor table for each site ---
      if (ltab /= 0) then
        call pshpr(iprint()-20)
        allocate(ntab(nbas+1))
        call pairs(nbas,nbas,1d0,plat,[range/2*(1+1d-6)],pos,
     .    [-1],3,-1,[0],nttab,ntab,iax,mxcsiz)
        call poppr
        allocate(rtab(3*nttab))
        call mkrtab(000,1d0,plat,pos,iax,nttab,pos,rtab)
C       if (ldisp == 1) ltab = 2
        call psho1(ldisp*10+ltab,lnn,iax,nbas,nttab,nlst1,lst1,
     .    nlst2,lst2,rtab,plat,pos,dist,mark,ipc,dclabl,ifi)
        deallocate(ntab,rtab,iax)
        goto 99
      endif

C --- Show shells for each class ---
      if (mxnbr0 == 0) then
        mxnbr = 2*(alat*range)**3
      else
        mxnbr = mxnbr0
      endif

      allocate(iax(niax*mxnbr),wk(mxnbr))
      avw = avwsr(plat,1d0,xx,nbas)
      avw = 1
c     call pshprt(50)
      do  ic = 1, nclass
        ib = iclbsj(ic,ipc,nbas,1)
        if (nlst1 > 0) then
          call hunti(lst1,nlst1,ib,0,j)
          if (ival(lst1,j+1) /= ib) cycle
        endif
        call nghbor(nbas,plat,pos,range,range,ib,mxnbr,npr,iax,wk)
        call xxsho(npr(1),nbas,plat,pos,iax,ipc,dclabl,nlstc,
     .    lstc,lves,ves,leula,eula,z)
      enddo
      deallocate(iax,wk)

   99 continue
      deallocate(lst1,lst2,mark,dist)
      return

  999 call rxs('shoshl: failed to parse ',sopts)
      end
      subroutine xxsho(npr,nbas,plat,pos,iax,ipc,dclabl,nlstc,lstc,lves,
     .  ves,leul,eula,z)
C- Kernel called by shoshl
C  nlstc,lstc:  a list of classes to include as pairs (nlstc>0)
      implicit none
      logical lves,leul
      integer npr,nbas,niax,ipc(*),nlstc,lstc(nlstc)
      parameter (niax=10)
      integer iax(niax,*)
      double precision plat(3,3),pos(3,1),dclabl(*),ves(*),eula(nbas,3),
     .  z(32)
      integer ih(2,120),scrwid
      parameter (scrwid=120)
      integer i,l,ishell,nshell,j,k,ii,kk,ic,jc,i1,lgunit,awrite,iclbsj,
     .  ib
      double precision dr(3),d,drr2,dshell,fuzz,z1(3),z2(3),alfa,beta,
     .  angle,pi,ddot
      character*8 clabl,outs1*25,outs2*(scrwid),outsv*(scrwid),
     .  outse*(scrwid)

      pi = 4*datan(1d0)
      fuzz = 1d-3
      dshell = 0
      nshell = 0
      ishell = 1
      if (leul) then
        alfa = eula(iax(1,1),1)
        beta = eula(iax(1,1),2)
        z1(1) = dcos(alfa)*dsin(beta)
        z1(2) = dsin(alfa)*dsin(beta)
        z1(3) = dcos(beta)
      endif
      ic = ipc(iax(1,1))
      call r8tos8(dclabl(ic),clabl)
      print 1,clabl,ic,nint(z(ic))
    1 format(/' Shell decomposition for class ',a,'  class',i4,'  z=',
     .        i2/' shell   d     nsh csiz  class ...')

      do  i = 1, npr
        d = dsqrt(drr2(plat,pos(1,iax(1,i)),pos(1,iax(2,i)),
     .    iax(3,i),iax(4,i),iax(5,i),dr))
C   ... new shell, or close of last shell
        if (dabs(d-dshell) > fuzz .or. i == npr) then
          i1 = i-1
          if (i == npr) i1 = i
          nshell = nshell+1
          write (outs1,2) nshell,dshell,i1+1-ishell,i1
    2     format(i4,f10.6,i4,i5,2x)
          call iinit(ih,2*(i-ishell))
C     ... ii is the number of different classes in this shell
          ii = 0
          do  j = ishell, i1
            ic = ipc(iax(2,j))
C       ... See whether already found one of these or if not in list
            kk = 0
            if (nlstc > 0) then
              kk = -1
              do  jc = 1, nlstc
                if (lstc(jc) > ic) exit
                if (lstc(jc) == ic) kk = 0
              enddo
            endif
            if (kk == 0) then
              do  k = 1, ii
                if (ih(2,k) /= ic) cycle
                kk = k
              enddo
            endif
C       ... We haven't --- increment ii and add this one
            if (kk == 0) then
              ii = ii+1
              kk = ii
              ih(2,kk) = ic
            endif
C       ... Increment number of occurrences of this species
            ih(1,kk) = ih(1,kk)+1
          enddo

C     ... Setup for printout
          outs2 = ' '
          outsv = ' '
          outse = ' '
          kk = 0
          do  k = 1, ii
            kk = kk+1
            call r8tos8(dclabl(ih(2,k)),clabl)
            l = awrite('%a  %np%i:'//clabl//
     .        '%a%?;n>1;(%i);%j;',outs2,len(outs2),0,
     .        (kk-1)*14,ih(2,k),ih(1,k),ih(1,k),0,0,0,0)
            if (lves) call awrit2('%np%d',outsv,len(outsv),0,
     .        (kk-1)*14,ves(ih(2,k)))
            ib = iclbsj(ih(2,k),ipc,-nbas,1)
            if (leul .and. ib > 0) then
              alfa = eula(ib,1)
              beta = eula(ib,2)
              z2(1) = dcos(alfa)*dsin(beta)
              z2(2) = dsin(alfa)*dsin(beta)
              z2(3) = dcos(beta)
              angle = dacos(max(-1d0,min(1d0,ddot(3,z1,1,z2,1))))
              if (angle > pi) angle = angle - 2*pi
              call awrit2('%np%d',outse,len(outse),0,(kk-1)*14,angle)
            endif
            if (l > scrwid-35) then
              call awrit0(outs1//outs2,' ',-scrwid,lgunit(1))
              if(lves) call awrit0('v%26f'//outsv,' ',-scrwid,lgunit(1))
              if(leul) call awrit0('e%26f'//outse,' ',-scrwid,lgunit(1))
              kk = 0
              outs2 = ' '
              outsv = ' '
            endif
          enddo
          if (outs2 /= ' ') then
            call awrit0(outs1//outs2,' ',-scrwid,lgunit(1))
            if (lves) call awrit0('v%26f'//outsv,' ',-scrwid,lgunit(1))
            if (leul) call awrit0('e%26f'//outse,' ',-scrwid,lgunit(1))
          endif
          outs1 = ' '

          ishell = i
          dshell = d
        endif
      enddo
      end
      subroutine shoang(sopts,nbas,pos,plat,mxnbr0,slabl,ips,zs)
C- Print bond angles
C ----------------------------------------------------------------
Ci Inputs
Ci   sopts :a set of modifiers, with the syntax
Ci         :  [:r=#][:spec=spec-list][:spec=spec-list]
Ci         :  :r=# sets range for shells
Ci         :  :sites=site-list collects angles only for sites within list
Ci         :  :spec=list       prints angles for bonds connecting to
Ci                             species in list
Ci   pos   :basis vectors, in units of alat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   mxnbr0
Ci   slabl :list of species labels
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   nbas  :size of basis
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   13 Sep 01 Added options sopts.  Altered argument list.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,mxnbr0
      character slabl(nbas)*8
      double precision pos(3,nbas),plat(3,3),zs(*)
      integer ips(nbas)
      character sopts*(*)
C ... Local parameters
      double precision avwsr,avw,range,xx
      integer npr(2),mxnbr,ib,nshell,nmxshl,niax,j,j1,j2,
     .  nspec,mxint,nsites,nbonds,
     .  ilst,parg,m,iv(10),ilast
      parameter (niax=10)
      character dc*1
      integer, allocatable :: slist(:),blist(:)
      integer, allocatable :: iax(:)
      real(8), allocatable :: wk(:)
      integer, allocatable :: num(:)
      real(8), allocatable :: ang(:)
      real(8), allocatable :: d(:)
C ... External calls
      external baslst,iinit,nghbor,nwordg,pvang1,pvang2,rxs

C ... Setup
      range = 2.5d0
      nspec = mxint(nbas,ips)
      nsites = 0
      nbonds = 0
      allocate(slist(nbas),blist(nbas))

C ... Switches
      dc = sopts(1:1)
      if (dc /= ' ') then
        j2 = 0
C   ... Return here to resume parsing for arguments
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) == dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 >= j1) then
          if (.false.) then

C         range
          elseif (sopts(j1:j1+1) == 'r=')  then
            m = 0
            j = parg('r=',4,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,range)
            if (j <= 0) goto 999

C         Site list
          elseif (sopts(j1:j1+4) == 'sites') then
            if (sopts(j1+5:j1+5) == '=') sopts(j1+5:j1+5) = dc
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zs,0,' ',
     .        xx,nsites,slist)
            j2 = j1+3+ilast

C         Bond list
          elseif (sopts(j1:j1+4) == 'bonds') then
            if (sopts(j1+5:j1+5) == '=') sopts(j1+5:j1+5) = dc
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zs,0,' ',
     .        xx,nbonds,blist)
            j2 = j1+3+ilast
          endif
          goto 10

        endif
      endif

      if (mxnbr0 == 0) then
        mxnbr = 2*range**3
      else
        mxnbr = mxnbr0
      endif
      allocate(iax(niax*mxnbr))
      allocate(wk(mxnbr))
      avw = avwsr(plat,1d0,xx,nbas)

C --- For every site in list, generate tables of bond angles ---
      ilst = 0
      do  ib = 1, nbas
        if (nsites > 0) then
          if (ilst+1 > nsites) cycle
          if (slist(ilst+1) /= ib) cycle
        endif
        ilst = ilst+1

C   ... Get neighbor lists
        call nghbor(nbas,plat,pos,range*avw,range*avw,ib,mxnbr,npr,iax,wk)

C   ... Get shell dimensions
        allocate(num(nspec*npr(1))); call iinit(num,nspec*npr(1))
        call pvang1(npr(1),nbas,plat,pos,iax,ips,num,nshell,nmxshl)
        deallocate(num)

C   ... Print bond angles
        j = nspec**2*nshell**2
        allocate(num(j)); call iinit(num,j)
        allocate(ang(j*nmxshl**2))
        allocate(d(nshell))
        call pvang2(npr(1),nbas,nspec,nshell,nmxshl,plat,pos,iax,
     .    ips,slabl,nbonds,blist,num,ang,d)
        deallocate(num,ang,d)
      enddo
      deallocate(iax,wk,slist,blist)

      return

  999 call rxs('shoang: failed to parse ',sopts)
      end

      subroutine pvang1(npr,nbas,plat,pos,iax,ips,num,nshell,nmxshl)
C- Help routine for shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of pairs in neighbor table
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information for one site.
Ci         :Table must be sorted by increasing distance from iax(1)
Ci   ips   :species table: site ib belongs to species ips(ib)
Co Outputs
Co   num   :(ishell,is) number of pairs in shell ishell of species is
Co   nshell:number of shells
Co   nmxshl:max value of num
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npr,nbas,nshell,nmxshl,niax,ips(*),num(npr,*)
      parameter (niax=10)
      integer iax(niax,*)
      double precision plat(3,3),pos(3,1)
C ... Local parameters
      integer i,is
      double precision d,wk(0:3),tol

      tol = 1d-6
C --- Get no. of shells and max no. of atoms in 1 shell and 1 class ---
      nshell = 0
      nmxshl = 0
      d = 0d0
      do  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,pos,iax(1,i),wk)
C       distance changed by more than tol ... new shell
        if (dabs(wk(0)-d) > tol) then
          nshell = nshell + 1
          d = wk(0)
        endif
        num(nshell,is) = num(nshell,is) + 1
        nmxshl = max0(nmxshl,num(nshell,is))
      enddo

      end

      subroutine pvang2(npr,nbas,nspec,nshell,nmxshl,plat,pos,iax,
     .  ips,slabl,nbonds,blist,num,ang,d)
C- Kernel called by shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of neighbors connecting site ib=iax(1,1)
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   nshell:number of shells
Ci   nmxshl:dimensions ang
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   slabl :struct containing global strings
Ci   num
Co Outputs
Co   ang   :table of angles
Co   d     :table of distances for each shell
Co   Angles and distances are printed out
Cu Updates
Cu   13 Sep 01
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer npr,nbas,nspec,nshell,nmxshl,niax,ips(*),
     .  num(nspec,nspec,nshell,nshell),nbonds,blist(nbonds)
      parameter (niax=10)
      integer iax(niax,*)
      character*8 slabl(*)
C ... Local parameters
      double precision plat(3,3),pos(3,1),d(nshell),
     .  ang(nmxshl**2,nspec,nspec,nshell,nshell)
      integer i,j,n,is,js,nsh1,nsh2,nmx2,k
      double precision rdtodg,d1,d2,dp,ddot,wk1(0:3),wk2(0:3)
C ... External calls
      external dlmn,hunti,rxx

      nmx2 = nmxshl**2
      rdtodg = 45d0 / datan(1.d0)

C --- Accumulate bond angles by shell and class ---
      nsh1 = 0
      d1 = 0d0
      do  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,pos,iax(1,i),wk1)
        if (dabs(wk1(0)-d1) > 1d-6) then
          nsh1 = nsh1 + 1
          d1 = wk1(0)
          d(nsh1) = d1
        endif
        nsh2 = nsh1
        d2 = d1
        if (nbonds > 0) then
          k = 0
          call hunti(blist,nbonds,iax(2,i),0,k)
          if (k >= nbonds) cycle
          if (blist(k+1) /= iax(2,i)) cycle
        endif
        do  j = i+1, npr
          js = ips(iax(2,j))
          call dlmn(nbas,plat,pos,iax(1,j),wk2)
          if (dabs(wk2(0)-d2) > 1d-6) then
            nsh2 = nsh2 + 1
            d2 = wk2(0)
          endif
          if (nbonds > 0) then
            k = 0
            call hunti(blist,nbonds,iax(2,j),0,k)
            if (k >= nbonds) cycle
            if (blist(k+1) /= iax(2,j)) cycle
          endif
          dp = ddot(3,wk1(1),1,wk2(1),1)
          if (dp > 1d0) dp =  1d0
          if (dp < -1d0) dp = -1d0
          if (nsh1 == nsh2 .and. js < is) then
            num(js,is,nsh1,nsh2) = num(js,is,nsh1,nsh2) + 1
            n = num(js,is,nsh1,nsh2)
            ang(n,js,is,nsh1,nsh2) = rdtodg*dacos(dp)
          else
            num(is,js,nsh1,nsh2) = num(is,js,nsh1,nsh2) + 1
            n = num(is,js,nsh1,nsh2)
            ang(n,is,js,nsh1,nsh2) = rdtodg*dacos(dp)
          endif
          call rxx(n > nmx2,'PVANG2: num gt nmx2')
        enddo
      enddo
      call rxx(nsh1 /= nshell,'PVANG2: nsh1 ne nshell')

C --- Printout ---
      call info2(1,1,0,
     .  ' Bond angles for site %i, species %i:'//trim(slabl(ips(iax(1,1))))//
     .  '%N shl1    d1    shl2    d2     cl1      cl2       angle(s) ...',
     .  iax(1,1), ips(iax(1,1)))

      do  nsh1 = 1, nshell
        do  nsh2 = nsh1, nshell
          do  is = 1, nspec
            do  js = 1, nspec
              n = num(is,js,nsh1,nsh2)
              if (n /= 0) then
                print 2, nsh1, d(nsh1), nsh2, d(nsh2), slabl(is),
     .            slabl(js), (ang(i,is,js,nsh1,nsh2), i = 1, n)
    2           format(2(1x,i3,1x,f9.6),1x,2(1x,a8),20(4(1x,f7.2):/47x))
              endif
            enddo
          enddo
        enddo
      enddo

      end
      subroutine psho1(ltab,lnn,iax,nbas,nttab,nlst1,lst1,nlst2,lst2,
     .  rtab,plat,pos,pos2,mark,ipc,dclabl,ifi)
C- Kernel called by supcel to displace pairs radially
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltab  :style which to print out neighbor table
Ci         :1s digit
Ci         :0 do nothing ;return.  Else print table as:
Ci         :1 (standard mode)
Ci         :    ib jb dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :2  (just the positions)
Ci         :   dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :3 (sparse matrix format)
Ci         :     1 jb dpos(1,jb)
Ci         :     2 jb dpos(2,jb)
Ci         :     3 jb dpos(3,jb)
Ci         :10s digit
Ci         :1 print out neighbor table for pos and
Ci         :  displaced pos2 as well
Ci   lnn   :If nonzero, only print first entry for a given pair site
Ci         :in neighbor table
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nbas  :size of basis
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   nlst1 :number of sites of "center" type
Ci   lst1  :list of sites of "center" type
Ci   nlst2 :number of sites of "neighbor" type
Ci   lst2  :list of sites of "neighbor" type
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   pos2  :displaced basis vectors (ltab >= 10)
Ci   mark  :work array of dimension nbas
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   dclabl:class name, packed as a real number
Ci   ifi   :file handle
Co Inputs/Outputs
Cu Updates
Cu   08 Aug 07 case ltab=12: allow for numerical imprecision in vector
Cu   19 Apr 03 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ltab,nlst1,lst1(nlst1),nlst2,lst2(nlst2),ifi,nbas,lnn
      integer niax,nttab,ipc(*),mark(nbas)
      double precision pos(3,*),pos2(3,*),rtab(3,nttab),dclabl(*),plat(9)
      parameter (niax=10)
      integer iax(niax,nttab)
C ... Dynamically allocated local arrays
      integer, allocatable :: iprm(:)
      real(8), allocatable :: d2(:,:)
C ... Local parameters
      integer iat,iap,low1,low2,i1,i2,ipr,stdo,i,oldi1,ic,jc
      character clabl*8, clabl2*8
C     character(len=120) :: buf(12)
      double precision d,ddot,dpos(3),qlat(9),tol
      procedure(integer) :: iprint,lgunit
      procedure(real(8)) :: dlength
C     logical latvec

      tol = 1d-6

      if (ltab == 0) return
      call mkqlat(plat,qlat,d)

      call info5(20,1,1,' ... shoshl: making neighbor list for'//
     .  ' %?#n#%-1j%i%j#%i# atom(s), style %i',nlst1,nbas,ltab,0,0)

      if (ltab == 3)
     .  call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)

      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
      call sanrg(.true.,mod(ltab,10),1,3,' shoshl:','tab')

      if (ltab == 3)
     .  call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)
      if (ltab == 13)
     .  call awrit1('%% rows 6 cols %i sparse',' ',80,ifi,nbas)

      low1 = 0
      low2 = 0
      oldi1 = 0

C     Reorder table within a site so pairs with approximately same d are ordered by iax(2)
      allocate(d2(6,nttab),iprm(nttab))
      do  iat = 1, nttab
        d2(1,iat) = iax(1,iat)
        d2(2,iat) = dlength(3,rtab(1,iat),1)
        d2(3,iat) = iax(2,iat)
        d2(4:6,iat) = rtab(1:3,iat)
      enddo
      call dvheap(6,nttab,d2,iprm,tol,101)

      do  iat = 1, nttab

        i1 = iax(1,iat)
        iap = iat        ! As given order
        iap = iprm(iat)  ! Sorted order order
        if (iax(1,iap) /= i1) call rx('bug in shoshl')

        i2 = iax(2,iap)

C   ... If site i1 isn't in the supplied list, skip this pair
        if (nlst1 /= 0) then
          call hunti(lst1,nlst1,i1,0,low1)
          if (low1 >= nlst1) cycle
          if (i1 /= lst1(low1+1)) cycle
        endif

C   ... New shell
        if (i1 /= oldi1) then
          call iinit(mark,nbas)
          oldi1 = i1
          ic = ipc(i1)
          call r8tos8(dclabl(ic),clabl)
          call awrit1('# neighbor list for site %i, class '//trim(clabl),
     .      ' ',80,ifi,i1)
          call info0(50,0,0,
     .      '#  ib  jb%10fCartesian coodinates%12fd/alat%33fMultiples of Plat')
     .
C         Don't print out on-site entry
C         goto 10
        endif

C   ... If site i2 isn't in the supplied lst2, skip this pair
        if (nlst2 /= 0) then
          call hunti(lst2,nlst2,i2,0,low2)
          if (i2 /= lst2(low2+1)) cycle
        endif

C   ... If i2 already marked, skip this pair
        if (mark(i2) /= 0) cycle

        if (mod(ltab,10) == 1) then
          jc = ipc(i2)
          call r8tos8(dclabl(jc),clabl2)
          d = dlength(3,rtab(1,iap),1)
          if (ltab < 10 .and. iprint() >= 50) then ! In multiples of plat
            call dgemm('T','N',3,1,3,1d0,qlat,3,rtab(1,iap),3,0d0,dpos,3)
            write (ifi,1) i1,i2,(rtab(i,iap),i=1,3),d,clabl,clabl2,dpos
          elseif (ltab < 10) then
            write (ifi,1) i1,i2,(rtab(i,iap),i=1,3),d,clabl,clabl2
    1       format(1x,2i4,3f12.7,2x,f12.7,2x,a,1x,a:2x,3f12.7)
          else
            do  i = 1, 3
              dpos(i) = pos2(i,i2)-pos(i,i2)
            enddo
            call shorbz(dpos,dpos,plat,qlat)
C           Only print out if displacement nonzero
            if (ddot(3,dpos,1,dpos,1) /= 0) then
                write (ifi,2) i1,i2,(rtab(i,iap),i=1,3),(dpos(i),i=1,3),d,clabl,clabl2
    2           format(1x,2i4,3f12.7,2x,3f12.7,2x,f12.7,2x,a,1x,a)
              mark(i2) = lnn
            endif
          endif
        elseif (ltab == 2) then
          write (ifi,3) (rtab(i,iap),i=1,3)
    3     format(1x,3f12.7)
          mark(i2) = lnn
        elseif (ltab == 12) then
          do  i = 1, 3
            dpos(i) = pos2(i,i2)-pos(i,i2)
          enddo
C         Only print out if displacement nonzero
          call shorbz(dpos,dpos,plat,qlat)
          if (ddot(3,dpos,1,dpos,1) > tol*tol) then
C         if (.not. latvec(1,tol,qlat,dpos)) then
            write (ifi,4) (rtab(i,iap),i=1,3),(dpos(i),i=1,3)
    4       format(1x,3f12.7,2x,3f12.7)
            mark(i2) = lnn
          endif
        elseif (ltab == 3) then
          do  i = 1, 3
            write (ifi,5) i,i2,rtab(i,iap)
    5       format(1x,i3,i5,3F12.7)
            mark(i2) = lnn
          enddo
        endif

      enddo
      deallocate(d2,iprm)

      end
