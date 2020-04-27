      subroutine angtab(sopts,nbas,bas,alat,rmax,qss,qlat,dlat,nkd,
     .  ipc,slabl,ips,zs,neul,eula,nl,qnu)
C- Print angles between spins on different sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :a set of modifiers, with the syntax
Ci         :  [~r=rmax[,rmin]][~sites=site-list]
Ci         :  :r=# sets range for shells
Ci         :  :sites=site-list collects angles only for sites within list
Ci         :  :spec=list       prints angles for bonds connecting to
Ci                             species in list
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   rmax  :augmentation radius, in a.u.,
Ci   qss   :Parameters defining spin-spiral
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   dlat  :direct lattice vectors
Ci   nkd   :number of dlat
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   slabl :list of species labels
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   eula  :Euler angles for noncollinear spins
Co Outputs
Cl Local variables
Cl         :
C  BUG: SS addition to Euler angle assumes along z, and no beta!
Cr Remarks
Cr
Cu Updates
Cu   15 Apr 11 New 'sign' switch
Cu   10 Feb 10 updated, new switches
Cu   26 Jan 03 got working again
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer nbas,nkd,neul,ipc(*),ips(nbas),nl
      double precision dlat(3,nkd),bas(3,nbas),rmax(*),alat,slabl(*),
     .  eula(nbas,neul,3),qss(4),qlat(3,3),qnu(3,nl,2,*),zs(*)
C ... Local parameters
      logical lsign
      double precision range(2),xx
      integer n1,i1,i2,j1,j2,irep,jrep,nrep,ifi,i,j,k,i1mach,ic,jc,ii,jj
      integer parg,m,iv(10),ilast,nsites,slist(nbas)
      integer NULLI
      parameter (NULLI=-99999)
      double precision z1(3),z2(3),angle(nbas),ddot,pi,
     .  dd0(3),dd(3),xd(3),dd1,sumrs,ovlprs,alfa,beta,dsum
      character outs*80
C      logical cmdopt,a2bin
      character dc*1
C     double precision rotm(3,3)

      if (neul /= 1) call rx('angtab not ready for l-dependent eula')

C ... Switches
      lsign = .false.
      range(1) = 0
      range(2) = 0
      nsites = nbas
      do  i = 1, nbas
        slist(i) = i
      enddo
      if (qss(1) == NULLI) then
        call dpzero(qss,4)
      endif
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

C         fold signs of moments into angle
          elseif (sopts(j1:j1+3) == 'sign')  then
            lsign = .true.

C         r<range and r>range(2)
          elseif (sopts(j1:j1+1) == 'r=')  then
            m = 0
            j = parg('r=',4,sopts(j1:),m,len(sopts(j1:)),
     .        ', '//dc,2,2,iv,range)
C           if (j == 2) call dswap(1,range(1),1,range(2),1)
            if (j <= 0) goto 999

C         Site list
          elseif (sopts(j1:j1+4) == 'sites') then
            if (sopts(j1+5:j1+5) == '=') sopts(j1+5:j1+5) = dc
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zs,0,' ',
     .        xx,nsites,slist)
            j2 = j1+3+ilast

C         Bond list
C          elseif (sopts(j1:j1+4) == 'bonds') then
C            if (sopts(j1+5:j1+5) == '=') sopts(j1+5:j1+5) = dc
C            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zs,0,' ',
C     .        xx,nbonds,blist)
C            j2 = j1+3+ilast
          else
            goto 999

          endif
          goto 10

        endif
      endif

C --- Print out Euler angles ---
      if (ddot(3,qss,1,qss,1) > 0) call
     .  awrit2('%x  Qss =%3:1;6d angle '//'%;6d',outs,80,-i1mach(2),qss,qss(4))
      print 344, (i, (eula(i,1,j), j=1,3), i=1,nbas)
  344 format('  ib      alpha       beta       gamma'/(i4,3f12.6))
      print *, ' '

C --- Print out Euler angles between each pair ---
      pi = 4*datan(1d0)
      ifi = i1mach(2)
      n1 = 8   ! Number of columns
      nrep = (nsites-1)/n1 + 1
      do  irep = 1, nrep
      i1 = (irep-1)*n1+1
      i2 = min0(irep*n1,nsites)
      do  jrep = 1, nrep
        j1 = (jrep-1)*n1+1
        j2 = min0(jrep*n1,nsites)
        write(ifi,121) (slist(j),j=j1,j2)
  121   format(3x,8i8)
        write(ifi,122) ('--------',j=j1,j2)
  122   format(7x,'---',8a8)
        do  ii = i1, i2
        i = slist(ii)
        ic = ipc(i)

C        call eua2rm(eula(i,1,1),eula(i,1,2),eula(i,1,3),rotm)
CC       print 335, ((rotm(k,jj),jj=1,3),k=1,3)
CC  335  format((3f15.9))
C        z1(1) = rotm(3,1)
C        z1(2) = rotm(3,2)
C        z1(3) = rotm(3,3)

        alfa = eula(i,1,1)
        beta = eula(i,1,2)
        z1(1) = dcos(alfa)*dsin(beta)
        z1(2) = dsin(alfa)*dsin(beta)
        z1(3) = dcos(beta)

        do  jj = j1, j2
          j = slist(jj)
          jc = ipc(j)

C     ... Add any rotation from qss(z).
          do  m = 1, 3
            dd0(m) = bas(m,j) - bas(m,i)
          enddo
          call shortn(dd0,dd,dlat,nkd)
C     ... This gives multiples of plat shortn took away
          do  k = 1, 3
            xd(k) =  (dd(1)-dd0(1))*qlat(1,k)
     .             + (dd(2)-dd0(2))*qlat(2,k)
     .             + (dd(3)-dd0(3))*qlat(3,k)
          enddo
C     ... Assume here qss is along z, qlat(3) also
          alfa = eula(j,1,1)
          if (qss(3) /= 0) then
            alfa = eula(j,1,1) - 2*pi*xd(3)*qss(3)/qlat(3,3)
          endif
          beta = eula(j,1,2)
          z2(1) = dcos(alfa)*dsin(beta)
          z2(2) = dsin(alfa)*dsin(beta)
          z2(3) = dcos(beta)
Cc          call eua2rm(alfa,eula(j,1,2),eula(j,1,3),rotm)
C          z2(1) = rotm(3,1)
C          z2(2) = rotm(3,2)
C          z2(3) = rotm(3,3)
          angle(jj) = dacos(max(-1d0,min(1d0,ddot(3,z1,1,z2,1))))
          if (lsign) then
C            print *, i,j,
C     .        (dsum(nl,qnu(1,1,1,ic),3) - dsum(nl,qnu(1,1,2,ic),3))*
C     .        (dsum(nl,qnu(1,1,1,jc),3) - dsum(nl,qnu(1,1,2,jc),3))
            if ((dsum(nl,qnu(1,1,1,ic),3) - dsum(nl,qnu(1,1,2,ic),3))*
     .          (dsum(nl,qnu(1,1,1,jc),3) - dsum(nl,qnu(1,1,2,jc),3))
     . < 0) angle(jj) = angle(jj) + pi
          endif
          if (angle(jj) > pi) angle(jj) = angle(jj) - 2*pi
          if (angle(jj) > pi) angle(jj) = angle(jj) - 2*pi
        enddo
        write(outs,130) i,(180/pi*angle(jj),jj=j1,j2)
  130   format(1x,i4,8f8.2)
        do  jj = j1, j2
          j = slist(jj)
          jc = ipc(j)
          do  m = 1, 3
            dd0(m) = bas(m,j) - bas(m,i)
          enddo
          call shortn(dd0,dd,dlat,nkd)
          do  m = 1, 3
            dd(m) = dd(m)*alat
          enddo
          dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          sumrs = rmax(ic) + rmax(jc)
          ovlprs = sumrs - dd1
          if (i == j) then
            outs(6+8*(jj-j1):5+8+8*(jj-j1)) = '   --'
          elseif (range(1) /= 0) then
            if (dd1 > range(1) .or. dd1 < range(2))
     .        outs(6+8*(jj-j1):5+8+8*(jj-j1)) = ' '
          elseif (100*ovlprs/dd1 < -1) then
            outs(6+8*(jj-j1):5+8+8*(jj-j1)) = ' '
          endif
        enddo
        call awrit0('%a',outs,-80,-i1mach(2))
      enddo
      enddo
      enddo

      return
  999 call rxs('angtab: failed to parse ',sopts)
      end
