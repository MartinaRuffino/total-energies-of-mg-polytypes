      subroutine ioxbs(ifi,nbas,nclass,alat,plat,range,rmax,dclabl,ipc,z,bas)
C- Write input file to xbs program
C ----------------------------------------------------------------
Ci Inputs
Ci   rmax : atom sizes
Ci   range: bonds drawn whose length < range(i)+range(j)
Cl Command-line arguments
Cl   -shift=x1,x2,x3 shifts all the coordinates by this amount,
Cl    in units of alat.
Cl   -dup=n1,n2,n3[,exprn] duplicates the unit cell n1 n2 n3 times
Cl    in the three lattice directions.  If 'exprn' is present, only
Cl    sites whose for which 'exprn' is satisfied are written out.
Cl    Variables that may be used in 'exprn' are:
Cl      alat,x1,x2,x3,p1,p2,p3,ic,ib,z,r
Cl      x1,x2,x3 are (dimensionless) site positions
Cl      p1,p2,p3 are the same positions, as projections onto plat
Cl      z is atomic number
Cl      r is radius
Cl   Example:
Cl     lmxbs "-dup=4,4,4,0<=p1&p1<.99&0<=p2&p2<.99&0<=p3&p3<.99"
Cl   picks out sites in first octant.
Cl     lmxbs "-dup=4,4,4,0<=x1&x1<1.01&0<=x2&x2<1.01&0<=x3&x3<1.01&z>0"
Cl   picks out sites in a cube of size alat, excluding empty spheres.
Cl
Cl   -sp:rule1[:rule2][...]
Cl    where 'rule' has the syntax  [exprn],[rmax],[red],[green],[blue]
Cl    Here exprn is a logical expression involving the class index ic
Cl    and atomic number z(ic).  Optional rmax, red, green, blue are
Cl    the ball size, and intensity of those colors for the ball.
Cl    Default values are taken for those numbers omitted.  For the
Cl    colors, they are random numbers; the default for rmax is whatever
Cl    is passed through to ioxbs.
Cl   Example:
Cl     lmxbs -sp:z==0,.5,.1,.1,.1:z==13&ic>3,,1,0,0:z==13,,0,1,0'
Cl    specifies that empty spheres are nearly black, with rmax=.5
Cl    and that Al atoms class index > 3 are red (rmax not spec'd)
Cl    and that all remaining Al atoms are green.
Cl   -scale=val
Cu Updates
Cu   24 Sep 18 Added --sites=list
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,nbas,nclass,ipc(nbas)
      double precision alat,plat(3,3),bas(3,nbas),z(nclass),
     .  range(nclass),rmax(nclass),dclabl(nclass)
C ... Dynamically allocated arrays
      integer, allocatable :: ilsts(:)
C ... Local parameters
      integer ib,ip,ic,jc,is,mode(3),it(3),ndup(3),parg,i,j,k,m,n,
     .  count,lgunit,ipr,ls,irul,nlsts
      double precision qlat(3,3),xx,r(3),r2(3),ddot,xxv(4),scale
C     double precision avw,avwsr,dum
      real ran1
      character*120 outs,cc*1,ctrm*3,clabl*8,clablj*8
      logical lcnst,lsw
      procedure(logical) :: cmdopt,a2bin
      procedure(integer) :: mkilsd,iinear

      call getpr(ipr)
      write (ifi,'(''* xbs file generate by program lmxbs'')')
      ls = len(outs)

      nlsts = 0
      if (cmdopt('-sites=',7,0,outs) .or. cmdopt('--sites=',8,0,outs)) then
        i = 8; if (outs(i:i) == '=') i=9
        nlsts = mkilsd(outs(i:),-1,it)
        if (nlsts <= 0) call rxs('ioxbs: bad or null site list : ',trim(outs))
        allocate(ilsts(nlsts))
        call mkilssr(11,outs(i:),nlsts,ilsts,[1,nbas])
      endif

      if (cmdopt('-shift=',7,0,outs)) then
        ip = 7
        call skipbl(outs,len(outs),ip)
        call dpzero(xxv,3)
        i = parg(' ',4,outs,ip,len(outs),', ',2,3,it,xxv)
        if (i < 0)
     .    call rxs2('ioxbs: failed to parse "',outs(1:ip+5),' ..."')
C       call prm('bas',' ',bas,3,3,nbas)
        call daxpy(nbas,1d0,xxv(1),0,bas(1,1),3)
        call daxpy(nbas,1d0,xxv(2),0,bas(2,1),3)
        call daxpy(nbas,1d0,xxv(3),0,bas(3,1),3)
C       call prm('bas',' ',bas,3,3,nbas)
      endif

C ... A random initializer for Knuth's random number generator
      call sectim(i,ic)
      call ran1in(ic/100)
C      if (rnge0 == 0) rnge0 = 2.7d0
C      avw = avwsr(plat,alat,dum,nbas)
C      dum = avw/alat*rnge0/2
C      do  10  ic = 1, nclass
C   10 if (range(ic) == 0d0) range(ic) = dum

C --- Write out site data ---
C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,xx)
      outs = ' '
      lcnst = .false.
      call iinit(ndup,4)
      if (cmdopt('-dup=',5,0,outs)) then
        ip = 5
        call skipbl(outs,len(outs),ip)
        i = parg(' ',2,outs,ip,len(outs),', ',2,3,it,ndup)
        mode(1) = 1
        mode(2) = 1
        mode(3) = 1
        lcnst = it(i) == 1
      endif

*     call shorps(nbas,plat,mode,bas,bas)

      if (ipr >= 40) then
        print 1
    1   format('   ic     ib',14x,'position',26x,'projection')
      endif
      count = 0
      call lodsyv('alat',1,alat,n)
      do  i = 0, ndup(1)
      do  j = 0, ndup(2)
      do  k = 0, ndup(3)
        do  ib = 1, nbas

          if (nlsts > 0) then
            ic = iinear(nlsts,ib,ilsts,1)
            if (ilsts(ic) /= ib) cycle
          endif

          ic = ipc(ib)
          call r8tos8(dclabl(ic),clabl)

C     ... Add bas + dup * lattice vector
          call dcopy(3,bas(1,ib),1,r,1)
          call daxpy(3,dble(i),plat(1,1),1,r,1)
          call daxpy(3,dble(j),plat(1,2),1,r,1)
          call daxpy(3,dble(k),plat(1,3),1,r,1)

          if (lcnst) then
C       ... projection of bas along plat(i)
            do  m = 1, 3
              r2(m) = ddot(3,r,1,qlat(1,m),1)
            enddo
            call lodsyv('ib',1,dble(ib),n)
            call lodsyv('ic',1,dble(ic),n)
            call lodsyv('z',1,z(ic),n)
            call lodsyv('x1',1,r(1),n)
            call lodsyv('x2',1,r(2),n)
            call lodsyv('x3',1,r(3),n)
            call lodsyv('p1',1,r2(1),n)
            call lodsyv('p2',1,r2(2),n)
            call lodsyv('p3',1,r2(3),n)
            call lodsyv('r',1,dsqrt(r(1)**2+r(2)**2+r(3)**2),n)
            n = ip
            call rxx(.not. a2bin(outs,lsw,0,0,' ',n,len(outs)),
     .               'iobxs: cannot parse exprn')
          else
            lsw = .true.
          endif
          if (lsw) then
            write (ifi,'(''atom  '',a,2x,3f12.6)')
     .        clabl, (alat*r(m),m=1,3)
            count = count+1
            if (ipr >= 40) print 363, ic,clabl,ib,r,r2
  363       format(i3,1x,a4,i4,2(f12.6,2f11.6))
          endif
        enddo
        enddo
        enddo
      enddo
      if (ipr >= 20) then
        n = (ndup(1)+1)*(ndup(2)+1)*(ndup(3)+1)*nbas
        call awrit3(' ioxbs:  wrote %i site positions%?#n# of %i#',
     .    ' ',80,lgunit(1),count,n-count,n)
      endif

C --- Write out spec data (random numbers for colors) ---
      write(ifi,'(a)')
      do  ic = 1, nclass
        call r8tos8(dclabl(ic),clabl)
        r(1) = ran1()
        r(2) = ran1()
        r(3) = ran1()
        call dscal(3,1/dsqrt(ddot(3,r,1,r,1)),r,1)
        xxv(1) = rmax(ic)
        call dcopy(3,r,1,xxv(2),1)
        if (cmdopt('-sp',3,0,outs)) then
          cc = outs(4:4)
          ip = 5
          irul = 0
C     ... Reentry for new sp
    4     continue
          irul = irul+1
          is = ip-1
          call chrpos(outs,',',ls,is)
          if (is >= ls) goto 39
C     j = fextg(ext)
C     call awrit0('%xls '//outs(ip:i)//'%a'//ext//
C     .      '%a | grep -s '//clabl,ss,len(ss),0)
C     call locase(ss)
C     call fsystm(ss,j)
            call lodsyv('ic',1,dble(ic),n)
            call lodsyv('z',1,z(ic),n)
            j = ip-1
            if (.not. a2bin(outs,lsw,0,0,',',j,is)) call
     .        rx('iobxs: failed to parse "'//outs(ip:is)//'"')
C       ... If it's this one
            if (lsw) then
              ctrm  = ', '//cc
              do  k = 0, 3
                is = is+1
                ip = is
                call chrps2(outs,ctrm,3,ls,is,n)
                if (ip == is) goto 35
                j = ip
                if (.not. a2bin(outs,xxv,4,k,ctrm(n:n),j,is)) call
     .            rx('iobxs: failed to parse "'//outs(ip:is)//'"')
   35           continue
C           ... The last argument
                if (n > 1) exit
              enddo
              print 3,ic,clabl,irul,xxv
    3         format(' spec',i3,':',a,2x,'rule',i2,'  r=',f8.5,
     .                '  color',3F7.3)
C       ... Otherwise, look for the next one
            else
              call chrps2(outs,' '//cc,2,ls,is,n)
              ip = is+2
              if (n == 2) goto 4
            endif
        endif
   39   continue

        write (ifi,'(''spec  '',a,f8.4,3f6.3)') clabl, xxv
      enddo

C --- Write out bond data ---
      write(ifi,'(a)')
      do  ic = 1, nclass
        call r8tos8(dclabl(ic),clabl)
        do  jc = ic, nclass
          call r8tos8(dclabl(jc),clablj)

          write (ifi,'(''bonds  '',a,2x,a,2f8.4,f8.3,2f6.2)')
     .      clabl, clablj, 0d0, (range(ic)+range(jc))/1, .1d0,
     .      0d0, 0d0
        enddo
      enddo

C --- Write out misc stuff ---
      call pioxbs('-scale=',20d0,scale)
      write(ifi,'(a)')
      call awrit1('scale  %d',' ',80,ifi,scale)
*     write (ifi,'(''scale  20'')')
      write (ifi,'(''dist   70'')')
      write (ifi,'(''inc     3'')')

      end
      subroutine pioxbs(strn,defval,parm)
      implicit none
      character*(*) strn
      double precision parm,defval
      character outs*120
      logical cmdopt,a2bin
      integer j

      j = len(strn)
      if (cmdopt(strn,j,0,outs)) then
        if (.not. a2bin(outs,parm,4,0,' ',j,len(outs)))
     .  call rx('ioxbs: failed to parse '//outs)
      else
        parm = defval
      endif
      end
