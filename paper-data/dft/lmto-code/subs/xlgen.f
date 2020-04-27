      subroutine xlgen(plat,rmax,rmax2,nvmax,opts,mode,nv,vecs)
C- Generate a list of lattice vectors, subject to constraints
C ----------------------------------------------------------------
Ci Inputs
Ci   plat  :dimensionless primitive lattice vectors
Ci   rmax  :largest radius for vector
Ci   rmax2 :if nonzero, largest radius for vector after
Ci         :multiples of plat are added
Ci   nvmax :maximum number of vectors allowed
Ci   opts  :1s digit:
Ci           1 add +/- any plat(j) to all other lattice vectors
Ci             found if plat(j) not in original list. (Ewald sums)
Ci          10s digit:
Ci           1 sort lattice vectors by increasing length
Ci           2 return nv only (or upper limit if 1s digit of opts is 1)
Ci           4 for padding, add a single vector and its reciprocal
Ci             instead of replicating entire original set.
Ci         100s digit:
Ci           1 return in vecs the multiples of plat(j) that
Ci             make up the lattice vector
Ci   mode  :vector of length 3 governing shifts along selected axes.
Ci         :0 suppresses shifts along plat(j)
Ci         :1 same as 0
Ci         :2 shifts to minimize length of pos
Co Outputs
Co   nv    :number of vectors found
Co   vecs  :list of vectors
Cu Updates
Cu  17 Jun 13 Replace f77 pointers with f90 ones
Cu  17 Mar 04 Bug fix for rmax2
Cu   2 Mar 04 New rmax2: truncate radius of lattice vectors to rmax2
Cu            when list has to be padded in order to include at
Cr            least one lattice vector.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nvmax,opts,mode(3)
      double precision plat(3,3),vecs(3,1),rmax,rmax2
C ... Local parameters
      double precision rsqr,v2,vj(3)
      integer i,imx(3),iprint,iv,ivck(3),j,jv,k,m,nglob,nv,stdo
C     integer oiwk,owk
      integer, allocatable :: iwk(:)
      real(8), allocatable :: wk(:)

C     call setpr(110)
      stdo = nglob('stdo')

C --- Setup ---
      call latlim(plat,rmax,imx(1),imx(2),imx(3))
      do  i = 1, 3
C   ... Switches flagging whether this plat in lattice vectors
        ivck(i) = 0
        if (mod(opts,10) == 1) ivck(i) = 1
        imx(i) = max(imx(i),1)
        if (mode(i) == 0) then
          imx(i) = 0
          ivck(i) = 0
        endif
      enddo
      rsqr = rmax*rmax
      nv = 0

C --- Loop over all triples, paring out those within rmax ---
      do  i = -imx(1), imx(1)
      do  j = -imx(2), imx(2)
      do  k = -imx(3), imx(3)
        v2 = 0
        do  m = 1, 3
          vj(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
          v2 = v2 + vj(m)**2
        enddo

C   --- A lattice vector found ---
        if (v2 > rsqr) cycle

C   ... Flag any plat in this vec as being present
        if (iabs(i) + iabs(j) + iabs(k) == 1) then
          if (i == 1) ivck(1) = 0
          if (j == 1) ivck(2) = 0
          if (k == 1) ivck(3) = 0
        endif

C   ... Increment nv and copy to vec(nv)
        nv = nv+1
        if (nv > nvmax .and. mod(opts/10,10) /= 2)
     .    call fexit(-1,111,' xlgen: too many vectors, n=%i',nv)
        if (mod(opts/10,10) == 2) then
        elseif (mod(opts/100,10) == 1) then
          vecs(1,nv) = i
          vecs(2,nv) = j
          vecs(3,nv) = k
        else
          vecs(1,nv) = vj(1)
          vecs(2,nv) = vj(2)
          vecs(3,nv) = vj(3)
        endif
      enddo
      enddo
      enddo

C --- Add plat if ivck ne 0 ---
      if (ivck(1)+ivck(2)+ivck(3) /= 0) then
        if (mod(opts/10,10) == 2) then
          nv = 3*nv
        elseif (mod(opts/10,10) == 4) then
          do  i = 1, 3
            if (ivck(i) == 1) then
              call dcopy(3,plat(1,i),1,vj,1)
              if (mod(opts/100,10) == 1) then
                vj(1) = 0
                vj(2) = 0
                vj(3) = 0
                vj(i) = ivck(i)
              endif
              vecs(1,nv+1) =  vj(1)
              vecs(2,nv+1) =  vj(2)
              vecs(3,nv+1) =  vj(3)
              vecs(1,nv+2) = -vj(1)
              vecs(2,nv+2) = -vj(2)
              vecs(3,nv+2) = -vj(3)
              nv = nv+2
            endif
          enddo
        else
          if (iprint() >= 20) write(stdo,1) ivck, rmax, rmax2
    1     format(/' xlgen: added missing plat: ivck=',3i2,'  rmax=',f8.3,'  rpad*rmax=',f8.3)
          if (3*nv > nvmax)
     .      call fexit(-1,111,' xlgen: too many vectors, n=%i',3*nv)
          if (ivck(1)+ivck(2)+ivck(3) /= 1)
     .      call rx('lgen: more than 1 missing plat ... try reducing EWALD_TOL')
          do  m = 1, 3
            v2 = ivck(1)*plat(m,1)+ivck(2)*plat(m,2)+ivck(3)*plat(m,3)
            if (mod(opts/100,10) == 1) v2 = ivck(m)
            call dcopy(nv,vecs(m,1),3,vecs(m,nv+1),3)
            call dcopy(nv,vecs(m,1),3,vecs(m,2*nv+1),3)
            call daxpy(nv, 1d0,v2,0,vecs(m,nv+1),3)
            call daxpy(nv,-1d0,v2,0,vecs(m,2*nv+1),3)
          enddo
          nv = 3*nv

C     ... Find and eliminate any replicas
          allocate(iwk(nv))
          call dvshel(3,nv,vecs,iwk,1)
*         call awrit2('%n:1i',' ',80,6,nv,iwk)
          k = 0
C     ... Mark any replica iv by iwk(i) -> -iv
          do  i = nv-1, 1, -1
            iv = iwk(i+1) + 1
            jv = iwk(i) + 1
            v2 = (vecs(1,iv)-vecs(1,jv))**2 +
     .           (vecs(2,iv)-vecs(2,jv))**2 +
     .           (vecs(3,iv)-vecs(3,jv))**2
            if (v2 < 1d-10) iwk(i+1) = -iv
          enddo

C     ... Flag vectors with radius > rmax2
          if (rmax2 > 0) then
            rsqr = rmax2*rmax2
            k = 0
            do  i = 0, nv-1
              if (iwk(i+1) >= 0) then
                iv = iwk(i+1) + 1
                v2 = vecs(1,iv)**2 + vecs(2,iv)**2 + vecs(3,iv)**2
                if (v2 > rsqr) then
                  if (iv < 0) call rx('bug in xlgen')
                  iwk(i+1) = -iv
                  k = k+1
                endif
              endif
            enddo
C           write(stdo,*) 'rpad reduced nv by',k
          endif
C         call prmx('unsorted vecs',vecs,3,3,nv)

C     ... Make a sorted list of replicas (any of iwk < 0)
          call ishell(nv,iwk)
C     ... For each replica, put lastmost vec into replica's place
          k = nv
          do  i = 0, nv-1
            iv = -iwk(i+1)
            if (iv <= 0) exit
            call dpcopy(vecs(1,k),vecs(1,iv),1,3,1d0)
            k = k-1
          enddo
          nv = k

          deallocate(iwk)
        endif
      endif
C     call prmx('after purging vecs',vecs,3,3,nv)

      if (mod(opts/10,10) == 2) return

C --- Sort vectors by increasing length ---
      if (mod(opts/10,10) == 1) then
        allocate(iwk(nv),wk(nv*3))
        call dvshel(3,nv,vecs,iwk,11)
*       call awrit2('%n:1i',' ',80,6,nv,iwk)
        call dvperm(3,nv,vecs,wk,iwk,.true.)
        deallocate(iwk,wk)
      endif

C --- Printout ---
      if (iprint() <= 70) return
      call awrit5(' xlgen: opts=%i  mode=%3:1i  rmax=%;4d  plat='//
     .  '%9:1;4d  nv=%i',' ',80,stdo,opts,mode,rmax,plat,nv)
      if (iprint() < 110) return
      write(stdo,2)
    2 format('  iv',6x,'px',8x,'py',8x,'pz',8x,'l')
      do  i = 1, nv
        v2 = vecs(1,i)**2 + vecs(2,i)**2 + vecs(3,i)**2
        write(stdo,3) i, (vecs(m,i),m=1,3), dsqrt(v2)
    3   format(i4,3F10.4,f10.3)
      enddo
      end
C      subroutine fmain
C
C      implicit none
C      integer wksize,nvmx,nv,opts,mode(3),ovecs
C      double precision plat(9),rmax
C      data plat /.5d0,.5d0,0d0, .5d0,-.5d0,0d0, 0d0,2d0,2d0/
C
C
C      nvmx = 500
C      rmax = 2
C
C      mode(1) = 0
C      mode(2) = 2
C      mode(3) = 2
C      opts = 0
C      call initqu(.true.)
C      call query('opts=?',2,opts)
C      call lgen(plat,rmax,nv,nvmx,w(ovecs))
C      write(stdo,*) 'old lgen found nv=', nv
C      write(stdo,*) 'call xlgen, opt=',opts
C      call xlgen(plat,rmax,nvmx,opts,mode,nv,w(ovecs))
C      end
