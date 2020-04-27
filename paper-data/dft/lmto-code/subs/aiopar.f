      logical function aiopar(alabl,lrel,pp,pprel,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl :class label
Ci   nl    :(global maximum l) + 1
Ci   lmax  :maximum l for a given site
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  pp    :scalar relativistic potential parameters (atomsr.f)
Cio        : pp(1) : enu
Cio        : pp(2) : calpha
Cio        : pp(3) : srdel = sqrt(delta) but with proper sign (that of phi-).
Cio        : pp(4) : palpha
Cio        : pp(5) : gamma, or Q in Varenna notes
Cio        : pp(6) : alpha, or qbar in Varenna notes
Cio        : Note: enu and c are relative to ves at MT boundary.  If ves shifts, enu and c should shift.
Cio  pprel :fpprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),ms1,ms2) fully relativistic potential parameters (atomsr.f)
Cio         :pprel(1,:,:,:) C
Cio         :pprel(2,:,:,:) gamma
Cio         :pprel(3,:,:,:) D+, where (D+) D = delta
Cio         :pprel(4,:,:,:) p
Cio         :pprel(5,:,:,:) enu
Cio  ves   :estat potential at rmax
Cio  aiopar:true unless read error or category not found.  In relativistic case, pprel must also be read
Cr Remarks
Cu Updates
Cu   25 Mar 15 Relativistic parameters include enu
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cr   11 Apr 94 Added convention: reading 1/p=0 => set p to 0
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      character alabl*8
      integer nl,lmax,nsp,ifi,lrel
      double precision ves,pp(6,nl,nsp),pprel(5,0:nl-1,2*nl,2,2)
C ... Local parameters
      logical lwarnp
      double precision xx,x2(0:8),mu
      equivalence (xx,x2)
      integer i,l,ip,ll,k,iprint,a2vec,ix2(9),nl2,ipr,nsp2,n,lmx
      integer mpipid,procid,imu,m,mumax
      logical scat,rdstrn,pars1v
      character s*100
C ... External calls
      external dpzero,getpr,info0,rx

      call getpr(ipr)
      aiopar = .false.
      lwarnp = .false. ! .true. => A warning message about p has already been printed
      procid = mpipid(1)

C --- File READ ---
      if (ifi > 0) then
      if (procid == 0) then

C   ... return unless file has PPAR: category
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (lmax /= nl2-1 .and. ipr >= 10)
     .    print *, 'aiopar (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp /= nsp2 .and. ipr >= 10)
     .    print *, 'aiopar (warning) mismatch in nsp, class '//alabl
        if (.not. pars1v(s,len(s),'ves=','=',4,ves))
     .    print *, 'aiopar (warning) failed to read ves, class '//alabl
        call dpzero(pp,6*nl*nsp)
        lmx = min(nl-1,nl2-1,lmax)
C   ... read ppar
        n = min(nsp,nsp2)
        read(ifi,'(a72)') s
        do  i = 1, n
        do  l = 0, nl2-1
          if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
          if (l <= lmx) then
            ip = 0
            k = a2vec(s,len(s),ip,4,' ',1,-2,-7,ix2,x2)
C    ...    Abort if failed to read 7 numbers
            if (k /= 7) call rx('AIOPAR, class'//alabl//'%a: failed to parse '//s)
            ll = x2(0)
C       ... no error if file has extra l's
            if (ll <= l) then
C         ... but pot pars must be available up to lmax
              if (ll /= l) call rx('AIOPAR: bad l quantum number')
C             Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
              x2(3) = dsign(1d0,x2(3))*dsqrt(dabs(x2(3)))
              if (x2(4) == 0) then
                if (.not. lwarnp) call info0(41,0,0,' AIOPAR: encountered 1/p=0 ... set p to 0')
                lwarnp = .true.
              else
                x2(4) = 1/x2(4)**2
              endif
              do  k = 1, 6
                pp(k,l+1,i) = x2(k)
                pp(k,l+1,nsp) = x2(k)
              enddo
            endif
          endif
        enddo
        enddo
      endif
C     Read relativistic ppar's
      if (procid == 0 .and. lrel == 2) then
C   ... return unless file has PPREL: category
        if (.not. scat(ifi,'PPREL:',':',.true.)) return
C        if (.not. scat(ifi,'PPREL:',':',.true.)) then
C          call dpzero(pprel,size(pprel))
C          call info0(20,0,0,' AIOPAR: missing PPREL ... copy enu from scalar rel')
C          do  l = 0, lmx
C            do  imu = 1, mumax
C              pprel(5,l,imu,1:2,1:2) = pp(1,l+1,1)/2 + pp(1,l+1,2)/2
C            enddo
C          enddo
C          return
C        endif
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (lmax /= nl2-1 .and. ipr >= 10)
     .    print *, 'aiopar (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp /= nsp2 .and. ipr >= 10)
     .    print *, 'aiopar (warning) mismatch in nsp, class '//alabl
        if (.not. pars1v(s,len(s),'ves=','=',4,ves))
     .    print *, 'aiopar (warning) failed to read ves, class '//alabl
        call dpzero(pprel,5*nl*2*nl*4)
        lmx = min(nl-1,nl2-1,lmax)
C   ... read ppar
        n = min(nsp,nsp2)
        read(ifi,'(a72)') s

        do  l = 0, lmx
          mumax = 2*(l+1)
          do  imu = 1, mumax
            do  i = 1, 2
              do  k = 1, 2
                if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
                if (l > lmx) goto 116
                ip = 0
                m = a2vec(s,len(s),ip,4,' ',1,-2,-9,ix2,x2)
                if (m == 8) then
                  call info0(2,0,0,' AIOPAR (warning) class '//alabl//'%a: substituting enu (rel)')
                  x2(8) = pp(1,l+1,1)/2 + pp(1,l+1,2)/2
C               Abort if failed to read 9 numbers
                elseif (m /= 9) then
                  call rx('AIOPAR (rel) class '//alabl//'%a: failed to parse '//s)
                endif
                mu = dble(imu-l) - 1.5d0
                if (dabs(mu-x2(1)) > 1d-12) goto 18
                ll = x2(0)
C           ... no error if file has extra l's
                if (ll > l) goto 116
C           ... but pot pars must be available up to lmax
                if (ll /= l) call rx('AIOPAR: bad l quantum number')
                do  m = 1, 5
                  pprel(m,l,imu,i,k) = x2(m+3)
                enddo
              enddo
            enddo
          enddo
        enddo
  116   continue

      endif

C     call mpibc1(pp,6*nl*nsp,4,.false.,'aiopar','pp')
C     call mpibc1(ves,1,4,.false.,'aiopar','ves')

C --- File WRITE ---
      else
        write(-ifi,21) alabl, lmax+1, nsp, ves
        do  i = 1, nsp
          do  l = 0, lmax
            xx = pp(3,l+1,i)
            if (dabs(pp(4,l+1,i)) < 1d-15) then
              write(-ifi,20) l, (pp(k,l+1,i), k=1,2), xx**2*dsign(1d0,xx), 0d0, (pp(k,l+1,i), k=5,6)
              lwarnp = .true.
            else
              write(-ifi,20) l, (pp(k,l+1,i), k=1,2), xx**2*dsign(1d0,xx),1/dsqrt(pp(4,l+1,i)), (pp(k,l+1,i), k=5,6)
            endif
          enddo
        enddo
        if (lwarnp) call info0(30,0,0,' AIOPAR: encountered p=0 ... wrote 1/p=0')
        aiopar = .true.

C       Write relativistic ppar's
        if (lrel == 2) then
          write(-ifi,22) alabl, lmax+1, nsp, ves
          do  l = 0, lmax
            mumax = 2*(l+1)
            do  imu = 1, mumax
              mu = dble(imu-l) - 1.5d0
              do  i = 1, 2
                do  k = 1, 2
                  write(-ifi,23) l,mu,i,k, (pprel(m,l,imu,i,k), m=1,5)
                enddo
              enddo
            enddo
          enddo
        endif

      endif

      aiopar = .true.
      return

   18 continue
      if (iprint() > 0) print *, 'aiopar: (input skipped) bad syntax, class '//alabl
      return

   20 format(i2,3f12.8,f12.7,3f12.8)
   23 format(i3,f5.1,2i4,5f13.8)
   21 format('PPAR:  ',a,'  nl=',i1,'  nsp=',i1,'  ves=',f12.8/' l',5x,
     .  'e_nu',10x,'C',8x,'+/-del',5x,'1/sqrt(p)',6x,'gam',9x,'alp')
   22 format('PPREL:  ',a,'  nl=',i1,'  nsp=',i1,'  ves=',f12.8/
     .  '  l  mu   ms1 ms2',5x,'C',11x,'gam',12x,'W',11x,'p',11x,'enu')

      end
