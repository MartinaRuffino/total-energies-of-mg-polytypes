      logical function aiomom(alabl,pl,ql,qlr,idmod,nl,lmax,nsp,z,rh,vrmax,bxc,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,lmax,nsp
Ci   nl:   dimensions pl,ql
Ci   ifi:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (atomsr.f)
Ci   ql:   moments q (atomsr.f) (ifi > 0)
Ci   qlr:  relativistic moments qr (atomsr.f) (ifi > 0)
Ci         qlr(1)=NULLI => array not allocated; do not attempt to read or write
Ci         qlr(1)=NULLI+1 => array allocated but not read.   Read but do not write
Ci     z:  nuclear charge. Used only on file read to get P for l between largest file lmax and lmax
Cio Inputs/Outputs
Cio   pl,ql:   moments q (atomsr.f)
Cio   aiomom:true unless read error or category not found
Cio   rh,vrmax density and Vxc at rmax
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   Bug in reading vrmax for nsp=2
Cu Updates
Cu   22 Feb 15 I/O of Fully relativistic moments, if qlr(1) > NULLI+1
Cu   12 Apr 11 I/O of direction of xc field defining moments
Cu   26 Apr 03 Added MPI calls
Cu   16 Apr 03 Read P,Q using a2vec
Cu   10 Apr 02 Set 'fill' pnu to free-electron value
Cu   16 May 01 On file read, do not rewind file
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,lmax,nsp,ifi
      character alabl*8
      integer idmod(0:nl-1)
      double precision pl(nl,nsp),ql(3,nl,nsp),qlr(4,nl,2*nl,2,2),rh,vrmax(2),bxc(3),z
C ... Local parameters
      integer,parameter :: NULLI=-99999
      logical scat,pars1v,parstr,rdstrn,lrel2
      integer i,l,ii,j,k,ipr,ix(6),nl0,nspf,a2vec,lmx,imu,mumax,stdo,nglob
      integer mpipid,procid
      double precision pi,dl,xx(6),scl,ddot,mu
      character s*100
CC#ifdefC AIX
CC      logical rdstrn
CC#endif
      call getpr(ipr)
      aiomom = .false.
      pi = 4*datan(1.d0)
      procid = mpipid(1)
      lrel2 = qlr(1,1,1,1,1) /= NULLI

C --- File READ ---
      if (ifi > 0) then
        if (procid /= 0) goto 99
        stdo = nglob('stdo')
C   ... Read scalar relativistic moments
        if (.not. scat(ifi,'MOMNTS:',':',.false.)) goto 101
C   ... read file nl and nsp
        backspace ifi
        read(ifi,'(a100)') s
        nl0 = nl; nspf = nsp; scl = 1d0
        if (.not. pars1v(s,len(s),'nl=','=',2,nl0)) call rx('AIOMOM: failed to find nl')
        if (lmax /= nl0-1 .and. ipr >= 10)
     .    call info2(10,0,0,' aiomom (warning) mismatch : lmax %i->%i, class '//trim(alabl),nl0-1,lmax)
        if (.not. pars1v(s,len(s),'nsp=','=',2,nspf)) call rx('AIOMOM: failed to find nsp')
        if (nsp /= nspf) then
          if (nsp > nspf) scl = dble(nspf)/dble(nsp)
          if (ipr >= 10) print *, 'aiomom (warning) mismatch in nsp, class '//alabl
        endif
        j = 0
        if (.not. parstr(s,'vrmax=',100,6,'=',j,i)) call rx('AIOMOM: failed to find rho,vrmax')
        if (a2vec(s,100,i,4,' ',1,j,3,ix,xx) /= 3) call rx('AIOMOM: failed to find rho,vrmax')
        rh = xx(1)
        vrmax(1) = xx(2)
        vrmax(2) = xx(3)
        j = 0
        if (parstr(s,'mag=',len(s),4,'=',j,i)) then
          if (a2vec(s,100,i,4,' ',1,j,3,ix,bxc) /= 3) call rx('AIOMOM: failed to read bxc')
        endif
        read(ifi,'(a100)') s
        lmx = min(nl0-1,lmax)
        call dpzero(pl,nl*nsp)
        call dpzero(ql,nl*nsp*3)
C       if (z > 0 .and. lmx < lmax) then
        if (lmx < lmax) then
          call defpq(1,z,lmax,1,pl,ql)
          if (nsp == 2) pl(:,2) = pl(:,1)
        endif
        call dpzero(ql,nl*nsp*3)
        do  i = 1, nspf
C        do  l = lmx+1, lmax
C          pl(l+1,i) = l+1 + .5d0 - datan(dble(l))/pi
C        enddo
        do  l = 0, nl0-1
          if (.not. rdstrn(ifi,s,len(s),.false.)) goto 12
          ii = 0
          if (a2vec(s,len(s),ii,4,' ',1,2,6,ix,xx) /= 6) goto 12
          ii = nint(xx(1))
          k = nint(xx(6))
C         read(ifi,100,err=12) ii, (xx(j), j=2,5), k
          if (ii /= l) call rx('AIOMOM: bad l quantum number')
          if (l <= lmx) then
            idmod(l) = k
            if (i <= nsp) then
              pl(l+1,i) = xx(2)
              pl(l+1,nsp) = xx(2)
            else
              pl(l+1,1) = (pl(l+1,1) + xx(2))/2
            endif
            do  ii = 1, 3
              if (i <= nsp) then
                ql(ii,l+1,nsp) = xx(2+ii)*scl
                ql(ii,l+1,i) = xx(2+ii)*scl
              else
                ql(ii,l+1,1) = xx(2+ii) + ql(ii,l+1,1)
              endif
            enddo
          endif
        enddo  ! loop over l
C ... Patch for bug in AIX err=
C#ifdefC AIX
C        if (i == 1 .and. nsp == 2) then
C          if (.not. rdstrn(ifi,s,len(s),.false.)) goto 19
C          if (s(1:5) /= '   0 ') goto 19
C          backspace ifi
C        endif
C#endif
        enddo ! loop over spin
        aiomom = .true.

C   ... Read fully relativistic moments
  101   continue
        if (.not. lrel2) goto 99
        if (.not. scat(ifi,'MOMREL:',':',.false.)) then
          do  l = 0, lmax
            do  imu = 1, 2*(l+1)
              qlr(4,l+1,imu,1,1) = NULLI ! Flag enu's not read
            enddo
          enddo
          goto 99
        endif
C       Future: check header line, as when reading MOMNTS; read via a2vec
        read(ifi,*) s
        do  l = 0, lmax
          mumax = 2*(l+1)
          do  imu = 1, mumax
            mu = dble(imu-l) - 1.5d0
            do  i = 1, 2
              do  k = 1, 2
                read(ifi,23) ix(1),xx(1),ix(3),ix(3),(qlr(ii,l+1,imu,i,k), ii=1,4)
              enddo
            enddo
          enddo
        enddo

        goto 99

C --- File WRITE ---
      elseif (procid == 0) then
        if (ddot(3,bxc,1,bxc,1) /= 0) then
          write(-ifi,333) alabl, rh, vrmax, lmax+1, nsp, bxc
        else
C         Use separate write statement to circumvent gnu compiler bug
          write(-ifi,332) alabl, rh, vrmax, lmax+1, nsp
        endif
  332   format('MOMNTS:  ',a4,'  rho,vrmax=',3f10.6,' nl=',i1,' nsp=',i1)
  333   format('MOMNTS:  ',a4,'  rho,vrmax=',3f10.6,' nl=',i1,' nsp=',i1:' mag=',3f9.6)
        write(-ifi,891)
  891   format('   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
        do  i = 1, nsp
          do  l = 0, lmax
            dl = dtan(pi*(.5d0 - pl(l+1,i)))
            if (dabs(dl) > 9999) dl = 0
C           write(-ifi,100) l,pl(l+1,i),(ql(ii,l+1,i),ii=1,3),idmod(l),dl
            call awrit5('%,4i%;13,7D%3;13,7D%,4i%;13,7D',' ',100,-ifi,
     .        l,pl(l+1,i),ql(1,l+1,i),idmod(l),dl)
          enddo
        enddo

C       Write the relativistic Q's
        lrel2 = lrel2 .and. qlr(1,1,1,1,1) /= NULLI+1
        if (.not. lrel2) return
        write(-ifi,22) alabl, rh, vrmax, lmax+1
        do  l = 0, lmax
          mumax = 2*(l+1)
          do  imu = 1, mumax
            mu = dble(imu-l) - 1.5d0
            do  i = 1, 2
              do  k = 1, 2
                write(-ifi,23) l,mu,i,k,(qlr(ii,l+1,imu,i,k), ii=1,4)
              enddo
            enddo
          enddo
        enddo
      endif
      return

   22 format('MOMREL:  ',a,'  rho,vrmax=',3f10.6,' nl=',i1/
     .  '  l  mu   ms1 ms2',5x,'q0',11x,'q1',11x,'q2',11x,'enu')
   23 format(i3,f5.1,2i4,4f13.8)


  100 format(i4,4f13.7,i4,f13.7)

C --- handle read exception ---
C#ifdefC AIX
C   19 l = 0
C      i = 2
C#endif
   12 continue
      if (l == 0 .and. i == 2) then
        if (ipr >= 20) print *, 'AIOMOM, ATOM=',alabl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  l = 0, lmax
          do  ii = 1, 3
            ql(ii,l+1,1) = ql(ii,l+1,2)
          enddo
        enddo
      else
        aiomom = .false.
      endif

   99 continue
C      call mpibc1(pl,nl*nsp,4,.false.,'aiomom','pl')
C      call mpibc1(ql,3*nl*nsp,4,.false.,'aiomom','ql')
C      call mpibc1(rh,1,4,.false.,'aiomom','ql')
C      call mpibc1(vrmax,2,4,.false.,'aiomom','ql')
      end
