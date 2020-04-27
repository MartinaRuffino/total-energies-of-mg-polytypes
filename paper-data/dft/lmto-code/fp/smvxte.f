      subroutine smvxte(job,s_lat,s_pot,n1,n2,n3,eps0,smpot,smrho)
C- Setup and perform analysis of response to an external potnetial
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  ng alat gv
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:kv
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  smvextc smcv0
Co     Stored:     smcv0
Co     Allocated:  *
Cio    Elts passed:smcv0 smvextc
Cio    Passed to:  *
Ci Inputs
Ci   job   :10s digit
Ci         : 0 smpot,smrho are in r space
Ci         : 1 smpot,smrho are in k space
Ci         : 1s digit
Ci         : 0 do nothing
Ci         :-1 count number of nonzero elements in s_pot%smvextc; return it in n1
Ci         :-2 Add estimate for screening charge to smrho
Ci         :-4 For each nonzero element in s_pot%smvextc, copy smpot(G) to s_pot%smcv0
Ci         :-6 Combination of -2 and -4
Ci         :-7 Combination of -2 and -9
Ci         :-8 For each nonzero element in s_pot%smvextc, copy smpot(G) to s_pot%smcv0
Ci         :   and write to file smpot0
Ci         :-9 For each nonzero element in s_pot%smvextc, read s_pot%smcv0 from file smpot0
Ci         :   if it exists; otherwise revert to job -4
Ci         : [Note: jobs -2, -6, -7 all estimate screening charge]
Ci         : 1 For each nonzero element, estimate eps^-1(G) as eps^-1 = dvext/dvtot
Ci         : 2 For each nonzero element, estimate response function as R = drho/dvext
Ci   n1,n2,n3 : (not applicable to job -1) number of divisions in smpot, smrho
Ci   eps0  : (jobs -2, -6, -7 only) dielectric constant use to estimate screening charge
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   smrho :smooth density on uniform mesh
Co Outputs
Co   n1    : (for job -1 only) number of nonzero Fourier components in external potential
Co         : (for job -8 only) file smpot0 is created
Co         : Other jobs: printout only
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   28 Jan 17 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer job,n1,n2,n3
      double precision eps0
      complex(8), target :: smpot(n1,n2,n3),smrho(n1,n2,n3)
C ... For structures
!       include 'structures.h'
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
C     complex(8), allocatable :: smfun(:,:,:)
      complex(8), allocatable :: cv0(:,:),cn(:)
C ... Local parameters
      logical lerr,lrspace
      integer k1,k2,k3,ig,ng,icv0,ncv0,job0,ifi
      double precision q2,tpiba,pi,pi8,fac,zx(2,5)
      complex(8) :: epsi,dphitot,drho
      procedure(integer) :: fxst,fopn
      character strn*80
      integer procid,mpipid
      integer, parameter :: master=0
      real(8), parameter :: tol = 1d-8

      lrspace = mod(iabs(job)/10,10) == 0
      job0 = isign(1,job)*mod(iabs(job),10)
      if (job0 == 0) return ! Do nothing
      k1 = n1; k2 = n2; k3 = n3
      ng = s_lat%ng
      procid = mpipid(1)
      pi   = 4d0*datan(1d0)
      pi8  = 8*pi
      tpiba=2*pi/s_lat%alat
      if (associated(s_pot%smcv0)) then; ncv0 = size(s_pot%smcv0,1); endif

C --- Gather ves(G) and smrho(G) before estimate for screening charged added to smrho
      if (job0 /= -1 .and. job0 /= -2) then
        if (lrspace) call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,-1)
        if (lrspace) call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)
        allocate(cv0(ng,2))
        call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smpot,cv0)
        call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smrho,cv0(1,2))
        if (lrspace) call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1) ! restore density to r space
      endif

C --- Add estimate for screening charge given eps0 ---
C     Add drho = R dphi where R = v^-1 (eps0^-1 - 1), dphi = s_pot%smvextc
C     Jobs -2, -6, -7 have 2's bit on; jobs -4, -6, -8 do not
      if (job0 < 0 .and. iand(iabs(job0),2) == 2 .and. eps0 /= 1) then
        call info2(15,1,0,' smvxte: add estimated screening charge using eps=%d',eps0,2)
        call info0(15,0,0,'    ig       kv%8fpi8/q2%13frho0%17fdrho%18fdphi')
        if (lrspace) call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)
        allocate(cn(ng))
        call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smrho,cn)
        do  ig  = 2, ng
          if (abs(s_pot%smvextc(ig)) < tol) cycle
          q2 = tpiba*tpiba*(s_lat%gv(ig,1)**2+s_lat%gv(ig,2)**2+s_lat%gv(ig,3)**2)
          fac = q2/pi8*(1/eps0-1)
C         (for prettified output)
          call dcopy(2,cn(ig),1,zx(1,1),1) ! rho0
          call dcopy(2,fac*s_pot%smvextc(ig),1,zx(1,3),1) ! drho
          call dcopy(2,s_pot%smvextc(ig),1,zx(1,2),1) ! external potential vext
          do  ifi = 1, 3
            if (abs(zx(1,ifi)) < tol) zx(1,ifi) = 0
            if (abs(zx(2,ifi)) < tol) zx(2,ifi) = 0
          enddo
          call info8(15,0,0,'%,6i%3,4i %;12,6D  %2:1;9F  %2:1;9F  %2:1;9F',ig,s_lat%kv(ig,:),
     .      pi8/q2, zx(1,1), zx(1,3), zx(1,2), 7, 8)
          cn(ig) = cn(ig) + fac*s_pot%smvextc(ig)
        enddo
        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cn,smrho)
        if (lrspace) call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
        deallocate(cn)
      endif
      if (job0 == -2) return

C --- Attempt to read s_pot%smcv0(:,1) from disk ---
      if (job0 == -7 .or. job0 == -9) then
        lerr = .true.  ! Will be set to .false. if read is successful
        if (fxst('smpot0') == 0) goto 20
        ifi = fopn('smpot0')
        rewind ifi
        read(ifi,"(6x,i6)",err=20) icv0
        if (icv0 /= ncv0) goto 20
        do  icv0  = 1, ncv0
          read(ifi,'(i6,2f18.12,2x,2f18.12)',err=20) ig, s_pot%smcv0(icv0,1), s_pot%smcv0(icv0,2)
          if (ig /= icv0) goto 20
        enddo
        lerr = .false.
        call info2(30,0,0,' smvxte: read %i 0-order potential elements from file smpot0',ncv0,2)
   20   continue ! Recovery point if file to read
        if (lerr) then
          call info0(30,0,0,' smvxte: failed to read 0-order potential from smpot0 ... use smpot')
          job0 = -4
        endif
        if (fxst('smpot0') /= 0) call fclose(ifi)
      endif

C --- This loop: ---
C    -job0=1     Counts number of nonzero points vext
C    -job0=4,6,8 Copies smpot(G) to s_pot%smcv0(:,1)
C     job0=1     Prints out estimate for eps
C     job0=2     Prints out estimate for for response function
      if (job0 == 1) then
        call info0(15,0,0,'%4fig%7fkv%7fpi8/q2%13fdphitot%15fphiext%15feps^-1')
      elseif (job0 == 2) then
        call info0(15,0,0,'%4fig%7fkv%7fpi8/q2%14fdrho%16fdphitot%16fphiext%16fR')
      elseif (job0 == -8) then
        call info0(15,0,0,'%4fig%7fkv%7fpi8/q2%14fphi0%18frho0%16fphiext')
      endif
      icv0 = 0
      do  ig  = 2, ng
        if (abs(s_pot%smvextc(ig)) < tol) cycle
        icv0 = icv0+1
        if (-job0 == 1) cycle
        if (icv0 > ncv0) call rx('s_pot%smcv0 wrongly dimensioned')
        if (-job0 == 4 .or. -job0 == 6 .or. -job0 == 8) then
          s_pot%smcv0(icv0,1) = cv0(ig,1)
          s_pot%smcv0(icv0,2) = cv0(ig,2)
        endif
        if (job0 < 1 .and. job0 /= -8) cycle
        dphitot = cv0(ig,1)-s_pot%smcv0(icv0,1)
        drho = cv0(ig,2)-s_pot%smcv0(icv0,2)
        epsi = dphitot/s_pot%smvextc(ig)
        q2 = tpiba*tpiba*(s_lat%gv(ig,1)**2+s_lat%gv(ig,2)**2+s_lat%gv(ig,3)**2)

C       (for prettified output)
        call dcopy(2,dphitot,1,zx(1,1),1) ! total phi
        call dcopy(2,s_pot%smvextc(ig),1,zx(1,2),1) ! external potential vext
        call dcopy(2,epsi,1,zx(1,3),1) ! inverse epsilon
        call dcopy(2,drho,1,zx(1,4),1) ! total rho
        call dcopy(2,drho/s_pot%smvextc(ig),1,zx(1,5),1) ! rho/vext
        if (job0 == -8) then
          call dcopy(2,cv0(ig,1),1,zx(1,1),1) ! smpot0
          call dcopy(2,cv0(ig,2),1,zx(1,4),1) ! smrho
        endif
        do  ifi = 1, 5
          if (abs(zx(1,ifi)) < tol) zx(1,ifi) = 0
          if (abs(zx(2,ifi)) < tol) zx(2,ifi) = 0
        enddo
        if (job0 == 1) then
        call info8(15,0,0,'%,6i%3,4i %;12,6D  %2:1;9F  %2:1;9F  %2:1;9F',ig,s_lat%kv(ig,:),
     .    pi8/q2, zx(1,1), zx(1,2), zx(1,3), 7, 8)
        elseif (job0 == 2) then
          call info8(15,0,0,'%,6i%3,4i %;12,6D  %2:1;9F  %2:1;9F  %2:1;9F  %2:1;9F',ig,s_lat%kv(ig,:),
     .      pi8/q2, zx(1,4), zx(1,1), zx(1,2), zx(1,5), 8)
        elseif (job0 == -8) then
          call info8(15,0,0,'%,6i%3,4i %;12,6D  %2:1;9F  %2:1;9F  %2:1;9F',ig,s_lat%kv(ig,:),
     .      pi8/q2, zx(1,1), zx(1,4), zx(1,2), 7, 8)
        endif
      enddo
      if (job0 == -1) then
        call info(40,1,0,' smvxt: found %i nonzero external V_G',icv0,0)
        n1 = icv0
        return
      endif
      if (icv0 /= ncv0) call rx('s_pot%smcv0 wrongly dimensioned')

C --- Write s_pot%smcv0(:,1) to disk ---
      if (job0 == -8) then
        ifi = fopn('smpot0')
        rewind ifi
        write(ifi,"('% rows',i6)") ncv0
        do  icv0  = 1, ncv0
          write(ifi,'(i6,2f18.12,2x,2f18.12)') icv0, s_pot%smcv0(icv0,1), s_pot%smcv0(icv0,2)
        enddo
        call fclose(ifi)
        call awrit1('smvxte: wrote %i potential elements to file smpot0',strn,len(strn),0,ncv0)
        call rx0(trim(strn))
      endif

C --- Restore potential to r space ---
      if (job0 /= -1 .and. job0 /= -2) then
        if (lrspace) call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
        deallocate(cv0)
      endif

      end
