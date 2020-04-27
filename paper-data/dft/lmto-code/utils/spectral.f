      subroutine fmain
C- Generate spectral functions A(omega) and optionally generate DOS
C Required inputs:
C     SEComg.UP (DN)
C     QIBZ      required for either of:
C               dos is to be integrated with proper q weights
C               --ws switch (save binary file for further processing)
C     TOTE.UP  (DN) : only for oneshot=1 case
Cl Local variables
Cl   sigqp : sigma(eqp)
Cl   sigw  : sigma(omega)
Cl   eps   : smearing (eV)
Co Outputs
Co   Files written to disk (switch --ws NOT used):
Co     Spectral functions sec_ib#_iq#.up (sec_ib#_iq#.dn)
Co     q-summed DOS dos.up, dos.dn
Co   If --ws is used, file se is generated.
Cr Remarks:
Cr  Examples:
Cr  Make spectral function files for 1st kp and QP levels above 10 eV
Cr    spectral --eps=.005 --domg=0.00340145 '--cnst:iq==1&eqp>-10'| tee out
Cr  Compare QP and spectrum DOS
Cr    fplg -y 0,2 -colsy 3 -lt 1,col=.5,.5,.5 dos.up -colsy 2 -lt 1,col=0,0,1 dos.up
Cr  Plot Real, Im parts of SE for one band and QP
Cr    fplot -disp -frme:xor=0:yab=0 0,1,0,1 -colsy 3,4 sec_ib2_iq1.up
Cr  Plot spectral function for G and G0
Cr    fplot -disp -frme:xor=0:yab=0 0,1,0,1 -y 0,2  -colsy 7,8 sec_ib2_iq1.up
Cr  Make file 'se' suitable for lmfgws
Cr    spectral --ws --nw=1 '--cnst:ib>=3&ib<=12'
Cu Updates
Cu   30 Nov 16 new --ws:x
C --------------------------------------------------
      implicit none
      character(3)  :: sss
      integer:: iarg,ifiq,nwmin,nwmax,ifidos,ifise,iw_r,is_r,
     .  iw0,iw,is,ib_r,iq_r,iq,ib,ndat,nqibz,iwx,iws,nwx,im,oneshot,iv0,
     .  ival,nbandx,nqibzx,iftote,ib_rx,fopnx,fopng,i1mach,stdo,nq,
     .  nband,garg,jx1,jx2,npoly=6,ix(2),nb0,isum,nomg,nsp,nw,NULLI
      integer nspinx,nqx,ntqx,ntq,isx
      parameter (NULLI=-99999)
      real(8):: qp(3),eig_r,se_r,omg_r,egw,sigqp,ein,s1,s2,sumchk,sumchk2,
     .  eee,egw1,Eqp,xx,ddot,dsum,ommin,ommax,dw,dw0,pi,aw,eps,
     .  dymx=1e-5,dy,sigwr,sigwi,range(2),ef,deltaw,
     .  q_rx(3)=-9999d0
      real(8),parameter :: tolq=1d-7
C     double precision polinta
      complex(8):: den
      logical :: a2bin, cmdstr, readqibz=.false., ldebug=.false.,
     .  lws=.false., lwsx=.false.
      integer,allocatable:: ibcount(:,:)
      real(8),allocatable:: omg(:),wibz(:),qibz(:,:),eig(:,:,:),
     .  ez1(:,:),dos(:),omgx(:),dosqp(:),sumag(:),sumag0(:),ag(:),
     .  ag0(:),qpf(:,:),sigqpi(:,:,:),sex(:,:,:)
      complex(8), allocatable:: sigm(:,:,:),sigw(:),sigmi(:,:,:,:)
      character fname*20, fnD*20, ext(2)*2, extl(2)*2
      character(130) fnA, first, outs, expr
      procedure(integer) :: iprint
      procedure(real(8)) :: dlength

C --- Parameter declarations and defaults ---
C     multi: spacing of given self-energy mesh is subdivided by
C            this factor, for finer resolution of A(w) and dos(w)
C
      integer :: multi=40
      data ext /'UP','DN'/
      data extl /'up','dn'/

C --- Parse command-line arguments ---
      pi = 4*datan(1d0)
      stdo = i1mach(2)
      oneshot = 0
      expr = ' '
      iarg = 1
      nw = 0
      dw = 0
      eps = 0
      range(1) = -9999
      range(2) =  9999
C     Evaluate a trivial expression to load the basic constants
      iw = 0
      if (.not. a2bin('1 ',ival,2,0,' ',iw,-1)) then
       call rx('spectral: problem with compilation')
      endif
      call setpr(30)

   15 continue
      first = ' '
      if (.not. cmdstr(iarg,first)) goto 100
      if (ldebug) call awrit0(' mc: parsing '//first,outs,-80,stdo)

      if (first(1:7) == '--1shot') then
        oneshot = 1

      else if (first(1:4) == '-nw=') then
        call rxx(garg('-nw=',iarg,2,' ',1,1,im,ib,nw) < 1,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:5) == '--nw=') then
        call rxx(garg('--nw=',iarg,2,' ',1,1,im,ib,nw) < 1,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:6) == '-domg=') then
        call rxx(garg('-domg=',iarg,4,' ',1,1,im,ib,dw) < 1,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:7) == '--domg=') then
        call rxx(garg('--domg=',iarg,4,' ',1,1,im,ib,dw) < 1,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:6) == '--eps=') then
        call rxx(garg('--eps=',iarg,4,' ',1,1,im,ix,eps) < 1,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:6) == '--ws:x') then
        lws = .true.
        lwsx = .true.

      else if (first(1:3) == '-pr' .or. first(1:4) == '--pr') then
        ib = 3; if (first(1:4) == '--pr') ib=4
        if (a2bin(first,iw,2,0,' ',ib,-1)) then
          call pshpr(iw)
        endif

      else if (first(1:4) == '--ws') then
        lws = .true.

      else if (first(1:8) == '--range=') then
        call rxx(garg('--range=',iarg,4,', ',3,2,im,ix,range) < 2,
     .    'spectral:  cannot parse expression '//trim(first))

      else if (first(1:2) == '-v') then
        iw = 2
        call parsyv(first,len(first),999,0,iw)

      else if (first(1:7) == '--cnst:') then
        expr = trim(first(7+verify(first(8:), ' ',.false.):))

      else
        print *, 'usage: spectral [switches]'
        print 333
  333   format(" Generates spectral and individual QP contributions",
     .    " to the Green's function."/
     .    " A separate file is generated for each QP, and a separate",
     .    " DOS file."//
     .    "   Switches:"//
     .  4x,"--1shot :",T19,"evals do not correspond with QP levels"//
     .  T6,"Choose either nw or domg."/
     .  T6,"Note: spectral function is only reliable on a fine energy",
     .      " mesh."/T6,"When --ws is used, --nw=1 is advised."/
     .  4x,"   -nw=#  :"/
     .  4x,"  --nw=#  :",T19,"Subdivide input energy mesh nw times"/
     .  4x," -domg=# :"/
     .  4x,"--domg=# :",T19,"Target energy spacing for A(omega), eV"/
     .  T19,"domg chooses the nw that approximately yields target",
     .      " spacing."//
     .  4x,"--range=#1,#2:",T19,"restrict generated files ",
     .      "to #1 < omega < #2"//
     .  4x,"--eps=# :",T19,"smearing width, eV"//
     .  4x,"--ws :",T19,"Write a self-energy file (se) for all q."/
     .              T19,"Individual files are not writen."/
     .  4x,"--ws:x :",T19,"Same as --ws, but also write exchange potential to file."//
     .  4x,"--cnst:expr : Exclude entries in SEComg.(UP,DN)",
     .     " for which expr is nonzero."/
     .  T19,"'expr' is an integer expression that can include",
     .  " the following variables:"/
     .  T21,"ib (band index)"/
     .  T21,"iq (k-point index)"/
     .  T21,"qx,qy,qz,q (components of q, and amplitude)"/
     .  T21,"eqp (quasiparticle level)"/
     .  T21,"spin (1 or 2)"/
     .  T19,"Example: Use only the first k-point ",
     .      "and exclude bands 1..3"/
     .  T19,"--cnst:iq==1&eqp>-10")
        call cexit(0,1)
      endif

C ... End of input
      iarg = iarg+1
      goto 15
  100 continue

C --- Read QIBZ ---
      if (fopnx('QIBZ',2,-1,-1) == 1) then
        ifiq = fopnx('QIBZ',2,1,-1)
C       open (ifiq, file='QIBZ',status='old',err=998) !q-points in IBZ.
        read(ifiq,*,err=998,end=998) nqibz
        allocate(qibz(1:3,nqibz),wibz(nqibz))
        do  iq = 1, nqibz
          read(ifiq,"(3d24.16,3x,d24.16)") qibz(1:3,iq),wibz(iq)
          wibz(iq) = wibz(iq)/2  ! so that sum of wibz is 1
        enddo
        call info2(10,0,0,' spectral: read %i qp from QIBZ',nqibz,0)
        if (abs(dsum(nqibz,wibz,1)-1) > 1d-10)
     .    call rx('BZ weights do not sum to 1')
        readqibz = .true.
      endif
  998 continue
      if (.not. readqibz) then
        call info0(10,0,0,' spectral: missing QP or data from QIBZ')
        call rxx(lws,'QIBZ must be present when using --ws')
      endif

C ... Get array dimensions; allocate arrays
      nwmin = 9999
      nwmax =-9999
      ommin = 9999
      ommax =-9999
      iw0   = 9999
      nq = 0
      nband = 0
      dw0 = 0
      nsp = 0
      call numsyv(iv0)
      do  is = 1, 2
        if (fopnx('SEComg.'//ext(is),2,-1,-1) == 1) then
          nsp = is
          ifise = fopnx('SEComg.'//ext(is),2,1,-1)
          rewind ifise
          do while (.true.)
            read(ifise,*,end=110)
     .        iw_r,ib_r,iq_r,is_r,qp(1:3),eig_r,omg_r

            if (expr(1:1) /= ' ') then
              call numsyv(iv0)
              call lodsyv('iq',1,dble(iq_r),ival)
              call lodsyv('ib',1,dble(ib_r),ival)
              call lodsyv('spin',1,dble(is_r),ival)
              call lodsyv('qx', 1,qp(1),ival)
              call lodsyv('qy', 1,qp(2),ival)
              call lodsyv('qz', 1,qp(3),ival)
              xx = sqrt(ddot(3,qp,1,qp,1))
              call lodsyv('q', 1,xx,ival)
              call lodsyv('eqp', 1,eig_r,ival)
C             call shosyv(0,0,0,6)
              iw = 0
              call rxx(.not. a2bin(expr,ival,2,0,' ',iw,len(expr)),
     .          'spectral:  cannot parse expression '//trim(expr))
              call clrsyv(iv0)
              if (ival == 0) cycle
            endif
            nwmin = min(nwmin,iw_r)
            nwmax = max(nwmax,iw_r)
            ommin = min(ommin,omg_r)
            ommax = max(ommax,omg_r)
            nq = max(nq,iq_r)
            nband = max(nband,ib_r)

C           Find spacing between energy points
            if (iw_r == nwmin+1) dw0 = omg_r - ommin
          enddo
  110     continue
        else
          call info0(10,0,0,' no file SEComg.'//ext(is))
          call rxx(is == 1.and.lws,
     .      ' file SEComg.UP required for --ws mode')
        endif
      enddo
      call info8(10,0,0,'%1fDimensions from file(s) SEComg.(UP,DN):%N'//
     .  ' nq=%i  nband=%i  nsp=%i  omega interval (%;6d,%;6d) eV'//
     .  ' with (%i,%i) points',nq,nband,nsp,ommin,ommax,nwmin,nwmax,0)
      if (nq == 0 .or. nband == 0 .or. nsp == 0) then
        print *, 'Nothing to calculate ... quitting'
        call cexit(0,1)
      endif

C     Set q-integration weights to unity, if not supplied
      if (.not. allocated(wibz)) then
        allocate(wibz(nq))
        call dvset(wibz,1,nq,1d0)
      endif

C ... Set mesh spacing
      if (nw /= 0) then
        multi = nw
        dw = dw0
        if (eps == 0) eps = .005d0
      else if (dw == 0) then
        if (eps == 0) eps = .005d0
        multi = max(2*dw0/eps,1d0)
      else
        multi = max(dw0/dw,1d0)
        if (eps == 0) eps = 2*dw
      endif
      call info5(10,0,0,'%1fEnergy mesh spacing = %;1d meV'//
     .  ' ... interpolate to target spacing %;1d meV.  '//
     .  'Broadening = %;1d meV',dw0*1e3,dw0/multi*1e3,eps*1000,0,0)

C --- For each  ib, iq, spin, read self-energy, make A ---
      allocate(omg(nwmin:nwmax),qpf(3,nq),eig(nband,nq,nsp),
     .  ez1(nband,nq),sigm(nwmin:nwmax,nband,nq))
      call dvset(eig,1,nband*nq*nsp,1d10) ! Flags evals not read
      allocate(ibcount(nband,2))
      nb0 = 0
      do  is = 1, nsp
        call iinit(ibcount(1,is),nband)
        if (fopnx('SEComg.'//ext(is),2,-1,-1) /= 1) cycle
        ifise = fopnx('SEComg.'//ext(is),2,1,-1)
        rewind ifise
C   ... Read until no more data
        do while (.true.)
          read(ifise,*,end=210)
     .      iw_r,ib_r,iq_r,is_r,qp(1:3),eig_r,omg_r,s1,s2
          if (is_r /= is) call rx('file spin does not match name')
          if (expr(1:1) /= ' ') then
            call numsyv(iv0)
            call lodsyv('iq',1,dble(iq_r),ival)
            call lodsyv('ib',1,dble(ib_r),ival)
            call lodsyv('spin',1,dble(is_r),ival)
            call lodsyv('qx', 1,qp(1),ival)
            call lodsyv('qy', 1,qp(2),ival)
            call lodsyv('qz', 1,qp(3),ival)
            xx = sqrt(ddot(3,qp,1,qp,1))
            call lodsyv('q', 1,xx,ival)
            call lodsyv('eqp', 1,eig_r,ival)
C           call shosyv(0,0,0,6)
            iw = 0
            call rxx(.not. a2bin(expr,ival,2,0,' ',iw,len(expr)),
     .        'spectral:  cannot parse expression '//trim(expr))
            call clrsyv(iv0)
            if (ival == 0) cycle
          endif

          ibcount(ib_r,is) = ibcount(ib_r,is)+1
          omg(iw_r) = omg_r
          eig(ib_r,iq_r,is) = eig_r
          sigm(iw_r,ib_r,iq_r) = dcmplx(s1,s2)
          qpf(1:3,iq_r) = qp
C         Confirm that qp(iq_r) matches QIBZ(iq_r); else reset
          if (readqibz) then
            if (sum(abs(qp - qibz(:,iq_r))) > 3d-6) then
              call info0(10,0,0, '%1f(warning) '//
     .          'file mismatch with QIBZ ... discard QIBZ weights')
              call info2(10,0,0, '%1ffile contains q ='//
     .          '%3;10,5D',qp,0)
              call dvset(wibz,1,size(wibz),1d0)
              readqibz = .false.
              call rxx(lws,'qp must synchronize when using --ws')
            endif
          endif

C     ... If one-shot, get eQP(Z=1) from TOTE
          if (oneshot == 1) then
            call rx('ib_rx not defined')
            if (sum(abs(q_rx-qp))<1d-6 .and. ib_rx==ib_r) then
              goto 120
            endif
            iftote = 1016
            fname = 'TOTE'//sss
            open(iftote, file=fname)
            read(iftote,*) nqibzx,nbandx
            do  iq = 1, nqibzx
            do  ib = 1, nbandx
              read(iftote,*) q_rx, ib_rx, iq_r, eee,egw,egw1
              if (sum(abs(q_rx-qp))<1d-6 .and. ib_rx==ib_r) then
                ez1(ib_r,iq_r) = egw1
                goto 119
              endif
            enddo
            enddo
            call rx('can not find Eqp in TOTE.*')
  119       continue
            close(iftote)
  120       continue
          endif
        enddo
  210   continue
        call fclose(ifise)
        ndat = nwmax-nwmin+1

C   ... Find first band for which there is data
        do  iq = 1, nq
        do  ib = 1, nband
          if (eig(ib,iq,is) /= 1d10) goto 220
          nb0 = ib
        enddo
        enddo
  220   continue
        iq = isum(nband,ibcount,1)
        if (mod(iq,ndat) /= 0) call info2(10,0,0,
     .        '%1f(warning) omega mesh not uniform for all QPs',0,0)
        call info5(10,0,0,'%1fSpectral functions starting from '//
     .    'band %i, spin %i, for %i QPs',nb0+1,is,iq/ndat,0,0)
        if (lws .and. iq /= ndat*nqibz*(nband-nb0)) then
          call info0(10,0,0,
     .    '%1f(warning) missing data for some qp, bands, or omega')
C         call rx('Incomplete SEComg file')
        elseif (iq /= ndat*nq*(nband-nb0)) then
          call info0(10,0,0,
     .    '%1f(warning) missing data for some qp, bands, or omega')
        elseif (is == 2 .and. iq /= isum(nband,ibcount(1,2),1)) then
          call rx(' spin mismatch in bands ... abort')
        endif

C   ... Number of energy mesh points, mesh range
        nomg = 0
        ommin = 9999
        ommax =-9999
        iwx = 0
        do  iw = nwmin, nwmax-1
          do  im = 0, multi
            if (im == multi .and. iw /= nwmax-1) cycle
            iwx = iwx+1
            xx = omg(iw) + (omg(iw+1)-omg(iw))*im/dble(multi)
            if (xx >= range(1) .and. xx <= range(2)) then
              ommin = min(ommin,xx)
              ommax = max(ommax,xx)
              nomg = nomg + 1
            endif
          enddo
        enddo

C   --- Generate and write dependent spectral function file for each ib,iq ---
C       Also accumulate total DOS
C       omgx: omega on finer energy mesh (sigma interpolated to it)
        if (.not. allocated(dos)) then
        allocate(dos(ndat*multi),dosqp(ndat*multi),omgx(ndat*multi),
     .           sigw(ndat*multi),sumag(ndat*multi),sumag0(ndat*multi))
        if (lws) then
         allocate(sigmi(nomg,nband-nb0,nq,nsp),sigqpi(nband-nb0,nq,nsp))
         if (lwsx) allocate(sex(nband-nb0,nq,nsp))
         if (.not. lwsx) allocate(sex(1,1,1))
        endif
        endif
        call dpzero(dos,ndat*multi)
        call dpzero(dosqp,ndat*multi)
        if (lws) then
          call info0(10,1,0,'    ib    iq%6fEqp%28pint A(G)   '//
     .      'int A(G0) rat[G] rat[G0]')
        else
          call info0(10,1,0,'%10ffile%26pEqp%35pint A(G)   '//
     .      'int A(G0) rat[G] rat[G0]')
        endif
C       call imxmn(nband,ibcount,1,iq,ib)
        do  iq = 1, nq
        do  ib = 1, nband
          if (eig(ib,iq,is) == 1d10) then
            if (expr(1:1) == ' ') call info2(10,0,0,
     .        '%1f(warning) nothing read for ib=%i iq=%i',ib,iq)
            cycle ! Nothing read for this band
          endif

          if (oneshot == 1) then
            Eqp = ez1(ib,iq)
          else
            Eqp = eig(ib,iq,is)
          endif

C         Z=1 case :  E_qp = eig + (sigm_x + sigm(eig) - vxc)
C     ... Make sigqp = sigm(eig)
          call polint(omg,dble(sigm(nwmin:nwmax,ib,iq)),ndat,
     .      npoly,eig(ib,iq,is),dymx,3,jx1,sigqp,dy)

C     ... Generate sumag = 1/pi int[Im G(w)] and 1/pi int[Im G0(w)]
          iwx = 0
          sumchk =  0
          sumchk2 = 0
          jx1 = 0
          do  iw = nwmin, nwmax-1
            do  im = 0, multi
            if (im == multi .and. iw /= nwmax-1) cycle
            iwx = iwx+1
            omgx(iwx) = omg(iw) + (omg(iw+1)-omg(iw))*im/dble(multi)
            dw = (omg(iw+1)-omg(iw))/dble(multi)
            if (oneshot == 1) then
              ein = omgx(iwx) + eig(ib,iq,is) - ez1(ib,iq)
            else
              ein = omgx(iwx)   ! omega is true QP energy
            endif

C       ... Make Sigma(omgx)
c           Use polinta to return f(ein) = polinta(ein,  x(1..n), f(1...n), n)
C            sigw = dcmplx(polinta(ein,omg,
C     .               dble(sigm(nwmin:nwmax,ib,iq)),nwmax-nwmin+1),
C     .                    polinta(ein,omg,
C     .              dimag(sigm(nwmin:nwmax,ib,iq)),nwmax-nwmin+1))
            call polint(omg,dble(sigm(nwmin:nwmax,ib,iq)),
     .         nwmax-nwmin+1,npoly,ein,dymx,3,jx1,sigwr,dy)
            call polint(omg,dimag(sigm(nwmin:nwmax,ib,iq)),
     .        nwmax-nwmin+1,npoly,ein,dymx,3,jx1,sigwi,dy)
            sigw(iwx) = dcmplx(sigwr,sigwi)

C       ... Denominator = omega - (E_qp + sigm(omega) - sigma(E_qp))
            den = dcmplx(dble(omgx(iwx) - (Eqp + sigw(iwx) - sigqp)),
     .                   abs(dimag(sigw(iwx)))+eps)! eps is smearing

C       ... Integrate by trapezoidal rule A(omgx) = abs(dimag(1d0/den))
C           and accumulate spectrum DOS
            aw = abs(dimag(1/pi/den))
            sumchk = sumchk + dw*aw
            sumag(iwx) = sumchk
C           dos(iwx) = dos(iwx) +  wibz(iq) * aw
C           if (iwx == 1) print *, iq,wibz(iq)

C       ... Integrate by trapezoidal rule QP DOS
            den = dcmplx(dble(omgx(iwx) - Eqp), eps  ) !eps = smearing
            aw = abs(dimag(1/pi/den))
C           dosqp(iwx) = dosqp(iwx) + wibz(iq) * aw
            sumchk2 = sumchk2 + dw*aw
            sumag0(iwx) = sumchk2

          enddo ! Frequency subdivision loop
          enddo ! Frequency loop
          nwx = iwx  ! dimension of fine frequency mesh

C     ... Generate A[G(w)] and A[G0(w)] by numerical differentiation
C         jx=0 => rational interpolation (attempt first), else polynomial
          allocate(ag(ndat*multi),ag0(ndat*multi))
          call poldvm(omgx,sumag,nwx,6,.true.,1d-6,jx1,ag)
          if (jx1 /= 0) then
            call poldvm(omgx,sumag,nwx,6,.false.,1d-6,iw,ag)
          endif
          call poldvm(omgx,sumag0,nwx,6,.true.,1d-6,jx2,ag0)
          if (jx2 /= 0) then
            call poldvm(omgx,sumag0,nwx,6,.false.,1d-6,iw,ag0)
          endif

C     ... Write sigma, int A[G(w)]  and A[G0(w)] to file; accumulate dos
          if (.not. lws) then
            call awrit2('sec_ib%i_iq%i.'//extl(is),fnA,len(fnA),0,ib,iq)
            ifise = fopng(trim(fnA),-1,0)
C           call fshow
            write(ifise,321) ib,iq,Eqp,qpf(1:3,iq)
  321       format('# ib=',i4,'  iq=',i4,'  Eqp=',f12.6:'  q=',3f11.6
     .       /'#',5x,'omega',9x,'omega-Eqp',5x,'Re sigm-vxc',4x,
     .        'Im sigm-vxc',6x,'int A(w)',6x,'int A0(w)',
     .        7x,'A(w)',11x,'A0(w)')

            iwx = 0
            do  iw = nwmin, nwmax-1
            do  im = 0, multi
              if (im == multi .and. iw /= nwmax-1) cycle
              iwx = iwx+1

              if (omgx(iwx) >= range(1).and.omgx(iwx) <= range(2)) then
               write(ifise,"(' ',15d15.7)") omgx(iwx), omgx(iwx)-Eqp,
     .         (sigw(iwx)-sigqp),sumag(iwx),sumag0(iwx),ag(iwx),ag0(iwx)
              endif
              dos(iwx) = dos(iwx) +  wibz(iq) * ag(iwx)
              dosqp(iwx) = dosqp(iwx) +  wibz(iq) * ag0(iwx)

            enddo ! Frequency subdivision loop
            enddo ! Frequency loop

            write(stdo,"(a20,f12.6,2f11.4,2L6)")
     .        trim(fnA), Eqp,sumchk,sumchk2, jx1 == 0, jx2 == 0
            call fclr(trim(fnA),ifise)

C         Accumulate interpolated sigma for se file, range of interpolated sigma
          else

C           Integrated data written to stdo
            if (iprint() > 30) then
              write(stdo,"(2i6,f12.6,2f11.4,2L6)")
     .          ib, iq, Eqp,sumchk,sumchk2, jx1 == 0, jx2 == 0
            endif

            iwx = 0
            iws = 0
            do  iw = nwmin, nwmax-1
              do  im = 0, multi
              if (im == multi .and. iw /= nwmax-1) cycle
              iwx = iwx+1
              if (omgx(iwx) >= range(1).and.omgx(iwx) <= range(2)) then
                iws = iws+1
C               sigmi(iws,ib-nb0,iq,is) = sigw(iwx)-sigqp
                sigmi(iws,ib-nb0,iq,is) = sigw(iwx)
                sigqpi(ib-nb0,iq,is) = sigqp
              endif
            enddo
            enddo
            call rxx(iws /= nomg,'dimension mismatch')
          endif
          deallocate(ag,ag0)
        enddo
        enddo

        if (lwsx) then
          fname = 'SEXU'; if (is == 2) fname = 'SEXD'
          ifise = fopng(fname,-1,1)  ! open as OLD
          call readx(ifise,50)
          read (ifise,*) nspinx,nqx,ntqx
          if (nsp /= nspinx)  call rx('spectral: wrong nspin SEXU')
          if (nq /= nqx)        call rx('spectral: wrong nq SEXU')
          if (is == 1) ntq = ntqx
          if (ntq /= ntqx)      call rx( 'hqpe: wrong ntq SEx')
          read (ifise,*)
          read (ifise,*) deltaw
          read (ifise,*)
          read (ifise,*) ef
          call readx(ifise,50) ! move file pointer to start of exchange
          read (ifise,*)

          do  iq = 1, nq
            do  ib  = 1, ntq
              read(ifise,"(3i5,3d24.16,3x,d24.16,3x,d24.16)") ib_r,iq_r,isx,q_rx,eig_r,se_r
              if (ib_r /= ib) call rx('band mismatch, file '//trim(fname))
              if (dlength(3,q_rx(:)-qpf(:,iq),1) > tolq) call rx('k-point mismatch, file '//trim(fname))
              if (ib <= nb0 .or. ib > nband) cycle
              sex(ib-nb0,iq,is) = se_r
            enddo
          enddo
        endif

C   --- Write spectral function file ---
        if (lws .and. is == nsp) then
          call info8(10,1,0,'%1fWriting SE file:  '//
     .      'nq=%i nband=%i nsp=%i, %i points in (%;6d,%;6d) eV window',
     .      nq,nband,nsp,nomg,ommin,ommax,0,0)
          ifise = fopng('se',-1,0)
          iw = 10; if (lwsx) iw = 30
          call ioseh(120,iw,nsp,nsp,nband-nb0,nq,nomg,ommin,ommax,dble(NULLI),0,[0],[0],-ifise)
C         Debugging check
C         call ioseh(120,0,nsp,nsp,nband-nb0,nq,nomg,ommin,ommax,dble(NULLI),4,[1,5,6,7],[3, 5,234, -1],-ifise)
          do  iq = 1, nq
            do  ib = 1, nband
              if (eig(ib,iq,1) == 1d10) eig(ib,iq,1) = NULLI
              if (eig(ib,iq,nsp) == 1d10) eig(ib,iq,nsp) = NULLI
            enddo
          enddo
          call iose2(20,iw,iw,nsp,nband-nb0,nq,nomg,qpf,
     .      eig(nb0+1:nband,:,:),sigqpi,sigmi,sex,sex,-ifise)
          call fclose(ifise)
        endif
        if (lws) cycle

C   --- Write DOS file ---
        fnD = 'dos.'//extl(is)
        if (readqibz) then
          call info0(10,1,0,' writing q-integrated dos to file '//
     .      trim(fnD)//' ...')
        else
         call info0(10,1,0,' writing partial q-summed dos to file '//
     .      trim(fnD)//' ...')
        endif
        ifidos = fopng(trim(fnD),-1,0)
        do  iwx = 1, nwx
          if (omgx(iwx) >= range(1) .and. omgx(iwx) <= range(2)) then
            write(ifidos,"(12d16.7)") omgx(iwx),dos(iwx),dosqp(iwx)
          endif
        enddo
        call fclr(trim(fnD),ifidos)
        deallocate(dos,dosqp,omgx,sigw,sumag,sumag0)
        cycle

  999   write(stdo,*) 'no file '// trim(fname)// ' ... nothing written'
      enddo

      end

c These routines are taken from Ferdi's rw.f
c-----------------------------------------------------------------
      subroutine readx(ifil,n)
      character*72 rchar
      do 10 i = 1,n
        read(ifil,5000)rchar
        j       = 0
c        call rmvbl (rchar,72,j)
c        rchar      = rchar(j+1:72)
        rchar=trim(adjustl(rchar))

        if(rchar(1:3) == '***')return
        if(rchar(1:3) == '###')return

   10 continue
Cstop2rx 2013.08.09 kino      stop 'readx: cannot find the string(rw.f)'
      call rx( 'readx: cannot find the string(rw.f)')
 5000 format(a72)
c     return
      end
