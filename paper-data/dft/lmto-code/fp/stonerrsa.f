C#define DEBUG
      subroutine stonerrsa(nq,nw,nmbas,qp,momsite,mmnorm,eiqrm,freq,
     .  x0et)
C- Transverse susceptibility matrix <e~| X |e~> from X^{-1}=X0^{-1}+U
C ----------------------------------------------------------------------
Ci Inputs
Ci   nq     :number of q points for which to calculate X(q,w)
Ci   nw     :number of frequencies where input x0et is given
Ci   nmbas  :number of magnetic sites
Ci   qp     :qp(3,nq)  vector of each q
Ci   momsite:magnetic moment m
Ci   mmnorm :|m|=sqrt(<m|B> <B|m> ) = norm(svec) ; this is not m.
Ci          :e~_a(r)=m_a(r) / sqrt(\int m_a(r)**2 dr)  then <e~_a|e~_a>=1
Ci   eiqrm  :eiqrm_i = <e~_a|e^{iqr}> =  <M_a|eiqr>/sqrt(<m_a|m_a>)
Ci   freq   :frequency mesh
Ci   x0et   :<e~|X0|e~>
Co Outputs:
Co   ... finish this
Cl Local variables
Cr Remarks
Cr Define:
Cr   e_a(r)  = m_a(r)/M_a, where M_a = \int m_a(r) dr (unit norm)
Cr   e~_a(r) = m_a(r)/sqrt<m_a|m_a>, where <m_a|m_a> = \int m_a(r)**2 dr
Cr Thus
Cr   M_a |e_a> = sqrt<m_a|m_a> |e~_a>
Cr Also define
Cr   |eb_a> = M_a/sqrt<m_a|m_a> |e~_a>
Cr We have the following normalizations:
Cr   <e~_a|e~_a> = 1   and   <eb_a|e_a> = 1
Cu Updates
Cu   07 Feb 09 (L. Ke) adapted from T. Kotani
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nq,nw,nmbas
      real(8) momsite(nmbas), mmnorm(nmbas),freq(nw),qp(3,nq)
      complex(8) eiqrm(nmbas,nq),x0et(nmbas,nmbas,nq,nw)
C ... Local parameters
      integer fopna,fopnx,i,ifi,ifi2,imax,imin,iq,iw,iwpola,iwpolf,ipr,
     .  ix,iy,j,jmax,jmin,nglob,nw_intp,nwx,stdo, iferro
      complex(8),allocatable :: qxq_intp(:,:),wk(:,:),qx0q_intp(:,:),
     .  x0inv(:,:,:,:),xinv(:,:,:,:),x2et(:,:,:,:),
     .  xinvh(:,:,:,:),dxidw(:,:,:,:),dxidw_eb(:,:,:,:)
      real(8),allocatable:: mxevl_xinvh(:,:),mxevl2_xinvh(:,:),e_intp(:)
      real(8) uub(nmbas,nw),uu(nmbas,nw),eval(nmbas),uu0_eb(nmbas)
      real(8) freqm(nw),mmnorm2(nmbas)
      real(8) emin_intp,dw_intp,emax_intp
      complex(8),allocatable ::qxq(:,:), qx0q(:,:)
      real(8),allocatable :: qxq_r(:,:),qxq_i(:,:),
     .     qxq_inv_r(:,:),qxq_inv_i(:,:)
      complex(8):: cxtmp(nmbas,nmbas),img=(0d0,1d0),
     .  meffi(nmbas,nmbas,nq),xinvh_w0(nmbas,nmbas,nq),
     .  x0inv_w0eb(nmbas,nmbas,nq),
     .  meffi_eb(nmbas,nmbas,nq),xinvh_w0eb(nmbas,nmbas,nq),
     .  meffi2(nmbas,nmbas,nq),meffi2_eb(nmbas,nmbas,nq)
C      complex(8):: uumtx(nmbas,nmbas)
      real(8) rydberg,polinta
      real(8) epole_fm(nq),epole_af(nq),vpole_fm(nq),vpole_af(nq)
      real(8) rtmp,rtmp1,rtmp2,rymev,epole
      external :: polinta
C     For calculating JMAT:
      complex(8) oo(nmbas,nmbas),evc(nmbas,nmbas),sqm(nmbas,nmbas),
     .  jjmat(nmbas,nmbas),jjmat2(nmbas,nmbas),mxi(nmbas,nmbas)
      integer(4):: nmx,nev,nw_wrX0qw
      parameter (rydberg=13.6058d0)
c$$$      write(*,*) 'ok, after '
C --- Setup ---
c      nwx = 400
      nw_wrX0qw= 800
      if (nw < nw_wrX0qw)
     .     nw_wrX0qw = nw    ! the cutoff mesh should smaller than the one without being cut

      nwx = nw
      rymev = rydberg*1d3

      stdo = nglob('stdo')
      call getpr(ipr)

      if (freq(1) /= 0d0)
     .  call rx1('nonzero 1st frequency, w(1) = %g',freq(1))

      allocate(qxq(nw,nq))
      allocate(qxq_r(nw,nq),qxq_i(nw,nq))
      allocate( qxq_inv_r(nw,nq),qxq_inv_i(nw,nq) )

      allocate(x0inv(nmbas,nmbas,nq,nwx))
      allocate(xinv(nmbas,nmbas,nq,nwx),x2et(nmbas,nmbas,nq,nwx) )
      allocate(xinvh(nmbas,nmbas,nq,nwx),mxevl_xinvh(nq,nwx) )
      allocate(dxidw(nmbas,nmbas,nq,nwx))
      allocate(dxidw_eb(nmbas,nmbas,nq,nwx))
      allocate(mxevl2_xinvh(nq,nwx) )

c$$$      write(*,*) 'ok, iam here'
c$$$C#ifdef DEBUG
C     open(111, file='x0qw.allq',form='unformatted')
c     nw_wrX0qw = nw;
      ifi = fopnx('x0qw.allq',2,2+4,-1)
      write (ifi) nq,nw_wrX0qw,nmbas
      write (ifi) qp
      write (ifi) freq(1:nw_wrX0qw)
      do  iq = 1, nq
      do  iw = 1, nw_wrX0qw
        write(ifi) dble(x0et(1:nmbas,1:nmbas,iq,iw))
        write(ifi) dimag(x0et(1:nmbas,1:nmbas,iq,iw))
      enddo
      enddo
      call fclose(ifi)
! C#endif


C     open(111, file='qm.allq',form='unformatted')
      ifi = fopnx('eiqrmm_allq',2,2+4,-1)
      write (ifi) nq, nmbas
      write (ifi) mmnorm, momsite
      do  iq = 1, nq
        write(ifi) dble( eiqrm(1:nmbas,iq))
        write(ifi) dimag(eiqrm(1:nmbas,iq))
      enddo
      call fclose(ifi)
! C#endif



C --- Determine <e~|U(w)|e~> from boundary condition at q=0 ---
C     Eq. 51, Liqin's notes
      x0inv(:,:,1,:) = x0et(:,:,1,:) ! only need x0(q=0,w)
      mmnorm2 = mmnorm**2
      do  iw = 1, nwx
        call matcinv(nmbas, x0inv(:,:,1,iw)) ! x0inv=<e~|x0^-1|e~>
        uub(1:nmbas,iw) =
     .    matmul(dble(x0inv(1:nmbas,1:nmbas,1,iw)),mmnorm(1:nmbas))
        do  i = 1, nmbas
          if (abs(momsite(i)) > 1d-3) then
            uu(i,iw) = freq(iw)*momsite(i)/mmnorm2(i)
     .               - uub(i,iw)/mmnorm(i)
          endif
        enddo
      write(stdo,"('UUet :',i4, f13.5, d15.7)")
     &  iw,rymev*freq(iw),uu(1,iw)
      enddo

C#ifdef DEBUG
C     open(112,file='etUet.allw')
C     ifi = fopnx('etUet.allw',2,2,-1)
      ifi = fopna('etUet',-1,0)
      do  iw = 1, nwx
        write(ifi,'(f21.13,3x,255d18.10)')
     .    freq(iw),(uu(ix,iw),ix=1,nmbas)
      enddo
      call fclose(ifi)
C#endif

C --- x2et=<e~|X(w)|e~>  &&  <eiqr|e~><e~|X|e~><eiqr|e~> ---
      call info0(20,1,0,' STONERRSA:  calculate full Xi for each q'//
     .  '%N Magnetic moments, by (magnetic) site:')
      if (ipr >= 20) write(stdo,"(1x,12f8.4)") momsite

C --- Determine ferromagnetic or anti-ferromagnetic ---
      iferro=0  ! duno
      if (nmbas == 1) iferro=1  !ferromagnetic
      if (nmbas == 2) then
        if (abs(sum(momsite)) < 1d-2) iferro=2  !anti-ferro
        if ( abs(momsite(1)-momsite(2)) < 1d-2)  iferro=1
      endif
      if (nmbas >= 3) then
         if (abs(sum(momsite)) < 1d-2)  iferro=2
         if ((sum(abs(momsite))-abs(sum(momsite))) < 1d-2)  iferro=1
      endif
      write(stdo,"(1x,'iferro:0(DUNO),1(FM),2(AF)',2x,'=',i3)") iferro
      write(stdo,*)
C     End of block to determine FM/AFM

C --- Full chiinv(q,omega) = chi0inv + U ---
C     Make qxq = <eiqr| X |eiqr> = <eiqr|e~> <e~|X|e~> <e~|eiqr>
C     where <e~|X|e~> is an nmbas x nmbas matrix
C     SWs occur at poles of qxq.

C! X^(-1)=X_0^{-1}+U      X=( I + X_{0}*U )^{-1} * X_{0}
C!                        X=X_{0} * ( I + U*X_{0} )^{-1}

c      x2et=0.d0
c      do  iq = 1, nq
c      do  iw = 1, nwx
c         cxtmp=0.d0
c         uumtx=0.d0
c         do  ix = 1, nmbas
c            if (iq == 1) then
c               uumtx(ix,ix) = uu(ix,iw)+  img*1d-30
c            else
c               uumtx(ix,ix) = uu(ix,iw)
c            endif
c         enddo
c         call matm(x0et(:,:,iq,iw), uumtx(:,:),cxtmp(:,: ),
c     .              nmbas,nmbas,nmbas  )
cc      call ZGEMM ( "N", "N", nmbas, nmbas, nmbas, 1d0,
cc     &             x0et(1,1,iq,iw), nmbas,
cc     &             uumtx, nmbas,
cc     &             0d0, cxtmp, nmbas )
c
c        do  ix = 1, nmbas
c          cxtmp(ix,ix) = cxtmp(ix,ix) + 1.d0
c        enddo
c        call matcinv(nmbas, cxtmp)
c        call matm(cxtmp(:,: ), x0et(:,:,iq,iw), x2et(:,:,iq,iw),
c     .              nmbas,nmbas,nmbas  )
cC! Full x_+-  ! x2et=<e~|X(w)|e~>
c      enddo
c      enddo


      x0inv = x0et
      do  iq = 1, nq
      do  iw = 1, nwx
        call matcinv(nmbas,x0inv(:,:,iq,iw))
        xinv(:,:,iq,iw) = x0inv(:,:,iq,iw)
        do  ix = 1, nmbas
          xinv(ix,ix,iq,iw) = x0inv(ix,ix,iq,iw) + uu(ix,iw)
        enddo
        cxtmp(:,:) = xinv(:,:,iq,iw)
        do  ix = 1, nmbas
c          cxtmp(ix,ix) = cxtmp(ix,ix) + img*1d-30
          cxtmp(ix,ix) = cxtmp(ix,ix) + img*1d-3
        enddo
C       Liqin: this is poor programming technique
        call matcinv(nmbas,cxtmp(:,:)) ! Full x_+-
        x2et(:,:,iq,iw) = cxtmp        ! x2et=<e~|X(w)|e~>
      enddo
      enddo

C     <eiqr|e~> <e~|X|e~> <e~|eiqr>
      do  iq = 1, nq
      do  iw = 1, nwx
        qxq(iw,iq) = sum(eiqrm(:,iq)
     .    *matmul(x2et(:,:,iq,iw),dconjg(eiqrm(:,iq) )))
        qxq_r(iw,iq) = dble(qxq(iw,iq) )
        qxq_i(iw,iq) = dimag(qxq(iw,iq) )

        qxq_inv_r(iw,iq) = dble( 1d0/qxq(iw,iq) ) !for interpolation
        qxq_inv_i(iw,iq) = dimag( 1d0/qxq(iw,iq) )

c      write(stdo,"( i4, f13.5, 2d15.7 )") iw,rymev*freq(iw)
c     .             , qxq_r(iw,iq),qxq_i(iw,iq)
      enddo
      enddo


C     open(111, file='qm.allq',form='unformatted')
c      ifi = fopnx('u_allq',2,2+4,-1)
c      write (ifi) nq,nwx, 1      ! 1 = dimension of qx0q for each q,w
c      write (ifi) qp
c      write (ifi) freq(1:nwx)
c      do  iw = 1, nwx
c        write(ifi) dble( uu(1:nmbas,iw))
c        write(ifi) dimag(uu(1:nmbas,iw))
c      enddo
c      call fclose(ifi)
! C#endif

      ifi  = fopnx('qxq_allq',2,2+4,-1)
      ifi2 = fopna('qxq_allq_ascii', -1, 2)
      write (ifi) nq,nwx, 1      ! 1 = dimension of qx0q for each q,w
      write (ifi) qp
      write (ifi) freq(1:nwx)
      do  iq = 1, nq
        do  iw = 1, nwx
          write(ifi) dble( qxq(iw,iq))
          write(ifi) dimag(qxq(iw,iq))

          write(ifi2,'(I5, 3F15.7, I7, E20.10, 2E20.10)')
     &          iq, qp(:,iq), iw, freq(iw), qxq(iw,iq)

        enddo
      enddo
      call fclose(ifi)
      call fclose(ifi2)

c      call info0(0,1,0,'%8pmatl invoke: xread(''qxq_allq'') ')

C --- Interpolate <eiqr|X|eiqr> ---
C     from original mesh to fine mesh to find pole in SW spectrum
C     Liqin: these should be passed as arguments
C      emin_intp = 0d0; dw_intp = 1d-1; emax_intp = 1000d0
      emin_intp = 0d0; dw_intp = 2d-1; emax_intp = 1000d0
      nw_intp = int((emax_intp-emin_intp)/dw_intp) + 1
      call info(20,0,0,' Make <q|X|q>: for %i energy points; '//
     .  'emax = %;5,5d eV',nwx,rydberg*freq(nwx))
      call info5(20,0,0,' Interpolate energy window with '//
     .  'emin:dw:emax = %;4d:%;4d:%;4d meV (%i points)',emin_intp,
     .  dw_intp,emax_intp,nw_intp,0)

      allocate(qxq_intp(nw_intp,nq),e_intp(nw_intp) )
      do  iq = 1, nq
      do  iw = 1, nw_intp
        if (iq == 1) then
          e_intp(iw) = (emin_intp+(iw-1)*dw_intp )/(rymev)
        endif
        rtmp1 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_r(1:nwx,iq),nwx)
        rtmp2 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_i(1:nwx,iq),nwx)
        qxq_intp(iw,iq) = 1d0/(rtmp1 + img*rtmp2) !for interpolation
      enddo
      enddo

      if (ipr >= 30) then
      call info0(30,1,0,' Data for pole search of <eiqr| X |eiqr>)%N'//
     .  '%5fq %13fqxq_r_max %14fqxq_r_min %14fqxq_i_max %14fqxq_i_min')
      do  iq = 1, nq
        jmax = maxloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
        jmin = minloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
        imax = maxloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
        imin = minloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
        write(stdo,101) iq,
     .    rymev*e_intp(jmax), -dble(qxq_intp(jmax,iq)),
     .    rymev*e_intp(jmin), -dble(qxq_intp(jmin,iq)),
     .    rymev*e_intp(imax), -dimag(qxq_intp(imax,iq)),
     .    rymev*e_intp(imin), -dimag(qxq_intp(imin,iq))
  101   format(3x,i3,1x, f7.2,'(',d15.7,')',1x
     .    ,f7.2,'(' ,d15.7,')',1x,f7.2,'(' ,d15.7,')'
     .    ,1x,f7.2,'(' ,d15.7,')' )
      enddo
      endif

cC#ifdef DEBUG
cC      open(106,file='qxqi.allq')
cC      open(107,file='qxqr.allq')
c      ifi  = fopna('qxqi',-1,0)
c      ifi2 = fopna('qxqr',-1,0)
c      do  iw = 1, nw_intp
c        write(ifi,"( f13.5, 100d15.7)") rymev*e_intp(iw),
c     .    (-dimag(qxq_intp(iw,iq)),iq=1,nq)
c        write(ifi2,"( f13.5, 100d15.7)") rymev*e_intp(iw),
c     .    (-dble(qxq_intp(iw,iq)),iq=1,nq)
c      enddo
c      call fclose(ifi); call fclose(ifi2)
cC#endif

C      qxq_intp(1,1) = 0.d0

      ifi  = fopnx('qxq_intp_allq',2,2+4,-1)
      ifi2 = fopna('qxq_intp_allq_ascii', -1, 2)

      write (ifi) nq,nw_intp, 1      ! 1 = dimension of qx0q for each q,w
      write (ifi2,*) "# nq, nw_intp = ", nq,nw_intp      ! 1 = dimension of qx0q for each q,w

      write (ifi) qp
      write (ifi) e_intp(1:nw_intp)
      do  iq = 1, nq
        do  iw = 1, nw_intp
          write(ifi) dble( qxq_intp(iw,iq))
          write(ifi) dimag(qxq_intp(iw,iq))

          write(ifi2,'(I5, 3F15.7, I7, E20.10, 2E20.10)') iq, qp(:,iq), iw, e_intp(iw), qxq_intp(iw,iq)

        enddo
      enddo
      call fclose(ifi)
      call fclose(ifi2)

c      call info0(0,1,0,'%8pmatl invoke: xread(''qxq_allq'') ')

C ... Finished interpolation of <eiqr|X|eiqr>

C --- Eigenvalue of hermitian xinvh(1:nmbas,1:nmbas) ---
C     Pole corresponds to maximum eval of hermitian part of chi
      do  iq = 1, nq
      do  iw = 1, nwx
        xinvh(:,:,iq,iw) = .5d0*( xinv(:,:,iq,iw)
     .    + transpose(dconjg(xinv(:,:,iq,iw))) )
        call zevl(nmbas,xinvh(1,1,iq,iw),eval)
        mxevl_xinvh(iq,iw) = maxval(eval)
      enddo
      enddo

C#ifdef DEBUG
      ifi = fopnx('evl_xh.allq',2,2,-1)
      do  iw = 1, nwx
        write(ifi,"( f13.5, 100d15.7)") rymev*freq(iw),
     .    (mxevl_xinvh(iq,iw),iq=1,nq)
      enddo
      call fclose(ifi)
C#endif
C ... Finished finding eigenvalues

C --- Pole search of xinvh ---
C Liqin: there doesn't seem to be a check whether FM or AFM
C answer:  updated  iferro=0(DUNO),1(FM),2(AF)'
      epole = 0
C ... Ferromagnetic case
      do  iq = 1, nq
        iwpolf = 1
C       Coarse search: bracket frequency where evl crosses zero
        do  iw = 1, nwx
          if (freq(iw) >= 0d0 ) then
            if (mxevl_xinvh(iq,iw) < 0d0 .and.
     .          mxevl_xinvh(iq,iw+1) > 0d0 ) then
              iwpolf = iw
              epole = freq(iw)
              exit
            endif
          endif
        enddo

        if (iq == 1) then
          rtmp = mxevl_xinvh(iq,iwpolf)
        elseif (iwpolf /= 1) then ! fine search
          do  ! Liqin ... infinite loop is dangerous
            epole = epole + 1d-7/rydberg
            rtmp  = polinta(epole,freq(iwpolf-1:iwpolf+2),
     .        mxevl_xinvh(iq,iwpolf-1:iwpolf+2),4)
            if (rtmp > 0) exit
          enddo
        endif
        epole_fm(iq) = epole
        vpole_fm(iq) = rtmp
      enddo

C ... Antiferromagnetic case
      do  iq = 1, nq
        iwpola = 1
        epole = freq(iwpola)
        vpole_af(iq) = mxevl_xinvh(iq,iwpola)
C       Coarse search on original mesh supplied by GWinput:
C       bracket frequency where evl crosses zero
        do  iw = 1, nwx-1
          if (freq(iw) >= 0d0 ) then
            if (mxevl_xinvh(iq,iw) > 0d0 .and.
     .          mxevl_xinvh(iq,iw+1) < 0d0 ) then
              iwpola = iw
              epole = freq(iw)
              exit
            endif
          endif
        enddo
C       Refine search where possible
        if (iq == 1) then
          rtmp = mxevl_xinvh(iq,iwpola)
        elseif (iwpola /= 1) then ! fine search
          do
            epole = epole + 1d-7/rydberg
            rtmp  = polinta(epole,freq(iwpola-1:iwpola+2),
     .        mxevl_xinvh(iq,iwpola-1:iwpola+2),4)
            if (rtmp < 0) exit
          enddo
        elseif (iwpola == 1) then
          epole = freq(iwpola)
          rtmp = mxevl_xinvh(iq,iwpola)
        endif
        epole_af(iq) = -epole
        vpole_af(iq) = rtmp
      enddo

      call info0(20,1,0,' Smallest evals of Hermitian part '//
     .  'of <q|X^-1(1:nmbas,1:nmbas)|q> (meV) %N'//
     .  '%5fq %9fFM pole %23fAFM pole')
      do  iq = 1, nq
         write(stdo,102) iq, rymev*epole_fm(iq), vpole_fm(iq)
     .    ,rymev*epole_af(iq), vpole_af(iq)
  102    format(3x,i3,3x,f8.2,3x, d15.7, 5x,f8.2,3x,d15.7)
      enddo
C ... End of pole search

C --- Determine effective moment meffi ---
C     meffi = [xinv(om=pole) - xinv(om=0)]/om
C     constructed so that Heisenberg J has pole at proper point
      do  iq = 1, nq
        rtmp1 = epole_fm(iq)
        rtmp2 = epole_af(iq)
      do  i = 1, nmbas
      do  j = 1, nmbas
        do  iw = 1, nwx-1
          dxidw(i,j,iq,iw)  =
c     .   (xinv(i,j,iq,iw+1) - xinv(i,j,iq,1))/(freq(iw+1)-freq(1))
     .      (xinvh(i,j,iq,iw+1) - xinvh(i,j,iq,1))/(freq(iw+1)-freq(1))
          freqm(iw) = 0.5d0*(freq(iw+1) + freq(iw))
        enddo
C       Interpolation near FM pole
        meffi(i,j,iq) =
     .      polinta(rtmp1,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8)
     .    + img* polinta(rtmp1,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
C       Interpolation near AFM pole
        meffi2(i,j,iq) =
     .      polinta(rtmp2,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8)
     .    + img* polinta(rtmp2,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
      enddo
      enddo
      enddo

C --- Xinvh(om=0) in ebar basis to establish corresondence with J ---
      do  iq = 1, nq
        xinvh_w0(:,:,iq) = xinvh(:,:,iq,1)
        x0inv_w0eb(:,:,iq) = x0inv(:,:,iq,1)
        do  i = 1, nmbas ! projected on |e-> instead of |e~>
           if (iq == 1) then
           uu0_eb(i) = uu(i,1)*mmnorm(i)/momsite(i)*mmnorm(i)/momsite(i)
           endif
          do  j = 1, nmbas
            x0inv_w0eb(i,j,iq) = x0inv_w0eb(i,j,iq)
     .        * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            xinvh_w0eb(i,j,iq) = xinvh_w0(i,j,iq)
     .        * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            meffi_eb(i,j,iq) = meffi(i,j,iq)
     .        * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            meffi2_eb(i,j,iq) = meffi2(i,j,iq)
     .        * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            do  iw = 1, nwx-1
              dxidw_eb(i,j,iq,iw) = dxidw(i,j,iq,iw)
     .          * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            enddo
          enddo
        enddo
      enddo

      call info0(20,1,0,' <eb|I_stoner|eb>_{q=0,w=0}')
      write(stdo,"(t33,12f10.4)")  uu0_eb

      call info0(20,1,0,' Inverse of effective mag. moment meffi(1,1)'//
     .  ' = [xinvh(1,1,iw)-xinvh(1,1,iw=0)]/w%N'//
     .  '%5fq%11fiw=1%21fFM pole%18fAFM pole%17fintp FM%18fintp AF')
      do  iq = 1, nq
         write(stdo,  "(3x,i3, 5(2d12.4,1x) )" ) iq,dxidw_eb(1,1,iq,1)
     .        ,dxidw_eb(1,1,iq,iwpolf),dxidw_eb(1,1,iq,iwpola)
     .        , meffi_eb(1,1,iq), meffi2_eb(1,1,iq)
      enddo

C#ifdef DEBUG
C     open(109,file='dxidw.allq')
C     ifi = fopna('dxidw',-1,0)
      ifi = fopnx('dxidw.allq',2,2,-1)
      do  iw = 1, nwx-1
        write(ifi,"(f8.2,100d15.7)") rymev*freqm(iw),
     .    (dxidw_eb(1,1,iq,iw), iq=1,nq)
      enddo
      call fclose(ifi)
C#endif

C --- Make Jmat file for FM and AFM cases ---
      call info0(20,1,0,'Writing files to disk: '//
     .  'Jmat.allq = J(q,w=intp FM pole), and  '//
     .  'Jmat_X2w0 = J(q,w=0)')
      nmx = nmbas
      nev = nmbas
C     open(110,file='Jmat.allq')
      ifi = fopnx('Jmat.allq',2,2,-1)
      do  iq = 1, nq
C       oo = 0d0
        mxi = 0d0
        do  i = 1, nmbas
C         oo(i,i) = 1d0
          mxi(i,i) = 1d0
        enddo
        allocate(wk(11,nmbas))
        call zhev(nmbas,meffi_eb(1,1,iq),wk,
     .    .false.,.true.,nmx,1d99,nev,wk,.false.,-1,eval,evc)
        deallocate(wk)
C       call diagcv(oo,meffi_eb(1,1,iq),evc, nmbas, eval,nmx,1d99,nev)

        oo = 0d0
        do  i = 1, nmbas
          if (eval(i) >= 0) then
            oo(i,i) = 1d0/sqrt(eval(i))
          else
            oo(i,i) = 1d0/csqrt( cmplx(eval(i)) )
          endif
C         oo(i,i)= 1d0/sqrt(eval(i))
        enddo
C       write(stdo,"( 'oo', 100f13.7)") (oo(i,i),i=1,nmbas)
        sqm = matmul(evc, matmul(oo, transpose(dconjg(evc))) )
        jjmat = matmul(sqm, matmul(xinvh_w0eb(:,:,iq),sqm))
        jjmat2 = matmul(sqm, matmul(mxi, sqm))
c       jjmat = xinvh_w0eb(:,:,iq)
c     bug--momsite could be negative
        do  ix = 1, nmbas
          do  iy = 1, nmbas
            jjmat(ix,iy) = jjmat(ix,iy)
     .            /csqrt(cmplx(momsite(ix)*momsite(iy)))
          enddo
        enddo

C        jjmat2 = xinvh_w0eb(:,:,iq)
C        do  ix = 1, nmbas
C          do  iy = 1, nmbas
C            jjmat2(ix,iy)= jjmat2(ix,iy)/sqrt(momsite(ix)*momsite(iy))
C          enddo
C        enddo

C         open(103,file='Jmat_X3w0.allq')
C         write(103,"('sqmIsqm: ',3d18.10, 3x, 255d18.10)")
C     .     qp(:,iq), (( jjmat2(ix,iy) ,ix=1,nmbas),iy=1,nmbas)


        write(ifi,"('JJMAT: ',3d18.10,3x,255d18.10)")
     .    qp(:,iq), ((jjmat(ix,iy), ix=1,nmbas),iy=1,nmbas)
        call zevl(nmbas,jjmat,eval)

C         write(stdo,"('e sw  ', 3f8.4,2x,255d15.7)")
C     .     qp(:,iq),(-momsite(1)*1d3*rydberg*eval(ix),ix=1,nmbas)
      enddo

C --- Save files: JJMAT: J=X0^(-1)  or J=X^(-1) at (w=0) ---
C     open(104,file='Jmat_X0w0.allq')
C     open(105,file='Jmat_X2w0.allq')
      ifi  = fopnx('Jmat_X0inv_w0eb.allq',2,2,-1)
      ifi2 = fopnx('Jmat_Xinv_w0eb.allq',2,2,-1)
      do  iq = 1, nq
        write(ifi,"(3d18.10, 3x, 255d18.10)")
     .    qp(:,iq), (( x0inv_w0eb(ix,iy,iq) ,ix=1,nmbas),iy=1,nmbas)   ! lmgf did change the sign
c$$$     .    qp(:,iq), (( -x0inv_w0eb(ix,iy,iq) ,ix=1,nmbas),iy=1,nmbas)
c        write(ifi2,"('xinvh_w0eb: ',3d18.10, 3x, 255d18.10)")
c     .    qp(:,iq), ((  xinvh_w0eb(ix,iy,iq) ,ix=1,nmbas),iy=1,nmbas)
        write(ifi2,"( 3d18.10, 3x, 255d18.10)")
c     .    qp(:,iq), ((  xinvh_w0eb(ix,iy,iq) ,ix=1,nmbas),iy=1,nmbas)
     .    qp(:,iq), ((  -xinvh_w0eb(ix,iy,iq) ,ix=1,nmbas),iy=1,nmbas)
      enddo

      deallocate(x0inv)
      deallocate(xinv,x2et)
      deallocate(xinvh,mxevl_xinvh)
      deallocate(dxidw)
      deallocate(dxidw_eb)
      deallocate(mxevl2_xinvh)

      end subroutine stonerrsa

      subroutine zevl(n,h,eval)
C- Return eigenvalues of a hermitian matrix h, leaving h unaltered
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :dimension of h
Ci   h     :small h = cnu-enu+sqrdel*S^beta*sqrdel
Co Outputs
Co   eval  :eigenvalues
Cl Local variables
Cr Remarks
Cu Updates
Cu   08 Feb 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      double complex h(n,n)
      double precision eval(n)
C ... Local parameters
      integer nev
      double precision xx
      complex(8),allocatable:: hloc(:,:),z(:,:),wk(:,:)

      allocate(wk(n,n),z(n,n),hloc(n,n))
      hloc = h
      call zhevx(n,n,hloc,xx,0,.true.,n,1d99,nev,wk,.false.,eval,n,z)
      deallocate(z,wk,hloc)

      end

      subroutine matcinv(n,a)
C- Inverse of a complex matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   n
Cio Inputs/Outputs
Ci   a     :matrix of dimension n.  On output, the inverse of a is returned
Cl Local variables
Cr Remarks
Cu Updates
Cu   09 Feb 09
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n
      complex(8):: a(n,n)
C ... Local parameters
      integer info,ipiv(n)
      complex(8),allocatable:: work(:)

      call zgetrf(n,n,a,n,ipiv,info)
      if (info /= 0) call rxi
     .  ('zinv failed to invert matrix: info =',info)
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if (info /= 0) call rxi
     .  ('zinv failed to invert matrix: info =',info)
      end

      double precision function polinta(x,xa,ya,n)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polint. If p(x) is the polynominal
c     of degree ndg such that p(xa(i))=ya(i), i=ns,..,ns+ndg then
c     the returned value polint=p(x). ns is obtained by hunting.
c     See Numerical Recipes
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ndgmx=4, nmx=ndgmx+1)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      logical ascnd
      save jlo
      data jlo/0/ , small/1d-30/
      ndg=min(ndgmx,n-1)
      ndt=ndg+1
      ascnd=xa(n) > xa(1)
      if (jlo <= 0 .or. jlo > n) then
      jlo=0
      jhi=n+1
      go to 30
      endif
      inc=1
      if (x > xa(jlo) .eqv. ascnd) then
   10 jhi=jlo+inc
      if (jhi > n) then
      jhi=n+1
      else if (x. gt. xa(jhi) .eqv. ascnd) then
      jlo=jhi
      inc=inc+inc
      go to 10
      endif
      else
      jhi=jlo
   20 jlo=jhi-inc
      if (jlo < 1) then
      jlo=0
      else if (x < xa(jlo) .eqv. ascnd) then
      jhi=jlo
      inc=inc+inc
      go to 20
      endif
      endif
   30 if (jhi-jlo /= 1) then
      jm=(jhi+jlo)/2
      if (x > xa(jm) .eqv. ascnd) then
      jlo=jm
      else
      jhi=jm
      endif
      go to 30
      endif
      nlo=max(1,jlo-ndg/2)
      nhi=min(n,nlo+ndg)
      nlo=nhi-ndg
      if (jlo == 0) then
      ns=1
      else if (jlo == n) then
      ns=ndt
      else if (abs(x-xa(jlo)) < abs(x-xa(jhi))) then
      ns=jlo-nlo+1
      else
      ns=jhi-nlo+1
      endif
      do  i = 1, ndt
      ii=nlo+i-1
      c(i)=ya(ii)
      d(i)=ya(ii)
      enddo
      polint=ya(nlo+ns-1)
      ns=ns-1
      do  m = 1, ndg
        do  i = 1, ndt-m
        ii=nlo+i-1
        ho=xa(ii)-x
        hp=xa(ii+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
c
c       an error can occur if two xa's are identical
        if (abs(den) < small) then
        write(6,1000)
 1000   format('   ***wrn in polint...data error')
        stop
        endif
c
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
        enddo
      if (2*ns < ndt-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
      polint=polint+dy
c takao
      enddo
      polinta=polint
      return
      end
