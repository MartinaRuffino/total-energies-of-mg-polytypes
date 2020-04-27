      subroutine mksigpade(s_pade,mode,nomg,ndsig,zomg,siginp,sigpade)
C- Interpolate Self-energy by Pade approximant
C ----------------------------------------------------------------------
Cio Structures
Cio  s_pade
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lpade npade zomega padcof
Cio    Passed to:  *
Ci Inputs
Ci   mode  :0 copy sigma from reference sigma
Ci         :1 If Pade has NOT been set up (s_pade%lpade): reverts to mode 0.  Otherwise:
Ci         :  Return sigpade = siginp + Pade[siginp(zomg)] - Pade[siginp(Im(zomg))]
Ci         :  Only correct if siginp evaluated at Im(zomg)
Ci         :  Idea is to exploit cancellation in errors of Pade approximants
Ci         :  It may be that Pade may not exactly match siginp at Im(zomg)
Ci         :2 Return Pade approximation for sigma
Ci   nomg  :leading dimension of siginp
Ci   ndsig :number of functions to interpolate
Ci   zomg  :Make sigma at frequency zomg
Ci   siginp:sigma on Matsubara axis, at frequency Im(zomg)
Co Outputs
Co  sigpade:either siginp (mode=0, possibly mode 1) or Pade approximation to siginp
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   14 Apr 17  Implemented mode 2
Cu   13 Feb 16  (MvS) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nomg,ndsig
      double precision zomg(2),siginp(nomg,2,ndsig),sigpade(2,ndsig)
C ... For structures
!      include 'structures.h'
      type (str_mpade):: s_pade
C ... Local parameters
      integer i
      double precision z(2)
      double precision xx(2,2),df(2)

      if (mode == 0 .or. mode == 1) then
        call dcopy(2*ndsig,siginp,nomg,sigpade,1)
        if (mode == 0) return
      endif
      if (mode == 2 .and. .not. s_pade%lpade)
     .  call rx('Pade not set up; cannot interpolate')
      if (.not. s_pade%lpade) return

      do  i = 1, ndsig

        call pade(1,zomg,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,1))
        if (mode == 2) then
          sigpade(:,i) = xx(:,1)
          cycle
        endif
        z(1) = 0; z(2) = zomg(2)
        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,2))
        df(1:2) = xx(:,1) - xx(:,2)
        sigpade(:,i) = sigpade(:,i) + df(:)

      enddo

      end
      subroutine supade(s_pade,mode,nomg,icut,ndsig,omegai,siginp)
C- Setup for Pade approximation to Sigma
C ----------------------------------------------------------------------
Cio Structures
Cio  s_pade
Ci     Elts read:  padcof zomega
Co     Stored:     npade
Co     Allocated:  padcof zomega*
Cio    Elts passed:padcof zomega
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 include points for positive Imaginary frequencies only
Ci         :1 include points for positive and negative Imaginary frequencies
Ci   nomg  :Number of Matsubara frequencies available; dimensions siginp
Ci   icut  :icut(1), if nonzero:
Ci         :for frequencies above icut, double spacing between successive points
Ci         :icut(2), if nonzero:
Ci         :Exclude frequencies above with index > icut(2)
Ci   ndsig :Number of independent matrix elements in the self-energy
Ci   omegai:Im(omega) (usually on the Matsubara frequencies)
Ci   siginp:Self energy for each ndsig and frequency in omegai
Co Outputs
Co  padcoff:Pade coefficients (mode=1)
Cs Command-line switches
Cl Local variables
Cl   deltai: spacing between matsubara frequencies
Cr Remarks
Cr
Cu Updates
Cu   12 Feb 16 (MvS) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nomg,ndsig,icut(2)
      double precision omegai(nomg),siginp(nomg,2,ndsig)
C ... For structures
!      include 'structures.h'
      type (str_mpade):: s_pade
C ... Dynamically allocated local arrays
      complex(8), pointer :: zomega(:),padcoff(:,:,:)
      integer, allocatable :: iprm(:)
      real(8), allocatable :: wk(:,:),omega(:,:)
C ... Local parameters
      integer i,k,iw,npade,deltai,mode0,nptsiw
      double precision rms1,rms2  !,z(2),xx(2,2)
      real(8), parameter :: ry2ev = 13.60569193d0
      procedure(integer) :: iprint
C      integer icut,ncuti(4)
C      double precision cuti(4)
C      data ncuti /1,5,10,30/
C      data cuti /5d0,20d0,50d0,9999d0/

      mode0 = mod(mode,10)


C ... Make a list of points on which to make Pade interpolation
C     icut = 1
      iw = 0; deltai=1
      npade = 0                 ! Loop determines npade = # points for Pade
      allocate(iprm(nomg))      ! Maps interpolating mesh points to original mesh
      do while (iw+deltai <= nomg)
        iw = iw + deltai
        if (icut(2) > 0 .and. iw > icut(2)) exit
        npade = npade+1
        iprm(npade) = iw
        if (icut(1) > 0 .and. iw >= icut(1)) then
          deltai = 2*deltai     !double the increment between points
        endif
C        if (omegai(iw)*ry2ev > cuti(icut)) then
C          icut = icut+1
C        endif
      enddo

      nptsiw = npade
      if (mode0 == 2) npade = npade*2           ! To include both sigma(iw) and sigma(-iw)
      s_pade%npade = npade
      allocate(s_pade%padcof(npade,npade+2,ndsig))
      allocate(s_pade%zomega(npade))
      padcoff => s_pade%padcof
      zomega  => s_pade%zomega
      do  iw = 1, nptsiw
        select case (mode0)
        case (0)
          zomega(iw) = dcmplx(0d0,omegai(iprm(iw)))
        case (1)
          zomega(iw) = dcmplx(0d0,-omegai(iprm(iw)))
        case (2)
          zomega(2*iw-1) = dcmplx(0d0,omegai(iprm(iw)))
          zomega(2*iw)   = dcmplx(0d0,-omegai(iprm(iw)))
        case default
          call rxi('supade : invalid mode',mode)
        end select
      enddo
      allocate(wk(2,npade))
      do  i = 1, ndsig
        do  iw = 1, nptsiw
        select case (mode0)
        case (0)
          wk(1:2,iw)   = siginp(iprm(iw),1:2,i)
        case (1)
          wk(1,iw)     = siginp(iprm(iw),1,i)
          wk(2,iw)     = -siginp(iprm(iw),2,i)
        case (2)
          wk(1:2,2*iw-1) = siginp(iprm(iw),1:2,i)
          wk(1,2*iw)     = siginp(iprm(iw),1,i)
          wk(2,2*iw)     = -siginp(iprm(iw),2,i)
        end select
        enddo
        call padcof(npade,zomega,wk,npade,padcoff(1,1,i))
      enddo

      rms1 = 0
      do  i = 1, ndsig
        call pade(npade,zomega,npade,zomega,npade,padcoff(1,1,i),wk)
        do  iw = 1, nptsiw
          k = iw; if (mode0 == 2) k = 2*iw-1
C          print *, dimag(zomega(k)),
C     .      wk(1,k)-siginp(iprm(iw),1,i),wk(2,k)-siginp(iprm(iw),2,i)
          rms1 = rms1 + (wk(1,k)-siginp(iprm(iw),1,i))**2
     .                + (wk(2,k)-siginp(iprm(iw),2,i))**2
        enddo
      enddo
      rms1 = sqrt(rms1/(npade*ndsig))
      deallocate(wk)

      rms2 = 0
      allocate(omega(2,nomg),wk(2,nomg))
      do  iw = 1, nomg
        omega(1,iw) = 0
        omega(2,iw) = omegai(iw)
      enddo
      do  i = 1, ndsig
        call pade(nomg,omega,npade,zomega,npade,padcoff(1,1,i),wk)
        do  iw = 1, nomg
C         print *, omega(2,iw),wk(1,iw)-siginp(iw,1,i),wk(2,iw)-siginp(iw,2,i)
          rms2 = rms2 + (wk(1,iw)-siginp(iw,1,i))**2
     .                + (wk(2,iw)-siginp(iw,2,i))**2
        enddo
      enddo
      s_pade%lpade = .true.
      rms2 = sqrt(rms2/(nomg*ndsig))
      call info5(30,1,0,
     .  ' Pade fit to sigma on %i points (icut=%i,%i): RMS err %;3g (fit pts) %;3g (all pts)',
     .  npade,icut(1),icut(2),rms1,rms2)
      if (iprint() >= 41-10) then
        call info0(10,0,0,'  point       omega')
        do  iw = 1, npade
          call info2(10,0,0,'%,5i %2;10,6D',iprm(iw),zomega(iw))
        enddo
      endif

C      do  i = 1, ndsig
C
C      do  iw = 1, nomg
C        z(1) = 0.00d0; z(2) = omegai(iw)
C        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,1))
C        z(1) = 0.001d0; z(2) = omegai(iw)
C        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,1))
C        z(1) =-0.001d0; z(2) = omegai(iw)
C        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,2))
C        print *, iw,i,j, sngl(xx(:,1)/.002d0 - xx(:,2)/.002d0)
C
CC        z(1) = 0.00d0; z(2) = omegai(iw)+.001d0
CC        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,1))
CC        z(1) = 0.00d0; z(2) = omegai(iw)-.001d0
CC        call pade(1,z,s_pade%npade,s_pade%zomega,s_pade%npade,s_pade%padcof(1,1,i),xx(1,2))
CC        print *, iw,i,j, xx(:,1)/.002d0 - xx(:,2)/.002d0
C
C      enddo
C        stop
C
C      enddo

      deallocate(omega,iprm,wk)

      end
