      subroutine modsig(mode,nomg,nband,nqf,nspse,ommin,ommax,sigma)
C- Modifications of/check for causality in self-energy sigma
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :1 check for causality in Im Sigma
Ci         :2 Abort if causality not satisfied
Ci   nomg  :
Ci   nband :number of bands
Ci   nqf   :
Ci   nspse :
Ci   ommin :
Ci   ommax :
Cio Inputs/Outputs
Ci   sigma :
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Aug 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nomg,nband,nqf,nspse
      real(8) :: ommin,ommax
      complex(8) :: sigma(nomg,nband,nqf,nspse)
C ... Local parameters
      integer iw,ib,iq,is,nerr
      real(8) :: omg,dw0,fac
      real(8), parameter:: tolw=1d-7

      dw0 = 0; if (nomg > 1) dw0 = (ommax-ommin)/max(nomg-1,1)

C ... Check causality
      if (mod(mode,10) == 0) return
      call info0(20,0,-1,' checking causality in sigma ...')
      nerr = 0
      do  iw = 1, max(nomg-1,1)
        omg = ommin + dw0*(iw-1)
        fac = 1
        if (omg > 0) fac = -1

        do  ib = 1, nband
        do  iq = 1, nqf
        do  is = 1, nspse

          if (abs(omg) <= tolw .and. abs(dimag(sigma(iw,ib,iq,is)))>tolw .or.
     .        abs(omg) >  tolw .and. dimag(sigma(iw,ib,iq,is))*fac < 0) then
            if (nerr == 0) call info5(20,1,0,
     .          ' modsig: Im(sigma) wrong sign at omg=%d  iw,ib,iq,is = %4:1i  Sig = %2:1;6F',
     .          omg,[iw,ib,iq,is],sigma(iw,ib,iq,is),4,5)
            nerr = nerr+1
          endif

        enddo
        enddo
        enddo ! ib,iq,is

      enddo ! iw

      if (nerr == 0) call info0(20,0,0,' ok')
      if (nerr > 0) call info2(20,0,0,' encountered %i noncausal values out of %i',nerr,max(nomg-1,1)*nband*nqf*nspse)
      if (nerr > 0 .and. mod(mode,10) > 1) call rx('sigma does not satisfy causality')

      end
