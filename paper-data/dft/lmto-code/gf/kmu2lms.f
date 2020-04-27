      subroutine kmu2lms1(opt,l,imu,nl,reskmu,reslms)
C- Rotates one (imu,l) block to lms representation
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :1 rotate kmu -> lms
Ci   l     :
Ci   imu   :
Ci   nl    :dimensions reslms
Cio Inputs/Outputs
Ci   reskmu:
Ci   reslms:
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   m are ordered l ... -l, and ms are ordered 1/2, -1/2
Cu Updates
Cu   15 May 18
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,l,imu,nl
      double complex reskmu(2,2),reslms(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2)
C ... Local parameters
      integer m1,m2,ms1,ms2,opt0,opt1,im1,im2
      double precision mu,u,clebsh(2,2)

C ... Setup for rotation of kappa-mu to lms rep'sn
      opt0 = mod(opt,10)
      if (opt0 == 0) return
      opt1 = mod(opt/10,10)

      mu = imu - l - 1.5d0
      u = mu/(l+0.5d0)
      clebsh(1,1) =  dsqrt((1+u)/2) ; clebsh(1,2) = -dsqrt((1-u)/2)
      clebsh(2,2) =  clebsh(1,1) ; clebsh(2,1) = - clebsh(1,2)

      if (opt0 == 1) then
        im1 = 0; im2 = 0
        call dpzero(reslms,2*size(reslms))
        do  ms1 = 1, 2
          do  ms2 = 1, 2
            m1 = int(mu - (ms1-1.5d0))
            m2 = int(mu - (ms2-1.5d0))
            im1 = 0; im2 = 0
            if (iabs(m1) <= l .and. iabs(m2) <= l) then
              reslms(m1,m2,ms1,ms2) =
     .          clebsh(1,ms1)*reskmu(1,1)*clebsh(1,ms2) +
     .          clebsh(1,ms1)*reskmu(1,2)*clebsh(2,ms2) +
     .          clebsh(2,ms1)*reskmu(2,1)*clebsh(1,ms2) +
     .          clebsh(2,ms1)*reskmu(2,2)*clebsh(2,ms2)
            endif
          enddo
        enddo
      else
        call rxi('kmu2lms not ready for opt=',opt)
      endif

      end
