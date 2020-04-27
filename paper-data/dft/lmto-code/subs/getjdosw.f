      subroutine getjdosw(mode,i0,strn,n1,n2,iv1,iv2)
C- Get list of channels for decomposition of DOS
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :0 parse nothing; return n1=n2=0
Ci           :1 parse for occ,unocc DOS strings, return n1,n2
Ci           :2 parse for occ,unocc DOS strings, return iv1,iv2
Ci           :3 same as 2, but:
Ci           :  check is made that n1,n2 match generated value
Ci   strn,i0 :  parse list from strn(i0:*)
Co Outputs
Co   n1    :number of channels in 1st string
Co   n2    :number of channels in 2nd string (0 if none)
Co   iv1   :vector of channels in 1st string
Co   iv2   :vector of channels in 2nd string
Cl Local variables
Cr Remarks
Cr   Strings denote a single list or double list , e.g.
Cr      --jdosw~5,6,8       (single list for channels 5,6,8)
Cr      --jdosw~5,6,8~7,9   (double list for chans 5,6,8 and 7,9)
Cr   The former is for DOS-like objects, the second for joint DOS.
Cu Updates
Cu   28 Nov 16 Use mkilsd to allow dup in strings
Cu   26 Sep 09 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,n1,n2,i0,iv1(n1),iv2(n2)
      character strn*(*)
C ... Local parameters
      character dc*(1)
      integer i,j,k,m
      procedure(integer) :: mkilsd

      dc = strn(i0:i0)
      if (mode < 2) then
        n1 = 0
        n2 = 0
      endif
      if (mode == 0 .or. dc == ' ') return

      i = i0
      call wordg(strn(i0+1:),0,dc//' ',1,i,j)
      if (i <= j) then
        m = mkilsd(strn(i+i0:j+i0),-1,k)
        if (m < 0) call rx('GETJDOSW: failed to parse '//trim(strn))
        if (mode == 1) then
          n1 = m
        else
          if (mode == 3 .and. n1 /= m)
     .      call rx('GETJDOSW: n1 wrongly dimensioned')
          call mkilss(11,strn(i+i0:j+i0),n1,iv1)
        endif
        i = j+2
      else
        i = i+1
      endif
      call nwordg(strn(i0+1:),0,dc//' ',1,i,j)
      if (i > j) return
C     call mkils0(strn(i+i0:j+i0),m,k)
      m = mkilsd(strn(i+i0:j+i0),-1,k)
      if (m < 0) call rx('GETJDOSW: failed to parse '//trim(strn))
        if (mode == 1) then
          n2 = m
        else
          if (mode == 3 .and. n2 /= m)
     .      call rx('GETJDOSW: n2 wrongly dimensioned')
          call mkilss(11,strn(i+i0:j+i0),n2,iv2)
        endif

      end
