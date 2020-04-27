      subroutine chcase(iopt,len,strn)
C- Changes the case of a string
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt  :=  -1 convert from uppercase to lowercase
Ci         :>   1 convert from lowercase to uppercase
Ci   len   :length of string
Co Outputs:
Co   strn  :string converted to specified case
Cu Updates
Cu   02 Nov 01 Adapted from Stuttgart lmto56
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iopt,len
      character strn(*)
C Local variables:
      integer i,nuplo
C Intrinsic functions:
      intrinsic  ichar,char

      nuplo = ichar('A') - ichar('a')
      do  i = 1, len
        if (iopt == -1) then
          if (strn(i) >= 'A' .and. strn(i) <= 'Z') then
            strn(i) = char(ichar(strn(i))-nuplo)
          endif
        else
          if (strn(i) >= 'a' .and. strn(i) <= 'z') then
            strn(i) = char(ichar(strn(i))+nuplo)
          endif
        endif
      enddo
      end
