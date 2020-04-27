
! The bitwise functions shall be replaced in the code one of these days.
   function bitand(i1, i2)
      implicit none
      integer, intent(in) :: i1, i2
      integer :: bitand
      bitand = iand(i1,i2)
   end function bitand

   function bitor(i1, i2)
      implicit none
      integer, intent(in) :: i1, i2
      integer :: bitor
      bitor = ior(i1,i2)
   end function bitor


   subroutine nlchar(i, s)
      use :: iso_c_binding
      implicit none
      integer, intent(in) :: i
      character, intent(out) :: s

      if (i == 1) s = c_new_line
   end subroutine nlchar

   subroutine locase(s)
      implicit none
      character(len=*), intent(inout) :: s
      integer,parameter :: ua = ichar('A'), uz=ichar('Z'), la = ichar('a'), lz = ichar('z')
      integer,parameter :: d = la - ua
      integer :: i, ic

      do i=1,len(s)
         ic = ichar(s(i:i))
         if ((ua <=ic) .and. (ic <= uz)) s(i:i) = char(ic + d)
      end do
   end subroutine locase



   subroutine flushs(l)
      use :: iso_c_binding
      use posix
      implicit none

      integer, intent(in) :: l
      logical :: opn
      integer :: r

      if (l >= 0) then
         inquire(unit=6,opened=opn)
         if (opn) flush(6)
      else
!         r = fflush(0_c_long)
         r = fflush(c_null_ptr)
      end if
   end subroutine flushs


   function ssizei(p1,p2) result(s)
   ! Returns the distance between two integer pointers
   ! This obviously works only on single, linear address space.
      implicit none
      interface
         function cssizei(p1,p2) bind(c) result(s)
            use iso_c_binding, only : c_int
            integer(c_int), intent(inout), dimension(*) :: p1, p2
            integer(c_int) :: s
         end function cssizei
      end interface
      integer, intent(inout), dimension(*) :: p1,p2
      integer(8) :: s
      s = cssizei(p1,p2)
   end function ssizei


   function ssized(p1,p2) result(s)
! Returns the distance between two double precision pointers
! This obviously works only on single, linear address space.
      implicit none
      interface
         function cssized(p1,p2) bind(c) result(s)
            use iso_c_binding, only : c_int, c_double
            real(c_double), intent(inout), dimension(*) :: p1, p2
            integer(c_int) :: s
         end function cssized
      end interface
      real(8), intent(inout), dimension(*) :: p1,p2
      integer(8) :: s
      s = cssized(p1,p2)
   end function ssized


   subroutine mkdcon(dmach, d1mach)
      implicit none
      integer, parameter :: p = 8
      real(p), intent(out) :: dmach(3), d1mach(5)
      real(p), parameter :: x = 1.0_p
      real(p) :: eps, tin, hug
      integer :: dgt, mxe, mne, rdx

      dgt = digits(x)
      mxe = maxexponent(x)
      mne = minexponent(x)
      rdx = radix(x)

      eps = 2.0_p**(1-dgt)
      tin = 2.0_p**(mne-1)
      hug = (2.0_p**dgt - 1.0_p)*2.0_p**(mxe-dgt)

      dmach = (/eps, tin, hug/)

      d1mach = (/tin, hug, eps/real(rdx,p), eps, log10(real(rdx,p)) /)

   end subroutine mkdcon


   subroutine gtenv(name,val)
      implicit none
! The passed variables need to be genuine 'string' in fortran's terms and not just an array of type character(len=1)
      character(len=*), intent(in) :: name
      character(len=*), intent(out) :: val
      call get_environment_variable(name,val)
   end subroutine gtenv


   subroutine cwrite( s, i1, i2, nwlin)
      character(len=*), intent(in) :: s
      integer, intent(in) :: i1, i2
      integer, intent(in) :: nwlin
      if (nwlin == 0) then
         write(6,'(a)',advance='no' ) s(i1+1:i2+1)
      else
         write(6,'(a)',advance='yes') s(i1+1:i2+1)
      endif
   end subroutine cwrite

   function frename(oldpath, newpath)
      use posix, only: rename
      use iso_c_binding, only : c_null_char
      implicit none
      character(len=*), intent(in) :: oldpath, newpath
      integer :: i,frename
      frename = rename(trim(oldpath)//c_null_char, trim(newpath)//c_null_char)
   end function frename

   subroutine fsystm(cmd, res)
      use :: iso_c_binding
      use posix
      implicit none

      character( kind = c_char, len=* ), intent(in) :: cmd
      character( kind = c_char, len=len(cmd)+1) :: ccmd
      integer(c_int), intent(out) :: res

      ccmd = trim(cmd)//c_null_char
      res = system(ccmd)
   end subroutine fsystm

   subroutine ptenv(pnam)
      use :: iso_c_binding
      use posix
      implicit none

      character(len=*,kind=c_char), intent(in) :: pnam
      integer(c_int) :: res
      character(len=len(pnam),kind=c_char) :: var, val, cvar, cval
      integer :: eqpos

      eqpos = index(pnam,'=')
      var = pnam(1:eqpos-1)
      val = pnam(eqpos+1:len(pnam))

      cvar = trim(var)//c_null_char
      cval = trim(val)//c_null_char

      res = setenv(cvar, cval, 1)

   end subroutine ptenv


   subroutine sectim(tsec,tusec)
      use :: iso_c_binding
      use posix
      implicit none

      integer, intent(out) :: tsec, tusec
      integer(c_long) :: t,tr
      integer :: dt(8)

      t = 0

      tsec = time(t)
      call date_and_time(values=dt)
      tusec = dt(8)
!     print *,tsec, tusec
   end subroutine sectim




!
!
!   program test
!
!     implicit none
!     integer :: r, ts,tu
!     character(len=20) :: val, varn, name
!     character(len=5) :: hi
!
!     call fsystm('ls /',r)
!     print *,'res: ', r
!
!     name = 'EXT'
!     val = 'trtscs'
!     varn(:) = trim(name)//"="//trim(val)
!
!     print *,'set var: ',varn
!     call ptenv(varn)
!     val = 'dbsv'
!     call gtenv(name,val)
!
! !     call get_environment_variable(name,val)
!     print '("val ot test: ",a)',trim(val)
!
!
!     call sectim(ts,tu)
!     print *,ts,tu
!
!     hi = 'hi'
!     call nlchar(1,hi(3:3))
!     write(6,'(a)',advance='no') trim(hi)
!
!     call flushs(-1)
!
!   end program test
