! Dimitar Pashov <d.pashov@gmail.com>

    module posix

!     use, intrinsic :: iso_c_binding

    implicit none
    private

    interface
        function chdir(path) bind(c)
        ! int chdir(const char *path)
            use iso_c_binding, only : c_char, c_int
            character(kind=c_char) :: path(*)
            integer(c_int) :: chdir
        end function chdir

        function getcwd(buf, size) bind(c)
        ! char *getcwd(char *buf, size_t size);
            use iso_c_binding, only : c_char, c_int, c_ptr
            character(kind=c_char) :: buf(*)
            integer(c_int), value :: size
            type(c_ptr) :: getcwd
        end function getcwd

        function fflush(stream) bind(c)
        ! int fflush(FILE *stream);
        ! Only to be used with c_null_ptr really!
            use :: iso_c_binding, only : c_int, c_ptr
            integer(c_int) :: fflush
            type(c_ptr), value :: stream
        end function fflush

        function system(command) bind(c)
        ! int system(const char *command);
            use :: iso_c_binding, only : c_int, c_char
            integer(c_int) :: system
            character(kind=c_char) :: command(*)
        end function system

        function setenv(name, value, overwrite) bind(c)
        ! int setenv(const char *name, const char *value, int overwrite);
            use :: iso_c_binding, only : c_char, c_int
            character(kind=c_char) :: name
            character(kind=c_char) :: value
            integer(kind=c_int),value :: overwrite
            integer(kind=c_int) :: setenv
        end function setenv

        function time(t) bind(c)
        ! time_t time(time_t *t);
            use :: iso_c_binding, only : c_long
            integer(c_long) :: time
            integer(c_long),value :: t
        end function time

        subroutine exit(rc) bind(c)
        ! void exit(int status);
            use :: iso_c_binding, only : c_int
            integer(c_int), value :: rc
        end subroutine exit

        function rename(oldpath, newpath) bind(c)
!     int rename(const char *oldpath, const char *newpath);
            use :: iso_c_binding, only : c_int, c_char
            character(kind=c_char) :: oldpath, newpath
            integer(c_int) :: rename
        end function rename

    end interface

    public :: chdir, getcwd, fflush, system, setenv, time, exit, rename

    end module posix
