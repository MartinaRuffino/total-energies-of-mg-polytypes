dnl aclocal.m4 
dnl Process this file with autoconf to produce a configure script.
dnl Revisions
dnl   3 Jul 01  Revised compiler environment variable from 'F77' to 'FC'
dnl             Added code to use find fftw, or use alternative cfft3

dnl Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY, to the extent permitted by law; without
dnl even the implied warranty of MERCHANTABILITY or FITNESS FOR A
dnl PARTICULAR PURPOSE.

# --- ACX_CHECK_CC_FLAGS ---
AC_DEFUN(ACX_CHECK_CC_FLAGS,
[
AC_REQUIRE([AC_PROG_CC])
AC_MSG_CHECKING([whether ${CC-cc} accepts $1])
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $1 -c conftest.c 2>&1`"; then
        ac_prog_cc_flags_ok=yes
        $2
else
        ac_prog_cc_flags_ok=no
        $3
fi
rm -r -f conftest*
AC_MSG_RESULT(${ac_prog_cc_flags_ok})
])

# --- AC_MSG_NONL ---
# analogous to AC_MSG_CHECKING, but without 'checking'
dnl AC_MSG_NONL(MESSAGE)
# define(AC_MSG_NONL,
# # for autoconf 2.13
# [echo $ac_n "$1"" $ac_c" 1>&AC_FD_MSG
# echo "configure:__oline__: $1" >&AC_FD_CC])
# AC_MSG_NONL(FEATURE)
# ------------------------
m4_define([AC_MSG_NONL],
[_AS_ECHO([$as_me:$LINENO: $1], AS_MESSAGE_LOG_FD)
_AS_ECHO_N([$1... ])[]dnl
])

# --- ACX_ZDIFF_KSH replaces /bin/sh with /bin/csh in directory/zdiff for selected hosts ---
# arguments :  directory, host
AC_DEFUN(ACX_ZDIFF_KSH,
[
( cd $1;
case "$2" in
  *aix* ) echo patch zdiff in directory $1 : replace /bin/sh with /bin/ksh
          cp -p zdiff zdiff~
          cat zdiff~ | sed s=/bin/sh=/bin/ksh= >zdiff
          chmod +x zdiff;;
esac
)
])

# --- AC_MSG_NOTE(message 
# writes message to file `config.notes'
# If message is "init" AC_MSG_NOTE removes the file
AC_DEFUN(AC_MSG_NOTES,
[
if test x"$1" = x"init"; then
  rm -f config.notes
  touch config.notes
else
  echo $ac_n "$1" >>config.notes
  echo "$2" >>config.notes
fi
]) # 

# --- ACX_VERSION_CHECK checks version ---
# arguments :  required-version-number, variable-to-set, path-of-file-with-version
AC_DEFUN(ACX_VERSION_CHECK,
[

AC_CHECK_FILES($3 , $2=yes , $2=no)
if test "${$2}" = "yes"; then
# echo in test $2 is ${$2}
  $2=`head -1 $3`
# echo in test now $2 is ${$2}
  if test "${$2}" = "$1"; then $2=yes;
  else
    AC_MSG_RESULT([Version ${$2} in file $3 does not match $1])
    $2=no;
  fi
fi
if test "${$2}" = no; then
  AC_MSG_RESULT([>> Setting `echo $2 | sed s/_/-/` to `echo $2 | sed s/_/-/`=no]); fi
])

# --- AC_FC_CMAIN ---
# -------------------
# AC_FC_CMAIN finds the name of the entry point for the fortran linker,
# so that the main program can be a C language program.
#

AC_DEFUN([AC_FC_CMAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_CACHE_CHECK([for name of entry point],
               ac_cv_fc_main,
[AC_LANG_PUSH(C)dnl
 ac_fc_m_save_LIBS=$LIBS
 LIBS="$LIBS $FLIBS"
 ac_cv_fc_main="main" # default entry point name

 list_of_candidates="main MAIN__ __main MAIN _MAIN __MAIN main_ main__ _main MAIN_ lastchoice"

 for ac_func in $list_of_candidates; do

  _AC_COMPILE_IFELSE([AC_LANG_PROGRAM([@%:@undef FC_DUMMY_MAIN
@%:@define main $ac_func])],
  [
   dnl mv conftest.$ac_objext cf77_test.$ac_objext
   ac_cv_fc_main=$ac_func
   ac_success=no
   AC_LANG_PUSH(Fortran)dnl
   AC_TRY_LINK_OBJ($ac_func,
         [ac_success=yes; break])
   AC_LANG_POP(Fortran)dnl])

  dnl echo did main $ac_func succeed ... $ac_success

 done

 ac_cv_fc_main=$ac_func

 dnl echo use \`$ac_cv_fc_main\' for the entry point from the linker
 AC_MSG_NONL([linker uses symbol ])

 rm -r -f conftest*
 LIBS=$ac_fc_m_save_LIBS
 AC_LANG_POP(C)dnl
])

AC_DEFINE_UNQUOTED([FC_CMAIN], $ac_cv_fc_maino,
                   [Define to name of a C language program that corresponds to
                    the main in a fortran program.])

])# AC_FC_CMAIN

# --- AC_FC_SUBPROG(Compilation flags,Fortran-code (optional))
#                            
# AC_FC_SUBPROG attempts to compile and fortran subroutine
#
# Inputs:
# $1 : optional compilation flags
# $2 : optional fortran code
#
# Macro AC_FC_SUBPROG also requires these shell variables to be set:
# $FC
#
# Outputs
# variable ac_success is set as follows:
#   If compilation is successful,
#     ac_success="ok"
#   If compilation is not successful,
#     ac_success="no"
# variable ac_error is returned as command for compilation
#
# File conftest.o is generated
#
AC_DEFUN(AC_FC_SUBPROG,
[
AC_REQUIRE([AC_PROG_FC])

# this is fortran routine
echo '      subroutine fmain' > conftest.f
echo "      $2" >> conftest.f
echo '      end' >> conftest.f

#cat conftest.f

rm -f conftest.o
ac_success=yes
$FC -c $1 conftest.f >/dev/null 2>&1
ac_error="$FC -c $1 conftest.f"
if ! (test -r conftest.o) ; then 
  ac_success=no
fi
])# AC_FC_SUBPROG

# --- AC_FC_PROG(1 : Fortran-code (optional),
#                2 : linker library or extra linker switches (optional),
#                3 : test program standard input (optional),
#                4 : test program command-line arguments (optional))
#                            
# AC_FC_PROG attempts to compile and link, and optionally execute, a fortran program
#
# Inputs:
# $1 : optional fortran code
# $2 : possible extra arguments to be passed to the linker
#      If either of the next two are set, test program is executed using these:
# $3 : standard input fed to test program, if link successful
# $4 : command-line arguments when (test program, if link successful
#
# Macro AC_FC_PROG also requires these shell variables to be set:
# $FC    $FFLAGS   $LK
#
# Outputs
# variable ac_success is set as follows:
#   If link is successful, and optional arguments 7 and 8 are not used,
#     ac_success="ok"
#   If link is successful, and optional arguments 7 or 8 are used,
#     ac_success= output of test program
#   If link is not successful,
#     ac_success="no"
#     ac_error= command that was unsuccessful
#
# Files conftest.o contest are generated
#
AC_DEFUN(AC_FC_PROG,
[
AC_REQUIRE([AC_PROG_FC])

# this is fortran routine, to be called by C main
echo '      program main' > conftest.f
echo "      $1" >> conftest.f
echo '      end' >> conftest.f

rm -f conftest.o contest
ac_success=no
$FC $FFLAGS -c conftest.f >/dev/null 2>&1
ac_error="$LK $FFLAGS $LDFLAGS conftest.o -o contest $2"

if (test -r conftest.o); then 
$LK $FFLAGS $LDFLAGS  conftest.o -o contest $2 >/dev/null 2>&1 
if test -x contest; then 
ac_success=ok
if (test x"$3" = x) && (test x"$4" = x); then :
else
ac_success=`echo $3 | contest $4`
fi
else
if ! (test -r contest) ; then 
  ac_error="$LK $FFLAGS $LDFLAGS  conftest.o -o contest $2"
fi
fi
else 
if ! (test -r conftest.o) ; then 
  ac_error="$FC $FFLAGS -c conftest.f"
fi
fi
])# AC_FC_PROG

# --- AC_CHANGE_ABSOLUTE_PATH_TO_SRCDIR(absolute-path,makefilestyle) ---
# Convert file specfied by absolute path to one relative to $srcdir,
# if the absolute path is within the reference path, or one level higher than it.
#
# $1 : name of variable to change
# $2 : if "makefile", write the altered variable in "Makefile" style: $(srcdir)/pathname
#      Also requires srcdir be set.
#
#  Example:
#    AC_CHANGE_ABSOLUTE_PATH_TO_SRCDIR(CCOMP)
#
AC_DEFUN(AC_CHANGE_ABSOLUTE_PATH_TO_SRCDIR,
[
AC_CHANGE_ABSOLUTE_PATH_TO_RELATIVE($srcdir,${$1},ac_tmp_path)

if ! test "${$1}" = "$ac_tmp_path" ; then
if test "$2" = "makefile" ; then
  AC_MSG_RESULT([convert $1 = ${$1} to \$(srcdir)/$ac_tmp_path])
  $1="\$(srcdir)/$ac_tmp_path"
else
  AC_MSG_RESULT([convert $1 = ${$1} to relative $ac_tmp_path])
  $1="$ac_tmp_path"
fi
else
  AC_MSG_RESULT([No conversion of path for $1 = ${$1}])
fi

])# AC_CHANGE_ABSOLUTE_PATH_TO_SRCDIR

# --- AC_CHANGE_ABSOLUTE_PATH_TO_RELATIVE(1 : reference path for relative path
#                                         2 : absolute-path
#                                         3 : var-for-rel-path
# If the name of absolute-path is within the reference path, or
# within one level higher than the reference path, the third argument
# is re-assigned to a relative path.
#
# Inputs:
# $1 : Reference path that defines starting point for relative path
# $2 : absolute pathname for file
# Output
# $3 : name of the shell variable which to assign
#
#  Example:
#    AC_CHANGE_ABSOLUTE_PATH_TO_RELATIVE($srcdir,$CCOMP,newccomp)
#
AC_DEFUN(AC_CHANGE_ABSOLUTE_PATH_TO_RELATIVE,
[
pathnow=$2
# substitution if pathnow is below reference path
#if ! test "${pathnow/$1/.}" = "$pathnow" ; then
#   $3="${pathnow/$1/.}"
# this one is more portable
if ! test "`echo ${pathnow} | sed s@$1@.@`" = "$pathnow" ; then
  $3="`echo ${pathnow} | sed s@$1@.@`"
# convert ./name to name
  pathnow=${$3}
  if ! test "${pathnow#./}" = "$pathnow" ; then
    $3="${pathnow#./}"
  fi
else
  updir="`(cd $1/..; pwd)`"
#   if ! test "${pathnow/$updir/..}" = "$pathnow" ; then
#     $3=${pathnow/$updir/..}
#   This one more portable
    if ! test "`echo ${pathnow} | sed s@$updir@..@`" = "$pathnow" ; then
      $3="`echo ${pathnow} | sed s@$updir@..@`"
  else
    $3="${pathnow}"
  fi
fi
])# AC_CHANGE_ABSOLUTE_PATH_TO_RELATIVE


# --- AC_FC_PROG_WITH_C_MAIN(1 : C-entry-point (required),
#                            2 : C-arg-decl (required if entry point has args),
#                            3 : C-code (optional),
#                            4 : Fortran-entry-point (required),
#                            5 : Fortran-code (optional),
#                            6 : linker library or extra linker switches (optional),
#                            7 : test program standard input (optional),
#                            8 : test program command-line arguments (optional))
#                            
# AC_FC_PROG_WITH_C_MAIN attempts to compile and link a fortran program called by a C main.
#
# Inputs:
# $1 : symbol for C main entry point
# $2 : symbol for possible declaration of arguments (optional)
# $3 : possible C code to be inserted before fortran call 
# $4 : symbol for fortran entry point corresponding to mangled symbol for 'fmain'
# $5 : optional fortran code
# $6 : possible extra arguments to be passed to the linker
#      If either of the next two are set, test program is executed using these:
# $7 : standard input fed to test program, if link successful
# $8 : command-line arguments when (test program, if link successful
#
# Macro AC_FC_PROG_WITH_C_MAIN also requires these shell variables to be set:
# ${CC-cc}   $FC    $FFLAGS   $LK
#
# Outputs
# variable ac_success is set as follows:
#   If link is successful, and optional arguments 7 and 8 are not used,
#     ac_success="ok"
#   If link is successful, and optional arguments 7 or 8 are used,
#     ac_success= output of test program
#   If link is not successful,
#     ac_success="no"
#     ac_error= command that was unsuccessful
#
# Files contest.o conftest.o contest are generated
#
AC_DEFUN(AC_FC_PROG_WITH_C_MAIN,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_FC])

# this is C main routine, calling fortran
echo "void $1" > contest.c
echo "$2" >>contest.c
echo '{' "$3" >>contest.c
echo "$4 ;" >>contest.c
echo '}' >>contest.c

# this is fortran routine, to be called by C main
echo '      subroutine fmain' > conftest.f
echo "      $5" >> conftest.f
echo '      end' >> conftest.f

rm -f contest.o conftest.o contest
ac_success=no
${CC-cc} $CFLAGS -c contest.c >/dev/null 2>&1
$FC $FFLAGS -c conftest.f >/dev/null 2>&1
ac_error="$LK $FFLAGS $LDFLAGS contest.o conftest.o -o contest $6"

if (test -r contest.o) && (test -r conftest.o); then 
$LK $FFLAGS $LDFLAGS contest.o conftest.o -o contest $6 >/dev/null 2>&1 
if test -x contest; then 
ac_success=ok
if (test x"$7" = x) && (test x"$8" = x); then :
else
ac_success=`echo $7 | contest $8`
fi
else
if ! (test -r contest) ; then 
  ac_error="$LK $FFLAGS $LDFLAGS contest.o conftest.o -o contest $6"
fi
fi
else 
if ! (test -r contest.o) ; then 
  ac_error="${CC-cc} -c contest.c"
fi
if ! (test -r conftest.o) ; then 
  ac_error="$FC $FFLAGS -c conftest.f"
fi
fi
])# AC_FC_PROG_WITH_C_MAIN

# --- AC_FC_F2C_LINKS sets variables for connecting fortran and C programs ---
# sets the following variables:
# ac_cv_fc_main names the entry point the linker assigns as the beginning
#   of the program.  Usually this variable is 'main' (as it would be if
#   you link a normal C program, but because usually the fortran compiler 
#   is used for the linking, some compilers use a different name.
#
# ac_cv_fc_mangling, which describes the fortran name-mangling scheme
#   This variable contains three fields, separated by commas:
#   lower case / upper case:
#      case translation of the Fortran symbols
#   underscore / no underscore:
#      whether the compiler appends "_" to symbol names
#   extra underscore / no extra underscore:
#      whether the compiler appends an extra "_" to symbol names already
#      containing at least one underscore

AC_DEFUN(AC_FC_F2C_LINKS,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_FC])
# echo FFLAGS is $FFLAGS

AC_FC_CMAIN
_AC_FC_NAME_MANGLING
])# AC_FC_F2C_LINKS


# --- ACX_MAKE_ARCHIVE_LIBRARY makes an archive library ---
# arguments :  directory-name
AC_DEFUN(ACX_MAKE_ARCHIVE_LIBRARY,
[

if test "$enable_$1" = "yes"; then
  thisdir=`pwd`
  AC_MSG_RESULT([creating soft link $1/Make.inc to Make.inc])
  (cd $1 && rm -f Make.inc && $LN_S $thisdir/Make.inc Make.inc)
  if test -r $1/test; then
    AC_MSG_RESULT([creating soft link $1/test/zdiff to testing/zdiff])
    (cd $1/test && rm -f zdiff && $LN_S $thisdir/testing/zdiff zdiff)
    AC_MSG_RESULT([creating soft link $1/test/add0 to testing/add0])
    (cd $1/test && rm -f add0 && $LN_S $thisdir/testing/add0 add0)
  fi
  makefiles="$makefiles $1/Makefile"
  AC_MSG_RESULT([Creating file $1/Makefile.in using startup/Maksmakfile:])
dnl  always call with initF90
dnl  MaksmakfileSW="--in"
dnl  if test "$FC_IS_F90" = "yes"; then MaksmakfileSW="$MaksmakfileSW --initF90"
dnl  fi
  MaksmakfileSW="--in --initF90"

  echo "(cd $1 && $thisdir/startup/Maksmakfile $MaksmakfileSW --mnemonic $MNEMONIC >Makefile.in)"
  (cd $1 && $thisdir/startup/Maksmakfile $MaksmakfileSW --mnemonic $MNEMONIC >Makefile.in)
fi

])

# same as ACX_MAKE_ARCHIVE_LIBRARY, but creates a split archive
AC_DEFUN(ACX_MAKE_ARCHIVE_LIBRARYS,
[

if test "$enable_$1" = "yes"; then
  thisdir=`pwd`
  AC_MSG_RESULT([creating soft link $1/Make.inc to Make.inc])
  (cd $1 && rm -f Make.inc && $LN_S $thisdir/Make.inc Make.inc)
#   AC_MSG_RESULT([creating soft link $1/subs2.a to $1/subs.a])
#   (cd $1 && rm -f subs2.a && $LN_S $thisdir/$1/subs.a subs2.a)
  if test -r $1/test; then
    AC_MSG_RESULT([creating soft link $1/test/zdiff to testing/zdiff])
    (cd $1/test && rm -f zdiff && $LN_S $thisdir/testing/zdiff zdiff)
    AC_MSG_RESULT([creating soft link $1/test/add0 to testing/add0])
    (cd $1/test && rm -f add0 && $LN_S $thisdir/testing/add0 add0)
  fi
  makefiles="$makefiles $1/Makefile"
  AC_MSG_RESULT([Creating file $1/Makefile.in using startup/Maksmakfile:])
  echo "(cd $1 && $thisdir/startup/Maksmakfile --in --split --mnemonic $MNEMONIC >Makefile.in)"
  (cd $1 && $thisdir/startup/Maksmakfile --in --split --mnemonic $MNEMONIC >Makefile.in)
fi

])

# --- AC_ADD_TO_CCOMP_SWITCHES(1 : string to compare against match-string; see arg 2
#                              2 : match-string.  If matches first argument,
#                                  ccomp switch will be included as 'defined'.
#                                  Otherwise switch will be 'undefined'.
#                              3 : if "quiet" then do nothing if ccomp switch is 'undefined'
#                              4 : Token to incorporate as ccomp switch
#                            
# AC_ADD_TO_CCOMP_SWITCHES appends to variable CCOMP_SWITCHES either:
# -dToken  or  -uToken   depending on whether first and second arguments match.
#
AC_DEFUN(AC_ADD_TO_CCOMP_SWITCHES,
[
  if test "$1" = "$2"; then
    CCOMP_SW="$CCOMP_SW -d$4"
  else
    if test x"$3" = x"quiet"; then : ; else
      CCOMP_SW="$CCOMP_SW -u$4"
    fi
  fi
])

# --- AC_CHECK_FOR_F90 checks if the fortran compiler allows f90 constructs ---
# Returns FC_IS_F90="yes" or "no"
# 27 Jul 09 bug fix Dimitar
AC_DEFUN(AC_CHECK_FOR_F90,
[

  echo '      program alloc' >conftest.f
  echo '      integer,allocatable:: igvx(:,:)' >>conftest.f
  echo '      allocate(igvx(3,4))' >>conftest.f
  echo '      end program alloc' >>conftest.f
  if `$FC -c $FFLAGS conftest.f >/dev/null 2>&1` ; then
    FC_IS_F90=yes
  else
    FC_IS_F90=no
  fi
  rm -f conftest.f conftest.o
])

# --- AC_CHECK_FOR_QUAD checks if the fortran compiler allows QUAD precision ---
# Returns FC_USES_QUAD="yes" or "no"
AC_DEFUN(AC_CHECK_FOR_QUAD,
[
  echo '      real*16 zn(2)' >conftest.f
  echo '      call zero(zn)' >>conftest.f
  echo '      zn(2) = zn(1)' >>conftest.f
  echo '      end' >>conftest.f
  if `$FC -c $FFLAGS conftest.f $LIBLOC >/dev/null 2>&1` ; then
    FC_USES_QUAD=yes
  else
    FC_USES_QUAD=no
  fi
  rm -f conftest.f conftest.o
])

# --- AC_CCOMP executes ccomp on specified fortran files ---
# If file with extension changed to `.fortran' does NOT exist,
#    AC_CCOMP renames existing file to one with that extension.
#    (If the file also doesn't exist, AC_CCOMP aborts with error)
# Next, AC_CCOMP invokes $CCOMP on the .fortran file, creating
#    a new file with specified name.
# $1 : list of files
# $2 : list of CCOMP directives to define or undefine
# $3 : list of CCOMP directives to undefine
# $4 : (optional) add to this name list of files renamed to .fortran
AC_DEFUN(AC_CCOMP,
[
CCOMP_UDEF=
CCOMP_DEF=
ifelse($3, , , 
[  
#  AC_MSG_RESULT([disabling $3 in files "$1" ...])
  for i in $3; do
  CCOMP_UDEF="$CCOMP_UDEF -u$i"
  done
])
ifelse($2, , , 
[  
#  AC_MSG_RESULT([enabling $2  in files "$1" ...])
  for i in $2; do
  CCOMP_DEF="$CCOMP_DEF -u$i -d$i"
  done
])
  for i in $1; do
    j=`echo $i | sed 's:[.]f$::'`.fortran
    dnl 9 Jul 03 always move .f file to .fortran
#    if test -r "$j"; then : ; else
     if mv $i $j ; then
       AC_MSG_RESULT([  mv $i $j])
     else
       AC_MSG_RESULT([configure aborting ... file $i could not be renamed to $j])
       exit
     fi
     ifelse($4, , : , [$4="[$]$4 $j"])
#   fi
    if test -r "$j"; then
     AC_MSG_RESULT($CCOMP $CCOMP_UDEF $CCOMP_DEF $j $i)
     $CCOMP $CCOMP_UDEF $CCOMP_DEF $j $i >/dev/null
    else
     AC_MSG_RESULT([configure aborting ... file $j was not readable])
    fi
  done
]
)

# --- AC_CHECK_FOR_AUTO_ARRAY checks if the fortran compiler allows automatic arrays ---
# Returns FC_AUTOARRAY="yes" or "no"
AC_DEFUN(AC_CHECK_FOR_AUTO_ARRAY,
[
cat > conftest.f <<EOF
      integer n
      double precision y
      n = 12
      call snot(n,y)
      end
      subroutine snot(n,y)
      double precision x(n)
      call dummy(x,n)
      y = x(n)
      end
EOF
if `$FC -c $FFLAGS conftest.f $LIBLOC >/dev/null 2>&1` ; then
  FC_AUTOARRAY=yes
else
  FC_AUTOARRAY=no
fi
rm -f conftest.f conftest.o
])

# --- AC_CHECK_FOR_FORTRAN_VARIABLE_STRINGS Check whether fortran compiler allows strings of variable length ---
# Returns FC_VSTRINGS="yes" or "no"
AC_DEFUN(AC_CHECK_FOR_FORTRAN_VARIABLE_STRINGS,
[
cat > conftest.f <<EOF
      subroutine snot(i,strn)
      character*(i) strn
      end
EOF
if `$FC -c $FFLAGS conftest.f >/dev/null 2>&1` ; then
  FC_VSTRINGS=yes
else
  FC_VSTRINGS=no
fi
rm -f conftest.f conftest.o
])

# --- ACX_FIND_ROUTINE_IN_LIB looks whether a particular routine
#     can be found when linking with supplied libraries
# $1 : routine name
# $2 : library in which to seek routine name
# $3 : If 'cmain', use C language main for entry point; else use Fortran 
# $4 : optional additional libraries
# $5 : action-if-found
# $6 : action-if-not-found
# (old)
# $1 : library in which to find routine
# $2 : optional additional libraries
# $3 : action-if-found
# $4 : action-if-not-found
# 
# Calls AC_FC_PROG_WITH_C_MAIN, which requires these variables:
# $FC    $FFLAGS   $LK
# Also if C main is to be used, ${CC-cc}, $ac_cv_fc_main, and  $fmain must be defined

AC_DEFUN(ACX_FIND_ROUTINE_IN_LIB,
[

  AC_MSG_CHECKING([for $1 in $2])
  if test x"$3" = x"cmain"; then
    AC_FC_PROG_WITH_C_MAIN([$ac_cv_fc_main(argc,argv)],,,
                           [$fmain()],[call $1],[$2 $4],,,)
  else
    AC_FC_PROG([call $1],[$2 $4],,,)
  fi

  if test "$ac_success" = "ok"; then
     AC_MSG_RESULT([ok])
     ifelse([$5], , :, [$5])
  else
    AC_MSG_RESULT([no])
    ifelse([$6], , , [$6])
  fi

# cat conftest.f
# echo $ac_success
# echo $ac_error

])# ACX_FIND_ROUTINE_IN_LIB

# --- ACX_WARN_LIB_MISSING checks if a library is present by trying to ---
#     link a call to a representative routine.
#     It is similar to ACX_FIND_ROUTINE_IN_LIB, but output is different
# $1: sets have_$1 to "yes" if link is successful, to "no" if not
# $2: routine name to find in link that verifies library is present
# $3: optional libraries (in addition to LIBLOC) to include in link
AC_DEFUN(ACX_WARN_LIB_MISSING,
[
  AC_FC_PROG_WITH_C_MAIN([$ac_cv_fc_main(argc,argv)],,,
  		         [$fmain()],[call $2],[$3 $LIBLOC],,,)

  if test "$ac_success" = "ok"; then
     AC_MSG_RESULT([LIBLOC apparently contains $1 library]);
     have_$1=yes
  else
    AC_MSG_RESULT([LIBLOC apparently does not contain $1 library]);
    have_$1=no
  fi
])

# AC_TRY_LINK_OBJ(OBJ, ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
# ------------------------------------------------------------
# Try to link a program that calls FUNC, handling GCC builtins.  If
# the link succeeds, execute ACTION-IF-FOUND; otherwise, execute
# ACTION-IF-NOT-FOUND.
AC_DEFUN([AC_TRY_LINK_OBJ],
[AC_LINK_OBJ_IFELSE([AC_LANG_CALL([], [$1])], [$2], [$3])])


# AC_LINK_OBJ_IFELSE(PROGRAM, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------
# Try to link PROGRAM.  Requires that the compiler for the current
# language was checked for, hence do not use this macro in macros looking
# for a compiler.
AC_DEFUN([AC_LINK_OBJ_IFELSE],
[AC_LANG_COMPILER_REQUIRE()dnl
_AC_LINK_OBJ_IFELSE($@)])

# _AC_LINK_OBJ_IFELSE(PROGRAM, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
# Try to link PROGRAM.
# This macro can be used during the selection of a compiler.
m4_define([_AC_LINK_OBJ_IFELSE],
[rm -f conftest$ac_exeext
 AS_IF([AC_TRY_EVAL(ac_linko) &&
        AC_TRY_COMMAND([test -s conftest$ac_exeext])],
      [$2],
      [echo "$as_me: failed program was:" >&AS_MESSAGE_LOG_FD
m4_ifvaln([$3], [$3])dnl])[]dnl
rm -f conftest.$ac_objext conftest$ac_exeext m4_ifval([$1], [conftest.$ac_ext])[]dnl
])# _AC_LINK_OBJ_IFELSE


# ------------------------ fortran.m4, with modifications ---------------
# This file is part of Autoconf.
# Fortran languages support.
# Copyright 2001
# Free Software Foundation, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
#
# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.


# _AC_LIST_MEMBER_IF(ELEMENT, LIST, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------------
#
# Processing the elements of a list is tedious in shell programming,
# as lists tend to be implemented as space delimited strings.
#
# This macro searches LIST for ELEMENT, and executes ACTION-IF-FOUND
# if ELEMENT is a member of LIST, otherwise it executes
# ACTION-IF-NOT-FOUND.
AC_DEFUN([_AC_LIST_MEMBER_IF],
[dnl Do some sanity checking of the arguments.
m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  ac_exists=false
  for ac_i in $2; do
    if test x"$1" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  AS_IF([test x"$ac_exists" = xtrue], [$3], [$4])[]dnl
])# _AC_LIST_MEMBER_IF


# _AC_LINKER_OPTION(LINKER-OPTIONS, SHELL-VARIABLE)
# -------------------------------------------------
#
# Specifying options to the compiler (whether it be the C, C++ or
# Fortran compiler) that are meant for the linker is compiler
# dependent.  This macro lets you give options to the compiler that
# are meant for the linker in a portable, compiler-independent way.
#
# This macro take two arguments, a list of linker options that the
# compiler should pass to the linker (LINKER-OPTIONS) and the name of
# a shell variable (SHELL-VARIABLE).  The list of linker options are
# appended to the shell variable in a compiler-dependent way.
#
# For example, if the selected language is C, then this:
#
#   _AC_LINKER_OPTION([-R /usr/local/lib/foo], foo_LDFLAGS)
#
# will expand into this if the selected C compiler is gcc:
#
#   foo_LDFLAGS="-Xlinker -R -Xlinker /usr/local/lib/foo"
#
# otherwise, it will expand into this:
#
#   foo_LDFLAGS"-R /usr/local/lib/foo"
#
# You are encouraged to add support for compilers that this macro
# doesn't currently support.
# FIXME: Get rid of this macro.
AC_DEFUN([_AC_LINKER_OPTION],
[if test "$ac_compiler_gnu" = yes; then
  for ac_link_opt in $1; do
    $2="[$]$2 -Xlinker $ac_link_opt"
  done
else
  $2="[$]$2 $1"
fi[]dnl
])# _AC_LINKER_OPTION



## ----------------------- ##
## 1. Language selection.  ##
## ----------------------- ##


# ----------------------------- #
# 1d. The Fortran language.  #
# ----------------------------- #


# AC_LANG(Fortran)
# -------------------
m4_define([AC_LANG(Fortran)],
[ac_ext=f
ac_compile='$FC -c $FFLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$FC -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_linko='$FC -o conftest$ac_exeext $FFLAGS $LDFLAGS conftest.$ac_objext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_f77_compiler_gnu
])


# AC_LANG_FORTRAN77
# -----------------
AU_DEFUN([AC_LANG_FORTRAN77], [AC_LANG(Fortran)])


# _AC_LANG_ABBREV(Fortran)
# ---------------------------
m4_define([_AC_LANG_ABBREV(Fortran)], [f77])



## ---------------------- ##
## 2.Producing programs.  ##
## ---------------------- ##


# ------------------------ #
# 2d. Fortran sources.  #
# ------------------------ #

# AC_LANG_SOURCE(Fortran)(BODY)
# --------------------------------
# FIXME: Apparently, according to former AC_TRY_COMPILER, the CPP
# directives must not be included.  But AC_TRY_RUN_NATIVE was not
# avoiding them, so?
m4_define([AC_LANG_SOURCE(Fortran)],
[$1])


# AC_LANG_PROGRAM(Fortran)([PROLOGUE], [BODY])
# -----------------------------------------------
# Yes, we discard the PROLOGUE.
m4_define([AC_LANG_PROGRAM(Fortran)],
[m4_ifval([$1],
       [m4_warn([syntax], [$0: ignoring PROLOGUE: $1])])dnl
      program main
$2
      end])


# AC_LANG_CALL(Fortran)(PROLOGUE, FUNCTION)
# --------------------------------------------
# FIXME: This is a guess, help!
m4_define([AC_LANG_CALL(Fortran)],
[AC_LANG_PROGRAM([$1],
[      call $2])])




## -------------------------------------------- ##
## 3. Looking for Compilers and Preprocessors.  ##
## -------------------------------------------- ##


# ----------------------------- #
# 3d. The Fortran compiler.  #
# ----------------------------- #


# AC_LANG_PREPROC(Fortran)
# ---------------------------
# Find the Fortran preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran)],
[m4_warn([syntax],
         [$0: No preprocessor defined for ]_AC_LANG)])


# AC_LANG_COMPILER(Fortran)
# ----------------------------
# Find the Fortran compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran)],
[AC_REQUIRE([AC_PROG_FC])])


# ac_cv_prog_g77
# --------------
# We used to name the cache variable this way.
AU_DEFUN([ac_cv_prog_g77],
[ac_cv_f77_compiler_gnu])


# AC_PROG_FC([COMPILERS...])
# ---------------------------
# COMPILERS is a space separated list of Fortran compilers to search
# for.  Fortran 95 isn't strictly backwards-compatible with Fortran,
# but `f95' is worth trying.
#
# Compilers are ordered by
#  1. F90, F95, F77
#  2. Good/tested native compilers, bad/untested native compilers
#  3. Wrappers around f2c go last.
#
# `fort77' and `fc' are wrappers around `f2c', `fort77' being better.
# It is believed that under HP-UX `fort77' is the name of the native
# compiler.  On some Cray systems, fort77 is a native compiler.
# cf77 and cft77 are (older) Cray F77 compilers.
# frt is the Fujitsu F77 compiler.
# pgf77 and pgf90 are the Portland Group F77 and F90 compilers.
# xlf/xlf90/xlf95 are IBM (AIX) F77/F90/F95 compilers.
# lf95 is the Lahey-Fujitsu compiler.
# fl32 is the Microsoft Fortran "PowerStation" compiler.
# af77 is the Apogee F77 compiler for Intergraph hardware running CLIX.
# epcf90 is the "Edinburgh Portable Compiler" F90.
# fort is the Compaq Fortran 90 (now 95) compiler for Tru64 and Linux/Alpha.
AC_DEFUN([AC_PROG_FC],
[AC_LANG_PUSH(Fortran)dnl
AC_ARG_VAR([FC],    [Fortran compiler command])dnl
AC_ARG_VAR([FFLAGS], [Fortran compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
dnl AC_CHECK_TOOLS(FC,
dnl    [m4_default([$1],
dnl              [g77 f77 xlf cf77 cft77 frt pgf77 fl32 af77 fort77 f90 xlf90 pgf90 epcf90 f95 fort xlf95 lf95 g95 ifc fc])])

AC_CHECK_TOOLS(FC,
   [m4_default([$1],
             [f90 cf90 cf77 fort xlf xlf90 ifc pgf90 epcf90 f95 xlf95 lf95 pgf77 fl32 af77 frt f77 g77 cf77 cft77 fort77 fc])])

# Provide some information about the compiler.
echo "$as_me:__oline__:" \
     "checking for _AC_LANG compiler version" >&AS_MESSAGE_LOG_FD
ac_compiler=`set X $ac_compile; echo $[2]`
_AC_EVAL([$ac_compiler --version </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -v </dev/null >&AS_MESSAGE_LOG_FD])
_AC_EVAL([$ac_compiler -V </dev/null >&AS_MESSAGE_LOG_FD])

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
# If we don't use `.F' as extension, the preprocessor is not run on the
# input file.
ac_save_ext=$ac_ext
ac_ext=F
_AC_LANG_COMPILER_GNU
ac_ext=$ac_save_ext
G77=`test $ac_compiler_gnu = yes && echo yes`
_AC_PROG_FC_G
AC_LANG_POP(Fortran)dnl
])# AC_PROG_FC


# _AC_PROG_FC_G
# --------------
# Check whether -g works, even if FFLAGS is set, in case the package
# plays around with FFLAGS (such as to build both debugging and normal
# versions of a library), tasteless as that idea is.
m4_define([_AC_PROG_FC_G],
[ac_test_FFLAGS=${FFLAGS+set}
ac_save_FFLAGS=$FFLAGS
FFLAGS=
AC_CACHE_CHECK(whether $FC accepts -g, ac_cv_prog_fc_g,
[FFLAGS=-g
_AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_fc_g=yes],
[ac_cv_prog_fc_g=no])
])
if test "$ac_test_FFLAGS" = set; then
  FFLAGS=$ac_save_FFLAGS
elif test $ac_cv_prog_fc_g = yes; then
  if test "$G77" = yes; then
    FFLAGS="-g -O2"
  else
    FFLAGS="-g"
  fi
else
  if test "$G77" = yes; then
    FFLAGS="-O2"
  else
    FFLAGS=
  fi
fi[]dnl
])# _AC_PROG_FC_G


# AC_PROG_FC_C_O
# ---------------
# Test if the Fortran compiler accepts the options `-c' and `-o'
# simultaneously, and define `FC_NO_MINUS_C_MINUS_O' if it does not.
#
# The usefulness of this macro is questionable, as I can't really see
# why anyone would use it.  The only reason I include it is for
# completeness, since a similar test exists for the C compiler.
AC_DEFUN([AC_PROG_FC_C_O],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_CACHE_CHECK([whether $FC understand -c and -o together],
               [ac_cv_prog_fc_c_o],
[AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
# We test twice because some compilers refuse to overwrite an existing
# `.o' file with `-o', although they will create one.
ac_try='$FC $FFLAGS -c conftest.$ac_ext -o conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
if AC_TRY_EVAL(ac_try) &&
     test -f conftest.$ac_objext &&
     AC_TRY_EVAL(ac_try); then
  ac_cv_prog_fc_c_o=yes
else
  ac_cv_prog_fc_c_o=no
fi
rm -r -f conftest*])
if test $ac_cv_prog_fc_c_o = no; then
  AC_DEFINE(FC_NO_MINUS_C_MINUS_O, 1,
            [Define to 1 if your Fortran compiler doesn't accept
             -c and -o together.])
fi
])# AC_PROG_FC_C_O





## ------------------------------- ##
## 4. Compilers' characteristics.  ##
## ------------------------------- ##


# ---------------------------------------- #
# 4d. Fortran compiler characteristics. #
# ---------------------------------------- #


# _AC_PROG_FC_V_OUTPUT([FLAG = $ac_cv_prog_fc_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (which default value, $ac_cv_prog_fc_v, is computed by
# _AC_PROG_FC_V), and return the output in $ac_fc_v_output.  This
# output is processed in the way expected by AC_FC_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_FC_V_OUTPUT],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_LANG_PUSH(Fortran)dnl

AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran compiler in order to get
# "verbose" output that we can then parse for the Fortran linker
# flags.
ac_save_FFLAGS=$FFLAGS
FFLAGS="$FFLAGS m4_default([$1], [$ac_cv_prog_fc_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_fc_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_fc_v_output" >&AS_MESSAGE_LOG_FD
FFLAGS=$ac_save_FFLAGS

rm -r -f conftest*
AC_LANG_POP(Fortran)dnl

# If we are using xlf then replace all the commas with spaces.
if echo $ac_fc_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_fc_v_output=`echo $ac_fc_v_output | sed 's/,/ /g'`
fi

# gfortran lib crt1.o conflicts with gcc crt1.10.6.o, darwin OS
if echo $ac_fc_v_output | grep i686-apple-darwin  >/dev/null 2>&1; then
  ac_fc_v_output=`echo $ac_fc_v_output | sed 's/-lcrt1.o //g'`
fi

# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ac_fc_v_output="`echo $ac_fc_v_output | 
	grep 'LPATH is:' |
	sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ac_fc_v_output"

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_fc_v_output | grep cft90 >/dev/null 2>&1; then
  ac_fc_v_output=`echo $ac_fc_v_output | sed "s/\"//g"`
fi[]dnl
])# _AC_PROG_FC_V_OUTPUT


# _AC_PROG_FC_V
# --------------
#
# Determine the flag that causes the Fortran compiler to print
# information of library and object files (normally -v)
# Needed for AC_FC_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_FC_V],
[AC_CACHE_CHECK([how to get verbose linking output from $FC],
                [ac_cv_prog_fc_v],
[AC_LANG_ASSERT(Fortran)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_fc_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_FC_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_fc_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_fc_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_fc_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $FC])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_FC_V


# AC_FC_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran
# intrinsic and run-time libraries that are required to successfully
# link a Fortran program or shared library.  The output variable
# FLIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran intrinsic and run-time libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran libraries.  Hence, the macro
# "AC_FC_LIBRARY_LDFLAGS" was created to determine these Fortran
# libraries.
#
# This macro was packaged in its current form by Matthew D. Langston.
# However, nearly all of this macro came from the "OCTAVE_FLIBS" macro
# in "octave-2.0.13/aclocal.m4", and full credit should go to John
# W. Eaton for writing this extremely useful macro.  Thank you John.
AC_DEFUN([AC_FC_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran)dnl
_AC_PROG_FC_V
AC_CACHE_CHECK([for Fortran libraries], ac_cv_flibs,
[if test "x$FLIBS" != "x"; then
  ac_cv_flibs="$FLIBS" # Let the user override the test.
else

_AC_PROG_FC_V_OUTPUT

ac_cv_flibs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_fc_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
              ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_flibs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_flibs="$ac_cv_flibs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
             "" | -*);;
             *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_flibs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_flibs="$ac_cv_flibs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_flibs, ,
                             ac_cv_flibs="$ac_cv_flibs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_fc_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_flibs)
      ;;
esac
fi # test "x$FLIBS" = "x"
])
FLIBS="$ac_cv_flibs"
AC_SUBST(FLIBS)
AC_LANG_POP(Fortran)dnl
])# AC_FC_LIBRARY_LDFLAGS


# AC_FC_DUMMY_MAIN([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------
#
# Detect name of dummy main routine required by the Fortran libraries,
# (if any) and define FC_DUMMY_MAIN to this name (which should be
# used for a dummy declaration, if it is defined).  On some systems,
# linking a C program to the Fortran library does not work unless you
# supply a dummy function called something like MAIN__.
#
# Execute ACTION-IF-NOT-FOUND if no way of successfully linking a C
# program with the FC libs is found; default to exiting with an error
# message.  Execute ACTION-IF-FOUND if a dummy routine name is needed
# and found or if it is not needed (default to defining FC_DUMMY_MAIN
# when needed).
#
# What is technically happening is that the Fortran libraries provide
# their own main() function, which usually initializes Fortran I/O and
# similar stuff, and then calls MAIN__, which is the entry point of
# your program.  Usually, a C program will override this with its own
# main() routine, but the linker sometimes complain if you don't
# provide a dummy (never-called) MAIN__ routine anyway.
#
# Of course, programs that want to allow Fortran subroutines to do
# I/O, etcetera, should call their main routine MAIN__() (or whatever)
# instead of main().  A separate autoconf test (AC_FC_MAIN) checks
# for the routine to use in this case (since the semantics of the test
# are slightly different).  To link to e.g. purely numerical
# libraries, this is normally not necessary, however, and most C/C++
# programs are reluctant to turn over so much control to Fortran.  =)
#
# The name variants we check for are (in order):
#   MAIN__ (g77, MAIN__ required on some systems; IRIX, MAIN__ optional)
#   MAIN_, __main (SunOS)
#   MAIN _MAIN __MAIN main_ main__ _main (we follow DDD and try these too)
AC_DEFUN([AC_FC_DUMMY_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
m4_define([_AC_LANG_PROGRAM_C_FC_HOOKS],
[#ifdef FC_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int FC_DUMMY_MAIN() { return 1; }
#endif
])
AC_CACHE_CHECK([for dummy main to link with Fortran libraries],
               ac_cv_fc_dummy_main,
[AC_LANG_PUSH(C)dnl
 ac_fc_dm_save_LIBS=$LIBS
 LIBS="$LIBS $FLIBS"

 # First, try linking without a dummy main:
 AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                [ac_cv_fc_dummy_main=none],
                [ac_cv_fc_dummy_main=unknown])

 if test $ac_cv_fc_dummy_main = unknown; then
   for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[@%:@define FC_DUMMY_MAIN $ac_func]])],
                    [ac_cv_fc_dummy_main=$ac_func; break])
   done
 fi
 rm -r -f conftest*
 LIBS=$ac_fc_dm_save_LIBS
 AC_LANG_POP(C)dnl
])
FC_DUMMY_MAIN=$ac_cv_fc_dummy_main
AS_IF([test "$FC_DUMMY_MAIN" != unknown],
      [m4_default([$1],
[if test $FC_DUMMY_MAIN != none; then
  AC_DEFINE_UNQUOTED([FC_DUMMY_MAIN], $FC_DUMMY_MAIN,
                     [Define to dummy `main' function (if any) required to
                      link to the Fortran libraries.])
fi])],
      [m4_default([$2],
                [AC_MSG_ERROR([linking to Fortran libraries from C fails])])])
])# AC_FC_DUMMY_MAIN


# AC_FC_MAIN
# -----------
# Define FC_MAIN to name of alternate main() function for use with
# the Fortran libraries.  (Typically, the libraries may define their
# own main() to initialize I/O, etcetera, that then call your own
# routine called MAIN__ or whatever.)  See AC_FC_DUMMY_MAIN, above.
# If no such alternate name is found, just define FC_MAIN to main.
#
AC_DEFUN([AC_FC_MAIN],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_CACHE_CHECK([for alternate main to link with Fortran libraries],
               ac_cv_fc_main,
[AC_LANG_PUSH(C)dnl
 ac_fc_m_save_LIBS=$LIBS
 LIBS="$LIBS $FLIBS"
 ac_cv_fc_main="main" # default entry point name

 for ac_func in MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ main__ _main; do
   AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@undef FC_DUMMY_MAIN
@%:@define main $ac_func])],
                  [ac_cv_fc_main=$ac_func; break])
 done
 rm -r -f conftest*
 LIBS=$ac_fc_m_save_LIBS
 AC_LANG_POP(C)dnl
])
AC_DEFINE_UNQUOTED([FC_MAIN], $ac_cv_fc_main,
                   [Define to alternate name for `main' routine that is
                    called from a `main' in the Fortran libraries.])
])# AC_FC_MAIN


# _AC_FC_NAME_MANGLING
# ---------------------
# Test for the name mangling scheme used by the Fortran compiler.
#
# Sets ac_cv_fc_mangling. The value contains three fields, separated
# by commas:
#
# lower case / upper case:
#    case translation of the Fortran symbols
# underscore / no underscore:
#    whether the compiler appends "_" to symbol names
# extra underscore / no extra underscore:
#    whether the compiler appends an extra "_" to symbol names already
#    containing at least one underscore
#
AC_DEFUN([_AC_FC_NAME_MANGLING],
[AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])dnl
AC_REQUIRE([AC_FC_DUMMY_MAIN])dnl
AC_CACHE_CHECK([for Fortran name-mangling scheme],
               ac_cv_fc_mangling,
[AC_LANG_PUSH(Fortran)dnl
# Compiles this routine without attempting to link it
AC_COMPILE_IFELSE(
[      subroutine foobar()
      return
      end
      subroutine foo_bar()
      return
      end],
[mv conftest.$ac_objext cf77_test.$ac_objext

  AC_LANG_PUSH(C)dnl

  ac_save_LIBS=$LIBS
  LIBS="cf77_test.$ac_objext $LIBS $FLIBS"

  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      AC_TRY_LINK_FUNC($ac_func,
         [ac_success=yes; break 2])
    done
  done

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
        foobar)
           ac_case=lower
           ac_foo_bar=foo_bar
           ;;
        FOOBAR)
           ac_case=upper
           ac_foo_bar=FOO_BAR
           ;;
     esac

     ac_success_extra=no
     for ac_extra in "" "_"; do
        ac_func="$ac_foo_bar$ac_underscore$ac_extra"
        AC_TRY_LINK_FUNC($ac_func,
        [ac_success_extra=yes; break])
     done

     if test "$ac_success_extra" = "yes"; then
	ac_cv_fc_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_fc_mangling="$ac_cv_fc_mangling, no underscore"
	else
           ac_cv_fc_mangling="$ac_cv_fc_mangling, underscore"
        fi
        if test -z "$ac_extra"; then
           ac_cv_fc_mangling="$ac_cv_fc_mangling, no extra underscore"
	else
           ac_cv_fc_mangling="$ac_cv_fc_mangling, extra underscore"
        fi
      else
	ac_cv_fc_mangling="unknown"
      fi
  else
     ac_cv_fc_mangling="unknown"
  fi

  LIBS=$ac_save_LIBS
  AC_LANG_POP(C)dnl
  rm -f cf77_test* conftest*])
AC_LANG_POP(Fortran)dnl
])
])# _AC_FC_NAME_MANGLING

# The replacement is empty.
AU_DEFUN([AC_FC_NAME_MANGLING], [])


# AC_FC_WRAPPERS
# ---------------
# Defines C macros FC_FUNC(name,NAME) and FC_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran compiler.
AC_DEFUN([AC_FC_WRAPPERS],
[AC_REQUIRE([_AC_FC_NAME_MANGLING])dnl
AH_TEMPLATE([FC_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran.])dnl
AH_TEMPLATE([FC_FUNC_],
    [As FC_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_fc_mangling in
  "lower case, no underscore, no extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [name])
          AC_DEFINE([FC_FUNC_(name,NAME)], [name]) ;;
  "lower case, no underscore, extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [name])
          AC_DEFINE([FC_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, no extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([FC_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([FC_FUNC_(name,NAME)], [name ## __]) ;;
  "upper case, no underscore, no extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([FC_FUNC_(name,NAME)], [NAME]) ;;
  "upper case, no underscore, extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([FC_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, no extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([FC_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, extra underscore")
          AC_DEFINE([FC_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([FC_FUNC_(name,NAME)], [NAME ## __]) ;;
  *)
          AC_MSG_WARN([unknown Fortran name-mangling scheme])
          ;;
esac
])# AC_FC_WRAPPERS


# AC_FC_FUNaC(NAME, [SHELLVAR = NAME])
# ------------------------------------
# For a Fortran subroutine of given NAME, define a shell variable
# $SHELLVAR to the Fortran-77 mangled name.  If the SHELLVAR
# argument is not supplied, it defaults to NAME.
AC_DEFUN([AC_FC_FUNC],
[AC_REQUIRE([_AC_FC_NAME_MANGLING])dnl
case $ac_cv_fc_mangling in
  upper*) ac_val="m4_toupper([$1])" ;;
  lower*) ac_val="m4_tolower([$1])" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_fc_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac
m4_if(m4_index([$1],[_]),-1,[],
[case $ac_cv_fc_mangling in *," extra underscore"*) ac_val="$ac_val"_ ;; esac
])
m4_default([$2],[$1])="$ac_val"
])# AC_FC_FUNC

