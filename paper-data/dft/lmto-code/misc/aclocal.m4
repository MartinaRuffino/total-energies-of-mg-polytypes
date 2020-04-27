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
rm -f conftest*
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

# --- ACX_MAKE_ARCHIVE_LIBRARY makes an archive library ---
# arguments :  directory-name
AC_DEFUN(ACX_MAKE_ARCHIVE_LIBRARY,
[

if test "$enable_$1" = "yes"; then
  thisdir=`pwd`
  AC_MSG_RESULT([creating soft link $1/Make.inc to Make.inc])
  (cd $1 && rm -f Make.inc && $LN_S $thisdir/Make.inc Make.inc)
  AC_MSG_RESULT([creating soft link $1/structures.h to structures.h])
  (cd $1 && rm -f structures.h && $LN_S $thisdir/structures.h structures.h)
  AC_MSG_RESULT([creating soft link $1/structures77.h to structures77.h])
  (cd $1 && rm -f structures77.h && $LN_S $thisdir/structures77.h structures77.h)
  AC_MSG_RESULT([creating soft link $1/interfaces.h to interfaces.h])
  (cd $1 && rm -f interfaces.h && $LN_S $thisdir/interfaces.h interfaces.h)
  AC_MSG_RESULT([creating soft link $1/interfaces77.h to interfaces77.h])
  (cd $1 && rm -f interfaces77.h && $LN_S $thisdir/interfaces77.h interfaces77.h)
  AC_MSG_RESULT([creating soft link $1/mod_ctx.mod to mod_ctx.mod])
  (cd $1 && rm -f mod_ctx.mod && $LN_S $thisdir/mods/mod_ctx.mod mod_ctx.mod)
dnl  if test -r $1/test; then
dnl    AC_MSG_RESULT([creating soft link $1/test/zdiff to testing/zdiff])
dnl    (cd $1/test && rm -f zdiff && $LN_S $thisdir/testing/zdiff zdiff)
dnl    AC_MSG_RESULT([creating soft link $1/test/add0 to testing/add0])
dnl    (cd $1/test && rm -f add0 && $LN_S $thisdir/testing/add0 add0)
dnl    AC_MSG_RESULT([creating soft link $1/test/poszer to testing/poszer])
dnl    (cd $1/test && rm -f poszer && $LN_S $thisdir/testing/poszer poszer)
dnl  fi
dnl  makefiles="$makefiles $1/Makefile"
dnl  AC_MSG_RESULT([Creating file $1/Makefile.in using misc/Maksmakfile:])
dnl  always call with initF90
dnl  MaksmakfileSW="--in"
dnl  if test "$FC_IS_F90" = "yes"; then MaksmakfileSW="$MaksmakfileSW --initF90"
dnl  fi
dnl  MaksmakfileSW="--in --initF90"
dnl  echo "(cd $1 && $thisdir/misc/Maksmakfile $MaksmakfileSW --mnemonic $MNEMONIC >Makefile.in)"
dnl  (cd $1 && $thisdir/misc/Maksmakfile $MaksmakfileSW --mnemonic $MNEMONIC >Makefile.in)

  AC_MSG_RESULT([Creating file $1/Make.patch using misc/Makepatchfile:])
  echo "(cd $1 && $thisdir/misc/Makepatchfile $MaksmakfileSW --mnemonic $MNEMONIC >Make.patch)"
  (cd $1 && $thisdir/misc/Makepatchfile  $MaksmakfileSW --mnemonic $MNEMONIC >Make.patch)
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
    AC_MSG_RESULT([creating soft link $1/test/poszer to testing/poszer])
    (cd $1/test && rm -f poszer && $LN_S $thisdir/testing/poszer poszer)
  fi
  makefiles="$makefiles $1/Makefile"
  AC_MSG_RESULT([Creating file $1/Makefile.in using misc/Maksmakfile:])
  echo "(cd $1 && $thisdir/misc/Maksmakfile --in --split --mnemonic $MNEMONIC >Makefile.in)"
  (cd $1 && $thisdir/misc/Maksmakfile --in --split --mnemonic $MNEMONIC >Makefile.in)
fi

])


# --- ACX_WARN_LIB_OVERLAP ---
AC_DEFUN(ACX_WARN_LIB_OVERLAP,
[

  if test "$1" = "$2" && test "${$4}" = yes; then
    echo "*** Warning: LIBLOC library $1 contains $3, yet $4=${$4}";
  fi

])

# --- ACX_CHECK_LIB_LINK checks if linker LK correctly links a file created by C compiler to a specified library ---
# $1 library name
# $2 C compiler
# $3 Fortran compiler
# $4 linker LK
# $5 action if found
# $6 action if not


AC_DEFUN(ACX_CHECK_LIB_LINK,
[
   AC_MSG_CHECKING([whether linker \"$4\" is able to link to LIBLOC])
  
  echo "void $ac_cv_fc_main(argc,argv)" > contest.c 
  echo '{printf("%d",argc);' >>contest.c
  echo "$FMAINTEST;}" >>contest.c
  echo '      subroutine fmain' > conftest.f
  echo '      end' >> conftest.f
  $2 -c $FCARGS contest.c >/dev/null 2>&1
  $3 -c $FCARGS conftest.f >/dev/null 2>&1 
#echo "void $ac_cv_fc_main(argc,argv)" > contest.c 
#echo '{printf("%d",argc);}' >>contest.c
#  $2 -c $FCARGS contest.c >/dev/null 2>&1
  $4 $FCARGS contest.o conftest.o $1 >/dev/null 2>&1 
 if test -x ./a.out; then 
      AC_MSG_RESULT([ok])
     ifelse([$5], , :, [$5])
  else
    AC_MSG_RESULT([no])
    ifelse([$6], , , [$6])
  fi
rm contest.c; rm contest.o
])

# --- ACX_CHECK_FOR_LIB checks if libraries are present by attempting to link with them ---
# $1 : library name
# $2 : action-if-found
# $3 : action-if-not-found
AC_DEFUN(ACX_CHECK_FOR_LIB,
[
  AC_MSG_CHECKING([whether library \"$1\" is in the link path])
  echo '      end' >conftest.f
  if `$FC $FCARGS conftest.f $1 >/dev/null 2>&1` ; then
     AC_MSG_RESULT([ok])
     ifelse([$2], , :, [$2])
  else
    AC_MSG_RESULT([no])
    ifelse([$3], , , [$3])
  fi
  rm -f conftest.f a.out
])



# --- ACX_FIND_ROUTINE_IN_LIB looks for a specific routine within a collection of libraries ---
# $1 : routine name
# $2 : libraries in which to find routine
# $3 : action-if-found
# $4 : action-if-not-found
AC_DEFUN(ACX_FIND_ROUTINE_IN_LIB,
[
  AC_MSG_CHECKING([for $1 in $2])
  echo '      call $1' >conftest.f
  echo '      end' >>conftest.f
  if `$FC $FCARGS conftest.f $2 >/dev/null 2>&1` ; then
     AC_MSG_RESULT([ok])
     ifelse([$3], , :, [$3])
  else
    AC_MSG_RESULT([no])
    ifelse([$4], , , [$4])
  fi
  rm -f conftest.f a.out

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
  if `$FC -c $FCARGS conftest.f >/dev/null 2>&1` ; then
    FC_IS_F90=yes
  else
    FC_IS_F90=no
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
if `$FC -c $FCARGS conftest.f $LIBLOC >/dev/null 2>&1` ; then
  FC_AUTOARRAY=yes
else
  FC_AUTOARRAY=no
fi
rm -f conftest.f conftest.o
])

# --- AC_CHECK_FOR_POINTER checks if the fortran compiler has POINTER extension ---
# Returns FC_POINTER="yes" or "no"
AC_DEFUN(AC_CHECK_FOR_POINTER,
[
cat > conftest.f <<EOF
      pointer (iarr, arr)
      double precision arr(20)
      end
EOF
if `$FC -c $FCARGS conftest.f >/dev/null 2>&1` ; then
  FC_POINTER=yes
else
  FC_POINTER=no
fi
rm -f conftest.f conftest.o
])


# --- AC_CHECK_FOR_COMPILER_EXTENSIONS checks whether the compiler is equippped ---
# for the enabled extensions.
# $1 package
# $2 compiler extensions
AC_DEFUN(AC_CHECK_FOR_COMPILER_EXTENSIONS,
[
echo $1
])

# AC_CHECK_FOR_GETARG checks if getarg is fortran-callable
# Returns FC_ALLOWS_GETARG="yes" or "no"
AC_DEFUN(AC_CHECK_FOR_GETARG,
[
cat > conftest.f <<EOF
      character*20 strn
      call getarg(1,strn)
      print '(a)', strn
      end
EOF
FC_ALLOWS_GETARG=no
if `$FC $FCARGS conftest.f >/dev/null 2>&1` ; then
  if (test `echo try | ./a.out try` = "try"); then FC_ALLOWS_GETARG=yes; fi
fi
rm -f conftest.f conftest.o
])

# AC_CHECK_ARG_LINK checks if either getarg or gtargc is fortran-callable through C main
# Returns FC_ARGS_LINK="yes" or "no"
# $1: the name of function/routine to check: "getarg" or "gtargc"
AC_DEFUN(AC_CHECK_ARG_LINK,
[
#first create test programs ctest.c and conftest.f, then compile and link
#ctest is a condensed version of fmain.c
cat > ctest.c << EOF1

/* Main and support programs to link unix shell with fortran programs.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI
#include "mpi.h"
#define MAXARGS 100
#endif

static int ret_val=0, aargc;
static char **pargv;


#if (LINUXI | LINUXA | LINUXF | LINUX_PGI | LINUX)
/* Note: switches should be done by the caller */
#define HAVE_IARGC 1
#endif

/* macro FC_FUNC converts a standard name to a fortranized version */
/*FC_FUNC must be bypassed for sun systems because c preprocessor doesn't recognize "##"*/
#if FC_UPPERCASE == 1
#  if FC_UNDERSCORE == 1
#    define FC_FUNC(x,X) X ## _
#  else
#    if FC_UNDERSCORE == 2
#      define FC_FUNC(x,X) X ## _  
#    else
#      define FC_FUNC(x,X) X
#    endif
#  endif
#else
#  if FC_UNDERSCORE == 1
#    define FC_FUNC(x,X) x ## _
#  else
#    if FC_UNDERSCORE == 2
#      define FC_FUNC(x,X) x ## _  
#    else
#      define FC_FUNC(x,X) x
#    endif
#  endif
#endif

/* --- Entry point --- */
#if NOPASSEDARGS == 0
CMAIN(argc,argv)
int argc; char *argv[];
#endif
#if NOPASSEDARGS == 1
CMAIN()
#endif

{
       
#if (FC_UNDERSCORE == 0 && FC_UPPERCASE == 1)
  void FMAIN();
#endif
#if (FC_UNDERSCORE == 1 && FC_UPPERCASE == 1)
  void FMAIN_();
#endif
#if (FC_UNDERSCORE == 2 && FC_UPPERCASE == 1)
  void FMAIN_();
#endif
#if (FC_UNDERSCORE == 0 && FC_UPPERCASE == 0)
  void fmain();
#endif
#if (FC_UNDERSCORE == 1 && FC_UPPERCASE == 0)
  void fmain_();
#endif
#if (FC_UNDERSCORE == 2 && FC_UPPERCASE == 0)
  void fmain_();
#endif
     

#if ! (NOPASSEDARGS == 1)
  aargc = argc;
  pargv = argv;
#endif

#ifdef AIX
   save_me();  /* so we won't get killed when page space is low */

#endif

/* --- Pass control to routine fmain --- */
/*fmain_();*/

FC_FUNC(fmain,FMAIN)();

 /* This is normal from fortran call */
   
 exit (ret_val);
   
}

/* --- function cexit: if *ps is nonzero, exit with retval pv --- */
/*void cexit_(pv,ps)*/

void FC_FUNC(cexit,CEXIT)(pv,ps)
int *pv,*ps;
{
  ret_val = *pv;
#ifdef CRAY2
  exit (ret_val);
#else
  if (*ps) exit (ret_val);
#endif
}

/* --- function nargc: retun the number of command-line arguments --- */
/*nargc_()*/

int FC_FUNC(nargc,NARGC)()
{
#if HAVE_IARGC == 1
  int i,iargc_();
  return(iargc_()+1);
#else
  return(aargc);
#endif
}


/* --- function gtargc: get command line arguments --- */
#if ! (NOPASSEDARGS == 1)
/*void gtargc_(iarg,ps,len)*/
void FC_FUNC(gtargc,GTARGC)(iarg,ps,len)
int *iarg; char *ps; short len;
{
  int i,maxlen; char *pps;

/* to handle fortran bug ... */
  len = (len < 0) ? -len : len;

  if (*iarg > aargc)
    { puts("getarg: request for nonexistent command line arg");
      exit(-1);
    }

/*copy string to ps, filling with blanks if passed string longer ...*/
  maxlen = strlen(*(pargv + *iarg));
  maxlen = (maxlen < len) ? maxlen : len;
  for (i = -1, pps=ps ; ++i<maxlen ;) *pps++ = *(*(pargv + *iarg) + i);
  while (i < len) {*pps++ = ' '; i++;}

}
#endif



#if DEC
void s_abort()
{
  exit(-1);
}
#endif

EOF1
# end of ctest.c

cat > conftest.f <<EOF2
      subroutine fmain()
      character*20 strn
      call $1(1,strn)
      print '(a)', strn
      end
EOF2
# end of conftest.f
$CC -c $FCARGS ctest.c -DCMAIN=$CMAIN -DNOPASSEDARGS=$NOPASSEDARGS -DFC_UNDERSCORE=$FC_UNDERSCORE -DFC_UPPERCASE=$FC_UPPERCASE >/dev/null 2>&1
FC_ARGS_LINK=no
$FC $FCARGS -c conftest.f >/dev/null 2>&1
$LK $FCARGS ctest.o conftest.o  #>/dev/null 2>&1
  if (test `echo try | ./a.out try` = "try"); then FC_ARGS_LINK=yes; fi

rm -f conftest.f conftest.o ctest.c ctest.o
])
 
# --- ACX_WARN_LIB_MISSING checks if a library is present by trying to ---
# link a call with a representative routine name.
# $1: part of -enable_$1
# $2: representative routine name
AC_DEFUN(ACX_WARN_LIB_MISSING,
[

# if test "${enable_$1}" = no; then
  echo '      call $2' >conftest.f
  echo '      end' >>conftest.f
  if `$FC $FCARGS conftest.f $LIBLOC >/dev/null 2>&1` ; then
    AC_MSG_RESULT([LIBLOC apparently contains $1 library]);
    have_$1=yes
  else
    AC_MSG_RESULT([LIBLOC apparently does not contain $1 library]);
  fi
  rm -f conftest.f a.out
# fi

])

# --- ACX_WARN_LINK_FAIL checks if a specified routine can be linked ---
# link a call with a representative routine name.
# $1 : routine name
# $2 : libraries to be included in the link
# $3 : if print, print message
# $4 : optional third argument: insert argument as first line of test routine
AC_DEFUN(ACX_WARN_LINK_FAIL,
[

  rm -f conftest.f
  ifelse([$4], , , [echo $4 >>conftest.f])
  echo '      call $1' >>conftest.f
  echo '      end' >>conftest.f
  ifelse([$3], "print", AC_MSG_CHECKING([whether $1 can be linked with libraries $2]), : )
  if `$FC $FCARGS conftest.f $2 >/dev/null 2>&1` ; then
    have_$1=yes
    ifelse([$3], "print", [echo yes], : )
  else
    ifelse([$3], "print", [echo no], : )
  fi
  rm -f conftest.f a.out

])

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

 rm -f conftest*
 LIBS=$ac_fc_m_save_LIBS
 AC_LANG_POP(C)dnl
])

AC_DEFINE_UNQUOTED([FC_CMAIN], $ac_cv_fc_maino,
                   [Define to name of a C langugate program that corresponds to
                    the main in a fortran program.])

])# AC_FC_CMAIN


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
rm -f conftest*])
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

rm -f conftest*
AC_LANG_POP(Fortran)dnl

# If we are using xlf then replace all the commas with spaces.
if echo $ac_fc_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_fc_v_output=`echo $ac_fc_v_output | sed 's/,/ /g'`
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
 rm -f conftest*
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
 rm -f conftest*
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


# AC_FC_FUNC(NAME, [SHELLVAR = NAME])
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